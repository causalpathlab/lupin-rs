//! Count-based SuSiE annotation: each cluster's pseudobulk counts regressed on
//! the panel's marker sets at once, so cell types compete for the marker
//! genes they share, and each type can only add expression. See [`model`].
//!
//! The clusters are fitted independently given two things they share:
//!
//! - a marker gene's background rate, its median rate across clusters (a
//!   marker is high in few clusters, so the median sits at the level of the
//!   rest), times the cluster's library size;
//! - a marker gene's NB dispersion φ_g. The first pass uses one φ, the median
//!   method-of-moments estimate over the genes no type claims (the model
//!   leaves those at their mean in every cluster). The refit estimates φ_g
//!   from every cluster's residuals around the first fit's means, so a type's
//!   effect is not mistaken for noise, and shrinks it towards the genes'
//!   median on the log scale.
//!
//! A cluster's pseudobulk sums many cells, so most types that enter at all
//! have a PIP near 1; each type's effect (the log fold of its markers) and
//! its share of the explained deviance rank them.

pub mod model;

use legume_numeric::matrix::rand_util::mix_seed;
use legume_numeric::matrix::utils::{median, quantiles};
use legume_numeric::mcmc::engine::{run_mcmc_parallel, split_rhat_chains, McmcConfig};
use model::{Cluster, ClusterCounts, Panel, SusieModel, SusiePrior, SusieSample};
use rayon::prelude::*;

/// Coverage of a credible set.
const COVERAGE: f32 = 0.95;
/// Smallest |correlation| between two marker sets in one credible set.
const MIN_PURITY: f32 = 0.5;
/// Draws (and warmup) of the single-chain first pass, whose means only
/// feed the dispersion estimate.
const FIRST_PASS_DRAWS: usize = 200;
/// Pseudo-clusters of weight pulling each gene's log φ_g to the median.
const DISPERSION_SHRINKAGE: f64 = 10.0;
/// The range an estimated φ is clamped to.
const PHI_RANGE: (f64, f64) = (1e-3, 10.0);

/// Sampler settings for every cluster.
#[derive(Clone, Copy, Debug, serde::Serialize, serde::Deserialize)]
pub struct SusieConfig {
    pub prior: SusiePrior,
    /// Independent chains per cluster, pooled; their disagreement is R̂.
    pub chains: usize,
    /// Posterior samples per chain.
    pub samples: usize,
    pub warmup: usize,
    pub thin: usize,
    pub seed: u64,
    /// A type needs at least this many matched marker genes with counts to
    /// enter, as the enrichment's `--min-markers`.
    #[serde(default = "default_min_markers")]
    pub min_markers: usize,
    /// One NB dispersion φ for every gene, fixed (0: Poisson); when `None`,
    /// gene-specific and estimated, with one refit.
    pub dispersion: Option<f32>,
}

impl Default for SusieConfig {
    fn default() -> Self {
        Self {
            prior: SusiePrior::default(),
            chains: 4,
            samples: 1000,
            warmup: 500,
            thin: 1,
            seed: 42,
            min_markers: default_min_markers(),
            dispersion: None,
        }
    }
}

fn default_min_markers() -> usize {
    3
}

impl SusieConfig {
    /// Refuse settings that leave no posterior.
    pub fn validate(&self) -> anyhow::Result<()> {
        anyhow::ensure!(
            self.prior.num_effects >= 1 && self.samples >= 1 && self.chains >= 1 && self.thin >= 1,
            "SuSiE needs at least one single effect, chain and posterior sample, and thin ≥ 1"
        );
        anyhow::ensure!(
            self.dispersion.is_none_or(|d| d >= 0.0),
            "SuSiE's NB dispersion is ≥ 0 (0: Poisson)"
        );
        Ok(())
    }
}

/// The posterior of one cluster.
pub struct ClusterFit {
    /// Posterior inclusion probability per cell type.
    pub pip: Vec<f32>,
    /// Posterior mean effect θ per cell type (log fold of its markers).
    pub theta: Vec<f32>,
    /// Each type's share of the deviance θ explains.
    pub explained: Vec<f32>,
    /// Credible sets: cell-type indices, most probable first.
    pub credible_sets: Vec<Vec<usize>>,
    /// The largest split R̂ over the chains (each split in half) of the types' θ.
    pub max_rhat: f32,
}

impl ClusterFit {
    /// The fit of a cluster without counts: nothing included, no call.
    fn empty(n_types: usize) -> Self {
        Self {
            pip: vec![0.0; n_types],
            theta: vec![0.0; n_types],
            explained: vec![0.0; n_types],
            credible_sets: Vec::new(),
            max_rhat: 1.0,
        }
    }

    /// The cluster's call: its most probable type among those in a credible
    /// set, PIPs within 0.01 of each other tied and broken by the share of
    /// the cluster each explains; `None` without a credible set.
    #[must_use]
    pub fn call(&self) -> Option<usize> {
        self.credible_sets.iter().flatten().copied().reduce(|a, b| {
            let (pa, pb) = (self.pip[a], self.pip[b]);
            let b_wins = if (pa - pb).abs() < 0.01 {
                self.explained[b] > self.explained[a]
            } else {
                pb > pa
            };
            if b_wins {
                b
            } else {
                a
            }
        })
    }
}

/// All clusters' fits plus the dispersions they used.
pub struct SusieFits {
    pub clusters: Vec<ClusterFit>,
    /// φ_g of each marker gene with counts (the likelihood's rows) in the
    /// last pass.
    pub dispersion: Vec<f32>,
}

/// Fit every cluster. `gene_sum_kg` is the G × K raw count sums, column-major
/// (cluster k is `[k·G .. (k+1)·G]`); `markers` lists each gene's cell types.
pub fn fit_all(
    gene_sum_kg: &[f64],
    n_genes: usize,
    n_clusters: usize,
    markers: &[Vec<usize>],
    n_types: usize,
    cfg: &SusieConfig,
) -> anyhow::Result<SusieFits> {
    anyhow::ensure!(
        gene_sum_kg.len() == n_genes * n_clusters && markers.len() == n_genes,
        "count sums and marker rows disagree on the genes"
    );
    cfg.validate()?;
    // A cluster id with no cells (a gap in the ids) has no counts: it is left
    // out of the shared background and dispersion, and gets no call.
    let column = |k: usize| &gene_sum_kg[k * n_genes..(k + 1) * n_genes];
    let present: Vec<usize> = (0..n_clusters)
        .filter(|&k| column(k).iter().any(|&v| v > 0.0))
        .collect();
    let at = |g: usize, k: usize| gene_sum_kg[k * n_genes + g];
    let expressed: Vec<usize> = (0..n_genes)
        .filter(|&g| present.iter().any(|&k| at(g, k) > 0.0))
        .collect();
    // A type with too few matched markers is left out, as the enrichment
    // leaves it out.
    let mut matched = vec![0usize; n_types];
    for &g in &expressed {
        markers[g].iter().for_each(|&c| matched[c] += 1);
    }
    let kept = |g: usize| -> Vec<usize> {
        markers[g]
            .iter()
            .copied()
            .filter(|&c| matched[c] >= cfg.min_markers)
            .collect()
    };
    let (marker_genes, free): (Vec<usize>, Vec<usize>) =
        expressed.iter().partition(|&&g| !kept(g).is_empty());
    anyhow::ensure!(!marker_genes.is_empty(), "no marker gene has counts");
    // A cluster's size: its counts of the genes no type claims, which a type's
    // markers rising cannot inflate (all its counts when no gene is free).
    let lib: Vec<f64> = present
        .iter()
        .map(|&k| {
            let free_sum: f64 = free.iter().map(|&g| at(g, k)).sum();
            if free_sum > 0.0 {
                free_sum
            } else {
                column(k).iter().sum()
            }
        })
        .collect();
    let total: f64 = lib.iter().sum();

    // The genes no type claims, at their pooled rate in every cluster.
    let global = cfg.dispersion.unwrap_or_else(|| {
        let est: Vec<f32> = free
            .iter()
            .filter_map(|&g| {
                let y = || present.iter().map(move |&k| at(g, k));
                let p = y().sum::<f64>() / total;
                // Genes too sparse to say anything are skipped.
                (p * total >= 10.0 * present.len() as f64)
                    .then(|| moment_phi(y(), lib.iter().map(|&n| n * p)) as f32)
            })
            .collect();
        median(&est)
    });

    // Each marker gene's rate in each present cluster, floored so a gene
    // silent in a cluster keeps a finite log.
    let floor = (0.5 / total) as f32;
    let rates: Vec<Vec<f32>> = marker_genes
        .iter()
        .map(|&g| {
            present
                .iter()
                .zip(&lib)
                .map(|(&k, &n)| ((at(g, k) / n) as f32).max(floor))
                .collect()
        })
        .collect();
    // The clusters on a background: each marker gene's rate × library size.
    let on = |background: &[f32]| -> Vec<Cluster> {
        present
            .iter()
            .zip(&lib)
            .map(|(&k, &n)| Cluster {
                y: marker_genes.iter().map(|&g| at(g, k) as f32).collect(),
                log_mu0: background
                    .iter()
                    .map(|&b| (f64::from(b) * n).ln() as f32)
                    .collect(),
            })
            .collect()
    };
    let types_of: Vec<Vec<usize>> = marker_genes.iter().map(|&g| kept(g)).collect();
    let corr = marker_correlation(&types_of, n_types);

    // Each present cluster's fit under `cfg`, seeded by its own id.
    let fit_pass = |panel: &Panel, data: &[Cluster], cfg: &SusieConfig| -> Vec<ClusterFit> {
        data.par_iter()
            .zip(&present)
            .map(|(cluster, &k)| {
                let counts = ClusterCounts { panel, cluster };
                let model = SusieModel {
                    data: counts,
                    prior: cfg.prior,
                };
                let config = McmcConfig {
                    n_samples: cfg.samples,
                    warmup: cfg.warmup,
                    thin: cfg.thin,
                    seed: mix_seed(cfg.seed, k as u64),
                };
                summarize(
                    &run_mcmc_parallel(&model, &config, cfg.chains),
                    counts,
                    &corr,
                )
            })
            .collect()
    };

    // A light first pass, on each marker's lower-quartile rate: a type in up
    // to three quarters of the clusters still shows above it.
    let quartile: Vec<f32> = rates.iter().map(|r| quantiles(r, &[0.25])[0]).collect();
    let first_data = on(&quartile);
    let mut dispersion = vec![global; marker_genes.len()];
    let mut panel = Panel::new(types_of.clone(), n_types, &dispersion);
    let first = fit_pass(
        &panel,
        &first_data,
        &SusieConfig {
            chains: 1,
            samples: FIRST_PASS_DRAWS,
            warmup: FIRST_PASS_DRAWS,
            thin: 1,
            ..*cfg
        },
    );
    // The background: a marker's median rate over the clusters the first
    // pass called none of its types (so a type in most clusters does not set
    // its own), else the lower quartile.
    let entered: Vec<Vec<bool>> = first
        .iter()
        .map(|f| {
            let mut e = vec![false; n_types];
            if let Some(c) = f.call() {
                e[c] = true;
            }
            e
        })
        .collect();
    let background: Vec<f32> = rates
        .iter()
        .zip(&types_of)
        .zip(&quartile)
        .map(|((r, types), &q)| {
            let absent: Vec<f32> = r
                .iter()
                .zip(&entered)
                .filter(|(_, e)| types.iter().all(|&c| !e[c]))
                .map(|(&v, _)| v)
                .collect();
            if absent.is_empty() {
                q
            } else {
                median(&absent)
            }
        })
        .collect();
    if cfg.dispersion.is_none() {
        // φ_g from every cluster's residuals around the first pass's means.
        let etas: Vec<Vec<f32>> = first_data
            .iter()
            .zip(&first)
            .map(|(cluster, f)| {
                ClusterCounts {
                    panel: &panel,
                    cluster,
                }
                .etas(&f.theta)
            })
            .collect();
        let raw: Vec<f64> = (0..marker_genes.len())
            .map(|i| {
                let y = first_data.iter().map(|c| f64::from(c.y[i]));
                let mu = etas.iter().map(|eta| f64::from(eta[i].exp()));
                moment_phi(y, mu)
            })
            .collect();
        dispersion = shrunk_dispersion(&raw, present.len());
        panel.set_dispersion(&dispersion);
    }
    let fitted = fit_pass(&panel, &on(&background), cfg);
    // Back to every id; one without counts gets no call.
    let mut fitted = fitted.into_iter();
    let clusters = (0..n_clusters)
        .map(|k| {
            if present.binary_search(&k).is_ok() {
                fitted.next().expect("one fit per present cluster")
            } else {
                ClusterFit::empty(n_types)
            }
        })
        .collect();
    Ok(SusieFits {
        clusters,
        dispersion,
    })
}

/// φ in Var = μ + φμ² by moments: Σ((y − μ)² − μ) / Σμ², clamped to
/// [`PHI_RANGE`] (its floor when nothing is expected).
fn moment_phi(y: impl Iterator<Item = f64>, mu: impl Iterator<Item = f64>) -> f64 {
    let (num, den) = y.zip(mu).fold((0.0, 0.0), |(n, d), (y, m)| {
        (n + (y - m).powi(2) - m, d + m * m)
    });
    if den > 0.0 {
        (num / den).clamp(PHI_RANGE.0, PHI_RANGE.1)
    } else {
        PHI_RANGE.0
    }
}

/// Per-gene φ shrunk on the log scale towards the genes' median, with
/// [`DISPERSION_SHRINKAGE`] pseudo-clusters of weight against `n_clusters`.
fn shrunk_dispersion(raw: &[f64], n_clusters: usize) -> Vec<f32> {
    let as_f32: Vec<f32> = raw.iter().map(|&d| d as f32).collect();
    let center = f64::from(median(&as_f32)).ln();
    let (k, w) = (n_clusters as f64, DISPERSION_SHRINKAGE);
    raw.iter()
        .map(|&d| {
            ((k * d.ln() + w * center) / (k + w))
                .exp()
                .clamp(PHI_RANGE.0, PHI_RANGE.1) as f32
        })
        .collect()
}

/// |Pearson correlation| between the types' marker indicators over the
/// marker genes.
fn marker_correlation(types_of: &[Vec<usize>], n_types: usize) -> Vec<Vec<f32>> {
    let n = types_of.len() as f64;
    let mut count = vec![0.0f64; n_types];
    let mut both = vec![vec![0.0f64; n_types]; n_types];
    for types in types_of {
        for &a in types {
            count[a] += 1.0;
            for &b in types {
                both[a][b] += 1.0;
            }
        }
    }
    (0..n_types)
        .map(|a| {
            (0..n_types)
                .map(|b| {
                    let (pa, pb) = (count[a] / n, count[b] / n);
                    let cov = both[a][b] / n - pa * pb;
                    let sd = (pa * (1.0 - pa) * pb * (1.0 - pb)).sqrt();
                    if sd > 0.0 {
                        (cov / sd).abs() as f32
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect()
}

/// Each single effect's mean choice probabilities over all draws. Draws
/// number the same effects differently (across chains, and within one as
/// effects trade types), so each draw's effects are matched to a reference,
/// greedily by overlap, before they are averaged: first to the first draw,
/// then to that average. PIP and θ do not depend on the numbering, an
/// effect's credible set does.
fn aligned_choice_means(chains: &[Vec<SusieSample>]) -> Vec<Vec<f32>> {
    let draws: Vec<&Vec<Vec<f32>>> = chains.iter().flatten().map(|s| &s.probs).collect();
    let Some(first) = draws.first() else {
        return Vec::new();
    };
    let (effects, choices) = (first.len(), first.first().map_or(0, Vec::len));
    let mut reference: Vec<Vec<f32>> = (*first).clone();
    for _ in 0..2 {
        let mut total = vec![vec![0.0f32; choices]; effects];
        for probs in &draws {
            for (r, c) in match_effects(&reference, probs) {
                for (v, &x) in total[r].iter_mut().zip(&probs[c]) {
                    *v += x / draws.len() as f32;
                }
            }
        }
        reference = total;
    }
    reference
}

/// Pairs (reference effect, draw effect), the most overlapping free pair
/// first.
fn match_effects(reference: &[Vec<f32>], draw: &[Vec<f32>]) -> Vec<(usize, usize)> {
    let dot = |a: &[f32], b: &[f32]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f32>();
    let n = reference.len().min(draw.len());
    let mut pairs: Vec<(usize, usize, f32)> = (0..n)
        .flat_map(|r| (0..n).map(move |c| (r, c)))
        .map(|(r, c)| (r, c, dot(&reference[r], &draw[c])))
        .collect();
    pairs.sort_by(|a, b| b.2.total_cmp(&a.2));
    let (mut used_r, mut used_c) = (vec![false; n], vec![false; n]);
    let mut out = Vec::with_capacity(n);
    for (r, c, _) in pairs {
        if !used_r[r] && !used_c[c] {
            (used_r[r], used_c[c]) = (true, true);
            out.push((r, c));
        }
    }
    out
}

/// PIP, mean θ, the share each type explains, credible sets and R̂ from one
/// cluster's chains, pooled.
fn summarize(
    chains: &[Vec<SusieSample>],
    counts: ClusterCounts<'_>,
    corr: &[Vec<f32>],
) -> ClusterFit {
    let samples: Vec<&SusieSample> = chains.iter().flatten().collect();
    let n_types = corr.len();
    let t = samples.len().max(1) as f32;
    let mut pip = vec![0.0f32; n_types];
    let mut theta = vec![0.0f32; n_types];
    for s in &samples {
        for c in 0..n_types {
            let excluded: f32 = s.probs.iter().map(|p| 1.0 - p[c]).product();
            pip[c] += (1.0 - excluded) / t;
            theta[c] += s.theta[c] / t;
        }
    }
    let alpha_bar = aligned_choice_means(chains);
    let mut credible_sets: Vec<Vec<usize>> = Vec::new();
    for bar in &alpha_bar {
        let null = bar[n_types];
        if null > 0.5 {
            continue;
        }
        // The effect's mass on the types, renormalized without the null.
        let on_types: Vec<f32> = bar[..n_types].iter().map(|v| v / (1.0 - null)).collect();
        let Some(set) = enrichment::consensus::credible_set(&on_types, COVERAGE, n_types) else {
            continue;
        };
        let pure = set
            .iter()
            .all(|&a| set.iter().all(|&b| a == b || corr[a][b] >= MIN_PURITY));
        if pure && !credible_sets.contains(&set) {
            credible_sets.push(set);
        }
    }
    let max_rhat = (0..n_types)
        .map(|c| {
            let per_chain: Vec<Vec<f32>> = chains
                .iter()
                .map(|chain| chain.iter().map(|s| s.theta[c]).collect())
                .collect();
            split_rhat_chains(&per_chain)
        })
        .fold(1.0f32, f32::max);
    ClusterFit {
        explained: counts.explained(&theta),
        pip,
        theta,
        credible_sets,
        max_rhat,
    }
}

#[cfg(test)]
#[path = "../tests/susie.rs"]
mod tests;
