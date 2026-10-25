//! SuSiE annotation on counts: each cluster's marker genes, as log-folds
//! against their background, regressed on the panel's marker sets at once, so
//! cell types compete for the marker genes they share, and each type can only
//! add expression. See [`model`]. SuSiE: Wang, Sarkar, Carbonetto & Stephens
//! (2020), "A simple new approach to variable selection in regression, with
//! application to genetic fine mapping", JRSS-B 82(5):1273–1300,
//! doi:10.1111/rssb.12388.
//!
//! The clusters are fitted independently given two things they share:
//!
//! - a marker gene's background rate, its median rate across the clusters the
//!   first pass called none of its types in (a marker is high in few
//!   clusters, so the median sits at the level of the rest), times the
//!   cluster's library size;
//! - a marker gene's misfit scale s²_g on the log scale. The first pass uses
//!   one s², the median robust (MAD) variance of the genes no type claims
//!   around their pooled rate — how far genes stray between clusters for
//!   reasons no type explains. The refit estimates s²_g from every cluster's
//!   residual log-folds around the first fit, and shrinks it towards the
//!   genes' median on the log scale.
//!
//! The Student-t likelihood bounds any one gene's pull, so a type's effect is
//! the typical fold of its markers; each type's effect and its share of the
//! explained log-likelihood rank the types that enter.

pub mod model;

use legume_numeric::matrix::rand_util::mix_seed;
use legume_numeric::matrix::utils::{median, quantiles};
use legume_numeric::mcmc::engine::{run_mcmc_parallel, split_rhat_chains, McmcConfig};
use model::{log_fold, Cluster, ClusterLogFold, Panel, SusieModel, SusiePrior, SusieSample};
use rayon::prelude::*;

/// Coverage of a credible set.
const COVERAGE: f32 = 0.95;
/// Smallest |correlation| between two marker sets in one credible set.
const MIN_PURITY: f32 = 0.5;
/// Draws (and warmup) of the single-chain first pass, whose fit only feeds
/// the background and the misfit scales.
const FIRST_PASS_DRAWS: usize = 200;
/// Pseudo-clusters of weight pulling each gene's log s²_g to the median.
const SCALE_SHRINKAGE: f64 = 10.0;
/// The range an estimated s² is clamped to (log-scale variance).
const SCALE2_RANGE: (f64, f64) = (1e-3, 25.0);

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
        Ok(())
    }
}

/// The posterior of one cluster.
pub struct ClusterFit {
    /// Posterior inclusion probability per cell type.
    pub pip: Vec<f32>,
    /// Posterior mean effect θ per cell type (log fold of its markers).
    pub theta: Vec<f32>,
    /// Each type's share of the log-likelihood θ gains; the shares sum to at
    /// most 1.
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
        let mut members: Vec<usize> = self.credible_sets.iter().flatten().copied().collect();
        call_order(&mut members, |&k| (self.pip[k], self.explained[k]));
        members.first().copied()
    }
}

/// PIPs closer than this are tied in SuSiE's call, the tie broken by the
/// share of the cluster each type explains.
const PIP_TIE: f32 = 0.01;

/// Put `items` in SuSiE's call order by their (PIP, explained share), best
/// first: by PIP, each run of PIPs within [`PIP_TIE`] below the run's first
/// taken as tied and ordered by the larger explained share. A total order
/// (the tie alone is not transitive); a NaN PIP comes last.
pub fn call_order<T>(items: &mut [T], key: impl Fn(&T) -> (f32, f32)) {
    let pip = |t: &T| {
        Some(key(t).0)
            .filter(|p| !p.is_nan())
            .unwrap_or(f32::NEG_INFINITY)
    };
    items.sort_by(|a, b| pip(b).total_cmp(&pip(a)));
    let mut start = 0;
    while start < items.len() {
        let lead = pip(&items[start]);
        // At least the lead itself, whose own gap is NaN when it is −∞.
        let end = start
            + items[start..]
                .partition_point(|t| lead - pip(t) < PIP_TIE)
                .max(1);
        items[start..end].sort_by(|a, b| key(b).1.total_cmp(&key(a).1));
        start = end;
    }
}

/// All clusters' fits plus the misfit scales they used.
pub struct SusieFits {
    pub clusters: Vec<ClusterFit>,
    /// s²_g (log-scale misfit variance) of each marker gene with counts (the
    /// likelihood's rows) in the last pass.
    pub scale2: Vec<f32>,
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
    // out of the shared background and scales, and gets no call.
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

    // How far the genes no type claims stray from their pooled rate.
    let global = {
        let est: Vec<f32> = free
            .iter()
            .filter_map(|&g| {
                let y: Vec<f64> = present.iter().map(|&k| at(g, k)).collect();
                let p = y.iter().sum::<f64>() / total;
                // Genes too sparse to say anything are skipped.
                (p * total >= 10.0 * present.len() as f64).then(|| {
                    let (x, v): (Vec<f32>, Vec<f32>) = y
                        .iter()
                        .zip(&lib)
                        .map(|(&y, &n)| log_fold(y, n * p))
                        .unzip();
                    misfit_scale2(&x, &v) as f32
                })
            })
            .collect();
        // Without a gene to measure it on: a unit log-scale variance.
        if est.is_empty() {
            1.0
        } else {
            median(&est)
        }
    };

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
            .map(|(&k, &n)| {
                Cluster::new(
                    marker_genes.iter().map(|&g| at(g, k)),
                    background.iter().map(|&b| f64::from(b) * n),
                )
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
                let lik = ClusterLogFold { panel, cluster };
                let model = SusieModel {
                    data: lik,
                    prior: cfg.prior,
                };
                let config = McmcConfig {
                    n_samples: cfg.samples,
                    warmup: cfg.warmup,
                    thin: cfg.thin,
                    seed: mix_seed(cfg.seed, k as u64),
                };
                summarize(&run_mcmc_parallel(&model, &config, cfg.chains), lik, &corr)
            })
            .collect()
    };

    // A light first pass, on each marker's lower-quartile rate: a type in up
    // to three quarters of the clusters still shows above it.
    let quartile: Vec<f32> = rates.iter().map(|r| quantiles(r, &[0.25])[0]).collect();
    let first_data = on(&quartile);
    let mut panel = Panel::new(types_of.clone(), n_types, &vec![global; marker_genes.len()]);
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
    // s²_g from every cluster's residual log-folds around the first fit. The
    // fit removes the types the first pass called; the MAD then keeps the
    // clusters it missed a type in (one present in most clusters, say) from
    // setting the spread.
    let fitted_etas: Vec<Vec<f32>> = first.iter().map(|f| panel.etas(&f.theta)).collect();
    let raw: Vec<f64> = (0..marker_genes.len())
        .map(|i| {
            let r: Vec<f32> = first_data
                .iter()
                .zip(&fitted_etas)
                .map(|(c, eta)| c.x[i] - eta[i])
                .collect();
            let v: Vec<f32> = first_data.iter().map(|c| c.sampling[i]).collect();
            misfit_scale2(&r, &v)
        })
        .collect();
    let scale2 = shrunk_scale2(&raw, present.len());
    panel.set_scale2(&scale2);
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
    Ok(SusieFits { clusters, scale2 })
}

/// A gene's misfit variance on the log scale: the squared MAD-based spread of
/// its log-folds `x` across clusters, less the median of their sampling
/// variances `v`, clamped to [`SCALE2_RANGE`]. Both sides robust: the MAD
/// keeps the clusters the gene is truly up in from setting its spread, the
/// median keeps a few clusters where it reads zero from cancelling it.
fn misfit_scale2(x: &[f32], v: &[f32]) -> f64 {
    let sd = f64::from(crate::annotate::geneset_units::mad_sd(x));
    (sd * sd - f64::from(median(v))).clamp(SCALE2_RANGE.0, SCALE2_RANGE.1)
}

/// Per-gene s² shrunk on the log scale towards the median of the genes
/// measured above the floor (genes too sparse to measure sit at it and say
/// nothing about the rest), with [`SCALE_SHRINKAGE`] pseudo-clusters of
/// weight against `n_clusters`.
fn shrunk_scale2(raw: &[f64], n_clusters: usize) -> Vec<f32> {
    let measured: Vec<f32> = raw
        .iter()
        .filter(|&&d| d > SCALE2_RANGE.0)
        .map(|&d| d as f32)
        .collect();
    let center = if measured.is_empty() {
        SCALE2_RANGE.0.ln()
    } else {
        f64::from(median(&measured)).ln()
    };
    let (k, w) = (n_clusters as f64, SCALE_SHRINKAGE);
    raw.iter()
        .map(|&d| {
            // In range: a log-scale weighted mean of in-range values.
            ((k * d.ln() + w * center) / (k + w)).exp() as f32
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
    // Summed in f64 and divided once: thousands of f32 terms of x/n drift.
    let n = draws.len() as f64;
    for _ in 0..2 {
        let mut total = vec![vec![0.0f64; choices]; effects];
        for probs in &draws {
            for (r, c) in match_effects(&reference, probs) {
                for (v, &x) in total[r].iter_mut().zip(&probs[c]) {
                    *v += f64::from(x);
                }
            }
        }
        reference = total
            .iter()
            .map(|row| row.iter().map(|&v| (v / n) as f32).collect())
            .collect();
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
    data: ClusterLogFold<'_>,
    corr: &[Vec<f32>],
) -> ClusterFit {
    let samples: Vec<&SusieSample> = chains.iter().flatten().collect();
    let n_types = corr.len();
    // Averaged in f64: thousands of f32 terms of 1/t drift past 1.
    let t = samples.len().max(1) as f64;
    let mut pip_sum = vec![0.0f64; n_types];
    let mut theta_sum = vec![0.0f64; n_types];
    for s in &samples {
        for c in 0..n_types {
            let excluded: f64 = s.probs.iter().map(|p| 1.0 - f64::from(p[c])).product();
            pip_sum[c] += 1.0 - excluded;
            theta_sum[c] += f64::from(s.theta[c]);
        }
    }
    let pip: Vec<f32> = pip_sum.iter().map(|&v| (v / t) as f32).collect();
    let theta: Vec<f32> = theta_sum.iter().map(|&v| (v / t) as f32).collect();
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
        explained: data.explained(&theta),
        pip,
        theta,
        credible_sets,
        max_rhat,
    }
}

#[cfg(test)]
#[path = "../tests/susie.rs"]
mod tests;
