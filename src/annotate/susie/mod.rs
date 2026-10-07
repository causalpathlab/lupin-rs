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
use legume_numeric::matrix::utils::median;
use legume_numeric::mcmc::engine::{run_mcmc, split_rhat, McmcConfig};
use model::{Cluster, ClusterCounts, Panel, SusieModel, SusiePrior, SusieSample};
use rayon::prelude::*;

/// Coverage of a credible set.
const COVERAGE: f32 = 0.95;
/// Smallest |correlation| between two marker sets in one credible set.
const MIN_PURITY: f32 = 0.5;
/// Pseudo-clusters of weight pulling each gene's log φ_g to the median.
const DISPERSION_SHRINKAGE: f64 = 10.0;
/// The range an estimated φ is clamped to.
const PHI_RANGE: (f64, f64) = (1e-3, 10.0);

/// Sampler settings for every cluster.
#[derive(Clone, Copy, Debug, serde::Serialize, serde::Deserialize)]
pub struct SusieConfig {
    pub prior: SusiePrior,
    pub samples: usize,
    pub warmup: usize,
    pub thin: usize,
    pub seed: u64,
    /// One NB dispersion φ for every gene, fixed (0: Poisson); when `None`,
    /// gene-specific and estimated, with one refit.
    pub dispersion: Option<f32>,
}

impl Default for SusieConfig {
    fn default() -> Self {
        Self {
            prior: SusiePrior::default(),
            samples: 1000,
            warmup: 500,
            thin: 1,
            seed: 42,
            dispersion: None,
        }
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
    /// The largest split-R̂ over the types' θ.
    pub max_rhat: f32,
}

impl ClusterFit {
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
    /// Marker genes with counts, the likelihood's rows.
    pub n_marker_genes: usize,
    /// φ_g of each of those genes in the last pass.
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
    let at = |g: usize, k: usize| gene_sum_kg[k * n_genes + g];
    let lib: Vec<f64> = (0..n_clusters)
        .map(|k| {
            gene_sum_kg[k * n_genes..(k + 1) * n_genes]
                .iter()
                .sum::<f64>()
                .max(1.0)
        })
        .collect();
    let total: f64 = lib.iter().sum();
    let expressed: Vec<usize> = (0..n_genes)
        .filter(|&g| (0..n_clusters).any(|k| at(g, k) > 0.0))
        .collect();
    let (marker_genes, free): (Vec<usize>, Vec<usize>) =
        expressed.iter().partition(|&&g| !markers[g].is_empty());
    anyhow::ensure!(!marker_genes.is_empty(), "no marker gene has counts");

    // The genes no type claims, at their pooled rate in every cluster.
    let global = cfg.dispersion.unwrap_or_else(|| {
        let est: Vec<f32> = free
            .iter()
            .filter_map(|&g| {
                let p = (0..n_clusters).map(|k| at(g, k)).sum::<f64>() / total;
                // Genes too sparse to say anything are skipped.
                (p * total >= 10.0 * n_clusters as f64).then(|| {
                    let mu = lib.iter().map(|&n| n * p);
                    moment_phi((0..n_clusters).map(|k| at(g, k)), mu) as f32
                })
            })
            .collect();
        median(&est)
    });

    // A marker gene's background rate: its median rate across clusters,
    // floored so a gene silent in most clusters keeps a finite log.
    let floor = (0.5 / total) as f32;
    let background: Vec<f32> = marker_genes
        .iter()
        .map(|&g| {
            let rates: Vec<f32> = (0..n_clusters)
                .map(|k| (at(g, k) / lib[k]) as f32)
                .collect();
            median(&rates).max(floor)
        })
        .collect();
    let data: Vec<Cluster> = (0..n_clusters)
        .map(|k| Cluster {
            y: marker_genes.iter().map(|&g| at(g, k) as f32).collect(),
            log_mu0: background
                .iter()
                .map(|&b| (f64::from(b) * lib[k]).ln() as f32)
                .collect(),
        })
        .collect();
    let types_of: Vec<Vec<usize>> = marker_genes.iter().map(|&g| markers[g].clone()).collect();
    let corr = marker_correlation(&types_of, n_types);

    let fit_pass = |panel: &Panel| -> Vec<ClusterFit> {
        data.par_iter()
            .enumerate()
            .map(|(k, cluster)| {
                let counts = ClusterCounts { panel, cluster };
                let model = SusieModel {
                    data: counts,
                    prior: cfg.prior,
                };
                let config = McmcConfig {
                    n_samples: cfg.samples,
                    warmup: cfg.warmup,
                    thin: cfg.thin.max(1),
                    seed: mix_seed(cfg.seed, k as u64),
                };
                summarize(&run_mcmc(&model, &config), counts, &corr)
            })
            .collect()
    };

    let mut dispersion = vec![global; marker_genes.len()];
    let mut panel = Panel::new(types_of.clone(), n_types, &dispersion);
    let mut clusters = fit_pass(&panel);
    if cfg.dispersion.is_none() {
        // φ_g from every cluster's residuals around the fitted means.
        let fitted: Vec<Vec<f32>> = data
            .iter()
            .zip(&clusters)
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
                let y = data.iter().map(|c| f64::from(c.y[i]));
                let mu = fitted.iter().map(|eta| f64::from(eta[i].exp()));
                moment_phi(y, mu)
            })
            .collect();
        dispersion = shrunk_dispersion(&raw, n_clusters);
        panel = Panel::new(types_of, n_types, &dispersion);
        clusters = fit_pass(&panel);
    }
    Ok(SusieFits {
        clusters,
        n_marker_genes: marker_genes.len(),
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

/// PIP, mean θ, the share each type explains, credible sets and R̂ from one
/// cluster's draws.
fn summarize(samples: &[SusieSample], counts: ClusterCounts<'_>, corr: &[Vec<f32>]) -> ClusterFit {
    let n_types = corr.len();
    let t = samples.len().max(1) as f32;
    let num_effects = samples.first().map_or(0, |s| s.probs.len());
    let mut pip = vec![0.0f32; n_types];
    let mut theta = vec![0.0f32; n_types];
    let mut alpha_bar = vec![vec![0.0f32; n_types + 1]; num_effects];
    for s in samples {
        for c in 0..n_types {
            let excluded: f32 = s.probs.iter().map(|p| 1.0 - p[c]).product();
            pip[c] += (1.0 - excluded) / t;
            theta[c] += s.theta[c] / t;
        }
        for (bar, p) in alpha_bar.iter_mut().zip(&s.probs) {
            for (b, &v) in bar.iter_mut().zip(p) {
                *b += v / t;
            }
        }
    }
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
            let x: Vec<f32> = samples.iter().map(|s| s.theta[c]).collect();
            split_rhat(&x)
        })
        .filter(|r| r.is_finite())
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
