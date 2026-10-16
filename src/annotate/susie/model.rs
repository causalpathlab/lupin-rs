//! Positive-effect SuSiE on one cluster's marker log-folds.
//!
//! ```text
//! x_g = ln((y_g + ½) / (μ⁰_g + ½))     marker gene g's count against its background
//! η_g = Σ_c A_gc θ_c                    its predicted log-fold
//! x_g ~ t_ν(η_g, v_g)                   v_g = s²_g + 1/(y_g + ½)
//! θ = Σ_l e_l · 1[γ_l]                  e_l = exp(b_l) > 0
//! γ_l ∈ {1..C, null}                    one cell type per single effect, or none
//! ```
//!
//! `μ⁰_g` is the gene's background in this cluster (library size × background
//! rate). `s²_g` is how far the gene strays from any fit across clusters, and
//! `1/(y + ½)` the sampling variance of the observed log count. (The variance
//! at the predicted count, which would weigh a silent marker fully against a
//! type predicting it high, was measured and lost accuracy on real data, so a
//! type whose own markers are mostly near zero can still be credited for a
//! few shared genes; a test records that gap.) The Student-t (ν =
//! [`NU`]) bounds what one gene can contribute (Lange, Little & Taylor 1989,
//! "Robust statistical modeling using the t distribution", JASA 84:881–896,
//! doi:10.1080/01621459.1989.10478852): a type's effect is then the typical
//! fold of its markers, so a handful of genes it shares with what the cluster
//! really is cannot carry it while the rest of its markers stay flat.
//!
//! Each sweep draws every effect's type γ_l exactly from its full conditional
//! (C + 1 choices, each a change to the likelihood over that type's markers
//! only), then its size by elliptical slice sampling on the standard-normal
//! raw `b_l = m_b + s_b·z`: a log-normal size, a positive fold. Effects can
//! only add expression, so a type cannot be credited for its markers being
//! low, and an effect that finds nothing rests on the null choice. The draws
//! keep each sweep's conditional choice probabilities, so a type's PIP is
//! Rao-Blackwellised.

use legume_numeric::mcmc::engine::{elliptical_slice_step, McmcModel};
use rand::rngs::SmallRng;
use rand::{Rng, RngExt};
use rand_distr::{Distribution, StandardNormal};

/// Degrees of freedom of the Student-t on the log-folds.
pub(super) const NU: f64 = 4.0;

/// The prior of the single effects.
#[derive(Clone, Copy, Debug, serde::Serialize, serde::Deserialize)]
pub struct SusiePrior {
    /// Number of single effects L.
    pub num_effects: usize,
    /// Prior probability that an effect picks no type.
    pub null_weight: f32,
    /// Mean of log e, e the log fold an effect adds to its type's markers.
    pub effect_mean: f32,
    /// Standard deviation of log e.
    pub effect_sd: f32,
}

impl Default for SusiePrior {
    fn default() -> Self {
        Self {
            num_effects: 5,
            null_weight: 0.5,
            effect_mean: 0.0,
            effect_sd: 1.0,
        }
    }
}

/// What every cluster shares: which types each marker gene marks and how far
/// each gene strays from a fit across clusters.
pub(super) struct Panel {
    types_of: Vec<Vec<usize>>,
    genes_of: Vec<Vec<usize>>,
    /// s²_g, the gene's misfit variance on the log scale.
    scale2: Vec<f32>,
}

impl Panel {
    /// The panel over `types_of` (marker gene → its types), with s²_g from
    /// `scale2`.
    pub(super) fn new(types_of: Vec<Vec<usize>>, n_types: usize, scale2: &[f32]) -> Self {
        let mut genes_of = vec![Vec::new(); n_types];
        for (g, types) in types_of.iter().enumerate() {
            for &c in types {
                genes_of[c].push(g);
            }
        }
        Self {
            types_of,
            genes_of,
            scale2: scale2.to_vec(),
        }
    }

    /// Set s²_g of every marker gene.
    pub(super) fn set_scale2(&mut self, scale2: &[f32]) {
        self.scale2 = scale2.to_vec();
    }

    pub(super) fn n_types(&self) -> usize {
        self.genes_of.len()
    }

    /// Every marker gene's predicted log-fold at θ.
    pub(super) fn etas(&self, theta: &[f32]) -> Vec<f32> {
        self.types_of
            .iter()
            .map(|types| types.iter().map(|&c| theta[c]).sum::<f32>())
            .collect()
    }
}

/// A count `y` against its background expectation `mu0`: the log-fold
/// ln((y + ½)/(μ⁰ + ½)) and the sampling variance of the log count, 1/(y + ½).
pub(super) fn log_fold(y: f64, mu0: f64) -> (f32, f32) {
    (
        ((y + 0.5) / (mu0 + 0.5)).ln() as f32,
        (1.0 / (y + 0.5)) as f32,
    )
}

/// One cluster's marker genes against their background.
pub(super) struct Cluster {
    /// ln((y + ½) / (μ⁰ + ½)) per marker gene.
    pub(super) x: Vec<f32>,
    /// 1/(y + ½): the sampling variance of each observed log count.
    pub(super) sampling: Vec<f32>,
}

impl Cluster {
    /// From counts `y` and background expectations `mu0`.
    pub(super) fn new(y: impl Iterator<Item = f64>, mu0: impl Iterator<Item = f64>) -> Self {
        let (x, sampling) = y.zip(mu0).map(|(y, m)| log_fold(y, m)).unzip();
        Self { x, sampling }
    }
}

/// One cluster on the panel: the likelihood the sampler targets.
#[derive(Clone, Copy)]
pub(super) struct ClusterLogFold<'a> {
    pub(super) panel: &'a Panel,
    pub(super) cluster: &'a Cluster,
}

impl ClusterLogFold<'_> {
    /// Gene g's Student-t log-likelihood at predictor η, up to constants.
    fn term(&self, g: usize, eta: f32) -> f64 {
        let v = f64::from(self.panel.scale2[g]) + f64::from(self.cluster.sampling[g]);
        let r = f64::from(self.cluster.x[g]) - f64::from(eta);
        -0.5 * (NU + 1.0) * (r * r / (NU * v)).ln_1p()
    }

    fn log_lik(&self, theta: &[f32]) -> f64 {
        self.panel
            .etas(theta)
            .iter()
            .enumerate()
            .map(|(g, &eta)| self.term(g, eta))
            .sum()
    }

    /// Each type's share of the log-likelihood θ gains over the background:
    /// the log-likelihood lost when its effect alone is removed, over the
    /// log-likelihood θ gains on θ = 0 (0 for types without an effect).
    pub(super) fn explained(&self, theta: &[f32]) -> Vec<f32> {
        let total = self.log_lik(theta) - self.log_lik(&vec![0.0; theta.len()]);
        let etas = self.panel.etas(theta);
        (0..theta.len())
            .map(|c| {
                if theta[c] <= 0.0 || total <= 0.0 {
                    return 0.0;
                }
                let lost: f64 = self.panel.genes_of[c]
                    .iter()
                    .map(|&g| self.term(g, etas[g]) - self.term(g, etas[g] - theta[c]))
                    .sum();
                (lost / total).max(0.0) as f32
            })
            .collect()
    }

    /// The change in log-likelihood from adding `effect` to type c's markers,
    /// at the predictors `etas` whose terms are `base`.
    fn gain(&self, etas: &[f32], base: &[f64], c: usize, effect: f32) -> f64 {
        self.panel.genes_of[c]
            .iter()
            .map(|&g| self.term(g, etas[g] + effect) - base[g])
            .sum()
    }
}

fn std_normal(rng: &mut impl Rng) -> f32 {
    let v: f64 = StandardNormal.sample(rng);
    v as f32
}

/// The sampler over one cluster.
pub(super) struct SusieModel<'a> {
    pub(super) data: ClusterLogFold<'a>,
    pub(super) prior: SusiePrior,
}

/// Sampler state: per effect, its type (`C` = null), its standard-normal raw
/// size and that size, the last conditional choice probabilities; θ, and the
/// marker genes' η at θ with their likelihood terms, kept in step with it.
pub(super) struct SusieState {
    gamma: Vec<usize>,
    z_eff: Vec<f32>,
    effect: Vec<f32>,
    probs: Vec<Vec<f32>>,
    theta: Vec<f32>,
    etas: Vec<f32>,
    base: Vec<f64>,
}

/// One posterior draw.
#[derive(Clone)]
pub(super) struct SusieSample {
    /// Each effect's conditional choice probabilities over the C types and
    /// the null choice (last).
    pub(super) probs: Vec<Vec<f32>>,
    /// θ over the C types.
    pub(super) theta: Vec<f32>,
}

impl SusieModel<'_> {
    fn effect(&self, z: f32) -> f32 {
        (self.prior.effect_mean + self.prior.effect_sd * z)
            .min(5.0)
            .exp()
    }

    /// Log prior of each choice: the null weight, the rest shared evenly.
    fn log_prior(&self) -> (f64, f64) {
        let w = f64::from(self.prior.null_weight).clamp(1e-6, 1.0 - 1e-6);
        let c = self.data.panel.n_types() as f64;
        (((1.0 - w) / c).ln(), w.ln())
    }

    /// Add `delta` to type c's effect, keeping η and the terms in step.
    fn shift(&self, s: &mut SusieState, c: usize, delta: f32) {
        s.theta[c] += delta;
        for &g in &self.data.panel.genes_of[c] {
            s.etas[g] += delta;
            s.base[g] = self.data.term(g, s.etas[g]);
        }
    }
}

impl McmcModel for SusieModel<'_> {
    type State = SusieState;
    type Sample = SusieSample;
    type Result = Vec<SusieSample>;

    fn init(&self, rng: &mut SmallRng) -> SusieState {
        let (c, l) = (self.data.panel.n_types(), self.prior.num_effects);
        let z_eff: Vec<f32> = (0..l).map(|_| std_normal(rng)).collect();
        let effect = z_eff.iter().map(|&z| self.effect(z)).collect();
        let theta = vec![0.0; c];
        let etas = self.data.panel.etas(&theta);
        let base = etas
            .iter()
            .enumerate()
            .map(|(g, &e)| self.data.term(g, e))
            .collect();
        SusieState {
            // Every effect starts on the null choice: "no type".
            gamma: vec![c; l],
            z_eff,
            effect,
            probs: vec![vec![0.0; c + 1]; l],
            theta,
            etas,
            base,
        }
    }

    fn sweep(&self, s: &mut SusieState, rng: &mut SmallRng) {
        let c = self.data.panel.n_types();
        let (log_type, log_null) = self.log_prior();
        // η and the terms from θ afresh: the shifts below add and take away
        // effects in f32, and their rounding would otherwise accumulate.
        s.etas = self.data.panel.etas(&s.theta);
        s.base = s
            .etas
            .iter()
            .enumerate()
            .map(|(g, &e)| self.data.term(g, e))
            .collect();
        for l in 0..self.prior.num_effects {
            if s.gamma[l] < c {
                self.shift(s, s.gamma[l], -s.effect[l]);
            }

            // The type, exactly from its full conditional.
            let effect = s.effect[l];
            let gains: Vec<f64> = (0..c)
                .map(|j| self.data.gain(&s.etas, &s.base, j, effect))
                .collect();
            // A type with no marker genes is no choice.
            let mut logp: Vec<f64> = (0..c)
                .map(|j| {
                    if self.data.panel.genes_of[j].is_empty() {
                        f64::NEG_INFINITY
                    } else {
                        log_type + gains[j]
                    }
                })
                .collect();
            logp.push(log_null);
            let m = logp.iter().fold(f64::NEG_INFINITY, |m, &v| m.max(v));
            let w: Vec<f64> = logp.iter().map(|&v| (v - m).exp()).collect();
            let total: f64 = w.iter().sum();
            let p: Vec<f32> = w.iter().map(|&v| (v / total) as f32).collect();
            let u: f32 = rng.random();
            let mut cum = 0.0;
            let j = p
                .iter()
                .position(|&v| {
                    cum += v;
                    u < cum
                })
                .unwrap_or(c);
            s.gamma[l] = j;
            s.probs[l] = p;

            // Its size: ESS on the likelihood of its type's markers (starting
            // from the gain just computed), or a prior draw when it picked none.
            if j < c {
                // ESS takes an f32 log-likelihood: the gain relative to the
                // current one keeps it near 0.
                let cur = gains[j];
                let ll =
                    |z: &f32| (self.data.gain(&s.etas, &s.base, j, self.effect(*z)) - cur) as f32;
                let nu = std_normal(rng);
                s.z_eff[l] = elliptical_slice_step(&s.z_eff[l], &nu, &ll, 0.0, rng).0;
                s.effect[l] = self.effect(s.z_eff[l]);
                self.shift(s, j, s.effect[l]);
            } else {
                s.z_eff[l] = std_normal(rng);
                s.effect[l] = self.effect(s.z_eff[l]);
            }
        }
    }

    fn collect(&self, s: &SusieState) -> SusieSample {
        SusieSample {
            probs: s.probs.clone(),
            theta: s.theta.clone(),
        }
    }

    fn summarize(&self, samples: Vec<SusieSample>) -> Vec<SusieSample> {
        samples
    }
}
