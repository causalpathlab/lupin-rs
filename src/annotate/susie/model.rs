//! Positive-effect SuSiE on one cluster's pseudobulk counts.
//!
//! ```text
//! y_g ~ NB(μ_g, φ_g)          the cluster's summed count of marker gene g
//! log μ_g = log μ⁰_g + Σ_c A_gc θ_c
//! θ = Σ_l e_l · 1[γ_l]        e_l = exp(b_l) > 0
//! γ_l ∈ {1..C, null}          one cell type per single effect, or none
//! ```
//!
//! `μ⁰_g` is the gene's background in this cluster (library size × background
//! rate). Each sweep draws every effect's type γ_l exactly from its full
//! conditional (C + 1 choices, each a change to the likelihood over that
//! type's markers only), then its size by elliptical slice sampling on the
//! standard-normal raw `b_l = m_b + s_b·z`: a log-normal size, the inverse
//! of the log link. Effects can only add expression, so a type cannot be
//! credited for its markers being low, and an effect that finds nothing rests
//! on the null choice. The draws keep each sweep's conditional choice
//! probabilities, so a type's PIP is Rao-Blackwellised.

use legume_numeric::mcmc::engine::{elliptical_slice_step, McmcModel};
use nalgebra::DVector;
use rand::rngs::SmallRng;
use rand::{Rng, RngExt};
use rand_distr::{Distribution, StandardNormal};

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

/// What every cluster shares: which types each marker gene marks and the
/// genes' NB sizes.
pub(super) struct Panel {
    types_of: Vec<Vec<usize>>,
    genes_of: Vec<Vec<usize>>,
    /// r_g = 1/φ_g (Var = μ + μ²/r); infinite is Poisson.
    nb_size: Vec<f32>,
    ln_r: Vec<f32>,
}

impl Panel {
    /// The panel over `types_of` (marker gene → its types), with φ_g from
    /// `dispersion` (φ = 0 is Poisson).
    pub(super) fn new(types_of: Vec<Vec<usize>>, n_types: usize, dispersion: &[f32]) -> Self {
        let mut genes_of = vec![Vec::new(); n_types];
        for (g, types) in types_of.iter().enumerate() {
            for &c in types {
                genes_of[c].push(g);
            }
        }
        let nb_size: Vec<f32> = dispersion
            .iter()
            .map(|&d| if d > 0.0 { 1.0 / d } else { f32::INFINITY })
            .collect();
        let ln_r = nb_size.iter().map(|r| r.ln()).collect();
        Self {
            types_of,
            genes_of,
            nb_size,
            ln_r,
        }
    }

    pub(super) fn n_types(&self) -> usize {
        self.genes_of.len()
    }
}

/// One cluster's counts over the panel's marker genes.
pub(super) struct Cluster {
    /// Summed count per marker gene.
    pub(super) y: Vec<f32>,
    /// log μ⁰ per marker gene.
    pub(super) log_mu0: Vec<f32>,
}

/// One cluster on the panel: the likelihood the sampler targets.
#[derive(Clone, Copy)]
pub(super) struct ClusterCounts<'a> {
    pub(super) panel: &'a Panel,
    pub(super) cluster: &'a Cluster,
}

impl ClusterCounts<'_> {
    /// Gene g's log-likelihood at linear predictor η, up to constants.
    fn term(&self, g: usize, eta: f32) -> f32 {
        let y = self.cluster.y[g];
        let r = self.panel.nb_size[g];
        if r.is_infinite() {
            // Poisson: y log μ − μ
            return y * eta - eta.exp();
        }
        // y log μ − (y + r) log(r + μ)
        y * eta - (y + r) * log_add_exp(self.panel.ln_r[g], eta)
    }

    /// Every marker gene's η at θ.
    pub(super) fn etas(&self, theta: &[f32]) -> Vec<f32> {
        self.cluster
            .log_mu0
            .iter()
            .zip(&self.panel.types_of)
            .map(|(&m0, types)| m0 + types.iter().map(|&c| theta[c]).sum::<f32>())
            .collect()
    }

    fn log_lik(&self, theta: &[f32]) -> f32 {
        self.etas(theta)
            .iter()
            .enumerate()
            .map(|(g, &eta)| self.term(g, eta))
            .sum()
    }

    /// Each type's share of the deviance θ explains over the background:
    /// the log-likelihood lost when its effect alone is removed, over the
    /// log-likelihood θ gains on θ = 0 (0 for types without an effect).
    pub(super) fn explained(&self, theta: &[f32]) -> Vec<f32> {
        let total = self.log_lik(theta) - self.log_lik(&vec![0.0; theta.len()]);
        let etas = self.etas(theta);
        (0..theta.len())
            .map(|c| {
                if theta[c] <= 0.0 || total <= 0.0 {
                    return 0.0;
                }
                let lost: f32 = self.panel.genes_of[c]
                    .iter()
                    .map(|&g| self.term(g, etas[g]) - self.term(g, etas[g] - theta[c]))
                    .sum();
                (lost / total).max(0.0)
            })
            .collect()
    }

    /// The change in log-likelihood from adding `effect` to type c's markers,
    /// at the predictors `etas` whose terms are `base`.
    fn gain(&self, etas: &[f32], base: &[f32], c: usize, effect: f32) -> f32 {
        self.panel.genes_of[c]
            .iter()
            .map(|&g| self.term(g, etas[g] + effect) - base[g])
            .sum()
    }
}

/// log(eᵃ + eᵇ).
fn log_add_exp(a: f32, b: f32) -> f32 {
    a.max(b) + (-(a - b).abs()).exp().ln_1p()
}

fn std_normal(rng: &mut impl Rng) -> f32 {
    let v: f64 = StandardNormal.sample(rng);
    v as f32
}

/// The sampler over one cluster.
pub(super) struct SusieModel<'a> {
    pub(super) data: ClusterCounts<'a>,
    pub(super) prior: SusiePrior,
}

/// Sampler state: per effect, its type (`C` = null), its standard-normal raw
/// size and that size, the last conditional choice probabilities; θ, and the
/// marker genes' η at θ with their likelihood terms, kept in step with it.
pub(super) struct SusieState {
    gamma: Vec<usize>,
    z_eff: Vec<DVector<f32>>,
    effect: Vec<f32>,
    probs: Vec<Vec<f32>>,
    theta: Vec<f32>,
    etas: Vec<f32>,
    base: Vec<f32>,
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
    fn effect(&self, z: &DVector<f32>) -> f32 {
        (self.prior.effect_mean + self.prior.effect_sd * z[0])
            .min(5.0)
            .exp()
    }

    /// Log prior of each choice: the null weight, the rest shared evenly.
    fn log_prior(&self) -> (f32, f32) {
        let w = self.prior.null_weight.clamp(1e-6, 1.0 - 1e-6);
        let c = self.data.panel.n_types() as f32;
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
        let z_eff: Vec<DVector<f32>> = (0..l)
            .map(|_| DVector::from_element(1, std_normal(rng)))
            .collect();
        let effect = z_eff.iter().map(|z| self.effect(z)).collect();
        let theta = vec![0.0; c];
        let etas = self.data.etas(&theta);
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
        for l in 0..self.prior.num_effects {
            if s.gamma[l] < c {
                self.shift(s, s.gamma[l], -s.effect[l]);
            }

            // The type, exactly from its full conditional.
            let effect = s.effect[l];
            let gains: Vec<f32> = (0..c)
                .map(|j| self.data.gain(&s.etas, &s.base, j, effect))
                .collect();
            let mut p: Vec<f32> = gains.iter().map(|&g| log_type + g).collect();
            p.push(log_null);
            let m = p.iter().fold(f32::NEG_INFINITY, |m, &v| m.max(v));
            p.iter_mut().for_each(|v| *v = (*v - m).exp());
            let total: f32 = p.iter().sum();
            p.iter_mut().for_each(|v| *v /= total);
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
            s.z_eff[l] = if j < c {
                let ll = |z: &DVector<f32>| self.data.gain(&s.etas, &s.base, j, self.effect(z));
                let nu = DVector::from_element(1, std_normal(rng));
                elliptical_slice_step(&s.z_eff[l], &nu, &ll, gains[j], rng).0
            } else {
                DVector::from_element(1, std_normal(rng))
            };
            s.effect[l] = self.effect(&s.z_eff[l]);
            if j < c {
                self.shift(s, j, s.effect[l]);
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
