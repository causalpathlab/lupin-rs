//! A fully Bayesian SuSiE over all clusters at once: the background and the
//! genes' misfit become parameters of one model instead of plug-in estimates.
//!
//! ```text
//! u_kg ~ N(β_g + η_kg, (σ²_g + δ_kg)/λ_kg)    η_k = Σ_c A_gc θ_kc (SuSiE, per cluster)
//! β_g  ~ N(m, τ²)                              the gene's baseline: its background
//! λ_kg ~ Gamma(ν/2, ν/2)                       integrated out, the Student-t
//! log σ²_g ~ N(μ_s, ω²),  μ_s ~ N(m_s, v_s),  ω² ~ InvGamma(a_ω, b_ω)
//! ```
//!
//! `u` is the log count over the library and `δ` its sampling variance,
//! 1/(y + ½). With β and σ² fixed this is exactly the plug-in model of
//! [`super::model::ClusterLogFold`], since a t is a Gaussian scale mixture over
//! a Gamma precision; freeing them is what the comparison with that model
//! measures.
//!
//! One sweep: every cluster's SuSiE step in parallel ([`super::model`],
//! Gaussian at the current β, σ², λ), then λ (conjugate Gamma), β
//! (conjugate Normal), log σ² (elliptical slice on its standardized value)
//! and the hyper-parameters (conjugate).

use super::model::{ClusterGauss, ClusterT, Panel, SusieModel, SusieSample, SusieState, NU};
use super::{summarize, ClusterFit, SusieConfig};
use legume_numeric::matrix::rand_util::mix_seed;
use legume_numeric::matrix::utils::median;
use legume_numeric::mcmc::engine::{
    elliptical_slice_step, run_mcmc_parallel, McmcConfig, McmcModel,
};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};
use rand_distr::{Distribution, Gamma, StandardNormal};
use rayon::prelude::*;

/// Prior variance of a gene's baseline β around the data's median log rate
/// (on the log scale: flat over any plausible background).
const BETA_VAR: f64 = 100.0;
/// Prior centre and variance of μ_s, the genes' mean log misfit variance.
const LOG_S2_MEAN: f64 = -2.3;
const LOG_S2_VAR: f64 = 9.0;
/// InvGamma(shape, scale) prior on ω², the spread of the genes' log σ².
const OMEGA_SHAPE: f64 = 2.0;
const OMEGA_SCALE: f64 = 1.0;

/// Parameters held at given values instead of sampled (the comparison's
/// stages).
pub(super) struct Fixed {
    pub(super) beta: Option<Vec<f32>>,
    pub(super) scale2: Option<Vec<f32>>,
}

/// Where the chains start: the plug-in's background (as a log baseline) and
/// misfit scales. Each chain starts at them plus a small seeded jitter, so the
/// chains still begin apart and R̂ can tell whether they converged.
pub(super) struct Start {
    pub(super) beta: Vec<f32>,
    pub(super) scale2: Vec<f32>,
}

/// Jitter of a chain's starting baseline and log misfit variance.
const START_JITTER: (f64, f64) = (0.1, 0.2);

/// Every cluster's posterior plus the genes' posterior means.
pub(super) struct JointFits {
    pub(super) clusters: Vec<ClusterFit>,
    /// Posterior mean baseline β_g per marker gene.
    pub(super) beta: Vec<f32>,
    /// Posterior mean σ²_g per marker gene.
    pub(super) scale2: Vec<f32>,
}

struct JointModel<'a> {
    panel: &'a Panel,
    u: &'a [Vec<f32>],
    delta: &'a [Vec<f32>],
    fixed: &'a Fixed,
    start: Option<&'a Start>,
    prior: super::SusiePrior,
    beta_mean: f64,
}

struct JointState {
    clusters: Vec<SusieState>,
    beta: Vec<f32>,
    /// log σ²_g = μ_s + ω·z_g.
    z: Vec<f32>,
    mu_s: f64,
    omega: f64,
    lambda: Vec<Vec<f32>>,
}

#[derive(Clone)]
struct JointSample {
    clusters: Vec<SusieSample>,
    beta: Vec<f32>,
    scale2: Vec<f32>,
}

fn normal(rng: &mut SmallRng) -> f64 {
    StandardNormal.sample(rng)
}

impl JointModel<'_> {
    fn n_genes(&self) -> usize {
        self.u.first().map_or(0, Vec::len)
    }

    fn scale2(&self, s: &JointState) -> Vec<f64> {
        match &self.fixed.scale2 {
            Some(v) => v.iter().map(|&a| f64::from(a)).collect(),
            None => {
                s.z.iter()
                    .map(|&z| (s.mu_s + s.omega * f64::from(z)).exp())
                    .collect()
            }
        }
    }

    /// Cluster k's Gaussian likelihood at the current β, σ², λ.
    fn gauss_inputs(
        &self,
        k: usize,
        beta: &[f32],
        s2: &[f64],
        lambda: &[f32],
    ) -> (Vec<f32>, Vec<f64>) {
        let x = self.u[k].iter().zip(beta).map(|(&u, &b)| u - b).collect();
        let w = (0..self.n_genes())
            .map(|g| (s2[g] + f64::from(self.delta[k][g])) / f64::from(lambda[g]))
            .collect();
        (x, w)
    }
}

impl JointModel<'_> {
    /// Every precision λ_kg from its conjugate Gamma at the current β, σ², η.
    fn draw_lambda(&self, s: &mut JointState, s2: &[f64], etas: &[Vec<f32>], rng: &mut SmallRng) {
        let shape = (NU + 1.0) / 2.0;
        for (k, eta) in etas.iter().enumerate() {
            for (g, &e) in eta.iter().enumerate() {
                let r = f64::from(self.u[k][g] - s.beta[g] - e);
                let rate = (NU + r * r / (s2[g] + f64::from(self.delta[k][g]))) / 2.0;
                s.lambda[k][g] = Gamma::new(shape, 1.0 / rate).unwrap().sample(rng) as f32;
            }
        }
    }
}

impl McmcModel for JointModel<'_> {
    type State = JointState;
    type Sample = JointSample;
    type Result = Vec<JointSample>;

    fn init(&self, rng: &mut SmallRng) -> JointState {
        let m = self.n_genes();
        let beta = match (&self.fixed.beta, self.start) {
            (Some(b), _) => b.clone(),
            (None, Some(st)) => st
                .beta
                .iter()
                .map(|&b| (f64::from(b) + START_JITTER.0 * normal(rng)) as f32)
                .collect(),
            (None, None) => (0..m)
                .map(|g| median(&self.u.iter().map(|u| u[g]).collect::<Vec<f32>>()))
                .collect(),
        };
        // log σ² = μ_s + ω·z: from the start's scales (jittered) when given.
        let (mu_s, omega, z) = match self.start {
            Some(st) if self.fixed.scale2.is_none() => {
                let logs: Vec<f64> = st
                    .scale2
                    .iter()
                    .map(|&v| f64::from(v).ln() + START_JITTER.1 * normal(rng))
                    .collect();
                let n = logs.len().max(1) as f64;
                let mu = logs.iter().sum::<f64>() / n;
                let sd = (logs.iter().map(|l| (l - mu).powi(2)).sum::<f64>() / n)
                    .sqrt()
                    .max(0.1);
                let z = logs.iter().map(|l| ((l - mu) / sd) as f32).collect();
                (mu, sd, z)
            }
            _ => (LOG_S2_MEAN, 1.0, vec![0.0; m]),
        };
        let lambda = vec![vec![1.0f32; m]; self.u.len()];
        let s2: Vec<f64> = match &self.fixed.scale2 {
            Some(v) => v.iter().map(|&a| f64::from(a)).collect(),
            None => z
                .iter()
                .map(|&z: &f32| (mu_s + omega * f64::from(z)).exp())
                .collect(),
        };
        let clusters = (0..self.u.len())
            .map(|k| {
                let (x, w) = self.gauss_inputs(k, &beta, &s2, &lambda[k]);
                let model = SusieModel {
                    data: ClusterGauss {
                        panel: self.panel,
                        x: &x,
                        w: &w,
                    },
                    prior: self.prior,
                };
                model.init(rng)
            })
            .collect();
        JointState {
            clusters,
            beta,
            z,
            mu_s,
            omega,
            lambda,
        }
    }

    fn sweep(&self, s: &mut JointState, rng: &mut SmallRng) {
        let m = self.n_genes();
        let s2 = self.scale2(s);

        // Every cluster's SuSiE step, in parallel, each on its own stream,
        // under the Student-t: λ integrated out, so the effects are not held
        // to the last draw of the precisions (that data augmentation mixes
        // slowly where a few genes are outlying). λ is drawn right after.
        let seeds: Vec<u64> = (0..self.u.len()).map(|_| rng.random()).collect();
        let JointState { clusters, beta, .. } = s;
        clusters
            .par_iter_mut()
            .enumerate()
            .zip(seeds)
            .for_each(|((k, st), seed)| {
                let x: Vec<f32> = self.u[k]
                    .iter()
                    .zip(beta.iter())
                    .map(|(&u, &b)| u - b)
                    .collect();
                let v: Vec<f64> = (0..m)
                    .map(|g| s2[g] + f64::from(self.delta[k][g]))
                    .collect();
                let model = SusieModel {
                    data: ClusterT {
                        panel: self.panel,
                        x: &x,
                        v: &v,
                    },
                    prior: self.prior,
                };
                model.sweep(st, &mut SmallRng::seed_from_u64(seed));
            });
        let etas: Vec<Vec<f32>> = s
            .clusters
            .iter()
            .map(|st| self.panel.etas(st.theta()))
            .collect();

        self.draw_lambda(s, &s2, &etas, rng);

        // β: conjugate Normal given λ, pooled over the clusters.
        if self.fixed.beta.is_none() {
            for g in 0..m {
                let (mut prec, mut num) = (1.0 / BETA_VAR, self.beta_mean / BETA_VAR);
                for (k, eta) in etas.iter().enumerate() {
                    let p = f64::from(s.lambda[k][g]) / (s2[g] + f64::from(self.delta[k][g]));
                    prec += p;
                    num += p * f64::from(self.u[k][g] - eta[g]);
                }
                s.beta[g] = (num / prec + normal(rng) / prec.sqrt()) as f32;
            }
        }

        // log σ²: elliptical slice on the standardized z_g under the Student-t
        // with λ integrated out (σ² and the gene's λs otherwise trade off, a
        // ridge the chain crawls along), then the hyper-parameters, conjugate
        // given log σ²; λ is drawn after, given the new σ².
        if self.fixed.scale2.is_none() {
            for g in 0..m {
                let terms = |z: f32| -> f64 {
                    let v = (s.mu_s + s.omega * f64::from(z)).exp();
                    etas.iter()
                        .enumerate()
                        .map(|(k, eta)| {
                            let total = v + f64::from(self.delta[k][g]);
                            let r = f64::from(self.u[k][g] - s.beta[g] - eta[g]);
                            -0.5 * total.ln() - 0.5 * (NU + 1.0) * (r * r / (NU * total)).ln_1p()
                        })
                        .sum()
                };
                let cur = terms(s.z[g]);
                let ll = |z: &f32| (terms(*z) - cur) as f32;
                let nu: f32 = normal(rng) as f32;
                s.z[g] = elliptical_slice_step(&s.z[g], &nu, &ll, 0.0, rng).0;
            }
            let logs: Vec<f64> =
                s.z.iter()
                    .map(|&z| s.mu_s + s.omega * f64::from(z))
                    .collect();
            let n = logs.len() as f64;
            let omega2 = s.omega * s.omega;
            let prec = 1.0 / LOG_S2_VAR + n / omega2;
            let mean = (LOG_S2_MEAN / LOG_S2_VAR + logs.iter().sum::<f64>() / omega2) / prec;
            s.mu_s = mean + normal(rng) / prec.sqrt();
            let ss: f64 = logs.iter().map(|l| (l - s.mu_s).powi(2)).sum();
            let draw = Gamma::new(OMEGA_SHAPE + n / 2.0, 1.0 / (OMEGA_SCALE + ss / 2.0))
                .unwrap()
                .sample(rng);
            s.omega = (1.0 / draw).sqrt();
            for (z, l) in s.z.iter_mut().zip(&logs) {
                *z = ((l - s.mu_s) / s.omega) as f32;
            }
        }

        // λ again, given the σ² just drawn with λ integrated out.
        let s2 = self.scale2(s);
        self.draw_lambda(s, &s2, &etas, rng);
    }

    fn collect(&self, s: &JointState) -> JointSample {
        JointSample {
            clusters: s.clusters.iter().map(SusieState::sample).collect(),
            beta: s.beta.clone(),
            scale2: self.scale2(s).iter().map(|&v| v as f32).collect(),
        }
    }

    fn summarize(&self, samples: Vec<JointSample>) -> Vec<JointSample> {
        samples
    }
}

/// Fit every cluster jointly. `u` and `delta` are clusters × marker genes:
/// each marker's log count over the cluster's library, and its sampling
/// variance.
pub(super) fn fit_joint(
    panel: &Panel,
    u: &[Vec<f32>],
    delta: &[Vec<f32>],
    fixed: &Fixed,
    start: Option<&Start>,
    cfg: &SusieConfig,
    corr: &[Vec<f32>],
) -> JointFits {
    let all: Vec<f32> = u.iter().flatten().copied().collect();
    let model = JointModel {
        panel,
        u,
        delta,
        fixed,
        start,
        prior: cfg.prior,
        beta_mean: f64::from(median(&all)),
    };
    let config = McmcConfig {
        n_samples: cfg.samples,
        warmup: cfg.warmup,
        thin: cfg.thin,
        seed: mix_seed(cfg.seed, 0),
    };
    let chains = run_mcmc_parallel(&model, &config, cfg.chains);
    let draws = chains.iter().map(Vec::len).sum::<usize>().max(1) as f64;
    let mean_of = |f: &dyn Fn(&JointSample) -> &[f32]| -> Vec<f32> {
        let m = model.n_genes();
        let mut acc = vec![0.0f64; m];
        for d in chains.iter().flatten() {
            for (a, &v) in acc.iter_mut().zip(f(d)) {
                *a += f64::from(v);
            }
        }
        acc.iter().map(|&a| (a / draws) as f32).collect()
    };
    let beta = mean_of(&|d| &d.beta);
    let scale2 = mean_of(&|d| &d.scale2);
    let clusters = (0..u.len())
        .map(|k| {
            let per_chain: Vec<Vec<SusieSample>> = chains
                .iter()
                .map(|c| c.iter().map(|d| d.clusters[k].clone()).collect())
                .collect();
            // The shares explained, at the posterior means.
            let x: Vec<f32> = u[k].iter().zip(&beta).map(|(&a, &b)| a - b).collect();
            let w: Vec<f64> = (0..x.len())
                .map(|g| f64::from(scale2[g]) + f64::from(delta[k][g]))
                .collect();
            let lik = ClusterGauss {
                panel,
                x: &x,
                w: &w,
            };
            summarize(&per_chain, lik, corr)
        })
        .collect();
    JointFits {
        clusters,
        beta,
        scale2,
    }
}
