//! The GO/GMT terms' test: each term's effect, as
//! [`enrichment::ontology_module_score`] reports it, against the effects of
//! random gene sets of the same size and expression level drawn from the
//! annotated genes. The test asks of a term what its effect says of it, so
//! the two agree in sign and order; a term's few cluster-specific genes do
//! not outweigh the rest.
//!
//! The null holds the term's genes per expression stratum (the cell types'
//! [`GeneStrata`]), so a set of highly expressed genes (a ribosome) meets
//! highly expressed decoys. Random genes
//! do not move together as a term's do, so the decoys alone would make any
//! coherent programme significant (a third of all terms in a monocyte
//! cluster of a PBMC run). The decoys' spread is therefore widened by how
//! much more the term varies over the clusters than they do, both measured
//! by the median absolute deviation, which one cluster standing out does
//! not move: a programme high in many clusters (ribosome, spliceosome) is
//! widened, one high in a single cluster is not. The p is the normal tail
//! of that z, as the units' decoys give theirs ([`super::geneset_units`]).
//! A cell type split over many clusters shares its programmes among them,
//! so they test weaker in each than their effect reads.

use super::geneset_units::{bh_where_finite, mad_sd, normal_upper_tail, Null};
use enrichment::gene_strata::GeneStrata;
use enrichment::Mat;
use legume_numeric::matrix::hypothesis::mean;
use legume_numeric::matrix::utils::median;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use rayon::prelude::*;

/// Random gene sets per term.
pub const TERM_NULL_DRAWS: usize = 200;
/// The decoys' spread (log-scale effect) below which they are the term
/// itself drawn again, and test nothing.
const MIN_DECOY_SPREAD: f32 = 1e-6;
/// The decoys' seed, so a round's test can be made again.
const SEED: u64 = 0x5EED_7E57;

/// Each cluster × term: the effect's z against its decoys, its p and the
/// BH q over the cluster's terms.
pub struct TermTest {
    pub z_kt: Mat,
    pub p_kt: Mat,
    pub q_kt: Mat,
}

/// Test every term of `terms` (rows of `profile_gk`) on every cluster, its
/// observed effect being `effect_kt` (from the same `profile_gk`, `terms`
/// and `universe`).
pub fn test_terms(
    profile_gk: &Mat,
    terms: &[(Box<str>, Vec<usize>)],
    universe: &[usize],
    effect_kt: &Mat,
) -> anyhow::Result<TermTest> {
    let (n_clusters, n_terms) = (profile_gk.ncols(), terms.len());
    anyhow::ensure!(
        effect_kt.shape() == (n_clusters, n_terms),
        "effects of {} × {}, for {n_clusters} clusters × {n_terms} terms",
        effect_kt.nrows(),
        effect_kt.ncols(),
    );
    let bg = Background::of(profile_gk, universe);
    let columns: Vec<(Vec<f32>, Vec<f32>)> = terms
        .par_iter()
        .enumerate()
        .map(|(t, (_, rows))| {
            let mut rng = SmallRng::seed_from_u64(SEED ^ (t as u64).wrapping_mul(0x9E37_79B9));
            let mut scratch = bg.strata.scratch();
            let mut drawn = Vec::new();
            let like = bg.strata.profile_of(&bg.members(rows));
            let null: Vec<Vec<f32>> = (0..TERM_NULL_DRAWS)
                .map(|_| {
                    bg.strata
                        .draw_matched(&like, &mut scratch, &mut drawn, &mut rng);
                    bg.effect(drawn.iter().map(|&(i, _)| i as usize))
                })
                .collect();
            let observed: Vec<f32> = effect_kt.column(t).iter().copied().collect();
            let spread: Vec<f32> = null.iter().map(|d| mad_sd(d)).collect();
            let widen = (mad_sd(&observed) / mean(&spread)).max(1.0);
            Ok((0..n_clusters)
                .map(|k| {
                    let scores: Vec<f32> = null.iter().map(|d| d[k]).collect();
                    let z = if mad_sd(&scores) > MIN_DECOY_SPREAD {
                        Null::of(&scores).z(observed[k]) / widen
                    } else {
                        // The decoys are the term again (its strata hold
                        // little else): nothing to test it against.
                        f32::NAN
                    };
                    (z, normal_upper_tail(z))
                })
                .unzip())
        })
        .collect::<anyhow::Result<_>>()?;
    let z_kt = Mat::from_fn(n_clusters, n_terms, |k, t| columns[t].0[k]);
    let p_kt = Mat::from_fn(n_clusters, n_terms, |k, t| columns[t].1[k]);
    let q_rows: Vec<Vec<f32>> = (0..n_clusters)
        .map(|k| bh_where_finite(&p_kt.row(k).iter().copied().collect::<Vec<_>>()))
        .collect();
    let q_kt = Mat::from_fn(n_clusters, n_terms, |k, t| q_rows[k][t]);
    Ok(TermTest { z_kt, p_kt, q_kt })
}

/// The annotated genes as [`enrichment::ontology_module_score`] sees them:
/// each cluster's `log1p(CP10K)` over them and its total, and their
/// expression strata.
struct Background {
    /// Each profile row's position among the annotated genes.
    at: Vec<Option<u32>>,
    /// `log1p(CP10K)` per cluster, over the annotated genes.
    lge: Vec<Vec<f64>>,
    total: Vec<f64>,
    strata: GeneStrata,
}

impl Background {
    fn of(profile_gk: &Mat, universe: &[usize]) -> Self {
        let mut genes = universe.to_vec();
        genes.sort_unstable();
        genes.dedup();
        let mut at = vec![None; profile_gk.nrows()];
        for (i, &g) in genes.iter().enumerate() {
            at[g] = Some(i as u32);
        }
        let lge: Vec<Vec<f64>> = (0..profile_gk.ncols())
            .map(|k| {
                let raw: Vec<f64> = genes
                    .iter()
                    .map(|&g| f64::from(profile_gk[(g, k)]))
                    .collect();
                let scale = 1e4 / raw.iter().sum::<f64>().max(1e-12);
                raw.iter().map(|&v| (v * scale).ln_1p()).collect()
            })
            .collect();
        let total = lge.iter().map(|l| l.iter().sum()).collect();
        let abundance: Vec<f32> = genes
            .iter()
            .map(|&g| profile_gk.row(g).iter().map(|v| v.max(0.0)).sum())
            .collect();
        let strata = GeneStrata::by_covariate(&abundance);
        Self {
            at,
            lge,
            total,
            strata,
        }
    }

    /// A term's rows among the annotated genes, for [`GeneStrata::profile_of`].
    fn members(&self, rows: &[usize]) -> Vec<(u32, f32)> {
        let mut m: Vec<u32> = rows
            .iter()
            .filter_map(|&g| self.at.get(g).copied().flatten())
            .collect();
        m.sort_unstable();
        m.dedup();
        m.into_iter().map(|i| (i, 1.0)).collect()
    }

    /// The effect of the annotated genes `members` on each cluster, as
    /// [`enrichment::ontology_module_score`] computes it: their mean against
    /// the other annotated genes', less its median over the clusters.
    fn effect(&self, members: impl Iterator<Item = usize> + Clone) -> Vec<f32> {
        let n = self.lge[0].len() as f64;
        let m = members.clone().count() as f64;
        let s: Vec<f32> = (self.lge.iter().zip(&self.total))
            .map(|(lge, &total)| {
                if m < 1.0 || n - m < 1.0 {
                    return 0.0;
                }
                let sum_in: f64 = members.clone().map(|i| lge[i]).sum();
                (sum_in / m - (total - sum_in) / (n - m)) as f32
            })
            .collect();
        let centre = median(&s);
        s.iter().map(|v| v - centre).collect()
    }
}

#[cfg(test)]
#[path = "tests/term_test.rs"]
mod tests;
