//! The GO/GMT terms' test: each term's effect, as
//! [`enrichment::ontology_module_score`] reports it, against the effects of
//! random gene sets of the same size and expression level drawn from the
//! annotated genes. The test asks of a term what its effect says of it, so
//! the two agree in sign and order; a term's few cluster-specific genes do
//! not outweigh the rest.
//!
//! The null holds the term's genes per expression bin, so a set of highly
//! expressed genes (a ribosome) meets highly expressed decoys. Random genes
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

use super::geneset_units::{bh_where_finite, Null};
use enrichment::Mat;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use rayon::prelude::*;

/// Random gene sets per term.
pub const TERM_NULL_DRAWS: usize = 200;
/// Expression bins the decoys keep a term's genes in.
const ABUNDANCE_BINS: usize = 20;
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
    draws: usize,
) -> anyhow::Result<TermTest> {
    let n_clusters = profile_gk.ncols();
    anyhow::ensure!(
        effect_kt.shape() == (n_clusters, terms.len()),
        "effects of {} × {}, for {n_clusters} clusters × {} terms",
        effect_kt.nrows(),
        effect_kt.ncols(),
        terms.len()
    );
    let bins = AbundanceBins::of(profile_gk, universe);
    let columns: Vec<(Vec<f32>, Vec<f32>)> = terms
        .par_iter()
        .enumerate()
        .map(|(t, (_, rows))| {
            let mut rng = SmallRng::seed_from_u64(SEED ^ (t as u64).wrapping_mul(0x9E37_79B9));
            let decoys: Vec<(Box<str>, Vec<usize>)> = (0..draws)
                .map(|_| (Box::from(""), bins.draw_like(rows, &mut rng)))
                .collect();
            let null = enrichment::ontology_module_score(profile_gk, &decoys, universe)?.effect_kt;
            let observed: Vec<f32> = effect_kt.column(t).iter().copied().collect();
            let decoy_spread: Vec<f32> = (0..draws)
                .map(|d| mad(&null.column(d).iter().copied().collect::<Vec<_>>()))
                .collect();
            let widen = (mad(&observed) / mean(&decoy_spread)).max(1.0);
            Ok((0..n_clusters)
                .map(|k| {
                    let scores: Vec<f32> = null.row(k).iter().copied().collect();
                    let z = Null::of(&scores).z(observed[k]) / widen;
                    (z, upper_tail(z))
                })
                .unzip())
        })
        .collect::<anyhow::Result<_>>()?;
    let mut z_kt = Mat::zeros(n_clusters, terms.len());
    let mut p_kt = Mat::zeros(n_clusters, terms.len());
    for (t, (z, p)) in columns.iter().enumerate() {
        for k in 0..n_clusters {
            (z_kt[(k, t)], p_kt[(k, t)]) = (z[k], p[k]);
        }
    }
    let mut q_kt = Mat::zeros(n_clusters, terms.len());
    for k in 0..n_clusters {
        let p: Vec<f32> = p_kt.row(k).iter().copied().collect();
        for (t, q) in bh_where_finite(&p).into_iter().enumerate() {
            q_kt[(k, t)] = q;
        }
    }
    Ok(TermTest { z_kt, p_kt, q_kt })
}

/// The median absolute deviation of `v`, scaled to a normal's sd.
fn mad(v: &[f32]) -> f32 {
    use legume_numeric::matrix::utils::median;
    let m = median(v);
    let dev: Vec<f32> = v.iter().map(|x| (x - m).abs()).collect();
    1.4826 * median(&dev)
}

fn mean(v: &[f32]) -> f32 {
    v.iter().sum::<f32>() / v.len().max(1) as f32
}

/// `P(Z ≥ z)`; NaN for a NaN `z`.
fn upper_tail(z: f32) -> f32 {
    use statrs::distribution::{ContinuousCDF, Normal};
    if z.is_nan() {
        return f32::NAN;
    }
    Normal::standard().sf(f64::from(z)) as f32
}

/// The annotated genes in [`ABUNDANCE_BINS`] equal bins by their mean
/// log expression over the clusters.
struct AbundanceBins {
    /// Each profile row's bin; `None` off the annotated genes.
    of_row: Vec<Option<usize>>,
    rows: Vec<Vec<usize>>,
}

impl AbundanceBins {
    fn of(profile_gk: &Mat, universe: &[usize]) -> Self {
        let mut genes: Vec<usize> = universe.to_vec();
        genes.sort_unstable();
        genes.dedup();
        let level = |g: usize| profile_gk.row(g).mean().ln_1p();
        genes.sort_by(|&a, &b| level(a).total_cmp(&level(b)));
        let mut of_row = vec![None; profile_gk.nrows()];
        let mut rows = vec![Vec::new(); ABUNDANCE_BINS];
        for (rank, &g) in genes.iter().enumerate() {
            let bin = rank * ABUNDANCE_BINS / genes.len();
            of_row[g] = Some(bin);
            rows[bin].push(g);
        }
        Self { of_row, rows }
    }

    /// Random rows with as many of each bin as `like` has.
    fn draw_like(&self, like: &[usize], rng: &mut SmallRng) -> Vec<usize> {
        let mut per_bin = [0usize; ABUNDANCE_BINS];
        for b in like
            .iter()
            .filter_map(|&g| self.of_row.get(g).copied().flatten())
        {
            per_bin[b] += 1;
        }
        per_bin
            .iter()
            .enumerate()
            .filter(|&(_, &n)| n > 0)
            .flat_map(|(b, &n)| {
                let pool = &self.rows[b];
                rand::seq::index::sample(rng, pool.len(), n.min(pool.len()))
                    .into_iter()
                    .map(|i| pool[i])
                    .collect::<Vec<_>>()
            })
            .collect()
    }
}

#[cfg(test)]
#[path = "tests/term_test.rs"]
mod tests;
