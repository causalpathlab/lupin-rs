//! Calibration of the enrichment nulls, on synthetic data where the truth is
//! known: with no signal the p-values must not over-call; planted markers
//! must be called, and graded by their strength rather than tied at a
//! floor. Scored as the pass scores ([`super::by_enrichment::null_config`]).
//! Slow, so ignored: `cargo test --release -- --ignored calibration`.

use super::by_enrichment::null_config;
use enrichment::{annotate, GroupInputs};
use legume_numeric::matrix::dense_mat_io::Mat;

const GENES: usize = 2000;
const CLUSTERS: usize = 10;
const TYPES: usize = 8;
const MARKERS_PER_TYPE: usize = 12;
const CELLS: usize = 3000;

struct Rng(u64);

impl Rng {
    fn uniform(&mut self) -> f32 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 40) as f32 + 0.5) / (1u64 << 24) as f32
    }

    fn exponential(&mut self) -> f32 {
        -self.uniform().ln()
    }
}

struct Summary {
    /// Share of null (cluster, type) p-values below 0.05.
    null_frac05: f32,
    /// Null pairs called at q < 0.1.
    null_calls: usize,
    /// q of each planted (cluster, type) pair.
    planted_q: Vec<f32>,
    /// p of each planted (cluster, type) pair.
    planted_p: Vec<f32>,
}

/// Score a synthetic run with lupin's null settings. Type `t` owns genes
/// `t * MARKERS_PER_TYPE ..`; each `(cluster, type, fold)` in `planted` raises
/// that type's markers `fold`-fold in that cluster. Cells spread over
/// `batches` batches at random, so each batch has its own cluster mix.
fn run(batches: usize, planted: &[(usize, usize, f32)], seed: u64) -> Summary {
    run_rows(batches, planted, seed, 1, false)
}

/// [`run`] on a matrix where each gene has `per_gene` rows that follow its
/// clusters alike (a faba matrix's channels and sites: the gene's pattern at
/// each row's own level, with a little noise). A type's markers are marked
/// on every row of the gene (`every_row`), else on its first alone.
fn run_rows(
    batches: usize,
    planted: &[(usize, usize, f32)],
    seed: u64,
    per_gene: usize,
    every_row: bool,
) -> Summary {
    let rows = GENES * per_gene;
    let mut rng = Rng(seed);
    let mut markers = Mat::zeros(rows, TYPES);
    for t in 0..TYPES {
        for m in 0..MARKERS_PER_TYPE {
            let gene = t * MARKERS_PER_TYPE + m;
            for r in 0..if every_row { per_gene } else { 1 } {
                markers[(gene * per_gene + r, t)] = 1.0;
            }
        }
    }
    let base: Vec<f32> = (0..GENES).map(|_| 0.2 + 3.0 * rng.exponential()).collect();
    let mut profile = Mat::zeros(rows, CLUSTERS);
    for (g, b) in base.iter().enumerate() {
        let pattern: Vec<f32> = (0..CLUSTERS).map(|_| 0.6 + 0.8 * rng.uniform()).collect();
        for r in 0..per_gene {
            // Each further row at a lower level: unspliced, then sites.
            let level = b / (1 + 2 * r) as f32;
            for (k, p) in pattern.iter().enumerate() {
                let noise = if r == 0 {
                    1.0
                } else {
                    0.95 + 0.1 * rng.uniform()
                };
                profile[(g * per_gene + r, k)] = level * p * noise;
            }
        }
    }
    for &(k, t, fold) in planted {
        for m in 0..MARKERS_PER_TYPE {
            for r in 0..per_gene {
                profile[((t * MARKERS_PER_TYPE + m) * per_gene + r, k)] *= fold;
            }
        }
    }

    let mut cells = Mat::zeros(CELLS, CLUSTERS);
    let mut counts = vec![vec![0f32; CLUSTERS]; batches];
    for n in 0..CELLS {
        let k = n % CLUSTERS;
        cells[(n, k)] = 1.0;
        let b = ((rng.uniform() * batches as f32) as usize).min(batches - 1);
        counts[b][k] += 1.0;
    }
    let mut membership = Mat::zeros(batches, CLUSTERS);
    let mut pseudobulk = Mat::zeros(rows, batches);
    for (b, row) in counts.iter().enumerate() {
        let total: f32 = row.iter().sum::<f32>().max(1.0);
        for (k, c) in row.iter().enumerate() {
            membership[(b, k)] = c / total;
        }
        for g in 0..rows {
            let mix: f32 = (0..CLUSTERS)
                .map(|k| membership[(b, k)] * profile[(g, k)])
                .sum();
            pseudobulk[(g, b)] = mix * (0.8 + 0.4 * rng.uniform());
        }
    }

    let group = GroupInputs {
        profile_gk: profile,
        pb_gene_gp: pseudobulk,
        pb_membership_pk: membership,
        cell_membership_nk: cells,
        gene_names: (0..rows).map(|g| format!("GENE{g}").into()).collect(),
        cell_names: (0..CELLS).map(|n| format!("cell{n}").into()).collect(),
    };
    let names: Vec<Box<str>> = (0..TYPES).map(|t| format!("CT{t}").into()).collect();
    let config = null_config(seed, 0.1, 1.0, 3, None);
    let out = annotate(&group, &markers, &names, &config).unwrap();

    let is_planted = |k: usize, t: usize| planted.iter().any(|&(pk, pt, _)| pk == k && pt == t);
    let mut null_p = Vec::new();
    let mut null_calls = 0;
    let (mut planted_q, mut planted_p) = (Vec::new(), Vec::new());
    for k in 0..CLUSTERS {
        for t in 0..TYPES {
            let q = out.qvalue_kc[(k, t)];
            if is_planted(k, t) {
                planted_q.push(q);
                planted_p.push(out.pvalue_kc[(k, t)]);
            } else {
                null_p.push(out.pvalue_kc[(k, t)]);
                null_calls += usize::from(q < 0.1);
            }
        }
    }
    let null_frac05 = null_p.iter().filter(|&&p| p < 0.05).count() as f32 / null_p.len() as f32;
    Summary {
        null_frac05,
        null_calls,
        planted_q,
        planted_p,
    }
}

#[test]
#[ignore]
fn calibration_null_does_not_over_call() {
    let s = run(6, &[], 7);
    assert!(
        s.null_frac05 <= 0.10,
        "{:.3} of null p below 0.05",
        s.null_frac05
    );
    assert!(
        s.null_calls <= 2,
        "{} null pairs called at q < 0.1",
        s.null_calls
    );
}

#[test]
#[ignore]
fn calibration_planted_markers_are_called() {
    let s = run(6, &[(0, 0, 3.0), (1, 1, 3.0), (2, 2, 2.0), (3, 3, 1.4)], 11);
    assert!(s.planted_q.iter().all(|&q| q < 0.05), "{:?}", s.planted_q);
}

#[test]
#[ignore]
fn calibration_a_gene_of_several_rows_is_marked_once() {
    // No signal; each gene has four rows that move together. Marked on
    // one row a gene, the null holds; marked on every row, a type's markers
    // are a quarter as many genes as the null counts them for.
    let once = run_rows(2, &[], 7, 4, false);
    assert!(
        once.null_frac05 <= 0.10,
        "one row a gene: {:.3} of null p below 0.05",
        once.null_frac05
    );
    let every = run_rows(2, &[], 7, 4, true);
    assert!(
        every.null_frac05 > 0.10,
        "every row: {:.3} of null p below 0.05 (expected to over-call)",
        every.null_frac05
    );
}

#[test]
#[ignore]
fn calibration_strong_types_are_graded_not_tied() {
    // Two types up in one cluster, one far more than the other: beyond what
    // the draws alone resolve, the stronger still tests smaller.
    let s = run(2, &[(0, 0, 4.0), (0, 1, 12.0)], 5);
    let (weak, strong) = (s.planted_p[0], s.planted_p[1]);
    assert!(strong < weak, "p {strong:e} against {weak:e}");
    let floor = 1.0 / (super::args::NUM_DRAWS as f32 + 1.0);
    assert!(weak < floor, "both past the draws' floor: {weak:e}");
}
