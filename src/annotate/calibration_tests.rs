//! Calibration of the enrichment nulls and the marker bootstrap, on
//! synthetic data where the truth is known: with no signal the p-values must
//! not over-call at any batch count (the sample permutation is skipped when
//! there are too few batches, see [`super::by_enrichment::sample_perm_draws`]);
//! planted markers must be called; the bootstrap must split an ambiguous
//! cluster. Slow, so ignored: `cargo test --release -- --ignored calibration`.

use super::args::{BOOT_NUM_DRAWS, MIN_CONFIDENCE, NUM_DRAWS};
use super::by_enrichment::sample_perm_draws;
use enrichment::consensus::Abstain;
use enrichment::marker_bootstrap::EnrichmentBootstrapConfig;
use enrichment::{annotate, AnnotateConfig, GroupInputs, SpecificityMode};
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
    /// Per cluster: bootstrap support and credible-set size.
    boot: Option<Vec<(f32, usize)>>,
}

/// Score a synthetic run with lupin's null settings. Type `t` owns genes
/// `t * MARKERS_PER_TYPE ..`; each `(cluster, type, fold)` in `planted` raises
/// that type's markers `fold`-fold in that cluster. Cells spread over
/// `batches` batches at random, so each batch has its own cluster mix.
fn run(batches: usize, planted: &[(usize, usize, f32)], seed: u64, n_boot: usize) -> Summary {
    let mut rng = Rng(seed);
    let mut markers = Mat::zeros(GENES, TYPES);
    for t in 0..TYPES {
        for m in 0..MARKERS_PER_TYPE {
            markers[(t * MARKERS_PER_TYPE + m, t)] = 1.0;
        }
    }
    let base: Vec<f32> = (0..GENES).map(|_| 0.2 + 3.0 * rng.exponential()).collect();
    let mut profile = Mat::zeros(GENES, CLUSTERS);
    for (g, b) in base.iter().enumerate() {
        for k in 0..CLUSTERS {
            profile[(g, k)] = b * (0.6 + 0.8 * rng.uniform());
        }
    }
    for &(k, t, fold) in planted {
        for m in 0..MARKERS_PER_TYPE {
            profile[(t * MARKERS_PER_TYPE + m, k)] *= fold;
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
    let mut pseudobulk = Mat::zeros(GENES, batches);
    for (b, row) in counts.iter().enumerate() {
        let total: f32 = row.iter().sum::<f32>().max(1.0);
        for (k, c) in row.iter().enumerate() {
            membership[(b, k)] = c / total;
        }
        for g in 0..GENES {
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
        gene_names: (0..GENES).map(|g| format!("GENE{g}").into()).collect(),
        cell_names: (0..CELLS).map(|n| format!("cell{n}").into()).collect(),
    };
    let names: Vec<Box<str>> = (0..TYPES).map(|t| format!("CT{t}").into()).collect();
    let config = AnnotateConfig {
        specificity: SpecificityMode::Simplex,
        num_row_randomization: NUM_DRAWS,
        num_sample_perm: sample_perm_draws(batches, 200),
        batch_labels: None,
        fdr_alpha: 0.1,
        q_softmax_temperature: 1.0,
        min_confidence: MIN_CONFIDENCE,
        seed,
        min_markers: 3,
        stratify_null: true,
        bootstrap: (n_boot > 0).then_some(EnrichmentBootstrapConfig {
            n_boot,
            abstain: Abstain::Support(0.5),
            set_coverage: 0.8,
            max_set_size: 3,
            boot_num_draws: BOOT_NUM_DRAWS,
        }),
    };
    let out = annotate(&group, &markers, &names, &config).unwrap();

    let is_planted = |k: usize, t: usize| planted.iter().any(|&(pk, pt, _)| pk == k && pt == t);
    let mut null_p = Vec::new();
    let mut null_calls = 0;
    let mut planted_q = Vec::new();
    for k in 0..CLUSTERS {
        for t in 0..TYPES {
            let q = out.qvalue_kc[(k, t)];
            if is_planted(k, t) {
                planted_q.push(q);
            } else {
                null_p.push(out.pvalue_kc[(k, t)]);
                null_calls += usize::from(q < 0.1);
            }
        }
    }
    let null_frac05 = null_p.iter().filter(|&&p| p < 0.05).count() as f32 / null_p.len() as f32;
    let boot = out.bootstrap.map(|b| {
        (0..CLUSTERS)
            .map(|k| (b.consensus.support[k], b.consensus.label_set[k].len()))
            .collect()
    });
    Summary {
        null_frac05,
        null_calls,
        planted_q,
        boot,
    }
}

#[test]
#[ignore]
fn calibration_null_does_not_over_call_at_any_batch_count() {
    for batches in [1usize, 2, 6] {
        let s = run(batches, &[], 7, 0);
        assert!(
            s.null_frac05 <= 0.10,
            "{batches} batch(es): {:.3} of null p below 0.05",
            s.null_frac05
        );
        assert!(
            s.null_calls <= 2,
            "{batches} batch(es): {} null pairs called at q < 0.1",
            s.null_calls
        );
    }
}

#[test]
#[ignore]
fn calibration_planted_markers_are_called_and_the_bootstrap_splits_ties() {
    let s = run(
        6,
        &[(0, 0, 3.0), (1, 1, 3.0), (2, 2, 2.0), (3, 3, 1.4)],
        11,
        50,
    );
    assert!(s.planted_q.iter().all(|&q| q < 0.05), "{:?}", s.planted_q);
    let boot = s.boot.unwrap();
    assert!(
        boot[..4].iter().all(|&(support, _)| support >= 0.9),
        "{boot:?}"
    );

    // Cluster 5 carries two types' markers equally.
    let s = run(6, &[(5, 4, 2.5), (5, 5, 2.5)], 13, 50);
    let (support, set) = s.boot.unwrap()[5];
    assert!(
        support < 0.8 && set == 2,
        "support {support}, set size {set}"
    );
}
