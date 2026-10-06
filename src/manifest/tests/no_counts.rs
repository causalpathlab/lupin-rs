//! A pass without the raw counts: an earlier pass's cache stands in exactly,
//! and the decoder's sums are taken from evenly spaced cells.

use super::*;
use crate::manifest::recalibrate::write_cache;
use crate::manifest::rounds::write_clusters;

#[test]
fn spaced_takes_at_most_the_cap_evenly_and_skips_unassigned_cells() {
    let mut labels = vec![0usize; 3 * DECODED_CELLS];
    labels.extend([1, usize::MAX, 1]);
    let picked = spaced(&labels, 2);
    assert_eq!(picked[0].0.len(), DECODED_CELLS);
    assert_eq!(picked[0].1, 3 * DECODED_CELLS, "the group's whole size");
    assert_eq!(picked[0].0[1] - picked[0].0[0], 3, "every third cell");
    assert_eq!(
        picked[1],
        (vec![3 * DECODED_CELLS, 3 * DECODED_CELLS + 2], 2)
    );
}

#[test]
fn a_modules_sum_goes_to_its_genes_by_their_shares() {
    // Two groups × two modules; genes 0, 1 in module 0, gene 2 in module 1.
    let sum = [10.0, 4.0, 20.0, 8.0];
    let gene_of = [(0, 0.25), (0, 0.75), (1, 1.0)];
    assert_eq!(
        to_genes(&sum, 2, &gene_of),
        vec![2.5, 7.5, 4.0, 5.0, 15.0, 8.0]
    );
}

#[test]
fn the_mixture_decoder_sums_beta_times_theta() {
    // Two genes × two topics; one cluster whose θ sums to (1, 3).
    let beta = Mat::from_row_slice(2, 2, &[0.9, 0.2, 0.1, 0.8]);
    let theta = Mat::from_row_slice(2, 1, &[1.0, 3.0]);
    let e = expected(&beta, &theta);
    let want = [0.9 + 0.6, 0.1 + 2.4].map(|v| v * NOMINAL_DEPTH);
    for (a, b) in e.iter().zip(want) {
        assert!((a - b).abs() < 1e-3, "{e:?}");
    }
}

/// A `vae` run whose counts are elsewhere, and a round beside it whose pass
/// cached its statistics over clusters 3 and 7.
fn run_with_cached_round(dir: &Path) -> (Loaded, Vec<f64>) {
    let at = |n: &str| dir.join(n).to_string_lossy().into_owned();
    std::fs::write(
        at("run.senna.json"),
        r#"{"version": 2, "kind": "vae", "prefix": "run",
            "data": {"input": ["elsewhere/a.zarr.zip"]}}"#,
    )
    .unwrap();
    let cells: Vec<Box<str>> = ["c0", "c1", "c2", "c3"].map(Into::into).to_vec();
    let genes: Vec<Box<str>> = ["G0", "G1"].map(Into::into).to_vec();
    // Cluster slots 0, 1 (ids 3, 7), gene-major per slot.
    let gene_sum_kg = vec![5.0, 1.0, 2.0, 8.0];
    let inputs = EnrichmentInputs {
        gene_names: genes.clone(),
        cell_names: cells.clone(),
        cluster_labels: vec![0, 0, 1, usize::MAX],
        n_clusters: 2,
        batch_labels: vec![0, 1, 1, 0],
        n_batches: 2,
        markers_gc: Mat::zeros(2, 0),
        celltype_names: Vec::new(),
        profile_gk: Mat::zeros(2, 2),
        pb_gene_gp: Mat::from_element(2, 2, 0.5),
        gene_sum_kg: gene_sum_kg.clone(),
        gene_weights: vec![1.0, 2.0],
        type_tree: None,
        cl_record: None,
        expression_source: None,
    };
    let cache = write_cache(&at("run.L1"), &inputs).unwrap();
    // The cache's sums are by slot (K0, K1); a pass writes them by cluster id.
    let mut m = Mat::zeros(2, 2);
    for c in 0..2 {
        for g in 0..2 {
            m[(g, c)] = gene_sum_kg[c * 2 + g] as f32;
        }
    }
    m.to_parquet_with_names(
        &at("run.L1.cluster_gene_sum.parquet"),
        (Some(&genes[..]), Some("gene")),
        Some(&["K3".into(), "K7".into()]),
    )
    .unwrap();
    write_clusters(
        &at("run.L1.clusters.parquet"),
        &cells,
        &[Some(3), Some(3), Some(7), None],
    )
    .unwrap();
    let rel = |p: &str| {
        Path::new(p)
            .file_name()
            .unwrap()
            .to_string_lossy()
            .into_owned()
    };
    let round = json!({
        "version": 2, "kind": "vae", "prefix": "run",
        "annotate": {
            "source": "run.senna.json",
            "expression_clusters": "run.L1.clusters.parquet",
            "stats_cache": {
                "gene_sum": rel(&cache.gene_sum),
                "batch_profile": rel(&cache.batch_profile),
                "gene_weight": rel(&cache.gene_weight),
                "cell_batch": rel(&cache.cell_batch),
            }
        }
    });
    std::fs::write(at("run.L1.senna.json"), round.to_string()).unwrap();
    (run::load(&at("run.senna.json")).unwrap(), gene_sum_kg)
}

#[test]
fn without_counts_a_cached_round_gives_its_clusters_and_sums() {
    let root = tempfile::tempdir().unwrap();
    let (loaded, gene_sum_kg) = run_with_cached_round(root.path());
    assert!(
        matches!(source(&loaded), Source::Cache(src) if src.file.ends_with("run.L1.senna.json")),
        "the round's cache stands in"
    );
    let args = crate::annotate_cmd::default_enrichment_args("x");
    let Source::Cache(src) = source(&loaded) else {
        panic!("a cache")
    };
    let e = from_cache(&args, &loaded, &src).unwrap();
    assert_eq!(e.n_clusters, 2);
    assert_eq!(e.cluster_labels, vec![0, 0, 1, usize::MAX]);
    assert_eq!(e.gene_sum_kg, gene_sum_kg, "ids 3, 7 land in slots 0, 1");
    assert_eq!(e.gene_weights, vec![1.0, 2.0]);
    assert_eq!(e.batch_labels, vec![0, 1, 1, 0]);
    assert_eq!(e.n_batches, 2);
}

#[test]
fn without_counts_or_cache_a_run_with_no_decoder_has_nothing() {
    let root = tempfile::tempdir().unwrap();
    let file = root.path().join("run.senna.json");
    std::fs::write(
        &file,
        r#"{"version": 2, "kind": "svd", "prefix": "run",
            "data": {"input": ["elsewhere/a.zarr.zip"]}}"#,
    )
    .unwrap();
    let loaded = run::load(&file.to_string_lossy()).unwrap();
    assert!(matches!(source(&loaded), Source::Nothing));
}
