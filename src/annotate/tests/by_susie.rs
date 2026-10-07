//! The SuSiE stage's outputs for [`super`].

use super::*;
use crate::annotate::outputs::{
    CLUSTER_CELLTYPE_EFFECT, CLUSTER_CELLTYPE_EXPLAINED, CLUSTER_CELLTYPE_PIP,
};

/// Six clusters of two cells over three types of 20 markers and 200 genes no
/// type claims; cluster k's type k % 3 markers at five times background.
fn inputs() -> EnrichmentInputs {
    let (n_free, n_types, per_type, n_clusters) = (200, 3, 20, 6);
    let n_genes = n_free + n_types * per_type;
    let mut markers_gc = Mat::zeros(n_genes, n_types);
    for c in 0..n_types {
        for g in 0..per_type {
            markers_gc[(n_free + c * per_type + g, c)] = 1.0;
        }
    }
    let mut gene_sum_kg = Vec::with_capacity(n_genes * n_clusters);
    for k in 0..n_clusters {
        for g in 0..n_genes {
            let base = 50.0 + (g % 7) as f64 * 10.0;
            let up = g >= n_free && (g - n_free) / per_type == k % n_types;
            gene_sum_kg.push(if up { 5.0 * base } else { base });
        }
    }
    let names =
        |p: &str, n: usize| -> Vec<Box<str>> { (0..n).map(|i| format!("{p}{i}").into()).collect() };
    EnrichmentInputs {
        gene_names: names("G", n_genes),
        cell_names: names("c", 2 * n_clusters),
        cluster_labels: (0..2 * n_clusters).map(|i| i / 2).collect(),
        n_clusters,
        batch_labels: vec![0; 2 * n_clusters],
        n_batches: 1,
        markers_gc,
        celltype_names: names("T", n_types),
        profile_gk: Mat::zeros(n_genes, n_clusters),
        pb_gene_gp: Mat::zeros(n_genes, 1),
        gene_sum_kg,
        gene_weights: vec![1.0; n_genes],
        type_tree: None,
        cl_record: None,
        expression_source: None,
    }
}

#[test]
fn the_stage_reports_its_tables_under_their_own_keys_and_leaves_q_alone() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("run").to_string_lossy().into_owned();
    let cfg = SusieConfig {
        samples: 300,
        warmup: 200,
        ..Default::default()
    };
    let mut outputs = AnnotationOutputs {
        cluster_celltype_q: Some("the enrichment's softmax Q".into()),
        ..Default::default()
    };
    run(&out, &cfg, &inputs(), &mut outputs).unwrap();

    assert_eq!(
        outputs.cluster_celltype_q.as_deref(),
        Some("the enrichment's softmax Q"),
        "Q stays the enrichment's"
    );
    for (path, suffix) in [
        (&outputs.cluster_celltype_pip, CLUSTER_CELLTYPE_PIP),
        (&outputs.cluster_celltype_effect, CLUSTER_CELLTYPE_EFFECT),
        (
            &outputs.cluster_celltype_explained,
            CLUSTER_CELLTYPE_EXPLAINED,
        ),
    ] {
        assert_eq!(path.as_deref(), Some(format!("{out}{suffix}").as_str()));
        assert!(std::path::Path::new(path.as_deref().unwrap()).exists());
    }
    let calls = crate::manifest::rounds::read_argmax(outputs.argmax.as_deref().unwrap()).unwrap();
    for i in 0..12 {
        assert_eq!(
            calls[&*format!("c{i}")].0,
            format!("T{}", (i / 2) % 3),
            "cell c{i}"
        );
    }
}
