//! `annotate --method projection` with the embedding marker score, end to end
//! on a planted toy run.

use super::*;
use crate::annotate::args::ProjectionMethod;
use crate::manifest::rounds::{read_argmax, write_clusters};
use crate::manifest::run;
use legume_numeric::matrix::traits::IoOps;

/// A bge run with ρ (g0..g2 on axis 0, g3..g5 on axis 1), twenty cells out on
/// the two axes, and the run's own clusters (0 for the first ten, 1 after).
fn toy_run(dir: &Path) -> (Loaded, Vec<Box<str>>) {
    let at = |n: &str| dir.join(n).to_string_lossy().into_owned();
    let genes: Vec<Box<str>> = (0..6).map(|g| format!("g{g}").into()).collect();
    let cells: Vec<Box<str>> = (0..20).map(|i| format!("c{i:02}").into()).collect();
    let mut rho = Mat::zeros(6, 2);
    for g in 0..6 {
        rho[(g, g / 3)] = 1.0 + 0.1 * g as f32;
        rho[(g, 1 - g / 3)] = 0.05 * g as f32;
    }
    let mut z = Mat::zeros(20, 2);
    for i in 0..20 {
        z[(i, i / 10)] = 2.0 + 0.1 * (i % 10) as f32;
        z[(i, 1 - i / 10)] = 0.2 * ((i % 3) as f32 - 1.0);
    }
    rho.to_parquet_with_names(
        &at("run.feature_embedding.parquet"),
        (Some(&genes[..]), Some("gene")),
        None,
    )
    .unwrap();
    z.to_parquet_with_names(
        &at("run.cell_embedding.parquet"),
        (Some(&cells[..]), Some("cell")),
        None,
    )
    .unwrap();
    let ids: Vec<Option<u32>> = (0..20).map(|i| Some(i / 10)).collect();
    write_clusters(&at("run.clusters.parquet"), &cells, &ids).unwrap();
    std::fs::write(
        at("panel.tsv"),
        "g0\tA\ng1\tA\ng2\tA\ng3\tB\ng4\tB\ng5\tB\n",
    )
    .unwrap();
    let manifest = serde_json::json!({
        "version": 2, "kind": "bge", "prefix": "run",
        "outputs": {"feature_embedding": "run.feature_embedding.parquet",
                    "cell_embedding": "run.cell_embedding.parquet"},
        "cluster": {"clusters": "run.clusters.parquet"}
    });
    std::fs::write(at("run.senna.json"), manifest.to_string()).unwrap();
    (run::load(&at("run.senna.json")).unwrap(), cells)
}

fn score_args(markers: &str, out: &str) -> AnnotateProjectionArgs {
    AnnotateProjectionArgs {
        projection: ProjectionMethod::Score,
        markers: markers.into(),
        out: out.into(),
        clusters: None,
        knn: 30,
        resolution: 1.0,
        num_clusters: None,
        min_cluster_size: 1,
        cluster_seed: None,
        num_perm: 0,
        seed: 42,
        min_markers: 3,
        no_idf: false,
        no_assign_qc: false,
        assign_mad: 2.5,
        fdr_alpha: 0.05,
        q_temperature: 1.0,
        obo: None,
        label_cl: None,
        ontology_fdr_q: 0.1,
        ontology_by: false,
        no_clean: false,
    }
}

#[test]
fn the_score_pass_calls_each_cell_and_records_its_round() {
    let dir = tempfile::tempdir().unwrap();
    let (loaded, cells) = toy_run(dir.path());
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    let panel = dir.path().join("panel.tsv").to_string_lossy().into_owned();
    annotate_by_projection(&score_args(&panel, &out), &loaded).unwrap();

    let argmax = read_argmax(&format!("{out}.argmax.tsv")).unwrap();
    for (i, c) in cells.iter().enumerate() {
        let (label, p) = &argmax[c.as_ref()];
        assert_eq!(label, if i < 10 { "A" } else { "B" }, "{c}");
        assert!(p.is_nan(), "the score has no probability");
    }
    let (_, nums) = read_table_columns(&format!("{out}{ANNOT_PARQUET}"), &[], &["margin"]).unwrap();
    assert!(
        nums[0].iter().all(|&m| m > 0.0),
        "every planted call has a margin"
    );

    let round = run::load(&annotated_path(&loaded.file, &out).to_string_lossy()).unwrap();
    let a = &round.manifest.annotate;
    assert!(a.argmax.is_some() && a.annotation.is_some());
    let share = a.cluster_celltype_share.as_deref().expect("a share table");
    let share = Mat::from_parquet(&resolve(&round.dir, share)).unwrap();
    assert_eq!(share.rows, ["K0".into(), "K1".into()] as [Box<str>; 2]);
    assert_eq!(share.mat[(0, 0)], 1.0);
    assert_eq!(share.mat[(1, 1)], 1.0);
    assert!(
        round.manifest.cluster.clusters.is_some(),
        "the round's cluster table"
    );
}

#[test]
fn a_score_round_lists_each_clusters_types_by_share_in_the_round_view() {
    let dir = tempfile::tempdir().unwrap();
    let (loaded, _) = toy_run(dir.path());
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    let panel = dir.path().join("panel.tsv").to_string_lossy().into_owned();
    annotate_by_projection(&score_args(&panel, &out), &loaded).unwrap();
    let view = crate::tui::round::RoundView::load(&annotated_path(&loaded.file, &out)).unwrap();
    for (c, want) in view.clusters.iter().zip(["A", "B"]) {
        let top = c.candidates.first().expect("a candidate");
        assert_eq!(
            (top.label.as_str(), top.share),
            (want, 1.0),
            "cluster {}",
            c.id
        );
    }
}
