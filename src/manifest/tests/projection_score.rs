//! `annotate --method projection --projection score`, end to end on the
//! planted toy run ([`crate::annotate::marker_score::toy`]).

use super::*;
use crate::annotate::marker_score::toy;
use crate::annotate::outputs::MARKER_SCORE_PARQUET;
use crate::manifest::rounds::read_argmax;
use crate::manifest::run;

#[test]
fn the_score_pass_calls_each_cell_and_records_its_round() {
    let dir = tempfile::tempdir().unwrap();
    let (manifest, panel) = toy::write_run(dir.path());
    let loaded = run::load(&manifest).unwrap();
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    annotate_by_marker_score(&toy::args(&panel, &out), &loaded).unwrap();

    let argmax = read_argmax(&format!("{out}.argmax.tsv")).unwrap();
    for i in 0..20 {
        let (label, p) = &argmax[&format!("c{i:02}")];
        assert_eq!(label, if i < 10 { "A" } else { "B" }, "c{i:02}");
        assert!(p.is_nan(), "the score has no probability");
    }
    let (_, nums) =
        read_table_columns(&format!("{out}{MARKER_SCORE_PARQUET}"), &[], &["margin"]).unwrap();
    assert!(
        nums[0].iter().all(|&m| m > 0.0),
        "every planted call has a margin"
    );

    let round = run::load(&annotated_path(&loaded.file, &out).to_string_lossy()).unwrap();
    let a = &round.manifest.annotate;
    assert!(a.argmax.is_some());
    assert!(
        a.marker_score
            .as_deref()
            .is_some_and(|p| p.ends_with(MARKER_SCORE_PARQUET)),
        "the score's own per-cell table: {:?}",
        a.marker_score
    );
    assert!(
        a.annotation.is_none(),
        "`annotation` is a cell × type posterior, which the score has none of"
    );
    assert!(
        round.manifest.cluster.clusters.is_some(),
        "the round's cluster table"
    );
    let recorded = &a.settings.as_ref().unwrap()["projection"];
    assert_eq!(recorded["projection"], "score");
    assert!(
        recorded.get("num_perm").is_none(),
        "no ORA settings on a score round"
    );
}

#[test]
fn a_score_round_clears_an_earlier_rounds_posterior() {
    let dir = tempfile::tempdir().unwrap();
    let (manifest, panel) = toy::write_run(dir.path());
    let mut m: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&manifest).unwrap()).unwrap();
    m["annotate"] = serde_json::json!({"annotation": "run.L0.annotation.parquet"});
    std::fs::write(&manifest, m.to_string()).unwrap();
    let loaded = run::load(&manifest).unwrap();
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    annotate_by_marker_score(&toy::args(&panel, &out), &loaded).unwrap();
    let round = run::load(&annotated_path(&loaded.file, &out).to_string_lossy()).unwrap();
    assert_eq!(
        round.manifest.annotate.annotation, None,
        "no earlier posterior beside the score's calls"
    );
}

/// A score round on the toy run, with an enrichment Q table added beside it
/// (cluster 0: A 0.7, B 0.3; cluster 1: B 1.0): two methods' readings.
pub(crate) fn round_with_both(dir: &std::path::Path) -> run::Loaded {
    use legume_numeric::matrix::dense_mat_io::Mat;
    use legume_numeric::matrix::traits::IoOps;
    let (manifest, panel) = toy::write_run(dir);
    let loaded = run::load(&manifest).unwrap();
    let out = dir.join("run.L1").to_string_lossy().into_owned();
    annotate_by_marker_score(&toy::args(&panel, &out), &loaded).unwrap();
    let round_path = annotated_path(&loaded.file, &out);
    let q = Mat::from_row_slice(2, 2, &[0.7, 0.3, 0.0, 1.0]);
    q.to_parquet_with_names(
        &format!("{out}.cluster_celltype_q.parquet"),
        (Some(&["K0".into(), "K1".into()][..]), Some("cluster")),
        Some(&["A".into(), "B".into()]),
    )
    .unwrap();
    let mut round = run::load(&round_path.to_string_lossy()).unwrap();
    round.manifest.annotate.cluster_celltype_q = Some("run.L1.cluster_celltype_q.parquet".into());
    round.manifest.save(&round.file).unwrap();
    run::load(&round_path.to_string_lossy()).unwrap()
}

#[test]
fn every_types_evidence_from_every_method_is_one_table() {
    let dir = tempfile::tempdir().unwrap();
    let mut round = round_with_both(dir.path());
    let cells = crate::manifest::rounds::read_cells(&round).unwrap();
    let rows = crate::manifest::rounds::type_evidence(&round, &cells).unwrap();
    let get = |k: u32, l: &str| {
        rows.iter()
            .find(|r| r.cluster == k && r.label == l)
            .unwrap()
    };
    assert_eq!(
        (get(0, "A").enrichment, get(0, "A").score),
        (Some(0.7), Some(1.0))
    );
    assert_eq!(
        get(0, "B").enrichment,
        Some(0.3),
        "the enrichment's runner-up"
    );
    assert_eq!(
        get(0, "B").score,
        Some(0.0),
        "the score called none of cluster 0 B"
    );
    assert_eq!(
        (get(1, "B").enrichment, get(1, "B").score),
        (Some(1.0), Some(1.0))
    );

    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    crate::manifest::rounds::write_summary(&mut round, &out).unwrap();
    let rel = round
        .manifest
        .annotate
        .cluster_evidence
        .as_deref()
        .expect("recorded");
    let path = resolve(&round.dir, rel);
    let fields = legume_numeric::matrix::parquet::peek_parquet_field_names(&path).unwrap();
    for f in [
        "cluster",
        "celltype",
        "share",
        "nes",
        "p",
        "q",
        "pip",
        "fold",
        "explained",
        "score",
    ] {
        assert!(fields.iter().any(|x| x.as_ref() == f), "{f} in {fields:?}");
    }
}

#[test]
fn beside_an_enrichment_the_score_takes_its_cells_and_clusters() {
    let dir = tempfile::tempdir().unwrap();
    let (manifest, panel) = toy::write_run(dir.path());
    let loaded = run::load(&manifest).unwrap();
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    // The enrichment's cells in its own order, one cell of the run left out,
    // in clusters of its own numbering.
    let names: Vec<Box<str>> = (0..20)
        .rev()
        .filter(|&i| i != 3)
        .map(|i| format!("c{i:02}").into())
        .collect();
    let groups: Vec<usize> = names
        .iter()
        .map(|n| usize::from(n.as_ref() < "c10") * 7)
        .collect();
    let table = score_beside(&loaded, &panel, &out, &names, &groups)
        .unwrap()
        .expect("the run decodes cells");
    let (strings, _) = read_table_columns(&table, &["cell", "call"], &[]).unwrap();
    assert_eq!(strings[0].len(), 19, "the enrichment's cells only");
    for (cell, call) in strings[0].iter().zip(&strings[1]) {
        let want = if cell.as_ref() < "c10" { "A" } else { "B" };
        assert_eq!(call.as_ref(), want, "{cell}");
    }
}

#[test]
fn a_run_that_does_not_decode_has_no_score_beside_its_enrichment() {
    let dir = tempfile::tempdir().unwrap();
    let (manifest, panel) = toy::write_run(dir.path());
    let mut loaded = run::load(&manifest).unwrap();
    loaded.manifest.outputs.feature_embedding = None;
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    let names: Vec<Box<str>> = vec!["c00".into()];
    assert!(score_beside(&loaded, &panel, &out, &names, &[0])
        .unwrap()
        .is_none());
}

#[test]
fn each_method_ranks_the_types_it_puts_forward() {
    use crate::annotate::rounds::Method;
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_both(dir.path());
    let cells = crate::manifest::rounds::read_cells(&round).unwrap();
    let rows = crate::manifest::rounds::type_evidence(&round, &cells).unwrap();
    let k0: Vec<_> = rows.iter().filter(|r| r.cluster == 0).cloned().collect();
    let labels =
        |m: Method| -> Vec<String> { m.ranked(&k0).iter().map(|r| r.label.clone()).collect() };
    assert_eq!(labels(Method::Enrichment), ["A", "B"]);
    assert_eq!(
        labels(Method::Score),
        ["A"],
        "a type no cell was called is not put forward"
    );
    assert!(labels(Method::Susie).is_empty(), "SuSiE did not run");
}

#[test]
fn review_json_gives_each_methods_top_and_every_types_evidence() {
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_both(dir.path());
    let v =
        serde_json::to_value(crate::manifest::rounds::review_json(&round, &[]).unwrap()).unwrap();
    let k0 = &v["0"];
    assert_eq!(k0["tops"]["enrichment"], "A");
    assert_eq!(k0["tops"]["score"], "A");
    assert!(k0["tops"].get("susie").is_none(), "{}", k0["tops"]);
    let a = k0["types"]
        .as_array()
        .unwrap()
        .iter()
        .find(|t| t["celltype"] == "A")
        .expect("A's evidence");
    assert!((a["share"].as_f64().unwrap() - 0.7).abs() < 1e-6);
    assert_eq!(a["score"], 1.0);
    assert!(a.get("pip").is_none(), "no SuSiE reading");
    assert!(k0.get("digest").is_some() && k0.get("history").is_some());
}

#[test]
fn the_scores_calls_take_the_panels_spelling() {
    use legume_numeric::matrix::parquet::{write_named_table, Column};
    let dir = tempfile::tempdir().unwrap();
    let (manifest, _) = toy::write_run(dir.path());
    let mut loaded = run::load(&manifest).unwrap();
    let call: Vec<Box<str>> = vec!["CD4  T, cell".into()];
    write_named_table(
        &dir.path()
            .join("run.marker_score.parquet")
            .to_string_lossy(),
        "cell",
        &["c00".into()],
        &[
            ("cluster".into(), Column::I32(&[0])),
            ("call".into(), Column::Str(&call)),
            ("margin".into(), Column::F32(&[0.5])),
        ],
    )
    .unwrap();
    loaded.manifest.annotate.marker_score = Some("run.marker_score.parquet".into());
    let cells = crate::manifest::rounds::read_cells(&loaded).unwrap();
    let score = crate::manifest::rounds::read_score_cells(&loaded, &cells)
        .unwrap()
        .unwrap();
    assert_eq!(score.types, ["CD4_T_cell".into()] as [Box<str>; 1]);
}

#[test]
fn a_rounds_method_is_the_one_it_records_making_its_call() {
    use crate::annotate::rounds::Method;
    use crate::manifest::rounds::called_by;
    let dir = tempfile::tempdir().unwrap();
    let mut round = round_with_both(dir.path());
    assert_eq!(
        called_by(&round.manifest.annotate),
        Some(Method::Score),
        "a score round with a Q table beside it is still the score's"
    );
    round.manifest.annotate.drop_marker_score();
    assert_eq!(round.manifest.annotate.called_by, None);
    assert_eq!(
        called_by(&round.manifest.annotate),
        Some(Method::Enrichment),
        "unrecorded, the tables tell"
    );
}

#[test]
fn review_reads_a_round_from_its_summary_when_its_cells_are_gone() {
    let dir = tempfile::tempdir().unwrap();
    let mut round = round_with_both(dir.path());
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    crate::manifest::rounds::write_summary(&mut round, &out).unwrap();
    round.manifest.cluster.clusters = None;
    let v = crate::manifest::rounds::review_json(&round, &[]).unwrap();
    assert!(v.contains_key(&0), "{v:?}");
    assert_eq!(v[&0]["tops"]["enrichment"], "A", "the tables need no cells");
}

#[test]
fn an_unreadable_table_costs_the_evidence_not_the_summary() {
    let dir = tempfile::tempdir().unwrap();
    let mut round = round_with_both(dir.path());
    std::fs::write(
        dir.path().join("run.L1.cluster_celltype_q.parquet"),
        "not parquet",
    )
    .unwrap();
    round.manifest.annotate.cluster_evidence = Some("an earlier round's".into());
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    crate::manifest::rounds::write_summary(&mut round, &out).unwrap();
    assert!(round.manifest.annotate.cluster_summary.is_some());
    assert_eq!(
        round.manifest.annotate.cluster_evidence, None,
        "no evidence, and none of another round's"
    );
}

#[test]
fn the_score_pass_writes_its_labels_as_parquet_too() {
    let dir = tempfile::tempdir().unwrap();
    let (manifest, panel) = toy::write_run(dir.path());
    let loaded = run::load(&manifest).unwrap();
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    annotate_by_marker_score(&toy::args(&panel, &out), &loaded).unwrap();
    let labels = read_argmax(&format!("{out}.argmax.parquet")).unwrap();
    assert_eq!(labels["c00"].0, "A");
}

#[test]
fn review_settles_tied_calls_names_clusters_k_and_places_them_on_the_tree() {
    use legume_numeric::matrix::dense_mat_io::Mat;
    use legume_numeric::matrix::traits::IoOps;
    let dir = tempfile::tempdir().unwrap();
    let mut round = round_with_both(dir.path());
    // Both types tied on q in cluster 1, whose share is all B's.
    let q = Mat::from_row_slice(2, 2, &[0.01, 0.5, 0.01, 0.01]);
    q.to_parquet_with_names(
        &dir.path()
            .join("run.L1.cluster_celltype_q_values.parquet")
            .to_string_lossy(),
        (Some(&["K0".into(), "K1".into()][..]), Some("cluster")),
        Some(&["A".into(), "B".into()]),
    )
    .unwrap();
    std::fs::write(
        dir.path().join("run.L1.celltype_tree.json"),
        r#"{"source":"cell_ontology","label_cl":{"A":"CL:9000001"},
            "groups":[{"name":"group_ab","cl_id":"CL:9000009","members":["A","B"]}]}"#,
    )
    .unwrap();
    let a = &mut round.manifest.annotate;
    a.cluster_celltype_q_values = Some("run.L1.cluster_celltype_q_values.parquet".into());
    a.celltype_tree = Some("run.L1.celltype_tree.json".into());
    // The summary as a pass with these tables writes it.
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    crate::manifest::rounds::write_summary(&mut round, &out).unwrap();
    round.manifest.save(&round.file).unwrap();
    let round = run::load(&round.file.to_string_lossy()).unwrap();

    let v = crate::manifest::rounds::review_json(&round, &[]).unwrap();
    let e = &v[&1]["digest"]["evidence"];
    assert_eq!(e["top"], "B", "{e}");
    assert_eq!(e["agrees"], true, "{e}");
    assert_eq!(v[&0]["cl"]["id"], "CL:9000001");
    assert_eq!(v[&0]["cl"]["via"], "label");
    assert_eq!(v[&1]["cl"]["id"], "CL:9000009");
    assert_eq!(v[&1]["cl"]["name"], "group_ab");

    let text = crate::manifest::rounds::review_text(&round, &[]).unwrap();
    assert!(text.starts_with("K0  "), "{text}");
    assert!(text.contains("\nK1  "), "{text}");
    assert!(text.contains("CL    CL:9000001 A\n"), "{text}");
    assert!(
        text.contains("CL    CL:9000009 group_ab (the group of B)\n"),
        "{text}"
    );
    let calls: Vec<&str> = text
        .lines()
        .skip_while(|l| !l.starts_with("K1"))
        .filter_map(|l| l.trim().strip_prefix("call"))
        .map(|l| l.split_whitespace().next().unwrap())
        .collect();
    assert_eq!(calls, ["B", "A"], "{text}");
}
