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
