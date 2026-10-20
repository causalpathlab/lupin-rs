//! The default `lupin annotate`: every method the run supports, their
//! readings in one round, end to end on the planted toy run with its counts
//! ([`crate::annotate::marker_score::toy::write_run_with_counts`]).

use super::*;
use crate::annotate::geneset_units::{
    MEMBERSHIP_PARQUET, TYPE_LINKS_TSV, UNITS_CALL_PARQUET, UNITS_PARQUET, UNIT_ANCHORS_PARQUET,
};
use crate::annotate::marker_score::toy;
use crate::annotate::rounds::Method;
use crate::manifest::rounds::read_argmax;
use crate::manifest::run;
use legume_numeric::matrix::parquet::read_table_columns;

/// The default pass's settings on the toy's panel, writing under `out`.
fn auto_args(panel: &str, out: &str) -> AnnotateArgs {
    AnnotateArgs {
        markers: panel.into(),
        ..crate::annotate_cmd::default_enrichment_args(out)
    }
}

/// Run `args` on the toy run with counts in `dir`; the round it recorded.
fn pass(dir: &std::path::Path, args: impl FnOnce(&str, &str) -> AnnotateArgs) -> run::Loaded {
    let (manifest, panel) = toy::write_run_with_counts(dir);
    pass_on(&manifest, &panel, dir, args)
}

fn pass_on(
    manifest: &str,
    panel: &str,
    dir: &std::path::Path,
    args: impl FnOnce(&str, &str) -> AnnotateArgs,
) -> run::Loaded {
    let loaded = run::load(manifest).unwrap();
    let out = dir.join("run.L1").to_string_lossy().into_owned();
    annotate_by_enrichment(&args(panel, &out), &loaded, None).unwrap();
    run::load(&annotated_path(&loaded.file, &out).to_string_lossy()).unwrap()
}

#[test]
fn the_default_pass_records_every_methods_reading_in_one_round() {
    let dir = tempfile::tempdir().unwrap();
    let round = pass(dir.path(), auto_args);
    let a = &round.manifest.annotate;
    let ends = |p: &Option<String>, tail: &str| p.as_deref().is_some_and(|p| p.ends_with(tail));
    assert!(a.cluster_celltype_q.is_some(), "the enrichment's");
    assert!(
        ends(&a.marker_score, ".marker_score.parquet"),
        "the score's"
    );
    assert!(
        ends(&a.units_call, UNITS_CALL_PARQUET),
        "{:?}",
        a.units_call
    );
    assert!(ends(&a.celltype_units, UNITS_PARQUET));
    assert!(ends(&a.celltype_links, TYPE_LINKS_TSV));
    assert!(ends(&a.celltype_membership, MEMBERSHIP_PARQUET));
    assert!(ends(&a.unit_anchors, UNIT_ANCHORS_PARQUET));
    assert!(
        matches!(a.called_by, Some(Method::Susie | Method::Enrichment)),
        "the call stays SuSiE's or the enrichment's: {:?}",
        a.called_by
    );

    // The units' per-cell call, beside the round's.
    let path = resolve(&round.dir, a.units_call.as_deref().unwrap());
    let (strings, nums) = read_table_columns(&path, &["cell", "call"], &["membership"]).unwrap();
    assert_eq!(strings[0].len(), 20);
    for ((cell, call), m) in strings[0].iter().zip(&strings[1]).zip(&nums[0]) {
        let want = if cell.as_ref() < "c10" { "A" } else { "B" };
        assert_eq!(call.as_ref(), want, "{cell}");
        assert!(*m > 0.5, "{cell}: {m}");
    }

    let cells = crate::manifest::rounds::read_cells(&round).unwrap();
    let rows = crate::manifest::rounds::type_evidence(&round, &cells).unwrap();
    let a0 = rows
        .iter()
        .find(|r| r.label == "A" && r.score.is_some())
        .expect("A's evidence");
    assert_eq!((a0.score, a0.units), (Some(1.0), Some(1.0)), "{a0:?}");
    assert!(a0.enrichment.is_some(), "{a0:?}");

    let evidence = resolve(&round.dir, a.cluster_evidence.as_deref().expect("recorded"));
    let fields = legume_numeric::matrix::parquet::peek_parquet_field_names(&evidence).unwrap();
    assert!(fields.iter().any(|f| f.as_ref() == "units"), "{fields:?}");

    let v = crate::manifest::rounds::review_json(&round, &[], true).unwrap();
    let tops = v
        .values()
        .map(|c| c["tops"]["units"].clone())
        .collect::<Vec<_>>();
    assert_eq!(tops, ["A", "B"], "each cluster's units top");
}

#[test]
fn the_call_goes_to_the_method_asked_for() {
    for (call, by) in [(Method::Units, "units"), (Method::Score, "score")] {
        let dir = tempfile::tempdir().unwrap();
        let round = pass(dir.path(), |p, o| AnnotateArgs {
            call: Some(call),
            ..auto_args(p, o)
        });
        let a = &round.manifest.annotate;
        assert_eq!(a.called_by, Some(call), "{by}");
        let argmax = resolve(&round.dir, a.argmax.as_deref().unwrap());
        let labels = read_argmax(&argmax).unwrap();
        for i in 0..20 {
            let (label, p) = &labels[&format!("c{i:02}")];
            assert_eq!(label, if i < 10 { "A" } else { "B" }, "{by}: c{i:02}");
            match call {
                Method::Units => assert!(*p > 0.5, "the units' membership"),
                _ => assert!(p.is_nan(), "the score has no probability"),
            }
        }
        assert!(
            a.cluster_celltype_q.is_some() && a.marker_score.is_some() && a.units_call.is_some(),
            "{by}: every reading is still the round's"
        );
        assert!(
            a.annotation.is_none(),
            "{by}: no per-cell posterior but the caller's: {:?}",
            a.annotation
        );
    }
}

#[test]
fn susie_reads_beside_a_call_it_does_not_make() {
    let dir = tempfile::tempdir().unwrap();
    let round = pass(dir.path(), |p, o| AnnotateArgs {
        call: Some(Method::Enrichment),
        ..auto_args(p, o)
    });
    let a = &round.manifest.annotate;
    assert_eq!(a.called_by, Some(Method::Enrichment));
    assert!(
        a.cluster_celltype_pip.is_some(),
        "SuSiE's tables are still the round's"
    );
}

#[test]
fn units_that_fail_cost_their_reading_not_the_pass() {
    let dir = tempfile::tempdir().unwrap();
    let (manifest, panel) = toy::write_run_with_counts(dir.path());
    let mut m: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&manifest).unwrap()).unwrap();
    m["outputs"]["feature_coembedding"] = "gone.parquet".into();
    std::fs::write(&manifest, m.to_string()).unwrap();
    let round = pass_on(&manifest, &panel, dir.path(), auto_args);
    let a = &round.manifest.annotate;
    assert!(a.units_call.is_none() && a.celltype_units.is_none());
    assert!(a.marker_score.is_some(), "the score still beside the call");
    assert!(a.cluster_summary.is_some());
}

#[test]
fn a_call_by_a_method_that_failed_fails_the_pass() {
    let dir = tempfile::tempdir().unwrap();
    let (manifest, panel) = toy::write_run_with_counts(dir.path());
    let mut m: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&manifest).unwrap()).unwrap();
    m["outputs"]["feature_coembedding"] = "gone.parquet".into();
    std::fs::write(&manifest, m.to_string()).unwrap();
    let loaded = run::load(&manifest).unwrap();
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    let args = AnnotateArgs {
        call: Some(Method::Units),
        ..auto_args(&panel, &out)
    };
    let e = annotate_by_enrichment(&args, &loaded, None).unwrap_err();
    assert!(format!("{e:#}").contains("units"), "{e:#}");
}

#[test]
fn an_explicit_enrichment_runs_only_the_enrichment() {
    let dir = tempfile::tempdir().unwrap();
    let round = pass(dir.path(), |p, o| AnnotateArgs {
        score_beside: false,
        units_beside: false,
        ..auto_args(p, o)
    });
    let a = &round.manifest.annotate;
    assert!(a.cluster_celltype_q.is_some());
    assert!(a.marker_score.is_none() && a.units_call.is_none() && a.celltype_units.is_none());
}
