//! [`super::load_feature_embedding`]: a run's decoder gene table.

use super::*;
use crate::manifest::run;
use legume_numeric::matrix::traits::IoOps;

/// A run of `kind` in `dir` with a `[2 × 2]` table per named output slot,
/// each filled with `fill` so a test can tell which one was read.
fn run_with(dir: &Path, kind: &str, tables: &[(&str, &[&str], f32)]) -> run::Loaded {
    let at = |n: &str| dir.join(n).to_string_lossy().into_owned();
    let mut outputs = serde_json::Map::new();
    for (slot, rows, fill) in tables {
        let file = format!("run.{slot}.parquet");
        let rows: Vec<Box<str>> = rows.iter().map(|r| (*r).into()).collect();
        Mat::from_element(rows.len(), 2, *fill)
            .to_parquet_with_names(&at(&file), (Some(&rows[..]), Some("gene")), None)
            .unwrap();
        outputs.insert((*slot).into(), file.into());
    }
    let manifest =
        serde_json::json!({"version": 2, "kind": kind, "prefix": "run", "outputs": outputs});
    std::fs::write(at("run.senna.json"), manifest.to_string()).unwrap();
    run::load(&at("run.senna.json")).unwrap()
}

#[test]
fn the_decoder_table_is_rho_never_the_coembedding() {
    let dir = tempfile::tempdir().unwrap();
    let r = run_with(
        dir.path(),
        "bge",
        &[
            ("feature_embedding", &["A", "B"], 1.0),
            ("feature_coembedding", &["A", "B"], 2.0),
        ],
    );
    let f = load_feature_embedding(&r.manifest, &r.dir, "run").unwrap();
    assert!(f.mat.iter().all(|&x| x == 1.0));
    assert_eq!(f.rows, ["A".into(), "B".into()] as [Box<str>; 2]);
}

#[test]
fn a_gem_table_keeps_its_spliced_rows_by_gene() {
    let dir = tempfile::tempdir().unwrap();
    let rows = [
        "A/count/spliced",
        "A/count/unspliced",
        "B/count/spliced",
        "B/count/unspliced",
    ];
    let r = run_with(dir.path(), "gem", &[("feature_embedding", &rows, 1.0)]);
    let f = load_feature_embedding(&r.manifest, &r.dir, "run").unwrap();
    assert_eq!(f.rows, ["A".into(), "B".into()] as [Box<str>; 2]);
}

#[test]
fn a_run_without_a_decoder_table_points_to_enrichment() {
    let dir = tempfile::tempdir().unwrap();
    let r = run_with(
        dir.path(),
        "svd",
        &[("feature_embedding", &["A", "B"], 1.0)],
    );
    let err = load_feature_embedding(&r.manifest, &r.dir, "run")
        .err()
        .expect("an error")
        .to_string();
    assert!(
        err.contains("svd") && err.contains("--method enrichment"),
        "{err}"
    );
}

#[test]
fn a_decoder_run_missing_its_table_says_which_slot() {
    let dir = tempfile::tempdir().unwrap();
    let r = run_with(
        dir.path(),
        "bge",
        &[("feature_coembedding", &["A", "B"], 2.0)],
    );
    let err = load_feature_embedding(&r.manifest, &r.dir, "run")
        .err()
        .expect("an error")
        .to_string();
    assert!(err.contains("outputs.feature_embedding"), "{err}");
}
