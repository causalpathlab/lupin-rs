//! A pass written under a staging prefix and promoted onto its own.

use super::*;
use std::fs;

fn write(p: &Path, text: &str) {
    fs::write(p, text).unwrap();
}

fn read(p: &Path) -> String {
    fs::read_to_string(p).unwrap()
}

#[test]
fn a_finished_pass_replaces_the_round_and_its_stale_files_only() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let out = d.join("run.L1").to_string_lossy().into_owned();
    let staging = staging_prefix(&out);
    // The old round, a file of its that the new one no longer writes, a
    // later round and an unrelated file.
    write(
        &d.join("run.L1.senna.json"),
        r#"{"annotate": {"argmax": "run.L1.argmax.tsv", "cluster_celltype_pip": "run.L1.pip.parquet"}}"#,
    );
    write(&d.join("run.L1.argmax.tsv"), "old");
    write(&d.join("run.L1.pip.parquet"), "old");
    write(&d.join("run.L1.r1.senna.json"), "{}");
    write(&d.join("run.L1.notes.txt"), "mine");
    // The new pass, staged.
    write(
        &d.join(".run.L1.staging.senna.json"),
        &format!(
            r#"{{"annotate": {{"argmax": ".run.L1.staging.argmax.tsv", "markers": "panel.tsv",
                 "settings": {{"enrichment": {{"out": "{staging}"}}}}}}}}"#
        ),
    );
    write(&d.join(".run.L1.staging.argmax.tsv"), "new");

    promote(Path::new("run.senna.json"), &staging, &out).unwrap();

    let m: serde_json::Value = serde_json::from_str(&read(&d.join("run.L1.senna.json"))).unwrap();
    assert_eq!(m["annotate"]["argmax"], "run.L1.argmax.tsv");
    assert_eq!(m["annotate"]["markers"], "panel.tsv", "other paths stay");
    assert_eq!(m["annotate"]["settings"]["enrichment"]["out"], out.as_str());
    assert_eq!(read(&d.join("run.L1.argmax.tsv")), "new");
    assert!(
        !d.join("run.L1.pip.parquet").exists(),
        "the old round's, unused now"
    );
    assert!(
        d.join("run.L1.r1.senna.json").exists(),
        "a later round is not the pass's"
    );
    assert!(
        d.join("run.L1.notes.txt").exists(),
        "nor is a file the round never named"
    );
    let left: Vec<String> = fs::read_dir(d)
        .unwrap()
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.contains("staging"))
        .collect();
    assert!(left.is_empty(), "{left:?}");
}

#[test]
fn a_stopped_pass_leaves_the_round_and_its_leftovers_are_discarded() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let out = d.join("run.L1").to_string_lossy().into_owned();
    write(
        &d.join("run.L1.senna.json"),
        r#"{"annotate": {"argmax": "run.L1.argmax.tsv"}}"#,
    );
    write(&d.join("run.L1.argmax.tsv"), "old");
    write(&d.join(".run.L1.staging.argmax.tsv"), "half");
    fs::create_dir(d.join(".run.L1.staging.susie-staging")).unwrap();
    discard(&staging_prefix(&out)).unwrap();
    assert_eq!(read(&d.join("run.L1.argmax.tsv")), "old");
    assert!(!d.join(".run.L1.staging.argmax.tsv").exists());
    assert!(!d.join(".run.L1.staging.susie-staging").exists());
}

#[test]
fn the_staging_prefix_sits_beside_its_round_as_written() {
    assert_eq!(staging_prefix("run.L1"), ".run.L1.staging");
    assert_eq!(staging_prefix("out/run.L1"), "out/.run.L1.staging");
}

#[test]
fn a_prefix_written_with_a_trailing_slash_stages_beside_its_round() {
    assert_eq!(staging_prefix("results/run/"), "results/.run.staging");
}
