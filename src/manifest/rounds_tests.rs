//! Round files and `relabel` end to end for [`super`].

use super::*;
use crate::manifest::run::RunManifest;
use serde_json::{json, Value};
use std::io::Write as _;
use std::path::PathBuf;

/// A first round in `{root}/r0`: four cells in clusters 0, 0, 1, 2.
fn first_round(root: &Path) -> PathBuf {
    let r0 = root.join("r0");
    fs::create_dir_all(&r0).unwrap();
    let cells: Vec<Box<str>> = ["a", "b", "c", "d"].iter().map(|s| Box::from(*s)).collect();
    let at = |f: &str| r0.join(f).to_string_lossy().into_owned();
    let ids = [Some(0), Some(0), Some(1), Some(2)];
    write_clusters(&at("run.clusters.parquet"), &cells, &ids).unwrap();
    let labels: Vec<Option<String>> = ["CT1", "CT1", "CT2", "CT3"]
        .iter()
        .map(|s| Some((*s).to_string()))
        .collect();
    write_argmax(
        &at("run.argmax.tsv"),
        &cells,
        &labels,
        &[0.9, 0.8, 0.7, 0.6],
    )
    .unwrap();
    let mut m = RunManifest::new(crate::manifest::run::RunKind::Topic, "run");
    m.cluster.clusters = Some("run.clusters.parquet".into());
    m.annotate.argmax = Some("run.argmax.tsv".into());
    let src = r0.join("run.senna.json");
    m.save(&src).unwrap();
    src
}

/// One decisions-file line.
fn line(decision: Value) -> String {
    format!("{decision}\n")
}

/// Several decisions-file lines.
fn lines(decisions: &[Value]) -> String {
    decisions.iter().map(|d| line(d.clone())).collect()
}

/// A keep on cluster 0, optionally naming the round it was made on.
fn keep0(round: Option<&str>) -> String {
    let mut d = json!({"cluster": 0, "action": "keep", "rationale": "r", "decided_by": "user"});
    if let Some(r) = round {
        d["round"] = json!(r);
    }
    line(d)
}

fn args(from: &Path, decisions: &Path, out: Option<&Path>) -> RelabelArgs {
    RelabelArgs {
        from: from.to_string_lossy().into(),
        decisions: decisions.to_string_lossy().into(),
        out: out.map(|o| o.to_string_lossy().into()),
        ..RelabelArgs::default()
    }
}

/// `--next` from `from` with `text` as the decisions.
fn next_from(from: &Path, text: &str) -> Result<PathBuf> {
    let source = run::load(&from.to_string_lossy())?;
    let chain = Chain::of(&source.file);
    let lock = chain.lock()?;
    let ds = parse_decisions(numbered(text), "test")?;
    chain.write_next(&lock, &source.file, ds, Path::new("."))
}

/// Start a watch on `src`'s chain with decisions at `decisions`.
fn watch_on(src: &Path, decisions: &Path) -> (WatchStatus, PathBuf, Chain) {
    begin_watch(&RelabelArgs {
        watch: true,
        ..args(src, decisions, None)
    })
    .unwrap()
}

fn append_to(path: &Path, text: &str) {
    let mut f = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .unwrap();
    f.write_all(text.as_bytes()).unwrap();
}

#[test]
fn cluster_tables_round_trip_with_unassigned_cells() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("c.parquet").to_string_lossy().into_owned();
    let cells: Vec<Box<str>> = ["a", "b", "c"].iter().map(|s| Box::from(*s)).collect();
    write_clusters(&path, &cells, &[Some(4), None, Some(0)]).unwrap();
    let (back_cells, ids) = read_clusters(&path).unwrap();
    assert_eq!(back_cells, cells);
    assert_eq!(ids, vec![Some(4), None, Some(0)]);
}

#[test]
fn relabel_writes_a_new_round_and_leaves_the_old_one() {
    let root = tempfile::tempdir().unwrap();
    let src = first_round(root.path());
    let before = fs::read_to_string(&src).unwrap();

    let decisions = root.path().join("d.jsonl");
    let text = lines(&[
        json!({"cluster": 0, "action": "label", "label": "CT4", "rationale": "why", "decided_by": "user"}),
        json!({"clusters": [1, 2], "action": "merge", "label": "CT2", "rationale": "same", "decided_by": "agent_proposed_user_accepted"}),
    ]);
    fs::write(&decisions, text).unwrap();
    let out = root.path().join("r1/next");
    run_relabel(&args(&src, &decisions, Some(&out))).unwrap();

    assert_eq!(
        fs::read_to_string(&src).unwrap(),
        before,
        "the old round is untouched"
    );
    let next = run::load(&format!("{}.senna.json", out.display())).unwrap();
    let cells1 = read_cells(&next).unwrap();
    assert_eq!(cells1.clusters, vec![Some(0), Some(0), Some(3), Some(3)]);
    let l: Vec<&str> = cells1
        .labels
        .iter()
        .map(|l| l.as_deref().unwrap())
        .collect();
    assert_eq!(l, ["CT4", "CT4", "CT2", "CT2"]);

    let a = &next.manifest.annotate;
    let from_next = |rel: &Option<String>| resolve(&next.dir, rel.as_deref().unwrap());
    assert!(same_file(Path::new(&from_next(&a.source)), &src));
    let history: History =
        serde_json::from_str(&fs::read_to_string(from_next(&a.history)).unwrap()).unwrap();
    assert_eq!(history[&3][0].merged_from, Some(vec![1, 2]));
    assert_eq!(history[&0][0].rationale, "why");
    let summary: Value =
        serde_json::from_str(&fs::read_to_string(from_next(&a.cluster_summary)).unwrap()).unwrap();
    assert_eq!(summary["3"]["size"], 2);
    assert_eq!(summary["0"]["label"], "CT4");
    let log = fs::read_to_string(from_next(&a.log)).unwrap();
    assert_eq!(log.lines().count(), 2);
    assert!(log.contains("timestamp"));
}

#[test]
fn an_explicit_output_never_overwrites_a_round() {
    let root = tempfile::tempdir().unwrap();
    let src = first_round(root.path());
    let d = root.path().join("d.jsonl");
    fs::write(&d, keep0(None)).unwrap();
    let out = root.path().join("o/x");
    run_relabel(&args(&src, &d, Some(&out))).unwrap();
    let err = run_relabel(&args(&src, &d, Some(&out)))
        .unwrap_err()
        .to_string();
    assert!(err.contains("never overwritten"), "{err}");
}

#[test]
fn marker_edits_write_the_round_its_own_panel() {
    let root = tempfile::tempdir().unwrap();
    let src = first_round(root.path());
    let panel0 = root.path().join("r0/markers.tsv");
    let panel0_text = "gene\tcelltype\nGENE1\tCT1\nGENE2\tCT1\nGENE3\tCT2\n";
    fs::write(&panel0, panel0_text).unwrap();
    let mut m = RunManifest::load(&src).unwrap().0;
    m.annotate.markers = Some("markers.tsv".into());
    m.save(&src).unwrap();

    let relabel_with = |from: &Path, text: &str, out: &str| -> PathBuf {
        let d = root.path().join(format!("{}.jsonl", out.replace('/', "_")));
        fs::write(&d, text).unwrap();
        let out = root.path().join(out);
        run_relabel(&args(from, &d, Some(&out))).unwrap();
        PathBuf::from(format!("{}.senna.json", out.display()))
    };

    let r1 = relabel_with(
        &src,
        &lines(&[
            json!({"action": "markers_drop", "label": "CT1", "features": ["GENE2"], "rationale": "weak", "decided_by": "user"}),
            json!({"action": "markers_add", "label": "CT3", "features": ["GENE4"], "rationale": "new", "decided_by": "user"}),
        ]),
        "r1/run",
    );
    let next = run::load(&r1.to_string_lossy()).unwrap();
    let a = &next.manifest.annotate;
    let panel1 = resolve(&next.dir, a.markers.as_deref().unwrap());
    assert!(panel1.ends_with("run.markers.tsv"));
    assert_eq!(
        fs::read_to_string(&panel1).unwrap(),
        "gene\tcelltype\nGENE1\tCT1\nGENE3\tCT2\nGENE4\tCT3\n"
    );
    assert_eq!(
        fs::read_to_string(&panel0).unwrap(),
        panel0_text,
        "the earlier round's panel is untouched"
    );
    let hist: rounds::MarkerHistory = serde_json::from_str(
        &fs::read_to_string(resolve(&next.dir, a.marker_history.as_deref().unwrap())).unwrap(),
    )
    .unwrap();
    assert_eq!(hist["CT1"][0].rationale, "weak");
    assert_eq!(hist["CT3"][0].features, Some(vec!["GENE4".to_string()]));

    // A round with no marker edits keeps the previous round's panel.
    let r2 = relabel_with(&r1, &keep0(None), "r2/run");
    let after = run::load(&r2.to_string_lossy()).unwrap();
    let panel2 = resolve(
        &after.dir,
        after.manifest.annotate.markers.as_deref().unwrap(),
    );
    assert!(same_file(Path::new(&panel2), Path::new(&panel1)));
}

#[test]
fn rounds_are_named_by_chain() {
    let c = Chain::of(Path::new("d/L0.senna.json"));
    assert_eq!(c.round_prefix(1), "d/L0.r1");
    let c = Chain::of(Path::new("d/L0.r3.senna.json"));
    assert_eq!(c.round_prefix(4), "d/L0.r4");
    // A trailing `.r` without digits is part of the name.
    let c = Chain::of(Path::new("d/a.rb.senna.json"));
    assert_eq!(c.round_prefix(1), "d/a.rb.r1");
    // A pinto run's rounds are lupin manifests.
    let c = Chain::of(Path::new("d/run.pinto.json"));
    assert_eq!(
        annotated_path(&c.base, &c.round_prefix(1)),
        PathBuf::from("d/run.r1.lupin.json")
    );
}

#[test]
fn next_grows_the_chain_and_refuses_a_stale_round() {
    let root = tempfile::tempdir().unwrap();
    let src = first_round(root.path());
    let r1 = next_from(&src, &keep0(None)).unwrap();
    assert!(r1.ends_with("r0/run.r1.senna.json"), "{}", r1.display());
    let r2 = next_from(&r1, &keep0(None)).unwrap();
    assert!(r2.ends_with("r0/run.r2.senna.json"));

    // Deciding again on r1 (a view that has not reloaded) is refused.
    let err = next_from(&r1, &keep0(None)).unwrap_err().to_string();
    assert!(err.contains("not the latest round"), "{err}");
    assert!(!root.path().join("r0/run.r3.senna.json").exists());
}

#[test]
fn a_chain_lists_rounds_after_its_base_and_skips_gaps() {
    let root = tempfile::tempdir().unwrap();
    let src = first_round(root.path());
    let r1 = next_from(&src, &keep0(None)).unwrap();
    let r2 = next_from(&r1, &keep0(None)).unwrap();
    // Starting from r1, the chain after it is only r2.
    let after: Vec<usize> = Chain::of(&r1).rounds().iter().map(|(k, _)| *k).collect();
    assert_eq!(after, [2]);
    // A missing r1 does not hide r2: latest is still r2.
    fs::remove_file(&r1).unwrap();
    let (k, latest) = Chain::of(&src).latest();
    assert_eq!(k, 2);
    assert!(same_file(&latest, &r2));
}

#[test]
fn a_held_lock_makes_other_writers_wait_then_fail() {
    let root = tempfile::tempdir().unwrap();
    let chain = Chain::of(&root.path().join("c.senna.json"));
    let held = chain.try_lock().unwrap();
    assert!(chain.try_lock().is_err());
    drop(held);
    assert!(chain.try_lock().is_ok());
}

#[test]
fn a_watch_announces_itself_before_any_decision() {
    let root = tempfile::tempdir().unwrap();
    let src = first_round(root.path());
    let (_, status_path, _) = watch_on(&src, &root.path().join("d.jsonl"));
    let status: WatchStatus =
        serde_json::from_str(&fs::read_to_string(&status_path).unwrap()).unwrap();
    assert_eq!(status.processed_lines, 0);
    assert_eq!(status.rounds, vec![status.base.clone()]);
    assert_eq!(status.latest, status.base);
    assert!(status.error.is_none());
    let base = resolve(&parent_dir(&status_path), &status.base);
    assert!(same_file(Path::new(&base), &src));
}

#[test]
fn watch_turns_each_appended_batch_into_a_round() {
    let root = tempfile::tempdir().unwrap();
    let src = first_round(root.path());
    let decisions = root.path().join("d.jsonl");
    let (mut status, status_path, chain) = watch_on(&src, &decisions);
    let append = |text: &str| append_to(&decisions, text);
    let step = |status: &mut WatchStatus| watch_step(status, &status_path, &chain).unwrap();

    // No file yet, then a line without its newline: nothing to do.
    assert!(!step(&mut status));
    let label = json!({"cluster": 0, "action": "label", "label": "CT4", "rationale": "r", "decided_by": "user", "round": "r0/run.senna.json"});
    append(&label.to_string());
    assert!(!step(&mut status));

    // The newline completes it: round 1.
    append("\n");
    assert!(step(&mut status));
    assert_eq!(status.rounds, ["run.senna.json", "run.r1.senna.json"]);
    assert_eq!(status.latest, "run.r1.senna.json");
    assert!(status.error.is_none());

    // A refused batch (no rationale) is recorded and skipped.
    append(&line(
        json!({"cluster": 1, "action": "keep", "decided_by": "user", "round": "r0/run.r1.senna.json"}),
    ));
    assert!(step(&mut status));
    assert_eq!(status.rounds.len(), 2);
    assert_eq!(status.error.as_ref().unwrap().lines, [2, 2]);

    // A watched decision that does not name its round is refused.
    let merge =
        json!({"clusters": [1, 2], "action": "merge", "rationale": "r", "decided_by": "user"});
    append(&line(merge.clone()));
    assert!(step(&mut status));
    let msg = &status.error.as_ref().unwrap().message;
    assert!(msg.contains("must name the round"), "{msg}");

    // Naming it, the batch applies to round 1 and clears the error.
    let mut named = merge;
    named["round"] = json!("r0/run.r1.senna.json");
    append(&line(named));
    assert!(step(&mut status));
    assert_eq!(status.rounds.len(), 3);
    assert!(status.error.is_none());

    // A decision made on round 1 after round 2 landed is refused.
    append(&keep0(Some("r0/run.r1.senna.json")));
    assert!(step(&mut status));
    assert_eq!(status.rounds.len(), 3);
    let msg = &status.error.as_ref().unwrap().message;
    assert!(msg.contains("not the latest round"), "{msg}");

    let latest = run::load(&resolve(&parent_dir(&status_path), &status.latest)).unwrap();
    let cells = read_cells(&latest).unwrap();
    assert_eq!(cells.clusters, vec![Some(0), Some(0), Some(3), Some(3)]);
    assert_eq!(
        cells.labels[0].as_deref(),
        Some("CT4"),
        "round 2 builds on round 1"
    );

    // A restarted watcher resumes rather than replaying.
    let resumed = start_watch(&status_path, &src, &decisions.to_string_lossy()).unwrap();
    assert_eq!(resumed.processed_lines, 5);
    assert_eq!(resumed.latest, status.latest);
}

#[test]
fn a_watcher_waits_out_a_held_lock_without_losing_decisions() {
    let root = tempfile::tempdir().unwrap();
    let src = first_round(root.path());
    let decisions = root.path().join("d.jsonl");
    let (mut status, status_path, chain) = watch_on(&src, &decisions);
    append_to(&decisions, &keep0(Some("r0/run.senna.json")));

    // Another writer holds the chain: the batch waits, nothing is consumed.
    let held = chain.try_lock().unwrap();
    assert!(!watch_step(&mut status, &status_path, &chain).unwrap());
    assert_eq!(status.processed_lines, 0);
    assert!(status.error.is_none());

    // Once it is free, the same batch applies.
    drop(held);
    assert!(watch_step(&mut status, &status_path, &chain).unwrap());
    assert_eq!(status.rounds.len(), 2);
    assert!(status.error.is_none());
}

#[test]
fn a_watcher_follows_rounds_written_by_direct_calls() {
    let root = tempfile::tempdir().unwrap();
    let src = first_round(root.path());
    let decisions = root.path().join("d.jsonl");
    let (mut status, status_path, chain) = watch_on(&src, &decisions);

    // A viewer writes r1 directly; the watcher's status picks it up.
    let r1 = next_from(&src, &keep0(None)).unwrap();
    assert!(!watch_step(&mut status, &status_path, &chain).unwrap());
    let on_disk: WatchStatus =
        serde_json::from_str(&fs::read_to_string(&status_path).unwrap()).unwrap();
    assert_eq!(on_disk.latest, "run.r1.senna.json");

    // A decision appended to the file builds on r1, not on the base.
    fs::write(&decisions, keep0(Some("r0/run.r1.senna.json"))).unwrap();
    assert!(watch_step(&mut status, &status_path, &chain).unwrap());
    assert!(status.error.is_none(), "{:?}", status.error);
    assert_eq!(status.latest, "run.r2.senna.json");
    let r2 = run::load(&root.path().join("r0/run.r2.senna.json").to_string_lossy()).unwrap();
    let source = resolve(&r2.dir, r2.manifest.annotate.source.as_deref().unwrap());
    assert!(same_file(Path::new(&source), &r1));
}

#[test]
fn review_prints_the_names_decisions_files_use() {
    use crate::annotate::rounds::{Action, DecidedBy};
    for a in [
        Action::Label,
        Action::Merge,
        Action::Keep,
        Action::Split,
        Action::MarkersAdd,
        Action::MarkersDrop,
    ] {
        assert_eq!(json!(a), json!(a.as_str()));
    }
    for d in [
        DecidedBy::User,
        DecidedBy::AgentProposedUserAccepted,
        DecidedBy::UserOverride,
    ] {
        assert_eq!(json!(d), json!(d.as_str()));
    }
}
