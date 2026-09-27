//! Round files and `relabel` end to end for [`super`].

use super::*;
use crate::manifest::run::RunManifest;
use std::io::Write as _;
use std::path::PathBuf;

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

/// A first round in `{root}/r0`: four cells in clusters 0, 0, 1, 2.
fn first_round(root: &Path) -> PathBuf {
    let r0 = root.join("r0");
    fs::create_dir_all(&r0).unwrap();
    let cells: Vec<Box<str>> = ["a", "b", "c", "d"].iter().map(|s| Box::from(*s)).collect();
    let at = |f: &str| r0.join(f).to_string_lossy().into_owned();
    write_clusters(
        &at("run.clusters.parquet"),
        &cells,
        &[Some(0), Some(0), Some(1), Some(2)],
    )
    .unwrap();
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

#[test]
fn relabel_writes_a_new_round_and_leaves_the_old_one() {
    let root = tempfile::tempdir().unwrap();
    let src = first_round(root.path());
    let before = fs::read_to_string(&src).unwrap();

    let decisions = root.path().join("d.jsonl");
    fs::write(
        &decisions,
        "{\"cluster\": 0, \"action\": \"label\", \"label\": \"CT4\", \"rationale\": \"why\", \"decided_by\": \"user\"}\n\
         {\"clusters\": [1, 2], \"action\": \"merge\", \"label\": \"CT2\", \"rationale\": \"same\", \"decided_by\": \"agent_proposed_user_accepted\"}\n",
    )
    .unwrap();
    let out = root.path().join("r1/next").to_string_lossy().into_owned();
    run_relabel(&RelabelArgs {
        from: src.to_string_lossy().into(),
        decisions: decisions.to_string_lossy().into(),
        out: out.clone().into(),
        watch: false,
    })
    .unwrap();

    assert_eq!(
        fs::read_to_string(&src).unwrap(),
        before,
        "the old round is untouched"
    );
    let next = run::load(&format!("{out}.senna.json")).unwrap();
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
    assert_eq!(
        Path::new(&from_next(&a.source)).canonicalize().unwrap(),
        src.canonicalize().unwrap()
    );
    let history: History =
        serde_json::from_str(&fs::read_to_string(from_next(&a.history)).unwrap()).unwrap();
    assert_eq!(history[&3][0].merged_from, Some(vec![1, 2]));
    assert_eq!(history[&0][0].rationale, "why");
    let summary: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(from_next(&a.cluster_summary)).unwrap()).unwrap();
    assert_eq!(summary["3"]["size"], 2);
    assert_eq!(summary["0"]["label"], "CT4");
    let log = fs::read_to_string(from_next(&a.log)).unwrap();
    assert_eq!(log.lines().count(), 2);
    assert!(log.contains("timestamp"));
}

#[test]
fn watch_turns_each_appended_batch_into_a_round() {
    let root = tempfile::tempdir().unwrap();
    let src = first_round(root.path());
    let decisions = root.path().join("d.jsonl");
    let out = root.path().join("w/run").to_string_lossy().into_owned();
    let status_path = PathBuf::from(format!("{out}{STATUS}"));
    fs::create_dir_all(root.path().join("w")).unwrap();
    let mut status = start_watch(
        &status_path,
        &src.to_string_lossy(),
        &decisions.to_string_lossy(),
    )
    .unwrap();
    let append = |text: &str| {
        let mut f = fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&decisions)
            .unwrap();
        f.write_all(text.as_bytes()).unwrap();
    };
    let step = |status: &mut WatchStatus| watch_step(status, &status_path, &out).unwrap();

    // No file yet, then a line without its newline: nothing to do.
    assert!(!step(&mut status));
    append(
        r#"{"cluster": 0, "action": "label", "label": "CT4", "rationale": "r", "decided_by": "user"}"#,
    );
    assert!(!step(&mut status));

    // The newline completes it: round 1.
    append("\n");
    assert!(step(&mut status));
    assert_eq!(status.rounds.len(), 1);
    assert!(status.error.is_none());

    // A refused batch is recorded and skipped; the round does not advance.
    append("{\"cluster\": 1, \"action\": \"keep\", \"decided_by\": \"user\"}\n");
    assert!(step(&mut status));
    assert_eq!(status.rounds.len(), 1);
    assert_eq!(status.error.as_ref().unwrap().lines, [2, 2]);

    // The next batch applies to round 1 and clears the error.
    append("{\"clusters\": [1, 2], \"action\": \"merge\", \"rationale\": \"r\", \"decided_by\": \"user\"}\n");
    assert!(step(&mut status));
    assert_eq!(status.rounds.len(), 2);
    assert!(status.error.is_none());

    let latest = run::load(&resolve(root.path().join("w").as_path(), &status.latest)).unwrap();
    let cells = read_cells(&latest).unwrap();
    assert_eq!(cells.clusters, vec![Some(0), Some(0), Some(3), Some(3)]);
    assert_eq!(
        cells.labels[0].as_deref(),
        Some("CT4"),
        "round 2 builds on round 1"
    );

    // A restarted watcher resumes rather than replaying.
    let resumed = start_watch(
        &status_path,
        &src.to_string_lossy(),
        &decisions.to_string_lossy(),
    )
    .unwrap();
    assert_eq!(resumed.processed_lines, 3);
    assert_eq!(resumed.latest, status.latest);
}
