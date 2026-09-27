//! Round files and `relabel` end to end for [`super`].

use super::*;
use crate::manifest::run::RunManifest;

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
    let r0 = root.path().join("r0");
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
