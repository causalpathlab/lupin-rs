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
fn superseding_moves_only_the_later_rounds_files() {
    let root = tempfile::tempdir().unwrap();
    let dir = root.path();
    for f in [
        "run.senna.json",
        "run.r1.senna.json",
        "run.r1.argmax.tsv",
        "run.r2.senna.json",
        "run.rest.keep",
        "run.rx.keep",
        "other.r1.senna.json",
    ] {
        fs::write(dir.join(f), "{}").unwrap();
    }
    let moved = supersede_later(&dir.join("run.senna.json")).unwrap();
    assert_eq!(moved, 2);
    let aside = dir.join("run.superseded");
    for f in [
        "run.r1.senna.json",
        "run.r1.argmax.tsv",
        "run.r2.senna.json",
    ] {
        assert!(aside.join(f).is_file(), "{f} set aside");
        assert!(!dir.join(f).exists(), "{f} gone from the chain");
    }
    for f in [
        "run.senna.json",
        "run.rest.keep",
        "run.rx.keep",
        "other.r1.senna.json",
    ] {
        assert!(dir.join(f).is_file(), "{f} left alone");
    }
    assert_eq!(supersede_later(&dir.join("run.senna.json")).unwrap(), 0);
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

/// `first_round` plus an expression profile (genes × K0..K2) and a panel:
/// GENE1 marks CT1 and peaks in K0, GENE2 marks CT2 and peaks in K1, GENE3
/// marks CT3 and peaks in K2; GENE4, on no panel, also peaks in K1.
fn round_with_profile(root: &Path) -> PathBuf {
    let src = first_round(root);
    let r0 = root.join("r0");
    let genes: Vec<Box<str>> = ["GENE1", "GENE2", "GENE3", "GENE4"]
        .iter()
        .map(|s| Box::from(*s))
        .collect();
    let cols: Vec<Box<str>> = ["K0", "K1", "K2"].iter().map(|s| Box::from(*s)).collect();
    let mut m = Mat::zeros(4, 3);
    for (g, k) in [(0, 0), (1, 1), (2, 2), (3, 1)] {
        m[(g, k)] = 50.0;
    }
    for g in 0..4 {
        for k in 0..3 {
            m[(g, k)] += 1.0;
        }
    }
    let profile = r0.join("run.cluster_expression.parquet");
    m.to_parquet_with_names(
        &profile.to_string_lossy(),
        (Some(&genes), Some("gene")),
        Some(&cols),
    )
    .unwrap();
    fs::write(
        r0.join("markers.tsv"),
        "gene\tcelltype\nGENE1\tCT1\nGENE2\tCT2\nGENE3\tCT3\n",
    )
    .unwrap();
    let mut man = RunManifest::load(&src).unwrap().0;
    man.annotate.cluster_expression = Some("run.cluster_expression.parquet".into());
    man.annotate.expression_clusters = Some("run.clusters.parquet".into());
    man.annotate.markers = Some("markers.tsv".into());
    man.save(&src).unwrap();
    src
}

fn preview_of(src: &Path, text: &str) -> Value {
    let source = run::load(&src.to_string_lossy()).unwrap();
    let ds = parse_decisions(numbered(text), "test").unwrap();
    preview(&source, ds, Path::new(".")).unwrap()
}

#[test]
fn preview_rescores_marker_edits_and_writes_nothing() {
    let root = tempfile::tempdir().unwrap();
    let src = round_with_profile(root.path());
    let before: Vec<_> = fs::read_dir(root.path().join("r0")).unwrap().collect();

    let p = preview_of(
        &src,
        &lines(&[
            json!({"cluster": 0, "action": "label", "label": "CT4", "rationale": "r", "decided_by": "user"}),
            json!({"action": "markers_drop", "label": "CT2", "features": ["GENE2"], "rationale": "r", "decided_by": "user"}),
            json!({"action": "markers_add", "label": "CT3", "features": ["GENE4"], "rationale": "r", "decided_by": "user"}),
        ]),
    );
    assert_eq!(p["rescored"], true);
    assert_eq!(p["cells_changed"], 2);
    // The named cluster: its label changes; its top call does not.
    assert_eq!(p["clusters"]["0"]["label_before"], "CT1");
    assert_eq!(p["clusters"]["0"]["label_after"], "CT4");
    assert_eq!(p["clusters"]["0"]["calls"][0]["label"], "CT1");
    // Not named, but the edited panel moves its top call.
    assert_eq!(p["clusters"]["1"]["top_before"], "CT2");
    assert_eq!(p["clusters"]["1"]["calls"][0]["label"], "CT3");
    assert_eq!(
        p["clusters"]["1"]["label_before"],
        p["clusters"]["1"]["label_after"]
    );
    // Untouched and unmoved: left out.
    assert!(p["clusters"].get("2").is_none(), "{p}");
    assert_eq!(p["markers"]["CT2"]["dropped"], json!(["GENE2"]));
    assert_eq!(p["markers"]["CT3"]["added"], json!(["GENE4"]));

    let after: Vec<_> = fs::read_dir(root.path().join("r0")).unwrap().collect();
    assert_eq!(before.len(), after.len(), "preview writes nothing");
}

#[test]
fn preview_shows_a_merge_under_its_new_id() {
    let root = tempfile::tempdir().unwrap();
    let src = round_with_profile(root.path());
    let p = preview_of(
        &src,
        &line(
            json!({"clusters": [1, 2], "action": "merge", "label": "CT2", "rationale": "r", "decided_by": "user"}),
        ),
    );
    let merged = &p["clusters"]["3"];
    assert_eq!(merged["label_after"], "CT2");
    assert!(
        merged["calls"].as_array().is_some_and(|c| !c.is_empty()),
        "{p}"
    );
    assert_eq!(p["cells_changed"], 1);
}

#[test]
fn preview_refuses_a_stale_round_and_reports_unscored_rounds() {
    let root = tempfile::tempdir().unwrap();
    let src = first_round(root.path());
    // No profile: labels only.
    let p = preview_of(&src, &keep0(None));
    assert_eq!(p["rescored"], false);
    assert!(
        p["clusters"].get("0").is_some(),
        "a named cluster is always shown"
    );

    next_from(&src, &keep0(None)).unwrap();
    let d = root.path().join("d.jsonl");
    fs::write(&d, keep0(None)).unwrap();
    let err = run_relabel(&RelabelArgs {
        preview: true,
        ..args(&src, &d, None)
    })
    .unwrap_err()
    .to_string();
    assert!(err.contains("not the latest round"), "{err}");
}

/// A round an enrichment pass wrote, with its statistics cached: 60 genes,
/// clusters 0..3 of 10 cells each over 4 batches; CT1, CT2 and CT3 own five
/// markers each, raised in clusters 0, 1 and 2.
fn enriched_round(root: &Path) -> PathBuf {
    use clap::Parser;
    #[derive(Parser)]
    struct Cli {
        #[command(flatten)]
        a: crate::annotate_cmd::AnnotateCliArgs,
    }
    let dir = root.join("e");
    fs::create_dir_all(&dir).unwrap();
    let at = |f: &str| dir.join(f).to_string_lossy().into_owned();
    let (g, k, per) = (60usize, 3usize, 10usize);
    let genes: Vec<Box<str>> = (0..g).map(|i| format!("GENE{i}").into()).collect();
    let cells: Vec<Box<str>> = (0..k * per).map(|i| format!("c{i}").into()).collect();
    let ids: Vec<Option<ClusterId>> = (0..k * per).map(|i| Some((i / per) as ClusterId)).collect();
    write_clusters(&at("run.clusters.parquet"), &cells, &ids).unwrap();
    let batches: Vec<Option<ClusterId>> =
        (0..k * per).map(|i| Some((i % 4) as ClusterId)).collect();
    write_clusters(&at("run.cell_batch.parquet"), &cells, &batches).unwrap();
    let labels: Vec<Option<String>> = ids
        .iter()
        .map(|id| Some(format!("CT{}", id.unwrap() + 1)))
        .collect();
    write_argmax(&at("run.argmax.tsv"), &cells, &labels, &vec![0.9; k * per]).unwrap();

    let mut sums = Mat::zeros(g, k);
    for i in 0..g {
        for c in 0..k {
            sums[(i, c)] = 10.0 + ((i * 7 + c * 3) % 5) as f32;
        }
    }
    for c in 0..k {
        for m in 0..5 {
            sums[(c * 5 + m, c)] *= 6.0;
        }
    }
    let kcols: Vec<Box<str>> = (0..k).map(|c| format!("K{c}").into()).collect();
    sums.to_parquet_with_names(
        &at("run.cluster_gene_sum.parquet"),
        (Some(&genes), Some("gene")),
        Some(&kcols),
    )
    .unwrap();
    sums.to_parquet_with_names(
        &at("run.cluster_expression.parquet"),
        (Some(&genes), Some("gene")),
        Some(&kcols),
    )
    .unwrap();
    let mut pb = Mat::zeros(g, 4);
    for i in 0..g {
        for b in 0..4 {
            pb[(i, b)] =
                (0..k).map(|c| sums[(i, c)]).sum::<f32>() / 300.0 * (1.0 + 0.05 * b as f32);
        }
    }
    let bcols: Vec<Box<str>> = (0..4).map(|b| format!("B{b}").into()).collect();
    pb.to_parquet_with_names(
        &at("run.batch_profile.parquet"),
        (Some(&genes), Some("gene")),
        Some(&bcols),
    )
    .unwrap();
    let mut w = Mat::zeros(g, 1);
    for i in 0..g {
        w[(i, 0)] = 1.0;
    }
    w.to_parquet_with_names(
        &at("run.gene_weight.parquet"),
        (Some(&genes), Some("gene")),
        Some(&["weight".into()]),
    )
    .unwrap();
    let mut panel = String::from("gene\tcelltype\n");
    for c in 0..k {
        for m in 0..5 {
            panel.push_str(&format!("GENE{}\tCT{}\n", c * 5 + m, c + 1));
        }
    }
    fs::write(dir.join("markers.tsv"), panel).unwrap();

    let cli = Cli::parse_from(["x", "-o", "o", "-m", "markers.tsv"]);
    let settings =
        serde_json::to_value(crate::annotate_cmd::build_enrichment_args(&cli.a)).unwrap();
    let mut m = RunManifest::new(crate::manifest::run::RunKind::Topic, "run");
    m.cluster.clusters = Some("run.clusters.parquet".into());
    let a = &mut m.annotate;
    a.argmax = Some("run.argmax.tsv".into());
    a.markers = Some("markers.tsv".into());
    a.expression_clusters = Some("run.clusters.parquet".into());
    a.cluster_expression = Some("run.cluster_expression.parquet".into());
    a.settings = Some(json!({ "enrichment": settings }));
    a.stats_cache = Some(crate::manifest::run::StatsCache {
        gene_sum: "run.cluster_gene_sum.parquet".into(),
        batch_profile: "run.batch_profile.parquet".into(),
        gene_weight: "run.gene_weight.parquet".into(),
        cell_batch: "run.cell_batch.parquet".into(),
    });
    let src = dir.join("run.senna.json");
    m.save(&src).unwrap();
    src
}

#[test]
fn a_relabel_round_is_rescored_on_its_merged_clusters() {
    let root = tempfile::tempdir().unwrap();
    let src = enriched_round(root.path());
    let d = root.path().join("d.jsonl");
    fs::write(
        &d,
        line(
            json!({"clusters": [1, 2], "action": "merge", "rationale": "r", "decided_by": "user"}),
        ),
    )
    .unwrap();
    run_relabel(&args(&src, &d, Some(&root.path().join("r1/run")))).unwrap();

    let r1 = run::load(&root.path().join("r1/run.senna.json").to_string_lossy()).unwrap();
    let a = &r1.manifest.annotate;
    let q = read_table(&resolve(
        &r1.dir,
        a.cluster_celltype_q_values.as_deref().unwrap(),
    ))
    .unwrap();
    assert_eq!(
        q.rows,
        vec![0, 3],
        "the merged cluster is scored under its new id"
    );
    let stats = a.stats.as_ref().unwrap();
    assert_eq!(stats["kind"], "post_selection");
    assert_eq!(stats["rounds_of_curation"], 1);

    let summary: Value = serde_json::from_str(
        &fs::read_to_string(resolve(&r1.dir, a.cluster_summary.as_deref().unwrap())).unwrap(),
    )
    .unwrap();
    assert_eq!(summary["0"]["evidence"]["top"], "CT1");
    assert_eq!(summary["0"]["evidence"]["agrees"], true);
    assert!(
        summary["0"]["evidence"]["q"].as_f64().unwrap() < 0.1,
        "{summary}"
    );
}

#[test]
fn a_marker_edit_rescores_the_round_and_its_preview() {
    let root = tempfile::tempdir().unwrap();
    let src = enriched_round(root.path());
    let edit = line(
        json!({"action": "markers_drop", "label": "CT3", "features": ["GENE14"], "rationale": "r", "decided_by": "user"}),
    );

    let p = preview_of(&src, &edit);
    assert_eq!(p["stats"], "recalibrated");

    let d = root.path().join("d.jsonl");
    fs::write(&d, edit).unwrap();
    run_relabel(&args(&src, &d, Some(&root.path().join("r1/run")))).unwrap();
    let r1 = run::load(&root.path().join("r1/run.senna.json").to_string_lossy()).unwrap();
    let a = &r1.manifest.annotate;
    let q = read_table(&resolve(
        &r1.dir,
        a.cluster_celltype_q_values.as_deref().unwrap(),
    ))
    .unwrap();
    assert_eq!(q.rows, vec![0, 1, 2]);
}

#[test]
fn a_preview_rescores_the_touched_types_as_a_save_does() {
    let root = tempfile::tempdir().unwrap();
    let src = enriched_round(root.path());
    let edit = line(
        json!({"action": "markers_drop", "label": "CT3", "features": ["GENE14"], "rationale": "r", "decided_by": "user"}),
    );
    let p = preview_of(&src, &edit);
    assert_eq!(p["stats"], "recalibrated");

    let d = root.path().join("d.jsonl");
    fs::write(&d, edit).unwrap();
    run_relabel(&args(&src, &d, Some(&root.path().join("r1/run")))).unwrap();
    let r1 = run::load(&root.path().join("r1/run.senna.json").to_string_lossy()).unwrap();
    let saved = read_table(&resolve(
        &r1.dir,
        r1.manifest.annotate.cluster_celltype_p.as_deref().unwrap(),
    ))
    .unwrap();
    let j = saved.cols.iter().position(|c| c == "CT3").unwrap();
    for id in &saved.rows {
        let calls = p["scores"][id.to_string()].as_array().unwrap();
        let ct3 = calls.iter().find(|c| c["label"] == "CT3").unwrap();
        let (preview_p, save_p) = (
            ct3["p"].as_f64().unwrap(),
            f64::from(saved.row(*id).unwrap()[j]),
        );
        assert!(
            (preview_p - save_p).abs() < 1e-6,
            "K{id}: preview p {preview_p} vs save {save_p}"
        );
        assert!(calls.iter().all(|c| c["share"].is_number()));
    }
}

/// `enriched_round` relabelled once, so it records p-values and NES.
fn rescored_round(root: &Path) -> Loaded {
    let src = enriched_round(root);
    let d = root.join("d0.jsonl");
    fs::write(
        &d,
        line(json!({"action": "markers_drop", "label": "CT3", "features": ["GENE14"], "rationale": "r", "decided_by": "user"})),
    )
    .unwrap();
    run_relabel(&args(&src, &d, Some(&root.join("r1/run")))).unwrap();
    run::load(&root.join("r1/run.senna.json").to_string_lossy()).unwrap()
}

#[test]
fn only_the_touched_types_are_rescored_and_they_match_a_full_rescore() {
    use super::super::recalibrate::{rescore, rescore_types};
    let root = tempfile::tempdir().unwrap();
    let r1 = rescored_round(root.path());
    let cells = read_cells(&r1).unwrap();
    let (before, _) = read_markers(&r1).unwrap();
    let after: Vec<(String, String)> = before
        .iter()
        .filter(|(g, t)| !(g == "GENE13" && t == "CT3"))
        .cloned()
        .collect();
    let touched = touched_types(&before, &after);
    assert_eq!(touched.iter().collect::<Vec<_>>(), ["CT3"]);

    let part = rescore_types(&r1, &cells, &after, &touched)
        .unwrap()
        .unwrap();
    let full = rescore(&r1, &cells, &after).unwrap().unwrap();
    let col = |name: &str| part.types.iter().position(|t| &**t == name).unwrap();
    let (ct3, ct1) = (col("CT3"), col("CT1"));
    for row in 0..part.ids.len() {
        assert_eq!(part.p_values[(row, ct3)], full.p_values[(row, ct3)]);
        assert_eq!(part.nes[(row, ct3)], full.nes[(row, ct3)]);
        assert!(part.z[(row, ct3)].is_finite(), "CT3 rescored");
        assert!(part.z[(row, ct1)].is_nan(), "CT1 kept, not rescored");
        assert!((part.p_values[(row, ct1)] - full.p_values[(row, ct1)]).abs() < 1e-6);
    }

    // A new type moves every type's IDF weight: all are rescored.
    let mut grown = after.clone();
    for g in ["GENE20", "GENE21", "GENE22"] {
        grown.push((g.into(), "CT9".into()));
    }
    let all = rescore_types(&r1, &cells, &grown, &touched_types(&after, &grown))
        .unwrap()
        .unwrap();
    assert!((0..all.types.len()).all(|t| all.z[(0, t)].is_finite()));
}

#[test]
fn an_edit_touches_its_type_and_the_types_sharing_its_genes() {
    let panel = |pairs: &[(&str, &str)]| -> Vec<(String, String)> {
        pairs
            .iter()
            .map(|(g, t)| ((*g).into(), (*t).into()))
            .collect()
    };
    let before = panel(&[("G1", "A"), ("G2", "A"), ("G2", "B"), ("G3", "C")]);
    let after = panel(&[("G1", "A"), ("G2", "B"), ("G3", "C"), ("G4", "New")]);
    let touched = touched_types(&before, &after);
    assert_eq!(
        touched.into_iter().collect::<Vec<_>>(),
        ["A", "B", "New"],
        "A lost G2, B shares G2, New is new; C is untouched"
    );
}

#[test]
fn a_round_needs_decisions() {
    let root = tempfile::tempdir().unwrap();
    let src = enriched_round(root.path());
    let d = root.path().join("empty.jsonl");
    fs::write(&d, "").unwrap();
    assert!(run_relabel(&args(&src, &d, Some(&root.path().join("x/run")))).is_err());
}

#[test]
fn a_coarse_label_agrees_with_a_top_call_inside_its_group() {
    let root = tempfile::tempdir().unwrap();
    let src = enriched_round(root.path());
    let dir = src.parent().unwrap();
    fs::write(
        dir.join("run.celltype_tree.json"),
        json!({"source": "marker_sharing", "groups": [
            {"name": "CT1/CT2", "members": ["CT1", "CT2"]},
            {"name": "CT3", "members": ["CT3"]}
        ]})
        .to_string(),
    )
    .unwrap();
    let mut m = RunManifest::load(&src).unwrap().0;
    m.annotate.celltype_tree = Some("run.celltype_tree.json".into());
    m.save(&src).unwrap();

    let d = root.path().join("d.jsonl");
    fs::write(&d, line(json!({"cluster": 0, "action": "label", "label": "CT1/CT2", "rationale": "r", "decided_by": "user"})))
        .unwrap();
    run_relabel(&args(&src, &d, Some(&root.path().join("r1/run")))).unwrap();
    let r1 = run::load(&root.path().join("r1/run.senna.json").to_string_lossy()).unwrap();
    let summary: Value = serde_json::from_str(
        &fs::read_to_string(resolve(
            &r1.dir,
            r1.manifest.annotate.cluster_summary.as_deref().unwrap(),
        ))
        .unwrap(),
    )
    .unwrap();
    assert_eq!(summary["0"]["evidence"]["top"], "CT1");
    assert_eq!(summary["0"]["evidence"]["agrees"], true);
}

#[test]
fn a_pasted_chat_answer_reads_as_decisions() {
    let text = "Here are my suggestions:\n```json\n{\"cluster\": 0, \"action\": \"keep\", \"rationale\": \"r\", \"decided_by\": \"agent_proposed_user_accepted\"}\n```\nHope this helps.\n";
    let ds = parse_decisions(numbered(text), "paste").unwrap();
    assert_eq!(ds.len(), 1);
    assert!(
        parse_decisions(numbered("{not json}\n"), "paste").is_err(),
        "a broken decision is still refused"
    );
}

#[test]
fn a_moved_runs_markers_are_read_from_where_the_run_now_is() {
    let tmp = tempfile::tempdir().unwrap();
    let proj = tmp.path().join("proj");
    std::fs::create_dir_all(proj.join("run")).unwrap();
    std::fs::create_dir_all(proj.join("data")).unwrap();
    std::fs::write(proj.join("data/panel.tsv"), "GENE1\tCT1\nGENE2\tCT2\n").unwrap();
    let manifest = proj.join("run/x.senna.json");
    std::fs::write(
        &manifest,
        r#"{"version":2,"kind":"topic","prefix":"/elsewhere/proj/run/x",
            "annotate":{"markers":"/elsewhere/proj/data/panel.tsv"}}"#,
    )
    .unwrap();
    let loaded = run::load(&manifest.to_string_lossy()).unwrap();
    let (pairs, _) = read_markers(&loaded).expect("the panel, found again");
    assert_eq!(pairs.len(), 2);
}
