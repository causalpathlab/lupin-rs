use super::*;
use crate::manifest::data_files::SearchPath;

fn types() -> Vec<(String, usize)> {
    vec![("A".into(), 10), ("B".into(), 5), ("C".into(), 3)]
}

#[test]
fn the_view_shows_the_project_statements_as_edges_and_records_new_ones() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    std::fs::create_dir_all(dir.join("lupin")).unwrap();
    std::fs::write(
        dir.join("lupin").join(PRECEDENCE),
        "A\tB\tprecedes\tknown\n",
    )
    .unwrap();
    let search = SearchPath::new(Some(dir));
    let v = OrderView::load(types(), None, &search, Vec::new(), 0);
    assert_eq!(
        v.file.as_deref(),
        Some(dir.join("lupin").join(PRECEDENCE).as_path())
    );
    let edges: Vec<(String, String, String)> = v
        .edges()
        .into_iter()
        .map(|e| (e.from, e.to, e.source))
        .collect();
    assert_eq!(edges, vec![("A".into(), "B".into(), "project".into())]);
    assert!(v.edges.is_empty() && v.error.is_none());

    v.record("B", "C", Relation::Precedes, "stated in the TUI")
        .unwrap();
    let v = OrderView::load(types(), None, &search, Vec::new(), 0);
    assert_eq!(v.edges().len(), 2);
    assert_eq!(v.edges()[1].to, "C");

    // A statement that would close a cycle is refused, and a repeated one
    // is not written twice.
    let file = v.file.clone().unwrap();
    let before = std::fs::read_to_string(&file).unwrap();
    let err = v
        .record("C", "A", Relation::Precedes, "a mistake")
        .unwrap_err();
    assert!(err.to_string().contains("cycle"), "{err}");
    let err = v.record("B", "C", Relation::Precedes, "again").unwrap_err();
    assert!(err.to_string().contains("already says"), "{err}");
    assert_eq!(std::fs::read_to_string(&file).unwrap(), before);

    // A cycle written by hand is shown as one.
    std::fs::write(&file, format!("{before}C\tA\tprecedes\ta mistake\n")).unwrap();
    let v = OrderView::load(types(), None, &search, Vec::new(), 0);
    assert!(v.prior.is_none());
    assert!(
        v.error.as_deref().is_some_and(|e| e.contains("cycle")),
        "{:?}",
        v.error
    );

    // A statement about the same pair replaces it.
    v.record("A", "C", Relation::Precedes, "undone").unwrap();
    let v = OrderView::load(types(), None, &search, Vec::new(), 0);
    assert_eq!(v.edges().len(), 2);
    // `unrelated` on a pair the rest of the prior orders is refused.
    let err = v.record("A", "C", Relation::Unrelated, "x").unwrap_err();
    assert!(err.to_string().contains("not recorded"), "{err}");
}

#[test]
fn without_a_project_directory_there_is_nowhere_to_write() {
    let search = SearchPath::new(None);
    let v = OrderView::load(types(), None, &search, Vec::new(), 0);
    assert!(v.file.is_none());
    assert!(v.record("A", "B", Relation::Precedes, "x").is_err());
}

#[test]
fn a_trajectory_run_is_offered_the_next_free_prefix() {
    let tmp = tempfile::tempdir().unwrap();
    let source = tmp.path().join("run.senna.json");
    std::fs::write(&source, "{}").unwrap();
    let stem = tmp.path().join("run").to_string_lossy().into_owned();
    assert_eq!(default_out(&source), format!("{stem}.T1"));
    std::fs::write(format!("{stem}.T1.senna.json"), "{}").unwrap();
    assert_eq!(default_out(&source), format!("{stem}.T2"));
}

#[test]
fn the_tui_carries_the_run_options_over_to_its_child() {
    use crate::trajectory::run::TrajectoryArgs;
    use clap::Parser;
    #[derive(Parser)]
    struct Cli {
        #[command(flatten)]
        t: TrajectoryArgs,
    }
    let a = Cli::parse_from([
        "lupin",
        "--root",
        "CT1",
        "--root",
        "CT2",
        "--labels",
        "l.tsv",
        "--prior",
        "p.tsv",
        "--prior-only",
        "--knn",
        "20",
        "--n-dcs",
        "8",
        "--min-cells",
        "5",
        "--min-connectivity",
        "0.2",
        "--check-only",
    ])
    .t;
    let mut argv = vec![
        "lupin".to_string(),
        "-f".into(),
        "r".into(),
        "-o".into(),
        "o".into(),
    ];
    argv.extend(a.child_argv());
    let b = Cli::parse_from(argv).t;
    assert_eq!(b.from.as_deref(), Some("r"));
    assert_eq!(b.out.as_deref(), Some("o"));
    let (a, b) = (
        TrajectoryArgs {
            from: None,
            out: None,
            ..a
        },
        TrajectoryArgs {
            from: None,
            out: None,
            ..b
        },
    );
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
}

#[test]
fn a_rerun_reuses_the_labels_file_the_last_trajectory_recorded() {
    use crate::tui::app::recorded_labels;
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let manifest = dir.join("run.T1.senna.json");
    std::fs::write(
        &manifest,
        r#"{"version":2,"kind":"topic","prefix":"run.T1",
            "trajectory":{"settings":{"labels":"labels.tsv"}}}"#,
    )
    .unwrap();
    let loaded = || crate::manifest::run::load(&manifest.to_string_lossy()).unwrap();
    assert_eq!(
        recorded_labels(&loaded()),
        None,
        "a missing file is not offered"
    );
    std::fs::write(dir.join("labels.tsv"), "c1\tCT1\n").unwrap();
    let got = recorded_labels(&loaded()).unwrap();
    assert!(std::path::Path::new(&got).ends_with("labels.tsv"), "{got}");
}

#[test]
fn a_statement_on_a_type_below_min_cells_makes_no_edge_and_is_set_aside() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    std::fs::create_dir_all(dir.join("lupin")).unwrap();
    std::fs::write(
        dir.join("lupin").join(PRECEDENCE),
        "A\tC\tprecedes\tknown\n",
    )
    .unwrap();
    let search = SearchPath::new(Some(dir));
    // C has 3 cells: below 4 it is no node, as `lupin trajectory` has it.
    let v = OrderView::load(types(), None, &search, Vec::new(), 4);
    assert!(v.edges().is_empty(), "no edge between nodes");
    let aside: Vec<(&str, &str)> = v
        .set_aside
        .iter()
        .map(|s| (s.from.as_str(), s.to.as_str()))
        .collect();
    assert_eq!(aside, [("A", "C")]);
    assert_eq!(
        OrderView::load(types(), None, &search, Vec::new(), 3)
            .edges()
            .len(),
        1
    );
    // With no cell counts (the panel's types), every type is a node.
    let unknown = vec![("A".to_string(), 0), ("C".to_string(), 0)];
    let v = OrderView::load(unknown, None, &search, Vec::new(), 4);
    assert_eq!(v.edges().len(), 1);
    assert!(v.set_aside.is_empty());
}

#[test]
fn statements_through_a_type_below_min_cells_that_order_nodes_are_not_set_aside() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    std::fs::create_dir_all(dir.join("lupin")).unwrap();
    std::fs::write(
        dir.join("lupin").join(PRECEDENCE),
        "A\tC\tprecedes\tknown\nC\tB\tprecedes\tknown\n",
    )
    .unwrap();
    let search = SearchPath::new(Some(dir));
    // C is no node, but A → C → B still orders A before B, as in the run.
    let v = OrderView::load(types(), None, &search, Vec::new(), 4);
    let edges: Vec<(String, String)> = v
        .edges()
        .iter()
        .map(|e| (e.from.clone(), e.to.clone()))
        .collect();
    assert_eq!(edges, [("A".to_string(), "B".to_string())]);
    assert!(v.set_aside.is_empty(), "{:?}", v.set_aside);
}

#[test]
fn a_run_given_no_min_cells_counts_types_as_the_run_does() {
    use super::TrajectoryRun;
    assert_eq!(
        TrajectoryRun::default().min_cells(),
        crate::trajectory::run::DEFAULT_MIN_CELLS
    );
    let given = TrajectoryRun {
        argv: vec!["--min-cells".into(), "5".into()],
        ..TrajectoryRun::default()
    };
    assert_eq!(given.min_cells(), 5);
}

// The trajectory shown in the TUI: its edges and its texts, written by
// code from the run's outputs, shown at once, and read again when they
// change.

use super::super::app::{App, Focus, Job, Running, TreeMode};
use super::super::drive::finish;
use crate::annotate_cmd::AnnotateCliArgs;
use crate::describe::trajectory::{self, tests::fixture_trajectory, OrderFacts};
use std::path::Path;

/// An app on a panel of three types in `dir`, with no round.
fn app_on(dir: &Path) -> App {
    #[derive(clap::Parser)]
    struct Cli {
        #[command(flatten)]
        annotate: AnnotateCliArgs,
    }
    let panel = dir.join("markers.tsv");
    std::fs::write(&panel, "gene\tcelltype\nG1\tCT1\nG2\tCT2\nG3\tCT3\n").unwrap();
    let argv = ["lupin", "-o", "x", "-m", panel.to_str().unwrap()];
    let args = <Cli as clap::Parser>::parse_from(argv).annotate;
    let types = crate::annotate::markers::read_panel(panel.to_str().unwrap()).unwrap();
    let tree = crate::manifest::ontology::panel_tree_on(None, &types);
    let mut app = App::new(args, PathBuf::new(), PathBuf::new(), tree);
    app.focus = Focus::Clusters;
    app
}

fn mkdir(p: PathBuf) -> PathBuf {
    std::fs::create_dir_all(&p).unwrap();
    p
}

/// The trajectory the order view shows.
fn shown(app: &App) -> Option<&ShownRun> {
    app.order_shown.as_ref()?.run.as_ref()
}

/// Trajectory `t`'s texts as code writes them from its facts now.
fn texts_of(t: &Path) -> Vec<trajectory::OrderText> {
    let loaded = crate::manifest::run::load(&t.to_string_lossy()).unwrap();
    trajectory::texts(&OrderFacts::read(&loaded).unwrap())
}

/// Trajectory `t`'s run paragraph now.
fn paragraph_of(t: &Path) -> String {
    texts_of(t)[0].text.clone()
}

/// Write trajectory `t`'s lineages as `paths` (one `a > b > c` per line).
fn relineage(t: &Path, paths: &[&str]) {
    let mut text = String::from("# lineage\tcomponent\tpath\n");
    for (i, p) in paths.iter().enumerate() {
        text += &format!("L{i}\t0\t{p}\n");
    }
    let file = t.with_file_name("t.trajectory_lineages.tsv");
    let before = std::fs::metadata(&file).and_then(|m| m.modified()).ok();
    std::fs::write(&file, text).unwrap();
    // Later than before, whatever the clock's grain.
    if let Some(before) = before {
        let later = std::time::SystemTime::now().max(before + std::time::Duration::from_secs(1));
        std::fs::File::options()
            .write(true)
            .open(&file)
            .unwrap()
            .set_modified(later)
            .unwrap();
    }
}

#[test]
fn the_shown_run_notices_a_rewrite_of_the_same_size() {
    let dir = tempfile::tempdir().unwrap();
    let t = fixture_trajectory(dir.path());
    relineage(&t, &["HSC_like > MPP_like > ERY"]);
    let mut shown = ShownRun::find(std::slice::from_ref(&t), None).0.unwrap();
    assert!(!shown.poll(), "nothing changed");
    relineage(&t, &["HSC_like > MPP_like > MON"]);
    assert!(shown.poll());
    assert_eq!(shown.texts, texts_of(&t));
}

#[test]
fn a_trajectory_found_nowhere_is_looked_for_again() {
    let dir = tempfile::tempdir().unwrap();
    let t = fixture_trajectory(dir.path());
    let edges = t.with_file_name("t.trajectory_edges.parquet");
    let away = t.with_file_name("away.parquet");
    std::fs::rename(&edges, &away).unwrap();
    let among = vec![t.clone()];
    let mut shown = Shown::find(among, None);
    assert!(shown.run.is_none(), "no edges yet");
    std::thread::sleep(POLL_EVERY);
    assert!(!shown.poll(), "nothing changed: not looked for again");
    let watched: Vec<&Path> = shown.watched.iter().map(|(f, _)| f.as_path()).collect();
    assert_eq!(
        watched,
        [t.as_path(), edges.as_path()],
        "what a look depends on"
    );
    std::fs::rename(&away, &edges).unwrap();
    std::thread::sleep(POLL_EVERY);
    assert!(shown.poll(), "written since");
    assert_eq!(shown.run.as_ref().map(|r| &r.manifest), Some(&t));
}

/// An app on a panel in `dir`, the order view up on trajectory `t`.
fn order_on(dir: &Path, t: &Path) -> App {
    let mut app = app_on(dir);
    app.source = t.to_path_buf();
    app.trajectory_tui = true;
    app.figures = None;
    app.toggle_order(false);
    app.figures = None;
    app
}

#[test]
fn the_order_view_shows_the_selected_types_line_or_the_runs_paragraph() {
    let dir = tempfile::tempdir().unwrap();
    let t = fixture_trajectory(dir.path());
    // Written from the run's outputs: there at once, no model run.
    let mut app = order_on(dir.path(), &t);
    assert!(app.child.is_none());
    app.focus = Focus::Clusters;
    let paragraph = paragraph_of(&t);
    assert_eq!(app.order_line().as_deref(), Some(paragraph.as_str()));
    app.focus = Focus::Order;
    let select = |app: &mut App, ty: &str| {
        let TreeMode::Order(v) = &mut app.tree_mode else {
            panic!("the order view")
        };
        if !v.types.iter().any(|(t, _)| t == ty) {
            v.types.push((ty.to_string(), 5));
        }
        v.sel = v.types.iter().position(|(t, _)| t == ty).unwrap();
    };
    // A type the trajectory does not order has no line: the paragraph.
    select(&mut app, "CT2");
    assert_eq!(app.order_line().as_deref(), Some(paragraph.as_str()));
    select(&mut app, "MPP_like");
    let mpp = texts_of(&t)
        .into_iter()
        .find(|x| x.subject.as_deref() == Some("MPP_like"))
        .unwrap();
    assert_eq!(app.order_line(), Some(mpp.text));
}

#[test]
fn another_trajectory_shown_shows_its_own_texts() {
    let dir = tempfile::tempdir().unwrap();
    let a = fixture_trajectory(&mkdir(dir.path().join("a")));
    let b = fixture_trajectory(&mkdir(dir.path().join("b")));
    relineage(&b, &["HSC_like > ERY"]);
    let mut app = order_on(dir.path(), &a);
    assert_eq!(shown(&app).map(|t| &t.manifest), Some(&a));
    app.source = b.clone();
    app.figures = None;
    app.tick();
    assert_eq!(shown(&app).map(|t| &t.manifest), Some(&b));
    app.focus = Focus::Clusters;
    assert_eq!(app.order_line(), Some(paragraph_of(&b)));
}

#[test]
fn opening_another_trajectory_shows_its_texts() {
    use crate::manifest::family::{Kind, Pick};
    let dir = tempfile::tempdir().unwrap();
    let a = fixture_trajectory(&mkdir(dir.path().join("a")));
    let b = fixture_trajectory(&mkdir(dir.path().join("b")));
    let mut app = order_on(dir.path(), &a);
    let pick = Pick {
        path: b.clone(),
        kind: Kind::Trajectory,
        labelled: false,
    };
    app.open_run(pick, true);
    app.tick();
    assert_eq!(shown(&app).map(|t| &t.manifest), Some(&b));
}

/// The order view's edges, as pairs.
fn order_edges(app: &App) -> Vec<(String, String)> {
    match &app.tree_mode {
        TreeMode::Order(v) => pairs(&v.edges),
        _ => panic!("the order view"),
    }
}

fn pairs(edges: &[crate::trajectory::edges::EdgeRow]) -> Vec<(String, String)> {
    (edges.iter())
        .map(|e| (e.a.to_string(), e.b.to_string()))
        .collect()
}

#[test]
fn a_new_output_shows_the_new_trajectorys_edges_with_its_texts() {
    use crate::trajectory::edges;
    let dir = tempfile::tempdir().unwrap();
    let a = fixture_trajectory(&mkdir(dir.path().join("a")));
    let b = fixture_trajectory(&mkdir(dir.path().join("b")));
    // B has fewer edges than A, and lineages of its own.
    let b_edges = b.with_file_name("t.trajectory_edges.parquet");
    let b_edges = b_edges.to_string_lossy();
    let fewer = edges::read(&b_edges).unwrap()[..2].to_vec();
    edges::write(&b_edges, &fewer).unwrap();
    let fewer = pairs(&fewer);
    relineage(&b, &["HSC_like > MPP_like"]);
    let (first, other) = (paragraph_of(&a), paragraph_of(&b));
    assert_ne!(first, other);
    let mut app = order_on(dir.path(), &a);
    app.focus = Focus::Clusters;
    app.tick();
    assert_eq!(app.order_line(), Some(first.clone()));
    assert_ne!(order_edges(&app).len(), fewer.len());
    let out = b.to_string_lossy().replace(".senna.json", "");
    app.set_output(&out);
    app.tick();
    assert_eq!(shown(&app).map(|t| &t.manifest), Some(&b));
    assert_eq!(app.order_line(), Some(other.clone()));
    assert_eq!(order_edges(&app), fewer, "the table is B's too");
    // Back to A, then to B while the view is away.
    let a_edges = pairs(
        &crate::trajectory::edges::read(
            &a.with_file_name("t.trajectory_edges.parquet")
                .to_string_lossy(),
        )
        .unwrap(),
    );
    app.set_output(&a.to_string_lossy().replace(".senna.json", ""));
    app.tick();
    assert_eq!(app.order_line(), Some(first));
    assert_eq!(order_edges(&app), a_edges);
    app.toggle_order(false);
    app.set_output(&out);
    app.toggle_order(false);
    app.figures = None;
    app.tick();
    assert_eq!(app.order_line(), Some(other));
    assert_eq!(order_edges(&app), fewer);
}

#[test]
fn the_order_view_entered_again_reads_the_texts_again() {
    let dir = tempfile::tempdir().unwrap();
    let a = fixture_trajectory(&mkdir(dir.path().join("a")));
    relineage(&a, &["HSC_like > MPP_like > ERY"]);
    let early = paragraph_of(&a);
    let mut app = order_on(dir.path(), &a);
    app.focus = Focus::Clusters;
    assert_eq!(app.order_line(), Some(early.clone()));
    app.toggle_order(false);
    // Run anew while the view is away, its lineages file the same size.
    relineage(&a, &["HSC_like > MPP_like > MON"]);
    let later = paragraph_of(&a);
    assert_ne!(later, early);
    app.toggle_order(false);
    app.focus = Focus::Clusters;
    assert_eq!(app.order_line(), Some(later));
    // So does opening the run again, the view up.
    relineage(&a, &["HSC_like > MPP_like > ERY"]);
    let pick = crate::manifest::family::Pick {
        path: a.clone(),
        kind: crate::manifest::family::Kind::Trajectory,
        labelled: false,
    };
    app.open_run(pick, true);
    assert_eq!(app.order_line(), Some(early));
}

#[test]
fn a_trajectory_run_again_on_its_output_leaves_no_stale_texts() {
    let dir = tempfile::tempdir().unwrap();
    let a = fixture_trajectory(&mkdir(dir.path().join("a")));
    let mut app = order_on(dir.path(), &a);
    // The run made again on the same output orders other types.
    relineage(&a, &["HSC_like > ERY"]);
    let child = std::process::Command::new("true").spawn().unwrap();
    app.child = Some(Running::new(child, Job::Trajectory(a.clone())));
    finish(&mut app);
    assert_eq!(shown(&app).map(|t| &t.texts), Some(&texts_of(&a)));
}

#[test]
fn the_order_views_edges_and_texts_come_from_one_trajectory() {
    use crate::trajectory::edges;
    let dir = tempfile::tempdir().unwrap();
    let a = fixture_trajectory(&mkdir(dir.path().join("a")));
    let b = fixture_trajectory(&mkdir(dir.path().join("b")));
    let b_edges = b.with_file_name("t.trajectory_edges.parquet");
    let b_edges = b_edges.to_string_lossy();
    edges::write(&b_edges, &edges::read(&b_edges).unwrap()[..2]).unwrap();
    relineage(&b, &["HSC_like > MPP_like"]);
    // A round made on A carries A's trajectory, borrowed; it is looked at
    // before B.
    let loaded = crate::manifest::run::load(&a.to_string_lossy()).unwrap();
    let borrowed = mkdir(dir.path().join("r2")).join("t.L1.senna.json");
    let copy = loaded.copy_to(borrowed.clone()).unwrap();
    copy.manifest.save(&copy.file).unwrap();
    let mut app = app_on(dir.path());
    app.target = borrowed.clone();
    app.source = b;
    app.trajectory_tui = true;
    app.figures = None;
    app.toggle_order(false);
    app.figures = None;
    app.tick();
    let a_edges = pairs(
        &edges::read(
            &a.with_file_name("t.trajectory_edges.parquet")
                .to_string_lossy(),
        )
        .unwrap(),
    );
    assert_eq!(order_edges(&app), a_edges);
    assert_eq!(shown(&app).map(|t| &t.manifest), Some(&borrowed));
    app.focus = Focus::Clusters;
    assert_eq!(app.order_line(), Some(paragraph_of(&a)));
}
