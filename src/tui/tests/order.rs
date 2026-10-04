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
    let v = OrderView::load(types(), None, &search, Vec::new());
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
    let v = OrderView::load(types(), None, &search, Vec::new());
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
    let v = OrderView::load(types(), None, &search, Vec::new());
    assert!(v.prior.is_none());
    assert!(
        v.error.as_deref().is_some_and(|e| e.contains("cycle")),
        "{:?}",
        v.error
    );

    // A statement about the same pair replaces it.
    v.record("A", "C", Relation::Precedes, "undone").unwrap();
    let v = OrderView::load(types(), None, &search, Vec::new());
    assert_eq!(v.edges().len(), 2);
    // `unrelated` on a pair the rest of the prior orders is refused.
    let err = v.record("A", "C", Relation::Unrelated, "x").unwrap_err();
    assert!(err.to_string().contains("not recorded"), "{err}");
}

#[test]
fn without_a_project_directory_there_is_nowhere_to_write() {
    let search = SearchPath::new(None);
    let v = OrderView::load(types(), None, &search, Vec::new());
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
