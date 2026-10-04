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

    // A statement that closes a cycle is kept in the file but the view says so.
    v.record("C", "A", Relation::Precedes, "a mistake").unwrap();
    let v = OrderView::load(types(), None, &search, Vec::new());
    assert!(v.prior.is_none());
    assert!(
        v.error.as_deref().is_some_and(|e| e.contains("cycle")),
        "{:?}",
        v.error
    );

    // `unrelated` withdraws it again.
    v.record("A", "C", Relation::Unrelated, "undone").unwrap();
    let v = OrderView::load(types(), None, &search, Vec::new());
    assert_eq!(v.edges().len(), 2);
    assert_eq!(v.unrelated().len(), 1);
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
        "HSC",
        "--root",
        "MPP",
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
