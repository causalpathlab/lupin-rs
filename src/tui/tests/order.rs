use super::*;
use crate::manifest::data_files::SearchPath;

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("lupin-order-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(dir.join("lupin")).unwrap();
    dir
}

fn types() -> Vec<(String, usize)> {
    vec![("A".into(), 10), ("B".into(), 5), ("C".into(), 3)]
}

#[test]
fn the_view_shows_the_project_statements_as_edges_and_records_new_ones() {
    let dir = scratch("edges");
    std::fs::write(
        dir.join("lupin").join(PRECEDENCE),
        "A\tB\tprecedes\tknown\n",
    )
    .unwrap();
    let search = SearchPath::new(Some(&dir));
    let v = OrderView::load(types(), None, &search, &[]);
    assert_eq!(
        v.file.as_deref(),
        Some(dir.join("lupin").join(PRECEDENCE).as_path())
    );
    let edges: Vec<(String, String, String)> = v
        .edges()
        .into_iter()
        .map(|(a, b, s, _)| (a, b, s))
        .collect();
    assert_eq!(edges, vec![("A".into(), "B".into(), "project".into())]);
    assert!(v.verdicts.is_empty() && v.error.is_none());

    v.record("B", "C", Relation::Precedes, "stated in the TUI")
        .unwrap();
    let v = OrderView::load(types(), None, &search, &[]);
    assert_eq!(v.edges().len(), 2);
    assert_eq!(v.edges()[1].1, "C");

    // A statement that closes a cycle is kept in the file but the view says so.
    v.record("C", "A", Relation::Precedes, "a mistake").unwrap();
    let v = OrderView::load(types(), None, &search, &[]);
    assert!(v.prior.is_none());
    assert!(
        v.error.as_deref().is_some_and(|e| e.contains("cycle")),
        "{:?}",
        v.error
    );

    // `unrelated` withdraws it again.
    v.record("A", "C", Relation::Unrelated, "undone").unwrap();
    let v = OrderView::load(types(), None, &search, &[]);
    assert_eq!(v.edges().len(), 2);
    assert_eq!(v.unrelated().len(), 1);
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn without_a_project_directory_there_is_nowhere_to_write() {
    let search = SearchPath::new(None);
    let v = OrderView::load(types(), None, &search, &[]);
    assert!(v.file.is_none());
    assert!(v.record("A", "B", Relation::Precedes, "x").is_err());
}
