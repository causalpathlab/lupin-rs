use super::*;
use std::time::Instant;

/// Take panels as they come until the scan ends.
fn wait(p: &mut Picker) {
    let start = Instant::now();
    while p.scan.is_some() {
        assert!(start.elapsed() < Duration::from_secs(10), "the scan ends");
        p.poll();
        std::thread::sleep(Duration::from_millis(5));
    }
}

fn names(p: &Picker) -> Vec<String> {
    p.shown().iter().map(|e| e.name.clone()).collect()
}

/// A directory with a panel, a one-type panel, a count table and a log.
fn dir() -> tempfile::TempDir {
    let root = tempfile::tempdir().unwrap();
    let at = |n: &str, s: &str| std::fs::write(root.path().join(n), s).unwrap();
    at("a_markers.tsv", "G1\tT1\nG2\tT2\nG3\tT2\n");
    at("one.tsv", "G1\tT1\n");
    at("counts.tsv", "G1\t3\nG2\t5\n");
    at("run.log", "G1\tT1\n");
    std::fs::create_dir(root.path().join("sub")).unwrap();
    root
}

#[test]
fn the_listing_comes_before_its_panels_are_read() {
    let root = dir();
    let mut p = Picker::new("t", root.path(), Want::Markers(None, 0));
    // Every file that may be a panel is listed while it is read; a log is not.
    let pending: Vec<&str> = p
        .entries
        .iter()
        .filter(|e| e.pending || e.panel.is_some())
        .map(|e| &*e.name)
        .collect();
    assert!(pending.contains(&"a_markers.tsv"));
    assert!(!names(&p).contains(&"run.log".to_string()));

    wait(&mut p);
    assert!(p.entries.iter().all(|e| !e.pending));
    let shown = names(&p);
    assert!(shown.contains(&"a_markers.tsv".to_string()));
    assert!(shown.contains(&"one.tsv".to_string()));
    assert!(!shown.contains(&"counts.tsv".to_string()), "a count table");
    assert!(shown.contains(&"sub".to_string()));
}

#[test]
fn panels_read_once_are_not_read_again() {
    let root = dir();
    let mut p = Picker::new("t", root.path(), Want::Markers(None, 0));
    wait(&mut p);
    p.open(root.path().join("sub"));
    p.open(root.path().to_path_buf());
    assert!(p.scan.is_none(), "nothing left to read");
    assert!(p.entries.iter().all(|e| !e.pending));
    assert!(names(&p).contains(&"a_markers.tsv".to_string()));
}

#[test]
fn the_best_panel_is_chosen_once_read_and_a_moved_selection_stays() {
    let root = dir();
    let genes = ["G1", "G2", "G3"].map(Box::<str>::from);
    let index = Arc::new(GeneRows::build(&genes));
    let mut p = Picker::new("t", root.path(), Want::Markers(Some(index), 3));
    wait(&mut p);
    assert_eq!(
        p.best.as_deref().and_then(Path::file_name),
        Some("a_markers.tsv".as_ref())
    );
    assert_eq!(p.selected_path(), p.best.clone());

    // Moved by hand, a later panel does not take the selection back.
    p.key(KeyCode::Up);
    let at = p.selected_path();
    p.settle(at.clone());
    assert_eq!(p.selected_path(), at);
}

#[test]
fn enter_on_a_file_not_read_yet_reads_it_then() {
    let root = dir();
    let mut p = Picker::new("t", root.path(), Want::Markers(None, 0));
    // Stop the scan: what is chosen is read on the key.
    p.scan = None;
    let one = p.shown().iter().position(|e| e.name == "one.tsv").unwrap();
    p.state.select(Some(one));
    assert!(matches!(p.key(KeyCode::Enter), Step::Stay));
    assert!(p.note.as_deref().unwrap().contains("one cell type"));
    let on = p.selected_path();
    assert_eq!(
        on.as_deref().and_then(Path::file_name),
        Some("one.tsv".as_ref()),
        "the note's file"
    );

    let counts = p
        .shown()
        .iter()
        .position(|e| e.name == "counts.tsv")
        .unwrap();
    p.state.select(Some(counts));
    assert!(matches!(p.key(KeyCode::Enter), Step::Stay));
    assert_eq!(p.note.as_deref(), Some("not a marker panel"));

    let panel = p
        .shown()
        .iter()
        .position(|e| e.name == "a_markers.tsv")
        .unwrap();
    p.state.select(Some(panel));
    assert!(matches!(p.key(KeyCode::Enter), Step::Chosen(_)));
}

#[test]
fn a_scan_that_ends_early_leaves_nothing_reading() {
    let root = dir();
    let mut p = Picker::new("t", root.path(), Want::Markers(None, 0));
    // The reader gone before it read a thing, as after a panic.
    let (tx, rx) = mpsc::channel();
    drop(tx);
    p.scan = Some(rx);
    p.poll();
    assert!(p.scan.is_none());
    assert!(p.entries.iter().all(|e| !e.pending), "none still reading");
    assert!(!p.hint("q").0.contains("reading"), "{}", p.hint("q").0);
}
