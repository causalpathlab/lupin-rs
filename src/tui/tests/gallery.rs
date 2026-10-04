use super::*;

fn figure_set(dir: &Path, base: &str) -> Vec<PathBuf> {
    let svg = dir.join(format!("{base}.svg"));
    let pdf = dir.join(format!("{base}.pdf"));
    std::fs::write(&svg, "<svg/>").unwrap();
    std::fs::write(&pdf, "%PDF-1.4 fake").unwrap();
    vec![svg, pdf]
}

#[test]
fn an_export_is_logged_checked_moved_and_removed() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().canonicalize().unwrap();
    let dir = dir.as_path();
    let log_dir = dir.join(DIR);
    let files = figure_set(dir, "run.trajectory.order");
    let pdf = files[1].clone();
    let mut g = Gallery::open(&log_dir);
    g.add(&files, "order · 7 in", "order", &dir.join("run.senna.json"))
        .unwrap();
    assert_eq!(g.entries.len(), 1);
    assert_eq!(g.entries[0].files.len(), 2);
    assert_eq!(g.check(), vec![Status::Intact]);

    // Same bytes under a new mtime: still ok; other bytes of the same size,
    // or another size: changed.
    std::fs::write(&files[1], "%PDF-1.4 fake").unwrap();
    assert_eq!(Gallery::open(&log_dir).check(), vec![Status::Intact]);
    let mut g2 = Gallery::open(&log_dir);
    g2.entries[0].files[1].mtime = 1;
    std::fs::write(&files[1], "%PDF-1.4 faKe").unwrap();
    assert_eq!(g2.check(), vec![Status::Changed]);
    std::fs::write(&files[1], "%PDF-1.4 other").unwrap();
    assert_eq!(Gallery::open(&log_dir).check(), vec![Status::Changed]);

    // A second set beside the manifest, not logged, is noticed.
    figure_set(dir, "run.trajectory.layout");
    let prefix = dir.join("run").to_string_lossy().into_owned();
    assert_eq!(g.not_listed(&prefix).len(), 1);

    // A changed file is neither moved nor deleted on the log's word.
    let mut g = Gallery::open(&log_dir);
    assert!(g.relocate(&pdf, &dir.join("figure1")).is_err());
    assert!(g.remove(&pdf, true).is_err());
    std::fs::write(&files[1], "%PDF-1.4 fake").unwrap();

    // Moving renames every file of the set and refuses to replace.
    let mut g = Gallery::open(&log_dir);
    g.relocate(&pdf, &dir.join("figure1")).unwrap();
    assert!(dir.join("figure1.svg").is_file() && dir.join("figure1.pdf").is_file());
    assert!(!files[0].exists());
    figure_set(dir, "taken");
    assert!(g
        .relocate(&dir.join("figure1.pdf"), &dir.join("taken"))
        .is_err());

    // A missing file is reported, not dropped.
    std::fs::remove_file(dir.join("figure1.pdf")).unwrap();
    let mut g = Gallery::open(&log_dir);
    assert_eq!(g.check(), vec![Status::Missing]);
    assert_eq!(g.entries.len(), 1);
    assert!(
        g.remove(&dir.join("figure1.pdf"), true).is_err(),
        "a missing set is not deleted blindly"
    );
    g.remove(&dir.join("figure1.pdf"), false).unwrap();
    assert!(g.entries.is_empty() && dir.join("figure1.svg").exists());
}

#[test]
fn an_unreadable_log_is_set_aside_not_emptied() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let log_dir = dir.join(DIR);
    std::fs::create_dir_all(&log_dir).unwrap();
    std::fs::write(log_dir.join("saved.json"), "{ not json").unwrap();
    let g = Gallery::open(&log_dir);
    assert!(g.entries.is_empty());
    assert!(g.set_aside.as_deref().is_some_and(|p| p.is_file()));
    assert!(!log_dir.join("saved.json").exists());
}
