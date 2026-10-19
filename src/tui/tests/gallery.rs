use super::*;
use legume_tui::gallery::LOG;
use sha2::{Digest, Sha256};

fn figure_set(dir: &Path, base: &str) -> Vec<PathBuf> {
    let svg = dir.join(format!("{base}.svg"));
    let pdf = dir.join(format!("{base}.pdf"));
    std::fs::write(&svg, "<svg/>").unwrap();
    std::fs::write(&pdf, "%PDF-1.4 fake").unwrap();
    vec![svg, pdf]
}

fn trace(dir: &Path, panel: &str) -> Trace {
    Trace::of(panel, &dir.join("run.senna.json"))
}

#[test]
fn an_export_is_logged_checked_moved_and_removed() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().canonicalize().unwrap();
    let dir = dir.as_path();
    let log_dir = dir.join(DIR);
    let files = figure_set(dir, "run.trajectory.order");
    let pdf = files[1].clone();
    let mut g = open(&log_dir);
    g.add(&files, "order · 7 in", None, trace(dir, "order"))
        .unwrap();
    assert_eq!(g.entries().len(), 1);
    assert_eq!(g.entries()[0].files.len(), 2);
    assert_eq!(g.entries()[0].extra.panel, "order");
    assert_eq!(g.check(), vec![Status::Intact]);

    // Same bytes under a new mtime: still ok; another size: changed.
    std::fs::write(&files[1], "%PDF-1.4 fake").unwrap();
    assert_eq!(open(&log_dir).check(), vec![Status::Intact]);
    std::fs::write(&files[1], "%PDF-1.4 other").unwrap();
    assert_eq!(open(&log_dir).check(), vec![Status::Changed]);

    // A second set beside the manifest, not logged, is noticed.
    figure_set(dir, "run.trajectory.layout");
    let prefix = dir.join("run.trajectory").to_string_lossy().into_owned();
    assert_eq!(not_listed(&g, &prefix).len(), 1);

    // A changed file is neither moved nor deleted on the log's word.
    let mut g = open(&log_dir);
    assert!(g.relocate(&pdf, &dir.join("figure1")).is_err());
    assert!(g.remove(&pdf, true).is_err());
    std::fs::write(&files[1], "%PDF-1.4 fake").unwrap();

    // Moving renames every file of the set and refuses to replace.
    let mut g = open(&log_dir);
    g.relocate(&pdf, &dir.join("figure1")).unwrap();
    assert!(dir.join("figure1.svg").is_file() && dir.join("figure1.pdf").is_file());
    assert!(!files[0].exists());
    figure_set(dir, "taken");
    assert!(g
        .relocate(&dir.join("figure1.pdf"), &dir.join("taken"))
        .is_err());

    // A missing file is reported, not dropped.
    std::fs::remove_file(dir.join("figure1.pdf")).unwrap();
    let mut g = open(&log_dir);
    assert_eq!(g.check(), vec![Status::Missing]);
    assert_eq!(g.entries().len(), 1);
    assert!(
        g.remove(&dir.join("figure1.pdf"), true).is_err(),
        "a missing set is not deleted blindly"
    );
    g.remove(&dir.join("figure1.pdf"), false).unwrap();
    assert!(g.entries().is_empty() && dir.join("figure1.svg").exists());
}

#[test]
fn a_log_lupin_wrote_before_loads_and_still_moves_its_files() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().canonicalize().unwrap();
    let files = figure_set(&dir, "run.trajectory.order");
    let record = |p: &Path| {
        let bytes = std::fs::read(p).unwrap();
        serde_json::json!({
            "path": p,
            "size": bytes.len(),
            "mtime": 0,
            "hash": crate::hex(&Sha256::digest(&bytes)),
        })
    };
    // As lupin's own gallery wrote it: no thumbnail, its trace inline.
    let old = serde_json::json!([{
        "path": files[1],
        "what": "order · 7 in · 200 dpi",
        "when": 1_700_000_000u64,
        "panel": "order",
        "manifest": dir.join("run.senna.json"),
        "files": [record(&files[0]), record(&files[1])],
    }]);
    let log_dir = dir.join(DIR);
    std::fs::create_dir(&log_dir).unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let private = std::fs::Permissions::from_mode(0o700);
        std::fs::set_permissions(&log_dir, private).unwrap();
    }
    std::fs::write(log_dir.join(LOG), old.to_string()).unwrap();

    let mut g = open(&log_dir);
    assert!(g.set_aside().is_none(), "read, not set aside");
    assert_eq!(g.entries().len(), 1);
    let e = &g.entries()[0];
    assert_eq!(e.extra, trace(&dir, "order"));
    assert_eq!((e.thumb.as_deref(), e.files.len()), (None, 2));
    assert_eq!(g.check(), vec![Status::Intact]);
    // Its hashes are the ones checked before a move.
    g.relocate(&files[1], &dir.join("figure1")).unwrap();
    assert!(dir.join("figure1.svg").is_file() && dir.join("figure1.pdf").is_file());
}

#[test]
fn an_unreadable_log_is_set_aside_not_emptied() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let log_dir = dir.join(DIR);
    std::fs::create_dir_all(&log_dir).unwrap();
    std::fs::write(log_dir.join(LOG), "{ not json").unwrap();
    let g = open(&log_dir);
    assert!(g.entries().is_empty());
    assert!(g.set_aside().is_some_and(Path::is_file));
    assert!(!log_dir.join(LOG).exists());
}

#[cfg(unix)]
#[test]
fn a_log_others_can_write_moves_and_deletes_nothing() {
    use std::os::unix::fs::PermissionsExt;
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().canonicalize().unwrap();
    let log_dir = dir.join(DIR);
    let files = figure_set(&dir, "run.trajectory.order");
    let mut g = open(&log_dir);
    g.add(&files, "order", None, trace(&dir, "order")).unwrap();
    let mode = |p: &Path| std::fs::metadata(p).unwrap().permissions().mode() & 0o777;
    assert_eq!(
        mode(&log_dir),
        0o700,
        "the log's directory is created private"
    );
    assert_eq!(mode(&log_dir.join(LOG)) & 0o077, 0, "so is the log");

    let anyone = std::fs::Permissions::from_mode(0o666);
    std::fs::set_permissions(log_dir.join(LOG), anyone).unwrap();
    let mut g = open(&log_dir);
    assert!(g.remove(&files[1], true).is_err());
    assert!(g.relocate(&files[1], &dir.join("elsewhere")).is_err());
    assert!(files.iter().all(|f| f.is_file()));
    // Nor is an entry taken off a log others can write, nor is the log
    // taken over by saving a new export over it.
    assert!(g.remove(&files[1], false).is_err());
    let more = figure_set(&dir, "run.trajectory.layout");
    assert!(g.add(&more, "layout", None, trace(&dir, "layout")).is_err());
    assert_eq!(mode(&log_dir.join(LOG)), 0o666, "the log is left as it was");
}

#[cfg(unix)]
#[test]
fn a_log_that_is_a_link_moves_and_deletes_nothing() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().canonicalize().unwrap();
    let files = figure_set(&dir, "run.trajectory.order");
    // A genuine log of this user's elsewhere, linked in as this one.
    let other = dir.join("elsewhere");
    let mut g = open(&other);
    g.add(&files, "order", None, trace(&dir, "order")).unwrap();
    let log_dir = dir.join(DIR);
    std::fs::create_dir(&log_dir).unwrap();
    std::os::unix::fs::symlink(other.join(LOG), log_dir.join(LOG)).unwrap();
    let mut g = open(&log_dir);
    assert!(g.remove(&files[1], true).is_err());
    assert!(files.iter().all(|f| f.is_file()));
}
