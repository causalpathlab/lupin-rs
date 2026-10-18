//! Suggested Cell Ontology terms for a panel type with none: offered with
//! their scores on `c`, written to the user's alias table only once
//! confirmed, and the tree rebuilt with the alias.

use super::super::app::{App, Focus, Job, Pending, Running};
use super::super::drive::{finish, press};
use super::*;
use crate::annotate::celltype_tree::ClTerms;
use crate::annotate::cl_rules::Aliases;
use crate::annotate::term_match::{Suggested, Suggestion};
use crate::annotate_cmd::AnnotateCliArgs;
use crate::manifest::data_files::SearchPath;
use ratatui::crossterm::event::KeyCode;
use std::path::{Path, PathBuf};

/// root ─┬─ group1 (class) ─┬─ CT1 cell
///       │                 └─ CT2 cell
///       └─ CT3 cell (not on the panel)
const OBO: &str = "format-version: 1.2

[Term]
id: CL:0
name: root cell

[Term]
id: CL:1
name: group1
subset: cellxgene_subset
is_a: CL:0

[Term]
id: CL:10
name: CT1 cell
is_a: CL:1

[Term]
id: CL:12
name: CT2 cell
is_a: CL:1

[Term]
id: CL:2
name: CT3 cell
is_a: CL:0
";

/// The panel: two types on terms and `XQ`, which matches none.
const PANEL: &str = "G1\tCT1 cells\nG2\tCT2 cells\nG3\tXQ\n";

/// `XQ`'s suggestions as `lupin data suggest-cl` writes them, in `dir`.
fn suggested(dir: &Path) -> PathBuf {
    let s = |id: &str, name: &str, score: f32| Suggestion {
        id: id.into(),
        name: name.into(),
        score,
    };
    let out = Suggested {
        labels: [(
            "XQ".to_string(),
            vec![
                s("CL:2", "CT3 cell", 0.91),
                s("CL:12", "CT2 cell", 0.62),
                s("CL:10", "CT1 cell", 0.58),
            ],
        )]
        .into(),
        ..Suggested::default()
    };
    let file = dir.join("suggested.json");
    std::fs::write(&file, serde_json::to_string(&out).unwrap()).unwrap();
    file
}

/// An app on the panel in `dir`, the user's config there too, its tree
/// focused on `XQ`.
fn app(dir: &Path) -> App {
    #[derive(clap::Parser)]
    struct Cli {
        #[command(flatten)]
        annotate: AnnotateCliArgs,
    }
    let panel = dir.join("panel.tsv");
    std::fs::write(&panel, PANEL).unwrap();
    let argv = ["lupin", "-o", "x", "-m", panel.to_str().unwrap()];
    let args = <Cli as clap::Parser>::parse_from(argv).annotate;
    let cl = ClTerms::parse(OBO, &crate::annotate::cl_rules::shipped());
    let types = crate::annotate::markers::read_panel(panel.to_str().unwrap()).unwrap();
    let mut app = App::new(args, PathBuf::new(), PathBuf::new(), types, Some(cl));
    app.data_search = SearchPath {
        user: Some(dir.join("config")),
        ..SearchPath::new(None)
    };
    app.focus = Focus::Tree;
    select(&mut app, "XQ");
    app
}

fn select(app: &mut App, label: &str) {
    let tree = &app.tree;
    app.tree_sel = tree
        .visible()
        .iter()
        .position(|&i| tree.label(i) == label)
        .expect("the label is on the tree");
}

#[test]
fn the_suggestions_file_reads_back_by_label() {
    let dir = tempfile::tempdir().unwrap();
    let s = Suggestions::read(&suggested(dir.path())).unwrap();
    let xq = s.of("XQ").unwrap();
    assert_eq!(xq.len(), 3);
    assert_eq!(xq[0].id, "CL:2");
    assert!(s.of("CT1_cells").is_none());
}

#[test]
fn the_types_without_a_term_are_the_ones_to_suggest_for() {
    let dir = tempfile::tempdir().unwrap();
    let app = app(dir.path());
    let cl = app.cl.as_ref().unwrap();
    assert_eq!(unmatched(cl, &app.tree), ["XQ"]);
    for t in ["CT1 cells", "CT2 cells", "XQ"] {
        assert_eq!(
            needs_term(cl, &app.tree, t),
            unmatched(cl, &app.tree).contains(&label_key(t)),
            "{t}: c and the job ask one question"
        );
    }
}

#[test]
fn a_type_with_its_own_term_needs_none_on_a_tree_built_from_markers() {
    // One type on a term is too few for the ontology's tree: the tree is
    // the markers' and puts no type on a term, yet CT1 has its own.
    let dir = tempfile::tempdir().unwrap();
    let panel = dir.path().join("few.tsv");
    std::fs::write(&panel, "G1\tCT1 cells\nG2\tXQ\nG3\tYZ\n").unwrap();
    let mut app = app(dir.path());
    app.panel = crate::annotate::markers::read_panel(panel.to_str().unwrap()).unwrap();
    app.rebuild_tree();
    assert!(app.tree.cl_of("CT1_cells").is_none(), "the markers' tree");
    let cl = app.cl.as_ref().unwrap();
    assert!(!needs_term(cl, &app.tree, "CT1_cells"));
    assert_eq!(unmatched(cl, &app.tree), ["XQ", "YZ"]);
    select(&mut app, "CT1_cells");
    app.suggestions = Some(Suggestions::read(&suggested(dir.path())).unwrap());
    press(&mut app, KeyCode::Char('c'));
    assert!(app.menu.is_none());
    assert!(app.status.contains("has its own"), "{}", app.status);
}

#[test]
fn c_offers_the_top_three_with_scores_and_a_confirmed_one_is_written_and_placed() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    assert!(app.tree.cl_of("XQ").is_none());
    app.suggestions = Some(Suggestions::read(&suggested(dir.path())).unwrap());
    press(&mut app, KeyCode::Char('c'));
    let m = app.menu.as_ref().expect("the suggestions");
    assert_eq!(m.choices.len(), 4, "three terms and none of them");
    assert!(
        m.choices[0].label.contains("CT3 cell"),
        "{}",
        m.choices[0].label
    );
    assert!(
        m.choices[0].label.contains("0.91"),
        "{}",
        m.choices[0].label
    );
    assert!(
        m.choices[0].label.contains("CL:2"),
        "{}",
        m.choices[0].label
    );
    // Taking one asks to confirm before anything is written.
    press(&mut app, KeyCode::Char('1'));
    let file = dir.path().join("config").join("cl_aliases.tsv");
    assert!(matches!(
        app.prompt.as_ref().map(|p| &p.pending),
        Some(Pending::Remember { id, file: f, .. }) if id == "CL:2" && *f == file
    ));
    assert!(!file.exists(), "nothing is written before the confirm");
    press(&mut app, KeyCode::Enter);
    let text = std::fs::read_to_string(&file).unwrap();
    let mut a = Aliases::default();
    assert_eq!(a.add_tsv(&text, "user"), 1);
    assert_eq!(a.get("XQ"), Some("CL:2"));
    // The tree is built again with the alias: XQ now sits on its term.
    assert_eq!(app.tree.cl_of("XQ"), Some("CL:2"));
    assert!(unmatched(app.cl.as_ref().unwrap(), &app.tree).is_empty());
}

#[test]
fn esc_at_the_confirm_or_none_of_them_writes_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("config").join("cl_aliases.tsv");
    let mut app = app(dir.path());
    app.suggestions = Some(Suggestions::read(&suggested(dir.path())).unwrap());
    press(&mut app, KeyCode::Char('c'));
    press(&mut app, KeyCode::Char('n'));
    assert!(app.menu.is_none() && app.prompt.is_none());
    assert!(app.status.contains("none"), "{}", app.status);
    press(&mut app, KeyCode::Char('c'));
    press(&mut app, KeyCode::Char('2'));
    press(&mut app, KeyCode::Esc);
    assert!(!file.exists());
    assert!(app.tree.cl_of("XQ").is_none());
}

#[test]
fn c_on_a_type_with_a_term_or_with_no_ontology_says_why_not() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    select(&mut app, "CT1_cells");
    press(&mut app, KeyCode::Char('c'));
    assert!(app.menu.is_none());
    assert!(app.status.contains("has its own"), "{}", app.status);
    let mut app = self::app(dir.path());
    app.cl = None;
    press(&mut app, KeyCode::Char('c'));
    assert!(app.status.contains("no Cell Ontology"), "{}", app.status);
}

#[test]
fn c_waits_for_a_running_job_before_ranking() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    let child = std::process::Command::new("sleep")
        .arg("30")
        .spawn()
        .unwrap();
    app.child = Some(Running::new(child, Job::Pass));
    press(&mut app, KeyCode::Char('c'));
    assert!(app.menu.is_none());
    assert!(app.status.contains("wait"), "{}", app.status);
    if let Some(mut r) = app.child.take() {
        let _ = r.child.kill();
        let _ = r.child.wait();
    }
}

#[test]
fn a_ranking_job_that_ends_opens_the_suggestions_for_the_type_asked() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    let file = suggested(dir.path());
    let child = std::process::Command::new("true").spawn().unwrap();
    app.child = Some(Running::new(
        child,
        Job::Suggest {
            file: file.clone(),
            label: "XQ".into(),
        },
    ));
    finish(&mut app);
    assert!(app.suggestions.as_ref().and_then(|s| s.of("XQ")).is_some());
    let m = app.menu.as_ref().expect("XQ's suggestions");
    assert!(m.choices[0].label.contains("CT3 cell"));
}

#[test]
fn without_the_model_the_tui_says_so_and_goes_on() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    let child = std::process::Command::new("false").spawn().unwrap();
    app.child = Some(Running::new(
        child,
        Job::Suggest {
            file: dir.path().join("none.json"),
            label: "XQ".into(),
        },
    ));
    finish(&mut app);
    assert!(app.failed.is_none(), "no failure popup: the rest works");
    assert!(app.menu.is_none());
    assert!(app.status.contains("suggestions"), "{}", app.status);
    assert!(app.status.contains("offline"), "{}", app.status);
}

#[test]
fn a_jobs_file_is_in_a_private_folder_of_the_apps_own() {
    let dir = tempfile::tempdir().unwrap();
    let mut a = app(dir.path());
    let mut b = app(dir.path());
    let fa = a.job_file("suggest.json").unwrap();
    let fb = b.job_file("suggest.json").unwrap();
    let folder = fa.parent().unwrap().to_path_buf();
    assert_ne!(
        folder,
        std::env::temp_dir(),
        "no name others can guess in the shared folder"
    );
    assert_ne!(fa, fb);
    assert_eq!(
        a.job_file("rescore.jsonl").unwrap().parent(),
        Some(folder.as_path())
    );
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mode = std::fs::metadata(&folder).unwrap().permissions().mode();
        assert_eq!(mode & 0o777, 0o700, "only the user reads or writes in it");
    }
    drop(a);
    assert!(!folder.exists(), "gone with the app");
}

#[test]
fn a_label_cl_file_taken_in_the_tui_rebuilds_the_tree_on_it() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    let file = dir.path().join("labels.tsv");
    std::fs::write(&file, "XQ\tCL:2\n").unwrap();
    let group = app.tree.node_of("group1").unwrap();
    app.tree.fold(group, true);
    select(&mut app, "XQ");
    app.set_label_cl(Some(&file));
    assert_eq!(app.tree.cl_of("XQ"), Some("CL:2"));
    assert_eq!(app.term_of("XQ").as_deref(), Some("CL:2"));
    let group = app.tree.node_of("group1").unwrap();
    assert!(app.tree.is_folded(group), "what was folded stays folded");
    let at = app.tree.visible()[app.tree_sel];
    assert_eq!(app.tree.label(at), "XQ", "the selection stays on its node");
}

#[test]
fn a_rebuild_on_a_label_no_longer_on_the_tree_keeps_the_selection_it_had() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    let before: Vec<String> = (app.tree.visible().iter())
        .map(|&i| app.tree.label(i).to_string())
        .collect();
    select(&mut app, "XQ");
    // A type gone from the panel moves the rows from where they were.
    app.panel.retain(|(_, t)| t != "CT1_cells");
    app.rebuild_tree_keeping(Some("gone"));
    let after: Vec<String> = (app.tree.visible().iter())
        .map(|&i| app.tree.label(i).to_string())
        .collect();
    assert_ne!(before, after, "{before:?}");
    let at = app.tree.visible()[app.tree_sel];
    assert_eq!(app.tree.label(at), "XQ", "{after:?}");
}

#[test]
fn a_rebuild_leaves_the_selection_on_a_row_of_the_new_tree() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    let last = app.tree.visible().len() - 1;
    app.tree_sel = last;
    app.panel.truncate(1);
    app.rebuild_tree();
    assert!(
        app.tree_sel < app.tree.visible().len(),
        "{} of {}",
        app.tree_sel,
        app.tree.visible().len()
    );
}
