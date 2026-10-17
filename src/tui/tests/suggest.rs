//! Suggested Cell Ontology terms for a panel type with none: offered with
//! their scores on `c`, written to the user's alias table only once
//! confirmed, and the tree rebuilt with the alias.

use super::super::app::{App, Focus, Job, Pending, Running};
use super::*;
use crate::annotate::celltype_tree::ClTerms;
use crate::annotate::cl_rules::Aliases;
use crate::annotate::term_match::{Suggested, Suggestion};
use crate::annotate_cmd::AnnotateCliArgs;
use crate::manifest::data_files::SearchPath;
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
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

fn press(app: &mut App, code: KeyCode) {
    app.key(KeyEvent::new(code, KeyModifiers::NONE));
}

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
    let tree = crate::manifest::ontology::panel_tree_on(Some(&cl), &types);
    let mut app = App::new(args, PathBuf::new(), PathBuf::new(), tree);
    app.cl = Some(cl);
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
    assert_eq!(unmatched(&app.tree), ["XQ"]);
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
    assert!(unmatched(&app.tree).is_empty());
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

/// Tick until the job is taken in.
fn finish(app: &mut App) {
    for _ in 0..200 {
        app.tick();
        if app.child.is_none() {
            return;
        }
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
    panic!("the job never ended");
}

#[test]
fn a_ranking_job_that_ends_opens_the_suggestions_for_the_type_asked() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    let file = suggested(dir.path());
    let child = std::process::Command::new("true").spawn().unwrap();
    app.child = Some(Running::new(child, Job::Suggest(file.clone())));
    app.suggest_for = Some("XQ".into());
    finish(&mut app);
    assert!(app.suggestions.as_ref().and_then(|s| s.of("XQ")).is_some());
    let m = app.menu.as_ref().expect("XQ's suggestions");
    assert!(m.choices[0].label.contains("CT3 cell"));
    assert!(!file.exists(), "the job's file is read and gone");
}

#[test]
fn without_the_model_the_tui_says_so_and_goes_on() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    let child = std::process::Command::new("false").spawn().unwrap();
    app.child = Some(Running::new(
        child,
        Job::Suggest(dir.path().join("none.json")),
    ));
    app.suggest_for = Some("XQ".into());
    finish(&mut app);
    assert!(app.failed.is_none(), "no failure popup: the rest works");
    assert!(app.menu.is_none());
    assert!(app.status.contains("suggestions"), "{}", app.status);
    assert!(app.status.contains("offline"), "{}", app.status);
}
