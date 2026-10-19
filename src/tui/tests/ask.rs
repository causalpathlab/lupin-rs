//! The questions lupin asks on a panel's load: one at a time, as a menu,
//! when the panel tree is idle; esc pauses them and `Q` takes them up again.

use super::super::app::{AfterPass, App, Focus};
use super::super::drive::press;
use super::super::suggest::Suggestions;
use crate::annotate::celltype_tree::ClTerms;
use crate::annotate::term_match::{Suggested, Suggestion};
use crate::annotate_cmd::AnnotateCliArgs;
use crate::manifest::data_files::SearchPath;
use ratatui::crossterm::event::KeyCode;
use std::path::{Path, PathBuf};

/// root ─── groupA (class) ─┬─ group1 (class) ─┬─ CT1 cell
///                          │                  └─ CT2 cell
///                          └─ CT4 cell
const OBO: &str = "format-version: 1.2

[Term]
id: CL:0
name: root cell

[Term]
id: CL:5
name: groupA
subset: cellxgene_subset
is_a: CL:0

[Term]
id: CL:1
name: group1
subset: cellxgene_subset
is_a: CL:5

[Term]
id: CL:10
name: CT1 cell
is_a: CL:1

[Term]
id: CL:12
name: CT2 cell
is_a: CL:1

[Term]
id: CL:14
name: CT4 cell
is_a: CL:5
";

/// CT4 stands alone (its siblings took group1); `XQ` matches no term.
const PANEL: &str = "G1\tCT1 cells\nG2\tCT2 cells\nG4\tCT4 cells\nG3\tXQ\n";

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
    let s = Suggested {
        labels: [(
            "XQ".to_string(),
            vec![Suggestion {
                id: "CL:14".into(),
                name: "CT4 cell".into(),
                score: 0.9,
            }],
        )]
        .into(),
        ..Suggested::default()
    };
    let file = dir.join("suggested.json");
    std::fs::write(&file, serde_json::to_string(&s).unwrap()).unwrap();
    app.suggestions = Some(Suggestions::read(&file).unwrap());
    app
}

#[test]
fn a_type_left_alone_is_asked_first_and_a_join_is_written_and_applied() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    app.tick();
    let m = app.menu.as_ref().expect("a question on load");
    assert!(m.question.contains("CT4"), "{}", m.question);
    assert_eq!(m.choices.len(), 2, "join or stand alone");
    press(&mut app, KeyCode::Char('1'));
    let file = dir.path().join("config").join("cl_groups.tsv");
    let text = std::fs::read_to_string(&file).unwrap();
    assert!(text.contains("CT4_cells\tjoin"), "{text}");
    assert_eq!(app.cl.as_ref().unwrap().joined("CT4 cells"), Some(true));
}

#[test]
fn a_label_without_a_term_is_asked_next_with_its_suggestions() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    app.tick();
    press(&mut app, KeyCode::Char('2'));
    app.tick();
    let m = app.menu.as_ref().expect("the next question");
    assert!(m.question.contains("XQ"), "{}", m.question);
    assert!(m.choices[0].label.contains("CT4 cell"));
}

#[test]
fn esc_pauses_the_questions_and_q_takes_them_up_again() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    app.tick();
    press(&mut app, KeyCode::Esc);
    app.tick();
    assert!(app.menu.is_none(), "paused");
    press(&mut app, KeyCode::Char('Q'));
    app.tick();
    assert!(app
        .menu
        .as_ref()
        .is_some_and(|m| m.question.contains("CT4")));
}

#[test]
fn esc_on_a_question_keeps_what_follows_a_running_pass() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    app.after_pass = AfterPass::Run;
    app.tick();
    press(&mut app, KeyCode::Esc);
    assert!(matches!(app.after_pass, AfterPass::Run));
}

#[test]
fn another_panel_asks_its_questions_though_the_last_were_paused() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app(dir.path());
    app.tick();
    press(&mut app, KeyCode::Esc);
    app.forget_questions();
    app.tick();
    assert!(app.menu.is_some(), "asked on the new panel");
}
