//! Remembering an alias picked in the TUI; the GO pane and setting;
//! the keys.

use super::*;
use crate::annotate::cl_rules::Aliases;

#[test]
fn a_remembered_alias_is_appended_and_reads_back_as_the_top_layer() {
    let root = tempfile::tempdir().unwrap();
    let file = root.path().join("run").join("lupin").join("cl_aliases.tsv");
    remember_alias(&file, "EMP", "CL:0000049", "picked\tin the TUI").unwrap();
    remember_alias(&file, "BaEoMa", "CL:0000767", "second").unwrap();
    let text = std::fs::read_to_string(&file).unwrap();
    assert!(text.starts_with('#'), "a new file starts with its header");
    assert_eq!(text.matches("Columns:").count(), 1, "one header");
    assert!(
        text.contains("EMP\tCL:0000049\tpicked in the TUI\n"),
        "tabs in the note are flattened"
    );
    let mut a = Aliases::default();
    assert_eq!(a.add_tsv(&text, "project"), 2);
    assert_eq!(a.get("emp"), Some("CL:0000049"));
}

/// `lupin -o x` plus `more`, as the command line parses them.
fn cli(more: &[&str]) -> AnnotateCliArgs {
    #[derive(clap::Parser)]
    struct Cli {
        #[command(flatten)]
        annotate: AnnotateCliArgs,
    }
    let argv = ["lupin", "-o", "x"].iter().chain(more);
    <Cli as clap::Parser>::parse_from(argv).annotate
}

fn app_with_terms(terms: bool) -> App {
    use super::super::round::{ClusterView, RoundView};
    use crate::annotate::rounds::Term;
    let args = cli(&[]);
    let panel = [("GENE1".to_string(), "CT1".to_string())];
    let tree = crate::manifest::ontology::panel_tree_on(None, &panel);
    let mut app = App::new(args, PathBuf::new(), PathBuf::new(), tree);
    let term = |t: &str| Term {
        source: "go".into(),
        term: t.into(),
        effect: 0.5,
        q: Some(0.01),
        p: Some(0.001),
        nes: Some(1.5),
    };
    app.round = Some(RoundView {
        manifest: PathBuf::new(),
        clusters: vec![ClusterView {
            id: 0,
            cells: 10,
            label: Some("CT1".into()),
            candidates: Vec::new(),
            shares: Vec::new(),
            genes: vec![("GENE1".into(), 1.0)],
            terms: if terms {
                vec![term("TERM1"), term("TERM2")]
            } else {
                Vec::new()
            },
        }],
        cell_names: Vec::new(),
        cell_clusters: Vec::new(),
        expression: None,
        markers: std::collections::BTreeMap::new(),
        loose_cells: 0,
        decided: std::collections::BTreeSet::new(),
        rescorable: false,
    });
    app.focus = Focus::Genes;
    app
}

#[test]
fn tab_reaches_the_go_pane_only_when_the_round_scored_terms() {
    let tabs = |app: &mut App| -> Vec<Focus> {
        app.focus = Focus::Clusters;
        (0..4)
            .map(|_| {
                press(app, KeyCode::Tab);
                app.focus
            })
            .collect()
    };
    let mut with = app_with_terms(true);
    assert!(tabs(&mut with) == [Focus::Genes, Focus::Tree, Focus::Go, Focus::Clusters]);
    let mut without = app_with_terms(false);
    assert!(tabs(&mut without) == [Focus::Genes, Focus::Tree, Focus::Clusters, Focus::Genes]);
}

#[test]
fn the_go_pane_scrolls_its_terms_and_m_keeps_to_genes() {
    let mut app = app_with_terms(true);
    app.focus = Focus::Go;
    app.pane_key(KeyCode::Right);
    assert_eq!(app.go_shift, 0, "a one-word name does not slide");
    app.pane_key(KeyCode::End);
    assert_eq!(app.go_sel, 1);
    app.pane_key(KeyCode::Esc);
    assert!(app.focus == Focus::Clusters);

    app.focus = Focus::Genes;
    app.pane_key(KeyCode::Char('m'));
    assert!(app.gene_view == GeneView::Markers);
    app.pane_key(KeyCode::Char('m'));
    assert!(app.gene_view == GeneView::Specific);
}

#[test]
fn left_and_right_slide_a_long_go_name_and_a_new_term_starts_over() {
    let mut app = app_with_terms(true);
    if let Some(r) = app.round.as_mut() {
        r.clusters[0].terms[0].term = "w1 w2 w3".into();
    }
    app.focus = Focus::Go;
    for _ in 0..5 {
        app.pane_key(KeyCode::Right);
    }
    assert_eq!(app.go_shift, 2, "stops at the last word");
    app.pane_key(KeyCode::Left);
    assert_eq!(app.go_shift, 1);
    app.pane_key(KeyCode::Down);
    assert_eq!(app.go_shift, 0);
}

#[test]
fn the_go_setting_toggles_go_unless_gene_sets_are_named() {
    let mut a = cli(&[]);
    assert_eq!(Setting::Go.value(&a), "off");
    Setting::Go.adjust(&mut a, true);
    assert!(a.go);
    assert_eq!(Setting::Go.value(&a), "on");
    let mut b = cli(&["--gaf", "f"]);
    Setting::Go.adjust(&mut b, true);
    assert!(!b.go);
    assert_eq!(Setting::Go.value(&b), "--gaf");
}

fn press(app: &mut App, code: KeyCode) {
    app.key(KeyEvent::new(code, KeyModifiers::NONE));
}

fn ctrl_c(app: &mut App) {
    app.key(KeyEvent::new(KeyCode::Char('c'), KeyModifiers::CONTROL));
}

/// A stand-in for a running pass.
fn running(app: &mut App) {
    let child = std::process::Command::new("sleep")
        .arg("30")
        .spawn()
        .unwrap();
    app.child = Some(Running::new(child, Job::Pass));
}

#[test]
fn page_keys_move_by_a_page_or_to_an_end() {
    let mut sel = 3;
    step(&mut sel, 25, KeyCode::PageDown);
    assert_eq!(sel, 13);
    step(&mut sel, 25, KeyCode::PageDown);
    step(&mut sel, 25, KeyCode::PageDown);
    assert_eq!(sel, 24);
    step(&mut sel, 25, KeyCode::Home);
    assert_eq!(sel, 0);
    step(&mut sel, 25, KeyCode::Up);
    assert_eq!(sel, 0);
    step(&mut sel, 25, KeyCode::End);
    assert_eq!(sel, 24);
    step(&mut sel, 25, KeyCode::Char('z'));
    assert_eq!(sel, 24);
    step(&mut sel, 0, KeyCode::End);
    assert_eq!(sel, 0);
}

#[test]
fn ctrl_c_asks_before_dropping_unsaved_edits() {
    let mut app = app_with_terms(false);
    ctrl_c(&mut app);
    assert!(app.quit, "nothing unsaved: quits at once");

    let mut app = app_with_terms(false);
    app.edits.push(super::super::round::Edit::Keep {
        cluster: 0,
        reason: "why".into(),
    });
    ctrl_c(&mut app);
    assert!(!app.quit);
    assert!(app.status.contains("unsaved"), "{}", app.status);
    ctrl_c(&mut app);
    assert!(app.quit);
}

#[test]
fn x_stops_a_pass_only_when_asked_twice_from_any_pane() {
    for focus in [Focus::Genes, Focus::Clusters] {
        let mut app = app_with_terms(false);
        running(&mut app);
        app.focus = focus;
        press(&mut app, KeyCode::Char('x'));
        assert!(app.child.is_some());
        assert_eq!(app.status, "x again stops it");
        press(&mut app, KeyCode::Char('x'));
        assert!(app.child.is_none(), "{focus:?}");
    }
}

#[test]
fn h_hides_a_gene_and_x_does_not() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app_with_terms(false);
    app.data_search = crate::manifest::data_files::SearchPath::new(Some(dir.path()));
    press(&mut app, KeyCode::Char('x'));
    assert!(!app.status.starts_with("hiding"), "x hides nothing");
    press(&mut app, KeyCode::Char('h'));
    assert!(app.status.starts_with("hiding GENE1"), "{}", app.status);
}

#[test]
fn a_candidate_number_past_the_list_says_so() {
    let mut app = app_with_terms(false);
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Char('7'));
    assert!(app.edits.is_empty());
    assert!(app.status.contains("0 candidate"), "{}", app.status);
}

/// An app on the order view of `types`, its project folder in `dir`.
fn order_app(dir: &std::path::Path, types: &[&str]) -> App {
    let mut app = app_with_terms(false);
    // No round: the order view lists the panel's types.
    app.round = None;
    app.original = types
        .iter()
        .map(|t| ((*t).to_string(), Default::default()))
        .collect();
    app.data_search = crate::manifest::data_files::SearchPath::new(Some(dir));
    app.reload_order();
    // Labelled from a file, as a run with no round is.
    app.trajectory.labels = Some("labels.tsv".into());
    app.focus = Focus::Order;
    app
}

fn menu(app: &App) -> Option<&super::super::menu::Menu> {
    app.menu.as_ref()
}

#[test]
fn with_nothing_to_order_r_asks_how_and_a_root_goes_on_to_the_run() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = order_app(dir.path(), &["CT1", "CT2", "CT3"]);
    press(&mut app, KeyCode::Char('r'));
    let m = menu(&app).expect("the order menu");
    assert_eq!(m.choices.len(), 4);
    // `1` starts from one type: the types, then Enter on the second.
    press(&mut app, KeyCode::Char('1'));
    let m = menu(&app).expect("the root menu");
    assert_eq!(m.choices[1].label, "CT2");
    press(&mut app, KeyCode::Down);
    press(&mut app, KeyCode::Enter);
    assert!(app
        .trajectory
        .argv
        .windows(2)
        .any(|w| w == ["--root", "CT2"]));
    assert!(matches!(
        app.prompt.as_ref().map(|p| &p.pending),
        Some(Pending::TrajectoryOut { .. })
    ));
    // The root is kept: the next `r` asks only for the prefix.
    app.prompt = None;
    press(&mut app, KeyCode::Char('r'));
    assert!(menu(&app).is_none());
}

#[test]
fn the_order_menu_stays_away_when_something_orders_the_types() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = order_app(dir.path(), &["CT1", "CT2"]);
    app.trajectory.argv = vec!["--prior".into(), "p.tsv".into()];
    press(&mut app, KeyCode::Char('r'));
    assert!(menu(&app).is_none(), "a --prior file orders them");

    std::fs::create_dir_all(dir.path().join("lupin")).unwrap();
    std::fs::write(
        dir.path().join("lupin").join("precedence.tsv"),
        "CT1\tCT2\tprecedes\n",
    )
    .unwrap();
    let mut app = order_app(dir.path(), &["CT1", "CT2"]);
    press(&mut app, KeyCode::Char('r'));
    assert!(menu(&app).is_none(), "a statement orders them");
}

#[test]
fn state_the_order_myself_closes_the_menu_on_the_table_and_esc_cancels() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = order_app(dir.path(), &["CT1", "CT2"]);
    press(&mut app, KeyCode::Char('r'));
    press(&mut app, KeyCode::Down);
    press(&mut app, KeyCode::Enter);
    assert!(app.prompt.is_none());
    assert!(app.focus == Focus::Order);
    assert!(app.status.contains("space marks two types"));
    press(&mut app, KeyCode::Char('r'));
    press(&mut app, KeyCode::Esc);
    assert!(app.prompt.is_none());
    assert!(app.trajectory.argv.is_empty());
}

#[test]
fn in_the_order_view_tab_follows_the_columns_clusters_ordering_ontology() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = order_app(dir.path(), &["CT1", "CT2"]);
    app.focus = Focus::Clusters;
    let seen: Vec<Focus> = (0..4)
        .map(|_| {
            press(&mut app, KeyCode::Tab);
            app.focus
        })
        .collect();
    assert!(seen == [Focus::Genes, Focus::Order, Focus::Tree, Focus::Clusters]);
}

#[test]
fn the_ontology_beside_the_table_follows_its_selection() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = order_app(dir.path(), &["CT1", "CT2", "CT3"]);
    press(&mut app, KeyCode::Down);
    press(&mut app, KeyCode::Down);
    let TreeMode::Order(v) = &app.tree_mode else {
        panic!("the order view")
    };
    let tree = &v.tree;
    let at = tree.visible()[v.tree_sel];
    assert_eq!(tree.label(at), "CT3");
}

#[test]
fn with_no_labels_r_asks_where_they_come_from() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = order_app(dir.path(), &["CT1", "CT2"]);
    app.trajectory.labels = None;
    press(&mut app, KeyCode::Char('r'));
    let m = menu(&app).expect("the labels menu");
    assert_eq!(m.choices.len(), 2);
    press(&mut app, KeyCode::Char('2'));
    assert!(app.prompt.is_none());
    assert!(app.want_file == Some(FileWant::Labels));
}

/// top ─ group ─┬─ CT1
///              └─ CT2
/// top ─ elsewhere
const TYPES_OBO: &str = "format-version: 1.2

[Term]
id: CL:100
name: top

[Term]
id: CL:101
name: group
is_a: CL:100

[Term]
id: CL:104
name: CT1
is_a: CL:101

[Term]
id: CL:105
name: CT2
is_a: CL:101

[Term]
id: CL:106
name: elsewhere
is_a: CL:100
";

#[test]
fn annotates_ontology_view_opens_on_the_rounds_labels_and_toggles_with_d() {
    use super::super::ontology::Scope;
    let cl = crate::annotate::celltype_tree::ClTerms::parse(
        TYPES_OBO,
        &crate::annotate::cl_rules::shipped(),
    );
    let mut app = app_with_terms(false);
    let panel = [
        ("GENE1".to_string(), "CT1".to_string()),
        ("GENE2".to_string(), "CT2".to_string()),
    ];
    app.tree = crate::manifest::ontology::panel_tree_on(Some(&cl), &panel);
    app.cl = Some(cl);
    app.focus = Focus::Tree;
    press(&mut app, KeyCode::Char('o'));
    let TreeMode::Ontology(v) = &app.tree_mode else {
        panic!("the ontology view");
    };
    assert_eq!(v.scope, Scope::Data);
    assert_eq!(
        v.data.keys().map(String::as_str).collect::<Vec<_>>(),
        ["CL:104"],
        "the round labels its cluster CT1 only"
    );
    // Marking works on the data's tree as on the whole ontology.
    press(&mut app, KeyCode::Char(' '));
    assert_eq!(app.tree_marked, ["CT1"]);
    press(&mut app, KeyCode::Char('d'));
    let TreeMode::Ontology(v) = &app.tree_mode else {
        panic!("the ontology view");
    };
    assert_eq!(v.scope, Scope::All);
    assert_eq!(v.selected().unwrap().id, "CL:104", "on the same term");
    // A search still runs.
    press(&mut app, KeyCode::Char('/'));
    for c in "elsewhere".chars() {
        press(&mut app, KeyCode::Char(c));
    }
    press(&mut app, KeyCode::Enter);
    let TreeMode::Ontology(v) = &app.tree_mode else {
        panic!("the ontology view");
    };
    assert_eq!(v.selected().unwrap().id, "CL:106");
}

#[test]
fn in_the_order_view_r_runs_the_trajectory_from_any_pane() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = order_app(dir.path(), &["CT1", "CT2", "CT3"]);
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Char('r'));
    assert!(app.form.is_none(), "no cluster settings on an order view");
    assert!(
        menu(&app).is_some(),
        "the trajectory's own question instead"
    );
}

#[test]
fn a_failed_job_keeps_its_popup_until_a_key_which_does_nothing_else() {
    let mut app = app_with_terms(false);
    app.failed = Some("trajectory failed (exit status: 1); see the log".into());
    let focus = app.focus;
    press(&mut app, KeyCode::Tab);
    assert!(app.failed.is_none());
    assert!(app.focus == focus, "the key only closed the popup");
}

#[test]
fn b_hides_the_running_job_and_brings_it_back() {
    let mut app = app_with_terms(false);
    running(&mut app);
    press(&mut app, KeyCode::Char('b'));
    assert!(app.progress_hidden);
    press(&mut app, KeyCode::Char('b'));
    assert!(!app.progress_hidden);
}

/// An order view on `{dir}/run.senna.json` with no output chosen yet.
fn form_app(dir: &std::path::Path) -> App {
    let mut app = order_app(dir, &["CT1", "CT2"]);
    app.source = dir.join("run.senna.json");
    app.out_chosen = false;
    app.args.markers = Default::default();
    app
}

#[test]
fn annotating_from_the_labels_menu_opens_the_form_on_the_panel() {
    use super::super::app::{setting_row, Setting};
    use super::super::menu::Action;
    let dir = tempfile::tempdir().unwrap();
    let mut app = form_app(dir.path());
    app.choose(Action::Annotate);
    assert!(app.form.is_some());
    let summary = app.form_summary().0;
    assert!(summary.contains("run.L1"), "{summary}");
    assert!(summary.contains("then the trajectory"), "{summary}");
    assert_eq!(
        app.form.as_ref().unwrap().row,
        setting_row(Setting::Markers)
    );
    assert!(app.child.is_none(), "the form starts nothing on its own");
}

#[test]
fn a_in_the_clusters_opens_the_form_and_offers_the_trajectory_after() {
    let dir = tempfile::tempdir().unwrap();
    // The first round is taken: the form offers the next free one.
    std::fs::write(dir.path().join("run.L1.senna.json"), "{}").unwrap();
    let mut app = form_app(dir.path());
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Char('A'));
    assert!(app.form.is_some());
    let summary = app.form_summary().0;
    assert!(summary.contains("run.L2"), "{summary}");
    assert!(summary.contains("replaces nothing"), "{summary}");
    assert!(
        summary.contains("a choice to run the trajectory"),
        "{summary}"
    );
    press(&mut app, KeyCode::Esc);
    assert!(app.form.is_none(), "esc closes without running");
    assert!(app.child.is_none());
}

#[test]
fn the_form_refuses_to_run_without_a_panel_and_enter_on_it_asks_for_one() {
    use super::super::app::{setting_row, FileWant, Setting};
    let dir = tempfile::tempdir().unwrap();
    let mut app = form_app(dir.path());
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Char('A'));
    app.form.as_mut().unwrap().row = setting_row(Setting::Run);
    press(&mut app, KeyCode::Enter);
    assert!(app.form.is_some() && app.child.is_none());
    let note = app.form.as_ref().unwrap().note.clone().unwrap_or_default();
    assert!(note.contains("marker panel"), "{note}");
    assert_eq!(
        app.form.as_ref().unwrap().row,
        setting_row(Setting::Markers),
        "back on the panel row"
    );
    // Enter there asks for a file and starts nothing.
    press(&mut app, KeyCode::Enter);
    assert_eq!(app.want_file, Some(FileWant::Markers));
    assert!(app.form.is_some() && app.child.is_none());
    // Enter on a value row only says how to run.
    app.form.as_mut().unwrap().row = setting_row(Setting::Knn);
    press(&mut app, KeyCode::Enter);
    assert!(app.form.is_some() && app.child.is_none());
}

#[test]
fn a_pass_that_replaces_a_round_asks_in_the_form_before_it_runs() {
    use super::super::app::{setting_row, Setting};
    let dir = tempfile::tempdir().unwrap();
    let mut app = form_app(dir.path());
    app.args.markers = "panel.tsv".into();
    app.out_chosen = true;
    app.args.out = dir.path().join("run.L1").to_string_lossy().into();
    app.target = dir.path().join("run.L1.senna.json");
    std::fs::write(&app.target, "{}").unwrap();
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Char('A'));
    assert_eq!(
        app.form.as_ref().unwrap().row,
        setting_row(Setting::Run),
        "a panel: on the run row"
    );
    assert!(app.form_summary().0.contains("replaces run.L1.senna.json"));
    press(&mut app, KeyCode::Enter);
    assert!(
        app.form.as_ref().unwrap().confirm && app.child.is_none(),
        "asked, not started"
    );
    let note = app.form.as_ref().unwrap().note.clone().unwrap_or_default();
    assert!(note.contains("run again to confirm"), "{note}");
}

#[test]
fn r_in_annotate_opens_the_same_form_without_a_trajectory() {
    let mut app = app_with_terms(false);
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Char('r'));
    assert!(app.form.is_some());
    assert!(!app.form_summary().0.contains("trajectory"));
    press(&mut app, KeyCode::Char('r'));
    assert!(app.form.is_none(), "r closes it again");
}

#[test]
fn after_a_pass_from_the_order_view_the_trajectory_is_offered() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = form_app(dir.path());
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Char('A'));
    press(&mut app, KeyCode::Esc);
    // As a run from the form leaves things: flags set, the form closed.
    press(&mut app, KeyCode::Char('A'));
    app.form = None;
    let mut child = std::process::Command::new("true").spawn().unwrap();
    child.wait().unwrap();
    app.child = Some(Running::new(child, Job::Pass));
    app.tick();
    let m = menu(&app).expect("the offer");
    assert_eq!(m.choices[0].label, "Run the trajectory on the new labels");
    press(&mut app, KeyCode::Char('2'));
    assert!(app.prompt.is_none() && app.status.contains("r runs"));
}

/// An order view with figures on screen and the ordering column focused.
fn figures_app(dir: &std::path::Path) -> App {
    let mut app = order_app(dir, &["CT1", "CT2"]);
    app.figures = Some(super::super::figure_pane::tests::pane());
    app.focus = Focus::Order;
    app
}

#[test]
fn in_the_order_view_t_restyles_and_capital_t_leaves() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = figures_app(dir.path());
    press(&mut app, KeyCode::Char('t'));
    assert!(app.in_order(), "t restyles the labels");
    assert!(app.status.starts_with("labels"), "{}", app.status);
    // With the table shown, t only says how.
    press(&mut app, KeyCode::Char('V'));
    press(&mut app, KeyCode::Char('t'));
    assert!(app.in_order());
    assert!(app.status.contains("T back to the tree"), "{}", app.status);
    press(&mut app, KeyCode::Char('T'));
    assert!(!app.in_order(), "T leaves the order view");
}

#[test]
fn the_figures_keys_wait_while_the_ontology_column_has_the_focus() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = figures_app(dir.path());
    app.focus = Focus::Tree;
    let style = format!("{:?}", app.figures.as_ref().unwrap().style);
    press(&mut app, KeyCode::Char('t'));
    press(&mut app, KeyCode::Char('c'));
    press(&mut app, KeyCode::Char('w'));
    let v = app.figures.as_ref().unwrap();
    assert_eq!(format!("{:?}", v.style), style, "no restyle");
    assert!(v.grid.is_none(), "no grid");
}

#[test]
fn esc_closes_the_exports_strip_then_the_grid_before_leaving() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = figures_app(dir.path());
    press(&mut app, KeyCode::Char('f'));
    assert!(app.figures.as_ref().unwrap().exports.open);
    press(&mut app, KeyCode::Esc);
    assert!(!app.figures.as_ref().unwrap().exports.open);
    assert_eq!(app.focus, Focus::Order, "only the strip closed");
    press(&mut app, KeyCode::Char('w'));
    press(&mut app, KeyCode::Esc);
    assert!(app.figures.as_ref().unwrap().grid.is_none());
    assert_eq!(app.focus, Focus::Order, "only the grid closed");
    press(&mut app, KeyCode::Esc);
    assert_eq!(app.focus, Focus::Clusters);
}

#[test]
fn a_grid_tile_opens_with_its_number_or_space() {
    use crate::trajectory::figures::Panel;
    let dir = tempfile::tempdir().unwrap();
    let mut app = figures_app(dir.path());
    press(&mut app, KeyCode::Char('w'));
    press(&mut app, KeyCode::Char('4'));
    let v = app.figures.as_ref().unwrap();
    assert!(v.grid.is_none());
    assert_eq!(v.current(), Panel::Order);
    press(&mut app, KeyCode::Char('w'));
    press(&mut app, KeyCode::Char('9'));
    assert!(app.status.contains("tile(s)"), "{}", app.status);
    press(&mut app, KeyCode::Char(' '));
    assert!(app.figures.as_ref().unwrap().grid.is_none());
}

#[test]
fn f_waits_for_the_figures() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = figures_app(dir.path());
    press(&mut app, KeyCode::Char('V'));
    press(&mut app, KeyCode::Char('f'));
    assert!(!app.figures.as_ref().unwrap().exports.open);
    assert!(
        app.status.starts_with("v shows the figures"),
        "{}",
        app.status
    );
}

#[test]
fn a_failed_job_notice_takes_the_key_before_an_open_prompt() {
    let mut app = app_with_terms(false);
    app.prompt = Some(Prompt {
        title: "why?".into(),
        text: String::new(),
        pending: Pending::Hide,
    });
    app.failed = Some("trajectory failed".into());
    press(&mut app, KeyCode::Char('a'));
    assert!(app.failed.is_none());
    assert_eq!(
        app.prompt.as_ref().unwrap().text,
        "",
        "the prompt got nothing"
    );
}

/// The family member `name` of a family made by `family_dir`, to show.
fn pick(dir: &std::path::Path, name: &str) -> crate::manifest::family::Pick {
    crate::manifest::family::family(&dir.join("X.senna.json"))
        .into_iter()
        .find(|m| m.name() == name)
        .expect("a member")
        .pick
}

/// An app on the run `X` of a family made by the family tests.
fn family_app(dir: &std::path::Path) -> App {
    crate::manifest::family::tests::family_dir(dir);
    let mut app = app_with_terms(false);
    app.round = None;
    app.source = dir.join("X.senna.json");
    app
}

#[test]
fn g_lists_the_family_and_a_round_opens_with_its_clusters() {
    let tmp = tempfile::tempdir().unwrap();
    let mut app = family_app(tmp.path());
    // No output chosen: switching offers the next free round.
    app.out_chosen = false;
    press(&mut app, KeyCode::Char('g'));
    let m = menu(&app).expect("the runs");
    assert_eq!(m.choices.len(), 5);
    assert_eq!(m.sel, 0, "on the run on screen");
    assert!(
        m.choices[0].label.starts_with("● X ·"),
        "{}",
        m.choices[0].label
    );
    // `2`: the round X.L0.
    press(&mut app, KeyCode::Char('2'));
    assert!(app.menu.is_none());
    let r = app.round.as_ref().expect("the round");
    assert_eq!(r.clusters.len(), 2);
    assert!(app.shown_manifest().ends_with("X.L0.senna.json"));
    // Passes start from the round shown, offered the next free round under
    // it, as when the TUI is opened on it.
    assert!(app.source.ends_with("X.L0.senna.json"));
    let (from, out) = crate::manifest::family::pass_origin(&tmp.path().join("X.L0.senna.json"));
    assert_eq!(app.source, from);
    assert_eq!(app.args.out.as_ref(), out);
    assert!(out.ends_with("X.L0.L2"), "L1 exists: {out}");
    // And back to the run: no round.
    press(&mut app, KeyCode::Char('g'));
    assert_eq!(menu(&app).unwrap().sel, 1, "the round is on screen now");
    press(&mut app, KeyCode::Char('1'));
    assert!(app.round.is_none());
}

#[test]
fn a_round_from_the_list_gives_the_order_view_its_labels() {
    let tmp = tempfile::tempdir().unwrap();
    let mut app = order_app(tmp.path(), &["CT9"]);
    let mut fam = family_app(tmp.path());
    std::mem::swap(&mut app.source, &mut fam.source);
    app.open_run(pick(tmp.path(), "X.L0"), false);
    let TreeMode::Order(v) = &app.tree_mode else {
        panic!("still the order view")
    };
    let types: Vec<&str> = v.types.iter().map(|(t, _)| t.as_str()).collect();
    assert_eq!(types, ["CT1", "CT2"]);
}

#[test]
fn unsaved_edits_are_asked_about_and_a_running_job_refuses() {
    let tmp = tempfile::tempdir().unwrap();
    let mut app = family_app(tmp.path());
    let round = pick(tmp.path(), "X.L0");
    app.open_run(round.clone(), false);
    app.edits.push(Edit::Keep {
        cluster: 0,
        reason: "kept".into(),
    });
    app.open_run(pick(tmp.path(), "X"), false);
    let m = menu(&app).expect("asked first");
    assert!(m.question.contains("1 unsaved edit"), "{}", m.question);
    assert!(app.round.is_some(), "nothing switched yet");
    // Esc keeps the edits; the one choice drops them and opens the run.
    assert_eq!(m.choices.len(), 1);
    press(&mut app, KeyCode::Esc);
    assert_eq!(app.edits.len(), 1);
    assert!(app.status.contains("kept"), "{}", app.status);
    app.open_run(pick(tmp.path(), "X"), true);
    assert!(app.round.is_none() && app.edits.is_empty());
    // A running job: no switching, no list.
    running(&mut app);
    app.open_run(round, false);
    assert!(app.round.is_none());
    assert!(
        app.status.contains("wait for the running job"),
        "{}",
        app.status
    );
    press(&mut app, KeyCode::Char('g'));
    assert!(app.menu.is_none());
    app.stop();
}

#[test]
fn switching_runs_gives_the_order_view_the_edges_of_what_is_shown() {
    let tmp = tempfile::tempdir().unwrap();
    let mut app = order_app(tmp.path(), &["CT9"]);
    let mut fam = family_app(tmp.path());
    std::mem::swap(&mut app.source, &mut fam.source);
    // The edges of an earlier trajectory, on screen.
    if let TreeMode::Order(v) = &mut app.tree_mode {
        v.edges = vec![crate::trajectory::edges::EdgeRow {
            a: "CT8".into(),
            b: "CT9".into(),
            connectivity: 0.5,
            in_prior: true,
            verdict: None,
            order_agreement: f32::NAN,
        }];
    }
    app.open_run(pick(tmp.path(), "X.L0"), false);
    let TreeMode::Order(v) = &app.tree_mode else {
        panic!("still the order view")
    };
    assert!(
        !v.edges.iter().any(|e| e.is_pair("CT8", "CT9")),
        "the earlier trajectory's edges are gone"
    );
}
