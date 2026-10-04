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
    app.child = Some((child, Instant::now(), Job::Pass));
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
fn x_stops_a_pass_only_when_asked_twice_and_never_from_the_genes_pane() {
    let mut app = app_with_terms(false);
    running(&mut app);
    app.focus = Focus::Genes;
    press(&mut app, KeyCode::Char('x'));
    press(&mut app, KeyCode::Char('x'));
    assert!(app.child.is_some(), "x hides genes there");

    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Char('x'));
    assert!(app.child.is_some());
    assert_eq!(app.status, "x again stops it");
    press(&mut app, KeyCode::Char('x'));
    assert!(app.child.is_none());
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
    match app.prompt.as_ref().map(|p| &p.pending) {
        Some(Pending::Choose(m)) => Some(m),
        _ => None,
    }
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
fn in_the_order_view_tab_visits_the_ontology_then_the_precedence_table() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = order_app(dir.path(), &["CT1", "CT2"]);
    app.focus = Focus::Clusters;
    let seen: Vec<Focus> = (0..4)
        .map(|_| {
            press(&mut app, KeyCode::Tab);
            app.focus
        })
        .collect();
    assert!(seen == [Focus::Genes, Focus::Tree, Focus::Order, Focus::Clusters]);
}

#[test]
fn the_ontology_beside_the_table_follows_its_selection() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = order_app(dir.path(), &["CT1", "CT2", "CT3"]);
    press(&mut app, KeyCode::Down);
    press(&mut app, KeyCode::Down);
    let tree = app.order_tree.as_ref().unwrap();
    let at = tree.visible()[app.order_tree_sel];
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
    assert!(!app.settings_open, "no cluster settings on an order view");
    assert!(
        menu(&app).is_some(),
        "the trajectory's own question instead"
    );
}
