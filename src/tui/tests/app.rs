//! Remembering an alias picked in the TUI; the GO pane and setting;
//! the keys.

use super::*;
use crate::annotate::cl_rules::Aliases;

#[test]
fn a_remembered_alias_is_appended_and_reads_back_as_the_top_layer() {
    let root = tempfile::tempdir().unwrap();
    let file = root.path().join("run").join("lupin").join("cl_aliases.tsv");
    remember_alias(&file, "CT1", "CL:9000001", "picked\tin the TUI").unwrap();
    remember_alias(&file, "CT2", "CL:9000002", "second").unwrap();
    let text = std::fs::read_to_string(&file).unwrap();
    assert!(text.starts_with('#'), "a new file starts with its header");
    assert_eq!(text.matches("Columns:").count(), 1, "one header");
    assert!(
        text.contains("CT1\tCL:9000001\tpicked in the TUI\n"),
        "tabs in the note are flattened"
    );
    let mut a = Aliases::default();
    assert_eq!(a.add_tsv(&text, "project"), 2);
    assert_eq!(a.get("ct1"), Some("CL:9000001"));
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
        by: Default::default(),
        figures: Default::default(),
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
    assert!(tabs(&mut with) == [Focus::Genes, Focus::Tree, Focus::Candidates, Focus::Go]);
    let mut without = app_with_terms(false);
    assert!(
        tabs(&mut without)
            == [
                Focus::Genes,
                Focus::Tree,
                Focus::Candidates,
                Focus::Clusters
            ]
    );
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
fn a_running_save_offers_no_stop_and_x_says_so_at_once() {
    let mut app = app_with_terms(false);
    let child = std::process::Command::new("sleep")
        .arg("30")
        .spawn()
        .unwrap();
    app.child = Some(Running::new(child, Job::Save(1)));
    app.tick();
    assert!(!app.status.contains("x: stop"), "{}", app.status);
    assert!(!drawn(&app).contains("x stop"), "the popup offers no stop");
    press(&mut app, KeyCode::Char('x'));
    assert!(app.status.contains("finishes on its own"), "{}", app.status);
    assert!(app.child.is_some());
    if let Some(mut r) = app.child.take() {
        let _ = r.child.kill();
        let _ = r.child.wait();
    }
}

#[test]
fn the_status_follows_the_latest_key_and_a_dropped_question_goes_with_it() {
    let mut app = app_with_terms(false);
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Backspace);
    assert!(app.status.contains("to undo"), "{}", app.status);
    press(&mut app, KeyCode::Down);
    assert!(!app.status.contains("to undo"), "{}", app.status);

    app.edits.push(super::super::round::Edit::Keep {
        cluster: 0,
        reason: "why".into(),
    });
    press(&mut app, KeyCode::Char('q'));
    assert!(app.status.contains("q again"), "{}", app.status);
    press(&mut app, KeyCode::Down);
    assert!(
        !app.status.contains("q again"),
        "the question is dropped: {}",
        app.status
    );

    running(&mut app);
    press(&mut app, KeyCode::Char('x'));
    assert!(app.status.contains("x again"), "{}", app.status);
    press(&mut app, KeyCode::Down);
    assert!(!app.status.contains("x again"), "{}", app.status);
    app.stop();
}

#[test]
fn a_last_stage_running_long_gives_no_eta() {
    let mut app = app_with_terms(false);
    running(&mut app);
    if let Some(r) = app.child.as_mut() {
        r.started = Instant::now() - Duration::from_secs(660);
        r.progress = JobProgress {
            report: Some(crate::progress::Report {
                done: 19.0,
                end: 20.0,
                total: 20.0,
                stage: "GO terms".into(),
            }),
            reported_at: Duration::from_secs(20),
        };
    }
    let screen = drawn(&app);
    assert!(screen.contains("elapsed 11m00s"), "{screen}");
    assert!(
        !screen.contains("ETA 0s"),
        "no time left to count on: {screen}"
    );
    assert!(screen.contains("no ETA: this stage runs past"), "{screen}");
    app.stop();
}

#[test]
fn a_running_job_redraws_on_the_clock_even_while_a_key_is_asked_again() {
    let mut app = app_with_terms(false);
    assert!(!app.redraw_due(Duration::from_secs(5)), "nothing runs");
    running(&mut app);
    press(&mut app, KeyCode::Char('x'));
    app.tick();
    assert_eq!(
        app.status, "x again stops it",
        "the question holds the status"
    );
    assert!(!app.redraw_due(Duration::from_millis(200)));
    assert!(
        app.redraw_due(Duration::from_secs(1)),
        "the popup's clock moves"
    );
    app.stop();
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
fn with_no_round_the_cluster_keys_say_there_are_no_clusters() {
    let mut app = app_with_terms(false);
    app.round = None;
    app.focus = Focus::Clusters;
    for key in ['1', 'S', 'a'] {
        press(&mut app, KeyCode::Char(key));
        assert_eq!(app.status, "no clusters yet", "{key}");
    }
    assert!(!app.cand_all, "nothing toggled for the next round");
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
    app.trajectory_tui = true;
    app.focus = Focus::Order;
    app
}

fn menu(app: &App) -> Option<&super::super::menu::Menu> {
    app.menu.as_ref()
}

#[test]
fn a_menus_long_options_are_shown_whole_when_the_screen_has_room() {
    use super::super::menu::{numbered, Action, Menu};
    let long = |n: usize| format!("{}end)", "x".repeat(n));
    let choices = |label: String| numbered(vec![(label, String::new(), Action::ListRuns)]);
    // The runs list takes the screen's width.
    let mut app = app_with_terms(false);
    let mut m = Menu::new("runs", choices(long(120)), "");
    m.wide = true;
    app.menu = Some(m);
    assert!(drawn(&app).contains("xend)"), "{}", drawn(&app));
    // A menu over the order view's clusters goes past the pane to fit.
    let dir = tempfile::tempdir().unwrap();
    let mut app = order_app(dir.path(), &["CT1", "CT2"]);
    app.menu = Some(Menu::new("how?", choices(long(70)), ""));
    assert!(drawn(&app).contains("xend)"), "{}", drawn(&app));
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
fn a_labels_file_gives_the_order_table_its_types_and_r_asks_how_to_order_them() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = order_app(dir.path(), &[]);
    let labels = dir.path().join("labels.tsv");
    std::fs::write(&labels, "c1\tCT1\nc2\tCT1\nc3\tCT2\n").unwrap();
    app.trajectory.labels = Some(labels.to_string_lossy().into());
    app.reload_order();
    let TreeMode::Order(v) = &app.tree_mode else {
        panic!("the order view")
    };
    assert_eq!(v.types, [("CT1".to_string(), 2), ("CT2".to_string(), 1)]);
    press(&mut app, KeyCode::Char('r'));
    assert!(
        menu(&app).is_some(),
        "nothing orders them: the menu asks how"
    );
    assert!(app.prompt.is_none());
}

#[test]
fn lupin_trajectorys_labels_reach_the_tui_once() {
    #[derive(clap::Parser)]
    struct Cli {
        #[command(flatten)]
        t: crate::trajectory::run::TrajectoryArgs,
    }
    let argv = ["lupin", "-f", "run.senna.json", "--labels", "l.tsv"];
    let t = <Cli as clap::Parser>::parse_from(argv).t;
    let run = super::super::order::TrajectoryRun::of(&t);
    assert_eq!(run.labels.as_deref(), Some("l.tsv"));
    assert!(!run.argv.iter().any(|a| a == "--labels"), "{:?}", run.argv);
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
fn in_the_order_view_tab_follows_the_columns_clusters_candidates_ordering_ontology() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = order_app(dir.path(), &["CT1", "CT2"]);
    app.focus = Focus::Clusters;
    let seen: Vec<Focus> = (0..5)
        .map(|_| {
            press(&mut app, KeyCode::Tab);
            app.focus
        })
        .collect();
    assert!(
        seen == [
            Focus::Genes,
            Focus::Candidates,
            Focus::Order,
            Focus::Tree,
            Focus::Clusters
        ]
    );
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
fn the_ontology_opens_on_a_term_the_datas_tree_does_not_list() {
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
    // `elsewhere` is in no type's lineage: the whole ontology, on it.
    let v = app.ontology_on(Some("CL:106".into())).unwrap();
    assert_eq!(v.selected().unwrap().id, "CL:106");
    // A term the data's tree lists keeps that tree.
    let v = app.ontology_on(Some("CL:104".into())).unwrap();
    assert_eq!(v.scope, super::super::ontology::Scope::Data);
    assert_eq!(v.selected().unwrap().id, "CL:104");
}

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
    (app.out_chosen, app.out_typed) = (false, false);
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
        setting_row(&app.args, Setting::Markers)
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
    app.form.as_mut().unwrap().row = setting_row(&app.args, Setting::Run);
    press(&mut app, KeyCode::Enter);
    assert!(app.form.is_some() && app.child.is_none());
    let note = app.form.as_ref().unwrap().note.clone().unwrap_or_default();
    assert!(note.contains("marker panel"), "{note}");
    assert_eq!(
        app.form.as_ref().unwrap().row,
        setting_row(&app.args, Setting::Markers),
        "back on the panel row"
    );
    // Enter there picks the panel in the file browser, and nothing starts.
    press(&mut app, KeyCode::Enter);
    assert!(app.menu.is_none());
    assert_eq!(app.want_file, Some(FileWant::Markers));
    assert!(app.form.is_some() && app.child.is_none());
    // Enter on a value row only says how to run.
    app.form.as_mut().unwrap().row = setting_row(&app.args, Setting::Knn);
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
        setting_row(&app.args, Setting::Run),
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
fn picking_a_new_panel_asks_again_before_a_pass_that_replaces_a_round() {
    use super::super::app::{setting_row, Setting};
    let dir = tempfile::tempdir().unwrap();
    let mut app = form_app(dir.path());
    let panel = dir.path().join("panel.tsv");
    std::fs::write(&panel, "G1\tCT1\nG2\tCT2\n").unwrap();
    app.args.markers = panel.to_string_lossy().into();
    app.out_chosen = true;
    app.args.out = dir.path().join("run.L1").to_string_lossy().into();
    app.target = dir.path().join("run.L1.senna.json");
    std::fs::write(&app.target, "{}").unwrap();
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Char('A'));
    press(&mut app, KeyCode::Enter);
    assert!(app.form.as_ref().unwrap().confirm, "asked once");
    app.set_markers(Some(&panel));
    assert_eq!(
        app.form.as_ref().unwrap().row,
        setting_row(&app.args, Setting::Run)
    );
    assert!(
        !app.form.as_ref().unwrap().confirm,
        "a new panel asks again, as a new output does"
    );
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
    let v = app.figures.as_ref().unwrap();
    assert_eq!(app.status, v.title(v.current()), "the grid's keys are gone");
    press(&mut app, KeyCode::Esc);
    assert_eq!(app.focus, Focus::Clusters);
}

#[test]
fn f_in_the_grid_says_the_strip_is_a_figures_and_opens_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = figures_app(dir.path());
    press(&mut app, KeyCode::Char('w'));
    press(&mut app, KeyCode::Char('f'));
    let v = app.figures.as_ref().unwrap();
    assert!(!v.exports.open, "no strip the grid would hide");
    assert!(v.grid.is_some(), "the grid stays");
    assert!(app.status.contains("opens the tile"), "{}", app.status);
    assert!(!super::super::figure_pane::GRID_KEYS.contains("f exports"));
}

#[test]
fn a_short_pane_gives_the_exports_strip_room_for_its_selected_row() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = map_app(dir.path(), true);
    app.focus = Focus::Tree;
    let x = &mut app.cell_map.as_mut().unwrap().exports;
    x.not_listed = vec!["run.cells.extra.pdf".into()];
    x.sel = 0;
    press(&mut app, KeyCode::Char('f'));
    let screen = drawn_in(&app, 100, 30).concat();
    assert!(screen.contains(" exports ("), "the strip");
    assert!(
        screen.contains("run.cells.extra.pdf"),
        "and its selected row"
    );
    assert!(app.cell_map.as_ref().unwrap().exports.fits.get());
    use super::super::ui::exports_height;
    assert_eq!(exports_height(30), Some(8), "room: the strip's own height");
    assert_eq!(exports_height(12), Some(6), "what the figure leaves");
    assert_eq!(exports_height(8), Some(3), "a row, taken from the figure");
    assert_eq!(exports_height(6), Some(3), "down to the figure's own row");
    assert_eq!(exports_height(5), None, "else not shown");
}

#[test]
fn a_grid_tile_opens_with_its_number_or_space() {
    use crate::figures::Panel;
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
fn the_status_follows_the_figures_and_t_c_say_only_the_scatter_takes_them() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = figures_app(dir.path());
    press(&mut app, KeyCode::Char('w'));
    assert!(app.status.contains("arrows choose"), "{}", app.status);
    press(&mut app, KeyCode::Char('V'));
    assert!(!app.status.contains("arrows choose"), "{}", app.status);
    press(&mut app, KeyCode::Char('v'));
    let p = app.figures.as_mut().unwrap();
    p.open_grid();
    p.grid.as_mut().unwrap().sel = 3;
    p.open_tile();
    assert!(!p.current().is_scatter());
    let style = p.style;
    for key in ['t', 'c'] {
        press(&mut app, KeyCode::Char(key));
        assert!(app.status.contains("scatter"), "{key}: {}", app.status);
    }
    assert!(
        app.figures.as_ref().unwrap().style == style,
        "nothing changed"
    );
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
    (app.out_chosen, app.out_typed) = (false, false);
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
fn a_pass_run_on_one_run_leaves_the_next_run_its_own_output() {
    let tmp = tempfile::tempdir().unwrap();
    let mut app = family_app(tmp.path());
    // As a pass from the form leaves things, with no output typed.
    (app.out_chosen, app.out_typed) = (true, false);
    app.args.out = "X.L1".into();
    app.open_run(pick(tmp.path(), "X.L0"), false);
    let (_, out) = crate::manifest::family::pass_origin(&tmp.path().join("X.L0.senna.json"));
    assert_eq!(app.args.out.as_ref(), out, "derived for the run shown");
    // A typed output stays.
    app.out_typed = true;
    app.args.out = "mine".into();
    app.open_run(pick(tmp.path(), "X"), false);
    assert_eq!(app.args.out.as_ref(), "mine");
}

#[test]
fn a_round_made_on_a_trajectory_run_lists_as_a_round() {
    let tmp = tempfile::tempdir().unwrap();
    let mut app = family_app(tmp.path());
    // A pass on X.T1 copies its trajectory section into the new round.
    std::fs::write(
        tmp.path().join("X.T1.L1.senna.json"),
        r#"{"version":2,"kind":"topic","prefix":"X.T1.L1",
            "cluster":{"clusters":"X.clusters.parquet"},
            "annotate":{"argmax":"X.L0.argmax.tsv","source":"X.T1.senna.json"},
            "trajectory":{"prior":"t.tsv"}}"#,
    )
    .unwrap();
    press(&mut app, KeyCode::Char('g'));
    let m = menu(&app).expect("the runs");
    let row = |name: &str| {
        m.choices
            .iter()
            .find(|c| c.label.contains(&format!("{name} ·")))
            .map(|c| c.label.clone())
            .unwrap()
    };
    assert!(row("X.T1.L1").contains("· round ·"), "{}", row("X.T1.L1"));
    assert!(row("X.T1").contains("· trajectory ·"), "{}", row("X.T1"));
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
fn a_save_opens_the_round_after_the_one_saved_whatever_the_pass_output() {
    let tmp = tempfile::tempdir().unwrap();
    let mut app = family_app(tmp.path());
    // The run opened is the round itself: its passes would go to X.L0.L2.
    app.source = tmp.path().join("X.L0.senna.json");
    app.target = tmp.path().join("X.L0.L2.senna.json");
    app.open_run(pick(tmp.path(), "X.L0"), false);
    // As `lupin relabel --next` leaves it: the next round beside it.
    let text = std::fs::read_to_string(tmp.path().join("X.L0.senna.json")).unwrap();
    std::fs::write(tmp.path().join("X.L0.r1.senna.json"), text).unwrap();
    let mut child = std::process::Command::new("true").spawn().unwrap();
    child.wait().unwrap();
    app.child = Some(Running::new(child, Job::Save(0)));
    app.tick();
    assert!(
        app.shown_manifest().ends_with("X.L0.r1.senna.json"),
        "{}: {}",
        app.shown_manifest().display(),
        app.status
    );
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

#[test]
fn opening_a_run_with_no_trajectory_drops_the_earlier_ones_figures() {
    let tmp = tempfile::tempdir().unwrap();
    let mut app = figures_app(tmp.path());
    app.trajectory_tui = true;
    let fam = family_app(tmp.path());
    app.source = fam.source.clone();
    app.open_run(pick(tmp.path(), "X.L0"), false);
    assert!(app.figures.is_none(), "X.L0 has no trajectory");
}

/// The form's marker panel, picked in the browser popup over it.
fn browsing_for_markers(dir: &Path) -> App {
    let mut app = app_with_terms(false);
    app.args.markers = Default::default();
    app.form = Some(Form::default());
    let want = super::super::picker::Want::Markers(None, 0);
    let b = super::super::picker::Picker::new(FileWant::Markers.title(), dir, want);
    app.browser = Some((FileWant::Markers, b));
    app
}

#[test]
fn the_browser_popup_picks_the_forms_marker_panel_and_draws_over_it() {
    let root = tempfile::tempdir().unwrap();
    std::fs::create_dir(root.path().join("sub")).unwrap();
    let panel = root.path().join("panel.tsv");
    std::fs::write(&panel, "GENE1\tCT1\nGENE2\tCT2\n").unwrap();
    let mut app = browsing_for_markers(root.path());

    let backend = ratatui::backend::TestBackend::new(120, 40);
    let mut term = ratatui::Terminal::new(backend).unwrap();
    term.draw(|f| super::super::ui::draw(f, &app)).unwrap();
    let screen: String = term
        .backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect();
    assert!(screen.contains("Pick a marker panel"), "the popup's title");
    assert!(screen.contains("panel.tsv"), "the panel is listed");
    assert!(screen.contains("sub/"), "directories are listed");

    // The browser starts on the panel; enter takes it, and the form stays.
    press(&mut app, KeyCode::Enter);
    assert!(app.browser.is_none());
    assert!(app.form.is_some(), "the form stays open to run");
    assert_eq!(
        Path::new(app.args.markers.as_ref()).file_name(),
        panel.file_name()
    );
}

#[test]
fn the_browser_shows_the_end_of_a_long_folder_and_esc_keeps_the_panel() {
    let root = tempfile::tempdir().unwrap();
    let deep = root
        .path()
        .join("a".repeat(60))
        .join("b".repeat(60))
        .join("annotate_here");
    std::fs::create_dir_all(&deep).unwrap();
    let mut app = browsing_for_markers(&deep);
    let backend = ratatui::backend::TestBackend::new(120, 40);
    let mut term = ratatui::Terminal::new(backend).unwrap();
    term.draw(|f| super::super::ui::draw(f, &app)).unwrap();
    let screen: String = term
        .backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect();
    assert!(
        screen.contains("annotate_here"),
        "the folder's own name: {screen}"
    );
    // Cancelled, an earlier panel stays and the status says so.
    app.args.markers = "earlier.tsv".into();
    press(&mut app, KeyCode::Esc);
    assert!(app.status.contains("earlier.tsv"), "{}", app.status);
    assert!(!app.status.contains("no marker panel"), "{}", app.status);
}

#[test]
fn esc_closes_the_browser_popup_and_leaves_the_form() {
    let root = tempfile::tempdir().unwrap();
    let mut app = browsing_for_markers(root.path());
    press(&mut app, KeyCode::Esc);
    assert!(app.browser.is_none());
    assert!(app.form.is_some(), "esc closes only the browser");
    assert!(app.args.markers.is_empty());
    assert_eq!(app.status, "no marker panel picked");
}

/// The screen `app` draws, as one string.
fn drawn(app: &App) -> String {
    drawn_in(app, 160, 45).concat()
}

/// The screen `app` draws on a `width` × `height` terminal, by row.
fn drawn_in(app: &App, width: u16, height: u16) -> Vec<String> {
    let backend = ratatui::backend::TestBackend::new(width, height);
    let mut term = ratatui::Terminal::new(backend).unwrap();
    term.draw(|f| super::super::ui::draw(f, app)).unwrap();
    let buffer = term.backend().buffer();
    buffer
        .content()
        .chunks(usize::from(width))
        .map(|row| row.iter().map(|c| c.symbol()).collect())
        .collect()
}

#[test]
fn a_narrow_cluster_list_leaves_out_the_cell_counts_rather_than_cut_them() {
    let dir = tempfile::tempdir().unwrap();
    let annotate = {
        let mut app = app_with_terms(false);
        app.focus = Focus::Clusters;
        app
    };
    let order = {
        let mut app = order_app(dir.path(), &["CT1", "CT2"]);
        app.round = app_with_terms(false).round.take();
        app.focus = Focus::Clusters;
        app
    };
    for (what, mut app) in [("annotate", annotate), ("order", order)] {
        app.round.as_mut().unwrap().clusters[0].cells = 1234567;
        let row = |screen: &[String]| {
            screen
                .iter()
                .find(|r| r.contains("K0 "))
                .cloned()
                .unwrap_or_default()
        };
        let wide = drawn_in(&app, 200, 50);
        assert!(row(&wide).contains("1234567"), "{what}: {}", row(&wide));
        let narrow = drawn_in(&app, 80, 24);
        let k0 = row(&narrow);
        assert!(k0.contains("CT1"), "{what}: {k0}");
        assert!(!k0.contains(" 1 ") && !k0.contains("123"), "{what}: {k0}");
        assert!(
            !narrow.concat().contains("cells   label"),
            "{what}: no header either"
        );
    }
}

#[test]
fn a_susie_candidate_shows_its_pip_and_fold() {
    use super::super::round::Candidate;
    let mut app = app_with_terms(false);
    let mut cand = Candidate::new("CT1", 0.85);
    let e = &mut cand.evidence;
    (e.pip, e.fold, e.nes, e.p, e.q) =
        (Some(1.0), Some(7.0), Some(2.16), Some(0.0001), Some(0.001));
    if let Some(r) = app.round.as_mut() {
        r.clusters[0].candidates = vec![cand.clone()];
    }
    let screen = drawn(&app);
    assert!(screen.contains("PIP"), "a PIP column");
    assert!(screen.contains("×7.0"), "the fold");
    assert!(evidence(&cand).contains("PIP 1.00") && evidence(&cand).contains("×7.0"));
}

#[test]
fn the_susie_setting_toggles_the_stage() {
    let mut a = cli(&[]);
    assert!(settings(&a).contains(&Setting::Susie));
    assert_eq!(Setting::Susie.value(&a), "on");
    Setting::Susie.adjust(&mut a, true);
    assert!(a.no_susie);
    assert_eq!(Setting::Susie.value(&a), "off");
    Setting::Susie.adjust(&mut a, false);
    assert!(!a.no_susie);
}

#[test]
fn a_susie_rounds_cluster_list_shows_the_calls_share_and_fold() {
    use super::super::round::Candidate;
    let mut app = app_with_terms(false);
    if let Some(r) = app.round.as_mut() {
        let mut c = Candidate::new("CT1", 0.85);
        (c.evidence.pip, c.evidence.fold) = (Some(1.0), Some(7.0));
        r.clusters[0].candidates = vec![c];
    }
    let screen = drawn(&app);
    // The screen is 160 cells a row; the list's header row names cells and label.
    let rows: Vec<String> = screen
        .chars()
        .collect::<Vec<_>>()
        .chunks(160)
        .map(|r| r.iter().collect())
        .collect();
    let header = rows
        .iter()
        .find(|r| r.contains("cells") && r.contains("label"))
        .expect("the cluster list's header");
    // The list pane: from its left border to its right one.
    let list: String = header.chars().skip(1).take_while(|&c| c != '│').collect();
    assert!(list.contains("explains") && list.contains("fold"), "{list}");
    assert!(
        !list.contains("PIP"),
        "a PIP beside the call reads as certainty: {list}"
    );
    let k0 = rows.iter().find(|r| r.contains("K0")).unwrap();
    let k0: String = k0.chars().skip(1).take_while(|&c| c != '│').collect();
    assert!(k0.contains("0.85") && k0.contains("×7.0"), "{k0}");
}

#[test]
fn arrows_on_a_form_row_that_changes_nothing_leave_the_round_current() {
    let mut app = app_with_terms(false);
    app.form = Some(Form::default());
    for s in [Setting::Run, Setting::Markers, Setting::Output] {
        app.form.as_mut().unwrap().row = setting_row(&app.args, s);
        press(&mut app, KeyCode::Left);
        press(&mut app, KeyCode::Right);
        assert!(!app.stale(), "{s:?} changed nothing");
    }
    app.form.as_mut().unwrap().row = setting_row(&app.args, Setting::Knn);
    press(&mut app, KeyCode::Right);
    assert!(app.stale(), "knn changed");
}

#[test]
fn enter_on_a_cluster_selects_its_type_only_in_a_tree_on_screen() {
    let panel = [
        ("GENE1".to_string(), "CT1".to_string()),
        ("GENE2".to_string(), "CT2".to_string()),
    ];
    let with_ct2 = || {
        let mut app = app_with_terms(false);
        app.tree = crate::manifest::ontology::panel_tree_on(None, &panel);
        app.round.as_mut().unwrap().clusters[0].label = Some("CT2".into());
        app.focus = Focus::Clusters;
        app
    };
    // On the panel tree, Enter's "on the tree" selects the cluster's type.
    let mut app = with_ct2();
    press(&mut app, KeyCode::Enter);
    press(&mut app, KeyCode::Char('t'));
    assert_eq!(app.tree.label(app.tree.visible()[app.tree_sel]), "CT2");
    assert!(app.focus == Focus::Tree);

    // In the order view the panel tree is hidden: its selection stays.
    let dir = tempfile::tempdir().unwrap();
    let mut app = with_ct2();
    app.data_search = crate::manifest::data_files::SearchPath::new(Some(dir.path()));
    app.reload_order();
    assert!(matches!(app.tree_mode, TreeMode::Order(_)));
    app.focus = Focus::Clusters;
    let before = app.tree_sel;
    press(&mut app, KeyCode::Enter);
    press(&mut app, KeyCode::Char('t'));
    assert_eq!(app.tree_sel, before, "the hidden panel tree is untouched");
    assert!(app.status.contains("panel"), "{}", app.status);
}

#[test]
fn capital_t_leaves_the_order_view_from_the_grid_too() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = figures_app(dir.path());
    press(&mut app, KeyCode::Char('w'));
    assert!(app.figures.as_ref().unwrap().grid.is_some());
    press(&mut app, KeyCode::Char('T'));
    assert!(!app.in_order(), "T leaves the order view");
}

#[test]
fn u_and_k_ask_why_then_backspace_undoes_and_e_waits_for_a_save() {
    use super::super::round::Edit;
    let mut app = app_with_terms(false);
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Char('u'));
    assert!(matches!(
        app.prompt.as_ref().map(|p| &p.pending),
        Some(Pending::Label { label, .. }) if label == UNASSIGNED_LABEL
    ));
    press(&mut app, KeyCode::Enter);
    assert!(
        matches!(&app.edits[..], [Edit::Label { cluster: 0, label, .. }] if label == UNASSIGNED_LABEL)
    );
    press(&mut app, KeyCode::Char('e'));
    assert!(
        app.status.starts_with("save the edits first"),
        "{}",
        app.status
    );

    press(&mut app, KeyCode::Backspace);
    assert!(app.edits.is_empty());
    assert_eq!(app.status, "undid K0's edit");
    press(&mut app, KeyCode::Backspace);
    assert_eq!(app.status, "no edit of K0 to undo");

    press(&mut app, KeyCode::Char('k'));
    assert!(matches!(
        app.prompt.as_ref().map(|p| &p.pending),
        Some(Pending::Keep { cluster: 0 })
    ));
    press(&mut app, KeyCode::Enter);
    assert!(matches!(&app.edits[..], [Edit::Keep { cluster: 0, .. }]));
}

#[test]
fn a_clusters_new_label_replaces_its_earlier_one() {
    use super::super::round::Edit;
    let mut app = app_with_terms(false);
    app.focus = Focus::Clusters;
    app.edits.push(Edit::Markers {
        label: "CT1".into(),
        genes: vec!["GENE1".into()],
        add: false,
        reason: "why".into(),
    });
    for label in ["CT1", "CT1+CT2"] {
        app.ask_label(label.into(), "why".into(), None);
        press(&mut app, KeyCode::Enter);
    }
    assert_eq!(app.edits.len(), 2, "the marker edit and one label");
    assert!(matches!(&app.edits[1], Edit::Label { label, .. } if label == "CT1+CT2"));
    press(&mut app, KeyCode::Backspace);
    assert_eq!(app.edits.len(), 1, "the marker edit stays");
}

#[test]
fn keeping_a_relabelled_cluster_keeps_its_new_label() {
    let mut app = app_with_terms(false);
    app.focus = Focus::Clusters;
    app.ask_label("CT2".into(), "why".into(), None);
    press(&mut app, KeyCode::Enter);
    press(&mut app, KeyCode::Char('k'));
    assert!(app.prompt.as_ref().unwrap().title.contains("keeps CT2"));
    press(&mut app, KeyCode::Enter);
    assert_eq!(app.current_label().as_deref(), Some("CT2"));
    let r = app.round.as_ref().unwrap();
    let saved = super::super::round::decisions(&app.edits, r, None);
    assert_eq!(saved.len(), 1);
    assert_eq!(
        saved[0].label.as_deref(),
        Some("CT2"),
        "the save records it"
    );
}

#[test]
fn close_bracket_visits_the_flagged_clusters_left_to_decide() {
    use super::super::round::{Candidate, ClusterView};
    let mut app = app_with_terms(false);
    let r = app.round.as_mut().unwrap();
    let cluster = |id, label: Option<&str>, candidates| ClusterView {
        id,
        cells: 10,
        label: label.map(Into::into),
        candidates,
        shares: Vec::new(),
        genes: vec![("GENE1".into(), 1.0)],
        terms: Vec::new(),
    };
    let sure = cluster(1, Some("CT1"), vec![Candidate::new("CT1", 1.0)]);
    let open = cluster(2, None, Vec::new());
    r.clusters.extend([sure, open]);
    app.focus = Focus::Clusters;
    // K0 (no candidates) and K2 (no label) are flagged; K1 is not.
    press(&mut app, KeyCode::Char(']'));
    assert_eq!(app.cluster_sel, 2);
    press(&mut app, KeyCode::Char(']'));
    assert_eq!(app.cluster_sel, 0);
    // Deciding both leaves nothing to visit.
    for _ in 0..2 {
        press(&mut app, KeyCode::Char('k'));
        press(&mut app, KeyCode::Enter);
        press(&mut app, KeyCode::Char(']'));
    }
    assert_eq!(app.status, "every flagged cluster is decided");
}

#[test]
fn question_mark_opens_the_guide_and_any_key_closes_it() {
    let mut app = app_with_terms(false);
    press(&mut app, KeyCode::Char('?'));
    assert!(app.help_open);
    press(&mut app, KeyCode::Char('q'));
    assert!(!app.help_open);
    assert!(!app.quit, "the closing key does nothing else");
}

#[test]
fn a_guide_taller_than_the_screen_scrolls_and_says_there_is_more() {
    for (trajectory, width, height) in [(true, 140, 45), (true, 200, 55), (false, 100, 30)] {
        let mut app = app_with_terms(false);
        app.trajectory_tui = trajectory;
        press(&mut app, KeyCode::Char('?'));
        let screen = drawn_in(&app, width, height).concat();
        assert!(
            screen.contains("↓ more"),
            "{width}x{height}: cut, and says so"
        );
        assert!(
            !screen.contains("a marker not enriched"),
            "the marks are below"
        );
        for key in [KeyCode::Down, KeyCode::PageDown, KeyCode::Up, KeyCode::End] {
            press(&mut app, key);
            assert!(app.help_open, "{key:?} scrolls, it does not close");
        }
        let screen = drawn_in(&app, width, height).concat();
        assert!(
            screen.contains("a marker not enriched"),
            "End shows the last"
        );
        assert!(!screen.contains("↓ more"), "nothing more below");
        press(&mut app, KeyCode::Home);
        assert_eq!(app.help_scroll, 0);
        press(&mut app, KeyCode::PageDown);
        assert!(app.help_scroll > 1, "a page");
        press(&mut app, KeyCode::Char('x'));
        assert!(!app.help_open, "any other key closes it");
        press(&mut app, KeyCode::Char('?'));
        assert_eq!(app.help_scroll, 0, "it opens at the top");
    }
}

#[test]
fn gene_keys_add_and_drop_markers_and_backspace_undoes_the_last() {
    use super::super::round::Edit;
    let mut app = app_with_terms(false);
    app.focus = Focus::Genes;
    press(&mut app, KeyCode::Char(' '));
    assert_eq!(app.marked, ["GENE1"]);
    press(&mut app, KeyCode::Char('d'));
    assert_eq!(app.status, "not markers of CT1");
    press(&mut app, KeyCode::Char('a'));
    assert!(matches!(
        app.prompt.as_ref().map(|p| &p.pending),
        Some(Pending::Markers { label, genes, add: true }) if label == "CT1" && genes == &["GENE1"]
    ));
    press(&mut app, KeyCode::Enter);
    assert!(app
        .edits
        .iter()
        .any(|e| matches!(e, Edit::Markers { add: true, .. })));
    // Esc on the "label it too?" question keeps the markers alone.
    app.prompt = None;
    app.focus = Focus::Genes;
    press(&mut app, KeyCode::Backspace);
    assert!(
        app.status
            .starts_with("undid adding GENE1 to CT1's markers"),
        "{}",
        app.status
    );
    assert!(app.edits.is_empty());
}

#[test]
fn on_the_panel_tree_enter_labels_and_plus_mixes_the_marked() {
    let panel = [
        ("GENE1".to_string(), "CT1".to_string()),
        ("GENE2".to_string(), "CT2".to_string()),
    ];
    let mut app = app_with_terms(false);
    app.tree = crate::manifest::ontology::panel_tree_on(None, &panel);
    app.focus = Focus::Tree;
    let at = |app: &App, t: &str| {
        let i = app.tree.node_of(t).unwrap();
        app.tree.visible().iter().position(|&v| v == i).unwrap()
    };
    app.tree_sel = at(&app, "CT2");
    press(&mut app, KeyCode::Enter);
    assert!(matches!(
        app.prompt.as_ref().map(|p| &p.pending),
        Some(Pending::Label { label, .. }) if label == "CT2"
    ));
    app.prompt = None;

    press(&mut app, KeyCode::Char('+'));
    assert!(app.status.starts_with("mark two or more"), "{}", app.status);
    press(&mut app, KeyCode::Char(' '));
    app.tree_sel = at(&app, "CT1");
    press(&mut app, KeyCode::Char(' '));
    press(&mut app, KeyCode::Char('+'));
    assert!(matches!(
        app.prompt.as_ref().map(|p| &p.pending),
        Some(Pending::MixName { parts }) if parts == &["CT1", "CT2"]
    ));
}

#[test]
fn in_the_order_table_two_marks_then_precedes_or_unrelated_ask_why() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = order_app(dir.path(), &["CT1", "CT2"]);
    press(&mut app, KeyCode::Char('>'));
    assert!(app.status.starts_with("mark exactly two"), "{}", app.status);
    press(&mut app, KeyCode::Char(' '));
    press(&mut app, KeyCode::Down);
    press(&mut app, KeyCode::Char(' '));
    assert_eq!(app.tree_marked, ["CT1", "CT2"]);
    for (key, relation) in [('>', Relation::Precedes), ('u', Relation::Unrelated)] {
        press(&mut app, KeyCode::Char(key));
        assert!(matches!(
            app.prompt.as_ref().map(|p| &p.pending),
            Some(Pending::Precedence { from, to, relation: r })
                if from == "CT1" && to == "CT2" && *r == relation
        ));
        app.prompt = None;
    }
    // - and + zoom the figures, as on the cell map: in the table they say so.
    for key in ['-', '+'] {
        press(&mut app, KeyCode::Char(key));
        assert!(app.prompt.is_none(), "{key} states nothing");
        assert!(app.status.contains("zoom"), "{key}: {}", app.status);
        assert!(app.status.contains("u unrelated"), "{key}: {}", app.status);
    }
    let bar = super::super::ui::help(&app);
    assert!(
        bar.contains("u unrelated") && !bar.contains("- unrelated"),
        "{bar}"
    );
    let order: String = super::super::ui::guide(true)
        .filter(|(s, _)| s.starts_with("order view"))
        .flat_map(|(_, keys)| keys)
        .map(|(key, what)| format!("{key} {what}\n"))
        .collect();
    assert!(
        order.contains("u the two marked types are unrelated"),
        "{order}"
    );
}

#[test]
fn in_the_form_enter_picks_the_panel_and_shift_enter_runs_from_any_row() {
    let mut app = app_with_terms(false);
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Char('r'));
    assert!(app.form.is_some());
    app.form.as_mut().unwrap().row = setting_row(&app.args, Setting::Markers);
    press(&mut app, KeyCode::Enter);
    assert!(matches!(app.want_file, Some(FileWant::Markers)));
    app.want_file = None;

    // Plain Enter on another row only says how to run.
    app.form.as_mut().unwrap().row = setting_row(&app.args, Setting::Susie);
    press(&mut app, KeyCode::Enter);
    assert!(app.status.ends_with("starts the pass"), "{}", app.status);
    // Shift+Enter runs from there, where the terminal tells it apart: with
    // no panel yet, it goes back to the panel's row.
    app.shift_enter = true;
    app.key(KeyEvent::new(KeyCode::Enter, KeyModifiers::SHIFT));
    let form = app.form.as_ref().unwrap();
    assert_eq!(form.row, setting_row(&app.args, Setting::Markers));
    assert!(form.note.as_deref().unwrap_or("").contains("marker panel"));
}

#[test]
fn the_annotate_tui_has_no_order_view_and_says_nothing_of_it() {
    let mut app = app_with_terms(false);
    app.focus = Focus::Tree;
    for key in ['t', 'T'] {
        press(&mut app, KeyCode::Char(key));
        assert!(!app.in_order(), "{key} stays on the tree");
    }
    for focus in [Focus::Clusters, Focus::Genes, Focus::Tree] {
        app.focus = focus;
        let bar = super::super::ui::help(&app);
        assert!(!bar.contains("order"), "{focus:?}: {bar}");
    }
    use crate::manifest::family::Kind;
    for kind in [Kind::Run, Kind::Round, Kind::Trajectory] {
        let what = kind.detail(false);
        assert!(
            !what.contains("order view") && !what.contains("figures with"),
            "{what}"
        );
    }
    for (section, keys) in super::super::ui::guide(false) {
        for (key, what) in keys {
            let text = format!("{section} {key} {what}");
            assert!(
                !text.contains("order view") && !text.contains("run the trajectory"),
                "{text}"
            );
        }
    }
}

#[test]
fn the_trajectory_tui_opens_the_order_view_from_the_tree() {
    let mut app = app_with_terms(false);
    app.trajectory_tui = true;
    app.focus = Focus::Tree;
    assert!(super::super::ui::help(&app).contains("t order"));
    press(&mut app, KeyCode::Char('t'));
    assert!(app.in_order());
    assert!(super::super::ui::guide(true).any(|(s, _)| s.starts_with("order view")));
}

#[test]
fn right_on_a_panel_leaf_goes_on_into_the_ontology_at_its_term() {
    use super::super::ontology::Scope;
    let parse = || {
        crate::annotate::celltype_tree::ClTerms::parse(
            TYPES_OBO,
            &crate::annotate::cl_rules::shipped(),
        )
    };
    let panel = [
        ("GENE1".to_string(), "CT1".to_string()),
        ("GENE2".to_string(), "CT2".to_string()),
    ];
    let leaf = |label: &str| {
        let mut app = app_with_terms(false);
        let cl = parse();
        app.tree = crate::manifest::ontology::panel_tree_on(Some(&cl), &panel);
        app.cl = Some(cl);
        app.focus = Focus::Tree;
        let i = app.tree.node_of(label).unwrap();
        app.tree.reveal(i);
        app.select_node(i);
        press(&mut app, KeyCode::Right);
        app
    };
    // CT1 labels the round's cluster: the data's tree, opened on it.
    let app = leaf("CT1");
    let TreeMode::Ontology(v) = &app.tree_mode else {
        panic!("the ontology view");
    };
    assert_eq!(v.scope, Scope::Data);
    assert_eq!(v.selected().unwrap().id, "CL:104");
    assert!(v.is_open("CL:104"));
    // CT2 labels none: the whole ontology around its term.
    let app = leaf("CT2");
    let TreeMode::Ontology(v) = &app.tree_mode else {
        panic!("the ontology view");
    };
    assert_eq!(v.scope, Scope::All);
    assert_eq!(v.selected().unwrap().id, "CL:105");
}

#[test]
fn right_on_an_unfolded_panel_branch_goes_on_into_the_ontology() {
    let cl = crate::annotate::celltype_tree::ClTerms::parse(
        TYPES_OBO,
        &crate::annotate::cl_rules::shipped(),
    );
    let panel = [
        ("GENE1".to_string(), "CT1".to_string()),
        ("GENE2".to_string(), "CT2".to_string()),
        // A type on the parent term too, so the panel tree branches.
        ("GENE3".to_string(), "group".to_string()),
    ];
    let mut app = app_with_terms(false);
    app.tree = crate::manifest::ontology::panel_tree_on(Some(&cl), &panel);
    app.cl = Some(cl);
    app.focus = Focus::Tree;
    let i = (0..app.tree.nodes.len())
        .find(|&i| !app.tree.nodes[i].children.is_empty())
        .expect("a branch");
    app.tree.fold(i, true);
    app.select_node(i);
    press(&mut app, KeyCode::Right);
    assert!(matches!(app.tree_mode, TreeMode::Panel), "first it unfolds");
    assert!(!app.tree.is_folded(i));
    press(&mut app, KeyCode::Right);
    let TreeMode::Ontology(v) = &app.tree_mode else {
        panic!("then on into the ontology");
    };
    let term = app.tree.nodes[i].cl_id.clone().unwrap();
    assert_eq!(v.selected().unwrap().id, term);
    assert!(v.is_open(&term));
}

#[test]
fn the_knn_row_steps_from_the_methods_default_and_back_to_it() {
    let mut a = cli(&[]);
    Setting::Knn.adjust(&mut a, true);
    assert_eq!(a.knn, Some(16));
    Setting::Knn.adjust(&mut a, false);
    assert_eq!(a.knn, None, "back on the default: no --knn");
    a.method = crate::annotate_cmd::AnnotateMethod::Projection;
    Setting::Knn.adjust(&mut a, false);
    assert_eq!(a.knn, Some(29), "projection's default is 30");
}

#[test]
fn a_running_jobs_progress_leaves_an_armed_question_on_screen() {
    let mut app = app_with_terms(false);
    running(&mut app);
    press(&mut app, KeyCode::Char('x'));
    app.tick();
    assert_eq!(app.status, "x again stops it");
    press(&mut app, KeyCode::Char('x'));
    assert!(app.child.is_none());
    // A stopped pass never reached the round.
    assert!(app.status.contains("as it was"), "{}", app.status);

    running(&mut app);
    press(&mut app, KeyCode::Char('q'));
    app.tick();
    assert!(app.status.contains("q again"), "{}", app.status);
    app.stop();
}

#[test]
fn a_long_pinned_label_keeps_its_note_whole() {
    use super::super::round::Edit;
    let mut app = app_with_terms(false);
    let long = "hematopoietic_oligopotent_progenitor_cell";
    let panel = [("GENE1".to_string(), long.to_string())];
    app.tree = crate::manifest::ontology::panel_tree_on(None, &panel);
    app.edits = vec![Edit::Label {
        cluster: 0,
        label: long.into(),
        reason: "why".into(),
    }];
    let screen = drawn(&app);
    assert!(screen.contains("not called (did not pass FDR)"), "{screen}");
}

#[test]
fn a_typed_name_off_the_panel_is_not_said_to_have_failed_fdr() {
    use super::super::round::Edit;
    let mut app = app_with_terms(false);
    let label = |app: &mut App, label: &str| {
        app.edits = vec![Edit::Label {
            cluster: 0,
            label: label.into(),
            reason: "why".into(),
        }];
    };
    label(&mut app, "Test_Custom_Type");
    let screen = drawn(&app);
    assert!(screen.contains("not a panel type"), "{screen}");
    assert!(!screen.contains("did not pass FDR"), "{screen}");
    // A panel type with no share was tested and not called.
    label(&mut app, "CT1");
    assert!(drawn(&app).contains("did not pass FDR"));
}

#[test]
fn a_round_with_no_cluster_expression_says_so_and_calls_no_marker_weak() {
    let mut app = app_with_terms(false);
    if let Some(r) = app.round.as_mut() {
        r.markers.insert("CT1".into(), ["GENE1".into()].into());
        r.clusters[0].genes.clear();
    }
    app.gene_view = GeneView::Markers;
    let screen = drawn(&app);
    assert!(screen.contains("1 markers of CT1"), "{screen}");
    assert!(
        !screen.contains("1 weak"),
        "unmeasured is not weak: {screen}"
    );
    assert!(!screen.contains("weak "), "{screen}");
    assert!(screen.contains("no cluster expression"), "{screen}");
    app.gene_view = GeneView::Specific;
    assert!(drawn(&app).contains("no cluster expression"));
}

#[test]
fn a_panel_leaf_picked_on_the_tree_gives_its_own_evidence_as_the_reason() {
    use super::super::round::Candidate;
    let mut app = app_with_terms(false);
    if let Some(r) = app.round.as_mut() {
        r.clusters[0].candidates = vec![Candidate::new("CT1", 0.18)];
    }
    let leaf = app.tree.node_of("CT1").unwrap();
    app.pick_node(leaf);
    let why = app.prompt.as_ref().unwrap().text.clone();
    assert!(!why.contains("ontology"), "a panel type: {why}");
    assert!(!why.contains("over CT1"), "not pooled with itself: {why}");
    assert!(why.contains("0.18"), "{why}");
}

#[test]
fn a_marker_edits_reason_says_specific_only_as_the_fold_change_does() {
    use legume_numeric::matrix::dense_mat_io::Mat;
    use legume_numeric::matrix::traits::MatWithNames;
    let mut app = app_with_terms(false);
    if let Some(r) = app.round.as_mut() {
        r.expression = Some(super::super::round::Expression::new(MatWithNames {
            rows: vec!["UP".into(), "DOWN".into()],
            cols: vec!["K0".into(), "K1".into()],
            mat: Mat::from_row_slice(2, 2, &[20.0, 1.0, 1.0, 20.0]),
        }));
        r.markers
            .insert("CT1".into(), ["UP".into(), "DOWN".into()].into());
    }
    let why = |app: &mut App, gene: &str, add: bool| {
        app.marked = vec![gene.into()];
        app.ask_markers_of("CT1".into(), add);
        app.prompt.take().map(|p| p.text).unwrap_or_default()
    };
    let up = why(&mut app, "UP", false);
    assert!(!up.contains("not specific"), "UP is up in K0: {up}");
    assert!(up.contains("UP +"), "{up}");
    let down = why(&mut app, "DOWN", false);
    assert!(down.contains("not specific"), "{down}");
}

#[test]
fn a_relabelled_cluster_shows_the_evidence_of_its_new_label() {
    use super::super::round::{Candidate, Edit};
    let cand = |label: &str, fold| {
        let mut c = Candidate::new(label, 0.5);
        (c.evidence.pip, c.evidence.fold) = (Some(1.0), Some(fold));
        c
    };
    let mut app = app_with_terms(false);
    if let Some(r) = app.round.as_mut() {
        r.clusters[0].candidates = vec![cand("CT1", 7.0), cand("CT_2", 3.0)];
    }
    let k0 = |app: &App| {
        drawn(app)
            .chars()
            .collect::<Vec<_>>()
            .chunks(160)
            .map(|r| r.iter().collect::<String>())
            .find(|r| r.contains("K0"))
            .map(|r| {
                r.chars()
                    .skip(1)
                    .take_while(|&c| c != '│')
                    .collect::<String>()
            })
            .unwrap()
    };
    let relabel = |app: &mut App, label: &str| {
        app.edits = vec![Edit::Label {
            cluster: 0,
            label: label.into(),
            reason: "why".into(),
        }];
    };
    relabel(&mut app, "CT_2");
    assert!(k0(&app).contains("×3.0"), "{}", k0(&app));
    // A label is its candidate's however its words are separated.
    relabel(&mut app, "CT 2");
    assert!(k0(&app).contains("×3.0"), "{}", k0(&app));
    relabel(&mut app, UNASSIGNED_LABEL);
    assert!(!k0(&app).contains('×'), "no call, no fold: {}", k0(&app));
}

#[test]
fn the_key_guide_wraps_a_long_description_instead_of_cutting_it() {
    let mut app = app_with_terms(false);
    app.help_open = true;
    let screen = drawn(&app);
    assert!(screen.contains("(◈ hiding)"), "the h row's end");
    assert!(
        screen.contains("pgup/dn home/end "),
        "the key column fits its keys"
    );
}

#[test]
fn a_label_is_recorded_in_its_canonical_form() {
    use super::super::round::Edit;
    let mut app = app_with_terms(false);
    app.prompt = Some(Prompt {
        title: String::new(),
        text: "why".into(),
        pending: Pending::Label {
            cluster: 0,
            label: "CT 2".into(),
            remember: None,
        },
    });
    press(&mut app, KeyCode::Enter);
    assert!(
        matches!(&app.edits[..], [Edit::Label { label, .. }] if label == "CT_2"),
        "{:?}",
        app.edits
    );
}

#[test]
fn the_form_lists_the_rows_of_its_method() {
    use crate::annotate_cmd::AnnotateMethod;
    let mut a = cli(&[]);
    a.method = AnnotateMethod::Enrichment;
    let rows = settings(&a);
    assert!(rows.contains(&Setting::Susie) && rows.contains(&Setting::Go));
    assert!(!rows.contains(&Setting::Projection));
    a.method = AnnotateMethod::Projection;
    a.projection = Some(crate::annotate::args::ProjectionMethod::Score);
    let rows = settings(&a);
    assert!(rows.contains(&Setting::Projection));
    for s in [Setting::Susie, Setting::Go, Setting::NumPerm] {
        assert!(!rows.contains(&s), "{s:?} is not the score's");
    }
    a.projection = Some(crate::annotate::args::ProjectionMethod::Ora);
    assert!(
        settings(&a).contains(&Setting::NumPerm),
        "the ORA's null permutes"
    );
    assert_eq!(
        rows.iter().position(|&s| s == Setting::Method),
        settings(&cli(&[]))
            .iter()
            .position(|&s| s == Setting::Method),
        "changing the method keeps the cursor on its row"
    );
}

#[test]
fn the_projection_row_cycles_the_runs_choice_the_score_and_the_ora() {
    use crate::annotate::args::ProjectionMethod;
    let mut a = cli(&[]);
    assert!(Setting::Projection.value(&a).starts_with("auto"));
    Setting::Projection.adjust(&mut a, true);
    assert_eq!(
        (a.projection, Setting::Projection.value(&a).as_str()),
        (Some(ProjectionMethod::Score), "score")
    );
    Setting::Projection.adjust(&mut a, true);
    assert_eq!(a.projection, Some(ProjectionMethod::Ora));
    Setting::Projection.adjust(&mut a, true);
    assert_eq!(a.projection, None, "back to the run's choice");
}

#[test]
fn the_counts_note_is_checked_again_each_time_the_form_opens() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = app_with_terms(false);
    app.focus = Focus::Clusters;
    app.source = dir.path().join("run.senna.json");
    std::fs::write(
        &app.source,
        r#"{"version": 2, "kind": "svd", "prefix": "run",
            "data": {"input": ["a.zarr.zip"]}}"#,
    )
    .unwrap();
    let note = |app: &App| app.form.as_ref().unwrap().expression.clone();
    press(&mut app, KeyCode::Char('r'));
    assert!(note(&app).is_some_and(|n| n.contains("will fail")));
    press(&mut app, KeyCode::Esc);
    // The counts are reachable now: the form, opened again, says nothing.
    std::fs::write(dir.path().join("a.zarr.zip"), "").unwrap();
    press(&mut app, KeyCode::Char('r'));
    assert_eq!(note(&app), None);
}

#[test]
fn the_counts_note_is_the_enrichments_and_not_shown_for_projection() {
    use crate::annotate_cmd::AnnotateMethod;
    let mut app = app_with_terms(false);
    app.form = Some(Form {
        expression: Some("The raw counts are not here".into()),
        ..Form::default()
    });
    assert!(app.form_summary().0.contains("raw counts"));
    app.args.method = AnnotateMethod::Projection;
    let (summary, _) = app.form_summary();
    assert!(!summary.contains("raw counts"), "{summary}");
}

#[test]
fn a_score_rounds_marker_edits_wait_for_a_new_pass() {
    use super::super::round::Edit;
    let mut app = app_with_terms(false);
    app.round.as_mut().unwrap().by = Some(crate::annotate::rounds::Method::Score);
    app.edits.push(Edit::Markers {
        label: "CT1".into(),
        genes: vec!["GENE1".into()],
        add: true,
        reason: String::new(),
    });
    app.tick();
    assert!(app.status.contains("new pass"), "{}", app.status);
}

#[test]
fn the_cluster_list_shows_a_susie_calls_share_not_its_pip() {
    use super::super::round::Candidate;
    use super::super::ui::{cluster_list_evidence, cluster_list_header};
    // A PIP is a type's inclusion in the cluster's model, near 1 for every
    // type a mixed cluster holds: next to the call it reads as certainty.
    let mut cand = Candidate::new("CT1", 0.20);
    (cand.evidence.pip, cand.evidence.fold) = (Some(1.0), Some(12.0));
    assert_eq!(cluster_list_evidence(true, Some(&cand)), ["0.20", "×12.0"]);
    assert_eq!(
        cluster_list_header(true),
        ["", "cells", "label", "explains", "fold", ""],
        "named apart from the enrichment's share"
    );
    let mut no_susie = cand.clone();
    (no_susie.evidence.pip, no_susie.evidence.fold) = (None, None);
    assert_eq!(
        cluster_list_evidence(true, Some(&no_susie)),
        ["0.20", "—"],
        "a known share is shown; only the fold is missing"
    );
}

/// The app on a four-cell round of two clusters, with a UMAP when `layout`.
fn map_app(dir: &std::path::Path, layout: bool) -> App {
    let l = crate::figures::tests::round(dir, layout);
    let panel = [("GENE1".to_string(), "A".to_string())];
    let tree = crate::manifest::ontology::panel_tree_on(None, &panel);
    let mut app = App::new(cli(&[]), PathBuf::new(), PathBuf::new(), tree);
    app.open(&l.file);
    app
}

/// The cluster the map lights.
fn map_focus(app: &App) -> Option<i32> {
    app.cell_map.as_ref().expect("a map").style.focus
}

#[test]
fn a_round_with_a_layout_shows_its_cells_lit_on_the_selected_cluster() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = map_app(dir.path(), true);
    assert!(
        matches!(app.tree_mode, TreeMode::Cells),
        "the map in the tree pane"
    );
    assert_eq!(map_focus(&app), Some(0));
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Down);
    assert_eq!(map_focus(&app), Some(1), "the light follows the selection");
}

#[test]
fn the_maps_keys_act_from_the_tree_pane_and_capital_v_shows_the_tree() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = map_app(dir.path(), true);
    app.focus = Focus::Tree;
    press(&mut app, KeyCode::Char('c'));
    let style = app.cell_map.as_ref().unwrap().style;
    assert_eq!(
        style.colouring,
        crate::figures::Colouring::Cluster,
        "{}",
        app.status
    );
    press(&mut app, KeyCode::Char('V'));
    assert!(matches!(app.tree_mode, TreeMode::Panel), "V shows the tree");
    press(&mut app, KeyCode::Char('v'));
    assert!(
        matches!(app.tree_mode, TreeMode::Cells),
        "v the cells again"
    );
}

#[test]
fn o_from_the_maps_ontology_goes_back_to_the_map() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = map_app(dir.path(), true);
    app.cl = Some(crate::annotate::celltype_tree::ClTerms::parse(
        TYPES_OBO,
        &crate::annotate::cl_rules::shipped(),
    ));
    app.focus = Focus::Tree;
    press(&mut app, KeyCode::Char('o'));
    assert!(
        matches!(app.tree_mode, TreeMode::Ontology(_)),
        "{}",
        app.status
    );
    let bar = super::super::ui::help(&app);
    assert!(bar.contains("o back to the cells"), "{bar}");
    press(&mut app, KeyCode::Char('o'));
    assert!(
        matches!(app.tree_mode, TreeMode::Cells),
        "the map, not the tree"
    );
    // From the panel tree, o still goes back to it.
    press(&mut app, KeyCode::Char('V'));
    press(&mut app, KeyCode::Char('o'));
    let bar = super::super::ui::help(&app);
    assert!(bar.contains("o panel"), "{bar}");
    press(&mut app, KeyCode::Char('o'));
    assert!(matches!(app.tree_mode, TreeMode::Panel));
    let map = super::super::ui::guide(false)
        .filter(|(s, _)| *s == "tree")
        .flat_map(|(_, keys)| keys)
        .find(|(key, _)| *key == "v / V")
        .map(|(_, what)| what)
        .unwrap();
    assert!(map.contains("o the Cell Ontology"), "{map}");
}

#[test]
fn each_panes_key_line_names_its_keys_and_a_narrow_one_ends_with_more() {
    use super::super::ui::{help, help_in};
    let dir = tempfile::tempdir().unwrap();
    let mut app = map_app(dir.path(), true);
    let wants = [
        (Focus::Clusters, &["u unassign", "⌫ undo", "r form"][..]),
        (
            Focus::Genes,
            &[
                "space mark",
                "X hide by pattern",
                "⌫ undo marker",
                "esc clusters",
            ][..],
        ),
        (Focus::Tree, &["o ontology"][..]),
    ];
    for (focus, want) in wants {
        app.focus = focus;
        let bar = help(&app);
        for w in want {
            assert!(bar.contains(w), "{focus:?} {w}: {bar}");
        }
        assert!(bar.starts_with(" ? keys · "), "{bar}");
        let narrow = help_in(&app, 80);
        assert!(
            ratatui::text::Line::from(narrow.as_str()).width() <= 80,
            "{narrow}"
        );
        assert!(narrow.ends_with(" · ? more"), "{focus:?}: {narrow}");
        let shown: Vec<&str> = narrow.trim().split(" · ").collect();
        let all: Vec<&str> = bar.trim().split(" · ").skip(1).collect();
        assert_eq!(
            shown[..shown.len() - 1],
            all[..shown.len() - 1],
            "the first keys, whole"
        );
    }
    app.focus = Focus::Clusters;
    let wide = help(&app);
    assert_eq!(help_in(&app, 400), wide, "all fit: no more to say");
}

#[test]
fn a_round_without_a_layout_keeps_the_tree_and_v_names_the_senna_command() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = map_app(dir.path(), false);
    assert!(matches!(app.tree_mode, TreeMode::Panel));
    assert!(app.cell_map.is_none());
    app.focus = Focus::Tree;
    press(&mut app, KeyCode::Char('v'));
    assert!(app.status.contains("senna layout"), "{}", app.status);
}

#[test]
fn on_the_map_t_restyles_even_in_the_trajectory_tui_and_esc_leaves() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = map_app(dir.path(), true);
    app.trajectory_tui = true;
    app.focus = Focus::Tree;
    press(&mut app, KeyCode::Char('t'));
    assert!(
        matches!(app.tree_mode, TreeMode::Cells),
        "t stays on the map"
    );
    assert!(app.status.starts_with("labels"), "{}", app.status);
    press(&mut app, KeyCode::Esc);
    assert_eq!(app.focus, Focus::Clusters);
}

#[test]
fn a_new_round_keeps_the_tree_chosen_over_the_map_and_the_maps_look() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = map_app(dir.path(), true);
    app.focus = Focus::Tree;
    press(&mut app, KeyCode::Char('c'));
    let file = app.round.as_ref().unwrap().manifest.clone();
    app.open(&file);
    let style = app.cell_map.as_ref().unwrap().style;
    assert_eq!(
        style.colouring,
        crate::figures::Colouring::Cluster,
        "the colouring kept"
    );
    press(&mut app, KeyCode::Char('V'));
    app.open(&file);
    assert!(matches!(app.tree_mode, TreeMode::Panel), "the tree kept");
}

#[test]
fn capital_s_sorts_the_candidates_by_the_next_method_and_a_number_picks_what_is_shown() {
    use super::super::round::Candidate;
    let mut app = app_with_terms(false);
    let c = &mut app.round.as_mut().unwrap().clusters[0];
    c.candidates = vec![
        Candidate {
            ranks: [Some(0), None, Some(1)],
            ..Candidate::new("CT1", 0.6)
        },
        Candidate {
            ranks: [Some(1), None, Some(0)],
            ..Candidate::new("CT2", 0.3)
        },
    ];
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Char('S'));
    assert!(app.status.contains("score"), "{}", app.status);
    press(&mut app, KeyCode::Char('1'));
    assert!(matches!(
        app.prompt.as_ref().map(|p| &p.pending),
        Some(Pending::Label { label, .. }) if label == "CT2"
    ));
    app.prompt = None;
    press(&mut app, KeyCode::Char('a'));
    assert!(app.status.contains("all"), "{}", app.status);
}

/// The app with cluster 0's candidates ranked by two methods: CT1 the
/// enrichment's top, CT2 the score's, CT3 the enrichment's third.
fn app_with_two_readings() -> App {
    use super::super::round::Candidate;
    let mut app = app_with_terms(false);
    let c = &mut app.round.as_mut().unwrap().clusters[0];
    let cand = |label: &str, share, ranks: [Option<usize>; 3]| {
        let mut c = Candidate {
            ranks,
            ..Candidate::new(label, share)
        };
        (c.evidence.enrichment, c.evidence.score) = (Some(share), ranks[2].map(|_| 0.5));
        c
    };
    c.candidates = vec![
        cand("CT1", 0.6, [Some(0), None, Some(1)]),
        cand("CT3", 0.3, [Some(1), None, None]),
        cand("CT2", 0.1, [Some(2), None, Some(0)]),
    ];
    app.focus = Focus::Clusters;
    app
}

#[test]
fn enter_offers_each_methods_top_choices_first_then_the_tree_and_a_name() {
    let mut app = app_with_two_readings();
    press(&mut app, KeyCode::Enter);
    let labels: Vec<String> = app
        .menu
        .as_ref()
        .expect("the label choices")
        .choices
        .iter()
        .map(|c| c.label.clone())
        .collect();
    assert_eq!(
        labels,
        ["CT1", "CT2", "CT3", "on the tree…", "a name…"],
        "the methods' tops first"
    );
    let screen = drawn(&app);
    assert!(
        screen.contains("on the tree") && screen.contains("a name"),
        "{screen}"
    );
    assert!(
        screen.contains("score"),
        "the highlighted choice's evidence"
    );

    press(&mut app, KeyCode::Char('2'));
    assert!(app.menu.is_none());
    assert!(matches!(
        app.prompt.as_ref().map(|p| &p.pending),
        Some(Pending::Label { label, .. }) if label == "CT2"
    ));
}

#[test]
fn a_name_typed_in_the_choices_labels_the_cluster() {
    let mut app = app_with_two_readings();
    press(&mut app, KeyCode::Enter);
    press(&mut app, KeyCode::Char('n'));
    assert!(matches!(
        app.prompt.as_ref().map(|p| &p.pending),
        Some(Pending::LabelName { .. })
    ));
    app.prompt.as_mut().unwrap().text = "Novel".into();
    press(&mut app, KeyCode::Enter);
    assert!(matches!(
        app.prompt.as_ref().map(|p| &p.pending),
        Some(Pending::Label { label, .. }) if label == "Novel"
    ));
}

#[test]
fn esc_closes_the_choices_and_arrows_move_in_them() {
    let mut app = app_with_two_readings();
    press(&mut app, KeyCode::Enter);
    press(&mut app, KeyCode::Down);
    assert_eq!(app.menu.as_ref().unwrap().sel, 1);
    press(&mut app, KeyCode::Esc);
    assert!(app.menu.is_none());
}

#[test]
fn on_the_map_the_trajectory_tuis_bar_names_capital_t_for_the_order() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = map_app(dir.path(), true);
    app.trajectory_tui = true;
    app.focus = Focus::Tree;
    let bar = super::super::ui::help(&app);
    assert!(
        !bar.contains("t order"),
        "t restyles the labels here: {bar}"
    );
    assert!(bar.contains("T order"), "{bar}");
    press(&mut app, KeyCode::Char('T'));
    assert!(app.in_order(), "T opens the order view");
}

#[test]
fn the_annotate_guide_names_the_maps_grid_and_exports_strip() {
    let text: String = super::super::ui::guide(false)
        .filter(|(section, _)| *section == "tree")
        .flat_map(|(_, keys)| keys.into_iter())
        .map(|(key, what)| format!("{key} {what}\n"))
        .collect();
    for want in [
        "an open branch",
        "w all figures",
        "f the exports strip",
        "M move",
        "D D delete",
    ] {
        assert!(text.contains(want), "{want}: {text}");
    }
}

#[test]
fn the_forms_keys_fit_its_popup_with_or_without_shift_enter() {
    for shift_enter in [false, true] {
        let mut app = app_with_terms(false);
        app.shift_enter = shift_enter;
        press(&mut app, KeyCode::Char('r'));
        assert!(app.form.is_some());
        let screen = drawn(&app);
        assert!(screen.contains("esc/q close"), "shift+enter {shift_enter}");
    }
}

#[test]
fn the_order_views_cluster_keys_name_s_save() {
    let dir = tempfile::tempdir().unwrap();
    let mut app = order_app(dir.path(), &["CT1", "CT2"]);
    app.focus = Focus::Clusters;
    let bar = super::super::ui::help(&app);
    assert!(bar.contains("s save"), "{bar}");
}

#[test]
fn with_the_guide_open_b_closes_it_and_leaves_the_jobs_popup() {
    let mut app = app_with_terms(false);
    running(&mut app);
    app.help_open = true;
    press(&mut app, KeyCode::Char('b'));
    assert!(!app.help_open, "any key closes the guide");
    assert!(!app.progress_hidden, "and does nothing else");
}

#[test]
fn the_trajectory_guide_agrees_with_itself_and_the_keys() {
    let entry = |section: &str, key: &str| {
        super::super::ui::guide(true)
            .find(|(s, _)| s.starts_with(section))
            .and_then(|(_, keys)| keys.into_iter().find(|(k, _)| *k == key))
            .map(|(_, what)| what)
            .unwrap()
    };
    let r = entry("anywhere", "r");
    assert!(r.contains("runs the trajectory"), "{r}");
    let c = entry("figures", "c");
    for colouring in ["cluster", "score call", "margin"] {
        assert!(c.contains(colouring), "{colouring}: {c}");
    }
    assert!(entry("figures", "+ - 0").contains("[ ]"));
}

#[test]
fn the_guide_says_a_save_is_not_stopped_and_t_on_the_map_restyles() {
    let entry = |section: &str, key: &str| {
        super::super::ui::guide(true)
            .find(|(s, _)| s.starts_with(section))
            .and_then(|(_, keys)| keys.into_iter().find(|(k, _)| *k == key))
            .map(|(_, what)| what)
            .unwrap()
    };
    let x = entry("anywhere", "x");
    assert!(x.contains("a save"), "x does not stop a save: {x}");
    let order = super::super::ui::guide(true)
        .map(|(s, _)| s)
        .find(|s| s.starts_with("order view"))
        .unwrap();
    assert!(order.contains("T on the cell map"), "{order}");
}

#[test]
fn the_log_pane_leaves_out_the_svg_fonts_fallback_notes() {
    use log::{Level, LevelFilter, Log, Metadata};
    fn at(target: &str, level: Level) -> Metadata<'_> {
        Metadata::builder().target(target).level(level).build()
    }
    let pane = crate::tui::PaneLogger(LevelFilter::Info);
    assert!(!pane.enabled(&at("usvg::text", Level::Warn)), "font noise");
    assert!(pane.enabled(&at("lupin::figures", Level::Warn)));
    assert!(pane.enabled(&at("lupin::tui", Level::Info)));
    // With -v, everything.
    let verbose = crate::tui::PaneLogger(LevelFilter::Debug);
    assert!(verbose.enabled(&at("usvg::text", Level::Warn)));
}

/// A stand-in for a running `job`.
fn running_job(app: &mut App, job: Job) {
    let child = std::process::Command::new("sleep")
        .arg("30")
        .spawn()
        .unwrap();
    app.child = Some(Running::new(child, job));
}

#[test]
fn a_keys_reply_holds_the_status_against_the_progress_for_a_while() {
    let mut app = app_with_terms(false);
    running_job(&mut app, Job::Save(1));
    press(&mut app, KeyCode::Char('x'));
    app.tick();
    assert!(app.status.contains("finishes on its own"), "{}", app.status);
    // Once the reply has been up its while, the progress takes over again.
    app.replied = app.replied.map(|t| t - KEY_REPLY_HOLD);
    app.tick();
    assert!(app.status.starts_with("save…"), "{}", app.status);
    // A key that says nothing leaves the progress at once.
    press(&mut app, KeyCode::Down);
    app.tick();
    assert!(app.status.starts_with("save…"), "{}", app.status);
    app.stop();
    if let Some(mut r) = app.child.take() {
        let _ = r.child.kill();
        let _ = r.child.wait();
    }

    running(&mut app);
    for code in [KeyCode::Char('r'), KeyCode::Char('A')] {
        app.focus = Focus::Clusters;
        press(&mut app, code);
        app.tick();
        assert!(app.status.contains("wait for"), "{code:?}: {}", app.status);
    }
    app.stop();
}

#[test]
fn the_settings_note_does_not_offer_r_while_a_job_refuses_it() {
    let mut app = app_with_terms(false);
    app.form = Some(Form::default());
    app.form.as_mut().unwrap().row = setting_row(&app.args, Setting::Knn);
    press(&mut app, KeyCode::Right);
    press(&mut app, KeyCode::Esc);
    assert!(drawn(&app).contains("r to re-run"));
    running_job(&mut app, Job::Save(0));
    let screen = drawn(&app);
    assert!(screen.contains("settings changed"), "{screen}");
    assert!(!screen.contains("r to re-run"), "{screen}");
    if let Some(mut r) = app.child.take() {
        let _ = r.child.kill();
        let _ = r.child.wait();
    }
}

#[test]
fn a_question_asked_again_lapses_and_any_other_key_drops_it() {
    let mut app = app_with_terms(false);
    running(&mut app);
    press(&mut app, KeyCode::Char('q'));
    assert!(app.status.contains("q again"), "{}", app.status);
    press(&mut app, KeyCode::Char('b'));
    press(&mut app, KeyCode::Char('q'));
    assert!(!app.quit, "b dropped the question: q asks again");
    assert!(app.child.is_some());
    assert!(app.status.contains("q again"), "{}", app.status);

    // Asked long ago, the question has lapsed: q asks again.
    let age = |app: &mut App| {
        if let Some((_, at)) = &mut app.armed {
            *at -= KEY_REPLY_HOLD;
        }
        if let Some(t) = &mut app.replied {
            *t -= KEY_REPLY_HOLD;
        }
    };
    age(&mut app);
    press(&mut app, KeyCode::Char('q'));
    assert!(!app.quit, "the question had lapsed");
    assert!(app.child.is_some());
    // And with nothing pressed the progress takes the status back.
    age(&mut app);
    app.tick();
    assert!(app.armed.is_none());
    assert!(app.status.starts_with("pass…"), "{}", app.status);
    press(&mut app, KeyCode::Char('x'));
    age(&mut app);
    press(&mut app, KeyCode::Char('x'));
    assert!(app.child.is_some(), "x's question lapses as well");
    app.stop();

    // With no job, a lapsed question leaves the status.
    app.edits.push(super::super::round::Edit::Keep {
        cluster: 0,
        reason: "why".into(),
    });
    press(&mut app, KeyCode::Char('q'));
    assert!(app.status.contains("q again"), "{}", app.status);
    age(&mut app);
    app.tick();
    assert!(!app.status.contains("q again"), "{}", app.status);
}

#[test]
fn q_closes_a_menu_as_esc_does_and_the_hints_say_so() {
    let tmp = tempfile::tempdir().unwrap();
    let mut app = family_app(tmp.path());
    press(&mut app, KeyCode::Char('g'));
    assert!(app.menu.is_some());
    assert!(drawn(&app).contains("esc or q cancels"));
    press(&mut app, KeyCode::Char('q'));
    assert!(app.menu.is_none(), "q closes the runs list");
    assert!(!app.quit, "and does not quit");

    let mut app = app_with_two_readings();
    press(&mut app, KeyCode::Enter);
    press(&mut app, KeyCode::Char('q'));
    assert!(app.menu.is_none(), "q closes the label choices");
    assert!(app.prompt.is_none() && !app.quit);

    press(&mut app, KeyCode::Char('r'));
    assert!(drawn(&app).contains("esc/q close"), "the form's hint");
    let q = super::super::ui::guide(false)
        .find(|(s, _)| s.starts_with("anywhere"))
        .and_then(|(_, keys)| keys.into_iter().find(|(k, _)| k.starts_with('q')))
        .map(|(_, what)| what)
        .unwrap();
    assert!(q.contains("closes a menu"), "{q}");
}

#[test]
fn the_settings_note_goes_once_a_pass_runs_with_them_or_they_are_set_back() {
    let tmp = tempfile::tempdir().unwrap();
    let mut app = family_app(tmp.path());
    app.open_run(pick(tmp.path(), "X.L0"), false);
    let knn = |app: &mut App, right: bool| {
        app.form = Some(Form::default());
        app.form.as_mut().unwrap().row = setting_row(&app.args, Setting::Knn);
        press(app, if right { KeyCode::Right } else { KeyCode::Left });
        app.form = None;
    };
    knn(&mut app, true);
    assert!(app.stale());
    knn(&mut app, false);
    assert!(!app.stale(), "set back: the round is current again");

    knn(&mut app, true);
    running(&mut app);
    assert!(!app.stale(), "the pass running has them");
    assert!(!drawn(&app).contains("settings changed"));
    app.stop();
    assert!(app.stale(), "stopped, it never reached the round");

    // A save's round was made with the settings of the round before it.
    let text = std::fs::read_to_string(tmp.path().join("X.L0.senna.json")).unwrap();
    std::fs::write(tmp.path().join("X.L0.r1.senna.json"), text).unwrap();
    let done = |app: &mut App, job: Job| {
        let mut child = std::process::Command::new("true").spawn().unwrap();
        child.wait().unwrap();
        app.child = Some(Running::new(child, job));
        app.tick();
    };
    done(&mut app, Job::Save(0));
    assert!(app.shown_manifest().ends_with("X.L0.r1.senna.json"));
    assert!(app.stale(), "a save keeps the note");
    app.target = tmp.path().join("X.L9.senna.json");
    done(&mut app, Job::Pass);
    assert!(!app.stale(), "a pass done ran with them");
}

#[test]
fn annotate_opened_on_a_round_shows_that_round() {
    let tmp = tempfile::tempdir().unwrap();
    let mut app = family_app(tmp.path());
    app.trajectory_tui = false;
    app.source = tmp.path().join("X.L0.senna.json");
    app.target = tmp.path().join("X.L0.L2.senna.json");
    app.open_start(true);
    assert!(
        app.shown_manifest().ends_with("X.L0.senna.json"),
        "{}",
        app.status
    );
    assert_eq!(app.round.as_ref().map(|r| r.clusters.len()), Some(2));
    // The run itself, with a round left under the output: that round.
    let mut app = family_app(tmp.path());
    app.trajectory_tui = false;
    app.target = tmp.path().join("X.L0.senna.json");
    app.open_start(false);
    assert!(app.shown_manifest().ends_with("X.L0.senna.json"));
}

#[test]
fn after_a_trajectory_run_g_marks_the_new_trajectory() {
    let tmp = tempfile::tempdir().unwrap();
    let mut app = family_app(tmp.path());
    app.trajectory_tui = true;
    app.open_run(pick(tmp.path(), "X.L0"), false);
    // As a trajectory run leaves its figures, read from X.T1.
    let mut pane = super::super::figure_pane::tests::pane();
    pane.data.manifest = tmp.path().join("X.T1.senna.json");
    app.figures = Some(pane);
    press(&mut app, KeyCode::Char('g'));
    let m = menu(&app).expect("the runs");
    let marked: Vec<&str> = m
        .choices
        .iter()
        .filter(|c| c.label.starts_with('●'))
        .map(|c| c.label.as_str())
        .collect();
    assert_eq!(marked.len(), 1, "{marked:?}");
    assert!(marked[0].starts_with("● X.T1 ·"), "{marked:?}");
}

#[test]
fn clusters_from_the_manifest_leave_the_leiden_rows_fixed_and_passed_over() {
    use super::super::app::{setting_row, Setting};
    let mut app = app_with_terms(false);
    app.fixed_clusters = true;
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Char('r'));
    // The three sit together: down from the row above them lands on the
    // row below, and up comes back.
    let above = setting_row(&app.args, Setting::Knn) - 1;
    let after = setting_row(&app.args, Setting::NumClusters) + 1;
    assert_eq!(after - above, 4);
    app.form.as_mut().unwrap().row = above;
    press(&mut app, KeyCode::Down);
    assert_eq!(app.form.as_ref().unwrap().row, after);
    press(&mut app, KeyCode::Up);
    assert_eq!(app.form.as_ref().unwrap().row, above);
    // On one of them, left and right change nothing.
    let before = format!("{:?}", app.args);
    app.form.as_mut().unwrap().row = setting_row(&app.args, Setting::Resolution);
    press(&mut app, KeyCode::Right);
    press(&mut app, KeyCode::Left);
    assert_eq!(format!("{:?}", app.args), before);
    assert!(!app.stale(), "no settings changed");
    let screen = drawn(&app);
    assert_eq!(screen.matches("from the manifest").count(), 4, "{screen}");
}

#[test]
fn another_runs_trajectory_is_named_over_its_figures_and_in_g() {
    let tmp = tempfile::tempdir().unwrap();
    let mut app = family_app(tmp.path());
    // A pass on X.T1 copies its trajectory section into the new round.
    std::fs::write(
        tmp.path().join("X.T1.L1.senna.json"),
        r#"{"version":2,"kind":"topic","prefix":"X.T1.L1",
            "cluster":{"clusters":"X.clusters.parquet"},
            "annotate":{"argmax":"X.L0.argmax.tsv","source":"X.T1.senna.json"},
            "trajectory":{"prior":"t.tsv","pseudotime":"X.T1.cell_pseudotime.parquet"}}"#,
    )
    .unwrap();
    app.trajectory_tui = true;
    press(&mut app, KeyCode::Char('g'));
    let m = menu(&app).expect("the runs");
    let detail = |name: &str| {
        m.choices
            .iter()
            .find(|c| c.label.contains(&format!("{name} ·")))
            .map(|c| c.detail.clone())
            .unwrap()
    };
    let note = "X.T1's trajectory, on its labels: r runs one for this round";
    assert!(detail("X.T1.L1").contains(note), "{}", detail("X.T1.L1"));
    assert!(!detail("X.T1").contains(note), "{}", detail("X.T1"));

    let tmp = tempfile::tempdir().unwrap();
    let mut app = figures_app(tmp.path());
    app.figures.as_mut().unwrap().data.trajectory_of = Some("X.T1".into());
    let screen = drawn(&app);
    assert!(
        screen.contains("figures from X.T1's trajectory"),
        "{screen}"
    );
}

/// The app with cluster 0's candidates `CT1` … `CT{n}`, the enrichment
/// ranking them in that order.
fn app_with_candidates(n: usize) -> App {
    use super::super::round::Candidate;
    let mut app = app_with_terms(false);
    let c = &mut app.round.as_mut().unwrap().clusters[0];
    c.candidates = (0..n)
        .map(|k| Candidate {
            ranks: [Some(k), None, None],
            ..Candidate::new(&format!("CT{}", k + 1), 0.5 / (k + 1) as f32)
        })
        .collect();
    app.focus = Focus::Clusters;
    app
}

#[test]
fn the_candidates_show_their_first_ten_where_the_screen_allows_and_say_how_many_more() {
    let app = app_with_candidates(14);
    let screen = drawn(&app);
    assert!(screen.contains("CT10"), "the tenth fits: {screen}");
    assert!(!screen.contains("CT11"), "{screen}");
    assert!(screen.contains("4 more"), "{screen}");
}

#[test]
fn tab_reaches_the_candidates_which_scroll_with_their_selection_and_enter_takes_one() {
    let mut app = app_with_candidates(30);
    press(&mut app, KeyCode::Char('a'));
    app.focus = Focus::Tree;
    press(&mut app, KeyCode::Tab);
    assert!(app.focus == Focus::Candidates);
    press(&mut app, KeyCode::End);
    let screen = drawn(&app);
    assert!(screen.contains("CT30"), "scrolled to the last: {screen}");
    assert!(
        !screen.contains("CT2 "),
        "the first rows scrolled off: {screen}"
    );
    press(&mut app, KeyCode::Up);
    press(&mut app, KeyCode::Enter);
    assert!(matches!(
        app.prompt.as_ref().map(|p| &p.pending),
        Some(Pending::Label { label, .. }) if label == "CT29"
    ));
    app.prompt = None;
    // A number takes that row here too; esc goes back to the clusters.
    press(&mut app, KeyCode::Char('2'));
    assert!(matches!(
        app.prompt.as_ref().map(|p| &p.pending),
        Some(Pending::Label { label, .. }) if label == "CT2"
    ));
    app.prompt = None;
    press(&mut app, KeyCode::Esc);
    assert!(app.focus == Focus::Clusters);
    // Another cluster's list starts at its top.
    app.cand_sel = 5;
    app.select_cluster(0);
    assert_eq!(app.cand_sel, 0);
}

#[test]
fn each_candidate_a_number_takes_shows_it_and_the_rest_show_none() {
    let mut app = app_with_candidates(12);
    press(&mut app, KeyCode::Char('a'));
    let screen = drawn(&app);
    assert!(screen.contains("9 CT9 "), "{screen}");
    // No key takes the tenth: it shows no digit, not a cut-off `1`.
    assert!(
        screen.contains("  CT10 ") && !screen.contains("1 CT10"),
        "{screen}"
    );
    // The label choices show their own numbers.
    press(&mut app, KeyCode::Enter);
    assert!(drawn(&app).contains(" 1  CT1 "), "{}", drawn(&app));
}

#[test]
fn the_candidates_value_column_follows_the_method_they_are_sorted_by() {
    use super::super::round::Candidate;
    let mut app = app_with_terms(false);
    let c = &mut app.round.as_mut().unwrap().clusters[0];
    let mut cand = Candidate {
        ranks: [Some(0), Some(0), None],
        ..Candidate::new("CT1", 0.61)
    };
    (
        cand.evidence.enrichment,
        cand.evidence.pip,
        cand.evidence.explained,
    ) = (Some(0.61), Some(1.0), Some(0.33));
    c.candidates = vec![cand];
    app.focus = Focus::Clusters;
    let screen = drawn(&app);
    assert!(
        screen.contains("enrich") && !screen.contains("0.33"),
        "{screen}"
    );
    press(&mut app, KeyCode::Char('S'));
    assert!(app.status.contains("susie"), "{}", app.status);
    let screen = drawn(&app);
    assert!(
        screen.contains("susie") && screen.contains("0.33"),
        "{screen}"
    );
}

#[test]
fn a_rescored_susie_round_keeps_its_call_and_says_what_is_from_before_the_edit() {
    use crate::annotate::rounds::{Method, TypeEvidence};
    let mut app = app_with_terms(false);
    let r = app.round.as_mut().unwrap();
    r.by = Some(Method::Susie);
    let e = TypeEvidence {
        label: "CT1".into(),
        enrichment: Some(0.7),
        pip: Some(1.0),
        fold: Some(6.0),
        explained: Some(0.6),
        ..TypeEvidence::default()
    };
    (r.clusters[0].candidates, r.clusters[0].shares) =
        super::super::round::ranked_view(&[e], Method::Susie);
    let flagged = r.clusters[0].flagged();
    app.show_scores(
        r#"{"scores": {"0": [{"label": "CT1", "share": 0.2, "nes": 1.0, "p": 0.01, "q": 0.05}]}}"#,
    );
    let c = &app.round.as_ref().unwrap().clusters[0];
    assert_eq!(
        (c.candidates[0].pip, c.candidates[0].enrichment),
        (Some(1.0), Some(0.2))
    );
    assert_eq!(
        c.flagged(),
        flagged,
        "no flag for what the rescore did not redo"
    );
    assert!(drawn(&app).contains("before the edit"), "{}", drawn(&app));
    // A label's reason says which of its numbers are from before the edit.
    app.focus = Focus::Clusters;
    press(&mut app, KeyCode::Char('1'));
    let why = &app.prompt.as_ref().expect("asks why").text;
    assert!(why.contains("SuSiE's from before"), "{why}");
}

#[test]
fn a_saved_round_called_by_another_method_says_so() {
    let tmp = tempfile::tempdir().unwrap();
    let mut app = family_app(tmp.path());
    app.source = tmp.path().join("X.L0.senna.json");
    app.open_run(pick(tmp.path(), "X.L0"), false);
    // A SuSiE round whose save, with a marker edit, has no SuSiE call.
    app.round.as_mut().unwrap().by = Some(crate::annotate::rounds::Method::Susie);
    let text = std::fs::read_to_string(tmp.path().join("X.L0.senna.json")).unwrap();
    std::fs::write(tmp.path().join("X.L0.r1.senna.json"), text).unwrap();
    let mut child = std::process::Command::new("true").spawn().unwrap();
    child.wait().unwrap();
    app.child = Some(Running::new(child, Job::Save(0)));
    app.tick();
    assert!(app.shown_manifest().ends_with("X.L0.r1.senna.json"));
    assert!(
        app.status.contains("a save does not redo susie"),
        "{}",
        app.status
    );
}

#[test]
fn a_tree_group_picked_as_a_label_is_not_said_to_have_failed_fdr() {
    use super::super::round::Edit;
    let mut app = app_with_terms(false);
    let panel =
        [("G1", "CT1"), ("G1", "CT2"), ("G2", "CT3")].map(|(g, t)| (g.to_string(), t.to_string()));
    app.tree = crate::manifest::ontology::panel_tree_on(None, &panel);
    let group = (0..app.tree.nodes.len())
        .find(|&i| !app.tree.nodes[i].children.is_empty())
        .expect("CT1 and CT2 share a marker");
    app.edits = vec![Edit::Label {
        cluster: 0,
        label: app.tree.label(group).to_string(),
        reason: "why".into(),
    }];
    let screen = drawn(&app);
    assert!(screen.contains("a tree group"), "{screen}");
    assert!(!screen.contains("did not pass FDR"), "{screen}");
}

#[test]
fn a_statement_on_a_type_below_min_cells_leaves_the_order_menu_to_ask_and_says_so() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::create_dir_all(dir.path().join("lupin")).unwrap();
    std::fs::write(
        dir.path().join("lupin").join("precedence.tsv"),
        "CT1\tCT2\tprecedes\n",
    )
    .unwrap();
    let mut app = order_app(dir.path(), &[]);
    let labels = dir.path().join("labels.tsv");
    std::fs::write(&labels, "c1\tCT1\nc2\tCT1\nc3\tCT2\n").unwrap();
    app.trajectory.labels = Some(labels.to_string_lossy().into());
    app.trajectory.argv = vec!["--min-cells".into(), "2".into()];
    app.reload_order();
    press(&mut app, KeyCode::Char('r'));
    let m = menu(&app).expect("CT2 has 1 cell: nothing orders the nodes");
    assert!(m.question.contains("--min-cells 2"), "{}", m.question);
    assert!(m.question.contains("CT1 precedes CT2"), "{}", m.question);
}
