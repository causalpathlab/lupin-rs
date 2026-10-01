//! Remembering an alias picked in the TUI; the GO terms view and setting.

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
        q: None,
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
fn m_reaches_the_go_terms_only_when_the_round_scored_them() {
    let mut app = app_with_terms(true);
    let views: Vec<GeneView> = (0..3)
        .map(|_| {
            app.pane_key(KeyCode::Char('m'));
            app.gene_view
        })
        .collect();
    assert!(views == [GeneView::Markers, GeneView::Terms, GeneView::Specific]);

    let mut app = app_with_terms(false);
    app.pane_key(KeyCode::Char('m'));
    app.pane_key(KeyCode::Char('m'));
    assert!(app.gene_view == GeneView::Specific);
}

#[test]
fn go_terms_are_read_not_edited() {
    let mut app = app_with_terms(true);
    app.gene_view = GeneView::Terms;
    assert_eq!(app.listed_genes(), ["TERM1", "TERM2"]);
    app.pane_key(KeyCode::Char(' '));
    app.pane_key(KeyCode::Char('x'));
    assert!(app.marked.is_empty());
    assert!(!app.hidden.hides("TERM1"));
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
