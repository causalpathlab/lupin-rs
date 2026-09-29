//! Browsing a made-up ontology.

use super::*;

/// root ─┬─ lymph (class) ─┬─ T cell (abbrev "TC") ── T helper
///       │                 └─ B cell
///       └─ stem cell (not on the panel)
const OBO: &str = "format-version: 1.2

[Term]
id: CL:0
name: root cell

[Term]
id: CL:1
name: lymph
subset: cellxgene_subset
is_a: CL:0

[Term]
id: CL:10
name: T cell
synonym: \"TC\" RELATED OMO:0003000 []
subset: cellxgene_subset
is_a: CL:1

[Term]
id: CL:11
name: T helper
is_a: CL:10

[Term]
id: CL:12
name: B cell
is_a: CL:1

[Term]
id: CL:2
name: stem cell
is_a: CL:0
";

fn setup() -> (ClTerms, PanelTree) {
    let cl = ClTerms::parse(OBO, &crate::annotate::cl_rules::shipped());
    let panel: Vec<(String, String)> = [("G1", "T cells"), ("G2", "T helper"), ("G3", "B cells")]
        .iter()
        .map(|(g, t)| ((*g).to_string(), (*t).to_string()))
        .collect();
    let tree = PanelTree::from_ontology(&cl, &panel).unwrap();
    (cl, tree)
}

fn ids(v: &OntologyView) -> Vec<(&str, Role)> {
    v.rows.iter().map(|r| (r.id.as_str(), r.role)).collect()
}

#[test]
fn a_term_shows_its_parents_itself_and_its_children() {
    let (cl, _) = setup();
    let v = OntologyView::at(&cl, "CL:10");
    assert_eq!(
        ids(&v),
        [
            ("CL:1", Role::Parent),
            ("CL:10", Role::Focus),
            ("CL:11", Role::Child)
        ]
    );
    assert_eq!(
        v.selected().unwrap().id,
        "CL:10",
        "the cursor starts on the focus"
    );
}

#[test]
fn walking_down_and_up_keeps_the_way_back() {
    let (cl, _) = setup();
    let mut v = OntologyView::at(&cl, "CL:0");
    v.sel = v.rows.iter().position(|r| r.id == "CL:2").unwrap();
    v.enter(&cl);
    assert_eq!(v.focus, "CL:2", "into a term off the panel");
    v.up(&cl);
    assert_eq!(v.focus, "CL:0");
    assert_eq!(
        v.selected().unwrap().id,
        "CL:2",
        "the cursor is back where we came from"
    );
}

#[test]
fn a_search_lists_hits_and_leaving_it_returns_to_the_focus() {
    let (cl, _) = setup();
    let mut v = OntologyView::at(&cl, "CL:1");
    v.search(&cl, "tc");
    assert_eq!(ids(&v), [("CL:10", Role::Hit)], "abbreviations are found");
    v.up(&cl);
    assert_eq!(v.focus, "CL:1");
    assert!(v.query.is_none());
}

#[test]
fn a_term_is_labelled_by_the_panel_type_on_it_else_by_its_name() {
    let (cl, tree) = setup();
    assert_eq!(term_label(&cl, &tree, "CL:10"), "T_cells");
    assert_eq!(term_label(&cl, &tree, "CL:2"), "stem_cell");
    assert_eq!(term_of(&cl, &tree, "T_cells").as_deref(), Some("CL:10"));
    assert_eq!(term_of(&cl, &tree, "stem_cell").as_deref(), Some("CL:2"));
    assert_eq!(term_of(&cl, &tree, "nothing"), None);
    let lymph: Vec<String> = type_ancestry(&cl, &tree)
        .into_iter()
        .filter(|(_, a)| a.contains("CL:1"))
        .map(|(l, _)| l)
        .collect();
    assert_eq!(
        lymph.len(),
        3,
        "every panel type sits under lymph: {lymph:?}"
    );
}

#[test]
fn only_a_listed_mix_is_a_mix() {
    let root = tempfile::tempdir().unwrap();
    let search = crate::manifest::data_files::SearchPath {
        install: None,
        source: None,
        cache: None,
        user: None,
        project: Some(root.path().join("lupin")),
    };
    let mut m = Mixed::load(&search).unwrap();
    // A `+` is part of many panel names: never read as a mix.
    assert_eq!(m.parts("EMP+HSC"), None);
    assert_eq!(m.parts("CD14+ monocyte"), None);
    assert_eq!(m.parts("HSC"), None);
    assert_eq!(m.parts("HSPC mix"), None, "a name not yet known");
    let file = root.path().join("lupin").join(MIXED);
    m.add("HSPC mix", &["EMP".into(), "HSC".into()], &file)
        .unwrap();
    assert_eq!(m.parts("hspc_mix"), Some(vec!["EMP".into(), "HSC".into()]));
    let again = Mixed::load(&search).unwrap();
    assert_eq!(
        again.parts("HSPC mix"),
        Some(vec!["EMP".into(), "HSC".into()])
    );
    let mut again = again;
    again
        .add("hspc MIX", &["EMP".into(), "HSC".into()], &file)
        .unwrap();
    let text = std::fs::read_to_string(&file).unwrap();
    assert_eq!(
        text.lines()
            .filter(|l| !l.starts_with('#'))
            .collect::<Vec<_>>(),
        ["HSPC mix\tEMP\tHSC"],
        "tab-separated, recorded once"
    );
}

#[test]
fn a_term_with_several_panel_types_is_labelled_the_same_either_way() {
    let cl = ClTerms::parse(OBO, &crate::annotate::cl_rules::shipped());
    let label = |pairs: &[(&str, &str)]| {
        let panel: Vec<(String, String)> = pairs
            .iter()
            .map(|(g, t)| ((*g).to_string(), (*t).to_string()))
            .collect();
        term_label(
            &cl,
            &PanelTree::from_ontology(&cl, &panel).unwrap(),
            "CL:10",
        )
    };
    let a = label(&[("G1", "T cells"), ("G2", "TC"), ("G3", "B cells")]);
    let b = label(&[("G2", "TC"), ("G1", "T cells"), ("G3", "B cells")]);
    assert_eq!(a, b);
}
