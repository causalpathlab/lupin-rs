//! [`super::PanelTree`] on a made-up ontology.

use super::*;

/// root ── blood (class) ─┬─ lymph (class) ─┬─ T cell ── T helper
///                        │                 └─ B cell
///                        └─ myeloid (class) ─┬─ mono
///                                            └─ DC
const OBO: &str = "format-version: 1.2

[Term]
id: CL:0
name: root cell

[Term]
id: CL:1
name: blood
subset: cellxgene_subset
is_a: CL:0

[Term]
id: CL:2
name: lymph
subset: cellxgene_subset
is_a: CL:1

[Term]
id: CL:3
name: myeloid
subset: cellxgene_subset
is_a: CL:1

[Term]
id: CL:10
name: T cell
subset: cellxgene_subset
is_a: CL:2

[Term]
id: CL:11
name: T helper
is_a: CL:10

[Term]
id: CL:12
name: B cell
is_a: CL:2

[Term]
id: CL:20
name: mono
is_a: CL:3

[Term]
id: CL:21
name: DC
is_a: CL:3
";

fn tree() -> PanelTree {
    let panel: Vec<(String, String)> = [
        ("G1", "T cells"),
        ("G2", "T helper"),
        ("G3", "B cells"),
        ("G4", "mono"),
        ("G5", "DC"),
        ("G6", "Osteoblast"),
    ]
    .iter()
    .map(|(g, t)| ((*g).to_string(), (*t).to_string()))
    .collect();
    PanelTree::from_ontology(
        &ClTerms::parse(OBO, &crate::annotate::cl_rules::shipped()),
        &panel,
    )
    .unwrap()
}

fn names(t: &PanelTree, ids: &[usize]) -> Vec<String> {
    ids.iter().map(|&i| t.nodes[i].name.clone()).collect()
}

#[test]
fn types_hang_under_shared_classes_nearest_first() {
    let t = tree();
    let helper = t.node_of("T helper").unwrap();
    assert_eq!(
        names(&t, &t.path(helper)),
        ["blood", "lymph", "T_cell", "T_helper"]
    );
    assert_eq!(t.nodes[helper].depth, 3);
    // The unmatched type stands on its own, after the tree.
    let visible = names(&t, &t.visible());
    assert_eq!(visible.first().map(String::as_str), Some("blood"));
    assert_eq!(visible.last().map(String::as_str), Some("Osteoblast"));
}

#[test]
fn a_node_is_called_by_its_panel_type_else_its_name() {
    let t = tree();
    let t_cell = t.node_of("T_cells").unwrap();
    assert_eq!(t.label(t_cell), "T_cells", "the panel's label, not CL's");
    let lymph = t.path(t_cell)[1];
    assert_eq!(t.label(lymph), "lymph");
    assert_eq!(
        t.node_of("lymph"),
        Some(lymph),
        "a class is found by its label"
    );
    let mut under = t.labels_under(lymph);
    under.sort_unstable();
    assert_eq!(under, ["B_cells", "T_cells", "T_helper"]);
}

#[test]
fn folding_hides_a_subtree_and_revealing_opens_the_way_back() {
    let mut t = tree();
    let blood = t.node_of("blood").unwrap();
    let dc = t.node_of("DC").unwrap();
    let all = t.visible().len();
    t.fold(blood, true);
    assert_eq!(t.visible().len(), 2, "blood and Osteoblast");
    t.fold(dc, true);
    assert!(!t.is_folded(dc), "a leaf does not fold");
    t.reveal(dc);
    assert!(t.visible().contains(&dc));
    assert_eq!(t.visible().len(), all);
}

#[test]
fn every_type_gets_a_childless_leaf_for_treebh() {
    let t = tree();
    let types: Vec<Box<str>> = ["T_cells", "T_helper", "B_cells", "Osteoblast", "absent"]
        .iter()
        .map(|&s| s.into())
        .collect();
    let tb = t.treebh(&types);
    assert_eq!(tb.leaf[4], None, "a type off the tree hangs under the root");
    let leaves: Vec<usize> = tb.leaf[..4].iter().map(|l| l.unwrap()).collect();
    for &l in &leaves {
        assert!(tb.children[l].is_empty(), "leaf {l} has children");
    }
    // T_cells' node has T_helper under it, so T_cells gets a leaf of its own beside it.
    let t_node = t.node_of("T_cells").unwrap() + 1;
    assert!(tb.children[t_node].contains(&leaves[0]));
    assert!(tb.children[t_node].contains(&(t.node_of("T_helper").unwrap() + 1)));
    // A childless node holding one type is that type's leaf.
    assert_eq!(leaves[1], t.node_of("T_helper").unwrap() + 1);
    assert!(tb.completed(&vec![true; types.len()]).is_ok());
}

#[test]
fn a_type_off_the_ontology_joins_the_class_it_shares_markers_with() {
    // "Pre T" matches no term but shares G2 with T helper; "Loner" shares nothing.
    let panel: Vec<(String, String)> = [
        ("G1", "T cells"),
        ("G2", "T helper"),
        ("G3", "B cells"),
        ("G4", "mono"),
        ("G5", "DC"),
        ("G2", "Pre T"),
        ("G9", "Loner"),
    ]
    .iter()
    .map(|(g, t)| ((*g).to_string(), (*t).to_string()))
    .collect();
    let t = PanelTree::from_ontology(
        &ClTerms::parse(OBO, &crate::annotate::cl_rules::shipped()),
        &panel,
    )
    .unwrap();
    let pre = t.node_of("Pre T").unwrap();
    let parent = t.nodes[pre].parent.map(|p| t.nodes[p].name.as_str());
    assert!(parent.is_some(), "Pre T is placed inside the tree");
    assert_ne!(parent, Some("Loner"));
    assert_eq!(t.nodes[t.node_of("Loner").unwrap()].parent, None);
}
