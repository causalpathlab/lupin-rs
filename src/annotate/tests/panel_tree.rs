//! [`super::PanelTree`] on a made-up ontology.

use super::*;

/// root ── group1 (class) ─┬─ group2 (class) ─┬─ CT1 cell ── CT1 helper
///                        │                 └─ CT0 cell
///                        └─ group3 (class) ─┬─ CT3
///                                            └─ CT4
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
id: CL:2
name: group2
subset: cellxgene_subset
is_a: CL:1

[Term]
id: CL:3
name: group3
subset: cellxgene_subset
is_a: CL:1

[Term]
id: CL:10
name: CT1 cell
subset: cellxgene_subset
is_a: CL:2

[Term]
id: CL:11
name: CT1 helper
is_a: CL:10

[Term]
id: CL:12
name: CT0 cell
is_a: CL:2

[Term]
id: CL:20
name: CT3
is_a: CL:3

[Term]
id: CL:21
name: CT4
is_a: CL:3
";

fn tree() -> PanelTree {
    let panel: Vec<(String, String)> = [
        ("G1", "CT1 cells"),
        ("G2", "CT1 helper"),
        ("G3", "CT0 cells"),
        ("G4", "CT3"),
        ("G5", "CT4"),
        ("G6", "stray1"),
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
    let helper = t.node_of("CT1 helper").unwrap();
    assert_eq!(
        names(&t, &t.path(helper)),
        ["group1", "group2", "CT1_cell", "CT1_helper"]
    );
    assert_eq!(t.nodes[helper].depth, 3);
    // The unmatched type stands on its own, after the tree.
    let visible = names(&t, &t.visible());
    assert_eq!(visible.first().map(String::as_str), Some("group1"));
    assert_eq!(visible.last().map(String::as_str), Some("stray1"));
}

#[test]
fn a_node_is_called_by_its_panel_type_else_its_name() {
    let t = tree();
    let ct1_cell = t.node_of("CT1_cells").unwrap();
    assert_eq!(
        t.label(ct1_cell),
        "CT1_cells",
        "the panel's label, not CL's"
    );
    let group2 = t.path(ct1_cell)[1];
    assert_eq!(t.label(group2), "group2");
    assert_eq!(
        t.node_of("group2"),
        Some(group2),
        "a class is found by its label"
    );
    let mut under = t.labels_under(group2);
    under.sort_unstable();
    assert_eq!(under, ["CT0_cells", "CT1_cells", "CT1_helper"]);
}

#[test]
fn folding_hides_a_subtree_and_revealing_opens_the_way_back() {
    let mut t = tree();
    let group1 = t.node_of("group1").unwrap();
    let dc = t.node_of("CT4").unwrap();
    let all = t.visible().len();
    t.fold(group1, true);
    assert_eq!(t.visible().len(), 2, "group1 and stray1");
    t.fold(dc, true);
    assert!(!t.is_folded(dc), "a leaf does not fold");
    t.reveal(dc);
    assert!(t.visible().contains(&dc));
    assert_eq!(t.visible().len(), all);
}

#[test]
fn every_type_gets_a_childless_leaf_for_treebh() {
    let t = tree();
    let types: Vec<Box<str>> = ["CT1_cells", "CT1_helper", "CT0_cells", "stray1", "absent"]
        .iter()
        .map(|&s| s.into())
        .collect();
    let tb = t.treebh(&types);
    assert_eq!(tb.leaf[4], None, "a type off the tree hangs under the root");
    let leaves: Vec<usize> = tb.leaf[..4].iter().map(|l| l.unwrap()).collect();
    for &l in &leaves {
        assert!(tb.children[l].is_empty(), "leaf {l} has children");
    }
    // CT1_cells' node has CT1_helper under it, so CT1_cells gets a leaf of its own beside it.
    let t_node = t.node_of("CT1_cells").unwrap() + 1;
    assert!(tb.children[t_node].contains(&leaves[0]));
    assert!(tb.children[t_node].contains(&(t.node_of("CT1_helper").unwrap() + 1)));
    // A childless node holding one type is that type's leaf.
    assert_eq!(leaves[1], t.node_of("CT1_helper").unwrap() + 1);
    assert!(tb.completed(&vec![true; types.len()]).is_ok());
}

#[test]
fn a_type_off_the_ontology_joins_the_class_it_shares_markers_with() {
    // "Pre CT1" matches no term but shares G2 with CT1 helper; "Loner" shares nothing.
    let panel: Vec<(String, String)> = [
        ("G1", "CT1 cells"),
        ("G2", "CT1 helper"),
        ("G3", "CT0 cells"),
        ("G4", "CT3"),
        ("G5", "CT4"),
        ("G2", "Pre CT1"),
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
    let pre = t.node_of("Pre CT1").unwrap();
    let parent = t.nodes[pre].parent.map(|p| t.nodes[p].name.as_str());
    assert!(parent.is_some(), "Pre CT1 is placed inside the tree");
    assert_ne!(parent, Some("Loner"));
    assert_eq!(t.nodes[t.node_of("Loner").unwrap()].parent, None);
}
