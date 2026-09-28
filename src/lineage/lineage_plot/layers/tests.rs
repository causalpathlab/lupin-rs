use super::{medoid, select_labeled_nodes, NodeLabels};
use std::collections::HashMap;

/// `(node names, roles, cell types, node → pixel position)`.
type Fixture = (
    Vec<Box<str>>,
    Vec<Box<str>>,
    Vec<Box<str>>,
    HashMap<Box<str>, (f32, f32)>,
);

/// 6 nodes: root (CT_root), 2 terminal CT1 in a tight pair + 1 far outlier CT1,
/// 1 internal CT2, 1 unassigned.
fn fixture() -> Fixture {
    let names: Vec<Box<str>> = ["n0", "n1", "n2", "n3", "n4", "n5"]
        .iter()
        .map(|s| (*s).into())
        .collect();
    let roles: Vec<Box<str>> = [
        "root", "terminal", "terminal", "terminal", "internal", "internal",
    ]
    .iter()
    .map(|s| (*s).into())
    .collect();
    let types: Vec<Box<str>> = ["CT_root", "CT1", "CT1", "CT1", "CT2", super::UNASSIGNED]
        .iter()
        .map(|s| (*s).into())
        .collect();
    let pos: HashMap<Box<str>, (f32, f32)> = [
        ("n0", (0.0, 0.0)),
        ("n1", (10.0, 0.0)),
        ("n2", (11.0, 0.0)),
        ("n3", (99.0, 0.0)), // far outlier
        ("n4", (0.0, 5.0)),
        ("n5", (1.0, 1.0)),
    ]
    .iter()
    .map(|(k, v)| ((*k).into(), *v))
    .collect();
    (names, roles, types, pos)
}

#[test]
fn per_type_labels_one_node_per_type_and_never_unassigned() {
    let (n, r, t, p) = fixture();
    let keep = select_labeled_nodes(&n, &r, &t, &p, NodeLabels::PerType);
    // root + CT1 rep + CT2 rep = 3; `unassigned` excluded.
    assert_eq!(keep.len(), 3, "got {keep:?}");
    assert!(keep.contains(&0), "root must always be labeled");
    assert!(
        !keep.contains(&5),
        "an `unassigned` node must never be labeled"
    );
}

/// The root already names its type — do not print it twice.
#[test]
fn per_type_does_not_duplicate_the_roots_type() {
    let (n, mut r, mut t, p) = fixture();
    // Make n4 a second CT_root node (internal).
    t[4] = "CT_root".into();
    r[4] = "internal".into();
    let keep = select_labeled_nodes(&n, &r, &t, &p, NodeLabels::PerType);
    let root_type: Vec<_> = keep.iter().filter(|&&i| &*t[i] == "CT_root").collect();
    assert_eq!(root_type.len(), 1, "CT_root labeled twice: {keep:?}");
    assert!(
        keep.contains(&0),
        "the surviving root-type label is the root's"
    );
}

/// The representative is the medoid of the type's TERMINAL nodes, so a far
/// outlier never wins over a tight cluster.
#[test]
fn per_type_representative_is_the_medoid_not_an_outlier() {
    let (n, r, t, p) = fixture();
    let keep = select_labeled_nodes(&n, &r, &t, &p, NodeLabels::PerType);
    let ct1: Vec<usize> = keep.iter().copied().filter(|&i| &*t[i] == "CT1").collect();
    assert_eq!(ct1.len(), 1);
    assert!(
        ct1[0] == 1 || ct1[0] == 2,
        "picked the far outlier n3 instead of the tight pair: {ct1:?}"
    );
}

#[test]
fn other_modes_scale_as_documented() {
    let (n, r, t, p) = fixture();
    let count = |m| select_labeled_nodes(&n, &r, &t, &p, m).len();
    assert_eq!(count(NodeLabels::None), 0);
    assert_eq!(count(NodeLabels::Root), 1);
    assert_eq!(count(NodeLabels::Terminal), 4); // 3 terminals + root
}

#[test]
fn medoid_picks_the_central_member() {
    let (n, _, _, p) = fixture();
    // n1(10,0), n2(11,0), n3(99,0) -> medoid is n1 or n2, never the outlier n3.
    let m = medoid(&[1, 2, 3], &n, &p).unwrap();
    assert!(m == 1 || m == 2, "medoid chose the outlier: {m}");
    assert_eq!(medoid(&[4], &n, &p), Some(4), "singleton is its own medoid");
    assert_eq!(medoid(&[], &n, &p), None);
}
