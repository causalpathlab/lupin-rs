//! Completing an incomplete prior by connectivity.

use super::*;

/// Symmetric connectivity over `n` types from `(a, b, value)` pairs.
fn conn(n: usize, pairs: &[(usize, usize, f64)]) -> DMatrix<f64> {
    let mut m = DMatrix::zeros(n, n);
    for &(a, b, v) in pairs {
        m[(a, b)] = v;
        m[(b, a)] = v;
    }
    m
}

#[test]
fn an_edgeless_type_joins_along_its_strongest_connection() {
    // 0 → 1 is the prior; 2 has no edge and links to 1 more than to 0.
    let c = conn(3, &[(0, 1, 0.9), (1, 2, 0.4), (0, 2, 0.1)]);
    let (host, edges) = hosts(&[true; 3], &[Some(0), Some(0), None], &c);
    assert_eq!(host, vec![None, None, Some(1)]);
    assert_eq!(edges, vec![(1, 2)]);
}

#[test]
fn a_chain_of_edgeless_types_takes_the_placed_type_it_hangs_from() {
    // 2 hangs from 0; 3 only reaches the prior through 2.
    let c = conn(4, &[(0, 1, 0.9), (0, 2, 0.3), (2, 3, 0.5)]);
    let (host, edges) = hosts(&[true; 4], &[Some(0), Some(0), None, None], &c);
    assert_eq!(host, vec![None, None, Some(0), Some(0)]);
    assert_eq!(edges, vec![(0, 2), (2, 3)]);
}

#[test]
fn two_prior_components_are_never_joined_and_unconnected_types_stay_out() {
    // Components 0 and 1; 2 bridges them, more to 3's side; 4 links nowhere;
    // 5 is not a node type.
    let c = conn(6, &[(0, 2, 0.2), (2, 3, 0.6), (2, 5, 0.9)]);
    let is_node = [true, true, true, true, true, false];
    let component = [Some(0), Some(0), None, Some(1), None, None];
    let (host, edges) = hosts(&is_node, &component, &c);
    assert_eq!(host, vec![None, None, Some(3), None, None, None]);
    assert_eq!(edges, vec![(3, 2)], "only the stronger side, one edge");
}

/// A round whose labels are coarse groups (`G1`, `G2`) and which keeps its
/// fine labels (`CT1` to `CT3`), or keeps none.
fn coarse_round(root: &std::path::Path, keeps_fine: bool) -> Loaded {
    use crate::manifest::run::{RunKind, RunManifest};
    let write = |file: &str, labels: [&str; 3]| {
        let text: String = (["a", "b", "c"].iter().zip(labels))
            .map(|(c, l)| format!("{c}\t{l}\n"))
            .collect();
        std::fs::write(root.join(file), text).unwrap();
    };
    write("run.argmax.tsv", ["G1", "G1", "G2"]);
    write("fine.argmax.tsv", ["CT1", "CT2", "CT3"]);
    let mut m = RunManifest::new(RunKind::Topic, "run");
    m.annotate.argmax = Some("run.argmax.tsv".into());
    if keeps_fine {
        m.annotate.fine_argmax = Some("fine.argmax.tsv".into());
    }
    let file = root.join("run.senna.json");
    m.save(&file).unwrap();
    load(&file.to_string_lossy()).unwrap()
}

#[test]
fn a_trajectory_orders_the_fine_labels_when_asked_or_when_a_root_is_one() {
    let dir = tempfile::tempdir().unwrap();
    let loaded = coarse_round(dir.path(), true);
    let roots = |rs: &[&str]| rs.iter().map(|r| Box::from(*r)).collect::<Vec<Box<str>>>();
    let level_of = |level, rs: &[&str]| label_source(&loaded, level, None, &roots(rs)).unwrap();
    let (path, level) = level_of(None, &[]);
    assert!(path.ends_with("run.argmax.tsv"), "{path}");
    assert_eq!(level, Some(Level::Coarse));
    assert_eq!(level_of(None, &["G1"]).1, Some(Level::Coarse));
    let (path, level) = level_of(None, &["CT1"]);
    assert!(path.ends_with("fine.argmax.tsv"), "a fine root: {path}");
    assert_eq!(level, Some(Level::Fine));
    assert_eq!(level_of(Some(Level::Fine), &[]).1, Some(Level::Fine));
    assert_eq!(
        level_of(Some(Level::Coarse), &["CT1"]).1,
        Some(Level::Coarse),
        "a level asked for stands; its root check says the rest"
    );
    let given = label_source(&loaded, None, Some("own.tsv"), &roots(&["CT1"])).unwrap();
    assert_eq!(given, ("own.tsv".to_string(), None));
}

#[test]
fn fine_labels_a_round_does_not_keep_are_refused_by_name() {
    let dir = tempfile::tempdir().unwrap();
    let loaded = coarse_round(dir.path(), false);
    let e = label_source(&loaded, Some(Level::Fine), None, &[]).unwrap_err();
    assert!(format!("{e:#}").contains("fine_argmax"), "{e:#}");
    let (path, level) = label_source(&loaded, None, None, &[Box::from("CT1")]).unwrap();
    assert!(path.ends_with("run.argmax.tsv"), "{path}");
    assert_eq!(level, Some(Level::Coarse));
}

#[test]
fn coarse_roots_never_read_the_fine_labels() {
    let dir = tempfile::tempdir().unwrap();
    let loaded = coarse_round(dir.path(), true);
    std::fs::remove_file(dir.path().join("fine.argmax.tsv")).unwrap();
    let (path, level) = label_source(&loaded, None, None, &[Box::from("G1")]).unwrap();
    assert!(path.ends_with("run.argmax.tsv"), "{path}");
    assert_eq!(level, Some(Level::Coarse));
}
