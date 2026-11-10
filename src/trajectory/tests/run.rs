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

/// Fine labels over six cells (one with none) and the coarse label of
/// each: `CT1` and the unlabelled cell under `G1`; `CT2` and `CT3` under
/// `G2`; `CT4` mostly under `G2`, one cell under `G1`.
struct FineOverCoarse {
    names: Vec<Box<str>>,
    group: Vec<usize>,
    is_node: Vec<bool>,
    coarse: Vec<&'static str>,
}

fn fine_over_coarse() -> FineOverCoarse {
    let fine = [
        "CT1",
        "CT1",
        enrichment::UNASSIGNED_LABEL,
        "CT2",
        "CT3",
        "CT3",
        "CT4",
        "CT4",
        "CT4",
    ];
    let coarse = ["G1", "G1", "G1", "G2", "G2", "G2", "G1", "G2", "G2"];
    let (names, group) = encode_groups(&fine);
    let is_node = names
        .iter()
        .map(|n| n.as_ref() != enrichment::UNASSIGNED_LABEL)
        .collect();
    FineOverCoarse {
        names,
        group,
        is_node,
        coarse: coarse.to_vec(),
    }
}

#[test]
fn on_fine_labels_a_coarse_root_stands_for_the_fine_types_mostly_under_it() {
    let FineOverCoarse {
        names,
        group,
        is_node,
        coarse,
    } = fine_over_coarse();
    let roots = |rs: &[&str]| rs.iter().map(|r| Box::from(*r)).collect::<Vec<Box<str>>>();
    let on = |rs: &[&str]| fine_roots(&names, &group, &is_node, &coarse, &roots(rs));
    assert_eq!(
        on(&["G1"]).unwrap(),
        roots(&["CT1"]),
        "the unlabelled cell takes no part"
    );
    assert_eq!(
        on(&["G2"]).unwrap(),
        roots(&["CT2", "CT3", "CT4"]),
        "CT4 is mostly under G2"
    );
    assert_eq!(
        on(&["CT2"]).unwrap(),
        roots(&["CT2"]),
        "a fine root stands as it is"
    );
    assert_eq!(
        on(&["G2", "CT3"]).unwrap(),
        roots(&["CT2", "CT3", "CT4"]),
        "each once"
    );
    let e = on(&["X"]).unwrap_err().to_string();
    assert!(e.contains("X") && e.contains("G1"), "neither label: {e}");
    assert!(
        !e.contains(enrichment::UNASSIGNED_LABEL),
        "never offered as a label: {e}"
    );
    assert!(on(&[enrichment::UNASSIGNED_LABEL]).is_err(), "never a root");
}

#[test]
fn a_coarse_root_over_no_fine_node_is_refused_saying_why() {
    let FineOverCoarse {
        names,
        group,
        mut is_node,
        coarse,
    } = fine_over_coarse();
    let ct1 = names.iter().position(|n| n.as_ref() == "CT1").unwrap();
    is_node[ct1] = false;
    let e = fine_roots(&names, &group, &is_node, &coarse, &[Box::from("G1")])
        .unwrap_err()
        .to_string();
    assert!(e.contains("G1") && e.contains("--min-cells"), "{e}");
}

#[test]
fn a_candidate_pair_runs_from_the_earlier_type_to_the_later() {
    let candidate = || EdgeCheck {
        from: 0,
        to: 1,
        verdict: Verdict::Candidate,
        order_agreement: f32::NAN,
    };
    // Type 0 is late (median 0.8), type 1 early (median 0.2).
    let by_type = vec![vec![0.7, 0.8, 0.9], vec![0.1, 0.2, 0.75]];
    let mut c = candidate();
    let same = [Some(0), Some(0)];
    direct(&mut c, &by_type, &same);
    assert_eq!((c.from, c.to), (1, 0), "the earlier median first");
    assert_eq!(c.order_agreement, 1.0, "all of the late type beyond 0.2");

    // A type no root reaches has no pseudotime: the pair keeps no direction.
    let mut c = candidate();
    direct(&mut c, &[vec![0.5], vec![]], &same);
    assert_eq!((c.from, c.to), (0, 1));
    assert!(c.order_agreement.is_nan());

    // Two components measure pseudotime from their own roots, and equal
    // medians say nothing: neither pair is given a direction.
    let mut c = candidate();
    direct(&mut c, &by_type, &[Some(0), Some(1)]);
    assert_eq!((c.from, c.to), (0, 1));
    assert!(c.order_agreement.is_nan());
    let mut c = candidate();
    direct(&mut c, &[vec![1.0, 1.0], vec![1.0, 1.0]], &same);
    assert!(c.order_agreement.is_nan());
}

#[test]
fn a_prior_edge_is_supported_weak_or_unsupported_by_its_edges() {
    use crate::trajectory::type_connectivity::Connectivity;
    let verdict = |edges, expected| verdict_of(&Connectivity::pair(edges, expected), 0, 1, 0.1);
    assert_eq!(verdict(400.0, 1450.0), Verdict::Supported);
    assert_eq!(verdict(72.0, 1450.0), Verdict::Weak, "joined, thinly");
    assert_eq!(verdict(18.0, 908.0), Verdict::Unsupported);
    assert_eq!(
        verdict(2.0, 9.3),
        Verdict::Unsupported,
        "a ratio of 0.2 on two edges is no evidence"
    );
    assert_eq!(Verdict::parse("weak"), Some(Verdict::Weak));
    assert_eq!(Verdict::Weak.as_str(), "weak");
}

#[test]
fn a_fine_runs_coarse_root_is_read_from_the_rounds_own_labels_by_cell() {
    let dir = tempfile::tempdir().unwrap();
    let loaded = coarse_round(dir.path(), true);
    // The run's cells in another order than the label files'.
    let cells: Vec<Box<str>> = ["c", "a", "b"].map(Box::from).into();
    let (names, group) = encode_groups(&["CT3", "CT1", "CT2"]);
    let inputs = Inputs {
        level: Some(Level::Fine),
        geometry: Mat::zeros(3, 1),
        cells_of: Vec::new(),
        is_node: vec![true; names.len()],
        cells,
        names,
        group,
    };
    let root = |r: &str| roots_of(&loaded, &inputs, &[Box::from(r)]).unwrap();
    assert_eq!(root("G1"), [Box::from("CT1"), Box::from("CT2")]);
    assert_eq!(root("G2"), [Box::from("CT3")]);
    assert_eq!(
        root("CT2"),
        [Box::from("CT2")],
        "a fine root reads nothing more"
    );
    let coarse = Inputs {
        level: Some(Level::Coarse),
        ..inputs
    };
    assert_eq!(
        roots_of(&loaded, &coarse, &[Box::from("G1")]).unwrap(),
        [Box::from("G1")]
    );
}
