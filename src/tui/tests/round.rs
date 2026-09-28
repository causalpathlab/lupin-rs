//! [`super::RoundView`]'s pieces on made-up rounds.

use super::*;

fn candidate(label: &str, share: f32) -> Candidate {
    Candidate {
        label: label.into(),
        share,
        nes: Some(1.0),
        p: Some(0.01),
        q: Some(0.05),
    }
}

fn cluster(id: ClusterId, cells: usize, label: Option<&str>) -> ClusterView {
    ClusterView {
        id,
        cells,
        label: label.map(String::from),
        candidates: vec![candidate("T", 0.6), candidate("B", 0.3)],
        shares: vec![("T".into(), 0.6), ("B".into(), 0.3), ("mono".into(), 0.1)],
        genes: Vec::new(),
    }
}

fn round() -> RoundView {
    RoundView {
        manifest: PathBuf::new(),
        clusters: vec![
            cluster(0, 30, Some("T")),
            cluster(1, 20, Some("B")),
            cluster(2, 5, None),
        ],
        cell_names: Vec::new(),
        cell_clusters: Vec::new(),
        expression: None,
        markers: BTreeMap::from([("T".into(), BTreeSet::from(["CD3E".into()]))]),
        loose_cells: 1,
        decided: BTreeSet::new(),
    }
}

fn relabel(cluster: ClusterId, label: &str) -> Edit {
    Edit::Label {
        cluster,
        label: label.into(),
        reason: "why".into(),
    }
}

#[test]
fn the_latest_edit_of_a_cluster_is_its_label() {
    let r = round();
    let edits = [
        relabel(1, "T"),
        relabel(1, "lymph"),
        relabel(0, UNASSIGNED_LABEL),
    ];
    assert_eq!(r.label_of(1, &edits).as_deref(), Some("lymph"));
    assert_eq!(r.label_of(0, &edits), None, "unassigned by hand");
    assert_eq!(r.label_of(2, &edits), None);
    assert_eq!(
        r.summary(&edits),
        [
            ("lymph".to_string(), 20),
            (UNASSIGNED_LABEL.to_string(), 36)
        ]
    );
    assert_eq!(
        r.summary(&[]),
        [
            ("T".to_string(), 30),
            ("B".to_string(), 20),
            (UNASSIGNED_LABEL.to_string(), 6)
        ]
    );
}

#[test]
fn edits_become_one_decision_per_cluster_after_the_marker_edits() {
    let r = round();
    let edits = [
        relabel(1, "T"),
        Edit::Markers {
            label: "B".into(),
            genes: vec!["MS4A1".into()],
            add: true,
            reason: "specific in K1".into(),
        },
        relabel(1, "lymph"),
    ];
    let d = decisions(&edits, &r);
    assert_eq!(d.len(), 2);
    assert_eq!(d[0].action.as_str(), "markers_add");
    assert_eq!(d[0].features, ["MS4A1"]);
    assert_eq!(d[1].cluster, [1]);
    assert_eq!(d[1].label.as_deref(), Some("lymph"));
    assert_eq!(d[1].evidence.len(), 2, "the candidates it was weighed on");
    assert_eq!(d[1].evidence[0]["nes"], 1.0);
    let q = d[1].evidence[0]["q"].as_f64().unwrap();
    assert!((q - 0.05).abs() < 1e-6, "{q}");
}

#[test]
fn specific_genes_rank_high_and_exclusive_first() {
    // G1 high but everywhere; G2 lower but only in K0; G3 absent.
    let expr = MatWithNames {
        rows: vec!["G1".into(), "G2".into(), "G3".into()],
        cols: vec!["K0".into(), "K1".into()],
        mat: Mat::from_row_slice(3, 2, &[10.0, 10.0, 8.0, 0.0, 0.0, 0.0]),
    };
    let g = specific_genes(&Expression::new(expr));
    let k0: Vec<&str> = g[&0].iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(k0, ["G2", "G1"]);
    assert_eq!(g[&1][0].0, "G1");
}

#[test]
fn a_gene_names_the_types_it_marks() {
    assert_eq!(r_marker(), ["T"]);
}

fn r_marker() -> Vec<String> {
    round()
        .marker_of("CD3E")
        .into_iter()
        .map(String::from)
        .collect()
}

#[test]
fn a_labels_markers_follow_the_edits_in_order() {
    let r = round();
    let markers = |add: bool, genes: &[&str]| Edit::Markers {
        label: "T".into(),
        genes: genes.iter().map(|g| (*g).to_string()).collect(),
        add,
        reason: "why".into(),
    };
    let edits = [
        markers(true, &["CD2", "CD5"]),
        markers(false, &["CD3E", "CD5"]),
    ];
    assert_eq!(r.markers_of("T", &edits), [("CD2".to_string(), true)]);
    let d = decisions(&edits, &r);
    assert_eq!(d[1].action.as_str(), "markers_drop");
    assert_eq!(d[1].features, ["CD3E", "CD5"]);
}

#[test]
fn a_genes_fold_change_is_its_cluster_over_the_rest() {
    let r = RoundView {
        expression: Some(Expression::new(MatWithNames {
            rows: vec!["G".into()],
            cols: vec!["K0".into(), "K1".into()],
            mat: Mat::from_row_slice(1, 2, &[3.0, 1.0]),
        })),
        ..round()
    };
    // The table's mean (2) as pseudo-count: log2((3 + 2) / (1 + 2)).
    let fc = r.fold_change(0, "G").unwrap();
    assert!((fc - (5.0f32 / 3.0).log2()).abs() < 1e-6, "{fc}");
    assert!(r.fold_change(1, "G").unwrap() < 0.0);
    assert_eq!(r.fold_change(0, "absent"), None);
}

#[test]
fn keeping_a_label_is_a_decision_that_undoes_an_earlier_relabel() {
    let r = round();
    let keep = Edit::Keep {
        cluster: 0,
        reason: "the markers agree".into(),
    };
    let edits = [relabel(0, "B"), keep.clone()];
    assert_eq!(
        r.label_of(0, &edits).as_deref(),
        Some("T"),
        "back to the round's label"
    );
    assert_eq!(keep.cluster(), Some(0));
    let d = decisions(&edits, &r);
    assert_eq!(d.len(), 1, "one decision per cluster: the last");
    assert_eq!(d[0].action.as_str(), "keep");
    assert_eq!(d[0].cluster, [0]);
    assert_eq!(d[0].label, None);
    assert_eq!(d[0].rationale, "the markers agree");
}

#[test]
fn a_rounds_history_says_which_clusters_are_decided() {
    let root = tempfile::tempdir().unwrap();
    let p = root.path().join("h.json");
    std::fs::write(
        &p,
        r#"{"1": [{"action": "keep"}], "3": [{"action": "label"}], "4": [], "x": [{}]}"#,
    )
    .unwrap();
    let d = decided_in(&p.to_string_lossy()).unwrap();
    assert_eq!(d, BTreeSet::from([1, 3]));
}
