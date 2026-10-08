//! [`super::RoundView`]'s pieces on made-up rounds.

use super::*;

fn candidate(label: &str, share: f32) -> Candidate {
    Candidate {
        label: label.into(),
        share,
        susie: None,
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
        candidates: vec![candidate("CT1", 0.6), candidate("CT2", 0.3)],
        shares: vec![
            ("CT1".into(), 0.6),
            ("CT2".into(), 0.3),
            ("CT3".into(), 0.1),
        ],
        genes: Vec::new(),
        terms: Vec::new(),
    }
}

fn round() -> RoundView {
    RoundView {
        manifest: PathBuf::new(),
        clusters: vec![
            cluster(0, 30, Some("CT1")),
            cluster(1, 20, Some("CT2")),
            cluster(2, 5, None),
        ],
        cell_names: Vec::new(),
        cell_clusters: Vec::new(),
        expression: None,
        markers: BTreeMap::from([("CT1".into(), BTreeSet::from(["GENE1".into()]))]),
        loose_cells: 1,
        decided: BTreeSet::new(),
        rescorable: false,
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
        relabel(1, "CT1"),
        relabel(1, "CT5"),
        relabel(0, UNASSIGNED_LABEL),
    ];
    assert_eq!(r.label_of(1, &edits).as_deref(), Some("CT5"));
    assert_eq!(r.label_of(0, &edits), None, "unassigned by hand");
    assert_eq!(r.label_of(2, &edits), None);
    assert_eq!(
        r.summary(&edits),
        [("CT5".to_string(), 20), (UNASSIGNED_LABEL.to_string(), 36)]
    );
    assert_eq!(
        r.summary(&[]),
        [
            ("CT1".to_string(), 30),
            ("CT2".to_string(), 20),
            (UNASSIGNED_LABEL.to_string(), 6)
        ]
    );
}

#[test]
fn edits_become_one_decision_per_cluster_after_the_marker_edits() {
    let r = round();
    let edits = [
        relabel(1, "CT1"),
        Edit::Markers {
            label: "CT2".into(),
            genes: vec!["GENE4".into()],
            add: true,
            reason: "specific in K1".into(),
        },
        relabel(1, "CT5"),
    ];
    let d = decisions(&edits, &r, None);
    assert_eq!(d.len(), 2);
    assert_eq!(d[0].action.as_str(), "markers_add");
    assert_eq!(d[0].features, ["GENE4"]);
    assert_eq!(d[1].cluster, [1]);
    assert_eq!(d[1].label.as_deref(), Some("CT5"));
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
    assert_eq!(r_marker(), ["CT1"]);
}

fn r_marker() -> Vec<String> {
    round()
        .marker_of("GENE1")
        .into_iter()
        .map(String::from)
        .collect()
}

#[test]
fn a_labels_markers_follow_the_edits_in_order() {
    let r = round();
    let markers = |add: bool, genes: &[&str]| Edit::Markers {
        label: "CT1".into(),
        genes: genes.iter().map(|g| (*g).to_string()).collect(),
        add,
        reason: "why".into(),
    };
    let edits = [
        markers(true, &["GENE2", "GENE3"]),
        markers(false, &["GENE1", "GENE3"]),
    ];
    assert_eq!(r.markers_of("CT1", &edits), [("GENE2".to_string(), true)]);
    let d = decisions(&edits, &r, None);
    assert_eq!(d[1].action.as_str(), "markers_drop");
    assert_eq!(d[1].features, ["GENE1", "GENE3"]);
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
    let edits = [relabel(0, "CT2"), keep.clone()];
    assert_eq!(
        r.label_of(0, &edits).as_deref(),
        Some("CT1"),
        "back to the round's label"
    );
    assert_eq!(keep.cluster(), Some(0));
    let d = decisions(&edits, &r, None);
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

#[test]
fn adding_a_gene_again_moves_it_from_an_unsaved_addition() {
    let add = |label: &str, genes: &[&str]| Edit::Markers {
        label: label.into(),
        genes: genes.iter().map(|g| (*g).to_string()).collect(),
        add: true,
        reason: "why".into(),
    };
    let mut edits = vec![add("X", &["G1", "G2"]), add("Z", &["G3"])];
    let from = take_added(&mut edits, 0, &["G1".into(), "G3".into()], "Y");
    assert_eq!(from, ["X", "Z"]);
    assert_eq!(edits.len(), 1, "Z's addition is left empty and dropped");
    let r = round();
    assert_eq!(
        r.markers_of("X", &edits)
            .into_iter()
            .filter(|(_, added)| *added)
            .map(|(g, _)| g)
            .collect::<Vec<_>>(),
        ["G2"]
    );
    // Adding to the same type again leaves its own addition alone.
    assert!(take_added(&mut edits, 0, &["G2".into()], "X").is_empty());
    assert_eq!(edits.len(), 1);
}

#[test]
fn rescored_scores_replace_a_clusters_and_swap_back() {
    let preview = serde_json::json!({"scores": {"0": [
        {"label": "CT4", "share": 0.8, "nes": 2.1, "p": 1e-4, "q": 1e-3},
        {"label": "CT1", "share": 0.2, "nes": 1.2, "p": 0.01, "q": 0.04},
        {"label": "CT2", "share": 0.0, "nes": null, "p": 0.9, "q": 1.0}
    ]}});
    let scores = parse_scores(&preview).unwrap();
    let (candidates, shares) = &scores[&0];
    assert_eq!(
        candidates
            .iter()
            .map(|c| c.label.as_str())
            .collect::<Vec<_>>(),
        ["CT4", "CT1"],
        "a type with no share is no candidate"
    );
    assert_eq!(candidates[0].nes, Some(2.1));
    assert_eq!(shares.len(), 3, "every type keeps its share");

    let mut r = round();
    let old = r.swap_scores(scores);
    assert_eq!(r.clusters[0].candidates[0].label, "CT4");
    assert_eq!(r.clusters[1].candidates[0].label, "CT1", "K1 not rescored");
    r.swap_scores(old);
    assert_eq!(r.clusters[0].candidates[0].label, "CT1");
    assert_eq!(r.clusters[0].shares.len(), 3);
    assert!(parse_scores(&serde_json::json!({"scores": null})).is_none());
}

#[test]
fn a_gene_being_saved_is_moved_by_a_drop_not_by_rewriting_the_save() {
    let add = |label: &str, genes: &[&str]| Edit::Markers {
        label: label.into(),
        genes: genes.iter().map(|g| (*g).to_string()).collect(),
        add: true,
        reason: "why".into(),
    };
    // The first edit is being saved.
    let mut edits = vec![add("X", &["G1"]), add("Z", &["G1", "G2"])];
    let from = take_added(&mut edits, 1, &["G1".into()], "Y");
    assert_eq!(from, ["X", "Z"]);
    assert_eq!(edits[0], add("X", &["G1"]), "the save's edit untouched");
    assert_eq!(edits[1], add("Z", &["G2"]));
    assert!(matches!(
        &edits[2],
        Edit::Markers { label, genes, add: false, .. } if label == "X" && genes == &["G1"]
    ));
}

#[test]
fn a_decisions_evidence_is_the_rounds_own_when_rescored_scores_are_shown() {
    let mut r = round();
    let recorded = r.swap_scores(Scores::from([(
        0,
        (vec![candidate("CT4", 0.9)], vec![("CT4".into(), 0.9)]),
    )]));
    let d = decisions(&[relabel(0, "CT1")], &r, Some(&recorded));
    assert_eq!(d[0].evidence[0]["term"], "CT1", "the round's own, not CT4");
}

#[test]
fn a_susie_rounds_candidates_follow_its_call() {
    let row = |v: [f32; 3]| -> Vec<(String, f32)> {
        ["A", "B", "C"]
            .iter()
            .map(|t| t.to_string())
            .zip(v)
            .collect()
    };
    // A and B both certain; B explains more. C barely enters.
    let pip = row([1.0, 1.0, 0.01]);
    let effect = row([0.95, 2.6, 0.1]);
    let explained = row([0.2, 0.5, 0.0]);
    let c = susie_candidates(&pip, &effect, &explained);
    let labels: Vec<&str> = c.iter().map(|c| c.label.as_str()).collect();
    assert_eq!(
        labels,
        ["B", "A"],
        "PIPs tie: the larger share leads; C is out"
    );
    assert_eq!(c[0].share, 0.5, "the share is the deviance explained");
    let s = c[0].susie.as_ref().unwrap();
    assert_eq!(s.pip, 1.0);
    assert!((s.fold - 2.6f32.exp()).abs() < 1e-4);
}

/// A SuSiE candidate with a PIP and an explained share.
fn susie_candidate(label: &str, pip: f32, share: f32) -> Candidate {
    Candidate {
        susie: Some(SusieEvidence { pip, fold: 4.0 }),
        ..candidate(label, share)
    }
}

#[test]
fn a_susie_cluster_is_flagged_when_its_call_is_unsure_or_close() {
    let with = |c: Vec<Candidate>| ClusterView {
        candidates: c,
        ..cluster(0, 10, Some("A"))
    };
    let clear = with(vec![
        susie_candidate("A", 1.0, 0.8),
        susie_candidate("B", 1.0, 0.3),
    ]);
    assert!(!clear.flagged(), "a sure call well ahead");
    let unsure = with(vec![susie_candidate("A", 0.9, 0.8)]);
    assert!(unsure.flagged(), "PIP under 0.95");
    let close = with(vec![
        susie_candidate("A", 1.0, 0.8),
        susie_candidate("B", 1.0, 0.45),
    ]);
    assert!(
        close.flagged(),
        "a runner-up explains at least half as much"
    );
    let none = ClusterView {
        label: None,
        ..with(vec![susie_candidate("A", 1.0, 0.8)])
    };
    assert!(none.flagged(), "no call");
}
