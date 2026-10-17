//! [`super::RoundView`]'s pieces on made-up rounds.

use super::*;

fn candidate(label: &str, share: f32) -> Candidate {
    let mut c = Candidate::new(label, share);
    (c.evidence.nes, c.evidence.p, c.evidence.q) = (Some(1.0), Some(0.01), Some(0.05));
    c
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
        by: Default::default(),
        figures: Default::default(),
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
        {"label": "CT2", "share": 0.0, "nes": null, "p": 0.9, "q": 1.0},
        {"label": "unassigned", "share": 0.5, "nes": null, "p": null, "q": null}
    ]}});
    let mut r = round();
    let scores = r.rescored(None, &parse_scores(&preview).unwrap());
    let (candidates, shares) = &scores[&0];
    assert_eq!(
        candidates
            .iter()
            .map(|c| c.label.as_str())
            .collect::<Vec<_>>(),
        ["CT4", "CT1"],
        "a type with no share is no candidate, nor unassigned"
    );
    assert_eq!(candidates[0].nes, Some(2.1));
    assert_eq!(shares.len(), 3, "every rescored type has its share");
    assert_eq!(
        candidates[0].ranks,
        [Some(0), None, None],
        "an enrichment round is ranked by the rescored enrichment"
    );

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
    use crate::annotate::rounds::{Method, TypeEvidence};
    let e = |label: &str, pip: f32, effect: f32, explained: f32| TypeEvidence {
        label: label.into(),
        pip: Some(pip),
        fold: Some(effect.exp()),
        explained: Some(explained),
        ..TypeEvidence::default()
    };
    // A and B both certain; B explains more. C barely enters.
    let rows = [
        e("A", 1.0, 0.95, 0.2),
        e("B", 1.0, 2.6, 0.5),
        e("C", 0.01, 0.1, 0.0),
    ];
    let (c, _) = ranked_view(&rows, Method::Susie);
    let labels: Vec<&str> = c.iter().map(|c| c.label.as_str()).collect();
    assert_eq!(
        labels,
        ["B", "A"],
        "PIPs tie: the larger share leads; C is out"
    );
    assert_eq!(c[0].share, 0.5, "the share is the log-likelihood explained");
    let s = c[0].susie().unwrap();
    assert_eq!(s.pip, 1.0);
    assert!((s.fold - 2.6f32.exp()).abs() < 1e-4);
}

/// A SuSiE candidate with a PIP and an explained share.
fn susie_candidate(label: &str, pip: f32, share: f32) -> Candidate {
    let mut c = candidate(label, share);
    (c.evidence.pip, c.evidence.fold) = (Some(pip), Some(4.0));
    c
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

#[test]
fn a_score_round_lists_each_clusters_types_by_share() {
    use crate::annotate::marker_score::toy;
    let dir = tempfile::tempdir().unwrap();
    let (manifest, panel) = toy::write_run(dir.path());
    let loaded = run::load(&manifest).unwrap();
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    crate::manifest::annotate::annotate_by_marker_score(&toy::args(&panel, &out), &loaded).unwrap();
    let view = RoundView::load(&crate::manifest::run::annotated_path(&loaded.file, &out)).unwrap();
    for (c, want) in view.clusters.iter().zip(["A", "B"]) {
        let top = c.candidates.first().expect("a candidate");
        assert_eq!(
            (top.label.as_str(), top.share),
            (want, 1.0),
            "cluster {}",
            c.id
        );
    }
}

#[test]
fn a_merged_score_round_lists_its_types_from_the_cells_calls() {
    use crate::annotate::marker_score::toy;
    use crate::manifest::rounds::{run_relabel, RelabelArgs};
    let dir = tempfile::tempdir().unwrap();
    let (manifest, panel) = toy::write_run(dir.path());
    let loaded = run::load(&manifest).unwrap();
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    crate::manifest::annotate::annotate_by_marker_score(&toy::args(&panel, &out), &loaded).unwrap();
    let round = crate::manifest::run::annotated_path(&loaded.file, &out);
    let decisions = dir.path().join("d.jsonl");
    std::fs::write(
        &decisions,
        r#"{"clusters": [0, 1], "action": "merge", "label": "AB", "rationale": "one", "decided_by": "user"}"#,
    )
    .unwrap();
    let next = dir.path().join("run.L2");
    run_relabel(&RelabelArgs {
        from: round.to_string_lossy().into(),
        decisions: decisions.to_string_lossy().into(),
        out: Some(next.to_string_lossy().into()),
        ..RelabelArgs::default()
    })
    .unwrap();
    let view = RoundView::load(&dir.path().join("run.L2.senna.json")).unwrap();
    let merged = &view.clusters[0];
    let mut shares: Vec<(&str, f32)> = merged
        .candidates
        .iter()
        .map(|c| (c.label.as_str(), c.share))
        .collect();
    shares.sort_by(|a, b| a.0.cmp(b.0));
    assert_eq!(
        shares,
        [("A", 0.5), ("B", 0.5)],
        "the merged cluster's cells, by their calls"
    );
}

#[test]
fn a_score_round_whose_table_is_gone_opens_without_shares() {
    use crate::annotate::marker_score::toy;
    let dir = tempfile::tempdir().unwrap();
    let (manifest, panel) = toy::write_run(dir.path());
    let loaded = run::load(&manifest).unwrap();
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    crate::manifest::annotate::annotate_by_marker_score(&toy::args(&panel, &out), &loaded).unwrap();
    std::fs::remove_file(format!(
        "{out}{}",
        crate::annotate::outputs::MARKER_SCORE_PARQUET
    ))
    .unwrap();
    let view = RoundView::load(&crate::manifest::run::annotated_path(&loaded.file, &out)).unwrap();
    assert!(view.clusters.iter().all(|c| c.candidates.is_empty()));
}

#[test]
fn call_shares_table_each_clusters_calls_with_uncalled_cells_in_its_size() {
    let clusters = [Some(3), Some(3), Some(3), Some(3), None, Some(5)];
    let calls = [
        Some("A"),
        Some("A"),
        Some("unassigned"),
        None,
        Some("A"),
        Some("B"),
    ];
    let t = crate::annotate::rounds::call_shares(&clusters, &calls);
    assert_eq!(t.rows, [3, 5]);
    assert_eq!(t.cols, ["A", "B"], "unassigned is not a call");
    assert_eq!(t.row(3).unwrap(), [0.5, 0.0]);
    assert_eq!(t.row(5).unwrap(), [0.0, 1.0]);
}

#[test]
fn a_rounds_candidates_carry_the_scores_share_beside_the_enrichments() {
    use crate::annotate::marker_score::toy;
    let dir = tempfile::tempdir().unwrap();
    let (manifest, panel) = toy::write_run(dir.path());
    let loaded = run::load(&manifest).unwrap();
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    crate::manifest::annotate::annotate_by_marker_score(&toy::args(&panel, &out), &loaded).unwrap();
    let view = RoundView::load(&crate::manifest::run::annotated_path(&loaded.file, &out)).unwrap();
    let top = &view.clusters[0].candidates[0];
    assert_eq!((top.label.as_str(), top.score), ("A", Some(1.0)));
}

fn ranked(label: &str, share: f32, ranks: [Option<usize>; 3]) -> Candidate {
    Candidate {
        ranks,
        ..candidate(label, share)
    }
}

#[test]
fn the_list_shows_ten_then_any_methods_top_below_them_and_all_on_asking() {
    use crate::annotate::rounds::Method;
    let mut cands: Vec<Candidate> = (0..12)
        .map(|i| {
            ranked(
                &format!("T{i:02}"),
                0.5 - 0.01 * i as f32,
                [Some(i), None, None],
            )
        })
        .collect();
    cands[11].ranks[2] = Some(0);
    let names = |by: Method, all: bool| -> Vec<String> {
        listed(&cands, by, all)
            .iter()
            .map(|c| c.label.clone())
            .collect()
    };
    let short = names(Method::Enrichment, false);
    let first_ten: Vec<String> = (0..10).map(|i| format!("T{i:02}")).collect();
    assert_eq!(short[..10], first_ten[..]);
    assert_eq!(
        short[10..],
        ["T11".to_string()],
        "the score's top is always listed"
    );
    assert_eq!(names(Method::Enrichment, true).len(), 12);
    assert_eq!(
        names(Method::Score, false)[0],
        "T11",
        "by the score, its top first"
    );
}

#[test]
fn a_rounds_candidates_know_where_each_method_ranks_them() {
    let dir = tempfile::tempdir().unwrap();
    let round = crate::manifest::annotate::projection_score_tests::round_with_both(dir.path());
    let view = RoundView::load(&round.file).unwrap();
    let k0 = &view.clusters[0].candidates;
    let ranks = |l: &str| k0.iter().find(|c| c.label == l).unwrap().ranks;
    assert_eq!(ranks("A"), [Some(0), None, Some(0)]);
    assert_eq!(ranks("B"), [Some(1), None, None]);
}

#[test]
fn the_choices_are_every_methods_top_three_firsts_first() {
    let cands = vec![
        ranked("A", 0.5, [Some(0), None, Some(1)]),
        ranked("C", 0.2, [Some(2), None, None]),
        ranked("B", 0.3, [Some(1), None, Some(0)]),
        ranked("D", 0.1, [Some(3), None, Some(2)]),
        ranked("E", 0.0, [Some(4), None, Some(3)]),
    ];
    let labels: Vec<&str> = choices(&cands).iter().map(|c| c.label.as_str()).collect();
    assert_eq!(labels, ["A", "B", "C", "D"]);
}

#[test]
fn a_rescore_keeps_the_rounds_method_and_the_readings_it_did_not_redo() {
    use crate::annotate::rounds::Method;
    let mut r = round();
    r.by = Some(Method::Susie);
    // SuSiE calls CT1, sure and well ahead; the score agrees.
    let mut ct1 = susie_candidate("CT1", 1.0, 0.6);
    (ct1.evidence.explained, ct1.evidence.score) = (Some(0.6), Some(0.7));
    let mut ct2 = susie_candidate("CT2", 0.99, 0.1);
    ct2.evidence.explained = Some(0.1);
    let (candidates, shares) = ranked_view(&[ct1.evidence, ct2.evidence], Method::Susie);
    (r.clusters[0].candidates, r.clusters[0].shares) = (candidates, shares.clone());
    let flagged = r.clusters[0].flagged();

    // The edit moves the enrichment to CT2.
    let preview = serde_json::json!({"scores": {"0": [
        {"label": "CT2", "share": 0.9, "nes": 2.0, "p": 1e-4, "q": 1e-3},
        {"label": "CT1", "share": 0.1, "nes": 0.5, "p": 0.2, "q": 0.4}
    ]}});
    let scores = r.rescored(None, &parse_scores(&preview).unwrap());
    let (c, rescored_shares) = &scores[&0];
    let labels: Vec<&str> = c.iter().map(|c| c.label.as_str()).collect();
    assert_eq!(labels, ["CT1", "CT2"], "still in SuSiE's order");
    assert_eq!((c[0].pip, c[0].score), (Some(1.0), Some(0.7)), "kept");
    assert_eq!(
        (c[0].enrichment, c[0].nes),
        (Some(0.1), Some(0.5)),
        "rescored"
    );
    assert_eq!(c[0].share, 0.6, "the share is still SuSiE's");
    assert_eq!(c[1].ranks[Method::Enrichment.index()], Some(0));
    assert_eq!(rescored_shares, &shares, "SuSiE's shares stand");
    assert!(!scores.contains_key(&1), "K1 was not rescored");

    // Rescored again, the kept readings are still the round's own.
    let recorded = r.swap_scores(scores);
    let again = r.rescored(Some(&recorded), &parse_scores(&preview).unwrap());
    assert_eq!(again[&0].0[0].pip, Some(1.0));
    assert_eq!(
        r.clusters[0].flagged(),
        flagged,
        "no flag for what was not redone"
    );
}
