//! Digest, decisions and history for [`super`].

use super::*;

fn some(v: &[u32]) -> Vec<Option<ClusterId>> {
    v.iter().map(|&x| Some(x)).collect()
}

fn labels(v: &[&str]) -> Vec<Option<String>> {
    v.iter().map(|s| Some((*s).to_string())).collect()
}

fn decisions(jsonl: &str) -> Vec<Decision> {
    jsonl
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).unwrap())
        .collect()
}

#[test]
fn digest_takes_the_majority_label_and_ranks_calls() {
    let clusters = some(&[0, 0, 0, 1]);
    let labs = labels(&["CT1", "CT1", "CT2", "CT2"]);
    let ev = Evidence {
        q: Some(Table {
            rows: vec![0, 1],
            cols: vec!["CT1".into(), "CT2".into(), "unassigned".into()],
            values: vec![0.01, 0.2, 0.0, 0.5, 0.001, 0.0],
        }),
        support: Some(Table {
            rows: vec![0],
            cols: vec!["CT1".into(), "CT2".into()],
            values: vec![0.9, 0.1],
        }),
        ..Evidence::default()
    };
    let d = digest(&clusters, &labs, &ev);
    assert_eq!(d[&0].size, 3);
    assert_eq!(d[&0].label.as_deref(), Some("CT1"));
    let top: Vec<&str> = d[&0].calls.iter().map(|c| c.label.as_str()).collect();
    assert_eq!(top, ["CT1", "CT2"], "unassigned is not a call");
    assert_eq!(d[&0].calls[0].support, Some(0.9));
    // No support for cluster 1: ranked by q alone.
    assert_eq!(d[&1].calls[0].label, "CT2");
    assert_eq!(d[&1].calls[0].support, None);
}

#[test]
fn the_digest_file_is_keyed_by_bare_integer_strings() {
    let d = digest(&some(&[12]), &labels(&["CT1"]), &Evidence::default());
    let v = serde_json::to_value(&d).unwrap();
    assert!(v.get("12").is_some(), "{v}");
}

#[test]
fn cluster_ids_parse_with_or_without_a_prefix() {
    assert_eq!(parse_cluster_id("12"), Some(12));
    assert_eq!(parse_cluster_id("K12"), Some(12));
    assert_eq!(parse_cluster_id("C3"), Some(3));
    assert_eq!(parse_cluster_id("CT"), None);
}

#[test]
fn label_and_merge_update_cells_and_record_history() {
    let mut clusters = some(&[0, 0, 1, 2, 2]);
    let mut labs = labels(&["CT1", "CT1", "CT2", "CT3", "CT2"]);
    let mut ds = decisions(
        r#"
        {"cluster": "1", "action": "label", "label": "CT4", "rationale": "r1", "decided_by": "user", "evidence": [{"kind": "marker", "term": "CT4"}, 2, 3, 4]}
        {"clusters": [0, "K2"], "action": "merge", "label": "CT1", "rationale": "r2", "decided_by": "agent_proposed_user_accepted"}
        "#,
    );
    let next = next_cluster_id(&clusters, &History::new());
    assert_eq!(next, 3);
    let h = apply(
        &mut ds,
        &mut clusters,
        &mut labs,
        next,
        "r1.senna.json",
        "T",
    )
    .unwrap();

    assert_eq!(clusters, some(&[3, 3, 1, 3, 3]));
    assert_eq!(labs, labels(&["CT1", "CT1", "CT4", "CT1", "CT1"]));
    assert_eq!(h[&1][0].label.as_deref(), Some("CT4"));
    assert_eq!(
        h[&1][0].evidence.len(),
        3,
        "evidence is cut to the top items"
    );
    assert_eq!(h[&3][0].merged_from, Some(vec![0, 2]));
    assert_eq!(h[&0][0].merged_into, Some(3));
    assert_eq!(h[&3][0].timestamp, "T");
    assert_eq!(
        ds[0].timestamp.as_deref(),
        Some("T"),
        "the log keeps the stamp"
    );
}

#[test]
fn a_merged_id_is_never_handed_out_again() {
    let mut history = History::new();
    history.insert(
        7,
        vec![HistoryEntry {
            round: "r".into(),
            action: Action::Merge,
            label: None,
            rationale: "r".into(),
            decided_by: DecidedBy::User,
            timestamp: "T".into(),
            evidence: vec![],
            merged_from: Some(vec![4, 9]),
            merged_into: None,
            split_from: None,
            features: None,
        }],
    );
    assert_eq!(next_cluster_id(&some(&[7, 1]), &history), 10);
}

#[test]
fn decisions_without_reasons_or_on_unknown_clusters_are_refused() {
    let run = |jsonl: &str| {
        let mut clusters = some(&[0, 1]);
        let mut labs = labels(&["CT1", "CT2"]);
        apply(&mut decisions(jsonl), &mut clusters, &mut labs, 2, "r", "T")
    };
    let no_reason = r#"{"cluster": 0, "action": "keep", "decided_by": "user"}"#;
    assert!(run(no_reason).is_err());
    let unknown = r#"{"cluster": 5, "action": "keep", "rationale": "r", "decided_by": "user"}"#;
    assert!(run(unknown).is_err());
    let twice = "{\"cluster\": 0, \"action\": \"keep\", \"rationale\": \"r\", \"decided_by\": \"user\"}\n\
                 {\"cluster\": 0, \"action\": \"label\", \"label\": \"CT3\", \"rationale\": \"r\", \"decided_by\": \"user\"}";
    assert!(run(twice).is_err());
    let split = r#"{"cluster": 0, "action": "split", "rationale": "r", "decided_by": "user"}"#;
    assert!(run(split).is_err());
    let unlabelled = r#"{"cluster": 0, "action": "label", "rationale": "r", "decided_by": "user"}"#;
    assert!(run(unlabelled).is_err());
}

#[test]
fn newer_history_goes_first() {
    let e = |r: &str| HistoryEntry {
        round: r.into(),
        action: Action::Keep,
        label: None,
        rationale: "r".into(),
        decided_by: DecidedBy::User,
        timestamp: "T".into(),
        evidence: vec![],
        merged_from: None,
        merged_into: None,
        split_from: None,
        features: None,
    };
    let older = History::from([(0, vec![e("r1")])]);
    let newer = History::from([(0, vec![e("r2")])]);
    let rounds: Vec<String> = prepend_history(older, newer)[&0]
        .iter()
        .map(|e| e.round.clone())
        .collect();
    assert_eq!(rounds, ["r2", "r1"]);
}

#[test]
fn timestamps_are_utc_calendar_dates() {
    assert_eq!(utc_timestamp(0), "1970-01-01T00:00:00Z");
    assert_eq!(utc_timestamp(951_782_400), "2000-02-29T00:00:00Z");
    assert_eq!(utc_timestamp(1_790_000_000), "2026-09-21T14:13:20Z");
}

#[test]
fn marker_edits_add_new_types_and_refuse_dropping_absent_markers() {
    let mut clusters = some(&[0]);
    let mut labs = labels(&["CT1"]);
    let mut ds = decisions(
        r#"{"action": "markers_add", "label": "CT 5", "features": ["GENE1", "GENE2"], "rationale": "r", "decided_by": "user"}
{"action": "markers_drop", "label": "CT1", "features": ["GENE3"], "rationale": "r", "decided_by": "user"}
{"cluster": 0, "action": "label", "label": "CT5", "rationale": "r", "decided_by": "user"}"#,
    );
    apply(&mut ds, &mut clusters, &mut labs, 1, "r1", "T").unwrap();
    assert_eq!(labs, labels(&["CT5"]), "cluster decisions still apply");

    let mut markers = vec![
        ("GENE3".to_string(), "CT1".to_string()),
        ("GENE4".to_string(), "CT1".to_string()),
        ("GENE1".to_string(), "CT_5".to_string()),
    ];
    let (edited, h) = apply_markers(&ds, &mut markers, "r1").unwrap();
    assert!(edited);
    assert_eq!(
        markers,
        vec![
            ("GENE4".to_string(), "CT1".to_string()),
            ("GENE1".to_string(), "CT_5".to_string()),
            ("GENE2".to_string(), "CT 5".to_string()),
        ],
        "an existing pair is not added twice, whatever the type's spacing"
    );
    assert_eq!(
        h["CT_5"][0].features,
        Some(vec!["GENE1".into(), "GENE2".into()])
    );
    assert_eq!(h["CT1"][0].timestamp, "T");

    let mut again = decisions(
        r#"{"action": "markers_drop", "label": "CT1", "features": ["GENE3"], "rationale": "r", "decided_by": "user"}"#,
    );
    apply(
        &mut again,
        &mut some(&[0]),
        &mut labels(&["CT1"]),
        1,
        "r2",
        "T",
    )
    .unwrap();
    assert!(apply_markers(&again, &mut markers, "r2").is_err());
}

#[test]
fn marker_edits_name_a_type_and_features_not_clusters() {
    let run = |line: &str| {
        apply(
            &mut decisions(line),
            &mut some(&[0]),
            &mut labels(&["CT1"]),
            1,
            "r",
            "T",
        )
    };
    assert!(run(
        r#"{"action": "markers_add", "label": "CT1", "rationale": "r", "decided_by": "user"}"#
    )
    .is_err());
    assert!(run(r#"{"action": "markers_add", "features": ["GENE1"], "rationale": "r", "decided_by": "user"}"#).is_err());
    assert!(run(r#"{"cluster": 0, "action": "markers_add", "label": "CT1", "features": ["GENE1"], "rationale": "r", "decided_by": "user"}"#).is_err());
}
