use super::*;

fn sample_sig() -> ClusterEvidence {
    ClusterEvidence {
        id: "K0".into(),
        coarse_label: "B_cells".into(),
        best_label: "B_cells".into(),
        best_q: Some(0.01),
        coarse_q: Some(0.01),
        label_support: None,
        best_significant: true,
        second: None,
        markers: vec!["CD19".into(), "MS4A1".into(), "CD79A".into()],
        neighbour_genes: vec![],
        incidence_words: vec![],
    }
}

#[test]
fn citation_accepts_template_with_best_label() {
    let c = ClusterEvidence {
        id: "K0".into(),
        coarse_label: "unassigned".into(),
        best_label: "NK".into(),
        best_q: Some(0.2),
        coarse_q: None,
        label_support: None,
        best_significant: false,
        second: None,
        markers: vec!["NCAM1".into(), "NKG7".into()],
        neighbour_genes: vec![],
        incidence_words: vec!["killer".into()],
    };
    let allowed = allowed_entities(std::slice::from_ref(&c));
    let s = template_sentence(&c);
    assert!(citation_check(&s, &allowed, &c).is_ok());
    assert!(s.contains("NK"));
    assert!(s.contains("no significant"));
    assert!(s.contains("NCAM1"));
    assert!(s.contains("Closest-type markers"));
    assert!(!s.contains("second significant"));
}

#[test]
fn significant_call_names_markers_and_q() {
    let c = sample_sig();
    let allowed = allowed_entities(std::slice::from_ref(&c));
    let s = template_sentence(&c);
    assert!(citation_check(&s, &allowed, &c).is_ok());
    assert!(s.contains("annotated as B_cells"));
    assert!(s.contains("supported by markers CD19, MS4A1, CD79A"));
    assert!(s.contains("best q=0.010"));
}

#[test]
fn runner_up_only_when_significant() {
    let mut c = sample_sig();
    c.second = Some(("NK".into(), 0.04));
    let allowed = allowed_entities(std::slice::from_ref(&c));
    let s = template_sentence(&c);
    assert!(citation_check(&s, &allowed, &c).is_ok());
    assert!(s.contains("second significant contender is NK (q=0.040)"));

    c.best_significant = false;
    c.coarse_label = "unassigned".into();
    c.second = None;
    let s2 = template_sentence(&c);
    assert!(!s2.contains("second significant"));
}

#[test]
fn citation_rejects_foreign_panel_type() {
    let c = sample_sig();
    let mut allowed = allowed_entities(std::slice::from_ref(&c));
    allowed.insert("Tcell".into());
    let bad = "K0 is annotated as Tcell.";
    assert!(citation_check(bad, &allowed, &c).is_err());
}

#[test]
fn attach_markers_uses_best_label_when_nonsig() {
    let mut clusters = vec![ClusterEvidence {
        id: "K1".into(),
        coarse_label: "unassigned".into(),
        best_label: "NK".into(),
        best_q: Some(0.4),
        coarse_q: None,
        label_support: None,
        best_significant: false,
        second: None,
        markers: vec![],
        neighbour_genes: vec![],
        incidence_words: vec![],
    }];
    let mut by_type = BTreeMap::new();
    by_type.insert("NK".into(), vec!["NKG7".into(), "GNLY".into()]);
    attach_markers(&mut clusters, &by_type);
    assert_eq!(clusters[0].markers, vec!["NKG7", "GNLY"]);
}

#[test]
fn incidence_words_aggregate_marker_edges() {
    let mut clusters = vec![ClusterEvidence {
        id: "K0".into(),
        coarse_label: "B_cells".into(),
        best_label: "B_cells".into(),
        best_q: Some(0.01),
        coarse_q: None,
        label_support: None,
        best_significant: true,
        second: None,
        markers: vec!["CD19".into(), "MS4A1".into()],
        neighbour_genes: vec!["PAX5".into()],
        incidence_words: vec![],
    }];
    let dir = tempfile_dir();
    let prefix = dir.join("text");
    let edges = format!("{}.feature_word.edges.tsv", prefix.display());
    std::fs::write(
        &edges,
        "gene\tCD19\tword\tlymphocyte\t2.0\n\
         gene\tCD19\tword\tbcell\t1.5\n\
         gene\tMS4A1\tword\tlymphocyte\t1.0\n\
         gene\tPAX5\tword\ttranscription\t3.0\n\
         gene\tOTHER\tword\tirrelevant\t9.0\n",
    )
    .unwrap();
    attach_incidence_words(&mut clusters, Some(prefix.to_str().unwrap())).unwrap();
    // lymphocyte 3.0, transcription 3.0, bcell 1.5 — OTHER ignored
    assert!(clusters[0].incidence_words.contains(&"lymphocyte".into()));
    assert!(clusters[0]
        .incidence_words
        .contains(&"transcription".into()));
    assert!(!clusters[0].incidence_words.contains(&"irrelevant".into()));
}

fn tempfile_dir() -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "lupin-describe-{}",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn enrichment_q_table_calls_each_cluster_by_its_lowest_q() {
    use legume_numeric::matrix::traits::IoOps;
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("run").to_string_lossy().into_owned();
    let path = format!("{prefix}{CLUSTER_CELLTYPE_Q_VALUES}");
    let q =
        legume_numeric::matrix::dense_mat_io::Mat::from_row_slice(2, 2, &[0.01, 0.05, 0.5, 0.3]);
    let rows: Vec<Box<str>> = vec!["K0".into(), "K1".into()];
    let cols: Vec<Box<str>> = vec!["B_cells".into(), "T_cells".into()];
    q.to_parquet_with_names(&path, (Some(&rows), Some("cluster")), Some(&cols))
        .unwrap();

    let files = EvidenceFiles::locate(&prefix);
    assert_eq!(files.enrichment_q.as_deref(), Some(path.as_str()));
    let mut ev = load_evidence(&files, 0.1).unwrap();
    assert_eq!(ev[0].coarse_label, "B_cells");
    assert!(ev[0].best_significant);
    assert_eq!(ev[1].coarse_label, UNASSIGNED_LABEL);
    assert_eq!(ev[1].best_label, "T_cells");

    attach_second_best(&mut ev, &files, 0.1).unwrap();
    assert_eq!(
        ev[0].second.as_ref().map(|(l, _)| l.as_str()),
        Some("T_cells")
    );
}

#[test]
fn a_manifest_points_describe_at_the_latest_pass_not_stale_files() {
    use crate::run_manifest::{RunKind, RunManifest};
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("run").to_string_lossy().into_owned();
    // A stale projection table from an earlier pass on the same prefix...
    std::fs::write(format!("{prefix}{ANNOT_PARQUET}"), b"stale").unwrap();
    // ...while the manifest records the latest (enrichment) pass.
    let mut m = RunManifest::new(RunKind::Topic, &prefix);
    m.annotate.argmax = Some("run.argmax.tsv".into());
    m.annotate.cluster_celltype_q_values = Some("run.cluster_celltype_q_values.parquet".into());
    m.save(std::path::Path::new(&format!("{prefix}.senna.json")))
        .unwrap();

    let files = EvidenceFiles::locate(&format!("{prefix}.senna.json"));
    assert!(
        files.annot.is_none(),
        "stale projection table must be ignored"
    );
    assert!(files
        .enrichment_q
        .as_deref()
        .is_some_and(|p| p.ends_with("run.cluster_celltype_q_values.parquet")));
    assert!(files.term_q.is_none());
}
