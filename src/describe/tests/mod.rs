use super::*;

fn check(
    s: &str,
    allowed: &BTreeSet<String>,
    c: &ClusterEvidence,
) -> Result<String, CitationError> {
    citation_check(s, allowed, &cluster_entities(c), &own_entities(c))
}

fn sample_sig() -> ClusterEvidence {
    ClusterEvidence {
        id: "K0".into(),
        coarse_label: "CT_A".into(),
        best_label: "CT_A".into(),
        best_q: Some(0.01),
        coarse_q: Some(0.01),
        label_support: None,
        best_significant: true,
        second: None,
        markers: vec!["GENE1".into(), "GENE2".into(), "GENE3".into()],
        neighbour_genes: vec![],
        incidence_words: vec![],
    }
}

#[test]
fn citation_accepts_template_with_best_label() {
    let c = ClusterEvidence {
        id: "K0".into(),
        coarse_label: "unassigned".into(),
        best_label: "CT_B".into(),
        best_q: Some(0.2),
        coarse_q: None,
        label_support: None,
        best_significant: false,
        second: None,
        markers: vec!["GENE4".into(), "GENE5".into()],
        neighbour_genes: vec![],
        incidence_words: vec!["wordk".into()],
    };
    let allowed = allowed_entities(std::slice::from_ref(&c));
    let s = template_sentence(&c);
    assert!(check(&s, &allowed, &c).is_ok());
    assert!(s.contains("CT_B"));
    assert!(s.contains("no significant"));
    assert!(s.contains("GENE4"));
    assert!(s.contains("Closest-type markers"));
    assert!(!s.contains("second significant"));
}

#[test]
fn significant_call_names_markers_and_q() {
    let c = sample_sig();
    let allowed = allowed_entities(std::slice::from_ref(&c));
    let s = template_sentence(&c);
    assert!(check(&s, &allowed, &c).is_ok());
    assert!(s.contains("annotated as CT_A"));
    assert!(s.contains("supported by markers GENE1, GENE2, GENE3"));
    assert!(s.contains("best q=0.010"));
}

#[test]
fn runner_up_only_when_significant() {
    let mut c = sample_sig();
    c.second = Some(("CT_B".into(), 0.04));
    let allowed = allowed_entities(std::slice::from_ref(&c));
    let s = template_sentence(&c);
    assert!(check(&s, &allowed, &c).is_ok());
    assert!(s.contains("second significant contender is CT_B (q=0.040)"));

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
    allowed.insert("CT_C".into());
    let bad = "K0 is annotated as CT_C.";
    assert!(check(bad, &allowed, &c).is_err());
}

#[test]
fn attach_markers_uses_best_label_when_nonsig() {
    let mut clusters = vec![ClusterEvidence {
        id: "K1".into(),
        coarse_label: "unassigned".into(),
        best_label: "CT_B".into(),
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
    by_type.insert("CT_B".into(), vec!["GENE5".into(), "GENE6".into()]);
    attach_markers(&mut clusters, &by_type);
    assert_eq!(clusters[0].markers, vec!["GENE5", "GENE6"]);
}

#[test]
fn incidence_words_aggregate_marker_edges() {
    let mut clusters = vec![ClusterEvidence {
        id: "K0".into(),
        coarse_label: "CT_A".into(),
        best_label: "CT_A".into(),
        best_q: Some(0.01),
        coarse_q: None,
        label_support: None,
        best_significant: true,
        second: None,
        markers: vec!["GENE1".into(), "GENE2".into()],
        neighbour_genes: vec!["GENE7".into()],
        incidence_words: vec![],
    }];
    let dir = tempfile_dir();
    let prefix = dir.join("text");
    let edges = format!("{}.feature_word.edges.tsv", prefix.display());
    std::fs::write(
        &edges,
        "gene\tGENE1\tword\tworda\t2.0\n\
         gene\tGENE1\tword\twordb\t1.5\n\
         gene\tGENE2\tword\tworda\t1.0\n\
         gene\tGENE7\tword\ttranscription\t3.0\n\
         gene\tOTHER\tword\tirrelevant\t9.0\n",
    )
    .unwrap();
    attach_incidence_words(&mut clusters, Some(prefix.to_str().unwrap())).unwrap();
    // worda 3.0, transcription 3.0, wordb 1.5 — OTHER ignored
    assert!(clusters[0].incidence_words.contains(&"worda".into()));
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
    let cols: Vec<Box<str>> = vec!["CT_A".into(), "T_cells".into()];
    q.to_parquet_with_names(&path, (Some(&rows), Some("cluster")), Some(&cols))
        .unwrap();

    let files = EvidenceFiles::locate(&prefix);
    assert_eq!(files.enrichment_q.as_deref(), Some(path.as_str()));
    let mut ev = load_evidence(&files, 0.1).unwrap();
    assert_eq!(ev[0].coarse_label, "CT_A");
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
    use crate::manifest::run::{RunKind, RunManifest};
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

#[test]
fn an_output_prefix_finds_its_labels_as_parquet_or_as_an_older_tsv() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("run").to_string_lossy().into_owned();
    let tsv = format!("{prefix}.argmax.tsv");
    std::fs::write(&tsv, "cell\tcell_type\nc0\tCT_A\n").unwrap();
    assert_eq!(
        EvidenceFiles::locate(&prefix).argmax.as_deref(),
        Some(tsv.as_str())
    );
    let ev = load_evidence(&EvidenceFiles::locate(&prefix), 0.1).unwrap();
    assert_eq!(ev[0].coarse_label, "CT_A");

    let cells: Vec<Box<str>> = vec!["c0".into()];
    graph_embedding_util::type_annotation::write_cell_labels(
        &prefix,
        &cells,
        &["CT_B".into()],
        &[0.9],
    )
    .unwrap();
    let parquet = format!("{prefix}.argmax.parquet");
    let files = EvidenceFiles::locate(&prefix);
    assert_eq!(files.argmax.as_deref(), Some(parquet.as_str()));
    assert_eq!(load_evidence(&files, 0.1).unwrap()[0].coarse_label, "CT_B");
}

/// `describe`'s arguments as the command line gives them.
fn describe_args(argv: &[&str]) -> DescribeArgs {
    use clap::Parser;
    #[derive(Parser)]
    struct Cli {
        #[command(flatten)]
        args: DescribeArgs,
    }
    Cli::try_parse_from(std::iter::once("describe").chain(argv.iter().copied()))
        .unwrap()
        .args
}

#[test]
fn the_flags_only_o_acts_on_are_named_without_it() {
    assert!(ignored_without_out(&describe_args(&["-f", "x"])).is_empty());
    let args = describe_args(&["-f", "x", "--markers", "m", "--lm", "--cache-only"]);
    assert_eq!(ignored_without_out(&args), ["--markers", "--cache-only"]);
}

#[test]
fn a_round_made_from_a_trajectory_is_described_as_a_round() {
    let dir = tempfile::tempdir().unwrap();
    let t = super::trajectory::tests::fixture_trajectory(dir.path());
    let loaded = crate::manifest::run::load(&t.to_string_lossy()).unwrap();
    let round = dir.path().join("t.L1.senna.json");
    let copy = loaded.copy_to(round.clone()).unwrap();
    copy.manifest.save(&copy.file).unwrap();
    let out = dir.path().join("t.L1").to_string_lossy().into_owned();
    let from = round.to_string_lossy().into_owned();
    // A round with no evidence of its own: nothing to describe.
    assert!(run_describe(&describe_args(&["-f", &from, "-o", &out])).is_err());
    assert!(
        !super::trajectory::md_path(&out).exists(),
        "not the trajectory it was made from"
    );
}

#[test]
fn cache_only_writes_the_models_texts_and_leaves_the_description() {
    use super::trajectory::{self, tests::fixture_trajectory, OrderFacts};
    use crate::describe::consensus::tests::Stub;
    let dir = tempfile::tempdir().unwrap();
    let t = fixture_trajectory(dir.path());
    let out = dir.path().join("t").to_string_lossy().into_owned();
    let md = trajectory::md_path(&out);
    std::fs::write(&md, "mine\n").unwrap();
    // The texts already cached, so no model is loaded.
    let loaded = crate::manifest::run::load(&t.to_string_lossy()).unwrap();
    let f = OrderFacts::read(&loaded).unwrap();
    let words: Vec<String> = std::iter::once("The order runs from HSC like.".to_string())
        .chain(f.types().iter().map(|t| format!("{t} sits here.")))
        .collect();
    let args = describe_args(&[
        "-f",
        &t.to_string_lossy(),
        "-o",
        &out,
        "--lm",
        "--lm-model",
        "qwen3-0.6b",
        "--device",
        "cpu",
        "--cache-only",
    ]);
    let cache = consensus::cache_path(&out);
    trajectory::describe(&f, Some(&args.lm.spec()), &cache, || {
        let refs: Vec<&str> = words.iter().map(String::as_str).collect();
        Ok(Box::new(Stub::new(&[], &refs)) as Box<dyn lm::Lm>)
    })
    .unwrap();
    run_describe(&args).unwrap();
    assert_eq!(std::fs::read_to_string(&md).unwrap(), "mine\n");
    assert!(cache.is_file());
}
