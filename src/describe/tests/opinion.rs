//! The language model's opinion on made-up evidence, with a stub model
//! that needs no download.

use super::*;
use crate::annotate::rounds::TypeEvidence;
use crate::describe::lm::Bans;
use crate::gene_text::cli::ComputeDevice;
use std::cell::Cell;
use std::collections::VecDeque;

/// A model that scores each choice from a table and writes the sentences
/// it is given in turn, recording what it was asked.
#[derive(Default)]
pub(crate) struct Stub {
    pub(crate) scores: HashMap<String, f32>,
    pub(crate) sentences: VecDeque<String>,
    /// Per sentence asked for: the prompt and what was banned.
    pub(crate) asked: Vec<(String, Bans)>,
}

impl Stub {
    pub(crate) fn new(scores: &[(&str, f32)], sentences: &[&str]) -> Self {
        Self {
            scores: scores.iter().map(|(c, s)| ((*c).to_string(), *s)).collect(),
            sentences: sentences.iter().map(|s| (*s).to_string()).collect(),
            asked: Vec::new(),
        }
    }
}

impl Lm for Stub {
    fn score(&mut self, _prompt: &str, continuations: &[String]) -> Result<Vec<f32>> {
        Ok(continuations
            .iter()
            .map(|c| self.scores.get(c.trim()).copied().unwrap_or(-9.0))
            .collect())
    }

    fn generate(&mut self, prompt: &str, _max: usize, bans: &Bans) -> Result<String> {
        self.asked.push((prompt.to_string(), bans.clone()));
        Ok(self.sentences.pop_front().unwrap_or_default())
    }
}

fn term(t: &str, effect: f32) -> Term {
    Term {
        source: "go".into(),
        term: t.into(),
        effect,
        q: Some(0.01),
        p: None,
        nes: None,
    }
}

/// A panel of CT_A, CT_B and CD8_T.
fn panel() -> Vec<PanelType> {
    [
        ("CT_A", &["GENE1", "GENE2"][..]),
        ("CT_B", &["GENE3"][..]),
        ("CD8_T", &["CD8A"][..]),
    ]
    .iter()
    .map(|(l, m)| PanelType {
        label: (*l).into(),
        markers: m.iter().map(|g| (*g).to_string()).collect(),
    })
    .collect()
}

/// Cluster 3: GENE1, GENE2 and CD8A among its top genes; CT_A's markers
/// GENE1 and GENE2 here.
fn input() -> ClusterInput {
    ClusterInput {
        cluster: 3,
        size: 120,
        genes: ["GENE1", "GENE2", "CD8A", "GENE9"]
            .iter()
            .map(|g| (*g).to_string())
            .collect(),
        here: BTreeMap::from([
            (
                "CT_A".to_string(),
                vec!["GENE1".to_string(), "GENE2".to_string()],
            ),
            ("CD8_T".to_string(), vec!["CD8A".to_string()]),
        ]),
        terms: vec![
            term("T cell activation", 0.5),
            term("B cell receptor signaling", -0.2),
        ],
    }
}

/// Every type and term of the run: CT_C and a third term are another
/// cluster's.
fn round() -> Round {
    let mut entities: BTreeSet<String> =
        ["CT_A", "CT_B", "CD8_T", "CT_C", "myeloid cell migration"]
            .iter()
            .map(|s| (*s).to_string())
            .collect();
    entities.extend(input().terms.iter().map(|t| t.term.clone()));
    let mut b = input();
    b.cluster = 4;
    b.genes.reverse();
    Round {
        clusters: vec![input(), b],
        panel: panel(),
        entities,
    }
}

fn opine_with(stub: &mut Stub) -> Opinion {
    opine(stub, "m", "r1", &input(), &round()).unwrap()
}

#[test]
fn the_opinion_prompt_states_the_panel_and_the_clusters_genes_and_terms_alone() {
    let p = opinion_prompt(&input(), &panel());
    assert!(p.starts_with("<|im_start|>system\n"), "{p}");
    assert!(
        p.ends_with("<think>\n\n</think>\n\nLABEL:"),
        "thinking off: {p}"
    );
    assert!(
        p.contains("TYPES:\n- CT A: GENE1, GENE2\n- CT B: GENE3\n- CD8 T: CD8A\n- unassigned: no type fits\n"),
        "{p}"
    );
    assert!(p.contains("CLUSTER: 120 cells.\n"), "{p}");
    assert!(p.contains("TOP GENES: GENE1, GENE2, CD8A, GENE9\n"), "{p}");
    assert!(
        p.contains("GENE SETS: T cell activation (up); B cell receptor signaling (down)"),
        "{p}"
    );
    // Nothing of another method's, the round's label or the ontology.
    for word in [
        "share",
        "q ",
        "pip",
        "rank",
        "enrichment",
        "susie",
        "score",
        "CALLS",
        "ONTOLOGY",
        "label ",
    ] {
        assert!(
            !p.to_lowercase().contains(&word.to_lowercase()),
            "{word}: {p}"
        );
    }
}

#[test]
fn the_biology_prompt_states_the_type_its_markers_here_the_genes_and_the_terms() {
    let p = biology_prompt(&input(), "CT_A");
    assert!(
        p.contains("CELL TYPE: CT A\nITS MARKERS HERE: GENE1, GENE2\nTOP GENES: GENE1"),
        "{p}"
    );
    assert!(p.contains("GENE SETS: T cell activation (up)"), "{p}");
    assert!(p.ends_with("<think>\n\n</think>\n\n"), "{p}");
    assert!(
        !biology_prompt(&input(), "CT_B").contains("MARKERS HERE"),
        "none of CT_B's markers are here"
    );
    for word in [
        "rank",
        "enrichment",
        "susie",
        "score",
        "before",
        "after",
        "first",
    ] {
        assert!(!p.to_lowercase().contains(word), "{word}: {p}");
    }
}

#[test]
fn the_label_is_the_best_scored_panel_type_and_every_type_has_its_probability() {
    let mut stub = Stub::new(
        &[("CT A", -1.0), ("CT B", -0.2), ("CD8 T", -3.0)],
        &[" express GENE1."],
    );
    let o = opine_with(&mut stub);
    assert_eq!(o.label, "CT_B");
    assert_eq!(o.p.len(), 4, "the three types and unassigned");
    let sum: f32 = o.p.values().sum();
    assert!((sum - 1.0).abs() < 1e-5, "{sum}");
    assert!(o.p_label() > o.p["CT_A"]);
    // The sentence is asked for after the type taken, its biology alone.
    let (asked, bans) = &stub.asked[0];
    assert!(
        asked.ends_with("<think>\n\n</think>\n\nCT B cells"),
        "{asked}"
    );
    assert!(asked.contains("CELL TYPE: CT B\n"), "{asked}");
    assert!(bans.numbers);
    for w in [
        "first",
        "before",
        "enrichment",
        "SuSiE",
        "supported",
        "root",
    ] {
        assert!(bans.words.contains(&w.to_string()), "{w}: {:?}", bans.words);
    }
    assert!(
        bans.keep.contains(&"CD8A".to_string()),
        "a name's digits: {:?}",
        bans.keep
    );
    assert_eq!(o.check, Check::Ok);
    assert_eq!(o.biology.as_deref(), Some("CT B cells express GENE1."));
    assert_eq!(o.cited, ["CT_B", "GENE1"]);
    assert_eq!((o.model.as_str(), o.revision.as_str()), ("m", "r1"));
    assert_eq!(o.key, key(&prompt(&input(), &panel()), "m", "r1"));
}

/// The grounding test: a sentence may name the type, its markers, the top
/// genes and the cluster's terms, no other type or term of the run, and no
/// number, sign or word that states what code states.
#[test]
fn a_biology_sentence_names_only_its_inputs_and_states_no_fact_of_the_evidence() {
    let r = round();
    let named = BTreeSet::from(["CT_A".to_string()]);
    let mut own: BTreeSet<String> = input().genes.into_iter().collect();
    own.extend(input().terms.iter().map(|t| t.term.clone()));
    let check = |t: &str| check(t, &r.entities, &named, &own);
    assert_eq!(check("CT A cells express GENE1 and CD8A."), Ok(()));
    assert_eq!(
        check("CT As, with T cell activation up."),
        Ok(()),
        "a plural names it"
    );
    assert_eq!(check("Its cells look odd."), Err(CitationError::NamesNone));
    assert_eq!(
        check("CT A, unlike CT C, expresses GENE1."),
        Err(CitationError::Foreign("CT_C".into()))
    );
    assert_eq!(
        check("CT A, with myeloid cell migration."),
        Err(CitationError::Foreign("myeloid cell migration".into()))
    );
    assert_eq!(
        check("CT A with GENE1 at 7."),
        Err(CitationError::Digits("7".into()))
    );
    assert_eq!(
        check("CT A has q = $-e$."),
        Err(CitationError::Markup("=".into()))
    );
    for (text, word) in [
        ("CT A expresses GENE1 first.", "first"),
        ("CT A gives rise to its neighbours.", "gives rise to"),
        ("CT A, as the enrichment shows.", "enrichment"),
        ("CT A sits at the root.", "root"),
        ("CT A has a supported edge.", "supported"),
    ] {
        assert_eq!(
            check(text),
            Err(CitationError::Vocabulary(word.into())),
            "{text}"
        );
    }
    assert_eq!(
        check(
            "CT A cells: GENE1, GENE2, CD8A, GENE9, T cell activation, B cell receptor signaling."
        ),
        Err(CitationError::List(6)),
        "a list, not a sentence"
    );
    assert_eq!(check("  "), Err(CitationError::Empty));
}

#[test]
fn a_biology_sentence_begins_with_the_types_cells() {
    assert_eq!(lead("CT_A"), "CT A cells");
    assert_eq!(lead("Plasma_Cell"), "Plasma Cell");
    assert_eq!(lead("NK_cells"), "NK cells");
}

#[test]
fn a_failing_sentence_is_tried_again_with_the_offender_banned_then_dropped() {
    let mut stub = Stub::new(
        &[("CT A", -0.1)],
        &[" ranks it by GENE1.", " express GENE1."],
    );
    let o = opine_with(&mut stub);
    assert_eq!(stub.asked.len(), 2);
    assert!(stub.asked[1].1.words.contains(&"ranks".to_string()));
    assert_eq!(
        (o.check, o.biology.as_deref()),
        (Check::Ok, Some("CT A cells express GENE1."))
    );
    let mut stub = Stub::new(&[("CT A", -0.1)], &[" or CT C.", " or CT C again."]);
    let o = opine_with(&mut stub);
    assert_eq!(
        stub.asked[1].1.words.last().map(String::as_str),
        Some("CT_C")
    );
    assert_eq!((o.check, o.biology), (Check::Template, None));
    assert!(o.cited.is_empty());
    // A number bans every digit, the names' too.
    let mut stub = Stub::new(&[("CT A", -0.1)], &[" at 7.", " at 8."]);
    opine_with(&mut stub);
    assert!(stub.asked[1].1.keep.is_empty(), "{:?}", stub.asked[1].1);
}

#[test]
fn unassigned_has_no_sentence_and_only_the_first_is_kept() {
    let mut stub = Stub::new(&[("unassigned", -0.1)], &["Odd cells."]);
    let o = opine_with(&mut stub);
    assert_eq!(o.label, UNASSIGNED_LABEL);
    assert!(stub.asked.is_empty(), "nothing asked");
    assert_eq!((o.check, o.biology), (Check::Template, None));
    let mut stub = Stub::new(
        &[("CT A", -0.1)],
        &[" express GENE1. They also express GENE2."],
    );
    assert_eq!(
        opine_with(&mut stub).biology.as_deref(),
        Some("CT A cells express GENE1.")
    );
}

fn candidate(label: &str, ranks: [Option<usize>; 4], ev: TypeEvidence) -> Candidate {
    Candidate {
        evidence: TypeEvidence {
            label: label.into(),
            ..ev
        },
        share: 0.0,
        ranks,
    }
}

/// CT_A first by the enrichment and the score, SuSiE's second; CT_B
/// SuSiE's first, the enrichment's second.
fn candidates() -> Vec<Candidate> {
    vec![
        candidate(
            "CT_A",
            [Some(0), Some(1), Some(0), None],
            TypeEvidence {
                enrichment: Some(0.7),
                q: Some(0.001),
                pip: Some(0.4),
                fold: Some(2.0),
                score: Some(0.55),
                ..TypeEvidence::default()
            },
        ),
        candidate(
            "CT_B",
            [Some(1), Some(0), None, Some(0)],
            TypeEvidence {
                enrichment: Some(0.2),
                q: Some(0.04),
                pip: Some(0.98),
                fold: Some(3.08),
                llm: Some(0.9),
                ..TypeEvidence::default()
            },
        ),
    ]
}

#[test]
fn the_facts_say_where_each_other_method_puts_the_type_with_its_readings() {
    let c = candidates();
    assert_eq!(
        facts("CT_B", &c),
        "CT_B: second by enrichment (share 0.20, q 0.04); first by SuSiE (PIP 0.98, ×3.1); \
         not ranked by the score."
    );
    assert_eq!(
        facts("CT_A", &c),
        "CT_A: first by enrichment (share 0.70, q 1e-3); second by SuSiE (PIP 0.40, ×2.0); \
         first by the score (share 0.55)."
    );
    assert_eq!(
        facts("CD8_T", &c),
        "CD8_T: not ranked by enrichment; not ranked by SuSiE; not ranked by the score."
    );
    assert_eq!(
        facts(UNASSIGNED_LABEL, &c),
        "unassigned; first are CT_A by enrichment, CT_B by SuSiE, CT_A by the score."
    );
    assert_eq!(
        facts("CT_A", &[]),
        "No other method ranks this cluster's types."
    );
    let mut stub = Stub::new(&[("CT B", -0.1)], &[" express GENE1."]);
    let o = opine_with(&mut stub);
    assert_eq!(
        o.text(&c),
        format!("{} CT B cells express GENE1.", facts("CT_B", &c))
    );
}

#[test]
fn the_key_is_the_prompts_the_model_and_its_revision() {
    let a = prompt(&input(), &panel());
    assert!(a.contains(&opinion_prompt(&input(), &panel())));
    assert!(a.contains(&biology_prompt(&input(), "CD8_T")));
    let mut other = input();
    other.genes.pop();
    assert_eq!(key(&a, "m", "r1"), key(&a, "m", "r1"));
    assert_ne!(
        key(&a, "m", "r1"),
        key(&prompt(&other, &panel()), "m", "r1")
    );
    assert_ne!(
        key(&a, "m", "r1"),
        key(&prompt(&input(), &panel()[1..]), "m", "r1")
    );
    assert_ne!(key(&a, "m", "r1"), key(&a, "m", "r2"));
    assert_ne!(key(&a, "m", "r1"), key(&a, "n", "r1"));
}

#[test]
fn an_opinion_reads_back_from_its_line_and_tables_what_the_model_puts_forward() {
    let mut stub = Stub::new(&[("CT A", -0.1), ("CT B", -1.0)], &[" express GENE1."]);
    let o = opine_with(&mut stub);
    let line = serde_json::to_string(&o).unwrap();
    assert!(line.contains(r#""check":"ok""#), "{line}");
    assert_eq!(serde_json::from_str::<Opinion>(&line).unwrap(), o);
    let t = table(&BTreeMap::from([(3, o.clone())]));
    assert_eq!(t.rows, [3]);
    assert_eq!(t.cols, ["CD8_T", "CT_A", "CT_B"], "unassigned left out");
    let at = |l: &str| t.values[t.cols.iter().position(|c| c == l).unwrap()];
    assert_eq!(at("CT_A"), o.p["CT_A"]);
    assert!(
        o.p["CD8_T"] < MIN_LLM_P && at("CD8_T").is_nan(),
        "below the floor: none"
    );
}

fn spec() -> LmSpec {
    LmSpec {
        model: super::super::lm::LmModel::Qwen3_0_6b,
        local_dir: None,
        device: ComputeDevice::Cpu,
    }
}

#[test]
fn a_round_is_opined_once_then_read_from_its_cache_which_only_grows() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_path(&dir.path().join("r").to_string_lossy());
    let loads = Cell::new(0);
    let load = || -> Result<Box<dyn Lm>> {
        loads.set(loads.get() + 1);
        Ok(Box::new(Stub::new(
            &[("CT A", -0.1)],
            &[" cells.", " cells too."],
        )))
    };
    let first = opine_round(&round(), &spec(), &cache, load).unwrap();
    assert_eq!(loads.get(), 1);
    assert_eq!(first.len(), 2);
    let text = std::fs::read_to_string(&cache).unwrap();
    assert_eq!(text.lines().count(), 2);
    // Again: every key is cached, so no model is loaded.
    let again = opine_round(&round(), &spec(), &cache, || -> Result<Box<dyn Lm>> {
        panic!("the cache has them all")
    })
    .unwrap();
    assert_eq!(again, first);
    assert_eq!(
        std::fs::read_to_string(&cache).unwrap(),
        text,
        "never rewritten"
    );
    // One cluster's genes changed: only it is asked again, appended.
    let mut changed = round();
    changed.clusters[1].genes.pop();
    let third = opine_round(&changed, &spec(), &cache, load).unwrap();
    assert_eq!(loads.get(), 2);
    assert_eq!(third[0], first[0]);
    let grown = std::fs::read_to_string(&cache).unwrap();
    assert!(grown.starts_with(&text) && grown.lines().count() == 3);
    // What is current for each cluster follows its inputs now.
    let now = current(&prompts(&changed), &cache);
    assert_eq!(now[&3], first[0]);
    assert_eq!(now[&4], third[1]);
    assert_eq!(current(&prompts(&round()), &cache)[&4], first[1]);
}

#[test]
fn a_failing_model_load_fails_the_round_only_when_a_cluster_needs_it() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_path(&dir.path().join("r").to_string_lossy());
    let err = opine_round(&round(), &spec(), &cache, || -> Result<Box<dyn Lm>> {
        anyhow::bail!(super::super::lm::OFFLINE)
    })
    .unwrap_err();
    assert!(format!("{err:#}").contains("run once online or give --lm-dir"));
}

/// A round of three clusters of ten cells in `dir`: CT1–CT3 each lead their
/// own cluster on the enrichment's share, their five markers rising there;
/// cluster 1 has a decision in its history.
pub(crate) fn fixture_round(dir: &Path) -> PathBuf {
    use crate::manifest::rounds::write_clusters;
    use crate::manifest::run::RunManifest;
    use legume_numeric::matrix::dense_mat_io::Mat;
    use legume_numeric::matrix::traits::IoOps;
    let at = |f: &str| dir.join(f).to_string_lossy().into_owned();
    let (g, k, per) = (30usize, 3usize, 10usize);
    let genes: Vec<Box<str>> = (0..g).map(|i| format!("GENE{i}").into()).collect();
    let cells: Vec<Box<str>> = (0..k * per).map(|i| format!("c{i}").into()).collect();
    let ids: Vec<Option<ClusterId>> = (0..k * per).map(|i| Some((i / per) as ClusterId)).collect();
    write_clusters(&at("run.clusters.parquet"), &cells, &ids).unwrap();
    let labels: Vec<Box<str>> = (0..k * per)
        .map(|i| format!("CT{}", i / per + 1).into())
        .collect();
    graph_embedding_util::type_annotation::write_cell_labels(
        &at("run"),
        &cells,
        &labels,
        &vec![0.9; k * per],
    )
    .unwrap();
    let mut expr = Mat::zeros(g, k);
    for i in 0..g {
        for c in 0..k {
            expr[(i, c)] = 1.0;
        }
    }
    for c in 0..k {
        for m in 0..5 {
            expr[(c * 5 + m, c)] = 20.0 - m as f32;
        }
    }
    let kcols: Vec<Box<str>> = (0..k).map(|c| format!("K{c}").into()).collect();
    expr.to_parquet_with_names(
        &at("run.cluster_expression.parquet"),
        (Some(&genes), Some("gene")),
        Some(&kcols),
    )
    .unwrap();
    let mut panel = String::from("gene\tcelltype\n");
    for c in 0..k {
        for m in 0..5 {
            panel += &format!("GENE{}\tCT{}\n", c * 5 + m, c + 1);
        }
    }
    std::fs::write(dir.join("markers.tsv"), panel).unwrap();
    let types: Vec<Box<str>> = (0..k).map(|c| format!("CT{}", c + 1).into()).collect();
    let table = |own: f32, other: f32, file: &str| {
        let mut m = Mat::zeros(k, k);
        for r in 0..k {
            for c in 0..k {
                m[(r, c)] = if r == c { own } else { other };
            }
        }
        m.to_parquet_with_names(&at(file), (Some(&kcols), Some("cluster")), Some(&types))
            .unwrap();
    };
    table(0.7, 0.15, "run.cluster_celltype_q.parquet");
    table(0.001, 0.5, "run.cluster_celltype_q_values.parquet");
    std::fs::write(
        dir.join("run.annotation_history.json"),
        r#"{"1": [{"round": "run.senna.json", "action": "label", "label": "CT2",
            "rationale": "its markers rise", "decided_by": "user",
            "timestamp": "2026-01-01T00:00:00Z"}]}"#,
    )
    .unwrap();
    let mut m = RunManifest::new(crate::manifest::run::RunKind::Topic, "run");
    m.cluster.clusters = Some("run.clusters.parquet".into());
    let a = &mut m.annotate;
    a.argmax = Some("run.argmax.parquet".into());
    a.markers = Some("markers.tsv".into());
    a.cluster_expression = Some("run.cluster_expression.parquet".into());
    a.cluster_celltype_q = Some("run.cluster_celltype_q.parquet".into());
    a.cluster_celltype_q_values = Some("run.cluster_celltype_q_values.parquet".into());
    a.history = Some("run.annotation_history.json".into());
    a.called_by = Some(Method::Enrichment);
    let file = dir.join("run.senna.json");
    m.save(&file).unwrap();
    file
}

#[test]
fn a_rounds_inputs_are_its_sizes_specific_genes_and_panel() {
    let dir = tempfile::tempdir().unwrap();
    let file = fixture_round(dir.path());
    let loaded = crate::manifest::run::load(&file.to_string_lossy()).unwrap();
    let r = inputs(&loaded).unwrap();
    assert_eq!(r.clusters.len(), 3);
    let labels: Vec<&str> = r.panel.iter().map(|t| t.label.as_str()).collect();
    assert_eq!(labels, ["CT1", "CT2", "CT3"], "in the panel's order");
    assert_eq!(
        r.panel[1].markers,
        ["GENE5", "GENE6", "GENE7", "GENE8", "GENE9"]
    );
    let c = &r.clusters[0];
    assert_eq!((c.cluster, c.size), (0, 10));
    assert_eq!(
        c.genes,
        ["GENE0", "GENE1", "GENE2", "GENE3", "GENE4"],
        "its specific genes, most specific first"
    );
    assert_eq!(
        c.here.keys().collect::<Vec<_>>(),
        ["CT1"],
        "only CT1's markers rise here"
    );
    for t in ["CT1", "CT2", "CT3"] {
        assert!(r.entities.contains(t), "{t} is a type of the run");
    }
    assert_eq!(cache_of(&loaded), dir.path().join("run.llm.jsonl"));
}

/// The model's opinions, cached beside the round, are evidence its review
/// reads; a round without them reviews as before.
#[test]
fn a_rounds_review_has_the_models_top_and_readings_only_where_it_ran() {
    let dir = tempfile::tempdir().unwrap();
    let file = fixture_round(dir.path());
    let loaded = crate::manifest::run::load(&file.to_string_lossy()).unwrap();
    let before = crate::manifest::rounds::review_json(&loaded, &[]).unwrap();
    assert!(before.values().all(|v| v["tops"].get("llm").is_none()));
    assert!(!serde_json::to_string(&before).unwrap().contains("\"llm\""));
    let r = inputs(&loaded).unwrap();
    opine_round(
        &r,
        &spec(),
        &cache_of(&loaded),
        || -> Result<Box<dyn Lm>> { Ok(Box::new(Stub::new(&[("CT3", -0.1), ("CT2", -1.5)], &[]))) },
    )
    .unwrap();
    let after = crate::manifest::rounds::review_json(&loaded, &[]).unwrap();
    for (id, v) in &after {
        assert_eq!(v["tops"]["llm"], "CT3", "K{id}");
        assert_eq!(
            v["digest"], before[id]["digest"],
            "its calls stay the round's"
        );
        let types = v["types"].as_array().unwrap();
        let ct3 = types.iter().find(|t| t["celltype"] == "CT3").unwrap();
        assert!(ct3["llm"].as_f64().unwrap() > 0.5);
    }
}
