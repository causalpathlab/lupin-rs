//! The consensus proposal on made-up evidence, with a stub model that
//! needs no download.

use super::*;
use crate::annotate::rounds::Term;
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

fn choice(label: &str, share: f32, q: f32, rank: usize, markers: &[&str]) -> ChoiceInput {
    ChoiceInput {
        evidence: TypeEvidence {
            label: label.into(),
            enrichment: Some(share),
            q: Some(q),
            ..TypeEvidence::default()
        },
        ranks: [Some(rank), None, None],
        markers: markers.iter().map(|m| (*m).to_string()).collect(),
    }
}

/// Cluster 3: CT_A leads CT_B and CD8_T on the enrichment.
fn input() -> ClusterInput {
    ClusterInput {
        cluster: 3,
        size: 120,
        label: Some("CT_A".into()),
        calls: vec![
            Call {
                label: "CT_A".into(),
                q: Some(0.001),
            },
            Call {
                label: "CT_B".into(),
                q: Some(0.04),
            },
        ],
        terms: vec![Term {
            source: "go".into(),
            term: "T cell activation".into(),
            effect: 0.5,
            q: Some(0.01),
            p: None,
            nes: None,
        }],
        cl: Some(ClPlacement {
            id: "CL:1".into(),
            name: "alpha cell".into(),
            abstained: false,
        }),
        tops: vec![
            (Method::Enrichment, "CT_A".into()),
            (Method::Susie, "CT_B".into()),
        ],
        choices: vec![
            choice("CT_A", 0.7, 0.001, 0, &["GENE1", "GENE2"]),
            choice("CT_B", 0.2, 0.04, 1, &["GENE3"]),
            choice("CD8_T", 0.1, 0.3, 2, &[]),
            ChoiceInput::unassigned(),
        ],
        rationale: Some("an earlier reason".into()),
    }
}

/// Every type and term of the run: CT_C and a second term are another
/// cluster's.
fn run() -> BTreeSet<String> {
    [
        "CT_A",
        "CT_B",
        "CD8_T",
        "CT_C",
        "T cell activation",
        "B cell receptor signaling",
    ]
    .iter()
    .map(|s| (*s).to_string())
    .collect()
}

fn propose_with(stub: &mut Stub) -> Proposal {
    propose(stub, "m", "r1", &input(), &run()).unwrap()
}

#[test]
fn the_label_prompt_states_the_evidence_with_its_numbers() {
    let p = label_prompt(&input());
    assert!(
        p.contains("CLUSTER K3, 120 cells, current label CT A."),
        "{p}"
    );
    assert!(
        p.contains("1 CT A — enrichment rank 1, share 0.70, q 1e-3; markers here: GENE1, GENE2\n"),
        "{p}"
    );
    assert!(p.contains("4 unassigned — no type fits\n"), "{p}");
    assert!(p.contains("CALLS: CT A (q 1e-3); CT B (q 0.04)"), "{p}");
    assert!(p.contains("TOPS: enrichment CT A; susie CT B"), "{p}");
    assert!(
        p.contains("TERMS: go T cell activation (effect 0.50, q 0.01)"),
        "{p}"
    );
    assert!(
        p.contains("ONTOLOGY: alpha cell\nEARLIER: an earlier reason"),
        "{p}"
    );
    assert!(p.ends_with("<think>\n\n</think>\n\nLABEL:"), "{p}");
    assert_eq!(prompt(&input()), format!("{p}\n{}", why_prompt(&input())));
}

#[test]
fn the_why_prompt_states_the_same_evidence_in_words_then_asks_for_the_label() {
    let p = why_prompt(&input());
    assert!(p.starts_with("<|im_start|>system\n"), "{p}");
    assert!(
        p.ends_with("<think>\n\n</think>\n\nLABEL:"),
        "thinking off: {p}"
    );
    assert!(p.contains("THE CLUSTER: its current label is CT A."), "{p}");
    assert!(
        p.contains(
            "- CT A: enrichment ranks it first with most of the evidence, strongly significant; \
             its markers here: GENE1, GENE2\n"
        ),
        "{p}"
    );
    assert!(
        p.contains(
            "- CD8 T: enrichment ranks it third with some of the evidence, not significant\n"
        ),
        "{p}"
    );
    assert!(p.contains("- unassigned: no type fits\n"), "{p}");
    assert!(
        p.contains("CALLS: CT A (strongly significant); CT B (significant)"),
        "{p}"
    );
    assert!(
        p.contains("EACH METHOD'S TOP: enrichment CT A; susie CT B"),
        "{p}"
    );
    assert!(
        p.contains("GENE SETS: T cell activation (up, significant)"),
        "{p}"
    );
    assert!(p.contains("ONTOLOGY: alpha cell\n"), "{p}");
    assert!(p.contains("EARLIER: an earlier reason"), "{p}");
    // The user's part has no numbers for the sentence to copy, but the names'.
    let user = &p[p.find("THE CLUSTER").unwrap()..];
    let rest = ["CD8", "GENE1", "GENE2", "GENE3"]
        .iter()
        .fold(user.to_string(), |u, n| u.replace(n, ""));
    assert!(!rest.chars().any(|c| c.is_ascii_digit()), "{user}");
}

#[test]
fn the_label_is_the_best_scored_choice_and_every_choice_has_its_probability() {
    let mut stub = Stub::new(
        &[("CT A", -1.0), ("CT B", -0.2), ("CD8 T", -3.0)],
        &[", by its marker GENE3."],
    );
    let p = propose_with(&mut stub);
    assert_eq!(p.label, "CT_B");
    assert_eq!(p.p.len(), 4, "the three types and unassigned");
    let sum: f32 = p.p.values().sum();
    assert!((sum - 1.0).abs() < 1e-5, "{sum}");
    assert!(p.p_label() > p.p["CT_A"]);
    // The rationale is asked for after the label taken.
    assert!(
        stub.asked[0].0.ends_with("LABEL: CT B\nWHY: CT B"),
        "{}",
        stub.asked[0].0
    );
    assert_eq!(p.check, Check::Ok);
    assert_eq!(p.rationale, "CT B, by its marker GENE3.");
    assert_eq!(p.cited, ["CT_B", "GENE3"]);
    assert_eq!((p.model.as_str(), p.revision.as_str()), ("m", "r1"));
    assert_eq!(p.key, key(&prompt(&input()), "m", "r1"));
}

/// The grounding test: a sentence may name only the types and terms of the
/// cluster's evidence, and no numbers but those in the names it cites.
#[test]
fn a_sentence_names_only_the_types_and_terms_of_its_evidence() {
    let c = input();
    let check = |t: &str| check_lm(t, &run(), &c.named(), &c.own());
    assert_eq!(check("CT A fits: GENE1 and GENE2 rise here."), Ok(()));
    assert_eq!(
        check("CT As, by their markers."),
        Ok(()),
        "a plural names it"
    );
    assert_eq!(
        check("CT Ab fits."),
        Err(CitationError::NamesNone),
        "another word does not"
    );
    assert_eq!(
        check("CT A, with T cell activation among its terms."),
        Ok(())
    );
    assert_eq!(
        check("CD8 T cells, as the ontology's alpha cell says."),
        Ok(())
    );
    assert_eq!(
        check("CT C fits better than CT A."),
        Err(CitationError::Foreign("CT_C".into()))
    );
    assert_eq!(
        check("CT A, as B cell receptor signaling shows."),
        Err(CitationError::Foreign("B cell receptor signaling".into()))
    );
    assert_eq!(
        check("CT A, with a share of 0.7."),
        Err(CitationError::Digits("0".into()))
    );
    assert_eq!(check("Its cells look odd."), Err(CitationError::NamesNone));
    assert_eq!(
        check("CT A with markers MS₄A₃."),
        Err(CitationError::Digits("₄".into())),
        "any digit"
    );
    assert_eq!(
        check("CT A has q = $-e$."),
        Err(CitationError::Markup("=".into()))
    );
    assert_eq!(check("  "), Err(CitationError::Empty));
}

#[test]
fn a_failing_sentence_is_tried_again_with_the_offender_banned() {
    let mut stub = Stub::new(
        &[("CT A", -0.1)],
        &[" or CT C fits better.", ", with GENE1 high."],
    );
    let p = propose_with(&mut stub);
    assert_eq!(stub.asked.len(), 2);
    assert!(stub.asked[0].1.numbers, "numbers are never written");
    assert_eq!(
        stub.asked[0].1.keep,
        ["CD8 T", "GENE1", "GENE2", "GENE3"],
        "but the names' digits"
    );
    assert_eq!(
        stub.asked[1].1.words,
        ["CT_C"],
        "the foreign type is banned"
    );
    assert_eq!(p.check, Check::Ok);
    assert_eq!(p.rationale, "CT A, with GENE1 high.");
}

#[test]
fn a_sentence_failing_twice_or_naming_nothing_gives_the_template() {
    for sentences in [
        &[" or CT C.", " or CT C again."][..],
        &[" with a share of 0.7.", " is odd 3."][..],
    ] {
        let mut stub = Stub::new(&[("CT A", -0.1)], sentences);
        let p = propose_with(&mut stub);
        assert_eq!(p.check, Check::Template);
        assert_eq!(stub.asked.len(), 2, "one retry");
        assert!(
            p.rationale.starts_with("K3 is annotated as CT_A."),
            "{}",
            p.rationale
        );
        assert!(p.rationale.contains("GENE1, GENE2"), "{}", p.rationale);
        assert!(p.cited.contains(&"CT_A".to_string()));
    }
}

#[test]
fn a_number_bans_every_digit_and_a_sentence_naming_nothing_is_not_retried() {
    let mut stub = Stub::new(
        &[("CT A", -0.1)],
        &[" with a share of 0.7.", " ranks first. It is also CT B."],
    );
    let p = propose_with(&mut stub);
    assert!(
        stub.asked[1].1.keep.is_empty(),
        "every digit banned: {:?}",
        stub.asked[1].1
    );
    assert_eq!(p.rationale, "CT A ranks first.", "one sentence");
    // Unassigned has no label to lead with: a sentence may name nothing.
    let mut stub = Stub::new(&[("unassigned", -0.1)], &["Odd cells."]);
    let p = propose_with(&mut stub);
    assert_eq!(p.label, UNASSIGNED_LABEL);
    assert_eq!(stub.asked.len(), 1, "no retry without an offender");
    assert_eq!(p.check, Check::Template);
}

#[test]
fn the_key_is_the_prompt_the_model_and_its_revision() {
    let a = prompt(&input());
    let mut other = input();
    other.label = Some("CT_B".into());
    assert_eq!(key(&a, "m", "r1"), key(&a, "m", "r1"));
    assert_ne!(key(&a, "m", "r1"), key(&prompt(&other), "m", "r1"));
    assert_ne!(key(&a, "m", "r1"), key(&a, "m", "r2"));
    assert_ne!(key(&a, "m", "r1"), key(&a, "n", "r1"));
}

#[test]
fn a_proposal_reads_back_from_its_line() {
    let mut stub = Stub::new(&[("CT A", -0.1)], &[" fits."]);
    let p = propose_with(&mut stub);
    let line = serde_json::to_string(&p).unwrap();
    assert!(line.contains(r#""check":"ok""#), "{line}");
    assert_eq!(serde_json::from_str::<Proposal>(&line).unwrap(), p);
}

fn spec() -> LmSpec {
    LmSpec {
        model: super::super::lm::LmModel::Qwen3_0_6b,
        local_dir: None,
        device: ComputeDevice::Cpu,
    }
}

fn round() -> Round {
    let mut b = input();
    b.cluster = 4;
    b.label = Some("CT_B".into());
    Round {
        clusters: vec![input(), b],
        entities: run(),
    }
}

#[test]
fn a_round_is_proposed_once_then_read_from_its_cache_which_only_grows() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_path(&dir.path().join("r").to_string_lossy());
    let loads = Cell::new(0);
    let load = || -> Result<Box<dyn Lm>> {
        loads.set(loads.get() + 1);
        Ok(Box::new(Stub::new(
            &[("CT A", -0.1)],
            &["CT A.", "CT A too."],
        )))
    };
    let first = propose_round(&round(), &spec(), &cache, load).unwrap();
    assert_eq!(loads.get(), 1);
    assert_eq!(first.len(), 2);
    let text = std::fs::read_to_string(&cache).unwrap();
    assert_eq!(text.lines().count(), 2);
    // Again: every key is cached, so no model is loaded.
    let again = propose_round(&round(), &spec(), &cache, || -> Result<Box<dyn Lm>> {
        panic!("the cache has them all")
    })
    .unwrap();
    assert_eq!(again, first);
    assert_eq!(
        std::fs::read_to_string(&cache).unwrap(),
        text,
        "never rewritten"
    );
    // One cluster's evidence changed: only it is asked again, appended.
    let mut changed = round();
    changed.clusters[1].calls.clear();
    let third = propose_round(&changed, &spec(), &cache, load).unwrap();
    assert_eq!(loads.get(), 2);
    assert_eq!(third[0], first[0]);
    let grown = std::fs::read_to_string(&cache).unwrap();
    assert!(grown.starts_with(&text) && grown.lines().count() == 3);
    // What is current for each cluster follows its evidence now.
    let now = current(&prompts(&changed), &cache);
    assert_eq!(now[&3], first[0]);
    assert_eq!(now[&4], third[1]);
    assert_eq!(current(&prompts(&round()), &cache)[&4], first[1]);
}

#[test]
fn a_failing_model_load_fails_the_round_only_when_a_cluster_needs_it() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_path(&dir.path().join("r").to_string_lossy());
    let err = propose_round(&round(), &spec(), &cache, || -> Result<Box<dyn Lm>> {
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
fn a_rounds_inputs_are_its_review_with_the_tuis_choices_and_their_markers_here() {
    let dir = tempfile::tempdir().unwrap();
    let file = fixture_round(dir.path());
    let loaded = crate::manifest::run::load(&file.to_string_lossy()).unwrap();
    let r = inputs(&loaded).unwrap();
    assert_eq!(r.clusters.len(), 3);
    let c = &r.clusters[0];
    assert_eq!(
        (c.cluster, c.size, c.label.as_deref()),
        (0, 10, Some("CT1"))
    );
    let labels: Vec<&str> = c.choices.iter().map(ChoiceInput::label).collect();
    assert_eq!(labels, ["CT1", "CT2", "CT3", UNASSIGNED_LABEL]);
    assert_eq!(c.choices[0].ranks[0], Some(0));
    assert_eq!(c.choices[0].evidence.enrichment, Some(0.7));
    assert_eq!(
        c.choices[0].markers,
        ["GENE0", "GENE1", "GENE2", "GENE3", "GENE4"],
        "its markers among the cluster's specific genes, most specific first"
    );
    assert!(
        c.choices[1].markers.is_empty(),
        "CT2's markers do not rise here"
    );
    assert_eq!(c.tops, [(Method::Enrichment, "CT1".to_string())]);
    assert_eq!(c.calls[0].label, "CT1");
    assert_eq!(r.clusters[1].rationale.as_deref(), Some("its markers rise"));
    assert!(c.rationale.is_none());
    for t in ["CT1", "CT2", "CT3"] {
        assert!(r.entities.contains(t), "{t} is a type of the run");
    }
}

#[test]
fn the_choices_offered_stop_at_nine_with_unassigned_last() {
    use crate::tui::round::Candidate;
    let candidates: Vec<Candidate> = (0..12)
        .map(|i| Candidate {
            evidence: TypeEvidence {
                label: format!("CT{i}"),
                ..TypeEvidence::default()
            },
            share: 0.1,
            ranks: [Some(i % 3), Some(i / 3 % 3), None],
        })
        .collect();
    let offered = offered(&candidates, &BTreeMap::new(), &[]);
    assert_eq!(offered.len(), CHOICES);
    assert_eq!(offered[CHOICES - 1].label(), UNASSIGNED_LABEL);
}
