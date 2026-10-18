//! A proposed label for each cluster of a round, with one grounded sentence
//! on why: the consensus `lupin describe --lm` writes and the annotate TUI
//! offers among its label choices. Nothing is applied: a proposal is taken
//! only when the user confirms it, as any other label.
//!
//! A cluster's evidence is what `review --json` gives (its size, label,
//! top calls and terms, its Cell Ontology placement, each method's top
//! type) with the TUI's label choices ([`crate::tui::round::choices`]),
//! each with its readings, ranks and the panel markers among the cluster's
//! specific genes; `unassigned` is always a choice. The model
//! ([`super::lm`]) scores each choice as the continuation of `LABEL:`, so
//! the label is one of the choices by construction; then, from the same
//! evidence stated in words (a small model copies the numbers it is shown),
//! it continues `WHY: {label}` greedily, and the sentence passes the citation check and
//! the check of what it claims of the methods' ranks ([`super::claims`]) or is replaced by
//! the template ([`super::template_led`]).
//!
//! Prior art for language-model cell typing: Hou & Ji, *Assessing GPT-4 for
//! cell type annotation in single-cell RNA-seq analysis*, Nat Methods
//! 21:1462 (2024; GPTCelltype).

use super::lm::{Bans, Lm, LmSpec, MARKUP};
use super::{citation_check, mentions, template_led, CitationError, ClusterEvidence};
use crate::annotate::markers::label_key;
use crate::annotate::rounds::{Call, ClPlacement, ClusterId, Method, Term, TypeEvidence};
use crate::manifest::run::Loaded;
use anyhow::{Context, Result};
use enrichment::UNASSIGNED_LABEL;
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use sha2::{Digest as _, Sha256};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::Write;
use std::path::{Path, PathBuf};

/// Choices offered to the model, `unassigned` the last: as many as the
/// TUI keys `1`–`9`.
pub const CHOICES: usize = 9;
/// The longest rationale, in tokens.
pub const RATIONALE_TOKENS: usize = 60;
/// Calls, terms and markers per choice the prompt quotes.
const TOP: usize = 5;
/// How the model is asked; part of the cache key.
const DECODE: &str = "label: mean log-probability of each choice, softmax; why: from the label, \
                      greedy to a newline, 60 tokens, no numbers or markup but the names', \
                      one retry with the offender banned; claims checked against the evidence";

/// One label choice as the model is shown it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChoiceInput {
    pub evidence: TypeEvidence,
    /// Where each method ([`Method::ALL`]) ranks it; 0 is its top.
    pub ranks: [Option<usize>; 3],
    /// Its panel markers among the cluster's specific genes, most specific
    /// first.
    pub markers: Vec<String>,
}

impl ChoiceInput {
    #[must_use]
    pub fn label(&self) -> &str {
        &self.evidence.label
    }

    fn unassigned() -> Self {
        Self {
            evidence: TypeEvidence {
                label: UNASSIGNED_LABEL.into(),
                ..TypeEvidence::default()
            },
            ..Self::default()
        }
    }
}

/// One cluster's evidence, as the prompt states it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ClusterInput {
    pub cluster: ClusterId,
    pub size: usize,
    pub label: Option<String>,
    pub calls: Vec<Call>,
    pub terms: Vec<Term>,
    pub cl: Option<ClPlacement>,
    /// Each method's first type.
    pub tops: Vec<(Method, String)>,
    pub choices: Vec<ChoiceInput>,
    /// The rationale of the cluster's last decision.
    pub rationale: Option<String>,
}

impl ClusterInput {
    /// What a sentence about the cluster must name one of: its choices and
    /// its terms.
    #[must_use]
    pub fn named(&self) -> BTreeSet<String> {
        let choices = self.choices.iter().map(|c| c.label().to_string());
        choices
            .chain(self.terms.iter().map(|t| t.term.clone()))
            .collect()
    }

    /// What it may name besides: the choices' markers, the label, the
    /// calls and the ontology's placement.
    #[must_use]
    pub fn own(&self) -> BTreeSet<String> {
        let mut s: BTreeSet<String> = self
            .choices
            .iter()
            .flat_map(|c| c.markers.iter().cloned())
            .collect();
        s.extend(self.label.clone());
        s.extend(self.calls.iter().map(|c| c.label.clone()));
        s.extend(self.cl.as_ref().map(|c| c.name.clone()));
        s
    }
}

/// A round as the model reads it: its clusters, and every type and term
/// in it, which a cluster's sentence may name only when they are its own.
#[derive(Clone, Debug, Default)]
pub struct Round {
    pub clusters: Vec<ClusterInput>,
    pub entities: BTreeSet<String>,
}

/// The round `loaded`, read as `review --json` reads it, with its label
/// choices, their markers among each cluster's specific genes, and the
/// panel's types.
pub fn inputs(loaded: &Loaded) -> Result<Round> {
    use crate::tui::round::{panel_sets, ranked_view, specific_genes, Expression};
    use legume_numeric::matrix::dense_mat_io::Mat;
    use legume_numeric::matrix::traits::IoOps;
    let r = crate::manifest::rounds::Reviewed::read(loaded, true)?;
    let a = &loaded.manifest.annotate;
    let by = crate::manifest::rounds::called_by(a).unwrap_or_default();
    let panel = match a.markers.as_deref().and_then(|m| loaded.recorded(m)) {
        Some(p) => panel_sets(&p)?,
        None => BTreeMap::new(),
    };
    let specific = match a
        .cluster_expression
        .as_deref()
        .and_then(|e| loaded.recorded(e))
    {
        Some(p) => specific_genes(&Expression::new(
            Mat::from_parquet_with_row_names(&p, Some(0))
                .with_context(|| format!("reading {p}"))?,
        )),
        None => BTreeMap::new(),
    };
    let mut entities: BTreeSet<String> = panel.keys().cloned().collect();
    let mut clusters = Vec::new();
    for (&id, d) in &r.digests {
        let types = r.tables.evidence(id);
        entities.extend(types.iter().map(|t| t.label.clone()));
        entities.extend(d.terms.iter().map(|t| t.term.clone()));
        let tops = Method::ALL
            .into_iter()
            .filter_map(|m| Some((m, m.ranked(&types).first()?.label.clone())))
            .collect();
        let genes = specific.get(&id).map_or(&[][..], Vec::as_slice);
        let (candidates, _) = ranked_view(&types, by);
        let offered = offered(&candidates, &panel, genes);
        let cl = d.cl.clone().or_else(|| {
            r.placement(d).map(|p| ClPlacement {
                id: p.id,
                name: p.name,
                abstained: false,
            })
        });
        let rationale = r
            .history
            .get(&id)
            .and_then(|h| h.first())
            .map(|h| h.rationale.clone())
            .filter(|s| !s.is_empty());
        clusters.push(ClusterInput {
            cluster: id,
            size: d.size,
            label: d.label.clone(),
            calls: d.calls.iter().take(TOP).cloned().collect(),
            terms: d.terms.iter().take(TOP).cloned().collect(),
            cl,
            tops,
            choices: offered,
            rationale,
        });
    }
    Ok(Round { clusters, entities })
}

/// The label choices offered for a cluster of `candidates`: the TUI's
/// ([`crate::tui::round::choices`]), as many as there are keys for with
/// `unassigned` last, each with its `panel` markers among the cluster's
/// specific `genes`, most specific first.
pub(crate) fn offered(
    candidates: &[crate::tui::round::Candidate],
    panel: &BTreeMap<String, BTreeSet<String>>,
    genes: &[(String, f32)],
) -> Vec<ChoiceInput> {
    let mut out: Vec<ChoiceInput> = crate::tui::round::choices(candidates)
        .into_iter()
        .filter(|c| c.label != UNASSIGNED_LABEL)
        .take(CHOICES - 1)
        .map(|c| {
            let markers = panel.get(&label_key(&c.label));
            ChoiceInput {
                evidence: c.evidence.clone(),
                ranks: c.ranks,
                markers: genes
                    .iter()
                    .filter(|(g, _)| markers.is_some_and(|m| m.contains(g)))
                    .take(TOP)
                    .map(|(g, _)| g.clone())
                    .collect(),
            }
        })
        .collect();
    out.push(ChoiceInput::unassigned());
    out
}

/// A label as the prompt writes it: words, not the panel's underscores.
fn shown(label: &str) -> String {
    label.replace('_', " ")
}

/// A rank in words.
fn ordinal(rank: usize) -> String {
    const WORDS: [&str; 10] = [
        "first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth",
        "tenth",
    ];
    WORDS
        .get(rank)
        .map_or_else(|| "lower".into(), |w| (*w).to_string())
}

/// A share of the cluster's evidence (or cells), in words.
fn share_words(v: f32) -> &'static str {
    match v {
        v if v >= 0.5 => "most",
        v if v >= 0.2 => "much",
        v if v >= 0.05 => "some",
        _ => "little",
    }
}

/// An FDR q-value, in words.
fn q_words(q: f32) -> &'static str {
    match q {
        q if q < 0.01 => "strongly significant",
        q if q < 0.1 => "significant",
        _ => "not significant",
    }
}

/// SuSiE's posterior inclusion probability, in words.
fn pip_words(pip: f32) -> &'static str {
    match pip {
        p if p >= 0.95 => "certain",
        p if p >= 0.5 => "likely",
        _ => "unlikely",
    }
}

/// The fold a type's markers rise by, in words.
fn fold_words(fold: f32) -> &'static str {
    match fold {
        f if f >= 4.0 => "a strong rise",
        f if f >= 2.0 => "a clear rise",
        _ => "a weak rise",
    }
}

/// Qwen3's chat with thinking switched off: `system`, `user`, and the
/// assistant's answer begun with `answer`.
pub(crate) fn chat(system: &str, user: &str, answer: &str) -> String {
    format!(
        "<|im_start|>system\n{system}<|im_end|>\n<|im_start|>user\n{user}<|im_end|>\n\
         <|im_start|>assistant\n<think>\n\n</think>\n\n{answer}"
    )
}

const SYSTEM: &str = "You label one single-cell cluster. Choose only from CHOICES.\n\
                      Cite only names that appear below. One sentence in words, no numbers.\n\
                      Answer in two lines: LABEL: the choice, then WHY: the sentence.";

/// What the model is asked about cluster `c`, for the cache key: both its
/// prompts.
#[must_use]
pub fn prompt(c: &ClusterInput) -> String {
    format!("{}\n{}", label_prompt(c), why_prompt(c))
}

/// A number for the label prompt, in two significant digits.
fn num(v: f32) -> String {
    if v != 0.0 && v.abs() < 0.01 {
        format!("{v:.0e}")
    } else {
        format!("{v:.2}")
    }
}

/// The prompt the label is chosen after, up to `LABEL:`: cluster `c`'s
/// evidence with its numbers.
#[must_use]
pub fn label_prompt(c: &ClusterInput) -> String {
    let label = c.label.as_deref().map_or(UNASSIGNED_LABEL.into(), shown);
    let mut user = format!(
        "CLUSTER K{}, {} cells, current label {label}.\nCHOICES:\n",
        c.cluster, c.size
    );
    for (i, ch) in c.choices.iter().enumerate() {
        user += &format!("{} {}", i + 1, shown(ch.label()));
        if ch.label() == UNASSIGNED_LABEL {
            user += " — no type fits\n";
            continue;
        }
        let e = &ch.evidence;
        let mut parts = Vec::new();
        if let Some(r) = ch.ranks[Method::Enrichment.index()] {
            let mut s = format!("enrichment rank {}", r + 1);
            if let Some(v) = e.enrichment {
                s += &format!(", share {}", num(v));
            }
            if let Some(q) = e.q {
                s += &format!(", q {}", num(q));
            }
            parts.push(s);
        }
        if let Some(pip) = e.pip {
            let mut s = format!("susie pip {}", num(pip));
            if let Some(f) = e.fold {
                s += &format!(", fold {f:.1}");
            }
            parts.push(s);
        }
        if let Some(v) = e.score {
            parts.push(format!("score share {}", num(v)));
        }
        if !ch.markers.is_empty() {
            parts.push(format!("markers here: {}", ch.markers.join(", ")));
        }
        user += &format!(" — {}\n", parts.join("; "));
    }
    let calls: Vec<String> = (c.calls.iter())
        .map(|k| match k.q {
            Some(q) => format!("{} (q {})", shown(&k.label), num(q)),
            None => shown(&k.label),
        })
        .collect();
    let terms: Vec<String> = (c.terms.iter())
        .map(|t| {
            let q = t.q.map(|q| format!(", q {}", num(q))).unwrap_or_default();
            format!("{} {} (effect {}{q})", t.source, t.term, num(t.effect))
        })
        .collect();
    user += &tail(c, &calls, "TOPS", &terms, "TERMS");
    chat(SYSTEM, user.trim_end(), "LABEL:")
}

/// The lines after the choices: the calls, each method's top (under
/// `tops`), the gene sets (under `sets`), the ontology and the last reason.
fn tail(c: &ClusterInput, calls: &[String], tops: &str, terms: &[String], sets: &str) -> String {
    let mut s = String::new();
    if !calls.is_empty() {
        s += &format!("CALLS: {}\n", calls.join("; "));
    }
    if !c.tops.is_empty() {
        let t: Vec<String> = (c.tops.iter())
            .map(|(m, t)| format!("{} {}", m.name(), shown(t)))
            .collect();
        s += &format!("{tops}: {}\n", t.join("; "));
    }
    if !terms.is_empty() {
        s += &format!("{sets}: {}\n", terms.join("; "));
    }
    if let Some(cl) = &c.cl {
        let abst = if cl.abstained { " (abstained)" } else { "" };
        s += &format!("ONTOLOGY: {}{abst}\n", cl.name);
    }
    if let Some(r) = &c.rationale {
        s += &format!("EARLIER: {r}\n");
    }
    s
}

/// The prompt the sentence is written after, up to `LABEL:`: the same
/// evidence in words, so the sentence has no numbers to copy (the TUI shows
/// them).
#[must_use]
pub fn why_prompt(c: &ClusterInput) -> String {
    let label = c.label.as_deref().map_or(UNASSIGNED_LABEL.into(), shown);
    let mut user = format!("THE CLUSTER: its current label is {label}.\nCHOICES:\n");
    for ch in &c.choices {
        user += &format!("- {}", shown(ch.label()));
        if ch.label() == UNASSIGNED_LABEL {
            user += ": no type fits\n";
            continue;
        }
        let e = &ch.evidence;
        let mut parts = Vec::new();
        if let Some(r) = ch.ranks[Method::Enrichment.index()] {
            let mut s = format!("enrichment ranks it {}", ordinal(r));
            if let Some(v) = e.enrichment {
                s += &format!(" with {} of the evidence", share_words(v));
            }
            if let Some(q) = e.q {
                s += &format!(", {}", q_words(q));
            }
            parts.push(s);
        }
        if let Some(pip) = e.pip {
            let mut s = format!(
                "SuSiE ranks it {}",
                ch.ranks[Method::Susie.index()].map_or("lower".into(), ordinal)
            );
            s += &format!(", {}", pip_words(pip));
            if let Some(f) = e.fold {
                s += &format!(", its markers with {}", fold_words(f));
            }
            parts.push(s);
        }
        if let Some(v) = e.score {
            parts.push(format!(
                "the score calls {} of the cells it",
                share_words(v)
            ));
        }
        if !ch.markers.is_empty() {
            parts.push(format!("its markers here: {}", ch.markers.join(", ")));
        }
        user += &format!(": {}\n", parts.join("; "));
    }
    let calls: Vec<String> = (c.calls.iter())
        .map(|k| match k.q {
            Some(q) => format!("{} ({})", shown(&k.label), q_words(q)),
            None => shown(&k.label),
        })
        .collect();
    let terms: Vec<String> = (c.terms.iter())
        .map(|t| {
            let way = if t.effect >= 0.0 { "up" } else { "down" };
            let q = t.q.map(|q| format!(", {}", q_words(q))).unwrap_or_default();
            format!("{} ({way}{q})", t.term)
        })
        .collect();
    user += &tail(c, &calls, "EACH METHOD'S TOP", &terms, "GENE SETS");
    chat(SYSTEM, user.trim_end(), "LABEL:")
}

/// The cache key of a prompt for a model at a revision, under [`DECODE`].
#[must_use]
pub fn key(prompt: &str, model: &str, revision: &str) -> String {
    let mut h = Sha256::new();
    for part in [prompt, model, revision, DECODE] {
        h.update(part.as_bytes());
        h.update([0]);
    }
    h.finalize().iter().map(|b| format!("{b:02x}")).collect()
}

/// Whether a sentence passed the check or is the template.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Check {
    Ok,
    Template,
}

/// One cluster's proposal: a line of the cache.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Proposal {
    pub cluster: ClusterId,
    pub label: String,
    /// Each choice's probability.
    pub p: BTreeMap<String, f32>,
    pub rationale: String,
    /// The cluster's entities the rationale names.
    pub cited: Vec<String>,
    pub check: Check,
    pub model: String,
    pub revision: String,
    pub key: String,
}

impl Proposal {
    /// The probability of the label proposed.
    #[must_use]
    pub fn p_label(&self) -> f32 {
        self.p.get(&self.label).copied().unwrap_or(0.0)
    }
}

/// The entities of `named` and `own` that `text` names, in their order.
pub(crate) fn cited(text: &str, named: &BTreeSet<String>, own: &BTreeSet<String>) -> Vec<String> {
    let t = super::norm(text);
    named
        .iter()
        .chain(own.iter().filter(|o| !named.contains(*o)))
        .filter(|e| mentions(&t, &super::norm(e)))
        .cloned()
        .collect()
}

/// The check a model's sentence passes beyond [`citation_check`]: no
/// numbers (nor formulas or markup standing in for them), other than in the
/// names it may cite (`CD8`), and some text.
pub fn check_lm(
    text: &str,
    run: &BTreeSet<String>,
    named: &BTreeSet<String>,
    own: &BTreeSet<String>,
) -> Result<(), CitationError> {
    if text.trim().is_empty() {
        return Err(CitationError::Empty);
    }
    let mut rest = super::norm(text);
    let mut names: Vec<String> = named.iter().chain(own).map(|e| super::norm(e)).collect();
    // The longest first, so a name inside another goes with it.
    names.sort_by_key(|n| std::cmp::Reverse(n.len()));
    for n in names.iter().filter(|n| !n.is_empty()) {
        rest = rest.replace(n.as_str(), " ");
    }
    if let Some(d) = rest
        .split(|c: char| !c.is_numeric())
        .find(|w| !w.is_empty())
    {
        return Err(CitationError::Digits(d.to_string()));
    }
    // Formulas and markup stand in for the numbers a small model was told
    // to leave out.
    if let Some(m) = rest.chars().find(|c| MARKUP.contains(*c)) {
        return Err(CitationError::Markup(m.to_string()));
    }
    citation_check(text, run, named, own).map(|_| ())
}

/// A sentence from `lm` begun with `lead` after `prompt`, decoded without
/// numbers or markup (but the digits of the names in `keep`) and checked by
/// `check`; retried once with the offending word's first token banned (or,
/// after a number, the names' digits too), `None` when it fails again or
/// names no one word to ban.
pub(crate) fn checked(
    lm: &mut dyn Lm,
    prompt: &str,
    lead: &str,
    max_tokens: usize,
    keep: Vec<String>,
    check: impl Fn(&str) -> Result<(), CitationError>,
) -> Result<Option<String>> {
    let asked = format!("{prompt}{lead}");
    let text = |cont: String| format!("{lead}{cont}").trim().to_string();
    let mut bans = Bans {
        words: Vec::new(),
        numbers: true,
        keep,
    };
    let draft = text(lm.generate(&asked, max_tokens, &bans)?);
    match check(&draft) {
        Ok(()) => return Ok(Some(draft)),
        Err(CitationError::Foreign(w) | CitationError::Markup(w)) => bans.words.push(w),
        Err(CitationError::Digits(_)) => bans.keep.clear(),
        Err(e) => {
            log::debug!("describe: {e}: {draft:?}");
            return Ok(None);
        }
    }
    let again = text(lm.generate(&asked, max_tokens, &bans)?);
    match check(&again) {
        Ok(()) => Ok(Some(again)),
        Err(e) => {
            log::debug!("describe: {e}, again: {again:?}");
            Ok(None)
        }
    }
}

/// The names of `named` and `own` with a digit, whose tokens a sentence
/// may use.
pub(crate) fn numbered(named: &BTreeSet<String>, own: &BTreeSet<String>) -> Vec<String> {
    (named.iter().chain(own))
        .filter(|n| n.chars().any(char::is_numeric))
        .map(|n| shown(n))
        .collect()
}

/// The first sentence of `text`: the model was asked for one.
fn first_sentence(text: &str) -> String {
    match text.find(". ") {
        Some(i) => text[..=i].to_string(),
        None => text.to_string(),
    }
}

/// The template sentence for `label` among `c`'s choices.
fn template(c: &ClusterInput, label: &str, p: &BTreeMap<String, f32>) -> String {
    let ch = c.choices.iter().find(|x| x.label() == label);
    let q = ch.and_then(|x| x.evidence.q);
    let second = c
        .choices
        .iter()
        .filter(|x| x.label() != label && x.label() != UNASSIGNED_LABEL)
        .filter_map(|x| Some((x, p.get(x.label()).copied()?, x.evidence.q?)))
        .filter(|(.., q)| *q < 0.1)
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(x, _, q)| (x.label().to_string(), q));
    let id = format!("K{}", c.cluster);
    // The cluster is annotated so only when the label is already its own.
    let lead = (c.label.as_deref() != Some(label))
        .then(|| format!("The consensus suggests {label} for {id}."));
    template_led(
        &ClusterEvidence {
            id,
            coarse_label: label.to_string(),
            best_label: label.to_string(),
            best_q: q,
            coarse_q: None,
            label_support: p.get(label).copied(),
            best_significant: label != UNASSIGNED_LABEL,
            second,
            markers: ch.map(|x| x.markers.clone()).unwrap_or_default(),
            neighbour_genes: Vec::new(),
            incidence_words: Vec::new(),
        },
        lead.as_deref(),
    )
}

/// The proposal for cluster `c` of a round whose types and terms are
/// `run`, by `lm` (named `model` at `revision`).
pub fn propose(
    lm: &mut dyn Lm,
    model: &str,
    revision: &str,
    c: &ClusterInput,
    run: &BTreeSet<String>,
) -> Result<Proposal> {
    let conts: Vec<String> = c
        .choices
        .iter()
        .map(|ch| format!(" {}\n", shown(ch.label())))
        .collect();
    let scores = lm.score(&label_prompt(c), &conts)?;
    anyhow::ensure!(scores.len() == conts.len(), "a score per choice");
    let top = scores.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let w: Vec<f32> = scores.iter().map(|s| (s - top).exp()).collect();
    let z: f32 = w.iter().sum();
    let mut best = 0;
    for (i, s) in scores.iter().enumerate() {
        if *s > scores[best] {
            best = i;
        }
    }
    let label = c.choices[best].label().to_string();
    let p: BTreeMap<String, f32> = c
        .choices
        .iter()
        .zip(&w)
        .map(|(ch, w)| (ch.label().to_string(), w / z))
        .collect();
    let (named, own) = (c.named(), c.own());
    // The sentence starts at the label, so it names it.
    let why = format!("{} {}\nWHY:", why_prompt(c), shown(&label));
    let lead = if label == UNASSIGNED_LABEL {
        String::new()
    } else {
        format!(" {}", shown(&label))
    };
    let sentence = checked(
        lm,
        &why,
        &lead,
        RATIONALE_TOKENS,
        numbered(&named, &own),
        |t| {
            check_lm(t, run, &named, &own)?;
            super::claims::rank_claims(t, c, &label).map_err(CitationError::Claim)
        },
    )?;
    let (rationale, check) = match sentence {
        Some(s) => (first_sentence(&s), Check::Ok),
        None => (template(c, &label, &p), Check::Template),
    };
    Ok(Proposal {
        cluster: c.cluster,
        cited: cited(&rationale, &named, &own),
        label,
        p,
        rationale,
        check,
        model: model.to_string(),
        revision: revision.to_string(),
        key: key(&prompt(c), model, revision),
    })
}

/// `{stem}.consensus.jsonl`: every proposal made for the round (or the
/// trajectory) at `stem`, one per line, appended and never rewritten.
#[must_use]
pub fn cache_path(stem: &str) -> PathBuf {
    PathBuf::from(format!("{stem}.consensus.jsonl"))
}

/// The records of `path` that read as `T`; none when it is missing.
pub fn read_cache<T: DeserializeOwned>(path: &Path) -> Vec<T> {
    std::fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .filter_map(|l| serde_json::from_str(l).ok())
        .collect()
}

/// Append `rec` to `path` as one line, in one write.
pub fn append<T: Serialize>(path: &Path, rec: &T) -> Result<()> {
    let line = serde_json::to_string(rec)? + "\n";
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .with_context(|| format!("opening {}", path.display()))?;
    f.write_all(line.as_bytes())
        .with_context(|| format!("appending to {}", path.display()))
}

/// Each cluster's prompt.
#[must_use]
pub fn prompts(round: &Round) -> BTreeMap<ClusterId, String> {
    round
        .clusters
        .iter()
        .map(|c| (c.cluster, prompt(c)))
        .collect()
}

/// Each cluster's latest proposal in `cache` made on the evidence its
/// prompt in `prompts` states, by any model.
#[must_use]
pub fn current(
    prompts: &BTreeMap<ClusterId, String>,
    cache: &Path,
) -> BTreeMap<ClusterId, Proposal> {
    read_cache::<Proposal>(cache)
        .into_iter()
        .filter(|p| {
            prompts
                .get(&p.cluster)
                .is_some_and(|pr| key(pr, &p.model, &p.revision) == p.key)
        })
        .map(|p| (p.cluster, p))
        .collect()
}

const STAGES: crate::progress::Stages = crate::progress::Stages::new(&[("describe", 1.0)]);

/// For each of `keys` in turn, the record of `cache` with that key, else
/// `make(lm, k)` by the model `load` gives, loaded at the first miss only;
/// `what` the items are, for the progress.
pub(crate) fn cached_or<T: DeserializeOwned + Clone>(
    cache: &Path,
    keys: impl ExactSizeIterator<Item = String>,
    what: &str,
    load: impl FnOnce() -> Result<Box<dyn Lm>>,
    mut make: impl FnMut(&mut dyn Lm, usize) -> Result<T>,
    key_of: impl Fn(&T) -> String,
) -> Result<Vec<T>> {
    let done: HashMap<String, T> = read_cache::<T>(cache)
        .into_iter()
        .map(|t| (key_of(&t), t))
        .collect();
    let mut load = Some(load);
    let mut lm: Option<Box<dyn Lm>> = None;
    let n = keys.len();
    let mut out = Vec::with_capacity(n);
    STAGES.start(0);
    for (k, key) in keys.enumerate() {
        let t = match done.get(&key) {
            Some(t) => t.clone(),
            None => {
                if lm.is_none() {
                    lm = Some(load.take().context("the model is loaded once")?()?);
                }
                make(lm.as_deref_mut().context("a model")?, k)?
            }
        };
        out.push(t);
        STAGES.within(0, k + 1, n, Some(&format!("{} of {n} {what}", k + 1)));
    }
    STAGES.finish();
    Ok(out)
}

/// A proposal for every cluster of `round` by the model `spec` names: from
/// `cache` when it has one made on the same evidence, else by the model,
/// loaded (by `load`) only then, and appended to `cache`.
pub fn propose_round(
    round: &Round,
    spec: &LmSpec,
    cache: &Path,
    load: impl FnOnce() -> Result<Box<dyn Lm>>,
) -> Result<Vec<Proposal>> {
    let (model, revision) = (spec.name(), spec.revision());
    let keys: Vec<String> = (round.clusters.iter())
        .map(|c| key(&prompt(c), &model, &revision))
        .collect();
    let out = cached_or(
        cache,
        keys.into_iter(),
        "clusters",
        load,
        |lm, k| {
            let p = propose(lm, &model, &revision, &round.clusters[k], &round.entities)?;
            append(cache, &p)?;
            Ok(p)
        },
        |p: &Proposal| p.key.clone(),
    )?;
    for p in &out {
        let template = if p.check == Check::Template {
            " · template"
        } else {
            ""
        };
        log::info!("K{}: {} ({:.2}){template}", p.cluster, p.label, p.p_label());
    }
    Ok(out)
}

#[cfg(test)]
#[path = "tests/consensus.rs"]
pub(crate) mod tests;
