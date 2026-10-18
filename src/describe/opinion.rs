//! A small local language model's opinion of each cluster of a round: a
//! method beside the enrichment (with SuSiE) and the projection score
//! ([`Method::Llm`]), run by `lupin describe --lm` (in the TUI, `L`) and
//! never by a pass.
//!
//! The model ([`super::lm`]) reads what a person would look at first: the
//! cluster's size, its top specific genes (as the TUI's genes pane ranks
//! them, [`crate::tui::round::specific_genes`], on the panel or not) and its
//! gene-set terms; never another method's reading, the round's label or the
//! Cell Ontology's placement. It scores every type of the panel, and
//! `unassigned`, as the continuation of `LABEL:`: a softmax of the scores
//! is its probability of each, so it names no type off the panel. Its top
//! is offered in the TUI (◆) as a suggestion the user confirms.
//!
//! The facts about that type (where the other methods put it, with their
//! readings) are written by code from the evidence ([`facts`]), never by
//! the model. The model adds one sentence on what the cluster's genes and
//! terms indicate for the type (its biology), from those alone, decoded
//! without numbers, formula signs or the words that state a rank, an order,
//! a role, a verdict or a method ([`super::vocabulary`]); the sentence may
//! name the type, its markers, the top genes and the cluster's terms, and
//! no other type of the run ([`super::citation_check`]). A sentence that
//! fails is asked once more with the offender banned, then dropped: the
//! facts stand alone.
//!
//! Prior art for language-model cell typing: Hou & Ji, *Assessing GPT-4 for
//! cell type annotation in single-cell RNA-seq analysis*, Nat Methods
//! 21:1462 (2024; GPTCelltype).

use super::lm::{Bans, Lm, LmSpec, MARKUP};
use super::{citation_check, mentions, CitationError};
use crate::annotate::markers::{label_key, read_marker_pairs};
use crate::annotate::rounds::{ClusterId, Method, Table, Term, MIN_LLM_P};
use crate::manifest::run::{derive_out_prefix, Loaded};
use crate::tui::round::Candidate;
use anyhow::{Context, Result};
use enrichment::UNASSIGNED_LABEL;
use serde::{Deserialize, Serialize};
use sha2::{Digest as _, Sha256};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::Write;
use std::path::{Path, PathBuf};

/// The cluster's specific genes the model reads.
pub const TOP_GENES: usize = 50;
/// Terms, and each type's markers, the prompts quote.
const TOP: usize = 5;
/// The longest biology sentence, in tokens.
pub const BIOLOGY_TOKENS: usize = 60;
/// How the model is asked; part of the cache key.
const DECODE: &str = "label: log-probability of each panel type and unassigned, softmax; \
                      biology: from the label, greedy to a newline, 60 tokens, no numbers, markup \
                      or evidence words but the names', one retry with the offender banned";

/// A type of the round's panel, with its first markers in the panel's
/// order.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PanelType {
    pub label: String,
    pub markers: Vec<String>,
}

/// One cluster as the model reads it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ClusterInput {
    pub cluster: ClusterId,
    pub size: usize,
    /// Its first [`TOP_GENES`] specific genes, most specific first.
    pub genes: Vec<String>,
    /// Each panel type's markers among all the cluster's specific genes,
    /// most specific first, up to [`TOP`]; types with none left out.
    pub here: BTreeMap<String, Vec<String>>,
    pub terms: Vec<Term>,
}

/// A round as the model reads it: its clusters, its panel, and every type
/// and term in it, which a cluster's sentence may name only when they are
/// its own.
#[derive(Clone, Debug, Default)]
pub struct Round {
    pub clusters: Vec<ClusterInput>,
    pub panel: Vec<PanelType>,
    pub entities: BTreeSet<String>,
}

/// The round `loaded`: each cluster's size, specific genes and terms, as
/// `review` and the TUI read them, and the panel's types.
pub fn inputs(loaded: &Loaded) -> Result<Round> {
    use crate::tui::round::{specific_genes, Expression};
    use legume_numeric::matrix::dense_mat_io::Mat;
    use legume_numeric::matrix::traits::IoOps;
    let r = crate::manifest::rounds::Reviewed::read(loaded, false)?;
    let a = &loaded.manifest.annotate;
    let pairs = match a.markers.as_deref().and_then(|m| loaded.recorded(m)) {
        Some(p) => read_marker_pairs(&p)?,
        None => Vec::new(),
    };
    anyhow::ensure!(
        !pairs.is_empty(),
        "the model chooses among the round's panel types, and the round records no panel that reads"
    );
    let mut panel: Vec<PanelType> = Vec::new();
    let mut sets: HashMap<String, BTreeSet<String>> = HashMap::new();
    for (gene, label) in &pairs {
        let at = match panel.iter().position(|t| *t.label == **label) {
            Some(i) => i,
            None => {
                panel.push(PanelType {
                    label: label.to_string(),
                    markers: Vec::new(),
                });
                panel.len() - 1
            }
        };
        let markers = &mut panel[at].markers;
        if markers.len() < TOP && !markers.iter().any(|m| **m == **gene) {
            markers.push(gene.to_string());
        }
        sets.entry(label.to_string())
            .or_default()
            .insert(gene.to_string());
    }
    let specific = match (a.cluster_expression.as_deref()).and_then(|e| loaded.recorded(e)) {
        Some(p) => specific_genes(&Expression::new(
            Mat::from_parquet_with_row_names(&p, Some(0))
                .with_context(|| format!("reading {p}"))?,
        )),
        None => BTreeMap::new(),
    };
    let mut entities: BTreeSet<String> = panel.iter().map(|t| t.label.clone()).collect();
    let mut clusters = Vec::new();
    for (&id, d) in &r.digests {
        entities.extend(d.terms.iter().map(|t| t.term.clone()));
        let genes = specific.get(&id).map_or(&[][..], Vec::as_slice);
        let here = (panel.iter())
            .filter_map(|t| {
                let set = sets.get(&t.label)?;
                let m: Vec<String> = (genes.iter())
                    .filter(|(g, _)| set.contains(g))
                    .take(TOP)
                    .map(|(g, _)| g.clone())
                    .collect();
                (!m.is_empty()).then(|| (t.label.clone(), m))
            })
            .collect();
        clusters.push(ClusterInput {
            cluster: id,
            size: d.size,
            genes: genes
                .iter()
                .take(TOP_GENES)
                .map(|(g, _)| g.clone())
                .collect(),
            here,
            terms: d.terms.iter().take(TOP).cloned().collect(),
        });
    }
    Ok(Round {
        clusters,
        panel,
        entities,
    })
}

/// A label as the prompts write it: words, not the panel's underscores.
fn shown(label: &str) -> String {
    label.replace('_', " ")
}

/// Qwen3's chat with thinking switched off: `system`, `user`, and the
/// assistant's answer begun with `answer`.
pub(crate) fn chat(system: &str, user: &str, answer: &str) -> String {
    format!(
        "<|im_start|>system\n{system}<|im_end|>\n<|im_start|>user\n{user}<|im_end|>\n\
         <|im_start|>assistant\n<think>\n\n</think>\n\n{answer}"
    )
}

const SYSTEM: &str = "You name the cell type of one single-cell cluster from its genes.\n\
                      Choose only from TYPES. Answer with the type alone.";

/// The lines on a cluster's genes and terms both prompts share.
fn genes_and_terms(c: &ClusterInput) -> String {
    let mut s = String::new();
    if !c.genes.is_empty() {
        s += &format!("TOP GENES: {}\n", c.genes.join(", "));
    }
    if !c.terms.is_empty() {
        let terms: Vec<String> = (c.terms.iter())
            .map(|t| {
                let way = if t.effect >= 0.0 { "up" } else { "down" };
                format!("{} ({way})", t.term)
            })
            .collect();
        s += &format!("GENE SETS: {}\n", terms.join("; "));
    }
    s
}

/// The prompt the type is chosen after, up to `LABEL:`: the panel's types,
/// each with its first markers, then cluster `c`'s size, top genes and
/// terms.
#[must_use]
pub fn opinion_prompt(c: &ClusterInput, panel: &[PanelType]) -> String {
    let mut user = String::from("TYPES:\n");
    for t in panel {
        user += &format!("- {}", shown(&t.label));
        if !t.markers.is_empty() {
            user += &format!(": {}", t.markers.join(", "));
        }
        user += "\n";
    }
    user += &format!("- {UNASSIGNED_LABEL}: no type fits\n");
    user += &format!("CLUSTER: {} cells.\n", c.size);
    user += &genes_and_terms(c);
    chat(SYSTEM, user.trim_end(), "LABEL:")
}

const BIOLOGY: &str = "You explain in one plain sentence what a single-cell cluster's genes say \
                       about its cell type: what the genes do, and why that fits the type.\n\
                       Write a sentence, not a list: name at most three genes. No numbers. Name \
                       only the cell type, genes and gene sets below.";

/// The prompt the biology sentence on `label` for cluster `c` is written
/// after: the type, its markers among the cluster's genes, the top genes
/// and the terms; no method's reading, rank or order.
#[must_use]
pub fn biology_prompt(c: &ClusterInput, label: &str) -> String {
    let mut user = format!("CELL TYPE: {}\n", shown(label));
    if let Some(m) = c.here.get(label) {
        user += &format!("ITS MARKERS HERE: {}\n", m.join(", "));
    }
    user += &genes_and_terms(c);
    chat(BIOLOGY, user.trim_end(), "")
}

/// What the model is asked about cluster `c`, for the cache key: the
/// opinion's prompt and the biology's for every type it may choose.
#[must_use]
pub fn prompt(c: &ClusterInput, panel: &[PanelType]) -> String {
    std::iter::once(opinion_prompt(c, panel))
        .chain(panel.iter().map(|t| biology_prompt(c, &t.label)))
        .collect::<Vec<_>>()
        .join("\n")
}

/// The cache key of a prompt for a model at a revision, under [`DECODE`]
/// and the words it bans ([`super::vocabulary`]).
#[must_use]
pub fn key(prompt: &str, model: &str, revision: &str) -> String {
    let mut h = Sha256::new();
    let words = super::vocabulary::all().collect::<Vec<_>>().join(",");
    for part in [prompt, model, revision, DECODE, &words] {
        h.update(part.as_bytes());
        h.update([0]);
    }
    h.finalize().iter().map(|b| format!("{b:02x}")).collect()
}

/// Whether the biology sentence passed its check (`ok`), or none stands
/// beside the facts (`template`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Check {
    Ok,
    Template,
}

/// One cluster's opinion: a line of the cache.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Opinion {
    pub cluster: ClusterId,
    /// The type of highest probability, or `unassigned`.
    pub label: String,
    /// Each panel type's probability, and `unassigned`'s.
    pub p: BTreeMap<String, f32>,
    /// What the cluster's genes and terms indicate for the label, by the
    /// model; none when its sentence failed, or for `unassigned`.
    #[serde(default)]
    pub biology: Option<String>,
    /// The entities the sentence names.
    pub cited: Vec<String>,
    pub check: Check,
    pub model: String,
    pub revision: String,
    pub key: String,
}

impl Opinion {
    /// The probability of the label.
    #[must_use]
    pub fn p_label(&self) -> f32 {
        self.p.get(&self.label).copied().unwrap_or(0.0)
    }

    /// The text beside the label: the facts by code on the candidates of
    /// its cluster, then the model's sentence.
    #[must_use]
    pub fn text(&self, candidates: &[Candidate]) -> String {
        let facts = facts(&self.label, candidates);
        match &self.biology {
            Some(b) => format!("{facts} {b}"),
            None => facts,
        }
    }
}

/// A method in the facts' words.
fn who(m: Method) -> &'static str {
    match m {
        Method::Enrichment => "enrichment",
        Method::Susie => "SuSiE",
        Method::Score => "the score",
        Method::Llm => "the model",
    }
}

/// A rank in words.
fn ordinal(rank: usize) -> String {
    const WORDS: [&str; 10] = [
        "first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth",
        "tenth",
    ];
    WORDS
        .get(rank)
        .map_or_else(|| format!("#{}", rank + 1), |w| (*w).to_string())
}

/// A number for the facts, in two significant digits.
fn num(v: f32) -> String {
    if v != 0.0 && v.abs() < 0.01 {
        format!("{v:.0e}")
    } else {
        format!("{v:.2}")
    }
}

/// Method `m`'s readings of `c`, in brackets.
fn readings(m: Method, c: &Candidate) -> String {
    let mut parts: Vec<String> = Vec::new();
    match m {
        Method::Enrichment => {
            parts.extend(c.enrichment.map(|v| format!("share {}", num(v))));
            parts.extend(c.q.map(|q| format!("q {}", num(q))));
        }
        Method::Susie => {
            parts.extend(c.pip.map(|v| format!("PIP {}", num(v))));
            parts.extend(c.fold.map(|f| format!("×{f:.1}")));
        }
        Method::Score => parts.extend(c.score.map(|v| format!("share {}", num(v)))),
        Method::Llm => {}
    }
    if parts.is_empty() {
        String::new()
    } else {
        format!(" ({})", parts.join(", "))
    }
}

/// The facts on `label` in a cluster of `candidates`, written by code from
/// the evidence: where each other method that ranks the cluster's types
/// puts it, with its readings (for `unassigned`, each one's first). True of
/// the evidence whatever the model thought.
#[must_use]
pub fn facts(label: &str, candidates: &[Candidate]) -> String {
    let ran: Vec<Method> = (Method::ALL.into_iter())
        .filter(|&m| m != Method::Llm)
        .filter(|m| candidates.iter().any(|c| c.ranks[m.index()].is_some()))
        .collect();
    if ran.is_empty() {
        return "No other method ranks this cluster's types.".into();
    }
    if label == UNASSIGNED_LABEL {
        let firsts: Vec<String> = (ran.iter())
            .filter_map(|&m| {
                let top = candidates.iter().find(|c| c.ranks[m.index()] == Some(0))?;
                Some(format!("{} by {}", top.label, who(m)))
            })
            .collect();
        return format!("{UNASSIGNED_LABEL}; first are {}.", firsts.join(", "));
    }
    let c = (candidates.iter()).find(|c| label_key(&c.label) == label_key(label));
    let parts: Vec<String> = (ran.iter())
        .map(|&m| match c.and_then(|c| Some((c, c.ranks[m.index()]?))) {
            Some((c, r)) => format!("{} by {}{}", ordinal(r), who(m), readings(m, c)),
            None => format!("not ranked by {}", who(m)),
        })
        .collect();
    format!("{label}: {}.", parts.join("; "))
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

/// The check a biology sentence passes beyond [`citation_check`]: some
/// text, and outside the names it may cite (`CD8`) no number, formula or
/// markup sign, nor a word of the [`super::vocabulary`]; a sentence, not a
/// list of more than [`MAX_NAMED`] of its names.
pub fn check(
    text: &str,
    run: &BTreeSet<String>,
    named: &BTreeSet<String>,
    own: &BTreeSet<String>,
) -> Result<(), CitationError> {
    if text.trim().is_empty() {
        return Err(CitationError::Empty);
    }
    // Another's name first: a digit in it (`CT2`) is no number, and the
    // retry is to ban the name.
    let cites = citation_check(text, run, named, own);
    if let Err(e @ CitationError::Foreign(_)) = cites {
        return Err(e);
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
    if let Some(w) = super::vocabulary::used(&rest) {
        return Err(CitationError::Vocabulary(w.to_string()));
    }
    // A sentence, not a list of the genes it was given.
    let own_named = cited(text, &BTreeSet::new(), own).len();
    if own_named > MAX_NAMED {
        return Err(CitationError::List(own_named));
    }
    cites.map(|_| ())
}

/// The most genes, terms and markers a biology sentence names: more is a
/// list of what it was given.
const MAX_NAMED: usize = 5;

/// How a biology sentence on `label` begins: the type's cells, so it reads
/// as a sentence about them.
pub(crate) fn lead(label: &str) -> String {
    let shown = shown(label);
    let last = shown.rsplit(' ').next().unwrap_or_default().to_lowercase();
    if last == "cell" || last == "cells" {
        shown
    } else {
        format!("{shown} cells")
    }
}

/// A sentence from `lm` begun with `lead` after `prompt`, decoded under
/// `bans` and checked by `check`; retried once with the offending word's
/// first token banned (or, after a number, the names' digits too), `None`
/// when it fails again or names no one word to ban.
pub(crate) fn checked(
    lm: &mut dyn Lm,
    prompt: &str,
    lead: &str,
    mut bans: Bans,
    check: impl Fn(&str) -> Result<(), CitationError>,
) -> Result<Option<String>> {
    let asked = format!("{prompt}{lead}");
    let text = |cont: String| format!("{lead}{cont}").trim().to_string();
    let draft = text(lm.generate(&asked, BIOLOGY_TOKENS, &bans)?);
    match check(&draft) {
        Ok(()) => return Ok(Some(draft)),
        Err(
            CitationError::Foreign(w) | CitationError::Markup(w) | CitationError::Vocabulary(w),
        ) => {
            // A phrase by its first word, as the bans start it.
            let w = w.split(' ').next().unwrap_or_default().to_string();
            bans.words.push(w);
        }
        Err(CitationError::Digits(_)) => bans.keep.clear(),
        Err(e) => {
            log::debug!("describe: {e}: {draft:?}");
            return Ok(None);
        }
    }
    let again = text(lm.generate(&asked, BIOLOGY_TOKENS, &bans)?);
    match check(&again) {
        Ok(()) => Ok(Some(again)),
        Err(e) => {
            log::debug!("describe: {e}, again: {again:?}");
            Ok(None)
        }
    }
}

/// The first sentence of `text`: the model was asked for one.
pub(crate) fn first_sentence(text: &str) -> String {
    match text.find(". ") {
        Some(i) => text[..=i].to_string(),
        None => text.to_string(),
    }
}

/// The opinion on cluster `c` of `round` by `lm` (named `model` at
/// `revision`): every panel type and `unassigned` scored by its
/// log-probability as the answer, a softmax of those its probability; then
/// the biology sentence on the top.
pub fn opine(
    lm: &mut dyn Lm,
    model: &str,
    revision: &str,
    c: &ClusterInput,
    round: &Round,
) -> Result<Opinion> {
    let labels: Vec<&str> = (round.panel.iter())
        .map(|t| t.label.as_str())
        .chain([UNASSIGNED_LABEL])
        .collect();
    let conts: Vec<String> = labels.iter().map(|l| format!(" {}\n", shown(l))).collect();
    let scores = lm.score(&opinion_prompt(c, &round.panel), &conts)?;
    anyhow::ensure!(scores.len() == conts.len(), "a score per type");
    let top = scores.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let w: Vec<f32> = scores.iter().map(|s| (s - top).exp()).collect();
    let z: f32 = w.iter().sum();
    let mut best = 0;
    for (i, s) in scores.iter().enumerate() {
        if *s > scores[best] {
            best = i;
        }
    }
    let label = labels[best].to_string();
    let p: BTreeMap<String, f32> = (labels.iter().zip(&w))
        .map(|(l, w)| ((*l).to_string(), w / z))
        .collect();
    let (biology, cited) = if label == UNASSIGNED_LABEL {
        (None, Vec::new())
    } else {
        let named = BTreeSet::from([label.clone()]);
        let mut own: BTreeSet<String> = c.genes.iter().cloned().collect();
        own.extend(c.here.get(&label).into_iter().flatten().cloned());
        own.extend(c.terms.iter().map(|t| t.term.clone()));
        let names: Vec<String> = named.iter().chain(&own).map(|n| shown(n)).collect();
        let lower: Vec<String> = names.iter().map(|n| n.to_lowercase()).collect();
        let bans = Bans {
            words: super::vocabulary::bans(&lower),
            numbers: true,
            keep: names,
        };
        let sentence = checked(lm, &biology_prompt(c, &label), &lead(&label), bans, |t| {
            check(t, &round.entities, &named, &own)
        })?
        .map(|s| first_sentence(&s));
        let cited = sentence
            .as_deref()
            .map_or_else(Vec::new, |s| cited(s, &named, &own));
        (sentence, cited)
    };
    Ok(Opinion {
        cluster: c.cluster,
        check: if biology.is_some() {
            Check::Ok
        } else {
            Check::Template
        },
        biology,
        cited,
        label,
        p,
        model: model.to_string(),
        revision: revision.to_string(),
        key: key(&prompt(c, &round.panel), model, revision),
    })
}

/// `{stem}.llm.jsonl`: every opinion made on the round at `stem`, one per
/// line, appended and never rewritten. It lives beside the round, whose
/// review and TUI read it as the model's evidence.
#[must_use]
pub fn cache_path(stem: &str) -> PathBuf {
    PathBuf::from(format!("{stem}.llm.jsonl"))
}

/// The cache of the round `loaded`.
#[must_use]
pub fn cache_of(loaded: &Loaded) -> PathBuf {
    cache_path(&derive_out_prefix(&loaded.file.to_string_lossy()))
}

/// The opinions of `path` that read; none when it is missing.
fn read_cache(path: &Path) -> Vec<Opinion> {
    std::fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .filter_map(|l| serde_json::from_str(l).ok())
        .collect()
}

/// Append `o` to `path` as one line, in one write.
pub fn append(path: &Path, o: &Opinion) -> Result<()> {
    let line = serde_json::to_string(o)? + "\n";
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
    (round.clusters.iter())
        .map(|c| (c.cluster, prompt(c, &round.panel)))
        .collect()
}

/// Each cluster's latest opinion in `cache` made on what its prompt in
/// `prompts` states, by any model.
#[must_use]
pub fn current(
    prompts: &BTreeMap<ClusterId, String>,
    cache: &Path,
) -> BTreeMap<ClusterId, Opinion> {
    read_cache(cache)
        .into_iter()
        .filter(|o| {
            prompts
                .get(&o.cluster)
                .is_some_and(|pr| key(pr, &o.model, &o.revision) == o.key)
        })
        .map(|o| (o.cluster, o))
        .collect()
}

/// The opinions as the model's evidence, cluster × type: each type's
/// probability where it is at least [`MIN_LLM_P`] (what the model puts
/// forward), `unassigned` left out.
#[must_use]
pub fn table(opinions: &BTreeMap<ClusterId, Opinion>) -> Table {
    let cols: Vec<String> = (opinions.values())
        .flat_map(|o| o.p.keys())
        .filter(|l| *l != UNASSIGNED_LABEL)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .cloned()
        .collect();
    let values = (opinions.values())
        .flat_map(|o| {
            cols.iter().map(|l| match o.p.get(l) {
                Some(&p) if p >= MIN_LLM_P => p,
                _ => f32::NAN,
            })
        })
        .collect();
    Table {
        rows: opinions.keys().copied().collect(),
        cols,
        values,
    }
}

/// The model's evidence on the round `loaded`, from its cache beside it:
/// `None` without one.
pub fn read_table(loaded: &Loaded) -> Result<Option<Table>> {
    let cache = cache_of(loaded);
    if !cache.is_file() {
        return Ok(None);
    }
    let round = inputs(loaded)?;
    Ok(Some(table(&current(&prompts(&round), &cache))))
}

const STAGES: crate::progress::Stages = crate::progress::Stages::new(&[("describe", 1.0)]);

/// An opinion on every cluster of `round` by the model `spec` names: from
/// `cache` when it has one made on the same inputs, else by the model,
/// loaded (by `load`) at the first miss only, and appended to `cache`.
pub fn opine_round(
    round: &Round,
    spec: &LmSpec,
    cache: &Path,
    load: impl FnOnce() -> Result<Box<dyn Lm>>,
) -> Result<Vec<Opinion>> {
    let (model, revision) = (spec.name(), spec.revision());
    let done: HashMap<String, Opinion> = (read_cache(cache).into_iter())
        .map(|o| (o.key.clone(), o))
        .collect();
    let mut load = Some(load);
    let mut lm: Option<Box<dyn Lm>> = None;
    let n = round.clusters.len();
    let mut out = Vec::with_capacity(n);
    STAGES.start(0);
    for (k, c) in round.clusters.iter().enumerate() {
        let o = match done.get(&key(&prompt(c, &round.panel), &model, &revision)) {
            Some(o) => o.clone(),
            None => {
                if lm.is_none() {
                    lm = Some(load.take().context("the model is loaded once")?()?);
                }
                let o = opine(
                    lm.as_deref_mut().context("a model")?,
                    &model,
                    &revision,
                    c,
                    round,
                )?;
                append(cache, &o)?;
                o
            }
        };
        out.push(o);
        STAGES.within(0, k + 1, n, Some(&format!("{} of {n} clusters", k + 1)));
    }
    STAGES.finish();
    for o in &out {
        let template = if o.check == Check::Template {
            " · facts only"
        } else {
            ""
        };
        log::info!("K{}: {} ({:.2}){template}", o.cluster, o.label, o.p_label());
    }
    Ok(out)
}

#[cfg(test)]
#[path = "tests/opinion.rs"]
pub(crate) mod tests;
