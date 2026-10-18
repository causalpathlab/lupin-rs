//! A small local language model's critique of a round's two readings of
//! each cluster: the round's call (the enrichment, with SuSiE when it ran)
//! and the projection score's top type, or, in a round without the score,
//! the other of the enrichment's and SuSiE's tops. Run by `lupin describe
//! --lm` (in the TUI, `L`) and never by a pass; it decides nothing.
//!
//! Code states every fact: which method put forward which type, and
//! whether they agree ([`Critique::facts`]). The model is asked one closed
//! question on the cluster's genes and terms and the two types with their
//! first markers, never a method's reading or rank: when the two agree,
//! whether the genes fit the type (`yes` / `no`); when they differ, which
//! fits (`A`, `B`, `both`: a mix or a state between them, `neither`). Each
//! answer is scored by its whole log-probability, a softmax of those its
//! probability. The two types are named `A` and `B` in both orders and the
//! probabilities averaged, so neither wins by its place.
//!
//! When the answer favours a type, the model adds one sentence on what the
//! cluster's genes indicate for it, under [`super::opinion`]'s bans and
//! checks, naming only the type and its markers among the cluster's top
//! genes; once more with the offender banned, then dropped.
//!
//! Prior art for language-model cell typing: Hou & Ji, *Assessing GPT-4 for
//! cell type annotation in single-cell RNA-seq analysis*, Nat Methods
//! 21:1462 (2024; GPTCelltype).

use super::lm::{Bans, Lm, LmModel, LmSpec};
use super::opinion::{self, Check, ClusterInput, PanelType, Round};
use crate::annotate::markers::label_key;
use crate::annotate::rounds::{ClusterId, Method};
use crate::manifest::run::{derive_out_prefix, Loaded};
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use sha2::{Digest as _, Sha256};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::Write;
use std::path::{Path, PathBuf};

/// How the model is asked; part of the cache key.
const DECODE: &str = "critique: log-probability of each answer, softmax, the two types in both \
                      orders averaged; biology: greedy to a newline, 60 tokens, no numbers, markup \
                      or evidence words, one retry with the offender banned";

/// What the question shows of a cluster beyond the two readings.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Ask {
    /// How many of its top specific genes.
    pub genes: usize,
    /// Those genes in an order drawn from this seed, not by rank.
    pub shuffle: Option<u64>,
    /// The sample: its species and tissue.
    pub context: Option<String>,
}

impl Ask {
    /// The first [`Ask::genes`] of `genes`, shuffled when asked.
    #[must_use]
    pub fn genes_of(&self, genes: &[String]) -> Vec<String> {
        let mut g: Vec<String> = genes.iter().take(self.genes).cloned().collect();
        if let Some(seed) = self.shuffle {
            g.sort_by_cached_key(|x| {
                let mut h = Sha256::new();
                h.update(seed.to_le_bytes());
                h.update(x.as_bytes());
                h.finalize().to_vec()
            });
        }
        g
    }
}

/// One of the two readings.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Reading {
    pub method: Method,
    pub label: String,
}

/// One cluster as the critique reads it.
#[derive(Clone, Debug)]
pub struct Input {
    pub cluster: ClusterInput,
    pub call: Reading,
    pub other: Reading,
}

impl Input {
    #[must_use]
    pub fn agree(&self) -> bool {
        label_key(&self.call.label) == label_key(&self.other.label)
    }
}

/// The method whose top stands beside the call `by`: the score's, else the
/// other of the enrichment's and SuSiE's.
fn beside(by: Method, tops: &[(Method, String)]) -> Option<Method> {
    let order: &[Method] = match by {
        Method::Score => &[Method::Susie, Method::Enrichment],
        Method::Susie => &[Method::Score, Method::Enrichment],
        _ => &[Method::Score, Method::Susie],
    };
    order
        .iter()
        .copied()
        .find(|m| tops.iter().any(|(t, _)| t == m))
}

/// Every cluster of the round `loaded` with two readings, and the round as
/// [`opinion::inputs`] reads it.
pub fn inputs(loaded: &Loaded, ask: &Ask) -> Result<(Round, Vec<Input>)> {
    let mut round = opinion::inputs(loaded)?;
    for c in &mut round.clusters {
        c.genes = ask.genes_of(&c.genes);
    }
    let r = crate::manifest::rounds::Reviewed::read(loaded, true)?;
    let by = crate::manifest::rounds::called_by(&loaded.manifest.annotate)
        .context("the critique reads the round's call, and the round records none")?;
    let mut out = Vec::new();
    for c in &round.clusters {
        let Some(d) = r.digests.get(&c.cluster) else {
            continue;
        };
        let Some(call) = d.calls.first() else {
            continue;
        };
        let types = r.tables.evidence(c.cluster);
        let tops: Vec<(Method, String)> = [Method::Enrichment, Method::Susie, Method::Score]
            .into_iter()
            .filter_map(|m| Some((m, m.ranked(&types).first()?.label.clone())))
            .collect();
        let Some(m) = beside(by, &tops) else {
            continue;
        };
        out.push(Input {
            cluster: c.clone(),
            call: Reading {
                method: by,
                label: call.label.clone(),
            },
            other: Reading {
                method: m,
                label: (tops.iter().find(|(t, _)| *t == m))
                    .map(|(_, l)| l.clone())
                    .unwrap_or_default(),
            },
        });
    }
    Ok((round, out))
}

/// A label as the prompts write it: words, not the panel's underscores.
fn shown(label: &str) -> String {
    label.replace('_', " ")
}

/// A type's first markers in the panel.
fn markers<'a>(panel: &'a [PanelType], label: &str) -> &'a [String] {
    (panel.iter())
        .find(|t| label_key(&t.label) == label_key(label))
        .map_or(&[][..], |t| t.markers.as_slice())
}

fn type_line(name: &str, label: &str, panel: &[PanelType]) -> String {
    let m = markers(panel, label);
    if m.is_empty() {
        format!("{name}: {}\n", shown(label))
    } else {
        format!("{name}: {} (markers {})\n", shown(label), m.join(", "))
    }
}

const SYSTEM: &str = "You check cell type calls for one single-cell cluster against its genes.\n\
                      Answer with one of the CHOICES alone.";

/// The cluster's sample, size, top genes and terms.
fn cluster_lines(c: &ClusterInput, ask: &Ask) -> String {
    let mut s = String::new();
    if let Some(x) = &ask.context {
        s += &format!("SAMPLE: {x}\n");
    }
    s += &format!("CLUSTER: {} cells.\n", c.size);
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

/// The answers a question takes, as the prompt lists them.
const AGREE: [&str; 2] = ["yes", "no"];
const DIFFER: [&str; 4] = ["A", "B", "both", "neither"];

/// The question on `i` with `first` named A (`true`: the call), up to
/// `ANSWER:`, in `model`'s chat.
fn question(model: LmModel, ask: &Ask, i: &Input, panel: &[PanelType], call_first: bool) -> String {
    let mut user = cluster_lines(&i.cluster, ask);
    if i.agree() {
        user += &type_line("TYPE", &i.call.label, panel);
        user += "Do the cluster's genes fit this type?\nCHOICES: yes; no";
    } else {
        let (a, b) = if call_first {
            (&i.call, &i.other)
        } else {
            (&i.other, &i.call)
        };
        user += &type_line("TYPE A", &a.label, panel);
        user += &type_line("TYPE B", &b.label, panel);
        user += "Which type do the cluster's genes fit?\n\
                 CHOICES: A; B; both (a mix of the two, or a state between them); neither";
    }
    model.chat(SYSTEM, &user, "ANSWER:")
}

const BIOLOGY: &str = "You explain in one plain sentence what a single-cell cluster's genes say \
                       about its cell type: what the genes do, and why that fits the type.\n\
                       Write a sentence, not a list: name at most three genes. No numbers. Name \
                       only the cell type and the genes below.";

/// The genes of `i`'s top genes that mark `label` in the panel.
fn own_markers(i: &Input, label: &str) -> Vec<String> {
    let here: BTreeSet<&String> = (i.cluster.here.iter())
        .filter(|(l, _)| label_key(l) == label_key(label))
        .flat_map(|(_, m)| m)
        .collect();
    (i.cluster.genes.iter())
        .filter(|g| here.contains(g))
        .cloned()
        .collect()
}

fn biology_prompt(model: LmModel, i: &Input, labels: &[&str]) -> String {
    let mut user = String::new();
    for l in labels {
        user += &format!("CELL TYPE: {}\n", shown(l));
        user += &format!("ITS GENES HERE: {}\n", own_markers(i, l).join(", "));
    }
    model.chat(BIOLOGY, user.trim_end(), "")
}

/// Everything the model is asked about `i`, for the cache key.
#[must_use]
pub fn prompt(model: LmModel, ask: &Ask, i: &Input, panel: &[PanelType]) -> String {
    let mut parts = vec![question(model, ask, i, panel, true)];
    if !i.agree() {
        parts.push(question(model, ask, i, panel, false));
    }
    for l in [&i.call.label, &i.other.label] {
        parts.push(biology_prompt(model, i, &[l]));
    }
    parts.join("\n")
}

/// The cache key of a prompt for a model at a revision, under [`DECODE`]
/// and the words it bans.
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

/// What the model said of the two readings.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Verdict {
    /// They agree, and the genes fit the type.
    Fits,
    /// They agree, and the genes do not fit the type.
    DoesNotFit,
    /// They differ, and the genes favour the call.
    Call,
    /// They differ, and the genes favour the other reading.
    Other,
    /// They differ, and the cluster is a mix of the two or between them.
    Both,
    /// They differ, and neither fits.
    Neither,
}

/// One cluster's critique: a line of the cache.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Critique {
    pub cluster: ClusterId,
    pub call: Reading,
    pub other: Reading,
    pub verdict: Verdict,
    /// Each verdict's probability, both orders averaged.
    pub p: BTreeMap<String, f32>,
    /// The verdict of each order alone, the call named A first.
    #[serde(default)]
    pub by_order: Vec<Verdict>,
    #[serde(default)]
    pub biology: Option<String>,
    pub cited: Vec<String>,
    pub check: Check,
    pub model: String,
    pub revision: String,
    pub key: String,
}

fn method_words(m: Method) -> &'static str {
    match m {
        Method::Enrichment => "the enrichment",
        Method::Susie => "SuSiE",
        Method::Score => "the score",
        Method::Llm => "the model",
    }
}

impl Critique {
    /// The probability of the verdict.
    #[must_use]
    pub fn p_verdict(&self) -> f32 {
        self.p.get(name(self.verdict)).copied().unwrap_or(0.0)
    }

    /// The facts, by code: which method put forward which type.
    #[must_use]
    pub fn facts(&self) -> String {
        let (c, o) = (&self.call, &self.other);
        if label_key(&c.label) == label_key(&o.label) {
            format!(
                "{} and {} agree on {}",
                method_words(c.method),
                method_words(o.method),
                c.label
            )
        } else {
            format!(
                "{} calls {}, {} puts {} first",
                method_words(c.method),
                c.label,
                method_words(o.method),
                o.label
            )
        }
    }

    /// The model's verdict in words.
    #[must_use]
    pub fn verdict_words(&self) -> String {
        match self.verdict {
            Verdict::Fits => format!("the genes fit {}", self.call.label),
            Verdict::DoesNotFit => format!("the genes do not fit {}", self.call.label),
            Verdict::Call => format!("the genes favour {}", self.call.label),
            Verdict::Other => format!("the genes favour {}", self.other.label),
            Verdict::Both => format!("a mix of {} and {}", self.call.label, self.other.label),
            Verdict::Neither => "neither fits the genes".into(),
        }
    }

    /// The strip's text: facts, verdict with its probability, sentence.
    #[must_use]
    pub fn text(&self) -> String {
        let mut s = format!(
            "{} — {} {:.2}",
            self.facts(),
            self.verdict_words(),
            self.p_verdict()
        );
        if let Some(b) = &self.biology {
            s += " — ";
            s += b;
        }
        s
    }

    /// The labels the verdict favours: the sentence's subject.
    fn favoured(&self) -> Vec<&str> {
        match self.verdict {
            Verdict::Fits | Verdict::Call => vec![self.call.label.as_str()],
            Verdict::Other => vec![self.other.label.as_str()],
            Verdict::Both => vec![self.call.label.as_str(), self.other.label.as_str()],
            Verdict::DoesNotFit | Verdict::Neither => Vec::new(),
        }
    }
}

fn name(v: Verdict) -> &'static str {
    match v {
        Verdict::Fits => "fits",
        Verdict::DoesNotFit => "does-not-fit",
        Verdict::Call => "call",
        Verdict::Other => "other",
        Verdict::Both => "both",
        Verdict::Neither => "neither",
    }
}

fn softmax(scores: &[f32]) -> Vec<f32> {
    let top = scores.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let w: Vec<f32> = scores.iter().map(|s| (s - top).exp()).collect();
    let z: f32 = w.iter().sum();
    w.iter().map(|w| w / z).collect()
}

fn best(p: &BTreeMap<Verdict, f32>) -> Verdict {
    let mut it = p.iter();
    let (mut v, mut top) = it
        .next()
        .map(|(v, p)| (*v, *p))
        .unwrap_or((Verdict::Neither, 0.0));
    for (w, &q) in it {
        if q > top {
            (v, top) = (*w, q);
        }
    }
    v
}

/// The verdicts' probabilities from the closed question, each order apart.
fn ask(
    lm: &mut dyn Lm,
    model: LmModel,
    how: &Ask,
    i: &Input,
    panel: &[PanelType],
) -> Result<(BTreeMap<Verdict, f32>, Vec<Verdict>)> {
    if i.agree() {
        let conts: Vec<String> = AGREE.iter().map(|a| format!(" {a}\n")).collect();
        let p = softmax(&lm.score(&question(model, how, i, panel, true), &conts)?);
        let p = BTreeMap::from([(Verdict::Fits, p[0]), (Verdict::DoesNotFit, p[1])]);
        let v = best(&p);
        return Ok((p, vec![v]));
    }
    let conts: Vec<String> = DIFFER.iter().map(|a| format!(" {a}\n")).collect();
    let mut sum: BTreeMap<Verdict, f32> = BTreeMap::new();
    let mut orders = Vec::new();
    for call_first in [true, false] {
        let p = softmax(&lm.score(&question(model, how, i, panel, call_first), &conts)?);
        let (a, b) = if call_first {
            (Verdict::Call, Verdict::Other)
        } else {
            (Verdict::Other, Verdict::Call)
        };
        let one = BTreeMap::from([
            (a, p[0]),
            (b, p[1]),
            (Verdict::Both, p[2]),
            (Verdict::Neither, p[3]),
        ]);
        orders.push(best(&one));
        for (v, q) in one {
            *sum.entry(v).or_default() += q / 2.0;
        }
    }
    Ok((sum, orders))
}

/// The critique of `i` by `lm` (named `model` at `revision`).
pub fn critique(
    lm: &mut dyn Lm,
    spec: &LmSpec,
    how: &Ask,
    i: &Input,
    round: &Round,
) -> Result<Critique> {
    let model = spec.model;
    let (p, by_order) = ask(lm, model, how, i, &round.panel)?;
    let verdict = best(&p);
    let mut c = Critique {
        cluster: i.cluster.cluster,
        call: i.call.clone(),
        other: i.other.clone(),
        verdict,
        p: p.iter().map(|(v, p)| (name(*v).to_string(), *p)).collect(),
        by_order,
        biology: None,
        cited: Vec::new(),
        check: Check::Template,
        model: spec.name(),
        revision: spec.revision(),
        key: key(
            &prompt(model, how, i, &round.panel),
            &spec.name(),
            &spec.revision(),
        ),
    };
    {
        let favoured: Vec<String> = c.favoured().iter().map(|s| s.to_string()).collect();
        let own: BTreeSet<String> = (favoured.iter()).flat_map(|l| own_markers(i, l)).collect();
        if !own.is_empty() {
            let named: BTreeSet<String> = favoured.iter().cloned().collect();
            // Every other top gene and panel marker is another's to name.
            let mut run = round.entities.clone();
            run.extend(i.cluster.genes.iter().cloned());
            run.extend(round.panel.iter().flat_map(|t| t.markers.iter().cloned()));
            let names: Vec<String> = named.iter().chain(&own).map(|n| shown(n)).collect();
            let lower: Vec<String> = names.iter().map(|n| n.to_lowercase()).collect();
            let bans = Bans {
                words: super::vocabulary::bans(&lower),
                numbers: true,
                keep: names,
            };
            let refs: Vec<&str> = favoured.iter().map(String::as_str).collect();
            let lead = opinion::lead(&refs.join(" and "));
            let sentence =
                opinion::checked(lm, &biology_prompt(model, i, &refs), &lead, bans, |t| {
                    opinion::check(t, &run, &named, &own)
                })?
                .map(|s| opinion::first_sentence(&s));
            if let Some(s) = &sentence {
                c.cited = opinion::cited(s, &named, &own);
                c.check = Check::Ok;
            }
            c.biology = sentence;
        }
    }
    Ok(c)
}

/// `{stem}.critique.jsonl`: every critique made on the round at `stem`.
#[must_use]
pub fn cache_path(stem: &str) -> PathBuf {
    PathBuf::from(format!("{stem}.critique.jsonl"))
}

#[must_use]
pub fn cache_of(loaded: &Loaded) -> PathBuf {
    cache_path(&derive_out_prefix(&loaded.file.to_string_lossy()))
}

fn read_cache(path: &Path) -> Vec<Critique> {
    std::fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .filter_map(|l| serde_json::from_str(l).ok())
        .collect()
}

fn append(path: &Path, c: &Critique) -> Result<()> {
    let line = serde_json::to_string(c)? + "\n";
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .with_context(|| format!("opening {}", path.display()))?;
    f.write_all(line.as_bytes())
        .with_context(|| format!("appending to {}", path.display()))
}

const STAGES: crate::progress::Stages = crate::progress::Stages::new(&[("critique", 1.0)]);

/// A critique of every cluster in `inputs`: from `cache` when made on the
/// same inputs, else by the model, loaded at the first miss only.
pub fn critique_round(
    round: &Round,
    inputs: &[Input],
    spec: &LmSpec,
    how: &Ask,
    cache: &Path,
    load: impl FnOnce() -> Result<Box<dyn Lm>>,
) -> Result<Vec<Critique>> {
    let (model, revision) = (spec.name(), spec.revision());
    let done: HashMap<String, Critique> = (read_cache(cache).into_iter())
        .map(|c| (c.key.clone(), c))
        .collect();
    let mut load = Some(load);
    let mut lm: Option<Box<dyn Lm>> = None;
    let n = inputs.len();
    let mut out = Vec::with_capacity(n);
    STAGES.start(0);
    for (k, i) in inputs.iter().enumerate() {
        let key = key(&prompt(spec.model, how, i, &round.panel), &model, &revision);
        let c = match done.get(&key) {
            Some(c) => c.clone(),
            None => {
                if lm.is_none() {
                    lm = Some(load.take().context("the model is loaded once")?()?);
                }
                let c = critique(lm.as_deref_mut().context("a model")?, spec, how, i, round)?;
                append(cache, &c)?;
                c
            }
        };
        out.push(c);
        STAGES.within(0, k + 1, n, Some(&format!("{} of {n} clusters", k + 1)));
    }
    STAGES.finish();
    for c in &out {
        log::info!("K{}: {}", c.cluster, c.text());
    }
    Ok(out)
}
