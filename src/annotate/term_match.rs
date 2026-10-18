//! Cell Ontology terms suggested by text embedding for the panel labels the
//! matching rules and aliases leave unmatched (abbreviations, compound
//! names: `BFU-E`, `CD8_Effector_Memory`).
//!
//! Every counted name of every live term (names and the synonyms the
//! [`crate::annotate::cl_rules::MatchRules`] count) is encoded once; a label
//! ranks terms by the cosine of its encoding to their nearest name. The
//! model can only pick among the ontology's own terms, never name one, and
//! nothing is taken on its own: no score threshold told right from wrong
//! well enough, so a label's top terms are offered with their scores and the
//! user confirms one into their alias table (the annotate TUI's `c`).
//!
//! The encoder is PubMedBERT tuned as a sentence embedding (NeuML's
//! `pubmedbert-base-embeddings`, at a pinned revision), mean-pooled as
//! sentence-transformers pools it:
//! - Gu Y, Tinn R, Cheng H, Lucas M, Usuyama N, Liu X, Naumann T, Gao J,
//!   Poon H (2021). Domain-specific language model pretraining for
//!   biomedical natural language processing. ACM Transactions on Computing
//!   for Healthcare 3(1):2. doi:10.1145/3458754
//! - Reimers N, Gurevych I (2019). Sentence-BERT: sentence embeddings using
//!   Siamese BERT-networks. EMNLP-IJCNLP 2019. doi:10.18653/v1/D19-1410
//!
//! Encoding the ontology takes a while on a CPU, so the index is cached
//! under lupin's cache directory, keyed by the ontology's release and the
//! model's revision ([`TermIndex::cache_name`]).

use crate::annotate::celltype_tree::ClTerms;
use crate::gene_text::encoder::{Encoder, ModelSpec, Pooling};
use crate::manifest::data_files::{ClData, Fetch, SearchPath};
use crate::progress::Stages;
use anyhow::{Context, Result};
use legume_numeric::matrix::traits::IoOps;
use log::{info, warn};
use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::path::Path;

/// The text model suggestions are ranked by, and its pinned revision.
pub const MODEL: &str = "NeuML/pubmedbert-base-embeddings";
pub const REVISION: &str = "b79526d6ef3645e0df4530322e266f24c829f5ef";
/// Terms suggested per label.
pub const SUGGESTED: usize = 3;
/// Labels and term names are short.
const MAX_TOKENS: usize = 64;

/// Texts to L2-normalised vectors, one per text.
pub trait TextEncoder {
    fn encode(&self, texts: &[&str]) -> Result<Vec<Vec<f32>>>;
}

impl TextEncoder for Encoder {
    fn encode(&self, texts: &[&str]) -> Result<Vec<Vec<f32>>> {
        self.encode_pooled(texts)
    }
}

/// Texts per encoder call; sorted by length, so padding stays small.
const BATCH: usize = 64;

/// Every counted term name, encoded, the names of a term together.
pub struct TermIndex {
    /// Each term, by id.
    terms: Vec<String>,
    /// Where each term's names end in [`Self::vecs`].
    ends: Vec<usize>,
    vecs: Vec<Vec<f32>>,
}

impl TermIndex {
    /// Every counted name of `terms` with its term, shortest first: what
    /// [`Self::build`] encodes and [`Self::cache_name`] digests.
    #[must_use]
    pub fn named(terms: &ClTerms) -> Vec<(&str, &str)> {
        let mut named: Vec<(&str, &str)> = terms
            .named_texts()
            .flat_map(|(text, ids)| ids.iter().map(move |id| (text, id.as_str())))
            .collect();
        named.sort_by(|a, b| a.0.len().cmp(&b.0.len()).then(a.cmp(b)));
        named
    }

    /// Encode every counted name ([`Self::named`]); `progress` hears the
    /// names encoded so far of all, after each batch.
    pub fn build(
        named: &[(&str, &str)],
        encoder: &impl TextEncoder,
        mut progress: impl FnMut(usize, usize),
    ) -> Result<Self> {
        let mut vecs = Vec::with_capacity(named.len());
        for chunk in named.chunks(BATCH) {
            let texts: Vec<&str> = chunk.iter().map(|(t, _)| *t).collect();
            vecs.extend(encoder.encode(&texts)?);
            progress(vecs.len(), named.len());
        }
        let ids = named.iter().map(|(_, id)| id.to_string()).collect();
        Ok(Self::of_rows(ids, vecs))
    }

    /// One name's vector per row, `ids[i]` its term: the names of a term
    /// put together, terms by id.
    fn of_rows(ids: Vec<String>, vecs: Vec<Vec<f32>>) -> Self {
        let mut rows: Vec<(String, Vec<f32>)> = ids.into_iter().zip(vecs).collect();
        rows.sort_by(|a, b| a.0.cmp(&b.0));
        let (mut terms, mut ends, mut vecs) = (Vec::new(), Vec::new(), Vec::new());
        for (id, v) in rows {
            if terms.last() != Some(&id) {
                if !vecs.is_empty() {
                    ends.push(vecs.len());
                }
                terms.push(id);
            }
            vecs.push(v);
        }
        if !vecs.is_empty() {
            ends.push(vecs.len());
        }
        Self { terms, ends, vecs }
    }

    /// Each term's names' rows in [`Self::vecs`], with the term.
    fn spans(&self) -> impl Iterator<Item = (&str, std::ops::Range<usize>)> {
        let starts = std::iter::once(0).chain(self.ends.iter().copied());
        (self.terms.iter().zip(starts.zip(&self.ends))).map(|(t, (a, &b))| (t.as_str(), a..b))
    }

    /// The cache file's name for the index of an ontology at `release`
    /// (its counted names `named`) by `model` at `revision`: the release and
    /// the revision, then a digest of both with every counted name, so other
    /// matching rules (which count other synonyms) make another index.
    #[must_use]
    pub fn cache_name(
        release: Option<&str>,
        named: &[(&str, &str)],
        model: &str,
        revision: &str,
    ) -> String {
        let mut h = Sha256::new();
        for part in [model, revision] {
            h.update(part);
            h.update([0]);
        }
        for (text, id) in named {
            h.update(text);
            h.update([0]);
            h.update(id);
            h.update([0]);
        }
        let digest = crate::hex(&h.finalize()[..8]);
        let safe = |s: &str| -> String {
            s.split(|c: char| !(c.is_ascii_alphanumeric() || c == '-' || c == '.'))
                .filter(|w| !w.is_empty())
                .collect::<Vec<_>>()
                .join("_")
        };
        let release = release.map_or("unreleased".into(), safe);
        let rev: String = safe(revision).chars().take(12).collect();
        format!("{release}-{rev}-{digest}.parquet")
    }

    /// Write the index to `file` (a parquet table: the term, then the
    /// vector), whole or not at all: through a file of its own beside it,
    /// so another process writing the same index never shares it.
    pub fn save(&self, file: &Path) -> Result<()> {
        let dir = file.parent().filter(|d| !d.as_os_str().is_empty());
        let dir = dir.unwrap_or_else(|| Path::new("."));
        std::fs::create_dir_all(dir)?;
        let dim = self.vecs.first().map_or(0, Vec::len);
        let m =
            DMatrix::from_row_iterator(self.vecs.len(), dim, self.vecs.iter().flatten().copied());
        let rows: Vec<Box<str>> = (self.spans())
            .flat_map(|(t, r)| std::iter::repeat_n(Box::<str>::from(t), r.len()))
            .collect();
        let partial = tempfile::Builder::new()
            .suffix(".partial")
            .tempfile_in(dir)
            .with_context(|| format!("a file to write {} through", file.display()))?;
        m.to_parquet_with_names(
            &partial.path().to_string_lossy(),
            (Some(&rows), Some("term")),
            None,
        )?;
        partial
            .persist(file)
            .with_context(|| format!("writing {}", file.display()))?;
        Ok(())
    }

    /// The index [`Self::save`] wrote to `file`.
    pub fn load(file: &Path) -> Result<Self> {
        let m = DMatrix::<f32>::from_parquet(&file.to_string_lossy())
            .with_context(|| format!("reading {}", file.display()))?;
        Ok(Self::of_rows(
            m.rows.into_iter().map(String::from).collect(),
            m.mat
                .row_iter()
                .map(|r| r.iter().copied().collect())
                .collect(),
        ))
    }

    /// The `k` terms nearest `v`, each scored by its nearest name, best
    /// first, ties by id.
    #[must_use]
    pub fn rank(&self, v: &[f32], k: usize) -> Vec<(String, f32)> {
        let cos = |u: &Vec<f32>| -> f32 { u.iter().zip(v).map(|(a, b)| a * b).sum() };
        // `max` passes over a NaN cosine, so no score is NaN.
        let scores = (self.spans().enumerate()).map(|(t, (_, r))| {
            (
                t,
                self.vecs[r]
                    .iter()
                    .map(cos)
                    .fold(f32::NEG_INFINITY, f32::max),
            )
        });
        crate::gene_text::edges::top_k(scores, k)
            .into_iter()
            .map(|(t, s)| (self.terms[t].clone(), s))
            .collect()
    }
}

/// A term suggested for a label: its id, its name, and the cosine of the
/// label to the term's nearest name.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Suggestion {
    pub id: String,
    pub name: String,
    pub score: f32,
}

/// Each label's `k` nearest terms, best first; none is taken.
pub fn suggest(
    terms: &ClTerms,
    index: &TermIndex,
    encoder: &impl TextEncoder,
    labels: &[&str],
    k: usize,
) -> Result<BTreeMap<String, Vec<Suggestion>>> {
    if labels.is_empty() {
        return Ok(BTreeMap::new());
    }
    let normalised: Vec<String> = labels.iter().map(|l| terms.normalise(l)).collect();
    let texts: Vec<&str> = normalised.iter().map(String::as_str).collect();
    let mut found = BTreeMap::new();
    for (label, v) in labels.iter().zip(encoder.encode(&texts)?) {
        let ranked = index
            .rank(&v, k)
            .into_iter()
            .map(|(id, score)| Suggestion {
                name: terms.name(&id).unwrap_or(&id).to_string(),
                id,
                score,
            })
            .collect();
        found.insert((*label).to_string(), ranked);
    }
    Ok(found)
}

/// What `lupin data suggest-cl` writes: each label's nearest terms, and what
/// ranked them.
#[derive(Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Suggested {
    pub model: String,
    pub revision: String,
    pub ontology_release: Option<String>,
    pub labels: BTreeMap<String, Vec<Suggestion>>,
}

static STAGES: Stages = Stages::new(&[
    ("loading the text model", 5.0),
    ("encoding the ontology's names (once; cached)", 90.0),
    ("ranking the labels", 5.0),
]);

/// The pinned model, from the Hub's cache, else the Hub unless offline
/// ([`crate::hub`]).
fn load_encoder() -> Result<Encoder> {
    let spec = ModelSpec {
        model_id: MODEL.into(),
        revision: REVISION.into(),
        local_dir: None,
    };
    Encoder::load(&spec, MAX_TOKENS, Pooling::Mean, candle_core::Device::Cpu)
        .with_context(|| format!("loading the text model {MODEL}"))
}

/// `lupin data suggest-cl`: the nearest terms of each of `labels` on the
/// ontology a run in `run_dir` finds (`obo` first), as JSON to `out` (else
/// stdout). The ontology's index is read from the cache, else built there.
pub fn run(
    run_dir: Option<&Path>,
    obo: Option<&str>,
    out: Option<&str>,
    labels: &[String],
) -> Result<()> {
    let data = ClData::load(SearchPath::new(run_dir), obo, None, Fetch::Allowed)?;
    let cache = data.search.cache.clone();
    let terms = data.into_terms()?.context(
        "no Cell Ontology at hand (see `lupin data where`): nothing to suggest terms from",
    )?;
    STAGES.start(0);
    let encoder = load_encoder()?;
    let named = TermIndex::named(&terms);
    let file = cache.map(|c| {
        c.join("cl-index").join(TermIndex::cache_name(
            terms.release.as_deref(),
            &named,
            MODEL,
            REVISION,
        ))
    });
    let cached = file
        .as_deref()
        .filter(|f| f.is_file())
        .map(|f| TermIndex::load(f).inspect_err(|e| warn!("{e:#}; building the index again")));
    let index = match cached {
        Some(Ok(index)) => index,
        _ => {
            let stage = 1;
            STAGES.start(stage);
            let index = TermIndex::build(&named, &encoder, |k, n| {
                STAGES.within(stage, k, n, Some(&format!("{k} of {n}")));
            })?;
            if let Some(f) = &file {
                match index.save(f) {
                    Ok(()) => info!("Cell Ontology index cached at {}", f.display()),
                    Err(e) => warn!("could not cache the index: {e:#}"),
                }
            }
            index
        }
    };
    STAGES.start(2);
    let labels: Vec<&str> = labels.iter().map(String::as_str).collect();
    let found = Suggested {
        model: MODEL.into(),
        revision: REVISION.into(),
        ontology_release: terms.release.clone(),
        labels: suggest(&terms, &index, &encoder, &labels, SUGGESTED)?,
    };
    let json = serde_json::to_string_pretty(&found)?;
    match out {
        Some(f) => std::fs::write(f, json + "\n").with_context(|| format!("writing {f}"))?,
        None => println!("{json}"),
    }
    STAGES.finish();
    Ok(())
}

#[cfg(test)]
#[path = "term_match_tests.rs"]
mod tests;
