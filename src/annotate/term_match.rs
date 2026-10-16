//! Panel labels matched to Cell Ontology terms by text embedding, for the
//! labels the matching rules and aliases leave unmatched (abbreviations,
//! compound names: `BFU-E`, `CD8_Effector_Memory`).
//!
//! Every counted name of every live term (names and the synonyms the
//! [`crate::annotate::cl_rules::MatchRules`] count) is encoded once; a label
//! ranks terms by the cosine of its encoding to their nearest name. Only a
//! clear winner is taken: a top score of at least `min_score`, ahead of the
//! runner-up by at least `min_gap`. The model can only pick among the
//! ontology's own terms, never name one.

use crate::annotate::celltype_tree::ClTerms;
use anyhow::Result;
use std::collections::{BTreeMap, HashMap};

/// Texts to L2-normalised vectors, one per text.
pub trait TextEncoder {
    fn encode(&self, texts: &[&str]) -> Result<Vec<Vec<f32>>>;
}

impl TextEncoder for crate::gene_text::encoder::Encoder {
    fn encode(&self, texts: &[&str]) -> Result<Vec<Vec<f32>>> {
        Ok(self
            .encode_batch(texts)?
            .into_iter()
            .map(|d| d.pooled)
            .collect())
    }
}

/// Texts per encoder call; sorted by length, so padding stays small.
const BATCH: usize = 64;

/// Every counted term name, encoded.
pub struct TermIndex {
    ids: Vec<String>,
    vecs: Vec<Vec<f32>>,
}

impl TermIndex {
    pub fn build(terms: &ClTerms, encoder: &impl TextEncoder) -> Result<Self> {
        let mut named: Vec<(&str, &str)> = terms
            .named_texts()
            .flat_map(|(text, ids)| ids.iter().map(move |id| (text, id.as_str())))
            .collect();
        named.sort_by(|a, b| a.0.len().cmp(&b.0.len()).then(a.cmp(b)));
        let mut vecs = Vec::with_capacity(named.len());
        for chunk in named.chunks(BATCH) {
            let texts: Vec<&str> = chunk.iter().map(|(t, _)| *t).collect();
            vecs.extend(encoder.encode(&texts)?);
        }
        let ids = named.into_iter().map(|(_, id)| id.to_string()).collect();
        Ok(Self { ids, vecs })
    }

    /// The `k` terms nearest `v`, each scored by its nearest name, best first.
    #[must_use]
    pub fn rank(&self, v: &[f32], k: usize) -> Vec<(String, f32)> {
        let mut best: HashMap<&str, f32> = HashMap::new();
        for (id, u) in self.ids.iter().zip(&self.vecs) {
            let cos: f32 = u.iter().zip(v).map(|(a, b)| a * b).sum();
            let e = best.entry(id.as_str()).or_insert(f32::NEG_INFINITY);
            *e = e.max(cos);
        }
        let mut ranked: Vec<(String, f32)> = best
            .into_iter()
            .map(|(id, s)| (id.to_string(), s))
            .collect();
        ranked.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        ranked.truncate(k);
        ranked
    }
}

/// When a label's nearest term is taken.
#[derive(Clone, Copy, Debug)]
pub struct Acceptance {
    pub min_score: f32,
    pub min_gap: f32,
}

/// A label's accepted term, with the evidence for it.
#[derive(Clone, Debug, PartialEq)]
pub struct Match {
    pub id: String,
    pub score: f32,
    pub runner_up: f32,
}

/// The labels whose nearest term is a clear winner, each with its term.
pub fn match_labels(
    terms: &ClTerms,
    index: &TermIndex,
    encoder: &impl TextEncoder,
    labels: &[&str],
    accept: Acceptance,
) -> Result<BTreeMap<String, Match>> {
    if labels.is_empty() {
        return Ok(BTreeMap::new());
    }
    let normalised: Vec<String> = labels.iter().map(|l| terms.normalise(l)).collect();
    let texts: Vec<&str> = normalised.iter().map(String::as_str).collect();
    let mut found = BTreeMap::new();
    for (label, v) in labels.iter().zip(encoder.encode(&texts)?) {
        let ranked = index.rank(&v, 2);
        let Some((id, score)) = ranked.first().cloned() else {
            continue;
        };
        let runner_up = ranked.get(1).map_or(f32::NEG_INFINITY, |r| r.1);
        if score >= accept.min_score && score - runner_up >= accept.min_gap {
            found.insert(
                (*label).to_string(),
                Match {
                    id,
                    score,
                    runner_up,
                },
            );
        }
    }
    Ok(found)
}

#[cfg(test)]
#[path = "term_match_tests.rs"]
mod tests;
