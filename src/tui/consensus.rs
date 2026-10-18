//! What `lupin describe --lm` writes, in the TUI. In annotate, each
//! cluster's suggested label among its choices with the model's sentence on
//! why, under the candidates (`L` runs it, `P` opens the choices on it). It
//! is a suggestion, not a call: it flags no cluster, nothing is taken until
//! the user confirms it, and the decision then records whether the user
//! took it or overrode it. In the order view, the trajectory's
//! paragraph and each type's line. Both are read from the job's cache as it
//! appends to it, and only what was made on the evidence on screen shows.

use crate::annotate::rounds::ClusterId;
use crate::describe::consensus::{self, cache_path, Proposal};
use crate::describe::trajectory::{self, Described, OrderFacts};
use crate::manifest::run::{derive_out_prefix, load};
use anyhow::Result;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// The size of `file` now; `None` when it is missing.
fn size(file: &Path) -> Option<u64> {
    std::fs::metadata(file).ok().map(|m| m.len())
}

/// A round's proposals, followed in its cache: each cluster's prompt as the
/// round's evidence makes it, so a line made on other evidence is left out.
pub struct Following {
    pub round: PathBuf,
    prompts: BTreeMap<ClusterId, String>,
    pub cache: PathBuf,
    /// The cache's size when last read.
    read: Option<u64>,
}

impl Following {
    pub fn of(round: &Path) -> Result<Self> {
        let loaded = load(&round.to_string_lossy())?;
        let inputs = consensus::inputs(&loaded)?;
        Ok(Self {
            round: round.to_path_buf(),
            prompts: consensus::prompts(&inputs),
            cache: cache_path(&derive_out_prefix(&round.to_string_lossy())),
            read: None,
        })
    }

    /// The proposals now in the cache, when it changed since last read.
    pub fn poll(&mut self) -> Option<BTreeMap<ClusterId, Proposal>> {
        let now = size(&self.cache);
        if now == self.read {
            return None;
        }
        self.read = now;
        Some(consensus::current(&self.prompts, &self.cache))
    }
}

/// A trajectory's texts, followed in its cache: the templates where it has
/// none yet.
pub struct OrderTexts {
    pub manifest: PathBuf,
    facts: OrderFacts,
    cache: PathBuf,
    read: Option<u64>,
    pub texts: Vec<Described>,
}

impl OrderTexts {
    pub fn of(manifest: &Path) -> Result<Self> {
        let loaded = load(&manifest.to_string_lossy())?;
        let facts = OrderFacts::read(&loaded)?;
        let cache = cache_path(&derive_out_prefix(&manifest.to_string_lossy()));
        Ok(Self {
            manifest: manifest.to_path_buf(),
            texts: trajectory::current(&facts, &cache),
            read: size(&cache),
            facts,
            cache,
        })
    }

    /// Read the cache again when it changed; whether it did.
    pub fn poll(&mut self) -> bool {
        let now = size(&self.cache);
        if now == self.read {
            return false;
        }
        self.read = now;
        self.texts = trajectory::current(&self.facts, &self.cache);
        true
    }

    /// Type `ty`'s line, or the run's paragraph without one.
    #[must_use]
    pub fn line(&self, ty: Option<&str>) -> Option<&str> {
        let of = |s: Option<&str>| {
            self.texts
                .iter()
                .find(|d| d.subject.as_deref() == s)
                .map(|d| d.text.as_str())
        };
        ty.and_then(|t| of(Some(t))).or_else(|| of(None))
    }
}

/// The strip under the candidates: the suggestion, its probability and why.
#[must_use]
pub fn strip(p: &Proposal) -> String {
    format!(
        "◆ suggested {} {:.2} — {}",
        p.label,
        p.p_label(),
        p.rationale
    )
}

#[cfg(test)]
#[path = "tests/consensus.rs"]
mod tests;
