//! What `lupin describe` gives the TUI. In annotate, the language model's
//! opinion of each cluster, a method beside the others: its probabilities
//! join the candidates' evidence, and its top is suggested under the
//! candidates (`L` runs it, `P` opens the choices on it) with the facts on
//! that type by code and the model's sentence on its biology. It is a
//! suggestion, not a call: it flags no cluster, nothing is taken until the
//! user confirms it, and the decision then records whether the user took
//! it or overrode it. The opinions are read from the round's cache as the
//! job appends to it, and only those made on the genes and terms on screen
//! show. In the order view, the trajectory's paragraph and each type's
//! line, written by code from the run's outputs and read again when they
//! change.

use crate::annotate::rounds::ClusterId;
use crate::describe::opinion::{self, Opinion};
use crate::describe::trajectory::{self, OrderFacts, OrderText};
use crate::manifest::run::load;
use crate::tui::round::Candidate;
use anyhow::Result;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// A file's size and last change, to tell when it changed; `None` when it
/// is missing.
type Stamp = Option<(u64, Option<std::time::SystemTime>)>;

fn stamp(file: &Path) -> Stamp {
    let m = std::fs::metadata(file).ok()?;
    Some((m.len(), m.modified().ok()))
}

/// A round's opinions, followed in its cache: each cluster's prompt as the
/// round's genes and terms make it, so a line made on others is left out.
pub struct Following {
    pub round: PathBuf,
    prompts: BTreeMap<ClusterId, String>,
    pub cache: PathBuf,
    /// The cache's [`Stamp`] when last read.
    read: Stamp,
}

impl Following {
    pub fn of(round: &Path) -> Result<Self> {
        let loaded = load(&round.to_string_lossy())?;
        let inputs = opinion::inputs(&loaded)?;
        Ok(Self {
            round: round.to_path_buf(),
            prompts: opinion::prompts(&inputs),
            cache: opinion::cache_of(&loaded),
            read: None,
        })
    }

    /// The opinions now in the cache, when it changed since last read.
    pub fn poll(&mut self) -> Option<BTreeMap<ClusterId, Opinion>> {
        let now = stamp(&self.cache);
        if now == self.read {
            return None;
        }
        self.read = now;
        Some(opinion::current(&self.prompts, &self.cache))
    }
}

/// A trajectory's texts, from its facts: read again when the manifest or a
/// file of its facts changes.
pub struct OrderTexts {
    pub manifest: PathBuf,
    /// The manifest and its facts' files, as last read.
    files: Vec<PathBuf>,
    read: Vec<Stamp>,
    pub texts: Vec<OrderText>,
}

impl OrderTexts {
    pub fn of(manifest: &Path) -> Result<Self> {
        let loaded = load(&manifest.to_string_lossy())?;
        let files: Vec<PathBuf> = std::iter::once(manifest.to_path_buf())
            .chain(
                OrderFacts::files(&loaded)
                    .into_iter()
                    .flatten()
                    .map(PathBuf::from),
            )
            .collect();
        Ok(Self {
            manifest: manifest.to_path_buf(),
            read: files.iter().map(|f| stamp(f)).collect(),
            files,
            texts: trajectory::texts(&OrderFacts::read(&loaded)?),
        })
    }

    /// Read the texts again when their files changed; whether they did. A
    /// manifest or facts that no longer read keep the texts as they were.
    pub fn poll(&mut self) -> bool {
        let now: Vec<Stamp> = self.files.iter().map(|f| stamp(f)).collect();
        if now == self.read {
            return false;
        }
        match Self::of(&self.manifest) {
            Ok(again) => *self = again,
            Err(_) => self.read = now,
        }
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

/// The strip under the candidates: the model's suggestion, its
/// probability, the facts on it and its sentence, on a cluster of
/// `candidates`.
#[must_use]
pub fn strip(o: &Opinion, candidates: &[Candidate]) -> String {
    format!(
        "◆ suggested {} {:.2} — {}",
        o.label,
        o.p_label(),
        o.text(candidates)
    )
}

#[cfg(test)]
#[path = "tests/describe.rs"]
mod tests;
