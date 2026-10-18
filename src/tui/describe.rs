//! What `lupin describe` gives the TUI: in the order view, the
//! trajectory's paragraph and each type's line, written by code from the
//! run's outputs and read again when they change.

use crate::describe::trajectory::{self, OrderFacts, OrderText};
use crate::manifest::run::load;
use anyhow::Result;
use std::path::{Path, PathBuf};

/// A file's size and last change, to tell when it changed; `None` when it
/// is missing.
type Stamp = Option<(u64, Option<std::time::SystemTime>)>;

fn stamp(file: &Path) -> Stamp {
    let m = std::fs::metadata(file).ok()?;
    Some((m.len(), m.modified().ok()))
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

#[cfg(test)]
#[path = "tests/describe.rs"]
mod tests;
