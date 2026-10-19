//! What the TUI exported: a log of each figure set (a `.pdf`, with an `.svg`
//! beside it but for a cell map, sharing a base name) and a thumbnail of
//! it, in `.lupin-view/` in the directory lupin runs from, kept by
//! `legume_tui::gallery` under its traced policy: every file recorded with
//! its size and hash and checked against the disk, nothing dropped, an
//! unreadable log set aside, and files moved or deleted only on the word of
//! a private log. What lupin keeps beside each entry to trace a figure back
//! is its [`Trace`] (`docs/trajectory-plan.md` §6). The runs' own
//! directories are not touched. The logs lupin wrote before load as they
//! are.

use crate::manifest::data_files::absolute;
use legume_tui::gallery::Policy;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

pub use legume_tui::gallery::{ago, now, Status};

/// Where the log lives, under the working directory.
pub const DIR: &str = ".lupin-view";

/// What a figure was drawn from, kept in its entry.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct Trace {
    /// The figure (`layout_umap`, `diffusion_dc1_dc2`, …).
    #[serde(default)]
    pub panel: String,
    /// The manifest the figure was drawn from, absolute.
    #[serde(default)]
    pub manifest: PathBuf,
}

impl Trace {
    pub fn of(panel: &str, manifest: &Path) -> Self {
        Self {
            panel: panel.into(),
            manifest: absolute(manifest),
        }
    }
}

/// The exports logged in one directory, newest first.
pub type Gallery = legume_tui::gallery::Gallery<Trace>;

/// The log in `dir` (`.lupin-view`).
pub fn open(dir: &Path) -> Gallery {
    Gallery::open(dir, Policy::traced())
}

/// The log under the working directory.
pub fn here() -> Gallery {
    open(Path::new(DIR))
}

/// Figures named `{figures}.*.pdf` (`figures` is the run's prefix and the
/// figures' kind, `run.trajectory`) that no entry of `g` lists.
pub fn not_listed(g: &Gallery, figures: &str) -> Vec<PathBuf> {
    let Some(dir) = Path::new(figures).parent() else {
        return Vec::new();
    };
    let stem = figure_stem(figures);
    let listed: std::collections::BTreeSet<&Path> =
        g.entries().iter().map(|e| e.path.as_path()).collect();
    let Ok(read) = std::fs::read_dir(if dir.as_os_str().is_empty() {
        Path::new(".")
    } else {
        dir
    }) else {
        return Vec::new();
    };
    let mut out: Vec<PathBuf> = read
        .filter_map(Result::ok)
        .map(|d| d.path())
        .filter(|p| {
            p.extension().is_some_and(|e| e == "pdf")
                && p.file_name()
                    .is_some_and(|n| n.to_string_lossy().starts_with(&stem))
                && !listed.contains(absolute(p).as_path())
        })
        .collect();
    out.sort();
    out
}

/// The start of a figure file's name for `figures` (the run's prefix and
/// the figures' kind, `dir/run.trajectory`): `run.trajectory.`.
pub fn figure_stem(figures: &str) -> String {
    Path::new(figures)
        .file_name()
        .map_or(String::new(), |s| format!("{}.", s.to_string_lossy()))
}

#[cfg(test)]
#[path = "tests/gallery.rs"]
mod tests;
