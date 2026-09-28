//! Genes kept out of the specific-genes view: exact names or `*` patterns
//! (`MT-*`, `RP[SL]` as `RPL*` and `RPS*`), one per line in
//! `hidden_genes.txt`, read from the user's config and the project's `lupin/`
//! folder (see [`crate::manifest::data_files`]) and added to from the TUI.
//! Data, not code: nothing is hidden unless a file says so.

use crate::manifest::data_files::{append_line, SearchPath};
use anyhow::Result;
use std::path::{Path, PathBuf};

/// The file the hidden genes are kept in.
pub const HIDDEN: &str = "hidden_genes.txt";

#[derive(Debug, Clone, Default, PartialEq)]
pub struct GeneFilter {
    patterns: Vec<String>,
}

impl GeneFilter {
    /// The patterns in `text`, one per line; `#` comments and blanks skipped.
    #[cfg(test)]
    #[must_use]
    pub fn parse(text: &str) -> Self {
        let mut f = Self::default();
        f.extend(text);
        f
    }

    fn extend(&mut self, text: &str) {
        self.patterns.extend(
            text.lines()
                .map(str::trim)
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
                .map(String::from),
        );
    }

    /// The user's and the project's hidden genes, both applied.
    pub fn load(search: &SearchPath) -> Result<Self> {
        let mut f = Self::default();
        for text in search.user_and_project(HIDDEN)? {
            f.extend(&text);
        }
        Ok(f)
    }

    /// Where the TUI adds patterns: the project's file, else the user's.
    #[must_use]
    pub fn file(search: &SearchPath) -> Option<PathBuf> {
        search.amend(HIDDEN)
    }

    /// Hide `pattern` from now on, and append it to `file`.
    pub fn add(&mut self, pattern: &str, file: &Path) -> Result<()> {
        let pattern = pattern.trim();
        if pattern.is_empty() || self.patterns.iter().any(|p| p == pattern) {
            return Ok(());
        }
        append_line(
            file,
            "Genes lupin's TUI keeps out of the specific-genes view: names or `*` patterns.",
            pattern,
        )?;
        self.patterns.push(pattern.to_string());
        Ok(())
    }

    #[must_use]
    pub fn hides(&self, gene: &str) -> bool {
        self.patterns.iter().any(|p| glob(p, gene))
    }

    #[cfg(test)]
    #[must_use]
    pub fn len(&self) -> usize {
        self.patterns.len()
    }
}

/// A pattern to hide genes like `gene`: its stem before the first `-` or
/// `.` with a `*`, else the name itself (`MT-ND3` → `MT-*`).
#[must_use]
pub fn suggest(gene: &str) -> String {
    match gene.find(['-', '.']) {
        Some(i) if i > 0 => format!("{}*", &gene[..=i]),
        _ => gene.to_string(),
    }
}

/// Whether `name` matches `pattern`, where `*` stands for any run of
/// characters; letters compare as written (gene symbols are case-bearing).
fn glob(pattern: &str, name: &str) -> bool {
    let parts: Vec<&str> = pattern.split('*').collect();
    if parts.len() == 1 {
        return pattern == name;
    }
    let (first, last) = (parts[0], parts[parts.len() - 1]);
    if !name.starts_with(first) || !name[first.len()..].ends_with(last) {
        return false;
    }
    let mut rest = &name[first.len()..name.len() - last.len()];
    for mid in &parts[1..parts.len() - 1] {
        match rest.find(mid) {
            Some(i) => rest = &rest[i + mid.len()..],
            None => return false,
        }
    }
    true
}

#[cfg(test)]
#[path = "tests/genes.rs"]
mod tests;
