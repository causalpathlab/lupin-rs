//! The order view: which cell types precede which, as the prior lupin would
//! build for the labels on screen (`docs/trajectory-plan.md` §6). The user
//! marks two types and states that one precedes the other, or that they are
//! unrelated; the statement goes to the project's `precedence.tsv`.

use crate::annotate::celltype_tree::ClTerms;
use crate::manifest::data_files::{append_row, SearchPath, PRECEDENCE};
use crate::manifest::run::{load, resolve};
use crate::trajectory::edges::{self, EdgeRow, Verdict};
use crate::trajectory::prior::{self, Prior, Relation, Statement};
use std::path::{Path, PathBuf};

/// A direct edge of the prior as the table shows it.
pub struct EdgeView<'a> {
    pub from: String,
    pub to: String,
    /// Where the statement came from, or the types it goes through.
    pub source: String,
    /// The trajectory run's row for the pair, when it has one.
    pub data: Option<&'a EdgeRow>,
}

pub struct OrderView {
    /// The labels on screen with their cell counts.
    pub types: Vec<(String, usize)>,
    /// The combined statements about them.
    pub statements: Vec<Statement>,
    /// The prior they make; `None` when they do not form a DAG (see `error`).
    pub prior: Option<Prior>,
    /// Why the statements do not form a DAG.
    pub error: Option<String>,
    /// Precedence files that could not be read.
    pub problems: Vec<String>,
    /// The latest trajectory run's pairs, with their verdicts.
    pub edges: Vec<EdgeRow>,
    /// The selected type.
    pub sel: usize,
    /// Where a statement is written: the project's `precedence.tsv`.
    pub file: Option<PathBuf>,
}

impl OrderView {
    /// The view for `types`: the Cell Ontology's statements about them and the
    /// user's and project's `precedence.tsv`, combined and reduced, with the
    /// `edges` of a trajectory run.
    pub fn load(
        types: Vec<(String, usize)>,
        cl: Option<&ClTerms>,
        search: &SearchPath,
        edges: Vec<EdgeRow>,
    ) -> Self {
        let names: Vec<Box<str>> = types.iter().map(|(t, _)| t.as_str().into()).collect();
        let refs: Vec<&str> = names.iter().map(AsRef::as_ref).collect();
        let (layers, problems) = prior::standing_layers(cl, &refs, search);
        let problems: Vec<String> = problems.iter().map(|e| format!("{e:#}")).collect();
        let statements = prior::combine(&layers);
        let is_node = vec![true; names.len()];
        let (prior, error) = match prior::build(&names, &is_node, statements.clone()) {
            Ok(p) => (Some(p), None),
            Err(e) => (None, Some(e.to_string())),
        };
        Self {
            types,
            statements,
            prior,
            error,
            problems,
            edges,
            sel: 0,
            file: search.amend(PRECEDENCE),
        }
    }

    pub fn selected(&self) -> Option<&str> {
        self.types.get(self.sel).map(|(t, _)| t.as_str())
    }

    /// The prior's direct edges, with the run's row for each when known.
    pub fn edges(&self) -> Vec<EdgeView<'_>> {
        let Some(p) = &self.prior else {
            return Vec::new();
        };
        p.edges
            .iter()
            .map(|e| {
                let from = self.types[e.from].0.clone();
                let to = self.types[e.to].0.clone();
                let source = e.source.map_or_else(
                    || {
                        format!(
                            "via {}",
                            e.via
                                .iter()
                                .map(|&i| self.types[i].0.as_str())
                                .collect::<Vec<_>>()
                                .join(", ")
                        )
                    },
                    |s| s.as_str().to_string(),
                );
                let data = self.row(&from, &to);
                EdgeView {
                    from,
                    to,
                    source,
                    data,
                }
            })
            .collect()
    }

    /// The run's row for a pair, whichever way it was written.
    pub fn row(&self, a: &str, b: &str) -> Option<&EdgeRow> {
        self.edges.iter().find(|e| e.is_pair(a, b))
    }

    /// Pairs the trajectory run found connected but the prior does not
    /// order, strongest first.
    pub fn candidates(&self) -> Vec<&EdgeRow> {
        let mut out: Vec<&EdgeRow> = self
            .edges
            .iter()
            .filter(|e| e.verdict == Some(Verdict::Candidate))
            .collect();
        out.sort_by(|p, q| q.connectivity.total_cmp(&p.connectivity));
        out
    }

    /// The `unrelated` statements, which show no edge.
    pub fn unrelated(&self) -> Vec<&Statement> {
        self.statements
            .iter()
            .filter(|s| s.relation == Relation::Unrelated)
            .collect()
    }

    /// Append a statement to the project's `precedence.tsv`.
    pub fn record(
        &self,
        from: &str,
        to: &str,
        relation: Relation,
        note: &str,
    ) -> anyhow::Result<PathBuf> {
        let file = self.file.clone().ok_or_else(|| {
            anyhow::anyhow!("no project or user directory to write precedence.tsv in")
        })?;
        append_row(
            &file,
            "Which cell types precede which, layered over lupin's own data files\n\
             (see `lupin data where`). Columns: from, to, relation (precedes | unrelated), note.",
            &[from, to, relation.as_str(), note],
        )?;
        Ok(file)
    }
}

/// The edge table of the first of `manifests` whose `trajectory` section
/// records one that can be read; empty otherwise.
pub fn run_edges(manifests: &[&Path]) -> Vec<EdgeRow> {
    manifests
        .iter()
        .find_map(|m| {
            let loaded = load(&m.to_string_lossy()).ok()?;
            let rel = loaded.manifest.trajectory.edges.as_deref()?;
            edges::read(&resolve(&loaded.dir, rel)).ok()
        })
        .unwrap_or_default()
}

#[cfg(test)]
#[path = "tests/order.rs"]
mod tests;
