//! The order view: which cell types precede which, as the prior lupin would
//! build for the labels on screen (`docs/trajectory-plan.md` §6). The user
//! marks two types and states that one precedes the other, or that they are
//! unrelated; the statement goes to the project's `precedence.tsv`.

use crate::annotate::celltype_tree::ClTerms;
use crate::annotate::panel_tree::PanelTree;
use crate::manifest::data_files::{append_row, SearchPath, PRECEDENCE};
use crate::manifest::run::{annotated_path, derive_out_prefix, load, resolve};
use crate::trajectory::edges::{self, EdgeRow, Verdict};
use crate::trajectory::prior::{self, Prior, Relation, Source, Statement};
use std::path::{Path, PathBuf};

/// How the order view runs `lupin trajectory`: the options it carries over
/// and the output prefix last used.
#[derive(Debug, Default)]
pub struct TrajectoryRun {
    pub argv: Vec<String>,
    pub out: Option<String>,
    /// A labels file typed in the TUI for a run with no annotation.
    pub labels: Option<String>,
}

/// The first `{stem}.T{k}` (k = 1, 2, …) beside `source` with no manifest yet.
pub fn default_out(source: &Path) -> String {
    let stem = derive_out_prefix(&source.to_string_lossy());
    (1..)
        .map(|k| format!("{stem}.T{k}"))
        .find(|o| !annotated_path(source, o).exists())
        .expect("some name is free")
}

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
    /// The types on the Cell Ontology, as annotate's tree pane draws them,
    /// and the row selected in it.
    pub tree: PanelTree,
    pub tree_sel: usize,
    /// The full ontology in place of `tree` (`o`).
    pub ontology: Option<super::ontology::OntologyView>,
    /// How many types map to a term, and the ontology's develops-from
    /// links among them.
    pub mapped: usize,
    pub links: usize,
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
        let panel: Vec<(String, String)> = types
            .iter()
            .map(|(t, _)| (String::new(), t.clone()))
            .collect();
        let mapped = cl.map_or(0, |cl| cl.map_labels(refs.iter().copied()).0.len());
        let links = statements
            .iter()
            .filter(|s| {
                matches!(s.source, Source::Cl | Source::ClInherited)
                    && s.relation == Relation::Precedes
            })
            .count();
        Self {
            tree: crate::manifest::ontology::panel_tree_on(cl, &panel),
            tree_sel: 0,
            ontology: None,
            mapped,
            links,
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

    /// Why `from relation to` should not be written: the file's last
    /// statement about the pair already says it, or the prior it makes
    /// would not build (a cycle, or `unrelated` on a pair the rest orders).
    /// On a prior that is already broken any statement goes, so the one that
    /// broke it can be replaced.
    pub fn refuse(&self, from: &str, to: &str, relation: Relation) -> Option<String> {
        let stated = self.file.as_deref().and_then(|f| {
            let text = std::fs::read_to_string(f).ok()?;
            let st = prior::parse_statements(&text, Source::Project, "").ok()?;
            st.into_iter()
                .rev()
                .find(|s| (s.from == from && s.to == to) || (s.from == to && s.to == from))
        });
        if let Some(s) = stated {
            let same_way = relation == Relation::Unrelated || s.from == from;
            if s.relation == relation && same_way {
                return Some(format!(
                    "{} already says {from} {} {to}",
                    self.file.as_deref().unwrap_or(Path::new("")).display(),
                    relation.as_str()
                ));
            }
        }
        self.prior.as_ref()?;
        let new = Statement {
            from: from.into(),
            to: to.into(),
            relation,
            source: Source::Project,
            note: String::new(),
        };
        let statements = prior::combine(&[self.statements.clone(), vec![new]]);
        let names: Vec<Box<str>> = self.types.iter().map(|(t, _)| t.as_str().into()).collect();
        let is_node = vec![true; names.len()];
        prior::build(&names, &is_node, statements)
            .err()
            .map(|e| format!("not recorded: {e}"))
    }

    /// Append a statement to the project's `precedence.tsv`, unless
    /// [`Self::refuse`] gives a reason not to.
    pub fn record(
        &self,
        from: &str,
        to: &str,
        relation: Relation,
        note: &str,
    ) -> anyhow::Result<PathBuf> {
        if let Some(why) = self.refuse(from, to, relation) {
            anyhow::bail!(why);
        }
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
