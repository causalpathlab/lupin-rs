//! The order view: which cell types precede which, as the prior lupin would
//! build for the labels on screen (`docs/trajectory-plan.md` §6). The user
//! marks two types and states that one precedes the other, or that they are
//! unrelated; the statement goes to the project's `precedence.tsv`.

use crate::annotate::celltype_tree::ClTerms;
use crate::manifest::data_files::{append_line, SearchPath, PRECEDENCE};
use crate::manifest::run::{load, resolve};
use crate::trajectory::prior::{self, Prior, Relation, Source, Statement};
use legume_numeric::matrix::parquet::read_table_columns;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// A pair's data verdict from a `trajectory` run.
pub struct Verdict {
    pub connectivity: f32,
    pub verdict: String,
    pub order_agreement: f32,
}

pub struct OrderView {
    /// The labels on screen with their cell counts.
    pub types: Vec<(String, usize)>,
    /// The combined statements about them.
    pub statements: Vec<Statement>,
    /// The prior they make; `None` when they do not form a DAG (see `error`).
    pub prior: Option<Prior>,
    pub error: Option<String>,
    /// Data verdicts by (from, to) from the latest trajectory run, when any.
    pub verdicts: BTreeMap<(String, String), Verdict>,
    /// The selected type.
    pub sel: usize,
    /// Where a statement is written: the project's `precedence.tsv`.
    pub file: Option<PathBuf>,
}

impl OrderView {
    /// The view for `types`: the Cell Ontology's statements about them and the
    /// user's and project's `precedence.tsv`, combined and reduced, with the
    /// verdicts of a trajectory run recorded in one of `manifests`.
    pub fn load(
        types: Vec<(String, usize)>,
        cl: Option<&ClTerms>,
        search: &SearchPath,
        manifests: &[&Path],
    ) -> Self {
        let names: Vec<Box<str>> = types.iter().map(|(t, _)| t.as_str().into()).collect();
        let mut layers = Vec::new();
        if let Some(cl) = cl {
            let (mapped, _) = cl.map_labels(names.iter().map(AsRef::as_ref));
            layers.push(prior::from_ontology(cl, &mapped));
        }
        let mut problems = Vec::new();
        for (layer, path) in search.user_and_project_files(PRECEDENCE) {
            let source = if layer == "user" {
                Source::User
            } else {
                Source::Project
            };
            let read = std::fs::read_to_string(&path)
                .map_err(|e| format!("reading {}: {e}", path.display()))
                .and_then(|text| {
                    prior::parse_statements(&text, source, &path.display().to_string())
                        .map_err(|e| format!("{e}"))
                });
            match read {
                Ok(st) => layers.push(st),
                Err(e) => problems.push(e),
            }
        }
        let statements = prior::combine(&layers);
        let is_node = vec![true; names.len()];
        let (prior, mut error) = match prior::build(&names, &is_node, statements.clone()) {
            Ok(p) => (Some(p), None),
            Err(e) => (None, Some(e.to_string())),
        };
        if !problems.is_empty() {
            let text = problems.join("; ");
            error = Some(error.map_or(text.clone(), |e| format!("{e}; {text}")));
        }
        let verdicts = manifests
            .iter()
            .find_map(|m| edges_table(m))
            .map(read_verdicts)
            .unwrap_or_default();
        Self {
            types,
            statements,
            prior,
            error,
            verdicts,
            sel: 0,
            file: search.amend(PRECEDENCE),
        }
    }

    pub fn selected(&self) -> Option<&str> {
        self.types.get(self.sel).map(|(t, _)| t.as_str())
    }

    /// The direct edges as `(from, to, source)`, with the verdict when known.
    pub fn edges(&self) -> Vec<(String, String, String, Option<&Verdict>)> {
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
                let v = self.verdict(&from, &to);
                (from, to, source, v)
            })
            .collect()
    }

    /// A pair's verdict from the trajectory run, whichever way it was written.
    pub fn verdict(&self, a: &str, b: &str) -> Option<&Verdict> {
        self.verdicts
            .get(&(a.to_string(), b.to_string()))
            .or_else(|| self.verdicts.get(&(b.to_string(), a.to_string())))
    }

    /// Pairs the trajectory run found connected but the prior does not order.
    pub fn candidates(&self) -> Vec<(&str, &str, &Verdict)> {
        let mut out: Vec<(&str, &str, &Verdict)> = self
            .verdicts
            .iter()
            .filter(|(_, v)| v.verdict == "candidate")
            .map(|((a, b), v)| (a.as_str(), b.as_str(), v))
            .collect();
        out.sort_by(|p, q| q.2.connectivity.total_cmp(&p.2.connectivity));
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
        append_line(
            &file,
            "Which cell types precede which, layered over lupin's own data files\n\
             (see `lupin data where`). Columns: from, to, relation (precedes | unrelated), note.",
            &format!(
                "{from}\t{to}\t{}\t{}",
                relation.as_str(),
                note.replace(['\t', '\n'], " ")
            ),
        )?;
        Ok(file)
    }
}

/// `trajectory.edges` of the manifest at `path`, resolved, when it has one.
fn edges_table(path: &Path) -> Option<String> {
    let loaded = load(&path.to_string_lossy()).ok()?;
    loaded
        .manifest
        .trajectory
        .edges
        .as_deref()
        .map(|rel| resolve(&loaded.dir, rel))
        .filter(|p| Path::new(p).is_file())
}

fn read_verdicts(path: String) -> BTreeMap<(String, String), Verdict> {
    let Ok((strings, numbers)) = read_table_columns(
        &path,
        &["a", "b", "verdict"],
        &["connectivity", "order_agreement"],
    ) else {
        return BTreeMap::new();
    };
    let (Some(a), Some(b), Some(v), Some(c), Some(o)) = (
        strings.first(),
        strings.get(1),
        strings.get(2),
        numbers.first(),
        numbers.get(1),
    ) else {
        return BTreeMap::new();
    };
    (0..a.len())
        .map(|i| {
            (
                (a[i].to_string(), b[i].to_string()),
                Verdict {
                    connectivity: c[i] as f32,
                    verdict: v[i].to_string(),
                    order_agreement: o[i] as f32,
                },
            )
        })
        .collect()
}

#[cfg(test)]
#[path = "tests/order.rs"]
mod tests;
