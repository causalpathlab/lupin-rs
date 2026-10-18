//! A trajectory run's order in words: one paragraph on the run (its
//! components, roots, the types where its lineages part and where they
//! end) and one line per type (its lineages, the types either side of it,
//! the verdicts of its edges). Both are written by code from the run's own
//! outputs (`trajectory_lineages.tsv` and `trajectory_edges.parquet`), so
//! they say only what the run found; no model words them.

use crate::manifest::run::{resolve, Loaded};
use crate::trajectory::edges::{self, EdgeRow};
use anyhow::{Context, Result};
use std::collections::BTreeSet;
use std::path::PathBuf;

/// One root-to-leaf path.
#[derive(Clone, Debug, PartialEq)]
pub struct Lineage {
    pub name: String,
    pub component: String,
    pub path: Vec<String>,
}

/// What a trajectory run found, as its outputs record it.
#[derive(Clone, Debug, Default)]
pub struct OrderFacts {
    pub lineages: Vec<Lineage>,
    pub edges: Vec<EdgeRow>,
}

/// The lineages of `trajectory_lineages.tsv`'s `text`.
#[must_use]
pub fn parse_lineages(text: &str) -> Vec<Lineage> {
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let mut it = l.split('\t');
            let (name, component, path) = (it.next()?, it.next()?, it.next()?);
            Some(Lineage {
                name: name.to_string(),
                component: component.to_string(),
                path: path.split(" > ").map(str::to_string).collect(),
            })
        })
        .collect()
}

impl OrderFacts {
    /// The facts of the trajectory run `loaded` records.
    pub fn read(loaded: &Loaded) -> Result<Self> {
        let [lineages, edges] = Self::files(loaded);
        let edges = match edges {
            Some(p) => edges::read(&p)?,
            None => Vec::new(),
        };
        let lineages = match lineages {
            Some(p) => parse_lineages(
                &std::fs::read_to_string(&p).with_context(|| format!("reading {p}"))?,
            ),
            None => Vec::new(),
        };
        Ok(Self { lineages, edges })
    }

    /// The files [`Self::read`] reads: the lineages and the edges, where
    /// `loaded` records them.
    #[must_use]
    pub fn files(loaded: &Loaded) -> [Option<String>; 2] {
        let t = &loaded.manifest.trajectory;
        [&t.lineages, &t.edges].map(|rel| rel.as_deref().map(|r| resolve(&loaded.dir, r)))
    }

    /// Where the lineages start, in order.
    #[must_use]
    pub fn roots(&self) -> Vec<String> {
        unique(self.lineages.iter().filter_map(|l| l.path.first().cloned()))
    }

    /// Where they end, in order.
    #[must_use]
    pub fn ends(&self) -> Vec<String> {
        unique(self.lineages.iter().filter_map(|l| l.path.last().cloned()))
    }

    /// The types on more than one lineage of a component where two of
    /// them part: the last type their paths share.
    #[must_use]
    pub fn branch_points(&self) -> Vec<String> {
        let mut out = Vec::new();
        for (i, a) in self.lineages.iter().enumerate() {
            for b in &self.lineages[i + 1..] {
                if a.component != b.component {
                    continue;
                }
                let shared = a
                    .path
                    .iter()
                    .zip(&b.path)
                    .take_while(|(x, y)| x == y)
                    .count();
                if shared > 0 && shared < a.path.len() && shared < b.path.len() {
                    out.push(a.path[shared - 1].clone());
                }
            }
        }
        unique(out.into_iter())
    }

    /// How many components the lineages span.
    #[must_use]
    pub fn components(&self) -> usize {
        self.lineages
            .iter()
            .map(|l| l.component.as_str())
            .collect::<BTreeSet<_>>()
            .len()
    }

    /// Every type on a lineage, then every other one an edge names, in
    /// order.
    #[must_use]
    pub fn types(&self) -> Vec<String> {
        let on = self.lineages.iter().flat_map(|l| l.path.iter().cloned());
        let edged = self
            .edges
            .iter()
            .filter(|e| e.verdict.is_some())
            .flat_map(|e| [e.a.to_string(), e.b.to_string()]);
        unique(on.chain(edged))
    }

    fn lineages_of(&self, t: &str) -> Vec<&Lineage> {
        self.lineages
            .iter()
            .filter(|l| l.path.iter().any(|p| p == t))
            .collect()
    }

    /// The edges of `t` with a verdict.
    fn edges_of(&self, t: &str) -> Vec<&EdgeRow> {
        self.edges
            .iter()
            .filter(|e| e.verdict.is_some() && (e.a.as_ref() == t || e.b.as_ref() == t))
            .collect()
    }

    /// The types either side of `t` on its lineages, before and after.
    fn sides(&self, t: &str) -> (Vec<String>, Vec<String>) {
        let mut before = Vec::new();
        let mut after = Vec::new();
        for l in self.lineages_of(t) {
            let i = l.path.iter().position(|p| p == t).unwrap_or(0);
            before.extend(i.checked_sub(1).map(|j| l.path[j].clone()));
            after.extend(l.path.get(i + 1).cloned());
        }
        (unique(before.into_iter()), unique(after.into_iter()))
    }
}

/// `items` without repeats, in order.
fn unique(items: impl Iterator<Item = String>) -> Vec<String> {
    let mut seen = BTreeSet::new();
    items.filter(|t| seen.insert(t.clone())).collect()
}

/// `a`, `a and b`, `a, b and c`.
fn and_list(items: &[String]) -> String {
    match items {
        [] => String::new(),
        [a] => a.clone(),
        [rest @ .., last] => format!("{} and {last}", rest.join(", ")),
    }
}

fn edge_text(e: &EdgeRow) -> String {
    let verdict = e.verdict.map_or("", |v| v.as_str());
    format!("{} → {} {verdict}", e.a, e.b)
}

/// The run paragraph when no model words it.
#[must_use]
pub fn run_template(f: &OrderFacts) -> String {
    if f.lineages.is_empty() {
        return "The run ordered no lineage: no root reaches its types.".into();
    }
    let mut s = format!(
        "{} lineage(s) in {} component(s). The order starts at {}.",
        f.lineages.len(),
        f.components(),
        and_list(&f.roots())
    );
    let bps = f.branch_points();
    if !bps.is_empty() {
        s += &format!(" It branches after {}.", and_list(&bps));
    }
    s += &format!(" Its lineages end at {}.", and_list(&f.ends()));
    let paths: Vec<String> = f
        .lineages
        .iter()
        .map(|l| format!("{}: {}", l.name, l.path.join(" > ")))
        .collect();
    format!("{s} {}.", paths.join("; "))
}

/// Type `t`'s line when no model words it.
#[must_use]
pub fn type_template(f: &OrderFacts, t: &str) -> String {
    let on: Vec<String> = f.lineages_of(t).iter().map(|l| l.name.clone()).collect();
    let mut s = if on.is_empty() {
        format!("{t} is on no lineage.")
    } else {
        let (before, after) = f.sides(t);
        let mut place = Vec::new();
        if !before.is_empty() {
            place.push(format!("after {}", and_list(&before)));
        }
        if !after.is_empty() {
            place.push(format!("before {}", and_list(&after)));
        }
        let place = if place.is_empty() {
            String::new()
        } else {
            format!(", {}", place.join(" and "))
        };
        format!("{t} is on {}{place}.", and_list(&on))
    };
    let edges: Vec<String> = f.edges_of(t).into_iter().map(edge_text).collect();
    if !edges.is_empty() {
        s += &format!(" Its edges: {}.", edges.join("; "));
    }
    s
}

/// One text: the run's paragraph (`subject` none) or a type's line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OrderText {
    pub subject: Option<String>,
    pub text: String,
}

/// The run's paragraph, then each type's line, in [`OrderFacts::types`]'
/// order.
#[must_use]
pub fn texts(f: &OrderFacts) -> Vec<OrderText> {
    let run = OrderText {
        subject: None,
        text: run_template(f),
    };
    let types = f.types().into_iter().map(|t| OrderText {
        text: type_template(f, &t),
        subject: Some(t),
    });
    std::iter::once(run).chain(types).collect()
}

/// `{out}.trajectory.describe.md`.
#[must_use]
pub fn md_path(out: &str) -> PathBuf {
    PathBuf::from(format!("{out}.trajectory.describe.md"))
}

/// The texts as Markdown: the paragraph, then a section per type.
#[must_use]
pub fn markdown(from: &str, texts: &[OrderText]) -> String {
    let mut s = format!("# Trajectory description\n\nFrom `{from}`.\n\n");
    for d in texts {
        if let Some(t) = &d.subject {
            s += &format!("## {t}\n\n");
        }
        s += &format!("{}\n\n", d.text);
    }
    s
}

#[cfg(test)]
#[path = "tests/trajectory.rs"]
pub(crate) mod tests;
