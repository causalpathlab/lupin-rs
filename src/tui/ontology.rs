//! Browsing the Cell Ontology around a term, to label a cluster with any
//! term, not only the panel's: the term's parents, the term, its children,
//! or the terms a search finds.

use crate::annotate::celltype_tree::ClTerms;
use crate::annotate::markers::label_key;
use crate::annotate::panel_tree::PanelTree;
use std::collections::BTreeSet;

/// Search hits listed.
const HITS: usize = 200;

/// Where a row stands relative to the focused term.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Role {
    Parent,
    Focus,
    Child,
    Hit,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Row {
    pub id: String,
    pub role: Role,
}

/// The ontology pane's state.
pub struct OntologyView {
    pub focus: String,
    pub rows: Vec<Row>,
    pub sel: usize,
    /// A search being typed.
    pub typing: Option<String>,
    /// The query whose hits are listed, when they are.
    pub query: Option<String>,
}

impl OntologyView {
    /// Focused on `focus`, the cursor on it.
    #[must_use]
    pub fn at(cl: &ClTerms, focus: &str) -> Self {
        let mut v = Self {
            focus: focus.to_string(),
            rows: Vec::new(),
            sel: 0,
            typing: None,
            query: None,
        };
        v.refocus(cl, focus);
        v
    }

    /// Focus on `id`: its parents, itself and its children, the cursor on it.
    pub fn refocus(&mut self, cl: &ClTerms, id: &str) {
        self.focus = id.to_string();
        self.query = None;
        let mut rows: Vec<Row> = cl
            .parents(id)
            .iter()
            .map(|p| Row {
                id: p.clone(),
                role: Role::Parent,
            })
            .collect();
        self.sel = rows.len();
        rows.push(Row {
            id: id.to_string(),
            role: Role::Focus,
        });
        rows.extend(cl.children(id).iter().map(|c| Row {
            id: c.clone(),
            role: Role::Child,
        }));
        self.rows = rows;
    }

    /// List the terms `query` finds.
    pub fn search(&mut self, cl: &ClTerms, query: &str) {
        self.rows = cl
            .search(query, HITS)
            .into_iter()
            .map(|id| Row {
                id,
                role: Role::Hit,
            })
            .collect();
        self.sel = 0;
        self.query = Some(query.to_string());
    }

    #[must_use]
    pub fn selected(&self) -> Option<&Row> {
        self.rows.get(self.sel)
    }

    /// Descend into the selected row's term.
    pub fn enter(&mut self, cl: &ClTerms) {
        if let Some(id) = self.selected().map(|r| r.id.clone()) {
            self.refocus(cl, &id);
        }
    }

    /// Up to the focus's first parent; out of a search, back to the focus.
    pub fn up(&mut self, cl: &ClTerms) {
        if self.query.is_some() {
            let focus = self.focus.clone();
            self.refocus(cl, &focus);
        } else if let Some(p) = cl.parents(&self.focus).first().cloned() {
            let from = self.focus.clone();
            self.refocus(cl, &p);
            // Keep the cursor on the term we came from.
            if let Some(i) = self.rows.iter().position(|r| r.id == from) {
                self.sel = i;
            }
        }
    }
}

/// The label a cluster gets from CL term `id`: the panel type sitting on
/// exactly that term when there is one, so a pick agrees with the panel,
/// else the term's name.
#[must_use]
pub fn term_label(cl: &ClTerms, tree: &PanelTree, id: &str) -> String {
    // The first by name when several panel types sit on one term, so a pick
    // names the same label every session.
    tree.typed_terms()
        .filter(|(_, t)| *t == id)
        .map(|(l, _)| l)
        .min()
        .map_or_else(|| label_key(cl.name(id).unwrap_or(id)), String::from)
}

/// The CL term a label names: a panel type's, else the term it matches
/// by name, synonym or abbreviation.
#[must_use]
pub fn term_of(cl: &ClTerms, tree: &PanelTree, label: &str) -> Option<String> {
    if let Some(id) = tree.cl_of(label) {
        return Some(id.to_string());
    }
    let (mapped, _) = cl.map_labels([label]);
    mapped.into_values().next()
}

/// Each panel type on a CL term with every term above it: which types a
/// term covers.
#[must_use]
pub fn type_ancestry(cl: &ClTerms, tree: &PanelTree) -> Vec<(String, BTreeSet<String>)> {
    tree.typed_terms()
        .map(|(l, id)| (l.to_string(), cl.ancestors_or_self(id)))
        .collect()
}

/// Mixed labels: a label that stands for several cell types a cluster's
/// evidence cannot separate, one per line of `mixed_labels.tsv`
/// (`name<TAB>type<TAB>type…`), read from the user's config and the project's
/// `lupin/` folder and added to from the TUI. Only a listed label is a mix:
/// a `+` inside a panel type's name (`CD14+ Monocytes`) means nothing.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Mixed {
    by_name: std::collections::BTreeMap<String, Vec<String>>,
}

/// The file mixed labels are kept in.
pub const MIXED: &str = "mixed_labels.tsv";

/// How mixed-label names compare: as lupin labels, ignoring case.
fn name_key(name: &str) -> String {
    label_key(&name.trim().to_lowercase())
}

impl Mixed {
    fn extend(&mut self, text: &str) {
        for line in text.lines().map(str::trim) {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut cols = line.split('\t').map(str::trim);
            let name = cols.next().unwrap_or_default();
            let parts: Vec<String> = cols.filter(|p| !p.is_empty()).map(label_key).collect();
            if !name.is_empty() && parts.len() > 1 {
                self.by_name.insert(name_key(name), parts);
            }
        }
    }

    /// The user's and the project's mixed labels.
    pub fn load(search: &crate::manifest::data_files::SearchPath) -> anyhow::Result<Self> {
        let mut m = Self::default();
        for text in search.user_and_project(MIXED)? {
            m.extend(&text);
        }
        Ok(m)
    }

    /// Name `parts` `name` from now on, and append it to `file` (once).
    pub fn add(
        &mut self,
        name: &str,
        parts: &[String],
        file: &std::path::Path,
    ) -> anyhow::Result<()> {
        if self
            .by_name
            .get(&name_key(name))
            .is_some_and(|p| p == parts)
        {
            return Ok(());
        }
        crate::manifest::data_files::append_line(
            file,
            "Mixed cell-type labels: the name, then the types it stands for, tab-separated.",
            &std::iter::once(name)
                .chain(parts.iter().map(String::as_str))
                .collect::<Vec<_>>()
                .join("\t"),
        )?;
        self.by_name.insert(name_key(name), parts.to_vec());
        Ok(())
    }

    /// The cell types `label` stands for, when it is a listed mix.
    #[must_use]
    pub fn parts(&self, label: &str) -> Option<Vec<String>> {
        self.by_name.get(&name_key(label)).cloned()
    }
}

#[cfg(test)]
#[path = "tests/ontology.rs"]
mod tests;
