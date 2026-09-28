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
    tree.typed_terms()
        .find(|(_, t)| *t == id)
        .map(|(l, _)| l.to_string())
        .unwrap_or_else(|| label_key(cl.name(id).unwrap_or(id)))
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

#[cfg(test)]
#[path = "tests/ontology.rs"]
mod tests;
