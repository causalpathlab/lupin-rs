//! Cell Ontology terms suggested for the panel types the ontology has no
//! term for ([`crate::annotate::term_match`]): `c` on such a type in the
//! tree pane offers its nearest terms with their scores; one is taken only
//! once the user confirms it, into their own alias table, and the tree is
//! built again with it. The ranking is a child `lupin data suggest-cl` job,
//! which fetches the text model and encodes the ontology the first time.

use super::menu::{numbered_but, Action, Choice, Menu};
use crate::annotate::celltype_tree::ClTerms;
use crate::annotate::markers::label_key;
use crate::annotate::panel_tree::PanelTree;
use crate::annotate::term_match::{Suggested, Suggestion};
use anyhow::{Context, Result};
use std::collections::BTreeMap;
use std::path::Path;

/// Each panel type's suggested terms, by its [`label_key`].
#[derive(Debug, Default)]
pub struct Suggestions {
    by_label: BTreeMap<String, Vec<Suggestion>>,
}

impl Suggestions {
    /// What a `lupin data suggest-cl` job wrote to `file`.
    pub fn read(file: &Path) -> Result<Self> {
        let text =
            std::fs::read_to_string(file).with_context(|| format!("reading {}", file.display()))?;
        let s: Suggested = serde_json::from_str(&text)
            .with_context(|| format!("{} is not the suggestions' JSON", file.display()))?;
        Ok(Self {
            by_label: (s.labels.into_iter())
                .map(|(l, v)| (label_key(&l), v))
                .collect(),
        })
    }

    #[must_use]
    pub fn of(&self, label: &str) -> Option<&[Suggestion]> {
        self.by_label.get(&label_key(label)).map(Vec::as_slice)
    }
}

/// Whether `label` is a panel type the ontology has no term for: what `c`
/// asks before ranking, and [`unmatched`] of every type. Asked of the
/// label's own reading ([`super::ontology::term_of`]), not of where the
/// tree put it: a tree built from the markers puts no type on a term.
#[must_use]
pub fn needs_term(cl: &ClTerms, tree: &PanelTree, label: &str) -> bool {
    tree.is_type(label) && super::ontology::term_of(cl, tree, label).is_none()
}

/// The panel types with no Cell Ontology term, by their [`label_key`]s, the
/// keys [`Suggestions`] gives them back by.
#[must_use]
pub fn unmatched(cl: &ClTerms, tree: &PanelTree) -> Vec<String> {
    let mut out: Vec<String> = (tree.nodes.iter())
        .flat_map(|n| &n.labels)
        .filter(|l| needs_term(cl, tree, l))
        .map(|l| label_key(l))
        .collect();
    out.sort();
    out.dedup();
    out
}

/// The choice among `label`'s suggested terms, or none of them; `file` is
/// where a confirmed one is written.
#[must_use]
pub fn menu(label: &str, suggested: &[Suggestion], file: &Path) -> Menu {
    let items = suggested
        .iter()
        .map(|s| {
            (
                format!("{} ({}) · {:.2}", s.name, s.id, s.score),
                format!(
                    "cosine {:.2} to its nearest name; enter, then confirm, writes {label} → {} to {}",
                    s.score,
                    s.id,
                    super::app::file_name(file)
                ),
                Action::Suggested {
                    label: label.to_string(),
                    id: s.id.clone(),
                    note: format!(
                        "suggested by text embedding ({:.2}), confirmed in the lupin TUI",
                        s.score
                    ),
                    file: file.to_path_buf(),
                },
            )
        })
        .collect();
    Menu::new(
        format!(
            "{label} has no Cell Ontology term. Its nearest by text embedding (a suggestion: you decide):"
        ),
        with_search_and_none(label, items, file),
        format!("no term taken for {label}"),
    )
}

/// The terms a search for `query` found for `label` (`(id, name)`), to
/// confirm one into `file` as a suggested one is.
#[must_use]
pub fn found_menu(label: &str, query: &str, found: &[(String, String)], file: &Path) -> Menu {
    let items = (found.iter())
        .map(|(id, name)| {
            (
                format!("{name} ({id})"),
                format!(
                    "enter, then confirm, writes {label} → {id} to {}",
                    super::app::file_name(file)
                ),
                Action::Suggested {
                    label: label.to_string(),
                    id: id.clone(),
                    note: format!("found by searching \"{query}\", confirmed in the lupin TUI"),
                    file: file.to_path_buf(),
                },
            )
        })
        .collect();
    Menu::new(
        format!("Cell Ontology terms for \"{query}\": which is {label}?"),
        with_search_and_none(label, items, file),
        format!("no term taken for {label}"),
    )
}

/// `items` numbered, then a search of the ontology and none of them.
fn with_search_and_none(
    label: &str,
    items: Vec<(String, String, Action)>,
    file: &Path,
) -> Vec<Choice> {
    let mut choices = numbered_but(items, "sn");
    choices.push(Choice {
        key: Some('s'),
        label: "search the ontology…".into(),
        detail: "type words of the term's name to find it".into(),
        action: Action::SearchTerm {
            label: label.to_string(),
            file: file.to_path_buf(),
        },
    });
    choices.push(Choice {
        key: Some('n'),
        label: "none of these".into(),
        detail: "no term fits: kept with the answers, so it is not asked again".into(),
        action: Action::NoneSuggested {
            label: label.to_string(),
        },
    });
    choices
}

#[cfg(test)]
#[path = "tests/suggest.rs"]
mod tests;
