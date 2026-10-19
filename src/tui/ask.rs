//! The questions lupin asks the user on a panel's load, one at a time as a
//! menu while the panel tree is idle: first whether a type left alone (its
//! siblings went under a nearer class) joins them, then a term for each type
//! the ontology has none for, from its suggestions. Answers are written where
//! the next run reads them (`cl_groups.tsv`, the user's alias table), so each
//! is asked once and every pass records it ([`crate::annotate::provenance`]);
//! "none of these" for a term is kept with the answers too.

use super::menu::{Action, Choice, Menu};
use crate::annotate::celltype_tree::{ClTerms, JoinQuestion, TypeTree};
use crate::annotate::panel_tree::PanelTree;
use std::collections::VecDeque;
use std::path::Path;

pub enum Question {
    Join(JoinQuestion),
    Term(String),
}

/// A panel's questions while they are asked.
#[derive(Default)]
pub struct Questions {
    /// Those left, in the order they are asked ([`queue`]).
    pub left: VecDeque<Question>,
    /// The one the open menu asks.
    pub asking: Option<Question>,
    /// A ranking was started for their terms: one per panel.
    pub ranked: bool,
}

/// The panel's questions, in the order they are asked.
#[must_use]
pub fn queue(cl: &ClTerms, tree: &PanelTree, panel: &[(String, String)]) -> VecDeque<Question> {
    (TypeTree::join_questions(cl, panel).into_iter())
        .map(Question::Join)
        .chain(
            (super::suggest::unmatched(cl, tree).into_iter())
                .filter(|l| !cl.declined(l))
                .map(Question::Term),
        )
        .collect()
}

/// Join `q`'s siblings under `class_name`, or stand alone; written to `file`.
#[must_use]
pub fn join_menu(q: &JoinQuestion, class_name: &str, file: &Path) -> Menu {
    let choice = |key, label: String, detail: String, join| Choice {
        key: Some(key),
        label,
        detail,
        action: Action::Join {
            label: q.label.clone(),
            join,
            file: file.to_path_buf(),
        },
    };
    let to = super::app::file_name(file);
    Menu::new(
        format!(
            "{} stands alone: the types it shares {class_name} with went under a narrower class ({}). Group it with them?",
            q.label, q.group
        ),
        vec![
            choice('1', format!("join them under {class_name}"), format!("writes {}<TAB>join to {to}", q.label), true),
            choice('2', "stand alone".into(), format!("writes {}<TAB>alone to {to}", q.label), false),
        ],
        format!("no answer for {}", q.label),
    )
}

#[cfg(test)]
#[path = "tests/ask.rs"]
mod tests;
