//! The record of how a panel's labels were put on the Cell Ontology: what
//! placed each label (an alias, the matching rules, the word the panel
//! leaves off, the run of its words, or its markers), the questions left
//! for the user and the answers given; and, drafted from it by code, the
//! paragraph a methods section needs. The same record gives the same text.

use crate::annotate::celltype_tree::{panel_types, ClTerms, How, JoinQuestion, TypeTree};
use serde::Serialize;

/// One label: what placed it, and on which term.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LabelRecord {
    pub label: String,
    pub how: How,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub term: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub term_name: Option<String>,
}

/// How a panel's labels were put on the ontology.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Record {
    pub ontology_release: Option<String>,
    /// Whether the groups came from the ontology; else from shared markers,
    /// too few labels being on terms.
    pub groups_from_ontology: bool,
    pub labels: Vec<LabelRecord>,
    /// The user's answers on types left alone: joined (`true`) or alone.
    pub answers: Vec<(String, bool)>,
    /// The questions still open.
    pub questions: Vec<JoinQuestion>,
}

impl Record {
    /// The record of `panel` on `terms`.
    #[must_use]
    pub fn of(terms: &ClTerms, panel: &[(String, String)]) -> Self {
        let types = panel_types(panel);
        let placed = terms.place_labels(types.iter().map(String::as_str));
        let labels = (types.iter())
            .map(|l| {
                let term = placed.term(l).cloned();
                LabelRecord {
                    label: l.clone(),
                    how: placed.how(l),
                    term_name: term
                        .as_deref()
                        .and_then(|id| terms.name(id))
                        .map(String::from),
                    term,
                }
            })
            .collect();
        let (answered, open): (Vec<_>, Vec<_>) =
            (TypeTree::left_alone(terms, &placed).into_iter()).partition(|(_, a)| a.is_some());
        Self {
            ontology_release: terms.release.clone(),
            // As `TypeTree::from_placed` decides.
            groups_from_ontology: placed.on_terms().is_some(),
            labels,
            answers: (answered.into_iter())
                .filter_map(|(q, a)| Some((q.label, a?)))
                .collect(),
            questions: open.into_iter().map(|(q, _)| q).collect(),
        }
    }
}

/// A methods paragraph drafted from `r`, by code.
#[must_use]
pub fn methods(r: &Record) -> String {
    let (mut names, mut aliases, mut by_word, mut under, mut kept, mut nowhere) =
        (0, 0, 0, 0, 0, 0);
    let mut word = None;
    for l in &r.labels {
        match &l.how {
            How::Name => names += 1,
            How::Alias => aliases += 1,
            How::PanelWord(w) => {
                by_word += 1;
                word = Some(w.as_str());
            }
            How::Under => under += 1,
            How::Kept => kept += 1,
            How::Nowhere => nowhere += 1,
        }
    }
    let was = |n: usize| if n == 1 { "was" } else { "were" };
    let release = (r.ontology_release.as_deref()).map_or_else(
        || "release not stated".to_string(),
        |v| format!("release {v}"),
    );
    let mut parts = Vec::new();
    if names > 0 {
        parts.push(format!("{names} by a term's name or synonym"));
    }
    if aliases > 0 {
        parts.push(format!("{aliases} by a curated alias"));
    }
    if let Some(w) = word {
        parts.push(format!(
            "{by_word} by completing the word \"{w}\" the panel leaves off its labels"
        ));
    }
    let n = r.labels.len();
    let matched = if parts.is_empty() {
        "no label matched a term".to_string()
    } else {
        parts.join(", ")
    };
    let mut text = format!(
        "The marker panel's {n} cell-type labels were matched to the Cell Ontology ({release}) \
         with lupin's matching rules (cl_matching.json): {matched}."
    );
    if !r.groups_from_ontology {
        text += " With too few labels on terms to group them by the ontology, the coarse groups were formed by shared marker genes.";
    }
    if under > 0 {
        text += &format!(
            " Of the labels no term matched, {under} {} placed under the term named by the longest run of its words, the words left off being ones the ontology uses for kinds of a term.",
            was(under)
        );
    }
    if kept > 0 {
        text += &format!(
            " {kept} {} kept to the lineage of the term its longest run of words names and grouped by shared markers.",
            was(kept)
        );
    }
    if nowhere > 0 {
        text += &format!(
            " {nowhere} matched no term and {} grouped by shared marker genes.",
            was(nowhere)
        );
    }
    if !r.answers.is_empty() {
        let said = |join: bool| -> Vec<&str> {
            (r.answers.iter())
                .filter(|(_, j)| *j == join)
                .map(|(l, _)| l.as_str())
                .collect()
        };
        let (joined, alone) = (said(true), said(false));
        text += " Where a type stood alone because its siblings fell under a narrower class, the user was asked whether it joins them:";
        if !joined.is_empty() {
            text += &format!(" joined: {}.", joined.join(", "));
        }
        if !alone.is_empty() {
            text += &format!(" kept alone: {}.", alone.join(", "));
        }
    }
    if !r.questions.is_empty() {
        let open: Vec<&str> = r.questions.iter().map(|q| q.label.as_str()).collect();
        text += &format!(
            " Left alone because their siblings fell under a narrower class, and not yet settled by the user, these stood as their own groups: {}.",
            open.join(", ")
        );
    }
    text
}
