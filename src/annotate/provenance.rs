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
                let term = (placed.on())
                    .chain(&placed.kept)
                    .find(|(m, _)| *m == l)
                    .map(|(_, id)| id.clone());
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
        Self {
            ontology_release: terms.release.clone(),
            labels,
            answers: (types.iter())
                .filter_map(|l| Some((l.clone(), terms.joined(l)?)))
                .collect(),
            questions: TypeTree::join_questions(terms, panel),
        }
    }
}

/// A methods paragraph drafted from `r`, by code.
#[must_use]
pub fn methods(r: &Record) -> String {
    let count = |f: &dyn Fn(&How) -> bool| r.labels.iter().filter(|l| f(&l.how)).count();
    let names = count(&|h| *h == How::Name);
    let aliases = count(&|h| *h == How::Alias);
    let word = r.labels.iter().find_map(|l| match &l.how {
        How::PanelWord(w) => Some(w.clone()),
        _ => None,
    });
    let by_word = count(&|h| matches!(h, How::PanelWord(_)));
    let under = count(&|h| *h == How::Under);
    let kept = count(&|h| *h == How::Kept);
    let nowhere = count(&|h| *h == How::Nowhere);
    let release = r
        .ontology_release
        .as_deref()
        .unwrap_or("of unknown release");
    let mut parts = vec![format!("{names} by a term's name or synonym")];
    if aliases > 0 {
        parts.push(format!("{aliases} by a curated alias"));
    }
    if let Some(w) = word {
        parts.push(format!(
            "{by_word} by completing the word \"{w}\" the panel leaves off its labels"
        ));
    }
    let mut text = format!(
        "The marker panel's {n} cell-type labels were matched to the Cell Ontology (release {release}) \
         with lupin's matching rules (cl_matching.json): {}.",
        parts.join(", "),
        n = r.labels.len()
    );
    if under > 0 {
        text += &format!(
            " Of the labels no term matched, {under} {} placed under the term named by the longest run of its words, the words left off being ones the ontology uses for kinds of a term.",
            if under == 1 { "was" } else { "were" }
        );
    }
    if kept > 0 {
        text += &format!(
            " {kept} {} kept to the lineage of the term its longest run of words names and grouped by shared markers.",
            if kept == 1 { "was" } else { "were" }
        );
    }
    if nowhere > 0 {
        text += &format!(
            " {nowhere} matched no term and {} grouped by shared marker genes.",
            if nowhere == 1 { "was" } else { "were" }
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
