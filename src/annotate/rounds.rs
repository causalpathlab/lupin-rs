//! Annotation rounds: the per-cluster digest a reviewer reads, the decisions
//! they make, and the history that keeps the reasoning behind each one.
//!
//! A round is a manifest. `lupin annotate` writes the first; each
//! `lupin relabel` reads a round plus a decisions file and writes the next,
//! whose `annotate.source` points back. Cluster ids are the integers in
//! `cluster.clusters` and are never reused: a merge takes a fresh id, so a
//! cluster's history stays attached to one id across rounds.
//!
//! Nothing here reads or writes files; [`crate::manifest::rounds`] does.

use crate::annotate::markers::label_key;
use anyhow::{bail, ensure, Result};
use serde::{Deserialize, Deserializer, Serialize};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};

pub type ClusterId = u32;

/// How many candidates of each kind a digest entry keeps.
const TOP: usize = 5;
/// How many evidence items a history entry keeps from its decision.
const TOP_EVIDENCE: usize = 3;

////////////
// digest //
////////////

/// What a reviewer needs to decide on one cluster. Keyed by cluster id in
/// `{out}.cluster_summary.json`; every list is sorted best first.
#[derive(Serialize, Deserialize, Default, Debug, Clone, PartialEq)]
pub struct Digest {
    pub size: usize,
    /// The label most of the cluster's cells carry in this round.
    pub label: Option<String>,
    pub calls: Vec<Call>,
    pub terms: Vec<Term>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cl: Option<ClPlacement>,
    /// The best call beside the cluster's label, to flag a label the
    /// evidence does not back.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub evidence: Option<EvidenceCheck>,
}

/// The top call and whether the cluster's label agrees with it.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct EvidenceCheck {
    pub top: String,
    pub q: Option<f32>,
    pub agrees: bool,
}

/// A candidate label with its FDR q-value, where known.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Call {
    pub label: String,
    pub q: Option<f32>,
}

/// A gene-set term (GO or GMT) and its effect on the cluster.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Term {
    pub source: String,
    pub term: String,
    pub effect: f32,
    pub q: Option<f32>,
    /// The term's test as the cell types are tested, when the pass ran it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub p: Option<f32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nes: Option<f32>,
}

/// Where the Cell Ontology walk placed the cluster.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct ClPlacement {
    pub id: String,
    pub name: String,
    pub abstained: bool,
}

/// A cluster × label table, rows already parsed to cluster ids.
pub struct Table {
    pub rows: Vec<ClusterId>,
    pub cols: Vec<String>,
    /// Row-major, `rows.len() × cols.len()`.
    pub values: Vec<f32>,
}

impl Table {
    pub(crate) fn row(&self, id: ClusterId) -> Option<&[f32]> {
        let r = self.rows.iter().position(|&x| x == id)?;
        let w = self.cols.len();
        Some(&self.values[r * w..(r + 1) * w])
    }
}

/// Everything a digest is built from; each part is optional.
#[derive(Default)]
pub struct Evidence {
    /// Cluster × label FDR q-values.
    pub q: Option<Table>,
    /// Per cluster, terms in rank order.
    pub terms: BTreeMap<ClusterId, Vec<Term>>,
    pub cl: BTreeMap<ClusterId, ClPlacement>,
}

/// Labels a table may carry that are not a call.
const NOT_A_CALL: &[&str] = &["unassigned"];

/// One digest entry per cluster id present in `clusters`.
#[must_use]
pub fn digest(
    clusters: &[Option<ClusterId>],
    labels: &[Option<String>],
    ev: &Evidence,
) -> BTreeMap<ClusterId, Digest> {
    let mut out: BTreeMap<ClusterId, Digest> = BTreeMap::new();
    let mut votes: BTreeMap<ClusterId, BTreeMap<&str, usize>> = BTreeMap::new();
    for (i, id) in clusters.iter().enumerate() {
        let Some(id) = *id else { continue };
        out.entry(id).or_default().size += 1;
        if let Some(Some(l)) = labels.get(i) {
            *votes.entry(id).or_default().entry(l.as_str()).or_default() += 1;
        }
    }
    for (id, d) in &mut out {
        d.label = votes.get(id).and_then(|v| {
            v.iter()
                .max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0)))
                .map(|(l, _)| (*l).to_string())
        });
        d.calls = calls_for(*id, ev);
        d.terms = ev
            .terms
            .get(id)
            .map(|t| t.iter().take(TOP).cloned().collect())
            .unwrap_or_default();
        d.cl = ev.cl.get(id).cloned();
        d.evidence = d.calls.first().map(|c| EvidenceCheck {
            top: c.label.clone(),
            q: c.q,
            agrees: d
                .label
                .as_deref()
                .is_some_and(|l| label_key(l) == label_key(&c.label)),
        });
    }
    out
}

fn calls_for(id: ClusterId, ev: &Evidence) -> Vec<Call> {
    let mut by_label: BTreeMap<String, Call> = BTreeMap::new();
    let mut take = |t: &Option<Table>, set: fn(&mut Call, f32)| {
        let Some(t) = t else { return };
        let Some(row) = t.row(id) else { return };
        for (c, &v) in t.cols.iter().zip(row) {
            if NOT_A_CALL.contains(&c.as_str()) || !v.is_finite() {
                continue;
            }
            let call = by_label.entry(c.clone()).or_insert_with(|| Call {
                label: c.clone(),
                q: None,
            });
            set(call, v);
        }
    };
    take(&ev.q, |c, v| c.q = Some(v));
    let mut calls: Vec<Call> = by_label.into_values().collect();
    // Smallest q first; unknowns last.
    calls.sort_by(|a, b| {
        let q = |c: &Call| c.q.unwrap_or(f32::INFINITY);
        q(a).total_cmp(&q(b))
    });
    calls.truncate(TOP);
    calls
}

///////////////
// decisions //
///////////////

#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Action {
    /// Give one cluster a label.
    Label,
    /// Join clusters into one with a fresh id, optionally labelling it.
    Merge,
    /// Leave a cluster as it is, recording why.
    Keep,
    /// Not supported yet; refused with an explanation.
    Split,
    /// Add features to a cell type's markers (a new type if it has none).
    MarkersAdd,
    /// Remove features from a cell type's markers.
    MarkersDrop,
}

impl Action {
    /// The name decisions files use.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Action::Label => "label",
            Action::Merge => "merge",
            Action::Keep => "keep",
            Action::Split => "split",
            Action::MarkersAdd => "markers_add",
            Action::MarkersDrop => "markers_drop",
        }
    }

    /// Edits the marker panel rather than cluster labels.
    #[must_use]
    pub fn edits_markers(self) -> bool {
        matches!(self, Action::MarkersAdd | Action::MarkersDrop)
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum DecidedBy {
    User,
    AgentProposedUserAccepted,
    UserOverride,
}

impl DecidedBy {
    /// The name decisions files use.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            DecidedBy::User => "user",
            DecidedBy::AgentProposedUserAccepted => "agent_proposed_user_accepted",
            DecidedBy::UserOverride => "user_override",
        }
    }
}

/// One line of a decisions file.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Decision {
    /// The cluster(s) acted on: an id or a list, as integers or strings.
    /// Empty for marker edits.
    #[serde(
        default,
        alias = "clusters",
        deserialize_with = "de_ids",
        skip_serializing_if = "Vec::is_empty"
    )]
    pub cluster: Vec<ClusterId>,
    pub action: Action,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    /// Marker edits: the features added to or dropped from `label`.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub features: Vec<String>,
    /// What the decision rests on, typically items quoted from the digest.
    #[serde(default)]
    pub evidence: Vec<Value>,
    /// Options weighed and why each was not taken.
    #[serde(default)]
    pub alternatives: Vec<Value>,
    #[serde(default)]
    pub rationale: String,
    pub decided_by: DecidedBy,
    /// Filled in when the decision is applied, if absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub timestamp: Option<String>,
    /// The round whose cluster ids this decision names, relative to the
    /// decisions file. When given, it must be the round being relabelled, so a
    /// decision made on a stale view cannot land on renumbered ids.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub round: Option<String>,
}

fn de_ids<'de, D: Deserializer<'de>>(d: D) -> Result<Vec<ClusterId>, D::Error> {
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum Id {
        N(ClusterId),
        S(String),
    }
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum OneOrMany {
        One(Id),
        Many(Vec<Id>),
    }
    let ids = match OneOrMany::deserialize(d)? {
        OneOrMany::One(id) => vec![id],
        OneOrMany::Many(ids) => ids,
    };
    ids.into_iter()
        .map(|id| match id {
            Id::N(n) => Ok(n),
            Id::S(s) => parse_cluster_id(&s)
                .ok_or_else(|| serde::de::Error::custom(format!("not a cluster id: `{s}`"))),
        })
        .collect()
}

/// A cluster id from `12`, `K12` or `C12`.
#[must_use]
pub fn parse_cluster_id(s: &str) -> Option<ClusterId> {
    s.trim_start_matches(|c: char| c.is_ascii_alphabetic())
        .parse()
        .ok()
}

/////////////
// history //
/////////////

/// One step in a cluster's history, newest first in the history file.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct HistoryEntry {
    /// The round's manifest, relative to the history file.
    pub round: String,
    pub action: Action,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    pub rationale: String,
    pub decided_by: DecidedBy,
    pub timestamp: String,
    #[serde(default)]
    pub evidence: Vec<Value>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub merged_from: Option<Vec<ClusterId>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub merged_into: Option<ClusterId>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub split_from: Option<ClusterId>,
    /// Marker edits: the features added or dropped.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub features: Option<Vec<String>>,
}

pub type History = BTreeMap<ClusterId, Vec<HistoryEntry>>;

/// `newer` entries go in front of each cluster's `older` ones.
#[must_use]
pub fn prepend_history<K: Ord>(
    older: BTreeMap<K, Vec<HistoryEntry>>,
    newer: BTreeMap<K, Vec<HistoryEntry>>,
) -> BTreeMap<K, Vec<HistoryEntry>> {
    let mut out = newer;
    for (id, mut entries) in older {
        out.entry(id).or_default().append(&mut entries);
    }
    out
}

/// The first id no round has used: past every current id and every id the
/// history mentions.
#[must_use]
pub fn next_cluster_id(clusters: &[Option<ClusterId>], history: &History) -> ClusterId {
    let current = clusters.iter().flatten().copied().max();
    let past = history
        .iter()
        .flat_map(|(id, es)| {
            std::iter::once(*id).chain(es.iter().flat_map(|e| {
                e.merged_from
                    .iter()
                    .flatten()
                    .copied()
                    .chain(e.merged_into)
                    .chain(e.split_from)
            }))
        })
        .max();
    current.max(past).map_or(0, |m| m + 1)
}

//////////////
// applying //
//////////////

/// Apply `decisions` to one round's per-cell cluster ids and labels, in
/// place, and return the history entries they make, tagged with `round`.
///
/// Every decision needs a rationale, may name only clusters present in this
/// round, and no cluster may be named twice. `now` stamps decisions that
/// carry no timestamp.
pub fn apply(
    decisions: &mut [Decision],
    clusters: &mut [Option<ClusterId>],
    labels: &mut [Option<String>],
    mut next_id: ClusterId,
    round: &str,
    now: &str,
) -> Result<History> {
    let present: BTreeSet<ClusterId> = clusters.iter().flatten().copied().collect();
    let mut seen = BTreeSet::new();
    for (n, d) in decisions.iter().enumerate() {
        let line = n + 1;
        ensure!(
            !d.rationale.trim().is_empty(),
            "decision {line}: a rationale is required"
        );
        if d.action.edits_markers() {
            ensure!(
                d.cluster.is_empty(),
                "decision {line}: a marker edit names a cell type, not clusters"
            );
            ensure!(
                d.label.as_deref().is_some_and(|l| !l.trim().is_empty()),
                "decision {line}: a marker edit needs the cell type as `label`"
            );
            ensure!(
                d.features.iter().any(|f| !f.trim().is_empty()),
                "decision {line}: a marker edit needs `features`"
            );
            continue;
        }
        ensure!(!d.cluster.is_empty(), "decision {line}: no cluster named");
        for id in &d.cluster {
            ensure!(
                present.contains(id),
                "decision {line}: cluster {id} is not in this round"
            );
            ensure!(
                seen.insert(*id),
                "decision {line}: cluster {id} is already decided in this round"
            );
        }
        match d.action {
            Action::Label => {
                ensure!(
                    d.cluster.len() == 1,
                    "decision {line}: `label` takes one cluster; use `merge` to join several"
                );
                ensure!(
                    d.label.as_deref().is_some_and(|l| !l.trim().is_empty()),
                    "decision {line}: `label` needs a label"
                );
            }
            Action::Merge => ensure!(
                d.cluster.len() >= 2,
                "decision {line}: `merge` needs at least two clusters"
            ),
            Action::Keep => {}
            Action::Split => bail!(
                "decision {line}: `split` is not supported yet; re-cluster that cluster's \
                 cells and pass the result to `lupin annotate --clusters`"
            ),
            Action::MarkersAdd | Action::MarkersDrop => unreachable!("handled above"),
        }
    }

    let mut history = History::new();
    for d in decisions.iter_mut() {
        let timestamp = d.timestamp.get_or_insert_with(|| now.to_string()).clone();
        let entry = |label: Option<String>| HistoryEntry {
            round: round.to_string(),
            action: d.action,
            label,
            rationale: d.rationale.clone(),
            decided_by: d.decided_by,
            timestamp: timestamp.clone(),
            evidence: d.evidence.iter().take(TOP_EVIDENCE).cloned().collect(),
            merged_from: None,
            merged_into: None,
            split_from: None,
            features: None,
        };
        match d.action {
            Action::Label | Action::Keep => {
                let id = d.cluster[0];
                if let Some(l) = &d.label {
                    relabel_cells(clusters, labels, id, &label_key(l));
                }
                for &id in &d.cluster {
                    history.entry(id).or_default().push(entry(d.label.clone()));
                }
            }
            Action::Merge => {
                let new_id = next_id;
                next_id += 1;
                for c in clusters.iter_mut() {
                    if c.is_some_and(|id| d.cluster.contains(&id)) {
                        *c = Some(new_id);
                    }
                }
                if let Some(l) = &d.label {
                    relabel_cells(clusters, labels, new_id, &label_key(l));
                }
                for &old in &d.cluster {
                    history.entry(old).or_default().push(HistoryEntry {
                        merged_into: Some(new_id),
                        ..entry(d.label.clone())
                    });
                }
                history.entry(new_id).or_default().push(HistoryEntry {
                    merged_from: Some(d.cluster.clone()),
                    ..entry(d.label.clone())
                });
            }
            Action::Split => unreachable!("refused above"),
            // Applied by [`apply_markers`].
            Action::MarkersAdd | Action::MarkersDrop => {}
        }
    }
    Ok(history)
}

/// Marker history: every round's marker edits per cell type, newest first.
pub type MarkerHistory = BTreeMap<String, Vec<HistoryEntry>>;

/// Apply the marker edits among `decisions` (already validated and stamped
/// by [`apply`]) to `(feature, cell type)` pairs, in place. Adding a pair
/// already there changes nothing; dropping one that is not there is refused,
/// since it is most likely a typo. Returns whether anything was edited, and
/// the history entries per cell type.
pub fn apply_markers(
    decisions: &[Decision],
    markers: &mut Vec<(String, String)>,
    round: &str,
) -> Result<(bool, MarkerHistory)> {
    let mut history = MarkerHistory::new();
    let mut edited = false;
    for (n, d) in decisions.iter().enumerate() {
        if !d.action.edits_markers() {
            continue;
        }
        let ty = d.label.as_deref().unwrap_or_default();
        let key = label_key(ty);
        let features: Vec<&str> = d
            .features
            .iter()
            .map(|f| f.trim())
            .filter(|f| !f.is_empty())
            .collect();
        let has =
            |m: &[(String, String)], f: &str| m.iter().any(|(g, t)| g == f && label_key(t) == key);
        if d.action == Action::MarkersAdd {
            for f in &features {
                if !has(markers, f) {
                    markers.push(((*f).to_string(), key.clone()));
                }
            }
        } else {
            let missing: Vec<&&str> = features.iter().filter(|f| !has(markers, f)).collect();
            ensure!(
                missing.is_empty(),
                "decision {}: {missing:?} are not markers of `{ty}`",
                n + 1
            );
            markers.retain(|(g, t)| !(label_key(t) == key && features.contains(&g.as_str())));
        }
        edited = true;
        history.entry(key.clone()).or_default().push(HistoryEntry {
            round: round.to_string(),
            action: d.action,
            label: Some(key.clone()),
            rationale: d.rationale.clone(),
            decided_by: d.decided_by,
            timestamp: d.timestamp.clone().unwrap_or_default(),
            evidence: d.evidence.iter().take(TOP_EVIDENCE).cloned().collect(),
            merged_from: None,
            merged_into: None,
            split_from: None,
            features: Some(features.iter().map(|f| (*f).to_string()).collect()),
        });
    }
    Ok((edited, history))
}

fn relabel_cells(
    clusters: &[Option<ClusterId>],
    labels: &mut [Option<String>],
    id: ClusterId,
    label: &str,
) {
    for (c, l) in clusters.iter().zip(labels.iter_mut()) {
        if *c == Some(id) {
            *l = Some(label.to_string());
        }
    }
}

/// `YYYY-MM-DDTHH:MM:SSZ` for a Unix time in seconds.
#[must_use]
pub fn utc_timestamp(unix_secs: u64) -> String {
    let days = i64::try_from(unix_secs / 86_400).unwrap_or(0);
    let secs = unix_secs % 86_400;
    // Days since 1970-01-01 to a civil date (Howard Hinnant's algorithm).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z",
        secs / 3600,
        (secs % 3600) / 60,
        secs % 60
    )
}

#[cfg(test)]
#[path = "rounds_tests.rs"]
mod tests;
