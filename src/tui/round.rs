//! One annotation round as the screen works on it: each cluster's label,
//! the cell types competing for it, and the genes that set it apart; plus
//! the edits made here, which become the `relabel` decisions of the next
//! round when saved.

use crate::annotate::markers::label_key;
use crate::annotate::markers::read_marker_pairs;
use crate::annotate::rounds::{
    digest, parse_cluster_id, Action, ClusterId, DecidedBy, Decision, Evidence, Table,
};
use crate::manifest::rounds::{read_cells, read_table};
use crate::manifest::run::{self, resolve};
use anyhow::Result;
use enrichment::UNASSIGNED_LABEL;
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::traits::{IoOps, MatWithNames};
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

/// Candidates listed per cluster.
const CANDIDATES: usize = 6;

/// A cell type competing for a cluster, with its enrichment statistics.
#[derive(Clone, Debug, PartialEq)]
pub struct Candidate {
    pub label: String,
    /// Its share of the cluster's evidence (Q): a softmax of z = Φ⁻¹(1 − p)
    /// over the types that pass FDR.
    pub share: f32,
    /// fgsea's normalized enrichment score, the effect size.
    pub nes: Option<f32>,
    pub p: Option<f32>,
    /// BH q-value within the cluster.
    pub q: Option<f32>,
}

/// A genes × clusters expression table with what fold changes over it need,
/// computed once: each gene's row, each cluster's column, each gene's total,
/// and the pseudo-count.
pub struct Expression {
    pub table: MatWithNames<Mat>,
    row: HashMap<String, usize>,
    col: HashMap<ClusterId, usize>,
    totals: Vec<f32>,
    /// The table's mean: a pseudo-count on its scale, so a gene near zero
    /// elsewhere does not win on a vanishing denominator.
    eps: f32,
}

impl Expression {
    /// Index `table` (columns `K{id}`).
    #[must_use]
    pub fn new(table: MatWithNames<Mat>) -> Self {
        let row = table
            .rows
            .iter()
            .enumerate()
            .map(|(r, g)| (g.to_string(), r))
            .collect();
        let col = table
            .cols
            .iter()
            .enumerate()
            .filter_map(|(k, c)| parse_cluster_id(c).map(|id| (id, k)))
            .collect();
        let totals = (0..table.mat.nrows())
            .map(|g| table.mat.row(g).sum())
            .collect();
        let eps = table.mat.mean().max(f32::MIN_POSITIVE);
        Self {
            table,
            row,
            col,
            totals,
            eps,
        }
    }

    #[must_use]
    pub fn row(&self, gene: &str) -> Option<usize> {
        self.row.get(gene).copied()
    }

    /// The columns of clusters `ids`.
    #[must_use]
    pub fn cols(&self, ids: &BTreeSet<ClusterId>) -> Vec<usize> {
        ids.iter()
            .filter_map(|id| self.col.get(id).copied())
            .collect()
    }

    /// Gene `r`'s mean over columns `inside`.
    #[must_use]
    pub fn mean_in(&self, r: usize, inside: &[usize]) -> f32 {
        inside.iter().map(|&k| self.table.mat[(r, k)]).sum::<f32>() / inside.len().max(1) as f32
    }

    /// Gene `r`'s log2 fold change over columns `inside` against the other
    /// columns' mean; `None` without both.
    #[must_use]
    pub fn log2fc(&self, r: usize, inside: &[usize]) -> Option<f32> {
        let rest = self
            .table
            .mat
            .ncols()
            .checked_sub(inside.len())
            .filter(|&n| n > 0)?;
        if inside.is_empty() {
            return None;
        }
        let sum_in: f32 = inside.iter().map(|&k| self.table.mat[(r, k)]).sum();
        let mean_out = (self.totals[r] - sum_in) / rest as f32;
        Some(((sum_in / inside.len() as f32 + self.eps) / (mean_out + self.eps)).log2())
    }
}

/// One cluster of the round.
pub struct ClusterView {
    pub id: ClusterId,
    pub cells: usize,
    /// The round's label; `None` for unassigned.
    pub label: Option<String>,
    /// Cell types by their share of the cluster's evidence, largest first.
    pub candidates: Vec<Candidate>,
    /// Every cell type's share of the cluster's evidence.
    pub shares: Vec<(String, f32)>,
    /// Genes by [`specific_genes`], most specific first.
    pub genes: Vec<(String, f32)>,
}

impl ClusterView {
    /// The top candidate's share: low means contested.
    #[must_use]
    pub fn top_share(&self) -> f32 {
        self.candidates.first().map_or(0.0, |c| c.share)
    }

    /// Unassigned, or its top candidate holds too little of its evidence.
    #[must_use]
    pub fn flagged(&self) -> bool {
        self.label.is_none() || self.top_share() < super::app::CONTESTED
    }
}

/// An edit made on screen, not yet saved.
#[derive(Clone, Debug, PartialEq)]
pub enum Edit {
    Label {
        cluster: ClusterId,
        label: String,
        reason: String,
    },
    /// Keep the cluster's label, recording why.
    Keep { cluster: ClusterId, reason: String },
    /// Add `genes` to `label`'s markers, or drop them (`add` false).
    Markers {
        label: String,
        genes: Vec<String>,
        add: bool,
        reason: String,
    },
}

impl Edit {
    /// The cluster a label or keep decides.
    #[must_use]
    pub fn cluster(&self) -> Option<ClusterId> {
        match self {
            Self::Label { cluster, .. } | Self::Keep { cluster, .. } => Some(*cluster),
            Self::Markers { .. } => None,
        }
    }
}

pub struct RoundView {
    pub manifest: PathBuf,
    pub clusters: Vec<ClusterView>,
    pub cell_names: Vec<Box<str>>,
    pub cell_clusters: Vec<Option<ClusterId>>,
    /// Genes × clusters, as the pass aggregated them.
    pub expression: Option<Expression>,
    /// The round's marker panel: label key → genes.
    pub markers: BTreeMap<String, BTreeSet<String>>,
    /// Cells no cluster holds.
    pub loose_cells: usize,
    /// Clusters an earlier round decided (labelled or kept), from its history.
    pub decided: BTreeSet<ClusterId>,
    /// The pass cached its statistics, so marker edits can be rescored.
    pub rescorable: bool,
}

impl RoundView {
    pub fn load(manifest: &Path) -> Result<Self> {
        let loaded = run::load(&manifest.to_string_lossy())?;
        let cells = read_cells(&loaded)?;
        let a = &loaded.manifest.annotate;
        let at = |rel: &Option<String>| rel.as_deref().map(|r| resolve(&loaded.dir, r));

        let table = |rel: &Option<String>| at(rel).map(|p| read_table(&p)).transpose();
        // Q: how the cluster's evidence splits over the types; and the
        // statistics it came from.
        let share = table(&a.cluster_celltype_q)?;
        let (nes_table, p_table, q_table) = (
            table(&a.cluster_celltype_nes)?,
            table(&a.cluster_celltype_p)?,
            table(&a.cluster_celltype_q_values)?,
        );
        let expression = at(&a.cluster_expression)
            .map(|p| Mat::from_parquet_with_row_names(&p, Some(0)))
            .transpose()?
            .map(Expression::new);
        let markers = match at(&a.markers) {
            Some(p) if Path::new(&p).is_file() => panel_sets(&p)?,
            _ => BTreeMap::new(),
        };

        // Each cluster's size and label, as `review` and `relabel` see them.
        let digests = digest(&cells.clusters, &cells.labels, &Evidence::default());
        let loose_cells = cells.clusters.iter().filter(|c| c.is_none()).count();
        let decided = match at(&a.history) {
            Some(p) => decided_in(&p)?,
            None => BTreeSet::new(),
        };
        let genes = expression.as_ref().map(specific_genes).unwrap_or_default();
        let clusters = digests
            .into_iter()
            .map(|(id, d)| {
                let shares: Vec<(String, f32)> = share
                    .as_ref()
                    .and_then(|t| t.row(id).map(|r| (&t.cols, r)))
                    .map(|(cols, r)| cols.iter().cloned().zip(r.iter().copied()).collect())
                    .unwrap_or_default();
                let stat = |t: &Option<Table>, label: &str| {
                    let t = t.as_ref()?;
                    let j = t.cols.iter().position(|c| c == label)?;
                    t.row(id).map(|r| r[j]).filter(|v| v.is_finite())
                };
                let candidates = top_by(shares.iter().cloned(), CANDIDATES)
                    .into_iter()
                    .filter(|(l, s)| *s > 0.0 && l != UNASSIGNED_LABEL)
                    .map(|(label, share)| Candidate {
                        nes: stat(&nes_table, &label),
                        p: stat(&p_table, &label),
                        q: stat(&q_table, &label),
                        label,
                        share,
                    })
                    .collect();
                ClusterView {
                    id,
                    cells: d.size,
                    label: d.label,
                    candidates,
                    shares,
                    genes: genes.get(&id).cloned().unwrap_or_default(),
                }
            })
            .collect();
        Ok(Self {
            manifest: manifest.to_path_buf(),
            clusters,
            cell_names: cells.names,
            cell_clusters: cells.clusters,
            expression,
            markers,
            loose_cells,
            decided,
            rescorable: a.stats_cache.is_some() && a.expression_clusters.is_some(),
        })
    }

    /// Cluster `id`'s label with the edits applied; `None` for unassigned.
    #[must_use]
    pub fn label_of(&self, id: ClusterId, edits: &[Edit]) -> Option<String> {
        // The cluster's last decision: a label, or keeping the round's.
        let edited = edits.iter().rev().find(|e| e.cluster() == Some(id));
        match edited {
            Some(Edit::Label { label, .. }) => (label != UNASSIGNED_LABEL).then(|| label.clone()),
            _ => self
                .clusters
                .iter()
                .find(|c| c.id == id)
                .and_then(|c| c.label.clone()),
        }
    }

    /// Cells per label with the edits applied, most first; unassigned last.
    #[must_use]
    pub fn summary(&self, edits: &[Edit]) -> Vec<(String, usize)> {
        let mut n: BTreeMap<String, usize> = BTreeMap::new();
        for c in &self.clusters {
            let l = self.label_of(c.id, edits);
            *n.entry(l.unwrap_or_else(|| UNASSIGNED_LABEL.into()))
                .or_default() += c.cells;
        }
        if self.loose_cells > 0 {
            *n.entry(UNASSIGNED_LABEL.into()).or_default() += self.loose_cells;
        }
        let mut v: Vec<(String, usize)> = n.into_iter().collect();
        v.sort_by(|a, b| {
            (a.0 == UNASSIGNED_LABEL)
                .cmp(&(b.0 == UNASSIGNED_LABEL))
                .then(b.1.cmp(&a.1))
                .then_with(|| a.0.cmp(&b.0))
        });
        v
    }

    /// `gene`'s log2 fold change in cluster `id` over the other clusters'
    /// mean, as [`specific_genes`] scores it; `None` when not measured.
    #[must_use]
    pub fn fold_change(&self, id: ClusterId, gene: &str) -> Option<f32> {
        let e = self.expression.as_ref()?;
        e.log2fc(e.row(gene)?, &e.cols(&BTreeSet::from([id])))
    }

    /// `label`'s markers with the edits applied, each with whether an edit
    /// added it (dropped ones are gone), by name.
    #[must_use]
    pub fn markers_of(&self, label: &str, edits: &[Edit]) -> Vec<(String, bool)> {
        let key = label_key(label);
        let mut out: BTreeMap<String, bool> = self
            .markers
            .get(&key)
            .into_iter()
            .flatten()
            .map(|g| (g.clone(), false))
            .collect();
        for e in edits {
            if let Edit::Markers {
                label, genes, add, ..
            } = e
            {
                if label_key(label) != key {
                    continue;
                }
                for g in genes {
                    if *add {
                        out.entry(g.clone()).or_insert(true);
                    } else {
                        out.remove(g);
                    }
                }
            }
        }
        out.into_iter().collect()
    }

    /// The types `gene` is a marker of in the round's panel.
    #[must_use]
    pub fn marker_of(&self, gene: &str) -> Vec<&str> {
        self.markers
            .iter()
            .filter(|(_, g)| g.contains(gene))
            .map(|(l, _)| l.as_str())
            .collect()
    }
}

/// Per cluster, its candidates and shares, as [`RoundView::load`] reads them
/// or a rescoring gave them.
pub type Scores = BTreeMap<ClusterId, (Vec<Candidate>, Vec<(String, f32)>)>;

/// A `lupin relabel --preview`'s `scores` (per cluster id, each type's
/// share, NES, p and q) as candidates and shares; `None` when it has none.
#[must_use]
pub fn parse_scores(preview: &serde_json::Value) -> Option<Scores> {
    let num = |v: &serde_json::Value| v.as_f64().map(|x| x as f32);
    preview["scores"]
        .as_object()?
        .iter()
        .map(|(id, calls)| {
            let id: ClusterId = id.parse().ok()?;
            let all: Vec<Candidate> = calls
                .as_array()?
                .iter()
                .map(|c| {
                    Some(Candidate {
                        label: c["label"].as_str()?.to_string(),
                        share: num(&c["share"]).unwrap_or(0.0),
                        nes: num(&c["nes"]),
                        p: num(&c["p"]),
                        q: num(&c["q"]),
                    })
                })
                .collect::<Option<_>>()?;
            let shares = all.iter().map(|c| (c.label.clone(), c.share)).collect();
            let candidates = all
                .into_iter()
                .filter(|c| c.share > 0.0 && c.label != UNASSIGNED_LABEL)
                .take(CANDIDATES)
                .collect();
            Some((id, (candidates, shares)))
        })
        .collect()
}

impl RoundView {
    /// Show `scores` in place of the clusters' own; the ones replaced, to
    /// put back.
    pub fn swap_scores(&mut self, mut scores: Scores) -> Scores {
        let mut old = Scores::new();
        for c in &mut self.clusters {
            if let Some((candidates, shares)) = scores.remove(&c.id) {
                old.insert(
                    c.id,
                    (
                        std::mem::replace(&mut c.candidates, candidates),
                        std::mem::replace(&mut c.shares, shares),
                    ),
                );
            }
        }
        old
    }
}

/// Take `genes` out of the unsaved additions to types other than `label`,
/// so adding a gene again moves it rather than giving it a second type;
/// marker edits left empty are dropped. The first `from` edits are being
/// saved and stay as they are: a gene they add elsewhere is dropped from
/// that type by a new edit instead. The types the genes came from, by name.
pub fn take_added(
    edits: &mut Vec<Edit>,
    from: usize,
    genes: &[String],
    label: &str,
) -> Vec<String> {
    let key = label_key(label);
    let mut moved = BTreeSet::new();
    let mut drops: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for (i, e) in edits.iter_mut().enumerate() {
        let Edit::Markers {
            label: other,
            genes: added,
            add: true,
            ..
        } = e
        else {
            continue;
        };
        if label_key(other) == key {
            continue;
        }
        let hit: Vec<String> = added
            .iter()
            .filter(|g| genes.contains(g))
            .cloned()
            .collect();
        if hit.is_empty() {
            continue;
        }
        moved.insert(other.clone());
        if i < from {
            drops.entry(other.clone()).or_default().extend(hit);
        } else {
            added.retain(|g| !genes.contains(g));
        }
    }
    let from = from.min(edits.len());
    let mut tail = edits.split_off(from);
    tail.retain(|e| !matches!(e, Edit::Markers { genes, .. } if genes.is_empty()));
    edits.extend(tail);
    for (other, genes) in drops {
        edits.push(Edit::Markers {
            label: other,
            genes,
            add: false,
            reason: format!("moved to {label}"),
        });
    }
    moved.into_iter().collect()
}

/// The last edit per cluster, and every marker edit, as the next
/// round's decisions; marker edits first, as they change what the labels
/// are scored on.
#[must_use]
///
/// A label's evidence is the round's own candidates: `recorded`'s, when a
/// rescoring is shown in their place.
pub fn decisions(edits: &[Edit], round: &RoundView, recorded: Option<&Scores>) -> Vec<Decision> {
    let decision = |action, cluster, label, features, evidence, rationale: &str| Decision {
        cluster,
        action,
        label,
        features,
        evidence,
        alternatives: Vec::new(),
        rationale: rationale.to_string(),
        decided_by: DecidedBy::User,
        timestamp: None,
        round: None,
    };
    let mut last: BTreeMap<ClusterId, &Edit> = BTreeMap::new();
    let mut out = Vec::new();
    for e in edits {
        match e {
            Edit::Label { cluster, .. } | Edit::Keep { cluster, .. } => {
                last.insert(*cluster, e);
            }
            Edit::Markers {
                label,
                genes,
                add,
                reason,
            } => {
                let action = if *add {
                    Action::MarkersAdd
                } else {
                    Action::MarkersDrop
                };
                out.push(decision(
                    action,
                    Vec::new(),
                    Some(label.clone()),
                    genes.clone(),
                    Vec::new(),
                    reason,
                ));
            }
        }
    }
    for (id, e) in last {
        let evidence = recorded
            .and_then(|r| r.get(&id))
            .map(|(candidates, _)| &candidates[..])
            .or_else(|| {
                round
                    .clusters
                    .iter()
                    .find(|c| c.id == id)
                    .map(|c| &c.candidates[..])
            })
            .unwrap_or_default()
            .iter()
            .map(|c| {
                json!({"kind": "marker", "term": c.label, "share": c.share, "nes": c.nes, "p": c.p, "q": c.q})
            })
            .collect();
        let (action, label, reason) = match e {
            Edit::Label { label, reason, .. } => (Action::Label, Some(label.clone()), reason),
            Edit::Keep { reason, .. } => (Action::Keep, None, reason),
            Edit::Markers { .. } => continue,
        };
        out.push(decision(
            action,
            vec![id],
            label,
            Vec::new(),
            evidence,
            reason,
        ));
    }
    out
}

/// The clusters a round's history (`annotation_history.json`: cluster id →
/// decisions) records a decision for.
fn decided_in(path: &str) -> Result<BTreeSet<ClusterId>> {
    let text = std::fs::read_to_string(path)?;
    let history: BTreeMap<String, serde_json::Value> = serde_json::from_str(&text)?;
    Ok(history
        .into_iter()
        .filter(|(_, v)| v.as_array().is_some_and(|a| !a.is_empty()))
        .filter_map(|(k, _)| k.parse().ok())
        .collect())
}

fn top_by(items: impl Iterator<Item = (String, f32)>, k: usize) -> Vec<(String, f32)> {
    let mut v: Vec<(String, f32)> = items.filter(|(_, s)| s.is_finite()).collect();
    v.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    v.truncate(k);
    v
}

/// Per cluster, the genes that set it apart, all of them: the genes above
/// the cluster's average, by [`Expression::log2fc`] over the others.
#[must_use]
pub fn specific_genes(expr: &Expression) -> BTreeMap<ClusterId, Vec<(String, f32)>> {
    let m = &expr.table.mat;
    expr.col
        .iter()
        .map(|(&id, &k)| {
            let floor = m.column(k).mean();
            let scored = (0..m.nrows())
                .filter(|&g| m[(g, k)] > floor)
                .filter_map(|g| {
                    expr.log2fc(g, &[k])
                        .map(|fc| (expr.table.rows[g].to_string(), fc))
                });
            (id, top_by(scored, usize::MAX))
        })
        .collect()
}

/// A marker panel as label key → genes.
pub fn panel_sets(path: &str) -> Result<BTreeMap<String, BTreeSet<String>>> {
    let mut out: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for (g, t) in read_marker_pairs(path)? {
        out.entry(label_key(&t))
            .or_default()
            .insert(g.into_string());
    }
    Ok(out)
}

#[cfg(test)]
#[path = "tests/round.rs"]
mod tests;
