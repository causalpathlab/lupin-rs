//! One annotation round as the screen works on it: each cluster's label,
//! the cell types competing for it, and the genes that set it apart; plus
//! the edits made here, which become the `relabel` decisions of the next
//! round when saved.

use crate::annotate::markers::{label_key, read_marker_pairs};
use crate::annotate::outputs::{
    CLUSTER_CELLTYPE_NES, CLUSTER_CELLTYPE_P, CLUSTER_CELLTYPE_Q_VALUES,
};
use crate::annotate::rounds::{ClusterId, Decision, Table};
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
/// Specific genes listed per cluster.
const TOP_GENES: usize = 20;

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
}

/// An edit made on screen, not yet saved.
#[derive(Clone, Debug, PartialEq)]
pub enum Edit {
    Label {
        cluster: ClusterId,
        label: String,
        reason: String,
    },
    /// Add `genes` to `label`'s markers, or drop them (`add` false).
    Markers {
        label: String,
        genes: Vec<String>,
        add: bool,
        reason: String,
    },
}

pub struct RoundView {
    pub manifest: PathBuf,
    pub clusters: Vec<ClusterView>,
    pub cell_names: Vec<Box<str>>,
    pub cell_clusters: Vec<Option<ClusterId>>,
    /// Genes × clusters, as the pass aggregated them.
    pub expression: Option<MatWithNames<Mat>>,
    /// Each gene's row of `expression`.
    pub gene_row: HashMap<String, usize>,
    /// The round's marker panel: label key → genes.
    pub markers: BTreeMap<String, BTreeSet<String>>,
    /// Cells no cluster holds.
    pub loose_cells: usize,
}

impl RoundView {
    pub fn load(manifest: &Path) -> Result<Self> {
        let loaded = run::load(&manifest.to_string_lossy())?;
        let cells = read_cells(&loaded)?;
        let a = &loaded.manifest.annotate;
        let at = |rel: &Option<String>| rel.as_deref().map(|r| resolve(&loaded.dir, r));

        // Q: how the cluster's evidence splits over the types; and beside
        // the q-values, the p-values and z they came from.
        let share = at(&a.cluster_celltype_q)
            .map(|p| read_table(&p))
            .transpose()?;
        let q_values = at(&a.cluster_celltype_q_values);
        let beside = |suffix: &str| {
            q_values
                .as_deref()
                .and_then(|q| q.strip_suffix(CLUSTER_CELLTYPE_Q_VALUES))
                .map(|stem| format!("{stem}{suffix}"))
                .filter(|p| Path::new(p).is_file())
                .map(|p| read_table(&p))
                .transpose()
        };
        let (nes_table, p_table) = (beside(CLUSTER_CELLTYPE_NES)?, beside(CLUSTER_CELLTYPE_P)?);
        let q_table = q_values.map(|p| read_table(&p)).transpose()?;
        let expression = at(&a.cluster_expression)
            .map(|p| Mat::from_parquet_with_row_names(&p, Some(0)))
            .transpose()?;
        let markers = match at(&a.markers) {
            Some(p) if Path::new(&p).is_file() => panel_sets(&p)?,
            _ => BTreeMap::new(),
        };

        let mut members: BTreeMap<ClusterId, Vec<usize>> = BTreeMap::new();
        let mut loose_cells = 0;
        for (i, c) in cells.clusters.iter().enumerate() {
            match c {
                Some(c) => members.entry(*c).or_default().push(i),
                None => loose_cells += 1,
            }
        }
        let genes = expression.as_ref().map(specific_genes).unwrap_or_default();
        let clusters = members
            .into_iter()
            .map(|(id, idx)| {
                let labels = idx.iter().filter_map(|&i| cells.labels[i].as_deref());
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
                    cells: idx.len(),
                    label: majority(labels),
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
            gene_row: expression
                .as_ref()
                .map(|e| {
                    e.rows
                        .iter()
                        .enumerate()
                        .map(|(r, g)| (g.to_string(), r))
                        .collect()
                })
                .unwrap_or_default(),
            expression,
            markers,
            loose_cells,
        })
    }

    /// Cluster `id`'s label with the edits applied; `None` for unassigned.
    #[must_use]
    pub fn label_of(&self, id: ClusterId, edits: &[Edit]) -> Option<String> {
        let edited = edits.iter().rev().find_map(|e| match e {
            Edit::Label { cluster, label, .. } if *cluster == id => Some(label.clone()),
            _ => None,
        });
        match edited {
            Some(l) => (l != UNASSIGNED_LABEL).then_some(l),
            None => self
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
        let r = *self.gene_row.get(gene)?;
        let n = e.cols.len();
        let k =
            (0..n).find(|&k| crate::annotate::rounds::parse_cluster_id(&e.cols[k]) == Some(id))?;
        if n < 2 {
            return None;
        }
        let eps = e.mat.mean().max(f32::MIN_POSITIVE);
        let x = e.mat[(r, k)];
        let rest = (e.mat.row(r).sum() - x) / (n - 1) as f32;
        Some(((x + eps) / (rest + eps)).log2())
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

/// The last edit per cluster, and every marker edit, as the next
/// round's decisions.
pub fn decisions(edits: &[Edit], round: &RoundView) -> Result<Vec<Decision>> {
    let mut last: BTreeMap<ClusterId, &Edit> = BTreeMap::new();
    let mut out = Vec::new();
    for e in edits {
        match e {
            Edit::Label { cluster, .. } => {
                last.insert(*cluster, e);
            }
            Edit::Markers {
                label,
                genes,
                add,
                reason,
            } => out.push(json!({
                "action": if *add { "markers_add" } else { "markers_drop" },
                "label": label,
                "features": genes,
                "rationale": reason,
                "decided_by": "user",
            })),
        }
    }
    for (id, e) in last {
        let Edit::Label { label, reason, .. } = e else {
            continue;
        };
        let c = round.clusters.iter().find(|c| c.id == id);
        let evidence: Vec<_> = c
            .map(|c| &c.candidates[..])
            .unwrap_or_default()
            .iter()
            .map(|c| json!({"kind": "marker", "term": c.label, "share": c.share, "nes": c.nes, "p": c.p, "q": c.q}))
            .collect();
        out.push(json!({
            "cluster": id,
            "action": "label",
            "label": label,
            "evidence": evidence,
            "rationale": reason,
            "decided_by": "user",
        }));
    }
    // Marker edits first, as they change what the labels are scored on.
    out.into_iter()
        .map(|v| Ok(serde_json::from_value(v)?))
        .collect()
}

/// The most frequent of `labels`; ties to the smallest; `None` when empty.
fn majority<'a>(labels: impl Iterator<Item = &'a str>) -> Option<String> {
    let mut n: BTreeMap<&str, usize> = BTreeMap::new();
    for l in labels {
        *n.entry(l).or_default() += 1;
    }
    n.into_iter()
        .max_by(|a, b| a.1.cmp(&b.1).then_with(|| b.0.cmp(a.0)))
        .map(|(l, _)| l.to_string())
}

fn top_by(items: impl Iterator<Item = (String, f32)>, k: usize) -> Vec<(String, f32)> {
    let mut v: Vec<(String, f32)> = items.filter(|(_, s)| s.is_finite()).collect();
    v.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    v.truncate(k);
    v
}

/// Per cluster, the genes that set it apart, from a genes × clusters
/// expression table (columns `K{id}`): among the genes above the cluster's
/// average, the largest log2 fold change over the other clusters' mean.
pub fn specific_genes(expr: &MatWithNames<Mat>) -> BTreeMap<ClusterId, Vec<(String, f32)>> {
    let m = &expr.mat;
    let (n_genes, n_clusters) = (m.nrows(), m.ncols());
    if n_genes == 0 || n_clusters < 2 {
        return BTreeMap::new();
    }
    let totals: Vec<f32> = (0..n_genes).map(|g| m.row(g).sum()).collect();
    // A pseudo-count on the scale of the table, so a gene near zero
    // elsewhere does not win on a vanishing denominator.
    let eps = m.mean().max(f32::MIN_POSITIVE);
    (0..n_clusters)
        .filter_map(|k| {
            let id = crate::annotate::rounds::parse_cluster_id(&expr.cols[k])?;
            let floor = m.column(k).mean();
            let scored = (0..n_genes).filter(|&g| m[(g, k)] > floor).map(|g| {
                let x = m[(g, k)];
                let rest = (totals[g] - x) / (n_clusters - 1) as f32;
                (expr.rows[g].to_string(), ((x + eps) / (rest + eps)).log2())
            });
            Some((id, top_by(scored, TOP_GENES)))
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
