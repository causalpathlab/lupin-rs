//! A round's genes × clusters expression profile and its marker panel, read
//! for the genes behind a cluster's calls: the genes that set it apart, and
//! how a type's markers fare in it. The TUI's gene views and `lupin mcp`'s
//! `annotate_genes` both read them here.

use crate::annotate::markers::{label_key, read_marker_pairs};
use crate::annotate::rounds::{parse_cluster_id, ClusterId};
use crate::manifest::run::{resolve, Loaded};
use anyhow::{Context, Result};
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::traits::{IoOps, MatWithNames};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;

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

    /// Gene `gene`'s log2 fold change in cluster `id` over the others;
    /// `None` when either is not in the table.
    #[must_use]
    pub fn fold_change(&self, id: ClusterId, gene: &str) -> Option<f32> {
        self.log2fc(self.row(gene)?, &self.cols(&BTreeSet::from([id])))
    }
}

/// `items` with a finite value, largest first, ties by name.
fn top_by(items: impl Iterator<Item = (String, f32)>) -> Vec<(String, f32)> {
    let mut v: Vec<(String, f32)> = items.filter(|(_, s)| s.is_finite()).collect();
    v.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
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
            (id, top_by(scored))
        })
        .collect()
}

/// A marker panel as label key → genes.
pub type Panel = BTreeMap<String, BTreeSet<String>>;

/// A marker panel as label key → genes, read from `path`.
pub fn panel_sets(path: &str) -> Result<Panel> {
    let mut out = Panel::new();
    for (g, t) in read_marker_pairs(path)? {
        out.entry(label_key(&t))
            .or_default()
            .insert(g.into_string());
    }
    Ok(out)
}

/// The expression profile `loaded` records, and its marker panel as label
/// key → genes (empty when it records none, or the file is gone).
pub fn profile_and_panel(loaded: &Loaded) -> Result<(Option<Expression>, Panel)> {
    let a = &loaded.manifest.annotate;
    let at = |rel: &Option<String>| rel.as_deref().map(|r| resolve(&loaded.dir, r));
    let expression = at(&a.cluster_expression)
        .map(|p| {
            Mat::from_parquet_with_row_names(&p, Some(0)).with_context(|| format!("reading {p}"))
        })
        .transpose()?
        .map(Expression::new);
    let markers = match at(&a.markers) {
        Some(p) if Path::new(&p).is_file() => panel_sets(&p)?,
        _ => BTreeMap::new(),
    };
    Ok((expression, markers))
}
