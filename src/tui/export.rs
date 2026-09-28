//! A round's final annotation, for use outside lupin:
//!
//! - `{round}.cell_annotation.parquet`: per cell, its cluster, label, the
//!   label's CL term and lineage (every level above it, top first), the
//!   cluster's evidence share for the label, and the cluster's top type;
//! - `{round}.celltype_markers.tsv`: per label, its marker genes, each from
//!   the original `panel` or `added` since, and the `suggested` genes that
//!   set its clusters apart; each gene with its log2 fold change in the
//!   label's clusters over the rest.

use super::round::{Edit, RoundView};
use crate::annotate::celltype_tree::ClTerms;
use crate::annotate::markers::label_key;
use crate::annotate::panel_tree::PanelTree;
use crate::annotate::rounds::ClusterId;
use anyhow::Result;
use enrichment::UNASSIGNED_LABEL;
use legume_numeric::matrix::common_io::write_lines;
use legume_numeric::matrix::parquet::{write_named_table, Column};
use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

/// Suggested genes listed per label.
const SUGGESTED: usize = 10;

/// A cluster's final call.
#[derive(Debug, Clone, PartialEq)]
pub struct Call {
    pub label: String,
    pub cl_id: String,
    pub lineage: String,
    pub share: f32,
    pub top: String,
}

/// Each cluster's call, the edits applied. A label off the panel's tree is
/// looked up in the Cell Ontology `cl`, when given, for its term, lineage and
/// the share of the panel types under it.
#[must_use]
pub fn cluster_calls(
    round: &RoundView,
    tree: &PanelTree,
    cl: Option<&ClTerms>,
    edits: &[Edit],
) -> BTreeMap<ClusterId, Call> {
    round
        .clusters
        .iter()
        .map(|c| {
            let label = round
                .label_of(c.id, edits)
                .unwrap_or_else(|| UNASSIGNED_LABEL.into());
            let node = (label != UNASSIGNED_LABEL)
                .then(|| tree.node_of(&label))
                .flatten();
            // Off the panel's tree: the CL term the label names, if any.
            let term = node
                .is_none()
                .then(|| cl.zip(Some(&label)))
                .flatten()
                .filter(|(_, l)| l.as_str() != UNASSIGNED_LABEL)
                .and_then(|(cl, l)| super::ontology::term_of(cl, tree, l).map(|t| (cl, t)));
            let under: BTreeSet<String> = match (node, &term) {
                (Some(i), _) => tree.labels_under(i).into_iter().map(String::from).collect(),
                (None, Some((cl, t))) => tree
                    .typed_terms()
                    .filter(|(_, id)| cl.ancestors_or_self(id).contains(t))
                    .map(|(l, _)| l.to_string())
                    .chain([label_key(&label)])
                    .collect(),
                (None, None) => BTreeSet::from([label_key(&label)]),
            };
            let share = c
                .shares
                .iter()
                .filter(|(t, _)| under.contains(&label_key(t)))
                .map(|(_, s)| s)
                .sum();
            let lineage = node.map_or_else(
                || match &term {
                    Some((cl, t)) => {
                        let mut names: Vec<String> = cl
                            .lineage(t)
                            .iter()
                            .map(|id| label_key(cl.name(id).unwrap_or(id)))
                            .collect();
                        names.pop();
                        names.push(label.clone());
                        names.join(" > ")
                    }
                    None => label.clone(),
                },
                |i| {
                    let path = tree.path(i);
                    let mut names: Vec<&str> = path[..path.len() - 1]
                        .iter()
                        .map(|&n| tree.nodes[n].name.as_str())
                        .collect();
                    names.push(&label);
                    names.join(" > ")
                },
            );
            let call = Call {
                cl_id: node
                    .and_then(|i| tree.nodes[i].cl_id.clone())
                    .or_else(|| term.as_ref().map(|(_, t)| t.clone()))
                    .unwrap_or_default(),
                lineage,
                share,
                top: c
                    .candidates
                    .first()
                    .map(|t| t.label.clone())
                    .unwrap_or_default(),
                label,
            };
            (c.id, call)
        })
        .collect()
}

/// One row of the marker table.
#[derive(Debug, Clone, PartialEq)]
pub struct MarkerRow {
    pub celltype: String,
    pub gene: String,
    /// `panel`, `added` or `suggested`.
    pub source: &'static str,
    pub log2fc: f32,
}

/// Per final label: its markers (the round's panel for it and every type
/// under it on the tree), from `original` or added since, and the genes
/// most specific to its clusters that are not markers yet.
#[must_use]
pub fn marker_rows(
    round: &RoundView,
    tree: &PanelTree,
    calls: &BTreeMap<ClusterId, Call>,
    original: &BTreeMap<String, BTreeSet<String>>,
) -> Vec<MarkerRow> {
    let mut clusters_of: BTreeMap<&str, BTreeSet<ClusterId>> = BTreeMap::new();
    for (id, c) in calls {
        if c.label != UNASSIGNED_LABEL {
            clusters_of.entry(c.label.as_str()).or_default().insert(*id);
        }
    }
    let expr = round.expression.as_ref();
    let mut out = Vec::new();
    for (label, ids) in clusters_of {
        let types: Vec<String> = match tree.node_of(label) {
            Some(i) => tree.labels_under(i).into_iter().map(String::from).collect(),
            None => vec![label_key(label)],
        };
        let in_panel = |set: &BTreeMap<String, BTreeSet<String>>, g: &str| {
            types
                .iter()
                .any(|t| set.get(t).is_some_and(|s| s.contains(g)))
        };
        let genes: BTreeSet<&str> = types
            .iter()
            .filter_map(|t| round.markers.get(t))
            .flatten()
            .map(String::as_str)
            .collect();
        let inside = expr.map(|e| e.cols(&ids)).unwrap_or_default();
        let fc_of = |g: &str| {
            expr.and_then(|e| e.log2fc(e.row(g)?, &inside))
                .unwrap_or(f32::NAN)
        };
        let mut rows: Vec<MarkerRow> = genes
            .iter()
            .map(|g| MarkerRow {
                celltype: label.to_string(),
                gene: (*g).to_string(),
                source: if in_panel(original, g) {
                    "panel"
                } else {
                    "added"
                },
                log2fc: fc_of(g),
            })
            .collect();
        rows.sort_by(|a, b| {
            a.source
                .cmp(b.source)
                .reverse()
                .then(b.log2fc.total_cmp(&a.log2fc))
        });
        if let Some(e) = expr {
            let names = &e.table.rows;
            // A suggestion must be above the label's clusters' average.
            let floor = (0..names.len()).map(|r| e.mean_in(r, &inside)).sum::<f32>()
                / names.len().max(1) as f32;
            let mut suggested: Vec<(usize, f32)> = (0..names.len())
                .filter(|&r| !genes.contains(&*names[r]) && e.mean_in(r, &inside) > floor)
                .filter_map(|r| e.log2fc(r, &inside).map(|f| (r, f)))
                .filter(|(_, f)| f.is_finite() && *f > 0.0)
                .collect();
            suggested.sort_by(|a, b| b.1.total_cmp(&a.1));
            rows.extend(
                suggested
                    .into_iter()
                    .take(SUGGESTED)
                    .map(|(r, f)| MarkerRow {
                        celltype: label.to_string(),
                        gene: names[r].to_string(),
                        source: "suggested",
                        log2fc: f,
                    }),
            );
        }
        out.extend(rows);
    }
    out
}

/// Write both files next to the round; returns their paths.
pub fn write(
    round: &RoundView,
    tree: &PanelTree,
    cl: Option<&ClTerms>,
    original: &BTreeMap<String, BTreeSet<String>>,
) -> Result<(PathBuf, PathBuf)> {
    let stem = crate::manifest::run::derive_out_prefix(&round.manifest.to_string_lossy());
    let calls = cluster_calls(round, tree, cl, &[]);

    let n = round.cell_names.len();
    let (mut cluster, mut label, mut cl_id, mut lineage, mut top) = (
        Vec::with_capacity(n),
        Vec::with_capacity(n),
        Vec::with_capacity(n),
        Vec::with_capacity(n),
        Vec::with_capacity(n),
    );
    let mut share = Vec::with_capacity(n);
    for c in &round.cell_clusters {
        let call = c.and_then(|c| calls.get(&c));
        cluster.push(c.map_or(-1, i64::from));
        let text = |f: fn(&Call) -> &str, none: &str| -> Box<str> { call.map_or(none, f).into() };
        label.push(text(|c| &c.label, UNASSIGNED_LABEL));
        cl_id.push(text(|c| &c.cl_id, ""));
        lineage.push(text(|c| &c.lineage, UNASSIGNED_LABEL));
        top.push(text(|c| &c.top, ""));
        share.push(call.map_or(f32::NAN, |c| c.share));
    }
    let cells = PathBuf::from(format!("{stem}.cell_annotation.parquet"));
    write_named_table(
        &cells.to_string_lossy(),
        "cell",
        &round.cell_names,
        &[
            ("cluster".into(), Column::I64(&cluster)),
            ("label".into(), Column::Str(&label)),
            ("cl_id".into(), Column::Str(&cl_id)),
            ("lineage".into(), Column::Str(&lineage)),
            ("share".into(), Column::F32(&share)),
            ("top_type".into(), Column::Str(&top)),
        ],
    )?;

    let markers = PathBuf::from(format!("{stem}.celltype_markers.tsv"));
    let mut lines: Vec<Box<str>> = vec!["celltype\tgene\tsource\tlog2fc".into()];
    lines.extend(
        marker_rows(round, tree, &calls, original)
            .into_iter()
            .map(|r| format!("{}\t{}\t{}\t{:.3}", r.celltype, r.gene, r.source, r.log2fc).into()),
    );
    write_lines(&lines, &markers.to_string_lossy())?;
    Ok((cells, markers))
}

#[cfg(test)]
#[path = "tests/export.rs"]
mod tests;
