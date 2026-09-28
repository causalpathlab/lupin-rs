//! The first round of a marker annotation calls clusters at a coarse level.
//!
//! Before the pass, the panel's cell types are grouped ([`TypeTree`]) and,
//! when the groups come from the Cell Ontology, the ontology and the
//! label → CL map are handed to the pass so its ontology walk runs by default.
//! After it, each cluster is called by the group its evidence adds up to;
//! the fine per-cell labels stay in `annotate.fine_argmax` and the fine calls
//! in the cluster summary, to refine in later rounds.

use crate::annotate::celltype_tree::{coarse_call, TypeTree};
use crate::annotate::rounds::ClusterId;
use crate::manifest::rounds::{read_cells, read_table, write_argmax, write_summary};
use crate::manifest::run::{self, rel_to_manifest, resolve, Loaded};
use anyhow::Result;
use legume_numeric::matrix::common_io::{mkdir_parent, write_lines};
use log::info;
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

const TREE: &str = ".celltype_tree.json";

const LABEL_CL: &str = ".label_cl.tsv";
const FINE_ARGMAX: &str = ".fine_argmax.tsv";

/// What [`prepare`] found for a marker pass.
pub struct Prepared {
    pub tree: TypeTree,
    /// The ontology file and the label → CL map to hand the pass, when the
    /// groups come from the Cell Ontology and the user named no map.
    pub ontology: Option<(String, String)>,
}

/// Whether `source` is a run no annotation round has been made from yet:
/// such a round is called at the coarse level.
#[must_use]
pub fn is_first_round(source: &Loaded) -> bool {
    let a = &source.manifest.annotate;
    a.celltype_tree.is_none() && a.history.is_none()
}

/// Group the panel at `markers` before a marker pass writing under `out`,
/// and write the groups to `{out}.celltype_tree.json`. `obo` and `label_cl`
/// are the user's; a label map is written for the pass when the groups come
/// from the Cell Ontology, every panel type matched a term, and the user
/// named no map.
pub fn prepare(
    markers: &str,
    out: &str,
    data: &super::data_files::ClData,
    label_cl: Option<&str>,
) -> Result<Prepared> {
    let panel = crate::annotate::markers::read_panel(markers)?;
    let tree = super::ontology::type_tree(data, &panel)?;
    let found = data.ontology.clone();
    mkdir_parent(out)?;
    fs::write(format!("{out}{TREE}"), serde_json::to_string_pretty(&tree)?)?;

    // The ontology walk needs every panel type on a term; with any left
    // unmatched it would refuse the whole pass, so it is not wired in then.
    let unmatched: Vec<&str> = tree
        .groups
        .iter()
        .flat_map(|g| g.members.iter())
        .filter(|m| !tree.label_cl.contains_key(*m))
        .map(String::as_str)
        .collect();
    if found.is_some() && label_cl.is_none() && !tree.label_cl.is_empty() && !unmatched.is_empty() {
        log::warn!(
            "{} panel type(s) match no Cell Ontology term ({}); the ontology walk is skipped. \
             Map them with --label-cl to run it",
            unmatched.len(),
            unmatched.join(", ")
        );
    }
    let ontology = match (found, label_cl) {
        (Some(obo), None) if !tree.label_cl.is_empty() && unmatched.is_empty() => {
            // Labels are already in the form the pass scores them under.
            let map = format!("{out}{LABEL_CL}");
            let lines: Vec<Box<str>> = tree
                .label_cl
                .iter()
                .map(|(label, id)| format!("{label}\t{id}").into_boxed_str())
                .collect();
            write_lines(&lines, &map)?;
            Some((obo.to_string_lossy().into_owned(), map))
        }
        _ => None,
    };
    Ok(Prepared { tree, ontology })
}

/// After the pass: record the groups and, when `coarse`, call each cluster
/// by the group its evidence adds up to, keeping the fine per-cell labels in
/// `annotate.fine_argmax`.
pub fn finish(manifest: &Path, tree: &TypeTree, coarse: bool) -> Result<()> {
    let mut loaded = run::load(&manifest.to_string_lossy())?;
    let out = run::derive_out_prefix(&loaded.file.to_string_lossy());
    loaded.manifest.annotate.celltype_tree =
        Some(rel_to_manifest(&loaded.dir, &format!("{out}{TREE}")));
    if !coarse {
        return loaded.manifest.save(&loaded.file);
    }
    let cells = read_cells(&loaded)?;
    // Per cluster × type, the pass's shares (enrichment): its Q probabilities.
    let probs = loaded
        .manifest
        .annotate
        .cluster_celltype_q
        .as_deref()
        .map(|rel| read_table(&resolve(&loaded.dir, rel)))
        .transpose()?;
    let min_share = 0.0;
    let mut votes: BTreeMap<ClusterId, Vec<&str>> = BTreeMap::new();
    for (id, label) in cells.clusters.iter().zip(&cells.labels) {
        if let Some(id) = id {
            let v = votes.entry(*id).or_default();
            v.extend(label.as_deref());
        }
    }
    let index = tree.index();
    let calls: BTreeMap<ClusterId, Option<String>> = votes
        .iter()
        .map(|(id, labels)| {
            let row = probs
                .as_ref()
                .and_then(|t| t.row(*id).map(|r| (&t.cols[..], r)));
            (*id, coarse_call(&index, row, labels, min_share))
        })
        .collect();

    let argmax_path = resolve(
        &loaded.dir,
        loaded
            .manifest
            .annotate
            .argmax
            .as_deref()
            .unwrap_or_default(),
    );
    let fine_path = format!("{out}{FINE_ARGMAX}");
    fs::copy(&argmax_path, &fine_path)?;
    // Enrichment calls whole clusters, so the cluster's group call is each
    // cell's; projection calls cells, so each keeps the group of its own fine
    // label, and a cell the pass abstained on stays unassigned.
    let labels: Vec<Option<String>> = if probs.is_some() {
        cells
            .clusters
            .iter()
            .map(|id| id.and_then(|id| calls.get(&id).cloned().flatten()))
            .collect()
    } else {
        cells
            .labels
            .iter()
            .map(|l| {
                l.as_deref()
                    .and_then(|l| index.get(&crate::annotate::markers::label_key(l)))
                    .map(|g| (*g).to_string())
            })
            .collect()
    };
    write_argmax(&argmax_path, &cells.names, &labels, &cells.probs)?;
    let named = calls.values().flatten().count();
    info!("first round: {named} cluster(s) called at the coarse level; fine labels in {fine_path}");

    loaded.manifest.annotate.fine_argmax = Some(rel_to_manifest(&loaded.dir, &fine_path));
    write_summary(&mut loaded, &out)?;
    loaded.manifest.save(&loaded.file)
}

#[cfg(test)]
#[path = "first_round_tests.rs"]
mod tests;
