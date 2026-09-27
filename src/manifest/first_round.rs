//! The first round of a marker annotation calls clusters at a coarse level.
//!
//! Before the pass, the panel's cell types are grouped ([`TypeTree`]) and,
//! when the groups come from the Cell Ontology, the ontology and the
//! label → CL map are handed to the pass so its ontology walk runs by default.
//! After it, each cluster is called by the group its evidence adds up to;
//! the fine per-cell labels stay in `annotate.fine_argmax` and the fine calls
//! in the cluster summary, to refine in later rounds.

use crate::annotate::celltype_tree::{coarse_call, TreeSource, TypeTree};
use crate::manifest::rounds::{
    read_argmax, read_clusters, read_table, write_argmax, write_summary,
};
use crate::manifest::run::{self, rel_to_manifest, resolve};
use anyhow::{Context, Result};
use legume_numeric::matrix::common_io::{mkdir_parent, write_lines};
use log::info;
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

pub const TREE: &str = ".celltype_tree.json";
pub const LABEL_CL: &str = ".label_cl.tsv";
pub const FINE_ARGMAX: &str = ".fine_argmax.tsv";

/// What [`prepare`] found: the groups, written at `{out}.celltype_tree.json`,
/// and the ontology and label map the pass should use, when there are any.
pub struct Prepared {
    pub tree: TypeTree,
    pub tree_path: String,
    pub obo: Option<String>,
    pub label_cl: Option<String>,
}

/// Group the panel at `markers` before a marker pass writing under `out`.
/// `obo` and `label_cl` are the user's; when absent and the Cell Ontology is
/// found, they are filled in.
pub fn prepare(
    markers: &str,
    out: &str,
    obo: Option<&str>,
    label_cl: Option<&str>,
) -> Result<Prepared> {
    let panel: Vec<(String, String)> = crate::annotate::markers::read_marker_pairs(markers)?
        .into_iter()
        .map(|(g, t)| (g.into_string(), t.into_string()))
        .collect();
    let (tree, found) = super::ontology::type_tree(obo, &panel)?;
    mkdir_parent(out)?;
    let tree_path = format!("{out}{TREE}");
    fs::write(&tree_path, serde_json::to_string_pretty(&tree)?)?;

    let obo = found.map(|p| p.to_string_lossy().into_owned());
    let label_cl = match label_cl {
        Some(l) => Some(l.to_string()),
        None if tree.source == TreeSource::CellOntology => {
            // Labels are already in the form the pass scores them under.
            let path = format!("{out}{LABEL_CL}");
            let lines: Vec<Box<str>> = tree
                .label_cl
                .iter()
                .map(|(label, id)| format!("{label}\t{id}").into_boxed_str())
                .collect();
            write_lines(&lines, &path)?;
            Some(path)
        }
        None => None,
    };
    Ok(Prepared {
        tree,
        tree_path,
        obo: obo.filter(|_| label_cl.is_some()),
        label_cl,
    })
}

/// After the pass: record the groups and, unless `fine`, turn the round's
/// labels into each cluster's coarse call, keeping the fine ones beside them.
pub fn finish(manifest: &Path, prepared: &Prepared, fine: bool) -> Result<()> {
    let mut loaded = run::load(&manifest.to_string_lossy())?;
    let out = run::derive_out_prefix(&loaded.file.to_string_lossy());
    loaded.manifest.annotate.celltype_tree =
        Some(rel_to_manifest(&loaded.dir, &prepared.tree_path));
    if fine {
        return loaded.manifest.save(&loaded.file);
    }
    let a = &loaded.manifest.annotate;
    let argmax_rel = a
        .argmax
        .clone()
        .context("the pass wrote no per-cell labels")?;
    let clusters_rel = loaded
        .manifest
        .cluster
        .clusters
        .clone()
        .context("the pass recorded no clusters")?;
    let (cells, ids) = read_clusters(&resolve(&loaded.dir, &clusters_rel))?;
    let argmax_path = resolve(&loaded.dir, &argmax_rel);
    let fine_labels = read_argmax(&argmax_path)?;

    // Per cluster: its probabilities over the panel's types, when the pass
    // wrote them (enrichment), and its cells' fine labels.
    let probs: BTreeMap<u32, Vec<(String, f32)>> = match a.cluster_celltype_q.as_deref() {
        Some(rel) => {
            let t = read_table(&resolve(&loaded.dir, rel))?;
            let w = t.cols.len();
            t.rows
                .iter()
                .enumerate()
                .map(|(r, id)| {
                    let row = t
                        .cols
                        .iter()
                        .cloned()
                        .zip(t.values[r * w..(r + 1) * w].iter().copied());
                    (*id, row.collect())
                })
                .collect()
        }
        None => BTreeMap::new(),
    };
    let mut votes: BTreeMap<u32, Vec<&str>> = BTreeMap::new();
    for (cell, id) in cells.iter().zip(&ids) {
        if let (Some(id), Some((label, _))) = (id, fine_labels.get(cell.as_ref())) {
            votes.entry(*id).or_default().push(label.as_str());
        }
    }
    let calls: BTreeMap<u32, Option<String>> = votes
        .keys()
        .chain(probs.keys())
        .map(|id| {
            let call = coarse_call(
                &prepared.tree,
                probs.get(id).map(Vec::as_slice),
                votes.get(id).map_or(&[][..], Vec::as_slice),
            );
            (*id, call)
        })
        .collect();

    let fine_path = format!("{out}{FINE_ARGMAX}");
    fs::copy(&argmax_path, &fine_path)?;
    let (labels, probs_out): (Vec<Option<String>>, Vec<f32>) = cells
        .iter()
        .zip(&ids)
        .map(|(cell, id)| {
            let p = fine_labels.get(cell.as_ref()).map_or(f32::NAN, |(_, p)| *p);
            (id.and_then(|id| calls.get(&id).cloned().flatten()), p)
        })
        .unzip();
    write_argmax(&argmax_path, &cells, &labels, &probs_out)?;
    let named = calls.values().flatten().count();
    info!("first round: {named} cluster(s) called at the coarse level; fine labels in {fine_path}");

    loaded.manifest.annotate.fine_argmax = Some(rel_to_manifest(&loaded.dir, &fine_path));
    write_summary(&mut loaded, &out)?;
    loaded.manifest.save(&loaded.file)
}

#[cfg(test)]
#[path = "first_round_tests.rs"]
mod tests;
