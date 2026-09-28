//! The Cell Ontology for annotation, with the rules and aliases that place
//! a panel on it, found by [`super::data_files`] for a run. With no ontology
//! at hand (offline, say), a panel's types are grouped by the marker genes
//! they share instead.

use super::data_files::{ClData, SearchPath};
use crate::annotate::celltype_tree::{ClTerms, TypeTree};
use crate::annotate::panel_tree::PanelTree;
use anyhow::Result;
use log::{info, warn};
use std::path::{Path, PathBuf};

/// The ontology data for a run whose manifest sits in `run_dir`: `obo`
/// (`--obo`) and `label_cl` (`--label-cl`) are the run's own files, the most
/// specific layer.
pub fn load(run_dir: Option<&Path>, obo: Option<&str>, label_cl: Option<&str>) -> Result<ClData> {
    ClData::load(SearchPath::new(run_dir), obo, label_cl)
}

/// The coarse level over `panel`: on the Cell Ontology when there is one and
/// at least two of the panel's types match its terms, else by shared markers.
/// Also the ontology file used, if any.
pub fn type_tree(data: &ClData, panel: &[(String, String)]) -> Result<(TypeTree, Option<PathBuf>)> {
    if let Some(terms) = data.terms()? {
        if let Some(tree) = TypeTree::from_ontology(&terms, panel) {
            info!(
                "cell-type groups from the Cell Ontology ({}): {} group(s), {} type(s) matched",
                terms.release.as_deref().unwrap_or("release unknown"),
                tree.groups.len(),
                tree.label_cl.len()
            );
            return Ok((tree, data.ontology.clone()));
        }
        warn!("fewer than two panel types match a Cell Ontology term; grouping by shared markers");
    }
    let tree = TypeTree::from_markers(panel);
    info!(
        "cell-type groups by shared markers: {} group(s)",
        tree.groups.len()
    );
    Ok((tree, data.ontology.clone()))
}

/// The panel's types on the Cell Ontology, else grouped by the markers they
/// share: the tree a round's labels are picked from and its q-values
/// adjusted over.
pub fn panel_tree(data: &ClData, panel: &[(String, String)]) -> Result<PanelTree> {
    Ok(panel_tree_on(data.terms()?.as_ref(), panel))
}

/// [`panel_tree`] on an ontology already read.
#[must_use]
pub fn panel_tree_on(terms: Option<&ClTerms>, panel: &[(String, String)]) -> PanelTree {
    terms
        .and_then(|t| PanelTree::from_ontology(t, panel))
        .unwrap_or_else(|| PanelTree::from_type_tree(&TypeTree::from_markers(panel)))
}
