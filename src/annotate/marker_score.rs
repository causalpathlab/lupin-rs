//! The embedding marker score on a run's matrices: each cell's best panel
//! type and its margin over the runner-up, and each cluster's share of cells
//! called each type. The score is legume-graph-embedding's
//! `type_annotation::embedding_score`, a module score (Tirosh et al. 2016,
//! *Science* 352:189; SCANPY, Wolf et al. 2018, *Genome Biology* 19:15) read
//! in closed form on the run's feature and cell embeddings, centred at the
//! balanced centre of its clusters.

use anyhow::Result;
pub use enrichment::consensus::UNASSIGNED;
use graph_embedding_util::type_annotation::embedding_score::{
    balanced_centre, calls, cell_covariance, scores, FeatureSpace,
};
use legume_numeric::matrix::dense_mat_io::{Mat, MatWithNames};
use std::collections::{BTreeSet, HashMap};

/// Per cell: the index of its best type in `types` (`None` without a
/// finite score) and the margin over the runner-up (NaN without one).
pub struct MarkerScore {
    pub types: Vec<Box<str>>,
    pub calls: Vec<Option<usize>>,
    pub margins: Vec<f64>,
}

/// Score `cells` (`[N × H]`) against the `markers` panel on the feature
/// embedding `feat` (`[G × H]`, rows named by gene), centred at the
/// balanced centre of `groups` (cell → cluster, [`UNASSIGNED`] left out).
pub fn score_cells(
    feat: &MatWithNames<Mat>,
    cells: &Mat,
    groups: &[usize],
    markers: &str,
) -> Result<MarkerScore> {
    let rho = feat.mat.map(f64::from);
    let z = cells.map(f64::from);
    let space = FeatureSpace::new(&rho, &cell_covariance(&z)?, &feat.rows)?;
    let panel = space.marker_sets(markers, 0.0)?;
    let dirs = space.directions(&panel.sets);
    let s = scores(&z, &balanced_centre(&z, groups)?, &dirs);
    let (calls, margins) = calls(&s, &dirs)
        .into_iter()
        .map(|c| (c.call, c.margin))
        .unzip();
    Ok(MarkerScore {
        types: panel.type_names,
        calls,
        margins,
    })
}

/// The clusters (sorted ids, [`UNASSIGNED`] left out) and, per cluster, the
/// share of its cells called each of `n_types` types; an uncalled cell
/// counts toward its cluster's size.
#[must_use]
pub fn share_by_cluster(
    calls: &[Option<usize>],
    groups: &[usize],
    n_types: usize,
) -> (Vec<usize>, Mat) {
    let ids: Vec<usize> = groups
        .iter()
        .copied()
        .filter(|&k| k != UNASSIGNED)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let row: HashMap<usize, usize> = ids.iter().enumerate().map(|(r, &k)| (k, r)).collect();
    let mut share = Mat::zeros(ids.len(), n_types);
    let mut size = vec![0_u32; ids.len()];
    for (call, k) in calls.iter().zip(groups) {
        let Some(&r) = row.get(k) else { continue };
        size[r] += 1;
        if let Some(t) = call {
            share[(r, *t)] += 1.0;
        }
    }
    for (r, &n) in size.iter().enumerate() {
        share.row_mut(r).scale_mut(1.0 / n as f32);
    }
    (ids, share)
}

#[cfg(test)]
#[path = "marker_score_tests.rs"]
mod tests;
