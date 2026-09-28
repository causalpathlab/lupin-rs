//! Marker-based cell-type calls attached to the trajectory's nodes.

use anyhow::{Context, Result};
use legume_numeric::matrix::dense_mat_io::axis_id_names;
use log::info;

use graph_embedding_util::type_annotation::{
    annotate_with_communities, CommunityCalls, InputEmbeddings, TermOraConfig,
};
use legume_numeric::matrix::branching::Branching;
use legume_numeric::matrix::dmatrix_io::DMatrix;
use legume_numeric::matrix::parquet::{write_named_table, Column};
use legume_numeric::matrix::traits::MatWithNames;

/// Inputs for [`annotate_trajectory`] — bundled to keep the fan-in a struct.
pub(super) struct AnnotateTrajArgs<'a> {
    /// The co-embedded gene vectors `[G × H]` the node calls score against.
    /// Loading these needs the producing run's manifest, so the caller hands
    /// them in already loaded — see [`crate::manifest::lineage`].
    pub(super) feature_embedding: &'a MatWithNames<DMatrix<f32>>,
    pub(super) out: &'a str,
    pub(super) markers: &'a str,
    /// Raw θ `[N × H]` — the same latent space `lupin annotate --method projection` scores in.
    pub(super) raw_theta: &'a DMatrix<f32>,
    pub(super) cell_names: &'a [Box<str>],
    /// Per-cell MST-node id (the k-means `labels`) — the annotation clustering.
    pub(super) labels: &'a [usize],
    /// Number of MST nodes.
    pub(super) k: usize,
    pub(super) num_perm: usize,
    pub(super) obo: Option<&'a str>,
    pub(super) label_cl: Option<&'a str>,
    pub(super) seed: u64,
}

/// Name each trajectory node by cell type: run the `lupin annotate --method projection` term-ORA core over the
/// MST-node grouping (raw θ vs the gem β dictionary), giving every node a permutation-
/// calibrated call. Writes `{out}.lineage_annot.*` and returns the per-node
/// [`CommunityCalls`]. Run BEFORE root selection (it doesn't depend on the root) so
/// `--root-type` can pick the root from these calls; the caller writes
/// `{out}.trajectory_annotation.parquet` afterwards via [`write_trajectory_annotation`].
pub(super) fn compute_node_calls(a: &AnnotateTrajArgs) -> Result<CommunityCalls> {
    // The co-embedded feature vectors, not β — see `senna::marker_embedding` for why a
    // Euclidean nearest-centroid call against β is not a well-posed question.
    let beta = a.feature_embedding;
    let cfg = TermOraConfig {
        n_perm: a.num_perm,
        // `--seed` drives the whole fit; it should drive the annotation's randomness too. It was
        // silently falling through to `TermOraConfig`'s own default of 42, so varying `--seed`
        // moved the centroids but left the permutation null untouched.
        min_markers: 3,
        seed: a.seed,
        obo: a.obo.map(str::to_owned),
        label_cl: a.label_cl.map(str::to_owned),
        panel_perm: 0,
        support_perm: 0,
        bootstrap: None,
        ..TermOraConfig::default()
    };
    let input = InputEmbeddings {
        feature_emb: &beta.mat,
        gene_names: &beta.rows,
        cell_emb: a.raw_theta,
        cell_names: a.cell_names,
    };

    annotate_with_communities(
        &input,
        a.markers,
        &format!("{}.lineage_annot", a.out),
        true, // IDF-weight markers, as `lupin annotate --method projection` does by default
        a.labels,
        a.k,
        None,
        &cfg,
    )
}

/// Cross the per-node calls with the rooted forest → the labeled trajectory: one row per
/// node — `role` (root | terminal | internal), `cell_type`, `confidence`. Terminals are
/// derived from the rooted children (a node with no children), not from the orientation,
/// so abstained edges cannot misclassify a leaf.
pub(super) fn write_trajectory_annotation(
    calls: &CommunityCalls,
    br: &Branching,
    path: &str,
) -> Result<()> {
    let k = br.parent.len();
    let mut has_child = vec![false; k];
    for p in br.parent.iter().flatten() {
        if *p < k {
            has_child[*p] = true;
        }
    }
    let node_names = axis_id_names("node_", k);
    let roles: Vec<Box<str>> = (0..k)
        .map(|node| {
            if br.parent[node].is_none() {
                Box::from("root")
            } else if !has_child[node] {
                Box::from("terminal")
            } else {
                Box::from("internal")
            }
        })
        .collect();
    write_named_table(
        path,
        "node",
        &node_names,
        &[
            (Box::from("role"), Column::Str(&roles)),
            (Box::from("cell_type"), Column::Str(&calls.labels)),
            (Box::from("confidence"), Column::F32(&calls.confidence)),
        ],
    )
    .with_context(|| format!("writing {path}"))?;
    info!("wrote {path} ({k} nodes; {} root(s))", br.roots.len());
    Ok(())
}

/////////////////////
// PHATE 2D layout //
/////////////////////
