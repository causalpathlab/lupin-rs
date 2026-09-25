//! Lineage family: velocity-oriented lineage inference over a `senna gem`
//! embedding (driven by [`run::run_lineage`]), plus `lineage_plot`, `assoc`
//! and `pseudotime`. Everything here takes plain file paths and loaded
//! matrices; [`crate::lineage_manifest`] resolves a `run.senna.json` into the
//! input structs.
//!
//! The generic numeric primitives live in `legume_numeric::matrix`: seeded k-means centroids
//! ([`legume_numeric::matrix::principal_graph::kmeans_centroids_seeded`]), the K×K distance matrix
//! ([`legume_numeric::matrix::principal_graph::pairwise_sqdist_rows_to_rows`]) + MST
//! ([`legume_numeric::matrix::principal_graph::mst_from_sqdist`]), and the Slingshot curves
//! ([`legume_numeric::matrix::principal_curve`]). This module holds two velocity-specific pieces on top of
//! that generic tree: the velocity [`orient`]ation of edges (δ from gem) and the
//! velocity [`forest`].

pub mod assoc;
pub mod lineage_plot;
pub mod pseudotime;

/// The `lupin lineage` command-line surface.
pub mod args;
mod cluster;
pub mod forest;
mod input;
mod layout;
pub mod orient;
mod root;
/// The `lupin lineage` run. Binary entry: [`run::run_lineage`].
pub mod run;
mod traj_annotation;
mod velocity_grid;
mod write;

/// What the producing run promised about its per-cell tables. Read from a run
/// manifest by [`crate::lineage_manifest`].
pub use input::LatentContract;
