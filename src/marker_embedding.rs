//! A run's gene table for marker matching, located through its manifest.
//!
//! Markers must be matched against genes placed on the cell manifold:
//! `outputs.feature_coembedding` for the embedding kinds, or ρ
//! (`outputs.feature_embedding`) for kinds whose ρ already shares the cells'
//! space. A `gem` table is keyed by feature row (spliced and unspliced per
//! gene), so only its spliced rows are kept, re-keyed by gene.

use anyhow::{Context, Result};
use data_beans::aux::feature_rows::split_count_row;
use legume_numeric::matrix::dense_mat_io::{Mat, MatWithNames};
use legume_numeric::matrix::traits::IoOps;
use log::info;
use std::path::Path;

use crate::run_manifest::{resolve, RunKind, RunManifest};

/// Load the marker-matching gene table. `prefix` only names the run in errors.
pub fn load_marker_feature_embedding(
    manifest: &RunManifest,
    dir: &Path,
    prefix: &str,
) -> Result<MatWithNames<Mat>> {
    let (slot, rel) = match (
        manifest.outputs.feature_coembedding.as_deref(),
        manifest.outputs.feature_embedding.as_deref(),
    ) {
        (Some(rel), _) => ("feature_coembedding", rel),
        (None, Some(_)) if manifest.kind.coembeds() => anyhow::bail!(
            "{prefix}: a {} run with no `outputs.feature_coembedding` — the co-embed was not \
             written (an interrupted run?), and its raw gene embedding ρ is not on the cell \
             manifold. Re-run the fit to completion.",
            manifest.kind
        ),
        (None, Some(rel)) => ("feature_embedding", rel),
        (None, None) => anyhow::bail!(
            "{prefix}: manifest has neither `outputs.feature_coembedding` nor \
             `outputs.feature_embedding` — this needs a gene embedding (a senna `gem` / `bge` / \
             `fne` / `resolve-embedding-space` run)"
        ),
    };
    let path = resolve(dir, rel);
    let feat = Mat::from_parquet(&path)
        .with_context(|| format!("reading gene embedding {path} (`outputs.{slot}`)"))?;
    if manifest.kind == RunKind::Gem {
        select_spliced_rows(feat, &path)
    } else {
        Ok(feat)
    }
}

fn select_spliced_rows(feat: MatWithNames<Mat>, path: &str) -> Result<MatWithNames<Mat>> {
    let (keep, rows): (Vec<usize>, Vec<Box<str>>) = feat
        .rows
        .iter()
        .enumerate()
        .filter_map(|(i, name)| match split_count_row(name) {
            Some((gene, false)) => Some((i, Box::from(gene))),
            _ => None,
        })
        .unzip();
    anyhow::ensure!(
        !keep.is_empty(),
        "{path} has no spliced count rows (found {} rows, e.g. `{}`). A spliced-only gem \
         run has no unspliced program to annotate.",
        feat.rows.len(),
        feat.rows.first().map_or("", |s| s.as_ref())
    );
    let mat = feat.mat.select_rows(&keep);
    info!(
        "gene embedding: {} of {} feature rows are spliced counts → {} genes [{} × {}]",
        keep.len(),
        feat.rows.len(),
        rows.len(),
        mat.nrows(),
        mat.ncols()
    );
    Ok(MatWithNames {
        mat,
        rows,
        cols: feat.cols,
    })
}
