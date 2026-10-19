//! A run's gene table for projection, located through its manifest
//! ([`projection_gene_table`]). The embedding marker score reads the decoder
//! table ρ (`outputs.feature_embedding`) of a kind that decodes cells
//! through it. The gene-set units are trained on genes placed on the cell
//! manifold: `outputs.feature_coembedding` for the kinds that
//! co-embed, or ρ for kinds whose ρ already shares the cells' space. A `gem`
//! table is keyed by feature row (spliced and unspliced per gene), so only
//! its spliced rows are kept, re-keyed by gene.

use anyhow::{Context, Result};
use data_beans::aux::feature_rows::split_count_row;
use legume_numeric::matrix::dense_mat_io::{Mat, MatWithNames};
use legume_numeric::matrix::traits::IoOps;
use log::info;
use std::path::Path;

use crate::annotate::args::ProjectionMethod;
use crate::manifest::run::{resolve, RunKind, RunManifest};

/// The manifest slot and relative path of the gene table `method` reads, or
/// why the run has none. The score reads the decoder table ρ
/// (`outputs.feature_embedding`) of a kind that decodes cells through it; the
/// units read the co-embedding, or the ρ of a kind whose ρ shares the cells'
/// space (a kind that co-embeds needs its co-embedding: its ρ is not on the
/// cell manifold).
pub fn projection_gene_table(
    manifest: &RunManifest,
    method: ProjectionMethod,
) -> Result<(&'static str, &str)> {
    let o = &manifest.outputs;
    let kind = &manifest.kind;
    match method {
        ProjectionMethod::Score => {
            anyhow::ensure!(
                kind.coembeds(),
                "a {kind} run has no decoder gene table to score markers on (bge, gem, \
                 simba and resolve-embedding-space runs do); use `lupin annotate --method \
                 enrichment`"
            );
            let rel = o.feature_embedding.as_deref().ok_or_else(|| {
                anyhow::anyhow!(
                    "a {kind} run with no `outputs.feature_embedding` (an interrupted fit?)"
                )
            })?;
            Ok(("feature_embedding", rel))
        }
        // Units align the cells with sets through the genes placed among
        // them.
        ProjectionMethod::Units => match (
            o.feature_coembedding.as_deref(),
            o.feature_embedding.as_deref(),
        ) {
            (Some(rel), _) => Ok(("feature_coembedding", rel)),
            (None, Some(_)) if kind.coembeds() => anyhow::bail!(
                "a {kind} run with no `outputs.feature_coembedding`: the co-embed was not \
                 written (an interrupted run?), and its raw gene embedding ρ is not on the \
                 cell manifold. Re-run the fit to completion, or use `--projection score`."
            ),
            (None, Some(rel)) => Ok(("feature_embedding", rel)),
            (None, None) => anyhow::bail!(
                "manifest has neither `outputs.feature_coembedding` nor \
                 `outputs.feature_embedding` — this needs a gene embedding (a senna `gem` / \
                 `bge` / `fne` / `resolve-embedding-space` run)"
            ),
        },
    }
}

/// The gene table the units are trained on: the genes placed among the cells. `prefix` only names the
/// run in errors.
pub fn load_marker_feature_embedding(
    manifest: &RunManifest,
    dir: &Path,
    prefix: &str,
) -> Result<MatWithNames<Mat>> {
    read_gene_table(manifest, dir, prefix, ProjectionMethod::Units)
}

/// The run's decoder gene table ρ, which the embedding marker score reads:
/// the table the model multiplies a cell's embedding by, `log rate = ρ_g · z
/// + a_g`. Never the co-embedding. `prefix` only names the run in errors.
pub fn load_feature_embedding(
    manifest: &RunManifest,
    dir: &Path,
    prefix: &str,
) -> Result<MatWithNames<Mat>> {
    read_gene_table(manifest, dir, prefix, ProjectionMethod::Score)
}

/// Read the gene table `method` reads; a table keyed by count row (`gem`, a
/// splice-channelized pinto run) keeps its spliced rows, re-keyed by gene.
fn read_gene_table(
    manifest: &RunManifest,
    dir: &Path,
    prefix: &str,
    method: ProjectionMethod,
) -> Result<MatWithNames<Mat>> {
    let (slot, rel) =
        projection_gene_table(manifest, method).map_err(|e| anyhow::anyhow!("{prefix}: {e}"))?;
    let path = resolve(dir, rel);
    let feat = Mat::from_parquet(&path)
        .with_context(|| format!("reading gene embedding {path} (`outputs.{slot}`)"))?;
    if manifest.kind == RunKind::Gem || feat.rows.iter().any(|r| split_count_row(r).is_some()) {
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

#[cfg(test)]
#[path = "marker_embedding_tests.rs"]
mod tests;
