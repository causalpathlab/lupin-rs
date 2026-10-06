//! An enrichment pass whose raw counts are not on this machine (a run trained
//! elsewhere and copied without its data). The cluster expression then comes
//! from what the run left beside it:
//!
//! 1. an earlier pass's cached statistics ([`super::recalibrate::write_cache`]):
//!    its clusters and their gene sums, exact, so the pass keeps those
//!    clusters whatever clustering it was asked for;
//! 2. else the model's decoder, at each cluster's cells: `β · θ̄` for the
//!    topic kinds and `masked-vae`, the mean of `softmax(z·W + b)` for `vae`.
//!    This is the model's expected expression, not the data's.
//!
//! Pseudobulk counts are not used: a pseudobulk mixes the clusters a cell
//! type is told apart by.

use super::annotate::{panel_inputs, resolve_clusters};
use super::family::family;
use super::rounds::read_clusters;
use super::run::{self, resolve, Loaded, RunKind};
use crate::annotate::aggregate::weighted_mean_profile;
use crate::annotate::args::AnnotateArgs;
use crate::annotate::inputs::EnrichmentInputs;
use crate::annotate::rounds::parse_cluster_id;
use anyhow::{Context, Result};
use legume_numeric::matrix::dense_mat_io::{read_mat, Mat};
use legume_numeric::matrix::traits::IoOps;
use log::{info, warn};
use rayon::prelude::*;
use serde_json::json;
use std::path::{Path, PathBuf};

/// Counts per cell the decoder's expected proportions are scaled to, so the
/// gene sums read like counts.
const NOMINAL_DEPTH: f64 = 1e4;

/// What a stand-in for the counts gives: the per-cluster and per-batch
/// expression over the cells, and where it came from.
struct Expression {
    gene_names: Vec<Box<str>>,
    cell_names: Vec<Box<str>>,
    cluster_labels: Vec<usize>,
    n_clusters: usize,
    batch_labels: Vec<usize>,
    n_batches: usize,
    /// Row-major `k · g` gene sums.
    gene_sum_kg: Vec<f64>,
    pb_gene_gp: Mat,
    gene_weights: Vec<f32>,
    source: serde_json::Value,
}

/// [`EnrichmentInputs`] for a run whose count files `missing` are not here.
pub(super) fn inputs(
    args: &AnnotateArgs,
    loaded: &Loaded,
    data: Option<&super::data_files::ClData>,
    missing: &[&str],
) -> Result<EnrichmentInputs> {
    warn!(
        "{} of the run's count file(s) are not here (first: {}); the cluster expression comes from what the run left beside it",
        missing.len(),
        missing[0]
    );
    let e = match from_cache(args, loaded)? {
        Some(e) => e,
        None => from_decoder(args, loaded)?.with_context(|| {
            format!(
                "the run's raw counts are not here ({} file(s) missing, first {}), no earlier pass \
                 cached its cluster sums beside it, and a {} run has no decoder to read expression from",
                missing.len(),
                missing[0],
                loaded.manifest.kind.as_str()
            )
        })?,
    };
    anyhow::ensure!(
        e.n_clusters >= 2,
        "annotate needs ≥ 2 clusters, found {}",
        e.n_clusters
    );
    let g = e.gene_names.len();
    let panel = panel_inputs(args, loaded, data, &e.gene_names)?;
    let profile_gk = weighted_mean_profile(&e.gene_sum_kg, e.n_clusters, g, &e.gene_weights);
    Ok(EnrichmentInputs {
        gene_names: e.gene_names,
        cell_names: e.cell_names,
        cluster_labels: e.cluster_labels,
        n_clusters: e.n_clusters,
        batch_labels: e.batch_labels,
        n_batches: e.n_batches,
        markers_gc: panel.markers_gc,
        celltype_names: panel.celltype_names,
        profile_gk,
        pb_gene_gp: e.pb_gene_gp,
        gene_sum_kg: e.gene_sum_kg,
        gene_weights: e.gene_weights,
        type_tree: panel.type_tree,
        cl_record: panel.cl_record,
        expression_source: Some(e.source),
    })
}

/// The run's manifest first, then the rest of its family.
fn family_manifests(loaded: &Loaded) -> Vec<PathBuf> {
    let same = |a: &Path, b: &Path| match (a.canonicalize(), b.canonicalize()) {
        (Ok(a), Ok(b)) => a == b,
        _ => a == b,
    };
    let mut out = vec![loaded.file.clone()];
    for m in family(&loaded.file) {
        if !out.iter().any(|p| same(p, &m.pick.path)) {
            out.push(m.pick.path);
        }
    }
    out
}

/// Where a pass on a run gets its cluster expression.
pub enum Source {
    /// The raw counts are here.
    Counts,
    /// Missing: an earlier pass's cache, this manifest's.
    Cache(PathBuf),
    /// Missing, with no cache: the decoder of a run of this kind.
    Decoder(&'static str),
    /// Missing, and nothing stands in.
    Nothing,
}

/// Where a pass on `loaded` gets its cluster expression.
pub fn source(loaded: &Loaded) -> Source {
    let files = loaded.manifest.data_inputs(&loaded.dir);
    if files.iter().all(|f| Path::new(f).exists()) {
        return Source::Counts;
    }
    if let Some(src) = newest_cache(loaded) {
        return Source::Cache(src.file);
    }
    match decoder_of(loaded) {
        Some(Decoder::Mixture) => Source::Decoder("β · θ̄"),
        Some(Decoder::Softmax) => Source::Decoder("mean softmax(z·W + b)"),
        None => Source::Nothing,
    }
}

/// The newest manifest in the run's family whose pass cached its
/// statistics, with every cached file here.
fn newest_cache(loaded: &Loaded) -> Option<Loaded> {
    let mut found: Vec<(std::time::SystemTime, Loaded)> = family_manifests(loaded)
        .into_iter()
        .filter_map(|p| run::load(&p.to_string_lossy()).ok())
        .filter(|src| {
            let a = &src.manifest.annotate;
            let (Some(c), Some(ids)) = (&a.stats_cache, &a.expression_clusters) else {
                return false;
            };
            [
                &c.gene_sum,
                &c.gene_weight,
                &c.batch_profile,
                &c.cell_batch,
                ids,
            ]
            .iter()
            .all(|r| Path::new(&resolve(&src.dir, r)).is_file())
        })
        .map(|src| {
            let t = src.file.metadata().and_then(|m| m.modified());
            (t.unwrap_or(std::time::UNIX_EPOCH), src)
        })
        .collect();
    found.sort_by_key(|(t, _)| std::cmp::Reverse(*t));
    found.into_iter().next().map(|(_, src)| src)
}

/// The newest cached statistics in the run's family, with their clusters.
fn from_cache(args: &AnnotateArgs, loaded: &Loaded) -> Result<Option<Expression>> {
    let Some(src) = newest_cache(loaded) else {
        return Ok(None);
    };
    let a = &src.manifest.annotate;
    let (Some(cache), Some(ids_rel)) = (&a.stats_cache, &a.expression_clusters) else {
        return Ok(None);
    };
    let at = |rel: &str| resolve(&src.dir, rel);
    if args.clusters.is_some() || loaded.manifest.cluster.clusters.is_none() {
        warn!(
            "without the counts the pass keeps the clusters {} summed; the clustering asked for is not used",
            src.file.display()
        );
    }

    // The cached clusters, numbered 0.. in id order.
    let (cell_names, ids) = read_clusters(&at(ids_rel))?;
    let mut order: Vec<u32> = ids.iter().flatten().copied().collect();
    order.sort_unstable();
    order.dedup();
    let slot = |id: u32| order.binary_search(&id).ok();
    let cluster_labels: Vec<usize> = ids
        .iter()
        .map(|id| id.and_then(slot).unwrap_or(usize::MAX))
        .collect();
    let k = order.len();

    let sums = read_mat(&at(&cache.gene_sum))?;
    let g = sums.rows.len();
    let mut gene_sum_kg = vec![0f64; g * k];
    for (j, col) in sums.cols.iter().enumerate() {
        let Some(dest) = parse_cluster_id(col).and_then(slot) else {
            continue;
        };
        for i in 0..g {
            gene_sum_kg[dest * g + i] = f64::from(sums.mat[(i, j)]);
        }
    }
    let weights = read_mat(&at(&cache.gene_weight))?;
    anyhow::ensure!(
        weights.rows == sums.rows,
        "{}: its genes are not the gene sums'",
        cache.gene_weight
    );
    let gene_weights: Vec<f32> = (0..g).map(|i| weights.mat[(i, 0)]).collect();
    let pb = read_mat(&at(&cache.batch_profile))?;
    let (batch_cells, batch_ids) = read_clusters(&at(&cache.cell_batch))?;
    let batch_of: std::collections::HashMap<&str, usize> = batch_cells
        .iter()
        .zip(&batch_ids)
        .filter_map(|(n, b)| b.map(|b| (n.as_ref(), b as usize)))
        .collect();
    let batch_labels = cell_names
        .iter()
        .map(|n| batch_of.get(n.as_ref()).copied().unwrap_or(0))
        .collect();
    info!(
        "cluster expression: the gene sums {} cached over its {k} clusters",
        src.file.display()
    );
    Ok(Some(Expression {
        gene_names: sums.rows,
        cell_names,
        cluster_labels,
        n_clusters: k,
        batch_labels,
        n_batches: pb.mat.ncols(),
        gene_sum_kg,
        pb_gene_gp: pb.mat,
        gene_weights,
        source: json!({
            "from": "cache",
            "manifest": src.file.file_name().map(|n| n.to_string_lossy().into_owned()),
        }),
    }))
}

/// How a run's decoder turns a cell's latent row into a gene distribution.
enum Decoder {
    /// `π = β · θ`: `θ = softmax(latent)` over the topics (`log θ` or raw
    /// `z`), `β` a gene distribution per topic.
    Mixture,
    /// `π = softmax_d(z · W + b)` over the decoder's features (scVI).
    Softmax,
}

/// The decoder of `loaded`'s kind, when it wrote a latent and a dictionary.
fn decoder_of(loaded: &Loaded) -> Option<Decoder> {
    let o = &loaded.manifest.outputs;
    o.latent.as_ref()?;
    o.softmax_dictionary.as_ref().or(o.dictionary.as_ref())?;
    match loaded.manifest.kind {
        RunKind::Topic | RunKind::Itopic | RunKind::JointTopic | RunKind::MaskedVae => {
            Some(Decoder::Mixture)
        }
        RunKind::Vae => Some(Decoder::Softmax),
        _ => None,
    }
}

/// The run's expected expression per cluster and batch, from its decoder.
fn from_decoder(args: &AnnotateArgs, loaded: &Loaded) -> Result<Option<Expression>> {
    let m = &loaded.manifest;
    let Some(decoder) = decoder_of(loaded) else {
        return Ok(None);
    };
    let o = &m.outputs;
    let (Some(latent_rel), Some(dict_rel)) = (
        o.latent.as_deref(),
        o.softmax_dictionary.as_deref().or(o.dictionary.as_deref()),
    ) else {
        return Ok(None);
    };
    let at = |rel: &str| resolve(&loaded.dir, rel);
    let latent = Mat::from_parquet_with_row_names(&at(latent_rel), Some(0))
        .map_err(|e| anyhow::anyhow!("reading {latent_rel}: {e}"))?;
    let dict = Mat::from_parquet_with_row_names(&at(dict_rel), Some(0))
        .map_err(|e| anyhow::anyhow!("reading {dict_rel}: {e}"))?;
    anyhow::ensure!(
        latent.mat.ncols() == dict.mat.ncols(),
        "the latent has {} factors but the dictionary {}",
        latent.mat.ncols(),
        dict.mat.ncols()
    );
    let cell_names = latent.rows;
    let (cluster_labels, k) = resolve_clusters(args, m, &loaded.dir, &cell_names)?;
    let (batch_labels, n_batches) = batch_labels(loaded, &cell_names);
    let g = dict.rows.len();

    // Expected gene distributions summed over each cluster's and batch's
    // cells, row-major `groups · g`, scaled to counts.
    let groups = Groups {
        clusters: &cluster_labels,
        k,
        batches: &batch_labels,
        n_batches,
    };
    let (sum_kg, sum_pg) = match decoder {
        Decoder::Mixture => {
            let beta = softmax_columns(&dict.mat);
            let (tk, tp) = theta_sums(&latent.mat, &groups);
            (expected(&beta, &tk), expected(&beta, &tp))
        }
        Decoder::Softmax => {
            let (features, gene_of) = vae_features(loaded, &dict.rows, &dict.mat)?;
            let bias = vae_bias(loaded, features.nrows())?;
            let (fk, fp) = decoded_sums(&latent.mat, &features, &bias, &groups);
            (to_genes(&fk, &gene_of), to_genes(&fp, &gene_of))
        }
    };
    let gene_sum_kg = sum_kg;
    let gene_weights = fisher_weights(loaded, &dict.rows);
    let pb_gene_gp = weighted_mean_profile(&sum_pg, n_batches, g, &gene_weights);
    info!(
        "cluster expression: the {} decoder's expected expression at each cluster's cells (no counts)",
        m.kind.as_str()
    );
    Ok(Some(Expression {
        gene_names: dict.rows,
        cell_names,
        cluster_labels,
        n_clusters: k,
        batch_labels,
        n_batches,
        gene_sum_kg,
        pb_gene_gp,
        gene_weights,
        source: json!({ "from": "decoder", "kind": m.kind.as_str(), "depth": NOMINAL_DEPTH }),
    }))
}

/// Cells per cluster (and per batch) the `vae` decoder is run at: their
/// mean stands for the group's.
const DECODED_CELLS: usize = 1000;

/// Each cell's cluster (`usize::MAX` for none) and batch.
struct Groups<'a> {
    clusters: &'a [usize],
    k: usize,
    batches: &'a [usize],
    n_batches: usize,
}

/// `θ = softmax(latent row)` summed over each cluster's and batch's cells:
/// `K × k` and `K × p`.
fn theta_sums(latent: &Mat, groups: &Groups) -> (Mat, Mat) {
    let n_topics = latent.ncols();
    let mut tk = Mat::zeros(n_topics, groups.k);
    let mut tp = Mat::zeros(n_topics, groups.n_batches);
    for (n, row) in latent.row_iter().enumerate() {
        let z: Vec<f32> = row.iter().copied().collect();
        for (t, th) in softmax(&z).into_iter().enumerate() {
            if groups.clusters[n] < groups.k {
                tk[(t, groups.clusters[n])] += th as f32;
            }
            tp[(t, groups.batches[n])] += th as f32;
        }
    }
    (tk, tp)
}

/// `β · θsums` (`g × K` by `K × groups`) as row-major `groups · g`, at
/// [`NOMINAL_DEPTH`] counts per cell.
fn expected(beta: &Mat, theta_sums: &Mat) -> Vec<f64> {
    let e = beta * theta_sums;
    e.column_iter()
        .flat_map(|c| {
            c.iter()
                .map(|&v| NOMINAL_DEPTH * f64::from(v))
                .collect::<Vec<_>>()
        })
        .collect()
}

/// Up to [`DECODED_CELLS`] of each group's cells, evenly spaced.
fn spaced(labels: &[usize], n_groups: usize) -> Vec<Vec<usize>> {
    let mut members = vec![Vec::new(); n_groups];
    for (n, &l) in labels.iter().enumerate() {
        if l < n_groups {
            members[l].push(n);
        }
    }
    members
        .into_iter()
        .map(|m| {
            let take = m.len().min(DECODED_CELLS);
            (0..take).map(|i| m[i * m.len() / take]).collect()
        })
        .collect()
}

/// `softmax_d(z·W + b)` over each cluster's and batch's cells, as sums over
/// all their cells (`groups · d`), estimated from [`spaced`] cells.
fn decoded_sums(
    latent: &Mat,
    features: &Mat,
    bias: &[f32],
    groups: &Groups,
) -> (Vec<f64>, Vec<f64>) {
    let d = features.nrows();
    let sums = |labels: &[usize], n_groups: usize| -> Vec<f64> {
        let picked = spaced(labels, n_groups);
        let size = |c: usize| labels.iter().filter(|&&l| l == c).count() as f64;
        let wt = features.transpose();
        let mut out = vec![0f64; n_groups * d];
        for (c, cells) in picked.iter().enumerate() {
            if cells.is_empty() {
                continue;
            }
            let row_sum = cells
                .par_chunks(256)
                .map(|chunk| {
                    let z = Mat::from_fn(chunk.len(), latent.ncols(), |r, t| latent[(chunk[r], t)]);
                    let logits = &z * &wt;
                    let mut acc = vec![0f64; d];
                    for row in logits.row_iter() {
                        let l: Vec<f32> = row.iter().zip(bias).map(|(v, b)| v + b).collect();
                        acc.iter_mut().zip(softmax(&l)).for_each(|(a, p)| *a += p);
                    }
                    acc
                })
                .reduce(
                    || vec![0f64; d],
                    |mut a, b| {
                        a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
                        a
                    },
                );
            let scale = NOMINAL_DEPTH * size(c) / cells.len() as f64;
            for (o, v) in out[c * d..][..d].iter_mut().zip(row_sum) {
                *o = v * scale;
            }
        }
        out
    };
    (
        sums(groups.clusters, groups.k),
        sums(groups.batches, groups.n_batches),
    )
}

/// Feature sums (`groups · d`) to gene sums (`groups · g`) by each gene's
/// feature and share of it.
fn to_genes(sum: &[f64], gene_of: &[(usize, f64)]) -> Vec<f64> {
    let g = gene_of.len();
    let d = gene_of.iter().map(|&(f, _)| f + 1).max().unwrap_or(0);
    let groups = sum.len().checked_div(d).unwrap_or(0);
    let mut out = vec![0f64; groups * g];
    for c in 0..groups {
        for (i, &(f, share)) in gene_of.iter().enumerate() {
            out[c * g + i] = sum[c * d + f] * share;
        }
    }
    out
}

fn softmax(x: &[f32]) -> Vec<f64> {
    let max = x.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let e: Vec<f64> = x.iter().map(|&v| f64::from(v - max).exp()).collect();
    let s: f64 = e.iter().sum();
    e.into_iter().map(|v| v / s).collect()
}

/// Each column of `m` (log weights or logits over genes) as a distribution.
fn softmax_columns(m: &Mat) -> Mat {
    let mut out = m.clone();
    for mut col in out.column_iter_mut() {
        let p = softmax(col.as_slice());
        col.iter_mut().zip(p).for_each(|(c, v)| *c = v as f32);
    }
    out
}

/// A `vae`'s decoder features: its loadings per feature, and each gene's
/// feature with the gene's share of it. A run that coarsened its genes into
/// modules decodes the modules, and its dictionary repeats a module's
/// loadings for each of its genes; a gene's share of its module is its mean
/// rate's share (equal shares without `feature_mean`).
fn vae_features(
    loaded: &Loaded,
    genes: &[Box<str>],
    dict: &Mat,
) -> Result<(Mat, Vec<(usize, f64)>)> {
    let prefix = loaded.run_prefix();
    let g = genes.len();
    let modules = std::fs::read_to_string(format!("{prefix}.coarsening.json"))
        .ok()
        .map(|text| -> Result<Vec<usize>> {
            let v: serde_json::Value = serde_json::from_str(&text)?;
            let finest = v["levels"]
                .as_array()
                .and_then(|l| l.last())
                .context("coarsening.json has no levels")?;
            let map: Vec<usize> = serde_json::from_value(finest["fine_to_coarse"].clone())?;
            anyhow::ensure!(
                map.len() == g,
                "coarsening covers {} genes, the dictionary {g}",
                map.len()
            );
            Ok(map)
        })
        .transpose()?;
    let Some(module_of) = modules else {
        return Ok((dict.clone(), (0..g).map(|i| (i, 1.0)).collect()));
    };
    let d = module_of.iter().max().map_or(0, |m| m + 1);
    let mut features = Mat::zeros(d, dict.ncols());
    for (i, &f) in module_of.iter().enumerate() {
        features.row_mut(f).copy_from(&dict.row(i));
    }
    let mean = read_mat(&format!("{prefix}.feature_mean.parquet"))
        .ok()
        .filter(|m| m.rows == genes)
        .map(|m| {
            (0..g)
                .map(|i| f64::from(m.mat[(i, 0)]).max(0.0))
                .collect::<Vec<_>>()
        });
    let mut total = vec![0f64; d];
    let mut size = vec![0f64; d];
    for (i, &f) in module_of.iter().enumerate() {
        total[f] += mean.as_ref().map_or(0.0, |m| m[i]);
        size[f] += 1.0;
    }
    let gene_of = module_of
        .iter()
        .enumerate()
        .map(|(i, &f)| {
            let share = match &mean {
                Some(m) if total[f] > 0.0 => m[i] / total[f],
                _ => 1.0 / size[f],
            };
            (f, share)
        })
        .collect();
    Ok((features, gene_of))
}

/// The `vae` decoder's per-feature offset `b`, which only its checkpoint
/// holds: the finest level's `gauss_decoder.bias` (`d` long).
fn vae_bias(loaded: &Loaded, d: usize) -> Result<Vec<f32>> {
    let path = loaded
        .manifest
        .outputs
        .extra
        .get("model")
        .and_then(|v| v.as_str())
        .map(|r| resolve(&loaded.dir, r))
        .unwrap_or_else(|| format!("{}.safetensors", loaded.run_prefix()));
    let tensors = candle_core::safetensors::load(&path, &candle_core::Device::Cpu)
        .with_context(|| format!("reading the vae decoder {path}"))?;
    let bias = tensors
        .iter()
        .filter(|(name, t)| name.ends_with("gauss_decoder.bias") && t.dims() == [d])
        .map(|(_, t)| t)
        .next()
        .with_context(|| format!("{path} has no decoder bias over {d} features"))?;
    Ok(bias.to_vec1::<f32>()?)
}

/// The cells' batches: the run's batch files, when they are here and cover
/// the cells in the latent's order; else the `@sample` the cell names carry;
/// else one batch.
fn batch_labels(loaded: &Loaded, cells: &[Box<str>]) -> (Vec<usize>, usize) {
    let files = loaded.manifest.data_batches(&loaded.dir);
    let from_files: Option<Vec<Box<str>>> = (!files.is_empty())
        .then(|| {
            files
                .iter()
                .map(|f| legume_numeric::matrix::common_io::read_lines(f).ok())
                .collect::<Option<Vec<_>>>()
        })
        .flatten()
        .map(|v| v.into_iter().flatten().collect());
    let labels: Vec<&str> = match &from_files {
        Some(l) if l.len() == cells.len() => l.iter().map(AsRef::as_ref).collect(),
        _ => cells
            .iter()
            .map(|c| c.rsplit_once('@').map_or("", |(_, s)| s))
            .collect(),
    };
    data_beans::alg::dc_poisson::compact_labels(&labels)
}

/// NB-Fisher weights cached at training, when they cover these genes; else
/// equal weights.
fn fisher_weights(loaded: &Loaded, genes: &[Box<str>]) -> Vec<f32> {
    match data_beans::alg::gene_weighting::load_fisher_weights(&loaded.run_prefix()) {
        Ok(Some((names, w))) if names == genes => w,
        _ => {
            info!("no cached NB-Fisher weights for these genes: equal weights");
            vec![1.0; genes.len()]
        }
    }
}

#[cfg(test)]
#[path = "tests/no_counts.rs"]
mod tests;
