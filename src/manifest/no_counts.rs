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
use super::recalibrate::{complete_cache, read_cache};
use super::rounds::read_clusters;
use super::run::{self, resolve, Loaded, RunKind};
use crate::annotate::aggregate::weighted_mean_profile;
use crate::annotate::args::AnnotateArgs;
use crate::annotate::inputs::EnrichmentInputs;
use anyhow::{Context, Result};
use legume_numeric::matrix::dense_mat_io::{read_mat, Mat};
use legume_numeric::matrix::traits::{IoOps, MatOps};
use log::{info, warn};
use rayon::prelude::*;
use serde_json::json;
use std::path::Path;

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

/// Where a pass on a run gets its cluster expression.
pub enum Source {
    /// The raw counts are here.
    Counts,
    /// Missing: the cache this manifest's pass wrote.
    Cache(Box<Loaded>),
    /// Missing, with no cache: the decoder of the run's kind.
    Decoder(Decoder),
    /// Missing, and nothing stands in.
    Nothing,
}

/// Where a pass on `loaded` gets its cluster expression.
pub fn source(loaded: &Loaded) -> Source {
    if missing_counts(loaded).is_empty() {
        return Source::Counts;
    }
    if let Some(src) = newest_cache(loaded) {
        return Source::Cache(Box::new(src));
    }
    match decoder_of(loaded) {
        Some((d, _, _)) => Source::Decoder(d),
        None => Source::Nothing,
    }
}

/// The run's count files that are not here.
fn missing_counts(loaded: &Loaded) -> Vec<String> {
    let mut files = loaded.manifest.data_inputs(&loaded.dir);
    files.retain(|f| !Path::new(f).exists());
    files
}

/// [`EnrichmentInputs`] from `source`, a stand-in for the counts
/// [`self::source`] chose.
pub(super) fn inputs(
    args: &AnnotateArgs,
    loaded: &Loaded,
    data: Option<&super::data_files::ClData>,
    source: Source,
) -> Result<EnrichmentInputs> {
    let missing = missing_counts(loaded);
    let first = missing.first().map_or("", String::as_str);
    warn!(
        "{} of the run's count file(s) are not here (first: {first}); the cluster expression comes from what the run left beside it",
        missing.len()
    );
    let e = match source {
        Source::Cache(src) => from_cache(args, loaded, &src)?,
        Source::Decoder(_) => from_decoder(args, loaded)?,
        Source::Nothing => anyhow::bail!(
            "the run's raw counts are not here ({} file(s) missing, first {first}), no earlier pass \
             cached its cluster sums beside it, and a {} run has no decoder to read expression from",
            missing.len(),
            loaded.manifest.kind.as_str()
        ),
        Source::Counts => anyhow::bail!("the run's raw counts are here: read them instead"),
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

/// The newest manifest in the run's family whose pass cached its
/// statistics, with every cached file here.
fn newest_cache(loaded: &Loaded) -> Option<Loaded> {
    let mut members = family(&loaded.file);
    members.sort_by_key(|m| std::cmp::Reverse(m.modified));
    members
        .into_iter()
        .filter_map(|m| run::load(&m.pick.path.to_string_lossy()).ok())
        .find(|src| complete_cache(src).is_some())
}

/// The cached statistics `src`'s pass wrote, with its clusters.
fn from_cache(args: &AnnotateArgs, loaded: &Loaded, src: &Loaded) -> Result<Expression> {
    let (cache, ids_rel) = complete_cache(src)
        .with_context(|| format!("{}: its cache is gone", src.file.display()))?;
    // The cached clusters, numbered 0.. in id order.
    let (cell_names, ids) = read_clusters(&resolve(&src.dir, ids_rel))?;
    if !asks_for(args, loaded, &cell_names, &ids) {
        warn!(
            "without the counts the pass keeps the clusters {} summed, not the clustering asked for",
            src.file.display()
        );
    }
    let mut order: Vec<u32> = ids.iter().flatten().copied().collect();
    order.sort_unstable();
    order.dedup();
    let slot = |id: u32| order.binary_search(&id).ok();
    let cluster_labels: Vec<usize> = ids
        .iter()
        .map(|id| id.and_then(slot).unwrap_or(usize::MAX))
        .collect();
    let k = order.len();
    let stats = read_cache(src, cache, slot, k, &cell_names)?;
    info!(
        "cluster expression: the gene sums {} cached over its {k} clusters",
        src.file.display()
    );
    Ok(Expression {
        gene_names: stats.gene_names,
        cell_names,
        cluster_labels,
        n_clusters: k,
        batch_labels: stats.batch_labels,
        n_batches: stats.pb_gene_gp.ncols(),
        gene_sum_kg: stats.gene_sum_kg,
        pb_gene_gp: stats.pb_gene_gp,
        gene_weights: stats.gene_weights,
        source: json!({
            "from": "cache",
            "manifest": src.file.file_name().map(|n| n.to_string_lossy().into_owned()),
        }),
    })
}

/// Whether the cached clusters `ids` (of `cells`) are the ones the pass was
/// asked for: the run's own cluster file, with no `--clusters`. Leiden (no
/// cluster file) or a file that assigns any cell otherwise is not.
fn asks_for(args: &AnnotateArgs, loaded: &Loaded, cells: &[Box<str>], ids: &[Option<u32>]) -> bool {
    let Some(rel) = loaded.manifest.cluster.clusters.as_deref() else {
        return false;
    };
    if args.clusters.is_some() {
        return false;
    }
    let Ok((names, own)) = read_clusters(&resolve(&loaded.dir, rel)) else {
        return false;
    };
    let own: std::collections::HashMap<&str, Option<u32>> =
        names.iter().map(AsRef::as_ref).zip(own).collect();
    cells
        .iter()
        .zip(ids)
        .all(|(c, id)| own.get(c.as_ref()).is_some_and(|o| o == id))
}

/// How a run's decoder turns a cell's latent row into a gene distribution.
#[derive(Clone, Copy)]
pub enum Decoder {
    /// `π = β · θ`: `θ = softmax(latent)` over the topics (`log θ` or raw
    /// `z`), `β` a gene distribution per topic.
    Mixture,
    /// `π = softmax_d(z · W + b)` over the decoder's features (scVI).
    Softmax,
}

impl Decoder {
    /// What it averages over a cluster's cells, in symbols.
    pub fn formula(self) -> &'static str {
        match self {
            Self::Mixture => "β · θ̄",
            Self::Softmax => "mean softmax(z·W + b)",
        }
    }
}

/// The decoder of `loaded`'s kind, with the latent and dictionary it decodes
/// (manifest-relative), when the run wrote both. A topic kind's β is its
/// empirical dictionary when it wrote one: per gene, from the counts, where
/// the model's own β is per gene module, split evenly over its genes when
/// the run coarsened them. A `vae`'s is its decoder's loadings.
fn decoder_of(loaded: &Loaded) -> Option<(Decoder, &str, &str)> {
    let o = &loaded.manifest.outputs;
    let latent = o.latent.as_deref()?;
    let model = o.softmax_dictionary.as_deref().or(o.dictionary.as_deref());
    let (d, dict) = match loaded.manifest.kind {
        RunKind::Topic | RunKind::Itopic | RunKind::JointTopic | RunKind::MaskedVae => (
            Decoder::Mixture,
            o.dictionary_empirical.as_deref().or(model)?,
        ),
        RunKind::Vae => (Decoder::Softmax, model?),
        _ => return None,
    };
    Some((d, latent, dict))
}

/// Each column of `dict` as a distribution over the genes: weights (the
/// empirical dictionary) normalised, log weights (the model's `log β`)
/// exponentiated first.
fn gene_distributions(dict: &Mat) -> Mat {
    if dict.iter().all(|&v| v >= 0.0) {
        let mut beta = dict.clone();
        for mut col in beta.column_iter_mut() {
            let total = col.sum();
            if total > 0.0 {
                col /= total;
            }
        }
        beta
    } else {
        dict.normalize_exp_logits_columns()
    }
}

/// The run's expected expression per cluster and batch, from its decoder.
fn from_decoder(args: &AnnotateArgs, loaded: &Loaded) -> Result<Expression> {
    let m = &loaded.manifest;
    let (decoder, latent_rel, dict_rel) =
        decoder_of(loaded).context("the run has no decoder to read expression from")?;
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
    let (gene_sum_kg, sum_pg) = match decoder {
        Decoder::Mixture => {
            let beta = gene_distributions(&dict.mat);
            let (tk, tp) = theta_sums(&latent.mat, &groups);
            (expected(&beta, &tk), expected(&beta, &tp))
        }
        Decoder::Softmax => {
            let (features, gene_of) = vae_features(loaded, &dict.rows, &dict.mat)?;
            let bias = vae_bias(loaded, features.nrows())?;
            let (fk, fp) = decoded_sums(&latent.mat, &features, &bias, &groups);
            let d = features.nrows();
            (to_genes(&fk, d, &gene_of), to_genes(&fp, d, &gene_of))
        }
    };
    let gene_weights = fisher_weights(loaded, &dict.rows);
    let pb_gene_gp = weighted_mean_profile(&sum_pg, n_batches, g, &gene_weights);
    info!(
        "cluster expression: the {} decoder's expected expression at each cluster's cells (no counts)",
        m.kind.as_str()
    );
    Ok(Expression {
        gene_names: dict.rows,
        cell_names,
        cluster_labels,
        n_clusters: k,
        batch_labels,
        n_batches,
        gene_sum_kg,
        pb_gene_gp,
        gene_weights,
        source: json!({
            "from": "decoder",
            "kind": m.kind.as_str(),
            "dictionary": dict_rel,
            "depth": NOMINAL_DEPTH,
        }),
    })
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
    let theta = latent.transpose().normalize_exp_logits_columns();
    let n_topics = theta.nrows();
    let mut tk = Mat::zeros(n_topics, groups.k);
    let mut tp = Mat::zeros(n_topics, groups.n_batches);
    for (n, col) in theta.column_iter().enumerate() {
        if groups.clusters[n] < groups.k {
            let mut c = tk.column_mut(groups.clusters[n]);
            c += &col;
        }
        let mut p = tp.column_mut(groups.batches[n]);
        p += &col;
    }
    (tk, tp)
}

/// `β · θsums` (`g × K` by `K × groups`) as row-major `groups · g` (a
/// column-major matrix's own order), at [`NOMINAL_DEPTH`] counts per cell.
fn expected(beta: &Mat, theta_sums: &Mat) -> Vec<f64> {
    (beta * theta_sums)
        .iter()
        .map(|&v| NOMINAL_DEPTH * f64::from(v))
        .collect()
}

/// Up to [`DECODED_CELLS`] of each group's cells, evenly spaced, and how many
/// cells each group has.
fn spaced(labels: &[usize], n_groups: usize) -> Vec<(Vec<usize>, usize)> {
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
            ((0..take).map(|i| m[i * m.len() / take]).collect(), m.len())
        })
        .collect()
}

/// `softmax_d(z·W + b)` over each cluster's and batch's cells, as sums over
/// all their cells (`groups · d`), estimated from [`spaced`] cells: every
/// group's chunks are decoded in parallel, a chunk's cells as the columns of
/// one matrix.
fn decoded_sums(
    latent: &Mat,
    features: &Mat,
    bias: &[f32],
    groups: &Groups,
) -> (Vec<f64>, Vec<f64>) {
    let d = features.nrows();
    let bias = Mat::from_column_slice(d, 1, bias);
    let picked: Vec<(bool, usize, Vec<usize>, usize)> = spaced(groups.clusters, groups.k)
        .into_iter()
        .enumerate()
        .map(|(c, (cells, n))| (true, c, cells, n))
        .chain(
            spaced(groups.batches, groups.n_batches)
                .into_iter()
                .enumerate()
                .map(|(c, (cells, n))| (false, c, cells, n)),
        )
        .collect();
    let tasks: Vec<(usize, &[usize])> = picked
        .iter()
        .enumerate()
        .flat_map(|(t, (_, _, cells, _))| cells.chunks(256).map(move |ch| (t, ch)))
        .collect();
    let partial: Vec<(usize, Vec<f64>)> = tasks
        .par_iter()
        .map(|&(t, chunk)| {
            let z = Mat::from_fn(latent.ncols(), chunk.len(), |k, r| latent[(chunk[r], k)]);
            let mut logits = features * z;
            for mut col in logits.column_iter_mut() {
                col += &bias.column(0);
            }
            logits.normalize_exp_logits_columns_inplace();
            let mut acc = vec![0f64; d];
            for col in logits.column_iter() {
                acc.iter_mut()
                    .zip(col.iter())
                    .for_each(|(a, &p)| *a += f64::from(p));
            }
            (t, acc)
        })
        .collect();
    let mut out_k = vec![0f64; groups.k * d];
    let mut out_p = vec![0f64; groups.n_batches * d];
    for (t, acc) in partial {
        let (is_cluster, c, cells, n) = &picked[t];
        let out = if *is_cluster { &mut out_k } else { &mut out_p };
        let scale = NOMINAL_DEPTH * *n as f64 / cells.len() as f64;
        for (o, v) in out[c * d..][..d].iter_mut().zip(acc) {
            *o += v * scale;
        }
    }
    (out_k, out_p)
}

/// Feature sums (`groups · d`) to gene sums (`groups · g`) by each gene's
/// feature and share of it.
fn to_genes(sum: &[f64], d: usize, gene_of: &[(usize, f64)]) -> Vec<f64> {
    let g = gene_of.len();
    let groups = sum.len().checked_div(d).unwrap_or(0);
    let mut out = vec![0f64; groups * g];
    for c in 0..groups {
        for (i, &(f, share)) in gene_of.iter().enumerate() {
            out[c * g + i] = sum[c * d + f] * share;
        }
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
    let prefix = loaded.model_prefix();
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
        .unwrap_or_else(|| format!("{}.safetensors", loaded.model_prefix()));
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
    use legume_numeric::matrix::common_io::read_lines;
    let files = loaded.manifest.data_batches(&loaded.dir);
    let from_files: Option<Vec<Box<str>>> = if files.is_empty() {
        None
    } else {
        files
            .iter()
            .map(|f| read_lines(f).ok())
            .collect::<Option<Vec<_>>>()
            .map(|v| v.concat())
    };
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
    match data_beans::alg::gene_weighting::load_fisher_weights(&loaded.model_prefix()) {
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
