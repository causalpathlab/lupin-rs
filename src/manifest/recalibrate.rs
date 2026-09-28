//! Rescoring a relabel round. An enrichment pass caches its sufficient
//! statistics ([`write_cache`]); a later round, whose clusters are merges of
//! the pass's clusters and whose marker panel may be edited, reruns the same
//! scoring on them ([`rescore`]) without re-reading the counts: merged
//! clusters' gene sums add, the batch × cluster membership is rebuilt from
//! the cells, and the panel is rebuilt from the round's markers.

use crate::annotate::aggregate::weighted_mean_profile;
use crate::annotate::args::AnnotateArgs;
use crate::annotate::by_enrichment;
use crate::annotate::inputs::EnrichmentInputs;
use crate::annotate::rounds::{parse_cluster_id, ClusterId};
use crate::manifest::rounds::{read_clusters, write_clusters, Cells};
use crate::manifest::run::{resolve, Loaded, StatsCache};
use anyhow::{Context, Result};
use legume_numeric::matrix::dense_mat_io::{read_mat, Mat};
use legume_numeric::matrix::traits::IoOps;
use log::info;
use std::collections::{BTreeMap, HashMap};

const GENE_SUM: &str = ".cluster_gene_sum.parquet";
const BATCH_PROFILE: &str = ".batch_profile.parquet";
const GENE_WEIGHT: &str = ".gene_weight.parquet";
const CELL_BATCH: &str = ".cell_batch.parquet";

/// Write an enrichment pass's sufficient statistics beside its outputs under
/// `out`; the paths are as written, for [`crate::manifest::annotate`] to
/// record. `None` for a pass with no per-batch profile (the GO/GMT pass),
/// which later rounds cannot rescore.
pub fn write_cache(out: &str, inputs: &EnrichmentInputs) -> Result<Option<StatsCache>> {
    let Some(pb) = &inputs.pb_gene_gp else {
        return Ok(None);
    };
    let g = inputs.gene_names.len();
    let k = inputs.n_clusters;
    let mut sums = Mat::zeros(g, k);
    for c in 0..k {
        for (i, v) in inputs.gene_sum_kg[c * g..(c + 1) * g].iter().enumerate() {
            sums[(i, c)] = *v as f32;
        }
    }
    let cols: Vec<Box<str>> = (0..k).map(|c| format!("K{c}").into()).collect();
    let rows = Some(&inputs.gene_names[..]);
    let gene_sum = format!("{out}{GENE_SUM}");
    sums.to_parquet_with_names(&gene_sum, (rows, Some("gene")), Some(&cols))?;
    let batch_cols: Vec<Box<str>> = (0..pb.ncols()).map(|b| format!("B{b}").into()).collect();
    let batch_profile = format!("{out}{BATCH_PROFILE}");
    pb.to_parquet_with_names(&batch_profile, (rows, Some("gene")), Some(&batch_cols))?;
    let mut w = Mat::zeros(g, 1);
    for (i, v) in inputs.gene_weights.iter().enumerate() {
        w[(i, 0)] = *v;
    }
    let gene_weight = format!("{out}{GENE_WEIGHT}");
    w.to_parquet_with_names(&gene_weight, (rows, Some("gene")), Some(&["weight".into()]))?;
    let batches: Vec<Option<ClusterId>> = inputs
        .batch_labels
        .iter()
        .map(|&b| u32::try_from(b).ok())
        .collect();
    let cell_batch = format!("{out}{CELL_BATCH}");
    write_clusters(&cell_batch, &inputs.cell_names, &batches)?;
    info!("cached the pass's statistics for rescoring later rounds");
    Ok(Some(StatsCache {
        gene_sum,
        batch_profile,
        gene_weight,
        cell_batch,
    }))
}

/// A round rescored: per cluster (rows, by `ids`) × cell type (`types`).
pub struct Rescored {
    pub ids: Vec<ClusterId>,
    pub types: Vec<Box<str>>,
    /// Softmaxed Q: probabilities over types.
    pub q_probs: Mat,
    /// FDR q-values.
    pub q_values: Mat,
    /// Bootstrap support per type plus a trailing `unassigned` column, when
    /// the bootstrap ran.
    pub support: Option<Mat>,
}

impl Rescored {
    /// Row names `K{id}`, as every cluster table uses.
    #[must_use]
    pub fn row_names(&self) -> Vec<Box<str>> {
        self.ids.iter().map(|id| format!("K{id}").into()).collect()
    }
}

/// Rescore the clusters of `after` against `panel` from `source`'s cached
/// statistics, with the settings `source`'s enrichment pass recorded and
/// `n_boot` bootstrap resamples (0 for none). `None` when `source` carries no
/// cache (a projection run, or a round from before caching).
pub(super) fn rescore(
    source: &Loaded,
    after: &Cells,
    panel: &[(String, String)],
    n_boot: usize,
) -> Result<Option<Rescored>> {
    let a = &source.manifest.annotate;
    let (Some(cache), Some(ids_rel)) = (&a.stats_cache, &a.expression_clusters) else {
        return Ok(None);
    };
    if panel.is_empty() {
        return Ok(None);
    }
    let mut args: AnnotateArgs = a
        .settings
        .as_ref()
        .and_then(|s| s.get("enrichment"))
        .cloned()
        .map(serde_json::from_value)
        .transpose()
        .context("reading the enrichment pass's settings")?
        .context("the round records no enrichment settings")?;
    args.n_boot = n_boot;
    let at = |rel: &str| resolve(&source.dir, rel);

    // The pass's cluster of each cell, and which new cluster each of those
    // became (merges only, so each goes to one).
    let (names, orig) = read_clusters(&at(ids_rel))?;
    let orig_of: HashMap<&str, ClusterId> = names
        .iter()
        .zip(&orig)
        .filter_map(|(n, id)| id.map(|id| (n.as_ref(), id)))
        .collect();
    let new_ids: Vec<ClusterId> = {
        let mut v: Vec<ClusterId> = after.clusters.iter().flatten().copied().collect();
        v.sort_unstable();
        v.dedup();
        v
    };
    let slot: HashMap<ClusterId, usize> =
        new_ids.iter().enumerate().map(|(i, id)| (*id, i)).collect();
    let mut votes: BTreeMap<ClusterId, BTreeMap<ClusterId, usize>> = BTreeMap::new();
    for (cell, id) in after.names.iter().zip(&after.clusters) {
        if let (Some(id), Some(&o)) = (id, orig_of.get(cell.as_ref())) {
            *votes.entry(o).or_default().entry(*id).or_default() += 1;
        }
    }
    let became: HashMap<ClusterId, usize> = votes
        .into_iter()
        .filter_map(|(o, to)| {
            let (id, _) = to.into_iter().max_by_key(|(_, n)| *n)?;
            Some((o, slot[&id]))
        })
        .collect();

    // Gene sums of the new clusters, then their profile.
    let sums = read_mat(&at(&cache.gene_sum))?;
    let g = sums.rows.len();
    let k = new_ids.len();
    let mut gene_sum_kg = vec![0f64; g * k];
    for (j, col) in sums.cols.iter().enumerate() {
        let Some(&dest) = parse_cluster_id(col).and_then(|o| became.get(&o)) else {
            continue;
        };
        for i in 0..g {
            gene_sum_kg[dest * g + i] += f64::from(sums.mat[(i, j)]);
        }
    }
    let weights = read_mat(&at(&cache.gene_weight))?;
    let weights: Vec<f32> = (0..g).map(|i| weights.mat[(i, 0)]).collect();
    let profile_gk = weighted_mean_profile(&gene_sum_kg, k, g, &weights);

    // Each cell's new cluster and its batch, in `after`'s cell order.
    let (batch_cells, batch_ids) = read_clusters(&at(&cache.cell_batch))?;
    let batch_of: HashMap<&str, usize> = batch_cells
        .iter()
        .zip(&batch_ids)
        .filter_map(|(n, b)| b.map(|b| (n.as_ref(), b as usize)))
        .collect();
    let pb = read_mat(&at(&cache.batch_profile))?;
    let n_batches = pb.mat.ncols();
    let cluster_labels: Vec<usize> = after
        .clusters
        .iter()
        .map(|id| {
            id.and_then(|id| slot.get(&id).copied())
                .unwrap_or(usize::MAX)
        })
        .collect();
    let batch_labels: Vec<usize> = after
        .names
        .iter()
        .map(|n| batch_of.get(n.as_ref()).copied().unwrap_or(0))
        .collect();

    let pairs: Vec<(Box<str>, Box<str>)> = panel
        .iter()
        .map(|(g, t)| (g.as_str().into(), t.as_str().into()))
        .collect();
    let annot = crate::annotate::markers::annotation_matrix_from_pairs(&pairs, &sums.rows)?;

    let inputs = EnrichmentInputs {
        gene_names: sums.rows.clone(),
        cell_names: after.names.clone(),
        cluster_labels,
        n_clusters: k,
        batch_labels,
        n_batches,
        markers_gc: annot.membership_ga,
        celltype_names: annot.annot_names.clone(),
        profile_gk,
        pb_gene_gp: Some(pb.mat),
        gene_sum_kg: Vec::new(),
        gene_weights: Vec::new(),
    };
    info!(
        "rescoring {k} cluster(s) against {} cell type(s){}",
        annot.annot_names.len(),
        if n_boot > 0 {
            format!(", bootstrap {n_boot}")
        } else {
            String::new()
        }
    );
    let out = by_enrichment::score(&args, &inputs)?;
    let support = out.bootstrap.map(|b| {
        let w = b.c + 1;
        let mut m = Mat::zeros(k, w);
        for r in 0..k {
            for c in 0..w {
                m[(r, c)] = b.consensus.post[r * w + c];
            }
        }
        m
    });
    Ok(Some(Rescored {
        ids: new_ids,
        types: annot.annot_names,
        q_probs: out.q_kc,
        q_values: out.qvalue_kc,
        support,
    }))
}
