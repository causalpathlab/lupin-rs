//! Rescoring a relabel round. An enrichment pass caches its sufficient
//! statistics ([`write_cache`]); a later round, whose clusters are merges of
//! the pass's clusters and whose marker panel may be edited, reruns the same
//! scoring on them ([`rescore`]) without re-reading the counts: merged
//! clusters' gene sums add, the batch × cluster membership is rebuilt from
//! the cells, and the panel is rebuilt from the round's markers.

use super::data_files::{ClData, Fetch, SearchPath};
use crate::annotate::aggregate::weighted_mean_profile;
use crate::annotate::args::AnnotateArgs;
use crate::annotate::by_enrichment;
use crate::annotate::inputs::EnrichmentInputs;
use crate::annotate::markers::label_key;
use crate::annotate::rounds::{parse_cluster_id, ClusterId};
use crate::manifest::rounds::{read_clusters, write_clusters, Cells};
use crate::manifest::run::{resolve, Loaded, StatsCache};
use anyhow::{Context, Result};
use legume_numeric::matrix::dense_mat_io::{read_mat, Mat};
use legume_numeric::matrix::traits::IoOps;
use log::info;
use std::collections::{BTreeMap, BTreeSet, HashMap};

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
    pub p_values: Mat,
    /// The ES restandardized by the gene-set null.
    pub z: Mat,
    /// fgsea's normalized enrichment score.
    pub nes: Mat,
    /// z = Φ⁻¹(1 − p), which the Q probabilities are a softmax of.
    pub probit_z: Mat,
}

impl Rescored {
    /// Row names `K{id}`, as every cluster table uses.
    #[must_use]
    pub fn row_names(&self) -> Vec<Box<str>> {
        self.ids.iter().map(|id| format!("K{id}").into()).collect()
    }
}

/// The Cell Ontology data `source`'s pass used, read from the files it
/// recorded; failing that, found beside the run without downloading.
fn pass_cl_data(source: &Loaded, args: &AnnotateArgs) -> Result<ClData> {
    let search = || SearchPath::new(Some(&source.dir));
    let recorded = source
        .manifest
        .annotate
        .settings
        .as_ref()
        .and_then(|s| s.pointer("/enrichment/cell_ontology"));
    if let Some(r) = recorded {
        if let Some(d) = ClData::from_record(r, search())? {
            return Ok(d);
        }
    }
    ClData::load(
        search(),
        args.obo.as_deref(),
        args.label_cl.as_deref(),
        Fetch::Never,
    )
}

/// Rescore the clusters of `after` against `panel` from `source`'s cached
/// statistics, with the settings `source`'s enrichment pass recorded. `None`
/// when `source` carries no
/// cache (a projection run, or a round from before caching).
pub(super) fn rescore(
    source: &Loaded,
    after: &Cells,
    panel: &[(String, String)],
) -> Result<Option<Rescored>> {
    let Some((args, inputs, ids)) = rescore_inputs(source, after, panel)? else {
        return Ok(None);
    };
    score_all(&args, inputs, ids).map(Some)
}

/// Every cluster × type of `inputs` scored, rows by `ids`.
fn score_all(
    args: &AnnotateArgs,
    inputs: EnrichmentInputs,
    ids: Vec<ClusterId>,
) -> Result<Rescored> {
    info!(
        "rescoring {} cluster(s) against {} cell type(s)",
        ids.len(),
        inputs.celltype_names.len()
    );
    let out = by_enrichment::score(args, &inputs)?;
    Ok(Rescored {
        ids,
        types: inputs.celltype_names,
        q_probs: out.q_kc,
        q_values: out.qvalue_kc,
        p_values: out.pvalue_kc,
        z: out.es_restandardized_kc,
        nes: out.nes_kc,
        probit_z: out.z_kc,
    })
}

/// [`rescore`] of a round whose clusters are `source`'s own, where only the
/// cell types in `touched` (label keys) scored differently: those are scored
/// again, every other type keeps `source`'s p-value and NES, and the q-values
/// and Q are redone over the whole row. The types a marker edit touches are
/// the ones it edits and the ones sharing an edited gene. A panel that gains
/// or loses a type moves every type's IDF weight, so it is rescored in full,
/// as it is when `source` records no p-values or NES for a cluster. `z` is
/// left NaN for the types not rescored: this is for previews, not for
/// writing.
pub(super) fn rescore_types(
    source: &Loaded,
    after: &Cells,
    panel: &[(String, String)],
    touched: &BTreeSet<String>,
) -> Result<Option<Rescored>> {
    let a = &source.manifest.annotate;
    let table = |rel: &Option<String>| {
        rel.as_deref()
            .map(|r| crate::manifest::rounds::read_table(&resolve(&source.dir, r)))
            .transpose()
    };
    let (Some(p_old), Some(nes_old)) = (
        table(&a.cluster_celltype_p)?,
        table(&a.cluster_celltype_nes)?,
    ) else {
        return rescore(source, after, panel);
    };
    let Some((args, inputs, ids)) = rescore_inputs(source, after, panel)? else {
        return Ok(None);
    };
    let names = &inputs.celltype_names;
    let keys = |cols: &[String]| cols.iter().map(|c| label_key(c)).collect::<BTreeSet<_>>();
    let same_types = keys(&p_old.cols) == names.iter().map(|n| label_key(n)).collect()
        && keys(&nes_old.cols) == keys(&p_old.cols);
    let rows_known = ids
        .iter()
        .all(|id| p_old.row(*id).is_some() && nes_old.row(*id).is_some());
    if !(same_types && rows_known) {
        return score_all(&args, inputs, ids).map(Some);
    }
    let col = |t: &crate::annotate::rounds::Table, name: &str| {
        t.cols.iter().position(|c| label_key(c) == label_key(name))
    };
    let (k, c) = (ids.len(), names.len());
    let mut p = Mat::zeros(k, c);
    let mut nes = Mat::zeros(k, c);
    let mut z = Mat::from_element(k, c, f32::NAN);
    for (i, id) in ids.iter().enumerate() {
        for (t, name) in names.iter().enumerate() {
            if let (Some(jp), Some(jn)) = (col(&p_old, name), col(&nes_old, name)) {
                p[(i, t)] = p_old.row(*id).map_or(1.0, |r| r[jp]);
                nes[(i, t)] = nes_old.row(*id).map_or(0.0, |r| r[jn]);
            }
        }
    }
    let types: Vec<usize> = (0..c)
        .filter(|&t| touched.contains(&label_key(&names[t])))
        .collect();
    let (markers_gc, config) = if types.is_empty() {
        let (_, m, cfg) = by_enrichment::prepare(&args, &inputs)?;
        (m, cfg)
    } else {
        let (scores, m, cfg) = by_enrichment::score_types(&args, &inputs, &types)?;
        for (j, &t) in scores.types.iter().enumerate() {
            for i in 0..k {
                p[(i, t)] = scores.pvalue_kc[(i, j)];
                nes[(i, t)] = scores.nes_kc[(i, j)];
                z[(i, t)] = scores.es_restandardized_kc[(i, j)];
            }
        }
        (m, cfg)
    };
    let adjusted = enrichment::adjust(&p, &markers_gc, &config)?;
    Ok(Some(Rescored {
        ids,
        types: inputs.celltype_names,
        q_probs: adjusted.q_kc,
        q_values: adjusted.qvalue_kc,
        p_values: adjusted.pvalue_kc,
        z,
        nes,
        probit_z: adjusted.z_kc,
    }))
}

/// What [`rescore`] scores: the enrichment settings of `source`'s pass, its
/// cached statistics regrouped by `after`'s clusters against `panel`, and the
/// clusters' ids by row. `None` when the pass cached nothing or the panel is
/// empty.
fn rescore_inputs(
    source: &Loaded,
    after: &Cells,
    panel: &[(String, String)],
) -> Result<Option<(AnnotateArgs, EnrichmentInputs, Vec<ClusterId>)>> {
    let a = &source.manifest.annotate;
    let (Some(cache), Some(ids_rel)) = (&a.stats_cache, &a.expression_clusters) else {
        return Ok(None);
    };
    if panel.is_empty() {
        return Ok(None);
    }
    let args: AnnotateArgs = a
        .settings
        .as_ref()
        .and_then(|s| s.get("enrichment"))
        .cloned()
        .map(serde_json::from_value)
        .transpose()
        .context("reading the enrichment pass's settings")?
        .context("the round records no enrichment settings")?;
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
        type_tree: Some(
            super::ontology::panel_tree(&pass_cl_data(source, &args)?, panel)?
                .treebh(&annot.annot_names),
        ),
        cl_record: None,
    };
    Ok(Some((args, inputs, new_ids)))
}
