//! Run-manifest glue for the annotation commands: resolve `-f run.senna.json`
//! into the loaded matrices [`crate::annotate`] takes, then write the artifact
//! paths it returns back into `manifest.annotate.*`. This is the only place
//! that touches [`RunManifest`] on behalf of an annotation command.
//!
//! Re-opening the raw counts and aggregating them per cluster also lives here,
//! over the sparse count backend; [`crate::annotate`] receives the resulting
//! cluster × gene expression.

use crate::annotate::args::{
    fixed_settings, AnnotateArgs, AnnotateOntologyArgs, MarkerScoreArgs, ProjectionMethod,
    UnitsArgs, BLOCK_SIZE,
};
use crate::annotate::by_enrichment;
use crate::annotate::inputs::{load_cluster_labels, EnrichmentInputs};
use crate::annotate::marker_score;
use crate::annotate::ontology;
use crate::annotate::outputs::AnnotationOutputs;
use crate::annotate::rounds::Method;
use crate::manifest::rounds;
use crate::manifest::run::{annotated_path, resolve, Loaded};

use crate::annotate::aggregate::{accumulate_gene_sum_pair, weighted_mean_profile};
use crate::annotate::markers::build_annotation_matrix;
use crate::manifest::run::{CellSpace, RunManifest};
use crate::marker_embedding::{load_feature_embedding, load_marker_feature_embedding};
use data_beans::aux::data_loading::{
    read_data_on_shared_rows, ReadSharedRowsArgs, SparseDataWithBatch,
};
use legume_numeric::matrix::dense_mat_io::Mat;

use anyhow::{anyhow, Context, Result};
use legume_numeric::matrix::traits::{IoOps, MatWithNames};
use log::info;
use rayon::prelude::*;
use rustc_hash::FxHashMap as HashMap;
use std::path::Path;

/// CLI passthrough for the internal Leiden fallback (used only when neither
/// `--clusters` nor `manifest.cluster.clusters` is provided).
#[derive(Debug, Clone, serde::Serialize)]
pub struct LeidenArgs {
    pub knn: usize,
    pub resolution: f64,
    pub num_clusters: Option<usize>,
    pub min_cluster_size: usize,
    pub seed: Option<u64>,
}

/// The `lupin annotate` enrichment defaults (`--knn 15 --resolution 1`), seeded.
impl Default for LeidenArgs {
    fn default() -> Self {
        Self {
            knn: 15,
            resolution: 1.0,
            num_clusters: None,
            min_cluster_size: 2,
            seed: Some(42),
        }
    }
}

/// `lupin annotate --method enrichment`: re-open the raw counts the manifest
/// points at, aggregate them per cluster, run the enrichment, and record what
/// it wrote.
/// `data` is the run's Cell Ontology data when the caller already loaded it.
pub fn annotate_by_enrichment(
    args: &AnnotateArgs,
    loaded: &Loaded,
    data: Option<&super::data_files::ClData>,
) -> Result<()> {
    let plan = by_enrichment::plan(args)?;
    let stages = by_enrichment::pass_stages(plan.gene_sets_too);
    stages.start(0);
    let inputs = load_enrichment_inputs(args, loaded, data)?;
    let mut outputs = by_enrichment::run(args, &plan, &inputs)?;
    // SuSiE makes the marker call over the enrichment's statistics, on raw
    // count sums (a decoder's expectation is none). NON-FATAL, like the
    // ontology walk: the enrichment's outputs are already written.
    let from_decoder = inputs
        .expression_source
        .as_ref()
        .is_some_and(|e| e["from"] == "decoder");
    // Unless another method's reading is asked to make the call (`--call`).
    let susie_calls = matches!(args.call, None | Some(Method::Susie));
    match (&args.susie, plan.ontology_mode, from_decoder) {
        (Some(cfg), false, false) => {
            let fit =
                crate::annotate::by_susie::run(&args.out, cfg, &inputs, &mut outputs, susie_calls);
            match fit {
                Err(e) if args.call == Some(Method::Susie) => {
                    return Err(
                        e.context("--call susie: the SuSiE stage, which makes the call, failed")
                    );
                }
                Err(e) => log::error!(
                    "SuSiE stage failed ({e:#}); the call is the enrichment's softmax share"
                ),
                Ok(()) => {}
            }
        }
        (Some(_), false, true) => {
            info!("no raw counts (decoder expression): the call is the enrichment's softmax share");
        }
        _ => {}
    }
    // The other methods' readings beside the enrichment's, as auto asks,
    // where the run has the gene table each reads: the score where it
    // decodes cells, the units where its genes are placed among the cells.
    // NON-FATAL, as the enrichment's call stands without them, unless the
    // method failing is the one asked to make the call.
    let beside = !plan.ontology_mode && !args.markers.is_empty();
    let (cells, groups) = (&inputs.cell_names, &inputs.cluster_labels);
    if args.score_beside && beside {
        match score_beside(loaded, &args.markers, &args.out, cells, groups) {
            Ok(table) => outputs.marker_score = table,
            Err(e) if args.call == Some(Method::Score) => {
                return Err(
                    e.context("--call score: the projection score, which makes the call, failed")
                );
            }
            Err(e) => log::warn!("no projection score beside the enrichment: {e:#}"),
        }
    }
    if args.units_beside && beside {
        match units_beside(loaded, &args.units_beside(), cells, groups) {
            Ok(Some(u)) => {
                outputs.units_call = u.units_call;
                outputs.celltype_units = u.celltype_units;
                outputs.celltype_links = u.celltype_links;
                outputs.celltype_membership = u.celltype_membership;
                outputs.term_alignment = u.term_alignment;
                outputs.term_units = u.term_units;
                outputs.unit_anchors = u.unit_anchors;
            }
            Ok(None) => {}
            Err(e) if args.call == Some(Method::Units) => {
                return Err(
                    e.context("--call units: the gene-set units, which make the call, failed")
                );
            }
            Err(e) => log::warn!("no gene-set units beside the enrichment: {e:#}"),
        }
    }
    take_call(args.call, &mut outputs, &args.out, cells)?;
    // The ids the cluster tables' `K{id}` rows refer to, so later rounds and
    // viewers key on the same clusters.
    let clusters_path = format!("{}{}", args.out, rounds::CLUSTERS);
    let ids = rounds::cluster_ids(&inputs.cluster_labels);
    rounds::write_clusters(&clusters_path, &inputs.cell_names, &ids)?;
    outputs.clusters = Some(clusters_path);
    // The decoder's expected expression is no count sum: a later pass or
    // rescoring must not take it for one (`from_decoder`, above).
    if !plan.ontology_mode && !from_decoder {
        outputs.stats_cache = Some(crate::manifest::recalibrate::write_cache(
            &args.out, &inputs,
        )?);
    }
    let pass = if plan.ontology_mode {
        Pass::GeneSets
    } else {
        let called = args.call.or(if outputs.cluster_celltype_pip.is_some() {
            Some(Method::Susie)
        } else {
            outputs
                .cluster_celltype_q
                .as_ref()
                .map(|_| Method::Enrichment)
        });
        Pass::Markers(&args.markers, called)
    };
    // Which nulls actually ran: the sample permutation is skipped when there
    // are too few batches to shuffle.
    let mut used = settings(args)?;
    if let serde_json::Value::Object(m) = &mut used {
        // The SuSiE settings that made the call, or none when it was skipped
        // or failed and the softmax share made it.
        if outputs.cluster_celltype_pip.is_none() {
            m.insert("susie".into(), serde_json::Value::Null);
        }
        let draws = by_enrichment::sample_perm_draws(inputs.n_batches, args.num_perm);
        if let Some(r) = &inputs.cl_record {
            m.insert("cell_ontology".into(), r.clone());
        }
        // Not the raw counts: an earlier pass's cache or the model's decoder.
        if let Some(e) = &inputs.expression_source {
            m.insert("expression".into(), e.clone());
        }
        m.insert(
            "null".into(),
            serde_json::json!({
                "gene_set_randomization": crate::annotate::args::NUM_DRAWS,
                "sample_permutation": draws,
                "batches": inputs.n_batches,
            }),
        );
    }
    stages.start(stages.named(crate::annotate::by_enrichment::STAGE_RECORD));
    record(loaded, &args.out, pass, &outputs, "enrichment", used)?;
    stages.finish();
    Ok(())
}

/// `lupin annotate --method projection --projection score`: the embedding
/// marker score ([`marker_score`]) on the run's decoder gene table and cell
/// embedding, centred on its clusters (`--clusters`, else the run's, else
/// internal Leiden), which also become the round's cluster table.
pub fn annotate_by_marker_score(args: &MarkerScoreArgs, loaded: &Loaded) -> Result<()> {
    let (feat, cells) = score_tables(loaded)?;
    let (groups, _) = resolve_clusters(
        args.clusters.as_deref(),
        &args.leiden,
        &loaded.manifest,
        &loaded.dir,
        &cells.rows,
    )?;
    let mut outputs = marker_score::run(args, &feat, &cells, &groups)?;
    let clusters = format!("{}{}", args.out, rounds::CLUSTERS);
    rounds::write_clusters(&clusters, &cells.rows, &rounds::cluster_ids(&groups))?;
    outputs.clusters = Some(clusters);
    record(
        loaded,
        &args.out,
        Pass::Markers(&args.markers, Some(Method::Score)),
        &outputs,
        "projection",
        projection_settings(args, ProjectionMethod::Score)?,
    )
}

/// `lupin annotate --method projection --projection units`: the panel's
/// cell types trained as gene-set units on the genes placed among the cells
/// ([`crate::annotate::geneset_units`]), each cell called by its region, the
/// round tabled by `--clusters`, else the run's, else internal Leiden.
pub fn annotate_by_units(args: &UnitsArgs, loaded: &Loaded) -> Result<()> {
    let run = loaded.file.display().to_string();
    let feat = load_marker_feature_embedding(&loaded.manifest, &loaded.dir, &run)?;
    let cell_path = loaded.geometry_latent_path()?;
    let cells = Mat::from_parquet(&cell_path)
        .with_context(|| format!("reading cell embedding {cell_path}"))?;
    let (groups, _) = resolve_clusters(
        args.clusters.as_deref(),
        &args.leiden,
        &loaded.manifest,
        &loaded.dir,
        &cells.rows,
    )?;
    let mut outputs = crate::annotate::geneset_units::run(args, &feat, &cells, &groups)?;
    let clusters = format!("{}{}", args.out, rounds::CLUSTERS);
    rounds::write_clusters(&clusters, &cells.rows, &rounds::cluster_ids(&groups))?;
    outputs.clusters = Some(clusters);
    record(
        loaded,
        &args.out,
        Pass::Markers(&args.markers, Some(Method::Units)),
        &outputs,
        "projection",
        projection_settings(args, ProjectionMethod::Units)?,
    )
}

/// What the score reads: the run's decoder gene table and its cell embedding.
fn score_tables(loaded: &Loaded) -> Result<(MatWithNames<Mat>, MatWithNames<Mat>)> {
    let run = loaded.file.display().to_string();
    let feat = load_feature_embedding(&loaded.manifest, &loaded.dir, &run)?;
    let cell_path = loaded.geometry_latent_path()?;
    let cells = Mat::from_parquet(&cell_path)
        .with_context(|| format!("reading cell embedding {cell_path}"))?;
    Ok((feat, cells))
}

/// The projection score as a second reading beside another pass: when the
/// run decodes cells, score those of `cells` (in `groups`, that pass's
/// clusters) and write only the score's per-cell table under `out`; its
/// path, or `None` when the run has no decoder gene table.
pub fn score_beside(
    loaded: &Loaded,
    markers: &str,
    out: &str,
    cells: &[Box<str>],
    groups: &[usize],
) -> Result<Option<String>> {
    if crate::marker_embedding::projection_gene_table(&loaded.manifest, ProjectionMethod::Score)
        .is_err()
    {
        return Ok(None);
    }
    let (feat, embedding) = score_tables(loaded)?;
    let (picked, kept) = pass_rows(embedding, cells, groups)?;
    marker_score::alongside(markers, out, &feat, &picked, &kept).map(Some)
}

/// The gene-set units as a reading beside another pass: when the run has
/// the gene table they read, train them on `args`' panel (and terms) and fit
/// them to those of `cells` the cell embedding has, the terms aligned with
/// `groups` (that pass's clusters); every file of theirs but the per-cell
/// labels is written under `args.out`. `None` when the run has no such table.
pub fn units_beside(
    loaded: &Loaded,
    args: &UnitsArgs,
    cells: &[Box<str>],
    groups: &[usize],
) -> Result<Option<AnnotationOutputs>> {
    if crate::marker_embedding::projection_gene_table(&loaded.manifest, ProjectionMethod::Units)
        .is_err()
    {
        return Ok(None);
    }
    let run = loaded.file.display().to_string();
    let feat = load_marker_feature_embedding(&loaded.manifest, &loaded.dir, &run)?;
    let cell_path = loaded.geometry_latent_path()?;
    let embedding = Mat::from_parquet(&cell_path)
        .with_context(|| format!("reading cell embedding {cell_path}"))?;
    let (picked, kept) = pass_rows(embedding, cells, groups)?;
    let terms = args.terms(&feat.rows)?;
    crate::annotate::geneset_units::alongside(
        &args.markers,
        &args.out,
        &feat,
        &picked,
        &kept,
        &terms,
        &args.config()?,
    )
    .map(Some)
}

/// The rows of `embedding` that are the pass's `cells`, in its order, with
/// their clusters of `groups`; the embedding as it is when they are its
/// rows already.
fn pass_rows(
    embedding: MatWithNames<Mat>,
    cells: &[Box<str>],
    groups: &[usize],
) -> Result<(MatWithNames<Mat>, Vec<usize>)> {
    if embedding.rows.as_slice() == cells {
        return Ok((embedding, groups.to_vec()));
    }
    let row: HashMap<&str, usize> = embedding
        .rows
        .iter()
        .enumerate()
        .map(|(i, c)| (c.as_ref(), i))
        .collect();
    let (rows, kept): (Vec<usize>, Vec<usize>) = cells
        .iter()
        .zip(groups)
        .filter_map(|(c, &g)| Some((*row.get(c.as_ref())?, g)))
        .unzip();
    anyhow::ensure!(
        !rows.is_empty(),
        "the cell embedding shares no cell with the pass's cells"
    );
    let picked = MatWithNames {
        mat: embedding.mat.select_rows(&rows),
        rows: rows.iter().map(|&i| embedding.rows[i].clone()).collect(),
        cols: embedding.cols,
    };
    Ok((picked, kept))
}

/// Make the round's call `call`'s when it is a reading beside the
/// enrichment's (`--call score` or `units`): its per-cell calls become the
/// round's labels over the pass's `cells` (a cell it has no row for is
/// unassigned), with the units' membership as the probability (the score
/// has none). The enrichment's per-cell posterior goes from the record, as
/// it is not the call's. SuSiE's call was made by its stage; `--call susie`
/// only checks it was.
fn take_call(
    call: Option<Method>,
    outputs: &mut AnnotationOutputs,
    out: &str,
    cells: &[Box<str>],
) -> Result<()> {
    use legume_numeric::matrix::parquet::read_table_columns;
    let (table, column) = match call {
        Some(Method::Score) => (&outputs.marker_score, "margin"),
        Some(Method::Units) => (&outputs.units_call, "membership"),
        Some(Method::Susie) => {
            anyhow::ensure!(
                outputs.cluster_celltype_pip.is_some(),
                "--call susie: the SuSiE stage did not run on this pass (decoder expression, \
                 or --no-susie)"
            );
            return Ok(());
        }
        Some(Method::Enrichment) | None => return Ok(()),
    };
    let name = call.map_or("", Method::name);
    let path = table
        .as_deref()
        .with_context(|| format!("--call {name}: the {name} gave no reading on this run"))?;
    let (strings, nums) = read_table_columns(path, &["cell", "call"], &[column])
        .with_context(|| format!("reading {path}"))?;
    let by: HashMap<&str, (&str, f32)> = (strings[0].iter().zip(&strings[1]).zip(&nums[0]))
        .map(|((c, l), &v)| (c.as_ref(), (l.as_ref(), v as f32)))
        .collect();
    let (labels, probs): (Vec<Box<str>>, Vec<f32>) = cells
        .iter()
        .map(|c| match by.get(c.as_ref()) {
            Some(&(l, v)) => (
                l.into(),
                if call == Some(Method::Units) {
                    v
                } else {
                    f32::NAN
                },
            ),
            None => (enrichment::UNASSIGNED_LABEL.into(), f32::NAN),
        })
        .unzip();
    graph_embedding_util::type_annotation::write_cell_labels(out, cells, &labels, &probs)?;
    info!("the round's call is the {name}'s");
    outputs.argmax = Some(crate::annotate::outputs::labels_path(out));
    outputs.annotation = None;
    Ok(())
}

/// A projection pass's settings, naming the way it called.
fn projection_settings(
    args: &impl serde::Serialize,
    projection: ProjectionMethod,
) -> Result<serde_json::Value> {
    let mut v = settings(args)?;
    if let serde_json::Value::Object(m) = &mut v {
        m.insert("projection".into(), projection.as_str().into());
    }
    Ok(v)
}

/// `lupin annotate` without markers: walk the CL tree over the cluster ×
/// cell-type matrix an earlier enrichment run recorded.
pub fn annotate_ontology(args: &AnnotateOntologyArgs, loaded: &Loaded) -> Result<()> {
    let q_rel = loaded
        .manifest
        .annotate
        .cluster_celltype_q
        .as_deref()
        .ok_or_else(|| {
            anyhow!(
            "{} has no `annotate.cluster_celltype_q` — run `lupin annotate --method enrichment \
             -m <markers>` on it first",
            loaded.file.display()
        )
        })?;
    let q_abs = resolve(&loaded.dir, q_rel);
    let outputs = ontology::run(args, &q_abs)?;
    record(
        loaded,
        &args.out,
        Pass::Ontology,
        &outputs,
        "ontology",
        settings(args)?,
    )
}

/// A pass's effective settings: its own arguments plus the fixed constants.
fn settings(args: &impl serde::Serialize) -> Result<serde_json::Value> {
    let mut v = serde_json::to_value(args)?;
    if let serde_json::Value::Object(m) = &mut v {
        m.insert("fixed".into(), fixed_settings());
    }
    Ok(v)
}

////////////////////////////////////
// annotation outputs → manifest  //
////////////////////////////////////

/// Which annotation pass ran, which decides the `manifest.annotate.*` slots it
/// owns. A pass writes every slot it owns — `None` clears a stale pointer
/// left by an earlier run — and leaves the rest alone.
enum Pass<'a> {
    /// Per-cell cell-type labels from a marker panel (enrichment or
    /// projection), and the method whose reading made the call.
    Markers(&'a str, Option<Method>),
    /// GO/GMT gene-set signature per cluster.
    GeneSets,
    /// CL ontology walk over an earlier enrichment.
    Ontology,
}

/// Record a pass's artifacts in the manifest, relative to its directory, plus
/// the settings it ran with under `annotate.settings.{method}`, and save it
/// as a new manifest (see [`annotated_path`]); the one it was read from is left as is.
fn record(
    loaded: &Loaded,
    out_prefix: &str,
    pass: Pass<'_>,
    out: &AnnotationOutputs,
    method: &str,
    settings: serde_json::Value,
) -> Result<()> {
    let mut loaded = loaded.copy_to(annotated_path(&loaded.file, out_prefix))?;
    let dir = &loaded.dir;
    let rel = |p: &Option<String>| {
        p.as_deref()
            .map(|abs| crate::manifest::run::rel_to_manifest(dir, abs))
    };
    let a = &mut loaded.manifest.annotate;
    match pass {
        Pass::Markers(markers, called) => {
            // Stored like every other path here, relative to the manifest,
            // so it still resolves when read from another directory.
            a.markers = Some(crate::manifest::run::rel_to_manifest(dir, markers));
            a.argmax = rel(&out.argmax);
            a.annotation = rel(&out.annotation);
            a.marker_score = rel(&out.marker_score);
            a.units_call = rel(&out.units_call);
            a.cluster_celltype_q = rel(&out.cluster_celltype_q);
            a.cluster_celltype_es = rel(&out.cluster_celltype_es);
            a.cluster_celltype_pip = rel(&out.cluster_celltype_pip);
            a.cluster_celltype_effect = rel(&out.cluster_celltype_effect);
            a.cluster_celltype_explained = rel(&out.cluster_celltype_explained);
            a.called_by = called;
            a.cluster_expression = rel(&out.cluster_expression);
            a.ontology_assignment = rel(&out.ontology_assignment);
            a.ontology_node_mass = rel(&out.ontology_node_mass);
            a.cluster_celltype_q_values = rel(&out.cluster_celltype_q_values);
            a.cluster_celltype_p = rel(&out.cluster_celltype_p);
            a.cluster_celltype_nes = rel(&out.cluster_celltype_nes);
            // The retired term ORA's tables are no later round's.
            a.cluster_term_q = None;
            a.marker_embedding = None;
            a.celltype_units = rel(&out.celltype_units);
            a.celltype_links = rel(&out.celltype_links);
            a.celltype_membership = rel(&out.celltype_membership);
            a.term_alignment = rel(&out.term_alignment);
            a.term_units = rel(&out.term_units);
            a.unit_anchors = rel(&out.unit_anchors);
            a.ontology_signature = rel(&out.ontology_signature);
            a.ontology_term_effect = rel(&out.ontology_term_effect);
            loaded.manifest.defaults.colour_by = Some("annotation".into());
        }
        Pass::GeneSets => {
            a.cluster_expression = rel(&out.cluster_expression);
            a.ontology_signature = rel(&out.ontology_signature);
            a.ontology_term_effect = rel(&out.ontology_term_effect);
        }
        Pass::Ontology => {
            a.ontology_assignment = rel(&out.ontology_assignment);
            a.ontology_node_mass = rel(&out.ontology_node_mass);
        }
    }
    if let Some(c) = &out.clusters {
        let c = crate::manifest::run::rel_to_manifest(dir, c);
        if out.cluster_expression.is_some() {
            loaded.manifest.annotate.expression_clusters = Some(c.clone());
        }
        loaded.manifest.cluster.clusters = Some(c);
    }
    // The previous round's decisions are in its own log; this round made none.
    loaded.manifest.annotate.log = None;
    // A fresh pass: its statistics are not post-selection.
    loaded.manifest.annotate.stats = None;
    if matches!(pass, Pass::Markers(..)) {
        // Enrichment caches its statistics; projection has none to cache.
        loaded.manifest.annotate.stats_cache = out.stats_cache.clone().map(|c| {
            let rel = |p: &str| crate::manifest::run::rel_to_manifest(&loaded.dir, p);
            crate::manifest::run::StatsCache {
                gene_sum: rel(&c.gene_sum),
                batch_profile: rel(&c.batch_profile),
                gene_weight: rel(&c.gene_weight),
                cell_batch: rel(&c.cell_batch),
            }
        });
    }
    let recorded = loaded
        .manifest
        .annotate
        .settings
        .get_or_insert_with(|| serde_json::json!({}));
    if let serde_json::Value::Object(m) = recorded {
        m.insert(method.into(), settings);
    }
    // Derived from what was just written; a failure here should not cost the pass.
    if let Err(e) = rounds::write_summary(&mut loaded, out_prefix) {
        log::warn!("no cluster summary for this round: {e:#}");
    }
    loaded.manifest.save(&loaded.file)
}

/////////////////////////////////////
// manifest → enrichment inputs    //
/////////////////////////////////////

/// How to re-open the raw counts a run trained on: its `data.input` and
/// `data.batch`, found from the manifest's directory ([`RunManifest::data_file`]), under the
/// multiome layout it recorded — a multiome run's files are modalities of one
/// cell set, glued by barcode and namespaced as training did. No cell QC:
/// annotation maps onto the run's existing cells.
fn raw_counts_load(manifest: &RunManifest, manifest_dir: &Path) -> Result<ReadSharedRowsArgs> {
    let boxed =
        |v: Vec<String>| -> Vec<Box<str>> { v.into_iter().map(String::into_boxed_str).collect() };
    let data_files = boxed(manifest.data_inputs(manifest_dir));
    let batch_files =
        (!manifest.data.batch.is_empty()).then(|| boxed(manifest.data_batches(manifest_dir)));
    let mut args = ReadSharedRowsArgs {
        data_files,
        batch_files,
        keep_empty_barcodes: true,
        ..Default::default()
    };
    // Replay a multiome load: cells glue by barcode, features are namespaced
    // `{name}/{modality}`, and (across sample groups) barcodes `{barcode}@{group}`.
    // The record is positional against `data.input`.
    if let Some(r) = manifest.data.multiome.as_ref() {
        let n = args.data_files.len();
        anyhow::ensure!(
            r.modality.len() == n && r.group.len() == n,
            "the run's multiome layout covers {} file(s) but data.input lists {n}",
            r.modality.len()
        );
        args.column_alignment = data_beans::sparse_io_vector::ColumnAlignment::Union;
        args.feature_kind = Some(data_beans::aux::feature_names::FeatureNameKind::Mixed);
        args.per_file_feature_suffix = Some(r.modality.iter().map(|s| s.as_str().into()).collect());
        args.per_file_barcode_suffix = r
            .barcode_tagged
            .then(|| r.group.iter().map(|s| Some(s.as_str().into())).collect());
    }
    Ok(args)
}

/// Relabel clusters smaller than `min_size` to `usize::MAX` (unassigned) and
/// compact the rest to `0..k`; returns `k`.
fn drop_small_clusters(labels: &mut [usize], min_size: usize) -> usize {
    let min_size = min_size.max(1);
    let n = labels.iter().copied().max().map_or(0, |m| m + 1);
    let mut sizes = vec![0usize; n];
    for &l in labels.iter() {
        sizes[l] += 1;
    }
    let mut remap = vec![usize::MAX; n];
    let mut k = 0;
    for (old, &sz) in sizes.iter().enumerate() {
        if sz >= min_size {
            remap[old] = k;
            k += 1;
        }
    }
    if k < n {
        let dropped: usize = sizes.iter().filter(|&&s| s < min_size).sum();
        info!(
            "Removed {} cluster(s) with < {min_size} cells ({dropped} cells unassigned)",
            n - k
        );
        for l in labels.iter_mut() {
            *l = remap[*l];
        }
    }
    k
}

/// Re-open the raw counts, resolve the clustering, parse the marker TSV, and
/// aggregate the NB-Fisher-weighted cluster and per-batch expression the
/// enrichment pass scores.
pub(super) fn load_enrichment_inputs(
    args: &AnnotateArgs,
    loaded: &Loaded,
    data: Option<&super::data_files::ClData>,
) -> Result<EnrichmentInputs> {
    let (manifest, manifest_dir) = (&loaded.manifest, loaded.dir.as_path());
    anyhow::ensure!(
        !manifest.data.input.is_empty(),
        "manifest.data.input is empty; cannot re-open raw counts for cluster aggregation"
    );

    match super::no_counts::source(loaded) {
        super::no_counts::Source::Counts => {}
        other => return super::no_counts::inputs(args, loaded, data, other),
    }
    let load = raw_counts_load(manifest, manifest_dir)?;
    info!("Re-opening raw counts: {} file(s)", load.data_files.len());
    let SparseDataWithBatch {
        data: data_vec,
        batch,
        ..
    } = read_data_on_shared_rows(load)?;
    let data_vec = &data_vec;
    let n_cells = data_vec.num_columns();
    let n_genes = data_vec.num_rows();
    let cell_names = data_vec.column_names()?;
    let gene_names = data_vec.row_names()?;
    info!("Raw counts: {n_genes} genes × {n_cells} cells");

    let (cluster_labels, n_clusters) = resolve_clusters(
        args.clusters.as_deref(),
        &args.leiden(),
        manifest,
        manifest_dir,
        &cell_names,
    )?;
    anyhow::ensure!(
        n_clusters >= 2,
        "annotate needs ≥ 2 clusters, found {n_clusters}"
    );

    // Build per-cell batch ids for permutation null.
    let (batch_labels, n_batches) = build_batch_labels(batch, n_cells)?;
    info!("Batches: {n_batches}");

    let panel = panel_inputs(args, loaded, data, &gene_names)?;

    let nb_fisher = nb_fisher_weights(&loaded.model_prefix(), data_vec, &gene_names)?;
    let (profile_gk, pb_gene_gp, gene_sum_kg) = aggregate_expression(
        data_vec,
        &cluster_labels,
        n_clusters,
        &batch_labels,
        n_batches,
        gene_names.len(),
        &nb_fisher,
    )?;

    Ok(EnrichmentInputs {
        gene_names,
        cell_names,
        cluster_labels,
        n_clusters,
        batch_labels,
        n_batches,
        markers_gc: panel.markers_gc,
        marker_support: panel.marker_support,
        celltype_names: panel.celltype_names,
        profile_gk,
        pb_gene_gp,
        gene_sum_kg,
        gene_weights: nb_fisher,
        type_tree: panel.type_tree,
        cl_record: panel.cl_record,
        expression_source: None,
    })
}

/// The marker panel's side of [`EnrichmentInputs`], aligned to `gene_names`.
pub(super) struct PanelInputs {
    pub markers_gc: Mat,
    pub marker_support: Vec<Vec<usize>>,
    pub celltype_names: Vec<Box<str>>,
    pub type_tree: Option<enrichment::treebh::TypeTree>,
    pub cl_record: Option<serde_json::Value>,
}

/// The marker matrix over `gene_names`, its cell types and their tree.
pub(super) fn panel_inputs(
    args: &AnnotateArgs,
    loaded: &Loaded,
    data: Option<&super::data_files::ClData>,
    gene_names: &[Box<str>],
) -> Result<PanelInputs> {
    // Markers aligned to data row order. Optional: GO/GMT ontology mode supplies
    // gene-sets instead of a curated marker TSV, so an empty path yields an empty
    // marker matrix (the marker-enrichment path is skipped by the caller).
    let mut cl_record = None;
    let (markers_gc, marker_support, celltype_names, type_tree) = if args.markers.is_empty() {
        info!("No marker TSV (ontology gene-set mode); skipping marker matrix");
        (
            Mat::zeros(gene_names.len(), 0),
            vec![Vec::new(); gene_names.len()],
            Vec::new(),
            None,
        )
    } else {
        let annot = build_annotation_matrix(&args.markers, gene_names)?;
        info!(
            "Marker matrix: {} genes × {} celltypes",
            annot.membership_ga.nrows(),
            annot.membership_ga.ncols()
        );
        // The tree each cluster's q-values are TreeBH-adjusted over.
        let panel = crate::annotate::markers::read_panel(&args.markers)?;
        let own;
        let data = match data {
            Some(d) => d,
            None => {
                own = super::ontology::load(
                    Some(&loaded.dir),
                    &args.markers,
                    args.obo.as_deref(),
                    args.label_cl.as_deref(),
                    super::data_files::Fetch::Never,
                )?;
                &own
            }
        };
        let tree = super::ontology::panel_tree(data, &panel)?;
        cl_record = Some(data.record(tree.release.as_deref()));
        let type_tree = Some(tree.treebh(&annot.annot_names));
        (
            annot.membership_ga,
            annot.support,
            annot.annot_names,
            type_tree,
        )
    };

    Ok(PanelInputs {
        markers_gc,
        marker_support,
        celltype_names,
        type_tree,
        cl_record,
    })
}

/// Cluster source priority:
///   1. `--clusters <path>` — user-explicit, hard-fail on missing.
///   2. `manifest.cluster.clusters` if it exists on disk — softly
///      fall back to internal Leiden when stale / relocated, since
///      we can re-cluster from the same latent the manifest points at.
///   3. Internal Leiden on `manifest.outputs.latent`.
pub(super) fn resolve_clusters(
    clusters: Option<&str>,
    leiden: &LeidenArgs,
    manifest: &RunManifest,
    manifest_dir: &Path,
    cell_names: &[Box<str>],
) -> Result<(Vec<usize>, usize)> {
    let manifest_cluster_path = manifest
        .cluster
        .clusters
        .as_deref()
        .map(|rel| resolve(manifest_dir, rel));
    let resolved_path = clusters.map(String::from).or_else(|| {
        manifest_cluster_path.filter(|p| {
            let exists = Path::new(p).is_file();
            if !exists {
                log::warn!(
                    "Cluster parquet {p} not found — falling back to internal Leiden \
                         on the manifest's latent."
                );
            }
            exists
        })
    });
    if let Some(path) = resolved_path {
        info!("Resolving cluster source: parquet {path}");
        return load_cluster_labels(&path, cell_names)
            .with_context(|| format!("failed to load cluster parquet {path}"));
    }
    info!(
        "Resolving cluster source: internal Leiden on manifest's cell embedding ({:?})",
        manifest.outputs.geometry_latent()
    );
    compute_clusters_from_latent(manifest, manifest_dir, cell_names, leiden)
}

/// Per-gene NB-Fisher weights: the cached parquet from training first, falling
/// back to recomputing them off the counts.
fn nb_fisher_weights(
    fisher_prefix: &str,
    data_vec: &data_beans::sparse_io_vector::SparseIoVec,
    gene_names: &[Box<str>],
) -> Result<Vec<f32>> {
    use data_beans::alg::gene_weighting::{compute_nb_fisher_weights, load_fisher_weights};

    let nb_fisher: Vec<f32> = match load_fisher_weights(fisher_prefix)? {
        Some((cached_genes, cached_w)) if cached_genes == gene_names => {
            info!(
                "Loaded {} NB-Fisher weights from {fisher_prefix}.fisher_weights.parquet",
                cached_w.len()
            );
            cached_w
        }
        Some((cached_genes, _)) => {
            info!(
                "Cached fisher_weights gene names ({}) don't match data ({}); recomputing",
                cached_genes.len(),
                gene_names.len()
            );
            compute_nb_fisher_weights(data_vec, Some(BLOCK_SIZE))?
        }
        None => compute_nb_fisher_weights(data_vec, Some(BLOCK_SIZE))?,
    };
    let (w_min, w_max, w_sum) = nb_fisher.par_iter().map(|&w| (w, w, w)).reduce(
        || (f32::INFINITY, 0.0f32, 0.0f32),
        |(lo, hi, s), (a, b, c)| (lo.min(a), hi.max(b), s + c),
    );
    info!(
        "NB-Fisher weights: min={:.4}, max={:.4}, mean={:.4}",
        w_min,
        w_max,
        w_sum / nb_fisher.len() as f32
    );
    Ok(nb_fisher)
}

/// One fused sweep over the counts for the per-cluster gene sums, plus the
/// per-batch axis the sample-permutation null needs (cell types and GO terms
/// alike).
fn aggregate_expression(
    data_vec: &data_beans::sparse_io_vector::SparseIoVec,
    cluster_labels: &[usize],
    n_clusters: usize,
    batch_labels: &[usize],
    n_batches: usize,
    g: usize,
    nb_fisher: &[f32],
) -> Result<(Mat, Mat, Vec<f64>)> {
    let (gene_sum_kg, gene_sum_pg) = accumulate_gene_sum_pair(
        data_vec,
        cluster_labels,
        n_clusters,
        batch_labels,
        n_batches,
        g,
        BLOCK_SIZE,
    )?;
    Ok((
        weighted_mean_profile(&gene_sum_kg, n_clusters, g, nb_fisher),
        weighted_mean_profile(&gene_sum_pg, n_batches, g, nb_fisher),
        gene_sum_kg,
    ))
}

/// Run Leiden on the manifest's latent matrix. For topic kinds, log-space
/// latents are exponentiated to probabilities first. Aligns labels back to
/// the data column order.
pub fn compute_clusters_from_latent(
    manifest: &RunManifest,
    manifest_dir: &Path,
    cell_names: &[Box<str>],
    args: &LeidenArgs,
) -> Result<(Vec<usize>, usize)> {
    // Prepared as every kNN graph of a run is; Leiden's own normalisation
    // below then changes nothing (unit rows stay unit, z-scores stay z-scores).
    let latent = crate::manifest::run::prepare_geometry(manifest, manifest_dir)?;

    // An embedding is trained with a dot-product / cosine-style objective, so the
    // kNN graph fed to Leiden should reflect ANGULAR distance. Plain Euclidean on
    // a raw embedding is dominated by the depth axis — senna documents that for gem
    // about its own output — which is why the choice follows the space rather
    // than a hand-listed pair of kinds.
    let cosine = manifest.kind.cell_space() == CellSpace::Embedding;
    let metric = if cosine {
        "cosine"
    } else {
        "z-scored euclidean"
    };

    info!(
        "Internal Leiden: knn={}, resolution={:.3}, target_k={:?}, min_cluster_size={}, metric={metric}",
        args.knn,
        args.resolution,
        args.num_clusters,
        args.min_cluster_size
    );

    let mut leiden = legume_numeric::matrix::clustering::leiden_clustering(
        &latent.mat,
        args.knn,
        args.resolution,
        args.num_clusters,
        args.seed,
        cosine,
    )?;
    let n_clusters = drop_small_clusters(&mut leiden, args.min_cluster_size);

    // Align by name. The training pipeline writes latents in data column
    // order, so name-mismatched cases are rare but worth diagnosing
    // loudly — silent misalignment makes every downstream cluster
    // statistic empty.
    let mut idx: HashMap<&str, usize> = HashMap::default();
    idx.reserve(latent.rows.len());
    for (i, name) in latent.rows.iter().enumerate() {
        idx.insert(name.as_ref(), i);
    }
    let mut labels = Vec::with_capacity(cell_names.len());
    let mut missing = 0usize;
    let mut pruned = 0usize;
    for cell in cell_names {
        if let Some(&i) = idx.get(cell.as_ref()) {
            let lab = leiden[i];
            if lab == usize::MAX {
                pruned += 1;
            }
            labels.push(lab);
        } else {
            labels.push(usize::MAX);
            missing += 1;
        }
    }
    let assigned = cell_names.len() - missing - pruned;
    info!(
        "Cluster alignment: {assigned} assigned, {pruned} pruned by min_cluster_size, \
         {missing} missing from latent (n_clusters={n_clusters})"
    );
    if missing > 0 {
        let preview_data: Vec<&str> = cell_names
            .iter()
            .take(3)
            .map(std::convert::AsRef::as_ref)
            .collect();
        let preview_latent: Vec<&str> = latent
            .rows
            .iter()
            .take(3)
            .map(std::convert::AsRef::as_ref)
            .collect();
        log::warn!(
            "{missing}/{} cells missing from latent (data examples: {:?}; latent examples: {:?}). \
             Likely cause: latent.parquet was written by a different pipeline run with \
             different barcode formatting (e.g. with vs. without `@<basename>` suffix).",
            cell_names.len(),
            preview_data,
            preview_latent,
        );
    }
    if assigned == 0 {
        anyhow::bail!(
            "Internal Leiden produced 0 assigned cells (missing={missing}, pruned={pruned}). \
             Check that the latent.parquet matches the data backend's cell barcodes."
        );
    }
    Ok((labels, n_clusters))
}

/// Build per-cell batch ids as compact u32 indices in [0, `n_batches`).
fn build_batch_labels(labels: Vec<Box<str>>, n_cells: usize) -> Result<(Vec<usize>, usize)> {
    if labels.is_empty() {
        return Ok((vec![0; n_cells], 1));
    }
    anyhow::ensure!(
        labels.len() == n_cells,
        "batch labels {} ≠ cell count {}",
        labels.len(),
        n_cells
    );

    let names: Vec<&str> = labels.iter().map(AsRef::as_ref).collect();
    Ok(data_beans::alg::dc_poisson::compact_labels(&names))
}

#[cfg(test)]
#[path = "tests/all_methods.rs"]
pub(crate) mod all_methods_tests;
#[cfg(test)]
#[path = "tests/projection_score.rs"]
pub(crate) mod projection_score_tests;
#[cfg(test)]
#[path = "tests/projection_units.rs"]
pub(crate) mod projection_units_tests;
