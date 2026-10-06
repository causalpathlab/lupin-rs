//! Run-manifest glue for the annotation commands: resolve `-f run.senna.json`
//! into the loaded matrices [`crate::annotate`] takes, then write the artifact
//! paths it returns back into `manifest.annotate.*`. This is the only place
//! that touches [`RunManifest`] on behalf of an annotation command.
//!
//! Re-opening the raw counts and aggregating them per cluster also lives here,
//! over the sparse count backend; [`crate::annotate`] receives the resulting
//! cluster × gene expression.

use crate::annotate::args::{
    fixed_settings, AnnotateArgs, AnnotateOntologyArgs, AnnotateProjectionArgs, BLOCK_SIZE,
};
use crate::annotate::by_enrichment;
use crate::annotate::by_projection::{self, ProjectionInputs};
use crate::annotate::inputs::{load_cluster_labels, EnrichmentInputs};
use crate::annotate::ontology;
use crate::annotate::outputs::{AnnotationOutputs, ANNOT_PARQUET};
use crate::manifest::rounds;
use crate::manifest::run::{annotated_path, resolve, Loaded};
use legume_numeric::matrix::parquet::read_table_columns;

use crate::annotate::aggregate::{accumulate_gene_sum_pair, weighted_mean_profile};
use crate::annotate::markers::build_annotation_matrix;
use crate::manifest::run::{CellSpace, RunManifest};
use crate::marker_embedding::load_marker_feature_embedding;
use data_beans::aux::data_loading::{
    read_data_on_shared_rows, ReadSharedRowsArgs, SparseDataWithBatch,
};
use legume_numeric::matrix::dense_mat_io::Mat;

use anyhow::{anyhow, Context, Result};
use legume_numeric::matrix::dmatrix_io::DMatrix;
use legume_numeric::matrix::traits::IoOps;
use log::info;
use rayon::prelude::*;
use rustc_hash::FxHashMap as HashMap;
use std::path::Path;

/// CLI passthrough for the internal Leiden fallback (used only when neither
/// `--clusters` nor `manifest.cluster.clusters` is provided).
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
    // The ids the cluster tables' `K{id}` rows refer to, so later rounds and
    // viewers key on the same clusters.
    let clusters_path = format!("{}{}", args.out, rounds::CLUSTERS);
    let ids: Vec<Option<u32>> = inputs
        .cluster_labels
        .iter()
        .map(|&k| u32::try_from(k).ok())
        .collect();
    rounds::write_clusters(&clusters_path, &inputs.cell_names, &ids)?;
    outputs.clusters = Some(clusters_path);
    // The decoder's expected expression is no count sum: a later pass or
    // rescoring must not take it for one.
    let from_decoder = inputs
        .expression_source
        .as_ref()
        .is_some_and(|e| e["from"] == "decoder");
    if !plan.ontology_mode && !from_decoder {
        outputs.stats_cache = Some(crate::manifest::recalibrate::write_cache(
            &args.out, &inputs,
        )?);
    }
    let pass = if plan.ontology_mode {
        Pass::GeneSets
    } else {
        Pass::Markers(&args.markers)
    };
    // Which nulls actually ran: the sample permutation is skipped when there
    // are too few batches to shuffle.
    let mut used = settings(args)?;
    if let serde_json::Value::Object(m) = &mut used {
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

/// `lupin annotate --method projection`: score the run's co-embedded gene
/// space and cell embedding against the marker panel, and record what it wrote.
pub fn annotate_by_projection(args: &AnnotateProjectionArgs, loaded: &Loaded) -> Result<()> {
    let run = loaded.file.display().to_string();
    // Genes on the cell manifold; for a `gem` run only the spliced rows,
    // re-keyed by gene — see [`crate::marker_embedding`].
    let feat = load_marker_feature_embedding(&loaded.manifest, &loaded.dir, &run).with_context(
        || {
            "projection needs a co-embedded gene space (a senna `gem` / `bge` / `fne` / \
             `resolve-embedding-space` run). For topic/svd runs use `lupin annotate --method enrichment`."
        },
    )?;
    let cell_path = loaded.geometry_latent_path()?;
    let cell = DMatrix::<f32>::from_parquet(&cell_path)
        .with_context(|| format!("reading cell embedding {cell_path}"))?;

    let mut outputs = by_projection::run(
        args,
        &ProjectionInputs {
            feature_embedding: &feat,
            cell_embedding: &cell,
        },
    )?;
    outputs.clusters = Some(projection_clusters(&args.out)?);
    record(
        loaded,
        &args.out,
        Pass::Markers(&args.markers),
        &outputs,
        "projection",
        settings(args)?,
    )
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
    /// Per-cell cell-type labels from a marker panel (enrichment or projection).
    Markers(&'a str),
    /// GO/GMT gene-set signature per cluster.
    GeneSets,
    /// CL ontology walk over an earlier enrichment.
    Ontology,
}

/// Projection clusters the cells itself; its per-cell `community` becomes
/// the round's cluster table.
fn projection_clusters(out: &str) -> Result<String> {
    let annot = format!("{out}{ANNOT_PARQUET}");
    let (strings, nums) = read_table_columns(&annot, &["cell"], &["community"])
        .with_context(|| format!("reading {annot}"))?;
    let ids: Vec<Option<u32>> = nums[0]
        .iter()
        .map(|&k| (k.is_finite() && k >= 0.0).then_some(k as u32))
        .collect();
    let path = format!("{out}{}", rounds::CLUSTERS);
    rounds::write_clusters(&path, &strings[0], &ids)?;
    Ok(path)
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
        Pass::Markers(markers) => {
            // Stored like every other path here, relative to the manifest,
            // so it still resolves when read from another directory.
            a.markers = Some(crate::manifest::run::rel_to_manifest(dir, markers));
            a.argmax = rel(&out.argmax);
            a.annotation = rel(&out.annotation);
            a.cluster_celltype_q = rel(&out.cluster_celltype_q);
            a.cluster_celltype_es = rel(&out.cluster_celltype_es);
            a.cluster_expression = rel(&out.cluster_expression);
            a.ontology_assignment = rel(&out.ontology_assignment);
            a.ontology_node_mass = rel(&out.ontology_node_mass);
            a.cluster_celltype_q_values = rel(&out.cluster_celltype_q_values);
            a.cluster_celltype_p = rel(&out.cluster_celltype_p);
            a.cluster_celltype_nes = rel(&out.cluster_celltype_nes);
            a.cluster_term_q = rel(&out.cluster_term_q);
            a.marker_embedding = rel(&out.marker_embedding);
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
    if matches!(pass, Pass::Markers(_)) {
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

    let (cluster_labels, n_clusters) = resolve_clusters(args, manifest, manifest_dir, &cell_names)?;
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
    let (markers_gc, celltype_names, type_tree) = if args.markers.is_empty() {
        info!("No marker TSV (ontology gene-set mode); skipping marker matrix");
        (Mat::zeros(gene_names.len(), 0), Vec::new(), None)
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
        (annot.membership_ga, annot.annot_names, type_tree)
    };

    Ok(PanelInputs {
        markers_gc,
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
    args: &AnnotateArgs,
    manifest: &RunManifest,
    manifest_dir: &Path,
    cell_names: &[Box<str>],
) -> Result<(Vec<usize>, usize)> {
    let manifest_cluster_path = manifest
        .cluster
        .clusters
        .as_deref()
        .map(|rel| resolve(manifest_dir, rel));
    let resolved_path = args.clusters.as_deref().map(String::from).or_else(|| {
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
    let leiden_args = LeidenArgs {
        knn: args.knn,
        resolution: args.resolution,
        num_clusters: args.num_clusters,
        min_cluster_size: args.min_cluster_size,
        seed: args.cluster_seed,
    };
    compute_clusters_from_latent(manifest, manifest_dir, cell_names, &leiden_args)
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
