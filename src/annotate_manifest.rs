//! Run-manifest glue for the annotation commands: resolve `-f run.senna.json`
//! into the loaded matrices [`crate::annotate`] takes, then write the artifact
//! paths it returns back into `manifest.annotate.*`. This is the only place
//! that touches [`RunManifest`] on behalf of an annotation command.
//!
//! Re-opening the raw counts and aggregating them per cluster also lives here,
//! over the sparse count backend; [`crate::annotate`] receives the resulting
//! cluster × gene expression.

use crate::annotate::args::{AnnotateArgs, AnnotateOntologyArgs, AnnotateProjectionArgs};
use crate::annotate::by_enrichment::{self, EnrichmentPlan};
use crate::annotate::by_projection::{self, ProjectionInputs};
use crate::annotate::inputs::{load_cluster_labels, EnrichmentInputs};
use crate::annotate::ontology;
use crate::annotate::outputs::AnnotationOutputs;
use crate::run_manifest::{self, resolve, Loaded};

use crate::annotate::aggregate::{
    accumulate_gene_sum, accumulate_gene_sum_pair, weighted_mean_profile,
};
use crate::annotate::markers::build_annotation_matrix;
use crate::marker_embedding::load_marker_feature_embedding;
use crate::run_manifest::{CellSpace, RunManifest};
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
pub fn annotate_by_enrichment(args: &AnnotateArgs, loaded: &mut Loaded) -> Result<()> {
    let plan = by_enrichment::plan(args)?;
    let inputs = load_enrichment_inputs(args, &plan, loaded)?;
    let outputs = by_enrichment::run(args, &plan, &inputs)?;
    let pass = if plan.ontology_mode {
        Pass::GeneSets
    } else {
        Pass::Markers(&args.markers)
    };
    record(loaded, pass, &outputs)
}

/// `lupin annotate --method projection`: score the run's co-embedded gene
/// space and cell embedding against the marker panel, and record what it wrote.
pub fn annotate_by_projection(args: &AnnotateProjectionArgs, loaded: &mut Loaded) -> Result<()> {
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

    let outputs = by_projection::run(
        args,
        &ProjectionInputs {
            feature_embedding: &feat,
            cell_embedding: &cell,
        },
    )?;
    record(loaded, Pass::Markers(&args.markers), &outputs)
}

/// `lupin annotate` without markers: walk the CL tree over the cluster ×
/// cell-type matrix an earlier enrichment run recorded.
pub fn annotate_ontology(args: &AnnotateOntologyArgs, loaded: &mut Loaded) -> Result<()> {
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
    record(loaded, Pass::Ontology, &outputs)
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

/// Record a pass's artifacts in the manifest, relative to its directory, and
/// save it back to the file it was read from.
fn record(loaded: &mut Loaded, pass: Pass<'_>, out: &AnnotationOutputs) -> Result<()> {
    let dir = &loaded.dir;
    let rel = |p: &Option<String>| {
        p.as_deref()
            .map(|abs| run_manifest::rel_to_manifest(dir, abs))
    };
    let a = &mut loaded.manifest.annotate;
    match pass {
        Pass::Markers(markers) => {
            a.markers = Some(markers.to_string());
            a.argmax = rel(&out.argmax);
            a.annotation = rel(&out.annotation);
            a.cluster_celltype_q = rel(&out.cluster_celltype_q);
            a.cluster_celltype_es = rel(&out.cluster_celltype_es);
            a.cluster_expression = rel(&out.cluster_expression);
            a.ontology_assignment = rel(&out.ontology_assignment);
            a.ontology_node_mass = rel(&out.ontology_node_mass);
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
    loaded.manifest.save(&loaded.file)
}

/////////////////////////////////////
// manifest → enrichment inputs    //
/////////////////////////////////////

/// How to re-open the raw counts a run trained on: its `data.input` and
/// `data.batch`, resolved against the manifest's directory, under the
/// multiome layout it recorded — a multiome run's files are modalities of one
/// cell set, glued by barcode and namespaced as training did. No cell QC:
/// annotation maps onto the run's existing cells.
fn raw_counts_load(
    manifest: &RunManifest,
    manifest_dir: &Path,
    preload: bool,
) -> Result<ReadSharedRowsArgs> {
    let to_box = |s: &String| resolve(manifest_dir, s).into_boxed_str();
    let data_files: Vec<Box<str>> = manifest.data.input.iter().map(to_box).collect();
    let batch_files =
        (!manifest.data.batch.is_empty()).then(|| manifest.data.batch.iter().map(to_box).collect());
    let mut args = ReadSharedRowsArgs {
        data_files,
        batch_files,
        preload,
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
/// aggregate the NB-Fisher-weighted cluster (and, for the marker path,
/// per-batch) expression the enrichment pass scores.
fn load_enrichment_inputs(
    args: &AnnotateArgs,
    plan: &EnrichmentPlan,
    loaded: &Loaded,
) -> Result<EnrichmentInputs> {
    let (manifest, manifest_dir) = (&loaded.manifest, loaded.dir.as_path());
    anyhow::ensure!(
        !manifest.data.input.is_empty(),
        "manifest.data.input is empty; cannot re-open raw counts for cluster aggregation"
    );

    let load = raw_counts_load(manifest, manifest_dir, args.preload_data)?;
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

    // Markers aligned to data row order. Optional: GO/GMT ontology mode supplies
    // gene-sets instead of a curated marker TSV, so an empty path yields an empty
    // marker matrix (the marker-enrichment path is skipped by the caller).
    let (markers_gc, celltype_names) = if args.markers.is_empty() {
        info!("No marker TSV (ontology gene-set mode); skipping marker matrix");
        (Mat::zeros(gene_names.len(), 0), Vec::new())
    } else {
        let annot = build_annotation_matrix(&args.markers, &gene_names)?;
        info!(
            "Marker matrix: {} genes × {} celltypes",
            annot.membership_ga.nrows(),
            annot.membership_ga.ncols()
        );
        (annot.membership_ga, annot.annot_names)
    };

    let nb_fisher = nb_fisher_weights(args, &loaded.run_prefix(), data_vec, &gene_names)?;
    let (profile_gk, pb_gene_gp) = aggregate_expression(
        args,
        plan,
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
        markers_gc,
        celltype_names,
        profile_gk,
        pb_gene_gp,
    })
}

/// Cluster source priority:
///   1. `--clusters <path>` — user-explicit, hard-fail on missing.
///   2. `manifest.cluster.clusters` if it exists on disk — softly
///      fall back to internal Leiden when stale / relocated, since
///      we can re-cluster from the same latent the manifest points at.
///   3. Internal Leiden on `manifest.outputs.latent`.
fn resolve_clusters(
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
    args: &AnnotateArgs,
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
            compute_nb_fisher_weights(data_vec, Some(args.block_size))?
        }
        None => compute_nb_fisher_weights(data_vec, Some(args.block_size))?,
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
/// per-batch axis the marker path's sample-permutation null needs. The GO/GMT
/// ontology path scores the per-cluster profile directly, so it asks for only
/// the cluster sums.
#[allow(clippy::too_many_arguments)]
fn aggregate_expression(
    args: &AnnotateArgs,
    plan: &EnrichmentPlan,
    data_vec: &data_beans::sparse_io_vector::SparseIoVec,
    cluster_labels: &[usize],
    n_clusters: usize,
    batch_labels: &[usize],
    n_batches: usize,
    g: usize,
    nb_fisher: &[f32],
) -> Result<(Mat, Option<Mat>)> {
    if plan.ontology_mode {
        let gene_sum_kg =
            accumulate_gene_sum(data_vec, cluster_labels, n_clusters, g, args.block_size)?;
        // μ[g, c] = w_NBF[g] · (Σ counts[g, n ∈ c]) / size_sum[c]; Simplex
        // specificity downstream supplies the cross-cluster housekeeping
        // suppression.
        return Ok((
            weighted_mean_profile(&gene_sum_kg, n_clusters, g, nb_fisher),
            None,
        ));
    }

    let (gene_sum_kg, gene_sum_pg) = accumulate_gene_sum_pair(
        data_vec,
        cluster_labels,
        n_clusters,
        batch_labels,
        n_batches,
        g,
        args.block_size,
    )?;
    Ok((
        weighted_mean_profile(&gene_sum_kg, n_clusters, g, nb_fisher),
        Some(weighted_mean_profile(&gene_sum_pg, n_batches, g, nb_fisher)),
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
    let latent_rel = manifest
        .outputs
        .geometry_latent()
        .ok_or_else(|| anyhow!("manifest missing outputs.cell_embedding and outputs.latent"))?;
    let latent_path = resolve(manifest_dir, latent_rel);
    let mut latent = Mat::from_parquet_with_row_names(&latent_path, Some(0))
        .with_context(|| format!("failed to load latent {latent_path}"))?;
    info!(
        "Loaded latent {latent_path}: {}×{}",
        latent.mat.nrows(),
        latent.mat.ncols()
    );

    // Gate on the space of the table just READ, not on what `latent.parquet`
    // would hold. A kind's `cell_space()` and its `latent_is_log_simplex()` are
    // separate declarations that do not have to agree, so keying this on
    // `latent_is_log_simplex` risks exponentiating a Euclidean embedding that
    // `geometry_latent` handed back instead of the simplex.
    //
    // This used to sniff `max <= 0.0` instead, which mis-handles masked-vae (raw
    // Gaussian `z`, usually max > 0 so it skipped by luck) and is undefined on an
    // all-NaN latent, where `max()`'s partial ordering decides the branch.
    if manifest.kind.cell_space() == CellSpace::LogSimplex {
        info!(
            "Log-simplex latent ({}); exponentiating to probabilities",
            manifest.kind
        );
        latent.mat.apply(|x| *x = x.exp());
    }

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
