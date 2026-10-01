//! Cluster-based annotation: marker-set enrichment on the per-cluster
//! expression matrix (NB-Fisher adjusted, re-aggregated from raw counts by the
//! caller — see [`crate::manifest::annotate`]).

use super::args::{AnnotateArgs, KEEP_IEA, MIN_CONFIDENCE, NUM_DRAWS};
use super::inputs::EnrichmentInputs;
use super::outputs::{clean_outputs, AnnotationOutputs, ENRICHMENT_OUTPUT_SUFFIXES};
use super::outputs::{
    write_cluster_tables, ARGMAX_TSV, CLUSTER_CELLTYPE_ES, CLUSTER_CELLTYPE_ES_STD,
    CLUSTER_CELLTYPE_NES, CLUSTER_CELLTYPE_P, CLUSTER_CELLTYPE_Q, CLUSTER_CELLTYPE_Q_VALUES,
    CLUSTER_CELLTYPE_Z, CLUSTER_TERM_NES, CLUSTER_TERM_P, CLUSTER_TERM_Q_VALUES,
};
use enrichment::{
    annotate, annotate_types, AnnotateConfig, AnnotateOutputs, GroupInputs, SpecificityMode,
    TypeScores,
};
use legume_numeric::matrix::common_io::mkdir_parent;
use legume_numeric::matrix::dense_mat_io::{axis_id_names, Mat};
use legume_numeric::matrix::traits::IoOps;
use log::info;
use rayon::prelude::*;

/// What the argument surface decides before any data is read: where the outputs
/// go, and which of the two gene-set pipelines runs. [`plan`] settles both, so
/// the caller can aggregate the right axes for [`run`].
pub struct EnrichmentPlan {
    pub out: Box<str>,
    /// `--go`/`--gaf`/`--gmt` without markers: the GO/GMT signature alone,
    /// no cell-type annotation.
    pub ontology_mode: bool,
    /// `--gaf`/`--gmt` with markers: the GO/GMT signature is scored too, on
    /// the marker pass's cluster profile.
    pub gene_sets_too: bool,
}

/// Validate the gene-set source flags and erase a previous run's artifacts.
pub fn plan(args: &AnnotateArgs) -> anyhow::Result<EnrichmentPlan> {
    let out = args.out.clone();
    mkdir_parent(&out)?;
    // --markers → curated cell-type annotation (+ optional inline CL ontology
    // via --obo/--label-cl); --gaf/--gmt → ontology gene-set signature
    // (cross-cluster-contrasted module score per cluster); both → both, on one
    // cluster profile.
    let markers = !args.markers.is_empty();
    let gene_sets = args.go || args.gaf.is_some() || args.gmt.is_some();
    anyhow::ensure!(
        markers || gene_sets,
        "a gene-set source is required: --markers, --go, --gaf, or --gmt"
    );
    anyhow::ensure!(
        !(args.gaf.is_some() && args.gmt.is_some()),
        "--gaf and --gmt are alternatives; pass one"
    );
    let ontology_mode = gene_sets && !markers;
    if gene_sets {
        anyhow::ensure!(
            args.go_obo.is_some(),
            "--go/--gaf/--gmt need the Gene Ontology (--go-obo, or the cached download) to name GO terms"
        );
    }
    if ontology_mode {
        anyhow::ensure!(
            args.label_cl.is_none(),
            "--label-cl is for --markers (curated CL); GO/GMT term ids are ontology ids already"
        );
    } else {
        anyhow::ensure!(
            args.obo.is_some() == args.label_cl.is_some(),
            "--obo and --label-cl must be given together to run inline ontology annotation \
             (got only one); omit both to skip it"
        );
    }
    if !args.no_clean {
        clean_outputs(&out, ENRICHMENT_OUTPUT_SUFFIXES);
    }
    Ok(EnrichmentPlan {
        out,
        ontology_mode,
        gene_sets_too: gene_sets && markers,
    })
}

/// Score the aggregated cluster expression against the marker panel (or, in
/// GO/GMT mode, against the gene sets) and write the artifacts under
/// `plan.out`. Returns their paths; recording them in a run manifest is the
/// caller's job.
/// Fewest batches the sample-permutation null may shuffle. It permutes whole
/// batches, so `P` batches allow only `P!` distinct orderings: with one or
/// two, null (cluster, type) pairs came out "significant" far above the
/// nominal rate on synthetic data, while four or more were calibrated.
pub const MIN_PERM_BATCHES: usize = 4;

/// The sample permutations to run for `n_batches`: the requested number, or
/// none when there are too few batches to permute. The gene-set
/// randomization null, which stays calibrated at any batch count, then
/// carries the test alone.
#[must_use]
pub fn sample_perm_draws(n_batches: usize, requested: usize) -> usize {
    if n_batches < MIN_PERM_BATCHES {
        0
    } else {
        requested
    }
}

/// The marker-path scoring of one enrichment pass: marker enrichment per
/// cluster against the gene-set null (and the sample-permutation null with
/// enough batches), and FDR.
/// Writes nothing; [`run`] writes the outputs, and a relabel round rescores
/// through this on its merged clusters and edited panel.
pub fn score(args: &AnnotateArgs, inputs: &EnrichmentInputs) -> anyhow::Result<AnnotateOutputs> {
    let (group, markers_gc, config) = prepare(args, inputs)?;
    info!(
        "Running cluster × marker enrichment: {} clusters × {} celltypes, \
         row-rand B={}, sample-perm B={}",
        inputs.n_clusters,
        inputs.celltype_names.len(),
        NUM_DRAWS,
        args.num_perm,
    );
    annotate(&group, &markers_gc, &inputs.celltype_names, &config)
}

/// [`score`] for the cell types `types` (indices into
/// `inputs.celltype_names`) alone: their columns exactly as `score` gives
/// them, with the weighted panel and config to finish the rows with
/// [`enrichment::adjust`] once the other types' p-values are known.
pub fn score_types(
    args: &AnnotateArgs,
    inputs: &EnrichmentInputs,
    types: &[usize],
) -> anyhow::Result<(TypeScores, Mat, AnnotateConfig)> {
    let (group, markers_gc, config) = prepare(args, inputs)?;
    info!(
        "Rescoring {} of {} celltype(s) over {} clusters",
        types.len(),
        inputs.celltype_names.len(),
        inputs.n_clusters,
    );
    let scores = annotate_types(&group, &markers_gc, &inputs.celltype_names, &config, types)?;
    Ok((scores, markers_gc, config))
}

/// The enrichment's inputs, weighted panel and config for `inputs`.
pub fn prepare(
    args: &AnnotateArgs,
    inputs: &EnrichmentInputs,
) -> anyhow::Result<(GroupInputs, Mat, AnnotateConfig)> {
    prepare_with(
        args,
        inputs,
        inputs.markers_gc.clone(),
        inputs.type_tree.clone(),
    )
}

/// [`prepare`] for any gene × set membership `sets_gc` (a marker panel, or
/// GO terms), its q-values TreeBH-adjusted over `type_tree` when given and
/// flat BH within the cluster otherwise.
fn prepare_with(
    args: &AnnotateArgs,
    inputs: &EnrichmentInputs,
    sets_gc: Mat,
    type_tree: Option<enrichment::treebh::TypeTree>,
) -> anyhow::Result<(GroupInputs, Mat, AnnotateConfig)> {
    let n_clusters = inputs.n_clusters;
    let n_batches = inputs.n_batches;
    let profile_gk = &inputs.profile_gk;
    // pb_membership[batch, cluster] = (# cells in batch with cluster id) / batch_size.
    let pb_membership_pk = build_pb_membership(
        &inputs.batch_labels,
        &inputs.cluster_labels,
        n_batches,
        n_clusters,
    );

    // One-hot cell membership (N × nClusters).
    let n_cells = inputs.cell_names.len();
    let mut cell_membership_nk = Mat::zeros(n_cells, n_clusters);
    for (n, &c) in inputs.cluster_labels.iter().enumerate() {
        if c < n_clusters {
            cell_membership_nk[(n, c)] = 1.0;
        }
    }

    // Data-aware specificity re-weighting: scale each marker by an empirical
    // specificity score derived from the actual cluster expression matrix.
    // Markers that light up broadly (shared by several types) get
    // attenuated; cluster-exclusive markers keep full weight. This
    // complements the IDF that already runs on the marker TSV.
    let mut markers_gc = sets_gc;
    apply_empirical_specificity_weights(&mut markers_gc, profile_gk);

    // The per-batch β̃ profile backs the sample-permutation null.
    let pb_gene_gp = inputs.pb_gene_gp.clone();
    let group = GroupInputs {
        profile_gk: profile_gk.clone(),
        pb_gene_gp,
        pb_membership_pk,
        cell_membership_nk,
        gene_names: inputs.gene_names.clone(),
        cell_names: inputs.cell_names.clone(),
    };

    let perm_draws = sample_perm_draws(n_batches, args.num_perm);
    if perm_draws < args.num_perm {
        log::warn!(
            "{n_batches} batch(es): too few for the sample-permutation null (needs \
             {MIN_PERM_BATCHES}); testing against the gene-set null alone"
        );
    }
    let config = AnnotateConfig {
        specificity: SpecificityMode::Simplex,
        num_row_randomization: NUM_DRAWS,
        num_sample_perm: perm_draws,
        // pb_membership_pk's rows ARE batches (one pseudobulk per batch),
        // so the sample-permutation null shuffles batches directly with no
        // inner stratification. Cell-level labels would be the wrong length
        // (caused a panic in `permute_indices`).
        batch_labels: None,
        fdr_alpha: args.fdr_alpha,
        q_softmax_temperature: args.q_temperature,
        min_confidence: MIN_CONFIDENCE,
        seed: args.seed,
        min_markers: args.min_markers,
        stratify_null: true,
        bootstrap: None,
        // The plain permutation p of the gene-set null's draws (floor 1 / (NUM_DRAWS + 1)).
        multilevel: None,
        // TreeBH over the panel's cell-type tree; flat BH within the cluster without one.
        type_tree,
    };

    Ok((group, markers_gc, config))
}

pub fn run(
    args: &AnnotateArgs,
    plan: &EnrichmentPlan,
    inputs: &EnrichmentInputs,
) -> anyhow::Result<AnnotationOutputs> {
    let out = plan.out.as_ref();
    let g = inputs.gene_names.len();
    let n_clusters = inputs.n_clusters;
    let profile_gk = &inputs.profile_gk;
    let cluster_names = axis_id_names("K", n_clusters);
    info!("Cluster expression: {g} genes × {n_clusters} clusters");
    let profile_max = profile_gk.iter().fold(0f32, |m, &v| m.max(v));
    if profile_max <= 1e-12 {
        anyhow::bail!(
            "Cluster expression matrix is all zero — every cell-axis cluster_label is \
             out of range (>= n_clusters). Check the cluster file's barcodes match the \
             data backend, or that the manifest's `clusters` path resolves correctly."
        );
    }

    ///////////////////////////////////
    // GO/GMT ontology gene-set mode //
    ///////////////////////////////////
    // The GO/GMT signature alone: no cell-type labels, so none of the marker
    // path's per-cell outputs.
    if plan.ontology_mode {
        return run_ontology_gene_sets(args, out, inputs, &cluster_names);
    }

    let AnnotateOutputs {
        q_kc,
        es_kc,
        es_restandardized_kc,
        nes_kc,
        z_kc,
        p_log2err_kc: _,
        perm_z_kc,
        pvalue_kc,
        qvalue_kc,
        cell_annotation_nc,
        argmax_labels,
        bootstrap: _,
    } = score(args, inputs)?;

    /////////////
    // Outputs //
    /////////////
    let cell_expr_path = format!("{out}.cluster_expression.parquet");
    profile_gk.to_parquet_with_names(
        &cell_expr_path,
        (Some(&inputs.gene_names), Some("gene")),
        Some(&cluster_names),
    )?;
    info!("wrote {cell_expr_path}");

    let annotation_path = format!("{out}.annotation.parquet");
    cell_annotation_nc.to_parquet_with_names(
        &annotation_path,
        (Some(&inputs.cell_names), Some("cell")),
        Some(&inputs.celltype_names),
    )?;
    info!("wrote {annotation_path}");

    // Per-cell label files via the SHARED writer (also emits `membership.tsv`),
    // so this pass and `annotate --method projection` produce an identical contract.
    let argmax_path = format!("{out}{ARGMAX_TSV}");
    {
        let cells: Vec<Box<str>> = argmax_labels
            .iter()
            .map(|l| Box::from(l.cell_name.as_ref()))
            .collect();
        let labels: Vec<Box<str>> = argmax_labels
            .iter()
            .map(|l| Box::from(l.label.as_ref()))
            .collect();
        let probs: Vec<f32> = argmax_labels.iter().map(|l| l.confidence).collect();
        graph_embedding_util::type_annotation::write_label_tsvs(out, &cells, &labels, &probs)?;
    }

    let written = write_cluster_tables(
        out,
        &cluster_names,
        &inputs.celltype_names,
        &[
            (&q_kc, CLUSTER_CELLTYPE_Q),
            (&es_kc, CLUSTER_CELLTYPE_ES),
            (&es_restandardized_kc, CLUSTER_CELLTYPE_ES_STD),
            (&nes_kc, CLUSTER_CELLTYPE_NES),
            (&z_kc, CLUSTER_CELLTYPE_Z),
            (&pvalue_kc, CLUSTER_CELLTYPE_P),
            (&qvalue_kc, CLUSTER_CELLTYPE_Q_VALUES),
        ],
    )?;
    let [q_path, es_path, _, nes_path, _, p_path, q_val_path] =
        <[String; 7]>::try_from(written).map_err(|_| anyhow::anyhow!("seven tables written"))?;
    // Correlation-preserving sample-permutation z (when num_perm > 0): the
    // preferred ontology input — graded, unlike the pooled p-value.
    if let Some(pz) = &perm_z_kc {
        write_cluster_tables(
            out,
            &cluster_names,
            &inputs.celltype_names,
            &[(pz, ".cluster_celltype_perm_z.parquet")],
        )?;
    }

    display_annotation_histogram(&cell_annotation_nc, &inputs.celltype_names);

    // Optional inline ontology annotation (TreeBH) — reuses the freshly computed
    // restandardized-ES z-matrix in memory, no parquet round-trip. NON-FATAL: a
    // bad label→CL map / OBO must not discard the already-written enrichment
    // outputs, so errors are logged and the run still finalizes.
    let mut ontology_assign: Option<String> = None;
    let mut ontology_mass: Option<String> = None;
    if let (Some(obo), Some(label_cl)) = (args.obo.as_deref(), args.label_cl.as_deref()) {
        match super::ontology::annotate_ontology_with_obo(
            out,
            label_cl,
            obo,
            args.ontology_fdr_q,
            args.ontology_by,
            // Prefer the correlation-preserving permutation z; fall back to the
            // row-randomization restandardized ES when no sample permutations ran.
            super::ontology::OntologyScore::Z(perm_z_kc.as_ref().unwrap_or(&es_restandardized_kc)),
            Some(&q_kc),
            &cluster_names,
            &inputs.celltype_names,
        ) {
            Ok((a, m)) => {
                ontology_assign = Some(a);
                ontology_mass = Some(m);
            }
            Err(e) => log::error!(
                "inline ontology annotation failed ({e}); enrichment outputs are intact"
            ),
        }
    }

    // GO/GMT terms on the same cluster profile. NON-FATAL, like the ontology
    // walk: the cell-type outputs above are already written.
    let mut gene_set_outputs = None;
    if plan.gene_sets_too {
        match gene_set_signature(args, out, inputs, &cluster_names) {
            Ok(paths) => gene_set_outputs = Some(paths),
            Err(e) => log::error!("GO term scoring failed ({e:#}); cell-type outputs are intact"),
        }
    }
    let (ontology_signature, ontology_term_effect) = gene_set_outputs.unzip();

    info!("annotate --method enrichment complete");
    Ok(AnnotationOutputs {
        ontology_signature,
        ontology_term_effect,
        cluster_celltype_q_values: Some(q_val_path),
        cluster_celltype_p: Some(p_path),
        cluster_celltype_nes: Some(nes_path),
        argmax: Some(argmax_path),
        annotation: Some(annotation_path),
        cluster_celltype_q: Some(q_path),
        cluster_celltype_es: Some(es_path),
        cluster_expression: Some(cell_expr_path),
        ontology_assignment: ontology_assign,
        ontology_node_mass: ontology_mass,
        ..AnnotationOutputs::default()
    })
}

/// GO/GMT ontology gene-set mode: the cluster profile and the GO/GMT signature
/// ([`gene_set_signature`]), with no cell-level labels (GO terms aren't cell
/// types).
fn run_ontology_gene_sets(
    args: &AnnotateArgs,
    out: &str,
    inputs: &EnrichmentInputs,
    cluster_names: &[Box<str>],
) -> anyhow::Result<AnnotationOutputs> {
    let (profile_gk, gene_names) = (&inputs.profile_gk, &inputs.gene_names);
    let cell_expr_path = format!("{out}.cluster_expression.parquet");
    profile_gk.to_parquet_with_names(
        &cell_expr_path,
        (Some(gene_names), Some("gene")),
        Some(cluster_names),
    )?;
    info!("wrote {cell_expr_path}");

    let (sig_path, effect_path) = gene_set_signature(args, out, inputs, cluster_names)?;

    info!("annotate --method enrichment (ontology gene-set mode) complete");
    Ok(AnnotationOutputs {
        cluster_expression: Some(cell_expr_path),
        ontology_signature: Some(sig_path),
        ontology_term_effect: Some(effect_path),
        ..AnnotationOutputs::default()
    })
}

/// Score the `--gaf`/`--gmt` terms on the cluster profile and write the
/// per-cluster signature and its `K × T` effect matrix; returns their paths.
/// Each term is also tested as the cell types are (fgsea NES, p, and q by BH
/// over the terms within each cluster), written as `K × T` tables beside it.
fn gene_set_signature(
    args: &AnnotateArgs,
    out: &str,
    inputs: &EnrichmentInputs,
    cluster_names: &[Box<str>],
) -> anyhow::Result<(String, String)> {
    use enrichment::ontology_module_score;
    let (profile_gk, gene_names) = (&inputs.profile_gk, &inputs.gene_names);

    let obo = args
        .go_obo
        .as_deref()
        .expect("gene-set scoring validated to have an ontology");
    // `--go`: the GO Consortium's annotations for the species the gene names
    // belong to.
    let gaf = match (&args.gaf, args.go && args.gmt.is_none()) {
        (Some(g), _) => Some(g.to_string()),
        (None, true) => {
            use super::go_signature::Species;
            let species = anyhow::Context::context(
                Species::detect(gene_names),
                "--go: cannot tell the species from the gene names; pass --gaf instead",
            )?;
            let gaf = crate::manifest::data_files::go_annotations(species, None)?;
            info!("GO annotations for {species:?} genes: {}", gaf.display());
            Some(gaf.to_string_lossy().into_owned())
        }
        (None, false) => None,
    };
    let gs = super::go_signature::load_go_gene_sets(
        obo,
        gaf.as_deref(),
        args.gmt.as_deref(),
        !KEEP_IEA,
        args.go_min_overlap,
        args.go_max_overlap,
        gene_names,
    )?;

    ////////////////////////////////////////////////////////////
    // Descriptive module-score signature (the GO/GMT scorer) //
    ////////////////////////////////////////////////////////////
    // Per (cluster, term): mean_in − mean_out of log1p(CP10K) on the cluster
    // profile, cross-cluster-contrasted. The top positive-effect terms per
    // cluster ARE the GO signature. This plain effect-size ranking recovers
    // cluster lineage (lineage and cell-cycle programs) more cleanly than a
    // permutation-z + TreeBH walk, whose ÷sd reweighting rewards small,
    // stable-null terms and whose depth preference descends to narrow processes.
    let ms = ontology_module_score(profile_gk, &gs.terms, &gs.universe)?;

    // The cell types' test on the same terms: their NES, p and q stand beside
    // the effect the signature ranks by.
    let mut sets_gt = Mat::zeros(gene_names.len(), gs.terms.len());
    for (t, (_, rows)) in gs.terms.iter().enumerate() {
        for &g in rows {
            sets_gt[(g, t)] = 1.0;
        }
    }
    // Taken by value: weighted in place, so the dense matrix is held once.
    let (group, weighted_gt, config) = prepare_with(args, inputs, sets_gt, None)?;
    info!(
        "Testing {} GO terms over {} clusters as the cell types are tested",
        gs.terms.len(),
        inputs.n_clusters
    );
    let all: Vec<usize> = (0..gs.terms.len()).collect();
    let tested = enrichment::annotate_types(&group, &weighted_gt, &ms.term_ids, &config, &all)?;
    let adjusted = enrichment::adjust(&tested.pvalue_kc, &weighted_gt, &config)?;
    write_cluster_tables(
        out,
        cluster_names,
        &ms.term_ids,
        &[
            (&tested.nes_kc, CLUSTER_TERM_NES),
            (&adjusted.pvalue_kc, CLUSTER_TERM_P),
            (&adjusted.qvalue_kc, CLUSTER_TERM_Q_VALUES),
        ],
    )?;
    let stats = super::go_signature::TermStats {
        nes_kt: &tested.nes_kc,
        p_kt: &adjusted.pvalue_kc,
        q_kt: &adjusted.qvalue_kc,
    };

    let sig_path = format!("{out}.ontology_signature.tsv");
    super::go_signature::write_go_signature(
        &sig_path,
        &gs.onto,
        &ms.effect_kt,
        &stats,
        &ms.term_ids,
        &gs.terms,
        "cluster",
        cluster_names,
    )?;
    let effect_path = format!("{out}.ontology_term_effect.parquet");
    ms.effect_kt.to_parquet_with_names(
        &effect_path,
        (Some(cluster_names), Some("cluster")),
        Some(&ms.term_ids),
    )?;
    info!("wrote {effect_path}");
    Ok((sig_path, effect_path))
}

/// Multiply each gene's marker entries by an empirical specificity score
/// derived from the cluster expression matrix.
///
/// For gene `g`, score = `(max_c μ[g,c]/Σ_c μ[g,c] − 1/K) / (1 − 1/K)` ∈
/// `[0, 1]`. A gene that fires evenly across all K clusters sits at the
/// uniform floor `1/K`, mapping to score 0. A gene exclusive to one
/// cluster has max simplex value 1, mapping to score 1.
fn apply_empirical_specificity_weights(markers_gc: &mut Mat, profile_gk: &Mat) {
    let g = markers_gc.nrows();
    let c = markers_gc.ncols();
    let k = profile_gk.ncols();
    debug_assert_eq!(profile_gk.nrows(), g);
    if k < 2 {
        return;
    }
    let inv_k = 1.0 / k as f32;
    let denom = (1.0 - inv_k).max(1e-8);

    let scores: Vec<f32> = (0..g)
        .into_par_iter()
        .map(|gi| {
            let row = profile_gk.row(gi);
            let sum: f32 = row.iter().sum();
            if sum <= 1e-12 {
                return 0.0;
            }
            let max = row.iter().fold(0.0f32, |m, &v| m.max(v));
            (((max / sum) - inv_k) / denom).clamp(0.0, 1.0)
        })
        .collect();

    let (mn, mx, sm) = scores
        .iter()
        .fold((f32::INFINITY, 0.0f32, 0.0f32), |(lo, hi, s), &x| {
            (lo.min(x), hi.max(x), s + x)
        });
    info!(
        "Empirical specificity weights: min={:.3}, max={:.3}, mean={:.3}",
        mn,
        mx,
        sm / g as f32
    );

    // Defensive guard: if the cluster expression matrix has no specificity
    // signal at all (every gene's row is uniform across clusters or all-zero
    // — typically because cluster_labels are misaligned and gene_sum_kg
    // ended up empty), multiplying by these all-zero scores would silently
    // zero the marker matrix, killing every enrichment downstream and
    // leaving every cell unassigned. Warn loudly and leave `markers_gc`
    // untouched so the caller can still get IDF-weighted enrichment.
    if mx <= 1e-6 {
        log::warn!(
            "Empirical specificity scores are all ~0 — likely cluster_labels misaligned or \
             gene_sum_kg is empty. Falling back to IDF-weighted markers without empirical \
             reweighting. Re-check that the cluster file's cell barcodes match the data, \
             or pass --no-empirical-specificity to silence this fallback."
        );
        return;
    }

    for gi in 0..g {
        let s = scores[gi];
        for ci in 0..c {
            markers_gc[(gi, ci)] *= s;
        }
    }
}

/// `pb_membership[b, c]` = fraction of batch `b`'s cells assigned to cluster `c`.
/// Rows that observe no cells stay all-zero, which the enrichment crate handles.
fn build_pb_membership(
    batch_labels: &[usize],
    cluster_labels: &[usize],
    n_batches: usize,
    n_clusters: usize,
) -> Mat {
    let mut out = Mat::zeros(n_batches, n_clusters);
    let mut batch_count = vec![0u64; n_batches];
    for (n, &b) in batch_labels.iter().enumerate() {
        if b >= n_batches {
            continue;
        }
        batch_count[b] += 1;
        let c = cluster_labels[n];
        if c < n_clusters {
            out[(b, c)] += 1.0;
        }
    }
    for b in 0..n_batches {
        let s = batch_count[b].max(1) as f32;
        for c in 0..n_clusters {
            out[(b, c)] /= s;
        }
    }
    out
}

fn display_annotation_histogram(annot: &Mat, annot_names: &[Box<str>]) {
    let n_cells = annot.nrows();
    let n_types = annot.ncols();

    // Per-cell argmax: rows are independent, so this fans out cleanly.
    let per_cell: Vec<(f32, Option<usize>)> = (0..n_cells)
        .into_par_iter()
        .map(|i| {
            let row = annot.row(i);
            let sum: f32 = row.iter().sum();
            if sum < 1e-12 {
                return (0.0, None);
            }
            let (idx, val) = row
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.total_cmp(b))
                .unwrap();
            (*val, Some(idx))
        })
        .collect();

    let mut type_counts = vec![0usize; n_types];
    let mut type_prob_sum = vec![0.0f32; n_types];
    let mut unassigned = 0usize;
    for &(prob, assign) in &per_cell {
        match assign {
            Some(c) => {
                type_counts[c] += 1;
                type_prob_sum[c] += prob;
            }
            None => unassigned += 1,
        }
    }
    let mut sorted_types: Vec<usize> = (0..n_types).collect();
    sorted_types.sort_by(|&a, &b| type_counts[b].cmp(&type_counts[a]));

    let max_count = *type_counts.iter().max().unwrap_or(&1).max(&unassigned);
    const MAX_BAR: usize = 20;

    let assigned_cells = n_cells - unassigned;
    let assigned_prob_sum: f32 = per_cell
        .iter()
        .filter_map(|(p, a)| a.map(|_| *p))
        .sum::<f32>();
    let mean_prob = if assigned_cells > 0 {
        assigned_prob_sum / assigned_cells as f32
    } else {
        0.0
    };
    let above_50 = per_cell.iter().filter(|(p, _)| *p > 0.5).count();
    let above_70 = per_cell.iter().filter(|(p, _)| *p > 0.7).count();

    eprintln!();
    eprintln!("Annotation Summary ({n_cells} cells)");
    eprintln!(
        "  Mean max-prob (assigned): {:.3}  >0.5: {} ({:.1}%)  >0.7: {} ({:.1}%)",
        mean_prob,
        above_50,
        100.0 * above_50 as f32 / n_cells as f32,
        above_70,
        100.0 * above_70 as f32 / n_cells as f32
    );
    if unassigned > 0 {
        let bar_len = (unassigned * MAX_BAR) / max_count.max(1);
        eprintln!(
            "  {:24} {:5} ({:5.1}%)      {}",
            "unassigned",
            unassigned,
            100.0 * unassigned as f32 / n_cells as f32,
            "▒".repeat(bar_len)
        );
    }
    eprintln!();

    for &ct in &sorted_types {
        if type_counts[ct] == 0 {
            continue;
        }
        let bar_len = (type_counts[ct] * MAX_BAR) / max_count.max(1);
        let bar: String = "█".repeat(bar_len);
        let avg_prob = type_prob_sum[ct] / type_counts[ct] as f32;
        eprintln!(
            "  {:24} {:5} ({:5.1}%) {:.2} {}",
            annot_names[ct],
            type_counts[ct],
            100.0 * type_counts[ct] as f32 / n_cells as f32,
            avg_prob,
            bar
        );
    }
    eprintln!();
}

#[cfg(test)]
mod perm_tests {
    use super::*;

    #[test]
    fn too_few_batches_skip_the_sample_permutation() {
        assert_eq!(sample_perm_draws(1, 500), 0);
        assert_eq!(sample_perm_draws(MIN_PERM_BATCHES - 1, 500), 0);
        assert_eq!(sample_perm_draws(MIN_PERM_BATCHES, 500), 500);
        assert_eq!(sample_perm_draws(1, 0), 0, "none requested, none run");
    }
}
