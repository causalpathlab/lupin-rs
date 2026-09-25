//! Unified `lupin annotate` — enrichment, or projection from a run manifest or embedding files.

use crate::annotate::args::{AnnotateArgs, AnnotateOntologyArgs, AnnotateProjectionArgs};
use crate::annotate::by_projection::{self, ProjectionInputs};
use crate::annotate_manifest::{annotate_by_enrichment, annotate_by_projection, annotate_ontology};
use crate::run_manifest::{self, RunManifest};
use anyhow::{Context, Result};
use clap::{Args, ValueEnum};
use legume_numeric::matrix::dmatrix_io::DMatrix;
use legume_numeric::matrix::traits::IoOps;

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, ValueEnum)]
pub enum AnnotateMethod {
    #[default]
    Auto,
    Enrichment,
    Projection,
}

#[derive(Args, Debug)]
pub struct AnnotateCliArgs {
    #[arg(
        long,
        value_enum,
        default_value_t = AnnotateMethod::Auto,
        help = "Annotation backend: enrichment (topic/svd), projection (co-embed), or auto"
    )]
    pub method: AnnotateMethod,

    #[arg(
        long,
        short = 'f',
        help = "Run manifest (`run.senna.json`) or its output prefix"
    )]
    pub from: Option<Box<str>>,

    #[arg(
        long,
        help = "Feature × D embedding parquet (projection from files; needs --cell-embedding)"
    )]
    pub feature_embedding: Option<Box<str>>,

    #[arg(
        long,
        help = "Cell × D embedding parquet (projection from files; needs --feature-embedding)"
    )]
    pub cell_embedding: Option<Box<str>>,

    #[arg(
        long,
        short = 'm',
        default_value = "",
        help = "Marker TSV (`gene<TAB>celltype`); omit for ontology-only follow-up on enrichment output"
    )]
    pub markers: Box<str>,

    #[arg(long, short = 'o', help = "Output prefix")]
    pub out: Option<Box<str>>,

    // ── shared clustering / stats ──
    // Default is method-specific when omitted: enrichment 15, projection/ORA 30.
    #[arg(long)]
    pub knn: Option<usize>,
    #[arg(long, default_value_t = 1.0)]
    pub resolution: f64,
    #[arg(long, default_value_t = 42)]
    pub seed: u64,
    #[arg(long = "num-perm", default_value_t = 500)]
    pub num_perm: usize,
    #[arg(long = "min-markers", default_value_t = 3)]
    pub min_markers: usize,
    #[arg(long = "fdr-alpha", default_value_t = 0.1)]
    pub fdr_alpha: f32,
    #[arg(long = "q-temperature", default_value_t = 1.0)]
    pub q_temperature: f32,
    #[arg(long = "no-clean")]
    pub no_clean: bool,

    // ── enrichment-only ──
    #[arg(long = "clusters")]
    pub clusters: Option<Box<str>>,
    #[arg(long = "num-clusters")]
    pub num_clusters: Option<usize>,
    #[arg(long = "min-cluster-size", default_value_t = 2)]
    pub min_cluster_size: usize,
    #[arg(long = "cluster-seed")]
    pub cluster_seed: Option<u64>,
    #[arg(long = "gaf")]
    pub gaf: Option<Box<str>>,
    #[arg(long = "gmt")]
    pub gmt: Option<Box<str>>,

    // ── projection / ORA ──
    #[arg(long = "no-idf")]
    pub no_idf: bool,
    #[arg(long = "no-assign-qc")]
    pub no_assign_qc: bool,
    #[arg(long = "assign-mad", default_value_t = 2.5)]
    pub assign_mad: f64,

    // ── marker stability bootstrap (enrichment and projection) ──
    #[arg(
        long = "n-boot",
        default_value_t = 200,
        help = "Marker-panel bootstrap resamples behind each call's support"
    )]
    pub n_boot: usize,
    #[arg(
        long = "min-support",
        default_value_t = 0.5,
        help = "Minimum fraction of resamples the top label must win to be called",
        long_help = "Minimum fraction of resamples the top label must win.\n\
                     Below this bar the cell is not called at all.\n\
                     \n\
                     Not scale-free: with C types chance agreement is 1/C,\n\
                     so 0.5 is ~3x chance on a 6-type panel and ~12x on a 24-type one.\n\
                     --abstain-separable uses a sign test instead."
    )]
    pub min_support: f32,
    #[arg(
        long = "abstain-separable",
        help = "Abstain by a sign test (top vs runner-up) instead of --min-support",
        long_help = "Keep the top label only if it beat the runner-up by more than\n\
                     resampling noise: an exact binomial sign test at --abstain-alpha.\n\
                     Means the same thing at any number of types, which --min-support does not."
    )]
    pub abstain_separable: bool,
    #[arg(
        long = "abstain-alpha",
        default_value_t = 0.05,
        help = "[--abstain-separable] Significance level of the sign test"
    )]
    pub abstain_alpha: f64,
    #[arg(
        long = "set-coverage",
        default_value_t = 0.8,
        help = "Coverage of the reported `label_set`: the smallest label set covering this share of resamples"
    )]
    pub set_coverage: f32,
    #[arg(
        long = "max-set-size",
        default_value_t = 3,
        help = "Largest `label_set` reported; a cell needing more to reach --set-coverage is unassigned"
    )]
    pub max_set_size: usize,

    // ── ontology (inline on annotate, or follow-up without markers) ──
    #[arg(long = "obo")]
    pub obo: Option<Box<str>>,
    #[arg(long = "label-cl")]
    pub label_cl: Option<Box<str>>,
    #[arg(long = "ontology-fdr-q", default_value_t = 0.1)]
    pub ontology_fdr_q: f64,
    #[arg(long = "ontology-by", help = "TreeBH Benjamini–Yekutieli")]
    pub ontology_by: bool,
    #[arg(long = "use-perm-p")]
    pub use_perm_p: bool,
}

pub fn run_annotate(args: &AnnotateCliArgs) -> Result<()> {
    if is_ontology_followup(args) {
        return run_ontology_followup(args);
    }

    match route(args)? {
        Route::Enrichment => run_enrichment(args),
        Route::Projection => run_projection(args),
        Route::EmbeddingFiles { feat, cell, prefix } => {
            run_projection_from_files(args, &feat, &cell, &prefix)
        }
    }
}

fn is_ontology_followup(args: &AnnotateCliArgs) -> bool {
    args.markers.is_empty()
        && args.gaf.is_none()
        && args.gmt.is_none()
        && args.from.is_some()
        && args.obo.is_some()
        && args.label_cl.is_some()
}

fn run_ontology_followup(args: &AnnotateCliArgs) -> Result<()> {
    let from = args
        .from
        .as_ref()
        .context("--from required for ontology follow-up")?;
    let obo = args.obo.as_ref().context("--obo required")?;
    let label_cl = args.label_cl.as_ref().context("--label-cl required")?;
    annotate_ontology(&AnnotateOntologyArgs {
        from: from.clone(),
        label_cl: label_cl.clone(),
        obo: obo.clone(),
        out: args.out.clone(),
        fdr_q: args.ontology_fdr_q,
        by: args.ontology_by,
        use_perm_p: args.use_perm_p,
    })
}

enum Route {
    Enrichment,
    /// Co-embed projection through the run manifest.
    Projection,
    /// Projection over an explicit `--feature-embedding` / `--cell-embedding` pair.
    EmbeddingFiles {
        feat: String,
        cell: String,
        prefix: String,
    },
}

/// Pick the backend. An explicit embedding pair always wins; otherwise the
/// run manifest decides between co-embed projection and enrichment.
fn route(args: &AnnotateCliArgs) -> Result<Route> {
    if let (Some(feat), Some(cell)) = (
        args.feature_embedding.as_deref(),
        args.cell_embedding.as_deref(),
    ) {
        if args.method == AnnotateMethod::Enrichment {
            return Ok(Route::Enrichment);
        }
        let prefix = args
            .out
            .as_deref()
            .or(args.from.as_deref())
            .map(run_manifest::derive_out_prefix)
            .unwrap_or_else(|| "annot".into());
        return Ok(Route::EmbeddingFiles {
            feat: feat.to_string(),
            cell: cell.to_string(),
            prefix,
        });
    }
    let projection = match args.method {
        AnnotateMethod::Enrichment => false,
        AnnotateMethod::Projection => true,
        AnnotateMethod::Auto => args
            .from
            .as_deref()
            .and_then(|f| run_manifest::load(f).ok())
            .is_some_and(|l| manifest_prefers_projection(&l.manifest)),
    };
    Ok(if projection {
        Route::Projection
    } else {
        Route::Enrichment
    })
}

/// A run with a co-embedded gene space annotates by projection.
fn manifest_prefers_projection(manifest: &RunManifest) -> bool {
    manifest.outputs.feature_coembedding.is_some()
        || (manifest.kind.coembeds() && manifest.outputs.feature_embedding.is_some())
}

fn run_enrichment(args: &AnnotateCliArgs) -> Result<()> {
    let from = args
        .from
        .clone()
        .context("--from is required for enrichment annotation")?;
    anyhow::ensure!(
        !args.markers.is_empty() || args.gaf.is_some() || args.gmt.is_some(),
        "enrichment needs --markers, --gaf, or --gmt"
    );
    annotate_by_enrichment(&build_enrichment_args(args, from))
}

fn run_projection(args: &AnnotateCliArgs) -> Result<()> {
    let from = args
        .from
        .clone()
        .context("--from is required for projection annotation")?;
    anyhow::ensure!(!args.markers.is_empty(), "projection needs --markers");
    annotate_by_projection(&build_projection_args(args, from))
}

/// Explicit embedding pair: the same projection pass as the manifest route
/// (bootstrap and abstention included), fed from the two parquets directly.
fn run_projection_from_files(
    args: &AnnotateCliArgs,
    feat_path: &str,
    cell_path: &str,
    prefix: &str,
) -> Result<()> {
    anyhow::ensure!(!args.markers.is_empty(), "projection needs --markers");

    let feat = DMatrix::<f32>::from_parquet(feat_path)
        .with_context(|| format!("reading feature embedding {feat_path}"))?;
    let cell = DMatrix::<f32>::from_parquet(cell_path)
        .with_context(|| format!("reading cell embedding {cell_path}"))?;

    by_projection::run(
        &build_projection_args(args, prefix.into()),
        prefix,
        &ProjectionInputs {
            feature_embedding: &feat,
            cell_embedding: &cell,
        },
    )?;
    Ok(())
}

fn build_enrichment_args(args: &AnnotateCliArgs, from: Box<str>) -> AnnotateArgs {
    AnnotateArgs {
        from,
        clusters: args.clusters.clone(),
        knn: args.knn.unwrap_or(15),
        resolution: args.resolution,
        num_clusters: args.num_clusters,
        min_cluster_size: args.min_cluster_size,
        cluster_seed: args.cluster_seed,
        markers: args.markers.clone(),
        gaf: args.gaf.clone(),
        gmt: args.gmt.clone(),
        no_iea: false,
        min_gene_set: 15,
        max_gene_set: 500,
        out: args.out.clone(),
        block_size: 1024,
        num_draws: 1000,
        num_perm: args.num_perm,
        min_markers: args.min_markers,
        fdr_alpha: args.fdr_alpha,
        q_temperature: args.q_temperature,
        min_confidence: 0.0,
        seed: args.seed,
        no_clean: args.no_clean,
        preload_data: false,
        no_empirical_specificity: false,
        obo: args.obo.clone(),
        label_cl: args.label_cl.clone(),
        ontology_fdr_q: args.ontology_fdr_q,
        ontology_by: args.ontology_by,
        no_gene_strata: false,
        no_bootstrap_markers: false,
        n_boot: args.n_boot,
        boot_num_draws: 100,
        min_support: args.min_support,
        abstain_separable: args.abstain_separable,
        abstain_alpha: args.abstain_alpha,
        set_coverage: args.set_coverage,
        max_set_size: args.max_set_size,
    }
}

fn build_projection_args(args: &AnnotateCliArgs, from: Box<str>) -> AnnotateProjectionArgs {
    AnnotateProjectionArgs {
        from,
        markers: args.markers.clone(),
        out: args.out.clone(),
        knn: args.knn.unwrap_or(30),
        resolution: args.resolution,
        num_perm: args.num_perm,
        seed: args.seed,
        min_markers: args.min_markers,
        no_idf: args.no_idf,
        no_assign_qc: args.no_assign_qc,
        assign_mad: args.assign_mad,
        fdr_alpha: args.fdr_alpha,
        q_temperature: args.q_temperature,
        obo: args.obo.clone(),
        label_cl: args.label_cl.clone(),
        ontology_fdr_q: args.ontology_fdr_q,
        ontology_by: args.ontology_by,
        panel_perm: 0,
        support_perm: 0,
        no_bootstrap_markers: false,
        n_boot: args.n_boot,
        no_recluster: false,
        min_support: args.min_support,
        abstain_separable: args.abstain_separable,
        abstain_alpha: args.abstain_alpha,
        set_coverage: args.set_coverage,
        max_set_size: args.max_set_size,
        no_clean: args.no_clean,
    }
}
