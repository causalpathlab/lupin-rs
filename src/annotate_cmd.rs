//! Unified `lupin annotate` — enrichment, or projection from a run manifest or embedding files.

use crate::annotate::args::{AnnotateArgs, AnnotateOntologyArgs, AnnotateProjectionArgs};
use crate::annotate::by_projection::{self, ProjectionInputs};
use crate::manifest::annotate::{
    annotate_by_enrichment, annotate_by_projection, annotate_ontology,
};
use crate::manifest::run::RunManifest;
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

#[derive(Args, Debug, Clone)]
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
        help = "Run manifest (`run.senna.json` or pinto's `run.pinto.json`) or its output prefix"
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
        help = "Marker TSV (`gene<TAB>celltype`); defaults to the round's `annotate.markers`. Omit for ontology-only follow-up on enrichment output"
    )]
    pub markers: Box<str>,

    #[arg(
        long,
        short = 'o',
        default_value = "",
        help = "Output prefix for every file this command writes. Without it, annotate opens its TUI, which asks for it (the run's prefix + `.L1` offered)"
    )]
    pub out: Box<str>,

    // ── shared clustering / stats ──
    // Default is method-specific when omitted: enrichment 15, projection/ORA 30.
    #[arg(long)]
    pub knn: Option<usize>,
    #[arg(long, default_value_t = 1.0)]
    pub resolution: f64,
    #[arg(long, default_value_t = 42)]
    pub seed: u64,
    #[arg(long = "num-perm", default_value_t = 1000)]
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
    #[arg(
        long,
        conflicts_with = "clusters",
        help = "A pinto run's cascade level to annotate (a `levels[].tag`, e.g. L2); default: its final level"
    )]
    pub level: Option<Box<str>>,
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
    #[arg(
        long,
        conflicts_with_all = ["gaf", "gmt"],
        help = "Also score GO terms per cluster, with the GO annotations of the run's species (human or mouse, from its gene names), downloaded once into the cache"
    )]
    pub go: bool,
    #[arg(
        long = "go-obo",
        help = "Gene Ontology .obo for --go/--gaf/--gmt; default: found or downloaded like the Cell Ontology"
    )]
    pub go_obo: Option<Box<str>>,
    #[arg(
        long = "go-min-overlap",
        default_value_t = 20,
        help = "Report GO/GMT terms sharing at least this many genes with the data's features (the overlap, not the term's own size); every term is tested"
    )]
    pub go_min_overlap: usize,
    #[arg(
        long = "go-max-overlap",
        default_value_t = 500,
        help = "Report GO/GMT terms sharing at most this many genes with the data's features; every term is tested"
    )]
    pub go_max_overlap: usize,

    // ── enrichment's SuSiE stage ──
    #[arg(
        long = "no-susie",
        help = "Call clusters by the enrichment's softmax share, without the SuSiE stage \
                (cluster counts regressed on the panel, so types compete for shared markers)"
    )]
    pub no_susie: bool,
    #[arg(
        long = "susie-effects",
        default_value_t = 5,
        help = "SuSiE: single effects per cluster (cell types one cluster can be)"
    )]
    pub susie_effects: usize,
    #[arg(
        long = "susie-dispersion",
        help = "SuSiE: one NB dispersion for every gene (0: Poisson); default: gene-specific, \
                estimated across the clusters"
    )]
    pub susie_dispersion: Option<f32>,
    #[arg(long = "mcmc-samples", default_value_t = 1000)]
    pub mcmc_samples: usize,
    #[arg(long = "mcmc-warmup", default_value_t = 500)]
    pub mcmc_warmup: usize,
    #[arg(long = "mcmc-thin", default_value_t = 1)]
    pub mcmc_thin: usize,

    // ── projection / ORA ──
    #[arg(long = "no-idf")]
    pub no_idf: bool,
    #[arg(long = "no-assign-qc")]
    pub no_assign_qc: bool,
    #[arg(long = "assign-mad", default_value_t = 2.5)]
    pub assign_mad: f64,

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

    #[arg(
        long,
        help = "Call a first round at the panel's fine cell types instead of its coarse groups \
                (later rounds are always fine)"
    )]
    pub fine: bool,

    #[arg(
        long,
        value_enum,
        default_value_t,
        help = "In the TUI, how figures reach the terminal"
    )]
    pub graphics: crate::tui::Graphics,
}

impl AnnotateCliArgs {
    /// Without an output prefix, annotate opens its TUI, which asks for one.
    #[must_use]
    pub fn opens_tui(&self) -> bool {
        self.out.is_empty()
    }

    /// The `annotate` arguments that give back `self`, for re-running a pass
    /// as a child process (with `--out`, so the child runs without the TUI).
    #[must_use]
    pub fn to_argv(&self) -> Vec<String> {
        let mut v: Vec<String> = Vec::new();
        let mut val = |flag: &str, x: String| {
            v.push(format!("--{flag}"));
            v.push(x);
        };
        let method = match self.method {
            AnnotateMethod::Auto => "auto",
            AnnotateMethod::Enrichment => "enrichment",
            AnnotateMethod::Projection => "projection",
        };
        val("method", method.into());
        val("markers", self.markers.to_string());
        val("out", self.out.to_string());
        val("resolution", self.resolution.to_string());
        val("seed", self.seed.to_string());
        val("num-perm", self.num_perm.to_string());
        val("min-markers", self.min_markers.to_string());
        val("fdr-alpha", self.fdr_alpha.to_string());
        val("q-temperature", self.q_temperature.to_string());
        val("min-cluster-size", self.min_cluster_size.to_string());
        val("assign-mad", self.assign_mad.to_string());
        val("ontology-fdr-q", self.ontology_fdr_q.to_string());
        val("susie-effects", self.susie_effects.to_string());
        val("mcmc-samples", self.mcmc_samples.to_string());
        val("mcmc-warmup", self.mcmc_warmup.to_string());
        val("mcmc-thin", self.mcmc_thin.to_string());
        let opts: [(&str, Option<String>); 16] = [
            ("from", self.from.as_deref().map(String::from)),
            (
                "feature-embedding",
                self.feature_embedding.as_deref().map(String::from),
            ),
            (
                "cell-embedding",
                self.cell_embedding.as_deref().map(String::from),
            ),
            ("knn", self.knn.map(|x| x.to_string())),
            // `--level` resolves to its clusters again in the child.
            (
                "clusters",
                self.clusters
                    .as_deref()
                    .filter(|_| self.level.is_none())
                    .map(String::from),
            ),
            ("level", self.level.as_deref().map(String::from)),
            ("num-clusters", self.num_clusters.map(|x| x.to_string())),
            ("cluster-seed", self.cluster_seed.map(|x| x.to_string())),
            ("gaf", self.gaf.as_deref().map(String::from)),
            ("gmt", self.gmt.as_deref().map(String::from)),
            ("obo", self.obo.as_deref().map(String::from)),
            ("go-obo", self.go_obo.as_deref().map(String::from)),
            ("go-min-overlap", Some(self.go_min_overlap.to_string())),
            ("go-max-overlap", Some(self.go_max_overlap.to_string())),
            ("label-cl", self.label_cl.as_deref().map(String::from)),
            (
                "susie-dispersion",
                self.susie_dispersion.map(|x| x.to_string()),
            ),
        ];
        for (flag, x) in opts {
            if let Some(x) = x {
                val(flag, x);
            }
        }
        let flags = [
            ("no-clean", self.no_clean),
            ("no-idf", self.no_idf),
            ("no-assign-qc", self.no_assign_qc),
            ("ontology-by", self.ontology_by),
            ("use-perm-p", self.use_perm_p),
            ("fine", self.fine),
            ("go", self.go),
            ("no-susie", self.no_susie),
        ];
        v.extend(
            flags
                .iter()
                .filter(|(_, on)| *on)
                .map(|(f, _)| format!("--{f}")),
        );
        v
    }
}

pub fn run_annotate(args: &AnnotateCliArgs) -> Result<()> {
    if args.opens_tui() {
        return crate::tui::run(args, None);
    }
    // One manifest load per invocation; every route below reuses it.
    let loaded = args
        .from
        .as_deref()
        .map(crate::manifest::run::load)
        .transpose()?;
    let mut args = args.clone();

    // A round carries its marker panel (revised by `relabel`); use it when
    // no other source of labels is given.
    if let Some(markers) = loaded.as_ref().and_then(|l| round_markers(&args, l)) {
        log::info!("No -m given: using the round's marker panel {markers}");
        args.markers = markers.into_boxed_str();
    }

    if let Some(l) = &loaded {
        apply_level(&mut args, &l.file)?;
    }
    // Ask before anything under -o is erased, not when the manifest is saved.
    if let Some(l) = &loaded {
        crate::manifest::run::may_replace(&crate::manifest::run::annotated_path(
            &l.file, &args.out,
        ))?;
    }

    if is_ontology_followup(&args) {
        let loaded = loaded
            .as_ref()
            .context("--from required for ontology follow-up")?;
        return annotate_ontology(&build_ontology_args(&args)?, loaded);
    }

    // A marker pass on a run: group the panel's cell types (finding the Cell
    // Ontology, whose walk then runs by default), and call a first round at
    // the coarse level.
    // The run's Cell Ontology data, read once for the grouping and the pass.
    let cl_data = match &loaded {
        Some(l) if !args.markers.is_empty() => Some(crate::manifest::ontology::load(
            Some(&l.dir),
            &args.markers,
            args.obo.as_deref(),
            args.label_cl.as_deref(),
            crate::manifest::data_files::Fetch::Allowed,
        )?),
        _ => None,
    };
    let prepared = cl_data
        .as_ref()
        .map(|d| {
            crate::manifest::first_round::prepare(
                &args.markers,
                &args.out,
                d,
                args.label_cl.as_deref(),
            )
        })
        .transpose()?;
    if let Some((obo, label_cl)) = prepared.as_ref().and_then(|p| p.ontology.clone()) {
        args.obo = Some(obo.into_boxed_str());
        args.label_cl = Some(label_cl.into_boxed_str());
    }
    resolve_go(&mut args, loaded.as_ref())?;
    let args = &args;

    match route(args, loaded.as_ref()) {
        Route::EmbeddingFiles { feat, cell } => run_projection_from_files(args, feat, cell)?,
        Route::Enrichment => {
            anyhow::ensure!(
                !args.markers.is_empty() || args.go || args.gaf.is_some() || args.gmt.is_some(),
                "enrichment needs --markers, --go, --gaf, or --gmt"
            );
            let loaded = loaded
                .as_ref()
                .context("--from is required for enrichment annotation")?;
            annotate_by_enrichment(&build_enrichment_args(args), loaded, cl_data.as_ref())?;
        }
        Route::Projection => {
            anyhow::ensure!(!args.markers.is_empty(), "projection needs --markers");
            anyhow::ensure!(
                !args.go && args.gaf.is_none() && args.gmt.is_none(),
                "GO terms are scored by the enrichment pass; add --method enrichment"
            );
            let loaded = loaded
                .as_ref()
                .context("--from is required for projection annotation")?;
            annotate_by_projection(&build_projection_args(args), loaded)?;
        }
    }
    if let (Some(p), Some(l)) = (&prepared, &loaded) {
        let coarse = !args.fine && crate::manifest::first_round::is_first_round(l);
        let manifest = crate::manifest::run::annotated_path(&l.file, &args.out);
        crate::manifest::first_round::finish(&manifest, &p.tree, coarse)?;
    }
    Ok(())
}

/// `--level`: annotate that level of a pinto run by taking its propensity as
/// `--clusters`; the tag stays in `args` for the round's settings. A no-op
/// without `--level`.
pub fn apply_level(args: &mut AnnotateCliArgs, manifest: &std::path::Path) -> Result<()> {
    let Some(tag) = args.level.clone() else {
        return Ok(());
    };
    anyhow::ensure!(
        crate::manifest::pinto::is_pinto(manifest),
        "--level is for a pinto run (.pinto.json); {} is not one",
        manifest.display()
    );
    let path = crate::manifest::pinto::level_propensity(manifest, &tag)?;
    log::info!("level {tag}: clusters from {path}");
    args.clusters = Some(path.into_boxed_str());
    Ok(())
}

/// The round's `annotate.markers`, resolved, when `-m`, `--gaf` and `--gmt`
/// are all absent and this is not an ontology follow-up.
fn round_markers(args: &AnnotateCliArgs, loaded: &crate::manifest::run::Loaded) -> Option<String> {
    if !args.markers.is_empty()
        || args.gaf.is_some()
        || args.gmt.is_some()
        || is_ontology_followup(args)
    {
        return None;
    }
    recorded_markers(loaded)
}

/// The marker panel `loaded` recorded, found again if the run moved (a path
/// from another machine is looked up by its tail); `None` when it is gone.
pub(crate) fn recorded_markers(loaded: &crate::manifest::run::Loaded) -> Option<String> {
    loaded.recorded(loaded.manifest.annotate.markers.as_deref()?)
}

/// Settle a GO pass's ontology: `--go-obo`, else (in a pass without markers,
/// where `--obo` has no Cell Ontology to name) `--obo`, else the Gene
/// Ontology found or downloaded like the Cell Ontology. `--go`'s annotations
/// are settled by the pass, which reads the data's gene names.
fn resolve_go(
    args: &mut AnnotateCliArgs,
    loaded: Option<&crate::manifest::run::Loaded>,
) -> Result<()> {
    let dir = loaded.map(|l| l.dir.as_path());
    if (args.go || args.gaf.is_some() || args.gmt.is_some()) && args.go_obo.is_none() {
        args.go_obo = match args.obo.take() {
            Some(obo) if args.markers.is_empty() => Some(obo),
            cl => {
                args.obo = cl;
                Some(
                    crate::manifest::data_files::go_ontology(dir)?
                        .to_string_lossy()
                        .into(),
                )
            }
        };
    }
    Ok(())
}

fn is_ontology_followup(args: &AnnotateCliArgs) -> bool {
    args.markers.is_empty()
        && args.gaf.is_none()
        && args.gmt.is_none()
        && args.from.is_some()
        && args.obo.is_some()
        && args.label_cl.is_some()
}

fn build_ontology_args(args: &AnnotateCliArgs) -> Result<AnnotateOntologyArgs> {
    Ok(AnnotateOntologyArgs {
        label_cl: args.label_cl.clone().context("--label-cl required")?,
        obo: args.obo.clone().context("--obo required")?,
        out: args.out.clone(),
        fdr_q: args.ontology_fdr_q,
        by: args.ontology_by,
        use_perm_p: args.use_perm_p,
    })
}

enum Route<'a> {
    Enrichment,
    /// Co-embed projection through the run manifest.
    Projection,
    /// Projection over an explicit `--feature-embedding` / `--cell-embedding` pair.
    EmbeddingFiles {
        feat: &'a str,
        cell: &'a str,
    },
}

/// Pick the backend. An explicit embedding pair wins; otherwise the run
/// manifest decides between co-embed projection and enrichment.
fn route<'a>(
    args: &'a AnnotateCliArgs,
    loaded: Option<&crate::manifest::run::Loaded>,
) -> Route<'a> {
    if let (Some(feat), Some(cell)) = (
        args.feature_embedding.as_deref(),
        args.cell_embedding.as_deref(),
    ) {
        if args.method != AnnotateMethod::Enrichment {
            return Route::EmbeddingFiles { feat, cell };
        }
    }
    let projection = match args.method {
        AnnotateMethod::Enrichment => false,
        AnnotateMethod::Projection => true,
        AnnotateMethod::Auto => loaded.is_some_and(|l| manifest_prefers_projection(&l.manifest)),
    };
    if projection {
        Route::Projection
    } else {
        Route::Enrichment
    }
}

/// A run with a co-embedded gene space annotates by projection. A kind that
/// co-embeds needs that space itself: its raw gene embedding is not on the
/// cell manifold, so a run written before the co-embedding existed falls back
/// to enrichment.
fn manifest_prefers_projection(manifest: &RunManifest) -> bool {
    manifest.outputs.feature_coembedding.is_some()
        || (!manifest.kind.coembeds() && manifest.outputs.feature_embedding.is_some())
}

/// Explicit embedding pair: the same projection pass as the manifest route
/// fed from the two parquets directly.
/// No manifest is touched.
fn run_projection_from_files(
    args: &AnnotateCliArgs,
    feat_path: &str,
    cell_path: &str,
) -> Result<()> {
    anyhow::ensure!(!args.markers.is_empty(), "projection needs --markers");
    let feat = DMatrix::<f32>::from_parquet(feat_path)
        .with_context(|| format!("reading feature embedding {feat_path}"))?;
    let cell = DMatrix::<f32>::from_parquet(cell_path)
        .with_context(|| format!("reading cell embedding {cell_path}"))?;
    by_projection::run(
        &build_projection_args(args),
        &ProjectionInputs {
            feature_embedding: &feat,
            cell_embedding: &cell,
        },
    )?;
    Ok(())
}

/// The enrichment settings `lupin annotate` uses by default, writing under
/// `out`: parsed from the command's own defaults, so the two cannot drift.
pub(crate) fn default_enrichment_args(out: &str) -> AnnotateArgs {
    #[derive(clap::Parser)]
    struct Defaults {
        #[command(flatten)]
        annotate: AnnotateCliArgs,
    }
    let d = <Defaults as clap::Parser>::parse_from(["lupin", "-o", out]);
    build_enrichment_args(&d.annotate)
}

pub(crate) fn build_enrichment_args(args: &AnnotateCliArgs) -> AnnotateArgs {
    AnnotateArgs {
        clusters: args.clusters.clone(),
        level: args.level.clone(),
        knn: args.knn.unwrap_or(15),
        resolution: args.resolution,
        num_clusters: args.num_clusters,
        min_cluster_size: args.min_cluster_size,
        cluster_seed: args.cluster_seed,
        markers: args.markers.clone(),
        gaf: args.gaf.clone(),
        gmt: args.gmt.clone(),
        go: args.go,
        go_obo: args.go_obo.clone(),
        go_min_overlap: args.go_min_overlap,
        go_max_overlap: args.go_max_overlap,
        out: args.out.clone(),
        num_perm: args.num_perm,
        min_markers: args.min_markers,
        fdr_alpha: args.fdr_alpha,
        q_temperature: args.q_temperature,
        seed: args.seed,
        no_clean: args.no_clean,
        obo: args.obo.clone(),
        label_cl: args.label_cl.clone(),
        ontology_fdr_q: args.ontology_fdr_q,
        ontology_by: args.ontology_by,
        susie: (!args.no_susie).then(|| build_susie_config(args)),
    }
}

fn build_susie_config(args: &AnnotateCliArgs) -> crate::annotate::susie::SusieConfig {
    let mut cfg = crate::annotate::susie::SusieConfig {
        samples: args.mcmc_samples,
        warmup: args.mcmc_warmup,
        thin: args.mcmc_thin,
        seed: args.seed,
        dispersion: args.susie_dispersion,
        ..Default::default()
    };
    cfg.prior.num_effects = args.susie_effects;
    cfg
}

fn build_projection_args(args: &AnnotateCliArgs) -> AnnotateProjectionArgs {
    AnnotateProjectionArgs {
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
        no_clean: args.no_clean,
    }
}

#[cfg(test)]
mod route_tests {
    use super::*;
    use crate::manifest::run::RunKind;

    #[test]
    fn a_co_embedding_kind_without_its_co_embedding_falls_back_to_enrichment() {
        let mut m = RunManifest::new(RunKind::Bge, "run");
        m.outputs.feature_embedding = Some("run.feature_embedding.parquet".into());
        assert!(!manifest_prefers_projection(&m));
        m.outputs.feature_coembedding = Some("run.feature_coembedding.parquet".into());
        assert!(manifest_prefers_projection(&m));
    }
}

#[cfg(test)]
#[path = "tests/annotate_cmd.rs"]
mod argv_tests;
