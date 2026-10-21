//! Unified `lupin annotate` — enrichment, or projection from a run manifest or embedding files.

use crate::annotate::args::{
    AnnotateArgs, AnnotateOntologyArgs, MarkerScoreArgs, ProjectionMethod,
};
use crate::annotate::susie::SusieConfig;
use crate::manifest::annotate::{
    annotate_by_enrichment, annotate_by_marker_score, annotate_ontology, LeidenArgs,
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
        help = "Annotation backend: enrichment (raw counts), projection (the run's embeddings; see --projection), or auto (every method the run supports, their readings in one round: the enrichment with SuSiE, the projection score where the run decodes cells, the gene-set units where its genes are placed among the cells; see --call)"
    )]
    pub method: AnnotateMethod,

    #[arg(
        long,
        value_enum,
        help = "Which method's reading makes the round's call (its labels) when several run; \
                default: SuSiE's where it ran, else the enrichment's. A method that does not \
                run on this pass is refused"
    )]
    pub call: Option<crate::annotate::rounds::Method>,

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
    // Default is method-specific when omitted: enrichment 15, projection 30.
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
        default_value_t = crate::annotate::args::GO_MIN_OVERLAP,
        help = "GO/GMT terms sharing at least this many genes with the data's features (the overlap, not the term's own size): the enrichment tests every term and reports these; the units train only these"
    )]
    pub go_min_overlap: usize,
    #[arg(
        long = "go-max-overlap",
        default_value_t = crate::annotate::args::GO_MAX_OVERLAP,
        help = "GO/GMT terms sharing at most this many genes with the data's features: the enrichment tests every term and reports these; the units train only these"
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
        default_value_t = SusieConfig::default().prior.num_effects,
        help = "SuSiE: single effects per cluster (cell types one cluster can be)"
    )]
    pub susie_effects: usize,
    /// Retired with the NB likelihood (0.5.0 shipped it); accepted so saved
    /// commands keep working, and ignored with a warning.
    #[arg(long = "susie-dispersion", hide = true)]
    pub susie_dispersion: Option<f32>,
    #[arg(
        long = "mcmc-chains",
        default_value_t = SusieConfig::default().chains,
        help = "SuSiE: independent chains per cluster, pooled (R̂ compares them)"
    )]
    pub mcmc_chains: usize,
    #[arg(
        long = "mcmc-samples",
        default_value_t = SusieConfig::default().samples,
        help = "SuSiE: samples per chain"
    )]
    pub mcmc_samples: usize,
    #[arg(long = "mcmc-warmup", default_value_t = SusieConfig::default().warmup)]
    pub mcmc_warmup: usize,
    #[arg(long = "mcmc-thin", default_value_t = SusieConfig::default().thin)]
    pub mcmc_thin: usize,

    // ── projection ──
    #[arg(
        long,
        value_enum,
        help = "Projection: score (the embedding marker score on the run's decoder gene table, \
                per cell with its margin) or units (cell types and GO terms trained as gene sets on \
                the run's fixed genes); default: the score where the run decodes cells, else the units"
    )]
    pub projection: Option<ProjectionMethod>,
    #[arg(
        long,
        value_enum,
        default_value_t = crate::gene_text::cli::ComputeDevice::Cpu,
        help = "Device training the gene-set units (--projection units, and beside the enrichment under auto)"
    )]
    pub device: crate::gene_text::cli::ComputeDevice,
    #[arg(long, default_value_t = 0, help = "Device ordinal (for cuda/metal)")]
    pub device_no: usize,
    /// Retired with the term-ORA projection; accepted so saved commands
    /// keep working, and ignored with a warning.
    #[arg(long = "no-idf", hide = true)]
    pub no_idf: bool,
    #[arg(long = "no-assign-qc", hide = true)]
    pub no_assign_qc: bool,
    #[arg(long = "assign-mad", hide = true)]
    pub assign_mad: Option<f64>,

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
    /// Gene sets asked for: `--go`, `--gaf` or `--gmt`.
    #[must_use]
    pub fn wants_gene_sets(&self) -> bool {
        self.go || self.gaf.is_some() || self.gmt.is_some()
    }

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
        if let Some(p) = self.projection {
            val("projection", p.as_str().into());
        }
        if let Some(c) = self.call {
            val("call", c.name().into());
        }
        val("markers", self.markers.to_string());
        val("out", self.out.to_string());
        val("resolution", self.resolution.to_string());
        val("seed", self.seed.to_string());
        val("num-perm", self.num_perm.to_string());
        val("min-markers", self.min_markers.to_string());
        val("fdr-alpha", self.fdr_alpha.to_string());
        val("q-temperature", self.q_temperature.to_string());
        val("min-cluster-size", self.min_cluster_size.to_string());
        val("ontology-fdr-q", self.ontology_fdr_q.to_string());
        val("susie-effects", self.susie_effects.to_string());
        val("mcmc-chains", self.mcmc_chains.to_string());
        val("mcmc-samples", self.mcmc_samples.to_string());
        val("mcmc-warmup", self.mcmc_warmup.to_string());
        val("mcmc-thin", self.mcmc_thin.to_string());
        val("device", self.device.as_str().into());
        val("device-no", self.device_no.to_string());
        let opts: [(&str, Option<String>); 15] = [
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
        ];
        for (flag, x) in opts {
            if let Some(x) = x {
                val(flag, x);
            }
        }
        let flags = [
            ("no-clean", self.no_clean),
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

/// `lupin annotate`: a pass with `-o`, else the TUI, in a session as
/// `session` asks.
pub fn run_annotate(args: &AnnotateCliArgs, session: &legume_tui::cli::SessionArgs) -> Result<()> {
    if args.opens_tui() {
        return crate::tui::run(args, None, session, false);
    }
    if let Some(cfg) = build_enrichment_args(args).susie {
        cfg.validate()?;
    }
    if args.no_idf || args.no_assign_qc || args.assign_mad.is_some() {
        log::warn!(
            "--no-idf, --no-assign-qc and --assign-mad are retired with the term-ORA \
             projection; ignoring them"
        );
    }
    // Refused before any prompt or download, as SuSiE's settings are.
    match args.projection {
        Some(ProjectionMethod::Score) => drop(build_score_args(args)?),
        Some(ProjectionMethod::Units) => drop(build_units_args(args)?),
        None => {}
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

    check_call(&args, loaded.as_ref())?;
    let Some(l) = &loaded else {
        return run_pass(args, None);
    };
    apply_level(&mut args, &l.file)?;
    // Ask before anything under -o is erased, not when the manifest is saved.
    crate::manifest::run::may_replace(&crate::manifest::run::annotated_path(&l.file, &args.out))?;
    // A pass on a run writes under a staging prefix and takes the round's
    // place only once it is done: stopped or failed, it leaves the old one.
    let out = args.out.to_string();
    let staging = crate::manifest::staging::staging_prefix(&out);
    crate::manifest::staging::discard(&staging)?;
    args.out = staging.as_str().into();
    match run_pass(args, Some(l)) {
        Ok(()) => crate::manifest::staging::promote(&l.file, &staging, &out),
        Err(e) => {
            if let Err(d) = crate::manifest::staging::discard(&staging) {
                log::warn!("could not remove the failed pass's files: {d:#}");
            }
            Err(e)
        }
    }
}

/// The pass `args` asks for, on `loaded` when it is a run's.
fn run_pass(
    mut args: AnnotateCliArgs,
    loaded: Option<&crate::manifest::run::Loaded>,
) -> Result<()> {
    if is_ontology_followup(&args) {
        let loaded = loaded.context("--from required for ontology follow-up")?;
        return annotate_ontology(&build_ontology_args(&args)?, loaded);
    }

    // A marker pass on a run: group the panel's cell types (finding the Cell
    // Ontology, whose walk then runs by default), and call a first round at
    // the coarse level.
    // The run's Cell Ontology data, read once for the grouping and the pass.
    let cl_data = match loaded {
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
    resolve_go(&mut args, loaded)?;
    let args = &args;

    match route(args, loaded) {
        Route::EmbeddingFiles { feat, cell } => run_projection_from_files(args, feat, cell)?,
        Route::Enrichment => {
            anyhow::ensure!(
                !args.markers.is_empty() || args.wants_gene_sets(),
                "enrichment needs --markers, --go, --gaf, or --gmt"
            );
            let loaded = loaded.context("--from is required for enrichment annotation")?;
            annotate_by_enrichment(&build_enrichment_args(args), loaded, cl_data.as_ref())?;
        }
        Route::Projection => {
            anyhow::ensure!(!args.markers.is_empty(), "projection needs --markers");
            let loaded = loaded.context("--from is required for projection annotation")?;
            let projection = resolve_projection(args, Some(&loaded.manifest));
            check_projection_gene_sets(args, projection)?;
            match projection {
                ProjectionMethod::Score => {
                    annotate_by_marker_score(&build_score_args(args)?, loaded)?
                }
                ProjectionMethod::Units => {
                    crate::manifest::annotate::annotate_by_units(&build_units_args(args)?, loaded)?
                }
            }
        }
    }
    if let (Some(p), Some(l)) = (&prepared, loaded) {
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
    if args.wants_gene_sets() && args.go_obo.is_none() {
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
/// decides between projection and enrichment ([`resolve_method_on`]).
fn route<'a>(
    args: &'a AnnotateCliArgs,
    loaded: Option<&crate::manifest::run::Loaded>,
) -> Route<'a> {
    if let (Some(feat), Some(cell)) = (
        args.feature_embedding.as_deref(),
        args.cell_embedding.as_deref(),
    ) {
        // Auto with GO terms or no panel is the enrichment's, files or not.
        let auto_enrichment = args.method == AnnotateMethod::Auto && enrichment_only(args);
        if args.method != AnnotateMethod::Enrichment && !auto_enrichment {
            return Route::EmbeddingFiles { feat, cell };
        }
    }
    if resolve_method_on(args, loaded) == AnnotateMethod::Projection {
        Route::Projection
    } else {
        Route::Enrichment
    }
}

/// The method a pass runs: the one chosen; for `auto`, the enrichment, with
/// every other method the run supports beside it (the projection score, the
/// units), unless `--projection` asks for that projection's own call on a
/// run with its gene table (GO terms or no panel staying the enrichment's,
/// unless the units are asked to align them).
#[must_use]
pub fn resolve_method(args: &AnnotateCliArgs, manifest: Option<&RunManifest>) -> AnnotateMethod {
    let projects =
        |m: &RunManifest| manifest_prefers_projection(m, resolve_projection(args, Some(m)));
    match args.method {
        AnnotateMethod::Auto
            if args.projection.is_some()
                && !enrichment_only(args)
                && manifest.is_some_and(projects) =>
        {
            AnnotateMethod::Projection
        }
        AnnotateMethod::Auto => AnnotateMethod::Enrichment,
        chosen => chosen,
    }
}

/// [`resolve_method`] on a run as it is on disk: `auto` on a run whose
/// enrichment has no expression to read (no raw counts here, no earlier
/// pass's cache, no decoder to stand in) projects instead, as the run's gene
/// table allows: the score where it decodes cells, else the units (which
/// align the GO terms asked for; the score cannot take them).
#[must_use]
pub fn resolve_method_on(
    args: &AnnotateCliArgs,
    loaded: Option<&crate::manifest::run::Loaded>,
) -> AnnotateMethod {
    let method = resolve_method(args, loaded.map(|l| &l.manifest));
    let Some(l) = loaded else { return method };
    let projection = resolve_projection(args, Some(&l.manifest));
    let projects = args.method == AnnotateMethod::Auto
        && method == AnnotateMethod::Enrichment
        && !args.markers.is_empty()
        && (projection == ProjectionMethod::Units || !args.wants_gene_sets())
        && manifest_prefers_projection(&l.manifest, projection)
        && !enrichment_reads(l);
    if projects {
        log::info!(
            "no expression for the enrichment (no raw counts, cache or decoder): \
             the {} alone",
            projection.as_str()
        );
        AnnotateMethod::Projection
    } else {
        method
    }
}

/// Whether the enrichment has expression to read on `l`: its raw counts, or
/// what stands in for them ([`crate::manifest::no_counts`]).
fn enrichment_reads(l: &crate::manifest::run::Loaded) -> bool {
    use crate::manifest::no_counts::{source, Source};
    !l.manifest.data.input.is_empty() && !matches!(source(l), Source::Nothing)
}

/// The methods whose readings a pass of `args` on `loaded` records, in
/// [`Method::ALL`]'s order: a projection pass its own; an enrichment pass
/// the enrichment, SuSiE (on count sums, the raw counts or an earlier
/// pass's cache of them, unless `--no-susie`) and, under
/// `auto`, each projection the run has the gene table for.
///
/// [`Method::ALL`]: crate::annotate::rounds::Method::ALL
#[must_use]
pub fn methods_run(
    args: &AnnotateCliArgs,
    loaded: Option<&crate::manifest::run::Loaded>,
) -> Vec<crate::annotate::rounds::Method> {
    use crate::annotate::rounds::Method;
    let manifest = loaded.map(|l| &l.manifest);
    let of = |p: ProjectionMethod| match p {
        ProjectionMethod::Score => Method::Score,
        ProjectionMethod::Units => Method::Units,
    };
    if let Route::EmbeddingFiles { .. } | Route::Projection = route(args, loaded) {
        return vec![of(resolve_projection(args, manifest))];
    }
    if args.markers.is_empty() {
        return vec![Method::Enrichment];
    }
    // Count sums: the raw counts, or an earlier pass's cache of them (a
    // decoder's expectation is none).
    let counts = loaded.is_some_and(|l| {
        use crate::manifest::no_counts::{source, Source};
        matches!(source(l), Source::Counts | Source::Cache(_))
    });
    let beside = |p: ProjectionMethod| {
        args.method == AnnotateMethod::Auto
            && manifest.is_some_and(|m| manifest_prefers_projection(m, p))
    };
    Method::ALL
        .into_iter()
        .filter(|m| match m {
            Method::Enrichment => true,
            Method::Susie => !args.no_susie && counts,
            Method::Score => beside(ProjectionMethod::Score),
            Method::Units => beside(ProjectionMethod::Units),
        })
        .collect()
}

/// Refuse a `--call` naming a method that does not run on this pass
/// ([`methods_run`]), before anything is read or written.
pub fn check_call(
    args: &AnnotateCliArgs,
    loaded: Option<&crate::manifest::run::Loaded>,
) -> Result<()> {
    let Some(call) = args.call else {
        return Ok(());
    };
    let run = methods_run(args, loaded);
    let names: Vec<&str> = run.iter().map(|m| m.name()).collect();
    anyhow::ensure!(
        run.contains(&call),
        "--call {}: {} does not run on this pass (it runs {}); drop --call, or name one of those",
        call.name(),
        call.name(),
        names.join(", ")
    );
    Ok(())
}

/// Gene sets (`--go`, `--gaf`, `--gmt`) on a projection pass: trained as
/// units by `--projection units`, scored by the enrichment pass otherwise.
fn check_projection_gene_sets(args: &AnnotateCliArgs, projection: ProjectionMethod) -> Result<()> {
    anyhow::ensure!(
        projection == ProjectionMethod::Units || !args.wants_gene_sets(),
        "GO terms are scored by the enrichment pass, or trained as units by \
         --projection units; add --method enrichment or --projection units"
    );
    Ok(())
}

/// No panel, or GO terms the units are not asked to align: what only the
/// enrichment does.
fn enrichment_only(args: &AnnotateCliArgs) -> bool {
    args.markers.is_empty()
        || (args.wants_gene_sets() && args.projection != Some(ProjectionMethod::Units))
}

/// How a projection pass calls: the `--projection` given, else the score
/// where the run decodes cells through its gene table and the units
/// otherwise (also without a run: embedding files do not say which table
/// they are).
#[must_use]
pub fn resolve_projection(
    args: &AnnotateCliArgs,
    manifest: Option<&RunManifest>,
) -> ProjectionMethod {
    args.projection.unwrap_or_else(|| match manifest {
        Some(m) if manifest_prefers_projection(m, ProjectionMethod::Score) => {
            ProjectionMethod::Score
        }
        _ => ProjectionMethod::Units,
    })
}

/// A run annotates by projection when it has the gene table its projection
/// method reads ([`crate::marker_embedding::projection_gene_table`]).
fn manifest_prefers_projection(manifest: &RunManifest, projection: ProjectionMethod) -> bool {
    crate::marker_embedding::projection_gene_table(manifest, projection).is_ok()
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
    // Without a run there are no clusters but `--clusters`: else one group,
    // centred at the mean cell.
    let groups = match args.clusters.as_deref() {
        Some(p) => crate::annotate::inputs::load_cluster_labels(p, &cell.rows)?.0,
        None => vec![0; cell.rows.len()],
    };
    match resolve_projection(args, None) {
        ProjectionMethod::Score => {
            crate::annotate::marker_score::run(&build_score_args(args)?, &feat, &cell, &groups)?;
        }
        ProjectionMethod::Units => {
            crate::annotate::geneset_units::run(&build_units_args(args)?, &feat, &cell, &groups)?;
        }
    }
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

/// Neighbours per cell when `--knn` is not given: the enrichment pass's
/// clustering graph, and projection's.
const ENRICHMENT_KNN: usize = 15;
const PROJECTION_KNN: usize = 30;

/// The `--knn` a pass by `method` uses when none is given.
pub(crate) fn default_knn(method: AnnotateMethod) -> usize {
    if method == AnnotateMethod::Projection {
        PROJECTION_KNN
    } else {
        ENRICHMENT_KNN
    }
}

pub(crate) fn build_enrichment_args(args: &AnnotateCliArgs) -> AnnotateArgs {
    AnnotateArgs {
        clusters: args.clusters.clone(),
        level: args.level.clone(),
        knn: args.knn.unwrap_or(ENRICHMENT_KNN),
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
        score_beside: args.method == AnnotateMethod::Auto,
        units_beside: args.method == AnnotateMethod::Auto,
        device: args.device,
        device_no: args.device_no,
        call: args.call,
    }
}

fn build_susie_config(args: &AnnotateCliArgs) -> crate::annotate::susie::SusieConfig {
    if args.susie_dispersion.is_some() {
        log::warn!(
            "--susie-dispersion is retired: the SuSiE stage no longer uses an NB dispersion \
             (it now estimates each marker's log-fold misfit); ignoring it"
        );
    }
    let mut cfg = crate::annotate::susie::SusieConfig {
        chains: args.mcmc_chains,
        min_markers: args.min_markers,
        samples: args.mcmc_samples,
        warmup: args.mcmc_warmup,
        thin: args.mcmc_thin,
        seed: args.seed,
        ..Default::default()
    };
    cfg.prior.num_effects = args.susie_effects;
    cfg
}

/// The score's arguments; a `--min-markers` other than the score's own is
/// refused, as the score cannot honour it.
fn build_score_args(args: &AnnotateCliArgs) -> Result<MarkerScoreArgs> {
    use graph_embedding_util::type_annotation::embedding_score::MIN_MARKERS;
    anyhow::ensure!(
        args.min_markers == MIN_MARKERS,
        "the marker score keeps types with at least {MIN_MARKERS} live markers; drop --min-markers"
    );
    Ok(MarkerScoreArgs {
        markers: args.markers.clone(),
        out: args.out.clone(),
        clusters: args.clusters.clone(),
        leiden: projection_leiden(args),
        no_clean: args.no_clean,
    })
}

/// Internal Leiden for a projection pass, from the command line's settings.
fn projection_leiden(args: &AnnotateCliArgs) -> LeidenArgs {
    LeidenArgs {
        knn: args.knn.unwrap_or(PROJECTION_KNN),
        resolution: args.resolution,
        num_clusters: args.num_clusters,
        min_cluster_size: args.min_cluster_size,
        seed: args.cluster_seed,
    }
}

/// The units' arguments; a `--min-markers` other than the default is
/// refused, as the units keep every type with two markers in the table.
fn build_units_args(args: &AnnotateCliArgs) -> Result<crate::annotate::args::UnitsArgs> {
    use graph_embedding_util::type_annotation::embedding_score::MIN_MARKERS;
    anyhow::ensure!(
        args.min_markers == MIN_MARKERS,
        "the units train every type with at least two markers in the gene table; drop --min-markers"
    );
    Ok(crate::annotate::args::UnitsArgs {
        clusters: args.clusters.clone(),
        leiden: projection_leiden(args),
        seed: args.seed,
        device: args.device,
        device_no: args.device_no,
        go: args.go,
        gaf: args.gaf.clone(),
        gmt: args.gmt.clone(),
        go_obo: args.go_obo.clone(),
        go_min: args.go_min_overlap,
        go_max: args.go_max_overlap,
        no_clean: args.no_clean,
        ..crate::annotate::args::UnitsArgs::new(&args.markers, &args.out)
    })
}

#[cfg(test)]
pub(crate) mod route_tests {
    use super::*;
    use crate::manifest::run::RunKind;

    pub(crate) fn cli(more: &[&str]) -> AnnotateCliArgs {
        #[derive(clap::Parser)]
        struct Cli {
            #[command(flatten)]
            annotate: AnnotateCliArgs,
        }
        let argv = ["lupin", "-o", "x", "-m", "m.tsv"].iter().chain(more);
        <Cli as clap::Parser>::parse_from(argv).annotate
    }

    fn run_with(kind: RunKind, rho: bool, coembedding: bool) -> RunManifest {
        let mut m = RunManifest::new(kind, "run");
        m.outputs.feature_embedding = rho.then(|| "run.feature_embedding.parquet".into());
        m.outputs.feature_coembedding =
            coembedding.then(|| "run.feature_coembedding.parquet".into());
        m
    }

    #[test]
    fn auto_resolves_to_the_method_the_run_and_the_request_allow() {
        let bge = run_with(RunKind::Bge, true, false);
        let svd = run_with(RunKind::Svd, false, false);
        assert_eq!(
            resolve_method(&cli(&[]), Some(&bge)),
            AnnotateMethod::Enrichment,
            "a decoder run's call is the enrichment's, the score beside it"
        );
        assert_eq!(
            resolve_method(&cli(&["--projection", "score"]), Some(&bge)),
            AnnotateMethod::Projection,
            "the score's call, asked for"
        );
        assert_eq!(
            resolve_method(&cli(&[]), Some(&svd)),
            AnnotateMethod::Enrichment
        );
        assert_eq!(
            resolve_method(&cli(&[]), None),
            AnnotateMethod::Enrichment,
            "no run"
        );
        assert_eq!(
            resolve_method(&cli(&["--go"]), Some(&bge)),
            AnnotateMethod::Enrichment,
            "GO terms are the enrichment's"
        );
        let mut no_panel = cli(&[]);
        no_panel.markers = "".into();
        assert_eq!(
            resolve_method(&no_panel, Some(&bge)),
            AnnotateMethod::Enrichment
        );
        assert_eq!(
            resolve_method(&cli(&["--projection", "units"]), Some(&bge)),
            AnnotateMethod::Enrichment,
            "the units asked for need the co-embedding"
        );
        assert_eq!(
            resolve_method(&cli(&["--method", "enrichment"]), Some(&bge)),
            AnnotateMethod::Enrichment,
            "a choice stands"
        );
    }

    #[test]
    fn gene_sets_go_to_the_units_when_they_are_asked_for() {
        let fne = run_with(RunKind::Fne, true, false);
        let units_go = cli(&["--projection", "units", "--go"]);
        assert_eq!(
            resolve_method(&units_go, Some(&fne)),
            AnnotateMethod::Projection,
            "the units align GO terms themselves"
        );
        assert!(
            matches!(route(&units_go, None), Route::Enrichment),
            "no run"
        );
        assert_eq!(
            resolve_method(&cli(&["--go"]), Some(&fne)),
            AnnotateMethod::Enrichment,
            "unasked, GO terms stay the enrichment's"
        );
    }

    #[test]
    fn the_command_lines_marker_floor_is_the_scores() {
        use graph_embedding_util::type_annotation::embedding_score::MIN_MARKERS;
        assert_eq!(
            cli(&[]).min_markers,
            MIN_MARKERS,
            "a default pass must not be refused"
        );
    }

    #[test]
    fn auto_with_go_terms_does_not_take_the_embedding_files_route() {
        let a = cli(&[
            "--go",
            "--feature-embedding",
            "f.parquet",
            "--cell-embedding",
            "c.parquet",
        ]);
        assert!(matches!(route(&a, None), Route::Enrichment));
    }

    #[test]
    fn projection_scores_where_the_run_decodes_and_trains_units_otherwise() {
        let bge = run_with(RunKind::Bge, true, true);
        let fne = run_with(RunKind::Fne, true, false);
        let pinto = run_with(RunKind::Other("pinto-cage".into()), false, true);
        let p = |more: &[&str], m: Option<&RunManifest>| resolve_projection(&cli(more), m);
        assert_eq!(p(&[], Some(&bge)), ProjectionMethod::Score);
        assert_eq!(
            p(&[], Some(&fne)),
            ProjectionMethod::Units,
            "no decoder table"
        );
        assert_eq!(p(&[], Some(&pinto)), ProjectionMethod::Units);
        assert_eq!(
            p(&["--min-markers", "5"], Some(&bge)),
            ProjectionMethod::Score,
            "settings no longer steer the choice"
        );
        assert_eq!(
            p(&[], None),
            ProjectionMethod::Units,
            "files: the table is unknown"
        );
        assert_eq!(
            p(&["--projection", "score"], None),
            ProjectionMethod::Score,
            "asked"
        );
        assert_eq!(
            resolve_method(&cli(&[]), Some(&pinto)),
            AnnotateMethod::Enrichment,
            "auto runs every method, the units beside the enrichment"
        );
    }

    #[test]
    fn auto_runs_every_method_and_a_chosen_one_runs_alone() {
        let auto = build_enrichment_args(&cli(&[]));
        assert!(auto.score_beside && auto.units_beside);
        assert_eq!(
            auto.call, None,
            "SuSiE's call where it ran, else the enrichment's"
        );
        let go = build_enrichment_args(&cli(&["--go"]));
        assert!(
            go.units_beside,
            "the units align the GO terms the enrichment tests"
        );
        let alone = build_enrichment_args(&cli(&["--method", "enrichment"]));
        assert!(!alone.score_beside && !alone.units_beside);
        let call = build_enrichment_args(&cli(&["--call", "units"]));
        assert_eq!(call.call, Some(crate::annotate::rounds::Method::Units));
    }

    /// A run in a fresh directory: `kind`, its gene tables as `outputs`
    /// says, its counts `here` or not.
    fn loaded_run(
        dir: &std::path::Path,
        kind: &str,
        outputs: serde_json::Value,
        here: bool,
    ) -> crate::manifest::run::Loaded {
        let at = |n: &str| dir.join(n).to_string_lossy().into_owned();
        if here {
            std::fs::write(at("counts.zarr.zip"), "").unwrap();
        }
        let m = serde_json::json!({
            "version": 2, "kind": kind, "prefix": "run",
            "data": {"input": ["counts.zarr.zip"]},
            "outputs": outputs,
        });
        std::fs::write(at("run.senna.json"), m.to_string()).unwrap();
        crate::manifest::run::load(&at("run.senna.json")).unwrap()
    }

    #[test]
    fn auto_without_counts_projects_where_the_run_has_a_gene_table() {
        let dir = tempfile::tempdir().unwrap();
        let fne = serde_json::json!({"feature_embedding": "run.feature_embedding.parquet"});
        let gone = loaded_run(dir.path(), "fne", fne.clone(), false);
        assert_eq!(
            resolve_method_on(&cli(&[]), Some(&gone)),
            AnnotateMethod::Projection,
            "the units alone: no counts for the enrichment"
        );
        assert_eq!(
            resolve_method_on(&cli(&["--method", "enrichment"]), Some(&gone)),
            AnnotateMethod::Enrichment,
            "a choice stands"
        );
        let here = tempfile::tempdir().unwrap();
        let counted = loaded_run(here.path(), "fne", fne, true);
        assert_eq!(
            resolve_method_on(&cli(&[]), Some(&counted)),
            AnnotateMethod::Enrichment
        );
        let bare = tempfile::tempdir().unwrap();
        let svd = loaded_run(bare.path(), "svd", serde_json::json!({}), false);
        assert_eq!(
            resolve_method_on(&cli(&[]), Some(&svd)),
            AnnotateMethod::Enrichment,
            "nothing to project on: the enrichment says why it cannot run"
        );
    }

    #[test]
    fn a_call_by_a_method_that_does_not_run_is_refused() {
        use crate::annotate::rounds::Method;
        let dir = tempfile::tempdir().unwrap();
        // Decodes cells (the score), but no co-embedding (no units).
        let bge = loaded_run(
            dir.path(),
            "bge",
            serde_json::json!({"feature_embedding": "run.feature_embedding.parquet"}),
            true,
        );
        let ok = |more: &[&str]| check_call(&cli(more), Some(&bge));
        assert!(ok(&["--call", "score"]).is_ok());
        assert!(ok(&["--call", "susie"]).is_ok());
        assert!(ok(&["--call", "enrichment"]).is_ok());
        assert!(ok(&[]).is_ok());
        for (more, named) in [
            (&["--call", "units"][..], "units"),
            (&["--call", "susie", "--no-susie"][..], "susie"),
            (&["--call", "score", "--method", "enrichment"][..], "score"),
            (
                &["--call", "enrichment", "--projection", "score"][..],
                "enrichment",
            ),
        ] {
            let e = ok(more).unwrap_err();
            let e = format!("{e:#}");
            assert!(e.contains(&format!("--call {named}")), "{more:?}: {e}");
        }
        assert_eq!(
            methods_run(&cli(&[]), Some(&bge)),
            [Method::Enrichment, Method::Susie, Method::Score]
        );
    }
}

#[cfg(test)]
#[path = "tests/annotate_cmd.rs"]
mod argv_tests;
