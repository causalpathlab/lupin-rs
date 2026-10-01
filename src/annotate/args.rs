//! Settings for the three annotation passes, built from the `lupin annotate`
//! flags. Values nobody tunes per run are the constants below; every pass
//! records its settings (these included) in `manifest.annotate.settings`.

/// Keep IEA (electronic) GAF annotations.
pub const KEEP_IEA: bool = true;
/// Cells per block when streaming the raw counts.
pub const BLOCK_SIZE: usize = 1024;
/// Random gene sets in the gene-set null: its p-values go no lower than `1 / (NUM_DRAWS + 1)`.
pub const NUM_DRAWS: usize = 10_000;
/// Minimum per-cell label confidence (0 = keep every call).
pub const MIN_CONFIDENCE: f32 = 0.0;

/// The constants above, for the manifest record.
#[must_use]
pub fn fixed_settings() -> serde_json::Value {
    serde_json::json!({
        "keep_iea": KEEP_IEA,
        "block_size": BLOCK_SIZE,
        "num_draws": NUM_DRAWS,
        "min_confidence": MIN_CONFIDENCE,
        "empirical_specificity": true,
        "gene_strata": true,
        "recluster": true,
        "panel_perm": 0,
        "support_perm": 0,
        "effect": "nes",
        "share": "softmax of z = probit(1 - p) over the FDR survivors",
        "p": "gene-set null permutation p (fgsea sign-aware)",
        "q": "treebh over the panel's cell-type tree",
    })
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct AnnotateArgs {
    /// Cluster parquet (cells × 1 cluster column); overrides `manifest.cluster.clusters`
    pub clusters: Option<Box<str>>,

    /// The pinto cascade level `clusters` came from (`--level`)
    pub level: Option<Box<str>>,

    /// Nearest neighbors for the internal Leiden cosine-KNN graph
    pub knn: usize,

    /// Modularity resolution (CPM) for the internal Leiden clustering
    pub resolution: f64,

    /// Optional target cluster count for Leiden auto-resolution
    pub num_clusters: Option<usize>,

    /// Minimum cluster size; smaller clusters become unassigned
    pub min_cluster_size: usize,

    /// Seed for internal Leiden clustering (deterministic)
    pub cluster_seed: Option<u64>,

    /// Marker-gene TSV, `gene<TAB>celltype` per line; empty in gene-set mode
    pub markers: Box<str>,

    /// GO annotation file (.gaf/.gaf.gz)
    pub gaf: Option<Box<str>>,

    /// MSigDB GMT gene-sets (`term<TAB>desc<TAB>genes…`); gene-set mode, like `gaf`
    pub gmt: Option<Box<str>>,

    /// Score GO terms with the GO annotations of the data's species (told
    /// from its gene names); `gaf` names the file instead
    pub go: bool,

    /// Gene Ontology .obo naming the terms of `go`/`gaf`/`gmt`
    pub go_obo: Option<Box<str>>,

    /// GO/GMT terms sharing fewer genes with the data's features are dropped
    /// (the overlap, not the term's own size)
    pub go_min_overlap: usize,

    /// GO/GMT terms sharing more genes with the data's features are dropped
    pub go_max_overlap: usize,

    /// Output prefix for annotation artifacts
    pub out: Box<str>,

    /// Number of PB-level sample permutations for the correlation-preserving null
    pub num_perm: usize,

    /// Drop a cell type with fewer than this many matched markers
    pub min_markers: usize,

    /// FDR α for the Q-matrix threshold
    pub fdr_alpha: f32,

    /// Softmax temperature used when row-normalizing Q over significant entries
    pub q_temperature: f32,

    /// RNG seed (deterministic; affects row randomization)
    pub seed: u64,

    /// Keep existing {out}.* annotation outputs.
    /// By default the explicit annotation set is erased first,
    /// never the embedding or manifest, for a fresh re-run.
    pub no_clean: bool,

    // ── optional inline ontology annotation (TreeBH) ──
    /// Cell Ontology .obo, e.g. cl-basic.obo
    pub obo: Option<Box<str>>,

    /// Curated `label<TAB>CL:id` map, one row per marker celltype.
    /// Needed with `obo`
    pub label_cl: Option<Box<str>>,

    /// Ontology TreeBH per-level FDR target (lower → descends less, abstains more)
    pub ontology_fdr_q: f64,

    /// Ontology TreeBH: Benjamini–Yekutieli within families (any dependence; more conservative)
    pub ontology_by: bool,
}

/// `lupin annotate --method projection` — firm marker-set annotation by projection
/// onto a co-embedded feature space (bge / fne / resolve-embedding-space).
/// Embedding-grounded (no raw-count re-read), complementary to
/// `annotate --method enrichment`. Drives the shared firm term-ORA core.
#[derive(Debug, serde::Serialize)]
pub struct AnnotateProjectionArgs {
    /// Marker-gene TSV: `gene<TAB>celltype` per line (tab/comma/space delimited)
    pub markers: Box<str>,

    /// Output prefix
    pub out: Box<str>,

    /// k for the cosine cell kNN graph fed to Leiden clustering
    pub knn: usize,

    /// Leiden resolution for cell clustering (higher → more, finer clusters)
    pub resolution: f64,

    /// Permutation draws calibrating the over-representation null (0 = analytic p only)
    pub num_perm: usize,

    /// RNG seed (clustering + permutation null)
    pub seed: u64,

    /// Drop a cell type with fewer than this many usable markers
    pub min_markers: usize,

    /// Disable IDF down-weighting of markers shared across many types
    pub no_idf: bool,

    /// Keep every cell→term assignment (skip the distance-outlier prune)
    pub no_assign_qc: bool,

    /// Outlier gate:
    /// prune a cell whose distance to its centroid exceeds median + k·MAD
    pub assign_mad: f64,

    /// FDR α for the per-cluster term call + Q sparsity (BH on the permutation p)
    pub fdr_alpha: f32,

    /// Softmax temperature when row-normalizing Q over significant terms
    pub q_temperature: f32,

    /// Cell Ontology .obo, e.g. cl-basic.obo
    pub obo: Option<Box<str>>,

    /// Curated `label<TAB>CL:id` map, one row per marker celltype.
    /// Needed with `obo`
    pub label_cl: Option<Box<str>>,

    /// Ontology TreeBH per-level FDR target (lower → descends less, abstains more)
    pub ontology_fdr_q: f64,

    /// Ontology TreeBH: Benjamini–Yekutieli within families (any dependence; more conservative)
    pub ontology_by: bool,

    /// Keep existing {out}.* projection outputs (default: erase the explicit set first)
    pub no_clean: bool,
}

/// Ontology follow-up of `lupin annotate` — hierarchical multi-resolution cell-type calling
/// (TreeBH) on the Cell Ontology, post-processing an `annotate --method enrichment`
/// run's cluster × celltype matrix.
#[derive(Debug, serde::Serialize)]
pub struct AnnotateOntologyArgs {
    /// Curated `label<TAB>CL:id` TSV mapping celltypes to Cell Ontology terms
    pub label_cl: Box<str>,

    /// Cell Ontology OBO file (e.g. cl-basic.obo)
    pub obo: Box<str>,

    /// Output prefix
    pub out: Box<str>,

    /// Per-level selective-FDR target (TreeBH).
    /// Lower → descends less, abstains more
    pub fdr_q: f64,

    /// Benjamini–Yekutieli within families (any dependence; more conservative)
    pub by: bool,

    /// Force the (saturated) permutation p-values instead of the default z→p
    pub use_perm_p: bool,
}
