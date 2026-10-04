//! Reader/writer for the `{prefix}.senna.json` run manifest a senna fit writes.
//!
//! Lupin models only the fields it reads or records. Every section keeps the
//! rest in an `extra` map, so a manifest round-trips through `load` → `save`
//! without losing anything senna (or a newer senna) put there. A run kind this
//! build does not know still loads, as [`RunKind::Other`].
//!
//! `--from` may name the manifest file itself or the run's output prefix; every
//! command resolves it through [`load`], which also remembers the file so
//! updates are saved back where they were read. Annotation instead writes a
//! new manifest through [`Loaded::copy_to`].

use std::fs;
use std::path::{Path, PathBuf};

use legume_numeric::matrix::dense_mat_io::{l2_normalize_rows_inplace, Mat, MatWithNames};
use legume_numeric::matrix::traits::{IoOps, MatOps};
use log::info;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};

/// Newest manifest layout this build understands.
pub const MANIFEST_VERSION: u32 = 2;

type Extra = Map<String, Value>;

/// Which senna command produced the run. Serialized as its kebab-case name.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(from = "String", into = "String")]
pub enum RunKind {
    Topic,
    Itopic,
    MaskedVae,
    JointTopic,
    Vae,
    Svd,
    JointSvd,
    Bge,
    Fne,
    ResolveEmbeddingSpace,
    Gem,
    Simba,
    /// A kind this build does not know. Treated conservatively: signed cell
    /// space, no log-simplex latent, no frozen gene table.
    Other(String),
}

/// What the table [`RunOutputs::geometry_latent`] holds.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum CellSpace {
    /// Euclidean coordinates where magnitude carries signal; angular distance fits.
    Embedding,
    /// Rows are `log θ`: `exp()` gives a probability vector.
    LogSimplex,
    /// Signed scores (loadings, a Gaussian `z`): no simplex, no magnitude semantics.
    Signed,
}

impl RunKind {
    #[must_use]
    pub fn as_str(&self) -> &str {
        match self {
            RunKind::Topic => "topic",
            RunKind::Itopic => "itopic",
            RunKind::MaskedVae => "masked-vae",
            RunKind::JointTopic => "joint-topic",
            RunKind::Vae => "vae",
            RunKind::Svd => "svd",
            RunKind::JointSvd => "joint-svd",
            RunKind::Bge => "bge",
            RunKind::Fne => "fne",
            RunKind::ResolveEmbeddingSpace => "resolve-embedding-space",
            RunKind::Gem => "gem",
            RunKind::Simba => "simba",
            RunKind::Other(s) => s,
        }
    }

    #[must_use]
    pub fn cell_space(&self) -> CellSpace {
        match self {
            RunKind::Bge
            | RunKind::Fne
            | RunKind::ResolveEmbeddingSpace
            | RunKind::Gem
            | RunKind::Simba => CellSpace::Embedding,
            RunKind::Topic | RunKind::Itopic | RunKind::JointTopic => CellSpace::LogSimplex,
            RunKind::MaskedVae
            | RunKind::Vae
            | RunKind::Svd
            | RunKind::JointSvd
            | RunKind::Other(_) => CellSpace::Signed,
        }
    }

    /// The whole gene-side model is a frozen `gene × H` table in `feature_embedding`.
    #[must_use]
    pub fn has_frozen_gene_table(&self) -> bool {
        matches!(self, RunKind::Bge | RunKind::Simba | RunKind::Gem)
    }

    /// Kinds that write a co-embedded gene table (`feature_coembedding`) beside ρ.
    #[must_use]
    pub fn coembeds(&self) -> bool {
        matches!(
            self,
            RunKind::Bge | RunKind::Gem | RunKind::Simba | RunKind::ResolveEmbeddingSpace
        )
    }
}

impl From<String> for RunKind {
    fn from(s: String) -> Self {
        match s.as_str() {
            "topic" => RunKind::Topic,
            "itopic" => RunKind::Itopic,
            "masked-vae" => RunKind::MaskedVae,
            "joint-topic" => RunKind::JointTopic,
            "vae" => RunKind::Vae,
            "svd" => RunKind::Svd,
            "joint-svd" => RunKind::JointSvd,
            "bge" => RunKind::Bge,
            "fne" => RunKind::Fne,
            "resolve-embedding-space" => RunKind::ResolveEmbeddingSpace,
            "gem" => RunKind::Gem,
            "simba" => RunKind::Simba,
            _ => RunKind::Other(s),
        }
    }
}

impl From<RunKind> for String {
    fn from(k: RunKind) -> Self {
        match k {
            RunKind::Other(s) => s,
            k => k.as_str().to_string(),
        }
    }
}

impl std::fmt::Display for RunKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunManifest {
    pub version: u32,
    pub kind: RunKind,
    /// The `--out` prefix the training command was run with.
    pub prefix: String,
    #[serde(default)]
    pub data: RunData,
    #[serde(default)]
    pub outputs: RunOutputs,
    #[serde(default)]
    pub layout: RunLayout,
    #[serde(default)]
    pub cluster: RunCluster,
    #[serde(default)]
    pub annotate: RunAnnotate,
    #[serde(default)]
    pub defaults: RunDefaults,
    #[serde(default)]
    pub trajectory: RunTrajectory,
    #[serde(flatten)]
    pub extra: Extra,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RunData {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub input: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub batch: Vec<String>,
    /// The multiome layout `input` was loaded under, positional against `input`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub multiome: Option<RunMultiome>,
    #[serde(flatten)]
    pub extra: Extra,
}

/// A multiome run's per-file layout, positional against `data.input`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RunMultiome {
    /// Modality label per input file; features become `{name}/{modality}`.
    pub modality: Vec<String>,
    /// Sample-group label per input file.
    pub group: Vec<String>,
    /// Whether barcodes were tagged `{barcode}@{group}` at load.
    pub barcode_tagged: bool,
    #[serde(flatten)]
    pub extra: Extra,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RunOutputs {
    /// Cell × K: log θ for topic runs, component scores for SVD runs.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub latent: Option<String>,
    /// Gene × K: signed loadings (SVD family) or a pre-split topic dictionary.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dictionary: Option<String>,
    /// Gene × K topic dictionary β, log space.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub softmax_dictionary: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dictionary_empirical: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anchor_labels: Option<String>,
    /// Gene × H raw gene embedding ρ.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub feature_embedding: Option<String>,
    /// Gene × H genes placed on the cell manifold (embedding kinds).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub feature_coembedding: Option<String>,
    /// Cell × H Euclidean embedding Z.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cell_embedding: Option<String>,
    #[serde(flatten)]
    pub extra: Extra,
}

impl RunOutputs {
    /// The topic dictionary or, failing that, the SVD loadings.
    #[must_use]
    pub fn gene_dictionary(&self) -> Option<&str> {
        self.softmax_dictionary
            .as_deref()
            .or(self.dictionary.as_deref())
    }

    /// The cell table for GEOMETRY (kNN, layout, clustering):
    /// `cell_embedding`, else `latent`.
    #[must_use]
    pub fn geometry_latent(&self) -> Option<&str> {
        self.cell_embedding.as_deref().or(self.latent.as_deref())
    }

    /// The cell table for a COMPOSITION view (structure bars, topic colour):
    /// `latent`, else `cell_embedding`.
    #[must_use]
    pub fn structure_latent(&self) -> Option<&str> {
        self.latent.as_deref().or(self.cell_embedding.as_deref())
    }

    /// The gene table paired with [`Self::structure_latent`].
    #[must_use]
    pub fn structure_dictionary(&self) -> Option<&str> {
        self.dictionary_empirical
            .as_deref()
            .or_else(|| self.gene_dictionary())
            .or(self.feature_embedding.as_deref())
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RunLayout {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cell_coords: Option<String>,
    #[serde(flatten)]
    pub extra: Extra,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RunCluster {
    /// Cells × 1 cluster id.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub clusters: Option<String>,
    #[serde(flatten)]
    pub extra: Extra,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RunAnnotate {
    /// N × C cell posterior.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annotation: Option<String>,
    /// Per-cell label + max probability.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub argmax: Option<String>,
    /// nClusters × C FDR-sparse softmax Q (probabilities, not q-values).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cluster_celltype_q: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cluster_celltype_es: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cluster_expression: Option<String>,
    /// Input marker TSV (provenance).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub markers: Option<String>,
    /// The manifest this one was copied from when annotate wrote it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
    /// Per-cluster digest for review, keyed by cluster id.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cluster_summary: Option<String>,
    /// Every round's decisions per cluster id, newest first.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub history: Option<String>,
    /// The decisions that made this round (JSONL).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub log: Option<String>,
    /// The coarse cell-type groups the first round was called with.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub celltype_tree: Option<String>,
    /// The fine per-cell labels behind a coarse first round.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fine_argmax: Option<String>,
    /// The cluster table `cluster_expression`'s columns refer to: the ids of
    /// the pass that wrote it, which later rounds' merges do not change.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expression_clusters: Option<String>,
    /// What an enrichment pass cached so later rounds can be rescored
    /// without re-reading counts.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stats_cache: Option<StatsCache>,
    /// How this round's statistics were made (e.g. post-selection after
    /// curation), for readers to caveat them.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stats: Option<Value>,
    /// Every round's marker edits per cell type, newest first.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub marker_history: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ontology_assignment: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ontology_node_mass: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ontology_signature: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ontology_term_effect: Option<String>,
    /// Enrichment: nClusters × C FDR q-values.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cluster_celltype_q_values: Option<String>,
    /// Enrichment: nClusters × C p-values.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cluster_celltype_p: Option<String>,
    /// Enrichment: nClusters × C NES (fgsea's normalized enrichment score).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cluster_celltype_nes: Option<String>,
    /// Projection: nClusters × term FDR q-values.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cluster_term_q: Option<String>,
    /// Projection: the marker panel's gene embedding per type.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub marker_embedding: Option<String>,
    /// Effective settings of the latest run of each annotate method
    /// (`enrichment` / `projection` / `ontology`), constants included.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub settings: Option<Value>,
    #[serde(flatten)]
    pub extra: Extra,
}

/// An enrichment pass's sufficient statistics, manifest-relative: the raw
/// per-cluster gene sums (columns `K{id}` of `annotate.expression_clusters`),
/// the per-batch profile, the per-gene weights and each cell's batch.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StatsCache {
    pub gene_sum: String,
    pub batch_profile: String,
    pub gene_weight: String,
    pub cell_batch: String,
}

/// What `lupin trajectory` wrote, manifest-relative.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RunTrajectory {
    /// The combined prior statements and direct edges (TSV).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub prior: Option<String>,
    /// Type pairs: connectivity, in-prior, verdict, order agreement.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub edges: Option<String>,
    /// Per cell: pseudotime, type, component, lineage weights.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pseudotime: Option<String>,
    /// Cells × diffusion components.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub diffusion: Option<String>,
    /// The root-to-leaf paths (TSV).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lineages: Option<String>,
    /// The settings and inputs the run used.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub settings: Option<Value>,
    #[serde(flatten)]
    pub extra: Extra,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RunDefaults {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub colour_by: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub palette: Option<String>,
    #[serde(flatten)]
    pub extra: Extra,
}

impl RunManifest {
    /// A bare manifest stating only the kind.
    #[must_use]
    pub fn new(kind: RunKind, prefix: &str) -> Self {
        Self {
            version: MANIFEST_VERSION,
            kind,
            prefix: prefix.into(),
            data: RunData::default(),
            outputs: RunOutputs::default(),
            layout: RunLayout::default(),
            cluster: RunCluster::default(),
            annotate: RunAnnotate::default(),
            trajectory: RunTrajectory::default(),
            defaults: RunDefaults::default(),
            extra: Extra::default(),
        }
    }

    /// The run's input files (`data.input`), found from `manifest_dir` (see
    /// [`Self::data_file`]).
    #[must_use]
    pub fn data_inputs(&self, manifest_dir: &Path) -> Vec<String> {
        self.data
            .input
            .iter()
            .map(|p| self.data_file(manifest_dir, p))
            .collect()
    }

    /// The run's batch files (`data.batch`), found as [`Self::data_inputs`].
    #[must_use]
    pub fn data_batches(&self, manifest_dir: &Path) -> Vec<String> {
        self.data
            .batch
            .iter()
            .map(|p| self.data_file(manifest_dir, p))
            .collect()
    }

    /// A data file the run recorded, found here even when the run was
    /// trained on another machine or the tree has moved: as recorded
    /// ([`resolve`]); else at the same place relative to the run's
    /// directory as it was relative to the training `prefix`'s; else under
    /// the manifest's directory or an ancestor of it, by the longest tail of
    /// the recorded path that exists there. Unfound, as recorded.
    #[must_use]
    pub fn data_file(&self, manifest_dir: &Path, recorded: &str) -> String {
        let given = resolve(manifest_dir, recorded);
        if Path::new(&given).exists() {
            return given;
        }
        // Lexically absolute, not canonical: the recorded paths are as
        // written, so resolving symlinks here alone would skew the rebase.
        let here = normalize(
            &std::path::absolute(manifest_dir).unwrap_or_else(|_| manifest_dir.to_path_buf()),
        );
        let recorded_abs = normalize(Path::new(recorded));
        let trained_in = Path::new(&self.prefix)
            .is_absolute()
            .then(|| normalize(Path::new(&self.prefix)))
            .and_then(|p| p.parent().map(Path::to_path_buf));
        let rebased = trained_in
            .filter(|_| recorded_abs.is_absolute())
            .map(|old| normalize(&here.join(relative_to(&recorded_abs, &old))))
            .filter(|p| p.exists());
        let found = rebased.or_else(|| {
            let parts: Vec<_> = recorded_abs
                .components()
                .filter(|c| matches!(c, std::path::Component::Normal(_)))
                .collect();
            // At least the file and its folder must match, so an unrelated
            // file of the same name elsewhere is not taken for the data.
            let min = parts.len().min(2);
            (0..=parts.len().checked_sub(min.max(1))?).find_map(|skip| {
                let tail: PathBuf = parts[skip..].iter().collect();
                here.ancestors().map(|a| a.join(&tail)).find(|p| p.exists())
            })
        });
        match found {
            Some(p) => {
                let p = p.to_string_lossy().into_owned();
                log::warn!("{recorded} is not here; using {p}");
                p
            }
            None => given,
        }
    }

    /// Read a manifest file; returns it with its directory, against which the
    /// relative paths inside resolve.
    pub fn load(path: &Path) -> anyhow::Result<(Self, PathBuf)> {
        let raw = fs::read_to_string(path)
            .map_err(|e| anyhow::anyhow!("read {}: {e}", path.display()))?;
        let mut m: Self = serde_json::from_str(&raw)
            .map_err(|e| anyhow::anyhow!("parse {}: {e}", path.display()))?;
        m.lift_v1_feature_slots();
        if m.version > MANIFEST_VERSION {
            log::warn!(
                "manifest {} is v{} but lupin reads up to v{MANIFEST_VERSION}; proceeding",
                path.display(),
                m.version
            );
        }
        // A `pinto-*` kind is a pinto run lupin annotated, not a stranger.
        if let RunKind::Other(k) = &m.kind {
            if !k.starts_with("pinto-") {
                log::warn!(
                    "manifest {}: run kind `{k}` is unknown to lupin; treating its cells as signed scores",
                    path.display()
                );
            }
        }
        Ok((m, parent_dir(path)))
    }

    /// v1 embedding runs stored the co-embed as `feature_embedding` and ρ as
    /// `feature_loading`; move them onto the v2 slots.
    fn lift_v1_feature_slots(&mut self) {
        if self.version >= 2 {
            return;
        }
        if self.kind.coembeds() && self.outputs.feature_coembedding.is_none() {
            self.outputs.feature_coembedding = self.outputs.feature_embedding.take();
        }
        if let Some(Value::String(rho)) = self.outputs.extra.remove("feature_loading") {
            self.outputs.feature_embedding = Some(rho);
        }
    }

    /// Write the manifest to `path`. An existing file is replaced only as
    /// [`may_replace`] allows.
    pub fn save(&self, path: &Path) -> anyhow::Result<()> {
        may_replace(path)?;
        let s = serde_json::to_string_pretty(self)?;
        fs::write(path, s).map_err(|e| anyhow::anyhow!("write {}: {e}", path.display()))?;
        wrote(path);
        info!("wrote {}", path.display());
        Ok(())
    }
}

/// Turn a path written relative to the working directory into the
/// manifest-relative form for storage; paths outside the manifest directory
/// stay absolute.
#[must_use]
pub fn rel_to_manifest(manifest_dir: &Path, written_path: &str) -> String {
    let Ok(cwd) = std::env::current_dir() else {
        return written_path.to_string();
    };
    let abs = cwd.join(written_path);
    let manifest_dir = cwd.join(manifest_dir);
    let manifest_abs = manifest_dir.canonicalize().unwrap_or(manifest_dir);
    // A stem or a file not yet written: canonicalize its directory instead.
    let written_abs = abs.canonicalize().unwrap_or_else(|_| {
        match (
            abs.parent().and_then(|d| d.canonicalize().ok()),
            abs.file_name(),
        ) {
            (Some(d), Some(name)) => d.join(name),
            _ => abs.clone(),
        }
    });
    match written_abs.strip_prefix(&manifest_abs) {
        Ok(rel) => rel.to_string_lossy().into_owned(),
        Err(_) => written_abs.to_string_lossy().into_owned(),
    }
}

/// The directory a file is in; `.` for a bare file name.
#[must_use]
pub fn parent_dir(p: &Path) -> PathBuf {
    p.parent()
        .filter(|d| !d.as_os_str().is_empty())
        .map_or_else(|| PathBuf::from("."), Path::to_path_buf)
}

/// Whether two paths name the same file. Paths that cannot be resolved (a
/// file not written yet) are compared as written.
#[must_use]
pub fn same_file(a: &Path, b: &Path) -> bool {
    let canon = |p: &Path| p.canonicalize().unwrap_or_else(|_| p.to_path_buf());
    canon(a) == canon(b)
}

/// A cell-coordinates table: the cell names (the parquet's row labels, in
/// order) and each column by name.
pub type CellCoords = (Vec<Box<str>>, rustc_hash::FxHashMap<String, Vec<f32>>);

/// The cells × columns table at `path` (a `senna layout` output), by column.
pub fn read_cell_coords(path: &str) -> anyhow::Result<CellCoords> {
    let MatWithNames { rows, cols, mat } = Mat::from_parquet(path)?;
    let by_name = cols
        .iter()
        .enumerate()
        .map(|(j, name)| (name.to_string(), mat.column(j).iter().copied().collect()))
        .collect();
    Ok((rows, by_name))
}

/// Resolve a manifest-relative path against the manifest's directory.
/// Absolute paths pass through.
#[must_use]
pub fn resolve(manifest_dir: &Path, rel: &str) -> String {
    let p = Path::new(rel);
    if p.is_absolute() {
        rel.to_string()
    } else {
        manifest_dir.join(p).to_string_lossy().into_owned()
    }
}

/// Set by `--overwrite`: replace existing manifests without asking.
static OVERWRITE: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// Manifests this process wrote or was allowed to replace; it may write them again.
static WRITTEN: std::sync::Mutex<std::collections::BTreeSet<PathBuf>> =
    std::sync::Mutex::new(std::collections::BTreeSet::new());

/// Let every later [`may_replace`] replace existing manifests (`--overwrite`).
pub fn allow_overwrite() {
    OVERWRITE.store(true, std::sync::atomic::Ordering::Relaxed);
}

/// `path` canonical, or as given when it cannot be (it does not exist).
fn key(path: &Path) -> PathBuf {
    path.canonicalize().unwrap_or_else(|_| path.to_path_buf())
}

fn wrote(path: &Path) {
    if let Ok(mut w) = WRITTEN.lock() {
        w.insert(key(path));
    }
}

/// Whether `path` may be written without asking: it does not exist, this
/// process wrote it (or was allowed to) already, or `--overwrite` is on.
#[must_use]
pub fn may_replace_unasked(path: &Path) -> bool {
    !path.exists()
        || OVERWRITE.load(std::sync::atomic::Ordering::Relaxed)
        || WRITTEN.lock().is_ok_and(|w| w.contains(&key(path)))
}

/// Whether `path` may be written: yes when it does not exist, when this
/// process wrote it (or was allowed to) already, or under `--overwrite`.
/// Otherwise ask on the terminal; with no terminal to ask on, refuse.
pub fn may_replace(path: &Path) -> anyhow::Result<()> {
    use std::io::{IsTerminal, Write};
    if may_replace_unasked(path) {
        return Ok(());
    }
    let refuse = || {
        anyhow::anyhow!(
            "{} exists; not replacing it (pass --overwrite to replace it, or choose another -o)",
            path.display()
        )
    };
    if !(std::io::stdin().is_terminal() && std::io::stderr().is_terminal()) {
        return Err(refuse());
    }
    eprint!("lupin: {} exists. Replace it? [y/N] ", path.display());
    std::io::stderr().flush()?;
    let mut answer = String::new();
    std::io::stdin().read_line(&mut answer)?;
    if matches!(answer.trim().to_lowercase().as_str(), "y" | "yes") {
        wrote(path);
        Ok(())
    } else {
        Err(refuse())
    }
}

/// `p` with `.` and `..` folded away, without touching the file system.
fn normalize(p: &Path) -> PathBuf {
    use std::path::Component;
    let mut out = PathBuf::new();
    for c in p.components() {
        match c {
            Component::CurDir => {}
            // `..` cancels a name before it; after nothing, a root or
            // another `..` it stays.
            Component::ParentDir => match out.components().next_back() {
                Some(Component::Normal(_)) => {
                    out.pop();
                }
                Some(Component::RootDir | Component::Prefix(_)) => {}
                _ => out.push(".."),
            },
            c => out.push(c),
        }
    }
    out
}

/// `p` as a path from `base`, both absolute and normalized.
fn relative_to(p: &Path, base: &Path) -> PathBuf {
    let (pc, bc): (Vec<_>, Vec<_>) = (p.components().collect(), base.components().collect());
    let common = pc.iter().zip(&bc).take_while(|(a, b)| a == b).count();
    let mut out: PathBuf = bc[common..].iter().map(|_| "..").collect();
    out.extend(&pc[common..]);
    out
}

/// Rewrite every relative path in `value` that names an existing file or
/// directory from `from_dir` so it names the same one from `to_dir`. The walk
/// covers keys lupin does not model, so another tool's paths survive a move;
/// strings that resolve to nothing are left alone.
pub fn rebase_paths(value: &mut Value, from_dir: &Path, to_dir: &Path) {
    match value {
        Value::String(s) => {
            if Path::new(s.as_str()).is_absolute() || s.is_empty() {
                return;
            }
            let old = from_dir.join(s.as_str());
            if old.exists() {
                *s = rel_to_manifest(to_dir, &old.to_string_lossy());
            }
        }
        Value::Array(items) => items
            .iter_mut()
            .for_each(|v| rebase_paths(v, from_dir, to_dir)),
        Value::Object(map) => map
            .values_mut()
            .for_each(|v| rebase_paths(v, from_dir, to_dir)),
        _ => {}
    }
}

/// The run's output prefix from `--from`: strip `.senna.json` or a trailing `.json`.
#[must_use]
pub fn derive_out_prefix(from: &str) -> String {
    from.strip_suffix(".senna.json")
        .or_else(|| from.strip_suffix(LUPIN_SUFFIX))
        .or_else(|| from.strip_suffix(super::pinto::SUFFIX))
        .or_else(|| from.strip_suffix(".json"))
        .unwrap_or(from)
        .to_string()
}

/// Default manifest filename for a run `--out` prefix.
#[must_use]
pub fn default_path(prefix: &str) -> String {
    format!("{prefix}.senna.json")
}

/// The suffix of a manifest lupin writes for a run senna did not train (a
/// pinto run). Same layout as a senna manifest; senna does not read it.
pub const LUPIN_SUFFIX: &str = ".lupin.json";

/// Where annotation writes its manifest for `-o prefix`: `.senna.json` for a
/// senna run, `.lupin.json` for anything else (a pinto run, or a round of one).
#[must_use]
pub fn annotated_path(source: &Path, prefix: &str) -> PathBuf {
    let name = source.to_string_lossy();
    if name.ends_with(LUPIN_SUFFIX) || super::pinto::is_pinto(source) {
        PathBuf::from(format!("{prefix}{LUPIN_SUFFIX}"))
    } else {
        PathBuf::from(default_path(prefix))
    }
}

/// The manifest file `--from` names: the path itself when it is a file,
/// otherwise `{prefix}.senna.json`, else `{prefix}.lupin.json` or
/// `{prefix}.pinto.json`, whichever exists first.
#[must_use]
pub fn manifest_file(from: &str) -> PathBuf {
    let direct = Path::new(from);
    if direct.is_file() {
        return direct.to_path_buf();
    }
    let prefix = derive_out_prefix(from);
    let senna = PathBuf::from(default_path(&prefix));
    [
        PathBuf::from(format!("{prefix}{LUPIN_SUFFIX}")),
        PathBuf::from(format!("{prefix}{}", super::pinto::SUFFIX)),
    ]
    .into_iter()
    .find(|p| !senna.is_file() && p.is_file())
    .unwrap_or(senna)
}

/// A loaded manifest, its directory (for resolving relative paths), and the
/// file it came from (for saving back).
pub struct Loaded {
    pub manifest: RunManifest,
    pub dir: PathBuf,
    pub file: PathBuf,
}

/// The geometry table of the run `manifest` describes, read and prepared for
/// kNN distances by the run's cell space: rows L2-normalised for an embedding
/// (Euclidean distance then ranks as cosine), `exp` then column z-scores for
/// `log θ`, column z-scores otherwise. Annotation's Leiden step and
/// `lupin trajectory` both feed their kNN graphs this.
pub fn prepare_geometry(manifest: &RunManifest, dir: &Path) -> anyhow::Result<MatWithNames<Mat>> {
    let rel = manifest.outputs.geometry_latent().ok_or_else(|| {
        anyhow::anyhow!("the manifest has neither `outputs.cell_embedding` nor `outputs.latent`")
    })?;
    let path = resolve(dir, rel);
    let mut x = Mat::from_parquet_with_row_names(&path, Some(0))
        .map_err(|e| anyhow::anyhow!("reading {path}: {e}"))?;
    match manifest.kind.cell_space() {
        CellSpace::Embedding => l2_normalize_rows_inplace(&mut x.mat),
        CellSpace::LogSimplex => {
            x.mat.apply(|v| *v = v.exp());
            x.mat.scale_columns_inplace();
        }
        CellSpace::Signed => x.mat.scale_columns_inplace(),
    }
    info!(
        "{} cells × {} dims from {path}",
        x.rows.len(),
        x.mat.ncols()
    );
    Ok(x)
}

impl Loaded {
    /// `{dir}/{name}` for a manifest at `{dir}/{name}.senna.json`: where the
    /// run's artifacts that the manifest does not record (NB-Fisher weights)
    /// sit. Derived from where the manifest IS, never from its
    /// `prefix` field, which holds the training machine's absolute path.
    #[must_use]
    pub fn run_prefix(&self) -> String {
        derive_out_prefix(&self.file.to_string_lossy())
    }

    /// The cell table for geometry (`cell_embedding`, else `latent`), resolved.
    pub fn geometry_latent_path(&self) -> anyhow::Result<String> {
        let rel = self.manifest.outputs.geometry_latent().ok_or_else(|| {
            anyhow::anyhow!(
                "{} has neither `outputs.cell_embedding` nor `outputs.latent`",
                self.file.display()
            )
        })?;
        Ok(resolve(&self.dir, rel))
    }

    /// The run's geometry table, read and prepared for kNN distances
    /// ([`prepare_geometry`]).
    pub fn prepared_geometry(&self) -> anyhow::Result<MatWithNames<Mat>> {
        prepare_geometry(&self.manifest, &self.dir)
    }
}

impl Loaded {
    /// This run as a new manifest at `file`, which must not be the one it was
    /// read from. Paths are rebased onto `file`'s directory; `prefix`, a stem
    /// rather than a file, keeps naming the original run, whose tables live
    /// there. Nothing is written.
    pub fn copy_to(&self, file: PathBuf) -> anyhow::Result<Loaded> {
        anyhow::ensure!(
            !same_file(&self.file, &file),
            "{} is the manifest this command reads from; choose a different --out",
            file.display()
        );
        let dir = parent_dir(&file);
        let mut value = serde_json::to_value(&self.manifest)?;
        let prefix = value.as_object_mut().and_then(|m| m.remove("prefix"));
        rebase_paths(&mut value, &self.dir, &dir);
        if let (Some(Value::String(p)), Some(m)) = (prefix, value.as_object_mut()) {
            let p = if Path::new(&p).is_absolute() {
                p
            } else {
                rel_to_manifest(&dir, &self.dir.join(&p).to_string_lossy())
            };
            m.insert("prefix".into(), Value::String(p));
        }
        let mut manifest: RunManifest = serde_json::from_value(value)?;
        manifest.annotate.source = Some(rel_to_manifest(&dir, &self.file.to_string_lossy()));
        Ok(Loaded {
            manifest,
            dir,
            file,
        })
    }
}

/// Load `--from`, given as a manifest path or a bare prefix.
pub fn load(from: &str) -> anyhow::Result<Loaded> {
    let file = manifest_file(from);
    if super::pinto::is_pinto(&file) {
        return super::pinto::load(&file);
    }
    let (manifest, dir) = RunManifest::load(&file).map_err(|e| {
        anyhow::anyhow!(
            "{e}\n\n`{from}` is neither a readable run manifest nor the prefix of one \
             (looked for `{}`)",
            file.display()
        )
    })?;
    info!(
        "Loaded run manifest {} (kind: {})",
        file.display(),
        manifest.kind
    );
    Ok(Loaded {
        manifest,
        dir,
        file,
    })
}

/// Load `--from` when given; otherwise no manifest and paths resolve against `.`.
pub fn load_optional(from: Option<&str>) -> anyhow::Result<(Option<RunManifest>, PathBuf)> {
    let Some(from) = from else {
        return Ok((None, PathBuf::from(".")));
    };
    let Loaded { manifest, dir, .. } = load(from)?;
    Ok((Some(manifest), dir))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_existing_manifest_is_not_replaced_without_asking() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("run.senna.json");
        fs::write(&path, "{}").unwrap();
        let m: RunManifest =
            serde_json::from_str(r#"{"version":2,"kind":"topic","prefix":"run"}"#).unwrap();
        // No terminal under `cargo test`: refused, and the file is untouched.
        assert!(m.save(&path).is_err());
        assert_eq!(fs::read_to_string(&path).unwrap(), "{}");
        // What this process wrote itself it may write again.
        let fresh = dir.path().join("fresh.senna.json");
        m.save(&fresh).unwrap();
        m.save(&fresh).unwrap();
    }

    #[test]
    fn normalize_folds_dots_without_losing_leading_parents() {
        assert_eq!(
            normalize(Path::new("../../a/b")),
            PathBuf::from("../../a/b")
        );
        assert_eq!(normalize(Path::new("a/./b/../c")), PathBuf::from("a/c"));
        assert_eq!(normalize(Path::new("/../a")), PathBuf::from("/a"));
        assert_eq!(normalize(Path::new("a/../../b")), PathBuf::from("../b"));
    }

    #[test]
    fn a_moved_runs_absolute_data_paths_are_found_again() {
        let tmp = tempfile::tempdir().unwrap();
        let proj = tmp.path().join("moved/proj");
        std::fs::create_dir_all(proj.join("run")).unwrap();
        std::fs::create_dir_all(proj.join("data")).unwrap();
        std::fs::write(proj.join("data/a.zarr"), "").unwrap();
        std::fs::create_dir_all(proj.join("place")).unwrap();
        std::fs::write(proj.join("place/b.tsv"), "").unwrap();
        std::fs::write(proj.join("c.tsv"), "").unwrap();
        let m = RunManifest {
            prefix: "/elsewhere/proj/run/x".into(),
            ..serde_json::from_str(r#"{"version":2,"kind":"topic","prefix":""}"#).unwrap()
        };
        let here = proj.join("run").canonicalize().unwrap();
        let root = proj.canonicalize().unwrap();
        // Rebased from the training prefix's directory.
        let a = m.data_file(&proj.join("run"), "/elsewhere/proj/data/a.zarr");
        assert_eq!(Path::new(&a), root.join("data/a.zarr"));
        // Found by its tail (file and folder) under an ancestor.
        let b = m.data_file(&proj.join("run"), "/other/place/b.tsv");
        assert_eq!(Path::new(&b), root.join("place/b.tsv"));
        // A same-named file alone is not the data.
        assert_eq!(
            m.data_file(&here, "/other/where/c.tsv"),
            "/other/where/c.tsv"
        );
        // Unfound: as recorded.
        assert_eq!(m.data_file(&here, "/nowhere/c.tsv"), "/nowhere/c.tsv");
    }

    #[test]
    fn unknown_fields_and_kinds_survive_a_round_trip() {
        let raw = r#"{
            "version": 2,
            "kind": "some-future-kind",
            "prefix": "out/run",
            "train_args": {"senna_version": "9.9", "args": {"k": 7}},
            "outputs": {"cell_embedding": "run.cell_embedding.parquet", "pb_tree": "run.pb_tree.parquet"},
            "annotate": {"argmax": "run.argmax.tsv", "future_slot": 3}
        }"#;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("run.senna.json");
        fs::write(&path, raw).unwrap();

        let (m, _) = RunManifest::load(&path).unwrap();
        assert_eq!(m.kind, RunKind::Other("some-future-kind".into()));
        assert_eq!(m.kind.cell_space(), CellSpace::Signed);
        let saved = dir.path().join("back.senna.json");
        m.save(&saved).unwrap();

        let back: Value = serde_json::from_str(&fs::read_to_string(&saved).unwrap()).unwrap();
        assert_eq!(back["kind"], "some-future-kind");
        assert_eq!(back["train_args"]["args"]["k"], 7);
        assert_eq!(back["outputs"]["pb_tree"], "run.pb_tree.parquet");
        assert_eq!(back["annotate"]["future_slot"], 3);
    }

    #[test]
    fn a_copy_in_a_sibling_directory_resolves_the_same_files() {
        let root = tempfile::tempdir().unwrap();
        let (a, b) = (root.path().join("a"), root.path().join("b"));
        fs::create_dir_all(a.join("layouts")).unwrap();
        fs::create_dir_all(&b).unwrap();
        for f in [
            "counts.zarr",
            "run.latent.parquet",
            "layouts/run.umap.cells.parquet",
        ] {
            fs::write(a.join(f), "").unwrap();
        }
        let raw = r#"{
            "version": 2,
            "kind": "topic",
            "prefix": "run",
            "data": {"input": ["counts.zarr"]},
            "outputs": {"latent": "run.latent.parquet"},
            "layout": {
                "current": "umap",
                "methods": {"umap": {"cell_coords": "layouts/run.umap.cells.parquet"}}
            }
        }"#;
        let file = a.join("run.senna.json");
        fs::write(&file, raw).unwrap();
        let src = load(&file.to_string_lossy()).unwrap();

        let copy = src.copy_to(b.join("out.senna.json")).unwrap();
        copy.manifest.save(&copy.file).unwrap();
        let v: Value = serde_json::from_str(&fs::read_to_string(&copy.file).unwrap()).unwrap();

        let same_file = |from_b: &Value, in_a: &str| {
            let got = Path::new(&resolve(&b, from_b.as_str().unwrap())).canonicalize();
            assert_eq!(
                got.unwrap(),
                a.join(in_a).canonicalize().unwrap(),
                "{from_b}"
            );
        };
        same_file(&v["data"]["input"][0], "counts.zarr");
        same_file(&v["outputs"]["latent"], "run.latent.parquet");
        same_file(
            &v["layout"]["methods"]["umap"]["cell_coords"],
            "layouts/run.umap.cells.parquet",
        );
        same_file(&v["annotate"]["source"], "run.senna.json");
        assert_eq!(v["layout"]["current"], "umap");

        let prefix = PathBuf::from(resolve(&b, v["prefix"].as_str().unwrap()));
        assert_eq!(
            prefix.parent().unwrap().canonicalize().unwrap(),
            a.canonicalize().unwrap()
        );
        assert_eq!(prefix.file_name().unwrap(), "run");

        // The input is untouched.
        assert_eq!(fs::read_to_string(&file).unwrap(), raw);
    }

    #[test]
    fn a_copy_onto_its_own_source_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let file = dir.path().join("run.senna.json");
        fs::write(&file, r#"{"version": 2, "kind": "topic", "prefix": "run"}"#).unwrap();
        let src = load(&file.to_string_lossy()).unwrap();
        assert!(src.copy_to(dir.path().join("run.senna.json")).is_err());
    }

    #[test]
    fn annotating_a_non_senna_run_writes_a_lupin_manifest() {
        let to = |src: &str| annotated_path(Path::new(src), "o/out");
        assert_eq!(to("r/run.senna.json"), PathBuf::from("o/out.senna.json"));
        assert_eq!(to("r/run.pinto.json"), PathBuf::from("o/out.lupin.json"));
        assert_eq!(to("r/round1.lupin.json"), PathBuf::from("o/out.lupin.json"));
    }

    #[test]
    fn known_kinds_keep_their_names() {
        for k in ["topic", "masked-vae", "resolve-embedding-space", "gem"] {
            let kind = RunKind::from(k.to_string());
            assert!(!matches!(kind, RunKind::Other(_)), "{k}");
            assert_eq!(String::from(kind), k);
        }
    }
}
