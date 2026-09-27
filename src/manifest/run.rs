//! Reader/writer for the `{prefix}.senna.json` run manifest a senna fit writes.
//!
//! Lupin models only the fields it reads or records. Every section keeps the
//! rest in an `extra` map, so a manifest round-trips through `load` → `save`
//! without losing anything senna (or a newer senna) put there. A run kind this
//! build does not know still loads, as [`RunKind::Other`].
//!
//! `--from` may name the manifest file itself or the run's output prefix; every
//! command resolves it through [`load`], which also remembers the file so
//! updates are saved back where they were read.

use std::fs;
use std::path::{Path, PathBuf};

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

    /// `{out}.latent.parquet` holds `log θ` on the probability simplex.
    #[must_use]
    pub fn latent_is_log_simplex(&self) -> bool {
        self.cell_space() == CellSpace::LogSimplex
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
    pub pseudotime: RunPseudotime,
    #[serde(default)]
    pub defaults: RunDefaults,
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

    /// The cell table for GEOMETRY (kNN, layout, clustering, trajectory):
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
    /// Projection: nClusters × term FDR q-values.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cluster_term_q: Option<String>,
    /// Projection: per-marker bootstrap support (live markers per type).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub marker_support: Option<String>,
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

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RunPseudotime {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pseudotime: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nodes_latent: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nodes_2d: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub edges: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub root_node: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tree_cell_coords: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tree_nodes_2d: Option<String>,
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
    /// A bare manifest stating only the kind (test fixtures).
    #[cfg(test)]
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
            pseudotime: RunPseudotime::default(),
            defaults: RunDefaults::default(),
            extra: Extra::default(),
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
        if let RunKind::Other(k) = &m.kind {
            log::warn!(
                "manifest {}: run kind `{k}` is unknown to lupin; treating its cells as signed scores",
                path.display()
            );
        }
        let dir = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .map_or_else(|| PathBuf::from("."), Path::to_path_buf);
        Ok((m, dir))
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

    pub fn save(&self, path: &Path) -> anyhow::Result<()> {
        let s = serde_json::to_string_pretty(self)?;
        fs::write(path, s).map_err(|e| anyhow::anyhow!("write {}: {e}", path.display()))?;
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
    let written_abs = abs.canonicalize().unwrap_or(abs);
    match written_abs.strip_prefix(&manifest_abs) {
        Ok(rel) => rel.to_string_lossy().into_owned(),
        Err(_) => written_abs.to_string_lossy().into_owned(),
    }
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

/// The run's output prefix from `--from`: strip `.senna.json` or a trailing `.json`.
#[must_use]
pub fn derive_out_prefix(from: &str) -> String {
    from.strip_suffix(".senna.json")
        .or_else(|| from.strip_suffix(".json"))
        .unwrap_or(from)
        .to_string()
}

/// Default manifest filename for a run `--out` prefix.
#[must_use]
pub fn default_path(prefix: &str) -> String {
    format!("{prefix}.senna.json")
}

/// The manifest file `--from` names: the path itself when it is a file,
/// otherwise `{prefix}.senna.json`.
#[must_use]
pub fn manifest_file(from: &str) -> PathBuf {
    let direct = Path::new(from);
    if direct.is_file() {
        direct.to_path_buf()
    } else {
        PathBuf::from(default_path(&derive_out_prefix(from)))
    }
}

/// A loaded manifest, its directory (for resolving relative paths), and the
/// file it came from (for saving back).
pub struct Loaded {
    pub manifest: RunManifest,
    pub dir: PathBuf,
    pub file: PathBuf,
}

impl Loaded {
    /// `{dir}/{name}` for a manifest at `{dir}/{name}.senna.json`: where the
    /// run's artifacts that the manifest does not record (NB-Fisher weights,
    /// velocity) sit. Derived from where the manifest IS, never from its
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
}

/// Load `--from`, given as a manifest path or a bare prefix.
pub fn load(from: &str) -> anyhow::Result<Loaded> {
    let file = manifest_file(from);
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
        m.save(&path).unwrap();

        let back: Value = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();
        assert_eq!(back["kind"], "some-future-kind");
        assert_eq!(back["train_args"]["args"]["k"], 7);
        assert_eq!(back["outputs"]["pb_tree"], "run.pb_tree.parquet");
        assert_eq!(back["annotate"]["future_slot"], 3);
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
