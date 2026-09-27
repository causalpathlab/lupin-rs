//! Read a pinto run (`{prefix}.pinto.json`) as a lupin [`RunManifest`].
//!
//! Pinto stores paths exactly as they were typed at fit time, relative to the
//! directory pinto ran in. Each one is found the way pinto's own viewer finds
//! it: as written if that exists, otherwise its file name beside the manifest.
//! The result is manifest-relative like any other run, so annotation can copy
//! it to a new `{out}.senna.json`. The `.pinto.json` itself is never written.
//!
//! What maps where:
//! - `data_files` → `data.input` (pinto records no batch files)
//! - `outputs.propensity` → `cluster.clusters`: its `cluster` column is the
//!   final level's hard label
//! - `outputs.cell_embedding` / `outputs.feature_embedding` (cage) →
//!   `outputs.cell_embedding` / `outputs.feature_coembedding`: one shared
//!   space scored by dot product, so the run annotates by projection

use super::run::{rel_to_manifest, Loaded, RunKind, RunManifest};
use anyhow::{Context, Result};
use log::{info, warn};
use serde::Deserialize;
use std::path::{Path, PathBuf};

/// The suffix pinto gives its manifest.
pub const SUFFIX: &str = ".pinto.json";

/// The part of pinto's manifest lupin reads; every other field is ignored.
#[derive(Deserialize)]
struct PintoManifest {
    command: String,
    prefix: String,
    #[serde(default)]
    data_files: Option<Vec<String>>,
    #[serde(default)]
    outputs: PintoOutputs,
}

#[derive(Deserialize, Default)]
struct PintoOutputs {
    #[serde(default)]
    propensity: Option<String>,
    #[serde(default)]
    cell_embedding: Option<String>,
    #[serde(default)]
    feature_embedding: Option<String>,
}

/// Whether `path` names a pinto manifest.
#[must_use]
pub fn is_pinto(path: &Path) -> bool {
    path.to_string_lossy().ends_with(SUFFIX)
}

/// Load `file` and translate it into a run manifest whose paths are relative
/// to `file`'s directory.
pub fn load(file: &Path) -> Result<Loaded> {
    let raw = std::fs::read_to_string(file).with_context(|| format!("read {}", file.display()))?;
    let p: PintoManifest =
        serde_json::from_str(&raw).with_context(|| format!("parse {}", file.display()))?;
    let dir = file
        .parent()
        .filter(|d| !d.as_os_str().is_empty())
        .map_or_else(|| PathBuf::from("."), Path::to_path_buf);

    let find = |written: &str| -> Option<String> {
        let found = locate(written, &dir);
        if found.is_none() {
            warn!("{}: `{written}` not found", file.display());
        }
        found.map(|abs| rel_to_manifest(&dir, &abs.to_string_lossy()))
    };
    let slot = |o: &Option<String>| o.as_deref().and_then(find);

    // The run's stem is where its manifest is (`{dir}/{name}.pinto.json`),
    // not the `prefix` pinto recorded, which is relative to where it ran.
    let name = file.file_name().map(|n| n.to_string_lossy().into_owned());
    let stem = name
        .as_deref()
        .and_then(|n| n.strip_suffix(SUFFIX))
        .unwrap_or(&p.prefix);
    let mut m = RunManifest::new(RunKind::Other(format!("pinto-{}", p.command)), stem);
    m.data.input = p
        .data_files
        .iter()
        .flatten()
        .filter_map(|f| find(f))
        .collect();
    m.cluster.clusters = slot(&p.outputs.propensity);
    m.outputs.cell_embedding = slot(&p.outputs.cell_embedding);
    m.outputs.feature_coembedding = slot(&p.outputs.feature_embedding);

    info!(
        "Loaded pinto run {} (command: {})",
        file.display(),
        p.command
    );
    Ok(Loaded {
        manifest: m,
        dir,
        file: file.to_path_buf(),
    })
}

/// A path pinto wrote: as written (from the working directory) if it exists,
/// else its file name in `dir`.
fn locate(written: &str, dir: &Path) -> Option<PathBuf> {
    let direct = PathBuf::from(written);
    if direct.exists() {
        return Some(direct);
    }
    let beside = dir.join(direct.file_name()?);
    beside.exists().then_some(beside)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::manifest::run::resolve;
    use std::fs;

    #[test]
    fn a_cage_run_maps_onto_projection_inputs() {
        let root = tempfile::tempdir().unwrap();
        let dir = root.path().join("fit");
        fs::create_dir_all(&dir).unwrap();
        for f in [
            "counts.zarr",
            "run.propensity.parquet",
            "run.cell_embedding.parquet",
            "run.feature_embedding.parquet",
        ] {
            fs::write(dir.join(f), "").unwrap();
        }
        // Paths as pinto writes them: relative to a working directory that
        // is not the manifest's.
        let raw = r#"{
            "command": "cage", "version": "0.1", "timestamp": "t",
            "prefix": "elsewhere/run", "n_cells": 3, "n_features": 2,
            "data_files": ["elsewhere/counts.zarr"],
            "outputs": {
                "propensity": "elsewhere/run.propensity.parquet",
                "cell_embedding": "elsewhere/run.cell_embedding.parquet",
                "feature_embedding": "elsewhere/run.feature_embedding.parquet",
                "cell_bias": "elsewhere/run.cell_bias.parquet"
            }
        }"#;
        let file = dir.join("run.pinto.json");
        fs::write(&file, raw).unwrap();

        let loaded = load(&file).unwrap();
        let m = &loaded.manifest;
        assert_eq!(m.kind, RunKind::Other("pinto-cage".into()));
        let at = |rel: &str| {
            Path::new(&resolve(&loaded.dir, rel))
                .canonicalize()
                .unwrap()
        };
        let want = |f: &str| dir.join(f).canonicalize().unwrap();
        assert_eq!(at(&m.data.input[0]), want("counts.zarr"));
        assert_eq!(
            at(m.cluster.clusters.as_deref().unwrap()),
            want("run.propensity.parquet")
        );
        assert_eq!(
            at(m.outputs.cell_embedding.as_deref().unwrap()),
            want("run.cell_embedding.parquet")
        );
        assert_eq!(
            at(m.outputs.feature_coembedding.as_deref().unwrap()),
            want("run.feature_embedding.parquet")
        );
        assert_eq!(m.prefix, "run");

        // Annotation's copy, in another directory, still finds every table
        // and names the pinto manifest it came from.
        fs::create_dir_all(root.path().join("out")).unwrap();
        let copy = loaded
            .copy_to(root.path().join("out/o.senna.json"))
            .unwrap();
        copy.manifest.save(&copy.file).unwrap();
        let (back, out_dir) = RunManifest::load(&copy.file).unwrap();
        let from_out = |rel: &str| Path::new(&resolve(&out_dir, rel)).canonicalize().unwrap();
        assert_eq!(
            from_out(back.cluster.clusters.as_deref().unwrap()),
            want("run.propensity.parquet")
        );
        assert_eq!(
            from_out(back.annotate.source.as_deref().unwrap()),
            file.canonicalize().unwrap()
        );

        // The source is never touched.
        assert_eq!(fs::read_to_string(&file).unwrap(), raw);
    }
}
