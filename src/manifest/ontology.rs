//! The Cell Ontology for annotation, found without being asked for: an
//! explicit `--obo`, else the copy lupin cached, else a fresh download into
//! that cache. With none of those (offline, say), the panel's coarse level is
//! inferred from the marker genes its types share instead.

use crate::annotate::celltype_tree::{ClTerms, TypeTree};
use anyhow::{Context, Result};
use log::{info, warn};
use std::fs;
use std::path::PathBuf;
use std::time::Duration;

/// Where the Cell Ontology's basic edition is published.
const CL_URL: &str = "https://purl.obolibrary.org/obo/cl/cl-basic.obo";
/// Set to skip the download (and use only `--obo` or the cache).
pub const OFFLINE_ENV: &str = "LUPIN_OFFLINE";

/// `{user cache}/lupin/cl-basic.obo`.
fn cache_path() -> Option<PathBuf> {
    dirs::cache_dir().map(|d| d.join("lupin").join("cl-basic.obo"))
}

/// The ontology file to use: `explicit` when given; else the cached copy;
/// else, unless offline, a download into the cache. `None` when there is
/// none to be had, which is not an error.
#[must_use]
pub fn resolve_obo(explicit: Option<&str>) -> Option<PathBuf> {
    if let Some(p) = explicit {
        return Some(PathBuf::from(p));
    }
    let cached = cache_path()?;
    if cached.is_file() {
        return Some(cached);
    }
    if std::env::var_os(OFFLINE_ENV).is_some() {
        info!("{OFFLINE_ENV} is set: not downloading the Cell Ontology");
        return None;
    }
    match download(&cached) {
        Ok(()) => Some(cached),
        Err(e) => {
            warn!("could not download the Cell Ontology ({e:#}); grouping cell types by shared markers instead");
            None
        }
    }
}

fn download(to: &std::path::Path) -> Result<()> {
    info!("downloading the Cell Ontology from {CL_URL}");
    let response = ureq::AgentBuilder::new()
        .timeout(Duration::from_secs(60))
        .build()
        .get(CL_URL)
        .call()
        .context("request failed")?;
    let dir = to.parent().context("no cache directory")?;
    fs::create_dir_all(dir)?;
    // Streamed to a side file and renamed once whole, so a half-finished
    // download is never read.
    let tmp = to.with_extension("obo.part");
    let mut file = fs::File::create(&tmp)?;
    std::io::copy(&mut response.into_reader(), &mut file).context("reading the response")?;
    drop(file);
    let head = fs::read_to_string(&tmp).unwrap_or_default();
    anyhow::ensure!(head.contains("[Term]"), "the response is not an OBO file");
    fs::rename(&tmp, to)?;
    info!("cached {}", to.display());
    Ok(())
}

/// The coarse level over `panel`: on the Cell Ontology when one is found and
/// at least two of the panel's types match its terms, else by shared markers.
pub fn type_tree(
    explicit_obo: Option<&str>,
    panel: &[(String, String)],
) -> Result<(TypeTree, Option<PathBuf>)> {
    let obo = resolve_obo(explicit_obo);
    if let Some(path) = &obo {
        let text =
            fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        let terms = ClTerms::parse(&text);
        if let Some(tree) = TypeTree::from_ontology(&terms, panel) {
            info!(
                "cell-type groups from the Cell Ontology ({}): {} group(s), {} type(s) matched",
                terms.release.as_deref().unwrap_or("release unknown"),
                tree.groups.len(),
                tree.label_cl.len()
            );
            return Ok((tree, obo));
        }
        warn!("fewer than two panel types match a Cell Ontology term; grouping by shared markers");
    }
    let tree = TypeTree::from_markers(panel);
    info!(
        "cell-type groups by shared markers: {} group(s)",
        tree.groups.len()
    );
    Ok((tree, obo))
}
