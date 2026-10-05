//! Marker panels bundled as data in the `lupin-panels` crate, named on the
//! command line as `panel:<name>`. A bundled panel is written into lupin's
//! cache as ordinary files (the table and its Cell Ontology sidecar), so the
//! rest of lupin reads it as it reads any panel, and the run records the
//! panel's name and the crate's version.

use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// How a bundled panel is named on the command line.
pub const PREFIX: &str = "panel:";

/// A bundled panel, as a run records it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Bundled {
    pub name: String,
    /// The `lupin-panels` version the panel came from.
    pub version: String,
}

/// The bundled panels' names; none when lupin is built without them.
#[must_use]
pub fn names() -> Vec<&'static str> {
    #[cfg(feature = "panels")]
    {
        lupin_panels::names().to_vec()
    }
    #[cfg(not(feature = "panels"))]
    {
        Vec::new()
    }
}

/// `markers` as a file path: a `panel:<name>` is written out under
/// [`panels_dir`] first, with the record of where it came from.
pub fn resolve(markers: &str) -> Result<(String, Option<Bundled>)> {
    match markers.strip_prefix(PREFIX) {
        Some(name) => {
            let (path, bundled) = write_out(name, &panels_dir())?;
            Ok((path.to_string_lossy().into_owned(), Some(bundled)))
        }
        None => Ok((markers.to_string(), None)),
    }
}

/// Where bundled panels are written: lupin's cache, else the temporary
/// directory.
#[must_use]
pub fn panels_dir() -> PathBuf {
    if cfg!(test) {
        // Tests write nothing into the machine's cache.
        return std::env::temp_dir().join("lupin-test").join("panels");
    }
    dirs::cache_dir()
        .map(|d| d.join("lupin"))
        .unwrap_or_else(|| std::env::temp_dir().join("lupin"))
        .join("panels")
}

/// Write bundled panel `name` under `root/<version>/<name>/` (its table and
/// sidecar, left as they are when already written); the table's path.
pub fn write_out(name: &str, root: &Path) -> Result<(PathBuf, Bundled)> {
    #[cfg(feature = "panels")]
    {
        use anyhow::Context;
        let Some(p) = lupin_panels::get(name) else {
            anyhow::bail!("no bundled panel {name}; bundled: {}", names().join(", "));
        };
        let dir = root.join(lupin_panels::VERSION).join(name);
        std::fs::create_dir_all(&dir).with_context(|| format!("creating {}", dir.display()))?;
        let table = dir.join(p.table_file);
        put(&table, p.table)?;
        if let (Some(bytes), Some(side)) = (p.sidecar, super::data_files::sidecar_of(&table)) {
            put(&side, bytes)?;
        }
        put(&dir.join("README.md"), p.readme)?;
        Ok((
            table,
            Bundled {
                name: name.to_string(),
                version: lupin_panels::VERSION.to_string(),
            },
        ))
    }
    #[cfg(not(feature = "panels"))]
    {
        let _ = root;
        anyhow::bail!(
            "{PREFIX}{name}: this lupin was built without bundled panels (feature `panels`)"
        )
    }
}

/// The bundled panel a table at `path` was written from by [`write_out`]
/// (`…/panels/<version>/<name>/<table>`), else `None`.
#[must_use]
pub fn bundled_of(path: &Path) -> Option<Bundled> {
    let dir = path.parent()?;
    let name = dir.file_name()?.to_str()?;
    let version_dir = dir.parent()?;
    let version = version_dir.file_name()?.to_str()?;
    (version_dir.parent()?.file_name()? == "panels" && names().contains(&name)).then(|| Bundled {
        name: name.to_string(),
        version: version.to_string(),
    })
}

/// Write `bytes` to `path` unless it already holds them, through a side file
/// renamed once whole.
#[cfg(feature = "panels")]
fn put(path: &Path, bytes: &[u8]) -> Result<()> {
    use anyhow::Context;
    use std::fs;
    if fs::read(path).is_ok_and(|b| b == bytes) {
        return Ok(());
    }
    let tmp = path.with_extension("part");
    fs::write(&tmp, bytes).with_context(|| format!("writing {}", tmp.display()))?;
    fs::rename(&tmp, path)?;
    Ok(())
}

#[cfg(all(test, feature = "panels"))]
#[path = "panels_tests.rs"]
mod tests;
