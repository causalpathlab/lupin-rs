//! What the TUI exported: a log of each figure set (`.svg` + `.pdf` sharing a
//! base name) in `.lupin-view/saved.json` in the directory lupin runs from, in
//! the gallery format `senna view` keeps, plus what lupin needs to trace a
//! figure back and to check it against the files on disk
//! (`docs/trajectory-plan.md` §6). The runs' own directories are not touched.

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

/// Where the log lives, under the working directory.
pub const DIR: &str = ".lupin-view";
const LOG: &str = "saved.json";
/// Saves remembered; older ones go.
const KEEP: usize = 200;

/// One file of an export, as it was written.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct FileRecord {
    pub path: PathBuf,
    pub size: u64,
    /// SHA-256 of the contents, hex.
    pub hash: String,
}

impl FileRecord {
    pub fn of(path: &Path) -> Result<Self> {
        let bytes = std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
        Ok(Self {
            path: path.to_path_buf(),
            size: bytes.len() as u64,
            hash: hex(&Sha256::digest(&bytes)),
        })
    }
}

/// One export: the PDF's path (as senna's gallery keys entries), what and
/// when, then lupin's trace.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Entry {
    /// The PDF, absolute.
    pub path: PathBuf,
    /// What was saved, for the gallery row.
    pub what: String,
    /// Seconds since the Unix epoch.
    pub when: u64,
    /// Kept for senna's gallery format; lupin draws from the SVG instead.
    #[serde(default)]
    pub thumb: String,
    /// The figure (`layout`, `diffusion_dc1_dc2`, …).
    pub panel: String,
    /// The manifest the figure was drawn from.
    pub manifest: PathBuf,
    /// Every file of the set, with size and hash.
    pub files: Vec<FileRecord>,
}

/// What a check finds about an entry.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    Ok,
    Changed,
    Missing,
}

impl Status {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Ok => "ok",
            Self::Changed => "changed since export",
            Self::Missing => "missing",
        }
    }
}

/// The exports logged in one directory, newest first.
pub struct Gallery {
    dir: PathBuf,
    pub entries: Vec<Entry>,
    /// A log that could not be read was moved aside to this file.
    pub set_aside: Option<PathBuf>,
}

impl Gallery {
    /// The log in `dir` (`.lupin-view`). An unreadable log is set aside as
    /// `saved.json.bad`, never treated as empty.
    pub fn open(dir: &Path) -> Self {
        let log = dir.join(LOG);
        let mut set_aside = None;
        let entries = match std::fs::read(&log) {
            Ok(bytes) => match serde_json::from_slice::<Vec<Entry>>(&bytes) {
                Ok(e) => e,
                Err(_) => {
                    let bad = dir.join(format!("{LOG}.bad"));
                    if std::fs::rename(&log, &bad).is_ok() {
                        set_aside = Some(bad);
                    }
                    Vec::new()
                }
            },
            Err(_) => Vec::new(),
        };
        Self {
            dir: dir.to_path_buf(),
            entries,
            set_aside,
        }
    }

    /// The log under the working directory.
    pub fn here() -> Self {
        Self::open(Path::new(DIR))
    }

    fn reload(&mut self) {
        let fresh = Self::open(&self.dir);
        self.entries = fresh.entries;
        if fresh.set_aside.is_some() {
            self.set_aside = fresh.set_aside;
        }
    }

    /// Write the log: to a temporary file, then into place.
    fn save(&self) -> Result<()> {
        std::fs::create_dir_all(&self.dir)?;
        let tmp = self.dir.join(format!("{LOG}.{}.tmp", std::process::id()));
        let text = serde_json::to_string_pretty(&self.entries)?;
        std::fs::write(&tmp, text)?;
        std::fs::rename(&tmp, self.dir.join(LOG))?;
        Ok(())
    }

    /// Log an export of `files` (the PDF among them). An entry for the same
    /// PDF is replaced. The log is reloaded first, so two viewers in one
    /// directory do not overwrite each other's saves.
    pub fn add(
        &mut self,
        files: &[PathBuf],
        what: &str,
        panel: &str,
        manifest: &Path,
    ) -> Result<()> {
        self.reload();
        let pdf = files
            .iter()
            .find(|f| f.extension().is_some_and(|e| e == "pdf"))
            .or(files.first())
            .context("nothing was exported")?;
        let records = files
            .iter()
            .map(|f| FileRecord::of(f))
            .collect::<Result<Vec<_>>>()?;
        let entry = Entry {
            path: absolute(pdf),
            what: what.into(),
            when: SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map_or(0, |d| d.as_secs()),
            thumb: String::new(),
            panel: panel.into(),
            manifest: absolute(manifest),
            files: records,
        };
        self.entries.retain(|e| e.path != entry.path);
        self.entries.insert(0, entry);
        self.entries.truncate(KEEP);
        self.save()
    }

    /// Each entry against the files on disk: the hash is compared only when
    /// size or modification time differ from a quick look, so unchanged
    /// files are not read.
    pub fn check(&self) -> Vec<Status> {
        self.entries
            .iter()
            .map(|e| {
                let mut status = Status::Ok;
                for f in &e.files {
                    match std::fs::metadata(&f.path) {
                        Err(_) => return Status::Missing,
                        Ok(m) => {
                            if m.len() != f.size {
                                status = Status::Changed;
                            } else if let Ok(now) = FileRecord::of(&f.path) {
                                if now.hash != f.hash {
                                    status = Status::Changed;
                                }
                            }
                        }
                    }
                }
                status
            })
            .collect()
    }

    /// Figures beside `manifest_prefix` named `{prefix}.trajectory.*.pdf`
    /// that no entry lists.
    pub fn not_listed(&self, manifest_prefix: &str) -> Vec<PathBuf> {
        let Some(dir) = Path::new(manifest_prefix).parent() else {
            return Vec::new();
        };
        let stem = Path::new(manifest_prefix)
            .file_name()
            .map_or(String::new(), |s| {
                format!("{}.trajectory.", s.to_string_lossy())
            });
        let listed: std::collections::BTreeSet<PathBuf> =
            self.entries.iter().map(|e| e.path.clone()).collect();
        let Ok(read) = std::fs::read_dir(if dir.as_os_str().is_empty() {
            Path::new(".")
        } else {
            dir
        }) else {
            return Vec::new();
        };
        let mut out: Vec<PathBuf> = read
            .filter_map(Result::ok)
            .map(|d| d.path())
            .filter(|p| {
                p.extension().is_some_and(|e| e == "pdf")
                    && p.file_name()
                        .is_some_and(|n| n.to_string_lossy().starts_with(&stem))
                    && !listed.contains(&absolute(p))
            })
            .collect();
        out.sort();
        out
    }

    /// Take an entry off the list, deleting its files when `delete`.
    pub fn remove(&mut self, index: usize, delete: bool) -> Result<()> {
        self.reload();
        if index >= self.entries.len() {
            return Ok(());
        }
        if delete {
            // Only the files that were exported, unchanged, are deleted.
            verify_set(&self.entries[index])?;
        }
        let e = self.entries.remove(index);
        if delete {
            for f in &e.files {
                std::fs::remove_file(&f.path)
                    .with_context(|| format!("deleting {}", f.path.display()))?;
            }
        }
        self.save()
    }

    /// Move an entry's set to a new base name (no extension); nothing is
    /// replaced.
    pub fn relocate(&mut self, index: usize, new_base: &Path) -> Result<()> {
        self.reload();
        let Some(e) = self.entries.get_mut(index) else {
            return Ok(());
        };
        verify_set(e)?;
        let targets: Vec<PathBuf> = e
            .files
            .iter()
            .map(|f| {
                let ext = f
                    .path
                    .extension()
                    .map_or(String::new(), |x| format!(".{}", x.to_string_lossy()));
                PathBuf::from(format!("{}{ext}", new_base.to_string_lossy()))
            })
            .collect();
        for t in &targets {
            anyhow::ensure!(!t.exists(), "{} exists", t.display());
        }
        for (f, t) in e.files.iter_mut().zip(&targets) {
            std::fs::rename(&f.path, t)
                .with_context(|| format!("moving {} to {}", f.path.display(), t.display()))?;
            f.path = absolute(t);
            if t.extension().is_some_and(|x| x == "pdf") {
                e.path = absolute(t);
            }
        }
        self.save()
    }
}

/// Every file of `e` is the export it was logged as: beside the PDF with the
/// same base name, an `svg`, `pdf` or `png`, and unchanged since. The log is
/// a plain file anyone could edit, so nothing is deleted or moved on its word
/// alone.
fn verify_set(e: &Entry) -> Result<()> {
    let base = std::fs::canonicalize(&e.path)
        .with_context(|| format!("{} is not there", e.path.display()))?
        .with_extension("");
    for f in &e.files {
        let c = std::fs::canonicalize(&f.path)
            .with_context(|| format!("{} is not there", f.path.display()))?;
        anyhow::ensure!(
            c.with_extension("") == base
                && matches!(
                    c.extension().and_then(|x| x.to_str()),
                    Some("svg" | "pdf" | "png")
                ),
            "{} is not a file of the export {}",
            c.display(),
            e.path.display()
        );
        anyhow::ensure!(
            FileRecord::of(&c)?.hash == f.hash,
            "{} has changed since it was exported; not touching it",
            c.display()
        );
    }
    Ok(())
}

/// `p` made absolute against the working directory (canonical when it exists).
fn absolute(p: &Path) -> PathBuf {
    std::fs::canonicalize(p).unwrap_or_else(|_| {
        std::env::current_dir()
            .map(|d| d.join(p))
            .unwrap_or_else(|_| p.to_path_buf())
    })
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// "3 min ago", for a gallery row.
pub fn ago(when: u64) -> String {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs());
    let d = now.saturating_sub(when);
    match d {
        0..=59 => "just now".into(),
        60..=3599 => format!("{} min ago", d / 60),
        3600..=86399 => format!("{} h ago", d / 3600),
        _ => format!("{} d ago", d / 86400),
    }
}

#[cfg(test)]
#[path = "tests/gallery.rs"]
mod tests;
