//! What the TUI exported: a log of each figure set (`.svg` + `.pdf` sharing a
//! base name) in `.lupin-view/saved.json` in the directory lupin runs from,
//! with what lupin needs to trace a figure back and to check it against the
//! files on disk (`docs/trajectory-plan.md` §6). The runs' own directories
//! are not touched.

use crate::manifest::data_files::absolute;
use anyhow::{ensure, Context, Result};
use ratatui::style::Color;
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
    /// Modification time, seconds since the Unix epoch; a quick check
    /// before the hash.
    #[serde(default)]
    pub mtime: u64,
    /// SHA-256 of the contents, hex.
    pub hash: String,
}

impl FileRecord {
    pub fn of(path: &Path) -> Result<Self> {
        let bytes = std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
        Ok(Self {
            path: absolute(path),
            size: bytes.len() as u64,
            mtime: std::fs::metadata(path).map_or(0, |m| mtime_of(&m)),
            hash: hex(&Sha256::digest(&bytes)),
        })
    }

    /// The file on disk is still this record: same size and mtime, else the
    /// same hash.
    fn unchanged(&self) -> Option<bool> {
        let m = std::fs::metadata(&self.path).ok()?;
        if m.len() != self.size {
            return Some(false);
        }
        if mtime_of(&m) == self.mtime && self.mtime != 0 {
            return Some(true);
        }
        Some(Self::of(&self.path).is_ok_and(|now| now.hash == self.hash))
    }
}

fn mtime_of(m: &std::fs::Metadata) -> u64 {
    m.modified()
        .ok()
        .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
        .map_or(0, |d| d.as_secs())
}

/// One export: the PDF's path, what and when, and its trace.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Entry {
    /// The PDF, absolute.
    pub path: PathBuf,
    /// What was saved, for the gallery row.
    pub what: String,
    /// Seconds since the Unix epoch.
    pub when: u64,
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
    Intact,
    Changed,
    Missing,
}

impl Status {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Intact => "ok",
            Self::Changed => "changed since export",
            Self::Missing => "missing",
        }
    }

    /// The row's colour; `None` for the default.
    pub fn color(self) -> Option<Color> {
        match self {
            Self::Intact => None,
            Self::Changed => Some(Color::Yellow),
            Self::Missing => Some(Color::Red),
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
    /// Write the log: to a fresh temporary file (created exclusively, under
    /// a random name, readable by this user only), then into place.
    fn save(&self) -> Result<()> {
        // Created private; an existing directory is left as it is.
        let mut make = std::fs::DirBuilder::new();
        make.recursive(true);
        #[cfg(unix)]
        std::os::unix::fs::DirBuilderExt::mode(&mut make, 0o700);
        make.create(&self.dir)?;
        let mut tmp = tempfile::NamedTempFile::new_in(&self.dir)?;
        serde_json::to_writer_pretty(&mut tmp, &self.entries)?;
        tmp.persist(self.dir.join(LOG))?;
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
        // A log this user does not own alone is not taken over by saving
        // over it: its entries would then be trusted for `D` and `m`.
        if std::fs::symlink_metadata(self.dir.join(LOG)).is_ok() {
            self.entries = read_private(&self.dir)?;
        }
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
            when: now(),
            panel: panel.into(),
            manifest: absolute(manifest),
            files: records,
        };
        self.entries.retain(|e| e.path != entry.path);
        self.entries.insert(0, entry);
        self.entries.truncate(KEEP);
        self.save()
    }

    /// Each entry against the files on disk: a file whose size and mtime are
    /// as recorded is taken as unchanged without reading it.
    pub fn check(&self) -> Vec<Status> {
        self.entries
            .iter()
            .map(|e| {
                let mut status = Status::Intact;
                for f in &e.files {
                    match f.unchanged() {
                        None => return Status::Missing,
                        Some(false) => status = Status::Changed,
                        Some(true) => {}
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

    /// The entry for the PDF `pdf`, from the log as it is now on disk
    /// (another viewer may have changed it since the row on screen was
    /// drawn), read only from a log this user alone can write: files are
    /// deleted or moved on its word.
    fn find(&mut self, pdf: &Path) -> Result<usize> {
        self.entries = read_private(&self.dir)?;
        self.entries
            .iter()
            .position(|e| e.path == pdf)
            .with_context(|| format!("{} is no longer listed", pdf.display()))
    }

    /// Take the export of `pdf` off the list, deleting its files when
    /// `delete`.
    pub fn remove(&mut self, pdf: &Path, delete: bool) -> Result<()> {
        let index = self.find(pdf)?;
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

    /// Move the export of `pdf` to a new base name (no extension); nothing is
    /// replaced, and a failed move puts back the files already moved.
    pub fn relocate(&mut self, pdf: &Path, new_base: &Path) -> Result<()> {
        let index = self.find(pdf)?;
        let e = &mut self.entries[index];
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
            ensure!(!t.exists(), "{} exists", t.display());
        }
        for (k, (f, t)) in e.files.iter().zip(&targets).enumerate() {
            if let Err(err) = std::fs::rename(&f.path, t) {
                for (g, u) in e.files.iter().zip(&targets).take(k) {
                    let _ = std::fs::rename(u, &g.path);
                }
                return Err(err)
                    .with_context(|| format!("moving {} to {}", f.path.display(), t.display()));
            }
        }
        for (f, t) in e.files.iter_mut().zip(&targets) {
            f.path = absolute(t);
            if t.extension().is_some_and(|x| x == "pdf") {
                e.path = absolute(t);
            }
        }
        self.save()
    }
}

/// The entries of the log in `dir`, only when the directory and the log are
/// this user's, no one else can write them, and neither is a link: a log
/// another user can edit, or swap for another, could name any of this
/// user's files, and files are deleted or moved on its word. The log is
/// checked on the file that was opened, and that file must be the one at
/// the log's path.
fn read_private(dir: &Path) -> Result<Vec<Entry>> {
    use std::io::Read;
    let log = dir.join(LOG);
    let mut f = std::fs::File::open(&log).with_context(|| format!("opening {}", log.display()))?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        // This process's user: the owner of a file it has just created.
        let me = tempfile::tempfile()?.metadata()?.uid();
        let mine = |m: &std::fs::Metadata| m.uid() == me && m.mode() & 0o022 == 0;
        let refuse = |p: &Path| {
            anyhow::anyhow!(
                "{} is not this user's alone (a link, or others can write it): not deleting \
                 or moving files on its word (chmod go-w, or remove it to start a new log)",
                p.display()
            )
        };
        let d = std::fs::symlink_metadata(dir)?;
        if d.file_type().is_symlink() || !d.is_dir() || !mine(&d) {
            return Err(refuse(dir));
        }
        let (at, opened) = (std::fs::symlink_metadata(&log)?, f.metadata()?);
        let same = at.dev() == opened.dev() && at.ino() == opened.ino();
        if at.file_type().is_symlink() || !opened.is_file() || !same || !mine(&opened) {
            return Err(refuse(&log));
        }
    }
    let mut bytes = Vec::new();
    f.read_to_end(&mut bytes)?;
    serde_json::from_slice(&bytes).with_context(|| format!("reading {}", log.display()))
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
        ensure!(
            c.with_extension("") == base
                && matches!(
                    c.extension().and_then(|x| x.to_str()),
                    Some("svg" | "pdf" | "png")
                ),
            "{} is not a file of the export {}",
            c.display(),
            e.path.display()
        );
        ensure!(
            FileRecord::of(&c)?.hash == f.hash,
            "{} has changed since it was exported; not touching it",
            c.display()
        );
    }
    Ok(())
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// Seconds since the Unix epoch.
pub fn now() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs())
}

/// "3 min ago", for a gallery row.
pub fn ago(when: u64, now: u64) -> String {
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
