//! A run's family: the run itself, the annotation rounds made from it
//! (`X.L0`, `X.L0.L1`, `X.r2`, …) and the trajectories made from those
//! (`X.T1`, …), found beside the manifest opened.
//!
//! The names come from one table of tags ([`Tag`]): the producers format
//! them from it and [`stem`] parses them with it, so the two cannot drift.
//! Members are the manifests (senna, lupin or pinto) in the opened
//! manifest's directory whose name grows the family's stem (`X.` …), and
//! any other manifest there whose `annotate.source` names a member (a
//! trajectory written under another prefix). Paths keep the spelling of the
//! directory listed, so a run opened through a link shows as it was opened.

use super::run::{
    annotated_path, derive_out_prefix, load, parent_dir, resolve, same_file, RunManifest,
    LUPIN_SUFFIX,
};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::path::{Path, PathBuf};

/// Every manifest suffix a member can have.
const SUFFIXES: [&str; 3] = [".senna.json", LUPIN_SUFFIX, super::pinto::SUFFIX];

/// A tag a derived manifest's name adds to the name it was made from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tag {
    /// An annotation round written by a pass (`X.L1`).
    Round,
    /// A trajectory's output (`X.T1`).
    Trajectory,
    /// A round of a relabel chain (`X.r2`).
    Chain,
}

impl Tag {
    const ALL: [Tag; 3] = [Tag::Round, Tag::Trajectory, Tag::Chain];

    fn letter(self) -> char {
        match self {
            Tag::Round => 'L',
            Tag::Trajectory => 'T',
            Tag::Chain => 'r',
        }
    }

    /// `{prefix}.{letter}{k}`.
    pub fn name(self, prefix: &str, k: usize) -> String {
        format!("{prefix}.{}{k}", self.letter())
    }

    /// `(prefix, k)` when `name` ends in this tag.
    pub fn parse(self, name: &str) -> Option<(&str, usize)> {
        let (prefix, k) = name.rsplit_once(&format!(".{}", self.letter()))?;
        (!k.is_empty() && k.bytes().all(|b| b.is_ascii_digit()))
            .then(|| k.parse().ok().map(|k| (prefix, k)))
            .flatten()
    }

    /// The first `{stem}.{tag}{k}` (k = 1, 2, …) beside `source` with no
    /// manifest yet, `stem` being `source`'s prefix.
    pub fn next_free(self, source: &Path) -> String {
        let stem = derive_out_prefix(&source.to_string_lossy());
        (1..)
            .map(|k| self.name(&stem, k))
            .find(|o| !annotated_path(source, o).exists())
            .expect("some k is free")
    }
}

/// Where passes from the manifest shown at `path` start, and the output a
/// new pass is offered: the next free round under it, never an existing one.
pub fn pass_origin(path: &Path) -> (PathBuf, String) {
    (path.to_path_buf(), Tag::Round.next_free(path))
}

/// What a member of the family is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// A run with no annotation (the senna or pinto run itself).
    Run,
    /// An annotation round: clusters with labels.
    Round,
    /// A trajectory's manifest copy.
    Trajectory,
}

impl Kind {
    pub fn of(m: &RunManifest) -> Self {
        let t = &m.trajectory;
        if t.prior.is_some() || t.pseudotime.is_some() {
            Kind::Trajectory
        } else if m.annotate.argmax.is_some() {
            Kind::Round
        } else {
            Kind::Run
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Kind::Run => "run",
            Kind::Round => "round",
            Kind::Trajectory => "trajectory",
        }
    }

    /// What showing a member of this kind does, in a line.
    /// What showing it does; `figures` when the TUI shows trajectories'
    /// figures (`lupin trajectory`'s).
    pub fn detail(self, figures: bool) -> &'static str {
        match (self, figures) {
            (Kind::Run, _) => "the run itself: no clusters on the left until a pass (A); passes start from it",
            (Kind::Round, true) => "show this round's clusters and labels; the order view takes its labels, and passes start from it",
            (Kind::Round, false) => "show this round's clusters and labels; passes start from it",
            (Kind::Trajectory, true) => "show this trajectory's figures with the round it carries; passes start from it",
            (Kind::Trajectory, false) => "show the round this trajectory carries (its figures are lupin trajectory's); passes start from it",
        }
    }
}

/// A member to show: what [`Member`] knows that showing it needs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Pick {
    pub path: PathBuf,
    pub kind: Kind,
    /// It carries an annotation (a round, or a trajectory's copy of one).
    pub labelled: bool,
}

pub struct Member {
    pub pick: Pick,
    /// What it holds, in a few words.
    pub holds: String,
    /// Seconds since the epoch of its last change.
    pub modified: u64,
}

impl Member {
    /// The file name without the manifest suffix.
    pub fn name(&self) -> String {
        derive_out_prefix(&file_name(&self.pick.path))
    }
}

/// A path's file name, or the whole path when it has none.
pub fn file_name(p: &Path) -> String {
    p.file_name().map_or_else(
        || p.display().to_string(),
        |n| n.to_string_lossy().into_owned(),
    )
}

/// The family's stem: `X` for `X`, `X.L0`, `X.L0.L1`, `X.T2`, `X.r3`.
pub fn stem(manifest: &Path) -> String {
    let mut s = derive_out_prefix(&file_name(manifest));
    while let Some(head) = Tag::ALL.iter().find_map(|t| t.parse(&s).map(|(h, _)| h)) {
        s = head.to_string();
    }
    s
}

/// The manifests in `opened`'s directory, by file name.
fn manifests_beside(opened: &Path) -> (PathBuf, Vec<String>) {
    let dir = parent_dir(opened);
    let names = std::fs::read_dir(&dir)
        .map(|entries| {
            entries
                .flatten()
                .map(|e| e.file_name().to_string_lossy().into_owned())
                .filter(|n| !n.starts_with('.') && SUFFIXES.iter().any(|s| n.ends_with(s)))
                // A pass not yet promoted onto its round (`staging`).
                .filter(|n| !n.contains(".staging."))
                .collect()
        })
        .unwrap_or_default();
    (dir, names)
}

/// Whether manifest `name` grows `stem`.
fn grows(name: &str, stem: &str) -> bool {
    let s = derive_out_prefix(name);
    s == stem || s.strip_prefix(stem).is_some_and(|r| r.starts_with('.'))
}

/// The members of `opened`'s family, the run first, then by name.
pub fn family(opened: &Path) -> Vec<Member> {
    let (dir, names) = manifests_beside(opened);
    let stem = stem(opened);
    // Every manifest beside it, read once.
    let mut all: BTreeMap<String, RunManifest> = BTreeMap::new();
    for name in names {
        if let Ok(l) = load(&dir.join(&name).to_string_lossy()) {
            all.insert(name, l.manifest);
        }
    }
    let mut kin: HashSet<String> = all.keys().filter(|n| grows(n, &stem)).cloned().collect();
    // Manifests made from a member under another prefix join it.
    loop {
        let more: Vec<String> = all
            .iter()
            .filter(|(n, m)| !kin.contains(*n) && source_name(m).is_some_and(|s| kin.contains(&s)))
            .map(|(n, _)| n.clone())
            .collect();
        if more.is_empty() {
            break;
        }
        kin.extend(more);
    }
    let mut out: Vec<Member> = all
        .into_iter()
        .filter(|(n, _)| kin.contains(n))
        .map(|(n, m)| member(&dir, dir.join(n), &m))
        .collect();
    out.sort_by(|a, b| {
        (a.pick.kind != Kind::Run)
            .cmp(&(b.pick.kind != Kind::Run))
            .then_with(|| a.name().cmp(&b.name()))
    });
    out
}

/// Whether `opened`'s family has an annotation round: the first family
/// manifest with labels answers, with nothing else read.
pub fn has_round(opened: &Path) -> bool {
    let (dir, names) = manifests_beside(opened);
    let stem = stem(opened);
    names.iter().filter(|n| grows(n, &stem)).any(|n| {
        load(&dir.join(n).to_string_lossy()).is_ok_and(|l| Kind::of(&l.manifest) == Kind::Round)
    })
}

/// The file name `annotate.source` names.
fn source_name(m: &RunManifest) -> Option<String> {
    let s = m.annotate.source.as_deref()?;
    Path::new(s)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
}

/// A cluster summary's labels, and nothing else of it.
#[derive(serde::Deserialize)]
struct Labelled {
    label: Option<String>,
}

fn member(dir: &Path, path: PathBuf, m: &RunManifest) -> Member {
    let kind = Kind::of(m);
    let clusters = || {
        let rel = m.annotate.cluster_summary.as_deref()?;
        let text = std::fs::read_to_string(resolve(dir, rel)).ok()?;
        let d: BTreeMap<String, Labelled> = serde_json::from_str(&text).ok()?;
        let labels: BTreeSet<&str> = d.values().filter_map(|x| x.label.as_deref()).collect();
        Some(format!("{} clusters, {} labels", d.len(), labels.len()))
    };
    let holds = match kind {
        Kind::Run => "no annotation".into(),
        Kind::Round => clusters().unwrap_or_else(|| "annotated".into()),
        Kind::Trajectory => {
            let roots = m
                .trajectory
                .settings
                .as_ref()
                .and_then(|s| s.get("roots"))
                .and_then(|r| r.as_array())
                .map(|r| r.iter().filter_map(|v| v.as_str()).collect::<Vec<_>>())
                .unwrap_or_default();
            let mut s = "trajectory".to_string();
            if let Some(f) = source_name(m) {
                s += &format!(" on {}", derive_out_prefix(&f));
            }
            if !roots.is_empty() {
                s += &format!(", from {}", roots.join(", "));
            }
            s
        }
    };
    let modified = std::fs::metadata(&path).map_or(0, |md| mtime_of(&md));
    Member {
        pick: Pick {
            labelled: m.annotate.argmax.is_some(),
            path,
            kind,
        },
        holds,
        modified,
    }
}

/// A file's last change, in seconds since the epoch (0 when unknown).
pub fn mtime_of(m: &std::fs::Metadata) -> u64 {
    m.modified()
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map_or(0, |d| d.as_secs())
}

/// The member that is `shown`, compared as files (a run opened through a
/// link or another spelling of its directory is still found).
pub fn current(members: &[Member], shown: &Path) -> Option<usize> {
    members.iter().position(|m| same_file(&m.pick.path, shown))
}

#[cfg(test)]
#[path = "family_tests.rs"]
pub(crate) mod tests;
