//! A run's family: the run itself, the annotation rounds made from it
//! (`X.L0`, `X.L0.L1`, `X.r2`, …) and the trajectories made from those
//! (`X.T1`, …), found beside the manifest opened so the TUI can move
//! between them (`g`).
//!
//! Members are the manifests in the opened manifest's directory whose name
//! grows the family's stem (`X.` …), and any other manifest there whose
//! `annotate.source` names a member (a trajectory written under another
//! prefix). Paths keep the spelling of the directory listed, so a run
//! opened through a link shows as it was opened.

use crate::manifest::run::{annotated_path, derive_out_prefix, resolve, same_file, RunManifest};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// What a member of the family is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// A run with no annotation (the senna run itself).
    Run,
    /// An annotation round: clusters with labels.
    Round,
    /// A trajectory's manifest copy.
    Trajectory,
}

impl Kind {
    pub fn name(self) -> &'static str {
        match self {
            Kind::Run => "run",
            Kind::Round => "round",
            Kind::Trajectory => "trajectory",
        }
    }
}

pub struct Member {
    pub path: PathBuf,
    pub kind: Kind,
    /// What it holds, in a few words.
    pub holds: String,
    /// Seconds since the epoch of its last change.
    pub modified: u64,
}

impl Member {
    /// The file name without the manifest suffix.
    pub fn name(&self) -> String {
        derive_out_prefix(&file_name(&self.path))
    }
}

fn file_name(p: &Path) -> String {
    p.file_name().map_or_else(
        || p.to_string_lossy().into_owned(),
        |n| n.to_string_lossy().into_owned(),
    )
}

/// The family's stem: `X` for `X`, `X.L0`, `X.L0.L1`, `X.T2`, `X.r3`.
pub fn stem(manifest: &Path) -> String {
    let mut s = derive_out_prefix(&file_name(manifest));
    while let Some((head, last)) = s.rsplit_once('.') {
        let tagged = last.len() > 1
            && matches!(last.as_bytes()[0], b'L' | b'T' | b'r')
            && last[1..].bytes().all(|b| b.is_ascii_digit());
        if !tagged {
            break;
        }
        s = head.to_string();
    }
    s
}

/// The members of `opened`'s family, the run first, then by name.
pub fn family(opened: &Path) -> Vec<Member> {
    let dir = match opened.parent() {
        Some(d) if !d.as_os_str().is_empty() => d.to_path_buf(),
        _ => PathBuf::from("."),
    };
    let suffix = annotated_path(opened, "").to_string_lossy().into_owned();
    let stem = stem(opened);
    let Ok(entries) = std::fs::read_dir(&dir) else {
        return Vec::new();
    };
    // Every manifest beside it, read once.
    let mut all: BTreeMap<String, (PathBuf, RunManifest)> = BTreeMap::new();
    for e in entries.flatten() {
        let name = e.file_name().to_string_lossy().into_owned();
        if !name.ends_with(&suffix) || name.starts_with('.') {
            continue;
        }
        let path = dir.join(&name);
        if let Ok((m, _)) = RunManifest::load(&path) {
            all.insert(name, (path, m));
        }
    }
    let grows = |name: &str| {
        let s = name.strip_suffix(&suffix).unwrap_or(name);
        s == stem || s.starts_with(&format!("{stem}."))
    };
    let mut kin: Vec<String> = all.keys().filter(|n| grows(n)).cloned().collect();
    // Manifests made from a member under another prefix join it.
    loop {
        let more: Vec<String> = all
            .iter()
            .filter(|(n, _)| !kin.contains(n))
            .filter(|(_, (_, m))| source_name(m).is_some_and(|s| kin.contains(&s)))
            .map(|(n, _)| n.clone())
            .collect();
        if more.is_empty() {
            break;
        }
        kin.extend(more);
    }
    let mut out: Vec<Member> = kin
        .into_iter()
        .filter_map(|n| all.remove(&n))
        .map(|(path, m)| member(&dir, path, &m))
        .collect();
    out.sort_by(|a, b| {
        (a.kind != Kind::Run)
            .cmp(&(b.kind != Kind::Run))
            .then_with(|| a.name().cmp(&b.name()))
    });
    out
}

/// The file name `annotate.source` names.
fn source_name(m: &RunManifest) -> Option<String> {
    let s = m.annotate.source.as_deref()?;
    Path::new(s)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
}

fn member(dir: &Path, path: PathBuf, m: &RunManifest) -> Member {
    let t = &m.trajectory;
    let kind = if t.prior.is_some() || t.pseudotime.is_some() {
        Kind::Trajectory
    } else if m.annotate.argmax.is_some() {
        Kind::Round
    } else {
        Kind::Run
    };
    let clusters = || {
        let rel = m.annotate.cluster_summary.as_deref()?;
        let text = std::fs::read_to_string(resolve(dir, rel)).ok()?;
        let d: BTreeMap<String, crate::annotate::rounds::Digest> =
            serde_json::from_str(&text).ok()?;
        let labels: std::collections::BTreeSet<&str> =
            d.values().filter_map(|x| x.label.as_deref()).collect();
        Some(format!("{} clusters, {} labels", d.len(), labels.len()))
    };
    let from = source_name(m);
    let holds = match kind {
        Kind::Run => "no annotation".into(),
        Kind::Round => clusters().unwrap_or_else(|| "annotated".into()),
        Kind::Trajectory => {
            let roots = t
                .settings
                .as_ref()
                .and_then(|s| s.get("roots"))
                .and_then(|r| r.as_array())
                .map(|r| r.iter().filter_map(|v| v.as_str()).collect::<Vec<_>>())
                .unwrap_or_default();
            let mut s = "trajectory".to_string();
            if let Some(f) = &from {
                s += &format!(" on {}", derive_out_prefix(f));
            }
            if !roots.is_empty() {
                s += &format!(", from {}", roots.join(", "));
            }
            s
        }
    };
    let modified = std::fs::metadata(&path)
        .and_then(|md| md.modified())
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map_or(0, |d| d.as_secs());
    Member {
        path,
        kind,
        holds,
        modified,
    }
}

/// The member that is `shown`, compared as files (a run opened through a
/// link or another spelling of its directory is still found).
pub fn current(members: &[Member], shown: &Path) -> Option<usize> {
    members.iter().position(|m| same_file(&m.path, shown))
}

/// Whether the family has an annotation round.
pub fn has_rounds(members: &[Member]) -> bool {
    members.iter().any(|m| m.kind == Kind::Round)
}

#[cfg(test)]
#[path = "tests/runs.rs"]
pub(super) mod tests;
