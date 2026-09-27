//! Manifest glue for annotation rounds ([`crate::annotate::rounds`]): read a
//! round's tables, write the next round, and the `review` / `relabel`
//! commands.
//!
//! A round writes, beside its manifest:
//! - `{out}.clusters.parquet`: per-cell cluster id (`cluster.clusters`)
//! - `{out}.argmax.tsv`: per-cell label (`annotate.argmax`)
//! - `{out}.cluster_summary.json`: the digest (`annotate.cluster_summary`)
//! - `{out}.annotation_log.jsonl`: this round's decisions (`annotate.log`)
//! - `{out}.annotation_history.json`: every round's decisions per cluster,
//!   newest first (`annotate.history`)
//! - when a round edits markers, `{out}.markers.tsv` (`annotate.markers`, the
//!   previous round's panel with the edits applied) and
//!   `{out}.marker_history.json` (`annotate.marker_history`, per cell type)

use crate::annotate::rounds::{
    self, apply, apply_markers, digest, next_cluster_id, parse_cluster_id, prepend_history,
    ClPlacement, ClusterId, Decision, Digest, Evidence, History, MarkerHistory, Table, Term,
};
use crate::manifest::run::{
    self, annotated_path, parent_dir, rel_to_manifest, resolve, same_file, Loaded,
};
use anyhow::{Context, Result};
use clap::Args;
use enrichment::UNASSIGNED_LABEL;
use legume_numeric::matrix::common_io::{mkdir_parent, write_lines};
use legume_numeric::matrix::dense_mat_io::{read_mat, Mat};
use legume_numeric::matrix::traits::IoOps;
use log::info;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

pub const CLUSTERS: &str = ".clusters.parquet";
pub const SUMMARY: &str = ".cluster_summary.json";
pub const LOG: &str = ".annotation_log.jsonl";
pub const HISTORY: &str = ".annotation_history.json";
pub const MARKERS: &str = ".markers.tsv";
pub const MARKER_HISTORY: &str = ".marker_history.json";

//////////////////
// cell tables  //
//////////////////

/// Cell names and, per cell, its cluster id (`None` when unassigned).
pub type CellClusters = (Vec<Box<str>>, Vec<Option<ClusterId>>);

/// Per-cell cluster ids from a cluster table: its `cluster` column (else the
/// first), with NaN, negative ids, and rows whose `entropy` is not finite
/// read as unassigned.
pub fn read_clusters(path: &str) -> Result<CellClusters> {
    let m = read_mat(path).with_context(|| format!("reading clusters {path}"))?;
    anyhow::ensure!(m.mat.ncols() >= 1, "{path}: no cluster column");
    let col = |name: &str| m.cols.iter().position(|c| c.as_ref() == name);
    let label = col("cluster").unwrap_or(0);
    let entropy = col("entropy");
    let ids = (0..m.mat.nrows())
        .map(|i| {
            let v = m.mat[(i, label)];
            let empty = entropy.is_some_and(|e| !m.mat[(i, e)].is_finite());
            (v.is_finite() && v >= 0.0 && !empty).then_some(v as ClusterId)
        })
        .collect();
    Ok((m.rows, ids))
}

pub fn write_clusters(path: &str, cells: &[Box<str>], ids: &[Option<ClusterId>]) -> Result<()> {
    let mut m = Mat::zeros(cells.len(), 1);
    for (i, id) in ids.iter().enumerate() {
        m[(i, 0)] = id.map_or(f32::NAN, |x| x as f32);
    }
    let cols: Vec<Box<str>> = vec!["cluster".into()];
    m.to_parquet_with_names(path, (Some(cells), Some("cell")), Some(&cols))?;
    info!("wrote {path}");
    Ok(())
}

/// `cell⇥cell_type⇥probability`, header first.
fn read_argmax(path: &str) -> Result<HashMap<String, (String, f32)>> {
    let raw = fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
    Ok(raw
        .lines()
        .skip(1)
        .filter_map(|l| {
            let mut f = l.split('\t');
            let (cell, label) = (f.next()?, f.next()?);
            let p = f.next().and_then(|p| p.parse().ok()).unwrap_or(f32::NAN);
            Some((cell.to_string(), (label.to_string(), p)))
        })
        .collect())
}

fn write_argmax(
    path: &str,
    cells: &[Box<str>],
    labels: &[Option<String>],
    probs: &[f32],
) -> Result<()> {
    let mut lines: Vec<Box<str>> = vec!["cell\tcell_type\tprobability".into()];
    for ((c, l), p) in cells.iter().zip(labels).zip(probs) {
        let l = l.as_deref().unwrap_or(UNASSIGNED_LABEL);
        lines.push(format!("{c}\t{l}\t{p:.4}").into_boxed_str());
    }
    write_lines(&lines, path)?;
    info!("wrote {path}");
    Ok(())
}

/// One round's cells: ids from `cluster.clusters`, labels and their
/// probabilities from `annotate.argmax`, joined by cell name.
struct Cells {
    names: Vec<Box<str>>,
    clusters: Vec<Option<ClusterId>>,
    labels: Vec<Option<String>>,
    probs: Vec<f32>,
}

fn read_cells(loaded: &Loaded) -> Result<Cells> {
    let a = &loaded.manifest.annotate;
    let clusters_rel = loaded
        .manifest
        .cluster
        .clusters
        .as_deref()
        .with_context(|| {
            format!(
                "{} records no `cluster.clusters`; annotate it with this build of lupin first",
                loaded.file.display()
            )
        })?;
    let (names, clusters) = read_clusters(&resolve(&loaded.dir, clusters_rel))?;
    let argmax = match a.argmax.as_deref() {
        Some(rel) => read_argmax(&resolve(&loaded.dir, rel))?,
        None => HashMap::new(),
    };
    let (labels, probs) = names
        .iter()
        .map(|c| match argmax.get(c.as_ref()) {
            Some((l, p)) if l != UNASSIGNED_LABEL => (Some(l.clone()), *p),
            Some((_, p)) => (None, *p),
            None => (None, f32::NAN),
        })
        .unzip();
    Ok(Cells {
        names,
        clusters,
        labels,
        probs,
    })
}

//////////////
// evidence //
//////////////

fn read_table(path: &str) -> Result<Table> {
    let m = Mat::from_parquet_with_row_names(path, Some(0))
        .with_context(|| format!("reading {path}"))?;
    let (keep, rows): (Vec<usize>, Vec<ClusterId>) = m
        .rows
        .iter()
        .enumerate()
        .filter_map(|(i, r)| parse_cluster_id(r).map(|id| (i, id)))
        .unzip();
    let values = keep
        .iter()
        .flat_map(|&i| (0..m.mat.ncols()).map(move |j| (i, j)))
        .map(|(i, j)| m.mat[(i, j)])
        .collect();
    Ok(Table {
        rows,
        cols: m.cols.iter().map(ToString::to_string).collect(),
        values,
    })
}

/// A header-first TSV as rows of `column → value`.
fn read_tsv(path: &str) -> Result<Vec<HashMap<String, String>>> {
    let raw = fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
    let mut lines = raw.lines();
    let header: Vec<&str> = lines.next().unwrap_or_default().split('\t').collect();
    Ok(lines
        .map(|l| {
            header
                .iter()
                .zip(l.split('\t'))
                .map(|(h, v)| ((*h).to_string(), v.to_string()))
                .collect()
        })
        .collect())
}

/// The GO/GMT signature TSV: its first column is the cluster, rows in rank order.
fn read_terms(path: &str, source: &str) -> Result<BTreeMap<ClusterId, Vec<Term>>> {
    let raw = fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
    let group = raw.lines().next().and_then(|h| h.split('\t').next());
    let group = group.unwrap_or("cluster").to_string();
    let mut out: BTreeMap<ClusterId, Vec<Term>> = BTreeMap::new();
    for row in read_tsv(path)? {
        let Some(id) = row.get(&group).and_then(|g| parse_cluster_id(g)) else {
            continue;
        };
        let name = row.get("term_name").filter(|n| !n.is_empty());
        let term = name
            .or_else(|| row.get("term_id"))
            .cloned()
            .unwrap_or_default();
        let effect = row
            .get("effect")
            .and_then(|e| e.parse().ok())
            .unwrap_or(f32::NAN);
        out.entry(id).or_default().push(Term {
            source: source.to_string(),
            term,
            effect,
            q: None,
        });
    }
    Ok(out)
}

fn read_cl(path: &str) -> Result<BTreeMap<ClusterId, ClPlacement>> {
    Ok(read_tsv(path)?
        .into_iter()
        .filter_map(|row| {
            let id = parse_cluster_id(row.get("cluster")?)?;
            Some((
                id,
                ClPlacement {
                    id: row.get("assigned_cl")?.clone(),
                    name: row.get("assigned_name").cloned().unwrap_or_default(),
                    abstained: row.get("abstained").is_some_and(|a| a == "true"),
                },
            ))
        })
        .collect())
}

/// Whatever evidence this round records; a missing table is skipped.
fn read_evidence(loaded: &Loaded) -> Result<Evidence> {
    let a = &loaded.manifest.annotate;
    let at = |rel: &Option<String>| rel.as_deref().map(|r| resolve(&loaded.dir, r));
    let mut ev = Evidence::default();
    if let Some(p) = at(&a.cluster_celltype_q_values).or_else(|| at(&a.cluster_term_q)) {
        ev.q = Some(read_table(&p)?);
    }
    if let Some(p) = at(&a.cluster_celltype_support) {
        ev.support = Some(read_table(&p)?);
    }
    if let Some(p) = at(&a.ontology_signature) {
        let gmt = a
            .settings
            .as_ref()
            .and_then(|s| s.pointer("/enrichment/gmt"))
            .is_some_and(|g| !g.is_null());
        ev.terms = read_terms(&p, if gmt { "gmt" } else { "go" })?;
    }
    if let Some(p) = at(&a.ontology_assignment) {
        ev.cl = read_cl(&p)?;
    }
    Ok(ev)
}

/// Build this round's digest, write `{out}.cluster_summary.json` and record
/// it in the manifest (which the caller saves).
pub fn write_summary(loaded: &mut Loaded, out_prefix: &str) -> Result<BTreeMap<ClusterId, Digest>> {
    let cells = read_cells(loaded)?;
    let d = digest(&cells.clusters, &cells.labels, &read_evidence(loaded)?);
    let path = format!("{out_prefix}{SUMMARY}");
    fs::write(&path, serde_json::to_string_pretty(&d)?)?;
    info!("wrote {path}");
    loaded.manifest.annotate.cluster_summary = Some(rel_to_manifest(&loaded.dir, &path));
    Ok(d)
}

/// The round's marker panel as `(feature, cell type)` pairs, and its marker
/// history; both empty when the round records none.
fn read_markers(loaded: &Loaded) -> Result<(Vec<(String, String)>, MarkerHistory)> {
    let a = &loaded.manifest.annotate;
    let pairs = match a.markers.as_deref().filter(|m| !m.is_empty()) {
        Some(rel) => data_beans::aux::gene_sets::read_membership_pairs(&resolve(&loaded.dir, rel))?
            .into_iter()
            .map(|(g, t)| (g.into_string(), t.into_string()))
            .collect(),
        None => Vec::new(),
    };
    let history = match a.marker_history.as_deref() {
        Some(rel) => {
            let path = resolve(&loaded.dir, rel);
            let raw = fs::read_to_string(&path).with_context(|| format!("reading {path}"))?;
            serde_json::from_str(&raw).with_context(|| format!("parsing {path}"))?
        }
        None => MarkerHistory::new(),
    };
    Ok((pairs, history))
}

fn read_history(loaded: &Loaded) -> Result<History> {
    let Some(rel) = loaded.manifest.annotate.history.as_deref() else {
        return Ok(History::new());
    };
    let path = resolve(&loaded.dir, rel);
    let raw = fs::read_to_string(&path).with_context(|| format!("reading {path}"))?;
    serde_json::from_str(&raw).with_context(|| format!("parsing {path}"))
}

/////////////
// relabel //
/////////////

#[derive(Args, Debug, Default)]
pub struct RelabelArgs {
    #[arg(
        long,
        short = 'f',
        help = "The round to start from (an annotated manifest or its prefix)"
    )]
    pub from: Box<str>,

    #[arg(
        long,
        short = 'd',
        help = "Decisions, one JSON object per line, or `-` for stdin (see `lupin review --help`)"
    )]
    pub decisions: Box<str>,

    #[arg(long, short = 'o', help = "Output prefix for the new round")]
    pub out: Option<Box<str>>,

    #[arg(
        long,
        conflicts_with = "out",
        help = "Write the next round of -f's chain ({chain}.r<N+1>) and print its path; \
                refused unless -f is the chain's latest round"
    )]
    pub next: bool,

    #[arg(
        long,
        conflicts_with_all = ["next", "out"],
        help = "Keep running: each batch of lines appended to the decisions file becomes the next \
                round of -f's chain; every decision must name its `round`"
    )]
    pub watch: bool,
}

/// Apply decisions to a round and write the next one; with `--watch`, keep
/// applying whatever is appended to the decisions file.
///
/// `--next` and `--watch` grow `-f`'s chain through [`Chain::write_next`]:
/// one lock, one rule that decisions are made on the latest round. `-o`
/// writes a round wherever it is told, and like every path never
/// overwrites one.
pub fn run_relabel(args: &RelabelArgs) -> Result<()> {
    if args.watch {
        anyhow::ensure!(
            &*args.decisions != "-",
            "--watch needs a decisions file, not stdin"
        );
        return watch(args);
    }
    let (raw, decisions_dir) = if &*args.decisions == "-" {
        (
            std::io::read_to_string(std::io::stdin())?,
            PathBuf::from("."),
        )
    } else {
        let raw = fs::read_to_string(&*args.decisions)
            .with_context(|| format!("reading {}", args.decisions))?;
        (raw, parent_dir(Path::new(&*args.decisions)))
    };
    let decisions = parse_decisions(numbered(&raw), &args.decisions)?;
    let source = run::load(&args.from)?;
    let written = match &args.out {
        Some(out) => relabel(&source, decisions, &decisions_dir, out)?,
        None if args.next => {
            let chain = Chain::of(&source.file);
            let lock = chain.lock()?;
            chain.write_next(&lock, &source.file, decisions, &decisions_dir)?
        }
        None => anyhow::bail!("give -o <prefix>, or --next to continue -f's chain"),
    };
    // The caller (a viewer, say) opens what was written.
    println!("{}", written.display());
    Ok(())
}

/// A text's lines, numbered from 1.
fn numbered(raw: &str) -> impl Iterator<Item = (usize, &str)> {
    raw.lines().enumerate().map(|(n, l)| (n + 1, l))
}

/// Numbered lines to decisions; blank lines are skipped.
fn parse_decisions<'a>(
    lines: impl IntoIterator<Item = (usize, &'a str)>,
    file: &str,
) -> Result<Vec<Decision>> {
    let decisions: Vec<Decision> = lines
        .into_iter()
        .filter(|(_, l)| !l.trim().is_empty())
        .map(|(n, l)| {
            serde_json::from_str(l).with_context(|| format!("{file} line {n}: not a decision"))
        })
        .collect::<Result<_>>()?;
    anyhow::ensure!(!decisions.is_empty(), "{file}: no decisions");
    Ok(decisions)
}

fn now() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs());
    rounds::utc_timestamp(secs)
}

///////////
// chain //
///////////

/// How long a one-shot writer waits for another to finish its round.
const LOCK_WAIT: Duration = Duration::from_secs(10);

/// Rounds `{prefix}.r1`, `{prefix}.r2`, ... grown from `base`. The latest
/// round is the highest one on disk, so every writer agrees on it.
pub struct Chain {
    prefix: String,
    base: PathBuf,
    /// The base's own round number when it is a round of this chain, else 0.
    base_k: usize,
}

impl Chain {
    /// The chain `round` belongs to: `X.rK.senna.json` is round K of `X`;
    /// any other manifest `X.*.json` starts chain `X`.
    #[must_use]
    pub fn of(round: &Path) -> Self {
        let stem = run::derive_out_prefix(&round.to_string_lossy());
        let (prefix, base_k) = split_round(&stem).unwrap_or((stem.as_str(), 0));
        Self {
            prefix: prefix.to_string(),
            base: round.to_path_buf(),
            base_k,
        }
    }

    fn round_prefix(&self, k: usize) -> String {
        format!("{}.r{k}", self.prefix)
    }

    /// Every round manifest on disk after the base, in order. Found by
    /// listing the directory, so a missing round does not hide later ones.
    fn rounds(&self) -> Vec<(usize, PathBuf)> {
        // `annotated_path(base, "")` is just the manifest suffix.
        let suffix = annotated_path(&self.base, "")
            .to_string_lossy()
            .into_owned();
        let prefix = Path::new(&self.prefix);
        let name = prefix.file_name().map(|n| n.to_string_lossy().into_owned());
        let Ok(entries) = fs::read_dir(parent_dir(prefix)) else {
            return Vec::new();
        };
        let mut found: Vec<(usize, PathBuf)> = entries
            .flatten()
            .filter_map(|e| {
                let file = e.file_name().to_string_lossy().into_owned();
                let (p, k) = split_round(file.strip_suffix(suffix.as_str())?)?;
                (Some(p) == name.as_deref() && k > self.base_k).then(|| (k, e.path()))
            })
            .collect();
        found.sort_unstable_by_key(|(k, _)| *k);
        found
    }

    /// The highest round on disk and its manifest; the base before any.
    #[must_use]
    pub fn latest(&self) -> (usize, PathBuf) {
        self.rounds()
            .pop()
            .unwrap_or_else(|| (self.base_k, self.base.clone()))
    }

    /// Write the next round from `made_on`, the round the decisions were made
    /// on, which must still be the latest: if another writer has moved on,
    /// the decisions' cluster ids may no longer mean what the decider saw.
    /// Holding `_lock` is what makes "latest" stay true while writing.
    fn write_next(
        &self,
        _lock: &ChainLock,
        made_on: &Path,
        decisions: Vec<Decision>,
        decisions_dir: &Path,
    ) -> Result<PathBuf> {
        let (k, latest) = self.latest();
        anyhow::ensure!(
            same_file(made_on, &latest),
            "{} is not the latest round (that is {}); reload and decide again",
            made_on.display(),
            latest.display()
        );
        let source = run::load(&latest.to_string_lossy())?;
        relabel(&source, decisions, decisions_dir, &self.round_prefix(k + 1))
    }

    /// Hold the chain's lock, waiting up to [`LOCK_WAIT`] for another writer.
    pub fn lock(&self) -> Result<ChainLock> {
        self.lock_within(LOCK_WAIT)
    }

    /// Hold the chain's lock only if it is free now.
    pub fn try_lock(&self) -> Result<ChainLock> {
        self.lock_within(Duration::ZERO)
    }

    /// `{prefix}.relabel.lock`, as an OS file lock: a writer that dies
    /// releases it, and the file itself stays.
    fn lock_within(&self, wait: Duration) -> Result<ChainLock> {
        let path = format!("{}.relabel.lock", self.prefix);
        mkdir_parent(&path)?;
        let file = fs::OpenOptions::new()
            .create(true)
            .truncate(false)
            .write(true)
            .open(&path)
            .with_context(|| format!("opening {path}"))?;
        let start = std::time::Instant::now();
        loop {
            match file.try_lock() {
                Ok(()) => return Ok(ChainLock { _file: file }),
                Err(fs::TryLockError::WouldBlock) => {
                    anyhow::ensure!(
                        start.elapsed() < wait,
                        "another relabel is writing this chain ({path} held)"
                    );
                    std::thread::sleep(Duration::from_millis(50));
                }
                Err(fs::TryLockError::Error(e)) => {
                    return Err(e).with_context(|| format!("locking {path}"))
                }
            }
        }
    }
}

/// `X.rK` as `(X, K)`.
fn split_round(stem: &str) -> Option<(&str, usize)> {
    let (prefix, k) = stem.rsplit_once(".r")?;
    let k = k
        .parse()
        .ok()
        .filter(|_| k.bytes().all(|b| b.is_ascii_digit()))?;
    Some((prefix, k))
}

/// The chain's lock, released when dropped (the OS releases it too if the
/// process dies first).
pub struct ChainLock {
    _file: fs::File,
}

/// Apply `decisions` to the round `source` and write the next round at
/// `out`; returns its manifest. Nothing is written if a decision is refused,
/// and an existing round is never overwritten.
fn relabel(
    source: &Loaded,
    mut decisions: Vec<Decision>,
    decisions_dir: &Path,
    out: &str,
) -> Result<PathBuf> {
    let manifest_path = annotated_path(&source.file, out);
    anyhow::ensure!(
        !manifest_path.exists(),
        "{} exists; rounds are never overwritten",
        manifest_path.display()
    );
    check_round(source, &decisions, decisions_dir)?;
    mkdir_parent(out)?;
    let round = manifest_path
        .file_name()
        .map_or_else(String::new, |n| n.to_string_lossy().into_owned());

    let mut cells = read_cells(source)?;
    let older = read_history(source)?;
    let next_id = next_cluster_id(&cells.clusters, &older);
    let newer = apply(
        &mut decisions,
        &mut cells.clusters,
        &mut cells.labels,
        next_id,
        &round,
        &now(),
    )?;
    let (mut markers, older_markers) = read_markers(source)?;
    let (markers_edited, newer_markers) = apply_markers(&decisions, &mut markers, &round)?;
    info!("applied {} decision(s)", decisions.len());

    let clusters_path = format!("{out}{CLUSTERS}");
    write_clusters(&clusters_path, &cells.names, &cells.clusters)?;
    let argmax_path = format!("{out}.argmax.tsv");
    write_argmax(&argmax_path, &cells.names, &cells.labels, &cells.probs)?;
    let log_path = format!("{out}{LOG}");
    let log: Vec<Box<str>> = decisions
        .iter()
        .map(|d| serde_json::to_string(d).map(String::into_boxed_str))
        .collect::<Result<_, _>>()?;
    write_lines(&log, &log_path)?;
    let history_path = format!("{out}{HISTORY}");
    let history = prepend_history(older, newer);
    fs::write(&history_path, serde_json::to_string_pretty(&history)?)?;
    info!("wrote {history_path}");

    let marker_paths = if markers_edited {
        let panel = format!("{out}{MARKERS}");
        let mut lines: Vec<Box<str>> = vec!["gene\tcelltype".into()];
        lines.extend(
            markers
                .iter()
                .map(|(g, t)| format!("{g}\t{t}").into_boxed_str()),
        );
        write_lines(&lines, &panel)?;
        info!("wrote {panel}");
        let hist = format!("{out}{MARKER_HISTORY}");
        let merged = prepend_history(older_markers, newer_markers);
        fs::write(&hist, serde_json::to_string_pretty(&merged)?)?;
        info!("wrote {hist}");
        Some((panel, hist))
    } else {
        None
    };

    let mut next = source.copy_to(manifest_path)?;
    let rel = |p: &str| Some(rel_to_manifest(&next.dir, p));
    if let Some((panel, hist)) = &marker_paths {
        next.manifest.annotate.markers = rel(panel);
        next.manifest.annotate.marker_history = rel(hist);
    }
    next.manifest.cluster.clusters = rel(&clusters_path);
    next.manifest.annotate.argmax = rel(&argmax_path);
    next.manifest.annotate.log = rel(&log_path);
    next.manifest.annotate.history = rel(&history_path);
    write_summary(&mut next, out)?;
    next.manifest.save(&next.file)?;
    Ok(next.file)
}

/// Refuse decisions that name a round other than `source`, the one they are
/// being applied to.
fn check_round(source: &Loaded, decisions: &[Decision], decisions_dir: &Path) -> Result<()> {
    let current = source
        .file
        .canonicalize()
        .unwrap_or_else(|_| source.file.clone());
    for d in decisions {
        let Some(named) = &d.round else { continue };
        anyhow::ensure!(
            same_file(Path::new(&resolve(decisions_dir, named)), &current),
            "ids refer to {named}, but the round being relabelled is {}; reload and decide again",
            source.file.display()
        );
    }
    Ok(())
}

///////////
// watch //
///////////

pub const STATUS: &str = ".relabel_status.json";
const POLL: Duration = Duration::from_millis(500);

/// What `relabel --watch` has done, rewritten after every batch at
/// `{chain}.relabel_status.json`. Rounds are named relative to that file,
/// which sits beside them.
#[derive(Serialize, Deserialize, Debug, Default)]
pub struct WatchStatus {
    /// The decisions file being watched.
    pub decisions: String,
    /// Lines of it already handled, applied or refused.
    pub processed_lines: usize,
    /// The round the watcher started from (`-f`).
    #[serde(default)]
    pub base: String,
    /// Every round in order, starting with `base`; filled in when written.
    pub rounds: Vec<String>,
    /// The chain's latest round; filled in when written.
    pub latest: String,
    /// Why the last batch was refused; `None` once one succeeds.
    pub error: Option<WatchError>,
    pub updated: String,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct WatchError {
    /// First and last line of the refused batch.
    pub lines: [usize; 2],
    pub message: String,
}

fn watch(args: &RelabelArgs) -> Result<()> {
    let (mut status, status_path, chain) = begin_watch(args)?;
    info!(
        "watching {} (from line {}); rounds go to {}.r<N>; status in {}",
        args.decisions,
        status.processed_lines + 1,
        chain.prefix,
        status_path.display()
    );
    loop {
        watch_step(&mut status, &status_path, &chain)?;
        std::thread::sleep(POLL);
    }
}

/// Start or resume a watch and write its status before any decision, so a
/// viewer can find the watcher at once. The chain is `-f`'s own.
fn begin_watch(args: &RelabelArgs) -> Result<(WatchStatus, PathBuf, Chain)> {
    let chain = Chain::of(&run::load(&args.from)?.file);
    mkdir_parent(&chain.prefix)?;
    let status_path = PathBuf::from(format!("{}{STATUS}", chain.prefix));
    let mut status = start_watch(&status_path, &chain.base, &args.decisions)?;
    write_status(&mut status, &status_path, &chain)?;
    Ok((status, status_path, chain))
}

/// Resume from an existing status for the same decisions file, else start
/// at `base` with nothing processed.
fn start_watch(status_path: &Path, base: &Path, decisions: &str) -> Result<WatchStatus> {
    let decisions_rel = rel_to_manifest(&parent_dir(status_path), decisions);
    if let Ok(raw) = fs::read_to_string(status_path) {
        let prev: WatchStatus = serde_json::from_str(&raw)
            .with_context(|| format!("parsing {}", status_path.display()))?;
        if prev.decisions == decisions_rel {
            return Ok(prev);
        }
        log::warn!(
            "{} watched {}; starting over for {decisions}",
            status_path.display(),
            prev.decisions
        );
    }
    Ok(WatchStatus {
        decisions: decisions_rel,
        base: file_name(base),
        ..WatchStatus::default()
    })
}

fn file_name(p: &Path) -> String {
    p.file_name()
        .map_or_else(String::new, |n| n.to_string_lossy().into_owned())
}

/// Apply the complete lines appended since the last step as one round.
/// Returns whether there was anything new. A refused batch is recorded in
/// the status and skipped, so the next append is not stuck behind it.
pub fn watch_step(status: &mut WatchStatus, status_path: &Path, chain: &Chain) -> Result<bool> {
    // Follow rounds other writers added since the last step.
    if file_name(&chain.latest().1) != status.latest {
        write_status(status, status_path, chain)?;
    }
    let decisions = resolve(&parent_dir(status_path), &status.decisions);
    let raw = match fs::read_to_string(&decisions) {
        Ok(r) => r,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(false),
        Err(e) => return Err(e).with_context(|| format!("reading {decisions}")),
    };
    // A line is complete once its newline is written.
    let complete = raw.matches('\n').count();
    if complete <= status.processed_lines {
        return Ok(false);
    }
    let batch: Vec<(usize, &str)> = numbered(&raw)
        .skip(status.processed_lines)
        .take(complete - status.processed_lines)
        .collect();
    let first_line = status.processed_lines + 1;
    if batch.iter().all(|(_, l)| l.trim().is_empty()) {
        status.processed_lines = complete;
        return Ok(false);
    }
    // Contention is not a verdict on the decisions: leave them for the next
    // poll rather than blocking it.
    let lock = match chain.try_lock() {
        Ok(lock) => lock,
        Err(e) => {
            log::warn!("lines {first_line}-{complete} wait: {e:#}");
            return Ok(false);
        }
    };
    status.processed_lines = complete;

    match apply_batch(chain, &lock, batch, &decisions) {
        Ok(file) => {
            info!("round written: {}", file.display());
            status.error = None;
        }
        Err(e) => {
            log::warn!("lines {first_line}-{complete} refused: {e:#}");
            status.error = Some(WatchError {
                lines: [first_line, complete],
                message: format!("{e:#}"),
            });
        }
    }
    write_status(status, status_path, chain)?;
    Ok(true)
}

/// One watched batch as the next round. With no `-f` per batch, the round
/// the decisions were made on comes from their `round`, which is required.
fn apply_batch(
    chain: &Chain,
    lock: &ChainLock,
    batch: Vec<(usize, &str)>,
    decisions: &str,
) -> Result<PathBuf> {
    let ds = parse_decisions(batch, decisions)?;
    let dir = parent_dir(Path::new(decisions));
    let made_on = ds
        .iter()
        .map(|d| d.round.as_deref())
        .collect::<Option<Vec<_>>>()
        .and_then(|r| r.first().map(|r| resolve(&dir, r)))
        .context("every watched decision must name the round it was made on (`round`)")?;
    chain.write_next(lock, Path::new(&made_on), ds, &dir)
}

/// Stamp the status with the chain as it is on disk and write it whole,
/// then rename it into place, so a reader never sees half a file.
fn write_status(status: &mut WatchStatus, status_path: &Path, chain: &Chain) -> Result<()> {
    status.rounds = std::iter::once(status.base.clone())
        .chain(chain.rounds().iter().map(|(_, r)| file_name(r)))
        .collect();
    status.latest = status.rounds.last().cloned().unwrap_or_default();
    status.updated = now();
    let tmp = status_path.with_extension("json.tmp");
    fs::write(&tmp, serde_json::to_string_pretty(status)?)?;
    fs::rename(&tmp, status_path)?;
    Ok(())
}

////////////
// review //
////////////

#[derive(Args, Debug)]
#[command(after_long_help = DECISIONS_HELP)]
pub struct ReviewArgs {
    #[arg(
        long,
        short = 'f',
        help = "An annotated manifest (a round) or its prefix"
    )]
    pub from: Box<str>,

    #[arg(long, short = 'c', help = "Only these clusters (repeatable)")]
    pub cluster: Vec<ClusterId>,

    #[arg(
        long,
        help = "Print JSON (digest and history per cluster) instead of text"
    )]
    pub json: bool,
}

const DECISIONS_HELP: &str = "\
Decisions file for `lupin relabel -d`: one JSON object per line.

  {\"cluster\": 3, \"action\": \"label\", \"label\": \"CT1\",
   \"evidence\": [{\"kind\": \"marker\", \"term\": \"CT1\", \"q\": 0.001}],
   \"alternatives\": [{\"label\": \"CT2\", \"why_not\": \"weaker support\"}],
   \"rationale\": \"...\", \"decided_by\": \"user\"}

action:     label (one cluster), merge (\"clusters\": [..], fresh id), keep,
            markers_add / markers_drop (\"label\": cell type, \"features\": [..];
            edits the round's marker panel, which the next annotate uses)
decided_by: user | agent_proposed_user_accepted | user_override
round:      optional; the round whose ids the decision names, relative to the
            decisions file. A decision on any other round is refused.
Every decision needs a rationale; it is kept in the round's history.
With several writers (a viewer and an agent), append each decision as one
write of one newline-terminated line under 4 KB, so lines never interleave.";

/// Print a round's digest and history, for a person or an agent deciding.
pub fn run_review(args: &ReviewArgs) -> Result<()> {
    let loaded = run::load(&args.from)?;
    let digests: BTreeMap<ClusterId, Digest> = match &loaded.manifest.annotate.cluster_summary {
        Some(rel) if Path::new(&resolve(&loaded.dir, rel)).is_file() => {
            serde_json::from_str(&fs::read_to_string(resolve(&loaded.dir, rel))?)?
        }
        _ => {
            let cells = read_cells(&loaded)?;
            digest(&cells.clusters, &cells.labels, &read_evidence(&loaded)?)
        }
    };
    let history = read_history(&loaded)?;
    let wanted = |id: &ClusterId| args.cluster.is_empty() || args.cluster.contains(id);

    if args.json {
        let out: BTreeMap<ClusterId, serde_json::Value> = digests
            .iter()
            .filter(|(id, _)| wanted(id))
            .map(|(id, d)| {
                let h = history.get(id).cloned().unwrap_or_default();
                (*id, serde_json::json!({ "digest": d, "history": h }))
            })
            .collect();
        println!("{}", serde_json::to_string_pretty(&out)?);
        return Ok(());
    }
    for (id, d) in digests.iter().filter(|(id, _)| wanted(id)) {
        print!(
            "{}",
            render(*id, d, history.get(id).map_or(&[][..], Vec::as_slice))
        );
    }
    Ok(())
}

fn render(id: ClusterId, d: &Digest, history: &[rounds::HistoryEntry]) -> String {
    use std::fmt::Write;
    let num = |v: Option<f32>| v.map_or_else(|| "-".into(), |v| format!("{v:.3}"));
    let mut s = String::new();
    let label = d.label.as_deref().unwrap_or("unassigned");
    let _ = writeln!(s, "C{id}  n={}  label={label}", d.size);
    for c in &d.calls {
        let _ = writeln!(
            s,
            "  call  {:<24} q={} support={}",
            c.label,
            num(c.q),
            num(c.support)
        );
    }
    for t in &d.terms {
        let _ = writeln!(
            s,
            "  {:<4}  {:<40} effect={:.3}",
            t.source, t.term, t.effect
        );
    }
    if let Some(cl) = &d.cl {
        let abst = if cl.abstained { " (abstained)" } else { "" };
        let _ = writeln!(s, "  CL    {} {}{abst}", cl.id, cl.name);
    }
    for h in history {
        let label = h.label.as_deref().unwrap_or("-");
        let _ = writeln!(
            s,
            "  hist  {} {} -> {label} by {}: {}",
            h.round,
            h.action.as_str(),
            h.decided_by.as_str(),
            h.rationale
        );
    }
    s
}

#[cfg(test)]
#[path = "rounds_tests.rs"]
mod tests;
