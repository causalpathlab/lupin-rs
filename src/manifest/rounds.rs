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

use crate::annotate::outputs::{
    write_cluster_tables, CLUSTER_CELLTYPE_ES_STD, CLUSTER_CELLTYPE_NES, CLUSTER_CELLTYPE_P,
    CLUSTER_CELLTYPE_Q, CLUSTER_CELLTYPE_Q_VALUES, CLUSTER_CELLTYPE_Z,
};
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
use legume_numeric::matrix::traits::MatWithNames;
use log::info;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet, HashMap};
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
pub(crate) fn read_argmax(path: &str) -> Result<HashMap<String, (String, f32)>> {
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

pub(super) fn write_argmax(
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
pub(crate) struct Cells {
    pub(crate) names: Vec<Box<str>>,
    pub(crate) clusters: Vec<Option<ClusterId>>,
    pub(crate) labels: Vec<Option<String>>,
    pub(crate) probs: Vec<f32>,
}

pub(crate) fn read_cells(loaded: &Loaded) -> Result<Cells> {
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

pub(crate) fn read_table(path: &str) -> Result<Table> {
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
        let num = |c: &str| row.get(c).and_then(|v| v.parse::<f32>().ok());
        out.entry(id).or_default().push(Term {
            source: source.to_string(),
            term,
            effect: num("effect").unwrap_or(f32::NAN),
            q: num("q"),
            p: num("p"),
            nes: num("nes"),
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

/// Each cluster's top GO (or GMT) terms, when the round scored any.
pub(crate) fn read_gene_set_terms(loaded: &Loaded) -> Result<BTreeMap<ClusterId, Vec<Term>>> {
    let a = &loaded.manifest.annotate;
    let Some(rel) = a.ontology_signature.as_deref() else {
        return Ok(BTreeMap::new());
    };
    let gmt = a
        .settings
        .as_ref()
        .and_then(|s| s.pointer("/enrichment/gmt"))
        .is_some_and(|g| !g.is_null());
    read_terms(&resolve(&loaded.dir, rel), if gmt { "gmt" } else { "go" })
}

/// Whatever evidence this round records; a missing table is skipped.
fn read_evidence(loaded: &Loaded) -> Result<Evidence> {
    let a = &loaded.manifest.annotate;
    let at = |rel: &Option<String>| rel.as_deref().map(|r| resolve(&loaded.dir, r));
    let mut ev = Evidence::default();
    if let Some(p) = at(&a.cluster_celltype_q_values).or_else(|| at(&a.cluster_term_q)) {
        ev.q = Some(read_table(&p)?);
    }
    ev.terms = read_gene_set_terms(loaded)?;
    if let Some(p) = at(&a.ontology_assignment) {
        ev.cl = read_cl(&p)?;
    }
    Ok(ev)
}

/// Build this round's digest, write `{out}.cluster_summary.json` and record
/// it in the manifest (which the caller saves).
pub fn write_summary(loaded: &mut Loaded, out_prefix: &str) -> Result<BTreeMap<ClusterId, Digest>> {
    let cells = read_cells(loaded)?;
    let mut d = digest(&cells.clusters, &cells.labels, &read_evidence(loaded)?);
    // A cluster labelled with a coarse group agrees with a top call that is
    // one of the group's types.
    if let Some(rel) = loaded.manifest.annotate.celltype_tree.as_deref() {
        let tree: crate::annotate::celltype_tree::TypeTree =
            serde_json::from_str(&fs::read_to_string(resolve(&loaded.dir, rel))?)?;
        let index = tree.index();
        for digest in d.values_mut() {
            let label = digest.label.clone();
            if let (Some(e), Some(label)) = (digest.evidence.as_mut(), label) {
                let group = index
                    .get(&crate::annotate::markers::label_key(&e.top))
                    .copied();
                e.agrees |= group.is_some_and(|g| g == label);
            }
        }
    }
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
        // Found again if the run moved; else where it was, for the error.
        Some(rel) => crate::annotate::markers::read_marker_pairs(
            &loaded
                .recorded(rel)
                .unwrap_or_else(|| resolve(&loaded.dir, rel)),
        )?
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
        conflicts_with_all = ["next", "out", "watch"],
        help = "Validate as --next would, write nothing, and print what the decisions would \
                change as JSON (marker edits re-rank cell types approximately)"
    )]
    pub preview: bool,

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
    anyhow::ensure!(!decisions.is_empty(), "{}: no decisions", args.decisions);
    let source = run::load(&args.from)?;
    if args.preview {
        let (_, latest) = Chain::of(&source.file).latest();
        ensure_latest(&source.file, &latest)?;
        println!(
            "{}",
            serde_json::to_string(&preview(&source, decisions, &decisions_dir)?)?
        );
        return Ok(());
    }
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

/// The latest round of the chain `round` starts, and the rounds after it
/// on disk (none for a fresh pass).
pub(crate) fn chain_rounds(round: &Path) -> (PathBuf, Vec<PathBuf>) {
    let chain = Chain::of(round);
    let later: Vec<PathBuf> = chain.rounds().into_iter().map(|(_, p)| p).collect();
    (chain.latest().1, later)
}

/// Set aside the rounds after `round` in its chain, under the chain's lock:
/// each round's files (`X.rK.*`) move into `X.superseded/`, so a new pass
/// that rewrote `round` starts a fresh chain. Returns how many rounds moved.
pub(crate) fn supersede_later(round: &Path) -> Result<usize> {
    let chain = Chain::of(round);
    let _lock = chain.lock()?;
    let later = chain.rounds();
    let prefix = Path::new(&chain.prefix);
    let dir = parent_dir(prefix);
    let chain_name = prefix
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let aside = dir.join(format!("{chain_name}.superseded"));
    for (k, _) in &later {
        let round_name = format!("{chain_name}.r{k}.");
        fs::create_dir_all(&aside)?;
        for e in fs::read_dir(&dir)?.flatten() {
            if e.file_name().to_string_lossy().starts_with(&round_name) {
                fs::rename(e.path(), aside.join(e.file_name()))?;
            }
        }
    }
    Ok(later.len())
}

/// A text's lines, numbered from 1.
fn numbered(raw: &str) -> impl Iterator<Item = (usize, &str)> {
    raw.lines().enumerate().map(|(n, l)| (n + 1, l))
}

/// Numbered lines to decisions. Only lines starting with `{` are decisions,
/// so an answer pasted from a chat, with its prose and code fences, reads as
/// is; each of those lines must parse.
fn parse_decisions<'a>(
    lines: impl IntoIterator<Item = (usize, &'a str)>,
    file: &str,
) -> Result<Vec<Decision>> {
    let decisions: Vec<Decision> = lines
        .into_iter()
        .filter(|(_, l)| l.trim_start().starts_with('{'))
        .map(|(n, l)| {
            serde_json::from_str(l).with_context(|| format!("{file} line {n}: not a decision"))
        })
        .collect::<Result<_>>()?;
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
        super::family::Tag::Chain.name(&self.prefix, k)
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
        ensure_latest(made_on, &latest)?;
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

/// Refuse decisions made on `made_on` when the chain has moved past it.
fn ensure_latest(made_on: &Path, latest: &Path) -> Result<()> {
    anyhow::ensure!(
        same_file(made_on, latest),
        "{} is not the latest round (that is {}); reload and decide again",
        made_on.display(),
        latest.display()
    );
    Ok(())
}

/// `X.rK` as `(X, K)`.
fn split_round(stem: &str) -> Option<(&str, usize)> {
    super::family::Tag::Chain.parse(stem)
}

/// The chain's lock, released when dropped (the OS releases it too if the
/// process dies first).
pub struct ChainLock {
    _file: fs::File,
}

/// A round with decisions applied, in memory.
struct Planned {
    /// Stamped with their timestamps, as the log keeps them.
    decisions: Vec<Decision>,
    cells: Cells,
    history: History,
    markers: Vec<(String, String)>,
    /// Every round's marker edits, when this one made any.
    marker_history: Option<MarkerHistory>,
}

/// Apply `decisions` to the round `source` in memory, as round `round`.
/// Every refusal happens here, before anything is written.
fn plan(
    source: &Loaded,
    mut decisions: Vec<Decision>,
    decisions_dir: &Path,
    round: &str,
) -> Result<Planned> {
    check_round(source, &decisions, decisions_dir)?;
    let mut cells = read_cells(source)?;
    let older = read_history(source)?;
    let next_id = next_cluster_id(&cells.clusters, &older);
    let newer = apply(
        &mut decisions,
        &mut cells.clusters,
        &mut cells.labels,
        next_id,
        round,
        &now(),
    )?;
    let (mut markers, older_markers) = read_markers(source)?;
    let (edited, newer_markers) = apply_markers(&decisions, &mut markers, round)?;
    Ok(Planned {
        decisions,
        cells,
        history: prepend_history(older, newer),
        markers,
        marker_history: edited.then(|| prepend_history(older_markers, newer_markers)),
    })
}

/// Apply `decisions` to the round `source` and write the next round at
/// `out`; returns its manifest. Nothing is written if a decision is refused,
/// and an existing round is never overwritten.
fn relabel(
    source: &Loaded,
    decisions: Vec<Decision>,
    decisions_dir: &Path,
    out: &str,
) -> Result<PathBuf> {
    let manifest_path = annotated_path(&source.file, out);
    anyhow::ensure!(
        !manifest_path.exists(),
        "{} exists; rounds are never overwritten",
        manifest_path.display()
    );
    let round = file_name(&manifest_path);
    let p = plan(source, decisions, decisions_dir, &round)?;
    info!("applied {} decision(s)", p.decisions.len());
    mkdir_parent(out)?;

    let clusters_path = format!("{out}{CLUSTERS}");
    write_clusters(&clusters_path, &p.cells.names, &p.cells.clusters)?;
    let argmax_path = format!("{out}.argmax.tsv");
    write_argmax(
        &argmax_path,
        &p.cells.names,
        &p.cells.labels,
        &p.cells.probs,
    )?;
    let log_path = format!("{out}{LOG}");
    let log: Vec<Box<str>> = p
        .decisions
        .iter()
        .map(|d| serde_json::to_string(d).map(String::into_boxed_str))
        .collect::<Result<_, _>>()?;
    write_lines(&log, &log_path)?;
    let history_path = format!("{out}{HISTORY}");
    fs::write(&history_path, serde_json::to_string_pretty(&p.history)?)?;
    info!("wrote {history_path}");

    let marker_paths = match &p.marker_history {
        Some(merged) => {
            let panel = format!("{out}{MARKERS}");
            let mut lines: Vec<Box<str>> = vec!["gene\tcelltype".into()];
            lines.extend(
                p.markers
                    .iter()
                    .map(|(g, t)| format!("{g}\t{t}").into_boxed_str()),
            );
            write_lines(&lines, &panel)?;
            info!("wrote {panel}");
            let hist = format!("{out}{MARKER_HISTORY}");
            fs::write(&hist, serde_json::to_string_pretty(merged)?)?;
            info!("wrote {hist}");
            Some((panel, hist))
        }
        None => None,
    };

    // Rescore the new clusters and panel when the pass cached its statistics.
    let rescored = super::recalibrate::rescore(source, &p.cells, &p.markers)?;
    let tables = rescored
        .as_ref()
        .map(|r| write_rescored(r, out))
        .transpose()?;

    let mut next = source.copy_to(manifest_path)?;
    // SuSiE's tables describe the source round's clusters on its panel: kept
    // while both stand, dropped once a merge regroups the clusters or the
    // markers change (the call is then the decisions over the enrichment's
    // rescored statistics).
    let regrouped = p
        .decisions
        .iter()
        .any(|d| matches!(d.action, crate::annotate::rounds::Action::Merge));
    if regrouped || p.marker_history.is_some() {
        next.manifest.annotate.drop_susie_tables();
    }
    let rel = |p: &str| Some(rel_to_manifest(&next.dir, p));
    if let Some(t) = &tables {
        let a = &mut next.manifest.annotate;
        a.cluster_celltype_q = rel(&t.q);
        a.cluster_celltype_q_values = rel(&t.q_values);
        a.cluster_celltype_p = rel(&t.p);
        a.cluster_celltype_nes = rel(&t.nes);
        let rounds = source
            .manifest
            .annotate
            .stats
            .as_ref()
            .and_then(|s| s.get("rounds_of_curation"))
            .and_then(Value::as_u64)
            .unwrap_or(0)
            + 1;
        a.stats = Some(json!({
            "kind": "post_selection",
            "rounds_of_curation": rounds,
        }));
    }
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

/// The tables a rescored round writes under its prefix.
struct RescoredTables {
    q: String,
    q_values: String,
    p: String,
    nes: String,
}

/// Write a rescored round's cluster × type tables under `out`, as a pass
/// writes them: the Q probabilities, q-values, p-values, z and NES.
fn write_rescored(r: &super::recalibrate::Rescored, out: &str) -> Result<RescoredTables> {
    let written = write_cluster_tables(
        out,
        &r.row_names(),
        &r.types,
        &[
            (&r.q_probs, CLUSTER_CELLTYPE_Q),
            (&r.q_values, CLUSTER_CELLTYPE_Q_VALUES),
            (&r.p_values, CLUSTER_CELLTYPE_P),
            (&r.nes, CLUSTER_CELLTYPE_NES),
            (&r.z, CLUSTER_CELLTYPE_ES_STD),
            (&r.probit_z, CLUSTER_CELLTYPE_Z),
        ],
    )?;
    let Ok([q, q_values, p, nes, _, _]) = <[String; 6]>::try_from(written) else {
        anyhow::bail!("rescoring wrote an unexpected number of tables");
    };
    Ok(RescoredTables {
        q,
        q_values,
        p,
        nes,
    })
}

/////////////
// preview //
/////////////

/// How many ranked calls a preview shows per cluster.
const PREVIEW_CALLS: usize = 5;

/// What `decisions` would do to `source`, written nowhere: per cluster the
/// label before and after, and the cell types re-ranked against the edited
/// marker panel. When the pass cached its statistics they are rescored as a
/// save would (only the types the edits touch, when no cluster changes; see
/// [`super::recalibrate::rescore_types`]) and `scores` holds every cluster's
/// share, NES, p and q; else, from the cluster expression profile, ranked
/// approximately (see [`approx_calls`]).
fn preview(source: &Loaded, decisions: Vec<Decision>, decisions_dir: &Path) -> Result<Value> {
    let before = read_cells(source)?;
    let (panel_before, _) = read_markers(source)?;
    let named: BTreeSet<ClusterId> = decisions.iter().flat_map(|d| d.cluster.clone()).collect();
    let after = plan(source, decisions, decisions_dir, "preview")?;

    // Grouped the way the new round would be, so a merged cluster's
    // "before" is what its cells were called before.
    let label_before = digest(&after.cells.clusters, &before.labels, &Evidence::default());
    let label_after = digest(
        &after.cells.clusters,
        &after.cells.labels,
        &Evidence::default(),
    );
    let touched: BTreeSet<ClusterId> = before
        .clusters
        .iter()
        .zip(&after.cells.clusters)
        .filter_map(|(b, a)| b.filter(|b| named.contains(b)).and(*a))
        .collect();

    // The round's recorded calls: the "before" of a rescoring, and the only
    // calls there are when nothing can be rescored.
    let recorded: Calls = digest(&before.clusters, &before.labels, &read_evidence(source)?)
        .into_iter()
        .map(|(id, d)| {
            let calls = d
                .calls
                .into_iter()
                .map(|c| CallView {
                    label: c.label,
                    score: None,
                    q: c.q,
                })
                .collect();
            (id, calls)
        })
        .collect();

    // Best first: rescored with the round's own scoring when the pass cached
    // its statistics, else ranked approximately on the expression profile,
    // else as recorded.
    // With the clusters as they were, only the cell types the marker edits
    // touch are scored again.
    let rescored = if after.cells.clusters == before.clusters {
        let touched = touched_types(&panel_before, &after.markers);
        super::recalibrate::rescore_types(source, &after.cells, &after.markers, &touched)?
    } else {
        super::recalibrate::rescore(source, &after.cells, &after.markers)?
    };
    let scores = rescored.as_ref().map(rescored_scores);
    let (stats, calls_before, calls_after): (&str, Option<Calls>, Calls) =
        if let Some(r) = &rescored {
            ("recalibrated", Some(recorded), rescored_calls(r))
        } else if let Some((table, groups)) = expression_profile(source, &after.cells)? {
            let [b, a] = approx_calls(&table, &groups, [&panel_before, &after.markers])?;
            ("approximate", Some(scored_calls(b)), scored_calls(a))
        } else {
            ("recorded", None, recorded)
        };
    let top = |calls: &Calls, id: &ClusterId| {
        calls
            .get(id)
            .and_then(|c| c.first())
            .map(|c| c.label.clone())
    };

    let mut clusters = serde_json::Map::new();
    for (id, after_d) in &label_after {
        let before_label = label_before.get(id).and_then(|d| d.label.clone());
        let top_after = top(&calls_after, id);
        let top_before = calls_before
            .as_ref()
            .map_or_else(|| top_after.clone(), |c| top(c, id));
        if !(touched.contains(id) || before_label != after_d.label || top_before != top_after) {
            continue;
        }
        let calls: Vec<Value> = calls_after
            .get(id)
            .into_iter()
            .flatten()
            .take(PREVIEW_CALLS)
            .map(|c| json!({"label": c.label, "score": c.score, "q": c.q}))
            .collect();
        clusters.insert(
            id.to_string(),
            json!({
                "label_before": before_label,
                "label_after": after_d.label,
                "top_before": top_before,
                "calls": calls,
            }),
        );
    }
    let cells_changed = before
        .labels
        .iter()
        .zip(&after.cells.labels)
        .filter(|(b, a)| b != a)
        .count();
    Ok(json!({
        "rescored": stats != "recorded",
        "stats": stats,
        "clusters": clusters,
        "cells_changed": cells_changed,
        "markers": marker_diff(&panel_before, &after.markers),
        "scores": scores,
    }))
}

/// The cell types (label keys) whose scores a panel edit changes: those
/// whose markers it edits, and those sharing an edited gene (whose IDF
/// weight moves).
fn touched_types(before: &[(String, String)], after: &[(String, String)]) -> BTreeSet<String> {
    use crate::annotate::markers::label_key;
    let pairs = |p: &[(String, String)]| -> BTreeSet<(String, String)> {
        p.iter().map(|(g, t)| (label_key(t), g.clone())).collect()
    };
    let (b, a) = (pairs(before), pairs(after));
    let edited: BTreeSet<&(String, String)> = b.symmetric_difference(&a).collect();
    let genes: BTreeSet<&str> = edited.iter().map(|(_, g)| g.as_str()).collect();
    edited
        .iter()
        .map(|(t, _)| t.clone())
        .chain(
            a.iter()
                .filter(|(_, g)| genes.contains(g.as_str()))
                .map(|(t, _)| t.clone()),
        )
        .collect()
}

/// A rescored round's every cluster × type, as a preview reports it: per
/// cluster id, each type's share (Q), NES, p and q, highest share first.
fn rescored_scores(r: &super::recalibrate::Rescored) -> Value {
    let finite = |v: f32| v.is_finite().then_some(v);
    r.ids
        .iter()
        .enumerate()
        .map(|(row, id)| {
            let mut types: Vec<usize> = (0..r.types.len()).collect();
            types.sort_by(|&a, &b| {
                r.q_probs[(row, b)]
                    .total_cmp(&r.q_probs[(row, a)])
                    .then(r.p_values[(row, a)].total_cmp(&r.p_values[(row, b)]))
            });
            let calls: Vec<Value> = types
                .into_iter()
                .map(|t| {
                    json!({
                        "label": &*r.types[t],
                        "share": r.q_probs[(row, t)],
                        "nes": finite(r.nes[(row, t)]),
                        "p": finite(r.p_values[(row, t)]),
                        "q": finite(r.q_values[(row, t)]),
                    })
                })
                .collect();
            (id.to_string(), Value::Array(calls))
        })
        .collect::<serde_json::Map<_, _>>()
        .into()
}

/// One ranked call as a preview reports it.
struct CallView {
    label: String,
    /// The approximate module score, when that is how it was ranked.
    score: Option<f32>,
    q: Option<f32>,
}

/// Per cluster, its calls best first.
type Calls = BTreeMap<ClusterId, Vec<CallView>>;

fn scored_calls(ranked: RankedCalls) -> Calls {
    ranked
        .into_iter()
        .map(|(id, v)| {
            let calls = v
                .into_iter()
                .map(|(label, score)| CallView {
                    label,
                    score,
                    q: None,
                })
                .collect();
            (id, calls)
        })
        .collect()
}

/// A rescored round's calls, smallest q first.
fn rescored_calls(r: &super::recalibrate::Rescored) -> Calls {
    r.ids
        .iter()
        .enumerate()
        .map(|(row, id)| {
            let mut calls: Vec<CallView> = r
                .types
                .iter()
                .enumerate()
                .map(|(t, label)| CallView {
                    label: label.to_string(),
                    score: None,
                    q: Some(r.q_values[(row, t)]).filter(|q| q.is_finite()),
                })
                .collect();
            calls.sort_by(|a, b| {
                let q = |c: &CallView| c.q.unwrap_or(f32::INFINITY);
                q(a).total_cmp(&q(b))
            });
            (*id, calls)
        })
        .collect()
}

/// Per cell type, the features the edits add and drop.
fn marker_diff(before: &[(String, String)], after: &[(String, String)]) -> Value {
    let key = |(g, t): &(String, String)| (crate::annotate::markers::label_key(t), g.clone());
    let b: BTreeSet<_> = before.iter().map(key).collect();
    let a: BTreeSet<_> = after.iter().map(key).collect();
    let mut out: BTreeMap<String, (Vec<String>, Vec<String>)> = BTreeMap::new();
    for (t, g) in a.difference(&b) {
        out.entry(t.clone()).or_default().0.push(g.clone());
    }
    for (t, g) in b.difference(&a) {
        out.entry(t.clone()).or_default().1.push(g.clone());
    }
    out.into_iter()
        .map(|(t, (added, dropped))| (t, json!({"added": added, "dropped": dropped})))
        .collect::<serde_json::Map<_, _>>()
        .into()
}

/// The round's gene × cluster expression profile, and for each cluster of
/// `after` how many of its cells come from each of the profile's clusters.
/// `None` when the round records no profile.
type Groups = BTreeMap<ClusterId, BTreeMap<ClusterId, usize>>;
fn expression_profile(
    source: &Loaded,
    after: &Cells,
) -> Result<Option<(MatWithNames<Mat>, Groups)>> {
    let a = &source.manifest.annotate;
    let Some(profile_rel) = a.cluster_expression.as_deref() else {
        return Ok(None);
    };
    let ids_rel = a
        .expression_clusters
        .as_deref()
        .or(source.manifest.cluster.clusters.as_deref())
        .context("no cluster table for the expression profile")?;
    let (names, ids) = read_clusters(&resolve(&source.dir, ids_rel))?;
    let original: HashMap<&str, ClusterId> = names
        .iter()
        .zip(&ids)
        .filter_map(|(n, id)| id.map(|id| (n.as_ref(), id)))
        .collect();
    let mut groups = Groups::new();
    for (cell, id) in after.names.iter().zip(&after.clusters) {
        if let (Some(id), Some(&from)) = (id, original.get(cell.as_ref())) {
            *groups.entry(*id).or_default().entry(from).or_default() += 1;
        }
    }
    let path = resolve(&source.dir, profile_rel);
    let table = Mat::from_parquet_with_row_names(&path, Some(0))
        .with_context(|| format!("reading {path}"))?;
    Ok(Some((table, groups)))
}

/// Per cluster, cell types best first with their score (`None` when the
/// ranking comes from recorded evidence rather than a score).
type RankedCalls = BTreeMap<ClusterId, Vec<(String, Option<f32>)>>;

/// A cell type and its marker rows with their weights.
type TypeMarkers = (String, Vec<(usize, f32)>);

/// Cell types ranked per cluster by a marker module score, for each of
/// `panels`: each cluster's profile is the cell-weighted mean of the profile
/// columns its cells come from; genes are `log1p`, z-scored across clusters;
/// a type scores the IDF-weighted mean z of its markers. Only the panels'
/// marker genes are read, and the profiles are built once for all panels. It
/// ranks as enrichment would, in milliseconds, but carries no q or support.
fn approx_calls<const N: usize>(
    table: &MatWithNames<Mat>,
    groups: &Groups,
    panels: [&[(String, String)]; N],
) -> Result<[RankedCalls; N]> {
    // Per panel: each type's marker rows and weights.
    let mut sparse: Vec<Vec<TypeMarkers>> = Vec::with_capacity(N);
    for panel in panels {
        let types = if panel.is_empty() {
            Vec::new()
        } else {
            let pairs: Vec<(Box<str>, Box<str>)> = panel
                .iter()
                .map(|(g, t)| (g.as_str().into(), t.as_str().into()))
                .collect();
            let annot =
                crate::annotate::markers::annotation_matrix_from_pairs(&pairs, &table.rows)?;
            let w = &annot.membership_ga;
            annot
                .annot_names
                .iter()
                .enumerate()
                .map(|(t, name)| {
                    let rows = (0..w.nrows())
                        .filter(|&r| w[(r, t)] > 0.0)
                        .map(|r| (r, w[(r, t)]))
                        .collect();
                    (name.to_string(), rows)
                })
                .collect()
        };
        sparse.push(types);
    }
    let marker_rows: BTreeSet<usize> = sparse
        .iter()
        .flatten()
        .flat_map(|(_, rows)| rows.iter().map(|&(r, _)| r))
        .collect();

    // Per cluster: log1p of its cell-weighted mean profile, marker rows only.
    let col: HashMap<ClusterId, usize> = table
        .cols
        .iter()
        .enumerate()
        .filter_map(|(j, c)| parse_cluster_id(c).map(|id| (id, j)))
        .collect();
    let profiles: Vec<(ClusterId, HashMap<usize, f32>)> = groups
        .iter()
        .filter_map(|(id, from)| {
            let parts: Vec<(usize, f32)> = from
                .iter()
                .filter_map(|(f, n)| col.get(f).map(|&j| (j, *n as f32)))
                .collect();
            let total: f32 = parts.iter().map(|(_, n)| n).sum();
            (total > 0.0).then(|| {
                let v = marker_rows
                    .iter()
                    .map(|&r| {
                        let m: f32 = parts.iter().map(|&(j, n)| n * table.mat[(r, j)]).sum();
                        (r, (m / total).max(0.0).ln_1p())
                    })
                    .collect();
                (*id, v)
            })
        })
        .collect();
    // Per marker row: its z-score in each cluster.
    let k = profiles.len().max(1) as f32;
    let z: HashMap<usize, (f32, f32)> = marker_rows
        .iter()
        .map(|&r| {
            let mean = profiles.iter().map(|(_, v)| v[&r]).sum::<f32>() / k;
            let var = profiles
                .iter()
                .map(|(_, v)| (v[&r] - mean).powi(2))
                .sum::<f32>()
                / k;
            (r, (mean, var.sqrt()))
        })
        .collect();

    Ok(std::array::from_fn(|p| {
        profiles
            .iter()
            .map(|(id, v)| {
                let mut scored: Vec<(String, Option<f32>)> = sparse[p]
                    .iter()
                    .filter_map(|(name, rows)| {
                        let (num, den) = rows.iter().fold((0.0f32, 0.0f32), |(n, d), &(r, wt)| {
                            let (mean, sd) = z[&r];
                            let zr = if sd > 0.0 { (v[&r] - mean) / sd } else { 0.0 };
                            (n + wt * zr, d + wt)
                        });
                        (den > 0.0).then(|| (name.clone(), Some(num / den)))
                    })
                    .collect();
                scored.sort_by(|a, b| b.1.unwrap_or(0.0).total_cmp(&a.1.unwrap_or(0.0)));
                (*id, scored)
            })
            .collect()
    }))
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
        let _ = writeln!(s, "  call  {:<24} q={}", c.label, num(c.q));
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
