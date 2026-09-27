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
use crate::manifest::run::{self, annotated_path, rel_to_manifest, resolve, Loaded};
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

#[derive(Args, Debug)]
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
        help = "Decisions, one JSON object per line (see `lupin review --help`)"
    )]
    pub decisions: Box<str>,

    #[arg(
        long,
        short = 'o',
        help = "Output prefix for the new round (with --watch: rounds are {out}.r1, {out}.r2, ...)"
    )]
    pub out: Box<str>,

    #[arg(
        long,
        help = "Keep running: each batch of lines appended to the decisions file becomes the next round"
    )]
    pub watch: bool,
}

/// Apply a decisions file to a round and write the next one, or with
/// `--watch`, keep applying whatever is appended to it.
pub fn run_relabel(args: &RelabelArgs) -> Result<()> {
    if args.watch {
        return watch(args);
    }
    let source = run::load(&args.from)?;
    let raw = fs::read_to_string(&*args.decisions)
        .with_context(|| format!("reading {}", args.decisions))?;
    let lines: Vec<(usize, &str)> = raw.lines().enumerate().map(|(n, l)| (n + 1, l)).collect();
    let decisions = parse_decisions(&lines, &args.decisions)?;
    relabel(
        &source,
        decisions,
        &parent_dir(Path::new(&*args.decisions)),
        &args.out,
    )?;
    Ok(())
}

/// `(line number, line)` pairs to decisions; blank lines are skipped.
fn parse_decisions(lines: &[(usize, &str)], file: &str) -> Result<Vec<Decision>> {
    let decisions: Vec<Decision> = lines
        .iter()
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

/// Apply `decisions` to the round `source` and write the next round at
/// `out`; returns its manifest. Nothing is written if a decision is refused.
fn relabel(
    source: &Loaded,
    mut decisions: Vec<Decision>,
    decisions_dir: &Path,
    out: &str,
) -> Result<PathBuf> {
    check_round(source, &decisions, decisions_dir)?;
    mkdir_parent(out)?;
    let manifest_path = annotated_path(&source.file, out);
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

/// Refuse decisions that name a round other than `source`.
fn check_round(source: &Loaded, decisions: &[Decision], decisions_dir: &Path) -> Result<()> {
    let canon = |p: &Path| p.canonicalize().unwrap_or_else(|_| p.to_path_buf());
    let current = canon(&source.file);
    for d in decisions {
        let Some(named) = &d.round else { continue };
        let named_path = PathBuf::from(resolve(decisions_dir, named));
        anyhow::ensure!(
            canon(&named_path) == current,
            "ids refer to {named}, latest is {}; review the latest round and decide again",
            source.file.display()
        );
    }
    Ok(())
}

///////////
// watch //
///////////

pub const STATUS: &str = ".relabel_status.json";
const POLL: std::time::Duration = std::time::Duration::from_millis(500);

/// What `relabel --watch` has done, rewritten after every batch at
/// `{out}.relabel_status.json`. Paths are relative to that file.
#[derive(Serialize, Deserialize, Debug, Default)]
pub struct WatchStatus {
    /// The decisions file being watched.
    pub decisions: String,
    /// Lines of it already handled, applied or refused.
    pub processed_lines: usize,
    /// The round the watcher started from (`-f`).
    #[serde(default)]
    pub base: String,
    /// Every round in order, starting with `base`.
    pub rounds: Vec<String>,
    /// The round the next batch applies to.
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
    let out = args.out.to_string();
    let (mut status, status_path) = begin_watch(args)?;
    info!(
        "watching {} (from line {}); rounds go to {out}.r<N>; status in {}",
        args.decisions,
        status.processed_lines + 1,
        status_path.display()
    );
    loop {
        watch_step(&mut status, &status_path, &out)?;
        std::thread::sleep(POLL);
    }
}

/// Start or resume a watch and write its status before any decision, so a
/// viewer can find the watcher at once.
fn begin_watch(args: &RelabelArgs) -> Result<(WatchStatus, PathBuf)> {
    mkdir_parent(&args.out)?;
    let status_path = PathBuf::from(format!("{}{STATUS}", args.out));
    let mut status = start_watch(&status_path, &args.from, &args.decisions)?;
    write_status(&mut status, &status_path)?;
    Ok((status, status_path))
}

/// Resume from an existing status for the same decisions file, else start
/// at `from` with nothing processed.
fn start_watch(status_path: &Path, from: &str, decisions: &str) -> Result<WatchStatus> {
    let dir = parent_dir(status_path);
    let decisions_rel = rel_to_manifest(&dir, decisions);
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
    let first = run::load(from)?;
    let base = rel_to_manifest(&dir, &first.file.to_string_lossy());
    Ok(WatchStatus {
        decisions: decisions_rel,
        rounds: vec![base.clone()],
        latest: base.clone(),
        base,
        updated: now(),
        ..WatchStatus::default()
    })
}

fn parent_dir(p: &Path) -> PathBuf {
    p.parent()
        .filter(|d| !d.as_os_str().is_empty())
        .map_or_else(|| PathBuf::from("."), Path::to_path_buf)
}

/// Apply the complete lines appended since the last step as one round.
/// Returns whether there was anything new. A refused batch is recorded in
/// the status and skipped, so the next append is not stuck behind it.
pub fn watch_step(status: &mut WatchStatus, status_path: &Path, out: &str) -> Result<bool> {
    let dir = parent_dir(status_path);
    let decisions = resolve(&dir, &status.decisions);
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
    let batch: Vec<(usize, &str)> = raw
        .lines()
        .enumerate()
        .skip(status.processed_lines)
        .take(complete - status.processed_lines)
        .map(|(n, l)| (n + 1, l))
        .collect();
    let first_line = status.processed_lines + 1;
    status.processed_lines = complete;
    if batch.iter().all(|(_, l)| l.trim().is_empty()) {
        return Ok(false);
    }

    // `rounds` starts with the base, so its length numbers the next round.
    let round_prefix = format!("{out}.r{}", status.rounds.len().max(1));
    let result = run::load(&resolve(&dir, &status.latest)).and_then(|source| {
        let ds = parse_decisions(&batch, &decisions)?;
        relabel(
            &source,
            ds,
            &parent_dir(Path::new(&decisions)),
            &round_prefix,
        )
    });
    match result {
        Ok(file) => {
            let rel = rel_to_manifest(&dir, &file.to_string_lossy());
            info!(
                "round {} written: {}",
                status.rounds.len() + 1,
                file.display()
            );
            status.rounds.push(rel.clone());
            status.latest = rel;
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
    write_status(status, status_path)?;
    Ok(true)
}

/// Stamp and write the status whole, then rename it into place, so a reader
/// never sees half a file.
fn write_status(status: &mut WatchStatus, status_path: &Path) -> Result<()> {
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
            "  hist  {} {:?} -> {label} by {:?}: {}",
            h.round, h.action, h.decided_by, h.rationale
        );
    }
    s
}

#[cfg(test)]
#[path = "rounds_tests.rs"]
mod tests;
