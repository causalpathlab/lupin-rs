//! Manifest glue for annotation rounds ([`crate::annotate::rounds`]): read a
//! round's tables, write the next round, and the `review` / `relabel`
//! commands.
//!
//! A round writes, beside its manifest:
//! - `{out}.clusters.parquet`: per-cell cluster id (`cluster.clusters`)
//! - `{out}.argmax.parquet`: per-cell label (`annotate.argmax`; older rounds
//!   record the TSV), with `{out}.argmax.tsv` and `{out}.membership.tsv`
//!   beside it
//! - `{out}.cluster_summary.json`: the digest (`annotate.cluster_summary`)
//! - `{out}.cluster_evidence.parquet`: every cluster's types with every
//!   method's statistic, one row each (`annotate.cluster_evidence`)
//! - `{out}.annotation_log.jsonl`: this round's decisions (`annotate.log`)
//! - `{out}.annotation_history.json`: every round's decisions per cluster,
//!   newest first (`annotate.history`)
//! - when a round edits markers, `{out}.markers.tsv` (`annotate.markers`, the
//!   previous round's panel with the edits applied) and
//!   `{out}.marker_history.json` (`annotate.marker_history`, per cell type)

use crate::annotate::celltype_tree::{TreePlacement, TypeTree, Via};
use crate::annotate::markers::label_key;
use crate::annotate::outputs::{
    write_cluster_tables, CLUSTER_CELLTYPE_ES_STD, CLUSTER_CELLTYPE_NES, CLUSTER_CELLTYPE_P,
    CLUSTER_CELLTYPE_Q, CLUSTER_CELLTYPE_Q_VALUES, CLUSTER_CELLTYPE_Z,
};
use crate::annotate::rounds::{
    self, apply, apply_markers, digest, next_cluster_id, parse_cluster_id, prepend_history,
    ClPlacement, ClusterId, Decision, Digest, Evidence, History, MarkerHistory, Method, Stat,
    Table, Term, TieBreak, TypeEvidence, TypeTables,
};
use crate::manifest::recalibrate::Terms;
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
/// `annotate.stats`' list of what a round's readers should know of its
/// statistics.
pub(crate) const CAVEATS: &str = "caveats";
/// The caveat of a merge whose GO/GMT terms could not be scored again.
pub(crate) const TERMS_NOT_RESCORED: &str = "the GO/GMT terms were not scored again: the merged \
     clusters have none until a new pass scores them; the other clusters keep theirs";
pub const SUMMARY: &str = ".cluster_summary.json";
pub const EVIDENCE: &str = ".cluster_evidence.parquet";
pub const LOG: &str = ".annotation_log.jsonl";
pub const HISTORY: &str = ".annotation_history.json";
pub const MARKERS: &str = ".markers.tsv";
pub const MARKER_HISTORY: &str = ".marker_history.json";

//////////////////
// cell tables  //
//////////////////

/// Cell names and, per cell, its cluster id (`None` when unassigned).
pub type CellClusters = (Vec<Box<str>>, Vec<Option<ClusterId>>);

/// Cluster labels (cell → cluster, `usize::MAX` for none) as the ids
/// [`write_clusters`] takes.
#[must_use]
pub fn cluster_ids(labels: &[usize]) -> Vec<Option<ClusterId>> {
    labels
        .iter()
        .map(|&k| ClusterId::try_from(k).ok())
        .collect()
}

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

/// A round's per-cell labels, `cell`, `cell_type` and (when present)
/// `probability`, from its parquet or its header-first TSV.
pub(crate) fn read_argmax(path: &str) -> Result<HashMap<String, (String, f32)>> {
    Ok(crate::cell_labels::read_label_rows(path)?
        .into_iter()
        .map(|(cell, label, p)| (cell.into(), (label.into(), p)))
        .collect())
}

/// Write a round's per-cell labels under `out` through the shared writer
/// (`{out}.argmax.parquet`, the TSVs beside it); the file to record,
/// [`crate::annotate::outputs::labels_path`].
pub(super) fn write_argmax(
    out: &str,
    cells: &[Box<str>],
    labels: &[Option<String>],
    probs: &[f32],
) -> Result<String> {
    let labels: Vec<Box<str>> = labels
        .iter()
        .map(|l| l.as_deref().unwrap_or(UNASSIGNED_LABEL).into())
        .collect();
    graph_embedding_util::type_annotation::write_cell_labels(out, cells, &labels, probs)?;
    Ok(crate::annotate::outputs::labels_path(out))
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
    read_terms(&resolve(&loaded.dir, rel), term_source(loaded))
}

/// Where a round's terms come from: `gmt` when its pass read a GMT, else `go`.
fn term_source(loaded: &Loaded) -> &'static str {
    match (loaded.manifest.annotate.settings.as_ref()).and_then(|s| s.pointer("/enrichment/gmt")) {
        Some(g) if !g.is_null() => "gmt",
        _ => "go",
    }
}

/// Whatever evidence this round records; a missing table is skipped.
fn read_evidence(loaded: &Loaded) -> Result<Evidence<'static>> {
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

/// A per-cell reading (the projection score's, the gene-set units'),
/// aligned to a round's cells.
#[derive(Default)]
pub(crate) struct ScoreCells {
    /// The types it called, in order of first call.
    pub(crate) types: Vec<Box<str>>,
    /// Per cell: its call, an index into `types`; `None` for no row or no call.
    pub(crate) call: Vec<Option<usize>>,
    /// Per cell: the score's lead of its call over the runner-up, or the
    /// cell's membership of the units' region it is called by; NaN for no
    /// row.
    pub(crate) margin: Vec<f32>,
}

impl ScoreCells {
    /// Per cell, the type called.
    pub(crate) fn calls(&self) -> Vec<Option<&str>> {
        self.call
            .iter()
            .map(|c| c.map(|k| self.types[k].as_ref()))
            .collect()
    }

    /// Each of `cells`' current clusters' share of its cells called each type.
    pub(crate) fn shares(&self, cells: &Cells) -> Table {
        rounds::call_shares(&cells.clusters, &self.calls())
    }
}

/// The score's per-cell table the round records, aligned to `cells`; `None`
/// when it records none or the file is gone.
pub(crate) fn read_score_cells(loaded: &Loaded, cells: &Cells) -> Result<Option<ScoreCells>> {
    let rel = loaded.manifest.annotate.marker_score.as_deref();
    read_call_cells(loaded, rel, "margin", cells)
}

/// The units' per-cell table the round records (`annotate.units_call`),
/// aligned to `cells`; `None` when it records none or the file is gone.
pub(crate) fn read_units_cells(loaded: &Loaded, cells: &Cells) -> Result<Option<ScoreCells>> {
    let rel = loaded.manifest.annotate.units_call.as_deref();
    read_call_cells(loaded, rel, "membership", cells)
}

/// A per-cell table (`cell`, `call` and the number `column`) at `rel`,
/// aligned to `cells`.
fn read_call_cells(
    loaded: &Loaded,
    rel: Option<&str>,
    column: &str,
    cells: &Cells,
) -> Result<Option<ScoreCells>> {
    use legume_numeric::matrix::parquet::read_table_columns;
    let Some(path) = rel.and_then(|rel| loaded.recorded(rel)) else {
        return Ok(None);
    };
    let (strings, nums) = read_table_columns(&path, &["cell", "call"], &[column])
        .with_context(|| format!("reading {path}"))?;
    let row: HashMap<&str, usize> = strings[0]
        .iter()
        .enumerate()
        .map(|(i, c)| (c.as_ref(), i))
        .collect();
    let mut out = ScoreCells::default();
    let mut code: HashMap<&str, usize> = HashMap::new();
    for c in &cells.names {
        let Some(&i) = row.get(c.as_ref()) else {
            out.call.push(None);
            out.margin.push(f32::NAN);
            continue;
        };
        // In the panel's own spelling, as the other methods' tables are.
        let call = strings[1][i].as_ref();
        out.call.push((call != UNASSIGNED_LABEL).then(|| {
            *code.entry(call).or_insert_with(|| {
                out.types.push(label_key(call).into());
                out.types.len() - 1
            })
        }));
        out.margin.push(nums[0][i] as f32);
    }
    Ok(Some(out))
}

/// The per-cell readings' shares of each of `cells`' clusters
/// ([`ScoreCells::shares`]): the score's and the units', where recorded.
#[derive(Default)]
pub(crate) struct CellShares {
    pub(crate) score: Option<Table>,
    pub(crate) units: Option<Table>,
}

impl CellShares {
    /// Each per-cell reading `loaded` records, tabled over `cells`'
    /// clusters.
    pub(crate) fn read(loaded: &Loaded, cells: &Cells) -> Result<Self> {
        Ok(Self {
            score: read_score_cells(loaded, cells)?.map(|s| s.shares(cells)),
            units: read_units_cells(loaded, cells)?.map(|s| s.shares(cells)),
        })
    }
}

/// Every cluster × type table the round records, read once: the
/// enrichment's and SuSiE's, with the per-cell readings' `shares` beside
/// them. A table recorded but gone from disk is skipped.
pub(crate) fn read_type_tables(loaded: &Loaded, mut shares: CellShares) -> Result<TypeTables> {
    let a = &loaded.manifest.annotate;
    let mut tables = Vec::new();
    for stat in Stat::ALL {
        let rel = match stat {
            Stat::Share => &a.cluster_celltype_q,
            Stat::Nes => &a.cluster_celltype_nes,
            Stat::P => &a.cluster_celltype_p,
            Stat::Q => &a.cluster_celltype_q_values,
            Stat::Pip => &a.cluster_celltype_pip,
            Stat::Fold => &a.cluster_celltype_effect,
            Stat::Explained => &a.cluster_celltype_explained,
            Stat::Score => {
                tables.extend(shares.score.take().map(|t| (stat, t)));
                continue;
            }
            Stat::Units => {
                tables.extend(shares.units.take().map(|t| (stat, t)));
                continue;
            }
        };
        match rel.as_deref().map(|r| (r, loaded.recorded(r))) {
            Some((_, Some(path))) => tables.push((stat, read_table(&path)?)),
            Some((r, None)) => log::warn!("{r} is recorded but gone: no {} evidence", stat.name()),
            None => {}
        }
    }
    Ok(TypeTables { tables })
}

/// The method whose reading made the round's call, as the round records it
/// (`annotate.called_by`); for a round from before that, by the tables it
/// records: SuSiE's, else the enrichment's share, else the score's per-cell
/// table, else the units'. `None` for a round with none of them.
pub(crate) fn called_by(a: &super::run::RunAnnotate) -> Option<Method> {
    if a.called_by.is_some() {
        return a.called_by;
    }
    let susie = a.cluster_celltype_pip.is_some()
        && a.cluster_celltype_effect.is_some()
        && a.cluster_celltype_explained.is_some();
    if susie {
        Some(Method::Susie)
    } else if a.cluster_celltype_q.is_some() {
        Some(Method::Enrichment)
    } else {
        (a.marker_score.as_ref().map(|_| Method::Score))
            .or_else(|| a.units_call.as_ref().map(|_| Method::Units))
    }
}

/// Every cluster's types with every method's evidence, cluster by cluster.
#[cfg(test)]
pub(crate) fn type_evidence(loaded: &Loaded, cells: &Cells) -> Result<Vec<TypeEvidence>> {
    Ok(evidence_rows(
        &tie_break(loaded, Some(cells))?.tables,
        cells,
    ))
}

/// `tables`' evidence for every cluster of `cells`, cluster by cluster.
fn evidence_rows(tables: &TypeTables, cells: &Cells) -> Vec<TypeEvidence> {
    let ids: BTreeSet<ClusterId> = cells.clusters.iter().flatten().copied().collect();
    ids.into_iter().flat_map(|id| tables.evidence(id)).collect()
}

/// Write `rows` as the long table `{out}.cluster_evidence.parquet`: cluster
/// (`K{id}`), celltype, then each [`Stat`], NaN where it has none.
fn write_type_evidence(path: &str, rows: &[TypeEvidence]) -> Result<()> {
    use legume_numeric::matrix::parquet::{write_table, Column};
    let cluster: Vec<Box<str>> = rows
        .iter()
        .map(|r| crate::figures::cluster_name(r.cluster).into())
        .collect();
    let celltype: Vec<Box<str>> = rows.iter().map(|r| r.label.as_str().into()).collect();
    let stats: Vec<Vec<f32>> = Stat::ALL
        .iter()
        .map(|&s| rows.iter().map(|r| r.get(s).unwrap_or(f32::NAN)).collect())
        .collect();
    let mut columns: Vec<(Box<str>, Column)> = vec![
        ("cluster".into(), Column::Str(&cluster)),
        ("celltype".into(), Column::Str(&celltype)),
    ];
    columns.extend(
        Stat::ALL
            .iter()
            .zip(&stats)
            .map(|(s, v)| (s.name().into(), Column::F32(v))),
    );
    write_table(path, &columns)?;
    info!("wrote {path}");
    Ok(())
}

/// What orders the round's calls tied on q: its type tables (with the
/// per-cell readings' shares over `cells`' clusters, when given) under the
/// method that made its call.
fn tie_break(loaded: &Loaded, cells: Option<&Cells>) -> Result<TieBreak> {
    let shares = match cells {
        Some(c) => CellShares::read(loaded, c)?,
        None => CellShares::default(),
    };
    Ok(TieBreak {
        tables: read_type_tables(loaded, shares)?,
        method: called_by(&loaded.manifest.annotate).unwrap_or_default(),
    })
}

/// The round's type tree (`annotate.celltype_tree`), when it records one
/// that is still there.
fn read_tree(loaded: &Loaded) -> Result<Option<TypeTree>> {
    let Some(path) =
        (loaded.manifest.annotate.celltype_tree.as_deref()).and_then(|rel| loaded.recorded(rel))
    else {
        return Ok(None);
    };
    let text = fs::read_to_string(&path).with_context(|| format!("reading {path}"))?;
    Ok(Some(serde_json::from_str(&text)?))
}

/// Build this round's digest, write `{out}.cluster_summary.json` and record
/// it in the manifest (which the caller saves).
pub fn write_summary(loaded: &mut Loaded, out_prefix: &str) -> Result<BTreeMap<ClusterId, Digest>> {
    let cells = read_cells(loaded)?;
    // A table that cannot be read costs the order of tied calls, never the
    // summary.
    let ties = tie_break(loaded, Some(&cells));
    if let Err(e) = &ties {
        log::warn!("tied calls by name: {e:#}");
    }
    let mut ev = read_evidence(loaded)?;
    ev.ties = ties.as_ref().ok();
    let mut d = digest(&cells.clusters, &cells.labels, &ev);
    // A cluster labelled with a coarse group agrees with a top call that is
    // one of the group's types.
    if let Some(tree) = read_tree(loaded)? {
        let index = tree.index();
        for digest in d.values_mut() {
            let label = digest.label.clone();
            let group_agrees = label.is_some_and(|label| {
                digest.top_calls().iter().any(|c| {
                    let key = label_key(&c.label);
                    index.get(&key).is_some_and(|g| *g == label)
                })
            });
            if let Some(e) = digest.evidence.as_mut() {
                e.agrees |= group_agrees;
            }
        }
    }
    let path = format!("{out_prefix}{SUMMARY}");
    fs::write(&path, serde_json::to_string_pretty(&d)?)?;
    info!("wrote {path}");
    loaded.manifest.annotate.cluster_summary = Some(rel_to_manifest(&loaded.dir, &path));
    let evidence = format!("{out_prefix}{EVIDENCE}");
    // A view of the tables just recorded: one that cannot be read costs the
    // evidence table, never the round; the parent's is never kept.
    loaded.manifest.annotate.cluster_evidence = None;
    let written = match &ties {
        Ok(t) => write_type_evidence(&evidence, &evidence_rows(&t.tables, &cells)),
        Err(e) => Err(anyhow::anyhow!("{e:#}")),
    };
    match written {
        Ok(()) => {
            loaded.manifest.annotate.cluster_evidence =
                Some(rel_to_manifest(&loaded.dir, &evidence));
        }
        Err(e) => log::warn!("no evidence table for this round: {e:#}"),
    }
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
        requires = "next",
        value_name = "PREFIX",
        help = "With --next: refuse unless the next round is PREFIX ({chain}.r<N+1>), so a \
                replayed step writes the round it wrote before"
    )]
    pub expect: Option<Box<str>>,

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
            if let Some(want) = &args.expect {
                let next = chain.next_prefix();
                anyhow::ensure!(
                    std::path::absolute(&**want)? == std::path::absolute(&next)?,
                    "the next round of the chain is {next}, not {want}"
                );
            }
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

    /// The prefix `--next` would write now: one past the latest round.
    pub(crate) fn next_prefix(&self) -> String {
        self.round_prefix(self.latest().0 + 1)
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
        "{} is not the latest round (that is {}): review that one and decide again",
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
    if super::recalibrate::rescores_by_contrast(source, &markers)
        && decisions.iter().any(|d| d.action.regroups())
    {
        let mut left = cells.clusters.iter().flatten();
        let first = left.next();
        anyhow::ensure!(
            left.any(|c| Some(c) != first),
            "the decisions leave one cluster: each cluster is scored against the others, so \
             nothing would be; to call them one type, label each cluster with it instead"
        );
    }
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
    let argmax_path = write_argmax(out, &p.cells.names, &p.cells.labels, &p.cells.probs)?;
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

    // Rescore the new clusters and panel when the pass cached its statistics;
    // the GO/GMT terms too when the clusters changed (they do not read the
    // panel).
    let regrouped = p.decisions.iter().any(|d| d.action.regroups());
    let wanted = if regrouped {
        Terms::Tested
    } else {
        Terms::None
    };
    let (rescored, terms) =
        match super::recalibrate::rescore_regrouped(source, &p.cells, &p.markers, wanted)? {
            Some((r, t)) => {
                let terms = t.map(|t| t.write(out, &r.row_names())).transpose()?;
                (Some(r), terms)
            }
            None => (None, None),
        };
    let unscored =
        regrouped && terms.is_none() && source.manifest.annotate.ontology_signature.is_some();
    if unscored {
        log::warn!("{TERMS_NOT_RESCORED}");
    }
    let tables = rescored
        .as_ref()
        .map(|r| write_rescored(r, out))
        .transpose()?;

    let mut next = source.copy_to(manifest_path)?;
    // SuSiE's tables describe the source round's clusters on its panel: kept
    // while both stand, dropped once a decision regroups the clusters or the
    // markers change (the call is then the decisions over the enrichment's
    // rescored statistics).
    if regrouped || p.marker_history.is_some() {
        next.manifest.annotate.drop_susie_tables();
    }
    // The source's terms describe its clusters, and a merge's own under
    // fresh ids: rescored, they describe every new cluster.
    if let Some((sig, effect)) = &terms {
        let a = &mut next.manifest.annotate;
        a.ontology_signature = Some(rel_to_manifest(&next.dir, sig));
        a.ontology_term_effect = Some(rel_to_manifest(&next.dir, effect));
    }
    // The score's and the units' per-cell calls hold over new clusters, but
    // not over a new panel: they were made on the old one.
    if p.marker_history.is_some() {
        next.manifest.annotate.drop_marker_score();
        next.manifest.annotate.drop_units_call();
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
    // Terms missing for merged clusters stay missing until a round scores
    // them all again.
    let stale =
        terms.is_none() && (unscored || caveats(source).iter().any(|c| c == TERMS_NOT_RESCORED));
    let stats = &mut next.manifest.annotate.stats;
    if stale {
        let s = stats.get_or_insert_with(|| json!({}));
        if !s.is_object() {
            *s = json!({ "recorded": s.take() });
        }
        s[CAVEATS] = json!([TERMS_NOT_RESCORED]);
    } else if let Some(s) = stats.as_mut().and_then(Value::as_object_mut) {
        s.remove(CAVEATS);
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
/// How many of a regrouped cluster's GO/GMT terms a preview shows.
const PREVIEW_TERMS: usize = 5;

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
    // New clusters' GO/GMT terms are scored too, so a merge is weighed with
    // what the merged cluster's genes do.
    let (rescored, terms) = if after.cells.clusters == before.clusters {
        let touched = touched_types(&panel_before, &after.markers);
        let r = super::recalibrate::rescore_types(source, &after.cells, &after.markers, &touched)?;
        (r, BTreeMap::new())
    } else {
        let rescore = super::recalibrate::rescore_regrouped;
        match rescore(source, &after.cells, &after.markers, Terms::Effects)? {
            Some((r, t)) => {
                // The clusters a decision made: the others keep the round's terms.
                let made = |id: &ClusterId| !before.clusters.contains(&Some(*id));
                let terms: BTreeMap<ClusterId, Vec<Term>> = t
                    .map(|t| {
                        let src = term_source(source);
                        (r.ids.iter().enumerate())
                            .filter(|(_, id)| made(id))
                            .map(|(k, id)| (*id, t.top_terms(k, src)))
                            .collect()
                    })
                    .unwrap_or_default();
                (Some(r), terms)
            }
            None => (None, BTreeMap::new()),
        }
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
        let mut entry = json!({
            "label_before": before_label,
            "label_after": after_d.label,
            "top_before": top_before,
            "calls": calls,
        });
        if let Some(t) = terms.get(id) {
            entry["terms"] = serde_json::to_value(&t[..t.len().min(PREVIEW_TERMS)])?;
        }
        clusters.insert(id.to_string(), entry);
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
    let key = |(g, t): &(String, String)| (label_key(t), g.clone());
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

/// The round's genes × clusters expression profile over its own clusters:
/// the pass's profile ([`expression_profile`]), each of the round's clusters
/// the cell-weighted mean of the pass clusters its cells came from, so a
/// merge's profile is its parts'. `None` when the round records no profile.
pub(crate) fn round_profile(loaded: &Loaded, cells: &Cells) -> Result<Option<MatWithNames<Mat>>> {
    let Some((table, groups)) = expression_profile(loaded, cells)? else {
        return Ok(None);
    };
    Ok(Some(regrouped_profile(&table, &groups)))
}

/// `table`'s columns (`K{id}`, pass clusters) mixed into `groups`' clusters,
/// each the mean of its pass clusters weighted by the cells it took from
/// each; a cluster none of whose pass clusters is in `table` is left out.
fn regrouped_profile(table: &MatWithNames<Mat>, groups: &Groups) -> MatWithNames<Mat> {
    let col: HashMap<ClusterId, usize> = (table.cols.iter().enumerate())
        .filter_map(|(j, c)| parse_cluster_id(c).map(|id| (id, j)))
        .collect();
    let parts: Vec<(ClusterId, Vec<(usize, f32)>)> = groups
        .iter()
        .filter_map(|(id, from)| {
            let parts: Vec<(usize, f32)> = (from.iter())
                .filter_map(|(f, n)| col.get(f).map(|&j| (j, *n as f32)))
                .collect();
            (!parts.is_empty()).then_some((*id, parts))
        })
        .collect();
    let mut mat = Mat::zeros(table.mat.nrows(), parts.len());
    for (k, (_, from)) in parts.iter().enumerate() {
        let total: f32 = from.iter().map(|(_, n)| n).sum();
        for &(j, n) in from {
            mat.column_mut(k).axpy(n / total, &table.mat.column(j), 1.0);
        }
    }
    MatWithNames {
        rows: table.rows.clone(),
        cols: parts
            .iter()
            .map(|(id, _)| format!("K{id}").into())
            .collect(),
        mat,
    }
}
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
    let merged = regrouped_profile(table, groups);
    let profiles: Vec<(ClusterId, HashMap<usize, f32>)> = (merged.cols.iter().enumerate())
        .filter_map(|(k, c)| {
            let v = (marker_rows.iter())
                .map(|&r| (r, merged.mat[(r, k)].max(0.0).ln_1p()))
                .collect();
            Some((parse_cluster_id(c)?, v))
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
            "ids refer to {named}, but the round being relabelled is {}: review that one and decide again",
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
        help = "Print JSON (per cluster: digest, history, each method's top type, every type's evidence) instead of text"
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
decided_by: user | agent_proposed_user_accepted | user_override | agent
round:      optional; the round whose ids the decision names, relative to the
            decisions file. A decision on any other round is refused.
Every decision needs a rationale; it is kept in the round's history.
With several writers (a viewer and an agent), append each decision as one
write of one newline-terminated line under 4 KB, so lines never interleave.";

/// The round's digests: its summary when written, else built from `cells`
/// (read when needed and not given); either way their calls tied on q
/// settled by `ties`.
fn read_digests(
    loaded: &Loaded,
    cells: Option<&Cells>,
    ties: &TieBreak,
) -> Result<BTreeMap<ClusterId, Digest>> {
    if let Some(path) =
        (loaded.manifest.annotate.cluster_summary.as_deref()).and_then(|rel| loaded.recorded(rel))
    {
        let mut d: BTreeMap<ClusterId, Digest> = serde_json::from_str(&fs::read_to_string(path)?)?;
        for (id, digest) in &mut d {
            digest.settle(*id, ties);
        }
        return Ok(d);
    }
    let read;
    let cells = match cells {
        Some(c) => c,
        None => {
            read = read_cells(loaded)?;
            &read
        }
    };
    let mut ev = read_evidence(loaded)?;
    ev.ties = Some(ties);
    Ok(digest(&cells.clusters, &cells.labels, &ev))
}

/// What `review` shows of a round: per cluster its digest, its history and
/// its place on the round's type tree; and the type tables.
struct Reviewed {
    digests: BTreeMap<ClusterId, Digest>,
    history: History,
    tree: Option<TypeTree>,
    tables: TypeTables,
}

impl Reviewed {
    /// The round as review shows it; `every_type` when the types' evidence
    /// is wanted too, `tops` when each method's top type is.
    fn read(loaded: &Loaded, every_type: bool, tops: bool) -> Result<Self> {
        // The cells give the score's and the units' shares over the current
        // clusters: for the types' evidence, the digests when no summary was
        // written, the order of tied calls when either made the call, and
        // their tops. A round whose cells cannot be read is still reviewed
        // from its summary.
        let a = &loaded.manifest.annotate;
        let summary = (a.cluster_summary.as_deref()).and_then(|rel| loaded.recorded(rel));
        let by_cells = matches!(called_by(a), Some(Method::Score | Method::Units));
        let cell_tops = tops && (a.marker_score.is_some() || a.units_call.is_some());
        let wants_cells = every_type || summary.is_none() || by_cells || cell_tops;
        let cells = wants_cells
            .then(|| {
                read_cells(loaded)
                    .inspect_err(|e| log::warn!("review without the cells: {e:#}"))
                    .ok()
            })
            .flatten();
        // A table or tree that cannot be read costs what it adds (the order
        // of tied calls, the types' evidence, the tree's placement), never
        // the review.
        let ties = tie_break(loaded, cells.as_ref()).unwrap_or_else(|e| {
            log::warn!("review without the type tables: {e:#}");
            TieBreak {
                tables: TypeTables::default(),
                method: called_by(&loaded.manifest.annotate).unwrap_or_default(),
            }
        });
        let tree = read_tree(loaded)
            .inspect_err(|e| log::warn!("review without the type tree: {e:#}"))
            .ok()
            .flatten();
        Ok(Self {
            digests: read_digests(loaded, cells.as_ref(), &ties)?,
            history: read_history(loaded)?,
            tree,
            tables: ties.tables,
        })
    }

    /// The Cell Ontology term the tree puts cluster `d`'s label on.
    fn placement(&self, d: &Digest) -> Option<TreePlacement> {
        self.tree.as_ref()?.placement(d.label.as_deref()?)
    }
}

/// `review --json`: per cluster (all, or those in `clusters`), its digest,
/// its history, each method's top type (`tops`, the first of
/// [`Method::ranked`], for the methods that put one forward), every type's
/// evidence (`types`, when `with_types`) and, when the round has a
/// type tree, the Cell Ontology term it puts the label on (`cl`).
pub(crate) fn review_json(
    loaded: &Loaded,
    clusters: &[ClusterId],
    with_types: bool,
) -> Result<BTreeMap<ClusterId, Value>> {
    let r = Reviewed::read(loaded, with_types, true)?;
    let (digests, history, tables) = (&r.digests, &r.history, &r.tables);
    Ok(digests
        .iter()
        .filter(|(id, _)| clusters.is_empty() || clusters.contains(id))
        .map(|(id, d)| {
            let types = tables.evidence(*id);
            let tops: serde_json::Map<String, Value> = Method::ALL
                .into_iter()
                .filter_map(|m| {
                    let top = m.ranked(&types).first()?.label.clone();
                    Some((m.name().to_string(), top.into()))
                })
                .collect();
            let h = history.get(id).cloned().unwrap_or_default();
            let mut v = json!({ "digest": d, "history": h, "tops": tops });
            if with_types {
                v["types"] = json!(types);
            }
            if let Some(cl) = r.placement(d) {
                v["cl"] = json!(cl);
            }
            (*id, v)
        })
        .collect())
}

/// `review`'s text: each cluster (all, or those in `clusters`) in turn.
pub(crate) fn review_text(loaded: &Loaded, clusters: &[ClusterId]) -> Result<String> {
    let r = Reviewed::read(loaded, false, false)?;
    let caveats = caveats(loaded)
        .into_iter()
        .map(|c| format!("caveat: {c}\n"));
    let clusters = r
        .digests
        .iter()
        .filter(|(id, _)| clusters.is_empty() || clusters.contains(id))
        .map(|(id, d)| {
            let history = r.history.get(id).map_or(&[][..], Vec::as_slice);
            render(*id, d, r.placement(d).as_ref(), history)
        });
    Ok(caveats.chain(clusters).collect())
}

/// What `loaded`'s readers should know of its statistics
/// (`annotate.stats.caveats`).
pub(crate) fn caveats(loaded: &Loaded) -> Vec<String> {
    let stats = loaded.manifest.annotate.stats.as_ref();
    stats
        .and_then(|s| s.get(CAVEATS))
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(|c| c.as_str().map(str::to_string))
        .collect()
}

pub fn run_review(args: &ReviewArgs) -> Result<()> {
    let loaded = run::load(&args.from)?;
    if args.json {
        let out = review_json(&loaded, &args.cluster, true)?;
        println!("{}", serde_json::to_string_pretty(&out)?);
        return Ok(());
    }
    print!("{}", review_text(&loaded, &args.cluster)?);
    Ok(())
}

fn render(
    id: ClusterId,
    d: &Digest,
    placed: Option<&TreePlacement>,
    history: &[rounds::HistoryEntry],
) -> String {
    use std::fmt::Write;
    let num = |v: Option<f32>| v.map_or_else(|| "-".into(), |v| format!("{v:.3}"));
    let mut s = String::new();
    let label = d.label.as_deref().unwrap_or("unassigned");
    let name = crate::figures::cluster_name(id);
    let _ = writeln!(s, "{name}  n={}  label={label}", d.size);
    for c in &d.calls {
        let _ = writeln!(s, "  call  {:<24} q={}", c.label, num(c.q));
    }
    for t in &d.terms {
        let _ = writeln!(
            s,
            "  {:<4}  {:<40} effect={:.3} q={}",
            t.source,
            t.term,
            t.effect,
            num(t.q)
        );
    }
    // The ontology walk's placement, else the type tree's.
    match (&d.cl, placed) {
        (Some(cl), _) => {
            let abst = if cl.abstained { " (abstained)" } else { "" };
            let _ = writeln!(s, "  CL    {} {}{abst}", cl.id, cl.name);
        }
        (None, Some(p)) => {
            let of = match p.via {
                Via::Group if p.name != label => format!(" (the group of {label})"),
                _ => String::new(),
            };
            let _ = writeln!(s, "  CL    {} {}{of}", p.id, p.name);
        }
        (None, None) => {}
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
pub(crate) mod tests;
