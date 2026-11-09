//! A round's marker genes × clusters as a heatmap and as numbers (`lupin
//! figure heatmap`, the MCP's `annotate_heatmap`): for each panel type its
//! top markers, by their largest log2 fold change in any cluster, and their
//! value in every cluster: log2 fold change over the other clusters, mean
//! level, that level z-scored across the clusters, or the share of the
//! cluster's cells the gene is detected in. The level is the round's
//! observed profile, or the topic model's expected level. Markers match the
//! profile's genes as every gene view of a round matches them
//! ([`crate::annotate::gene_rows::LevelRows`]). The figure is one SVG,
//! written as `.svg`, `.pdf` or `.png`, and the numbers go beside it as
//! `.json`.

use super::Figure;
use crate::annotate::expression::{label_marker_set, round_panel, Expression};
use crate::annotate::gene_rows::CHANNELS;
use crate::annotate::rounds::ClusterId;
use crate::manifest::no_counts::{self, Decoder};
use crate::manifest::rounds::{cluster_labels, read_cells, read_tree, Cells};
use crate::manifest::run::{self, resolve, Loaded};
use anyhow::{Context, Result};
use data_beans::aux::feature_rows::SPLICED;
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::traits::{IoOps, MatOps, MatWithNames};
use legume_plot::palette::sample_blue_red;
use legume_plot::svg_emit::escape_xml;
use legume_tui::chart::heatmap::zscore_rows;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fmt::Write as _;
use std::ops::AddAssign;
use std::path::{Path, PathBuf};

/// What a heatmap cell shows of a gene in a cluster.
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum HeatValue {
    /// log2 fold change over the other clusters' mean.
    Log2fc,
    /// The cluster's mean level.
    Mean,
    /// The mean level z-scored across the clusters.
    Zscore,
    /// The share of the cluster's cells with a count above zero.
    Frac,
}

impl HeatValue {
    fn title(self) -> &'static str {
        match self {
            Self::Log2fc => "log2 fold change over the other clusters",
            Self::Mean => "mean level",
            Self::Zscore => "mean level, z across clusters",
            Self::Frac => "share of cells detected",
        }
    }

    /// Whether its colours run either side of zero.
    fn diverging(self) -> bool {
        matches!(self, Self::Log2fc | Self::Zscore)
    }
}

/// Which level the markers are measured on.
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum HeatSource {
    /// The round's cluster profile, from the observed counts.
    Observed,
    /// The topic model's expected level: each cell's topic shares times
    /// the topics' gene distributions, per 10,000, averaged over the
    /// cluster.
    Expected,
}

impl HeatSource {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Observed => "observed",
            Self::Expected => "expected",
        }
    }
}

/// The order of the clusters, left to right.
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ClusterOrder {
    /// By label, in the order of the types' rows, so a type's markers and
    /// its clusters meet on the diagonal.
    Label,
    /// By cluster id, as lupin numbers them.
    Id,
}

#[derive(clap::Args, Debug)]
pub struct HeatmapArgs {
    #[arg(short, long, help = "A round (its manifest or prefix)")]
    pub from: Box<str>,
    #[arg(
        short,
        long,
        help = "The figure, ending .pdf, .svg or .png; its numbers go beside it as .json"
    )]
    pub out: Box<str>,
    #[arg(
        long = "type",
        help = "A cell type whose markers to show (repeatable); every panel type when none"
    )]
    pub types: Vec<Box<str>>,
    #[arg(long, default_value_t = 5, help = "Markers shown per type")]
    pub top: usize,
    #[arg(long, value_enum, default_value_t = HeatValue::Zscore)]
    pub value: HeatValue,
    #[arg(long, value_enum, default_value_t = HeatSource::Observed)]
    pub source: HeatSource,
    #[arg(long, value_enum, default_value_t = ClusterOrder::Label)]
    pub order: ClusterOrder,
    #[arg(
        long,
        default_value = SPLICED,
        value_parser = clap::builder::PossibleValuesParser::new(CHANNELS),
        help = "Where a gene has several rows (a faba file's), the count channel its level is read from"
    )]
    pub channel: String,
}

/// A heatmap's numbers.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct HeatData {
    pub round: String,
    pub value: HeatValue,
    pub source: HeatSource,
    /// Left to right.
    pub clusters: Vec<HeatCluster>,
    /// Top to bottom, grouped by type.
    pub rows: Vec<HeatRow>,
    /// Types none of whose markers is measured, so they have no rows.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub unmeasured: Vec<String>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct HeatCluster {
    pub id: ClusterId,
    pub label: Option<String>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct HeatRow {
    /// The marker as the panel names it.
    pub gene: String,
    /// The profile's row it was read from.
    pub row: String,
    #[serde(rename = "type")]
    pub type_: String,
    /// Per cluster, in [`HeatData::clusters`] order; `None` where it has
    /// no value.
    pub values: Vec<Option<f32>>,
}

/// `lupin figure heatmap`.
pub fn run_heatmap(args: &HeatmapArgs) -> Result<()> {
    let (base, formats) = figure_base(&args.out)?;
    let loaded = run::load(&args.from)?;
    let data = heat_data(&loaded, args)?;
    let fig = svg(&data);
    legume_numeric::matrix::common_io::mkdir_parent(&args.out)?;
    legume_plot::write_figure(&fig.svg, fig.w, fig.h, &base, formats)?;
    let json = numbers_path(&args.out);
    std::fs::write(&json, serde_json::to_string_pretty(&data)?)
        .with_context(|| format!("writing {}", json.display()))?;
    log::info!("wrote {}", json.display());
    Ok(())
}

/// The figure `out` without its extension, and the one format it names.
pub fn figure_base(out: &str) -> Result<(String, legume_plot::FigureFormats)> {
    let path = Path::new(out);
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    let formats = legume_plot::FigureFormats {
        svg: ext == "svg",
        png: ext == "png",
        pdf: ext == "pdf",
    };
    anyhow::ensure!(
        formats.svg || formats.png || formats.pdf,
        "{out}: a figure ends .pdf, .svg or .png"
    );
    Ok((
        path.with_extension("").to_string_lossy().into_owned(),
        formats,
    ))
}

/// The numbers file beside the figure `out`.
pub fn numbers_path(out: &str) -> PathBuf {
    Path::new(out).with_extension("json")
}

/// The heatmap's numbers for round `loaded`.
pub fn heat_data(loaded: &Loaded, args: &HeatmapArgs) -> Result<HeatData> {
    let cells = read_cells(loaded)?;
    let panel = round_panel(loaded)?;
    anyhow::ensure!(!panel.is_empty(), "the round records no marker panel");
    let tree = read_tree(loaded).ok().flatten();
    let table = match args.source {
        HeatSource::Observed => crate::manifest::rounds::round_profile(loaded, &cells)?
            .context("the round records no expression profile (annotate.cluster_expression)")?,
        HeatSource::Expected => expected_profile(loaded, &cells)?,
    };
    let expr = Expression::with_channel(table, &args.channel);
    let labels = cluster_labels(loaded)?;
    let types: Vec<String> = if args.types.is_empty() {
        panel.keys().cloned().collect()
    } else {
        args.types.iter().map(ToString::to_string).collect()
    };
    let all_cols: Vec<usize> = labels.keys().filter_map(|&id| expr.col(id)).collect();
    // Each type's measured markers, by their largest fold change.
    let mut chosen: Vec<(String, String, usize)> = Vec::new();
    let mut unmeasured = Vec::new();
    for t in &types {
        let set = label_marker_set(tree.as_ref(), &panel, t);
        anyhow::ensure!(
            !set.is_empty(),
            "{t:?} is not a type of the round's panel, nor a group of its types, with markers"
        );
        let mut genes: Vec<(String, usize, f32)> = set
            .into_iter()
            .filter_map(|g| {
                let r = expr.row(&g)?;
                let best = (all_cols.iter())
                    .filter_map(|&k| expr.log2fc(r, &[k]))
                    .fold(f32::NEG_INFINITY, f32::max);
                Some((g, r, best))
            })
            .collect();
        if genes.is_empty() {
            unmeasured.push(t.clone());
            continue;
        }
        genes.sort_by(|a, b| b.2.total_cmp(&a.2).then_with(|| a.0.cmp(&b.0)));
        chosen.extend(
            genes
                .into_iter()
                .take(args.top)
                .map(|(g, r, _)| (t.clone(), g, r)),
        );
    }
    let clusters = ordered(&labels, &types, tree.as_ref(), &panel, args.order);
    let cols: Vec<Option<usize>> = clusters.iter().map(|c| expr.col(c.id)).collect();
    let level = |r: usize| -> Vec<Option<f32>> {
        (cols.iter())
            .map(|k| k.map(|k| expr.table.mat[(r, k)]))
            .collect()
    };
    let frac = match args.value {
        HeatValue::Frac => {
            let cluster_of: HashMap<&str, ClusterId> = (cells.names.iter())
                .zip(&cells.clusters)
                .filter_map(|(n, c)| Some((&**n, (*c)?)))
                .collect();
            let ids: Vec<ClusterId> = clusters.iter().map(|c| c.id).collect();
            let rows: Vec<&str> = (chosen.iter())
                .map(|&(_, _, r)| &*expr.table.rows[r])
                .collect();
            crate::manifest::annotate::detected_share(loaded, &cluster_of, &ids, &rows)?
        }
        _ => Vec::new(),
    };
    let values = (chosen.iter().enumerate()).map(|(i, &(_, _, r))| match args.value {
        HeatValue::Mean => level(r),
        HeatValue::Zscore => zscored(&level(r)),
        HeatValue::Log2fc => (cols.iter())
            .map(|k| k.and_then(|k| expr.log2fc(r, &[k])))
            .collect(),
        HeatValue::Frac => frac[i].clone(),
    });
    let rows = (chosen.iter().cloned())
        .zip(values)
        .map(|((type_, gene, r), values)| HeatRow {
            gene,
            row: expr.table.rows[r].to_string(),
            type_,
            values,
        })
        .collect();
    Ok(HeatData {
        round: loaded.file.to_string_lossy().into_owned(),
        value: args.value,
        source: args.source,
        clusters,
        rows,
        unmeasured,
    })
}

/// `v` z-scored over its values, as the screen's heatmaps z-score a row.
fn zscored(v: &[Option<f32>]) -> Vec<Option<f32>> {
    let mut z: Vec<f32> = v.iter().map(|x| x.unwrap_or(f32::NAN)).collect();
    let cols = z.len();
    zscore_rows(&mut z, cols, f32::INFINITY);
    z.into_iter().map(|x| x.is_finite().then_some(x)).collect()
}

/// The clusters in `order`: by label, the place of a cluster's label (or,
/// for a coarse label, of its first member) among `types`, then the rest by
/// label; else by id.
fn ordered(
    labels: &BTreeMap<ClusterId, Option<String>>,
    types: &[String],
    tree: Option<&crate::annotate::celltype_tree::TypeTree>,
    panel: &crate::annotate::expression::Panel,
    order: ClusterOrder,
) -> Vec<HeatCluster> {
    let key = crate::annotate::markers::label_key;
    let place = |label: &Option<String>| -> usize {
        let Some(l) = label else { return usize::MAX };
        let at = |t: &str| types.iter().position(|u| key(u) == key(t));
        at(l)
            .or_else(|| {
                let members = tree.and_then(|t| t.coarse_members(l, panel))?;
                members.iter().filter_map(|m| at(m)).min()
            })
            .unwrap_or(usize::MAX - 1)
    };
    let mut out: Vec<HeatCluster> = labels
        .iter()
        .map(|(&id, label)| HeatCluster {
            id,
            label: label.clone(),
        })
        .collect();
    if order == ClusterOrder::Label {
        out.sort_by(|a, b| {
            (place(&a.label).cmp(&place(&b.label)))
                .then_with(|| a.label.cmp(&b.label))
                .then_with(|| a.id.cmp(&b.id))
        });
    }
    out
}

/// The topic model's expected level per cluster (columns `K{id}`), as a
/// pass without counts decodes it ([`no_counts`]): the topics' gene
/// distributions times the cluster's mean topic shares, per
/// [`no_counts::NOMINAL_DEPTH`] counts.
fn expected_profile(loaded: &Loaded, cells: &Cells) -> Result<MatWithNames<Mat>> {
    let Some((Decoder::Mixture, latent, dictionary)) = no_counts::decoder_of(loaded) else {
        anyhow::bail!(
            "a {} run records no topic model (a latent and a dictionary) whose expected \
             level `source: expected` reads: use `source: observed`",
            loaded.manifest.kind
        );
    };
    let read = |rel: &str| {
        let path = resolve(&loaded.dir, rel);
        Mat::from_parquet_with_row_names(&path, Some(0)).with_context(|| format!("reading {path}"))
    };
    let (theta, beta) = (read(latent)?, read(dictionary)?);
    let k = theta.mat.ncols();
    anyhow::ensure!(
        beta.mat.ncols() == k,
        "the latent has {k} topics and the dictionary {}",
        beta.mat.ncols()
    );
    let row: HashMap<&str, usize> = (theta.rows.iter().enumerate())
        .map(|(i, n)| (&**n, i))
        .collect();
    let ids: Vec<ClusterId> = (cells.clusters.iter().flatten().copied())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    // Each cluster's mean topic shares, topics × clusters.
    let theta_kn = theta.mat.transpose().normalize_exp_logits_columns();
    let mut shares = Mat::zeros(k, ids.len());
    let mut n = vec![0usize; ids.len()];
    for (name, c) in cells.names.iter().zip(&cells.clusters) {
        let (Some(c), Some(&i)) = (c, row.get(&**name)) else {
            continue;
        };
        let j = ids.binary_search(c).expect("its cluster is listed");
        shares.column_mut(j).add_assign(&theta_kn.column(i));
        n[j] += 1;
    }
    for (j, &n) in n.iter().enumerate() {
        shares.column_mut(j).unscale_mut(n.max(1) as f32);
    }
    let depth = no_counts::NOMINAL_DEPTH as f32;
    let mat = no_counts::gene_distributions(&beta.mat) * shares * depth;
    Ok(MatWithNames {
        rows: beta.rows,
        cols: ids.iter().map(|id| format!("K{id}").into()).collect(),
        mat,
    })
}

/// Pixels of a heatmap cell, and of its text.
const CELL: u32 = 16;
const FONT: f32 = 11.0;
/// Width of a character of the text, roughly.
const CHAR_W: f32 = 6.4;
/// A z-score's colours stop here.
const Z_CLIP: f32 = 3.0;
const MISSING: &str = "#e6e6e6";

/// The heatmap of `data` as an SVG: rows `GENE · type` grouped by type,
/// clusters `K{id} · label` across the top, a colour key on the right.
pub(crate) fn svg(data: &HeatData) -> Figure {
    let (nr, nc) = (data.rows.len(), data.clusters.len());
    let row_labels: Vec<String> = (data.rows.iter())
        .map(|r| format!("{} · {}", r.gene, r.type_))
        .collect();
    let col_labels: Vec<String> = (data.clusters.iter())
        .map(|c| match &c.label {
            Some(l) => format!("{} · {l}", super::cluster_name(c.id)),
            None => super::cluster_name(c.id),
        })
        .collect();
    let longest = |v: &[String]| v.iter().map(|s| s.chars().count()).max().unwrap_or(0) as f32;
    let left = (longest(&row_labels) * CHAR_W + 12.0).ceil() as u32;
    // Column labels at 45°.
    let top = (longest(&col_labels) * CHAR_W * 0.71 + 40.0).ceil() as u32;
    // Right of the cells: the colour key, and as far as the last columns'
    // labels reach.
    let reach = (col_labels.iter().enumerate())
        .map(|(j, l)| {
            let past = (nc - j) as f32 * CELL as f32 - CELL as f32 / 2.0;
            l.chars().count() as f32 * CHAR_W * 0.71 - past
        })
        .fold(0.0f32, f32::max);
    let key_w = 110.max((reach + 12.0).ceil() as u32);
    let w = left + nc as u32 * CELL + key_w;
    let h = (top + nr.max(1) as u32 * CELL + 20).max(top + 150);
    let lim = colour_limit(data);
    let colour = |v: Option<f32>| -> String {
        let Some(v) = v.filter(|v| v.is_finite()) else {
            return MISSING.into();
        };
        let low = if data.value.diverging() { -1.0 } else { 0.0 };
        let (r, g, b) = sample_blue_red(0.5 + 0.5 * (v / lim).clamp(low, 1.0));
        format!("#{}", crate::hex(&[r, g, b]))
    };
    let mut s = String::new();
    let _ = write!(
        s,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" font-family="Helvetica, Arial, sans-serif" font-size="{FONT}"><rect width="100%" height="100%" fill="white"/>"#
    );
    let title = format!(
        "markers × clusters · {} · {}",
        data.value.title(),
        data.source.as_str()
    );
    let _ = write!(
        s,
        r#"<text x="8" y="16" font-weight="bold">{}</text>"#,
        escape_xml(&title)
    );
    for (j, label) in col_labels.iter().enumerate() {
        let x = left as f32 + (j as f32 + 0.5) * CELL as f32;
        let y = top as f32 - 6.0;
        let _ = write!(
            s,
            r#"<text x="{x:.1}" y="{y:.1}" transform="rotate(-45 {x:.1} {y:.1})">{}</text>"#,
            escape_xml(label)
        );
    }
    for (i, (row, label)) in data.rows.iter().zip(&row_labels).enumerate() {
        let y = top + i as u32 * CELL;
        let _ = write!(
            s,
            r#"<text x="{}" y="{:.1}" text-anchor="end">{}</text>"#,
            left - 6,
            y as f32 + CELL as f32 * 0.72,
            escape_xml(label)
        );
        for (j, v) in row.values.iter().enumerate() {
            let _ = write!(
                s,
                r#"<rect x="{}" y="{y}" width="{CELL}" height="{CELL}" fill="{}"/>"#,
                left + j as u32 * CELL,
                colour(*v)
            );
        }
        // A line between one type's rows and the next's.
        if i > 0 && data.rows[i - 1].type_ != row.type_ {
            let _ = write!(
                s,
                r#"<line x1="{}" y1="{y}" x2="{}" y2="{y}" stroke="white" stroke-width="2"/>"#,
                left,
                left + nc as u32 * CELL
            );
        }
    }
    // The colour key: low to high, bottom to top.
    let (kx, ky, kh) = (left + nc as u32 * CELL + 20, top, 120u32);
    let low = if data.value.diverging() { -lim } else { 0.0 };
    for step in 0..kh {
        let v = low + (lim - low) * (1.0 - step as f32 / kh as f32);
        let _ = write!(
            s,
            r#"<rect x="{kx}" y="{}" width="14" height="1.5" fill="{}"/>"#,
            ky + step,
            colour(Some(v))
        );
    }
    for (v, y) in [(lim, ky + 8), (low, ky + kh)] {
        let _ = write!(s, r#"<text x="{}" y="{y}">{v:.2}</text>"#, kx + 20);
    }
    s.push_str("</svg>");
    Figure { svg: s, w, h }
}

/// The value at which the colours stop: a share's 1, a z-score's
/// [`Z_CLIP`], else the largest magnitude shown.
fn colour_limit(data: &HeatData) -> f32 {
    let largest = (data.rows.iter())
        .flat_map(|r| r.values.iter().flatten())
        .map(|v| v.abs())
        .filter(|v| v.is_finite())
        .fold(0.0f32, f32::max);
    match data.value {
        HeatValue::Frac => 1.0,
        HeatValue::Zscore => largest.min(Z_CLIP),
        _ => largest,
    }
    .max(f32::EPSILON)
}

/// The heatmap of `data` as PNG bytes.
pub(crate) fn png(data: &HeatData) -> Result<Vec<u8>> {
    let img = super::render_svg(&svg(data))?;
    let mut bytes = Vec::new();
    img.write_to(
        &mut std::io::Cursor::new(&mut bytes),
        image::ImageFormat::Png,
    )?;
    Ok(bytes)
}

#[cfg(test)]
#[path = "heatmap_tests.rs"]
mod tests;
