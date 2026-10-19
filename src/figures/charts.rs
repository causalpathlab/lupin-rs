//! A round's charts, drawn in place of its cell map as senna draws its own
//! through `legume_tui::chart`: a structure plot of each cell's mixture
//! (its posterior over the cell types, its membership of a units round's
//! regions, or a topic run's topics) stacked from the bottom in a panel per
//! cluster, and heatmaps of features × clusters (the panel's markers, the
//! GO terms a units pass aligned, the cell types), z-scored by row. Each
//! reads the table its manifest records on first use; without it the chart
//! says which instead of drawing.

use super::{cell_index, cluster_name, distinct_colours, Colouring, Panel, Size, THEME};
use crate::manifest::rounds::read_table;
use crate::manifest::run::{Loaded, RunKind};
use anyhow::{Context, Result};
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::parquet::read_table_columns;
use legume_numeric::matrix::traits::IoOps;
use legume_tui::chart::heatmap::{draw_heat, zscore_rows, Heat, HeatStyle, CLIP};
use legume_tui::chart::note;
use legume_tui::chart::structure::{figure, BarStyle, Figure, Grouping, Layout, Names, Weights};
use legume_tui::export::{write_pdf, Page};
use legume_tui::map::colour::{to_linear, Lin};
use legume_tui::map::raster::Frame;
use legume_tui::map::text::Font;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;
use std::sync::{Arc, OnceLock};

use super::units::{STARS, STAR_Q};

/// What a structure plot stacks for each cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Mixture {
    /// The cell types: the pass's per-cell posterior, else a units round's
    /// region membership.
    Types,
    /// A topic run's topic shares.
    Topics,
}

impl Mixture {
    fn slug(self) -> &'static str {
        match self {
            Self::Types => "types",
            Self::Topics => "topics",
        }
    }

    pub(crate) fn title(self) -> &'static str {
        match self {
            Self::Types => "cell types per cell, by cluster",
            Self::Topics => "topics per cell, by cluster",
        }
    }
}

/// Which features × clusters a heatmap shows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum HeatKind {
    Markers,
    Terms,
    Types,
}

impl HeatKind {
    /// In the order `c` steps through them.
    pub(crate) const ALL: [Self; 3] = [Self::Markers, Self::Terms, Self::Types];

    fn slug(self) -> &'static str {
        match self {
            Self::Markers => "markers",
            Self::Terms => "terms",
            Self::Types => "types",
        }
    }

    pub(crate) fn title(self) -> &'static str {
        match self {
            Self::Markers => "markers × clusters",
            Self::Terms => "GO terms × clusters",
            Self::Types => "cell types × clusters",
        }
    }
}

impl Panel {
    /// Whether it is a chart, drawn in place of the map.
    pub(crate) fn is_chart(self) -> bool {
        matches!(self, Self::Structure { .. } | Self::Heatmap { .. })
    }

    pub(super) fn chart_slug(self) -> Option<String> {
        match self {
            Self::Structure { of } => Some(format!("structure_{}", of.slug())),
            Self::Heatmap { of } => Some(format!("heatmap_{}", of.slug())),
            _ => None,
        }
    }
}

/// Each cell's mixture: the components' names and `cells × k` weights,
/// row-major, NaN for a cell the table lacks.
#[derive(Debug)]
pub(crate) struct Mix {
    pub(crate) names: Vec<String>,
    pub(crate) w: Vec<f32>,
}

/// A heatmap's features × clusters, z-scored by row; NaN where missing.
#[derive(Debug)]
pub(crate) struct HeatTable {
    pub(crate) title: String,
    pub(crate) rows: Vec<String>,
    pub(crate) cols: Vec<String>,
    pub(crate) z: Vec<f32>,
}

/// A table the round records, found; or why there is none.
type Found = Result<String, String>;

/// Why a trajectory's figures have no charts.
const NOT_A_ROUND: &str = "charts are a round's";
/// Why a round has no cell types per cell.
const NO_TYPES: &str = "the round records no per-cell posterior (annotate.annotation) \
                        or units membership (annotate.celltype_membership)";
/// Why a round has no cell types per cluster.
const NO_TYPE_Q: &str = "the round records no enrichment Q (annotate.cluster_celltype_q) \
                         or units membership (annotate.celltype_membership)";

/// The tables a round's charts read, and each chart's data once read.
pub(crate) struct Charts {
    /// The round's cells by name, which the per-cell tables are read onto.
    cells: Vec<Box<str>>,
    posterior: Found,
    membership: Found,
    topics: Found,
    expression: Found,
    markers: Found,
    alignment: Found,
    type_q: Found,
    mixtures: [OnceLock<Result<Arc<Mix>, String>>; 2],
    heats: [OnceLock<Result<Arc<HeatTable>, String>>; 3],
}

impl Default for Charts {
    fn default() -> Self {
        let none = || Err(NOT_A_ROUND.to_string());
        Self {
            cells: Vec::new(),
            posterior: none(),
            membership: none(),
            topics: none(),
            expression: none(),
            markers: none(),
            alignment: none(),
            type_q: none(),
            mixtures: Default::default(),
            heats: Default::default(),
        }
    }
}

impl Charts {
    /// The tables `loaded` records for the charts of its `cells`.
    pub(crate) fn of(loaded: &Loaded, cells: &[Box<str>]) -> Self {
        let a = &loaded.manifest.annotate;
        let found = |rel: &Option<String>, field: &str| -> Found {
            let rel = rel
                .as_deref()
                .ok_or_else(|| format!("the round records no {field}"))?;
            loaded
                .recorded(rel)
                .ok_or_else(|| format!("{rel} ({field}) is recorded but gone"))
        };
        let m = &loaded.manifest;
        let topic_run = matches!(
            m.kind,
            RunKind::Topic | RunKind::Itopic | RunKind::JointTopic
        );
        Self {
            cells: cells.to_vec(),
            posterior: found(&a.annotation, "annotate.annotation"),
            membership: found(&a.celltype_membership, "annotate.celltype_membership"),
            topics: if topic_run {
                found(&m.outputs.latent, "outputs.latent")
            } else {
                Err(format!("a {} run has no topics", m.kind))
            },
            expression: found(&a.cluster_expression, "annotate.cluster_expression"),
            markers: found(&a.markers, "annotate.markers"),
            alignment: found(&a.term_alignment, "annotate.term_alignment"),
            type_q: found(&a.cluster_celltype_q, "annotate.cluster_celltype_q"),
            ..Self::default()
        }
    }
}

/// The first of `options` found; else why not, a table recorded but gone
/// before one never recorded, which `none` says.
fn first_found<'a>(options: &[&'a Found], none: &str) -> Result<&'a str, String> {
    if let Some(Ok(p)) = options.iter().find(|f| f.is_ok()) {
        return Ok(p);
    }
    match options
        .iter()
        .find_map(|f| f.as_ref().err().filter(|e| e.contains("gone")))
    {
        Some(e) => Err(e.clone()),
        None => Err(none.into()),
    }
}

impl super::CellFigures {
    /// Whether the round records the table `of` stacks.
    pub(crate) fn has_mixture(&self, of: Mixture) -> bool {
        let c = &self.charts;
        match of {
            Mixture::Types => c.posterior.is_ok() || c.membership.is_ok(),
            Mixture::Topics => c.topics.is_ok(),
        }
    }

    /// Each cell's mixture `of`, read once; why not.
    pub(crate) fn mixture(&self, of: Mixture) -> Result<Arc<Mix>, String> {
        let c = &self.charts;
        let slot = match of {
            Mixture::Types => 0,
            Mixture::Topics => 1,
        };
        c.mixtures[slot]
            .get_or_init(|| {
                let path = match of {
                    Mixture::Types => first_found(&[&c.posterior, &c.membership], NO_TYPES)?,
                    Mixture::Topics => c.topics.as_deref().map_err(Clone::clone)?,
                };
                let mut mix = read_cells(path, &c.cells).map_err(|e| format!("{e:#}"))?;
                if of == Mixture::Topics {
                    shares_of_logs(&mut mix.w, mix.names.len());
                }
                Ok(Arc::new(mix))
            })
            .clone()
    }

    /// `mix`'s structure plot: a panel per cluster in id order, cells no
    /// cluster holds last.
    pub(crate) fn structure(&self, mix: &Mix) -> Layout {
        let (ids, of) = self.groups();
        let names: Vec<String> = ids.iter().map(|&k| cluster_name(k)).collect();
        Layout::build(
            &Weights {
                k: mix.names.len(),
                w: &mix.w,
                valid: None,
            },
            &Grouping {
                of: &of,
                names: &names,
                order: None,
            },
        )
    }

    /// The round's clusters in id order, and each cell's place among them
    /// (`u32::MAX` for none).
    fn groups(&self) -> (Vec<i32>, Vec<u32>) {
        let ids: Vec<i32> = (self.cluster.iter())
            .filter(|&&k| k >= 0)
            .copied()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let of = (self.cluster.iter())
            .map(|k| ids.binary_search(k).map_or(u32::MAX, |i| i as u32))
            .collect();
        (ids, of)
    }

    /// The heatmap `of`, read once; why not.
    pub(crate) fn heat(&self, of: HeatKind) -> Result<Arc<HeatTable>, String> {
        let slot = HeatKind::ALL.iter().position(|&h| h == of).unwrap_or(0);
        self.charts.heats[slot]
            .get_or_init(|| {
                let made = match of {
                    HeatKind::Markers => self.markers_heat(),
                    HeatKind::Terms => self.terms_heat(),
                    HeatKind::Types => self.types_heat(),
                };
                made.map(Arc::new)
            })
            .clone()
    }

    /// The panel's markers (rows `GENE · type`, grouped by type) over the
    /// clusters of `cluster_expression`, as that pass made them.
    fn markers_heat(&self) -> Result<HeatTable, String> {
        let c = &self.charts;
        let path = c.expression.as_deref().map_err(Clone::clone)?;
        let markers = c.markers.as_deref().map_err(Clone::clone)?;
        let t = Mat::from_parquet_with_row_names(path, Some(0))
            .with_context(|| format!("reading {path}"))
            .map_err(|e| format!("{e:#}"))?;
        let panel = read_panel(markers).map_err(|e| format!("{e:#}"))?;
        // A marker matches its row as written, else by case.
        let exact: HashMap<&str, usize> = (t.rows.iter().enumerate())
            .map(|(r, g)| (&**g, r))
            .collect();
        let folded: HashMap<String, usize> = (t.rows.iter().enumerate())
            .map(|(r, g)| (g.to_ascii_uppercase(), r))
            .collect();
        let cols = t.cols.len();
        let (mut rows, mut z) = (Vec::new(), Vec::new());
        for (ty, genes) in &panel {
            for g in genes {
                let found = exact.get(g.as_str()).copied();
                let Some(r) = found.or_else(|| folded.get(&g.to_ascii_uppercase()).copied()) else {
                    continue;
                };
                rows.push(format!("{g} · {ty}"));
                z.extend((0..cols).map(|k| t.mat[(r, k)]));
            }
        }
        if rows.is_empty() {
            return Err(format!("none of the panel's markers is in {path}"));
        }
        zscore_rows(&mut z, cols, CLIP);
        Ok(HeatTable {
            title: "markers × clusters · mean expression, z by row".into(),
            rows,
            cols: t.cols.iter().map(ToString::to_string).collect(),
            z,
        })
    }

    /// The GO terms some cluster stars (as the map stars them: q below
    /// [`STAR_Q`], the best [`STARS`]), their scores z-scored over every
    /// cluster, missing where q is not below [`STAR_Q`].
    fn terms_heat(&self) -> Result<HeatTable, String> {
        let path = self.charts.alignment.as_deref().map_err(Clone::clone)?;
        let (s, v) = read_table_columns(path, &["cluster", "term"], &["score", "q"])
            .with_context(|| format!("reading {path}"))
            .map_err(|e| format!("{e:#}"))?;
        let mut clusters: Vec<u32> = (s[0].iter())
            .filter_map(|k| crate::annotate::rounds::parse_cluster_id(k))
            .collect();
        clusters.sort_unstable();
        clusters.dedup();
        let mut terms: Vec<&str> = Vec::new();
        for t in &s[1] {
            if !terms.contains(&&**t) {
                terms.push(t);
            }
        }
        let (n, m) = (terms.len(), clusters.len());
        let (mut score, mut q) = (vec![f32::NAN; n * m], vec![f32::NAN; n * m]);
        for i in 0..s[0].len() {
            let Some(k) = crate::annotate::rounds::parse_cluster_id(&s[0][i]) else {
                continue;
            };
            let (Ok(k), Some(t)) = (
                clusters.binary_search(&k),
                terms.iter().position(|t| **t == *s[1][i]),
            ) else {
                continue;
            };
            score[t * m + k] = v[0][i] as f32;
            q[t * m + k] = v[1][i] as f32;
        }
        // Each cluster's stars, cluster by cluster.
        let mut shown: Vec<usize> = Vec::new();
        for k in 0..m {
            let mut best: Vec<usize> = (0..n)
                .filter(|&t| q[t * m + k] < STAR_Q && score[t * m + k].is_finite())
                .collect();
            best.sort_by(|&a, &b| score[b * m + k].total_cmp(&score[a * m + k]));
            for t in best.into_iter().take(STARS) {
                if !shown.contains(&t) {
                    shown.push(t);
                }
            }
        }
        if shown.is_empty() {
            return Err(format!(
                "no cluster has a GO term with q below {STAR_Q} in {path}"
            ));
        }
        let mut z: Vec<f32> = (shown.iter())
            .flat_map(|&t| score[t * m..(t + 1) * m].to_vec())
            .collect();
        zscore_rows(&mut z, m, CLIP);
        for (r, &t) in shown.iter().enumerate() {
            for k in 0..m {
                let q = q[t * m + k];
                if q.is_nan() || q >= STAR_Q {
                    z[r * m + k] = f32::NAN;
                }
            }
        }
        // Named as the map lists them, else by id.
        let names: HashMap<&str, &str> = (self.units.iter())
            .flat_map(|u| u.by_cluster.values().flatten())
            .map(|t| (&*t.id, t.name.as_str()))
            .collect();
        Ok(HeatTable {
            title: format!("GO terms × clusters · score, z by row; q ≥ {STAR_Q} blank"),
            rows: (shown.iter())
                .map(|&t| names.get(terms[t]).unwrap_or(&terms[t]).to_string())
                .collect(),
            cols: clusters.iter().map(|&k| cluster_name(k)).collect(),
            z,
        })
    }

    /// The cell types over the clusters: the enrichment's Q, else each
    /// cluster's mean membership of a units round's regions.
    fn types_heat(&self) -> Result<HeatTable, String> {
        let c = &self.charts;
        let path = first_found(&[&c.type_q, &c.membership], NO_TYPE_Q)?;
        let fail = |e: anyhow::Error| format!("{e:#}");
        let (title, rows, cols, mut z) = if c.type_q.is_ok() {
            let t = read_table(path).map_err(fail)?;
            let (n, m) = (t.cols.len(), t.rows.len());
            let z: Vec<f32> = (0..n)
                .flat_map(|j| (0..m).map(move |k| (j, k)))
                .map(|(j, k)| t.values[k * n + j])
                .collect();
            let cols: Vec<String> = t.rows.iter().map(|&k| cluster_name(k)).collect();
            ("enrichment Q", t.cols, cols, z)
        } else {
            let mix = read_cells(path, &c.cells).map_err(fail)?;
            let (ids, of) = self.groups();
            let n = mix.names.len();
            let (mut sum, mut count) = (vec![0f64; n * ids.len()], vec![0usize; ids.len()]);
            for (i, &g) in of.iter().enumerate() {
                let row = &mix.w[i * n..(i + 1) * n];
                if g == u32::MAX || !row.iter().all(|v| v.is_finite()) {
                    continue;
                }
                let g = g as usize;
                count[g] += 1;
                for (j, &v) in row.iter().enumerate() {
                    sum[j * ids.len() + g] += f64::from(v);
                }
            }
            let z: Vec<f32> = (sum.iter().enumerate())
                .map(|(i, &s)| match count[i % ids.len()] {
                    0 => f32::NAN,
                    c => (s / c as f64) as f32,
                })
                .collect();
            let cols: Vec<String> = ids.iter().map(|&k| cluster_name(k)).collect();
            ("units membership share", mix.names, cols, z)
        };
        zscore_rows(&mut z, cols.len(), CLIP);
        Ok(HeatTable {
            title: format!("cell types × clusters · {title}, z by row"),
            rows,
            cols,
            z,
        })
    }

    /// Why `panel`, a chart, has no table to draw from; `None` when the
    /// round records its tables, which are read only when it is drawn.
    pub(crate) fn chart_missing(&self, panel: Panel) -> Option<String> {
        let c = &self.charts;
        let found = |f: &Found| f.as_ref().map(drop).map_err(Clone::clone);
        let why = match panel {
            Panel::Structure { of: Mixture::Types } => {
                first_found(&[&c.posterior, &c.membership], NO_TYPES).map(drop)
            }
            Panel::Structure {
                of: Mixture::Topics,
            } => found(&c.topics),
            Panel::Heatmap {
                of: HeatKind::Markers,
            } => found(&c.expression).and_then(|()| found(&c.markers)),
            Panel::Heatmap {
                of: HeatKind::Terms,
            } => found(&c.alignment),
            Panel::Heatmap {
                of: HeatKind::Types,
            } => first_found(&[&c.type_q, &c.membership], NO_TYPE_Q).map(drop),
            _ => Err(format!("{} is not a chart", panel.title())),
        };
        why.err()
    }

    /// The chart `panel` on the paper page at `size`; without its table, a
    /// line saying so. With `keep_text`, its text is kept as text, for a
    /// PDF.
    pub(crate) fn chart_frame(&self, panel: Panel, size: Size, keep_text: bool) -> Frame {
        let mut frame = Frame::blank(size.w as usize, size.h as usize, THEME.background);
        if keep_text {
            frame.keep_text();
        }
        let font = Font::for_cell_height(size.text_px, false);
        let bold = Font::for_cell_height(size.text_px, true);
        let missing = |frame: &mut Frame, e: &str| {
            note(frame, &format!("no {}: {e}", panel.title()), bold, &THEME);
        };
        match panel {
            Panel::Structure { of } => match self.mixture(of) {
                Ok(mix) => {
                    let palette = self.palette(of, &mix.names);
                    let layout = self.structure(&mix);
                    let weights = Weights {
                        k: mix.names.len(),
                        w: &mix.w,
                        valid: None,
                    };
                    let style = BarStyle {
                        palette: &palette,
                        focus: None,
                        gap: legume_tui::chart::structure::GAP,
                        hits: false,
                    };
                    let fig = Figure {
                        title: Some(of.title()),
                        key: Some(&mix.names),
                        names: Names::Auto,
                        font,
                        bold,
                        theme: &THEME,
                    };
                    figure(&mut frame, &layout, &weights, &style, &fig);
                }
                Err(e) => missing(&mut frame, &e),
            },
            Panel::Heatmap { of } => match self.heat(of) {
                Ok(h) => {
                    let heat = Heat {
                        title: &h.title,
                        rows: &h.rows,
                        cols: &h.cols,
                        col_swatch: None,
                        z: &h.z,
                        clip: CLIP,
                        bold_row: None,
                    };
                    let style = HeatStyle {
                        row_chars: 28,
                        ..HeatStyle::for_theme(&THEME)
                    };
                    // A heatmap with no room says so itself.
                    let _ = draw_heat(&mut frame, &heat, font, bold, &THEME, &style);
                }
                Err(e) => missing(&mut frame, &e),
            },
            _ => missing(&mut frame, "not a chart"),
        }
        frame
    }

    /// Colours of the components `names` stacks: a type in its colour on
    /// the map coloured by type, the rest past those.
    fn palette(&self, of: Mixture, names: &[String]) -> Vec<Lin> {
        let lin = |(r, g, b): legume_plot::palette::Rgb| to_linear([r, g, b]);
        let on_map = match of {
            Mixture::Types => &self.cell_colours(Colouring::Type).1[..],
            Mixture::Topics => &[],
        };
        let more = distinct_colours(on_map.len() + names.len());
        (names.iter().enumerate())
            .map(|(j, n)| match on_map.iter().find(|(m, _)| m == n) {
                Some(&(_, c)) => lin(c),
                None => lin(more[on_map.len() + j]),
            })
            .collect()
    }

    /// The chart `panel` as a PDF at `path`, as the map is saved: the
    /// chart as one image of `size` at `dpi`, its text set over it as
    /// text. The image, for a thumbnail.
    pub(crate) fn write_chart(
        &self,
        panel: Panel,
        size: Size,
        dpi: f32,
        path: &Path,
    ) -> Result<image::RgbaImage> {
        let mut frame = self.chart_frame(panel, size, true);
        let img = frame.to_image();
        let page = Page {
            img: img.clone(),
            text: frame.take_text(),
        };
        write_pdf(&[page], dpi, path).with_context(|| format!("writing {}", path.display()))?;
        Ok(img)
    }
}

/// The marker panel at `path` as each type's genes, types and genes in
/// order, as the round reads it.
fn read_panel(path: &str) -> Result<BTreeMap<String, BTreeSet<String>>> {
    let mut out: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for (g, t) in crate::annotate::markers::read_marker_pairs(path)? {
        out.entry(t.into_string())
            .or_default()
            .insert(g.into_string());
    }
    Ok(out)
}

/// A cells × components table at `path` on the round's `cells`, by name.
fn read_cells(path: &str, cells: &[Box<str>]) -> Result<Mix> {
    let t = Mat::from_parquet_with_row_names(path, Some(0))
        .with_context(|| format!("reading {path}"))?;
    let on = super::aligned(&t.rows, &t.mat, &cell_index(cells), cells.len());
    let k = on.ncols();
    let w = (0..on.nrows())
        .flat_map(|i| (0..k).map(move |j| (i, j)))
        .map(|(i, j)| on[(i, j)])
        .collect();
    Ok(Mix {
        names: t.cols.iter().map(ToString::to_string).collect(),
        w,
    })
}

/// Rows of `k` logs (a topic run's log θ) as shares summing to one; a row
/// with a missing value stays missing.
fn shares_of_logs(w: &mut [f32], k: usize) {
    for row in w.chunks_mut(k.max(1)) {
        let top = row.iter().copied().fold(f32::MIN, f32::max);
        if !row.iter().all(|v| v.is_finite()) {
            row.fill(f32::NAN);
            continue;
        }
        row.iter_mut().for_each(|v| *v = (*v - top).exp());
        let sum: f32 = row.iter().sum();
        row.iter_mut().for_each(|v| *v /= sum);
    }
}

#[cfg(test)]
#[path = "charts_tests.rs"]
pub(crate) mod tests;
