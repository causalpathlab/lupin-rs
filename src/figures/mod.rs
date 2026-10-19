//! Figures of a run's cells, from the outputs its manifest records: its 2D
//! layouts coloured by what the run says of each cell, and, for a trajectory
//! run (`docs/trajectory-plan.md` §6), the prior's edges as arrows, the
//! diffusion map, the types in pseudotime order and the connectivity between
//! types. A scatter of cells is a map drawn as senna draws its own, through
//! `legume_tui::map` on senna's paper page, and exported as a PDF of its
//! raster with the text set over it; the order and connectivity panels are
//! SVGs that legume-plot writes as figure files and resvg rasterises for the
//! terminal. The TUI's figure pane is where they are shown and exported.

use crate::manifest::rounds::{Cells, ScoreCells};
use crate::manifest::run::{derive_out_prefix, resolve, Loaded, RunManifest};
use crate::trajectory::edges::{self, EdgeRow, Verdict};
use anyhow::{Context, Result};
use enrichment::UNASSIGNED_LABEL;
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::parquet::{peek_parquet_field_names, read_table_columns};
use legume_numeric::matrix::traits::IoOps;
use legume_numeric::matrix::utils::{median, quantiles};
use legume_plot::hinton::{hinton_size, render_hinton, HintonOpts};
use legume_plot::palette::{sample_blue_red, Rgb};
use legume_plot::svg_emit::escape_xml;
use legume_plot::FigureFormats;
use legume_tui::export::{write_pdf, Page};
use legume_tui::map::colour::{to_linear, Lin, Srgb};
use legume_tui::map::labels::{place, Label, Placed, Rect};
use legume_tui::map::legend::draw_legend;
use legume_tui::map::offsets::Offsets;
use legume_tui::map::points::{Job, Mark, PointLayer};
use legume_tui::map::raster::{edge, Frame, Shape};
use legume_tui::map::text::{Font, TextRun};
use legume_tui::map::theme::Theme;
use legume_tui::map::viewport::{Bounds, Limits, Viewport, YAxis};
use rustc_hash::FxHashMap;
use std::borrow::Cow;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};

/// Pairs outside the prior on the connectivity panel, and the pairs it
/// orders.
const GREY: Rgb = (190, 190, 190);
const INK: Rgb = (40, 40, 40);
/// The page a cell map is drawn on, and its palette: senna's.
const THEME: Theme = Theme::PAPER;
/// How much of the frame a whole cell map fills, as senna's.
const MAP_MARGIN: f32 = 0.92;
/// A prior edge's arrow on a cell map, and an inferred edge's: the data's
/// edge, not the prior's.
const ARROW: Srgb = [40, 40, 40];
const INFERRED: Srgb = [115, 115, 115];
/// Entries a categorical legend lists before `+k more`.
const LEGEND_MAX: usize = 20;
/// Figure width in inches and dots per inch for an export.
const EXPORT_WIDTH_IN: f32 = 7.0;
const EXPORT_DPI: f32 = 200.0;

/// Label sizes `t` steps through, as `senna view` does: small, medium,
/// large, largest, then off.
pub(crate) const TEXT_SCALES: [f32; 4] = [1.0, 1.4, 1.8, 2.4];
const TEXT_SIZES: [&str; 4] = ["small", "medium", "large", "largest"];
/// The label states `t` cycles through: each size, then off.
const LABELS: [Option<usize>; 5] = [Some(0), Some(1), Some(2), Some(3), None];

/// The value after `cur` in `all`, wrapping round; the first when `cur` is
/// not in it.
pub(crate) fn next_in<T: Copy + PartialEq>(all: &[T], cur: T) -> T {
    let at = all.iter().position(|&v| v == cur);
    all[at.map_or(0, |k| (k + 1) % all.len())]
}

/// What the cells of a scatter are coloured by.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Colouring {
    Pseudotime,
    Type,
    Lineage,
    Component,
    /// A round's clusters.
    Cluster,
    /// The projection score's per-cell call.
    Call,
    /// The projection score's margin: how far its call leads the runner-up.
    Margin,
}

impl Colouring {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Pseudotime => "pseudotime",
            Self::Type => "cell type",
            Self::Lineage => "lineage",
            Self::Component => "component",
            Self::Cluster => "cluster",
            Self::Call => "score call",
            Self::Margin => "score margin",
        }
    }

    /// Its slot in [`CellFigures`]'s colour cache.
    fn slot(self) -> usize {
        self as usize
    }
}

/// How the scatters are drawn: the label size (an index into
/// [`TEXT_SCALES`], `None` for no labels), the colouring, and the cluster in
/// focus, whose cells keep their colours over the others in grey.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Style {
    pub(crate) labels: Option<usize>,
    pub(crate) colouring: Colouring,
    pub(crate) focus: Option<i32>,
}

impl Default for Style {
    fn default() -> Self {
        Self {
            labels: Some(1),
            colouring: Colouring::Pseudotime,
            focus: None,
        }
    }
}

impl Style {
    /// The style a figure pane opens in: the first colouring `data` has.
    pub(crate) fn opening(data: &CellFigures) -> Self {
        Self {
            colouring: data.colourings()[0],
            ..Self::default()
        }
    }

    /// Labels in turn: small, medium, large, largest, off (`t`, as
    /// `senna view` does); the note for the status line.
    pub(crate) fn cycle_labels(&mut self) -> String {
        self.labels = next_in(&LABELS, self.labels);
        match self.labels {
            Some(k) => {
                let then = next_in(&LABELS, self.labels).map_or("off", |j| TEXT_SIZES[j]);
                format!("labels {} · t for {then}", TEXT_SIZES[k])
            }
            None => "labels off · t shows them small".into(),
        }
    }

    /// The next of the colourings `all` (`c`); the note for the status line.
    pub(crate) fn cycle_colouring(&mut self, all: &[Colouring]) -> String {
        self.colouring = next_in(all, self.colouring);
        format!(
            "coloured by {} · c for {}",
            self.colouring.name(),
            next_in(all, self.colouring).name()
        )
    }
}

/// How far in `+` can zoom a scatter, in steps of ×√2 (level 2 is ×2, 12
/// is ×64), and how far an arrow pans, as a share of the part on screen.
const MAX_ZOOM_LEVEL: u8 = 12;
const PAN_STEP: f32 = 0.2;

/// The part of a scatter on screen: the whole extent zoomed `level` steps
/// about the centre `(cx, cy)`, both in shares of the whole extent (y up).
/// The level is a count, so zooming in and back out lands on whole exactly.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct View {
    pub(crate) level: u8,
    pub(crate) cx: f32,
    pub(crate) cy: f32,
}

impl Default for View {
    fn default() -> Self {
        Self {
            level: 0,
            cx: 0.5,
            cy: 0.5,
        }
    }
}

impl View {
    pub(crate) fn is_whole(&self) -> bool {
        self.level == 0
    }

    /// How many times the whole extent is magnified.
    pub(crate) fn factor(&self) -> f32 {
        2f32.powf(f32::from(self.level) / 2.0)
    }

    /// Zoom in (`inward`) or out by one step about the centre.
    pub(crate) fn zoom(&mut self, inward: bool) {
        self.level = if inward {
            (self.level + 1).min(MAX_ZOOM_LEVEL)
        } else {
            self.level.saturating_sub(1)
        };
        self.clamp();
    }

    /// Move the part on screen by `(dx, dy)` steps (right and up positive).
    pub(crate) fn pan(&mut self, dx: f32, dy: f32) {
        let span = 1.0 / self.factor();
        self.cx += dx * PAN_STEP * span;
        self.cy += dy * PAN_STEP * span;
        self.clamp();
    }

    /// Keep the part on screen inside the whole extent.
    fn clamp(&mut self) {
        let half = 0.5 / self.factor();
        self.cx = self.cx.clamp(half, 1.0 - half);
        self.cy = self.cy.clamp(half, 1.0 - half);
    }

    /// The camera on `whole` for a `w × h` frame: fitted as senna fits a
    /// map, zoomed and centred as the view says, and kept from showing past
    /// the data (`Limits.keep_inside`), so a zoomed map never pans off it.
    pub(crate) fn viewport(&self, whole: Bounds, w: usize, h: usize) -> Viewport {
        let mut vp = Viewport::fit(whole, w, h, MAP_MARGIN, YAxis::Up);
        let [x0, y0, x1, y1] = whole;
        vp.scale *= self.factor();
        vp.cx = x0 + self.cx * (x1 - x0);
        vp.cy = y0 + self.cy * (y1 - y0);
        let keep = Limits {
            scale: None,
            keep_inside: Some(whole),
        };
        vp.pan_px(0.0, 0.0, &keep);
        vp
    }
}

/// Which figure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Panel {
    /// One of the run's 2D layouts (an index into `CellFigures::layouts`),
    /// the prior's edges as arrows.
    Layout { k: usize },
    /// Two diffusion components.
    Diffusion { x: usize, y: usize },
    /// Median pseudotime per type, with the middle half of its cells.
    Order,
    /// PAGA connectivity between the node types, in pseudotime order.
    Connectivity,
}

impl Panel {
    /// Whether it is a scatter of cells, which the style colours and labels
    /// and which pans and zooms.
    pub(crate) fn is_scatter(self) -> bool {
        matches!(self, Self::Layout { .. } | Self::Diffusion { .. })
    }

    /// The figure's name in file names and the export log.
    pub(crate) fn slug(self) -> String {
        match self {
            Self::Layout { .. } => "layout".into(),
            Self::Diffusion { x, y } => format!("diffusion_dc{x}_dc{y}"),
            Self::Order => "order".into(),
            Self::Connectivity => "connectivity".into(),
        }
    }

    pub(crate) fn title(self) -> String {
        match self {
            Self::Layout { .. } => "layout".into(),
            Self::Diffusion { x, y } => format!("diffusion map · DC{x} × DC{y}"),
            Self::Order => "types by median pseudotime".into(),
            Self::Connectivity => "connectivity between types".into(),
        }
    }
}

/// A drawn figure: the SVG and its size, which a Hinton diagram sets itself.
pub(crate) struct Figure {
    pub(crate) svg: String,
    pub(crate) w: u32,
    pub(crate) h: u32,
}

/// A figure as drawn: a scatter as a cell map, the other panels as SVG.
pub(crate) enum Picture {
    Map(Box<CellMap>),
    Svg(Figure),
}

/// The size a figure is drawn at: its pixels, and the height of a line of
/// its text (a terminal cell's on screen, as senna sizes its labels).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Size {
    pub(crate) w: u32,
    pub(crate) h: u32,
    pub(crate) text_px: f32,
}

/// A scatter drawn as a map: the frame, the camera it was drawn with, and
/// where its labels and legend went (which the tests read, and a click on a
/// label will).
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) struct CellMap {
    pub(crate) frame: Frame,
    pub(crate) vp: Viewport,
    /// The labels placed, keyed by their place in `names`.
    pub(crate) placed: Vec<Placed<u32>>,
    pub(crate) names: Vec<Box<str>>,
    pub(crate) legend: Option<Rect>,
    /// The legend's and the labels' text, when it was kept as text (for a
    /// PDF) rather than drawn into the frame.
    pub(crate) text: Vec<TextRun>,
}

/// The cells of a scatter as the map draws them: each in its colour, grey
/// without one, in the order the figure draws them.
pub(crate) struct CellLayer<'a> {
    x: &'a [f32],
    y: &'a [f32],
    /// The cells with coordinates, in a fixed random order.
    shuffled: Vec<u32>,
    order: Arc<Vec<u32>>,
    colour: Cow<'a, [Option<Rgb>]>,
}

impl PointLayer for CellLayer<'_> {
    fn xy(&self, i: usize) -> [f32; 2] {
        [self.x[i], self.y[i]]
    }

    fn shuffled(&self) -> &[u32] {
        &self.shuffled
    }

    fn order(&self) -> Arc<Vec<u32>> {
        Arc::clone(&self.order)
    }

    fn mark(&self, i: usize) -> Option<Mark> {
        let colour = match self.colour.get(i).copied().flatten() {
            Some((r, g, b)) => to_linear([r, g, b]),
            None => to_linear(THEME.no_group),
        };
        Some(Mark {
            colour,
            alpha: 1.0,
            shape: Shape::Circle,
            size: 1.0,
        })
    }
}

/// A type's pseudotime summary, one row of the order panel.
pub(crate) struct TypeOrder {
    pub(crate) name: Box<str>,
    /// Cells with a finite pseudotime.
    pub(crate) cells: usize,
    pub(crate) q1: f32,
    pub(crate) median: f32,
    pub(crate) q3: f32,
}

/// An export: its base name, the files written and a line for the log.
pub(crate) struct Export {
    pub(crate) base: PathBuf,
    pub(crate) files: Vec<PathBuf>,
    pub(crate) what: String,
}

/// A run's per-cell outputs, aligned to its cells: a trajectory run's
/// pseudotime, types, lineages, components, diffusion map and edges, or a
/// round's labels (as `types`) and clusters. What a run lacks stays empty.
#[derive(Default)]
pub(crate) struct CellFigures {
    /// What made the outputs, which names the figures' exports.
    pub(crate) source: Source,
    /// The manifest the outputs were read from; export names derive from it.
    pub(crate) manifest: PathBuf,
    /// The prefix of the trajectory the outputs are of, when another
    /// manifest's (a round keeps the section of the run it was made on).
    pub(crate) trajectory_of: Option<String>,
    /// Per cell; NaN for cells no root reaches.
    pub(crate) pseudotime: Vec<f32>,
    pub(crate) types: Vec<Box<str>>,
    /// Per cell: its prior component, -1 for none.
    pub(crate) component: Vec<i32>,
    /// Per cell: the one lineage its type lies on, -1 for a type shared by
    /// several lineages or on none.
    pub(crate) lineage: Vec<i32>,
    /// Cells × diffusion components.
    pub(crate) diffusion: Option<Mat>,
    /// The run's layouts that share cells with it, best first.
    pub(crate) layouts: Vec<Layout>,
    /// Per cell: a round's cluster, -1 for none.
    pub(crate) cluster: Vec<i32>,
    /// A round's projection score per cell (its call and margin), when the
    /// round carries one.
    pub(crate) score: Option<ScoreCells>,
    /// What to run for the layouts of [`WANTED_LAYOUTS`] it lacks, if any.
    pub(crate) layout_hint: Option<String>,
    pub(crate) edges: Vec<EdgeRow>,
    /// Each colouring's cell colours and legend, made on first use.
    pub(crate) colours: [OnceLock<CellColours>; 7],
    /// Every cell, grey first, then by pseudotime: the draw order, made on
    /// first use.
    pub(crate) by_pseudotime: OnceLock<Vec<usize>>,
    /// The share of the cells a diffusion component is spread over below
    /// which the run left it out of pseudotime; the figure opens on the
    /// first pair it kept.
    pub(crate) dc_min_share: f32,
    /// Each scatter's label medians over all its cells, made on first use.
    pub(crate) medians: Mutex<Vec<(MediansKey, Arc<Medians>)>>,
}

/// Each cell's colour (`None` for grey) and the legend of a categorical
/// colouring.
type CellColours = (Vec<Option<Rgb>>, Vec<(String, Rgb)>);
/// What made a [`CellFigures`]' outputs.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum Source {
    #[default]
    Trajectory,
    Round,
}

/// A cluster's name on figures, as the TUI and `review` name it.
pub(crate) fn cluster_name(k: impl std::fmt::Display) -> String {
    format!("K{k}")
}

/// Each cell's colour on a scatter and the legend of those colours.
type Shown<'a> = (Cow<'a, [Option<Rgb>]>, Cow<'a, [(String, Rgb)]>);

/// A scatter's x and y per cell, borrowed from a layout or made from the
/// diffusion map.
type Coords<'a> = (Cow<'a, [f32]>, Cow<'a, [f32]>);

/// Each label's median position on a scatter, and its count of cells.
type Medians = BTreeMap<Box<str>, ((f32, f32), usize)>;

/// The scatter and labelling a set of medians is of.
type MediansKey = (Panel, LabelBy);

/// What a scatter's labels name: the cells' types or their clusters.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum LabelBy {
    Type,
    Cluster,
    Call,
}

/// A 2D layout of the run's cells, NaN where a cell has none.
pub(crate) struct Layout {
    /// Which of senna's layouts it is (`phate`, `umap`, …), when known.
    pub(crate) method: Option<String>,
    pub(crate) x: Vec<f32>,
    pub(crate) y: Vec<f32>,
}

impl CellFigures {
    /// The outputs `manifest` (read from `file`, in `dir`) records; `None`
    /// when it has no `trajectory.pseudotime`.
    pub(crate) fn load(manifest: &RunManifest, dir: &Path, file: &Path) -> Result<Option<Self>> {
        let t = &manifest.trajectory;
        let Some(pt_rel) = t.pseudotime.as_deref() else {
            return Ok(None);
        };
        let at = |rel: &str| resolve(dir, rel);
        let pt_path = at(pt_rel);
        let (cells, types, pseudotime, component, lineage) =
            read_pseudotime(&pt_path).with_context(|| format!("reading {pt_path}"))?;
        let index = cell_index(&cells);

        let diffusion = match t.diffusion.as_deref() {
            Some(rel) => {
                let d = Mat::from_parquet_with_row_names(&at(rel), Some(0))?;
                Some(aligned(&d.rows, &d.mat, &index, cells.len()))
            }
            None => None,
        };
        let (all, layout_hint) = read_layouts(manifest, dir, file, &index, cells.len());
        let edges = match t.edges.as_deref() {
            Some(rel) => edges::read(&at(rel))?,
            None => Vec::new(),
        };
        Ok(Some(Self {
            manifest: file.to_path_buf(),
            trajectory_of: crate::manifest::family::borrowed_trajectory(manifest, file),
            pseudotime,
            types,
            component,
            lineage,
            diffusion,
            layouts: all,
            layout_hint,
            // What the run left out with; runs from before it was recorded
            // kept every component, but open on the spread ones all the same.
            dc_min_share: t
                .settings
                .as_ref()
                .and_then(|s| s.get("dc_min_share")?.as_f64())
                .map_or(crate::trajectory::diffusion::DC_MIN_SHARE, |v| v as f32),
            edges,
            ..Self::default()
        }))
    }

    /// A round's labels and clusters (`cells`, read from `loaded`) on its
    /// run's layouts; none when the run has no layout, with what to run for
    /// one in `layout_hint`.
    pub(crate) fn of_round(loaded: &Loaded, cells: &Cells, score: Option<ScoreCells>) -> Self {
        let index = cell_index(&cells.names);
        let (layouts, layout_hint) = read_layouts(
            &loaded.manifest,
            &loaded.dir,
            &loaded.file,
            &index,
            cells.names.len(),
        );
        Self {
            source: Source::Round,
            manifest: loaded.file.clone(),
            types: cells
                .labels
                .iter()
                .map(|l| l.as_deref().unwrap_or(UNASSIGNED_LABEL).into())
                .collect(),
            cluster: cells
                .clusters
                .iter()
                .map(|k| k.map_or(-1, |k| k as i32))
                .collect(),
            layouts,
            layout_hint,
            score,
            ..Self::default()
        }
    }

    /// What the figures are of, in export names: `trajectory` for a
    /// trajectory run's, `cells` for a round's.
    pub(crate) fn kind(&self) -> &'static str {
        match self.source {
            Source::Trajectory => "trajectory",
            Source::Round => "cells",
        }
    }

    /// The run's output prefix, which export names start from.
    pub(crate) fn prefix(&self) -> String {
        derive_out_prefix(&self.manifest.to_string_lossy())
    }

    /// A panel's title, naming the layout drawn and, on a scatter, the
    /// colouring and the cluster in focus.
    pub(crate) fn title(&self, panel: Panel, style: &Style) -> String {
        let focus = style
            .focus
            .map_or(String::new(), |k| format!(" · {}", cluster_name(k)));
        let colouring = format!("{}{focus}", style.colouring.name());
        match panel {
            Panel::Layout { k } => {
                let method = self.layouts.get(k).and_then(|l| l.method.as_deref());
                format!(
                    "{} · {colouring}",
                    method.map_or("layout".into(), str::to_uppercase)
                )
            }
            Panel::Diffusion { .. } => format!("{} · {colouring}", panel.title()),
            _ => panel.title(),
        }
    }

    /// The panel's name in file names and the export log, naming the layout.
    pub(crate) fn slug(&self, panel: Panel) -> String {
        match panel {
            Panel::Layout { k } => match self.layouts.get(k).and_then(|l| l.method.as_deref()) {
                Some(m) => format!("layout_{m}"),
                None => panel.slug(),
            },
            _ => panel.slug(),
        }
    }

    /// The coordinates a scatter can show, in the order `m` steps through
    /// them: the layouts, then the diffusion map at its first pair spread
    /// over the cells.
    pub(crate) fn scatters(&self) -> Vec<Panel> {
        let mut out: Vec<Panel> = (0..self.layouts.len())
            .map(|k| Panel::Layout { k })
            .collect();
        if let Some(d) = self.diffusion.as_ref().filter(|d| d.ncols() >= 3) {
            out.push(first_spread_pair(d, self.dc_min_share));
        }
        out
    }

    /// The scatter after `panel` in [`Self::scatters`]; `None` when there is
    /// only one, or `panel` is not a scatter.
    pub(crate) fn next_scatter(&self, panel: Panel) -> Option<Panel> {
        let all = self.scatters();
        if all.len() < 2 {
            return None;
        }
        let at = match panel {
            Panel::Layout { k } => k,
            Panel::Diffusion { .. } => all.len() - 1,
            _ => return None,
        };
        Some(all[(at + 1) % all.len()])
    }

    /// The colourings the outputs support, in the order `c` steps through
    /// them; a lineage or component colouring needs at least two.
    pub(crate) fn colourings(&self) -> Vec<Colouring> {
        let several = |v: &[i32]| {
            let mut seen = v.iter().filter(|&&c| c >= 0);
            seen.next().is_some_and(|&a| seen.any(|&b| b != a))
        };
        let mut out = Vec::new();
        if !self.pseudotime.is_empty() {
            out.push(Colouring::Pseudotime);
        }
        out.push(Colouring::Type);
        if several(&self.cluster) {
            out.push(Colouring::Cluster);
        }
        if self.score.is_some() {
            out.extend([Colouring::Call, Colouring::Margin]);
        }
        if several(&self.lineage) {
            out.push(Colouring::Lineage);
        }
        if several(&self.component) {
            out.push(Colouring::Component);
        }
        out
    }

    /// The figures these outputs allow, in display order: one scatter (the
    /// best layout, else the diffusion map), the order and the connectivity.
    pub(crate) fn panels(&self) -> Vec<Panel> {
        let mut out: Vec<Panel> = self.scatters().into_iter().take(1).collect();
        if !self.pseudotime.is_empty() {
            out.push(Panel::Order);
        }
        if !self.edges.is_empty() {
            out.push(Panel::Connectivity);
        }
        out
    }

    /// The next diffusion pair after `(x, y)`: the y axis (or the x axis when
    /// `x_axis`) steps through the components after the trivial first one,
    /// skipping the other axis's.
    pub(crate) fn next_pair(
        &self,
        x: usize,
        y: usize,
        forward: bool,
        x_axis: bool,
    ) -> (usize, usize) {
        let n = self.diffusion.as_ref().map_or(0, Mat::ncols);
        if n < 3 {
            return (x, y);
        }
        let step = |v: usize| {
            if forward {
                if v + 1 < n {
                    v + 1
                } else {
                    1
                }
            } else if v > 1 {
                v - 1
            } else {
                n - 1
            }
        };
        let (moving, other) = if x_axis { (x, y) } else { (y, x) };
        let mut next = step(moving);
        if next == other {
            next = step(next);
        }
        if x_axis {
            (next, y)
        } else {
            (x, next)
        }
    }

    /// The figure, drawn to fit `size`; a scatter in `style` at `view` (a
    /// thumbnail shows it whole).
    pub(crate) fn figure(
        &self,
        panel: Panel,
        size: Size,
        style: &Style,
        view: View,
    ) -> Result<Picture> {
        Ok(match panel {
            Panel::Layout { .. } | Panel::Diffusion { .. } => {
                Picture::Map(Box::new(self.cell_map(panel, size, style, view, false)?))
            }
            Panel::Order => Picture::Svg(self.order_figure(size.w, size.h)),
            Panel::Connectivity => Picture::Svg(self.connectivity_figure(size.w, size.h)),
        })
    }

    /// The coordinates of the scatter `panel`.
    fn coords(&self, panel: Panel) -> Result<Coords<'_>> {
        match panel {
            Panel::Layout { k } => {
                let l = self.layouts.get(k).context("the run has no layout")?;
                Ok((Cow::Borrowed(&l.x[..]), Cow::Borrowed(&l.y[..])))
            }
            Panel::Diffusion { x, y } => {
                let d = self.diffusion.as_ref().context("no diffusion map")?;
                let column = |j: usize| d.column(j).iter().copied().collect::<Vec<f32>>();
                Ok((Cow::Owned(column(x)), Cow::Owned(column(y))))
            }
            _ => anyhow::bail!("{} is not a scatter", panel.title()),
        }
    }

    /// Write `panel` as `{prefix}.{kind}.{panel}`, moved past existing
    /// files as `-2`, `-3` …: a scatter as it is on screen, as a PDF of its
    /// raster with its text set over it (as senna exports its map); the
    /// other panels as `.svg` and `.pdf`.
    pub(crate) fn export(&self, panel: Panel, style: &Style, view: View) -> Result<Export> {
        let base = free_base(&format!(
            "{}.{}.{}",
            self.prefix(),
            self.kind(),
            self.slug(panel)
        ));
        let w = (EXPORT_WIDTH_IN * EXPORT_DPI) as u32;
        // 11 pt text at the export's resolution.
        let size = Size {
            w,
            h: w * 3 / 4,
            text_px: EXPORT_DPI * 11.0 / 72.0,
        };
        let b = base.to_string_lossy();
        let files = if panel.is_scatter() {
            let mut m = self.cell_map(panel, size, style, view, true)?;
            let pdf = PathBuf::from(format!("{b}.pdf"));
            let page = Page {
                img: m.frame.to_image(),
                text: std::mem::take(&mut m.text),
            };
            write_pdf(&[page], EXPORT_DPI, &pdf)
                .with_context(|| format!("writing {}", pdf.display()))?;
            vec![pdf]
        } else {
            let Picture::Svg(fig) = self.figure(panel, size, style, view)? else {
                unreachable!("only a scatter is a map")
            };
            let formats = FigureFormats {
                svg: true,
                png: false,
                pdf: true,
            };
            legume_plot::write_figure(&fig.svg, fig.w, fig.h, &b, formats)?;
            vec![
                PathBuf::from(format!("{b}.svg")),
                PathBuf::from(format!("{b}.pdf")),
            ]
        };
        let zoomed = if panel.is_scatter() && !view.is_whole() {
            format!(" · zoomed ×{:.1}", view.factor())
        } else {
            String::new()
        };
        Ok(Export {
            base,
            files,
            what: format!(
                "{}{zoomed} · {EXPORT_WIDTH_IN:.0} in · {EXPORT_DPI:.0} dpi",
                self.title(panel, style),
            ),
        })
    }

    /// The scatter `panel` as a map, as senna draws its own: the cells in
    /// the style's colouring (grey when they have no value) by the shared
    /// point job on the paper page; with the prior's direct edges as arrows
    /// between the type medians on a layout; the legend in the top-left
    /// corner; and the labels at the medians of their cells (of the cells on
    /// screen when zoomed), placed so none covers another. With `keep_text`
    /// the text is kept as text, for a PDF, instead of drawn.
    pub(crate) fn cell_map(
        &self,
        panel: Panel,
        size: Size,
        style: &Style,
        view: View,
        keep_text: bool,
    ) -> Result<CellMap> {
        let (x, y) = self.coords(panel)?;
        let (x, y) = (&x[..], &y[..]);
        let layer = self.cell_layer(x, y, style);
        anyhow::ensure!(!layer.shuffled.is_empty(), "no cell has coordinates");
        let whole = bounds_of(x, y, layer.shuffled.iter().map(|&i| i as usize));
        let vp = view.viewport(whole, size.w as usize, size.h as usize);
        let layers = std::slice::from_ref(&layer);
        let mut job = Job::new(vp, THEME.background, layers);
        while !job.step(layers) {}
        let mut frame = job.finish();
        let finite: Vec<usize> = layer.shuffled.iter().map(|&i| i as usize).collect();
        // The type medians over all cells anchor the arrows; only made when
        // they are drawn.
        if matches!(panel, Panel::Layout { .. }) && !self.edges.is_empty() {
            let all = self.medians(panel, LabelBy::Type, x, y, &finite);
            let stroke = (size.w.min(size.h) as f32 / 400.0).clamp(1.0, 3.0) * 1.2;
            self.draw_arrows(&mut frame, &all, &vp, stroke);
        }
        if keep_text {
            frame.keep_text();
        }
        let (_, legend) = self.shown_colours(style);
        let legend = self.draw_legend(&mut frame, &legend, size.text_px);
        let mut names = Vec::new();
        let mut labels = Vec::new();
        if let Some(k) = style.labels {
            let font = Font::for_cell_height(size.text_px * TEXT_SCALES[k], true);
            // Coloured by cluster, the labels name the clusters; else the
            // types. The whole scatter's are kept; a zoomed one's are of the
            // cells on screen.
            let by = match style.colouring {
                Colouring::Cluster => LabelBy::Cluster,
                Colouring::Call => LabelBy::Call,
                _ => LabelBy::Type,
            };
            let medians = if view.is_whole() {
                self.medians(panel, by, x, y, &finite)
            } else {
                let [x0, y0, x1, y1] = vp.window();
                let shown: Vec<usize> = finite
                    .iter()
                    .copied()
                    .filter(|&i| (x0..=x1).contains(&x[i]) && (y0..=y1).contains(&y[i]))
                    .collect();
                Arc::new(self.label_medians(by, x, y, &shown))
            };
            // A label in its group's ink when the cells are coloured by its
            // groups, as senna inks them; else in the page's.
            let (_, all) = self.cell_colours(style.colouring);
            let ink_of = |name: &str| {
                all.iter()
                    .find(|(n, _)| n == name)
                    .map_or(to_linear(THEME.ink), |&(_, (r, g, b))| {
                        THEME.ink_for(to_linear([r, g, b]))
                    })
            };
            for (t, &((mx, my), cells)) in medians.iter() {
                if t.as_ref() == UNASSIGNED_LABEL {
                    continue;
                }
                let (lx, ly) = keep_inside(vp.to_px([mx, my]), font, t, &frame);
                labels.push(Label {
                    text: t.to_string(),
                    x: lx,
                    y: ly,
                    ink: ink_of(t),
                    tier: 0,
                    priority: cells as f32,
                    font,
                    key: Some(names.len() as u32),
                });
                names.push(t.clone());
            }
        }
        let reserved: Vec<Rect> = legend.into_iter().collect();
        let placed = place(&mut frame, labels, &reserved, &Offsets::default());
        let text = frame.take_text();
        Ok(CellMap {
            frame,
            vp,
            placed,
            names,
            legend,
            text,
        })
    }

    /// The cells of `(x, y)` with coordinates as a map layer in the style's
    /// colouring and draw order.
    pub(crate) fn cell_layer<'a>(
        &'a self,
        x: &'a [f32],
        y: &'a [f32],
        style: &Style,
    ) -> CellLayer<'a> {
        let finite = (0..x.len()).filter(|&i| x[i].is_finite() && y[i].is_finite());
        let shuffled = shuffled(finite);
        let (colour, _) = self.shown_colours(style);
        let cells: Vec<usize> = shuffled.iter().map(|&i| i as usize).collect();
        let order = self
            .draw_order(&cells, style.colouring, &colour)
            .into_iter()
            .map(|i| i as u32)
            .collect();
        CellLayer {
            x,
            y,
            shuffled,
            order: Arc::new(order),
            colour,
        }
    }

    /// A categorical colouring's legend in the top-left corner, as senna
    /// draws it: a dot and the name in its group's ink per entry, at most
    /// [`LEGEND_MAX`] rows and no more than fill half the frame; none in a
    /// frame too small for two rows (a grid thumbnail).
    fn draw_legend(
        &self,
        frame: &mut Frame,
        legend: &[(String, Rgb)],
        text_px: f32,
    ) -> Option<Rect> {
        let font = Font::for_cell_height(text_px * 0.8, false);
        let rows = legend_rows(frame.h, font.line_height())?;
        let entries: Vec<(String, Lin, Lin)> = legend
            .iter()
            .map(|(name, (r, g, b))| {
                let c = to_linear([*r, *g, *b]);
                (name.clone(), c, THEME.ink_for(c))
            })
            .collect();
        draw_legend(frame, &entries, font, Some((rows, to_linear(THEME.ink))))
    }

    /// The order `cells` are drawn in, later over earlier: grey cells first,
    /// so coloured ones draw over them; by pseudotime, low to high, so the
    /// few late cells are not buried under the many early ones.
    fn draw_order(
        &self,
        cells: &[usize],
        colouring: Colouring,
        colour_of: &[Option<Rgb>],
    ) -> Vec<usize> {
        if colouring != Colouring::Pseudotime || self.pseudotime.is_empty() {
            return grey_first(cells.iter().copied(), colour_of);
        }
        // Every cell sorted once; each figure keeps the cells it shows.
        let all = self.by_pseudotime.get_or_init(|| {
            let t = &self.pseudotime;
            let mut all: Vec<usize> = (0..t.len()).collect();
            all.sort_by(|&a, &b| {
                (t[a].is_finite().cmp(&t[b].is_finite())).then(t[a].total_cmp(&t[b]))
            });
            all
        });
        let mut shown = vec![false; self.pseudotime.len()];
        for &i in cells {
            shown[i] = true;
        }
        grey_first(all.iter().copied().filter(|&i| shown[i]), colour_of)
    }

    /// Each cell's colour in the style's colouring and the legend: outside
    /// the cluster in focus, no colour (grey), and the legend of the colours
    /// left. A categorical colouring gives each entry its own colour.
    pub(crate) fn shown_colours(&self, style: &Style) -> Shown<'_> {
        let (all, legend) = self.cell_colours(style.colouring);
        let Some(k) = style.focus.filter(|_| !self.cluster.is_empty()) else {
            return (all.into(), legend.into());
        };
        let lit: Vec<Option<Rgb>> = all
            .iter()
            .zip(&self.cluster)
            .map(|(&c, &at)| c.filter(|_| at == k))
            .collect();
        // The distinct lit colours, until every entry's is found.
        let mut colours: Vec<Rgb> = Vec::new();
        for &c in lit.iter().flatten() {
            if colours.len() == legend.len() {
                break;
            }
            if !colours.contains(&c) {
                colours.push(c);
            }
        }
        let kept = legend
            .iter()
            .filter(|(_, c)| colours.contains(c))
            .cloned()
            .collect();
        (lit.into(), Cow::Owned(kept))
    }

    /// Each cell's colour in `colouring` and its legend, made once.
    fn cell_colours(&self, colouring: Colouring) -> &CellColours {
        self.colours[colouring.slot()].get_or_init(|| match colouring {
            // A round has no pseudotime: its cells are grey.
            Colouring::Pseudotime if self.pseudotime.is_empty() => {
                (vec![None; self.types.len()], Vec::new())
            }
            Colouring::Pseudotime => (
                self.pseudotime
                    .iter()
                    .map(|&t| t.is_finite().then(|| sample_blue_red(t)))
                    .collect(),
                Vec::new(),
            ),
            Colouring::Type => {
                let keys: Vec<Option<&str>> = self
                    .types
                    .iter()
                    .map(|t| Some(t.as_ref()).filter(|&t| t != UNASSIGNED_LABEL))
                    .collect();
                categorical(&keys, str::to_string)
            }
            Colouring::Lineage => categorical(&codes(&self.lineage), |c| format!("lineage {c}")),
            Colouring::Component => {
                categorical(&codes(&self.component), |c| format!("component {c}"))
            }
            Colouring::Cluster => categorical(&codes(&self.cluster), cluster_name),
            Colouring::Call | Colouring::Margin if self.score.is_none() => {
                (vec![None; self.types.len()], Vec::new())
            }
            Colouring::Call => categorical(
                &self
                    .score
                    .iter()
                    .flat_map(ScoreCells::calls)
                    .collect::<Vec<_>>(),
                str::to_string,
            ),
            // Blue for a call that barely leads, red for the widest margin.
            Colouring::Margin => {
                let margin = self.score.as_ref().map_or(&[][..], |s| &s.margin[..]);
                let widest = margin
                    .iter()
                    .copied()
                    .filter(|m| m.is_finite())
                    .fold(f32::MIN_POSITIVE, f32::max);
                let colours = margin
                    .iter()
                    .map(|&m| m.is_finite().then(|| sample_blue_red(m / widest)))
                    .collect();
                (colours, Vec::new())
            }
        })
    }

    /// The label medians (`by`) over all of `cells` on the scatter `panel`,
    /// kept for the next figure of it. They are made outside the lock, so figures
    /// drawn in parallel do not wait on each other.
    fn medians(
        &self,
        panel: Panel,
        by: LabelBy,
        x: &[f32],
        y: &[f32],
        cells: &[usize],
    ) -> Arc<Medians> {
        let key = (panel, by);
        let kept = |m: &Vec<(MediansKey, Arc<Medians>)>| {
            m.iter().find(|(k, _)| *k == key).map(|(_, m)| m.clone())
        };
        if let Some(m) = kept(&self.medians.lock().expect("not poisoned")) {
            return m;
        }
        let made = Arc::new(self.label_medians(by, x, y, cells));
        let mut all = self.medians.lock().expect("not poisoned");
        if let Some(m) = kept(&all) {
            return m;
        }
        all.push((key, made.clone()));
        made
    }

    /// The prior's direct edges as arrows between the type `medians`:
    /// supported ones solid, the rest faded, inferred ones in grey.
    fn draw_arrows(&self, frame: &mut Frame, medians: &Medians, vp: &Viewport, stroke: f32) {
        let drawn = |e: &&EdgeRow| e.in_prior || e.verdict == Some(Verdict::Inferred);
        for e in self.edges.iter().filter(drawn) {
            if let (Some(&(a, _)), Some(&(b, _))) = (medians.get(&e.a), medians.get(&e.b)) {
                let (ink, alpha) = match e.verdict {
                    Some(Verdict::Supported) => (ARROW, 0.9),
                    Some(Verdict::Inferred) => (INFERRED, 0.9),
                    _ => (ARROW, 0.3),
                };
                let (a, b) = (vp.to_px([a.0, a.1]), vp.to_px([b.0, b.1]));
                arrow(frame, a, b, stroke, to_linear(ink), alpha);
            }
        }
    }

    /// Each type's, or each cluster's (`by`), median position among `cells`.
    fn label_medians(&self, by: LabelBy, x: &[f32], y: &[f32], cells: &[usize]) -> Medians {
        fn medians_by<K: Ord>(
            cells: &[usize],
            key: impl Fn(usize) -> Option<K>,
            x: &[f32],
            y: &[f32],
        ) -> BTreeMap<K, ((f32, f32), usize)> {
            let mut by: BTreeMap<K, (Vec<f32>, Vec<f32>)> = BTreeMap::new();
            for &i in cells {
                if let Some(k) = key(i) {
                    let e = by.entry(k).or_default();
                    e.0.push(x[i]);
                    e.1.push(y[i]);
                }
            }
            by.into_iter()
                .map(|(k, (xs, ys))| (k, ((median(&xs), median(&ys)), xs.len())))
                .collect()
        }
        match by {
            LabelBy::Type | LabelBy::Call => {
                let calls = match (by, &self.score) {
                    (LabelBy::Call, Some(s)) => s.calls(),
                    (LabelBy::Call, None) => Vec::new(),
                    _ => self.types.iter().map(|t| Some(t.as_ref())).collect(),
                };
                medians_by(cells, |i| calls.get(i).copied().flatten(), x, y)
                    .into_iter()
                    .map(|(t, at)| (t.into(), at))
                    .collect()
            }
            LabelBy::Cluster => {
                let codes = codes(&self.cluster);
                medians_by(cells, |i| codes[i], x, y)
                    .into_iter()
                    .map(|(c, at)| (cluster_name(c).into(), at))
                    .collect()
            }
        }
    }

    /// The types with a finite pseudotime, in order of their median.
    pub(crate) fn type_order(&self) -> Vec<TypeOrder> {
        let mut by_type: BTreeMap<&str, Vec<f32>> = BTreeMap::new();
        for (t, &p) in self.types.iter().zip(&self.pseudotime) {
            if p.is_finite() {
                by_type.entry(t.as_ref()).or_default().push(p);
            }
        }
        let mut rows: Vec<TypeOrder> = by_type
            .into_iter()
            .map(|(t, v)| {
                let q = quantiles(&v, &[0.25, 0.5, 0.75]);
                TypeOrder {
                    name: t.into(),
                    cells: v.len(),
                    q1: q[0],
                    median: q[1],
                    q3: q[2],
                }
            })
            .collect();
        rows.sort_by(|a, b| a.median.total_cmp(&b.median));
        rows
    }

    fn order_figure(&self, w: u32, h: u32) -> Figure {
        let rows = self.type_order();
        let (left, right, top, bottom) = (w as f32 * 0.36, w as f32 - 16.0, 28.0, h as f32 - 12.0);
        let rh = ((bottom - top) / rows.len().max(1) as f32).min(26.0);
        let x = |v: f32| left + v * (right - left);
        let font = rh.clamp(8.0, 13.0) * 0.9;
        let mut s = format!(
            "<?xml version='1.0' encoding='UTF-8'?>\n<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {w} {h}' width='{w}' height='{h}' font-family='Helvetica, Arial, sans-serif'>\
             <rect width='{w}' height='{h}' fill='white'/>"
        );
        for tick in [0.0, 0.25, 0.5, 0.75, 1.0] {
            s += &format!(
                "<line x1='{0:.1}' y1='{1}' x2='{0:.1}' y2='{2:.1}' stroke='#ddd'/>\
                 <text x='{0:.1}' y='{3}' font-size='{4:.1}' text-anchor='middle' fill='#666'>{5}</text>",
                x(tick),
                top - 6.0,
                bottom,
                top - 10.0,
                font,
                tick
            );
        }
        for (k, r) in rows.iter().enumerate() {
            let y = top + rh * (k as f32 + 0.5);
            s += &format!(
                "<text x='{:.1}' y='{:.1}' font-size='{font:.1}' text-anchor='end' fill='#222'>{} ({})</text>\
                 <line x1='{:.1}' y1='{y:.1}' x2='{:.1}' y2='{y:.1}' stroke='#9ec5f4' stroke-width='{:.1}' stroke-linecap='round'/>\
                 <circle cx='{:.1}' cy='{y:.1}' r='{:.1}' fill='#2a78d6' stroke='white' stroke-width='1.5'/>",
                left - 8.0,
                y + font * 0.35,
                escape_xml(&r.name),
                r.cells,
                x(r.q1),
                x(r.q3),
                (rh * 0.28).max(2.0),
                x(r.median),
                (rh * 0.22).clamp(2.5, 5.0)
            );
        }
        s += "</svg>";
        Figure { svg: s, w, h }
    }

    /// Connectivity as a Hinton diagram (box area ∝ connectivity), the types
    /// in pseudotime order; pairs the prior orders are in ink, the rest grey.
    fn connectivity_figure(&self, w: u32, h: u32) -> Figure {
        // The types in pseudotime order, then any no root reaches.
        let mut nodes: Vec<Box<str>> = self.type_order().into_iter().map(|r| r.name).collect();
        for e in &self.edges {
            for t in [&e.a, &e.b] {
                if !nodes.contains(t) {
                    nodes.push(t.clone());
                }
            }
        }
        nodes.retain(|t| self.edges.iter().any(|e| e.a == *t || e.b == *t));
        let n = nodes.len().max(1);
        let at = |t: &str| nodes.iter().position(|x| x.as_ref() == t);
        let mut mat = vec![0.0f32; n * n];
        let mut colors = vec![GREY; n * n];
        for e in &self.edges {
            if let (Some(i), Some(j)) = (at(&e.a), at(&e.b)) {
                for k in [i * n + j, j * n + i] {
                    mat[k] = e.connectivity;
                    colors[k] = if e.in_prior { INK } else { GREY };
                }
            }
        }
        // Row labels take 9 font widths, the legend 6 cells (legume-plot's
        // layout); fit the grid into what is left.
        let font_px = (w.min(h) as f32 / 40.0).clamp(6.0, 12.0);
        let cell_px = ((w as f32 - font_px * 9.0) / (n as f32 + 6.0))
            .min((h as f32 - font_px * 8.0) / n as f32)
            .clamp(4.0, 40.0);
        let opts = HintonOpts {
            row_labels: Some(&nodes),
            col_labels: Some(&nodes),
            cell_colors: Some(&colors),
            cell_px,
            font_px,
            grid_stroke_px: 0.4,
            title: Some("connectivity between types"),
            ..HintonOpts::default()
        };
        let size = hinton_size(n, n, &opts);
        Figure {
            svg: render_hinton(&mat, n, n, &opts),
            w: size.width_px,
            h: size.height_px,
        }
    }
}

/// The bounding box of the points `idx` of `(x, y)`.
fn bounds_of(x: &[f32], y: &[f32], idx: impl Iterator<Item = usize>) -> Bounds {
    let (mut x0, mut x1, mut y0, mut y1) = (f32::MAX, f32::MIN, f32::MAX, f32::MIN);
    for i in idx {
        x0 = x0.min(x[i]);
        x1 = x1.max(x[i]);
        y0 = y0.min(y[i]);
        y1 = y1.max(y[i]);
    }
    [x0, y0, x1, y1]
}

/// `cells` in a fixed random order, the same on every draw: each prefix a
/// uniform subsample, as the map's point job wants them.
fn shuffled(cells: impl Iterator<Item = usize>) -> Vec<u32> {
    // splitmix64 of the cell's index: a permutation without a generator.
    let mix = |i: usize| {
        let mut z = (i as u64).wrapping_add(0x9e37_79b9_7f4a_7c15);
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    };
    let mut out: Vec<u32> = cells.map(|i| i as u32).collect();
    out.sort_by_cached_key(|&i| mix(i as usize));
    out
}

/// An arrow from `a` to `b`: a shaft `stroke` pixels wide and a two-stroke
/// head at `b`, each as thin edges side by side.
fn arrow(frame: &mut Frame, a: (f32, f32), b: (f32, f32), stroke: f32, ink: Lin, alpha: f32) {
    let (dx, dy) = (b.0 - a.0, b.1 - a.1);
    let len = dx.hypot(dy);
    if !len.is_finite() || len < 1.0 {
        return;
    }
    let (ux, uy) = (dx / len, dy / len);
    let lines = stroke.round().max(1.0) as usize;
    let mut thick = |p: (f32, f32), q: (f32, f32)| {
        let (vx, vy) = (q.0 - p.0, q.1 - p.1);
        let l = vx.hypot(vy).max(1e-6);
        let (nx, ny) = (-vy / l, vx / l);
        for k in 0..lines {
            let o = k as f32 - 0.5 * (lines - 1) as f32;
            let at = |(x, y): (f32, f32)| (x + o * nx, y + o * ny);
            edge(frame, at(p), at(q), (0.0, 0.0), ink, alpha);
        }
    };
    thick(a, b);
    let head = (stroke * 5.0).min(0.5 * len);
    for side in [-1.0, 1.0] {
        let back = (
            b.0 - head * ux + side * 0.5 * head * -uy,
            b.1 - head * uy + side * 0.5 * head * ux,
        );
        thick(b, back);
    }
}

/// The layouts senna makes that the trajectory's figures look for.
pub(crate) const WANTED_LAYOUTS: [&str; 2] = ["umap", "phate"];

/// The layouts of `manifest` (read from `file` in `dir`) and of the runs it
/// was made from, along `annotate.source`, as paths, PHATE first: a layout
/// senna adds to the run after a round or trajectory was made is found too.
/// Also the first run of that chain, which `senna layout` adds to.
pub(crate) fn layouts_along(
    manifest: &RunManifest,
    dir: &Path,
    file: &Path,
) -> (Vec<(Option<String>, String)>, PathBuf) {
    let mut out: Vec<(Option<String>, String)> = Vec::new();
    let mut add = |m: &RunManifest, dir: &Path| {
        for (method, rel) in layouts(m) {
            let path = resolve(dir, &rel);
            if !out.iter().any(|(_, p)| *p == path) {
                out.push((method, path));
            }
        }
    };
    add(manifest, dir);
    let mut run = file.to_path_buf();
    for l in crate::manifest::run::source_chain(manifest, dir) {
        add(&l.manifest, &l.dir);
        run = l.file;
    }
    out.sort_by_key(|(m, _)| m.as_deref() != Some("phate"));
    (out, run)
}

/// What to run for the [`WANTED_LAYOUTS`] missing from `found` (the
/// layouts' methods), on `run`.
pub(crate) fn layout_hint<'a>(
    found: impl IntoIterator<Item = Option<&'a str>>,
    run: &Path,
) -> Option<String> {
    let found: Vec<Option<&str>> = found.into_iter().collect();
    let missing: Vec<&str> = WANTED_LAYOUTS
        .into_iter()
        .filter(|w| !found.contains(&Some(*w)))
        .collect();
    if missing.is_empty() {
        return None;
    }
    let name = run
        .file_name()
        .map_or(run.to_string_lossy(), |n| n.to_string_lossy());
    let which = match missing[..] {
        [one] => one.to_string(),
        _ => format!("{{{}}}", missing.join("|")),
    };
    let shown: Vec<String> = missing.iter().map(|m| m.to_uppercase()).collect();
    Some(format!(
        "no {} for this run: `senna layout {which} --from {name}` adds it, and the figures pick it up when the run is opened again",
        shown.join(" or ")
    ))
}

/// The layouts recorded along `file`'s run that read onto its `n` cells
/// (`index`), and what to run for the ones it lacks.
fn read_layouts(
    manifest: &RunManifest,
    dir: &Path,
    file: &Path,
    index: &CellIndex,
    n: usize,
) -> (Vec<Layout>, Option<String>) {
    let mut all = Vec::new();
    let (found, run) = layouts_along(manifest, dir, file);
    for (method, path) in &found {
        match read_layout(path, index, n) {
            Ok((x, y)) => all.push(Layout {
                method: method.clone(),
                x,
                y,
            }),
            Err(e) => log::warn!("layout {path}: {e:#}"),
        }
    }
    // Shown on `m`; a run says it once, when it ends.
    let hint = layout_hint(all.iter().map(|l| l.method.as_deref()), &run);
    (all, hint)
}

/// The layouts to draw pseudotime on, best first: senna's PHATE (`senna
/// layout phate`, recorded under `layout.methods.phate`), as PHATE is built
/// to show trajectories, then the run's current layout, then the other
/// layouts it records. Each with its method's name, when known, and its
/// manifest-relative path.
fn layouts(manifest: &RunManifest) -> Vec<(Option<String>, String)> {
    let l = &manifest.layout;
    let methods = l.extra.get("methods").and_then(|m| m.as_object());
    let coords = |m: &str| {
        methods
            .and_then(|ms| ms.get(m)?.get("cell_coords")?.as_str())
            .map(str::to_string)
    };
    let current = l.extra.get("current").and_then(|v| v.as_str());
    let mut out: Vec<(Option<String>, String)> = Vec::new();
    let mut add = |method: Option<&str>, path: Option<String>| {
        if let Some(p) = path {
            if !out.iter().any(|(_, q)| *q == p) {
                out.push((method.map(str::to_string), p));
            }
        }
    };
    add(Some("phate"), coords("phate"));
    add(current, l.cell_coords.clone());
    for m in methods.into_iter().flat_map(|ms| ms.keys()) {
        add(Some(m), coords(m));
    }
    out
}

/// The pseudotime table at `path`, in one read: per cell its name, type,
/// pseudotime, component (-1 when the table has none) and the one lineage
/// its type lies on (-1 when shared or none), from the `L0`, `L1`, …
/// weights.
#[allow(clippy::type_complexity)]
fn read_pseudotime(
    path: &str,
) -> Result<(Vec<Box<str>>, Vec<Box<str>>, Vec<f32>, Vec<i32>, Vec<i32>)> {
    let fields = peek_parquet_field_names(path)?;
    let has = |name: &str| fields.iter().any(|f| f.as_ref() == name);
    let lineages: Vec<String> = (0..)
        .map(|k| format!("L{k}"))
        .take_while(|c| has(c))
        .collect();
    let with_component = has("component");
    let mut numeric = vec!["pseudotime"];
    if with_component {
        numeric.push("component");
    }
    numeric.extend(lineages.iter().map(String::as_str));
    let (strings, numbers) = read_table_columns(path, &["cell", "type"], &numeric)?;
    let [cells, types]: [Vec<Box<str>>; 2] = strings
        .try_into()
        .map_err(|_| anyhow::anyhow!("no cell and type columns"))?;
    let n = cells.len();
    let pseudotime = numbers[0].iter().map(|&v| v as f32).collect();
    let component = if with_component {
        numbers[1].iter().map(|&v| v as i32).collect()
    } else {
        vec![-1; n]
    };
    let weights = &numbers[1 + usize::from(with_component)..];
    let lineage = (0..n).map(|i| sole_lineage(weights, i)).collect();
    Ok((cells, types, pseudotime, component, lineage))
}

/// The lineage cell `i` weighs most on, -1 when it ties or weighs on none.
fn sole_lineage(weights: &[Vec<f64>], i: usize) -> i32 {
    let (mut best, mut at, mut tied) = (0.0, -1, false);
    for (k, w) in weights.iter().enumerate() {
        if w[i] > best {
            (best, at, tied) = (w[i], k as i32, false);
        } else if w[i] == best && best > 0.0 {
            tied = true;
        }
    }
    if tied {
        -1
    } else {
        at
    }
}

/// Each cell name's cell, and every cell of a name the cells repeat
/// (merged samples can repeat a barcode), apart: most names are one cell's.
pub(crate) struct CellIndex<'a> {
    first: FxHashMap<&'a str, usize>,
    repeated: FxHashMap<&'a str, Vec<usize>>,
}

impl CellIndex<'_> {
    fn contains(&self, name: &str) -> bool {
        self.first.contains_key(name)
    }
}

pub(crate) fn cell_index(cells: &[Box<str>]) -> CellIndex<'_> {
    let mut first: FxHashMap<&str, usize> = FxHashMap::default();
    let mut repeated: FxHashMap<&str, Vec<usize>> = FxHashMap::default();
    for (i, c) in cells.iter().enumerate() {
        match first.get(c.as_ref()) {
            None => {
                first.insert(c.as_ref(), i);
            }
            Some(&f) => repeated
                .entry(c.as_ref())
                .or_insert_with(|| vec![f])
                .push(i),
        }
    }
    CellIndex { first, repeated }
}

/// The layout at `path` as `(x, y)` for each of the `n` cells `index` maps
/// names onto, NaN for a cell the layout lacks; an error when it has none
/// of them, so the next layout is tried.
fn read_layout(path: &str, index: &CellIndex, n: usize) -> Result<(Vec<f32>, Vec<f32>)> {
    let t = Mat::from_parquet(path)?;
    let col = |name: &str| t.cols.iter().position(|c| c.as_ref() == name);
    let (Some(x), Some(y)) = (col("x"), col("y")) else {
        anyhow::bail!("{path} has no x and y columns");
    };
    if !t.rows.iter().any(|r| index.contains(r)) {
        anyhow::bail!("{path} shares no cell with the run's outputs");
    }
    let m = aligned(&t.rows, &t.mat, index, n);
    Ok((
        m.column(x).iter().copied().collect(),
        m.column(y).iter().copied().collect(),
    ))
}

/// `mat`'s rows reordered to the `n` cells by name; a cell the matrix
/// lacks is NaN. A name the cells repeat takes its rows in order when the
/// matrix repeats it as often, else its one row goes to each of its cells.
fn aligned(rows: &[Box<str>], mat: &Mat, index: &CellIndex, n: usize) -> Mat {
    let mut out = Mat::from_element(n, mat.ncols(), f32::NAN);
    let mut times: FxHashMap<&str, usize> = FxHashMap::default();
    for r in rows
        .iter()
        .filter(|r| index.repeated.contains_key(r.as_ref()))
    {
        *times.entry(r.as_ref()).or_default() += 1;
    }
    let mut seen: FxHashMap<&str, usize> = FxHashMap::default();
    for (i, r) in rows.iter().enumerate() {
        let r = r.as_ref();
        match index.repeated.get(r) {
            Some(cells) if times.get(r) == Some(&cells.len()) => {
                let k = seen.entry(r).or_default();
                out.set_row(cells[*k], &mat.row(i));
                *k += 1;
            }
            Some(cells) => cells.iter().for_each(|&j| out.set_row(j, &mat.row(i))),
            None => {
                if let Some(&j) = index.first.get(r) {
                    out.set_row(j, &mat.row(i));
                }
            }
        }
    }
    out
}

/// Colours for `keys` (one per cell, `None` for grey), a palette entry per
/// distinct key in order, and the legend naming them.
fn categorical<K: Ord + Copy>(keys: &[Option<K>], name: impl Fn(K) -> String) -> CellColours {
    let mut used: Vec<K> = keys.iter().flatten().copied().collect();
    used.sort_unstable();
    used.dedup();
    let colour = distinct_colours(used.len());
    let colours = keys
        .iter()
        .map(|k| {
            k.and_then(|k| used.binary_search(&k).ok())
                .map(|j| colour[j])
        })
        .collect();
    let legend = used
        .iter()
        .zip(&colour)
        .map(|(&k, &c)| (name(k), c))
        .collect();
    (colours, legend)
}

/// `n` different colours, from senna's palette: golden-angle hues at equal
/// lightness, alternating lightness past a dozen.
fn distinct_colours(n: usize) -> Vec<Rgb> {
    THEME
        .palette
        .srgb_all(n)
        .into_iter()
        .map(Rgb::from)
        .collect()
}

/// `cells` in order, the grey ones (no colour in `colour_of`) first.
fn grey_first(cells: impl Iterator<Item = usize> + Clone, colour_of: &[Option<Rgb>]) -> Vec<usize> {
    let coloured = |i: &usize| colour_of.get(*i).is_some_and(Option::is_some);
    let mut order: Vec<usize> = cells.clone().filter(|i| !coloured(i)).collect();
    order.extend(cells.filter(coloured));
    order
}

/// Lineage, component or cluster codes as keys, -1 (none) as `None`.
fn codes(v: &[i32]) -> Vec<Option<i32>> {
    v.iter().map(|&c| (c >= 0).then_some(c)).collect()
}

/// `xy` moved so a label of `text` in `font` centred on it stays inside
/// `frame`, instead of being dropped off its edge.
fn keep_inside(xy: (f32, f32), font: Font, text: &str, frame: &Frame) -> (f32, f32) {
    // The label's box, padded as `place` pads it.
    let half_w = 0.5 * font.width(text) as f32 + 3.0;
    let half_h = 0.5 * font.line_height() as f32 + 1.0;
    let inside = |v: f32, half: f32, len: f32| {
        if 2.0 * half >= len {
            0.5 * len
        } else {
            v.clamp(half, len - half)
        }
    };
    (
        inside(xy.0, half_w, frame.w as f32),
        inside(xy.1, half_h, frame.h as f32),
    )
}

/// Rows a legend may take in a frame `h` pixels high, with lines `line`
/// pixels high: at most [`LEGEND_MAX`] and half the frame, counting the
/// `+k more` line; `None` when fewer than two fit.
fn legend_rows(h: usize, line: usize) -> Option<usize> {
    let fits = (h / 2).saturating_sub(line / 2) / line.max(1);
    (fits >= 2).then_some(fits.min(LEGEND_MAX))
}

/// `base`, or `base-2`, `base-3`, … while any file of the set exists.
fn free_base(base: &str) -> PathBuf {
    let taken = |b: &str| {
        ["svg", "pdf", "png"]
            .iter()
            .any(|ext| Path::new(&format!("{b}.{ext}")).exists())
    };
    if !taken(base) {
        return PathBuf::from(base);
    }
    (2..)
        .map(|k| format!("{base}-{k}"))
        .find(|b| !taken(b))
        .map(PathBuf::from)
        .expect("some suffix is free")
}

/// `picture` as pixels: a cell map's frame as it is, an SVG through resvg.
pub(crate) fn render(picture: &Picture) -> Result<image::RgbaImage> {
    match picture {
        Picture::Map(m) => Ok(m.frame.to_image()),
        Picture::Svg(fig) => render_svg(fig),
    }
}

/// `fig` rendered to pixels by resvg, with the system's fonts.
fn render_svg(fig: &Figure) -> Result<image::RgbaImage> {
    use resvg::{tiny_skia, usvg};
    // The system fonts are loaded once, not for every figure drawn.
    static FONTS: std::sync::OnceLock<std::sync::Arc<usvg::fontdb::Database>> =
        std::sync::OnceLock::new();
    let fontdb = FONTS
        .get_or_init(|| {
            let mut db = usvg::fontdb::Database::new();
            db.load_system_fonts();
            std::sync::Arc::new(db)
        })
        .clone();
    let options = usvg::Options {
        fontdb,
        ..usvg::Options::default()
    };
    let tree = usvg::Tree::from_str(&fig.svg, &options).context("parsing the figure")?;
    let mut pixmap = tiny_skia::Pixmap::new(fig.w, fig.h)
        .with_context(|| format!("a {} × {} figure", fig.w, fig.h))?;
    resvg::render(
        &tree,
        tiny_skia::Transform::identity(),
        &mut pixmap.as_mut(),
    );
    let rgba = pixmap
        .pixels()
        .iter()
        .flat_map(|p| {
            let c = p.demultiply();
            [c.red(), c.green(), c.blue(), c.alpha()]
        })
        .collect();
    image::RgbaImage::from_raw(fig.w, fig.h, rgba).context("the figure's pixels")
}

/// The first two diffusion components after the trivial one spread over at
/// least `min_share` of the cells, as the run kept them: the components of a
/// barely attached group show that group against a line of every other cell.
/// DC1 × DC2 when fewer than two are.
fn first_spread_pair(d: &Mat, min_share: f32) -> Panel {
    use crate::trajectory::diffusion::participation;
    let spread: Vec<usize> = (1..d.ncols())
        .filter(|&j| {
            participation(d.column(j).iter().map(|&v| f64::from(v))) >= f64::from(min_share)
        })
        .take(2)
        .collect();
    match spread[..] {
        [x, y] => Panel::Diffusion { x, y },
        _ => Panel::Diffusion { x: 1, y: 2 },
    }
}

#[cfg(test)]
#[path = "tests.rs"]
pub(crate) mod tests;
