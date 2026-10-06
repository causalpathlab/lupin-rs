//! Figures of a trajectory run, from the outputs its manifest records
//! (`docs/trajectory-plan.md` §6): the layout coloured by pseudotime with the
//! prior's edges as arrows, the diffusion map, the types in pseudotime order,
//! and the connectivity between types. Each is an SVG that legume-plot
//! renders to PNG for the terminal or writes as a figure file; the TUI's
//! figure pane is where they are shown and exported.

use super::edges::{self, EdgeRow, Verdict};
use crate::manifest::run::{derive_out_prefix, resolve, RunManifest};
use anyhow::{Context, Result};
use enrichment::UNASSIGNED_LABEL;
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::parquet::{peek_parquet_field_names, read_table_columns};
use legume_numeric::matrix::traits::IoOps;
use legume_numeric::matrix::utils::{median, quantiles};
use legume_plot::hinton::{hinton_size, render_hinton, HintonOpts};
use legume_plot::palette::{self, sample_blue_red, Palette, Rgb};
use legume_plot::rasterize::{
    rasterize_arrow_layer_png, rasterize_per_point_png, DataBounds, Extent, PointShape,
};
use legume_plot::svg_emit::{emit_svg, escape_xml, SvgOpts, TopicLayer};
use legume_plot::FigureFormats;
use rustc_hash::FxHashMap;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};

/// Cells no root reaches, pairs outside the prior, and the frame.
const GREY: Rgb = (190, 190, 190);
const INK: Rgb = (40, 40, 40);
/// An inferred edge's arrow: the data's edge, not the prior's.
const INFERRED: Rgb = (115, 115, 115);
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
}

impl Colouring {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Pseudotime => "pseudotime",
            Self::Type => "cell type",
            Self::Lineage => "lineage",
            Self::Component => "component",
        }
    }

    /// Its slot in [`TrajectoryData`]'s colour cache.
    fn slot(self) -> usize {
        self as usize
    }
}

/// How the scatters are drawn: the label size (an index into
/// [`TEXT_SCALES`], `None` for no labels) and the colouring.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Style {
    pub(crate) labels: Option<usize>,
    pub(crate) colouring: Colouring,
}

impl Default for Style {
    fn default() -> Self {
        Self {
            labels: Some(1),
            colouring: Colouring::Pseudotime,
        }
    }
}

impl Style {
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

    /// The part of `whole` on screen.
    pub(crate) fn of(&self, whole: &DataBounds) -> DataBounds {
        let half = 0.5 / self.factor();
        let (w, h) = (whole.xmax - whole.xmin, whole.ymax - whole.ymin);
        DataBounds {
            xmin: whole.xmin + (self.cx - half) * w,
            xmax: whole.xmin + (self.cx + half) * w,
            ymin: whole.ymin + (self.cy - half) * h,
            ymax: whole.ymin + (self.cy + half) * h,
        }
    }
}

/// Which figure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Panel {
    /// One of the run's 2D layouts (an index into `TrajectoryData::layouts`),
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

/// A trajectory run's outputs, aligned to its cells.
#[derive(Default)]
pub(crate) struct TrajectoryData {
    /// The manifest the outputs were read from; export names derive from it.
    pub(crate) manifest: PathBuf,
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
    /// What to run for the layouts of [`WANTED_LAYOUTS`] it lacks, if any.
    pub(crate) layout_hint: Option<String>,
    pub(crate) edges: Vec<EdgeRow>,
    /// Each colouring's cell colours and legend, made on first use.
    pub(crate) colours: [OnceLock<CellColours>; 4],
    /// Each scatter's type medians over all its cells, made on first use.
    pub(crate) medians: Mutex<Vec<(Panel, Arc<Medians>)>>,
}

/// Each cell's colour (`None` for grey) and the legend of a categorical
/// colouring.
type CellColours = (Vec<Option<Rgb>>, Vec<(String, Rgb)>);
/// Each type's median position on a scatter.
type Medians = BTreeMap<Box<str>, (f32, f32)>;

/// A 2D layout of the run's cells, NaN where a cell has none.
pub(crate) struct Layout {
    /// Which of senna's layouts it is (`phate`, `umap`, …), when known.
    pub(crate) method: Option<String>,
    pub(crate) x: Vec<f32>,
    pub(crate) y: Vec<f32>,
}

impl TrajectoryData {
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
        let index: FxHashMap<&str, usize> = cells
            .iter()
            .enumerate()
            .map(|(i, c)| (c.as_ref(), i))
            .collect();

        let diffusion = match t.diffusion.as_deref() {
            Some(rel) => {
                let d = Mat::from_parquet_with_row_names(&at(rel), Some(0))?;
                Some(aligned(&d.rows, &d.mat, &index, cells.len()))
            }
            None => None,
        };
        let mut all = Vec::new();
        let (found, run) = layouts_along(manifest, dir, file);
        for (method, path) in &found {
            match read_layout(path, &index) {
                Ok((x, y)) => all.push(Layout {
                    method: method.clone(),
                    x,
                    y,
                }),
                Err(e) => log::warn!("layout {path}: {e:#}"),
            }
        }
        let methods: Vec<Option<String>> = all.iter().map(|l| l.method.clone()).collect();
        let layout_hint = layout_hint(&methods, &run);
        if let Some(h) = &layout_hint {
            log::warn!("{h}");
        }
        let edges = match t.edges.as_deref() {
            Some(rel) => edges::read(&at(rel))?,
            None => Vec::new(),
        };
        Ok(Some(Self {
            manifest: file.to_path_buf(),
            pseudotime,
            types,
            component,
            lineage,
            diffusion,
            layouts: all,
            layout_hint,
            edges,
            ..Self::default()
        }))
    }

    /// The run's output prefix, which export names start from.
    pub(crate) fn prefix(&self) -> String {
        derive_out_prefix(&self.manifest.to_string_lossy())
    }

    /// A panel's title, naming the layout drawn and, on a scatter, the
    /// colouring.
    pub(crate) fn title(&self, panel: Panel, colouring: Colouring) -> String {
        let colouring = colouring.name();
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
            out.push(first_spread_pair(d));
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
        let mut out = vec![Colouring::Pseudotime, Colouring::Type];
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
        out.push(Panel::Order);
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

    /// The figure, drawn to fit `w × h` pixels; a scatter in `style` at
    /// `view` (a thumbnail shows it whole).
    pub(crate) fn figure(
        &self,
        panel: Panel,
        w: u32,
        h: u32,
        style: &Style,
        view: View,
    ) -> Result<Figure> {
        match panel {
            Panel::Layout { k } => {
                let l = self.layouts.get(k).context("the run has no layout")?;
                self.scatter(panel, &l.x, &l.y, (w, h), true, style, view)
            }
            Panel::Diffusion { x, y } => {
                let d = self.diffusion.as_ref().context("no diffusion map")?;
                let cx: Vec<f32> = d.column(x).iter().copied().collect();
                let cy: Vec<f32> = d.column(y).iter().copied().collect();
                self.scatter(panel, &cx, &cy, (w, h), false, style, view)
            }
            Panel::Order => Ok(self.order_figure(w, h)),
            Panel::Connectivity => Ok(self.connectivity_figure(w, h)),
        }
    }

    /// Write `panel` as `{prefix}.trajectory.{panel}.svg` and `.pdf`, moved
    /// past existing files as `-2`, `-3` …, a scatter as it is on screen.
    pub(crate) fn export(&self, panel: Panel, style: &Style, view: View) -> Result<Export> {
        let base = free_base(&format!(
            "{}.trajectory.{}",
            self.prefix(),
            self.slug(panel)
        ));
        let w = (EXPORT_WIDTH_IN * EXPORT_DPI) as u32;
        let fig = self.figure(panel, w, w * 3 / 4, style, view)?;
        let b = base.to_string_lossy();
        let formats = FigureFormats {
            svg: true,
            png: false,
            pdf: true,
        };
        legume_plot::write_figure(&fig.svg, fig.w, fig.h, &b, formats)?;
        let files = vec![
            PathBuf::from(format!("{b}.svg")),
            PathBuf::from(format!("{b}.pdf")),
        ];
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
                self.title(panel, style.colouring),
            ),
        })
    }

    /// The cells of `(x, y)` inside `view` in the style's colouring (grey
    /// when they have no value), the type labels at the medians of their
    /// cells on screen, with the prior's direct edges as arrows between type
    /// medians when `arrows`.
    #[allow(clippy::too_many_arguments)]
    fn scatter(
        &self,
        panel: Panel,
        x: &[f32],
        y: &[f32],
        (w, h): (u32, u32),
        arrows: bool,
        style: &Style,
        view: View,
    ) -> Result<Figure> {
        let ext = Extent { w, h };
        let finite: Vec<usize> = (0..x.len())
            .filter(|&i| x[i].is_finite() && y[i].is_finite())
            .collect();
        anyhow::ensure!(!finite.is_empty(), "no cell has coordinates");
        let bounds = view.of(&bounds_of(x, y, &finite));
        let inside = |i: &usize| {
            (bounds.xmin..=bounds.xmax).contains(&x[*i])
                && (bounds.ymin..=bounds.ymax).contains(&y[*i])
        };
        let shown: Vec<usize> = finite.iter().copied().filter(inside).collect();
        let (colour_of, legend) = self.cell_colours(style.colouring);
        let order = self.draw_order(&shown, style.colouring);
        let pts: Vec<(f32, f32)> = order
            .iter()
            .map(|&i| to_pixel((x[i], y[i]), &bounds, ext))
            .collect();
        let colors: Vec<Rgb> = order
            .iter()
            .map(|&i| colour_of[i].unwrap_or(GREY))
            .collect();
        let radius = (w.min(h) as f32 / 400.0).clamp(1.0, 3.0);
        let mut layers = vec![raster_layer(rasterize_per_point_png(
            &pts,
            &colors,
            ext,
            radius,
            0.85,
            PointShape::Circle,
        )?)];
        // The medians over all cells anchor the arrows, and the labels of
        // the whole scatter; only made when one of them is drawn.
        let whole_labels = style.labels.is_some() && view.is_whole();
        let all = (arrows || whole_labels).then(|| self.medians(panel, x, y, &finite));
        if let (true, Some(all)) = (arrows, &all) {
            layers.extend(self.arrow_layers(all, &bounds, ext, radius)?);
        }
        let mut svg = emit_svg(
            &layers,
            &SvgOpts {
                width_px: w,
                height_px: h,
                frame_stroke_px: 1.0,
                background: Some((255, 255, 255)),
                ..SvgOpts::default()
            },
        );
        // Text over the layers, as legume-plot draws its labels. Spliced in:
        // legume-plot's labels come with an image per layer.
        let base_font = (w.min(h) as f32 / 48.0).max(6.0);
        let mut text = legend_svg(legend, base_font, h);
        if let Some(k) = style.labels {
            let font = base_font * TEXT_SCALES[k];
            let on_screen;
            let medians = match &all {
                Some(all) if view.is_whole() => &**all,
                _ => {
                    on_screen = self.type_medians(x, y, &shown);
                    &on_screen
                }
            };
            for (t, &at) in medians {
                if t.as_ref() != UNASSIGNED_LABEL {
                    let xy = keep_inside(to_pixel(at, &bounds, ext), font, t, ext);
                    text += &svg_text(xy, font, Some("middle"), true, t);
                }
            }
        }
        let end = svg.rfind("</svg>").context("an SVG without its end")?;
        svg.insert_str(end, &text);
        Ok(Figure { svg, w, h })
    }

    /// Each cell's colour in `colouring` and its legend, made once.
    /// The order `cells` are drawn in, later over earlier: grey cells first,
    /// so coloured ones draw over them; by pseudotime, low to high, so the
    /// few late cells are not buried under the many early ones.
    fn draw_order(&self, cells: &[usize], colouring: Colouring) -> Vec<usize> {
        let (colour_of, _) = self.cell_colours(colouring);
        let mut order = cells.to_vec();
        match colouring {
            Colouring::Pseudotime => {
                let t = |i: usize| self.pseudotime[i];
                order.sort_by(|&a, &b| {
                    (t(a).is_finite().cmp(&t(b).is_finite())).then(t(a).total_cmp(&t(b)))
                });
            }
            _ => order.sort_by_key(|&i| colour_of[i].is_some()),
        }
        order
    }

    fn cell_colours(&self, colouring: Colouring) -> &CellColours {
        self.colours[colouring.slot()].get_or_init(|| match colouring {
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
        })
    }

    /// The type medians over all of `cells` on the scatter `panel`, kept for
    /// the next figure of it. They are made outside the lock, so figures
    /// drawn in parallel do not wait on each other.
    fn medians(&self, panel: Panel, x: &[f32], y: &[f32], cells: &[usize]) -> Arc<Medians> {
        let kept = |m: &Vec<(Panel, Arc<Medians>)>| {
            m.iter().find(|(p, _)| *p == panel).map(|(_, m)| m.clone())
        };
        if let Some(m) = kept(&self.medians.lock().expect("not poisoned")) {
            return m;
        }
        let made = Arc::new(self.type_medians(x, y, cells));
        let mut all = self.medians.lock().expect("not poisoned");
        if let Some(m) = kept(&all) {
            return m;
        }
        all.push((panel, made.clone()));
        made
    }

    /// The prior's direct edges as arrows between the type `medians`:
    /// supported ones solid, the rest faded.
    fn arrow_layers(
        &self,
        medians: &Medians,
        bounds: &DataBounds,
        ext: Extent,
        radius: f32,
    ) -> Result<Vec<TopicLayer>> {
        let mut strong = Vec::new();
        let mut weak = Vec::new();
        let mut inferred = Vec::new();
        let drawn = |e: &&EdgeRow| e.in_prior || e.verdict == Some(Verdict::Inferred);
        for e in self.edges.iter().filter(drawn) {
            if let (Some(&a), Some(&b)) = (medians.get(&e.a), medians.get(&e.b)) {
                let seg = (to_pixel(a, bounds, ext), to_pixel(b, bounds, ext));
                match e.verdict {
                    Some(Verdict::Supported) => strong.push(seg),
                    Some(Verdict::Inferred) => inferred.push(seg),
                    _ => weak.push(seg),
                }
            }
        }
        let stroke = radius * 1.2;
        let mut layers = Vec::new();
        for (segs, colour, alpha) in [
            (strong, INK, 0.9),
            (weak, INK, 0.3),
            (inferred, INFERRED, 0.9),
        ] {
            if !segs.is_empty() {
                layers.push(raster_layer(rasterize_arrow_layer_png(
                    &segs,
                    ext,
                    stroke,
                    stroke * 5.0,
                    colour,
                    alpha,
                )?));
            }
        }
        Ok(layers)
    }

    /// Each type's median position among `cells`.
    fn type_medians(&self, x: &[f32], y: &[f32], cells: &[usize]) -> Medians {
        let mut by_type: BTreeMap<&str, (Vec<f32>, Vec<f32>)> = BTreeMap::new();
        for &i in cells {
            let e = by_type.entry(self.types[i].as_ref()).or_default();
            e.0.push(x[i]);
            e.1.push(y[i]);
        }
        by_type
            .into_iter()
            .map(|(t, (xs, ys))| (t.into(), (median(&xs), median(&ys))))
            .collect()
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

/// Data → pixel with y pointing up (larger data-y higher on screen).
fn to_pixel(p: (f32, f32), bounds: &DataBounds, ext: Extent) -> (f32, f32) {
    let (x, y) = bounds.to_pixel(p, ext);
    (x, ext.h as f32 - y)
}

/// The bounding box of the points `idx` of `(x, y)`.
fn bounds_of(x: &[f32], y: &[f32], idx: &[usize]) -> DataBounds {
    let (mut x0, mut x1, mut y0, mut y1) = (f32::MAX, f32::MIN, f32::MAX, f32::MIN);
    for &i in idx {
        x0 = x0.min(x[i]);
        x1 = x1.max(x[i]);
        y0 = y0.min(y[i]);
        y1 = y1.max(y[i]);
    }
    DataBounds::from_minmax(x0, x1, y0, y1)
}

/// The layouts to draw pseudotime on, best first: senna's PHATE (`senna
/// layout phate`, recorded under `layout.methods.phate`), as PHATE is built
/// to show trajectories, then the run's current layout, then the other
/// layouts it records. Each with its method's name, when known, and its
/// manifest-relative path.
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
    let mut source = manifest.annotate.source.as_deref().map(|s| resolve(dir, s));
    // A chain of rounds is short; the bound stops a cycle.
    for _ in 0..32 {
        let Some(src) = source.take() else { break };
        let Ok(l) = crate::manifest::run::load(&src) else {
            break;
        };
        add(&l.manifest, &l.dir);
        source = l
            .manifest
            .annotate
            .source
            .as_deref()
            .map(|s| resolve(&l.dir, s));
        run = l.file;
    }
    out.sort_by_key(|(m, _)| m.as_deref() != Some("phate"));
    (out, run)
}

/// What to run for the [`WANTED_LAYOUTS`] missing from `found` (the
/// layouts' methods), on `run`.
pub(crate) fn layout_hint(found: &[Option<String>], run: &Path) -> Option<String> {
    let missing: Vec<&str> = WANTED_LAYOUTS
        .into_iter()
        .filter(|w| !found.iter().any(|m| m.as_deref() == Some(*w)))
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
        "no {} for this run: `senna layout {which} --from {name}` adds it, and the figures pick it up when they reopen",
        shown.join(" or ")
    ))
}

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

/// The layout at `path` as `(x, y)` per cell of `index`, NaN for a cell the
/// layout lacks; an error when it has none of them, so the next layout is
/// tried.
fn read_layout(path: &str, index: &FxHashMap<&str, usize>) -> Result<(Vec<f32>, Vec<f32>)> {
    let t = Mat::from_parquet(path)?;
    let col = |name: &str| t.cols.iter().position(|c| c.as_ref() == name);
    let (Some(x), Some(y)) = (col("x"), col("y")) else {
        anyhow::bail!("{path} has no x and y columns");
    };
    if !t.rows.iter().any(|r| index.contains_key(r.as_ref())) {
        anyhow::bail!("{path} shares no cell with the trajectory");
    }
    let m = aligned(&t.rows, &t.mat, index, index.len());
    Ok((
        m.column(x).iter().copied().collect(),
        m.column(y).iter().copied().collect(),
    ))
}

/// `mat`'s rows reordered to `cells` by name; a cell the matrix lacks is NaN.
fn aligned(rows: &[Box<str>], mat: &Mat, index: &FxHashMap<&str, usize>, n: usize) -> Mat {
    let mut out = Mat::from_element(n, mat.ncols(), f32::NAN);
    for (i, r) in rows.iter().enumerate() {
        if let Some(&j) = index.get(r.as_ref()) {
            out.set_row(j, &mat.row(i));
        }
    }
    out
}

/// A raster-only layer: no hull, no label.
fn raster_layer(png: Vec<u8>) -> TopicLayer {
    TopicLayer {
        label: String::new(),
        png,
        hull_px: Vec::new(),
        label_xy_px: (f32::NAN, f32::NAN),
        color: INK,
    }
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

/// `n` different colours: legume-plot's palette for `n` while it has that
/// many, then hues a golden angle apart in two lightnesses, so no two
/// groups share a colour however many there are.
fn distinct_colours(n: usize) -> Vec<Rgb> {
    let pal = palette::resolve(&Palette::Auto, n);
    // The palette cycles; its own colours are those before the first repeat.
    let own = (1..n)
        .find(|&j| palette::color(&pal, j) == palette::color(&pal, 0))
        .unwrap_or(n);
    (0..n)
        .map(|j| {
            if j < own {
                return palette::color(&pal, j);
            }
            let k = j - own;
            let hue = (k as f32 * 137.508) % 360.0;
            let light = if k % 2 == 0 { 0.42 } else { 0.62 };
            hsl(hue, 0.65, light)
        })
        .collect()
}

/// An HSL colour (hue in degrees, saturation and lightness in 0..=1).
fn hsl(h: f32, s: f32, l: f32) -> Rgb {
    let c = (1.0 - (2.0 * l - 1.0).abs()) * s;
    let x = c * (1.0 - ((h / 60.0) % 2.0 - 1.0).abs());
    let (r, g, b) = match (h / 60.0) as u32 {
        0 => (c, x, 0.0),
        1 => (x, c, 0.0),
        2 => (0.0, c, x),
        3 => (0.0, x, c),
        4 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };
    let m = l - c / 2.0;
    let byte = |v: f32| ((v + m) * 255.0).round().clamp(0.0, 255.0) as u8;
    (byte(r), byte(g), byte(b))
}

/// Lineage or component codes as keys, -1 (none) as `None`.
fn codes(v: &[i32]) -> Vec<Option<i32>> {
    v.iter().map(|&c| (c >= 0).then_some(c)).collect()
}

/// `text` at `xy` in [`INK`], centred vertically and anchored at `anchor`
/// (the start when `None`), with a white halo as legume-plot draws a label
/// when `halo`.
/// `xy` moved so a centred label of `text` stays inside the figure, its
/// width taken as about 0.55 font sizes per character.
fn keep_inside(xy: (f32, f32), font: f32, text: &str, ext: Extent) -> (f32, f32) {
    let half = 0.275 * font * text.chars().count() as f32 + 2.0;
    let (w, h) = (ext.w as f32, ext.h as f32);
    let x = if 2.0 * half >= w {
        w / 2.0
    } else {
        xy.0.clamp(half, w - half)
    };
    (x, xy.1.clamp(font, (h - font).max(font)))
}

fn svg_text(xy: (f32, f32), font: f32, anchor: Option<&str>, halo: bool, text: &str) -> String {
    let (x, y) = xy;
    let (r, g, b) = INK;
    let anchor = anchor.map_or(String::new(), |a| format!(" text-anchor='{a}'"));
    let halo = if halo {
        format!(
            " paint-order='stroke' stroke='white' stroke-width='{:.2}' stroke-linejoin='round'",
            (font * 0.35).max(1.5)
        )
    } else {
        String::new()
    };
    format!(
        "<text x='{x:.2}' y='{y:.2}' font-family='Helvetica, Arial, sans-serif' font-size='{font:.2}'{anchor} \
         dominant-baseline='central'{halo} fill='rgb({r},{g},{b})'>{}</text>\n",
        escape_xml(text)
    )
}

/// A categorical colouring's legend in the top-left corner of a figure
/// `height` pixels high: a swatch and a name per entry, at most
/// [`LEGEND_MAX`] of them and no more than fill half the figure.
fn legend_svg(entries: &[(String, Rgb)], font: f32, height: u32) -> String {
    let mut s = String::new();
    let row = font * 1.4;
    // At most half the figure's height, counting the "+k more" line; none
    // in a figure too small for two rows (a grid thumbnail).
    let fits = ((height as f32 * 0.5 - font) / row).floor().max(0.0) as usize;
    if fits < 2 {
        return s;
    }
    let shown = if entries.len() <= fits.min(LEGEND_MAX) {
        entries.len()
    } else {
        fits.min(LEGEND_MAX) - 1
    };
    for (k, (name, (r, g, b))) in entries.iter().take(shown).enumerate() {
        let y = font + row * k as f32;
        s += &format!(
            "<rect x='{:.2}' y='{:.2}' width='{font:.2}' height='{font:.2}' fill='rgb({r},{g},{b})'/>",
            font * 0.6,
            y - font * 0.5,
        );
        s += &svg_text((font * 1.9, y), font, None, true, name);
    }
    if entries.len() > shown {
        let more = format!("+{} more", entries.len() - shown);
        s += &svg_text(
            (font * 0.6, font + row * shown as f32),
            font,
            None,
            false,
            &more,
        );
    }
    s
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

/// `fig` rendered to pixels, through legume-plot's renderer (which writes a
/// file, here a private temporary one) and the `image` crate.
pub(crate) fn render(fig: &Figure) -> Result<image::RgbaImage> {
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

/// The first two diffusion components after the trivial one that are spread
/// over the cells ([`DC_MIN_SHARE`]): the components of a barely attached
/// group show that group against a line of every other cell. DC1 × DC2 when
/// fewer than two are.
fn first_spread_pair(d: &Mat) -> Panel {
    use super::diffusion::{participation, DC_MIN_SHARE};
    let spread: Vec<usize> = (1..d.ncols())
        .filter(|&j| {
            participation(d.column(j).iter().map(|&v| f64::from(v))) >= f64::from(DC_MIN_SHARE)
        })
        .take(2)
        .collect();
    match spread[..] {
        [x, y] => Panel::Diffusion { x, y },
        _ => Panel::Diffusion { x: 1, y: 2 },
    }
}

#[cfg(test)]
#[path = "figures_tests.rs"]
mod tests;
