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
use legume_numeric::matrix::parquet::{read_table_columns, ParquetReader};
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

/// Cells no root reaches, pairs outside the prior, and the frame.
const GREY: Rgb = (190, 190, 190);
const INK: Rgb = (40, 40, 40);
/// Entries a categorical legend lists before `+k more`.
const LEGEND_MAX: usize = 20;
/// Figure width in inches and dots per inch for an export.
const EXPORT_WIDTH_IN: f32 = 7.0;
const EXPORT_DPI: f32 = 200.0;

/// Label sizes `t` steps through, as `senna view` does: small, medium,
/// large, largest, then off.
pub(crate) const TEXT_SCALES: [f32; 4] = [1.0, 1.4, 1.8, 2.4];
const TEXT_SIZES: [&str; 4] = ["small", "medium", "large", "largest"];

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

/// How far `+` / `-` zoom a scatter, how far in it can go, and how far an
/// arrow pans, as a share of the part on screen.
const ZOOM_STEP: f32 = 1.4;
const MAX_ZOOM: f32 = 64.0;
const PAN_STEP: f32 = 0.2;

/// The part of a scatter on screen: the whole extent zoomed `zoom` times
/// about the centre `(cx, cy)`, both in shares of the whole extent (y up).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct View {
    pub(crate) zoom: f32,
    pub(crate) cx: f32,
    pub(crate) cy: f32,
}

impl Default for View {
    fn default() -> Self {
        Self {
            zoom: 1.0,
            cx: 0.5,
            cy: 0.5,
        }
    }
}

impl View {
    pub(crate) fn is_whole(&self) -> bool {
        self.zoom <= 1.0
    }

    /// Zoom in (`inward`) or out by one step about the centre.
    pub(crate) fn zoom(&mut self, inward: bool) {
        let f = if inward { ZOOM_STEP } else { 1.0 / ZOOM_STEP };
        self.zoom = (self.zoom * f).clamp(1.0, MAX_ZOOM);
        self.clamp();
    }

    /// Move the part on screen by `(dx, dy)` steps (right and up positive).
    pub(crate) fn pan(&mut self, dx: f32, dy: f32) {
        let span = 1.0 / self.zoom;
        self.cx += dx * PAN_STEP * span;
        self.cy += dy * PAN_STEP * span;
        self.clamp();
    }

    /// Keep the part on screen inside the whole extent.
    fn clamp(&mut self) {
        let half = 0.5 / self.zoom;
        self.cx = self.cx.clamp(half, 1.0 - half);
        self.cy = self.cy.clamp(half, 1.0 - half);
    }

    /// The part of `whole` on screen.
    pub(crate) fn of(&self, whole: &DataBounds) -> DataBounds {
        let half = 0.5 / self.zoom;
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
    pub(crate) edges: Vec<EdgeRow>,
    pub(crate) style: Style,
    /// The part of the scatter on screen (`+` `-` and the arrows).
    pub(crate) view: View,
}

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
        let (strings, numbers) = read_table_columns(&pt_path, &["cell", "type"], &["pseudotime"])
            .with_context(|| format!("reading {pt_path}"))?;
        let cells: Vec<Box<str>> = strings[0].clone();
        let types: Vec<Box<str>> = strings[1].clone();
        let pseudotime: Vec<f32> = numbers[0].iter().map(|&v| v as f32).collect();
        let (component, lineage) = components_and_lineages(&pt_path, cells.len())?;
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
        for (method, rel) in layouts(manifest) {
            match read_layout(&at(&rel), &index) {
                Ok((x, y)) => all.push(Layout { method, x, y }),
                Err(e) => log::warn!("layout {rel}: {e:#}"),
            }
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
            edges,
            style: Style::default(),
            view: View::default(),
        }))
    }

    /// The run's output prefix, which export names start from.
    pub(crate) fn prefix(&self) -> String {
        derive_out_prefix(&self.manifest.to_string_lossy())
    }

    /// A panel's title, naming the layout drawn and the colouring.
    pub(crate) fn title(&self, panel: Panel) -> String {
        let colouring = self.style.colouring.name();
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
    /// them: the layouts, then the diffusion map.
    pub(crate) fn scatters(&self) -> Vec<Panel> {
        let mut out: Vec<Panel> = (0..self.layouts.len())
            .map(|k| Panel::Layout { k })
            .collect();
        if self.diffusion.as_ref().is_some_and(|d| d.ncols() >= 3) {
            out.push(Panel::Diffusion { x: 1, y: 2 });
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

    /// Labels in turn: small, medium, large, largest, off; the note for the
    /// status line.
    pub(crate) fn cycle_labels(&mut self) -> String {
        let next = match self.style.labels {
            Some(k) if k + 1 < TEXT_SCALES.len() => Some(k + 1),
            Some(_) => None,
            None => Some(0),
        };
        self.style.labels = next;
        match next {
            Some(k) => {
                let then = TEXT_SIZES.get(k + 1).copied().unwrap_or("off");
                format!("labels {} · t for {then}", TEXT_SIZES[k])
            }
            None => "labels off · t shows them small".into(),
        }
    }

    /// The next colouring; the note for the status line.
    pub(crate) fn cycle_colouring(&mut self) -> String {
        let all = self.colourings();
        let at = all
            .iter()
            .position(|&c| c == self.style.colouring)
            .unwrap_or(0);
        self.style.colouring = all[(at + 1) % all.len()];
        let then = all[(at + 2) % all.len()];
        format!(
            "coloured by {} · c for {}",
            self.style.colouring.name(),
            then.name()
        )
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

    /// The next diffusion pair after `(x, y)`, cycling through the components
    /// after the trivial first one.
    pub(crate) fn next_pair(&self, x: usize, y: usize, forward: bool) -> (usize, usize) {
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
        let ny = step(y);
        if ny == x {
            (x, step(ny))
        } else {
            (x, ny)
        }
    }

    /// The figure, drawn to fit `w × h` pixels, a scatter at the view on
    /// screen.
    pub(crate) fn figure(&self, panel: Panel, w: u32, h: u32) -> Result<Figure> {
        self.figure_at(panel, w, h, self.view)
    }

    /// The figure with a scatter at `view` (a thumbnail shows it whole).
    pub(crate) fn figure_at(&self, panel: Panel, w: u32, h: u32, view: View) -> Result<Figure> {
        match panel {
            Panel::Layout { k } => {
                let l = self.layouts.get(k).context("the run has no layout")?;
                self.scatter(&l.x, &l.y, w, h, true, view)
            }
            Panel::Diffusion { x, y } => {
                let d = self.diffusion.as_ref().context("no diffusion map")?;
                let cx: Vec<f32> = d.column(x).iter().copied().collect();
                let cy: Vec<f32> = d.column(y).iter().copied().collect();
                self.scatter(&cx, &cy, w, h, false, view)
            }
            Panel::Order => Ok(self.order_figure(w, h)),
            Panel::Connectivity => Ok(self.connectivity_figure(w, h)),
        }
    }

    /// Write `panel` as `{prefix}.trajectory.{panel}.svg` and `.pdf`, moved
    /// past existing files as `-2`, `-3` ….
    pub(crate) fn export(&self, panel: Panel) -> Result<Export> {
        let base = free_base(&format!(
            "{}.trajectory.{}",
            self.prefix(),
            self.slug(panel)
        ));
        let w = (EXPORT_WIDTH_IN * EXPORT_DPI) as u32;
        let fig = self.figure(panel, w, w * 3 / 4)?;
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
        Ok(Export {
            base,
            files,
            what: format!(
                "{}{} · {EXPORT_WIDTH_IN:.0} in · {EXPORT_DPI:.0} dpi",
                self.title(panel),
                match panel {
                    Panel::Layout { .. } | Panel::Diffusion { .. } if !self.view.is_whole() =>
                        format!(" · zoomed ×{:.1}", self.view.zoom),
                    _ => String::new(),
                }
            ),
        })
    }

    /// The cells of `(x, y)` inside `view` in the style's colouring (grey
    /// when they have no value), the type labels at the medians of their
    /// cells on screen, with the prior's direct edges as arrows between type
    /// medians when `arrows`.
    fn scatter(
        &self,
        x: &[f32],
        y: &[f32],
        w: u32,
        h: u32,
        arrows: bool,
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
        let (colour_of, legend) = self.cell_colours();
        // Grey cells first, so coloured ones draw over them.
        let mut order = shown.clone();
        order.sort_by_key(|&i| colour_of[i].is_some());
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
        if arrows {
            layers.extend(self.arrow_layers(x, y, &finite, &bounds, ext, radius)?);
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
        // Text over the layers, as legume-plot draws its labels.
        let base_font = (w.min(h) as f32 / 48.0).max(6.0);
        let mut text = legend_svg(&legend, base_font);
        if let Some(k) = self.style.labels {
            let font = base_font * TEXT_SCALES[k];
            for (t, at) in self.type_medians(x, y, &shown) {
                if t != UNASSIGNED_LABEL {
                    text += &halo_text(to_pixel(at, &bounds, ext), font, INK, t);
                }
            }
        }
        let end = svg.rfind("</svg>").context("an SVG without its end")?;
        svg.insert_str(end, &text);
        Ok(Figure { svg, w, h })
    }

    /// Each cell's colour in the style's colouring (`None` for grey: no
    /// pseudotime, unassigned, or no single lineage or component), and the
    /// legend of a categorical colouring.
    fn cell_colours(&self) -> (Vec<Option<Rgb>>, Vec<(String, Rgb)>) {
        let categorical = |codes: &[i32], name: &dyn Fn(i32) -> String| {
            let mut used: Vec<i32> = codes.iter().copied().filter(|&c| c >= 0).collect();
            used.sort_unstable();
            used.dedup();
            let pal = palette::resolve(&Palette::Auto, used.len());
            let colour = |c: i32| used.binary_search(&c).ok().map(|k| palette::color(&pal, k));
            let legend = used
                .iter()
                .map(|&c| (name(c), colour(c).expect("used")))
                .collect();
            (codes.iter().map(|&c| colour(c)).collect(), legend)
        };
        match self.style.colouring {
            Colouring::Pseudotime => (
                self.pseudotime
                    .iter()
                    .map(|&t| t.is_finite().then(|| sample_blue_red(t)))
                    .collect(),
                Vec::new(),
            ),
            Colouring::Type => {
                let mut names: Vec<&str> = self
                    .types
                    .iter()
                    .map(AsRef::as_ref)
                    .filter(|&t| t != UNASSIGNED_LABEL)
                    .collect();
                names.sort_unstable();
                names.dedup();
                let codes: Vec<i32> = self
                    .types
                    .iter()
                    .map(|t| names.binary_search(&t.as_ref()).map_or(-1, |k| k as i32))
                    .collect();
                categorical(&codes, &|c| names[c as usize].to_string())
            }
            Colouring::Lineage => categorical(&self.lineage, &|c| format!("lineage {c}")),
            Colouring::Component => categorical(&self.component, &|c| format!("component {c}")),
        }
    }

    /// The prior's direct edges as arrows between type medians: supported
    /// ones solid, the rest faded.
    fn arrow_layers(
        &self,
        x: &[f32],
        y: &[f32],
        cells: &[usize],
        bounds: &DataBounds,
        ext: Extent,
        radius: f32,
    ) -> Result<Vec<TopicLayer>> {
        let medians = self.type_medians(x, y, cells);
        let mut strong = Vec::new();
        let mut weak = Vec::new();
        for e in self.edges.iter().filter(|e| e.in_prior) {
            if let (Some(&a), Some(&b)) = (medians.get(e.a.as_ref()), medians.get(e.b.as_ref())) {
                let seg = (to_pixel(a, bounds, ext), to_pixel(b, bounds, ext));
                if e.verdict == Some(Verdict::Supported) {
                    strong.push(seg);
                } else {
                    weak.push(seg);
                }
            }
        }
        let stroke = radius * 1.2;
        let mut layers = Vec::new();
        for (segs, alpha) in [(strong, 0.9), (weak, 0.3)] {
            if !segs.is_empty() {
                layers.push(raster_layer(rasterize_arrow_layer_png(
                    &segs,
                    ext,
                    stroke,
                    stroke * 5.0,
                    INK,
                    alpha,
                )?));
            }
        }
        Ok(layers)
    }

    /// Each type's median position among `cells`.
    fn type_medians(&self, x: &[f32], y: &[f32], cells: &[usize]) -> BTreeMap<&str, (f32, f32)> {
        let mut by_type: BTreeMap<&str, (Vec<f32>, Vec<f32>)> = BTreeMap::new();
        for &i in cells {
            let e = by_type.entry(self.types[i].as_ref()).or_default();
            e.0.push(x[i]);
            e.1.push(y[i]);
        }
        by_type
            .into_iter()
            .map(|(t, (xs, ys))| (t, (median(&xs), median(&ys))))
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

/// Per cell of the pseudotime table at `path` (`n` rows): its component
/// and the one lineage its type lies on (-1 when shared or none), from the
/// `component` and `L0`, `L1`, … columns; all -1 when they are absent.
fn components_and_lineages(path: &str, n: usize) -> Result<(Vec<i32>, Vec<i32>)> {
    let t =
        ParquetReader::new(path, None, None, None).with_context(|| format!("reading {path}"))?;
    let ncols = t.column_names.len();
    let col = |name: &str| t.column_names.iter().position(|c| c.as_ref() == name);
    let at = |i: usize, j: usize| t.row_major_data[i * ncols + j];
    let rows = t.row_major_data.len().checked_div(ncols).unwrap_or(0);
    anyhow::ensure!(rows == n, "{path}: {rows} rows for {n} cells");
    let component = match col("component") {
        Some(j) => (0..n).map(|i| at(i, j) as i32).collect(),
        None => vec![-1; n],
    };
    let lineages: Vec<usize> = (0..).map_while(|k| col(&format!("L{k}"))).collect();
    let lineage = (0..n)
        .map(|i| {
            let w: Vec<f64> = lineages.iter().map(|&j| at(i, j)).collect();
            let best = w.iter().copied().fold(0.0, f64::max);
            let on: Vec<usize> = (0..w.len()).filter(|&k| w[k] == best).collect();
            match on.as_slice() {
                [k] if best > 0.0 => *k as i32,
                _ => -1,
            }
        })
        .collect();
    Ok((component, lineage))
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

/// `text` centred on `xy` with a white halo, as legume-plot draws a label.
fn halo_text(xy: (f32, f32), font: f32, colour: Rgb, text: &str) -> String {
    let (x, y) = xy;
    let (r, g, b) = colour;
    format!(
        "<text x='{x:.2}' y='{y:.2}' font-family='Helvetica, Arial, sans-serif' font-size='{font:.2}' \
         text-anchor='middle' dominant-baseline='central' paint-order='stroke' stroke='white' \
         stroke-width='{:.2}' stroke-linejoin='round' fill='rgb({r},{g},{b})'>{}</text>\n",
        (font * 0.35).max(1.5),
        escape_xml(text)
    )
}

/// A categorical colouring's legend in the top-left corner: a swatch and a
/// name per entry, at most [`LEGEND_MAX`] of them.
fn legend_svg(entries: &[(String, Rgb)], font: f32) -> String {
    let mut s = String::new();
    let row = font * 1.4;
    let shown = entries.len().min(LEGEND_MAX);
    for (k, (name, (r, g, b))) in entries.iter().take(shown).enumerate() {
        let y = font + row * k as f32;
        s += &format!(
            "<rect x='{:.2}' y='{:.2}' width='{font:.2}' height='{font:.2}' fill='rgb({r},{g},{b})'/>\
             <text x='{:.2}' y='{:.2}' font-family='Helvetica, Arial, sans-serif' font-size='{font:.2}' \
             dominant-baseline='central' paint-order='stroke' stroke='white' stroke-width='{:.2}' \
             fill='rgb(40,40,40)'>{}</text>\n",
            font * 0.6,
            y - font * 0.5,
            font * 1.9,
            y,
            (font * 0.35).max(1.5),
            escape_xml(name)
        );
    }
    if entries.len() > shown {
        s += &format!(
            "<text x='{:.2}' y='{:.2}' font-family='Helvetica, Arial, sans-serif' font-size='{font:.2}' \
             dominant-baseline='central' fill='rgb(40,40,40)'>+{} more</text>\n",
            font * 0.6,
            font + row * shown as f32,
            entries.len() - shown
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
    let tmp = tempfile::Builder::new()
        .prefix("lupin-figure-")
        .suffix(".png")
        .tempfile()
        .context("creating a temporary file for the figure")?;
    legume_plot::render_png(&fig.svg, fig.w, fig.h, tmp.path())?;
    let img = image::open(tmp.path())
        .with_context(|| format!("reading {}", tmp.path().display()))?
        .to_rgba8();
    Ok(img)
}

#[cfg(test)]
#[path = "figures_tests.rs"]
mod tests;
