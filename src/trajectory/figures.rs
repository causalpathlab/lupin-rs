//! Figures of a trajectory run, from the outputs its manifest records
//! (`docs/trajectory-plan.md` §6): the layout coloured by pseudotime with the
//! prior's edges as arrows, the diffusion map, the types in pseudotime order,
//! and the connectivity between types. Each is an SVG that legume-plot
//! renders to PNG for the terminal or writes as a figure file; the TUI's
//! figure pane is where they are shown and exported.

use super::edges::{self, EdgeRow, Verdict};
use crate::manifest::run::{derive_out_prefix, resolve, RunManifest};
use anyhow::{Context, Result};
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::parquet::read_table_columns;
use legume_numeric::matrix::traits::IoOps;
use legume_numeric::matrix::utils::{median, quantiles};
use legume_plot::hinton::{hinton_size, render_hinton, HintonOpts};
use legume_plot::palette::{sample_blue_red, Rgb};
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
/// Figure width in inches and dots per inch for an export.
const EXPORT_WIDTH_IN: f32 = 7.0;
const EXPORT_DPI: f32 = 200.0;

/// Which figure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Panel {
    /// The run's 2D layout (senna's PHATE when it made one), cells coloured
    /// by pseudotime, prior edges as arrows.
    Layout,
    /// Two diffusion components, cells coloured by pseudotime.
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
            Self::Layout => "layout".into(),
            Self::Diffusion { x, y } => format!("diffusion_dc{x}_dc{y}"),
            Self::Order => "order".into(),
            Self::Connectivity => "connectivity".into(),
        }
    }

    pub(crate) fn title(self) -> String {
        match self {
            Self::Layout => "layout · pseudotime".into(),
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
    /// Cells × diffusion components.
    pub(crate) diffusion: Option<Mat>,
    /// The run's layout, `(x, y)` per cell, NaN where a cell has none.
    pub(crate) layout: Option<(Vec<f32>, Vec<f32>)>,
    /// Which of senna's layouts it is (`phate`, `umap`, …), when known.
    pub(crate) layout_method: Option<String>,
    pub(crate) edges: Vec<EdgeRow>,
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
        let (mut layout_method, mut layout) = (None, None);
        for (method, rel) in layouts(manifest) {
            match read_layout(&at(&rel), &index) {
                Ok(l) => {
                    (layout_method, layout) = (method, Some(l));
                    break;
                }
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
            diffusion,
            layout,
            layout_method,
            edges,
        }))
    }

    /// The run's output prefix, which export names start from.
    pub(crate) fn prefix(&self) -> String {
        derive_out_prefix(&self.manifest.to_string_lossy())
    }

    /// A panel's title, naming the layout drawn.
    pub(crate) fn title(&self, panel: Panel) -> String {
        match (panel, &self.layout_method) {
            (Panel::Layout, Some(m)) => format!("{} · pseudotime", m.to_uppercase()),
            _ => panel.title(),
        }
    }

    /// The figures these outputs allow, in display order.
    pub(crate) fn panels(&self) -> Vec<Panel> {
        let mut out = Vec::new();
        if self.layout.is_some() {
            out.push(Panel::Layout);
        }
        if self.diffusion.as_ref().is_some_and(|d| d.ncols() >= 3) {
            out.push(Panel::Diffusion { x: 1, y: 2 });
        }
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

    /// The figure, drawn to fit `w × h` pixels.
    pub(crate) fn figure(&self, panel: Panel, w: u32, h: u32) -> Result<Figure> {
        match panel {
            Panel::Layout => {
                let (x, y) = self.layout.as_ref().context("the run has no layout")?;
                self.scatter(x, y, w, h, true)
            }
            Panel::Diffusion { x, y } => {
                let d = self.diffusion.as_ref().context("no diffusion map")?;
                let cx: Vec<f32> = d.column(x).iter().copied().collect();
                let cy: Vec<f32> = d.column(y).iter().copied().collect();
                self.scatter(&cx, &cy, w, h, false)
            }
            Panel::Order => Ok(self.order_figure(w, h)),
            Panel::Connectivity => Ok(self.connectivity_figure(w, h)),
        }
    }

    /// Write `panel` as `{prefix}.trajectory.{panel}.svg` and `.pdf`, moved
    /// past existing files as `-2`, `-3` ….
    pub(crate) fn export(&self, panel: Panel) -> Result<Export> {
        let base = free_base(&format!("{}.trajectory.{}", self.prefix(), panel.slug()));
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
                "{} · {EXPORT_WIDTH_IN:.0} in · {EXPORT_DPI:.0} dpi",
                self.title(panel)
            ),
        })
    }

    /// Cells at `(x, y)` coloured by pseudotime (grey when none), with the
    /// prior's direct edges as arrows between type medians when `arrows`.
    fn scatter(&self, x: &[f32], y: &[f32], w: u32, h: u32, arrows: bool) -> Result<Figure> {
        let ext = Extent { w, h };
        let finite: Vec<usize> = (0..x.len())
            .filter(|&i| x[i].is_finite() && y[i].is_finite())
            .collect();
        anyhow::ensure!(!finite.is_empty(), "no cell has coordinates");
        let bounds = bounds_of(x, y, &finite);
        // Unreached cells first, so coloured ones draw over them.
        let mut order = finite.clone();
        order.sort_by_key(|&i| self.pseudotime[i].is_finite());
        let pts: Vec<(f32, f32)> = order
            .iter()
            .map(|&i| to_pixel((x[i], y[i]), &bounds, ext))
            .collect();
        let colors: Vec<Rgb> = order
            .iter()
            .map(|&i| {
                let t = self.pseudotime[i];
                if t.is_finite() {
                    sample_blue_red(t)
                } else {
                    GREY
                }
            })
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
        let svg = emit_svg(
            &layers,
            &SvgOpts {
                width_px: w,
                height_px: h,
                frame_stroke_px: 1.0,
                background: Some((255, 255, 255)),
                ..SvgOpts::default()
            },
        );
        Ok(Figure { svg, w, h })
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
/// to show trajectories, then the run's current layout. Each with its
/// method's name, when known, and its manifest-relative path.
fn layouts(manifest: &RunManifest) -> Vec<(Option<String>, String)> {
    let l = &manifest.layout;
    let phate = l
        .extra
        .get("methods")
        .and_then(|m| m.get("phate")?.get("cell_coords")?.as_str());
    let current = l.extra.get("current").and_then(|v| v.as_str());
    let mut out: Vec<(Option<String>, String)> = Vec::new();
    if let Some(p) = phate {
        out.push((Some("phate".into()), p.into()));
    }
    if let Some(p) = l.cell_coords.as_deref().filter(|p| Some(*p) != phate) {
        out.push((current.map(str::to_string), p.into()));
    }
    out
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
