//! Figures of a trajectory run, from the outputs its manifest records
//! (`docs/trajectory-plan.md` §6): the layout coloured by pseudotime with the
//! prior's edges as arrows, the diffusion map, the types in pseudotime order,
//! and the connectivity between types. Each is an SVG that legume-plot
//! renders to PNG for the terminal or writes as a figure file.

use crate::manifest::run::{derive_out_prefix, resolve, Loaded};
use anyhow::{Context, Result};
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::parquet::read_table_columns;
use legume_numeric::matrix::traits::IoOps;
use legume_numeric::matrix::utils::median;
use legume_plot::palette::{sample_blue_red, Rgb};
use legume_plot::rasterize::{
    rasterize_arrow_layer_png, rasterize_per_point_png, DataBounds, Extent, PointShape,
};
use legume_plot::svg_emit::{emit_svg, escape_xml, SvgOpts, TopicLayer};
use rustc_hash::FxHashMap;
use std::collections::BTreeMap;
use std::path::Path;

/// Cells no root reaches, and the frame.
const GREY: Rgb = (190, 190, 190);
const INK: Rgb = (40, 40, 40);

/// Which figure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Panel {
    /// The run's 2D layout, cells coloured by pseudotime, prior edges as arrows.
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

/// A type pair from `{out}.trajectory_edges.parquet`.
pub(crate) struct EdgeRow {
    pub(crate) a: Box<str>,
    pub(crate) b: Box<str>,
    pub(crate) connectivity: f32,
    pub(crate) in_prior: bool,
    pub(crate) verdict: Box<str>,
}

/// A trajectory run's outputs, aligned to its cells.
pub(crate) struct TrajectoryData {
    /// The manifest the outputs were read from; export names derive from it.
    pub(crate) prefix: String,
    /// Per cell; NaN for cells no root reaches.
    pub(crate) pseudotime: Vec<f32>,
    pub(crate) types: Vec<Box<str>>,
    /// Cells × diffusion components.
    pub(crate) diffusion: Option<Mat>,
    /// The run's layout, `(x, y)` per cell, NaN where a cell has none.
    pub(crate) layout: Option<(Vec<f32>, Vec<f32>)>,
    pub(crate) edges: Vec<EdgeRow>,
}

impl TrajectoryData {
    /// The outputs the manifest `loaded` records; `None` when it has no
    /// `trajectory.pseudotime`.
    pub(crate) fn load(loaded: &Loaded) -> Result<Option<Self>> {
        let t = &loaded.manifest.trajectory;
        let Some(pt_rel) = t.pseudotime.as_deref() else {
            return Ok(None);
        };
        let at = |rel: &str| resolve(&loaded.dir, rel);
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
        let layout = match loaded.manifest.layout.cell_coords.as_deref() {
            Some(rel) => {
                let path = at(rel);
                let (names, cols) = crate::plot::scatter::read_cell_coords(&path)?;
                let (Some(x), Some(y)) = (cols.get("x"), cols.get("y")) else {
                    anyhow::bail!("{path} has no x and y columns");
                };
                let mut lx = vec![f32::NAN; cells.len()];
                let mut ly = vec![f32::NAN; cells.len()];
                for (i, name) in names.iter().enumerate() {
                    if let Some(&j) = index.get(name.as_ref()) {
                        lx[j] = x[i];
                        ly[j] = y[i];
                    }
                }
                Some((lx, ly))
            }
            None => None,
        };
        let edges = match t.edges.as_deref() {
            Some(rel) => {
                let (s, n) = read_table_columns(
                    &at(rel),
                    &["a", "b", "in_prior", "verdict"],
                    &["connectivity"],
                )?;
                (0..s[0].len())
                    .map(|i| EdgeRow {
                        a: s[0][i].clone(),
                        b: s[1][i].clone(),
                        connectivity: n[0][i] as f32,
                        in_prior: s[2][i].as_ref() == "true",
                        verdict: s[3][i].clone(),
                    })
                    .collect()
            }
            None => Vec::new(),
        };
        Ok(Some(Self {
            prefix: derive_out_prefix(&loaded.file.to_string_lossy()),
            pseudotime,
            types,
            diffusion,
            layout,
            edges,
        }))
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

    /// The figure as SVG, `w × h` pixels.
    pub(crate) fn svg(&self, panel: Panel, w: u32, h: u32) -> Result<String> {
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
            Panel::Order => Ok(self.order_svg(w, h)),
            Panel::Connectivity => Ok(self.connectivity_svg(w, h)),
        }
    }

    /// Cells at `(x, y)` coloured by pseudotime (grey when none), with the
    /// prior's direct edges as arrows between type medians when `arrows`.
    fn scatter(&self, x: &[f32], y: &[f32], w: u32, h: u32, arrows: bool) -> Result<String> {
        let ext = Extent { w, h };
        let finite: Vec<usize> = (0..x.len())
            .filter(|&i| x[i].is_finite() && y[i].is_finite())
            .collect();
        anyhow::ensure!(!finite.is_empty(), "no cell has coordinates");
        let (mut x0, mut x1, mut y0, mut y1) = (f32::MAX, f32::MIN, f32::MAX, f32::MIN);
        for &i in &finite {
            x0 = x0.min(x[i]);
            x1 = x1.max(x[i]);
            y0 = y0.min(y[i]);
            y1 = y1.max(y[i]);
        }
        let bounds = DataBounds::from_minmax(x0, x1, y0, y1);
        let px = |i: usize| crate::plot::to_pixel((x[i], y[i]), &bounds, ext);
        // Unreached cells first, so coloured ones draw over them.
        let mut order = finite.clone();
        order.sort_by_key(|&i| self.pseudotime[i].is_finite());
        let pts: Vec<(f32, f32)> = order.iter().map(|&i| px(i)).collect();
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
        let mut layers = vec![TopicLayer {
            label: String::new(),
            png: rasterize_per_point_png(&pts, &colors, ext, radius, 0.85, PointShape::Circle)?,
            hull_px: Vec::new(),
            label_xy_px: (f32::NAN, f32::NAN),
            color: INK,
        }];
        if arrows {
            let medians = self.type_medians(x, y, &finite);
            let mut strong = Vec::new();
            let mut weak = Vec::new();
            for e in self.edges.iter().filter(|e| e.in_prior) {
                if let (Some(&a), Some(&b)) = (medians.get(e.a.as_ref()), medians.get(e.b.as_ref()))
                {
                    let seg = (
                        crate::plot::to_pixel(a, &bounds, ext),
                        crate::plot::to_pixel(b, &bounds, ext),
                    );
                    if e.verdict.as_ref() == "supported" {
                        strong.push(seg);
                    } else {
                        weak.push(seg);
                    }
                }
            }
            let stroke = radius * 1.2;
            for (segs, alpha) in [(strong, 0.9), (weak, 0.3)] {
                if !segs.is_empty() {
                    layers.push(TopicLayer {
                        label: String::new(),
                        png: rasterize_arrow_layer_png(
                            &segs,
                            ext,
                            stroke,
                            stroke * 5.0,
                            INK,
                            alpha,
                        )?,
                        hull_px: Vec::new(),
                        label_xy_px: (f32::NAN, f32::NAN),
                        color: INK,
                    });
                }
            }
        }
        Ok(emit_svg(
            &layers,
            &SvgOpts {
                width_px: w,
                height_px: h,
                frame_stroke_px: 1.0,
                ..SvgOpts::default()
            },
        ))
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

    /// Per type (with a finite pseudotime): `(type, cells, q1, median, q3)`,
    /// in order of the median.
    pub(crate) fn type_order(&self) -> Vec<(Box<str>, usize, f32, f32, f32)> {
        let mut by_type: BTreeMap<&str, Vec<f32>> = BTreeMap::new();
        for (t, &p) in self.types.iter().zip(&self.pseudotime) {
            if p.is_finite() {
                by_type.entry(t.as_ref()).or_default().push(p);
            }
        }
        let mut rows: Vec<(Box<str>, usize, f32, f32, f32)> = by_type
            .into_iter()
            .map(|(t, mut v)| {
                v.sort_by(f32::total_cmp);
                let q = |f: f32| v[((v.len() - 1) as f32 * f).round() as usize];
                (t.into(), v.len(), q(0.25), median(&v), q(0.75))
            })
            .collect();
        rows.sort_by(|a, b| a.3.total_cmp(&b.3));
        rows
    }

    fn order_svg(&self, w: u32, h: u32) -> String {
        let rows = self.type_order();
        let (left, right, top, bottom) = (w as f32 * 0.36, w as f32 - 16.0, 28.0, h as f32 - 12.0);
        let rh = ((bottom - top) / rows.len().max(1) as f32).min(26.0);
        let x = |v: f32| left + v * (right - left);
        let font = rh.clamp(8.0, 13.0) * 0.9;
        let mut s = svg_head(w, h);
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
        for (k, (t, n, q1, m, q3)) in rows.iter().enumerate() {
            let y = top + rh * (k as f32 + 0.5);
            s += &format!(
                "<text x='{:.1}' y='{:.1}' font-size='{font:.1}' text-anchor='end' fill='#222'>{} ({n})</text>\
                 <line x1='{:.1}' y1='{y:.1}' x2='{:.1}' y2='{y:.1}' stroke='#9ec5f4' stroke-width='{:.1}' stroke-linecap='round'/>\
                 <circle cx='{:.1}' cy='{y:.1}' r='{:.1}' fill='#2a78d6' stroke='white' stroke-width='1.5'/>",
                left - 8.0,
                y + font * 0.35,
                escape_xml(t),
                x(*q1),
                x(*q3),
                (rh * 0.28).max(2.0),
                x(*m),
                (rh * 0.22).clamp(2.5, 5.0)
            );
        }
        s + "</svg>"
    }

    fn connectivity_svg(&self, w: u32, h: u32) -> String {
        let order: Vec<Box<str>> = self.type_order().into_iter().map(|r| r.0).collect();
        let nodes: Vec<&str> = order
            .iter()
            .map(AsRef::as_ref)
            .filter(|t| {
                self.edges
                    .iter()
                    .any(|e| e.a.as_ref() == *t || e.b.as_ref() == *t)
            })
            .collect();
        let n = nodes.len().max(1);
        let (left, top) = (w as f32 * 0.3, 8.0);
        let cell = ((w as f32 - left - 8.0) / n as f32).min((h as f32 - top - 8.0) / n as f32);
        let font = (cell * 0.75).clamp(6.0, 12.0);
        let mut s = svg_head(w, h);
        let at = |t: &str| nodes.iter().position(|&x| x == t);
        let mut m = vec![vec![None::<(f32, bool)>; n]; n];
        for e in &self.edges {
            if let (Some(i), Some(j)) = (at(&e.a), at(&e.b)) {
                m[i][j] = Some((e.connectivity, e.in_prior));
                m[j][i] = Some((e.connectivity, e.in_prior));
            }
        }
        for (i, t) in nodes.iter().enumerate() {
            s += &format!(
                "<text x='{:.1}' y='{:.1}' font-size='{font:.1}' text-anchor='end' fill='#222'>{}</text>",
                left - 6.0,
                top + cell * (i as f32 + 0.5) + font * 0.35,
                escape_xml(t)
            );
            for (j, cell_value) in m[i].iter().enumerate() {
                let (fill, stroke) = match *cell_value {
                    _ if i == j => ("#ffffff".to_string(), "none"),
                    Some((c, in_prior)) => (blue(c), if in_prior { "#222" } else { "none" }),
                    None => ("#f2f4f7".to_string(), "none"),
                };
                s += &format!(
                    "<rect x='{:.1}' y='{:.1}' width='{:.1}' height='{:.1}' fill='{fill}' stroke='{stroke}' stroke-width='1.5'/>",
                    left + cell * j as f32 + 0.5,
                    top + cell * i as f32 + 0.5,
                    cell - 1.0,
                    cell - 1.0
                );
            }
        }
        s + "</svg>"
    }
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

fn svg_head(w: u32, h: u32) -> String {
    format!(
        "<?xml version='1.0' encoding='UTF-8'?>\n<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {w} {h}' width='{w}' height='{h}' font-family='Helvetica, Arial, sans-serif'>\
         <rect width='{w}' height='{h}' fill='white'/>"
    )
}

/// White to blue by `t` in [0, 1].
fn blue(t: f32) -> String {
    let t = t.clamp(0.0, 1.0);
    let lerp = |a: f32, b: f32| (a + (b - a) * t).round() as u8;
    format!(
        "#{:02x}{:02x}{:02x}",
        lerp(255.0, 13.0),
        lerp(255.0, 54.0),
        lerp(255.0, 107.0)
    )
}

/// `svg` rendered to pixels, through legume-plot's renderer (which writes a
/// file) and the `image` crate.
pub(crate) fn render(svg: &str, w: u32, h: u32) -> Result<image::RgbaImage> {
    let tmp = std::env::temp_dir().join(format!("lupin-figure-{}.png", std::process::id()));
    legume_plot::render_png(svg, w, h, &tmp)?;
    let img = image::open(&tmp)
        .with_context(|| format!("reading {}", tmp.display()))?
        .to_rgba8();
    let _ = std::fs::remove_file(&tmp);
    Ok(img)
}

/// Write `svg` as `{base}.svg` and `{base}.pdf`; returns the files written.
pub(crate) fn export(svg: &str, w: u32, h: u32, base: &Path) -> Result<Vec<std::path::PathBuf>> {
    let base_s = base.to_string_lossy();
    legume_plot::write_figure(
        svg,
        w,
        h,
        &base_s,
        legume_plot::FigureFormats {
            svg: true,
            png: false,
            pdf: true,
        },
    )?;
    Ok(vec![
        std::path::PathBuf::from(format!("{base_s}.svg")),
        std::path::PathBuf::from(format!("{base_s}.pdf")),
    ])
}

#[cfg(test)]
#[path = "figures_tests.rs"]
mod tests;
