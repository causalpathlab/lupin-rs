//! `lupin plot` — publication-quality rasterized scatter with vector
//! labels over transparent background.
//!
//! Preferred invocation is `lupin plot --from {prefix}.senna.json`
//! (manifest produced by `senna topic` + enriched by `senna layout`); it
//! fills in cell-coords / topics / labels / colour-by / palette from
//! the manifest's `viz{}`, `outputs{}`, and `defaults{}` sections.
//! Explicit CLI flags still override whatever the manifest provides.
//!
//! See `plot::` module docs for the overall SVG→PNG/PDF pipeline. This
//! file is the clap entry point and glue: resolves manifest + overrides
//! into a `ResolvedInputs`, buckets cells by group, dispatches per-group
//! rasterization via rayon, emits SVG, then renders PNG + PDF.

use crate::annotate::inputs::load_cluster_labels;
use crate::manifest::annotate::{compute_clusters_from_latent, LeidenArgs};
use crate::manifest::run::RunManifest;
use crate::manifest::run::{load_optional, manifest_file, resolve};
use crate::plot::hull::{convex_hull, hull_centroid, median_xy, trim_outliers_by_median, Pt};
use crate::plot::palette::{self, Palette};
use crate::plot::rasterize::{rasterize_group_png, DataBounds, Extent, PointShape};
use crate::plot::svg_emit::{emit_svg, SvgOpts, TopicLayer};
use crate::plot::to_pixel;
use clap::{Args, ValueEnum};
use legume_numeric::matrix::common_io::mkdir_parent;
use legume_numeric::matrix::dense_mat_io::{Mat, MatWithNames};
use legume_numeric::matrix::traits::*;
use log::info;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fs;
use std::path::{Path, PathBuf};

#[cfg(test)]
#[path = "scatter_tests.rs"]
mod tests;

const PT_PER_INCH: f32 = 72.0;

/// Source of the per-cell group ID used for coloring.
#[derive(ValueEnum, Clone, Debug)]
#[clap(rename_all = "kebab-case")]
pub enum ColorBy {
    /// Use the `cluster` column in `cell_coords.parquet`.
    Cluster,
    /// Use the `pb_id` column in `cell_coords.parquet`.
    #[value(alias = "pb_id")]
    PbId,
    /// Argmax of a separate topic-proportions parquet (`--topics`).
    Topic,
    /// Per-cell celltype label from `lupin annotate --method enrichment`'s argmax TSV
    /// (`--annotation`, defaults to `manifest.annotate.argmax`).
    Annotation,
    /// Pseudotime from `lupin trajectory` (the manifest's `trajectory.pseudotime`),
    /// on a blue→red ramp, with the prior's direct edges as arrows between
    /// the types' medians; unsupported edges faded.
    Pseudotime,
}

/// Label placement strategy per group.
#[derive(ValueEnum, Clone, Debug, PartialEq, Eq)]
#[clap(rename_all = "kebab-case")]
pub enum LabelPosition {
    /// Coordinate-wise median of the group's points (robust default).
    Median,
    /// Area-weighted centroid of the group's convex hull.
    HullCentroid,
}

#[derive(Args, Debug)]
pub struct PlotArgs {
    #[arg(
        long,
        short = 'f',
        help = "Run manifest JSON (+ updated by `senna layout`)",
        long_help = "Run manifest JSON from any senna embedding run (+ updated by `senna layout`).\n\
                     \n\
                     If set, fills in --cell-coords, --topics, --labels, --colour-by,\n\
                     and --palette from the manifest's viz/outputs/defaults sections.\n\
                     Any explicit CLI flag still overrides the manifest value.\n\
                     Paths inside the manifest resolve relative to its own directory,\n\
                     so you can move a run directory around freely."
    )]
    pub from: Option<Box<str>>,

    #[arg(
        long,
        short = 'c',
        help = "Cell coordinates parquet (from `senna layout`)",
        long_help = "Cell coordinates parquet (from `senna layout`);\n\
                     required unless --from provides it."
    )]
    pub cell_coords: Option<Box<str>>,

    #[arg(
        long,
        short = 'o',
        help = "Output prefix",
        long_help = "Writes {out}.plot.pdf.\n\
                     Pass --svg / --png to additionally emit {out}.plot.svg / {out}.plot.png;\n\
                     --no-pdf suppresses the PDF."
    )]
    pub out: Box<str>,

    #[arg(
        long = "colour-by",
        alias = "color-by",
        value_enum,
        help = "Colour source (default: manifest's `defaults.colour_by`, else `cluster`)"
    )]
    pub colour_by: Option<ColorBy>,

    #[arg(
        long,
        help = "Topic proportions parquet (cells × K); required with --colour-by topic"
    )]
    pub topics: Option<Box<str>>,

    #[arg(
        long,
        help = "Annotation argmax TSV from `lupin annotate --method enrichment`",
        long_help = "Annotation argmax TSV from `lupin annotate --method enrichment`.\n\
                     Columns are cell\\tcell_type\\tprobability.\n\
                     Defaults to the manifest's annotate.argmax."
    )]
    pub annotation: Option<Box<str>>,

    #[arg(
        long,
        default_value_t = false,
        help = "Preload data when auto-running `senna layout`",
        long_help = "Preload data when auto-running `senna layout`,\n\
                     for a manifest missing layout.cell_coords.\n\
                     No-op if cell_coords already exists.",
        hide = true
    )]
    pub preload_data: bool,

    #[arg(
        long,
        help = "TSV mapping group_id<TAB>display_name (one per line)",
        long_help = "TSV mapping group_id<TAB>display_name (one per line).\n\
                     Missing IDs fall back to T{id}."
    )]
    pub labels: Option<Box<str>>,

    #[arg(
        long,
        help = "Drop groups with fewer than N cells (0 = keep all)",
        long_help = "Filter out small/dead groups before rendering. When unset,\n\
                     defaults to max(50, n_cells / 200) for --colour-by topic,\n\
                     which kills argmax ghosts on dead topics, and 0 otherwise.\n\
                     Pass --min-topic-cells 0 to opt out of the auto threshold."
    )]
    pub min_topic_cells: Option<usize>,

    #[arg(long, default_value_t = 6.0, help = "Plot width (inches)")]
    pub width: f32,

    #[arg(long, default_value_t = 6.0, help = "Plot height (inches)")]
    pub height: f32,

    #[arg(long, default_value_t = 300, help = "Output DPI (raster layers)")]
    pub dpi: u32,

    #[arg(long, default_value_t = 2.0, help = "Point size (pt)")]
    pub point_size: f32,

    #[arg(long, default_value_t = 0.6, help = "Point alpha (0..=1)")]
    pub alpha: f32,

    #[arg(
        long,
        value_enum,
        default_value_t = PointShape::Circle,
        help = "Marker shape"
    )]
    pub point_shape: PointShape,

    #[arg(
        long,
        default_value_t = false,
        help = "Cycle marker shape per group (circle→triangle→square→diamond)"
    )]
    pub point_shape_cycle: bool,

    #[arg(
        long,
        value_enum,
        help = "Qualitative palette (default: manifest's `defaults.palette`, else `auto`)"
    )]
    pub palette: Option<Palette>,

    #[arg(
        long,
        value_enum,
        default_value_t = LabelPosition::Median,
        help = "Label placement strategy"
    )]
    pub label_position: LabelPosition,

    #[arg(long, default_value_t = 10.0, help = "Label font size (pt)")]
    pub label_font_size: f32,

    #[arg(long, default_value_t = false, help = "Suppress vector text labels")]
    pub no_labels: bool,

    #[arg(
        long,
        default_value_t = false,
        help = "Draw convex hull polygons around each group",
        long_help = "Draw convex hull polygons around each group. Off by default:\n\
                     scRNA groups are rarely separable in 2D, and hulls overstate that."
    )]
    pub hull: bool,

    #[arg(
        long,
        default_value_t = 0.95,
        help = "Fraction of closest-to-median points used for each hull (1.0 = all)",
        long_help = "Only applies when --hull is enabled.\n\
                     For each group, keep only the points nearest the median.\n\
                     The median is coordinate-wise; the distance is Euclidean.\n\
                     {coverage} is the fraction kept, before computing the hull.\n\
                     Strips a few fringe cells so one outlier can't drag the polygon.\n\
                     Set to 1.0 to use every point."
    )]
    pub hull_coverage: f32,

    #[arg(
        long,
        default_value_t = 0.0,
        help = "Hull fill opacity (0..=1; 0 = outline only)"
    )]
    pub hull_fill_alpha: f32,

    #[arg(
        long,
        default_value_t = false,
        help = "Also emit SVG (default: PDF only)"
    )]
    pub svg: bool,

    #[arg(
        long,
        default_value_t = false,
        help = "Also emit flattened PNG (default: PDF only)"
    )]
    pub png: bool,

    #[arg(long, default_value_t = false, help = "Skip PDF output")]
    pub no_pdf: bool,
}

pub fn fit_plot(args: &PlotArgs) -> anyhow::Result<()> {
    let mut resolved = resolve_inputs(args)?;
    mkdir_parent(&resolved.out)?;
    if matches!(resolved.colour_by, ColorBy::Pseudotime) {
        match trajectory_figure(&resolved)? {
            Some(data) => return plot_pseudotime(args, &resolved, &data),
            None if args.colour_by.is_some() => anyhow::bail!(
                "--colour-by pseudotime needs a manifest with a `trajectory` section: run \
                 `lupin trajectory` first"
            ),
            None => {
                log::warn!(
                    "the manifest's default colour is pseudotime but it has no trajectory section; \
                     colouring by annotation or cluster instead"
                );
                resolved.colour_by = if resolved.annotation.is_some() {
                    ColorBy::Annotation
                } else {
                    ColorBy::Cluster
                };
            }
        }
    }
    let (cell_names, coords_by_name) = read_cell_coords(&resolved.cell_coords)?;
    let x = coords_by_name
        .get("x")
        .ok_or_else(|| anyhow::anyhow!("cell_coords parquet missing 'x' column"))?;
    let y = coords_by_name
        .get("y")
        .ok_or_else(|| anyhow::anyhow!("cell_coords parquet missing 'y' column"))?;
    let n_cells = x.len();
    info!("Loaded {n_cells} cells from {}", resolved.cell_coords);

    // Auto-built id → display name map for `ColorBy::Annotation` (replaces
    // the empty default below). The `--labels` TSV, when present, layers
    // on top so users can still rename specific celltypes.
    let mut auto_label_map: FxHashMap<i64, String> = FxHashMap::default();

    let group_ids = match resolved.colour_by {
        // Drawn by the trajectory's own figure above, never grouped.
        ColorBy::Pseudotime => anyhow::bail!("pseudotime colouring is not a grouping"),
        ColorBy::Cluster => {
            let ids = match coords_by_name.get("cluster") {
                Some(col) => col
                    .iter()
                    .map(|&v| if v.is_nan() || v < 0.0 { -1 } else { v as i64 })
                    .collect::<Vec<_>>(),
                None => resolve_cluster_ids_for_plot(&mut resolved, &cell_names)?,
            };
            // Label cluster groups as `C{g}` so a downstream `--labels` TSV
            // can still rename them by id, but the default isn't the
            // topic-style `T{g}`.
            for &g in &ids {
                if g >= 0 {
                    auto_label_map.entry(g).or_insert_with(|| format!("C{g}"));
                }
            }
            ids
        }
        ColorBy::PbId => coords_by_name
            .get("pb_id")
            .ok_or_else(|| anyhow::anyhow!("cell_coords missing 'pb_id' column"))?
            .iter()
            .map(|&v| if v.is_nan() { -1 } else { v as i64 })
            .collect::<Vec<_>>(),
        ColorBy::Annotation => {
            let path = resolved.annotation.as_ref().ok_or_else(|| {
                anyhow::anyhow!(
                    "--colour-by annotation requires --annotation PATH \
                     (or `lupin annotate --method enrichment` must have populated manifest.annotate.argmax)"
                )
            })?;
            let (ids, label_map) = argmax_annotation(path, &cell_names)?;
            auto_label_map = label_map;
            ids
        }
        ColorBy::Topic => {
            let topics_path = resolved.topics.as_ref().ok_or_else(|| {
                anyhow::anyhow!(
                    "--colour-by topic requires --topics PATH (or manifest outputs.topics)"
                )
            })?;
            argmax_topics(topics_path, n_cells)?
        }
    };

    // Single pass: per-group bucketing + global bounds.
    let mut pts_by_group: FxHashMap<i64, Vec<Pt>> = FxHashMap::default();
    let (mut xmin, mut xmax) = (f32::INFINITY, f32::NEG_INFINITY);
    let (mut ymin, mut ymax) = (f32::INFINITY, f32::NEG_INFINITY);
    for i in 0..n_cells {
        let g = group_ids[i];
        if g < 0 || !x[i].is_finite() || !y[i].is_finite() {
            continue;
        }
        pts_by_group.entry(g).or_default().push((x[i], y[i]));
        if x[i] < xmin {
            xmin = x[i];
        }
        if x[i] > xmax {
            xmax = x[i];
        }
        if y[i] < ymin {
            ymin = y[i];
        }
        if y[i] > ymax {
            ymax = y[i];
        }
    }
    if pts_by_group.is_empty() {
        anyhow::bail!("no valid group assignments found");
    }

    // Auto-threshold for `--colour-by topic` kills argmax ghosts on dead
    // topics. Floor of 50 protects small-N runs.
    let min_topic_cells: usize = args.min_topic_cells.unwrap_or_else(|| {
        if matches!(resolved.colour_by, ColorBy::Topic) {
            (n_cells / 200).max(50)
        } else {
            0
        }
    });
    if min_topic_cells > 0 {
        let before = pts_by_group.len();
        pts_by_group.retain(|_g, pts| pts.len() >= min_topic_cells);
        let dropped = before - pts_by_group.len();
        if dropped > 0 {
            info!("Dropped {dropped} groups with fewer than {min_topic_cells} cells");
        }
        if pts_by_group.is_empty() {
            anyhow::bail!("--min-topic-cells {min_topic_cells} dropped every group");
        }
    }

    let mut unique: Vec<i64> = pts_by_group.keys().copied().collect();
    unique.sort_unstable();
    info!("Plotting {} groups", unique.len());

    let width_px = (args.width * args.dpi as f32).round() as u32;
    let height_px = (args.height * args.dpi as f32).round() as u32;
    let ext = Extent {
        w: width_px,
        h: height_px,
    };
    let bounds = DataBounds::from_minmax(xmin, xmax, ymin, ymax);

    let radius_px = args.point_size * args.dpi as f32 / PT_PER_INCH / 2.0;
    let label_font_px = args.label_font_size * args.dpi as f32 / PT_PER_INCH;
    let palette = palette::resolve(&resolved.palette, unique.len());

    // Start from the auto-built map (populated for ColorBy::Annotation,
    // empty otherwise), then layer the user-supplied --labels TSV on
    // top so per-id renames win over the auto names.
    let mut label_map = auto_label_map;
    if let Some(p) = resolved.labels.as_deref() {
        for (id, name) in read_labels_tsv(p)? {
            label_map.insert(id, name);
        }
    }

    // Skip hull/trim computation when nothing needs it — hull polygons
    // are opt-in now and median is the default label position.
    let need_hull_geometry = args.hull || args.label_position == LabelPosition::HullCentroid;

    let layers: Vec<TopicLayer> = unique
        .par_iter()
        .map(|g| -> anyhow::Result<TopicLayer> {
            let pts = pts_by_group.get(g).expect("group present");
            let pts_px: Vec<(f32, f32)> = pts.iter().map(|&p| to_pixel(p, &bounds, ext)).collect();

            // Key palette by group id, not enumerate index, so colors
            // stay aligned across every view that maps id → color.
            let color_idx = (*g).max(0) as usize;
            let color = palette::color(&palette, color_idx);
            let shape = if args.point_shape_cycle {
                PointShape::cycle_nth(color_idx)
            } else {
                args.point_shape
            };
            let png = rasterize_group_png(
                &pts_px,
                ext,
                legume_plot::RadiusSpec::Scalar(radius_px),
                color,
                args.alpha,
                shape,
            )?;

            let (hull_px, hull_data) = if need_hull_geometry {
                let trimmed = trim_outliers_by_median(pts, args.hull_coverage);
                let h = convex_hull(&trimmed);
                let px: Vec<Pt> = h.iter().map(|&p| to_pixel(p, &bounds, ext)).collect();
                (px, h)
            } else {
                (Vec::new(), Vec::new())
            };

            let label_xy_data = match args.label_position {
                LabelPosition::Median => median_xy(pts),
                LabelPosition::HullCentroid => hull_centroid(&hull_data),
            };
            let label_xy_px = to_pixel(label_xy_data, &bounds, ext);

            let label = label_map.get(g).cloned().unwrap_or_else(|| format!("T{g}"));

            Ok(TopicLayer {
                label,
                png,
                hull_px,
                label_xy_px,
                color,
            })
        })
        .collect::<anyhow::Result<Vec<_>>>()?;

    let svg = emit_svg(
        &layers,
        &SvgOpts {
            width_px,
            height_px,
            draw_hulls: args.hull,
            draw_labels: !args.no_labels,
            label_font_size_px: label_font_px,
            hull_stroke_px: (radius_px * 0.8).max(1.0),
            hull_fill_alpha: args.hull_fill_alpha,
            ..Default::default()
        },
    );

    let base = resolved.out.clone();
    legume_plot::write_figure(
        &svg,
        width_px,
        height_px,
        &format!("{base}.plot"),
        legume_plot::FigureFormats {
            svg: args.svg,
            png: args.png,
            pdf: !args.no_pdf,
        },
    )?;

    Ok(())
}

/// Resolved inputs after merging `--from` manifest defaults with
/// explicit CLI flags. CLI wins; missing-and-not-in-manifest yields a
/// clear error. All paths here are absolute (or at least resolved
/// relative to the manifest's directory) so downstream readers don't
/// have to guess what they're relative to.
struct ResolvedInputs {
    cell_coords: String,
    topics: Option<String>,
    annotation: Option<String>,
    /// Cluster assignments parquet path (cells × 1, "cluster" column).
    /// Used by the `ColorBy::Cluster` fallback when `cell_coords` lacks
    /// a `cluster` column. Populated from `manifest.cluster.clusters`
    /// when `--from` is in effect.
    clusters: Option<String>,
    labels: Option<String>,
    out: String,
    colour_by: ColorBy,
    palette: Palette,
    /// Loaded manifest (when `--from` is in effect). The cluster fallback
    /// may write back to it (`cluster.clusters` field) and persist to
    /// `manifest_path`.
    manifest: Option<RunManifest>,
    manifest_path: Option<String>,
    manifest_dir: PathBuf,
}

fn resolve_inputs(args: &PlotArgs) -> anyhow::Result<ResolvedInputs> {
    let (manifest, manifest_dir) = load_optional(args.from.as_deref())?;

    let resolve_opt = |s: &str| resolve(&manifest_dir, s);

    let cell_coords = args
        .cell_coords
        .as_deref()
        .map(String::from)
        .or_else(|| {
            manifest
                .as_ref()
                .and_then(|m| m.layout.cell_coords.as_deref())
                .map(resolve_opt)
        })
        .ok_or_else(|| {
            anyhow::anyhow!(
                "no --cell-coords given and manifest {} has no layout.cell_coords; \
                 run `senna layout umap --from …` (or phate/tsne) first",
                args.from.as_deref().unwrap_or("(none)")
            )
        })?;

    // Log θ when the run resolved topics, else the cell embedding: argmax
    // over its axes is the natural "dominant direction" colouring.
    let topics = args.topics.as_deref().map(String::from).or_else(|| {
        manifest
            .as_ref()
            .and_then(|m| m.outputs.structure_latent())
            .map(resolve_opt)
    });

    let annotation = args.annotation.as_deref().map(String::from).or_else(|| {
        manifest
            .as_ref()
            .and_then(|m| m.annotate.argmax.as_deref())
            .map(resolve_opt)
    });

    let clusters = manifest
        .as_ref()
        .and_then(|m| m.cluster.clusters.as_deref())
        .map(resolve_opt);

    let labels = args.labels.as_deref().map(String::from).or_else(|| {
        manifest
            .as_ref()
            .and_then(|m| m.outputs.anchor_labels.as_deref())
            .map(resolve_opt)
    });

    let out = args.out.to_string();

    // Colour-by precedence:
    //   1. explicit CLI flag
    //   2. manifest `defaults.colour_by`
    //   3. Annotation, if `manifest.annotate.argmax` is populated
    //   4. Cluster (existing fallback)
    // Step 3 makes annotation the natural default for the
    // train → annotate → plot workflow without forcing users to set it.
    let colour_by = args
        .colour_by
        .clone()
        .or_else(|| {
            manifest
                .as_ref()
                .and_then(|m| m.defaults.colour_by.as_deref())
                .and_then(parse_colour_by)
        })
        .or_else(|| annotation.as_ref().map(|_| ColorBy::Annotation))
        .unwrap_or(ColorBy::Cluster);

    let palette = args
        .palette
        .clone()
        .or_else(|| {
            manifest
                .as_ref()
                .and_then(|m| m.defaults.palette.as_deref())
                .and_then(parse_palette)
        })
        .unwrap_or(Palette::Auto);

    Ok(ResolvedInputs {
        cell_coords,
        topics,
        annotation,
        clusters,
        labels,
        out,
        colour_by,
        palette,
        manifest,
        manifest_path: args
            .from
            .as_deref()
            .map(|f| manifest_file(f).to_string_lossy().into_owned()),
        manifest_dir,
    })
}

/// A manifest default lupin does not know (hand-edited, or a mode since
/// removed, e.g. `pseudotime`) is reported and skipped, not silently ignored.
/// The trajectory outputs the manifest records, when `--from` gave one with a
/// `trajectory` section.
fn trajectory_figure(
    resolved: &ResolvedInputs,
) -> anyhow::Result<Option<crate::trajectory::figures::TrajectoryData>> {
    let (Some(manifest), Some(path)) = (&resolved.manifest, &resolved.manifest_path) else {
        return Ok(None);
    };
    let loaded = crate::manifest::run::Loaded {
        manifest: manifest.clone(),
        dir: resolved.manifest_dir.clone(),
        file: PathBuf::from(path),
    };
    crate::trajectory::figures::TrajectoryData::load(&loaded)
}

/// The layout coloured by pseudotime, through the trajectory's own figure.
fn plot_pseudotime(
    args: &PlotArgs,
    resolved: &ResolvedInputs,
    data: &crate::trajectory::figures::TrajectoryData,
) -> anyhow::Result<()> {
    let width_px = (args.width * args.dpi as f32).round() as u32;
    let height_px = (args.height * args.dpi as f32).round() as u32;
    let svg = data.svg(
        crate::trajectory::figures::Panel::Layout,
        width_px,
        height_px,
    )?;
    legume_plot::write_figure(
        &svg,
        width_px,
        height_px,
        &format!("{}.plot", resolved.out),
        legume_plot::FigureFormats {
            svg: args.svg,
            png: args.png,
            pdf: !args.no_pdf,
        },
    )?;
    Ok(())
}

fn parse_colour_by(s: &str) -> Option<ColorBy> {
    use clap::ValueEnum;
    let parsed = ColorBy::from_str(s, true).ok();
    if parsed.is_none() {
        log::warn!("ignoring unknown defaults.colour_by '{s}' in the manifest");
    }
    parsed
}

fn parse_palette(s: &str) -> Option<Palette> {
    use clap::ValueEnum;
    Palette::from_str(s, true).ok()
}

pub(crate) type CellCoords = (Vec<Box<str>>, FxHashMap<String, Vec<f32>>);

/// Returns `(cell_names, columns_by_name)`. Cell names are the parquet
/// row labels (in data column order), needed when matching against an
/// annotation TSV by cell name.
pub(crate) fn read_cell_coords(path: &str) -> anyhow::Result<CellCoords> {
    let MatWithNames { rows, cols, mat } = Mat::from_parquet(path)?;
    let mut by_name: FxHashMap<String, Vec<f32>> = FxHashMap::default();
    for (j, name) in cols.iter().enumerate() {
        let col: Vec<f32> = (0..mat.nrows()).map(|i| mat[(i, j)]).collect();
        by_name.insert(name.to_string(), col);
    }
    Ok((rows, by_name))
}

/// Argmax over cells × K, returning the **axis ID** per row rather than
/// the column position: a topic table names its columns `T{c}`, so a
/// downstream `T5` legend swatch always means topic 5 even if the columns
/// aren't in 0..K-1 order on disk. An embedding table (`h{c}`) has no IDs of
/// its own, so its axes are numbered by position.
fn argmax_topics(path: &str, n_cells_expected: usize) -> anyhow::Result<Vec<i64>> {
    let MatWithNames { cols, mat, .. } = Mat::from_parquet(path)?;
    if mat.nrows() != n_cells_expected {
        anyhow::bail!(
            "topics parquet has {} rows but cell_coords has {}",
            mat.nrows(),
            n_cells_expected
        );
    }
    // `T{c}` carries an ID; `h{c}` is an embedding axis, numbered by position;
    // anything else is refused rather than argmaxed as if it were a composition.
    let topic_ids = crate::plot::try_parse_axis_ids(&cols, "T")
        .or_else(|| crate::plot::try_parse_axis_ids(&cols, "h"))
        .ok_or_else(|| {
            anyhow::anyhow!(
                "topics parquet at {path} has columns that are neither topic IDs (`T{{c}}`) \
                 nor embedding axes (`h{{c}}`)"
            )
        })?;
    let mut out = Vec::with_capacity(mat.nrows());
    for i in 0..mat.nrows() {
        let mut best_j = 0usize;
        let mut best_v = f32::NEG_INFINITY;
        for j in 0..mat.ncols() {
            let v = mat[(i, j)];
            if v > best_v {
                best_v = v;
                best_j = j;
            }
        }
        out.push(topic_ids[best_j]);
    }
    Ok(out)
}

/// Cluster fallback for `ColorBy::Cluster` when the `cell_coords`
/// parquet has no `cluster` column. Mirrors `lupin annotate --method enrichment`'s strategy:
///   1. If `manifest.cluster.clusters` is set, load that parquet.
///   2. Otherwise run Leiden against the manifest's cell embedding
///      (`outputs.cell_embedding`, else `outputs.latent`), write
///      `{out}.clusters.parquet`, and update the manifest in place.
///
/// Returns per-cell ids aligned to `cell_names` (`cell_coords` row order),
/// with `-1` for unassigned/missing cells.
fn resolve_cluster_ids_for_plot(
    resolved: &mut ResolvedInputs,
    cell_names: &[Box<str>],
) -> anyhow::Result<Vec<i64>> {
    let Some(manifest) = resolved.manifest.as_mut() else {
        anyhow::bail!(
            "--colour-by cluster but cell_coords has no 'cluster' column \
             and no manifest is loaded; re-run `senna layout` with --clusters \
             or pass --from <manifest>"
        );
    };

    // Path 1: manifest already has a cluster parquet.
    if let Some(path) = resolved.clusters.as_deref() {
        info!("Loading clusters from {path}");
        let (labels_usize, n_clusters) = load_cluster_labels(path, cell_names)?;
        info!("Loaded {n_clusters} clusters from manifest.cluster.clusters");
        return Ok(usize_to_signed(&labels_usize));
    }

    // Path 2: leiden on the manifest's latent, with the annotate defaults.
    info!(
        "No 'cluster' column in cell_coords and manifest.cluster.clusters is unset; \
         running internal Leiden on the manifest latent"
    );
    let (labels_usize, n_clusters) = compute_clusters_from_latent(
        manifest,
        &resolved.manifest_dir,
        cell_names,
        &LeidenArgs::default(),
    )?;
    info!("Internal Leiden produced {n_clusters} clusters");

    // Persist `{out}.clusters.parquet` next to the plot output and patch
    // the manifest so future runs (annotate/plot) reuse it.
    let parquet_path = format!("{}.clusters.parquet", resolved.out);
    write_cluster_assignments_parquet(&parquet_path, cell_names, &labels_usize)?;
    info!("Wrote {parquet_path}");

    if let Some(manifest_path) = resolved.manifest_path.as_deref() {
        let rel = crate::manifest::run::rel_to_manifest(&resolved.manifest_dir, &parquet_path);
        manifest.cluster.clusters = Some(rel);
        // Only a convenience for later runs: not worth asking about.
        if crate::manifest::run::may_replace_unasked(Path::new(manifest_path)) {
            manifest.save(Path::new(manifest_path))?;
        } else {
            log::info!("clusters not recorded in {manifest_path} (--overwrite records them)");
        }
    }
    resolved.clusters = Some(parquet_path);

    Ok(usize_to_signed(&labels_usize))
}

fn usize_to_signed(labels: &[usize]) -> Vec<i64> {
    labels
        .iter()
        .map(|&v| if v == usize::MAX { -1 } else { v as i64 })
        .collect()
}

fn write_cluster_assignments_parquet(
    path: &str,
    cell_names: &[Box<str>],
    labels: &[usize],
) -> anyhow::Result<()> {
    use legume_numeric::matrix::traits::IoOps;
    let mut data = Mat::zeros(cell_names.len(), 1);
    for (i, &c) in labels.iter().enumerate() {
        data[(i, 0)] = if c == usize::MAX { f32::NAN } else { c as f32 };
    }
    let cols: Vec<Box<str>> = vec!["cluster".into()];
    data.to_parquet_with_names(path, (Some(cell_names), Some("cell")), Some(&cols))?;
    Ok(())
}

/// Read `lupin annotate --method enrichment`'s argmax TSV and produce per-cell integer
/// group ids + a stable id → celltype-name label map.
///
/// Format: `cell\tcell_type\tprobability` with optional header. Cells
/// absent from the TSV map to `-1` (filtered downstream by the same
/// `g < 0` skip used for unassigned clusters / NaN `pb_ids`). Celltype
/// strings are sorted alphabetically before id assignment so the same
/// celltype gets the same colour across reruns and across sibling
/// commands (e.g. `plot-topic --group-by annotation`, which sorts
/// celltype panels alphabetically too — see plot/topic/mod.rs:533).
fn argmax_annotation(
    path: &str,
    cell_names: &[Box<str>],
) -> anyhow::Result<(Vec<i64>, FxHashMap<i64, String>)> {
    let mem = crate::cell_labels::read_cell_labels(path)?;

    // Stable id assignment: sort unique celltype strings, assign 0..N.
    let mut unique: Vec<Box<str>> = mem.unique_groups();
    // Push "unassigned" to the very end of the colour cycle when present,
    // so the argmax-thresholded cells don't steal the leading palette
    // slot.
    if let Some(pos) = unique
        .iter()
        .position(|s| s.as_ref() == enrichment::UNASSIGNED_LABEL)
    {
        let last = unique.remove(pos);
        unique.push(last);
    }
    let name_to_id: FxHashMap<Box<str>, i64> = unique
        .iter()
        .enumerate()
        .map(|(i, n)| (n.clone(), i as i64))
        .collect();

    let mut group_ids = Vec::with_capacity(cell_names.len());
    let mut n_missing = 0usize;
    for c in cell_names {
        if let Some(label) = mem.get(c) {
            group_ids.push(*name_to_id.get(label).unwrap_or(&-1))
        } else {
            group_ids.push(-1);
            n_missing += 1;
        }
    }
    if n_missing > 0 {
        info!(
            "annotation: {n_missing}/{} cells absent from {path} → dropped from plot",
            cell_names.len()
        );
    }

    let label_map: FxHashMap<i64, String> = name_to_id
        .into_iter()
        .map(|(name, id)| (id, String::from(name)))
        .collect();
    Ok((group_ids, label_map))
}

/// Parse a two-column TSV of `group_id<TAB>display_name`. Blank lines
/// and lines starting with `#` are skipped.
fn read_labels_tsv(path: &str) -> anyhow::Result<FxHashMap<i64, String>> {
    let content = fs::read_to_string(Path::new(path))?;
    let mut map: FxHashMap<i64, String> = FxHashMap::default();
    for (line_no, raw) in content.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.splitn(2, '\t');
        let id_str = parts
            .next()
            .ok_or_else(|| anyhow::anyhow!("labels TSV line {}: missing id", line_no + 1))?;
        let name = parts
            .next()
            .ok_or_else(|| anyhow::anyhow!("labels TSV line {}: missing name", line_no + 1))?;
        let id: i64 = id_str
            .parse()
            .map_err(|e| anyhow::anyhow!("labels TSV line {}: bad id: {e}", line_no + 1))?;
        map.insert(id, name.to_string());
    }
    Ok(map)
}
