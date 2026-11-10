//! `lupin figure round`: a round's charts and cell maps as files, without
//! the screen, as the TUI's figure pane draws them: the heatmaps of its top
//! GO/GMT terms, its cell types and its markers × clusters (each with its
//! numbers beside it, as JSON), the structure plots of its cells' cell
//! types or topics, and its cells on any layout in any colouring.

use super::charts::{HeatKind, Mixture};
use super::look::Look;
use super::maps::{once_each, Maps};
use super::units::{STARS, STAR_Q};
use super::{CellFigures, Drawing, Panel};
use anyhow::{Context, Result};
use clap::ValueEnum;
use serde::Serialize;
use std::path::PathBuf;

/// A round's chart, by name.
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChartKind {
    /// Each cluster's top GO/GMT terms × clusters.
    Terms,
    /// The cell types × clusters.
    Types,
    /// The panel's markers × clusters.
    Markers,
    /// Each cell's cell types, stacked, a panel per cluster.
    Structure,
    /// Each cell's topics (a topic run's), stacked, a panel per cluster.
    Topics,
}

impl ChartKind {
    /// Its name, as the command line takes it.
    pub(crate) fn name(self) -> String {
        (self.to_possible_value())
            .expect("no variant is skipped")
            .get_name()
            .to_string()
    }

    fn panel(self) -> Panel {
        match self {
            Self::Terms => Panel::Heatmap {
                of: HeatKind::Terms,
            },
            Self::Types => Panel::Heatmap {
                of: HeatKind::Types,
            },
            Self::Markers => Panel::Heatmap {
                of: HeatKind::Markers,
            },
            Self::Structure => Panel::Structure { of: Mixture::Types },
            Self::Topics => Panel::Structure {
                of: Mixture::Topics,
            },
        }
    }
}

/// What `lupin figure round` draws: its charts, then its maps.
#[derive(Debug)]
pub(crate) struct RoundDrawing {
    pub(crate) charts: Vec<ChartKind>,
    pub(crate) maps: Drawing,
}

/// A chart a round offers, or why it cannot draw it.
#[derive(Debug, Serialize)]
pub(crate) struct ChartOffered {
    pub(crate) chart: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) missing: Option<String>,
}

impl CellFigures {
    /// Each chart, and why the round cannot draw it when it cannot.
    pub(crate) fn charts_offered(&self) -> Vec<ChartOffered> {
        (ChartKind::value_variants().iter())
            .map(|&c| ChartOffered {
                chart: c.name(),
                missing: self.chart_missing(c.panel()),
            })
            .collect()
    }

    /// What to draw: the charts `charts` names, which the round must
    /// record, and the maps `maps` names in `colour`; neither, every chart
    /// the round records.
    pub(crate) fn round_drawing<S: AsRef<str>>(
        &self,
        charts: &[ChartKind],
        maps: &[S],
        colour: Option<&str>,
        look: &Look,
    ) -> Result<RoundDrawing> {
        let charts: Vec<ChartKind> = if charts.is_empty() && maps.is_empty() {
            let have: Vec<ChartKind> = (ChartKind::value_variants().iter().copied())
                .filter(|c| self.chart_missing(c.panel()).is_none())
                .collect();
            anyhow::ensure!(
                !have.is_empty(),
                "the round records no chart: {}",
                (self.charts_offered().iter())
                    .map(|o| format!("{}: {}", o.chart, o.missing.as_deref().unwrap_or("")))
                    .collect::<Vec<_>>()
                    .join("; ")
            );
            have
        } else {
            if let Some((c, why)) =
                (charts.iter()).find_map(|&c| Some((c, self.chart_missing(c.panel())?)))
            {
                anyhow::bail!("no {} chart: {why}", c.name());
            }
            once_each(charts.to_vec())
        };
        let maps = Drawing {
            maps: Maps::These(self.maps_named(maps)?),
            style: self.style_in(colour, look)?,
            look: look.clone(),
        };
        Ok(RoundDrawing { charts, maps })
    }

    /// Each chart's files under `base`: `{base}.{chart}.png` and `.pdf`, a
    /// heatmap's numbers as `.json` too.
    fn chart_files(panel: Panel, base: &str) -> Vec<PathBuf> {
        let b = format!("{base}.{}", panel.slug());
        let mut exts = vec!["png", "pdf"];
        if matches!(panel, Panel::Heatmap { .. }) {
            exts.push("json");
        }
        exts.into_iter()
            .map(|e| PathBuf::from(format!("{b}.{e}")))
            .collect()
    }

    /// The files [`Self::write_round`] writes under `base`: the charts',
    /// then the maps' ([`Self::drawing_files`]).
    pub(crate) fn round_files(&self, d: &RoundDrawing, base: &str) -> Vec<PathBuf> {
        let mut out: Vec<PathBuf> = (d.charts.iter())
            .flat_map(|c| Self::chart_files(c.panel(), base))
            .collect();
        out.extend(self.drawing_files(&d.maps, base));
        out
    }

    /// A round's charts and maps as files under `base`
    /// ([`Self::round_files`]).
    pub(crate) fn write_round(&self, d: &RoundDrawing, base: &str) -> Result<Vec<PathBuf>> {
        legume_numeric::matrix::common_io::mkdir_parent(base)?;
        for panel in d.charts.iter().map(|c| c.panel()) {
            let files = Self::chart_files(panel, base);
            let (img, _) = self.write_chart_fitted(panel, &d.maps.look, &files[1])?;
            (img.save(&files[0])).with_context(|| format!("writing {}", files[0].display()))?;
            if let (Panel::Heatmap { of }, Some(json)) = (panel, files.get(2)) {
                let heat = self.heat(of).map_err(|e| anyhow::anyhow!(e))?;
                std::fs::write(json, serde_json::to_string_pretty(&*heat)?)
                    .with_context(|| format!("writing {}", json.display()))?;
            }
        }
        self.write_drawing(&d.maps, base)?;
        Ok(self.round_files(d, base))
    }

    /// The charts `d` turns onto their page, being wider than it, by name.
    pub(crate) fn round_turned(&self, d: &RoundDrawing) -> Result<Vec<String>> {
        let turned = |c: &ChartKind| {
            Ok(self
                .chart_fit(c.panel(), &d.maps.look)?
                .turn
                .map(|_| c.name()))
        };
        (d.charts.iter())
            .filter_map(|c| turned(c).transpose())
            .collect()
    }

    /// Each heatmap `d` draws under `base`, by name, and the file its
    /// numbers go to.
    pub(crate) fn round_numbers(d: &RoundDrawing, base: &str) -> Vec<(String, PathBuf)> {
        (d.charts.iter())
            .filter_map(|c| Some((c.name(), Self::chart_files(c.panel(), base).get(2)?.clone())))
            .collect()
    }
}

#[derive(clap::Args, Debug)]
pub struct RoundFigureArgs {
    #[arg(short, long, help = "A round (its manifest or prefix)")]
    pub from: Box<str>,
    #[arg(
        short,
        long,
        help = "Output prefix: each chart as {out}.heatmap_{terms|types|markers} or {out}.structure_{types|topics}, .png and .pdf (a heatmap's numbers as .json too); each map as {out}.layout_{name}.{colour}.png and .pdf"
    )]
    pub out: Box<str>,
    #[arg(
        long = "chart",
        value_enum,
        help = "Draw these charts (every one the round records when neither --chart nor --map is given); repeat for several"
    )]
    pub charts: Vec<ChartKind>,
    #[arg(
        long = "map",
        value_name = "MAP",
        help = "Draw the round's cells on these layouts (phate, umap, …, best or all); repeat for several"
    )]
    pub maps: Vec<Box<str>>,
    #[arg(
        long,
        visible_alias = "color",
        help = "Colour the map's cells by type (default), fine (a coarse first round's fine labels), cluster, call or margin"
    )]
    pub colour: Option<Box<str>>,
    #[command(flatten)]
    pub look: Look,
    #[arg(long, default_value_t = STARS, help = "The terms heatmap: each cluster's top terms")]
    pub top: usize,
    #[arg(long, default_value_t = STAR_Q, help = "The terms heatmap: terms with q below this")]
    pub q: f32,
}

/// The figures of the round `from`.
pub(crate) fn round_figures(from: &str) -> Result<CellFigures> {
    round_figures_of(&crate::manifest::run::load(from)?)
}

/// The figures of the round `loaded`.
pub(crate) fn round_figures_of(loaded: &crate::manifest::run::Loaded) -> Result<CellFigures> {
    let cells = crate::manifest::rounds::read_cells(loaded)?;
    let score = crate::manifest::rounds::read_score_cells(loaded, &cells)?;
    Ok(CellFigures::of_round(loaded, &cells, score))
}

/// `lupin figure round`: a round's charts and maps as files.
pub fn run_round_figures(args: &RoundFigureArgs) -> Result<()> {
    let mut figures = round_figures(&args.from)?;
    figures.set_term_cut(args.top, args.q)?;
    let d = figures.round_drawing(&args.charts, &args.maps, args.colour.as_deref(), &args.look)?;
    for f in figures.write_round(&d, &args.out)? {
        log::info!("wrote {}", f.display());
    }
    Ok(())
}
