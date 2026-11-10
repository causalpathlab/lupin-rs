//! How `lupin figure` draws a map, as lupin's screen sets it: the labels'
//! size (`t`), the cluster lit (a click), a units round's GO stars, links
//! and rings (`g`, `l`, `r`), and the page (the save dialog's width and dpi).

use super::{cluster_name, CellFigures, Size, Style, EXPORT_DPI, EXPORT_WIDTH_IN};
use anyhow::Result;
use clap::ValueEnum;
use legume_tui::export::PageSize;

/// A figure's height over its width, as the screen's exports draw it.
pub(crate) const FIGURE_ASPECT: f32 = 3.0 / 4.0;

/// The type labels' size on a map.
#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum LabelSize {
    Small,
    Medium,
    Large,
    Largest,
    Off,
}

/// What a units round's map draws over its cells.
#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum Overlay {
    /// The lit cluster's GO terms as numbered stars.
    Go,
    /// The strongest links between the types.
    Links,
    /// Each type's fitted radius as a ring.
    Rings,
    /// None of them.
    None,
}

#[derive(clap::Args, Clone, Debug, PartialEq)]
pub struct Look {
    #[arg(
        long,
        value_enum,
        help = "The type labels on a map: small, medium (default), large, largest or off"
    )]
    pub labels: Option<LabelSize>,
    #[arg(
        long,
        value_name = "CLUSTER",
        help = "Light this cluster's cells, the others grey (on a units round's map, with its GO terms)"
    )]
    pub focus: Option<i32>,
    #[arg(
        long = "overlay",
        value_enum,
        value_delimiter = ',',
        help = "A units round's map: go, links, rings or none (default go,links)"
    )]
    pub overlays: Vec<Overlay>,
    #[arg(long, default_value_t = EXPORT_WIDTH_IN, help = "Page width in inches; the height follows")]
    pub width: f32,
    #[arg(
        long,
        default_value_t = EXPORT_DPI as u32,
        help = "Dots per inch of the cells (a PDF's text stays text)"
    )]
    pub dpi: u32,
}

impl Default for Look {
    fn default() -> Self {
        Self {
            labels: None,
            focus: None,
            overlays: Vec::new(),
            width: EXPORT_WIDTH_IN,
            dpi: EXPORT_DPI as u32,
        }
    }
}

impl Look {
    /// The page it asks for.
    pub(crate) fn page(&self) -> PageSize {
        PageSize {
            width_in: self.width,
            dpi: self.dpi,
        }
    }

    /// `style` as it asks, for `data`: a cluster lit must be one of its
    /// clusters, and overlays are a units round's.
    pub(crate) fn apply(&self, style: &mut Style, data: &CellFigures) -> Result<()> {
        if let Some(size) = self.labels {
            style.labels = match size {
                LabelSize::Small => Some(0),
                LabelSize::Medium => Some(1),
                LabelSize::Large => Some(2),
                LabelSize::Largest => Some(3),
                LabelSize::Off => None,
            };
        }
        if let Some(k) = self.focus {
            if !data.cluster.contains(&k) {
                let mut ks: Vec<i32> = data.cluster.clone();
                ks.sort_unstable();
                ks.dedup();
                let ks: Vec<String> = ks.into_iter().map(cluster_name).collect();
                anyhow::bail!(
                    "no cluster {} to light: {}",
                    cluster_name(k),
                    if ks.is_empty() {
                        "the cells have no clusters".into()
                    } else {
                        format!("the clusters are {}", ks.join(", "))
                    }
                );
            }
            style.focus = Some(k);
        }
        if !self.overlays.is_empty() {
            anyhow::ensure!(
                data.units.is_some(),
                "go, links and rings are drawn on a units round's map; this one has no units"
            );
            let on = |o| self.overlays.contains(&o);
            (style.go, style.links, style.rings) =
                (on(Overlay::Go), on(Overlay::Links), on(Overlay::Rings));
        }
        self.figure_size()?;
        Ok(())
    }

    /// The pixels a figure is drawn at ([`FIGURE_ASPECT`]).
    pub(crate) fn figure_size(&self) -> Result<Size> {
        self.size(FIGURE_ASPECT)
    }

    /// The pixels a figure `aspect` (height over width) tall is drawn at,
    /// its text 11 pt; or why the page cannot be drawn.
    pub(crate) fn size(&self, aspect: f32) -> Result<Size> {
        let page = self.page();
        let (w, h) = (page.check(aspect))
            .map_err(|e| anyhow::anyhow!("--width {} --dpi {}: {e}", self.width, self.dpi))?;
        Ok(Size {
            w: w as u32,
            h: h as u32,
            text_px: page.cell_px(),
        })
    }
}
