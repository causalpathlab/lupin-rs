//! How `lupin figure` draws a map, as lupin's screen sets it: the labels'
//! size (`t`), the cluster lit (a click), a units round's GO stars, links
//! and rings (`g`, `l`, `r`), and the page (the save dialog's width and dpi).
//!
//! A figure is drawn for its page. On paper (a journal's page) it is the
//! page's width and as tall as it needs, and a chart wider than the page is
//! drawn along the page's length and turned onto it. On a slide (the figure
//! box of a 16:9 deck, as the `legume` report sets one: a title bar above,
//! a caption below) it is the box's height and as wide as it needs, and a
//! chart taller than the box is drawn across the slide and turned onto it.
//! Either way the text keeps its size on the page: 11 pt.

use super::turn::Turn;
use super::{cluster_name, CellFigures, Size, Style, EXPORT_DPI, EXPORT_WIDTH_IN};
use anyhow::Result;
use clap::ValueEnum;
use legume_tui::export::PageSize;

/// A figure's height over its width, as the screen's exports draw it.
pub(crate) const FIGURE_ASPECT: f32 = 3.0 / 4.0;

/// A page's height over its width, at most, before a figure wider than its
/// page is turned onto it: a journal page's.
pub(crate) const MAX_TALL: f32 = 1.3;

/// A slide's figure box, in inches: its height, and its widest (the 16:9
/// deck's text width).
pub(crate) const SLIDE_HEIGHT_IN: f32 = 3.5;
pub(crate) const SLIDE_WIDTH_IN: f32 = 12.2;

/// What a figure is set on.
#[derive(ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PageKind {
    /// A journal's page: fixed width, as tall as the figure needs.
    #[default]
    Paper,
    /// A 16:9 slide's figure box: fixed height, as wide as the figure needs.
    Slide,
}

/// The frame a figure is drawn in, and how it is turned onto its page, if
/// it is (drawn along the page's other side).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Fit {
    pub(crate) size: Size,
    pub(crate) turn: Option<Turn>,
}

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
    #[arg(
        long,
        value_enum,
        default_value_t,
        help = "What the figures are set on: paper (a journal page, 7 in wide, as tall as needed) or slide (a 16:9 deck's figure box, 3.5 in tall, as wide as needed); a chart too wide for paper, or too tall for a slide, is turned onto it"
    )]
    pub page: PageKind,
    #[arg(
        long,
        help = "Paper: the page's width in inches (7); slide: the widest a figure gets (12.2)"
    )]
    pub width: Option<f32>,
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
            page: PageKind::Paper,
            width: None,
            dpi: EXPORT_DPI as u32,
        }
    }
}

impl Look {
    /// Its page as legume-tui sizes one: the width asked for at its dpi.
    pub(crate) fn page(&self) -> PageSize {
        PageSize {
            width_in: self.width_in(),
            dpi: self.dpi,
        }
    }

    /// The pixels its text is a line of: 11 pt.
    pub(crate) fn text_px(&self) -> f32 {
        self.page().cell_px()
    }

    /// The width it asks for, in inches: paper's page, or a slide's widest.
    pub(crate) fn width_in(&self) -> f32 {
        self.width.unwrap_or(match self.page {
            PageKind::Paper => EXPORT_WIDTH_IN,
            PageKind::Slide => SLIDE_WIDTH_IN,
        })
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

    /// The pixels a figure that asks for no room of its own (a map) is
    /// drawn at: paper's width at [`FIGURE_ASPECT`], or a slide's height at
    /// its inverse.
    pub(crate) fn figure_size(&self) -> Result<Size> {
        Ok(self.fit((0, 0))?.size)
    }

    /// The frame for a figure whose content needs `natural` pixels across
    /// and down at this page's text size, and whether it is turned onto
    /// the page.
    ///
    /// Paper: the page's width, as tall as the content needs (at least
    /// [`FIGURE_ASPECT`]); content wider than the page and no deeper than
    /// the page is wide is drawn along the page's length (a full page, or
    /// longer) and turned left. Slide: the box's height, as wide as the content
    /// needs (from [`FIGURE_ASPECT`] to the slide's widest); content taller
    /// than the box and no wider than it is tall is drawn across the slide
    /// and turned right, so nothing on a screen reads upside down.
    pub(crate) fn fit(&self, (nw, nh): (usize, usize)) -> Result<Fit> {
        let (width, page) = (self.width_in(), self.page());
        let (w, short) = (page.check(FIGURE_ASPECT))
            .map_err(|e| anyhow::anyhow!("--width {width} --dpi {}: {e}", self.dpi))?;
        let px = |inches: f32| (inches * self.dpi as f32).round() as usize;
        let by = |n: usize, f: f32| (n as f32 * f).round() as usize;
        let (size, turn) = match self.page {
            PageKind::Paper => {
                if nw > w && nh <= w {
                    ((nw.max(by(w, MAX_TALL)), w), Some(Turn::Left))
                } else {
                    ((w, nh.max(short)), None)
                }
            }
            PageKind::Slide => {
                let h = px(SLIDE_HEIGHT_IN);
                let narrow = by(h, 1.0 / FIGURE_ASPECT);
                let wide = px(width).max(narrow);
                if nh > h && nw <= h {
                    ((h, nh.clamp(narrow, wide)), Some(Turn::Right))
                } else {
                    ((nw.clamp(narrow, wide), h), None)
                }
            }
        };
        let n = size.0.checked_mul(size.1);
        anyhow::ensure!(
            size.0 > 0 && size.1 > 0 && n.is_some_and(|n| n <= legume_tui::export::MAX_PIXELS),
            "--width {width} --dpi {}: a {} × {} px figure cannot be drawn; lower the dpi",
            self.dpi,
            size.0,
            size.1
        );
        Ok(Fit {
            size: Size {
                w: size.0 as u32,
                h: size.1 as u32,
                text_px: page.cell_px(),
            },
            turn,
        })
    }
}
