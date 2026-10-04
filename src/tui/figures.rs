//! The figure pane of the order view: a trajectory run's figures drawn in the
//! terminal, exported as SVG + PDF, and logged in the gallery
//! (`docs/trajectory-plan.md` §6).

use super::gallery::{Gallery, Status};
use crate::manifest::run::load;
use crate::trajectory::figures::{self, Panel, TrajectoryData};
use anyhow::{Context, Result};
use clap::ValueEnum;
use image::DynamicImage;
use ratatui::layout::{Rect, Size};
use ratatui_image::picker::cap_parser::QueryStdioOptions;
use ratatui_image::picker::{Picker, ProtocolType};
use ratatui_image::protocol::Protocol;
use ratatui_image::{FilterType, Resize};
use std::cell::RefCell;
use std::path::{Path, PathBuf};
use std::time::Duration;

/// How figures reach the terminal, as `senna view` offers it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum Graphics {
    /// Ask the terminal; half-block characters when it does not answer.
    #[default]
    Auto,
    Kitty,
    Sixel,
    Iterm2,
    /// Unicode half-blocks, which every terminal shows.
    Blocks,
}

/// Figure width in inches and dots per inch for an export.
const EXPORT_WIDTH_IN: f32 = 7.0;
const EXPORT_DPI: f32 = 200.0;
/// Pixels per terminal cell assumed when sizing the raster for the pane.
const PX_PER_CELL: (u32, u32) = (10, 20);

/// A rendered figure, kept for the pane's size and panel.
struct Rendered {
    panel: Panel,
    size: (u16, u16),
    protocol: Protocol,
}

pub struct FiguresView {
    pub data: TrajectoryData,
    /// The manifest the figures come from.
    pub manifest: PathBuf,
    pub panels: Vec<Panel>,
    pub panel: usize,
    /// Show the figure instead of the order table.
    pub shown: bool,
    picker: Option<Picker>,
    rendered: RefCell<Option<Rendered>>,
    /// Drawing failed; shown in the pane.
    pub error: RefCell<Option<String>>,
    pub gallery: Gallery,
    pub gallery_open: bool,
    pub gallery_sel: usize,
    pub statuses: Vec<Status>,
    pub not_listed: Vec<PathBuf>,
}

impl FiguresView {
    /// The figures of the first of `manifests` with a `trajectory` section
    /// whose outputs can be read; `None` when there is none.
    pub fn load(manifests: &[&Path], graphics: Graphics) -> Option<Self> {
        let (manifest, data) = manifests.iter().find_map(|m| {
            let loaded = load(&m.to_string_lossy()).ok()?;
            TrajectoryData::load(&loaded)
                .ok()
                .flatten()
                .map(|d| (m.to_path_buf(), d))
        })?;
        let panels = data.panels();
        let picker = picker(graphics);
        let gallery = Gallery::here();
        let mut view = Self {
            data,
            manifest,
            panels,
            panel: 0,
            shown: false,
            picker,
            rendered: RefCell::new(None),
            error: RefCell::new(None),
            gallery,
            gallery_open: false,
            gallery_sel: 0,
            statuses: Vec::new(),
            not_listed: Vec::new(),
        };
        view.refresh();
        Some(view)
    }

    pub fn current(&self) -> Panel {
        self.panels[self.panel]
    }

    /// Show the figures, or move to the next panel when they are shown.
    pub fn next_panel(&mut self) {
        if self.shown {
            self.panel = (self.panel + 1) % self.panels.len();
        }
        self.shown = true;
    }

    /// Another pair of diffusion components, when that panel is shown.
    pub fn step_pair(&mut self, forward: bool) {
        if let Panel::Diffusion { x, y } = self.current() {
            let (nx, ny) = self.data.next_pair(x, y, forward);
            self.panels[self.panel] = Panel::Diffusion { x: nx, y: ny };
        }
    }

    /// The figure ready to draw into `area`, rendered on demand and kept
    /// while the panel and the pane's size stay the same.
    pub fn protocol(&self, area: Rect) -> Option<std::cell::Ref<'_, Protocol>> {
        let picker = self.picker.as_ref()?;
        let size = (area.width, area.height);
        let stale = self
            .rendered
            .borrow()
            .as_ref()
            .is_none_or(|r| r.panel != self.current() || r.size != size);
        if stale {
            let (w, h) = (
                u32::from(area.width) * PX_PER_CELL.0,
                u32::from(area.height) * PX_PER_CELL.1,
            );
            match self.render(picker, w.max(64), h.max(64), area) {
                Ok(protocol) => {
                    *self.rendered.borrow_mut() = Some(Rendered {
                        panel: self.current(),
                        size,
                        protocol,
                    });
                    *self.error.borrow_mut() = None;
                }
                Err(e) => {
                    *self.error.borrow_mut() = Some(format!("{e:#}"));
                    return None;
                }
            }
        }
        std::cell::Ref::filter_map(self.rendered.borrow(), |r| r.as_ref().map(|r| &r.protocol)).ok()
    }

    fn render(&self, picker: &Picker, w: u32, h: u32, area: Rect) -> Result<Protocol> {
        let svg = self.data.svg(self.current(), w, h)?;
        let img = figures::render(&svg, w, h)?;
        Ok(picker.new_protocol(
            DynamicImage::ImageRgba8(img),
            Size::new(area.width, area.height),
            Resize::Fit(Some(FilterType::Triangle)),
        )?)
    }

    /// Export the current panel as `{prefix}.trajectory.{panel}.svg` and
    /// `.pdf`, moved past existing files as `-2`, `-3` …, and log it.
    pub fn export(&mut self) -> Result<String> {
        let panel = self.current();
        let base = free_base(&format!("{}.trajectory.{}", self.data.prefix, panel.slug()));
        let (w, h) = (
            (EXPORT_WIDTH_IN * EXPORT_DPI) as u32,
            (EXPORT_WIDTH_IN * EXPORT_DPI * 0.75) as u32,
        );
        let svg = self.data.svg(panel, w, h)?;
        let files = figures::export(&svg, w, h, &base)?;
        let what = format!(
            "{} · {EXPORT_WIDTH_IN:.0} in · {EXPORT_DPI:.0} dpi",
            panel.title()
        );
        let logged = self
            .gallery
            .add(&files, &what, &panel.slug(), &self.manifest)
            .map(|()| String::new())
            .unwrap_or_else(|e| format!(" (not listed: {e:#})"));
        self.refresh();
        Ok(format!("wrote {}.svg and .pdf{logged}", base.display()))
    }

    /// Re-read the gallery and check every entry against its files.
    pub fn refresh(&mut self) {
        self.gallery = Gallery::open(Path::new(super::gallery::DIR));
        self.statuses = self.gallery.check();
        self.not_listed = self.gallery.not_listed(&self.data.prefix);
        self.gallery_sel = self.gallery_sel.min(self.rows().saturating_sub(1));
    }

    /// Rows of the gallery strip: the entries, then the files not listed.
    pub fn rows(&self) -> usize {
        self.gallery.entries.len() + self.not_listed.len()
    }

    /// Add the selected not-listed file set to the gallery.
    pub fn adopt_selected(&mut self) -> Result<String> {
        let i = self
            .gallery_sel
            .checked_sub(self.gallery.entries.len())
            .context("the selected row is already listed")?;
        let pdf = self.not_listed.get(i).context("nothing selected")?.clone();
        let mut files = vec![pdf.clone()];
        let svg = pdf.with_extension("svg");
        if svg.is_file() {
            files.push(svg);
        }
        let name = pdf
            .file_stem()
            .map_or(String::new(), |s| s.to_string_lossy().into_owned());
        let panel = name.rsplit(".trajectory.").next().unwrap_or("").to_string();
        self.gallery.add(&files, &name, &panel, &self.manifest)?;
        self.refresh();
        Ok(format!("listed {}", pdf.display()))
    }

    /// Take the selected entry off the list, deleting its files when `delete`.
    pub fn remove_selected(&mut self, delete: bool) -> Result<String> {
        let i = self.gallery_sel;
        anyhow::ensure!(
            i < self.gallery.entries.len(),
            "the selected row is not a listed export"
        );
        let name = self.gallery.entries[i].path.display().to_string();
        self.gallery.remove(i, delete)?;
        self.refresh();
        Ok(if delete {
            format!("deleted {name}")
        } else {
            format!("unlisted {name}")
        })
    }

    /// Move the selected entry's files to `new_base` (no extension).
    pub fn relocate_selected(&mut self, new_base: &str) -> Result<String> {
        let i = self.gallery_sel;
        anyhow::ensure!(
            i < self.gallery.entries.len(),
            "the selected row is not a listed export"
        );
        self.gallery.relocate(i, Path::new(new_base))?;
        self.refresh();
        Ok(format!("moved to {new_base}.svg and .pdf"))
    }

    /// The selected entry's base name, for a move prompt.
    pub fn selected_base(&self) -> Option<String> {
        let e = self.gallery.entries.get(self.gallery_sel)?;
        Some(e.path.with_extension("").to_string_lossy().into_owned())
    }
}

/// The terminal's picture protocol, asked for unless `graphics` names one;
/// half-blocks when the terminal does not answer.
fn picker(graphics: Graphics) -> Option<Picker> {
    let mut picker = if graphics == Graphics::Blocks {
        Picker::halfblocks()
    } else {
        Picker::from_query_stdio_with_options(QueryStdioOptions {
            timeout: Duration::from_millis(500),
            ..Default::default()
        })
        .unwrap_or_else(|_| Picker::halfblocks())
    };
    match graphics {
        Graphics::Auto | Graphics::Blocks => {}
        Graphics::Kitty => picker.set_protocol_type(ProtocolType::Kitty),
        Graphics::Sixel => picker.set_protocol_type(ProtocolType::Sixel),
        Graphics::Iterm2 => picker.set_protocol_type(ProtocolType::Iterm2),
    }
    Some(picker)
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
