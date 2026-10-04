//! The figure pane of the order view: a trajectory run's figures drawn in the
//! terminal, exported as SVG + PDF, and logged in the gallery
//! (`docs/trajectory-plan.md` §6).

use super::gallery::{Gallery, Status};
use crate::manifest::run::load;
use crate::trajectory::figures::{self, Panel, TrajectoryData};
use anyhow::{Context, Result};
use clap::ValueEnum;
use image::DynamicImage;
use ratatui_image::picker::cap_parser::QueryStdioOptions;
use ratatui_image::picker::{Picker, ProtocolType};
use ratatui_image::protocol::StatefulProtocol;
use std::cell::{RefCell, RefMut};
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

/// The raster a figure is drawn at once; ratatui-image scales it to the pane.
const RASTER: (u32, u32) = (1200, 900);

/// The exports strip: the gallery, checked against the files, and the
/// figure files beside the run that it does not list.
pub struct Exports {
    /// The run's output prefix, whose figures are looked for.
    prefix: String,
    pub gallery: Gallery,
    pub open: bool,
    pub sel: usize,
    pub statuses: Vec<Status>,
    pub not_listed: Vec<PathBuf>,
}

impl Exports {
    fn new(prefix: String) -> Self {
        let mut e = Self {
            prefix,
            gallery: Gallery::here(),
            open: false,
            sel: 0,
            statuses: Vec::new(),
            not_listed: Vec::new(),
        };
        e.refresh();
        e
    }

    /// Re-read the gallery and check every entry against its files.
    pub fn refresh(&mut self) {
        self.gallery = Gallery::here();
        self.statuses = self.gallery.check();
        self.not_listed = self.gallery.not_listed(&self.prefix);
        self.sel = self.sel.min(self.rows().saturating_sub(1));
    }

    /// Rows of the strip: the entries, then the files not listed.
    pub fn rows(&self) -> usize {
        self.gallery.entries.len() + self.not_listed.len()
    }

    /// The selected row as an index into the gallery's entries.
    fn listed_index(&self) -> Result<usize> {
        anyhow::ensure!(
            self.sel < self.gallery.entries.len(),
            "the selected row is not a listed export"
        );
        Ok(self.sel)
    }

    /// Add the selected not-listed file set to the gallery, as a figure of
    /// `manifest`.
    pub fn adopt(&mut self, manifest: &Path) -> Result<String> {
        let i = self
            .sel
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
        let panel = name.rsplit(".trajectory.").next().unwrap_or("");
        self.gallery.add(&files, &name, panel, manifest)?;
        self.refresh();
        Ok(format!("listed {}", pdf.display()))
    }

    /// Take the selected entry off the list, deleting its files when `delete`.
    pub fn remove(&mut self, delete: bool) -> Result<String> {
        let i = self.listed_index()?;
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
    pub fn relocate(&mut self, new_base: &str) -> Result<String> {
        let i = self.listed_index()?;
        self.gallery.relocate(i, Path::new(new_base))?;
        self.refresh();
        Ok(format!("moved to {new_base}.svg and .pdf"))
    }

    /// The selected entry's base name, for a move prompt.
    pub fn selected_base(&self) -> Option<String> {
        let e = self.gallery.entries.get(self.sel)?;
        Some(e.path.with_extension("").to_string_lossy().into_owned())
    }
}

pub struct FigurePane {
    pub data: TrajectoryData,
    pub panels: Vec<Panel>,
    pub sel: usize,
    /// Show the figure instead of the order table.
    pub shown: bool,
    picker: Picker,
    /// The current panel, encoded once and resized to the pane by
    /// ratatui-image on each draw; or why it could not be drawn.
    rendered: RefCell<Option<(Panel, Result<StatefulProtocol, String>)>>,
    pub exports: Exports,
}

impl FigurePane {
    /// The figures of the first of `manifests` with a `trajectory` section;
    /// `None` when there is none, an error when its outputs cannot be read.
    pub fn load(manifests: &[&Path], graphics: Graphics) -> Result<Option<Self>> {
        for m in manifests {
            let Ok(loaded) = load(&m.to_string_lossy()) else {
                continue;
            };
            if let Some(data) = TrajectoryData::load(&loaded.manifest, &loaded.dir, &loaded.file)? {
                let exports = Exports::new(data.prefix());
                return Ok(Some(Self {
                    panels: data.panels(),
                    data,
                    sel: 0,
                    shown: false,
                    picker: picker(graphics),
                    rendered: RefCell::new(None),
                    exports,
                }));
            }
        }
        Ok(None)
    }

    pub fn current(&self) -> Panel {
        self.panels[self.sel]
    }

    /// Show the figures, or move to the next panel when they are shown.
    pub fn next_panel(&mut self) {
        if self.shown {
            self.sel = (self.sel + 1) % self.panels.len();
        }
        self.shown = true;
    }

    /// Another pair of diffusion components, when that panel is shown.
    pub fn step_pair(&mut self, forward: bool) {
        if let Panel::Diffusion { x, y } = self.current() {
            let (nx, ny) = self.data.next_pair(x, y, forward);
            self.panels[self.sel] = Panel::Diffusion { x: nx, y: ny };
        }
    }

    /// The current figure, ready to draw, rendered on first use.
    pub fn protocol(&self) -> RefMut<'_, Result<StatefulProtocol, String>> {
        let panel = self.current();
        let mut slot = self.rendered.borrow_mut();
        if slot.as_ref().is_none_or(|(p, _)| *p != panel) {
            *slot = Some((panel, self.render(panel).map_err(|e| format!("{e:#}"))));
        }
        RefMut::map(slot, |s| &mut s.as_mut().expect("filled above").1)
    }

    fn render(&self, panel: Panel) -> Result<StatefulProtocol> {
        let fig = self.data.figure(panel, RASTER.0, RASTER.1)?;
        let img = figures::render(&fig)?;
        Ok(self
            .picker
            .new_resize_protocol(DynamicImage::ImageRgba8(img)))
    }

    /// Export the current panel and log it.
    pub fn export(&mut self) -> Result<String> {
        let panel = self.current();
        let e = self.data.export(panel)?;
        let logged = self
            .exports
            .gallery
            .add(&e.files, &e.what, &panel.slug(), &self.data.manifest)
            .map(|()| String::new())
            .unwrap_or_else(|e| format!(" (not listed: {e:#})"));
        self.exports.refresh();
        Ok(format!("wrote {}.svg and .pdf{logged}", e.base.display()))
    }
}

/// The terminal's picture protocol, asked for unless `graphics` names one;
/// half-blocks when the terminal does not answer.
fn picker(graphics: Graphics) -> Picker {
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
    picker
}
