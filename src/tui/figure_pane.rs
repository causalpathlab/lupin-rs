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
use std::cell::{Ref, RefCell};
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
        let set_aside = self.gallery.set_aside.take();
        self.gallery = Gallery::here();
        self.gallery.set_aside = self.gallery.set_aside.take().or(set_aside);
        self.statuses = self.gallery.check();
        self.not_listed = self.gallery.not_listed(&self.prefix);
        self.sel = self.sel.min(self.rows().saturating_sub(1));
    }

    /// Rows of the strip: the entries, then the files not listed.
    pub fn rows(&self) -> usize {
        self.gallery.entries.len() + self.not_listed.len()
    }

    /// The selected listed export's PDF.
    fn selected_pdf(&self) -> Result<PathBuf> {
        self.gallery
            .entries
            .get(self.sel)
            .map(|e| e.path.clone())
            .context("the selected row is not a listed export")
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
        let pdf = self.selected_pdf()?;
        self.gallery.remove(&pdf, delete)?;
        let name = pdf.display();
        self.refresh();
        Ok(if delete {
            format!("deleted {name}")
        } else {
            format!("unlisted {name}")
        })
    }

    /// Move the export of `pdf` to `new_base` (no extension).
    pub fn relocate(&mut self, pdf: &Path, new_base: &str) -> Result<String> {
        self.gallery.relocate(pdf, Path::new(new_base))?;
        self.refresh();
        Ok(format!("moved to {new_base}.svg and .pdf"))
    }
}

/// A panel at a pane size, in terminal cells.
type Shown = (Panel, u16, u16);

pub struct FigurePane {
    pub data: TrajectoryData,
    pub panels: Vec<Panel>,
    pub sel: usize,
    /// Show the figure instead of the order table.
    pub shown: bool,
    picker: Picker,
    /// The figure on screen for (panel, pane size), or why it could not be
    /// drawn.
    rendered: RefCell<Option<(Shown, Result<Protocol, String>)>>,
    /// The current panel rasterised once; a new pane size only re-fits it.
    image: RefCell<Option<(Panel, Result<DynamicImage, String>)>>,
    pub exports: Exports,
}

impl FigurePane {
    /// The figures of the first of `manifests` with a `trajectory` section;
    /// `None` when there is none, an error when its outputs cannot be read.
    pub fn load(manifests: &[&Path], picker: &Picker) -> Result<Option<Self>> {
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
                    picker: picker.clone(),
                    image: RefCell::new(None),
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

    /// The current figure drawn for `area`, as senna view draws its
    /// figures. The panel is rasterised once (at the pane's pixel size when
    /// first shown) and only fitted again when the pane changes size; a
    /// failure is kept too, so it is not retried every frame.
    pub fn protocol(&self, area: Rect) -> Ref<'_, Result<Protocol, String>> {
        let key = (self.current(), area.width, area.height);
        if self
            .rendered
            .borrow()
            .as_ref()
            .is_none_or(|(k, _)| *k != key)
        {
            let drawn = self.image(key.0, area).and_then(|img| {
                self.picker
                    .new_protocol(
                        img,
                        Size::new(area.width, area.height),
                        Resize::Fit(Some(FilterType::Triangle)),
                    )
                    .map_err(|e| format!("{e:#}"))
            });
            *self.rendered.borrow_mut() = Some((key, drawn));
        }
        Ref::map(self.rendered.borrow(), |s| {
            &s.as_ref().expect("filled above").1
        })
    }

    /// `panel` rasterised, from the cache when it is the one last drawn.
    fn image(&self, panel: Panel, area: Rect) -> Result<DynamicImage, String> {
        if let Some((p, img)) = self.image.borrow().as_ref() {
            if *p == panel {
                return img.clone();
            }
        }
        let px = self.picker.font_size();
        let w = (u32::from(area.width) * u32::from(px.width)).max(64);
        let h = (u32::from(area.height) * u32::from(px.height)).max(64);
        let img = self
            .data
            .figure(panel, w, h)
            .and_then(|f| figures::render(&f))
            .map(DynamicImage::ImageRgba8)
            .map_err(|e| format!("{e:#}"));
        *self.image.borrow_mut() = Some((panel, img.clone()));
        img
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
pub fn picker(graphics: Graphics) -> Picker {
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
