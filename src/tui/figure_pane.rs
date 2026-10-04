//! The figure pane of the order view: a trajectory run's figures drawn in the
//! terminal, exported as SVG + PDF, and logged in the gallery
//! (`docs/trajectory-plan.md` §6).

use super::gallery::{Gallery, Status};
use crate::manifest::run::load;
use crate::trajectory::figures::{self, Panel, TrajectoryData, View};
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
    /// The export a first `D` asked to delete.
    delete_armed: Option<PathBuf>,
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
            delete_armed: None,
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
    /// A first `D` on the selected export: its name, to confirm with a
    /// second `D` on the same one; `None` when this is that second `D`.
    pub fn arm_delete(&mut self) -> Option<String> {
        let pdf = self.selected_pdf().ok()?;
        if self.delete_armed.take().as_ref() == Some(&pdf) {
            return None;
        }
        let name = pdf.with_extension("").display().to_string();
        self.delete_armed = Some(pdf);
        Some(name)
    }

    pub fn disarm(&mut self) {
        self.delete_armed = None;
    }

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

/// A thumbnail drawn at a tile size, in terminal cells, or why it could not
/// be.
type Tile = Option<((u16, u16), Result<Protocol, String>)>;

/// The thumbnail grid (`w`, as in `senna view`): every figure as a tile.
pub struct Grid {
    pub tiles: Vec<Panel>,
    pub sel: usize,
    drawn: RefCell<Vec<Tile>>,
}

impl Grid {
    /// Columns of the grid for its tiles: as square as it gets.
    pub fn cols(&self) -> usize {
        grid_cols(self.tiles.len())
    }

    /// Move the selection by an arrow key; `false` for any other key.
    pub fn step(&mut self, code: ratatui::crossterm::event::KeyCode) -> bool {
        use ratatui::crossterm::event::KeyCode;
        let (n, cols) = (self.tiles.len(), self.cols());
        self.sel = match code {
            KeyCode::Left => (self.sel + n - 1) % n,
            KeyCode::Right => (self.sel + 1) % n,
            KeyCode::Up if self.sel >= cols => self.sel - cols,
            KeyCode::Down if self.sel + cols < n => self.sel + cols,
            KeyCode::Up | KeyCode::Down => self.sel,
            _ => return false,
        };
        true
    }
}

/// Columns for `n` tiles: the smallest square that holds them.
pub fn grid_cols(n: usize) -> usize {
    (1..=n.max(1)).find(|c| c * c >= n).unwrap_or(1)
}

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
    /// The diffusion pair `m` comes back to.
    pair: (usize, usize),
    /// The thumbnail grid, while it is open.
    pub grid: Option<Grid>,
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
                return Ok(Some(Self::new(data, picker)));
            }
        }
        Ok(None)
    }

    fn new(data: TrajectoryData, picker: &Picker) -> Self {
        let exports = Exports::new(data.prefix());
        Self {
            panels: data.panels(),
            data,
            sel: 0,
            shown: false,
            picker: picker.clone(),
            image: RefCell::new(None),
            rendered: RefCell::new(None),
            pair: (1, 2),
            grid: None,
            exports,
        }
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
            self.set_view(View::default());
        }
    }

    /// Whether the panel on screen is a scatter, which pans and zooms.
    fn on_scatter(&self) -> bool {
        matches!(
            self.current(),
            Panel::Layout { .. } | Panel::Diffusion { .. }
        )
    }

    /// Zoom the scatter in or out by a step (`+` / `-`).
    pub fn zoom(&mut self, inward: bool) -> String {
        if !self.on_scatter() {
            return "only the scatter zooms (m and v choose it)".into();
        }
        let mut view = self.data.view;
        view.zoom(inward);
        self.set_view(view);
        if view.is_whole() {
            "the whole scatter".into()
        } else {
            format!("zoom ×{:.1} · arrows pan · 0 shows it whole", view.zoom)
        }
    }

    /// Pan the scatter by a step (the arrows); `dx`, `dy` right and up.
    pub fn pan(&mut self, dx: f32, dy: f32) -> String {
        if !self.on_scatter() {
            return "only the scatter pans (m and v choose it)".into();
        }
        let mut view = self.data.view;
        view.pan(dx, dy);
        if view == self.data.view {
            return if view.is_whole() {
                "+ zooms in; the arrows then pan".into()
            } else {
                "at the edge of the scatter".into()
            };
        }
        self.set_view(view);
        String::new()
    }

    /// The whole scatter again (`0`).
    pub fn reset_view(&mut self) -> String {
        if !self.on_scatter() {
            return "only the scatter zooms (m and v choose it)".into();
        }
        self.set_view(View::default());
        "the whole scatter".into()
    }

    fn set_view(&mut self, view: View) {
        if self.data.view != view {
            self.data.view = view;
            self.redraw();
        }
    }

    /// Open the thumbnail grid on the figure on screen (`w`).
    pub fn open_grid(&mut self) {
        let mut tiles = self.data.scatters();
        let pair = match self.current() {
            Panel::Diffusion { x, y } => (x, y),
            _ => self.pair,
        };
        for t in &mut tiles {
            if let Panel::Diffusion { .. } = t {
                *t = Panel::Diffusion {
                    x: pair.0,
                    y: pair.1,
                };
            }
        }
        tiles.extend(
            self.panels
                .iter()
                .copied()
                .filter(|p| !matches!(p, Panel::Layout { .. } | Panel::Diffusion { .. })),
        );
        let cur = self.current();
        let sel = tiles.iter().position(|&t| t == cur).unwrap_or(0);
        self.grid = Some(Grid {
            drawn: RefCell::new(vec![None; tiles.len()]),
            tiles,
            sel,
        });
    }

    /// Show the grid's selected tile and close the grid (Enter).
    pub fn open_tile(&mut self) {
        let Some(g) = self.grid.take() else { return };
        let Some(&tile) = g.tiles.get(g.sel) else {
            return;
        };
        let scatter = |p: &Panel| matches!(p, Panel::Layout { .. } | Panel::Diffusion { .. });
        let at = if scatter(&tile) {
            self.panels.iter().position(scatter)
        } else {
            self.panels.iter().position(|&p| p == tile)
        };
        let Some(at) = at else { return };
        if self.panels[at] != tile {
            if let Panel::Diffusion { x, y } = self.panels[at] {
                self.pair = (x, y);
            }
            self.panels[at] = tile;
            self.set_view(View::default());
        }
        self.sel = at;
        self.shown = true;
    }

    /// Tile `i` of the grid drawn for `area`, the scatter whole; kept until
    /// the tile's size or a setting changes.
    pub fn tile(&self, i: usize, area: Rect) -> Option<Ref<'_, Result<Protocol, String>>> {
        let g = self.grid.as_ref()?;
        let panel = *g.tiles.get(i)?;
        let size = (area.width, area.height);
        if g.drawn.borrow()[i].as_ref().is_none_or(|(s, _)| *s != size) {
            let px = self.picker.font_size();
            let w = (u32::from(area.width) * u32::from(px.width)).max(32);
            let h = (u32::from(area.height) * u32::from(px.height)).max(32);
            let drawn = self
                .data
                .figure_at(panel, w, h, View::default())
                .and_then(|f| figures::render(&f))
                .map_err(|e| format!("{e:#}"))
                .and_then(|img| {
                    self.picker
                        .new_protocol(
                            DynamicImage::ImageRgba8(img),
                            Size::new(area.width, area.height),
                            Resize::Fit(Some(FilterType::Triangle)),
                        )
                        .map_err(|e| format!("{e:#}"))
                });
            g.drawn.borrow_mut()[i] = Some((size, drawn));
        }
        Some(Ref::map(g.drawn.borrow(), |d| {
            &d[i].as_ref().expect("filled above").1
        }))
    }

    /// Labels on the scatters in turn (`t`), as `senna view` does.
    pub fn cycle_labels(&mut self) -> String {
        let note = self.data.cycle_labels();
        self.redraw();
        note
    }

    /// The next colouring of the scatters (`c`).
    pub fn cycle_colouring(&mut self) -> String {
        let note = self.data.cycle_colouring();
        self.redraw();
        note
    }

    /// The scatter's next coordinates (`m`): the run's layouts, then the
    /// diffusion map, which comes back on the pair last shown.
    pub fn next_layout(&mut self) -> String {
        let cur = self.current();
        let Some(next) = self.data.next_scatter(cur) else {
            return match cur {
                Panel::Layout { .. } | Panel::Diffusion { .. } => {
                    "only one layout in this run".into()
                }
                _ => "m switches the layout of the scatter (v shows it)".into(),
            };
        };
        if let Panel::Diffusion { x, y } = cur {
            self.pair = (x, y);
        }
        self.set_view(View::default());
        self.panels[self.sel] = match next {
            Panel::Diffusion { .. } => Panel::Diffusion {
                x: self.pair.0,
                y: self.pair.1,
            },
            p => p,
        };
        self.data.title(self.current())
    }

    /// Draw the current panel and the thumbnails again, after the style or
    /// the view changed.
    fn redraw(&self) {
        self.image.borrow_mut().take();
        self.rendered.borrow_mut().take();
        if let Some(g) = &self.grid {
            g.drawn.borrow_mut().iter_mut().for_each(|t| *t = None);
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
            .add(
                &e.files,
                &e.what,
                &self.data.slug(panel),
                &self.data.manifest,
            )
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

#[cfg(test)]
#[path = "tests/figure_pane.rs"]
mod tests;
