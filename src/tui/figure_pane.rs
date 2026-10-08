//! The figure pane of the order view: a trajectory run's figures drawn in the
//! terminal, exported as SVG + PDF, and logged in the gallery
//! (`docs/trajectory-plan.md` §6).

use super::gallery::{Gallery, Status};
use crate::manifest::run::load;
use crate::trajectory::figures::{self, Panel, Style, TrajectoryData, View};
use anyhow::{Context, Result};
use clap::ValueEnum;
use image::DynamicImage;
use ratatui::layout::{Rect, Size};
use ratatui_image::picker::{Picker, ProtocolType};
use ratatui_image::protocol::Protocol;
use ratatui_image::{FilterType, Resize};
use rayon::prelude::*;
use std::cell::{Ref, RefCell};
use std::path::{Path, PathBuf};

/// How figures reach the terminal, as `senna view` offers it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum Graphics {
    /// The terminal's own, as its environment names it (kitty, Ghostty,
    /// iTerm2, WezTerm, foot, mlterm); half-block characters otherwise.
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

/// What a figure on screen depends on: the panel and, for a scatter, the
/// style and the view (the other panels ignore both).
#[derive(Clone, Copy, Debug, PartialEq)]
struct Look {
    panel: Panel,
    style: Option<Style>,
    view: View,
}

/// A figure drawn for the terminal, or why it could not be.
type Drawn = Result<Protocol, String>;

/// A figure on screen: what it shows, at a size in terminal cells.
type Sized = (Look, (u16, u16));

/// Thumbnails kept across openings of the grid before they are dropped.
const TILES_KEPT: usize = 32;

/// The thumbnail grid (`w`, as in `senna view`): every figure as a tile.
pub struct Grid {
    pub tiles: Vec<Panel>,
    pub sel: usize,
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
    /// How the scatters are drawn (`t`, `c`).
    pub style: Style,
    /// The part of the scatter on screen (`+` `-` and the arrows).
    view: View,
    /// The diffusion pair the diffusion map shows (`,` `.`).
    pair: (usize, usize),
    picker: Picker,
    /// The figure on screen, drawn for the pane.
    rendered: RefCell<Option<(Sized, Drawn)>>,
    /// The figure on screen rasterised once; a new pane size only re-fits it.
    image: RefCell<Option<(Look, Result<DynamicImage, String>)>>,
    /// Thumbnails drawn, kept while the grid closes and opens again.
    tiles: RefCell<Vec<(Sized, Drawn)>>,
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
            style: Style::default(),
            view: View::default(),
            pair: (1, 2),
            picker: picker.clone(),
            image: RefCell::new(None),
            rendered: RefCell::new(None),
            tiles: RefCell::new(Vec::new()),
            grid: None,
            exports,
        }
    }

    /// The panel on screen; the diffusion map at the pair chosen.
    pub fn current(&self) -> Panel {
        self.at_pair(self.panels[self.sel])
    }

    /// `panel`, a diffusion map at the pair chosen.
    fn at_pair(&self, panel: Panel) -> Panel {
        match panel {
            Panel::Diffusion { .. } => Panel::Diffusion {
                x: self.pair.0,
                y: self.pair.1,
            },
            p => p,
        }
    }

    /// A panel's title, naming the layout and, on a scatter, the colouring.
    pub fn title(&self, panel: Panel) -> String {
        self.data.title(panel, self.style.colouring)
    }

    /// Show the figures, or move to the next panel when they are shown.
    pub fn next_panel(&mut self) {
        if self.shown {
            self.sel = (self.sel + 1) % self.panels.len();
            self.view = View::default();
        }
        self.shown = true;
    }

    /// Another pair of diffusion components, when that panel is shown.
    pub fn step_pair(&mut self, forward: bool, x_axis: bool) -> String {
        let Panel::Diffusion { x, y } = self.current() else {
            return "the axes step on the diffusion map (m reaches it)".into();
        };
        self.pair = self.data.next_pair(x, y, forward, x_axis);
        self.view = View::default();
        self.title(self.current())
    }

    /// Change the scatter's view by `f`; `None` when the panel on screen is
    /// not a scatter, else the view before and after.
    fn with_view(&mut self, f: impl FnOnce(&mut View)) -> Option<(View, View)> {
        if !self.current().is_scatter() {
            return None;
        }
        let before = self.view;
        f(&mut self.view);
        Some((before, self.view))
    }

    /// Zoom the scatter in or out by a step (`+` / `-`).
    pub fn zoom(&mut self, inward: bool) -> String {
        match self.with_view(|v| v.zoom(inward)) {
            None => "only the scatter zooms (m and v choose it)".into(),
            Some((_, v)) if v.is_whole() => "the whole scatter".into(),
            Some((_, v)) => format!("zoom ×{:.1} · arrows pan · 0 shows it whole", v.factor()),
        }
    }

    /// Pan the scatter by a step (the arrows); `dx`, `dy` right and up.
    pub fn pan(&mut self, dx: f32, dy: f32) -> String {
        match self.with_view(|v| v.pan(dx, dy)) {
            None => "only the scatter pans (m and v choose it)".into(),
            Some((a, b)) if a != b => String::new(),
            Some((_, v)) if v.is_whole() => "+ zooms in; the arrows then pan".into(),
            Some(_) => "at the edge of the scatter".into(),
        }
    }

    /// The whole scatter again (`0`).
    pub fn reset_view(&mut self) -> String {
        match self.with_view(|v| *v = View::default()) {
            None => "only the scatter zooms (m and v choose it)".into(),
            Some(_) => "the whole scatter".into(),
        }
    }

    /// What `panel` looks like drawn now, a scatter at `view`.
    fn look(&self, panel: Panel, view: View) -> Look {
        let scatter = panel.is_scatter();
        Look {
            panel,
            style: scatter.then_some(self.style),
            view: if scatter { view } else { View::default() },
        }
    }

    /// Open the thumbnail grid on the figure on screen (`w`).
    pub fn open_grid(&mut self) {
        let mut tiles: Vec<Panel> = self
            .data
            .scatters()
            .into_iter()
            .map(|p| self.at_pair(p))
            .collect();
        tiles.extend(self.panels.iter().copied().filter(|p| !p.is_scatter()));
        let cur = self.current();
        let sel = tiles.iter().position(|&t| t == cur).unwrap_or(0);
        self.grid = Some(Grid { tiles, sel });
    }

    /// Show the grid's selected tile and close the grid (Enter).
    pub fn open_tile(&mut self) {
        let Some(g) = self.grid.take() else { return };
        let Some(&tile) = g.tiles.get(g.sel) else {
            return;
        };
        let at = if tile.is_scatter() {
            self.panels.iter().position(|p| p.is_scatter())
        } else {
            self.panels.iter().position(|&p| p == tile)
        };
        let Some(at) = at else { return };
        self.panels[at] = tile;
        // A figure opened from the grid starts whole, as its tile shows it.
        self.view = View::default();
        self.sel = at;
        self.shown = true;
    }

    /// The grid's selected tile, while the grid is open.
    pub fn selected_tile(&self) -> Option<Panel> {
        let g = self.grid.as_ref()?;
        g.tiles.get(g.sel).copied()
    }

    /// Draw the grid's tiles `wanted` (index, area) that are not drawn yet:
    /// the figures in parallel, then their terminal pictures.
    pub fn draw_tiles(&self, wanted: &[(usize, Rect)]) {
        let Some(g) = &self.grid else { return };
        let missing: Vec<(Sized, Rect)> = wanted
            .iter()
            .filter_map(|&(i, area)| {
                let key = (
                    self.look(*g.tiles.get(i)?, View::default()),
                    (area.width, area.height),
                );
                let drawn = self.tiles.borrow().iter().any(|(k, _)| *k == key);
                (!drawn).then_some((key, area))
            })
            .collect();
        if missing.is_empty() {
            return;
        }
        let (data, px) = (&self.data, self.picker.font_size());
        let images: Vec<Result<image::RgbaImage, String>> = missing
            .par_iter()
            .map(|((look, _), area)| {
                let w = (u32::from(area.width) * u32::from(px.width)).max(32);
                let h = (u32::from(area.height) * u32::from(px.height)).max(32);
                let style = look.style.unwrap_or_default();
                data.figure(look.panel, w, h, &style, View::default())
                    .and_then(|f| figures::render(&f))
                    .map_err(|e| format!("{e:#}"))
            })
            .collect();
        let mut tiles = self.tiles.borrow_mut();
        if tiles.len() + missing.len() > TILES_KEPT {
            tiles.clear();
        }
        for ((key, area), img) in missing.into_iter().zip(images) {
            let drawn = img.and_then(|img| {
                self.picker
                    .new_protocol(
                        DynamicImage::ImageRgba8(img),
                        Size::new(area.width, area.height),
                        Resize::Fit(Some(FilterType::Triangle)),
                    )
                    .map_err(|e| format!("{e:#}"))
            });
            tiles.push((key, drawn));
        }
    }

    /// Tile `i` of the grid drawn for `area`, the scatter whole.
    pub fn tile(&self, i: usize, area: Rect) -> Option<Ref<'_, Drawn>> {
        let g = self.grid.as_ref()?;
        let key = (
            self.look(*g.tiles.get(i)?, View::default()),
            (area.width, area.height),
        );
        let find = || self.tiles.borrow().iter().position(|(k, _)| *k == key);
        if find().is_none() {
            self.draw_tiles(&[(i, area)]);
        }
        let at = find()?;
        Some(Ref::map(self.tiles.borrow(), |t| &t[at].1))
    }

    /// The scatter's next coordinates (`m`): the run's layouts, then the
    /// diffusion map at the pair last shown.
    pub fn next_layout(&mut self) -> String {
        let cur = self.current();
        let Some(next) = self.data.next_scatter(cur) else {
            return if cur.is_scatter() {
                match &self.data.layout_hint {
                    Some(h) => format!("only one layout in this run; {h}"),
                    None => "only one layout in this run".into(),
                }
            } else {
                "m switches the layout of the scatter (v shows it)".into()
            };
        };
        self.panels[self.sel] = next;
        self.view = View::default();
        self.title(self.current())
    }

    /// The current figure drawn for `area`, as senna view draws its
    /// figures. The figure is rasterised once (at the pane's pixel size when
    /// first shown) and only fitted again when the pane changes size; a
    /// failure is kept too, so it is not retried every frame.
    pub fn protocol(&self, area: Rect) -> Ref<'_, Drawn> {
        let look = self.look(self.current(), self.view);
        let key = (look, (area.width, area.height));
        if self
            .rendered
            .borrow()
            .as_ref()
            .is_none_or(|(k, _)| *k != key)
        {
            let drawn = self.image(look, area).and_then(|img| {
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

    /// `look` rasterised, from the cache when it is the one last drawn.
    fn image(&self, look: Look, area: Rect) -> Result<DynamicImage, String> {
        if let Some((l, img)) = self.image.borrow().as_ref() {
            if *l == look {
                return img.clone();
            }
        }
        let px = self.picker.font_size();
        let w = (u32::from(area.width) * u32::from(px.width)).max(64);
        let h = (u32::from(area.height) * u32::from(px.height)).max(64);
        let img = self
            .data
            .figure(look.panel, w, h, &self.style, look.view)
            .and_then(|f| figures::render(&f))
            .map(DynamicImage::ImageRgba8)
            .map_err(|e| format!("{e:#}"));
        *self.image.borrow_mut() = Some((look, img.clone()));
        img
    }

    /// Export the current panel as it is on screen and log it; with the
    /// grid open, its selected tile as the tile shows it (whole).
    pub fn export(&mut self) -> Result<String> {
        let (panel, view) = match self.selected_tile() {
            Some(tile) => (tile, View::default()),
            None => (self.current(), self.view),
        };
        let e = self.data.export(panel, &self.style, view)?;
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

/// The terminal's picture protocol: the one `graphics` names, else the one
/// the environment says the terminal has; half-blocks when it says none or
/// the cell's pixel size is unknown.
///
/// The terminal is not asked. ratatui-image's query reads stdin on a thread
/// that outlives its timeout and restores the terminal mode after it says it
/// is done, so beside the TUI's own key reader it eats keys and can leave
/// the terminal out of raw mode.
pub fn picker(graphics: Graphics) -> Picker {
    let cell = ratatui::crossterm::terminal::window_size()
        .ok()
        .and_then(|w| cell_pixels(&w));
    let protocol = match graphics {
        Graphics::Blocks => ProtocolType::Halfblocks,
        Graphics::Kitty => ProtocolType::Kitty,
        Graphics::Sixel => ProtocolType::Sixel,
        Graphics::Iterm2 => ProtocolType::Iterm2,
        Graphics::Auto if cell.is_none() => ProtocolType::Halfblocks,
        Graphics::Auto => protocol_from_env(|k| std::env::var(k).ok()),
    };
    if protocol == ProtocolType::Halfblocks {
        return Picker::halfblocks();
    }
    // A protocol named outright, with the cell's size unknown: a common one.
    #[allow(deprecated)] // in favour of the stdin query, which is what is avoided
    let mut picker = Picker::from_fontsize(cell.unwrap_or((10, 20)).into());
    picker.set_protocol_type(protocol);
    picker
}

/// The picture protocol the terminal's environment variables (`var`) name.
pub(crate) fn protocol_from_env(var: impl Fn(&str) -> Option<String>) -> ProtocolType {
    let v = |k: &str| var(k).unwrap_or_default();
    let (term, program) = (v("TERM"), v("TERM_PROGRAM"));
    // Through tmux a picture needs passthrough, which `--graphics` can ask for.
    if var("TMUX").is_some() {
        ProtocolType::Halfblocks
    } else if var("KITTY_WINDOW_ID").is_some()
        || term.contains("kitty")
        || term.contains("ghostty")
        || program.eq_ignore_ascii_case("ghostty")
    {
        ProtocolType::Kitty
    } else if program == "iTerm.app" || program == "WezTerm" || v("LC_TERMINAL") == "iTerm2" {
        ProtocolType::Iterm2
    } else if term.starts_with("foot") || term.starts_with("mlterm") {
        ProtocolType::Sixel
    } else {
        ProtocolType::Halfblocks
    }
}

/// A cell's size in pixels, from the window's when the terminal reports it.
pub(crate) fn cell_pixels(w: &ratatui::crossterm::terminal::WindowSize) -> Option<(u16, u16)> {
    (w.columns > 0 && w.rows > 0 && w.width > 0 && w.height > 0)
        .then(|| (w.width / w.columns, w.height / w.rows))
}

#[cfg(test)]
#[path = "tests/figure_pane.rs"]
pub(crate) mod tests;
