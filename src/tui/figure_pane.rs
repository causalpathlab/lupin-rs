//! The figure pane of the order view: a trajectory run's figures drawn in the
//! terminal, a scatter as senna draws its map, exported (a scatter as a PDF
//! through senna's save dialog, the other panels at once as SVG + PDF) and
//! logged in the gallery with a thumbnail (`docs/trajectory-plan.md` §6). The round's cell map is one too, with the
//! mouse: a drag pans, the wheel zooms, a click picks the cluster under it,
//! and a label pressed is dragged. Labels moved (by the mouse or `E`) are kept
//! beside the run in lupin's own `{prefix}.lupin_labels.json`. On a round,
//! `H` steps from the map to its structure plot and heatmap and back, as
//! senna steps through its charts, each saved through the same dialog.

use super::app::step;
use super::gallery::{Gallery, Status, Trace};
use crate::figures::charts::{HeatKind, Mixture};
use crate::figures::{
    self, CellFigures, CellMap, Colouring, Panel, Picture, Size, Source, Style, View,
};
use crate::manifest::run::load;
use anyhow::{Context, Result};
use image::RgbaImage;
use legume_tui::export::{PageSize, PdfArgs};
use legume_tui::look::Look as Chrome;
use legume_tui::map::display::fit_protocol;
use legume_tui::map::labels::{Hits, Placed};
use legume_tui::map::offsets::{OffsetBook, Offsets};
use legume_tui::map::pointer::{EditKey, LabelEdit, MapAction, MapArea, Pointer, PointerOpts};
use legume_tui::map::viewport::Viewport;
use legume_tui::save::{Numbering, SaveDialog, SaveKey, SaveRequest, SaveTarget};
use ratatui::crossterm::event::{KeyCode, KeyEvent, MouseEvent};
use ratatui::layout::Rect;
use ratatui_image::picker::Picker;
use ratatui_image::protocol::Protocol;
use rayon::prelude::*;
use std::cell::{Cell, Ref, RefCell};
use std::path::{Path, PathBuf};
use std::sync::mpsc::{channel, Receiver, TryRecvError};
use std::sync::Arc;

/// The exports strip: the gallery, checked against the files, and the
/// figure files beside the run that it does not list.
pub struct Exports {
    /// The run's output prefix and the figures' kind (`run.trajectory`),
    /// whose figures are looked for.
    prefix: String,
    pub gallery: Gallery,
    pub open: bool,
    pub sel: usize,
    pub statuses: Vec<Status>,
    pub not_listed: Vec<PathBuf>,
    /// The export a first `D` asked to delete.
    delete_armed: Option<PathBuf>,
    /// The pane, as last drawn, has room for the strip's selected row.
    pub fits: std::cell::Cell<bool>,
}

impl Exports {
    fn new(prefix: String) -> Self {
        let mut e = Self {
            prefix,
            gallery: super::gallery::here(),
            open: false,
            sel: 0,
            statuses: Vec::new(),
            not_listed: Vec::new(),
            delete_armed: None,
            fits: std::cell::Cell::new(true),
        };
        e.refresh();
        e
    }

    /// Re-read the gallery and check every entry against its files.
    pub fn refresh(&mut self) {
        self.gallery.reload();
        self.statuses = self.gallery.check();
        self.not_listed = super::gallery::not_listed(&self.gallery, &self.prefix);
        self.sel = self.sel.min(self.rows().saturating_sub(1));
    }

    /// Rows of the strip: the entries, then the files not listed.
    pub fn rows(&self) -> usize {
        self.gallery.entries().len() + self.not_listed.len()
    }

    /// The selected listed export's PDF.
    fn selected_pdf(&self) -> Result<PathBuf> {
        self.gallery
            .entries()
            .get(self.sel)
            .map(|e| e.path.clone())
            .context("the selected row is not a listed export")
    }

    /// Add the selected not-listed file set to the gallery, as a figure of
    /// `manifest`.
    pub fn adopt(&mut self, manifest: &Path) -> Result<String> {
        let i = self
            .sel
            .checked_sub(self.gallery.entries().len())
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
        let stem = super::gallery::figure_stem(&self.prefix);
        let panel = name.strip_prefix(&stem).unwrap_or(&name);
        self.gallery
            .add(&files, &name, None, Trace::of(panel, manifest))?;
        self.refresh();
        // The list is newest first: the selection follows the file.
        let listed = crate::manifest::data_files::absolute(&pdf);
        if let Some(i) = self.gallery.entries().iter().position(|e| e.path == listed) {
            self.sel = i;
        }
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
        Ok(format!("moved to {new_base}"))
    }
}

/// What a figure on screen depends on: the panel and, for a scatter, the
/// style, the view and the labels moved so far (the other panels ignore
/// them).
#[derive(Clone, Copy, Debug, PartialEq)]
struct Look {
    panel: Panel,
    style: Option<Style>,
    view: View,
    moves: u32,
}

/// A figure drawn for the terminal, or why it could not be.
type Drawn = Result<Protocol, String>;

/// A figure on screen: what it shows, at a size in terminal cells.
type Sized = (Look, (u16, u16));

/// A figure rasterised, with where its labels went when it is a scatter.
type Raster = Result<(RgbaImage, Option<Arc<OnMap>>), String>;

/// A look being rasterised, and where its pixels will come.
type Rasterising = (Look, Receiver<Raster>);

/// A scatter as drawn: its camera and its labels, which the pointer and
/// label editing find.
pub(crate) struct OnMap {
    pub(crate) panel: Panel,
    pub(crate) vp: Viewport,
    pub(crate) placed: Vec<Placed<u32>>,
    pub(crate) names: Vec<Box<str>>,
    /// The labels' line height in pixels: an offset's unit.
    pub(crate) line: f32,
    hits: Hits<u32>,
}

impl OnMap {
    fn of(panel: Panel, m: &CellMap) -> Self {
        let mut hits = Hits::default();
        hits.record(m.vp, &m.placed);
        Self {
            panel,
            vp: m.vp,
            placed: m.placed.clone(),
            names: m.names.clone(),
            line: m.line,
            hits,
        }
    }

    /// The moves `moved` keeps by name, by the labels' places here.
    fn by_place(&self, moved: &Offsets<String>) -> Offsets<u32> {
        self.names
            .iter()
            .enumerate()
            .filter_map(|(k, n)| Some((k as u32, moved.get(&n.to_string())?)))
            .collect()
    }
}

/// What the pointer did on the map.
#[derive(Debug, PartialEq, Eq)]
pub enum MapMouse {
    /// Nothing to show.
    None,
    /// The map or a label moved: draw it again.
    Redraw,
    /// The pointer is over this cluster's cells, or no cluster's.
    Hover(Option<i32>),
    /// A click on this cluster's cells, or on no cluster's.
    Click(Option<i32>),
    /// Taken; the status line says this.
    Say(String),
}

/// The file beside the run that keeps its moved labels, lupin's own.
const LABELS_SUFFIX: &str = "lupin_labels.json";

/// Why a figure rasterised off the screen's thread never came.
const NOT_DRAWN: &str = "the figure was not drawn";

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

/// What a key did on the figure pane.
#[derive(Debug, PartialEq, Eq)]
pub enum FigureKey {
    /// Not one of its keys.
    Ignored,
    /// Taken, with nothing to say.
    Took,
    /// Taken; the status line says this.
    Say(String),
    /// Ask where to move the export of this PDF.
    Relocate(PathBuf),
}

/// What the thumbnail grid's keys do.
pub const GRID_KEYS: &str =
    "arrows choose · Enter, space or 1-9 opens · p exports it · t c restyle · V leaves the figures · esc or w closes";

/// Why the exports strip's keys do nothing: the pane is too short for it.
const STRIP_NO_ROOM: &str =
    "the exports strip has no room in this pane: a taller terminal shows it (f closes it)";

pub struct FigurePane {
    /// Shared with the thread rasterising a figure off the screen's.
    pub data: Arc<CellFigures>,
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
    image: RefCell<Option<(Look, Raster)>>,
    /// A figure being rasterised off the screen's thread, for a quick frame.
    pending: RefCell<Option<Rasterising>>,
    /// Thumbnails drawn, kept while the grid closes and opens again.
    tiles: RefCell<Vec<(Sized, Drawn)>>,
    /// The thumbnail grid, while it is open.
    pub grid: Option<Grid>,
    pub exports: Exports,
    /// The labels moved, beside the run; why not when the file there cannot
    /// be read (it is then left as it is, and the labels where they are).
    book: Result<OffsetBook, String>,
    /// Label moves so far, so a figure is drawn again after each.
    moves: u32,
    /// A label dragged since the book was last written.
    unsaved: bool,
    pointer: Pointer<u32>,
    /// Moving labels from the keys (`E`), while it is on.
    edit: Option<LabelEdit<u32>>,
    /// The scatter on screen, as drawn, and where on the screen.
    pub(crate) on_screen: RefCell<Option<Arc<OnMap>>>,
    pub(crate) area: Cell<Option<Rect>>,
    /// The save dialog, while it is open.
    save: Option<Saving>,
    /// A round's chart shown in place of the panel (`H`).
    chart: Option<Panel>,
    /// What the structure plot and the heatmap show when `H` reaches them
    /// (`c` on them).
    mixture: Mixture,
    heat: HeatKind,
}

impl FigurePane {
    /// The figures of the first of `manifests` with a `trajectory` section;
    /// `None` when there is none, an error when its outputs cannot be read.
    pub fn load(manifests: &[&Path], picker: &Picker) -> Result<Option<Self>> {
        for m in manifests {
            let Ok(loaded) = load(&m.to_string_lossy()) else {
                continue;
            };
            if let Some(data) = CellFigures::load(&loaded.manifest, &loaded.dir, &loaded.file)? {
                return Ok(Some(Self::new(data, picker)));
            }
        }
        Ok(None)
    }

    fn new(data: CellFigures, picker: &Picker) -> Self {
        let exports = Exports::new(format!("{}.{}", data.prefix(), data.kind()));
        let style = Style::opening(&data);
        let path = labels_file(&data);
        let book = OffsetBook::load(&path).map_err(|e| {
            format!(
                "the moved labels in {} cannot be read ({e}): the labels stay where they are and the file is left as it is",
                path.display()
            )
        });
        Self {
            panels: data.panels(),
            data: Arc::new(data),
            sel: 0,
            shown: false,
            style,
            view: View::default(),
            pair: (1, 2),
            picker: picker.clone(),
            image: RefCell::new(None),
            pending: RefCell::new(None),
            rendered: RefCell::new(None),
            tiles: RefCell::new(Vec::new()),
            grid: None,
            exports,
            book,
            moves: 0,
            unsaved: false,
            // One wheel step is one of the view's zoom steps.
            pointer: Pointer::new(PointerOpts {
                wheel: std::f32::consts::SQRT_2,
                ..PointerOpts::default()
            }),
            edit: None,
            on_screen: RefCell::new(None),
            area: Cell::new(None),
            save: None,
            chart: None,
            mixture: Mixture::Types,
            heat: HeatKind::Markers,
        }
    }

    /// Draw as `before` did: its labels, colouring (when this data has it),
    /// view and exports strip; a new round of the same run looks the same.
    pub fn keep_look(&mut self, before: &FigurePane) {
        let colouring = before.style.colouring;
        self.style.labels = before.style.labels;
        (self.style.go, self.style.links, self.style.rings) =
            (before.style.go, before.style.links, before.style.rings);
        if self.data.colourings().contains(&colouring) {
            self.style.colouring = colouring;
        }
        if before.data.layouts.len() == self.data.layouts.len() {
            self.panels = before.panels.clone();
            self.sel = before.sel.min(self.panels.len().saturating_sub(1));
        }
        self.view = before.view;
        self.exports.open = before.exports.open;
        (self.chart, self.mixture, self.heat) = (before.chart, before.mixture, before.heat);
    }

    /// Draw for the terminal `picker` is for: every terminal picture made
    /// again, and with `changed` (another protocol or cell size) every
    /// figure rasterised again too.
    pub fn repick(&mut self, picker: &Picker, changed: bool) {
        self.picker = picker.clone();
        *self.rendered.get_mut() = None;
        *self.on_screen.get_mut() = None;
        self.tiles.get_mut().clear();
        if changed {
            *self.image.get_mut() = None;
            *self.pending.get_mut() = None;
        }
    }

    /// `data`'s figures, shown at once.
    pub fn showing(data: CellFigures, picker: &Picker) -> Self {
        let mut p = Self::new(data, picker);
        p.shown = true;
        p
    }

    /// The pane's keys and its exports strip's.
    pub fn key(&mut self, code: KeyCode) -> FigureKey {
        use FigureKey::{Ignored, Relocate, Say, Took};
        if code != KeyCode::Char('D') {
            self.exports.disarm();
        }
        // The thumbnail grid is modal: arrows choose, Enter opens, p exports
        // the tile chosen, and only the app's own keys (quit, run, stop,
        // save, export the round, panes, the tree) go by.
        if let Some(g) = &mut self.grid {
            if g.step(code) {
                return Took;
            }
            // As in senna view, space or a tile's number opens it too.
            if let KeyCode::Char(c @ '1'..='9') = code {
                let k = c as usize - '1' as usize;
                if k >= g.tiles.len() {
                    return Say(format!("{} tile(s)", g.tiles.len()));
                }
                g.sel = k;
            }
            return match code {
                KeyCode::Enter | KeyCode::Char(' ' | '1'..='9') => {
                    self.open_tile();
                    Say(self.title(self.current()))
                }
                KeyCode::Esc | KeyCode::Char('w') => {
                    self.grid = None;
                    Say(self.title(self.current()))
                }
                KeyCode::Char('V') => {
                    self.grid = None;
                    self.shown = false;
                    Took
                }
                KeyCode::Char('t') => Say(self.style.cycle_labels()),
                KeyCode::Char('c') => Say(self.style.cycle_colouring(&self.data.colourings())),
                KeyCode::Char('p') => match self.selected_tile() {
                    Some(tile) if tile.is_scatter() => self.ask_save(tile, true),
                    _ => Say(self
                        .export()
                        .unwrap_or_else(|e| format!("export failed: {e:#}"))),
                },
                // The strip is drawn under a figure, never over the grid.
                KeyCode::Char('f') => {
                    Say("f lists the exports under a figure: enter opens the tile".into())
                }
                KeyCode::Char(c) if !"qrxseT".contains(c) => Say(GRID_KEYS.into()),
                _ => Ignored,
            };
        }
        // The strip's keys act only while it is on screen; the arrows it
        // does not take pan. Open in a pane too short for it, they say so.
        let open = self.shown && self.exports.open;
        let strip = open && self.exports.fits.get();
        if open && !strip && matches!(code, KeyCode::Enter | KeyCode::Char('d' | 'D' | 'M')) {
            return Say(STRIP_NO_ROOM.into());
        }
        if strip {
            let rows = self.exports.rows();
            if step(&mut self.exports.sel, rows, code) {
                return Took;
            }
        }
        let shown = self.shown;
        match code {
            // Esc closes the exports strip before it leaves the pane.
            KeyCode::Esc if strip => {
                self.exports.open = false;
                Took
            }
            // A round's charts, as senna's: H steps, Esc the map again.
            KeyCode::Char('H') if shown => Say(self.step_chart()),
            KeyCode::Esc if self.chart.is_some() => {
                self.chart = None;
                Say("the map · H its charts".into())
            }
            KeyCode::Char('c') if self.chart.is_some() => Say(self.step_shown()),
            KeyCode::Char('p') if self.chart.is_some() => self.ask_save(self.current(), false),
            // The strip and its check belong to the figures.
            KeyCode::Char('f' | 'R') if !shown => {
                Say("v shows the figures, then f lists their exports".into())
            }
            KeyCode::Char('v') => {
                self.next_panel();
                Took
            }
            KeyCode::Char('V') => {
                self.shown = false;
                self.grid = None;
                Took
            }
            KeyCode::Char(',') if shown => Say(self.step_pair(false, false)),
            KeyCode::Char('.') if shown => Say(self.step_pair(true, false)),
            KeyCode::Char('[') if shown => Say(self.step_pair(false, true)),
            KeyCode::Char(']') if shown => Say(self.step_pair(true, true)),
            // senna view's keys: labels, colouring, layout.
            // Labels and colours are the scatter's.
            KeyCode::Char('t' | 'c') if shown && !self.current().is_scatter() => {
                Say("t and c restyle the scatter (v reaches it)".into())
            }
            KeyCode::Char('t') if shown => Say(self.style.cycle_labels()),
            KeyCode::Char('c') if shown => Say(self.style.cycle_colouring(&self.data.colourings())),
            KeyCode::Char('m') if shown => Say(self.next_layout()),
            // A units round's GO stars, links and rings.
            KeyCode::Char(c @ ('g' | 'l' | 'r')) if self.units_keys() => {
                Say(self.style.toggle_units(c))
            }
            // senna view's navigation: the grid, zoom and pan.
            KeyCode::Char('w') if shown => {
                self.open_grid();
                Say(GRID_KEYS.into())
            }
            KeyCode::Char('+' | '=') if shown => Say(self.zoom(true)),
            KeyCode::Char('-' | '_') if shown => Say(self.zoom(false)),
            KeyCode::Char('0') if shown => Say(self.reset_view()),
            KeyCode::Left | KeyCode::Right | KeyCode::Up | KeyCode::Down if shown => {
                let (dx, dy) = match code {
                    KeyCode::Left => (-1.0, 0.0),
                    KeyCode::Right => (1.0, 0.0),
                    KeyCode::Up => (0.0, 1.0),
                    _ => (0.0, -1.0),
                };
                let note = self.pan(dx, dy);
                if note.is_empty() {
                    Took
                } else {
                    Say(note)
                }
            }
            KeyCode::Char('p') if shown && self.current().is_scatter() => {
                self.ask_save(self.current(), false)
            }
            KeyCode::Char('p') if shown => Say(self
                .export()
                .unwrap_or_else(|e| format!("export failed: {e:#}"))),
            KeyCode::Char('f') => {
                self.exports.open = !self.exports.open;
                if self.exports.open && !self.exports.fits.get() {
                    Say(STRIP_NO_ROOM.into())
                } else {
                    Took
                }
            }
            KeyCode::Char('R') => {
                let x = &mut self.exports;
                x.refresh();
                Say(format!(
                    "{} export(s) checked, {} not listed",
                    x.gallery.entries().len(),
                    x.not_listed.len()
                ))
            }
            KeyCode::Enter if strip => Say(self
                .exports
                .adopt(&self.data.manifest)
                .unwrap_or_else(|e| format!("{e:#}"))),
            KeyCode::Char('d') if strip => Say(self
                .exports
                .remove(false)
                .unwrap_or_else(|e| format!("{e:#}"))),
            KeyCode::Char('D') if strip => {
                let x = &mut self.exports;
                Say(match x.arm_delete() {
                    Some(name) => format!("D again deletes {name} and its files"),
                    None => x.remove(true).unwrap_or_else(|e| format!("{e:#}")),
                })
            }
            KeyCode::Char('M') if strip => {
                match self.exports.gallery.entries().get(self.exports.sel) {
                    Some(e) => Relocate(e.path.clone()),
                    None => Took,
                }
            }
            _ => Ignored,
        }
    }

    /// Whether `g`, `l` and `r` are the units' here: a units round's map
    /// on screen, not its grid.
    pub fn units_keys(&self) -> bool {
        self.shown
            && self.grid.is_none()
            && self.data.units.is_some()
            && self.current().is_scatter()
    }

    /// The panel on screen: a chart, else the panel chosen, the diffusion
    /// map at the pair chosen.
    pub fn current(&self) -> Panel {
        self.chart
            .unwrap_or_else(|| self.at_pair(self.panels[self.sel]))
    }

    /// Whether a chart is shown in place of the map.
    pub fn charting(&self) -> bool {
        self.chart.is_some()
    }

    /// `H`: the map, then the structure plot, then the heatmap, then the
    /// map again; what is shown.
    fn step_chart(&mut self) -> String {
        if self.data.source != Source::Round {
            return "charts are a round's: H on its cell map".into();
        }
        self.chart = match self.chart {
            None => {
                // Topics when the cell types are not there to stack.
                let other = match self.mixture {
                    Mixture::Types => Mixture::Topics,
                    Mixture::Topics => Mixture::Types,
                };
                if !self.data.has_mixture(self.mixture) && self.data.has_mixture(other) {
                    self.mixture = other;
                }
                Some(Panel::Structure { of: self.mixture })
            }
            Some(Panel::Structure { .. }) => Some(Panel::Heatmap { of: self.heat }),
            Some(_) => None,
        };
        match self.chart {
            Some(c) => self.chart_note(c),
            None => "the map · H its charts".into(),
        }
    }

    /// `c` on a chart: the other mixture, or the next heatmap.
    fn step_shown(&mut self) -> String {
        match self.chart {
            Some(Panel::Structure { of }) => {
                let other = match of {
                    Mixture::Types => Mixture::Topics,
                    Mixture::Topics => Mixture::Types,
                };
                if let Some(why) = self.data.chart_missing(Panel::Structure { of: other }) {
                    return format!("{} only: no {} ({why})", of.title(), other.title());
                }
                self.mixture = other;
            }
            Some(Panel::Heatmap { of }) => {
                let all = HeatKind::ALL;
                let at = all.iter().position(|&h| h == of).unwrap_or(0);
                self.heat = all[(at + 1) % all.len()];
            }
            _ => return String::new(),
        }
        let c = match self.chart {
            Some(Panel::Structure { .. }) => Panel::Structure { of: self.mixture },
            _ => Panel::Heatmap { of: self.heat },
        };
        self.chart = Some(c);
        self.chart_note(c)
    }

    /// What the status says of chart `c`: its keys, or why it has no
    /// table to draw.
    fn chart_note(&self, c: Panel) -> String {
        let keys = "H next · esc the map";
        match self.data.chart_missing(c) {
            Some(why) => format!("no {}: {why} · {keys}", c.title()),
            None => {
                let what = match c {
                    Panel::Structure { .. } => "c topics or types",
                    _ => "c the next heatmap",
                };
                format!("{} · {what} · p save · {keys}", c.title())
            }
        }
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
        self.data.title(panel, &self.style)
    }

    /// Show the figures, or move to the next panel when they are shown.
    pub fn next_panel(&mut self) {
        self.chart = None;
        if self.shown {
            self.sel = (self.sel + 1) % self.panels.len();
            self.view = View::default();
        }
        self.shown = true;
    }

    /// Another pair of diffusion components, when that panel is shown.
    pub fn step_pair(&mut self, forward: bool, x_axis: bool) -> String {
        let Panel::Diffusion { x, y } = self.current() else {
            return if self.data.diffusion.is_some() {
                "the axes step on the diffusion map (m reaches it)".into()
            } else {
                "no diffusion map here: the axes step on a trajectory's".into()
            };
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
        let reach = self.reach();
        match self.with_view(|v| v.zoom(inward, reach)) {
            None => "only the scatter zooms (m and v choose it)".into(),
            Some((_, v)) if v.is_whole() => "the whole scatter".into(),
            Some((_, v)) => format!("zoom ×{:.1} · arrows pan · 0 shows it whole", v.factor()),
        }
    }

    /// Pan the scatter by a step (the arrows); `dx`, `dy` right and up.
    pub fn pan(&mut self, dx: f32, dy: f32) -> String {
        let reach = self.reach();
        match self.with_view(|v| v.pan(dx, dy, reach)) {
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
            moves: if scatter { self.moves } else { 0 },
        }
    }

    /// Open the thumbnail grid on the figure on screen (`w`).
    pub fn open_grid(&mut self) {
        self.chart = None;
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
        let (data, px) = (&self.data, self.cell_px());
        let moved: Vec<Offsets<String>> = missing
            .iter()
            .map(|((look, _), _)| self.moved(look.panel))
            .collect();
        let images: Vec<Result<RgbaImage, String>> = missing
            .par_iter()
            .zip(&moved)
            .map(|(((look, _), area), moved)| {
                let size = Size {
                    w: (u32::from(area.width) * px.0).max(32),
                    h: (u32::from(area.height) * px.1).max(32),
                    text_px: px.1 as f32,
                };
                let style = look.style.unwrap_or_default();
                data.figure(look.panel, size, &style, View::default(), moved)
                    .and_then(|f| figures::render(&f))
                    .map_err(|e| format!("{e:#}"))
            })
            .collect();
        let mut tiles = self.tiles.borrow_mut();
        if tiles.len() + missing.len() > TILES_KEPT {
            tiles.clear();
        }
        for ((key, area), img) in missing.into_iter().zip(images) {
            tiles.push((key, self.fit(img, area)));
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
        self.area.set(Some(area));
        if self
            .rendered
            .borrow()
            .as_ref()
            .is_none_or(|(k, _)| *k != key)
        {
            self.show(key, self.image(look, area), area);
        }
        Ref::map(self.rendered.borrow(), |s| {
            &s.as_ref().expect("filled above").1
        })
    }

    /// The figure as last drawn, without waiting: for a frame that goes
    /// out at once, the current figure for `area` once it is rasterised
    /// (which starts off the screen's thread), else whatever is drawn so
    /// far (none when nothing was), and whether the figure is still due.
    pub fn drawn_so_far(&self, area: Rect) -> (Option<Ref<'_, Drawn>>, bool) {
        let look = self.look(self.current(), self.view);
        let key = (look, (area.width, area.height));
        let current = |r: &Option<(Sized, Drawn)>| r.as_ref().is_some_and(|(k, _)| *k == key);
        self.area.set(Some(area));
        if !current(&self.rendered.borrow()) {
            if let Some(img) = self.image_ready(look, area) {
                self.show(key, img, area);
            }
        }
        let rendered = self.rendered.borrow();
        let behind = !current(&rendered);
        let drawn = Ref::filter_map(rendered, |s| s.as_ref().map(|(_, d)| d)).ok();
        (drawn, behind)
    }

    /// Put `img` on screen for `area` as `key`, keeping where its labels
    /// went.
    fn show(&self, key: Sized, img: Raster, area: Rect) {
        let on = img.as_ref().ok().and_then(|(_, on)| on.clone());
        *self.rendered.borrow_mut() = Some((key, self.fit(img.map(|(img, _)| img), area)));
        *self.on_screen.borrow_mut() = on;
    }

    /// `img` fitted to `area` as the terminal draws it, as senna fits its
    /// map.
    fn fit(&self, img: Result<RgbaImage, String>, area: Rect) -> Drawn {
        img.and_then(|img| fit_protocol(&self.picker, img, area).map_err(|e| format!("{e:#}")))
    }

    /// `look` rasterised from `data` for `area`, at its pixel size, its
    /// text a cell high, its labels moved by `moved`.
    fn rasterised(
        data: &CellFigures,
        look: Look,
        style: &Style,
        moved: &Offsets<String>,
        px: (u32, u32),
        area: Rect,
    ) -> Raster {
        let fig = data
            .figure(look.panel, frame_size(px, area), style, look.view, moved)
            .map_err(|e| format!("{e:#}"))?;
        let on = match &fig {
            Picture::Map(m) => Some(Arc::new(OnMap::of(look.panel, m))),
            Picture::Svg(_) | Picture::Chart(_) => None,
        };
        let img = figures::render(&fig).map_err(|e| format!("{e:#}"))?;
        Ok((img, on))
    }

    /// The scatter's frame in pixels, as last drawn on screen; before
    /// that, an export's.
    pub(crate) fn frame(&self) -> (u32, u32) {
        match self.rendered.borrow().as_ref() {
            Some(((_, (w, h)), _)) => {
                let s = frame_size(self.cell_px(), Rect::new(0, 0, *w, *h));
                (s.w, s.h)
            }
            None => {
                let s = figures::export_size(None);
                (s.w, s.h)
            }
        }
    }

    /// How much of the scatter on screen its frame shows when whole
    /// ([`View::reach`]); the whole extent when it has no cells to fit.
    fn reach(&self) -> [f32; 2] {
        let (w, h) = self.frame();
        match self.data.whole(self.current()) {
            Ok(whole) => View::reach(whole, w as usize, h as usize),
            Err(_) => [1.0, 1.0],
        }
    }

    /// The terminal cell's size in frame pixels.
    fn cell_px(&self) -> (u32, u32) {
        let (w, h) = legume_tui::map::display::cell(&self.picker);
        (w as u32, h as u32)
    }

    /// `look` rasterised, when it is the one last drawn.
    fn cached(&self, look: Look) -> Option<Raster> {
        let cache = self.image.borrow();
        cache
            .as_ref()
            .filter(|(l, _)| *l == look)
            .map(|(_, img)| img.clone())
    }

    /// `look` rasterised, from the cache when it is the one last drawn;
    /// else waiting for the one being rasterised, or drawing it here.
    fn image(&self, look: Look, area: Rect) -> Raster {
        if let Some(img) = self.cached(look) {
            return img;
        }
        let waiting = self.pending.borrow_mut().take();
        let img = match waiting {
            Some((l, rx)) if l == look => rx.recv().unwrap_or_else(|_| Err(NOT_DRAWN.into())),
            _ => {
                let moved = self.moved(look.panel);
                Self::rasterised(&self.data, look, &self.style, &moved, self.cell_px(), area)
            }
        };
        *self.image.borrow_mut() = Some((look, img.clone()));
        img
    }

    /// `look` rasterised when it is, without waiting: from the cache, or
    /// from the one being rasterised; else that is started (off the
    /// screen's thread) and `None` said.
    fn image_ready(&self, look: Look, area: Rect) -> Option<Raster> {
        if let Some(img) = self.cached(look) {
            return Some(img);
        }
        // One figure is rasterised at a time: a look moved past is let
        // finish before the latest starts, so keys held down start few.
        let mut pending = self.pending.borrow_mut();
        if let Some((l, rx)) = pending.as_ref() {
            let img = match rx.try_recv() {
                Ok(img) => img,
                Err(TryRecvError::Empty) => return None,
                Err(TryRecvError::Disconnected) => Err(NOT_DRAWN.into()),
            };
            let l = *l;
            *pending = None;
            if l == look {
                *self.image.borrow_mut() = Some((look, img.clone()));
                return Some(img);
            }
        }
        let (tx, rx) = channel();
        let (data, style, px) = (Arc::clone(&self.data), self.style, self.cell_px());
        let moved = self.moved(look.panel);
        std::thread::spawn(move || {
            let _ = tx.send(Self::rasterised(&data, look, &style, &moved, px, area));
        });
        *pending = Some((look, rx));
        None
    }

    /// Export the current panel as it is on screen and log it; with the
    /// grid open, its selected tile as the tile shows it (whole).
    pub fn export(&mut self) -> Result<String> {
        let (panel, view, frame) = match self.selected_tile() {
            Some(tile) => (tile, View::default(), None),
            None => (self.current(), self.view, Some(self.frame())),
        };
        let moved = self.moved(panel);
        let e = self.data.export(panel, &self.style, view, frame, &moved)?;
        let logged = self.log(&e.files, &e.what, e.picture.as_ref(), panel);
        let exts: Vec<String> = e
            .files
            .iter()
            .filter_map(|f| Some(format!(".{}", f.extension()?.to_string_lossy())))
            .collect();
        Ok(format!(
            "wrote {} {}{logged}",
            e.base.display(),
            exts.join(" and ")
        ))
    }
}

impl FigurePane {
    /// Log `files`, an export of `panel` showing `picture`, in the gallery
    /// and check the strip again; why it is not listed, if it is not.
    fn log(
        &mut self,
        files: &[PathBuf],
        what: &str,
        picture: Option<&RgbaImage>,
        panel: Panel,
    ) -> String {
        let trace = Trace::of(&self.data.slug(panel), &self.data.manifest);
        let logged = self
            .exports
            .gallery
            .add(files, what, picture, trace)
            .map(|()| String::new())
            .unwrap_or_else(|e| format!(" (not listed: {e:#})"));
        self.exports.refresh();
        logged
    }

    /// Whether the save dialog is open.
    pub fn saving(&self) -> bool {
        self.save.is_some()
    }

    /// `p` on a scatter: the save dialog, for the map on screen or, from
    /// the grid, its tile's map.
    fn ask_save(&mut self, panel: Panel, tile: bool) -> FigureKey {
        let target = MapTarget {
            panel,
            view: if tile { View::default() } else { self.view },
            frame: self.frame(),
            tile,
            chart: panel.is_chart(),
            stem: format!(
                "{}.{}.{}",
                self.data.prefix(),
                self.data.kind(),
                self.data.slug(panel)
            ),
            title: self.title(panel),
        };
        let first = target.scopes()[0];
        let dialog = SaveDialog::new(&target, first, Numbering::Suffix);
        self.save = Some(Saving { target, dialog });
        FigureKey::Took
    }

    /// A key for the save dialog: it takes every key while it is open.
    pub fn save_key(&mut self, k: &KeyEvent) -> FigureKey {
        let Some(mut s) = self.save.take() else {
            return FigureKey::Ignored;
        };
        match s.dialog.key(*k, &s.target) {
            SaveKey::Open => {
                self.save = Some(s);
                FigureKey::Took
            }
            SaveKey::Cancelled => FigureKey::Took,
            SaveKey::Save(req) => match self.save_map(&s.target, &req) {
                Ok(said) => FigureKey::Say(said),
                // The dialog stays, saying why, for another name or size.
                Err(e) => {
                    s.dialog.say(format!("not saved: {e:#}"));
                    self.save = Some(s);
                    FigureKey::Took
                }
            },
        }
    }

    /// Write the map `req` asks for and log it with its thumbnail.
    fn save_map(&mut self, t: &MapTarget, req: &SaveRequest<MapScope>) -> Result<String> {
        let page = PageSize {
            width_in: req.width_in,
            dpi: req.dpi,
        };
        let size = Size {
            w: req.px.0 as u32,
            h: req.px.1 as u32,
            text_px: page.cell_px(),
        };
        if t.chart {
            let img = (self.data).write_chart(t.panel, size, req.dpi as f32, &req.path)?;
            let what = format!("{} · {} in · {} dpi", t.title, req.width_in, req.dpi);
            let logged = self.log(&req.files, &what, Some(&img), t.panel);
            return Ok(format!("wrote {}{logged}", req.path.display()));
        }
        let view = match req.scope {
            MapScope::Shown => t.view,
            MapScope::Whole => View::default(),
        };
        let moved = self.moved(t.panel);
        let img = self.data.write_map(
            t.panel,
            size,
            req.dpi as f32,
            &self.style,
            view,
            &moved,
            &req.path,
        )?;
        let zoomed = if view.is_whole() {
            String::new()
        } else {
            format!(" · zoomed ×{:.1}", view.factor())
        };
        let what = format!(
            "{}{zoomed} · {} in · {} dpi",
            self.title(t.panel),
            req.width_in,
            req.dpi
        );
        let logged = self.log(&req.files, &what, Some(&img), t.panel);
        Ok(format!("wrote {}{logged}", req.path.display()))
    }

    /// Draw the save dialog, when it is open, in the middle of `area`.
    pub fn draw_save(&self, f: &mut ratatui::Frame, area: Rect) {
        if let Some(s) = &self.save {
            s.dialog.render(f, area, &s.target, &Chrome::terminal());
        }
    }

    /// The labels moved on `panel`, by name: none when the book cannot be
    /// read.
    fn moved(&self, panel: Panel) -> Offsets<String> {
        match &self.book {
            Ok(b) => b.layout(&self.layout_id(panel)),
            Err(_) => Offsets::default(),
        }
    }

    fn layout_id(&self, panel: Panel) -> String {
        layout_id(&self.data, panel)
    }

    /// The scatter on screen as drawn, when it is the panel shown.
    fn shown_map(&self) -> Option<Arc<OnMap>> {
        let on = self.on_screen.borrow().clone()?;
        (self.shown && self.grid.is_none() && on.panel == self.current()).then_some(on)
    }

    /// What mouse event `m` does on the scatter on screen: a drag pans it,
    /// the wheel zooms it at the pointer, a press on a label drags the
    /// label, and a click or a hover finds the cluster under the pointer.
    pub fn mouse(&mut self, m: &MouseEvent) -> MapMouse {
        let (Some(rect), Some(on)) = (self.area.get(), self.shown_map()) else {
            return self.let_go().map_or(MapMouse::None, MapMouse::Say);
        };
        let cell = legume_tui::map::display::cell(&self.picker);
        let act = self.pointer.handle(&MapArea { rect, cell }, m, |px, py| {
            on.hits.at(&on.vp, px, py)
        });
        let (w, h) = (on.vp.w, on.vp.h);
        let reply = match act {
            MapAction::None => MapMouse::None,
            MapAction::Hover { px, py } => MapMouse::Hover(self.cluster_at(&on, px, py, cell)),
            MapAction::Away => MapMouse::Hover(None),
            MapAction::Click { px, py, .. } => MapMouse::Click(self.cluster_at(&on, px, py, cell)),
            MapAction::Pan { dx, dy } => self.move_view(|v, whole| v.pan_px(dx, dy, whole, w, h)),
            MapAction::Zoom { factor, at } => {
                self.move_view(|v, whole| v.zoom_at(factor > 1.0, at, whole, w, h))
            }
            MapAction::MoveLabel { key, dx, dy } if on.line > 0.0 => {
                self.move_label(&on, key, (dx / on.line, dy / on.line))
            }
            MapAction::MoveLabel { .. } => MapMouse::None,
        };
        match self.keep_labels() {
            Some(e) => MapMouse::Say(e),
            None => reply,
        }
    }

    /// Let go of a press and stop moving labels, as when the map goes
    /// away. A label dragged so far is kept; what went wrong keeping it.
    pub fn let_go(&mut self) -> Option<String> {
        self.pointer.cancel();
        self.edit = None;
        self.keep_labels()
    }

    /// Write the labels moved by a drag once no press is held, however the
    /// press ended; what went wrong.
    fn keep_labels(&mut self) -> Option<String> {
        (self.unsaved && !self.pointer.dragging())
            .then(|| self.save_labels().err())
            .flatten()
    }

    /// Move the view by `f` on the whole extent of the scatter.
    fn move_view(
        &mut self,
        f: impl FnOnce(&mut View, legume_tui::map::viewport::Bounds),
    ) -> MapMouse {
        let Ok(whole) = self.data.whole(self.current()) else {
            return MapMouse::None;
        };
        let before = self.view;
        f(&mut self.view, whole);
        if self.view == before {
            MapMouse::None
        } else {
            MapMouse::Redraw
        }
    }

    /// Move label `key` of `on` by `d_em` lines.
    fn move_label(&mut self, on: &OnMap, key: u32, d_em: (f32, f32)) -> MapMouse {
        let Some(name) = on.names.get(key as usize) else {
            return MapMouse::None;
        };
        let id = self.layout_id(on.panel);
        match &mut self.book {
            Err(why) => MapMouse::Say(why.clone()),
            Ok(book) => {
                let mut o = book.layout(&id);
                o.nudge(name.to_string(), d_em);
                book.set_layout(&id, &o);
                self.moves += 1;
                self.unsaved = true;
                MapMouse::Redraw
            }
        }
    }

    /// Write the moved labels beside the run.
    fn save_labels(&mut self) -> Result<(), String> {
        self.unsaved = false;
        match &self.book {
            Ok(b) => b.save().map_err(|e| {
                format!(
                    "the moved labels were not kept in {}: {e}",
                    b.path().display()
                )
            }),
            Err(why) => Err(why.clone()),
        }
    }

    /// The cluster under pixel `(px, py)` of `on`: a cluster's label names
    /// it, else the cell drawn nearest within a terminal cell.
    fn cluster_at(&self, on: &OnMap, px: f32, py: f32, cell: (f32, f32)) -> Option<i32> {
        if self.data.cluster.is_empty() {
            return None;
        }
        if self.style.colouring == Colouring::Cluster {
            let named = on
                .hits
                .at(&on.vp, px, py)
                .and_then(|k| on.names.get(k as usize)?.strip_prefix('K')?.parse().ok());
            if named.is_some() {
                return named;
            }
        }
        let i = self
            .data
            .cell_at(on.panel, &on.vp, px, py, cell.0.max(cell.1))?;
        let k = *self.data.cluster.get(i)?;
        (k >= 0).then_some(k)
    }

    /// Moving labels from the keys is on.
    pub fn editing(&self) -> bool {
        self.edit.is_some()
    }

    /// Key `k` for moving labels: `E` starts it; then Tab and BackTab
    /// choose a label, the arrows move it (Shift: further), Backspace puts
    /// it back on its cells and Esc ends it. Other keys are ignored.
    pub fn edit_key(&mut self, k: &KeyEvent) -> FigureKey {
        use FigureKey::{Ignored, Say, Took};
        if self.grid.is_some() {
            return Ignored;
        }
        let Some(mut edit) = self.edit else {
            if !LabelEdit::<u32>::starts(k) {
                return Ignored;
            }
            return self.start_editing();
        };
        let Some(on) = self.shown_map() else {
            self.edit = None;
            return Say("the labels move on the scatter on screen".into());
        };
        let id = self.layout_id(on.panel);
        let Ok(book) = &mut self.book else {
            self.edit = None;
            return Ignored;
        };
        let mut by_name = book.layout(&id);
        let mut o = on.by_place(&by_name);
        let did = edit.key(k, &on.placed, &mut o);
        self.edit = Some(edit);
        match did {
            EditKey::Pass => Ignored,
            EditKey::Leave => {
                self.edit = None;
                Took
            }
            EditKey::Taken => Say(editing_note(&on, edit)),
            EditKey::Changed => {
                for (k, n) in on.names.iter().enumerate() {
                    let n = n.to_string();
                    by_name.reset(&n);
                    if let Some(d) = o.get(&(k as u32)) {
                        by_name.nudge(n, d);
                    }
                }
                book.set_layout(&id, &by_name);
                self.moves += 1;
                Say(match self.save_labels() {
                    Ok(()) => editing_note(&on, edit),
                    Err(e) => e,
                })
            }
        }
    }

    /// `E`: moving labels from the keys, on the first label as read.
    fn start_editing(&mut self) -> FigureKey {
        if let Err(why) = &self.book {
            return FigureKey::Say(why.clone());
        }
        if !self.current().is_scatter() {
            return FigureKey::Say("labels move on the scatter (m and v choose it)".into());
        }
        let Some(on) = self.shown_map().filter(|on| !on.placed.is_empty()) else {
            return FigureKey::Say("no labels on screen to move: t shows them".into());
        };
        let edit = LabelEdit::enter(&on.placed);
        self.edit = Some(edit);
        FigureKey::Say(editing_note(&on, edit))
    }
}

/// `panel`'s key in the labels' book: the round and the layout, so a label
/// moved on one round's UMAP stays off the next round's and the PHATE.
fn layout_id(data: &CellFigures, panel: Panel) -> String {
    let prefix = data.prefix();
    let round = Path::new(&prefix)
        .file_name()
        .map_or(prefix.clone(), |n| n.to_string_lossy().into_owned());
    format!("{round}:{}", data.slug(panel))
}

/// The file beside the run of `data` that keeps its moved labels.
fn labels_file(data: &CellFigures) -> PathBuf {
    PathBuf::from(format!("{}.{LABELS_SUFFIX}", data.prefix()))
}

/// What a save of a scatter takes in: the map as it is on screen (zoomed
/// and panned, in the screen's shape), or the whole map.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MapScope {
    Shown,
    Whole,
}

/// The scatter a save is of, as the save dialog asks about it: its view
/// and the shape on screen when `p` was pressed.
struct MapTarget {
    panel: Panel,
    view: View,
    frame: (u32, u32),
    /// From the grid, whose tile shows the scatter whole.
    tile: bool,
    /// A round's chart, saved as it is drawn.
    chart: bool,
    /// `{prefix}.{kind}.{panel}`, numbered past saved files.
    stem: String,
    title: String,
}

impl SaveTarget for MapTarget {
    type Scope = MapScope;

    fn scopes(&self) -> Vec<MapScope> {
        if self.chart {
            vec![MapScope::Shown]
        } else if self.tile {
            vec![MapScope::Whole]
        } else {
            vec![MapScope::Shown, MapScope::Whole]
        }
    }

    fn scope_name(&self, s: MapScope) -> String {
        match s {
            MapScope::Shown if self.chart => self.title.clone(),
            MapScope::Shown if self.view.is_whole() => format!("{} as on screen", self.title),
            MapScope::Shown => format!(
                "{} as on screen, zoomed ×{:.1}",
                self.title,
                self.view.factor()
            ),
            MapScope::Whole => format!("{}, the whole map", self.title),
        }
    }

    fn default_stem(&self, _: MapScope) -> String {
        self.stem.clone()
    }

    /// A PDF only, as senna saves its map.
    fn exts(&self, _: MapScope) -> Vec<String> {
        vec!["pdf".into()]
    }

    /// The screen's shape, so the page shows the data the screen shows;
    /// the whole map at 4:3.
    fn aspect(&self, s: MapScope) -> f32 {
        match (s, self.frame) {
            (MapScope::Shown, (w, h)) if w > 0 && h > 0 => w as f32 / h as f32,
            _ => 4.0 / 3.0,
        }
    }
}

/// The save dialog open on a scatter.
struct Saving {
    target: MapTarget,
    dialog: SaveDialog<MapScope>,
}

/// `lupin annotate -f ROUND --pdf FILE`: the round's cell map as the TUI
/// opens on it (its first layout, whole, in its first colouring, the
/// labels where they were moved), saved as a PDF at `pdf`'s size without
/// the TUI. Where it went. Never over the round's manifest, only to a
/// `.pdf`, and over a file already there only as `--overwrite` allows.
pub fn round_pdf(from: &str, pdf: &PdfArgs) -> Result<PathBuf> {
    let path = pdf.pdf.clone().context("no --pdf file")?;
    let loaded = load(from)?;
    anyhow::ensure!(
        !crate::manifest::run::same_file(&path, &loaded.file),
        "--pdf {} is the round's manifest: choose another file",
        path.display()
    );
    anyhow::ensure!(
        path.extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("pdf")),
        "--pdf {} is not a .pdf",
        path.display()
    );
    crate::manifest::run::may_replace(&path)?;
    let cells = crate::manifest::rounds::read_cells(&loaded)?;
    let score = crate::manifest::rounds::read_score_cells(&loaded, &cells)?;
    let data = CellFigures::of_round(&loaded, &cells, score);
    if data.layouts.is_empty() {
        anyhow::bail!(
            "{} has no layout to draw its cells on{}",
            loaded.file.display(),
            data.layout_hint
                .as_deref()
                .map_or(String::new(), |h| format!(": {h}"))
        );
    }
    let page = pdf.size();
    // Height over width, as the page sizes it.
    let (w, h) = page
        .check(3.0 / 4.0)
        .map_err(|e| anyhow::anyhow!("--width {} --dpi {}: {e}", pdf.width, pdf.dpi))?;
    let panel = data.panels()[0];
    let moved = match OffsetBook::load(&labels_file(&data)) {
        Ok(b) => b.layout(&layout_id(&data, panel)),
        Err(e) => {
            log::warn!("the moved labels cannot be read ({e}): drawn where they fall");
            Offsets::default()
        }
    };
    let size = Size {
        w: w as u32,
        h: h as u32,
        text_px: page.cell_px(),
    };
    data.write_map(
        panel,
        size,
        page.dpi as f32,
        &Style::opening(&data),
        View::default(),
        &moved,
        &path,
    )?;
    Ok(path)
}

/// The status while labels are moved from the keys: the one chosen, and
/// the keys.
fn editing_note(on: &OnMap, edit: LabelEdit<u32>) -> String {
    let name = edit
        .at
        .and_then(|k| on.names.get(k as usize))
        .map_or("no label", AsRef::as_ref);
    format!("moving {name} · arrows move it (shift: further) · tab next · ⌫ back on its cells · esc done")
}

/// The pixels a figure is rasterised at for `area`, `px` per cell, its
/// text a cell high.
fn frame_size(px: (u32, u32), area: Rect) -> Size {
    Size {
        w: (u32::from(area.width) * px.0).max(64),
        h: (u32::from(area.height) * px.1).max(64),
        text_px: px.1 as f32,
    }
}

#[cfg(test)]
#[path = "tests/figure_pane.rs"]
pub(crate) mod tests;
