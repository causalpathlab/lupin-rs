//! Choosing the files the TUI was not given on the command line: a run
//! manifest, then a marker panel, from a browser over the file system; the
//! same browser opens in a popup for the files the TUI's forms ask for.
//!
//! A directory lists at once, from its entries' names; marker panels are read
//! and scored against the run's genes on a thread behind it, likeliest first,
//! and remembered for as long as the browser is open.

use crate::annotate::gene_rows::GeneRows;
use crate::annotate::markers::read_marker_pairs;
use crate::manifest::{pinto, run};
use anyhow::Result;
use ratatui::crossterm::event::{self, Event, KeyCode, KeyEventKind};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::style::{Color, Modifier, Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, List, ListItem, ListState, Paragraph};
use ratatui::{DefaultTerminal, Frame};
use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::sync::{mpsc, Arc};
use std::time::Duration;

/// Files larger than this are not read as marker panels.
const MAX_PANEL_BYTES: u64 = 8 << 20;

/// What the browser is choosing.
pub enum Want {
    /// A run manifest (`.senna.json`, `.lupin.json`, `.pinto.json`).
    Manifest,
    /// A marker panel, scored against the run's genes when they are known.
    Markers(Option<Arc<GeneRows>>, usize),
    /// A `cell<TAB>type` labels table (`.tsv`, `.txt`, `.csv`, or gzipped).
    Labels,
}

/// A marker panel's size and how much of it the run has.
#[derive(Clone, Copy)]
struct Panel {
    types: usize,
    genes: usize,
    matched: Option<usize>,
}

struct Entry {
    name: String,
    path: PathBuf,
    dir: bool,
    panel: Option<Panel>,
    /// A file still to be read as a marker panel.
    pending: bool,
}

impl Entry {
    fn fits(&self, want: &Want) -> bool {
        match want {
            Want::Manifest => is_manifest(&self.path),
            Want::Markers(..) => self.panel.is_some() || self.pending,
            Want::Labels => is_table(&self.path),
        }
    }
}

fn is_table(p: &Path) -> bool {
    let n = p.to_string_lossy();
    let n = n.strip_suffix(".gz").unwrap_or(&n);
    [".tsv", ".txt", ".csv"].iter().any(|ext| n.ends_with(ext))
}

pub struct Picker {
    title: &'static str,
    want: Want,
    dir: PathBuf,
    entries: Vec<Entry>,
    state: ListState,
    /// List every file, not only those that fit.
    all: bool,
    /// The panel covering the most of the run's genes.
    best: Option<PathBuf>,
    /// Why the last choice was refused.
    note: Option<String>,
    /// Files already read as marker panels (or found not to be one).
    read: HashMap<PathBuf, Option<Panel>>,
    /// Panels read behind the listing, as they come.
    scan: Option<mpsc::Receiver<(PathBuf, Option<Panel>)>>,
    /// The selection was moved by hand: arriving panels leave it be.
    moved: bool,
}

fn is_manifest(p: &Path) -> bool {
    let n = p.to_string_lossy();
    n.ends_with(".senna.json") || n.ends_with(run::LUPIN_SUFFIX) || pinto::is_pinto(p)
}

/// Whether a file by this name may be a marker panel, worth reading.
fn may_be_panel(path: &Path) -> bool {
    let skip = [
        "parquet", "zarr", "h5", "h5ad", "json", "bam", "bai", "png", "pdf", "log", "zip",
    ];
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    !skip.contains(&ext)
}

/// `path` read as a marker panel, if it is one: at least one pair, not every
/// label a number (a count table), and, when the run's genes are known, at
/// least one of them (not a log or some other table).
fn read_panel(path: &Path, index: Option<&GeneRows>) -> Option<Panel> {
    if path.metadata().ok()?.len() > MAX_PANEL_BYTES {
        return None;
    }
    let pairs = read_marker_pairs(path.to_str()?).ok()?;
    if pairs.iter().all(|(_, t)| t.parse::<f64>().is_ok()) {
        return None;
    }
    let types: BTreeSet<&str> = pairs.iter().map(|(_, t)| &**t).collect();
    let genes: BTreeSet<&str> = pairs.iter().map(|(g, _)| &**g).collect();
    let matched = index.map(|ix| genes.iter().filter(|g| ix.match_rows(g).is_some()).count());
    if matched == Some(0) {
        return None;
    }
    Some(Panel {
        types: types.len(),
        genes: genes.len(),
        matched,
    })
}

impl Picker {
    fn open(&mut self, dir: PathBuf) {
        let dir = dir.canonicalize().unwrap_or(dir);
        let markers = matches!(self.want, Want::Markers(..));
        let mut entries: Vec<Entry> = std::fs::read_dir(&dir)
            .into_iter()
            .flatten()
            .flatten()
            .filter(|e| !e.file_name().to_string_lossy().starts_with('.'))
            .map(|e| {
                let path = e.path();
                // The entry's own type, without a stat; a link is followed.
                let is_dir = match e.file_type() {
                    Ok(t) if !t.is_symlink() => t.is_dir(),
                    _ => path.is_dir(),
                };
                // A `.zarr` store is a directory, but never one to browse into here.
                let dir = is_dir && path.extension().is_none_or(|x| x != "zarr");
                let (panel, pending) = match self.read.get(&path) {
                    _ if dir || !markers => (None, false),
                    Some(known) => (*known, false),
                    None => (None, !is_dir && may_be_panel(&path)),
                };
                Entry {
                    name: e.file_name().to_string_lossy().into_owned(),
                    path,
                    dir,
                    panel,
                    pending,
                }
            })
            .collect();
        entries.sort_by(|a, b| b.dir.cmp(&a.dir).then_with(|| a.name.cmp(&b.name)));
        if let Some(up) = dir.parent() {
            entries.insert(
                0,
                Entry {
                    name: "..".into(),
                    path: up.to_path_buf(),
                    dir: true,
                    panel: None,
                    pending: false,
                },
            );
        }
        self.dir = dir;
        self.entries = entries;
        self.moved = false;
        self.start_scan();
        self.settle(None);
    }

    /// Read the listing's pending files as panels on a thread, the names that
    /// mention markers first. A scan of the directory left behind stops at
    /// its next file.
    fn start_scan(&mut self) {
        self.scan = None;
        let Want::Markers(index, _) = &self.want else {
            return;
        };
        let mut todo: Vec<PathBuf> = self
            .entries
            .iter()
            .filter(|e| e.pending)
            .map(|e| e.path.clone())
            .collect();
        if todo.is_empty() {
            return;
        }
        todo.sort_by_key(|p| {
            let name = p.file_name().map(|n| n.to_string_lossy().to_lowercase());
            let name = name.unwrap_or_default();
            (!name.contains("marker"), name.len())
        });
        let index = index.clone();
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            for p in todo {
                let panel = read_panel(&p, index.as_deref());
                if tx.send((p, panel)).is_err() {
                    break;
                }
            }
        });
        self.scan = Some(rx);
    }

    /// Take the panels read since the last call; whether any came.
    pub fn poll(&mut self) -> bool {
        let Some(rx) = &self.scan else {
            return false;
        };
        let mut came = Vec::new();
        loop {
            match rx.try_recv() {
                Ok(r) => came.push(r),
                Err(mpsc::TryRecvError::Empty) => break,
                Err(mpsc::TryRecvError::Disconnected) => {
                    // The reader is done, or gone early: what it left is
                    // unread, no longer reading.
                    self.scan = None;
                    let left: Vec<PathBuf> = self
                        .entries
                        .iter()
                        .filter(|e| e.pending && !came.iter().any(|(p, _)| *p == e.path))
                        .map(|e| e.path.clone())
                        .collect();
                    came.extend(left.into_iter().map(|p| (p, None)));
                    break;
                }
            }
        }
        if came.is_empty() {
            return false;
        }
        let at = self.selected_path();
        for (path, panel) in came {
            self.learn(path, panel);
        }
        self.settle(at);
        true
    }

    /// Record what `path` turned out to be.
    fn learn(&mut self, path: PathBuf, panel: Option<Panel>) {
        if let Some(e) = self.entries.iter_mut().find(|e| e.path == path) {
            (e.panel, e.pending) = (panel, false);
        }
        self.read.insert(path, panel);
    }

    fn selected_path(&self) -> Option<PathBuf> {
        let at = self.state.selected()?;
        self.shown().get(at).map(|e| e.path.clone())
    }

    /// The best panel so far, and the selection: on it, else the first file
    /// that fits, until moved by hand; then on `keep` where it is still shown.
    fn settle(&mut self, keep: Option<PathBuf>) {
        self.best = self
            .entries
            .iter()
            .filter_map(|e| {
                e.panel
                    .as_ref()
                    .filter(|p| p.types > 1)
                    .map(|p| (p, &e.path))
            })
            .filter_map(|(p, path)| Some((p.matched?, path)))
            .max_by_key(|(m, _)| *m)
            .filter(|(m, _)| *m > 0)
            .map(|(_, p)| p.clone());
        let shown = self.shown();
        let at = if self.moved {
            keep.and_then(|k| shown.iter().position(|e| e.path == k))
                .unwrap_or_else(|| {
                    let at = self.state.selected().unwrap_or(0);
                    at.min(shown.len().saturating_sub(1))
                })
        } else {
            shown
                .iter()
                .position(|e| Some(&e.path) == self.best.as_ref())
                .or_else(|| shown.iter().position(|e| e.fits(&self.want)))
                .unwrap_or(0)
        };
        self.state.select(Some(at));
    }

    fn shown(&self) -> Vec<&Entry> {
        self.entries
            .iter()
            .filter(|e| self.all || e.dir || e.fits(&self.want))
            .collect()
    }

    fn row(&self, e: &Entry) -> ListItem<'static> {
        if e.dir {
            return ListItem::new(Line::from(format!("  {}/", e.name)).fg(Color::Cyan));
        }
        let star = if Some(&e.path) == self.best.as_ref() {
            "★ "
        } else {
            "  "
        };
        let mut spans = vec![Span::raw(format!("{star}{}", e.name))];
        if e.pending {
            spans.push(Span::raw("   …").fg(Color::Indexed(244)));
        }
        if let Some(p) = &e.panel {
            let cover = match (p.matched, &self.want) {
                (Some(m), Want::Markers(_, n)) => format!(" · {m}/{} in the run's {n}", p.genes),
                _ => String::new(),
            };
            spans.push(
                Span::raw(format!("   {} types · {} genes{cover}", p.types, p.genes))
                    .fg(Color::Indexed(250)),
            );
        }
        let line = Line::from(spans);
        ListItem::new(if e.fits(&self.want) { line } else { line.dim() })
    }

    /// The directory shown.
    pub fn dir(&self) -> &Path {
        &self.dir
    }

    pub fn title(&self) -> &'static str {
        self.title
    }

    /// The key line: why the last choice was refused, else the keys.
    pub fn hint(&self, cancel: &str) -> (String, bool) {
        if let Some(note) = &self.note {
            return (format!(" {note} "), true);
        }
        let left = self.entries.iter().filter(|e| e.pending).count();
        let reading = if left > 0 {
            format!("reading {left} file{} · ", if left == 1 { "" } else { "s" })
        } else {
            String::new()
        };
        let all = if self.all {
            "a fitting only"
        } else {
            "a all files"
        };
        (
            format!(
                " {reading}↑↓ pgup/dn move · enter open/choose · ← up · {all} · ~ home · {cancel} "
            ),
            false,
        )
    }

    /// The entries, unframed, in `area`.
    pub fn draw_list(&self, f: &mut Frame, area: Rect) {
        let items: Vec<ListItem> = self.shown().into_iter().map(|e| self.row(e)).collect();
        let list =
            List::new(items).highlight_style(Style::default().add_modifier(Modifier::REVERSED));
        f.render_stateful_widget(list, area, &mut self.state.clone());
    }

    /// The browser on the whole screen, before the TUI starts.
    fn draw(&self, f: &mut Frame) {
        let [head, body, keys] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Min(3),
            Constraint::Length(1),
        ])
        .areas(f.area());
        f.render_widget(
            Paragraph::new(format!(" {} — {}", self.title, self.dir.display())).bold(),
            head,
        );
        let block = Block::default().borders(Borders::ALL);
        self.draw_list(f, block.inner(body));
        f.render_widget(block, body);
        let (hint, note) = self.hint("q cancel");
        let hint = Paragraph::new(hint);
        f.render_widget(
            if note {
                hint.fg(Color::LightYellow)
            } else {
                hint.dim()
            },
            keys,
        );
    }

    /// Browse from `start` for a file of the kind `want`.
    pub fn new(title: &'static str, start: &Path, want: Want) -> Self {
        let mut p = Picker {
            title,
            want,
            dir: PathBuf::new(),
            entries: Vec::new(),
            state: ListState::default(),
            all: false,
            best: None,
            note: None,
            read: HashMap::new(),
            scan: None,
            moved: false,
        };
        p.open(start.to_path_buf());
        p
    }

    /// Take a key: move, open a directory, choose a file or cancel.
    pub fn key(&mut self, code: KeyCode) -> Step {
        self.note = None;
        let at = self.state.selected().unwrap_or(0);
        let mut sel = at;
        if super::app::step(&mut sel, self.shown().len(), code) {
            self.state.select(Some(sel));
            self.moved = true;
            return Step::Stay;
        }
        match code {
            KeyCode::Char('q') | KeyCode::Esc => return Step::Cancelled,
            KeyCode::Left | KeyCode::Backspace => {
                if let Some(up) = self.dir.parent().map(Path::to_path_buf) {
                    self.open(up);
                }
            }
            KeyCode::Char('~') => {
                if let Some(home) = std::env::var_os("HOME") {
                    self.open(PathBuf::from(home));
                }
            }
            KeyCode::Char('a') => {
                self.all = !self.all;
                self.state.select(Some(0));
                self.moved = true;
            }
            KeyCode::Enter | KeyCode::Right => {
                let Some((dir, path, pending)) = self
                    .shown()
                    .get(at)
                    .map(|e| (e.dir, e.path.clone(), e.pending))
                else {
                    return Step::Stay;
                };
                if dir {
                    self.open(path);
                    return Step::Stay;
                }
                // A file not read yet is read now, rather than waited for;
                // the selection stays on it, which any note is about.
                if pending {
                    if let Want::Markers(index, _) = &self.want {
                        let panel = read_panel(&path, index.as_deref());
                        self.learn(path.clone(), panel);
                        self.moved = true;
                        self.settle(Some(path.clone()));
                    }
                }
                let panel = self.read.get(&path).copied().flatten();
                match panel {
                    Some(p) if p.types < 2 => {
                        self.note =
                            Some("one cell type: annotation needs types to choose between".into());
                    }
                    None if pending && !self.all => {
                        self.note = Some("not a marker panel".into());
                    }
                    _ => return Step::Chosen(path),
                }
            }
            _ => {}
        }
        Step::Stay
    }
}

/// What a key did in the browser.
pub enum Step {
    Stay,
    Chosen(PathBuf),
    Cancelled,
}

/// Browse from `start` for a file of the kind `want` on the whole screen;
/// `None` when cancelled.
pub fn pick(title: &'static str, start: &Path, want: Want) -> Result<Option<PathBuf>> {
    let mut p = Picker::new(title, start, want);
    let mut terminal: DefaultTerminal = ratatui::init();
    let chosen = (|| -> Result<Option<PathBuf>> {
        let mut dirty = true;
        loop {
            dirty |= p.poll();
            if dirty {
                terminal.draw(|f| p.draw(f))?;
                dirty = false;
            }
            // Panels still being read redraw as they come; else wait for a key.
            if p.scan.is_some() && !event::poll(Duration::from_millis(100))? {
                continue;
            }
            dirty = true;
            let Event::Key(k) = event::read()? else {
                continue;
            };
            if k.kind != KeyEventKind::Press {
                continue;
            }
            match p.key(k.code) {
                Step::Stay => {}
                Step::Chosen(file) => return Ok(Some(file)),
                Step::Cancelled => return Ok(None),
            }
        }
    })();
    ratatui::restore();
    chosen
}

#[cfg(test)]
#[path = "tests/picker.rs"]
mod tests;
