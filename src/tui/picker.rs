//! Choosing the files the TUI was not given on the command line: a run
//! manifest, then a marker panel, from a browser over the file system; the
//! same browser opens in a popup for the files the TUI's forms ask for.

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
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

/// Files larger than this are not read as marker panels.
const MAX_PANEL_BYTES: u64 = 8 << 20;

/// What the browser is choosing.
pub enum Want {
    /// A run manifest (`.senna.json`, `.lupin.json`, `.pinto.json`).
    Manifest,
    /// A marker panel, scored against the run's genes when they are known.
    Markers(Option<Box<GeneRows>>, usize),
    /// A `cell<TAB>type` labels table (`.tsv`, `.txt`, `.csv`, or gzipped).
    Labels,
}

/// A marker panel's size and how much of it the run has.
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
}

impl Entry {
    fn fits(&self, want: &Want) -> bool {
        match want {
            Want::Manifest => is_manifest(&self.path),
            Want::Markers(..) => self.panel.is_some(),
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
}

fn is_manifest(p: &Path) -> bool {
    let n = p.to_string_lossy();
    n.ends_with(".senna.json") || n.ends_with(run::LUPIN_SUFFIX) || pinto::is_pinto(p)
}

/// `path` read as a marker panel, if it is one: at least one pair, not every
/// label a number (a count table), and, when the run's genes are known, at
/// least one of them (not a log or some other table).
fn read_panel(path: &Path, want: &Want) -> Option<Panel> {
    let Want::Markers(index, _) = want else {
        return None;
    };
    let skip = [
        "parquet", "zarr", "h5", "h5ad", "json", "bam", "bai", "png", "pdf", "log", "zip",
    ];
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    if skip.contains(&ext) || path.metadata().ok()?.len() > MAX_PANEL_BYTES {
        return None;
    }
    let pairs = read_marker_pairs(path.to_str()?).ok()?;
    if pairs.iter().all(|(_, t)| t.parse::<f64>().is_ok()) {
        return None;
    }
    let types: BTreeSet<&str> = pairs.iter().map(|(_, t)| &**t).collect();
    let genes: BTreeSet<&str> = pairs.iter().map(|(g, _)| &**g).collect();
    let matched = index
        .as_ref()
        .map(|ix| genes.iter().filter(|g| ix.match_rows(g).is_some()).count());
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
        let mut entries: Vec<Entry> = std::fs::read_dir(&dir)
            .into_iter()
            .flatten()
            .flatten()
            .filter(|e| !e.file_name().to_string_lossy().starts_with('.'))
            .map(|e| {
                let path = e.path();
                // A `.zarr` store is a directory, but never one to browse into here.
                let dir = path.is_dir() && path.extension().is_none_or(|x| x != "zarr");
                let panel = if dir {
                    None
                } else {
                    read_panel(&path, &self.want)
                };
                Entry {
                    name: e.file_name().to_string_lossy().into_owned(),
                    path,
                    dir,
                    panel,
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
                },
            );
        }
        self.best = entries
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
        self.dir = dir;
        self.entries = entries;
        // Start on the best panel, else the first file that fits.
        let shown = self.shown();
        let at = shown
            .iter()
            .position(|e| Some(&e.path) == self.best.as_ref())
            .or_else(|| shown.iter().position(|e| e.fits(&self.want)))
            .unwrap_or(0);
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
        let all = if self.all {
            "a fitting only"
        } else {
            "a all files"
        };
        (
            format!(" ↑↓ pgup/dn move · enter open/choose · ← up · {all} · ~ home · {cancel} "),
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
            }
            KeyCode::Enter | KeyCode::Right => match self
                .shown()
                .get(at)
                .map(|e| (e.dir, e.path.clone(), e.panel.as_ref().map(|x| x.types)))
            {
                None => {}
                Some((true, dir, _)) => self.open(dir),
                Some((false, _, Some(types))) if types < 2 => {
                    self.note =
                        Some("one cell type: annotation needs types to choose between".into());
                }
                Some((false, file, _)) => return Step::Chosen(file),
            },
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
        loop {
            terminal.draw(|f| p.draw(f))?;
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
