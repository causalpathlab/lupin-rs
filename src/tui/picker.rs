//! Choosing the files `--tui` was not given on the command line: a run
//! manifest, then a marker panel, from a browser over the file system.

use crate::annotate::gene_rows::GeneRows;
use crate::annotate::markers::read_marker_pairs;
use crate::manifest::{pinto, run};
use anyhow::Result;
use ratatui::crossterm::event::{self, Event, KeyCode, KeyEventKind};
use ratatui::layout::{Constraint, Layout};
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
    Markers(Option<GeneRows>, usize),
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
        }
    }
}

struct Picker {
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
        let items: Vec<ListItem> = self.shown().into_iter().map(|e| self.row(e)).collect();
        let list = List::new(items)
            .block(Block::default().borders(Borders::ALL))
            .highlight_style(Style::default().add_modifier(Modifier::REVERSED));
        f.render_stateful_widget(list, body, &mut self.state.clone());
        if let Some(note) = &self.note {
            f.render_widget(
                Paragraph::new(format!(" {note}")).fg(Color::LightYellow),
                keys,
            );
            return;
        }
        let all = if self.all {
            "a fitting only"
        } else {
            "a all files"
        };
        f.render_widget(
            Paragraph::new(format!(
                " ↑↓ move · enter open/choose · ← up · {all} · ~ home · q cancel"
            ))
            .dim(),
            keys,
        );
    }
}

/// Browse from `start` for a file of the kind `want`; `None` when cancelled.
pub fn pick(title: &'static str, start: &Path, want: Want) -> Result<Option<PathBuf>> {
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
            p.note = None;
            let at = p.state.selected().unwrap_or(0);
            let (n, at_entry) = {
                let shown = p.shown();
                let e = shown.get(at);
                (
                    shown.len(),
                    e.map(|e| (e.dir, e.path.clone(), e.panel.as_ref().map(|x| x.types))),
                )
            };
            match k.code {
                KeyCode::Char('q') | KeyCode::Esc => return Ok(None),
                KeyCode::Up | KeyCode::Char('k') => p.state.select(Some(at.saturating_sub(1))),
                KeyCode::Down | KeyCode::Char('j') => {
                    p.state.select(Some((at + 1).min(n.saturating_sub(1))));
                }
                KeyCode::PageUp => p.state.select(Some(at.saturating_sub(20))),
                KeyCode::PageDown => p.state.select(Some((at + 20).min(n.saturating_sub(1)))),
                KeyCode::Left | KeyCode::Backspace | KeyCode::Char('h') => {
                    if let Some(up) = p.dir.parent().map(Path::to_path_buf) {
                        p.open(up);
                    }
                }
                KeyCode::Char('~') => {
                    if let Some(home) = std::env::var_os("HOME") {
                        p.open(PathBuf::from(home));
                    }
                }
                KeyCode::Char('a') => {
                    p.all = !p.all;
                    p.state.select(Some(0));
                }
                KeyCode::Enter | KeyCode::Right | KeyCode::Char('l') => match at_entry {
                    None => {}
                    Some((true, dir, _)) => p.open(dir),
                    Some((false, _, Some(types))) if types < 2 => {
                        p.note =
                            Some("one cell type: annotation needs types to choose between".into());
                    }
                    Some((false, file, _)) => return Ok(Some(file)),
                },
                _ => {}
            }
        }
    })();
    ratatui::restore();
    chosen
}
