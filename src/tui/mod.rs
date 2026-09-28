//! `lupin annotate --tui`: run a pass, then go cluster by cluster, giving
//! each a label (a candidate, or any node of the Cell Ontology over the
//! panel), with the genes that set it apart at hand to add as markers.
//!
//! A pass is a child `lupin annotate` ([`runner`]); a round is read back from
//! the manifest it records ([`round`]); saving writes the edits as the next
//! round, as `lupin relabel --next` does.

mod app;
mod export;
mod ontology;
mod round;
mod runner;
mod ui;

use crate::annotate::markers::read_marker_pairs;
use crate::annotate_cmd::AnnotateCliArgs;
use crate::manifest::run::{self, annotated_path, resolve};
use anyhow::{Context, Result};
use app::App;
use ratatui::crossterm::event::{self, Event, KeyEventKind};
use std::sync::Mutex;
use std::time::Duration;

/// This process's own log lines while the screen is up.
static OWN_LOG: Mutex<Vec<String>> = Mutex::new(Vec::new());

struct PaneLogger(log::LevelFilter);

impl log::Log for PaneLogger {
    fn enabled(&self, m: &log::Metadata) -> bool {
        m.level() <= self.0
    }
    fn log(&self, r: &log::Record) {
        if self.enabled(r.metadata()) {
            if let Ok(mut v) = OWN_LOG.lock() {
                v.push(format!("[{}] {}", r.level(), r.args()));
            }
        }
    }
    fn flush(&self) {}
}

/// Send this process's `log` output to the log pane instead of stderr.
pub fn init_logger(verbose: bool) {
    let level = if verbose {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    };
    if log::set_boxed_logger(Box::new(PaneLogger(level))).is_ok() {
        log::set_max_level(level);
    }
}

fn drain_own_log() -> Vec<String> {
    OWN_LOG
        .lock()
        .map(|mut v| std::mem::take(&mut *v))
        .unwrap_or_default()
}

pub fn run(args: &AnnotateCliArgs) -> Result<()> {
    let from = args
        .from
        .as_deref()
        .context("--tui needs a run manifest (-f)")?;
    let loaded = run::load(from)?;
    let mut args = args.clone();
    args.from = Some(loaded.file.to_string_lossy().into());
    if args.markers.is_empty() {
        let rel = loaded
            .manifest
            .annotate
            .markers
            .as_deref()
            .context("--tui needs a marker panel (-m)")?;
        args.markers = resolve(&loaded.dir, rel).into_boxed_str();
    }

    eprintln!("lupin: placing the panel on the Cell Ontology…");
    let panel: Vec<(String, String)> = read_marker_pairs(&args.markers)?
        .into_iter()
        .map(|(g, t)| (g.into_string(), t.into_string()))
        .collect();
    let data = crate::manifest::ontology::load(
        Some(&loaded.dir),
        args.obo.as_deref(),
        args.label_cl.as_deref(),
    )?;
    let terms = data.terms()?;
    let tree = crate::manifest::ontology::panel_tree_on(terms.as_ref(), &panel);

    let target = annotated_path(&loaded.file, &args.out);
    let mut app = App::new(args, loaded.file.clone(), target.clone(), tree);
    app.fixed_clusters = loaded.manifest.cluster.clusters.is_some();
    app.original = round::panel_sets(&app.args.markers)?;
    if let Some(cl) = &terms {
        app.panel_ancestry = ontology::type_ancestry(cl, &app.tree);
    }
    app.cl = terms;
    app.data_search = data.search.clone();
    // Pick up where an earlier session left this prefix: its latest round.
    if target.is_file() {
        let (latest, _) = crate::manifest::rounds::chain_rounds(&target);
        app.open(&latest);
    }

    let mut terminal = ratatui::init();
    let result = (|| -> Result<()> {
        while !app.quit {
            app.tick();
            terminal.draw(|f| ui::draw(f, &app))?;
            if event::poll(Duration::from_millis(150))? {
                if let Event::Key(k) = event::read()? {
                    if k.kind == KeyEventKind::Press {
                        app.key(k);
                    }
                }
            }
        }
        Ok(())
    })();
    ratatui::restore();
    if let Some(r) = &app.round {
        eprintln!("lupin: latest round {}", r.manifest.display());
    }
    result
}
