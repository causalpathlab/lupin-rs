//! `lupin annotate --tui`: run a pass, then go cluster by cluster, giving
//! each a label (a candidate, or any node of the Cell Ontology over the
//! panel), with the genes that set it apart at hand to add as markers.
//!
//! A pass is a child `lupin annotate` ([`runner`]); a round is read back from
//! the manifest it records ([`round`]); saving writes the edits as the next
//! round, as `lupin relabel --next` does.

mod app;
mod export;
mod genes;
mod ontology;
mod picker;
mod round;
mod runner;
mod ui;

use crate::annotate::gene_rows::GeneRows;
use crate::annotate_cmd::AnnotateCliArgs;
use crate::manifest::run::{self, annotated_path, resolve};
use anyhow::Result;
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
    let cwd = std::env::current_dir()?;
    let from = match args.from.as_deref() {
        Some(f) => f.to_string(),
        None => match picker::pick("Pick a run manifest", &cwd, picker::Want::Manifest)? {
            Some(f) => f.to_string_lossy().into_owned(),
            None => return Ok(()),
        },
    };
    let loaded = run::load(&from)?;
    let mut args = args.clone();
    crate::annotate_cmd::apply_level(&mut args, &loaded.file)?;
    // The screen's genes, NES, p and q and live rescoring all come from an
    // enrichment pass; projection writes none of them.
    if args.method == crate::annotate_cmd::AnnotateMethod::Auto {
        args.method = crate::annotate_cmd::AnnotateMethod::Enrichment;
    }
    args.from = Some(loaded.file.to_string_lossy().into());
    if args.markers.is_empty() {
        args.markers = match loaded.manifest.annotate.markers.as_deref() {
            Some(rel) => resolve(&loaded.dir, rel).into_boxed_str(),
            None => {
                let index = run_genes(&loaded).map(|g| GeneRows::build(&g));
                let n = index.as_ref().map_or(0, GeneRows::n_genes);
                let want = picker::Want::Markers(index, n);
                match picker::pick("Pick a marker panel", &loaded.dir, want)? {
                    Some(p) => p.to_string_lossy().into(),
                    None => return Ok(()),
                }
            }
        };
    }
    if args.out.is_empty() {
        // The run's prefix, next to its manifest, one level down: reopening
        // the same run picks up this session's rounds.
        let stem = run::derive_out_prefix(&loaded.file.to_string_lossy());
        args.out = format!("{stem}.L1").into_boxed_str();
        eprintln!("lupin: writing under -o {}", args.out);
    }

    eprintln!("lupin: placing the panel on the Cell Ontology…");
    let panel = crate::annotate::markers::read_panel(&args.markers)?;
    let data = crate::manifest::ontology::load(
        Some(&loaded.dir),
        args.obo.as_deref(),
        args.label_cl.as_deref(),
        crate::manifest::data_files::Fetch::Allowed,
    )?;
    let search = data.search.clone();
    let terms = data.into_terms()?;
    let tree = crate::manifest::ontology::panel_tree_on(terms.as_ref(), &panel);

    let target = annotated_path(&loaded.file, &args.out);
    let mut app = App::new(args, loaded.file.clone(), target.clone(), tree);
    app.fixed_clusters = loaded.manifest.cluster.clusters.is_some();
    app.original = round::panel_sets(&app.args.markers)?;
    if let Some(cl) = &terms {
        app.panel_ancestry = ontology::type_ancestry(cl, &app.tree);
    }
    app.cl = terms;
    app.hidden = genes::GeneFilter::load(&search)?;
    app.mixed = ontology::Mixed::load(&search)?;
    app.data_search = search;
    // Pick up where an earlier session left this prefix: its latest round.
    if target.is_file() {
        let (latest, _) = crate::manifest::rounds::chain_rounds(&target);
        app.open(&latest);
    }

    let mut terminal = ratatui::init();
    let result = (|| -> Result<()> {
        // Drawn when something changed: a key, a resize, the log, the status
        // or a rescoring starting or ending.
        let mut dirty = true;
        while !app.quit {
            let (logged, status) = (app.log.len(), app.status.clone());
            let rescoring = app.rescoring.is_some();
            app.tick();
            dirty |= app.log.len() != logged
                || app.status != status
                || app.rescoring.is_some() != rescoring;
            if dirty {
                terminal.draw(|f| ui::draw(f, &app))?;
                dirty = false;
            }
            if event::poll(Duration::from_millis(150))? {
                match event::read()? {
                    Event::Key(k) if k.kind == KeyEventKind::Press => {
                        app.key(k);
                        dirty = true;
                    }
                    Event::Resize(..) => dirty = true,
                    _ => {}
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

/// The genes the run was trained on, from its dictionary (else its feature
/// embedding), to score marker panels against; `None` when it has neither.
fn run_genes(loaded: &run::Loaded) -> Option<Vec<Box<str>>> {
    use legume_numeric::matrix::parquet::read_parquet_string_column;
    let o = &loaded.manifest.outputs;
    let rel = o.dictionary.as_deref().or(o.feature_embedding.as_deref())?;
    // Only the row-name column: the numbers are not needed.
    read_parquet_string_column(&resolve(&loaded.dir, rel), 0).ok()
}
