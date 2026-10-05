//! `lupin annotate` without `-o` (and `lupin trajectory` without `-f` and
//! `-o`): run a pass, then go cluster by cluster, giving
//! each a label (a candidate, or any node of the Cell Ontology over the
//! panel), with the genes that set it apart at hand to add as markers.
//!
//! A pass is a child `lupin annotate` ([`runner`]); a round is read back from
//! the manifest it records ([`round`]); saving writes the edits as the next
//! round, as `lupin relabel --next` does.

mod app;
mod export;
mod figure_pane;
mod gallery;
mod genes;
mod menu;
mod ontology;
mod order;
mod picker;
mod round;
mod runner;
mod ui;

use crate::annotate::gene_rows::GeneRows;
use crate::annotate_cmd::AnnotateCliArgs;
use crate::manifest::run::{self, annotated_path, resolve};
use anyhow::{Context, Result};
use app::App;
pub use figure_pane::Graphics;
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

/// `lupin trajectory` without both `--from` and `--out`: the TUI on the
/// order view, which picks the manifest when none is given and asks for the
/// output prefix when it runs.
pub fn run_trajectory(t: &crate::trajectory::run::TrajectoryArgs) -> Result<()> {
    use clap::Parser;
    #[derive(Parser)]
    struct Annotate {
        #[command(flatten)]
        args: AnnotateCliArgs,
    }
    let mut a = Annotate::try_parse_from(["lupin annotate"])
        .context("building the TUI's arguments")?
        .args;
    a.from.clone_from(&t.from);
    a.obo.clone_from(&t.obo);
    a.label_cl.clone_from(&t.label_cl);
    a.graphics = t.graphics;
    let run_with = order::TrajectoryRun {
        argv: t.child_argv(),
        out: t.out.as_deref().map(str::to_string),
        labels: None,
    };
    run(&a, Some(run_with))
}

/// Run the TUI on `args`; with `trajectory`, it opens on the order view and
/// runs `lupin trajectory` with those options.
pub fn run(args: &AnnotateCliArgs, trajectory: Option<order::TrajectoryRun>) -> Result<()> {
    use std::io::IsTerminal;
    anyhow::ensure!(
        std::io::stdin().is_terminal() && std::io::stdout().is_terminal(),
        "no terminal for the TUI: give the output prefix (-o, and -f for trajectory) to run without it"
    );
    let start_in_order = trajectory.is_some();
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
        // The recorded panel, found again if the run moved; one that is gone
        // is picked anew.
        args.markers = match crate::annotate_cmd::recorded_markers(&loaded) {
            Some(p) => p.into_boxed_str(),
            // The order view needs no marker panel: open on the run's types.
            None if start_in_order => Default::default(),
            None => {
                let index = run_genes(&loaded).map(|g| Box::new(GeneRows::build(&g)));
                let n = index.as_deref().map_or(0, GeneRows::n_genes);
                let want = picker::Want::Markers(index, n);
                match picker::pick("Pick a marker panel", &loaded.dir, want)? {
                    Some(p) => p.to_string_lossy().into(),
                    None => return Ok(()),
                }
            }
        };
    }
    // Passes start from the manifest opened. Without `-o`, annotate offers
    // its first round (`.L1`), so reopening the run picks its rounds up
    // again; the trajectory's view offers the next free round under it,
    // never an existing one.
    let out_chosen = !args.out.is_empty();
    if !out_chosen {
        args.out = if start_in_order {
            crate::manifest::family::pass_origin(&loaded.file).1
        } else {
            let stem = run::derive_out_prefix(&loaded.file.to_string_lossy());
            crate::manifest::family::Tag::Round.name(&stem, 1)
        }
        .into_boxed_str();
    }

    eprintln!("lupin: placing the panel on the Cell Ontology…");
    let panel = if args.markers.is_empty() {
        Vec::new()
    } else {
        match crate::annotate::markers::read_panel(&args.markers) {
            Ok(p) => p,
            // The order view needs no panel.
            Err(e) if start_in_order => {
                eprintln!("lupin: no marker panel ({e:#})");
                args.markers = Default::default();
                Vec::new()
            }
            Err(e) => return Err(e),
        }
    };
    let data = crate::manifest::ontology::load(
        Some(&loaded.dir),
        &args.markers,
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
    app.out_chosen = out_chosen;
    app.trajectory = trajectory.unwrap_or_default();
    if !app.args.markers.is_empty() {
        app.original = round::panel_sets(&app.args.markers)?;
    }
    if let Some(cl) = &terms {
        app.panel_ancestry = ontology::type_ancestry(cl, &app.tree);
    }
    app.cl = terms;
    app.hidden = genes::GeneFilter::load(&search)?;
    app.mixed = ontology::Mixed::load(&search)?;
    app.data_search = search;
    // Pick up where an earlier session left this prefix: its latest round.
    // The trajectory TUI starts from the manifest named instead, as the
    // batch run would.
    if target.is_file() && !start_in_order {
        let (latest, _) = crate::manifest::rounds::chain_rounds(&target);
        app.open(&latest);
    }
    if start_in_order {
        // A run that is already annotated shows its round in the cluster
        // panes; the trajectory then reads the same labels.
        if loaded.manifest.annotate.argmax.is_some() {
            app.open(&loaded.file);
        }
        app.focus = app::Focus::Tree;
        app.toggle_order(true);
        // The ordering is this view's main job: start there.
        app.focus = app::Focus::Order;
    }
    // A run with no annotation whose rounds sit beside it: say where they are.
    if app.round.is_none() && crate::manifest::family::has_round(&loaded.file) {
        app.status = "this run has annotated rounds: g lists them".into();
    }

    let mut terminal = ratatui::init();
    app.shift_enter = push_keys();
    let result = (|| -> Result<()> {
        // Drawn when something changed: a key, a resize, the log, the status
        // or a rescoring starting or ending.
        let mut dirty = true;
        while !app.quit {
            let (logged, status) = (app.log.len(), app.status.clone());
            let rescoring = app.rescoring.is_some();
            let progress = app.child.as_ref().map(|r| r.progress.reported_at);
            app.tick();
            dirty |= app.log.len() != logged
                || app.status != status
                || app.rescoring.is_some() != rescoring
                || app.child.as_ref().map(|r| r.progress.reported_at) != progress;
            if dirty {
                terminal.draw(|f| ui::draw(f, &app))?;
                dirty = false;
            }
            if let Some(want) = app.want_file.take() {
                // The file browser takes the screen, then gives it back.
                if app.shift_enter {
                    pop_keys();
                }
                let kind = match want {
                    app::FileWant::Markers => {
                        let index = run_genes(&loaded).map(|g| Box::new(GeneRows::build(&g)));
                        let n = index.as_deref().map_or(0, GeneRows::n_genes);
                        picker::Want::Markers(index, n)
                    }
                    // Tab-separated text files.
                    _ => picker::Want::Labels,
                };
                let picked = picker::pick(want.title(), &loaded.dir, kind)?;
                // A fresh terminal redraws every cell on its first draw.
                terminal = ratatui::init();
                if app.shift_enter {
                    push_keys();
                }
                match want {
                    app::FileWant::Markers => app.set_markers(picked.as_deref()),
                    app::FileWant::Labels => app.set_labels(picked.as_deref()),
                    app::FileWant::Prior => app.set_prior(picked.as_deref()),
                    app::FileWant::LabelCl => app.set_label_cl(picked.as_deref()),
                }
                dirty = true;
            }
            // Every event already waiting is handled before the next draw, so
            // a held key draws only where it ends.
            let mut wait = Duration::from_millis(150);
            while !app.quit && app.want_file.is_none() && event::poll(wait)? {
                wait = Duration::ZERO;
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
    if app.shift_enter {
        pop_keys();
    }
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

/// Ask the terminal to tell Shift+Enter from Enter (the kitty keyboard
/// protocol) when it can; whether it does.
fn push_keys() -> bool {
    use ratatui::crossterm::event::{KeyboardEnhancementFlags, PushKeyboardEnhancementFlags};
    let can = ratatui::crossterm::terminal::supports_keyboard_enhancement().unwrap_or(false);
    can && ratatui::crossterm::execute!(
        std::io::stdout(),
        PushKeyboardEnhancementFlags(KeyboardEnhancementFlags::DISAMBIGUATE_ESCAPE_CODES)
    )
    .is_ok()
}

/// Undo [`push_keys`]; called only after it succeeded.
fn pop_keys() {
    let _ = ratatui::crossterm::execute!(
        std::io::stdout(),
        ratatui::crossterm::event::PopKeyboardEnhancementFlags
    );
}
