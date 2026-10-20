//! `lupin annotate` without `-o` (and `lupin trajectory` without `-f` and
//! `-o`): run a pass, then go cluster by cluster, giving
//! each a label (a candidate, or any node of the Cell Ontology over the
//! panel), with the genes that set it apart at hand to add as markers.
//!
//! A pass is a child `lupin annotate` ([`runner`]); a round is read back from
//! the manifest it records ([`round`]); saving writes the edits as the next
//! round, as `lupin relabel --next` does.
//!
//! The screens run as a detachable session ([`legume_tui`]): the pickers and
//! the main screen share one session and terminal, which outlive a dropped
//! connection; `lupin attach` comes back, and `lupin kill` quits as `q` does.

mod app;
mod ask;
mod export;
mod figure_pane;
mod gallery;
mod genes;
mod menu;
mod ontology;
mod order;
mod picker;
pub(crate) mod round;
mod runner;
mod suggest;
mod ui;

#[cfg(test)]
#[path = "tests/drive.rs"]
mod drive;

use crate::annotate::gene_rows::GeneRows;
use crate::annotate_cmd::AnnotateCliArgs;
use crate::manifest::run::{self, annotated_path, resolve};
use anyhow::{Context, Result};
use app::App;
pub use figure_pane::round_pdf;
pub use legume_tui::graphics::Graphics;
use legume_tui::{cli::SessionArgs, Keyboard, Modes, Mouse, Round, Session};
use ratatui::crossterm::event::{Event, KeyEventKind, KeyboardEnhancementFlags};
use ratatui::DefaultTerminal;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Duration;

/// The tool lupin's TUI sessions belong to: they are named
/// `lupin-annotate-1`, `lupin-trajectory-2`, and `lupin attach` finds them.
pub const TOOL: legume_tui::Tool = legume_tui::Tool::new("lupin");

/// The detach and switch keys for the `?` guide, each with what it does,
/// in lower case as the guide shows keys; none when the screen runs plain.
pub(crate) fn session_keys() -> Vec<(String, String)> {
    legume_tui::help_keys()
        .into_iter()
        .map(|k| (k.key.to_lowercase(), k.long))
        .collect()
}

/// This process's own log lines, for the log pane.
static OWN_LOG: Mutex<Vec<String>> = Mutex::new(Vec::new());

/// The log pane's logger, at its level, and the stderr logger behind the
/// gate, which stays quiet while the screen is up ([`run`]).
struct PaneLogger(log::LevelFilter, &'static legume_tui::logs::Gate);

impl log::Log for PaneLogger {
    fn enabled(&self, m: &log::Metadata) -> bool {
        // The SVG renderer's notes (a font it falls back from) are noise
        // in the pane, unless asked for with -v.
        let svg_note = m.target().starts_with("usvg") && self.0 < log::LevelFilter::Debug;
        m.level() <= self.0 && !svg_note
    }
    fn log(&self, r: &log::Record) {
        if self.enabled(r.metadata()) {
            if let Ok(mut v) = OWN_LOG.lock() {
                v.push(format!("[{}] {}", r.level(), r.args()));
            }
            self.1.log(r);
        }
    }
    fn flush(&self) {
        self.1.flush();
    }
}

/// Send this process's `log` output to the log pane, and to `stderr` (as
/// the other subcommands log) while the screen is not up.
pub fn init_logger(verbose: bool, stderr: env_logger::Logger) {
    let level = if verbose {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    };
    let pane = PaneLogger(level, legume_tui::logs::gate(stderr));
    if log::set_boxed_logger(Box::new(pane)).is_ok() {
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
pub fn run_trajectory(
    t: &crate::trajectory::run::TrajectoryArgs,
    session: &SessionArgs,
) -> Result<()> {
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
    run(&a, Some(order::TrajectoryRun::of(t)), session)
}

/// Run the TUI on `args`; with `trajectory`, it opens on the order view and
/// runs `lupin trajectory` with those options. It runs as a detachable
/// session (`session`): a new one when the run is named (`--from`), else
/// the live ones are offered first.
pub fn run(
    args: &AnnotateCliArgs,
    trajectory: Option<order::TrajectoryRun>,
    session: &SessionArgs,
) -> Result<()> {
    use std::io::IsTerminal;
    anyhow::ensure!(
        std::io::stdin().is_terminal() && std::io::stdout().is_terminal(),
        "no terminal for the TUI: give the output prefix (-o, and -f for trajectory) to run without it"
    );
    // First: hosted, this process only relays the terminal from here on,
    // and the program the host starts does the rest.
    let screen = if trajectory.is_some() {
        "trajectory"
    } else {
        "annotate"
    };
    let session = TOOL.session(screen, &session.options(args.from.is_some()))?;
    // stderr's copy of the log is dropped while the screen is up; the pane
    // keeps every line.
    let quiet = legume_tui::logs::quiet_stderr();
    let mut terminal = session.terminal_with(modes(false))?;
    let latest = screens(&session, &mut terminal, args, trajectory);
    // The latest round is left: said on the terminal attached as the
    // session ends, else kept for whoever looks next.
    match &latest {
        Ok(Some(m)) => session.leave(&format!("lupin: latest round {}", m.display())),
        Ok(None) => {}
        // Seen before the session ends, then said as main says it.
        Err(e) => {
            let _ = session.hold_error(&mut terminal, e);
        }
    }
    session.restore();
    drop(quiet);
    latest.map(|_| ())
}

/// The screen's terminal modes: Shift+Enter told from Enter where the
/// terminal can, and with `mouse` every mouse event reported (the cell map's
/// pointer), which takes the terminal's own text selection away.
fn modes(mouse: bool) -> Modes {
    let keys = Keyboard::IfSupported(KeyboardEnhancementFlags::DISAMBIGUATE_ESCAPE_CODES);
    Modes::new()
        .keyboard_enhancement(keys)
        .mouse(if mouse { Mouse::Full } else { Mouse::Off })
}

/// A message on the whole screen while the TUI loads, and the session's
/// status.
fn say(session: &Session, terminal: &mut DefaultTerminal, msg: &str) -> Result<()> {
    session.status(msg);
    session.draw(terminal, |f| {
        use ratatui::widgets::{Paragraph, Wrap};
        let text = Paragraph::new(format!(" lupin: {msg}")).wrap(Wrap { trim: false });
        f.render_widget(text, f.area());
    })?;
    Ok(())
}

/// The pickers, then the main screen, on `terminal`; the latest round on
/// screen when it is left, if any.
fn screens(
    session: &Session,
    terminal: &mut DefaultTerminal,
    args: &AnnotateCliArgs,
    trajectory: Option<order::TrajectoryRun>,
) -> Result<Option<PathBuf>> {
    let start_in_order = trajectory.is_some();
    let cwd = std::env::current_dir()?;
    let from = match args.from.as_deref() {
        Some(f) => f.to_string(),
        None => {
            let want = picker::Want::Manifest;
            match picker::pick(session, terminal, "Pick a run manifest", &cwd, want)? {
                Some(f) => f.to_string_lossy().into_owned(),
                None => return Ok(None),
            }
        }
    };
    say(session, terminal, "reading the run…")?;
    let loaded = run::load(&from)?;
    let mut args = args.clone();
    crate::annotate_cmd::apply_level(&mut args, &loaded.file)?;
    // The screen's genes, NES, p and q, GO terms and live rescoring all come
    // from an enrichment pass, so auto opens on it (every other method the
    // run supports beside it) unless projection was asked for;
    // projection is a choice in the form.
    use crate::annotate_cmd::{resolve_method, AnnotateMethod};
    let asked_projection = args.projection.is_some() || args.feature_embedding.is_some();
    if args.method == AnnotateMethod::Auto
        && !asked_projection
        && resolve_method(&args, Some(&loaded.manifest)) == AnnotateMethod::Projection
    {
        args.method = AnnotateMethod::Enrichment;
    }
    args.from = Some(loaded.file.to_string_lossy().into());
    // The run's genes, read once for every marker panel browse.
    let mut genes: Option<Option<Arc<GeneRows>>> = None;
    if args.markers.is_empty() {
        // The recorded panel, found again if the run moved; one that is gone
        // is picked anew.
        args.markers = match crate::annotate_cmd::recorded_markers(&loaded) {
            Some(p) => p.into_boxed_str(),
            // The order view needs no marker panel: open on the run's types.
            None if start_in_order => Default::default(),
            None => {
                let index = genes.get_or_insert_with(|| gene_index(&loaded)).clone();
                let n = index.as_deref().map_or(0, GeneRows::n_genes);
                let want = picker::Want::Markers(index, n);
                match picker::pick(session, terminal, "Pick a marker panel", &loaded.dir, want)? {
                    Some(p) => p.to_string_lossy().into(),
                    None => return Ok(None),
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

    say(session, terminal, "placing the panel on the Cell Ontology…")?;
    // Why the order view opens without the recorded panel, for its status.
    let mut no_panel = None;
    let panel = if args.markers.is_empty() {
        Vec::new()
    } else {
        match crate::annotate::markers::read_panel(&args.markers) {
            Ok(p) => p,
            // The order view needs no panel.
            Err(e) if start_in_order => {
                no_panel = Some(format!("no marker panel ({e:#})"));
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
    let target = annotated_path(&loaded.file, &args.out);
    let mut app = App::new(args, loaded.file.clone(), target.clone(), panel, terms);
    app.fixed_clusters = loaded.manifest.cluster.clusters.is_some();
    (app.out_chosen, app.out_typed) = (out_chosen, out_chosen);
    app.trajectory = trajectory.unwrap_or_default();
    if !app.args.markers.is_empty() {
        app.original = round::panel_sets(&app.args.markers)?;
    }
    app.hidden = genes::GeneFilter::load(&search)?;
    app.mixed = ontology::Mixed::load(&search)?;
    app.data_search = search;
    // Annotate picks up where an earlier session left this prefix; else a
    // run that is already annotated shows its round in the cluster panes
    // (the trajectory then reads the same labels).
    app.trajectory_tui = start_in_order;
    app.open_start(loaded.manifest.annotate.argmax.is_some());
    if start_in_order {
        app.focus = app::Focus::Tree;
        app.toggle_order(true);
        // The ordering is this view's main job: start there.
        app.focus = app::Focus::Order;
    }
    // A run with no annotation whose rounds sit beside it: say where they are.
    if app.round.is_none() && crate::manifest::family::has_round(&loaded.file) {
        app.status = "this run has annotated rounds: g lists them".into();
    }
    if let Some(why) = no_panel {
        app.status = why;
    }

    // The screen is the pickers' still: cleared, and the attached terminal
    // asked what it can show, before the main screen is drawn, so no key
    // typed on it meets the query. Taken by the loop's first round.
    session.request_repaint();
    let mut first = Some(session.poll_with(terminal, Duration::ZERO, || {
        fit(session, &mut app);
        Ok(())
    })?);
    // Drawn when something changed: a key, a resize, an attach, the log, the
    // status or a rescoring starting or ending; and each second a job runs.
    // Every frame goes out at once, the figures as last drawn while a new one
    // is rasterised off this thread; the screen is drawn again shortly until
    // it is in (`behind`).
    app.quick_frame.set(true);
    let mut dirty = true;
    let mut behind = false;
    let mut drawn = std::time::Instant::now();
    // The mouse is the screen's only while the cell map is shown.
    let mut mouse = false;
    while !app.quit {
        if app.wants_mouse() != mouse {
            mouse = !mouse;
            session.set_modes(modes(mouse))?;
        }
        let (logged, status) = (app.log.len(), app.status.clone());
        let rescoring = app.rescoring.is_some();
        let progress = app.child.as_ref().map(|r| r.progress.reported_at);
        app.tick();
        dirty |= app.browser.as_mut().is_some_and(|(_, b)| b.poll());
        dirty |= app.log.len() != logged
            || app.status != status
            || app.rescoring.is_some() != rescoring
            || app.child.as_ref().map(|r| r.progress.reported_at) != progress
            || app.redraw_due(drawn.elapsed());
        if dirty || behind {
            app.frame_behind.set(false);
            // Detached, nothing is drawn and the jobs go on.
            session.draw(terminal, |f| ui::draw(f, &app))?;
            behind = app.frame_behind.get();
            dirty = false;
            drawn = std::time::Instant::now();
        }
        session.status(&app.session_status());
        if let Some(want) = app.want_file.take() {
            // The file browser opens in a popup over what asked for it.
            let kind = match want {
                app::FileWant::Markers => {
                    let index = genes.get_or_insert_with(|| gene_index(&loaded)).clone();
                    let n = index.as_deref().map_or(0, GeneRows::n_genes);
                    picker::Want::Markers(index, n)
                }
                // Tab-separated text files.
                _ => picker::Want::Labels,
            };
            let start = match want {
                // Beside the panel already chosen, else the run.
                app::FileWant::Markers if !app.args.markers.is_empty() => {
                    Path::new(app.args.markers.as_ref())
                        .parent()
                        .filter(|d| d.is_dir())
                        .unwrap_or(&loaded.dir)
                        .to_path_buf()
                }
                _ => loaded.dir.clone(),
            };
            app.browser = Some((want, picker::Picker::new(want.title(), &start, kind)));
            dirty = true;
        }
        // Every event already waiting is handled before the next draw, so
        // a held key draws only where it ends.
        // A job that has ended is read before the next key.
        if app.job_ended() {
            continue;
        }
        let mut wait = Duration::from_millis(if behind { 30 } else { 150 });
        while !app.quit && app.want_file.is_none() {
            // After an attach: the screen is cleared and the new terminal
            // asked what it can show, before its modes are on again.
            let round = match first.take() {
                Some(round) => round,
                None => session.poll_with(terminal, wait, || {
                    fit(session, &mut app);
                    Ok(())
                })?,
            };
            wait = Duration::ZERO;
            match round {
                Round::Event(Event::Key(k)) if k.kind == KeyEventKind::Press => {
                    app.key(k);
                    dirty = true;
                }
                Round::Event(Event::Mouse(m)) => dirty |= app.mouse(m),
                Round::Event(Event::Resize(..)) => dirty = true,
                Round::Event(_) => {}
                Round::Repainted => {
                    app.shift_enter = session.keyboard_enhanced();
                    dirty = true;
                    break;
                }
                // `lupin kill` quits as q does, without asking.
                Round::Stopped => {
                    if app.saving_now() {
                        app.status = "lupin kill: the save finishes, then lupin quits".into();
                        session.draw(terminal, |f| ui::draw(f, &app))?;
                    }
                    app.kill();
                    break;
                }
                _ => break,
            }
        }
    }
    Ok(app.round.as_ref().map(|r| r.manifest.clone()))
}

/// Fit the figures to the terminal attached now: its picture protocol as
/// `--graphics` names it, else asked of the terminal, else guessed from its
/// environment. Asked with the terminal raw and before mouse reports are on.
fn fit(session: &Session, app: &mut App) {
    use legume_tui::graphics::{reprobe, Fallback};
    let old = app.picker().clone();
    app.repick(reprobe(
        session,
        &old,
        app.args.graphics,
        Fallback::ProbeThenEnv,
    ));
}

/// The run's genes, indexed to score marker panels against.
fn gene_index(loaded: &run::Loaded) -> Option<Arc<GeneRows>> {
    run_genes(loaded).map(|g| Arc::new(GeneRows::build(&g)))
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
