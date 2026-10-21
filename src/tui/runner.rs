//! Passes, saves and rescoring previews run as child `lupin annotate` / `lupin relabel`, so
//! their progress bars (which hide themselves off a terminal) and logs stay
//! out of the screen, the screen stays live, and stopping is a kill. Their
//! stderr comes back line by line. In a hosted session each runs in a
//! process group of its own: a signal meant for the TUI does not reach it,
//! and the TUI stops it itself (`x`, quitting, `lupin kill`). Each is
//! [`Held`] where a lost host's hangup reaches it too, a save aside.

use crate::annotate_cmd::AnnotateCliArgs;
use crate::child::{lupin_command, relabel_argv};
use anyhow::{Context, Result};
use legume_tui::child::Stopper;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;
use std::process::{Child, ExitStatus, Stdio};
use std::sync::mpsc::{channel, Receiver, Sender};

/// A child the TUI started, held in a [`Stopper`]: stopping it reaches the
/// process group it leads (hosted), and so does a hangup of the TUI (its
/// host lost), unless it saves. How it ended is kept once seen.
pub struct Held {
    stopper: Stopper,
    ended: Option<ExitStatus>,
}

impl Held {
    /// `child`, stopped as `stopper` says ([`Stopper::saving`] for a save).
    pub fn new(child: Child, stopper: Stopper) -> Self {
        stopper.hold(child);
        Self {
            stopper,
            ended: None,
        }
    }

    /// How it ended, once it has; None while it runs.
    pub fn try_wait(&mut self) -> io::Result<Option<ExitStatus>> {
        if self.ended.is_none() {
            self.ended = self.stopper.try_wait()?;
        }
        Ok(self.ended)
    }

    /// Whether it saves ([`Stopper::saving`]): quitting lets it finish.
    pub fn saves(&self) -> bool {
        self.stopper.saves()
    }

    /// Kill it now, with the group it leads.
    pub fn kill(&mut self) {
        self.stopper.kill();
    }

    /// Wait for it to end.
    pub fn wait(&mut self) -> io::Result<ExitStatus> {
        if let Some(st) = self.ended {
            return Ok(st);
        }
        let mut child = self
            .stopper
            .release()
            .ok_or_else(|| io::Error::other("the child is no longer held"))?;
        let st = child.wait()?;
        self.ended = Some(st);
        Ok(st)
    }

    /// Its process id, while it runs.
    #[cfg(test)]
    pub fn id(&self) -> u32 {
        self.stopper.id().expect("a child held")
    }
}

/// A stand-in child, stopped as a pass is.
#[cfg(test)]
impl From<Child> for Held {
    fn from(child: Child) -> Self {
        Self::new(child, Stopper::new())
    }
}

/// Start `args` as a child pass. `overwrite` once the user has agreed to
/// replace the pass's existing round.
pub fn spawn_pass(args: &AnnotateCliArgs, overwrite: bool, log: Sender<String>) -> Result<Held> {
    let mut argv = vec!["annotate".to_string()];
    argv.extend(args.to_argv());
    if overwrite {
        argv.push("--overwrite".into());
    }
    spawn(&argv, log, Stopper::new())
}

/// Start `lupin trajectory` with `argv` (`-f`, `-o` and the run's options)
/// as a child.
pub fn spawn_trajectory(argv: &[String], log: Sender<String>) -> Result<Held> {
    let mut v = vec!["trajectory".to_string()];
    v.extend_from_slice(argv);
    spawn(&v, log, Stopper::new())
}

/// Start `lupin data suggest-cl` as a child: the nearest Cell Ontology terms
/// of `labels` on the ontology run `from` finds (`obo` first), as JSON in
/// `out`.
pub fn spawn_suggest(
    from: &Path,
    obo: Option<&str>,
    labels: &[String],
    out: &Path,
    log: Sender<String>,
) -> Result<Held> {
    let mut argv = vec!["data".to_string(), "suggest-cl".into()];
    if from.is_file() {
        argv.extend(["-f".into(), from.to_string_lossy().into_owned()]);
    }
    if let Some(obo) = obo {
        argv.extend(["--obo".into(), obo.to_string()]);
    }
    argv.extend(["-o".into(), out.to_string_lossy().into_owned(), "--".into()]);
    argv.extend_from_slice(labels);
    spawn(&argv, log, Stopper::new())
}

/// Start `lupin relabel -f <round> -d <decisions> --next` as a child: a
/// save, which quitting lets finish and a hangup does not reach.
pub fn spawn_relabel(round: &Path, decisions: &Path, log: Sender<String>) -> Result<Held> {
    let argv = relabel_argv(round, decisions, &["--next"]);
    spawn(&argv, log, Stopper::saving())
}

/// Start `lupin relabel -f <round> -d <decisions> --preview` as a child:
/// what the decisions would change, rescored, as JSON on the receiver once
/// it is done.
pub fn spawn_preview(
    round: &Path,
    decisions: &Path,
    log: Sender<String>,
) -> Result<(Held, Receiver<String>)> {
    let argv = relabel_argv(round, decisions, &["--preview"]);
    let (tx, rx) = channel();
    let held = start(&argv, log, Some(tx), Stopper::new())?;
    Ok((held, rx))
}

/// Start this `lupin` with `argv`, its stderr sent to `log` line by line.
/// It reports its progress there too, as `@progress` lines (see
/// [`crate::progress`]), for the popup. It is held in `stopper`.
fn spawn(argv: &[String], log: Sender<String>, stopper: Stopper) -> Result<Held> {
    start(argv, log, None, stopper)
}

/// [`spawn`], or, with `out`, a child whose whole stdout is sent there once
/// it is done, and which reports no progress. Every child the TUI starts is
/// started and held here.
fn start(
    argv: &[String],
    log: Sender<String>,
    out: Option<Sender<String>>,
    stopper: Stopper,
) -> Result<Held> {
    let mut cmd = lupin_command(argv)?;
    // With the session's opt-out (`--no-detach`) passed on; hosted, in a
    // process group of its own.
    if let Some(session) = legume_tui::current() {
        session.command(&mut cmd, legume_tui::Group::OwnWhenHosted);
    }
    if out.is_none() {
        cmd.env(crate::progress::ENV, "1");
    }
    let stdout = if out.is_some() {
        Stdio::piped()
    } else {
        Stdio::null()
    };
    let mut child = cmd
        .stdin(Stdio::null())
        .stdout(stdout)
        .stderr(Stdio::piped())
        .spawn()
        .with_context(|| format!("starting lupin {}", argv[0]))?;
    let stderr = child.stderr.take().context("no stderr from the child")?;
    std::thread::spawn(move || {
        for line in BufReader::new(stderr).lines() {
            let Ok(line) = line else { break };
            if log.send(line).is_err() {
                break;
            }
        }
    });
    if let Some(out) = out {
        let mut stdout = child.stdout.take().context("no stdout from the child")?;
        std::thread::spawn(move || {
            let mut text = String::new();
            if stdout.read_to_string(&mut text).is_ok() {
                let _ = out.send(text);
            }
        });
    }
    Ok(Held::new(child, stopper))
}
