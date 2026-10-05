//! Passes, saves and rescoring previews run as child `lupin annotate` / `lupin relabel`, so
//! their progress bars (which hide themselves off a terminal) and logs stay
//! out of the screen, the screen stays live, and stopping is a kill. Their
//! stderr comes back line by line.

use crate::annotate_cmd::AnnotateCliArgs;
use anyhow::{Context, Result};
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::process::{Child, Command, Stdio};
use std::sync::mpsc::{channel, Receiver, Sender};

/// Start `args` as a child pass. Its first round is left at the fine types
/// (`--fine`): labels are chosen here.
/// `overwrite` once the user has agreed to replace the pass's existing round.
pub fn spawn_pass(args: &AnnotateCliArgs, overwrite: bool, log: Sender<String>) -> Result<Child> {
    let args = AnnotateCliArgs {
        fine: true,
        ..args.clone()
    };
    let mut argv = vec!["annotate".to_string()];
    argv.extend(args.to_argv());
    if overwrite {
        argv.push("--overwrite".into());
    }
    spawn(&argv, log)
}

/// Start `lupin trajectory` with `argv` (`-f`, `-o` and the run's options)
/// as a child.
pub fn spawn_trajectory(argv: &[String], log: Sender<String>) -> Result<Child> {
    let mut v = vec!["trajectory".to_string()];
    v.extend_from_slice(argv);
    spawn(&v, log)
}

/// Start `lupin relabel -f <round> -d <decisions> --next` as a child.
pub fn spawn_relabel(round: &Path, decisions: &Path, log: Sender<String>) -> Result<Child> {
    let argv = [
        "relabel".to_string(),
        "-f".into(),
        round.to_string_lossy().into_owned(),
        "-d".into(),
        decisions.to_string_lossy().into_owned(),
        "--next".into(),
    ];
    spawn(&argv, log)
}

/// Start `lupin relabel -f <round> -d <decisions> --preview` as a child:
/// what the decisions would change, rescored, as JSON on the receiver once
/// it is done.
pub fn spawn_preview(
    round: &Path,
    decisions: &Path,
    log: Sender<String>,
) -> Result<(Child, Receiver<String>)> {
    let argv = [
        "relabel".to_string(),
        "-f".into(),
        round.to_string_lossy().into_owned(),
        "-d".into(),
        decisions.to_string_lossy().into_owned(),
        "--preview".into(),
    ];
    let mut child = start(&argv, log, Stdio::piped(), false)?;
    let mut stdout = child.stdout.take().context("no stdout from the child")?;
    let (tx, rx) = channel();
    std::thread::spawn(move || {
        let mut out = String::new();
        if stdout.read_to_string(&mut out).is_ok() {
            let _ = tx.send(out);
        }
    });
    Ok((child, rx))
}

/// Start this `lupin` with `argv`, its stderr sent to `log` line by line.
/// It reports its progress there too, as `@progress` lines (see
/// [`crate::progress`]), for the popup.
fn spawn(argv: &[String], log: Sender<String>) -> Result<Child> {
    start(argv, log, Stdio::null(), true)
}

/// [`spawn`], with the child's stdout as `stdout`; `progress` asks it for
/// `@progress` lines.
fn start(argv: &[String], log: Sender<String>, stdout: Stdio, progress: bool) -> Result<Child> {
    let exe = std::env::current_exe().context("locating the lupin executable")?;
    let mut cmd = Command::new(exe);
    if progress {
        cmd.env(crate::progress::ENV, "1");
    }
    let mut child = cmd
        .args(argv)
        .env("RUST_LOG", "info")
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
    Ok(child)
}
