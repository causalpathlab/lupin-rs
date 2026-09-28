//! Passes and saves run as child `lupin annotate` / `lupin relabel`, so
//! their progress bars (which hide themselves off a terminal) and logs stay
//! out of the screen, the screen stays live, and stopping is a kill. Their
//! stderr comes back line by line.

use crate::annotate_cmd::AnnotateCliArgs;
use anyhow::{Context, Result};
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::{Child, Command, Stdio};
use std::sync::mpsc::Sender;

/// Start `args` as a child pass. Its first round is left at the fine types
/// (`--fine`): labels are chosen here.
pub fn spawn_pass(args: &AnnotateCliArgs, log: Sender<String>) -> Result<Child> {
    let args = AnnotateCliArgs {
        fine: true,
        tui: false,
        ..args.clone()
    };
    let mut argv = vec!["annotate".to_string()];
    argv.extend(args.to_argv());
    spawn(&argv, log)
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

/// Start this `lupin` with `argv`, its stderr sent to `log` line by line.
fn spawn(argv: &[String], log: Sender<String>) -> Result<Child> {
    let exe = std::env::current_exe().context("locating the lupin executable")?;
    let mut child = Command::new(exe)
        .args(argv)
        .env("RUST_LOG", "info")
        .stdin(Stdio::null())
        .stdout(Stdio::null())
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
