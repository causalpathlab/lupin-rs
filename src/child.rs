//! This `lupin` run as a child: the TUI's passes and saves and `lupin mcp`'s
//! tools start their work this way, so it can be followed line by line and
//! stopped with a kill.

use anyhow::Result;
use std::path::Path;
use std::process::Command;

/// `lupin` with `argv`, logging at `info`: the binary this process runs,
/// even when the file was replaced since it started (see
/// [`legume_mcp::own_command`]).
pub fn lupin_command(argv: &[String]) -> Result<Command> {
    let mut cmd = legume_mcp::own_command(argv)?;
    cmd.env("RUST_LOG", "info");
    Ok(cmd)
}

/// `relabel -f <round> -d <decisions>` and then `tail` (`--next`,
/// `--preview`, …).
pub fn relabel_argv(round: &Path, decisions: &Path, tail: &[&str]) -> Vec<String> {
    let mut v = vec![
        "relabel".to_string(),
        "-f".into(),
        round.to_string_lossy().into_owned(),
        "-d".into(),
        decisions.to_string_lossy().into_owned(),
    ];
    v.extend(tail.iter().map(|t| t.to_string()));
    v
}
