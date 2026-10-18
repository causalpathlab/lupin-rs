//! lupin's TUI as a detachable session, driven through pseudo-terminals
//! standing in for the user's: it starts hosted, a lost client leaves it
//! running, an attach from another terminal repaints it, a TUI started
//! under a session's marker runs plain, `lupin kill` stops it through its
//! own quit path, leaving no process, socket or label behind, and
//! `--no-detach` runs it plain.

#![cfg(unix)]

use legume_tui::testing::{alive, Private, Term};
use legume_tui::Size;
use std::path::Path;
use std::time::Duration;

/// The run browser's title: drawn on every full paint.
const PICKER: &str = "Pick a run manifest";

/// The size every test terminal has.
const SIZE: Size = Size {
    rows: 40,
    cols: 120,
    width: 1200,
    height: 800,
};

/// A private directory for lupin's sessions; a session a failed test left
/// running is killed with it.
fn private() -> Private {
    Private::new(env!("CARGO_BIN_EXE_lupin")).patience(Duration::from_secs(30))
}

/// lupin with `args` on a new terminal, offline, its user files in `home`.
fn term(dir: &Private, home: &Path, args: &[&str], env: &[(&str, &str)]) -> Term {
    let home = home.to_string_lossy().into_owned();
    let mut all = vec![
        ("LUPIN_OFFLINE", "1"),
        ("HOME", home.as_str()),
        ("XDG_CACHE_HOME", home.as_str()),
        ("XDG_DATA_HOME", home.as_str()),
        ("XDG_CONFIG_HOME", home.as_str()),
    ];
    all.extend_from_slice(env);
    dir.term_sized(args, &all, SIZE)
}

#[test]
fn annotate_detaches_attaches_and_is_killed_through_its_quit_path() {
    let dir = private();
    let home = tempfile::tempdir().unwrap();

    // Started in a terminal with no run named, it is hosted, on the run
    // browser: this terminal is its first client.
    let mut first = term(&dir, home.path(), &["annotate"], &[]);
    first.wait_for(PICKER);
    assert!(!first.cooked(), "the client's terminal is raw");
    let listed = dir.wait_listed("lupin-annotate-1\tattached");
    assert!(listed.contains("pick a run manifest"), "{listed:?}");
    let pid = dir.program_pid();
    assert_ne!(
        pid,
        first.child.id(),
        "the program runs apart from its client"
    );

    // The client lost without a word (a dropped connection): the program
    // runs on, detached.
    first.child.kill().unwrap();
    first.child.wait().unwrap();
    drop(first);
    dir.wait_listed("lupin-annotate-1\tdetached");
    assert!(alive(pid));

    // Attached from another terminal: the browser is painted whole there,
    // and keys reach it (`a` lists every file).
    let mut second = term(&dir, home.path(), &["attach"], &[]);
    second.wait_for(PICKER);
    second.type_(b"a");
    second.wait_for("fitting only");
    dir.wait_listed("lupin-annotate-1\tattached");

    // A TUI started under the session's marker, as from a shell inside it,
    // runs plain: no second host, no second socket.
    let mut nested = term(
        &dir,
        home.path(),
        &["trajectory"],
        &[(legume_tui::SESSION_VAR, "lupin-annotate-1")],
    );
    nested.wait_for(PICKER);
    assert_eq!(dir.sockets(), ["lupin-annotate-1"]);
    nested.type_(b"q");
    assert_eq!(nested.exit_code(), 0);
    assert_eq!(dir.sockets(), ["lupin-annotate-1"]);

    // Killed while attached: it quits as q does, the client ends with it,
    // and nothing is left: no process, socket or file.
    let said = dir.run(&["kill"]);
    assert!(said.contains("lupin-annotate-1 stopped"), "{said:?}");
    assert_eq!(second.exit_code(), 0);
    assert!(second.cooked(), "terminal restored after the program ended");
    assert!(dir.sockets().is_empty());
    assert!(dir.is_empty(), "left behind in {:?}", dir.path());
    dir.gone(pid);
    assert!(dir.run(&["sessions"]).contains("no live lupin sessions"));
}

#[test]
fn no_detach_runs_plain_in_this_terminal() {
    let dir = private();
    let home = tempfile::tempdir().unwrap();
    // Opted out, by the flag or the variable: in this process, no client
    // between, no socket, and ctrl+c ends it with the terminal restored.
    for (args, env) in [
        (vec!["annotate", "--no-detach"], vec![]),
        (vec!["trajectory"], vec![(legume_tui::NO_DETACH_VAR, "1")]),
    ] {
        let mut plain = term(&dir, home.path(), &args, &env);
        plain.wait_for(PICKER);
        assert!(!plain.cooked());
        assert!(dir.sockets().is_empty(), "{args:?} {env:?}");
        assert!(!plain.text().contains("\x1b_Ga=d,d=A"), "no client's reset");
        plain.type_(b"\x03");
        assert_eq!(plain.exit_code(), 0);
        assert!(plain.cooked(), "terminal restored");
    }
}

/// A run with nothing annotated yet, for the trajectory's order view: it
/// opens on the main screen without a marker panel.
fn bare_run(dir: &Path) -> String {
    let run = dir.join("Z.senna.json");
    std::fs::write(&run, r#"{"version":2,"kind":"topic","prefix":"Z"}"#).unwrap();
    run.to_string_lossy().into_owned()
}

/// The main screen's header: drawn on every full paint.
const MAIN: &str = "lupin trajectory";

#[test]
fn the_main_screen_comes_back_on_attach_and_stops_on_kill() {
    let dir = private();
    let home = tempfile::tempdir().unwrap();
    let work = tempfile::tempdir().unwrap();
    let run = bare_run(work.path());

    // A run named: a new session, straight to the main screen.
    let mut first = term(&dir, home.path(), &["trajectory", "-f", &run], &[]);
    first.wait_for(MAIN);
    let listed = dir.wait_listed("lupin-trajectory-1\tattached");
    assert!(listed.contains("Z.senna.json"), "{listed:?}");
    let pid = dir.program_pid();

    // Detached by the key: asked for a name, Enter keeps the id.
    first.type_(b"\x1b\x04");
    first.answer(1, b"\r");
    assert_eq!(first.exit_code(), 0);
    assert!(first.cooked(), "terminal restored after detach");
    first.wait_for("Detached lupin-trajectory-1.");
    dir.wait_listed("lupin-trajectory-1\tdetached");
    assert!(alive(pid));

    // Back on another terminal: cleared and painted whole there, and its
    // `?` guide lists the session keys.
    let mut second = term(&dir, home.path(), &["attach", "trajectory-1"], &[]);
    second.wait_for("\x1b[?1049h\x1b[H\x1b[2J");
    let reset = second.text().find("\x1b[?1049h\x1b[H\x1b[2J").unwrap();
    second.wait_for_after(reset, MAIN);
    second.type_(b"?");
    second.wait_for("any other key closes");
    // At its end.
    second.type_(b"\x1b[F");
    second.wait_for("ctrl+alt+d");
    second.wait_for("ctrl+alt+s");
    second.wait_for("pick another session");

    // Killed: it quits as q does, without asking, and leaves nothing.
    let said = dir.run(&["kill", "trajectory-1"]);
    assert!(said.contains("lupin-trajectory-1 stopped"), "{said:?}");
    assert_eq!(second.exit_code(), 0);
    assert!(second.cooked());
    assert!(dir.sockets().is_empty());
    assert!(dir.is_empty(), "left behind in {:?}", dir.path());
    dir.gone(pid);
}

#[test]
fn an_error_while_detached_waits_for_an_attach_to_be_seen() {
    let dir = private();
    let home = tempfile::tempdir().unwrap();
    let work = tempfile::tempdir().unwrap();
    // A run that is read only when written to: the load waits on it.
    let run = work.path().join("Z.senna.json");
    let made = std::process::Command::new("mkfifo")
        .arg(&run)
        .status()
        .unwrap();
    assert!(made.success());
    let run = run.to_string_lossy().into_owned();

    let mut first = term(&dir, home.path(), &["trajectory", "-f", &run], &[]);
    first.wait_for("reading");
    let pid = dir.program_pid();
    first.type_(b"\x1b\x04");
    first.answer(1, b"\r");
    assert_eq!(first.exit_code(), 0);
    dir.wait_listed("lupin-trajectory-1\tdetached");

    // The load fails with nobody there: the session stays, saying why.
    std::fs::write(&run, "not a manifest").unwrap();
    dir.wait_listed("stopped: ");
    assert!(alive(pid));

    // Attached, the error is on screen; a key ends the session, and the
    // error is said on the terminal after it.
    let mut second = term(&dir, home.path(), &["attach", "trajectory-1"], &[]);
    second.wait_for("quits");
    let shown = second.seen_len();
    second.type_(b"q");
    assert_ne!(second.exit_code(), 0);
    assert!(second.cooked());
    second.wait_for_after(shown, "neither a readable run manifest");
    assert!(dir.sockets().is_empty());
    assert!(dir.is_empty(), "left behind in {:?}", dir.path());
    dir.gone(pid);
}
