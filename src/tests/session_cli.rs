//! The session flags and subcommands on lupin's command line, and the help
//! that describes them.

use super::*;
use clap::CommandFactory;

fn parse(argv: &[&str]) -> Cli {
    Cli::try_parse_from(std::iter::once("lupin").chain(argv.iter().copied())).unwrap()
}

/// The long help of subcommand `name`.
fn long_help(name: &str) -> String {
    let mut cli = Cli::command();
    cli.build();
    cli.find_subcommand_mut(name)
        .unwrap()
        .render_long_help()
        .to_string()
}

#[test]
fn the_session_flags_sit_beside_annotate_and_stay_out_of_its_child_passes() {
    let cli = parse(&["annotate", "--no-detach", "--session", "work", "-o", "x"]);
    let Commands::Annotate { args, session, .. } = cli.cmd else {
        panic!("annotate");
    };
    assert!(session.no_detach);
    assert_eq!(session.session.as_deref(), Some("work"));
    // A pass started from the TUI is `lupin annotate` with these: no
    // session flag goes with them.
    let argv = args.to_argv();
    assert!(!argv
        .iter()
        .any(|a| a.contains("detach") || a == "--session"));
}

#[test]
fn trajectory_takes_the_session_flags_too() {
    let cli = parse(&["trajectory", "-f", "run", "--no-detach"]);
    let Commands::Trajectory { args, session, .. } = cli.cmd else {
        panic!("trajectory");
    };
    assert!(session.no_detach && session.session.is_none());
    assert_eq!(args.from.as_deref(), Some("run"));
    assert!(!args.child_argv().iter().any(|a| a.contains("detach")));
}

#[test]
fn the_session_subcommands_parse() {
    assert!(matches!(parse(&["sessions"]).cmd, Commands::Sessions));
    let Commands::Attach(a) = parse(&["attach", "work"]).cmd else {
        panic!("attach");
    };
    assert_eq!(a.name.as_deref(), Some("work"));
    let Commands::Kill(k) = parse(&["kill", "--force"]).cmd else {
        panic!("kill");
    };
    assert!(k.force && k.name.is_none());
    let Commands::Rename(r) = parse(&["rename", "annotate-1", "work"]).cmd else {
        panic!("rename");
    };
    assert_eq!(
        (r.session.as_str(), r.label.as_str()),
        ("annotate-1", "work")
    );
}

#[test]
fn the_tui_commands_help_ends_with_their_session_in_lupins_words() {
    for (name, goes_on, resume) in [
        (
            "annotate",
            "a pass or a save goes on",
            "`lupin annotate` with no --from",
        ),
        (
            "trajectory",
            "a pass or a trajectory run goes on",
            "`lupin trajectory` with no --from",
        ),
    ] {
        let help = long_help(name);
        // Wrapped to lines, so read as one.
        let flat = help.split_whitespace().collect::<Vec<_>>().join(" ");
        for said in [
            goes_on,
            resume,
            "ctrl+alt+d detaches",
            "ctrl+alt+s",
            "`lupin attach`",
            "--no-detach",
            "`?` in the TUI shows the keys",
        ] {
            assert!(flat.contains(said), "{name}: {said:?} in {help}");
        }
    }
    let kill = long_help("kill");
    let flat = kill.split_whitespace().collect::<Vec<_>>().join(" ");
    assert!(flat.contains("unsaved edits are dropped"), "{kill}");
}

#[test]
fn annotate_with_pdf_draws_the_rounds_map_without_a_terminal() {
    let dir = tempfile::tempdir().unwrap();
    let loaded = crate::figures::tests::round(dir.path(), true);
    let out = dir.path().join("map.pdf");
    let (from, to) = (loaded.file.to_string_lossy(), out.to_string_lossy());
    let cli = parse(&["annotate", "-f", &from, "--pdf", &to, "--width", "3.5"]);
    assert!(!cli.cmd.opens_tui(), "no screen for a PDF");
    let Commands::Annotate {
        args, session, pdf, ..
    } = cli.cmd
    else {
        panic!("annotate");
    };
    assert_eq!(
        (pdf.width, pdf.dpi),
        (3.5, 300),
        "senna's flags and default"
    );
    // A pass started from the TUI takes none of them.
    assert!(!args.to_argv().iter().any(|a| a == "--pdf" || a == "--dpi"));
    annotate_command(&args, &session, &pdf).unwrap();
    assert!(std::fs::read(&out).unwrap().starts_with(b"%PDF-"));

    // A PDF is of a round: one to draw, and no pass beside it.
    let cli = parse(&["annotate", "--pdf", &to]);
    let Commands::Annotate {
        args, session, pdf, ..
    } = cli.cmd
    else {
        panic!("annotate");
    };
    assert!(annotate_command(&args, &session, &pdf).is_err(), "no -f");
    let cli = parse(&["annotate", "-f", &from, "-o", "x", "--pdf", &to]);
    let Commands::Annotate {
        args, session, pdf, ..
    } = cli.cmd
    else {
        panic!("annotate");
    };
    assert!(annotate_command(&args, &session, &pdf).is_err(), "with -o");
}

#[test]
fn a_job_opens_the_tui_on_the_step_given_whole() {
    let cli = parse(&["annotate", "--job", "--from", "r.senna.json", "--out", "p"]);
    assert!(cli.cmd.opens_tui(), "-o given, and still the screen");
    let Commands::Annotate { job, args, .. } = cli.cmd else {
        panic!("not annotate");
    };
    assert!(job);
    assert_eq!(args.out.as_ref(), "p");
    let cli = parse(&["trajectory", "--job", "-f", "r.senna.json", "-o", "t"]);
    assert!(cli.cmd.opens_tui());
    // Without it, a step given whole runs without the screen, as before.
    let cli = parse(&["annotate", "--from", "r.senna.json", "--out", "p"]);
    assert!(!cli.cmd.opens_tui());
}
