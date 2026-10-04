//! The TUI's passes re-run as child processes: [`super::AnnotateCliArgs::to_argv`].

use super::*;

#[derive(clap::Parser)]
struct Cli {
    #[command(flatten)]
    annotate: AnnotateCliArgs,
}

fn parse(argv: &[&str]) -> AnnotateCliArgs {
    let argv = std::iter::once("lupin").chain(argv.iter().copied());
    <Cli as clap::Parser>::try_parse_from(argv)
        .unwrap()
        .annotate
}

#[test]
fn the_argv_parses_back_to_the_same_arguments() {
    let a = parse(&[
        "-f",
        "run.senna.json",
        "-m",
        "m.tsv.gz",
        "-o",
        "out/x",
        "--method",
        "enrichment",
        "--knn",
        "20",
        "--resolution",
        "0.7",
        "--num-clusters",
        "12",
        "--cluster-seed",
        "3",
        "--obo",
        "cl.obo",
        "--fine",
        "--no-idf",
    ]);
    let argv = a.to_argv();
    let argv: Vec<&str> = argv.iter().map(String::as_str).collect();
    let b = parse(&argv);
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
}

#[test]
fn a_moved_runs_marker_panel_is_found_by_its_path_tail() {
    let tmp = tempfile::tempdir().unwrap();
    let proj = tmp.path().join("proj");
    std::fs::create_dir_all(proj.join("run")).unwrap();
    std::fs::create_dir_all(proj.join("data")).unwrap();
    std::fs::write(proj.join("data/panel.tsv"), "GENE1\tCT1\n").unwrap();
    let manifest = proj.join("run/x.senna.json");
    std::fs::write(
        &manifest,
        r#"{"version":2,"kind":"topic","prefix":"/elsewhere/proj/run/x",
            "annotate":{"markers":"/elsewhere/proj/data/panel.tsv"}}"#,
    )
    .unwrap();
    let loaded = crate::manifest::run::load(&manifest.to_string_lossy()).unwrap();
    let found = crate::annotate_cmd::recorded_markers(&loaded).expect("found again");
    assert!(
        std::path::Path::new(&found).ends_with("data/panel.tsv"),
        "{found}"
    );
}
