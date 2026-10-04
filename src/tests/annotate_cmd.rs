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
