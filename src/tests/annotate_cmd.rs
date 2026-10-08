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
fn the_susie_stage_parses_back_to_the_same_arguments() {
    let a = parse(&[
        "-f",
        "run.senna.json",
        "-m",
        "m.tsv",
        "-o",
        "out/x",
        "--clusters",
        "c.parquet",
        "--susie-effects",
        "3",
        "--susie-dispersion",
        "0",
        "--mcmc-chains",
        "2",
        "--mcmc-samples",
        "200",
        "--mcmc-warmup",
        "100",
        "--mcmc-thin",
        "2",
    ]);
    assert!(build_enrichment_args(&parse(&["-o", "x", "--no-susie"]))
        .susie
        .is_none());
    assert!(build_enrichment_args(&a).susie.is_some());
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

#[test]
fn susie_samples_a_thousand_draws_per_chain_by_default() {
    let cfg = build_enrichment_args(&parse(&["-o", "x"])).susie.unwrap();
    assert_eq!(cfg.samples, 1000);
    assert_eq!(crate::annotate::susie::SusieConfig::default().samples, 1000);
}

#[test]
fn susie_settings_that_leave_no_posterior_fail_the_command() {
    for bad in [
        &["--susie-effects", "0"][..],
        &["--mcmc-chains", "0"],
        &["--mcmc-samples", "0"],
        &["--mcmc-thin", "0"],
        &["--susie-dispersion=-1"],
    ] {
        let mut argv = vec!["-o", "x", "-m", "m.tsv"];
        argv.extend_from_slice(bad);
        let err = run_annotate(&parse(&argv)).unwrap_err().to_string();
        assert!(err.contains("SuSiE"), "{bad:?}: {err}");
    }
}
