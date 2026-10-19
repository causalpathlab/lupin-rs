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
    ] {
        let mut argv = vec!["-o", "x", "-m", "m.tsv"];
        argv.extend_from_slice(bad);
        let err = run_annotate(&parse(&argv), &Default::default())
            .unwrap_err()
            .to_string();
        assert!(err.contains("SuSiE"), "{bad:?}: {err}");
    }
}

#[test]
fn projection_leaves_its_method_to_the_run_unless_given_and_the_choice_parses_back() {
    let base = [
        "-f",
        "run.senna.json",
        "-m",
        "m.tsv",
        "-o",
        "out/x",
        "--method",
        "projection",
    ];
    assert_eq!(parse(&base).projection, None);
    let mut ora = base.to_vec();
    ora.extend(["--projection", "ora"]);
    let a = parse(&ora);
    assert_eq!(a.projection, Some(ProjectionMethod::Ora));
    let argv = a.to_argv();
    let argv: Vec<&str> = argv.iter().map(String::as_str).collect();
    assert_eq!(format!("{a:?}"), format!("{:?}", parse(&argv)));
}

#[test]
fn projection_from_embedding_files_runs_the_ora_unless_the_score_is_asked_for() {
    use legume_numeric::matrix::traits::IoOps;
    let dir = tempfile::tempdir().unwrap();
    let (feat, cells, _) = crate::annotate::marker_score::toy::matrices();
    let at = |n: &str| dir.path().join(n).to_string_lossy().into_owned();
    feat.mat
        .to_parquet_with_names(
            &at("feat.parquet"),
            (Some(&feat.rows[..]), Some("gene")),
            None,
        )
        .unwrap();
    cells
        .mat
        .to_parquet_with_names(
            &at("cells.parquet"),
            (Some(&cells.rows[..]), Some("cell")),
            None,
        )
        .unwrap();
    let panel = crate::annotate::marker_score::toy::panel(dir.path());
    let out = at("files");
    let (feat_path, cells_path) = (at("feat.parquet"), at("cells.parquet"));
    let argv = |more: &[&str]| {
        let mut v = vec![
            "-m",
            &panel,
            "-o",
            &out,
            "--method",
            "projection",
            "--feature-embedding",
            &feat_path,
            "--cell-embedding",
            &cells_path,
        ];
        v.extend(more);
        parse(&v)
    };
    assert_eq!(
        resolve_projection(&argv(&[]), None),
        ProjectionMethod::Ora,
        "files say no table"
    );
    run_projection_from_files(
        &argv(&["--projection", "score"]),
        &at("feat.parquet"),
        &at("cells.parquet"),
    )
    .unwrap();
    let table = format!("{out}{}", crate::annotate::outputs::MARKER_SCORE_PARQUET);
    let fields = legume_numeric::matrix::parquet::peek_parquet_field_names(&table).unwrap();
    assert!(fields.iter().any(|f| f.as_ref() == "margin"), "{fields:?}");
}

#[test]
fn the_score_refuses_settings_it_cannot_honour() {
    let base = [
        "-f",
        "run.senna.json",
        "-m",
        "m.tsv",
        "-o",
        "out/x",
        "--projection",
        "score",
    ];
    let with = |more: &[&str]| {
        let mut v = base.to_vec();
        v.extend(more);
        build_score_args(&parse(&v))
    };
    assert!(with(&[]).is_ok());
    let err = with(&["--no-idf"]).expect_err("refused").to_string();
    assert!(err.contains("IDF"), "{err}");
    let err = with(&["--min-markers", "5"])
        .expect_err("refused")
        .to_string();
    assert!(err.contains("min-markers"), "{err}");
}

#[test]
fn the_retired_susie_dispersion_flag_is_still_accepted() {
    // 0.5.0 shipped it; saved commands that pass it must keep working.
    let a = parse(&["-o", "x", "-m", "m.tsv", "--susie-dispersion", "0"]);
    assert!(!a.to_argv().iter().any(|s| s.contains("susie-dispersion")));
}
