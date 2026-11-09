use super::*;
use crate::manifest::rounds::tests::{round_with_profile, round_with_profile_rows};

/// The command's arguments for `round` and `value`, the rest its defaults.
fn args(round: &Path, value: HeatValue) -> HeatmapArgs {
    HeatmapArgs {
        from: round.to_string_lossy().into(),
        out: "unused.pdf".into(),
        types: Vec::new(),
        top: 5,
        value,
        source: HeatSource::Observed,
        order: ClusterOrder::Label,
        channel: "spliced".into(),
    }
}

/// The heatmap of `round` with `value`, the rest as the command defaults
/// them, as `adjust` leaves them.
fn heat(round: &Path, value: HeatValue, adjust: impl FnOnce(&mut HeatmapArgs)) -> Result<HeatData> {
    let mut args = args(round, value);
    adjust(&mut args);
    heat_data(&run::load(&args.from).unwrap(), &args)
}

/// Row `gene`'s values.
fn row<'a>(d: &'a HeatData, gene: &str) -> &'a [Option<f32>] {
    &d.rows.iter().find(|r| r.gene == gene).unwrap().values
}

#[test]
fn each_types_markers_are_measured_in_every_cluster() {
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path());
    let d = heat(&round, HeatValue::Log2fc, |_| {}).unwrap();
    let ids: Vec<ClusterId> = d.clusters.iter().map(|c| c.id).collect();
    assert_eq!(ids, [0, 1, 2]);
    assert_eq!(d.clusters[1].label.as_deref(), Some("CT2"));
    let rows: Vec<(&str, &str)> = (d.rows.iter())
        .map(|r| (r.gene.as_str(), r.type_.as_str()))
        .collect();
    assert_eq!(rows, [("GENE1", "CT1"), ("GENE2", "CT2"), ("GENE3", "CT3")]);
    let g1 = row(&d, "GENE1");
    assert!(g1[0].unwrap() > 0.0 && g1[1].unwrap() < 0.0, "{g1:?}");

    let z = heat(&round, HeatValue::Zscore, |_| {}).unwrap();
    let sum: f32 = row(&z, "GENE2").iter().flatten().sum();
    assert!(sum.abs() < 1e-4, "z sums to zero across clusters: {sum}");
    let mean = heat(&round, HeatValue::Mean, |_| {}).unwrap();
    assert_eq!(row(&mean, "GENE1")[0], Some(51.0));
}

#[test]
fn faba_rows_are_read_as_the_panels_genes() {
    let dir = tempfile::tempdir().unwrap();
    let rows = [
        "GENE1/count/spliced",
        "gene2/count/spliced",
        "ENSG3_GENE3",
        "GENE4",
    ];
    let round = round_with_profile_rows(dir.path(), rows);
    let d = heat(&round, HeatValue::Mean, |_| {}).unwrap();
    assert_eq!(d.rows.len(), 3, "{d:?}");
    assert_eq!(d.rows[0].row, "GENE1/count/spliced");
    assert!(d.unmeasured.is_empty());
}

#[test]
fn types_with_no_measured_marker_are_named_and_clusters_follow_the_types() {
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path());
    // CT9 is a panel type whose one marker the profile does not measure.
    std::fs::write(
        dir.path().join("r0").join("markers.tsv"),
        "gene\tcelltype\nGENE1\tCT1\nGENE2\tCT2\nGENE3\tCT3\nGENE9\tCT9\n",
    )
    .unwrap();
    let types = |a: &mut HeatmapArgs| {
        a.types = vec!["CT3".into(), "CT9".into(), "CT1".into()];
    };
    let d = heat(&round, HeatValue::Log2fc, types).unwrap();
    assert_eq!(d.unmeasured, ["CT9"]);
    let ids: Vec<ClusterId> = d.clusters.iter().map(|c| c.id).collect();
    assert_eq!(ids, [2, 0, 1], "CT3's cluster, CT1's, then the rest");
    let by_id = heat(&round, HeatValue::Log2fc, |a| {
        types(a);
        a.order = ClusterOrder::Id;
    })
    .unwrap();
    let ids: Vec<ClusterId> = by_id.clusters.iter().map(|c| c.id).collect();
    assert_eq!(ids, [0, 1, 2]);
}

#[test]
fn the_expected_level_needs_a_topic_model() {
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path());
    let e = heat(&round, HeatValue::Mean, |a| a.source = HeatSource::Expected).unwrap_err();
    assert!(format!("{e:#}").contains("source: observed"), "{e:#}");
}

#[test]
fn the_expected_level_is_the_clusters_topic_shares_times_the_topics() {
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path());
    let r0 = dir.path().join("r0");
    let names = |v: &[&str]| v.iter().map(|s| Box::from(*s)).collect::<Vec<Box<str>>>();
    // Cells a, b (cluster 0) all topic 0; c, d (clusters 1, 2) all topic 1.
    let big = 30.0f32;
    let theta = Mat::from_row_slice(4, 2, &[big, 0.0, big, 0.0, 0.0, big, 0.0, big]);
    let at = |f: &str| r0.join(f).to_string_lossy().into_owned();
    theta
        .to_parquet_with_names(
            &at("run.latent.parquet"),
            (Some(&names(&["a", "b", "c", "d"])), Some("cell")),
            None,
        )
        .unwrap();
    // Topic 0 is GENE1; topic 1 is GENE2 and GENE3 alike.
    let lo = -40.0f32;
    let half = 0.5f32.ln();
    let beta = Mat::from_row_slice(4, 2, &[0.0, lo, lo, half, lo, half, lo, lo]);
    beta.to_parquet_with_names(
        &at("run.dictionary.parquet"),
        (
            Some(&names(&["GENE1", "GENE2", "GENE3", "GENE4"])),
            Some("gene"),
        ),
        None,
    )
    .unwrap();
    let (mut m, _) = run::RunManifest::load(&round).unwrap();
    m.outputs.latent = Some("run.latent.parquet".into());
    m.outputs.softmax_dictionary = Some("run.dictionary.parquet".into());
    m.save(&round).unwrap();
    let d = heat(&round, HeatValue::Mean, |a| a.source = HeatSource::Expected).unwrap();
    let near = |v: Option<f32>, want: f32| (v.unwrap() - want).abs() < 1.0;
    assert!(near(row(&d, "GENE1")[0], 1e4), "{:?}", row(&d, "GENE1"));
    assert!(near(row(&d, "GENE1")[1], 0.0));
    assert!(near(row(&d, "GENE2")[1], 5e3), "{:?}", row(&d, "GENE2"));
}

#[test]
fn the_share_detected_is_counted_from_the_raw_counts() {
    use data_beans::sparse_io::{create_sparse_from_triplets, SparseIoBackend};
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path());
    let r0 = dir.path().join("r0");
    let zarr = r0.join("data.zarr").to_string_lossy().into_owned();
    // GENE1 in a only (cluster 0 has a and b); GENE2 in c (cluster 1).
    let triplets = [(0u64, 0u64, 3.0f32), (1, 2, 1.0), (2, 1, 0.0)];
    let mut data = create_sparse_from_triplets(
        &triplets,
        (4, 4, 3),
        Some(&zarr),
        Some(&SparseIoBackend::Zarr),
    )
    .unwrap();
    let names = |v: &[&str]| v.iter().map(|s| Box::from(*s)).collect::<Vec<Box<str>>>();
    data.register_row_names_vec(&names(&["GENE1", "GENE2", "GENE3", "GENE4"]));
    data.register_column_names_vec(&names(&["a", "b", "c", "d"]));
    drop(data);
    let (mut m, _) = run::RunManifest::load(&round).unwrap();
    m.data.input = vec!["data.zarr".into()];
    m.save(&round).unwrap();
    // GENE1 marks two types, so it has two rows.
    std::fs::write(
        r0.join("markers.tsv"),
        "gene\tcelltype\nGENE1\tCT1\nGENE1\tCT2\nGENE2\tCT2\nGENE3\tCT3\n",
    )
    .unwrap();
    let d = heat(&round, HeatValue::Frac, |_| {}).unwrap();
    let gene1: Vec<_> = d.rows.iter().filter(|r| r.gene == "GENE1").collect();
    assert_eq!(gene1.len(), 2);
    for r in gene1 {
        assert_eq!(r.values, [Some(0.5), Some(0.0), Some(0.0)], "{}", r.type_);
    }
    assert_eq!(row(&d, "GENE2"), [Some(0.0), Some(1.0), Some(0.0)]);
    assert_eq!(
        row(&d, "GENE3"),
        [Some(0.0), Some(0.0), Some(0.0)],
        "a zero is not a count"
    );
}

#[test]
fn the_figure_draws_every_row_and_cluster_and_renders() {
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path());
    let d = heat(&round, HeatValue::Zscore, |_| {}).unwrap();
    let fig = svg(&d);
    for said in ["GENE1 · CT1", "GENE3 · CT3", "K1 · CT2"] {
        assert!(fig.svg.contains(said), "{said}");
    }
    assert!(png(&d).unwrap().starts_with(b"\x89PNG"));
}

#[test]
fn the_command_writes_the_figure_its_name_asks_for_and_its_numbers() {
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path());
    for ext in ["svg", "png", "pdf"] {
        let out = dir
            .path()
            .join(format!("h.{ext}"))
            .to_string_lossy()
            .into_owned();
        run_heatmap(&HeatmapArgs {
            out: out.clone().into(),
            top: 2,
            ..args(&round, HeatValue::Log2fc)
        })
        .unwrap();
        assert!(Path::new(&out).is_file(), "{out}");
    }
    let numbers: HeatData =
        serde_json::from_str(&std::fs::read_to_string(dir.path().join("h.json")).unwrap()).unwrap();
    assert_eq!(numbers.value, HeatValue::Log2fc);
    assert!(figure_base("h.txt").is_err());
}

#[test]
fn a_long_last_cluster_label_fits_inside_the_figure() {
    let label = "a cell type with a rather long name indeed".to_string();
    let d = HeatData {
        round: "r".into(),
        value: HeatValue::Zscore,
        source: HeatSource::Observed,
        clusters: vec![HeatCluster {
            id: 0,
            label: Some(label.clone()),
        }],
        rows: vec![HeatRow {
            gene: "G".into(),
            row: "G".into(),
            type_: "T".into(),
            values: vec![Some(1.0)],
        }],
        unmeasured: Vec::new(),
    };
    let fig = svg(&d);
    let reach = format!("K0 · {label}").chars().count() as f32 * CHAR_W * 0.71;
    let left = ("G · T".chars().count() as f32 * CHAR_W + 12.0).ceil();
    assert!(
        fig.w as f32 >= left + CELL as f32 / 2.0 + reach,
        "{} < {}",
        fig.w,
        left + reach
    );
}

#[test]
fn a_type_the_panel_does_not_have_is_refused_by_name() {
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path());
    let e = heat(&round, HeatValue::Zscore, |a| a.types = vec!["CT1x".into()]).unwrap_err();
    assert!(format!("{e:#}").contains("CT1x"), "{e:#}");
}
