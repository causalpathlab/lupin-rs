//! The coarse first round for [`super`].

use super::*;
use crate::annotate::celltype_tree::TreeSource;
use crate::manifest::rounds::read_argmax;
use crate::manifest::rounds::write_clusters;
use crate::manifest::run::RunManifest;
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::traits::IoOps;

/// root ─┬─ group a ─┬─ CT1
///       │           └─ CT 2
///       └─ group b ─── CT3
const OBO: &str = "data-version: test/1
[Term]
id: CL:9000000
name: root cell
[Term]
id: CL:9000001
name: group a
subset: cellxgene_subset
is_a: CL:9000000
[Term]
id: CL:9000002
name: group b
subset: cellxgene_subset
is_a: CL:9000000
[Term]
id: CL:9000011
name: CT1
is_a: CL:9000001
[Term]
id: CL:9000012
name: CT 2
is_a: CL:9000001
[Term]
id: CL:9000021
name: CT3
is_a: CL:9000002
";

fn setup(root: &Path) -> (String, String) {
    let obo = root.join("cl.obo");
    fs::write(&obo, OBO).unwrap();
    let markers = root.join("markers.tsv");
    fs::write(
        &markers,
        "gene\tcelltype\nGENE1\tCT1\nGENE2\tCT 2\nGENE3\tCT3\n",
    )
    .unwrap();
    (
        obo.to_string_lossy().into_owned(),
        markers.to_string_lossy().into_owned(),
    )
}

/// The fixture ontology under the shipped rules, no aliases.
fn data(obo: &str) -> crate::manifest::data_files::ClData {
    crate::manifest::data_files::ClData {
        rules: crate::annotate::cl_rules::shipped(),
        aliases: crate::annotate::cl_rules::Aliases::default(),
        ontology: Some(obo.into()),
        rule_files: Vec::new(),
        alias_files: Vec::new(),
        search: crate::manifest::data_files::SearchPath::new(None),
        parsed: std::sync::OnceLock::new(),
    }
}

#[test]
fn prepare_groups_the_panel_on_the_ontology_and_maps_its_labels() {
    let root = tempfile::tempdir().unwrap();
    let (obo, markers) = setup(root.path());
    let out = root.path().join("o/run").to_string_lossy().into_owned();
    let p = prepare(&markers, &out, &data(&obo), None).unwrap();
    assert_eq!(p.tree.source, TreeSource::CellOntology);
    let (wired_obo, map_path) = p.ontology.clone().unwrap();
    assert_eq!(wired_obo, obo);
    let map = fs::read_to_string(map_path).unwrap();
    // One spelling per type: the panel's labels are already keyed.
    let ct2: Vec<&str> = map.lines().filter(|l| l.contains("CL:9000012")).collect();
    assert_eq!(ct2, ["CT_2\tCL:9000012"]);
    let tree: TypeTree =
        serde_json::from_str(&fs::read_to_string(format!("{out}.celltype_tree.json")).unwrap())
            .unwrap();
    assert_eq!(tree.group_of("CT_2"), Some("group_a"));

    // The user's own map is kept, not replaced.
    let own = prepare(&markers, &out, &data(&obo), Some("mine.tsv")).unwrap();
    assert!(
        own.ontology.is_none(),
        "the user's map stands; nothing is wired in"
    );
}

#[test]
fn finish_calls_each_cluster_by_its_group_and_keeps_the_fine_labels() {
    let root = tempfile::tempdir().unwrap();
    let (obo, markers) = setup(root.path());
    let out = root.path().join("run").to_string_lossy().into_owned();
    let p = prepare(&markers, &out, &data(&obo), None).unwrap();

    // Cluster 0: CT1 and CT 2 cells; cluster 1: CT3. Q says cluster 0 leans
    // CT3 by the single best type, but group a by the sum.
    let cells: Vec<Box<str>> = ["a", "b", "c", "d"].iter().map(|s| Box::from(*s)).collect();
    write_clusters(
        &format!("{out}.clusters.parquet"),
        &cells,
        &[Some(0), Some(0), Some(1), Some(1)],
    )
    .unwrap();
    let fine: Vec<Option<String>> = ["CT1", "CT_2", "CT3", "CT3"]
        .iter()
        .map(|s| Some((*s).to_string()))
        .collect();
    write_argmax(
        &format!("{out}.argmax.tsv"),
        &cells,
        &fine,
        &[0.9, 0.8, 0.7, 0.6],
    )
    .unwrap();
    let mut q = Mat::zeros(2, 3);
    for (r, row) in [[0.3, 0.3, 0.4], [0.0, 0.1, 0.9]].iter().enumerate() {
        for (c, v) in row.iter().enumerate() {
            q[(r, c)] = *v;
        }
    }
    let rows: Vec<Box<str>> = vec!["K0".into(), "K1".into()];
    let cols: Vec<Box<str>> = vec!["CT1".into(), "CT_2".into(), "CT3".into()];
    q.to_parquet_with_names(
        &format!("{out}.q.parquet"),
        (Some(&rows), Some("cluster")),
        Some(&cols),
    )
    .unwrap();
    let mut m = RunManifest::new(crate::manifest::run::RunKind::Topic, "run");
    m.cluster.clusters = Some("run.clusters.parquet".into());
    m.annotate.argmax = Some("run.argmax.tsv".into());
    m.annotate.cluster_celltype_q = Some("run.q.parquet".into());
    let manifest = root.path().join("run.senna.json");
    m.save(&manifest).unwrap();

    finish(&manifest, &p.tree, true).unwrap();

    let loaded = run::load(&manifest.to_string_lossy()).unwrap();
    let a = &loaded.manifest.annotate;
    let coarse = read_argmax(&resolve(&loaded.dir, a.argmax.as_deref().unwrap())).unwrap();
    assert_eq!(coarse["a"].0, "group_a");
    assert_eq!(coarse["b"].0, "group_a");
    assert_eq!(coarse["c"].0, "CT3", "alone in its class: named for itself");
    assert_eq!(coarse["a"].1, 0.9, "each cell keeps its probability");
    let kept = read_argmax(&resolve(&loaded.dir, a.fine_argmax.as_deref().unwrap())).unwrap();
    assert_eq!(kept["b"].0, "CT_2");
    let summary: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(resolve(&loaded.dir, a.cluster_summary.as_deref().unwrap())).unwrap(),
    )
    .unwrap();
    assert_eq!(summary["0"]["label"], "group_a");
    assert!(a.celltype_tree.is_some());
}

#[test]
fn a_fine_first_round_only_records_the_groups() {
    let root = tempfile::tempdir().unwrap();
    let (obo, markers) = setup(root.path());
    let out = root.path().join("run").to_string_lossy().into_owned();
    let p = prepare(&markers, &out, &data(&obo), None).unwrap();
    let mut m = RunManifest::new(crate::manifest::run::RunKind::Topic, "run");
    m.annotate.argmax = Some("run.argmax.tsv".into());
    let manifest = root.path().join("run.senna.json");
    m.save(&manifest).unwrap();
    finish(&manifest, &p.tree, false).unwrap();
    let loaded = run::load(&manifest.to_string_lossy()).unwrap();
    assert!(loaded.manifest.annotate.celltype_tree.is_some());
    assert!(loaded.manifest.annotate.fine_argmax.is_none());
}

#[test]
fn an_unmatched_panel_type_keeps_the_ontology_walk_off() {
    let root = tempfile::tempdir().unwrap();
    let (obo, _) = setup(root.path());
    let markers = root.path().join("m2.tsv");
    // CT_x matches no term: the walk would refuse the pass, so it is not wired.
    fs::write(
        &markers,
        "gene\tcelltype\nGENE1\tCT1\nGENE2\tCT 2\nGENE9\tCT x\n",
    )
    .unwrap();
    let out = root.path().join("o/run").to_string_lossy().into_owned();
    let p = prepare(&markers.to_string_lossy(), &out, &data(&obo), None).unwrap();
    assert!(p.ontology.is_none());
    assert_eq!(
        p.tree.source,
        TreeSource::CellOntology,
        "the groups still come from the ontology"
    );
}

#[test]
fn without_cluster_probabilities_each_cell_keeps_its_own_groups_call() {
    let root = tempfile::tempdir().unwrap();
    let (obo, markers) = setup(root.path());
    let out = root.path().join("run").to_string_lossy().into_owned();
    let p = prepare(&markers, &out, &data(&obo), None).unwrap();
    // One cluster: a CT1 cell, a CT3 cell and an abstained cell (projection).
    let cells: Vec<Box<str>> = ["a", "b", "c"].iter().map(|s| Box::from(*s)).collect();
    write_clusters(
        &format!("{out}.clusters.parquet"),
        &cells,
        &[Some(0), Some(0), Some(0)],
    )
    .unwrap();
    let fine = vec![Some("CT1".to_string()), Some("CT3".to_string()), None];
    write_argmax(
        &format!("{out}.argmax.tsv"),
        &cells,
        &fine,
        &[0.9, 0.8, 0.1],
    )
    .unwrap();
    let mut m = RunManifest::new(crate::manifest::run::RunKind::Topic, "run");
    m.cluster.clusters = Some("run.clusters.parquet".into());
    m.annotate.argmax = Some("run.argmax.tsv".into());
    let manifest = root.path().join("run.senna.json");
    m.save(&manifest).unwrap();

    finish(&manifest, &p.tree, true).unwrap();

    let loaded = run::load(&manifest.to_string_lossy()).unwrap();
    let coarse = read_argmax(&resolve(
        &loaded.dir,
        loaded.manifest.annotate.argmax.as_deref().unwrap(),
    ))
    .unwrap();
    assert_eq!(coarse["a"].0, "group_a");
    assert_eq!(coarse["b"].0, "CT3");
    assert_eq!(
        coarse["c"].0, "unassigned",
        "an abstained cell stays unassigned"
    );
}
