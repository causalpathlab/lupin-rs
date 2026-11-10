//! [`super`] on a made-up round over a small Cell Ontology.

use super::*;
use crate::manifest::rounds::tests::round_with_profile;
use legume_numeric::matrix::parquet::read_parquet_string_columns_by_name;
use std::path::Path;

/// root ─┬─ group1 (class) ─┬─ CT1
///       │                 └─ CT2
///       └─ other
/// CT3, on the panel, has no term.
pub(crate) const OBO: &str = "format-version: 1.2

[Term]
id: CL:0
name: root cell

[Term]
id: CL:1
name: group1
subset: cellxgene_subset
is_a: CL:0

[Term]
id: CL:10
name: CT1
synonym: \"type one\" EXACT []
is_a: CL:1

[Term]
id: CL:12
name: CT2
is_a: CL:1

[Term]
id: CL:2
name: other
is_a: CL:0
";

/// A round (cells a–d in K0, K0, K1, K2 labelled CT1, CT1, CT2, CT3; the
/// panel GENE1 CT1, GENE2 CT2, GENE3 CT3) and the ontology [`OBO`] beside
/// it.
pub(crate) fn round_on_ontology(root: &Path) -> (PathBuf, PathBuf) {
    let round = round_with_profile(root);
    let obo = root.join("r0").join("cl.obo");
    std::fs::write(&obo, OBO).unwrap();
    (round, obo)
}

fn state(root: &Path) -> PanelState {
    let (round, obo) = round_on_ontology(root);
    PanelState::load(&round.to_string_lossy(), Some(&obo.to_string_lossy())).unwrap()
}

#[test]
fn a_round_is_loaded_with_its_panel_on_the_ontology() {
    let dir = tempfile::tempdir().unwrap();
    let s = state(dir.path());
    assert_eq!(s.panel.len(), 3);
    let cl = s.ontology().unwrap();
    assert_eq!(cl.name("CL:10"), Some("CT1"));
    assert!(s.tree.node_of("CT1").is_some());
    assert!(s.original.contains_key("CT1"), "{:?}", s.original.keys());
}

#[test]
fn the_export_writes_each_cells_call_and_each_labels_markers_under_the_prefix_asked() {
    let dir = tempfile::tempdir().unwrap();
    let s = state(dir.path());
    let out = dir.path().join("out").join("x");
    let (cells, markers) = s.export(Some(&out.to_string_lossy())).unwrap();
    assert_eq!(cells, dir.path().join("out/x.cell_annotation.parquet"));
    assert_eq!(markers, dir.path().join("out/x.celltype_markers.tsv"));
    let cols =
        read_parquet_string_columns_by_name(&cells.to_string_lossy(), &["label", "cl_id"]).unwrap();
    let labels: Vec<&str> = cols[0].iter().map(AsRef::as_ref).collect();
    assert_eq!(labels, ["CT1", "CT1", "CT2", "CT3"]);
    assert_eq!(&*cols[1][0], "CL:10", "CT1's term");
    assert_eq!(&*cols[1][3], "", "CT3 has none");
    let tsv = std::fs::read_to_string(&markers).unwrap();
    assert!(tsv.starts_with("celltype\tgene\tsource\tlog2fc\n"), "{tsv}");
    assert!(tsv.contains("CT1\tGENE1\tpanel\t"), "{tsv}");
    let beside: Vec<_> = std::fs::read_dir(dir.path().join("r0"))
        .unwrap()
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.contains("cell_annotation"))
        .collect();
    assert!(beside.is_empty(), "nothing beside the round: {beside:?}");
}

#[test]
fn the_export_goes_beside_the_round_by_default() {
    let dir = tempfile::tempdir().unwrap();
    let s = state(dir.path());
    let (cells, _) = s.export(None).unwrap();
    assert_eq!(cells, dir.path().join("r0/run.cell_annotation.parquet"));
}
