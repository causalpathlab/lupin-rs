//! [`super`]'s calls and marker table on a made-up round.

use super::*;
use crate::annotate::celltype_tree::{TreeSource, TypeGroup, TypeTree};
use crate::tui::round::{Candidate, ClusterView};
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::traits::MatWithNames;

/// `group1` over `CT1` and `CT2`.
fn tree() -> PanelTree {
    PanelTree::from_type_tree(&TypeTree {
        source: TreeSource::MarkerSharing,
        release: None,
        label_cl: BTreeMap::new(),
        groups: vec![TypeGroup {
            name: "group1".into(),
            cl_id: None,
            members: vec!["CT1".into(), "CT2".into()],
        }],
    })
}

fn cluster(id: ClusterId, label: &str, shares: &[(&str, f32)]) -> ClusterView {
    let shares: Vec<(String, f32)> = shares.iter().map(|(t, s)| ((*t).into(), *s)).collect();
    ClusterView {
        id,
        cells: 1,
        label: Some(label.into()),
        candidates: shares
            .iter()
            .map(|(label, share)| Candidate {
                label: label.clone(),
                share: *share,
                susie: None,
                nes: None,
                p: None,
                q: None,
            })
            .collect(),
        shares,
        genes: Vec::new(),
        terms: Vec::new(),
    }
}

/// K0 is CT1, K1 is CT2; GENE1 is up in K0, GENE3 in K1, NEW in K0 only.
fn round() -> RoundView {
    let set = |g: &[&str]| g.iter().map(|g| (*g).to_string()).collect::<BTreeSet<_>>();
    RoundView {
        manifest: PathBuf::new(),
        clusters: vec![
            cluster(0, "CT1", &[("CT1", 0.7), ("CT2", 0.2)]),
            cluster(1, "CT2", &[("CT2", 0.6), ("CT1", 0.1)]),
        ],
        cell_names: vec!["c0".into(), "c1".into()],
        cell_clusters: vec![Some(0), Some(1)],
        expression: Some(crate::tui::round::Expression::new(MatWithNames {
            rows: vec!["GENE1".into(), "GENE3".into(), "NEW".into(), "GENE4".into()],
            cols: vec!["K0".into(), "K1".into()],
            mat: Mat::from_row_slice(4, 2, &[8.0, 1.0, 1.0, 8.0, 12.0, 0.0, 20.0, 20.0]),
        })),
        markers: BTreeMap::from([
            ("CT1".into(), set(&["GENE1", "GENE2"])),
            ("CT2".into(), set(&["GENE3"])),
        ]),
        loose_cells: 0,
        decided: BTreeSet::new(),
        rescorable: false,
    }
}

fn relabel(cluster: ClusterId, label: &str) -> Edit {
    Edit::Label {
        cluster,
        label: label.into(),
        reason: "why".into(),
    }
}

#[test]
fn a_class_label_pools_its_types_share_and_names_its_lineage() {
    let calls = cluster_calls(
        &round(),
        &tree(),
        None,
        &crate::tui::ontology::Mixed::default(),
        &[relabel(1, "group1")],
    );
    assert_eq!(calls[&0].label, "CT1");
    assert_eq!(calls[&0].lineage, "group1 > CT1");
    assert!((calls[&0].share - 0.7).abs() < 1e-6);
    assert_eq!(calls[&1].label, "group1");
    assert_eq!(calls[&1].lineage, "group1");
    assert!((calls[&1].share - 0.7).abs() < 1e-6, "CT2 0.6 + CT1 0.1");
    assert_eq!(calls[&1].top, "CT2");
}

#[test]
fn markers_are_panel_or_added_and_the_rest_suggested() {
    let r = round();
    let t = tree();
    let calls = cluster_calls(&r, &t, None, &crate::tui::ontology::Mixed::default(), &[]);
    // CD2 was added this session: not in the original panel.
    let original = BTreeMap::from([
        ("CT1".to_string(), BTreeSet::from(["GENE1".to_string()])),
        ("CT2".to_string(), BTreeSet::from(["GENE3".to_string()])),
    ]);
    let rows = marker_rows(
        &r,
        &t,
        &crate::tui::ontology::Mixed::default(),
        &calls,
        &original,
    );
    let of = |ct: &str| -> Vec<(String, &str)> {
        rows.iter()
            .filter(|m| m.celltype == ct)
            .map(|m| (m.gene.clone(), m.source))
            .collect()
    };
    assert_eq!(
        of("CT1"),
        [
            ("GENE1".to_string(), "panel"),
            ("GENE2".to_string(), "added"),
            ("NEW".to_string(), "suggested")
        ]
    );
    let gene1 = rows.iter().find(|m| m.gene == "GENE1").unwrap();
    assert!(gene1.log2fc > 0.5, "up in the CT1 cluster");
    assert!(
        rows.iter()
            .find(|m| m.gene == "GENE2")
            .unwrap()
            .log2fc
            .is_nan(),
        "not measured"
    );
    assert!(
        !rows.iter().any(|m| m.gene == "GENE4"),
        "high but everywhere: not suggested"
    );
    assert_eq!(of("CT2"), [("GENE3".to_string(), "panel")]);
}

#[test]
fn a_term_off_the_panel_gets_its_id_and_lineage_from_the_ontology() {
    let obo = "[Term]\nid: CL:0\nname: root cell\n\n\
               [Term]\nid: CL:1\nname: group1\nsubset: cellxgene_subset\nis_a: CL:0\n\n\
               [Term]\nid: CL:10\nname: CT1\nis_a: CL:1\n\n\
               [Term]\nid: CL:12\nname: CT2\nis_a: CL:1\n\n\
               [Term]\nid: CL:3\nname: CT3 cell\nis_a: CL:0\n";
    let cl = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let panel: Vec<(String, String)> = [("G1", "CT1"), ("G2", "CT2")]
        .iter()
        .map(|(g, t)| ((*g).to_string(), (*t).to_string()))
        .collect();
    let t = PanelTree::from_ontology(&cl, &panel).unwrap();
    let calls = cluster_calls(
        &round(),
        &t,
        Some(&cl),
        &crate::tui::ontology::Mixed::default(),
        &[relabel(0, "CT3_cell"), relabel(1, "group1")],
    );
    assert_eq!(calls[&0].cl_id, "CL:3");
    assert_eq!(calls[&0].lineage, "root_cell > CT3_cell");
    assert_eq!(calls[&1].cl_id, "CL:1", "a class the panel's tree has");
    assert!(
        (calls[&1].share - 0.7).abs() < 1e-6,
        "CT2 0.6 + CT1 0.1, both under group1"
    );
}

#[test]
fn a_mixed_label_sums_its_parts_and_sits_under_their_common_ancestor() {
    let obo = "[Term]\nid: CL:0\nname: root cell\n\n\
               [Term]\nid: CL:1\nname: group1\nsubset: cellxgene_subset\nis_a: CL:0\n\n\
               [Term]\nid: CL:10\nname: CT1\nis_a: CL:1\n\n\
               [Term]\nid: CL:12\nname: CT2\nis_a: CL:1\n\n\
               [Term]\nid: CL:3\nname: CT3 cell\nis_a: CL:0\n";
    let cl = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let panel: Vec<(String, String)> = [("G1", "CT1"), ("G2", "CT2")]
        .iter()
        .map(|(g, t)| ((*g).to_string(), (*t).to_string()))
        .collect();
    let t = PanelTree::from_ontology(&cl, &panel).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let mut mixed = crate::tui::ontology::Mixed::default();
    mixed
        .add(
            "CT2+CT1",
            &["CT2".into(), "CT1".into()],
            &dir.path().join("m.tsv"),
        )
        .unwrap();
    let calls = cluster_calls(&round(), &t, Some(&cl), &mixed, &[relabel(0, "CT2+CT1")]);
    assert_eq!(calls[&0].label, "CT2+CT1");
    assert!((calls[&0].share - 0.9).abs() < 1e-6, "CT1 0.7 + CT2 0.2");
    assert_eq!(calls[&0].cl_id, "CL:12|CL:10");
    assert_eq!(calls[&0].lineage, "root_cell > group1 > CT2+CT1");
    // Without the ontology: the parts on the panel's tree still add up.
    let alone = cluster_calls(&round(), &tree(), None, &mixed, &[relabel(0, "CT2+CT1")]);
    assert!((alone[&0].share - 0.9).abs() < 1e-6);
    // Unlisted, `CT2+CT1` is just a label off the panel.
    let plain = cluster_calls(
        &round(),
        &tree(),
        None,
        &crate::tui::ontology::Mixed::default(),
        &[relabel(0, "CT2+CT1")],
    );
    assert!(plain[&0].share < 0.9);
}

#[test]
fn a_mixed_labels_markers_are_its_parts_markers() {
    let r = round();
    let t = tree();
    let dir = tempfile::tempdir().unwrap();
    let mut mixed = crate::tui::ontology::Mixed::default();
    mixed
        .add(
            "mix1",
            &["CT2".into(), "CT1".into()],
            &dir.path().join("m.tsv"),
        )
        .unwrap();
    let calls = cluster_calls(&r, &t, None, &mixed, &[relabel(0, "mix1")]);
    let rows = marker_rows(&r, &t, &mixed, &calls, &r.markers);
    let genes: BTreeSet<&str> = rows
        .iter()
        .filter(|m| m.celltype == "mix1")
        .map(|m| m.gene.as_str())
        .collect();
    for g in ["GENE1", "GENE2", "GENE3"] {
        assert!(genes.contains(g), "{g} in {genes:?}");
    }
}

#[test]
fn unmeasured_markers_come_last() {
    let mut r = round();
    // CD2 sorts before GENE1 by name but is not measured.
    r.markers.insert(
        "CT1".into(),
        ["GENE2", "GENE1", "ZZZ"].map(String::from).into(),
    );
    let t = tree();
    let calls = cluster_calls(&r, &t, None, &crate::tui::ontology::Mixed::default(), &[]);
    let rows = marker_rows(
        &r,
        &t,
        &crate::tui::ontology::Mixed::default(),
        &calls,
        &r.markers,
    );
    let ct1_rows: Vec<&MarkerRow> = rows
        .iter()
        .filter(|m| m.celltype == "CT1" && m.source == "panel")
        .collect();
    assert_eq!(ct1_rows[0].gene, "GENE1");
    let first_nan = ct1_rows.iter().position(|m| m.log2fc.is_nan()).unwrap();
    assert!(
        ct1_rows[first_nan..].iter().all(|m| m.log2fc.is_nan()),
        "measured first: {:?}",
        ct1_rows.iter().map(|m| &m.gene).collect::<Vec<_>>()
    );
}

#[test]
fn a_named_mix_is_placed_by_its_parts() {
    let mut mixed = crate::tui::ontology::Mixed::default();
    let root = tempfile::tempdir().unwrap();
    mixed
        .add(
            "mix1",
            &["CT2".into(), "CT1".into()],
            &root.path().join("m.tsv"),
        )
        .unwrap();
    let calls = cluster_calls(&round(), &tree(), None, &mixed, &[relabel(0, "mix1")]);
    assert_eq!(calls[&0].label, "mix1");
    assert!((calls[&0].share - 0.9).abs() < 1e-6, "CT1 0.7 + CT2 0.2");
}
