//! [`super`]'s calls and marker table on a made-up round.

use super::*;
use crate::annotate::celltype_tree::{TreeSource, TypeGroup, TypeTree};
use crate::tui::round::{Candidate, ClusterView};
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::traits::MatWithNames;

/// `lymph` over `T` and `B`.
fn tree() -> PanelTree {
    PanelTree::from_type_tree(&TypeTree {
        source: TreeSource::MarkerSharing,
        release: None,
        label_cl: BTreeMap::new(),
        groups: vec![TypeGroup {
            name: "lymph".into(),
            cl_id: None,
            members: vec!["T".into(), "B".into()],
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
                nes: None,
                p: None,
                q: None,
            })
            .collect(),
        shares,
        genes: Vec::new(),
    }
}

/// K0 is T, K1 is B; CD3E is up in K0, MS4A1 in K1, NEW in K0 only.
fn round() -> RoundView {
    let set = |g: &[&str]| g.iter().map(|g| (*g).to_string()).collect::<BTreeSet<_>>();
    RoundView {
        manifest: PathBuf::new(),
        clusters: vec![
            cluster(0, "T", &[("T", 0.7), ("B", 0.2)]),
            cluster(1, "B", &[("B", 0.6), ("T", 0.1)]),
        ],
        cell_names: vec!["c0".into(), "c1".into()],
        cell_clusters: vec![Some(0), Some(1)],
        expression: Some(crate::tui::round::Expression::new(MatWithNames {
            rows: vec!["CD3E".into(), "MS4A1".into(), "NEW".into(), "ACTB".into()],
            cols: vec!["K0".into(), "K1".into()],
            mat: Mat::from_row_slice(4, 2, &[8.0, 1.0, 1.0, 8.0, 12.0, 0.0, 20.0, 20.0]),
        })),
        markers: BTreeMap::from([
            ("T".into(), set(&["CD3E", "CD2"])),
            ("B".into(), set(&["MS4A1"])),
        ]),
        loose_cells: 0,
        decided: BTreeSet::new(),
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
        &[relabel(1, "lymph")],
    );
    assert_eq!(calls[&0].label, "T");
    assert_eq!(calls[&0].lineage, "lymph > T");
    assert!((calls[&0].share - 0.7).abs() < 1e-6);
    assert_eq!(calls[&1].label, "lymph");
    assert_eq!(calls[&1].lineage, "lymph");
    assert!((calls[&1].share - 0.7).abs() < 1e-6, "B 0.6 + T 0.1");
    assert_eq!(calls[&1].top, "B");
}

#[test]
fn markers_are_panel_or_added_and_the_rest_suggested() {
    let r = round();
    let t = tree();
    let calls = cluster_calls(&r, &t, None, &crate::tui::ontology::Mixed::default(), &[]);
    // CD2 was added this session: not in the original panel.
    let original = BTreeMap::from([
        ("T".to_string(), BTreeSet::from(["CD3E".to_string()])),
        ("B".to_string(), BTreeSet::from(["MS4A1".to_string()])),
    ]);
    let rows = marker_rows(&r, &t, &calls, &original);
    let of = |ct: &str| -> Vec<(String, &str)> {
        rows.iter()
            .filter(|m| m.celltype == ct)
            .map(|m| (m.gene.clone(), m.source))
            .collect()
    };
    assert_eq!(
        of("T"),
        [
            ("CD3E".to_string(), "panel"),
            ("CD2".to_string(), "added"),
            ("NEW".to_string(), "suggested")
        ]
    );
    let cd3e = rows.iter().find(|m| m.gene == "CD3E").unwrap();
    assert!(cd3e.log2fc > 0.5, "up in the T cluster");
    assert!(
        rows.iter()
            .find(|m| m.gene == "CD2")
            .unwrap()
            .log2fc
            .is_nan(),
        "not measured"
    );
    assert!(
        !rows.iter().any(|m| m.gene == "ACTB"),
        "high but everywhere: not suggested"
    );
    assert_eq!(of("B"), [("MS4A1".to_string(), "panel")]);
}

#[test]
fn a_term_off_the_panel_gets_its_id_and_lineage_from_the_ontology() {
    let obo = "[Term]\nid: CL:0\nname: root cell\n\n\
               [Term]\nid: CL:1\nname: lymph\nsubset: cellxgene_subset\nis_a: CL:0\n\n\
               [Term]\nid: CL:10\nname: T\nis_a: CL:1\n\n\
               [Term]\nid: CL:12\nname: B\nis_a: CL:1\n\n\
               [Term]\nid: CL:3\nname: stem cell\nis_a: CL:0\n";
    let cl = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let panel: Vec<(String, String)> = [("G1", "T"), ("G2", "B")]
        .iter()
        .map(|(g, t)| ((*g).to_string(), (*t).to_string()))
        .collect();
    let t = PanelTree::from_ontology(&cl, &panel).unwrap();
    let calls = cluster_calls(
        &round(),
        &t,
        Some(&cl),
        &crate::tui::ontology::Mixed::default(),
        &[relabel(0, "stem_cell"), relabel(1, "lymph")],
    );
    assert_eq!(calls[&0].cl_id, "CL:3");
    assert_eq!(calls[&0].lineage, "root_cell > stem_cell");
    assert_eq!(calls[&1].cl_id, "CL:1", "a class the panel's tree has");
    assert!(
        (calls[&1].share - 0.7).abs() < 1e-6,
        "B 0.6 + T 0.1, both under lymph"
    );
}

#[test]
fn a_mixed_label_sums_its_parts_and_sits_under_their_common_ancestor() {
    let obo = "[Term]\nid: CL:0\nname: root cell\n\n\
               [Term]\nid: CL:1\nname: lymph\nsubset: cellxgene_subset\nis_a: CL:0\n\n\
               [Term]\nid: CL:10\nname: T\nis_a: CL:1\n\n\
               [Term]\nid: CL:12\nname: B\nis_a: CL:1\n\n\
               [Term]\nid: CL:3\nname: stem cell\nis_a: CL:0\n";
    let cl = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let panel: Vec<(String, String)> = [("G1", "T"), ("G2", "B")]
        .iter()
        .map(|(g, t)| ((*g).to_string(), (*t).to_string()))
        .collect();
    let t = PanelTree::from_ontology(&cl, &panel).unwrap();
    let calls = cluster_calls(
        &round(),
        &t,
        Some(&cl),
        &crate::tui::ontology::Mixed::default(),
        &[relabel(0, "B+T")],
    );
    assert_eq!(calls[&0].label, "B+T");
    assert!((calls[&0].share - 0.9).abs() < 1e-6, "T 0.7 + B 0.2");
    assert_eq!(calls[&0].cl_id, "CL:12|CL:10");
    assert_eq!(calls[&0].lineage, "root_cell > lymph > B+T");
    // Without the ontology: the parts on the panel's tree still add up.
    let alone = cluster_calls(
        &round(),
        &tree(),
        None,
        &crate::tui::ontology::Mixed::default(),
        &[relabel(0, "B+T")],
    );
    assert!((alone[&0].share - 0.9).abs() < 1e-6);
}

#[test]
fn a_named_mix_is_placed_by_its_parts() {
    let mut mixed = crate::tui::ontology::Mixed::default();
    let root = tempfile::tempdir().unwrap();
    mixed
        .add(
            "lymphs",
            &["B".into(), "T".into()],
            &root.path().join("m.tsv"),
        )
        .unwrap();
    let calls = cluster_calls(&round(), &tree(), None, &mixed, &[relabel(0, "lymphs")]);
    assert_eq!(calls[&0].label, "lymphs");
    assert!((calls[&0].share - 0.9).abs() < 1e-6, "T 0.7 + B 0.2");
}
