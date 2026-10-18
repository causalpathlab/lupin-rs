//! Coarse groups for [`super`], on a made-up ontology and panel.

use super::*;

/// root ─┬─ group a ─┬─ CT1          (group a, group b: analysis classes)
///       │           └─ CT2
///       └─ group b ─── CT3 (exact synonym "CT3 alias")
const OBO: &str = "format-version: 1.2
data-version: test/2026-01-01

[Term]
id: CL:9000000
name: root cell

[Term]
id: CL:9000001
name: group a
subset: cellxgene_subset
is_a: CL:9000000 ! root cell

[Term]
id: CL:9000002
name: group b
subset: blood_and_immune_upper_slim
is_a: CL:9000000 ! root cell

[Term]
id: CL:9000011
name: CT1
synonym: \"C1X\" RELATED OMO:0003000 []
is_a: CL:9000001 ! group a

[Term]
id: CL:9000012
name: CT2
is_a: CL:9000001 ! group a

[Term]
id: CL:9000021
name: CT3
synonym: \"CT3 alias\" EXACT []
synonym: \"CT3 loose\" RELATED []
is_a: CL:9000002 ! group b

[Term]
id: CL:9000099
name: CT9
is_obsolete: true

[Typedef]
id: part_of
name: part of
";

fn panel(pairs: &[(&str, &str)]) -> Vec<(String, String)> {
    pairs
        .iter()
        .map(|(g, t)| ((*g).to_string(), (*t).to_string()))
        .collect()
}

fn group(name: &str, members: &[&str]) -> TypeGroup {
    TypeGroup {
        name: name.to_string(),
        cl_id: Some(format!("CL:{name}")),
        members: members.iter().map(|m| (*m).to_string()).collect(),
    }
}

fn members(tree: &TypeTree, group: &str) -> Vec<String> {
    tree.groups
        .iter()
        .find(|g| g.name == group)
        .map(|g| g.members.clone())
        .unwrap_or_default()
}

#[test]
fn obo_terms_keep_names_exact_synonyms_parents_and_release() {
    let t = ClTerms::parse(OBO, &crate::annotate::cl_rules::shipped());
    assert_eq!(t.release.as_deref(), Some("test/2026-01-01"));
    assert_eq!(t.term_count(), 6, "the obsolete term is skipped");
    let (mapped, unmapped) = t.map_labels(["CT1", "ct3_alias", "CT3 loose", "CT9", "CT7", "CT1s"]);
    assert_eq!(mapped["CT1"], "CL:9000011");
    assert_eq!(
        mapped["ct3_alias"], "CL:9000021",
        "case and `_` do not matter"
    );
    assert_eq!(
        unmapped,
        ["CT3 loose", "CT9", "CT7"],
        "only exact synonyms, no obsolete terms"
    );
}

#[test]
fn labels_match_terms_whose_words_come_in_another_order() {
    let obo = "[Term]\nid: CL:1\nname: memory CT1 cell\n\n\
               [Term]\nid: CL:2\nname: xi-omicron CT2 cell\n\n\
               [Term]\nid: CL:3\nname: cell CT1 memory\nis_obsolete: true\n\n\
               [Term]\nid: CL:4\nname: red blue cell\n\n\
               [Term]\nid: CL:5\nname: blue red cell\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let (mapped, unmapped) = t.map_labels([
        "CT1 cells memory",
        "Xi omicron CT2 cells",
        "Red blue cells",
        "Cells red blue",
    ]);
    assert_eq!(mapped["CT1 cells memory"], "CL:1");
    assert_eq!(mapped["Xi omicron CT2 cells"], "CL:2", "hyphens are spaces");
    assert_eq!(mapped["Red blue cells"], "CL:4", "an exact name wins");
    assert_eq!(unmapped, ["Cells red blue"], "two terms share its words");
}

#[test]
fn a_label_that_leaves_off_cell_matches_the_term_that_ends_in_it() {
    let obo = "[Term]\nid: CL:1\nname: immature CT1 cell\n\n\
               [Term]\nid: CL:2\nname: pro-CT2 cell\n\n\
               [Term]\nid: CL:3\nname: CT3\n\n\
               [Term]\nid: CL:4\nname: CT3 cell\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let (mapped, unmapped) = t.map_labels(["Immature_CT1", "Pro-CT2", "CT3", "CT9"]);
    assert_eq!(mapped["Immature_CT1"], "CL:1");
    assert_eq!(mapped["Pro-CT2"], "CL:2", "a hyphen parts two words");
    assert_eq!(mapped["CT3"], "CL:3", "a name as written wins");
    assert_eq!(unmapped, ["CT9"]);
    let names_only = "[Term]\nid: CL:5\nname: alpha retinal ganglion cell\nsynonym: \"alpha type cell\" EXACT []\n\n\
                      [Term]\nid: CL:6\nname: endothelial cell\n";
    let t1 = ClTerms::parse(names_only, &crate::annotate::cl_rules::shipped());
    let (m1, u1) = t1.map_labels(["Alpha_Type", "Endothelial"]);
    assert_eq!(m1["Endothelial"], "CL:6", "one word completes a name");
    assert_eq!(u1, ["Alpha_Type"], "but not a synonym");
    let literal = ClTerms::parse(obo, &crate::annotate::cl_rules::MatchRules::default());
    assert_eq!(
        literal.map_labels(["Immature_CT1"]).1,
        ["Immature_CT1"],
        "only when the rules say so"
    );
}

#[test]
fn an_unmatched_type_joins_the_group_holding_its_specific_markers() {
    // X shares three markers most types hold with B and C, and two specific
    // markers with A: A explains it, though raw overlap favours B and C.
    let pairs = panel(&[
        ("S1", "A1"),
        ("S2", "A1"),
        ("a1", "A1"),
        ("P1", "B1"),
        ("P2", "B1"),
        ("P3", "B1"),
        ("b1", "B1"),
        ("P1", "B2"),
        ("P2", "B2"),
        ("P3", "B2"),
        ("b2", "B2"),
        ("P1", "C1"),
        ("P2", "C1"),
        ("P3", "C1"),
        ("c1", "C1"),
        ("P1", "X"),
        ("P2", "X"),
        ("P3", "X"),
        ("S1", "X"),
        ("S2", "X"),
    ]);
    let mut groups = vec![
        group("a", &["A1"]),
        group("b", &["B1", "B2"]),
        group("c", &["C1"]),
    ];
    attach_by_sharing(&mut groups, vec!["X".to_string()], &pairs, |_, _| true);
    assert_eq!(groups[0].members, ["A1", "X"]);
}

#[test]
fn a_broad_group_does_not_outbid_a_near_copy() {
    // X's two markers are A1's; B1 and B2 hold one each between them.
    let pairs = panel(&[
        ("x1", "A1"),
        ("x2", "A1"),
        ("a", "A1"),
        ("x1", "B1"),
        ("b1", "B1"),
        ("x2", "B2"),
        ("b2", "B2"),
        ("x1", "X"),
        ("x2", "X"),
    ]);
    let mut groups = vec![group("a", &["A1"]), group("b", &["B1", "B2"])];
    attach_by_sharing(&mut groups, vec!["X".to_string()], &pairs, |_, _| true);
    assert_eq!(groups[0].members, ["A1", "X"]);
}

#[test]
fn a_marker_almost_every_type_holds_does_not_stop_the_split() {
    let tree = TypeTree::from_markers(&panel(&[
        ("U", "A"),
        ("P", "A"),
        ("x", "A"),
        ("U", "B"),
        ("P", "B"),
        ("x", "B"),
        ("U", "C"),
        ("P", "C"),
        ("y", "C"),
        ("U", "D"),
        ("P", "D"),
        ("y", "D"),
        ("U", "E"),
        ("e", "E"),
    ]));
    let groups: Vec<Vec<String>> = tree.groups.iter().map(|g| g.members.clone()).collect();
    assert!(
        groups.contains(&vec!["A".to_string(), "B".to_string()]),
        "{groups:?}"
    );
    assert!(
        groups.contains(&vec!["C".to_string(), "D".to_string()]),
        "{groups:?}"
    );
}

#[test]
fn types_linked_only_by_a_marker_all_hold_stay_apart() {
    let tree = TypeTree::from_markers(&panel(&[
        ("U", "A"),
        ("a", "A"),
        ("U", "B"),
        ("b", "B"),
        ("U", "C"),
        ("c", "C"),
    ]));
    assert_eq!(tree.groups.len(), 3, "{:?}", tree.groups);
}

#[test]
fn every_component_is_split_not_only_a_whole_panel() {
    // E holds only its own marker, so the panel is not one component; the
    // other four are still split by what they specifically share.
    let tree = TypeTree::from_markers(&panel(&[
        ("P", "A"),
        ("x", "A"),
        ("P", "B"),
        ("x", "B"),
        ("P", "C"),
        ("y", "C"),
        ("P", "D"),
        ("y", "D"),
        ("e", "E"),
    ]));
    let groups: Vec<Vec<String>> = tree.groups.iter().map(|g| g.members.clone()).collect();
    assert!(
        groups.contains(&vec!["A".to_string(), "B".to_string()]),
        "{groups:?}"
    );
    assert!(
        groups.contains(&vec!["C".to_string(), "D".to_string()]),
        "{groups:?}"
    );
}

#[test]
fn a_dilute_broad_group_still_takes_a_near_copy_of_its_member() {
    // X is nearly T1, whose group is otherwise unlike X; N1 shares one marker.
    let mut rows = vec![
        ("x1", "T1"),
        ("x2", "T1"),
        ("x3", "T1"),
        ("t", "T1"),
        ("x1", "X"),
        ("x2", "X"),
        ("x3", "X"),
        ("t", "X"),
        ("n1", "X"),
        ("n1", "N1"),
        ("m", "N1"),
    ];
    let others = ["T2", "T3", "T4", "T5", "T6"];
    for (o, own) in others.iter().zip(["o2", "o3", "o4", "o5", "o6"]) {
        rows.push(("t", o));
        rows.push((own, o));
    }
    let pairs = panel(&rows);
    let mut members = vec!["T1"];
    members.extend(others);
    let mut groups = vec![group("t", &members), group("n", &["N1"])];
    attach_by_sharing(&mut groups, vec!["X".to_string()], &pairs, |_, _| true);
    assert!(groups[0].members.iter().any(|m| m == "X"), "{groups:?}");
}

#[test]
fn the_original_overlap_case_still_joins_its_specific_group() {
    let pairs = panel(&[
        ("S1", "A1"),
        ("a1", "A1"),
        ("a2", "A1"),
        ("P1", "B1"),
        ("P2", "B1"),
        ("P3", "B1"),
        ("b1", "B1"),
        ("P1", "B2"),
        ("P2", "B2"),
        ("P3", "B2"),
        ("b2", "B2"),
        ("P1", "C1"),
        ("P2", "C1"),
        ("P3", "C1"),
        ("c1", "C1"),
        ("P1", "X"),
        ("P2", "X"),
        ("P3", "X"),
        ("S1", "X"),
    ]);
    let mut groups = vec![
        group("a", &["A1"]),
        group("b", &["B1", "B2"]),
        group("c", &["C1"]),
    ];
    attach_by_sharing(&mut groups, vec!["X".to_string()], &pairs, |_, _| true);
    assert_eq!(groups[0].members, ["A1", "X"]);
}
#[test]
fn marker_groups_are_not_made_by_markers_most_types_hold() {
    // P1..P4 are held by six of seven types (HLA-like); B1 shares only `b`
    // with B2, but four P markers with each D and E type.
    let mut rows = vec![
        ("b", "B1"),
        ("b", "B2"),
        ("x1", "B2"),
        ("d", "D1"),
        ("d", "D2"),
        ("d1", "D1"),
        ("d2", "D2"),
    ];
    for t in ["B1", "D1", "D2", "E1", "E2", "E3"] {
        for p in ["P1", "P2", "P3", "P4"] {
            rows.push((p, t));
        }
    }
    rows.extend([("e1", "E1"), ("e2", "E2"), ("e3", "E3")]);
    let tree = TypeTree::from_markers(&panel(&rows));
    let group_of = |t: &str| {
        tree.groups
            .iter()
            .position(|g| g.members.iter().any(|m| m == t))
            .unwrap()
    };
    assert_eq!(group_of("B1"), group_of("B2"), "{:?}", tree.groups);
    assert_ne!(group_of("B1"), group_of("D1"), "{:?}", tree.groups);
}

#[test]
fn an_implied_word_needs_a_word_to_complete() {
    let obo = "[Term]\nid: CL:0\nname: cell\n\n\
               [Term]\nid: CL:1\nname: progenitor cell\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let (mapped, unmapped) = t.map_labels(["", "-", "Progenitor"]);
    assert_eq!(
        mapped["Progenitor"], "CL:1",
        "a general label, a general term"
    );
    assert_eq!(unmapped, ["", "-"], "no word to complete");
}
#[test]
fn a_terms_own_name_beats_another_terms_synonym_in_every_tier() {
    let obo = "[Term]\nid: CL:1\nname: pro CT8 cell\n\n\
               [Term]\nid: CL:2\nname: early CT8 cell\nsynonym: \"pro CT8 cell\" EXACT []\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let labels = ["Pro_CT8", "Pro_CT8_cell", "Pros CT8 cell", "CT8 cell pro"];
    let (mapped, _) = t.map_labels(labels);
    for l in labels {
        assert_eq!(mapped.get(l).map(String::as_str), Some("CL:1"), "{l}");
    }
}

#[test]
fn an_ambiguous_implied_word_gives_way_to_the_next() {
    use crate::annotate::cl_rules::MatchRules;
    let rules = MatchRules::from_layers(&[serde_json::json!({
        "synonyms": {"EXACT": []},
        "implied": {"suffixes": ["cell", "lymphocyte"]}
    })])
    .unwrap();
    let obo = "[Term]\nid: CL:1\nname: memory CT9 cell\n\n\
               [Term]\nid: CL:2\nname: memory CT9 cell\n\n\
               [Term]\nid: CL:3\nname: memory CT9 lymphocyte\n";
    let t = ClTerms::parse(obo, &rules);
    assert_eq!(t.map_labels(["Memory_CT9"]).0["Memory_CT9"], "CL:3");
}

#[test]
fn ontology_groups_are_the_nearest_shared_classes() {
    let t = ClTerms::parse(OBO, &crate::annotate::cl_rules::shipped());
    let p = panel(&[
        ("GENE1", "CT1"),
        ("GENE2", "CT2"),
        ("GENE3", "CT3"),
        ("GENE3", "CT_4"),
        ("GENE9", "CT5"),
    ]);
    let tree = TypeTree::from_ontology(&t, &p).unwrap();
    assert_eq!(tree.source, TreeSource::CellOntology);
    assert_eq!(members(&tree, "group_a"), ["CT1", "CT2"]);
    // CT3 shares group b with no other type, so it stands as itself; CT_4
    // matches no term but shares a marker with CT3.
    assert_eq!(members(&tree, "CT3"), ["CT3", "CT_4"]);
    assert_eq!(
        members(&tree, "CT5"),
        ["CT5"],
        "no term, no shared marker: its own group"
    );
    assert_eq!(tree.group_of("CT 4"), Some("CT3"));
    assert_eq!(tree.label_cl["CT3"], "CL:9000021");
}

#[test]
fn ontology_groups_need_an_analysis_class_shared_by_two_types() {
    let t = ClTerms::parse(OBO, &crate::annotate::cl_rules::shipped());
    // CT1 and CT3 share only the root, which is no analysis class.
    let tree = TypeTree::from_ontology(&t, &panel(&[("GENE1", "CT1"), ("GENE3", "CT3")])).unwrap();
    assert_eq!(members(&tree, "CT1"), ["CT1"]);
    assert_eq!(members(&tree, "CT3"), ["CT3"]);
    // Fewer than two matched types: no ontology tree.
    assert!(TypeTree::from_ontology(&t, &panel(&[("GENE1", "CT1"), ("GENE7", "CT7")])).is_none());
}

#[test]
fn marker_sharing_groups_types_with_common_markers() {
    let tree = TypeTree::from_markers(&panel(&[
        ("GENE1", "CT1"),
        ("GENE2", "CT1"),
        ("GENE2", "CT2"),
        ("GENE3", "CT3"),
    ]));
    assert_eq!(tree.source, TreeSource::MarkerSharing);
    assert_eq!(members(&tree, "CT1/CT2"), ["CT1", "CT2"]);
    assert_eq!(members(&tree, "CT3"), ["CT3"]);
}

#[test]
fn one_connected_panel_is_split_at_its_widest_gap() {
    // CT1–CT2 and CT3–CT4 overlap heavily; CT2–CT3 barely.
    let tree = TypeTree::from_markers(&panel(&[
        ("GENE1", "CT1"),
        ("GENE2", "CT1"),
        ("GENE1", "CT2"),
        ("GENE2", "CT2"),
        ("GENE9", "CT2"),
        ("GENE9", "CT3"),
        ("GENE3", "CT3"),
        ("GENE4", "CT3"),
        ("GENE3", "CT4"),
        ("GENE4", "CT4"),
    ]));
    let mut groups: Vec<Vec<String>> = tree
        .groups
        .iter()
        .map(|g| {
            let mut m = g.members.clone();
            m.sort();
            m
        })
        .collect();
    groups.sort();
    assert_eq!(groups, [vec!["CT1", "CT2"], vec!["CT3", "CT4"]]);
}

#[test]
fn a_coarse_call_sums_evidence_over_a_group() {
    let tree = TypeTree::from_markers(&panel(&[
        ("GENE1", "CT1"),
        ("GENE1", "CT2"),
        ("GENE3", "CT3"),
    ]));
    let index = tree.index();
    // CT3 is the single best type, but CT1 and CT2 together outweigh it.
    let types = ["CT1".to_string(), "CT2".to_string(), "CT3".to_string()];
    let probs = [0.3, 0.3, 0.4];
    assert_eq!(
        coarse_call(&index, Some((&types, &probs)), &[]).as_deref(),
        Some("CT1/CT2")
    );
    // No evidence at all falls back to the cells' labels.
    assert_eq!(
        coarse_call(&index, Some((&types, &[0.0; 3])), &["CT3"]).as_deref(),
        Some("CT3")
    );
    // Without probabilities, the cells' fine labels vote.
    assert_eq!(
        coarse_call(&index, None, &["CT3", "CT3", "CT1"]).as_deref(),
        Some("CT3")
    );
    assert_eq!(coarse_call(&index, None, &[]), None);
}

#[test]
fn abbreviations_match_but_other_related_synonyms_do_not() {
    let t = ClTerms::parse(OBO, &crate::annotate::cl_rules::shipped());
    let (mapped, unmapped) = t.map_labels(["C1X", "C1Xs", "CT3 loose"]);
    assert_eq!(
        mapped["C1X"], "CL:9000011",
        "an abbreviation (OMO:0003000) matches"
    );
    assert_eq!(mapped["C1Xs"], "CL:9000011", "in the plural too");
    assert_eq!(
        unmapped,
        ["CT3 loose"],
        "a plain RELATED synonym still does not"
    );
}

#[test]
fn the_ontology_can_be_walked_and_searched() {
    let t = ClTerms::parse(OBO, &crate::annotate::cl_rules::shipped());
    assert!(
        t.has("CL:9000001") && !t.has("CL:9000099"),
        "obsolete terms are gone"
    );
    assert_eq!(t.parents("CL:9000011"), ["CL:9000001"]);
    assert_eq!(t.children("CL:9000001"), ["CL:9000011", "CL:9000012"]);
    assert!(t.children("CL:9000011").is_empty());
    assert_eq!(
        t.lineage("CL:9000011"),
        ["CL:9000000", "CL:9000001", "CL:9000011"]
    );
    // Whole-name hits first, then shorter names.
    assert_eq!(t.search("group", 10), ["CL:9000001", "CL:9000002"]);
    assert_eq!(
        t.search("c1x", 10),
        ["CL:9000011"],
        "abbreviations are searchable"
    );
    assert_eq!(t.search("ct", 2).len(), 2, "capped at the limit");
    assert!(t.search("  ", 10).is_empty());
}

#[test]
fn a_search_ranks_terms_named_by_the_query_above_synonym_hits() {
    let obo = "format-version: 1.2

[Term]
id: CL:1
name: T cell
synonym: \"T lymphocyte\" EXACT []

[Term]
id: CL:2
name: lymphocyte

[Term]
id: CL:3
name: innate lymphoid cell
";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    assert_eq!(
        t.search("lympho", 10),
        ["CL:2", "CL:3", "CL:1"],
        "a name starting with it, then a name with it, then a synonym"
    );
}

#[test]
fn a_label_is_placed_on_its_own_term_else_its_groups() {
    let tree = TypeTree {
        source: TreeSource::CellOntology,
        release: None,
        label_cl: [("CT1".to_string(), "CL:9000001".to_string())].into(),
        label_under: BTreeMap::new(),
        groups: vec![TypeGroup {
            name: "group_a".into(),
            cl_id: Some("CL:9000009".into()),
            members: vec!["CT1".into(), "CT2".into()],
        }],
    };
    let at = |l: &str| tree.placement(l).map(|p| (p.id, p.name, p.via));
    assert_eq!(
        at("CT1"),
        Some(("CL:9000001".into(), "CT1".into(), Via::Label))
    );
    assert_eq!(
        at("CT2"),
        Some(("CL:9000009".into(), "group_a".into(), Via::Group))
    );
    assert_eq!(
        at("group_a"),
        Some(("CL:9000009".into(), "group_a".into(), Via::Group))
    );
    assert_eq!(at("CT9"), None);
}

#[test]
fn a_text_naming_a_term_and_its_specialisations_takes_the_most_general() {
    let obo = "[Term]\nid: CL:1\nname: CT1 progenitor cell\nsynonym: \"CTP\" RELATED OMO:0003000 []\n\n\
               [Term]\nid: CL:2\nname: CT1 progenitor cell, mammalian\nsynonym: \"CTP\" RELATED OMO:0003000 []\nis_a: CL:1\n\n\
               [Term]\nid: CL:3\nname: CT2 cell\nsynonym: \"C2\" EXACT []\n\n\
               [Term]\nid: CL:4\nname: CT4 cell\nsynonym: \"C2\" EXACT []\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let (mapped, unmapped) = t.map_labels(["CTP", "C2"]);
    assert_eq!(
        mapped["CTP"], "CL:1",
        "the general term, never a guess below it"
    );
    assert_eq!(unmapped, ["C2"], "two unrelated terms stay a tie");
}

#[test]
fn a_name_is_also_known_without_its_trailing_qualifier() {
    let obo = "[Term]\nid: CL:1\nname: CT5 neuron (sensu Vertebrata)\n\n\
               [Term]\nid: CL:2\nname: CT6 cell, human\n\n\
               [Term]\nid: CL:3\nname: CT7 cell (Mmus)\n\n\
               [Term]\nid: CL:4\nname: CT7 cell (Primate)\n\n\
               [Term]\nid: CL:5\nname: CT8 cell, CD9-positive\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let (mapped, unmapped) = t.map_labels(["CT5 neurons", "CT6", "CT7", "CT8"]);
    assert_eq!(mapped["CT5 neurons"], "CL:1");
    assert_eq!(mapped["CT6"], "CL:2", "an implied word completes it too");
    assert_eq!(
        unmapped,
        ["CT7", "CT8"],
        "two species' terms are a tie; a marker is no qualifier"
    );
    let literal = ClTerms::parse(obo, &crate::annotate::cl_rules::MatchRules::default());
    assert_eq!(literal.map_labels(["CT5 neuron"]).1, ["CT5 neuron"]);
}

#[test]
fn plural_forms_match_as_the_rules_count_them() {
    let obo = "[Term]\nid: CL:1\nname: CT8ial cell\nsynonym: \"CT8ia\" RELATED OMO:0003004 []\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    assert_eq!(t.map_labels(["CT8ia"]).0["CT8ia"], "CL:1");
}

#[test]
fn separators_may_be_ignored_as_a_last_resort() {
    let obo = "[Term]\nid: CL:1\nname: proxyblast\n\n\
               [Term]\nid: CL:2\nname: prezed cell\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let (mapped, unmapped) = t.map_labels(["Pro-Xyblasts", "Pre-Zed"]);
    assert_eq!(mapped["Pro-Xyblasts"], "CL:1");
    assert_eq!(mapped["Pre-Zed"], "CL:2", "with an implied word");
    assert!(unmapped.is_empty());
    let spaced = format!("{obo}\n[Term]\nid: CL:3\nname: pre zed cell\n");
    let t = ClTerms::parse(&spaced, &crate::annotate::cl_rules::shipped());
    assert_eq!(
        t.map_labels(["Pre-Zed"]).0["Pre-Zed"],
        "CL:3",
        "the words as written first"
    );
    let literal = ClTerms::parse(obo, &crate::annotate::cl_rules::MatchRules::default());
    assert_eq!(literal.map_labels(["Pro-Xyblast"]).1, ["Pro-Xyblast"]);
}

/// A lineage written out in a term's name, with its specialisations under it.
const LINEAGE: &str = "[Term]\nid: CL:1\nname: CT4-positive, alpha-beta T cell\n\n\
    [Term]\nid: CL:2\nname: central memory CT4-positive, alpha-beta T cell\nis_a: CL:1\n\n\
    [Term]\nid: CL:3\nname: central memory CT4-positive, alpha-beta T cell, human\nis_a: CL:2\n\n\
    [Term]\nid: CL:4\nname: central memory CT4-negative, alpha-beta T cell\n\n\
    [Term]\nid: CL:5\nname: effector memory CT4-positive, alpha-beta T cell\nis_a: CL:1\n\n\
    [Term]\nid: CL:6\nname: central nervous system cell\n\n\
    [Term]\nid: CL:7\nname: effector memory CT4-positive, alpha-beta T cell, terminally differentiated\nis_a: CL:1\n";

#[test]
fn a_label_whose_words_a_term_spells_out_takes_the_most_general_such_term() {
    let t = ClTerms::parse(LINEAGE, &crate::annotate::cl_rules::shipped());
    let (mapped, unmapped) = t.map_labels([
        "CT4 Central Memory",
        "CT4+ central memory",
        "CT4 Effector Memory",
        "CT4 Memory",
        "Central",
    ]);
    assert_eq!(mapped["CT4 Central Memory"], "CL:2", "not CT4-negative's");
    assert_eq!(mapped["CT4+ central memory"], "CL:2");
    assert_eq!(
        mapped["CT4 Effector Memory"], "CL:5",
        "with no term above both, the name adding the fewest words"
    );
    assert_eq!(
        unmapped,
        ["CT4 Memory", "Central"],
        "names adding as many words, with nothing named above them; one word is too little"
    );
    let literal = ClTerms::parse(LINEAGE, &crate::annotate::cl_rules::MatchRules::default());
    assert_eq!(literal.map_labels(["CT4 Central Memory"]).0.len(), 0);
}

#[test]
fn a_word_the_panel_leaves_off_its_labels_completes_them() {
    let obo = "[Term]\nid: CL:1\nname: organx A cell\nsynonym: \"organx alpha cell\" EXACT []\n\n\
               [Term]\nid: CL:2\nname: organx beta cell\n\n\
               [Term]\nid: CL:3\nname: organx gamma cell\n\n\
               [Term]\nid: CL:4\nname: eyey ganglion cell\nsynonym: \"alpha cell\" EXACT []\n\n\
               [Term]\nid: CL:5\nname: brainz delta cell\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let (mapped, unmapped) = t.map_labels(["Alpha", "Beta", "Gamma", "Delta", "Duct"]);
    assert_eq!(
        mapped["Alpha"], "CL:1",
        "the panel's word, not another tissue's synonym"
    );
    assert_eq!(mapped["Beta"], "CL:2");
    assert_eq!(mapped["Gamma"], "CL:3");
    assert_eq!(
        unmapped,
        ["Delta", "Duct"],
        "another word, for one label only"
    );
    let (few, _) = t.map_labels(["Alpha", "Beta", "Delta"]);
    assert!(
        few.is_empty(),
        "two labels are too few to say what the panel leaves off"
    );
}

#[test]
fn a_known_context_word_completes_by_lookup_as_the_scan_does() {
    let obo = "[Term]\nid: CL:1\nname: organx A cell\nsynonym: \"organx alpha cell\" EXACT []\n\n\
               [Term]\nid: CL:2\nname: organx beta cell\n\n\
               [Term]\nid: CL:3\nname: organx gamma cell\n\n\
               [Term]\nid: CL:4\nname: organx duct\n\n\
               [Term]\nid: CL:5\nname: brainz delta cell\n\n\
               [Term]\nid: CL:6\nname: organx cell\n";
    let mut t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let panel: Vec<(String, String)> = ["Alpha", "Beta", "Gamma", "Duct", "Delta", "Pi"]
        .map(|l| (String::new(), l.to_string()))
        .into();
    let (_, unmapped) = t.map_labels(panel.iter().map(|(_, l)| l.as_str()));
    let labels: Vec<String> = [
        "Alpha",
        "Beta",
        "Duct",
        "Delta",
        "Pi",
        "organx beta cell",
        "-",
    ]
    .map(String::from)
    .into();
    let scanned = t.completions(&labels).remove("organx").unwrap_or_default();
    assert_eq!(t.completed_by(&labels, "organx"), scanned);
    assert_eq!(
        scanned.keys().collect::<Vec<_>>(),
        ["Alpha", "Beta", "Duct"],
        "with or without the implied word; not another word's label"
    );
    assert_eq!(unmapped, ["Delta", "Pi"]);
    assert!(
        t.map_labels(["Gamma"]).0.is_empty(),
        "one label is too few to tell"
    );
    t.set_panel(&panel);
    assert_eq!(t.panel_context, Some(Some("organx".into())));
    let (one, _) = t.map_labels(["Gamma"]);
    assert_eq!(one["Gamma"], "CL:3", "one label read with the panel's word");
    assert!(
        t.completed_by(&labels, "cell").is_empty(),
        "an implied word is no context"
    );
}

/// Two lineages (analysis classes) of two types each, and kinds of three
/// types the ontology names with a word before them: `early X` is a kind of
/// `X`; `pre X` is not.
const MODIFIED: &str = "[Term]\nid: CL:0\nname: root cell\n\n\
    [Term]\nid: CL:100\nname: lin A cell\nsubset: cellxgene_subset\nis_a: CL:0\n\n\
    [Term]\nid: CL:200\nname: lin B cell\nsubset: cellxgene_subset\nis_a: CL:0\n\n\
    [Term]\nid: CL:1\nname: CT1 cell\nis_a: CL:100\n\n\
    [Term]\nid: CL:2\nname: CT2 cell\nis_a: CL:100\n\n\
    [Term]\nid: CL:3\nname: CT3 cell\nis_a: CL:0\n\n\
    [Term]\nid: CL:4\nname: CT4 cell\nis_a: CL:0\n\n\
    [Term]\nid: CL:5\nname: CT5 cell\nis_a: CL:200\n\n\
    [Term]\nid: CL:6\nname: CT6 cell\nis_a: CL:200\n\n\
    [Term]\nid: CL:11\nname: early CT1 cell\nis_a: CL:1\n\n\
    [Term]\nid: CL:12\nname: early CT2 cell\nis_a: CL:2\n\n\
    [Term]\nid: CL:13\nname: early CT3 cell\nis_a: CL:3\n\n\
    [Term]\nid: CL:21\nname: pre CT1 cell\nis_a: CL:0\n\n\
    [Term]\nid: CL:22\nname: pre CT2 cell\nis_a: CL:0\n\n\
    [Term]\nid: CL:23\nname: pre CT3 cell\nis_a: CL:0\n";

#[test]
fn a_label_is_placed_under_a_term_by_words_the_ontology_uses_for_its_kinds() {
    let t = ClTerms::parse(MODIFIED, &crate::annotate::cl_rules::shipped());
    let placed = t.place_labels([
        "Early CT5",
        "Pre CT5",
        "Late CT5",
        "Early",
        "CT5 x CT6",
        "CT1",
    ]);
    assert_eq!(placed.on["CT1"], "CL:1");
    assert_eq!(placed.under["Early CT5"], "CL:5", "an early CT5 is a CT5");
    assert_eq!(
        placed.hint["Pre CT5"], "CL:5",
        "a pre CT5 is only of its lineage"
    );
    assert_eq!(
        placed.hint["Late CT5"], "CL:5",
        "a word the ontology never uses"
    );
    assert_eq!(
        placed.none,
        ["Early", "CT5 x CT6"],
        "nothing left, or two terms"
    );
    let literal = ClTerms::parse(MODIFIED, &crate::annotate::cl_rules::MatchRules::default());
    let placed = literal.place_labels(["Early CT5 cell"]);
    assert!(placed.under.is_empty() && placed.hint.is_empty());
}

#[test]
fn a_type_under_a_term_groups_by_it_and_a_hinted_one_keeps_to_its_lineage() {
    let t = ClTerms::parse(MODIFIED, &crate::annotate::cl_rules::shipped());
    // Pre CT5 shares two markers with CT1 and one with CT6: by markers alone
    // it would join lineage A.
    let p = panel(&[
        ("G1", "CT1"),
        ("G2", "CT2"),
        ("G5", "CT5"),
        ("G6", "CT6"),
        ("G1", "Pre CT5"),
        ("G2", "Pre CT5"),
        ("G6", "Pre CT5"),
        ("G7", "Early CT6"),
        ("G8", "Early CT4"),
    ]);
    let tree = TypeTree::from_ontology(&t, &p).unwrap();
    let alone = tree
        .groups
        .iter()
        .find(|g| g.members == ["Early CT4"])
        .unwrap();
    assert_eq!(
        (alone.name.as_str(), alone.cl_id.as_deref()),
        ("Early_CT4", Some("CL:4")),
        "alone, named for itself, not for the term it is a kind of"
    );
    assert_eq!(members(&tree, "lin_A_cell"), ["CT1", "CT2"]);
    assert_eq!(
        members(&tree, "lin_B_cell"),
        ["CT5", "CT6", "Early CT6", "Pre CT5"]
    );
    assert!(
        !tree.label_cl.contains_key("Early CT6"),
        "placed, not matched"
    );
}

#[test]
fn a_group_is_not_one_of_its_members_types() {
    // CT1low is a kind of CT1 in the ontology; as a group, CT1 would make it
    // look like one.
    let obo = "[Term]\nid: CL:0\nname: root cell\n\n\
               [Term]\nid: CL:1\nname: lin M cell\nsubset: cellxgene_subset\nis_a: CL:0\n\n\
               [Term]\nid: CL:2\nname: CT1 cell\nsynonym: \"C1X\" EXACT []\nsubset: cellxgene_subset\nis_a: CL:1\n\n\
               [Term]\nid: CL:3\nname: CT1low cell\nis_a: CL:2\n\n\
               [Term]\nid: CL:4\nname: CT5 cell\nsubset: cellxgene_subset\nis_a: CL:0\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let tree = TypeTree::from_ontology(
        &t,
        &panel(&[("G1", "CT1"), ("G2", "CT1low"), ("G5", "CT5")]),
    )
    .unwrap();
    assert_eq!(members(&tree, "lin_M_cell"), ["CT1", "CT1low"]);
    let same = TypeTree::from_ontology(&t, &panel(&[("G1", "CT1"), ("G2", "C1X"), ("G5", "CT5")]))
        .unwrap();
    assert_eq!(
        members(&same, "CT1_cell"),
        ["C1X", "CT1"],
        "both types are that term"
    );
}

#[test]
fn a_group_inside_another_joins_it() {
    // Naive and Reg are nearest CT4 T; CM and EM share a class below it.
    let obo = "[Term]\nid: CL:0\nname: root cell\n\n\
               [Term]\nid: CL:1\nname: lin T cell\nsubset: cellxgene_subset\nis_a: CL:0\n\n\
               [Term]\nid: CL:2\nname: CT4 T cell\nsubset: cellxgene_subset\nis_a: CL:1\n\n\
               [Term]\nid: CL:3\nname: CT4 memory T cell\nsubset: cellxgene_subset\nis_a: CL:2\n\n\
               [Term]\nid: CL:11\nname: naive CT4 T cell\nis_a: CL:2\n\n\
               [Term]\nid: CL:12\nname: reg CT4 T cell\nis_a: CL:2\n\n\
               [Term]\nid: CL:13\nname: CM CT4 T cell\nis_a: CL:3\n\n\
               [Term]\nid: CL:14\nname: EM CT4 T cell\nis_a: CL:3\n\n\
               [Term]\nid: CL:20\nname: lin B cell\nsubset: cellxgene_subset\nis_a: CL:0\n\n\
               [Term]\nid: CL:21\nname: CT8 cell\nis_a: CL:20\n\n\
               [Term]\nid: CL:22\nname: CT9 cell\nis_a: CL:20\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let p = panel(&[
        ("G1", "Naive CT4 T"),
        ("G2", "Reg CT4 T"),
        ("G3", "CM CT4 T"),
        ("G4", "EM CT4 T"),
        ("G8", "CT8"),
        ("G9", "CT9"),
    ]);
    let tree = TypeTree::from_ontology(&t, &p).unwrap();
    assert_eq!(
        members(&tree, "CT4_T_cell"),
        ["CM CT4 T", "EM CT4 T", "Naive CT4 T", "Reg CT4 T"]
    );
    assert_eq!(
        members(&tree, "lin_B_cell"),
        ["CT8", "CT9"],
        "beside it, not in it"
    );
    assert_eq!(tree.groups.len(), 2);
}

#[test]
fn groups_below_one_curated_lineage_are_one_group() {
    // CT4 and CT8 are both T; lin T and lin B both sit in one general system,
    // which is no lineage.
    let obo = "[Term]\nid: CL:0\nname: root cell\n\n\
               [Term]\nid: CL:1\nname: sys cell\nsubset: general_cell_types_upper_slim\nis_a: CL:0\n\n\
               [Term]\nid: CL:2\nname: lin T cell\nsubset: blood_and_immune_upper_slim\nis_a: CL:1\n\n\
               [Term]\nid: CL:3\nname: CT4 T cell\nsubset: cellxgene_subset\nis_a: CL:2\n\n\
               [Term]\nid: CL:4\nname: CT8 T cell\nsubset: cellxgene_subset\nis_a: CL:2\n\n\
               [Term]\nid: CL:11\nname: naive CT4 T cell\nis_a: CL:3\n\n\
               [Term]\nid: CL:12\nname: CM CT4 T cell\nis_a: CL:3\n\n\
               [Term]\nid: CL:13\nname: naive CT8 T cell\nis_a: CL:4\n\n\
               [Term]\nid: CL:14\nname: CM CT8 T cell\nis_a: CL:4\n\n\
               [Term]\nid: CL:20\nname: lin B cell\nsubset: blood_and_immune_upper_slim\nis_a: CL:1\n\n\
               [Term]\nid: CL:21\nname: CT5 cell\nis_a: CL:20\n\n\
               [Term]\nid: CL:22\nname: CT6 cell\nis_a: CL:20\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let p = panel(&[
        ("G1", "Naive CT4 T"),
        ("G2", "CM CT4 T"),
        ("G3", "Naive CT8 T"),
        ("G4", "CM CT8 T"),
        ("G5", "CT5"),
        ("G6", "CT6"),
    ]);
    let tree = TypeTree::from_ontology(&t, &p).unwrap();
    assert_eq!(
        members(&tree, "lin_T_cell"),
        ["CM CT4 T", "CM CT8 T", "Naive CT4 T", "Naive CT8 T"]
    );
    assert_eq!(members(&tree, "lin_B_cell"), ["CT5", "CT6"]);
    assert_eq!(tree.groups.len(), 2);
}

#[test]
fn a_lone_type_s_group_that_takes_another_is_named_for_its_term() {
    // Early CT4 stands alone under CT4, named for itself; Stray matches no
    // term and shares its markers. Together they are CT4's group, not one of
    // them: a group named for a member reads as that member's type.
    let t = ClTerms::parse(MODIFIED, &crate::annotate::cl_rules::shipped());
    let p = panel(&[
        ("G1", "CT1"),
        ("G2", "CT2"),
        ("G8", "Early CT4"),
        ("G9", "Early CT4"),
        ("G8", "Stray"),
        ("G9", "Stray"),
    ]);
    let tree = TypeTree::from_ontology(&t, &p).unwrap();
    let both = (tree.groups.iter())
        .find(|g| g.members.iter().any(|m| m == "Stray"))
        .unwrap();
    assert_eq!(both.members, ["Early CT4", "Stray"]);
    assert_eq!(
        (both.name.as_str(), both.cl_id.as_deref()),
        ("CT4_cell", Some("CL:4"))
    );
}

/// Terms named for markers, with and without their sign, and hyphenated
/// names whose hyphens part words of letters alone.
const SIGNED: &str = "[Term]\nid: CL:1\nname: CD4 T cell\n\n\
    [Term]\nid: CL:2\nname: CD34 cell\n\n\
    [Term]\nid: CL:3\nname: CD4 CD8 thymocyte\n\n\
    [Term]\nid: CL:4\nname: CD16-negative monocyte\nsynonym: \"CD16- monocyte\" EXACT []\n\n\
    [Term]\nid: CL:5\nname: pro-B cell\n\n\
    [Term]\nid: CL:6\nname: T cell\n";

#[test]
fn a_marker_s_sign_is_kept_however_the_label_is_compared() {
    let t = ClTerms::parse(SIGNED, &crate::annotate::cl_rules::shipped());
    let (mapped, unmapped) = t.map_labels([
        "CD4- T",
        "CD4-T",
        "CD34-",
        "CD34- cells",
        "CD4-CD8- thymocyte",
        "CD16 Monocytes",
        "Pro-B",
        "Pro B",
        "T-cell",
        "Tcells",
    ]);
    assert_eq!(
        mapped["Pro-B"], "CL:5",
        "a hyphen between letters parts words"
    );
    assert_eq!(mapped["Pro B"], "CL:5");
    assert_eq!(mapped["T-cell"], "CL:6");
    assert_eq!(
        mapped["Tcells"], "CL:6",
        "separators ignored as a last resort"
    );
    assert_eq!(
        unmapped,
        [
            "CD4- T",
            "CD4-T",
            "CD34-",
            "CD34- cells",
            "CD4-CD8- thymocyte",
            "CD16 Monocytes"
        ],
        "a negative marker is not the marker, and the marker is not its negative"
    );
}

#[test]
fn a_lineage_over_another_group_merges_nothing() {
    // lin M holds phagocytes and lin AB's T cells; a lone regulatory type
    // sits under lin M by another way down, beside lin AB, not in it.
    let obo = "[Term]\nid: CL:0\nname: root cell\n\n\
               [Term]\nid: CL:1\nname: lin M cell\nsubset: blood_and_immune_upper_slim\nis_a: CL:0\n\n\
               [Term]\nid: CL:2\nname: phago cell\nsubset: cellxgene_subset\nis_a: CL:1\n\n\
               [Term]\nid: CL:3\nname: mono A cell\nis_a: CL:2\n\n\
               [Term]\nid: CL:4\nname: mono B cell\nis_a: CL:2\n\n\
               [Term]\nid: CL:5\nname: mature T cell\nsubset: cellxgene_subset\nis_a: CL:1\n\n\
               [Term]\nid: CL:6\nname: lin AB cell\nsubset: blood_and_immune_upper_slim\nis_a: CL:5\n\n\
               [Term]\nid: CL:7\nname: ab one cell\nis_a: CL:6\n\n\
               [Term]\nid: CL:8\nname: ab two cell\nis_a: CL:6\n\n\
               [Term]\nid: CL:9\nname: reg T cell\nis_a: CL:5\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let p = panel(&[
        ("G1", "Mono A"),
        ("G2", "Mono B"),
        ("G3", "AB one"),
        ("G4", "AB two"),
        ("G5", "Reg T"),
    ]);
    let tree = TypeTree::from_ontology(&t, &p).unwrap();
    assert_eq!(members(&tree, "phago_cell"), ["Mono A", "Mono B"]);
    assert_eq!(members(&tree, "lin_AB_cell"), ["AB one", "AB two"]);
    assert_eq!(members(&tree, "reg_T_cell"), ["Reg T"], "on its own");
    assert_eq!(tree.groups.len(), 3, "one level: {:?}", tree.groups);
}

#[test]
fn a_negative_marker_s_word_is_no_panel_context() {
    let obo = "[Term]\nid: CL:1\nname: CD2-negative ya cell\n\n\
               [Term]\nid: CL:2\nname: CD3-negative za cell\n\n\
               [Term]\nid: CL:3\nname: CD5-negative wa cell\n";
    let mut t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let labels = ["CD2 ya", "CD3 za", "CD5 wa"];
    let (mapped, _) = t.map_labels(labels);
    assert!(mapped.is_empty(), "{mapped:?}");
    t.set_panel(&panel(&[
        ("G1", "CD2 ya"),
        ("G2", "CD3 za"),
        ("G3", "CD5 wa"),
    ]));
    assert_eq!(t.panel_context, Some(None));
}

#[test]
fn an_alias_for_one_label_leaves_the_panel_s_word_to_the_others() {
    let obo = "[Term]\nid: CL:1\nname: organx A cell\nsynonym: \"organx alpha cell\" EXACT []\n\n\
               [Term]\nid: CL:2\nname: organx beta cell\n\n\
               [Term]\nid: CL:3\nname: organx gamma cell\n";
    let mut t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    t.add_aliases("Alpha\tCL:1\n", "test");
    let labels = ["Alpha", "Beta", "Gamma"];
    let (mapped, unmapped) = t.map_labels(labels);
    assert!(unmapped.is_empty(), "{unmapped:?}");
    assert_eq!(mapped["Beta"], "CL:2");
    t.set_panel(&panel(&[("G1", "Alpha"), ("G2", "Beta"), ("G3", "Gamma")]));
    assert_eq!(t.panel_context, Some(Some("organx".into())));
}

#[test]
fn a_kind_of_a_type_the_panel_s_word_completes_is_placed_under_it() {
    let obo = "[Term]\nid: CL:1\nname: organx A cell\nsynonym: \"organx alpha cell\" EXACT []\n\n\
               [Term]\nid: CL:2\nname: organx beta cell\n\n\
               [Term]\nid: CL:3\nname: organx gamma cell\n\n\
               [Term]\nid: CL:11\nname: foo cell\n\n\
               [Term]\nid: CL:12\nname: mature foo cell\nis_a: CL:11\n\n\
               [Term]\nid: CL:21\nname: bar cell\n\n\
               [Term]\nid: CL:22\nname: mature bar cell\nis_a: CL:21\n\n\
               [Term]\nid: CL:31\nname: baz cell\n\n\
               [Term]\nid: CL:32\nname: mature baz cell\nis_a: CL:31\n";
    let labels = ["Alpha", "Beta", "Gamma", "Mature Beta"];
    let mut t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let placed = t.place_labels(labels);
    assert_eq!(placed.on["Beta"], "CL:2");
    assert_eq!(placed.under["Mature Beta"], "CL:2", "told from the labels");
    t.set_panel(&panel(&labels.map(|l| ("G", l))));
    let placed = t.place_labels(["Mature Beta"]);
    assert_eq!(
        placed.under["Mature Beta"], "CL:2",
        "read with the panel's word"
    );
}

#[test]
fn the_panel_s_word_completes_labels_with_word_order_kept() {
    let obo = "[Term]\nid: CL:1\nname: organx alpha cell\n\n\
               [Term]\nid: CL:2\nname: organx beta cell\n\n\
               [Term]\nid: CL:3\nname: organx gamma cell\n\n\
               [Term]\nid: CL:4\nname: cell blue\n";
    let mut rules = crate::annotate::cl_rules::shipped();
    rules.any_word_order = false;
    let t = ClTerms::parse(obo, &rules);
    let (mapped, unmapped) = t.map_labels(["Alpha", "Beta", "Gamma", "Blue cell"]);
    assert_eq!(mapped["Alpha"], "CL:1");
    assert_eq!(mapped.len(), 3);
    assert_eq!(unmapped, ["Blue cell"], "no other word order");
}

#[test]
fn a_type_under_a_term_is_placed_on_that_term_not_its_group_s() {
    let t = ClTerms::parse(MODIFIED, &crate::annotate::cl_rules::shipped());
    let p = panel(&[("G1", "CT5"), ("G2", "CT6"), ("G3", "Early CT5")]);
    let tree = TypeTree::from_ontology(&t, &p).unwrap();
    assert_eq!(members(&tree, "lin_B_cell"), ["CT5", "CT6", "Early CT5"]);
    let at = |l: &str| tree.placement(l).map(|p| (p.id, p.name, p.via));
    assert_eq!(
        at("Early CT5"),
        Some(("CL:5".into(), "Early CT5".into(), Via::Under))
    );
    assert_eq!(at("CT6"), Some(("CL:6".into(), "CT6".into(), Via::Label)));
}

#[test]
fn only_a_taxon_after_a_comma_is_a_qualifier() {
    let obo = "[Term]\nid: CL:1\nname: CT6 cell, human\n\n\
               [Term]\nid: CL:2\nname: CTQ cell, glycinergic\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let (mapped, unmapped) = t.map_labels(["CT6", "CTQ"]);
    assert_eq!(mapped["CT6"], "CL:1");
    assert_eq!(unmapped, ["CTQ"], "a kind is no qualifier");
}

#[test]
fn a_marker_s_positive_or_negative_word_reads_as_a_word() {
    let obo = "[Term]\nid: CL:1\nname: CD4-positive helper T cell\n\n\
               [Term]\nid: CL:2\nname: CD14-positive monocyte\n\n\
               [Term]\nid: CL:3\nname: CD16-negative monocyte\nsynonym: \"CD16- monocyte\" EXACT []\n\n\
               [Term]\nid: CL:4\nname: proCT9blast\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let (mapped, unmapped) = t.map_labels([
        "CD4 positive helper T cells",
        "CD4-positive helper T cells",
        "CD14positive monocyte",
        "CD16 negative monocytes",
        "CD16 Monocytes",
        "Pro-CT9blast",
    ]);
    assert_eq!(mapped["CD4 positive helper T cells"], "CL:1");
    assert_eq!(mapped["CD4-positive helper T cells"], "CL:1");
    assert_eq!(
        mapped["CD14positive monocyte"], "CL:2",
        "run together as a last resort"
    );
    assert_eq!(mapped["CD16 negative monocytes"], "CL:3");
    assert_eq!(
        unmapped,
        ["CD16 Monocytes", "Pro-CT9blast"],
        "the marker is not its negative; a hyphen beside a digit is kept, as in CD4-T"
    );
}

#[test]
fn an_alias_elsewhere_counts_for_no_panel_word() {
    let obo = "[Term]\nid: CL:1\nname: organx beta cell\n\n\
               [Term]\nid: CL:2\nname: organx gamma cell\n\n\
               [Term]\nid: CL:3\nname: organx delta cell\n\n\
               [Term]\nid: CL:11\nname: tissuey xa cell\n\n\
               [Term]\nid: CL:12\nname: tissuey xb cell\n\n\
               [Term]\nid: CL:13\nname: tissuey xc cell\n\n\
               [Term]\nid: CL:21\nname: first cell\n\n\
               [Term]\nid: CL:22\nname: second cell\n\n\
               [Term]\nid: CL:23\nname: third cell\n";
    let labels = ["Beta", "Gamma", "Delta", "Xa", "Xb", "Xc"];
    let mut t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    t.add_aliases("Xa\tCL:21\nXb\tCL:22\nXc\tCL:23\n", "test");
    let (mapped, unmapped) = t.map_labels(labels);
    assert!(unmapped.is_empty(), "{unmapped:?}");
    assert_eq!(
        mapped["Beta"], "CL:1",
        "the panel's word, not the aliased labels'"
    );
    assert_eq!(mapped["Xa"], "CL:21", "the alias");
    t.set_panel(&panel(&labels.map(|l| ("G", l))));
    assert_eq!(t.panel_context, Some(Some("organx".into())));
    // Aliased to the very terms the other word completes them to, they count.
    let mut t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    t.add_aliases("Xa\tCL:11\nXb\tCL:12\nXc\tCL:13\n", "test");
    let (mapped, _) = t.map_labels(labels);
    assert_eq!(mapped["Xa"], "CL:11");
    assert!(
        !mapped.contains_key("Beta"),
        "a tie tells no word: {mapped:?}"
    );
}

#[test]
fn a_run_of_words_a_name_matches_comes_before_one_the_panel_s_word_completes() {
    let obo = "[Term]\nid: CL:1\nname: organx alpha cell\n\n\
               [Term]\nid: CL:2\nname: organx beta cell\n\n\
               [Term]\nid: CL:3\nname: organx gamma cell\n\n\
               [Term]\nid: CL:4\nname: foo cell\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let placed = t.place_labels(["Alpha", "Beta", "Gamma", "Beta Foo"]);
    assert_eq!(placed.on["Beta"], "CL:2");
    assert_eq!(
        placed.hint.get("Beta Foo").map(String::as_str),
        Some("CL:4"),
        "{placed:?}"
    );
}
