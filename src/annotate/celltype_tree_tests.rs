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
    let (_, u1) = t1.map_labels(["Alpha_Type", "Endothelial"]);
    assert_eq!(
        u1,
        ["Alpha_Type", "Endothelial"],
        "not a synonym; not one word"
    );
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
    attach_by_sharing(&mut groups, vec!["X".to_string()], &pairs);
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
    attach_by_sharing(&mut groups, vec!["X".to_string()], &pairs);
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
fn two_types_sharing_a_marker_stay_one_group() {
    let tree = TypeTree::from_markers(&panel(&[("g", "A"), ("g", "B")]));
    assert_eq!(tree.groups.len(), 1);
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
fn an_implied_word_needs_a_label_of_two_or_more_words() {
    let obo = "[Term]\nid: CL:0\nname: cell\n\n\
               [Term]\nid: CL:1\nname: progenitor cell\n\n\
               [Term]\nid: CL:2\nname: early progenitor cell\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let (mapped, unmapped) = t.map_labels(["", "-", "Progenitor", "Early_Progenitor"]);
    assert_eq!(mapped["Early_Progenitor"], "CL:2");
    assert_eq!(
        unmapped,
        ["", "-", "Progenitor"],
        "no word, or one generic word"
    );
}

#[test]
fn a_terms_own_name_beats_another_terms_synonym_in_every_tier() {
    let obo = "[Term]\nid: CL:1\nname: pro CT8 cell\n\n\
               [Term]\nid: CL:2\nname: early CT8 cell\nsynonym: \"pro CT8 cell\" EXACT []\n";
    let t = ClTerms::parse(obo, &crate::annotate::cl_rules::shipped());
    let (mapped, _) = t.map_labels(["Pro_CT8", "Pro_CT8_cell", "Pro CT8 cells"]);
    for l in ["Pro_CT8", "Pro_CT8_cell", "Pro CT8 cells"] {
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
