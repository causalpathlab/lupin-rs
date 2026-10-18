//! [`super`]: the rules as data, layered.

use super::*;
use serde_json::json;

#[test]
fn without_a_rules_file_matching_is_literal() {
    let r = MatchRules::from_layers(&[]).unwrap();
    assert!(r.counts("EXACT", &[]));
    assert!(!r.counts("RELATED", &["OMO:0003000"]), "no abbreviations");
    assert_eq!(r.normalise("T Cells"), "t cells", "no plural rule");
    assert_eq!(r.singular("t cells"), "t cells");
    assert!(!r.any_word_order);
    assert!(r.implied.suffixes.is_empty(), "no word is implied");
    assert!(!r.classes.is_class(&["cellxgene_subset".into()]));
}

#[test]
fn the_shipped_rules_reproduce_the_matching_lupin_had() {
    let shipped: Value =
        serde_json::from_str(include_str!("../../../data/cl_matching.json")).unwrap();
    let r = MatchRules::from_layers(&[shipped]).unwrap();
    assert!(r.counts("RELATED", &["OMO:0003000"]));
    assert!(
        r.counts("RELATED", &["OMO:0003004"]),
        "a plural form (OMO:0003004) counts"
    );
    assert!(
        !r.counts("RELATED", &[]),
        "a plain RELATED synonym does not count"
    );
    assert!(!r.counts("NARROW", &[]));
    assert_eq!(r.normalise("CT1 Cells"), "ct1 cell");
    assert_eq!(r.singular("ct1 cells memory"), "ct1 cell memory");
    assert_eq!(
        r.singular("class ms"),
        "class ms",
        "not `ss`, not short words"
    );
    assert!(r.any_word_order);
    assert_eq!(r.implied.suffixes, ["cell"]);
    assert!(r.classes.is_class(&["blood_and_immune_upper_slim".into()]));
    assert!(r
        .classes
        .is_lineage(&["blood_and_immune_upper_slim".into()]));
    assert!(
        !r.classes.is_lineage(&[
            "general_cell_types_upper_slim".into(),
            "cellxgene_subset".into()
        ]),
        "a whole system is no lineage"
    );
    assert!(!r
        .classes
        .is_class(&["cellxgene_subset".into(), "upper_level".into()]));
}

#[test]
fn a_later_layer_overrides_only_what_it_sets() {
    let base = json!({"any_word_order": true, "plural": {"singular_min_len": 4, "suffixes": [["cells", "cell"]]}});
    let user = json!({"plural": {"singular_min_len": null}, "_doc": "mine"});
    let r = MatchRules::from_layers(&[base, user]).unwrap();
    assert!(r.any_word_order, "kept");
    assert_eq!(r.plural.singular_min_len, None, "null removes");
    assert_eq!(
        r.plural.suffixes,
        [("cells".to_string(), "cell".to_string())]
    );
}

#[test]
fn aliases_layer_and_compare_as_labels() {
    let mut a = Aliases::default();
    a.add_tsv(
        "# note\nlabel\tcl_id\tnote\nPro B\tCL:1\tbuilt-in\nNK\tCL:2\n",
        "built-in",
    );
    a.add_tsv("pro_b\tCL:9\tmine\n", "project");
    assert_eq!(a.get("PRO B"), Some("CL:9"), "the later file wins");
    assert_eq!(a.get("nk"), Some("CL:2"));
    assert_eq!(a.get("Pro-B"), Some("CL:9"), "a hyphen parts words too");
    assert_eq!(a.len(), 2);
    // As `--label-cl` maps have been written: commas, a header, an NA row.
    let n = a.add_tsv(
        "cell_type,cl_id\nCT1 x,CL:9000001\nCT9\tNA\nbroken line\n\"CT2, y\",CL:9000002\n",
        "x",
    );
    assert_eq!(n, 2, "only the CL pairs are taken");
    assert_eq!(a.get("ct1_x"), Some("CL:9000001"));
    assert_eq!(
        a.get("CT2 y"),
        Some("CL:9000002"),
        "a quoted label keeps its comma"
    );
    a.add_tsv("CT3\tCL:9000003, a note\n", "x");
    assert_eq!(
        a.get("CT3"),
        Some("CL:9000003"),
        "a note after the id is not part of it"
    );
}

#[test]
fn a_negative_markers_hyphen_is_kept_in_an_alias() {
    let mut a = Aliases::default();
    a.add_tsv(
        "CD34\tCL:1\nCD34-\tCL:2\nCD4 T\tCL:3\nCD4- T\tCL:4\nPro-B\tCL:5\n",
        "project",
    );
    assert_eq!(a.get("cd34"), Some("CL:1"));
    assert_eq!(a.get("CD34-"), Some("CL:2"), "not CD34's");
    assert_eq!(a.get("CD4_T"), Some("CL:3"));
    assert_eq!(a.get("CD4-_T"), Some("CL:4"), "not CD4 T's");
    assert_eq!(
        a.get("pro b"),
        Some("CL:5"),
        "a hyphen between two words parts them"
    );
}

#[test]
fn a_negative_markers_hyphen_before_the_next_marker_is_kept() {
    let mut a = Aliases::default();
    a.add_tsv(
        "CD34 CD38-\tCL:1\nCD34-CD38-\tCL:2\nCD4 CD8-\tCL:3\nCD4-CD8-\tCL:4\n\
         CD45RA-CD38-\tCL:5\npre B\tCL:6\n",
        "project",
    );
    assert_eq!(a.get("CD34 CD38-"), Some("CL:1"));
    assert_eq!(a.get("CD34-CD38-"), Some("CL:2"), "not CD34 CD38-'s");
    assert_eq!(a.get("CD4 CD8-"), Some("CL:3"));
    assert_eq!(a.get("CD4-CD8-"), Some("CL:4"), "not CD4 CD8-'s");
    assert_eq!(
        a.get("CD45RA CD38-"),
        None,
        "a marker with letters after its digits"
    );
    assert_eq!(
        a.get("pre-B"),
        Some("CL:6"),
        "a word without digits is no marker"
    );
}

#[test]
fn a_hyphen_parts_two_words_of_letters_only() {
    let same = |x: &str, y: &str| assert_eq!(key(x), key(y), "{x} is {y}");
    let differ = |x: &str, y: &str| assert_ne!(key(x), key(y), "{x} is not {y}");
    same("Pro-B", "Pro B");
    same("pre-B", "pre B");
    same("T-cell", "T cell");
    same("alpha-beta T cell", "alpha beta T cell");
    same("Pre-Pro-B", "Pre Pro B");
    // A marker's own words for its sign, as the ontology writes them.
    same("CD4-positive T cell", "CD4 positive T cell");
    same("CD16-negative monocyte", "CD16 negative monocyte");
    same("CD4-CD8-negative", "CD4-CD8 negative");
    // A side with a digit, a sign at a word's end: kept as written.
    for literal in [
        "CD34-",
        "CD34-CD38-",
        "CD4- T",
        "NK-CD56bright",
        "CD27-IgD-",
        "B-1a",
        "CD8-T",
    ] {
        assert_eq!(
            key(literal),
            label_key(&literal.to_lowercase()),
            "{literal}"
        );
    }
    differ("CD34-", "CD34");
    differ("CD34-CD38-", "CD34 CD38-");
    differ("CD4- T", "CD4 T");
    differ("CD4-CD8-", "CD4 CD8");
    differ("CD4-CD8", "CD4 CD8");
    differ("NK-CD56bright", "NK CD56bright");
    differ("B-1a", "B 1a");
    differ("CD16- monocyte", "CD16 monocyte");
    differ("CD4-positive T cell", "CD4 T cell");
}
