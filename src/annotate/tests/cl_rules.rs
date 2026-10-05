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
    assert!(!r.classes.is_class(&["cellxgene_subset".into()]));
}

#[test]
fn the_shipped_rules_reproduce_the_matching_lupin_had() {
    let shipped: Value =
        serde_json::from_str(include_str!("../../../data/cl_matching.json")).unwrap();
    let r = MatchRules::from_layers(&[shipped]).unwrap();
    assert!(r.counts("RELATED", &["OMO:0003000"]));
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
    assert!(r.classes.is_class(&["blood_and_immune_upper_slim".into()]));
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
    assert_eq!(a.len(), 2);
    // As `--label-cl` maps have been written: commas, a header, an NA row.
    let n = a.add_tsv(
        "cell_type,cl_id\nCT1 x,CL:9000001\nCT9\tNA\nbroken line\n",
        "x",
    );
    assert_eq!(n, 1, "only the CL pair is taken");
    assert_eq!(a.get("ct1_x"), Some("CL:9000001"));
}

#[test]
fn no_alias_table_ships_with_lupin() {
    // Label aliases are data about a marker panel: they travel with it.
    let data = concat!(env!("CARGO_MANIFEST_DIR"), "/data");
    assert!(!std::path::Path::new(data).join("cl_aliases.tsv").exists());
}
