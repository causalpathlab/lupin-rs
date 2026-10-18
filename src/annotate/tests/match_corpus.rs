//! The regression corpus for matching panel labels to terms and grouping
//! them: labels as panels write them, each with the term it must reach (or
//! none), on a made-up ontology shaped as the Cell Ontology is (marker signs,
//! comma qualifiers, species-qualified names, abbreviations, plural
//! synonyms, state classes over several lineages). Every matching rule is
//! checked against all of it, so a rule that fixes one label and breaks
//! another fails here.

use super::*;
use crate::annotate::cl_rules::{shipped, Aliases};

/// root ─┬─ CTH lineage (class) ─┬─ CT1 cell ─┬─ MK1-positive, alpha-beta CT1 cell ─┬─ MK1-positive helper CT1 cell
///       │                      │            │                                      └─ MK1-positive, alpha-beta CT1 cell, human
///       │                      │            ├─ MK5-positive, alpha-beta CT1 cell
///       │                      │            ├─ MK3-positive, MK4-negative CT1 cell
///       │                      │            ├─ MK3-negative, MK4-positive CT1 cell
///       │                      │            └─ regulatory CT1 cell ── MK1-positive, MK6-positive, alpha-beta regulatory CT1 cell
///       │                      ├─ CT2 cell ── memory CT2 cell ── MK2-negative memory CT2 cell
///       │                      └─ CS1 cell (class; also under CTS)
///       ├─ CTC lineage (class) ─┬─ CT7 cell
///       │                      └─ CS2 cell (also under CTS)
///       ├─ CTS cell (class: a state of several lineages)
///       ├─ CT8 precursor cell ── CT8 cell, human
///       └─ CT9 cell, mouse; CT9 cell, human
const OBO: &str = "format-version: 1.2
data-version: corpus/2026-01-01

[Term]
id: CL:9100000
name: root cell

[Term]
id: CL:9100001
name: CTH lineage cell
subset: general_cell_types_upper_slim
is_a: CL:9100000

[Term]
id: CL:9100002
name: CTC lineage cell
subset: general_cell_types_upper_slim
is_a: CL:9100000

[Term]
id: CL:9100003
name: CTS cell
subset: general_cell_types_upper_slim
subset: cellxgene_subset
is_a: CL:9100000

[Term]
id: CL:9100010
name: CT1 cell
subset: cellxgene_subset
is_a: CL:9100001

[Term]
id: CL:9100011
name: MK1-positive, alpha-beta CT1 cell
is_a: CL:9100010

[Term]
id: CL:9100018
name: MK1-positive, alpha-beta CT1 cell, human
is_a: CL:9100011

[Term]
id: CL:9100012
name: MK1-positive helper CT1 cell
synonym: \"MK1-positive CT1-helper cell\" EXACT []
is_a: CL:9100011

[Term]
id: CL:9100013
name: MK5-positive, alpha-beta CT1 cell
is_a: CL:9100010

[Term]
id: CL:9100014
name: MK3-positive, MK4-negative CT1 cell
is_a: CL:9100010

[Term]
id: CL:9100015
name: MK3-negative, MK4-positive CT1 cell
is_a: CL:9100010

[Term]
id: CL:9100016
name: regulatory CT1 cell
synonym: \"CTR\" EXACT []
is_a: CL:9100010

[Term]
id: CL:9100017
name: MK1-positive, MK6-positive, alpha-beta regulatory CT1 cell
synonym: \"CTR\" BROAD []
is_a: CL:9100016

[Term]
id: CL:9100020
name: CT2 cell
is_a: CL:9100001

[Term]
id: CL:9100021
name: memory CT2 cell
is_a: CL:9100020

[Term]
id: CL:9100022
name: MK2-negative memory CT2 cell
is_a: CL:9100021

[Term]
id: CL:9100030
name: CS1 cell
subset: blood_and_immune_upper_slim
synonym: \"CS1X\" RELATED OMO:0003000 []
is_a: CL:9100001
is_a: CL:9100003

[Term]
id: CL:9100031
name: CS2 cell
synonym: \"CS2X\" RELATED OMO:0003000 []
synonym: \"CS2 stromal cells\" RELATED OMO:0003004 []
is_a: CL:9100002
is_a: CL:9100003

[Term]
id: CL:9100040
name: CT7 cell
is_a: CL:9100002

[Term]
id: CL:9100019
name: MK7-A-positive CT1 cell
is_a: CL:9100010

[Term]
id: CL:9100060
name: CT10 cell
is_a: CL:9100000

[Term]
id: CL:9100061
name: CT10, human
is_a: CL:9100060

[Term]
id: CL:9100050
name: CT8 cell, human
is_a: CL:9100051

[Term]
id: CL:9100051
name: CT8 precursor cell
is_a: CL:9100000

[Term]
id: CL:9100052
name: CT9 cell, mouse
is_a: CL:9100000

[Term]
id: CL:9100053
name: CT9 cell, human
is_a: CL:9100000
";

/// Each label, the term it must reach (`None`: none), and why.
const CASES: &[(&str, Option<&str>, &str)] = &[
    // A term's name as written, and as panels space, pluralise and sign it.
    (
        "MK1-positive helper CT1 cell",
        Some("MK1-positive helper CT1 cell"),
        "as written",
    ),
    (
        "MK1 positive helper CT1 cells",
        Some("MK1-positive helper CT1 cell"),
        "a spaced sign is the hyphenated one",
    ),
    (
        "MK1+ helper CT1 cells",
        Some("MK1-positive helper CT1 cell"),
        "`+` is positive",
    ),
    (
        "helper CT1 cells MK1+",
        Some("MK1-positive helper CT1 cell"),
        "words in another order",
    ),
    (
        "MK1-positive CT1-helper cells",
        Some("MK1-positive helper CT1 cell"),
        "an exact synonym",
    ),
    (
        "MK2-negative memory CT2 cell",
        Some("MK2-negative memory CT2 cell"),
        "a negative sign as written",
    ),
    (
        "MK2 negative memory CT2 cells",
        Some("MK2-negative memory CT2 cell"),
        "a spaced negative sign",
    ),
    (
        "MK2- memory CT2 cells",
        Some("MK2-negative memory CT2 cell"),
        "a trailing `-` is negative",
    ),
    (
        "MK2-positive memory CT2 cell",
        None,
        "the opposite sign is another type",
    ),
    (
        "memory CT2 cells",
        Some("memory CT2 cell"),
        "no sign, no signed term",
    ),
    // Signs stay with their markers.
    (
        "MK4-positive MK3-negative CT1 cell",
        Some("MK3-negative, MK4-positive CT1 cell"),
        "each sign stays with its marker",
    ),
    (
        "MK3+ MK4- CT1 cells",
        Some("MK3-positive, MK4-negative CT1 cell"),
        "`+`/`-` stay with their markers",
    ),
    (
        "CT1 cells MK7-A+",
        Some("MK7-A-positive CT1 cell"),
        "a marker with a hyphen keeps it, its sign `+`",
    ),
    (
        "MK7-A positive CT1 cells",
        Some("MK7-A-positive CT1 cell"),
        "a marker with a hyphen, its sign spaced",
    ),
    // Abbreviations and synonyms.
    (
        "CTR",
        Some("regulatory CT1 cell"),
        "an exact synonym; a broad one does not count",
    ),
    (
        "MK1-positive CTR",
        None,
        "a signed abbreviation names no term; dropping the sign would be a guess",
    ),
    ("CS1X", Some("CS1 cell"), "an abbreviation"),
    ("CS2X", Some("CS2 cell"), "an abbreviation"),
    // Comma qualifiers: the most general term the label fits.
    (
        "MK5-positive CT1 cells",
        Some("MK5-positive, alpha-beta CT1 cell"),
        "the label leaves off a qualifier",
    ),
    (
        "MK1-positive CT1 cells",
        Some("MK1-positive, alpha-beta CT1 cell"),
        "of the terms it fits, the most general",
    ),
    (
        "MK5+ CT1",
        Some("MK5-positive, alpha-beta CT1 cell"),
        "a qualifier and the trailing word both left off",
    ),
    (
        "MK1-positive alpha beta CT1 cells",
        Some("MK1-positive, alpha-beta CT1 cell"),
        "a qualifier spaced as a label writes it",
    ),
    // Species-qualified names.
    (
        "CT8 cells",
        Some("CT8 cell, human"),
        "the only term the label fits, but for its species",
    ),
    (
        "CT9 cell",
        None,
        "two species' terms, neither above the other",
    ),
    (
        "CT10",
        Some("CT10 cell"),
        "completing a term's own name comes before leaving off a qualifier",
    ),
    // Plural synonyms.
    (
        "CS2 stromal cells",
        Some("CS2 cell"),
        "a plural-form synonym",
    ),
    // What must stay unmatched.
    ("CT6", None, "no such term"),
    ("cell", None, "a word every term has"),
];

fn terms() -> ClTerms {
    ClTerms::parse(OBO, &shipped())
}

#[test]
fn every_corpus_label_reaches_its_term() {
    let t = terms();
    let mut wrong = Vec::new();
    for &(label, want, why) in CASES {
        let (mapped, _) = t.map_labels([label]);
        let got = mapped.get(label).and_then(|id| t.name(id));
        if got != want {
            wrong.push(format!("{label:?}: got {got:?}, want {want:?} ({why})"));
        }
    }
    assert!(wrong.is_empty(), "\n{}", wrong.join("\n"));
}

#[test]
fn a_marker_sign_is_one_word_however_it_is_written() {
    let r = shipped();
    for spelled in ["MK2-negative", "MK2 negative", "MK2-", "mk2-NEGATIVE"] {
        assert_eq!(r.normalise(spelled), "mk2-negative", "{spelled:?}");
    }
    for spelled in ["MK1-positive", "MK1 positive", "MK1+"] {
        assert_eq!(r.normalise(spelled), "mk1-positive", "{spelled:?}");
    }
    assert_eq!(
        r.normalise("MK1+MK6+ CT1"),
        "mk1-positive mk6-positive ct1",
        "signs written together part"
    );
    assert_eq!(
        r.normalise("MK1+, MK6- CT1"),
        "mk1-positive mk6-negative ct1"
    );
    for label in [
        "MK2-negative memory CT2 cells",
        "MK1+ helper",
        "pro-CT2",
        "CT1-helper cell",
    ] {
        let once = r.normalise(label);
        assert_eq!(
            r.normalise(&once),
            once,
            "normalising twice changes nothing: {label:?}"
        );
    }
    assert_eq!(
        r.normalise("pro-CT2"),
        "pro-ct2",
        "a hyphen before a word is no sign"
    );
    for spelled in ["MK7-A+", "MK7-A positive", "MK7-A-positive"] {
        assert_eq!(r.normalise(spelled), "mk7-a-positive", "{spelled:?}");
    }
    assert_eq!(r.normalise("MK7-A-"), "mk7-a-negative");
}

#[test]
fn an_alias_keeps_its_term_over_every_rule() {
    let mut aliases = Aliases::default();
    aliases.add_tsv("MK1-positive CTR\tCL:9100017\n", "test");
    let t = terms().with_aliases(aliases);
    let (mapped, _) = t.map_labels(["MK1-positive CTR", "MK1+ CTR", "CTR"]);
    assert_eq!(mapped["MK1-positive CTR"], "CL:9100017");
    assert_eq!(
        mapped["CTR"], "CL:9100016",
        "the alias is for its own label"
    );
    assert_eq!(
        mapped["MK1+ CTR"], "CL:9100017",
        "an alias is compared as names are, its sign however written"
    );
}

#[test]
fn a_later_alias_table_wins_over_an_earlier_spelling_of_its_label() {
    let mut t = terms();
    t.add_aliases("MK1-positive CTR\tCL:9100017\n", "earlier");
    t.add_aliases("MK1+ CTR\tCL:9100016\n", "later");
    let (mapped, _) = t.map_labels(["MK1-positive CTR"]);
    assert_eq!(mapped["MK1-positive CTR"], "CL:9100016");
    t.add_aliases("MK1 positive CTR\tCL:9100017\n", "latest");
    let (mapped, _) = t.map_labels(["MK1+ CTR"]);
    assert_eq!(mapped["MK1+ CTR"], "CL:9100017");
}

#[test]
fn a_search_for_a_word_being_typed_finds_signed_terms() {
    let t = terms();
    let found: Vec<&str> = (t.search("MK1-", 20).iter())
        .filter_map(|id| t.name(id))
        .collect();
    assert!(
        found.contains(&"MK1-positive helper CT1 cell"),
        "a trailing `-` is not yet a sign: {found:?}"
    );
    assert!(t.search("+", 5).is_empty());
}

#[test]
fn rule_words_are_compared_as_names_are() {
    let mut rules = shipped();
    rules.implied.suffixes = vec!["Cells".into()];
    rules.qualifiers = vec!["Alpha-Beta".into()];
    let t = ClTerms::parse(OBO, &rules);
    let (mapped, unmapped) = t.map_labels(["memory CT2", "MK5-positive CT1 cells"]);
    assert!(unmapped.is_empty(), "{unmapped:?}");
    assert_eq!(t.name(&mapped["memory CT2"]), Some("memory CT2 cell"));
    assert_eq!(
        t.name(&mapped["MK5-positive CT1 cells"]),
        Some("MK5-positive, alpha-beta CT1 cell")
    );
}

/// The shipped rules, with the corpus's state class.
fn tree(types: &[&str]) -> TypeTree {
    let mut rules = shipped();
    rules.classes.states = vec!["CL:9100003".into()];
    let panel: Vec<(String, String)> = types
        .iter()
        .enumerate()
        .map(|(i, t)| (format!("GENE{i}"), (*t).to_string()))
        .collect();
    TypeTree::from_ontology(&ClTerms::parse(OBO, &rules), &panel).expect("the types match terms")
}

#[test]
fn a_state_class_does_not_group_types_of_two_lineages() {
    let t = tree(&["CS1X", "CS2X", "CT2 cell", "CT7 cell"]);
    assert_eq!(t.group_of("CS1X"), t.group_of("CT2 cell"), "by lineage");
    assert_eq!(t.group_of("CS2X"), t.group_of("CT7 cell"), "by lineage");
    assert_ne!(t.group_of("CS1X"), t.group_of("CS2X"));
    let alone = tree(&["CS1X", "CS2X"]);
    assert_ne!(
        alone.group_of("CS1X"),
        alone.group_of("CS2X"),
        "not even when nothing else shares a class with them"
    );
}
