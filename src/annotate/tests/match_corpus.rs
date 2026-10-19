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
///       ├─ CT9 cell, mouse; CT9 cell, human
///       ├─ TISA K1 cell ("TISA omega cell"), TISA K2 cell ("TISA sigma cell"), TISA tau cell
///       ├─ omega CTG cell ("omega cell")
///       ├─ CT11/CT12/CT13 cell, each with "mature …" under it; CT14 cell, and "pre CT14 cell" not under it
///       └─ (under CTH lineage) CT16 cell, CT17 cell
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
id: CL:9100070
name: TISA K1 cell
synonym: \"TISA omega cell\" EXACT []
is_a: CL:9100000

[Term]
id: CL:9100071
name: TISA K2 cell
synonym: \"TISA sigma cell\" EXACT []
is_a: CL:9100000

[Term]
id: CL:9100072
name: TISA tau cell
is_a: CL:9100000

[Term]
id: CL:9100073
name: omega CTG cell
synonym: \"omega cell\" EXACT []
is_a: CL:9100000

[Term]
id: CL:9100080
name: CT11 cell
is_a: CL:9100000

[Term]
id: CL:9100081
name: mature CT11 cell
is_a: CL:9100080

[Term]
id: CL:9100082
name: CT12 cell
is_a: CL:9100000

[Term]
id: CL:9100083
name: mature CT12 cell
is_a: CL:9100082

[Term]
id: CL:9100084
name: CT13 cell
is_a: CL:9100000

[Term]
id: CL:9100085
name: mature CT13 cell
is_a: CL:9100084

[Term]
id: CL:9100086
name: CT14 cell
is_a: CL:9100000

[Term]
id: CL:9100087
name: pre CT14 cell
is_a: CL:9100000

[Term]
id: CL:9100088
name: CT16 cell
is_a: CL:9100001

[Term]
id: CL:9100089
name: CT17 cell
is_a: CL:9100001

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

#[test]
fn a_word_the_panel_leaves_off_completes_its_labels() {
    let t = terms();
    let (mapped, unmapped) = t.map_labels(["Omega", "Sigma", "Tau", "CT1 cell"]);
    assert!(unmapped.is_empty(), "{unmapped:?}");
    let name = |l: &str| t.name(&mapped[l]);
    assert_eq!(name("Omega"), Some("TISA K1 cell"), "through a synonym");
    assert_eq!(name("Sigma"), Some("TISA K2 cell"));
    assert_eq!(name("Tau"), Some("TISA tau cell"));
    let (_, unmapped) = t.map_labels(["Omega", "Sigma"]);
    assert_eq!(
        unmapped,
        ["Omega", "Sigma"],
        "too few labels to tell the word"
    );
}

/// A panel whose types each have their own marker, and `shared` markers.
fn tree_of(types: &[&str], shared: &[(&str, &str)]) -> TypeTree {
    let mut panel: Vec<(String, String)> = (types.iter().enumerate())
        .map(|(i, t)| (format!("GENE{i}"), (*t).to_string()))
        .collect();
    panel.extend(
        shared
            .iter()
            .map(|(g, t)| ((*g).to_string(), (*t).to_string())),
    );
    TypeTree::from_ontology(&terms(), &panel).expect("the types match terms")
}

#[test]
fn a_label_nothing_matches_goes_under_the_term_its_longest_run_names() {
    let t = tree_of(&["Mature CT16", "CT17 cell", "CT2 cell", "CT7 cell"], &[]);
    assert_eq!(
        t.group_of("Mature CT16"),
        t.group_of("CT17 cell"),
        "`mature` is a word the ontology uses for kinds of a term"
    );
    assert!(
        !t.label_cl.contains_key("Mature CT16"),
        "placed, not matched"
    );
}

#[test]
fn a_run_with_other_words_left_keeps_the_label_to_its_lineage() {
    // `Pre CT16` shares more markers with CT7 (another lineage) than with CT17.
    let shared = [
        ("GENEA", "Pre CT16"),
        ("GENEA", "CT7 cell"),
        ("GENEB", "Pre CT16"),
        ("GENEB", "CT7 cell"),
        ("GENEC", "Pre CT16"),
        ("GENEC", "CT17 cell"),
    ];
    let t = tree_of(&["Pre CT16", "CT17 cell", "CT2 cell", "CT7 cell"], &shared);
    assert_eq!(
        t.group_of("Pre CT16"),
        t.group_of("CT17 cell"),
        "`pre` is no kind word, but CT16 keeps the label in its lineage"
    );
}

#[test]
fn a_word_left_off_may_be_the_implied_one_in_any_place() {
    let t = tree_of(&["CT16 mature cells", "CT17 cell", "CT2 cell"], &[]);
    assert_eq!(
        t.group_of("CT16 mature cells"),
        t.group_of("CT17 cell"),
        "`cell` left off is no word against placing"
    );
    let obo = "[Term]\nid: CL:1\nname: cell\n\n[Term]\nid: CL:2\nname: CT1 cell\nis_a: CL:1\n";
    let t = ClTerms::parse(obo, &shipped());
    assert_eq!(
        t.place_with("CT1 odd cells", None)
            .map(|(id, _)| id)
            .as_deref(),
        Some("CL:2"),
        "a run of implied words alone names nothing"
    );
}

#[test]
fn the_panel_word_completes_runs_too() {
    let t = tree_of(&["Omega", "Sigma", "Tau", "Mature Omega", "CT2 cell"], &[]);
    assert_eq!(t.group_of("Mature Omega"), t.group_of("Omega"));
}

#[test]
fn placed_labels_count_toward_a_tree() {
    let t = tree_of(&["CT2 cell", "Mature CT16", "Mature CT17"], &[]);
    assert_eq!(t.group_of("Mature CT16"), t.group_of("Mature CT17"));
}

#[test]
fn a_label_kept_to_a_lineage_and_alone_keeps_its_term() {
    let shared = [("GENEA", "Pre CT16"), ("GENEA", "CT7 cell")];
    let t = tree_of(&["Pre CT16", "CT17 cell", "CT2 cell", "CT7 cell"], &shared);
    let own = t.groups.iter().find(|g| g.members == ["Pre CT16"]);
    assert_eq!(
        own.and_then(|g| g.cl_id.as_deref()),
        Some("CL:9100088"),
        "{:?}",
        t.groups
    );
}

#[test]
fn the_panel_tree_puts_a_placed_label_with_its_lineage() {
    let panel: Vec<(String, String)> = ["Mature CT16", "CT17 cell", "CT7 cell", "CT10 cell"]
        .iter()
        .enumerate()
        .map(|(i, t)| (format!("GENE{i}"), (*t).to_string()))
        .collect();
    let p = crate::annotate::panel_tree::PanelTree::from_ontology(&terms(), &panel).unwrap();
    let parent = |i: usize| p.nodes[i].parent;
    let up = |l: &str| {
        let mut at = p.node_of(l);
        std::iter::from_fn(move || {
            at = at.and_then(parent);
            at
        })
        .collect::<Vec<_>>()
    };
    let lineage = parent(p.node_of("CT17 cell").unwrap()).expect("CT17 under its lineage");
    assert!(
        up("Mature CT16").contains(&lineage),
        "under its term, in the lineage it shares with CT17: {:?}",
        up("Mature CT16")
    );
}

#[test]
fn a_panel_reading_says_how_each_label_was_placed() {
    let mut aliases = Aliases::default();
    aliases.add_tsv("CTR\tCL:9100016\n", "test");
    let t = terms().with_aliases(aliases);
    let p = t.place_labels([
        "CTR",
        "CT2 cell",
        "Omega",
        "Sigma",
        "Tau",
        "Mature CT16",
        "Pre CT16",
        "CT6",
    ]);
    let how = |l: &str| p.how(l);
    assert_eq!(how("CTR"), How::Alias);
    assert_eq!(how("CT2 cell"), How::Name);
    assert_eq!(how("Omega"), How::PanelWord("tisa".into()));
    assert_eq!(how("Mature CT16"), How::Under);
    assert_eq!(how("Pre CT16"), How::Kept);
    assert_eq!(how("CT6"), How::Nowhere);
}

/// Two types under CT1 cell (a class) and CT2 cell, which shares only the
/// broader CTH lineage with them.
const LONE: [&str; 3] = [
    "MK1-positive, alpha-beta CT1 cell",
    "MK5-positive, alpha-beta CT1 cell",
    "CT2 cell",
];

#[test]
fn a_type_its_siblings_left_alone_is_a_question() {
    let t = tree_of(&LONE, &[]);
    assert_ne!(t.group_of("CT2 cell"), t.group_of(LONE[0]));
    let panel: Vec<(String, String)> = (LONE.iter().enumerate())
        .map(|(i, l)| (format!("GENE{i}"), (*l).to_string()))
        .collect();
    let terms = terms();
    let asked = TypeTree::join_questions(&terms, &panel);
    assert_eq!(
        asked,
        [JoinQuestion {
            label: "CT2 cell".into(),
            class: "CL:9100001".into(),
            group: t.group_of(LONE[0]).unwrap().to_string(),
        }]
    );
}

#[test]
fn an_answer_joins_the_lone_type_under_the_class_it_shares() {
    let panel: Vec<(String, String)> = (LONE.iter().enumerate())
        .map(|(i, l)| (format!("GENE{i}"), (*l).to_string()))
        .collect();
    let mut terms = terms();
    terms.add_joins("CT2 cell\tjoin\n");
    let t = TypeTree::from_ontology(&terms, &panel).unwrap();
    assert_eq!(t.group_of("CT2 cell"), t.group_of(LONE[0]));
    let g = t.groups.iter().find(|g| g.members.len() == 3).unwrap();
    assert_eq!(
        g.cl_id.as_deref(),
        Some("CL:9100001"),
        "named for the class they share"
    );
    assert!(
        TypeTree::join_questions(&terms, &panel).is_empty(),
        "answered"
    );
    let mut terms = self::terms();
    terms.add_joins("CT2 cell\talone\n");
    let t = TypeTree::from_ontology(&terms, &panel).unwrap();
    assert_ne!(t.group_of("CT2 cell"), t.group_of(LONE[0]));
    assert!(
        TypeTree::join_questions(&terms, &panel).is_empty(),
        "answered"
    );
}

#[test]
fn the_record_drafts_a_methods_paragraph() {
    let mut terms = terms();
    terms.add_aliases("CTR\tCL:9100016\n", "aliases");
    terms.add_joins("CT2 cell\tjoin\n");
    let labels = [
        "CTR",
        "MK1-positive, alpha-beta CT1 cell",
        "MK5-positive, alpha-beta CT1 cell",
        "CT2 cell",
        "Omega",
        "Sigma",
        "Tau",
        "Mature CT14",
        "CT6",
    ];
    let panel: Vec<(String, String)> = (labels.iter().enumerate())
        .map(|(i, l)| (format!("GENE{i}"), (*l).to_string()))
        .collect();
    let r = crate::annotate::provenance::Record::of(&terms, &panel);
    assert_eq!(r.labels.len(), 9);
    assert_eq!(r.answers, [("CT2 cell".to_string(), true)]);
    let text = crate::annotate::provenance::methods(&r);
    for want in [
        "9 cell-type labels",
        "release corpus/2026-01-01",
        "1 by a curated alias",
        "3 by completing the word \"tisa\"",
        "1 was placed under",
        "1 matched no term",
        "CT2 cell",
    ] {
        assert!(text.contains(want), "{want:?} not in:\n{text}");
    }
    assert_eq!(
        text,
        crate::annotate::provenance::methods(&r),
        "the same record, the same text"
    );
}

#[test]
fn an_answer_is_read_in_any_case_and_only_answers_that_apply_are_recorded() {
    let mut terms = terms();
    assert_eq!(
        terms.add_joins("CT2 cell\tJoin\nCT7 cell\tjoin\nCT1 cell\tmaybe\n"),
        2
    );
    let panel: Vec<(String, String)> = (LONE.iter().chain(["CT7 cell"].iter()).enumerate())
        .map(|(i, l)| (format!("GENE{i}"), (*l).to_string()))
        .collect();
    let r = crate::annotate::provenance::Record::of(&terms, &panel);
    assert_eq!(
        r.answers,
        [("CT2 cell".to_string(), true)],
        "CT7 is not left alone here: its answer did nothing"
    );
}

#[test]
fn the_methods_paragraph_says_when_groups_came_from_markers_and_reads_without_a_release() {
    let obo = "[Term]\nid: CL:1\nname: CT1 cell\n";
    let terms = ClTerms::parse(obo, &shipped());
    let panel: Vec<(String, String)> = ["CT1 cell", "XQ", "YZ"]
        .iter()
        .enumerate()
        .map(|(i, l)| (format!("GENE{i}"), (*l).to_string()))
        .collect();
    let text = crate::annotate::provenance::methods(&crate::annotate::provenance::Record::of(
        &terms, &panel,
    ));
    assert!(!text.contains("release of unknown"), "{text}");
    assert!(text.contains("shared marker genes"), "{text}");
    assert!(text.contains("too few"), "{text}");
}
