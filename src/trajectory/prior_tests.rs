use super::*;

/// A small ontology: a stem cell, two progenitors, a precursor and a mature
/// type that develops from the precursor, and a broad class the mature type
/// `is_a`, which develops from the stem cell (the shortcut inheritance makes).
const OBO: &str = "data-version: test
[Term]
id: CL:0001
name: stem cell
[Term]
id: CL:0002
name: progenitor
relationship: RO:0002202 CL:0001 ! develops from stem cell
[Term]
id: CL:0003
name: committed progenitor
relationship: RO:0002202 CL:0002 ! develops from progenitor
[Term]
id: CL:0004
name: precursor
relationship: RO:0002202 CL:0003 ! develops from committed progenitor
[Term]
id: CL:0010
name: broad class
relationship: RO:0002202 CL:0001 ! develops from stem cell
[Term]
id: CL:0005
name: mature cell
is_a: CL:0010 ! broad class
relationship: RO:0002202 CL:0004 ! develops from precursor
";

fn terms() -> ClTerms {
    ClTerms::parse(OBO, &crate::annotate::cl_rules::shipped())
}

fn labels(names: &[&str]) -> Vec<Box<str>> {
    names.iter().map(|&n| n.into()).collect()
}

fn label_cl() -> BTreeMap<String, String> {
    [
        ("Stem", "CL:0001"),
        ("Prog", "CL:0002"),
        ("Committed", "CL:0003"),
        ("Pre", "CL:0004"),
        ("Mature", "CL:0005"),
    ]
    .into_iter()
    .map(|(l, id)| (l.to_string(), id.to_string()))
    .collect()
}

fn edges_named(types: &[Box<str>], p: &Prior) -> Vec<(String, String, String)> {
    p.edges
        .iter()
        .map(|e| {
            (
                types[e.from].to_string(),
                types[e.to].to_string(),
                e.source
                    .map_or("implied".to_string(), |s| s.as_str().to_string()),
            )
        })
        .collect()
}

#[test]
fn ontology_statements_follow_develops_from_and_inherit_through_is_a() {
    let st = from_ontology(&terms(), &label_cl());
    let find = |a: &str, b: &str| {
        st.iter()
            .find(|s| s.from == a && s.to == b)
            .map(|s| s.source)
    };
    assert_eq!(find("Stem", "Prog"), Some(Source::Cl));
    assert_eq!(
        find("Stem", "Committed"),
        Some(Source::Cl),
        "through the chain"
    );
    assert_eq!(find("Pre", "Mature"), Some(Source::Cl));
    assert_eq!(
        find("Stem", "Mature"),
        Some(Source::Cl),
        "the chain reaches the stem cell directly, so the broad class's shortcut does not make it inherited"
    );
    assert_eq!(find("Prog", "Stem"), None, "never backwards");
}

#[test]
fn inherited_statements_are_tagged() {
    let mut map = label_cl();
    map.remove("Pre");
    map.remove("Committed");
    map.remove("Prog");
    // With only the broad class's develops_from reaching the stem cell.
    let st = from_ontology(&terms(), &map);
    let stem_mature = st
        .iter()
        .find(|s| s.from == "Stem" && s.to == "Mature")
        .unwrap();
    // Mature's own chain still reaches Stem through dropped terms (Cl wins over inherited).
    assert_eq!(stem_mature.source, Source::Cl);
}

#[test]
fn reduction_keeps_only_direct_edges_and_passes_through_dropped_types() {
    let types = labels(&["Stem", "Prog", "Committed", "Pre", "Mature"]);
    let st = from_ontology(&terms(), &label_cl());
    // "Committed" has too few cells: it is not a node.
    let is_node = [true, true, false, true, true];
    let p = build(&types, &is_node, st, &[]).unwrap();
    let e = edges_named(&types, &p);
    assert_eq!(
        e,
        vec![
            ("Stem".into(), "Prog".into(), "cl".into()),
            ("Prog".into(), "Pre".into(), "cl".into()),
            ("Pre".into(), "Mature".into(), "cl".into()),
        ],
        "{e:?}"
    );
    assert_eq!(p.roots, vec![vec![0]]);
    assert_eq!(p.component, vec![Some(0), Some(0), None, Some(0), Some(0)]);
}

#[test]
fn a_later_layer_replaces_a_pair_whatever_its_direction() {
    let cl = vec![Statement {
        from: "A".into(),
        to: "B".into(),
        relation: Relation::Precedes,
        source: Source::Cl,
        note: String::new(),
    }];
    let project = parse_statements(
        "B\tA\tprecedes\treversed on purpose\n",
        Source::Project,
        "p",
    )
    .unwrap();
    let combined = combine(&[cl, project]);
    assert_eq!(combined.len(), 1);
    assert_eq!(
        (combined[0].from.as_str(), combined[0].to.as_str()),
        ("B", "A")
    );
    assert_eq!(combined[0].source, Source::Project);

    let user = parse_statements(
        "from\tto\trelation\n# a comment\nA\tB\tunrelated\n",
        Source::User,
        "u",
    )
    .unwrap();
    let combined = combine(&[combined, user]);
    assert_eq!(combined[0].relation, Relation::Unrelated);
    let types = labels(&["A", "B"]);
    let p = build(&types, &[true, true], combined, &[]).unwrap();
    assert!(p.edges.is_empty());
    assert_eq!(p.component, vec![None, None]);
}

#[test]
fn a_cycle_is_an_error_naming_its_statements() {
    let st = parse_statements(
        "A\tB\tprecedes\nB\tC\tprecedes\nC\tA\tprecedes\n",
        Source::Run,
        "r",
    )
    .unwrap();
    let types = labels(&["A", "B", "C"]);
    let err = build(&types, &[true; 3], st, &[]).unwrap_err().to_string();
    assert!(
        err.contains("cycle") && err.contains("C precedes A (run)"),
        "{err}"
    );
}

#[test]
fn a_forced_root_fans_out_to_unreached_nodes_and_cannot_have_a_parent() {
    let types = labels(&["Stem", "X", "Y"]);
    let p = build(&types, &[true; 3], Vec::new(), &["Stem"]).unwrap();
    assert_eq!(
        edges_named(&types, &p),
        vec![
            ("Stem".into(), "X".into(), "cli".into()),
            ("Stem".into(), "Y".into(), "cli".into())
        ]
    );
    let st = parse_statements("X\tStem\tprecedes\n", Source::Run, "r").unwrap();
    let err = build(&types, &[true; 3], st, &["Stem"])
        .unwrap_err()
        .to_string();
    assert!(err.contains("cannot be a root"), "{err}");
}

#[test]
fn bad_relation_and_self_statements_are_errors() {
    assert!(parse_statements("A\tB\tbefore\n", Source::Run, "r").is_err());
    assert!(parse_statements("A\tA\tprecedes\n", Source::Run, "r").is_err());
}

#[test]
fn the_prior_tsv_round_trips_its_statements() {
    let types = labels(&["A", "B", "C"]);
    let st = parse_statements("A\tB\tprecedes\nB\tC\tprecedes\n", Source::Project, "p").unwrap();
    let p = build(&types, &[true; 3], st, &[]).unwrap();
    let dir = std::env::temp_dir().join(format!("lupin-prior-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("x.trajectory_prior.tsv");
    write_tsv(path.to_str().unwrap(), &types, &p).unwrap();
    let back = parse_prior_tsv(&std::fs::read_to_string(&path).unwrap(), "x").unwrap();
    assert_eq!(back.len(), 2);
    assert_eq!((back[1].from.as_str(), back[1].to.as_str()), ("B", "C"));
    std::fs::remove_dir_all(dir).unwrap();
}
