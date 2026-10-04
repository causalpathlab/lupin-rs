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
[Term]
id: CL:0006
name: other mature cell
is_a: CL:0010 ! broad class
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
        ("Other", "CL:0006"),
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

fn statements(text: &str) -> Vec<Statement> {
    parse_statements(text, Source::Run, "test").unwrap()
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
        "its own chain reaches the stem cell, so the broad class's shortcut does not make it inherited"
    );
    assert_eq!(
        find("Stem", "Other"),
        Some(Source::ClInherited),
        "only the broad class develops from the stem cell"
    );
    assert_eq!(find("Prog", "Stem"), None, "never backwards");
}

#[test]
fn reduction_keeps_only_direct_edges_and_passes_through_dropped_types() {
    let types = labels(&["Stem", "Prog", "Committed", "Pre", "Mature"]);
    let mut map = label_cl();
    map.remove("Other");
    let st = from_ontology(&terms(), &map);
    // "Committed" has too few cells: it is not a node.
    let is_node = [true, true, false, true, true];
    let p = build(&types, &is_node, st).unwrap();
    assert_eq!(
        edges_named(&types, &p),
        vec![
            ("Stem".into(), "Prog".into(), "cl".into()),
            ("Prog".into(), "Pre".into(), "cl".into()),
            ("Pre".into(), "Mature".into(), "cl".into()),
        ]
    );
    assert!(
        p.edges[1].via.is_empty(),
        "the ontology states Prog → Pre outright through the chain"
    );
    assert_eq!(p.roots, vec![vec![0]]);
    assert_eq!(p.component, vec![Some(0), Some(0), None, Some(0), Some(0)]);
    assert!(p.related(0, 4) && p.related(4, 0) && !p.related(2, 2));
    assert_eq!(p.lineages(), vec![(0, vec![0, 1, 3, 4])]);
}

#[test]
fn an_edge_through_a_dropped_type_records_the_way() {
    let types = labels(&["A", "M", "B"]);
    let st = statements("A\tM\tprecedes\nM\tB\tprecedes\n");
    let p = build(&types, &[true, false, true], st).unwrap();
    assert_eq!(
        edges_named(&types, &p),
        vec![("A".into(), "B".into(), "implied".into())]
    );
    assert_eq!(p.edges[0].via, vec![1]);
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
    let p = build(&types, &[true, true], combined).unwrap();
    assert!(p.edges.is_empty());
    assert_eq!(p.component, vec![None, None]);
    assert!(p.roots.is_empty());
}

#[test]
fn a_cycle_is_an_error_naming_its_statements() {
    let st = statements("A\tB\tprecedes\nB\tC\tprecedes\nC\tA\tprecedes\n");
    let err = build(&labels(&["A", "B", "C"]), &[true; 3], st)
        .unwrap_err()
        .to_string();
    assert!(
        err.contains("cycle") && err.contains("C precedes A (run)"),
        "{err}"
    );
}

#[test]
fn forced_roots_fan_out_to_unreached_nodes_but_not_to_each_other() {
    let types = labels(&["Stem", "X", "Y", "Stem2"]);
    let layer = root_layer(&types, &[true; 4], &[], &["Stem", "Stem2"]).unwrap();
    let pairs: Vec<(&str, &str)> = layer
        .iter()
        .map(|s| (s.from.as_str(), s.to.as_str()))
        .collect();
    assert_eq!(
        pairs,
        vec![("Stem", "X"), ("Stem", "Y"), ("Stem2", "X"), ("Stem2", "Y")]
    );
    assert!(layer.iter().all(|s| s.source == Source::Cli));
    let p = build(&types, &[true; 4], combine(&[layer])).unwrap();
    assert_eq!(p.roots, vec![vec![0, 3]]);

    let st = statements("X\tStem\tprecedes\n");
    let err = root_layer(&types, &[true; 4], &st, &["Stem"])
        .unwrap_err()
        .to_string();
    assert!(err.contains("cannot be a root"), "{err}");
    let err = root_layer(&types, &[true; 4], &[], &["Nope"])
        .unwrap_err()
        .to_string();
    assert!(err.contains("not a type of this run"), "{err}");
}

#[test]
fn bad_relation_and_self_statements_are_errors() {
    assert!(parse_statements("A\tB\tbefore\n", Source::Run, "r").is_err());
    assert!(parse_statements("A\tA\tprecedes\n", Source::Run, "r").is_err());
}

#[test]
fn the_prior_tsv_reads_back_as_statements() {
    let types = labels(&["A", "B", "C"]);
    let p = build(
        &types,
        &[true; 3],
        statements("A\tB\tprecedes\tfirst\nB\tC\tprecedes\n"),
    )
    .unwrap();
    let dir = std::env::temp_dir().join(format!("lupin-prior-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("x.trajectory_prior.tsv");
    write_tsv(path.to_str().unwrap(), &types, &p).unwrap();
    let back =
        parse_statements(&std::fs::read_to_string(&path).unwrap(), Source::Run, "x").unwrap();
    assert_eq!(back.len(), 2, "edge rows are skipped");
    assert_eq!(
        (
            back[0].from.as_str(),
            back[0].to.as_str(),
            back[0].note.as_str()
        ),
        ("A", "B", "first")
    );
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn an_unrelated_pair_the_prior_still_orders_is_an_error() {
    let mut st = statements("A\tM\tprecedes\nM\tB\tprecedes\n");
    st.extend(parse_statements("A\tB\tunrelated\n", Source::Project, "p").unwrap());
    let err = build(
        &labels(&["A", "M", "B"]),
        &[true, false, true],
        combine(&[st]),
    )
    .unwrap_err()
    .to_string();
    assert!(
        err.contains("A and B are stated unrelated (project)")
            && err.contains("A → M → B")
            && err.contains("M precedes B (run)"),
        "{err}"
    );
}

#[test]
fn only_a_node_before_it_keeps_a_type_from_being_a_root() {
    let types = labels(&["Small", "R", "X"]);
    let st = statements("Small\tR\tprecedes\n");
    let layer = root_layer(&types, &[false, true, true], &st, &["R"]).unwrap();
    assert_eq!(layer.len(), 1, "R → X");
    let err = root_layer(&types, &[true, true, true], &st, &["R"])
        .unwrap_err()
        .to_string();
    assert!(err.contains("Small → R"), "{err}");
}

#[test]
fn opposite_ontology_statements_are_a_cycle_not_a_choice() {
    let cl = |from: &str, to: &str| Statement {
        from: from.into(),
        to: to.into(),
        relation: Relation::Precedes,
        source: Source::Cl,
        note: String::new(),
    };
    let combined = combine(&[vec![cl("A", "B"), cl("B", "A")]]);
    assert_eq!(combined.len(), 2);
    let err = build(&labels(&["A", "B"]), &[true; 2], combined)
        .unwrap_err()
        .to_string();
    assert!(err.contains("cycle"), "{err}");
    // In a user's file the later line wins, as before.
    let user = parse_statements("A\tB\tprecedes\nB\tA\tprecedes\n", Source::User, "u").unwrap();
    assert_eq!(combine(&[user]).len(), 1);
}

#[test]
fn kind_rows_are_read_only_under_the_kind_header_and_indented_comments_are_skipped() {
    let st = parse_statements(
        "  # an indented comment\nedge\tB\tprecedes\nstatement\tC\tprecedes\n",
        Source::Run,
        "r",
    )
    .unwrap();
    assert_eq!(
        st.iter()
            .map(|s| (s.from.as_str(), s.to.as_str()))
            .collect::<Vec<_>>(),
        vec![("edge", "B"), ("statement", "C")]
    );
    let err = parse_statements("# kind\tfrom\tto\nA\tB\tprecedes\n", Source::Run, "r")
        .unwrap_err()
        .to_string();
    assert!(err.contains("`edge` or `statement`"), "{err}");
}
