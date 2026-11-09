//! A trajectory's order in words, written by code, on made-up facts.

use super::*;
use crate::trajectory::edges::Verdict;
use std::path::Path;

fn edge(a: &str, b: &str, verdict: Option<Verdict>) -> EdgeRow {
    EdgeRow {
        a: a.into(),
        b: b.into(),
        connectivity: 0.8,
        in_prior: matches!(verdict, Some(Verdict::Supported | Verdict::Unsupported)),
        verdict,
        order_agreement: 0.9,
    }
}

fn lineage(name: &str, component: &str, path: &[&str]) -> lineages::Lineage {
    lineages::Lineage {
        name: name.into(),
        component: component.into(),
        path: path.iter().map(|t| t.to_string()).collect(),
    }
}

/// TY3 → TY5, which parts toward TY2 and toward TY4; TY1 is a
/// component of its own; TY6 pairs with TY4 as a candidate only.
fn facts() -> OrderFacts {
    OrderFacts {
        lineages: vec![
            lineage("L0", "0", &["TY3", "TY5", "TY2"]),
            lineage("L1", "0", &["TY3", "TY5", "TY4"]),
            lineage("L2", "1", &["TY1"]),
        ],
        edges: vec![
            edge("TY3", "TY5", Some(Verdict::Supported)),
            edge("TY5", "TY2", Some(Verdict::Supported)),
            edge("TY5", "TY4", Some(Verdict::Unsupported)),
            edge("TY4", "TY6", Some(Verdict::Candidate)),
            edge("TY2", "TY1", None),
        ],
    }
}

#[test]
fn the_facts_give_roots_branch_points_and_ends() {
    let f = facts();
    assert_eq!(f.roots(), ["TY3", "TY1"]);
    assert_eq!(f.branch_points(), ["TY5"], "where L0 and L1 part");
    assert_eq!(f.ends(), ["TY2", "TY4", "TY1"]);
    assert_eq!(f.components(), 2);
    assert_eq!(
        f.types(),
        ["TY3", "TY5", "TY2", "TY4", "TY1", "TY6"],
        "the lineages' types, then those only an edge names"
    );
}

#[test]
fn the_texts_say_the_order_and_each_types_place() {
    let f = facts();
    let run = run_text(&f);
    assert!(run.contains("The order starts at TY3 and TY1."), "{run}");
    assert!(run.contains("It branches after TY5."), "{run}");
    assert!(run.contains("end at TY2, TY4 and TY1."), "{run}");
    assert!(run.contains("L0: TY3 > TY5 > TY2"), "{run}");
    let ty5 = type_text(&f, "TY5");
    assert!(
        ty5.starts_with("TY5 is on L0 and L1, after TY3 and before TY2 and TY4."),
        "{ty5}"
    );
    assert!(ty5.contains("TY5 → TY4 unsupported"), "{ty5}");
    let ty6 = type_text(&f, "TY6");
    assert!(ty6.starts_with("TY6 is on no lineage."), "{ty6}");
    assert!(ty6.contains("TY4 → TY6 candidate"), "{ty6}");
}

#[test]
fn the_texts_are_the_runs_paragraph_then_each_types_line_and_read_as_markdown() {
    let f = facts();
    let t = texts(&f);
    assert_eq!(t.len(), 1 + f.types().len());
    assert_eq!(
        t[0],
        OrderText {
            subject: None,
            text: run_text(&f)
        }
    );
    let ty5 = t
        .iter()
        .find(|t| t.subject.as_deref() == Some("TY5"))
        .unwrap();
    assert_eq!(ty5.text, type_text(&f, "TY5"));
    let md = markdown("t.senna.json", &t);
    assert!(md.starts_with("# Trajectory description\n\nFrom `t.senna.json`.\n\n"));
    assert!(md.contains(&format!("## TY5\n\n{}\n", ty5.text)), "{md}");
}

/// A trajectory manifest in `dir` with the facts' lineages and edges and a
/// pseudotime per cell.
pub(crate) fn fixture_trajectory(dir: &Path) -> PathBuf {
    use crate::manifest::run::RunManifest;
    use legume_numeric::matrix::parquet::{write_named_table, Column};
    let f = facts();
    let at = |name: &str| dir.join(name).to_string_lossy().into_owned();
    lineages::write(&at("t.trajectory_lineages.tsv"), &f.lineages).unwrap();
    edges::write(&at("t.trajectory_edges.parquet"), &f.edges).unwrap();
    let cells: Vec<Box<str>> = (0..6).map(|i| format!("c{i}").into()).collect();
    let types: Vec<Box<str>> = ["TY3", "TY3", "TY5", "TY2", "TY2", "TY4"]
        .iter()
        .map(|s| Box::from(*s))
        .collect();
    let pt = [0.0f32, 0.02, 0.2, 0.6, 0.8, 0.6];
    write_named_table(
        &at("t.cell_pseudotime.parquet"),
        "cell",
        &cells,
        &[
            ("pseudotime".into(), Column::F32(&pt)),
            ("type".into(), Column::Str(&types)),
        ],
    )
    .unwrap();
    let mut m = RunManifest::new(crate::manifest::run::RunKind::Topic, "t");
    m.trajectory.lineages = Some("t.trajectory_lineages.tsv".into());
    m.trajectory.edges = Some("t.trajectory_edges.parquet".into());
    m.trajectory.pseudotime = Some("t.cell_pseudotime.parquet".into());
    let file = dir.join("t.senna.json");
    m.save(&file).unwrap();
    file
}

#[test]
fn the_facts_read_from_a_trajectory_run_with_the_files_they_come_from() {
    let dir = tempfile::tempdir().unwrap();
    let file = fixture_trajectory(dir.path());
    let loaded = crate::manifest::run::load(&file.to_string_lossy()).unwrap();
    let f = OrderFacts::read(&loaded).unwrap();
    assert_eq!(f.lineages, facts().lineages);
    assert_eq!(f.edges.len(), 5);
    let [lineages, edges] = OrderFacts::files(&loaded);
    assert!(lineages.unwrap().ends_with("t.trajectory_lineages.tsv"));
    assert!(edges.unwrap().ends_with("t.trajectory_edges.parquet"));
}
