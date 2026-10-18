//! The claims a model's text makes, checked against made-up facts.

use super::*;
use crate::annotate::rounds::TypeEvidence;
use crate::describe::consensus::ChoiceInput;
use crate::describe::trajectory::parse_lineages;
use crate::trajectory::edges::{EdgeRow, Verdict};

fn edge(a: &str, b: &str, verdict: Verdict) -> EdgeRow {
    EdgeRow {
        a: a.into(),
        b: b.into(),
        connectivity: 0.5,
        in_prior: true,
        verdict: Some(verdict),
        order_agreement: 0.9,
    }
}

/// stem → wide progenitor → narrow progenitor, which parts toward CT1,
/// CT2 and CT3; CT3 → CT4 mature. Narrow → CT3 is unsupported, CT3 → CT4
/// supported.
fn order() -> OrderFacts {
    use Verdict::{Supported, Unsupported};
    OrderFacts {
        lineages: parse_lineages(
            "L0\t0\tstem_cell > wide_progenitor_cell > narrow_progenitor_cell > CT1_cell\n\
             L1\t0\tstem_cell > wide_progenitor_cell > narrow_progenitor_cell > CT2_cell\n\
             L2\t0\tstem_cell > wide_progenitor_cell > narrow_progenitor_cell > CT3_cell > CT4_mature_cell\n",
        ),
        edges: vec![
            edge("stem_cell", "wide_progenitor_cell", Supported),
            edge("wide_progenitor_cell", "narrow_progenitor_cell", Supported),
            edge("narrow_progenitor_cell", "CT1_cell", Supported),
            edge("narrow_progenitor_cell", "CT2_cell", Supported),
            edge("narrow_progenitor_cell", "CT3_cell", Unsupported),
            edge("CT3_cell", "CT4_mature_cell", Supported),
        ],
        medians: Default::default(),
    }
}

fn ok(text: &str, subject: Option<&str>) {
    let f = order();
    assert_eq!(order_claims(text, &f, subject), Ok(()), "{text}");
}

fn wrong(text: &str, subject: Option<&str>) {
    let f = order();
    assert!(order_claims(text, &f, subject).is_err(), "{text}");
}

#[test]
fn an_edges_verdict_is_the_one_between_the_types_it_names() {
    let ct3 = Some("CT3_cell");
    wrong(
        "CT3 cell sits between the narrow progenitor cell and the CT4 mature cell, with a \
         supported edge from the narrow progenitor cell to CT3 cell and an unsupported edge \
         from CT3 cell to the CT4 mature cell.",
        ct3,
    );
    ok(
        "CT3 cell sits between the narrow progenitor cell and the CT4 mature cell, with an \
         unsupported edge from the narrow progenitor to CT3 cell and a supported edge from \
         CT3 cell to the mature cell.",
        ct3,
    );
    // Named by no type: the subject's own edges.
    ok("CT3 cell has supported and unsupported connections.", ct3);
    wrong("CT1 cell has an unsupported edge.", Some("CT1_cell"));
}

#[test]
fn between_names_a_type_on_either_side() {
    wrong(
        "CT1 cell sits between the stem cell and the wide progenitor cell.",
        Some("CT1_cell"),
    );
    wrong(
        "The CT4 mature cell sits between CT3 cell and the stem cell in the order.",
        Some("CT4_mature_cell"),
    );
    ok(
        "The narrow progenitor cell sits between the wide progenitor cell and CT2 cell.",
        Some("narrow_progenitor_cell"),
    );
}

#[test]
fn an_order_word_runs_from_the_earlier_type_to_the_later() {
    let ct1 = Some("CT1_cell");
    wrong("CT1 cell gives rise to the narrow progenitor cell.", ct1);
    ok(
        "CT1 cell comes after the narrow progenitor cell and ends its lineage.",
        ct1,
    );
    ok(
        "The stem cell, which gives rise to the wide progenitor cell, starts the order.",
        None,
    );
    wrong(
        "The order runs from CT2 cell to the wide progenitor cell.",
        None,
    );
}

#[test]
fn a_role_is_held_by_a_type_the_clause_names_else_the_subject() {
    wrong(
        "The wide progenitor cell is a root.",
        Some("wide_progenitor_cell"),
    );
    ok(
        "The stem cell is at the base of the order.",
        Some("stem_cell"),
    );
    wrong("Lineages branch at CT1 cell.", None);
    ok(
        "The order starts at the stem cell and branches at the narrow progenitor cell into \
         CT1 cell, CT2 cell and CT3 cell.",
        None,
    );
    wrong("Its lineages end at the wide progenitor cell.", None);
    wrong("It is where the lineages branch.", Some("CT2_cell"));
}

/// CT_A tops the enrichment and is SuSiE's second; CT_B is SuSiE's first.
fn cluster() -> ClusterInput {
    let choice = |label: &str, ranks: [Option<usize>; 3]| ChoiceInput {
        evidence: TypeEvidence {
            label: label.into(),
            ..TypeEvidence::default()
        },
        ranks,
        markers: Vec::new(),
    };
    ClusterInput {
        choices: vec![
            choice("CT_A", [Some(0), Some(1), Some(1)]),
            choice("CT_B", [Some(1), Some(0), None]),
            choice("unassigned", [None; 3]),
        ],
        tops: vec![
            (Method::Enrichment, "CT_A".into()),
            (Method::Susie, "CT_B".into()),
            (Method::Score, "CT_C".into()),
        ],
        ..ClusterInput::default()
    }
}

#[test]
fn a_methods_rank_is_the_one_it_gives_the_type_named() {
    let c = cluster();
    let check = |t: &str| rank_claims(t, &c, "CT_A");
    assert!(check("CT A has the highest enrichment score and SuSiE ranks it first").is_err());
    assert_eq!(
        check("CT A has the highest enrichment score and SuSiE ranks it second"),
        Ok(())
    );
    assert_eq!(
        check("CT A ranks first by enrichment, while SuSiE ranks CT B first"),
        Ok(())
    );
    assert!(check("CT A is the top type by SuSiE").is_err());
    assert!(check("CT A, which the score ranks first, fits").is_err());
    assert_eq!(check("CT A is the score's second-highest type"), Ok(()));
    assert_eq!(
        check("CT A fits its markers best"),
        Ok(()),
        "no method named"
    );
}

#[test]
fn a_rank_word_runs_over_a_list_of_methods_but_not_past_another_claim() {
    let c = cluster();
    let check = |t: &str, label: &str| rank_claims(t, &c, label);
    assert!(check("CT A has the highest enrichment and SuSiE ranking", "CT_A").is_err());
    assert!(check(
        "CT A has the highest enrichment, SuSiE ranking, and score",
        "CT_A"
    )
    .is_err());
    assert!(check(
        "CT B has the strongest evidence from enrichment and SuSiE",
        "CT_B"
    )
    .is_err());
    assert!(check("CT A is the top result across all methods", "CT_A").is_err());
    assert_eq!(
        check(
            "CT A has the highest enrichment rank and is supported by SuSiE",
            "CT_A"
        ),
        Ok(()),
        "support is no rank"
    );
    assert_eq!(
        check(
            "CT A has the highest enrichment rank and SuSiE ranks it as likely",
            "CT_A"
        ),
        Ok(())
    );
}
