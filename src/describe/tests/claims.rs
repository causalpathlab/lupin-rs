//! The claims a model's text makes, checked against made-up facts.

use super::*;
use crate::annotate::rounds::TypeEvidence;
use crate::describe::consensus::ChoiceInput;

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
    assert!(
        check(
            "CT A has the top enrichment and all methods ranking",
            "CT_A"
        )
        .is_err(),
        "every method of `all methods` takes the rank, SuSiE's second too"
    );
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
    assert_eq!(
        check(
            "CT A is ranked first by enrichment, and all methods support it",
            "CT_A"
        ),
        Ok(()),
        "a clause saying more than the methods takes no rank from the list before it"
    );
    assert!(check("CT A has the top enrichment, or SuSiE ranking", "CT_A").is_err());
}

#[test]
fn a_listed_method_takes_the_rank_whatever_noun_follows_it() {
    let c = cluster();
    let check = |t: &str| rank_claims(t, &c, "CT_A");
    for tail in ["support", "ranks", "rankings here", "scores"] {
        let t = format!("CT A has the highest enrichment and SuSiE {tail}");
        assert!(check(&t).is_err(), "{t}");
    }
    assert_eq!(
        check("CT A is ranked first by enrichment, and all methods support it"),
        Ok(())
    );
    // Other joins and nouns: the list goes on.
    for t in [
        "CT A has the highest enrichment and SuSiE scores here",
        "CT A has the highest enrichment, and SuSiE support",
        "CT A has the highest enrichment or SuSiE support",
        "CT A has the highest enrichment scores and SuSiE support",
        "CT A has the highest enrichment score and the SuSiE ranking",
        "CT A has the highest enrichment and SuSiE's ranking",
        "CT A has the highest enrichment and SuSiE agrees",
        "CT A has the highest enrichment and all methods agree",
        "CT A has the highest enrichment and SuSiE, which puts it first",
    ] {
        assert!(check(t).is_err(), "{t}");
    }
    // A predicate of its own: an object, a helper, a type, its own rank.
    for t in [
        "CT A has the highest enrichment and all methods support it",
        "CT A has the highest enrichment and SuSiE also supports it",
        "CT A has the highest enrichment and SuSiE is less sure",
        "CT A has the highest enrichment and SuSiE does not",
        "CT A has the highest enrichment and SuSiE favours CT B",
        "CT A has the highest enrichment and SuSiE puts it second",
        "CT A has the highest enrichment and SuSiE ranks second",
        "CT A has the highest enrichment and SuSiE, which puts it second",
    ] {
        assert_eq!(check(t), Ok(()), "{t}");
    }
}
