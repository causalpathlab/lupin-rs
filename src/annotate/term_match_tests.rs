//! [`super`] with a fake encoder: fixed unit vectors per text.

use super::*;
use crate::annotate::celltype_tree::ClTerms;
use std::collections::HashMap;

const OBO: &str = "[Term]\nid: CL:1\nname: burst forming unit erythroid\nsynonym: \"BFU-E\" RELATED OMO:0003000 []\n\n\
                   [Term]\nid: CL:2\nname: colony forming unit erythroid\n\n\
                   [Term]\nid: CL:3\nname: pancreatic A cell\nsynonym: \"alpha cell\" EXACT []\n\n\
                   [Term]\nid: CL:4\nname: retired cell\nis_obsolete: true\n";

/// Texts as unit vectors on named axes; an unknown text is the zero vector.
struct Fake(HashMap<&'static str, Vec<f32>>);

impl TextEncoder for Fake {
    fn encode(&self, texts: &[&str]) -> anyhow::Result<Vec<Vec<f32>>> {
        Ok(texts
            .iter()
            .map(|t| self.0.get(t).cloned().unwrap_or_else(|| vec![0.0; 3]))
            .collect())
    }
}

fn unit(v: [f32; 3]) -> Vec<f32> {
    let n = v.iter().map(|x| x * x).sum::<f32>().sqrt();
    v.iter().map(|x| x / n).collect()
}

fn fake() -> Fake {
    Fake(HashMap::from([
        ("burst forming unit erythroid", unit([1.0, 0.0, 0.0])),
        ("bfu e", unit([1.0, 0.1, 0.0])),
        ("colony forming unit erythroid", unit([0.9, 0.44, 0.0])),
        ("pancreatic a cell", unit([0.0, 0.0, 1.0])),
        ("alpha cell", unit([0.0, 0.2, 1.0])),
        ("retired cell", unit([0.0, 1.0, 0.0])),
        // labels, as normalised
        ("bfue", unit([1.0, 0.05, 0.0])),
        ("alpha", unit([0.0, 0.15, 1.0])),
        ("erythroid unit", unit([0.95, 0.31, 0.0])),
        ("retired", unit([0.0, 1.0, 0.0])),
    ]))
}

fn index() -> (ClTerms, TermIndex) {
    let terms = ClTerms::parse(OBO, &crate::annotate::cl_rules::shipped());
    let idx = TermIndex::build(&terms, &fake()).unwrap();
    (terms, idx)
}

#[test]
fn a_label_ranks_terms_by_their_closest_counted_text() {
    let (_, idx) = index();
    let v = &fake().encode(&["alpha"]).unwrap()[0];
    let ranked = idx.rank(v, 2);
    assert_eq!(ranked[0].0, "CL:3", "its synonym is nearest");
    assert!(ranked[0].1 > ranked[1].1);
}

#[test]
fn an_obsolete_term_is_never_a_candidate() {
    let (_, idx) = index();
    let v = &fake().encode(&["retired"]).unwrap()[0];
    assert!(idx.rank(v, 4).iter().all(|(id, _)| id != "CL:4"));
}

#[test]
fn only_a_clear_winner_is_accepted() {
    let (terms, idx) = index();
    let accept = Acceptance {
        min_score: 0.9,
        min_gap: 0.05,
    };
    let found = match_labels(
        &terms,
        &idx,
        &fake(),
        &["Alpha", "BFUE", "Erythroid unit", "Unknown"],
        accept,
    )
    .unwrap();
    assert_eq!(found["Alpha"].id, "CL:3");
    assert_eq!(found["BFUE"].id, "CL:1");
    assert!(found["BFUE"].score > found["BFUE"].runner_up);
    assert!(
        !found.contains_key("Erythroid unit"),
        "two terms are about as close"
    );
    assert!(!found.contains_key("Unknown"), "nothing is close");
}

/// An encoder that may not be called with nothing to encode.
struct NoEmpty;

impl TextEncoder for NoEmpty {
    fn encode(&self, texts: &[&str]) -> anyhow::Result<Vec<Vec<f32>>> {
        anyhow::ensure!(!texts.is_empty(), "empty batch");
        Ok(texts.iter().map(|_| vec![1.0, 0.0, 0.0]).collect())
    }
}

#[test]
fn no_labels_need_no_encoding() {
    let (terms, idx) = index();
    let accept = Acceptance {
        min_score: 0.9,
        min_gap: 0.05,
    };
    assert!(match_labels(&terms, &idx, &NoEmpty, &[], accept)
        .unwrap()
        .is_empty());
}
