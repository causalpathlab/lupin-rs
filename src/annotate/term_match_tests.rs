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
    let idx = TermIndex::build(&TermIndex::named(&terms), &fake(), |_, _| {}).unwrap();
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
fn a_label_gets_its_nearest_terms_by_name_and_score_and_none_is_taken() {
    let (terms, idx) = index();
    let found = suggest(
        &terms,
        &idx,
        &fake(),
        &["Alpha", "BFUE", "Erythroid unit"],
        2,
    )
    .unwrap();
    let alpha = &found["Alpha"];
    assert_eq!(alpha[0].id, "CL:3");
    assert_eq!(alpha[0].name, "pancreatic A cell");
    assert_eq!(found["BFUE"][0].id, "CL:1");
    assert!(found["BFUE"][0].score > found["BFUE"][1].score);
    // Two terms about as close are both offered: the user decides.
    let ids: Vec<&str> = found["Erythroid unit"]
        .iter()
        .map(|s| s.id.as_str())
        .collect();
    assert_eq!(ids.len(), 2);
    assert!(ids.contains(&"CL:1") && ids.contains(&"CL:2"));
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
    assert!(suggest(&terms, &idx, &NoEmpty, &[], 3).unwrap().is_empty());
}

#[test]
fn building_reports_each_batch_up_to_every_name() {
    let terms = ClTerms::parse(OBO, &crate::annotate::cl_rules::shipped());
    let mut seen = Vec::new();
    TermIndex::build(&TermIndex::named(&terms), &fake(), |done, total| {
        seen.push((done, total))
    })
    .unwrap();
    let (done, total) = *seen.last().unwrap();
    assert_eq!(done, total, "the last report is the whole");
    assert!(total >= 5, "every counted name is encoded");
}

#[test]
fn the_index_reads_back_from_its_cache_file() {
    let (_, idx) = index();
    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("sub").join("index.parquet");
    idx.save(&file).unwrap();
    let back = TermIndex::load(&file).unwrap();
    let v = &fake().encode(&["alpha"]).unwrap()[0];
    let (a, b) = (idx.rank(v, 3), back.rank(v, 3));
    assert_eq!(
        a.iter().map(|r| &r.0).collect::<Vec<_>>(),
        b.iter().map(|r| &r.0).collect::<Vec<_>>()
    );
    assert!((a[0].1 - b[0].1).abs() < 1e-6);
}

#[test]
fn the_cache_is_keyed_by_the_ontology_release_and_the_model_revision() {
    let rules = crate::annotate::cl_rules::shipped();
    let v1 = ClTerms::parse(&format!("data-version: r1\n\n{OBO}"), &rules);
    let v2 = ClTerms::parse(&format!("data-version: r2\n\n{OBO}"), &rules);
    let key = |t: &ClTerms, rev: &str| {
        TermIndex::cache_name(t.release.as_deref(), &TermIndex::named(t), "org/model", rev)
    };
    assert_eq!(key(&v1, "abc"), key(&v1, "abc"));
    assert!(key(&v1, "abc").starts_with("r1-"));
    assert_ne!(key(&v1, "abc"), key(&v2, "abc"), "another release");
    assert_ne!(key(&v1, "abc"), key(&v1, "abd"), "another revision");
}
