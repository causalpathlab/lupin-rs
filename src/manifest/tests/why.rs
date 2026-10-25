//! [`super`]: contrasts asked of a small synthetic round.

use super::*;
use crate::annotate::celltype_tree::{TreeSource, TypeGroup, TypeTree};
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::traits::MatWithNames;
use rand::rngs::SmallRng;
use rand::SeedableRng;

const GENES: usize = 400;
const CLUSTERS: usize = 9;

fn gene(g: usize) -> String {
    format!("G{g}")
}

/// 400 genes over 9 clusters at four expression levels; genes 0..16 are up
/// in K0.
fn expression() -> Expression {
    let mut rng = SmallRng::seed_from_u64(1);
    let mut x = Mat::zeros(GENES, CLUSTERS);
    for g in 0..GENES {
        let level = [5.0, 20.0, 80.0, 300.0][g % 4];
        for k in 0..CLUSTERS {
            let noise: f32 = rand::RngExt::random_range(&mut rng, 0.8..1.25);
            x[(g, k)] = level * noise * if g < 16 && k == 0 { 6.0 } else { 1.0 };
        }
    }
    Expression::new(MatWithNames {
        rows: (0..GENES).map(|g| gene(g).into()).collect(),
        cols: (0..CLUSTERS).map(|k| format!("K{k}").into()).collect(),
        mat: x,
    })
}

fn evidence(label: &str, share: f32, q: f32) -> TypeEvidence {
    TypeEvidence {
        cluster: 0,
        label: label.into(),
        enrichment: Some(share),
        q: Some(q),
        ..TypeEvidence::default()
    }
}

/// Type `a` marks genes 0..8 (up in K0), type `b` genes 6, 7 and 300..304
/// (flat); a term holds gene 3 and genes 8..16, which no panel lists. K0 is
/// labelled `a`, the enrichment's first over `b`.
fn round() -> Round {
    let expr = expression();
    let specific = specific_genes(&expr);
    let panel: Panel = [
        ("a".to_string(), (0..8).map(gene).collect()),
        ("b".to_string(), [6, 7, 300, 301, 302, 303].map(gene).into()),
    ]
    .into();
    let digest = Digest {
        size: 100,
        label: Some("a".into()),
        ..Digest::default()
    };
    let types = vec![evidence("a", 0.8, 0.01), evidence("b", 0.2, 0.04)];
    let ids: Vec<ClusterId> = (0..CLUSTERS as ClusterId).collect();
    let column = |k0: f32| -> Table {
        Table {
            rows: ids.clone(),
            cols: vec!["T:1".into()],
            values: ids.iter().map(|&k| if k == 0 { k0 } else { 0.0 }).collect(),
        }
    };
    let terms = Terms {
        terms: vec![TermGenes {
            id: "T:1".into(),
            name: "a programme".into(),
            genes: std::iter::once(3).chain(8..16).map(gene).collect(),
        }],
        effect: column(1.0),
        z: Some(column(5.0)),
        p: Some(column(1e-6)),
        q: Some(column(1e-4)),
        col: [("T:1".to_string(), 0)].into(),
    };
    Round {
        expr: Some(expr),
        panel,
        tree: None,
        specific,
        evidence: [(0, (digest, types))].into(),
        terms: Some(terms),
        method: Method::Enrichment,
        caveats: Vec::new(),
    }
}

#[test]
fn a_gene_list_is_named_by_its_kind() {
    assert_eq!(GeneList::parse("specific").unwrap(), GeneList::Specific);
    assert_eq!(
        GeneList::parse("markers:NK").unwrap(),
        GeneList::Markers("NK".into())
    );
    assert_eq!(
        GeneList::parse("term:GO:0002228").unwrap(),
        GeneList::Term("GO:0002228".into())
    );
    assert!(GeneList::parse("markers:").is_err());
    assert!(GeneList::parse("top").is_err());
}

#[test]
fn readings_favour_a_side_only_when_none_leans_the_other_way() {
    use Favours::*;
    assert_eq!(Favours::of([Fact, Neither, CannotTell]), Fact);
    assert_eq!(Favours::of([Fact, Foil]), Split);
    assert_eq!(Favours::of([Neither, CannotTell]), Neither);
    assert_eq!(Favours::of([CannotTell]), CannotTell);
    assert_eq!(Favours::of([]), CannotTell);
}

#[test]
fn a_gene_shown_is_placed_and_one_not_shown_says_why() {
    let r = round();
    let a = why_gene(&r, 0, "G3", Some("G300"), &GeneList::Specific).unwrap();
    assert_eq!(a.favours, Favours::Fact);
    assert!(a.fact["rank"].as_u64().is_some(), "{}", a.fact);
    assert_eq!(a.foil["shown"], json!(false));
    let why = a.foil["why_not"].as_str().unwrap();
    assert!(why.contains("at or below K0's mean"), "{why}");
    assert_eq!(a.foil["marker_of"], json!(["b"]));

    // Without a foil, the gene listed next.
    let next = why_gene(&r, 0, "G3", None, &GeneList::Specific).unwrap();
    assert!(next.question.unwrap().foil.is_some());
    assert!(next.notes[0].contains("the gene listed next"));

    let m = why_gene(&r, 0, "G8", None, &GeneList::Markers("a".into())).unwrap();
    assert_eq!(m.fact["why_not"], json!("not a marker of a"));
}

#[test]
fn a_label_weighs_the_panel_as_one_line_and_terms_without_its_genes() {
    let r = round();
    let a = why_label(&r, 0, None, None, 2).unwrap();
    let q = a.question.as_ref().unwrap();
    assert_eq!((q.fact.as_str(), q.foil.as_deref()), ("a", Some("b")));
    let v = a.verdict.as_ref().unwrap();
    assert_eq!((v.favours, v.fact, v.foil), (Favours::Fact, 2, 0), "{a:#?}");

    let panel = &a.children[0];
    assert_eq!(panel.line.as_deref(), Some("panel"));
    assert!(panel.fact["markers"]["p"].as_f64().unwrap() < 0.05);
    assert!(
        panel.notes.iter().any(|n| n.contains("2 shared markers")),
        "{:?}",
        panel.notes
    );
    // Its gene question: a's top distinct marker against b's.
    assert_eq!(
        panel.children[0].question.as_ref().unwrap().kind,
        Kind::Gene
    );

    let terms = &a.children[1];
    assert_eq!(terms.favours, Favours::Fact);
    let t = &terms.fact["terms"][0];
    assert_eq!(t["genes_without_panel"], json!(8));
    assert_eq!(t["independent"], json!(true), "{t}");
    assert_eq!(terms.foil["terms"], json!([]));

    let text = render(&a);
    assert!(text.starts_with("why a rather than b, given K0?"), "{text}");
    assert!(
        text.contains("verdict: favours the fact (2 lines for a, 0 for b"),
        "{text}"
    );
    assert!(text.contains("  panel line"), "{text}");
}

#[test]
fn a_term_backs_no_side_once_the_panel_genes_are_all_it_has() {
    let mut r = round();
    // The term now holds a's markers only.
    r.terms.as_mut().unwrap().terms[0].genes = (0..6).map(gene).collect();
    let a = why_label(&r, 0, Some("a"), Some("b"), 1).unwrap();
    let terms = &a.children[1];
    assert_eq!(terms.favours, Favours::Neither, "{:#?}", terms.fact);
    assert_eq!(terms.fact["terms"][0]["independent"], json!(false));
    assert!(a.children.iter().all(|c| c.children.is_empty()), "depth 1");
}

#[test]
fn a_coarse_label_with_no_rival_outside_asks_between_its_types() {
    let mut r = round();
    r.tree = Some(TypeTree {
        source: TreeSource::CellOntology,
        release: None,
        label_cl: Default::default(),
        groups: vec![TypeGroup {
            name: "grp".into(),
            cl_id: None,
            members: vec!["a".into(), "b".into()],
        }],
    });
    let a = why_label(&r, 0, Some("grp"), None, 1).unwrap();
    let q = a.question.as_ref().unwrap();
    assert_eq!((q.fact.as_str(), q.foil.as_deref()), ("a", Some("b")));
    assert!(
        a.notes[0].contains("between its types a and b"),
        "{:?}",
        a.notes
    );
}

#[test]
fn a_term_ranks_by_its_effect() {
    let r = round();
    let a = why_term(&r, 0, "a programme", None);
    assert!(
        a.unwrap_err().to_string().contains("no other term"),
        "one term has no rival"
    );
    let err = why_term(&r, 0, "programme", None).unwrap_err().to_string();
    assert!(err.contains("a programme"), "names the near term: {err}");
}
