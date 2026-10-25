//! [`super::test_terms`]: a term's effect against its decoys.

use super::*;

/// 400 genes over 9 clusters at four expression levels; genes 0..8 are up
/// in cluster 0, genes 8..16 vary at random.
fn profile() -> Mat {
    let mut rng = SmallRng::seed_from_u64(1);
    let mut x = Mat::zeros(400, 9);
    for g in 0..400 {
        let level = [5.0, 20.0, 80.0, 300.0][g % 4];
        for k in 0..9 {
            let noise: f32 = rand::RngExt::random_range(&mut rng, 0.8..1.25);
            x[(g, k)] = level * noise * if g < 8 && k == 0 { 6.0 } else { 1.0 };
        }
    }
    x
}

#[test]
fn a_terms_test_agrees_with_its_effect() {
    let x = profile();
    let terms: Vec<(Box<str>, Vec<usize>)> = vec![
        ("up".into(), (0..8).collect()),
        ("flat".into(), (8..16).collect()),
    ];
    let universe: Vec<usize> = (0..400).collect();
    let effect = enrichment::ontology_module_score(&x, &terms, &universe)
        .unwrap()
        .effect_kt;
    let t = test_terms(&x, &terms, &universe, &effect).unwrap();
    assert!(t.z_kt[(0, 0)] > 3.0, "z {}", t.z_kt[(0, 0)]);
    assert!(t.q_kt[(0, 0)] < 0.01, "q {}", t.q_kt[(0, 0)]);
    let others = (1..9).map(|k| t.z_kt[(k, 0)]).fold(f32::MIN, f32::max);
    assert!(
        t.z_kt[(0, 0)] > 2.0 * others,
        "only its own cluster stands out: {} vs {others}",
        t.z_kt[(0, 0)]
    );
    assert!(t.q_kt[(0, 1)] > 0.05, "a flat term is not called");
    let again = test_terms(&x, &terms, &universe, &effect).unwrap();
    assert_eq!(again.p_kt, t.p_kt, "the decoys are seeded");
}

#[test]
fn the_decoys_effect_is_the_module_scores() {
    let x = profile();
    let universe: Vec<usize> = (0..400).collect();
    let terms: Vec<(Box<str>, Vec<usize>)> = vec![
        ("a".into(), (0..8).collect()),
        ("b".into(), (100..160).step_by(3).collect()),
    ];
    let want = enrichment::ontology_module_score(&x, &terms, &universe)
        .unwrap()
        .effect_kt;
    let bg = Background::of(&x, &universe);
    for (t, (_, rows)) in terms.iter().enumerate() {
        let members: Vec<usize> = bg.members(rows).iter().map(|&(i, _)| i as usize).collect();
        let got = bg.effect(members.iter().copied());
        for (k, g) in got.iter().enumerate() {
            assert!(
                (g - want[(k, t)]).abs() < 1e-5,
                "{t} {k}: {g} vs {}",
                want[(k, t)]
            );
        }
    }
}

#[test]
fn a_term_its_strata_cannot_hide_is_not_tested() {
    let x = profile();
    // Ten annotated genes in one stratum, nine of them the term's.
    let universe: Vec<usize> = (0..10).collect();
    let terms: Vec<(Box<str>, Vec<usize>)> = vec![("most".into(), (0..10).collect())];
    let effect = enrichment::ontology_module_score(&x, &terms, &universe)
        .unwrap()
        .effect_kt;
    let t = test_terms(&x, &terms, &universe, &effect).unwrap();
    assert!(
        (0..9).all(|k| t.z_kt[(k, 0)].is_nan() && t.q_kt[(k, 0)].is_nan()),
        "z {:?}",
        t.z_kt.column(0).as_slice()
    );
}
