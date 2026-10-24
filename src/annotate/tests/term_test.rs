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
    let t = test_terms(&x, &terms, &universe, &effect, 200).unwrap();
    assert!(t.z_kt[(0, 0)] > 3.0, "z {}", t.z_kt[(0, 0)]);
    assert!(t.q_kt[(0, 0)] < 0.01, "q {}", t.q_kt[(0, 0)]);
    let others = (1..9).map(|k| t.z_kt[(k, 0)]).fold(f32::MIN, f32::max);
    assert!(
        t.z_kt[(0, 0)] > 2.0 * others,
        "only its own cluster stands out: {} vs {others}",
        t.z_kt[(0, 0)]
    );
    assert!(t.q_kt[(0, 1)] > 0.05, "a flat term is not called");
    let again = test_terms(&x, &terms, &universe, &effect, 200).unwrap();
    assert_eq!(again.p_kt, t.p_kt, "the decoys are seeded");
}

#[test]
fn decoys_keep_a_terms_expression_levels() {
    let x = profile();
    let universe: Vec<usize> = (0..400).collect();
    let bins = AbundanceBins::of(&x, &universe);
    let mut rng = SmallRng::seed_from_u64(2);
    let like = [3usize, 7, 11, 15];
    let drawn = bins.draw_like(&like, &mut rng);
    let bin = |g: usize| bins.of_row[g].unwrap();
    let mut want: Vec<usize> = like.iter().map(|&g| bin(g)).collect();
    let mut got: Vec<usize> = drawn.iter().map(|&g| bin(g)).collect();
    want.sort_unstable();
    got.sort_unstable();
    assert_eq!(got, want);
}
