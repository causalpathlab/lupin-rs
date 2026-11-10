use super::*;
use crate::trajectory::encode_groups;
use crate::trajectory::reference::{ScanpyReference, CONNECTIVITY, DIR};
use legume_numeric::matrix::common_io::read_lines_of_words_delim;

/// scanpy's PAGA on scanpy's reference subset matches, pair by pair.
#[test]
fn connectivity_matches_scanpy_on_the_reference_subset() {
    let path = format!("{DIR}/{CONNECTIVITY}");
    let Some(g) = ScanpyReference::read() else {
        return;
    };
    if !std::path::Path::new(&path).exists() {
        eprintln!("skipping: {path} is absent");
        return;
    }
    let (names, group) = encode_groups(&g.labels);
    let conn = connectivity(&g.neighbours(), &group, names.len()).ratio;
    let code = |n: &str| names.binary_search_by(|x| x.as_ref().cmp(n)).unwrap();

    let want = read_lines_of_words_delim(&path, "\t", 0).unwrap();
    let mut worst = 0.0f64;
    for r in &want.lines {
        let v: f64 = r[2].parse().unwrap();
        worst = worst.max((conn[(code(&r[0]), code(&r[1]))] - v).abs());
    }
    eprintln!(
        "PAGA: {} pairs, worst gap to scanpy {worst:.2e}",
        want.lines.len()
    );
    assert!(worst < 1e-6, "worst connectivity gap to scanpy {worst}");
}

#[test]
fn the_poisson_tail_is_exact_at_small_and_large_means() {
    assert_eq!(poisson_tail(0, 3.0), 1.0);
    // P(X >= 1 | 0.5) = 1 - e^-0.5.
    assert!((poisson_tail(1, 0.5) - (1.0 - (-0.5f64).exp())).abs() < 1e-12);
    // P(X >= 3 | 0.09) = 1 - e^-0.09 (1 + 0.09 + 0.09^2 / 2).
    let want = 1.0 - (-0.09f64).exp() * (1.0 + 0.09 + 0.09 * 0.09 / 2.0);
    assert!((poisson_tail(3, 0.09) - want).abs() < 1e-12);
    // Far in the tail of a large mean, and at its middle.
    assert!(poisson_tail(600, 400.0) < 1e-15);
    assert!((poisson_tail(400, 400.0) - 0.5).abs() < 0.02);
}

#[test]
fn a_pair_is_joined_when_its_edges_are_more_than_noise_gives() {
    // Two types of some 250 cells: 72 edges where chance would give 1450,
    // a ratio of 0.05, are far more than 1% of chance.
    assert!(Connectivity::pair(72.0, 1450.0).joined(0, 1));
    // 18 where chance would give 908: within what 1% noise gives.
    assert!(!Connectivity::pair(18.0, 908.0).joined(0, 1));
    // Two small types: 2 edges of 9.3 expected is a ratio of 0.2 and no
    // evidence; 3 is.
    assert!(!Connectivity::pair(2.0, 9.3).joined(0, 1));
    assert!(Connectivity::pair(3.0, 9.3).joined(0, 1));
    assert!(!Connectivity::pair(0.0, 9.3).joined(0, 1));
    // The ratio of the pairs not joined is no connection at all.
    assert_eq!(Connectivity::pair(2.0, 9.3).joined_ratio()[(0, 1)], 0.0);
    assert!(Connectivity::pair(3.0, 9.3).joined_ratio()[(0, 1)] > 0.3);
}
