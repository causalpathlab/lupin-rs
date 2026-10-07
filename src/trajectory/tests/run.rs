//! Completing an incomplete prior by connectivity.

use super::*;

/// Symmetric connectivity over `n` types from `(a, b, value)` pairs.
fn conn(n: usize, pairs: &[(usize, usize, f64)]) -> DMatrix<f64> {
    let mut m = DMatrix::zeros(n, n);
    for &(a, b, v) in pairs {
        m[(a, b)] = v;
        m[(b, a)] = v;
    }
    m
}

#[test]
fn an_edgeless_type_joins_along_its_strongest_connection() {
    // 0 → 1 is the prior; 2 has no edge and links to 1 more than to 0.
    let c = conn(3, &[(0, 1, 0.9), (1, 2, 0.4), (0, 2, 0.1)]);
    let (host, edges) = hosts(&[true; 3], &[Some(0), Some(0), None], &c);
    assert_eq!(host, vec![None, None, Some(1)]);
    assert_eq!(edges, vec![(1, 2)]);
}

#[test]
fn a_chain_of_edgeless_types_takes_the_placed_type_it_hangs_from() {
    // 2 hangs from 0; 3 only reaches the prior through 2.
    let c = conn(4, &[(0, 1, 0.9), (0, 2, 0.3), (2, 3, 0.5)]);
    let (host, edges) = hosts(&[true; 4], &[Some(0), Some(0), None, None], &c);
    assert_eq!(host, vec![None, None, Some(0), Some(0)]);
    assert_eq!(edges, vec![(0, 2), (2, 3)]);
}

#[test]
fn two_prior_components_are_never_joined_and_unconnected_types_stay_out() {
    // Components 0 and 1; 2 bridges them, more to 3's side; 4 links nowhere;
    // 5 is not a node type.
    let c = conn(6, &[(0, 2, 0.2), (2, 3, 0.6), (2, 5, 0.9)]);
    let is_node = [true, true, true, true, true, false];
    let component = [Some(0), Some(0), None, Some(1), None, None];
    let (host, edges) = hosts(&is_node, &component, &c);
    assert_eq!(host, vec![None, None, Some(3), None, None, None]);
    assert_eq!(edges, vec![(3, 2)], "only the stronger side, one edge");
}
