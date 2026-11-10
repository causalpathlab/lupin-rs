//! Connectivity between cell groups on the kNN graph, as scanpy 1.10.3's
//! `tl.paga` computes it (connectivity model v1.2).
//!
//! Method: F. A. Wolf, F. K. Hamey, M. Plass, J. Solana, J. S. Dahlin,
//! B. Göttgens, N. Rajewsky, L. Simon, F. J. Theis (2019), "PAGA: graph
//! abstraction reconciles clustering with trajectory inference through a
//! topology preserving map of single cells", *Genome Biology* 20, 59.
//! doi:10.1186/s13059-019-1663-x. Reference implementation replicated here:
//! scanpy 1.10.3.
//!
//! The statistic counts each cell's own directed neighbour list: the edges
//! from group `a` into group `b` plus those from `b` into `a`, against the
//! number expected if every group's edges landed on random cells,
//! `(e_a n_b + e_b n_a) / (n − 1)`, where `e_a` is the number of edges leaving
//! `a`'s cells and `n_a` its size; capped at 1.
//!
//! The ratio is on the scale of random mixing, so it falls as two groups
//! grow (their interface is a shrinking share of their edges) and it is
//! large on a handful of edges between two small groups. A pair therefore
//! also carries its edge count, and is *joined* ([`Connectivity::joined`])
//! only when that count is more than noise gives: a Poisson test against
//! [`NOISE`] of the random expectation, what cells labelled at random at
//! that rate would put between the two groups.

use super::diffusion::Neighbours;
use nalgebra::DMatrix;

/// The share of random mixing that noise may give a pair of groups that do
/// not touch: stray edges and stray labels.
pub(crate) const NOISE: f64 = 0.01;
/// A pair is joined when its edges are this unlikely under [`NOISE`].
pub(crate) const JOINED_P: f64 = 1e-3;

/// The kNN edges between every pair of groups: symmetric `n_groups ×
/// n_groups` matrices, zero on the diagonal.
pub(crate) struct Connectivity {
    /// PAGA's ratio of the edges to their random expectation, capped at 1;
    /// zero between groups no edge joins.
    pub(crate) ratio: DMatrix<f64>,
    /// The edges between the two groups, each cell's own list counted.
    pub(crate) edges: DMatrix<f64>,
    /// The edges expected if every group's edges landed on random cells.
    pub(crate) expected: DMatrix<f64>,
}

impl Connectivity {
    /// Whether groups `a` and `b` are joined: their edges are more than
    /// [`NOISE`] of the random expectation gives, at [`JOINED_P`].
    pub(crate) fn joined(&self, a: usize, b: usize) -> bool {
        poisson_tail(self.edges[(a, b)] as u64, NOISE * self.expected[(a, b)]) < JOINED_P
    }

    /// Two groups with `edges` between them where `expected` would land by
    /// chance.
    #[cfg(test)]
    pub(crate) fn pair(edges: f64, expected: f64) -> Self {
        let m = |v: f64| DMatrix::from_row_slice(2, 2, &[0.0, v, v, 0.0]);
        Self {
            ratio: m((edges / expected).min(1.0)),
            edges: m(edges),
            expected: m(expected),
        }
    }

    /// The ratio of the pairs that are joined, zero for the rest.
    pub(crate) fn joined_ratio(&self) -> DMatrix<f64> {
        DMatrix::from_fn(self.ratio.nrows(), self.ratio.ncols(), |a, b| {
            if self.joined(a, b) {
                self.ratio[(a, b)]
            } else {
                0.0
            }
        })
    }
}

/// `P(X ≥ x)` for `X` Poisson of mean `mean`.
pub(crate) fn poisson_tail(x: u64, mean: f64) -> f64 {
    use statrs::distribution::{DiscreteCDF, Poisson};
    if x == 0 {
        return 1.0;
    }
    Poisson::new(mean).map_or(0.0, |p| p.sf(x - 1))
}

/// [`Connectivity`] of the groups on graph `nb`. `group[c]` is cell `c`'s
/// group; over the graph's representatives, each with its own group.
pub(crate) fn connectivity(nb: &Neighbours, group: &[usize], n_groups: usize) -> Connectivity {
    let n = nb.n_cells();
    let group: Vec<usize> = nb.reps.iter().map(|&c| group[c]).collect();
    let mut size = vec![0.0; n_groups];
    let mut out_edges = vec![0.0; n_groups];
    let mut between = DMatrix::<f64>::zeros(n_groups, n_groups);
    for (c, row) in nb.idx.iter().enumerate() {
        let a = group[c];
        size[a] += 1.0;
        // Always `k − 1` per cell on exact kNN lists; counted as scanpy
        // counts it, from each cell's own list.
        out_edges[a] += row.len() as f64;
        for &j in row {
            let b = group[j];
            if a != b {
                between[(a, b)] += 1.0;
            }
        }
    }
    let edges = &between + between.transpose();
    let expected = DMatrix::from_fn(n_groups, n_groups, |a, b| {
        if a == b {
            0.0
        } else {
            (out_edges[a] * size[b] + out_edges[b] * size[a]) / (n as f64 - 1.0)
        }
    });
    let ratio = DMatrix::from_fn(n_groups, n_groups, |a, b| {
        let (e, x) = (edges[(a, b)], expected[(a, b)]);
        if e == 0.0 {
            0.0
        } else if x == 0.0 {
            1.0
        } else {
            (e / x).min(1.0)
        }
    });
    Connectivity {
        ratio,
        edges,
        expected,
    }
}

#[cfg(test)]
#[path = "type_connectivity_tests.rs"]
mod tests;
