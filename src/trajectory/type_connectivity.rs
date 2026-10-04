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

use super::diffusion::Neighbours;
use nalgebra::DMatrix;

/// Symmetric `n_groups × n_groups` connectivity, zero on the diagonal and
/// between groups no edge joins. `group[c]` is cell `c`'s group.
pub(crate) fn connectivity(nb: &Neighbours, group: &[usize], n_groups: usize) -> DMatrix<f64> {
    let n = nb.n_cells();
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
    let pair = &between + between.transpose();
    DMatrix::from_fn(n_groups, n_groups, |a, b| {
        let edges = pair[(a, b)];
        if edges == 0.0 {
            return 0.0;
        }
        let expected = (out_edges[a] * size[b] + out_edges[b] * size[a]) / (n as f64 - 1.0);
        if expected == 0.0 {
            1.0
        } else {
            (edges / expected).min(1.0)
        }
    })
}

#[cfg(test)]
#[path = "type_connectivity_tests.rs"]
mod tests;
