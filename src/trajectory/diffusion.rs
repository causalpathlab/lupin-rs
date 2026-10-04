//! Diffusion map and diffusion pseudotime, as scanpy 1.10.3 computes them
//! (`neighbors(method='gauss', knn=True)`, `diffmap`, `dpt`).
//!
//! Methods:
//! - Diffusion pseudotime: L. Haghverdi, M. Büttner, F. A. Wolf, F. Buettner,
//!   F. J. Theis (2016), "Diffusion pseudotime robustly reconstructs lineage
//!   branching", *Nature Methods* 13, 845–848. doi:10.1038/nmeth.3971
//! - Diffusion maps and the density normalisation: R. R. Coifman, S. Lafon
//!   (2006), "Diffusion maps", *Applied and Computational Harmonic Analysis*
//!   21(1), 5–30. doi:10.1016/j.acha.2006.04.006
//! - The Gaussian kernel with per-cell bandwidths: L. Haghverdi, F. Buettner,
//!   F. J. Theis (2015), "Diffusion maps for high-dimensional single-cell
//!   analysis of differentiation data", *Bioinformatics* 31(18), 2989–2998.
//!   doi:10.1093/bioinformatics/btv325
//! - Reference implementation replicated here: scanpy 1.10.3; F. A. Wolf,
//!   P. Angerer, F. J. Theis (2018), "SCANPY: large-scale single-cell gene
//!   expression data analysis", *Genome Biology* 19, 15.
//!   doi:10.1186/s13059-017-1382-0
//! - The eigenpairs come from a randomised SVD (legume-numeric `rsvd_with`):
//!   N. Halko, P.-G. Martinsson, J. A. Tropp (2011), "Finding structure with
//!   randomness", *SIAM Review* 53(2), 217–288. doi:10.1137/090771806
//!
//! The kNN lists, sparse matrices, graph components, median and randomised
//! SVD all come from legume-numeric; what is here is scanpy's own arithmetic
//! on top of them: the Gaussian kernel with per-cell bandwidths, the density
//! and symmetric normalisation, signed eigenpairs, and the DPT distance.

use anyhow::{bail, ensure, Result};
use legume_numeric::matrix::graph::{connected_components, AdjListGraph};
use legume_numeric::matrix::knn::{knn_rows, ALL_PAIRS_THRESHOLD};
use legume_numeric::matrix::traits::{MatTriplets, RandomizedAlgs, RsvdArgs};
use legume_numeric::matrix::utils::median;
use log::warn;
use nalgebra::DMatrix;
use nalgebra_sparse::CscMatrix;
use rustc_hash::FxHashMap;

/// Diffusion components computed, as `tl.diffmap(n_comps=15)`; DPT uses the
/// first `n_dcs` of them, as `tl.dpt(n_dcs=10)` does by default.
pub(crate) const DIFFMAP_COMPS: usize = 15;

/// Eigenvalues at or above this enter the DPT distance unweighted (scanpy
/// compares in float32).
const UNWEIGHTED_FROM: f32 = 0.9994;

/// The randomised SVD's starting effort: the leading eigenvalues of a
/// diffusion operator sit close together, and the defaults (5, 5) do not
/// separate them.
const EIGEN_ARGS: RsvdArgs = RsvdArgs {
    power_iters: 20,
    oversample: 10,
};

/// Doubling the power iterations stops here.
const MAX_POWER_ITERS: usize = 160;

/// Largest residual `‖Su − λu‖` accepted for a unit eigenvector.
const RESIDUAL_TOL: f64 = 1e-3;

/// A root type with more cells than this has its medoid searched among the
/// cells nearest its centroid, so the search stays quadratic in this many.
const MEDOID_POOL: usize = 2000;

/// Each representative cell's `k − 1` nearest other representatives and
/// their distances: the one kNN graph the diffusion map and PAGA share. `k`
/// counts the cell itself, as scanpy's `n_neighbors` does. Cells that
/// coincide exactly are collapsed onto one representative, the first of them.
pub(crate) struct Neighbours {
    pub(crate) idx: Vec<Vec<usize>>,
    dist: Vec<Vec<f32>>,
    /// The cell each representative is.
    pub(crate) reps: Vec<usize>,
    /// Each cell's representative, as an index into `reps`.
    pub(crate) rep_of: Vec<usize>,
}

impl Neighbours {
    /// The kNN lists of `geometry` (cells × dims, already prepared), over one
    /// representative per set of exactly coinciding cells (with a warning):
    /// scanpy's kernel is not defined at a zero distance.
    pub(crate) fn new(geometry: &DMatrix<f32>, k: usize) -> Result<Self> {
        let n = geometry.nrows();
        ensure!(k >= 2, "k counts the cell itself, so it must be at least 2");
        let bad = (0..n)
            .filter(|&i| geometry.row(i).iter().any(|v| !v.is_finite()))
            .count();
        ensure!(
            bad == 0,
            "{bad} cells have a non-finite coordinate in the geometry"
        );
        let (reps, rep_of) = collapse_duplicates(geometry);
        if reps.len() < n {
            warn!(
                "{} of {n} cells coincide exactly with another cell and share its \
                 neighbours, pseudotime and diffusion coordinates ({} distinct cells)",
                n - reps.len(),
                reps.len()
            );
        }
        let m = reps.len();
        ensure!(m > k, "{m} distinct cells are too few for {k} neighbours");
        if m > ALL_PAIRS_THRESHOLD {
            warn!(
                "{m} cells: neighbours come from the approximate inverted-file search \
                 (exact up to {ALL_PAIRS_THRESHOLD}), so the result no longer replicates scanpy"
            );
        }
        let distinct = if m < n {
            DMatrix::from_fn(m, geometry.ncols(), |i, j| geometry[(reps[i], j)])
        } else {
            geometry.clone()
        };
        let (idx, dist) = knn_rows(&distinct, k - 1);
        Ok(Self {
            idx,
            dist,
            reps,
            rep_of,
        })
    }

    /// Representatives in the graph.
    pub(crate) fn n_cells(&self) -> usize {
        self.idx.len()
    }

    /// Connected component of each cell in the (symmetrised) kNN graph.
    fn components(&self) -> Vec<usize> {
        let edges: Vec<(usize, usize)> = self
            .idx
            .iter()
            .enumerate()
            .flat_map(|(i, row)| row.iter().map(move |&j| (i, j)))
            .collect();
        connected_components(&AdjListGraph::from_unweighted_edges(self.n_cells(), &edges))
    }
}

/// The first cell of each set of cells with identical coordinates, and each
/// cell's set as an index into that list.
fn collapse_duplicates(geometry: &DMatrix<f32>) -> (Vec<usize>, Vec<usize>) {
    let mut first: FxHashMap<Vec<u32>, usize> = FxHashMap::default();
    let mut reps = Vec::new();
    let rep_of = (0..geometry.nrows())
        .map(|i| {
            // `+ 0.0` makes −0 and +0 one key.
            let key: Vec<u32> = geometry
                .row(i)
                .iter()
                .map(|v| (v + 0.0).to_bits())
                .collect();
            *first.entry(key).or_insert_with(|| {
                reps.push(i);
                reps.len() - 1
            })
        })
        .collect();
    (reps, rep_of)
}

/// The diffusion components of a cell kNN graph, per cell (copies of a
/// representative take its rows).
pub(crate) struct DiffusionMap {
    /// Eigenvalues, largest first, rounded to float32 as scanpy keeps them.
    pub(crate) evals: Vec<f32>,
    /// Eigenvectors of the symmetric transition matrix, one column each.
    pub(crate) evecs: DMatrix<f64>,
    /// `evecs` with each column weighted as in the DPT distance: `λ / (1 − λ)`
    /// below [`UNWEIGHTED_FROM`], 1 at or above it. DPT distance is the
    /// Euclidean distance between rows of this.
    coords: DMatrix<f64>,
    /// Connected component of each cell in the kNN graph.
    component: Vec<usize>,
}

impl DiffusionMap {
    /// The diffusion map on `nb` with `n_comps` components, of which the
    /// DPT distance uses the first `n_dcs` (all of them if `n_dcs` is larger).
    pub(crate) fn new(nb: &Neighbours, n_comps: usize, n_dcs: usize) -> Result<Self> {
        let n = nb.n_cells();
        ensure!(n_dcs >= 1, "at least one diffusion component is needed");
        let kernel = gauss_kernel(&nb.idx, &nb.dist)?;
        let s = symmetric_transitions(n, kernel)?;
        let (evals, evecs) = transition_eigenpairs(&s, n_comps.max(n_dcs).min(n - 1))?;
        // Every cell gets its representative's row.
        let evecs = DMatrix::from_fn(nb.rep_of.len(), evecs.ncols(), |i, j| {
            evecs[(nb.rep_of[i], j)]
        });
        let n_dcs = n_dcs.min(evals.len());
        for &l in &evals[..n_dcs] {
            if (l - UNWEIGHTED_FROM).abs() < 1e-4 {
                warn!(
                    "eigenvalue {l} sits within 1e-4 of {UNWEIGHTED_FROM}, where scanpy switches \
                     it between weighted and unweighted; pseudotime may differ from scanpy's"
                );
            }
        }
        let mut coords = evecs.columns(0, n_dcs).into_owned();
        for (mut col, &l) in coords.column_iter_mut().zip(&evals) {
            if l < UNWEIGHTED_FROM {
                col *= f64::from(l / (1.0 - l));
            }
        }
        Ok(Self {
            evals,
            evecs,
            coords,
            component: {
                let c = nb.components();
                nb.rep_of.iter().map(|&r| c[r]).collect()
            },
        })
    }

    /// DPT distance between cells `i` and `c`; infinite across components.
    pub(crate) fn distance(&self, i: usize, c: usize) -> f64 {
        if self.component[i] != self.component[c] {
            return f64::INFINITY;
        }
        self.coords
            .row(i)
            .iter()
            .zip(self.coords.row(c).iter())
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt()
    }

    /// The medoid of `cells` in DPT distance, among those in the graph
    /// component holding most of them. Beyond [`MEDOID_POOL`] cells the exact
    /// search runs over that many cells nearest the type's centroid.
    pub(crate) fn medoid(&self, cells: &[usize]) -> Option<usize> {
        let mut count: FxHashMap<usize, usize> = FxHashMap::default();
        for &c in cells {
            *count.entry(self.component[c]).or_default() += 1;
        }
        let (&main, _) = count
            .iter()
            .max_by_key(|(comp, n)| (**n, std::cmp::Reverse(**comp)))?;
        let mut pool: Vec<usize> = cells
            .iter()
            .copied()
            .filter(|&c| self.component[c] == main)
            .collect();
        if pool.len() > MEDOID_POOL {
            let d = self.coords.ncols();
            let mut centroid = vec![0.0; d];
            for &c in &pool {
                for (s, v) in centroid.iter_mut().zip(self.coords.row(c).iter()) {
                    *s += v;
                }
            }
            centroid.iter_mut().for_each(|s| *s /= pool.len() as f64);
            let to_centroid = |c: usize| -> f64 {
                self.coords
                    .row(c)
                    .iter()
                    .zip(&centroid)
                    .map(|(a, b)| (a - b) * (a - b))
                    .sum()
            };
            let mut keyed: Vec<(f64, usize)> = pool.iter().map(|&c| (to_centroid(c), c)).collect();
            keyed.sort_by(|a, b| a.0.total_cmp(&b.0));
            pool = keyed
                .into_iter()
                .take(MEDOID_POOL)
                .map(|(_, c)| c)
                .collect();
        }
        pool.iter()
            .map(|&c| (c, pool.iter().map(|&o| self.distance(c, o)).sum::<f64>()))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(c, _)| c)
    }

    /// Each cell's DPT distance to the nearest of `roots`.
    pub(crate) fn distances_from(&self, roots: &[usize]) -> Vec<f64> {
        (0..self.coords.nrows())
            .map(|c| {
                roots
                    .iter()
                    .map(|&r| self.distance(r, c))
                    .fold(f64::INFINITY, f64::min)
            })
            .collect()
    }

    /// Pseudotime from `root`: DPT distance divided by its largest finite
    /// value, as scanpy scales it.
    #[cfg(test)]
    pub(crate) fn pseudotime(&self, root: usize) -> Vec<f32> {
        let d = self.distances_from(&[root]);
        let top = d
            .iter()
            .copied()
            .filter(|x| x.is_finite())
            .fold(0.0, f64::max);
        d.into_iter()
            .map(|x| (if top > 0.0 { x / top } else { x }) as f32)
            .collect()
    }
}

/// scanpy's Gaussian kernel over each cell's `k − 1` nearest other cells,
/// as `(i, j, w)` sorted by `(i, j)` with both directions present: a pair
/// only one cell lists takes the value that cell gives it.
fn gauss_kernel(idx: &[Vec<usize>], dist: &[Vec<f32>]) -> Result<Vec<(usize, usize, f64)>> {
    let sigma_sq: Vec<f64> = dist
        .iter()
        .map(|d| f64::from(median(&d.iter().map(|x| x * x).collect::<Vec<f32>>())))
        .collect();
    let zero = sigma_sq.iter().filter(|&&s| s <= 0.0).count();
    if zero > 0 {
        bail!(
            "{zero} cells have a zero kernel bandwidth (most of their neighbours coincide \
             with them); scanpy's kernel is undefined there, so use another embedding"
        );
    }
    let sigma: Vec<f64> = sigma_sq.iter().map(|s| s.sqrt()).collect();

    // Each listed pair goes in both ways, flagged 0 where the cell listed it
    // and 1 where it is the mirror; sorting by (i, j, flag) and keeping the
    // first of each (i, j) keeps a cell's own value wherever it has one.
    let mut w: Vec<(usize, usize, u8, f64)> = Vec::with_capacity(2 * idx.len() * idx[0].len());
    for (i, (row, d)) in idx.iter().zip(dist).enumerate() {
        for (&j, &dij) in row.iter().zip(d) {
            let den = sigma_sq[i] + sigma_sq[j];
            let v = (2.0 * sigma[i] * sigma[j] / den).sqrt() * (-f64::from(dij * dij) / den).exp();
            w.push((i, j, 0, v));
            w.push((j, i, 1, v));
        }
    }
    w.sort_unstable_by_key(|&(i, j, f, _)| (i, j, f));
    w.dedup_by_key(|&mut (i, j, _, _)| (i, j));
    Ok(w.into_iter().map(|(i, j, _, v)| (i, j, v)).collect())
}

/// `S = Z⁻¹ K Z⁻¹` with `K = Q⁻¹ W Q⁻¹`, `q` the kernel's degrees and `z`
/// the square roots of `K`'s, as scanpy's `compute_transitions` builds it.
fn symmetric_transitions(n: usize, mut kernel: Vec<(usize, usize, f64)>) -> Result<CscMatrix<f64>> {
    let mut q = vec![0.0; n];
    for &(i, _, v) in &kernel {
        q[i] += v;
    }
    let mut z = vec![0.0; n];
    for e in &mut kernel {
        e.2 /= q[e.0] * q[e.1];
        z[e.0] += e.2;
    }
    let z: Vec<f64> = z.into_iter().map(f64::sqrt).collect();
    for e in &mut kernel {
        e.2 /= z[e.0] * z[e.1];
    }
    CscMatrix::from_nonzero_triplets(n, n, &kernel)
}

/// The `n_dcs` eigenpairs of the symmetric `s` largest in magnitude, sorted
/// by signed eigenvalue, largest first. The randomised SVD gives the
/// vectors; each eigenvalue is its vector's Rayleigh quotient `uᵀSu`, which
/// carries the sign. Every pair must pass the residual check; until they do
/// the power iterations double, up to [`MAX_POWER_ITERS`].
fn transition_eigenpairs(s: &CscMatrix<f64>, n_dcs: usize) -> Result<(Vec<f32>, DMatrix<f64>)> {
    let mut args = EIGEN_ARGS;
    loop {
        let (u, _, _) = s.rsvd_with(n_dcs, &args)?;
        let su = s * &u;
        let mut pairs: Vec<(f64, usize)> = Vec::with_capacity(u.ncols());
        let mut worst: f64 = 0.0;
        for c in 0..u.ncols() {
            let lambda = u.column(c).dot(&su.column(c));
            worst = worst.max((su.column(c) - u.column(c) * lambda).norm());
            pairs.push((lambda, c));
        }
        // A NaN residual must fail, not pass.
        if worst.is_nan() || worst > RESIDUAL_TOL {
            ensure!(
                worst.is_finite(),
                "the transition matrix has non-finite entries (residual {worst})"
            );
        } else {
            pairs.sort_by(|a, b| b.0.total_cmp(&a.0));
            let evals = pairs.iter().map(|&(l, _)| l as f32).collect();
            let evecs =
                DMatrix::from_columns(&pairs.iter().map(|&(_, c)| u.column(c)).collect::<Vec<_>>());
            return Ok((evals, evecs));
        }
        // Another try: say so on the progress popup, halfway through the
        // stage at most.
        let tries = (MAX_POWER_ITERS / EIGEN_ARGS.power_iters).ilog2() as usize + 1;
        let tried = (args.power_iters / EIGEN_ARGS.power_iters).ilog2() as usize + 1;
        super::run::STAGES.within(
            super::run::STAGE_DIFFUSION,
            tried,
            2 * tries,
            Some(&format!("{} power iterations", 2 * args.power_iters)),
        );
        ensure!(
            args.power_iters < MAX_POWER_ITERS,
            "diffusion components did not converge (residual {worst:.2e} after {} power \
             iterations); the spectrum may hold a ± pair of equal size",
            args.power_iters
        );
        args.power_iters *= 2;
        args.oversample += EIGEN_ARGS.oversample;
    }
}

#[cfg(test)]
#[path = "diffusion_tests.rs"]
mod tests;
