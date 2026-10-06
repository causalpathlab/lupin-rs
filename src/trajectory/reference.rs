//! scanpy's results on a reference subset of the bench data, the external
//! reference the tests check lupin against. The fixtures are local files
//! (`src/tests/fixtures/`, not in the repository) read at run time; when they
//! are absent the tests skip.

use super::diffusion::{DiffusionMap, Neighbours};
use legume_numeric::matrix::common_io::{read_lines, read_lines_of_words_delim, ReadLinesOut};
use nalgebra::DMatrix;

pub(crate) const DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/src/tests/fixtures");
const PSEUDOTIME: &str = "scanpy_pseudotime_reference.tsv.gz";
pub(crate) const CONNECTIVITY: &str = "scanpy_type_connectivity_reference.tsv";

/// The column named `name` in `table`'s header.
pub(crate) fn col(table: &ReadLinesOut<Box<str>>, name: &str) -> usize {
    table
        .header
        .iter()
        .position(|h| h.as_ref() == name)
        .unwrap_or_else(|| panic!("no column {name}"))
}

/// Prepared geometry, labels, root, pseudotime and eigenvalues of the subset.
pub(crate) struct ScanpyReference {
    pub(crate) labels: Vec<Box<str>>,
    pub(crate) geometry: DMatrix<f32>,
    pub(crate) root: usize,
    pub(crate) pseudotime: Vec<f32>,
    pub(crate) evals: Vec<f32>,
}

impl ScanpyReference {
    /// The fixture; `None` (and a note) when it is absent.
    pub(crate) fn read() -> Option<Self> {
        let path = format!("{DIR}/{PSEUDOTIME}");
        if !std::path::Path::new(&path).exists() {
            eprintln!("skipping: {path} is absent");
            return None;
        }
        let comment = read_lines(&path).unwrap().into_iter().next().unwrap();
        let evals = comment
            .split("evals=")
            .nth(1)
            .unwrap()
            .split(',')
            .map(|v| v.trim().parse().unwrap())
            .collect();
        let table = read_lines_of_words_delim(&path, "\t", 0).unwrap();
        let dims: Vec<usize> = (0..)
            .map_while(|j| {
                table
                    .header
                    .iter()
                    .position(|h| h.as_ref() == format!("x{j}"))
            })
            .collect();
        let (label, is_root, pt) = (
            col(&table, "label"),
            col(&table, "is_root"),
            col(&table, "pseudotime"),
        );
        let rows = &table.lines;
        Some(Self {
            labels: rows.iter().map(|r| r[label].clone()).collect(),
            geometry: DMatrix::from_fn(rows.len(), dims.len(), |i, j| {
                rows[i][dims[j]].parse().unwrap()
            }),
            root: rows
                .iter()
                .position(|r| r[is_root].as_ref() == "true")
                .unwrap(),
            pseudotime: rows.iter().map(|r| r[pt].parse().unwrap()).collect(),
            evals,
        })
    }

    /// The kNN lists at scanpy's `n_neighbors = 15`.
    pub(crate) fn neighbours(&self) -> Neighbours {
        Neighbours::new(&self.geometry, 15).unwrap()
    }

    /// The diffusion map with as many components as the fixture has, all of
    /// them in the DPT distance, as the fixture's pseudotime was made.
    pub(crate) fn diffusion_map(&self) -> DiffusionMap {
        let n = self.evals.len();
        DiffusionMap::new(&self.neighbours(), n, n, 0.0).unwrap()
    }
}
