use super::*;
use crate::trajectory::reference::{col, ScanpyReference};
use legume_numeric::matrix::agreement::spearman;
use legume_numeric::matrix::common_io::read_lines_of_words_delim;

#[test]
fn eigenvalues_match_scanpy_on_the_reference_subset() {
    let Some(g) = ScanpyReference::read() else {
        return;
    };
    let dm = g.diffusion_map();
    let worst = dm
        .evals
        .iter()
        .zip(&g.evals)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f32, f32::max);
    eprintln!("worst eigenvalue gap to scanpy: {worst:.2e}");
    assert!(
        worst < 1e-3,
        "worst eigenvalue gap to scanpy {worst}: {:?}",
        dm.evals
    );
}

#[test]
fn pseudotime_matches_scanpy_on_the_reference_subset() {
    let Some(g) = ScanpyReference::read() else {
        return;
    };
    let rho = spearman(&g.diffusion_map().pseudotime(g.root), &g.pseudotime);
    eprintln!("Spearman with scanpy's pseudotime: {rho:.5}");
    assert!(rho >= 0.99, "Spearman with scanpy's pseudotime {rho}");
}

#[test]
fn medoid_root_is_scanpys_root() {
    let Some(g) = ScanpyReference::read() else {
        return;
    };
    let root_type = &g.labels[g.root];
    let cells: Vec<usize> = (0..g.labels.len())
        .filter(|&i| &g.labels[i] == root_type)
        .collect();
    assert_eq!(g.diffusion_map().medoid(&cells), Some(g.root));
}

#[test]
fn exact_duplicate_cells_share_their_pseudotime() {
    let mut x = DMatrix::<f32>::from_fn(40, 3, |i, j| {
        (i * 7 + j * 3) as f32 % 11.0 + i as f32 * 0.01
    });
    let first = x.row(0).into_owned();
    x.set_row(1, &first);
    x.set_row(2, &first);
    let nb = Neighbours::new(&x, 5).unwrap();
    assert_eq!(nb.n_cells(), 38);
    assert_eq!(nb.rep_of[1], nb.rep_of[0]);
    assert_eq!(nb.rep_of[2], nb.rep_of[0]);
    let dm = DiffusionMap::new(&nb, 5, 5, 0.0).unwrap();
    assert_eq!(dm.evecs.nrows(), 40);
    let pt = dm.pseudotime(10);
    assert!(pt[0].is_finite());
    assert_eq!(pt[0], pt[1]);
    assert_eq!(pt[0], pt[2]);
    assert_eq!(dm.evecs.row(0), dm.evecs.row(2));
}

#[test]
fn dpt_uses_the_first_n_dcs_of_the_computed_components() {
    let x = DMatrix::<f32>::from_fn(60, 3, |i, j| {
        ((i * 13 + j * 5) % 17) as f32 + i as f32 * 0.03
    });
    let nb = Neighbours::new(&x, 6).unwrap();
    let all = DiffusionMap::new(&nb, 8, 8, 0.0).unwrap();
    let first = DiffusionMap::new(&nb, 8, 3, 0.0).unwrap();
    assert_eq!(first.evals.len(), 8, "every computed component is kept");
    assert_eq!(first.evals, all.evals);
    assert_ne!(first.pseudotime(0), all.pseudotime(0));
}

/// The full bench run against scanpy's full-data pseudotime, on local data
/// named by `LUPIN_TRAJECTORY_BENCH=<prefix>` (`{prefix}.senna.json` and
/// `{prefix}.scanpy_dpt.tsv`); runs only on request:
/// `LUPIN_TRAJECTORY_BENCH=… cargo test full_bench -- --ignored --nocapture`.
#[test]
#[ignore]
fn full_bench_matches_scanpy() {
    use crate::manifest::run::load;
    use rustc_hash::FxHashMap;
    let Ok(prefix) = std::env::var("LUPIN_TRAJECTORY_BENCH") else {
        eprintln!("skipping: LUPIN_TRAJECTORY_BENCH is unset");
        return;
    };
    let manifest = format!("{prefix}.senna.json");
    let reference = format!("{prefix}.scanpy_dpt.tsv");
    let x = load(&manifest).unwrap().prepared_geometry().unwrap();
    let table = read_lines_of_words_delim(&reference, "\t", 0).unwrap();
    let (cell, pt) = (col(&table, "cell"), col(&table, "pseudotime"));
    let order: FxHashMap<&str, usize> = x
        .rows
        .iter()
        .enumerate()
        .map(|(i, r)| (r.as_ref(), i))
        .collect();
    let mut want = vec![f32::NAN; x.rows.len()];
    for r in &table.lines {
        want[order[r[cell].as_ref()]] = r[pt].parse().unwrap();
    }
    let root = table
        .lines
        .iter()
        .find(|r| r[pt].parse::<f32>().unwrap() == 0.0)
        .map(|r| order[r[cell].as_ref()])
        .unwrap();

    let t0 = std::time::Instant::now();
    // All 15 components in the distance, as this bench has always run; it
    // must match the n_dcs of the scanpy run it is compared with.
    let dm = DiffusionMap::new(&Neighbours::new(&x.mat, 15).unwrap(), 15, 15, 0.0).unwrap();
    let got = dm.pseudotime(root);
    let secs = t0.elapsed().as_secs_f64();
    let rho = spearman(&got, &want);
    eprintln!(
        "full bench: {} cells, Spearman {rho:.5}, {secs:.2} s",
        x.rows.len()
    );
    assert!(rho >= 0.98, "Spearman with scanpy on the full bench {rho}");
}

#[test]
fn participation_is_one_spread_evenly_and_one_over_n_on_one_cell() {
    assert!((participation([1.0, -1.0, 1.0, -1.0].into_iter()) - 1.0).abs() < 1e-12);
    assert!((participation([0.0, 3.0, 0.0, 0.0].into_iter()) - 0.25).abs() < 1e-12);
    assert_eq!(participation(std::iter::empty()), 0.0);
}

/// 600 cells along a line, and 6 cells (1%) in a clump far off it.
fn line_and_clump() -> DMatrix<f32> {
    DMatrix::<f32>::from_fn(606, 2, |i, j| match (i < 600, j) {
        (true, 0) => i as f32 * 0.05,
        (true, _) => ((i * 7) % 5) as f32 * 0.01,
        (false, 0) => 60.0 + (i % 3) as f32 * 0.01,
        (false, _) => (i % 4) as f32 * 0.01,
    })
}

#[test]
fn a_clumps_own_components_are_left_out_of_the_distance() {
    let nb = Neighbours::new(&line_and_clump(), 6).unwrap();
    let kept = DiffusionMap::new(&nb, 10, 10, 0.0).unwrap();
    assert!(kept.left_out.is_empty(), "0 keeps every component");
    let dm = DiffusionMap::new(&nb, 10, 10, f64::from(DC_MIN_SHARE)).unwrap();
    assert!(!dm.left_out.is_empty(), "the clump's components go");
    for &j in &dm.left_out {
        let share = participation(dm.evecs.column(j).iter().copied());
        assert!(share < f64::from(DC_MIN_SHARE), "DC{j} spread over {share}");
    }
    let main: Vec<usize> = (0..dm.evecs.ncols())
        .filter(|j| !dm.left_out.contains(j))
        .collect();
    assert!(main.len() >= 8, "the line's own components stay: {main:?}");
}
