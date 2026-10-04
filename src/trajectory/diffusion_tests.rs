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
    let hsc: Vec<usize> = (0..g.labels.len())
        .filter(|&i| g.labels[i].as_ref() == "HSC")
        .collect();
    assert_eq!(g.diffusion_map().medoid(&hsc), Some(g.root));
}

#[test]
fn exact_duplicate_cells_are_refused() {
    let mut x = DMatrix::<f32>::from_fn(40, 3, |i, j| {
        (i * 7 + j * 3) as f32 % 11.0 + i as f32 * 0.01
    });
    let first = x.row(0).into_owned();
    x.set_row(1, &first);
    let err = Neighbours::new(&x, 5)
        .err()
        .expect("duplicates must be refused");
    assert!(err.to_string().contains("coincide"), "{err}");
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
    let dm = DiffusionMap::new(&Neighbours::new(&x.mat, 15).unwrap(), 15).unwrap();
    let got = dm.pseudotime(root);
    let secs = t0.elapsed().as_secs_f64();
    let rho = spearman(&got, &want);
    eprintln!(
        "full bench: {} cells, Spearman {rho:.5}, {secs:.2} s",
        x.rows.len()
    );
    assert!(rho >= 0.98, "Spearman with scanpy on the full bench {rho}");
}
