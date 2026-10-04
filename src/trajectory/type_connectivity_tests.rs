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
    let conn = connectivity(&g.neighbours(), &group, names.len());
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
