//! The marker matrix's support for [`super`].

use super::*;

#[test]
fn the_support_keeps_a_gene_every_type_claims() {
    let pairs: Vec<(Box<str>, Box<str>)> = [("G0", "A"), ("G1", "A"), ("G1", "B"), ("G2", "B")]
        .iter()
        .map(|(g, t)| (Box::from(*g), Box::from(*t)))
        .collect();
    let rows: Vec<Box<str>> = ["G0", "G1", "G2", "G3"]
        .iter()
        .map(|g| Box::from(*g))
        .collect();
    let info = annotation_matrix_from_pairs(&pairs, &rows).unwrap();
    // G1 marks both types: its IDF weight is ln(2/2) = 0, but it is a marker.
    assert_eq!(info.membership_ga[(1, 0)], 0.0);
    assert_eq!(info.support, vec![vec![0], vec![0, 1], vec![1], vec![]]);
}

#[test]
fn a_gene_of_several_rows_is_marked_on_its_level_row_alone() {
    // Marked on every row, a gene counts for as many markers as it has rows,
    // and its rows move together where the null draws rows one by one.
    let pairs: Vec<(Box<str>, Box<str>)> = [("GENE1", "A"), ("GENE2", "B")]
        .map(|(g, t)| (Box::from(g), Box::from(t)))
        .into();
    let rows: Vec<Box<str>> = [
        "GENE1/m6a/chr1:5/methylated",
        "GENE1/count/unspliced",
        "GENE1/count/spliced",
        "GENE2/m6a/methylated",
        "GENE3/count/spliced",
    ]
    .map(Box::from)
    .into();
    let info = annotation_matrix_from_pairs(&pairs, &rows).unwrap();
    let marked: Vec<usize> = (0..rows.len())
        .filter(|&r| !info.support[r].is_empty())
        .collect();
    assert_eq!(marked, [2, 3], "GENE1's spliced counts; GENE2's only row");
    assert!(info.membership_ga[(2, 0)] > 0.0);
}
