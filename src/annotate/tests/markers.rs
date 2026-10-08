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
