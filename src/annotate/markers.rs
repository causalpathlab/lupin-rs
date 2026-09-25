//! The `gene<TAB>celltype` marker panel as a gene × cell-type matrix aligned
//! to the data's gene order.

use data_beans::utilities::name_matching::GeneIndex;
use legume_numeric::matrix::dense_mat_io::Mat;

/// IDF-weighted gene × cell-type membership plus the sorted cell-type names
/// indexing its columns.
pub struct AnnotInfo {
    pub membership_ga: Mat,
    pub annot_names: Vec<Box<str>>,
}

/// Read a marker TSV and match its genes to `row_names` (exact → symbol →
/// flexible); unmatched markers are logged and dropped. Cell-type names have
/// spaces replaced by `_`.
pub fn build_annotation_matrix(
    marker_gene_path: &str,
    row_names: &[Box<str>],
) -> anyhow::Result<AnnotInfo> {
    let marker_pairs = data_beans::aux::gene_sets::read_membership_pairs(marker_gene_path)?;
    anyhow::ensure!(
        !marker_pairs.is_empty(),
        "empty/invalid marker gene information"
    );

    let normalized: Vec<Box<str>> = marker_pairs
        .iter()
        .map(|(_, t)| t.replace(' ', "_").into_boxed_str())
        .collect();
    let mut annot_names = normalized.clone();
    annot_names.sort_unstable();
    annot_names.dedup();

    let mut membership = Mat::zeros(row_names.len(), annot_names.len());
    let mut matched = 0;
    let mut unmatched = Vec::new();
    let gene_index = GeneIndex::build(row_names);
    for ((gene, _), ty) in marker_pairs.iter().zip(&normalized) {
        let a = annot_names
            .binary_search(ty)
            .expect("every type is in annot_names");
        if let Some(g) = gene_index.match_gene(gene) {
            membership[(g, a)] = 1.0;
            matched += 1;
        } else {
            unmatched.push(gene.clone());
        }
    }

    if !unmatched.is_empty() && unmatched.len() <= 10 {
        log::info!("Unmatched marker genes: {unmatched:?}");
    } else if !unmatched.is_empty() {
        log::info!("{} marker genes not found in dictionary", unmatched.len());
    }

    // w_g = ln(C / c_g): genes every type claims drop out of the score.
    let max_idf = enrichment::markers::apply_idf_weights(&mut membership);
    log::info!(
        "Matched {matched}/{} marker genes to {} cell types (IDF max ln(C) = {max_idf:.3})",
        marker_pairs.len(),
        annot_names.len(),
    );
    Ok(AnnotInfo {
        membership_ga: membership,
        annot_names,
    })
}
