//! The `gene<TAB>celltype` marker panel as a gene × cell-type matrix aligned
//! to the data's gene order.

use crate::annotate::gene_rows::GeneRows;
use anyhow::Context;
use data_beans::aux::gene_sets::read_membership_pairs;
use legume_numeric::matrix::dense_mat_io::Mat;

/// IDF-weighted gene × cell-type membership plus the sorted cell-type names
/// indexing its columns.
pub struct AnnotInfo {
    /// G × C IDF-weighted membership (0 for a gene every type claims).
    pub membership_ga: Mat,
    pub annot_names: Vec<Box<str>>,
    /// Each row's cell types (indices into `annot_names`), unweighted.
    pub support: Vec<Vec<usize>>,
}

/// A cell-type label's canonical form: words split on whitespace, `,` and
/// `_`, joined by `_`. Labels that differ only in those separators name one
/// type (`CT 1`, `CT_1` and `CT, 1` are all `CT_1`).
#[must_use]
pub fn label_key(label: &str) -> String {
    label
        .split(|c: char| c.is_whitespace() || c == ',' || c == '_')
        .filter(|w| !w.is_empty())
        .collect::<Vec<_>>()
        .join("_")
}

/// A marker panel as `(gene, cell type)` pairs, each label in its
/// [`label_key`] form, read as senna reads it ([`read_membership_pairs`]): a
/// line with a tab splits on its tabs, one without is CSV (a quoted label
/// keeps its commas), and cells keep their columns, unquoted. Blank lines,
/// `#`/`%` comments and a `gene`/`symbol` header are skipped, and skipped rows
/// are counted in a warning. A header that names the type column first
/// (`celltype<TAB>gene`, as a round's `celltype_markers.tsv` export has it)
/// reads the two columns the other way round. Reads through gzip.
pub fn read_marker_pairs(path: &str) -> anyhow::Result<Vec<(Box<str>, Box<str>)>> {
    Ok(key_labels(
        read_membership_pairs(path).context("reading marker panel")?,
    ))
}

/// A panel's rows as `(gene, label)` pairs, the label in [`label_key`] form,
/// dropping a label that keys to nothing (`_`, `,`). Rows under a
/// type-then-gene header are turned round, the header dropped.
pub fn key_labels(mut pairs: Vec<(Box<str>, Box<str>)>) -> Vec<(Box<str>, Box<str>)> {
    let type_first = pairs
        .first()
        .is_some_and(|(a, b)| names_type_then_gene(a, b));
    if type_first {
        pairs.remove(0);
    }
    pairs
        .into_iter()
        .map(|(a, b)| if type_first { (b, a) } else { (a, b) })
        .filter_map(|(gene, label)| {
            let label = label_key(&label);
            (!label.is_empty()).then(|| (gene, label.into_boxed_str()))
        })
        .collect()
}

/// Whether a header row names a cell-type column, then a gene column.
fn names_type_then_gene(first: &str, second: &str) -> bool {
    let is = |cell: &str, names: &[&str]| names.iter().any(|n| cell.eq_ignore_ascii_case(n));
    is(first, &["celltype", "cell_type", "type", "label"]) && is(second, &["gene", "symbol"])
}

/// [`read_marker_pairs`] as owned `(gene, cell type)` strings.
pub fn read_panel(path: &str) -> anyhow::Result<Vec<(String, String)>> {
    Ok(read_marker_pairs(path)?
        .into_iter()
        .map(|(g, t)| (g.into_string(), t.into_string()))
        .collect())
}

/// Read a marker TSV and match its genes to `row_names` (exact → symbol →
/// flexible); a matched gene marks every row of that gene
/// ([`GeneRows`]), and unmatched markers are logged and dropped. Cell-type names are
/// in their [`label_key`] form.
pub fn build_annotation_matrix(
    marker_gene_path: &str,
    row_names: &[Box<str>],
) -> anyhow::Result<AnnotInfo> {
    let marker_pairs = read_marker_pairs(marker_gene_path)?;
    annotation_matrix_from_pairs(&marker_pairs, row_names)
}

/// [`build_annotation_matrix`] for `(gene, cell type)` pairs already in memory.
pub fn annotation_matrix_from_pairs(
    marker_pairs: &[(Box<str>, Box<str>)],
    row_names: &[Box<str>],
) -> anyhow::Result<AnnotInfo> {
    anyhow::ensure!(
        !marker_pairs.is_empty(),
        "empty/invalid marker gene information"
    );

    let normalized: Vec<Box<str>> = marker_pairs
        .iter()
        .map(|(_, t)| label_key(t).into_boxed_str())
        .collect();
    let mut annot_names = normalized.clone();
    annot_names.sort_unstable();
    annot_names.dedup();

    let mut membership = Mat::zeros(row_names.len(), annot_names.len());
    let mut matched = 0;
    let mut marked = 0;
    let mut unmatched = Vec::new();
    let gene_rows = GeneRows::build(row_names);
    if let Some(m) = gene_rows.modality_summary() {
        log::info!(
            "Feature rows by modality ({m}); each marker marks every row of its gene ({} genes)",
            gene_rows.n_genes()
        );
    }
    for ((gene, _), ty) in marker_pairs.iter().zip(&normalized) {
        let a = annot_names
            .binary_search(ty)
            .expect("every type is in annot_names");
        if let Some(rows) = gene_rows.match_rows(gene) {
            for &g in rows {
                membership[(g, a)] = 1.0;
            }
            marked += rows.len();
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

    let support: Vec<Vec<usize>> = (0..membership.nrows())
        .map(|g| {
            (0..membership.ncols())
                .filter(|&a| membership[(g, a)] > 0.0)
                .collect()
        })
        .collect();
    // w_g = ln(C / c_g): genes every type claims drop out of the score.
    let max_idf = enrichment::markers::apply_idf_weights(&mut membership);
    log::info!(
        "Matched {matched}/{} marker genes ({marked} rows) to {} cell types (IDF max ln(C) = {max_idf:.3})",
        marker_pairs.len(),
        annot_names.len(),
    );
    Ok(AnnotInfo {
        membership_ga: membership,
        annot_names,
        support,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labels_differing_only_in_separators_are_one_type() {
        for l in ["CT 1, a", "CT_1_a", "  CT,1 a ", "CT__1 ,a"] {
            assert_eq!(label_key(l), "CT_1_a", "{l}");
        }
    }

    #[test]
    fn a_panel_reads_as_senna_reads_it_and_keys_its_labels() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("p.tsv");
        std::fs::write(
            &tsv,
            "gene\tcelltype\n# note\n% note\nGENE1\tCT1, sub a\nGENE2\tCT 2\tignored\n\n\
             GENE3,CT3, sub b\n\"GENE4\",\"CT 4\"\n\tCT5, sub c\n",
        )
        .unwrap();
        let pairs = read_marker_pairs(&tsv.to_string_lossy()).unwrap();
        let pairs: Vec<(&str, &str)> = pairs.iter().map(|(g, t)| (&**g, &**t)).collect();
        assert_eq!(
            pairs,
            [
                ("GENE1", "CT1_sub_a"),
                ("GENE2", "CT_2"),
                ("GENE3", "CT3"),
                ("GENE4", "CT_4")
            ],
            "a tab line splits on its tabs, a comma line on its commas; cells keep their \
             columns, unquoted, and an empty gene skips its row"
        );
    }

    #[test]
    fn a_gzipped_panel_reads_like_the_plain_one() {
        use legume_numeric::matrix::common_io::write_lines;
        let dir = tempfile::tempdir().unwrap();
        let lines: Vec<Box<str>> = ["gene\tcelltype", "GENE1\tCT1", "GENE2\tCT 2"]
            .into_iter()
            .map(Box::from)
            .collect();
        let plain = dir.path().join("p.tsv").to_string_lossy().into_owned();
        let gz = dir.path().join("p.tsv.gz").to_string_lossy().into_owned();
        write_lines(&lines, &plain).unwrap();
        write_lines(&lines, &gz).unwrap();
        assert!(!std::fs::read(&gz).unwrap().starts_with(b"gene"));
        let got = read_marker_pairs(&gz).unwrap();
        assert_eq!(got, read_marker_pairs(&plain).unwrap());
        assert_eq!(got.len(), 2);
    }
}

#[cfg(test)]
#[path = "tests/markers.rs"]
mod support_tests;
