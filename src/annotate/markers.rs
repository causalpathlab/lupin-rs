//! The `gene<TAB>celltype` marker panel as a gene × cell-type matrix aligned
//! to the data's gene order.

use crate::annotate::gene_rows::GeneRows;
use legume_numeric::matrix::common_io::read_lines;
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
/// [`label_key`] form. Lines split on the tab; a line with none splits at its
/// first comma, so the label keeps anything after it. Blank lines, `#`
/// comments and a `gene`/`symbol` header are skipped. Reads through gzip.
pub fn read_marker_pairs(path: &str) -> anyhow::Result<Vec<(Box<str>, Box<str>)>> {
    let lines =
        read_lines(path).map_err(|e| anyhow::anyhow!("reading marker panel {path}: {e}"))?;
    Ok(lines
        .iter()
        .filter_map(|line| {
            let line = line.trim();
            let (gene, label) = line.split_once('\t').or_else(|| line.split_once(','))?;
            let (gene, label) = (gene.trim(), label.split('\t').next().unwrap_or(label));
            let skip = gene.is_empty()
                || gene.starts_with('#')
                || matches!(gene.to_lowercase().as_str(), "gene" | "symbol");
            let label = label_key(label);
            (!skip && !label.is_empty()).then(|| (Box::from(gene), label.into_boxed_str()))
        })
        .collect())
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
    fn labels_differing_only_in_separators_are_one_type() {
        for l in ["CT 1, a", "CT_1_a", "  CT,1 a ", "CT__1 ,a"] {
            assert_eq!(label_key(l), "CT_1_a", "{l}");
        }
    }

    #[test]
    fn a_panel_splits_on_tabs_or_the_first_comma_and_keys_its_labels() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("p.tsv");
        std::fs::write(
            &tsv,
            "gene\tcelltype\n# note\nGENE1\tCT1, sub a\nGENE2\tCT 2\tignored\n\nGENE3,CT3, sub b\n",
        )
        .unwrap();
        let pairs = read_marker_pairs(&tsv.to_string_lossy()).unwrap();
        let pairs: Vec<(&str, &str)> = pairs.iter().map(|(g, t)| (&**g, &**t)).collect();
        assert_eq!(
            pairs,
            [
                ("GENE1", "CT1_sub_a"),
                ("GENE2", "CT_2"),
                ("GENE3", "CT3_sub_b")
            ]
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
