//! Per-cell label files: `cell<TAB>cell_type[<TAB>…]`, as in the
//! `{out}.argmax.parquet` / `argmax.tsv` / `membership.tsv` that annotate
//! writes, or a user-supplied override. Every consumer reads them here so
//! header handling and barcode matching agree.

use anyhow::{Context, Result};
use legume_numeric::matrix::common_io::{read_lines_of_words_delim, ReadLinesOut};
use legume_numeric::matrix::membership::{detect_delimiter, Membership};

/// One cell's label and its probability (NaN when the table has none).
pub type LabelRow = (Box<str>, Box<str>, f32);

/// A labels table with a header, read by column name: `cell`, `cell_type`
/// and, when present, `probability`; parquet or delimited text (`.gz`
/// fine). Rows with an empty cell or label are left out.
pub fn read_label_rows(path: &str) -> Result<Vec<LabelRow>> {
    use legume_numeric::matrix::table::TableReader;
    let t = TableReader::open(path).with_context(|| format!("reading cell labels {path}"))?;
    let mut cols = t.select(&["cell", "cell_type"])?;
    cols.extend(t.column_index("probability"));
    let mut out = Vec::new();
    for row in t.rows(&cols)? {
        let mut row = row?;
        if row[0].is_empty() || row[1].is_empty() {
            continue;
        }
        let p = row.get(2).and_then(|p| p.parse().ok()).unwrap_or(f32::NAN);
        out.push((std::mem::take(&mut row[0]), std::mem::take(&mut row[1]), p));
    }
    Ok(out)
}

/// Read a cell → label file. A parquet one is read by column name
/// ([`read_label_rows`]). Text: the delimiter follows the extension (tab
/// unless `.csv`), `.gz` is fine, blank and `#` lines are skipped, and a
/// leading `cell…` header row is dropped. Lookups fall back from
/// `{barcode}@{sample}` to the bare barcode, so either naming joins.
pub fn read_cell_labels(path: &str) -> Result<Membership> {
    if legume_numeric::matrix::table::is_parquet_path(path) {
        let pairs: Vec<(Box<str>, Box<str>)> = read_label_rows(path)?
            .into_iter()
            .map(|(cell, label, _)| (cell, label))
            .collect();
        return Ok(Membership::from_pairs(pairs, false).with_delimiter('@'));
    }
    let ReadLinesOut { lines, .. } = read_lines_of_words_delim(path, detect_delimiter(path), -1)
        .with_context(|| format!("reading cell labels {path}"))?;
    let mut rows = lines
        .into_iter()
        .filter(|w| !w.is_empty() && !w[0].is_empty() && !w[0].starts_with('#'))
        .peekable();
    if rows
        .peek()
        .is_some_and(|w| w[0].eq_ignore_ascii_case("cell"))
    {
        rows.next();
    }
    let mut n_short = 0usize;
    let pairs: Vec<(Box<str>, Box<str>)> = rows
        .filter_map(|mut w| {
            if w.len() < 2 {
                n_short += 1;
                return None;
            }
            let label = w.swap_remove(1);
            Some((w.swap_remove(0), label))
        })
        .collect();
    if n_short > 0 {
        log::warn!("{path}: skipped {n_short} line(s) with no label column");
    }
    Ok(Membership::from_pairs(pairs, false).with_delimiter('@'))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn skips_header_and_joins_tagged_barcodes() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.argmax.tsv");
        std::fs::write(
            &path,
            "cell\tcell_type\tprobability\nAAAC\tB\t0.9\n# note\n\nGGGT\tT\t0.8\n",
        )
        .unwrap();
        let mem = read_cell_labels(path.to_str().unwrap()).unwrap();
        assert_eq!(mem.get("AAAC@s1"), Some("B"));
        assert_eq!(mem.get("GGGT"), Some("T"));
        assert_eq!(mem.get("cell"), None);
        assert_eq!(mem.unique_groups(), vec![Box::from("B"), Box::from("T")]);
    }

    #[test]
    fn a_parquet_label_table_reads_as_its_text_does() {
        use legume_numeric::matrix::parquet::{write_named_table, Column};
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.argmax.parquet");
        let types: Vec<Box<str>> = vec!["B".into(), "T".into()];
        write_named_table(
            path.to_str().unwrap(),
            "cell",
            &["AAAC".into(), "GGGT".into()],
            &[
                ("cell_type".into(), Column::Str(&types)),
                ("probability".into(), Column::F32(&[0.9, 0.8])),
            ],
        )
        .unwrap();
        let mem = read_cell_labels(path.to_str().unwrap()).unwrap();
        assert_eq!(mem.get("AAAC@s1"), Some("B"));
        assert_eq!(mem.get("GGGT"), Some("T"));
    }

    #[test]
    fn a_parquet_label_table_is_read_by_column_name_without_empty_labels() {
        use legume_numeric::matrix::parquet::{write_named_table, Column};
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.annot.parquet");
        let community: Vec<Box<str>> = vec!["7".into(), "8".into()];
        let types: Vec<Box<str>> = vec!["B".into(), "".into()];
        write_named_table(
            path.to_str().unwrap(),
            "cell",
            &["AAAC".into(), "GGGT".into()],
            &[
                ("community".into(), Column::Str(&community)),
                ("cell_type".into(), Column::Str(&types)),
            ],
        )
        .unwrap();
        let mem = read_cell_labels(path.to_str().unwrap()).unwrap();
        assert_eq!(
            mem.get("AAAC"),
            Some("B"),
            "cell_type, not the second column"
        );
        assert_eq!(mem.get("GGGT"), None, "an empty label is no label");
    }
}
