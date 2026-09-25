//! Per-cell label files: `cell<TAB>cell_type[<TAB>…]`, as in the
//! `{out}.argmax.tsv` / `membership.tsv` that annotate and lineage write, or a
//! user-supplied override. Every consumer reads them here so header handling
//! and barcode matching agree.

use anyhow::{Context, Result};
use legume_numeric::matrix::common_io::{read_lines_of_words_delim, ReadLinesOut};
use legume_numeric::matrix::membership::{detect_delimiter, Membership};

/// Read a cell → label file. The delimiter follows the extension (tab unless
/// `.csv`), `.gz` is fine, blank and `#` lines are skipped, and a leading
/// `cell…` header row is dropped. Lookups fall back from `{barcode}@{sample}`
/// to the bare barcode, so either naming joins.
pub fn read_cell_labels(path: &str) -> Result<Membership> {
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
}
