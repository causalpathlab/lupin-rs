//! Marker genes resolved to every data row that measures the gene.
//!
//! A plain expression matrix has one row per gene, so matching a marker to a
//! row name is enough. A faba matrix names its rows
//! `{gene}/{modality}[/{subunit}]/{channel}` (`data_beans::aux::feature_rows`),
//! so one gene owns several rows: two count channels, gene-level modification
//! rows, per-site rows. Matching the name alone picks whichever row came
//! first, often a modification site rather than the gene's expression.
//!
//! Here each row is read for its gene and its data type, a marker matches the
//! gene, and it marks every row of that gene: each modality, channel and
//! site. Allele rows are keyed by locus, not gene, and never carry a marker.

use data_beans::aux::feature_rows::{self as fr, parse_feature_row, FeatureRow};
use data_beans::utilities::name_matching::GeneIndex;
use std::collections::{BTreeMap, HashMap};

/// `name` as a feature row of a known modality. A name that carries one more
/// trailing field (a per-file modality suffix from a multiome load) is read
/// without it.
fn parse(name: &str) -> Option<FeatureRow<'_>> {
    let known = |r: &FeatureRow| fr::channels(r.modality).is_some();
    parse_feature_row(name).filter(known).or_else(|| {
        let (head, _) = name.rsplit_once('/')?;
        parse_feature_row(head).filter(known)
    })
}

/// The gene `name` measures: the unit of a feature row, else the name
/// itself; `None` for an allele row, whose unit is a locus.
#[must_use]
pub fn gene_of(name: &str) -> Option<&str> {
    match parse(name) {
        Some(row) if row.modality == fr::BAF => None,
        Some(row) => Some(row.gene),
        None => Some(name),
    }
}

/// Every row of each gene, behind a name index over the genes.
pub struct GeneRows {
    index: GeneIndex,
    rows: Vec<Vec<usize>>,
    /// Rows per modality, empty when no row parsed as a feature row.
    modalities: BTreeMap<Box<str>, usize>,
}

impl GeneRows {
    #[must_use]
    pub fn build(row_names: &[Box<str>]) -> Self {
        let mut genes: Vec<Box<str>> = Vec::new();
        let mut rows: Vec<Vec<usize>> = Vec::new();
        let mut at: HashMap<Box<str>, usize> = HashMap::new();
        let mut modalities = BTreeMap::new();
        for (i, name) in row_names.iter().enumerate() {
            if let Some(row) = parse(name) {
                *modalities.entry(row.modality.into()).or_insert(0) += 1;
            }
            let Some(gene) = gene_of(name) else {
                continue;
            };
            let g = *at.entry(gene.into()).or_insert_with(|| {
                genes.push(gene.into());
                rows.push(Vec::new());
                genes.len() - 1
            });
            rows[g].push(i);
        }
        Self {
            index: GeneIndex::build(&genes),
            rows,
            modalities,
        }
    }

    /// The rows marker `gene` marks, or `None` when no gene matches.
    #[must_use]
    pub fn match_rows(&self, gene: &str) -> Option<&[usize]> {
        self.index.match_gene(gene).map(|g| &*self.rows[g])
    }

    /// Number of distinct genes.
    #[must_use]
    pub fn n_genes(&self) -> usize {
        self.rows.len()
    }

    /// `modality: rows` for each modality seen, or `None` for a matrix with
    /// no feature rows.
    #[must_use]
    pub fn modality_summary(&self) -> Option<String> {
        (!self.modalities.is_empty()).then(|| {
            self.modalities
                .iter()
                .map(|(m, n)| format!("{m}: {n}"))
                .collect::<Vec<_>>()
                .join(", ")
        })
    }
}

/// The count channels a faba gene's level can be read from.
pub const CHANNELS: [&str; 3] = [fr::SPLICED, fr::UNSPLICED, fr::TOTAL];

/// Each gene's one row to read its level from, where a gene owns several
/// (a faba matrix): its gene-level count row of the channel asked for,
/// else of [`fr::SPLICED`], [`fr::TOTAL`], then any other channel, else
/// any gene-level row, else a per-site row (ties by row order); a plain
/// row is its own.
/// Markers match genes as [`GeneRows`] matches them: by case, `ENSG…_SYMBOL`
/// either way, and without the `/{modality}/{channel}` suffix.
pub struct LevelRows {
    index: GeneIndex,
    /// Per gene, its level row.
    level: Vec<usize>,
    /// Each name asked for, matched: a name that matches no gene is
    /// matched by a scan of every gene, which a screen redrawn many times
    /// a second should not repeat.
    asked: std::sync::Mutex<HashMap<Box<str>, Option<usize>>>,
}

impl LevelRows {
    #[must_use]
    pub fn build(row_names: &[Box<str>], channel: &str) -> Self {
        let genes = GeneRows::build(row_names);
        let rank = |r: usize| match parse(&row_names[r]) {
            None => 0,
            Some(row) if row.subunit.is_some() => 6,
            Some(row) if row.modality != fr::COUNT => 5,
            Some(row) if row.channel == channel => 0,
            Some(row) if row.channel == fr::SPLICED => 1,
            Some(row) if row.channel == fr::TOTAL => 2,
            Some(_) => 3,
        };
        let level = genes
            .rows
            .iter()
            .map(|rows| {
                *rows
                    .iter()
                    .min_by_key(|&&r| rank(r))
                    .expect("a gene has a row")
            })
            .collect();
        Self {
            index: genes.index,
            level,
            asked: std::sync::Mutex::default(),
        }
    }

    /// The row gene `gene`'s level is read from; `None` when no gene
    /// matches.
    #[must_use]
    pub fn row(&self, gene: &str) -> Option<usize> {
        let mut asked = self.asked.lock().unwrap_or_else(|e| e.into_inner());
        if let Some(&r) = asked.get(gene) {
            return r;
        }
        let r = self.index.match_gene(gene).map(|g| self.level[g]);
        asked.insert(gene.into(), r);
        r
    }

    /// `row_names` with each gene's level row (by [`fr::SPLICED`]) kept
    /// and its other rows blank: a [`GeneIndex`] built over them, which
    /// takes a key's first row, then finds each gene at its level row.
    #[must_use]
    pub fn level_names(row_names: &[Box<str>]) -> Vec<Box<str>> {
        let mut out: Vec<Box<str>> = vec![Box::from(""); row_names.len()];
        for &r in Self::build(row_names, fr::SPLICED).rows() {
            out[r] = row_names[r].clone();
        }
        out
    }

    /// Every gene's level row, in the order the genes first appear.
    #[must_use]
    pub fn rows(&self) -> &[usize] {
        &self.level
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(names: &[&str]) -> Vec<Box<str>> {
        names.iter().map(|s| Box::from(*s)).collect()
    }

    #[test]
    fn a_marker_marks_every_row_of_its_gene() {
        let names = rows(&[
            "ENSG1_GENE1/m6a/chr1:100/methylated",
            "ENSG1_GENE1/m6a/methylated",
            "ENSG1_GENE1/count/spliced",
            "ENSG1_GENE1/count/unspliced",
            "ENSG1_GENE1/count/total",
            "chr1:500/baf/alt",
        ]);
        let gr = GeneRows::build(&names);
        assert_eq!(gr.match_rows("GENE1"), Some(&[0, 1, 2, 3, 4][..]));
        assert_eq!(gr.match_rows("ENSG1"), Some(&[0, 1, 2, 3, 4][..]));
        assert_eq!(gr.match_rows("chr1:500"), None);
        assert_eq!(gr.n_genes(), 1);
        assert_eq!(
            gr.modality_summary().as_deref(),
            Some("baf: 1, count: 3, m6a: 2")
        );
    }

    #[test]
    fn genes_without_counts_match_their_modification_rows() {
        let names = rows(&[
            "GENE1/atoi/chr1:5/edited",
            "GENE1/atoi/edited",
            "GENE1/atoi/unedited",
            "GENE2/m6a/chr2:9/methylated",
            "GENE3/apa/proximal",
        ]);
        let gr = GeneRows::build(&names);
        assert_eq!(gr.match_rows("GENE1"), Some(&[0, 1, 2][..]));
        assert_eq!(gr.match_rows("GENE2"), Some(&[3][..]));
        assert_eq!(gr.match_rows("GENE3"), Some(&[4][..]));
    }

    #[test]
    fn plain_and_suffixed_names_resolve_too() {
        let names = rows(&[
            "ENSG1_GENE1",
            "GENE2",
            "ENSG3_GENE3/count/spliced/rna",
            "GENE4/rna",
        ]);
        let gr = GeneRows::build(&names);
        assert_eq!(gr.match_rows("GENE1"), Some(&[0][..]));
        assert_eq!(gr.match_rows("gene2"), Some(&[1][..]));
        assert_eq!(gr.match_rows("GENE3"), Some(&[2][..]));
        assert_eq!(gr.match_rows("GENE4"), Some(&[3][..]));
        assert_eq!(gr.modality_summary().as_deref(), Some("count: 1"));
    }

    #[test]
    fn a_genes_level_is_read_from_its_spliced_count_row_by_default() {
        let names = rows(&[
            "ENSG1_GENE1/m6a/chr1:100/methylated",
            "ENSG1_GENE1/count/unspliced",
            "ENSG1_GENE1/m6a/methylated",
            "ENSG1_GENE1/count/spliced",
            "GENE2/count/unspliced",
            "GENE3/m6a/chr2:9/methylated",
            "GENE3/m6a/methylated",
            "GENE4",
        ]);
        let spliced = LevelRows::build(&names, fr::SPLICED);
        assert_eq!(spliced.row("GENE1"), Some(3));
        assert_eq!(spliced.row("gene1"), Some(3), "by case too");
        assert_eq!(spliced.row("ENSG1"), Some(3));
        assert_eq!(spliced.row("ENSG1_GENE1"), Some(3));
        assert_eq!(spliced.row("GENE2"), Some(4), "the count row it has");
        assert_eq!(spliced.row("GENE3"), Some(6), "gene level before a site");
        assert_eq!(spliced.row("GENE4"), Some(7));
        assert_eq!(spliced.row("GENE5"), None);
        assert_eq!(spliced.rows(), [3, 4, 6, 7], "one row per gene");

        // A gene index over the level names finds the spliced row, not the
        // first row of the gene.
        let index = GeneIndex::build(&LevelRows::level_names(&names));
        assert_eq!(index.match_gene("GENE1"), Some(3));
        assert_eq!(index.match_gene("GENE3"), Some(6));

        let unspliced = LevelRows::build(&names, fr::UNSPLICED);
        assert_eq!(unspliced.row("GENE1"), Some(1));
        // A channel the gene lacks: its spliced row, whatever the order.
        let total = LevelRows::build(&names, fr::TOTAL);
        assert_eq!(total.row("GENE1"), Some(3));
    }

    #[test]
    fn a_gene_symbol_with_a_slash_keeps_its_rows_together() {
        let names = rows(&[
            "ENSG1_GENE1/B/count/spliced",
            "ENSG1_GENE1/B/count/unspliced",
        ]);
        let gr = GeneRows::build(&names);
        assert_eq!(gr.match_rows("ENSG1_GENE1/B"), Some(&[0, 1][..]));
    }
}
