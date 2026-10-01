//! GO/GMT signature helpers for `annotate --method enrichment`: load + reconcile
//! gene-sets against a gene dictionary, and write a per-group top-N signature
//! TSV from a (group × term) effect matrix.

use data_beans::aux::gene_sets::{read_gaf, read_gmt, GafOpts};
use data_beans::aux::ontology::Ontology;
use data_beans::utilities::name_matching::GeneIndex;
use log::info;
use std::io::Write;

/// Coverage floor below which enrichment is meaningless — fail loudly rather
/// than emit an empty signature.
const MIN_COVERAGE_FRAC: f32 = 0.1;
const MIN_COVERAGE_TERMS: usize = 5;
/// Top terms reported per group in the signature TSV.
const TOP_N: usize = 10;

/// A species whose GO annotations lupin can download.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Species {
    Human,
    Mouse,
}

impl Species {
    /// The species of a run's genes, by vote over its rows: an Ensembl id
    /// prefix (`ENSG`, `ENSMUSG`) when the name has one, else the symbol's
    /// case (human symbols are upper case, mouse symbols capitalized).
    /// `None` when neither side wins.
    #[must_use]
    pub fn detect(row_names: &[Box<str>]) -> Option<Self> {
        let (mut human, mut mouse) = (0usize, 0usize);
        for name in row_names {
            let Some(gene) = super::gene_rows::gene_of(name) else {
                continue;
            };
            let upper = gene.to_ascii_uppercase();
            if upper.starts_with("ENSMUSG") {
                mouse += 1;
                continue;
            }
            if upper.starts_with("ENSG") && upper[4..].starts_with(|c: char| c.is_ascii_digit()) {
                human += 1;
                continue;
            }
            let symbol = gene.rsplit('_').next().unwrap_or(gene);
            let mut letters = symbol.chars().filter(char::is_ascii_alphabetic);
            let Some(first) = letters.next() else {
                continue;
            };
            let rest: Vec<char> = letters.collect();
            if rest.is_empty() {
                continue;
            }
            if first.is_ascii_uppercase() && rest.iter().all(char::is_ascii_uppercase) {
                human += 1;
            } else if first.is_ascii_uppercase() && rest.iter().all(char::is_ascii_lowercase) {
                mouse += 1;
            }
        }
        match human.cmp(&mouse) {
            std::cmp::Ordering::Greater => Some(Self::Human),
            std::cmp::Ordering::Less => Some(Self::Mouse),
            std::cmp::Ordering::Equal => None,
        }
    }

    /// The GO Consortium's annotation file for this species.
    #[must_use]
    pub fn gaf_file(self) -> &'static str {
        match self {
            Self::Human => "goa_human.gaf.gz",
            Self::Mouse => "mgi.gaf.gz",
        }
    }

    /// Where [`Self::gaf_file`] is published.
    #[must_use]
    pub fn gaf_url(self) -> String {
        format!(
            "https://current.geneontology.org/annotations/{}",
            self.gaf_file()
        )
    }
}

/// Reconciled GO/GMT gene-sets ready for scoring.
pub struct GeneSetInputs {
    pub onto: Ontology,
    /// `(term id, member rows into the supplied `gene_names`)`, size-windowed,
    /// sorted by id.
    pub terms: Vec<(Box<str>, Vec<usize>)>,
    /// All matched annotated rows (the background universe).
    pub universe: Vec<usize>,
}

/// Load GO/GMT gene-sets, reconcile to `gene_names` (+ size window + coverage
/// gate). Exactly one of `gaf`/`gmt` must be `Some` (the caller validates this).
#[allow(clippy::too_many_arguments)]
pub fn load_go_gene_sets(
    obo: &str,
    gaf: Option<&str>,
    gmt: Option<&str>,
    no_iea: bool,
    min_gene_set: usize,
    max_gene_set: usize,
    gene_names: &[Box<str>],
) -> anyhow::Result<GeneSetInputs> {
    let onto = Ontology::load_obo(obo)?;
    info!("loaded ontology: {} terms from {obo}", onto.len());

    let gene_sets = if let Some(gaf) = gaf {
        info!("reading GAF gene-sets from {gaf} (no_iea={no_iea})");
        read_gaf(gaf, &GafOpts { no_iea })?.into_gene_sets(Some(&onto))
    } else {
        let gmt = gmt.expect("exactly one of --gaf/--gmt is required");
        info!("reading GMT gene-sets from {gmt}");
        read_gmt(gmt)?
    };
    info!(
        "gene-sets: {} terms, {} genes, {} annotations",
        gene_sets.n_terms(),
        gene_sets.n_genes(),
        gene_sets.n_annotations()
    );

    let idx = GeneIndex::build(gene_names);
    let rec = gene_sets.reconcile(&idx, min_gene_set, Some(max_gene_set));
    rec.log_coverage();
    rec.ensure_coverage(MIN_COVERAGE_FRAC, MIN_COVERAGE_TERMS)?;
    let universe = rec.universe;
    let mut terms: Vec<(Box<str>, Vec<usize>)> = rec.term_rows.into_iter().collect();
    terms.sort_by(|a, b| a.0.cmp(&b.0));
    info!("scoring {} size-windowed terms", terms.len());

    Ok(GeneSetInputs {
        onto,
        terms,
        universe,
    })
}

/// Write a per-group top-`N` GO signature TSV from a `group × term` effect
/// matrix, each group's terms ranked by descending positive effect. `term_ids`
/// indexes the matrix columns and aligns 1:1 with `terms` (for the gene count).
pub fn write_go_signature(
    path: &str,
    onto: &Ontology,
    effect_kt: &enrichment::Mat,
    term_ids: &[Box<str>],
    terms: &[(Box<str>, Vec<usize>)],
    group_axis: &str,
    group_names: &[Box<str>],
) -> anyhow::Result<()> {
    let n_terms = term_ids.len();
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "{group_axis}\trank\tterm_id\tterm_name\teffect\tn_genes")?;
    for (k, gname) in group_names.iter().enumerate() {
        let mut ranked: Vec<(usize, f32)> = (0..n_terms)
            .map(|t| (t, effect_kt[(k, t)]))
            .filter(|&(_, e)| e > 0.0)
            .collect();
        ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        for (rank, &(t, e)) in ranked.iter().take(TOP_N).enumerate() {
            let id = &term_ids[t];
            let name = onto.name(id).unwrap_or(id.as_ref());
            writeln!(
                f,
                "{gname}\t{}\t{id}\t{name}\t{e:.4}\t{}",
                rank + 1,
                terms[t].1.len()
            )?;
        }
    }
    info!("wrote {path}");
    Ok(())
}

#[cfg(test)]
mod species_tests {
    use super::Species;

    fn names(v: &[&str]) -> Vec<Box<str>> {
        v.iter().map(|s| Box::from(*s)).collect()
    }

    #[test]
    fn the_species_comes_from_ensembl_ids_else_symbol_case() {
        let human = names(&["ENSG00000000001_GENE1", "GENE2", "chr1:5/baf/alt"]);
        assert_eq!(Species::detect(&human), Some(Species::Human));
        let mouse = names(&["ENSMUSG00000000001_Gene1", "Gene2/count/spliced", "Gene3"]);
        assert_eq!(Species::detect(&mouse), Some(Species::Mouse));
        assert_eq!(Species::detect(&names(&["GENE1", "Gene2"])), None);
        assert_eq!(Species::detect(&names(&["1", "chr1:5/baf/alt"])), None);
    }

    #[test]
    fn each_species_has_its_annotation_file() {
        assert!(Species::Human.gaf_url().ends_with("/goa_human.gaf.gz"));
        assert!(Species::Mouse.gaf_url().ends_with("/mgi.gaf.gz"));
    }
}
