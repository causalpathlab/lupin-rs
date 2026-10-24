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
/// The fewest of a term's genes among the data's features for it to count
/// towards the coverage check.
const MIN_MATCHED: usize = 2;
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
    max_gene_set: Option<usize>,
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

    // Coverage is of the gene ids, over every term the data can test; the
    // size window then picks the terms scored, and says so when it leaves none.
    let idx = GeneIndex::build(gene_names);
    let rec = gene_sets.reconcile(&idx, MIN_MATCHED, None);
    rec.log_coverage();
    rec.ensure_coverage(MIN_COVERAGE_FRAC, MIN_COVERAGE_TERMS)?;
    let universe = rec.universe;
    let window = |n: usize| n >= min_gene_set && max_gene_set.is_none_or(|mx| n <= mx);
    let mut terms: Vec<(Box<str>, Vec<usize>)> = (rec.term_rows.into_iter())
        .filter(|(_, rows)| window(rows.len()))
        .collect();
    anyhow::ensure!(
        !terms.is_empty(),
        "no gene-set term shares {min_gene_set}..={} genes with the data's features; \
         widen --go-min-overlap/--go-max-overlap",
        max_gene_set.map_or_else(|| "∞".to_string(), |m| m.to_string()),
    );
    terms.sort_by(|a, b| a.0.cmp(&b.0));
    info!("scoring {} size-windowed terms", terms.len());

    Ok(GeneSetInputs {
        onto,
        terms,
        universe,
    })
}

/// The gene sets a pass asks for: `--gaf`, else (`go`) the GO Consortium's
/// annotations for the species `gene_names` belong to, else `--gmt`; read
/// against `obo` and reconciled to `gene_names`.
/// Terms are kept with `sizes` of their genes among `gene_names`.
pub fn requested_gene_sets(
    obo: &str,
    go: bool,
    gaf: Option<&str>,
    gmt: Option<&str>,
    gene_names: &[Box<str>],
    overlap: std::ops::RangeInclusive<usize>,
) -> anyhow::Result<GeneSetInputs> {
    let gaf = match (gaf, go && gmt.is_none()) {
        (Some(g), _) => Some(g.to_string()),
        (None, true) => {
            let species = anyhow::Context::context(
                Species::detect(gene_names),
                "--go: cannot tell the species from the gene names; pass --gaf instead",
            )?;
            let gaf = crate::manifest::data_files::go_annotations(species, None)?;
            info!("GO annotations for {species:?} genes: {}", gaf.display());
            Some(gaf.to_string_lossy().into_owned())
        }
        (None, false) => None,
    };
    load_go_gene_sets(
        obo,
        gaf.as_deref(),
        gmt,
        !crate::annotate::args::KEEP_IEA,
        *overlap.start(),
        Some(*overlap.end()),
        gene_names,
    )
}

/// Each cluster × term's test, as the cell types are tested.
pub struct TermStats<'a> {
    pub nes_kt: &'a enrichment::Mat,
    pub p_kt: &'a enrichment::Mat,
    pub q_kt: &'a enrichment::Mat,
}

/// Row `k`'s top-[`TOP_N`] terms of a `group × term` effect matrix, by
/// descending positive effect: `(term column, effect)`.
pub fn ranked_terms(effect_kt: &enrichment::Mat, k: usize) -> Vec<(usize, f32)> {
    let mut ranked: Vec<(usize, f32)> = (0..effect_kt.ncols())
        .map(|t| (t, effect_kt[(k, t)]))
        .filter(|&(_, e)| e > 0.0)
        .collect();
    ranked.sort_by(|a, b| b.1.total_cmp(&a.1));
    ranked.truncate(TOP_N);
    ranked
}

/// Write a per-group top-`N` GO signature TSV from a `group × term` effect
/// matrix ([`ranked_terms`]), with their `stats`. `term_ids` indexes the
/// matrix columns and aligns 1:1 with `terms` (for the gene count).
#[allow(clippy::too_many_arguments)]
pub fn write_go_signature(
    path: &str,
    onto: &Ontology,
    effect_kt: &enrichment::Mat,
    stats: &TermStats,
    term_ids: &[Box<str>],
    terms: &[(Box<str>, Vec<usize>)],
    group_axis: &str,
    group_names: &[Box<str>],
) -> anyhow::Result<()> {
    let mut f = std::fs::File::create(path)?;
    writeln!(
        f,
        "{group_axis}\trank\tterm_id\tterm_name\teffect\tn_genes\tnes\tp\tq"
    )?;
    for (k, gname) in group_names.iter().enumerate() {
        for (rank, (t, e)) in ranked_terms(effect_kt, k).into_iter().enumerate() {
            let id = &term_ids[t];
            let name = onto.name(id).unwrap_or(id.as_ref());
            writeln!(
                f,
                "{gname}\t{}\t{id}\t{name}\t{e:.4}\t{}\t{:.3}\t{:.3e}\t{:.3e}",
                rank + 1,
                terms[t].1.len(),
                stats.nes_kt[(k, t)],
                stats.p_kt[(k, t)],
                stats.q_kt[(k, t)],
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

#[cfg(test)]
mod report_tests {
    use super::*;

    #[test]
    fn a_signature_lists_a_groups_positive_terms_by_effect() {
        let dir = tempfile::tempdir().unwrap();
        let obo = dir.path().join("t.obo");
        std::fs::write(
            &obo,
            "[Term]\nid: GO:1\nname: weak term\n\n[Term]\nid: GO:2\nname: strong term\n\n\
             [Term]\nid: GO:3\nname: depleted term\n",
        )
        .unwrap();
        let onto = Ontology::load_obo(obo.to_str().unwrap()).unwrap();
        let ids: Vec<Box<str>> = vec!["GO:1".into(), "GO:2".into(), "GO:3".into()];
        let terms: Vec<(Box<str>, Vec<usize>)> = vec![
            ("GO:1".into(), vec![0, 1]),
            ("GO:2".into(), vec![0, 1, 2]),
            ("GO:3".into(), vec![2, 3]),
        ];
        let effect = enrichment::Mat::from_row_slice(1, 3, &[0.2, 0.9, -0.5]);
        assert_eq!(ranked_terms(&effect, 0), vec![(1, 0.9), (0, 0.2)]);

        let stat = enrichment::Mat::from_row_slice(1, 3, &[0.01, 0.02, 0.5]);
        let stats = TermStats {
            nes_kt: &stat,
            p_kt: &stat,
            q_kt: &stat,
        };
        let tsv = dir.path().join("sig.tsv");
        let groups: Vec<Box<str>> = vec!["K0".into()];
        write_go_signature(
            tsv.to_str().unwrap(),
            &onto,
            &effect,
            &stats,
            &ids,
            &terms,
            "cluster",
            &groups,
        )
        .unwrap();
        let text = std::fs::read_to_string(&tsv).unwrap();
        assert!(text.contains("K0\t1\tGO:2\tstrong term"), "{text}");
        assert!(text.contains("K0\t2\tGO:1\tweak term"), "{text}");
        assert!(!text.contains("depleted"), "{text}");
    }
}
