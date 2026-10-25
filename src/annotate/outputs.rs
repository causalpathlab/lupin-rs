//! What an annotation pass erases before it runs and what it reports having
//! written, shared by both subcommands (`annotate --method projection`,
//! `annotate --method enrichment`) so their I/O contract stays in lock-step.
//!
//! Wiring these paths into `manifest.annotate.*` is the caller's job — see
//! [`crate::manifest::annotate`].

use std::path::Path;

/// `{out}.argmax.tsv`: per-cell `cell⇥cell_type⇥probability`, written by every
/// pass beside the parquet, and what older rounds record as `annotate.argmax`.
pub const ARGMAX_TSV: &str = ".argmax.tsv";

/// `name` and, when it is a round's labels file (`{out}.argmax.tsv` or
/// `.argmax.parquet`), every label file the shared writer puts beside it, so
/// the files a manifest names cover the ones it does not.
#[must_use]
pub fn with_label_siblings(name: &str) -> Vec<String> {
    use graph_embedding_util::type_annotation::CELL_LABEL_FILE_SUFFIXES;
    match labels_prefix(name) {
        Some(out) => CELL_LABEL_FILE_SUFFIXES
            .iter()
            .map(|s| format!("{out}{s}"))
            .collect(),
        None => vec![name.to_string()],
    }
}

/// The prefix of a round's labels file (`{out}.argmax.tsv` or
/// `.argmax.parquet`), as [`labels_path`] names it; `None` for another name.
#[must_use]
pub fn labels_prefix(path: &str) -> Option<&str> {
    use graph_embedding_util::type_annotation::CELL_LABELS_SUFFIX;
    path.strip_suffix(ARGMAX_TSV)
        .or_else(|| path.strip_suffix(CELL_LABELS_SUFFIX))
}

/// The per-cell labels file a pass or round under `out` records as
/// `annotate.argmax`: the parquet, with the TSVs still written beside it.
/// Every reader takes either, so an older round's recorded TSV still reads.
#[must_use]
pub fn labels_path(out: &str) -> String {
    use graph_embedding_util::type_annotation::CELL_LABELS_SUFFIX;
    format!("{out}{CELL_LABELS_SUFFIX}")
}

/// The per-cell labels file found under `out`, in the forms a round may
/// record, newest first: the parquet, else an older pass's TSV.
#[must_use]
pub fn labels_under(out: &str) -> Option<String> {
    [labels_path(out), format!("{out}{ARGMAX_TSV}")]
        .into_iter()
        .find(|p| Path::new(p).exists())
}
/// `{out}.annot.parquet`: the retired term-ORA projection's per-cell table
/// (`community`, `coarse_label`, …), read back from older rounds.
pub const ANNOT_PARQUET: &str = ".annot.parquet";
/// Projection score pass: each cell's cluster, call and margin.
pub const MARKER_SCORE_PARQUET: &str = ".marker_score.parquet";
/// `{out}.cluster_term_q.parquet`: the retired term-ORA projection's cluster
/// × term FDR q, read back from older rounds.
pub const CLUSTER_TERM_Q: &str = ".cluster_term_q.parquet";
/// `{out}.cluster_celltype_q_values.parquet`: enrichment's cluster × cell-type FDR q.
pub const CLUSTER_CELLTYPE_Q_VALUES: &str = ".cluster_celltype_q_values.parquet";
/// Enrichment: cluster × cell type Q (each row a softmax over the types that pass FDR).
pub const CLUSTER_CELLTYPE_Q: &str = ".cluster_celltype_q.parquet";
/// Enrichment: cluster × cell type raw enrichment scores.
pub const CLUSTER_CELLTYPE_ES: &str = ".cluster_celltype_es.parquet";
/// Enrichment: cluster × cell type p-values.
pub const CLUSTER_CELLTYPE_P: &str = ".cluster_celltype_p.parquet";
/// Enrichment: cluster × cell type NES (fgsea's normalized enrichment score), the effect size.
pub const CLUSTER_CELLTYPE_NES: &str = ".cluster_celltype_nes.parquet";
/// Enrichment: cluster × cell type z = Φ⁻¹(1 − p), which Q (the share) is a softmax of.
pub const CLUSTER_CELLTYPE_Z: &str = ".cluster_celltype_z.parquet";
/// Enrichment: cluster × cell type ES restandardized by the gene-set null.
pub const CLUSTER_CELLTYPE_ES_STD: &str = ".cluster_celltype_es_std.parquet";
/// Enrichment: cluster × GO/GMT term z of its effect against matched random gene sets.
pub const CLUSTER_TERM_Z: &str = ".cluster_term_z.parquet";
/// Enrichment: cluster × GO/GMT term p-values.
pub const CLUSTER_TERM_P: &str = ".cluster_term_p.parquet";
/// Enrichment: cluster × GO/GMT term FDR q (BH over the terms within each
/// cluster), named as [`CLUSTER_CELLTYPE_Q_VALUES`] is.
pub const CLUSTER_TERM_Q_VALUES: &str = ".cluster_term_q_values.parquet";
/// GO/GMT terms: cluster × term effect, which the term tables sit beside.
pub const ONTOLOGY_TERM_EFFECT: &str = ".ontology_term_effect.parquet";

/// The `{out}` prefix a round's term effect matrix and its tables share.
#[must_use]
pub fn term_prefix(effect_path: &str) -> &str {
    effect_path
        .strip_suffix(ONTOLOGY_TERM_EFFECT)
        .unwrap_or(effect_path)
}

/// GO/GMT terms: each scored term's measured genes (term id, name, genes
/// comma-joined), beside its tables, for asking why a term or a gene.
pub const TERM_GENES: &str = ".term_genes.tsv";

/// SuSiE stage: cluster × cell type posterior inclusion probability.
pub const CLUSTER_CELLTYPE_PIP: &str = ".cluster_celltype_pip.parquet";
/// SuSiE stage: cluster × cell type posterior mean effect (log fold of the markers).
pub const CLUSTER_CELLTYPE_EFFECT: &str = ".cluster_celltype_effect.parquet";
/// SuSiE stage: cluster × cell type share of the log-likelihood its effect gains.
pub const CLUSTER_CELLTYPE_EXPLAINED: &str = ".cluster_celltype_explained.parquet";
/// SuSiE stage: one row per credible set.
pub const CLUSTER_CREDIBLE_SETS: &str = ".cluster_credible_sets.tsv";

/// `{prefix}{suffix}` files written by `annotate --method enrichment` (relative to its
/// bare `{out}` prefix), besides the per-cell label files, which the shared
/// writer's `CELL_LABEL_FILE_SUFFIXES` lists. NOTE: this prefix is shared with the training run's
/// artifacts (`{out}.cell_embedding.parquet`, `{out}.senna.json`, …), so the
/// list is EXPLICIT — never a glob — to avoid deleting the embedding/manifest.
pub const ENRICHMENT_OUTPUT_SUFFIXES: &[&str] = &[
    ".annotation.parquet",
    CLUSTER_CELLTYPE_Q,
    CLUSTER_CELLTYPE_ES,
    CLUSTER_CELLTYPE_ES_STD,
    CLUSTER_CELLTYPE_NES,
    CLUSTER_CELLTYPE_P,
    CLUSTER_CELLTYPE_Z,
    CLUSTER_CELLTYPE_Q_VALUES,
    ".cluster_celltype_perm_z.parquet",
    CLUSTER_CELLTYPE_PIP,
    CLUSTER_CELLTYPE_EFFECT,
    CLUSTER_CELLTYPE_EXPLAINED,
    CLUSTER_CREDIBLE_SETS,
    ".cluster_expression.parquet",
    CLUSTER_TERM_Z,
    CLUSTER_TERM_P,
    CLUSTER_TERM_Q_VALUES,
    TERM_GENES,
    ".ontology_signature.tsv",
    ONTOLOGY_TERM_EFFECT,
    ".ontology_assignment.tsv",
    ".ontology_node_mass.parquet",
    // written by the marker bootstrap of older builds; erased so none outlive it
    ".cluster_celltype_support.parquet",
    ".cluster_qc.tsv",
    ".type_qc.tsv",
];

/// Erase the exact `{prefix}{suffix}` output files (if present) for a fresh
/// re-run. Only the listed files are removed — sibling artifacts (the
/// embedding, the manifest) are never touched. Best-effort: a removal error is
/// logged, not fatal.
pub fn clean_outputs(prefix: &str, suffixes: &[&str]) {
    let mut removed = 0usize;
    for s in suffixes {
        let path = format!("{prefix}{s}");
        if Path::new(&path).exists() {
            match legume_numeric::matrix::common_io::remove_file(&path) {
                Ok(()) => removed += 1,
                Err(e) => log::warn!("--clean: could not remove {path}: {e}"),
            }
        }
    }
    if removed > 0 {
        log::info!("--clean: removed {removed} existing output file(s) under {prefix}");
    }
}

/// Absolute paths of the artifacts an annotation pass produced. Which fields
/// are filled is method-specific: projection emits cluster × term rather than
/// cluster × celltype enrichment, and GO/GMT gene-set mode emits a signature
/// instead of per-cell labels.
///
/// A `None` is meaningful to the manifest adapter: it CLEARS any stale pointer
/// left by an earlier run of a different annotate method.
#[derive(Default)]
pub struct AnnotationOutputs {
    pub argmax: Option<String>,
    pub annotation: Option<String>,
    /// The projection score's per-cell cluster, call and margin.
    pub marker_score: Option<String>,
    /// The gene-set units' per-cell cluster, call and membership.
    pub units_call: Option<String>,
    pub cluster_celltype_q: Option<String>,
    pub cluster_celltype_es: Option<String>,
    /// The SuSiE stage's PIP, effect and explained-share tables, when it made the call.
    pub cluster_celltype_pip: Option<String>,
    pub cluster_celltype_effect: Option<String>,
    pub cluster_celltype_explained: Option<String>,
    pub cluster_expression: Option<String>,
    pub ontology_assignment: Option<String>,
    pub ontology_node_mass: Option<String>,
    pub ontology_signature: Option<String>,
    pub ontology_term_effect: Option<String>,
    pub cluster_celltype_q_values: Option<String>,
    pub cluster_celltype_p: Option<String>,
    pub cluster_celltype_nes: Option<String>,
    /// The cell types as gene-set units, and the cells each pair shares.
    pub celltype_units: Option<String>,
    pub celltype_links: Option<String>,
    /// Each cell's membership of each type's region (cells × types).
    pub celltype_membership: Option<String>,
    /// Each cluster against each ontology term trained as a unit.
    pub term_alignment: Option<String>,
    /// Each of those terms' centre, beside the cell types'.
    pub term_units: Option<String>,
    /// Each cell type's and term's nearest cells, where a layout draws it.
    pub unit_anchors: Option<String>,
    /// Per-cell cluster ids whose values match the cluster tables' `K{id}` rows.
    pub clusters: Option<String>,
    /// The pass's sufficient statistics, for rescoring later rounds.
    pub stats_cache: Option<crate::manifest::run::StatsCache>,
}

/// Write cluster × cell type `tables` under `out`, each as `{out}{suffix}`
/// with rows `rows` and columns `cols`; returns the paths in order.
pub fn write_cluster_tables(
    out: &str,
    rows: &[Box<str>],
    cols: &[Box<str>],
    tables: &[(&legume_numeric::matrix::dense_mat_io::Mat, &str)],
) -> anyhow::Result<Vec<String>> {
    use legume_numeric::matrix::traits::IoOps;
    tables
        .iter()
        .map(|(m, suffix)| {
            let path = format!("{out}{suffix}");
            m.to_parquet_with_names(&path, (Some(rows), Some("cluster")), Some(cols))?;
            log::info!("wrote {path}");
            Ok(path)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_round_records_its_labels_as_the_parquet() {
        assert_eq!(labels_path("dir/run.L1"), "dir/run.L1.argmax.parquet");
    }

    #[test]
    fn a_labels_prefix_is_read_from_either_recorded_form() {
        for name in ["dir/run.L1.argmax.parquet", "dir/run.L1.argmax.tsv"] {
            assert_eq!(labels_prefix(name), Some("dir/run.L1"), "{name}");
            let mut siblings = with_label_siblings(name);
            siblings.sort();
            assert_eq!(
                siblings,
                [
                    "dir/run.L1.argmax.parquet",
                    "dir/run.L1.argmax.tsv",
                    "dir/run.L1.membership.tsv"
                ],
                "{name}"
            );
        }
        assert_eq!(labels_prefix("dir/run.L1.membership.tsv"), None);
        for other in [
            "dir/run.L1.fine_argmax.tsv",
            "dir/run.L1.fine_argmax.parquet",
        ] {
            assert_eq!(labels_prefix(other), None, "{other}");
        }
    }
}
