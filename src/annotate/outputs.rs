//! What an annotation pass erases before it runs and what it reports having
//! written, shared by both subcommands (`annotate --method projection`,
//! `annotate --method enrichment`) so their I/O contract stays in lock-step.
//!
//! Wiring these paths into `manifest.annotate.*` is the caller's job — see
//! [`crate::manifest::annotate`].

use std::path::Path;

/// `{out}.argmax.tsv`: per-cell `cell⇥cell_type⇥probability`, written by every pass.
pub const ARGMAX_TSV: &str = ".argmax.tsv";
/// `{out}.annot.parquet`: projection's per-cell table (`community`, `coarse_label`, …).
pub const ANNOT_PARQUET: &str = ".annot.parquet";
/// `{out}.cluster_term_q.parquet`: projection's cluster × term FDR q.
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
/// Enrichment: cluster × GO/GMT term NES, tested as the cell types are.
pub const CLUSTER_TERM_NES: &str = ".cluster_term_nes.parquet";
/// Enrichment: cluster × GO/GMT term p-values.
pub const CLUSTER_TERM_P: &str = ".cluster_term_p.parquet";
/// Enrichment: cluster × GO/GMT term FDR q (BH over the terms within each
/// cluster), named as [`CLUSTER_CELLTYPE_Q_VALUES`] is.
pub const CLUSTER_TERM_Q_VALUES: &str = ".cluster_term_q_values.parquet";

/// `{prefix}{suffix}` files written by `annotate --method enrichment` (relative to its
/// bare `{out}` prefix). NOTE: this prefix is shared with the training run's
/// artifacts (`{out}.cell_embedding.parquet`, `{out}.senna.json`, …), so the
/// list is EXPLICIT — never a glob — to avoid deleting the embedding/manifest.
pub const ENRICHMENT_OUTPUT_SUFFIXES: &[&str] = &[
    ".annotation.parquet",
    ARGMAX_TSV,
    ".membership.tsv",
    CLUSTER_CELLTYPE_Q,
    CLUSTER_CELLTYPE_ES,
    CLUSTER_CELLTYPE_ES_STD,
    CLUSTER_CELLTYPE_NES,
    CLUSTER_CELLTYPE_P,
    CLUSTER_CELLTYPE_Z,
    CLUSTER_CELLTYPE_Q_VALUES,
    ".cluster_celltype_perm_z.parquet",
    ".cluster_expression.parquet",
    CLUSTER_TERM_NES,
    CLUSTER_TERM_P,
    CLUSTER_TERM_Q_VALUES,
    ".ontology_signature.tsv",
    ".ontology_term_effect.parquet",
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
    pub cluster_celltype_q: Option<String>,
    pub cluster_celltype_es: Option<String>,
    pub cluster_expression: Option<String>,
    pub ontology_assignment: Option<String>,
    pub ontology_node_mass: Option<String>,
    pub ontology_signature: Option<String>,
    pub ontology_term_effect: Option<String>,
    pub cluster_celltype_q_values: Option<String>,
    pub cluster_celltype_p: Option<String>,
    pub cluster_celltype_nes: Option<String>,
    pub cluster_term_q: Option<String>,
    pub marker_embedding: Option<String>,
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
