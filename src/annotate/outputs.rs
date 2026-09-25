//! What an annotation pass erases before it runs and what it reports having
//! written, shared by both subcommands (`annotate --method projection`,
//! `annotate --method enrichment`) so their I/O contract stays in lock-step.
//!
//! Wiring these paths into `manifest.annotate.*` is the caller's job — see
//! [`crate::annotate_manifest`].

use std::path::Path;

/// `{out}.argmax.tsv`: per-cell `cell⇥cell_type⇥probability`, written by every pass.
pub const ARGMAX_TSV: &str = ".argmax.tsv";
/// `{out}.annot.parquet`: projection's per-cell table (`community`, `coarse_label`, …).
pub const ANNOT_PARQUET: &str = ".annot.parquet";
/// `{out}.cluster_term_q.parquet`: projection's cluster × term FDR q.
pub const CLUSTER_TERM_Q: &str = ".cluster_term_q.parquet";
/// `{out}.cluster_celltype_q_values.parquet`: enrichment's cluster × cell-type FDR q.
pub const CLUSTER_CELLTYPE_Q_VALUES: &str = ".cluster_celltype_q_values.parquet";

/// `{prefix}{suffix}` files written by `annotate --method enrichment` (relative to its
/// bare `{out}` prefix). NOTE: this prefix is shared with the training run's
/// artifacts (`{out}.cell_embedding.parquet`, `{out}.senna.json`, …), so the
/// list is EXPLICIT — never a glob — to avoid deleting the embedding/manifest.
pub const ENRICHMENT_OUTPUT_SUFFIXES: &[&str] = &[
    ".annotation.parquet",
    ARGMAX_TSV,
    ".membership.tsv",
    ".cluster_celltype_q.parquet",
    ".cluster_celltype_es.parquet",
    ".cluster_celltype_es_std.parquet",
    ".cluster_celltype_p.parquet",
    CLUSTER_CELLTYPE_Q_VALUES,
    ".cluster_celltype_perm_z.parquet",
    ".cluster_expression.parquet",
    ".ontology_assignment.tsv",
    ".ontology_node_mass.parquet",
    // marker-panel stability bootstrap
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
}
