//! The enrichment pass's SuSiE stage: each cluster's marker genes, as
//! log-folds against their background, regressed on the marker panel with
//! positive single effects ([`super::susie`]), so the types compete for the
//! markers they share.
//!
//! The enrichment's statistics (ES, NES, p, its softmax Q, q-values, the
//! ontology walk) stay as it wrote them; SuSiE makes the call. It rewrites the
//! per-cell labels and annotation (each cell its cluster's PIPs), and adds
//! each type's posterior inclusion probability (PIP), its effect (the log fold
//! of its markers), its share of the log-likelihood the effects gain, and the
//! credible sets.

use super::inputs::EnrichmentInputs;
use super::outputs::{
    write_cluster_tables, AnnotationOutputs, CLUSTER_CELLTYPE_EFFECT, CLUSTER_CELLTYPE_EXPLAINED,
    CLUSTER_CELLTYPE_PIP, CLUSTER_CREDIBLE_SETS,
};
use super::susie::{fit_all, SusieConfig};
use legume_numeric::matrix::common_io::write_lines;
use legume_numeric::matrix::dense_mat_io::{axis_id_names, Mat};
use legume_numeric::matrix::traits::IoOps;
use legume_numeric::matrix::utils::quantiles;
use log::{info, warn};
use std::path::Path;

/// Fit the clusters in `inputs` (raw count sums, not a decoder's expectation)
/// and write the stage's outputs under `out`, pointing its PIP, effect and
/// share tables at theirs; when it makes the call (`calls`), `outputs`' label
/// and annotation too, at what it rewrote. Otherwise its reading stands
/// beside another method's call, whose labels it leaves alone.
pub fn run(
    out: &str,
    cfg: &SusieConfig,
    inputs: &EnrichmentInputs,
    outputs: &mut AnnotationOutputs,
    calls_too: bool,
) -> anyhow::Result<()> {
    anyhow::ensure!(!inputs.celltype_names.is_empty(), "no marker panel");
    let (n_genes, n_clusters) = (inputs.gene_names.len(), inputs.n_clusters);
    let n_types = inputs.celltype_names.len();
    // The panel's unweighted support: the regression credits shared markers
    // itself, so the IDF weights do not enter.
    let fits = fit_all(
        &inputs.gene_sum_kg,
        n_genes,
        n_clusters,
        &inputs.marker_support,
        n_types,
        cfg,
    )?;
    let s: Vec<f32> = quantiles(&fits.scale2, &[0.0, 0.5, 1.0])
        .into_iter()
        .map(f32::sqrt)
        .collect();
    info!(
        "SuSiE: {n_clusters} clusters × {n_types} cell types over {} marker genes; \
         log-fold misfit scale median {:.3} (range {:.3}–{:.3})",
        fits.scale2.len(),
        s[1],
        s[0],
        s[2]
    );

    let cluster_names = axis_id_names("K", n_clusters);
    let mut pip = Mat::zeros(n_clusters, n_types);
    let mut effect = Mat::zeros(n_clusters, n_types);
    let mut explained = Mat::zeros(n_clusters, n_types);
    // Each cluster's call (`ClusterFit::call`).
    let mut calls: Vec<Option<usize>> = Vec::with_capacity(n_clusters);
    let mut sets: Vec<Box<str>> =
        vec!["cluster\tset\tcell_types\tpip\teffect\tfold\texplained".into()];
    for (k, f) in fits.clusters.iter().enumerate() {
        if f.max_rhat > 1.1 {
            warn!(
                "SuSiE: {}'s chains disagree (R̂ {:.2}); raise --mcmc-samples",
                cluster_names[k], f.max_rhat
            );
        }
        for c in 0..n_types {
            pip[(k, c)] = f.pip[c];
            effect[(k, c)] = f.theta[c];
            explained[(k, c)] = f.explained[c];
        }
        calls.push(f.call());
        // The sets by how much of the cluster their lead type explains.
        let mut order: Vec<&Vec<usize>> = f.credible_sets.iter().collect();
        order.sort_by(|a, b| f.explained[b[0]].total_cmp(&f.explained[a[0]]));
        for (i, set) in order.into_iter().enumerate() {
            let names: Vec<&str> = set.iter().map(|&c| &*inputs.celltype_names[c]).collect();
            let lead = set[0];
            sets.push(
                format!(
                    "{}\t{i}\t{}\t{:.4}\t{:.3}\t{:.2}\t{:.3}",
                    cluster_names[k],
                    names.join("|"),
                    f.pip[lead],
                    f.theta[lead],
                    f.theta[lead].exp(),
                    f.explained[lead]
                )
                .into(),
            );
        }
    }

    // The tables first: the labels below replace the enrichment's, so they
    // are written only once everything else is.
    let written = write_cluster_tables(
        out,
        &cluster_names,
        &inputs.celltype_names,
        &[
            (&pip, CLUSTER_CELLTYPE_PIP),
            (&effect, CLUSTER_CELLTYPE_EFFECT),
            (&explained, CLUSTER_CELLTYPE_EXPLAINED),
        ],
    )?;
    let sets_path = format!("{out}{CLUSTER_CREDIBLE_SETS}");
    write_lines(&sets, &sets_path)?;
    info!("wrote {sets_path}");

    // Pointed at once the files it stands for are all written.
    let point = |outputs: &mut AnnotationOutputs| {
        let mut written = written.into_iter();
        outputs.cluster_celltype_pip = written.next();
        outputs.cluster_celltype_effect = written.next();
        outputs.cluster_celltype_explained = written.next();
    };
    if !calls_too {
        point(outputs);
        return Ok(());
    }

    // Each cell takes its cluster's PIPs and call.
    let n_cells = inputs.cell_names.len();
    let mut annotation = Mat::zeros(n_cells, n_types);
    let (mut labels, mut probs) = (Vec::with_capacity(n_cells), Vec::with_capacity(n_cells));
    for (i, &k) in inputs.cluster_labels.iter().enumerate() {
        let call = calls.get(k).copied();
        if call.is_some() {
            for t in 0..n_types {
                annotation[(i, t)] = pip[(k, t)];
            }
        }
        match call.flatten() {
            Some(c) => {
                labels.push(inputs.celltype_names[c].clone());
                probs.push(pip[(k, c)]);
            }
            None => {
                labels.push(enrichment::UNASSIGNED_LABEL.into());
                probs.push(0.0);
            }
        }
    }
    let mut called: std::collections::BTreeMap<&str, usize> = Default::default();
    for c in &calls {
        let name = c.map_or(enrichment::UNASSIGNED_LABEL, |c| &*inputs.celltype_names[c]);
        *called.entry(name).or_default() += 1;
    }
    let summary: Vec<String> = called.iter().map(|(t, n)| format!("{t} {n}")).collect();
    info!("SuSiE calls (clusters): {}", summary.join(", "));
    // The per-cell files replace the enrichment's: written under a staging
    // prefix, then moved into place, so a failure while writing them leaves
    // the enrichment's call whole.
    let staged = format!("{out}.susie-staging");
    let annotation_path = format!("{out}.annotation.parquet");
    annotation.to_parquet_with_names(
        &format!("{staged}.annotation.parquet"),
        (Some(&inputs.cell_names), Some("cell")),
        Some(&inputs.celltype_names),
    )?;
    graph_embedding_util::type_annotation::write_cell_labels(
        &staged,
        &inputs.cell_names,
        &labels,
        &probs,
    )?;
    let labels_files = graph_embedding_util::type_annotation::CELL_LABEL_FILE_SUFFIXES;
    // The cleanup list keeps retired names: move only what was written.
    for suffix in std::iter::once(&".annotation.parquet").chain(labels_files) {
        let from = format!("{staged}{suffix}");
        if Path::new(&from).is_file() {
            std::fs::rename(from, format!("{out}{suffix}"))?;
        }
    }
    info!("wrote {annotation_path} and the per-cell labels");

    outputs.argmax = Some(super::outputs::labels_path(out));
    outputs.annotation = Some(annotation_path);
    point(outputs);
    Ok(())
}

#[cfg(test)]
#[path = "tests/by_susie.rs"]
mod tests;
