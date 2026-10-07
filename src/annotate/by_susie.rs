//! The enrichment pass's SuSiE stage: each cluster's pseudobulk counts
//! regressed on the marker panel with positive single effects
//! ([`super::susie`]), so the types compete for the markers they share.
//!
//! The enrichment's statistics (ES, NES, p, its softmax Q, q-values, the
//! ontology walk) stay as it wrote them; SuSiE makes the call. It rewrites the
//! per-cell labels and annotation (each cell its cluster's PIPs), and adds
//! each type's posterior inclusion probability (PIP), its effect (the log fold
//! of its markers), its share of the explained deviance, and the credible
//! sets.
//! PIPs near 1 are common on pseudobulk counts; the effect and the share
//! rank the types that are in.

use super::inputs::EnrichmentInputs;
use super::outputs::{
    write_cluster_tables, AnnotationOutputs, ARGMAX_TSV, CLUSTER_CELLTYPE_EFFECT,
    CLUSTER_CELLTYPE_EXPLAINED, CLUSTER_CELLTYPE_PIP, CLUSTER_CREDIBLE_SETS,
};
use super::susie::{fit_all, SusieConfig};
use legume_numeric::matrix::common_io::write_lines;
use legume_numeric::matrix::dense_mat_io::{axis_id_names, Mat};
use legume_numeric::matrix::traits::IoOps;
use legume_numeric::matrix::utils::quantiles;
use log::{info, warn};

/// Fit the clusters in `inputs` (raw count sums, not a decoder's expectation)
/// and write the stage's outputs under `out`, pointing `outputs`' label and
/// annotation at what it rewrote and its PIP, effect and share tables at theirs.
pub fn run(
    out: &str,
    cfg: &SusieConfig,
    inputs: &EnrichmentInputs,
    outputs: &mut AnnotationOutputs,
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
    let phi = quantiles(&fits.dispersion, &[0.0, 0.5, 1.0]);
    info!(
        "SuSiE: {n_clusters} clusters × {n_types} cell types over {} marker genes; \
         dispersion median {:.3} (range {:.3}–{:.3})",
        fits.dispersion.len(),
        phi[1],
        phi[0],
        phi[2]
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
    let annotation_path = format!("{out}.annotation.parquet");
    annotation.to_parquet_with_names(
        &annotation_path,
        (Some(&inputs.cell_names), Some("cell")),
        Some(&inputs.celltype_names),
    )?;
    info!("wrote {annotation_path}");
    graph_embedding_util::type_annotation::write_label_tsvs(
        out,
        &inputs.cell_names,
        &labels,
        &probs,
    )?;

    outputs.argmax = Some(format!("{out}{ARGMAX_TSV}"));
    outputs.annotation = Some(annotation_path);
    let mut written = written.into_iter();
    outputs.cluster_celltype_pip = written.next();
    outputs.cluster_celltype_effect = written.next();
    outputs.cluster_celltype_explained = written.next();
    Ok(())
}

#[cfg(test)]
#[path = "tests/by_susie.rs"]
mod tests;
