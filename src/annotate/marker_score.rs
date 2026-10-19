//! The embedding marker score on a run's matrices: each cell's best panel
//! type and its margin over the runner-up. The score is
//! legume-graph-embedding's `type_annotation::embedding_score`, a module
//! score (Tirosh et al. 2016, *Science* 352:189; SCANPY, Wolf et al. 2018,
//! *Genome Biology* 19:15) read in closed form on the run's feature and cell
//! embeddings, centred at the balanced centre of its clusters. A round tables
//! the per-cell calls by cluster ([`crate::annotate::rounds::call_shares`]).

use crate::annotate::args::MarkerScoreArgs;
use crate::annotate::outputs::{
    clean_outputs, labels_path, AnnotationOutputs, MARKER_SCORE_PARQUET,
};
use anyhow::Result;
use enrichment::UNASSIGNED_LABEL;
use graph_embedding_util::type_annotation::embedding_score::{
    balanced_centre, calls, cell_covariance, scores, FeatureSpace,
};
use legume_numeric::matrix::dense_mat_io::{Mat, MatWithNames};
use legume_numeric::matrix::parquet::{write_named_table, Column};

/// Per cell: the index of its best type in `types` (`None` without a
/// finite score) and the margin over the runner-up (NaN without one).
pub struct MarkerScore {
    pub types: Vec<Box<str>>,
    pub calls: Vec<Option<usize>>,
    pub margins: Vec<f64>,
}

/// Score `cells` (`[N × H]`) against the `markers` panel on the feature
/// embedding `feat` (`[G × H]`, rows named by gene), centred at the
/// balanced centre of `groups` (cell → cluster, `usize::MAX` left out).
pub fn score_cells(
    feat: &MatWithNames<Mat>,
    cells: &Mat,
    groups: &[usize],
    markers: &str,
) -> Result<MarkerScore> {
    let rho = feat.mat.map(f64::from);
    let z = cells.map(f64::from);
    let space = FeatureSpace::new(&rho, &cell_covariance(&z)?, &feat.rows)?;
    let panel = space.marker_sets(markers, 0.0)?;
    let dirs = space.directions(&panel.sets);
    let s = scores(&z, &balanced_centre(&z, groups)?, &dirs);
    let (calls, margins) = calls(&s, &dirs)
        .into_iter()
        .map(|c| (c.call, c.margin))
        .unzip();
    Ok(MarkerScore {
        types: panel.type_names,
        calls,
        margins,
    })
}

/// Score `cells` against `args.markers` on `feat`, centred on `groups`, and
/// write the pass's files under `args.out`: the shared per-cell label
/// contract (`argmax.parquet` beside `argmax.tsv` and `membership.tsv`; the
/// probability is NaN, as the
/// score has none) and each cell's cluster, call and margin
/// (`marker_score.parquet`, not a posterior: never `annotation`). The
/// round's cluster table is the caller's; the
/// round view tables each cluster's calls from the per-cell ones.
pub fn run(
    args: &MarkerScoreArgs,
    feat: &MatWithNames<Mat>,
    cells: &MatWithNames<Mat>,
    groups: &[usize],
) -> Result<AnnotationOutputs> {
    let out = args.out.as_ref();
    legume_numeric::matrix::common_io::mkdir_parent(out)?;
    if !args.no_clean {
        clean_outputs(
            out,
            graph_embedding_util::type_annotation::CELL_LABEL_FILE_SUFFIXES,
        );
        clean_outputs(out, &[MARKER_SCORE_PARQUET]);
    }
    let (table, call) = score_and_write(&args.markers, out, feat, cells, groups)?;
    let probs = vec![f32::NAN; cells.rows.len()];
    graph_embedding_util::type_annotation::write_cell_labels(out, &cells.rows, &call, &probs)?;

    Ok(AnnotationOutputs {
        argmax: Some(labels_path(out)),
        marker_score: Some(table),
        ..AnnotationOutputs::default()
    })
}

/// The score as a second reading beside another pass (an enrichment round
/// carries it): score `cells` and write only the per-cell table
/// (`marker_score.parquet`) under `out`, leaving the round's call alone.
/// Returns the table's path.
pub fn alongside(
    markers: &str,
    out: &str,
    feat: &MatWithNames<Mat>,
    cells: &MatWithNames<Mat>,
    groups: &[usize],
) -> Result<String> {
    Ok(score_and_write(markers, out, feat, cells, groups)?.0)
}

/// Score the cells and write their cluster, call and margin; the table's
/// path and the calls.
fn score_and_write(
    markers: &str,
    out: &str,
    feat: &MatWithNames<Mat>,
    cells: &MatWithNames<Mat>,
    groups: &[usize],
) -> Result<(String, Vec<Box<str>>)> {
    let score = score_cells(feat, &cells.mat, groups, markers)?;
    let call: Vec<Box<str>> = score
        .calls
        .iter()
        .map(|c| {
            c.map_or(UNASSIGNED_LABEL, |t| score.types[t].as_ref())
                .into()
        })
        .collect();
    let table = format!("{out}{MARKER_SCORE_PARQUET}");
    let cluster: Vec<i32> = groups
        .iter()
        .map(|&k| i32::try_from(k).unwrap_or(-1))
        .collect();
    let margin: Vec<f32> = score.margins.iter().map(|&m| m as f32).collect();
    write_named_table(
        &table,
        "cell",
        &cells.rows,
        &[
            ("cluster".into(), Column::I32(&cluster)),
            ("call".into(), Column::Str(&call)),
            ("margin".into(), Column::F32(&margin)),
        ],
    )?;
    Ok((table, call))
}

/// A planted toy run for tests: six genes (g0..g2 load on axis 0, g3..g5 on
/// axis 1), twenty cells out on the two axes in clusters 0 and 1, and a
/// panel of A (g0..g2) and B (g3..g5).
#[cfg(test)]
pub(crate) mod toy {
    use super::*;

    pub(crate) fn matrices() -> (MatWithNames<Mat>, MatWithNames<Mat>, Vec<usize>) {
        let genes: Vec<Box<str>> = (0..6).map(|g| format!("g{g}").into()).collect();
        let mut rho = Mat::zeros(6, 2);
        for g in 0..6 {
            rho[(g, g / 3)] = 1.0 + 0.1 * g as f32;
            rho[(g, 1 - g / 3)] = 0.05 * g as f32;
        }
        let cells: Vec<Box<str>> = (0..20).map(|i| format!("c{i:02}").into()).collect();
        let mut z = Mat::zeros(20, 2);
        for i in 0..20 {
            z[(i, i / 10)] = 2.0 + 0.1 * (i % 10) as f32;
            z[(i, 1 - i / 10)] = 0.2 * ((i % 3) as f32 - 1.0);
        }
        let h: Vec<Box<str>> = vec!["h0".into(), "h1".into()];
        (
            MatWithNames {
                mat: rho,
                rows: genes,
                cols: h.clone(),
            },
            MatWithNames {
                mat: z,
                rows: cells,
                cols: h,
            },
            (0..20).map(|i| i / 10).collect(),
        )
    }

    /// The panel, written in `dir`; its path.
    pub(crate) fn panel(dir: &std::path::Path) -> String {
        let path = dir.join("panel.tsv");
        std::fs::write(&path, "g0\tA\ng1\tA\ng2\tA\ng3\tB\ng4\tB\ng5\tB\n").unwrap();
        path.to_string_lossy().into_owned()
    }

    /// The toy as a bge run in `dir`: ρ, the cell embedding, the run's own
    /// clusters, the panel and the manifest. Returns the manifest's path and
    /// the panel's.
    pub(crate) fn write_run(dir: &std::path::Path) -> (String, String) {
        use legume_numeric::matrix::traits::IoOps;
        let at = |n: &str| dir.join(n).to_string_lossy().into_owned();
        let (feat, cells, groups) = matrices();
        feat.mat
            .to_parquet_with_names(
                &at("run.feature_embedding.parquet"),
                (Some(&feat.rows[..]), Some("gene")),
                None,
            )
            .unwrap();
        cells
            .mat
            .to_parquet_with_names(
                &at("run.cell_embedding.parquet"),
                (Some(&cells.rows[..]), Some("cell")),
                None,
            )
            .unwrap();
        let ids = crate::manifest::rounds::cluster_ids(&groups);
        crate::manifest::rounds::write_clusters(&at("run.clusters.parquet"), &cells.rows, &ids)
            .unwrap();
        let manifest = serde_json::json!({
            "version": 2, "kind": "bge", "prefix": "run",
            "outputs": {"feature_embedding": "run.feature_embedding.parquet",
                        "cell_embedding": "run.cell_embedding.parquet"},
            "cluster": {"clusters": "run.clusters.parquet"}
        });
        std::fs::write(at("run.senna.json"), manifest.to_string()).unwrap();
        (at("run.senna.json"), panel(dir))
    }

    /// The toy run of [`write_run`] with its raw counts beside it (the
    /// toy's genes and six background genes g6..g11 over its cells, each
    /// cluster's three markers high) and its genes doubling as the ones
    /// placed among the cells: a run every method reads. Returns the
    /// manifest's path and the panel's.
    pub(crate) fn write_run_with_counts(dir: &std::path::Path) -> (String, String) {
        use data_beans::sparse_io::{create_sparse_from_dmatrix, SparseIoBackend};
        let (manifest, panel) = write_run(dir);
        let (_, cells, _) = matrices();
        let genes: Vec<Box<str>> = (0..12).map(|g| format!("g{g}").into()).collect();
        let counts = nalgebra::DMatrix::<f32>::from_fn(12, 20, |g, i| match g {
            0..=5 if g / 3 == i / 10 => 20.0 + ((g + i) % 5) as f32,
            0..=5 => ((g + i) % 2) as f32,
            _ => 5.0 + ((g * i) % 3) as f32,
        });
        let path = dir.join("counts.zarr").to_string_lossy().into_owned();
        let mut backend =
            create_sparse_from_dmatrix(&counts, Some(&path), Some(&SparseIoBackend::Zarr)).unwrap();
        backend.register_row_names_vec(&genes);
        backend.register_column_names_vec(&cells.rows);
        drop(backend);
        let mut m: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&manifest).unwrap()).unwrap();
        m["data"] = serde_json::json!({"input": ["counts.zarr"]});
        m["outputs"]["feature_coembedding"] = "run.feature_embedding.parquet".into();
        std::fs::write(&manifest, m.to_string()).unwrap();
        (manifest, panel)
    }

    /// Score arguments writing under `out`, clustered as the run is.
    pub(crate) fn args(markers: &str, out: &str) -> MarkerScoreArgs {
        MarkerScoreArgs {
            markers: markers.into(),
            out: out.into(),
            clusters: None,
            leiden: crate::manifest::annotate::LeidenArgs::default(),
            no_clean: false,
        }
    }
}

#[cfg(test)]
#[path = "marker_score_tests.rs"]
mod tests;
