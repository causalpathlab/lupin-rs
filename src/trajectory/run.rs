//! The `lupin trajectory` command: diffusion pseudotime from a named root
//! type and PAGA connectivity between the run's cell types.

use super::diffusion::{DiffusionMap, Neighbours};
use super::encode_groups;
use super::type_connectivity::connectivity;
use crate::cell_labels::read_cell_labels;
use crate::manifest::run::{load, resolve};
use anyhow::{Context, Result};
use clap::Args;
use legume_numeric::matrix::common_io::mkdir_parent;
use legume_numeric::matrix::dense_mat_io::{axis_id_names, Mat};
use legume_numeric::matrix::parquet::{write_named_table, Column};
use legume_numeric::matrix::traits::IoOps;
use legume_numeric::matrix::utils::median;
use log::info;
use nalgebra::DMatrix;

#[derive(Args, Debug)]
pub struct TrajectoryArgs {
    #[arg(short = 'f', long = "from", help = "Run manifest (or its prefix)")]
    pub from: Box<str>,

    #[arg(short = 'o', long, help = "Output prefix")]
    pub out: Box<str>,

    #[arg(
        long,
        help = "Cell type the trajectory starts from; its medoid is the root cell"
    )]
    pub root: Box<str>,

    #[arg(
        long,
        help = "Per-cell labels (cell<TAB>type); defaults to the run's annotate.argmax"
    )]
    pub labels: Option<Box<str>>,

    #[arg(
        long,
        default_value_t = 15,
        help = "Neighbours per cell, counting the cell itself"
    )]
    pub knn: usize,

    #[arg(long, default_value_t = 15, help = "Diffusion components")]
    pub n_dcs: usize,

    #[arg(
        long,
        default_value_t = 20,
        help = "Types with fewer cells are reported but not given a PAGA row"
    )]
    pub min_cells: usize,
}

/// Settings a recomputation needs, apart from where the data came from.
pub(crate) struct Params {
    pub(crate) knn: usize,
    pub(crate) n_dcs: usize,
    pub(crate) min_cells: usize,
}

/// The cells, their prepared geometry and their types.
pub(crate) struct Inputs {
    pub(crate) cells: Vec<Box<str>>,
    pub(crate) geometry: Mat,
    pub(crate) labels: Vec<Box<str>>,
    pub(crate) names: Vec<Box<str>>,
    pub(crate) group: Vec<usize>,
}

/// What a run computes.
pub(crate) struct Trajectory {
    pub(crate) root: usize,
    pub(crate) pseudotime: Vec<f32>,
    pub(crate) map: DiffusionMap,
    /// PAGA connectivity between every pair of types.
    pub(crate) connectivity: DMatrix<f64>,
    /// Types with at least `min_cells` cells, other than unassigned.
    pub(crate) kept: Vec<usize>,
    /// Kept types in order of median pseudotime: `(median, type)`; a type no
    /// root reaches has an infinite median.
    pub(crate) order: Vec<(f32, usize)>,
}

pub fn run_trajectory(args: &TrajectoryArgs) -> Result<()> {
    let inputs = load_inputs(args)?;
    let params = Params {
        knn: args.knn,
        n_dcs: args.n_dcs,
        min_cells: args.min_cells,
    };
    let t = compute(&inputs, &args.root, &params)?;
    info!(
        "root: {} ({} medoid); eigenvalues {:?}",
        inputs.cells[t.root], args.root, t.map.evals
    );
    log_order(&inputs, &t, params.min_cells);
    let written = write(&inputs, &t, &args.out)?;
    info!("wrote {}", written.join(", "));
    Ok(())
}

pub(crate) fn load_inputs(args: &TrajectoryArgs) -> Result<Inputs> {
    let loaded = load(&args.from)?;
    let x = loaded.prepared_geometry()?;
    let label_path = match args.labels.as_deref() {
        Some(p) => p.to_string(),
        None => loaded
            .manifest
            .annotate
            .argmax
            .as_deref()
            .map(|rel| resolve(&loaded.dir, rel))
            .context("the run has no annotation; annotate it or pass --labels")?,
    };
    let membership = read_cell_labels(&label_path)?;
    let mut missing = 0usize;
    let labels: Vec<Box<str>> = x
        .rows
        .iter()
        .map(|c| {
            membership.get(c).map(Into::into).unwrap_or_else(|| {
                missing += 1;
                enrichment::UNASSIGNED_LABEL.into()
            })
        })
        .collect();
    if missing > 0 {
        info!(
            "{missing}/{} cells are not in {label_path}; they count as {}",
            x.rows.len(),
            enrichment::UNASSIGNED_LABEL
        );
    }
    let (names, group) = encode_groups(&labels);
    Ok(Inputs {
        cells: x.rows,
        geometry: x.mat,
        labels,
        names,
        group,
    })
}

pub(crate) fn compute(inputs: &Inputs, root_type: &str, params: &Params) -> Result<Trajectory> {
    let g = inputs
        .names
        .binary_search_by(|n| n.as_ref().cmp(root_type))
        .ok()
        .with_context(|| {
            format!(
                "no cell is labelled {root_type:?}; the labels are: {}",
                inputs.names.join(", ")
            )
        })?;
    let root_cells: Vec<usize> = (0..inputs.group.len())
        .filter(|&i| inputs.group[i] == g)
        .collect();

    let nb = Neighbours::new(&inputs.geometry, params.knn)?;
    let map = DiffusionMap::new(&nb, params.n_dcs)?;
    let root = map.medoid(&root_cells).context("no root cell")?;
    let pseudotime = map.pseudotime(root);

    // PAGA groups every label, as scanpy requires.
    let connectivity = connectivity(&nb, &inputs.group, inputs.names.len());
    let mut by_type: Vec<Vec<f32>> = vec![Vec::new(); inputs.names.len()];
    for (&g, &t) in inputs.group.iter().zip(&pseudotime) {
        by_type[g].push(t);
    }
    let kept: Vec<usize> = (0..inputs.names.len())
        .filter(|&g| {
            by_type[g].len() >= params.min_cells
                && inputs.names[g].as_ref() != enrichment::UNASSIGNED_LABEL
        })
        .collect();
    let mut order: Vec<(f32, usize)> = kept
        .iter()
        .map(|&g| {
            let finite: Vec<f32> = by_type[g]
                .iter()
                .copied()
                .filter(|t| t.is_finite())
                .collect();
            let m = if finite.is_empty() {
                f32::INFINITY
            } else {
                median(&finite)
            };
            (m, g)
        })
        .collect();
    order.sort_by(|p, q| p.0.total_cmp(&q.0));
    Ok(Trajectory {
        root,
        pseudotime,
        map,
        connectivity,
        kept,
        order,
    })
}

fn log_order(inputs: &Inputs, t: &Trajectory, min_cells: usize) {
    let mut size = vec![0usize; inputs.names.len()];
    inputs.group.iter().for_each(|&g| size[g] += 1);
    info!("median pseudotime by type (≥ {min_cells} cells):");
    for &(m, g) in &t.order {
        info!("  {m:.3}  {} ({})", inputs.names[g], size[g]);
    }
}

/// Writes `{out}.cell_pseudotime.parquet`, `{out}.diffusion.parquet` and
/// `{out}.type_connectivity.parquet`; returns their paths.
pub(crate) fn write(inputs: &Inputs, t: &Trajectory, out: &str) -> Result<Vec<String>> {
    mkdir_parent(out)?;
    let pt_path = format!("{out}.cell_pseudotime.parquet");
    write_named_table(
        &pt_path,
        "cell",
        &inputs.cells,
        &[
            ("pseudotime".into(), Column::F32(&t.pseudotime)),
            ("type".into(), Column::Str(&inputs.labels)),
        ],
    )?;

    let dc = t.map.evecs.map(|v| v as f32);
    let dc_path = format!("{out}.diffusion.parquet");
    dc.to_parquet_with_names(
        &dc_path,
        (Some(&inputs.cells), Some("cell")),
        Some(&axis_id_names("DC", dc.ncols())),
    )?;

    let (mut a, mut b, mut v) = (Vec::new(), Vec::new(), Vec::new());
    for (i, &p) in t.kept.iter().enumerate() {
        for &q in &t.kept[i + 1..] {
            a.push(inputs.names[p].clone());
            b.push(inputs.names[q].clone());
            v.push(t.connectivity[(p, q)] as f32);
        }
    }
    let conn_path = format!("{out}.type_connectivity.parquet");
    write_named_table(
        &conn_path,
        "a",
        &a,
        &[
            ("b".into(), Column::Str(&b)),
            ("connectivity".into(), Column::F32(&v)),
        ],
    )?;
    Ok(vec![pt_path, dc_path, conn_path])
}
