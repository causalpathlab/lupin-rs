//! Streaming per-group gene sums over the raw counts.
//!
//! Walks count columns in parallel blocks, accumulating
//! `T[k, g] = Σ_{n ∈ k} y[g, n]` row-major as `Vec<f64>` of length `k · m`.
//! Cells whose label is outside `0..k` are skipped.

use data_beans::sparse_io_vector::SparseIoVec;
use legume_numeric::matrix::dense_mat_io::Mat;
use rayon::prelude::*;

/// Peak group-sum memory allowed across all workers: every worker holds one
/// full set of accumulators, so parallelism is bounded by memory, not cores.
const ACCUMULATOR_BUDGET_BYTES: usize = 1 << 30;

/// Per-group gene sums for one grouping.
pub fn accumulate_gene_sum(
    data_vec: &SparseIoVec,
    labels: &[usize],
    k: usize,
    m: usize,
    block_size: usize,
) -> anyhow::Result<Vec<f64>> {
    let mut out = accumulate_gene_sum_multi(data_vec, &[(labels, k)], m, block_size)?;
    Ok(out.pop().expect("one grouping in, one out"))
}

/// Two groupings (e.g. cluster and batch) in a single sweep over the columns.
pub fn accumulate_gene_sum_pair(
    data_vec: &SparseIoVec,
    labels_a: &[usize],
    k_a: usize,
    labels_b: &[usize],
    k_b: usize,
    m: usize,
    block_size: usize,
) -> anyhow::Result<(Vec<f64>, Vec<f64>)> {
    let mut out =
        accumulate_gene_sum_multi(data_vec, &[(labels_a, k_a), (labels_b, k_b)], m, block_size)?;
    let sum_b = out.pop().expect("two groupings in, two out");
    let sum_a = out.pop().expect("two groupings in, two out");
    Ok((sum_a, sum_b))
}

fn accumulate_gene_sum_multi(
    data_vec: &SparseIoVec,
    groupings: &[(&[usize], usize)],
    m: usize,
    block_size: usize,
) -> anyhow::Result<Vec<Vec<f64>>> {
    let n = groupings[0].0.len();
    for (labels, _) in groupings {
        anyhow::ensure!(
            labels.len() == n,
            "grouping length mismatch: {} vs {n}",
            labels.len()
        );
    }
    let blocks: Vec<(usize, usize)> = (0..n)
        .step_by(block_size.max(1))
        .map(|lb| (lb, (lb + block_size).min(n)))
        .collect();

    let acc_bytes = groupings.iter().map(|&(_, k)| k * m).sum::<usize>() * size_of::<f64>();
    let affordable = (ACCUMULATOR_BUDGET_BYTES / acc_bytes.max(1)).max(1);
    let workers = rayon::current_num_threads()
        .min(affordable)
        .min(blocks.len().max(1));
    if workers < rayon::current_num_threads() {
        log::info!(
            "gene-sum aggregation: {workers} workers (not {}) — each holds {} MiB of group sums",
            rayon::current_num_threads(),
            acc_bytes >> 20
        );
    }
    let zeros = || -> Vec<Vec<f64>> {
        groupings
            .iter()
            .map(|&(_, k)| vec![0.0f64; k * m])
            .collect()
    };

    let chunk = blocks.len().div_ceil(workers.max(1)).max(1);
    blocks
        .par_chunks(chunk)
        .map(|chunk| {
            let mut acc = zeros();
            for &(lb, ub) in chunk {
                add_block(&mut acc, data_vec, groupings, lb, ub, m)?;
            }
            anyhow::Ok(acc)
        })
        .try_reduce_with(|mut a, b| {
            for (x, y) in a.iter_mut().zip(b) {
                for (xi, yi) in x.iter_mut().zip(y) {
                    *xi += yi;
                }
            }
            anyhow::Ok(a)
        })
        .unwrap_or_else(|| Ok(zeros()))
}

fn add_block(
    sums: &mut [Vec<f64>],
    data_vec: &SparseIoVec,
    groupings: &[(&[usize], usize)],
    lb: usize,
    ub: usize,
    m: usize,
) -> anyhow::Result<()> {
    let csc = data_vec.read_columns_csc(lb..ub)?;
    for j in 0..csc.ncols() {
        let col = csc.col(j);
        for (sum, &(labels, k)) in sums.iter_mut().zip(groupings) {
            let kk = labels[lb + j];
            if kk >= k {
                continue;
            }
            let row = &mut sum[kk * m..(kk + 1) * m];
            for (&g, &v) in col.row_indices().iter().zip(col.values()) {
                row[g] += f64::from(v);
            }
        }
    }
    Ok(())
}

/// Row-major `k · m` gene sums → `m × k` mean profiles with optional per-gene
/// weights: `μ[g, c] = w[g] · sum[c, g] / Σ_g sum[c, g]`.
pub fn weighted_mean_profile(
    gene_sum: &[f64],
    n_groups: usize,
    n_genes: usize,
    weights: &[f32],
) -> Mat {
    debug_assert_eq!(gene_sum.len(), n_groups * n_genes);
    debug_assert!(weights.is_empty() || weights.len() == n_genes);
    let mut out = Mat::zeros(n_genes, n_groups);
    out.as_mut_slice()
        .par_chunks_mut(n_genes)
        .enumerate()
        .for_each(|(c, col)| {
            let row = &gene_sum[c * n_genes..(c + 1) * n_genes];
            let inv = (1.0 / row.iter().sum::<f64>().max(1.0)) as f32;
            for (gi, slot) in col.iter_mut().enumerate() {
                let w = weights.get(gi).copied().unwrap_or(1.0);
                *slot = row[gi] as f32 * inv * w;
            }
        });
    out
}
