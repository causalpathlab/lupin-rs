//! Gene sets as units on a run's gene table: each cell type (its markers)
//! and ontology term (its annotated genes) is trained by SGD as a node of a
//! gene–set graph whose gene rows are the run's own, frozen, so the sets are
//! placed in the space the cells share with the genes and cannot bend it.
//! The engine is PyTorch-BigGraph's (Lerer et al. 2019), as
//! `graph_embedding_util::fne` implements it; gene rows are pinned with its
//! preset rows, as `senna fne --freeze-feature-embedding` does.

use crate::annotate::gene_rows::GeneRows;
use crate::annotate::markers::label_key;
use crate::annotate::outputs::{labels_path, AnnotationOutputs};
use anyhow::Result;
use graph_embedding_util::fne::{
    train, FneConfig, NodeTypeTable, PresetMode, PresetRows, Relation, RelationPolarity,
    RelationTable, TypedEdgeList,
};
use legume_numeric::matrix::dense_mat_io::{Mat, MatWithNames};
use legume_numeric::matrix::parquet::{write_named_table, Column};
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet};

/// Each cell type's genes in the table and given, radius, share of the
/// cells and centre (`h0…`), by type.
pub const UNITS_PARQUET: &str = ".celltype_units.parquet";
/// Each pair of cell types and the cells they share (`Regions::links`).
pub const TYPE_LINKS_TSV: &str = ".celltype_links.tsv";
/// Each cluster against each ontology term: score, p and q ([`align`]).
pub const TERM_ALIGNMENT_PARQUET: &str = ".term_alignment.parquet";
/// Each ontology term's genes in the table and given, and centre (`h0…`).
pub const TERM_UNITS_PARQUET: &str = ".term_units.parquet";
/// What a pass writes beside the per-cell labels, cleared before it runs so
/// no file of an earlier pass (its terms, say) is read as this one's.
pub const UNITS_OUTPUT_SUFFIXES: [&str; 4] = [
    UNITS_PARQUET,
    TYPE_LINKS_TSV,
    TERM_ALIGNMENT_PARQUET,
    TERM_UNITS_PARQUET,
];

/// A run's genes and their rows (row-major, `dim` wide), in the space its
/// cells share with them.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneTable {
    pub names: Vec<String>,
    pub rows: Vec<f32>,
    pub dim: usize,
}

/// What a gene set stands for: a cell type (its markers) or an ontology
/// term (its annotated genes). Each kind is its own node type, so a set is
/// contrasted only with sets of its kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SetKind {
    CellType,
    Term,
}

impl SetKind {
    const ALL: [Self; 2] = [Self::CellType, Self::Term];

    fn name(self) -> &'static str {
        match self {
            Self::CellType => "cell_type",
            Self::Term => "term",
        }
    }
}

/// Ontology terms for a pass, and the genes they were annotated from: the
/// pool their decoys are drawn from.
#[derive(Debug, Clone, Default)]
pub struct Terms {
    pub sets: Vec<GeneSet>,
    pub pool: Vec<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GeneSet {
    pub name: String,
    pub kind: SetKind,
    pub genes: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct UnitConfig {
    pub epochs: usize,
    /// The fewest of a set's genes the table must have to train it.
    pub min_genes: usize,
    pub seed: u64,
    /// EM rounds fitting the regions to the cells.
    pub region_iters: usize,
    /// Random gene sets per size bin of the terms, trained as their null.
    pub decoys: usize,
    /// Where training runs.
    pub device: legume_numeric::candle::candle_core::Device,
}

impl Default for UnitConfig {
    fn default() -> Self {
        Self {
            epochs: 50,
            min_genes: 2,
            seed: 1,
            region_iters: 100,
            decoys: 0,
            device: legume_numeric::candle::candle_core::Device::Cpu,
        }
    }
}

/// The trained sets: each one's centre, and the genes the table had of it.
#[derive(Debug, Clone)]
pub struct Units {
    pub names: Vec<String>,
    pub kinds: Vec<SetKind>,
    /// Row-major, `dim` wide, one row per name.
    pub centres: Vec<f32>,
    pub dim: usize,
    /// Per set, the genes in the table and the genes given.
    pub covered: Vec<(usize, usize)>,
    /// Per set, the size its decoys have ([`size_bin`]).
    pub bins: Vec<usize>,
    /// The sets with too few genes in the table: name, found, given.
    pub dropped: Vec<(String, usize, usize)>,
    /// Random gene sets trained beside the sets, as their null.
    pub decoys: Vec<Decoy>,
}

/// A random gene set trained as a term is: its gene rows and centre.
#[derive(Debug, Clone)]
pub struct Decoy {
    /// The size bin it stands in for (its genes number fewer when the pool
    /// is smaller than the bin).
    pub bin: usize,
    /// Its gene rows, kept for the tests that check where they came from.
    #[cfg(test)]
    pub genes: Vec<u32>,
    pub centre: Vec<f32>,
}

impl Units {
    /// Set `i`'s centre.
    #[must_use]
    pub fn centre(&self, i: usize) -> &[f32] {
        &self.centres[i * self.dim..(i + 1) * self.dim]
    }

    /// The sets of `kind`, by index.
    fn of_kind(&self, kind: SetKind) -> Vec<usize> {
        (0..self.names.len())
            .filter(|&i| self.kinds[i] == kind)
            .collect()
    }

    #[cfg(test)]
    pub fn coverage(&self, name: &str) -> Option<(usize, usize)> {
        let i = self.names.iter().position(|n| n == name)?;
        Some(self.covered[i])
    }
}

/// Train `sets` against `genes`, whose rows stay as they are. A set's genes
/// are found as the other passes find markers ([`GeneRows`]). The terms get
/// `cfg.decoys` random gene sets per size bin ([`size_bin`]), trained beside
/// them as nodes of their type: the null [`align`] reads. Decoys are drawn
/// from `term_pool`, the genes terms could be annotated from (real sets are
/// not random genes, so a null drawn from every gene sits apart from all of
/// them), or from the whole table when it is empty.
pub fn train_units(
    genes: &GeneTable,
    sets: &[GeneSet],
    term_pool: &[String],
    cfg: &UnitConfig,
) -> Result<Units> {
    use rand::rngs::SmallRng;
    use rand::SeedableRng;
    let names: Vec<Box<str>> = genes.names.iter().map(|g| g.as_str().into()).collect();
    let index = GeneRows::build(&names);
    let rows_of = |given: &[String]| -> Vec<u32> {
        let mut rows: Vec<u32> = (given.iter())
            .filter_map(|g| index.match_rows(g))
            .flatten()
            .map(|&r| r as u32)
            .collect();
        rows.sort_unstable();
        rows.dedup();
        rows
    };
    let n_genes = genes.names.len();
    let mut dropped = Vec::new();
    let mut rng = SmallRng::seed_from_u64(cfg.seed);
    // Every node to train, kind by kind in node-type order: the kept sets,
    // then the terms' decoys (`None`), each with its size bin.
    let mut nodes: Vec<(SetKind, Vec<u32>, Option<&GeneSet>, usize)> = Vec::new();
    for kind in SetKind::ALL {
        let mut bins = std::collections::BTreeSet::new();
        for s in sets.iter().filter(|s| s.kind == kind) {
            let found = rows_of(&s.genes);
            if found.len() < cfg.min_genes {
                dropped.push((s.name.clone(), found.len(), s.genes.len()));
            } else {
                let bin = size_bin(found.len()).min(n_genes);
                bins.insert(bin);
                nodes.push((kind, found, Some(s), bin));
            }
        }
        if kind != SetKind::Term {
            continue;
        }
        let pool: Vec<u32> = if term_pool.is_empty() {
            (0..n_genes as u32).collect()
        } else {
            rows_of(term_pool)
        };
        for bin in bins {
            for _ in 0..cfg.decoys {
                let pick = rand::seq::index::sample(&mut rng, pool.len(), bin.min(pool.len()));
                let rows: Vec<u32> = pick.iter().map(|i| pool[i]).collect();
                nodes.push((kind, rows, None, bin));
            }
        }
    }
    anyhow::ensure!(
        nodes.iter().any(|n| n.2.is_some()),
        "no gene set has {} genes in the table",
        cfg.min_genes
    );
    // Checked before training: the regions need a cell type to fit.
    anyhow::ensure!(
        !sets.iter().any(|s| s.kind == SetKind::CellType)
            || nodes.iter().any(|n| n.0 == SetKind::CellType),
        "no cell type has {} of its markers in the gene table",
        cfg.min_genes
    );

    let present: Vec<SetKind> = (SetKind::ALL.into_iter())
        .filter(|k| nodes.iter().any(|n| n.0 == *k))
        .collect();
    let mut table = vec![("gene", n_genes)];
    table.extend((present.iter()).map(|k| (k.name(), nodes.iter().filter(|n| n.0 == *k).count())));
    let types = NodeTypeTable::new(&table)?;
    let rels = RelationTable::new(
        (present.iter().enumerate())
            .map(|(i, k)| Relation {
                name: format!("gene:{}", k.name()).into(),
                lhs_type: 0,
                rhs_type: (i + 1) as u16,
                weight: 1.0,
                undirected: false,
                polarity: RelationPolarity::Friend,
            })
            .collect(),
        &types,
    )?;
    let mut edges = TypedEdgeList::default();
    for (id, (kind, rows, ..)) in nodes.iter().enumerate() {
        let rel = present.iter().position(|k| k == kind).unwrap_or_default();
        for &g in rows {
            edges.lhs.push(g);
            edges.rhs.push((n_genes + id) as u32);
            edges.rel.push(rel as u16);
        }
    }
    let config = FneConfig {
        dim: genes.dim,
        epochs: cfg.epochs,
        num_batch_negs: 50.min(edges.lhs.len()),
        eval_fraction: 0.0,
        preset: Some(PresetRows {
            ids: (0..n_genes as u32).collect(),
            rows: genes.rows.clone(),
            mode: PresetMode::Freeze,
        }),
        seed: cfg.seed,
        device: cfg.device.clone(),
        ..FneConfig::default()
    };
    let out = train(edges, types, rels, &config)?;
    let all: Vec<f32> = out.embedding.flatten_all()?.to_vec1()?;
    let (gene_rows, trained) = all.split_at(n_genes * genes.dim);
    anyhow::ensure!(
        gene_rows == genes.rows.as_slice(),
        "training moved the gene table it was to keep fixed"
    );
    let row = |id: usize| trained[id * genes.dim..(id + 1) * genes.dim].to_vec();
    let mut units = Units {
        names: Vec::new(),
        kinds: Vec::new(),
        centres: Vec::new(),
        dim: genes.dim,
        covered: Vec::new(),
        bins: Vec::new(),
        dropped,
        decoys: Vec::new(),
    };
    for (id, (kind, rows, set, bin)) in nodes.into_iter().enumerate() {
        match set {
            Some(s) => {
                units.names.push(s.name.clone());
                units.kinds.push(kind);
                units.centres.extend(row(id));
                units.covered.push((rows.len(), s.genes.len()));
                units.bins.push(bin);
            }
            None => units.decoys.push(Decoy {
                bin,
                #[cfg(test)]
                genes: rows,
                centre: row(id),
            }),
        }
    }
    Ok(units)
}

/// The size a set of `n` genes is compared at: the nearest of 2, 3, 4, 6,
/// 8, 11, 16, … (`round(2^(i/2))`), so decoys are drawn for a few sizes, not
/// every one.
#[must_use]
pub fn size_bin(n: usize) -> usize {
    let i = (2.0 * (n.max(2) as f64).log2()).round();
    2f64.powf(i / 2.0).round() as usize
}

/// The sets of one kind as regions of the cells' space: each a ball around
/// its centre with a radius fitted to the cells, the cells' membership of
/// each, and how strongly each pair of regions shares cells.
#[derive(Debug, Clone)]
pub struct Regions {
    pub names: Vec<String>,
    /// The root mean square distance of a region's cells from its centre,
    /// on the unit sphere.
    pub radius: Vec<f32>,
    /// The share of the cells in each region.
    pub weight: Vec<f32>,
    /// `[n_cells × K]` row-major: each cell's membership, summing to 1.
    pub membership: Vec<f32>,
    /// `[K × K]` row-major, symmetric: the cells two regions share,
    /// `Σ_c w_ck w_cl / √(n_k n_l)` with `n_k = Σ_c w_ck`.
    pub links: Vec<f32>,
}

/// Fit the `kind` sets of `units` to `cells` (`[n × dim]` row-major, in the
/// genes' space as [`whitened`] puts both: trained by contrast with other
/// genes, a centre is a direction from the genes' mean) as a
/// mixture of isotropic balls on the unit sphere, with their centres fixed
/// where training put them: the radii and the regions'
/// shares by maximum likelihood (EM, exact for fixed centres), `iters`
/// rounds.
pub fn fit_regions(
    units: &Units,
    kind: SetKind,
    cells: &[f32],
    n: usize,
    iters: usize,
) -> Result<Regions> {
    let dim = units.dim;
    anyhow::ensure!(cells.len() == n * dim, "cells are not {n} rows of {dim}");
    let picked = units.of_kind(kind);
    let k = picked.len();
    anyhow::ensure!(k > 0, "no {} sets to fit", kind.name());
    let centres: Vec<f32> = picked.iter().flat_map(|&i| unit(units.centre(i))).collect();
    // Squared chord distances, cell by region: 2 − 2 u·m between unit
    // vectors (a zero cell is at chord² 1 from every centre).
    let dist: Vec<f32> = cells
        .par_chunks(dim)
        .flat_map_iter(|cell| {
            let norm = cell.iter().map(|x| x * x).sum::<f32>().sqrt();
            let scale = if norm > 0.0 { 1.0 / norm } else { 0.0 };
            let zero = norm == 0.0;
            centres.chunks(dim).map(move |m| {
                if zero {
                    1.0
                } else {
                    let dot: f32 = cell.iter().zip(m).map(|(a, b)| a * b).sum();
                    (2.0 - 2.0 * dot * scale).max(0.0)
                }
            })
        })
        .collect();
    // The sphere's own dimension.
    let d_eff = (dim.max(2) - 1) as f32;
    let mut var = vec![dist.iter().sum::<f32>() / (dist.len() as f32 * d_eff); k];
    let mut weight = vec![1.0 / k as f32; k];
    let mut w = vec![0.0f32; n * k];
    for _ in 0..iters.max(1) {
        // Each region's log prior and normaliser, and the scale of its
        // distances, once per round.
        let base: Vec<f32> = (0..k)
            .map(|j| weight[j].ln() - 0.5 * d_eff * var[j].ln())
            .collect();
        let inv: Vec<f32> = var.iter().map(|v| 0.5 / v).collect();
        let (mass, spread) = w
            .par_chunks_mut(k)
            .zip(dist.par_chunks(k))
            .fold(
                || (vec![0.0f32; k], vec![0.0f32; k]),
                |(mut mass, mut spread), (row, d)| {
                    for j in 0..k {
                        row[j] = base[j] - d[j] * inv[j];
                    }
                    let top = row.iter().copied().fold(f32::NEG_INFINITY, f32::max);
                    let mut z = 0.0;
                    for x in row.iter_mut() {
                        *x = (*x - top).exp();
                        z += *x;
                    }
                    for j in 0..k {
                        row[j] /= z;
                        mass[j] += row[j];
                        spread[j] += row[j] * d[j];
                    }
                    (mass, spread)
                },
            )
            .reduce(|| (vec![0.0f32; k], vec![0.0f32; k]), add_pair);
        for j in 0..k {
            weight[j] = (mass[j] / n as f32).max(f32::MIN_POSITIVE);
            var[j] = (spread[j] / (d_eff * mass[j].max(f32::MIN_POSITIVE))).max(1e-8);
        }
    }
    // Each region's mass and Wᵀ W, in one pass over the cells.
    let (mass, shared) = w
        .par_chunks(k)
        .fold(
            || (vec![0.0f32; k], vec![0.0f32; k * k]),
            |(mut mass, mut shared), row| {
                for a in 0..k {
                    mass[a] += row[a];
                    for b in a..k {
                        shared[a * k + b] += row[a] * row[b];
                    }
                }
                (mass, shared)
            },
        )
        .reduce(|| (vec![0.0f32; k], vec![0.0f32; k * k]), add_pair);
    let mut links = vec![0.0f32; k * k];
    for a in 0..k {
        for b in a..k {
            let v = shared[a * k + b] / (mass[a] * mass[b]).sqrt().max(f32::MIN_POSITIVE);
            links[a * k + b] = v;
            links[b * k + a] = v;
        }
    }
    Ok(Regions {
        names: picked.iter().map(|&i| units.names[i].clone()).collect(),
        radius: var.iter().map(|v| (v * d_eff).sqrt()).collect(),
        weight,
        membership: w,
        links,
    })
}

/// `lupin annotate --projection units` on loaded tables: `args`' panel and
/// ontology terms trained on the gene table `feat`, its outputs under
/// `args.out` cleared first unless `--no-clean` ([`run_pass`]).
pub fn run(
    args: &crate::annotate::args::UnitsArgs,
    feat: &MatWithNames<Mat>,
    cells: &MatWithNames<Mat>,
    groups: &[usize],
) -> Result<AnnotationOutputs> {
    if !args.no_clean {
        crate::annotate::outputs::clean_outputs(
            &args.out,
            graph_embedding_util::type_annotation::CELL_LABEL_FILE_SUFFIXES,
        );
        crate::annotate::outputs::clean_outputs(&args.out, &UNITS_OUTPUT_SUFFIXES);
    }
    let terms = args.terms(&feat.rows)?;
    run_pass(
        &args.markers,
        &args.out,
        feat,
        cells,
        groups,
        &terms,
        &args.config()?,
    )
}

/// A projection pass by gene-set units: the panel's cell types, and the
/// ontology `terms` when given, trained on the gene table `feat` (frozen).
/// The cell types are fitted to `cells` as regions, each cell called by the
/// region it most belongs to (its membership the probability), and written
/// with their links beside the labels; the terms are aligned with the
/// clusters of `groups` against their decoys.
pub fn run_pass(
    markers: &str,
    out: &str,
    feat: &MatWithNames<Mat>,
    cells: &MatWithNames<Mat>,
    groups: &[usize],
    terms: &Terms,
    cfg: &UnitConfig,
) -> Result<AnnotationOutputs> {
    legume_numeric::matrix::common_io::mkdir_parent(out)?;
    let mut panel: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for (gene, t) in crate::annotate::markers::read_panel(markers)? {
        panel.entry(label_key(&t)).or_default().insert(gene);
    }
    let mut sets: Vec<GeneSet> = (panel.into_iter())
        .map(|(name, genes)| GeneSet {
            name,
            kind: SetKind::CellType,
            genes: genes.into_iter().collect(),
        })
        .collect();
    sets.extend_from_slice(&terms.sets);
    // Both sides read per unit of the cells' spread, so an axis every type
    // spreads along does not outweigh the ones that tell types apart.
    let (gene_rows, flat) = whitened(feat, cells, groups)?;
    let units = train_units(&gene_rows, &sets, &terms.pool, cfg)?;
    for (name, found, given) in &units.dropped {
        log::warn!("{name}: {found} of its {given} markers are in the gene table; not trained");
    }
    let n = cells.rows.len();
    let r = fit_regions(&units, SetKind::CellType, &flat, n, cfg.region_iters)?;
    let k = r.names.len();
    let (labels, probs): (Vec<Box<str>>, Vec<f32>) = (0..n)
        .map(|c| {
            let row = &r.membership[c * k..(c + 1) * k];
            let (best, p) =
                (row.iter().enumerate())
                    .fold((0, f32::MIN), |b, (j, &p)| if p > b.1 { (j, p) } else { b });
            (r.names[best].as_str().into(), p)
        })
        .unzip();
    graph_embedding_util::type_annotation::write_cell_labels(out, &cells.rows, &labels, &probs)?;
    write_sets(
        &format!("{out}{UNITS_PARQUET}"),
        &units,
        &units.of_kind(SetKind::CellType),
        &[("radius", &r.radius), ("weight", &r.weight)],
    )?;
    write_links(&format!("{out}{TYPE_LINKS_TSV}"), &r)?;
    let (term_alignment, term_units) = if units.kinds.contains(&SetKind::Term) {
        let path = format!("{out}{TERM_ALIGNMENT_PARQUET}");
        let a = align(&units, SetKind::Term, &flat, groups)?;
        write_alignment(&path, &a)?;
        let sets = format!("{out}{TERM_UNITS_PARQUET}");
        write_sets(&sets, &units, &units.of_kind(SetKind::Term), &[])?;
        (Some(path), Some(sets))
    } else {
        (None, None)
    };
    Ok(AnnotationOutputs {
        argmax: Some(labels_path(out)),
        celltype_units: Some(format!("{out}{UNITS_PARQUET}")),
        celltype_links: Some(format!("{out}{TYPE_LINKS_TSV}")),
        term_alignment,
        term_units,
        ..AnnotationOutputs::default()
    })
}

/// The genes (from their mean) and the cells (from the balanced centre of
/// `groups`, each cluster's mean weighed alike) of `feat` and `cells`, both
/// whitened by the cells' covariance within their clusters (`Σ_w^{-1/2}`,
/// as linear discriminant analysis reads them; axes the cells barely spread
/// along are floored at a millionth of the widest): spread every cluster shares then
/// weighs no more than the axes that tell clusters apart. Genes as a [`GeneTable`], cells row-major.
fn whitened(
    feat: &MatWithNames<Mat>,
    cells: &MatWithNames<Mat>,
    groups: &[usize],
) -> Result<(GeneTable, Vec<f32>)> {
    use nalgebra::{DMatrix, DVector};
    let (n, dim) = cells.mat.shape();
    anyhow::ensure!(feat.mat.ncols() == dim, "genes and cells differ in width");
    anyhow::ensure!(groups.len() == n, "a cluster per cell");
    let row = |c: usize| {
        DVector::<f64>::from_iterator(dim, cells.mat.row(c).iter().map(|&v| f64::from(v)))
    };
    let mut by: BTreeMap<usize, (DVector<f64>, usize)> = BTreeMap::new();
    for (c, &g) in groups.iter().enumerate().filter(|(_, g)| **g != usize::MAX) {
        let x = row(c);
        anyhow::ensure!(
            x.iter().all(|v| v.is_finite()),
            "cell {c} has a non-finite coordinate"
        );
        let e = by.entry(g).or_insert_with(|| (DVector::zeros(dim), 0));
        e.0 += x;
        e.1 += 1;
    }
    anyhow::ensure!(!by.is_empty(), "no cell is in a cluster");
    let means: BTreeMap<usize, DVector<f64>> = by
        .into_iter()
        .map(|(g, (s, k))| (g, s / k as f64))
        .collect();
    // The balanced centre: each cluster's mean weighed alike.
    let centre = means.values().sum::<DVector<f64>>() / means.len() as f64;
    // Pooled within the clusters: the spread that tells no cluster apart,
    // summed cell by cell into one dim × dim matrix.
    let mut within = DMatrix::<f64>::zeros(dim, dim);
    let mut d = DVector::<f64>::zeros(dim);
    let mut counted = 0usize;
    for (c, g) in groups.iter().enumerate() {
        if let Some(m) = means.get(g) {
            d.copy_from(&row(c));
            d -= m;
            within.ger(1.0, &d, &d, 1.0);
            counted += 1;
        }
    }
    let cov = within / counted.saturating_sub(means.len()).max(1) as f64;
    let eig = cov.symmetric_eigen();
    let top = eig.eigenvalues.max();
    anyhow::ensure!(
        top.is_finite() && top > 0.0,
        "the cells do not spread within their clusters: nothing to read the genes by"
    );
    // Axes the cells barely spread along are floored at a millionth of the
    // widest, so they are scaled up at most a thousandfold.
    let floor = top * 1e-6;
    let inv_sqrt = DVector::from_iterator(
        dim,
        eig.eigenvalues.iter().map(|&l| 1.0 / l.max(floor).sqrt()),
    );
    let w = &eig.eigenvectors * DMatrix::from_diagonal(&inv_sqrt) * eig.eigenvectors.transpose();
    // Each row from `origin`, through `w`, row-major f32: one row at a time.
    let project = |m: &Mat, origin: &DVector<f64>| -> Vec<f32> {
        let mut out = Vec::with_capacity(m.nrows() * dim);
        let mut x = DVector::<f64>::zeros(dim);
        for r in 0..m.nrows() {
            for (j, v) in m.row(r).iter().enumerate() {
                x[j] = f64::from(*v) - origin[j];
            }
            out.extend(w.tr_mul(&x).iter().map(|&v| v as f32));
        }
        out
    };
    let gene_mean = DVector::<f64>::from_iterator(
        dim,
        (0..dim).map(|j| {
            feat.mat
                .column(j)
                .iter()
                .map(|&v| f64::from(v))
                .sum::<f64>()
                / feat.mat.nrows().max(1) as f64
        }),
    );
    Ok((
        GeneTable {
            names: feat.rows.iter().map(ToString::to_string).collect(),
            rows: project(&feat.mat, &gene_mean),
            dim,
        },
        project(&cells.mat, &centre),
    ))
}

/// The sets `at` of `units` as a table: genes found and given, the `extra`
/// columns, then the centre (`h0…`), one row per set.
fn write_sets(path: &str, units: &Units, at: &[usize], extra: &[(&str, &[f32])]) -> Result<()> {
    let found: Vec<f32> = at.iter().map(|&i| units.covered[i].0 as f32).collect();
    let given: Vec<f32> = at.iter().map(|&i| units.covered[i].1 as f32).collect();
    let centre: Vec<Vec<f32>> = (0..units.dim)
        .map(|d| (at.iter()).map(|&i| units.centre(i)[d]).collect())
        .collect();
    let dims: Vec<String> = (0..units.dim).map(|d| format!("h{d}")).collect();
    let mut cols = vec![
        ("genes_found".into(), Column::F32(&found)),
        ("genes_given".into(), Column::F32(&given)),
    ];
    cols.extend(extra.iter().map(|(n, v)| ((*n).into(), Column::F32(v))));
    cols.extend((dims.iter().zip(&centre)).map(|(n, v)| (n.as_str().into(), Column::F32(v))));
    let rows: Vec<Box<str>> = at.iter().map(|&i| units.names[i].as_str().into()).collect();
    write_named_table(path, "set", &rows, &cols)?;
    Ok(())
}

/// `a` as a long table: cluster (`K{id}`), term, score, p, q.
fn write_alignment(path: &str, a: &Alignment) -> Result<()> {
    let t = a.terms.len();
    let cluster: Vec<Box<str>> = (a.clusters.iter())
        .flat_map(|k| std::iter::repeat_n(format!("K{k}").into(), t))
        .collect();
    let term: Vec<Box<str>> = (0..a.clusters.len())
        .flat_map(|_| a.terms.iter().map(|x| x.as_str().into()))
        .collect();
    let rows: Vec<Box<str>> = (0..cluster.len()).map(|i| i.to_string().into()).collect();
    write_named_table(
        path,
        "row",
        &rows,
        &[
            ("cluster".into(), Column::Str(&cluster)),
            ("term".into(), Column::Str(&term)),
            ("score".into(), Column::F32(&a.score)),
            ("p".into(), Column::F32(&a.p)),
            ("q".into(), Column::F32(&a.q)),
        ],
    )?;
    Ok(())
}

fn write_links(path: &str, r: &Regions) -> Result<()> {
    let k = r.names.len();
    let mut text = String::from("type_a\ttype_b\tlink\n");
    for a in 0..k {
        for b in (a + 1)..k {
            text += &format!("{}\t{}\t{}\n", r.names[a], r.names[b], r.links[a * k + b]);
        }
    }
    std::fs::write(path, text)?;
    Ok(())
}

/// How each cluster's cells align with each set of one kind: `[clusters ×
/// sets]` row-major.
#[derive(Debug, Clone)]
pub struct Alignment {
    pub terms: Vec<String>,
    pub clusters: Vec<usize>,
    /// The mean cosine of a cluster's cells (as [`whitened`] puts them) to the
    /// set's centre.
    pub score: Vec<f32>,
    /// The normal upper tail of the score among its decoys' scores (their
    /// mean and spread), so it is not floored at one over their count; NaN
    /// without two decoys.
    pub p: Vec<f32>,
    /// `p` adjusted by Benjamini–Hochberg within each cluster.
    pub q: Vec<f32>,
}

/// Align the clusters of `cells` (`[n × dim]`, `groups` per cell,
/// `usize::MAX` left out) with the `kind` sets of `units`, each against the
/// decoys of its size.
pub fn align(units: &Units, kind: SetKind, cells: &[f32], groups: &[usize]) -> Result<Alignment> {
    let dim = units.dim;
    let n = groups.len();
    anyhow::ensure!(cells.len() == n * dim, "cells are not {n} rows of {dim}");
    let clusters: Vec<usize> = {
        let set: std::collections::BTreeSet<usize> = groups
            .iter()
            .copied()
            .filter(|&g| g != usize::MAX)
            .collect();
        set.into_iter().collect()
    };
    // Each cluster's mean unit vector, in one pass over the cells.
    let at: BTreeMap<usize, usize> = clusters.iter().enumerate().map(|(i, &k)| (k, i)).collect();
    let mut means = vec![vec![0.0f32; dim]; clusters.len()];
    let mut counts = vec![0usize; clusters.len()];
    for (cell, g) in cells.chunks(dim).zip(groups) {
        if let Some(&i) = at.get(g) {
            means[i]
                .iter_mut()
                .zip(unit(cell))
                .for_each(|(a, b)| *a += b);
            counts[i] += 1;
        }
    }
    for (m, &c) in means.iter_mut().zip(&counts) {
        m.iter_mut().for_each(|a| *a /= c.max(1) as f32);
    }
    let dot = |a: &[f32], b: &[f32]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f32>();
    let picked = units.of_kind(kind);
    let dirs: Vec<Vec<f32>> = picked.iter().map(|&i| unit(units.centre(i))).collect();
    let decoy_dirs: Vec<(usize, Vec<f32>)> = (units.decoys.iter())
        .map(|d| (d.bin, unit(&d.centre)))
        .collect();
    let t = picked.len();
    let mut score = Vec::with_capacity(clusters.len() * t);
    let mut p = Vec::with_capacity(clusters.len() * t);
    let mut q = Vec::with_capacity(clusters.len() * t);
    for m in &means {
        // The cluster's null per decoy size: its scores' mean and spread.
        let mut by_bin: BTreeMap<usize, Vec<f32>> = BTreeMap::new();
        for (bin, dir) in &decoy_dirs {
            by_bin.entry(*bin).or_default().push(dot(m, dir));
        }
        let nulls: BTreeMap<usize, Null> =
            by_bin.into_iter().map(|(b, v)| (b, Null::of(&v))).collect();
        let row_p: Vec<f32> = (picked.iter().zip(&dirs))
            .map(|(&i, dir)| {
                let s = dot(m, dir);
                score.push(s);
                nulls
                    .get(&units.bins[i])
                    .map_or(f32::NAN, |n| n.upper_tail(s))
            })
            .collect();
        q.extend(bh_where_finite(&row_p));
        p.extend(row_p);
    }
    Ok(Alignment {
        terms: picked.iter().map(|&i| units.names[i].clone()).collect(),
        clusters,
        score,
        p,
        q,
    })
}

/// Benjamini–Hochberg over the finite `p`; a NaN p (a term without a null)
/// stays NaN and is not counted.
pub(crate) fn bh_where_finite(p: &[f32]) -> Vec<f32> {
    let at: Vec<usize> = (0..p.len()).filter(|&i| p[i].is_finite()).collect();
    let finite: Vec<f32> = at.iter().map(|&i| p[i]).collect();
    let adjusted = legume_numeric::matrix::hypothesis::benjamini_hochberg(&finite);
    let mut q = vec![f32::NAN; p.len()];
    for (&i, v) in at.iter().zip(adjusted) {
        q[i] = v;
    }
    q
}

/// Decoy scores as a normal: their mean and sample standard deviation;
/// `None` for the spread of fewer than two.
struct Null {
    mean: f64,
    sd: Option<f64>,
}

impl Null {
    fn of(scores: &[f32]) -> Self {
        let n = scores.len() as f64;
        let mean = scores.iter().map(|&x| f64::from(x)).sum::<f64>() / n.max(1.0);
        let sd = (n >= 2.0).then(|| {
            let ss: f64 = scores.iter().map(|&x| (f64::from(x) - mean).powi(2)).sum();
            (ss / (n - 1.0)).sqrt()
        });
        Self { mean, sd }
    }

    /// `P(Z ≥ z)`, `s` standardised by the null; NaN without a spread.
    fn upper_tail(&self, s: f32) -> f32 {
        use statrs::distribution::{ContinuousCDF, Normal};
        let Some(sd) = self.sd else {
            return f32::NAN;
        };
        let z = (f64::from(s) - self.mean) / sd.max(f64::MIN_POSITIVE);
        Normal::standard().sf(z) as f32
    }
}

/// Two pairs of sums, added.
fn add_pair(
    (mut a, mut b): (Vec<f32>, Vec<f32>),
    (c, d): (Vec<f32>, Vec<f32>),
) -> (Vec<f32>, Vec<f32>) {
    a.iter_mut().zip(c).for_each(|(x, y)| *x += y);
    b.iter_mut().zip(d).for_each(|(x, y)| *x += y);
    (a, b)
}

/// `v` scaled to unit length (zero stays zero).
fn unit(v: &[f32]) -> Vec<f32> {
    let n = v.iter().map(|x| x * x).sum::<f32>().sqrt();
    if n > 0.0 {
        v.iter().map(|x| x / n).collect()
    } else {
        v.to_vec()
    }
}

#[cfg(test)]
#[path = "tests/geneset_units.rs"]
mod tests;
