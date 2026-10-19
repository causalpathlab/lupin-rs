//! Gene sets as units on a run's gene table: each cell type (its markers)
//! and ontology term (its annotated genes) is trained by SGD as a node of a
//! gene–set graph whose gene rows are the run's own, frozen, so the sets are
//! placed in the space the cells share with the genes and cannot bend it.
//! The engine is PyTorch-BigGraph's (Lerer et al. 2019), as
//! `graph_embedding_util::fne` implements it; gene rows are pinned with its
//! preset rows, as `senna fne --freeze-feature-embedding` does.

use crate::annotate::markers::label_key;
use crate::annotate::outputs::{labels_path, AnnotationOutputs};
use anyhow::Result;
use graph_embedding_util::fne::{
    train, FneConfig, NodeTypeTable, PresetMode, PresetRows, Relation, RelationPolarity,
    RelationTable, TypedEdgeList,
};
use legume_numeric::matrix::dense_mat_io::{Mat, MatWithNames};
use legume_numeric::matrix::parquet::{write_named_table, Column};
use std::collections::HashMap;

/// Each cell type's genes in the table and given, radius, share of the
/// cells and centre (`h0…`), by type.
pub const UNITS_PARQUET: &str = ".celltype_units.parquet";
/// Each pair of cell types and the cells they share (`Regions::links`).
pub const TYPE_LINKS_TSV: &str = ".celltype_links.tsv";

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
            device: legume_numeric::candle::candle_core::Device::Cpu,
        }
    }
}

/// The trained sets: each one's centre, and the genes the table had of it.
#[derive(Debug, Clone)]
pub struct Units {
    /// The genes' mean row. Trained by contrast with other genes, a set's
    /// centre is a direction from here, not a point among its genes; cells
    /// are read from here too ([`fit_regions`]).
    pub origin: Vec<f32>,
    pub names: Vec<String>,
    pub kinds: Vec<SetKind>,
    /// Row-major, `dim` wide, one row per name.
    pub centres: Vec<f32>,
    pub dim: usize,
    /// Per set, the genes in the table and the genes given.
    pub covered: Vec<(usize, usize)>,
    /// The sets with too few genes in the table: name, found, given.
    pub dropped: Vec<(String, usize, usize)>,
}

impl Units {
    fn index(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }

    #[cfg(test)]
    pub fn centre(&self, name: &str) -> Option<&[f32]> {
        let i = self.index(name)?;
        Some(&self.centres[i * self.dim..(i + 1) * self.dim])
    }

    #[cfg(test)]
    pub fn coverage(&self, name: &str) -> Option<(usize, usize)> {
        self.index(name).map(|i| self.covered[i])
    }
}

/// Train `sets` against `genes`, whose rows stay as they are.
pub fn train_units(genes: &GeneTable, sets: &[GeneSet], cfg: &UnitConfig) -> Result<Units> {
    let row_of: HashMap<&str, u32> = (genes.names.iter().enumerate())
        .map(|(i, g)| (g.as_str(), i as u32))
        .collect();
    let mut dropped = Vec::new();
    // The sets kept, by kind in node-type order, each with its gene rows.
    let mut kept: Vec<(&GeneSet, Vec<u32>)> = Vec::new();
    for kind in SetKind::ALL {
        for s in sets.iter().filter(|s| s.kind == kind) {
            let mut found: Vec<u32> = (s.genes.iter())
                .filter_map(|g| row_of.get(g.as_str()).copied())
                .collect();
            found.sort_unstable();
            found.dedup();
            if found.len() < cfg.min_genes {
                dropped.push((s.name.clone(), found.len(), s.genes.len()));
            } else {
                kept.push((s, found));
            }
        }
    }
    anyhow::ensure!(
        !kept.is_empty(),
        "no gene set has {} genes in the table",
        cfg.min_genes
    );

    let n_genes = genes.names.len();
    let present: Vec<SetKind> = (SetKind::ALL.into_iter())
        .filter(|k| kept.iter().any(|(s, _)| s.kind == *k))
        .collect();
    let mut table = vec![("gene", n_genes)];
    table.extend(
        present
            .iter()
            .map(|k| (k.name(), kept.iter().filter(|(s, _)| s.kind == *k).count())),
    );
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
    for (set_id, (s, rows)) in kept.iter().enumerate() {
        let rel = present
            .iter()
            .position(|k| *k == s.kind)
            .unwrap_or_default();
        for &g in rows {
            edges.lhs.push(g);
            edges.rhs.push((n_genes + set_id) as u32);
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
    let (gene_rows, centres) = all.split_at(n_genes * genes.dim);
    let origin: Vec<f32> = (0..genes.dim)
        .map(|d| {
            (0..n_genes)
                .map(|g| genes.rows[g * genes.dim + d])
                .sum::<f32>()
                / n_genes as f32
        })
        .collect();
    anyhow::ensure!(
        gene_rows == genes.rows.as_slice(),
        "training moved the gene table it was to keep fixed"
    );
    Ok(Units {
        origin,
        names: kept.iter().map(|(s, _)| s.name.clone()).collect(),
        kinds: kept.iter().map(|(s, _)| s.kind).collect(),
        centres: centres.to_vec(),
        dim: genes.dim,
        covered: kept.iter().map(|(s, r)| (r.len(), s.genes.len())).collect(),
        dropped,
    })
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
/// genes' space, read from the genes' origin as the centres are) as a
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
    let picked: Vec<usize> = (0..units.names.len())
        .filter(|&i| units.kinds[i] == kind)
        .collect();
    let k = picked.len();
    anyhow::ensure!(k > 0, "no {} sets to fit", kind.name());
    let centres: Vec<Vec<f32>> = (picked.iter())
        .map(|&i| unit(&units.centres[i * dim..(i + 1) * dim]))
        .collect();
    // Squared chord distances, cell by region.
    let dist: Vec<f32> = (0..n)
        .flat_map(|c| {
            let from_origin: Vec<f32> = (cells[c * dim..(c + 1) * dim].iter())
                .zip(&units.origin)
                .map(|(x, o)| x - o)
                .collect();
            let u = unit(&from_origin);
            (centres.iter())
                .map(move |m| u.iter().zip(m).map(|(a, b)| (a - b) * (a - b)).sum())
                .collect::<Vec<f32>>()
        })
        .collect();
    // The sphere's own dimension.
    let d_eff = (dim.max(2) - 1) as f32;
    let mut var = vec![dist.iter().sum::<f32>() / (dist.len() as f32 * d_eff); k];
    let mut weight = vec![1.0 / k as f32; k];
    let mut w = vec![0.0f32; n * k];
    for _ in 0..iters.max(1) {
        for c in 0..n {
            let row = &mut w[c * k..(c + 1) * k];
            for j in 0..k {
                row[j] =
                    weight[j].ln() - 0.5 * d_eff * var[j].ln() - dist[c * k + j] / (2.0 * var[j]);
            }
            let top = row.iter().copied().fold(f32::NEG_INFINITY, f32::max);
            let z: f32 = row.iter().map(|x| (x - top).exp()).sum();
            row.iter_mut().for_each(|x| *x = (*x - top).exp() / z);
        }
        for j in 0..k {
            let mass: f32 = (0..n).map(|c| w[c * k + j]).sum();
            let spread: f32 = (0..n).map(|c| w[c * k + j] * dist[c * k + j]).sum();
            weight[j] = (mass / n as f32).max(f32::MIN_POSITIVE);
            var[j] = (spread / (d_eff * mass.max(f32::MIN_POSITIVE))).max(1e-8);
        }
    }
    let mass: Vec<f32> = (0..k).map(|j| (0..n).map(|c| w[c * k + j]).sum()).collect();
    let mut links = vec![0.0f32; k * k];
    for a in 0..k {
        for b in 0..k {
            let shared: f32 = (0..n).map(|c| w[c * k + a] * w[c * k + b]).sum();
            links[a * k + b] = shared / (mass[a] * mass[b]).sqrt().max(f32::MIN_POSITIVE);
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

/// A projection pass by gene-set units: the panel's cell types trained on
/// the gene table `feat` (frozen), fitted to `cells` as regions, each cell
/// called by the region it most belongs to (its membership the
/// probability); the units and their links written beside the labels.
pub fn run_pass(
    markers: &str,
    out: &str,
    feat: &MatWithNames<Mat>,
    cells: &MatWithNames<Mat>,
    cfg: &UnitConfig,
) -> Result<AnnotationOutputs> {
    legume_numeric::matrix::common_io::mkdir_parent(out)?;
    let mut by_type: Vec<GeneSet> = Vec::new();
    for (gene, t) in crate::annotate::markers::read_panel(markers)? {
        let name = label_key(&t);
        match by_type.iter_mut().find(|s| s.name == name) {
            Some(s) => s.genes.push(gene),
            None => by_type.push(GeneSet {
                name,
                kind: SetKind::CellType,
                genes: vec![gene],
            }),
        }
    }
    let units = train_units(&table_of(feat), &by_type, cfg)?;
    for (name, found, given) in &units.dropped {
        log::warn!("{name}: {found} of its {given} markers are in the gene table; not trained");
    }
    let n = cells.rows.len();
    let flat = table_of(cells).rows;
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
    write_units(&format!("{out}{UNITS_PARQUET}"), &units, &r)?;
    write_links(&format!("{out}{TYPE_LINKS_TSV}"), &r)?;
    Ok(AnnotationOutputs {
        argmax: Some(labels_path(out)),
        celltype_units: Some(format!("{out}{UNITS_PARQUET}")),
        celltype_links: Some(format!("{out}{TYPE_LINKS_TSV}")),
        ..AnnotationOutputs::default()
    })
}

/// `m`'s rows as a [`GeneTable`] (row-major).
fn table_of(m: &MatWithNames<Mat>) -> GeneTable {
    GeneTable {
        names: m.rows.iter().map(ToString::to_string).collect(),
        rows: (0..m.mat.nrows())
            .flat_map(|i| (0..m.mat.ncols()).map(move |j| (i, j)))
            .map(|(i, j)| m.mat[(i, j)])
            .collect(),
        dim: m.mat.ncols(),
    }
}

fn write_units(path: &str, units: &Units, r: &Regions) -> Result<()> {
    let at: Vec<usize> = r.names.iter().filter_map(|n| units.index(n)).collect();
    let found: Vec<f32> = at.iter().map(|&i| units.covered[i].0 as f32).collect();
    let given: Vec<f32> = at.iter().map(|&i| units.covered[i].1 as f32).collect();
    let centre: Vec<Vec<f32>> = (0..units.dim)
        .map(|d| {
            at.iter()
                .map(|&i| units.centres[i * units.dim + d])
                .collect()
        })
        .collect();
    let names: Vec<String> = (0..units.dim).map(|d| format!("h{d}")).collect();
    let mut cols = vec![
        ("genes_found".into(), Column::F32(&found)),
        ("genes_given".into(), Column::F32(&given)),
        ("radius".into(), Column::F32(&r.radius)),
        ("weight".into(), Column::F32(&r.weight)),
    ];
    cols.extend((names.iter().zip(&centre)).map(|(n, v)| (n.as_str().into(), Column::F32(v))));
    let rows: Vec<Box<str>> = r.names.iter().map(|n| n.as_str().into()).collect();
    write_named_table(path, "type", &rows, &cols)?;
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
