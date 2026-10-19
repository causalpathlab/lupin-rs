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
/// Each cluster against each ontology term: score, p and q ([`align`]).
pub const TERM_ALIGNMENT_PARQUET: &str = ".term_alignment.parquet";

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
    /// Random gene sets per kind and size bin, trained as the null.
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
    /// Per set, the size its decoys have ([`size_bin`]).
    pub bins: Vec<usize>,
    /// The sets with too few genes in the table: name, found, given.
    pub dropped: Vec<(String, usize, usize)>,
    /// Random gene sets trained beside the sets, as their null.
    pub decoys: Vec<Decoy>,
}

/// A random gene set trained as a `kind` set is: its gene rows and centre.
#[derive(Debug, Clone)]
pub struct Decoy {
    pub kind: SetKind,
    pub genes: Vec<u32>,
    pub centre: Vec<f32>,
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

/// Train `sets` against `genes`, whose rows stay as they are, with
/// `cfg.decoys` random gene sets per kind and size bin ([`size_bin`]) trained
/// beside them as nodes of the same type, the null [`align`] reads.
#[cfg(test)]
pub fn train_units(genes: &GeneTable, sets: &[GeneSet], cfg: &UnitConfig) -> Result<Units> {
    train_units_from(genes, sets, &[], cfg)
}

/// [`train_units`], each kind's decoys drawn from its `pools` of genes (the
/// genes its sets could be annotated from: real sets are not random genes,
/// so a null drawn from every gene sits apart from all of them); a kind
/// without a pool draws from the whole table.
pub fn train_units_from(
    genes: &GeneTable,
    sets: &[GeneSet],
    pools: &[(SetKind, Vec<String>)],
    cfg: &UnitConfig,
) -> Result<Units> {
    use rand::rngs::SmallRng;
    use rand::SeedableRng;
    let row_of: HashMap<&str, u32> = (genes.names.iter().enumerate())
        .map(|(i, g)| (g.as_str(), i as u32))
        .collect();
    let n_genes = genes.names.len();
    let mut dropped = Vec::new();
    let mut rng = SmallRng::seed_from_u64(cfg.seed);
    // Every node to train, kind by kind in node-type order: the kept sets,
    // then their decoys. `None` for a decoy.
    let mut nodes: Vec<(SetKind, Vec<u32>, Option<&GeneSet>)> = Vec::new();
    for kind in SetKind::ALL {
        let mut bins = std::collections::BTreeSet::new();
        for s in sets.iter().filter(|s| s.kind == kind) {
            let mut found: Vec<u32> = (s.genes.iter())
                .filter_map(|g| row_of.get(g.as_str()).copied())
                .collect();
            found.sort_unstable();
            found.dedup();
            if found.len() < cfg.min_genes {
                dropped.push((s.name.clone(), found.len(), s.genes.len()));
            } else {
                bins.insert(size_bin(found.len()).min(n_genes));
                nodes.push((kind, found, Some(s)));
            }
        }
        let pool: Vec<u32> = match pools.iter().find(|(k, _)| *k == kind) {
            Some((_, names)) => {
                let mut rows: Vec<u32> = (names.iter())
                    .filter_map(|g| row_of.get(g.as_str()).copied())
                    .collect();
                rows.sort_unstable();
                rows.dedup();
                rows
            }
            None => (0..n_genes as u32).collect(),
        };
        for bin in bins {
            for _ in 0..cfg.decoys {
                let pick = rand::seq::index::sample(&mut rng, pool.len(), bin.min(pool.len()));
                let rows: Vec<u32> = pick.iter().map(|i| pool[i]).collect();
                nodes.push((kind, rows, None));
            }
        }
    }
    anyhow::ensure!(
        nodes.iter().any(|n| n.2.is_some()),
        "no gene set has {} genes in the table",
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
    for (id, (kind, rows, _)) in nodes.iter().enumerate() {
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
        origin: (0..genes.dim)
            .map(|d| {
                (0..n_genes)
                    .map(|g| genes.rows[g * genes.dim + d])
                    .sum::<f32>()
                    / n_genes as f32
            })
            .collect(),
        names: Vec::new(),
        kinds: Vec::new(),
        centres: Vec::new(),
        dim: genes.dim,
        covered: Vec::new(),
        bins: Vec::new(),
        dropped,
        decoys: Vec::new(),
    };
    for (id, (kind, rows, set)) in nodes.iter().enumerate() {
        match set {
            Some(s) => {
                units.names.push(s.name.clone());
                units.kinds.push(*kind);
                units.centres.extend(row(id));
                units.covered.push((rows.len(), s.genes.len()));
                units.bins.push(size_bin(rows.len()).min(n_genes));
            }
            None => units.decoys.push(Decoy {
                kind: *kind,
                genes: rows.clone(),
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
    by_type.extend_from_slice(&terms.sets);
    let pools = [(SetKind::Term, terms.pool.clone())];
    let pools = if terms.pool.is_empty() {
        &[][..]
    } else {
        &pools[..]
    };
    let units = train_units_from(&table_of(feat), &by_type, pools, cfg)?;
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
    let term_alignment = if units.kinds.contains(&SetKind::Term) {
        let path = format!("{out}{TERM_ALIGNMENT_PARQUET}");
        write_alignment(&path, &align(&units, SetKind::Term, &flat, groups)?)?;
        Some(path)
    } else {
        None
    };
    Ok(AnnotationOutputs {
        argmax: Some(labels_path(out)),
        celltype_units: Some(format!("{out}{UNITS_PARQUET}")),
        celltype_links: Some(format!("{out}{TYPE_LINKS_TSV}")),
        term_alignment,
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
    /// The mean cosine of a cluster's cells (from the genes' origin) to the
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
    // Each cluster's mean unit vector, read from the genes' origin.
    let means: Vec<Vec<f32>> = (clusters.iter())
        .map(|&k| {
            let mut m = vec![0.0f32; dim];
            let members = (0..n).filter(|&c| groups[c] == k);
            let mut count = 0usize;
            for c in members {
                let from: Vec<f32> = (cells[c * dim..(c + 1) * dim].iter())
                    .zip(&units.origin)
                    .map(|(x, o)| x - o)
                    .collect();
                m.iter_mut().zip(unit(&from)).for_each(|(a, b)| *a += b);
                count += 1;
            }
            m.iter_mut().for_each(|a| *a /= count.max(1) as f32);
            m
        })
        .collect();
    let dot = |a: &[f32], b: &[f32]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f32>();
    let picked: Vec<usize> = (0..units.names.len())
        .filter(|&i| units.kinds[i] == kind)
        .collect();
    let decoys: Vec<&Decoy> = units.decoys.iter().filter(|d| d.kind == kind).collect();
    let decoy_dirs: Vec<Vec<f32>> = decoys.iter().map(|d| unit(&d.centre)).collect();
    let t = picked.len();
    let mut score = Vec::with_capacity(clusters.len() * t);
    let mut p = Vec::with_capacity(clusters.len() * t);
    let mut q = Vec::with_capacity(clusters.len() * t);
    for m in &means {
        let mut row_p = Vec::with_capacity(t);
        for &i in &picked {
            let s = dot(m, &unit(&units.centres[i * dim..(i + 1) * dim]));
            let null: Vec<f32> = (decoys.iter().zip(&decoy_dirs))
                .filter(|(d, _)| d.genes.len() == units.bins[i])
                .map(|(_, dir)| dot(m, dir))
                .collect();
            score.push(s);
            row_p.push(upper_tail(s, &null));
        }
        q.extend(legume_numeric::matrix::hypothesis::benjamini_hochberg(
            &row_p,
        ));
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

/// `P(Z ≥ z)` for `s` among `null`, `z` standardised by the null's mean and
/// sample standard deviation.
fn upper_tail(s: f32, null: &[f32]) -> f32 {
    use statrs::distribution::{ContinuousCDF, Normal};
    let n = null.len() as f64;
    if n < 2.0 {
        return f32::NAN;
    }
    let mean = null.iter().map(|&x| f64::from(x)).sum::<f64>() / n;
    let var = null
        .iter()
        .map(|&x| (f64::from(x) - mean).powi(2))
        .sum::<f64>()
        / (n - 1.0);
    let z = (f64::from(s) - mean) / var.sqrt().max(f64::MIN_POSITIVE);
    Normal::standard().sf(z) as f32
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
