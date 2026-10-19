//! A units round's gene-set units on its cell map: each cell type at the
//! median of its anchors (the cells nearest its centre, `unit_anchors`) on
//! the layout drawn, the strongest links between types as edges, a ring of
//! its fitted radius, and the selected cluster's significant GO terms as
//! numbered stars at their own anchors. Any layout the run has reads them,
//! as the anchors are cells.

use super::{median, CellIndex};
use crate::manifest::run::{resolve, Loaded};
use anyhow::{Context, Result};
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::parquet::read_table_columns;
use legume_numeric::matrix::traits::IoOps;
use std::collections::BTreeMap;
use std::path::Path;

/// A cluster's terms starred: at most this many, best score first.
pub(crate) const STARS: usize = 6;
/// The q below which a term is starred.
pub(crate) const STAR_Q: f32 = 0.05;

/// A unit's anchors: the round's cells nearest it and their cosines.
#[derive(Debug, Clone, Default)]
pub(crate) struct Unit {
    pub(crate) name: Box<str>,
    pub(crate) cells: Vec<usize>,
    pub(crate) weight: Vec<f32>,
    /// A cell type's radius on the unit sphere; NaN for a term.
    pub(crate) radius: f32,
}

/// A term aligned with a cluster: its id, its name (the id without the
/// ontology), score and q.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct GoTerm {
    pub(crate) id: Box<str>,
    pub(crate) name: String,
    pub(crate) score: f32,
    pub(crate) q: f32,
}

/// What a units round says of its units, read once.
#[derive(Debug, Default)]
pub(crate) struct UnitMap {
    pub(crate) types: Vec<Unit>,
    pub(crate) terms: BTreeMap<Box<str>, Unit>,
    /// Each linked pair of types (indices into `types`) and its link.
    pub(crate) links: Vec<(usize, usize, f32)>,
    /// Each cluster's starred terms ([`top_terms`]).
    pub(crate) by_cluster: BTreeMap<i32, Vec<GoTerm>>,
}

/// A cell type on one layout: where it sits, and its ring's radius in the
/// layout's units, when its anchors spread.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TypeMark {
    pub(crate) name: Box<str>,
    pub(crate) at: Option<(f32, f32)>,
    pub(crate) ring: Option<f32>,
}

/// The units on one layout: the types, the links kept (their strength as
/// a share of the strongest), and the cluster's stars, numbered in order.
#[derive(Debug, Clone, Default)]
pub(crate) struct Marks {
    pub(crate) types: Vec<TypeMark>,
    pub(crate) links: Vec<(usize, usize, f32)>,
    pub(crate) stars: Vec<(GoTerm, (f32, f32))>,
}

impl UnitMap {
    /// The units `loaded` records onto its cells (`index`); `None` when it
    /// records no anchors.
    pub(crate) fn read(loaded: &Loaded, index: &CellIndex) -> Result<Option<Self>> {
        let a = &loaded.manifest.annotate;
        let at = |rel: &Option<String>| rel.as_deref().map(|r| resolve(&loaded.dir, r));
        let Some(path) = at(&a.unit_anchors) else {
            return Ok(None);
        };
        let (s, w) = read_table_columns(&path, &["unit", "kind", "cell"], &["weight"])
            .with_context(|| format!("reading {path}"))?;
        let mut types: Vec<Unit> = Vec::new();
        let mut terms: BTreeMap<Box<str>, Unit> = BTreeMap::new();
        for i in 0..s[0].len() {
            let Some(&cell) = index.first.get(&*s[2][i]) else {
                continue;
            };
            let name = &s[0][i];
            let unit = if &*s[1][i] == "term" {
                terms.entry(name.clone()).or_insert_with(|| Unit {
                    name: name.clone(),
                    radius: f32::NAN,
                    ..Unit::default()
                })
            } else {
                match types.iter().position(|t| t.name == *name) {
                    Some(k) => &mut types[k],
                    None => {
                        types.push(Unit {
                            name: name.clone(),
                            radius: f32::NAN,
                            ..Unit::default()
                        });
                        types.last_mut().expect("just pushed")
                    }
                }
            };
            unit.cells.push(cell);
            unit.weight.push(w[0][i] as f32);
        }
        if let Some(p) = at(&a.celltype_units) {
            let t = Mat::from_parquet_with_row_names(&p, Some(0))
                .with_context(|| format!("reading {p}"))?;
            if let Some(j) = t.cols.iter().position(|c| &**c == "radius") {
                for (r, name) in t.rows.iter().enumerate() {
                    if let Some(u) = types.iter_mut().find(|u| u.name == *name) {
                        u.radius = t.mat[(r, j)];
                    }
                }
            }
        }
        let links = match at(&a.celltype_links) {
            Some(p) => read_links(&p, &types)?,
            None => Vec::new(),
        };
        let by_cluster = match at(&a.term_alignment) {
            Some(p) => {
                let go = go_ontology(loaded);
                read_alignment(&p, |id| Some(go.as_ref()?.name(id)?.to_string()))?
            }
            None => BTreeMap::new(),
        };
        Ok(Some(Self {
            types,
            terms,
            links,
            by_cluster,
        }))
    }

    /// Cluster `k`'s starred terms; none without a cluster.
    pub(crate) fn terms_of(&self, k: Option<i32>) -> &[GoTerm] {
        k.and_then(|k| self.by_cluster.get(&k))
            .map_or(&[], Vec::as_slice)
    }

    /// The units on the layout `(x, y)`, with the stars of cluster `focus`.
    pub(crate) fn marks(&self, x: &[f32], y: &[f32], focus: Option<i32>) -> Marks {
        let types = (self.types.iter())
            .map(|u| {
                let at = anchor_median(u, x, y);
                TypeMark {
                    name: u.name.clone(),
                    at,
                    ring: at.and_then(|m| ring(u, m, x, y)),
                }
            })
            .collect();
        // The strongest links, as many as there are types.
        let mut links: Vec<(usize, usize, f32)> = (self.links.iter())
            .copied()
            .filter(|l| l.2.is_finite() && l.2 > 0.0)
            .collect();
        links.sort_by(|a, b| b.2.total_cmp(&a.2));
        links.truncate(self.types.len().max(1));
        let top = links.first().map_or(1.0, |l| l.2);
        let links = links.into_iter().map(|(a, b, v)| (a, b, v / top)).collect();
        let stars = (self.terms_of(focus).iter())
            .filter_map(|t| {
                let at = anchor_median(self.terms.get(&t.id)?, x, y)?;
                Some((t.clone(), at))
            })
            .collect();
        Marks {
            types,
            links,
            stars,
        }
    }
}

/// The median of `u`'s anchors with coordinates on `(x, y)`.
fn anchor_median(u: &Unit, x: &[f32], y: &[f32]) -> Option<(f32, f32)> {
    let (xs, ys): (Vec<f32>, Vec<f32>) = (u.cells.iter())
        .filter(|&&c| x[c].is_finite() && y[c].is_finite())
        .map(|&c| (x[c], y[c]))
        .unzip();
    (!xs.is_empty()).then(|| (median(&xs), median(&ys)))
}

/// `u`'s radius on the layout: its radius on the sphere times its anchors'
/// spread on the layout (root mean square from `at`) over their spread on
/// the sphere (the root mean square chord `√(2 − 2w)` from its centre).
fn ring(u: &Unit, at: (f32, f32), x: &[f32], y: &[f32]) -> Option<f32> {
    let rms = |v: &mut dyn Iterator<Item = f32>| {
        let (s, n) = v.fold((0.0, 0usize), |(s, n), d| (s + d * d, n + 1));
        (n > 0).then(|| (s / n as f32).sqrt())
    };
    let on = (u.cells.iter()).filter(|&&c| x[c].is_finite() && y[c].is_finite());
    let layout = rms(&mut on.map(|&c| (x[c] - at.0).hypot(y[c] - at.1)))?;
    let sphere = rms(&mut u.weight.iter().map(|w| (2.0 - 2.0 * w).max(0.0).sqrt()))?;
    let r = u.radius * layout / sphere;
    (r.is_finite() && r > 0.0).then_some(r)
}

/// The terms of `rows` with `q` below [`STAR_Q`], best score first, at
/// most [`STARS`].
pub(crate) fn top_terms(mut rows: Vec<GoTerm>) -> Vec<GoTerm> {
    rows.retain(|t| t.q < STAR_Q && t.score.is_finite());
    rows.sort_by(|a, b| b.score.total_cmp(&a.score));
    rows.truncate(STARS);
    rows
}

/// The pairs of `types` linked in the table at `path` (`type_a`, `type_b`,
/// `link`).
fn read_links(path: &str, types: &[Unit]) -> Result<Vec<(usize, usize, f32)>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
    let at = |n: &str| types.iter().position(|t| &*t.name == n);
    Ok(text
        .lines()
        .skip(1)
        .filter_map(|l| {
            let mut f = l.split('\t');
            let (a, b) = (at(f.next()?)?, at(f.next()?)?);
            Some((a, b, f.next()?.trim().parse().ok()?))
        })
        .collect())
}

/// The alignment table at `path`, as each cluster's starred terms, each
/// by its `name` where it has one.
fn read_alignment(
    path: &str,
    name: impl Fn(&str) -> Option<String>,
) -> Result<BTreeMap<i32, Vec<GoTerm>>> {
    let (s, v) = read_table_columns(path, &["cluster", "term"], &["score", "q"])
        .with_context(|| format!("reading {path}"))?;
    let mut by: BTreeMap<i32, Vec<GoTerm>> = BTreeMap::new();
    for i in 0..s[0].len() {
        let Some(k) = crate::annotate::rounds::parse_cluster_id(&s[0][i]) else {
            continue;
        };
        let id = s[1][i].clone();
        by.entry(k as i32).or_default().push(GoTerm {
            name: name(&id).unwrap_or_else(|| id.to_string()),
            id,
            score: v[0][i] as f32,
            q: v[1][i] as f32,
        });
    }
    Ok(by
        .into_iter()
        .map(|(k, rows)| (k, top_terms(rows)))
        .collect())
}

/// The ontology the pass read its terms from (`go_obo` in its settings,
/// as given or beside the manifest); none when it is gone.
fn go_ontology(loaded: &Loaded) -> Option<data_beans::aux::ontology::Ontology> {
    let settings = loaded.manifest.annotate.settings.as_ref();
    let obo = settings.and_then(|s| s.get("go_obo")?.as_str())?;
    let path = if Path::new(obo).is_file() {
        obo.to_string()
    } else {
        resolve(&loaded.dir, obo)
    };
    data_beans::aux::ontology::Ontology::load_obo(&path)
        .map_err(|e| log::warn!("GO term names: {e:#}; the terms go by their ids"))
        .ok()
}

#[cfg(test)]
#[path = "units_tests.rs"]
pub(crate) mod tests;
