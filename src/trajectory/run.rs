//! The `lupin trajectory` command: the prior over the run's cell types,
//! checked against the kNN graph, then diffusion pseudotime from the prior's
//! roots, written with a `trajectory` section into a copy of the manifest.

use super::diffusion::{DiffusionMap, Neighbours, DIFFMAP_COMPS};
use super::edges::{self, EdgeRow, Verdict};
use super::encode_groups;
use super::prior::{self, Prior, Source, Statement};
use super::type_connectivity::connectivity;
use crate::cell_labels::read_cell_labels;
use crate::manifest::data_files::Fetch;
use crate::manifest::run::{
    annotated_path, load, may_replace, rel_to_manifest, resolve, same_file, Loaded,
};
use crate::progress::Stages;
use anyhow::{bail, Context, Result};
use clap::Args;
use legume_numeric::matrix::common_io::mkdir_parent;
use legume_numeric::matrix::dense_mat_io::{axis_id_names, Mat};
use legume_numeric::matrix::parquet::{write_named_table, Column};
use legume_numeric::matrix::traits::IoOps;
use legume_numeric::matrix::utils::median;
use log::{info, warn};
use nalgebra::DMatrix;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Args, Debug)]
pub struct TrajectoryArgs {
    #[arg(
        short = 'f',
        long = "from",
        help = "Run manifest (or its prefix); the TUI asks for it when omitted"
    )]
    pub from: Option<Box<str>>,

    #[arg(
        short = 'o',
        long,
        help = "Output prefix; also receives a copy of the manifest with a `trajectory` section. \
                With --from it runs without the TUI; otherwise the TUI asks for it"
    )]
    pub out: Option<Box<str>>,

    #[arg(
        long,
        help = "A cell type the trajectory starts from (repeatable); it gets a `cli` edge to every node it cannot otherwise reach"
    )]
    pub root: Vec<Box<str>>,

    #[arg(
        long,
        help = "Per-cell labels (cell<TAB>type); defaults to the run's annotate.argmax"
    )]
    pub labels: Option<Box<str>>,

    #[arg(
        long,
        help = "A precedence file for this run (from<TAB>to<TAB>precedes|unrelated[<TAB>note]), or a {out}.trajectory_prior.tsv to replay"
    )]
    pub prior: Option<Box<str>>,

    #[arg(
        long,
        requires = "prior",
        help = "Use --prior alone: no Cell Ontology and no user or project precedence.tsv"
    )]
    pub prior_only: bool,

    #[arg(long, help = "Cell Ontology .obo (else the data files' search path)")]
    pub obo: Option<Box<str>>,

    #[arg(long, help = "label<TAB>CL:id aliases on top of the shipped ones")]
    pub label_cl: Option<Box<str>>,

    #[arg(
        long,
        default_value_t = 15,
        help = "Neighbours per cell, counting the cell itself"
    )]
    pub knn: usize,

    #[arg(
        long,
        default_value_t = 10,
        value_parser = clap::value_parser!(u16).range(1..),
        help = "Diffusion components in the pseudotime distance, of the 15 computed (as scanpy's dpt)"
    )]
    pub n_dcs: u16,

    #[arg(
        long,
        default_value_t = 20,
        help = "Types with fewer cells are not nodes of the prior; their cells follow the nearest root"
    )]
    pub min_cells: usize,

    #[arg(
        long,
        default_value_t = 0.1,
        help = "PAGA connectivity at or above which a prior edge is supported and a pair outside the prior is a candidate"
    )]
    pub min_connectivity: f32,

    #[arg(
        long,
        help = "Stop after the prior and its connectivity check; write no pseudotime"
    )]
    pub check_only: bool,

    #[arg(
        long,
        value_enum,
        default_value_t = crate::tui::Graphics::Auto,
        help = "In the TUI, how figures reach the terminal"
    )]
    pub graphics: crate::tui::Graphics,
}

impl TrajectoryArgs {
    /// Without both `--from` and `--out`, the command opens its TUI, which
    /// asks for what is missing.
    #[must_use]
    pub fn opens_tui(&self) -> bool {
        self.from.is_none() || self.out.is_none()
    }

    /// The options a run started from the TUI carries over, less `--from`,
    /// `--out` and `--graphics`.
    pub(crate) fn child_argv(&self) -> Vec<String> {
        let mut v: Vec<String> = Vec::new();
        let mut val = |flag: &str, x: String| {
            v.push(format!("--{flag}"));
            v.push(x);
        };
        for r in &self.root {
            val("root", r.to_string());
        }
        for (flag, x) in [
            ("labels", &self.labels),
            ("prior", &self.prior),
            ("obo", &self.obo),
            ("label-cl", &self.label_cl),
        ] {
            if let Some(x) = x {
                val(flag, x.to_string());
            }
        }
        val("knn", self.knn.to_string());
        val("n-dcs", self.n_dcs.to_string());
        val("min-cells", self.min_cells.to_string());
        val("min-connectivity", self.min_connectivity.to_string());
        if self.prior_only {
            v.push("--prior-only".into());
        }
        if self.check_only {
            v.push("--check-only".into());
        }
        v
    }
}

/// What the manifest at `dir` records about the run; input files
/// manifest-relative, as every path in a manifest is.
fn settings(args: &TrajectoryArgs, dir: &std::path::Path) -> serde_json::Value {
    let rel = |p: &Option<Box<str>>| p.as_deref().map(|p| rel_to_manifest(dir, p));
    serde_json::json!({
        "knn": args.knn, "n_dcs": args.n_dcs, "min_cells": args.min_cells,
        "min_connectivity": args.min_connectivity,
        "roots": args.root, "prior": rel(&args.prior), "prior_only": args.prior_only,
        "labels": rel(&args.labels), "obo": rel(&args.obo), "label_cl": rel(&args.label_cl),
    })
}

/// The cells, their prepared geometry and their types.
struct Inputs {
    pub(crate) cells: Vec<Box<str>>,
    pub(crate) geometry: Mat,
    /// The distinct labels, sorted.
    pub(crate) names: Vec<Box<str>>,
    /// Each cell's index into `names`.
    pub(crate) group: Vec<usize>,
    /// Each type's cells.
    pub(crate) cells_of: Vec<Vec<usize>>,
    /// Types with at least `min_cells` cells, other than unassigned.
    pub(crate) is_node: Vec<bool>,
}

impl Inputs {
    fn label(&self, cell: usize) -> &str {
        self.names[self.group[cell]].as_ref()
    }
}

/// A direct prior edge with its data verdict.
#[derive(Debug, Clone)]
struct EdgeCheck {
    pub(crate) from: usize,
    pub(crate) to: usize,
    pub(crate) verdict: Verdict,
    /// Fraction of `to`'s cells beyond `from`'s median pseudotime; NaN before
    /// pseudotime is computed.
    pub(crate) order_agreement: f32,
}

/// What a run computes.
struct Trajectory {
    pub(crate) prior: Prior,
    /// PAGA connectivity between every pair of types.
    pub(crate) connectivity: DMatrix<f64>,
    pub(crate) edges: Vec<EdgeCheck>,
    /// Pairs of nodes the prior does not order whose connectivity reaches the
    /// threshold, strongest first.
    pub(crate) candidates: Vec<(usize, usize)>,
    /// `None` before pseudotime is computed (`--check-only`).
    pub(crate) ordering: Option<Ordering>,
}

/// Pseudotime and lineages.
struct Ordering {
    pub(crate) map: DiffusionMap,
    /// Per cell; NaN for cells no root reaches.
    pub(crate) pseudotime: Vec<f32>,
    /// Per cell, the prior component (`None` when no root reaches it).
    pub(crate) component: Vec<Option<usize>>,
    /// Root-to-leaf paths over node indices, per component.
    pub(crate) lineages: Vec<(usize, Vec<usize>)>,
    /// Types × lineages: a type's weight is uniform over the paths through it.
    pub(crate) type_weights: Vec<Vec<f32>>,
    /// Each type's finite pseudotimes.
    pub(crate) by_type: Vec<Vec<f32>>,
    /// Node types in order of median pseudotime: `(median, type)`; a type no
    /// root reaches has an infinite median.
    pub(crate) order: Vec<(f32, usize)>,
}

/// With both `--from` and `--out`, the run; otherwise the TUI, which asks for
/// what is missing.
pub fn run_trajectory(args: &TrajectoryArgs) -> Result<()> {
    match (args.from.as_deref(), args.out.as_deref()) {
        (Some(from), Some(out)) if !args.opens_tui() => run_batch(args, from, out),
        _ => crate::tui::run_trajectory(args),
    }
}

/// The batch run's stages for the TUI's progress popup, weighted by their
/// rough share of the time: the kNN graph and the eigensolve dominate.
pub(crate) const STAGES: Stages = Stages::new(&[
    ("reading the run", 1.0),
    ("building the prior", 0.5),
    ("kNN graph", 4.0),
    ("connectivity check", 0.5),
    ("diffusion components", 4.0),
    ("pseudotime and lineages", 1.0),
    ("writing outputs", 1.0),
]);
/// [`STAGES`] by name.
pub(crate) const STAGE_DIFFUSION: usize = 4;

fn run_batch(args: &TrajectoryArgs, from: &str, out: &str) -> Result<()> {
    STAGES.start(0);
    let loaded = load(from)?;
    let manifest_out = annotated_path(&loaded.file, out);
    // Refused before anything is written: the outputs would replace the
    // input's own tables and the manifest copy could not be made.
    anyhow::ensure!(
        !same_file(&loaded.file, &manifest_out),
        "{} is the manifest this command reads from; choose a different --out",
        manifest_out.display()
    );
    may_replace(&manifest_out)?;
    mkdir_parent(out)?;
    let inputs = load_inputs(&loaded, args)?;
    STAGES.start(1);
    let mut layers = gather_statements(&loaded, args, &inputs)?;
    let roots: Vec<&str> = args.root.iter().map(AsRef::as_ref).collect();
    if !roots.is_empty() {
        let so_far = prior::combine(&layers);
        layers.push(prior::root_layer(
            &inputs.names,
            &inputs.is_node,
            &so_far,
            &roots,
        )?);
    }
    let prior = prior::build(&inputs.names, &inputs.is_node, prior::combine(&layers))?;
    report_prior(&inputs, &prior);
    if prior.edges.is_empty() {
        bail!(
            "the prior orders no pair of types, so there is nothing to check or to order from: \
             name a root (--root TYPE), add statements to precedence.tsv or --prior, or map the \
             types to Cell Ontology terms (--label-cl)"
        );
    }

    STAGES.start(2);
    let nb = Neighbours::new(&inputs.geometry, args.knn)?;
    STAGES.start(3);
    let mut t = check(&inputs, &nb, prior, args.min_connectivity);
    report_check(&inputs, &t);
    if !args.check_only {
        t.ordering = Some(order(&inputs, &nb, &t.prior, usize::from(args.n_dcs))?);
        agreement(&mut t);
        report_order(&inputs, &t);
    }
    STAGES.start(6);
    let written = write(&inputs, &t, out)?;
    record(&loaded, &manifest_out, &written, args)?;
    STAGES.finish();
    info!(
        "wrote {}",
        written.values().cloned().collect::<Vec<_>>().join(", ")
    );
    Ok(())
}

fn load_inputs(loaded: &Loaded, args: &TrajectoryArgs) -> Result<Inputs> {
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
    let mut missing = Vec::new();
    let labels: Vec<&str> = x
        .rows
        .iter()
        .map(|c| {
            membership.get(c).unwrap_or_else(|| {
                missing.push(c.as_ref());
                enrichment::UNASSIGNED_LABEL
            })
        })
        .collect();
    if missing.len() == x.rows.len() {
        bail!(
            "none of the {} cells is in {label_path} (e.g. {}); do the cell ids match?",
            x.rows.len(),
            missing
                .iter()
                .take(3)
                .copied()
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
    if !missing.is_empty() {
        warn!(
            "{}/{} cells are not in {label_path}; they count as {}",
            missing.len(),
            x.rows.len(),
            enrichment::UNASSIGNED_LABEL
        );
    }
    let (names, group) = encode_groups(&labels);
    let mut cells_of = vec![Vec::new(); names.len()];
    for (i, &g) in group.iter().enumerate() {
        cells_of[g].push(i);
    }
    let is_node: Vec<bool> = (0..names.len())
        .map(|g| {
            cells_of[g].len() >= args.min_cells && names[g].as_ref() != enrichment::UNASSIGNED_LABEL
        })
        .collect();
    let small: Vec<String> = (0..names.len())
        .filter(|&g| !is_node[g] && names[g].as_ref() != enrichment::UNASSIGNED_LABEL)
        .map(|g| format!("{} ({})", names[g], cells_of[g].len()))
        .collect();
    if !small.is_empty() {
        info!(
            "{} type(s) below --min-cells {} are not nodes: {}",
            small.len(),
            args.min_cells,
            small.join(", ")
        );
    }
    Ok(Inputs {
        cells: x.rows,
        geometry: x.mat,
        names,
        group,
        cells_of,
        is_node,
    })
}

/// The prior's statements, layer by layer: the Cell Ontology, the user's and
/// the project's `precedence.tsv`, then `--prior`.
fn gather_statements(
    loaded: &Loaded,
    args: &TrajectoryArgs,
    inputs: &Inputs,
) -> Result<Vec<Vec<Statement>>> {
    let mut layers: Vec<Vec<Statement>> = Vec::new();
    if !args.prior_only {
        // The run's marker panel, for its alias sidecar.
        let panel = crate::annotate_cmd::recorded_markers(loaded).unwrap_or_default();
        let data = crate::manifest::ontology::load(
            Some(&loaded.dir),
            &panel,
            args.obo.as_deref(),
            args.label_cl.as_deref(),
            Fetch::Allowed,
        )?;
        let node_names: Vec<&str> = (0..inputs.names.len())
            .filter(|&g| inputs.is_node[g])
            .map(|g| inputs.names[g].as_ref())
            .collect();
        let terms = data.terms()?;
        if terms.is_none() {
            warn!("no Cell Ontology at hand; the prior comes from precedence files alone");
        }
        let (standing, problems) = prior::standing_layers(terms, &node_names, &data.search);
        if let Some(e) = problems.into_iter().next() {
            return Err(e);
        }
        layers = standing;
    }
    if let Some(p) = args.prior.as_deref() {
        let text = std::fs::read_to_string(p).with_context(|| format!("reading --prior {p}"))?;
        let st = prior::parse_statements(&text, Source::Run, p)?;
        info!("--prior {p}: {} statement(s)", st.len());
        layers.push(st);
    }
    Ok(layers)
}

fn report_prior(inputs: &Inputs, prior: &Prior) {
    info!(
        "prior: {} statement(s), {} direct edge(s), {} component(s)",
        prior.statements.len(),
        prior.edges.len(),
        prior.roots.len()
    );
    for (c, roots) in prior.roots.iter().enumerate() {
        let names: Vec<&str> = roots.iter().map(|&r| inputs.names[r].as_ref()).collect();
        info!("  component {c}: root(s) {}", names.join(", "));
    }
    let loose: Vec<&str> = (0..inputs.names.len())
        .filter(|&g| inputs.is_node[g] && prior.component[g].is_none())
        .map(|g| inputs.names[g].as_ref())
        .collect();
    if !loose.is_empty() {
        warn!(
            "{} node type(s) have no prior edge and get no pseudotime: {}",
            loose.len(),
            loose.join(", ")
        );
    }
}

/// PAGA connectivity over all types, verdicts on the direct edges, and
/// candidate pairs the prior does not order.
fn check(inputs: &Inputs, nb: &Neighbours, prior: Prior, min_connectivity: f32) -> Trajectory {
    let relabelled = (0..inputs.cells.len())
        .filter(|&c| inputs.group[c] != inputs.group[nb.reps[nb.rep_of[c]]])
        .count();
    if relabelled > 0 {
        warn!(
            "{relabelled} cells coincide with a cell of another type; connectivity counts \
             them under that cell's type"
        );
    }
    let conn = connectivity(nb, &inputs.group, inputs.names.len());
    let n = inputs.names.len();
    let edges = prior
        .edges
        .iter()
        .map(|e| EdgeCheck {
            from: e.from,
            to: e.to,
            verdict: if conn[(e.from, e.to)] as f32 >= min_connectivity {
                Verdict::Supported
            } else {
                Verdict::Unsupported
            },
            order_agreement: f32::NAN,
        })
        .collect();
    let mut candidates: Vec<(usize, usize)> = (0..n)
        .flat_map(|a| (a + 1..n).map(move |b| (a, b)))
        .filter(|&(a, b)| {
            inputs.is_node[a]
                && inputs.is_node[b]
                && !prior.related(a, b)
                && conn[(a, b)] as f32 >= min_connectivity
        })
        .collect();
    candidates.sort_by(|&p, &q| conn[q].total_cmp(&conn[p]));
    Trajectory {
        prior,
        connectivity: conn,
        edges,
        candidates,
        ordering: None,
    }
}

fn report_check(inputs: &Inputs, t: &Trajectory) {
    let n_sup = t
        .edges
        .iter()
        .filter(|e| e.verdict == Verdict::Supported)
        .count();
    info!(
        "connectivity check: {n_sup}/{} prior edge(s) supported, {} candidate pair(s) the prior does not order",
        t.edges.len(),
        t.candidates.len()
    );
    for e in t.edges.iter().filter(|e| e.verdict != Verdict::Supported) {
        info!(
            "  unsupported: {} → {} (connectivity {:.3})",
            inputs.names[e.from],
            inputs.names[e.to],
            t.connectivity[(e.from, e.to)]
        );
    }
    for &(a, b) in t.candidates.iter().take(10) {
        info!(
            "  candidate: {} — {} (connectivity {:.3})",
            inputs.names[a],
            inputs.names[b],
            t.connectivity[(a, b)]
        );
    }
}

/// Diffusion pseudotime from each component's roots (the root types'
/// medoids), scaled to [0, 1] per component, and the lineages.
fn order(inputs: &Inputs, nb: &Neighbours, prior: &Prior, n_dcs: usize) -> Result<Ordering> {
    STAGES.start(STAGE_DIFFUSION);
    let map = DiffusionMap::new(nb, DIFFMAP_COMPS, n_dcs)?;
    STAGES.start(STAGE_DIFFUSION + 1);
    let n = inputs.cells.len();
    // Distance to the nearest root of each component.
    let mut dist: Vec<Vec<f64>> = Vec::with_capacity(prior.roots.len());
    for (c, roots) in prior.roots.iter().enumerate() {
        let mut cells = Vec::with_capacity(roots.len());
        for &g in roots {
            let r = map.medoid(&inputs.cells_of[g]).context("no root cell")?;
            info!(
                "component {c}: root {} is the medoid of {}",
                inputs.cells[r], inputs.names[g]
            );
            cells.push(r);
        }
        dist.push(map.distances_from(&cells));
    }
    // A node type's cell belongs to its type's component; any other cell to
    // the component whose root is nearest.
    let component: Vec<Option<usize>> = (0..n)
        .map(|i| {
            let g = inputs.group[i];
            if inputs.is_node[g] {
                prior.component[g]
            } else {
                (0..dist.len())
                    .filter(|&c| dist[c][i].is_finite())
                    .min_by(|&a, &b| dist[a][i].total_cmp(&dist[b][i]))
            }
        })
        .collect();
    let mut top = vec![0.0f64; dist.len()];
    for i in 0..n {
        if let Some(c) = component[i] {
            if dist[c][i].is_finite() {
                top[c] = top[c].max(dist[c][i]);
            }
        }
    }
    let pseudotime: Vec<f32> = (0..n)
        .map(|i| match component[i] {
            Some(c) if dist[c][i].is_finite() => {
                (dist[c][i] / top[c].max(f64::MIN_POSITIVE)) as f32
            }
            _ => f32::NAN,
        })
        .collect();

    let lineages = prior.lineages();
    let type_weights: Vec<Vec<f32>> = (0..inputs.names.len())
        .map(|g| {
            let through: Vec<bool> = lineages.iter().map(|(_, p)| p.contains(&g)).collect();
            let k = through.iter().filter(|&&t| t).count() as f32;
            through
                .iter()
                .map(|&t| if t { 1.0 / k } else { 0.0 })
                .collect()
        })
        .collect();

    let mut by_type: Vec<Vec<f32>> = vec![Vec::new(); inputs.names.len()];
    for (&g, &pt) in inputs.group.iter().zip(&pseudotime) {
        if pt.is_finite() {
            by_type[g].push(pt);
        }
    }
    let mut order: Vec<(f32, usize)> = (0..inputs.names.len())
        .filter(|&g| inputs.is_node[g])
        .map(|g| {
            let m = if by_type[g].is_empty() {
                f32::INFINITY
            } else {
                median(&by_type[g])
            };
            (m, g)
        })
        .collect();
    order.sort_by(|p, q| p.0.total_cmp(&q.0));
    Ok(Ordering {
        map,
        pseudotime,
        component,
        lineages,
        type_weights,
        by_type,
        order,
    })
}

/// For each direct edge A → B, the fraction of B's cells beyond A's median.
fn agreement(t: &mut Trajectory) {
    let Some(o) = &t.ordering else { return };
    for e in &mut t.edges {
        let (from, to) = (&o.by_type[e.from], &o.by_type[e.to]);
        if from.is_empty() || to.is_empty() {
            continue;
        }
        let m = median(from);
        e.order_agreement = to.iter().filter(|&&p| p > m).count() as f32 / to.len() as f32;
    }
}

fn report_order(inputs: &Inputs, t: &Trajectory) {
    let Some(o) = &t.ordering else { return };
    let unreached = o.pseudotime.iter().filter(|p| !p.is_finite()).count();
    if unreached > 0 {
        warn!(
            "{unreached}/{} cells are reached by no root (another kNN component, or a type \
             without a prior edge) and get no pseudotime",
            o.pseudotime.len()
        );
    }
    info!("eigenvalues {:?}", o.map.evals);
    info!("median pseudotime by node type:");
    for &(m, g) in &o.order {
        info!(
            "  {m:.3}  {} ({})",
            inputs.names[g],
            inputs.cells_of[g].len()
        );
    }
    info!("order agreement of the prior's edges:");
    for e in &t.edges {
        info!(
            "  {} → {}: {:.2} ({})",
            inputs.names[e.from],
            inputs.names[e.to],
            e.order_agreement,
            e.verdict.as_str()
        );
    }
    let paths: Vec<String> = o
        .lineages
        .iter()
        .map(|(_, p)| {
            p.iter()
                .map(|&g| inputs.names[g].as_ref())
                .collect::<Vec<_>>()
                .join(" > ")
        })
        .collect();
    info!("{} lineage(s): {}", o.lineages.len(), paths.join("; "));
}

/// The outputs, keyed by their manifest slot.
fn write(inputs: &Inputs, t: &Trajectory, out: &str) -> Result<BTreeMap<&'static str, String>> {
    let mut written = BTreeMap::new();
    let prior_path = format!("{out}.trajectory_prior.tsv");
    prior::write_tsv(&prior_path, &inputs.names, &t.prior)?;
    written.insert("prior", prior_path);
    written.insert("edges", write_edges(inputs, t, out)?);
    match &t.ordering {
        Some(o) => written.extend(write_ordering(inputs, o, out)?),
        // A check onto a prefix that held a full run: its ordering no longer
        // goes with the prior written above.
        None => {
            for stale in [
                format!("{out}.cell_pseudotime.parquet"),
                format!("{out}.diffusion.parquet"),
                format!("{out}.trajectory_lineages.tsv"),
            ] {
                if std::path::Path::new(&stale).is_file() {
                    std::fs::remove_file(&stale).with_context(|| format!("removing {stale}"))?;
                }
            }
        }
    }
    Ok(written)
}

/// `{out}.trajectory_edges.parquet`: every node pair, so the table is the
/// connectivity matrix as well as the verdicts.
fn write_edges(inputs: &Inputs, t: &Trajectory, out: &str) -> Result<String> {
    let edge_at: BTreeMap<(usize, usize), &EdgeCheck> = t
        .edges
        .iter()
        .flat_map(|e| [((e.from, e.to), e), ((e.to, e.from), e)])
        .collect();
    let candidate: BTreeSet<(usize, usize)> = t.candidates.iter().copied().collect();
    let nodes: Vec<usize> = (0..inputs.names.len())
        .filter(|&g| inputs.is_node[g])
        .collect();
    let mut rows = Vec::new();
    for (i, &p) in nodes.iter().enumerate() {
        for &q in &nodes[i + 1..] {
            let edge = edge_at.get(&(p, q)).copied();
            let (from, to) = edge.map_or((p, q), |e| (e.from, e.to));
            rows.push(EdgeRow {
                a: inputs.names[from].clone(),
                b: inputs.names[to].clone(),
                connectivity: t.connectivity[(p, q)] as f32,
                in_prior: edge.is_some(),
                verdict: match edge {
                    Some(e) => Some(e.verdict),
                    None if candidate.contains(&(p, q)) => Some(Verdict::Candidate),
                    None => None,
                },
                order_agreement: edge.map_or(f32::NAN, |e| e.order_agreement),
            });
        }
    }
    let path = format!("{out}.trajectory_edges.parquet");
    edges::write(&path, &rows)?;
    Ok(path)
}

/// `{out}.cell_pseudotime.parquet`, `{out}.diffusion.parquet` and
/// `{out}.trajectory_lineages.tsv`.
fn write_ordering(
    inputs: &Inputs,
    o: &Ordering,
    out: &str,
) -> Result<BTreeMap<&'static str, String>> {
    let mut written = BTreeMap::new();
    let n = inputs.cells.len();
    let labels: Vec<Box<str>> = (0..n).map(|i| inputs.label(i).into()).collect();
    let component: Vec<i32> = o
        .component
        .iter()
        .map(|c| c.map_or(-1, |c| c as i32))
        .collect();
    let lineage_cols: Vec<Vec<f32>> = (0..o.lineages.len())
        .map(|k| (0..n).map(|i| o.type_weights[inputs.group[i]][k]).collect())
        .collect();
    let mut cols: Vec<(Box<str>, Column)> = vec![
        ("pseudotime".into(), Column::F32(&o.pseudotime)),
        ("type".into(), Column::Str(&labels)),
        ("component".into(), Column::I32(&component)),
    ];
    for (k, c) in lineage_cols.iter().enumerate() {
        cols.push((format!("L{k}").into(), Column::F32(c)));
    }
    let pt_path = format!("{out}.cell_pseudotime.parquet");
    write_named_table(&pt_path, "cell", &inputs.cells, &cols)?;
    written.insert("pseudotime", pt_path);

    let dc = o.map.evecs.map(|v| v as f32);
    let dc_path = format!("{out}.diffusion.parquet");
    dc.to_parquet_with_names(
        &dc_path,
        (Some(&inputs.cells), Some("cell")),
        Some(&axis_id_names("DC", dc.ncols())),
    )?;
    written.insert("diffusion", dc_path);

    let mut text = String::from("# lineage\tcomponent\tpath\n");
    for (k, (c, p)) in o.lineages.iter().enumerate() {
        let path: Vec<&str> = p.iter().map(|&g| inputs.names[g].as_ref()).collect();
        text += &format!("L{k}\t{c}\t{}\n", path.join(" > "));
    }
    let lin_path = format!("{out}.trajectory_lineages.tsv");
    std::fs::write(&lin_path, text).with_context(|| format!("writing {lin_path}"))?;
    written.insert("lineages", lin_path);
    Ok(written)
}

/// A copy of the run's manifest at `manifest_out` with the `trajectory`
/// section filled in.
fn record(
    loaded: &Loaded,
    manifest_out: &std::path::Path,
    written: &BTreeMap<&'static str, String>,
    args: &TrajectoryArgs,
) -> Result<()> {
    let mut copy = loaded.copy_to(manifest_out.to_path_buf())?;
    let settings = settings(args, &copy.dir);
    let rel = |k: &str| written.get(k).map(|p| rel_to_manifest(&copy.dir, p));
    let t = &mut copy.manifest.trajectory;
    t.prior = rel("prior");
    t.edges = rel("edges");
    t.pseudotime = rel("pseudotime");
    t.diffusion = rel("diffusion");
    t.lineages = rel("lineages");
    t.settings = Some(settings);
    copy.manifest.save(&copy.file)
}
