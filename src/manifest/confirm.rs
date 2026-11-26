//! A round's calls checked on data that did not make its clusters: Gaussian
//! data fission of the cell embedding (Leiner, Duan, Wasserman & Ramdas,
//! "Data fission: splitting a single data point", *JASA*, 2025;
//! arXiv:2112.11079), and a held-out half of each type's markers. A line of
//! evidence beside the round's others, which it leaves as they are;
//! `lupin why` weighs it as its own line. `lupin docs fission` has the
//! method in full.
//!
//! # Fission
//!
//! With `w ~ N(0, Σ)` drawn per cell and `τ > 0`, the embedding `z` is split
//! into `f = z + τw` and `g = z − w/τ`. When `z_i ~ N(μ_i, Σ)`, `f` and `g`
//! are independent, so the clusters chosen on `f` (by legume-numeric's
//! Leiden, as the run's own) are fixed in advance as far as `g` is
//! concerned. The cells are not Gaussian: the split is independent only as
//! far as the embedding is asymptotically normal, a central-limit argument
//! and no more. `τ = 1` splits evenly; a smaller `τ` keeps more for the
//! clusters and less for the test.
//!
//! `Σ` is the cells' noise, which `z` alone does not tell apart from their
//! spread of means. A wrong `Σ` leaves `Cov(f, g) = Σ_noise − Σ`: too small
//! (the round's over-split clusters) and a split confirms itself; too large
//! (the cells' total spread, where types differ) and `g` is shifted away
//! from wherever the noise split a type. It is fitted instead to the
//! condition that holds when it is right: `f` and `g` do not co-vary within
//! any cluster of `f` ([`fitted_split`]). What dependence is left, in units
//! of what independent halves show by chance, is reported.
//!
//! # The ranking and its two tests
//!
//! A call needs only each cluster's ranking of the genes. Cluster `k`'s is
//! read from `g` through the decoder's gene table `ρ`:
//!
//! ```text
//! t_gk = ρ_g · (ḡ_k − m) / √(ρ_gᵀ Σ_g ρ_g / n_k)
//! ```
//!
//! with `ḡ_k` the cluster's mean of `g`, `m` the mean of the cluster means
//! and `Σ_g` the covariance of `g` within the clusters: the cluster's mean
//! logit for the gene against the balanced centre, in standard errors. Its
//! sufficient statistics are the clusters' sizes and means of `g` and one
//! pooled covariance. Each type is tested twice on it, and stands only on
//! both:
//!
//! - against random gene sets of its size, by the enrichment pass's own
//!   engine (weighted KS, multilevel tail, BH within the cluster), with no
//!   specificity transform, since the ranking is already a contrast. The
//!   ranking is a projection of each gene's `ρ_g` on one direction, so
//!   markers with like rows move together, which random sets do not; alone,
//!   this test is anti-conservative;
//! - by its markers' mean, `u_T · (ḡ_k − m)`, against its exact variance
//!   under Gaussian `g` ([`Ranking::up_p`]): the markers are up in the
//!   cluster, their correlation counted.
//!
//! The counts are not split: they made the embedding, and the gene table
//! with it, so a test that reads them (the round's enrichment and SuSiE)
//! stays post-selection. Fission removes the double dip of clustering and
//! testing on one embedding, given that embedding.
//!
//! # Held-out markers
//!
//! Each type's markers are split in two halves by a hash of the gene name.
//! One half picks the cluster's call on the ranking from `f`; the other
//! half tests it on the ranking from `g`, both ways. A call chosen on `g`
//! would carry over to the other half through the markers' correlation.
//!
//! # What still reads the round
//!
//! The round's within-cluster covariance starts the fit of `Σ` (which the
//! fit's condition then corrects). The clusters of `f` are matched to the
//! round's by the cells they share, and a call is said to agree with the
//! round's label: both describe, and change no statistic.

use super::annotate::{leiden_on, score_tables, LeidenArgs};
use super::rounds::{labels_and_tree, read_clusters};
use super::run::{prepare_cells, CellSpace};
use crate::annotate::markers::{annotation_matrix_from_pairs, label_key, read_marker_pairs};
use crate::annotate::rounds::ClusterId;
use anyhow::{Context, Result};
use enrichment::{annotate, AnnotateOutputs, GroupInputs, SpecificityMode};
use legume_numeric::matrix::dense_mat_io::Mat;
use log::info;
use nalgebra::{DMatrix, DVector};
use rand::SeedableRng;
use rand_distr::{Distribution, StandardNormal};
use serde::{Deserialize, Serialize};
use sha2::{Digest as _, Sha256};
use std::collections::BTreeMap;

/// A type stands, and a held-out half holds, below this q (or p).
const ALPHA: f32 = 0.05;
/// `lupin confirm`'s defaults, said once for the command and [`ConfirmArgs::new`].
const TAU: f64 = 1.0;
const SEED: u64 = 42;
const KNN: usize = 15;
const RESOLUTION: f64 = 1.0;
const MIN_CLUSTER_SIZE: usize = 10;
const MIN_MARKERS: usize = 3;
const FDR_ALPHA: f32 = 0.1;
/// What the confirming files end in, after the output prefix.
pub const CONFIRM_TSV: &str = ".confirm.tsv";
pub const CONFIRM_JSON: &str = ".confirm.json";
pub const CONFIRM_CLUSTERS: &str = ".confirm.clusters.parquet";

#[derive(clap::Args, Debug, Clone)]
pub struct ConfirmArgs {
    #[arg(
        long,
        short = 'f',
        help = "An annotated manifest (a round) or its prefix"
    )]
    pub from: Box<str>,

    #[arg(
        long,
        default_value_t = TAU,
        help = "Fission's τ: f = z + τw clusters, g = z − w/τ tests; 1 splits evenly"
    )]
    pub tau: f64,

    #[arg(
        long,
        value_enum,
        default_value_t = Sigma::Fitted,
        help = "The noise's covariance: fitted (so that f and g do not co-vary within the clusters of f), total (the cells' own: right only when no cell differs but by noise) or within (the round's clusters': right only when they are the cells' types)"
    )]
    pub sigma: Sigma,

    #[arg(
        long,
        default_value_t = SEED,
        help = "Seed of the noise, Leiden and the nulls"
    )]
    pub seed: u64,

    #[arg(
        long,
        default_value_t = KNN,
        help = "Leiden's k nearest neighbours on f"
    )]
    pub knn: usize,

    #[arg(long, default_value_t = RESOLUTION, help = "Leiden's resolution on f")]
    pub resolution: f64,

    #[arg(
        long = "num-clusters",
        help = "Aim Leiden on f at this many clusters (default: by --resolution alone; the round's own count was chosen on data that holds g)"
    )]
    pub num_clusters: Option<usize>,

    #[arg(
        long = "min-cluster-size",
        default_value_t = MIN_CLUSTER_SIZE,
        help = "Clusters of f with fewer cells are left out"
    )]
    pub min_cluster_size: usize,

    #[arg(long = "min-markers", default_value_t = MIN_MARKERS)]
    pub min_markers: usize,

    #[arg(long = "fdr-alpha", default_value_t = FDR_ALPHA)]
    pub fdr_alpha: f32,

    #[arg(long, help = "Print the result as JSON")]
    pub json: bool,
}

/// The covariance fission's noise is drawn with.
#[derive(clap::ValueEnum, Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq, Default)]
#[serde(rename_all = "lowercase")]
pub enum Sigma {
    /// Fitted so that `f` and `g` do not co-vary within the clusters of `f`
    /// ([`fitted_split`]), from the cells' covariance within the round's
    /// clusters.
    #[default]
    Fitted,
    /// The cells' covariance: right when no cell differs from another but by
    /// noise. Where some do, `Cov(f, g)` is minus the spread between them,
    /// and a cluster of `f` split by the noise shifts `g` the other way.
    Total,
    /// The cells' covariance within the round's clusters: right when they
    /// are the cells' types; where the round over-split, too small, and
    /// `Cov(f, g)` is positive.
    Within,
}

/// Fitting [`Sigma::Fitted`] stops once the cross-covariance is within this
/// many times what independent halves would show by chance
/// ([`cross_within`]), or after [`FIT_ROUNDS`].
const FIT_TOL: f64 = 1.5;
const FIT_ROUNDS: usize = 5;

impl ConfirmArgs {
    /// The command's defaults for the round `from`, as `lupin confirm -f
    /// from` takes them.
    pub fn new(from: &str) -> Self {
        Self {
            from: from.into(),
            tau: TAU,
            sigma: Sigma::default(),
            seed: SEED,
            knn: KNN,
            resolution: RESOLUTION,
            num_clusters: None,
            min_cluster_size: MIN_CLUSTER_SIZE,
            min_markers: MIN_MARKERS,
            fdr_alpha: FDR_ALPHA,
            json: false,
        }
    }

    /// Leiden's settings on `f`.
    fn leiden(&self) -> LeidenArgs {
        LeidenArgs {
            knn: self.knn,
            resolution: self.resolution,
            num_clusters: self.num_clusters,
            min_cluster_size: self.min_cluster_size,
            seed: Some(self.seed),
        }
    }
}

/// One cluster of `f`, its call on `g`, and the round's label it meets.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Confirmed {
    pub cluster: ClusterId,
    pub cells: usize,
    /// The round's cluster most of its cells are in, and their share.
    pub round_cluster: Option<ClusterId>,
    pub overlap: f32,
    pub round_label: Option<String>,
    /// The call on `g` with the whole panel: the type with the largest
    /// share of its Q, none when no type is significant.
    pub call: Option<String>,
    /// The call is the round's label, or one of the types a coarse label
    /// stands for.
    pub agrees: bool,
    /// Every type's test on `g` with the whole panel, by label key.
    pub types: BTreeMap<String, TypeTest>,
    /// The call each half of the markers makes, tested by the other half.
    pub held_out: [HeldOut; 2],
}

/// A type's enrichment in a cluster's ranking; `None` where it has none
/// (too few markers in the gene table).
#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Default)]
pub struct TypeTest {
    pub nes: Option<f32>,
    /// Against random gene sets of the same size.
    pub p: Option<f32>,
    pub q: Option<f32>,
    /// That the markers are up in the cluster: their mean logit against the
    /// centre, its variance from `Σ_g`, so markers that move together are
    /// counted as such, which random gene sets do not do.
    #[serde(default)]
    pub p_up: Option<f32>,
}

impl TypeTest {
    /// Below 0.05 on both tests, by its q against random gene sets.
    #[must_use]
    pub fn stands(&self) -> bool {
        Self::both(self.q, self.p_up)
    }

    /// [`Self::stands`] by its p: a held-out half is one test, not a panel's.
    #[must_use]
    pub fn stands_by_p(&self) -> bool {
        Self::both(self.p, self.p_up)
    }

    fn both(a: Option<f32>, b: Option<f32>) -> bool {
        a.is_some_and(|v| v < ALPHA) && b.is_some_and(|v| v < ALPHA)
    }
}

/// The call one half of the markers makes, and the other half's test of it.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct HeldOut {
    /// Made on the ranking from `f` with this half.
    pub call: Option<String>,
    /// The held-out half's p for that call on `g`, against random gene sets
    /// and against no cluster differing (none when it has too few markers
    /// to test).
    pub p: Option<f32>,
    #[serde(default)]
    pub p_up: Option<f32>,
    /// The call is the cluster's, and stands on `g` by the other half: both
    /// p below 0.05.
    pub holds: bool,
}

impl HeldOut {
    /// In a line: the call and its test by the other half.
    #[must_use]
    pub fn summary(&self) -> String {
        let num = |v: Option<f32>| v.map_or("-".into(), |v| format!("{v:.3}"));
        match (&self.call, self.p) {
            (Some(c), Some(p)) => format!(
                "{c} p {p:.3}, up p {}{}",
                num(self.p_up),
                if self.holds { "" } else { " ✗" }
            ),
            (Some(c), None) => format!("{c} untested"),
            (None, _) => "no call".into(),
        }
    }
}

/// What `lupin confirm` found, and how.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Confirmation {
    pub round: String,
    pub tau: f64,
    #[serde(default)]
    pub sigma: Sigma,
    /// What is left of `f` and `g`'s dependence: their covariance within the
    /// clusters of `f`, in units of what independent halves show by chance
    /// (about 1 when they are independent).
    #[serde(default)]
    pub dependence: f64,
    pub seed: u64,
    pub clusters: Vec<Confirmed>,
    /// Clusters of `f` whose call on `g` is the round's label.
    pub agree: usize,
    /// Clusters whose call held out in both directions: each half picked
    /// it, and the other half confirmed it.
    pub held: usize,
    pub files: BTreeMap<String, String>,
    /// The round's clusters file this was made from, by content: a
    /// confirmation of other clusters is not read as this round's.
    #[serde(default)]
    pub clusters_digest: String,
}

/// The round's clusters file, by content (a short SHA-256).
fn clusters_digest(loaded: &super::run::Loaded) -> Result<String> {
    let rel = (loaded.manifest.cluster.clusters.as_deref())
        .context("the round records no `cluster.clusters`")?;
    let path = super::run::resolve(&loaded.dir, rel);
    let bytes = std::fs::read(&path).with_context(|| format!("reading {path}"))?;
    Ok(crate::hex(&Sha256::digest(&bytes)[..8]))
}

/// The confirmation `lupin confirm` wrote beside the round `loaded`, when
/// there is one and it was made from the round's clusters as they are.
pub fn read(loaded: &super::run::Loaded) -> Result<Option<Confirmation>> {
    let path = format!("{}{CONFIRM_JSON}", loaded.run_prefix());
    let text = match std::fs::read_to_string(&path) {
        Ok(text) => text,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(e) => return Err(e).with_context(|| format!("reading {path}")),
    };
    let c: Confirmation = serde_json::from_str(&text).with_context(|| format!("reading {path}"))?;
    if clusters_digest(loaded).ok().as_deref() != Some(c.clusters_digest.as_str()) {
        log::warn!("{path} was made from other clusters than the round's: confirm it again");
        return Ok(None);
    }
    Ok(Some(c))
}

/// `lupin confirm`.
pub fn run_confirm(args: &ConfirmArgs) -> Result<()> {
    let c = confirm(args)?;
    if args.json {
        println!("{}", serde_json::to_string_pretty(&c)?);
    } else {
        print!("{}", render(&c));
    }
    Ok(())
}

/// Split the round's embedding, cluster one half and test the round's
/// panel on the other; write the table, its record and the clusters of `f`.
pub fn confirm(args: &ConfirmArgs) -> Result<Confirmation> {
    anyhow::ensure!(
        args.tau.is_finite() && args.tau > 0.0,
        "--tau must be positive"
    );
    let loaded = super::run::load(&args.from)?;
    let markers = (loaded.manifest.annotate.markers.as_deref())
        .map(|r| super::run::resolve(&loaded.dir, r))
        .context("the round records no marker panel (`annotate.markers`)")?;
    let (feat, emb) = score_tables(&loaded)?;
    anyhow::ensure!(
        feat.mat.ncols() == emb.mat.ncols(),
        "the gene table has {} dimensions but the cell embedding {}",
        feat.mat.ncols(),
        emb.mat.ncols()
    );
    let clusters_rel = (loaded.manifest.cluster.clusters.as_deref())
        .context("the round records no `cluster.clusters`")?;
    let (names, ids) = read_clusters(&super::run::resolve(&loaded.dir, clusters_rel))?;
    let round = by_name(&emb.rows, &names, &ids);
    // What the report reads of the round, before the split's Leiden runs.
    let (labels, tree) = labels_and_tree(&loaded)?;
    let pairs = read_marker_pairs(&markers)?;
    let panel: BTreeMap<String, ()> = pairs.iter().map(|(_, t)| (label_key(t), ())).collect();

    let z = emb.mat.map(f64::from);
    let space = loaded.manifest.kind.cell_space();
    let Split {
        f,
        g,
        groups,
        k,
        dependence,
    } = split(&z, &round, space, args)?;
    anyhow::ensure!(
        k >= 2,
        "Leiden on f gave {k} cluster(s); a call needs at least 2"
    );
    info!(
        "fission: {k} clusters of f, tested on g (τ = {}, dependence {dependence:.3})",
        args.tau
    );

    // The tests read g's ranking; the held-out halves pick their call on
    // f's, so the half that tests a call did not choose it on the same data.
    let rho = feat.mat.map(f64::from);
    let (on_g, on_f) = rayon::join(
        || Ranking::of(&rho, &g, &groups, k),
        || Ranking::of(&rho, &f, &groups, k),
    );
    let (on_g, on_f) = (on_g?, on_f?);
    let halves = halves(&pairs);
    // The five panels' tests are independent: run side by side.
    let runs: Vec<(&Ranking, &[Pair])> = vec![
        (&on_g, &pairs),
        (&on_g, &halves[0]),
        (&on_g, &halves[1]),
        (&on_f, &halves[0]),
        (&on_f, &halves[1]),
    ];
    let mut scored = {
        use rayon::prelude::*;
        (runs.into_par_iter())
            .map(|(r, p)| test_panel(r, p, &feat.rows, args))
            .collect::<Result<Vec<Scored>>>()?
            .into_iter()
    };
    let mut next = || scored.next().expect("five panels scored");
    let whole = next();
    let tested = [next(), next()];
    let picked = [next(), next()];

    let agrees = |label: &str, call: &str| {
        let types = (tree.as_ref())
            .and_then(|t| t.coarse_members(label, &panel))
            .map_or_else(|| vec![label.to_string()], <[String]>::to_vec);
        types.iter().any(|t| label_key(t) == label_key(call))
    };
    let overlap = overlap(&groups, &round, k);
    let mut clusters = Vec::with_capacity(k);
    for (kk, (cells, round_cluster, share)) in overlap.into_iter().enumerate() {
        let round_label = round_cluster.and_then(|r| labels.get(&r).cloned().flatten());
        let call = whole.call(kk);
        // A half holds when its call (made on f) is the cluster's call and
        // stands on g by the other half.
        let held = |by: &Scored, on: &Scored| {
            let picked = by.call(kk);
            let t = picked.as_deref().and_then(|c| on.test(kk, c));
            let same =
                matches!((&picked, &call), (Some(a), Some(b)) if label_key(a) == label_key(b));
            HeldOut {
                holds: same && t.is_some_and(|t| t.stands_by_p()),
                call: picked,
                p: t.and_then(|t| t.p),
                p_up: t.and_then(|t| t.p_up),
            }
        };
        let held_out = [held(&picked[0], &tested[1]), held(&picked[1], &tested[0])];
        clusters.push(Confirmed {
            cluster: kk as ClusterId,
            cells,
            round_cluster,
            overlap: share,
            agrees: matches!((&round_label, &call), (Some(l), Some(c)) if agrees(l, c)),
            round_label,
            call,
            types: whole.tests(kk),
            held_out,
        });
    }

    // Beside the round, where `lupin why` reads it, while no other writer
    // moves the round's files.
    let out = loaded.run_prefix();
    let _lock = super::rounds::Chain::of(&loaded.file).lock()?;
    let mut files = BTreeMap::new();
    let path = format!("{out}{CONFIRM_CLUSTERS}");
    super::rounds::write_clusters(&path, &emb.rows, &super::rounds::cluster_ids(&groups))?;
    files.insert("clusters".into(), path);
    let path = format!("{out}{CONFIRM_TSV}");
    std::fs::write(&path, tsv(&clusters)).with_context(|| format!("writing {path}"))?;
    files.insert("table".into(), path);
    let c = Confirmation {
        round: loaded.file.display().to_string(),
        tau: args.tau,
        sigma: args.sigma,
        dependence,
        seed: args.seed,
        agree: clusters.iter().filter(|c| c.agrees).count(),
        held: (clusters.iter())
            .filter(|c| c.held_out.iter().all(|h| h.holds))
            .count(),
        clusters,
        files,
        clusters_digest: clusters_digest(&loaded)?,
    };
    let path = format!("{out}{CONFIRM_JSON}");
    std::fs::write(&path, serde_json::to_string_pretty(&c)?)
        .with_context(|| format!("writing {path}"))?;
    info!("wrote {path}");
    Ok(c)
}

impl Confirmed {
    /// The whole panel's test of the call.
    #[must_use]
    pub fn call_test(&self) -> Option<TypeTest> {
        self.types.get(&label_key(self.call.as_deref()?)).copied()
    }
}

/// Each of `rows`' round cluster, found by cell name in `names`.
fn by_name(
    rows: &[Box<str>],
    names: &[Box<str>],
    clusters: &[Option<ClusterId>],
) -> Vec<Option<ClusterId>> {
    let at: BTreeMap<&str, Option<ClusterId>> = names
        .iter()
        .map(AsRef::as_ref)
        .zip(clusters.iter().copied())
        .collect();
    rows.iter()
        .map(|r| at.get(r.as_ref()).copied().flatten())
        .collect()
}

/// A marker and its cell type, as the panel reader gives them.
type Pair = (Box<str>, Box<str>);

/// Each cluster's size and mean, by cluster.
type Means<G> = BTreeMap<G, (usize, DVector<f64>)>;

/// Each cluster's size and mean of `x` (`[N × H]`) over the cells of
/// `groups` (`usize::MAX`, or `None`, left out), and the covariance within
/// them, pooled (ddof: one per cluster).
fn within<G: Copy + Ord>(
    x: &DMatrix<f64>,
    groups: impl Iterator<Item = Option<G>> + Clone,
) -> Result<(Means<G>, DMatrix<f64>)> {
    let h = x.ncols();
    let mut by: Means<G> = BTreeMap::new();
    for (row, gr) in x.row_iter().zip(groups.clone()) {
        if let Some(gr) = gr {
            let e = by.entry(gr).or_insert_with(|| (0, DVector::zeros(h)));
            e.0 += 1;
            e.1 += row.transpose();
        }
    }
    for (n, s) in by.values_mut() {
        *s /= *n as f64;
    }
    let (centred, dof) = centred(x, groups, &by);
    anyhow::ensure!(
        dof > 0,
        "too few cells in clusters for a covariance within them"
    );
    Ok((by, centred.tr_mul(&centred) / dof as f64))
}

/// The rows of `x` in a cluster of `groups`, less their cluster's mean in
/// `by` (`[n × H]`), and the degrees of freedom left (one spent per
/// cluster).
fn centred<G: Copy + Ord>(
    x: &DMatrix<f64>,
    groups: impl Iterator<Item = Option<G>>,
    by: &Means<G>,
) -> (DMatrix<f64>, usize) {
    let rows: Vec<(usize, &DVector<f64>)> = (groups.enumerate())
        .filter_map(|(i, gr)| Some((i, &by.get(&gr?)?.1)))
        .collect();
    let out = DMatrix::from_fn(rows.len(), x.ncols(), |r, j| {
        let (i, m) = rows[r];
        x[(i, j)] - m[j]
    });
    (out, rows.len().saturating_sub(by.len()))
}

/// Each cell's cluster, `None` for the cells left out (`usize::MAX`).
fn assigned(groups: &[usize]) -> impl Iterator<Item = Option<usize>> + Clone + '_ {
    groups.iter().map(|&c| (c != usize::MAX).then_some(c))
}

/// The two halves, the clusters of `f`, and what is left of their
/// dependence.
struct Split {
    f: DMatrix<f64>,
    g: DMatrix<f64>,
    groups: Vec<usize>,
    /// The clusters of `f`.
    k: usize,
    dependence: f64,
}

/// Split `z` by `args.sigma` and cluster `f`.
fn split(
    z: &DMatrix<f64>,
    round: &[Option<ClusterId>],
    space: CellSpace,
    args: &ConfirmArgs,
) -> Result<Split> {
    let all = vec![Some(0); z.nrows()];
    // The fit starts where the round's clusters put the noise: under no
    // difference it corrects that in a round; from the cells' total spread
    // f would be noisy, its clusters mixed, and their spread of means read
    // as dependence.
    let start = match args.sigma {
        Sigma::Within | Sigma::Fitted => within(z, round.iter().copied())?.1,
        Sigma::Total => within(z, all.iter().copied())?.1,
    };
    let rounds = if args.sigma == Sigma::Fitted {
        FIT_ROUNDS
    } else {
        1
    };
    fitted_split(z, start, rounds, space, args)
}

/// Fission with `Σ` fitted to the condition that holds when it is right:
/// `f` and `g` independent, so not co-varying within any cluster of `f`.
/// Their pooled cross-covariance within the clusters, `C`, is zero then;
/// otherwise each round adds the shortfall `Σ_noise − Σ` it implies
/// ([`shortfall`]) to `Σ` (kept positive semi-definite), splits
/// again with the same draw of noise and clusters `f` again, until `C` is
/// within [`FIT_TOL`] of its chance size or `rounds` are spent. The last
/// split is returned, with that ratio.
fn fitted_split(
    z: &DMatrix<f64>,
    mut sigma: DMatrix<f64>,
    rounds: usize,
    space: CellSpace,
    args: &ConfirmArgs,
) -> Result<Split> {
    let tau = args.tau;
    let eps = noise(z.nrows(), z.ncols(), args.seed);
    let mut at = 0;
    loop {
        let (f, g) = fission(z, &eps, &sigma, tau);
        let mut x: Mat = f.map(|v| v as f32);
        prepare_cells(space, &mut x);
        let (groups, k) = leiden_on(&x, space, &args.leiden())?;
        let cross = cross_within(&f, &g, &groups)?;
        let sym = (&cross.fg + cross.fg.transpose()) / 2.0;
        let dependence = sym.norm() / cross.chance.max(f64::MIN_POSITIVE);
        at += 1;
        info!("fission round {at}: dependence {dependence:.3}");
        if dependence <= FIT_TOL || at >= rounds {
            return Ok(Split {
                f,
                g,
                groups,
                k,
                dependence,
            });
        }
        let d = shortfall(&cross, &sigma, tau).unwrap_or(sym);
        sigma = psd_map(&(sigma + d), |l| l);
    }
}

/// `f` and `g` within the clusters of `f`.
struct Cross {
    /// `[H × H]`: `Cov(f, g)`, pooled; ddof one per cluster.
    fg: DMatrix<f64>,
    /// `Cov(f)`, likewise.
    ff: DMatrix<f64>,
    /// The Frobenius norm `fg`'s symmetric part has by chance when `f` and
    /// `g` are independent: entry `(a, b)` then has variance
    /// `(F_aa G_bb + F_bb G_aa + 2 F_ab G_ab) / 4n`, `F` and `G` their own
    /// covariances within the clusters.
    chance: f64,
}

/// The shortfall `D = Σ_noise − Σ` that `cross` implies. Within a cluster
/// of `f` the spread of `f` is cut short by the clustering, and their
/// covariance with it; the slope of `g` on `f` is not: `M = Cov(g, f)
/// Cov(f)⁻¹` estimates `D V⁻¹`, `V = D + (1 + τ²) Σ` being `f`'s spread
/// about the cells' means, so `D = (I − M)⁻¹ (1 + τ²) M Σ` (its symmetric
/// part). `None` when `Cov(f)` or `I − M` cannot be inverted.
fn shortfall(cross: &Cross, sigma: &DMatrix<f64>, tau: f64) -> Option<DMatrix<f64>> {
    let m = cross.fg.transpose() * cross.ff.clone().try_inverse()?;
    let i = DMatrix::<f64>::identity(m.nrows(), m.ncols());
    let d = (&i - &m).try_inverse()? * m * sigma * (1.0 + tau * tau);
    Some((&d + d.transpose()) / 2.0)
}

/// [`Cross`] for `f` and `g` within the clusters of `groups`.
fn cross_within(f: &DMatrix<f64>, g: &DMatrix<f64>, groups: &[usize]) -> Result<Cross> {
    let (fm, fc) = within(f, assigned(groups))?;
    let (gm, gc) = within(g, assigned(groups))?;
    let ((cf, dof), (cg, _)) = (
        centred(f, assigned(groups), &fm),
        centred(g, assigned(groups), &gm),
    );
    let h = f.ncols();
    let dof = dof.max(1) as f64;
    let c = cf.tr_mul(&cg) / dof;
    let chance = (0..h)
        .flat_map(|a| (0..h).map(move |b| (a, b)))
        .map(|(a, b)| {
            (fc[(a, a)] * gc[(b, b)] + fc[(b, b)] * gc[(a, a)] + 2.0 * fc[(a, b)] * gc[(a, b)])
                / (4.0 * dof)
        })
        .sum::<f64>()
        .sqrt();
    Ok(Cross {
        fg: c,
        ff: fc,
        chance,
    })
}

/// `[n × h]` standard normals drawn from `seed`: the noise before `Σ`.
fn noise(n: usize, h: usize, seed: u64) -> DMatrix<f64> {
    let mut rng = rand::rngs::SmallRng::seed_from_u64(seed);
    DMatrix::<f64>::from_fn(n, h, |_, _| StandardNormal.sample(&mut rng))
}

/// `f = z + τw` and `g = z − w/τ`, `w = Σ^{1/2} ε` per cell: `N(0, Σ)`
/// for standard normal `eps` ([`noise`]).
fn fission(
    z: &DMatrix<f64>,
    eps: &DMatrix<f64>,
    sigma: &DMatrix<f64>,
    tau: f64,
) -> (DMatrix<f64>, DMatrix<f64>) {
    let w = eps * psd_map(sigma, f64::sqrt);
    (z + &w * tau, z - &w / tau)
}

/// Each cluster's ranking of the genes from one half, and what the test of
/// its markers' mean reads.
struct Ranking<'a> {
    rho: &'a DMatrix<f64>,
    /// `[G × K]`: each gene's mean logit in each cluster against the
    /// balanced centre, in standard errors (module docs).
    t: Mat,
    /// `[H × K]`: each cluster's mean less the balanced centre.
    d: DMatrix<f64>,
    /// The half's covariance within the clusters.
    sigma: DMatrix<f64>,
    /// `√(ρ_gᵀ Σ ρ_g)` per gene.
    sd: DVector<f64>,
    size: Vec<usize>,
}

impl<'a> Ranking<'a> {
    fn of(rho: &'a DMatrix<f64>, x: &DMatrix<f64>, groups: &[usize], k: usize) -> Result<Self> {
        let (by, sigma) = within(x, assigned(groups))?;
        let mut quad = rho * &sigma;
        quad.component_mul_assign(rho);
        let mut d = DMatrix::<f64>::zeros(x.ncols(), k);
        let mut size = vec![0usize; k];
        for (&c, (n, m)) in &by {
            d.set_column(c, m);
            size[c] = *n;
        }
        let centre = d.column_mean();
        for mut col in d.column_iter_mut() {
            col -= &centre;
        }
        let sd = quad.column_sum().map(|v| v.max(0.0).sqrt());
        let contrast = rho * &d;
        // Every cluster of `f` has cells: Leiden's are numbered densely.
        let t = Mat::from_fn(rho.nrows(), k, |gi, c| {
            let se = sd[gi] / (size[c] as f64).sqrt();
            if se > 0.0 {
                (contrast[(gi, c)] / se) as f32
            } else {
                0.0
            }
        });
        Ok(Self {
            rho,
            t,
            d,
            sigma,
            sd,
            size,
        })
    }

    /// `[K × C]` one-sided p that each type's markers are up in each
    /// cluster: their weighted mean standardised logit against the centre,
    /// `u_T · d_k` with `u_T = Σ_g w_g ρ_g / sd_g / Σ_g w_g`, over its
    /// standard error `√(u_Tᵀ Σ u_T · v_k)`. `v_k = (1 − 1/K)²/n_k +
    /// Σ_{j≠k} 1/(K² n_j)` is the variance of a cluster's mean less the
    /// mean of all the cluster means, per unit of `Σ`. The markers'
    /// correlation is in `Σ`, as random gene sets do not have it. NaN for a
    /// type with no hit weight.
    fn up_p(&self, weights: &Mat) -> Mat {
        use crate::annotate::geneset_units::normal_upper_tail;
        let (h, k) = self.d.shape();
        let kf = k as f64;
        let inv: Vec<f64> = self.size.iter().map(|&n| 1.0 / n as f64).collect();
        let total: f64 = inv.iter().sum();
        let v: Vec<f64> = (0..k)
            .map(|kk| (1.0 - 1.0 / kf).powi(2) * inv[kk] + (total - inv[kk]) / (kf * kf))
            .collect();
        // Each type's direction and its variance once, then every cluster.
        let c = weights.ncols();
        let mut u = DMatrix::<f64>::zeros(h, c);
        for (cc, mut col) in u.column_iter_mut().enumerate() {
            let mut mass = 0.0;
            for g in 0..weights.nrows() {
                let w = f64::from(weights[(g, cc)]);
                if w > 0.0 && self.sd[g] > 0.0 {
                    col += self.rho.row(g).transpose() * (w / self.sd[g]);
                    mass += w;
                }
            }
            col /= if mass > 0.0 { mass } else { f64::NAN };
        }
        let spread = (&self.sigma * &u).component_mul(&u).row_sum();
        let shift = u.tr_mul(&self.d);
        Mat::from_fn(k, c, |kk, cc| {
            let se = (spread[cc] * v[kk]).sqrt();
            if se > 0.0 {
                normal_upper_tail((shift[(cc, kk)] / se) as f32)
            } else {
                f32::NAN
            }
        })
    }

    /// The enrichment's inputs for this ranking, per cluster: no
    /// pseudobulk, as no sample permutation is run, and each cluster its own
    /// "cell", as only the cluster × type tables are read.
    fn inputs(&self, genes: &[Box<str>]) -> GroupInputs {
        let k = self.size.len();
        GroupInputs {
            pb_gene_gp: Mat::zeros(genes.len(), 1),
            pb_membership_pk: Mat::zeros(1, k),
            profile_gk: self.t.clone(),
            cell_membership_nk: Mat::identity(k, k),
            gene_names: genes.to_vec(),
            cell_names: (0..k).map(|c| format!("F{c}").into()).collect(),
        }
    }
}

/// `V f(Λ) Vᵀ` for symmetric `Σ = V Λ Vᵀ`, negative eigenvalues (from
/// rounding, or a fitting step) read as 0: `f = sqrt` gives `Σ^{1/2}`, so
/// `w = ε Σ^{1/2}` is `N(0, Σ)` for standard normal rows `ε`; the identity
/// projects `Σ` onto the positive semi-definite matrices.
fn psd_map(sigma: &DMatrix<f64>, f: impl Fn(f64) -> f64) -> DMatrix<f64> {
    let eig = sigma.clone().symmetric_eigen();
    &eig.eigenvectors
        * DMatrix::from_diagonal(&eig.eigenvalues.map(|l| f(l.max(0.0))))
        * eig.eigenvectors.transpose()
}

/// The panel `pairs` tested on the rankings of `group`.
struct Scored {
    names: Vec<Box<str>>,
    out: AnnotateOutputs,
    /// `[K × C]` p that the type's markers are up ([`Ranking::up_p`]).
    up: Mat,
    /// The types the enrichment tested (enough markers, not every gene).
    live: Vec<bool>,
}

impl Scored {
    /// Cluster `k`'s call: of the types that stand on both tests, the one
    /// with the largest share of its Q; none when no type stands.
    fn call(&self, k: usize) -> Option<String> {
        let row = self.out.q_kc.row(k);
        (0..row.len())
            .filter(|&c| self.at(k, c).stands())
            .max_by(|&a, &b| row[a].total_cmp(&row[b]))
            .map(|c| self.names[c].to_string())
    }

    /// Type `c`'s test in cluster `k`; nothing for a type not tested.
    fn at(&self, k: usize, c: usize) -> TypeTest {
        if !self.live[c] {
            return TypeTest::default();
        }
        let finite = |v: f32| v.is_finite().then_some(v);
        TypeTest {
            nes: finite(self.out.nes_kc[(k, c)]),
            p: finite(self.out.pvalue_kc[(k, c)]),
            q: finite(self.out.qvalue_kc[(k, c)]),
            p_up: finite(self.up[(k, c)]),
        }
    }

    /// Type `name`'s test in cluster `k`, when this panel has it.
    fn test(&self, k: usize, name: &str) -> Option<TypeTest> {
        let key = label_key(name);
        (self.names.iter())
            .position(|n| n.as_ref() == key)
            .map(|c| self.at(k, c))
    }

    /// Every type's test in cluster `k`, by label key.
    fn tests(&self, k: usize) -> BTreeMap<String, TypeTest> {
        (self.names.iter().enumerate())
            .map(|(c, n)| (n.to_string(), self.at(k, c)))
            .collect()
    }
}

/// The panel `pairs` tested on ranking `r` of the clusters `groups`: by
/// the enrichment against random gene sets, and by its markers' mean.
fn test_panel(
    r: &Ranking,
    pairs: &[Pair],
    genes: &[Box<str>],
    args: &ConfirmArgs,
) -> Result<Scored> {
    let annot = annotation_matrix_from_pairs(pairs, genes)?;
    let mut config = crate::annotate::by_enrichment::null_config(
        args.seed,
        args.fdr_alpha,
        1.0,
        args.min_markers,
        None,
    );
    // The ranking is already a contrast between the clusters, on a scale
    // of standard errors: no specificity transform, and no abundance
    // strata (an unmeasured gene sits at 0, mid-ranking, not at the foot).
    config.specificity = SpecificityMode::Raw;
    config.stratify_null = false;
    let group = r.inputs(genes);
    let out = annotate(&group, &annot.membership_ga, &annot.annot_names, &config)?;
    // Weighted as the enrichment weighs them, the types it leaves out left
    // out.
    let live = enrichment::live_types(&annot.membership_ga, &config);
    let mut weights = annot.membership_ga.clone();
    for (mut col, &on) in weights.column_iter_mut().zip(&live) {
        if !on {
            col.fill(0.0);
        }
    }
    Ok(Scored {
        names: annot.annot_names,
        out,
        up: r.up_p(&weights),
        live,
    })
}

/// Each type's markers in two halves: ordered by a hash of the gene name,
/// alternately.
fn halves(pairs: &[Pair]) -> [Vec<Pair>; 2] {
    let mut by: BTreeMap<String, Vec<&Box<str>>> = BTreeMap::new();
    for (gene, t) in pairs {
        by.entry(label_key(t)).or_default().push(gene);
    }
    let mut out = [Vec::new(), Vec::new()];
    for (t, mut genes) in by {
        genes.sort_by_cached_key(|g| Sha256::digest(g.as_bytes()));
        genes.dedup();
        for (i, g) in genes.into_iter().enumerate() {
            out[i % 2].push((g.clone(), t.as_str().into()));
        }
    }
    out
}

/// Per cluster of `groups`: its cells, the round cluster most of them are
/// in and their share.
fn overlap(
    groups: &[usize],
    round: &[Option<ClusterId>],
    k: usize,
) -> Vec<(usize, Option<ClusterId>, f32)> {
    let mut counts: Vec<BTreeMap<ClusterId, usize>> = vec![BTreeMap::new(); k];
    let mut cells = vec![0usize; k];
    for (&c, r) in groups.iter().zip(round) {
        if c < k {
            cells[c] += 1;
            if let Some(r) = r {
                *counts[c].entry(*r).or_default() += 1;
            }
        }
    }
    (counts.into_iter().zip(cells))
        .map(|(m, n)| {
            let best = m
                .into_iter()
                .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0)));
            match best {
                Some((r, c)) => (n, Some(r), c as f32 / n as f32),
                None => (n, None, 0.0),
            }
        })
        .collect()
}

fn tsv(rows: &[Confirmed]) -> String {
    let mut s = String::from(
        "cluster\tcells\tround_cluster\toverlap\tround_label\tcall\tnes\tp\tq\tp_up\tagrees\tcall_a\tp_b\tp_up_b\tcall_b\tp_a\tp_up_a\n",
    );
    let opt = |o: &Option<String>| o.clone().unwrap_or_default();
    let num = |v: Option<f32>| v.map(|v| v.to_string()).unwrap_or_default();
    for r in rows {
        let [a, b] = &r.held_out;
        let t = r.call_test().unwrap_or_default();
        s.push_str(&format!(
            "{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.cluster,
            r.cells,
            r.round_cluster.map(|c| c.to_string()).unwrap_or_default(),
            r.overlap,
            opt(&r.round_label),
            opt(&r.call),
            num(t.nes),
            num(t.p),
            num(t.q),
            num(t.p_up),
            r.agrees,
            opt(&a.call),
            num(a.p),
            num(a.p_up),
            opt(&b.call),
            num(b.p),
            num(b.p_up),
        ));
    }
    s
}

/// The result as text: a line per cluster of `f`, then the tally.
pub fn render(c: &Confirmation) -> String {
    let mut s = format!(
        "Fission of {} (τ = {}): clusters of f, called on g\n",
        c.round, c.tau
    );
    for r in &c.clusters {
        let round = match (r.round_cluster, &r.round_label) {
            (Some(k), Some(l)) => format!("K{k} {l} ({:.0}%)", 100.0 * r.overlap),
            (Some(k), None) => format!("K{k} ({:.0}%)", 100.0 * r.overlap),
            _ => "-".into(),
        };
        let t = r.call_test();
        let num = |v: Option<f32>| v.map_or("-".into(), |v| format!("{v:.3}"));
        s.push_str(&format!(
            "  F{} {} cells ~ {round}: {} (NES {}, q {}, up p {}){}; held out: {} / {}\n",
            r.cluster,
            r.cells,
            r.call.as_deref().unwrap_or("no call"),
            num(t.and_then(|t| t.nes)),
            num(t.and_then(|t| t.q)),
            num(t.and_then(|t| t.p_up)),
            if r.agrees { "" } else { " ≠ round" },
            r.held_out[0].summary(),
            r.held_out[1].summary(),
        ));
    }
    s.push_str(&format!(
        "{} of {} clusters of f called as the round labels them; {} held out both ways (p < {ALPHA})\n",
        c.agree,
        c.clusters.len(),
        c.held
    ));
    s
}

#[cfg(test)]
#[path = "tests/confirm.rs"]
pub(crate) mod tests;
