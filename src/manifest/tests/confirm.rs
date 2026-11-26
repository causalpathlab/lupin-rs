//! [`super`] on a planted run: five clusters, each along its own axis of
//! the embedding, the first four with a type whose markers are the genes
//! read most along it; the fifth is no type of the panel's.

use super::*;
use crate::manifest::annotate::annotate_by_marker_score;
use crate::manifest::run::{self, annotated_path};
use legume_numeric::matrix::traits::IoOps;
use rand::rngs::SmallRng;

const H: usize = 6;
const K: usize = 4;
/// The planted clusters: the typed ones and one of no type.
const CLUSTERS: usize = K + 1;
const PER: usize = 80;
const GENES: usize = 240;
const MARKERS: usize = 10;
const TYPES: [&str; K] = ["A", "B", "C", "D"];

/// The run in `dir`, each type with `markers[t]` markers: its manifest's
/// path and its panel's.
fn planted(dir: &std::path::Path, markers: [usize; K]) -> (String, String) {
    let mut rng = SmallRng::seed_from_u64(7);
    let normal = |rng: &mut SmallRng| -> f32 { StandardNormal.sample(rng) };
    let rho = Mat::from_fn(GENES, H, |_, _| normal(&mut rng));
    let n = CLUSTERS * PER;
    let z = Mat::from_fn(n, H, |i, h| {
        let k = i / PER;
        (if h == k { 3.0 } else { 0.0 }) + 0.5 * normal(&mut rng)
    });
    let genes: Vec<Box<str>> = (0..GENES).map(|g| format!("g{g}").into()).collect();
    let cells: Vec<Box<str>> = (0..n).map(|i| format!("c{i:03}").into()).collect();
    let at = |name: &str| dir.join(name).to_string_lossy().into_owned();
    rho.to_parquet_with_names(
        &at("run.feature_embedding.parquet"),
        (Some(&genes), Some("gene")),
        None,
    )
    .unwrap();
    z.to_parquet_with_names(
        &at("run.cell_embedding.parquet"),
        (Some(&cells), Some("cell")),
        None,
    )
    .unwrap();
    let groups: Vec<usize> = (0..n).map(|i| i / PER).collect();
    crate::manifest::rounds::write_clusters(
        &at("run.clusters.parquet"),
        &cells,
        &crate::manifest::rounds::cluster_ids(&groups),
    )
    .unwrap();
    // Each type's markers: the genes read most along its cluster's axis.
    let mut panel = String::new();
    for (k, t) in TYPES.iter().enumerate() {
        let mut by: Vec<usize> = (0..GENES).collect();
        by.sort_by(|&a, &b| rho[(b, k)].total_cmp(&rho[(a, k)]));
        for &g in &by[..markers[k]] {
            panel.push_str(&format!("g{g}\t{t}\n"));
        }
    }
    std::fs::write(at("panel.tsv"), panel).unwrap();
    let manifest = serde_json::json!({
        "version": 2, "kind": "bge", "prefix": "run",
        "outputs": {"feature_embedding": "run.feature_embedding.parquet",
                    "cell_embedding": "run.cell_embedding.parquet"},
        "cluster": {"clusters": "run.clusters.parquet"}
    });
    std::fs::write(at("run.senna.json"), manifest.to_string()).unwrap();
    (at("run.senna.json"), at("panel.tsv"))
}

/// The planted run's score round in `dir`: its manifest's path.
pub(crate) fn round(dir: &std::path::Path) -> String {
    round_with(dir, [MARKERS; K])
}

/// [`round`] with `markers[t]` markers for each type.
fn round_with(dir: &std::path::Path, markers: [usize; K]) -> String {
    let (manifest, panel) = planted(dir, markers);
    let loaded = run::load(&manifest).unwrap();
    let out = dir.join("run.L1").to_string_lossy().into_owned();
    annotate_by_marker_score(
        &crate::annotate::marker_score::toy::args(&panel, &out),
        &loaded,
    )
    .unwrap();
    annotated_path(&loaded.file, &out)
        .to_string_lossy()
        .into_owned()
}

fn args(from: &str) -> ConfirmArgs {
    ConfirmArgs::new(from)
}

#[test]
fn the_planted_types_are_called_on_the_half_the_clusters_were_not_made_from() {
    let dir = tempfile::tempdir().unwrap();
    let from = round(dir.path());
    let c = confirm(&args(&from)).unwrap();
    assert!(
        c.dependence < 3.0,
        "f and g still co-vary: {}",
        c.dependence
    );
    assert!(
        c.clusters.len() >= CLUSTERS,
        "{} clusters of f",
        c.clusters.len()
    );
    let (typed, typeless): (Vec<_>, Vec<_>) =
        (c.clusters.iter()).partition(|x| x.round_cluster.is_some_and(|r| (r as usize) < K));
    // The cluster of no type: no type of the panel stands in it.
    for x in &typeless {
        let stand: Vec<_> = (x.types.iter())
            .filter(|(_, t)| t.stands())
            .map(|(n, _)| n)
            .collect();
        assert!(
            stand.is_empty(),
            "F{} of no type: {stand:?} stand",
            x.cluster
        );
    }
    for x in &typed {
        assert!(
            x.agrees,
            "F{}: {:?} against {:?}",
            x.cluster, x.call, x.round_label
        );
        assert!(x.overlap > 0.9, "F{} overlaps by {}", x.cluster, x.overlap);
        for h in &x.held_out {
            assert!(h.holds, "F{}: {h:?}", x.cluster);
            assert_eq!(h.call, x.call, "F{}: a half calls another type", x.cluster);
        }
    }
    assert_eq!(c.agree, typed.len());
    assert_eq!(c.held, typed.len());
    for f in c.files.values() {
        assert!(std::path::Path::new(f).is_file(), "{f}");
    }
    let tsv = std::fs::read_to_string(&c.files["table"]).unwrap();
    assert_eq!(tsv.lines().count(), c.clusters.len() + 1);
}

#[test]
fn the_noise_is_drawn_again_for_another_seed_and_alike_for_the_same() {
    let z = DMatrix::<f64>::from_fn(40, 3, |i, h| ((i * 7 + h * 3) % 11) as f64);
    let round: Vec<Option<ClusterId>> = (0..40).map(|i| Some((i / 20) as ClusterId)).collect();
    let sigma = within(&z, round.iter().copied()).unwrap().1;
    let (f1, g1) = fission(&z, &noise(40, 3, 1), &sigma, 1.0);
    let (f2, _) = fission(&z, &noise(40, 3, 1), &sigma, 1.0);
    let (f3, _) = fission(&z, &noise(40, 3, 2), &sigma, 1.0);
    assert_eq!(f1, f2);
    assert_ne!(f1, f3);
    // z is the precision-weighted mean of the halves: (f + τ²g) / (1 + τ²).
    let back = (&f1 + &g1) / 2.0;
    assert!((back - &z).abs().max() < 1e-9);
}

#[test]
fn a_types_markers_are_halved_and_none_is_lost() {
    let pairs: Vec<(Box<str>, Box<str>)> = (0..7)
        .map(|g| (format!("g{g}").into(), "T cell".into()))
        .chain([("x".into(), "B".into())])
        .collect();
    let [a, b] = halves(&pairs);
    assert_eq!(a.len() + b.len(), pairs.len());
    let of = |h: &[(Box<str>, Box<str>)], t: &str| h.iter().filter(|p| p.1.as_ref() == t).count();
    assert_eq!((of(&a, "T_cell"), of(&b, "T_cell")), (4, 3));
    assert!(a.iter().all(|p| !b.contains(p)));
}

#[test]
fn why_weighs_the_fission_line_beside_the_others() {
    let dir = tempfile::tempdir().unwrap();
    let from = round(dir.path());
    let loaded = run::load(&from).unwrap();
    let r = crate::manifest::why::Round::read(&loaded).unwrap();
    let before = crate::manifest::why::why_label(&r, 0, None, Some("B"), 1).unwrap();
    assert!(before
        .children
        .iter()
        .all(|l| l.line.as_deref() != Some("fission")));

    confirm(&args(&from)).unwrap();
    let r = crate::manifest::why::Round::read(&loaded).unwrap();
    let a = crate::manifest::why::why_label(&r, 0, None, Some("B"), 1).unwrap();
    let line = (a.children.iter())
        .find(|l| l.line.as_deref() == Some("fission"))
        .expect("a fission line");
    assert_eq!(
        line.favours,
        crate::manifest::why::Favours::Fact,
        "{line:?}"
    );
    assert!(
        line.notes.iter().any(|n| n.contains("held-out")),
        "{:?}",
        line.notes
    );
    assert_eq!(a.children.len(), before.children.len() + 1);
}

/// One homogeneous population (no cluster, no type): every call is false.
/// Clustered and tested on the same `z`, the double dip; or clustered on
/// `f` and tested on `g`. The rate of (cluster, type) pairs whose markers
/// stand on both tests, and of those whose markers' mean is up at 0.05.
fn false_calls(fissioned: bool, seed: u64) -> (f64, f64) {
    let mut rng = SmallRng::seed_from_u64(seed);
    let normal = |rng: &mut SmallRng| -> f64 { StandardNormal.sample(rng) };
    let (n, genes) = (1500, GENES);
    let rho = DMatrix::<f64>::from_fn(genes, H, |_, _| normal(&mut rng));
    let z = DMatrix::<f64>::from_fn(n, H, |_, _| normal(&mut rng));
    // Each type's markers read most along an axis the clusters may split on.
    let mut pairs: Vec<(Box<str>, Box<str>)> = Vec::new();
    for (k, t) in TYPES.iter().enumerate() {
        let mut by: Vec<usize> = (0..genes).collect();
        by.sort_by(|&a, &b| rho[(b, k)].total_cmp(&rho[(a, k)]));
        pairs.extend(
            by[..MARKERS]
                .iter()
                .map(|g| (format!("g{g}").into(), (*t).into())),
        );
    }
    let names: Vec<Box<str>> = (0..genes).map(|g| format!("g{g}").into()).collect();
    let mut a = args("");
    a.seed = seed;
    // The double dip's clusters: Leiden on z itself, as a round's are.
    let on_z = || {
        let mut x: Mat = z.map(|v| v as f32);
        prepare_cells(CellSpace::Embedding, &mut x);
        leiden_on(&x, CellSpace::Embedding, &a.leiden()).unwrap()
    };
    let (test_on, groups, k) = if fissioned {
        let round: Vec<_> = on_z()
            .0
            .iter()
            .map(|&c| ClusterId::try_from(c).ok())
            .collect();
        let s = split(&z, &round, CellSpace::Embedding, &a).unwrap();
        assert!(
            s.dependence < 2.0,
            "seed {seed}: dependence {}",
            s.dependence
        );
        (s.g, s.groups, s.k)
    } else {
        let (groups, k) = on_z();
        (z.clone(), groups, k)
    };
    let r = Ranking::of(&rho, &test_on, &groups, k).unwrap();
    let s = test_panel(&r, &pairs, &names, &a).unwrap();
    let (mut stand, mut up, mut all) = (0usize, 0usize, 0usize);
    for kk in 0..k {
        for t in s.tests(kk).values() {
            all += 1;
            stand += usize::from(t.stands());
            up += usize::from(t.p_up.is_some_and(|p| p < ALPHA));
        }
    }
    (stand as f64 / all as f64, up as f64 / all as f64)
}

#[test]
fn fission_keeps_false_calls_at_their_level_where_the_double_dip_does_not() {
    let runs = |fissioned| -> (f64, f64) {
        let r: Vec<(f64, f64)> = (1..=6).map(|s| false_calls(fissioned, s)).collect();
        let n = r.len() as f64;
        (
            r.iter().map(|x| x.0).sum::<f64>() / n,
            r.iter().map(|x| x.1).sum::<f64>() / n,
        )
    };
    let (dip_stand, dip_up) = runs(false);
    let (fis_stand, fis_up) = runs(true);
    eprintln!("double dip: stands {dip_stand:.3}, up {dip_up:.3}; fission: stands {fis_stand:.3}, up {fis_up:.3}");
    // At 0.05: the double dip finds markers up in near half of the pairs;
    // fission near the level, a little above it from Σ's fit at this size.
    assert!(
        dip_up > 0.3,
        "the double dip should find markers up: {dip_up}"
    );
    assert!(
        fis_up < 0.09,
        "fission's markers-up test should hold its level: {fis_up}"
    );
    assert!(fis_stand <= fis_up && fis_stand < 0.1, "{fis_stand}");
}

#[test]
fn a_half_with_too_few_markers_is_untested_not_failed() {
    let dir = tempfile::tempdir().unwrap();
    // A's five markers halve into three and two: the two cannot be tested.
    let from = round_with(dir.path(), [5, MARKERS, MARKERS, MARKERS]);
    let c = confirm(&args(&from)).unwrap();
    let a = (c.clusters.iter())
        .find(|x| x.call.as_deref() == Some("A"))
        .expect("a cluster called A");
    // The half of three calls A on f, and the half of two cannot test it on
    // g; the half of two makes no call. Neither is a failed test.
    let [x, y] = &a.held_out;
    let (calls, none) = if x.call.is_some() { (x, y) } else { (y, x) };
    assert_eq!(calls.call.as_deref(), Some("A"), "{:?}", a.held_out);
    assert_eq!((calls.p, calls.p_up), (None, None), "{calls:?}");
    assert_eq!(none.call, None, "{none:?}");
    assert!(!calls.holds && !none.holds);
    assert!(render(&c).contains("A untested"), "{}", render(&c));
}

#[test]
fn a_confirmation_of_other_clusters_is_not_read_as_the_rounds() {
    let dir = tempfile::tempdir().unwrap();
    let from = round(dir.path());
    confirm(&args(&from)).unwrap();
    let loaded = run::load(&from).unwrap();
    assert!(read(&loaded).unwrap().is_some());
    // The round's clusters change under it.
    let path = run::resolve(
        &loaded.dir,
        loaded.manifest.cluster.clusters.as_deref().unwrap(),
    );
    let (cells, mut ids) = crate::manifest::rounds::read_clusters(&path).unwrap();
    ids.reverse();
    crate::manifest::rounds::write_clusters(&path, &cells, &ids).unwrap();
    assert!(read(&loaded).unwrap().is_none());
}
