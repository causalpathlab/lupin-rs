use super::*;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use rand_distr::{Distribution, Gamma, LogNormal, Poisson};

/// Simulated pseudobulks: `n_free` genes no type claims, then each type's
/// markers (`sets`, gene indices after the free ones); cluster k's markers of
/// the types in `planted[k]` scaled by the given fold. NB with φ = 0.05.
struct Sim {
    sums: Vec<f64>,
    markers: Vec<Vec<usize>>,
    n_genes: usize,
    n_clusters: usize,
}

fn simulate(n_free: usize, sets: &[Vec<usize>], planted: &[Vec<(usize, f64)>], seed: u64) -> Sim {
    let mut rng = SmallRng::seed_from_u64(seed);
    let n_marker = sets.iter().flatten().max().map_or(0, |m| m + 1);
    let n_genes = n_free + n_marker;
    let mut markers = vec![Vec::new(); n_genes];
    for (c, set) in sets.iter().enumerate() {
        for &g in set {
            markers[n_free + g].push(c);
        }
    }
    let rate = LogNormal::new(0.0, 1.0).unwrap();
    let base: Vec<f64> = (0..n_genes).map(|_| rate.sample(&mut rng)).collect();
    let z: f64 = base.iter().sum();
    let phi = 0.05;
    let mut sums = Vec::with_capacity(n_genes * planted.len());
    for fold_of in planted {
        let lib = 2e5;
        for g in 0..n_genes {
            let mut mu = lib * base[g] / z;
            for &(c, fold) in fold_of {
                if markers[g].contains(&c) {
                    mu *= fold;
                }
            }
            let lam = Gamma::new(1.0 / phi, mu * phi).unwrap().sample(&mut rng);
            sums.push(Poisson::new(lam.max(1e-9)).unwrap().sample(&mut rng));
        }
    }
    Sim {
        sums,
        markers,
        n_genes,
        n_clusters: planted.len(),
    }
}

fn fit(sim: &Sim, n_types: usize) -> SusieFits {
    let cfg = SusieConfig {
        samples: 400,
        warmup: 300,
        ..Default::default()
    };
    fit_all(
        &sim.sums,
        sim.n_genes,
        sim.n_clusters,
        &sim.markers,
        n_types,
        &cfg,
    )
    .unwrap()
}

/// Four disjoint 30-gene marker sets.
fn disjoint() -> Vec<Vec<usize>> {
    (0..4).map(|c| (c * 30..(c + 1) * 30).collect()).collect()
}

#[test]
fn a_planted_type_is_found_in_each_cluster() {
    let planted: Vec<Vec<(usize, f64)>> = (0..8).map(|k| vec![(k % 4, 4.0)]).collect();
    let sim = simulate(400, &disjoint(), &planted, 1);
    let fits = fit(&sim, 4);
    for (k, f) in fits.clusters.iter().enumerate() {
        let top = (0..4)
            .max_by(|&a, &b| f.pip[a].total_cmp(&f.pip[b]))
            .unwrap();
        assert_eq!(top, k % 4, "cluster {k}: pip {:?}", f.pip);
        assert!(f.pip[top] > 0.9, "cluster {k}: pip {:?}", f.pip);
        assert!(f.credible_sets.contains(&vec![top]), "cluster {k}");
    }
}

#[test]
fn shared_markers_are_credited_to_the_type_that_has_them_all() {
    // A and B share 18 of their 50 markers; only A is up, in cluster 0.
    let a: Vec<usize> = (0..50).collect();
    let b: Vec<usize> = (32..82).collect();
    let others: Vec<Vec<usize>> = (0..3)
        .map(|c| (82 + c * 30..82 + (c + 1) * 30).collect())
        .collect();
    let mut sets = vec![a, b];
    sets.extend(others);
    let mut planted: Vec<Vec<(usize, f64)>> = vec![vec![(0, 4.0)]];
    planted.extend((0..6).map(|k| vec![(2 + k % 3, 4.0)]));
    let sim = simulate(400, &sets, &planted, 2);
    let f = &fit(&sim, 5).clusters[0];
    assert!(f.pip[0] > 0.9, "A: pip {:?}", f.pip);
    assert!(f.pip[1] < 0.2, "B: pip {:?}", f.pip);
}

#[test]
fn markers_planted_down_get_no_effect() {
    // Cluster 0: type 0's markers at a fifth of their background.
    let mut planted: Vec<Vec<(usize, f64)>> = vec![vec![(0, 0.2)]];
    planted.extend((0..6).map(|k| vec![(1 + k % 3, 4.0)]));
    let sim = simulate(400, &disjoint(), &planted, 3);
    let f = &fit(&sim, 4).clusters[0];
    assert!(f.theta[0] < 0.1, "theta {:?}", f.theta);
    assert!(
        f.credible_sets.iter().all(|s| !s.contains(&0)),
        "{:?}",
        f.credible_sets
    );
}

#[test]
fn a_cluster_with_no_type_has_no_credible_set() {
    let mut planted: Vec<Vec<(usize, f64)>> = vec![vec![]];
    planted.extend((0..8).map(|k| vec![(k % 4, 4.0)]));
    let sim = simulate(400, &disjoint(), &planted, 4);
    let f = &fit(&sim, 4).clusters[0];
    assert!(
        f.credible_sets.is_empty(),
        "{:?} pip {:?}",
        f.credible_sets,
        f.pip
    );
}

#[test]
fn the_poisson_limit_finds_the_planted_types_too() {
    let planted: Vec<Vec<(usize, f64)>> = (0..8).map(|k| vec![(k % 4, 4.0)]).collect();
    let sim = simulate(400, &disjoint(), &planted, 5);
    let cfg = SusieConfig {
        samples: 400,
        warmup: 300,
        dispersion: Some(0.0),
        ..Default::default()
    };
    let fits = fit_all(
        &sim.sums,
        sim.n_genes,
        sim.n_clusters,
        &sim.markers,
        4,
        &cfg,
    )
    .unwrap();
    assert!(fits.dispersion.iter().all(|&d| d == 0.0));
    for (k, f) in fits.clusters.iter().enumerate() {
        assert!(f.pip[k % 4] > 0.9, "cluster {k}: pip {:?}", f.pip);
    }
}

#[test]
fn the_planted_type_explains_most_of_its_cluster() {
    // Cluster 0: type 0 up strongly, type 1 up a little.
    let mut planted: Vec<Vec<(usize, f64)>> = vec![vec![(0, 6.0), (1, 1.5)]];
    planted.extend((0..8).map(|k| vec![(k % 4, 4.0)]));
    let sim = simulate(400, &disjoint(), &planted, 8);
    let f = &fit(&sim, 4).clusters[0];
    assert!(f.explained[0] > 0.5, "explained {:?}", f.explained);
    assert!(
        f.explained[0] > 2.0 * f.explained[1],
        "explained {:?}",
        f.explained
    );
}

/// A fit with the given PIPs, effects and shares, every type in its own
/// credible set except those with PIP 0.
fn fit_of(pip: &[f32], theta: &[f32], explained: &[f32]) -> ClusterFit {
    ClusterFit {
        pip: pip.to_vec(),
        theta: theta.to_vec(),
        explained: explained.to_vec(),
        credible_sets: (0..pip.len())
            .filter(|&c| pip[c] > 0.0)
            .map(|c| vec![c])
            .collect(),
        max_rhat: 1.0,
    }
}

#[test]
fn a_cluster_is_called_by_pip_then_by_the_share_it_explains() {
    // PIPs tie at 1: the type explaining more of the cluster wins, not the
    // one whose markers rose most.
    let f = fit_of(&[1.0, 1.0, 0.0], &[2.0, 2.6, 0.0], &[0.163, 0.157, 0.0]);
    assert_eq!(f.call(), Some(0));
    // A clearly higher PIP still wins.
    let f = fit_of(&[0.6, 1.0, 0.0], &[2.0, 0.5, 0.0], &[0.5, 0.1, 0.0]);
    assert_eq!(f.call(), Some(1));
    // No credible set, no call.
    let f = fit_of(&[0.0, 0.0, 0.0], &[0.0; 3], &[0.0; 3]);
    assert_eq!(f.call(), None);
}

#[test]
fn empty_clusters_neither_shift_the_background_nor_get_a_call() {
    // Clusters 0..9 are real (cluster 0 no type); 10..21 have no cells, as
    // gaps in cluster ids leave them.
    let mut planted: Vec<Vec<(usize, f64)>> = vec![vec![]];
    planted.extend((0..9).map(|k| vec![(k % 4, 4.0)]));
    let mut sim = simulate(400, &disjoint(), &planted, 9);
    sim.sums.extend(std::iter::repeat_n(0.0, sim.n_genes * 12));
    sim.n_clusters += 12;
    let fits = fit(&sim, 4);
    assert!(
        fits.clusters[0].credible_sets.is_empty(),
        "the null cluster: {:?}",
        fits.clusters[0].pip
    );
    for f in &fits.clusters[10..] {
        assert!(f.call().is_none() && f.pip.iter().all(|&p| p == 0.0));
    }
    for k in 1..10 {
        assert_eq!(fits.clusters[k].call(), Some((k - 1) % 4), "cluster {k}");
    }
}

#[test]
fn settings_that_leave_no_posterior_are_refused() {
    let sim = simulate(100, &disjoint(), &[vec![(0, 4.0)], vec![(1, 4.0)]], 10);
    let run = |cfg: SusieConfig| {
        fit_all(
            &sim.sums,
            sim.n_genes,
            sim.n_clusters,
            &sim.markers,
            4,
            &cfg,
        )
    };
    let mut no_effects = SusieConfig::default();
    no_effects.prior.num_effects = 0;
    assert!(run(no_effects).is_err());
    assert!(run(SusieConfig {
        samples: 0,
        ..Default::default()
    })
    .is_err());
    assert!(run(SusieConfig {
        dispersion: Some(-0.1),
        ..Default::default()
    })
    .is_err());
}

#[test]
fn huge_pseudobulks_still_separate_types_that_share_markers() {
    // As shared_markers_are_credited_…, at ~10⁶ counts per marker gene,
    // where f32 terms lose about a nat each.
    let a: Vec<usize> = (0..50).collect();
    let b: Vec<usize> = (32..82).collect();
    let others: Vec<Vec<usize>> = (0..3)
        .map(|c| (82 + c * 30..82 + (c + 1) * 30).collect())
        .collect();
    let mut sets = vec![a, b];
    sets.extend(others);
    let mut planted: Vec<Vec<(usize, f64)>> = vec![vec![(0, 4.0)]];
    planted.extend((0..6).map(|k| vec![(2 + k % 3, 4.0)]));
    let mut sim = simulate(400, &sets, &planted, 11);
    sim.sums.iter_mut().for_each(|v| *v *= 2_000.0);
    let f = &fit(&sim, 5).clusters[0];
    assert!(f.pip[0] > 0.9, "A: pip {:?}", f.pip);
    assert!(f.pip[1] < 0.2, "B: pip {:?}", f.pip);
}

#[test]
fn chains_that_settle_apart_are_flagged() {
    let same = vec![vec![1.0f32, 1.1, 0.9, 1.0, 1.05, 0.95]; 4];
    assert!(between_chain_rhat(&same) < 1.05);
    let apart = vec![
        vec![0.0f32, 0.1, -0.1, 0.0, 0.05, -0.05],
        vec![2.0f32, 2.1, 1.9, 2.0, 2.05, 1.95],
    ];
    assert!(between_chain_rhat(&apart) > 2.0);
}

#[test]
fn chains_that_number_the_same_effects_differently_agree_on_the_credible_sets() {
    use model::{Cluster, ClusterCounts, Panel, SusieSample};
    let panel = Panel::new(vec![vec![0], vec![1], vec![]], 2, &[0.1, 0.1, 0.1]);
    let cluster = Cluster {
        y: vec![10.0, 10.0, 10.0],
        log_mu0: vec![0.0, 0.0, 0.0],
    };
    let counts = ClusterCounts {
        panel: &panel,
        cluster: &cluster,
    };
    let draw = |first: usize| SusieSample {
        probs: vec![
            (0..3).map(|j| f32::from(u8::from(j == first))).collect(),
            (0..3)
                .map(|j| f32::from(u8::from(j == 1 - first)))
                .collect(),
        ],
        theta: vec![1.0, 1.0],
    };
    // Chain 0 numbers the types' effects (0, 1); chain 1, (1, 0).
    let chains = vec![vec![draw(0); 10], vec![draw(1); 10]];
    let corr = marker_correlation(&[vec![0], vec![1], vec![]], 2);
    let mut sets = summarize(&chains, counts, &corr).credible_sets;
    sets.sort();
    assert_eq!(sets, vec![vec![0], vec![1]]);
}
