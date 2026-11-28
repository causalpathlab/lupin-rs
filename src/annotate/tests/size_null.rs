//! [`super::fit`]: a type's statistic read against panels of its size.

use super::*;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use rand_distr::{Distribution, StandardNormal};

const CLUSTERS: usize = 35;
/// The panel's sizes, as unequal as paper1's (VLMC 18 … Exc deep-layer 381).
const SIZES: [f64; 8] = [18.0, 19.0, 34.0, 40.0, 51.0, 193.0, 256.0, 381.0];

/// Each type's statistic over the clusters: standard normal noise widened by
/// `inflate(size)` (what co-varying markers do to a set's null), lifted by
/// `truth` where `true_in(type, cluster)`.
fn round(
    inflate: impl Fn(f64) -> f64,
    truth: f64,
    true_in: impl Fn(usize, usize) -> bool,
) -> Vec<TypeValues> {
    let mut rng = SmallRng::seed_from_u64(7);
    (SIZES.iter().enumerate())
        .map(|(c, &size)| TypeValues {
            size,
            values: (0..CLUSTERS)
                .map(|k| {
                    let e: f64 = StandardNormal.sample(&mut rng);
                    e * inflate(size) + if true_in(c, k) { truth } else { 0.0 }
                })
                .collect(),
        })
        .collect()
}

/// Each cluster's best type by `s`.
fn calls(types: &[TypeValues], s: impl Fn(usize, f64) -> f64) -> Vec<usize> {
    (0..CLUSTERS)
        .map(|k| {
            (0..types.len())
                .max_by(|&a, &b| s(a, types[a].values[k]).total_cmp(&s(b, types[b].values[k])))
                .unwrap()
        })
        .collect()
}

/// CAMERA's inflation for markers that co-vary at `rho`.
fn vif(rho: f64) -> impl Fn(f64) -> f64 {
    move |m| (1.0 + rho * (m - 1.0)).sqrt()
}

/// Each type is true in clusters 4c..4c+4, so 32 of the 35 have a truth.
fn own(c: usize, k: usize) -> bool {
    k / 4 == c
}

#[test]
fn a_panel_without_a_size_trend_is_read_as_it_is() {
    let types = round(|_| 1.0, 4.0, own);
    let trend = fit(&types);
    for t in &trend.types {
        assert!(
            (t.scale - 1.0).abs() < 0.2,
            "{} markers read at {}",
            t.size,
            t.scale
        );
    }
    let (raw, cal) = (
        calls(&types, |_, s| s),
        calls(&types, |c, s| trend.calibrate(c, s)),
    );
    let same = (0..CLUSTERS).filter(|&k| raw[k] == cal[k]).count();
    assert!(same >= 33, "calls kept on {same} of {CLUSTERS}");
}

#[test]
fn a_large_panel_no_longer_wins_on_its_size() {
    // Markers that co-vary at ρ 0.1: 381 widen the null ×6.2, 18 ×1.6.
    let types = round(vif(0.1), 5.0, own);
    let raw = calls(&types, |_, s| s);
    let largest = SIZES.len() - 1;
    let won_raw = raw.iter().filter(|&&c| c == largest).count();
    assert!(won_raw >= 12, "the raw call follows size: {won_raw}");

    let trend = fit(&types);
    assert_eq!(trend.status, "fitted");
    assert!(
        trend.types[largest].scale > 2.0,
        "{:?}",
        trend.types[largest]
    );
    let cal = calls(&types, |c, s| trend.calibrate(c, s));
    // Read by the true inflation, the best any size null could do.
    let oracle = calls(&types, |c, s| s / vif(0.1)(SIZES[c]));
    let right = |tops: &[usize]| (0..32).filter(|&k| tops[k] == k / 4).count();
    let (raw_right, cal_right, best) = (right(&raw), right(&cal), right(&oracle));
    assert!(
        cal_right > raw_right && cal_right + 2 >= best,
        "calls right on raw {raw_right}, calibrated {cal_right}, oracle {best} of 32"
    );
    // Beyond its own clusters (28..32), the largest type wins by chance alone.
    let stray = |tops: &[usize]| {
        (0..CLUSTERS)
            .filter(|&k| tops[k] == largest && k / 4 != largest)
            .count()
    };
    assert!(
        stray(&cal) * 2 <= stray(&raw),
        "the largest type wins {} clusters not its own, {} raw: {cal:?}",
        stray(&cal),
        stray(&raw)
    );
    // The 18-marker type's own clusters are its again.
    assert!((0..4).all(|k| cal[k] == 0), "{cal:?}");
}

#[test]
fn a_type_in_most_clusters_keeps_them() {
    // The 256-marker type (6) is true in 28 of 35 clusters, as a tumour's
    // neoplastic cells are; the other types sit at chance.
    let types = round(vif(0.02), 10.0, |c, k| c == 6 && k < 28);
    let trend = fit(&types);
    let cal = calls(&types, |c, s| trend.calibrate(c, s));
    let kept = (0..28).filter(|&k| cal[k] == 6).count();
    assert!(kept >= 25, "the dominant type keeps {kept} of 28: {cal:?}");
}

#[test]
fn too_few_types_or_clusters_fit_nothing() {
    let few = |n: usize, k: usize| -> Vec<TypeValues> {
        (0..n)
            .map(|c| TypeValues {
                size: 10.0 * (c + 1) as f64,
                values: (0..k).map(|i| (i * c) as f64).collect(),
            })
            .collect()
    };
    for t in [fit(&few(3, 30)), fit(&few(8, 4))] {
        assert!(t.status.starts_with("too few"), "{}", t.status);
        assert!(!t.moves());
    }
}

#[test]
fn the_effective_size_counts_low_weights_less() {
    assert_eq!(effective_size([1.0; 10]), 10.0);
    let uneven = effective_size([1.0, 1.0, 0.1, 0.1, 0.0]);
    assert!(uneven > 2.0 && uneven < 3.0, "{uneven}");
}

#[test]
fn the_size_correlation_reads_a_monotone_trend_as_one() {
    let trend = SizeTrend::none("off", &SIZES);
    let y: Vec<f64> = SIZES.iter().map(|m| m.sqrt()).collect();
    assert!((trend.size_correlation(&y).unwrap() - 1.0).abs() < 1e-12);
}
