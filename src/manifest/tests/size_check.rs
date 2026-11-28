//! [`super::moved`] and [`super::follow`]: what the size null did to a round.

use super::*;
use crate::annotate::size_null::{fit, TypeValues};

/// A type's evidence: its raw p, and its share under the size null.
fn ev(label: &str, p: f32, share: f32) -> TypeEvidence {
    TypeEvidence {
        label: label.into(),
        p: Some(p),
        q: Some(p),
        enrichment: Some(share),
        ..TypeEvidence::default()
    }
}

/// A fit that moves its types: chance spread rising with size over four
/// types of 10, 20, 40 and 400 markers.
fn moving_fit() -> SizeTrend {
    let sizes = [10.0, 20.0, 40.0, 400.0];
    let types: Vec<TypeValues> = (sizes.iter())
        .map(|&m: &f64| TypeValues {
            size: m,
            values: (0..20).map(|k| f64::from(k) * 0.1 * m.sqrt()).collect(),
        })
        .collect();
    let trend = fit(&types).named(&["A", "B", "C", "Big"]);
    assert!(trend.moves(), "{}", trend.status);
    trend
}

#[test]
fn a_cluster_whose_top_the_size_null_moved_says_so() {
    let fit = moving_fit();
    let types = vec![ev("Big", 1e-6, 0.2), ev("A", 1e-3, 0.8)];
    let m = moved(&fit, &types).expect("moved");
    assert_eq!(m.raw_top, "Big");
    assert_eq!(m.top, "A");
    // Not moved: both orders agree.
    let same = vec![ev("Big", 1e-6, 0.8), ev("A", 1e-3, 0.2)];
    assert_eq!(moved(&fit, &same), None);
}

#[test]
fn a_fit_that_moves_nothing_flags_nothing() {
    let none = SizeTrend::none("no trend", &[10.0, 400.0]).named(&["A", "Big"]);
    let types = vec![ev("Big", 1e-6, 0.2), ev("A", 1e-3, 0.8)];
    assert_eq!(moved(&none, &types), None);
}

#[test]
fn the_wins_that_follow_size_are_counted_before_and_after() {
    let fit = moving_fit();
    // Raw: the 400-marker type tops every cluster; after: each its own.
    let clusters: BTreeMap<ClusterId, Vec<TypeEvidence>> = ["A", "B", "C", "Big"]
        .iter()
        .enumerate()
        .map(|(k, own)| {
            let types = ["A", "B", "C", "Big"]
                .iter()
                .map(|t| {
                    let p = if *t == "Big" { 1e-8 } else { 1e-2 };
                    ev(t, p, if t == own { 1.0 } else { 0.0 })
                })
                .collect();
            (k as ClusterId, types)
        })
        .collect();
    let f = follow(&fit, &clusters);
    assert!(f.raw.unwrap() > 0.7, "{f:?}");
    assert!(
        f.calibrated.is_none(),
        "every type wins once: no spread, {f:?}"
    );
    assert_eq!(f.moved.len(), 3);
    assert!(f.line().contains("3 cluster(s) moved"));
}
