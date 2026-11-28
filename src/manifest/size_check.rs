//! What the size null did to a round, as review shows it: the clusters whose
//! enrichment top it moved, and how far the panel's wins follow its sizes
//! before and after.

use super::rounds::SIZE_NULL_FIT;
use super::run::Loaded;
use crate::annotate::markers::label_key;
use crate::annotate::rounds::{ClusterId, Method, TypeEvidence};
use crate::annotate::size_null::SizeTrend;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// The enrichment's size-null fit the round records.
pub(crate) fn recorded_fit(loaded: &Loaded) -> Option<SizeTrend> {
    let s = loaded.manifest.annotate.settings.as_ref()?;
    SizeTrend::deserialize(s.pointer(&format!("/enrichment/{SIZE_NULL_FIT}"))?).ok()
}

/// The tested types (raw p below 1) by the random-set test's p, before the
/// size null, smallest first; ties go to the larger NES.
pub(crate) fn raw_order(types: &[TypeEvidence]) -> Vec<&TypeEvidence> {
    let high = |v: Option<f32>| v.unwrap_or(f32::NEG_INFINITY);
    let mut tested: Vec<&TypeEvidence> = (types.iter())
        .filter(|e| e.p.is_some_and(|p| p < 1.0))
        .collect();
    tested.sort_by(|a, b| {
        (a.p.unwrap().total_cmp(&b.p.unwrap())).then(high(b.nes).total_cmp(&high(a.nes)))
    });
    tested
}

/// The first of [`raw_order`].
pub(crate) fn raw_top(types: &[TypeEvidence]) -> Option<&TypeEvidence> {
    raw_order(types).first().copied()
}

/// A cluster whose enrichment top the size null moved.
#[derive(Serialize, Debug, Clone, PartialEq)]
pub(crate) struct Moved {
    /// Its top by the raw p, and by the size null.
    pub raw_top: String,
    pub top: String,
}

/// The enrichment's called top type, by the size null.
fn called_top(types: &[TypeEvidence]) -> Option<String> {
    (Method::Enrichment.ranked(types).first()).map(|e| e.label.clone())
}

/// Where `fit` moved a cluster's enrichment top: `None` when it did not, no
/// fit moved anything, or the cluster has no call (the round keeps no raw
/// q-values, so whether the raw top would have been called is not known).
pub(crate) fn moved(fit: &SizeTrend, types: &[TypeEvidence]) -> Option<Moved> {
    if !fit.moves() {
        return None;
    }
    between(called_top(types)?, raw_top(types)?)
}

/// A move from the raw top `raw` to the called `top`, if they differ.
fn between(top: String, raw: &TypeEvidence) -> Option<Moved> {
    (label_key(&top) != label_key(&raw.label)).then(|| Moved {
        raw_top: raw.label.clone(),
        top,
    })
}

/// How far the panel's wins follow its sizes: the Spearman ρ between a
/// type's log effective size and the called clusters it tops, by the raw p
/// and by the size null.
#[derive(Serialize, Debug, Clone, PartialEq)]
pub(crate) struct SizeFollow {
    pub status: String,
    pub raw: Option<f64>,
    pub calibrated: Option<f64>,
    /// The clusters whose top the size null moved.
    pub moved: BTreeMap<ClusterId, Moved>,
}

/// [`SizeFollow`] over `clusters`' evidence.
pub(crate) fn follow(
    fit: &SizeTrend,
    clusters: &BTreeMap<ClusterId, Vec<TypeEvidence>>,
) -> SizeFollow {
    let slot: BTreeMap<String, usize> = (fit.types.iter().enumerate())
        .map(|(i, t)| (label_key(&t.name), i))
        .collect();
    let n = fit.types.len();
    let (mut raw, mut cal) = (vec![0.0; n], vec![0.0; n]);
    let count = |wins: &mut [f64], label: Option<&str>| {
        if let Some(&i) = label.and_then(|l| slot.get(&label_key(l))) {
            wins[i] += 1.0;
        }
    };
    let mut moved_at = BTreeMap::new();
    // Over the called clusters alone, as `moved` is.
    for (&id, types) in clusters {
        let (Some(top), Some(raw_top)) = (called_top(types), raw_top(types)) else {
            continue;
        };
        count(&mut raw, Some(&raw_top.label));
        count(&mut cal, Some(&top));
        if let Some(m) = between(top, raw_top).filter(|_| fit.moves()) {
            moved_at.insert(id, m);
        }
    }
    SizeFollow {
        status: fit.status.clone(),
        raw: fit.size_correlation(&raw),
        calibrated: fit.size_correlation(&cal),
        moved: moved_at,
    }
}

impl SizeFollow {
    /// Review's line for it.
    pub(crate) fn line(&self) -> String {
        let rho = |v: Option<f64>| v.map_or_else(|| "-".into(), |v| format!("{v:+.2}"));
        format!(
            "size null: {}; clusters won ~ panel size: raw ρ {}, calibrated ρ {}; {} cluster(s) moved\n",
            self.status,
            rho(self.raw),
            rho(self.calibrated),
            self.moved.len()
        )
    }
}

#[cfg(test)]
#[path = "tests/size_check.rs"]
mod tests;
