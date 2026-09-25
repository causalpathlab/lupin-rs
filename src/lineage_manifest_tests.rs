//! What a run manifest tells the lineage θ resolver.
//!
//! The resolver itself is tested in `lineage::input`; here it is the
//! mapping from a stamped `kind` to the contract.

use super::*;
use crate::run_manifest::{default_path, RunKind, RunManifest};
use std::path::Path;

/// Minimal `{prefix}.senna.json` stating only the kind.
fn write_kind_only(prefix: &str, kind: RunKind) -> anyhow::Result<()> {
    RunManifest::new(kind, prefix).save(Path::new(&default_path(prefix)))
}

/// A unique scratch prefix per test, so the manifests written here cannot collide.
fn scratch(tag: &str) -> String {
    std::env::temp_dir()
        .join(format!(
            "lupin_lineage_manifest_{}_{tag}",
            std::process::id()
        ))
        .to_string_lossy()
        .into_owned()
}

#[test]
fn a_topic_run_promises_a_log_simplex_latent() {
    let p = scratch("topic");
    write_kind_only(&p, RunKind::Topic).unwrap();
    let c = load_contract(&p).0;
    assert!(c.latent_is_log_simplex);
    assert!(!c.is_gem());
    assert_eq!(c.kind.as_deref(), Some("topic"));
}

#[test]
fn a_gem_run_promises_a_euclidean_co_embedding_and_no_velocity() {
    let p = scratch("gem");
    write_kind_only(&p, RunKind::Gem).unwrap();
    let c = load_contract(&p).0;
    assert!(
        !c.latent_is_log_simplex,
        "gem writes no latent.parquet at all"
    );
    assert!(c.is_gem(), "the caller uses this to explain the missing δ");
}

#[test]
fn an_unreadable_manifest_yields_an_unknown_contract() {
    // Nothing stamped at this prefix: the resolver must be told that, not
    // handed a default that reads as a real claim.
    let c = load_contract(&scratch("nothing_here")).0;
    assert!(c.kind.is_none());
    assert!(!c.latent_is_log_simplex);
    assert!(
        c.source.contains("senna.json"),
        "the fallback warning names where a manifest was looked for, got {}",
        c.source
    );
}
