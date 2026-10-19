//! `annotate --method projection --projection units`, end to end on the
//! planted toy run ([`crate::annotate::marker_score::toy`]), its gene table
//! read as the genes placed among the cells.

use super::*;
use crate::annotate::geneset_units::{TYPE_LINKS_TSV, UNITS_PARQUET};
use crate::annotate::marker_score::toy;
use crate::manifest::rounds::read_argmax;
use crate::manifest::run;

#[test]
fn the_units_pass_calls_each_cell_and_records_its_round() {
    let dir = tempfile::tempdir().unwrap();
    let (manifest, panel) = toy::write_run(dir.path());
    // The toy's genes double as the ones placed on the cell manifold.
    let mut m: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&manifest).unwrap()).unwrap();
    m["outputs"]["feature_coembedding"] = "run.feature_embedding.parquet".into();
    std::fs::write(&manifest, m.to_string()).unwrap();
    let loaded = run::load(&manifest).unwrap();
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    let args = UnitsArgs {
        epochs: 60,
        seed: 7,
        ..UnitsArgs::new(&panel, &out)
    };
    annotate_by_units(&args, &loaded).unwrap();

    let argmax = read_argmax(&format!("{out}.argmax.tsv")).unwrap();
    for i in 0..20 {
        let (label, p) = &argmax[&format!("c{i:02}")];
        assert_eq!(label, if i < 10 { "A" } else { "B" }, "c{i:02}");
        assert!(*p > 0.5, "a call is the region the cell most belongs to");
    }
    let round = run::load(&annotated_path(&loaded.file, &out).to_string_lossy()).unwrap();
    let a = &round.manifest.annotate;
    assert!(a
        .celltype_units
        .as_deref()
        .is_some_and(|p| p.ends_with(UNITS_PARQUET)));
    assert!(a
        .celltype_links
        .as_deref()
        .is_some_and(|p| p.ends_with(TYPE_LINKS_TSV)));
    assert!(
        round.manifest.cluster.clusters.is_some(),
        "the round's cluster table"
    );
    assert_eq!(
        a.settings.as_ref().unwrap()["projection"]["projection"],
        "units"
    );
}
