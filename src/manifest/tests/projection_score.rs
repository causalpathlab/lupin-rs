//! `annotate --method projection --projection score`, end to end on the
//! planted toy run ([`crate::annotate::marker_score::toy`]).

use super::*;
use crate::annotate::marker_score::toy;
use crate::manifest::rounds::read_argmax;
use crate::manifest::run;
use legume_numeric::matrix::traits::IoOps;

#[test]
fn the_score_pass_calls_each_cell_and_records_its_round() {
    let dir = tempfile::tempdir().unwrap();
    let (manifest, panel) = toy::write_run(dir.path());
    let loaded = run::load(&manifest).unwrap();
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    annotate_by_marker_score(&toy::args(&panel, &out), &loaded).unwrap();

    let argmax = read_argmax(&format!("{out}.argmax.tsv")).unwrap();
    for i in 0..20 {
        let (label, p) = &argmax[&format!("c{i:02}")];
        assert_eq!(label, if i < 10 { "A" } else { "B" }, "c{i:02}");
        assert!(p.is_nan(), "the score has no probability");
    }
    let (_, nums) = read_table_columns(&format!("{out}{ANNOT_PARQUET}"), &[], &["margin"]).unwrap();
    assert!(
        nums[0].iter().all(|&m| m > 0.0),
        "every planted call has a margin"
    );

    let round = run::load(&annotated_path(&loaded.file, &out).to_string_lossy()).unwrap();
    let a = &round.manifest.annotate;
    assert!(a.argmax.is_some() && a.annotation.is_some());
    let share = a.cluster_celltype_share.as_deref().expect("a share table");
    let share = Mat::from_parquet(&resolve(&round.dir, share)).unwrap();
    assert_eq!(share.rows, ["K0".into(), "K1".into()] as [Box<str>; 2]);
    assert_eq!((share.mat[(0, 0)], share.mat[(1, 1)]), (1.0, 1.0));
    assert!(
        round.manifest.cluster.clusters.is_some(),
        "the round's cluster table"
    );
    let recorded = &a.settings.as_ref().unwrap()["projection"];
    assert_eq!(recorded["projection"], "score");
    assert!(
        recorded.get("num_perm").is_none(),
        "no ORA settings on a score round"
    );
}
