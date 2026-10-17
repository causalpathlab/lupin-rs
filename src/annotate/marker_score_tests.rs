//! [`super`] on the planted toy run ([`super::toy`]).

use super::toy;
use super::*;

#[test]
fn planted_cells_get_their_type_with_a_positive_margin() {
    let dir = tempfile::tempdir().unwrap();
    let (feat, cells, groups) = toy::matrices();
    let s = score_cells(&feat, &cells.mat, &groups, &toy::panel(dir.path())).unwrap();
    assert_eq!(s.types, ["A".into(), "B".into()] as [Box<str>; 2]);
    for i in 0..20 {
        assert_eq!(s.calls[i], Some(i / 10), "cell {i}");
        assert!(s.margins[i] > 0.0, "cell {i}: {}", s.margins[i]);
    }
}

#[test]
fn a_one_type_panel_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let (feat, cells, groups) = toy::matrices();
    let panel = dir.path().join("one.tsv");
    std::fs::write(&panel, "g0\tA\ng1\tA\ng2\tA\n").unwrap();
    let err = score_cells(&feat, &cells.mat, &groups, &panel.to_string_lossy())
        .err()
        .expect("an error")
        .to_string();
    assert!(err.contains('A'), "{err}");
}

#[test]
fn the_score_beside_another_pass_writes_its_table_and_leaves_the_call_alone() {
    let dir = tempfile::tempdir().unwrap();
    let (feat, cells, groups) = toy::matrices();
    let out = dir.path().join("run.L1").to_string_lossy().into_owned();
    std::fs::write(format!("{out}.argmax.tsv"), "the other pass's call").unwrap();
    let table = alongside(&toy::panel(dir.path()), &out, &feat, &cells, &groups).unwrap();
    assert!(table.ends_with(crate::annotate::outputs::MARKER_SCORE_PARQUET));
    let (strings, _) =
        legume_numeric::matrix::parquet::read_table_columns(&table, &["call"], &[]).unwrap();
    assert_eq!(strings[0][0].as_ref(), "A");
    assert_eq!(
        std::fs::read_to_string(format!("{out}.argmax.tsv")).unwrap(),
        "the other pass's call"
    );
}
