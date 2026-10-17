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
