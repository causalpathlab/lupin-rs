//! [`super`] on the planted toy run ([`super::toy`]).

use super::toy;
use super::*;
use enrichment::consensus::UNASSIGNED;

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
fn a_clusters_share_counts_every_cell_it_holds_in_rows_k0_up() {
    let calls = [Some(0), Some(0), Some(1), None, Some(1)];
    let groups = [3, 3, 3, 3, UNASSIGNED];
    let share = share_by_cluster(&calls, &groups, 4, 2);
    assert_eq!(share.nrows(), 4, "K0..K3, as every cluster table");
    assert_eq!(share[(3, 0)], 0.5);
    assert_eq!(
        share[(3, 1)],
        0.25,
        "an uncalled cell counts in its cluster"
    );
    assert_eq!(share[(0, 0)], 0.0, "an empty cluster has no share");
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
