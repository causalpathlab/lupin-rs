//! [`super`] on a planted toy run.

use super::*;

/// Six genes: g0..g2 load on axis 0, g3..g5 on axis 1. Ten cells sit out on
/// each axis; clusters 0 and 1 hold them, and a seventh cell is unassigned.
fn toy(dir: &std::path::Path) -> (MatWithNames<Mat>, Mat, Vec<usize>, String) {
    let genes: Vec<Box<str>> = (0..6).map(|g| format!("g{g}").into()).collect();
    let mut rho = Mat::zeros(6, 2);
    for g in 0..6 {
        rho[(g, g / 3)] = 1.0 + 0.1 * g as f32;
        rho[(g, 1 - g / 3)] = 0.05 * g as f32;
    }
    let mut cells = Mat::zeros(21, 2);
    for i in 0..20 {
        let axis = i / 10;
        cells[(i, axis)] = 2.0 + 0.1 * (i % 10) as f32;
        cells[(i, 1 - axis)] = 0.2 * ((i % 3) as f32 - 1.0);
    }
    let mut groups: Vec<usize> = (0..20).map(|i| i / 10).collect();
    groups.push(UNASSIGNED);
    let panel = dir.join("panel.tsv");
    std::fs::write(&panel, "g0\tA\ng1\tA\ng2\tA\ng3\tB\ng4\tB\ng5\tB\n").unwrap();
    let feat = MatWithNames {
        mat: rho,
        rows: genes,
        cols: vec!["h0".into(), "h1".into()],
    };
    (feat, cells, groups, panel.to_string_lossy().into_owned())
}

#[test]
fn planted_cells_get_their_type_with_a_positive_margin() {
    let dir = tempfile::tempdir().unwrap();
    let (feat, cells, groups, panel) = toy(dir.path());
    let s = score_cells(&feat, &cells, &groups, &panel).unwrap();
    assert_eq!(s.types, ["A".into(), "B".into()] as [Box<str>; 2]);
    for i in 0..20 {
        assert_eq!(s.calls[i], Some(i / 10), "cell {i}");
        assert!(s.margins[i] > 0.0, "cell {i}: {}", s.margins[i]);
    }
}

#[test]
fn a_clusters_share_counts_every_cell_it_holds() {
    let calls = [Some(0), Some(0), Some(1), None, Some(1)];
    let groups = [3, 3, 3, 3, UNASSIGNED];
    let (ids, share) = share_by_cluster(&calls, &groups, 2);
    assert_eq!(ids, [3]);
    assert_eq!(share[(0, 0)], 0.5);
    assert_eq!(
        share[(0, 1)],
        0.25,
        "an uncalled cell counts in its cluster"
    );
}

#[test]
fn a_one_type_panel_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let (feat, cells, groups, _) = toy(dir.path());
    let panel = dir.path().join("one.tsv");
    std::fs::write(&panel, "g0\tA\ng1\tA\ng2\tA\n").unwrap();
    let err = score_cells(&feat, &cells, &groups, &panel.to_string_lossy())
        .err()
        .expect("an error")
        .to_string();
    assert!(err.contains('A'), "{err}");
}
