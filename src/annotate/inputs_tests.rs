//! Cluster-parquet alignment for [`super::load_cluster_labels`].

use super::load_cluster_labels;
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::traits::IoOps;

fn write_clusters(dir: &std::path::Path, name: &str, cells: &[&str], labels: &[f32]) -> String {
    let mut m = Mat::zeros(cells.len(), 1);
    for (i, &lab) in labels.iter().enumerate() {
        m[(i, 0)] = lab;
    }
    let rows: Vec<Box<str>> = cells.iter().map(|s| Box::from(*s)).collect();
    let cols: Vec<Box<str>> = vec!["cluster".into()];
    let path = dir.join(name).to_string_lossy().into_owned();
    m.to_parquet_with_names(&path, (Some(&rows), Some("cell")), Some(&cols))
        .unwrap();
    path
}

#[test]
fn aligns_cluster_parquet_to_data_cell_order() {
    let dir = tempfile::tempdir().unwrap();
    let path = write_clusters(dir.path(), "c.parquet", &["b", "a", "c"], &[1.0, 0.0, 2.0]);
    let cell_names: Vec<Box<str>> = ["a", "b", "c"].iter().map(|s| Box::from(*s)).collect();
    let (labels, n_clusters) = load_cluster_labels(&path, &cell_names).unwrap();
    assert_eq!(labels, vec![0, 1, 2]);
    assert_eq!(n_clusters, 3);
}

#[test]
fn rejects_parquet_with_no_overlap_to_data_cells() {
    let dir = tempfile::tempdir().unwrap();
    let path = write_clusters(dir.path(), "c.parquet", &["x", "y"], &[0.0, 1.0]);
    let cell_names: Vec<Box<str>> = vec!["a".into(), "b".into()];
    assert!(load_cluster_labels(&path, &cell_names).is_err());
}

#[test]
fn a_membership_row_with_undefined_entropy_is_unassigned() {
    // Propensity-style table: C0, C1, cluster, entropy. Cell "z" has no
    // membership, so its cluster 0 is an argmax of zeros.
    let dir = tempfile::tempdir().unwrap();
    let rows = [
        ("a", [0.9, 0.1, 0.0, 0.3]),
        ("b", [0.2, 0.8, 1.0, 0.5]),
        ("z", [0.0, 0.0, 0.0, f32::NAN]),
    ];
    let mut m = Mat::zeros(rows.len(), 4);
    for (i, (_, v)) in rows.iter().enumerate() {
        for (j, &x) in v.iter().enumerate() {
            m[(i, j)] = x;
        }
    }
    let names: Vec<Box<str>> = rows.iter().map(|(n, _)| Box::from(*n)).collect();
    let cols: Vec<Box<str>> = ["C0", "C1", "cluster", "entropy"]
        .iter()
        .map(|s| Box::from(*s))
        .collect();
    let path = dir.path().join("p.parquet").to_string_lossy().into_owned();
    m.to_parquet_with_names(&path, (Some(&names), Some("cell")), Some(&cols))
        .unwrap();

    let (labels, n_clusters) = load_cluster_labels(&path, &names).unwrap();
    assert_eq!(labels, vec![0, 1, usize::MAX]);
    assert_eq!(n_clusters, 2);
}
