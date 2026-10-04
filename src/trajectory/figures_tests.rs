use super::*;

fn data() -> TrajectoryData {
    let n = 40;
    let types: Vec<Box<str>> = (0..n)
        .map(|i| if i < 20 { "A" } else { "B" }.into())
        .collect();
    let mut pseudotime: Vec<f32> = (0..n).map(|i| i as f32 / (n - 1) as f32).collect();
    pseudotime[0] = f32::NAN;
    let x: Vec<f32> = (0..n).map(|i| i as f32).collect();
    let y: Vec<f32> = (0..n).map(|i| (i % 7) as f32).collect();
    let mut d = Mat::zeros(n, 3);
    for i in 0..n {
        d[(i, 1)] = x[i];
        d[(i, 2)] = y[i];
    }
    TrajectoryData {
        manifest: "run.senna.json".into(),
        pseudotime,
        types,
        diffusion: Some(d),
        layout: Some((x, y)),
        layout_method: None,
        edges: vec![EdgeRow {
            a: "A".into(),
            b: "B".into(),
            connectivity: 0.8,
            in_prior: true,
            verdict: Some(Verdict::Supported),
            order_agreement: 0.9,
        }],
    }
}

#[test]
fn every_panel_renders_to_svg_and_pixels() {
    let t = data();
    let panels = t.panels();
    assert_eq!(panels.len(), 4);
    for p in panels {
        let fig = t.figure(p, 320, 240).unwrap();
        assert!(fig.svg.starts_with("<?xml"), "{}", p.slug());
        let img = render(&fig).unwrap();
        assert_eq!((img.width(), img.height()), (fig.w, fig.h));
        if p != Panel::Connectivity {
            assert_eq!((fig.w, fig.h), (320, 240), "{}", p.slug());
        }
    }
}

#[test]
fn the_order_panel_lists_types_by_median_and_skips_unreached_cells() {
    let rows = data().type_order();
    assert_eq!(rows[0].name.as_ref(), "A");
    assert_eq!(rows[0].cells, 19, "the NaN cell is left out");
    assert!(rows[0].median < rows[1].median);
}

#[test]
fn phate_is_drawn_when_senna_made_one() {
    let mut m: RunManifest = serde_json::from_str(
        r#"{"version": 2, "kind": "topic", "prefix": "r", "layout": {"cell_coords": "r.umap.parquet", "current": "umap",
            "methods": {"umap": {"cell_coords": "r.umap.parquet"}}}}"#,
    )
    .unwrap();
    assert_eq!(
        layouts(&m),
        vec![(Some("umap".into()), "r.umap.parquet".to_string())]
    );
    m.layout.extra["methods"]["phate"] = serde_json::json!({"cell_coords": "r.phate.parquet"});
    assert_eq!(
        layouts(&m),
        vec![
            (Some("phate".into()), "r.phate.parquet".to_string()),
            (Some("umap".into()), "r.umap.parquet".to_string())
        ]
    );
    let mut t = data();
    t.layout_method = Some("phate".into());
    assert_eq!(t.title(Panel::Layout), "PHATE · pseudotime");
}

#[test]
fn diffusion_pairs_cycle_past_the_trivial_component() {
    let t = data();
    assert_eq!(
        t.next_pair(1, 2, true),
        (1, 2),
        "with three components the only pair after the trivial one is (1, 2)"
    );
    assert_eq!(
        t.next_pair(1, 2, false),
        (1, 2),
        "only two non-trivial components"
    );
}

#[test]
fn export_writes_svg_and_pdf_past_existing_files() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let mut t = data();
    t.manifest = dir.join("x.senna.json");
    let first = t.export(Panel::Order).unwrap();
    assert_eq!(first.base, dir.join("x.trajectory.order"));
    assert!(first.files.iter().all(|f| f.is_file()), "{:?}", first.files);
    let second = t.export(Panel::Order).unwrap();
    assert_eq!(second.base, dir.join("x.trajectory.order-2"));
}
