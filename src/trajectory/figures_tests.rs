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
        component: (0..n).map(|i| i32::from(i >= 30)).collect(),
        lineage: (0..n).map(|i| if i < 20 { -1 } else { 0 }).collect(),
        diffusion: Some(d),
        layouts: vec![Layout { method: None, x, y }],
        edges: vec![EdgeRow {
            a: "A".into(),
            b: "B".into(),
            connectivity: 0.8,
            in_prior: true,
            verdict: Some(Verdict::Supported),
            order_agreement: 0.9,
        }],
        style: Style::default(),
        view: View::default(),
    }
}

#[test]
fn every_panel_renders_to_svg_and_pixels() {
    let t = data();
    let panels = t.panels();
    assert_eq!(panels.len(), 3, "one scatter, the order, the connectivity");
    for p in panels.into_iter().chain(t.scatters()) {
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
    t.layouts[0].method = Some("phate".into());
    assert_eq!(t.title(Panel::Layout { k: 0 }), "PHATE · pseudotime");
    assert_eq!(t.slug(Panel::Layout { k: 0 }), "layout_phate");
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

#[test]
fn every_recorded_layout_is_offered_after_phate_and_the_current_one() {
    let m: RunManifest = serde_json::from_str(
        r#"{"version": 2, "kind": "topic", "prefix": "r", "layout": {"cell_coords": "r.umap.parquet", "current": "umap",
            "methods": {"umap": {"cell_coords": "r.umap.parquet"}, "tsne": {"cell_coords": "r.tsne.parquet"},
                        "phate": {"cell_coords": "r.phate.parquet"}}}}"#,
    )
    .unwrap();
    let names: Vec<Option<String>> = layouts(&m).into_iter().map(|(m, _)| m).collect();
    assert_eq!(
        names,
        vec![
            Some("phate".into()),
            Some("umap".into()),
            Some("tsne".into())
        ]
    );
}

#[test]
fn m_steps_through_the_layouts_then_the_diffusion_map() {
    let mut t = data();
    let second = Layout {
        method: Some("umap".into()),
        x: t.layouts[0].y.clone(),
        y: t.layouts[0].x.clone(),
    };
    t.layouts.push(second);
    let a = Panel::Layout { k: 0 };
    let b = t.next_scatter(a).unwrap();
    assert_eq!(b, Panel::Layout { k: 1 });
    let c = t.next_scatter(b).unwrap();
    assert!(matches!(c, Panel::Diffusion { .. }));
    assert_eq!(t.next_scatter(c), Some(a));
    assert_eq!(t.next_scatter(Panel::Order), None);
    t.layouts.clear();
    t.diffusion = None;
    assert_eq!(t.next_scatter(a), None, "nothing to switch to");
}

#[test]
fn type_labels_follow_the_text_size_and_turn_off() {
    let mut t = data();
    let labels = |t: &TrajectoryData| {
        let svg = t.figure(Panel::Layout { k: 0 }, 320, 240).unwrap().svg;
        let a = svg.matches(">A</text>").count();
        let b = svg.matches(">B</text>").count();
        (a, b, svg)
    };
    let (a, b, svg) = labels(&t);
    assert_eq!((a, b), (1, 1), "one label per type, on at medium");
    let medium = svg
        .find("font-size='")
        .map(|i| svg[i..].to_string())
        .unwrap();
    assert_eq!(t.cycle_labels(), "labels large · t for largest");
    let (_, _, svg) = labels(&t);
    let large = svg
        .find("font-size='")
        .map(|i| svg[i..].to_string())
        .unwrap();
    assert_ne!(medium[..20], large[..20], "the size changed");
    t.cycle_labels();
    assert_eq!(t.cycle_labels(), "labels off · t shows them small");
    assert_eq!(labels(&t).0 + labels(&t).1, 0);
    assert_eq!(t.cycle_labels(), "labels small · t for medium");
}

#[test]
fn the_unassigned_type_gets_no_label() {
    let mut t = data();
    t.types[39] = enrichment::UNASSIGNED_LABEL.into();
    let svg = t.figure(Panel::Layout { k: 0 }, 320, 240).unwrap().svg;
    assert!(!svg.contains(&format!(">{}</text>", enrichment::UNASSIGNED_LABEL)));
}

#[test]
fn colourings_cycle_over_what_the_outputs_support() {
    let mut t = data();
    assert_eq!(
        t.colourings(),
        vec![Colouring::Pseudotime, Colouring::Type, Colouring::Component],
        "a single lineage is no colouring"
    );
    let by_time = t.figure(Panel::Layout { k: 0 }, 160, 120).unwrap().svg;
    assert_eq!(
        t.cycle_colouring(),
        "coloured by cell type · c for component"
    );
    let by_type = t.figure(Panel::Layout { k: 0 }, 160, 120).unwrap().svg;
    assert_ne!(by_time, by_type);
    assert!(
        by_type.contains(">A</text>") && by_type.contains("<rect"),
        "a legend"
    );
    assert_eq!(t.title(Panel::Layout { k: 0 }), "layout · cell type");
    assert_eq!(
        t.cycle_colouring(),
        "coloured by component · c for pseudotime"
    );
    assert_eq!(
        t.cycle_colouring(),
        "coloured by pseudotime · c for cell type"
    );
}

#[test]
fn the_view_zooms_pans_and_stays_inside_the_scatter() {
    let mut v = View::default();
    assert!(v.is_whole());
    v.pan(1.0, 0.0);
    assert_eq!(v, View::default(), "the whole scatter does not pan");
    v.zoom(false);
    assert_eq!(v.zoom, 1.0, "no further out than the whole");
    v.zoom(true);
    v.zoom(true);
    assert!((v.zoom - ZOOM_STEP * ZOOM_STEP).abs() < 1e-5);
    let half = 0.5 / v.zoom;
    for _ in 0..50 {
        v.pan(1.0, -1.0);
    }
    assert!(
        (v.cx - (1.0 - half)).abs() < 1e-6,
        "stops at the right edge"
    );
    assert!((v.cy - half).abs() < 1e-6, "and at the bottom");
    for _ in 0..50 {
        v.zoom(true);
    }
    assert_eq!(v.zoom, MAX_ZOOM);
    let whole = DataBounds {
        xmin: 0.0,
        xmax: 10.0,
        ymin: -5.0,
        ymax: 5.0,
    };
    let part = View {
        zoom: 2.0,
        cx: 0.25,
        cy: 0.75,
    }
    .of(&whole);
    assert_eq!(
        (part.xmin, part.xmax, part.ymin, part.ymax),
        (0.0, 5.0, 0.0, 5.0)
    );
    let all = View::default().of(&whole);
    assert_eq!((all.xmin, all.xmax), (0.0, 10.0));
}

#[test]
fn a_zoomed_scatter_labels_only_the_types_on_screen() {
    let mut t = data();
    let svg = |t: &TrajectoryData| t.figure(Panel::Layout { k: 0 }, 320, 240).unwrap().svg;
    let whole = svg(&t);
    assert!(whole.contains(">A</text>") && whole.contains(">B</text>"));
    // The left half holds only type A's cells.
    t.view = View {
        zoom: 2.0,
        cx: 0.25,
        cy: 0.5,
    };
    let left = svg(&t);
    assert!(left.contains(">A</text>"));
    assert!(!left.contains(">B</text>"), "B is off screen");
    assert_ne!(whole, left);
    // A thumbnail shows the scatter whole whatever the view.
    let thumb = t
        .figure_at(Panel::Layout { k: 0 }, 320, 240, View::default())
        .unwrap()
        .svg;
    assert_eq!(thumb, whole);
    // The order panel ignores the view.
    assert_eq!(
        t.figure(Panel::Order, 320, 240).unwrap().svg,
        t.figure_at(Panel::Order, 320, 240, View::default())
            .unwrap()
            .svg
    );
}
