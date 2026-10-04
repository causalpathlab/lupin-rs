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
        ..TrajectoryData::default()
    }
}

/// `t`'s figure of `panel` at 320 × 240 in `style`, whole.
fn draw(t: &TrajectoryData, panel: Panel, style: &Style) -> Figure {
    t.figure(panel, 320, 240, style, View::default()).unwrap()
}

#[test]
fn every_panel_renders_to_svg_and_pixels() {
    let t = data();
    let panels = t.panels();
    assert_eq!(panels.len(), 3, "one scatter, the order, the connectivity");
    for p in panels.into_iter().chain(t.scatters()) {
        let fig = draw(&t, p, &Style::default());
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
    assert_eq!(
        t.title(Panel::Layout { k: 0 }, Colouring::Pseudotime),
        "PHATE · pseudotime"
    );
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
    let first = t
        .export(Panel::Order, &Style::default(), View::default())
        .unwrap();
    assert_eq!(first.base, dir.join("x.trajectory.order"));
    assert!(first.files.iter().all(|f| f.is_file()), "{:?}", first.files);
    let second = t
        .export(Panel::Order, &Style::default(), View::default())
        .unwrap();
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
    let t = data();
    let mut style = Style::default();
    let labels = |style: &Style| {
        let svg = draw(&t, Panel::Layout { k: 0 }, style).svg;
        let a = svg.matches(">A</text>").count();
        let b = svg.matches(">B</text>").count();
        (a, b, svg)
    };
    let (a, b, svg) = labels(&style);
    assert_eq!((a, b), (1, 1), "one label per type, on at medium");
    let medium = svg
        .find("font-size='")
        .map(|i| svg[i..].to_string())
        .unwrap();
    assert_eq!(style.cycle_labels(), "labels large · t for largest");
    let (_, _, svg) = labels(&style);
    let large = svg
        .find("font-size='")
        .map(|i| svg[i..].to_string())
        .unwrap();
    assert_ne!(medium[..20], large[..20], "the size changed");
    style.cycle_labels();
    assert_eq!(style.cycle_labels(), "labels off · t shows them small");
    assert_eq!(labels(&style).0 + labels(&style).1, 0);
    assert_eq!(style.cycle_labels(), "labels small · t for medium");
}

#[test]
fn the_unassigned_type_gets_no_label() {
    let mut t = data();
    t.types[39] = enrichment::UNASSIGNED_LABEL.into();
    let svg = draw(&t, Panel::Layout { k: 0 }, &Style::default()).svg;
    assert!(!svg.contains(&format!(">{}</text>", enrichment::UNASSIGNED_LABEL)));
}

#[test]
fn colourings_cycle_over_what_the_outputs_support() {
    let t = data();
    let all = t.colourings();
    assert_eq!(
        all,
        vec![Colouring::Pseudotime, Colouring::Type, Colouring::Component],
        "a single lineage is no colouring"
    );
    let mut style = Style::default();
    let by_time = draw(&t, Panel::Layout { k: 0 }, &style).svg;
    assert_eq!(
        style.cycle_colouring(&all),
        "coloured by cell type · c for component"
    );
    let by_type = draw(&t, Panel::Layout { k: 0 }, &style).svg;
    assert_ne!(by_time, by_type);
    assert!(
        by_type.contains(">A</text>") && by_type.contains("<rect"),
        "a legend"
    );
    assert_eq!(
        t.title(Panel::Layout { k: 0 }, style.colouring),
        "layout · cell type"
    );
    assert_eq!(
        style.cycle_colouring(&all),
        "coloured by component · c for pseudotime"
    );
    assert_eq!(
        style.cycle_colouring(&all),
        "coloured by pseudotime · c for cell type"
    );
}

#[test]
fn a_cell_takes_the_one_lineage_it_weighs_most_on() {
    let w = vec![vec![0.2, 0.0, 0.5], vec![0.7, 0.0, 0.5]];
    assert_eq!(sole_lineage(&w, 0), 1);
    assert_eq!(sole_lineage(&w, 1), -1, "on none");
    assert_eq!(sole_lineage(&w, 2), -1, "a tie");
}

#[test]
fn the_view_zooms_pans_and_stays_inside_the_scatter() {
    let mut v = View::default();
    assert!(v.is_whole());
    v.pan(1.0, 0.0);
    assert_eq!(v, View::default(), "the whole scatter does not pan");
    v.zoom(false);
    assert_eq!(v.level, 0, "no further out than the whole");
    v.zoom(true);
    v.zoom(true);
    assert_eq!(v.factor(), 2.0);
    let half = 0.5 / v.factor();
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
    assert_eq!(v.level, MAX_ZOOM_LEVEL);
    assert_eq!(v.factor(), 64.0);
    let whole = DataBounds {
        xmin: 0.0,
        xmax: 10.0,
        ymin: -5.0,
        ymax: 5.0,
    };
    let part = View {
        level: 2,
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
fn zooming_in_and_back_out_shows_the_whole_scatter() {
    let mut v = View::default();
    for _ in 0..3 {
        v.zoom(true);
    }
    v.pan(1.0, 1.0);
    for _ in 0..3 {
        v.zoom(false);
    }
    assert!(v.is_whole());
    assert_eq!(v.factor(), 1.0);
    assert_eq!(v, View::default(), "and centred again");
}

#[test]
fn a_zoomed_scatter_labels_only_the_types_on_screen() {
    let t = data();
    let style = Style::default();
    let svg = |view: View| {
        t.figure(Panel::Layout { k: 0 }, 320, 240, &style, view)
            .unwrap()
            .svg
    };
    let whole = svg(View::default());
    assert!(whole.contains(">A</text>") && whole.contains(">B</text>"));
    // The left half holds only type A's cells.
    let left = svg(View {
        level: 2,
        cx: 0.25,
        cy: 0.5,
    });
    assert!(left.contains(">A</text>"));
    assert!(!left.contains(">B</text>"), "B is off screen");
    assert_ne!(whole, left);
    assert_eq!(
        svg(View::default()),
        whole,
        "the medians kept draw the same"
    );
    // The order panel ignores the view.
    let zoomed = View {
        level: 2,
        ..View::default()
    };
    assert_eq!(
        t.figure(Panel::Order, 320, 240, &style, zoomed)
            .unwrap()
            .svg,
        draw(&t, Panel::Order, &style).svg
    );
}

#[test]
fn every_type_gets_its_own_colour_however_many_there_are() {
    let names: Vec<String> = (1..=25).map(|k| format!("CT{k}")).collect();
    let keys: Vec<Option<&str>> = names.iter().map(|n| Some(n.as_str())).collect();
    let (colours, legend) = categorical(&keys, str::to_string);
    assert_eq!(legend.len(), 25);
    let mut seen: Vec<Rgb> = colours.iter().flatten().copied().collect();
    seen.sort_unstable();
    seen.dedup();
    assert_eq!(seen.len(), 25, "no two types share a colour");
}

#[test]
fn the_legend_fits_the_figure_or_stays_out_of_a_thumbnail() {
    let entries: Vec<(String, Rgb)> = (1..=25).map(|k| (format!("CT{k}"), INK)).collect();
    assert_eq!(legend_svg(&entries, 6.0, 20), "", "no room for two rows");
    let rows = |svg: &str| svg.matches("<rect").count();
    let small = legend_svg(&entries, 6.0, 100);
    assert!(rows(&small) <= 5, "half of 100 px holds few rows");
    assert!(small.contains("more"), "and says how many are left out");
    let large = legend_svg(&entries, 6.0, 2000);
    assert_eq!(rows(&large), LEGEND_MAX - 1);
    assert!(large.contains("+6 more"));
    assert_eq!(rows(&legend_svg(&entries[..3], 6.0, 2000)), 3);
}

#[test]
fn a_label_near_the_edge_stays_inside_the_figure() {
    let ext = Extent { w: 400, h: 300 };
    let name = "CT1_with_a_long_name";
    let half = 0.275 * 12.0 * name.len() as f32;
    let (x, y) = keep_inside((0.0, 0.0), 12.0, name, ext);
    assert!(x - half >= 0.0 && y >= 12.0, "{x} {y}");
    let (x, _) = keep_inside((400.0, 150.0), 12.0, name, ext);
    assert!(x + half <= 400.0, "{x}");
    assert_eq!(keep_inside((200.0, 150.0), 12.0, name, ext), (200.0, 150.0));
}
