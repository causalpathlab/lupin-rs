use super::*;

fn data() -> CellFigures {
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
    CellFigures {
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
        ..CellFigures::default()
    }
}

/// The size the tests draw at: 320 × 240, text 16 px.
const SIZE: Size = Size {
    w: 320,
    h: 240,
    text_px: 16.0,
};

/// `t`'s figure of `panel` at [`SIZE`] in `style`, whole.
fn draw(t: &CellFigures, panel: Panel, style: &Style) -> Picture {
    t.figure(panel, SIZE, style, View::default(), &Offsets::default())
        .unwrap()
}

/// `t`'s first layout as a cell map at [`SIZE`] in `style` at `view`, its
/// text kept as text.
fn map_at(t: &CellFigures, style: &Style, view: View) -> CellMap {
    t.cell_map(
        Panel::Layout { k: 0 },
        SIZE,
        style,
        view,
        true,
        &Offsets::default(),
    )
    .unwrap()
}

fn map(t: &CellFigures, style: &Style) -> CellMap {
    map_at(t, style, View::default())
}

/// The names of the labels `m` placed on its cells.
fn labels(m: &CellMap) -> Vec<&str> {
    m.placed
        .iter()
        .map(|p| m.names[p.key as usize].as_ref())
        .collect()
}

/// Every text `m` kept: its legend and its labels.
fn texts(m: &CellMap) -> Vec<&str> {
    m.text.iter().map(|r| r.text.as_str()).collect()
}

#[test]
fn every_panel_renders_to_pixels() {
    let t = data();
    let panels = t.panels();
    assert_eq!(panels.len(), 3, "one scatter, the order, the connectivity");
    for p in panels.into_iter().chain(t.scatters()) {
        let fig = draw(&t, p, &Style::default());
        let img = render(&fig).unwrap();
        match &fig {
            Picture::Map(m) => {
                assert!(p.is_scatter(), "{}", p.slug());
                assert_eq!((img.width(), img.height()), (320, 240), "{}", p.slug());
                assert_eq!((m.frame.w, m.frame.h), (320, 240));
            }
            Picture::Svg(f) => {
                assert!(!p.is_scatter(), "{}", p.slug());
                assert!(f.svg.starts_with("<?xml"), "{}", p.slug());
                assert_eq!((img.width(), img.height()), (f.w, f.h));
                if p != Panel::Connectivity {
                    assert_eq!((f.w, f.h), (320, 240), "{}", p.slug());
                }
            }
            Picture::Chart(_) => panic!("a trajectory's panels are no chart"),
        }
    }
}

#[test]
fn the_cell_map_is_the_shared_jobs_frame_on_the_paper_page() {
    use legume_tui::map::points::Job;
    use legume_tui::map::theme::Theme;
    let mut t = data();
    // No arrows, labels or legend: only the points.
    t.edges.clear();
    let style = Style {
        labels: None,
        ..Style::default()
    };
    let m = map(&t, &style);
    let l = &t.layouts[0];
    let layers = [t.cell_layer(&l.x, &l.y, &style)];
    let mut job = Job::new(m.vp, Theme::PAPER.background, &layers);
    while !job.step(&layers) {}
    let img = render(&Picture::Map(Box::new(m))).unwrap();
    assert_eq!(img, job.finish().to_image(), "the shared job's pixels");
    let [r, g, b] = Theme::PAPER.background;
    assert_eq!(img.get_pixel(0, 0).0, [r, g, b, 255], "on the paper page");
}

#[test]
fn the_cell_map_colours_groups_from_sennas_palette_and_greys_the_rest() {
    use legume_tui::map::colour::{to_linear, Categorical};
    use legume_tui::map::theme::Theme;
    let mut t = data();
    t.types = (0..t.types.len())
        .map(|i| if i < 20 { "CT1" } else { "CT2" }.into())
        .collect();
    let style = Style {
        colouring: Colouring::Type,
        ..Style::default()
    };
    let (colours, _) = t.shown_colours(&style);
    let senna = |i| Some(Categorical::SENNA.srgb(i, 2).into());
    assert_eq!((colours[0], colours[39]), (senna(0), senna(1)));
    let l = &t.layouts[0];
    let layer = t.cell_layer(&l.x, &l.y, &style);
    let mark = layer.mark(39).unwrap();
    assert_eq!(mark.colour, to_linear(Categorical::SENNA.srgb(1, 2)));
    t.types[39] = enrichment::UNASSIGNED_LABEL.into();
    t.colours = Default::default();
    let layer = t.cell_layer(&l.x, &l.y, &style);
    assert_eq!(
        layer.mark(39).unwrap().colour,
        to_linear(Theme::PAPER.no_group),
        "a cell of no type in the page's grey"
    );
}

#[test]
fn labels_are_placed_at_the_type_medians() {
    let mut t = data();
    t.types = (0..t.types.len())
        .map(|i| if i < 20 { "CT1" } else { "CT2" }.into())
        .collect();
    let style = Style {
        colouring: Colouring::Type,
        ..Style::default()
    };
    let m = map(&t, &style);
    assert_eq!(labels(&m).len(), 2, "{:?}", labels(&m));
    let l = &t.layouts[0];
    let cells: Vec<usize> = (0..t.types.len()).collect();
    let medians = t.label_medians(LabelBy::Type, &l.x, &l.y, &cells);
    for p in &m.placed {
        let name = &m.names[p.key as usize];
        let ((x, y), _) = medians[name];
        let at = m.vp.to_px([x, y]);
        assert!(
            (p.anchor.0 - at.0).abs() < 1e-3 && (p.anchor.1 - at.1).abs() < 1e-3,
            "{name}: {:?} at {at:?}",
            p.anchor
        );
        let centre = 0.5 * (p.rect[0] + p.rect[2]);
        assert!((centre - at.0).abs() <= 1.0, "{name} centred on it");
    }
}

#[test]
fn the_cell_map_keeps_inside_the_data_when_zoomed_to_its_edge() {
    let t = data();
    let style = Style::default();
    let reach = View::reach(t.whole(Panel::Layout { k: 0 }).unwrap(), 320, 240);
    let mut v = View::default();
    for _ in 0..4 {
        v.zoom(true, reach);
    }
    for _ in 0..50 {
        v.pan(1.0, 1.0, reach);
    }
    let m = map_at(&t, &style, v);
    // Across, the window is narrower than the cells: stopped at their
    // right edge. Down, it is wider: centred on them.
    let [x0, y0, x1, y1] = m.vp.window();
    assert!((x1 - 39.0).abs() < 1e-3 && x0 > 0.0, "{:?}", m.vp.window());
    assert!((0.5 * (y0 + y1) - 3.0).abs() < 1e-3, "{:?}", m.vp.window());
}

#[test]
fn the_prior_edges_are_drawn_as_arrows_between_type_medians() {
    let mut t = data();
    let style = Style {
        labels: None,
        ..Style::default()
    };
    let with = render(&Picture::Map(Box::new(map(&t, &style)))).unwrap();
    t.edges.clear();
    let without = render(&Picture::Map(Box::new(map(&t, &style)))).unwrap();
    assert_ne!(with, without, "the arrow is drawn");
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
        t.title(Panel::Layout { k: 0 }, &Style::default()),
        "PHATE · pseudotime"
    );
    assert_eq!(t.slug(Panel::Layout { k: 0 }), "layout_phate");
}

#[test]
fn diffusion_pairs_cycle_past_the_trivial_component() {
    let t = data();
    assert_eq!(
        t.next_pair(1, 2, true, false),
        (1, 2),
        "with three components the only pair after the trivial one is (1, 2)"
    );
    assert_eq!(
        t.next_pair(1, 2, false, false),
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
        .export(
            Panel::Order,
            &Style::default(),
            View::default(),
            None,
            &Offsets::default(),
        )
        .unwrap();
    assert_eq!(first.base, dir.join("x.trajectory.order"));
    assert!(first.files.iter().all(|f| f.is_file()), "{:?}", first.files);
    let second = t
        .export(
            Panel::Order,
            &Style::default(),
            View::default(),
            None,
            &Offsets::default(),
        )
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
    let m = map(&t, &style);
    let mut names = labels(&m);
    names.sort_unstable();
    assert_eq!(names, ["A", "B"], "one label per type, on at medium");
    let height = |m: &CellMap| m.placed[0].rect[3] - m.placed[0].rect[1];
    let medium = height(&m);
    assert_eq!(style.cycle_labels(), "labels large · t for largest");
    assert!(height(&map(&t, &style)) > medium, "the size changed");
    style.cycle_labels();
    assert_eq!(style.cycle_labels(), "labels off · t shows them small");
    let off = map(&t, &style);
    assert!(off.placed.is_empty() && off.text.is_empty());
    assert_eq!(style.cycle_labels(), "labels small · t for medium");
}

#[test]
fn the_unassigned_type_gets_no_label() {
    let mut t = data();
    t.types[39] = enrichment::UNASSIGNED_LABEL.into();
    let m = map(&t, &Style::default());
    assert!(!texts(&m).contains(&enrichment::UNASSIGNED_LABEL));
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
    let by_time = map(&t, &style);
    assert!(by_time.legend.is_none(), "no legend for pseudotime");
    assert_eq!(
        style.cycle_colouring(&all),
        "coloured by cell type · c for component"
    );
    let by_type = map(&t, &style);
    assert_ne!(by_time.frame.to_image(), by_type.frame.to_image());
    assert!(by_type.legend.is_some(), "a legend");
    assert_eq!(
        texts(&by_type).iter().filter(|&&s| s == "A").count(),
        2,
        "A in the legend and as a label"
    );
    assert_eq!(
        t.title(Panel::Layout { k: 0 }, &style),
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
    v.pan(1.0, 0.0, ALL);
    assert_eq!(v, View::default(), "the whole scatter does not pan");
    v.zoom(false, ALL);
    assert_eq!(v.level, 0, "no further out than the whole");
    v.zoom(true, ALL);
    v.zoom(true, ALL);
    assert_eq!(v.factor(), 2.0);
    let half = 0.5 / v.factor();
    for _ in 0..50 {
        v.pan(1.0, -1.0, ALL);
    }
    assert!(
        (v.cx - (1.0 - half)).abs() < 1e-6,
        "stops at the right edge"
    );
    assert!((v.cy - half).abs() < 1e-6, "and at the bottom");
    for _ in 0..50 {
        v.zoom(true, ALL);
    }
    assert_eq!(v.level, MAX_ZOOM_LEVEL);
    assert_eq!(v.factor(), 64.0);
    // Zoomed ×2 into the top-left quarter of a square: the window is half
    // the fitted one, pushed back inside the data at its left and top.
    let whole = [0.0, -5.0, 10.0, 5.0];
    let vp = View {
        level: 2,
        cx: 0.25,
        cy: 0.75,
    }
    .viewport(whole, 100, 100);
    let [x0, _, x1, y1] = vp.window();
    let near = |a: f32, b: f32| (a - b).abs() < 1e-3;
    assert!(near(x0, 0.0) && near(y1, 5.0), "{:?}", vp.window());
    assert!(near(x1 - x0, 10.0 / 0.92 / 2.0), "{:?}", vp.window());
    let all = View::default().viewport(whole, 100, 100);
    assert_eq!((all.cx, all.cy), (5.0, 0.0), "whole, centred");
}

/// What a frame showing exactly the whole extent reaches on each axis.
const ALL: [f32; 2] = [1.0, 1.0];

/// The view's centre as the frame draws it, in shares of `whole`.
fn drawn_centre(v: View, whole: Bounds, w: usize, h: usize) -> (f32, f32) {
    let vp = v.viewport(whole, w, h);
    let [x0, y0, x1, y1] = whole;
    ((vp.cx - x0) / (x1 - x0), (vp.cy - y0) / (y1 - y0))
}

#[test]
fn the_view_is_where_the_frame_draws_it_so_every_pan_moves_the_picture() {
    // A square of data in a 3:2 frame: across, the frame shows more than
    // the data at ×1.4, so it does not pan that way.
    let whole = [0.0, 0.0, 10.0, 10.0];
    let (w, h) = (1200, 800);
    let reach = View::reach(whole, w, h);
    assert!(reach[0] > reach[1], "{reach:?}");
    let mut v = View::default();
    v.zoom(true, reach);
    for _ in 0..3 {
        v.pan(1.0, 0.0, reach);
    }
    assert_eq!(v.cx, 0.5, "across, the data is all on screen");
    // ×2: as far right as the frame goes, the view is what is drawn, and
    // the first step back moves the picture.
    v.zoom(true, reach);
    for _ in 0..50 {
        v.pan(1.0, 1.0, reach);
    }
    let (dx, dy) = drawn_centre(v, whole, w, h);
    assert!((v.cx - dx).abs() < 1e-5 && (v.cy - dy).abs() < 1e-5);
    assert!(v.cx < 1.0 - 0.5 / v.factor(), "{v:?}");
    let before = drawn_centre(v, whole, w, h);
    v.pan(-1.0, 0.0, reach);
    let after = drawn_centre(v, whole, w, h);
    assert!(after.0 < before.0 - 1e-3, "{before:?} {after:?}");
    // Zooming out keeps it so.
    v.zoom(false, reach);
    let (dx, dy) = drawn_centre(v, whole, w, h);
    assert!((v.cx - dx).abs() < 1e-5 && (v.cy - dy).abs() < 1e-5);
}

#[test]
fn a_zoomed_export_shows_the_data_the_screen_shows() {
    let whole = [0.0, 0.0, 10.0, 10.0];
    let v = View {
        level: 3,
        cx: 0.4,
        cy: 0.6,
    };
    for (w, h) in [(1200, 800), (600, 900), (640, 480)] {
        let size = export_size(Some((w, h)));
        let screen = v.viewport(whole, w as usize, h as usize).window();
        let paper = v.viewport(whole, size.w as usize, size.h as usize).window();
        for (a, b) in screen.iter().zip(&paper) {
            assert!((a - b).abs() < 1e-2, "{w}×{h}: {screen:?} {paper:?}");
        }
    }
    let size = export_size(None);
    assert_eq!(size.h, size.w * 3 / 4, "4:3 with no screen to follow");
}

#[test]
fn zooming_in_and_back_out_shows_the_whole_scatter() {
    let mut v = View::default();
    for _ in 0..3 {
        v.zoom(true, ALL);
    }
    v.pan(1.0, 1.0, ALL);
    for _ in 0..3 {
        v.zoom(false, ALL);
    }
    assert!(v.is_whole());
    assert_eq!(v.factor(), 1.0);
    assert_eq!(v, View::default(), "and centred again");
}

#[test]
fn a_zoomed_scatter_labels_only_the_types_on_screen() {
    let t = data();
    let style = Style::default();
    let whole = map(&t, &style);
    let mut names = labels(&whole);
    names.sort_unstable();
    assert_eq!(names, ["A", "B"]);
    // The left quarter holds only type A's cells.
    let left = map_at(
        &t,
        &style,
        View {
            level: 4,
            cx: 0.125,
            cy: 0.5,
        },
    );
    assert_eq!(labels(&left), ["A"], "B is off screen");
    assert_ne!(whole.frame.to_image(), left.frame.to_image());
    assert_eq!(
        map(&t, &style).frame.to_image(),
        whole.frame.to_image(),
        "the medians kept draw the same"
    );
    // The order panel ignores the view.
    let zoomed = View {
        level: 2,
        ..View::default()
    };
    let svg = |p: Picture| match p {
        Picture::Svg(f) => f.svg,
        Picture::Map(_) | Picture::Chart(_) => panic!("the order panel is an SVG"),
    };
    assert_eq!(
        svg(t
            .figure(Panel::Order, SIZE, &style, zoomed, &Offsets::default())
            .unwrap()),
        svg(draw(&t, Panel::Order, &style))
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
    assert_eq!(legend_rows(40, 20), None, "no room for two rows");
    let small = legend_rows(200, 20).unwrap();
    assert!(small <= 5, "half of 200 px holds few rows: {small}");
    assert_eq!(legend_rows(4000, 20), Some(LEGEND_MAX));
    // Drawn: 25 types list LEGEND_MAX - 1 and count the rest.
    let mut t = data();
    t.types = (0..t.types.len())
        .map(|i| format!("CT{}", i % 25 + 1).into())
        .collect();
    let style = Style {
        colouring: Colouring::Type,
        labels: None,
        ..Style::default()
    };
    let size = Size {
        w: 1600,
        h: 1200,
        text_px: 16.0,
    };
    let m = t
        .cell_map(
            Panel::Layout { k: 0 },
            size,
            &style,
            View::default(),
            true,
            &Offsets::default(),
        )
        .unwrap();
    let text = texts(&m);
    assert_eq!(text.len(), LEGEND_MAX, "{text:?}");
    assert_eq!(text.last(), Some(&"+6 more"));
    let thumb = t
        .cell_map(
            Panel::Layout { k: 0 },
            Size {
                w: 60,
                h: 40,
                text_px: 16.0,
            },
            &style,
            View::default(),
            true,
            &Offsets::default(),
        )
        .unwrap();
    assert!(thumb.legend.is_none() && thumb.text.is_empty());
}

#[test]
fn a_label_near_the_edge_stays_inside_the_figure() {
    let mut t = data();
    // One cell of a long-named type at the far right, at the top, clear
    // of the other type's label.
    t.types = (0..t.types.len())
        .map(|i| {
            if i == 39 {
                "CT1_with_a_long_name"
            } else {
                "CT2"
            }
            .into()
        })
        .collect();
    t.layouts[0].y[39] = 6.0;
    let style = Style {
        colouring: Colouring::Type,
        ..Style::default()
    };
    let m = map(&t, &style);
    let placed = m
        .placed
        .iter()
        .find(|p| m.names[p.key as usize].as_ref() == "CT1_with_a_long_name")
        .expect("kept, not dropped off the edge");
    assert!(
        placed.rect[0] >= 0.0 && placed.rect[2] <= 320.0,
        "{:?}",
        placed.rect
    );
    let at = m.vp.to_px([39.0, 6.0]);
    assert!(placed.anchor.0 < at.0, "moved in from its cell");
}

#[test]
fn either_axis_of_the_diffusion_map_steps_and_skips_the_other() {
    let mut t = data();
    t.diffusion = Some(Mat::from_element(t.types.len(), 5, 0.0));
    // y steps past x; x steps past y; both wrap past the trivial first.
    assert_eq!(t.next_pair(1, 2, true, false), (1, 3));
    assert_eq!(t.next_pair(1, 4, true, false), (1, 2));
    assert_eq!(t.next_pair(1, 3, true, true), (2, 3));
    assert_eq!(t.next_pair(2, 3, true, true), (4, 3));
    assert_eq!(t.next_pair(1, 3, false, true), (4, 3));
}

#[test]
fn cells_draw_grey_first_then_by_pseudotime_low_to_high() {
    let mut t = data();
    t.pseudotime = vec![0.9, f32::NAN, 0.1, 0.5];
    t.types.truncate(4);
    let order = t.draw_order(&[0, 1, 2, 3], Colouring::Pseudotime, &[]);
    assert_eq!(order, vec![1, 2, 3, 0]);
}

#[test]
fn layouts_come_from_the_run_a_round_was_made_from_too() {
    let root = tempfile::tempdir().unwrap();
    let at = |n: &str| root.path().join(n);
    std::fs::write(
        at("r.senna.json"),
        r#"{"version": 2, "kind": "topic", "prefix": "r", "layout": {"cell_coords": "r.umap.parquet",
            "current": "umap", "methods": {"umap": {"cell_coords": "r.umap.parquet"},
            "phate": {"cell_coords": "r.phate.parquet"}}}}"#,
    )
    .unwrap();
    // A round made before PHATE was added: it copied only the UMAP.
    std::fs::write(
        at("r.L1.senna.json"),
        r#"{"version": 2, "kind": "topic", "prefix": "r", "annotate": {"source": "r.senna.json"},
            "layout": {"cell_coords": "r.umap.parquet", "current": "umap"}}"#,
    )
    .unwrap();
    let l = crate::manifest::run::load(&at("r.L1.senna.json").to_string_lossy()).unwrap();
    let (found, run) = layouts_along(&l.manifest, &l.dir, &l.file);
    let methods: Vec<Option<&str>> = found.iter().map(|(m, _)| m.as_deref()).collect();
    assert_eq!(methods, vec![Some("phate"), Some("umap")]);
    assert!(found[0].1.ends_with("r.phate.parquet"));
    assert!(
        run.ends_with("r.senna.json"),
        "senna layout goes on the run"
    );
}

#[test]
fn a_missing_umap_or_phate_names_the_senna_command() {
    let run = Path::new("/x/r.senna.json");
    let has = |m: &[&'static str]| m.iter().map(|&s| Some(s)).collect::<Vec<_>>();
    assert_eq!(layout_hint(has(&["phate", "umap"]), run), None);
    let one = layout_hint(has(&["umap"]), run).unwrap();
    assert!(
        one.contains("no PHATE") && one.contains("`senna layout phate --from r.senna.json`"),
        "{one}"
    );
    let both = layout_hint([], run).unwrap();
    assert!(
        both.contains("no UMAP or PHATE") && both.contains("senna layout {umap|phate}"),
        "{both}"
    );
}

#[test]
fn a_layout_is_read_onto_every_cell_even_with_a_name_twice() {
    use legume_numeric::matrix::traits::IoOps as _;
    // Merged samples can repeat a barcode: three cells, two names.
    let cells: Vec<Box<str>> = ["a", "b", "a"].iter().map(|s| Box::from(*s)).collect();
    let index = cell_index(&cells);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("layout.parquet");
    let rows: Vec<Box<str>> = ["a", "b"].iter().map(|s| Box::from(*s)).collect();
    let cols: Vec<Box<str>> = ["x", "y"].iter().map(|s| Box::from(*s)).collect();
    let m = Mat::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
    m.to_parquet_with_names(
        &path.to_string_lossy(),
        (Some(&rows), Some("cell")),
        Some(&cols),
    )
    .unwrap();
    let (x, y) = read_layout(&path.to_string_lossy(), &index, cells.len()).unwrap();
    assert_eq!(x.len(), 3, "one per cell");
    assert_eq!((x[1], y[1]), (3.0, 4.0));
    assert_eq!((x[0], y[0]), (1.0, 2.0), "each cell of the name");
    assert_eq!((x[2], y[2]), (1.0, 2.0));
}

#[test]
fn a_table_that_repeats_a_name_gives_each_cell_its_own_row() {
    let cells: Vec<Box<str>> = ["a", "b", "a"].iter().map(|s| Box::from(*s)).collect();
    let index = cell_index(&cells);
    // The table names its rows as the cells do: the k-th `a` row is the
    // k-th `a` cell's.
    let rows: Vec<Box<str>> = ["a", "b", "a"].iter().map(|s| Box::from(*s)).collect();
    let m = Mat::from_row_slice(3, 1, &[1.0, 2.0, 9.0]);
    let out = aligned(&rows, &m, &index, 3);
    assert_eq!(
        out.column(0).iter().copied().collect::<Vec<_>>(),
        [1.0, 2.0, 9.0]
    );
}

/// The figures of the round `l`.
fn of_round(l: &crate::manifest::run::Loaded) -> CellFigures {
    let cells = crate::manifest::rounds::read_cells(l).unwrap();
    let score = crate::manifest::rounds::read_score_cells(l, &cells).unwrap();
    CellFigures::of_round(l, &cells, score)
}

/// A round on four cells in two clusters, labelled A and B, with a UMAP.
fn round_with_layout(dir: &Path) -> crate::manifest::run::Loaded {
    round(dir, true)
}

/// A round on four cells in two clusters, labelled A and B; a UMAP when
/// `layout`.
pub(crate) fn round(dir: &Path, layout: bool) -> crate::manifest::run::Loaded {
    use legume_numeric::matrix::traits::IoOps as _;
    let at = |n: &str| dir.join(n).to_string_lossy().into_owned();
    let names: Vec<Box<str>> = ["c0", "c1", "c2", "c3"]
        .iter()
        .map(|s| Box::from(*s))
        .collect();
    crate::manifest::rounds::write_clusters(
        &at("r.clusters.parquet"),
        &names,
        &[Some(0), Some(0), Some(1), Some(1)],
    )
    .unwrap();
    std::fs::write(
        at("r.argmax.tsv"),
        "cell\tcell_type\tprobability\nc0\tA\t0.9\nc1\tA\t0.8\nc2\tB\t0.7\nc3\tunassigned\t0.1\n",
    )
    .unwrap();
    let cols: Vec<Box<str>> = ["x", "y"].iter().map(|s| Box::from(*s)).collect();
    Mat::from_row_slice(4, 2, &[0.0, 0.0, 1.0, 0.0, 5.0, 5.0, 6.0, 5.0])
        .to_parquet_with_names(
            &at("r.umap.parquet"),
            (Some(&names), Some("cell")),
            Some(&cols),
        )
        .unwrap();
    let layout = if layout {
        r#", "layout": {"cell_coords": "r.umap.parquet", "current": "umap"}"#
    } else {
        ""
    };
    std::fs::write(
        at("r.senna.json"),
        format!(
            r#"{{"version": 2, "kind": "bge", "prefix": "r",
                "cluster": {{"clusters": "r.clusters.parquet"}},
                "annotate": {{"argmax": "r.argmax.tsv"}}{layout}}}"#
        ),
    )
    .unwrap();
    crate::manifest::run::load(&at("r.senna.json")).unwrap()
}

#[test]
fn a_rounds_cells_come_with_its_labels_clusters_and_layout() {
    let dir = tempfile::tempdir().unwrap();
    let f = of_round(&round_with_layout(dir.path()));
    assert_eq!(
        f.types.iter().map(AsRef::as_ref).collect::<Vec<_>>(),
        ["A", "A", "B", UNASSIGNED_LABEL]
    );
    assert_eq!(f.cluster, vec![0, 0, 1, 1]);
    assert_eq!((f.layouts.len(), f.layouts[0].x[2]), (1, 5.0));
    assert_eq!(
        f.colourings(),
        vec![Colouring::Type, Colouring::Cluster],
        "no pseudotime"
    );
    assert_eq!(
        f.panels(),
        vec![Panel::Layout { k: 0 }],
        "no trajectory panels"
    );
    let fig = draw(&f, Panel::Layout { k: 0 }, &Style::opening(&f));
    assert!(matches!(fig, Picture::Map(_)));
}

#[test]
fn a_round_without_a_layout_has_no_figures_but_names_the_senna_command() {
    let dir = tempfile::tempdir().unwrap();
    let f = of_round(&round(dir.path(), false));
    assert!(f.layouts.is_empty());
    assert!(
        f.layout_hint
            .as_deref()
            .is_some_and(|h| h.contains("senna layout")),
        "{:?}",
        f.layout_hint
    );
}

#[test]
fn a_focused_cluster_keeps_its_colours_and_draws_over_the_greyed_rest() {
    let dir = tempfile::tempdir().unwrap();
    let f = of_round(&round_with_layout(dir.path()));
    let style = Style {
        focus: Some(0),
        ..Style::opening(&f)
    };
    let (colours, _) = f.shown_colours(&style);
    assert!(
        colours[0].is_some() && colours[1].is_some(),
        "cluster 0 in colour"
    );
    assert!(
        colours[2].is_none() && colours[3].is_none(),
        "the rest grey"
    );
    let order = f.draw_order(&[0, 1, 2, 3], style.colouring, &colours);
    assert_eq!(&order[..2], &[2, 3], "grey first, under the focus");
    assert!(
        f.title(Panel::Layout { k: 0 }, &style).contains("K0"),
        "the title names the focus"
    );
}

#[test]
fn a_rounds_figures_export_as_cells_not_trajectory() {
    let dir = tempfile::tempdir().unwrap();
    let f = of_round(&round_with_layout(dir.path()));
    assert_eq!(f.kind(), "cells");
    assert_eq!(data().kind(), "trajectory");
    let e = f
        .export(
            Panel::Layout { k: 0 },
            &Style::opening(&f),
            View::default(),
            None,
            &Offsets::default(),
        )
        .unwrap();
    assert!(
        e.base.to_string_lossy().contains("r.cells.layout_umap"),
        "{}",
        e.base.display()
    );
    // The cell map is a PDF only, as senna's: its raster with the text set
    // over it.
    let pdf = PathBuf::from(format!("{}.pdf", e.base.display()));
    assert_eq!(e.files, vec![pdf.clone()]);
    assert!(std::fs::read(&pdf).unwrap().starts_with(b"%PDF"));
    assert!(!Path::new(&format!("{}.svg", e.base.display())).exists());
}

#[test]
fn coloured_by_cluster_the_labels_name_the_clusters() {
    let dir = tempfile::tempdir().unwrap();
    let f = of_round(&round_with_layout(dir.path()));
    let by = |colouring| {
        let style = Style {
            colouring,
            ..Style::opening(&f)
        };
        let m = map(&f, &style);
        texts(&m)
            .into_iter()
            .map(str::to_string)
            .collect::<Vec<_>>()
    };
    let has = |t: &[String], s: &str| t.iter().filter(|x| *x == s).count();
    let types = by(Colouring::Type);
    assert!(has(&types, "A") > 0 && has(&types, "K0") == 0, "{types:?}");
    let clusters = by(Colouring::Cluster);
    assert_eq!(has(&clusters, "A"), 0, "no type names over the clusters");
    assert_eq!(has(&clusters, "K0"), 2, "legend and label");
}

#[test]
fn with_a_cluster_in_focus_the_legend_lists_only_its_colours() {
    let dir = tempfile::tempdir().unwrap();
    let f = of_round(&round_with_layout(dir.path()));
    let style = Style {
        labels: None,
        focus: Some(1),
        ..Style::opening(&f)
    };
    let m = map(&f, &style);
    assert!(texts(&m).contains(&"B"), "cluster 1's label");
    assert!(!texts(&m).contains(&"A"), "not cluster 0's");
}

#[test]
fn a_round_drawn_by_pseudotime_it_lacks_is_grey_not_a_panic() {
    let dir = tempfile::tempdir().unwrap();
    let f = of_round(&round_with_layout(dir.path()));
    let fig = draw(&f, Panel::Layout { k: 0 }, &Style::default());
    assert!(matches!(fig, Picture::Map(_)));
}

/// [`round`] with the projection score's per-cell table beside it: c0, c1
/// called B with margins 0.1 and 0.9, c2 called A; c3 has no row.
fn round_with_score(dir: &Path) -> crate::manifest::run::Loaded {
    use legume_numeric::matrix::parquet::{write_named_table, Column};
    let mut l = round(dir, true);
    let cells: Vec<Box<str>> = ["c0", "c1", "c2"].map(Box::from).to_vec();
    let call: Vec<Box<str>> = ["B", "B", "A"].map(Box::from).to_vec();
    let table = dir.join("r.marker_score.parquet");
    write_named_table(
        &table.to_string_lossy(),
        "cell",
        &cells,
        &[
            ("cluster".into(), Column::I32(&[0, 0, 1])),
            ("call".into(), Column::Str(&call)),
            ("margin".into(), Column::F32(&[0.1, 0.9, 0.5])),
        ],
    )
    .unwrap();
    l.manifest.annotate.marker_score = Some("r.marker_score.parquet".into());
    l
}

#[test]
fn a_round_with_the_score_colours_by_its_calls_and_their_margins() {
    let dir = tempfile::tempdir().unwrap();
    let f = of_round(&round_with_score(dir.path()));
    let colourings = f.colourings();
    assert!(colourings.contains(&Colouring::Call), "{colourings:?}");
    assert!(colourings.contains(&Colouring::Margin), "{colourings:?}");

    let (calls, legend) = f.cell_colours(Colouring::Call);
    assert_eq!(calls[0], calls[1], "both called B");
    assert_ne!(calls[0], calls[2], "A apart");
    assert!(
        calls[3].is_none(),
        "a cell the score has no row for is grey"
    );
    let names: Vec<&str> = legend.iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(names, ["A", "B"]);

    let (margins, _) = f.cell_colours(Colouring::Margin);
    assert!(margins[..3].iter().all(Option::is_some) && margins[3].is_none());
    assert_ne!(margins[0], margins[1], "a wide margin from a narrow one");

    let style = Style {
        colouring: Colouring::Call,
        ..Style::opening(&f)
    };
    let m = map(&f, &style);
    assert!(
        labels(&m).contains(&"B"),
        "coloured by call, the labels name the calls"
    );
}

#[test]
fn a_round_without_the_score_has_no_call_or_margin_colouring() {
    let dir = tempfile::tempdir().unwrap();
    let colourings = of_round(&round_with_layout(dir.path())).colourings();
    assert!(!colourings.contains(&Colouring::Call) && !colourings.contains(&Colouring::Margin));
}

#[test]
fn a_drag_pans_the_view_and_the_wheel_zooms_about_the_pointer() {
    let t = data();
    let whole = t.whole(Panel::Layout { k: 0 }).unwrap();
    let (w, h) = (320, 240);
    let mut v = View::default();
    v.pan_px(10.0, 0.0, whole, w, h);
    assert_eq!(v, View::default(), "the whole scatter does not pan");
    // A step in, the data under the pointer kept under it.
    let at = (80.0, 120.0);
    let before = v.viewport(whole, w, h).to_data(at.0, at.1);
    v.zoom_at(true, at, whole, w, h);
    assert_eq!(v.level, 1);
    let after = v.viewport(whole, w, h).to_data(at.0, at.1);
    assert!((before[0] - after[0]).abs() < 1e-3, "{before:?} {after:?}");
    // Dragged right, the picture follows the pointer.
    let cell = [10.0, 3.0];
    let was = v.viewport(whole, w, h).to_px(cell);
    v.pan_px(10.0, 0.0, whole, w, h);
    let now = v.viewport(whole, w, h).to_px(cell);
    assert!((now.0 - was.0 - 10.0).abs() < 1e-2, "{was:?} {now:?}");
    v.zoom_at(false, at, whole, w, h);
    assert_eq!(v, View::default(), "a step out shows it whole again");
}

#[test]
fn a_moved_label_sits_its_offset_in_lines_from_its_cells() {
    let t = data();
    let style = Style {
        colouring: Colouring::Type,
        ..Style::default()
    };
    let at = |m: &CellMap| {
        *m.placed
            .iter()
            .find(|p| m.names[p.key as usize].as_ref() == "A")
            .expect("A placed")
    };
    let still = map(&t, &style);
    assert!(still.line > 0.0);
    let mut o = Offsets::default();
    o.nudge("A".to_string(), (0.0, 2.0));
    let moved = t
        .cell_map(
            Panel::Layout { k: 0 },
            SIZE,
            &style,
            View::default(),
            true,
            &o,
        )
        .unwrap();
    let (a, b) = (at(&still), at(&moved));
    assert!(b.moved && !a.moved);
    let dy = b.rect[1] - a.rect[1];
    assert!(
        (dy - 2.0 * moved.line).abs() < 2.0,
        "{dy} for {}",
        moved.line
    );
}

#[test]
fn the_cell_under_a_pixel_is_the_nearest_within_reach() {
    let t = data();
    let vp = View::default().viewport(t.whole(Panel::Layout { k: 0 }).unwrap(), 320, 240);
    let (px, py) = vp.to_px([5.0, 5.0]);
    let panel = Panel::Layout { k: 0 };
    assert_eq!(t.cell_at(panel, &vp, px + 1.0, py, 4.0), Some(5));
    assert_eq!(t.cell_at(panel, &vp, px, py - 60.0, 4.0), None, "none near");
}
