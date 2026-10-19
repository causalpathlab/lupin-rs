use super::*;
use crate::figures::Layout;
use crate::trajectory::edges::{EdgeRow, Verdict};
use legume_numeric::matrix::dense_mat_io::Mat;
use ratatui::crossterm::event::KeyCode;

/// What `p` has made for the terminal: its picture protocol, and whether a
/// terminal picture and a rasterised figure are cached.
pub(crate) fn made(p: &FigurePane) -> (ratatui_image::picker::ProtocolType, bool, bool) {
    let rendered = p.rendered.borrow().is_some();
    (
        p.picker.protocol_type(),
        rendered,
        p.image.borrow().is_some(),
    )
}

/// Two layouts, a diffusion map and one edge, on placeholder types.
pub(crate) fn pane() -> FigurePane {
    let n = 30;
    let x: Vec<f32> = (0..n).map(|i| i as f32).collect();
    let y: Vec<f32> = (0..n).map(|i| (i % 5) as f32).collect();
    let mut d = Mat::zeros(n, 3);
    for i in 0..n {
        d[(i, 1)] = x[i];
        d[(i, 2)] = y[i];
    }
    let data = CellFigures {
        manifest: "run.senna.json".into(),
        pseudotime: (0..n).map(|i| i as f32 / n as f32).collect(),
        types: (0..n)
            .map(|i| if i < 15 { "CT1" } else { "CT2" }.into())
            .collect(),
        component: vec![0; n],
        lineage: vec![0; n],
        diffusion: Some(d),
        layouts: vec![
            Layout {
                method: Some("phate".into()),
                x: x.clone(),
                y: y.clone(),
            },
            Layout {
                method: Some("umap".into()),
                x: y,
                y: x,
            },
        ],
        edges: vec![EdgeRow {
            a: "CT1".into(),
            b: "CT2".into(),
            connectivity: 0.5,
            in_prior: true,
            verdict: Some(Verdict::Supported),
            order_agreement: 0.9,
        }],
        ..CellFigures::default()
    };
    let mut p = FigurePane::new(data, &Picker::halfblocks());
    p.shown = true;
    p
}

#[test]
fn the_grid_holds_every_figure_and_opens_the_chosen_one() {
    let mut p = pane();
    p.open_grid();
    let g = p.grid.as_ref().unwrap();
    assert_eq!(
        g.tiles,
        vec![
            Panel::Layout { k: 0 },
            Panel::Layout { k: 1 },
            Panel::Diffusion { x: 1, y: 2 },
            Panel::Order,
            Panel::Connectivity,
        ]
    );
    assert_eq!(g.sel, 0, "it opens on the figure on screen");
    assert_eq!(g.cols(), 3);
    let g = p.grid.as_mut().unwrap();
    assert!(g.step(KeyCode::Down));
    assert_eq!(g.sel, 3);
    assert!(g.step(KeyCode::Down), "no row below: it stays");
    assert_eq!(g.sel, 3);
    assert!(g.step(KeyCode::Up));
    assert_eq!(g.sel, 0);
    assert!(g.step(KeyCode::Up), "no row above: it stays");
    assert_eq!(g.sel, 0);
    assert!(g.step(KeyCode::Left), "left wraps");
    assert_eq!(g.sel, 4);
    assert!(!g.step(KeyCode::Enter));
    // Open the UMAP: it takes the scatter's place.
    g.sel = 1;
    p.open_tile();
    assert!(p.grid.is_none());
    assert_eq!(p.current(), Panel::Layout { k: 1 });
    // And the order panel.
    p.open_grid();
    p.grid.as_mut().unwrap().sel = 3;
    p.open_tile();
    assert_eq!(p.current(), Panel::Order);
    assert_eq!(p.panels[0], Panel::Layout { k: 1 }, "the scatter is kept");
}

#[test]
fn tiles_are_drawn_once_and_again_after_a_setting_changes() {
    let mut p = pane();
    p.open_grid();
    let area = Rect::new(0, 0, 20, 8);
    let all: Vec<(usize, Rect)> = (0..5).map(|i| (i, area)).collect();
    p.draw_tiles(&all);
    assert_eq!(p.tiles.borrow().len(), 5);
    for i in 0..5 {
        assert!(p.tile(i, area).unwrap().is_ok(), "tile {i}");
    }
    assert!(p.tile(5, area).is_none());
    assert_eq!(p.tiles.borrow().len(), 5, "drawn once");
    p.style.cycle_labels();
    p.draw_tiles(&all);
    assert_eq!(
        p.tiles.borrow().len(),
        8,
        "the three scatters again; the order and connectivity are kept"
    );
    p.grid = None;
    p.open_grid();
    p.draw_tiles(&all);
    assert_eq!(p.tiles.borrow().len(), 8, "kept while the grid was closed");
}

#[test]
fn only_the_scatter_zooms_and_a_new_panel_shows_it_whole() {
    let mut p = pane();
    assert!(p.zoom(true).starts_with("zoom ×1.4"));
    assert!(!p.view.is_whole());
    assert_eq!(p.pan(1.0, 0.0), "");
    assert!(p.view.cx > 0.5);
    p.next_layout();
    assert!(p.view.is_whole(), "m shows the new layout whole");
    p.zoom(true);
    assert_eq!(p.reset_view(), "the whole scatter");
    assert!(p.view.is_whole());
    p.zoom(true);
    p.next_panel();
    assert_eq!(p.current(), Panel::Order);
    assert!(p.view.is_whole(), "v shows the next panel whole");
    assert!(p.zoom(true).starts_with("only the scatter"));
    assert!(p.view.is_whole());
}

#[test]
fn the_diffusion_map_keeps_its_pair_across_layouts() {
    let mut p = pane();
    p.next_layout();
    p.next_layout();
    assert_eq!(p.current(), Panel::Diffusion { x: 1, y: 2 });
    p.pair = (2, 1);
    p.next_layout();
    p.next_layout();
    p.next_layout();
    assert_eq!(
        p.current(),
        Panel::Diffusion { x: 2, y: 1 },
        "m comes back to it"
    );
    p.open_grid();
    assert!(p
        .grid
        .as_ref()
        .unwrap()
        .tiles
        .contains(&Panel::Diffusion { x: 2, y: 1 }));
}

#[test]
fn a_scatter_opened_from_the_grid_shows_it_whole() {
    let mut p = pane();
    p.zoom(true);
    assert!(!p.view.is_whole());
    // To the order panel through the grid, then back to the same scatter.
    p.open_grid();
    let order = p
        .grid
        .as_ref()
        .unwrap()
        .tiles
        .iter()
        .position(|&t| t == Panel::Order);
    p.grid.as_mut().unwrap().sel = order.unwrap();
    p.open_tile();
    assert_eq!(p.current(), Panel::Order);
    p.open_grid();
    p.grid.as_mut().unwrap().sel = 0;
    p.open_tile();
    assert_eq!(p.current(), Panel::Layout { k: 0 });
    assert!(p.view.is_whole(), "as its tile shows it");
}

#[test]
fn the_grid_names_the_tile_it_would_export() {
    let mut p = pane();
    assert_eq!(p.selected_tile(), None);
    p.open_grid();
    p.grid.as_mut().unwrap().sel = 1;
    assert_eq!(p.selected_tile(), Some(Panel::Layout { k: 1 }));
}

#[test]
fn a_new_terminal_gets_new_pictures_and_new_figures_only_when_it_differs() {
    let mut p = pane();
    let area = Rect::new(0, 0, 30, 10);
    assert!(p.protocol(area).is_ok());
    assert!(p.rendered.borrow().is_some() && p.image.borrow().is_some());
    // The same kind of terminal: the picture is made again, from the figure
    // already drawn.
    p.repick(&Picker::halfblocks(), false);
    assert!(p.rendered.borrow().is_none(), "the picture goes");
    assert!(p.image.borrow().is_some(), "the figure stays");
    assert!(p.protocol(area).is_ok());
    // Another protocol or cell size: the figure is drawn again too.
    #[allow(deprecated)]
    let mut kitty = Picker::from_fontsize((10, 20).into());
    kitty.set_protocol_type(ratatui_image::picker::ProtocolType::Kitty);
    p.repick(&kitty, true);
    assert!(p.rendered.borrow().is_none() && p.image.borrow().is_none());
    assert_eq!(
        p.picker.protocol_type(),
        ratatui_image::picker::ProtocolType::Kitty
    );
    assert!(p.protocol(area).is_ok());
}

#[test]
fn the_axis_keys_say_where_they_work_and_what_they_chose() {
    let mut p = pane();
    assert!(
        p.step_pair(true, true).contains("diffusion map"),
        "only there"
    );
    while !matches!(p.current(), Panel::Diffusion { .. }) {
        p.next_layout();
    }
    let said = p.step_pair(true, false);
    assert_eq!(said, p.title(p.current()), "the new pair, as titled");
}

#[test]
fn a_file_just_listed_stays_selected_when_the_strip_resorts() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().canonicalize().unwrap();
    let write = |base: &str| {
        let pdf = dir.join(format!("{base}.pdf"));
        std::fs::write(&pdf, "%PDF-1.4 fake").unwrap();
        pdf
    };
    let listed = write("run.trajectory.a");
    write("run.trajectory.b");
    let mut x = Exports {
        prefix: dir.join("run.trajectory").to_string_lossy().into_owned(),
        gallery: Gallery::open(&dir.join(super::super::gallery::DIR)),
        open: true,
        sel: 0,
        statuses: Vec::new(),
        not_listed: Vec::new(),
        delete_armed: None,
        fits: std::cell::Cell::new(true),
    };
    x.gallery
        .add(&[listed], "a", "a", &dir.join("run.senna.json"))
        .unwrap();
    x.refresh();
    assert_eq!((x.gallery.entries.len(), x.not_listed.len()), (1, 1));
    // The not-listed row is the second; listing it puts it first.
    x.sel = 1;
    x.adopt(&dir.join("run.senna.json")).unwrap();
    let on = &x.gallery.entries[x.sel].path;
    assert!(on.ends_with("run.trajectory.b.pdf"), "{}", on.display());
}

#[test]
fn the_axes_keys_say_when_there_is_no_diffusion_map_to_step() {
    let mut p = pane();
    Arc::get_mut(&mut p.data).unwrap().diffusion = None;
    let said = p.step_pair(true, false);
    assert!(!said.contains("m reaches it"), "{said}");
    assert!(said.contains("no diffusion map"), "{said}");
    assert!(pane().step_pair(true, false).contains("m reaches it"));
}

#[test]
fn the_strips_keys_say_why_when_it_has_no_room() {
    let mut p = pane();
    p.exports.not_listed = vec!["run.trajectory.extra.pdf".into()];
    assert_eq!(p.key(KeyCode::Char('f')), FigureKey::Took);
    p.exports.fits.set(false);
    for key in [
        KeyCode::Enter,
        KeyCode::Char('d'),
        KeyCode::Char('D'),
        KeyCode::Char('M'),
    ] {
        match p.key(key) {
            FigureKey::Say(s) => assert!(s.contains("no room"), "{key:?}: {s}"),
            k => panic!("{key:?}: {k:?}"),
        }
    }
    assert!(p.exports.not_listed.len() == 1, "nothing listed or removed");
}

#[test]
fn a_pan_the_frame_cannot_show_says_so_and_the_export_follows_the_frame() {
    let mut p = pane();
    let area = Rect::new(0, 0, 30, 10);
    assert!(p.protocol(area).is_ok());
    let (cw, ch) = p.cell_px();
    assert_eq!(p.frame(), (30 * cw, 10 * ch), "the pane's pixels");
    p.zoom(true);
    // The first layout is wide and flat: its height is all on screen.
    assert_eq!(p.pan(0.0, 1.0), "at the edge of the scatter");
    assert_eq!(p.view.cy, 0.5);
    assert_eq!(p.pan(1.0, 0.0), "");
}

#[test]
fn a_quick_frame_shows_the_figure_last_drawn_while_the_new_one_is_drawn_apart() {
    let mut p = pane();
    let area = Rect::new(0, 0, 30, 10);
    assert!(
        matches!(p.drawn_so_far(area), (None, true)),
        "nothing drawn yet"
    );
    // A full frame waits for it.
    assert!(p.protocol(area).is_ok());
    assert!(matches!(p.drawn_so_far(area), (Some(_), false)));
    // A new look: the old drawing stands in until the new one is in.
    p.style.cycle_labels();
    assert!(matches!(p.drawn_so_far(area), (Some(_), true)));
    let start = std::time::Instant::now();
    while p.drawn_so_far(area).1 {
        assert!(start.elapsed().as_secs() < 30, "never drawn");
        std::thread::sleep(std::time::Duration::from_millis(5));
    }
    assert!(p.drawn_so_far(area).0.is_some_and(|d| d.is_ok()));
}

#[test]
fn the_scatter_on_screen_is_the_cell_map_on_the_paper_page() {
    let p = pane();
    let area = Rect::new(0, 0, 30, 10);
    let look = p.look(p.current(), p.view);
    let (cw, ch) = p.cell_px();
    let img = FigurePane::rasterised(&p.data, look, &p.style, (cw, ch), area).unwrap();
    assert_eq!((img.width(), img.height()), (30 * cw, 10 * ch));
    let [r, g, b] = legume_tui::map::theme::Theme::PAPER.background;
    assert_eq!(
        img.get_pixel(img.width() - 1, 0).0,
        [r, g, b, 255],
        "the page, not white"
    );
}
