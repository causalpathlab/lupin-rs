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
    let mut p = FigurePane::new(figures(), &Picker::halfblocks());
    p.shown = true;
    p
}

/// [`pane`]'s figures.
fn figures() -> CellFigures {
    let n = 30;
    let x: Vec<f32> = (0..n).map(|i| i as f32).collect();
    let y: Vec<f32> = (0..n).map(|i| (i % 5) as f32).collect();
    let mut d = Mat::zeros(n, 3);
    for i in 0..n {
        d[(i, 1)] = x[i];
        d[(i, 2)] = y[i];
    }
    CellFigures {
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
    }
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
        gallery: super::super::gallery::open(&dir.join(super::super::gallery::DIR)),
        open: true,
        sel: 0,
        statuses: Vec::new(),
        not_listed: Vec::new(),
        delete_armed: None,
        fits: std::cell::Cell::new(true),
    };
    x.gallery
        .add(
            &[listed],
            "a",
            None,
            Trace::of("a", &dir.join("run.senna.json")),
        )
        .unwrap();
    x.refresh();
    assert_eq!((x.gallery.entries().len(), x.not_listed.len()), (1, 1));
    // The not-listed row is the second; listing it puts it first.
    x.sel = 1;
    x.adopt(&dir.join("run.senna.json")).unwrap();
    let on = &x.gallery.entries()[x.sel].path;
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
    let moved = Offsets::default();
    let (img, _) = FigurePane::rasterised(&p.data, look, &p.style, &moved, (cw, ch), area).unwrap();
    assert_eq!((img.width(), img.height()), (30 * cw, 10 * ch));
    let [r, g, b] = legume_tui::map::theme::Theme::PAPER.background;
    assert_eq!(
        img.get_pixel(img.width() - 1, 0).0,
        [r, g, b, 255],
        "the page, not white"
    );
}

use legume_tui::map::offsets::OffsetBook;
use ratatui::crossterm::event::{KeyEvent, KeyModifiers, MouseButton, MouseEvent, MouseEventKind};
use std::path::Path;

/// Where the pointer tests draw the map.
const MAP: Rect = Rect::new(2, 1, 60, 24);

/// [`figures`] as a round's cells, in two clusters, its run `run.L1` in
/// `dir`; labelled by type and drawn in [`MAP`].
fn map_in(dir: &Path) -> FigurePane {
    let mut data = figures();
    data.source = crate::figures::Source::Round;
    data.manifest = dir.join("run.L1.senna.json");
    data.cluster = (0..30).map(|i| i32::from(i >= 15)).collect();
    let mut p = FigurePane::showing(data, &Picker::halfblocks());
    p.style.colouring = crate::figures::Colouring::Type;
    drawn(&p);
    p
}

fn drawn(p: &FigurePane) {
    assert!(p.protocol(MAP).is_ok());
}

fn at(kind: MouseEventKind, (column, row): (u16, u16)) -> MouseEvent {
    MouseEvent {
        kind,
        column,
        row,
        modifiers: KeyModifiers::NONE,
    }
}

fn key(code: KeyCode, modifiers: KeyModifiers) -> KeyEvent {
    KeyEvent::new(code, modifiers)
}

/// The screen cell over frame pixel `(x, y)` of the map.
fn cell_of(p: &FigurePane, (x, y): (f32, f32)) -> (u16, u16) {
    let (cw, ch) = p.cell_px();
    (
        MAP.x + (x / cw as f32) as u16,
        MAP.y + (y / ch as f32) as u16,
    )
}

/// Where label `name` was placed on the map as last drawn, and whether it
/// was moved.
fn label(p: &FigurePane, name: &str) -> ([f32; 4], bool) {
    let on = p.on_screen.borrow().clone().expect("drawn");
    let k = on
        .names
        .iter()
        .position(|n| n.as_ref() == name)
        .expect("named");
    let q = on
        .placed
        .iter()
        .find(|q| q.key as usize == k)
        .expect("placed");
    (q.rect, q.moved)
}

/// The screen cell over label `name`.
fn label_cell(p: &FigurePane, name: &str) -> (u16, u16) {
    let (r, _) = label(p, name);
    cell_of(p, (0.5 * (r[0] + r[2]), 0.5 * (r[1] + r[3])))
}

/// The screen cell over cell `i` of the first layout.
fn cell_cell(p: &FigurePane, i: usize) -> (u16, u16) {
    let on = p.on_screen.borrow().clone().expect("drawn");
    let l = &p.data.layouts[0];
    cell_of(p, on.vp.to_px([l.x[i], l.y[i]]))
}

#[test]
fn the_pointer_shows_and_picks_a_cluster_drags_the_map_and_the_wheel_zooms() {
    use MouseButton::Left;
    use MouseEventKind::{Down, Drag, Moved, ScrollUp, Up};
    let dir = tempfile::tempdir().unwrap();
    let mut p = map_in(dir.path());
    p.style.labels = None;
    drawn(&p);
    let c = cell_cell(&p, 20);
    assert_eq!(p.mouse(&at(Moved, c)), MapMouse::Hover(Some(1)));
    assert_eq!(p.mouse(&at(Down(Left), c)), MapMouse::None);
    assert_eq!(p.mouse(&at(Up(Left), c)), MapMouse::Click(Some(1)));
    let off = (MAP.x, MAP.y + MAP.height - 1);
    assert_eq!(
        p.mouse(&at(Moved, off)),
        MapMouse::Hover(None),
        "no cell near"
    );
    // Whole, the map does not pan; the wheel zooms it, then it does.
    let mid = (MAP.x + 30, MAP.y + 12);
    assert_eq!(p.mouse(&at(ScrollUp, mid)), MapMouse::Redraw);
    assert_eq!(p.view.level, 1);
    drawn(&p);
    let cx = p.view.cx;
    p.mouse(&at(Down(Left), mid));
    assert_eq!(
        p.mouse(&at(Drag(Left), (mid.0 + 2, mid.1))),
        MapMouse::Redraw
    );
    p.mouse(&at(Up(Left), (mid.0 + 2, mid.1)));
    assert!(p.view.cx < cx, "dragged right, the view goes left");
    // Off the map, nothing.
    assert_eq!(p.mouse(&at(ScrollUp, (0, 0))), MapMouse::None);
}

#[test]
fn a_dragged_label_is_kept_beside_the_run_and_read_again() {
    use MouseButton::Left;
    use MouseEventKind::{Down, Drag, Up};
    let dir = tempfile::tempdir().unwrap();
    let mut p = map_in(dir.path());
    assert!(!label(&p, "CT1").1);
    let from = label_cell(&p, "CT1");
    let to = (from.0, from.1 + 2);
    p.mouse(&at(Down(Left), from));
    assert_eq!(p.mouse(&at(Drag(Left), to)), MapMouse::Redraw);
    p.mouse(&at(Up(Left), to));
    assert_eq!(p.view, View::default(), "the label moved, not the map");
    let file = dir.path().join("run.L1.lupin_labels.json");
    let o = OffsetBook::load(&file)
        .unwrap()
        .layout("run.L1:layout_phate");
    let (dx, dy) = o.get(&"CT1".to_string()).expect("kept");
    assert!(dx == 0.0 && dy > 0.0, "{dx} {dy}");
    let written: Vec<_> = std::fs::read_dir(dir.path())
        .unwrap()
        .map(|e| e.unwrap().file_name())
        .collect();
    assert_eq!(written, ["run.L1.lupin_labels.json"], "nothing else");
    // Opened again, the label is where it was put.
    let again = map_in(dir.path());
    assert!(label(&again, "CT1").1);
    assert!(!label(&again, "CT2").1);
}

#[test]
fn a_label_let_go_without_a_release_is_kept_too() {
    use MouseButton::Left;
    use MouseEventKind::{Down, Drag};
    let dir = tempfile::tempdir().unwrap();
    let mut p = map_in(dir.path());
    let from = label_cell(&p, "CT1");
    p.mouse(&at(Down(Left), from));
    assert_eq!(
        p.mouse(&at(Drag(Left), (from.0, from.1 + 2))),
        MapMouse::Redraw
    );
    assert_eq!(p.let_go(), None);
    let file = dir.path().join("run.L1.lupin_labels.json");
    let o = OffsetBook::load(&file)
        .unwrap()
        .layout("run.L1:layout_phate");
    assert!(o.get(&"CT1".to_string()).is_some(), "kept");
}

#[test]
fn a_labels_file_that_cannot_be_read_is_left_as_it_is_and_said() {
    use MouseButton::Left;
    use MouseEventKind::{Down, Drag, Up};
    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("run.L1.lupin_labels.json");
    std::fs::write(&file, "{ not a book").unwrap();
    let mut p = map_in(dir.path());
    let from = label_cell(&p, "CT1");
    let to = (from.0, from.1 + 2);
    p.mouse(&at(Down(Left), from));
    match p.mouse(&at(Drag(Left), to)) {
        MapMouse::Say(s) => assert!(s.contains("run.L1.lupin_labels.json"), "{s}"),
        r => panic!("{r:?}"),
    }
    p.mouse(&at(Up(Left), to));
    match p.edit_key(&key(KeyCode::Char('E'), KeyModifiers::SHIFT)) {
        FigureKey::Say(s) => assert!(s.contains("run.L1.lupin_labels.json"), "{s}"),
        r => panic!("{r:?}"),
    }
    assert!(!p.editing());
    drawn(&p);
    assert!(!label(&p, "CT1").1, "the labels stay where they are");
    assert_eq!(std::fs::read_to_string(&file).unwrap(), "{ not a book");
}

#[test]
fn capital_e_moves_the_labels_from_the_keys_and_keeps_them() {
    let dir = tempfile::tempdir().unwrap();
    let mut p = map_in(dir.path());
    let none = KeyModifiers::NONE;
    assert_eq!(
        p.edit_key(&key(KeyCode::Char('t'), none)),
        FigureKey::Ignored
    );
    assert!(matches!(
        p.edit_key(&key(KeyCode::Char('E'), KeyModifiers::SHIFT)),
        FigureKey::Say(_)
    ));
    assert!(p.editing());
    assert!(matches!(
        p.edit_key(&key(KeyCode::Right, none)),
        FigureKey::Say(_)
    ));
    let file = dir.path().join("run.L1.lupin_labels.json");
    let o = OffsetBook::load(&file)
        .unwrap()
        .layout("run.L1:layout_phate");
    let moved: Vec<_> = o.iter().map(|(k, d)| (k.clone(), d)).collect();
    assert_eq!(
        moved,
        [("CT1".to_string(), (0.5, 0.0))],
        "the first as read"
    );
    // Tab goes on to the next; Backspace puts it back.
    p.edit_key(&key(KeyCode::Tab, none));
    p.edit_key(&key(KeyCode::Down, KeyModifiers::SHIFT));
    p.edit_key(&key(KeyCode::BackTab, none));
    p.edit_key(&key(KeyCode::Backspace, none));
    let o = OffsetBook::load(&file)
        .unwrap()
        .layout("run.L1:layout_phate");
    let moved: Vec<_> = o.iter().map(|(k, d)| (k.clone(), d)).collect();
    assert_eq!(moved, [("CT2".to_string(), (0.0, 2.0))]);
    assert_eq!(p.edit_key(&key(KeyCode::Esc, none)), FigureKey::Took);
    assert!(!p.editing());
    // Out of the mode, the arrows pan again.
    assert_eq!(p.edit_key(&key(KeyCode::Right, none)), FigureKey::Ignored);
}

/// The screen cell over data point `xy` of `p`'s scatter, as last drawn
/// wherever it was drawn.
pub(crate) fn screen_at(p: &FigurePane, xy: [f32; 2]) -> (u16, u16) {
    let on = p.on_screen.borrow().clone().expect("drawn");
    let area = p.area.get().expect("drawn");
    let (cw, ch) = p.cell_px();
    let (x, y) = on.vp.to_px(xy);
    (
        area.x + (x / cw as f32) as u16,
        area.y + (y / ch as f32) as u16,
    )
}

/// [`map_in`], its exports logged in `dir`'s gallery.
fn saving_in(dir: &Path) -> FigurePane {
    let mut p = map_in(dir);
    p.exports.gallery = super::super::gallery::open(&dir.join(super::super::gallery::DIR));
    p
}

#[test]
fn p_on_the_map_asks_how_to_save_it_and_esc_saves_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let dir = dir.path().canonicalize().unwrap();
    let mut p = saving_in(&dir);
    assert!(!p.saving());
    assert_ne!(p.key(KeyCode::Char('p')), FigureKey::Ignored);
    assert!(p.saving(), "the save dialog is open");
    assert_eq!(
        p.save_key(&key(KeyCode::Esc, KeyModifiers::NONE)),
        FigureKey::Took
    );
    assert!(!p.saving());
    assert!(std::fs::read_dir(&dir).unwrap().all(|f| f
        .unwrap()
        .path()
        .extension()
        .is_none_or(|e| e != "pdf")));
}

#[test]
fn a_saved_map_is_a_pdf_whose_text_is_text_logged_with_a_thumbnail() {
    let dir = tempfile::tempdir().unwrap();
    let dir = dir.path().canonicalize().unwrap();
    let mut p = saving_in(&dir);
    p.key(KeyCode::Char('p'));
    let FigureKey::Say(said) = p.save_key(&key(KeyCode::Enter, KeyModifiers::NONE)) else {
        panic!("a save says where it went");
    };
    assert!(!p.saving(), "{said}");
    let pdf = dir.join("run.L1.cells.layout_phate.pdf");
    assert!(said.contains(&pdf.display().to_string()), "{said}");
    let bytes = std::fs::read(&pdf).unwrap();
    let s = String::from_utf8_lossy(&bytes);
    assert!(s.starts_with("%PDF-"));
    // The label CT1 set as text over the points' image.
    assert!(s.contains("(CT1) Tj"), "the labels are text: {}", tj(&s));
    assert!(s.contains("/Subtype /Image"));
    // 7 in wide, as the dialog opens.
    assert!(s.contains("/MediaBox [0 0 504 "), "7 in is 504 pt");
    assert!(
        !pdf.with_extension("svg").exists(),
        "a PDF only, as senna's"
    );

    let g = &p.exports.gallery;
    let e = &g.entries()[0];
    assert_eq!(e.path, pdf);
    assert_eq!(e.extra.panel, "layout_phate");
    assert_eq!(e.extra.manifest, dir.join("run.L1.senna.json"));
    assert_eq!(e.files.len(), 1, "the file traced");
    let thumb = image::open(g.thumb(e).expect("a thumbnail")).unwrap();
    assert!((1..=240).contains(&thumb.width()));
    assert_eq!(p.exports.statuses, vec![Status::Intact]);

    // The next save is offered a free name.
    p.key(KeyCode::Char('p'));
    p.save_key(&key(KeyCode::Enter, KeyModifiers::NONE));
    assert!(dir.join("run.L1.cells.layout_phate-2.pdf").is_file());
    assert_eq!(p.exports.gallery.entries().len(), 2);
}

#[test]
fn the_other_panels_still_export_at_once_as_svg_and_pdf_with_a_thumbnail() {
    let dir = tempfile::tempdir().unwrap();
    let dir = dir.path().canonicalize().unwrap();
    let mut data = figures();
    data.manifest = dir.join("run.senna.json");
    let mut p = FigurePane::showing(data, &Picker::halfblocks());
    p.exports.gallery = super::super::gallery::open(&dir.join(super::super::gallery::DIR));
    while p.current() != Panel::Order {
        p.next_panel();
    }
    assert!(matches!(p.key(KeyCode::Char('p')), FigureKey::Say(_)));
    assert!(!p.saving(), "no dialog");
    let file = |ext: &str| dir.join(format!("run.trajectory.order.{ext}"));
    assert!(file("svg").is_file() && file("pdf").is_file());
    let g = &p.exports.gallery;
    assert!(g.thumb(&g.entries()[0]).is_some_and(|t| t.is_file()));
}

#[test]
fn a_rounds_cell_map_is_saved_without_the_tui() {
    let dir = tempfile::tempdir().unwrap();
    let loaded = crate::figures::tests::round(dir.path(), true);
    let out = dir.path().join("map.pdf");
    let pdf = legume_tui::export::PdfArgs {
        pdf: Some(out.clone()),
        width: 3.5,
        dpi: 150,
    };
    let before = std::fs::read(&loaded.file).unwrap();
    assert_eq!(
        round_pdf(&loaded.file.to_string_lossy(), &pdf).unwrap(),
        out
    );
    let bytes = std::fs::read(&out).unwrap();
    let s = String::from_utf8_lossy(&bytes);
    assert!(s.starts_with("%PDF-"));
    // 3.5 in at 150 dpi: 525 px, 252 pt.
    assert!(s.contains("/MediaBox [0 0 252 "), "3.5 in wide");
    // Its labels, the round's types, are text.
    assert!(s.contains("(A) Tj"), "the labels are text: {}", tj(&s));
    assert_eq!(
        std::fs::read(&loaded.file).unwrap(),
        before,
        "the manifest is left"
    );
}

#[test]
fn a_rounds_cell_map_never_replaces_its_manifest_or_another_file() {
    let dir = tempfile::tempdir().unwrap();
    let loaded = crate::figures::tests::round(dir.path(), true);
    let from = loaded.file.to_string_lossy().to_string();
    let to = |out: PathBuf| legume_tui::export::PdfArgs {
        pdf: Some(out),
        width: 3.5,
        dpi: 150,
    };
    let before = std::fs::read(&loaded.file).unwrap();
    // The manifest itself, and through a path that only resolves to it.
    let roundabout = dir.path().join(".").join(loaded.file.file_name().unwrap());
    for at in [loaded.file.clone(), roundabout] {
        let e = round_pdf(&from, &to(at.clone())).unwrap_err();
        assert!(e.to_string().contains("manifest"), "{e}");
    }
    assert_eq!(
        std::fs::read(&loaded.file).unwrap(),
        before,
        "the manifest is left"
    );
    // Not a .pdf.
    let png = dir.path().join("map.png");
    assert!(round_pdf(&from, &to(png.clone())).is_err());
    assert!(!png.exists());
    // A PDF already there: kept without --overwrite (no terminal to ask on).
    let taken = dir.path().join("taken.pdf");
    std::fs::write(&taken, b"mine").unwrap();
    let e = round_pdf(&from, &to(taken.clone())).unwrap_err();
    assert!(e.to_string().contains("--overwrite"), "{e}");
    assert_eq!(std::fs::read(&taken).unwrap(), b"mine");
}

/// The text a PDF sets, as its `Tj` lines.
fn tj(pdf: &str) -> String {
    pdf.lines()
        .filter(|l| l.ends_with("Tj"))
        .collect::<Vec<_>>()
        .join(" | ")
}

/// A round of four cells in K0 and K1 with a posterior per cell, its
/// exports logged in `dir`'s gallery.
fn charts_in(dir: &Path) -> FigurePane {
    use crate::figures::charts::tests::{charted, posterior};
    posterior(dir);
    let data = charted(dir, |m| {
        m["annotate"]["annotation"] = "r.annotation.parquet".into();
    });
    let mut p = FigurePane::showing(data, &Picker::halfblocks());
    p.exports.gallery = super::super::gallery::open(&dir.join(super::super::gallery::DIR));
    drawn(&p);
    p
}

#[test]
fn h_steps_from_the_map_to_the_structure_plot_and_the_heatmap_and_esc_returns() {
    use crate::figures::charts::{HeatKind, Mixture};
    let dir = tempfile::tempdir().unwrap();
    let mut p = charts_in(dir.path());
    let map = p.current();
    assert!(map.is_scatter());
    let FigureKey::Say(said) = p.key(KeyCode::Char('H')) else {
        panic!("H says what it shows");
    };
    assert_eq!(p.current(), Panel::Structure { of: Mixture::Types });
    assert!(said.contains("esc"), "{said}");
    drawn(&p);
    // The heatmap next; this round records no expression, which it says.
    let FigureKey::Say(said) = p.key(KeyCode::Char('H')) else {
        panic!("H says what it shows");
    };
    assert_eq!(
        p.current(),
        Panel::Heatmap {
            of: HeatKind::Markers
        }
    );
    assert!(said.contains("cluster_expression"), "{said}");
    drawn(&p);
    // Then the map again.
    p.key(KeyCode::Char('H'));
    assert_eq!(p.current(), map);
    // Esc leaves a chart for the map, then is not the pane's.
    p.key(KeyCode::Char('H'));
    assert!(p.current().is_chart());
    assert_ne!(p.key(KeyCode::Esc), FigureKey::Ignored);
    assert_eq!(p.current(), map);
    assert_eq!(p.key(KeyCode::Esc), FigureKey::Ignored);
}

#[test]
fn h_while_moving_labels_ends_the_moving_so_esc_on_the_chart_is_the_maps() {
    let dir = tempfile::tempdir().unwrap();
    let mut p = charts_in(dir.path());
    let map = p.current();
    let shift = KeyModifiers::SHIFT;
    assert!(matches!(
        p.edit_key(&key(KeyCode::Char('E'), shift)),
        FigureKey::Say(_)
    ));
    assert!(p.editing());
    // H is not the labels': it passes to the pane, which opens a chart.
    assert_eq!(
        p.edit_key(&key(KeyCode::Char('H'), shift)),
        FigureKey::Ignored
    );
    p.key(KeyCode::Char('H'));
    assert!(p.current().is_chart());
    assert!(!p.editing(), "the labels move on the map, not a chart");
    // So Esc goes to the pane and returns to the map at once.
    let esc = key(KeyCode::Esc, KeyModifiers::NONE);
    assert_eq!(p.edit_key(&esc), FigureKey::Ignored);
    p.key(KeyCode::Esc);
    assert_eq!(p.current(), map);
}

#[test]
fn c_on_a_chart_steps_what_it_shows() {
    use crate::figures::charts::{HeatKind, Mixture};
    let dir = tempfile::tempdir().unwrap();
    let mut p = charts_in(dir.path());
    p.key(KeyCode::Char('H'));
    // Not a topic run: the cell types are all a structure plot stacks.
    let FigureKey::Say(said) = p.key(KeyCode::Char('c')) else {
        panic!("c says why");
    };
    assert_eq!(p.current(), Panel::Structure { of: Mixture::Types });
    assert!(said.contains("topics"), "{said}");
    p.key(KeyCode::Char('H'));
    let mut seen = vec![p.current()];
    for _ in 0..3 {
        p.key(KeyCode::Char('c'));
        seen.push(p.current());
    }
    let heat = |of| Panel::Heatmap { of };
    assert_eq!(
        seen,
        [
            heat(HeatKind::Markers),
            heat(HeatKind::Terms),
            heat(HeatKind::Types),
            heat(HeatKind::Markers)
        ]
    );
    // The map, H later, opens on the heatmap chosen.
    p.key(KeyCode::Char('c'));
    p.key(KeyCode::Esc);
    p.key(KeyCode::Char('H'));
    p.key(KeyCode::Char('H'));
    assert_eq!(p.current(), heat(HeatKind::Terms));
}

#[test]
fn a_saved_chart_is_a_pdf_whose_text_is_text_logged_with_a_thumbnail() {
    let dir = tempfile::tempdir().unwrap();
    let dir = dir.path().canonicalize().unwrap();
    let mut p = charts_in(&dir);
    p.key(KeyCode::Char('H'));
    assert_ne!(p.key(KeyCode::Char('p')), FigureKey::Ignored);
    assert!(p.saving(), "the save dialog is open");
    let FigureKey::Say(said) = p.save_key(&key(KeyCode::Enter, KeyModifiers::NONE)) else {
        panic!("a save says where it went");
    };
    let pdf = dir.join("r.cells.structure_types.pdf");
    assert!(said.contains(&pdf.display().to_string()), "{said}");
    let s = String::from_utf8_lossy(&std::fs::read(&pdf).unwrap()).into_owned();
    assert!(s.starts_with("%PDF-"));
    assert!(
        s.contains("(K0) Tj"),
        "the panels' names are text: {}",
        tj(&s)
    );
    let g = &p.exports.gallery;
    let e = &g.entries()[0];
    assert_eq!(e.path, pdf);
    assert_eq!(e.extra.panel, "structure_types");
    assert!(g.thumb(e).is_some_and(|t| t.is_file()));
}

#[test]
fn a_trajectorys_figures_have_no_charts() {
    let mut p = pane();
    let before = p.current();
    let FigureKey::Say(said) = p.key(KeyCode::Char('H')) else {
        panic!("H says why not");
    };
    assert!(said.contains("round"), "{said}");
    assert_eq!(p.current(), before);
}
