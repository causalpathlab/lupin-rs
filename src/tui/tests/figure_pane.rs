use super::*;
use crate::trajectory::edges::{EdgeRow, Verdict};
use crate::trajectory::figures::Layout;
use legume_numeric::matrix::dense_mat_io::Mat;
use ratatui::crossterm::event::KeyCode;

/// Two layouts, a diffusion map and one edge, on placeholder types.
fn pane() -> FigurePane {
    let n = 30;
    let x: Vec<f32> = (0..n).map(|i| i as f32).collect();
    let y: Vec<f32> = (0..n).map(|i| (i % 5) as f32).collect();
    let mut d = Mat::zeros(n, 3);
    for i in 0..n {
        d[(i, 1)] = x[i];
        d[(i, 2)] = y[i];
    }
    let data = TrajectoryData {
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
        ..TrajectoryData::default()
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
