use super::*;
use crate::trajectory::edges::{EdgeRow, Verdict};
use crate::trajectory::figures::{Layout, Style};
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
        style: Style::default(),
        view: View::default(),
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
    for i in 0..5 {
        assert!(p.tile(i, area).unwrap().is_ok(), "tile {i}");
    }
    assert!(p.tile(5, area).is_none());
    let drawn = |p: &FigurePane| {
        p.grid
            .as_ref()
            .unwrap()
            .drawn
            .borrow()
            .iter()
            .filter(|t| t.is_some())
            .count()
    };
    assert_eq!(drawn(&p), 5);
    p.cycle_labels();
    assert_eq!(drawn(&p), 0, "a new label size draws them again");
}

#[test]
fn only_the_scatter_zooms_and_a_new_layout_shows_it_whole() {
    let mut p = pane();
    assert!(p.zoom(true).starts_with("zoom ×1.4"));
    assert!(!p.data.view.is_whole());
    assert_eq!(p.pan(1.0, 0.0), "");
    assert!(p.data.view.cx > 0.5);
    p.next_layout();
    assert!(p.data.view.is_whole(), "m shows the new layout whole");
    p.zoom(true);
    assert_eq!(p.reset_view(), "the whole scatter");
    assert!(p.data.view.is_whole());
    p.next_panel();
    assert_eq!(p.current(), Panel::Order);
    assert!(p.zoom(true).starts_with("only the scatter"));
    assert!(p.data.view.is_whole());
}
