//! A units round's cell types and GO terms on its cell map: placed at their
//! anchors' medians, linked, ringed and starred, each toggled.

use super::*;
use crate::figures::tests::round;
use crate::figures::{CellFigures, CellMap, Colouring, Panel, Size, Style, View};
use legume_numeric::matrix::parquet::{write_named_table, Column};
use legume_tui::map::offsets::Offsets;
use std::path::Path;

const SIZE: Size = Size {
    w: 320,
    h: 240,
    text_px: 12.0,
};

fn strs(v: &[&str]) -> Vec<Box<str>> {
    v.iter().map(|s| Box::from(*s)).collect()
}

/// [`round`] (c0 at (0, 0), c1 (1, 0), c2 (5, 5), c3 (6, 5); A and B) as a
/// units pass leaves it: A anchored on c0, c1 and c2, B on c2 and c3; two
/// GO terms, the first named by the ontology the pass read; K0 aligned
/// with GO:0000001 and K1 with GO:0000002.
pub(crate) fn units_round(dir: &Path, obo: bool) -> crate::manifest::run::Loaded {
    let mut l = round(dir, true);
    let at = |n: &str| dir.join(n).to_string_lossy().into_owned();
    let anchors = [
        ("A", "cell_type", "c0", 0.9),
        ("A", "cell_type", "c1", 0.8),
        ("A", "cell_type", "c2", 0.1),
        ("B", "cell_type", "c2", 0.9),
        ("B", "cell_type", "c3", 0.7),
        ("GO:0000001", "term", "c0", 0.6),
        ("GO:0000001", "term", "c3", 0.5),
        ("GO:0000002", "term", "c2", 0.6),
    ];
    let col = |k: usize| {
        (anchors.iter())
            .map(|a| Box::from([a.0, a.1, a.2][k]))
            .collect::<Vec<Box<str>>>()
    };
    let weight: Vec<f32> = anchors.iter().map(|a| a.3 as f32).collect();
    let rows: Vec<Box<str>> = (0..anchors.len()).map(|i| i.to_string().into()).collect();
    write_named_table(
        &at("r.unit_anchors.parquet"),
        "row",
        &rows,
        &[
            ("unit".into(), Column::Str(&col(0))),
            ("kind".into(), Column::Str(&col(1))),
            ("cell".into(), Column::Str(&col(2))),
            ("weight".into(), Column::F32(&weight)),
        ],
    )
    .unwrap();
    write_named_table(
        &at("r.celltype_units.parquet"),
        "set",
        &strs(&["A", "B"]),
        &[
            ("genes_found".into(), Column::F32(&[3.0, 3.0])),
            ("radius".into(), Column::F32(&[0.5, 0.25])),
        ],
    )
    .unwrap();
    std::fs::write(
        at("r.celltype_links.tsv"),
        "type_a\ttype_b\tlink\nA\tB\t0.4\n",
    )
    .unwrap();
    write_named_table(
        &at("r.term_alignment.parquet"),
        "row",
        &strs(&["0", "1", "2", "3"]),
        &[
            (
                "cluster".into(),
                Column::Str(&strs(&["K0", "K0", "K1", "K1"])),
            ),
            (
                "term".into(),
                Column::Str(&strs(&[
                    "GO:0000001",
                    "GO:0000002",
                    "GO:0000001",
                    "GO:0000002",
                ])),
            ),
            ("score".into(), Column::F32(&[0.5, 0.3, 0.1, 0.4])),
            ("p".into(), Column::F32(&[0.001, 0.2, 0.5, 0.001])),
            ("q".into(), Column::F32(&[0.002, 0.2, 0.5, 0.002])),
        ],
    )
    .unwrap();
    let mut settings = serde_json::json!({ "projection": "units" });
    if obo {
        std::fs::write(
            at("go.obo"),
            "format-version: 1.2\n\n[Term]\nid: GO:0000001\nname: first process\n",
        )
        .unwrap();
        settings["go_obo"] = at("go.obo").into();
    }
    let a = &mut l.manifest.annotate;
    a.unit_anchors = Some("r.unit_anchors.parquet".into());
    a.celltype_units = Some("r.celltype_units.parquet".into());
    a.celltype_links = Some("r.celltype_links.tsv".into());
    a.term_alignment = Some("r.term_alignment.parquet".into());
    a.settings = Some(settings);
    // Written too, for what opens the round from its file.
    std::fs::write(&l.file, serde_json::to_string(&l.manifest).unwrap()).unwrap();
    l
}

pub(crate) fn of_round(l: &crate::manifest::run::Loaded) -> CellFigures {
    let cells = crate::manifest::rounds::read_cells(l).unwrap();
    CellFigures::of_round(l, &cells, None)
}

fn map(f: &CellFigures, style: &Style) -> CellMap {
    f.cell_map(
        Panel::Layout { k: 0 },
        SIZE,
        style,
        View::default(),
        true,
        &Offsets::default(),
    )
    .unwrap()
}

/// The labels `m` placed, by name, with their anchors.
fn placed(m: &CellMap) -> Vec<(&str, (f32, f32))> {
    (m.placed.iter())
        .map(|p| (m.names[p.key as usize].as_ref(), p.anchor))
        .collect()
}

fn layout(f: &CellFigures) -> (&[f32], &[f32]) {
    (&f.layouts[0].x, &f.layouts[0].y)
}

#[test]
fn each_unit_sits_at_the_median_of_its_anchors_on_the_layout() {
    let dir = tempfile::tempdir().unwrap();
    let f = of_round(&units_round(dir.path(), true));
    let u = f.units.as_ref().expect("the round has anchors");
    let (x, y) = layout(&f);
    let m = u.marks(x, y, Some(0));
    let at = |n: &str| m.types.iter().find(|t| &*t.name == n).unwrap().at;
    assert_eq!(at("A"), Some((1.0, 0.0)), "median of c0, c1, c2");
    assert_eq!(at("B"), Some((5.5, 5.0)), "median of c2, c3");
    assert_eq!(m.links, vec![(0, 1, 1.0)], "the one link, the strongest");
    // K0's one significant term, at its anchors' median and named.
    assert_eq!(m.stars.len(), 1);
    let (term, star) = &m.stars[0];
    assert_eq!(
        (&*term.id, term.name.as_str()),
        ("GO:0000001", "first process")
    );
    assert_eq!(*star, (3.0, 2.5));
    assert!((term.score - 0.5).abs() < 1e-6 && (term.q - 0.002).abs() < 1e-6);
    // K1's is the other term; no cluster in focus, no stars.
    assert_eq!(&*u.marks(x, y, Some(1)).stars[0].0.id, "GO:0000002");
    assert!(u.marks(x, y, None).stars.is_empty());
}

#[test]
fn a_ring_is_the_types_radius_read_at_its_anchors_spread() {
    let dir = tempfile::tempdir().unwrap();
    let f = of_round(&units_round(dir.path(), true));
    let (x, y) = layout(&f);
    let m = f.units.as_ref().unwrap().marks(x, y, None);
    let b = m.types.iter().find(|t| &*t.name == "B").unwrap();
    // B's anchors are 0.5 from their median on the layout, and their chord
    // from its centre is √(2 − 2w) on the sphere.
    let chord = |w: f32| (2.0 - 2.0 * w).sqrt();
    let sphere = ((chord(0.9).powi(2) + chord(0.7).powi(2)) / 2.0).sqrt();
    let want = 0.25 * 0.5 / sphere;
    assert!((b.ring.unwrap() - want).abs() < 1e-4, "{:?} {want}", b.ring);
}

#[test]
fn a_clusters_stars_are_its_six_best_significant_terms() {
    let rows: Vec<GoTerm> = (0..9)
        .map(|i| GoTerm {
            id: format!("GO:{i}").into(),
            name: format!("GO:{i}"),
            score: i as f32,
            q: if i == 8 { 0.5 } else { 0.01 },
        })
        .collect();
    let top = top_terms(rows);
    let ids: Vec<&str> = top.iter().map(|t| &*t.id).collect();
    assert_eq!(ids, ["GO:7", "GO:6", "GO:5", "GO:4", "GO:3", "GO:2"]);
}

#[test]
fn without_the_ontology_a_term_goes_by_its_id() {
    let dir = tempfile::tempdir().unwrap();
    let f = of_round(&units_round(dir.path(), false));
    assert_eq!(f.go_terms(Some(0))[0].name, "GO:0000001");
}

#[test]
fn the_map_draws_the_types_at_their_anchors_as_labels_that_move() {
    let dir = tempfile::tempdir().unwrap();
    let f = of_round(&units_round(dir.path(), true));
    let style = Style {
        focus: Some(0),
        ..Style::opening(&f)
    };
    assert_eq!(style.colouring, Colouring::Type);
    let m = map(&f, &style);
    let labels = placed(&m);
    let a: Vec<_> = labels.iter().filter(|(n, _)| *n == "A").collect();
    assert_eq!(
        a.len(),
        1,
        "A once, at its unit, not its cells too: {labels:?}"
    );
    let (dx, _) = m.vp.to_px([1.0, 0.0]);
    assert!(
        (a[0].1 .0 - dx).abs() < 1.0,
        "over the median of its anchors"
    );
    assert!(
        labels.iter().any(|(n, _)| *n == "1"),
        "K0's star is numbered"
    );
    // A label moved by name stays moved.
    let mut moved = Offsets::default();
    moved.nudge("A".to_string(), (0.0, 2.0));
    let m = f
        .cell_map(
            Panel::Layout { k: 0 },
            SIZE,
            &style,
            View::default(),
            true,
            &moved,
        )
        .unwrap();
    assert!(m
        .placed
        .iter()
        .any(|p| &*m.names[p.key as usize] == "A" && p.moved));
}

#[test]
fn g_l_and_r_take_the_stars_links_and_rings_off_and_on() {
    let dir = tempfile::tempdir().unwrap();
    let f = of_round(&units_round(dir.path(), true));
    let on = Style {
        focus: Some(0),
        rings: true,
        ..Style::opening(&f)
    };
    assert!(on.go && on.links, "stars and links open on");
    let pixels = |s: &Style| map(&f, s).frame.to_image().into_raw();
    let mut s = on;
    s.toggle_units('g');
    assert!(!s.go);
    assert!(!placed(&map(&f, &s)).iter().any(|(n, _)| *n == "1"));
    for k in ['l', 'r'] {
        let mut s = on;
        let note = s.toggle_units(k);
        assert_ne!(s, on, "{k}: {note}");
        assert_ne!(pixels(&s), pixels(&on), "{k} changes the map");
        s.toggle_units(k);
        assert_eq!(s, on);
    }
}

#[test]
fn a_round_without_anchors_draws_no_units() {
    let dir = tempfile::tempdir().unwrap();
    let f = of_round(&round(dir.path(), true));
    assert!(f.units.is_none());
    assert!(f.go_terms(Some(0)).is_empty());
    let plain = Style::opening(&f);
    let off = Style {
        go: false,
        links: false,
        rings: false,
        ..plain
    };
    assert_eq!(
        map(&f, &plain).frame.to_image().into_raw(),
        map(&f, &off).frame.to_image().into_raw()
    );
}
