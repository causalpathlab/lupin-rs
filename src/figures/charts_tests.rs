//! A round's charts from the tables its manifest records: the structure
//! plot's per-cell mixtures, the heatmaps' features × clusters, and what is
//! said when a table is not there.

use super::*;
use crate::figures::tests::round;
use crate::figures::{CellFigures, Panel, Picture, Size, Style, View};
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::parquet::{write_named_table, Column};
use legume_numeric::matrix::traits::IoOps;
use legume_tui::map::offsets::Offsets;
use std::path::Path;

const SIZE: Size = Size {
    w: 480,
    h: 320,
    text_px: 12.0,
};

fn strs(v: &[&str]) -> Vec<Box<str>> {
    v.iter().map(|s| Box::from(*s)).collect()
}

/// `values` (row-major) as a table at `dir/name`, rows named `rows` (under
/// `head`) and columns `cols`.
fn table(dir: &Path, name: &str, head: &str, rows: &[&str], cols: &[&str], values: &[f32]) {
    Mat::from_row_slice(rows.len(), cols.len(), values)
        .to_parquet_with_names(
            &dir.join(name).to_string_lossy(),
            (Some(&strs(rows)), Some(head)),
            Some(&strs(cols)),
        )
        .unwrap();
}

/// [`round`] (c0 and c1 in K0, c2 and c3 in K1) with `edit` made to its
/// manifest, as the TUI reads it.
pub(crate) fn charted(dir: &Path, edit: impl FnOnce(&mut serde_json::Value)) -> CellFigures {
    let l = round(dir, true);
    let mut m: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&l.file).unwrap()).unwrap();
    edit(&mut m);
    std::fs::write(&l.file, m.to_string()).unwrap();
    let l = crate::manifest::run::load(&l.file.to_string_lossy()).unwrap();
    let cells = crate::manifest::rounds::read_cells(&l).unwrap();
    CellFigures::of_round(&l, &cells, None)
}

const CELLS: [&str; 4] = ["c0", "c1", "c2", "c3"];

/// Each cell's posterior over CT1 and CT2, as enrichment writes it.
pub(crate) fn posterior(dir: &Path) {
    let v = [0.9, 0.1, 0.8, 0.2, 0.3, 0.7, 0.1, 0.9];
    table(
        dir,
        "r.annotation.parquet",
        "cell",
        &CELLS,
        &["CT1", "CT2"],
        &v,
    );
}

/// Each cell's membership of the CT1 and CT2 regions, as a units pass
/// writes it.
fn membership(dir: &Path) {
    let v = [0.6, 0.4, 0.4, 0.6, 0.2, 0.8, 0.0, 1.0];
    let name = "r.celltype_membership.parquet";
    table(dir, name, "cell", &CELLS, &["CT1", "CT2"], &v);
}

fn types_of(f: &CellFigures) -> Arc<Mix> {
    f.mixture(Mixture::Types).expect("a mixture")
}

fn picture(f: &CellFigures, panel: Panel) -> image::RgbaImage {
    let fig = f
        .figure(
            panel,
            SIZE,
            &Style::default(),
            View::default(),
            &Offsets::default(),
        )
        .unwrap();
    assert!(matches!(fig, Picture::Chart(_)), "a chart is drawn as one");
    crate::figures::render(&fig).unwrap()
}

/// Pixels of `img` in a colour, not a grey: bars and heatmap cells, never
/// the page or its text.
fn coloured(img: &image::RgbaImage) -> usize {
    img.pixels()
        .filter(|p| {
            let [r, g, b, _] = p.0;
            r.max(g).max(b) - r.min(g).min(b) > 40
        })
        .count()
}

/// The lines of text a chart sets, kept as text.
fn lines(f: &CellFigures, panel: Panel) -> Vec<String> {
    let mut frame = f.chart_frame(panel, SIZE, true);
    frame.take_text().into_iter().map(|r| r.text).collect()
}

/// What a chart says in its text.
fn said(f: &CellFigures, panel: Panel) -> String {
    lines(f, panel).join(" | ")
}

/// `z` is `want`, to rounding, missing where `want` is.
fn near(z: &[f32], want: &[f32]) {
    assert_eq!(z.len(), want.len(), "{z:?}");
    for (a, b) in z.iter().zip(want) {
        assert!(
            (a.is_nan() && b.is_nan()) || (a - b).abs() < 1e-4,
            "{z:?} against {want:?}"
        );
    }
}

#[test]
fn the_structure_plot_stacks_each_cells_posterior_in_a_panel_per_cluster() {
    let dir = tempfile::tempdir().unwrap();
    posterior(dir.path());
    membership(dir.path());
    let f = charted(dir.path(), |m| {
        m["annotate"]["annotation"] = "r.annotation.parquet".into();
        m["annotate"]["celltype_membership"] = "r.celltype_membership.parquet".into();
    });
    // The posterior, when the pass wrote one, over the units' regions.
    let mix = types_of(&f);
    assert_eq!(mix.names, ["CT1", "CT2"]);
    assert_eq!(&mix.w[..4], &[0.9, 0.1, 0.8, 0.2]);
    let l = f.structure(&mix);
    let panels: Vec<(&str, &[u32])> = (l.panels.iter())
        .map(|p| (p.name.as_str(), &p.rows[..]))
        .collect();
    assert_eq!(panels, [("K0", &[0, 1][..]), ("K1", &[3, 2][..])]);
    let img = picture(&f, Panel::Structure { of: Mixture::Types });
    assert_eq!((img.width(), img.height()), (SIZE.w, SIZE.h));
    assert!(coloured(&img) > 100, "bars in the types' colours");
    let text = said(&f, Panel::Structure { of: Mixture::Types });
    for name in ["K0", "K1", "CT1", "CT2"] {
        assert!(text.contains(name), "{name} in {text}");
    }
}

#[test]
fn a_units_round_stacks_each_cells_membership_of_the_regions() {
    let dir = tempfile::tempdir().unwrap();
    membership(dir.path());
    let f = charted(dir.path(), |m| {
        m["annotate"]["celltype_membership"] = "r.celltype_membership.parquet".into();
    });
    let mix = types_of(&f);
    assert_eq!(mix.names, ["CT1", "CT2"]);
    assert_eq!(&mix.w[..2], &[0.6, 0.4]);
    assert!(f.has_mixture(Mixture::Types));
    assert!(!f.has_mixture(Mixture::Topics), "not a topic run");
}

#[test]
fn a_topic_run_stacks_its_cells_topic_shares() {
    let dir = tempfile::tempdir().unwrap();
    // log θ, as a topic run records its latent.
    let theta = [0.75f32, 0.25, 0.5, 0.5, 0.2, 0.8, 0.1, 0.9].map(f32::ln);
    table(
        dir.path(),
        "r.latent.parquet",
        "cell",
        &CELLS,
        &["T0", "T1"],
        &theta,
    );
    let f = charted(dir.path(), |m| {
        m["kind"] = "topic".into();
        m["outputs"]["latent"] = "r.latent.parquet".into();
    });
    assert!(f.has_mixture(Mixture::Topics));
    let mix = f.mixture(Mixture::Topics).unwrap();
    assert_eq!(mix.names, ["T0", "T1"]);
    for (got, want) in mix.w.iter().zip([0.75, 0.25, 0.5, 0.5, 0.2, 0.8, 0.1, 0.9]) {
        assert!((got - want).abs() < 1e-5, "{:?}", mix.w);
    }
    assert!(
        coloured(&picture(
            &f,
            Panel::Structure {
                of: Mixture::Topics
            }
        )) > 100
    );
}

#[test]
fn the_markers_heatmap_has_the_panels_genes_by_type_z_scored_over_the_clusters() {
    let dir = tempfile::tempdir().unwrap();
    let genes = ["GENE1", "GENE2", "GENE3", "GENE4"];
    let v = [5.0, 1.0, 2.0, 2.0, 0.0, 4.0, 3.0, 3.0];
    table(
        dir.path(),
        "r.cluster_expression.parquet",
        "gene",
        &genes,
        &["K0", "K1"],
        &v,
    );
    std::fs::write(
        dir.path().join("r.markers.tsv"),
        "gene\tcell_type\nGENE3\tCT2\nGENE1\tCT1\nGENE2\tCT1\nGENE9\tCT2\n",
    )
    .unwrap();
    let f = charted(dir.path(), |m| {
        m["annotate"]["cluster_expression"] = "r.cluster_expression.parquet".into();
        m["annotate"]["markers"] = "r.markers.tsv".into();
    });
    let h = f.heat(HeatKind::Markers).unwrap();
    assert_eq!(h.cols, ["K0", "K1"]);
    // Grouped by type; a marker the table lacks is left out.
    assert_eq!(h.rows, ["GENE1 · CT1", "GENE2 · CT1", "GENE3 · CT2"]);
    near(&h.z, &[1.0, -1.0, 0.0, 0.0, -1.0, 1.0]);
    let img = picture(
        &f,
        Panel::Heatmap {
            of: HeatKind::Markers,
        },
    );
    assert!(coloured(&img) > 100, "cells in the diverging colours");
    assert!(said(
        &f,
        Panel::Heatmap {
            of: HeatKind::Markers
        }
    )
    .contains("GENE1"));
}

#[test]
fn the_terms_heatmap_shows_the_starred_terms_scores_missing_where_q_is_not_below_005() {
    let dir = tempfile::tempdir().unwrap();
    let rows = [
        ("K0", "GO:0000001", 0.8, 0.01),
        ("K0", "GO:0000002", 0.1, 0.5),
        ("K0", "GO:0000003", 0.3, 0.9),
        ("K1", "GO:0000001", 0.2, 0.05),
        ("K1", "GO:0000002", 0.7, 0.02),
        ("K1", "GO:0000003", 0.1, 0.9),
    ];
    let col = |k: usize| {
        (rows.iter())
            .map(|r| Box::from([r.0, r.1][k]))
            .collect::<Vec<Box<str>>>()
    };
    let score: Vec<f32> = rows.iter().map(|r| r.2).collect();
    let q: Vec<f32> = rows.iter().map(|r| r.3).collect();
    let names: Vec<Box<str>> = (0..rows.len()).map(|i| i.to_string().into()).collect();
    write_named_table(
        &dir.path()
            .join("r.term_alignment.parquet")
            .to_string_lossy(),
        "row",
        &names,
        &[
            ("cluster".into(), Column::Str(&col(0))),
            ("term".into(), Column::Str(&col(1))),
            ("score".into(), Column::F32(&score)),
            ("p".into(), Column::F32(&q)),
            ("q".into(), Column::F32(&q)),
        ],
    )
    .unwrap();
    let f = charted(dir.path(), |m| {
        m["annotate"]["term_alignment"] = "r.term_alignment.parquet".into();
    });
    let h = f.heat(HeatKind::Terms).unwrap();
    assert_eq!(h.cols, ["K0", "K1"]);
    // Only the terms a cluster stars; GO:0000003 is below no cluster's q.
    assert_eq!(h.rows, ["GO:0000001", "GO:0000002"]);
    // z over every cluster, then missing where q = 0.05 is not below it.
    near(&h.z, &[1.0, f32::NAN, f32::NAN, 1.0]);
    assert!(
        coloured(&picture(
            &f,
            Panel::Heatmap {
                of: HeatKind::Terms
            }
        )) > 100
    );
}

#[test]
fn the_types_heatmap_reads_the_enrichments_q_else_the_units_membership_share() {
    let dir = tempfile::tempdir().unwrap();
    let q = [0.9, 0.1, 0.2, 0.8];
    let name = "r.cluster_celltype_q.parquet";
    table(
        dir.path(),
        name,
        "cluster",
        &["K0", "K1"],
        &["CT1", "CT2"],
        &q,
    );
    membership(dir.path());
    let f = charted(dir.path(), |m| {
        m["annotate"]["cluster_celltype_q"] = "r.cluster_celltype_q.parquet".into();
        m["annotate"]["celltype_membership"] = "r.celltype_membership.parquet".into();
    });
    let h = f.heat(HeatKind::Types).unwrap();
    assert_eq!(h.rows, ["CT1", "CT2"]);
    assert_eq!(h.cols, ["K0", "K1"]);
    near(&h.z, &[1.0, -1.0, -1.0, 1.0]);
    assert!(h.title.contains('Q'), "{}", h.title);

    // A units round: each cluster's mean membership, CT1 0.5 in K0 and
    // 0.1 in K1.
    let dir = tempfile::tempdir().unwrap();
    membership(dir.path());
    let f = charted(dir.path(), |m| {
        m["annotate"]["celltype_membership"] = "r.celltype_membership.parquet".into();
    });
    let h = f.heat(HeatKind::Types).unwrap();
    assert_eq!(h.rows, ["CT1", "CT2"]);
    near(&h.z, &[1.0, -1.0, -1.0, 1.0]);
    assert!(h.title.contains("membership"), "{}", h.title);
    assert!(
        coloured(&picture(
            &f,
            Panel::Heatmap {
                of: HeatKind::Types
            }
        )) > 100
    );
}

#[test]
fn a_chart_without_its_table_says_so_instead_of_drawing() {
    let dir = tempfile::tempdir().unwrap();
    let f = charted(dir.path(), |_| {});
    let cases = [
        (
            Panel::Structure { of: Mixture::Types },
            "annotate.annotation",
        ),
        (
            Panel::Structure {
                of: Mixture::Topics,
            },
            "topic",
        ),
        (
            Panel::Heatmap {
                of: HeatKind::Markers,
            },
            "cluster_expression",
        ),
        (
            Panel::Heatmap {
                of: HeatKind::Terms,
            },
            "term_alignment",
        ),
        (
            Panel::Heatmap {
                of: HeatKind::Types,
            },
            "cluster_celltype_q",
        ),
    ];
    for (panel, why) in cases {
        let e = f.chart_missing(panel).expect("missing");
        assert!(e.contains(why), "{panel:?}: {e}");
        // One line on the page saying so, and nothing drawn.
        let text = lines(&f, panel);
        assert_eq!(text.len(), 1, "{panel:?}: {text:?}");
        assert!(text[0].contains(why), "{panel:?} says so: {text:?}");
        assert_eq!(coloured(&picture(&f, panel)), 0, "{panel:?} draws nothing");
    }
    // A table recorded but gone says that.
    let f = charted(dir.path(), |m| {
        m["annotate"]["cluster_expression"] = "r.gone.parquet".into();
    });
    let e = f
        .chart_missing(Panel::Heatmap {
            of: HeatKind::Markers,
        })
        .unwrap();
    assert!(e.contains("gone"), "{e}");
}

#[test]
fn a_chart_is_written_as_a_pdf_whose_text_is_text() {
    let dir = tempfile::tempdir().unwrap();
    posterior(dir.path());
    let f = charted(dir.path(), |m| {
        m["annotate"]["annotation"] = "r.annotation.parquet".into();
    });
    let pdf = dir.path().join("chart.pdf");
    let img = f
        .write_chart(Panel::Structure { of: Mixture::Types }, SIZE, 150.0, &pdf)
        .unwrap();
    assert_eq!((img.width(), img.height()), (SIZE.w, SIZE.h));
    let s = String::from_utf8_lossy(&std::fs::read(&pdf).unwrap()).into_owned();
    assert!(s.starts_with("%PDF-"));
    assert!(
        s.contains("(K0) Tj") && s.contains("(CT1) Tj"),
        "names are text"
    );
}
