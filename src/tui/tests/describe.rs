//! The trajectory's texts in the order view: written by code from the
//! run's outputs, shown at once, and read again when they change.

use super::super::app::{App, Focus, Job, Running};
use crate::annotate::rounds::{ClusterId, Method};
use crate::annotate_cmd::AnnotateCliArgs;
use std::path::{Path, PathBuf};

/// A round of three clusters of ten cells in `dir`: CT1–CT3 each lead their
/// own cluster on the enrichment's share, their five markers rising there;
/// cluster 1 has a decision in its history.
fn fixture_round(dir: &Path) -> PathBuf {
    use crate::manifest::rounds::write_clusters;
    use crate::manifest::run::RunManifest;
    use legume_numeric::matrix::dense_mat_io::Mat;
    use legume_numeric::matrix::traits::IoOps;
    let at = |f: &str| dir.join(f).to_string_lossy().into_owned();
    let (g, k, per) = (30usize, 3usize, 10usize);
    let genes: Vec<Box<str>> = (0..g).map(|i| format!("GENE{i}").into()).collect();
    let cells: Vec<Box<str>> = (0..k * per).map(|i| format!("c{i}").into()).collect();
    let ids: Vec<Option<ClusterId>> = (0..k * per).map(|i| Some((i / per) as ClusterId)).collect();
    write_clusters(&at("run.clusters.parquet"), &cells, &ids).unwrap();
    let labels: Vec<Box<str>> = (0..k * per)
        .map(|i| format!("CT{}", i / per + 1).into())
        .collect();
    graph_embedding_util::type_annotation::write_cell_labels(
        &at("run"),
        &cells,
        &labels,
        &vec![0.9; k * per],
    )
    .unwrap();
    let mut expr = Mat::zeros(g, k);
    for i in 0..g {
        for c in 0..k {
            expr[(i, c)] = 1.0;
        }
    }
    for c in 0..k {
        for m in 0..5 {
            expr[(c * 5 + m, c)] = 20.0 - m as f32;
        }
    }
    let kcols: Vec<Box<str>> = (0..k).map(|c| format!("K{c}").into()).collect();
    expr.to_parquet_with_names(
        &at("run.cluster_expression.parquet"),
        (Some(&genes), Some("gene")),
        Some(&kcols),
    )
    .unwrap();
    let mut panel = String::from("gene\tcelltype\n");
    for c in 0..k {
        for m in 0..5 {
            panel += &format!("GENE{}\tCT{}\n", c * 5 + m, c + 1);
        }
    }
    std::fs::write(dir.join("markers.tsv"), panel).unwrap();
    let types: Vec<Box<str>> = (0..k).map(|c| format!("CT{}", c + 1).into()).collect();
    let table = |own: f32, other: f32, file: &str| {
        let mut m = Mat::zeros(k, k);
        for r in 0..k {
            for c in 0..k {
                m[(r, c)] = if r == c { own } else { other };
            }
        }
        m.to_parquet_with_names(&at(file), (Some(&kcols), Some("cluster")), Some(&types))
            .unwrap();
    };
    table(0.7, 0.15, "run.cluster_celltype_q.parquet");
    table(0.001, 0.5, "run.cluster_celltype_q_values.parquet");
    std::fs::write(
        dir.join("run.annotation_history.json"),
        r#"{"1": [{"round": "run.senna.json", "action": "label", "label": "CT2",
            "rationale": "its markers rise", "decided_by": "user",
            "timestamp": "2026-01-01T00:00:00Z"}]}"#,
    )
    .unwrap();
    let mut m = RunManifest::new(crate::manifest::run::RunKind::Topic, "run");
    m.cluster.clusters = Some("run.clusters.parquet".into());
    let a = &mut m.annotate;
    a.argmax = Some("run.argmax.parquet".into());
    a.markers = Some("markers.tsv".into());
    a.cluster_expression = Some("run.cluster_expression.parquet".into());
    a.cluster_celltype_q = Some("run.cluster_celltype_q.parquet".into());
    a.cluster_celltype_q_values = Some("run.cluster_celltype_q_values.parquet".into());
    a.history = Some("run.annotation_history.json".into());
    a.called_by = Some(Method::Enrichment);
    let file = dir.join("run.senna.json");
    m.save(&file).unwrap();
    file
}

/// An app on the fixture round in `dir`.
fn app_on(round: &Path) -> App {
    #[derive(clap::Parser)]
    struct Cli {
        #[command(flatten)]
        annotate: AnnotateCliArgs,
    }
    let panel = round.with_file_name("markers.tsv");
    let argv = ["lupin", "-o", "x", "-m", panel.to_str().unwrap()];
    let args = <Cli as clap::Parser>::parse_from(argv).annotate;
    let types = crate::annotate::markers::read_panel(panel.to_str().unwrap()).unwrap();
    let tree = crate::manifest::ontology::panel_tree_on(None, &types);
    let mut app = App::new(args, round.to_path_buf(), round.to_path_buf(), tree);
    app.open(round);
    app.focus = Focus::Clusters;
    app
}

trait TapMkdir {
    fn tap_mkdir(self) -> Self;
}

impl TapMkdir for PathBuf {
    fn tap_mkdir(self) -> Self {
        std::fs::create_dir_all(&self).unwrap();
        self
    }
}

use super::super::app::TreeMode;
use crate::describe::trajectory::{self, tests::fixture_trajectory, OrderFacts};

/// Trajectory `t`'s texts as code writes them from its facts now.
fn texts_of(t: &Path) -> Vec<trajectory::OrderText> {
    let loaded = crate::manifest::run::load(&t.to_string_lossy()).unwrap();
    trajectory::texts(&OrderFacts::read(&loaded).unwrap())
}

/// Trajectory `t`'s run paragraph now.
fn paragraph_of(t: &Path) -> String {
    texts_of(t)[0].text.clone()
}

/// Write trajectory `t`'s lineages as `paths` (one `a > b > c` per line).
fn relineage(t: &Path, paths: &[&str]) {
    let mut text = String::from("# lineage\tcomponent\tpath\n");
    for (i, p) in paths.iter().enumerate() {
        text += &format!("L{i}\t0\t{p}\n");
    }
    let file = t.with_file_name("t.trajectory_lineages.tsv");
    let before = std::fs::metadata(&file).and_then(|m| m.modified()).ok();
    std::fs::write(&file, text).unwrap();
    // Later than before, whatever the clock's grain.
    if let Some(before) = before {
        let later = std::time::SystemTime::now().max(before + std::time::Duration::from_secs(1));
        std::fs::File::options()
            .write(true)
            .open(&file)
            .unwrap()
            .set_modified(later)
            .unwrap();
    }
}

#[test]
fn the_order_texts_notice_a_rewrite_of_the_same_size() {
    let dir = tempfile::tempdir().unwrap();
    let t = fixture_trajectory(dir.path());
    relineage(&t, &["HSC_like > MPP_like > ERY"]);
    let mut texts = super::super::describe::OrderTexts::of(&t).unwrap();
    assert!(!texts.poll(), "nothing changed");
    relineage(&t, &["HSC_like > MPP_like > MON"]);
    assert!(texts.poll());
    assert_eq!(texts.texts, texts_of(&t));
}

/// An app on a round in `dir/r`, the order view up on trajectory `t`.
fn order_on(dir: &Path, t: &Path) -> App {
    let mut app = app_on(&fixture_round(&dir.join("r").tap_mkdir()));
    app.source = t.to_path_buf();
    app.trajectory_tui = true;
    app.figures = None;
    app.toggle_order(false);
    app.figures = None;
    app
}

#[test]
fn the_order_view_shows_the_selected_types_line_or_the_runs_paragraph() {
    let dir = tempfile::tempdir().unwrap();
    let t = fixture_trajectory(dir.path());
    // Written from the run's outputs: there at once, no model run.
    let mut app = order_on(dir.path(), &t);
    assert!(app.child.is_none());
    app.focus = Focus::Clusters;
    let paragraph = paragraph_of(&t);
    assert_eq!(app.order_line().as_deref(), Some(paragraph.as_str()));
    app.focus = Focus::Order;
    let select = |app: &mut App, ty: &str| {
        let TreeMode::Order(v) = &mut app.tree_mode else {
            panic!("the order view")
        };
        if !v.types.iter().any(|(t, _)| t == ty) {
            v.types.push((ty.to_string(), 5));
        }
        v.sel = v.types.iter().position(|(t, _)| t == ty).unwrap();
    };
    // A type the trajectory does not order has no line: the paragraph.
    select(&mut app, "CT2");
    assert_eq!(app.order_line().as_deref(), Some(paragraph.as_str()));
    select(&mut app, "MPP_like");
    let mpp = texts_of(&t)
        .into_iter()
        .find(|x| x.subject.as_deref() == Some("MPP_like"))
        .unwrap();
    assert_eq!(app.order_line(), Some(mpp.text));
}

#[test]
fn another_trajectory_shown_shows_its_own_texts() {
    let dir = tempfile::tempdir().unwrap();
    let a = fixture_trajectory(&dir.path().join("a").tap_mkdir());
    let b = fixture_trajectory(&dir.path().join("b").tap_mkdir());
    relineage(&b, &["HSC_like > ERY"]);
    let mut app = order_on(dir.path(), &a);
    assert_eq!(app.order_texts.as_ref().map(|t| &t.manifest), Some(&a));
    app.source = b.clone();
    app.figures = None;
    app.load_order_texts(false);
    assert_eq!(app.order_texts.as_ref().map(|t| &t.manifest), Some(&b));
    app.focus = Focus::Clusters;
    assert_eq!(app.order_line(), Some(paragraph_of(&b)));
}

#[test]
fn opening_another_trajectory_shows_its_texts() {
    use crate::manifest::family::{Kind, Pick};
    let dir = tempfile::tempdir().unwrap();
    let a = fixture_trajectory(&dir.path().join("a").tap_mkdir());
    let b = fixture_trajectory(&dir.path().join("b").tap_mkdir());
    let mut app = order_on(dir.path(), &a);
    let pick = Pick {
        path: b.clone(),
        kind: Kind::Trajectory,
        labelled: false,
    };
    app.open_run(pick, true);
    app.tick();
    assert_eq!(app.order_texts.as_ref().map(|t| &t.manifest), Some(&b));
}

/// Tick until the job is taken in.
fn finish(app: &mut App) {
    for _ in 0..200 {
        app.tick();
        if app.child.is_none() {
            return;
        }
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
    panic!("the job never ended");
}

/// The order view's edges, as pairs.
fn order_edges(app: &App) -> Vec<(String, String)> {
    match &app.tree_mode {
        TreeMode::Order(v) => pairs(&v.edges),
        _ => panic!("the order view"),
    }
}

fn pairs(edges: &[crate::trajectory::edges::EdgeRow]) -> Vec<(String, String)> {
    (edges.iter())
        .map(|e| (e.a.to_string(), e.b.to_string()))
        .collect()
}

#[test]
fn a_new_output_shows_the_new_trajectorys_edges_with_its_texts() {
    use crate::trajectory::edges;
    let dir = tempfile::tempdir().unwrap();
    let a = fixture_trajectory(&dir.path().join("a").tap_mkdir());
    let b = fixture_trajectory(&dir.path().join("b").tap_mkdir());
    // B has fewer edges than A, and lineages of its own.
    let b_edges = b.with_file_name("t.trajectory_edges.parquet");
    let b_edges = b_edges.to_string_lossy();
    let fewer = edges::read(&b_edges).unwrap()[..2].to_vec();
    edges::write(&b_edges, &fewer).unwrap();
    let fewer = pairs(&fewer);
    relineage(&b, &["HSC_like > MPP_like"]);
    let (first, other) = (paragraph_of(&a), paragraph_of(&b));
    assert_ne!(first, other);
    let mut app = order_on(dir.path(), &a);
    app.focus = Focus::Clusters;
    app.tick();
    assert_eq!(app.order_line(), Some(first.clone()));
    assert_ne!(order_edges(&app).len(), fewer.len());
    let out = b.to_string_lossy().replace(".senna.json", "");
    app.set_output(&out);
    app.tick();
    assert_eq!(app.order_texts.as_ref().map(|t| &t.manifest), Some(&b));
    assert_eq!(app.order_line(), Some(other.clone()));
    assert_eq!(order_edges(&app), fewer, "the table is B's too");
    // Back to A, then to B while the view is away.
    let a_edges = pairs(
        &crate::trajectory::edges::read(
            &a.with_file_name("t.trajectory_edges.parquet")
                .to_string_lossy(),
        )
        .unwrap(),
    );
    app.set_output(&a.to_string_lossy().replace(".senna.json", ""));
    app.tick();
    assert_eq!(app.order_line(), Some(first));
    assert_eq!(order_edges(&app), a_edges);
    app.toggle_order(false);
    app.set_output(&out);
    app.toggle_order(false);
    app.figures = None;
    app.tick();
    assert_eq!(app.order_line(), Some(other));
    assert_eq!(order_edges(&app), fewer);
}

#[test]
fn the_order_view_entered_again_reads_the_texts_again() {
    let dir = tempfile::tempdir().unwrap();
    let a = fixture_trajectory(&dir.path().join("a").tap_mkdir());
    relineage(&a, &["HSC_like > MPP_like > ERY"]);
    let early = paragraph_of(&a);
    let mut app = order_on(dir.path(), &a);
    app.focus = Focus::Clusters;
    assert_eq!(app.order_line(), Some(early.clone()));
    app.toggle_order(false);
    // Run anew while the view is away, its lineages file the same size.
    relineage(&a, &["HSC_like > MPP_like > MON"]);
    let later = paragraph_of(&a);
    assert_ne!(later, early);
    app.toggle_order(false);
    app.focus = Focus::Clusters;
    assert_eq!(app.order_line(), Some(later));
    // So does opening the run again, the view up.
    relineage(&a, &["HSC_like > MPP_like > ERY"]);
    let pick = crate::manifest::family::Pick {
        path: a.clone(),
        kind: crate::manifest::family::Kind::Trajectory,
        labelled: false,
    };
    app.open_run(pick, true);
    assert_eq!(app.order_line(), Some(early));
}

#[test]
fn a_trajectory_run_again_on_its_output_leaves_no_stale_texts() {
    let dir = tempfile::tempdir().unwrap();
    let a = fixture_trajectory(&dir.path().join("a").tap_mkdir());
    let mut app = order_on(dir.path(), &a);
    // The run made again on the same output orders other types.
    relineage(&a, &["HSC_like > ERY"]);
    let child = std::process::Command::new("true").spawn().unwrap();
    app.child = Some(Running::new(child, Job::Trajectory(a.clone())));
    finish(&mut app);
    assert_eq!(
        app.order_texts.as_ref().map(|t| &t.texts),
        Some(&texts_of(&a))
    );
}

#[test]
fn the_order_views_edges_and_texts_come_from_one_trajectory() {
    use crate::trajectory::edges;
    let dir = tempfile::tempdir().unwrap();
    let a = fixture_trajectory(&dir.path().join("a").tap_mkdir());
    let b = fixture_trajectory(&dir.path().join("b").tap_mkdir());
    let b_edges = b.with_file_name("t.trajectory_edges.parquet");
    let b_edges = b_edges.to_string_lossy();
    edges::write(&b_edges, &edges::read(&b_edges).unwrap()[..2]).unwrap();
    relineage(&b, &["HSC_like > MPP_like"]);
    // A round made on A carries A's trajectory, borrowed; it is looked at
    // before B.
    let loaded = crate::manifest::run::load(&a.to_string_lossy()).unwrap();
    let borrowed = dir.path().join("r2").tap_mkdir().join("t.L1.senna.json");
    let copy = loaded.copy_to(borrowed.clone()).unwrap();
    copy.manifest.save(&copy.file).unwrap();
    let mut app = app_on(&fixture_round(&dir.path().join("r").tap_mkdir()));
    app.target = borrowed.clone();
    app.source = b;
    app.trajectory_tui = true;
    app.figures = None;
    app.toggle_order(false);
    app.figures = None;
    app.tick();
    let a_edges = pairs(
        &edges::read(
            &a.with_file_name("t.trajectory_edges.parquet")
                .to_string_lossy(),
        )
        .unwrap(),
    );
    assert_eq!(order_edges(&app), a_edges);
    assert_eq!(
        app.order_texts.as_ref().map(|t| &t.manifest),
        Some(&borrowed)
    );
    app.focus = Focus::Clusters;
    assert_eq!(app.order_line(), Some(paragraph_of(&a)));
}
