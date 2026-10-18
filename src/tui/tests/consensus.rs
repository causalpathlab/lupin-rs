//! The consensus suggestion in the annotate TUI: shown under the
//! candidates and among the label choices, flagging no cluster, taken only
//! once confirmed, and recorded as the user's call on it; and
//! the trajectory's description in the order view.

use super::super::app::{App, Focus, Job, Running};
use super::super::round::decisions;
use crate::annotate::rounds::DecidedBy;
use crate::annotate_cmd::AnnotateCliArgs;
use crate::describe::consensus::tests::{fixture_round, Stub};
use crate::describe::consensus::{self, cache_path, Proposal};
use crate::describe::lm::{Lm, LmModel, LmSpec};
use crate::gene_text::cli::ComputeDevice;
use crate::manifest::run::derive_out_prefix;
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
use std::path::{Path, PathBuf};

fn press(app: &mut App, code: KeyCode) {
    app.key(KeyEvent::new(code, KeyModifiers::NONE));
}

const WHY: &str = "CT2 fits by its scores.";

fn spec() -> LmSpec {
    LmSpec {
        model: LmModel::Qwen3_0_6b,
        local_dir: None,
        device: ComputeDevice::Cpu,
    }
}

/// The round's cache.
fn cache_of(round: &Path) -> PathBuf {
    cache_path(&derive_out_prefix(&round.to_string_lossy()))
}

/// CT2 proposed for every cluster of `round`, cached as the job writes it.
fn propose_ct2(round: &Path) -> Vec<Proposal> {
    let loaded = crate::manifest::run::load(&round.to_string_lossy()).unwrap();
    let r = consensus::inputs(&loaded).unwrap();
    let load = || -> anyhow::Result<Box<dyn Lm>> {
        let cont = WHY.trim_start_matches("CT2");
        Ok(Box::new(Stub::new(&[("CT2", -0.1)], &[cont, cont, cont])))
    };
    consensus::propose_round(&r, &spec(), &cache_of(round), load).unwrap()
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

/// Cluster `id` has a suggestion on screen.
fn has(app: &App, id: u32) -> bool {
    app.round
        .as_ref()
        .is_some_and(|r| r.proposals.contains_key(&id))
}

#[test]
fn a_rounds_cached_suggestions_show_under_each_clusters_candidates_and_flag_none() {
    let dir = tempfile::tempdir().unwrap();
    let round = fixture_round(dir.path());
    propose_ct2(&round);
    let mut app = app_on(&round);
    for (sel, id) in [(0usize, 0u32), (1, 1), (2, 2)] {
        app.cluster_sel = sel;
        assert!(has(&app, id));
        let line = app.proposal_line().unwrap();
        assert!(line.starts_with("◆ suggested CT2 "), "{line}");
    }
    // A suggestion is no call: the clusters' flags stay the evidence's.
    let r = app.round.as_ref().unwrap();
    let flagged: Vec<bool> = r.clusters.iter().map(|c| c.flagged()).collect();
    assert_eq!(flagged, [false, false, false]);
}

#[test]
fn proposals_made_on_other_evidence_are_not_shown() {
    let dir = tempfile::tempdir().unwrap();
    let round = fixture_round(dir.path());
    let mut stale = propose_ct2(&round);
    std::fs::remove_file(cache_of(&round)).unwrap();
    stale[0].key = "another prompt's".into();
    consensus::append(&cache_of(&round), &stale[0]).unwrap();
    let app = app_on(&round);
    assert!(!has(&app, 0));
    assert!(app.proposal_line().is_none());
}

#[test]
fn the_strip_and_the_choices_show_the_proposal_and_p_takes_it_with_its_reason() {
    let dir = tempfile::tempdir().unwrap();
    let round = fixture_round(dir.path());
    propose_ct2(&round);
    let mut app = app_on(&round);
    let line = app.proposal_line().expect("K0's proposal");
    assert!(line.starts_with("◆ suggested CT2 "), "{line}");
    assert!(line.ends_with(&format!("— {WHY}")), "{line}");
    // The choices mark it, keyed as before, its reason in its detail.
    press(&mut app, KeyCode::Enter);
    let m = app.menu.as_ref().expect("the label choices");
    let at = m
        .choices
        .iter()
        .position(|c| c.label.starts_with("◆ CT2"))
        .expect("the proposal marked");
    assert_eq!(m.choices[at].key, Some('2'), "CT2 keeps its key");
    assert!(
        m.choices[at].detail.contains(WHY),
        "{}",
        m.choices[at].detail
    );
    assert_eq!(m.sel, 0, "Enter opens on the first choice");
    press(&mut app, KeyCode::Esc);
    // P opens them on the proposal; Enter takes it, its reason the model's.
    press(&mut app, KeyCode::Char('P'));
    assert_eq!(app.menu.as_ref().unwrap().sel, at);
    press(&mut app, KeyCode::Enter);
    assert_eq!(app.prompt.as_ref().unwrap().text, WHY);
    assert!(
        app.edits.is_empty(),
        "nothing is applied before the confirm"
    );
    press(&mut app, KeyCode::Enter);
    let r = app.round.as_ref().unwrap();
    let d = decisions(&app.edits, r, None);
    assert_eq!(d.len(), 1);
    assert_eq!(d[0].label.as_deref(), Some("CT2"));
    assert_eq!(d[0].decided_by, DecidedBy::AgentProposedUserAccepted);
    assert_eq!(d[0].rationale, WHY);
    let c = &d[0].evidence[0]["consensus"];
    assert_eq!(c["model"], "Qwen/Qwen3-0.6B");
    assert_eq!(c["revision"], LmModel::Qwen3_0_6b.revision());
    assert!(c["p"].as_f64().unwrap() > 0.5);
    assert!(
        d[0].evidence.iter().any(|e| e["term"] == "CT2"),
        "the cited type's row: {:?}",
        d[0].evidence
    );
    let alt: Vec<&str> = d[0]
        .alternatives
        .iter()
        .filter_map(|a| a["label"].as_str())
        .collect();
    assert!(alt.contains(&"CT1") && !alt.contains(&"CT2"), "{alt:?}");
    assert!(d[0].alternatives.iter().all(|a| a["p"].is_number()));
}

#[test]
fn another_label_while_a_proposal_is_shown_is_the_users_override() {
    let dir = tempfile::tempdir().unwrap();
    let round = fixture_round(dir.path());
    propose_ct2(&round);
    let mut app = app_on(&round);
    // The first candidate, CT1, over the proposal.
    press(&mut app, KeyCode::Char('1'));
    assert_ne!(app.prompt.as_ref().unwrap().text, WHY);
    press(&mut app, KeyCode::Enter);
    let d = decisions(&app.edits, app.round.as_ref().unwrap(), None);
    assert_eq!(d[0].label.as_deref(), Some("CT1"));
    assert_eq!(d[0].decided_by, DecidedBy::UserOverride);
    let first = &d[0].alternatives[0];
    assert_eq!(first["label"], "CT2");
    assert_eq!(first["proposed_by"], "consensus");
    assert_eq!(first["rationale"], WHY);
    // A round without proposals decides as before.
    let plain = fixture_round(&dir.path().join("plain").tap_mkdir());
    let mut app = app_on(&plain);
    press(&mut app, KeyCode::Char('1'));
    press(&mut app, KeyCode::Enter);
    let d = decisions(&app.edits, app.round.as_ref().unwrap(), None);
    assert_eq!(d[0].decided_by, DecidedBy::User);
    assert!(d[0].alternatives.is_empty());
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

#[test]
fn p_without_a_proposal_says_how_to_get_one() {
    let dir = tempfile::tempdir().unwrap();
    let round = fixture_round(dir.path());
    let mut app = app_on(&round);
    press(&mut app, KeyCode::Char('P'));
    assert!(app.menu.is_none());
    assert!(app.status.contains("L"), "{}", app.status);
}

#[test]
fn l_waits_for_a_running_job() {
    let dir = tempfile::tempdir().unwrap();
    let round = fixture_round(dir.path());
    let mut app = app_on(&round);
    let child = std::process::Command::new("sleep")
        .arg("30")
        .spawn()
        .unwrap();
    app.child = Some(Running::new(child, Job::Pass));
    press(&mut app, KeyCode::Char('L'));
    assert!(app.status.contains("wait"), "{}", app.status);
    assert!(matches!(
        app.child.as_ref().map(|r| &r.job),
        Some(Job::Pass)
    ));
    if let Some(mut r) = app.child.take() {
        let _ = r.child.kill();
        let _ = r.child.wait();
    }
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

#[test]
fn proposals_fill_in_as_the_running_job_appends_them() {
    let dir = tempfile::tempdir().unwrap();
    let round = fixture_round(dir.path());
    let made = propose_ct2(&round);
    let cache = cache_of(&round);
    std::fs::remove_file(&cache).unwrap();
    let mut app = app_on(&round);
    assert!(!has(&app, 0));
    app.follow_consensus().unwrap();
    let child = std::process::Command::new("sleep")
        .arg("30")
        .spawn()
        .unwrap();
    app.child = Some(Running::new(child, Job::Consensus(round.clone())));
    consensus::append(&cache, &made[0]).unwrap();
    app.tick();
    assert!(has(&app, 0), "K0's line is in");
    assert!(!has(&app, 2), "K2's is not yet");
    consensus::append(&cache, &made[2]).unwrap();
    app.tick();
    assert!(has(&app, 2));
    if let Some(r) = &mut app.child {
        let _ = r.child.kill();
    }
    finish(&mut app);
    assert!(app.failed.is_none());
}

#[test]
fn a_consensus_job_without_the_model_says_so_and_the_rest_works() {
    let dir = tempfile::tempdir().unwrap();
    let round = fixture_round(dir.path());
    let mut app = app_on(&round);
    app.follow_consensus().unwrap();
    app.push_log(format!("Error: {}", crate::describe::lm::OFFLINE));
    let child = std::process::Command::new("false").spawn().unwrap();
    app.child = Some(Running::new(child, Job::Consensus(round.clone())));
    finish(&mut app);
    assert!(app.failed.is_none(), "no failure popup: the rest works");
    assert!(
        app.status.contains("run once online or give --lm-dir"),
        "{}",
        app.status
    );
    // Labels still go on as before.
    press(&mut app, KeyCode::Char('1'));
    press(&mut app, KeyCode::Enter);
    assert_eq!(app.edits.len(), 1);
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
    std::fs::write(t.with_file_name("t.trajectory_lineages.tsv"), text).unwrap();
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
fn l_in_the_order_view_runs_no_model() {
    let dir = tempfile::tempdir().unwrap();
    let t = fixture_trajectory(dir.path());
    let mut app = order_on(dir.path(), &t);
    app.focus = Focus::Order;
    press(&mut app, KeyCode::Char('L'));
    assert!(app.child.is_none());
    assert!(app.status.contains("no model"), "{}", app.status);
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

#[test]
fn the_tuis_describe_job_writes_its_cache_and_no_description() {
    let lm = crate::describe::LmArgs::default();
    let argv = super::super::runner::describe_argv(Path::new("d/r.senna.json"), "d/r", &lm);
    assert_eq!(argv[..5], ["describe", "-f", "d/r.senna.json", "-o", "d/r"]);
    assert!(argv.iter().any(|a| a == "--lm"));
    assert!(
        argv.iter().any(|a| a == "--cache-only"),
        "never writes over a description the user made: {argv:?}"
    );
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
