//! The screen's state and what each key does to it.
//!
//! Work goes cluster by cluster: pick a cluster, weigh its candidates and
//! genes, and give it a label (a candidate, any node of the ontology tree,
//! or none), each with a reason. Genes can be added to a cell type's
//! markers. Saving writes the edits as the next round.

use super::menu::{numbered, Action, Menu, Outcome};
use super::round::{decisions, Edit, RoundView};
use super::runner;
use crate::annotate::markers::label_key;
use crate::annotate::panel_tree::PanelTree;
use crate::annotate::rounds::ClusterId;
use crate::annotate_cmd::{AnnotateCliArgs, AnnotateMethod};
use crate::manifest::rounds::chain_rounds;
use crate::manifest::run::annotated_path;
use crate::trajectory::prior::Relation;
use enrichment::UNASSIGNED_LABEL;
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
use std::path::{Path, PathBuf};
use std::process::Child;
use std::sync::mpsc::{channel, Receiver, Sender};
use std::time::Instant;

/// Log lines kept for the log pane.
const LOG_KEEP: usize = 1000;
/// A cluster whose top candidate has less of its evidence than this is flagged.
pub const CONTESTED: f32 = 0.5;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Focus {
    Clusters,
    Genes,
    Tree,
    /// The order view's precedence table (the tree pane then shows the
    /// ontology beside it).
    Order,
    /// The selected cluster's GO terms, when the round scored them.
    Go,
}

/// The settings the screen edits, in the order it lists them.
#[derive(Clone, Copy)]
pub enum Setting {
    Method,
    Knn,
    Resolution,
    NumClusters,
    NumPerm,
    Go,
    /// The output prefix of passes; Enter edits it.
    Output,
}

pub const SETTINGS: [Setting; 7] = [
    Setting::Method,
    Setting::Knn,
    Setting::Resolution,
    Setting::NumClusters,
    Setting::NumPerm,
    Setting::Go,
    Setting::Output,
];

impl Setting {
    pub fn name(self) -> &'static str {
        match self {
            Self::Output => "output",
            Self::Method => "method",
            Self::Knn => "knn",
            Self::Resolution => "resolution",
            Self::NumClusters => "clusters (k)",
            Self::NumPerm => "permutations",
            Self::Go => "GO terms",
        }
    }

    pub fn value(self, a: &AnnotateCliArgs) -> String {
        match self {
            Self::Output => a.out.to_string(),
            Self::Method => match a.method {
                AnnotateMethod::Auto => "auto".into(),
                AnnotateMethod::Enrichment => "enrichment".into(),
                AnnotateMethod::Projection => "projection".into(),
            },
            Self::Knn => a.knn.map_or("default".into(), |k| k.to_string()),
            Self::Resolution => format!("{:.2}", a.resolution),
            Self::NumClusters => a.num_clusters.map_or("auto".into(), |k| k.to_string()),
            Self::NumPerm => a.num_perm.to_string(),
            Self::Go => match (&a.gaf, &a.gmt) {
                (Some(_), _) => "--gaf".into(),
                (_, Some(_)) => "--gmt".into(),
                _ if a.go => "on".into(),
                _ => "off".into(),
            },
        }
    }

    /// One step up (`up`) or down.
    fn adjust(self, a: &mut AnnotateCliArgs, up: bool) {
        let step = |x: usize, by: usize, min: usize| {
            if up {
                x + by
            } else {
                x.saturating_sub(by).max(min)
            }
        };
        match self {
            Self::Method => {
                use AnnotateMethod::{Auto, Enrichment, Projection};
                a.method = match (a.method, up) {
                    (Auto, true) | (Projection, false) => Enrichment,
                    (Enrichment, true) | (Auto, false) => Projection,
                    (Projection, true) | (Enrichment, false) => Auto,
                };
            }
            Self::Knn => a.knn = Some(step(a.knn.unwrap_or(15), 1, 2)),
            Self::Resolution => {
                let r = a.resolution + if up { 0.1 } else { -0.1 };
                a.resolution = (r.max(0.1) * 100.0).round() / 100.0;
            }
            Self::NumClusters => {
                a.num_clusters = match (a.num_clusters, up) {
                    (None, true) => Some(2),
                    (None, false) => None,
                    (Some(k), false) if k <= 2 => None,
                    (Some(k), _) => Some(step(k, 1, 2)),
                };
            }
            Self::NumPerm => a.num_perm = step(a.num_perm, 100, 0),
            // Named gene sets are the command line's; only `--go` toggles.
            Self::Go if a.gaf.is_none() && a.gmt.is_none() => a.go = !a.go,
            Self::Go | Self::Output => {}
        }
    }
}

/// What the tree pane shows.
pub enum TreeMode {
    /// The panel's types on the tree lupin built.
    Panel,
    /// The Cell Ontology around a term.
    Ontology(super::ontology::OntologyView),
    /// Which cell types precede which: the prior `lupin trajectory` builds.
    Order(super::order::OrderView),
}

/// What the genes pane lists.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum GeneView {
    /// The cluster's most specific genes.
    Specific,
    /// The markers of the cluster's label.
    Markers,
}

/// What an open prompt will do with its reason.
pub enum Pending {
    Label {
        cluster: ClusterId,
        label: String,
        /// A panel type with no CL term, and the term chosen for it: offer
        /// to remember the alias once the label is taken.
        remember: Option<(String, String)>,
    },
    /// Keep cluster `cluster`'s label; the text is why.
    Keep { cluster: ClusterId },
    /// Hide the genes the text's pattern matches from the specific genes.
    Hide,
    /// The name (the text) of a mixed label for `parts`.
    MixName { parts: Vec<String> },
    /// The cell type (the text) to add or drop `genes` as markers of.
    MarkerType { genes: Vec<String>, add: bool },
    /// Append `label → id` to the alias file `file`; the text is the note.
    Remember {
        label: String,
        id: String,
        file: PathBuf,
    },
    Markers {
        label: String,
        genes: Vec<String>,
        add: bool,
    },
    /// State that `from` precedes `to` (or that they are unrelated) in the
    /// project's precedence file; the text is why.
    Precedence {
        from: String,
        to: String,
        relation: Relation,
    },
    /// Move the export of `pdf` to the base name typed.
    Relocate { pdf: PathBuf },
    /// The output prefix (the text) for passes; then start one when
    /// `then_start`.
    Output { then_start: bool },
    /// Run `lupin trajectory` with the output prefix typed; `replace` is the
    /// prefix whose existing manifest the user has agreed to replace.
    TrajectoryOut { replace: Option<String> },
    /// A choice among options (see [`super::menu`]). Takes keys, not text.
    Choose(super::menu::Menu),
}

/// A file the main loop picks in the file browser, which takes the screen.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FileWant {
    /// A marker panel, to annotate a run that has no annotation.
    Markers,
    /// A `cell<TAB>type` labels file for a trajectory run.
    Labels,
    /// A precedence file for a trajectory run (`--prior`).
    Prior,
    /// A `label<TAB>CL:id` file mapping the labels to Cell Ontology terms.
    LabelCl,
}

/// A one-line prompt for a decision's reason, prefilled.
pub struct Prompt {
    pub title: String,
    pub text: String,
    pub pending: Pending,
}

/// A running job's progress: the latest `@progress` line and log line.
#[derive(Default, Clone)]
pub struct JobProgress {
    pub done: f64,
    pub total: f64,
    pub stage: String,
    /// The job's latest log line.
    pub last: String,
}

impl JobProgress {
    /// The share done, once the job has reported any.
    pub fn fraction(&self) -> Option<f64> {
        (self.total > 0.0).then(|| (self.done / self.total).clamp(0.0, 1.0))
    }
}

/// Time left after `elapsed` with `f` done, assuming the rest goes at the
/// same pace; `None` until there is enough to go on (2% and 2 s).
pub fn eta(f: f64, elapsed: std::time::Duration) -> Option<std::time::Duration> {
    ((0.02..1.0).contains(&f) && elapsed.as_secs_f64() >= 2.0)
        .then(|| elapsed.mul_f64((1.0 - f) / f))
}

/// `ETA 1m10s`, or `estimating…`.
pub fn eta_text(f: f64, elapsed: std::time::Duration) -> String {
    eta(f, elapsed).map_or("estimating…".into(), |d| {
        format!("ETA {}", duration_text(d))
    })
}

/// `45s`, `3m05s`, `1h02m`.
pub fn duration_text(d: std::time::Duration) -> String {
    let s = d.as_secs();
    match s {
        0..60 => format!("{s}s"),
        60..3600 => format!("{}m{:02}s", s / 60, s % 60),
        _ => format!("{}h{:02}m", s / 3600, (s % 3600) / 60),
    }
}

/// What a child process is doing.
#[derive(Clone, PartialEq, Eq)]
pub enum Job {
    Pass,
    /// Writing this many edits as the next round.
    Save(usize),
    /// `lupin trajectory`, writing this manifest.
    Trajectory(PathBuf),
}

impl Job {
    /// The job's name in the status line.
    fn name(&self) -> &'static str {
        match self {
            Self::Pass => "pass",
            Self::Save(_) => "save",
            Self::Trajectory(_) => "trajectory",
        }
    }
}

/// A `lupin relabel --preview` rescoring the round against `edits`.
pub struct Rescoring {
    child: Child,
    out: Receiver<String>,
    edits: Vec<Edit>,
    file: PathBuf,
}

impl Rescoring {
    fn stop(mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
        let _ = std::fs::remove_file(&self.file);
    }
}

/// A key that needs pressing twice, and what it will do.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Armed {
    Quit,
    Run,
    Stop,
}

pub struct App {
    /// The running job's progress, as its `@progress` lines report it.
    pub progress: Option<JobProgress>,
    /// The progress popup is hidden (`b`); the status line keeps the gist.
    pub progress_hidden: bool,
    /// A job that failed: its popup stays, red, until a key.
    pub failed: Option<String>,
    /// The terminal tells Shift+Enter from Enter (the kitty keyboard
    /// protocol): Shift+Enter then runs a pass and Enter only edits.
    pub shift_enter: bool,
    pub args: AnnotateCliArgs,
    /// The run manifest passes start from.
    pub source: PathBuf,
    /// Where a pass records itself (see [`crate::manifest::run::annotated_path`]).
    pub target: PathBuf,
    /// The source names its own clusters, which the Leiden knobs then do not touch.
    pub fixed_clusters: bool,
    pub tree: PanelTree,
    pub tree_mode: TreeMode,
    /// The Cell Ontology, when there is one.
    pub cl: Option<crate::annotate::celltype_tree::ClTerms>,
    /// Each panel type on a CL term with every term above it.
    pub panel_ancestry: Vec<(String, std::collections::BTreeSet<String>)>,
    /// Where the run's data files are found, and aliases remembered.
    pub data_search: crate::manifest::data_files::SearchPath,
    /// The marker panel the session started from: label key → genes.
    pub original: std::collections::BTreeMap<String, std::collections::BTreeSet<String>>,
    pub round: Option<RoundView>,
    pub edits: Vec<Edit>,
    pub focus: Focus,
    pub cluster_sel: usize,
    pub gene_sel: usize,
    /// The selected row of the GO pane.
    pub go_sel: usize,
    /// Words of the selected GO term's name slid out of view (← →).
    pub go_shift: usize,
    pub gene_view: GeneView,
    /// Genes kept out of the specific-genes view.
    pub hidden: super::genes::GeneFilter,
    /// Show the hidden genes anyway, dimmed.
    pub show_hidden: bool,
    /// Genes marked for adding or dropping as markers.
    pub marked: Vec<String>,
    /// Tree nodes' labels marked for a mixed label.
    pub tree_marked: Vec<String>,
    /// Mixed labels by name.
    pub mixed: super::ontology::Mixed,
    pub tree_sel: usize,
    /// The order view's types on the Cell Ontology, drawn beside the
    /// precedence table as the tree pane draws the panel.
    pub order_tree: Option<PanelTree>,
    pub order_tree_sel: usize,
    /// The full ontology in place of `order_tree` (`o`).
    pub order_cl: Option<super::ontology::OntologyView>,
    pub setting: usize,
    /// The clustering and pass settings popup is open.
    pub settings_open: bool,
    /// The key guide is open.
    pub help_open: bool,
    pub prompt: Option<Prompt>,
    pub log: Vec<String>,
    log_tx: Sender<String>,
    log_rx: Receiver<String>,
    pub child: Option<(Child, Instant, Job)>,
    /// A preview rescoring the round against the marker edits, running.
    pub rescoring: Option<Rescoring>,
    /// The marker edits the scores on screen reflect (none: the round's own).
    scored_for: Vec<Edit>,
    /// The round's own scores while rescored ones are shown.
    pub recorded: Option<super::round::Scores>,
    pub status: String,
    /// Settings changed since the round on screen was made.
    pub stale: bool,
    armed: Option<Armed>,
    pub quit: bool,
    /// The trajectory figures of the manifests on screen, when they have any.
    pub figures: Option<super::figure_pane::FigurePane>,
    /// The output prefix was given or confirmed; until then the first pass
    /// asks for it.
    pub out_chosen: bool,
    /// How the order view runs `lupin trajectory`.
    pub trajectory: super::order::TrajectoryRun,
    /// A file to pick in the file browser, which the main loop opens.
    pub want_file: Option<FileWant>,
    /// Once the pass running now is done, run the trajectory.
    trajectory_after_pass: bool,
    /// The terminal's picture protocol, asked for once.
    picker: Option<ratatui_image::picker::Picker>,
}

impl App {
    pub fn new(args: AnnotateCliArgs, source: PathBuf, target: PathBuf, tree: PanelTree) -> Self {
        let (log_tx, log_rx) = channel();
        Self {
            args,
            source,
            target,
            fixed_clusters: false,
            tree,
            tree_mode: TreeMode::Panel,
            order_tree: None,
            order_tree_sel: 0,
            order_cl: None,
            cl: None,
            panel_ancestry: Vec::new(),
            data_search: crate::manifest::data_files::SearchPath::new(None),
            original: std::collections::BTreeMap::new(),
            round: None,
            edits: Vec::new(),
            focus: Focus::Clusters,
            cluster_sel: 0,
            gene_sel: 0,
            go_sel: 0,
            go_shift: 0,
            gene_view: GeneView::Specific,
            hidden: super::genes::GeneFilter::default(),
            show_hidden: false,
            marked: Vec::new(),
            tree_marked: Vec::new(),
            mixed: super::ontology::Mixed::default(),
            tree_sel: 0,
            setting: 0,
            settings_open: false,
            help_open: false,
            prompt: None,
            log: Vec::new(),
            log_tx,
            log_rx,
            child: None,
            rescoring: None,
            scored_for: Vec::new(),
            recorded: None,
            status: "r: run a pass".into(),
            stale: false,
            armed: None,
            quit: false,
            shift_enter: false,
            figures: None,
            out_chosen: true,
            trajectory: super::order::TrajectoryRun::default(),
            want_file: None,
            trajectory_after_pass: false,
            picker: None,
            progress: None,
            progress_hidden: false,
            failed: None,
        }
    }

    pub fn push_log(&mut self, line: String) {
        self.log.push(line);
        if self.log.len() > LOG_KEEP {
            self.log.drain(..self.log.len() - LOG_KEEP);
        }
    }

    /// Take in finished work: new log lines, and a pass that has ended.
    pub fn tick(&mut self) {
        while let Ok(line) = self.log_rx.try_recv() {
            // A running job's progress lines feed the popup, not the log.
            if let Some((done, total, stage)) = crate::progress::parse(&line) {
                if self.child.is_some() {
                    let p = self.progress.get_or_insert_with(JobProgress::default);
                    (p.done, p.total, p.stage) = (done, total, stage);
                }
                continue;
            }
            if self.child.is_some() {
                self.progress.get_or_insert_with(JobProgress::default).last = line.clone();
            }
            self.push_log(line);
        }
        for line in super::drain_own_log() {
            self.push_log(line);
        }
        self.keep_scores_current();
        let Some((child, started, job)) = &mut self.child else {
            return;
        };
        let (secs, job) = (started.elapsed().as_secs(), job.clone());
        match child.try_wait() {
            Ok(None) => {
                let p = self.progress.get_or_insert_with(JobProgress::default);
                let elapsed = started.elapsed();
                self.status = match p.fraction() {
                    Some(f) if self.progress_hidden => format!(
                        "{} {:.0}% · {} · b shows  (x: stop)",
                        job.name(),
                        100.0 * f,
                        eta_text(f, elapsed)
                    ),
                    _ => format!("{}… {secs}s  (x: stop)", job.name()),
                };
            }
            Ok(Some(st)) if st.success() => {
                self.child = None;
                self.progress = None;
                if let Job::Trajectory(m) = job {
                    self.trajectory_done(&m, secs);
                    return;
                }
                let then_trajectory =
                    job == Job::Pass && std::mem::take(&mut self.trajectory_after_pass);
                // A new pass rewrote the base round: rounds made on the old
                // one no longer apply.
                if job == Job::Pass {
                    if let Err(e) = crate::manifest::rounds::supersede_later(&self.target) {
                        self.push_log(format!("could not set the old rounds aside: {e:#}"));
                    }
                }
                let (latest, _) = chain_rounds(&self.target);
                // Edits made while a save ran were not in it: carry them into
                // the new round (a save keeps the clusters' ids).
                let later = match job {
                    Job::Save(n) => self.edits.split_off(n.min(self.edits.len())),
                    _ => Vec::new(),
                };
                self.open(&latest);
                if matches!(self.tree_mode, TreeMode::Order(_)) {
                    self.reload_order();
                }
                let kept = later.len();
                self.edits = later;
                let done = match job {
                    Job::Save(n) if kept > 0 => format!(
                        "saved {n} edit(s), {kept} made since still unsaved; {}",
                        self.export()
                    ),
                    Job::Save(n) => format!("saved {n} edit(s); {}", self.export()),
                    _ => format!("pass done in {secs}s"),
                };
                self.status = format!("{done}. {}", self.status);
                if then_trajectory {
                    if self.prompt.is_none() {
                        self.ask_trajectory_out();
                    } else {
                        self.status += " r in the order view runs the trajectory.";
                    }
                }
            }
            Ok(Some(st)) => {
                self.child = None;
                self.progress = None;
                self.trajectory_after_pass = false;
                self.status = format!("{} failed ({st}); see the log", job.name());
                self.failed = Some(self.status.clone());
            }
            Err(e) => {
                self.child = None;
                self.progress = None;
                self.trajectory_after_pass = false;
                self.status = format!("lost the {}: {e}", job.name());
                self.failed = Some(self.status.clone());
            }
        }
    }

    /// The running job for its popup: what it is and how long it has run.
    pub fn running_job(&self) -> Option<(String, std::time::Duration)> {
        let (_, started, job) = self.child.as_ref()?;
        let what = match job {
            Job::Pass => format!("pass → {}", file_name(&self.target)),
            Job::Save(n) => format!("saving {n} edit(s)"),
            Job::Trajectory(m) => format!("trajectory → {}", file_name(m)),
        };
        Some((what, started.elapsed()))
    }

    /// Show round `manifest`, dropping unsaved edits.
    pub fn open(&mut self, manifest: &Path) {
        match RoundView::load(manifest) {
            Ok(r) => {
                let flagged = r.clusters.iter().filter(|c| c.flagged()).count();
                self.status = format!(
                    "{}: {} clusters, {flagged} to look at (]: next)",
                    file_name(manifest),
                    r.clusters.len()
                );
                self.stop_rescoring();
                self.recorded = None;
                self.scored_for.clear();
                self.round = Some(r);
                self.edits.clear();
                self.marked.clear();
                self.tree_marked.clear();
                self.gene_sel = 0;
                self.stale = false;
                self.cluster_sel = self.cluster_sel.min(self.n_clusters().saturating_sub(1));
            }
            Err(e) => self.status = format!("could not read {}: {e:#}", manifest.display()),
        }
    }

    /// The unsaved marker edits, in order.
    fn marker_edits(&self) -> Vec<Edit> {
        self.edits
            .iter()
            .filter(|e| matches!(e, Edit::Markers { .. }))
            .cloned()
            .collect()
    }

    /// Show scores that reflect the marker edits: take in a finished
    /// rescoring, restart one the edits have moved past, start one for new
    /// edits (not while a pass or save runs), and put the round's own
    /// scores back when no marker edit is left.
    fn keep_scores_current(&mut self) {
        let want = self.marker_edits();
        if let Some(r) = &mut self.rescoring {
            match r.child.try_wait() {
                Ok(None) if r.edits == want => return,
                Ok(None) => {}
                Ok(Some(st)) => {
                    let r = self.rescoring.take().expect("a rescoring");
                    let out = r.out.recv_timeout(std::time::Duration::from_secs(2));
                    let _ = std::fs::remove_file(&r.file);
                    if r.edits == want {
                        self.scored_for = want.clone();
                        match out.ok().filter(|_| st.success()) {
                            Some(json) => self.show_scores(&json),
                            None => self.rescore_failed("rescoring failed; see the log"),
                        }
                    }
                }
                Err(_) => {}
            }
            if self.rescoring.as_ref().is_some_and(|r| r.edits != want) {
                self.stop_rescoring();
            }
        }
        if self.rescoring.is_some() || want == self.scored_for || self.child.is_some() {
            return;
        }
        let Some(round) = &mut self.round else { return };
        if want.is_empty() {
            if let Some(old) = self.recorded.take() {
                round.swap_scores(old);
            }
            self.scored_for.clear();
            return;
        }
        if !round.rescorable {
            self.scored_for = want;
            self.status = "this round cached no statistics: marker edits are scored on save".into();
            return;
        }
        let file = std::env::temp_dir().join(format!(
            "lupin-rescore-{}-{}.jsonl",
            std::process::id(),
            self.edits.len()
        ));
        let started = (|| -> anyhow::Result<_> {
            let lines: Vec<String> = decisions(&want, round, None)
                .iter()
                .map(serde_json::to_string)
                .collect::<Result<_, _>>()?;
            std::fs::write(&file, lines.join("\n") + "\n")?;
            runner::spawn_preview(&round.manifest, &file, self.log_tx.clone())
        })();
        match started {
            // The table's title says it is rescoring; the status keeps what
            // the edit said.
            Ok((child, out)) => {
                self.rescoring = Some(Rescoring {
                    child,
                    out,
                    edits: want,
                    file,
                });
            }
            Err(e) => {
                self.scored_for = want;
                self.status = format!("could not rescore: {e:#}");
            }
        }
    }

    /// Show a preview's rescored candidates, keeping the round's own.
    fn show_scores(&mut self, json: &str) {
        let scores = match serde_json::from_str::<serde_json::Value>(json) {
            Ok(v) => super::round::parse_scores(&v),
            Err(_) => return self.rescore_failed("rescoring gave no readable result; see the log"),
        };
        let Some(scores) = scores else {
            return self
                .rescore_failed("this round cannot be rescored: marker edits are scored on save");
        };
        let Some(round) = &mut self.round else { return };
        let old = round.swap_scores(scores);
        // The round's own are what is kept, across several rescorings.
        if self.recorded.is_none() {
            self.recorded = Some(old);
        }
    }

    /// Say why, and show the round's own scores again rather than ones the
    /// edits have moved past.
    fn rescore_failed(&mut self, why: &str) {
        if let (Some(round), Some(old)) = (&mut self.round, self.recorded.take()) {
            round.swap_scores(old);
        }
        self.status = why.into();
    }

    fn stop_rescoring(&mut self) {
        if let Some(r) = self.rescoring.take() {
            r.stop();
        }
    }

    /// Whether the order view shows its figures in place of the table.
    fn figures_shown(&self) -> bool {
        self.figures.as_ref().is_some_and(|f| f.shown)
    }

    fn n_clusters(&self) -> usize {
        self.round.as_ref().map_or(0, |r| r.clusters.len())
    }

    pub fn selected(&self) -> Option<&super::round::ClusterView> {
        self.round.as_ref()?.clusters.get(self.cluster_sel)
    }

    /// The selected cluster's label with the edits applied.
    pub fn current_label(&self) -> Option<String> {
        let r = self.round.as_ref()?;
        r.label_of(self.selected()?.id, &self.edits)
    }

    /// The selected cluster's label, else its top candidate.
    fn label_or_top(&self) -> Option<String> {
        self.current_label()
            .or_else(|| self.selected()?.candidates.first().map(|c| c.label.clone()))
    }

    /// Select cluster `i`, clearing the gene selection.
    fn select_cluster(&mut self, i: usize) {
        self.cluster_sel = i;
        self.gene_sel = 0;
        self.go_sel = 0;
        self.go_shift = 0;
        self.marked.clear();
        self.tree_marked.clear();
    }

    /// Whether the round scored GO terms, so the GO pane is shown.
    #[must_use]
    pub fn has_go(&self) -> bool {
        self.round
            .as_ref()
            .is_some_and(|r| r.clusters.iter().any(|c| !c.terms.is_empty()))
    }

    /// The label of the tree pane's selected node: a panel node's, or a CL
    /// term's in the ontology view.
    fn tree_selected_label(&self) -> Option<String> {
        match (&self.tree_mode, &self.cl) {
            (TreeMode::Order(v), _) => v.selected().map(str::to_string),
            (TreeMode::Ontology(v), Some(cl)) => v
                .selected()
                .map(|r| super::ontology::term_label(cl, &self.tree, &r.id)),
            _ => self
                .tree
                .visible()
                .get(self.tree_sel)
                .map(|&i| self.tree.label(i).to_string()),
        }
    }

    /// Mark or unmark the selected node for a mixed label.
    fn toggle_tree_mark(&mut self) {
        let Some(l) = self.tree_selected_label() else {
            return;
        };
        match self.tree_marked.iter().position(|m| *m == l) {
            Some(i) => {
                self.tree_marked.remove(i);
            }
            None => self.tree_marked.push(l),
        }
    }

    /// Label the selected cluster by the marked nodes together: a mixed
    /// label, for a cluster the evidence cannot split between them. Asks
    /// what to call it first.
    fn mix_marked(&mut self) {
        if self.tree_marked.len() < 2 {
            self.status =
                "mark two or more nodes with space, then + gives the cluster a mixed label".into();
            return;
        }
        let mut parts = std::mem::take(&mut self.tree_marked);
        parts.sort();
        parts.dedup();
        self.prompt = Some(Prompt {
            title: format!(
                " name the mix of {} (e.g. an abbreviation): ",
                parts.join(", ")
            ),
            text: parts.join(MIX),
            pending: Pending::MixName { parts },
        });
    }

    /// Label the selected cluster `name`, standing for `parts`; a name that
    /// is not just the parts joined is remembered with them.
    fn ask_mixed(&mut self, name: String, parts: Vec<String>) {
        let said: Vec<String> = parts
            .iter()
            .map(|p| {
                let cand = self
                    .selected()
                    .and_then(|c| c.candidates.iter().find(|t| t.label == *p));
                match cand {
                    Some(t) => format!("{p} ({})", evidence(t)),
                    None => p.clone(),
                }
            })
            .collect();
        let reason = format!(
            "mixed: {}; the evidence does not separate them",
            said.join(" + ")
        );
        // Every mix is listed, so a `+` in a type's own name is never read as one.
        if let Some(file) = self.data_search.amend(super::ontology::MIXED) {
            if let Err(e) = self.mixed.add(&name, &parts, &file) {
                self.status = format!("could not record the mix: {e:#}");
            }
        }
        self.ask_label(name, reason, None);
    }

    /// Put the tree's cursor on node `i`, when it is visible.
    fn select_node(&mut self, i: usize) {
        self.tree_sel = self
            .tree
            .visible()
            .iter()
            .position(|&v| v == i)
            .unwrap_or(0);
    }

    /// The selected cluster's candidates `under` covers, as `type share`.
    fn pooled(&self, under: impl Fn(&str) -> bool) -> Vec<String> {
        self.selected()
            .map(|c| &c.candidates[..])
            .unwrap_or_default()
            .iter()
            .filter(|c| under(&c.label))
            .map(|c| format!("{} {:.2}", c.label, c.share))
            .collect()
    }

    /// Whether cluster `id` has an unsaved label or keep.
    pub fn edited(&self, id: ClusterId) -> bool {
        self.edits.iter().any(|e| e.cluster() == Some(id))
    }

    /// Whether cluster `id` is decided: now, or in an earlier round.
    pub fn decided(&self, id: ClusterId) -> bool {
        self.edited(id) || self.round.as_ref().is_some_and(|r| r.decided.contains(&id))
    }

    /// Keep the selected cluster's label, asking why.
    fn keep(&mut self) {
        let Some(c) = self.selected() else { return };
        let label = self
            .current_label()
            .unwrap_or_else(|| UNASSIGNED_LABEL.into());
        let reason = match c.candidates.iter().find(|t| t.label == label) {
            Some(t) => format!("keep {label} ({})", evidence(t)),
            None => format!("keep {label}"),
        };
        let id = c.id;
        let reason = self.rescored_note(reason);
        self.prompt = Some(Prompt {
            title: format!(" K{id} keeps {label}: why? "),
            text: reason,
            pending: Pending::Keep { cluster: id },
        });
    }

    fn start(&mut self) {
        if self.child.is_some() {
            return;
        }
        // A pass needs a marker panel: pick one, and the pass goes on
        // (set_markers starts it).
        if self.args.markers.is_empty() {
            self.status = "a pass needs a marker panel: pick one".into();
            self.want_file = Some(FileWant::Markers);
            return;
        }
        if !self.out_chosen {
            self.ask_output(true);
            return;
        }
        let (_, later) = chain_rounds(&self.target);
        let replaces = self.target.is_file();
        if (replaces || !self.edits.is_empty() || !later.is_empty())
            && self.armed != Some(Armed::Run)
        {
            self.armed = Some(Armed::Run);
            let mut what = Vec::new();
            if replaces {
                what.push(format!(
                    "replaces {} and its outputs",
                    self.target.display()
                ));
            }
            if !self.edits.is_empty() {
                what.push(format!("drops {} unsaved edit(s)", self.edits.len()));
            }
            if !later.is_empty() {
                what.push(format!("sets aside {} saved round(s)", later.len()));
            }
            self.status = format!("a new pass {}: r again to go on", what.join(", "));
            return;
        }
        self.push_log(format!(
            "── pass: lupin annotate {}",
            self.args.to_argv().join(" ")
        ));
        match runner::spawn_pass(&self.args, replaces, self.log_tx.clone()) {
            Ok(c) => {
                self.child = Some((c, Instant::now(), Job::Pass));
                self.status = "running…".into();
            }
            Err(e) => self.status = format!("{e:#}"),
        }
    }

    /// Stop a running pass or trajectory. A save is not stopped midway,
    /// which could leave a round half-written: it finishes on its own.
    fn stop(&mut self) {
        match self.child.as_ref().map(|c| c.2.clone()) {
            Some(Job::Pass | Job::Trajectory(_)) => {
                self.trajectory_after_pass = false;
                if let Some((mut c, _, _)) = self.child.take() {
                    let _ = c.kill();
                    let _ = c.wait();
                    self.status = "stopped".into();
                }
            }
            Some(Job::Save(_)) => self.status = "a save finishes on its own; wait for it".into(),
            None => {}
        }
    }

    /// How many of the edits a running save is writing: those stay as they
    /// are until it is done.
    fn saving(&self) -> usize {
        match self.child.as_ref().map(|c| c.2.clone()) {
            Some(Job::Save(n)) => n.min(self.edits.len()),
            _ => 0,
        }
    }

    /// Whether a pass is running, which is about to replace the clusters.
    fn pass_running(&self) -> bool {
        self.child.as_ref().is_some_and(|c| c.2 == Job::Pass)
    }

    fn save(&mut self) {
        if self.child.is_some() {
            return;
        }
        let Some(round) = &self.round else {
            self.status = "nothing to save: run a pass first".into();
            return;
        };
        if self.edits.is_empty() {
            self.status = "no edits to save".into();
            return;
        }
        let file = decisions_file(&round.manifest);
        let written = (|| -> anyhow::Result<()> {
            let lines: Vec<String> = decisions(&self.edits, round, self.recorded.as_ref())
                .iter()
                .map(serde_json::to_string)
                .collect::<Result<_, _>>()?;
            std::fs::write(&file, lines.join("\n") + "\n")?;
            Ok(())
        })();
        let started = written
            .and_then(|()| runner::spawn_relabel(&round.manifest, &file, self.log_tx.clone()));
        match started {
            Ok(c) => {
                self.push_log(format!(
                    "── save: lupin relabel -f {} -d {} --next",
                    round.manifest.display(),
                    file.display()
                ));
                self.child = Some((c, Instant::now(), Job::Save(self.edits.len())));
                self.status = "saving…".into();
            }
            Err(e) => self.status = format!("save failed: {e:#}"),
        }
    }

    /// Write the open round's cell annotation and marker table; says where.
    fn export(&mut self) -> String {
        let Some(r) = &self.round else {
            return "nothing to export".into();
        };
        match super::export::write(r, &self.tree, self.cl.as_ref(), &self.mixed, &self.original) {
            Ok((cells, markers)) => {
                self.push_log(format!(
                    "exported {} and {}",
                    cells.display(),
                    markers.display()
                ));
                format!("exported {}, {}", file_name(&cells), file_name(&markers))
            }
            Err(e) => format!("export failed: {e:#}"),
        }
    }

    /// Open a prompt to label the selected cluster `label`, then maybe
    /// offer to `remember` an alias.
    fn ask_label(&mut self, label: String, reason: String, remember: Option<(String, String)>) {
        let Some(c) = self.selected() else { return };
        self.prompt = Some(Prompt {
            title: format!(" K{} → {label}: why? ", c.id),
            text: reason,
            pending: Pending::Label {
                cluster: c.id,
                label,
                remember,
            },
        });
    }

    /// After adding markers of `label`, ask whether the selected cluster is
    /// a `label` too; the markers' reason is the label's.
    fn ask_label_too(&mut self, label: String, reason: String) {
        let Some(c) = self.selected() else { return };
        self.prompt = Some(Prompt {
            title: format!(
                " label K{} as {label} too? enter: yes (the text is why) · esc: markers only ",
                c.id
            ),
            text: reason,
            pending: Pending::Label {
                cluster: c.id,
                label,
                remember: None,
            },
        });
    }

    /// Label the selected cluster by its `k`-th candidate.
    fn pick_candidate(&mut self, k: usize) {
        let Some(c) = self.selected() else { return };
        let Some(pick) = c.candidates.get(k).cloned() else {
            return;
        };
        let reason = match c.candidates.first() {
            Some(top) if k > 0 => format!(
                "{} ({}) over the top candidate {} ({})",
                pick.label,
                evidence(&pick),
                top.label,
                evidence(top)
            ),
            _ => format!("top candidate {} ({})", pick.label, evidence(&pick)),
        };
        let reason = self.rescored_note(reason);
        self.ask_label(pick.label, reason, None);
    }

    /// `reason`, noting when the numbers it quotes are a rescoring's.
    fn rescored_note(&self, reason: String) -> String {
        if self.recorded.is_some() {
            format!("{reason}; rescored with the marker edits")
        } else {
            reason
        }
    }

    /// Label the selected cluster by tree node `i`.
    fn pick_node(&mut self, i: usize) {
        let label = self.tree.label(i).to_string();
        let under = self.tree.labels_under(i);
        let pooled = self.pooled(|t| under.contains(&label_key(t).as_str()));
        let reason = if pooled.is_empty() {
            format!("{label} from the ontology")
        } else {
            format!("{label} from the ontology, over {}", pooled.join(" + "))
        };
        self.ask_label(label, reason, None);
    }

    fn unassign(&mut self) {
        let Some(c) = self.selected() else { return };
        let reason = match c.candidates.first() {
            Some(top) => format!(
                "no convincing candidate (top {} {})",
                top.label,
                evidence(top)
            ),
            None => "no candidate".into(),
        };
        self.ask_label(UNASSIGNED_LABEL.into(), reason, None);
    }

    /// Drop the last marker edit not being saved.
    fn undo_marker(&mut self) {
        let saving = self.saving();
        let last = self.edits[saving..]
            .iter()
            .rposition(|e| matches!(e, Edit::Markers { .. }))
            .map(|i| i + saving);
        self.status = match last.map(|i| self.edits.remove(i)) {
            Some(Edit::Markers {
                label, genes, add, ..
            }) => format!(
                "undid {} {} {label}'s markers",
                if add { "adding" } else { "dropping" },
                genes.join(", ")
            ),
            _ => "no marker edit to undo".into(),
        };
    }

    /// Drop the selected cluster's edits.
    fn undo(&mut self) {
        let Some(id) = self.selected().map(|c| c.id) else {
            return;
        };
        let before = self.edits.len();
        self.edits.retain(|e| e.cluster() != Some(id));
        self.status = match before - self.edits.len() {
            0 => format!("no edit of K{id} to undo"),
            _ => format!("undid K{id}'s edit"),
        };
    }

    /// The genes the genes pane lists: the cluster's specific genes, or
    /// its label's markers with the edits applied.
    pub fn listed_genes(&self) -> Vec<String> {
        match self.gene_view {
            GeneView::Specific => self
                .selected()
                .map(|c| {
                    c.genes
                        .iter()
                        .filter(|g| self.show_hidden || !self.hidden.hides(&g.0))
                        .map(|g| g.0.clone())
                        .collect()
                })
                .unwrap_or_default(),
            GeneView::Markers => self
                .label_markers()
                .into_iter()
                .map(|(g, _, _)| g)
                .collect(),
        }
    }

    /// The selected cluster's label's markers (the edits applied), each with
    /// its log2 fold change in the cluster and whether an edit added it; the
    /// most specific first, unmeasured last.
    pub fn label_markers(&self) -> Vec<(String, Option<f32>, bool)> {
        let (Some(r), Some(c), Some(l)) = (&self.round, self.selected(), self.current_label())
        else {
            return Vec::new();
        };
        let mut v: Vec<(String, Option<f32>, bool)> = r
            .markers_of(&l, &self.edits)
            .into_iter()
            .map(|(g, added)| {
                let fc = r.fold_change(c.id, &g);
                (g, fc, added)
            })
            .collect();
        v.sort_by(|a, b| {
            let f = |x: &Option<f32>| x.unwrap_or(f32::NEG_INFINITY);
            f(&b.1).total_cmp(&f(&a.1)).then_with(|| a.0.cmp(&b.0))
        });
        v
    }

    /// The marked genes, else the selected one.
    fn chosen_genes(&self) -> Vec<String> {
        if self.marked.is_empty() {
            self.listed_genes()
                .get(self.gene_sel)
                .cloned()
                .into_iter()
                .collect()
        } else {
            self.marked.clone()
        }
    }

    /// Keep genes matching `pattern` out of the specific genes, remembered
    /// in the project's hidden-genes file.
    fn hide(&mut self, pattern: &str) {
        let Some(file) = super::genes::GeneFilter::file(&self.data_search) else {
            return;
        };
        self.status = match self.hidden.add(pattern, &file) {
            Ok(()) => format!(
                "hiding {pattern} (in {}; H shows hidden genes)",
                file.display()
            ),
            Err(e) => format!("could not hide {pattern}: {e:#}"),
        };
        let n = self.listed_genes().len();
        self.gene_sel = self.gene_sel.min(n.saturating_sub(1));
    }

    /// Open a prompt to add (`add`) or drop the marked genes, else the
    /// selected one, as markers of the selected cluster's label. Genes
    /// already in (for adding) or not in (for dropping) are left out.
    fn ask_markers(&mut self, add: bool) {
        // Markers belong to a cell type: the cluster's label, else its top
        // candidate, else one named here (a new one, for a type the panel lacks).
        let label = self
            .current_label()
            .or_else(|| self.selected()?.candidates.first().map(|c| c.label.clone()));
        match label {
            Some(l) => self.ask_markers_of(l, add),
            None => self.ask_marker_type(add),
        }
    }

    /// Ask which cell type (existing or new) the chosen genes are markers
    /// of, prefilled with the cluster's label.
    fn ask_marker_type(&mut self, add: bool) {
        let genes = self.chosen_genes();
        if genes.is_empty() {
            return;
        }
        self.prompt = Some(Prompt {
            title: format!(
                " {} {} as markers of which cell type? (a new name makes a new type) ",
                if add { "add" } else { "drop" },
                genes.join(", ")
            ),
            text: self.current_label().unwrap_or_default(),
            pending: Pending::MarkerType { genes, add },
        });
    }

    /// [`Self::ask_markers`] for cell type `label`.
    fn ask_markers_of(&mut self, label: String, add: bool) {
        let (Some(c), Some(r)) = (self.selected(), &self.round) else {
            return;
        };
        let chosen = self.chosen_genes();
        let current: Vec<String> = r
            .markers_of(&label, &self.edits)
            .into_iter()
            .map(|(g, _)| g)
            .collect();
        let genes: Vec<String> = chosen
            .into_iter()
            .filter(|g| current.contains(g) != add)
            .collect();
        if genes.is_empty() {
            self.status = if add {
                format!("already markers of {label}")
            } else {
                format!("not markers of {label}")
            };
            return;
        }
        let fc: Vec<String> = genes
            .iter()
            .map(|g| {
                r.fold_change(c.id, g)
                    .map_or_else(|| format!("{g} not measured"), |f| format!("{g} {f:+.1}"))
            })
            .collect();
        let (verb, reason) = if add {
            ("add", format!("specific to K{} ({})", c.id, fc.join(", ")))
        } else {
            (
                "drop",
                format!("not specific to K{} ({})", c.id, fc.join(", ")),
            )
        };
        self.prompt = Some(Prompt {
            title: format!(" {verb} {} as markers of {label}: why? ", genes.join(", ")),
            text: reason,
            pending: Pending::Markers { label, genes, add },
        });
    }

    fn prompt_key(&mut self, k: KeyEvent) {
        let Some(p) = &mut self.prompt else { return };
        if let Pending::Choose(menu) = &mut p.pending {
            match menu.key(k.code) {
                Outcome::Open => {}
                Outcome::Cancel => {
                    self.prompt = None;
                    self.trajectory_after_pass = false;
                    self.status = "no trajectory run started".into();
                }
                Outcome::Chosen(a) => {
                    self.prompt = None;
                    self.choose(a);
                }
            }
            return;
        }
        match k.code {
            KeyCode::Esc => {
                if matches!(
                    self.prompt.take().map(|p| p.pending),
                    Some(Pending::Output { then_start: true })
                ) {
                    self.trajectory_after_pass = false;
                }
            }
            KeyCode::Enter => {
                let edits_clusters = matches!(
                    p.pending,
                    Pending::Label { .. } | Pending::Keep { .. } | Pending::Markers { .. }
                );
                if edits_clusters && self.pass_running() {
                    self.status =
                        "a pass is running and will replace these clusters: wait for it".into();
                    return;
                }
                let Some(p) = self.prompt.take() else { return };
                let reason = p.text.trim().to_string();
                if reason.is_empty() {
                    self.status = "type something, or Esc to cancel".into();
                    self.prompt = Some(p);
                    return;
                }
                match p.pending {
                    Pending::Label {
                        cluster,
                        label,
                        remember,
                    } => {
                        self.edits.push(Edit::Label {
                            cluster,
                            label,
                            reason,
                        });
                        self.focus = Focus::Clusters;
                        if let (Some((alias, id)), Some(file)) = (
                            remember,
                            self.data_search.amend(crate::manifest::data_files::ALIASES),
                        ) {
                            self.prompt = Some(Prompt {
                                title: format!(
                                    " remember {alias} → {id} in {}? enter: yes (the text is the note) · esc: no ",
                                    file.display()
                                ),
                                text: "picked in the lupin annotate TUI".into(),
                                pending: Pending::Remember {
                                    label: alias,
                                    id,
                                    file,
                                },
                            });
                        }
                    }
                    Pending::Keep { cluster } => {
                        self.edits.push(Edit::Keep { cluster, reason });
                    }
                    Pending::Hide => self.hide(&reason),
                    Pending::MixName { parts } => self.ask_mixed(reason, parts),
                    Pending::MarkerType { genes, add } => {
                        self.marked = genes;
                        self.ask_markers_of(reason, add);
                    }
                    Pending::Precedence { from, to, relation } => {
                        let TreeMode::Order(v) = &self.tree_mode else {
                            return;
                        };
                        self.tree_marked.clear();
                        self.status = match v.record(&from, &to, relation, &reason) {
                            Ok(file) => format!(
                                "recorded {from} {} {to} in {}",
                                relation.as_str(),
                                file.display()
                            ),
                            Err(e) => format!("could not record the statement: {e:#}"),
                        };
                        self.reload_order();
                    }
                    Pending::Output { then_start } => {
                        self.set_output(&reason);
                        if then_start && self.out_chosen && self.args.out.as_ref() == reason {
                            self.start();
                            if self.child.is_some() {
                                self.settings_open = false;
                            }
                        }
                    }
                    Pending::TrajectoryOut { replace } => {
                        let exists = annotated_path(&self.source, &reason).exists();
                        if exists && replace.as_deref() != Some(reason.as_str()) {
                            // A name typed over the one offered: ask again.
                            self.prompt = Some(Prompt {
                                title: format!(
                                    " {reason}: its manifest exists: Enter replaces it "
                                ),
                                text: reason.clone(),
                                pending: Pending::TrajectoryOut {
                                    replace: Some(reason),
                                },
                            });
                        } else {
                            self.run_trajectory(&reason, exists);
                        }
                    }
                    // Handled above, key by key.
                    Pending::Choose(_) => {}
                    Pending::Relocate { pdf } => {
                        self.status = match &mut self.figures {
                            Some(v) => v
                                .exports
                                .relocate(&pdf, &reason)
                                .unwrap_or_else(|e| format!("could not move: {e:#}")),
                            None => "no figures".into(),
                        };
                    }
                    Pending::Remember { label, id, file } => {
                        self.status = match remember_alias(&file, &label, &id, &reason) {
                            Ok(()) => format!("remembered {label} → {id} in {}", file.display()),
                            Err(e) => format!("could not remember the alias: {e:#}"),
                        };
                    }
                    Pending::Markers { label, genes, add } => {
                        let moved = if add {
                            let saving = self.saving();
                            super::round::take_added(&mut self.edits, saving, &genes, &label)
                        } else {
                            Vec::new()
                        };
                        self.status = format!(
                            "{} {} {label}'s markers on save (the next round is rescored){}",
                            genes.len(),
                            if add { "added to" } else { "dropped from" },
                            if moved.is_empty() {
                                String::new()
                            } else {
                                format!("; moved from {}", moved.join(", "))
                            }
                        );
                        self.edits.push(Edit::Markers {
                            label: label.clone(),
                            genes,
                            add,
                            reason: reason.clone(),
                        });
                        self.marked.clear();
                        // Markers only change what the next round scores on:
                        // offer to label the cluster the genes came from too.
                        let other =
                            self.current_label().map(|l| label_key(&l)) != Some(label_key(&label));
                        if add && other {
                            self.ask_label_too(label, reason);
                        }
                    }
                }
            }
            KeyCode::Backspace => {
                p.text.pop();
            }
            KeyCode::Char(c) if !k.modifiers.contains(KeyModifiers::CONTROL) => p.text.push(c),
            _ => {}
        }
    }

    pub fn key(&mut self, k: KeyEvent) {
        if k.modifiers.contains(KeyModifiers::CONTROL) && k.code == KeyCode::Char('c') {
            let armed = self.armed.take();
            return self.request_quit(armed, "ctrl-c");
        }
        if self.prompt.is_some() {
            return self.prompt_key(k);
        }
        let beside = self.focus == Focus::Tree && matches!(self.tree_mode, TreeMode::Order(_));
        let view = match &mut self.tree_mode {
            TreeMode::Ontology(v) => Some(v),
            TreeMode::Order(_) if beside => self.order_cl.as_mut(),
            _ => None,
        };
        if let (Some(v), Some(cl)) = (view, &self.cl) {
            if let Some(text) = &mut v.typing {
                match k.code {
                    KeyCode::Esc => v.typing = None,
                    KeyCode::Enter => {
                        let q = std::mem::take(text);
                        v.typing = None;
                        let left = v.search(cl, &q);
                        self.status = if left {
                            format!(
                                "no term in the data matches {q:?}: {} in the whole ontology (←: back)",
                                v.rows.len()
                            )
                        } else {
                            format!("{} term(s) match {q:?} (←: back)", v.rows.len())
                        };
                    }
                    KeyCode::Backspace => {
                        text.pop();
                    }
                    KeyCode::Char(c) => text.push(c),
                    _ => {}
                }
                return;
            }
        }
        // A failed job's popup goes with any key, which does nothing else.
        if self.failed.take().is_some() {
            return;
        }
        // `b` hides the running job's popup, or brings it back.
        if k.code == KeyCode::Char('b') && self.child.is_some() {
            self.progress_hidden = !self.progress_hidden;
            return;
        }
        if self.help_open {
            // Any key closes the guide.
            self.help_open = false;
            return;
        }
        if k.code == KeyCode::Char('?') {
            self.help_open = true;
            return;
        }
        let armed = self.armed.take();
        if self.settings_open {
            return self.settings_key(k, armed);
        }
        // The focused pane's keys come first (in the genes pane, x hides).
        if self.pane_key(k.code) {
            return;
        }
        match k.code {
            KeyCode::Char('q') => self.request_quit(armed, "q"),
            // In the trajectory's view `r` runs the trajectory from any pane;
            // a pass (the cluster settings) is for a run with no annotation,
            // offered by that menu.
            KeyCode::Char('r') if matches!(self.tree_mode, TreeMode::Order(_)) => {
                self.ask_trajectory_out();
            }
            KeyCode::Char('r') => self.settings_open = true,
            // Stops a running pass or save, asked twice.
            KeyCode::Char('x') if self.child.is_some() => {
                if armed == Some(Armed::Stop) {
                    self.stop();
                } else {
                    self.armed = Some(Armed::Stop);
                    self.status = "x again stops it".into();
                }
            }
            KeyCode::Char('s') => self.save(),
            KeyCode::Char('e') => {
                self.status = if self.edits.is_empty() {
                    self.export()
                } else {
                    "save the edits first (s): the export is of a saved round".into()
                };
            }
            KeyCode::Tab => self.focus = self.cycle(true),
            KeyCode::BackTab => self.focus = self.cycle(false),
            _ => {}
        }
    }

    /// Quit on `key` (`q` or `ctrl-c`), asking again first while a pass runs
    /// or edits are unsaved; never during a save.
    fn request_quit(&mut self, armed: Option<Armed>, key: &str) {
        let busy = self.child.as_ref().map(|c| c.2.clone());
        let unsaved = match self.edits.len() {
            0 => String::new(),
            n => format!(" and {n} unsaved edit(s)"),
        };
        if matches!(busy, Some(Job::Save(_))) {
            self.status = format!("saving… {key} once it is done");
        } else if let (true, Some(job)) = (armed != Some(Armed::Quit), &busy) {
            self.armed = Some(Armed::Quit);
            self.status = format!(
                "a {} is running{unsaved}: {key} again stops it and quits",
                job.name()
            );
        } else if armed == Some(Armed::Quit) || self.edits.is_empty() {
            self.stop();
            self.stop_rescoring();
            self.quit = true;
        } else {
            self.armed = Some(Armed::Quit);
            self.status = format!(
                "{} unsaved edit(s): {key} again to quit, s to save",
                self.edits.len()
            );
        }
    }

    /// Keys in the clustering and pass settings popup.
    fn settings_key(&mut self, k: KeyEvent, armed: Option<Armed>) {
        let code = k.code;
        // Shift+Enter runs; plain Enter runs only where the terminal cannot
        // tell the two apart.
        let run = code == KeyCode::Enter
            && (k.modifiers.contains(KeyModifiers::SHIFT) || !self.shift_enter);
        match code {
            KeyCode::Esc | KeyCode::Char('r' | 'q') => self.settings_open = false,
            KeyCode::Up | KeyCode::Down => {
                step(&mut self.setting, SETTINGS.len(), code);
            }
            KeyCode::Left | KeyCode::Right => {
                SETTINGS[self.setting].adjust(&mut self.args, code == KeyCode::Right);
                self.stale = self.round.is_some();
            }
            KeyCode::Enter if !run || matches!(SETTINGS[self.setting], Setting::Output) => {
                if matches!(SETTINGS[self.setting], Setting::Output) {
                    self.ask_output(false);
                } else {
                    self.status = "shift+enter runs the pass".into();
                }
            }
            KeyCode::Enter => {
                self.armed = armed;
                self.start();
                // Closed once running; open while a confirmation is pending.
                self.settings_open = self.child.is_none();
            }
            _ => {}
        }
    }

    fn cycle(&self, forward: bool) -> Focus {
        use Focus::{Clusters, Genes, Go, Order, Tree};
        // The order view's columns, left to right: clusters (genes below
        // them), ordering, ontology.
        let mut order = if matches!(self.tree_mode, TreeMode::Order(_)) {
            vec![Clusters, Genes, Order, Tree]
        } else {
            vec![Clusters, Genes, Tree]
        };
        if self.has_go() && !matches!(self.tree_mode, TreeMode::Order(_)) {
            order.push(Go);
        }
        let i = order.iter().position(|f| *f == self.focus).unwrap_or(0);
        let n = order.len();
        order[if forward {
            (i + 1) % n
        } else {
            (i + n - 1) % n
        }]
    }

    /// The focused pane's binding for `code`, if it has one: `false` leaves
    /// the key to the global bindings. Moving in the pane's list never claims it.
    fn pane_key(&mut self, code: KeyCode) -> bool {
        // The order view's figures take their keys before the tree's.
        let order = matches!(self.tree_mode, TreeMode::Order(_));
        let on_figures =
            self.focus == Focus::Order || (self.focus == Focus::Tree && self.figures_shown());
        if order && on_figures && self.figure_key(code) {
            return true;
        }
        let move_in = |sel: &mut usize, n: usize| step(sel, n, code);
        match self.focus {
            Focus::Clusters => {
                let (before, mut sel) = (self.cluster_sel, self.cluster_sel);
                move_in(&mut sel, self.n_clusters());
                if sel != before {
                    self.select_cluster(sel);
                }
                match code {
                    KeyCode::Char(c @ '1'..='9') => {
                        let k = c as usize - '1' as usize;
                        let n = self.selected().map_or(0, |c| c.candidates.len());
                        if k < n {
                            self.pick_candidate(k);
                        } else {
                            self.status = format!("this cluster has {n} candidate(s)");
                        }
                    }
                    KeyCode::Char('u') => self.unassign(),
                    KeyCode::Char('k') => self.keep(),
                    KeyCode::Backspace | KeyCode::Delete => self.undo(),
                    KeyCode::Char(']') => self.next_flagged(),
                    KeyCode::Enter => self.jump_to_tree(),
                    _ => return false,
                }
            }
            Focus::Genes => {
                let listed = self.listed_genes();
                move_in(&mut self.gene_sel, listed.len());
                match code {
                    KeyCode::Char(' ') => {
                        if let Some(g) = listed.get(self.gene_sel).cloned() {
                            match self.marked.iter().position(|m| *m == g) {
                                Some(i) => {
                                    self.marked.remove(i);
                                }
                                None => self.marked.push(g),
                            }
                        }
                    }
                    KeyCode::Char('a') => self.ask_markers(true),
                    KeyCode::Char('A') => self.ask_marker_type(true),
                    KeyCode::Char('x') => {
                        let genes = self.chosen_genes();
                        for g in genes {
                            self.hide(&g);
                        }
                        self.marked.clear();
                    }
                    KeyCode::Char('X') => {
                        if let Some(g) = self.chosen_genes().first() {
                            self.prompt = Some(Prompt {
                                title: " hide genes matching (`*` = anything): ".into(),
                                text: super::genes::suggest(g),
                                pending: Pending::Hide,
                            });
                        }
                    }
                    KeyCode::Char('H') => {
                        self.show_hidden = !self.show_hidden;
                        self.gene_sel = 0;
                    }
                    KeyCode::Char('d') => self.ask_markers(false),
                    KeyCode::Backspace | KeyCode::Delete => self.undo_marker(),
                    KeyCode::Char('m') => {
                        self.gene_view = match self.gene_view {
                            GeneView::Specific => GeneView::Markers,
                            GeneView::Markers => GeneView::Specific,
                        };
                        self.gene_sel = 0;
                        self.marked.clear();
                    }
                    KeyCode::Esc => self.focus = Focus::Clusters,
                    _ => return false,
                }
            }
            Focus::Go => {
                let before = self.go_sel;
                let n = self.selected().map_or(0, |c| c.terms.len());
                step(&mut self.go_sel, n, code);
                if self.go_sel != before {
                    self.go_shift = 0;
                }
                let words = self
                    .selected()
                    .and_then(|c| c.terms.get(self.go_sel))
                    .map_or(0, |t| t.term.split_whitespace().count());
                match code {
                    KeyCode::Right if self.go_shift + 1 < words => self.go_shift += 1,
                    KeyCode::Left => self.go_shift = self.go_shift.saturating_sub(1),
                    KeyCode::Esc => self.focus = Focus::Clusters,
                    _ => return false,
                }
            }
            Focus::Order => match code {
                KeyCode::Char('t') => self.toggle_order(false),
                KeyCode::Char(' ') => self.toggle_tree_mark(),
                _ => {
                    let took = self.order_key(code);
                    self.sync_tree_to_order();
                    return took;
                }
            },
            Focus::Tree if matches!(self.tree_mode, TreeMode::Order(_)) => {
                return self.order_tree_key(code);
            }
            Focus::Tree if code == KeyCode::Char('t') => self.toggle_order(false),
            Focus::Tree if code == KeyCode::Char('o') => self.toggle_ontology(),
            // A search is of the ontology: open it there.
            Focus::Tree
                if code == KeyCode::Char('/') && matches!(self.tree_mode, TreeMode::Panel) =>
            {
                self.toggle_ontology();
                if let TreeMode::Ontology(v) = &mut self.tree_mode {
                    v.typing = Some(String::new());
                }
            }
            Focus::Tree if code == KeyCode::Char(' ') => self.toggle_tree_mark(),
            Focus::Tree if code == KeyCode::Char('+') => self.mix_marked(),
            Focus::Tree if matches!(self.tree_mode, TreeMode::Ontology(_)) => {
                return self.ontology_key(code);
            }
            Focus::Tree => {
                let visible = self.tree.visible();
                move_in(&mut self.tree_sel, visible.len());
                let at = visible.get(self.tree_sel).copied();
                match (code, at) {
                    (KeyCode::Enter, Some(i)) => self.pick_node(i),
                    (KeyCode::Left, Some(i)) => {
                        if self.tree.nodes[i].children.is_empty() || self.tree.is_folded(i) {
                            // Up to the parent, like a file tree.
                            if let Some(p) = self.tree.nodes[i].parent {
                                self.select_node(p);
                            }
                        } else {
                            self.tree.fold(i, true);
                        }
                    }
                    (KeyCode::Right, Some(i)) => self.tree.fold(i, false),
                    (KeyCode::Esc, _) => self.focus = Focus::Clusters,
                    _ => return false,
                }
            }
        }
        true
    }

    /// Switch the tree pane to the order view (which types precede which) and
    /// back; `show_figures` opens on the figures (`lupin trajectory`).
    pub(super) fn toggle_order(&mut self, show_figures: bool) {
        self.tree_marked.clear();
        if matches!(self.tree_mode, TreeMode::Order(_)) {
            self.tree_mode = TreeMode::Panel;
            if self.focus == Focus::Order {
                self.focus = Focus::Tree;
            }
            return;
        }
        self.focus = Focus::Order;
        if self.figures.is_none() {
            let manifests = self.manifests();
            let refs: Vec<&Path> = manifests.iter().map(PathBuf::as_path).collect();
            let picker = self.picker().clone();
            match super::figure_pane::FigurePane::load(&refs, &picker) {
                Ok(v) => self.figures = v,
                Err(e) => self.push_log(format!("[WARN] trajectory figures: {e:#}")),
            }
            if let Some(v) = &mut self.figures {
                v.shown = show_figures;
            }
        }
        self.reload_order();
        self.status = "order view: space marks a type, then > states the first precedes the second, - that they are unrelated".into();
    }

    /// The manifests whose trajectory outputs the order view shows: the
    /// trajectory run made here, the round's target, then the run opened.
    fn manifests(&self) -> Vec<PathBuf> {
        let mut v: Vec<PathBuf> = self
            .trajectory
            .out
            .iter()
            .map(|o| annotated_path(&self.source, o))
            .collect();
        v.extend([self.target.clone(), self.source.clone()]);
        v
    }

    /// Ask for the output prefix of passes, prefilled with the one in use;
    /// start a pass once given when `then_start`.
    fn ask_output(&mut self, then_start: bool) {
        self.prompt = Some(Prompt {
            title: " output prefix for this run's rounds: ".into(),
            text: self.args.out.to_string(),
            pending: Pending::Output { then_start },
        });
    }

    /// Take `out` as the output prefix of passes; its latest round, when it
    /// has one, is opened.
    fn set_output(&mut self, out: &str) {
        if self.child.is_some() {
            self.status = "wait for the running job before changing the output".into();
            return;
        }
        if out == self.args.out.as_ref() {
            self.out_chosen = true;
            return;
        }
        if !self.edits.is_empty() {
            self.status = format!(
                "{} unsaved edit(s): save (s) or drop them before changing the output",
                self.edits.len()
            );
            return;
        }
        self.args.out = out.into();
        self.target = annotated_path(&self.source, out);
        self.out_chosen = true;
        if self.target.is_file() {
            let (latest, _) = chain_rounds(&self.target);
            self.open(&latest);
        }
        self.status = format!("writing under {out}");
    }

    /// Ask for the trajectory run's output prefix: the last one used, else
    /// the next free `{stem}.T{k}`.
    fn ask_trajectory_out(&mut self) {
        if self.child.is_some() {
            self.status = "wait for the running job, or stop it with x".into();
            return;
        }
        if !self.edits.is_empty() {
            self.status = format!(
                "{} unsaved edit(s): save them (s) so the trajectory sees the labels on screen",
                self.edits.len()
            );
            return;
        }
        let loaded = crate::manifest::run::load(&self.source.to_string_lossy()).ok();
        let mut has_labels = self
            .trajectory
            .argv
            .iter()
            .any(|a| a == "--labels" || a.starts_with("--labels="))
            || self.trajectory.labels.is_some()
            || self.round.is_some()
            || loaded
                .as_ref()
                .is_some_and(|l| l.manifest.annotate.argmax.is_some());
        // A trajectory run made from a labels file recorded it: reuse it.
        if !has_labels {
            if let Some(path) = loaded.as_ref().and_then(recorded_labels) {
                self.push_log(format!("labels from the run's last trajectory: {path}"));
                self.trajectory.labels = Some(path);
                has_labels = true;
            }
        }
        if !has_labels {
            self.ask_labels(None);
            return;
        }
        if self.ask_root() {
            return;
        }
        let text = self
            .trajectory
            .out
            .clone()
            .unwrap_or_else(|| super::order::default_out(&self.source));
        let exists = annotated_path(&self.source, &text).exists();
        self.prompt = Some(Prompt {
            title: if exists {
                format!(" run the trajectory as {text} (its manifest exists: Enter replaces it): ")
            } else {
                " run the trajectory: output prefix ".into()
            },
            pending: Pending::TrajectoryOut {
                replace: exists.then(|| text.clone()),
            },
            text,
        });
    }

    /// Show `menu`.
    fn show(&mut self, menu: Menu) {
        self.prompt = Some(Prompt {
            title: String::new(),
            text: String::new(),
            pending: Pending::Choose(menu),
        });
    }

    /// A run with no labels: annotate it here, or take a labels file;
    /// `note` heads the question (a pass that asks to be confirmed).
    fn ask_labels(&mut self, note: Option<String>) {
        let question = match note {
            Some(n) => format!("{n}. Choose annotate again to go on."),
            None => "This run has no cell-type labels yet, and the trajectory needs them. \
                     Where should they come from?"
                .into(),
        };
        self.show(Menu::new(
            question,
            numbered(vec![
                (
                    "Annotate this run now".into(),
                    "pick a marker panel (or use the run's), run an annotation pass, then the trajectory on its labels".into(),
                    Action::Annotate,
                ),
                (
                    "Use a cell<TAB>type labels file".into(),
                    "pick the file in the file browser; the trajectory reads its labels".into(),
                    Action::LabelsFile,
                ),
            ]),
        ));
    }

    /// When nothing would order the types (no statement makes an edge, and
    /// no `--root` or `--prior` was given), ask how the order should come
    /// instead of starting a run that can only fail; `true` when it asked.
    fn ask_root(&mut self) -> bool {
        if self.order_given() {
            return false;
        }
        let TreeMode::Order(v) = &self.tree_mode else {
            return false;
        };
        if v.error.is_some() || v.types.is_empty() || !v.edges().is_empty() {
            return false;
        }
        self.ask_order(None);
        true
    }

    /// `--root` or `--prior` among the trajectory's options.
    fn order_given(&self) -> bool {
        self.trajectory.argv.iter().any(|a| {
            ["--root", "--prior"]
                .iter()
                .any(|f| a == f || a.starts_with(&format!("{f}=")))
        })
    }

    /// The menu of ways to give the types an order; `note` heads it.
    fn ask_order(&mut self, note: Option<&str>) {
        let TreeMode::Order(v) = &self.tree_mode else {
            return;
        };
        let names: Vec<&str> = v.types.iter().map(|(t, _)| t.as_str()).collect();
        let why = match &self.cl {
            Some(cl) => {
                let (mapped, _) = cl.map_labels(names.iter().copied());
                match mapped.len() {
                    0 => format!("none of the {} labels map to Cell Ontology terms", names.len()),
                    n => format!(
                        "{n} of {} labels map to Cell Ontology terms, with no develops-from link between them",
                        names.len()
                    ),
                }
            }
            None => "no Cell Ontology is at hand".into(),
        };
        let question = format!(
            "{}Nothing orders these types yet ({why}). How should the trajectory get its order?",
            note.map_or(String::new(), |n| format!("{n}. "))
        );
        self.show(Menu::new(
            question,
            numbered(vec![
                (
                    "Start from one type (--root)".into(),
                    "choose the type the trajectory starts from; every other type is ordered after it by the data".into(),
                    Action::PickRoot,
                ),
                (
                    "State the order myself".into(),
                    "in the order table: space marks two types, > says the first precedes the second, - that they are unrelated; then r".into(),
                    Action::StateOrder,
                ),
                (
                    "Use a precedence file (--prior)".into(),
                    "pick a from<TAB>to<TAB>precedes|unrelated file, or an earlier run's trajectory_prior.tsv".into(),
                    Action::PriorFile,
                ),
                (
                    "Map my labels to Cell Ontology terms (--label-cl)".into(),
                    "pick a label<TAB>CL:id file; the ontology's develops-from links then order the types".into(),
                    Action::LabelCl,
                ),
            ]),
        ));
    }

    /// The menu of types to start from, with their cell counts.
    fn ask_root_type(&mut self) {
        let TreeMode::Order(v) = &self.tree_mode else {
            return;
        };
        let items = v
            .types
            .iter()
            .map(|(t, n)| {
                (
                    t.clone(),
                    format!("start from {t} ({n} cells), kept for later runs"),
                    Action::Root(t.clone()),
                )
            })
            .collect();
        let mut menu = Menu::new(
            "Which type does the trajectory start from?",
            numbered(items),
        );
        menu.sel = v.sel.min(menu.choices.len().saturating_sub(1));
        self.show(menu);
    }

    /// Do what the chosen option says.
    pub(super) fn choose(&mut self, action: Action) {
        match action {
            Action::PickRoot => self.ask_root_type(),
            Action::Root(t) => {
                self.push_log(format!("trajectory root: {t} (kept for later runs)"));
                self.trajectory.argv.extend(["--root".into(), t]);
                self.ask_trajectory_out();
            }
            Action::StateOrder => {
                self.focus = Focus::Order;
                if let Some(f) = &mut self.figures {
                    f.shown = false;
                }
                self.status = "space marks two types, > says the first precedes the second, - that they are unrelated; then r runs".into();
            }
            Action::PriorFile => self.want_file = Some(FileWant::Prior),
            Action::LabelCl => self.want_file = Some(FileWant::LabelCl),
            Action::LabelsFile => self.want_file = Some(FileWant::Labels),
            Action::Annotate => {
                self.annotate_then_trajectory();
                // A pass that replaces a round asks first; the menu asks
                // again, and annotate a second time goes on.
                if self.armed == Some(Armed::Run) && self.child.is_none() {
                    let note = self
                        .status
                        .trim_end_matches(": r again to go on")
                        .to_string();
                    self.ask_labels(Some(note));
                }
            }
        }
    }

    /// Take `path` (picked in the file browser) as the trajectory's
    /// precedence file, then ask for the output prefix.
    pub fn set_prior(&mut self, path: Option<&Path>) {
        let Some(p) = path else {
            self.status = "no precedence file picked".into();
            return;
        };
        let p = p.to_string_lossy().into_owned();
        self.push_log(format!("trajectory prior: {p} (kept for later runs)"));
        self.trajectory.argv.extend(["--prior".into(), p]);
        self.ask_trajectory_out();
    }

    /// Take `path` as the labels' Cell Ontology terms: reload the ontology
    /// with it and the order view; go on to the run when it now orders
    /// something, else ask again.
    pub fn set_label_cl(&mut self, path: Option<&Path>) {
        let Some(p) = path else {
            self.status = "no label<TAB>CL:id file picked".into();
            return;
        };
        let p = p.to_string_lossy().into_owned();
        let dir = self.source.parent().map(Path::to_path_buf);
        let terms = crate::manifest::ontology::load(
            dir.as_deref(),
            self.args.obo.as_deref(),
            Some(&p),
            crate::manifest::data_files::Fetch::Allowed,
        )
        .and_then(crate::manifest::data_files::ClData::into_terms);
        match terms {
            Ok(t) => {
                self.cl = t;
                if let Some(cl) = &self.cl {
                    self.panel_ancestry = super::ontology::type_ancestry(cl, &self.tree);
                }
                self.args.label_cl = Some(p.clone().into());
                // One --label-cl: the new file replaces one given before.
                let argv = std::mem::take(&mut self.trajectory.argv);
                let mut it = argv.into_iter();
                while let Some(a) = it.next() {
                    if a == "--label-cl" {
                        it.next();
                    } else if !a.starts_with("--label-cl=") {
                        self.trajectory.argv.push(a);
                    }
                }
                self.trajectory
                    .argv
                    .extend(["--label-cl".into(), p.clone()]);
                self.push_log(format!("label to Cell Ontology terms: {p}"));
                self.reload_order();
                if !self.ask_root() {
                    self.ask_trajectory_out();
                } else {
                    self.ask_order(Some("Still no order from the ontology"));
                }
            }
            Err(e) => self.status = format!("{p}: {e:#}"),
        }
    }

    /// Start `lupin trajectory` from the round on screen (else the run
    /// opened) into `out`.
    fn run_trajectory(&mut self, out: &str, replace: bool) {
        let from = self
            .round
            .as_ref()
            .map_or_else(|| self.source.clone(), |r| r.manifest.clone());
        let mut argv = vec![
            "-f".to_string(),
            from.to_string_lossy().into_owned(),
            "-o".into(),
            out.into(),
        ];
        argv.extend(self.trajectory.argv.iter().cloned());
        // A labels file typed in the TUI stands in for a missing annotation:
        // a round on screen has its own.
        if let (None, Some(l)) = (&self.round, &self.trajectory.labels) {
            argv.extend(["--labels".into(), l.clone()]);
        }
        if replace {
            argv.push("--overwrite".into());
        }
        self.push_log(format!("── lupin trajectory {}", argv.join(" ")));
        match runner::spawn_trajectory(&argv, self.log_tx.clone()) {
            Ok(c) => {
                let manifest = annotated_path(&from, out);
                self.child = Some((c, Instant::now(), Job::Trajectory(manifest)));
                self.trajectory.out = Some(out.into());
                self.status = "running the trajectory…".into();
            }
            Err(e) => self.status = format!("{e:#}"),
        }
    }

    /// The terminal's picture protocol, asked for the first time figures are
    /// shown and reused after.
    fn picker(&mut self) -> &ratatui_image::picker::Picker {
        let graphics = self.args.graphics;
        self.picker
            .get_or_insert_with(|| super::figure_pane::picker(graphics))
    }

    /// Annotate the run, then run the trajectory on its labels: the marker
    /// panel is the manifest's, else picked in the file browser.
    fn annotate_then_trajectory(&mut self) {
        self.trajectory_after_pass = true;
        if self.args.markers.is_empty() {
            self.want_file = Some(FileWant::Markers);
            return;
        }
        self.start();
    }

    /// Take `path` as the marker panel (picked in the file browser), then go
    /// on with the pass a missing annotation asked for.
    pub fn set_markers(&mut self, path: Option<&Path>) {
        let Some(path) = path else {
            self.trajectory_after_pass = false;
            self.status = "no marker panel picked".into();
            return;
        };
        let p = path.to_string_lossy();
        let read = crate::annotate::markers::read_panel(&p)
            .and_then(|panel| Ok((panel, super::round::panel_sets(&p)?)));
        match read {
            Ok((panel, sets)) => {
                self.tree = crate::manifest::ontology::panel_tree_on(self.cl.as_ref(), &panel);
                if let Some(cl) = &self.cl {
                    self.panel_ancestry = super::ontology::type_ancestry(cl, &self.tree);
                }
                self.original = sets;
                self.args.markers = p.into();
                self.start();
            }
            Err(e) => {
                self.trajectory_after_pass = false;
                self.status = format!("{e:#}");
            }
        }
    }

    /// Take `path` (picked in the file browser) as the labels of the
    /// trajectory run, then ask for its output prefix.
    pub fn set_labels(&mut self, path: Option<&Path>) {
        match path {
            Some(p) => {
                self.trajectory.labels = Some(p.to_string_lossy().into_owned());
                self.ask_trajectory_out();
            }
            None => self.status = "no labels file picked".into(),
        }
    }

    /// A trajectory run wrote `manifest`: show its figures and verdicts.
    fn trajectory_done(&mut self, manifest: &Path, secs: u64) {
        let picker = self.picker().clone();
        self.figures = match super::figure_pane::FigurePane::load(&[manifest], &picker) {
            Ok(f) => f,
            Err(e) => {
                self.push_log(format!("[WARN] trajectory figures: {e:#}"));
                None
            }
        };
        let edges = match &self.figures {
            Some(f) => f.data.edges.clone(),
            None => super::order::run_edges(&[manifest]),
        };
        if let TreeMode::Order(v) = &mut self.tree_mode {
            v.edges = edges;
        }
        if let Some(f) = &mut self.figures {
            f.shown = true;
        }
        if matches!(self.tree_mode, TreeMode::Order(_)) {
            self.reload_order();
        }
        self.status = format!("trajectory done in {secs}s: {}", manifest.display());
    }

    /// The labels on screen with their cell counts: the round's labels, else
    /// the panel's types.
    fn order_types(&self) -> Vec<(String, usize)> {
        match &self.round {
            Some(r) => r
                .summary(&self.edits)
                .into_iter()
                .filter(|(l, _)| l != UNASSIGNED_LABEL)
                .collect(),
            // No round yet (`lupin trajectory`): the trajectory run's
            // types with their cells, else the panel's types.
            None => match &self.figures {
                Some(f) => {
                    let mut n: std::collections::BTreeMap<String, usize> = Default::default();
                    for t in &f.data.types {
                        *n.entry(t.to_string()).or_default() += 1;
                    }
                    n.into_iter()
                        .filter(|(l, _)| l != UNASSIGNED_LABEL)
                        .collect()
                }
                None => self.original.keys().map(|k| (k.clone(), 0)).collect(),
            },
        }
    }

    /// (Re)build the order view from the labels on screen and the data files.
    /// The trajectory run's edges are kept from the view on screen (or the
    /// figures), so only the precedence files are read again.
    fn reload_order(&mut self) {
        let (sel, edges) = match std::mem::replace(&mut self.tree_mode, TreeMode::Panel) {
            TreeMode::Order(v) => (v.sel, v.edges),
            _ => (
                0,
                match &self.figures {
                    Some(f) => f.data.edges.clone(),
                    None => {
                        let manifests = self.manifests();
                        let refs: Vec<&Path> = manifests.iter().map(PathBuf::as_path).collect();
                        super::order::run_edges(&refs)
                    }
                },
            ),
        };
        let mut v = super::order::OrderView::load(
            self.order_types(),
            self.cl.as_ref(),
            &self.data_search,
            edges,
        );
        v.sel = sel.min(v.types.len().saturating_sub(1));
        if let Some(e) = &v.error {
            self.status = format!("the prior is not a DAG: {e}");
        } else if !v.problems.is_empty() {
            self.status = format!("unreadable precedence file: {}", v.problems.join("; "));
        }
        let panel: Vec<(String, String)> = v
            .types
            .iter()
            .map(|(t, _)| (String::new(), t.clone()))
            .collect();
        self.order_tree = Some(crate::manifest::ontology::panel_tree_on(
            self.cl.as_ref(),
            &panel,
        ));
        self.tree_mode = TreeMode::Order(v);
        self.sync_tree_to_order();
    }

    /// Keys of the ontology beside the precedence table: the order view's
    /// types on the Cell Ontology, or the full ontology (`o`). Moving onto a
    /// type selects it in the table too.
    fn order_tree_key(&mut self, code: KeyCode) -> bool {
        if let (Some(v), Some(cl)) = (&mut self.order_cl, &self.cl) {
            if step(&mut v.sel, v.rows.len(), code) {
                return true;
            }
            match code {
                KeyCode::Right => v.enter(cl),
                KeyCode::Left => v.up(cl),
                KeyCode::Char('/') => v.typing = Some(String::new()),
                KeyCode::Char('d') => {
                    if !v.toggle_scope(cl) {
                        self.status = "no type of the data sits on a Cell Ontology term".into();
                    }
                }
                KeyCode::Char('o') => self.order_cl = None,
                KeyCode::Char('t') => self.toggle_order(false),
                KeyCode::Esc => self.focus = Focus::Clusters,
                _ => return false,
            }
            return true;
        }
        let Some(tree) = &mut self.order_tree else {
            return false;
        };
        let visible = tree.visible();
        if step(&mut self.order_tree_sel, visible.len(), code) {
            self.sync_order_to_tree();
            return true;
        }
        let at = visible.get(self.order_tree_sel).copied();
        match (code, at) {
            (KeyCode::Left, Some(i)) => {
                if tree.nodes[i].children.is_empty() || tree.is_folded(i) {
                    // Up to the parent, like a file tree.
                    if let Some(p) = tree.nodes[i].parent {
                        if let Some(k) = tree.visible().iter().position(|&j| j == p) {
                            self.order_tree_sel = k;
                        }
                    }
                } else {
                    tree.fold(i, true);
                }
                self.sync_order_to_tree();
            }
            (KeyCode::Right, Some(i)) => tree.fold(i, false),
            (KeyCode::Char(' '), Some(_)) => {
                // Only the order view's types can be ordered.
                if self.sync_order_to_tree() {
                    self.toggle_tree_mark();
                } else {
                    self.status = "not one of the run's types: only they can be ordered".into();
                }
            }
            (KeyCode::Char('o'), at) => {
                let Some(cl) = &self.cl else {
                    self.status = "no Cell Ontology at hand (see `lupin data where`)".into();
                    return true;
                };
                let focus = at
                    .and_then(|i| tree.nodes[i].cl_id.clone())
                    .filter(|id| cl.has(id))
                    .unwrap_or_else(|| ROOT_TERM.to_string());
                let data = self.ontology_data();
                let Some(cl) = &self.cl else {
                    return true;
                };
                self.order_cl = Some(super::ontology::OntologyView::in_data(cl, data, &focus));
            }
            (KeyCode::Char('t'), _) => self.toggle_order(false),
            (KeyCode::Esc, _) => self.focus = Focus::Clusters,
            _ => return false,
        }
        true
    }

    /// Select in the precedence table the type selected in the ontology
    /// beside it; `false` when that node is not one of the table's types.
    fn sync_order_to_tree(&mut self) -> bool {
        let (Some(tree), TreeMode::Order(v)) = (&self.order_tree, &mut self.tree_mode) else {
            return false;
        };
        let Some(&i) = tree.visible().get(self.order_tree_sel) else {
            return false;
        };
        let l = label_key(tree.label(i));
        match v.types.iter().position(|(t, _)| label_key(t) == l) {
            Some(k) => {
                v.sel = k;
                true
            }
            None => false,
        }
    }

    /// Select in the ontology beside the table the type selected in the
    /// precedence table, unfolding the way to it.
    fn sync_tree_to_order(&mut self) {
        let (Some(tree), TreeMode::Order(v)) = (&mut self.order_tree, &self.tree_mode) else {
            return;
        };
        let Some(i) = v.selected().and_then(|t| tree.node_of(t)) else {
            return;
        };
        tree.reveal(i);
        if let Some(k) = tree.visible().iter().position(|&j| j == i) {
            self.order_tree_sel = k;
        }
    }

    /// Keys of the order view: move, `>` / `-` on two marked types (the
    /// figure pane's keys come first, in `pane_key`).
    fn order_key(&mut self, code: KeyCode) -> bool {
        let shown = self.figures_shown();
        let TreeMode::Order(v) = &mut self.tree_mode else {
            return false;
        };
        // The table is hidden behind the figures: its keys wait for `V`.
        if !shown && step(&mut v.sel, v.types.len(), code) {
            return true;
        }
        let relation = match code {
            KeyCode::Char('r') => {
                self.ask_trajectory_out();
                return true;
            }
            KeyCode::Char('>') => Relation::Precedes,
            KeyCode::Char('-') => Relation::Unrelated,
            KeyCode::Esc => {
                self.focus = Focus::Clusters;
                return true;
            }
            _ => return false,
        };
        if shown {
            self.status = "V shows the order table to mark types".into();
            return true;
        }
        if self.tree_marked.len() != 2 {
            self.status = "mark exactly two types with space, in order, then > or -".into();
            return true;
        }
        let (from, to) = (self.tree_marked[0].clone(), self.tree_marked[1].clone());
        if let Some(why) = v.refuse(&from, &to, relation) {
            self.status = why;
            return true;
        }
        let title = match relation {
            Relation::Precedes => format!(" why does {from} precede {to}? "),
            Relation::Unrelated => format!(" why are {from} and {to} unrelated? "),
        };
        self.prompt = Some(Prompt {
            title,
            text: String::new(),
            pending: Pending::Precedence { from, to, relation },
        });
        true
    }

    /// Keys of the figure pane and its exports strip; `false` when the key is
    /// not one of them (or no figures are loaded).
    fn figure_key(&mut self, code: KeyCode) -> bool {
        let Some(v) = &mut self.figures else {
            if matches!(code, KeyCode::Char('v' | 'p')) {
                self.status = "no trajectory outputs for the manifests on screen: run `lupin trajectory` first".into();
                return true;
            }
            return false;
        };
        if code != KeyCode::Char('D') {
            v.exports.disarm();
        }
        // The thumbnail grid is modal: arrows choose, Enter opens, p exports
        // the tile chosen, and only the app's own keys (quit, run, stop,
        // save, export the round, panes) go by.
        if let Some(g) = &mut v.grid {
            if g.step(code) {
                return true;
            }
            match code {
                KeyCode::Enter => {
                    v.open_tile();
                    self.status = v.title(v.current());
                }
                KeyCode::Esc | KeyCode::Char('w') => v.grid = None,
                KeyCode::Char('V') => {
                    v.grid = None;
                    v.shown = false;
                }
                KeyCode::Char('t') => self.status = v.style.cycle_labels(),
                KeyCode::Char('c') => self.status = v.style.cycle_colouring(&v.data.colourings()),
                KeyCode::Char('p') => {
                    self.status = v
                        .export()
                        .unwrap_or_else(|e| format!("export failed: {e:#}"));
                }
                KeyCode::Char('f') => v.exports.open = !v.exports.open,
                KeyCode::Char(c) if !"qrxse".contains(c) => {
                    self.status = GRID_KEYS.into();
                }
                _ => return false,
            }
            return true;
        }
        let x = &mut v.exports;
        // The strip's keys act only while it is on screen; the arrows it
        // does not take pan.
        let strip = v.shown && x.open;
        if strip {
            let rows = x.rows();
            if step(&mut x.sel, rows, code) {
                return true;
            }
        }
        match code {
            KeyCode::Char('v') => v.next_panel(),
            KeyCode::Char('V') => {
                v.shown = false;
                v.grid = None;
            }
            KeyCode::Char(',') if v.shown => v.step_pair(false),
            KeyCode::Char('.') if v.shown => v.step_pair(true),
            // senna view's keys: labels, colouring, layout.
            KeyCode::Char('t') if v.shown => self.status = v.style.cycle_labels(),
            KeyCode::Char('c') if v.shown => {
                self.status = v.style.cycle_colouring(&v.data.colourings());
            }
            // The table is hidden behind the figures: its marks wait for `V`.
            KeyCode::Char(' ') if v.shown => {
                self.status = "V shows the order table to mark types".into();
            }
            KeyCode::Char('m') if v.shown => self.status = v.next_layout(),
            // senna view's navigation: the grid, zoom and pan.
            KeyCode::Char('w') if v.shown => {
                v.open_grid();
                self.status = GRID_KEYS.into();
            }
            KeyCode::Char('+' | '=') if v.shown => self.status = v.zoom(true),
            KeyCode::Char('-' | '_') if v.shown => self.status = v.zoom(false),
            KeyCode::Char('0') if v.shown => self.status = v.reset_view(),
            KeyCode::Left | KeyCode::Right | KeyCode::Up | KeyCode::Down if v.shown => {
                let (dx, dy) = match code {
                    KeyCode::Left => (-1.0, 0.0),
                    KeyCode::Right => (1.0, 0.0),
                    KeyCode::Up => (0.0, 1.0),
                    _ => (0.0, -1.0),
                };
                let note = v.pan(dx, dy);
                if !note.is_empty() {
                    self.status = note;
                }
            }
            KeyCode::Char('p') if v.shown => {
                self.status = v
                    .export()
                    .unwrap_or_else(|e| format!("export failed: {e:#}"));
            }
            KeyCode::Char('f') => x.open = !x.open,
            KeyCode::Char('R') => {
                x.refresh();
                self.status = format!(
                    "{} export(s) checked, {} not listed",
                    x.gallery.entries.len(),
                    x.not_listed.len()
                );
            }
            KeyCode::Enter if strip => {
                self.status = x
                    .adopt(&v.data.manifest)
                    .unwrap_or_else(|e| format!("{e:#}"));
            }
            KeyCode::Char('d') if strip => {
                self.status = x.remove(false).unwrap_or_else(|e| format!("{e:#}"));
            }
            KeyCode::Char('D') if strip => {
                self.status = match x.arm_delete() {
                    Some(name) => format!("D again deletes {name} (svg and pdf)"),
                    None => x.remove(true).unwrap_or_else(|e| format!("{e:#}")),
                };
                return true;
            }
            KeyCode::Char('M') if strip => {
                if let Some(e) = x.gallery.entries.get(x.sel) {
                    self.prompt = Some(Prompt {
                        title: " move the export to (base name, no extension): ".into(),
                        text: e.path.with_extension("").to_string_lossy().into_owned(),
                        pending: Pending::Relocate {
                            pdf: e.path.clone(),
                        },
                    });
                }
            }
            _ => return false,
        }
        true
    }

    /// Switch the tree pane between the panel's tree and the Cell Ontology,
    /// the ontology opening on the selected node's term, else the cluster's.
    fn toggle_ontology(&mut self) {
        if matches!(self.tree_mode, TreeMode::Ontology(_)) {
            self.tree_mode = TreeMode::Panel;
            return;
        }
        let Some(cl) = &self.cl else {
            self.status = "no Cell Ontology at hand (see `lupin data where`)".into();
            return;
        };
        let from_node = self
            .tree
            .visible()
            .get(self.tree_sel)
            .and_then(|&i| self.tree.nodes[i].cl_id.clone());
        let from_cluster = || super::ontology::term_of(cl, &self.tree, &self.label_or_top()?);
        let focus = from_node
            .or_else(from_cluster)
            .filter(|id| cl.has(id))
            .unwrap_or_else(|| ROOT_TERM.to_string());
        let data = self.ontology_data();
        self.tree_mode =
            TreeMode::Ontology(super::ontology::OntologyView::in_data(cl, data, &focus));
    }

    /// The types in the data by the Cell Ontology term each sits on, with
    /// their cells: the order view's types, else the round's labels, else
    /// the marker panel's types.
    fn ontology_data(&self) -> std::collections::BTreeMap<String, Vec<(String, usize)>> {
        let Some(cl) = &self.cl else {
            return Default::default();
        };
        let types: Vec<(String, usize)> = match (&self.tree_mode, &self.round) {
            (TreeMode::Order(v), _) => v.types.clone(),
            (_, Some(r)) => r
                .summary(&self.edits)
                .into_iter()
                .filter(|(l, _)| l != UNASSIGNED_LABEL)
                .collect(),
            _ => self
                .tree
                .typed_terms()
                .map(|(l, _)| (l.to_string(), 0))
                .collect(),
        };
        let mut data: std::collections::BTreeMap<String, Vec<(String, usize)>> = Default::default();
        for (t, n) in types {
            if let Some(id) = super::ontology::term_of(cl, &self.tree, &t).filter(|id| cl.has(id)) {
                data.entry(id).or_default().push((t, n));
            }
        }
        data
    }

    fn ontology_key(&mut self, code: KeyCode) -> bool {
        let (TreeMode::Ontology(v), Some(cl)) = (&mut self.tree_mode, &self.cl) else {
            return false;
        };
        step(&mut v.sel, v.rows.len(), code);
        match code {
            KeyCode::Right => v.enter(cl),
            KeyCode::Left => v.up(cl),
            KeyCode::Char('/') => v.typing = Some(String::new()),
            KeyCode::Char('d') => {
                if !v.toggle_scope(cl) {
                    self.status = "no type of the data sits on a Cell Ontology term".into();
                }
            }
            KeyCode::Enter => {
                if let Some(id) = v.selected().map(|r| r.id.clone()) {
                    self.pick_term(&id);
                }
            }
            KeyCode::Esc => self.focus = Focus::Clusters,
            _ => return false,
        }
        true
    }

    /// Label the selected cluster by CL term `id`.
    fn pick_term(&mut self, id: &str) {
        let Some(cl) = &self.cl else { return };
        let label = super::ontology::term_label(cl, &self.tree, id);
        let term_of = |l: &str| super::ontology::term_of(cl, &self.tree, l);
        let pooled =
            self.pooled(|t| term_of(t).is_some_and(|t| cl.ancestors_or_self(&t).contains(id)));
        let name = cl.name(id).unwrap_or(id);
        let reason = if pooled.is_empty() {
            format!("{name} ({id}) from the Cell Ontology")
        } else {
            format!(
                "{name} ({id}) from the Cell Ontology, over {}",
                pooled.join(" + ")
            )
        };
        // The top candidate, if the ontology has no term for it: offer to
        // remember that this is the term it means.
        let remember = self
            .selected()
            .and_then(|c| c.candidates.first())
            .filter(|t| term_of(&t.label).is_none())
            .map(|t| (t.label.clone(), id.to_string()));
        self.ask_label(label, reason, remember);
    }

    /// Select the next flagged cluster without an edit, wrapping around.
    fn next_flagged(&mut self) {
        let Some(r) = &self.round else { return };
        let n = r.clusters.len();
        let next = (1..=n)
            .map(|d| (self.cluster_sel + d) % n)
            .find(|&i| r.clusters[i].flagged() && !self.decided(r.clusters[i].id));
        match next {
            Some(i) => self.select_cluster(i),
            None => self.status = "every flagged cluster is decided".into(),
        }
    }

    /// Focus the tree on the selected cluster's label (or its top candidate).
    fn jump_to_tree(&mut self) {
        if let Some(i) = self.label_or_top().and_then(|l| self.tree.node_of(&l)) {
            self.tree.reveal(i);
            self.select_node(i);
        }
        self.focus = Focus::Tree;
    }
}

/// A TUI that ends on an error or a panic stops its children too, so no
/// pass or trajectory goes on writing behind it. A save is let finish, as
/// stopping it could leave a round half-written.
impl Drop for App {
    fn drop(&mut self) {
        if let Some((mut c, _, job)) = self.child.take() {
            if !matches!(job, Job::Save(_)) {
                let _ = c.kill();
            }
            let _ = c.wait();
        }
        self.stop_rescoring();
    }
}

/// What the thumbnail grid's keys do.
const GRID_KEYS: &str =
    "arrows choose · Enter opens · p exports it · t c restyle · V table · esc or w closes";

/// The labels file a trajectory run on `loaded` recorded in its settings,
/// resolved against the manifest's directory, when it is still there.
pub(super) fn recorded_labels(loaded: &crate::manifest::run::Loaded) -> Option<String> {
    let settings = loaded.manifest.trajectory.settings.as_ref()?;
    let rel = settings.get("labels")?.as_str()?;
    let path = crate::manifest::run::resolve(&loaded.dir, rel);
    Path::new(&path).is_file().then_some(path)
}

/// `X.decisions.jsonl` for round `X.senna.json`: what a save hands relabel.
fn decisions_file(round: &Path) -> PathBuf {
    let stem = crate::manifest::run::derive_out_prefix(&round.to_string_lossy());
    PathBuf::from(format!("{stem}.decisions.jsonl"))
}

/// Joins the parts of a mixed label (`Basophils+Mast_cells`).
pub const MIX: &str = "+";

/// The Cell Ontology's root, where browsing starts without a better place.
const ROOT_TERM: &str = "CL:0000000";

/// Append `label → id` to the alias table `file`, creating it (and its
/// directory) with a header the first time.
fn remember_alias(file: &Path, label: &str, id: &str, note: &str) -> anyhow::Result<()> {
    crate::manifest::data_files::append_row(
        file,
        "Cell-type labels mapped to Cell Ontology terms, layered over lupin's own\n\
         (see `lupin data where`). Columns: label, CL id, note.",
        &[label, id, note],
    )
}

/// A candidate's evidence for a reason: share, and NES and q when known.
fn evidence(c: &super::round::Candidate) -> String {
    let mut s = format!("share {:.2}", c.share);
    if let Some(nes) = c.nes {
        s += &format!(", NES {nes:.2}");
    }
    if let Some(q) = c.q {
        s += &format!(", q {q:.3}");
    }
    s
}

pub(super) fn file_name(p: &Path) -> String {
    p.file_name().map_or_else(
        || p.display().to_string(),
        |n| n.to_string_lossy().into_owned(),
    )
}

/// Rows a page key moves.
const PAGE: usize = 10;

/// Move selection `sel` in a list of `n` rows by `code`: a row (↑↓), a page
/// (PgUp/PgDn) or to an end (Home/End). Other keys leave it; `false` for
/// them.
pub(super) fn step(sel: &mut usize, n: usize, code: KeyCode) -> bool {
    let last = n.saturating_sub(1);
    *sel = match code {
        KeyCode::Up => sel.saturating_sub(1),
        KeyCode::Down => (*sel + 1).min(last),
        KeyCode::PageUp => sel.saturating_sub(PAGE),
        KeyCode::PageDown => (*sel + PAGE).min(last),
        KeyCode::Home => 0,
        KeyCode::End => last,
        _ => return false,
    };
    true
}

#[cfg(test)]
#[path = "tests/app.rs"]
mod tests;
