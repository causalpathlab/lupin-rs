//! The screen's state and what each key does to it.
//!
//! Work goes cluster by cluster: pick a cluster, weigh its candidates and
//! genes, and give it a label (a candidate, any node of the ontology tree,
//! or none), each with a reason. Genes can be added to a cell type's
//! markers. Saving writes the edits as the next round.

use super::round::{decisions, Edit, RoundView};
use super::runner;
use crate::annotate::markers::label_key;
use crate::annotate::panel_tree::PanelTree;
use crate::annotate::rounds::ClusterId;
use crate::annotate_cmd::{AnnotateCliArgs, AnnotateMethod};
use crate::manifest::rounds::chain_rounds;
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
}

/// The settings the screen edits, in the order it lists them.
#[derive(Clone, Copy)]
pub enum Setting {
    Method,
    Knn,
    Resolution,
    NumClusters,
    NumPerm,
}

pub const SETTINGS: [Setting; 5] = [
    Setting::Method,
    Setting::Knn,
    Setting::Resolution,
    Setting::NumClusters,
    Setting::NumPerm,
];

impl Setting {
    pub fn name(self) -> &'static str {
        match self {
            Self::Method => "method",
            Self::Knn => "knn",
            Self::Resolution => "resolution",
            Self::NumClusters => "clusters (k)",
            Self::NumPerm => "permutations",
        }
    }

    pub fn value(self, a: &AnnotateCliArgs) -> String {
        match self {
            Self::Method => match a.method {
                AnnotateMethod::Auto => "auto".into(),
                AnnotateMethod::Enrichment => "enrichment".into(),
                AnnotateMethod::Projection => "projection".into(),
            },
            Self::Knn => a.knn.map_or("default".into(), |k| k.to_string()),
            Self::Resolution => format!("{:.2}", a.resolution),
            Self::NumClusters => a.num_clusters.map_or("auto".into(), |k| k.to_string()),
            Self::NumPerm => a.num_perm.to_string(),
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
        }
    }
}

/// What the tree pane shows.
pub enum TreeMode {
    /// The panel's types on the tree lupin built.
    Panel,
    /// The Cell Ontology around a term.
    Ontology(super::ontology::OntologyView),
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
}

/// A one-line prompt for a decision's reason, prefilled.
pub struct Prompt {
    pub title: String,
    pub text: String,
    pub pending: Pending,
}

/// What a child process is doing.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Job {
    Pass,
    /// Writing this many edits as the next round.
    Save(usize),
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
}

pub struct App {
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
            cl: None,
            panel_ancestry: Vec::new(),
            data_search: crate::manifest::data_files::SearchPath::new(None),
            original: std::collections::BTreeMap::new(),
            round: None,
            edits: Vec::new(),
            focus: Focus::Clusters,
            cluster_sel: 0,
            gene_sel: 0,
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
            self.push_log(line);
        }
        for line in super::drain_own_log() {
            self.push_log(line);
        }
        self.keep_scores_current();
        let Some((child, started, job)) = &mut self.child else {
            return;
        };
        let (secs, job) = (started.elapsed().as_secs(), *job);
        match child.try_wait() {
            Ok(None) => {
                let what = match job {
                    Job::Pass => "pass",
                    Job::Save(_) => "saving",
                };
                self.status = format!("{what}… {secs}s  (x: stop)");
            }
            Ok(Some(st)) if st.success() => {
                self.child = None;
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
                    Job::Pass => Vec::new(),
                };
                self.open(&latest);
                let kept = later.len();
                self.edits = later;
                let done = match job {
                    Job::Pass => format!("pass done in {secs}s"),
                    Job::Save(n) if kept > 0 => format!(
                        "saved {n} edit(s), {kept} made since still unsaved; {}",
                        self.export()
                    ),
                    Job::Save(n) => format!("saved {n} edit(s); {}", self.export()),
                };
                self.status = format!("{done}. {}", self.status);
            }
            Ok(Some(st)) => {
                self.child = None;
                let what = match job {
                    Job::Pass => "pass",
                    Job::Save(_) => "save",
                };
                self.status = format!("{what} failed ({st}); see the log");
            }
            Err(e) => {
                self.child = None;
                self.status = format!("lost the pass: {e}");
            }
        }
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
                            None => self.status = "rescoring failed; see the log".into(),
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
        let file = std::env::temp_dir().join(format!(
            "lupin-rescore-{}-{}.jsonl",
            std::process::id(),
            self.edits.len()
        ));
        let started = (|| -> anyhow::Result<_> {
            let lines: Vec<String> = decisions(&want, round)
                .iter()
                .map(serde_json::to_string)
                .collect::<Result<_, _>>()?;
            std::fs::write(&file, lines.join("\n") + "\n")?;
            runner::spawn_preview(&round.manifest, &file, self.log_tx.clone())
        })();
        match started {
            Ok((child, out)) => {
                self.status = "rescoring the edited types…".into();
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
        let Some(round) = &mut self.round else { return };
        let scores = serde_json::from_str::<serde_json::Value>(json)
            .ok()
            .as_ref()
            .and_then(super::round::parse_scores);
        let Some(scores) = scores else {
            self.status =
                "this round cannot be rescored (no cached statistics): save to rescore".into();
            return;
        };
        let old = round.swap_scores(scores);
        // The round's own are what is kept, across several rescorings.
        if self.recorded.is_none() {
            self.recorded = Some(old);
        }
        self.status = "rescored with your marker edits".into();
    }

    fn stop_rescoring(&mut self) {
        if let Some(r) = self.rescoring.take() {
            r.stop();
        }
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
        self.marked.clear();
        self.tree_marked.clear();
    }

    /// The label of the tree pane's selected node: a panel node's, or a CL
    /// term's in the ontology view.
    fn tree_selected_label(&self) -> Option<String> {
        match (&self.tree_mode, &self.cl) {
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
        self.prompt = Some(Prompt {
            title: format!(" K{} keeps {label}: why? ", c.id),
            text: reason,
            pending: Pending::Keep { cluster: c.id },
        });
    }

    fn start(&mut self) {
        if self.child.is_some() {
            return;
        }
        let (_, later) = chain_rounds(&self.target);
        if (!self.edits.is_empty() || !later.is_empty()) && self.armed != Some(Armed::Run) {
            self.armed = Some(Armed::Run);
            self.status = format!(
                "a new pass drops {} unsaved edit(s) and sets aside {} saved round(s): r again to go on",
                self.edits.len(),
                later.len()
            );
            return;
        }
        self.push_log(format!(
            "── pass: lupin annotate {}",
            self.args.to_argv().join(" ")
        ));
        match runner::spawn_pass(&self.args, self.log_tx.clone()) {
            Ok(c) => {
                self.child = Some((c, Instant::now(), Job::Pass));
                self.status = "running…".into();
            }
            Err(e) => self.status = format!("{e:#}"),
        }
    }

    /// Stop a running pass. A save is not stopped midway, which could leave
    /// a round half-written: it finishes on its own.
    fn stop(&mut self) {
        match self.child.as_ref().map(|c| c.2) {
            Some(Job::Pass) => {
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
            let lines: Vec<String> = decisions(&self.edits, round)
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
        self.ask_label(pick.label, reason, None);
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

    /// Drop the selected cluster's edits.
    fn undo(&mut self) {
        let Some(id) = self.selected().map(|c| c.id) else {
            return;
        };
        self.edits.retain(|e| e.cluster() != Some(id));
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
        match k.code {
            KeyCode::Esc => self.prompt = None,
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
                    self.status = "a decision needs a reason".into();
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
                                text: "picked in lupin annotate --tui".into(),
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
                    Pending::Remember { label, id, file } => {
                        self.status = match remember_alias(&file, &label, &id, &reason) {
                            Ok(()) => format!("remembered {label} → {id} in {}", file.display()),
                            Err(e) => format!("could not remember the alias: {e:#}"),
                        };
                    }
                    Pending::Markers { label, genes, add } => {
                        let moved = if add {
                            super::round::take_added(&mut self.edits, &genes, &label)
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
            // A running save is left to finish; a pass is stopped.
            if self.pass_running() {
                self.stop();
            }
            self.stop_rescoring();
            self.quit = true;
            return;
        }
        if self.prompt.is_some() {
            return self.prompt_key(k);
        }
        if let (TreeMode::Ontology(v), Some(cl)) = (&mut self.tree_mode, &self.cl) {
            if let Some(text) = &mut v.typing {
                match k.code {
                    KeyCode::Esc => v.typing = None,
                    KeyCode::Enter => {
                        let q = std::mem::take(text);
                        v.typing = None;
                        v.search(cl, &q);
                        self.status = format!("{} term(s) match {q:?} (←: back)", v.rows.len());
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
            return self.settings_key(k.code, armed);
        }
        match k.code {
            KeyCode::Char('q') => {
                let busy = self.child.as_ref().map(|c| c.2);
                let unsaved = match self.edits.len() {
                    0 => String::new(),
                    n => format!(" and {n} unsaved edit(s)"),
                };
                if matches!(busy, Some(Job::Save(_))) {
                    self.status = "saving… q once it is done".into();
                } else if armed != Some(Armed::Quit) && busy.is_some() {
                    self.armed = Some(Armed::Quit);
                    self.status = format!("a pass is running{unsaved}: q again stops it and quits");
                } else if armed == Some(Armed::Quit) || self.edits.is_empty() {
                    self.stop();
                    self.stop_rescoring();
                    self.quit = true;
                } else {
                    self.armed = Some(Armed::Quit);
                    self.status = format!(
                        "{} unsaved edit(s): q again to quit, s to save",
                        self.edits.len()
                    );
                }
            }
            KeyCode::Char('r') => self.settings_open = true,
            // Stops a running pass or save; otherwise it is the pane's (hide, in genes).
            KeyCode::Char('x') if self.child.is_some() => self.stop(),
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
            code => self.pane_key(code),
        }
    }

    /// Keys in the clustering and pass settings popup.
    fn settings_key(&mut self, code: KeyCode, armed: Option<Armed>) {
        match code {
            KeyCode::Esc | KeyCode::Char('r' | 'q') => self.settings_open = false,
            KeyCode::Up => self.setting = self.setting.saturating_sub(1),
            KeyCode::Down if self.setting + 1 < SETTINGS.len() => self.setting += 1,
            KeyCode::Left | KeyCode::Right => {
                SETTINGS[self.setting].adjust(&mut self.args, code == KeyCode::Right);
                self.stale = self.round.is_some();
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
        use Focus::{Clusters, Genes, Tree};
        let order = [Clusters, Tree, Genes];
        let i = order.iter().position(|f| *f == self.focus).unwrap_or(0);
        let n = order.len();
        order[if forward {
            (i + 1) % n
        } else {
            (i + n - 1) % n
        }]
    }

    fn pane_key(&mut self, code: KeyCode) {
        let up = code == KeyCode::Up;
        let down = code == KeyCode::Down;
        let move_in = |sel: &mut usize, n: usize| {
            if up {
                *sel = sel.saturating_sub(1);
            } else if down && *sel + 1 < n {
                *sel += 1;
            }
        };
        match self.focus {
            Focus::Clusters => {
                let (before, mut sel) = (self.cluster_sel, self.cluster_sel);
                move_in(&mut sel, self.n_clusters());
                if sel != before {
                    self.select_cluster(sel);
                }
                match code {
                    KeyCode::Char(c @ '1'..='9') => self.pick_candidate(c as usize - '1' as usize),
                    KeyCode::Char('u') => self.unassign(),
                    KeyCode::Char('k') => self.keep(),
                    KeyCode::Backspace | KeyCode::Delete => self.undo(),
                    KeyCode::Char(']') => self.next_flagged(),
                    KeyCode::Enter => self.jump_to_tree(),
                    _ => {}
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
                    KeyCode::Char('m') => {
                        self.gene_view = match self.gene_view {
                            GeneView::Specific => GeneView::Markers,
                            GeneView::Markers => GeneView::Specific,
                        };
                        self.gene_sel = 0;
                        self.marked.clear();
                    }
                    KeyCode::Esc => self.focus = Focus::Clusters,
                    _ => {}
                }
            }
            Focus::Tree if code == KeyCode::Char('o') => self.toggle_ontology(),
            Focus::Tree if code == KeyCode::Char(' ') => self.toggle_tree_mark(),
            Focus::Tree if code == KeyCode::Char('+') => self.mix_marked(),
            Focus::Tree if matches!(self.tree_mode, TreeMode::Ontology(_)) => {
                self.ontology_key(code);
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
                    _ => {}
                }
            }
        }
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
        self.tree_mode = TreeMode::Ontology(super::ontology::OntologyView::at(cl, &focus));
    }

    fn ontology_key(&mut self, code: KeyCode) {
        let (TreeMode::Ontology(v), Some(cl)) = (&mut self.tree_mode, &self.cl) else {
            return;
        };
        match code {
            KeyCode::Up => v.sel = v.sel.saturating_sub(1),
            KeyCode::Down if v.sel + 1 < v.rows.len() => v.sel += 1,
            KeyCode::Right => v.enter(cl),
            KeyCode::Left => v.up(cl),
            KeyCode::Char('/') => v.typing = Some(String::new()),
            KeyCode::Enter => {
                if let Some(id) = v.selected().map(|r| r.id.clone()) {
                    self.pick_term(&id);
                }
            }
            KeyCode::Esc => self.focus = Focus::Clusters,
            _ => {}
        }
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
    crate::manifest::data_files::append_line(
        file,
        "Cell-type labels mapped to Cell Ontology terms, layered over lupin's own\n\
         (see `lupin data where`). Columns: label, CL id, note.",
        &format!("{label}\t{id}\t{}", note.replace(['\t', '\n'], " ")),
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

fn file_name(p: &Path) -> String {
    p.file_name().map_or_else(
        || p.display().to_string(),
        |n| n.to_string_lossy().into_owned(),
    )
}

#[cfg(test)]
#[path = "tests/app.rs"]
mod tests;
