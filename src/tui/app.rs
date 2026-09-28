//! The screen's state and what each key does to it.
//!
//! Work goes cluster by cluster: pick a cluster, weigh its candidates and
//! genes, and give it a label (a candidate, any node of the ontology tree,
//! or none), each with a reason. Genes can be added to a cell type's
//! markers. Saving writes the edits as the next round.

use super::round::{decisions, Edit, RoundView};
use super::runner;
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
    Settings,
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
    /// Genes marked for adding or dropping as markers.
    pub marked: Vec<String>,
    pub tree_sel: usize,
    pub setting: usize,
    pub prompt: Option<Prompt>,
    pub log: Vec<String>,
    log_tx: Sender<String>,
    log_rx: Receiver<String>,
    pub child: Option<(Child, Instant, Job)>,
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
            marked: Vec::new(),
            tree_sel: 0,
            setting: 0,
            prompt: None,
            log: Vec::new(),
            log_tx,
            log_rx,
            child: None,
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
                self.open(&latest);
                let done = match job {
                    Job::Pass => format!("pass done in {secs}s"),
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
                self.round = Some(r);
                self.edits.clear();
                self.marked.clear();
                self.stale = false;
                self.cluster_sel = self.cluster_sel.min(self.n_clusters().saturating_sub(1));
            }
            Err(e) => self.status = format!("could not read {}: {e:#}", manifest.display()),
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

    /// Whether cluster `id` has an unsaved label edit.
    pub fn edited(&self, id: ClusterId) -> bool {
        self.edits
            .iter()
            .any(|e| matches!(e, Edit::Label { cluster, .. } if *cluster == id))
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

    fn stop(&mut self) {
        if let Some((mut c, _, _)) = self.child.take() {
            let _ = c.kill();
            let _ = c.wait();
            self.status = "stopped".into();
        }
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
        match super::export::write(r, &self.tree, self.cl.as_ref(), &self.original) {
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
        let pooled =
            self.pooled(|t| under.contains(&crate::annotate::markers::label_key(t).as_str()));
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
        self.edits
            .retain(|e| !matches!(e, Edit::Label { cluster, .. } if *cluster == id));
    }

    /// The genes the genes pane lists: the cluster's specific genes, or
    /// its label's markers with the edits applied.
    pub fn listed_genes(&self) -> Vec<String> {
        match self.gene_view {
            GeneView::Specific => self
                .selected()
                .map(|c| c.genes.iter().map(|g| g.0.clone()).collect())
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

    /// Open a prompt to add (`add`) or drop the marked genes, else the
    /// selected one, as markers of the selected cluster's label. Genes
    /// already in (for adding) or not in (for dropping) are left out.
    fn ask_markers(&mut self, add: bool) {
        let (Some(c), Some(r)) = (self.selected(), &self.round) else {
            return;
        };
        let Some(label) = self.current_label() else {
            self.status = "label the cluster first: markers belong to a cell type".into();
            return;
        };
        let chosen = if self.marked.is_empty() {
            self.listed_genes()
                .get(self.gene_sel)
                .cloned()
                .into_iter()
                .collect()
        } else {
            self.marked.clone()
        };
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
                        if let (Some((alias, id)), Some(file)) =
                            (remember, self.data_search.amend_aliases())
                        {
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
                    Pending::Remember { label, id, file } => {
                        self.status = match remember_alias(&file, &label, &id, &reason) {
                            Ok(()) => format!("remembered {label} → {id} in {}", file.display()),
                            Err(e) => format!("could not remember the alias: {e:#}"),
                        };
                    }
                    Pending::Markers { label, genes, add } => {
                        self.status = format!(
                            "{} {} {label}'s markers on save (the next round is rescored)",
                            genes.len(),
                            if add { "added to" } else { "dropped from" }
                        );
                        self.edits.push(Edit::Markers {
                            label,
                            genes,
                            add,
                            reason,
                        });
                        self.marked.clear();
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
            self.stop();
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
        let armed = self.armed.take();
        match k.code {
            KeyCode::Char('q') => {
                let busy = self.child.as_ref().map(|c| c.2);
                if armed != Some(Armed::Quit) && busy.is_some() {
                    self.armed = Some(Armed::Quit);
                    self.status = match busy {
                        Some(Job::Save(_)) => "a save is running: q again stops it and loses it",
                        _ => "a pass is running: q again stops it",
                    }
                    .into();
                } else if armed == Some(Armed::Quit) || self.edits.is_empty() {
                    self.stop();
                    self.quit = true;
                } else {
                    self.armed = Some(Armed::Quit);
                    self.status = format!(
                        "{} unsaved edit(s): q again to quit, s to save",
                        self.edits.len()
                    );
                }
            }
            KeyCode::Char('r') => {
                self.armed = armed;
                self.start();
            }
            KeyCode::Char('x') => self.stop(),
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

    fn cycle(&self, forward: bool) -> Focus {
        use Focus::{Clusters, Genes, Settings, Tree};
        let order = [Clusters, Tree, Genes, Settings];
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
            Focus::Settings => {
                move_in(&mut self.setting, SETTINGS.len());
                let inc = code == KeyCode::Right;
                if inc || code == KeyCode::Left {
                    SETTINGS[self.setting].adjust(&mut self.args, inc);
                    self.stale = self.round.is_some();
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
            .find(|&i| r.clusters[i].flagged() && !self.edited(r.clusters[i].id));
        match next {
            Some(i) => self.select_cluster(i),
            None => self.status = "every flagged cluster has an edit".into(),
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

/// The Cell Ontology's root, where browsing starts without a better place.
const ROOT_TERM: &str = "CL:0000000";

/// Append `label → id` to the alias table `file`, creating it (and its
/// directory) with a header the first time.
fn remember_alias(file: &Path, label: &str, id: &str, note: &str) -> anyhow::Result<()> {
    use std::io::Write;
    if let Some(dir) = file.parent() {
        std::fs::create_dir_all(dir)?;
    }
    let new = !file.exists();
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(file)?;
    if new {
        writeln!(
            f,
            "# Cell-type labels mapped to Cell Ontology terms, layered over lupin's own\n\
             # (see `lupin data where`). Columns: label, CL id, note.\nlabel\tcl_id\tnote"
        )?;
    }
    writeln!(f, "{label}\t{id}\t{}", note.replace(['\t', '\n'], " "))?;
    Ok(())
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
