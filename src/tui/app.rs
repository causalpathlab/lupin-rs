//! The screen's state and what each key does to it.
//!
//! Work goes cluster by cluster: pick a cluster, weigh its candidates and
//! genes, and give it a label (a candidate, any node of the ontology tree,
//! or none), each with a reason. Genes can be added to a cell type's
//! markers. Saving writes the edits as the next round.

use super::figure_pane::FigureKey;
use super::menu::{numbered, Action, Menu, Outcome};
use super::ontology::ViewKey;
use super::round::{decisions, Edit, RoundView};
use super::runner::{self, Held};
use crate::annotate::args::ProjectionMethod;
use crate::annotate::celltype_tree::ClTerms;
use crate::annotate::markers::label_key;
use crate::annotate::panel_tree::PanelTree;
use crate::annotate::rounds::ClusterId;
use crate::annotate_cmd::{AnnotateCliArgs, AnnotateMethod};
use crate::manifest::rounds::chain_rounds;
use crate::manifest::run::annotated_path;
use crate::trajectory::prior::Relation;
use anyhow::Context as _;
use enrichment::UNASSIGNED_LABEL;
use legume_tui::graphics::Repick;
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
use std::path::{Path, PathBuf};
use std::process::ExitStatus;
use std::sync::mpsc::{channel, Receiver, Sender};
use std::time::{Duration, Instant};

/// Log lines kept for the log pane.
const LOG_KEEP: usize = 1000;
/// A cluster whose top candidate has less of its evidence than this is flagged.
pub const CONTESTED: f32 = 0.5;
/// How long a key's reply holds the status line while a job runs, before
/// the job's progress takes it back; and how long a key asked again (`x`,
/// `q`) waits for its second press.
pub const KEY_REPLY_HOLD: Duration = Duration::from_secs(4);

/// The annotation form's state: the selected row, and why it did not run
/// or what running it again will do (`confirm`: the user has been told the
/// pass replaces something, and the next run goes on).
#[derive(Debug, Default)]
pub struct Form {
    pub row: usize,
    pub note: Option<String>,
    pub confirm: bool,
    /// Where the pass gets its cluster expression, when not the counts.
    pub expression: Option<String>,
}

impl Form {
    /// The setting on the cursor's row of `args`'s form ([`settings`]).
    #[must_use]
    pub fn setting(&self, args: &AnnotateCliArgs) -> Setting {
        let rows = settings(args);
        rows[self.row.min(rows.len() - 1)]
    }
}

/// What follows a pass.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AfterPass {
    #[default]
    Nothing,
    /// Run the trajectory on the new labels.
    Run,
    /// Ask whether to (a pass started from the order view).
    Offer,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Focus {
    Clusters,
    Genes,
    Tree,
    /// The order view's precedence table (the tree pane then shows the
    /// ontology beside it).
    Order,
    /// The selected cluster's GO terms, when the round scored them.
    Go,
    /// The selected cluster's candidates, which scroll with their selection.
    Candidates,
}

/// The rows of the annotation form, in the order it lists them.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Setting {
    /// The marker panel; Enter picks it in the file browser.
    Markers,
    /// The output prefix of passes; Enter edits it.
    Output,
    Method,
    /// Projection's way of calling: the embedding marker score or the ORA.
    Projection,
    Knn,
    Resolution,
    NumClusters,
    NumPerm,
    Go,
    /// The enrichment's SuSiE stage, which makes the marker call.
    Susie,
    /// Start the pass (Enter here, or Shift+Enter anywhere).
    Run,
}

/// The form's rows for `a`'s method: the enrichment's own (GO, SuSiE) only
/// when it may run, projection's way of calling only when projection may,
/// and permutations when a null permutes (the enrichment's, the ORA's).
#[must_use]
pub fn settings(a: &AnnotateCliArgs) -> Vec<Setting> {
    use Setting::{
        Go, Knn, Markers, Method, NumClusters, NumPerm, Output, Projection, Resolution, Run, Susie,
    };
    let enrichment = a.method != AnnotateMethod::Projection;
    let projection = a.method != AnnotateMethod::Enrichment;
    let permutes = enrichment || a.projection != Some(ProjectionMethod::Score);
    [
        (Markers, true),
        (Output, true),
        (Method, true),
        (Projection, projection),
        (Knn, true),
        (Resolution, true),
        (NumClusters, true),
        (NumPerm, permutes),
        (Go, enrichment),
        (Susie, enrichment),
        (Run, true),
    ]
    .into_iter()
    .filter_map(|(s, shown)| shown.then_some(s))
    .collect()
}

/// Where `s` sits in `a`'s form ([`settings`]).
pub fn setting_row(a: &AnnotateCliArgs, s: Setting) -> usize {
    settings(a).iter().position(|&t| t == s).unwrap_or(0)
}

impl Setting {
    /// The Leiden knobs, which a source naming its own clusters leaves
    /// fixed.
    #[must_use]
    pub fn leiden(self) -> bool {
        matches!(self, Self::Knn | Self::Resolution | Self::NumClusters)
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Markers => "marker panel",
            Self::Run => "▶ run",
            Self::Output => "output",
            Self::Method => "method",
            Self::Projection => "projection",
            Self::Knn => "knn",
            Self::Resolution => "resolution",
            Self::NumClusters => "clusters (k)",
            Self::NumPerm => "permutations",
            Self::Go => "GO terms",
            Self::Susie => "SuSiE call",
        }
    }

    pub fn value(self, a: &AnnotateCliArgs) -> String {
        match self {
            Self::Markers if a.markers.is_empty() => "required: enter picks one".into(),
            Self::Markers => file_name(Path::new(a.markers.as_ref())),
            Self::Run => String::new(),
            // The prefix's name: the directory is the run's, in the summary.
            Self::Output => Path::new(a.out.as_ref())
                .file_name()
                .map_or_else(|| a.out.to_string(), |n| n.to_string_lossy().into_owned()),
            Self::Method => match a.method {
                AnnotateMethod::Auto => "auto".into(),
                AnnotateMethod::Enrichment => "enrichment".into(),
                AnnotateMethod::Projection => "projection".into(),
            },
            Self::Projection => a
                .projection
                .map_or(
                    "auto (score where the run decodes)",
                    ProjectionMethod::as_str,
                )
                .into(),
            Self::Knn => a.knn.map_or_else(
                || format!("default ({})", crate::annotate_cmd::default_knn(a.method)),
                |k| k.to_string(),
            ),
            Self::Resolution => format!("{:.2}", a.resolution),
            Self::NumClusters => a.num_clusters.map_or("auto".into(), |k| k.to_string()),
            Self::NumPerm => a.num_perm.to_string(),
            Self::Go => match (&a.gaf, &a.gmt) {
                (Some(_), _) => "--gaf".into(),
                (_, Some(_)) => "--gmt".into(),
                _ if a.go => "on".into(),
                _ => "off".into(),
            },
            Self::Susie => if a.no_susie { "off" } else { "on" }.into(),
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
            Self::Projection => {
                use ProjectionMethod::{Ora, Score, Units};
                a.projection = match (a.projection, up) {
                    (None, true) | (Some(Ora), false) => Some(Score),
                    (Some(Score), true) | (Some(Units), false) => Some(Ora),
                    (Some(Ora), true) | (None, false) => Some(Units),
                    (Some(Units), true) | (Some(Score), false) => None,
                };
            }
            Self::Knn => {
                // On the method's default it is left unset: no `--knn`.
                let d = crate::annotate_cmd::default_knn(a.method);
                let k = step(a.knn.unwrap_or(d), 1, 2);
                a.knn = (k != d).then_some(k);
            }
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
            Self::Susie => a.no_susie = !a.no_susie,
            // Named gene sets are the command line's; only `--go` toggles.
            Self::Go if a.gaf.is_none() && a.gmt.is_none() => a.go = !a.go,
            Self::Go | Self::Output | Self::Markers | Self::Run => {}
        }
    }

    /// One line on what the row sets, for the form's foot.
    pub fn explain(self) -> &'static str {
        match self {
            Self::Markers => "the marker panel the pass scores the clusters' cell types with; enter picks a file",
            Self::Output => "the prefix the new round is written under; enter edits it",
            Self::Method => "how clusters get their cell types: enrichment of the panel's markers, projection, or auto",
            Self::Projection => "projection calls each cell: the marker score on the run's embedding (with its margin), or the co-embedding ORA",
            Self::Knn => "neighbours per cell in the clustering graph",
            Self::Resolution => "Leiden resolution: higher gives more, smaller clusters",
            Self::NumClusters => "a fixed number of clusters, or auto from the resolution",
            Self::NumPerm => "permutations for the enrichment null: more is slower and finer",
            Self::Go => "also score GO terms for each cluster (slower)",
            Self::Susie => "call each cluster by SuSiE (its counts on the panel; types compete for shared markers), or by the enrichment's share",
            Self::Run => "start the pass with these settings",
        }
    }
}

/// What the tree pane shows.
pub enum TreeMode {
    /// The panel's types on the tree lupin built.
    Panel,
    /// The Cell Ontology around a term.
    Ontology(Box<super::ontology::OntologyView>),
    /// Which cell types precede which: the prior `lupin trajectory` builds.
    Order(Box<super::order::OrderView>),
    /// The round's cells on the run's layout, the selected cluster lit.
    Cells,
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
    /// The label (the text) to give cluster `cluster`, typed by hand.
    LabelName { cluster: ClusterId },
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
    /// Move the export of `pdf` to the base name typed; from the round's
    /// cell map when `cells`, else the order view's figures.
    Relocate { pdf: PathBuf, cells: bool },
    /// The output prefix (the text) for passes.
    Output,
    /// Run `lupin trajectory` with the output prefix typed; `replace` is the
    /// prefix whose existing manifest the user has agreed to replace.
    TrajectoryOut { replace: Option<String> },
}

/// Where a pass on the run `source` gets its cluster expression, in words,
/// when the raw counts are not here.
fn expression_note(source: &Path) -> Option<String> {
    use crate::manifest::no_counts::{self, Source};
    let loaded = crate::manifest::run::load(&source.to_string_lossy()).ok()?;
    Some(match no_counts::source(&loaded) {
        Source::Counts => return None,
        Source::Cache(src) => format!(
            "The raw counts are not here: it scores the clusters {} cached and their gene sums, whatever the clustering rows say",
            file_name(&src.file)
        ),
        Source::Decoder(d) => format!(
            "The raw counts are not here and no pass cached its sums: it scores the {} decoder's expected expression ({}), a rough stand-in",
            loaded.manifest.kind.as_str(),
            d.formula()
        ),
        Source::Nothing => format!(
            "The raw counts are not here, no pass cached its sums, and a {} run has no decoder to stand in: the pass will fail",
            loaded.manifest.kind.as_str()
        ),
    })
}

/// A file the main loop opens the file browser for, in a popup.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FileWant {
    /// A marker panel, for the annotation form.
    Markers,
    /// A `cell<TAB>type` labels file for a trajectory run.
    Labels,
    /// A precedence file for a trajectory run (`--prior`).
    Prior,
    /// A `label<TAB>CL:id` file mapping the labels to Cell Ontology terms.
    LabelCl,
}

impl FileWant {
    /// The file browser's title.
    pub fn title(self) -> &'static str {
        match self {
            Self::Markers => "Pick a marker panel",
            Self::Labels => "Pick a cell<TAB>type labels file",
            Self::Prior => "Pick a precedence file (from<TAB>to<TAB>precedes|unrelated)",
            Self::LabelCl => "Pick a label<TAB>CL:id file",
        }
    }
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
    pub report: Option<crate::progress::Report>,
    /// The job's time when `report` came.
    pub reported_at: Duration,
}

impl JobProgress {
    /// The estimate `elapsed` into the job.
    pub fn estimate(&self, elapsed: Duration) -> Option<crate::progress::Estimate> {
        let r = self.report.as_ref()?;
        Some(crate::progress::estimate(r, self.reported_at, elapsed))
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

    /// Whether `x` stops it: a save is not stopped midway.
    #[must_use]
    pub fn stoppable(&self) -> bool {
        !matches!(self, Self::Save(_))
    }
}

/// A `lupin relabel --preview` rescoring the round against `edits`.
pub struct Rescoring {
    child: Held,
    out: Receiver<String>,
    edits: Vec<Edit>,
    file: PathBuf,
}

impl Rescoring {
    fn stop(mut self) {
        self.child.kill();
        let _ = self.child.wait();
        let _ = std::fs::remove_file(&self.file);
    }
}

/// `lupin data suggest-cl` ranking Cell Ontology terms for the panel types
/// without one. It runs beside the other jobs, in a slot of its own: the
/// first ranking (the text model fetched, the ontology encoded) takes
/// minutes, and passes, saves and rescoring go on meanwhile.
pub struct Ranking {
    child: Held,
    /// Its log lines, apart from the other jobs'.
    log: Receiver<String>,
    started: Instant,
    /// Where it writes the suggestions.
    file: PathBuf,
    /// The type `c` was pressed on, whose suggestions open when it is done.
    label: String,
}

impl Ranking {
    #[must_use]
    pub fn new(child: Held, log: Receiver<String>, file: PathBuf, label: String) -> Self {
        Self {
            child,
            log,
            started: Instant::now(),
            file,
            label,
        }
    }

    fn stop(mut self) {
        self.child.kill();
        let _ = self.child.wait();
    }

    /// Its log lines not yet forwarded once the child has ended: its reader
    /// sends the last ones (often the reason it failed) after the exit, so
    /// they are waited for, a moment at most, until the reader is done.
    fn rest_of_log(&self) -> Vec<String> {
        let until = Instant::now() + RANKING_LOG_TAIL;
        std::iter::from_fn(|| {
            let left = until.checked_duration_since(Instant::now())?;
            self.log.recv_timeout(left).ok()
        })
        .collect()
    }

    #[cfg(test)]
    pub fn pid(&self) -> u32 {
        self.child.id()
    }
}

/// How long an ended ranking's last log lines are waited for.
const RANKING_LOG_TAIL: Duration = Duration::from_millis(500);

/// How the ranking's status line begins: it is renewed only over its own.
const RANKING_STATUS: &str = "ranking Cell Ontology terms";

/// A key that needs pressing twice, and what it will do.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Armed {
    Quit,
    Stop,
}

/// A job the TUI started, with the progress it reports.
pub struct Running {
    pub child: Held,
    pub started: Instant,
    pub job: Job,
    /// As its `@progress` lines report it.
    pub progress: JobProgress,
    /// It ended well and the status says so; the next tick reads what it
    /// wrote, once the screen has shown that.
    ended: bool,
}

impl Running {
    pub fn new(child: Held, job: Job) -> Self {
        Self {
            child,
            started: Instant::now(),
            job,
            progress: JobProgress::default(),
            ended: false,
        }
    }
}

pub struct App {
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
    pub cl: Option<ClTerms>,
    /// Each panel type on a CL term with every term above it.
    pub panel_ancestry: Vec<(String, std::collections::BTreeSet<String>)>,
    /// Where the run's data files are found, and aliases remembered.
    pub data_search: crate::manifest::data_files::SearchPath,
    /// Each label's Cell Ontology term ([`Self::term_of`]) as worked out
    /// since the ontology, its aliases or the panel's tree last changed.
    terms_of: std::cell::RefCell<std::collections::HashMap<String, Option<String>>>,
    /// The marker panel, as read: (gene, type) rows.
    pub panel: Vec<(String, String)>,
    /// The marker panel the session started from: label key → genes.
    pub original: std::collections::BTreeMap<String, std::collections::BTreeSet<String>>,
    pub round: Option<RoundView>,
    pub edits: Vec<Edit>,
    pub focus: Focus,
    pub cluster_sel: usize,
    pub gene_sel: usize,
    /// The selected row of the candidates as listed.
    pub cand_sel: usize,
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
    /// Cell Ontology terms suggested for the panel types without one, once
    /// ranked (`c` in the tree pane).
    pub suggestions: Option<super::suggest::Suggestions>,
    /// The ranking those come from, while it runs.
    pub ranking: Option<Ranking>,
    /// The files jobs hand to their child: a folder only the user can
    /// open, made on first use and removed with the app.
    jobs_dir: Option<tempfile::TempDir>,
    /// The trajectory the order view shows, whose edges its table and whose
    /// texts its strip are, as [`Self::keep_order_current`] keeps it.
    pub order_shown: Option<super::order::Shown>,
    pub tree_sel: usize,
    /// The order view's types on the Cell Ontology, drawn beside the
    /// precedence table as the tree pane draws the panel.
    /// The annotation form, while it is open.
    pub form: Option<Form>,
    /// The key guide is open.
    pub help_open: bool,
    /// The guide's first line on screen, and how far it scrolls as last
    /// drawn (0 when it fits).
    pub help_scroll: usize,
    pub help_scrolls: std::cell::Cell<usize>,
    pub prompt: Option<Prompt>,
    /// A choice among options (see [`super::menu`]), which takes keys.
    pub menu: Option<super::menu::Menu>,
    pub log: Vec<String>,
    log_tx: Sender<String>,
    log_rx: Receiver<String>,
    pub child: Option<Running>,
    /// A preview rescoring the round against the marker edits, running.
    pub rescoring: Option<Rescoring>,
    /// The marker edits the scores on screen reflect (none: the round's own).
    scored_for: Vec<Edit>,
    /// The round's own scores while rescored ones are shown.
    pub recorded: Option<super::round::Scores>,
    pub status: String,
    /// When a key last left a reply in the status: while it is recent the
    /// running job's progress does not replace it.
    replied: Option<Instant>,
    /// The pass settings before the first change made since the round on
    /// screen (see [`App::stale`]); none while unchanged.
    settings_base: Option<Vec<String>>,
    /// A key that asked to be pressed again, when, and the status that
    /// asks: the question lapses after [`KEY_REPLY_HOLD`].
    armed: Option<(Armed, Instant, String)>,
    pub quit: bool,
    /// The trajectory figures of the manifests on screen, when they have any.
    pub figures: Option<super::figure_pane::FigurePane>,
    /// The round's cells on its run's layout, when the run has one.
    pub cell_map: Option<super::figure_pane::FigurePane>,
    /// What to run for a layout, when the round's run has none.
    map_hint: Option<String>,
    /// The user put the tree in the cell map's place (`V`).
    map_off: bool,
    /// The Cell Ontology on screen was opened from the cell map, which `o`
    /// goes back to.
    pub ontology_from_map: bool,
    /// The method the candidates are listed by (`S`); the round's own when
    /// none was chosen.
    cand_by: Option<crate::annotate::rounds::Method>,
    /// Every candidate listed, not only the first ten and each method's top (`a`).
    pub cand_all: bool,
    /// The output prefix was given or confirmed; until then the first pass
    /// asks for it.
    pub out_chosen: bool,
    /// The output prefix was given (`-o`) or typed, not only run with: it
    /// stays when another run is opened.
    pub out_typed: bool,
    /// How the order view runs `lupin trajectory`.
    pub trajectory: super::order::TrajectoryRun,
    /// Started as `lupin trajectory`'s TUI: only it has the order view.
    pub trajectory_tui: bool,
    /// A file to pick in the file browser, which the main loop opens.
    pub want_file: Option<FileWant>,
    /// The file browser, open in a popup over the form or menu that asked.
    pub browser: Option<(FileWant, super::picker::Picker)>,
    /// What follows the pass running now (or the form about to start one).
    after_pass: AfterPass,
    /// The terminal's picture protocol, asked for once.
    picker: Option<ratatui_image::picker::Picker>,
    /// The frame being drawn goes out at once: the figures show as last
    /// drawn while a new one is rasterised, so a key's reply is not held
    /// up by drawing one anew. Off, a frame waits for its figures.
    pub quick_frame: std::cell::Cell<bool>,
    /// A quick frame showed a figure not yet drawn anew.
    pub frame_behind: std::cell::Cell<bool>,
}

impl App {
    /// The app on marker panel `panel` (gene, type rows), its tree placed
    /// on the ontology `cl` when there is one.
    pub fn new(
        args: AnnotateCliArgs,
        source: PathBuf,
        target: PathBuf,
        panel: Vec<(String, String)>,
        cl: Option<ClTerms>,
    ) -> Self {
        let (log_tx, log_rx) = channel();
        let mut app = Self {
            args,
            source,
            target,
            fixed_clusters: false,
            // Built below, from `panel`.
            tree: crate::manifest::ontology::panel_tree_on(None, &[]),
            tree_mode: TreeMode::Panel,
            cl,
            panel_ancestry: Vec::new(),
            panel,
            terms_of: Default::default(),
            data_search: crate::manifest::data_files::SearchPath::new(None),
            original: std::collections::BTreeMap::new(),
            round: None,
            edits: Vec::new(),
            focus: Focus::Clusters,
            cluster_sel: 0,
            gene_sel: 0,
            cand_sel: 0,
            go_sel: 0,
            go_shift: 0,
            gene_view: GeneView::Specific,
            hidden: super::genes::GeneFilter::default(),
            show_hidden: false,
            marked: Vec::new(),
            tree_marked: Vec::new(),
            mixed: super::ontology::Mixed::default(),
            suggestions: None,
            ranking: None,
            jobs_dir: None,
            order_shown: None,
            tree_sel: 0,
            form: None,
            help_open: false,
            help_scroll: 0,
            help_scrolls: Default::default(),
            prompt: None,
            menu: None,
            log: Vec::new(),
            log_tx,
            log_rx,
            child: None,
            rescoring: None,
            scored_for: Vec::new(),
            recorded: None,
            status: "r: run a pass".into(),
            replied: None,
            settings_base: None,
            armed: None,
            quit: false,
            shift_enter: false,
            figures: None,
            cell_map: None,
            map_hint: None,
            map_off: false,
            ontology_from_map: false,
            cand_by: None,
            cand_all: false,
            out_chosen: true,
            out_typed: true,
            trajectory: super::order::TrajectoryRun::default(),
            trajectory_tui: false,
            want_file: None,
            browser: None,
            after_pass: AfterPass::Nothing,
            picker: None,
            quick_frame: Default::default(),
            frame_behind: Default::default(),
            progress_hidden: false,
            failed: None,
        };
        app.rebuild_tree();
        app
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
            if let Some(report) = crate::progress::parse(&line) {
                if let Some(r) = &mut self.child {
                    let at = r.started.elapsed();
                    (r.progress.report, r.progress.reported_at) = (Some(report), at);
                }
                continue;
            }
            self.push_log(line);
        }
        for line in super::drain_own_log() {
            self.push_log(line);
        }
        self.keep_scores_current();
        self.keep_order_current();
        self.keep_ranking();
        // A question asked again lapses, and the status with it unless
        // something has said more since.
        if let Some((_, _, asked)) = self
            .armed
            .take_if(|(_, at, _)| at.elapsed() >= KEY_REPLY_HOLD)
        {
            if self.status == asked {
                self.status.clear();
            }
        }
        let Some(r) = &mut self.child else {
            return;
        };
        let (secs, job) = (r.started.elapsed().as_secs(), r.job.clone());
        match r.child.try_wait() {
            // The progress, unless a key waits to be pressed again or a
            // key's reply is still fresh.
            Ok(None)
                if self.armed.is_some()
                    || self.replied.is_some_and(|t| t.elapsed() < KEY_REPLY_HOLD) => {}
            Ok(None) => {
                let elapsed = r.started.elapsed();
                let stop = if job.stoppable() { " (x: stop)" } else { "" };
                self.status = match r.progress.estimate(elapsed) {
                    Some(e) if self.progress_hidden => format!(
                        "{} {:.0}% · {} · b shows{stop}",
                        job.name(),
                        100.0 * e.fraction,
                        crate::progress::eta_text(&e)
                    ),
                    _ => format!("{}… {secs}s{stop}", job.name()),
                };
            }
            // Said first: reading the output can take a moment.
            Ok(Some(st)) if st.success() && !r.ended => {
                r.ended = true;
                self.status = format!("{} done in {secs}s: reading what it wrote…", job.name());
            }
            Ok(Some(st)) if st.success() => {
                self.child = None;
                if let Job::Trajectory(m) = job {
                    self.trajectory_done(&m, secs);
                    return;
                }
                let after = std::mem::take(&mut self.after_pass);
                let after = if job == Job::Pass {
                    after
                } else {
                    AfterPass::Nothing
                };
                // A new pass rewrote the base round: rounds made on the old
                // one no longer apply.
                if job == Job::Pass {
                    if let Err(e) = crate::manifest::rounds::supersede_later(&self.target) {
                        self.push_log(format!("could not set the old rounds aside: {e:#}"));
                    }
                }
                // A save writes the round after the one on screen; a pass,
                // its output's.
                let saved = match (&job, &self.round) {
                    (Job::Save(_), Some(r)) => r.manifest.clone(),
                    _ => self.target.clone(),
                };
                let (latest, _) = chain_rounds(&saved);
                // A save's round was made with the settings of the one before
                // it, so a change to them still waits; a pass ran with them.
                let base = self.settings_base.take();
                // Edits made while a save ran were not in it: carry them into
                // the new round (a save keeps the clusters' ids).
                let later = match job {
                    Job::Save(n) => self.edits.split_off(n.min(self.edits.len())),
                    _ => Vec::new(),
                };
                let by_before = self.round.as_ref().and_then(|r| r.by);
                self.open(&latest);
                self.settings_base = if job == Job::Pass { None } else { base };
                if self.in_order() {
                    self.reload_order();
                }
                let kept = later.len();
                self.edits = later;
                // A marker edit's save drops SuSiE's and the score's
                // readings, made on the old panel: the new round's call is
                // another method's, which may flag other clusters.
                let by_after = self.round.as_ref().and_then(|r| r.by);
                let saved_other = |b: &crate::annotate::rounds::Method| {
                    matches!(job, Job::Save(_)) && by_after != Some(*b)
                };
                let called = match by_before.filter(saved_other) {
                    Some(before) => {
                        let after = by_after.unwrap_or_default().name();
                        self.push_log(format!(
                            "the new round is called by {after}: a save does not redo {}'s call (r re-runs the pass)",
                            before.name()
                        ));
                        format!(
                            ", called by {after}: {}'s call needs a new pass (r)",
                            before.name()
                        )
                    }
                    None => String::new(),
                };
                // The status keeps to what to do next: the files exported
                // are named in the log.
                let exported = |app: &mut Self| match app.export() {
                    Ok(_) => "labels and markers exported".to_string(),
                    Err(e) => e,
                };
                let done = match job {
                    Job::Save(n) if kept > 0 => format!(
                        "saved {n} edit(s){called}, {kept} made since still unsaved; {}",
                        exported(self)
                    ),
                    Job::Save(n) => format!("saved {n} edit(s){called}; {}", exported(self)),
                    _ => format!("pass done in {secs}s"),
                };
                // The round's counts, as opening it says them; the header
                // names it.
                let counts = self.round.as_ref().map(round_counts);
                self.status = match (&job, counts) {
                    (Job::Save(_), Some(c)) => format!("{done}; {c}"),
                    _ => format!("{done}. {}", self.status),
                };
                match after {
                    AfterPass::Run if !self.asking() => self.ask_trajectory_out(),
                    AfterPass::Run => self.status += " r in the order view runs the trajectory.",
                    AfterPass::Offer if !self.asking() => self.offer_trajectory_run(),
                    _ => {}
                }
            }
            Ok(Some(st)) => {
                self.child = None;
                self.after_pass = AfterPass::Nothing;
                self.status = format!("{} failed ({st}); see the log", job.name());
                self.failed = Some(self.status.clone());
            }
            Err(e) => {
                self.child = None;
                self.after_pass = AfterPass::Nothing;
                self.status = format!("lost the {}: {e}", job.name());
                self.failed = Some(self.status.clone());
            }
        }
    }

    /// Take in the ranking's log, and its suggestions once it has ended.
    /// Suggestions are extra: without them (offline, no model) the rest goes
    /// on, so a failed ranking raises no failure popup.
    fn keep_ranking(&mut self) {
        let Some(r) = &mut self.ranking else { return };
        let mut lines: Vec<String> = std::iter::from_fn(|| r.log.try_recv().ok()).collect();
        let secs = r.started.elapsed().as_secs();
        let ended = r.child.try_wait();
        if !matches!(ended, Ok(None)) {
            lines.extend(r.rest_of_log());
        }
        for line in lines
            .into_iter()
            .filter(|l| crate::progress::parse(l).is_none())
        {
            self.push_log(line);
        }
        self.ranking_ended(ended, secs);
    }

    /// What the ranking's child's state, `ended`, after `secs`, means.
    pub(super) fn ranking_ended(&mut self, ended: std::io::Result<Option<ExitStatus>>, secs: u64) {
        match ended {
            Ok(None) => {
                // The other job's progress, or a key's reply, says more.
                // A job's or a key's word stays: the ranking renews its own.
                let quiet = self.child.is_none()
                    && self.armed.is_none()
                    && self.replied.is_none_or(|t| t.elapsed() >= KEY_REPLY_HOLD)
                    && (self.status.is_empty() || self.status.starts_with(RANKING_STATUS));
                if quiet {
                    self.status = format!("{RANKING_STATUS}… {secs}s (x: stop)");
                }
            }
            Ok(Some(st)) if st.success() => {
                let Some(r) = self.ranking.take() else { return };
                self.suggestions_done(&r.file, &r.label, secs);
            }
            _ => {
                if let Some(r) = self.ranking.take() {
                    r.stop();
                }
                self.status =
                    "no Cell Ontology suggestions (offline, or no text model): see the log; the rest works"
                        .into();
                self.push_log(format!("[WARN] {}", self.status));
            }
        }
    }

    /// Drop the suggestions, and a ranking under way: both are of the
    /// terms or the panel as they were.
    fn forget_suggestions(&mut self) {
        self.suggestions = None;
        if let Some(r) = self.ranking.take() {
            r.stop();
            self.push_log(
                "── suggestions: the ranking stopped, the terms or the panel changed".into(),
            );
        }
    }

    /// Stop the ranking.
    fn stop_ranking(&mut self) {
        if let Some(r) = self.ranking.take() {
            r.stop();
            self.status = "stopped the term ranking".into();
        }
    }

    /// The job has ended and the next tick reads what it wrote: keys wait
    /// for that.
    #[must_use]
    pub fn job_ended(&self) -> bool {
        self.child.as_ref().is_some_and(|r| r.ended)
    }

    /// Whether the screen is due a redraw `since` the last one with nothing
    /// else changed: a running job's elapsed time and ETA move on the clock.
    #[must_use]
    pub fn redraw_due(&self, since: Duration) -> bool {
        (self.child.is_some() || self.ranking.is_some()) && since >= Duration::from_secs(1)
    }

    /// The running job for its popup: what it is and how long it has run.
    pub fn running_job(&self) -> Option<(String, std::time::Duration)> {
        let Running { started, job, .. } = self.child.as_ref()?;
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
                self.status = format!("{}: {}", file_name(manifest), round_counts(&r));
                self.set_round(Some(r));
            }
            Err(e) => {
                self.status = format!("could not read {}: {e:#}", manifest.display());
                self.push_log(format!("[WARN] {}", self.status));
            }
        }
    }

    /// Show the round the TUI opens on: in annotate, the latest round an
    /// earlier session left under the output; else the manifest opened,
    /// when it is `annotated`. The trajectory TUI starts from the manifest
    /// named, as the batch run would.
    pub fn open_start(&mut self, annotated: bool) {
        if self.target.is_file() && !self.trajectory_tui {
            let (latest, _) = chain_rounds(&self.target);
            self.open(&latest);
        } else if annotated {
            let opened = self.source.clone();
            self.open(&opened);
        }
    }

    /// Put `r` on screen (or none), dropping what belonged to the round
    /// shown before: edits, marks, scores and the stale note.
    fn set_round(&mut self, r: Option<RoundView>) {
        self.stop_rescoring();
        self.recorded = None;
        self.scored_for.clear();
        let mut r = r;
        self.load_cell_map(r.as_mut().map(|r| std::mem::take(&mut r.figures)));
        self.round = r;
        self.edits.clear();
        self.marked.clear();
        self.tree_marked.clear();
        self.gene_sel = 0;
        self.settings_base = None;
        self.cluster_sel = self.cluster_sel.min(self.n_clusters().saturating_sub(1));
        self.sync_map_focus();
    }

    /// The round's cells (`figures`) on its run's layout, shown in the tree
    /// pane; the panel's tree when the run has no layout.
    /// A new round of the same run keeps the map's look and the user's
    /// choice of the tree over it.
    fn load_cell_map(&mut self, figures: Option<crate::figures::CellFigures>) {
        let before = self.cell_map.take();
        self.map_hint = None;
        if let Some(data) = figures {
            if data.layouts.is_empty() {
                self.map_hint = data.layout_hint;
            } else {
                let picker = self.picker().clone();
                let mut map = super::figure_pane::FigurePane::showing(data, &picker);
                if let Some(b) = &before {
                    map.keep_look(b);
                }
                self.cell_map = Some(map);
            }
        }
        match (&self.tree_mode, self.cell_map.is_some()) {
            (TreeMode::Panel, true) if !self.map_off => self.tree_mode = TreeMode::Cells,
            (TreeMode::Cells, false) => self.tree_mode = TreeMode::Panel,
            _ => {}
        }
    }

    /// Light the selected cluster on the cell map.
    fn sync_map_focus(&mut self) {
        let k = self.selected().map(|c| c.id as i32);
        if let Some(m) = &mut self.cell_map {
            m.style.focus = k;
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
            self.status = if round.by == Some(crate::annotate::rounds::Method::Score) {
                "a projection score round has no statistics to rescore: a new pass scores the edited panel"
            } else {
                "this round cached no statistics: marker edits are scored on save"
            }
            .into();
            return;
        }
        let started = (|| -> anyhow::Result<_> {
            let file = job_file(
                &mut self.jobs_dir,
                &format!("rescore-{}.jsonl", self.edits.len()),
            )?;
            let lines: Vec<String> = decisions(&want, round, None)
                .iter()
                .map(serde_json::to_string)
                .collect::<Result<_, _>>()?;
            std::fs::write(&file, lines.join("\n") + "\n")?;
            let (child, out) = runner::spawn_preview(&round.manifest, &file, self.log_tx.clone())?;
            Ok((child, out, file))
        })();
        match started {
            // The table's title says it is rescoring; the status keeps what
            // the edit said.
            Ok((child, out, file)) => {
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
        let preview = match serde_json::from_str::<serde_json::Value>(json) {
            Ok(v) => v,
            Err(_) => return self.rescore_failed("rescoring gave no readable result; see the log"),
        };
        let rescored = super::round::parse_scores(&preview);
        let Some(rescored) = rescored else {
            return self
                .rescore_failed("this round cannot be rescored: marker edits are scored on save");
        };
        let Some(round) = &mut self.round else { return };
        let scores = round.rescored(self.recorded.as_ref(), &rescored);
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
    pub(super) fn figures_shown(&self) -> bool {
        self.figures.as_ref().is_some_and(|f| f.shown)
    }

    /// The keys that start a pass from the form, in words.
    pub(super) fn run_keys(&self) -> &'static str {
        if self.shift_enter {
            "shift+enter or enter on ▶ run"
        } else {
            "enter on ▶ run (the last row)"
        }
    }

    /// A prompt or menu waits for an answer.
    pub(super) fn asking(&self) -> bool {
        self.prompt.is_some() || self.menu.is_some()
    }

    /// Keys go to the panel tree, and to nothing on top of it: no popup,
    /// form, help or file picker, and no search being typed (the panel tree
    /// shows none).
    fn idle_on_panel_tree(&self) -> bool {
        self.focus == Focus::Tree
            && matches!(self.tree_mode, TreeMode::Panel)
            && !self.asking()
            && self.failed.is_none()
            && self.browser.is_none()
            && self.form.is_none()
            && !self.help_open
    }

    /// The trajectory's order view is up (the tree pane shows it).
    pub(super) fn in_order(&self) -> bool {
        matches!(self.tree_mode, TreeMode::Order(_))
    }

    fn n_clusters(&self) -> usize {
        self.round.as_ref().map_or(0, |r| r.clusters.len())
    }

    pub fn selected(&self) -> Option<&super::round::ClusterView> {
        self.round.as_ref()?.clusters.get(self.cluster_sel)
    }

    /// Select cluster `id` again, as a choice made for it may outlast the
    /// selection (a pass reloads the round); `false`, saying so, when the
    /// round no longer has it.
    fn reselect(&mut self, id: ClusterId) -> bool {
        let at = self
            .round
            .as_ref()
            .and_then(|r| r.clusters.iter().position(|c| c.id == id));
        match at {
            Some(i) => {
                if i != self.cluster_sel {
                    self.select_cluster(i);
                }
                true
            }
            None => {
                self.status = format!("K{id} is no longer in the round");
                false
            }
        }
    }

    /// The method the candidates are listed by: the one chosen (`S`), else
    /// the round's, while it ranks the selected cluster's candidates; else
    /// the first that does.
    pub fn cand_method(&self) -> crate::annotate::rounds::Method {
        use crate::annotate::rounds::Method;
        let ranks = |m: &Method| {
            self.selected()
                .is_some_and(|c| c.candidates.iter().any(|t| t.ranks[m.index()].is_some()))
        };
        let wanted = self.cand_by.or_else(|| self.round.as_ref()?.by);
        wanted
            .filter(ranks)
            .or_else(|| Method::ALL.into_iter().find(ranks))
            .unwrap_or_default()
    }

    /// The selected cluster's candidates as listed: by [`App::cand_method`],
    /// the first ten and each method's top, or all (`a`).
    pub fn listed_candidates(&self) -> Vec<&super::round::Candidate> {
        self.selected().map_or_else(Vec::new, |c| {
            super::round::listed(&c.candidates, self.cand_method(), self.cand_all)
        })
    }

    /// List the candidates by the next method that ranks any of them.
    fn next_cand_method(&mut self) {
        use crate::annotate::rounds::Method;
        let ranking: Vec<Method> = Method::ALL
            .into_iter()
            .filter(|m| {
                self.selected()
                    .is_some_and(|c| c.candidates.iter().any(|t| t.ranks[m.index()].is_some()))
            })
            .collect();
        let next =
            (!ranking.is_empty()).then(|| crate::figures::next_in(&ranking, self.cand_method()));
        self.status = match next {
            Some(m) => {
                self.cand_by = Some(m);
                format!("candidates by {}", m.name())
            }
            None => "no method ranks this cluster's candidates".into(),
        };
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
        self.cand_sel = 0;
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
        let own = c.candidates.iter().find(|t| t.label == label);
        let reason = match own {
            Some(t) => format!("keep {label} ({})", evidence(t)),
            None => format!("keep {label}"),
        };
        let id = c.id;
        let reason = self.rescored_note(reason, own);
        self.prompt = Some(Prompt {
            title: format!(" K{id} keeps {label}: why? "),
            text: reason,
            pending: Pending::Keep { cluster: id },
        });
    }

    /// Open the annotation form: the marker panel, the output and the pass
    /// settings, run from its last row; `after` follows the pass.
    pub(super) fn open_form(&mut self, after: AfterPass) {
        if self.child.is_some() {
            self.status = "wait for the running job, or stop it with x".into();
            return;
        }
        self.after_pass = after;
        // From the trajectory's view a pass writes a new round beside the
        // run unless one was chosen.
        if !self.out_chosen && self.in_order() {
            let (_, out) = crate::manifest::family::pass_origin(&self.source);
            self.target = annotated_path(&self.source, &out);
            self.args.out = out.into();
        }
        self.form = Some(Form {
            row: setting_row(
                &self.args,
                if self.args.markers.is_empty() {
                    Setting::Markers
                } else {
                    Setting::Run
                },
            ),
            // Asked each time: the counts may have come or gone since.
            expression: expression_note(&self.source),
            ..Form::default()
        });
    }

    /// Close the form without running.
    fn close_form(&mut self) {
        self.form = None;
        self.after_pass = AfterPass::Nothing;
    }

    /// What a pass from the form would replace or drop, in words.
    pub(super) fn pass_effects(&self) -> Vec<String> {
        let (_, later) = chain_rounds(&self.target);
        let mut what = Vec::new();
        if self.target.is_file() {
            what.push(format!(
                "replaces {} and its outputs",
                file_name(&self.target)
            ));
        }
        if !self.edits.is_empty() {
            what.push(format!("drops {} unsaved edit(s)", self.edits.len()));
        }
        if !later.is_empty() {
            what.push(format!("sets aside {} saved round(s)", later.len()));
        }
        what
    }

    /// The form's first line: what running it will do, and whether that
    /// replaces or drops anything.
    pub fn form_summary(&self) -> (String, bool) {
        let run = file_name(&self.source);
        let out = file_name(&self.target);
        let effects = self.pass_effects();
        let warns = !effects.is_empty();
        let what = if effects.is_empty() {
            format!("a new pass on {run} writes {out} and replaces nothing")
        } else {
            format!("a new pass on {run} {}", effects.join(", "))
        };
        let then = match self.after_pass {
            AfterPass::Run => "; then the trajectory on its labels",
            AfterPass::Offer => "; then a choice to run the trajectory",
            AfterPass::Nothing => "",
        };
        // The note on missing counts is about the enrichment's input.
        let note = self
            .form
            .as_ref()
            .and_then(|f| f.expression.as_deref())
            .filter(|_| self.args.method != AnnotateMethod::Projection);
        match note {
            Some(e) => (format!("{what}{then}. {e}"), true),
            None => (format!("{what}{then}"), warns),
        }
    }

    /// Run the form: a marker panel is needed, and a pass that replaces
    /// something runs on the second asking.
    fn run_form(&mut self) {
        let busy = self.child.is_some();
        let no_panel = self.args.markers.is_empty();
        let effects = self.pass_effects();
        let Some(form) = &mut self.form else { return };
        if busy {
            form.note = Some("wait for the running job, or stop it with x".into());
            return;
        }
        if no_panel {
            form.confirm = false;
            form.note =
                Some("a pass needs a marker panel: enter on the first row picks one".into());
            form.row = setting_row(&self.args, Setting::Markers);
            return;
        }
        self.out_chosen = true;
        if !effects.is_empty() && !form.confirm {
            form.confirm = true;
            form.note = Some(format!(
                "this pass {}: run again to confirm",
                effects.join(", ")
            ));
            return;
        }
        self.start();
        if self.child.is_some() {
            self.form = None;
        }
    }

    /// Start a pass with the form's settings (`run_form` has checked them).
    fn start(&mut self) {
        let replaces = self.target.is_file();
        self.push_log(format!(
            "── pass: lupin annotate {}",
            self.args.to_argv().join(" ")
        ));
        match runner::spawn_pass(&self.args, replaces, self.log_tx.clone()) {
            Ok(c) => {
                self.child = Some(Running::new(c, Job::Pass));
                self.status = "running…".into();
            }
            Err(e) => self.status = format!("{e:#}"),
        }
    }

    /// Stop a running pass or trajectory. A save is not stopped midway,
    /// which could leave a round half-written: it finishes on its own.
    fn stop(&mut self) {
        match self.child.as_ref().map(|r| r.job.clone()) {
            Some(job @ (Job::Pass | Job::Trajectory(_))) => {
                self.after_pass = AfterPass::Nothing;
                let Some(mut r) = self.child.take() else {
                    return;
                };
                r.child.kill();
                let _ = r.child.wait();
                if job != Job::Pass {
                    self.status = "stopped".into();
                    return;
                }
                // A pass writes under a staging prefix until it is done:
                // stopped, its half goes and the round stays as it was.
                let staging = crate::manifest::staging::staging_prefix(&self.args.out);
                if let Err(e) = crate::manifest::staging::discard(&staging) {
                    self.push_log(format!("[WARN] the stopped pass's files: {e:#}"));
                }
                self.status = format!("stopped: {} is as it was", self.args.out);
            }
            Some(Job::Save(_)) => self.status = "a save finishes on its own; wait for it".into(),
            None => {}
        }
    }

    /// How many of the edits a running save is writing: those stay as they
    /// are until it is done.
    fn saving(&self) -> usize {
        match self.child.as_ref().map(|r| r.job.clone()) {
            Some(Job::Save(n)) => n.min(self.edits.len()),
            _ => 0,
        }
    }

    /// Record a cluster's label or keep, in place of its earlier one (unless
    /// a running save is writing that).
    fn decide(&mut self, edit: Edit) {
        let saving = self.saving();
        let id = edit.cluster();
        let mut k = 0;
        self.edits.retain(|e| {
            k += 1;
            k <= saving || e.cluster() != id
        });
        self.edits.push(edit);
    }

    /// Whether the pass settings differ from those the round on screen was
    /// made with (as they stood when it was shown): not while a pass runs
    /// with them, nor once they are set back, nor with no round on screen.
    pub fn stale(&self) -> bool {
        !self.pass_running()
            && self.round.is_some()
            && self
                .settings_base
                .as_ref()
                .is_some_and(|b| *b != pass_settings(&self.args))
    }

    /// Whether a pass is running, which is about to replace the clusters.
    fn pass_running(&self) -> bool {
        self.child.as_ref().is_some_and(|r| r.job == Job::Pass)
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
                self.child = Some(Running::new(c, Job::Save(self.edits.len())));
                self.status = "saving…".into();
            }
            Err(e) => self.status = format!("save failed: {e:#}"),
        }
    }

    /// Write the open round's cell annotation and marker table, logging
    /// where; says where, or why not.
    fn export(&mut self) -> Result<String, String> {
        let Some(r) = &self.round else {
            return Err("nothing to export".into());
        };
        match super::export::write(r, &self.tree, self.cl.as_ref(), &self.mixed, &self.original) {
            Ok((cells, markers)) => {
                self.push_log(format!(
                    "exported {} and {}",
                    cells.display(),
                    markers.display()
                ));
                Ok(format!(
                    "exported {}, {}",
                    file_name(&cells),
                    file_name(&markers)
                ))
            }
            Err(e) => Err(format!("export failed: {e:#}")),
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

    /// Label the selected cluster by its `k`-th listed candidate.
    fn pick_candidate(&mut self, k: usize) {
        if let Some(pick) = self.listed_candidates().get(k).copied().cloned() {
            self.pick(pick);
        }
    }

    /// Open the label choices for the selected cluster: each method's top
    /// types ([`super::round::choices`]), then the tree and a typed name.
    fn open_choice(&mut self) {
        use super::menu::{numbered_but, Choice, Menu};
        let Some(c) = self.selected() else { return };
        let options = super::round::choices(&c.candidates)
            .into_iter()
            .map(|t| {
                let detail = self.choice_detail(c, t);
                let label = t.label.clone();
                let action = Action::Label {
                    cluster: c.id,
                    label: label.clone(),
                };
                (label, detail, action)
            })
            .collect();
        let mut choices = numbered_but(options, "tn");
        choices.extend([
            Choice {
                key: Some('t'),
                label: "on the tree…".into(),
                detail: "Pick the type, or a coarser node, on the panel tree.".into(),
                action: Action::LabelOnTree,
            },
            Choice {
                key: Some('n'),
                label: "a name…".into(),
                detail: "Type any name; it is scored once the edits are saved.".into(),
                action: Action::LabelByName { cluster: c.id },
            },
        ]);
        self.menu = Some(Menu::new(
            format!("Label K{}: each method's top types first", c.id),
            choices,
            "",
        ));
    }

    /// A label choice's evidence, one paragraph a line: each method's
    /// reading and its place in that method's ranking, the type's markers,
    /// and the cluster's gene-set terms.
    fn choice_detail(&self, c: &super::round::ClusterView, t: &super::round::Candidate) -> String {
        use crate::annotate::rounds::Method;
        let num = |v: Option<f32>| v.map_or_else(|| "—".into(), |v| format!("{v:.2}"));
        let place = |m: Method| match t.ranks[m.index()] {
            Some(0) => " · its top".to_string(),
            Some(r) => format!(" · its #{}", r + 1),
            None => String::new(),
        };
        let mut lines = Vec::new();
        if t.enrichment.is_some() || t.nes.is_some() {
            lines.push(format!(
                "enrichment: share {} · NES {} · q {}{}",
                num(t.enrichment),
                num(t.nes),
                num(t.q),
                place(Method::Enrichment)
            ));
        }
        if let Some(s) = t.susie() {
            lines.push(format!(
                "SuSiE: PIP {:.2} · ×{:.1} · explained {}{}",
                s.pip,
                s.fold,
                num(t.explained),
                place(Method::Susie)
            ));
        }
        if let Some(s) = t.score {
            lines.push(format!(
                "score: {s:.2} of the cluster's cells{}",
                place(Method::Score)
            ));
        }
        if let Some(r) = &self.round {
            let markers: Vec<String> = (r.markers_of(&t.label, &self.edits).into_iter())
                .map(|(g, _)| g)
                .take(12)
                .collect();
            if !markers.is_empty() {
                lines.push(format!("markers: {}", markers.join(" ")));
            }
        }
        if !c.terms.is_empty() {
            let terms: Vec<&str> = c.terms.iter().take(5).map(|t| t.term.as_str()).collect();
            lines.push(format!("the cluster's gene sets: {}", terms.join("; ")));
        }
        lines.join("\n")
    }

    /// Label the selected cluster by candidate `pick`.
    fn pick(&mut self, pick: super::round::Candidate) {
        let Some(c) = self.selected() else { return };
        let reason = match c.candidates.first() {
            Some(top) if top.label != pick.label => format!(
                "{} ({}) over the top candidate {} ({})",
                pick.label,
                evidence(&pick),
                top.label,
                evidence(top)
            ),
            _ => format!("top candidate {} ({})", pick.label, evidence(&pick)),
        };
        let reason = self.rescored_note(reason, Some(&pick));
        self.ask_label(pick.evidence.label, reason, None);
    }

    /// `reason`, noting when the numbers it quotes of candidate `t` are a
    /// rescoring's, and which a rescoring did not redo.
    fn rescored_note(&self, reason: String, t: Option<&super::round::Candidate>) -> String {
        if self.recorded.is_none() {
            return reason;
        }
        let old: Vec<&str> = t
            .into_iter()
            .flat_map(|t| {
                [
                    (t.pip.is_some(), "SuSiE's"),
                    (t.score.is_some(), "the score's"),
                ]
            })
            .filter_map(|(has, name)| has.then_some(name))
            .collect();
        if old.is_empty() {
            format!("{reason}; rescored with the marker edits")
        } else {
            format!(
                "{reason}; the enrichment's rescored with the marker edits, {} from before them",
                old.join(" and ")
            )
        }
    }

    /// Label the selected cluster by tree node `i`.
    fn pick_node(&mut self, i: usize) {
        let label = self.tree.label(i).to_string();
        let key = label_key(&label);
        let under = self.tree.labels_under(i);
        // A panel type is the label's own evidence; the types under it pool.
        let own = self
            .selected()
            .and_then(|c| c.candidates.iter().find(|t| label_key(&t.label) == key))
            .map(evidence);
        let pooled = self.pooled(|t| {
            let t = label_key(t);
            t != key && under.contains(&t.as_str())
        });
        let mut reason = if self.tree.nodes[i].labels.is_empty() {
            format!("{label} from the ontology")
        } else {
            format!("{label} on the panel tree")
        };
        if let Some(own) = own {
            reason += &format!(" ({own})");
        }
        if !pooled.is_empty() {
            reason += &format!(", over {}", pooled.join(" + "));
        }
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
            }) => {
                let (did, to) = if add {
                    ("adding", "to")
                } else {
                    ("dropping", "from")
                };
                format!("undid {did} {} {to} {label}'s markers", genes.join(", "))
            }
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
        let folds: Vec<Option<f32>> = genes.iter().map(|g| r.fold_change(c.id, g)).collect();
        let fc: Vec<String> = genes
            .iter()
            .zip(&folds)
            .map(|(g, f)| f.map_or_else(|| format!("{g} not measured"), |f| format!("{g} {f:+.1}")))
            .collect();
        // Specific (or not) only as every gene's fold change says.
        let up = |f: &Option<f32>| f.is_some_and(|f| f > 0.0);
        let down = |f: &Option<f32>| f.is_some_and(|f| f <= 0.0);
        let (id, fc) = (c.id, fc.join(", "));
        let (verb, reason) = match add {
            true if folds.iter().all(up) => ("add", format!("specific to K{id} ({fc})")),
            true => ("add", format!("markers of {label} (in K{id}: {fc})")),
            false if folds.iter().all(down) => ("drop", format!("not specific to K{id} ({fc})")),
            false => ("drop", format!("not markers of {label} (in K{id}: {fc})")),
        };
        self.prompt = Some(Prompt {
            title: format!(" {verb} {} as markers of {label}: why? ", genes.join(", ")),
            text: reason,
            pending: Pending::Markers { label, genes, add },
        });
    }

    /// Keys of the open menu: move, take an option, or cancel.
    fn menu_key(&mut self, code: KeyCode) {
        let Some(menu) = &mut self.menu else { return };
        match menu.key(code) {
            Outcome::Open => {}
            Outcome::Cancel => {
                let note = self.menu.take().map(|m| m.on_cancel).unwrap_or_default();
                // Nothing follows a pass the menu was asked about.
                self.after_pass = AfterPass::Nothing;
                self.status = note;
            }
            Outcome::Chosen(a) => {
                self.menu = None;
                self.choose(a);
            }
        }
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
                        // As the panel writes its types, so every pane matches it.
                        self.decide(Edit::Label {
                            cluster,
                            label: label_key(&label),
                            reason,
                        });
                        self.focus = Focus::Clusters;
                        if let (Some((alias, id)), Some(file)) =
                            (remember, self.data_search.alias_target())
                        {
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
                        // Keeping a label picked here keeps that label, not
                        // the round's.
                        let now = self.round.as_ref().and_then(|r| {
                            let own = &r.clusters.iter().find(|c| c.id == cluster)?.label;
                            let now = r.label_of(cluster, &self.edits);
                            (now != *own).then(|| now.unwrap_or_else(|| UNASSIGNED_LABEL.into()))
                        });
                        self.decide(match now {
                            Some(label) => Edit::Label {
                                cluster,
                                label,
                                reason,
                            },
                            None => Edit::Keep { cluster, reason },
                        });
                    }
                    Pending::Hide => self.hide(&reason),
                    Pending::LabelName { cluster } => {
                        if self.reselect(cluster) {
                            self.ask_label(reason, "named by hand".into(), None);
                        }
                    }
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
                    Pending::Output => {
                        self.set_output(&reason);
                        if let Some(f) = &mut self.form {
                            f.confirm = false;
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
                    Pending::Relocate { pdf, cells } => {
                        let pane = if cells {
                            &mut self.cell_map
                        } else {
                            &mut self.figures
                        };
                        self.status = match pane {
                            Some(v) => v
                                .exports
                                .relocate(&pdf, &reason)
                                .unwrap_or_else(|e| format!("could not move: {e:#}")),
                            None => "no figures".into(),
                        };
                    }
                    Pending::Remember { label, id, file } => {
                        self.status = match remember_alias(&file, &label, &id, &reason) {
                            Ok(()) => match self.learn_alias(&label, &id, &file) {
                                None => {
                                    format!("remembered {label} → {id} in {}", file.display())
                                }
                                Some(later) => format!(
                                    "wrote {label} → {id} to {}, but {} maps {label} too and is read after it: that reading stays (edit it there)",
                                    file.display(),
                                    later.display()
                                ),
                            },
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
        // The status says what the latest key did: an older one, or a
        // question the key has dropped (`x again`, `q again`), goes.
        self.status.clear();
        self.key_unsynced(k);
        self.replied = (!self.status.is_empty()).then(Instant::now);
        self.sync_map_focus();
    }

    fn key_unsynced(&mut self, k: KeyEvent) {
        // Any key but a timely second press drops a question asked again.
        let armed = self
            .armed
            .take()
            .filter(|(_, at, _)| at.elapsed() < KEY_REPLY_HOLD)
            .map(|(a, ..)| a);
        if k.modifiers.contains(KeyModifiers::CONTROL) && k.code == KeyCode::Char('c') {
            return self.request_quit(armed, "ctrl-c");
        }
        // The popups take keys from the top down, as `ui::draw_overlays`
        // stacks them: a failed job's notice goes with any key, which does
        // nothing else; then a prompt or menu.
        if self.failed.take().is_some() {
            return;
        }
        if self.browser.is_some() {
            return self.browser_key(k.code);
        }
        if self.menu.is_some() {
            return self.menu_key(k.code);
        }
        if self.prompt.is_some() {
            return self.prompt_key(k);
        }
        let beside = self.focus == Focus::Tree && self.in_order();
        let view = match &mut self.tree_mode {
            TreeMode::Ontology(v) => Some(v.as_mut()),
            TreeMode::Order(v) if beside => v.ontology.as_mut(),
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
        if self.help_open {
            // The arrows, pages and ends scroll the guide; any other key
            // closes it.
            if !step(&mut self.help_scroll, self.help_scrolls.get() + 1, k.code) {
                self.help_open = false;
            }
            return;
        }
        // `b` hides the running job's popup, or brings it back.
        if k.code == KeyCode::Char('b') && self.child.is_some() {
            self.progress_hidden = !self.progress_hidden;
            return;
        }
        if k.code == KeyCode::Char('?') {
            self.help_open = true;
            self.help_scroll = 0;
            return;
        }
        if self.form.is_some() {
            return self.settings_key(k);
        }
        // `g` lists the run's family from any pane.
        if k.code == KeyCode::Char('g') {
            return self.list_runs();
        }
        // `x` stops a running job from any pane, asked twice; a save, which
        // is not stopped, says so at once.
        if k.code == KeyCode::Char('x') {
            match self.child.as_ref().map(|r| r.job.stoppable()) {
                Some(false) => self.stop(),
                Some(true) if armed == Some(Armed::Stop) => self.stop(),
                Some(true) => self.ask_again(Armed::Stop, "x again stops it".into()),
                // With no other job, x is the ranking's.
                None if self.ranking.is_none() => {}
                None if armed == Some(Armed::Stop) => self.stop_ranking(),
                None => self.ask_again(Armed::Stop, "x again stops the term ranking".into()),
            }
            return;
        }
        // The focused pane's keys come first.
        if self.pane_key(k.code) {
            return;
        }
        match k.code {
            KeyCode::Char('q') => self.request_quit(armed, "q"),
            // In the trajectory's view `r` runs the trajectory from any pane;
            // a pass (the cluster settings) is for a run with no annotation,
            // offered by that menu.
            KeyCode::Char('r') if self.in_order() => {
                self.ask_trajectory_out();
            }
            KeyCode::Char('r') => self.open_form(AfterPass::Nothing),
            KeyCode::Char('s') => self.save(),
            KeyCode::Char('e') => {
                self.status = if self.edits.is_empty() {
                    self.export().unwrap_or_else(|e| e)
                } else {
                    "save the edits first (s): the export is of a saved round".into()
                };
            }
            KeyCode::Tab => self.focus = self.cycle(true),
            KeyCode::BackTab => self.focus = self.cycle(false),
            _ => {}
        }
    }

    /// Ask for `a`'s key again with `status`, for [`KEY_REPLY_HOLD`].
    fn ask_again(&mut self, a: Armed, status: String) {
        self.status.clone_from(&status);
        self.armed = Some((a, Instant::now(), status));
    }

    /// Quit on `key` (`q` or `ctrl-c`), asking again first while a pass runs
    /// or edits are unsaved; never during a save.
    fn request_quit(&mut self, armed: Option<Armed>, key: &str) {
        let busy = self.child.as_ref().map(|r| r.job.clone());
        let unsaved = match self.edits.len() {
            0 => String::new(),
            n => format!(" and {n} unsaved edit(s)"),
        };
        if matches!(busy, Some(Job::Save(_))) {
            self.status = format!("saving… {key} once it is done");
        } else if let (true, Some(job)) = (armed != Some(Armed::Quit), &busy) {
            let asked = format!(
                "a {} is running{unsaved}: {key} again stops it and quits",
                job.name()
            );
            self.ask_again(Armed::Quit, asked);
        } else if armed == Some(Armed::Quit) || self.edits.is_empty() {
            self.stop();
            self.stop_rescoring();
            self.quit = true;
        } else {
            let asked = format!(
                "{} unsaved edit(s): {key} again to quit, s to save",
                self.edits.len()
            );
            self.ask_again(Armed::Quit, asked);
        }
    }

    /// `lupin kill`: quit as `q` does, asking nothing. A save under way
    /// finishes first (its round taken in), a running pass or trajectory,
    /// the rescoring and the term ranking are stopped, and unsaved edits are
    /// dropped.
    pub fn kill(&mut self) {
        while self.saving_now() {
            if let Some(r) = &mut self.child {
                let _ = r.child.wait();
            }
            self.tick();
        }
        self.stop();
        self.stop_rescoring();
        if let Some(r) = self.ranking.take() {
            r.stop();
        }
        self.edits.clear();
        self.quit = true;
    }

    /// Whether a save is running.
    pub fn saving_now(&self) -> bool {
        self.child.as_ref().is_some_and(|r| r.child.saves())
    }

    /// One line for `lupin sessions` and the resume picker: the run or round
    /// on screen, the job running and the unsaved edits.
    pub fn session_status(&self) -> String {
        let shown = self.round.as_ref().map_or(&self.source, |r| &r.manifest);
        let mut line = crate::manifest::family::file_name(shown);
        if let Some(r) = &self.child {
            line += &format!(", {} running", r.job.name());
        }
        if !self.edits.is_empty() {
            line += &format!(", {} unsaved edit(s)", self.edits.len());
        }
        line
    }

    /// Keys in the annotation form: ↑↓ move, ←→ change a value, Enter picks
    /// the panel, edits the output or runs (on `▶ run`); Shift+Enter runs
    /// from any row where the terminal tells it from Enter.
    fn settings_key(&mut self, k: KeyEvent) {
        let code = k.code;
        let Some(form) = &mut self.form else { return };
        let rows = settings(&self.args);
        let at = form.setting(&self.args);
        match code {
            KeyCode::Esc | KeyCode::Char('r' | 'q') => self.close_form(),
            KeyCode::Up
            | KeyCode::Down
            | KeyCode::Home
            | KeyCode::End
            | KeyCode::PageUp
            | KeyCode::PageDown => {
                let from = form.row;
                step(&mut form.row, rows.len(), code);
                // Clusters from the manifest: the cursor passes over the
                // Leiden rows, which sit between others.
                while self.fixed_clusters && rows[form.row].leiden() {
                    if form.row > from {
                        form.row += 1;
                    } else {
                        form.row -= 1;
                    }
                }
            }
            KeyCode::Left | KeyCode::Right if self.fixed_clusters && at.leiden() => {}
            KeyCode::Left | KeyCode::Right => {
                let before = pass_settings(&self.args);
                at.adjust(&mut self.args, code == KeyCode::Right);
                if self.round.is_some() {
                    self.settings_base.get_or_insert(before);
                }
            }
            KeyCode::Enter if self.shift_enter && k.modifiers.contains(KeyModifiers::SHIFT) => {
                self.run_form();
            }
            KeyCode::Enter => match at {
                Setting::Markers => self.want_file = Some(FileWant::Markers),
                Setting::Output => self.ask_output(),
                Setting::Run => self.run_form(),
                _ => self.status = format!("{} starts the pass", self.run_keys()),
            },
            _ => {}
        }
    }

    fn cycle(&self, forward: bool) -> Focus {
        use Focus::{Candidates, Clusters, Genes, Go, Order, Tree};
        // The order view's columns, left to right: clusters (genes and
        // candidates below them), ordering, ontology; annotate's: clusters
        // and genes, then the tree and the candidates below it.
        let mut order = if self.in_order() {
            vec![Clusters, Genes, Candidates, Order, Tree]
        } else {
            vec![Clusters, Genes, Tree, Candidates]
        };
        if self.has_go() && !self.in_order() {
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

    /// The keys on the selected cluster's candidates, from the clusters or
    /// the candidates pane: a number takes that row, `S` lists them by the
    /// next method and `a` lists all; `false` for any other key.
    fn candidate_key(&mut self, code: KeyCode) -> bool {
        if !matches!(code, KeyCode::Char('1'..='9' | 'S' | 'a')) {
            return false;
        }
        if self.selected().is_none() {
            self.status = "no clusters yet".into();
            return true;
        }
        match code {
            KeyCode::Char(c @ '1'..='9') => {
                let k = c as usize - '1' as usize;
                let n = self.listed_candidates().len();
                if k < n {
                    self.pick_candidate(k);
                } else {
                    self.status = format!("this cluster has {n} candidate(s)");
                }
            }
            KeyCode::Char('S') => {
                self.next_cand_method();
                self.cand_sel = 0;
            }
            _ => {
                self.cand_all = !self.cand_all;
                self.cand_sel = 0;
                self.status = if self.cand_all {
                    "all candidates listed".into()
                } else {
                    format!(
                        "the first {} candidates and each method's top",
                        super::round::LISTED
                    )
                };
            }
        }
        true
    }

    /// The focused pane's binding for `code`, if it has one: `false` leaves
    /// the key to the global bindings. Moving in the pane's list never claims it.
    fn pane_key(&mut self, code: KeyCode) -> bool {
        // The ordering column's figures take their keys before its table's;
        // the ontology column keeps its own.
        if self.in_order() && self.focus == Focus::Order && self.figure_key(code) {
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
                // The keys that act on a cluster wait for a round.
                let on_cluster = matches!(
                    code,
                    KeyCode::Char('u' | 'k' | ']')
                        | KeyCode::Backspace
                        | KeyCode::Delete
                        | KeyCode::Enter
                );
                if on_cluster && self.selected().is_none() {
                    self.status = "no clusters yet".into();
                    return true;
                }
                if self.candidate_key(code) {
                    return true;
                }
                match code {
                    KeyCode::Char('u') => self.unassign(),
                    KeyCode::Char('k') => self.keep(),
                    KeyCode::Backspace | KeyCode::Delete => self.undo(),
                    KeyCode::Char(']') => self.next_flagged(),
                    KeyCode::Enter => self.open_choice(),
                    // Annotate again: the form, and from the order view a
                    // choice to run the trajectory after.
                    KeyCode::Char('A') => {
                        let order = self.in_order();
                        self.open_form(if order {
                            AfterPass::Offer
                        } else {
                            AfterPass::Nothing
                        });
                    }
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
                    KeyCode::Char('h') => {
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
            Focus::Candidates => {
                let n = self.listed_candidates().len();
                // A save or rescoring may have left the list shorter.
                self.cand_sel = self.cand_sel.min(n.saturating_sub(1));
                move_in(&mut self.cand_sel, n);
                match code {
                    KeyCode::Enter if self.selected().is_some() => {
                        self.pick_candidate(self.cand_sel);
                    }
                    KeyCode::Esc => self.focus = Focus::Clusters,
                    _ => return self.candidate_key(code),
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
                KeyCode::Char('T') => self.toggle_order(false),
                // `t` restyles the figures' labels (taken above when shown).
                KeyCode::Char('t') => self.status = "v figures · T back to the tree".into(),
                KeyCode::Char(' ') => self.toggle_tree_mark(),
                _ => {
                    let took = self.order_key(code);
                    self.sync_tree_to_order();
                    return took;
                }
            },
            Focus::Tree if self.in_order() => {
                return self.order_tree_key(code);
            }
            // The cell map takes its keys first; T still leaves for the order.
            Focus::Tree
                if matches!(self.tree_mode, TreeMode::Cells) && code != KeyCode::Char('T') =>
            {
                return self.cell_map_key(code);
            }
            Focus::Tree if self.trajectory_tui && matches!(code, KeyCode::Char('t' | 'T')) => {
                self.toggle_order(false);
            }
            Focus::Tree
                if code == KeyCode::Char('v') && matches!(self.tree_mode, TreeMode::Panel) =>
            {
                if let Some(m) = &mut self.cell_map {
                    m.exports.disarm();
                    self.map_off = false;
                    self.tree_mode = TreeMode::Cells;
                    self.status = "the round's cells · V shows the tree".into();
                } else {
                    self.status = self
                        .map_hint
                        .clone()
                        .unwrap_or_else(|| "no layout to draw the cells on".into());
                }
            }
            Focus::Tree if code == KeyCode::Char('o') => self.toggle_ontology(),
            Focus::Tree
                if code == KeyCode::Char('c') && matches!(self.tree_mode, TreeMode::Panel) =>
            {
                self.suggest_terms();
            }
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
                        if let Some(p) = self.tree.left(i) {
                            self.select_node(p);
                        }
                    }
                    // A folded branch unfolds; a leaf, or a branch already
                    // open, goes on into the Cell Ontology, beyond the panel.
                    (KeyCode::Right, Some(i)) => {
                        if self.tree.is_folded(i) {
                            self.tree.fold(i, false);
                        } else {
                            self.go_on_into_ontology(i);
                        }
                    }
                    (KeyCode::Esc, _) => self.focus = Focus::Clusters,
                    _ => return false,
                }
            }
        }
        true
    }

    /// The Cell Ontology view opened at panel node `i`'s term, else at the
    /// nearest term above it.
    fn go_on_into_ontology(&mut self, i: usize) {
        let label = self.tree.label(i).to_string();
        let own = self.node_term(i);
        let term = own.clone().or_else(|| {
            self.tree
                .path(i)
                .into_iter()
                .rev()
                .skip(1)
                .find_map(|a| self.node_term(a))
        });
        let Some(mut v) = self.ontology_at(term.clone()) else {
            return;
        };
        let Some(cl) = &self.cl else { return };
        // Without a term here or above, the view's own focus (the root).
        let id = term.unwrap_or_else(|| v.focus.clone());
        v.open_at(cl, &id);
        let name = cl.name(&id).unwrap_or(&id).to_string();
        self.status = if own.is_some() {
            format!("{label} in the Cell Ontology: → opens a term, h hides the unmatched, o back to the panel")
        } else {
            format!("{label} has no Cell Ontology term: at {name}, the nearest above it (o back to the panel)")
        };
        self.tree_mode = TreeMode::Ontology(Box::new(v));
        self.ontology_from_map = false;
    }

    /// Panel node `i`'s Cell Ontology term: the one it was placed on, else
    /// the term its label names.
    fn node_term(&self, i: usize) -> Option<String> {
        self.tree.nodes[i].cl_id.clone().or_else(|| {
            let cl = self.cl.as_ref()?;
            self.term_of(self.tree.label(i)).filter(|id| cl.has(id))
        })
    }

    /// Switch the tree pane to the order view (which types precede which) and
    /// back; `show_figures` opens on the figures (`lupin trajectory`).
    pub(super) fn toggle_order(&mut self, show_figures: bool) {
        self.tree_marked.clear();
        if self.in_order() {
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
        self.show_trajectory();
        self.status = "order view: space marks a type, then > states the first precedes the second, u that they are unrelated".into();
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

    /// Ask for the output prefix of passes, prefilled with the one in use.
    fn ask_output(&mut self) {
        self.prompt = Some(Prompt {
            title: " output prefix for this run's rounds: ".into(),
            text: self.args.out.to_string(),
            pending: Pending::Output,
        });
    }

    /// Take `out` as the output prefix of passes; its latest round, when it
    /// has one, is opened.
    pub(super) fn set_output(&mut self, out: &str) {
        if self.child.is_some() {
            self.status = "wait for the running job before changing the output".into();
            return;
        }
        if out == self.args.out.as_ref() {
            (self.out_chosen, self.out_typed) = (true, true);
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
        (self.out_chosen, self.out_typed) = (true, true);
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
            self.ask_labels();
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
        self.menu = Some(menu);
    }

    /// A run with no labels: annotate it here, or take a labels file.
    fn ask_labels(&mut self) {
        let mut items = vec![
            (
                "Annotate this run now".into(),
                "a form asks for the marker panel, the output and the clustering settings; the pass, then the trajectory on its labels".into(),
                Action::Annotate,
            ),
            (
                "Use a cell<TAB>type labels file".into(),
                "pick the file in the file browser; the trajectory reads its labels".into(),
                Action::LabelsFile,
            ),
        ];
        if crate::manifest::family::has_round(&self.source) {
            items.insert(
                0,
                (
                    "Open an annotated round of this run".into(),
                    "the run's rounds are beside it: pick one (g), and the trajectory takes its labels".into(),
                    Action::ListRuns,
                ),
            );
        }
        self.show(Menu::new(
            "This run has no cell-type labels yet, and the trajectory needs them. \
             Where should they come from?",
            numbered(items),
            NOT_STARTED,
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
        // Statements on types too small to be nodes order none of them.
        let aside = match &v.set_aside[..] {
            [] => String::new(),
            st => format!(
                " Set aside, naming a type below --min-cells {}: {}.",
                v.min_cells,
                st.iter()
                    .map(|s| format!("{} precedes {}", s.from, s.to))
                    .collect::<Vec<_>>()
                    .join("; ")
            ),
        };
        let question = format!(
            "{}Nothing orders these types yet ({why}).{aside} How should the trajectory get its order?",
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
                    "in the order table: space marks two types, > says the first precedes the second, u that they are unrelated; then r".into(),
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
            NOT_STARTED,
        ));
    }

    /// The menu of types to start from, with their cell counts: the nodes
    /// of the run's prior, as a root below `--min-cells` makes the run fail.
    fn ask_root_type(&mut self) {
        let TreeMode::Order(v) = &self.tree_mode else {
            return;
        };
        let (nodes, small): (Vec<usize>, Vec<usize>) =
            (0..v.types.len()).partition(|&i| v.is_node[i]);
        let items = nodes
            .iter()
            .map(|&i| {
                let (t, n) = &v.types[i];
                (
                    t.clone(),
                    format!("start from {t} ({n} cells), kept for later runs"),
                    Action::Root(t.clone()),
                )
            })
            .collect();
        const NAMED: usize = 8;
        let small: Vec<&str> = small.iter().map(|&i| v.types[i].0.as_str()).collect();
        let left_out = match small.len() {
            0 => String::new(),
            n => format!(
                " Left out, below --min-cells {}: {}{}.",
                v.min_cells,
                small[..n.min(NAMED)].join(", "),
                if n > NAMED {
                    format!(" and {} more", n - NAMED)
                } else {
                    String::new()
                }
            ),
        };
        let mut menu = Menu::new(
            format!("Which type does the trajectory start from?{left_out}"),
            numbered(items),
            NOT_STARTED,
        );
        menu.sel = nodes.iter().position(|&i| i >= v.sel).unwrap_or(0);
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
                self.status = "space marks two types, > says the first precedes the second, u that they are unrelated; then r runs".into();
            }
            Action::PriorFile => self.want_file = Some(FileWant::Prior),
            Action::LabelCl => self.want_file = Some(FileWant::LabelCl),
            Action::LabelsFile => self.want_file = Some(FileWant::Labels),
            Action::Annotate => self.open_form(AfterPass::Run),
            Action::RunTrajectory => self.ask_trajectory_out(),
            Action::NotNow => {
                self.status = "r runs the trajectory when you are ready".into();
            }
            Action::ListRuns => self.list_runs(),
            Action::OpenRun { pick, drop_edits } => self.open_run(pick, drop_edits),
            Action::Label { cluster, label } => {
                if !self.reselect(cluster) {
                    return;
                }
                let pick = self
                    .selected()
                    .and_then(|c| c.candidates.iter().find(|t| t.label == label))
                    .cloned();
                match pick {
                    Some(pick) => self.pick(pick),
                    None => {
                        self.status = format!("{label} is no longer a candidate of K{cluster}");
                    }
                }
            }
            Action::LabelOnTree => self.jump_to_tree(),
            Action::Suggested {
                label,
                id,
                note,
                file,
            } => {
                let name = self.cl.as_ref().and_then(|cl| cl.name(&id)).unwrap_or(&id);
                self.prompt = Some(Prompt {
                    title: format!(
                        " {label} → {name} ({id}): write it to {}? enter: yes (the text is the note) · esc: no ",
                        file.display()
                    ),
                    text: note,
                    pending: Pending::Remember { label, id, file },
                });
            }
            Action::NoneSuggested { label } => {
                self.status = format!("none taken: {label} stays without a Cell Ontology term");
            }
            Action::LabelByName { cluster } => {
                self.prompt = Some(Prompt {
                    title: format!(" label K{cluster} as (any name): "),
                    text: String::new(),
                    pending: Pending::LabelName { cluster },
                });
            }
        }
    }

    /// The manifest on screen: the trajectory whose figures the trajectory
    /// TUI shows (a run made here, until another run is opened), else the
    /// round shown, else the run passes start from.
    pub(super) fn shown_manifest(&self) -> PathBuf {
        if let Some(f) = self.figures.as_ref().filter(|_| self.trajectory_tui) {
            return f.data.manifest.clone();
        }
        self.round
            .as_ref()
            .map_or_else(|| self.source.clone(), |r| r.manifest.clone())
    }

    /// List the run's family (`g`): the run, its rounds and the
    /// trajectories made from them, the one on screen marked.
    pub(super) fn list_runs(&mut self) {
        use crate::manifest::family;
        if self.child.is_some() {
            self.status = "wait for the running job, or stop it with x".into();
            return;
        }
        let shown = self.shown_manifest();
        let members = family::family(&shown);
        if members.is_empty() {
            self.status = format!("no run manifests beside {}", shown.display());
            return;
        }
        let at = family::current(&members, &shown);
        let now = super::gallery::now();
        let items = members
            .into_iter()
            .enumerate()
            .map(|(i, m)| {
                let mark = if Some(i) == at { "● " } else { "  " };
                let label = format!(
                    "{mark}{} · {} · {} · {}",
                    m.name(),
                    m.pick.kind.name(),
                    m.holds,
                    super::gallery::ago(m.modified, now)
                );
                let mut detail = m.pick.kind.detail(self.trajectory_tui).to_string();
                if let Some(of) = m.borrowed.as_ref().filter(|_| self.trajectory_tui) {
                    detail += &format!(
                        "; it carries {of}'s trajectory, on its labels: r runs one for this round"
                    );
                }
                let action = Action::OpenRun {
                    pick: m.pick,
                    drop_edits: false,
                };
                (label, detail, action)
            })
            .collect();
        let mut menu = Menu::new(
            format!(
                "The runs of {} (● on screen). Which should the TUI show?",
                family::stem(&shown)
            ),
            numbered(items),
            "kept the run on screen",
        );
        menu.sel = at.unwrap_or(0);
        menu.wide = true;
        self.show(menu);
    }

    /// Show `pick`, a member of the run's family: a round (or a
    /// trajectory's copy of one) in the cluster panes, a trajectory's
    /// figures, the run with no round. Passes then start from it, as when
    /// the TUI is opened on it.
    pub(super) fn open_run(&mut self, pick: crate::manifest::family::Pick, drop_edits: bool) {
        use crate::manifest::family::{pass_origin, Kind};
        if self.child.is_some() {
            self.status = "wait for the running job, or stop it with x".into();
            return;
        }
        let path = pick.path.as_path();
        let name = file_name(path);
        if !drop_edits && !self.edits.is_empty() {
            let detail =
                "the edits are lost; esc keeps them, and s first saves them as a new round";
            let n = self.edits.len();
            return self.show(Menu::new(
                format!("{n} unsaved edit(s) on the round on screen. Open {name} anyway?"),
                numbered(vec![(
                    format!("Open {name}, dropping the edits"),
                    detail.into(),
                    Action::OpenRun {
                        pick: pick.clone(),
                        drop_edits: true,
                    },
                )]),
                "kept the run on screen (s saves the edits)",
            ));
        }
        // The round shown was made with the settings as they stood before
        // any change: a change made since still waits for a pass.
        let base = self.settings_base.take();
        if pick.labelled {
            self.open(path);
        } else {
            self.set_round(None);
        }
        self.settings_base = base;
        if !crate::manifest::run::same_file(path, &self.source) {
            let (source, out) = pass_origin(path);
            self.source = source;
            // The form's Leiden rows follow the new source: fixed when it
            // names its own clusters.
            match crate::manifest::run::load(&self.source.to_string_lossy()) {
                Ok(l) => self.fixed_clusters = l.manifest.cluster.clusters.is_some(),
                Err(e) => self.push_log(format!("[WARN] {e:#}")),
            }
            // A pass run on the old run chose its output for that run only.
            if !self.out_typed {
                self.args.out = out.into_boxed_str();
                self.out_chosen = false;
            }
            self.args.from = Some(self.source.to_string_lossy().into());
            self.target = annotated_path(&self.source, &self.args.out);
        }
        // Only `lupin trajectory`'s TUI shows a trajectory's figures.
        // Another run's figures (and the edges read from them) never stay.
        let trajectory = pick.kind == Kind::Trajectory && self.trajectory_tui;
        let shown = self.figures.take().is_none_or(|f| f.shown);
        if trajectory {
            let picker = self.picker().clone();
            match super::figure_pane::FigurePane::load(&[path], &picker) {
                Ok(Some(mut v)) => {
                    v.shown = shown;
                    self.figures = Some(v);
                }
                Ok(None) => {}
                Err(e) => self.push_log(format!("[WARN] trajectory figures: {e:#}")),
            }
        }
        if self.in_order() {
            // The order view's edges and texts follow the trajectory now shown.
            self.show_trajectory();
        }
        self.status = match (&self.round, trajectory) {
            (Some(r), true) => format!(
                "{name}: its figures, and its round's {} clusters",
                r.clusters.len()
            ),
            (Some(r), false) => format!("{name}: {} clusters (g lists the runs)", r.clusters.len()),
            (None, true) => format!("{name}: its figures; it carries no round (g lists the runs)"),
            (None, false) => {
                format!("{name}: the run, no annotation (A annotates it, g lists its rounds)")
            }
        };
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
        // The first map layers over the ontology at hand; a second one
        // replaces the first, so everything is read again.
        let layer = self.cl.is_some() && self.args.label_cl.is_none();
        let terms = if layer {
            std::fs::read_to_string(&p)
                .map_err(anyhow::Error::from)
                .map(|text| {
                    let mut cl = self.cl.take();
                    if let Some(cl) = &mut cl {
                        cl.add_aliases(&text, &p);
                    }
                    cl
                })
        } else {
            let dir = self.source.parent().map(Path::to_path_buf);
            crate::manifest::ontology::load(
                dir.as_deref(),
                &self.args.markers,
                self.args.obo.as_deref(),
                Some(&p),
                crate::manifest::data_files::Fetch::Allowed,
            )
            .and_then(crate::manifest::data_files::ClData::into_terms)
        };
        match terms {
            Ok(t) => {
                self.cl = t;
                // Ranked against the terms as they were.
                self.forget_suggestions();
                self.rebuild_tree_keeping(None);
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
                self.child = Some(Running::new(c, Job::Trajectory(manifest)));
                self.trajectory.out = Some(out.into());
                self.status = "running the trajectory…".into();
            }
            Err(e) => self.status = format!("{e:#}"),
        }
    }

    /// The attached terminal's picture protocol, as the screen last asked it
    /// ([`App::repick`]); block characters until it has.
    pub(super) fn picker(&mut self) -> &ratatui_image::picker::Picker {
        self.picker
            .get_or_insert_with(ratatui_image::picker::Picker::halfblocks)
    }

    /// The pictures for the terminal attached now, as `repick` says: kitty's
    /// sent again (they went with the last terminal), each figure's terminal
    /// picture made again from it; with another protocol or cell size, the
    /// figures drawn again too, with the new picker.
    pub fn repick(&mut self, repick: Repick) {
        let (picker, changed) = match repick {
            Repick::Changed(picker) => (picker, true),
            Repick::Resend => (self.picker().clone(), false),
            _ => return,
        };
        for pane in [&mut self.figures, &mut self.cell_map]
            .into_iter()
            .flatten()
        {
            pane.repick(&picker, changed);
        }
        self.picker = Some(picker);
    }

    /// After a pass started from the order view: run the trajectory on the
    /// new labels, or not now.
    fn offer_trajectory_run(&mut self) {
        self.show(Menu::new(
            "The new round is open. Run the trajectory on its labels?",
            numbered(vec![
                (
                    "Run the trajectory on the new labels".into(),
                    "asks for the output prefix, then runs lupin trajectory".into(),
                    Action::RunTrajectory,
                ),
                ("Not now".into(), "r runs it later".into(), Action::NotNow),
            ]),
            "r runs the trajectory when you are ready",
        ));
    }

    /// A key for the file browser; a choice or a cancel closes it and goes
    /// to what asked for the file.
    fn browser_key(&mut self, code: KeyCode) {
        let Some((want, b)) = &mut self.browser else {
            return;
        };
        let want = *want;
        let picked = match b.key(code) {
            super::picker::Step::Stay => return,
            super::picker::Step::Chosen(p) => Some(p),
            super::picker::Step::Cancelled => None,
        };
        self.browser = None;
        let picked = picked.as_deref();
        match want {
            FileWant::Markers => self.set_markers(picked),
            FileWant::Labels => self.set_labels(picked),
            FileWant::Prior => self.set_prior(picked),
            FileWant::LabelCl => self.set_label_cl(picked),
        }
    }

    /// Take `path` as the marker panel (picked in the file browser from the
    /// annotation form, which stays open to run).
    pub fn set_markers(&mut self, path: Option<&Path>) {
        let Some(path) = path else {
            self.status = if self.args.markers.is_empty() {
                "no marker panel picked".into()
            } else {
                format!(
                    "kept the marker panel {}",
                    file_name(Path::new(self.args.markers.as_ref()))
                )
            };
            return;
        };
        let p = path.to_string_lossy();
        let read = crate::annotate::markers::read_panel(&p)
            .and_then(|panel| Ok((panel, super::round::panel_sets(&p)?)));
        match read {
            Ok((panel, sets)) => {
                self.panel = panel;
                self.forget_suggestions();
                self.rebuild_tree();
                self.original = sets;
                self.status = format!("marker panel: {}", file_name(path));
                self.args.markers = p.into();
                if let Some(f) = &mut self.form {
                    f.note = None;
                    f.confirm = false;
                    f.row = setting_row(&self.args, Setting::Run);
                }
            }
            Err(e) => self.status = format!("{e:#}"),
        }
    }

    /// Take `path` (picked in the file browser) as the labels of the
    /// trajectory run, then ask for its output prefix.
    pub fn set_labels(&mut self, path: Option<&Path>) {
        match path {
            Some(p) => {
                self.trajectory.labels = Some(p.to_string_lossy().into_owned());
                if self.in_order() {
                    self.reload_order();
                }
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
        if let Some(f) = &mut self.figures {
            f.shown = true;
        }
        if self.in_order() {
            // Its manifest is written anew, maybe over the one shown: its
            // edges and texts are read again.
            self.show_trajectory();
        }
        // The header names the manifest whose figures are up.
        self.status = if self.figures.is_some() && self.trajectory_tui {
            format!("trajectory done in {secs}s: its figures (g lists the runs)")
        } else {
            format!("trajectory done in {secs}s: {}", file_name(manifest))
        };
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
            // types with their cells, else the labels file's, else the
            // panel's types.
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
                None => match self.labels_file_types() {
                    Some(types) => types,
                    None => self.original.keys().map(|k| (k.clone(), 0)).collect(),
                },
            },
        }
    }

    /// The trajectory's labels file's types with their cells.
    fn labels_file_types(&self) -> Option<Vec<(String, usize)>> {
        let path = self.trajectory.labels.as_deref()?;
        let labels = match crate::cell_labels::read_cell_labels(path) {
            Ok(l) => l,
            Err(e) => {
                log::warn!("the labels file: {e:#}");
                return None;
            }
        };
        let mut n: std::collections::BTreeMap<String, usize> = Default::default();
        for (_, t) in labels.iter() {
            *n.entry(t.to_string()).or_default() += 1;
        }
        Some(
            n.into_iter()
                .filter(|(l, _)| l != UNASSIGNED_LABEL)
                .collect(),
        )
    }

    /// (Re)build the order view from the labels on screen, the data files
    /// and the edges of the trajectory shown, already read.
    fn reload_order(&mut self) {
        let (sel, ontology) = match std::mem::replace(&mut self.tree_mode, TreeMode::Panel) {
            TreeMode::Order(v) => (v.sel, v.ontology),
            _ => (0, None),
        };
        let edges = (self.order_shown.as_ref())
            .and_then(|s| s.run.as_ref())
            .map(|r| r.edges.clone())
            .unwrap_or_default();
        let mut v = super::order::OrderView::load(
            self.order_types(),
            self.cl.as_ref(),
            &self.data_search,
            edges,
            self.trajectory.min_cells(),
        );
        v.sel = sel.min(v.types.len().saturating_sub(1));
        if let Some(e) = &v.error {
            self.status = format!("the prior is not a DAG: {e}");
        } else if !v.problems.is_empty() {
            self.status = format!("unreadable precedence file: {}", v.problems.join("; "));
        }
        v.ontology = ontology;
        self.tree_mode = TreeMode::Order(Box::new(v));
        self.sync_tree_to_order();
    }

    /// Keys of the ontology beside the precedence table: the order view's
    /// types on the Cell Ontology, or the full ontology (`o`). Moving onto a
    /// type selects it in the table too.
    fn order_tree_key(&mut self, code: KeyCode) -> bool {
        let (TreeMode::Order(v), Some(cl)) = (&mut self.tree_mode, &self.cl) else {
            return false;
        };
        if let Some(ov) = &mut v.ontology {
            match ov.key(cl, code) {
                ViewKey::Taken => {}
                ViewKey::NoData => {
                    self.status = "no type of the data sits on a Cell Ontology term".into();
                }
                ViewKey::Other => match code {
                    KeyCode::Char('o') => v.ontology = None,
                    KeyCode::Char('T') => self.toggle_order(false),
                    KeyCode::Esc => self.focus = Focus::Clusters,
                    _ => return false,
                },
            }
            return true;
        }
        let visible = v.tree.visible();
        if step(&mut v.tree_sel, visible.len(), code) {
            self.sync_order_to_tree();
            return true;
        }
        let at = visible.get(v.tree_sel).copied();
        match (code, at) {
            (KeyCode::Left, Some(i)) => {
                if let Some(p) = v.tree.left(i) {
                    if let Some(k) = v.tree.visible().iter().position(|&j| j == p) {
                        v.tree_sel = k;
                    }
                }
                self.sync_order_to_tree();
            }
            (KeyCode::Right, Some(i)) => v.tree.fold(i, false),
            (KeyCode::Char(' '), Some(_)) => {
                // Only the order view's types can be ordered.
                if self.sync_order_to_tree() {
                    self.toggle_tree_mark();
                } else {
                    self.status = "not one of the run's types: only they can be ordered".into();
                }
            }
            (KeyCode::Char('o'), at) => {
                let term = at.and_then(|i| v.tree.nodes[i].cl_id.clone());
                let view = self.ontology_on(term);
                if let TreeMode::Order(v) = &mut self.tree_mode {
                    v.ontology = view;
                }
            }
            (KeyCode::Char('T'), _) => self.toggle_order(false),
            (KeyCode::Esc, _) => self.focus = Focus::Clusters,
            _ => return false,
        }
        true
    }

    /// Select in the precedence table the type selected in the ontology
    /// beside it; `false` when that node is not one of the table's types.
    fn sync_order_to_tree(&mut self) -> bool {
        let TreeMode::Order(v) = &mut self.tree_mode else {
            return false;
        };
        let Some(&i) = v.tree.visible().get(v.tree_sel) else {
            return false;
        };
        let l = label_key(v.tree.label(i));
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
        let TreeMode::Order(v) = &mut self.tree_mode else {
            return;
        };
        let Some(i) = v.selected().and_then(|t| v.tree.node_of(t)) else {
            return;
        };
        v.tree.reveal(i);
        if let Some(k) = v.tree.visible().iter().position(|&j| j == i) {
            v.tree_sel = k;
        }
    }

    /// Keys of the order view: move, `>` / `u` on two marked types (the
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
            KeyCode::Char('u') => Relation::Unrelated,
            // `+` and `-` zoom the figures (taken above when shown), as on
            // the cell map.
            KeyCode::Char('+' | '=' | '-' | '_') => {
                self.status = "+ - zoom the figures (v shows them) · u unrelated".into();
                return true;
            }
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
            self.status = "mark exactly two types with space, in order, then > or u".into();
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
        // The table is hidden behind the figures: its marks wait for `V`.
        if v.shown && v.grid.is_none() && code == KeyCode::Char(' ') {
            v.exports.disarm();
            self.status = "V shows the order table to mark types".into();
            return true;
        }
        let k = v.key(code);
        self.took_figure_key(k, false)
    }

    /// Keys of the round's cell map, in the tree pane: `V` shows the tree
    /// again, `o` the Cell Ontology, Esc the clusters; the rest are the
    /// figure pane's.
    fn cell_map_key(&mut self, code: KeyCode) -> bool {
        let Some(v) = &mut self.cell_map else {
            return false;
        };
        if code == KeyCode::Char('V') {
            v.exports.disarm();
            v.grid = None;
            self.map_off = true;
            self.tree_mode = TreeMode::Panel;
            self.status = "the panel's tree · v shows the cells again".into();
            return true;
        }
        let k = v.key(code);
        if self.took_figure_key(k, true) {
            return true;
        }
        match code {
            KeyCode::Char('o') => self.toggle_ontology(),
            KeyCode::Esc => self.focus = Focus::Clusters,
            _ => return false,
        }
        true
    }

    /// Act on what a figure pane's key did, the round's cell map's when
    /// `cells`; whether it took the key.
    fn took_figure_key(&mut self, k: FigureKey, cells: bool) -> bool {
        match k {
            FigureKey::Ignored => return false,
            FigureKey::Took => {}
            FigureKey::Say(s) => self.status = s,
            FigureKey::Relocate(pdf) => {
                self.prompt = Some(Prompt {
                    title: " move the export to (base name, no extension): ".into(),
                    text: pdf.with_extension("").to_string_lossy().into_owned(),
                    pending: Pending::Relocate { pdf, cells },
                });
            }
        }
        true
    }

    /// Switch the tree pane between the panel's tree and the Cell Ontology,
    /// the ontology opening on the selected node's term, else the cluster's.
    fn toggle_ontology(&mut self) {
        if matches!(self.tree_mode, TreeMode::Ontology(_)) {
            let map = std::mem::take(&mut self.ontology_from_map) && self.cell_map.is_some();
            self.tree_mode = if map {
                TreeMode::Cells
            } else {
                TreeMode::Panel
            };
            return;
        }
        let from_map = matches!(self.tree_mode, TreeMode::Cells);
        let from_node = self
            .tree
            .visible()
            .get(self.tree_sel)
            .and_then(|&i| self.node_term(i));
        let from_cluster = self.label_or_top().and_then(|l| self.term_of(&l));
        if let Some(v) = self.ontology_on(from_node.or(from_cluster)) {
            self.tree_mode = TreeMode::Ontology(Box::new(v));
            self.ontology_from_map = from_map;
        }
    }

    /// The Cell Ontology view, "in the data", opened on `term` (else the
    /// root); `None`, said in the status, without an ontology at hand.
    fn ontology_at(&mut self, term: Option<String>) -> Option<super::ontology::OntologyView> {
        let Some(cl) = &self.cl else {
            self.status = "no Cell Ontology at hand (see `lupin data where`)".into();
            return None;
        };
        let focus = term
            .filter(|id| cl.has(id))
            .unwrap_or_else(|| ROOT_TERM.to_string());
        let data = self.ontology_data();
        let cl = self.cl.as_ref()?;
        Some(super::ontology::OntologyView::in_data(cl, data, &focus))
    }

    /// [`Self::ontology_at`] with the cursor on `term` (`o`): a term the
    /// data's tree does not list opens the whole ontology around it.
    fn ontology_on(&mut self, term: Option<String>) -> Option<super::ontology::OntologyView> {
        let mut v = self.ontology_at(term.clone())?;
        if let (Some(id), Some(cl)) = (term.filter(|id| v.focus == *id), &self.cl) {
            v.select(cl, &id);
        }
        Some(v)
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
            if let Some(id) = self.term_of(&t).filter(|id| cl.has(id)) {
                data.entry(id).or_default().push((t, n));
            }
        }
        data
    }

    fn ontology_key(&mut self, code: KeyCode) -> bool {
        let (TreeMode::Ontology(v), Some(cl)) = (&mut self.tree_mode, &self.cl) else {
            return false;
        };
        match v.key(cl, code) {
            ViewKey::Taken => return true,
            ViewKey::NoData => {
                self.status = "no type of the data sits on a Cell Ontology term".into();
                return true;
            }
            ViewKey::Other => {}
        }
        match code {
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
        let term_of = |l: &str| self.term_of(l);
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

    /// `c` on the panel tree: the selected type's suggested Cell Ontology
    /// terms when it has no term of its own, ranked first by a job (for
    /// every such type) when they are not yet.
    fn suggest_terms(&mut self) {
        let Some(i) = self.tree.visible().get(self.tree_sel).copied() else {
            return;
        };
        let label = self.tree.label(i).to_string();
        let Some(cl) = &self.cl else {
            self.status = "no Cell Ontology at hand: nothing to suggest terms from".into();
            return;
        };
        if !super::suggest::needs_term(cl, &self.tree, &label) {
            let what = if self.tree.is_type(&label) {
                "has its own Cell Ontology term"
            } else {
                "is a group, not a panel type"
            };
            self.status = format!("{label} {what}: c suggests terms for a type without one");
            return;
        }
        let labels = super::suggest::unmatched(cl, &self.tree);
        if self.offer_suggestions(&label) {
            return;
        }
        if self.ranking.is_some() {
            self.status =
                "the terms are being ranked: c offers them once they are (x stops it)".into();
            return;
        }
        let (tx, rx) = channel();
        let started = self.job_file("suggest.json").and_then(|file| {
            let obo = self.args.obo.as_deref();
            runner::spawn_suggest(&self.source, obo, &labels, &file, tx).map(|c| (c, file))
        });
        match started {
            Ok((c, file)) => {
                self.push_log(format!(
                    "── suggestions: lupin data suggest-cl for {} type(s) without a term",
                    labels.len()
                ));
                self.ranking = Some(Ranking::new(c, rx, file, label));
                self.status = format!("{RANKING_STATUS} (the first time fetches the text model and encodes the ontology: a minute or two)");
            }
            Err(e) => self.status = format!("{e:#}"),
        }
    }

    /// A file for a job's child to read or write, in the app's own folder.
    pub fn job_file(&mut self, name: &str) -> anyhow::Result<PathBuf> {
        job_file(&mut self.jobs_dir, name)
    }

    /// Offer `label`'s suggested terms, when they are ranked and it still
    /// has no term; `false` when they are not.
    fn offer_suggestions(&mut self, label: &str) -> bool {
        let Some(cl) = &self.cl else { return false };
        if !super::suggest::needs_term(cl, &self.tree, label) {
            return false;
        }
        let Some(s) = self.suggestions.as_ref().and_then(|s| s.of(label)) else {
            return false;
        };
        match self.data_search.user_aliases() {
            _ if s.is_empty() => self.status = format!("no term near {label}"),
            Some(file) => self.menu = Some(super::suggest::menu(label, s, &file)),
            None => self.status = "no user or project folder to write an alias to".into(),
        }
        true
    }

    /// A ranking job wrote `file`: keep its suggestions, and offer those of
    /// the type it was started on while the user waits on the panel tree;
    /// else the next `c` does, so no menu takes keys meant for elsewhere.
    fn suggestions_done(&mut self, file: &Path, label: &str, secs: u64) {
        match super::suggest::Suggestions::read(file) {
            Ok(s) => {
                self.suggestions = Some(s);
                self.status = format!(
                    "Cell Ontology terms ranked in {secs}s: c on a type without a term offers them"
                );
                if !self.idle_on_panel_tree() {
                    self.status = format!("suggestions for {label} ready: press c on it");
                } else {
                    self.offer_suggestions(label);
                }
            }
            Err(e) => self.status = format!("no suggestions: {e:#}"),
        }
    }

    /// Where the trajectory shown is looked for: the figures' manifest,
    /// then the order view's.
    fn trajectory_candidates(&self) -> Vec<PathBuf> {
        let figures = self.figures.as_ref().map(|f| f.data.manifest.clone());
        figures.into_iter().chain(self.manifests()).collect()
    }

    /// Show in the order view the trajectory found now among
    /// [`Self::trajectory_candidates`]: its edges in the table, its texts in
    /// the strip, read together.
    fn show_trajectory(&mut self) {
        let among = self.trajectory_candidates();
        self.order_shown = Some(super::order::Shown::find(among));
        self.reload_order();
    }

    /// In the order view, keep the table and the strip on the trajectory
    /// shown, whatever changed it: when where it is looked for has changed,
    /// it is found again; else, looked at every [`super::order::POLL_EVERY`],
    /// once a file it was read from changed.
    fn keep_order_current(&mut self) {
        if !self.in_order() {
            return;
        }
        let among = self.trajectory_candidates();
        match &mut self.order_shown {
            Some(s) if s.among == among => {
                if s.due() && s.refresh() {
                    self.reload_order();
                }
            }
            _ => self.show_trajectory(),
        }
    }

    /// The order view's strip: the selected type's line while the table has
    /// the focus, else (or for a type the run does not order) the run's
    /// paragraph.
    #[must_use]
    pub fn order_line(&self) -> Option<&str> {
        let t = self.order_shown.as_ref()?.run.as_ref()?;
        let ty = match &self.tree_mode {
            TreeMode::Order(v) if self.focus == Focus::Order => {
                v.types.get(v.sel).map(|(t, _)| t.as_str())
            }
            _ => None,
        };
        t.line(ty)
    }

    /// Take alias `label → id`, just written to `file`, in this session:
    /// the ontology matches by it and the panel's tree is built again, the
    /// cursor on `label`. Unless a table read after `file` (a later layer,
    /// or `--label-cl`) maps `label` too: it wins on every start, so the
    /// session keeps its reading, and that table comes back.
    fn learn_alias(&mut self, label: &str, id: &str, file: &Path) -> Option<PathBuf> {
        // As the ontology's loading reads it: as given, else beside the run.
        let run = (self.args.label_cl.as_deref()).and_then(|f| self.data_search.run_file(f));
        let cl = self.cl.as_ref();
        let later = (self.data_search.alias_layers_after(file).into_iter())
            .chain(run)
            .filter(|p| p.as_path() != file)
            .find(|p| {
                let mut a = crate::annotate::cl_rules::Aliases::default();
                let text = std::fs::read_to_string(p).unwrap_or_default();
                a.add_tsv(&text, &p.display().to_string());
                // As the ontology compares labels: any spelling of its signs.
                match cl {
                    Some(cl) => a.iter().any(|(l, _)| cl.alike(l, label)),
                    None => a.get(label).is_some(),
                }
            });
        if later.is_some() {
            return later;
        }
        let cl = self.cl.as_mut()?;
        cl.add_aliases(&format!("{label}\t{id}\n"), &file.display().to_string());
        self.rebuild_tree_keeping(Some(label));
        None
    }

    /// [`Self::rebuild_tree`], keeping what was folded, and the selection
    /// on node `on`, else on the one selected before, shown.
    pub(super) fn rebuild_tree_keeping(&mut self, on: Option<&str>) {
        let tree = &self.tree;
        let folded: Vec<String> = (0..tree.nodes.len())
            .filter(|&i| tree.is_folded(i))
            .map(|i| tree.label(i).to_string())
            .collect();
        let selected = (tree.visible().get(self.tree_sel)).map(|&i| tree.label(i).to_string());
        self.rebuild_tree();
        for l in &folded {
            if let Some(i) = self.tree.node_of(l) {
                self.tree.fold(i, true);
            }
        }
        let node = |l: Option<&str>| l.and_then(|l| self.tree.node_of(l));
        if let Some(i) = node(on).or_else(|| node(selected.as_deref())) {
            self.tree.reveal(i);
            self.select_node(i);
        }
        self.clamp_tree_sel();
    }

    /// Build the panel's tree again from [`Self::panel`], with each type's
    /// terms; the selection stays on a row of the new tree.
    pub(super) fn rebuild_tree(&mut self) {
        self.tree = crate::manifest::ontology::panel_tree_on(self.cl.as_ref(), &self.panel);
        self.clamp_tree_sel();
        self.terms_of.borrow_mut().clear();
        if let Some(cl) = &self.cl {
            self.panel_ancestry = super::ontology::type_ancestry(cl, &self.tree);
        }
    }

    /// Keep the panel tree's selection on a row it shows.
    fn clamp_tree_sel(&mut self) {
        self.tree_sel = (self.tree_sel).min(self.tree.visible().len().saturating_sub(1));
    }

    /// The Cell Ontology term `label` names ([`super::ontology::term_of`]);
    /// none without an ontology.
    pub fn term_of(&self, label: &str) -> Option<String> {
        let cl = self.cl.as_ref()?;
        if let Some(t) = self.terms_of.borrow().get(label) {
            return t.clone();
        }
        let t = super::ontology::term_of(cl, &self.tree, label);
        (self.terms_of.borrow_mut()).insert(label.to_string(), t.clone());
        t
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
        // Only the panel tree on screen is moved, the cell map giving way to
        // it: another view leaves its selection where the user left it.
        if let Some(back) = match self.tree_mode {
            TreeMode::Panel => None,
            TreeMode::Cells => {
                self.tree_mode = TreeMode::Panel;
                None
            }
            TreeMode::Ontology(_) => Some("o on the tree"),
            TreeMode::Order(_) => Some("T"),
        } {
            self.status = format!("the panel tree is not shown: {back} returns to it");
        } else if let Some(i) = self.label_or_top().and_then(|l| self.tree.node_of(&l)) {
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
        if let Some(mut r) = self.child.take() {
            if !r.child.saves() {
                r.child.kill();
            }
            let _ = r.child.wait();
        }
        if let Some(r) = self.ranking.take() {
            r.stop();
        }
        self.stop_rescoring();
    }
}

/// `name` in the app's jobs folder, made on first use: a name in the shared
/// temporary folder could be taken first by someone else, or a link there
/// send the write elsewhere.
fn job_file(dir: &mut Option<tempfile::TempDir>, name: &str) -> anyhow::Result<PathBuf> {
    let dir = match dir {
        Some(d) => d,
        empty => {
            let mut b = tempfile::Builder::new();
            b.prefix("lupin-tui-");
            #[cfg(unix)]
            b.permissions(std::os::unix::fs::PermissionsExt::from_mode(0o700));
            empty.insert(b.tempdir().context("making a folder for the jobs' files")?)
        }
    };
    Ok(dir.path().join(name))
}

/// A round's clusters and how many need a look, as the status says them.
fn round_counts(r: &RoundView) -> String {
    let flagged = r.clusters.iter().filter(|c| c.flagged()).count();
    format!(
        "{} clusters, {flagged} to look at (]: next)",
        r.clusters.len()
    )
}

/// A pass's arguments but its run and output, which are what it reads and
/// where it writes, not how.
fn pass_settings(args: &AnnotateCliArgs) -> Vec<String> {
    let argv = args.to_argv();
    let mut kept = Vec::with_capacity(argv.len());
    let mut it = argv.into_iter();
    while let Some(a) = it.next() {
        if a == "--out" || a == "--from" {
            it.next();
        } else {
            kept.push(a);
        }
    }
    kept
}

/// The status line when a trajectory menu is closed without a choice.
const NOT_STARTED: &str = "no trajectory run started";

/// The labels file a trajectory run on `loaded` recorded in its settings,
/// resolved against the manifest's directory, when it is still there.
pub(super) fn recorded_labels(loaded: &crate::manifest::run::Loaded) -> Option<String> {
    let settings = loaded.manifest.trajectory.settings.as_ref()?;
    loaded.recorded(settings.get("labels")?.as_str()?)
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

/// A candidate's evidence for a reason: share, SuSiE's PIP and fold when it
/// made the call, and NES and q when known.
fn evidence(c: &super::round::Candidate) -> String {
    let mut s = format!("share {:.2}", c.share);
    if let Some(e) = c.susie() {
        s += &format!(", PIP {:.2}, ×{:.1}", e.pip, e.fold);
    }
    if let Some(nes) = c.nes {
        s += &format!(", NES {nes:.2}");
    }
    if let Some(q) = c.q {
        s += &format!(", q {q:.3}");
    }
    if let Some(score) = c.score {
        s += &format!(", score {score:.2} of the cells");
    }
    s
}

pub(super) use crate::manifest::family::file_name;

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
