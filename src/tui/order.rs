//! The order view: which cell types precede which, as the prior lupin would
//! build for the labels on screen (`docs/trajectory-plan.md` §6). The user
//! marks two types and states that one precedes the other, or that they are
//! unrelated; the statement goes to the project's `precedence.tsv`.

use crate::annotate::celltype_tree::ClTerms;
use crate::annotate::panel_tree::PanelTree;
use crate::describe::trajectory::{self, OrderFacts, OrderText};
use crate::manifest::data_files::{append_row, SearchPath, PRECEDENCE};
use crate::manifest::run::load;
use crate::trajectory::edges::{self, EdgeRow, Verdict};
use crate::trajectory::lineages;
use crate::trajectory::prior::{self, Prior, Relation, Source, Statement};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

/// How the order view runs `lupin trajectory`: the options it carries over
/// and the output prefix last used.
#[derive(Debug, Default)]
pub struct TrajectoryRun {
    pub argv: Vec<String>,
    pub out: Option<String>,
    /// A labels file, given as `--labels` or picked in the TUI, for a run
    /// with no annotation.
    pub labels: Option<String>,
}

impl TrajectoryRun {
    /// The `--min-cells` the run is given: types with fewer cells are not
    /// nodes of its prior. Without one, the run's default.
    #[must_use]
    pub fn min_cells(&self) -> usize {
        self.argv
            .windows(2)
            .rev()
            .find(|w| w[0] == "--min-cells")
            .and_then(|w| w[1].parse().ok())
            .unwrap_or(crate::trajectory::run::DEFAULT_MIN_CELLS)
    }

    /// How the TUI runs `lupin trajectory` with `t`'s options: `--labels`
    /// is kept apart, so the order view lists its types and a round on
    /// screen stands in for it.
    #[must_use]
    pub fn of(t: &crate::trajectory::run::TrajectoryArgs) -> Self {
        let mut argv = t.child_argv();
        if let Some(i) = argv.iter().position(|a| a == "--labels") {
            argv.drain(i..(i + 2).min(argv.len()));
        }
        Self {
            argv,
            out: t.out.as_deref().map(str::to_string),
            labels: t.labels.as_deref().map(str::to_string),
        }
    }
}

/// The first `{stem}.T{k}` (k = 1, 2, …) beside `source` with no manifest yet.
pub fn default_out(source: &Path) -> String {
    crate::manifest::family::Tag::Trajectory.next_free(source)
}

/// A direct edge of the prior as the table shows it.
pub struct EdgeView<'a> {
    pub from: String,
    pub to: String,
    /// Where the statement came from, or the types it goes through.
    pub source: String,
    /// The trajectory run's row for the pair, when it has one.
    pub data: Option<&'a EdgeRow>,
}

pub struct OrderView {
    /// The labels on screen with their cell counts.
    pub types: Vec<(String, usize)>,
    /// The combined statements about them.
    pub statements: Vec<Statement>,
    /// Which types are nodes of the prior: at least the run's
    /// `--min-cells`, as `lupin trajectory` counts them.
    pub is_node: Vec<bool>,
    /// The run's `--min-cells`.
    pub min_cells: usize,
    /// The `precedes` statements naming a type below `min_cells`, which
    /// order no pair of nodes by themselves.
    pub set_aside: Vec<Statement>,
    /// The prior they make; `None` when they do not form a DAG (see `error`).
    pub prior: Option<Prior>,
    /// Why the statements do not form a DAG.
    pub error: Option<String>,
    /// Precedence files that could not be read.
    pub problems: Vec<String>,
    /// The latest trajectory run's pairs, with their verdicts.
    pub edges: Vec<EdgeRow>,
    /// The selected type.
    pub sel: usize,
    /// Where a statement is written: the project's `precedence.tsv`.
    pub file: Option<PathBuf>,
    /// The types on the Cell Ontology, as annotate's tree pane draws them,
    /// and the row selected in it.
    pub tree: PanelTree,
    pub tree_sel: usize,
    /// The full ontology in place of `tree` (`o`).
    pub ontology: Option<super::ontology::OntologyView>,
    /// How many types map to a term, and the ontology's develops-from
    /// links among them.
    pub mapped: usize,
    pub links: usize,
}

impl OrderView {
    /// The view for `types`: the Cell Ontology's statements about them and the
    /// user's and project's `precedence.tsv`, combined and reduced, with the
    /// `edges` of a trajectory run. As `lupin trajectory` builds its prior,
    /// a type with fewer than `min_cells` cells is no node (a count of 0 is
    /// unknown, and counts as one): the Cell Ontology speaks of the nodes
    /// only, and edges join nodes.
    pub fn load(
        types: Vec<(String, usize)>,
        cl: Option<&ClTerms>,
        search: &SearchPath,
        edges: Vec<EdgeRow>,
        min_cells: usize,
    ) -> Self {
        let names: Vec<Box<str>> = types.iter().map(|(t, _)| t.as_str().into()).collect();
        let refs: Vec<&str> = names.iter().map(AsRef::as_ref).collect();
        let is_node: Vec<bool> = types
            .iter()
            .map(|&(_, n)| n == 0 || n >= min_cells)
            .collect();
        let (nodes, small): (Vec<_>, Vec<_>) = refs.iter().zip(&is_node).partition(|(_, &n)| n);
        let nodes: Vec<&str> = nodes.into_iter().map(|(t, _)| *t).collect();
        let small: std::collections::HashSet<&str> = small.into_iter().map(|(t, _)| *t).collect();
        let (layers, problems) = prior::standing_layers(cl, &nodes, search);
        let problems: Vec<String> = problems.iter().map(|e| format!("{e:#}")).collect();
        let statements = prior::combine(&layers);
        let (prior, error) = match prior::build(&names, &is_node, statements.clone()) {
            Ok(p) => (Some(p), None),
            Err(e) => (None, Some(e.to_string())),
        };
        // A statement naming a small type is set aside when the nodes'
        // edges are the same without it: a chain through the type (A → s →
        // B) still orders A before B, as in the run.
        let pairs = |p: &Prior| p.edges.iter().map(|e| (e.from, e.to)).collect::<Vec<_>>();
        let set_aside = statements
            .iter()
            .enumerate()
            .filter(|(_, s)| {
                s.relation == Relation::Precedes
                    && (small.contains(s.from.as_str()) || small.contains(s.to.as_str()))
            })
            .filter(|&(i, _)| {
                prior.as_ref().is_none_or(|p| {
                    let mut without = statements.clone();
                    without.remove(i);
                    prior::build(&names, &is_node, without).is_ok_and(|q| pairs(&q) == pairs(p))
                })
            })
            .map(|(_, s)| s.clone())
            .collect();
        let panel: Vec<(String, String)> = types
            .iter()
            .map(|(t, _)| (String::new(), t.clone()))
            .collect();
        let mapped = cl.map_or(0, |cl| cl.map_labels(refs.iter().copied()).0.len());
        let links = statements
            .iter()
            .filter(|s| {
                matches!(s.source, Source::Cl | Source::ClInherited)
                    && s.relation == Relation::Precedes
            })
            .count();
        Self {
            tree: crate::manifest::ontology::panel_tree_on(cl, &panel),
            tree_sel: 0,
            ontology: None,
            mapped,
            links,
            types,
            statements,
            is_node,
            min_cells,
            set_aside,
            prior,
            error,
            problems,
            edges,
            sel: 0,
            file: search.amend(PRECEDENCE),
        }
    }

    pub fn selected(&self) -> Option<&str> {
        self.types.get(self.sel).map(|(t, _)| t.as_str())
    }

    /// The prior's direct edges, with the run's row for each when known.
    pub fn edges(&self) -> Vec<EdgeView<'_>> {
        let Some(p) = &self.prior else {
            return Vec::new();
        };
        p.edges
            .iter()
            .map(|e| {
                let from = self.types[e.from].0.clone();
                let to = self.types[e.to].0.clone();
                let source = e.source.map_or_else(
                    || {
                        format!(
                            "via {}",
                            e.via
                                .iter()
                                .map(|&i| self.types[i].0.as_str())
                                .collect::<Vec<_>>()
                                .join(", ")
                        )
                    },
                    |s| s.as_str().to_string(),
                );
                let data = self.row(&from, &to);
                EdgeView {
                    from,
                    to,
                    source,
                    data,
                }
            })
            .collect()
    }

    /// The run's row for a pair, whichever way it was written.
    pub fn row(&self, a: &str, b: &str) -> Option<&EdgeRow> {
        self.edges.iter().find(|e| e.is_pair(a, b))
    }

    /// Pairs the trajectory run found connected but the prior does not
    /// order, strongest first.
    pub fn candidates(&self) -> Vec<&EdgeRow> {
        self.with_verdict(Verdict::Candidate)
    }

    /// Edges the trajectory run added to place node types the prior left
    /// without one, strongest first.
    pub fn inferred(&self) -> Vec<&EdgeRow> {
        self.with_verdict(Verdict::Inferred)
    }

    /// The run's pairs with verdict `v`, strongest first.
    fn with_verdict(&self, v: Verdict) -> Vec<&EdgeRow> {
        let mut out: Vec<&EdgeRow> = self.edges.iter().filter(|e| e.verdict == Some(v)).collect();
        out.sort_by(|p, q| q.connectivity.total_cmp(&p.connectivity));
        out
    }

    /// The `unrelated` statements, which show no edge.
    pub fn unrelated(&self) -> Vec<&Statement> {
        self.statements
            .iter()
            .filter(|s| s.relation == Relation::Unrelated)
            .collect()
    }

    /// Why `from relation to` should not be written: the file's last
    /// statement about the pair already says it, or the prior it makes
    /// would not build (a cycle, or `unrelated` on a pair the rest orders).
    /// On a prior that is already broken any statement goes, so the one that
    /// broke it can be replaced.
    pub fn refuse(&self, from: &str, to: &str, relation: Relation) -> Option<String> {
        let stated = self.file.as_deref().and_then(|f| {
            let text = std::fs::read_to_string(f).ok()?;
            let st = prior::parse_statements(&text, Source::Project, "").ok()?;
            st.into_iter()
                .rev()
                .find(|s| (s.from == from && s.to == to) || (s.from == to && s.to == from))
        });
        if let Some(s) = stated {
            let same_way = relation == Relation::Unrelated || s.from == from;
            if s.relation == relation && same_way {
                return Some(format!(
                    "{} already says {from} {} {to}",
                    self.file.as_deref().unwrap_or(Path::new("")).display(),
                    relation.as_str()
                ));
            }
        }
        self.prior.as_ref()?;
        let new = Statement {
            from: from.into(),
            to: to.into(),
            relation,
            source: Source::Project,
            note: String::new(),
        };
        let statements = prior::combine(&[self.statements.clone(), vec![new]]);
        let names: Vec<Box<str>> = self.types.iter().map(|(t, _)| t.as_str().into()).collect();
        prior::build(&names, &self.is_node, statements)
            .err()
            .map(|e| format!("not recorded: {e}"))
    }

    /// Append a statement to the project's `precedence.tsv`, unless
    /// [`Self::refuse`] gives a reason not to.
    pub fn record(
        &self,
        from: &str,
        to: &str,
        relation: Relation,
        note: &str,
    ) -> anyhow::Result<PathBuf> {
        if let Some(why) = self.refuse(from, to, relation) {
            anyhow::bail!(why);
        }
        let file = self.file.clone().ok_or_else(|| {
            anyhow::anyhow!("no project or user directory to write precedence.tsv in")
        })?;
        append_row(
            &file,
            "Which cell types precede which, layered over lupin's own data files\n\
             (see `lupin data where`). Columns: from, to, relation (precedes | unrelated), note.",
            &[from, to, relation.as_str(), note],
        )?;
        Ok(file)
    }
}

/// A file's size and last change, to tell when it changed; `None` when it
/// is missing.
type Stamp = Option<(u64, Option<std::time::SystemTime>)>;

fn stamp(file: &Path) -> Stamp {
    let m = std::fs::metadata(file).ok()?;
    Some((m.len(), m.modified().ok()))
}

/// The trajectory the order view shows: its edges for the table and the
/// texts code writes from its outputs for the strip, read together and
/// read again together when the manifest or an output changes.
pub struct ShownRun {
    pub manifest: PathBuf,
    pub edges: Vec<EdgeRow>,
    pub texts: Vec<OrderText>,
    /// The manifest and its outputs, with their stamps as last read.
    files: Vec<(PathBuf, Stamp)>,
}

impl ShownRun {
    /// The first of `manifests` whose `trajectory` section records an edge
    /// table that can be read (its own or one a round borrowed from the run
    /// it was made on). `known` are edges already read from one of them
    /// (the figures'), taken rather than read again. When it finds none,
    /// the files that told so: each manifest, and the edge table one
    /// records.
    pub fn find(
        manifests: &[PathBuf],
        known: Option<(&Path, &[EdgeRow])>,
    ) -> (Option<Self>, Vec<PathBuf>) {
        let mut looked = Vec::new();
        for m in manifests {
            let known = known.filter(|(at, _)| at == m).map(|(_, e)| e);
            match Self::read(m, known) {
                Ok(run) => return (Some(run), Vec::new()),
                Err(files) => looked.extend(files),
            }
        }
        (None, looked)
    }

    /// The run `manifest` records, else the files that told it has none.
    fn read(manifest: &Path, known: Option<&[EdgeRow]>) -> Result<Self, Vec<PathBuf>> {
        let none = |more: Option<&str>| {
            let mut files = vec![manifest.to_path_buf()];
            files.extend(more.map(PathBuf::from));
            files
        };
        let loaded = load(&manifest.to_string_lossy()).map_err(|_| none(None))?;
        let [lineages_file, edges] = OrderFacts::files(&loaded);
        let edges_file = edges.ok_or_else(|| none(None))?;
        let edges = match known {
            Some(e) => e.to_vec(),
            None => edges::read(&edges_file).map_err(|_| none(Some(&edges_file)))?,
        };
        let lineages = match &lineages_file {
            Some(p) => lineages::read(p).unwrap_or_else(|e| {
                log::warn!("no trajectory texts: {e:#}");
                Vec::new()
            }),
            None => Vec::new(),
        };
        let facts = OrderFacts { lineages, edges };
        let files = [
            Some(manifest.to_path_buf()),
            lineages_file.map(PathBuf::from),
            Some(edges_file.into()),
        ]
        .into_iter()
        .flatten()
        .collect();
        Ok(Self {
            manifest: manifest.to_path_buf(),
            texts: trajectory::texts(&facts),
            edges: facts.edges,
            files: stamped(files),
        })
    }

    /// Read it again when its files changed; whether they did. A manifest
    /// or edges that no longer read keep it as it was.
    pub fn poll(&mut self) -> bool {
        if self.files.iter().all(|(f, s)| stamp(f) == *s) {
            return false;
        }
        match Self::read(&self.manifest, None) {
            Ok(again) => *self = again,
            Err(_) => {
                for (f, s) in &mut self.files {
                    *s = stamp(f);
                }
            }
        }
        true
    }

    /// Type `ty`'s line, or the run's paragraph without one.
    #[must_use]
    pub fn line(&self, ty: Option<&str>) -> Option<&str> {
        let of = |s: Option<&str>| {
            self.texts
                .iter()
                .find(|d| d.subject.as_deref() == s)
                .map(|d| d.text.as_str())
        };
        ty.and_then(|t| of(Some(t))).or_else(|| of(None))
    }
}

/// How often the order view looks at its trajectory's files for a change.
pub const POLL_EVERY: Duration = Duration::from_secs(1);

/// The trajectory the order view shows and where it was looked for, and
/// when its files were last looked at.
pub struct Shown {
    pub among: Vec<PathBuf>,
    pub run: Option<ShownRun>,
    /// With no run found, the files that told so, with their stamps then.
    watched: Vec<(PathBuf, Stamp)>,
    polled: Instant,
}

impl Shown {
    /// The run found among `among` ([`ShownRun::find`]).
    #[must_use]
    pub fn find(among: Vec<PathBuf>, known: Option<(&Path, &[EdgeRow])>) -> Self {
        let (run, looked) = ShownRun::find(&among, known);
        Self {
            among,
            run,
            watched: stamped(looked),
            polled: Instant::now(),
        }
    }

    /// [`ShownRun::poll`], at most every [`POLL_EVERY`], or with no run
    /// found yet a look among the candidates again once a file that told so
    /// changed; whether it changed.
    pub fn poll(&mut self) -> bool {
        if self.polled.elapsed() < POLL_EVERY {
            return false;
        }
        self.polled = Instant::now();
        if let Some(run) = &mut self.run {
            return run.poll();
        }
        if self.watched.iter().all(|(f, s)| stamp(f) == *s) {
            return false;
        }
        let (run, looked) = ShownRun::find(&self.among, None);
        self.run = run;
        self.watched = stamped(looked);
        self.run.is_some()
    }
}

/// `files` with their stamps now.
fn stamped(files: Vec<PathBuf>) -> Vec<(PathBuf, Stamp)> {
    files
        .into_iter()
        .map(|f| {
            let s = stamp(&f);
            (f, s)
        })
        .collect()
}

#[cfg(test)]
#[path = "tests/order.rs"]
mod tests;
