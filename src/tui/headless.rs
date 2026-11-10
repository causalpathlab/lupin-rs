//! What the screen builds when it opens a round, built without the screen:
//! the round's marker panel on the Cell Ontology, with the user's and the
//! project's mixed labels, hidden genes and data files. `lupin export` and
//! the MCP's ontology tools work on it as the TUI's panes do.

use super::genes::GeneFilter;
use super::ontology::Mixed;
use crate::annotate::celltype_tree::ClTerms;
use crate::annotate::celltype_tree::TypeTree;
use crate::annotate::markers::label_key;
use crate::annotate::panel_tree::PanelTree;
use crate::manifest::data_files::{append_row, ALIASES_HEADER, GROUPS_HEADER};
use crate::manifest::data_files::{Fetch, SearchPath};
use crate::manifest::run::{self, Loaded};
use anyhow::{Context, Result};
use serde_json::{json, Map, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

/// A round's panel as the TUI holds it.
pub struct PanelState {
    pub loaded: Loaded,
    pub panel: Vec<(String, String)>,
    /// Where its data files are read and the user's answers written.
    pub search: SearchPath,
    /// The Cell Ontology, when one is found.
    pub cl: Option<ClTerms>,
    /// The panel's types on the ontology, else grouped by shared markers.
    pub tree: PanelTree,
    /// The panel's own markers per type, before any edit.
    pub original: BTreeMap<String, BTreeSet<String>>,
    pub mixed: Mixed,
}

impl PanelState {
    /// The panel of round (or run) `from`, its ontology `obo` when given,
    /// else the one its data files find (downloaded once, as the TUI does).
    pub fn load(from: &str, obo: Option<&str>) -> Result<Self> {
        let loaded = run::load(from)?;
        let markers = crate::annotate_cmd::recorded_markers(&loaded).unwrap_or_default();
        let panel = if markers.is_empty() {
            Vec::new()
        } else {
            crate::annotate::markers::read_panel(&markers)?
        };
        // The panel's own markers per type, from the pairs just read.
        let mut original: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
        for (gene, ty) in &panel {
            original
                .entry(label_key(ty))
                .or_default()
                .insert(gene.clone());
        }
        let data = crate::manifest::ontology::load(
            Some(&loaded.dir),
            &markers,
            obo,
            None,
            Fetch::Allowed,
        )?;
        let search = data.search.clone();
        let cl = data.into_terms()?;
        let tree = crate::manifest::ontology::panel_tree_on(cl.as_ref(), &panel);
        Ok(Self {
            mixed: Mixed::load(&search)?,
            loaded,
            panel,
            search,
            cl,
            tree,
            original,
        })
    }

    /// The Cell Ontology, or why there is none.
    pub fn ontology(&self) -> Result<&ClTerms> {
        self.cl
            .as_ref()
            .context("no Cell Ontology found: give `obo`, or run `lupin data fetch` once online")
    }

    /// The round's final annotation as `{out}.cell_annotation.parquet` and
    /// `{out}.celltype_markers.tsv` (`out`: the round's own prefix when
    /// none), as `e` on the screen writes it.
    pub fn export(&self, out: Option<&str>) -> Result<(PathBuf, PathBuf)> {
        let round = super::round::RoundView::load(&self.loaded.file)?;
        let stem = match out {
            Some(o) => o.to_string(),
            None => run::derive_out_prefix(&self.loaded.file.to_string_lossy()),
        };
        super::export::write(
            &round,
            &self.tree,
            self.cl.as_ref(),
            &self.mixed,
            &self.original,
            &stem,
        )
    }
    /// The panel's types, each once, as the panel names them.
    fn types(&self) -> Vec<String> {
        crate::annotate::celltype_tree::panel_types(&self.panel)
    }

    /// Term `id` as a tool shows it.
    fn term(cl: &ClTerms, id: &str) -> Value {
        json!({ "id": id, "name": cl.name(id) })
    }

    fn terms<'a>(cl: &ClTerms, ids: impl IntoIterator<Item = &'a String>) -> Vec<Value> {
        ids.into_iter().map(|id| Self::term(cl, id)).collect()
    }

    /// The terms whose names, synonyms or abbreviations hold `query`
    /// (`/` on the screen), at most `limit`.
    pub fn search(&self, query: &str, limit: usize) -> Result<Value> {
        let cl = self.ontology()?;
        Ok(json!({ "terms": Self::terms(cl, &cl.search(query, limit)) }))
    }

    /// Term `id` and its place: its parents and children (← and → on the
    /// screen), its lineage from the top, what it develops from, and the
    /// panel's types at or under it.
    pub fn term_facts(&self, id: &str) -> Result<Value> {
        let cl = self.ontology()?;
        anyhow::ensure!(
            cl.has(id),
            "{id}: no such Cell Ontology term (`search` finds terms by name)"
        );
        let types: Vec<String> = (self.types().into_iter())
            .filter(|t| {
                super::ontology::term_of(cl, &self.tree, t)
                    .is_some_and(|own| cl.ancestors_or_self(&own).contains(id))
            })
            .collect();
        Ok(json!({
            "id": id,
            "name": cl.name(id),
            "parents": Self::terms(cl, cl.parents(id)),
            "children": Self::terms(cl, cl.children(id)),
            "lineage": Self::terms(cl, &cl.lineage(id)),
            "develops_from": Self::terms(cl, cl.develops_from(id)),
            "types": types,
        }))
    }

    /// The panel on the ontology (`d` on the screen, the terms in the
    /// data): each type's term and lineage, and the types with none.
    pub fn on_ontology(&self) -> Result<Value> {
        let cl = self.ontology()?;
        let (mut on, mut unmatched) = (Map::new(), Vec::new());
        for t in self.types() {
            match super::ontology::term_of(cl, &self.tree, &t) {
                Some(id) => {
                    let lineage = Self::terms(cl, &cl.lineage(&id));
                    on.insert(
                        t,
                        json!({ "term": id, "name": cl.name(&id), "lineage": lineage }),
                    );
                }
                None => unmatched.push(t),
            }
        }
        Ok(json!({ "types": on, "unmatched": unmatched }))
    }

    /// The questions the screen asks of the panel (`Q`): whether a type
    /// left alone (its siblings went under a nearer class) joins them, and
    /// which term a type with none is, each until answered.
    pub fn questions(&self) -> Result<Value> {
        let cl = self.ontology()?;
        let (mut join, mut terms) = (Vec::new(), Vec::new());
        for q in super::ask::queue(cl, &self.tree, &self.panel) {
            match q {
                super::ask::Question::Join(q) => join.push(json!({
                    "type": q.label,
                    "class": Self::term(cl, &q.class),
                    "group": q.group,
                })),
                // As the panel names it, as the other tools do.
                super::ask::Question::Term(key) => {
                    let named = (self.types().into_iter()).find(|t| label_key(t) == key);
                    terms.push(json!({ "type": named.unwrap_or(key) }));
                }
            }
        }
        Ok(json!({ "join": join, "terms": terms }))
    }

    /// Answer the panel's question on type `ty`: `join` or `alone` for a
    /// type left alone, `none` for one no term fits, else the Cell Ontology
    /// term it is. Kept, with `note`, where every later pass reads it (the
    /// user's answers, else the project's). The file written.
    pub fn answer(&self, ty: &str, answer: &str, note: &str) -> Result<PathBuf> {
        let cl = self.ontology()?;
        anyhow::ensure!(!note.trim().is_empty(), "an answer needs a `note`: why");
        anyhow::ensure!(
            self.tree.is_type(ty),
            "{ty} is not one of the panel's types: {}",
            self.types().join(", ")
        );
        let given = answer.trim().to_lowercase();
        // The answers on types go with the groups; a term, to the aliases.
        let (file, header, value) = match given.as_str() {
            "join" | "alone" => {
                let asked = TypeTree::join_questions(cl, &self.panel);
                anyhow::ensure!(
                    asked.iter().any(|q| label_key(&q.label) == label_key(ty)),
                    "{ty} is not a type left alone: join and alone answer whether one joins its siblings"
                );
                (self.search.user_groups(), GROUPS_HEADER, given.as_str())
            }
            "none" => {
                anyhow::ensure!(
                    super::suggest::needs_term(cl, &self.tree, ty),
                    "{ty} has a Cell Ontology term already"
                );
                (self.search.user_groups(), GROUPS_HEADER, "none")
            }
            _ => {
                let id = answer.trim();
                anyhow::ensure!(
                    cl.has(id),
                    "{id}: no such Cell Ontology term (answer join, alone, none or a term's id)"
                );
                (self.search.user_aliases(), ALIASES_HEADER, id)
            }
        };
        let file = file.context("no user or project folder to keep the answer in")?;
        append_row(&file, header, &[ty, value, note])?;
        Ok(file)
    }
}

/// Hide genes matching `patterns` (exact names or `PFX-*`) from the
/// specific genes of the run in `run_dir` from now on, as `h` and `X` on
/// the screen do. The file written.
pub fn hide_genes(run_dir: &Path, patterns: &[String]) -> Result<PathBuf> {
    anyhow::ensure!(
        !patterns.is_empty() && patterns.iter().all(|p| !p.trim().is_empty()),
        "give each gene or pattern (PFX-*) to hide"
    );
    let search = SearchPath::new(Some(run_dir));
    let file = GeneFilter::file(&search)
        .context("no user or project folder to keep the hidden genes in")?;
    let mut hidden = GeneFilter::load(&search)?;
    for p in patterns {
        hidden.add(p, &file)?;
    }
    Ok(file)
}

/// Name the mix of `parts` `name` in the project's mixed labels (else the
/// user's), as `+` on the screen does, so the export places it under their
/// common ancestor; the file written.
pub fn remember_mix(run_dir: &Path, name: &str, parts: &[String]) -> Result<PathBuf> {
    let search = SearchPath::new(Some(run_dir));
    let file =
        Mixed::file(&search).context("no user or project folder to keep the mixed label in")?;
    Mixed::load(&search)?.add(name, parts, &file)?;
    Ok(file)
}

/// `lupin export`: a round's final annotation, for use outside lupin.
#[derive(clap::Args, Debug)]
pub struct ExportArgs {
    #[arg(short, long, help = "A round (its manifest or prefix)")]
    pub from: Box<str>,
    #[arg(
        short,
        long,
        help = "Output prefix: {out}.cell_annotation.parquet and {out}.celltype_markers.tsv (the round's own prefix by default)"
    )]
    pub out: Option<Box<str>>,
    #[arg(long, help = "The Cell Ontology OBO file, as annotate's --obo")]
    pub obo: Option<Box<str>>,
}

pub fn run_export(args: &ExportArgs) -> Result<()> {
    let state = PanelState::load(&args.from, args.obo.as_deref())?;
    anyhow::ensure!(
        state.loaded.manifest.annotate.argmax.is_some(),
        "{} is not an annotated round",
        state.loaded.file.display()
    );
    let (cells, markers) = state.export(args.out.as_deref())?;
    for f in [cells, markers] {
        log::info!("wrote {}", f.display());
    }
    Ok(())
}

#[cfg(test)]
#[path = "tests/headless.rs"]
pub(crate) mod tests;
