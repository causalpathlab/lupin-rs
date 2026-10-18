//! A trajectory run's order in words: one paragraph on the run (its
//! components, roots, the types where its lineages part and where they
//! end) and one line per type (its lineages, the types either side of it,
//! the verdicts of its edges). The facts are the run's own outputs
//! (`trajectory_lineages.tsv`, `trajectory_edges.parquet` and each type's
//! median pseudotime); with `--lm` the model words them ([`super::lm`]),
//! greedily and deciding nothing, and each text names only the types its
//! facts allow or the template stands in.

use super::consensus::{append, cached_or, chat, check_lm, cited, key, Check};
use super::lm::{Lm, LmSpec};
use crate::manifest::run::{resolve, Loaded};
use crate::trajectory::edges::{self, EdgeRow};
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

/// The run paragraph's and a type line's longest text, in tokens.
const RUN_TOKENS: usize = 120;
const TYPE_TOKENS: usize = 60;

/// One root-to-leaf path.
#[derive(Clone, Debug, PartialEq)]
pub struct Lineage {
    pub name: String,
    pub component: String,
    pub path: Vec<String>,
}

/// What a trajectory run found, as its outputs record it.
#[derive(Clone, Debug, Default)]
pub struct OrderFacts {
    pub lineages: Vec<Lineage>,
    pub edges: Vec<EdgeRow>,
    /// Each type's median pseudotime over the cells a root reaches.
    pub medians: BTreeMap<String, f32>,
}

/// The lineages of `trajectory_lineages.tsv`'s `text`.
#[must_use]
pub fn parse_lineages(text: &str) -> Vec<Lineage> {
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let mut it = l.split('\t');
            let (name, component, path) = (it.next()?, it.next()?, it.next()?);
            Some(Lineage {
                name: name.to_string(),
                component: component.to_string(),
                path: path.split(" > ").map(str::to_string).collect(),
            })
        })
        .collect()
}

/// The median of `v`'s finite values; `None` without any.
fn median(mut v: Vec<f32>) -> Option<f32> {
    v.retain(|x| x.is_finite());
    v.sort_by(f32::total_cmp);
    let n = v.len();
    (n > 0).then(|| {
        if n % 2 == 1 {
            v[n / 2]
        } else {
            0.5 * (v[n / 2 - 1] + v[n / 2])
        }
    })
}

impl OrderFacts {
    /// The facts of the trajectory run `loaded` records.
    pub fn read(loaded: &Loaded) -> Result<Self> {
        let t = &loaded.manifest.trajectory;
        let at = |rel: &Option<String>| rel.as_deref().map(|r| resolve(&loaded.dir, r));
        let edges = match at(&t.edges) {
            Some(p) => edges::read(&p)?,
            None => Vec::new(),
        };
        let lineages = match at(&t.lineages) {
            Some(p) => parse_lineages(
                &std::fs::read_to_string(&p).with_context(|| format!("reading {p}"))?,
            ),
            None => Vec::new(),
        };
        let mut by_type: BTreeMap<String, Vec<f32>> = BTreeMap::new();
        if let Some(p) = at(&t.pseudotime) {
            let (s, n) =
                legume_numeric::matrix::parquet::read_table_columns(&p, &["type"], &["pseudotime"])
                    .with_context(|| format!("reading {p}"))?;
            for (ty, pt) in s[0].iter().zip(&n[0]) {
                by_type.entry(ty.to_string()).or_default().push(*pt as f32);
            }
        }
        let medians = by_type
            .into_iter()
            .filter_map(|(t, v)| Some((t, median(v)?)))
            .collect();
        Ok(Self {
            lineages,
            edges,
            medians,
        })
    }

    /// Where the lineages start, in order.
    #[must_use]
    pub fn roots(&self) -> Vec<String> {
        unique(self.lineages.iter().filter_map(|l| l.path.first().cloned()))
    }

    /// Where they end, in order.
    #[must_use]
    pub fn ends(&self) -> Vec<String> {
        unique(self.lineages.iter().filter_map(|l| l.path.last().cloned()))
    }

    /// The types on more than one lineage of a component where two of
    /// them part: the last type their paths share.
    #[must_use]
    pub fn branch_points(&self) -> Vec<String> {
        let mut out = Vec::new();
        for (i, a) in self.lineages.iter().enumerate() {
            for b in &self.lineages[i + 1..] {
                if a.component != b.component {
                    continue;
                }
                let shared = a
                    .path
                    .iter()
                    .zip(&b.path)
                    .take_while(|(x, y)| x == y)
                    .count();
                if shared > 0 && shared < a.path.len() && shared < b.path.len() {
                    out.push(a.path[shared - 1].clone());
                }
            }
        }
        unique(out.into_iter())
    }

    /// How many components the lineages span.
    #[must_use]
    pub fn components(&self) -> usize {
        self.lineages
            .iter()
            .map(|l| l.component.as_str())
            .collect::<BTreeSet<_>>()
            .len()
    }

    /// Every type on a lineage, then every other one an edge names, in
    /// order.
    #[must_use]
    pub fn types(&self) -> Vec<String> {
        let on = self.lineages.iter().flat_map(|l| l.path.iter().cloned());
        let edged = self
            .edges
            .iter()
            .filter(|e| e.verdict.is_some())
            .flat_map(|e| [e.a.to_string(), e.b.to_string()]);
        unique(on.chain(edged))
    }

    fn lineages_of(&self, t: &str) -> Vec<&Lineage> {
        self.lineages
            .iter()
            .filter(|l| l.path.iter().any(|p| p == t))
            .collect()
    }

    /// The edges of `t` with a verdict.
    fn edges_of(&self, t: &str) -> Vec<&EdgeRow> {
        self.edges
            .iter()
            .filter(|e| e.verdict.is_some() && (e.a.as_ref() == t || e.b.as_ref() == t))
            .collect()
    }

    /// The types either side of `t` on its lineages, before and after.
    fn sides(&self, t: &str) -> (Vec<String>, Vec<String>) {
        let mut before = Vec::new();
        let mut after = Vec::new();
        for l in self.lineages_of(t) {
            let i = l.path.iter().position(|p| p == t).unwrap_or(0);
            before.extend(i.checked_sub(1).map(|j| l.path[j].clone()));
            after.extend(l.path.get(i + 1).cloned());
        }
        (unique(before.into_iter()), unique(after.into_iter()))
    }

    /// What `t`'s line may name besides `t`: its edges' other ends, the
    /// types either side of it and its lineages' roots.
    #[must_use]
    pub fn allowed_for(&self, t: &str) -> BTreeSet<String> {
        let (before, after) = self.sides(t);
        let mut s: BTreeSet<String> = before.into_iter().chain(after).collect();
        for e in self.edges_of(t) {
            s.insert(e.a.to_string());
            s.insert(e.b.to_string());
        }
        s.extend(
            self.lineages_of(t)
                .iter()
                .filter_map(|l| l.path.first().cloned()),
        );
        s.remove(t);
        s
    }
}

/// `items` without repeats, in order.
fn unique(items: impl Iterator<Item = String>) -> Vec<String> {
    let mut seen = BTreeSet::new();
    items.filter(|t| seen.insert(t.clone())).collect()
}

/// `a`, `a and b`, `a, b and c`.
fn and_list(items: &[String]) -> String {
    match items {
        [] => String::new(),
        [a] => a.clone(),
        [rest @ .., last] => format!("{} and {last}", rest.join(", ")),
    }
}

fn edge_text(e: &EdgeRow) -> String {
    let verdict = e.verdict.map_or("", |v| v.as_str());
    format!("{} → {} {verdict}", e.a, e.b)
}

/// The run paragraph when no model words it.
#[must_use]
pub fn run_template(f: &OrderFacts) -> String {
    if f.lineages.is_empty() {
        return "The run ordered no lineage: no root reaches its types.".into();
    }
    let mut s = format!(
        "{} lineage(s) in {} component(s). The order starts at {}.",
        f.lineages.len(),
        f.components(),
        and_list(&f.roots())
    );
    let bps = f.branch_points();
    if !bps.is_empty() {
        s += &format!(" It branches after {}.", and_list(&bps));
    }
    s += &format!(" Its lineages end at {}.", and_list(&f.ends()));
    let paths: Vec<String> = f
        .lineages
        .iter()
        .map(|l| format!("{}: {}", l.name, l.path.join(" > ")))
        .collect();
    format!("{s} {}.", paths.join("; "))
}

/// Type `t`'s line when no model words it.
#[must_use]
pub fn type_template(f: &OrderFacts, t: &str) -> String {
    let on: Vec<String> = f.lineages_of(t).iter().map(|l| l.name.clone()).collect();
    let mut s = if on.is_empty() {
        format!("{t} is on no lineage.")
    } else {
        let (before, after) = f.sides(t);
        let mut place = Vec::new();
        if !before.is_empty() {
            place.push(format!("after {}", and_list(&before)));
        }
        if !after.is_empty() {
            place.push(format!("before {}", and_list(&after)));
        }
        let place = if place.is_empty() {
            String::new()
        } else {
            format!(", {}", place.join(" and "))
        };
        format!("{t} is on {}{place}.", and_list(&on))
    };
    let edges: Vec<String> = f.edges_of(t).into_iter().map(edge_text).collect();
    if !edges.is_empty() {
        s += &format!(" Its edges: {}.", edges.join("; "));
    }
    s
}

const SYSTEM: &str = "You describe the order of cell types along a trajectory.\n\
                      Cite only names that appear below. Plain sentences, no numbers.";

fn facts_text(f: &OrderFacts, types: Option<&BTreeSet<String>>) -> String {
    let keep = |t: &str| types.is_none_or(|s| s.contains(t));
    let mut s = String::new();
    let lineages: Vec<String> = f
        .lineages
        .iter()
        .filter(|l| l.path.iter().any(|p| keep(p)))
        .map(|l| l.path.join(" > "))
        .collect();
    if !lineages.is_empty() {
        s += &format!("LINEAGES:\n{}\n", lineages.join("\n"));
    }
    let edges: Vec<String> = f
        .edges
        .iter()
        .filter(|e| e.verdict.is_some() && keep(&e.a) && keep(&e.b))
        .map(|e| {
            let order = if e.order_agreement.is_finite() {
                format!(", order agreement {:.2}", e.order_agreement)
            } else {
                String::new()
            };
            format!(
                "{} (connectivity {:.2}{order})",
                edge_text(e),
                e.connectivity
            )
        })
        .collect();
    if !edges.is_empty() {
        s += &format!("EDGES: {}\n", edges.join("; "));
    }
    let medians: Vec<String> = f
        .medians
        .iter()
        .filter(|(t, _)| keep(t))
        .map(|(t, m)| format!("{t} {m:.2}"))
        .collect();
    if !medians.is_empty() {
        s += &format!("MEDIAN PSEUDOTIME: {}\n", medians.join("; "));
    }
    s
}

/// One text to write: the run's paragraph (`subject` none) or a type's
/// line, with what it may name and its template.
#[derive(Clone, Debug)]
pub struct Subject {
    pub subject: Option<String>,
    pub prompt: String,
    /// It names one of these.
    pub named: BTreeSet<String>,
    /// It may name these too.
    pub own: BTreeSet<String>,
    pub template: String,
    pub max_tokens: usize,
}

/// The run's paragraph, then each type's line.
#[must_use]
pub fn subjects(f: &OrderFacts) -> Vec<Subject> {
    let (roots, bps, ends) = (f.roots(), f.branch_points(), f.ends());
    let mut user = format!("ROOTS: {}\n", roots.join(", "));
    if !bps.is_empty() {
        user += &format!("BRANCH POINTS: {}\n", bps.join(", "));
    }
    user += &format!("ENDS: {}\n", ends.join(", "));
    user += &facts_text(f, None);
    user += "Describe the components, the roots, where the lineages branch and where they end, \
             in at most three sentences.";
    let mut out = vec![Subject {
        subject: None,
        prompt: chat(SYSTEM, &user, ""),
        named: roots.iter().chain(&bps).chain(&ends).cloned().collect(),
        own: BTreeSet::new(),
        template: run_template(f),
        max_tokens: RUN_TOKENS,
    }];
    for t in f.types() {
        let own = f.allowed_for(&t);
        let mut shown = own.clone();
        shown.insert(t.clone());
        let (before, after) = f.sides(&t);
        let mut user = format!("TYPE: {t}\n");
        if !before.is_empty() {
            user += &format!("BEFORE IT: {}\n", before.join(", "));
        }
        if !after.is_empty() {
            user += &format!("AFTER IT: {}\n", after.join(", "));
        }
        user += &facts_text(f, Some(&shown));
        user += &format!(
            "Describe where {t} sits in the order and what its edges' verdicts say, in one sentence."
        );
        out.push(Subject {
            subject: Some(t.clone()),
            prompt: chat(SYSTEM, &user, ""),
            named: BTreeSet::from([t.clone()]),
            own,
            template: type_template(f, &t),
            max_tokens: TYPE_TOKENS,
        });
    }
    out
}

/// One text, as the cache keeps it.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Described {
    /// The type; none for the run's paragraph.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub subject: Option<String>,
    pub text: String,
    pub cited: Vec<String>,
    pub check: Check,
    pub model: String,
    pub revision: String,
    pub key: String,
}

/// `s`'s text by `lm` (named `model` at `revision`), checked against the
/// run's `types`; the template when it fails.
pub fn word(
    lm: &mut dyn Lm,
    model: &str,
    revision: &str,
    s: &Subject,
    types: &BTreeSet<String>,
) -> Result<Described> {
    let text = super::consensus::checked(
        lm,
        &s.prompt,
        "",
        s.max_tokens,
        super::consensus::numbered(&s.named, &s.own),
        |t| check_lm(t, types, &s.named, &s.own),
    )?;
    let (text, check) = match text {
        Some(t) => (t, Check::Ok),
        None => (s.template.clone(), Check::Template),
    };
    Ok(Described {
        subject: s.subject.clone(),
        cited: cited(&text, &s.named, &s.own),
        text,
        check,
        model: model.to_string(),
        revision: revision.to_string(),
        key: key(&s.prompt, model, revision),
    })
}

/// Each subject's text: from `cache` when it has one on the same facts by
/// the model `spec` names, else worded by the model (loaded by `load`
/// only then) and appended to `cache`. Without a `spec`, the templates.
pub fn describe(
    f: &OrderFacts,
    spec: Option<&LmSpec>,
    cache: &Path,
    load: impl FnOnce() -> Result<Box<dyn Lm>>,
) -> Result<Vec<Described>> {
    let all = subjects(f);
    let types: BTreeSet<String> = f.types().into_iter().collect();
    let Some(spec) = spec else {
        return Ok(all
            .iter()
            .map(|s| Described {
                subject: s.subject.clone(),
                text: s.template.clone(),
                cited: cited(&s.template, &s.named, &s.own),
                check: Check::Template,
                model: String::new(),
                revision: String::new(),
                key: String::new(),
            })
            .collect());
    };
    let (model, revision) = (spec.name(), spec.revision());
    let keys: Vec<String> = all
        .iter()
        .map(|s| key(&s.prompt, &model, &revision))
        .collect();
    cached_or(
        cache,
        keys.into_iter(),
        "texts",
        load,
        |lm, k| {
            let d = word(lm, &model, &revision, &all[k], &types)?;
            append(cache, &d)?;
            Ok(d)
        },
        |d: &Described| d.key.clone(),
    )
}

/// The texts of the run's facts `f` the cache `cache` has, by any model, in
/// subject order: the template where it has none.
#[must_use]
pub fn current(f: &OrderFacts, cache: &Path) -> Vec<Described> {
    let recs = super::consensus::read_cache::<Described>(cache);
    subjects(f)
        .into_iter()
        .map(|s| {
            recs.iter()
                .rev()
                .find(|d| d.subject == s.subject && key(&s.prompt, &d.model, &d.revision) == d.key)
                .cloned()
                .unwrap_or_else(|| Described {
                    subject: s.subject.clone(),
                    cited: cited(&s.template, &s.named, &s.own),
                    text: s.template,
                    check: Check::Template,
                    model: String::new(),
                    revision: String::new(),
                    key: String::new(),
                })
        })
        .collect()
}

/// `{out}.trajectory.describe.md`.
#[must_use]
pub fn md_path(out: &str) -> PathBuf {
    PathBuf::from(format!("{out}.trajectory.describe.md"))
}

/// The texts as Markdown: the paragraph, then a section per type.
#[must_use]
pub fn markdown(from: &str, texts: &[Described]) -> String {
    let mut s = format!("# Trajectory description\n\nFrom `{from}`.\n\n");
    for d in texts {
        if let Some(t) = &d.subject {
            s += &format!("## {t}\n\n");
        }
        s += &format!("{}\n\n", d.text);
    }
    s
}

#[cfg(test)]
#[path = "tests/trajectory.rs"]
pub(crate) mod tests;
