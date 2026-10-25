//! Contrastive questions about a round: why X rather than Y, given Z. A
//! question names a fact (what the round shows or decided), a foil (the
//! rival the asker has in mind; the strongest one when left out) and what
//! is given (a cluster, and for a gene the list it is shown in). The answer
//! is the rule that decides it and the values on both sides, each from a
//! table of the round, and which side they favour.
//!
//! A label's answer triangulates: its lines of evidence are each a contrast,
//! and each raises the gene and term questions that explain it. The panel
//! (the methods' readings, the markers that tell the two types apart and the
//! ontology walk's placement, which scores the types' enrichment) is one
//! line however many readings it holds; the GO/GMT terms are the check from
//! outside it, counted only for terms that still stand out once the panel's
//! genes are taken out of them.

use super::rounds::{caveats, evidence_by_cluster, read_cells, read_table, ClusterEvidence};
use super::run::{resolve, Loaded};
use crate::annotate::celltype_tree::TypeTree;
use crate::annotate::expression::{profile_and_panel, specific_genes, Expression, Panel};
use crate::annotate::go_signature::{read_term_genes, TermGenes};
use crate::annotate::markers::label_key;
use crate::annotate::outputs::{CLUSTER_TERM_P, CLUSTER_TERM_Q_VALUES, CLUSTER_TERM_Z, TERM_GENES};
use crate::annotate::rounds::{ClusterId, Digest, Method, Table, TypeEvidence};
use anyhow::{bail, Context, Result};
use serde::Serialize;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};

/// A term or a set of markers stands out in a cluster below this q (or p).
const ALPHA: f32 = 0.05;
/// Two methods' readings (shares) closer than this are a tie.
const READING_TIE: f32 = 0.01;
/// Genes listed per side as a term's drivers.
const DRIVERS: usize = 5;
/// A term backs a side only from among the cluster's top terms by effect.
const TOP_TERMS: usize = 50;
/// Terms per side checked for support from outside the panel.
const BACKING_TERMS: usize = 3;
/// A term's effect at least this share of the cluster's counts as shared by
/// another cluster.
const SHARED_SHARE: f32 = 0.5;

/// What a question is about.
#[derive(Serialize, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum Kind {
    Label,
    Gene,
    Term,
}

/// The list a gene is shown in, within a cluster.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeneList {
    /// The cluster's specific genes.
    Specific,
    /// A cell type's markers (a coarse label's: its types').
    Markers(String),
    /// A GO/GMT term's genes, by id or name.
    Term(String),
}

impl GeneList {
    /// `specific`, `markers:<type>` or `term:<id or name>`.
    pub fn parse(s: &str) -> Result<Self> {
        match s.split_once(':') {
            None if s == "specific" => Ok(Self::Specific),
            Some(("markers", t)) if !t.is_empty() => Ok(Self::Markers(t.into())),
            Some(("term", t)) if !t.is_empty() => Ok(Self::Term(t.into())),
            _ => bail!(
                "a gene list is `specific`, `markers:<type>` or `term:<id or name>`, not {s:?}"
            ),
        }
    }

    fn name(&self) -> String {
        match self {
            Self::Specific => "specific".into(),
            Self::Markers(t) => format!("markers:{t}"),
            Self::Term(t) => format!("term:{t}"),
        }
    }
}

#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct Question {
    pub kind: Kind,
    pub fact: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub foil: Option<String>,
    pub cluster: ClusterId,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub list: Option<String>,
}

/// Which side an answer, or a line of evidence, favours.
#[derive(Serialize, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Favours {
    Fact,
    Foil,
    /// The evidence does not separate them.
    Neither,
    /// Its readings disagree.
    Split,
    /// Nothing to read it from.
    CannotTell,
}

impl Favours {
    /// What readings that each favour a side favour together: the side all
    /// that lean agree on, `Split` when they lean both ways.
    fn of(readings: impl IntoIterator<Item = Favours>) -> Favours {
        let (mut fact, mut foil, mut any) = (false, false, false);
        for f in readings {
            any |= f != Favours::CannotTell;
            fact |= f == Favours::Fact;
            foil |= f == Favours::Foil;
        }
        match (fact, foil, any) {
            (true, true, _) => Favours::Split,
            (true, false, _) => Favours::Fact,
            (false, true, _) => Favours::Foil,
            (false, false, true) => Favours::Neither,
            (false, false, false) => Favours::CannotTell,
        }
    }

    /// The side with the higher `fact`/`foil`, `Neither` when equal,
    /// `CannotTell` without both.
    fn higher(fact: Option<f32>, foil: Option<f32>) -> Favours {
        match (fact, foil) {
            (Some(a), Some(b)) if a > b => Favours::Fact,
            (Some(a), Some(b)) if a < b => Favours::Foil,
            (Some(_), Some(_)) => Favours::Neither,
            _ => Favours::CannotTell,
        }
    }
}

/// A label's verdict over its independent lines.
#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct Verdict {
    pub favours: Favours,
    /// How many independent lines favour each side.
    pub fact: usize,
    pub foil: usize,
    /// Lines that do not separate them, disagree within, or have nothing.
    pub neither: usize,
}

#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct Answer {
    /// The question; a line of a label's evidence has none of its own.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub question: Option<Question>,
    /// The line of evidence, for a line.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub line: Option<String>,
    /// What decides it, in one line.
    pub rule: String,
    pub fact: Value,
    pub foil: Value,
    pub favours: Favours,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub verdict: Option<Verdict>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub notes: Vec<String>,
    /// The round's tables the values come from.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub sources: Vec<&'static str>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub children: Vec<Answer>,
    /// The round's caveats, at the root.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub caveats: Vec<String>,
}

impl Answer {
    fn new(rule: impl Into<String>, fact: Value, foil: Value, favours: Favours) -> Self {
        Self {
            question: None,
            line: None,
            rule: rule.into(),
            fact,
            foil,
            favours,
            verdict: None,
            notes: Vec::new(),
            sources: Vec::new(),
            children: Vec::new(),
            caveats: Vec::new(),
        }
    }
}

/// The GO/GMT terms a round scored: each term's genes and its cluster ×
/// term effect, z, p and q.
pub struct Terms {
    pub terms: Vec<TermGenes>,
    effect: Table,
    z: Option<Table>,
    p: Option<Table>,
    q: Option<Table>,
    /// Each term id's column in the tables.
    col: BTreeMap<String, usize>,
}

impl Terms {
    fn read(loaded: &Loaded) -> Result<Option<Self>> {
        let Some(rel) = loaded.manifest.annotate.ontology_term_effect.as_deref() else {
            return Ok(None);
        };
        let path = resolve(&loaded.dir, rel);
        let effect = read_table(&path)?;
        let prefix = path
            .strip_suffix(".ontology_term_effect.parquet")
            .unwrap_or(&path);
        let at = |suffix: &str| {
            let p = format!("{prefix}{suffix}");
            std::path::Path::new(&p).is_file().then_some(p)
        };
        let table = |suffix| at(suffix).map(|p| read_table(&p)).transpose();
        let col = (effect.cols.iter().enumerate())
            .map(|(j, id)| (id.clone(), j))
            .collect();
        // A round from before terms kept their genes still has their tables.
        let terms = match at(TERM_GENES) {
            Some(p) => read_term_genes(&p)?,
            None => Vec::new(),
        };
        Ok(Some(Self {
            terms,
            z: table(CLUSTER_TERM_Z)?,
            p: table(CLUSTER_TERM_P)?,
            q: table(CLUSTER_TERM_Q_VALUES)?,
            effect,
            col,
        }))
    }

    fn has_genes(&self) -> bool {
        !self.terms.is_empty()
    }

    /// The term named `name` (an id, or a name in any case).
    fn find(&self, name: &str) -> Result<&TermGenes> {
        let lower = name.to_lowercase();
        if let Some(t) =
            (self.terms.iter()).find(|t| t.id == name || t.name.to_lowercase() == lower)
        {
            return Ok(t);
        }
        let near: Vec<&str> = (self.terms.iter())
            .filter(|t| t.name.to_lowercase().contains(&lower))
            .take(5)
            .map(|t| t.name.as_str())
            .collect();
        if near.is_empty() {
            bail!("no scored term is {name:?}");
        }
        bail!(
            "no scored term is {name:?}; terms that hold it: {}",
            near.join("; ")
        )
    }

    fn get(table: Option<&Table>, k: ClusterId, j: usize) -> Option<f32> {
        let t = table?;
        let i = t.rows.iter().position(|&r| r == k)?;
        Some(t.values[i * t.cols.len() + j]).filter(|v| v.is_finite())
    }

    fn effect(&self, k: ClusterId, id: &str) -> Option<f32> {
        Self::get(Some(&self.effect), k, *self.col.get(id)?)
    }

    fn stat(&self, k: ClusterId, id: &str) -> (Option<f32>, Option<f32>, Option<f32>) {
        let Some(&j) = self.col.get(id) else {
            return (None, None, None);
        };
        let g = |t: &Option<Table>| Self::get(t.as_ref(), k, j);
        (g(&self.z), g(&self.p), g(&self.q))
    }

    /// Cluster `k`'s terms with a positive effect, the largest first.
    fn ranked(&self, k: ClusterId) -> Vec<(&TermGenes, f32)> {
        let mut v: Vec<(&TermGenes, f32)> = (self.terms.iter())
            .filter_map(|t| Some((t, self.effect(k, &t.id)?)))
            .filter(|(_, e)| *e > 0.0)
            .collect();
        v.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.name.cmp(&b.0.name)));
        v
    }

    /// The clusters besides `k` where term `id`'s effect is at least
    /// [`SHARED_SHARE`] of `k`'s.
    fn shared_with(&self, k: ClusterId, id: &str) -> Vec<ClusterId> {
        let Some(e) = self.effect(k, id).filter(|e| *e > 0.0) else {
            return Vec::new();
        };
        (self.effect.rows.iter().copied())
            .filter(|&r| r != k)
            .filter(|&r| self.effect(r, id).is_some_and(|x| x >= SHARED_SHARE * e))
            .collect()
    }
}

/// What a round can answer from.
pub struct Round {
    expr: Option<Expression>,
    panel: Panel,
    tree: Option<TypeTree>,
    specific: BTreeMap<ClusterId, Vec<(String, f32)>>,
    evidence: ClusterEvidence,
    terms: Option<Terms>,
    method: Method,
    caveats: Vec<String>,
}

impl Round {
    pub fn read(loaded: &Loaded) -> Result<Self> {
        let cells = read_cells(loaded)?;
        let (expr, panel) = profile_and_panel(loaded, &cells)?;
        let specific = expr.as_ref().map(specific_genes).unwrap_or_default();
        let (evidence, tree) = evidence_by_cluster(loaded)?;
        Ok(Self {
            expr,
            panel,
            tree,
            specific,
            evidence,
            terms: Terms::read(loaded)?,
            method: super::rounds::called_by(&loaded.manifest.annotate).unwrap_or_default(),
            caveats: caveats(loaded),
        })
    }

    fn cluster(&self, k: ClusterId) -> Result<&(Digest, Vec<TypeEvidence>)> {
        self.evidence
            .get(&k)
            .with_context(|| format!("the round has no cluster K{k}"))
    }

    fn fc(&self, k: ClusterId, gene: &str) -> Option<f32> {
        self.expr.as_ref()?.fold_change(k, gene)
    }

    /// The panel types `label` stands for: its group's members for a coarse
    /// label, else itself.
    fn types_of(&self, label: &str) -> Vec<String> {
        match self
            .tree
            .as_ref()
            .and_then(|t| t.coarse_members(label, &self.panel))
        {
            Some(m) => m.to_vec(),
            None => vec![label.to_string()],
        }
    }

    /// `label`'s markers (a coarse label's: its types').
    fn markers(&self, label: &str) -> BTreeSet<String> {
        (self.types_of(label).iter())
            .filter_map(|t| self.panel.get(&label_key(t)))
            .flatten()
            .cloned()
            .collect()
    }

    /// The types `gene` is a marker of.
    fn marker_of(&self, gene: &str) -> Vec<String> {
        (self.panel.iter())
            .filter(|(_, g)| g.contains(gene))
            .map(|(t, _)| t.clone())
            .collect()
    }

    /// `genes` by their fold change in `k`, the largest first, unmeasured
    /// last.
    fn by_fc<'a>(&self, k: ClusterId, genes: impl IntoIterator<Item = &'a String>) -> Ranked {
        let mut v: Vec<(String, Option<f32>)> = genes
            .into_iter()
            .map(|g| (g.clone(), self.fc(k, g)))
            .collect();
        let f = |x: &Option<f32>| x.unwrap_or(f32::NEG_INFINITY);
        v.sort_by(|a, b| f(&b.1).total_cmp(&f(&a.1)).then_with(|| a.0.cmp(&b.0)));
        v
    }

    fn root(&self, mut a: Answer, q: Question) -> Answer {
        a.question = Some(q);
        a.caveats.clone_from(&self.caveats);
        a
    }
}

// ---- genes ----

/// Genes with their log2 fold change in a cluster, `None` when unmeasured.
type Ranked = Vec<(String, Option<f32>)>;

/// The genes `list` shows in cluster `k`, in its order, with the rule.
fn gene_list(r: &Round, k: ClusterId, list: &GeneList) -> Result<(Ranked, String, &'static str)> {
    Ok(match list {
        GeneList::Specific => (
            (r.specific.get(&k).into_iter().flatten())
                .map(|(g, f)| (g.clone(), Some(*f)))
                .collect(),
            format!("K{k}'s specific genes: those above the cluster's mean expression, by log2 fold change over the other clusters"),
            "cluster_expression",
        ),
        GeneList::Markers(t) => (
            r.by_fc(k, &r.markers(t)),
            format!("{t}'s panel markers, by log2 fold change in K{k}"),
            "markers",
        ),
        GeneList::Term(t) => {
            let terms = r.terms.as_ref().context("the round scored no GO/GMT terms")?;
            anyhow::ensure!(terms.has_genes(), "the round does not keep its terms' genes: annotate it again");
            let term = terms.find(t)?;
            (
                r.by_fc(k, &term.genes),
                format!("{}'s measured genes, by log2 fold change in K{k}", term.name),
                "term_genes",
            )
        }
    })
}

/// Where `gene` stands against `list` in cluster `k`.
fn gene_side(
    r: &Round,
    k: ClusterId,
    gene: &str,
    shown: &[(String, Option<f32>)],
    list: &GeneList,
) -> Value {
    let fc = r.fc(k, gene);
    let mut v = json!({ "gene": gene, "log2fc": fc });
    match shown.iter().position(|(g, _)| g == gene) {
        Some(i) => {
            v["rank"] = json!(i + 1);
            v["of"] = json!(shown.len());
        }
        None => {
            v["shown"] = json!(false);
            v["why_not"] = json!(why_not_shown(r, k, gene, list));
        }
    }
    let types = r.marker_of(gene);
    if !types.is_empty() {
        v["marker_of"] = json!(types);
    }
    if let Some(terms) = r.terms.as_ref() {
        let top: Vec<&str> = (terms.ranked(k).into_iter())
            .filter(|(t, _)| t.genes.iter().any(|g| g == gene))
            .take(3)
            .map(|(t, _)| t.name.as_str())
            .collect();
        if !top.is_empty() {
            v["in_terms"] = json!(top);
        }
    }
    v
}

/// Why `gene` is not in `list` for cluster `k`.
fn why_not_shown(r: &Round, k: ClusterId, gene: &str, list: &GeneList) -> String {
    let Some(e) = r.expr.as_ref() else {
        return "the round has no cluster expression".into();
    };
    let (Some(g), Some(c)) = (e.row(gene), e.col(k)) else {
        return "not measured".into();
    };
    match list {
        GeneList::Specific => {
            let m = &e.table.mat;
            let floor = m.column(c).mean();
            format!(
                "at or below K{k}'s mean expression ({:.3} against a mean of {floor:.3})",
                m[(g, c)]
            )
        }
        GeneList::Markers(t) => format!("not a marker of {t}"),
        GeneList::Term(t) => format!("not in {t}'s genes"),
    }
}

/// Why `fact` rather than `foil` is shown in `list` for cluster `k`; the
/// foil, when left out, is the gene listed just below the fact.
pub fn why_gene(
    r: &Round,
    k: ClusterId,
    fact: &str,
    foil: Option<&str>,
    list: &GeneList,
) -> Result<Answer> {
    r.cluster(k)?;
    let (shown, rule, source) = gene_list(r, k, list)?;
    let at = |g: &str| shown.iter().position(|(s, _)| s == g);
    let mut notes = Vec::new();
    let foil = match foil {
        Some(f) => Some(f.to_string()),
        None => {
            let next = at(fact)
                .and_then(|i| shown.get(i + 1))
                .map(|(g, _)| g.clone());
            if let Some(n) = &next {
                notes.push(format!("no foil given: {n}, the gene listed next"));
            }
            next
        }
    };
    let favours = match (at(fact), foil.as_deref().map(at)) {
        (_, None) => Favours::CannotTell,
        (Some(a), Some(Some(b))) => Favours::higher(Some(-(a as f32)), Some(-(b as f32))),
        (Some(_), Some(None)) => Favours::Fact,
        (None, Some(Some(_))) => Favours::Foil,
        (None, Some(None)) => Favours::Neither,
    };
    let side = |g: &str| gene_side(r, k, g, &shown, list);
    let mut a = Answer::new(
        rule,
        side(fact),
        foil.as_deref().map_or(Value::Null, side),
        favours,
    );
    a.notes = notes;
    a.sources = vec![source, "cluster_expression"];
    Ok(r.root(
        a,
        Question {
            kind: Kind::Gene,
            fact: fact.into(),
            foil,
            cluster: k,
            list: Some(list.name()),
        },
    ))
}

// ---- terms ----

fn term_side(r: &Round, terms: &Terms, k: ClusterId, t: &TermGenes) -> Value {
    let (z, p, q) = terms.stat(k, &t.id);
    let rank = (terms.ranked(k).iter())
        .position(|(x, _)| x.id == t.id)
        .map(|i| i + 1);
    let drivers: Vec<Value> = (r.by_fc(k, &t.genes).into_iter())
        .take(DRIVERS)
        .map(|(g, f)| json!([g, f]))
        .collect();
    let shared = terms.shared_with(k, &t.id);
    let mut v = json!({
        "term": t.name, "id": t.id, "effect": terms.effect(k, &t.id),
        "z": z, "p": p, "q": q, "rank": rank, "genes": t.genes.len(),
    });
    if !drivers.is_empty() {
        v["drivers"] = json!(drivers);
    }
    v["shared_with"] = json!(shared.iter().map(|c| format!("K{c}")).collect::<Vec<_>>());
    v
}

/// Why term `fact` rather than `foil` in cluster `k`; the foil, when left
/// out, is the term ranked just below the fact (or the top term, for a
/// fact not ranked).
pub fn why_term(r: &Round, k: ClusterId, fact: &str, foil: Option<&str>) -> Result<Answer> {
    r.cluster(k)?;
    let terms = r
        .terms
        .as_ref()
        .context("the round scored no GO/GMT terms")?;
    anyhow::ensure!(
        terms.has_genes(),
        "the round does not keep its terms' genes: annotate it again"
    );
    let a = terms.find(fact)?;
    let ranked = terms.ranked(k);
    let mut notes = Vec::new();
    let b = match foil {
        Some(f) => terms.find(f)?,
        None => {
            let i = ranked.iter().position(|(t, _)| t.id == a.id);
            let pick = match i {
                Some(i) => ranked.get(i + 1),
                None => ranked.first(),
            };
            let (b, _) = pick.context("no other term with a positive effect in the cluster")?;
            notes.push(format!("no foil given: {}, the term ranked next", b.name));
            b
        }
    };
    let (ea, eb) = (terms.effect(k, &a.id), terms.effect(k, &b.id));
    let (qa, qb) = (terms.stat(k, &a.id).2, terms.stat(k, &b.id).2);
    let favours = Favours::higher(ea, eb);
    // The order is by effect; the test may read them the other way when one
    // is a programme many clusters share, which widens its null.
    if favours == Favours::Fact && qa.zip(qb).is_some_and(|(x, y)| x > y) {
        notes.push(format!(
            "{} ranks higher by effect, but {} stands out more from matched gene sets (q {:.3} against {:.3}): {} is shared with {} other clusters",
            a.name,
            b.name,
            qb.unwrap_or(f32::NAN),
            qa.unwrap_or(f32::NAN),
            a.name,
            terms.shared_with(k, &a.id).len()
        ));
    }
    let mut ans = Answer::new(
        format!("K{k}'s terms rank by effect (their genes' mean log expression against the other annotated genes, less its median over the clusters); q tests that effect against matched random gene sets, widened for a programme high in many clusters"),
        term_side(r, terms, k, a),
        term_side(r, terms, k, b),
        favours,
    );
    ans.notes = notes;
    ans.sources = vec![
        "ontology_term_effect",
        "cluster_term_z",
        "cluster_term_q_values",
        "term_genes",
    ];
    Ok(r.root(
        ans,
        Question {
            kind: Kind::Term,
            fact: a.name.clone(),
            foil: Some(b.name.clone()),
            cluster: k,
            list: None,
        },
    ))
}

// ---- labels ----

/// Each set's effect in cluster `k` against matched random gene sets, as
/// the GO/GMT terms are tested ([`crate::annotate::term_test`]): its z and p,
/// `None` for a set with no measured gene.
fn set_test(
    expr: &Expression,
    k: ClusterId,
    sets: &[&BTreeSet<String>],
) -> Result<Vec<Option<(f32, f32)>>> {
    let profile = &expr.table.mat;
    let Some(c) = expr.col(k) else {
        return Ok(vec![None; sets.len()]);
    };
    let rows: Vec<Vec<usize>> = (sets.iter())
        .map(|s| s.iter().filter_map(|g| expr.row(g)).collect())
        .collect();
    let named: Vec<(Box<str>, Vec<usize>)> = (rows.iter().enumerate())
        .filter(|(_, r)| !r.is_empty())
        .map(|(i, r)| (i.to_string().into_boxed_str(), r.clone()))
        .collect();
    if named.is_empty() {
        return Ok(vec![None; sets.len()]);
    }
    let universe: Vec<usize> = (0..profile.nrows())
        .filter(|&g| profile.row(g).iter().any(|&v| v > 0.0))
        .collect();
    let ms = enrichment::ontology_module_score(profile, &named, &universe)?;
    let test = crate::annotate::term_test::test_terms(profile, &named, &universe, &ms.effect_kt)?;
    let mut out = vec![None; sets.len()];
    for (j, (i, _)) in named.iter().enumerate() {
        let i: usize = i.parse().expect("an index");
        out[i] = Some((test.z_kt[(c, j)], test.p_kt[(c, j)]));
    }
    Ok(out)
}

/// `label`'s evidence from `types`: its own, or for a coarse label its best
/// type's by `method`; with that type's name.
fn reading<'a>(
    r: &Round,
    types: &'a [TypeEvidence],
    method: Method,
    label: &str,
) -> Option<(&'a TypeEvidence, usize)> {
    let names: BTreeSet<String> = r.types_of(label).iter().map(|t| label_key(t)).collect();
    (method.ranked(types).into_iter().enumerate())
        .find(|(_, e)| names.contains(&label_key(&e.label)))
        .map(|(i, e)| (e, i))
}

/// The cluster's strongest rival to `fact`, by the round's method: its first
/// type that is neither the fact nor one of its types; else, when the
/// method puts forward no other, the other type of the smallest q.
fn rival(r: &Round, types: &[TypeEvidence], fact: &str) -> Option<(String, bool)> {
    let mine: BTreeSet<String> = r.types_of(fact).iter().map(|t| label_key(t)).collect();
    let other = |l: &String| {
        !mine.contains(&label_key(l)) && !r.types_of(l).iter().any(|t| mine.contains(&label_key(t)))
    };
    if let Some(l) = (r.method.ranked(types).into_iter())
        .map(|e| e.label.clone())
        .find(other)
    {
        return Some((l, true));
    }
    let low = |v: Option<f32>| v.unwrap_or(f32::INFINITY);
    (types.iter())
        .filter(|e| other(&e.label))
        .min_by(|a, b| {
            low(a.q)
                .total_cmp(&low(b.q))
                .then(low(b.nes.map(|x| -x)).total_cmp(&low(a.nes.map(|x| -x))))
        })
        .map(|e| (e.label.clone(), false))
}

/// Coarse label `label`'s two best types in the cluster: those the round's
/// method puts forward, best first, then the others in the group's order.
fn inner_pair(r: &Round, types: &[TypeEvidence], label: &str) -> Option<(String, String)> {
    let members = r.types_of(label);
    let mine: BTreeSet<String> = members.iter().map(|t| label_key(t)).collect();
    let ranked = (r.method.ranked(types).into_iter())
        .map(|e| e.label.clone())
        .filter(|l| mine.contains(&label_key(l)));
    let mut seen = BTreeSet::new();
    let mut it = ranked
        .chain(members.iter().cloned())
        .filter(|l| seen.insert(label_key(l)));
    Some((it.next()?, it.next()?))
}

/// The panel line: each method's reading, the distinguishing markers' test
/// and the ontology walk's placement.
fn panel_line(r: &Round, k: ClusterId, fact: &str, foil: &str, depth: usize) -> Result<Answer> {
    let (digest, types) = r.cluster(k)?;
    let mut readings = Vec::new();
    let (mut fa, mut fo) = (serde_json::Map::new(), serde_json::Map::new());
    for m in Method::ALL {
        let (a, b) = (reading(r, types, m, fact), reading(r, types, m, foil));
        if a.is_none() && b.is_none() {
            continue;
        }
        let side = |x: Option<(&TypeEvidence, usize)>| match x {
            Some((e, i)) => {
                json!({ "type": e.label, "rank": i + 1, "reading": m.reading(e), "q": e.q })
            }
            None => json!({ "rank": null }),
        };
        // The higher place in the method's order; one it does not put
        // forward loses to one it does.
        let f = Favours::higher(a.map(|(_, i)| -(i as f32)), b.map(|(_, i)| -(i as f32)));
        let f = match (a, b, f) {
            (Some(_), None, _) => Favours::Fact,
            (None, Some(_), _) => Favours::Foil,
            (Some((x, _)), Some((y, _)), _)
                if (m.reading(x).zip(m.reading(y)))
                    .is_some_and(|(x, y)| (x - y).abs() < READING_TIE) =>
            {
                Favours::Neither
            }
            (_, _, f) => f,
        };
        readings.push((m.name(), f));
        fa.insert(m.name().into(), side(a));
        fo.insert(m.name().into(), side(b));
    }

    // The markers that tell the two apart, tested on the cluster's profile.
    let (ma, mb) = (r.markers(fact), r.markers(foil));
    let only_a: BTreeSet<String> = ma.difference(&mb).cloned().collect();
    let only_b: BTreeSet<String> = mb.difference(&ma).cloned().collect();
    let shared: Vec<&String> = ma.intersection(&mb).collect();
    let mut notes = Vec::new();
    if let Some(expr) = r.expr.as_ref() {
        let tested = set_test(expr, k, &[&only_a, &only_b])?;
        let up = |s: &BTreeSet<String>| {
            r.by_fc(k, s)
                .iter()
                .filter(|(_, f)| f.is_some_and(|f| f >= 1.0))
                .count()
        };
        let side = |s: &BTreeSet<String>, t: Option<(f32, f32)>| {
            let top: Vec<Value> = r
                .by_fc(k, s)
                .into_iter()
                .take(DRIVERS)
                .map(|(g, f)| json!([g, f]))
                .collect();
            json!({ "distinct": s.len(), "up": up(s), "z": t.map(|t| t.0), "p": t.map(|t| t.1), "top": top })
        };
        let (ta, tb) = (tested[0], tested[1]);
        let stands = |t: Option<(f32, f32)>| t.is_some_and(|(_, p)| p < ALPHA);
        let f = match (stands(ta), stands(tb)) {
            (true, false) => Favours::Fact,
            (false, true) => Favours::Foil,
            (true, true) => Favours::higher(ta.map(|t| t.0), tb.map(|t| t.0)),
            (false, false) if ta.is_none() || tb.is_none() => Favours::CannotTell,
            (false, false) => Favours::Neither,
        };
        readings.push(("markers", f));
        fa.insert("markers".into(), side(&only_a, ta));
        fo.insert("markers".into(), side(&only_b, tb));
        if only_a.is_empty() || only_b.is_empty() {
            notes.push(
                "one side has no marker the other lacks: its markers cannot tell them apart".into(),
            );
        }
        if !shared.is_empty() {
            notes.push(format!(
                "{} shared markers say nothing of the contrast: {}",
                shared.len(),
                shared
                    .iter()
                    .map(|s| s.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
        }
    }

    // Where the ontology walk put the cluster: on either type's term (or
    // its group's), or elsewhere.
    if let (Some(cl), Some(tree)) = (digest.cl.as_ref(), r.tree.as_ref()) {
        let on = |l: &str| {
            let mut ids: Vec<String> = r
                .types_of(l)
                .iter()
                .filter_map(|t| tree.placement(t))
                .map(|p| p.id)
                .collect();
            ids.extend(tree.placement(l).map(|p| p.id));
            ids.contains(&cl.id)
        };
        let f = match (on(fact), on(foil)) {
            (true, false) => Favours::Fact,
            (false, true) => Favours::Foil,
            _ => Favours::Neither,
        };
        readings.push(("placement", f));
        fa.insert("placement".into(), json!(on(fact)));
        fo.insert("placement".into(), json!(on(foil)));
        notes.push(format!(
            "the ontology walk placed K{k} at {} {}",
            cl.id, cl.name
        ));
    }

    let favours = Favours::of(readings.iter().map(|(_, f)| *f));
    let mut line = Answer::new(
        "the panel: each method's order, a test of the markers only one type has against matched random gene sets, and the ontology walk's placement; one line, as all rest on the marker panel",
        Value::Object(fa),
        Value::Object(fo),
        favours,
    );
    line.line = Some("panel".into());
    line.notes = notes;
    line.notes.push(format!(
        "readings: {}",
        readings
            .iter()
            .map(|(n, f)| format!("{n} {}", favours_word(*f)))
            .collect::<Vec<_>>()
            .join(", ")
    ));
    line.sources = vec![
        "cluster_celltype_q",
        "cluster_expression",
        "markers",
        "ontology_assignment",
    ];
    if depth > 1 {
        let top = |s: &BTreeSet<String>| r.by_fc(k, s).into_iter().next().map(|(g, _)| g);
        if let (Some(a), Some(b)) = (top(&only_a), top(&only_b)) {
            let mut child = why_gene(r, k, &a, Some(&b), &GeneList::Specific)?;
            child.caveats.clear();
            line.children.push(child);
        }
    }
    Ok(line)
}

/// The terms backing one side: among the cluster's top terms by effect,
/// significant, and holding more of this side's distinct markers that are
/// up in the cluster (log2 FC ≥ 1) than of the other's; each tested again
/// without any panel gene.
fn backing(
    r: &Round,
    terms: &Terms,
    k: ClusterId,
    (distinct, other): (&BTreeSet<String>, &BTreeSet<String>),
    panel_genes: &BTreeSet<String>,
) -> Result<Vec<Value>> {
    let Some(expr) = r.expr.as_ref() else {
        return Ok(Vec::new());
    };
    let up = |t: &TermGenes, s: &BTreeSet<String>| {
        (t.genes.iter())
            .filter(|g| s.contains(*g) && r.fc(k, g).is_some_and(|f| f >= 1.0))
            .count()
    };
    let picked: Vec<&TermGenes> = (terms.ranked(k).into_iter())
        .take(TOP_TERMS)
        .filter(|(t, _)| terms.stat(k, &t.id).2.is_some_and(|q| q < ALPHA))
        .filter(|(t, _)| up(t, distinct) > up(t, other))
        .map(|(t, _)| t)
        .take(BACKING_TERMS)
        .collect();
    let rest: Vec<BTreeSet<String>> = (picked.iter())
        .map(|t| {
            t.genes
                .iter()
                .filter(|g| !panel_genes.contains(*g))
                .cloned()
                .collect()
        })
        .collect();
    let tested = set_test(expr, k, &rest.iter().collect::<Vec<_>>())?;
    Ok((picked.iter().zip(&rest).zip(tested))
        .map(|((t, rest), without)| {
            let overlap: Vec<&String> = t.genes.iter().filter(|g| distinct.contains(*g)).collect();
            json!({
                "term": t.name, "id": t.id, "q": terms.stat(k, &t.id).2, "effect": terms.effect(k, &t.id),
                "markers_in_it": overlap, "genes": t.genes.len(), "genes_without_panel": rest.len(),
                "without_panel": without.map(|(z, p)| json!({ "z": z, "p": p })),
                "independent": without.is_some_and(|(_, p)| p < ALPHA),
            })
        })
        .collect())
}

/// The term holding some of `distinct` with the largest effect in `k`.
fn holding<'a>(
    terms: &'a Terms,
    k: ClusterId,
    distinct: &BTreeSet<String>,
) -> Option<&'a TermGenes> {
    (terms.terms.iter())
        .filter(|t| t.genes.iter().any(|g| distinct.contains(g)))
        .filter_map(|t| Some((t, terms.effect(k, &t.id)?)))
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(t, _)| t)
}

/// The GO/GMT line: the terms that back each side, counted only where they
/// still stand out without the panel's genes.
fn terms_line(r: &Round, k: ClusterId, fact: &str, foil: &str, depth: usize) -> Result<Answer> {
    let rule = "GO/GMT terms among the cluster's top 50 by effect, significant, and holding more of one type's up distinct markers than of the other's; counted where the term still stands out once every panel gene is taken out of it";
    let Some(terms) = r.terms.as_ref().filter(|t| t.has_genes()) else {
        let mut a = Answer::new(rule, Value::Null, Value::Null, Favours::CannotTell);
        a.line = Some("terms".into());
        a.notes.push(match r.terms {
            None => "the round scored no GO/GMT terms (annotate with --go, --gaf or --gmt)".into(),
            Some(_) => "the round does not keep its terms' genes: annotate it again".into(),
        });
        return Ok(a);
    };
    let (ma, mb) = (r.markers(fact), r.markers(foil));
    let panel_genes: BTreeSet<String> = r.panel.values().flatten().cloned().collect();
    let only_a: BTreeSet<String> = ma.difference(&mb).cloned().collect();
    let only_b: BTreeSet<String> = mb.difference(&ma).cloned().collect();
    let ba = backing(r, terms, k, (&only_a, &only_b), &panel_genes)?;
    let bb = backing(r, terms, k, (&only_b, &only_a), &panel_genes)?;
    let backed = |v: &[Value]| v.iter().any(|t| t["independent"] == json!(true));
    let favours = match (backed(&ba), backed(&bb)) {
        (true, false) => Favours::Fact,
        (false, true) => Favours::Foil,
        _ => Favours::Neither,
    };
    let mut a = Answer::new(
        rule,
        json!({ "terms": ba }),
        json!({ "terms": bb }),
        favours,
    );
    a.line = Some("terms".into());
    a.sources = vec![
        "ontology_term_effect",
        "cluster_term_q_values",
        "term_genes",
        "cluster_expression",
    ];
    if ba.is_empty() && bb.is_empty() {
        a.notes
            .push("no top significant term holds either type's up distinct markers".into());
    }
    if depth > 1 {
        // Each side's best term: a backing one, else the term holding its
        // distinct markers with the largest effect in the cluster.
        let best = |backing: &[Value], distinct: &BTreeSet<String>| {
            (backing
                .first()
                .and_then(|t| t["id"].as_str())
                .map(str::to_string))
            .or_else(|| holding(terms, k, distinct).map(|t| t.id.clone()))
        };
        if let (Some(x), Some(y)) = (best(&ba, &only_a), best(&bb, &only_b)) {
            if x != y {
                let mut child = why_term(r, k, &x, Some(&y))?;
                child.caveats.clear();
                a.children.push(child);
            }
        }
    }
    Ok(a)
}

/// A label, with the panel types it stands for when it is coarse.
fn label_side(r: &Round, label: &str) -> Value {
    let types = r.types_of(label);
    match types.as_slice() {
        [t] if label_key(t) == label_key(label) => json!({ "label": label }),
        _ => json!({ "label": label, "types": types }),
    }
}

/// Why cluster `k` is `fact` rather than `foil`: the panel line and the
/// GO/GMT line, each with its gene or term questions to `depth` levels,
/// and the verdict over them. The fact defaults to the cluster's label; the
/// foil to its strongest rival by the round's method.
pub fn why_label(
    r: &Round,
    k: ClusterId,
    fact: Option<&str>,
    foil: Option<&str>,
    depth: usize,
) -> Result<Answer> {
    let (digest, types) = r.cluster(k)?;
    let fact = match fact {
        Some(f) => f.to_string(),
        None => digest
            .label
            .clone()
            .with_context(|| format!("K{k} is unassigned: name the label to ask about"))?,
    };
    let mut notes = Vec::new();
    let (fact, foil) = match foil {
        Some(f) => (fact, f.to_string()),
        None => {
            let (a, b, note) = match (rival(r, types, &fact), inner_pair(r, types, &fact)) {
                (Some((f, true)), _) => (fact, f.clone(), format!(
                    "no foil given: {f}, the strongest rival by the {} order",
                    r.method.name()
                )),
                // A coarse label with no rival outside it: the question it
                // leaves open is between its own types.
                (_, Some((a, b))) => (a.clone(), b.clone(), format!(
                    "no foil given and no rival outside {fact}: asking between its types {a} and {b}"
                )),
                (Some((f, false)), None) => (fact, f.clone(), format!(
                    "no foil given, and the {} order puts forward no other type: {f}, the closest by q",
                    r.method.name()
                )),
                (None, None) => bail!("no rival to {fact} in K{k}: name one"),
            };
            notes.push(note);
            (a, b)
        }
    };
    let lines = vec![
        panel_line(r, k, &fact, &foil, depth)?,
        terms_line(r, k, &fact, &foil, depth)?,
    ];
    let count = |f: Favours| lines.iter().filter(|l| l.favours == f).count();
    let (n_fact, n_foil) = (count(Favours::Fact), count(Favours::Foil));
    let favours = match (n_fact, n_foil) {
        (0, 0) => Favours::Neither,
        (_, 0) => Favours::Fact,
        (0, _) => Favours::Foil,
        _ => Favours::Split,
    };
    if favours == Favours::Neither {
        notes.push(format!("the evidence does not separate {fact} and {foil}: a label over both (their group) says no more than it holds"));
    }
    let mut a = Answer::new(
        "the lines of evidence that rest on different data: the marker panel, and the GO/GMT terms without the panel's genes",
        label_side(r, &fact),
        label_side(r, &foil),
        favours,
    );
    a.verdict = Some(Verdict {
        favours,
        fact: n_fact,
        foil: n_foil,
        neither: lines.len() - n_fact - n_foil,
    });
    a.notes = notes;
    a.children = lines;
    Ok(r.root(
        a,
        Question {
            kind: Kind::Label,
            fact,
            foil: Some(foil),
            cluster: k,
            list: None,
        },
    ))
}

// ---- the command ----

/// How to ask, said the same way to agents (the MCP instructions), to users
/// (the README, the TUI help) and here.
pub const ASKING: &str = "Ask why in contrasts: why X rather than Y, given Z. The answer weighs \
lines of evidence that rest on different data, says which side each favours, \
and holds the gene and term questions behind them. A tie or a disagreement is \
an answer: the evidence does not separate them.";

/// What a decision quotes in its `evidence` to keep the answer it rests on:
/// the question and the verdict (or what the answer favours).
#[must_use]
pub fn citation(a: &Answer) -> Value {
    json!({
        "why": a.question,
        "favours": a.favours,
        "verdict": a.verdict,
    })
}

/// The three slots of a question, as `lupin why` and the MCP tool take them.
pub const WHY_HELP: &str = "\
Ask why the round shows X rather than Y, given Z: a contrast.
  fact   what the round shows or decided (a label, a gene, a term)
  foil   the rival you have in mind; left out, the strongest one
  given  the cluster, and for a gene the list it is shown in
The answer is the rule that decides it, the values on both sides (each from
a table of the round) and which side they favour. A label's answer weighs its
independent lines of evidence (the marker panel; the GO/GMT terms without the
panel's genes), each with the gene and term questions behind it. A tie or a
disagreement is an answer: the evidence does not separate them.

Examples:
  lupin why -f run.r1 -c 0                       K0's label against its strongest rival
  lupin why -f run.r1 -c 0 --label NK --over gdT
  lupin why -f run.r1 -c 0 --gene NCAM1 --over CD3E --list markers:NK
  lupin why -f run.r1 -c 5 --term 'inflammatory response'";

#[derive(clap::Args, Debug)]
#[command(after_long_help = WHY_HELP)]
pub struct WhyArgs {
    #[arg(
        long,
        short = 'f',
        help = "An annotated manifest (a round) or its prefix"
    )]
    pub from: Box<str>,

    #[arg(long, short = 'c', help = "The cluster the question is about")]
    pub cluster: ClusterId,

    #[arg(
        long,
        group = "fact",
        help = "Why this label (the cluster's own when no fact is named)"
    )]
    pub label: Option<String>,

    #[arg(long, group = "fact", help = "Why this gene, in the list --list names")]
    pub gene: Option<String>,

    #[arg(long, group = "fact", help = "Why this GO/GMT term (an id or a name)")]
    pub term: Option<String>,

    #[arg(long, help = "Rather than this (the strongest rival when left out)")]
    pub over: Option<String>,

    #[arg(
        long,
        default_value = "specific",
        help = "For --gene: the list it is shown in: specific, markers:<type> or term:<id or name>"
    )]
    pub list: String,

    #[arg(
        long,
        default_value_t = 2,
        help = "For a label: levels of gene and term questions under its lines (1: the lines only)"
    )]
    pub depth: usize,

    #[arg(long, help = "Print the answer as JSON")]
    pub json: bool,
}

/// Answer the question `args` puts to its round.
pub fn ask(args: &WhyArgs) -> Result<Answer> {
    let loaded = super::run::load(&args.from)?;
    ask_round(&Round::read(&loaded)?, args)
}

/// Answer the question `args` puts to round `r`, already read.
pub fn ask_round(r: &Round, args: &WhyArgs) -> Result<Answer> {
    let (k, over) = (args.cluster, args.over.as_deref());
    match (&args.gene, &args.term) {
        (Some(g), _) => why_gene(r, k, g, over, &GeneList::parse(&args.list)?),
        (_, Some(t)) => why_term(r, k, t, over),
        _ => why_label(r, k, args.label.as_deref(), over, args.depth),
    }
}

pub fn run_why(args: &WhyArgs) -> Result<()> {
    let a = ask(args)?;
    if args.json {
        println!("{}", serde_json::to_string_pretty(&a)?);
    } else {
        print!("{}", render(&a));
    }
    Ok(())
}

// ---- text ----

fn favours_word(f: Favours) -> &'static str {
    match f {
        Favours::Fact => "favours the fact",
        Favours::Foil => "favours the foil",
        Favours::Neither => "does not separate them",
        Favours::Split => "split",
        Favours::CannotTell => "cannot tell",
    }
}

/// `v` on one line: an object's fields as `key=value` (a nested one in
/// braces), numbers to three significant places.
fn compact(v: &Value) -> String {
    compact_at(v, 0)
}

fn compact_at(v: &Value, depth: usize) -> String {
    match v {
        Value::Null => "-".into(),
        Value::Number(n) => n.as_f64().map_or_else(
            || n.to_string(),
            |x| {
                if x.fract() == 0.0 && x.abs() < 1e6 {
                    format!("{x}")
                } else if x.abs() >= 0.01 {
                    format!("{x:.3}")
                } else {
                    format!("{x:.1e}")
                }
            },
        ),
        Value::String(s) => s.clone(),
        Value::Array(a) => {
            // Items that hold spaces themselves (objects, term names) are
            // set apart.
            let spaced = |x: &Value| x.is_object() || x.as_str().is_some_and(|t| t.contains(' '));
            let sep = if a.iter().any(spaced) { "; " } else { " " };
            let items: Vec<String> = a.iter().map(|x| compact_at(x, depth + 1)).collect();
            format!("[{}]", items.join(sep))
        }
        Value::Object(o) => {
            let fields: Vec<String> = (o.iter())
                .filter(|(_, v)| !v.is_null())
                .map(|(k, v)| format!("{k}={}", compact_at(v, depth + 1)))
                .collect();
            match depth {
                0 => fields.join(", "),
                _ => format!("{{{}}}", fields.join(", ")),
            }
        }
        Value::Bool(b) => b.to_string(),
    }
}

/// The answer as indented text: the question, the rule, both sides and
/// what they favour, then its children.
#[must_use]
pub fn render(a: &Answer) -> String {
    let mut s = String::new();
    render_into(a, 0, ("fact", "foil"), &mut s);
    s
}

/// `a` at `depth`; a line without a question of its own names its sides as
/// its parent's do (`names`).
fn render_into(a: &Answer, depth: usize, names: (&str, &str), s: &mut String) {
    use std::fmt::Write;
    let pad = "  ".repeat(depth);
    let (fact, foil) = match &a.question {
        Some(q) => {
            let foil = q.foil.as_deref().unwrap_or("-");
            let given = q.list.as_ref().map_or_else(
                || format!("K{}", q.cluster),
                |l| format!("K{}, {l}", q.cluster),
            );
            let _ = writeln!(s, "{pad}why {} rather than {foil}, given {given}?", q.fact);
            (q.fact.clone(), foil.to_string())
        }
        None => {
            let _ = writeln!(s, "{pad}{} line", a.line.as_deref().unwrap_or("a"));
            (names.0.to_string(), names.1.to_string())
        }
    };
    let _ = writeln!(s, "{pad}  rule: {}", a.rule);
    let _ = writeln!(s, "{pad}  {fact}: {}", compact(&a.fact));
    let _ = writeln!(s, "{pad}  {foil}: {}", compact(&a.foil));
    match &a.verdict {
        Some(v) => {
            let _ = writeln!(
                s,
                "{pad}  verdict: {} ({} lines for {fact}, {} for {foil}, {} neither)",
                favours_word(v.favours),
                v.fact,
                v.foil,
                v.neither
            );
        }
        None => {
            let _ = writeln!(s, "{pad}  {}", favours_word(a.favours));
        }
    }
    for n in &a.notes {
        let _ = writeln!(s, "{pad}  note: {n}");
    }
    for c in &a.caveats {
        let _ = writeln!(s, "{pad}  caveat: {c}");
    }
    for c in &a.children {
        render_into(c, depth + 1, (&fact, &foil), s);
    }
}

#[cfg(test)]
#[path = "tests/why.rs"]
mod tests;
