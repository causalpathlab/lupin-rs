//! The prior: which cell types precede which (`docs/trajectory-plan.md` §3).
//!
//! Statements come from the Cell Ontology's `develops_from` relations and from
//! `precedence.tsv` files along lupin's data-file search path, a later layer's
//! statement about a pair of types replacing earlier ones. The `precedes`
//! statements are closed transitively over every type, restricted to the node
//! types and reduced to the direct edges, which must form a DAG.

use crate::annotate::celltype_tree::ClTerms;
use anyhow::{bail, ensure, Result};
use log::info;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// What a statement says about two types.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Relation {
    Precedes,
    Unrelated,
}

impl Relation {
    fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "precedes" => Some(Self::Precedes),
            "unrelated" => Some(Self::Unrelated),
            _ => None,
        }
    }

    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Precedes => "precedes",
            Self::Unrelated => "unrelated",
        }
    }
}

/// Where a statement came from, least to most specific.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Source {
    /// The Cell Ontology, by the type's own `develops_from`.
    Cl,
    /// The Cell Ontology, through an `is_a` ancestor of the type.
    ClInherited,
    /// The user's `precedence.tsv`.
    User,
    /// The project's `precedence.tsv`, beside the run.
    Project,
    /// The `--prior` file.
    Run,
    /// `--root`.
    Cli,
}

impl Source {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Cl => "cl",
            Self::ClInherited => "cl-inherited",
            Self::User => "user",
            Self::Project => "project",
            Self::Run => "run",
            Self::Cli => "cli",
        }
    }
}

/// `from` precedes `to`, or the two are unrelated.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Statement {
    pub(crate) from: String,
    pub(crate) to: String,
    pub(crate) relation: Relation,
    pub(crate) source: Source,
    pub(crate) note: String,
}

impl Statement {
    fn pair(&self) -> (String, String) {
        let (a, b) = (&self.from, &self.to);
        if a <= b {
            (a.clone(), b.clone())
        } else {
            (b.clone(), a.clone())
        }
    }
}

/// What the Cell Ontology says about the labels in `label_cl` (label → term):
/// B develops from A when B's term, or an `is_a` ancestor of it, develops from
/// A's term or from a term that is an `is_a` descendant of it, directly or
/// through intermediate terms. A statement found only through an ancestor of
/// B is tagged [`Source::ClInherited`].
pub(crate) fn from_ontology(
    terms: &ClTerms,
    label_cl: &BTreeMap<String, String>,
) -> Vec<Statement> {
    let mut labels_of: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for (label, id) in label_cl {
        labels_of.entry(id.as_str()).or_default().push(label);
    }
    let mut out = Vec::new();
    for (label_b, b) in label_cl {
        // Terms reached from `b` by develops_from steps, with is_a steps
        // allowed anywhere; `inherited` records an is_a step before the first
        // develops_from step. Keep the better (not inherited) flag per term.
        let mut reached: BTreeMap<String, bool> = BTreeMap::new();
        let mut seen: BTreeSet<(String, bool, bool)> = BTreeSet::new();
        let mut queue = VecDeque::from([(b.clone(), false, false)]);
        while let Some((x, inherited, developed)) = queue.pop_front() {
            if !seen.insert((x.clone(), inherited, developed)) {
                continue;
            }
            for d in terms.develops_from(&x) {
                let e = reached.entry(d.clone()).or_insert(inherited);
                *e = *e && inherited;
                queue.push_back((d.clone(), inherited, true));
            }
            for p in terms.parents(&x) {
                if developed {
                    // An ancestor of a reached term is reached too.
                    let e = reached.entry(p.clone()).or_insert(inherited);
                    *e = *e && inherited;
                }
                queue.push_back((p.clone(), inherited || !developed, developed));
            }
        }
        for (a, &inherited) in &reached {
            for label_a in labels_of.get(a.as_str()).into_iter().flatten() {
                if *label_a == label_b {
                    continue;
                }
                out.push(Statement {
                    from: (*label_a).to_string(),
                    to: label_b.clone(),
                    relation: Relation::Precedes,
                    source: if inherited {
                        Source::ClInherited
                    } else {
                        Source::Cl
                    },
                    note: format!("{} develops from {}", b, a),
                });
            }
        }
    }
    out
}

/// Statements from a `precedence.tsv` text: `from<TAB>to<TAB>relation[<TAB>note]`,
/// `#` comments and a `from…` header skipped. `origin` names the file in errors.
pub(crate) fn parse_statements(text: &str, source: Source, origin: &str) -> Result<Vec<Statement>> {
    let mut out = Vec::new();
    for (n, line) in text.lines().enumerate() {
        let line = line.trim_end();
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').map(str::trim).collect();
        if n == 0 && fields.first() == Some(&"from") {
            continue;
        }
        ensure!(
            fields.len() >= 3,
            "{origin}:{}: expected from<TAB>to<TAB>relation, got {line:?}",
            n + 1
        );
        let relation = Relation::parse(fields[2]).ok_or_else(|| {
            anyhow::anyhow!(
                "{origin}:{}: relation must be `precedes` or `unrelated`, got {:?}",
                n + 1,
                fields[2]
            )
        })?;
        ensure!(
            fields[0] != fields[1],
            "{origin}:{}: a type cannot precede itself ({})",
            n + 1,
            fields[0]
        );
        out.push(Statement {
            from: fields[0].to_string(),
            to: fields[1].to_string(),
            relation,
            source,
            note: fields[3..].join(" "),
        });
    }
    Ok(out)
}

/// One statement per pair of types: within a layer the last wins, and a later
/// layer's statement about a pair replaces an earlier layer's whatever its
/// direction. Sorted by the pair.
pub(crate) fn combine(layers: &[Vec<Statement>]) -> Vec<Statement> {
    let mut by_pair: BTreeMap<(String, String), Statement> = BTreeMap::new();
    for layer in layers {
        for s in layer {
            by_pair.insert(s.pair(), s.clone());
        }
    }
    by_pair.into_values().collect()
}

/// A direct edge of the reduced prior, over node indices.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Edge {
    pub(crate) from: usize,
    pub(crate) to: usize,
    /// The statement's source when the edge is stated outright.
    pub(crate) source: Option<Source>,
    /// Types the edge passes through when it is implied by a path (through a
    /// dropped type, or a shortcut the reduction removed).
    pub(crate) via: Vec<usize>,
}

/// The combined statements over a run's types and the direct edges they imply.
#[derive(Debug)]
pub(crate) struct Prior {
    pub(crate) statements: Vec<Statement>,
    pub(crate) edges: Vec<Edge>,
    /// Each node's prior component, `None` for a node with no edge.
    pub(crate) component: Vec<Option<usize>>,
    /// Nodes with an edge but no incoming one, per component.
    pub(crate) roots: Vec<Vec<usize>>,
}

/// Build the prior over `types` (every label in the run), of which `is_node`
/// marks the node types. `forced_roots` (`--root`) must be nodes with no
/// incoming statement; a forced root with no path to a node gets a `cli`
/// statement to it, so a run with no other prior still orders from it.
pub(crate) fn build(
    types: &[Box<str>],
    is_node: &[bool],
    statements: Vec<Statement>,
    forced_roots: &[&str],
) -> Result<Prior> {
    let index = |name: &str| types.iter().position(|t| t.as_ref() == name);
    let n = types.len();
    let mut statements: Vec<Statement> = statements
        .into_iter()
        .filter(|s| {
            let known = index(&s.from).is_some() && index(&s.to).is_some();
            if !known {
                info!(
                    "prior statement about a type not in this run is set aside: {} {} {}",
                    s.from,
                    s.relation.as_str(),
                    s.to
                );
            }
            known
        })
        .collect();

    let adjacency = |statements: &[Statement]| -> Vec<Vec<usize>> {
        let mut adj = vec![Vec::new(); n];
        for s in statements
            .iter()
            .filter(|s| s.relation == Relation::Precedes)
        {
            adj[index(&s.from).unwrap()].push(index(&s.to).unwrap());
        }
        adj
    };
    let mut adj = adjacency(&statements);
    let mut reach = reachability(&adj);

    for &r in forced_roots {
        let ri = index(r).ok_or_else(|| {
            anyhow::anyhow!(
                "--root {r:?} is not a type of this run; the types are: {}",
                types.join(", ")
            )
        })?;
        ensure!(is_node[ri], "--root {r:?} has too few cells to be a node");
        if let Some(s) = statements
            .iter()
            .find(|s| s.relation == Relation::Precedes && s.to == r)
        {
            bail!(
                "--root {r:?} cannot be a root: {} precedes it ({})",
                s.from,
                s.source.as_str()
            );
        }
        for t in (0..n).filter(|&t| t != ri && is_node[t] && !reach[ri].contains(&t)) {
            statements.push(Statement {
                from: r.to_string(),
                to: types[t].to_string(),
                relation: Relation::Precedes,
                source: Source::Cli,
                note: "--root".into(),
            });
        }
        adj = adjacency(&statements);
        reach = reachability(&adj);
    }

    for a in 0..n {
        for &b in &reach[a] {
            if b != a && reach[b].contains(&a) {
                let cycle = path(&adj, a, b)
                    .into_iter()
                    .chain(path(&adj, b, a).into_iter().skip(1))
                    .map(|i| types[i].to_string())
                    .collect::<Vec<_>>();
                let involved: Vec<String> = statements
                    .iter()
                    .filter(|s| {
                        s.relation == Relation::Precedes
                            && cycle.contains(&s.from)
                            && cycle.contains(&s.to)
                    })
                    .map(|s| format!("{} precedes {} ({})", s.from, s.to, s.source.as_str()))
                    .collect();
                bail!(
                    "the prior has a cycle: {}; from {}",
                    cycle.join(" → "),
                    involved.join("; ")
                );
            }
        }
    }

    let nodes: Vec<usize> = (0..n).filter(|&i| is_node[i]).collect();
    let mut edges = Vec::new();
    for &a in &nodes {
        for &b in &nodes {
            if a == b || !reach[a].contains(&b) {
                continue;
            }
            let implied = nodes
                .iter()
                .any(|&c| c != a && c != b && reach[a].contains(&c) && reach[c].contains(&b));
            if implied {
                continue;
            }
            let stated = statements
                .iter()
                .find(|s| {
                    s.relation == Relation::Precedes
                        && index(&s.from) == Some(a)
                        && index(&s.to) == Some(b)
                })
                .map(|s| s.source);
            let via = if stated.is_some() {
                Vec::new()
            } else {
                let p = path(&adj, a, b);
                p[1..p.len() - 1].to_vec()
            };
            edges.push(Edge {
                from: a,
                to: b,
                source: stated,
                via,
            });
        }
    }

    // Components over the direct edges, numbered in order of their lowest node.
    let mut component = vec![None; n];
    let mut next = 0;
    for &start in &nodes {
        if component[start].is_some() || !edges.iter().any(|e| e.from == start || e.to == start) {
            continue;
        }
        let mut queue = VecDeque::from([start]);
        while let Some(x) = queue.pop_front() {
            if component[x].is_some() {
                continue;
            }
            component[x] = Some(next);
            for e in &edges {
                if e.from == x {
                    queue.push_back(e.to);
                } else if e.to == x {
                    queue.push_back(e.from);
                }
            }
        }
        next += 1;
    }
    let mut roots = vec![Vec::new(); next];
    for &i in &nodes {
        if let Some(c) = component[i] {
            if !edges.iter().any(|e| e.to == i) {
                roots[c].push(i);
            }
        }
    }
    Ok(Prior {
        statements,
        edges,
        component,
        roots,
    })
}

/// Every node's set of nodes reachable by one or more steps.
fn reachability(adj: &[Vec<usize>]) -> Vec<BTreeSet<usize>> {
    (0..adj.len())
        .map(|start| {
            let mut seen = BTreeSet::new();
            let mut queue: VecDeque<usize> = adj[start].iter().copied().collect();
            while let Some(x) = queue.pop_front() {
                if seen.insert(x) {
                    queue.extend(adj[x].iter().copied());
                }
            }
            seen
        })
        .collect()
}

/// A shortest path from `a` to `b` (both ends included); empty if none.
fn path(adj: &[Vec<usize>], a: usize, b: usize) -> Vec<usize> {
    let mut parent = vec![None; adj.len()];
    let mut queue = VecDeque::from([a]);
    let mut seen = BTreeSet::from([a]);
    while let Some(x) = queue.pop_front() {
        for &y in &adj[x] {
            if seen.insert(y) {
                parent[y] = Some(x);
                if y == b {
                    let mut p = vec![b];
                    while let Some(q) = parent[*p.last().unwrap()] {
                        p.push(q);
                    }
                    p.reverse();
                    return p;
                }
                queue.push_back(y);
            }
        }
    }
    Vec::new()
}

/// `{out}.trajectory_prior.tsv`: the combined statements and the direct edges,
/// each with its source, a complete prior on its own (`--prior … --prior-only`
/// reads it back).
pub(crate) fn write_tsv(path: &str, types: &[Box<str>], prior: &Prior) -> Result<()> {
    let mut text = String::from("# kind\tfrom\tto\trelation\tsource\tnote\n");
    for s in &prior.statements {
        text += &format!(
            "statement\t{}\t{}\t{}\t{}\t{}\n",
            s.from,
            s.to,
            s.relation.as_str(),
            s.source.as_str(),
            s.note
        );
    }
    for e in &prior.edges {
        let (source, note) = match e.source {
            Some(s) => (s.as_str().to_string(), String::new()),
            None => (
                "implied".to_string(),
                format!(
                    "via {}",
                    e.via
                        .iter()
                        .map(|&i| types[i].as_ref())
                        .collect::<Vec<_>>()
                        .join(", ")
                ),
            ),
        };
        text += &format!(
            "edge\t{}\t{}\tprecedes\t{source}\t{note}\n",
            types[e.from], types[e.to]
        );
    }
    std::fs::write(path, text).map_err(|e| anyhow::anyhow!("writing {path}: {e}"))
}

/// The statements of a `{out}.trajectory_prior.tsv` (its `statement` rows),
/// so a prior can be read back exactly as it was used.
pub(crate) fn parse_prior_tsv(text: &str, origin: &str) -> Result<Vec<Statement>> {
    let mut out = Vec::new();
    for (n, line) in text.lines().enumerate() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.first() != Some(&"statement") {
            continue;
        }
        ensure!(f.len() >= 5, "{origin}:{}: malformed prior row", n + 1);
        let relation = Relation::parse(f[3])
            .ok_or_else(|| anyhow::anyhow!("{origin}:{}: bad relation {:?}", n + 1, f[3]))?;
        out.push(Statement {
            from: f[1].to_string(),
            to: f[2].to_string(),
            relation,
            source: Source::Run,
            note: f.get(5).unwrap_or(&"").to_string(),
        });
    }
    Ok(out)
}

#[cfg(test)]
#[path = "prior_tests.rs"]
mod tests;
