//! The prior: which cell types precede which (`docs/trajectory-plan.md` §3).
//!
//! Statements come from the Cell Ontology's `develops_from` relations and from
//! `precedence.tsv` files along lupin's data-file search path, a later layer's
//! statement about a pair of types replacing earlier ones. The `precedes`
//! statements are closed transitively over every type, restricted to the node
//! types and reduced to the direct edges, which must form a DAG.

use crate::annotate::celltype_tree::ClTerms;
use anyhow::{bail, ensure, Context, Result};
use legume_numeric::matrix::graph::{connected_components, AdjListGraph};
use log::info;
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::path::Path;

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
                    note: format!("{b} develops from {a}"),
                });
            }
        }
    }
    out
}

/// Statements from a `precedence.tsv` text, `from<TAB>to<TAB>relation[<TAB>note]`,
/// or from a `{out}.trajectory_prior.tsv` (whose rows carry a leading `kind`
/// and a source: its `statement` rows are read, its `edge` rows skipped). `#`
/// comments and a `from…` header are skipped; `origin` names the file in
/// errors.
pub(crate) fn parse_statements(text: &str, source: Source, origin: &str) -> Result<Vec<Statement>> {
    let mut out = Vec::new();
    for (n, line) in text.lines().enumerate() {
        let line = line.trim_end();
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields: Vec<&str> = line.split('\t').map(str::trim).collect();
        if fields.first() == Some(&"from")
            && fields
                .get(2)
                .is_some_and(|r| r.eq_ignore_ascii_case("relation"))
        {
            continue;
        }
        let mut note_from = 3;
        match fields.first() {
            Some(&"edge") => continue,
            Some(&"statement") => {
                fields.remove(0);
                // `from to relation source note`: the source column is dropped.
                note_from = 4;
            }
            _ => {}
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
            note: fields.get(note_from..).unwrap_or(&[]).join(" "),
        });
    }
    Ok(out)
}

/// The statements of a user- or project-layer `precedence.tsv` at `path`
/// (`layer` is `"user"` or `"project"`, as `SearchPath::user_and_project_files`
/// names them).
pub(crate) fn read_layer(path: &Path, layer: &str) -> Result<Vec<Statement>> {
    let text =
        std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let source = if layer == "user" {
        Source::User
    } else {
        Source::Project
    };
    parse_statements(&text, source, &path.display().to_string())
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

/// `--root`: each named type must be a node with no incoming statement, and
/// gets a `cli` statement to every node no statement lets it reach (other
/// forced roots aside), so a run with no other prior still orders every type
/// from it. Returns the layer to combine last.
pub(crate) fn root_layer(
    types: &[Box<str>],
    is_node: &[bool],
    statements: &[Statement],
    roots: &[&str],
) -> Result<Vec<Statement>> {
    let index = index_of(types);
    let statements: Vec<Statement> = statements
        .iter()
        .filter(|s| known(s, &index))
        .cloned()
        .collect();
    let statements = statements.as_slice();
    let adj = adjacency(types.len(), &index, statements);
    let reach = reachability(&adj);
    let mut root_idx = Vec::new();
    for &r in roots {
        let &ri = index.get(r).ok_or_else(|| {
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
        root_idx.push(ri);
    }
    // An explicit `unrelated` between the root and a type stands.
    let unrelated: BTreeSet<(usize, usize)> = statements
        .iter()
        .filter(|s| s.relation == Relation::Unrelated)
        .map(|s| {
            let (a, b) = (index[s.from.as_str()], index[s.to.as_str()]);
            (a.min(b), a.max(b))
        })
        .collect();
    let mut out = Vec::new();
    for &ri in &root_idx {
        for t in 0..types.len() {
            if t != ri
                && is_node[t]
                && !reach[ri][t]
                && !root_idx.contains(&t)
                && !unrelated.contains(&(ri.min(t), ri.max(t)))
            {
                out.push(Statement {
                    from: types[ri].to_string(),
                    to: types[t].to_string(),
                    relation: Relation::Precedes,
                    source: Source::Cli,
                    note: "--root".into(),
                });
            }
        }
    }
    Ok(out)
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
    /// `reach[a][b]`: a path of `precedes` statements leads from `a` to `b`.
    reach: Vec<Vec<bool>>,
    /// Each node's prior component, `None` for a node with no edge.
    pub(crate) component: Vec<Option<usize>>,
    /// Nodes with an edge but no incoming one, per component.
    pub(crate) roots: Vec<Vec<usize>>,
}

impl Prior {
    /// Whether the prior orders `a` and `b`, either way.
    pub(crate) fn related(&self, a: usize, b: usize) -> bool {
        self.reach[a][b] || self.reach[b][a]
    }

    /// Every root-to-leaf path of the direct edges, with its component,
    /// sorted.
    pub(crate) fn lineages(&self) -> Vec<(usize, Vec<usize>)> {
        let mut out = Vec::new();
        for (c, roots) in self.roots.iter().enumerate() {
            for &r in roots {
                let mut stack = vec![vec![r]];
                while let Some(path) = stack.pop() {
                    let last = *path.last().expect("a path has a start");
                    let next: Vec<usize> = self
                        .edges
                        .iter()
                        .filter(|e| e.from == last)
                        .map(|e| e.to)
                        .collect();
                    if next.is_empty() {
                        out.push((c, path));
                    } else {
                        for to in next {
                            let mut p = path.clone();
                            p.push(to);
                            stack.push(p);
                        }
                    }
                }
            }
        }
        out.sort();
        out
    }
}

/// Build the prior over `types` (every label in the run), of which `is_node`
/// marks the node types. Statements about types not in the run are set aside
/// with a note.
pub(crate) fn build(
    types: &[Box<str>],
    is_node: &[bool],
    statements: Vec<Statement>,
) -> Result<Prior> {
    let index = index_of(types);
    let n = types.len();
    let statements: Vec<Statement> = statements
        .into_iter()
        .filter(|s| {
            let k = known(s, &index);
            if !k {
                info!(
                    "prior statement about a type not in this run is set aside: {} {} {}",
                    s.from,
                    s.relation.as_str(),
                    s.to
                );
            }
            k
        })
        .collect();
    let adj = adjacency(n, &index, &statements);
    let reach = reachability(&adj);
    let stated: BTreeMap<(usize, usize), Source> = statements
        .iter()
        .filter(|s| s.relation == Relation::Precedes)
        .map(|s| ((index[s.from.as_str()], index[s.to.as_str()]), s.source))
        .collect();

    let mutual = (0..n)
        .flat_map(|a| (0..n).map(move |b| (a, b)))
        .find(|&(a, b)| a != b && reach[a][b] && reach[b][a]);
    if let Some((a, b)) = mutual {
        let cycle: Vec<&str> = path(&adj, a, b)
            .into_iter()
            .chain(path(&adj, b, a).into_iter().skip(1))
            .map(|i| types[i].as_ref())
            .collect();
        let involved: Vec<String> = statements
            .iter()
            .filter(|s| {
                s.relation == Relation::Precedes
                    && cycle.contains(&s.from.as_str())
                    && cycle.contains(&s.to.as_str())
            })
            .map(|s| format!("{} precedes {} ({})", s.from, s.to, s.source.as_str()))
            .collect();
        bail!(
            "the prior has a cycle: {}; from {}",
            cycle.join(" → "),
            involved.join("; ")
        );
    }

    let nodes: Vec<usize> = (0..n).filter(|&i| is_node[i]).collect();
    let mut edges = Vec::new();
    for &a in &nodes {
        for &b in nodes.iter().filter(|&&b| b != a && reach[a][b]) {
            let implied = nodes
                .iter()
                .any(|&c| c != a && c != b && reach[a][c] && reach[c][b]);
            if implied {
                continue;
            }
            let source = stated.get(&(a, b)).copied();
            let via = if source.is_some() {
                Vec::new()
            } else {
                let p = path(&adj, a, b);
                p[1..p.len() - 1].to_vec()
            };
            edges.push(Edge {
                from: a,
                to: b,
                source,
                via,
            });
        }
    }

    // Components over the direct edges, numbered densely in order of their
    // lowest node; a node with no edge is in none.
    let pairs: Vec<(usize, usize)> = edges.iter().map(|e| (e.from, e.to)).collect();
    let labels = connected_components(&AdjListGraph::from_unweighted_edges(n, &pairs));
    let mut dense: BTreeMap<usize, usize> = BTreeMap::new();
    let mut component = vec![None; n];
    for i in 0..n {
        if pairs.iter().any(|&(a, b)| a == i || b == i) {
            let next = dense.len();
            component[i] = Some(*dense.entry(labels[i]).or_insert(next));
        }
    }
    let n_components = dense.len();
    let mut roots = vec![Vec::new(); n_components];
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
        reach,
        component,
        roots,
    })
}

/// Whether both of a statement's types are in the run.
fn known(s: &Statement, index: &BTreeMap<&str, usize>) -> bool {
    index.contains_key(s.from.as_str()) && index.contains_key(s.to.as_str())
}

fn index_of(types: &[Box<str>]) -> BTreeMap<&str, usize> {
    types
        .iter()
        .enumerate()
        .map(|(i, t)| (t.as_ref(), i))
        .collect()
}

/// Successors of each type by the `precedes` statements.
fn adjacency(n: usize, index: &BTreeMap<&str, usize>, statements: &[Statement]) -> Vec<Vec<usize>> {
    let mut adj = vec![Vec::new(); n];
    for s in statements
        .iter()
        .filter(|s| s.relation == Relation::Precedes)
    {
        adj[index[s.from.as_str()]].push(index[s.to.as_str()]);
    }
    adj
}

/// `reach[a][b]`: `b` is reachable from `a` by one or more steps.
fn reachability(adj: &[Vec<usize>]) -> Vec<Vec<bool>> {
    (0..adj.len())
        .map(|start| {
            let mut seen = vec![false; adj.len()];
            let mut queue: VecDeque<usize> = adj[start].iter().copied().collect();
            while let Some(x) = queue.pop_front() {
                if !seen[x] {
                    seen[x] = true;
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
    let mut seen = vec![false; adj.len()];
    seen[a] = true;
    while let Some(x) = queue.pop_front() {
        for &y in &adj[x] {
            if !seen[y] {
                seen[y] = true;
                parent[y] = Some(x);
                if y == b {
                    let mut p = vec![b];
                    while let Some(q) = parent[*p.last().expect("non-empty")] {
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
/// reads its statements back).
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

#[cfg(test)]
#[path = "prior_tests.rs"]
mod tests;
