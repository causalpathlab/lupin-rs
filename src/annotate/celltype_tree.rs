//! A coarse level over a marker panel's cell types, for the first round of
//! annotation: each fine type belongs to one broad group, and a cluster's
//! first call is the group its evidence adds up to.
//!
//! The groups come from the Cell Ontology when there is one: the panel's
//! types are matched to CL terms by name or exact synonym, and each goes
//! under its nearest ancestor among the ontology's analysis classes (its
//! curated `*_upper_slim` subsets and `cellxgene_subset`) that it shares with
//! another panel type. A fine panel collapses to lineage classes; a broad
//! one stays about as it is. Without an ontology, or when too few types match, types
//! that share marker genes are grouped together. Nothing here names a cell
//! type; every group is derived from the panel and the ontology at run time.
//!
//! Pure: [`crate::manifest`] finds the ontology and the evidence.

use crate::annotate::markers::label_key;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap};

/// The Cell Ontology terms a panel is placed on: names, exact synonyms and
/// `is_a` parents, from an OBO file.
#[derive(Default)]
pub struct ClTerms {
    /// Normalised name or exact synonym → term ids.
    by_name: HashMap<String, BTreeSet<String>>,
    name_of: HashMap<String, String>,
    parents: HashMap<String, Vec<String>>,
    /// Terms in the ontology's analysis subsets: its curated upper-level
    /// slims (`*_upper_slim`) and `cellxgene_subset`, less the abstract
    /// `upper_level` terms.
    classes: BTreeSet<String>,
    /// The file's `data-version`, when it states one.
    pub release: Option<String>,
}

/// One `[Term]` stanza while it is read.
#[derive(Default)]
struct Stanza {
    id: String,
    /// The name first, then exact synonyms.
    names: Vec<String>,
    parents: Vec<String>,
    subsets: Vec<String>,
    obsolete: bool,
}

impl Stanza {
    fn add_to(self, terms: &mut ClTerms) {
        if self.id.is_empty() || self.obsolete {
            return;
        }
        for n in &self.names {
            terms
                .by_name
                .entry(normalise(n))
                .or_default()
                .insert(self.id.clone());
        }
        if let Some(first) = self.names.first() {
            terms.name_of.insert(self.id.clone(), first.clone());
        }
        let upper = self.subsets.iter().any(|s| s.ends_with("upper_level"));
        let class = |s: &String| s.ends_with("_upper_slim") || s == "cellxgene_subset";
        if !upper && self.subsets.iter().any(class) {
            terms.classes.insert(self.id.clone());
        }
        terms.parents.insert(self.id, self.parents);
    }
}

/// How labels and CL names are compared: case, the separators
/// [`label_key`] ignores (whitespace, `,`, `_`) and a plural `cells` do not
/// matter.
fn normalise(s: &str) -> String {
    let mut s = label_key(&s.to_lowercase()).replace('_', " ");
    if let Some(stem) = s.strip_suffix(" cells") {
        s = format!("{stem} cell");
    }
    s
}

impl ClTerms {
    /// Parse the `[Term]` stanzas of an OBO file: `id`, `name`, exact
    /// synonyms, `is_a` parents and subsets; obsolete terms are skipped.
    #[must_use]
    pub fn parse(obo: &str) -> Self {
        let mut terms = ClTerms::default();
        let mut stanza: Option<Stanza> = None;
        for line in obo.lines() {
            let line = line.trim();
            if line.starts_with('[') {
                if let Some(t) = stanza.take() {
                    t.add_to(&mut terms);
                }
                stanza = (line == "[Term]").then(Stanza::default);
                continue;
            }
            let Some(t) = stanza.as_mut() else {
                if let Some(v) = line.strip_prefix("data-version:") {
                    terms.release = Some(v.trim().to_string());
                }
                continue;
            };
            if let Some(v) = line.strip_prefix("id:") {
                t.id = v.trim().to_string();
            } else if let Some(v) = line.strip_prefix("name:") {
                // The name goes first, so it is the one reported.
                t.names.insert(0, v.trim().to_string());
            } else if let Some(v) = line.strip_prefix("synonym:") {
                let mut parts = v.trim().splitn(3, '"');
                if let (Some(""), Some(text), Some(rest)) =
                    (parts.next(), parts.next(), parts.next())
                {
                    if rest.trim_start().starts_with("EXACT") {
                        t.names.push(text.to_string());
                    }
                }
            } else if let Some(v) = line.strip_prefix("is_a:") {
                if let Some(p) = v.split_whitespace().next() {
                    t.parents.push(p.to_string());
                }
            } else if let Some(v) = line.strip_prefix("subset:") {
                t.subsets.push(v.trim().to_string());
            } else if line == "is_obsolete: true" {
                t.obsolete = true;
            }
        }
        if let Some(t) = stanza {
            t.add_to(&mut terms);
        }
        terms
    }

    #[cfg(test)]
    fn term_count(&self) -> usize {
        self.parents.len()
    }

    #[must_use]
    pub fn name(&self, id: &str) -> Option<&str> {
        self.name_of.get(id).map(String::as_str)
    }

    /// Match panel labels to terms by name or exact synonym. A label that
    /// matches no term, or several, stays unmapped rather than guessed.
    #[must_use]
    pub fn map_labels<'a>(
        &self,
        labels: impl IntoIterator<Item = &'a str>,
    ) -> (BTreeMap<String, String>, Vec<String>) {
        let mut mapped = BTreeMap::new();
        let mut unmapped = Vec::new();
        for label in labels {
            match self.by_name.get(&normalise(label)) {
                Some(ids) if ids.len() == 1 => {
                    mapped.insert(
                        label.to_string(),
                        ids.iter().next().cloned().unwrap_or_default(),
                    );
                }
                _ => unmapped.push(label.to_string()),
            }
        }
        (mapped, unmapped)
    }

    /// `id` and every term above it by `is_a`.
    fn ancestors_or_self(&self, id: &str) -> BTreeSet<String> {
        let mut seen = BTreeSet::new();
        let mut stack = vec![id.to_string()];
        while let Some(t) = stack.pop() {
            if seen.insert(t.clone()) {
                stack.extend(self.parents.get(&t).into_iter().flatten().cloned());
            }
        }
        seen
    }
}

/// Where a [`TypeTree`]'s groups came from.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum TreeSource {
    CellOntology,
    MarkerSharing,
}

/// One broad group of fine cell types.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct TypeGroup {
    pub name: String,
    /// The CL term the group is, for an ontology group.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cl_id: Option<String>,
    pub members: Vec<String>,
}

/// The coarse level over a panel's cell types.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct TypeTree {
    pub source: TreeSource,
    /// The ontology release, for an ontology tree.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub release: Option<String>,
    /// Panel labels matched to CL terms.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub label_cl: BTreeMap<String, String>,
    pub groups: Vec<TypeGroup>,
}

impl TypeTree {
    /// The group a fine type belongs to.
    #[must_use]
    pub fn group_of(&self, label: &str) -> Option<&str> {
        let key = label_key(label);
        self.groups
            .iter()
            .find(|g| g.members.iter().any(|m| label_key(m) == key))
            .map(|g| g.name.as_str())
    }

    /// The panel's types on the Cell Ontology. `None` when fewer than two
    /// types match a term, which leaves nothing to group.
    #[must_use]
    pub fn from_ontology(terms: &ClTerms, panel: &[(String, String)]) -> Option<Self> {
        let types = panel_types(panel);
        let (mapped, unmapped) = terms.map_labels(types.iter().map(String::as_str));
        if mapped.len() < 2 {
            return None;
        }
        let ancestry: BTreeMap<&str, BTreeSet<String>> = mapped
            .iter()
            .map(|(label, id)| (label.as_str(), terms.ancestors_or_self(id)))
            .collect();
        // Each type goes under its nearest ancestor (or itself) among the
        // ontology's analysis classes that it shares with another panel type.
        // A type sharing none stands as its own group.
        let depth = |t: &str| terms.ancestors_or_self(t).len();
        let shared = |t: &str| ancestry.values().filter(|a| a.contains(t)).count() >= 2;
        let mut groups: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for (label, a) in &ancestry {
            let node = a
                .iter()
                .filter(|t| terms.classes.contains(*t) && shared(t))
                .max_by(|x, y| depth(x).cmp(&depth(y)).then_with(|| y.cmp(x)))
                .cloned()
                .unwrap_or_else(|| mapped[*label].clone());
            groups.entry(node).or_default().push((*label).to_string());
        }
        // A class every sharer left for a nearer one groups nothing: its lone
        // type is named for itself, not for the broader class.
        let lone: Vec<(String, String)> = groups
            .iter()
            .filter(|(node, m)| m.len() == 1 && mapped[&m[0]] != **node)
            .map(|(node, m)| (node.clone(), m[0].clone()))
            .collect();
        for (node, label) in lone {
            groups.remove(&node);
            groups
                .entry(mapped[&label].clone())
                .or_default()
                .push(label);
        }
        let mut groups: Vec<TypeGroup> = groups
            .into_iter()
            .map(|(id, members)| TypeGroup {
                name: label_key(terms.name(&id).unwrap_or(&id)),
                cl_id: Some(id),
                members,
            })
            .collect();
        attach_by_sharing(&mut groups, unmapped, panel);
        Some(TypeTree {
            source: TreeSource::CellOntology,
            release: terms.release.clone(),
            label_cl: mapped,
            groups,
        })
    }

    /// Groups of types that share marker genes: the connected components of
    /// "has a marker in common". When everything connects, the component is
    /// split by average linkage on marker overlap at its widest gap.
    #[must_use]
    pub fn from_markers(panel: &[(String, String)]) -> Self {
        let types = panel_types(panel);
        let sets = marker_sets(panel);
        let n = types.len();
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(p: &mut [usize], i: usize) -> usize {
            if p[i] != i {
                let r = find(p, p[i]);
                p[i] = r;
            }
            p[i]
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if jaccard(&sets[&types[i]], &sets[&types[j]]) > 0.0 {
                    let (a, b) = (find(&mut parent, i), find(&mut parent, j));
                    parent[a] = b;
                }
            }
        }
        let mut comps: BTreeMap<usize, Vec<String>> = BTreeMap::new();
        for (i, t) in types.iter().enumerate() {
            let r = find(&mut parent, i);
            comps.entry(r).or_default().push(t.clone());
        }
        let mut members: Vec<Vec<String>> = comps.into_values().collect();
        if members.len() == 1 && n >= 3 {
            members = split_by_linkage(&types, &sets);
        }
        TypeTree {
            source: TreeSource::MarkerSharing,
            release: None,
            label_cl: BTreeMap::new(),
            groups: members
                .into_iter()
                .map(|m| TypeGroup {
                    name: group_name(&m),
                    cl_id: None,
                    members: m,
                })
                .collect(),
        }
    }
}

/// The panel's cell types, in first-seen order, one per scoring name.
fn panel_types(panel: &[(String, String)]) -> Vec<String> {
    let mut seen = BTreeSet::new();
    panel
        .iter()
        .filter(|(_, t)| seen.insert(label_key(t)))
        .map(|(_, t)| t.trim().to_string())
        .collect()
}

fn marker_sets(panel: &[(String, String)]) -> BTreeMap<String, BTreeSet<String>> {
    let names: BTreeMap<String, String> = panel_types(panel)
        .into_iter()
        .map(|t| (label_key(&t), t))
        .collect();
    let mut sets: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for (g, t) in panel {
        let t = names[&label_key(t)].clone();
        sets.entry(t).or_default().insert(g.trim().to_string());
    }
    sets
}

fn jaccard(a: &BTreeSet<String>, b: &BTreeSet<String>) -> f64 {
    let inter = a.intersection(b).count();
    let union = a.union(b).count();
    if union == 0 {
        0.0
    } else {
        inter as f64 / union as f64
    }
}

/// Average-linkage agglomeration on `1 - jaccard`, cut where successive
/// merge heights jump the most.
fn split_by_linkage(
    types: &[String],
    sets: &BTreeMap<String, BTreeSet<String>>,
) -> Vec<Vec<String>> {
    let mut clusters: Vec<Vec<usize>> = (0..types.len()).map(|i| vec![i]).collect();
    let dist = |a: &[usize], b: &[usize]| {
        let total: f64 = a
            .iter()
            .flat_map(|&i| b.iter().map(move |&j| (i, j)))
            .map(|(i, j)| 1.0 - jaccard(&sets[&types[i]], &sets[&types[j]]))
            .sum();
        total / (a.len() * b.len()) as f64
    };
    // Record every partition and the height of the merge that ends it.
    let mut history: Vec<(f64, Vec<Vec<usize>>)> = Vec::new();
    while clusters.len() > 1 {
        let mut best = (f64::INFINITY, 0, 1);
        for i in 0..clusters.len() {
            for j in (i + 1)..clusters.len() {
                let d = dist(&clusters[i], &clusters[j]);
                if d < best.0 {
                    best = (d, i, j);
                }
            }
        }
        history.push((best.0, clusters.clone()));
        let merged = clusters.remove(best.2);
        clusters[best.1].extend(merged);
    }
    // The partition just before the largest jump in merge height, keeping at
    // least two groups.
    let cut = history
        .windows(2)
        .enumerate()
        .max_by(|(_, a), (_, b)| (a[1].0 - a[0].0).total_cmp(&(b[1].0 - b[0].0)))
        .map_or(history.len().saturating_sub(1), |(i, _)| i + 1);
    let (_, partition) = &history[cut.min(history.len() - 1)];
    partition
        .iter()
        .map(|c| c.iter().map(|&i| types[i].clone()).collect())
        .collect()
}

/// A marker-sharing group's name: its members, the first few spelled out.
fn group_name(members: &[String]) -> String {
    if members.len() <= 3 {
        members.join("/")
    } else {
        format!("{}/{}/+{}", members[0], members[1], members.len() - 2)
    }
}

/// Put each unmapped type in the group whose members share most markers
/// with it, or in a group of its own.
fn attach_by_sharing(
    groups: &mut Vec<TypeGroup>,
    unmapped: Vec<String>,
    panel: &[(String, String)],
) {
    let sets = marker_sets(panel);
    for t in unmapped {
        let own = &sets[&t];
        let best = groups
            .iter()
            .enumerate()
            .map(|(i, g)| {
                let union: BTreeSet<String> = g
                    .members
                    .iter()
                    .filter_map(|m| sets.get(m))
                    .flatten()
                    .cloned()
                    .collect();
                (i, jaccard(own, &union))
            })
            .filter(|(_, j)| *j > 0.0)
            .max_by(|a, b| a.1.total_cmp(&b.1));
        match best {
            Some((i, _)) => groups[i].members.push(t),
            None => groups.push(TypeGroup {
                name: t.clone(),
                cl_id: None,
                members: vec![t],
            }),
        }
    }
}

/// A cluster's coarse call from its evidence: per group, the cluster's
/// probabilities for the group's types summed, when there are any; otherwise
/// the group most of its cells' fine labels fall in.
#[must_use]
pub fn coarse_call(
    tree: &TypeTree,
    probs: Option<&[(String, f32)]>,
    cell_labels: &[&str],
) -> Option<String> {
    let mut mass: BTreeMap<&str, f32> = BTreeMap::new();
    match probs.filter(|p| p.iter().any(|(_, v)| *v > 0.0)) {
        Some(p) => {
            for (t, v) in p {
                if let Some(g) = tree.group_of(t) {
                    *mass.entry(g).or_default() += v.max(0.0);
                }
            }
        }
        None => {
            for l in cell_labels {
                if let Some(g) = tree.group_of(l) {
                    *mass.entry(g).or_default() += 1.0;
                }
            }
        }
    }
    mass.into_iter()
        .filter(|(_, m)| *m > 0.0)
        .max_by(|a, b| a.1.total_cmp(&b.1).then_with(|| b.0.cmp(a.0)))
        .map(|(g, _)| g.to_string())
}

#[cfg(test)]
#[path = "celltype_tree_tests.rs"]
mod tests;
