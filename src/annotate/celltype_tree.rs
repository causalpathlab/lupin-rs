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

use crate::annotate::cl_rules::{Aliases, MatchRules};
use crate::annotate::markers::label_key;
use data_beans::alg::union_find::UnionFind;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap};

/// The Cell Ontology terms a panel is placed on: names, the synonyms the
/// [`MatchRules`] count, and `is_a` parents, from an OBO file; plus curated
/// [`Aliases`] checked before any name.
#[derive(Default)]
pub struct ClTerms {
    rules: MatchRules,
    aliases: Aliases,
    /// Normalised name or counted synonym → term ids.
    by_name: HashMap<String, BTreeSet<String>>,
    /// The same names by [`word_bag`], for labels whose words come in
    /// another order (`B cells memory` for `memory B cell`).
    by_words: HashMap<String, BTreeSet<String>>,
    name_of: HashMap<String, String>,
    parents: HashMap<String, Vec<String>>,
    /// `is_a` children, the reverse of `parents`.
    children: HashMap<String, Vec<String>>,
    /// `develops_from` (RO:0002202): the terms a term's cells develop from.
    develops_from: HashMap<String, Vec<String>>,
    /// Terms in the analysis subsets the rules name (`ClassRules`).
    classes: BTreeSet<String>,
    /// The file's `data-version`, when it states one.
    pub release: Option<String>,
}

/// One `[Term]` stanza while it is read.
#[derive(Default)]
struct Stanza {
    id: String,
    /// The name first, then the synonyms the rules count.
    names: Vec<String>,
    parents: Vec<String>,
    develops_from: Vec<String>,
    subsets: Vec<String>,
    obsolete: bool,
}

impl Stanza {
    fn add_to(self, terms: &mut ClTerms) {
        if self.id.is_empty() || self.obsolete {
            return;
        }
        for n in &self.names {
            let n = terms.rules.normalise(n);
            if terms.rules.any_word_order {
                terms
                    .by_words
                    .entry(word_bag(&terms.rules.singular(&n)))
                    .or_default()
                    .insert(self.id.clone());
            }
            terms.by_name.entry(n).or_default().insert(self.id.clone());
        }
        if let Some(first) = self.names.first() {
            terms.name_of.insert(self.id.clone(), first.clone());
        }
        if terms.rules.classes.is_class(&self.subsets) {
            terms.classes.insert(self.id.clone());
        }
        if !self.develops_from.is_empty() {
            terms
                .develops_from
                .insert(self.id.clone(), self.develops_from);
        }
        terms.parents.insert(self.id, self.parents);
    }
}

/// A normalised name's words, hyphenated ones split, in sorted order.
fn word_bag(normalised: &str) -> String {
    let mut words: Vec<&str> = normalised
        .split([' ', '-'])
        .filter(|w| !w.is_empty())
        .collect();
    words.sort_unstable();
    words.join(" ")
}

impl ClTerms {
    /// Parse the `[Term]` stanzas of an OBO file under `rules`: `id`, `name`,
    /// the synonyms the rules count, `is_a` parents and subsets; obsolete
    /// terms are skipped.
    #[must_use]
    pub fn parse(obo: &str, rules: &MatchRules) -> Self {
        let mut terms = ClTerms {
            rules: rules.clone(),
            ..ClTerms::default()
        };
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
                    // `SCOPE [TYPE…] [xrefs]`: the rules say which count.
                    let mut words = rest.split_whitespace();
                    let scope = words.next().unwrap_or_default();
                    let types: Vec<&str> = words.take_while(|w| !w.starts_with('[')).collect();
                    if terms.rules.counts(scope, &types) {
                        t.names.push(text.to_string());
                    }
                }
            } else if let Some(v) = line.strip_prefix("is_a:") {
                if let Some(p) = v.split_whitespace().next() {
                    t.parents.push(p.to_string());
                }
            } else if let Some(v) = line.strip_prefix("relationship:") {
                // `develops_from` is written by its RO id in cl-basic.
                let mut words = v.split_whitespace();
                if let (Some("RO:0002202" | "develops_from"), Some(origin)) =
                    (words.next(), words.next())
                {
                    t.develops_from.push(origin.to_string());
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
        for (child, parents) in &terms.parents {
            for p in parents {
                terms
                    .children
                    .entry(p.clone())
                    .or_default()
                    .push(child.clone());
            }
        }
        for kids in terms.children.values_mut() {
            kids.sort_by(|a, b| {
                let name = |id: &String| terms.name_of.get(id).cloned().unwrap_or_default();
                name(a).cmp(&name(b))
            });
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

    /// Match panel labels to terms by name or exact synonym, as written or
    /// with plural words made singular, else by the same words in any order
    /// (hyphens as spaces). A label that matches no term, or several, stays
    /// unmapped rather than guessed.
    #[must_use]
    pub fn map_labels<'a>(
        &self,
        labels: impl IntoIterator<Item = &'a str>,
    ) -> (BTreeMap<String, String>, Vec<String>) {
        let mut mapped = BTreeMap::new();
        let mut unmapped = Vec::new();
        for label in labels {
            // A curated alias first, when its term is in this ontology.
            if let Some(id) = self.aliases.get(label).filter(|id| self.has(id)) {
                mapped.insert(label.to_string(), id.to_string());
                continue;
            }
            // As written, else with plural words made singular.
            let exact = self.rules.normalise(label);
            let single = self.rules.singular(&exact);
            let found = self
                .by_name
                .get(&exact)
                .or_else(|| self.by_name.get(&single))
                .or_else(|| self.by_words.get(&word_bag(&single)));
            match found {
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

    /// Check `aliases` before any name when matching labels.
    #[must_use]
    pub fn with_aliases(mut self, aliases: Aliases) -> Self {
        self.aliases = aliases;
        self
    }

    /// Layer an alias table (`label<TAB>CL:id` rows, read from `source`)
    /// over the aliases, as a run's `--label-cl` is; how many rows it took.
    pub fn add_aliases(&mut self, text: &str, source: &str) -> usize {
        self.aliases.add_tsv(text, source)
    }

    /// How far below the ontology's top `id` is: its ancestors, itself
    /// included.
    #[must_use]
    pub fn depth(&self, id: &str) -> usize {
        self.ancestors_or_self(id).len()
    }

    /// Whether `id` is a term of the ontology.
    #[must_use]
    pub fn has(&self, id: &str) -> bool {
        self.name_of.contains_key(id)
    }

    /// The terms `id` develops from (`develops_from`, RO:0002202).
    #[must_use]
    pub fn develops_from(&self, id: &str) -> &[String] {
        self.develops_from.get(id).map_or(&[], Vec::as_slice)
    }

    /// `id`'s `is_a` parents.
    #[must_use]
    pub fn parents(&self, id: &str) -> &[String] {
        self.parents.get(id).map_or(&[], Vec::as_slice)
    }

    /// `id`'s `is_a` children, by name.
    #[must_use]
    pub fn children(&self, id: &str) -> &[String] {
        self.children.get(id).map_or(&[], Vec::as_slice)
    }

    /// `id` and a chain of first-listed parents up to a term with none, top first.
    #[must_use]
    pub fn lineage(&self, id: &str) -> Vec<String> {
        let mut out = vec![id.to_string()];
        while let Some(p) = self.parents(out.last().map_or(id, String::as_str)).first() {
            if out.contains(p) {
                break;
            }
            out.push(p.clone());
        }
        out.reverse();
        out
    }

    /// Terms whose name, exact synonym or abbreviation contains `query`
    /// (compared as [`normalise`]d), those matching whole first, then by
    /// name length; at most `limit`.
    #[must_use]
    pub fn search(&self, query: &str, limit: usize) -> Vec<String> {
        let q = self.rules.normalise(query);
        if q.is_empty() {
            return Vec::new();
        }
        let mut hits: BTreeMap<String, (bool, usize)> = BTreeMap::new();
        for (name, ids) in &self.by_name {
            if !name.contains(&q) {
                continue;
            }
            for id in ids {
                let len = self.name_of.get(id).map_or(usize::MAX, String::len);
                let e = hits.entry(id.clone()).or_insert((false, len));
                e.0 |= *name == q;
            }
        }
        let mut v: Vec<(String, (bool, usize))> = hits.into_iter().collect();
        v.sort_by(|a, b| {
            b.1 .0
                .cmp(&a.1 .0)
                .then(a.1 .1.cmp(&b.1 .1))
                .then(a.0.cmp(&b.0))
        });
        v.into_iter().take(limit).map(|(id, _)| id).collect()
    }

    /// Whether `id` is one of the ontology's analysis classes.
    #[must_use]
    pub(crate) fn is_class(&self, id: &str) -> bool {
        self.classes.contains(id)
    }

    /// `id` and every term above it by `is_a`.
    pub(crate) fn ancestors_or_self(&self, id: &str) -> BTreeSet<String> {
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
    #[cfg(test)]
    #[must_use]
    pub fn group_of(&self, label: &str) -> Option<&str> {
        self.index().get(&label_key(label)).copied()
    }

    /// Each member type's group, keyed by its [`label_key`].
    #[must_use]
    pub fn index(&self) -> HashMap<String, &str> {
        self.groups
            .iter()
            .flat_map(|g| {
                g.members
                    .iter()
                    .map(move |m| (label_key(m), g.name.as_str()))
            })
            .collect()
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
        // Per analysis class: how many panel types sit under it, and its depth.
        let mut sharers: HashMap<&str, usize> = HashMap::new();
        for a in ancestry.values() {
            for t in a.iter().filter(|t| terms.classes.contains(*t)) {
                *sharers.entry(t.as_str()).or_default() += 1;
            }
        }
        let depth: HashMap<&str, usize> = sharers.keys().map(|t| (*t, terms.depth(t))).collect();
        let mut groups: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for (label, a) in &ancestry {
            let node = a
                .iter()
                .filter(|t| sharers.get(t.as_str()).is_some_and(|&n| n >= 2))
                .max_by(|x, y| {
                    depth[x.as_str()]
                        .cmp(&depth[y.as_str()])
                        .then_with(|| y.cmp(x))
                })
                .cloned()
                .unwrap_or_else(|| mapped[*label].clone());
            groups.entry(node).or_default().push((*label).to_string());
        }
        let mut groups: Vec<TypeGroup> = groups
            .into_iter()
            .map(|(id, members)| {
                // A class every sharer left for a nearer one groups nothing:
                // its lone type is named for itself, not for the broader class.
                if members.len() == 1 {
                    (mapped[&members[0]].clone(), members)
                } else {
                    (id, members)
                }
            })
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
        let mut uf = UnionFind::new(n);
        for i in 0..n {
            for j in (i + 1)..n {
                if !sets[&types[i]].is_disjoint(&sets[&types[j]]) {
                    uf.union(i, j);
                }
            }
        }
        let mut comps: BTreeMap<usize, Vec<String>> = BTreeMap::new();
        for (i, t) in types.iter().enumerate() {
            comps.entry(uf.find(i)).or_default().push(t.clone());
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
pub(crate) fn panel_types(panel: &[(String, String)]) -> Vec<String> {
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
            Some((i, _)) => {
                groups[i].members.push(t);
                // A group no ontology term names is named for its members.
                if groups[i].cl_id.is_none() {
                    groups[i].name = group_name(&groups[i].members);
                }
            }
            None => groups.push(TypeGroup {
                name: t.clone(),
                cl_id: None,
                members: vec![t],
            }),
        }
    }
}

/// A cluster's coarse call from its evidence: per group (via
/// [`TypeTree::index`]), the cluster's shares for the group's types summed,
/// when there are any (`types` and their `values`), else the cells' fine
/// labels counted per group; the group with the most called.
#[must_use]
pub fn coarse_call(
    index: &HashMap<String, &str>,
    probs: Option<(&[String], &[f32])>,
    cell_labels: &[&str],
) -> Option<String> {
    let group = |t: &str| index.get(&label_key(t)).copied();
    let mut mass: BTreeMap<&str, f32> = BTreeMap::new();
    match probs.filter(|(_, v)| v.iter().any(|&v| v > 0.0)) {
        Some((types, values)) => {
            for (t, v) in types.iter().zip(values) {
                if let Some(g) = group(t) {
                    *mass.entry(g).or_default() += v.max(0.0);
                }
            }
        }
        None => {
            // Distinct labels first, so each is looked up once.
            let mut counts: BTreeMap<&str, f32> = BTreeMap::new();
            for l in cell_labels {
                *counts.entry(l).or_default() += 1.0;
            }
            for (l, n) in counts {
                if let Some(g) = group(l) {
                    *mass.entry(g).or_default() += n;
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
