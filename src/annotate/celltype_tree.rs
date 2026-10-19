//! A coarse level over a marker panel's cell types, for the first round of
//! annotation: each fine type belongs to one broad group, and a cluster's
//! first call is the group its evidence adds up to.
//!
//! The groups come from the Cell Ontology when there is one: the panel's
//! types are matched to CL terms under the matching rules, and each goes
//! under its nearest ancestor among the ontology's analysis classes (its
//! curated `*_upper_slim` subsets and `cellxgene_subset`) that it shares with
//! another panel type; a type that matches no term joins the group sharing
//! its most specific markers. A fine panel collapses to lineage classes; a
//! broad one stays about as it is. Without an ontology, or when too few
//! types match, types that share specific markers are grouped together.
//! Nothing here names a cell type; every group is derived from the panel and
//! the ontology at run time.
//!
//! Pure: [`crate::manifest`] finds the ontology and the evidence.

use crate::annotate::cl_rules::{Aliases, MatchRules};
use crate::annotate::markers::label_key;
use data_beans::alg::union_find::UnionFind;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::sync::OnceLock;

/// The Cell Ontology terms a panel is placed on: names, the synonyms the
/// [`MatchRules`] count, and `is_a` parents, from an OBO file; plus curated
/// [`Aliases`] checked before any name.
#[derive(Default)]
pub struct ClTerms {
    rules: MatchRules,
    aliases: Aliases,
    /// The aliases by label as names are compared ([`MatchRules::normalise`]).
    alias_of: HashMap<String, String>,
    /// The user's answers on types left alone (see [`JoinQuestion`]): join
    /// (`true`) or stand alone, by label as names are compared.
    joins: HashMap<String, bool>,
    /// The types the user found no fitting term for, likewise.
    declined: HashSet<String>,
    /// [`MatchRules::qualifier_words`].
    qualifiers: Vec<Vec<String>>,
    /// Normalised name or counted synonym → term ids.
    by_name: HashMap<String, BTreeSet<String>>,
    /// The same names by [`ClTerms::bag`]: for labels whose words come in
    /// another order (`CT1 cells memory` for `memory CT1 cell`), when the
    /// rules allow it, and for the panel's word.
    by_words: HashMap<String, BTreeSet<String>>,
    /// Each term's own name without its [`MatchRules::qualifiers`], when it
    /// has any, by [`ClTerms::key`].
    by_unqualified: HashMap<String, BTreeSet<String>>,
    /// Each term's own name (no synonyms) by [`ClTerms::key`]: what an
    /// implied trailing word may complete.
    by_term_name: HashMap<String, BTreeSet<String>>,
    name_of: HashMap<String, String>,
    parents: HashMap<String, Vec<String>>,
    /// `is_a` children, the reverse of `parents`.
    children: HashMap<String, Vec<String>>,
    /// `develops_from` (RO:0002202): the terms a term's cells develop from.
    develops_from: HashMap<String, Vec<String>>,
    /// Terms in the analysis subsets the rules name (`ClassRules`).
    classes: BTreeSet<String>,
    /// The words the ontology uses for kinds of a term (see `PartRules`),
    /// worked out when first needed.
    kind_words: OnceLock<BTreeSet<String>>,
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
        for (i, n) in self.names.iter().enumerate() {
            let n = terms.rules.normalise(n);
            let key = terms.key(&n);
            if i == 0 {
                if let Some(core) = terms.rules.unqualified(&n, &terms.qualifiers) {
                    let core = terms.key(&core);
                    terms
                        .by_unqualified
                        .entry(core)
                        .or_default()
                        .insert(self.id.clone());
                }
                (terms.by_term_name.entry(key.clone()).or_default()).insert(self.id.clone());
            }
            if terms.rules.any_word_order {
                terms
                    .by_words
                    .entry(key)
                    .or_default()
                    .insert(self.id.clone());
            }
            terms.by_name.entry(n).or_default().insert(self.id.clone());
        }
        if let Some(first) = self.names.first() {
            terms.name_of.insert(self.id.clone(), first.clone());
        }
        if terms.rules.classes.is_class(&self.id, &self.subsets) {
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

impl ClTerms {
    /// Parse the `[Term]` stanzas of an OBO file under `rules`: `id`, `name`,
    /// the synonyms the rules count, `is_a` parents and subsets; obsolete
    /// terms are skipped.
    #[must_use]
    pub fn parse(obo: &str, rules: &MatchRules) -> Self {
        let mut terms = ClTerms {
            rules: rules.clone(),
            qualifiers: rules.qualifier_words(),
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
    /// (hyphens as spaces, but a marker keeps its sign), else with an
    /// implied trailing word completing a term's name, else as the most
    /// general term whose name it is but for its qualifiers. Where a text names several terms, a term's own name wins
    /// over another's synonym; a label still matching several, or none,
    /// stays unmapped rather than guessed.
    #[must_use]
    pub fn map_labels<'a>(
        &self,
        labels: impl IntoIterator<Item = &'a str>,
    ) -> (BTreeMap<String, String>, Vec<String>) {
        let (mapped, unmapped, _) = self.match_panel(labels);
        (mapped, unmapped)
    }

    /// [`Self::map_labels`], and the word the panel leaves off when there
    /// is one: the labels left are completed by it.
    fn match_panel<'a>(
        &self,
        labels: impl IntoIterator<Item = &'a str>,
    ) -> (BTreeMap<String, String>, Vec<String>, Option<String>) {
        let (placed, unmapped) = self.read_panel(labels);
        (placed.mapped, unmapped, placed.context)
    }

    /// [`Self::match_panel`] as a [`Placed`] holding the matched labels and
    /// what matched each, and the labels left.
    fn read_panel<'a>(&self, labels: impl IntoIterator<Item = &'a str>) -> (Placed, Vec<String>) {
        let mut placed = Placed::default();
        let mut unmapped = Vec::new();
        for label in labels {
            match self.match_how(label) {
                Some((id, how)) => placed.put(label.to_string(), id, how),
                None => unmapped.push(label.to_string()),
            }
        }
        placed.context = self.context_word(&unmapped);
        if let Some(w) = placed.context.clone() {
            unmapped.retain(|l| match self.complete(l, &w) {
                Some(id) => {
                    placed.put(l.clone(), id, How::PanelWord(w.clone()));
                    false
                }
                None => true,
            });
        }
        (placed, unmapped)
    }

    /// A panel's labels on the ontology: [`Self::map_labels`], then each
    /// label left placed by the longest run of its words
    /// ([`Self::place_with`]), the panel's word completing runs too.
    #[must_use]
    pub fn place_labels<'a>(&self, labels: impl IntoIterator<Item = &'a str>) -> Placed {
        let (mut placed, unmapped) = self.read_panel(labels);
        for label in unmapped {
            match self.place_with(&label, placed.context.as_deref()) {
                Some((id, true)) => placed.put(label, id, How::Under),
                Some((id, false)) => placed.put(label, id, How::Kept),
                None => placed.none.push(label),
            }
        }
        placed
    }

    /// One label's term: a curated alias first, when its term is in this
    /// ontology, else by the names.
    fn match_one(&self, label: &str) -> Option<String> {
        self.match_how(label).map(|(id, _)| id)
    }

    /// [`Self::match_one`], and which of the two matched.
    fn match_how(&self, label: &str) -> Option<(String, How)> {
        (self.alias(label).map(|id| (id, How::Alias)))
            .or_else(|| self.match_name(label).map(|id| (id, How::Name)))
    }

    /// The term `label` is aliased to, when this ontology has it.
    fn alias(&self, label: &str) -> Option<String> {
        let exact = self.rules.normalise(label);
        (self.alias_of.get(&exact))
            .filter(|id| self.has(id))
            .cloned()
    }

    /// One label's term by the names alone.
    fn match_name(&self, label: &str) -> Option<String> {
        let exact = self.rules.normalise(label);
        // As written, else made singular, else in any word order, else
        // with an implied trailing word completing a term's own name
        // (see `ImpliedRules`), else as the most general term it fits
        // but for qualifiers, with or without the implied word. Where a
        // text names several terms, the one it is the own name of wins
        // over those it is a synonym of.
        let single = self.rules.singular(&exact);
        let key = self.key(&exact);
        let only = |ids: &BTreeSet<String>| -> Option<String> {
            (ids.len() == 1).then(|| ids.iter().next().cloned())?
        };
        let one = |ids: &BTreeSet<String>, text: &str| -> Option<String> {
            only(ids).or_else(|| {
                let named = self.by_term_name.get(&self.key(text))?;
                only(&ids.intersection(named).cloned().collect())
            })
        };
        if let Some(ids) = self.by_name.get(&exact) {
            one(ids, &exact)
        } else if let Some(ids) = self.by_name.get(&single) {
            one(ids, &single)
        } else if let Some(ids) = (self.by_words.get(&key)).filter(|_| self.rules.any_word_order) {
            one(ids, &single)
        } else if !key.is_empty() {
            // A term's own name completed first; the looser link, a name
            // its qualifiers left off, only then.
            let completed: Vec<String> = (self.rules.implied.suffixes.iter())
                .map(|w| self.key(&self.rules.normalise(&format!("{exact} {w}"))))
                .collect();
            (completed.iter())
                .find_map(|k| self.by_term_name.get(k).and_then(only))
                .or_else(|| {
                    (std::iter::once(&key).chain(&completed))
                        .find_map(|k| self.most_general(self.by_unqualified.get(k)?))
                })
        } else {
            None
        }
    }

    /// `label`'s words as names' are compared: normalised, made singular.
    fn label_words(&self, label: &str) -> Vec<String> {
        let single = self.rules.singular(&self.rules.normalise(label));
        self.rules.words(&single).map(String::from).collect()
    }

    /// The word the panel leaves off its `labels` (those no name matched):
    /// the one that completes the most of them ([`Self::complete`]), when it
    /// completes enough and more than any other word (see `ContextRules`).
    fn context_word(&self, labels: &[String]) -> Option<String> {
        let min = self.rules.context.min_labels?;
        if labels.len() < min.max(1) {
            return None;
        }
        let implied = &self.rules.implied.suffixes;
        // Candidate words: each a name has beside all of a label's words, with
        // at most an implied one.
        let own: Vec<Vec<String>> = (labels.iter())
            .map(|l| self.label_words(l))
            .filter(|w| !w.is_empty())
            .collect();
        let mut words: BTreeSet<&str> = BTreeSet::new();
        for bag in self.by_words.keys() {
            let all: Vec<&str> = bag.split(' ').collect();
            for own in &own {
                if all.len() > own.len() + 2 || all.len() <= own.len() {
                    continue;
                }
                let mut rest = all.clone();
                let all_there = own.iter().all(|w| {
                    let at = rest.iter().position(|r| r == w);
                    at.map(|i| rest.remove(i)).is_some()
                });
                if !all_there {
                    continue;
                }
                if rest.len() == 2 {
                    rest.retain(|r| !implied.iter().any(|w| w == r));
                }
                if let [w] = rest[..] {
                    words.insert(w);
                }
            }
        }
        let mut ranked: Vec<(usize, &str)> = (words.into_iter())
            .filter(|w| !implied.iter().any(|s| s == w) && !self.rules.has_sign(w))
            .map(|w| {
                (
                    labels
                        .iter()
                        .filter(|l| self.complete(l, w).is_some())
                        .count(),
                    w,
                )
            })
            .collect();
        ranked.sort_by_key(|a| std::cmp::Reverse(a.0));
        let (n, w) = *ranked.first()?;
        let unique = ranked.get(1).is_none_or(|r| r.0 < n);
        (n >= min.max(1) && unique).then(|| w.to_string())
    }

    /// `label` with word `w`, and with `w` and an implied word, among the
    /// names (and counted synonyms) in any order: the most general term
    /// they name.
    fn complete(&self, label: &str, w: &str) -> Option<String> {
        if self.label_words(label).is_empty() {
            return None;
        }
        let bag = |text: String| self.bag(&self.rules.normalise(&text));
        let ids: BTreeSet<String> = std::iter::once(bag(format!("{label} {w}")))
            .chain((self.rules.implied.suffixes.iter()).map(|s| bag(format!("{label} {w} {s}"))))
            .filter_map(|key| self.by_words.get(&key))
            .flatten()
            .cloned()
            .collect();
        self.most_general(&ids)
    }

    /// The term the longest run of `label`'s words names ([`Self::match_one`]),
    /// when the runs of that length name one term, and whether every word
    /// left off is a word the ontology uses for kinds of a term (see
    /// `PartRules`); an implied word names nothing alone and is never in
    /// the way.
    /// A run is also completed by the panel's word `context`.
    pub(crate) fn place_with(&self, label: &str, context: Option<&str>) -> Option<(String, bool)> {
        self.rules.part.min_pairs?;
        let implied = |w: &String| self.rules.implied.suffixes.contains(w);
        let kinds = self.kind_words.get_or_init(|| self.kind_words());
        let w = self.label_words(label);
        for len in (1..w.len()).rev() {
            let mut found: BTreeMap<String, bool> = BTreeMap::new();
            for start in 0..=w.len() - len {
                let run = &w[start..start + len];
                if run.iter().all(implied) {
                    continue;
                }
                let text = run.join(" ");
                let id = (self.match_one(&text))
                    .or_else(|| context.and_then(|c| self.complete(&text, c)));
                if let Some(id) = id {
                    let mut left = w[..start].iter().chain(&w[start + len..]);
                    let kind = left.all(|d| implied(d) || kinds.contains(d));
                    *found.entry(id).or_default() |= kind;
                }
            }
            match found.len() {
                0 => {}
                1 => return found.into_iter().next(),
                _ => return None,
            }
        }
        None
    }

    /// The words the ontology uses for kinds of a term (see `PartRules`):
    /// over every pair of own names `X` and `X` with one word more, those
    /// whose added word makes a kind of `X` in enough pairs and enough of
    /// its pairs.
    fn kind_words(&self) -> BTreeSet<String> {
        let Some(min) = self.rules.part.min_pairs else {
            return BTreeSet::new();
        };
        let mut by_bag: HashMap<Vec<String>, Vec<&str>> = HashMap::new();
        for (id, name) in &self.name_of {
            let mut bag = self.label_words(name);
            bag.sort_unstable();
            by_bag.entry(bag).or_default().push(id);
        }
        // Per added word: the pairs where it makes a kind of the shorter
        // name, and the others.
        let mut pairs: HashMap<&str, (usize, usize)> = HashMap::new();
        for (bag, ids) in &by_bag {
            let mut above = None;
            for (i, word) in bag.iter().enumerate() {
                // A word twice in a name is one pair.
                if i > 0 && bag[i - 1] == *word {
                    continue;
                }
                let mut base = bag.clone();
                base.remove(i);
                let Some(bases) = by_bag.get(&base) else {
                    continue;
                };
                let above: &Vec<BTreeSet<String>> = above.get_or_insert_with(|| {
                    ids.iter().map(|id| self.ancestors_or_self(id)).collect()
                });
                let e = pairs.entry(word).or_default();
                for a in above {
                    for b in bases {
                        if a.contains(*b) {
                            e.0 += 1;
                        } else {
                            e.1 += 1;
                        }
                    }
                }
            }
        }
        let share = self.rules.part.min_share;
        (pairs.into_iter())
            .filter(|(_, (kind, other))| {
                *kind >= min && *kind as f64 >= share * (kind + other) as f64
            })
            .map(|(w, _)| w.to_string())
            .collect()
    }

    /// The one of `ids` every other is under, when there is one.
    fn most_general(&self, ids: &BTreeSet<String>) -> Option<String> {
        let under: Vec<BTreeSet<String>> =
            ids.iter().map(|id| self.ancestors_or_self(id)).collect();
        ids.iter()
            .find(|g| under.iter().all(|a| a.contains(*g)))
            .cloned()
    }

    /// A normalised text as [`Self::by_words`] keys it: made singular, its
    /// words sorted.
    fn bag(&self, normalised: &str) -> String {
        self.rules.word_bag(&self.rules.singular(normalised))
    }

    /// A normalised text as the name indexes key it: made singular, and its
    /// words sorted when the rules allow any word order.
    fn key(&self, normalised: &str) -> String {
        if self.rules.any_word_order {
            self.bag(normalised)
        } else {
            self.rules.singular(normalised)
        }
    }

    /// Every counted name (as names are compared) with the terms it names.
    pub fn named_texts(&self) -> impl Iterator<Item = (&str, &BTreeSet<String>)> {
        self.by_name.iter().map(|(name, ids)| (name.as_str(), ids))
    }

    /// A label as names are compared.
    #[must_use]
    pub fn normalise(&self, label: &str) -> String {
        self.rules.normalise(label)
    }

    /// Check `aliases` before any name when matching labels.
    #[must_use]
    pub fn with_aliases(mut self, aliases: Aliases) -> Self {
        self.aliases = aliases;
        self.index_aliases();
        self
    }

    /// Layer an alias table (`label<TAB>CL:id` rows, read from `source`)
    /// over the aliases, as a run's `--label-cl` is; how many rows it took.
    pub fn add_aliases(&mut self, text: &str, source: &str) -> usize {
        let n = self.aliases.add_tsv(text, source);
        self.index_aliases();
        n
    }

    /// Layer the user's answers on the panel's types: `label<TAB>join` or
    /// `label<TAB>alone` for a type left alone, `label<TAB>none` for one no
    /// term fits; read as every two-column table is, a later row winning.
    /// How many rows it took.
    pub fn add_joins(&mut self, text: &str) -> usize {
        let rows = legume_numeric::matrix::common_io::TwoFieldRows::from_text(text);
        let mut n = 0;
        for (label, answer) in &rows.rows {
            let key = self.rules.normalise(label);
            let join = match answer.trim().to_lowercase().as_str() {
                "join" => true,
                "alone" => false,
                "none" => {
                    self.declined.insert(key);
                    n += 1;
                    continue;
                }
                other => {
                    log::warn!("{label}: answer {other:?} is not join, alone or none; skipped");
                    continue;
                }
            };
            self.joins.insert(key, join);
            n += 1;
        }
        n
    }

    /// The user's answer for type `label` left alone, if any.
    pub(crate) fn joined(&self, label: &str) -> Option<bool> {
        self.joins.get(&self.rules.normalise(label)).copied()
    }

    /// Whether the user found no fitting term for type `label`.
    #[must_use]
    pub fn declined(&self, label: &str) -> bool {
        self.declined.contains(&self.rules.normalise(label))
    }

    /// Whether labels `a` and `b` are one label as names are compared.
    #[must_use]
    pub fn alike(&self, a: &str, b: &str) -> bool {
        self.rules.normalise(a) == self.rules.normalise(b)
    }

    fn index_aliases(&mut self) {
        self.alias_of = (self.aliases.iter())
            .map(|(label, id)| (self.rules.normalise(label), id.to_string()))
            .collect();
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
        // A query is often a word being typed (`x-` on the way to
        // `x-positive`): it is looked for as written and with its signs read.
        let q = self.rules.normalise(query);
        let typed = self.rules.normalise_unsigned(query);
        if q.is_empty() && typed.is_empty() {
            return Vec::new();
        }
        let has = |name: &str| {
            (!q.is_empty() && name.contains(&q)) || (!typed.is_empty() && name.contains(&typed))
        };
        // Ranked: a whole-name hit (any name), then by where the term's own
        // name has the query (at its start, inside, not at all: a synonym's
        // hit), then shorter names.
        let mut hits: BTreeMap<String, (bool, u8, usize)> = BTreeMap::new();
        for (name, ids) in &self.by_name {
            if !has(name) {
                continue;
            }
            for id in ids {
                let e = hits.entry(id.clone()).or_insert_with(|| {
                    let own = self.name_of.get(id);
                    let at = match own.map(|n| self.rules.normalise(n)) {
                        Some(n) if n.starts_with(&q) || n.starts_with(&typed) => 0,
                        Some(n) if has(&n) => 1,
                        _ => 2,
                    };
                    (false, at, own.map_or(usize::MAX, String::len))
                });
                e.0 |= *name == q || *name == typed;
            }
        }
        // Stable over the ids' order, which breaks the remaining ties.
        let mut v: Vec<_> = hits.into_iter().collect();
        v.sort_by_key(|(_, (exact, at, len))| (!*exact, *at, *len));
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

/// A panel's labels on the Cell Ontology ([`ClTerms::place_labels`]).
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Placed {
    /// Matched to a term.
    pub mapped: BTreeMap<String, String>,
    /// Under a term, by the longest run of their words.
    pub under: BTreeMap<String, String>,
    /// Kept to a term's lineage when they join a group by their markers.
    pub kept: BTreeMap<String, String>,
    /// Nowhere on the ontology.
    pub none: Vec<String>,
    /// The word the panel leaves off.
    pub context: Option<String>,
    /// What put each label on a term, matched, under or kept.
    by: BTreeMap<String, How>,
}

/// The fewest labels on terms for the ontology to group a panel's types.
const FEWEST_ON_TERMS: usize = 2;

/// What put a label where it is ([`Placed::how`]), for the record.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case", tag = "by", content = "word")]
pub enum How {
    /// A curated alias.
    Alias,
    /// A term's name or counted synonym, under the matching rules.
    Name,
    /// The word the panel leaves off its labels.
    PanelWord(String),
    /// Under the term the longest run of its words names.
    Under,
    /// Kept to a term's lineage, grouped by its markers.
    Kept,
    /// Nowhere on the ontology: grouped by its markers.
    Nowhere,
}

impl Placed {
    /// Each label on a term: matched, or under it.
    pub fn on(&self) -> impl Iterator<Item = (&String, &String)> {
        self.mapped.iter().chain(&self.under)
    }

    /// What put `label` where it is.
    #[must_use]
    pub fn how(&self, label: &str) -> How {
        self.by.get(label).cloned().unwrap_or(How::Nowhere)
    }

    /// The term `label` is matched to, under, or kept to.
    #[must_use]
    pub fn term(&self, label: &str) -> Option<&String> {
        (self.mapped.get(label))
            .or_else(|| self.under.get(label))
            .or_else(|| self.kept.get(label))
    }

    /// `label` on term `id`, put there as `how` says.
    fn put(&mut self, label: String, id: String, how: How) {
        let to = match how {
            How::Alias | How::Name | How::PanelWord(_) => &mut self.mapped,
            How::Under => &mut self.under,
            How::Kept => &mut self.kept,
            How::Nowhere => unreachable!("a label on term {id} is somewhere"),
        };
        to.insert(label.clone(), id);
        self.by.insert(label, how);
    }

    /// Each label on a term by its term, when enough are for the ontology
    /// to group them.
    #[must_use]
    pub fn on_terms(&self) -> Option<BTreeMap<&str, &str>> {
        let on: BTreeMap<&str, &str> = (self.on())
            .map(|(l, id)| (l.as_str(), id.as_str()))
            .collect();
        (on.len() >= FEWEST_ON_TERMS).then_some(on)
    }

    /// The labels on no term: those kept to a lineage, then the others.
    #[must_use]
    pub fn loose(&self) -> Vec<String> {
        self.kept.keys().chain(&self.none).cloned().collect()
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

/// Where a round's type tree puts a label on the Cell Ontology.
#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct TreePlacement {
    pub id: String,
    /// The panel type or group the term was matched for.
    pub name: String,
    pub via: Via,
}

/// What placed a label: its own term, or its group's.
#[derive(Serialize, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum Via {
    Label,
    Group,
}

impl TypeTree {
    /// `label`'s place on the Cell Ontology: the term it matched, else the
    /// term of the group it is (a coarse label) or belongs to.
    #[must_use]
    pub fn placement(&self, label: &str) -> Option<TreePlacement> {
        let key = label_key(label);
        if let Some((name, id)) = self.label_cl.iter().find(|(l, _)| label_key(l) == key) {
            return Some(TreePlacement {
                id: id.clone(),
                name: name.clone(),
                via: Via::Label,
            });
        }
        let group = self
            .groups
            .iter()
            .find(|g| label_key(&g.name) == key)
            .or_else(|| {
                self.groups
                    .iter()
                    .find(|g| g.members.iter().any(|m| label_key(m) == key))
            })?;
        Some(TreePlacement {
            id: group.cl_id.clone()?,
            name: group.name.clone(),
            via: Via::Group,
        })
    }

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
        Self::from_placed(
            terms,
            panel,
            &terms.place_labels(types.iter().map(String::as_str)),
        )
    }

    /// [`Self::from_ontology`] on the panel's labels already placed
    /// ([`ClTerms::place_labels`]). A label under a term is grouped as its
    /// term is, without being recorded as matched; one kept to a term's
    /// lineage joins a group in it by its markers, else stands alone on
    /// that term.
    #[must_use]
    pub fn from_placed(
        terms: &ClTerms,
        panel: &[(String, String)],
        placed: &Placed,
    ) -> Option<Self> {
        let on = placed.on_terms()?;
        let (mut groups, lone) = class_groups(terms, &on);
        // A lone type the user said to join goes with the group its sharers
        // took, under the class they all share. A second joining the same
        // group goes where the first took it, under the broader of the two
        // classes when one is in the other; else it stays alone, the two
        // sharing no class to group under.
        let mut moved: BTreeMap<String, String> = BTreeMap::new();
        for l in &lone {
            let (Some(true), Some(target)) = (terms.joined(&l.label), &l.target) else {
                continue;
            };
            let (from, to) = match moved.get(target) {
                None => (target.clone(), l.class.clone()),
                Some(to) if terms.ancestors_or_self(to).contains(&l.class) => {
                    (to.clone(), l.class.clone())
                }
                Some(to) if terms.ancestors_or_self(&l.class).contains(to) => {
                    (l.class.clone(), to.clone())
                }
                Some(to) => {
                    log::warn!(
                        "{}: joins a group already under {to}, a class it is not in; left alone",
                        l.label
                    );
                    continue;
                }
            };
            if let Some(members) = groups.remove(&from) {
                groups.entry(to.clone()).or_default().extend(members);
            }
            moved.insert(target.clone(), to);
        }
        let mut groups: Vec<TypeGroup> = groups
            .into_iter()
            .map(|(id, members)| {
                let id = group_term(&on, id, &members);
                TypeGroup {
                    name: label_key(terms.name(&id).unwrap_or(&id)),
                    cl_id: Some(id),
                    members,
                }
            })
            .collect();
        // A group in a kept label's lineage: its term above or below the one
        // that keeps it.
        let in_lineage = |label: &str, g: &TypeGroup| {
            placed.kept.get(label).is_none_or(|kept| {
                g.cl_id.as_deref().is_some_and(|id| {
                    terms.ancestors_or_self(kept).contains(id)
                        || terms.ancestors_or_self(id).contains(kept)
                })
            })
        };
        let kept = |label: &str| placed.kept.get(label).cloned();
        attach_by_sharing(&mut groups, placed.loose(), panel, in_lineage, kept);
        Some(TypeTree {
            source: TreeSource::CellOntology,
            release: terms.release.clone(),
            label_cl: placed.mapped.clone(),
            groups,
        })
    }

    /// The panel's types left alone whose sharers all took one group
    /// ([`JoinQuestion`]), each with the user's answer when there is one.
    #[must_use]
    pub fn left_alone(terms: &ClTerms, placed: &Placed) -> Vec<(JoinQuestion, Option<bool>)> {
        let Some(on) = placed.on_terms() else {
            return Vec::new();
        };
        let (groups, lone) = class_groups(terms, &on);
        (lone.into_iter())
            .filter_map(|l| {
                let target = l.target?;
                // As the tree names it.
                let named = match groups.get(&target) {
                    Some(members) => group_term(&on, target, members),
                    None => target,
                };
                let answer = terms.joined(&l.label);
                Some((
                    JoinQuestion {
                        label: l.label,
                        group: label_key(terms.name(&named).unwrap_or(&named)),
                        class: l.class,
                    },
                    answer,
                ))
            })
            .collect()
    }

    /// The panel's types left alone that the user has not answered for.
    #[must_use]
    pub fn join_questions(terms: &ClTerms, panel: &[(String, String)]) -> Vec<JoinQuestion> {
        let types = panel_types(panel);
        let placed = terms.place_labels(types.iter().map(String::as_str));
        (Self::left_alone(terms, &placed).into_iter())
            .filter(|(_, a)| a.is_none())
            .map(|(q, _)| q)
            .collect()
    }

    /// Groups of types that share marker genes: the connected components of
    /// "has a marker in common". When everything connects, the component is
    /// split by average linkage on marker overlap at its widest gap.
    #[must_use]
    pub fn from_markers(panel: &[(String, String)]) -> Self {
        let types = panel_types(panel);
        let sets = marker_sets(panel);
        let n = types.len();
        // Types share a group through a marker some type lacks; a marker
        // every type holds (IDF 0) tells none apart. Each part is then split
        // where its types stop sharing specific markers.
        let idf = marker_idf(&sets);
        let mut uf = UnionFind::new(n);
        for i in 0..n {
            for j in (i + 1)..n {
                let (a, b) = (&sets[&types[i]], &sets[&types[j]]);
                if a.intersection(b).any(|g| idf[g.as_str()] > 0.0) {
                    uf.union(i, j);
                }
            }
        }
        let mut comps: BTreeMap<usize, Vec<String>> = BTreeMap::new();
        for (i, t) in types.iter().enumerate() {
            comps.entry(uf.find(i)).or_default().push(t.clone());
        }
        let members: Vec<Vec<String>> = comps
            .into_values()
            .flat_map(|c| {
                if c.len() >= 3 {
                    split_by_linkage(&c, &sets, &idf)
                } else {
                    vec![c]
                }
            })
            .collect();
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

/// The term a group of `members` under class `id` is named for: the class,
/// but for a group of one, which a class every sharer left for a nearer one
/// leaves; its type is named for itself, not for the broader class.
fn group_term(on: &BTreeMap<&str, &str>, id: String, members: &[String]) -> String {
    match members {
        [only] => on[only.as_str()].to_string(),
        _ => id,
    }
}

/// A type that shares an analysis class with other panel types, but stands
/// alone because they all went under a nearer class: the class it shares,
/// and the group they took (by its class) when they all took one.
struct Lone {
    label: String,
    class: String,
    target: Option<String>,
}

/// The panel's types on terms (`on`), each under its nearest ancestor (or
/// itself) among the ontology's analysis classes that it shares with another
/// panel type, by class; and the types that left alone ([`Lone`]).
fn class_groups(
    terms: &ClTerms,
    on: &BTreeMap<&str, &str>,
) -> (BTreeMap<String, Vec<String>>, Vec<Lone>) {
    let ancestry: BTreeMap<&str, BTreeSet<String>> = (on.iter())
        .map(|(label, id)| (*label, terms.ancestors_or_self(id)))
        .collect();
    // Per analysis class: how many panel types sit under it, and its depth.
    let mut sharers: HashMap<&str, usize> = HashMap::new();
    for a in ancestry.values() {
        for t in a.iter().filter(|t| terms.classes.contains(*t)) {
            *sharers.entry(t.as_str()).or_default() += 1;
        }
    }
    let depth: HashMap<&str, usize> = sharers.keys().map(|t| (*t, terms.depth(t))).collect();
    let mut groups: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut node_of: BTreeMap<&str, String> = BTreeMap::new();
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
            .unwrap_or_else(|| on[*label].to_string());
        node_of.insert(label, node.clone());
        groups.entry(node).or_default().push((*label).to_string());
    }
    let lone = (groups.iter())
        .filter(|(class, members)| members.len() == 1 && on[members[0].as_str()] != class.as_str())
        .map(|(class, members)| {
            let label = &members[0];
            let theirs: BTreeSet<&String> = (ancestry.iter())
                .filter(|(l, a)| **l != label.as_str() && a.contains(class))
                .map(|(l, _)| &node_of[*l])
                .collect();
            Lone {
                label: label.clone(),
                class: class.clone(),
                target: (theirs.len() == 1).then(|| (*theirs.iter().next().unwrap()).clone()),
            }
        })
        .collect();
    (groups, lone)
}

/// A question for the user: does `label`, left alone, join `group` under
/// the analysis class `class` they share, or stand alone? Answered in
/// `cl_groups.tsv` ([`ClTerms::add_joins`]).
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct JoinQuestion {
    pub label: String,
    pub class: String,
    pub group: String,
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

/// Each marker's IDF over the panel's types, `ln(types / types holding it)`,
/// the weighting `enrichment::markers::apply_idf_weights` gives the scoring
/// pass: a marker most types hold says little about which a type is.
fn marker_idf(sets: &BTreeMap<String, BTreeSet<String>>) -> HashMap<&str, f64> {
    let mut holders: HashMap<&str, usize> = HashMap::new();
    for g in sets.values().flatten() {
        *holders.entry(g.as_str()).or_default() += 1;
    }
    let n = sets.len() as f64;
    holders
        .into_iter()
        .map(|(g, k)| (g, (n / k as f64).ln()))
        .collect()
}

/// Jaccard with each marker weighed by its IDF.
fn weighted_jaccard(a: &BTreeSet<String>, b: &BTreeSet<String>, idf: &HashMap<&str, f64>) -> f64 {
    let w = |g: &String| idf[g.as_str()];
    let shared: f64 = a.intersection(b).map(w).sum();
    let either: f64 = a.union(b).map(w).sum();
    if either > 0.0 {
        shared / either
    } else {
        0.0
    }
}

/// Average-linkage agglomeration on `1 -` IDF-weighted Jaccard, cut where
/// successive merge heights jump the most; with no jump at all, the types
/// stay one group.
fn split_by_linkage(
    types: &[String],
    sets: &BTreeMap<String, BTreeSet<String>>,
    idf: &HashMap<&str, f64>,
) -> Vec<Vec<String>> {
    let mut clusters: Vec<Vec<usize>> = (0..types.len()).map(|i| vec![i]).collect();
    let dist = |a: &[usize], b: &[usize]| {
        let total: f64 = a
            .iter()
            .flat_map(|&i| b.iter().map(move |&j| (i, j)))
            .map(|(i, j)| 1.0 - weighted_jaccard(&sets[&types[i]], &sets[&types[j]], idf))
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
    // No jump at all: nothing in the component sets types apart.
    let jumps = history.windows(2).map(|w| w[1].0 - w[0].0);
    if jumps.clone().fold(0.0, f64::max) <= 0.0 {
        return vec![types.to_vec()];
    }
    let cut = jumps
        .enumerate()
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
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

/// Put each unmapped type in the group, of those `may_join` allows it,
/// with the member that holds most of its markers' IDF weight
/// ([`marker_idf`]; the first such group on a tie), or in a group of its
/// own, on the term `term_of` gives it.
fn attach_by_sharing(
    groups: &mut Vec<TypeGroup>,
    unmapped: Vec<String>,
    panel: &[(String, String)],
    may_join: impl Fn(&str, &TypeGroup) -> bool,
    term_of: impl Fn(&str) -> Option<String>,
) {
    let sets = marker_sets(panel);
    let idf = marker_idf(&sets);
    for t in unmapped {
        let own = &sets[&t];
        let best = groups
            .iter()
            .enumerate()
            .filter(|(_, grp)| may_join(&t, grp))
            .map(|(i, grp)| {
                let held = grp
                    .members
                    .iter()
                    .filter_map(|m| sets.get(m))
                    .map(|m| own.intersection(m).map(|g| idf[g.as_str()]).sum::<f64>())
                    .fold(0.0, f64::max);
                (i, held)
            })
            .filter(|(_, held)| *held > 0.0)
            .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0)));
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
                cl_id: term_of(&t),
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

#[cfg(test)]
#[path = "tests/match_corpus.rs"]
mod match_corpus;
