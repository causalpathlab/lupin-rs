//! A coarse level over a marker panel's cell types, for the first round of
//! annotation: each fine type belongs to one broad group, and a cluster's
//! first call is the group its evidence adds up to.
//!
//! The groups come from the Cell Ontology when there is one: the panel's
//! types are placed on CL terms under the matching rules
//! ([`ClTerms::place_labels`]), and each goes under its nearest ancestor
//! among the ontology's analysis classes (its curated `*_upper_slim` subsets
//! and `cellxgene_subset`) that it shares with another type's term, short of
//! one that is a panel type's own term over a kind of it; a group whose
//! class is a kind of another group's then joins that one, and groups below
//! one curated lineage (a domain `*_upper_slim` term) are one when that adds
//! no level: when no other group lies below the lineage, or when the
//! lineage is already one of the merging groups' class (a panel type's own
//! term). A type
//! the ontology has only a kind of term for (`Early X` for `X`) groups as
//! that term; one that matches no term joins the group sharing its most
//! specific markers, keeping to the lineage its words hint at. A fine panel
//! collapses to lineage classes; a broad one stays about as it is. Without
//! an ontology, or when too few types match, types that share specific
//! markers are grouped together.
//! Nothing here names a cell type; every group is derived from the panel and
//! the ontology at run time.
//!
//! Pure: [`crate::manifest`] finds the ontology and the evidence.

use crate::annotate::cl_rules::{part_hyphens, Aliases, MatchRules, SIGN_WORDS};
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
    /// The names and counted synonyms, as written.
    names: NameIndex,
    /// The same without a trailing qualifier (see `QualifierRules`), looked
    /// at only when the names as written say nothing about a label.
    unqualified: NameIndex,
    name_of: HashMap<String, String>,
    parents: HashMap<String, Vec<String>>,
    /// `is_a` children, the reverse of `parents`.
    children: HashMap<String, Vec<String>>,
    /// `develops_from` (RO:0002202): the terms a term's cells develop from.
    develops_from: HashMap<String, Vec<String>>,
    /// Terms in the analysis subsets the rules name (`ClassRules`).
    classes: BTreeSet<String>,
    /// Terms in the lineage subsets the rules name (`ClassRules`).
    lineages: BTreeSet<String>,
    /// The words the ontology uses for kinds of a term (see `PartRules`),
    /// found on first use.
    modifiers: std::sync::OnceLock<BTreeSet<String>>,
    /// The word a panel leaves off its labels, once [`Self::set_panel`]
    /// has read one from the panel: labels are then completed by it, however
    /// few are mapped. Without one, it is told from the labels mapped.
    panel_context: Option<String>,
    /// The file's `data-version`, when it states one.
    pub release: Option<String>,
}

/// Term names keyed every way a label is compared with them.
#[derive(Default)]
struct NameIndex {
    /// Normalised name or counted synonym → term ids.
    by_name: HashMap<String, BTreeSet<String>>,
    /// The same names by [`word_bag`]: for labels whose words come in
    /// another order (`CT1 cells memory` for `memory CT1 cell`) when the
    /// rules allow it, and for the panel's word in any case.
    by_words: HashMap<String, BTreeSet<String>>,
    /// Each term's own name (no synonyms) by [`ClTerms::name_key`]: what an
    /// implied trailing word may complete.
    by_term_name: HashMap<String, BTreeSet<String>>,
    /// Names and counted synonyms without spaces ([`compact`]), when the rules
    /// compare so ([`MatchRules::compact`]); and the own names alone.
    by_compact: HashMap<String, BTreeSet<String>>,
    by_compact_name: HashMap<String, BTreeSet<String>>,
    /// Each own name's [`marker_words`], for labels a name spells out among
    /// other words ([`crate::annotate::cl_rules::ContainedRules`]).
    phrases: Vec<(BTreeSet<String>, String)>,
}

impl NameIndex {
    /// Index `name` of term `id`; `own` when it is the term's own name.
    fn add(&mut self, rules: &MatchRules, id: &str, name: &str, own: bool) {
        let add = |map: &mut HashMap<String, BTreeSet<String>>, key: String| {
            map.entry(key).or_default().insert(id.to_string());
        };
        let n = rules.normalise(name);
        let single = rules.singular(&n);
        add(&mut self.by_words, word_bag(&single));
        if rules.compact {
            add(&mut self.by_compact, compact(&single));
            if own {
                add(&mut self.by_compact_name, compact(&single));
            }
        }
        // Own names only: a synonym is already a looser link (see
        // `ImpliedRules`).
        if own {
            add(&mut self.by_term_name, name_key(rules, &n));
            if rules.contained.min_words.is_some() {
                self.phrases.push((marker_words(&single), id.to_string()));
            }
        }
        add(&mut self.by_name, n);
    }
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
        for (i, name) in self.names.iter().enumerate() {
            terms.names.add(&terms.rules, &self.id, name, i == 0);
            if let Some(stem) = unqualified(&terms.rules, name) {
                terms.unqualified.add(&terms.rules, &self.id, stem, i == 0);
            }
        }
        if let Some(first) = self.names.first() {
            terms.name_of.insert(self.id.clone(), first.clone());
        }
        if terms.rules.classes.is_lineage(&self.subsets) {
            terms.lineages.insert(self.id.clone());
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

/// A normalised name's words: spaces part them, and so does a hyphen
/// between two words of letters alone or before a marker's written-out sign
/// (`CD4-positive`); any other hyphen is kept ([`part_hyphens`]).
fn words(normalised: &str) -> Vec<String> {
    (part_hyphens(normalised).split(' '))
        .filter(|w| !w.is_empty())
        .map(String::from)
        .collect()
}

/// Whether a normalised name has a word at all.
fn has_words(normalised: &str) -> bool {
    normalised.split(' ').any(|w| !w.is_empty())
}

/// A normalised name's [`words`] in sorted order.
fn word_bag(normalised: &str) -> String {
    let mut words = words(normalised);
    words.sort_unstable();
    words.join(" ")
}

/// A name as [`NameIndex::by_term_name`] keys it: normalised, made
/// singular, and its words sorted when the rules allow any word order.
fn name_key(rules: &MatchRules, name: &str) -> String {
    let single = rules.singular(&rules.normalise(name));
    if rules.any_word_order {
        word_bag(&single)
    } else {
        single
    }
}

/// The one id of `ids`, when there is exactly one.
fn only(ids: &BTreeSet<String>) -> Option<String> {
    (ids.len() == 1).then(|| ids.iter().next().cloned())?
}

/// A normalised name's [`words`] run together: no spaces, nor a hyphen
/// between words; a sign's hyphen is kept.
fn compact(normalised: &str) -> String {
    words(normalised).concat()
}

/// A normalised name's words as a label may give them: parted by spaces
/// only, so `CD4-negative` stays itself, and a marker's `-positive` or `+`
/// dropped, so `CD4-positive`, `CD4+` and `CD4` are one word.
fn marker_words(normalised: &str) -> BTreeSet<String> {
    normalised
        .split(' ')
        .map(|w| {
            w.strip_suffix("-positive")
                .or_else(|| w.strip_suffix('+'))
                .unwrap_or(w)
        })
        .filter(|w| !w.is_empty())
        .map(String::from)
        .collect()
}

/// `name` without the trailing qualifier the rules drop
/// ([`crate::annotate::cl_rules::QualifierRules`]): a parenthesis, or a
/// taxon the rules name after a comma.
fn unqualified<'a>(rules: &MatchRules, name: &'a str) -> Option<&'a str> {
    let name = name.trim_end();
    let stem = if rules.qualifiers.parenthetical && name.ends_with(')') {
        name.rfind(" (").map(|i| &name[..i])
    } else {
        None
    };
    stem.or_else(|| {
        let (stem, word) = name.rsplit_once(", ")?;
        let taxa = &rules.qualifiers.comma_words;
        taxa.iter()
            .any(|t| t.eq_ignore_ascii_case(word))
            .then_some(stem)
    })
    .map(str::trim_end)
    .filter(|s| !s.is_empty())
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

    /// The words the ontology uses for kinds of a term (see `PartRules`):
    /// over every pair of own names `X` and `X` with one word more, those
    /// whose added word makes a kind of `X` in enough pairs and in enough of
    /// its pairs.
    fn modifier_words(&self) -> BTreeSet<String> {
        let Some(min) = self.rules.part.min_pairs else {
            return BTreeSet::new();
        };
        let mut by_bag: HashMap<Vec<String>, Vec<&str>> = HashMap::new();
        for (id, name) in &self.name_of {
            let single = self.rules.singular(&self.rules.normalise(name));
            let bag: Vec<String> = (word_bag(&single).split(' '))
                .filter(|w| !w.is_empty())
                .map(String::from)
                .collect();
            by_bag.entry(bag).or_default().push(id);
        }
        // Per added word: pairs where it makes a kind of the shorter name, and not.
        let mut pairs: HashMap<&str, (usize, usize)> = HashMap::new();
        for (bag, ids) in &by_bag {
            let above: Vec<BTreeSet<String>> =
                ids.iter().map(|id| self.ancestors_or_self(id)).collect();
            for (i, word) in bag.iter().enumerate() {
                let mut base = bag.clone();
                base.remove(i);
                let Some(bases) = by_bag.get(&base) else {
                    continue;
                };
                let e = pairs.entry(word).or_default();
                for a in &above {
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
        pairs
            .into_iter()
            .filter(|(_, (kind, other))| {
                *kind >= min && *kind as f64 >= self.rules.part.min_share * (kind + other) as f64
            })
            .map(|(w, _)| w.to_string())
            .collect()
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
    /// (a hyphen between words as a space, a marker's sign kept, see
    /// [`part_hyphens`]), else with an implied
    /// trailing word completing a term's name, else as the rules allow
    /// beyond: without spaces, by words a term's name spells out among others, and with a
    /// word the panel leaves off its labels. Where a text names several
    /// terms, a term's own name wins over another's synonym, and a term
    /// every other is a kind of wins over its specialisations; a label still
    /// matching several, or none, stays unmapped rather than guessed.
    ///
    /// The word the panel leaves off is the one [`Self::set_panel`] read
    /// from the panel, so a few of its labels, or one, read as the panel's
    /// do; without one, it is told from `labels`.
    #[must_use]
    pub fn map_labels<'a>(&self, labels: impl IntoIterator<Item = &'a str>) -> Mapped {
        let labels: Vec<&str> = labels.into_iter().collect();
        let read = self.read(&labels);
        let context = self.context_word(&read);
        self.complete_read(read, context.as_deref())
    }

    /// Each label's term by alias or name ([`Self::match_plain`]), and the
    /// labels no name matches: one pass over `labels`, for the mapping and
    /// the panel's word alike.
    fn read(&self, labels: &[&str]) -> Read {
        let mut read = Read::default();
        for label in labels {
            let name = self.match_name(label);
            if name.is_none() {
                read.unnamed.push((*label).to_string());
            }
            match self.alias_of(label).or(name) {
                Some(id) => {
                    read.mapped.insert((*label).to_string(), id);
                }
                None => read.unmapped.push((*label).to_string()),
            }
        }
        read
    }

    /// `read` with the labels it left completed by the panel's word
    /// `context` ([`Self::complete`]).
    fn complete_read(&self, read: Read, context: Option<&str>) -> Mapped {
        let Read {
            mut mapped,
            mut unmapped,
            ..
        } = read;
        if let Some(w) = context {
            unmapped.retain(|l| match self.complete(l, w) {
                Some(id) => {
                    mapped.insert(l.clone(), id);
                    false
                }
                None => true,
            });
        }
        (mapped, unmapped)
    }

    /// The word the panel leaves off: the one [`Self::set_panel`] read,
    /// else the one told from the labels `read` saw.
    fn context_word(&self, read: &Read) -> Option<String> {
        (self.panel_context.clone()).or_else(|| self.context_of(&read.unnamed).map(|(w, _)| w))
    }

    /// Read the word marker panel `panel` leaves off its labels, for
    /// [`Self::map_labels`] to complete any of them by; a panel that leaves
    /// none (or an empty one) leaves the word to be told from the labels
    /// given.
    pub fn set_panel(&mut self, panel: &[(String, String)]) {
        let types = panel_types(panel);
        let types: Vec<&str> = types.iter().map(String::as_str).collect();
        let read = self.read(&types);
        self.panel_context = self.context_of(&read.unnamed).map(|(w, _)| w);
    }

    /// [`Self::map_labels`], then each label left by a run of its words
    /// (see `PartRules`): under the term that run matches, or with that term
    /// only as a hint of its lineage.
    #[must_use]
    pub fn place_labels<'a>(&self, labels: impl IntoIterator<Item = &'a str>) -> Placed {
        let labels: Vec<&str> = labels.into_iter().collect();
        let read = self.read(&labels);
        let context = self.context_word(&read);
        let (on, unmapped) = self.complete_read(read, context.as_deref());
        let mut placed = Placed {
            on,
            ..Placed::default()
        };
        for label in unmapped {
            match self.by_part(&label, context.as_deref()) {
                Some((id, true)) => {
                    placed.under.insert(label, id);
                }
                Some((id, false)) => {
                    placed.hint.insert(label, id);
                }
                None => placed.none.push(label),
            }
        }
        placed
    }

    /// The term the longest run of `label`'s words matches, when one run
    /// of that length matches one term, and whether every word left off is
    /// a modifier. At each length, the runs a name or alias matches
    /// ([`Self::match_plain`]) come first; only when none does are runs
    /// completed by the panel's word `context`.
    fn by_part(&self, label: &str, context: Option<&str>) -> Option<(String, bool)> {
        self.rules.part.min_pairs?;
        let w = self.label_words(label);
        for len in (1..w.len()).rev() {
            let runs = || (0..=w.len() - len).map(|start| (start, w[start..start + len].join(" ")));
            let found = |term: &dyn Fn(&str) -> Option<String>| {
                let mut found: BTreeMap<String, bool> = BTreeMap::new();
                for (start, run) in runs() {
                    if let Some(id) = term(&run) {
                        let left = w[..start].iter().chain(&w[start + len..]);
                        let modifiers = self.modifiers.get_or_init(|| self.modifier_words());
                        let kind = left.clone().all(|d| modifiers.contains(d));
                        *found.entry(id).or_default() |= kind;
                    }
                }
                found
            };
            let mut at = found(&|run| self.match_plain(run));
            if at.is_empty() {
                if let Some(c) = context {
                    at = found(&|run| self.complete(run, c));
                }
            }
            match at.len() {
                0 => {}
                1 => return at.into_iter().next(),
                _ => return None,
            }
        }
        None
    }

    /// One label's term by every tier but the panel's word: a curated alias
    /// first, when its term is in this ontology, else the names.
    fn match_plain(&self, label: &str) -> Option<String> {
        self.alias_of(label).or_else(|| self.match_name(label))
    }

    /// The term `label` is aliased to, when this ontology has it.
    fn alias_of(&self, label: &str) -> Option<String> {
        (self.aliases.get(label))
            .filter(|id| self.has(id))
            .map(str::to_string)
    }

    /// One label's term by the names alone, every tier but the aliases and
    /// the panel's context.
    fn match_name(&self, label: &str) -> Option<String> {
        let exact = self.rules.normalise(label);
        let single = self.rules.singular(&exact);
        // A label with no word names nothing.
        if !has_words(&single) {
            return None;
        }
        self.match_in(&self.names, &exact, &single)
            .or_else(|| self.match_in(&self.unqualified, &exact, &single))
            .flatten()
    }

    /// A label (`exact`, normalised; `single`, made singular) on `ix`, tier
    /// by tier: the first tier whose text `ix` has decides, matched or not
    /// (`Some(None)` for a text naming several terms); `None` when none has.
    fn match_in(&self, ix: &NameIndex, exact: &str, single: &str) -> Option<Option<String>> {
        // As written, else made singular, else in any word order. Where a
        // text names several terms, the one it is the own name of wins over
        // those it is a synonym of.
        let one = |ids: &BTreeSet<String>, text: &str| -> Option<String> {
            only(ids)
                .or_else(|| {
                    let named = ix.by_term_name.get(&self.name_key(text))?;
                    only(&ids.intersection(named).cloned().collect())
                })
                .or_else(|| self.most_general(ids))
        };
        if let Some(ids) = ix.by_name.get(exact) {
            return Some(one(ids, exact));
        }
        if let Some(ids) = ix.by_name.get(single) {
            return Some(one(ids, single));
        }
        if let Some(ids) = (self.rules.any_word_order)
            .then(|| ix.by_words.get(&word_bag(single)))
            .flatten()
        {
            return Some(one(ids, single));
        }
        // An implied trailing word completing a term's own name (see
        // `ImpliedRules`), the first that names one; then the same compared
        // without separators.
        let implied = |index: &HashMap<String, BTreeSet<String>>, key: &dyn Fn(&str) -> String| {
            let named: Vec<&BTreeSet<String>> = (self.rules.implied.suffixes.iter())
                .filter_map(|w| index.get(&key(&format!("{single} {w}"))))
                .collect();
            (!named.is_empty()).then(|| named.into_iter().find_map(|ids| self.most_general(ids)))
        };
        if let Some(found) = implied(&ix.by_term_name, &|t| self.name_key(t)) {
            return Some(found);
        }
        if let Some(ids) = ix.by_compact.get(&compact(single)) {
            return Some(self.most_general(ids));
        }
        if let Some(found) = implied(&ix.by_compact_name, &compact) {
            return Some(found);
        }
        self.contained(ix, single)
    }

    /// The most general term whose own name spells out `single`'s words
    /// among others (see `ContainedRules`), else the one whose name adds the
    /// fewest words; `None` when no name spells them out.
    fn contained(&self, ix: &NameIndex, single: &str) -> Option<Option<String>> {
        let min = self.rules.contained.min_words?;
        let want = marker_words(single);
        let counted = want
            .iter()
            .filter(|w| !self.rules.implied.suffixes.contains(w))
            .count();
        if counted < min {
            return None;
        }
        let found: Vec<(usize, &String)> = (ix.phrases.iter())
            .filter(|(ws, _)| want.is_subset(ws))
            .map(|(ws, id)| (ws.len(), id))
            .collect();
        let ids: BTreeSet<String> = found.iter().map(|(_, id)| (*id).clone()).collect();
        // Else the name that adds the fewest words, when only one does.
        let tightest = || {
            let fewest = found.iter().map(|(n, _)| *n).min()?;
            only(
                &(found.iter())
                    .filter(|(n, _)| *n == fewest)
                    .map(|(_, id)| (*id).clone())
                    .collect(),
            )
        };
        (!ids.is_empty()).then(|| self.most_general(&ids).or_else(tightest))
    }

    /// Per word that, with at most an implied one, turns labels of `labels`
    /// into a term's name or counted synonym: those labels, with their
    /// terms (see `ContextRules`). A scan of every name finds the words;
    /// [`Self::completed_by`] completes by each.
    fn completions(&self, labels: &[String]) -> BTreeMap<&str, BTreeMap<String, String>> {
        let implied = &self.rules.implied.suffixes;
        let mut found: BTreeSet<&str> = BTreeSet::new();
        for label in labels {
            let own = self.label_words(label);
            if own.is_empty() {
                continue;
            }
            for bag in self.names.by_words.keys() {
                let mut rest: Vec<&str> = bag.split(' ').collect();
                let all_there = own.iter().all(|w| {
                    rest.iter()
                        .position(|r| r == w)
                        .map(|i| rest.remove(i))
                        .is_some()
                });
                if !all_there || rest.len() > 2 {
                    continue;
                }
                if rest.len() == 2 {
                    let Some(i) = rest.iter().position(|r| implied.iter().any(|w| w == r)) else {
                        continue;
                    };
                    rest.remove(i);
                }
                if let [w] = rest[..] {
                    found.insert(w);
                }
            }
        }
        (found.into_iter())
            .map(|w| (w, self.completed_by(labels, w)))
            .filter(|(_, completed)| !completed.is_empty())
            .collect()
    }

    /// The labels of `labels` word `w` completes ([`Self::complete`]), with
    /// their terms.
    fn completed_by(&self, labels: &[String], w: &str) -> BTreeMap<String, String> {
        (labels.iter())
            .filter_map(|l| Some((l.clone(), self.complete(l, w)?)))
            .collect()
    }

    /// `label`'s words with word `w`, or with `w` and an implied word,
    /// looked up among the names: the most general term they name. An
    /// implied word or a marker's sign ([`SIGN_WORDS`]) completes nothing: a
    /// label leaving off `negative` is no negative marker's.
    fn complete(&self, label: &str, w: &str) -> Option<String> {
        let implied = &self.rules.implied.suffixes;
        if implied.iter().any(|s| s == w) || SIGN_WORDS.contains(&w) {
            return None;
        }
        let own = self.label_words(label);
        if own.is_empty() {
            return None;
        }
        let bag = |more: &[&str]| {
            let mut all: Vec<&str> = own.iter().map(String::as_str).chain([w]).collect();
            all.extend(more);
            all.sort_unstable();
            all.join(" ")
        };
        let ids: BTreeSet<String> = std::iter::once(bag(&[]))
            .chain(implied.iter().map(|s| bag(&[s.as_str()])))
            .filter_map(|key| self.names.by_words.get(&key))
            .flatten()
            .cloned()
            .collect();
        self.most_general(&ids)
    }

    /// `label`'s words, as names' are compared: normalised and made singular.
    fn label_words(&self, label: &str) -> Vec<String> {
        words(&self.rules.singular(&self.rules.normalise(label)))
    }

    /// The panel's context word among the completions of `labels` (the ones
    /// no name matches, [`Self::completions`]): the one that completes the
    /// most labels, when it completes enough of them and more than any other
    /// word, with the labels it completes and their terms. An aliased label
    /// counts for a word only when the word completes it to its alias's term.
    /// None when the rules leave context off.
    fn context_of(&self, labels: &[String]) -> Option<(String, BTreeMap<String, String>)> {
        let min = self.rules.context.min_labels?;
        let mut ranked: Vec<(&str, BTreeMap<String, String>)> = (self.completions(labels))
            .into_iter()
            .map(|(w, completed)| {
                let agree = |(l, id): &(String, String)| self.alias_of(l).is_none_or(|a| a == *id);
                (
                    w,
                    completed
                        .into_iter()
                        .filter(agree)
                        .collect::<BTreeMap<_, _>>(),
                )
            })
            .filter(|(_, completed)| !completed.is_empty())
            .collect();
        ranked.sort_by_key(|(_, m)| std::cmp::Reverse(m.len()));
        let mut ranked = ranked.into_iter();
        let (w, completed) = ranked.next()?;
        let n = completed.len();
        (n >= min.max(1) && ranked.next().is_none_or(|r| r.1.len() < n))
            .then(|| (w.to_string(), completed))
    }

    /// The one of `ids` every other is a kind of (`is_a`), if there is one.
    fn most_general(&self, ids: &BTreeSet<String>) -> Option<String> {
        if let Some(id) = only(ids) {
            return Some(id);
        }
        let above: Vec<BTreeSet<String>> =
            ids.iter().map(|id| self.ancestors_or_self(id)).collect();
        ids.iter()
            .find(|id| above.iter().all(|a| a.contains(*id)))
            .cloned()
    }

    /// A name as [`NameIndex::by_term_name`] keys it ([`name_key`]).
    fn name_key(&self, name: &str) -> String {
        name_key(&self.rules, name)
    }

    /// Every counted name (as names are compared) with the terms it names.
    pub fn named_texts(&self) -> impl Iterator<Item = (&str, &BTreeSet<String>)> {
        self.names
            .by_name
            .iter()
            .map(|(name, ids)| (name.as_str(), ids))
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
        // Ranked: a whole-name hit (any name), then by where the term's own
        // name has the query (at its start, inside, not at all: a synonym's
        // hit), then shorter names.
        let mut hits: BTreeMap<String, (bool, u8, usize)> = BTreeMap::new();
        for (name, ids) in &self.names.by_name {
            if !name.contains(&q) {
                continue;
            }
            for id in ids {
                let e = hits.entry(id.clone()).or_insert_with(|| {
                    let own = self.name_of.get(id);
                    let at = match own.map(|n| self.rules.normalise(n)) {
                        Some(n) if n.starts_with(&q) => 0,
                        Some(n) if n.contains(&q) => 1,
                        _ => 2,
                    };
                    (false, at, own.map_or(usize::MAX, String::len))
                });
                e.0 |= *name == q;
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

/// Labels with the terms they map to, and the labels left unmapped.
pub type Mapped = (BTreeMap<String, String>, Vec<String>);

/// What one pass over labels found ([`ClTerms::read`]): each label's term by
/// alias or name, the labels left, and the labels no name matches (aliased
/// or not), which the panel's word is told from.
#[derive(Default)]
struct Read {
    mapped: BTreeMap<String, String>,
    unmapped: Vec<String>,
    unnamed: Vec<String>,
}

/// Where [`ClTerms::place_labels`] put a panel's labels.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Placed {
    /// Label → the term it matched.
    pub on: BTreeMap<String, String>,
    /// Label → the term it is a kind of.
    pub under: BTreeMap<String, String>,
    /// Label → a term of its lineage.
    pub hint: BTreeMap<String, String>,
    /// Labels with nothing on the ontology.
    pub none: Vec<String>,
}

impl Placed {
    /// The term each label sits on or under.
    #[must_use]
    pub fn anchored(&self) -> BTreeMap<&str, &str> {
        (self.on.iter())
            .chain(&self.under)
            .map(|(l, id)| (l.as_str(), id.as_str()))
            .collect()
    }

    /// The labels on no term, hinted or not, hinted first by name.
    #[must_use]
    pub fn loose(&self) -> Vec<String> {
        self.hint
            .keys()
            .cloned()
            .chain(self.none.iter().cloned())
            .collect()
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
    /// Panel labels placed under a CL term, as kinds of it
    /// ([`ClTerms::place_labels`]).
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub label_under: BTreeMap<String, String>,
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

/// What placed a label: its own term, the term it is a kind of, or its
/// group's.
#[derive(Serialize, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum Via {
    Label,
    Under,
    Group,
}

impl TypeTree {
    /// `label`'s place on the Cell Ontology: the term it matched, else the
    /// one it is placed under, else the term of the group it is (a coarse
    /// label) or belongs to.
    #[must_use]
    pub fn placement(&self, label: &str) -> Option<TreePlacement> {
        let key = label_key(label);
        let own = [
            (&self.label_cl, Via::Label),
            (&self.label_under, Via::Under),
        ];
        for (terms, via) in own {
            if let Some((name, id)) = terms.iter().find(|(l, _)| label_key(l) == key) {
                return Some(TreePlacement {
                    id: id.clone(),
                    name: name.clone(),
                    via,
                });
            }
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
        let placed = terms.place_labels(types.iter().map(String::as_str));
        Self::from_placed(terms, panel, &placed)
    }

    /// [`Self::from_ontology`] with the panel's types already placed
    /// ([`ClTerms::place_labels`]).
    #[must_use]
    pub fn from_placed(
        terms: &ClTerms,
        panel: &[(String, String)],
        placed: &Placed,
    ) -> Option<Self> {
        if placed.on.len() < 2 {
            return None;
        }
        // A type under a term groups as that term would: the groups are
        // made of the terms types sit on or under, each with its types.
        let mut labels_on: BTreeMap<&str, Vec<String>> = BTreeMap::new();
        for (label, id) in placed.anchored() {
            labels_on.entry(id).or_default().push(label.to_string());
        }
        // Each term's terms above it, and so its depth, walked once.
        let walked: std::cell::RefCell<HashMap<String, std::rc::Rc<BTreeSet<String>>>> =
            Default::default();
        let above = |id: &str| {
            (walked.borrow_mut().entry(id.to_string()))
                .or_insert_with(|| std::rc::Rc::new(terms.ancestors_or_self(id)))
                .clone()
        };
        let depth_of = |id: &str| above(id).len();
        let ancestry: BTreeMap<&str, std::rc::Rc<BTreeSet<String>>> =
            labels_on.keys().map(|id| (*id, above(id))).collect();
        // Each term goes under its nearest ancestor (or itself) among the
        // ontology's analysis classes that it shares with another type's
        // term. A term sharing none stands as its own group.
        // Per analysis class: how many of the terms sit under it.
        let mut sharers: HashMap<&str, usize> = HashMap::new();
        for a in ancestry.values() {
            for t in a.iter().filter(|t| terms.classes.contains(*t)) {
                *sharers.entry(t.as_str()).or_default() += 1;
            }
        }
        let depth: HashMap<&str, usize> = sharers.keys().map(|t| (*t, depth_of(t))).collect();
        // A class that is one panel type's own term groups no other term
        // below it: over a kind of it, it would make that kind look like the
        // type.
        let on: BTreeSet<&str> = placed.on.values().map(String::as_str).collect();
        let groups_by = |t: &str| sharers.get(t).copied().unwrap_or(0) >= 2 && !on.contains(t);
        // Group class → the terms in it.
        let mut groups: BTreeMap<String, Vec<&str>> = BTreeMap::new();
        for (id, a) in &ancestry {
            let node = a
                .iter()
                .filter(|t| groups_by(t))
                .max_by(|x, y| {
                    depth[x.as_str()]
                        .cmp(&depth[y.as_str()])
                        .then_with(|| y.cmp(x))
                })
                .cloned()
                .unwrap_or_else(|| (*id).to_string());
            groups.entry(node).or_default().push(id);
        }
        // One coarse level: a group whose class is a kind of another's joins
        // it, the deepest first, so nested groups end in the outermost.
        let mut shared: Vec<(usize, String)> = (groups.iter())
            .filter(|(_, m)| m.len() >= 2)
            .map(|(id, _)| (depth_of(id), id.clone()))
            .collect();
        shared.sort_by_key(|(d, _)| std::cmp::Reverse(*d));
        for (_, inner) in &shared {
            let over = above(inner);
            let outer = (shared.iter())
                .filter(|(_, o)| o != inner && over.contains(o) && groups.contains_key(o))
                .max_by_key(|(d, _)| *d)
                .map(|(_, o)| o);
            if let Some(outer) = outer {
                let moved = groups.remove(inner).unwrap_or_default();
                groups.entry(outer.clone()).or_default().extend(moved);
            }
        }
        // Groups below one lineage are one group: each group's lineage is the
        // nearest above its class, or above its lone term. The groups stay
        // one level: a lineage over a group it does not take merges nothing,
        // unless it is already one merging group's class (a panel type's own
        // term), which adds no level.
        let class_of = |id: &str, members: &[&str]| match members {
            [m] => (*m).to_string(),
            _ => id.to_string(),
        };
        let classes: BTreeMap<String, String> = (groups.iter())
            .map(|(id, members)| (id.clone(), class_of(id, members)))
            .collect();
        let mut by_lineage: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for (id, class) in &classes {
            let lineage = (above(class).iter())
                .filter(|t| terms.lineages.contains(*t))
                .max_by_key(|t| depth_of(t))
                .cloned();
            if let Some(l) = lineage {
                by_lineage.entry(l).or_default().push(id.clone());
            }
        }
        // Every merging group is taken out before any lineage is put in, as
        // a lineage may be the class of a group merging elsewhere.
        let merged: Vec<(String, Vec<&str>)> = (by_lineage.into_iter())
            .filter(|(lineage, ids)| {
                let is_a_class = ids.iter().any(|id| classes[id] == *lineage);
                ids.len() >= 2
                    && (is_a_class
                        || (classes.iter())
                            .filter(|(id, _)| !ids.contains(id))
                            .all(|(_, class)| !above(class).contains(lineage)))
            })
            .map(|(lineage, ids)| {
                let members = (ids.iter())
                    .flat_map(|id| groups.remove(id).unwrap_or_default())
                    .collect();
                (lineage, members)
            })
            .collect();
        for (lineage, members) in merged {
            groups.entry(lineage).or_default().extend(members);
        }
        let mut groups: Vec<TypeGroup> = groups
            .into_iter()
            .map(|(id, ids)| {
                // A class every sharer left for a nearer one groups nothing:
                // its lone term's types are named for it, not for the
                // broader class.
                let id = match &ids[..] {
                    [one] => (*one).to_string(),
                    _ => id,
                };
                let mut members: Vec<String> = ids
                    .iter()
                    .flat_map(|t| labels_on[t].iter().cloned())
                    .collect();
                members.sort_unstable();
                TypeGroup {
                    name: String::new(),
                    cl_id: Some(id),
                    members,
                }
            })
            .collect();
        // A hinted type joins no group of another lineage: one whose term
        // neither is above its hint nor below it.
        let fits = |g: &TypeGroup, t: &str| {
            let (Some(hint), Some(class)) = (placed.hint.get(t), g.cl_id.as_deref()) else {
                return true;
            };
            above(hint).contains(class) || above(class).contains(hint)
        };
        attach_by_sharing(&mut groups, placed.loose(), panel, fits);
        // Named once every type is in: a lone type under a term is no more
        // that term than it was, but a group it shares is the term's.
        for g in &mut groups {
            if let Some(id) = &g.cl_id {
                g.name = match &g.members[..] {
                    [m] if placed.under.contains_key(m) => label_key(m),
                    _ => label_key(terms.name(id).unwrap_or(id)),
                };
            }
        }
        Some(TypeTree {
            source: TreeSource::CellOntology,
            release: terms.release.clone(),
            label_cl: placed.on.clone(),
            label_under: placed.under.clone(),
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
            label_under: BTreeMap::new(),
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

/// Put each unmapped type in the group, among those it `fits`, with the
/// member that holds most of its markers' IDF weight ([`marker_idf`]; the
/// first such group on a tie), or in a group of its own.
fn attach_by_sharing(
    groups: &mut Vec<TypeGroup>,
    unmapped: Vec<String>,
    panel: &[(String, String)],
    fits: impl Fn(&TypeGroup, &str) -> bool,
) {
    let sets = marker_sets(panel);
    let idf = marker_idf(&sets);
    for t in unmapped {
        let own = &sets[&t];
        let best = groups
            .iter()
            .enumerate()
            .filter(|(_, grp)| fits(grp, &t))
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
