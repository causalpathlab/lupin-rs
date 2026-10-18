//! How panel labels meet the Cell Ontology, as data: the matching rules
//! (`cl_matching.json`) and the label → term aliases a marker panel carries
//! (its sidecar, or a `cl_aliases.tsv`). Nothing here is decided in code: without a rules file
//! matching is literal (a term's name or exact synonym, as written), and every
//! heuristic (plurals, word order, abbreviations) comes from a file. Finding
//! and layering the files is [`crate::manifest::data_files`]'s job.

use crate::annotate::markers::label_key;
use legume_numeric::matrix::common_io::TwoFieldRows;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;

/// Which OBO subsets make a term an analysis class.
#[derive(Serialize, Deserialize, Clone, Debug, Default, PartialEq)]
#[serde(default)]
pub struct ClassRules {
    pub subsets: Vec<String>,
    pub subset_suffixes: Vec<String>,
    pub exclude_subset_suffixes: Vec<String>,
    /// Terms that are no class however they are subset: states cells of
    /// many lineages pass through, which would group types across them.
    pub states: Vec<String>,
}

impl ClassRules {
    /// Whether the term `id`, in `subsets`, is a class.
    #[must_use]
    pub fn is_class(&self, id: &str, subsets: &[String]) -> bool {
        !self.states.iter().any(|s| s == id) && self.in_class_subset(subsets)
    }

    fn in_class_subset(&self, subsets: &[String]) -> bool {
        let excluded = subsets.iter().any(|s| {
            self.exclude_subset_suffixes
                .iter()
                .any(|x| s.ends_with(x.as_str()))
        });
        !excluded
            && subsets.iter().any(|s| {
                self.subsets.contains(s)
                    || self.subset_suffixes.iter().any(|x| s.ends_with(x.as_str()))
            })
    }
}

/// Plural handling.
#[derive(Serialize, Deserialize, Clone, Debug, Default, PartialEq)]
#[serde(default)]
pub struct PluralRules {
    /// Trailing words replaced (`cells` → `cell`).
    pub suffixes: Vec<(String, String)>,
    /// A trailing `s` is dropped from words this long or longer (not `ss`);
    /// `None` keeps words as they are.
    pub singular_min_len: Option<usize>,
}

/// The matching rules.
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq)]
#[serde(default)]
pub struct MatchRules {
    /// Where to download the ontology from, when it is not at hand.
    pub ontology_url: Option<String>,
    /// Synonym scope (`EXACT`, `RELATED`, …) → the synonym types that count;
    /// an empty list counts every synonym of that scope.
    pub synonyms: BTreeMap<String, Vec<String>>,
    pub classes: ClassRules,
    pub plural: PluralRules,
    pub implied: ImpliedRules,
    /// Match a label whose words are a term's in another order.
    pub any_word_order: bool,
    /// Marker signs: each sign's word and the symbol that writes it
    /// (`positive` → `+`). However a label writes one (`X-positive`,
    /// `X positive`, `X+`), it is compared as one word, `x-positive`.
    pub marker_signs: BTreeMap<String, String>,
    /// Words a term's name may have that a label leaves off (`alpha-beta`,
    /// a species): a label that matches no name is matched to the most
    /// general of the terms whose names it is but for these words, when
    /// one is above all the others.
    pub qualifiers: Vec<String>,
    pub context: ContextRules,
    pub part: PartRules,
}

/// A word a whole panel leaves off its labels (a tissue word).
#[derive(Serialize, Deserialize, Clone, Debug, Default, PartialEq)]
#[serde(default)]
pub struct ContextRules {
    /// The one word that, with at most an implied one, completes a term's
    /// name or counted synonym for at least this many of a panel's labels no
    /// name matched, and for more than any other word, completes them;
    /// `None` turns it off.
    pub min_labels: Option<usize>,
}

/// Labels nothing matches, placed by the longest run of their words that
/// one term names.
#[derive(Serialize, Deserialize, Clone, Debug, Default, PartialEq)]
#[serde(default)]
pub struct PartRules {
    /// A label goes under its run's term when every word left off is one
    /// the ontology uses for kinds of a term: over the pairs of own names
    /// `X` and `X` with that word more, the longer a kind of `X` (`is_a`) in
    /// at least `min_pairs` pairs and in `min_share` of them. Otherwise the
    /// term only keeps the label to its lineage. `None` turns it off.
    pub min_pairs: Option<usize>,
    pub min_share: f64,
}

/// Words a label may leave off its end.
#[derive(Serialize, Deserialize, Clone, Debug, Default, PartialEq)]
#[serde(default)]
pub struct ImpliedRules {
    /// A label that matches no term is tried again with each of these
    /// appended (`Immature_B` → `immature B cell`), taking the first that
    /// completes exactly one term's name: not a synonym, already a looser
    /// link (`alpha cell` is a retinal ganglion cell's synonym; a pancreas
    /// panel's Alpha is not one). A label with no word completes nothing.
    pub suffixes: Vec<String>,
}

/// Literal matching: names and exact synonyms as written, no classes.
impl Default for MatchRules {
    fn default() -> Self {
        Self {
            ontology_url: None,
            synonyms: BTreeMap::from([("EXACT".to_string(), Vec::new())]),
            classes: ClassRules::default(),
            plural: PluralRules::default(),
            implied: ImpliedRules::default(),
            any_word_order: false,
            marker_signs: BTreeMap::new(),
            qualifiers: Vec::new(),
            context: ContextRules::default(),
            part: PartRules::default(),
        }
    }
}

impl MatchRules {
    /// The rules from layered JSON documents, each overriding the keys it
    /// sets (RFC 7396 merge patch; `_doc` keys are notes). No layers: the
    /// literal default.
    pub fn from_layers(layers: &[Value]) -> anyhow::Result<Self> {
        if layers.is_empty() {
            return Ok(Self::default());
        }
        let mut merged = Value::Object(serde_json::Map::new());
        for l in layers {
            merge_patch(&mut merged, l);
        }
        strip_docs(&mut merged);
        Ok(serde_json::from_value(merged)?)
    }

    /// Whether a synonym of `scope` with `types` counts as a name.
    #[must_use]
    pub fn counts(&self, scope: &str, types: &[&str]) -> bool {
        self.synonyms.get(scope).is_some_and(|want| {
            want.is_empty() || types.iter().any(|t| want.iter().any(|w| w == t))
        })
    }

    /// A label or name as compared: lower case, lupin's separators as one
    /// space, each marker sign one word ([`Self::marker_signs`]), trailing
    /// plural words replaced.
    #[must_use]
    pub fn normalise(&self, s: &str) -> String {
        self.singular_words(self.signed(&label_key(&s.to_lowercase()).replace('_', " ")))
    }

    /// [`Self::normalise`] but for marker signs: a text as written.
    #[must_use]
    pub fn normalise_unsigned(&self, s: &str) -> String {
        self.singular_words(label_key(&s.to_lowercase()).replace('_', " "))
    }

    /// Trailing plural words replaced.
    fn singular_words(&self, mut s: String) -> String {
        for (from, to) in &self.plural.suffixes {
            if let Some(stem) = s.strip_suffix(&format!(" {from}")) {
                s = format!("{stem} {to}");
            } else if s == *from {
                s.clone_from(to);
            }
        }
        s
    }

    /// Each marker sign written `marker-word`: a sign word after its marker
    /// (`x positive`) joins it, and a word ending in a sign's symbol is
    /// marker-sign pairs at every symbol (`x+y-` → `x-positive y-negative`),
    /// but a hyphen is a sign only at the word's end (`x-y+` →
    /// `x-y-positive`); a word that ends in no symbol (`pro-x`) has no sign.
    fn signed(&self, spaced: &str) -> String {
        if self.marker_signs.is_empty() {
            return spaced.to_string();
        }
        // Each one-character symbol and its sign.
        let symbols: Vec<(char, &str)> = (self.marker_signs.iter())
            .filter_map(|(sign, sym)| {
                let mut c = sym.chars();
                c.next()
                    .filter(|_| c.next().is_none())
                    .map(|c| (c, sign.as_str()))
            })
            .collect();
        let sign_of = |c: char| symbols.iter().find(|(s, _)| *s == c).map(|(_, w)| *w);
        // Sign the last word, unless it has a sign already.
        let sign_last = |out: &mut Vec<String>, sign: &str| match out.last_mut() {
            Some(last) if !self.has_sign(last) => {
                last.push('-');
                last.push_str(sign);
                true
            }
            _ => false,
        };
        let mut out: Vec<String> = Vec::new();
        for w in spaced.split(' ').filter(|w| !w.is_empty()) {
            if self.marker_signs.contains_key(w) && sign_last(&mut out, w) {
                continue;
            }
            if w.chars().last().and_then(sign_of).is_none() {
                out.push(w.to_string());
                continue;
            }
            // A hyphen symbol inside a word is the marker's own (`x-y+`);
            // only at the end is it a sign.
            let mut marker = String::new();
            let n = w.chars().count();
            for (i, c) in w.chars().enumerate() {
                match sign_of(c).filter(|_| c != '-' || i + 1 == n) {
                    None => marker.push(c),
                    // A symbol alone signs the word before it (`x +`).
                    Some(sign) if marker.is_empty() => {
                        sign_last(&mut out, sign);
                    }
                    Some(sign) => out.push(format!("{}-{sign}", std::mem::take(&mut marker))),
                }
            }
        }
        out.join(" ")
    }

    /// A normalised text's words: spaces part them, and hyphens too, but
    /// not in a marker with its sign (`x-y-positive` stays one word).
    pub fn words<'a>(&'a self, normalised: &'a str) -> impl Iterator<Item = &'a str> + 'a {
        normalised
            .split(' ')
            .filter(|w| !w.is_empty())
            .flat_map(move |w| {
                let whole = self.has_sign(w).then_some(w);
                let parts = whole
                    .is_none()
                    .then(|| w.split('-').filter(|p| !p.is_empty()));
                whole.into_iter().chain(parts.into_iter().flatten())
            })
    }

    /// Whether `w` is a marker with its sign (`x-positive`).
    #[must_use]
    pub fn has_sign(&self, w: &str) -> bool {
        (self.marker_signs.keys())
            .any(|s| (w.strip_suffix(s.as_str())).is_some_and(|m| m.len() > 1 && m.ends_with('-')))
    }

    /// [`Self::words`] in sorted order, for matching words in any order.
    #[must_use]
    pub fn word_bag(&self, normalised: &str) -> String {
        let mut words: Vec<&str> = self.words(normalised).collect();
        words.sort_unstable();
        words.join(" ")
    }

    /// Each of [`Self::qualifiers`] as its words, compared as names are.
    #[must_use]
    pub fn qualifier_words(&self) -> Vec<Vec<String>> {
        (self.qualifiers.iter())
            .map(|q| {
                self.words(&self.normalise(q))
                    .map(String::from)
                    .collect::<Vec<_>>()
            })
            .filter(|q| !q.is_empty())
            .collect()
    }

    /// A normalised name's [`Self::words`] without its `qualifiers` (from
    /// [`Self::qualifier_words`]), when it has any and other words besides.
    #[must_use]
    pub fn unqualified(&self, normalised: &str, qualifiers: &[Vec<String>]) -> Option<String> {
        let words: Vec<&str> = self.words(normalised).collect();
        let (mut core, mut rest) = (Vec::new(), &words[..]);
        while let Some(w) = rest.first() {
            match (qualifiers.iter())
                .find(|q| q.len() <= rest.len() && q.iter().zip(rest).all(|(a, b)| a == b))
            {
                Some(q) => rest = &rest[q.len()..],
                None => {
                    core.push(*w);
                    rest = &rest[1..];
                }
            }
        }
        (!core.is_empty() && core.len() < words.len()).then(|| core.join(" "))
    }

    /// Every word made singular, when the rules say so.
    #[must_use]
    pub fn singular(&self, normalised: &str) -> String {
        let Some(min) = self.plural.singular_min_len else {
            return normalised.to_string();
        };
        normalised
            .split(' ')
            .map(|w| match w.strip_suffix('s') {
                Some(stem) if w.len() >= min && !stem.ends_with('s') => stem,
                _ => w,
            })
            .collect::<Vec<_>>()
            .join(" ")
    }
}

/// `patch` onto `target`, RFC 7396: objects merge key by key, `null` removes,
/// anything else replaces.
pub fn merge_patch(target: &mut Value, patch: &Value) {
    match patch {
        Value::Object(p) => {
            if !target.is_object() {
                *target = Value::Object(serde_json::Map::new());
            }
            if let Value::Object(t) = target {
                for (k, v) in p {
                    if v.is_null() {
                        t.remove(k);
                    } else {
                        merge_patch(t.entry(k.clone()).or_insert(Value::Null), v);
                    }
                }
            }
        }
        other => *target = other.clone(),
    }
}

fn strip_docs(v: &mut Value) {
    if let Value::Object(m) = v {
        m.retain(|k, _| !k.starts_with('_'));
        m.values_mut().for_each(strip_docs);
    }
}

/// Curated label → CL term aliases.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Aliases {
    /// [`label_key`], lower case → term, and when the row was taken.
    map: BTreeMap<String, (String, usize)>,
    taken: usize,
}

impl Aliases {
    /// Add the rows of an alias table, overriding earlier ones:
    /// `label<TAB>CL:id[<TAB>note]`, or `label,CL:id` as `--label-cl` maps
    /// have been written, read as every two-column table is
    /// ([`TwoFieldRows`]). `#`/`%` comments and blank lines are skipped, and so,
    /// with a warning naming `source`, is any row without a CL-style id (a
    /// header, `NA`). Returns how many rows were taken.
    pub fn add_tsv(&mut self, text: &str, source: &str) -> usize {
        let rows = TwoFieldRows::from_text(text);
        let (mut n, mut skipped) = (0, rows.gapped);
        for (label, id) in &rows.rows {
            // The id's first word, as a note may follow it (`CL:x, note`).
            let id = id.split([',', ' ']).next().unwrap_or_default();
            if !id.contains(':') {
                skipped += 1;
                continue;
            }
            self.taken += 1;
            self.map.insert(key(label), (id.to_string(), self.taken));
            n += 1;
        }
        if skipped > 0 {
            log::warn!("{source}: skipped {skipped} row(s) without a `label<TAB>CL:id` pair");
        }
        n
    }

    /// The term `label` is aliased to.
    #[must_use]
    pub fn get(&self, label: &str) -> Option<&str> {
        self.map.get(&key(label)).map(|(id, _)| id.as_str())
    }

    /// Each aliased label (as [`label_key`] keys it, lower case) and its
    /// term, in the order the rows were taken: a later one overrides.
    pub fn iter(&self) -> impl Iterator<Item = (&str, &str)> {
        let mut rows: Vec<_> = self.map.iter().collect();
        rows.sort_by_key(|(_, (_, at))| *at);
        rows.into_iter()
            .map(|(l, (id, _))| (l.as_str(), id.as_str()))
    }

    /// How many labels are aliased.
    #[must_use]
    pub fn len(&self) -> usize {
        self.map.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }
}

/// How aliases compare labels: as lupin labels, ignoring case.
fn key(label: &str) -> String {
    label_key(&label.to_lowercase())
}

/// The rules the repository ships (`data/cl_matching.json`), for tests.
#[cfg(test)]
pub(crate) fn shipped() -> MatchRules {
    let v: Value = serde_json::from_str(include_str!("../../data/cl_matching.json"))
        .expect("the shipped rules are JSON");
    MatchRules::from_layers(&[v]).expect("the shipped rules parse")
}

#[cfg(test)]
#[path = "tests/cl_rules.rs"]
mod tests;
