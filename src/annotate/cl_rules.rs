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
    /// Subsets whose terms are lineages: groups of the first round below
    /// one lineage are one group.
    pub lineage_suffixes: Vec<String>,
    /// Subsets too general to be a lineage.
    pub lineage_exclude: Vec<String>,
}

impl ClassRules {
    /// Whether a term in `subsets` is a lineage.
    #[must_use]
    pub fn is_lineage(&self, subsets: &[String]) -> bool {
        subsets.iter().any(|s| {
            !self.lineage_exclude.contains(s)
                && self
                    .lineage_suffixes
                    .iter()
                    .any(|x| s.ends_with(x.as_str()))
        })
    }

    /// Whether a term in `subsets` is a class.
    #[must_use]
    pub fn is_class(&self, subsets: &[String]) -> bool {
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
    pub qualifiers: QualifierRules,
    /// As a last resort, compare without spaces, nor a hyphen between
    /// words (`Pro-X blast` for `proXblast`); a hyphen [`part_hyphens`]
    /// keeps stays.
    pub compact: bool,
    pub contained: ContainedRules,
    pub context: ContextRules,
    pub part: PartRules,
}

/// The trailing qualifiers a term is also known without: a taxon in
/// parentheses (`neuron (sensu Vertebrata)`, `cell (Mmus)`), or one of
/// `comma_words` after a comma (`cell, human`); no other word after a comma
/// is a qualifier (`amacrine cell, GABAergic` is a kind of amacrine cell).
/// A text that then names a term and its own specialisations takes the
/// most general of them.
#[derive(Serialize, Deserialize, Clone, Debug, Default, PartialEq)]
#[serde(default)]
pub struct QualifierRules {
    pub parenthetical: bool,
    /// The taxa a comma may qualify a name by, compared ignoring case; the
    /// older switch `comma_word` reads as the shipped taxa or none.
    pub comma_words: Vec<String>,
}

/// A label whose words a term's own name (not a synonym) spells out, in any order and among
/// others (`CD4 Central Memory` in `central memory CD4-positive, alpha-beta
/// T cell`; `X-positive` and `X+` are both `X`), takes the most general
/// such term: one every other such term is a kind of, else the one whose
/// name adds the fewest words, else none.
#[derive(Serialize, Deserialize, Clone, Debug, Default, PartialEq)]
#[serde(default)]
pub struct ContainedRules {
    /// The fewest words, not counting an implied one, a label needs for this;
    /// `None` turns it off.
    pub min_words: Option<usize>,
}

/// A label nothing matches, placed by the longest run of its words that
/// matches one term (`Early GMP` by `GMP`): under that term when every word
/// left off is a modifier, one the ontology uses for kinds of a term (it
/// names `X` and `word X` with the second a kind of the first in at least
/// `min_pairs` pairs, and in at least `min_share` of all such pairs:
/// `early`, `mature`, `cycling`; not `pre`, `progenitor`); else the term is
/// only a hint of its lineage, which a type joining a group by its markers
/// keeps to (`Pre-cDC` joins no group of another lineage). The panel's word
/// ([`ContextRules`]) completes a run only when no run of that length
/// matches as written.
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq)]
#[serde(default)]
pub struct PartRules {
    /// `None` turns placing by part of a label off.
    pub min_pairs: Option<usize>,
    pub min_share: f64,
}

impl Default for PartRules {
    fn default() -> Self {
        Self {
            min_pairs: None,
            min_share: 1.0,
        }
    }
}

/// A word the whole panel leaves off its labels (a tissue: `pancreatic`
/// for `Alpha`, `Beta`, `Ductal`), found as the one word that, with an
/// implied one, completes a term's name or counted synonym for the most
/// labels no name matched (an aliased one counting only when the word
/// completes it to its alias's term); it then completes the labels nothing
/// else matched.
#[derive(Serialize, Deserialize, Clone, Debug, Default, PartialEq)]
#[serde(default)]
pub struct ContextRules {
    /// How many labels the word must complete; `None` turns it off.
    pub min_labels: Option<usize>,
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
            qualifiers: QualifierRules::default(),
            compact: false,
            contained: ContainedRules::default(),
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
        comma_switch(&mut merged);
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
    /// space, trailing plural words replaced.
    #[must_use]
    pub fn normalise(&self, s: &str) -> String {
        let mut s = label_key(&s.to_lowercase()).replace('_', " ");
        for (from, to) in &self.plural.suffixes {
            if let Some(stem) = s.strip_suffix(&format!(" {from}")) {
                s = format!("{stem} {to}");
            } else if s == *from {
                s.clone_from(to);
            }
        }
        s
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

/// Rules written before the comma's taxa were listed (a pass's recorded
/// files, a user's layer) say `"comma_word": true` or `false`: the shipped
/// taxa, or none. The switch, where a layer sets it, decides.
fn comma_switch(rules: &mut Value) {
    let Some(q) = rules.get_mut("qualifiers").and_then(Value::as_object_mut) else {
        return;
    };
    let Some(on) = q.remove("comma_word").as_ref().and_then(Value::as_bool) else {
        return;
    };
    let taxa = if on {
        let shipped: Value = serde_json::from_str(SHIPPED).expect("the shipped rules are JSON");
        shipped["qualifiers"]["comma_words"].clone()
    } else {
        Value::Array(Vec::new())
    };
    q.insert("comma_words".into(), taxa);
}

/// The rules the repository ships, `data/cl_matching.json`.
const SHIPPED: &str = include_str!("../../data/cl_matching.json");

fn strip_docs(v: &mut Value) {
    if let Value::Object(m) = v {
        m.retain(|k, _| !k.starts_with('_'));
        m.values_mut().for_each(strip_docs);
    }
}

/// Curated label → CL term aliases.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Aliases {
    /// [`label_key`], lower case → term.
    map: BTreeMap<String, String>,
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
            self.map.insert(key(label), id.to_string());
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
        self.map.get(&key(label)).map(String::as_str)
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

/// How aliases compare labels: as lupin labels, ignoring case, with
/// hyphens read by [`part_hyphens`].
fn key(label: &str) -> String {
    label_key(&part_hyphens(&label.to_lowercase()))
}

/// `text` with each hyphen between two words of letters alone read as a
/// space (`Pro-B` is `Pro B`, `T-cell` is `T cell`, `alpha-beta` is `alpha
/// beta`), and so each before a marker's sign written out as the ontology
/// writes it ([`SIGN_WORDS`]: `CD4-positive` is `CD4 positive`); every other
/// hyphen is kept as written, so a negative marker or a name with digits
/// stays its literal text (`CD34-` is not `CD34`, `CD34-CD38-` is not
/// `CD34 CD38-`, `CD4- T` is not `CD4 T`, `NK-CD56bright` and `B-1a` are as
/// written). Aliases and term names are compared under this one rule.
pub(crate) fn part_hyphens(text: &str) -> String {
    let chars: Vec<char> = text.chars().collect();
    (chars.iter().enumerate())
        .map(|(i, &c)| {
            if c != '-' {
                return c;
            }
            let before = word(chars[..i].iter().rev().copied());
            let after = word(chars[i + 1..].iter().copied());
            let letters = |w: &str| !w.is_empty() && w.chars().all(char::is_alphabetic);
            let parts = (letters(&before) && letters(&after))
                || (!before.is_empty() && SIGN_WORDS.contains(&after.as_str()));
            if parts {
                ' '
            } else {
                c
            }
        })
        .collect()
}

/// The words the Cell Ontology writes a marker's sign with
/// (`CD4-positive`, `CD16-negative`), which a label may write with a space.
pub(crate) const SIGN_WORDS: [&str; 2] = ["positive", "negative"];

/// The word `chars` start with: up to the first character that is neither a
/// letter nor a digit.
fn word(chars: impl Iterator<Item = char>) -> String {
    chars.take_while(|c| c.is_alphanumeric()).collect()
}

/// The rules the repository ships (`data/cl_matching.json`), for tests.
#[cfg(test)]
pub(crate) fn shipped() -> MatchRules {
    let v: Value = serde_json::from_str(SHIPPED).expect("the shipped rules are JSON");
    MatchRules::from_layers(&[v]).expect("the shipped rules parse")
}

#[cfg(test)]
#[path = "tests/cl_rules.rs"]
mod tests;
