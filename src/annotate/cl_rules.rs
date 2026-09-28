//! How panel labels meet the Cell Ontology, as data: the matching rules
//! (`cl_matching.json`) and the curated label → term aliases
//! (`cl_aliases.tsv`). Nothing here is decided in code: without a rules file
//! matching is literal (a term's name or exact synonym, as written), and every
//! heuristic (plurals, word order, abbreviations) comes from a file. Finding
//! and layering the files is [`crate::manifest::data_files`]'s job.

use crate::annotate::markers::label_key;
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
}

impl ClassRules {
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
    /// Match a label whose words are a term's in another order.
    pub any_word_order: bool,
}

/// Literal matching: names and exact synonyms as written, no classes.
impl Default for MatchRules {
    fn default() -> Self {
        Self {
            ontology_url: None,
            synonyms: BTreeMap::from([("EXACT".to_string(), Vec::new())]),
            classes: ClassRules::default(),
            plural: PluralRules::default(),
            any_word_order: false,
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
    /// have been written. `#` comments and blank lines are skipped, and so,
    /// with a warning naming `source`, is any row without a CL-style id (a
    /// header, `NA`). Returns how many rows were taken.
    pub fn add_tsv(&mut self, text: &str, source: &str) -> usize {
        let (mut n, mut skipped) = (0, 0);
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let (label, rest) = line
                .split_once('\t')
                .or_else(|| line.split_once(','))
                .unwrap_or((line, ""));
            let id = rest.split(['\t', ',']).next().unwrap_or_default().trim();
            let label = label.trim();
            if label.is_empty() || !id.contains(':') {
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
