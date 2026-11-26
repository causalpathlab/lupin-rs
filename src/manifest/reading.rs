//! Readings: a reader's label for a cluster, kept as evidence beside the
//! methods' readings and never as a decision on the cluster. The reader is
//! typically an agent given only the cluster's most specific genes (`lupin
//! ask --blind`), which shares nothing with the marker panel the methods
//! read: its agreement with them is a check from outside the panel.
//!
//! A reading arrives as a `reading` decision through `lupin relabel`, so it
//! takes the chain's lock and lands in the next round with the round's
//! history. There it is resolved once: placed on the Cell Ontology (the term
//! the reader gave, when the ontology has it, else its label by the matching
//! rules), and the panel types it agrees with recorded. A reading agrees
//! with a type whose term is its own or above it (the reader is at least as
//! specific), and with a type below it within the type's own coarse group
//! (a broader reading names the group, not something above it). Readings of
//! clusters a later round merges are dropped: the reader never saw the
//! merged cells.
//!
//! A reading is blind when it cites the id of the blind prompt written for
//! the round it reads ([`BLIND_PROMPT`]); any other is kept as not blind.

use crate::annotate::celltype_tree::{ClTerms, TypeTree};
use crate::annotate::markers::label_key;
use crate::annotate::rounds::{Action, ClusterId, DecidedBy, Decision, Table};
use crate::manifest::run::Loaded;
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;

/// `{out}.agent_readings.jsonl`: a round's readings, one per line
/// (`annotate.agent_readings`).
pub const READINGS: &str = ".agent_readings.jsonl";

/// `{round}.blind_prompt.md`: the blind prompt written for a round, beside
/// it; its last line gives its id, [`prompt_id`] of the text above it.
pub const BLIND_PROMPT: &str = ".blind_prompt.md";

/// What starts a blind prompt's last line, its id.
pub const PROMPT_ID: &str = "\nPrompt id: ";

/// One reading, resolved.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Reading {
    pub cluster: ClusterId,
    /// The reader's label, as given.
    pub label: String,
    /// The Cell Ontology term it is placed on.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cl: Option<String>,
    /// The panel types it agrees with ([`label_key`]s).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub agrees_with: Vec<String>,
    /// Whether it answers the blind prompt written for the round it reads.
    pub blind: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub prompt: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub model: Option<String>,
    /// The genes the reader relied on.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub genes: Vec<String>,
    pub rationale: String,
    pub decided_by: DecidedBy,
    /// The round it was recorded in.
    pub round: String,
    pub timestamp: String,
}

/// A prompt's id: the start of its text's SHA-256, in hex, so a reading
/// made from a prompt can be checked against it later.
#[must_use]
pub fn prompt_id(text: &str) -> String {
    use sha2::{Digest, Sha256};
    crate::hex(&Sha256::digest(text.as_bytes())[..8])
}

/// Where the blind prompt for the round `loaded` is written.
#[must_use]
pub fn blind_prompt_path(loaded: &Loaded) -> String {
    let prefix = super::run::derive_out_prefix(&loaded.file.to_string_lossy());
    format!("{prefix}{BLIND_PROMPT}")
}

/// The id of the blind prompt written for the round `loaded`, when there is
/// one beside it whose text is still the text its id names.
fn blind_prompt_id(loaded: &Loaded) -> Option<String> {
    let text = fs::read_to_string(blind_prompt_path(loaded)).ok()?;
    let (body, id) = text.rsplit_once(PROMPT_ID)?;
    let id = id.trim();
    (prompt_id(body) == id).then(|| id.to_string())
}

/// The readings the round `loaded` records; none when it records none.
pub(crate) fn read(loaded: &Loaded) -> Result<Vec<Reading>> {
    let Some(rel) = loaded.manifest.annotate.agent_readings.as_deref() else {
        return Ok(Vec::new());
    };
    let Some(path) = loaded.recorded(rel) else {
        log::warn!("{rel} is recorded but gone: no readings");
        return Ok(Vec::new());
    };
    let text = fs::read_to_string(&path).with_context(|| format!("reading {path}"))?;
    super::rounds::numbered(&text)
        .filter(|(_, l)| !l.trim().is_empty())
        .map(|(n, l)| {
            serde_json::from_str(l).with_context(|| format!("{path}: line {n} is not a reading"))
        })
        .collect()
}

/// The readings that count among `rs`: the blind ones, when there are any.
pub(crate) fn counted<'a>(rs: impl IntoIterator<Item = &'a Reading>) -> Vec<&'a Reading> {
    let mut rs: Vec<&Reading> = rs.into_iter().collect();
    if rs.iter().any(|r| r.blind) {
        rs.retain(|r| r.blind);
    }
    rs
}

/// One panel type on the Cell Ontology.
struct TypeTerm {
    key: String,
    /// Its own term and every term above it.
    above: Option<(String, BTreeSet<String>)>,
    /// Its coarse group's term.
    group: Option<String>,
}

/// The panel's types on the Cell Ontology, to tell which a reading agrees
/// with.
struct PanelTerms<'a> {
    types: Vec<TypeTerm>,
    terms: Option<&'a ClTerms>,
}

impl<'a> PanelTerms<'a> {
    fn of(
        markers: &[(String, String)],
        tree: Option<&TypeTree>,
        terms: Option<&'a ClTerms>,
    ) -> Self {
        let keys: Vec<String> = crate::annotate::celltype_tree::panel_types(markers)
            .iter()
            .map(|t| label_key(t))
            .collect();
        let mut own: BTreeMap<String, String> = BTreeMap::new();
        let mut group: BTreeMap<String, String> = BTreeMap::new();
        if let Some(tree) = tree {
            own.extend(
                tree.label_cl
                    .iter()
                    .map(|(l, id)| (label_key(l), id.clone())),
            );
            for g in &tree.groups {
                let Some(id) = &g.cl_id else { continue };
                group.extend(g.members.iter().map(|m| (label_key(m), id.clone())));
            }
        }
        if let Some(t) = terms {
            let missing: Vec<&str> = (keys.iter())
                .filter(|k| !own.contains_key(*k))
                .map(String::as_str)
                .collect();
            let placed = t.place_labels(missing.iter().copied());
            for k in missing {
                if let Some(id) = placed.term(k) {
                    own.insert(k.to_string(), id.clone());
                }
            }
        }
        let types = keys
            .into_iter()
            .map(|key| TypeTerm {
                above: (own.get(&key).cloned())
                    .zip(terms)
                    .map(|(id, t)| (id.clone(), t.ancestors_or_self(&id))),
                group: group.get(&key).cloned(),
                key,
            })
            .collect();
        Self { types, terms }
    }

    /// The types a reading labelled `label`, on `cl`, agrees with.
    fn agreeing(&self, label: &str, cl: Option<&str>) -> Vec<String> {
        // A reader's spelling of a type's name, whatever its case.
        let name = label_key(label).to_lowercase();
        let above_c = cl.zip(self.terms).map(|(c, t)| (c, t.ancestors_or_self(c)));
        (self.types.iter())
            .filter(|t| {
                let by_term = match (&above_c, &t.above) {
                    (Some((c, above_c)), Some((own, above_own))) => {
                        // As specific as the type, or more; or broader,
                        // but within the type's own group.
                        above_c.contains(own)
                            || (above_own.contains(*c)
                                && t.group.as_ref().is_some_and(|g| above_c.contains(g)))
                    }
                    _ => false,
                };
                by_term || t.key.to_lowercase() == name
            })
            .map(|t| t.key.clone())
            .collect()
    }
}

/// The next round's readings: `source`'s and the `reading`s among
/// `decisions` (made in round `round`), of the clusters in `present`, each
/// placed on the Cell Ontology and found to agree with types of the panel
/// `markers`, which the round may have edited. The ontology is read once,
/// and only when there is a reading.
pub(crate) fn next_readings(
    source: &Loaded,
    decisions: &[Decision],
    markers: &[(String, String)],
    round: &str,
    present: &BTreeSet<ClusterId>,
) -> Result<Vec<Reading>> {
    let older = read(source)?;
    let new: Vec<&Decision> = (decisions.iter())
        .filter(|d| d.action == Action::Reading)
        .collect();
    if older.is_empty() && new.is_empty() {
        return Ok(Vec::new());
    }
    let data = super::recalibrate::round_cl_data(source)?;
    let terms = data.terms()?;
    let tree = super::rounds::read_tree(source)?;
    let panel = PanelTerms::of(markers, tree.as_ref(), terms);
    let blind_id = blind_prompt_id(source);
    let placed = new.into_iter().map(|d| {
        let label = d.label.clone().unwrap_or_default();
        let given = (d.reading.cl.as_deref())
            .map(str::trim)
            .filter(|c| !c.is_empty());
        let cl = match (given, terms) {
            (Some(c), Some(t)) if t.has(c) => Some(c.to_string()),
            (_, Some(t)) => t.place_labels([label.as_str()]).term(&label).cloned(),
            (given, None) => given.map(str::to_string),
        };
        Reading {
            cluster: d.cluster[0],
            label,
            cl,
            agrees_with: Vec::new(),
            blind: d.reading.prompt.is_some() && d.reading.prompt == blind_id,
            prompt: d.reading.prompt.clone(),
            model: d.reading.model.clone(),
            genes: d.features.clone(),
            rationale: d.rationale.clone(),
            decided_by: d.decided_by,
            round: round.to_string(),
            timestamp: d.timestamp.clone().unwrap_or_default(),
        }
    });
    let mut readings = carried(older, placed.collect(), present);
    for r in &mut readings {
        r.agrees_with = panel.agreeing(&r.label, r.cl.as_deref());
    }
    Ok(readings)
}

/// `older` readings, then `newer`, of clusters still in `present`: a reader
/// never saw a merged cluster's cells, even when the merge comes in the
/// batch that brings its reading.
fn carried(
    older: Vec<Reading>,
    newer: Vec<Reading>,
    present: &BTreeSet<ClusterId>,
) -> Vec<Reading> {
    older
        .into_iter()
        .chain(newer)
        .filter(|r| present.contains(&r.cluster))
        .collect()
}

/// Write `readings` at `{out}{READINGS}`; `None` (nothing written) when
/// there are none.
pub(crate) fn write(out: &str, readings: &[Reading]) -> Result<Option<String>> {
    if readings.is_empty() {
        return Ok(None);
    }
    let path = format!("{out}{READINGS}");
    let lines: Vec<Box<str>> = (readings.iter())
        .map(|r| serde_json::to_string(r).map(String::into_boxed_str))
        .collect::<Result<_, _>>()?;
    legume_numeric::matrix::common_io::write_lines(&lines, &path)?;
    log::info!("wrote {path}");
    Ok(Some(path))
}

/// Each read cluster's share of its readings that agree with each type: the
/// panel types they agree with, and for a reading that agrees with none, its
/// own label. Blind readings only, when a cluster has any.
#[must_use]
pub(crate) fn shares(readings: &[Reading]) -> Table {
    let mut by: BTreeMap<ClusterId, Vec<&Reading>> = BTreeMap::new();
    for r in readings {
        by.entry(r.cluster).or_default().push(r);
    }
    for rs in by.values_mut() {
        *rs = counted(rs.iter().copied());
    }
    let named = |r: &Reading| -> Vec<String> {
        if r.agrees_with.is_empty() {
            vec![label_key(&r.label)]
        } else {
            r.agrees_with.clone()
        }
    };
    let cols: Vec<String> = by
        .values()
        .flatten()
        .flat_map(|r| named(r))
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut values = Vec::with_capacity(by.len() * cols.len());
    for rs in by.values() {
        let n = rs.len() as f32;
        let mut count = vec![0usize; cols.len()];
        for r in rs {
            for c in named(r) {
                if let Ok(j) = cols.binary_search(&c) {
                    count[j] += 1;
                }
            }
        }
        values.extend(count.into_iter().map(|c| c as f32 / n));
    }
    Table {
        rows: by.keys().copied().collect(),
        cols,
        values,
    }
}

/// A reading of `cluster` labelled `label`, agreeing with `agrees`, for
/// tests.
#[cfg(test)]
pub(crate) fn reading(cluster: ClusterId, label: &str, agrees: &[&str], blind: bool) -> Reading {
    Reading {
        cluster,
        label: label.into(),
        cl: None,
        agrees_with: agrees.iter().map(|s| (*s).to_string()).collect(),
        blind,
        prompt: None,
        model: None,
        genes: Vec::new(),
        rationale: "r".into(),
        decided_by: DecidedBy::Agent,
        round: "r1".into(),
        timestamp: String::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const OBO: &str = "[Term]\nid: CL:4\nname: leukocyte\n\n\
[Term]\nid: CL:1\nname: T cell\nis_a: CL:4\n\n\
[Term]\nid: CL:2\nname: CD4 T cell\nis_a: CL:1\n\n\
[Term]\nid: CL:3\nname: CD8 T cell\nis_a: CL:1\n\n\
[Term]\nid: CL:5\nname: naive CD4 T cell\nis_a: CL:2\n";

    /// Panel types CD4 and CD8 on their terms, in the group `T cell`.
    fn t_cells() -> (Vec<(String, String)>, TypeTree, ClTerms) {
        use crate::annotate::celltype_tree::{TreeSource, TypeGroup};
        let tree = TypeTree {
            source: TreeSource::CellOntology,
            release: None,
            label_cl: [("CD4".into(), "CL:2".into()), ("CD8".into(), "CL:3".into())].into(),
            groups: vec![TypeGroup {
                name: "T cell".into(),
                cl_id: Some("CL:1".into()),
                members: vec!["CD4".into(), "CD8".into()],
            }],
        };
        let terms = ClTerms::parse(OBO, &crate::annotate::cl_rules::shipped());
        let markers = [("G1", "CD4"), ("G2", "CD8")]
            .map(|(g, t)| (g.to_string(), t.to_string()))
            .into();
        (markers, tree, terms)
    }

    #[test]
    fn a_reading_agrees_with_its_lineage_within_the_types_group() {
        let (markers, tree, terms) = t_cells();
        let panel = PanelTerms::of(&markers, Some(&tree), Some(&terms));
        let agree = |cl: &str| panel.agreeing("x", Some(cl));
        // As specific as a type, or more.
        assert_eq!(agree("CL:2"), ["CD4"]);
        assert_eq!(agree("CL:5"), ["CD4"]);
        // Broader, but the types' own group.
        assert_eq!(agree("CL:1"), ["CD4", "CD8"]);
        // Above the group: it names nothing the panel tells apart.
        assert!(agree("CL:4").is_empty());
    }

    #[test]
    fn without_the_ontology_a_reading_agrees_by_name_only() {
        let (markers, tree, _) = t_cells();
        let panel = PanelTerms::of(&markers, Some(&tree), None);
        assert_eq!(panel.agreeing("cd4", Some("CL:1")), ["CD4"]);
        assert!(panel.agreeing("T cell", Some("CL:1")).is_empty());
    }

    #[test]
    fn prompt_ids_are_stable_and_distinct() {
        assert_eq!(prompt_id("a"), prompt_id("a"));
        assert_ne!(prompt_id("a"), prompt_id("b"));
        assert_eq!(prompt_id("a").len(), 16);
    }

    #[test]
    fn shares_count_runs_and_prefer_blind_readings() {
        let rs = [
            reading(0, "CT1", &["CT1"], true),
            reading(0, "CT1", &["CT1"], true),
            reading(0, "Other", &[], true),
            reading(0, "CT2", &["CT2"], false),
            reading(1, "Broad", &["CT1", "CT2"], false),
        ];
        let t = shares(&rs);
        assert_eq!(t.rows, vec![0, 1]);
        assert_eq!(t.cols, vec!["CT1", "CT2", "Other"]);
        let k0 = t.row(0).unwrap();
        assert!((k0[0] - 2.0 / 3.0).abs() < 1e-6);
        assert_eq!(k0[1], 0.0);
        assert!((k0[2] - 1.0 / 3.0).abs() < 1e-6);
        assert_eq!(t.row(1).unwrap(), &[1.0, 1.0, 0.0]);
    }

    #[test]
    fn merged_clusters_lose_their_readings() {
        let present: BTreeSet<ClusterId> = [1].into();
        let kept = carried(
            vec![reading(0, "CT1", &[], true), reading(1, "CT2", &[], true)],
            vec![reading(1, "CT3", &[], true), reading(0, "CT4", &[], true)],
            &present,
        );
        let labels: Vec<&str> = kept.iter().map(|r| r.label.as_str()).collect();
        assert_eq!(labels, ["CT2", "CT3"]);
    }
}
