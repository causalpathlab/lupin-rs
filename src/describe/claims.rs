//! What a model's text says of its evidence, checked against it. The
//! citation check ([`super::citation_check`]) polices which names a text
//! uses; this polices the claims it makes of them, in the words its prompt
//! states the facts in: a method ranking a type first or second, a type
//! sitting between two others or before one, an edge's verdict, a root,
//! branch point or end. A claim is read from the clause it is in (split
//! at commas and `and`, `but`, `while`, `with`, `which`), about the types
//! that clause names, else, every kind alike, the type a relative clause
//! speaks of, the sentence's subject or the text's; a claim no word of the
//! vocabulary makes is not checked.

use super::consensus::ClusterInput;
use super::norm;
use super::trajectory::OrderFacts;
use crate::annotate::rounds::Method;
use crate::trajectory::edges::Verdict;
use std::collections::BTreeMap;

/// Ok, or the claim the facts do not bear out, in words.
pub type Claim = Result<(), String>;

/// Words that join clauses, a claim read in each on its own.
const JOINS: [&str; 7] = ["and", "but", "while", "whereas", "with", "which", "where"];
/// Those of [`JOINS`] that lead a relative clause, which speaks of a type
/// named before it.
const RELATIVE: [&str; 2] = ["which", "where"];
/// Words for the place a type holds: a lineage's root, a branch point or
/// an end.
const ROOT: [&str; 11] = [
    "root",
    "roots",
    "rooted",
    "start",
    "starts",
    "starting",
    "begin",
    "begins",
    "base",
    "origin",
    "originates",
];
const BRANCH: [&str; 9] = [
    "branch",
    "branches",
    "branching",
    "branched",
    "split",
    "splits",
    "diverge",
    "diverges",
    "bifurcates",
];
const END: [&str; 9] = [
    "end",
    "ends",
    "ending",
    "ended",
    "terminal",
    "terminates",
    "terminate",
    "endpoint",
    "final",
];
/// Words putting the type before them earlier than the one after.
const FORWARD: [&str; 18] = [
    "gives rise to",
    "give rise to",
    "giving rise to",
    "gave rise to",
    "leads to",
    "lead to",
    "leading to",
    "precedes",
    "precede",
    "preceding",
    "before",
    "upstream of",
    "differentiates into",
    "differentiate into",
    "develops into",
    "develop into",
    "becomes",
    "become",
];
/// Words putting it later.
const BACKWARD: [&str; 15] = [
    "after",
    "follows",
    "follow",
    "following",
    "downstream of",
    "derived from",
    "derives from",
    "arises from",
    "arise from",
    "comes from",
    "come from",
    "develops from",
    "develop from",
    "stems from",
    "stem from",
];
/// Words for a method's first.
const FIRST: [&str; 8] = [
    "first",
    "top",
    "highest",
    "strongest",
    "best",
    "leading",
    "greatest",
    "likeliest",
];
const ORDINALS: [&str; 9] = [
    "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth",
];

/// Where a text names one of its entities.
#[derive(Clone, Debug)]
struct Span {
    start: usize,
    end: usize,
    entity: String,
}

/// Whether `t[i..j]` stands as whole words.
fn whole(t: &str, i: usize, j: usize) -> bool {
    let edge = |c: Option<char>| c.is_none_or(|c| !c.is_alphanumeric());
    edge(t[..i].chars().next_back()) && edge(t[j..].chars().next())
}

/// Where `t` has `w` as whole words, a plural `s` after it allowed when
/// `plural`; each occurrence's start and end.
pub(crate) fn find(t: &str, w: &str, plural: bool) -> Vec<(usize, usize)> {
    t.match_indices(w)
        .filter_map(|(i, _)| {
            let j = i + w.len();
            if whole(t, i, j) {
                Some((i, j))
            } else if plural && t[j..].starts_with('s') && whole(t, i, j + 1) {
                Some((i, j + 1))
            } else {
                None
            }
        })
        .collect()
}

/// The ways a text names each of `entities`: its name; and, when they are
/// no other entity's, its trailing words (`monocyte` for `CT14-positive
/// monocyte`, not a lone `cell`) and either without a last `cell`
/// (`narrow progenitor` for `narrow progenitor cell`).
fn forms(entities: &[String]) -> Vec<(String, String)> {
    const GENERIC: [&str; 4] = ["cell", "cells", "type", "types"];
    let mut short: BTreeMap<String, Vec<&String>> = BTreeMap::new();
    let names: BTreeMap<String, &String> = entities.iter().map(|e| (norm(e), e)).collect();
    for (n, e) in &names {
        let words: Vec<&str> = n.split_whitespace().collect();
        let mut tails: Vec<String> = (1..words.len())
            .map(|k| words[k..].join(" "))
            .filter(|f| f.split(' ').any(|w| !GENERIC.contains(&w)))
            .filter(|f| f.contains(' ') || f.len() >= 4)
            .collect();
        tails.push(n.clone());
        let bare: Vec<String> = (tails.iter())
            .filter_map(|f| f.strip_suffix(" cell"))
            .filter(|f| f.len() >= 3)
            .map(str::to_string)
            .collect();
        for f in tails.into_iter().chain(bare) {
            if f != *n && !names.contains_key(&f) {
                short.entry(f).or_default().push(e);
            }
        }
    }
    let own = names.iter().map(|(n, e)| (n.clone(), (*e).clone()));
    let short = (short.into_iter())
        .filter(|(_, es)| es.len() == 1)
        .map(|(f, es)| (f, es[0].clone()));
    own.chain(short).collect()
}

/// A text ([`norm`]ed) with where it names its entities and where its
/// clauses and sentences part.
struct Text {
    t: String,
    spans: Vec<Span>,
    /// Clause breaks, then sentence ends, as positions.
    breaks: Vec<usize>,
    stops: Vec<usize>,
}

impl Text {
    fn new(text: &str, entities: &[String]) -> Self {
        let t = norm(text);
        let mut found: Vec<Span> = forms(entities)
            .into_iter()
            .flat_map(|(f, e)| {
                find(&t, &f, true)
                    .into_iter()
                    .map(move |(start, end)| Span {
                        start,
                        end,
                        entity: e.clone(),
                    })
            })
            .collect();
        // The longest of overlapping names.
        found.sort_by_key(|s| (s.start, std::cmp::Reverse(s.end)));
        let mut spans: Vec<Span> = Vec::new();
        for s in found {
            if spans.last().is_none_or(|l| s.start >= l.end) {
                spans.push(s);
            }
        }
        let mut text = Self {
            t,
            spans,
            breaks: Vec::new(),
            stops: Vec::new(),
        };
        text.stops = text.marks(&['.', '!', '?']);
        // `and` inside `between A and B` joins no clauses.
        let tied: Vec<usize> = text.betweens().iter().map(|b| b.3).collect();
        let mut breaks = text.marks(&[',', ';', ':']);
        for w in JOINS {
            breaks.extend(
                (find(&text.t, w, false).into_iter())
                    .map(|(i, _)| i)
                    .filter(|i| !text.inside(*i) && !tied.contains(i)),
            );
        }
        breaks.extend(&text.stops);
        breaks.sort_unstable();
        text.breaks = breaks;
        text
    }

    /// Whether `p` is inside a name.
    fn inside(&self, p: usize) -> bool {
        self.spans.iter().any(|s| s.start <= p && p < s.end)
    }

    /// Where the text has any of `chars` outside a name.
    fn marks(&self, chars: &[char]) -> Vec<usize> {
        (self.t.char_indices())
            .filter(|(i, c)| chars.contains(c) && !self.inside(*i))
            .map(|(i, _)| i)
            .collect()
    }

    /// The part around `p` between the nearest of `marks`.
    fn around(p: usize, marks: &[usize], len: usize) -> (usize, usize) {
        let a = marks.iter().rev().find(|&&m| m < p).map_or(0, |m| m + 1);
        let b = marks.iter().find(|&&m| m >= p).copied().unwrap_or(len);
        (a, b)
    }

    fn clause(&self, p: usize) -> (usize, usize) {
        Self::around(p, &self.breaks, self.t.len())
    }

    fn sentence(&self, p: usize) -> (usize, usize) {
        Self::around(p, &self.stops, self.t.len())
    }

    /// The entities named in `range`, in order.
    fn named(&self, (a, b): (usize, usize)) -> Vec<&str> {
        (self.spans.iter())
            .filter(|s| s.start >= a && s.end <= b)
            .map(|s| s.entity.as_str())
            .collect()
    }

    /// Who the text speaks of at `p`, among the order words at `said`: the
    /// last type `p`'s clause names before it; else, in a [`RELATIVE`]
    /// clause, the last its sentence names before it; else the sentence's
    /// subject, the last type it names before its first order word (`A
    /// comes after B and before C` says `before` of A).
    fn subject(&self, p: usize, said: &[usize]) -> Option<&str> {
        let (a, _) = self.clause(p);
        let (s, _) = self.sentence(p);
        let last = |from: usize, to: usize| self.named((from, to)).last().copied();
        let relative = a > 0 && self.words(&RELATIVE).iter().any(|&(q, _)| q + 1 == a);
        let lead = said.iter().copied().filter(|&q| q >= s).min().unwrap_or(p);
        last(a, p)
            .or_else(|| last(s, p).filter(|_| relative))
            .or_else(|| last(s, lead.min(p)))
    }

    /// The first named from `p` on in its sentence.
    fn after(&self, p: usize) -> Option<&Span> {
        let (_, b) = self.sentence(p);
        self.spans.iter().find(|s| s.start >= p && s.end <= b)
    }

    /// Where the text has any of `words`, with the word.
    fn words<'w>(&self, words: &[&'w str]) -> Vec<(usize, &'w str)> {
        let mut out: Vec<(usize, &str)> = (words.iter())
            .flat_map(|w| {
                find(&self.t, w, false)
                    .into_iter()
                    .map(move |(i, _)| (i, *w))
            })
            .collect();
        out.sort_unstable();
        out
    }

    /// Each `between A and B`: where `between` is, the two entities, and
    /// where the `and` is.
    fn betweens(&self) -> Vec<(usize, String, String, usize)> {
        let mut out = Vec::new();
        for (p, _) in find(&self.t, "between", false) {
            let Some(x) = self.after(p) else { continue };
            let Some(y) = self.after(x.end) else { continue };
            let gap = self.t[x.end..y.start].trim();
            let gap = gap.trim_end_matches(" the").trim_end_matches(" an");
            if gap == "and" {
                let and = x.end + self.t[x.end..].find("and").unwrap_or(0);
                out.push((p, x.entity.clone(), y.entity.clone(), and));
            }
        }
        out
    }
}

/// Whether `a` comes before `b` on a lineage.
fn earlier(f: &OrderFacts, a: &str, b: &str) -> bool {
    f.lineages.iter().any(|l| {
        let at = |t: &str| l.path.iter().position(|p| p == t);
        matches!((at(a), at(b)), (Some(i), Some(j)) if i < j)
    })
}

/// The verdicts of the edges between `a` and `b`, or of `a`'s edges.
fn verdicts<'f>(
    f: &'f OrderFacts,
    a: &'f str,
    b: Option<&'f str>,
) -> impl Iterator<Item = Verdict> + 'f {
    f.edges
        .iter()
        .filter(move |e| match b {
            Some(b) => e.is_pair(a, b),
            None => e.a.as_ref() == a || e.b.as_ref() == a,
        })
        .filter_map(|e| e.verdict)
}

/// The claims a trajectory text (its run paragraph, `subject` none, or a
/// type's line) makes of the order `f`.
pub fn order_claims(text: &str, f: &OrderFacts, subject: Option<&str>) -> Claim {
    let x = Text::new(text, &f.types());
    let betweens = x.betweens();
    let (forward, backward) = (x.words(&FORWARD), x.words(&BACKWARD));
    let mut said: Vec<usize> = betweens.iter().map(|b| b.0).collect();
    said.extend(forward.iter().chain(&backward).map(|(p, _)| *p));
    // Who a claim at `p` is about, else the subject: the same for every
    // kind of claim.
    let who = |p: usize| x.subject(p, &said).or(subject);
    for (p, a, b, _) in betweens {
        let Some(t) = who(p) else { continue };
        let fits = |u: &str, v: &str| earlier(f, u, t) && earlier(f, t, v);
        if !(fits(&a, &b) || fits(&b, &a)) {
            return Err(format!("{t} is not between {a} and {b}"));
        }
    }
    for (words, forward) in [(forward, true), (backward, false)] {
        for (p, w) in words {
            let (Some(l), Some(r)) = (who(p), x.after(p + w.len())) else {
                continue;
            };
            let (u, v) = if forward {
                (l, &*r.entity)
            } else {
                (&*r.entity, l)
            };
            if u != v && !earlier(f, u, v) {
                return Err(format!("{u} is not before {v}"));
            }
        }
    }
    for (p, _) in find(&x.t, "from", false) {
        let Some(a) = x.after(p) else { continue };
        let Some(b) = x.after(a.end) else { continue };
        let gap = x.t[a.end..b.start].trim();
        let to = ["to", "to the", "to a", "to an"].contains(&gap);
        let edge = f
            .edges
            .iter()
            .any(|e| *e.a == *a.entity && *e.b == *b.entity);
        if to && !earlier(f, &a.entity, &b.entity) && !edge {
            return Err(format!("{} does not lead to {}", a.entity, b.entity));
        }
    }
    let verdict_words = [
        Verdict::Supported,
        Verdict::Unsupported,
        Verdict::Candidate,
        Verdict::Inferred,
    ]
    .map(|v| (v.as_str(), v));
    let verdict_words = verdict_words
        .iter()
        .chain(&[("not supported", Verdict::Unsupported)]);
    for &(word, v) in verdict_words {
        for (p, _) in find(&x.t, word, false) {
            // `not supported` is read as one word.
            if word == "supported" && x.t[..p].ends_with("not ") {
                continue;
            }
            let named = x.named(x.clause(p));
            let (a, b) = match (named.as_slice(), who(p)) {
                ([a, b, ..], _) => (*a, Some(*b)),
                ([a], Some(s)) if *a != s => (s, Some(*a)),
                ([a], _) => (*a, None),
                ([], Some(s)) => (s, None),
                ([], None) => continue,
            };
            if !verdicts(f, a, b).any(|e| e == v) {
                let of = b.map_or_else(|| format!("{a}'s edges"), |b| format!("{a} – {b}"));
                return Err(format!("no {word} edge in {of}"));
            }
        }
    }
    let roles: [(&[&str], Vec<String>, &str); 3] = [
        (&ROOT, f.roots(), "a root"),
        (&BRANCH, f.branch_points(), "a branch point"),
        (&END, f.ends(), "an end"),
    ];
    for (words, holders, role) in &roles {
        for (p, _) in x.words(words) {
            let named = x.named(x.clause(p));
            let named = if named.is_empty() {
                who(p).into_iter().collect()
            } else {
                named
            };
            if !named.is_empty() && !named.iter().any(|t| holders.iter().any(|h| h == t)) {
                return Err(format!("{} is not {role}", named.join(" or ")));
            }
        }
    }
    Ok(())
}

/// Where `method` ranks `label` among cluster `c`'s evidence: first by its
/// top, else by the choice's rank.
fn rank_of(c: &ClusterInput, method: Method, label: &str) -> Option<usize> {
    if c.tops.iter().any(|(m, t)| *m == method && t == label) {
        return Some(0);
    }
    (c.choices.iter())
        .find(|ch| ch.label() == label)
        .and_then(|ch| ch.ranks[method.index()])
        .filter(|&r| r > 0)
}

/// Where `t` says a rank: each rank word's position and rank, 0 for a
/// first; a superlative after `second-` is that ordinal's.
fn ranks(t: &str) -> Vec<(usize, usize)> {
    let mut out: Vec<(usize, usize)> = Vec::new();
    for (i, w) in ORDINALS.iter().enumerate() {
        out.extend(find(t, w, false).into_iter().map(|(p, _)| (p, i + 1)));
    }
    for w in FIRST {
        out.extend(
            (find(t, w, false).into_iter())
                .filter(|(p, _)| !ORDINALS.iter().any(|o| t[..*p].ends_with(&format!("{o}-"))))
                .map(|(p, _)| (p, 0)),
        );
    }
    out.sort_unstable();
    out
}

/// The claims a cluster's sentence makes of cluster `c`'s methods: each
/// method named takes the nearest rank word of its clause, else, when it
/// is listed with the method before it (`the highest enrichment and SuSiE
/// ranking`), that one's; the claim is about the choice its clause names,
/// else `label` (the one the sentence proposes).
pub fn rank_claims(text: &str, c: &ClusterInput, label: &str) -> Claim {
    let mut entities: Vec<String> = c.choices.iter().map(|ch| ch.label().to_string()).collect();
    entities.extend(c.calls.iter().map(|k| k.label.clone()));
    entities.extend(c.tops.iter().map(|(_, t)| t.clone()));
    entities.sort();
    entities.dedup();
    // `enrichment score` is the enrichment, not the projection score.
    let text = norm(text)
        .replace("enrichment score", "enrichment")
        .replace("susie score", "susie");
    let x = Text::new(&text, &entities);
    let mut named: Vec<(usize, usize, Method)> = Vec::new();
    for m in Method::ALL {
        named.extend(
            find(&x.t, m.name(), false)
                .into_iter()
                .map(|(p, e)| (p, e, m)),
        );
        for all in [
            "all methods",
            "every method",
            "each method",
            "all three methods",
        ] {
            named.extend(find(&x.t, all, false).into_iter().map(|(p, e)| (p, e, m)));
        }
    }
    named.sort_by_key(|&(p, e, m)| (p, e, m.index()));
    let said = ranks(&x.t);
    // Words a list of methods may have between them.
    let listed = |gap: &str| {
        gap.split(|ch: char| ch == ',' || ch.is_whitespace())
            .all(|w| ["", "and", "the", "ranking", "rankings", "'s"].contains(&w))
    };
    let mut taken: Vec<Option<usize>> = Vec::new();
    for (i, &(p, _, m)) in named.iter().enumerate() {
        let clause = x.clause(p);
        let own = (said.iter())
            .filter(|(q, _)| (clause.0..clause.1).contains(q))
            .min_by_key(|(q, _)| q.abs_diff(p))
            .map(|&(_, r)| r);
        let rank = own.or_else(|| {
            // The methods of `all methods` share its place, and so its rank.
            let (q, prev_end, _) = *named.get(i.checked_sub(1)?)?;
            let joined = q == p || x.t.get(prev_end..p).is_some_and(listed);
            joined.then(|| taken[i - 1]).flatten()
        });
        taken.push(rank);
        let Some(rank) = rank else { continue };
        let about = x.named(clause);
        let about = if about.is_empty() { vec![label] } else { about };
        if !about.iter().any(|t| rank_of(c, m, t) == Some(rank)) {
            let place = rank.checked_sub(1).map_or("first", |r| ORDINALS[r]);
            return Err(format!(
                "{} does not rank {} {place}",
                m.name(),
                about.join(" or ")
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "tests/claims.rs"]
mod tests;
