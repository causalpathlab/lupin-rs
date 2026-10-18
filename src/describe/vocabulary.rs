//! The words a model's biology sentence may not use: those that state an
//! order, a role, a verdict, a rank or a method's reading. Code writes every
//! such fact from the evidence ([`super::opinion::facts`],
//! [`super::trajectory`]); a sentence by the model that uses one is asked
//! again without it, then dropped. Words are found whole, as
//! [`super::mentions`] finds names.

/// Words putting one type before another.
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
/// Words putting it after.
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
/// Words for a rank.
const RANK: [&str; 23] = [
    "first",
    "second",
    "third",
    "fourth",
    "fifth",
    "sixth",
    "seventh",
    "eighth",
    "ninth",
    "tenth",
    "last",
    "top",
    "highest",
    "strongest",
    "best",
    "leading",
    "greatest",
    "likeliest",
    "rank",
    "ranks",
    "ranked",
    "ranking",
    "rankings",
];
/// An edge's verdicts ([`crate::trajectory::edges::Verdict`]).
const VERDICT: [&str; 4] = ["supported", "unsupported", "candidate", "inferred"];
/// The methods and their readings.
const METHOD: [&str; 12] = [
    "enrichment",
    "susie",
    "score",
    "scores",
    "scored",
    "pip",
    "posterior",
    "probability",
    "method",
    "methods",
    "significant",
    "significance",
];

/// Every word and phrase of the vocabulary.
fn all() -> impl Iterator<Item = &'static str> {
    FORWARD
        .into_iter()
        .chain(BACKWARD)
        .chain(ROOT)
        .chain(BRANCH)
        .chain(END)
        .chain(RANK)
        .chain(VERDICT)
        .chain(METHOD)
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

/// The first word of the vocabulary `text` (lower case) uses as a whole
/// word, in the vocabulary's order.
#[must_use]
pub fn used(text: &str) -> Option<&'static str> {
    all().find(|w| !find(text, w, false).is_empty())
}

/// The words a generation is kept from starting, each phrase by its first
/// word, in lower case and capitalised, but those `names` (lower case)
/// have as a whole word: a name the sentence may cite keeps its words.
#[must_use]
pub fn bans(names: &[String]) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for w in all().filter_map(|w| w.split(' ').next()) {
        if names.iter().any(|n| !find(n, w, false).is_empty()) {
            continue;
        }
        let mut c = w.chars();
        let cap: String = c
            .next()
            .map(char::to_uppercase)
            .into_iter()
            .flatten()
            .chain(c)
            .collect();
        for form in [w.to_string(), cap] {
            if !out.contains(&form) {
                out.push(form);
            }
        }
    }
    out.push("SuSiE".into());
    out
}

#[cfg(test)]
#[path = "tests/vocabulary.rs"]
mod tests;
