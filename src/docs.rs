//! `lupin docs`: the method write-ups, compiled into the binary.
//!
//! `include_str!`, not paths read at runtime: the binary is often the only
//! thing on the machine that ran the analysis, and a doc you cannot reach from
//! there is a doc nobody reads. It also breaks the build if a file is moved or
//! deleted, which keeps the list honest about what exists.

use anyhow::Result;
use clap::builder::PossibleValue;
use clap::{Args, ValueEnum};

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum Topic {
    /// Marker cell-type annotation by projection, end to end.
    Annotation,
    /// Why the annotation pools cells into coarse clusters.
    Grouping,
    /// Annotation on a gene-annotated Cell Ontology (partly implemented).
    OntologyPlan,
    /// Supervised trajectory: an explicit prior, checked against the data (not implemented).
    TrajectoryPlan,
}

/// Every write-up with its one-line blurb, so the listing can never advertise
/// a topic the command cannot print.
const DOCS: &[(Topic, &str, &str)] = &[
    (
        Topic::Annotation,
        "METHOD  marker cell-type annotation by projection, end to end",
        include_str!("../docs/annotation-methods.md"),
    ),
    (
        Topic::Grouping,
        "METHOD  why the annotation pools cells into coarse clusters",
        include_str!("../docs/annotation-grouping.md"),
    ),
    (
        Topic::OntologyPlan,
        "PLAN    (partly implemented) annotation on a gene-annotated Cell Ontology",
        include_str!("../docs/annotation-ontology-plan.md"),
    ),
    (
        Topic::TrajectoryPlan,
        "PLAN    (not implemented) supervised trajectory: prior, check, diffusion pseudotime",
        include_str!("../docs/trajectory-plan.md"),
    ),
];

#[derive(Args, Debug)]
pub struct DocsArgs {
    #[arg(value_enum, help = "Which write-up to print (omit to list them)")]
    pub topic: Option<Topic>,
}

pub fn run_docs(args: &DocsArgs) -> Result<()> {
    let Some(want) = args.topic else {
        println!("lupin method write-ups (`lupin docs <TOPIC>` to read one):\n");
        // The name clap accepts, e.g. `ontology-plan`, not the Debug form.
        let slugs: Vec<String> = DOCS
            .iter()
            .map(|(topic, _, _)| {
                topic
                    .to_possible_value()
                    .as_ref()
                    .map(PossibleValue::get_name)
                    .unwrap_or_default()
                    .to_string()
            })
            .collect();
        let width = slugs.iter().map(String::len).max().unwrap_or(0);
        for (slug, (_, blurb, _)) in slugs.iter().zip(DOCS) {
            println!("  {slug:<width$} {blurb}");
        }
        return Ok(());
    };
    let text = DOCS
        .iter()
        .find(|(t, _, _)| *t == want)
        .map(|(_, _, text)| *text)
        .expect("every Topic has a row in DOCS");
    println!("{text}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_topic_has_a_write_up() {
        for t in Topic::value_variants() {
            assert!(
                DOCS.iter().any(|(d, _, text)| d == t && !text.is_empty()),
                "{t:?}"
            );
        }
    }
}
