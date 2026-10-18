//! `{out}.trajectory_lineages.tsv`: one root-to-leaf path per line, with
//! its name and component. The run writes it; `describe` and the order
//! view read it.

use anyhow::{Context, Result};

/// One root-to-leaf path.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Lineage {
    pub name: String,
    pub component: String,
    pub path: Vec<String>,
}

/// Between the types of a lineage's path.
const STEP: &str = " > ";

pub(crate) fn write(path: &str, lineages: &[Lineage]) -> Result<()> {
    let mut text = String::from("# lineage\tcomponent\tpath\n");
    for l in lineages {
        text += &format!("{}\t{}\t{}\n", l.name, l.component, l.path.join(STEP));
    }
    std::fs::write(path, text).with_context(|| format!("writing {path}"))
}

pub(crate) fn read(path: &str) -> Result<Vec<Lineage>> {
    Ok(parse(
        &std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?,
    ))
}

/// The lineages of the file's `text`.
#[must_use]
pub(crate) fn parse(text: &str) -> Vec<Lineage> {
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let mut it = l.split('\t');
            let (name, component, path) = (it.next()?, it.next()?, it.next()?);
            Some(Lineage {
                name: name.to_string(),
                component: component.to_string(),
                path: path.split(STEP).map(str::to_string).collect(),
            })
        })
        .collect()
}
