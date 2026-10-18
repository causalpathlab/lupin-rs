//! `{out}.trajectory_lineages.tsv`: one root-to-leaf path per line, with
//! its name and component. The run writes it; `describe` and the order
//! view read it.

use anyhow::{Context, Result};
use legume_numeric::matrix::common_io::{read_lines_of_words_delim, ReadLinesOut};

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
    let ReadLinesOut { lines, .. } =
        read_lines_of_words_delim(path, "\t", -1).with_context(|| format!("reading {path}"))?;
    Ok(lines
        .into_iter()
        .filter(|w| {
            w.first()
                .is_some_and(|f| !f.is_empty() && !f.starts_with('#'))
        })
        .filter_map(|w| {
            let [name, component, path, ..] = &w[..] else {
                return None;
            };
            Some(Lineage {
                name: name.to_string(),
                component: component.to_string(),
                path: path.split(STEP).map(str::to_string).collect(),
            })
        })
        .collect())
}
