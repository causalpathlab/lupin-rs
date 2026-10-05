//! Marker panels for lupin, as data: each panel is a `gene<TAB>cell type`
//! table, its Cell Ontology sidecar (`label<TAB>CL:id<TAB>note`) and a
//! README naming its source and licence, compiled in from `data/<name>/`.
//! lupin reads them with `--markers panel:<name>`.

/// One panel's files.
#[derive(Debug, Clone, Copy)]
pub struct Panel {
    pub name: &'static str,
    /// The marker table's file name, as lupin writes it out.
    pub table_file: &'static str,
    pub table: &'static [u8],
    /// The Cell Ontology sidecar, when the panel has one.
    pub sidecar: Option<&'static [u8]>,
    pub readme: &'static [u8],
}

/// This crate's version, for a run's record of the panel it used.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

static PANELS: &[Panel] = &[Panel {
    name: "example",
    table_file: "example.tsv",
    table: include_bytes!("../data/example/example.tsv"),
    sidecar: Some(include_bytes!("../data/example/example.cl.tsv")),
    readme: include_bytes!("../data/example/README.md"),
}];

static NAMES: &[&str] = &["example"];

/// The panels' names.
#[must_use]
pub fn names() -> &'static [&'static str] {
    NAMES
}

/// The panel called `name`.
#[must_use]
pub fn get(name: &str) -> Option<Panel> {
    PANELS.iter().find(|p| p.name == name).copied()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_name_has_a_panel_with_a_table_and_a_readme() {
        assert_eq!(names().len(), PANELS.len());
        for name in names() {
            let p = get(name).unwrap();
            assert!(!p.table.is_empty() && !p.readme.is_empty(), "{name}");
            assert!(p.table_file.starts_with(name), "{name}");
        }
        assert!(get("no such panel").is_none());
    }
}
