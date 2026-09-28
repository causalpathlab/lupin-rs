//! Remembering an alias picked in the TUI.

use super::*;
use crate::annotate::cl_rules::Aliases;

#[test]
fn a_remembered_alias_is_appended_and_reads_back_as_the_top_layer() {
    let root = tempfile::tempdir().unwrap();
    let file = root.path().join("run").join("lupin").join("cl_aliases.tsv");
    remember_alias(&file, "EMP", "CL:0000049", "picked\tin the TUI").unwrap();
    remember_alias(&file, "BaEoMa", "CL:0000767", "second").unwrap();
    let text = std::fs::read_to_string(&file).unwrap();
    assert!(text.starts_with('#'), "a new file starts with its header");
    assert_eq!(text.matches("label\tcl_id\tnote").count(), 1, "one header");
    assert!(
        text.contains("EMP\tCL:0000049\tpicked in the TUI\n"),
        "tabs in the note are flattened"
    );
    let mut a = Aliases::default();
    assert_eq!(a.add_tsv(&text, "project"), 2);
    assert_eq!(a.get("emp"), Some("CL:0000049"));
}
