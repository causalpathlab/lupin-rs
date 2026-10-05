use super::*;

#[test]
fn a_bundled_panel_is_written_out_with_its_sidecar() {
    let tmp = tempfile::tempdir().unwrap();
    let root = tmp.path().join("panels");
    let (table, b) = write_out("example", &root).unwrap();
    assert_eq!(b.name, "example");
    assert_eq!(b.version, lupin_panels::VERSION);
    let panel = crate::annotate::markers::read_panel(&table.to_string_lossy()).unwrap();
    assert!(panel.iter().any(|(g, t)| g == "GENE1" && t == "CT1"));
    let side = crate::manifest::data_files::sidecar_of(&table).unwrap();
    assert!(side.is_file(), "the sidecar sits beside the table");
    assert_eq!(bundled_of(&table), Some(b.clone()), "the run records it");
    assert_eq!(bundled_of(Path::new("d/x.tsv")), None);
    // Written again, nothing changes.
    assert_eq!(write_out("example", &root).unwrap().0, table);
}

#[test]
fn an_unknown_name_lists_the_bundled_ones() {
    let root = tempfile::tempdir().unwrap();
    let e = write_out("no-such", root.path()).unwrap_err().to_string();
    assert!(
        e.contains("no bundled panel no-such") && e.contains("example"),
        "{e}"
    );
}

#[test]
fn a_path_is_left_as_it_is() {
    assert_eq!(resolve("x.tsv").unwrap(), ("x.tsv".to_string(), None));
}

#[test]
fn the_bundled_sidecar_maps_the_panels_labels() {
    let root = tempfile::tempdir().unwrap();
    let (table, _) = write_out("example", root.path()).unwrap();
    let search =
        crate::manifest::data_files::SearchPath::new(None).with_panel(&table.to_string_lossy());
    let d = crate::manifest::data_files::ClData::load(
        search,
        None,
        None,
        crate::manifest::data_files::Fetch::Never,
    )
    .unwrap();
    assert_eq!(d.aliases.get("CT1"), Some("CL:9000001"));
    assert_eq!(d.aliases.get("CT3"), None);
}
