//! Layering data files along a search path, offline.

use super::*;

fn offline() {
    std::env::set_var(OFFLINE_ENV, "1");
}

fn search(root: &Path) -> SearchPath {
    SearchPath {
        install: Some(root.join("share")),
        source: None,
        cache: Some(root.join("cache")),
        user: Some(root.join("user")),
        project: Some(root.join("run").join("lupin")),
    }
}

fn write(p: PathBuf, text: &str) {
    fs::create_dir_all(p.parent().unwrap()).unwrap();
    fs::write(p, text).unwrap();
}

#[test]
fn nothing_found_means_literal_rules_and_no_aliases() {
    offline();
    let root = tempfile::tempdir().unwrap();
    let d = ClData::load(search(root.path()), None, None, Fetch::Never).unwrap();
    assert_eq!(d.rules, MatchRules::default());
    assert!(d.aliases.is_empty());
    assert!(d.ontology.is_none());
    assert!(d.rule_files.is_empty() && d.alias_files.is_empty());
}

#[test]
fn later_layers_win_and_every_file_read_is_recorded() {
    offline();
    let root = tempfile::tempdir().unwrap();
    let s = search(root.path());
    write(
        root.path().join("share").join(RULES),
        r#"{"any_word_order": true, "synonyms": {"EXACT": []}}"#,
    );
    write(
        root.path().join("user").join(RULES),
        r#"{"any_word_order": false}"#,
    );
    write(
        root.path().join("share").join(ALIASES),
        "NK\tCL:1\nMono\tCL:2\n",
    );
    write(root.path().join("run/lupin").join(ALIASES), "NK\tCL:9\n");
    let run = root.path().join("run.tsv");
    write(run.clone(), "Mono\tCL:8\n");
    write(
        root.path().join("cache").join(ONTOLOGY),
        "[Term]\nid: CL:1\nname: x\n",
    );

    let d = ClData::load(s, None, Some(&run.to_string_lossy()), Fetch::Never).unwrap();
    assert!(
        !d.rules.any_word_order,
        "the user layer overrides the install"
    );
    assert!(
        d.rules.counts("EXACT", &[]),
        "keys the user did not set are kept"
    );
    assert_eq!(
        d.aliases.get("nk"),
        Some("CL:9"),
        "the project beats the install"
    );
    assert_eq!(
        d.aliases.get("mono"),
        Some("CL:8"),
        "the run's file beats all"
    );
    assert_eq!(d.rule_files.len(), 2);
    assert_eq!(d.alias_files.len(), 3, "install, project and the run's own");
    assert!(d.ontology.as_ref().unwrap().ends_with(ONTOLOGY));
    assert_eq!(
        d.search.amend(ALIASES).unwrap(),
        root.path().join("run/lupin").join(ALIASES)
    );
}

#[test]
fn the_install_copy_is_used_before_the_cache() {
    offline();
    let root = tempfile::tempdir().unwrap();
    write(root.path().join("share").join(ALIASES), "NK\tCL:1\n");
    let s = search(root.path());
    write(s.cached(ALIASES).unwrap(), "NK\tCL:2\n");
    let d = ClData::load(search(root.path()), None, None, Fetch::Never).unwrap();
    assert_eq!(d.aliases.get("NK"), Some("CL:1"));
}

#[test]
fn this_releases_files_are_published_at_its_tag() {
    let url = release_url(ALIASES);
    assert!(url.starts_with("https://raw.githubusercontent.com/causalpathlab/lupin-rs/v"));
    assert!(url.ends_with("/data/cl_aliases.tsv"));
}

#[test]
fn a_rescore_reads_exactly_the_files_its_pass_recorded() {
    offline();
    let root = tempfile::tempdir().unwrap();
    write(
        root.path().join("share").join(RULES),
        r#"{"any_word_order": true}"#,
    );
    write(root.path().join("share").join(ALIASES), "NK\tCL:1\n");
    write(
        root.path().join("cache").join(ONTOLOGY),
        "[Term]\nid: CL:1\nname: x\n",
    );
    let pass = ClData::load(search(root.path()), None, None, Fetch::Never).unwrap();
    let record = pass.record(Some("test"));

    // Later, a user file appears: the rescore must not see it.
    write(root.path().join("user").join(ALIASES), "NK\tCL:7\n");
    let again = ClData::from_record(&record, search(root.path()))
        .unwrap()
        .unwrap();
    assert_eq!(again.aliases.get("NK"), Some("CL:1"));
    assert!(again.rules.any_word_order);
    assert!(again.ontology.as_ref().unwrap().is_absolute());

    // A recorded file gone: no record to trust.
    fs::remove_file(root.path().join("share").join(ALIASES)).unwrap();
    assert!(ClData::from_record(&record, search(root.path()))
        .unwrap()
        .is_none());
    assert!(
        ClData::from_record(&serde_json::json!({}), search(root.path()))
            .unwrap()
            .is_none()
    );
}

#[test]
fn a_run_file_is_found_as_given_else_beside_the_run() {
    let root = tempfile::tempdir().unwrap();
    let s = search(root.path());
    write(root.path().join("run").join("map.tsv"), "NK\tCL:1\n");
    assert_eq!(
        s.run_file("map.tsv").unwrap(),
        root.path().join("run").join("map.tsv")
    );
    assert!(s.run_file("absent.tsv").is_none());
}

#[test]
fn appended_lines_come_after_a_comment_header() {
    let root = tempfile::tempdir().unwrap();
    let f = root.path().join("a").join("x.tsv");
    append_line(&f, "what\ncolumns", "one").unwrap();
    append_line(&f, "what\ncolumns", "two").unwrap();
    assert_eq!(
        fs::read_to_string(&f).unwrap(),
        "# what\n# columns\none\ntwo\n"
    );
}

#[test]
fn the_gene_ontology_is_found_beside_the_run_and_never_fetched_offline() {
    offline();
    let root = tempfile::tempdir().unwrap();
    let run = root.path().join("run");
    let err = go_ontology(Some(&run)).unwrap_err().to_string();
    assert!(err.contains("--go-obo"), "{err}");
    let obo = run.join("lupin").join(GO_ONTOLOGY);
    write(obo.clone(), "[Term]\nid: GO:0000001\n");
    assert_eq!(go_ontology(Some(&run)).unwrap(), obo);
}

#[test]
fn every_ontology_is_cached_across_releases() {
    let root = tempfile::tempdir().unwrap();
    let s = search(root.path());
    let cache = root.path().join("cache");
    assert_eq!(s.cached(GO_ONTOLOGY).unwrap(), cache.join(GO_ONTOLOGY));
    assert_eq!(s.cached(ONTOLOGY).unwrap(), cache.join(ONTOLOGY));
    assert!(s.cached(RULES).unwrap().starts_with(cache.join("data")));
}

#[test]
fn go_annotations_are_found_beside_the_run_and_cached_across_releases() {
    use crate::annotate::go_signature::Species;
    offline();
    let root = tempfile::tempdir().unwrap();
    let run = root.path().join("run");
    let err = go_annotations(Species::Human, Some(&run))
        .unwrap_err()
        .to_string();
    assert!(err.contains("--gaf"), "{err}");
    let gaf = run.join("lupin").join(Species::Human.gaf_file());
    write(gaf.clone(), "!gaf-version: 2.2\n");
    assert_eq!(go_annotations(Species::Human, Some(&run)).unwrap(), gaf);
    let s = search(root.path());
    assert_eq!(
        s.cached(Species::Mouse.gaf_file()).unwrap(),
        root.path().join("cache").join(Species::Mouse.gaf_file())
    );
}
