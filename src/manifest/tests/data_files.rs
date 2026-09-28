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
    let d = ClData::load(search(root.path()), None, None).unwrap();
    assert_eq!(d.rules, MatchRules::default());
    assert!(d.aliases.is_empty());
    assert!(d.ontology.is_none());
    assert!(d.sources.is_empty());
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

    let d = ClData::load(s, None, Some(&run.to_string_lossy())).unwrap();
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
    assert_eq!(d.sources.len(), 6, "{:?}", d.sources);
    assert!(d.ontology.as_ref().unwrap().ends_with(ONTOLOGY));
    assert_eq!(
        d.search.amend_aliases().unwrap(),
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
    let d = ClData::load(search(root.path()), None, None).unwrap();
    assert_eq!(d.aliases.get("NK"), Some("CL:1"));
}

#[test]
fn this_releases_files_are_published_at_its_tag() {
    let url = release_url(ALIASES);
    assert!(url.starts_with("https://raw.githubusercontent.com/causalpathlab/lupin-rs/v"));
    assert!(url.ends_with("/data/cl_aliases.tsv"));
}
