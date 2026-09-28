//! [`super::GeneFilter`]: patterns, files and suggestions.

use super::*;

#[test]
fn patterns_match_names_and_star_runs() {
    let f = GeneFilter::parse("# mito\nMT-*\nRPL*\n*-AS1\nXIST\nRP*L1*\n\n");
    for g in ["MT-ND3", "RPL13A", "EPB41L4A-AS1", "XIST", "RP11-ML1Q"] {
        assert!(f.hides(g), "{g}");
    }
    for g in ["MTOR", "RPS6", "XIST2", "AS1-X", "CRHBP"] {
        assert!(!f.hides(g), "{g}");
    }
    assert_eq!(f.len(), 5);
}

#[test]
fn a_suggestion_is_the_stem_with_a_star() {
    assert_eq!(suggest("MT-ND3"), "MT-*");
    assert_eq!(suggest("RP11-111M22.2"), "RP11-*");
    assert_eq!(suggest("CTD-2090I13.1"), "CTD-*");
    assert_eq!(suggest("CRHBP"), "CRHBP");
}

#[test]
fn added_patterns_persist_and_load_with_the_users() {
    let root = tempfile::tempdir().unwrap();
    let search = SearchPath {
        install: None,
        source: None,
        cache: None,
        user: Some(root.path().join("user")),
        project: Some(root.path().join("run").join("lupin")),
    };
    std::fs::create_dir_all(root.path().join("user")).unwrap();
    std::fs::write(root.path().join("user").join(HIDDEN), "RPL*\n").unwrap();

    let mut f = GeneFilter::load(&search).unwrap();
    assert!(f.hides("RPL13A") && !f.hides("MT-ND3"));
    let file = GeneFilter::file(&search).unwrap();
    assert!(
        file.starts_with(root.path().join("run")),
        "the project's file"
    );
    f.add("MT-*", &file).unwrap();
    f.add("MT-*", &file).unwrap();
    assert!(f.hides("MT-ND3"));

    let again = GeneFilter::load(&search).unwrap();
    assert!(again.hides("MT-ND3") && again.hides("RPL13A"));
    assert_eq!(again.len(), 2, "no duplicate line");
}
