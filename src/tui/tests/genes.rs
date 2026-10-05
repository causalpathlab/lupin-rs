//! [`super::GeneFilter`]: patterns, files and suggestions.

use super::*;

#[test]
fn patterns_match_names_and_star_runs() {
    let f = GeneFilter::parse("# family\nPA-*\nPB*\n*-SX1\nGENE9\nPC*L1*\n\n");
    for g in ["PA-GENE3", "PB13A", "GENE4-SX1", "GENE9", "PC11-ML1Q"] {
        assert!(f.hides(g), "{g}");
    }
    for g in ["PAOR", "PD6", "GENE92", "SX1-X", "GENE5"] {
        assert!(!f.hides(g), "{g}");
    }
    assert_eq!(f.len(), 5);
}

#[test]
fn a_suggestion_is_the_stem_with_a_star() {
    assert_eq!(suggest("PA-GENE3"), "PA-*");
    assert_eq!(suggest("PC11-111M22.2"), "PC11-*");
    assert_eq!(suggest("PE-2090I13.1"), "PE-*");
    assert_eq!(suggest("GENE5"), "GENE5");
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
        panel: None,
    };
    std::fs::create_dir_all(root.path().join("user")).unwrap();
    std::fs::write(root.path().join("user").join(HIDDEN), "PB*\n").unwrap();

    let mut f = GeneFilter::load(&search).unwrap();
    assert!(f.hides("PB13A") && !f.hides("PA-GENE3"));
    let file = GeneFilter::file(&search).unwrap();
    assert!(
        file.starts_with(root.path().join("run")),
        "the project's file"
    );
    f.add("PA-*", &file).unwrap();
    f.add("PA-*", &file).unwrap();
    assert!(f.hides("PA-GENE3"));

    let again = GeneFilter::load(&search).unwrap();
    assert!(again.hides("PA-GENE3") && again.hides("PB13A"));
    assert_eq!(again.len(), 2, "no duplicate line");
}
