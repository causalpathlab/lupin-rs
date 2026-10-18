//! The vocabulary a biology sentence may not use.

use super::*;

#[test]
fn a_word_is_found_whole_and_a_phrase_as_its_words() {
    assert_eq!(used("its markers rise before the end."), Some("before"));
    assert_eq!(used("ct a gives rise to ct b."), Some("gives rise to"));
    assert_eq!(used("ct a ranks it second."), Some("second"));
    assert_eq!(used("ct a, as susie says."), Some("susie"));
    assert_eq!(used("a supported edge."), Some("supported"));
    assert_eq!(
        used("endothelial genes, a basophil, scoring nothing."),
        None
    );
    assert_eq!(used("ct a expresses gene1 and gene2."), None);
}

#[test]
fn the_bans_start_each_word_in_both_cases_but_a_names_own() {
    let b = bans(&[]);
    for w in [
        "before",
        "Before",
        "gives",
        "Gives",
        "second",
        "SuSiE",
        "end",
        "enrichment",
    ] {
        assert!(b.contains(&w.to_string()), "{w}: {b:?}");
    }
    assert_eq!(
        b.iter().filter(|w| *w == "gives").count(),
        1,
        "a word once, whatever phrases start with it"
    );
    let b = bans(&["t cell receptor signaling, first step".to_string()]);
    assert!(!b.contains(&"first".to_string()), "a name's own word");
    assert!(b.contains(&"second".to_string()));
}
