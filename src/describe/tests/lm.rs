//! What a generation may write, on a made-up vocabulary.

use super::*;
use std::str::FromStr;

/// The two stop tokens (special), a word, a markup sign, a name with a
/// digit and a digit.
fn tokenizer() -> Tokenizer {
    let special = |id: u32, t: &str| {
        serde_json::json!({"id": id, "content": t, "single_word": false, "lstrip": false,
            "rstrip": false, "normalized": false, "special": true})
    };
    let json = serde_json::json!({
        "version": "1.0",
        "truncation": null,
        "padding": null,
        "added_tokens": [special(0, "<|im_end|>"), special(1, "<|endoftext|>")],
        "normalizer": null,
        "pre_tokenizer": {"type": "WhitespaceSplit"},
        "post_processor": null,
        "decoder": null,
        "model": {
            "type": "WordLevel",
            "vocab": {"<|im_end|>": 0, "<|endoftext|>": 1, "cell": 2, "<": 3, "CD8": 4, "7": 5, "[UNK]": 6},
            "unk_token": "[UNK]"
        }
    });
    Tokenizer::from_str(&json.to_string()).unwrap()
}

#[test]
fn a_generation_without_numbers_or_markup_can_still_end_its_turn() {
    let t = tokenizer();
    let stops = stop_ids(&t);
    assert_eq!(stops, [0, 1]);
    let bans = Bans {
        words: vec!["cell".into()],
        numbers: true,
        keep: vec!["CD8".into()],
    };
    let ban = banned(&t, &numeric_ids(&t), &stops, &bans).unwrap();
    assert!(
        ban.contains(&3) && ban.contains(&5),
        "markup and digits: {ban:?}"
    );
    assert!(ban.contains(&2), "the word banned");
    assert!(!ban.contains(&4), "a name kept keeps its digit");
    let bans = Bans {
        keep: vec!["cell".into()],
        ..bans
    };
    let ban = banned(&t, &numeric_ids(&t), &stops, &bans).unwrap();
    assert!(!ban.contains(&2), "a name kept keeps its words too");
    assert!(
        !ban.contains(&0) && !ban.contains(&1),
        "the stop tokens are never banned, though their text has `<`, `|` and `>`: {ban:?}"
    );
}
