//! [`super`] on a made-up cache, never online.

use super::*;

const REV: &str = "0123456789abcdef0123456789abcdef01234567";

/// A cache holding `files` of `org/model` at [`REV`].
fn cache_with(files: &[&str]) -> (tempfile::TempDir, Cache) {
    let root = tempfile::tempdir().unwrap();
    let repo = root.path().join("models--org--model");
    std::fs::create_dir_all(repo.join("refs")).unwrap();
    std::fs::write(repo.join("refs").join(REV), REV).unwrap();
    let snap = repo.join("snapshots").join(REV);
    std::fs::create_dir_all(&snap).unwrap();
    for f in files {
        std::fs::write(snap.join(f), "x").unwrap();
    }
    let cache = Cache::new(root.path().to_path_buf());
    (root, cache)
}

#[test]
fn either_switch_rules_out_a_download() {
    let os = |s: &str| Some(std::ffi::OsString::from(s));
    let s = |s: &str| Some(s.to_string());
    assert_eq!(offline_switch(None, None), None);
    assert_eq!(offline_switch(os("1"), None), Some(OFFLINE_ENV));
    assert_eq!(offline_switch(os(""), s("1")), Some(OFFLINE_ENV));
    for on in ["1", "true", "YES", "On"] {
        assert_eq!(offline_switch(None, s(on)), Some(HF_OFFLINE_ENV), "{on}");
    }
    for off in ["0", "false", ""] {
        assert_eq!(offline_switch(None, s(off)), None, "{off}");
    }
}

#[test]
fn cached_files_are_read_offline() {
    let (_root, cache) = cache_with(&["config.json", "model.safetensors"]);
    let m = Model::in_cache("org/model", REV, cache, Some(HF_OFFLINE_ENV));
    let got = m
        .files(&["config.json", "model.safetensors"], false)
        .unwrap();
    assert!(got[1].ends_with("model.safetensors"));
    assert!(m.cached("tokenizer.json").is_none());
}

#[test]
fn offline_a_file_the_cache_lacks_is_an_error_naming_the_switch() {
    let (root, cache) = cache_with(&["config.json"]);
    let m = Model::in_cache("org/model", REV, cache, Some(HF_OFFLINE_ENV));
    let e = m
        .files(&["config.json", "model.safetensors"], false)
        .unwrap_err()
        .to_string();
    assert!(e.contains("org/model") && e.contains(HF_OFFLINE_ENV), "{e}");
    assert!(e.contains(&root.path().display().to_string()), "{e}");
}
