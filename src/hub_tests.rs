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
    let got = m.files(["config.json", "model.safetensors"]).unwrap();
    assert!(got[1].ends_with("model.safetensors"));
    assert!(m.cached("tokenizer.json").is_none());
}

#[test]
fn offline_a_file_the_cache_lacks_is_an_error_naming_the_switch() {
    let (root, cache) = cache_with(&["config.json"]);
    let m = Model::in_cache("org/model", REV, cache, Some(HF_OFFLINE_ENV));
    let e = m
        .files(["config.json", "model.safetensors"])
        .unwrap_err()
        .to_string();
    assert!(e.contains("org/model") && e.contains(HF_OFFLINE_ENV), "{e}");
    assert!(e.contains(&root.path().display().to_string()), "{e}");
}

#[test]
fn a_commits_snapshot_is_read_without_a_ref_to_it() {
    // The Hub's Python tools write no `refs/<commit>` for a commit revision.
    let (root, cache) = cache_with(&["config.json"]);
    let repo = root.path().join("models--org--model");
    std::fs::remove_file(repo.join("refs").join(REV)).unwrap();
    let m = Model::in_cache("org/model", REV, cache.clone(), Some(HF_OFFLINE_ENV));
    assert!(m
        .cached("config.json")
        .is_some_and(|p| p.ends_with("config.json")));
    assert!(m.cached("tokenizer.json").is_none());
    // A branch is no commit: only its ref says which snapshot it is.
    let main = repo.join("snapshots").join("main");
    std::fs::create_dir_all(&main).unwrap();
    std::fs::write(main.join("config.json"), "x").unwrap();
    let m = Model::in_cache("org/model", "main", cache, Some(HF_OFFLINE_ENV));
    assert!(m.cached("config.json").is_none());
}

#[test]
fn files_fetched_into_the_default_cache_before_are_found_when_hf_home_names_another() {
    // The cache HF_HOME names lacks the model; the default one, where an
    // earlier lupin fetched it whatever HF_HOME said, has it.
    let (_now, empty) = cache_with(&[]);
    let (_before, default) = cache_with(&["config.json"]);
    let m = Model::in_cache("org/model", REV, empty, Some(HF_OFFLINE_ENV)).reading_too(default);
    let [got] = m.files(["config.json"]).unwrap();
    assert!(got.ends_with("config.json"));
}
