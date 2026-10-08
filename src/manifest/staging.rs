//! Passes write under a staging prefix beside their own and are promoted
//! onto it when they finish, so a stopped or failed pass leaves the round
//! it would have replaced as it was.

use super::run::{annotated_path, parent_dir};
use anyhow::{Context, Result};
use serde_json::Value;
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

/// The prefix a pass to `out` writes under until it is promoted: beside
/// it, a dotfile, so no listing of runs or files shows a pass in progress.
#[must_use]
pub fn staging_prefix(out: &str) -> String {
    let (_, base) = split(out);
    // The directory as `out` writes it, so the prefix records the same way.
    match Path::new(out).parent().map(Path::to_string_lossy) {
        Some(dir) if !dir.is_empty() => format!("{dir}/.{base}.staging"),
        _ => format!(".{base}.staging"),
    }
}

/// `prefix`'s directory and file-name part.
fn split(prefix: &str) -> (PathBuf, String) {
    let p = Path::new(prefix);
    let base = p
        .file_name()
        .map_or_else(String::new, |n| n.to_string_lossy().into_owned());
    (parent_dir(p), base)
}

/// The names in `dir` that start `{base}.`.
fn named(dir: &Path, base: &str) -> Result<Vec<String>> {
    let lead = format!("{base}.");
    Ok(fs::read_dir(dir)
        .with_context(|| format!("reading {}", dir.display()))?
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.starts_with(&lead))
        .collect())
}

fn remove(path: &Path) -> std::io::Result<()> {
    if path.is_dir() {
        fs::remove_dir_all(path)
    } else {
        fs::remove_file(path)
    }
}

/// Remove what a pass left under `staging`: a stopped one's half.
pub fn discard(staging: &str) -> Result<()> {
    let (dir, base) = split(staging);
    if !dir.is_dir() {
        return Ok(());
    }
    for name in named(&dir, &base)? {
        let p = dir.join(&name);
        remove(&p).with_context(|| format!("removing {}", p.display()))?;
    }
    Ok(())
}

/// Every string in `v`, rewritten by `f`.
fn map_strings(v: &mut Value, f: &impl Fn(&str) -> String) {
    match v {
        Value::String(s) => *s = f(s),
        Value::Array(a) => a.iter_mut().for_each(|x| map_strings(x, f)),
        Value::Object(o) => o.values_mut().for_each(|x| map_strings(x, f)),
        _ => {}
    }
}

/// The files of the round `{base}` that manifest `v` names: its strings
/// that are `{base}.*` names beside it.
fn own_names(v: &Value, base: &str) -> BTreeSet<String> {
    fn walk(v: &Value, lead: &str, out: &mut BTreeSet<String>) {
        match v {
            Value::String(s) if s.starts_with(lead) && !s.contains('/') => {
                out.insert(s.clone());
            }
            Value::Array(a) => a.iter().for_each(|x| walk(x, lead, out)),
            Value::Object(o) => o.values().for_each(|x| walk(x, lead, out)),
            _ => {}
        }
    }
    let mut out = BTreeSet::new();
    walk(v, &format!("{base}."), &mut out);
    out
}

fn read_json(p: &Path) -> Result<Value> {
    let text = fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?;
    serde_json::from_str(&text).with_context(|| format!("parsing {}", p.display()))
}

/// Put the pass staged under `staging` (a run from `source`'s kind of
/// manifest) in the place of `out`'s round: its files renamed onto `out`,
/// the paths in its manifest with them, the manifest last; then the old
/// round's own files the new one no longer names are removed. Files the old
/// manifest never named, such as later rounds, are not the pass's to touch.
pub fn promote(source: &Path, staging: &str, out: &str) -> Result<()> {
    let (dir, staged_base) = split(staging);
    let (_, out_base) = split(out);
    let staged = annotated_path(source, staging);
    let target = annotated_path(source, out);

    let mut manifest = read_json(&staged)?;
    map_strings(&mut manifest, &|s| s.replace(&staged_base, &out_base));
    let new_names = own_names(&manifest, &out_base);
    let old_names = if target.is_file() {
        own_names(&read_json(&target)?, &out_base)
    } else {
        BTreeSet::new()
    };

    let staged_name = staged.file_name().map(|n| n.to_string_lossy().into_owned());
    for name in named(&dir, &staged_base)? {
        if Some(&name) == staged_name.as_ref() {
            continue;
        }
        let to = dir.join(format!("{out_base}{}", &name[staged_base.len()..]));
        if to.exists() {
            remove(&to).with_context(|| format!("replacing {}", to.display()))?;
        }
        fs::rename(dir.join(&name), &to).with_context(|| format!("moving {name} into place"))?;
    }
    fs::write(&staged, serde_json::to_string_pretty(&manifest)?)?;
    fs::rename(&staged, &target).with_context(|| format!("writing {}", target.display()))?;
    log::info!("the pass is in place: {}", target.display());

    for stale in old_names.difference(&new_names) {
        let p = dir.join(stale);
        if p.exists() {
            remove(&p).with_context(|| format!("removing {}", p.display()))?;
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "tests/staging.rs"]
mod tests;
