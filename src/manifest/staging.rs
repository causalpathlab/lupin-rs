//! Passes write under a staging prefix beside their own and are promoted
//! onto it when they finish, so a stopped or failed pass leaves the round
//! it would have replaced as it was.

use super::run::annotated_path;
use anyhow::{Context, Result};
use serde_json::Value;
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

/// The prefix a pass to `out` writes under until it is promoted.
#[must_use]
pub fn staging_prefix(out: &str) -> String {
    format!("{out}.staging")
}

/// `prefix`'s directory and file-name part.
fn split(prefix: &str) -> (PathBuf, String) {
    let p = Path::new(prefix);
    let dir = match p.parent() {
        Some(d) if !d.as_os_str().is_empty() => d.to_path_buf(),
        _ => PathBuf::from("."),
    };
    let base = p
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    (dir, base)
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

fn strings(v: &Value, out: &mut BTreeSet<String>) {
    match v {
        Value::String(s) => {
            out.insert(s.clone());
        }
        Value::Array(a) => a.iter().for_each(|x| strings(x, out)),
        Value::Object(o) => o.values().for_each(|x| strings(x, out)),
        _ => {}
    }
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
    let own = |names: BTreeSet<String>| -> BTreeSet<String> {
        let lead = format!("{out_base}.");
        names
            .into_iter()
            .filter(|n| n.starts_with(&lead) && !n.contains('/'))
            .collect()
    };
    let mut new_names = BTreeSet::new();
    strings(&manifest, &mut new_names);
    let new_names = own(new_names);
    let old_names = match target.is_file() {
        true => {
            let mut s = BTreeSet::new();
            strings(&read_json(&target)?, &mut s);
            own(s)
        }
        false => BTreeSet::new(),
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
