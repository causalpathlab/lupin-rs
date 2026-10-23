//! `LUPIN_VERSION`: the crate's version, and the commit it was built from
//! when built in a git checkout (a crates.io build has none), so two builds
//! of one version can be told apart.

use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/index");
    let version = env!("CARGO_PKG_VERSION");
    let commit = Command::new("git")
        .args(["describe", "--always", "--dirty", "--exclude", "*"])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty());
    let full = match commit {
        Some(c) => format!("{version} ({c})"),
        None => version.to_string(),
    };
    println!("cargo:rustc-env=LUPIN_VERSION={full}");
}
