//! `LUPIN_VERSION`: the crate's version, and the commit it was built from
//! when built in a git checkout (a crates.io build has none), so two builds
//! of one version can be told apart.

use std::process::Command;

fn main() {
    // HEAD moves on a checkout, its branch's ref on a commit (or
    // packed-refs), and the sources on an edit (`-dirty`). Only paths that
    // exist: cargo reruns for a missing one every build. Optional locks off,
    // so `describe` leaves the index alone.
    let git = |args: &[&str]| {
        Command::new("git")
            .env("GIT_OPTIONAL_LOCKS", "0")
            .args(args)
            .output()
            .ok()
            .filter(|o| o.status.success())
            .and_then(|o| String::from_utf8(o.stdout).ok())
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
    };
    let mut watched = vec!["src".to_string()];
    for path in ["HEAD", "packed-refs"] {
        watched.extend(git(&["rev-parse", "--git-path", path]));
    }
    if let Some(r) = git(&["symbolic-ref", "-q", "HEAD"]) {
        watched.extend(git(&["rev-parse", "--git-path", &r]));
    }
    for w in watched.iter().filter(|w| std::path::Path::new(w).exists()) {
        println!("cargo:rerun-if-changed={w}");
    }
    let version = env!("CARGO_PKG_VERSION");
    let commit = git(&["describe", "--always", "--dirty", "--exclude", "*"]);
    let full = commit.map_or_else(|| version.to_string(), |c| format!("{version} ({c})"));
    println!("cargo:rustc-env=LUPIN_VERSION={full}");
}
