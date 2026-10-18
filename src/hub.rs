//! The Hugging Face Hub, as every model lupin loads reaches it: each file is
//! read from the cache `HF_HOME` names (else `~/.cache/huggingface`) when it
//! is there, and downloaded only when neither `LUPIN_OFFLINE` nor
//! `HF_HUB_OFFLINE` (read as the Hub's own clients read it) is set.

use crate::manifest::data_files::OFFLINE_ENV;
use anyhow::{Context, Result};
use hf_hub::api::sync::{ApiBuilder, ApiRepo};
use hf_hub::{Cache, Repo, RepoType};
use log::info;
use std::ffi::OsString;
use std::path::PathBuf;

/// The Hub's own switch against downloads.
pub const HF_OFFLINE_ENV: &str = "HF_HUB_OFFLINE";

/// The switch that rules out a download, if one does: `LUPIN_OFFLINE` set
/// to anything, else `HF_HUB_OFFLINE` set to a true value (`1`, `true`,
/// `yes`, `on`, any case).
fn offline_switch(lupin: Option<OsString>, hf: Option<String>) -> Option<&'static str> {
    let on = |v: &str| ["1", "true", "yes", "on"].contains(&v.trim().to_lowercase().as_str());
    if lupin.is_some() {
        Some(OFFLINE_ENV)
    } else {
        hf.filter(|v| on(v)).map(|_| HF_OFFLINE_ENV)
    }
}

/// A model repository at a pinned revision, in the environment's cache.
pub struct Model {
    id: String,
    revision: String,
    cache: Cache,
    /// The switch set against downloads, if any.
    offline: Option<&'static str>,
}

impl Model {
    #[must_use]
    pub fn new(id: &str, revision: &str) -> Self {
        let offline = offline_switch(
            std::env::var_os(OFFLINE_ENV),
            std::env::var(HF_OFFLINE_ENV).ok(),
        );
        Self::in_cache(id, revision, Cache::from_env(), offline)
    }

    fn in_cache(id: &str, revision: &str, cache: Cache, offline: Option<&'static str>) -> Self {
        Self {
            id: id.to_string(),
            revision: revision.to_string(),
            cache,
            offline,
        }
    }

    fn repo(&self) -> Repo {
        Repo::with_revision(self.id.clone(), RepoType::Model, self.revision.clone())
    }

    /// File `name`, when the cache has it: under the snapshot the revision's
    /// ref names, or, for a commit, under its own snapshot, which the Hub's
    /// Python tools write with no ref to it.
    #[must_use]
    fn cached(&self, name: &str) -> Option<PathBuf> {
        let repo = self.repo();
        self.cache.repo(repo.clone()).get(name).or_else(|| {
            let commit =
                self.revision.len() == 40 && self.revision.bytes().all(|b| b.is_ascii_hexdigit());
            let path = (self.cache.path().join(repo.folder_name()))
                .join("snapshots")
                .join(&self.revision)
                .join(name);
            (commit && path.exists()).then_some(path)
        })
    }

    /// The Hub, to fetch what the cache lacks, drawing its bars unless the
    /// TUI reads this command's progress; offline, an error naming the
    /// switch.
    fn online(&self) -> Result<ApiRepo> {
        let at = self.cache.path().display();
        if let Some(switch) = self.offline {
            anyhow::bail!(
                "{}@{} is not cached under {at} and {switch} is set",
                self.id,
                self.revision
            );
        }
        info!(
            "fetching {}@{} from the Hub (cached under {at})",
            self.id, self.revision
        );
        // Not `ApiBuilder::from_env`, which would replace the cache (a test's).
        let mut api =
            ApiBuilder::from_cache(self.cache.clone()).with_progress(!crate::progress::enabled());
        if let Ok(endpoint) = std::env::var("HF_ENDPOINT") {
            api = api.with_endpoint(endpoint);
        }
        Ok(api
            .build()
            .context("Hugging Face Hub client")?
            .repo(self.repo()))
    }

    /// Files `names`: from the cache when it has every one, else from the
    /// Hub ([`Self::online`]).
    pub fn files<const N: usize>(&self, names: [&str; N]) -> Result<[PathBuf; N]> {
        let all = match names.iter().map(|n| self.cached(n)).collect() {
            Some(all) => all,
            None => {
                let hub = self.online()?;
                (names.iter())
                    .map(|n| {
                        hub.get(n)
                            .with_context(|| format!("fetching {n} of {}", self.id))
                    })
                    .collect::<Result<Vec<_>>>()?
            }
        };
        // One path per name, by construction.
        Ok(<[PathBuf; N]>::try_from(all).unwrap_or_else(|_| unreachable!()))
    }
}

#[cfg(test)]
#[path = "hub_tests.rs"]
mod tests;
