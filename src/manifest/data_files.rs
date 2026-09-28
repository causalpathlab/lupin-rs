//! lupin's data files, kept out of the binary: the Cell Ontology
//! (`cl-basic.obo`), the matching rules (`cl_matching.json`) and the curated
//! aliases (`cl_aliases.tsv`). Each is found on a search path and can be
//! amended without rebuilding:
//!
//! 1. **base**: the install's `share/lupin/` (`LUPIN_DATA_DIR` overrides),
//!    else the `data/` of the source this binary was built from (a checkout,
//!    or the crate `cargo install` unpacked) while it is still there, else
//!    the user cache, filled by download (the rules and aliases from this
//!    release's tag of the repository, the ontology from the rules'
//!    `ontology_url`) unless `LUPIN_OFFLINE` is set;
//! 2. **user**: `~/.config/lupin/` (`LUPIN_CONFIG_DIR` overrides);
//! 3. **project**: `lupin/` beside the run manifest;
//! 4. **run**: a file named on the command line (`--label-cl`, `--obo`).
//!
//! Rules layer key by key, aliases row by row, later layers winning; the
//! ontology is the most specific one found. [`ClData::sources`] says which
//! files were read, for the run's record.

use crate::annotate::celltype_tree::ClTerms;
use crate::annotate::cl_rules::{Aliases, MatchRules};
use anyhow::{Context, Result};
use log::{info, warn};
use serde_json::{json, Value};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

pub const RULES: &str = "cl_matching.json";
pub const ALIASES: &str = "cl_aliases.tsv";
pub const ONTOLOGY: &str = "cl-basic.obo";

/// Set to skip every download.
pub const OFFLINE_ENV: &str = "LUPIN_OFFLINE";
pub const DATA_DIR_ENV: &str = "LUPIN_DATA_DIR";
pub const CONFIG_DIR_ENV: &str = "LUPIN_CONFIG_DIR";

/// The directories searched, least specific first.
#[derive(Debug, Clone)]
pub struct SearchPath {
    /// The install's data (`share/lupin/`), when there is one.
    pub install: Option<PathBuf>,
    /// The `data/` of the source this binary was built from, if still there.
    pub source: Option<PathBuf>,
    /// lupin's download cache: the ontology at its top, the rules and
    /// aliases under `data/<release>/` (see [`Self::cached`]).
    pub cache: Option<PathBuf>,
    pub user: Option<PathBuf>,
    pub project: Option<PathBuf>,
}

impl SearchPath {
    /// The search path for a run whose manifest sits in `run_dir`.
    #[must_use]
    pub fn new(run_dir: Option<&Path>) -> Self {
        let env_dir = |k: &str| std::env::var_os(k).map(PathBuf::from);
        let install = env_dir(DATA_DIR_ENV).or_else(|| {
            let exe = std::env::current_exe().ok()?;
            let share = exe.parent()?.parent()?.join("share").join("lupin");
            share.is_dir().then_some(share)
        });
        let source =
            Some(Path::new(env!("CARGO_MANIFEST_DIR")).join("data")).filter(|d| d.is_dir());
        let cache = dirs::cache_dir().map(|d| d.join("lupin"));
        let user = env_dir(CONFIG_DIR_ENV)
            .or_else(|| dirs::home_dir().map(|h| h.join(".config").join("lupin")));
        let project = run_dir.map(|d| d.join("lupin"));
        Self {
            install,
            source,
            cache,
            user,
            project,
        }
    }

    /// Where the cache keeps `name`: the ontology whatever the release, the
    /// rules and aliases per release, as they are published per release.
    #[must_use]
    pub fn cached(&self, name: &str) -> Option<PathBuf> {
        let c = self.cache.as_ref()?;
        Some(if name == ONTOLOGY {
            c.join(name)
        } else {
            c.join("data").join(env!("CARGO_PKG_VERSION")).join(name)
        })
    }

    /// `name` in every layer that has it, least specific first. The base
    /// layer is the install's copy, else the cache's, downloaded when
    /// `download` says so.
    fn layers(&self, name: &str, download: bool) -> Vec<PathBuf> {
        let mut out = Vec::new();
        let shipped = [&self.install, &self.source]
            .into_iter()
            .flatten()
            .map(|d| d.join(name))
            .find(|p| p.is_file());
        match shipped {
            Some(p) => out.push(p),
            None => {
                if let Some(p) = self.cached(name) {
                    if p.is_file() || (download && fetch_into(name, &p, None).is_ok()) {
                        out.push(p);
                    }
                }
            }
        }
        out.extend(
            [&self.user, &self.project]
                .into_iter()
                .flatten()
                .map(|d| d.join(name))
                .filter(|p| p.is_file()),
        );
        out
    }

    /// Where the TUI remembers an alias: the project's file, else the user's.
    #[must_use]
    pub fn amend_aliases(&self) -> Option<PathBuf> {
        self.project
            .as_ref()
            .or(self.user.as_ref())
            .map(|d| d.join(ALIASES))
    }
}

/// The Cell Ontology data a run uses.
pub struct ClData {
    pub rules: MatchRules,
    pub aliases: Aliases,
    /// The ontology file, when one is at hand.
    pub ontology: Option<PathBuf>,
    /// Every file read, for the run's record.
    pub sources: Vec<String>,
    pub search: SearchPath,
}

impl ClData {
    /// Load the rules and aliases along `search`, then `run_aliases` (a
    /// `--label-cl` file), and find the ontology (`explicit_obo` first).
    pub fn load(
        search: SearchPath,
        explicit_obo: Option<&str>,
        run_aliases: Option<&str>,
    ) -> Result<Self> {
        let online = std::env::var_os(OFFLINE_ENV).is_none();
        let mut sources = Vec::new();

        let mut layers = Vec::new();
        for p in search.layers(RULES, online) {
            let text =
                fs::read_to_string(&p).with_context(|| format!("reading {}", p.display()))?;
            layers.push(
                serde_json::from_str::<Value>(&text)
                    .with_context(|| format!("{} is not JSON", p.display()))?,
            );
            sources.push(p.display().to_string());
        }
        if layers.is_empty() {
            warn!(
                "no {RULES} found (see `lupin data where`): matching cell-type labels to the \
                 Cell Ontology literally, by name and exact synonym only"
            );
        }
        let rules = MatchRules::from_layers(&layers)?;

        let mut aliases = Aliases::default();
        let mut files = search.layers(ALIASES, online);
        files.extend(run_aliases.map(PathBuf::from));
        for p in files {
            let text =
                fs::read_to_string(&p).with_context(|| format!("reading {}", p.display()))?;
            aliases.add_tsv(&text, &p.display().to_string())?;
            sources.push(p.display().to_string());
        }
        if !aliases.is_empty() {
            info!(
                "{} cell-type alias(es) to Cell Ontology terms",
                aliases.len()
            );
        }

        let ontology = match explicit_obo {
            Some(p) => Some(PathBuf::from(p)),
            None => search
                .layers(ONTOLOGY, false)
                .pop()
                .or_else(|| fetch_ontology(&search, &rules, online)),
        };
        if let Some(p) = &ontology {
            sources.push(p.display().to_string());
        }
        Ok(Self {
            rules,
            aliases,
            ontology,
            sources,
            search,
        })
    }

    /// The ontology parsed under the rules, with the aliases; `None` without
    /// an ontology.
    pub fn terms(&self) -> Result<Option<ClTerms>> {
        self.ontology
            .as_ref()
            .map(|p| {
                let text =
                    fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?;
                Ok(ClTerms::parse(&text, &self.rules).with_aliases(self.aliases.clone()))
            })
            .transpose()
    }

    /// What the run used, for its record.
    #[must_use]
    pub fn record(&self, release: Option<&str>) -> Value {
        json!({ "files": self.sources, "ontology_release": release })
    }
}

/// The ontology downloaded into the cache from the rules' `ontology_url`.
fn fetch_ontology(search: &SearchPath, rules: &MatchRules, online: bool) -> Option<PathBuf> {
    let to = search.cached(ONTOLOGY)?;
    let url = rules.ontology_url.as_deref()?;
    if !online {
        info!("{OFFLINE_ENV} is set: not downloading the Cell Ontology");
        return None;
    }
    match fetch_into(ONTOLOGY, &to, Some(url)) {
        Ok(()) => Some(to),
        Err(e) => {
            warn!("could not download the Cell Ontology ({e:#}); grouping cell types by shared markers instead");
            None
        }
    }
}

/// Where this release's copy of data file `name` is published: the
/// repository (from the package metadata) at this version's tag.
#[must_use]
pub fn release_url(name: &str) -> String {
    let repo = env!("CARGO_PKG_REPOSITORY").replace("https://github.com/", "");
    format!(
        "https://raw.githubusercontent.com/{repo}/v{}/data/{name}",
        env!("CARGO_PKG_VERSION")
    )
}

/// Download `name` (from `url`, else this release's copy) to `to`, streamed
/// to a side file and renamed once whole, so a half-finished download is
/// never read.
pub fn fetch_into(name: &str, to: &Path, url: Option<&str>) -> Result<()> {
    let url = url.map_or_else(|| release_url(name), String::from);
    info!("downloading {name} from {url}");
    let response = ureq::AgentBuilder::new()
        .timeout(Duration::from_secs(60))
        .build()
        .get(&url)
        .call()
        .context("request failed")?;
    fs::create_dir_all(to.parent().context("no cache directory")?)?;
    let tmp = to.with_extension("part");
    let mut file = fs::File::create(&tmp)?;
    std::io::copy(&mut response.into_reader(), &mut file).context("reading the response")?;
    drop(file);
    if name == ONTOLOGY {
        let head = fs::read_to_string(&tmp).unwrap_or_default();
        anyhow::ensure!(head.contains("[Term]"), "the response is not an OBO file");
    }
    fs::rename(&tmp, to)?;
    info!("cached {}", to.display());
    Ok(())
}

/// `lupin data`.
#[derive(clap::Args, Debug)]
pub struct DataArgs {
    #[command(subcommand)]
    pub cmd: DataCmd,
}

#[derive(clap::Subcommand, Debug)]
pub enum DataCmd {
    /// Show the search path and the files each layer contributes
    Where {
        /// A run manifest, for its project layer
        #[arg(long, short = 'f')]
        from: Option<Box<str>>,
    },
    /// Download this release's rules and aliases, and the Cell Ontology, into
    /// the cache (for machines that will run offline)
    Fetch {
        /// Replace what the cache already holds
        #[arg(long)]
        force: bool,
    },
}

pub fn run_data(args: &DataArgs) -> Result<()> {
    match &args.cmd {
        DataCmd::Where { from } => {
            let run_dir = from
                .as_deref()
                .map(|f| super::run::load(f).map(|l| l.dir))
                .transpose()?;
            let search = SearchPath::new(run_dir.as_deref());
            let show = |label: &str, dir: Option<PathBuf>| {
                let dir = dir.map_or_else(|| "-".to_string(), |d| d.display().to_string());
                println!("{label:<9} {dir}");
            };
            show("install", search.install.clone());
            show("source", search.source.clone());
            show("cache", search.cache.clone());
            show("user", search.user.clone());
            show("project", search.project.clone());
            println!();
            for name in [RULES, ALIASES, ONTOLOGY] {
                let found = search.layers(name, false);
                println!("{name}:");
                if found.is_empty() {
                    println!(
                        "  (none{})",
                        if name == ONTOLOGY {
                            "; `lupin data fetch` or --obo"
                        } else {
                            "; `lupin data fetch`"
                        }
                    );
                }
                for p in found {
                    println!("  {}", p.display());
                }
            }
            Ok(())
        }
        DataCmd::Fetch { force } => {
            let search = SearchPath::new(None);
            for name in [RULES, ALIASES] {
                let to = search.cached(name).context("no cache directory")?;
                if *force || !to.is_file() {
                    fetch_into(name, &to, None)?;
                }
                println!("{}", to.display());
            }
            let rules = ClData::load(search.clone(), None, None)?.rules;
            let to = search.cached(ONTOLOGY).context("no cache directory")?;
            if *force || !to.is_file() {
                let url = rules
                    .ontology_url
                    .as_deref()
                    .context("the matching rules name no ontology_url")?;
                fetch_into(ONTOLOGY, &to, Some(url))?;
            }
            println!("{}", to.display());
            Ok(())
        }
    }
}

#[cfg(test)]
#[path = "tests/data_files.rs"]
mod tests;
