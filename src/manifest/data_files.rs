//! lupin's data files, kept out of the binary: the Cell Ontology
//! (`cl-basic.obo`) and the matching rules (`cl_matching.json`). Each is found
//! on a search path and can be amended without rebuilding:
//!
//! 1. **base**: the install's `share/lupin/` (`LUPIN_DATA_DIR` overrides),
//!    else the `data/` of the source this binary was built from (a checkout,
//!    or the crate `cargo install` unpacked) while it is still there, else
//!    the user cache, filled by download (the rules from this release's tag
//!    of the repository, the ontology from the rules' `ontology_url`) unless
//!    `LUPIN_OFFLINE` is set;
//! 2. **user**: `~/.config/lupin/` (`LUPIN_CONFIG_DIR` overrides);
//! 3. **project**: `lupin/` beside the run manifest;
//! 4. **run**: a file named on the command line (`--label-cl`, `--obo`).
//!
//! Label aliases (`label<TAB>CL:id<TAB>note`, for labels name matching
//! cannot settle) are data about a marker panel, not about lupin: none ship
//! with it. They are read from the panel's sidecar (`x.tsv.gz` →
//! `x.cl.tsv`, see [`sidecar_of`]), then the user's and the project's
//! `cl_aliases.tsv`, then `--label-cl`, row by row, later layers winning.
//! Rules layer key by key; the ontology is the most specific one found.
//! [`ClData::record`] says which files were read, for the run's record.
//!
//! The Gene Ontology (`go-basic.obo`), which names the terms of a GO pass, and
//! a species' GO annotations (`goa_human.gaf.gz`, ...) are found the same way
//! ([`go_ontology`], [`go_annotations`]): the most specific layer that has
//! them, else downloaded into the cache.

use crate::annotate::celltype_tree::ClTerms;
use crate::annotate::cl_rules::{Aliases, MatchRules};
use crate::annotate::go_signature::Species;
use anyhow::{Context, Result};
use log::{info, warn};
use serde_json::{json, Value};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::time::Duration;

pub const RULES: &str = "cl_matching.json";
pub const ALIASES: &str = "cl_aliases.tsv";
/// Which cell types precede which, for `lupin trajectory` and the TUI's order view.
pub const PRECEDENCE: &str = "precedence.tsv";
pub const ONTOLOGY: &str = "cl-basic.obo";
pub const GO_ONTOLOGY: &str = "go-basic.obo";
pub const GO_ONTOLOGY_URL: &str = "https://purl.obolibrary.org/obo/go/go-basic.obo";

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
    /// The marker panel's alias sidecar ([`sidecar_of`]), whether or not it
    /// exists yet.
    pub panel: Option<PathBuf>,
}

/// The alias sidecar of marker panel `panel`: its table extensions
/// (`.gz`, then `.tsv`, `.txt` or `.csv`) replaced by `.cl.tsv`, beside it.
#[must_use]
pub fn sidecar_of(panel: &Path) -> Option<PathBuf> {
    let name = panel.file_name()?.to_str()?;
    let base = name.strip_suffix(".gz").unwrap_or(name);
    let base = [".tsv", ".txt", ".csv"]
        .iter()
        .find_map(|e| base.strip_suffix(e))
        .unwrap_or(base);
    (!base.is_empty()).then(|| panel.with_file_name(format!("{base}.cl.tsv")))
}

impl SearchPath {
    /// The search path for a run whose manifest sits in `run_dir`. Under
    /// test, only the project layer: nothing from the machine the tests run
    /// on, and no downloads.
    #[must_use]
    pub fn new(run_dir: Option<&Path>) -> Self {
        let project = run_dir.map(|d| d.join("lupin"));
        if cfg!(test) {
            return Self {
                install: None,
                source: None,
                cache: None,
                user: None,
                project,
                panel: None,
            };
        }
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
        Self {
            install,
            source,
            cache,
            user,
            project,
            panel: None,
        }
    }

    /// The search path with marker panel `panel`'s sidecar as the first
    /// alias layer; an empty name leaves it without one.
    #[must_use]
    pub fn with_panel(mut self, panel: &str) -> Self {
        self.panel = (!panel.is_empty())
            .then(|| sidecar_of(Path::new(panel)))
            .flatten();
        self
    }

    /// The alias files, least specific first: the panel's sidecar, then the
    /// user's and the project's `cl_aliases.tsv`, those that exist.
    #[must_use]
    pub fn alias_layers(&self) -> Vec<PathBuf> {
        [&self.panel]
            .into_iter()
            .flatten()
            .cloned()
            .chain(
                [&self.user, &self.project]
                    .into_iter()
                    .flatten()
                    .map(|d| d.join(ALIASES)),
            )
            .filter(|p| p.is_file())
            .collect()
    }

    /// Where a label's alias chosen in the TUI is written: the panel's
    /// sidecar when its directory takes a new file, else the project's
    /// (else the user's) `cl_aliases.tsv`.
    #[must_use]
    pub fn alias_target(&self) -> Option<PathBuf> {
        let writable = |p: &Path| {
            p.parent()
                .is_some_and(|d| d.metadata().is_ok_and(|m| !m.permissions().readonly()))
        };
        self.panel
            .clone()
            .filter(|p| writable(p))
            .or_else(|| self.amend(ALIASES))
    }

    /// Where the cache keeps `name`: an ontology whatever the release, the
    /// rules and aliases per release, as they are published per release.
    #[must_use]
    pub fn cached(&self, name: &str) -> Option<PathBuf> {
        let c = self.cache.as_ref()?;
        Some(if is_external(name) {
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

    /// A file named on the command line: as given, else beside the run
    /// (the project layer's parent); `None`, with a warning, when neither
    /// exists.
    #[must_use]
    pub fn run_file(&self, f: &str) -> Option<PathBuf> {
        let given = PathBuf::from(f);
        let beside = self
            .project
            .as_ref()
            .and_then(|p| p.parent())
            .map(|d| d.join(f));
        let found = [Some(given), beside]
            .into_iter()
            .flatten()
            .find(|p| p.is_file());
        if found.is_none() {
            warn!("{f} not found; going on without it");
        }
        found
    }

    /// Where the TUI writes data file `name`: the project's, else the user's.
    #[must_use]
    pub fn amend(&self, name: &str) -> Option<PathBuf> {
        self.project
            .as_ref()
            .or(self.user.as_ref())
            .map(|d| d.join(name))
    }

    /// `name` in the user's and then the project's layer, the files that
    /// exist, each with its layer's name.
    #[must_use]
    pub fn user_and_project_files(&self, name: &str) -> Vec<(&'static str, PathBuf)> {
        [("user", &self.user), ("project", &self.project)]
            .into_iter()
            .filter_map(|(layer, d)| d.as_ref().map(|d| (layer, d.join(name))))
            .filter(|(_, p)| p.is_file())
            .collect()
    }

    /// `name`'s text from the user's and then the project's layer, the files
    /// that exist.
    pub fn user_and_project(&self, name: &str) -> Result<Vec<String>> {
        self.user_and_project_files(name)
            .into_iter()
            .map(|(_, p)| {
                fs::read_to_string(&p).with_context(|| format!("reading {}", p.display()))
            })
            .collect()
    }
}

/// The Cell Ontology data a run uses.
pub struct ClData {
    pub rules: MatchRules,
    pub aliases: Aliases,
    /// The ontology file, when one is at hand.
    pub ontology: Option<PathBuf>,
    /// The rules and alias files read, in layer order.
    pub rule_files: Vec<PathBuf>,
    pub alias_files: Vec<PathBuf>,
    pub search: SearchPath,
    /// The ontology, parsed on first use.
    pub parsed: OnceLock<Option<ClTerms>>,
}

/// Whether [`ClData::load`] may download what it cannot find. Only an
/// annotate pass and the TUI's start do: a relabel, a preview or a rescore
/// reuses what its pass recorded ([`ClData::from_record`]).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Fetch {
    Allowed,
    Never,
}

impl ClData {
    /// Load the rules and aliases along `search`, then `run_aliases` (a
    /// `--label-cl` file), and find the ontology (`explicit_obo` first).
    /// Run files are taken as given, else beside the run.
    pub fn load(
        search: SearchPath,
        explicit_obo: Option<&str>,
        run_aliases: Option<&str>,
        fetch: Fetch,
    ) -> Result<Self> {
        let online = fetch == Fetch::Allowed && std::env::var_os(OFFLINE_ENV).is_none();
        let rule_files = search.layers(RULES, online);
        if rule_files.is_empty() {
            warn!(
                "no {RULES} found (see `lupin data where`): matching cell-type labels to the \
                 Cell Ontology literally, by name and exact synonym only"
            );
        }
        let mut alias_files = search.alias_layers();
        alias_files.extend(run_aliases.and_then(|f| search.run_file(f)));
        let ontology = match explicit_obo {
            Some(f) => search.run_file(f),
            None => search.layers(ONTOLOGY, false).pop(),
        };
        let mut data = Self::from_files(rule_files, alias_files, ontology, search)?;
        if data.ontology.is_none() && explicit_obo.is_none() {
            data.ontology = fetch_ontology(&data.search, &data.rules, online);
        }
        Ok(data)
    }

    /// The data a pass recorded ([`Self::record`]), read from exactly those
    /// files; `None` when there is no record or a file is gone.
    pub fn from_record(record: &Value, search: SearchPath) -> Result<Option<Self>> {
        let paths = |key: &str| -> Option<Vec<PathBuf>> {
            record
                .get(key)?
                .as_array()?
                .iter()
                .map(|v| v.as_str().map(PathBuf::from))
                .collect()
        };
        let (Some(rules), Some(aliases)) = (paths("rules"), paths("aliases")) else {
            return Ok(None);
        };
        let ontology = record
            .get("ontology")
            .and_then(Value::as_str)
            .map(PathBuf::from);
        let all_there = rules
            .iter()
            .chain(&aliases)
            .chain(&ontology)
            .all(|p| p.is_file());
        if !all_there {
            warn!("files the pass recorded are gone; finding the Cell Ontology data afresh");
            return Ok(None);
        }
        Self::from_files(rules, aliases, ontology, search).map(Some)
    }

    fn from_files(
        rule_files: Vec<PathBuf>,
        alias_files: Vec<PathBuf>,
        ontology: Option<PathBuf>,
        search: SearchPath,
    ) -> Result<Self> {
        let mut layers = Vec::new();
        for p in &rule_files {
            let text = fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?;
            layers.push(
                serde_json::from_str::<Value>(&text)
                    .with_context(|| format!("{} is not JSON", p.display()))?,
            );
        }
        let rules = MatchRules::from_layers(&layers)?;
        let mut aliases = Aliases::default();
        for p in &alias_files {
            let text = fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?;
            aliases.add_tsv(&text, &p.display().to_string());
        }
        if !aliases.is_empty() {
            info!(
                "{} cell-type alias(es) to Cell Ontology terms",
                aliases.len()
            );
        }
        Ok(Self {
            rules,
            aliases,
            ontology,
            rule_files,
            alias_files,
            search,
            parsed: OnceLock::new(),
        })
    }

    /// The ontology parsed under the rules, with the aliases, parsed once;
    /// `None` without an ontology.
    pub fn terms(&self) -> Result<Option<&ClTerms>> {
        if self.parsed.get().is_none() {
            let parsed = self.parse()?;
            let _ = self.parsed.set(parsed);
        }
        Ok(self.parsed.get().and_then(Option::as_ref))
    }

    /// [`Self::terms`], owned.
    pub fn into_terms(mut self) -> Result<Option<ClTerms>> {
        match self.parsed.take() {
            Some(t) => Ok(t),
            None => self.parse(),
        }
    }

    fn parse(&self) -> Result<Option<ClTerms>> {
        self.ontology
            .as_ref()
            .map(|p| {
                let text =
                    fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?;
                Ok(ClTerms::parse(&text, &self.rules).with_aliases(self.aliases.clone()))
            })
            .transpose()
    }

    /// What the run used, for its record: each file by its role, as an
    /// absolute path, so a later rescore reads the same ones
    /// ([`Self::from_record`]) wherever it runs from.
    #[must_use]
    pub fn record(&self, release: Option<&str>) -> Value {
        let abs = |p: &PathBuf| absolute(p).display().to_string();
        json!({
            "rules": self.rule_files.iter().map(abs).collect::<Vec<_>>(),
            "aliases": self.alias_files.iter().map(abs).collect::<Vec<_>>(),
            "ontology": self.ontology.as_ref().map(abs),
            "ontology_release": release,
        })
    }
}

/// Append `line` to `file`, creating it (and its directory) with the
/// comment `header` first.
pub fn append_line(file: &Path, header: &str, line: &str) -> Result<()> {
    use std::io::Write;
    if let Some(dir) = file.parent() {
        fs::create_dir_all(dir)?;
    }
    let new = !file.exists();
    let mut f = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(file)?;
    if new {
        for h in header.lines() {
            writeln!(f, "# {h}")?;
        }
    }
    writeln!(f, "{line}")?;
    Ok(())
}

/// Append a TSV row of `fields` to `file` as [`append_line`] does; a tab or a
/// newline inside a field becomes a space.
pub fn append_row(file: &Path, header: &str, fields: &[&str]) -> Result<()> {
    let row: Vec<String> = fields
        .iter()
        .map(|f| f.replace(['\t', '\n'], " "))
        .collect();
    append_line(file, header, &row.join("\t"))
}

/// `p` made absolute (and canonical when it exists).
pub(crate) fn absolute(p: &Path) -> PathBuf {
    p.canonicalize().unwrap_or_else(|_| {
        std::env::current_dir()
            .map(|d| d.join(p))
            .unwrap_or_else(|_| p.to_path_buf())
    })
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

/// Whether data file `name` is an ontology, checked to be OBO when
/// downloaded.
fn is_obo(name: &str) -> bool {
    name.ends_with(".obo")
}

/// Whether data file `name` is published by its source rather than with a
/// lupin release, so one cached copy serves every release.
fn is_external(name: &str) -> bool {
    is_obo(name) || name.contains(".gaf")
}

/// `name`: the most specific layer's copy, else the cache's, downloaded from
/// `url` unless [`OFFLINE_ENV`] is set. `flag` is how to name a file instead.
fn external_file(run_dir: Option<&Path>, name: &str, url: &str, flag: &str) -> Result<PathBuf> {
    let search = SearchPath::new(run_dir);
    if let Some(p) = search.layers(name, false).pop() {
        return Ok(p);
    }
    anyhow::ensure!(
        std::env::var_os(OFFLINE_ENV).is_none(),
        "{OFFLINE_ENV} is set and no {name} is at hand; pass {flag} or run `lupin data fetch`"
    );
    let to = search.cached(name).context("no cache directory")?;
    fetch_into(name, &to, Some(url))
        .with_context(|| format!("downloading {name}; pass {flag} instead"))?;
    Ok(to)
}

/// The Gene Ontology ([`external_file`] from [`GO_ONTOLOGY_URL`]).
pub fn go_ontology(run_dir: Option<&Path>) -> Result<PathBuf> {
    external_file(run_dir, GO_ONTOLOGY, GO_ONTOLOGY_URL, "--go-obo")
}

/// `species`' GO annotations ([`external_file`] from the GO Consortium).
pub fn go_annotations(species: Species, run_dir: Option<&Path>) -> Result<PathBuf> {
    external_file(run_dir, species.gaf_file(), &species.gaf_url(), "--gaf")
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
    if is_obo(name) {
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
    /// Download this release's matching rules, the Cell Ontology, the Gene
    /// Ontology and the human and mouse GO annotations into the cache (for
    /// machines that will run offline)
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
            let gafs = [Species::Human, Species::Mouse].map(Species::gaf_file);
            println!(
                "aliases   the panel's sidecar (x.tsv.gz → x.cl.tsv), then user and project {ALIASES}, then --label-cl"
            );
            for p in search.alias_layers() {
                println!("  {}", p.display());
            }
            for name in [RULES, ONTOLOGY, GO_ONTOLOGY].into_iter().chain(gafs) {
                let found = search.layers(name, false);
                println!("{name}:");
                if found.is_empty() {
                    println!(
                        "  (none; `lupin data fetch`{})",
                        match name {
                            ONTOLOGY => " or --obo",
                            GO_ONTOLOGY => " or --go-obo",
                            _ if is_external(name) => " or --gaf",
                            _ => "",
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
            {
                let name = RULES;
                let to = search.cached(name).context("no cache directory")?;
                if *force || !to.is_file() {
                    fetch_into(name, &to, None)?;
                }
                println!("{}", to.display());
            }
            let rules = ClData::load(search.clone(), None, None, Fetch::Never)?.rules;
            let to = search.cached(ONTOLOGY).context("no cache directory")?;
            if *force || !to.is_file() {
                let url = rules
                    .ontology_url
                    .as_deref()
                    .context("the matching rules name no ontology_url")?;
                fetch_into(ONTOLOGY, &to, Some(url))?;
            }
            println!("{}", to.display());
            let go = [(GO_ONTOLOGY, GO_ONTOLOGY_URL.to_string())]
                .into_iter()
                .chain([Species::Human, Species::Mouse].map(|s| (s.gaf_file(), s.gaf_url())));
            for (name, url) in go {
                let to = search.cached(name).context("no cache directory")?;
                if *force || !to.is_file() {
                    fetch_into(name, &to, Some(&url))?;
                }
                println!("{}", to.display());
            }
            Ok(())
        }
    }
}

#[cfg(test)]
#[path = "tests/data_files.rs"]
mod tests;
