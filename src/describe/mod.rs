//! `lupin describe` — compose a short sentence from annotate evidence.
//!
//! The template composer **never decides**: evidence comes from annotate artifacts (and
//! optional keyword incidence from a `word-graph` vocab). With `--lm` a small local model
//! ([`lm`]) proposes a label among each cluster's choices with one sentence on why
//! ([`consensus`]), and describes a trajectory's order ([`trajectory`]); nothing is
//! applied, the user confirms it in the TUI. Every sentence, template or model's, is
//! citation-checked against the entities its evidence allows ([`citation_check`]).

pub(crate) mod claims;
pub mod consensus;
pub mod lm;
pub mod trajectory;

use enrichment::UNASSIGNED_LABEL;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use crate::annotate::outputs::{
    ANNOT_PARQUET, ARGMAX_TSV, CLUSTER_CELLTYPE_Q_VALUES, CLUSTER_TERM_Q,
};
use crate::gene_text::vocab::Vocabulary;
use anyhow::{bail, Context, Result};
use clap::Args;
use graph_embedding_util::type_annotation::CELL_LABELS_SUFFIX;
use legume_numeric::matrix::common_io::open_buf_reader;
use legume_numeric::matrix::parquet::peek_parquet_field_names;
use log::info;
use serde_json::json;

#[derive(Args, Debug)]
pub struct DescribeArgs {
    #[arg(
        long,
        short = 'f',
        help = "Run manifest or annotate output prefix (reads projection, enrichment, or argmax artifacts)"
    )]
    pub from: Box<str>,

    #[arg(
        long,
        help = "Optional `word-graph` / `text-qc` prefix; fishes keywords from `{prefix}.feature_word.edges.tsv` via each cluster's markers"
    )]
    pub text_prefix: Option<Box<str>>,

    #[arg(
        long,
        help = "Optional marker TSV (`gene<TAB>celltype`) when annotate did not write marker tables"
    )]
    pub markers: Option<Box<str>>,

    #[arg(
        long,
        help = "Optional feature×H embedding parquet; expands each cluster's marker set by nearest neighbours before fishing keywords"
    )]
    pub feature_embedding: Option<Box<str>>,

    #[arg(
        long = "neighbour-k",
        default_value_t = 5,
        help = "How many embedding neighbours to add per cluster when `--feature-embedding` is set"
    )]
    pub neighbour_k: usize,

    #[arg(
        long = "fdr-alpha",
        default_value_t = 0.1,
        help = "FDR gate for calling a runner-up significant (must match annotate --fdr-alpha)"
    )]
    pub fdr_alpha: f32,

    #[arg(
        long,
        short = 'o',
        help = "Output prefix; without it, the TUI opens on the round (or the trajectory) and runs the consensus"
    )]
    pub out: Option<Box<str>>,

    #[arg(
        long,
        requires = "lm",
        help = "With --lm: write only the model's cache (`{out}.consensus.jsonl`), not the description files; what the TUI runs"
    )]
    pub cache_only: bool,

    #[command(flatten)]
    pub lm: LmArgs,
}

/// The small local model `describe` writes with ([`lm`]).
#[derive(Args, Debug, Clone, PartialEq, Eq)]
pub struct LmArgs {
    #[arg(
        long,
        help = "Propose a label per cluster among its choices, with a citation-checked sentence on why (or describe a trajectory's order), by a small local model; nothing is applied"
    )]
    pub lm: bool,

    #[arg(long, value_enum, default_value_t = lm::LmModel::default(), help = "The model (pinned revision, fetched once from the Hub)")]
    pub lm_model: lm::LmModel,

    #[arg(
        long,
        help = "Local model directory (config.json, tokenizer.json, *.safetensors); skips the Hub"
    )]
    pub lm_dir: Option<String>,

    #[arg(long, value_enum, default_value_t = lm::default_device(), help = "Compute device for the model")]
    pub device: crate::gene_text::cli::ComputeDevice,
}

impl Default for LmArgs {
    fn default() -> Self {
        Self {
            lm: false,
            lm_model: lm::LmModel::default(),
            lm_dir: None,
            device: lm::default_device(),
        }
    }
}

impl LmArgs {
    #[must_use]
    pub fn spec(&self) -> lm::LmSpec {
        lm::LmSpec {
            model: self.lm_model,
            local_dir: self.lm_dir.as_deref().map(std::path::PathBuf::from),
            device: self.device,
        }
    }

    /// The flags a child `lupin describe` is run with: `--lm` and these.
    #[must_use]
    pub fn to_argv(&self) -> Vec<String> {
        use clap::ValueEnum;
        let mut v = vec![
            "--lm".into(),
            "--lm-model".into(),
            self.lm_model.arg().into(),
        ];
        if let Some(d) = &self.lm_dir {
            v.extend(["--lm-dir".into(), d.clone()]);
        }
        if let Some(d) = self.device.to_possible_value() {
            v.extend(["--device".into(), d.get_name().to_string()]);
        }
        v
    }

    /// The model loader, run only when a text is not cached.
    fn loader(&self) -> impl FnOnce() -> Result<Box<dyn lm::Lm>> {
        let spec = self.spec();
        move || Ok(Box::new(lm::Qwen3::load(&spec)?) as Box<dyn lm::Lm>)
    }
}

#[derive(Debug, Clone)]
struct ClusterEvidence {
    id: String,
    coarse_label: String,
    best_label: String,
    best_q: Option<f32>,
    coarse_q: Option<f32>,
    label_support: Option<f32>,
    best_significant: bool,
    /// Second-lowest FDR-significant term `(label, q)`, when one exists beside the primary.
    second: Option<(String, f32)>,
    /// Marker genes/features that define the called (or closest) panel type.
    markers: Vec<String>,
    /// Extra genes pulled in via `--feature-embedding` neighbours of the markers.
    neighbour_genes: Vec<String>,
    incidence_words: Vec<String>,
}

impl ClusterEvidence {
    fn has_significant_call(&self) -> bool {
        self.best_significant && self.coarse_label != UNASSIGNED_LABEL
    }
}

/// `describe` on a trajectory run: the run's paragraph and each type's line
/// to `{out}.trajectory.describe.md`, by the model with `--lm`.
fn describe_trajectory(
    loaded: &crate::manifest::run::Loaded,
    args: &DescribeArgs,
    out: &str,
) -> Result<()> {
    let facts = trajectory::OrderFacts::read(loaded)?;
    let spec = args.lm.lm.then(|| args.lm.spec());
    let cache = consensus::cache_path(out);
    let texts = trajectory::describe(&facts, spec.as_ref(), &cache, args.lm.loader())?;
    if args.cache_only {
        return Ok(());
    }
    let md = trajectory::md_path(out);
    std::fs::write(&md, trajectory::markdown(&args.from, &texts))
        .with_context(|| format!("writing {}", md.display()))?;
    info!("wrote {} ({} text(s))", md.display(), texts.len());
    Ok(())
}

pub fn run_describe(args: &DescribeArgs) -> Result<()> {
    let Some(out) = args.out.as_deref() else {
        // The TUI writes the model's cache alone (`--cache-only`).
        if args.text_prefix.is_some() || args.markers.is_some() || args.feature_embedding.is_some()
        {
            log::warn!(
                "--text-prefix, --markers and --feature-embedding go into the description files, \
                 which only -o writes; the TUI shows the model's suggestions"
            );
        }
        return crate::tui::run_describe(&args.from, &args.lm);
    };
    legume_numeric::matrix::common_io::mkdir_parent(out)?;
    let loaded = crate::manifest::run::load(&args.from).ok();
    // A trajectory run: its order is what is described.
    if let Some(l) = loaded
        .as_ref()
        .filter(|l| crate::manifest::family::is_trajectory(&l.manifest, &l.file))
    {
        return describe_trajectory(l, args, out);
    }
    // The model's proposals come first: they need only the round.
    let proposals = if args.lm.lm {
        let loaded = loaded
            .as_ref()
            .context("--lm proposes labels for an annotated round: give its manifest")?;
        let round = consensus::inputs(loaded)?;
        let cache = consensus::cache_path(out);
        consensus::propose_round(&round, &args.lm.spec(), &cache, args.lm.loader())?
    } else {
        Vec::new()
    };
    if args.cache_only {
        return Ok(());
    }

    let files = EvidenceFiles::locate(&args.from);
    let mut enriched = match load_evidence(&files, args.fdr_alpha) {
        Ok(e) => e,
        // With proposals, the template sentences are extra.
        Err(e) if !proposals.is_empty() => {
            log::warn!("describe: no template sentences ({e:#})");
            Vec::new()
        }
        Err(e) => return Err(e),
    };
    attach_second_best(&mut enriched, &files, args.fdr_alpha)?;
    let markers_tsv = args.markers.as_deref().or(files.markers_tsv.as_deref());
    let markers_by_type = load_markers_by_type(&files, markers_tsv)?;
    attach_markers(&mut enriched, &markers_by_type);
    if let Some(emb) = args.feature_embedding.as_deref() {
        attach_embedding_neighbours(&mut enriched, emb, args.neighbour_k)?;
    }
    attach_incidence_words(&mut enriched, args.text_prefix.as_deref())?;

    let allowed = allowed_entities(&enriched);
    let mut sentences = Vec::new();
    for c in &enriched {
        let draft = template_sentence(c);
        let sentence = citation_check(&draft, &allowed, &cluster_entities(c), &own_entities(c))
            .map_err(|e| anyhow::anyhow!("citation check failed for {}: {e}", c.id))?;
        sentences.push((c.id.clone(), sentence));
    }

    let proposal_of = |id: &str| proposals.iter().find(|p| format!("K{}", p.cluster) == id);
    let mut evidence = json!({
        "from": args.from.as_ref(),
        "text_prefix": args.text_prefix.as_ref().map(|s| s.as_ref()),
        "markers": args.markers.as_ref().map(|s| s.as_ref()),
        "clusters": enriched.iter().map(|c| {
            let mut v = json!({
                "id": c.id,
                "coarse_label": c.coarse_label,
                "best_label": c.best_label,
                "best_q": c.best_q,
                "coarse_q": c.coarse_q,
                "label_support": c.label_support,
                "best_significant": c.best_significant,
                "second_label": c.second.as_ref().map(|(l, _)| l),
                "second_q": c.second.as_ref().map(|(_, q)| q),
                "markers": c.markers,
                "neighbour_genes": c.neighbour_genes,
                "incidence_words": c.incidence_words,
            });
            if let Some(p) = proposal_of(&c.id) {
                v["proposal"] = json!(p);
            }
            v
        }).collect::<Vec<_>>(),
        "allowed_entities": allowed.iter().cloned().collect::<Vec<_>>(),
    });
    if args.lm.lm {
        let spec = args.lm.spec();
        evidence["lm"] = json!({"model": spec.name(), "revision": spec.revision()});
        evidence["proposals"] = json!(proposals);
    }

    let json_path = format!("{out}.describe.json");
    let mut jw = File::create(&json_path).with_context(|| format!("create {json_path}"))?;
    serde_json::to_writer_pretty(&mut jw, &evidence)?;
    writeln!(jw)?;
    info!("wrote {json_path}");

    let md_path = format!("{out}.describe.md");
    let mut mw = File::create(&md_path).with_context(|| format!("create {md_path}"))?;
    writeln!(mw, "# Cluster descriptions")?;
    writeln!(mw)?;
    writeln!(mw, "From `{from}`.", from = args.from.as_ref())?;
    writeln!(mw)?;
    let proposed = |mw: &mut File, p: &consensus::Proposal| {
        writeln!(
            mw,
            "Proposed: **{}** ({:.2}) — {}",
            p.label,
            p.p_label(),
            p.rationale
        )
        .and_then(|()| writeln!(mw))
    };
    for (id, s) in &sentences {
        writeln!(mw, "## {id}")?;
        writeln!(mw)?;
        writeln!(mw, "{s}")?;
        writeln!(mw)?;
        if let Some(p) = proposal_of(id) {
            proposed(&mut mw, p)?;
        }
    }
    for p in &proposals {
        let id = format!("K{}", p.cluster);
        if !sentences.iter().any(|(s, _)| *s == id) {
            writeln!(mw, "## {id}")?;
            writeln!(mw)?;
            proposed(&mut mw, p)?;
        }
    }
    info!("wrote {md_path} ({} cluster(s))", sentences.len());
    Ok(())
}

/// Where each piece of annotate evidence lives.
#[derive(Default)]
struct EvidenceFiles {
    /// Per-cell table with `community` / `coarse_label` (projection).
    annot: Option<String>,
    /// Enrichment's cluster × cell-type FDR q-values.
    enrichment_q: Option<String>,
    /// Projection's cluster × term FDR q-values.
    term_q: Option<String>,
    argmax: Option<String>,
    marker_embedding: Option<String>,
    /// The marker panel annotate was run with.
    markers_tsv: Option<String>,
}

impl EvidenceFiles {
    /// From the manifest's `annotate.*` slots when `--from` is a run
    /// manifest (or its prefix) — those point at the latest pass, whatever its
    /// `-o` was. Otherwise `--from` is an annotate output prefix and
    /// the files are looked for by name.
    fn locate(from: &str) -> Self {
        match crate::manifest::run::load(from) {
            Ok(loaded) => {
                let a = &loaded.manifest.annotate;
                let at = |slot: &Option<String>| {
                    slot.as_deref()
                        .map(|rel| crate::manifest::run::resolve(&loaded.dir, rel))
                };
                let term_q = at(&a.cluster_term_q);
                Self {
                    // `annotation` holds projection's per-cell table only on a projection pass.
                    annot: term_q.is_some().then(|| at(&a.annotation)).flatten(),
                    enrichment_q: at(&a.cluster_celltype_q_values),
                    term_q,
                    argmax: at(&a.argmax),
                    marker_embedding: at(&a.marker_embedding),
                    // Manifest-relative; older manifests recorded the path as
                    // typed on the command line, so fall back to it as given.
                    markers_tsv: a
                        .markers
                        .as_deref()
                        .filter(|m| !m.is_empty())
                        .map(|m| loaded.recorded(m).unwrap_or_else(|| m.to_string())),
                }
            }
            Err(_) => {
                let found = |suffix: &str| {
                    let p = format!("{from}{suffix}");
                    Path::new(&p).exists().then_some(p)
                };
                Self {
                    annot: found(ANNOT_PARQUET),
                    enrichment_q: found(CLUSTER_CELLTYPE_Q_VALUES),
                    term_q: found(CLUSTER_TERM_Q),
                    // The parquet, else the TSV of an older pass.
                    argmax: found(CELL_LABELS_SUFFIX).or_else(|| found(ARGMAX_TSV)),
                    marker_embedding: found(".marker_embedding.parquet"),
                    markers_tsv: None,
                }
            }
        }
    }
}

fn load_evidence(files: &EvidenceFiles, fdr_alpha: f32) -> Result<Vec<ClusterEvidence>> {
    if let Some(p) = &files.annot {
        return load_from_annot_parquet(p);
    }
    if let Some(p) = &files.enrichment_q {
        return load_from_enrichment_q(p, fdr_alpha);
    }
    if let Some(p) = &files.argmax {
        return load_from_argmax(p);
    }
    bail!(
        "no annotate evidence found: need a projection annot table, an enrichment \
         q-value table, or per-cell labels (run `lupin annotate` first)"
    );
}

/// Enrichment's cluster × cell-type FDR table: each cluster's call is its
/// lowest-q cell type, significant when that q clears `fdr_alpha`.
fn load_from_enrichment_q(path: &str, fdr_alpha: f32) -> Result<Vec<ClusterEvidence>> {
    let q = read_q_table(path)?;
    Ok(q.rows
        .iter()
        .enumerate()
        .map(|(i, id)| {
            let best = (0..q.mat.ncols())
                .filter(|&j| q.mat[(i, j)].is_finite())
                .min_by(|&a, &b| q.mat[(i, a)].total_cmp(&q.mat[(i, b)]));
            let best_label = best.map_or(UNASSIGNED_LABEL, |j| q.cols[j].as_ref());
            let best_q = best.map(|j| q.mat[(i, j)]);
            let significant = best_q.is_some_and(|x| x < fdr_alpha);
            ClusterEvidence {
                id: id.to_string(),
                coarse_label: if significant {
                    best_label.to_string()
                } else {
                    UNASSIGNED_LABEL.to_string()
                },
                best_label: best_label.to_string(),
                best_q,
                coarse_q: best_q.filter(|_| significant),
                label_support: None,
                best_significant: significant,
                second: None,
                markers: Vec::new(),
                neighbour_genes: Vec::new(),
                incidence_words: Vec::new(),
            }
        })
        .collect())
}

fn read_q_table(
    path: &str,
) -> Result<
    legume_numeric::matrix::dense_mat_io::MatWithNames<legume_numeric::matrix::dense_mat_io::Mat>,
> {
    use legume_numeric::matrix::traits::IoOps;
    legume_numeric::matrix::dense_mat_io::Mat::from_parquet_with_row_names(path, Some(0))
        .with_context(|| format!("read {path}"))
}

fn load_from_annot_parquet(path: &str) -> Result<Vec<ClusterEvidence>> {
    let fields = peek_parquet_field_names(path).with_context(|| format!("peek {path}"))?;
    let has = |name: &str| fields.iter().any(|f| f.as_ref() == name);
    anyhow::ensure!(
        has("coarse_label") && has("community"),
        "{path}: need coarse_label and community columns"
    );

    let string_cols: Vec<&str> = if has("best_label") {
        vec!["coarse_label", "best_label"]
    } else {
        vec!["coarse_label"]
    };
    let mut numeric_cols: Vec<&str> = vec!["community"];
    if has("best_q") {
        numeric_cols.push("best_q");
    }
    if has("coarse_q") {
        numeric_cols.push("coarse_q");
    }
    if has("label_support") {
        numeric_cols.push("label_support");
    }
    if has("best_significant") {
        numeric_cols.push("best_significant");
    }

    let (strings, nums) =
        legume_numeric::matrix::parquet::read_table_columns(path, &string_cols, &numeric_cols)
            .with_context(|| format!("read {path}"))?;

    let coarse = &strings[0];
    let best = if strings.len() > 1 {
        &strings[1]
    } else {
        coarse
    };
    let community: Vec<i32> = nums[0].iter().map(|x| *x as i32).collect();
    let best_q_col = numeric_cols.iter().position(|&c| c == "best_q");
    let coarse_q_col = numeric_cols.iter().position(|&c| c == "coarse_q");
    let support_col = numeric_cols.iter().position(|&c| c == "label_support");
    let best_sig_col = numeric_cols.iter().position(|&c| c == "best_significant");

    let mut by_comm: BTreeMap<i32, ClusterEvidence> = BTreeMap::new();
    for i in 0..community.len() {
        let id = community[i];
        by_comm.entry(id).or_insert_with(|| {
            let coarse_s = coarse[i].to_string();
            let best_s = best[i].to_string();
            ClusterEvidence {
                id: format!("K{id}"),
                coarse_label: coarse_s.clone(),
                best_label: best_s,
                best_q: best_q_col.map(|j| nums[j][i] as f32),
                coarse_q: coarse_q_col.map(|j| nums[j][i] as f32),
                label_support: support_col.map(|j| nums[j][i] as f32),
                best_significant: best_sig_col
                    .map(|j| nums[j][i] as i32 != 0)
                    .unwrap_or(coarse_s != UNASSIGNED_LABEL),
                second: None,
                markers: Vec::new(),
                neighbour_genes: Vec::new(),
                incidence_words: Vec::new(),
            }
        });
    }
    Ok(by_comm.into_values().collect())
}

fn load_from_argmax(path: &str) -> Result<Vec<ClusterEvidence>> {
    let labels = crate::cell_labels::read_cell_labels(path)?.unique_groups();
    Ok(labels
        .into_iter()
        .enumerate()
        .map(|(i, lab)| {
            let sig = lab.as_ref() != UNASSIGNED_LABEL;
            ClusterEvidence {
                id: format!("L{i}"),
                coarse_label: lab.to_string(),
                best_label: lab.into(),
                best_q: None,
                coarse_q: None,
                label_support: None,
                best_significant: sig,
                second: None,
                markers: Vec::new(),
                neighbour_genes: Vec::new(),
                incidence_words: Vec::new(),
            }
        })
        .collect())
}

/// Load panel markers keyed by cell-type name.
///
/// Preference: projection's `marker_embedding`, then a marker TSV
/// (`--markers`, else the one annotate used).
fn load_markers_by_type(
    files: &EvidenceFiles,
    markers_tsv: Option<&str>,
) -> Result<BTreeMap<String, Vec<String>>> {
    if let Some(p) = &files.marker_embedding {
        return load_markers_from_embedding_parquet(p);
    }
    if let Some(path) = markers_tsv {
        return load_markers_from_tsv(path);
    }
    info!("describe: no marker_embedding; pass --markers for a gene<TAB>celltype TSV");
    Ok(BTreeMap::new())
}

fn load_markers_from_embedding_parquet(path: &str) -> Result<BTreeMap<String, Vec<String>>> {
    let fields = peek_parquet_field_names(path).with_context(|| format!("peek {path}"))?;
    let has = |name: &str| fields.iter().any(|f| f.as_ref() == name);
    anyhow::ensure!(has("gene") && has("type"), "{path}: need gene and type");

    let string_cols = ["gene", "type"];
    let numeric_cols: Vec<&str> = if has("weight") {
        vec!["weight"]
    } else {
        vec![]
    };
    let (strings, nums) =
        legume_numeric::matrix::parquet::read_table_columns(path, &string_cols, &numeric_cols)
            .with_context(|| format!("read {path}"))?;
    let genes = &strings[0];
    let types = &strings[1];
    let mut scored: BTreeMap<String, Vec<(i32, f32, String)>> = BTreeMap::new();
    for i in 0..genes.len() {
        let w = nums.first().map(|col| col[i] as f32).unwrap_or(1.0);
        scored
            .entry(types[i].to_string())
            .or_default()
            .push((1, w, genes[i].to_string()));
    }
    Ok(finalize_marker_lists(scored))
}

fn load_markers_from_tsv(path: &str) -> Result<BTreeMap<String, Vec<String>>> {
    let reader = open_buf_reader(path).with_context(|| format!("open {path}"))?;
    let mut scored: BTreeMap<String, Vec<(i32, f32, String)>> = BTreeMap::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split('\t');
        let Some(gene) = parts.next() else {
            continue;
        };
        let Some(ct) = parts.next() else {
            continue;
        };
        if gene.eq_ignore_ascii_case("gene") || gene.eq_ignore_ascii_case("feature") {
            continue;
        }
        scored
            .entry(ct.to_string())
            .or_default()
            .push((1, 1.0, gene.to_string()));
    }
    Ok(finalize_marker_lists(scored))
}

fn finalize_marker_lists(
    scored: BTreeMap<String, Vec<(i32, f32, String)>>,
) -> BTreeMap<String, Vec<String>> {
    let mut out = BTreeMap::new();
    for (ct, mut rows) in scored {
        rows.sort_by(|a, b| {
            b.0.cmp(&a.0)
                .then_with(|| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal))
                .then_with(|| a.2.cmp(&b.2))
        });
        let mut seen = BTreeSet::new();
        let mut genes = Vec::new();
        for (_, _, g) in rows {
            if seen.insert(g.clone()) {
                genes.push(g);
            }
            if genes.len() >= 12 {
                break;
            }
        }
        out.insert(ct, genes);
    }
    out
}

fn attach_markers(clusters: &mut [ClusterEvidence], by_type: &BTreeMap<String, Vec<String>>) {
    for c in clusters {
        let key = if c.has_significant_call() {
            c.coarse_label.as_str()
        } else if c.best_label != UNASSIGNED_LABEL {
            c.best_label.as_str()
        } else {
            continue;
        };
        if let Some(genes) = by_type.get(key) {
            c.markers = genes.iter().take(8).cloned().collect();
        }
    }
}

/// Fill `second` from the cluster × label FDR table: projection's
/// `cluster_term_q`, else enrichment's `cluster_celltype_q_values`.
///
/// Among terms with `q < fdr_alpha`, take the lowest-q term that is **not** the
/// primary call. Omitted when the runner-up is nonsignificant or the matrix is absent.
fn attach_second_best(
    clusters: &mut [ClusterEvidence],
    files: &EvidenceFiles,
    fdr_alpha: f32,
) -> Result<()> {
    let Some(path) = files.term_q.as_ref().or(files.enrichment_q.as_ref()) else {
        info!("describe: no cluster q table; skipping significant runner-up");
        return Ok(());
    };
    let loaded = read_q_table(path)?;
    let row_of: BTreeMap<&str, usize> = loaded
        .rows
        .iter()
        .enumerate()
        .map(|(i, r)| (r.as_ref(), i))
        .collect();

    for c in clusters.iter_mut() {
        if !c.has_significant_call() {
            continue;
        }
        let Some(&ri) = row_of.get(c.id.as_str()) else {
            continue;
        };
        let primary = c.coarse_label.as_str();
        let second = loaded
            .cols
            .iter()
            .enumerate()
            .map(|(j, name)| (loaded.mat[(ri, j)], name.as_ref()))
            .filter(|&(q, lab)| q.is_finite() && q < fdr_alpha && lab != primary)
            .min_by(|a, b| a.0.total_cmp(&b.0).then_with(|| a.1.cmp(b.1)));
        if let Some((q, lab)) = second {
            c.second = Some((lab.to_string(), q));
        }
    }
    Ok(())
}

fn load_incidence_words(prefix: &str) -> Result<Vec<String>> {
    let path = format!("{prefix}.vocab.tsv");
    if !Path::new(&path).exists() {
        info!("describe: no {path}; skipping keyword incidence");
        return Ok(Vec::new());
    }
    // n_docs is only used for df_frac display in Vocabulary; 1 is fine for incidence listing.
    let vocab = Vocabulary::read_tsv(&path, 1)?;
    Ok(vocab.kept.iter().take(32).map(|w| w.to_string()).collect())
}

/// `feature → [(word, weight), …]` from `{prefix}.feature_word.edges.tsv`.
/// Keys are lowercased feature names. Only `* → word` rows are kept.
fn load_feature_word_edges(prefix: &str) -> Result<BTreeMap<String, Vec<(String, f32)>>> {
    let path = format!("{prefix}.feature_word.edges.tsv");
    if !Path::new(&path).exists() {
        return Ok(BTreeMap::new());
    }
    let f = File::open(&path).with_context(|| format!("open {path}"))?;
    let mut out: BTreeMap<String, Vec<(String, f32)>> = BTreeMap::new();
    for line in BufReader::new(f).lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            continue;
        }
        // lhs_type, lhs, rhs_type, rhs, weight
        if parts[2] != "word" {
            continue;
        }
        let gene = parts[1].to_lowercase();
        let word = parts[3].to_string();
        let w: f32 = parts[4].parse().unwrap_or(0.0);
        if w > 0.0 && !word.is_empty() {
            out.entry(gene).or_default().push((word, w));
        }
    }
    Ok(out)
}

fn gene_key(g: &str) -> String {
    g.to_lowercase()
}

/// Per-cluster keywords: sum feature→word weights over that cluster's markers
/// (and embedding neighbours), ranked by total weight.
///
/// Falls back to a global vocab head only when no `feature_word.edges.tsv` exists.
fn attach_incidence_words(
    clusters: &mut [ClusterEvidence],
    text_prefix: Option<&str>,
) -> Result<()> {
    let Some(prefix) = text_prefix else {
        return Ok(());
    };
    let edges = load_feature_word_edges(prefix)?;
    if edges.is_empty() {
        let fallback = load_incidence_words(prefix)?;
        if fallback.is_empty() {
            return Ok(());
        }
        info!(
            "describe: no feature_word edges under `{prefix}`; \
             using global vocab head (not per-cluster)"
        );
        for c in clusters.iter_mut() {
            c.incidence_words = fallback.iter().take(8).cloned().collect();
        }
        return Ok(());
    }

    for c in clusters.iter_mut() {
        let mut score: BTreeMap<String, f32> = BTreeMap::new();
        let genes = c.markers.iter().chain(c.neighbour_genes.iter());
        for g in genes {
            if let Some(words) = edges.get(&gene_key(g)) {
                for (word, w) in words {
                    *score.entry(word.clone()).or_default() += *w;
                }
            }
        }
        let mut ranked: Vec<(String, f32)> = score.into_iter().collect();
        ranked.sort_by(|a, b| {
            b.1.partial_cmp(&a.1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.0.cmp(&b.0))
        });
        c.incidence_words = ranked.into_iter().take(8).map(|(w, _)| w).collect();
    }
    Ok(())
}

/// Expand each cluster's gene set by cosine neighbours of the marker centroid
/// in a feature×H embedding (prefer `feature_coembedding` when available).
fn attach_embedding_neighbours(
    clusters: &mut [ClusterEvidence],
    emb_path: &str,
    k: usize,
) -> Result<()> {
    if k == 0 {
        return Ok(());
    }
    use legume_numeric::matrix::dmatrix_io::DMatrix;
    use legume_numeric::matrix::traits::IoOps;
    let loaded = DMatrix::<f32>::from_parquet_with_row_names(emb_path, Some(0))
        .with_context(|| format!("read {emb_path}"))?;
    let (n, h) = (loaded.mat.nrows(), loaded.mat.ncols());
    if n == 0 || h == 0 {
        return Ok(());
    }

    let mut index: BTreeMap<String, usize> = BTreeMap::new();
    for (i, name) in loaded.rows.iter().enumerate() {
        index.insert(gene_key(name), i);
    }

    // Precompute L2 norms for cosine.
    let norms: Vec<f32> = (0..n)
        .map(|i| {
            let mut s = 0.0f32;
            for j in 0..h {
                let v = loaded.mat[(i, j)];
                s += v * v;
            }
            s.sqrt().max(1e-12)
        })
        .collect();

    for c in clusters.iter_mut() {
        if c.markers.is_empty() {
            continue;
        }
        let mut centroid = vec![0.0f32; h];
        let mut n_hit = 0usize;
        let mut used = BTreeSet::new();
        for g in &c.markers {
            let Some(&ri) = index.get(&gene_key(g)) else {
                continue;
            };
            n_hit += 1;
            used.insert(ri);
            for (j, slot) in centroid.iter_mut().enumerate() {
                *slot += loaded.mat[(ri, j)];
            }
        }
        if n_hit == 0 {
            continue;
        }
        for v in &mut centroid {
            *v /= n_hit as f32;
        }
        let cnorm = centroid
            .iter()
            .map(|v| v * v)
            .sum::<f32>()
            .sqrt()
            .max(1e-12);

        let mut scored: Vec<(f32, usize)> = Vec::with_capacity(n);
        for (i, &gnorm) in norms.iter().enumerate() {
            if used.contains(&i) {
                continue;
            }
            let mut dot = 0.0f32;
            for (j, &cv) in centroid.iter().enumerate() {
                dot += cv * loaded.mat[(i, j)];
            }
            scored.push((dot / (cnorm * gnorm), i));
        }
        if scored.len() > k {
            scored.select_nth_unstable_by(k, |a, b| {
                b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal)
            });
            scored.truncate(k);
        }
        scored.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
        c.neighbour_genes = scored
            .into_iter()
            .map(|(_, i)| loaded.rows[i].to_string())
            .collect();
    }
    info!("describe: expanded markers with up to {k} neighbours from {emb_path}");
    Ok(())
}

fn cluster_entities(c: &ClusterEvidence) -> BTreeSet<String> {
    let mut s = BTreeSet::new();
    if c.coarse_label != UNASSIGNED_LABEL {
        s.insert(c.coarse_label.clone());
    }
    if c.best_label != UNASSIGNED_LABEL {
        s.insert(c.best_label.clone());
    }
    if let Some((lab, _)) = &c.second {
        s.insert(lab.clone());
    }
    s.extend(c.incidence_words.iter().cloned());
    s.extend(c.markers.iter().cloned());
    // Neighbours are evidence for keyword fishing, not citation-checked panel types.
    s
}

fn allowed_entities(clusters: &[ClusterEvidence]) -> BTreeSet<String> {
    let mut s = BTreeSet::new();
    for c in clusters {
        s.extend(cluster_entities(c));
        s.extend(c.neighbour_genes.iter().cloned());
    }
    s.insert(UNASSIGNED_LABEL.into());
    s
}

fn significance_clause(c: &ClusterEvidence, omit_best_q: bool) -> String {
    let mut parts = Vec::new();
    if !omit_best_q {
        if let Some(q) = c.best_q.filter(|q| q.is_finite()) {
            parts.push(format!("best q={q:.3}"));
        }
    }
    if let Some(q) = c.coarse_q.filter(|q| q.is_finite()) {
        if c.best_q.map(|b| (b - q).abs() > 1e-6).unwrap_or(true) {
            parts.push(format!("coarse q={q:.3}"));
        }
    }
    if let Some(s) = c.label_support.filter(|s| s.is_finite()) {
        parts.push(format!("label support={s:.2}"));
    }
    if parts.is_empty() {
        String::new()
    } else {
        format!(" ({})", parts.join("; "))
    }
}

fn marker_sentence(c: &ClusterEvidence, omit_best_q: bool) -> String {
    let sig = significance_clause(c, omit_best_q);
    if c.markers.is_empty() {
        if sig.is_empty() {
            return String::new();
        }
        return if c.has_significant_call() {
            format!(" Significance:{sig}.")
        } else {
            format!(" Closest-type significance:{sig}.")
        };
    }
    let genes = c
        .markers
        .iter()
        .take(5)
        .cloned()
        .collect::<Vec<_>>()
        .join(", ");
    if c.has_significant_call() {
        format!(" The call is supported by markers {genes}{sig}.")
    } else {
        format!(" Closest-type markers include {genes}{sig}.")
    }
}

fn runner_up_clause(c: &ClusterEvidence) -> String {
    match &c.second {
        Some((lab, q)) if c.has_significant_call() => {
            format!(" A second significant contender is {lab} (q={q:.3}).")
        }
        _ => String::new(),
    }
}

fn template_sentence(c: &ClusterEvidence) -> String {
    template_led(c, None)
}

/// [`template_sentence`] with `lead` in place of its opening `{id} is
/// annotated as {label}.`, for a label the cluster does not have yet.
fn template_led(c: &ClusterEvidence, lead: Option<&str>) -> String {
    let words = if c.incidence_words.is_empty() {
        String::new()
    } else {
        format!(
            " Keywords from marker text: {}.",
            c.incidence_words
                .iter()
                .take(5)
                .cloned()
                .collect::<Vec<_>>()
                .join(", ")
        )
    };
    let runner = runner_up_clause(c);
    if c.has_significant_call() {
        let markers = marker_sentence(c, false);
        let lead = lead.map_or_else(
            || format!("{} is annotated as {}.", c.id, c.coarse_label),
            str::to_string,
        );
        format!("{lead}{markers}{runner}{words}")
    } else if c.best_label != UNASSIGNED_LABEL {
        let q = c
            .best_q
            .map(|q| format!(" (best q={q:.3})"))
            .unwrap_or_default();
        // best q already appears in the lead clause — omit it from the marker sentence.
        let markers = marker_sentence(c, c.best_q.is_some());
        format!(
            "{} has no significant call; closest panel type is {}{}.{}{}",
            c.id, c.best_label, q, markers, words
        )
    } else {
        let markers = marker_sentence(c, false);
        format!("{} is unassigned.{}{}", c.id, markers, words)
    }
}

/// Why a sentence failed the citation check.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CitationError {
    /// It names none of the entities it should name one of.
    NamesNone,
    /// It names this entity of the run, which is not its own.
    Foreign(String),
    Braces,
    /// (A model's text) a number outside the names it cites.
    Digits(String),
    /// (A model's text) a formula or markup sign.
    Markup(String),
    /// (A model's text) nothing.
    Empty,
    /// (A model's text) a claim its evidence does not bear out
    /// ([`claims`]).
    Claim(String),
}

impl std::fmt::Display for CitationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NamesNone => write!(f, "names none of its evidence"),
            Self::Foreign(e) => write!(f, "mentions `{e}` outside its evidence"),
            Self::Braces => write!(f, "braces in the sentence"),
            Self::Digits(d) => write!(f, "a number ({d})"),
            Self::Markup(m) => write!(f, "markup ({m})"),
            Self::Empty => write!(f, "no sentence"),
            Self::Claim(c) => write!(f, "a claim its evidence does not bear out: {c}"),
        }
    }
}

impl std::error::Error for CitationError {}

/// `s` as the check compares it: lowercase, the panel's underscores as
/// spaces.
pub(crate) fn norm(s: &str) -> String {
    s.to_lowercase().replace('_', " ")
}

/// `text` names `entity` (both [`norm`]ed) as a whole word or words, or
/// their plural (`cells`).
pub(crate) fn mentions(text: &str, entity: &str) -> bool {
    !entity.is_empty() && !claims::find(text, entity, true).is_empty()
}

/// Check `draft` against what it may say: it names at least one of `named`
/// (when there are any), and no entity of the `run` that is neither in
/// `named` nor in `own` (what the subject may name but need not, such as
/// its markers). Short names (under three letters) are not policed: gene
/// symbols collide with them.
pub fn citation_check(
    draft: &str,
    run: &BTreeSet<String>,
    named: &BTreeSet<String>,
    own: &BTreeSet<String>,
) -> std::result::Result<String, CitationError> {
    let text = norm(draft);
    let named_n: Vec<String> = named.iter().map(|e| norm(e)).collect();
    for ent in run {
        if ent.len() < 3
            || ent.as_str() == UNASSIGNED_LABEL
            || own.contains(ent)
            || named.contains(ent)
        {
            continue;
        }
        let e = norm(ent);
        // A name inside one of the subject's own (`B` in `Pre B`) is that one.
        if mentions(&text, &e)
            && !named_n
                .iter()
                .any(|a| a.contains(&e) || e.contains(a.as_str()))
        {
            return Err(CitationError::Foreign(ent.clone()));
        }
    }
    if !named_n.is_empty() && !named_n.iter().any(|e| mentions(&text, e)) {
        return Err(CitationError::NamesNone);
    }
    if draft.contains('{') || draft.contains('}') {
        return Err(CitationError::Braces);
    }
    Ok(draft.to_string())
}

/// What cluster `c`'s template sentence may name besides its own
/// entities: its markers, keywords and their neighbours.
fn own_entities(c: &ClusterEvidence) -> BTreeSet<String> {
    c.markers
        .iter()
        .chain(&c.incidence_words)
        .chain(&c.neighbour_genes)
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests;
