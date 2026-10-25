//! `lupin mcp`: lupin's annotation and trajectory loops as Model Context
//! Protocol tools, served over stdio with `legume-mcp`, so an agent can run
//! them: run a pass, review it, preview and apply decisions, order the
//! types and read the order back.
//!
//! lupin states the facts (calls, q, ranks, edges, verdicts, all written by
//! code); the agent brings the judgment, as `relabel` decisions with a
//! rationale or as precedence statements with a note, both kept with the
//! round or run they made. Passes, relabels, `ask` and trajectories run as
//! child `lupin` processes, as from the TUI: stdout carries the protocol,
//! and their log lines come back as progress. No tool replaces an existing
//! manifest: a new output takes the next free name, in today's results
//! folder unless the agent names one. `legume_mcp::project` keeps the
//! project's records after Noble (PLoS Comput Biol 5(7): e1000424, 2009,
//! doi:10.1371/journal.pcbi.1000424): every step that writes, with the
//! inputs made for it, goes into the day's `runall.sh`, so the work can be
//! redone without the agent.

use crate::annotate::expression::{label_marker_set, profile_and_panel, specific_genes, UP_LOG2FC};
use crate::annotate::rounds::{ClusterId, DecidedBy};
use crate::child::{lupin_command, relabel_argv};
use crate::describe::trajectory::{self as order_text, OrderFacts};
use crate::manifest::family::{file_name, Tag};
use crate::manifest::rounds::{review_json, Chain};
use crate::manifest::run::{self, annotated_path, derive_out_prefix, manifest_file};
use anyhow::{Context, Result};
use clap::ValueEnum;
use legume_mcp::{Call, Content, Handler, ServerInfo, ToolSpec};
use serde_json::{json, Map, Value};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;

const INSTRUCTIONS: &str = "\
lupin annotates single-cell clusters with cell types and orders the types \
into trajectories. Its tools state facts computed by code; you make the \
judgments, each with a written rationale that lupin keeps.

Annotation loop:
1. annotate_run on a run manifest (senna or pinto) writes a new round; or, \
without a marker panel, annotate_ask writes an unlabelled round and lists \
each cluster's specific genes for you to label. With one of `go` (the \
species' GO annotations), `gaf` or `gmt` (a pathway file), it also tests GO \
or pathway terms per cluster: evidence of what a cluster does, beside its \
cell type.
2. annotate_review shows each cluster: its calls with q, its top terms, \
each method's top type, its Cell Ontology placement and its history. \
annotate_genes adds the genes: what sets each cluster apart, and how the \
panel markers of its label and of each method's top type fare in it. Where \
the methods disagree, weigh their evidence, terms included; do not label \
beyond it. A cluster's terms rank by effect, what sets it apart (NK-cell \
killing, platelet clotting, MHC-II presentation); their q tests that effect \
against random gene sets and discounts programmes many clusters share \
(ribosome, spliceosome), so a cell type split over several clusters tests \
weaker in each than its effect reads: weigh effect and q together. A round's `caveats`, when it has any, say what its statistics \
lack: read them first.
   Ask why in contrasts: why X rather than Y, given Z. The answer weighs \
lines of evidence that rest on different data, says which side each favours, \
and holds the gene and term questions behind them. A tie or a disagreement is \
an answer: the evidence does not separate them. Before you label, relabel or \
merge, put the choice to annotate_why: why a cluster is A rather than B (left \
out: its label against its strongest rival), why a gene is shown rather than \
another, why one term rather than another. Base your rationale on its verdict \
and name the lines it rests on, and put its `cite` in the decision's \
`evidence`; where the evidence does not separate A and B, keep a label over \
both or say so.
3. annotate_preview with your decisions shows what they would change, \
rescored; then annotate_relabel applies them and returns the next round. \
For a merge, the preview lists the merged cluster's top terms with their \
effect and q, so weigh whether the parts share a function before merging. A merge that would leave one cluster is \
refused (each cluster is scored against the others): label the clusters \
the same type instead. annotate_relabel's `warnings` say what the new round \
lacks, such as terms not scored again. \
Rounds are never overwritten, so a decision can always be revisited.

Trajectory loop:
1. trajectory_run on an annotated round, with roots and any precedence \
statements (from, to, precedes|unrelated, note), writes a run.
2. trajectory_read returns its lineages, each edge's connectivity and \
verdict (supported, unsupported, candidate, inferred) and the order in words. \
Revise the statements where the data disagree and run again.

annotate_run and trajectory_run take `background`: on a large run set it, \
and the call returns a job id at once; the step runs on lupin's own screen in \
a session the user can watch with `lupin attach`, and goes on if this server \
ends.

docs prints the method write-ups. Paths are resolved against the server's \
working directory; results give the paths written.";

fn schema(properties: Value, required: &[&str]) -> Value {
    json!({ "type": "object", "properties": properties, "required": required })
}

fn tools() -> Vec<ToolSpec> {
    let options = json!({
        "type": "array", "items": { "type": "string" },
        "description": "Further command-line options, as `lupin annotate --help` lists them"
    });
    let decisions = json!({
        "type": "array",
        "description": "Decisions, one object each: {cluster or clusters, action: label|merge|keep|markers_add|markers_drop, label?, features?, evidence?: [..], alternatives?: [{label, why_not}], rationale (required)}. decided_by defaults to agent (your own decision); give agent_proposed_user_accepted only when the user accepted it.",
        "items": { "type": "object" }
    });
    let background = json!({
        "type": "boolean",
        "description": "Run as a background job: the call returns its id at once; follow it with job_wait or job_status, and use the manifest only once it is done"
    });
    vec![
        ToolSpec::new(
            "annotate_run",
            "Annotate a run's clusters with cell types (`lupin annotate`) and write a new round. Returns the round's path and each cluster's size, label and each method's top type. A first round's labels are coarse Cell Ontology groups unless `options` holds `--fine`.",
            schema(
                json!({
                    "from": { "type": "string", "description": "Run manifest or its prefix (or a round, to annotate it again)" },
                    "markers": { "type": "string", "description": "Marker panel (cell type, gene per line); a round's own panel is used when omitted" },
                    "method": { "type": "string", "enum": ["auto", "enrichment", "projection"] },
                    "go": { "type": "boolean", "description": "Also score GO terms per cluster, with the GO annotations of the run's species (human or mouse, from its gene names), downloaded once (one of go, gaf, gmt)" },
                    "gaf": { "type": "string", "description": "Also score GO terms from this GO annotation file (.gaf; one of go, gaf, gmt)" },
                    "gmt": { "type": "string", "description": "Also score these gene sets or pathways (.gmt: name, description, genes per line; one of go, gaf, gmt)" },
                    "out": { "type": "string", "description": "Output prefix; the next free {stem}.L{k} when omitted" },
                    "options": options,
                    "background": background.clone(),
                }),
                &["from"],
            ),
        ),
        ToolSpec::new(
            "annotate_review",
            "Show an annotation round cluster by cluster, as `lupin review --json`: digest (size, label, calls with q, top GO/GMT terms with effect and q, evidence check), history, each method's top type and Cell Ontology placement; with `types`, every type's evidence too (large: ask for it on a few clusters). When the round has caveats (e.g. merged clusters whose terms were not scored again), a second item lists them. For why one label, gene or term rather than another, use annotate_why.",
            schema(
                json!({
                    "round": { "type": "string", "description": "Annotated manifest (a round) or its prefix" },
                    "clusters": { "type": "array", "items": { "type": "integer" }, "description": "Only these clusters" },
                    "types": { "type": "boolean", "description": "Add every type's evidence per cluster" },
                }),
                &["round"],
            ),
        ),
        ToolSpec::new(
            "annotate_genes",
            "Per cluster of a round, the genes behind its calls: the genes that set it apart (log2 fold change over the other clusters), and for its label, each method's top type and any `types` given, those types' panel markers: how many, how many are up (log2 fold change ≥ 1) and the `top` highest with their log2 fold change. A merged cluster's genes are its parts' together (their cell-weighted mean profile). A coarse label (a Cell Ontology group of the panel's types) is measured on its members' markers, listed in `group_of`. For why one label, gene or term rather than another, use annotate_why.",
            schema(
                json!({
                    "round": { "type": "string" },
                    "clusters": { "type": "array", "items": { "type": "integer" }, "description": "Only these clusters" },
                    "types": { "type": "array", "items": { "type": "string" }, "description": "More cell types whose markers to show" },
                    "top": { "type": "integer", "description": "Genes listed per cluster and per type (20)" },
                }),
                &["round"],
            ),
        ),
        ToolSpec::new(
            "annotate_why",
            "Ask why a round shows X rather than Y, given Z (`lupin why --json`). A label: why cluster `cluster` is `label` rather than `over` (both left out: its label against its strongest rival by the round's method; a coarse label with no rival outside it is asked between its own types). Returns each line of evidence (the panel: each method's order, a test of the markers only one type has, the ontology walk's placement, counted as one line; the GO/GMT terms that still stand out once the panel's genes are taken out of them), which side each favours, the gene and term questions behind them, and the verdict over the lines. A gene: why `gene` rather than `over` in `list`, with each one's place, or why it is not shown. A term: why `term` rather than `over`, by effect, with z, q, the genes driving each and the clusters that share it. `cite` is what a decision's `evidence` quotes to keep the answer.",
            schema(
                json!({
                    "round": { "type": "string" },
                    "cluster": { "type": "integer" },
                    "label": { "type": "string", "description": "Why this label (one of label, gene, term; none: the cluster's label)" },
                    "gene": { "type": "string", "description": "Why this gene, in `list`" },
                    "term": { "type": "string", "description": "Why this GO/GMT term (an id or a name)" },
                    "over": { "type": "string", "description": "Rather than this; left out, the strongest rival" },
                    "list": { "type": "string", "description": "For a gene: specific (default), markers:<type> or term:<id or name>" },
                    "depth": { "type": "integer", "description": "For a label: levels of gene and term questions under its lines (2)" },
                }),
                &["round", "cluster"],
            ),
        ),
        ToolSpec::new(
            "annotate_preview",
            "What decisions would change in a round, rescored, without writing anything (`lupin relabel --preview`). A merge also lists the merged cluster's top GO/GMT terms, with effect and q, when the round scored terms.",
            schema(
                json!({ "round": { "type": "string" }, "decisions": decisions.clone() }),
                &["round", "decisions"],
            ),
        ),
        ToolSpec::new(
            "annotate_relabel",
            "Apply decisions to a round and write the next round of its chain (`lupin relabel --next`); the rationales are kept in its history. Returns the new round's path, and `warnings` when it lacks something (e.g. a merge whose GO/GMT terms could not be scored again).",
            schema(
                json!({ "round": { "type": "string" }, "decisions": decisions }),
                &["round", "decisions"],
            ),
        ),
        ToolSpec::new(
            "annotate_ask",
            "Without a marker panel: write an unlabelled first round and return each cluster's most specific genes (`lupin ask`). Label the clusters with annotate_relabel on the round it names.",
            schema(
                json!({
                    "from": { "type": "string", "description": "Run manifest or its prefix" },
                    "out": { "type": "string", "description": "Output prefix; the next free {stem}.L{k} when omitted" },
                    "context": { "type": "string", "description": "What the cells are: tissue, species, condition" },
                    "top": { "type": "integer", "description": "Specific genes listed per cluster (15)" },
                }),
                &["from"],
            ),
        ),
        ToolSpec::new(
            "trajectory_run",
            "Order an annotated round's cells along a prior over their cell types (`lupin trajectory`): Cell Ontology and precedence statements, checked against the kNN graph, then diffusion pseudotime. Returns what trajectory_read returns for the new run.",
            schema(
                json!({
                    "from": { "type": "string", "description": "Annotated round (or run) manifest or its prefix" },
                    "out": { "type": "string", "description": "Output prefix; the next free {stem}.T{k} when omitted" },
                    "roots": { "type": "array", "items": { "type": "string" }, "description": "Cell types the trajectory starts from" },
                    "statements": {
                        "type": "array",
                        "description": "Precedence statements for this run, written to {out}.precedence.tsv",
                        "items": {
                            "type": "object",
                            "properties": {
                                "from": { "type": "string" },
                                "to": { "type": "string" },
                                "relation": { "type": "string", "enum": ["precedes", "unrelated"] },
                                "note": { "type": "string", "description": "Why; kept with the run" },
                            },
                            "required": ["from", "to", "relation", "note"],
                        },
                    },
                    "prior": { "type": "string", "description": "A precedence file, or a run's {out}.trajectory_prior.tsv to replay (not with statements)" },
                    "options": { "type": "array", "items": { "type": "string" }, "description": "Further options, as `lupin trajectory --help` lists them" },
                    "background": background,
                }),
                &["from"],
            ),
        ),
        ToolSpec::new(
            "trajectory_read",
            "Read a trajectory run: its lineages (root to leaf), its roots, branch points and ends, each type pair's connectivity, place in the prior, verdict and order agreement, and the order in words.",
            schema(
                json!({ "run": { "type": "string", "description": "Trajectory manifest or its prefix" } }),
                &["run"],
            ),
        ),
        ToolSpec::new(
            "docs",
            "lupin's method write-ups; without a topic, the list of topics.",
            schema(
                json!({ "topic": { "type": "string", "enum": crate::docs::slugs() } }),
                &[],
            ),
        ),
    ]
}

/// `lupin mcp watch`: follow the server in `dir` (else the one running, or
/// one picked among several), read-only: a person sees what the agent
/// does, and stops nothing.
pub fn watch(dir: Option<PathBuf>) -> Result<()> {
    // The log stays off the screen while it is up.
    let _quiet = legume_tui::logs::quiet_stderr();
    legume_mcp::watch::view::run(&crate::tui::TOOL, dir, false, |_| None)
}

pub fn run_mcp() -> Result<()> {
    let info = ServerInfo {
        name: "lupin".into(),
        version: env!("LUPIN_VERSION").into(),
        instructions: Some(INSTRUCTIONS.into()),
    };
    legume_mcp::serve(info, Lupin)
}

struct Lupin;

impl Handler for Lupin {
    fn tools(&self) -> Vec<ToolSpec> {
        tools()
    }

    fn call(&mut self, name: &str, args: Value, cx: &Call) -> Result<Vec<Content>> {
        let a = Arguments(args);
        match name {
            "annotate_run" => annotate_run(&a, cx),
            "annotate_review" => {
                let clusters: Vec<ClusterId> = a.parse("clusters")?.unwrap_or_default();
                let with_types = a.parse::<bool>("types")?.unwrap_or(false);
                let loaded = run::load(&a.req_path("round")?)?;
                let review = review_json(&loaded, &clusters, with_types)?;
                let caveats = crate::manifest::rounds::caveats(&loaded);
                let mut out = vec![Content::Json(json!(review))];
                if !caveats.is_empty() {
                    out.push(Content::Json(json!({ "caveats": caveats })));
                }
                Ok(out)
            }
            "annotate_genes" => annotate_genes(&a),
            "annotate_preview" => Ok(vec![Content::Json(relabel_preview(&a, cx)?)]),
            "annotate_relabel" => {
                let (round, warnings) = relabel_apply(&a, cx)?;
                let mut out = json!({ "round": absolute(&round) });
                if !warnings.is_empty() {
                    out["warnings"] = json!(warnings);
                }
                Ok(vec![Content::Json(out)])
            }
            "annotate_ask" => annotate_ask(&a, cx),
            "annotate_why" => annotate_why(&a),
            "trajectory_run" => trajectory_run(&a, cx),
            "trajectory_read" => Ok(vec![Content::Json(trajectory_read(&a.req_path("run")?)?)]),
            "docs" => Ok(vec![Content::Text(match a.str("topic") {
                Some(t) => crate::docs::text(
                    crate::docs::Topic::from_str(t, true).map_err(|e| anyhow::anyhow!(e))?,
                )
                .to_string(),
                None => crate::docs::listing(),
            })]),
            _ => anyhow::bail!("no tool {name}"),
        }
    }
}

/// A tool's arguments.
struct Arguments(Value);

impl Arguments {
    /// A text argument; an empty one counts as not given.
    fn str(&self, key: &str) -> Option<&str> {
        self.0
            .get(key)
            .and_then(Value::as_str)
            .filter(|s| !s.is_empty())
    }

    /// A path argument, made absolute against the server's working
    /// directory: children and the records then mean the same file wherever
    /// the project's root is.
    fn path(&self, key: &str) -> Option<String> {
        self.str(key).map(|p| absolute(Path::new(p)))
    }

    fn req_path(&self, key: &str) -> Result<String> {
        self.path(key)
            .with_context(|| format!("`{key}` is required"))
    }

    fn parse<T: serde::de::DeserializeOwned>(&self, key: &str) -> Result<Option<T>> {
        match self.0.get(key) {
            None | Some(Value::Null) => Ok(None),
            Some(v) => serde_json::from_value(v.clone())
                .map(Some)
                .with_context(|| format!("reading `{key}`")),
        }
    }
}

/// `lupin` with `argv` as a child, refusing `--overwrite` anywhere in it:
/// no tool replaces a manifest, whatever the agent passes on. The command
/// line is parsed as lupin parses it first, so a bad one is refused in a
/// line, before anything runs.
fn lupin(argv: &[String]) -> Result<Command> {
    anyhow::ensure!(
        !argv.iter().any(|o| o == "--overwrite"),
        "--overwrite is not taken here: choose another `out`"
    );
    if let Err(e) = <crate::Cli as clap::Parser>::try_parse_from(
        std::iter::once("lupin").chain(argv.iter().map(String::as_str)),
    ) {
        use clap::error::ErrorKind::{DisplayHelp, DisplayVersion};
        anyhow::ensure!(
            !matches!(e.kind(), DisplayHelp | DisplayVersion),
            "--help and --version are not taken here: the `docs` tool describes lupin"
        );
        // The error and its tips, without the usage line and the pointer to --help.
        let text = e.to_string();
        let said: Vec<&str> = text
            .lines()
            .take_while(|l| !l.starts_with("Usage:"))
            .filter(|l| !l.trim().is_empty() && !l.starts_with("For more information"))
            .collect();
        anyhow::bail!("{}", said.join("\n").trim_start_matches("error: "));
    }
    lupin_command(argv)
}

fn absolute(p: &Path) -> String {
    std::path::absolute(p)
        .unwrap_or_else(|_| p.to_path_buf())
        .to_string_lossy()
        .into_owned()
}

/// Where a run writes: its output prefix and manifest, and whether the
/// prefix was chosen here rather than by the agent.
struct Output {
    out: String,
    manifest: PathBuf,
    chosen: bool,
}

/// `out`, or the next free `tag` name in today's results folder (beside the
/// manifest `from` names when the server keeps no records); its manifest
/// must not exist yet, nor be one a background job is still writing.
/// Checked here, before the tool writes anything beside it.
fn output(a: &Arguments, cx: &Call, from: &str, tag: Tag) -> Result<Output> {
    let source = manifest_file(from);
    let (out, chosen) = match a.path("out") {
        Some(o) => (o, false),
        None => {
            let beside = match (cx.results_dir(), source.file_name()) {
                (Some(dir), Some(name)) => dir.join(name),
                _ => source.clone(),
            };
            let held = |m: &Path| cx.ensure_free(&[m.to_path_buf()]).is_err();
            (tag.next_free_but(&beside, held), true)
        }
    };
    let manifest = annotated_path(&source, &out);
    anyhow::ensure!(
        !manifest.exists(),
        "{} exists; not replacing it: choose another `out`, or omit it for the next free name",
        manifest.display()
    );
    cx.ensure_free(std::slice::from_ref(&manifest))?;
    if let Some(dir) = Path::new(&out).parent() {
        std::fs::create_dir_all(dir).with_context(|| format!("making {}", dir.display()))?;
    }
    Ok(Output {
        out,
        manifest,
        chosen,
    })
}

/// `--{flag}={value}`: one word, so a value starting with `-` stays a value.
fn flag(flag: &str, value: &str) -> String {
    format!("--{flag}={value}")
}

/// Whether the agent asked for the step as a background job.
fn background(a: &Arguments) -> Result<bool> {
    Ok(a.parse::<bool>("background")?.unwrap_or(false))
}

/// `record` (`lupin {argv}`: an annotate pass or a trajectory, `argv[0]`)
/// in the background, writing `manifest` into `out`: on lupin's own screen
/// in a session of its own (`--job`), which the user watches with `lupin
/// attach` and which outlives this server; where there are no sessions, a
/// job of this server's. What it returns: the job's id, its session, and
/// the manifest (`key`), to use once the job is done.
fn background_job(
    cx: &Call,
    argv: &[String],
    record: Command,
    out: &str,
    manifest: &Path,
    summary: &str,
    key: &str,
) -> Result<Value> {
    #[cfg(unix)]
    {
        let mut words = vec![argv[0].clone(), "--job".into()];
        words.extend_from_slice(&argv[1..]);
        let started = cx.start_session_job(legume_mcp::SessionJob {
            tool: crate::tui::TOOL,
            screen: argv[0].clone(),
            label: Path::new(out)
                .file_name()
                .map(|n| n.to_string_lossy().into_owned()),
            watch: lupin_command(&words)?,
            record,
            outputs: vec![manifest.to_path_buf()],
            summary: summary.into(),
        })?;
        let session = started.session;
        Ok(json!({
            "job": started.job,
            key: absolute(manifest),
            "session": session,
            "note": format!(
                "running as a background job on lupin's screen in session {session} (`lupin attach {session}` shows it to the user), which goes on if this server ends: job_wait or job_status until it is done, then use it"
            ),
        }))
    }
    #[cfg(not(unix))]
    {
        let _ = (argv, out);
        let job = cx.spawn_job(
            record,
            std::slice::from_ref(&manifest.to_path_buf()),
            summary,
        )?;
        Ok(json!({
            "job": job,
            key: absolute(manifest),
            "note": "running as a background job: job_wait or job_status until it is done, then use it",
        }))
    }
}

fn annotate_run(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let from = a.req_path("from")?;
    let Output {
        out,
        manifest: round,
        ..
    } = output(a, cx, &from, Tag::Round)?;
    let mut argv = vec!["annotate".into(), flag("from", &from), flag("out", &out)];
    let mut summary = format!("annotate {}", file_name(Path::new(&from)));
    if let Some(m) = a.path("markers") {
        argv.push(flag("markers", &m));
        summary += &format!(" with {}", file_name(Path::new(&m)));
    }
    if let Some(m) = a.str("method") {
        argv.push(flag("method", m));
        summary += &format!(" by {m}");
    }
    if a.parse::<bool>("go")?.unwrap_or(false) {
        argv.push("--go".into());
        summary += " and GO terms";
    }
    for key in ["gaf", "gmt"] {
        if let Some(p) = a.path(key) {
            argv.push(flag(key, &p));
            summary += &format!(" and {}", file_name(Path::new(&p)));
        }
    }
    argv.extend(a.parse::<Vec<String>>("options")?.unwrap_or_default());
    if background(a)? {
        let job = background_job(cx, &argv, lupin(&argv)?, &out, &round, &summary, "round")?;
        return Ok(vec![Content::Json(job)]);
    }
    cx.run_making(lupin(&argv)?, std::slice::from_ref(&round), &summary)?;
    cx.state(&json!({ "round": absolute(&round) }));
    let loaded = run::load(&round.to_string_lossy())?;
    let clusters: Map<String, Value> = review_json(&loaded, &[], false)?
        .into_iter()
        .map(|(id, v)| {
            let d = &v["digest"];
            let brief = json!({ "size": d["size"], "label": d["label"], "tops": v["tops"] });
            (id.to_string(), brief)
        })
        .collect();
    Ok(vec![Content::Json(
        json!({ "round": absolute(&round), "clusters": clusters }),
    )])
}

fn annotate_why(a: &Arguments) -> Result<Vec<Content>> {
    use crate::manifest::why;
    let named = ["label", "gene", "term"]
        .iter()
        .filter(|k| a.str(k).is_some())
        .count();
    anyhow::ensure!(named <= 1, "name one of `label`, `gene` and `term`");
    let named = |k: &str| a.str(k).map(str::to_string);
    let mut args = why::WhyArgs::new(
        &a.req_path("round")?,
        a.parse("cluster")?.context("`cluster` is required")?,
    );
    (args.label, args.gene, args.term, args.over) =
        (named("label"), named("gene"), named("term"), named("over"));
    if let Some(l) = named("list") {
        args.list = l;
    }
    if let Some(d) = a.parse("depth")? {
        args.depth = d;
    }
    let answer = why::ask(&args)?;
    let mut v = serde_json::to_value(&answer)?;
    v["cite"] = why::citation(&answer);
    Ok(vec![Content::Json(v)])
}

fn annotate_genes(a: &Arguments) -> Result<Vec<Content>> {
    let loaded = run::load(&a.req_path("round")?)?;
    let clusters: Vec<ClusterId> = a.parse("clusters")?.unwrap_or_default();
    let more: Vec<String> = a.parse("types")?.unwrap_or_default();
    let top = a.parse::<usize>("top")?.unwrap_or(20);
    let cells = crate::manifest::rounds::read_cells(&loaded)?;
    let (expr, panel) = profile_and_panel(&loaded, &cells)?;
    let expr = expr.context("the round records no expression profile")?;
    let specific = specific_genes(&expr);
    // A coarse label is a Cell Ontology group of the panel's types: its
    // markers are its members'.
    let tree = crate::manifest::rounds::read_tree(&loaded).ok().flatten();
    let out: Map<String, Value> = review_json(&loaded, &clusters, false)?
        .into_iter()
        .map(|(id, v)| {
            let label = v["digest"]["label"].as_str().map(str::to_string);
            let tops = v["tops"].as_object().into_iter().flat_map(|m| m.values());
            let types: BTreeSet<String> = label
                .into_iter()
                .chain(tops.filter_map(|t| t.as_str().map(str::to_string)))
                .chain(more.iter().cloned())
                .collect();
            let markers: Map<String, Value> = types
                .into_iter()
                .map(|t| {
                    let group = tree.as_ref().and_then(|tr| tr.coarse_members(&t, &panel));
                    let set = label_marker_set(tree.as_ref(), &panel, &t);
                    let mut genes: Vec<(&String, f32)> = set
                        .iter()
                        .filter_map(|g| Some((g, expr.fold_change(id, g)?)))
                        .collect();
                    genes.sort_by(|x, y| y.1.total_cmp(&x.1));
                    let up = genes.iter().filter(|(_, f)| *f >= UP_LOG2FC).count();
                    let listed: Vec<Value> =
                        genes.iter().take(top).map(|(g, f)| json!([g, f])).collect();
                    let mut entry = json!({
                        "markers": set.len(), "measured": genes.len(), "up": up, "top": listed,
                    });
                    if let Some(m) = group {
                        entry["group_of"] = json!(m);
                    }
                    (t, entry)
                })
                .collect();
            let specific: Vec<Value> = specific
                .get(&id)
                .into_iter()
                .flatten()
                .take(top)
                .map(|(g, f)| json!([g, f]))
                .collect();
            let entry = json!({ "specific": specific, "markers": markers });
            (id.to_string(), entry)
        })
        .collect();
    Ok(vec![Content::Json(Value::Object(out))])
}

/// The round `a` names, which must exist.
fn round_arg(a: &Arguments) -> Result<PathBuf> {
    let round = manifest_file(&a.req_path("round")?);
    anyhow::ensure!(round.is_file(), "{}: no such round", round.display());
    Ok(round)
}

/// The decisions in `a`, one JSON line each, for a decisions file beside
/// `round`: each names the round (by its file name, so relative to that
/// file, whatever the agent wrote there), so it cannot land on renumbered
/// ids, and is the agent's unless it says otherwise.
fn decisions_text(a: &Arguments, round: &Path) -> Result<String> {
    let decisions: Vec<Map<String, Value>> = a.parse("decisions")?.unwrap_or_default();
    anyhow::ensure!(!decisions.is_empty(), "no decisions");
    let mut text = String::new();
    for mut d in decisions {
        d.insert("round".into(), json!(file_name(round)));
        d.entry("decided_by")
            .or_insert_with(|| json!(DecidedBy::Agent.as_str()));
        text += &format!("{}\n", Value::Object(d));
    }
    Ok(text)
}

/// `lupin relabel --preview` with the decisions in a temporary file beside
/// the round; the preview's JSON. A preview writes nothing, so it stays out
/// of the driver script.
fn relabel_preview(a: &Arguments, cx: &Call) -> Result<Value> {
    let round = round_arg(a)?;
    let text = decisions_text(a, &round)?;
    let tmp = tempfile::Builder::new()
        .prefix(".mcp-preview-")
        .suffix(".jsonl")
        .tempfile_in(run::parent_dir(&round))?;
    std::fs::write(tmp.path(), text)?;
    let argv = relabel_argv(&round, tmp.path(), &["--preview"]);
    let out = cx.run_reading(lupin(&argv)?)?;
    serde_json::from_str(&out).context("reading the preview")
}

/// Apply the decisions as the next round of the round's chain: `--next`,
/// told the round it must write (`--expect`), so a replay writes the same
/// one. The decisions go to a fresh `{round}.mcp_decisions.{k}.jsonl` beside
/// the round, and into the driver script. Also returns what the agent should
/// know of the new round that the child's log alone would say.
fn relabel_apply(a: &Arguments, cx: &Call) -> Result<(PathBuf, Vec<String>)> {
    let round = round_arg(a)?;
    let text = decisions_text(a, &round)?;
    // Refused here before anything is written; relabel checks again under
    // the chain's lock.
    let chain = Chain::of(&round);
    let (_, latest) = chain.latest();
    anyhow::ensure!(
        run::same_file(&round, &latest),
        "{} is not the latest round of its chain ({} is): review that one and decide again",
        round.display(),
        latest.display()
    );
    let file = claim_decisions_file(&round)?;
    cx.write_file(&file, &text)?;
    let prefix = chain.next_prefix();
    let next = annotated_path(&round, &prefix);
    let argv = relabel_argv(&round, &file, &["--next", "--expect", &prefix]);
    let summary = format!(
        "relabel {} with {} decision(s)",
        file_name(&round),
        text.lines().count()
    );
    cx.run_making(lupin(&argv)?, std::slice::from_ref(&next), &summary)?;
    cx.state(&json!({ "round": absolute(&next), "from": absolute(&round) }));
    // The round is written: what it says of itself is extra, not the result.
    let warnings = run::load(&next.to_string_lossy())
        .map(|l| crate::manifest::rounds::caveats(&l))
        .unwrap_or_default();
    Ok((next, warnings))
}

/// The next `{round}.mcp_decisions.{k}.jsonl` beside `round`, made empty so
/// two writers cannot take the same one.
fn claim_decisions_file(round: &Path) -> Result<PathBuf> {
    let stem = derive_out_prefix(&round.to_string_lossy());
    for k in 1.. {
        let file = PathBuf::from(format!("{stem}.mcp_decisions.{k}.jsonl"));
        match std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&file)
        {
            Ok(_) => return Ok(file),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(e).with_context(|| format!("making {}", file.display())),
        }
    }
    unreachable!("some k is free")
}

fn annotate_ask(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let from = a.req_path("from")?;
    let Output {
        out,
        manifest: round,
        ..
    } = output(a, cx, &from, Tag::Round)?;
    let mut argv = vec!["ask".into(), flag("from", &from), flag("out", &out)];
    if let Some(c) = a.str("context") {
        argv.push(flag("context", c));
    }
    if let Some(n) = a.parse::<usize>("top")? {
        argv.push(flag("top", &n.to_string()));
    }
    let summary = format!(
        "unlabelled first round of {} for ask",
        file_name(Path::new(&from))
    );
    let prompt = cx.run_making(lupin(&argv)?, std::slice::from_ref(&round), &summary)?;
    cx.state(&json!({ "round": absolute(&round) }));
    Ok(vec![
        Content::Json(json!({ "round": absolute(&round) })),
        Content::Text(prompt),
    ])
}

/// A precedence statement as `trajectory_run` takes it.
#[derive(serde::Deserialize)]
struct Statement {
    from: String,
    to: String,
    relation: String,
    note: String,
}

fn trajectory_run(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let from = a.req_path("from")?;
    let Output {
        out,
        manifest,
        chosen,
    } = output(a, cx, &from, Tag::Trajectory)?;
    let mut argv = vec!["trajectory".into(), flag("from", &from), flag("out", &out)];
    let mut summary = format!("trajectory of {}", file_name(Path::new(&from)));
    let roots: Vec<String> = a.parse("roots")?.unwrap_or_default();
    if !roots.is_empty() {
        summary += &format!(" from {}", roots.join(", "));
    }
    argv.extend(roots.iter().map(|r| flag("root", r)));
    let statements: Vec<Statement> = a.parse("statements")?.unwrap_or_default();
    let mut written = None;
    match (a.path("prior"), statements.is_empty()) {
        (Some(_), false) => anyhow::bail!("give `statements` or `prior`, not both"),
        (Some(p), true) => argv.push(flag("prior", &p)),
        (None, false) => {
            let path = format!("{out}.precedence.tsv");
            let mut text = String::new();
            for s in &statements {
                anyhow::ensure!(
                    matches!(s.relation.as_str(), "precedes" | "unrelated"),
                    "a relation is `precedes` or `unrelated`, not {:?}",
                    s.relation
                );
                for field in [&s.from, &s.to, &s.note] {
                    anyhow::ensure!(
                        !field.contains(['\t', '\n']),
                        "statement fields cannot hold tabs or newlines: {field:?}"
                    );
                }
                text += &format!("{}\t{}\t{}\t{}\n", s.from, s.to, s.relation, s.note);
            }
            // A name chosen here has no run yet, so a file there is a failed
            // attempt's and is replaced; at a name the agent gave, any other
            // file is not.
            if let (false, Ok(old)) = (chosen, std::fs::read_to_string(&path)) {
                anyhow::ensure!(
                    old == text,
                    "{path} exists; not replacing it: choose another `out`"
                );
            }
            argv.push(flag("prior", &path));
            summary += &format!(", {} statement(s)", statements.len());
            written = Some((path, text));
        }
        (None, true) => {}
    }
    argv.extend(a.parse::<Vec<String>>("options")?.unwrap_or_default());
    // Every argument is read and checked before the statements are written,
    // so a call refused writes nothing.
    let (cmd, in_background) = (lupin(&argv)?, background(a)?);
    if let Some((path, text)) = written {
        cx.write_file(Path::new(&path), &text)?;
    }
    if in_background {
        let job = background_job(cx, &argv, cmd, &out, &manifest, &summary, "run")?;
        return Ok(vec![Content::Json(job)]);
    }
    cx.run_making(cmd, std::slice::from_ref(&manifest), &summary)?;
    cx.state(&json!({ "trajectory": absolute(&manifest), "from": from }));
    Ok(vec![Content::Json(trajectory_read(
        &manifest.to_string_lossy(),
    )?)])
}

fn trajectory_read(path: &str) -> Result<Value> {
    let loaded = run::load(path)?;
    let facts = OrderFacts::read(&loaded)?;
    anyhow::ensure!(
        !facts.lineages.is_empty() || !facts.edges.is_empty(),
        "{} records no trajectory",
        loaded.file.display()
    );
    let lineages: Vec<Value> = facts
        .lineages
        .iter()
        .map(|l| json!({ "name": l.name, "component": l.component, "path": l.path }))
        .collect();
    let edges: Vec<Value> = facts
        .edges
        .iter()
        .map(|e| {
            json!({
                "a": e.a, "b": e.b,
                "connectivity": e.connectivity,
                "in_prior": e.in_prior,
                "verdict": e.verdict.map(|v| v.as_str()),
                "order_agreement": e.order_agreement,
            })
        })
        .collect();
    let texts: Vec<Value> = order_text::texts(&facts)
        .into_iter()
        .map(|t| json!({ "type": t.subject, "text": t.text }))
        .collect();
    Ok(json!({
        "run": absolute(&loaded.file),
        "roots": facts.roots(),
        "branch_points": facts.branch_points(),
        "ends": facts.ends(),
        "lineages": lineages,
        "edges": edges,
        "texts": texts,
    }))
}

#[cfg(test)]
#[path = "mcp_tests.rs"]
mod tests;
