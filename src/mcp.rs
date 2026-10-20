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

use crate::annotate::markers::label_key;
use crate::annotate::rounds::ClusterId;
use crate::describe::trajectory::{self as order_text, OrderFacts};
use crate::manifest::family::Tag;
use crate::manifest::rounds::{review_json, Chain};
use crate::manifest::run::{
    self, annotated_path, derive_out_prefix, manifest_file, resolve, same_file,
};
use crate::tui::round::{panel_sets, specific_genes, Expression};
use anyhow::{Context, Result};
use clap::{Args, ValueEnum};
use legume_mcp::{Call, Content, Handler, ServerInfo, ToolSpec};
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::traits::IoOps;
use serde_json::{json, Map, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Args, Debug)]
pub struct McpArgs {}

const INSTRUCTIONS: &str = "\
lupin annotates single-cell clusters with cell types and orders the types \
into trajectories. Its tools state facts computed by code; you make the \
judgments, each with a written rationale that lupin keeps.

Annotation loop:
1. annotate_run on a run manifest (senna or pinto) writes a new round; or, \
without a marker panel, annotate_ask writes an unlabelled round and lists \
each cluster's specific genes for you to label.
2. annotate_review shows each cluster: its calls with q, its top terms, \
each method's top type, its Cell Ontology placement and its history. \
annotate_genes adds the genes: what sets each cluster apart, and how the \
panel markers of its label and of each method's top type fare in it. Where \
the methods disagree, weigh their evidence; do not label beyond it.
3. annotate_preview with your decisions shows what they would change, \
rescored; then annotate_relabel applies them and returns the next round. \
Rounds are never overwritten, so a decision can always be revisited.

Trajectory loop:
1. trajectory_run on an annotated round, with roots and any precedence \
statements (from, to, precedes|unrelated, note), writes a run.
2. trajectory_read returns its lineages, each edge's connectivity and \
verdict (supported, unsupported, candidate, inferred) and the order in words. \
Revise the statements where the data disagree and run again.

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
        "description": "Decisions, one object each: {cluster or clusters, action: label|merge|keep|markers_add|markers_drop, label?, features?, evidence?: [..], alternatives?: [{label, why_not}], rationale (required)}. decided_by defaults to agent_proposed_user_accepted.",
        "items": { "type": "object" }
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
                    "out": { "type": "string", "description": "Output prefix; the next free {stem}.L{k} when omitted" },
                    "options": options,
                }),
                &["from"],
            ),
        ),
        ToolSpec::new(
            "annotate_review",
            "Show an annotation round cluster by cluster, as `lupin review --json`: digest (size, label, calls with q, terms, evidence check), history, each method's top type and Cell Ontology placement; with `types`, every type's evidence too (large: ask for it on a few clusters).",
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
            "Per cluster of a round, the genes behind its calls: the genes that set it apart (log2 fold change over the other clusters), and for its label, each method's top type and any `types` given, those types' panel markers: how many, how many are up (log2 fold change ≥ 1) and the `top` highest with their log2 fold change. Clusters merged since the pass list no genes.",
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
            "annotate_preview",
            "What decisions would change in a round, rescored, without writing anything (`lupin relabel --preview`).",
            schema(
                json!({ "round": { "type": "string" }, "decisions": decisions.clone() }),
                &["round", "decisions"],
            ),
        ),
        ToolSpec::new(
            "annotate_relabel",
            "Apply decisions to a round and write the next round of its chain (`lupin relabel --next`); the rationales are kept in its history. Returns the new round's path.",
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
            schema(json!({ "topic": { "type": "string" } }), &[]),
        ),
    ]
}

pub fn run_mcp(_args: &McpArgs) -> Result<()> {
    let info = ServerInfo {
        name: "lupin".into(),
        version: env!("CARGO_PKG_VERSION").into(),
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
                let every_type = a.parse::<bool>("types")?.unwrap_or(false);
                let loaded = run::load(&a.req("round")?)?;
                let mut review = review_json(&loaded, &clusters)?;
                if !every_type {
                    for v in review.values_mut() {
                        if let Some(m) = v.as_object_mut() {
                            m.remove("types");
                        }
                    }
                }
                Ok(vec![Content::Json(json!(review))])
            }
            "annotate_genes" => annotate_genes(&a),
            "annotate_preview" => Ok(vec![Content::Json(relabel_preview(&a, cx)?)]),
            "annotate_relabel" => {
                let round = relabel_apply(&a, cx)?;
                Ok(vec![Content::Json(json!({ "round": absolute(&round) }))])
            }
            "annotate_ask" => annotate_ask(&a, cx),
            "trajectory_run" => trajectory_run(&a, cx),
            "trajectory_read" => Ok(vec![Content::Json(trajectory_read(&a.req("run")?)?)]),
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
    fn str(&self, key: &str) -> Option<&str> {
        self.0.get(key).and_then(Value::as_str)
    }

    fn req(&self, key: &str) -> Result<String> {
        self.str(key)
            .map(str::to_string)
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

    /// `options`, refusing `--overwrite`: no tool replaces a manifest.
    fn options(&self) -> Result<Vec<String>> {
        let v: Vec<String> = self.parse("options")?.unwrap_or_default();
        anyhow::ensure!(
            !v.iter().any(|o| o == "--overwrite"),
            "--overwrite is not taken here: choose another `out`"
        );
        Ok(v)
    }
}

/// `lupin` itself, run with `argv`. A binary replaced while the server runs
/// (a `cargo install`) is gone from its path; the server's own still runs.
fn lupin(argv: &[String]) -> Result<Command> {
    let exe = std::env::current_exe().context("finding the lupin binary")?;
    let own = Path::new("/proc/self/exe");
    let mut cmd = Command::new(if !exe.exists() && own.exists() {
        own
    } else {
        &exe
    });
    cmd.args(argv);
    Ok(cmd)
}

fn absolute(p: &Path) -> String {
    std::path::absolute(p)
        .unwrap_or_else(|_| p.to_path_buf())
        .to_string_lossy()
        .into_owned()
}

/// `out`, or the next free `tag` name in today's results folder (beside the
/// manifest `from` names when the server keeps no records); and the
/// manifest a run from `from` to it writes, which must not exist yet.
fn output(a: &Arguments, cx: &Call, from: &str, tag: Tag) -> Result<(String, PathBuf)> {
    let source = manifest_file(from);
    let out = match a.str("out") {
        Some(o) => o.to_string(),
        None => match (cx.results_dir(), source.file_name()) {
            (Some(dir), Some(name)) => tag.next_free(&dir.join(name)),
            _ => tag.next_free(&source),
        },
    };
    let manifest = annotated_path(&source, &out);
    anyhow::ensure!(
        !manifest.exists(),
        "{} exists; not replacing it: choose another `out`, or omit it for the next free name",
        manifest.display()
    );
    Ok((out, manifest))
}

fn annotate_run(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let from = a.req("from")?;
    let (out, round) = output(a, cx, &from, Tag::Round)?;
    let mut argv: Vec<String> = ["annotate", "-f", &from, "-o", &out]
        .map(String::from)
        .into();
    if let Some(m) = a.str("markers") {
        argv.extend(["-m".into(), m.into()]);
    }
    if let Some(m) = a.str("method") {
        argv.extend(["--method".into(), m.into()]);
    }
    argv.extend(a.options()?);
    let summary = format!(
        "annotate {}{}{}",
        file_name(Path::new(&from)),
        a.str("markers").map_or(String::new(), |m| format!(
            " with {}",
            file_name(Path::new(m))
        )),
        a.str("method")
            .map_or(String::new(), |m| format!(" by {m}")),
    );
    cx.run_making(lupin(&argv)?, std::slice::from_ref(&round), &summary)?;
    let loaded = run::load(&round.to_string_lossy())?;
    let clusters: Map<String, Value> = review_json(&loaded, &[])?
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

fn annotate_genes(a: &Arguments) -> Result<Vec<Content>> {
    let loaded = run::load(&a.req("round")?)?;
    let clusters: Vec<ClusterId> = a.parse("clusters")?.unwrap_or_default();
    let more: Vec<String> = a.parse("types")?.unwrap_or_default();
    let top = a.parse::<usize>("top")?.unwrap_or(20);
    let at = |rel: &Option<String>| rel.as_deref().map(|r| resolve(&loaded.dir, r));
    let ann = &loaded.manifest.annotate;
    let profile = at(&ann.cluster_expression).context("the round records no expression profile")?;
    let expr = Expression::new(
        Mat::from_parquet_with_row_names(&profile, Some(0))
            .with_context(|| format!("reading {profile}"))?,
    );
    let specific = specific_genes(&expr);
    let panel = match at(&ann.markers) {
        Some(p) if Path::new(&p).is_file() => panel_sets(&p)?,
        _ => BTreeMap::new(),
    };
    let out: Map<String, Value> = review_json(&loaded, &clusters)?
        .into_iter()
        .map(|(id, v)| {
            let cols = expr.cols(&BTreeSet::from([id]));
            let fc = |g: &str| expr.row(g).and_then(|r| expr.log2fc(r, &cols));
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
                    let set = panel.get(&label_key(&t));
                    let mut genes: Vec<(&String, f32)> = set
                        .into_iter()
                        .flatten()
                        .filter_map(|g| Some((g, fc(g)?)))
                        .collect();
                    genes.sort_by(|x, y| y.1.total_cmp(&x.1));
                    let up = genes.iter().filter(|(_, f)| *f >= 1.0).count();
                    let listed: Vec<Value> =
                        genes.iter().take(top).map(|(g, f)| json!([g, f])).collect();
                    let n = set.map_or(0, BTreeSet::len);
                    let entry =
                        json!({ "markers": n, "measured": genes.len(), "up": up, "top": listed });
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

/// The decisions in `a`, one JSON line each. Each names the round (as
/// `round_ref`, relative to the decisions file), so it cannot land on
/// renumbered ids, and is the agent's unless it says otherwise.
fn decisions_text(a: &Arguments, round_ref: &str) -> Result<String> {
    let decisions: Vec<Map<String, Value>> = a.parse("decisions")?.unwrap_or_default();
    anyhow::ensure!(!decisions.is_empty(), "no decisions");
    let mut text = String::new();
    for mut d in decisions {
        d.entry("round").or_insert_with(|| json!(round_ref));
        d.entry("decided_by")
            .or_insert_with(|| json!("agent_proposed_user_accepted"));
        text += &format!("{}\n", Value::Object(d));
    }
    Ok(text)
}

/// The round `a` names, which must exist.
fn round_arg(a: &Arguments) -> Result<PathBuf> {
    let round = manifest_file(&a.req("round")?);
    anyhow::ensure!(round.is_file(), "{}: no such round", round.display());
    Ok(round)
}

/// `lupin relabel --preview` with the decisions in a temporary file; the
/// preview's JSON.
fn relabel_preview(a: &Arguments, cx: &Call) -> Result<Value> {
    let round = round_arg(a)?;
    let text = decisions_text(a, &absolute(&round))?;
    let tmp = tempfile::Builder::new().suffix(".jsonl").tempfile()?;
    std::fs::write(tmp.path(), text)?;
    let argv = [
        "relabel".to_string(),
        "-f".into(),
        round.to_string_lossy().into_owned(),
        "-d".into(),
        tmp.path().to_string_lossy().into_owned(),
        "--preview".into(),
    ];
    let out = cx.run_reading(lupin(&argv)?)?;
    serde_json::from_str(&out).context("reading the preview")
}

/// Apply the decisions as the next round of the round's chain, named here
/// (`-o`) rather than left to `--next`, so the step replays to the same
/// round. The decisions go to the next free `{round}.mcp_decisions.{k}.jsonl`
/// beside the round, and into the driver script.
fn relabel_apply(a: &Arguments, cx: &Call) -> Result<PathBuf> {
    let round = round_arg(a)?;
    let chain = Chain::of(&round);
    let (k, latest) = chain.latest();
    anyhow::ensure!(
        same_file(&round, &latest),
        "{} is not the latest round of its chain: review {} and decide on it",
        round.display(),
        latest.display()
    );
    let name = round.file_name().map(|n| n.to_string_lossy().into_owned());
    let text = decisions_text(a, &name.unwrap_or_default())?;
    let stem = derive_out_prefix(&round.to_string_lossy());
    let file = (1..)
        .map(|k| PathBuf::from(format!("{stem}.mcp_decisions.{k}.jsonl")))
        .find(|p| !p.exists())
        .expect("some k is free");
    cx.write_file(&file, &text)?;
    let prefix = chain.round_prefix(k + 1);
    let next = annotated_path(&round, &prefix);
    let argv = [
        "relabel".to_string(),
        "-f".into(),
        round.to_string_lossy().into_owned(),
        "-d".into(),
        file.to_string_lossy().into_owned(),
        "-o".into(),
        prefix,
    ];
    let summary = format!(
        "relabel {} with {} decision(s)",
        file_name(&round),
        text.lines().count()
    );
    cx.run_making(lupin(&argv)?, std::slice::from_ref(&next), &summary)?;
    Ok(next)
}

/// `p`'s file name, for summaries.
fn file_name(p: &Path) -> String {
    p.file_name().map_or_else(
        || p.display().to_string(),
        |n| n.to_string_lossy().into_owned(),
    )
}

fn annotate_ask(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let from = a.req("from")?;
    let (out, round) = output(a, cx, &from, Tag::Round)?;
    let mut argv: Vec<String> = ["ask", "-f", &from, "-o", &out].map(String::from).into();
    if let Some(c) = a.str("context") {
        argv.extend(["--context".into(), c.into()]);
    }
    if let Some(n) = a.parse::<usize>("top")? {
        argv.extend(["--top".into(), n.to_string()]);
    }
    let summary = format!(
        "unlabelled first round of {} for ask",
        file_name(Path::new(&from))
    );
    let prompt = cx.run_making(lupin(&argv)?, std::slice::from_ref(&round), &summary)?;
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
    let from = a.req("from")?;
    let (out, manifest) = output(a, cx, &from, Tag::Trajectory)?;
    let mut argv: Vec<String> = ["trajectory", "-f", &from, "-o", &out]
        .map(String::from)
        .into();
    for r in a.parse::<Vec<String>>("roots")?.unwrap_or_default() {
        argv.extend(["--root".into(), r]);
    }
    let statements: Vec<Statement> = a.parse("statements")?.unwrap_or_default();
    match (a.str("prior"), statements.is_empty()) {
        (Some(_), false) => anyhow::bail!("give `statements` or `prior`, not both"),
        (Some(p), true) => argv.extend(["--prior".into(), p.into()]),
        (None, false) => {
            let path = format!("{out}.precedence.tsv");
            let mut text = String::new();
            for s in &statements {
                for field in [&s.from, &s.to, &s.relation, &s.note] {
                    anyhow::ensure!(
                        !field.contains(['\t', '\n']),
                        "statement fields cannot hold tabs or newlines: {field:?}"
                    );
                }
                text += &format!("{}\t{}\t{}\t{}\n", s.from, s.to, s.relation, s.note);
            }
            // Left by a failed attempt at this run, it is replaced; any
            // other file there is not.
            if let Ok(old) = std::fs::read_to_string(&path) {
                anyhow::ensure!(
                    old == text,
                    "{path} exists; not replacing it: choose another `out`"
                );
            }
            cx.write_file(Path::new(&path), &text)?;
            argv.extend(["--prior".into(), path]);
        }
        (None, true) => {}
    }
    argv.extend(a.options()?);
    let roots: Vec<String> = a.parse("roots")?.unwrap_or_default();
    let summary = format!(
        "trajectory of {}{}{}",
        file_name(Path::new(&from)),
        if roots.is_empty() {
            String::new()
        } else {
            format!(" from {}", roots.join(", "))
        },
        if statements.is_empty() {
            String::new()
        } else {
            format!(", {} statement(s)", statements.len())
        },
    );
    cx.run_making(lupin(&argv)?, std::slice::from_ref(&manifest), &summary)?;
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
