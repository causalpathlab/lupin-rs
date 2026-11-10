//! `lupin mcp`: lupin's annotation and trajectory loops as Model Context
//! Protocol tools, served over stdio with `legume-mcp`, so an agent can run
//! them: run a pass, review it, preview and apply decisions, order the
//! types and read the order back.
//!
//! lupin states the facts (calls, q, ranks, edges, verdicts, all written by
//! code); the agent brings the judgment, as `relabel` decisions with a
//! rationale or as precedence statements with a note, both kept with the
//! round or run they made. A person follows the server with `lupin mcp
//! watch`, stops a call or leaves the agent a note there, and opens what a
//! call made on lupin's screen. Passes, relabels, `ask` and trajectories run as
//! child `lupin` processes, as from the TUI: stdout carries the protocol,
//! and their log lines come back as progress. No tool replaces an existing
//! manifest: a new output takes the next free name, in today's results
//! folder unless the agent names one. `legume_mcp::project` keeps the
//! project's records after Noble (PLoS Comput Biol 5(7): e1000424, 2009,
//! doi:10.1371/journal.pcbi.1000424): every step that writes, with the
//! inputs made for it, goes into the day's `runall.sh`, so the work can be
//! redone without the agent.

use crate::annotate::expression::{
    label_marker_set, profile_and_panel_on, specific_genes, Expression, Panel, UP_LOG2FC,
};
use crate::annotate::gene_rows::CHANNELS;
use crate::annotate::rounds::{ClusterId, DecidedBy};
use crate::child::{lupin_command, relabel_argv};
use crate::describe::trajectory::{self as order_text, OrderFacts};
use crate::figures::look::{LabelSize, Look, Overlay, PageKind};
use crate::manifest::family::{file_name, Tag};
use crate::manifest::rounds::{review_json, Chain};
use crate::manifest::run::{self, annotated_path, derive_out_prefix, manifest_file};
use crate::tui::headless::PanelState;
use anyhow::{Context, Result};
use clap::ValueEnum;
use data_beans::aux::feature_rows::SPLICED;
use legume_mcp::{Call, Content, Handler, ServerInfo, ToolSpec};
use serde_json::{json, Map, Value};
use std::collections::{BTreeMap, BTreeSet};
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
beyond it. A cluster's terms rank by effect, what sets it apart (killing, \
clotting, antigen presentation); their q tests that effect \
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
A first round's labels are coarse; give `level: fine` (or a fine type as a root) to order its fine labels. On fine labels a coarse root stands for the fine types mostly under it.
2. trajectory_read returns its lineages, each edge's connectivity and \
verdict (supported; weak: joined by the kNN graph, thinly; unsupported: not \
joined; candidate; inferred) and the order in words. \
Revise the statements where the data disagree and run again.

annotate_run and trajectory_run take `background`: on a large run set it, \
and the call returns a job id at once; the step runs on lupin's own screen in \
a session the user can watch with `lupin attach`, and goes on if this server \
ends.

Figures: annotate_heatmap draws a round's markers × clusters and returns its \
numbers too; annotate_figure draws a round's other charts (each cluster's top \
GO/GMT terms × clusters, the cell types × clusters, structure plots) and its \
cells on a layout, heatmaps with their numbers; trajectory_figure draws a \
trajectory run's map, order and connectivity, or with `maps` any of its layouts and diffusion component \
pairs (`list` first: it says how well each component follows pseudotime \
and parts the lineages, and suggests a pair). All three write into today's \
results folder like any step.

The figure tools take the look lupin's screen gives a map: `labels` (their \
size, or off), `focus` (a cluster lit, the rest grey), a units round's \
`overlays`, and the page's `width` and `dpi`; labels a person moved on the \
screen stay where they were put. Figures are drawn for where they go: \
`page: slide` (the default, for the day's report deck: the figure box's \
height, as wide as needed, a chart too tall turned right across it) or \
`page: paper` (a journal page's width, as tall as needed, a chart too wide \
turned left onto it); a result's `turned` names the charts turned.

Around the loops: annotate_family lists a run's rounds and trajectories; \
annotate_export writes a round's final annotation for use outside lupin. \
ontology browses the Cell Ontology over the panel; ontology_questions lists \
what lupin asks of it (types left alone that may join their siblings, types \
with no term), and ontology_answer answers one, with a note, for every later \
pass. A label for two types the evidence does not separate (A+B) gives them \
as `parts`. annotate_hide_genes keeps genes out of the specific genes; \
exports lists, moves and deletes the figures saved on lupin's screen.

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
        "description": "Decisions, one object each: {cluster or clusters, action: label|merge|keep|markers_add|markers_drop, label?, parts?, features?, evidence?: [..], alternatives?: [{label, why_not}], rationale (required)}. A label that stands for two or more types the evidence does not separate (A+B) gives them as `parts`: it is named so for the export, which places it under their common ancestor. decided_by defaults to agent (your own decision); give agent_proposed_user_accepted only when the user accepted it.",
        "items": { "type": "object" }
    });
    let obo = json!({
        "type": "string",
        "description": "The Cell Ontology OBO file; the one lupin's data files find when omitted"
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
            "Per cluster of a round, the genes behind its calls: the genes that set it apart (log2 fold change over the other clusters; those hidden by annotate_hide_genes left out and counted in `hidden`), and for its label, each method's top type and any `types` given, those types' panel markers: how many, how many are up (log2 fold change ≥ 1) and the `top` highest with their log2 fold change. A merged cluster's genes are its parts' together (their cell-weighted mean profile). A coarse label (a Cell Ontology group of the panel's types) is measured on its members' markers, listed in `group_of`. Markers match genes by case, ENSG..._SYMBOL either way and without a faba file's /count/spliced suffix; when under half the panel's markers are measured, a second item warns, with names that did not match. For why one label, gene or term rather than another, use annotate_why.",
            schema(
                json!({
                    "round": { "type": "string" },
                    "clusters": { "type": "array", "items": { "type": "integer" }, "description": "Only these clusters" },
                    "types": { "type": "array", "items": { "type": "string" }, "description": "More cell types whose markers to show" },
                    "top": { "type": "integer", "description": "Genes listed per cluster and per type (20)" },
                    "channel": { "type": "string", "enum": CHANNELS, "description": "Where a gene has several rows (a faba file's GENE/count/spliced and GENE/count/unspliced), the count channel its level is read from (spliced)" },
                    "show_hidden": { "type": "boolean", "description": "List the hidden genes (annotate_hide_genes) among the specific genes too" },
                }),
                &["round"],
            ),
        ),
        ToolSpec::new(
            "annotate_heatmap",
            "A round's marker genes × clusters as a heatmap (`lupin figure heatmap`): for each type its `top` markers, by their largest log2 fold change in any cluster, and their `value` in every cluster. Markers match genes as annotate_genes matches them. Writes the figure (.pdf, .svg or .png, by `out`) and its numbers beside it (.json) into today's results folder, recorded as any step; returns the numbers (clusters left to right, rows top to bottom, `unmeasured`: types none of whose markers is measured) and the image.",
            schema(
                json!({
                    "round": { "type": "string" },
                    "types": { "type": "array", "items": { "type": "string" }, "description": "Cell types whose markers to show (every panel type when omitted)" },
                    "top": { "type": "integer", "description": "Markers per type (5)" },
                    "value": { "type": "string", "enum": ["log2fc", "mean", "zscore", "frac"], "description": "log2 fold change over the other clusters, mean level, that level z-scored across clusters (default), or the share of the cluster's cells the gene is detected in (read from the run's raw counts)" },
                    "source": { "type": "string", "enum": ["observed", "expected"], "description": "The round's observed profile (default), or the topic model's expected level (a topic run's)" },
                    "order": { "type": "string", "enum": ["label", "id"], "description": "Clusters by label, following the types' rows (default), or by lupin's cluster id" },
                    "channel": { "type": "string", "enum": CHANNELS, "description": "Where a gene has several rows, the count channel its level is read from (spliced)" },
                    "out": { "type": "string", "description": "The figure, ending .pdf, .svg or .png; {round}.heatmap.pdf in today's results folder when omitted" },
                }),
                &["round"],
            ),
        ),
        ToolSpec::new(
            "annotate_figure",
            "A round's charts and cell maps (`lupin figure round`), as the screen draws them. Charts: `terms`, each cluster's top GO/GMT terms × clusters (by effect, z by row, blank where q is not below `q`; a round annotated with go, gaf or gmt); `types`, the cell types × clusters; `markers`, the panel's markers × clusters; `structure`, each cell's cell types stacked by cluster; `topics`, a topic run's topics so. Maps: the round's cells on a layout, coloured by type, fine type (a coarse first round's fine labels), cluster, call or margin. Each chart is sized to what it shows and to its `page`: on a slide as wide as it needs at the slide's height, and one taller than the slide (many rows) turned right across it; on paper as tall as it needs at the page's width, and one wider than the page (many clusters) turned left onto it; PNG and PDF alike, listed in `turned`. Each chart as {out}.{chart}.png and .pdf, a heatmap's numbers as .json too (rows, cols, z, the unscaled values and their q); each map as {out}.{map}.{colour}.png and .pdf. Written into today's results folder, recorded as any step; returns the files, each heatmap's numbers, the images, and what the round offers (each chart, or why it has none; its layouts and colourings). `list` returns that alone, drawing nothing. For the markers heatmap with its own options (types, value, source), use annotate_heatmap.",
            schema(
                json!({
                    "round": { "type": "string" },
                    "charts": { "type": "array", "items": { "type": "string", "enum": ["terms", "types", "markers", "structure", "topics"] }, "description": "The charts to draw; every one the round records when neither charts nor maps is given" },
                    "maps": { "type": "array", "items": { "type": "string" }, "description": "Layouts to draw the cells on: phate, umap, …, best or all" },
                    "colour": { "type": "string", "enum": ["type", "fine", "cluster", "call", "margin"], "description": "What colours the map's cells (type); fine: the fine labels a coarse first round keeps, by cell" },
                    "labels": { "type": "string", "enum": ["small", "medium", "large", "largest", "off"], "description": "The type labels on a map (medium)" },
                    "focus": { "type": "integer", "description": "Light this cluster's cells, the others grey" },
                    "page": { "type": "string", "enum": ["slide", "paper"], "description": "What the figures are set on: slide (default: the day's report deck, a figure box 3.5 in tall and as wide as needed; a chart too tall is turned onto it) or paper (a journal page 7 in wide, as tall as needed; a chart too wide is turned onto it)" },
                    "width": { "type": "number", "description": "Paper: the page's width in inches (7); slide: the widest a figure gets (12.2)" },
                    "dpi": { "type": "integer", "description": "Dots per inch of the cells (200); a PDF's text stays text" },
                    "overlays": { "type": "array", "items": { "type": "string", "enum": ["go", "links", "rings", "none"] }, "description": "A units round's map: the lit cluster's GO terms as stars, the strongest links between types, each type's fitted radius as a ring (go and links; none for none)" },
                    "top": { "type": "integer", "description": "The terms heatmap: each cluster's top terms by effect (6)" },
                    "q": { "type": "number", "description": "The terms heatmap: only terms with q below this (0.05)" },
                    "out": { "type": "string", "description": "Output prefix; {round}.figure in today's results folder when omitted" },
                    "list": { "type": "boolean", "description": "Only list the charts and maps the round offers" },
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
                    "level": { "type": "string", "enum": ["coarse", "fine"], "description": "Which of the round's labels to order: its own (coarse: a first round's are Cell Ontology groups) or the fine labels a coarse first round keeps. Left out: fine when a root is a fine label and not a coarse one; the result's `level` says which. On fine labels, a root that is one of the round's own labels stands for each fine type most of whose cells carry it" },
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
            "Read a trajectory run: the level of the round's labels it ordered (coarse or fine), its lineages (root to leaf), its roots, branch points and ends, each type pair's connectivity, place in the prior, verdict and order agreement, and the order in words.",
            schema(
                json!({ "run": { "type": "string", "description": "Trajectory manifest or its prefix" } }),
                &["run"],
            ),
        ),
        ToolSpec::new(
            "trajectory_figure",
            "A trajectory run's figures (`lupin figure trajectory`). Without `maps`: pseudotime on its best map (PHATE or another layout, else the diffusion map), as PNG and PDF; the types in pseudotime order and the connectivity between them (the lineage graph), as SVG and PDF. With `maps`: only those maps, each as {out}.{map}.{colour}.png and .pdf. A map has the prior's edges as arrows between the types; a categorical colouring draws its largest group first, so small ones stay visible. Written into today's results folder, recorded as any step; returns the files, each map's image, and `maps`: what the run offers (its layouts; per diffusion component its sd, share of cells, whether pseudotime kept it, |Spearman ρ| with pseudotime and η² between the lineages; the pair the diffusion map opens on; the colourings). `list` returns `maps` alone, drawing nothing.",
            schema(
                json!({
                    "run": { "type": "string", "description": "Trajectory manifest or its prefix" },
                    "out": { "type": "string", "description": "Output prefix; {run}.figure in today's results folder when omitted" },
                    "maps": { "type": "array", "items": { "type": "string" }, "description": "The maps to draw: a layout by name (phate, umap, …), best, all (every layout and the diffusion map), suggested or diffusion (the diffusion map on its suggested pair: of the kept components, the one following pseudotime best on x, the one parting the lineages best on y), or two diffusion components as dc1,dc3 (DC0 is trivial and refused)" },
                    "colour": { "type": "string", "enum": ["pseudotime", "type", "fine", "lineage", "component"], "description": "What colours the cells (pseudotime); type is the run's level of labels; fine, on a run of a coarse first round's labels, the fine labels it keeps, by cell" },
                    "labels": { "type": "string", "enum": ["small", "medium", "large", "largest", "off"], "description": "The type labels on a map (medium)" },
                    "focus": { "type": "integer", "description": "Light this cluster's cells, the others grey" },
                    "page": { "type": "string", "enum": ["slide", "paper"], "description": "What the figures are set on: slide (default: the day's report deck, a figure box 3.5 in tall and as wide as needed; a chart too tall is turned onto it) or paper (a journal page 7 in wide, as tall as needed; a chart too wide is turned onto it)" },
                    "width": { "type": "number", "description": "Paper: the page's width in inches (7); slide: the widest a figure gets (12.2)" },
                    "dpi": { "type": "integer", "description": "Dots per inch of the cells (200); a PDF's text stays text" },
                    "list": { "type": "boolean", "description": "Only list the maps and colourings the run offers, with each diffusion component's numbers" },
                }),
                &["run"],
            ),
        ),
        ToolSpec::new(
            "annotate_family",
            "The family of a run, round or trajectory (`g` on lupin's screen): the run itself, its annotated rounds and its trajectories, each with its path, kind (run, round, trajectory), whether it carries labels, what it holds and when it last changed (seconds since the epoch). Any member is what the other tools take as `round`, `from` or `run`.",
            schema(
                json!({ "from": { "type": "string", "description": "A run, round or trajectory manifest or its prefix" } }),
                &["from"],
            ),
        ),
        ToolSpec::new(
            "annotate_export",
            "A round's final annotation for use outside lupin (`lupin export`, `e` on lupin's screen): {out}.cell_annotation.parquet, per cell its cluster, label, the label's Cell Ontology term and lineage, the cluster's evidence share for the label and its top type; and {out}.celltype_markers.tsv, per label its markers (from the panel, or added since) and the genes that set its clusters apart, each with its log2 fold change. Written into today's results folder, recorded as any step; returns the files and each label's cells and term.",
            schema(
                json!({
                    "round": { "type": "string" },
                    "out": { "type": "string", "description": "Output prefix; {round} in today's results folder when omitted" },
                    "obo": obo.clone(),
                }),
                &["round"],
            ),
        ),
        ToolSpec::new(
            "ontology",
            "The Cell Ontology as lupin's screen browses it, over a round's marker panel. With `search`: the terms whose names, synonyms or abbreviations hold it. With `term` (an id): its name, parents, children, lineage from the top, what it develops from, and the panel's types at or under it. With neither: each panel type's term and lineage, and the types with none (`unmatched`; ontology_questions asks about them).",
            schema(
                json!({
                    "round": { "type": "string" },
                    "search": { "type": "string" },
                    "term": { "type": "string", "description": "A term's id (CL:…)" },
                    "limit": { "type": "integer", "description": "Terms a search returns (10)" },
                    "obo": obo.clone(),
                }),
                &["round"],
            ),
        ),
        ToolSpec::new(
            "ontology_questions",
            "The questions lupin's screen asks of a round's marker panel until they are answered: `join`, types left alone (their siblings went under a nearer class) that may join them under `class`; `terms`, types with no Cell Ontology term. With `suggest`, each such type's nearest terms by text embedding, with their scores (the first time fetches the text model and encodes the ontology: a minute or two; none offline without the model). Answer each with ontology_answer.",
            schema(
                json!({
                    "round": { "type": "string" },
                    "suggest": { "type": "boolean", "description": "Rank terms for the types with none" },
                    "obo": obo.clone(),
                }),
                &["round"],
            ),
        ),
        ToolSpec::new(
            "ontology_answer",
            "Answer one of ontology_questions for a panel type, kept with your note where every later pass reads it and records it (the user's answers, else the project's): `join` or `alone` for a type left alone; a term's id (CL:…) for a type with none, kept in the alias table; `none` when no term fits. The next pass groups and places the type so; nothing already written changes. Returns the file written and, for a term, the type's reading now.",
            schema(
                json!({
                    "round": { "type": "string" },
                    "type": { "type": "string", "description": "The panel type asked about" },
                    "answer": { "type": "string", "description": "join, alone, none, or a term's id (CL:…)" },
                    "note": { "type": "string", "description": "Why; kept with the answer" },
                    "obo": obo.clone(),
                }),
                &["round", "type", "answer", "note"],
            ),
        ),
        ToolSpec::new(
            "annotate_hide_genes",
            "Hide genes from the specific genes annotate_genes and lupin's screen list, from now on (`h` and `X` on the screen): exact names or patterns (PFX-* for every gene starting PFX-). Kept in the project's hidden_genes.txt (else the user's); the markers are never hidden. Returns the file written.",
            schema(
                json!({
                    "round": { "type": "string" },
                    "patterns": { "type": "array", "items": { "type": "string" } },
                }),
                &["round", "patterns"],
            ),
        ),
        ToolSpec::new(
            "exports",
            "The exports log lupin's screen keeps of the figures saved (`f` on the screen, in .lupin-view/ under the server's working directory), each figure set with its files, size and hash checked against the disk. `list`: each entry, newest first, with its status (ok, changed since export, missing) and the manifest and figure it was drawn from. `add`: log a figure file (its .pdf, with its .png or .svg beside it). `unlist`: take one off the log. `move`: move its files to `to` (a folder, or a new name), never over a file. `delete`: delete its files for good, only while they are as logged and with `confirm` (ask the user first, as the screen asks twice).",
            schema(
                json!({
                    "action": { "type": "string", "enum": ["list", "add", "unlist", "move", "delete"] },
                    "path": { "type": "string", "description": "The figure's file, as listed (not for list)" },
                    "to": { "type": "string", "description": "For move: a folder, or the new name" },
                    "what": { "type": "string", "description": "For add: what the figure shows" },
                    "confirm": { "type": "boolean", "description": "For delete: true, once the user agrees; the files cannot be brought back" },
                }),
                &["action"],
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
/// one picked among several). A person sees what the agent does, can stop
/// a call or job and leave the agent a note, and opens on lupin's own
/// screen (`lupin view`) the round or trajectory a call made or read.
pub fn watch(dir: Option<PathBuf>) -> Result<()> {
    // The log stays off the screen while it is up.
    let _quiet = legume_tui::logs::quiet_stderr();
    legume_mcp::watch::view::run(&crate::tui::TOOL, dir, true, |call| {
        let manifest = call_manifest(&call.args, &call.text)?;
        let mut view = legume_mcp::own_command(["view"]).ok()?;
        view.arg(manifest);
        Some(view)
    })
}

/// The round or trajectory run a call made or read: the one its result
/// names (`round`, `run`), else the one its arguments name; only once its
/// manifest is there (a background job's is not until it ends). The feed
/// cuts a result at 4 KB: past that, the arguments name it.
fn call_manifest(args: &Value, text: &str) -> Option<PathBuf> {
    // The first value: a note the person left may follow it.
    let result = serde_json::Deserializer::from_str(text)
        .into_iter::<Value>()
        .next()
        .and_then(Result::ok)
        .unwrap_or(Value::Null);
    let named = |v: &Value| {
        ["round", "run"]
            .into_iter()
            .find_map(|k| v.get(k)?.as_str().map(manifest_file))
    };
    let source = args.get("from").and_then(Value::as_str).map(manifest_file);
    [named(&result), named(args), source]
        .into_iter()
        .flatten()
        .find(|m| m.is_file())
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
            "annotate_heatmap" => annotate_heatmap(&a, cx),
            "annotate_figure" => annotate_figure(&a, cx),
            "trajectory_figure" => trajectory_figure(&a, cx),
            "trajectory_run" => trajectory_run(&a, cx),
            "trajectory_read" => Ok(vec![Content::Json(trajectory_read(&a.req_path("run")?)?)]),
            "annotate_family" => annotate_family(&a),
            "annotate_export" => annotate_export(&a, cx),
            "ontology" => ontology(&a),
            "ontology_questions" => ontology_questions(&a, cx),
            "ontology_answer" => ontology_answer(&a),
            "annotate_hide_genes" => {
                let patterns: Vec<String> = a.parse("patterns")?.unwrap_or_default();
                let file =
                    crate::tui::headless::hide_genes(&run::parent_dir(&round_arg(&a)?), &patterns)?;
                Ok(vec![Content::Json(json!({ "written": absolute(&file) }))])
            }
            "exports" => exports(Path::new(crate::tui::gallery::DIR), &a),
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

    fn req_str(&self, key: &str) -> Result<&str> {
        self.str(key)
            .with_context(|| format!("`{key}` is required"))
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
    none_exist(std::slice::from_ref(&manifest))?;
    cx.ensure_free(std::slice::from_ref(&manifest))?;
    legume_numeric::matrix::common_io::mkdir_parent(&out)?;
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
            label: Some(session_label(out)),
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

/// A background job's session label: the start of its output's name, to
/// know it by, and a short hash of the whole prefix, to tell it apart. The
/// label is in the session's socket name, whose whole path a Unix socket
/// holds to about a hundred bytes, so a long name cannot be all of it.
fn session_label(out: &str) -> String {
    use sha2::{Digest, Sha256};
    const NAME: usize = 6;
    let name: String = (Path::new(out).file_name())
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default()
        .trim_start_matches('.')
        .chars()
        .map(|c| match c {
            c if c.is_ascii_alphanumeric() || matches!(c, '-' | '_' | '.') => c,
            _ => '-',
        })
        .take(NAME)
        .collect();
    let hash = crate::hex(&Sha256::digest(out.as_bytes())[..4]);
    if name.is_empty() {
        hash
    } else {
        format!("{name}-{hash}")
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
    let channel = channel_arg(a)?;
    let cells = crate::manifest::rounds::read_cells(&loaded)?;
    let (expr, panel) = profile_and_panel_on(&loaded, &cells, channel)?;
    let expr = expr.context("the round records no expression profile")?;
    let warnings = unmatched_warning(&expr, &panel);
    let specific = specific_genes(&expr);
    // The genes the screen hides from them, unless asked for.
    let hidden = if a.parse::<bool>("show_hidden")?.unwrap_or(false) {
        None
    } else {
        let search = crate::manifest::data_files::SearchPath::new(Some(&loaded.dir));
        Some(crate::tui::genes::GeneFilter::load(&search)?)
    };
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
            let (mut listed, mut hid) = (Vec::new(), 0);
            for (g, f) in specific.get(&id).into_iter().flatten() {
                if hidden.as_ref().is_some_and(|h| h.hides(g)) {
                    hid += 1;
                } else if listed.len() < top {
                    listed.push(json!([g, f]));
                }
            }
            let mut entry = json!({ "specific": listed, "markers": markers });
            if hid > 0 {
                entry["hidden"] = json!(hid);
            }
            (id.to_string(), entry)
        })
        .collect();
    let mut content = vec![Content::Json(Value::Object(out))];
    if let Some(w) = warnings {
        content.push(Content::Json(json!({ "warnings": [w] })));
    }
    Ok(content)
}

/// The `channel` `a` names, spliced when none: one of [`CHANNELS`].
fn channel_arg(a: &Arguments) -> Result<&str> {
    let channel = a.str("channel").unwrap_or(SPLICED);
    anyhow::ensure!(
        CHANNELS.contains(&channel),
        "`channel` is one of {}, not {channel:?}",
        CHANNELS.join(", ")
    );
    Ok(channel)
}

/// Below this share of the panel's markers measured, the genes shown say
/// little of the panel, and a warning says so.
const MATCHED_SHARE: f32 = 0.5;

/// What to say when few of `panel`'s markers are rows of `expr`: how many,
/// with some that are not and some row names, so a naming mismatch shows.
fn unmatched_warning(expr: &Expression, panel: &Panel) -> Option<String> {
    let genes: BTreeSet<&String> = panel.values().flatten().collect();
    let (matched, unmatched): (Vec<&String>, Vec<&String>) =
        genes.iter().partition(|g| expr.row(g).is_some());
    let share = matched.len() as f32 / genes.len().max(1) as f32;
    if genes.is_empty() || share >= MATCHED_SHARE {
        return None;
    }
    let some = |names: Vec<&str>| names.into_iter().take(5).collect::<Vec<_>>().join(", ");
    Some(format!(
        "{} of the panel's {} markers are measured in the round's expression profile, \
         so its markers' counts here say little: not matched, e.g. {}; the profile's rows \
         are named e.g. {}",
        matched.len(),
        genes.len(),
        some(unmatched.iter().map(|g| g.as_str()).collect()),
        some(expr.table.rows.iter().map(|r| &**r).collect()),
    ))
}

/// A free name `base` in today's results folder, else beside `beside`
/// when the server keeps no records, with the outputs `plan` says it
/// writes.
fn free_name(
    cx: &Call,
    beside: &Path,
    base: &str,
    plan: impl Fn(&str) -> Result<Vec<PathBuf>>,
) -> Result<(String, Vec<PathBuf>)> {
    if cx.results_dir().is_some() {
        cx.free_prefix(base, plan, |_| false)
    } else {
        legume_mcp::free_prefix(&run::parent_dir(beside), base, plan, |_| false)
    }
}

/// The output prefix `a` names, none of whose files (`plan`) may exist, else
/// the next free `base` beside `beside` (in today's results folder when
/// records are kept); and its files.
fn out_prefix(
    a: &Arguments,
    cx: &Call,
    beside: &Path,
    base: &str,
    plan: impl Fn(&str) -> Result<Vec<PathBuf>>,
) -> Result<(String, Vec<PathBuf>)> {
    match a.path("out") {
        Some(o) => {
            let files = plan(&o)?;
            none_exist(&files)?;
            Ok((o, files))
        }
        None => free_name(cx, beside, base, plan),
    }
}

/// `outputs` at a name the agent gave: none may exist yet.
fn none_exist(outputs: &[PathBuf]) -> Result<()> {
    if let Some(p) = outputs.iter().find(|p| p.exists()) {
        anyhow::bail!(
            "{} exists; not replacing it: choose another `out`, or omit it for the next free name",
            p.display()
        );
    }
    Ok(())
}

/// The name of the manifest `file` without its `.{kind}.json`.
fn stem(file: &Path) -> String {
    file_name(Path::new(&derive_out_prefix(&file.to_string_lossy())))
}

fn annotate_heatmap(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    use crate::figures::heatmap::{self, figure_base, numbers_path};
    let round = round_arg(a)?;
    channel_arg(a)?;
    let plan = |figure: &str| -> Result<Vec<PathBuf>> {
        figure_base(figure)?;
        Ok(vec![PathBuf::from(figure), numbers_path(figure)])
    };
    let (out, outputs) = match a.path("out") {
        Some(o) => {
            let outputs = plan(&o)?;
            none_exist(&outputs)?;
            (o, outputs)
        }
        None => {
            let base = format!("{}.heatmap", stem(&round));
            let (b, outputs) = free_name(cx, &round, &base, |b| plan(&format!("{b}.pdf")))?;
            (format!("{b}.pdf"), outputs)
        }
    };
    let mut argv: Vec<String> = vec![
        "figure".into(),
        "heatmap".into(),
        flag("from", &round.to_string_lossy()),
        flag("out", &out),
    ];
    let types: Vec<String> = a.parse("types")?.unwrap_or_default();
    argv.extend(types.iter().map(|t| flag("type", t)));
    if let Some(n) = a.parse::<usize>("top")? {
        argv.push(flag("top", &n.to_string()));
    }
    for key in ["value", "source", "order", "channel"] {
        if let Some(v) = a.str(key) {
            argv.push(flag(key, v));
        }
    }
    let summary = format!(
        "marker heatmap ({}) of {}",
        a.str("value").unwrap_or("zscore"),
        file_name(&round)
    );
    cx.run_making(lupin(&argv)?, &outputs, &summary)?;
    let text = std::fs::read_to_string(&outputs[1])
        .with_context(|| format!("reading {}", outputs[1].display()))?;
    let data: heatmap::HeatData = serde_json::from_str(&text)?;
    let mut v = serde_json::to_value(&data)?;
    v["figure"] = json!(absolute(Path::new(&out)));
    v["numbers"] = json!(absolute(&outputs[1]));
    Ok(vec![Content::Json(v), Content::Png(heatmap::png(&data)?)])
}

fn annotate_figure(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    use crate::figures::round::{round_figures, ChartKind};
    use crate::figures::units::{STARS, STAR_Q};
    let from = a.req_path("round")?;
    let mut figures = round_figures(&from)?;
    let offered = json!({
        "charts": figures.charts_offered(),
        "maps": figures.maps_offered(),
    });
    if a.parse::<bool>("list")?.unwrap_or(false) {
        return Ok(vec![Content::Json(json!({ "offers": offered }))]);
    }
    let charts: Vec<String> = a.parse("charts")?.unwrap_or_default();
    let kinds: Vec<ChartKind> = (charts.iter())
        .map(|c| ChartKind::from_str(c, true).map_err(|e| anyhow::anyhow!("chart {c:?}: {e}")))
        .collect::<Result<_>>()?;
    let maps: Vec<String> = a.parse("maps")?.unwrap_or_default();
    let colour = a.str("colour");
    let (top, q) = (a.parse::<usize>("top")?, a.parse::<f32>("q")?);
    figures.set_term_cut(top.unwrap_or(STARS), q.unwrap_or(STAR_Q))?;
    // A chart the round lacks, a map or a colouring is refused here, saying
    // why, before anything runs.
    let look = look_arg(a)?;
    let drawing = figures.round_drawing(&kinds, &maps, colour, &look)?;
    let mut flags: Vec<String> = (drawing.charts.iter())
        .map(|k| flag("chart", &k.name()))
        .collect();
    flags.extend(
        figures
            .map_specs(&drawing.maps)
            .iter()
            .map(|m| flag("map", m)),
    );
    flags.extend(colour.map(|c| flag("colour", c)));
    flags.extend(look_flags(&look));
    flags.extend(top.map(|t| flag("top", &t.to_string())));
    flags.extend(q.map(|q| flag("q", &q.to_string())));
    let summary = format!("charts and maps of {}", file_name(Path::new(&from)));
    let (out, files) = draw_figures(a, cx, "round", &from, flags, &summary, |b| {
        figures.round_files(&drawing, b)
    })?;
    let mut numbers = Map::new();
    for (chart, path) in crate::figures::CellFigures::round_numbers(&drawing, &out) {
        let text = std::fs::read_to_string(&path)
            .with_context(|| format!("reading {}", path.display()))?;
        numbers.insert(chart, serde_json::from_str(&text)?);
    }
    let mut facts = json!({ "numbers": numbers, "offers": offered });
    let turned = figures.round_turned(&drawing)?;
    if !turned.is_empty() {
        facts["turned"] = json!(turned);
    }
    figure_content(facts, &files)
}

fn trajectory_figure(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let from = a.req_path("run")?;
    let figures = crate::figures::trajectory_figures(&from)?;
    let offered = serde_json::to_value(figures.maps_offered())?;
    if a.parse::<bool>("list")?.unwrap_or(false) {
        return Ok(vec![Content::Json(json!({ "maps": offered }))]);
    }
    let maps: Vec<String> = a.parse("maps")?.unwrap_or_default();
    let colour = a.str("colour");
    // A map or colouring the run lacks is refused here, naming what it has.
    let look = look_arg(a)?;
    let drawing = figures.drawing(&maps, colour, &look)?;
    // The maps as the child reads them back, `suggested` as its pair.
    let mut flags: Vec<String> = (figures.map_specs(&drawing).iter())
        .map(|m| flag("map", m))
        .collect();
    flags.extend(colour.map(|c| flag("colour", c)));
    flags.extend(look_flags(&look));
    let mut summary = format!("trajectory figures of {}", file_name(Path::new(&from)));
    if !maps.is_empty() {
        summary += &format!(" on {}", maps.join(", "));
    }
    let (_, files) = draw_figures(a, cx, "trajectory", &from, flags, &summary, |b| {
        figures.drawing_files(&drawing, b)
    })?;
    figure_content(json!({ "maps": offered }), &files)
}

/// How a figure tool's maps look and the page they are written on, as `a`
/// asks: what `lupin figure` takes as its look.
fn look_arg(a: &Arguments) -> Result<Look> {
    // A figure an agent draws goes to the day's report, a deck of slides,
    // unless asked for paper.
    let mut look = Look {
        page: PageKind::Slide,
        ..Look::default()
    };
    if let Some(p) = a.str("page") {
        look.page = PageKind::from_str(p, true).map_err(|e| anyhow::anyhow!("`page`: {e}"))?;
    }
    if let Some(l) = a.str("labels") {
        look.labels =
            Some(LabelSize::from_str(l, true).map_err(|e| anyhow::anyhow!("`labels`: {e}"))?);
    }
    look.focus = a.parse("focus")?;
    for o in a.parse::<Vec<String>>("overlays")?.unwrap_or_default() {
        look.overlays
            .push(Overlay::from_str(&o, true).map_err(|e| anyhow::anyhow!("`overlays`: {e}"))?);
    }
    if let Some(w) = a.parse::<f32>("width")? {
        look.width = Some(w);
    }
    if let Some(d) = a.parse::<u32>("dpi")? {
        look.dpi = d;
    }
    Ok(look)
}

/// `look` as `lupin figure`'s flags, those it changes from the default.
fn look_flags(look: &Look) -> Vec<String> {
    let name = |v: Option<clap::builder::PossibleValue>| v.map(|v| v.get_name().to_string());
    let base = Look::default();
    let mut flags = Vec::new();
    if let Some(l) = look.labels.and_then(|l| name(l.to_possible_value())) {
        flags.push(flag("labels", &l));
    }
    if let Some(k) = look.focus {
        flags.push(flag("focus", &k.to_string()));
    }
    let overlays: Vec<String> = (look.overlays.iter())
        .filter_map(|o| name(o.to_possible_value()))
        .collect();
    if !overlays.is_empty() {
        flags.push(flag("overlay", &overlays.join(",")));
    }
    if let Some(w) = look.width {
        flags.push(flag("width", &w.to_string()));
    }
    if look.page != base.page {
        flags.extend(name(look.page.to_possible_value()).map(|p| flag("page", &p)));
    }
    if look.dpi != base.dpi {
        flags.push(flag("dpi", &look.dpi.to_string()));
    }
    flags
}

/// `lupin figure {sub}` of `from` with `flags`, into the `out` `a` names
/// (none of whose files may exist) or the next free `{stem}.figure`, as
/// `plan` names the files it writes; recorded as `summary`. The output
/// prefix and the files.
fn draw_figures(
    a: &Arguments,
    cx: &Call,
    sub: &str,
    from: &str,
    flags: Vec<String>,
    summary: &str,
    plan: impl Fn(&str) -> Vec<PathBuf>,
) -> Result<(String, Vec<PathBuf>)> {
    let file = manifest_file(from);
    let (out, files) = out_prefix(a, cx, &file, &format!("{}.figure", stem(&file)), |b| {
        Ok(plan(b))
    })?;
    let mut argv = vec![
        "figure".into(),
        sub.into(),
        flag("from", from),
        flag("out", &out),
    ];
    argv.extend(flags);
    cx.run_making(lupin(&argv)?, &files, summary)?;
    Ok((out, files))
}

/// A figure tool's result: `facts` with the files written, then each PNG
/// among them as an image.
fn figure_content(mut facts: Value, files: &[PathBuf]) -> Result<Vec<Content>> {
    facts["files"] = json!(files.iter().map(|f| absolute(f)).collect::<Vec<_>>());
    let mut content = vec![Content::Json(facts)];
    for png in files
        .iter()
        .filter(|f| f.extension().is_some_and(|e| e == "png"))
    {
        let bytes = std::fs::read(png).with_context(|| format!("reading {}", png.display()))?;
        content.push(Content::Png(bytes));
    }
    Ok(content)
}

/// A round's panel as lupin's screen holds it, on the ontology `obo` names
/// when it names one.
fn panel_state(a: &Arguments) -> Result<PanelState> {
    let round = round_arg(a)?;
    PanelState::load(&round.to_string_lossy(), a.path("obo").as_deref())
}

fn annotate_family(a: &Arguments) -> Result<Vec<Content>> {
    let from = manifest_file(&a.req_path("from")?);
    anyhow::ensure!(from.is_file(), "{}: no such manifest", from.display());
    let members: Vec<Value> = (crate::manifest::family::family(&from).iter())
        .map(|m| {
            let mut v = json!({
                "name": m.name(),
                "path": absolute(&m.pick.path),
                "kind": m.pick.kind.name(),
                "labelled": m.pick.labelled,
                "holds": m.holds,
                "modified": m.modified,
            });
            if let Some(b) = &m.borrowed {
                v["borrowed"] = json!(b);
            }
            v
        })
        .collect();
    Ok(vec![Content::Json(json!({ "members": members }))])
}

fn annotate_export(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    use legume_numeric::matrix::parquet::read_parquet_string_columns_by_name;
    let round = round_arg(a)?;
    let (out, files) = out_prefix(a, cx, &round, &stem(&round), |b| {
        let (cells, markers) = crate::tui::export::files(b);
        Ok(vec![cells, markers])
    })?;
    let mut argv = vec![
        "export".into(),
        flag("from", &round.to_string_lossy()),
        flag("out", &out),
    ];
    if let Some(obo) = a.path("obo") {
        argv.push(flag("obo", &obo));
    }
    cx.run_making(
        lupin(&argv)?,
        &files,
        &format!("export of {}", file_name(&round)),
    )?;
    let cols = read_parquet_string_columns_by_name(
        &files[0].to_string_lossy(),
        &["label", "cl_id", "lineage"],
    )?;
    let [label, cl_id, lineage] = &cols[..] else {
        anyhow::bail!("{}: not the export's columns", files[0].display());
    };
    // Each label's cells, counted on the columns as read.
    let mut counted: BTreeMap<&str, (u64, &str, &str)> = BTreeMap::new();
    for ((l, id), line) in label.iter().zip(cl_id).zip(lineage) {
        counted.entry(l).or_insert((0, id, line)).0 += 1;
    }
    let labels: Map<String, Value> = (counted.into_iter())
        .map(|(l, (n, id, line))| {
            (
                l.into(),
                json!({ "cells": n, "cl_id": id, "lineage": line }),
            )
        })
        .collect();
    Ok(vec![Content::Json(json!({
        "cells": absolute(&files[0]),
        "markers": absolute(&files[1]),
        "labels": labels,
    }))])
}

fn ontology(a: &Arguments) -> Result<Vec<Content>> {
    let state = panel_state(a)?;
    let v = match (a.str("search"), a.str("term")) {
        (Some(_), Some(_)) => anyhow::bail!("give `search` or `term`, not both"),
        (Some(q), None) => state.search(q, a.parse::<usize>("limit")?.unwrap_or(10))?,
        (None, Some(id)) => state.term_facts(id)?,
        (None, None) => state.on_ontology()?,
    };
    Ok(vec![Content::Json(v)])
}

fn ontology_questions(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let state = panel_state(a)?;
    let mut v = state.questions()?;
    let types: Vec<String> = (v["terms"].as_array().into_iter().flatten())
        .filter_map(|t| t["type"].as_str().map(str::to_string))
        .collect();
    if a.parse::<bool>("suggest")?.unwrap_or(false) && !types.is_empty() {
        let obo = a.path("obo");
        let argv =
            crate::tui::runner::suggest_argv(&state.loaded.file, obo.as_deref(), &types, None);
        let ranked = crate::tui::suggest::Suggestions::parse(&cx.run_reading(lupin(&argv)?)?)?;
        for t in v["terms"].as_array_mut().into_iter().flatten() {
            let found = t["type"].as_str().and_then(|ty| ranked.of(ty));
            t["suggestions"] = json!(found.unwrap_or_default());
        }
    }
    Ok(vec![Content::Json(v)])
}

fn ontology_answer(a: &Arguments) -> Result<Vec<Content>> {
    let state = panel_state(a)?;
    let (ty, answer) = (a.req_str("type")?, a.req_str("answer")?);
    let file = state.answer(ty, answer, a.str("note").unwrap_or_default())?;
    let mut v = json!({ "written": absolute(&file) });
    // Read as the next pass reads it: a later table may map the type too.
    let now = panel_state(a)?.on_ontology()?;
    if let Some(on) = now["types"].get(ty) {
        v["now"] = on.clone();
    }
    Ok(vec![Content::Json(v)])
}

/// The `exports` tool on the log `g`.
fn exports(g: &Path, a: &Arguments) -> Result<Vec<Content>> {
    use crate::tui::gallery;
    let mut log = gallery::open(g);
    // As the log keeps it: resolved, once the file is there.
    let path = || {
        let p = PathBuf::from(a.req_path("path")?);
        Ok::<_, anyhow::Error>(p.canonicalize().unwrap_or(p))
    };
    let what = a.req_str("action")?;
    match what {
        "list" => {}
        "add" => {
            let pdf = path()?;
            anyhow::ensure!(pdf.is_file(), "{}: no such file", pdf.display());
            let files: Vec<PathBuf> = std::iter::once(pdf.clone())
                .chain(
                    ["png", "svg"]
                        .map(|e| pdf.with_extension(e))
                        .into_iter()
                        .filter(|f| f.is_file() && *f != pdf),
                )
                .collect();
            let note = a
                .str("what")
                .map_or_else(|| file_name(&pdf), str::to_string);
            log.add(&files, &note, None, gallery::Trace::default())?;
        }
        "unlist" => log.remove(&path()?, false)?,
        "delete" => {
            anyhow::ensure!(
                a.parse::<bool>("confirm")?.unwrap_or(false),
                "delete removes the figure's files for good: give `confirm: true` once the user agrees"
            );
            log.remove(&path()?, true)?;
        }
        "move" => {
            let to = a.req_path("to")?;
            log.relocate(&path()?, Path::new(&to))?;
        }
        other => anyhow::bail!("`action` is list, add, unlist, move or delete, not {other:?}"),
    }
    let now = gallery::now();
    let entries: Vec<Value> = (log.entries().iter().zip(log.check()))
        .map(|(e, status)| {
            json!({
                "path": e.path,
                "what": e.what,
                "when": e.when,
                "ago": gallery::ago(e.when, now),
                "status": status.as_str(),
                "figure": e.extra.panel,
                "manifest": e.extra.manifest,
                "files": e.files.iter().map(|f| &f.path).collect::<Vec<_>>(),
            })
        })
        .collect();
    Ok(vec![Content::Json(
        json!({ "log": absolute(log.dir()), "entries": entries }),
    )])
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
fn decisions_text(a: &Arguments, round: &Path) -> Result<(String, Vec<Map<String, Value>>)> {
    let decisions: Vec<Map<String, Value>> = a.parse("decisions")?.unwrap_or_default();
    anyhow::ensure!(!decisions.is_empty(), "no decisions");
    let mut text = String::new();
    for mut d in decisions.iter().cloned() {
        d.insert("round".into(), json!(file_name(round)));
        d.entry("decided_by")
            .or_insert_with(|| json!(DecidedBy::Agent.as_str()));
        text += &format!("{}\n", Value::Object(d));
    }
    Ok((text, decisions))
}

/// `lupin relabel --preview` with the decisions in a temporary file beside
/// the round; the preview's JSON. A preview writes nothing, so it stays out
/// of the driver script.
fn relabel_preview(a: &Arguments, cx: &Call) -> Result<Value> {
    let round = round_arg(a)?;
    let (text, _) = decisions_text(a, &round)?;
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
    let (text, decisions) = decisions_text(a, &round)?;
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
    // A mixed label, checked now, is named with its parts once the round is
    // written.
    let mut mixes = Vec::new();
    for d in &decisions {
        if let (Some(name), Some(parts)) = (d.get("label").and_then(Value::as_str), d.get("parts"))
        {
            let parts: Vec<String> = serde_json::from_value(parts.clone())
                .with_context(|| format!("reading `parts` of {name}"))?;
            crate::tui::ontology::Mixed::check(name, &parts)?;
            mixes.push((name, parts));
        }
    }
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
    for (name, parts) in &mixes {
        crate::tui::headless::remember_mix(&run::parent_dir(&round), name, parts)?;
    }
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
    if let Some(l) = a.str("level") {
        argv.push(flag("level", l));
    }
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
    let level = (loaded.manifest.trajectory.settings.as_ref())
        .and_then(|s| s.get("level"))
        .cloned()
        .unwrap_or(Value::Null);
    Ok(json!({
        "run": absolute(&loaded.file),
        "level": level,
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
