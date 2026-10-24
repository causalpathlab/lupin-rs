# lupin-rs

**L**abel, **U**nfold, **P**lace, **I**nterpret, **N**arrate — text graphs, cell-type
annotation, trajectories, and short descriptions. Figures are drawn in each
tool's own viewer (`senna view`, `pinto view`, lupin's TUI).

Everything ships as one crate: `lupin-rs` (binary `lupin`). Annotation
and gene-text logic live as internal modules — not separate crates.io packages.

Lupin reads the runs senna writes (`run.senna.json` and the tables it lists) but
does not link against senna. `-f/--from` takes the manifest file or its output
prefix. `lupin annotate` leaves that manifest untouched and writes a new one,
`{out}.senna.json`: a copy of the run with its paths rebased onto the new
location, every field lupin does not use preserved, and an `annotate` section
added, so viewers open the annotated run directly.

`lupin annotate` also reads a pinto run (`run.pinto.json`, or its prefix): its
count files, the final level's `cluster` column in the propensity table, and,
for `cage`, the shared cell and feature embeddings, which it annotates by
projection. Cell coordinates, when the run had them, become a `spatial`
layout. The pinto manifest is never written; the result is a new
`{out}.lupin.json`, with the same layout as a senna manifest. Senna does not
read it; lupin commands take it with `-f`. Every command that writes files
takes an explicit `-o/--out` prefix; lupin never writes to a location derived
from the manifest, so a run copied to another machine works as is.

## Installation

```sh
cargo install lupin-rs
```

Requires Rust 1.91+. Optional: `--features cuda` / `--features metal` / `--features hdf5`.

## Quick start

```sh
lupin text-qc --uniprot-tsv human.tsv --obo go-basic.obo -o run
lupin word-graph --uniprot-tsv human.tsv --obo go-basic.obo -o run
lupin annotate -f run.senna.json -m markers.tsv -o out   # writes out.senna.json
lupin trajectory -f out.senna.json -o out.T1           # prior, check, pseudotime
lupin trajectory                                       # the TUI: pick a run, state the order, run
lupin describe -f run.senna.json --text-prefix run -o out
lupin describe -f out.T1.senna.json -o out.T1          # a trajectory's order in words
```

## Sessions

In an interactive terminal the TUI (`lupin annotate` without `-o`, `lupin
trajectory` without both `-f` and `-o`) runs as a detachable session, so a
pass, a save or a trajectory run goes on when the connection drops.
ctrl+alt+d detaches (the first time, it asks for a name), ctrl+alt+s
switches to another session, and starting the TUI again with no `--from`
offers the live sessions to return to. `lupin view <manifest>` opens the
screen on an output: a trajectory run on the order view, a round on the
round view. `lupin attach`, `lupin sessions`,
`lupin rename` and `lupin kill` manage them; `kill` quits as `q` does,
without asking: a save under way finishes, running jobs stop and unsaved
edits are dropped. `--no-detach` (or `LEGUME_TUI_NO_DETACH=1`) runs it in
the terminal only.

## Cell-type labels

Spaces, commas and underscores in a cell-type label are interchangeable:
`CT 1, a`, `CT_1_a` and `CT,1 a` name one type, written `CT_1_a` in every
output. A marker panel line with a tab splits on its tabs, so a label may
contain commas; a line without one is CSV, where a label with a comma must be
quoted (`CD4,"CD4-positive, alpha-beta T cell"`).

## First round: high-level calls

A marker pass calls each cluster by a broad group of the panel's cell types,
with the evidence of the group's types added up; the fine calls stay in the
cluster summary (and each cell's fine label in `annotate.fine_argmax`) to
refine in later rounds. `--fine` calls the fine types directly.

The groups come from the Cell Ontology: `--obo`, else a copy cached under the
user cache directory, else a download into that cache (`LUPIN_OFFLINE=1`
skips it). Panel labels are matched to terms by name or exact synonym,
ignoring case, and each type goes under its nearest analysis class (the ontology's upper slims
and `cellxgene_subset`) shared with another panel type. With the ontology
found, its walk runs by default with the matched labels (`--label-cl`
overrides). Without it, types that share marker genes are grouped instead.
`{out}.celltype_tree.json` records the groups, their source and the
ontology release.

## Without a marker panel: ask an AI

```sh
lupin ask -f run.senna.json -o a0 --context "tissue, species"   # prints a prompt
# paste the prompt into any AI chat, then paste its answer back:
lupin relabel -f a0.senna.json -d - --next < answer.txt
```

`lupin ask` writes a first round with unlabelled clusters and a prompt listing
each cluster's most specific genes (log fold change of counts per 10k over the
other clusters, among the genes it expresses strongly). The prompt asks for
decision lines only: labels with a rationale and the genes relied on, and
marker sets for each label. `relabel` reads the answer as pasted (prose and
code fences are skipped), records who decided (`agent_proposed_user_accepted`)
and rescores the round against the suggested markers. lupin sends nothing
anywhere.

## Annotation rounds

An annotated manifest is a round. `lupin review` prints each cluster's evidence
(candidate labels with q and support, top GO/GMT terms, Cell Ontology placement)
and every decision made on it so far; `--json` gives the same to a program or an
agent, with `cl` the term the round's type tree puts the cluster's label on. Decisions go in a JSONL file, each with its evidence, the alternatives
weighed and a rationale; `lupin relabel` applies them and writes the next round:

```sh
lupin annotate -f run.senna.json -m markers.tsv -o r0
lupin review -f r0                       # or: --json, -c 3
lupin relabel -f r0 -d decisions.jsonl -o r1
```

Each round records `annotate.source` (the round before it),
`annotate.cluster_summary` (`{out}.cluster_summary.json`, keyed by cluster id),
`annotate.log` (this round's decisions) and `annotate.history`
(`{out}.annotation_history.json`: every round's decisions per cluster, newest
first). Merged clusters take fresh ids, so an id always names the same cells'
history. `lupin review --help` documents the decisions format.

Rounds after the first form a chain: `r0.r1`, `r0.r2`, ... beside `r0`, and
the latest is the highest one on disk. Two ways grow it, and they can run at
once; every writer takes the chain's lock, and a decision made on anything but
the latest round is refused ("reload and decide again"):

- `lupin relabel -f <round> -d - --next` reads decisions on stdin, writes the
  next round and prints its path. A viewer runs it once per decision.
- `lupin relabel --watch -f <round> -d decisions.jsonl` keeps running and turns
  each batch of lines appended to the file into the next round, so an agent or
  an editor can decide while it runs. Each watched decision names the round it
  was made on (`round`). `{chain}.relabel_status.json` lists the rounds, the
  latest (including rounds written by direct calls) and the last refused batch;
  a restarted watcher resumes from it.

```sh
echo '{"cluster": 3, "action": "label", "label": "CT1", "rationale": "...", "decided_by": "user"}' \
  | lupin relabel -f r0 -d - --next        # prints r0.r1.senna.json
lupin relabel --watch -f r0 -d decisions.jsonl
```

`--preview` takes the same input and checks as `--next` but writes nothing: it
prints, as JSON, each affected cluster's label before and after and the cells
that would change. When the round records a cluster expression profile
(enrichment), cell types are also re-ranked against the edited marker panel by
a marker module score, a quick approximation of what the next `lupin annotate`
would call.

Rounds of an enrichment run are rescored. The pass caches its sufficient
statistics (per-cluster gene sums, the per-batch profile, gene weights and
each cell's batch; `annotate.stats_cache`), and every relabel round reruns the
same scoring on its merged clusters and edited marker panel without reading
the counts: fresh q-values every round. These statistics come after curation
on the same data, so the round records
`annotate.stats = {kind: post_selection, rounds_of_curation}`; a fresh
`lupin annotate` pass is the confirmatory one. Each cluster's summary
entry carries `evidence` (the top call, its q, and whether the
cluster's label agrees: it is one of the calls tied at that q, which go by
the method's own reading, then the NES). `--preview` uses the same
rescoring when it can.

A round is never overwritten, whichever way it is written.

Decisions can also revise the marker panel: `markers_add` / `markers_drop`
name a cell type (`label`) and `features`. The round then writes its own
`{out}.markers.tsv` (the previous panel with the edits) as `annotate.markers`,
and `{out}.marker_history.json` keeps each cell type's edits with their
rationale. `lupin annotate -f <round>` without `-m` re-annotates from that
panel.

## GO terms

`--go` scores GO terms per cluster on the same cluster
expression profile, next to the cell-type calls:

```sh
lupin annotate -f run.senna.json -m markers.tsv --go -o r0   # cell types and GO terms
lupin annotate -f r0.senna.json --go -o r1                   # the round's panel, plus GO
```

The species is told from the gene names (Ensembl id prefixes, else symbol
case; human and mouse), and its GO annotations (`goa_human.gaf.gz`,
`mgi.gaf.gz`) and the Gene Ontology (`go-basic.obo`) are downloaded once from
the GO Consortium into the cache below, looked up in the same layers as the
Cell Ontology data. `--gaf` or `--gmt` names other gene sets, `--go-obo`
another ontology. A GAF's annotations by an author's untraced statement
(NAS) or with no data (ND) are left out, through a filtered copy kept in the
cache: one NAS paper puts the cytosolic ribosome under "negative regulation
of myoblast fusion", which would top every ribosome-rich cluster.
The terms sharing 20 to 500 genes with the data are tested
and reported (`--go-min-overlap`, `--go-max-overlap`; the overlap counts, not
the term's own size): the window is fixed before anything is scored, the
broadest terms tell no cluster apart, and q-values are over the terms in it.
Each cluster's top
terms, ranked by their effect (the genes' mean in the cluster against the
other clusters), are in `{out}.ontology_signature.tsv` and in its round
summary entry. Every term's effect is also tested, against 200 random gene
sets per term with as many genes from each of 20 expression bins (a ribosome
meets other highly expressed genes): z, p (normal tail) and q (BH over the
terms within each cluster) in `{out}.cluster_term_{z,p,q_values}.parquet`.
Random genes do not move together as a term's genes do, so the decoys'
spread is widened by how much more the term varies over the clusters than
they do (by the median absolute deviation, which one cluster standing out
does not move): a programme high in many clusters (ribosome, spliceosome,
oxidative phosphorylation) tests weak in each, one that sets a cluster apart
does not. A cell type split over many clusters shares its programmes among
them in the same way, so they test weaker in each than their effect reads.
On a PBMC run, 0.4% of the terms per cluster (at most 4%) clear q < 0.05,
the NK, platelet, B, dendritic and monocyte clusters by their own
programmes, in seconds.
A relabel round that merges clusters rescores the terms on the
new clusters from the pass's cached statistics; one that keeps the clusters
keeps the terms (a pass that cached no statistics leaves a merged cluster
without terms until a new pass, and says so). A preview of a merge (`lupin
relabel --preview`) lists the merged cluster's top terms with their effect
and q, so the merge is weighed with them. A merged cluster's
genes are read as its parts' cell-weighted mean.
GO terms need `--method enrichment`.
In the annotate TUI (`lupin annotate` without `-o`), `GO terms` in the annotation form (`r`) turns `--go` on
for the next pass; once a round has terms, a third column lists the selected
cluster's top terms with their effect, p and q beside its cell types (`tab`
reaches it).

## Cell Ontology data

lupin places panel labels on the Cell Ontology with two data files, kept out
of the binary so they can be updated and amended without rebuilding:

| file | what |
|---|---|
| `cl-basic.obo` | the Cell Ontology (downloaded once, then cached) |
| `cl_matching.json` | matching rules: which synonyms count (abbreviations such as CT1 for its full name, plural forms), plurals, word order, marker signs (`X+`, `X positive` and `X-positive` read alike), qualifiers a label may leave off (`alpha-beta`, a species), a word a whole panel leaves off its labels (a tissue word), placing a label nothing matches by the longest run of its words that a term names, which ontology subsets are classes and which terms are states that group nothing |

Each is looked up in layers, later ones winning (rules key by key): the
install's `share/lupin/` (or `LUPIN_DATA_DIR`) or the source's `data/`, else a
cache filled by download; your `~/.config/lupin/` (`$XDG_CONFIG_HOME/lupin/`
when that is set, `LUPIN_CONFIG_DIR` over both); a
project `lupin/` folder beside the run manifest; and `--obo` for one run.
Without a rules file, matching is literal (names and exact synonyms).

Labels that matching cannot settle are mapped by **aliases**,
`label<TAB>CL:id<TAB>note` rows, kept as data rather than code. They are
read row by row, later layers winning:

1. the shared abbreviations, `data/cl_aliases.tsv` in this repository,
   downloaded from the release's tag and cached (`lupin data fetch`); not
   packaged with the crate, and skipped offline when nothing is cached;
2. the panel's sidecar, `x.cl.tsv` beside a panel `x.tsv.gz`;
3. `cl_aliases.tsv` in `~/.config/lupin/`, then in the project's `lupin/`;
4. `--label-cl` for one run.

To change a shared row, map the label in any later layer. Each pass records
the files and ontology release it used under
`annotate.settings.enrichment.cell_ontology`.

```sh
lupin data where -f run.senna.json   # which file each layer contributes
lupin data fetch                     # cache everything (GO files too), for offline machines
lupin data suggest-cl CT1 CT2        # each label's nearest Cell Ontology terms, as JSON
```

In the annotate TUI, `o` in the tree pane switches to the Cell Ontology
itself, opening on the terms in the data (`d` toggles the whole ontology):
browse a term's parents and children, `/` to search names, synonyms and
abbreviations, `Enter` to label a cluster with any term. In the genes pane `h`
hides a gene and `X` hides by pattern; `x` always stops a running job. `g`
lists the run's family (the run, its annotated rounds and its trajectories)
to switch between them without restarting. When the cluster's top
candidate has no term, lupin offers to remember the pick in the panel's
sidecar (else the project's `lupin/cl_aliases.tsv`). `?` lists every key.

For a panel type with no term at all, `c` in the tree pane suggests the
three terms nearest it by text embedding, with their scores
(`lupin data suggest-cl`, PubMedBERT sentence embeddings at a pinned
revision, fetched once into the Hugging Face cache, `HF_HOME`, else
`~/.cache/huggingface`; `LUPIN_OFFLINE` or `HF_HUB_OFFLINE=1` rules the
download out). A suggestion is never taken on its own: the one you confirm is
written to your `~/.config/lupin/cl_aliases.tsv` and the tree is built again
with it. The first ranking fetches the model and encodes the ontology, which
takes a minute or two on a CPU; the index is cached under lupin's cache
directory per ontology release and model revision. Offline without the
model there are no suggestions, and the rest works as before.

What you curate in the TUI is kept as plain files in `lupin/` beside the run
manifest (read after `~/.config/lupin/`, which holds the same names for every
project): `cl_aliases.tsv` (label → CL term, after the panel's sidecar), `hidden_genes.txt` (genes or `*`
patterns such as `MT-*` kept out of the specific-genes view) and
`mixed_labels.tsv` (a name for a mixed label, e.g. `CT1 mix<TAB>CT2<TAB>CT3`, given
to a cluster the evidence cannot split).

## Agents (MCP)

`lupin mcp` serves lupin's annotation and trajectory loops over the
[Model Context Protocol](https://modelcontextprotocol.io) on stdin/stdout, so
any MCP client can drive them, with any model behind it (Claude, GPT, Gemini,
or a local one). lupin computes and states the facts; the agent supplies the
judgments, each with a rationale that lupin keeps.

### Connect a client

The server is one command, `lupin mcp`. A client needs that, and where the
project is: the server keeps its records (below) in its working directory, or
in `LEGUME_MCP_PROJECT` when set. Set it whenever the client may start servers
somewhere else (desktop apps usually do). If the client cannot find `lupin`,
give its full path (`which lupin`).

Most clients (Claude Desktop, Cursor's `.cursor/mcp.json`, Gemini CLI's
`settings.json`, and others) take this shape:

```json
{
  "mcpServers": {
    "lupin": {
      "command": "lupin",
      "args": ["mcp"],
      "env": { "LEGUME_MCP_PROJECT": "/path/to/project" }
    }
  }
}
```

Claude Code:

```sh
claude mcp add lupin -e LEGUME_MCP_PROJECT="$PWD" -- lupin mcp
```

VS Code (`.vscode/mcp.json`):

```json
{ "servers": { "lupin": { "type": "stdio", "command": "lupin", "args": ["mcp"],
                          "env": { "LEGUME_MCP_PROJECT": "${workspaceFolder}" } } } }
```

Codex CLI (`~/.codex/config.toml`):

```toml
[mcp_servers.lupin]
command = "lupin"
args = ["mcp"]
env = { LEGUME_MCP_PROJECT = "/path/to/project" }
```

Key names differ between clients; the server needs only the command, its
argument and, optionally, the environment. Nothing is sent anywhere by lupin:
the model sees what the tools return, and where that goes is up to the client
and model you choose (a local model keeps it on the machine).

### What the agent does

The server tells the agent how to work when it connects (the MCP
`instructions`). Some clients do not pass those to the model; with one of them,
put this in the agent's own instructions (an `AGENTS.md`, a system prompt):

```text
lupin annotates single-cell clusters with cell types and orders the types into trajectories. Its tools state facts computed by code; you make the judgments, each with a written rationale that lupin keeps.

Annotation loop:
1. annotate_run on a run manifest (senna or pinto) writes a new round; or, without a marker panel, annotate_ask writes an unlabelled round and lists each cluster's specific genes for you to label. With one of `go` (the species' GO annotations), `gaf` or `gmt` (a pathway file), it also tests GO or pathway terms per cluster: evidence of what a cluster does, beside its cell type.
2. annotate_review shows each cluster: its calls with q, its top terms, each method's top type, its Cell Ontology placement and its history. annotate_genes adds the genes: what sets each cluster apart, and how the panel markers of its label and of each method's top type fare in it. Where the methods disagree, weigh their evidence, terms included; do not label beyond it. A cluster's terms rank by effect, what sets it apart (NK-cell killing, platelet clotting, MHC-II presentation); their q tests that effect against random gene sets and discounts programmes many clusters share (ribosome, spliceosome), so a cell type split over several clusters tests weaker in each than its effect reads: weigh effect and q together. A round's `caveats`, when it has any, say what its statistics lack: read them first.
3. annotate_preview with your decisions shows what they would change, rescored; then annotate_relabel applies them and returns the next round. For a merge, the preview lists the merged cluster's top terms with their effect and q, so weigh whether the parts share a function before merging. A merge that would leave one cluster is refused (each cluster is scored against the others): label the clusters the same type instead. annotate_relabel's `warnings` say what the new round lacks, such as terms not scored again. Rounds are never overwritten, so a decision can always be revisited.

Trajectory loop:
1. trajectory_run on an annotated round, with roots and any precedence statements (from, to, precedes|unrelated, note), writes a run.
2. trajectory_read returns its lineages, each edge's connectivity and verdict (supported, unsupported, candidate, inferred) and the order in words. Revise the statements where the data disagree and run again.

annotate_run and trajectory_run take `background`: on a large run set it, and the call returns a job id at once; the step runs on lupin's own screen in a session the user can watch with `lupin attach`, and goes on if this server ends.

docs prints the method write-ups. Paths are resolved against the server's working directory; results give the paths written.
```

With records kept (the default), legume-mcp adds where they are and asks the
agent to open with a `notebook` plan (the question, the data, how the result
will be checked) and to close with its conclusions. Records or not, it also
tells the agent how to follow background jobs.

A large `annotate_run` or `trajectory_run` can outlast a client's wait on one
call: with `background: true` it returns `{job, round}` or `{job, run}` at
once, and the agent follows the job with `job_wait` or `job_status` (or stops
it with `job_cancel`), using the manifest only once the job is done. A job that
succeeds is recorded like any step; jobs end with their server. While a job
runs, its output is held: a call with no `out` skips its name, and one given
it is refused, naming the job, before anything is written.

| Tool | Does |
|---|---|
| `annotate_run` | a pass (`lupin annotate`) on a run manifest; returns the round and each cluster's label and each method's top type |
| `annotate_ask` | without a panel: an unlabelled round and each cluster's most specific genes |
| `annotate_review` | per cluster: calls with q, terms, each method's top type, Cell Ontology placement, history; `types` adds every type's evidence |
| `annotate_genes` | per cluster: the genes that set it apart, and how the panel markers of its label and of each method's top type fare |
| `annotate_preview` | what decisions would change, rescored; writes nothing |
| `annotate_relabel` | applies decisions as the next round of the chain; rounds are never overwritten |
| `trajectory_run` | orders an annotated round along a prior (roots, precedence statements) and the kNN graph |
| `trajectory_read` | lineages, each edge's connectivity and verdict, the order in words |
| `docs` | lupin's method write-ups |
| `notebook` | (from legume-mcp) an entry in the day's lab notebook |
| `job_status`, `job_wait`, `job_cancel` | (from legume-mcp) follow or stop a background job |

Decisions use the `relabel` format (`lupin review --help`). They are recorded
as `decided_by: agent` unless the agent says a person accepted them
(`agent_proposed_user_accepted`). Passes run as child `lupin` processes. No
tool replaces an existing manifest: an output the agent does not name takes the
next free `{stem}.L{k}` / `{stem}.T{k}` in the day's results folder.

### Records

The server keeps the project's records, laid out after Noble's guide to
organizing computational biology projects [1] (see `legume-mcp`):

- `results/<YYYY-MM-DD>/runall.sh`: every step that wrote something, in order,
  with paths relative to the project and steps skipped when their outputs
  exist, so `bash results/<day>/runall.sh` in a copy of the project redoes the
  work without the agent. A relabel runs `--next --expect {chain}.r<k>`, so it
  replays to the same round.
- `results/<YYYY-MM-DD>/notebook.md`: the agent's plan, observations and
  conclusions; the plan names the question, the data and the check, after [2].
- `results/INDEX.tsv`: one row per output with a one-line summary, for `grep`.

`LEGUME_MCP_RECORDS=off` keeps none.

### Watching

While an agent works, the server writes a live feed to
`$XDG_RUNTIME_DIR/legume-mcp/lupin-<pid>/`, and after each round or trajectory
it writes, it publishes it (`{"round": …}` or `{"trajectory": …}`) as the view
to open. A person follows and steers it (stop a call, leave the agent a note)
with legume-mcp's `watch` reader and view. `LEGUME_MCP_WATCH=off` writes none.

[1] W. S. Noble, "A quick guide to organizing computational biology projects",
PLoS Comput Biol 5(7): e1000424 (2009). doi:10.1371/journal.pcbi.1000424

[2] W. S. Noble, "Ten simple rules for defining a computational biology
project", PLoS Comput Biol 19(1): e1010786 (2023). doi:10.1371/journal.pcbi.1010786

## Method write-ups

`lupin docs` lists them; `lupin docs <topic>` prints one (compiled into the
binary): `annotation`, `grouping`, `ontology-plan`, `trajectory-plan`.

## Related crates

[`senna-rs`](https://crates.io/crates/senna-rs) (train / embed / layout; produces the runs lupin reads),
[`legume-plot`](https://crates.io/crates/legume-plot) (render primitives),
[`data-beans`](https://crates.io/crates/data-beans),
[`legume-numeric`](https://crates.io/crates/legume-numeric).
