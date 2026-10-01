# lupin-rs

**L**abel, **U**nfold, **P**lace, **I**nterpret, **N**arrate — text graphs, cell-type
annotation, lineage, pseudotime, association, plots, and short descriptions.

Everything ships as one crate: `lupin-rs` (binary `lupin`). Annotation, lineage,
and gene-text logic live as internal modules — not separate crates.io packages.

Lupin reads the runs senna writes (`run.senna.json` and the tables it lists) but
does not link against senna. `-f/--from` takes the manifest file or its output
prefix. `lupin annotate` leaves that manifest untouched and writes a new one,
`{out}.senna.json`: a copy of the run with its paths rebased onto the new
location, every field lupin does not use preserved, and an `annotate` section
added, so viewers open the annotated run directly. `lineage` and `pseudotime`
record their outputs back in the manifest they read.

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
lupin lineage -f out/gem -o out/lin
lupin pseudotime -f run.senna.json -o out
lupin plot --from run.senna.json -o out/plot
lupin describe -f run.senna.json --text-prefix run -o out
```

## Cell-type labels

Spaces, commas and underscores in a cell-type label are interchangeable:
`CT 1, a`, `CT_1_a` and `CT,1 a` name one type, written `CT_1_a` in every
output. Marker panels split on the tab (or, with none, on a line's first
comma), so labels may contain commas.

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
agent. Decisions go in a JSONL file, each with its evidence, the alternatives
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
cluster's label agrees). `--preview` uses the same rescoring when it can.

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
another ontology. A term is scored when 20 to 500 of its genes are features
of the data (`--go-min-overlap`, `--go-max-overlap`); the overlap counts, not
the term's own size. Each cluster's top
terms, ranked by their effect (the genes' mean in the cluster against the
other clusters), are in `{out}.ontology_signature.tsv` and in its round
summary entry. Every term is also tested as the cell types are: fgsea's NES,
its p-value against random gene sets and sample permutations, and q by BH
over the terms within each cluster (`{out}.cluster_term_{nes,p,q_values}.parquet`).
Testing thousands of terms takes minutes where the cell types take seconds.
GO terms need `--method enrichment`.
In `lupin annotate --tui`, `GO terms` in the settings (`r`) turns `--go` on
for the next pass; once a round has terms, a third column lists the selected
cluster's top terms with their effect, p and q beside its cell types (`tab`
reaches it).

## Cell Ontology data

lupin places panel labels on the Cell Ontology using three data files, kept out
of the binary so they can be updated and amended without rebuilding:

| file | what |
|---|---|
| `cl-basic.obo` | the Cell Ontology (downloaded once, then cached) |
| `cl_matching.json` | matching rules: which synonyms count (abbreviations such as HSC, GMP), plurals, word order, which ontology subsets are classes |
| `cl_aliases.tsv` | curated `label<TAB>CL:id<TAB>note` mappings for names matching cannot settle (Azimuth's `CD14 Mono`, `Prog Mk`, …) |

Each is looked up in layers, later ones winning (rules key by key, aliases row
by row): the install's `share/lupin/` (or `LUPIN_DATA_DIR`) or the source's
`data/`, else a cache filled by download; your `~/.config/lupin/`
(`LUPIN_CONFIG_DIR`); a project `lupin/` folder beside the run manifest; and
`--label-cl` / `--obo` for one run. Without a rules file, matching is literal
(names and exact synonyms). Each pass records the files and ontology release it
used under `annotate.settings.enrichment.cell_ontology`.

```sh
lupin data where -f run.senna.json   # which file each layer contributes
lupin data fetch                     # cache everything (GO files too), for offline machines
```

In `lupin annotate --tui`, `o` in the tree pane switches to the Cell Ontology
itself: browse a term's parents and children, `/` to search names, synonyms and
abbreviations, `Enter` to label a cluster with any term. When the cluster's top
candidate has no term, lupin offers to remember the pick in the project's
`lupin/cl_aliases.tsv`. `?` lists every key.

What you curate in the TUI is kept as plain files in `lupin/` beside the run
manifest (read after `~/.config/lupin/`, which holds the same names for every
project): `cl_aliases.tsv` (label → CL term), `hidden_genes.txt` (genes or `*`
patterns such as `MT-*` kept out of the specific-genes view) and
`mixed_labels.tsv` (a name for a mixed label, e.g. `HSPC mix<TAB>EMP<TAB>HSC`, given
to a cluster the evidence cannot split).

## Method write-ups

`lupin docs` lists them; `lupin docs <topic>` prints one (compiled into the
binary): `annotation`, `grouping`, `ontology-plan`, `rooting-plan`.

## Related crates

[`senna-rs`](https://crates.io/crates/senna-rs) (train / embed / layout; produces the runs lupin reads),
[`legume-plot`](https://crates.io/crates/legume-plot) (render primitives),
[`data-beans`](https://crates.io/crates/data-beans),
[`legume-numeric`](https://crates.io/crates/legume-numeric).
