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

`lupin relabel --watch` keeps running and turns each batch of lines appended
to the decisions file into the next round (`{out}.r1`, `{out}.r2`, ...), so
decisions can come from a viewer, an agent or an editor while it runs.
`{out}.relabel_status.json` names the round it started from (`base`), every
round since, the latest, and the last refused batch; a restarted watcher
resumes from it. A decision may name the round its ids refer to (`round`);
one naming any round but the latest is refused, so a decision made on a
stale view never lands on renumbered clusters.

```sh
lupin relabel --watch -f r0 -d decisions.jsonl -o rounds/run
```

## Related crates

[`senna-rs`](https://crates.io/crates/senna-rs) (train / embed / layout; produces the runs lupin reads),
[`legume-plot`](https://crates.io/crates/legume-plot) (render primitives),
[`data-beans`](https://crates.io/crates/data-beans),
[`legume-numeric`](https://crates.io/crates/legume-numeric).
