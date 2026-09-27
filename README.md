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
record their outputs back in the manifest they read. Every command that writes files
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

## Related crates

[`senna-rs`](https://crates.io/crates/senna-rs) (train / embed / layout; produces the runs lupin reads),
[`legume-plot`](https://crates.io/crates/legume-plot) (render primitives),
[`data-beans`](https://crates.io/crates/data-beans),
[`legume-numeric`](https://crates.io/crates/legume-numeric).
