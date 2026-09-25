# lupin-rs

**L**abel, **U**nfold, **P**lace, **I**nterpret, **N**arrate — text graphs, cell-type
annotation, lineage, pseudotime, association, plots, and short descriptions.

Everything ships as one crate: `lupin-rs` (binary `lupin`). Annotation, lineage,
and gene-text logic live as internal modules — not separate crates.io packages.

Lupin reads the runs senna writes (`run.senna.json` and the tables it lists) but
does not link against senna. `-f/--from` takes the manifest file or its output
prefix, and commands that produce artifacts record them back in the manifest,
preserving every field lupin does not use. Every command that writes files
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
lupin annotate -f run.senna.json -m markers.tsv -o out
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
