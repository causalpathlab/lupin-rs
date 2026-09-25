# lupin-rs

**L**abel, **U**nfold, **P**lace, **I**nterpret, **N**arrate — text graphs, cell-type
annotation, lineage, pseudotime, association, and short citation-checked descriptions.

Workspace crates: [`legume-annotate`](annotate/) (lib `annotate`), [`legume-lineage`](lineage/) (lib `lineage`), [`legume-gene-text`](gene-text/) (lib `gene_text`), and the `lupin` binary (`lupin-rs` on crates.io).

## Installation

```sh
cargo install lupin-rs
```

Requires Rust 1.91+ ([rustup](https://rustup.rs/)). Optional GPU backends: rebuild with `--features cuda` or `--features metal`. HDF5 inputs: `--features hdf5`.

## Quick start

```sh
lupin text-qc --uniprot-tsv human.tsv --obo go-basic.obo -o run
lupin word-graph --uniprot-tsv human.tsv --obo go-basic.obo -o run
lupin annotate -f run.senna.json -m markers.tsv -o out
lupin lineage -f out/gem -o out/lin
lupin pseudotime -f run.senna.json -o out
lupin describe -f out --text-prefix run -o out
```

See [`lupin/README.md`](lupin/README.md) and [`gene-text/README.md`](gene-text/README.md) for command details.

## Related crates

Plotting and embedding stacks come from crates.io: [`senna-rs`](https://crates.io/crates/senna-rs), [`legume-plot`](https://crates.io/crates/legume-plot), [`data-beans`](https://crates.io/crates/data-beans), [`legume-numeric`](https://crates.io/crates/legume-numeric).
