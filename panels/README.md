# lupin-panels

Marker panels for [lupin](https://github.com/causalpathlab/lupin-rs), shipped
as data in their own crate so that lupin's code carries no cell-type
vocabulary. lupin reads a panel by name: `lupin annotate --markers
panel:<name>`.

## A panel

Each panel is a directory `data/<name>/` with:

| File | Contents |
|---|---|
| `<name>.tsv` | the marker table: `gene<TAB>cell type`, one marker per line; a `gene` header and `#` lines are skipped |
| `<name>.cl.tsv` | optional: labels mapped to Cell Ontology terms, `label<TAB>CL:id<TAB>note`, for labels name matching cannot settle |
| `README.md` | where the markers come from, how they were built, and the licence or terms they are shared under |

## Adding a panel

1. Check that the source's licence lets the markers be redistributed, and
   how it must be credited; say so in the panel's README.
2. Add `data/<name>/` with the three files above.
3. Add the panel to `PANELS` and `NAMES` in `src/lib.rs`.
4. Bump the crate's version (a new panel or a changed one is a new version),
   and run `cargo test -p lupin-panels`.

A run records the panel's name and this crate's version, so it can be
reproduced with the same panel.
