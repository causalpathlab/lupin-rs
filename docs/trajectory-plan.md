# Plan: supervised trajectory — an explicit prior, checked against the data

## Status

Not implemented. lupin 0.3.0 removed the previous trajectory stack (`lineage`,
`pseudotime`, `dyn-assoc`, `lineage-plot`); `dev/trajectory-salvage.md` in the
repository records what it did and why it went. This document is the design of
its replacement, `lupin trajectory`.

---

## 1. Why the prior has to be explicit

The removed `lineage` rooted its trajectory at a node the user named, usually
through the marker annotation. With no velocity table (no senna command writes
one), every edge was undirected, so the root alone set pseudotime and branch
order: changing it gave nearly uncorrelated pseudotime on the same embedding.
The output restated the annotation it was rooted on, behind k-means, an MST and
principal curves that made it look data-driven.

The prior cannot be avoided: "which cell type comes first" is biology, not
something a transcriptome embedding settles on its own. So the redesign
**states the prior openly, lets the user supervise it, and reports how far the
data supports it**:

- **Nodes are the user's annotated cell types**, not k-means centroids.
- **Edges are a partial order** ("A precedes B"), assembled from the Cell
  Ontology and the user's own statements, every edge tagged with its source.
- **Every prior edge is tested** against the cell kNN graph, and strong data
  connections the prior lacks are reported as candidates.
- **Cells are ordered by diffusion pseudotime** within the prior's paths.

## 2. Command

```
lupin trajectory -f run.senna.json -o out [--prior FILE [--prior-only]] [--root TYPE]...
                 [--tui] [--label-cl FILE] [--obo FILE]
                 [--knn 15] [--n-dcs 15] [--min-cells 20] [--check-only]
```

- **Labels** are the given round's `annotate.argmax`, so curation in later
  rounds (merges, relabels) is what the trajectory sees.
- **Small types**: types with fewer than `--min-cells` cells are not nodes;
  their cells are treated like unassigned cells (below), and prior edges
  through them are joined across (§3).
- **Unassigned cells** stay in the cell graph — they are often the cells
  between two types — but carry no type: they are not tested in §4 and get
  pseudotime only when connected to a component's typed cells.
- `--root TYPE` names a source of the prior; naming a type that has an
  incoming edge is an error.
- `--prior-only` uses the `--prior` file alone: no ontology, no other layers.
- `--check-only` stops after the connectivity check (§4) and writes no
  pseudotime.

## 3. The prior

Statements come from two sources and are combined into one graph over the
run's types.

1. **Cell Ontology.** The ontology parser keeps `relationship: RO:0002202`
   (`develops_from`; a few hundred in cl-basic) besides `is_a`. Labels map to
   CL terms as in annotation (aliases, `--label-cl`). B develops from A when
   CL says so directly, through terms off the panel, or through an `is_a`
   ancestor of B (CL coverage is sparse: B cell has no `develops_from` of its
   own). Inherited statements are tagged `cl-inherited`.
2. **`precedence.tsv`**, `from<TAB>to<TAB>relation<TAB>note` with `relation`
   `precedes` or `unrelated`. It is a data file like `cl_aliases.tsv`, found
   along the same search path: the user's config, then the project's `lupin/`
   directory beside the run, then a file named by `--prior` for this run only
   (a named file that does not exist is an error, not a warning). The TUI
   (§6) appends to the project layer, as remembered aliases are.

**Combining.** Statements are about a pair of types. A later layer's
statement about {A, B} replaces every earlier statement about that pair,
whatever its direction: `B precedes A` in the project file reverses an
ontology `A → B`; `unrelated` removes it. The order is ontology, then user,
then project, then `--prior`.

**From statements to edges.** The combined statements are closed
transitively over all panel types (so `A → M → B` still orders A before B when
M is dropped for size), restricted to the run's node types, then reduced to
the **direct edges**: an edge is dropped when a longer path implies it. This
removes the shortcuts inheritance creates (CL's "leukocyte develops from
haematopoietic stem cell" would otherwise put HSC → B next to
HSC → MPP → CLP → pro-B → B). Only direct edges are tested in §4 and define
lineages in §5.

The result must be acyclic; a cycle is an error naming its statements and
their sources (`cl`, `cl-inherited`, `user`, `project`, `run`).

## 4. Checking the prior against the data

The connectivity check does not depend on the root.

- **One kNN graph** on the run's geometry table, prepared as annotation's
  Leiden step does (log-simplex exponentiated; cosine for an embedding,
  z-scored otherwise). It keeps each cell's own neighbour list and distances
  (`k` neighbours counting the cell itself, as scanpy's `n_neighbors` does),
  and its edge weights are the Gaussian kernel of §5. Both §4 and §5 use it.
- **Connectivity of types a, b**: the kernel weight between their cells as a
  fraction of the smaller type's total weight, so the value does not scale
  with the number of types on the panel.
- **Null**: random relabelling of cells with type sizes fixed. Its mean and
  variance for a pair are computed in closed form from the graph, giving a
  z-score and p-value without sampling; Benjamini–Hochberg over all pairs.
- **Verdicts**: a direct prior edge is `supported` when q < α and the fraction
  is at least `f_min`, `weak` when only one holds, `unsupported` otherwise; a
  pair outside the prior meeting both is a `candidate`. The defaults for α and
  `f_min` are set in phase 0 from the bone-marrow data and reported with every
  run.

## 5. Ordering: diffusion pseudotime

Diffusion pseudotime (Haghverdi et al. 2016), reproducing scanpy 1.10's
`neighbors(method='gauss', knn=True)`, `diffmap` and `dpt` so the two can be
compared directly:

- **Kernel**: `W(x, y) = sqrt(2σ_x σ_y / (σ_x² + σ_y²)) · exp(−d² / (σ_x² + σ_y²))`
  with σ_x² the median squared distance to x's neighbours, symmetrised;
  density-normalised `K = W / q qᵀ`; then `S = D^-1/2 K D^-1/2`. (destiny uses
  a k-th-neighbour bandwidth, so a small gap to destiny is expected.)
- **Diffusion components**: the `n_dcs` largest eigenpairs of the sparse,
  symmetric `S` by a Lanczos-type solver, as scanpy's `eigsh`; `S` need not be
  positive semi-definite, so a singular-value solver, which loses an
  eigenvalue's sign, is not used. The eigenvectors are those of `S`, not
  rescaled, as in scanpy.
- **Distance**, as scanpy computes it: over all components, eigenvalues below
  0.9994 contribute `(λ / (1 − λ))² (ψ(x) − ψ(y))²` and those at or above it
  `(ψ(x) − ψ(y))²` unweighted; the square root of the sum. Cells in another
  connected piece of the cell graph are at infinite distance.
- **Where it runs**: per connected component of the prior, on that
  component's typed and unassigned cells, then per connected piece of that
  cell subgraph — a prior component whose types never touch in the data is
  several pieces, and each gets its own pseudotime and a warning naming the
  break. A thin bridge (eigenvalue near 1) is reported, not silently weighted.
- **Root cell**: among the root type's cells, the one with the largest mean
  diffusion distance to the cells of the component's terminal types —
  farthest from every fate, and independent of eigenvector signs. With several
  sources in a component, each has a root cell and a cell's pseudotime is its
  distance to the nearest one. Pseudotime is divided by its largest finite
  value, as in scanpy.
- **Lineages**: each root-to-leaf path of the direct-edge graph; a cell's
  lineage weights are uniform over the paths through its type. Cells of types
  outside the prior get no pseudotime (NaN), counted in the log.
- **Order agreement**: for a prior edge A → B, the fraction of B's cells
  beyond A's median pseudotime. Unlike §4 this depends on the root, and is
  reported as such.

## 6. Supervising the prior in the TUI

`lupin trajectory --tui` (also reachable from `annotate --tui`) adds an
**order view** in the tree pane: the run's types with their cell counts, the
direct edges with their sources, and — once a check has run — each edge's
verdict. Mark type A, mark type B, then `>` for "A precedes B" or `-` for
"unrelated" (`x` already stops a running pass anywhere in the TUI); a short
reason is asked for, as for other edits. Saving appends to the project layer's
`precedence.tsv`; annotation rounds (`decisions.jsonl`) are not touched.

## 7. Outputs and the manifest

| output | contents |
|---|---|
| `{out}.trajectory_prior.tsv` | the combined statements and the direct edges, each with its source — a complete prior on its own |
| `{out}.trajectory_edges.parquet` | type pairs: weight, fraction, z, p, q, in-prior, verdict, order agreement |
| `{out}.cell_pseudotime.parquet` | per cell: `pseudotime`, `type`, `component`, `piece`, one weight column per lineage |
| `{out}.diffusion.parquet` | cells × diffusion components |

A `trajectory` section goes into the manifest annotation would write
(`{out}.senna.json`, or `{out}.lupin.json` for a pinto or lupin input), asked
about before any existing file is replaced: the outputs above, the settings,
and each `precedence.tsv` read with its layer and a content hash. Because
`precedence.tsv` files only grow, a rerun that should reproduce this prior
reads `{out}.trajectory_prior.tsv` itself (`--prior … --prior-only`), not the
recorded paths.

`lupin plot --colour-by pseudotime` colours cells from `trajectory.pseudotime`
through a continuous-colour path (the categorical one cannot), with an
optional overlay of the direct edges as arrows between type medians, width by
connectivity, unsupported edges faded. A manifest with no `trajectory` section
is an error when the flag is given, and a warning-and-default when it comes
from `defaults.colour_by` — pre-0.3.0 runs may still carry
`colour_by = "pseudotime"` and an old `pseudotime` block, which is never read.

## 8. Validation

lupin's diffusion pseudotime must match established implementations before it
is trusted.

- **References**: scanpy 1.10 `tl.dpt` (primary) and R `destiny::DPT` (the
  original authors'), run from scripts under `dev/trajectory-bench/` with their
  own environments.
- **Data**: 10x 10k bone-marrow mononuclear cells (`10k_BMMNC_5pv2`), with
  existing senna topic / VAE / SVD runs, annotated with the BoneMarrowMap
  marker panel. Haematopoiesis gives a known order.
- **Like for like**: same cells, same transformed geometry, same `k`
  (counting self), same number of diffusion components, same root cell.
  scanpy runs on each component's cells with its kNN graph rebuilt on them;
  lupin's subgraph keeps the full-data neighbours, so the comparison also
  reports lupin with the graph rebuilt per component, to separate the two
  effects.
- **Metrics**: per-cell Spearman ρ between lupin and scanpy (target ≥ 0.98),
  per-component correlation of the diffusion maps, per-type median-pseudotime
  ranks, run time.
- **Biology**: HSC → MPP → committed progenitors → mature types recovered; a
  deliberately wrong prior (reversed root, a false edge) flagged as
  `unsupported` or by low order agreement.
- **Golden test**: a ~500-cell subset with scanpy's pseudotime kept as a test
  fixture, so CI checks the match without Python.

Unit tests cover `develops_from` parsing, layer replacement and reversal,
transitive reduction (including through dropped types), the cycle error,
connectivity and its closed-form null on a synthetic two-branch set, and DPT
on a synthetic Y shape and on a disconnected one.

## 9. Not in the first version

- RNA velocity.
- Automatic root discovery.
- Association along the trajectory. It comes later as `trajectory assoc`,
  lifting the removed spline-GAM and Bayesian contrast/trend tests; of their
  results only a trend's sign depends on the root.

## 10. Phases

0. Reference baseline: annotate the 10k BMMNC run, run scanpy (and destiny)
   DPT through the bench harness, record the numbers, the golden fixture and
   the default α and `f_min`.
1. `develops_from` parsing, `precedence.tsv` along the search path, combining
   and reduction, and `--check-only` with the connectivity check.
2. Diffusion pseudotime, outputs, manifest section; validated against phase 0.
3. TUI order view.
4. `plot --colour-by pseudotime` and the edge overlay.
5. Later: association.
