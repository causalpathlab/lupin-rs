# Plan: supervised trajectory, an explicit prior checked against the data

## Status

Partly implemented. `lupin trajectory` builds the prior (§3), checks it
against the kNN graph (§4), orders cells by diffusion pseudotime from the
prior's roots (§5), writes its outputs and a `trajectory` manifest section
(§7); without `-f` and `-o` it opens the TUI on the order view, with its
figures and export log
(§6). lupin draws its figures only there: it has no plot commands, each
tool having its own viewer. lupin 0.3.0 removed the previous trajectory stack
(`lineage`, `pseudotime`, `dyn-assoc`, `lineage-plot`);
`dev/trajectory-salvage.md` in the repository records what it did and why it
went.

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

The numerical core is not new: lupin **replicates scanpy's unsupervised
routines** (`neighbors(method='gauss')`, `diffmap`, `dpt` and `paga`)
faithfully, and adds only the supervised, interactive layer on top: the groups
are the user's types, the root and the edges come from the prior, and the user
edits the prior where the data disagrees. Where scanpy's routine breaks on an
input, lupin stops with a clear message instead of a home-grown fix.

## 2. Command

```
lupin trajectory -f run.senna.json -o out [--prior FILE [--prior-only]] [--root TYPE]...
                 [--labels FILE] [--label-cl FILE] [--obo FILE]
                 [--knn 15] [--n-dcs 10] [--min-cells 20] [--min-connectivity 0.1]
                 [--check-only] [--graphics auto|kitty|sixel|iterm2|blocks]
lupin trajectory [-f run.senna.json] [-o out] [options]   # the TUI
```

With both `-f` and `-o` the command runs and exits. With either missing it
opens the TUI on the order view: a manifest is picked in its file browser
when `-f` is not given, and `r` runs the trajectory as a child process with
the options given, asking for the output prefix (the last one used, else the
next free `{stem}.T{k}`; an existing one is replaced only on a second
Enter). A run that is already annotated opens with its round in the
cluster panes. On a run with no annotation, `r` asks where the labels come
from, in a menu that explains each choice: annotate the run first, or a
`cell<TAB>type` labels file picked in the file browser. Annotating opens
the **annotation form**, which asks for what a pass needs before it starts:
the marker panel (the manifest's when it can be found, else required, Enter
picks it in the file browser), the output prefix (the next free
`{stem}.L{k}`), the method and the clustering and permutation settings. Its
first line says what running it will do (which round it writes, or what it
replaces, and that the trajectory follows); the last row, `▶ run`, starts
the pass (Shift+Enter does too where the terminal tells it from Enter). A
pass with no panel is refused with the reason in the form, and one that
replaces a round asks to be run again there. `A` in the clusters pane
re-annotates an annotated run through the same form; after that pass the new
round opens and a menu offers the trajectory on its labels. When nothing orders the types (no statement gives an
edge, and no `--root` or `--prior`), `r` asks how the order should come
instead of starting a run that can only fail: start from one type
(`--root`, picked from the types with their cell counts and kept for later
runs), state the order in the table, a precedence file (`--prior`), or a
`label<TAB>CL:id` file (`--label-cl`) so the ontology's develops-from links
order the types (asked again if they still do not). In every menu ↑↓
choose, Enter or the option's key takes it, and esc cancels. `lupin annotate` follows the same rule: without `-o` it opens its
TUI, where `r` (or `A` in the clusters pane) opens the same annotation
form.

- **Labels** are the given round's `annotate.argmax`, so curation in later
  rounds (merges, relabels) is what the trajectory sees; `--labels` names
  another `cell<TAB>type` file. A cell missing from the labels counts as
  unassigned (reported); labels matching no cell at all are an error.
- **Small types**: types with fewer than `--min-cells` cells are not nodes;
  their cells are treated like unassigned cells (below), and prior edges
  through them are joined across (§3).
- **Unassigned cells** stay in the cell graph (they are often the cells
  between two types) but carry no type: they are not tested in §4 and get
  pseudotime only when connected to a component's typed cells.
- `--root TYPE` (repeatable) names a source of the prior; a type that a
  node type precedes cannot be one (a type too small to be a node does not
  count). A root gets a `cli` statement to every
  node it cannot otherwise reach, so a run with no other prior still orders
  every type from it (the check then says which of those edges the data
  supports).
- `--prior FILE` is a precedence file for this run, or a
  `{out}.trajectory_prior.tsv` from an earlier run to replay; `--prior-only`
  uses it alone: no ontology, no user or project layer.
- `--check-only` stops after the connectivity check (§4) and writes no
  pseudotime.

## 3. The prior

Statements come from two sources and are combined into one graph over the
run's types.

1. **Cell Ontology.** The ontology parser keeps `relationship: RO:0002202`
   (`develops_from`; a few hundred in cl-basic) besides `is_a`. Labels map to
   CL terms as in annotation (aliases, `--label-cl`). B develops from A when
   CL says so directly, through terms off the panel, or through an `is_a`
   ancestor of B (CL coverage is sparse: CT3 may have no `develops_from` of
   its own). Inherited statements are tagged `cl-inherited`.
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
then project, then `--prior`. Within a file the later line wins; the
ontology's statements have no order, so if it says both `A → B` and `B → A`
both are kept and reported as a cycle. A `# kind` header marks a
`{out}.trajectory_prior.tsv`, whose rows start with `statement` or `edge`.

**From statements to edges.** The combined statements are closed
transitively over all panel types (so `A → M → B` still orders A before B when
M is dropped for size), restricted to the run's node types, then reduced to
the **direct edges**: an edge is dropped when a longer path implies it. This
removes the shortcuts inheritance creates (an inherited "CT5 develops from
CT1" would otherwise put CT1 → CT5 next to
CT1 → CT2 → CT3 → CT4 → CT5). Only direct edges are tested in §4 and define
lineages in §5.

The result must be acyclic; a cycle is an error naming its statements and
their sources (`cl`, `cl-inherited`, `user`, `project`, `run`). So is an
`unrelated` pair that other statements still order (A → M → B with
`A unrelated B`): it names the path and its statements, since `unrelated`
removes only a statement about that pair.

## 4. Checking the prior against the data

The connectivity check does not depend on the root.

- **One kNN graph** on the run's geometry table, prepared as annotation's
  Leiden step does (log-simplex exponentiated; rows L2-normalised for an
  embedding, so Euclidean distance is `√(2 − 2 cos)`: it ranks neighbours as
  cosine does and is the distance scanpy's kernel uses; z-scored otherwise).
  Each cell's neighbour list and distances come from legume-numeric's
  `knn_rows_l2` (exact) or `knn_rows_ivf` (IVF), switched at
  `knn_graph::ALL_PAIRS_THRESHOLD` (65,536 cells). Both leave the cell itself
  out, so lupin asks for `k − 1` and counts the cell in, as scanpy's
  `n_neighbors` does; the cell is not its own neighbour, and the kernel has no
  diagonal. Exact duplicate cells make scanpy drop a true neighbour for a
  self-loop; lupin instead collapses each set of identical cells onto one
  representative, builds the graph over those, and gives every copy its
  representative's pseudotime and diffusion coordinates, with a warning
  (how many cells, and how many coincide with a cell of another type; the
  representative's type is the one counted in connectivity). Edge weights are the Gaussian kernel of §5. The graph is built once over all cells, as scanpy
  builds it on a whole dataset, and both §4 and §5 use it.
- **Connectivity of types a, b**: scanpy's PAGA (`tl.paga`, connectivity
  model v1.2), which counts the kNN edges between groups against a random
  null. It counts each cell's own directed list of `k − 1` neighbours, so an
  edge both cells list counts twice and a group's edge total is the sum of
  its cells' list lengths. Every cell is in a group, as scanpy requires: the
  node types, and the unassigned and small-type cells as groups of their own,
  which take part in the statistic but get no verdict.
- **Verdicts**: a direct prior edge is `supported` when its connectivity
  reaches `--min-connectivity` and `unsupported` below it; a pair of nodes
  joined by no path of the prior whose connectivity reaches it is a
  `candidate`. PAGA's spanning tree is not used: with many pairs saturated at
  1 its tied edges follow the order of the type names.
- **What it can and cannot show**: on the bench data every edge of the known
  hierarchy sits above the default threshold and every planted
  edge between separated branches below it, which is how the default was
  set; but a planted edge between two types that touch (siblings under one
  progenitor, a stem cell and a small type beside it) saturates like a true
  one. Connectivity tells touching types from separated ones and no more;
  order agreement (§5) carries the rest, and the report says which check
  flagged an edge. Candidates are many for the same reason, and the log
  shows only the strongest.

## 5. Ordering: diffusion pseudotime

Diffusion pseudotime (Haghverdi et al. 2016), reproducing scanpy 1.10's
`neighbors(method='gauss', knn=True)`, `diffmap` and `dpt` so the two can be
compared directly:

- **Kernel**: `W(x, y) = sqrt(2σ_x σ_y / (σ_x² + σ_y²)) · exp(−d² / (σ_x² + σ_y²))`
  with σ_x² the median squared distance to x's `k − 1` neighbours, symmetrised;
  density-normalised `K = W / q qᵀ`; then `S = D^-1/2 K D^-1/2`.
- **Diffusion components**: 15 eigenpairs (`diffmap`'s `n_comps`) of the sparse, symmetric
  `S` largest in magnitude, as scanpy's `eigsh(which='LM')` finds them, here
  in float64 by legume-numeric's randomised SVD (`rsvd_with`, from its next
  release; 20 power iterations, 10 oversample columns), each eigenvalue's
  sign taken from its Rayleigh quotient `uᵀSu`, then sorted by signed value,
  largest first, as scanpy orders them. At the default 5 iterations the
  clustered leading eigenvalues come out wrong enough that pseudotime departs
  from scanpy's; at 20 and 10 the bench matches `eigsh`. Every pair is checked by its
  residual `‖Su − λu‖`; if any is too large (a flatter spectrum, or a ± pair
  of equal size that the SVD cannot separate), lupin doubles the iterations
  up to a cap and otherwise stops, saying so. The eigenvectors are those of
  `S`, not rescaled, and the eigenvalues are rounded to float32 before the
  0.9994 test below, as in scanpy; one within 1e-4 of 0.9994, where the
  solver's error could put it on the other side, is reported.
- **Distance**, as scanpy computes it: over the first `--n-dcs` components
  (default 10, as `tl.dpt`; all 15 are written and drawn), eigenvalues below
  0.9994 contribute `(λ / (1 − λ))² (ψ(x) − ψ(y))²` and those at or above it
  `(ψ(x) − ψ(y))²` unweighted; the square root of the sum. Cells in another
  connected piece of the cell graph are at infinite distance.
- **Where it runs**: one diffusion map over the whole §4 graph, as scanpy runs
  it. A cell of a prior type belongs to that type's prior component; an
  unassigned cell belongs to the component whose root is nearest. A type with
  no prior edge is in no component: its cells get no pseudotime. Each
  component is measured from its own roots and divided by its own largest
  finite distance, so every component spans 0 to 1. A cell no root can reach,
  in another connected component of the graph, is at infinite distance, as in
  scanpy, and is reported. A region joined by only a few edges gets an extra
  eigenvalue at or above 0.9994, unweighted (above); lupin names such regions
  in the report.
- **Degenerate geometry**: a cell whose neighbours mostly coincide with it has
  σ_x = 0 (after duplicates are collapsed, near-coincident cells can still do
  this), which turns scanpy's kernel NaN and fails its eigensolve (senna VAE
  and topic latents with saturated cells did this on the bench data). lupin
  checks for σ_x = 0 before the eigensolve and stops, naming how many cells
  are affected and suggesting another embedding.
- **Root cell**: the root type's medoid (among its cells in the graph
  component holding most of them, the one with the smallest mean diffusion
  distance to the others), independent of eigenvector signs. On the bench
  data the first draft's rule, the root-type cell farthest from every
  terminal type, picked an outlier and ordered the types far worse than the
  medoid on the same geometry.
  With several sources in a component, a cell's pseudotime is its distance
  to the nearest root.
- **Lineages**: each root-to-leaf path of the direct-edge graph, listed in
  `{out}.trajectory_lineages.tsv`; a cell's lineage weights (`L0`, `L1`, …)
  are uniform over the paths through its type, zero for a cell of no node
  type. Cells of types outside the prior, and cells no root reaches, get no
  pseudotime (NaN), counted in the log.
- **Order agreement**: for a prior edge A → B, the fraction of B's cells
  beyond A's median pseudotime. Unlike §4 this depends on the root, and is
  reported as such.

## 6. Supervising the prior in the TUI

`lupin trajectory` (also `t` in the annotate TUI's tree pane) adds an
**order view**: beside the clusters, the run's types on the Cell Ontology,
drawn as annotate's tree pane draws the panel (how many types map to terms
and how many develops-from links join them above it, the selected type's
place in the ontology, or that it has none, below), and the precedence table:
the run's types with their cell counts, the direct edges with their sources,
and, once a check has run, each edge's verdict. Tab goes clusters, genes,
ontology, precedence; selecting a type in either column selects it in the
other, and `o` in the ontology column opens the Cell Ontology there. Below
150 columns the two share one column, showing whichever has the focus.

The Cell Ontology view (`o`, here and in annotate's tree pane) opens on the
terms **in the data**: the terms the labels sit on (the order view's types,
else the round's labels, else the panel's types), under their lowest common
ancestor, each marked with its types and cells; a chain of terms that only
leads on to one term is folded into that term's line (→ opens it, ← folds
it). `d` switches to the whole ontology around the same term and back. A
search (`/`) lists the matches in the data's tree, or in the whole ontology
when none are there, and says so. Mark type A, mark type B, then `>` for "A precedes B" or `-` for
"unrelated" (`x` already stops a running pass anywhere in the TUI); a short
reason is asked for, as for other edits. Enter on the reason writes the
statement to the project layer's `precedence.tsv` at once (no save step);
annotation rounds (`decisions.jsonl`) are not touched. A new statement about a
pair replaces the earlier one, which is how a statement is corrected. A
statement is refused, with the reason shown and no reason asked, when the
file's last line about the pair already makes it, when it would close a
cycle, or when it says `unrelated` about a pair the rest of the prior orders.

**Figures and what was exported.**

The figure panels take the order table's place in the tree pane while they
are shown (`v` shows them and steps through them, `V` brings the table back);
space, `>` and `-` do nothing then, since types are marked in the table. The
panels are drawn from the run's trajectory outputs, the same set the bench summary shows:

- **Scatter**: the cells on one of the run's layouts, with the prior's
  direct edges as arrows between type medians (supported edges solid, the
  rest faded), or on two diffusion components (`,`/`.` step the y axis, `[`/`]` the x axis).
  `m` steps through the coordinates: senna's PHATE first when the run has
  one (`senna layout phate`, Moon et al. 2019, in legume-numeric's
  `matrix::layout`; `layout.methods.phate`), being made to show
  trajectories, then the current layout, then the other layouts the
  manifest records, then the diffusion map. lupin reads layouts and never
  writes one. As in senna view, `t` steps the cell-type labels at each
  type's median through small, medium (the default), large, largest and
  off, and `c` colours the cells by pseudotime, cell type, lineage (a
  type on one lineage only) or component, the last two only when there are
  several; cells without a value are grey, and a categorical colouring has
  a legend (no more than half the figure high, none in a thumbnail). Every
  group gets its own colour: legume-plot's palette while it lasts, then
  hues a golden angle apart. Unassigned cells get no label.
  `+`/`=` and `-`/`_` zoom the scatter in and out about its centre in
  steps of ×√2 (up to ×64, and back out to whole exactly), the
  arrows pan it (with the exports strip open, ↑ and ↓ move in the strip
  and ← → still pan), and `0` shows it whole again, as in senna view. The figure is
  drawn again from the data at the part on screen, so points and labels stay
  sharp; a label sits at the median of its type's cells on screen, and a
  type with none there has no label. A new layout (`m`), diffusion pair or
  figure (`v`, or opened from the grid) starts whole. The order and
  connectivity panels do not zoom.
- **Order by type**: median pseudotime and middle half per type, in order of
  the median.
- **Connectivity**: PAGA connectivity between the node types as a Hinton
  diagram (box area ∝ connectivity; legume-plot's `render_hinton`), ordered by
  pseudotime, pairs the prior orders in ink.
- **All figures** (`w`, as senna view's grid): each layout the run has, the
  diffusion map at the pair last shown, the order and the connectivity as
  thumbnails, drawn whole with the labels and colouring on screen. The
  arrows choose one, Enter opens it in the pane, `p` exports it as its tile
  shows it, `f` opens the exports strip, `V` goes back to the order table,
  and esc or `w` closes the grid on the figure it came from. `t` and `c`
  also work there.
- The direct prior edges with their connectivity, order agreement and verdict
  are in the order view's table, so an edit there shows its effect after the
  next run.

Images are drawn as senna view draws them, so both viewers behave alike:
`ratatui-image` (same major version, crossterm) with `--graphics
{auto,kitty,sixel,iterm2,blocks}`, querying the terminal on `auto` and falling
back to half-blocks.

**Export.** `p` on a panel writes it through legume-plot's `write_figure` as
SVG and PDF, 7 in wide at 200 dpi, with the labels, colouring, layout and
zoom on screen; the log notes a zoomed export. An export is a set of files sharing one base name, handled as a
unit: the name `{out}.trajectory.{panel}` (`layout_phate`, `layout_umap`, …
for a layout) moves to
`-2`, `-3` … while any file of the set exists.

**What was exported.** Each export is listed in `./.lupin-view/saved.json` in
the directory lupin runs from (newest first; an existing PDF path replaces
its entry; reloaded from disk before every change, and a change finds its
entry by path, not by row). The list is written to a temporary file and
renamed into place, and a list that cannot be read is set aside as
`saved.json.bad` and reported, never treated as empty. Each entry records the
panel, the manifest it came from, and each file's absolute path, size,
modification time and content hash. A strip in the TUI (`f`) shows the list;
an entry can be moved (`M`; nothing is replaced, and a failed move is undone),
removed from the list (`d`), or removed with its files (`D`, which deletes
only on a second `D` on the same export). Exporting (`p`), deleting (`D`) and
moving (`M`) act on the log only when it and its directory are this user's
alone: writable by no one else and not links (lupin creates them so); saving
over a log others can write is refused. Files are moved or deleted only when
they are still the export that was logged: beside the PDF, with the same base
name, and unchanged. A failed log write never fails the export.

**Refresh.** `R` re-reads the list and checks every entry against the files
on disk: `ok`, `changed since export` or `missing`. A file whose size and
modification time match the record is taken as unchanged without hashing;
otherwise it is hashed, so a touched file with the same bytes stays `ok`.
Missing entries stay listed, marked, until removed. Figures beside the run
named `{out}.trajectory.{panel}` but not on the list are shown as
`not listed`; Enter on one adds it. The same check runs when the TUI opens.

## 7. Outputs and the manifest

| output | contents |
|---|---|
| `{out}.trajectory_prior.tsv` | the combined statements and the direct edges, each with its source; a complete prior on its own |
| `{out}.trajectory_edges.parquet` | type pairs: PAGA connectivity, in-prior, verdict, order agreement |
| `{out}.cell_pseudotime.parquet` | per cell: `pseudotime`, `type`, `component`, one weight column per lineage |
| `{out}.diffusion.parquet` | cells × diffusion components |
| `{out}.trajectory_lineages.tsv` | each lineage's component and root-to-leaf path |

A `trajectory` section goes into the manifest annotation would write
(`{out}.senna.json`, or `{out}.lupin.json` for a pinto or lupin input), asked
about before any existing file is replaced: the outputs above and the
settings (`knn`, `n_dcs`, `min_cells`, `min_connectivity`, `roots`, `prior`,
`prior_only`). Because `precedence.tsv` files only grow, a rerun that should
reproduce this prior reads `{out}.trajectory_prior.tsv` itself
(`--prior … --prior-only`), which holds every statement that was used.

## 8. Validation

lupin's diffusion pseudotime must match established implementations before it
is trusted.

- **Reference**: scanpy 1.10 (`dpt` and `paga`, exact neighbours), run
  locally as a test-driven-development tool; none of the bench, its notes or
  its reference fixtures is committed.
- **Data**: a sample with published cell-type labels whose developmental
  order is known, embedded by `senna svd`; the sample itself, the bench
  outputs and the reference fixtures stay out of the repository. Baseline:
  scanpy's pseudotime from the root type's medoid agrees with the expected
  stage order.
- **Like for like**: same cells, same prepared geometry, same `k` (counting
  self), same number of diffusion components computed and used by DPT, same
  root cell; no exact duplicate cells (lupin collapses them, scanpy does
  not).
- **Metrics**: per-cell Spearman ρ between lupin's and scanpy's pseudotime on
  the full bench data (target ≥ 0.98; defined while both searches are exact,
  up to 65,536 cells; scanpy's default goes approximate from 8,192, so the
  bench runs it exact), the diffusion subspaces compared by
  principal angles (single components can rotate within a near-degenerate
  pair, as λ₂ and λ₃ nearly do here), per-type median-pseudotime ranks, run
  time.
- **Biology**: the known stage order recovered from the root type; a
  deliberately wrong prior (reversed root, a false edge) flagged as
  `unsupported` or by low order agreement.
- **Reference tests**: a stratified subset of the bench cells with its
  prepared geometry, scanpy's pseudotime and scanpy's PAGA on it, kept as
  local fixtures that the tests read at run time and skip when absent, so the
  match can be checked without Python and no data enters the repository.

Unit tests cover `develops_from` parsing, layer replacement and reversal,
transitive reduction (including through dropped types), the cycle error,
PAGA connectivity against the reference PAGA fixture, the degenerate-geometry
error, and DPT on a synthetic Y shape and on a disconnected one.

## 9. Not in the first version

- RNA velocity.
- Automatic root discovery.
- Association along the trajectory. It comes later as `trajectory assoc`,
  lifting the removed spline-GAM and Bayesian contrast/trend tests; of their
  results only a trend's sign depends on the root.

## 10. Phases

0. Reference baseline (done): scanpy DPT and PAGA on the bench sample, the
   reference fixtures and the baseline, all kept locally.
1. (done) `develops_from` parsing, `precedence.tsv` along the search path,
   combining and reduction, and `--check-only` with the connectivity check;
   its `--min-connectivity` default was set by planting false edges on the
   bench data.
2. (done) Diffusion pseudotime, lineages, outputs, manifest section; validated
   against phase 0.
3. (done) TUI order view, figures and export log.
4. (done) The layout figure with the edge overlay, in the TUI; lupin's
   plot commands removed (each tool has its own viewer).
5. Later: association.
