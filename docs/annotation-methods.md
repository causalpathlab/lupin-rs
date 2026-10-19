# Marker-based cell-type annotation by projection: methods

What `lupin annotate --method projection` does and why each part of it is there. Every default
quoted here is the shipped default; results are stated qualitatively.

Projection calls cell types on a run's embeddings rather than its raw counts, in one of two
ways (`--projection`):

- **`score`**: the embedding marker score, §0. The default for a run that decodes cells
  through its gene table (bge, gem, simba, resolve-embedding-space).
- **`units`**: the panel's cell types, and GO terms when asked for, trained as gene-set units
  on the run's genes, kept fixed, §1. The default for other runs with a gene embedding, and for
  embedding files given directly (`--feature-embedding`), whose table lupin cannot tell.

## Every method in one round (`--method auto`)

Under `--method auto` (the default) every method the run supports runs, and their readings go
into one round, on the enrichment's clusters:

- the enrichment, with its SuSiE stage where it has count sums (the raw counts here, or an
  earlier pass's cache of them; not a decoder's expectation);
- the projection score (§0) where the run decodes cells;
- the gene-set units (§1) where the run's genes are placed among the cells (its co-embedding,
  or a ρ that shares the cells' space). GO terms asked for (`--go`, `--gaf`, `--gmt`) are then
  both tested by the enrichment and aligned by the units.

A projection beside the enrichment is non-fatal: if it fails, the round goes on without its
reading, with a warning. The round records the score's per-cell table
(`annotate.marker_score`) and the units' (`{out}.units_call.parquet`, `annotate.units_call`:
each cell's cluster, call and membership of its region), with the units' other files (§1).
The round view and `review --json` show every reading, and every type's evidence from every
method goes to `{out}.cluster_evidence.parquet` (`annotate.cluster_evidence`), whose `score`
and `units` columns are each method's share of the cluster's cells called the type.

**Which reading makes the call.** `--call {enrichment, susie, score, units}` picks the method
whose reading becomes the round's per-cell labels (`annotate.argmax`), and so its cluster
labels; the others stay beside it as evidence. Under `--call score` or `units` the round
records no `annotate.annotation`: the enrichment's per-cell posterior is not the call's. By default it is SuSiE's where it ran, else the
enrichment's. A `--call` naming a method that does not run on the pass is refused before
anything is read, and a pass whose calling method fails is an error, not a fallback. The
round records the method as `annotate.called_by`. The TUI's form sets it as `call by`, under
auto.

**One method alone.** `--method enrichment` (what senna and pinto run) runs the enrichment
alone; `--method projection`, or `--projection` under auto, runs that projection alone, with
its own call. Under auto, a run whose enrichment has no expression to read (no raw counts here,
no earlier pass's cache, no decoder to stand in) is projected instead: by the score where it
decodes cells, else by the units.

---

## 0. The embedding marker score (`--projection score`)

Code: `type_annotation::embedding_score` in the `legume-graph-embedding` crate, driven by
`lupin annotate` (`annotate::marker_score`).

**The gene table is the decoder's.** For the runs that decode cells through it (bge, gem,
simba, resolve-embedding-space), a gene's log rate in a cell is `ρ_g · z_i + a_g`, with ρ the
run's `outputs.feature_embedding` and z its cell embedding. The co-embedding is never read: it
places a gene where its direction points, among the cells it is most aligned with, which is not
where its expressing cells are. For a `gem` run only the spliced rows are kept.

**A type's score is a module score read in closed form.** Standardise each gene's logit over the
cells and average a type's markers, weighted by IDF over the panel's types (a marker many types
share counts little). Every step is linear, so the mean moves inside the dot product:
`s_T(i) = (z_i − m) · u_T`, where `u_T` is the weighted mean of `ρ_g / σ_g` over the type's
markers and `σ_g` is the spread of gene g's logit over the cells. The bias `a_g` cancels, and no
gene × cell matrix is formed. This is the module score of Tirosh et al. (2016, *Science*
352:189) as SCANPY's `score_genes` implements it (Wolf et al. 2018, *Genome Biology* 19:15),
read here on the model's implied expression. The control sets are replaced by a reference point
`m`, as described next.

**Centred at the balanced centre.** `m` is the mean of the cluster means: `--clusters`, else the
run's `cluster.clusters`, else internal Leiden. That way a large cluster does not define
"average" for the others. These clusters are also the round's cluster table.

**The call is per cell, with its margin.** Each cell takes its best type. The margin over the
runner-up is kept as the call's confidence: across the bench datasets it is the signal that
separates confident cells from doubtful ones. A per-cell gene-set null and a hierarchical
walk were both tried, and neither raised the precision of what they kept.

**Outputs:**
- `argmax.parquet` (recorded as `annotate.argmax`, `argmax.tsv` beside it): the call; the
  probability is NaN, since the score has none;
- `marker_score.parquet`: each cell's cluster, call and margin, recorded as
  `annotate.marker_score` (not `annotate.annotation`, which holds a cell × type posterior). The TUI tables a cluster's candidates from these calls over the
  round's current clusters, so merges and splits keep their evidence. A marker edit needs a
  new pass: the score caches no statistics to rescore with.

**Guards:**
- IDF weighting and the three-marker floor are part of the score: another `--min-markers` is
  refused;
- a type with fewer than three live markers is dropped;
- fewer than two types left is an error naming the types lost;
- a non-finite embedding is an error.

---

## 1. Gene-set units (`--projection units`)

Code: `annotate::geneset_units`, on `fne` in the `legume-graph-embedding` crate.

**An alignment through the genes.** The cells and the genes already share a space: the
co-embedding for the kinds that co-embed, else ρ where it shares the cells' space. A cell type
(its markers) or a GO term (its annotated genes) is a gene set, and each set is placed in that
space by SGD as a node of a gene–set graph. The engine is PyTorch-BigGraph's dot-product softmax
with sampled negatives and row Adagrad (Lerer et al. 2019, *SysML*), as `fne` implements it. The
gene rows are the run's own and frozen (`senna fne --freeze-feature-embedding` does the same),
so the sets cannot bend the space and few genes cannot overfit it. Averaging the markers'
positions was replaced: an average lands near the middle of the cloud and needs several markers.

**Read per unit of the cells' spread.** Genes and cells are whitened by the cells' pooled
covariance within their clusters, as discriminant analysis does, with cells centred at the
balanced centre of the cluster means. An axis along which every cluster spreads alike (depth,
batch) then counts for little.

**A type is a region.** With the centres fixed, a mixture of isotropic balls on the unit sphere
is fitted to the cells by EM: each type's radius, its share of the cells, and each cell's
membership. A cell is called by the type it most belongs to. Two types are linked by the cells
they share, `Σ_c w_ck w_cl / √(n_k n_l)`.

**GO terms are aligned with clusters against decoys.** `--go`, `--gaf` or `--gmt` add terms
(sharing `--go-min-overlap`…`--go-max-overlap` genes with the table) as units of their own
kind. Random gene sets drawn from the annotated genes, binned by size, are trained in the same
run as decoys. A cluster's score for a term is its mean cell direction against the term's
centre; p is the normal tail of the decoys' scores of that size, and q is Benjamini–Hochberg
within the cluster.

**Outputs** (recorded under `annotate`):
- the per-cell call and membership, as every pass writes them (beside another method's call,
  only the units' own table below holds them);
- `{out}.units_call.parquet` (`units_call`): each cell's cluster, call and membership, which a
  round tables by cluster as it does the score's;
- `{out}.celltype_units.parquet` (`celltype_units`): each type's genes found and given, radius,
  share and centre;
- `{out}.celltype_membership.parquet` (`celltype_membership`): each cell's membership of each
  type's region, cells × types, which the TUI's structure plot stacks;
- `{out}.celltype_links.tsv` (`celltype_links`): the linked types;
- with terms, `{out}.term_alignment.parquet` (`term_alignment`), cluster × term score, p and q,
  and `{out}.term_units.parquet` (`term_units`), each term's centre;
- `{out}.unit_anchors.parquet` (`unit_anchors`): each type's and term's 30 nearest cells by
  cosine (unit, kind, cell, weight), which place them on the TUI's cell map.

`--device cuda` (with `--device-no`) trains on a GPU in a CUDA build.

---

## Retired: the co-embedding term ORA

`--projection ora` (Euclidean nearest-centroid on the marker average, then cluster × type
over-representation) was removed in favour of the units; `--no-idf`, `--no-assign-qc` and
`--assign-mad` are still accepted and ignored with a warning. What it taught, so it is not
re-derived:

- **Per-cell neighbourhoods do not replace clusters.** Over-representation ranks types by how
  surprising a count is, not how likely: in a small neighbourhood the two rankings invert. Full
  write-up: `lupin docs grouping`.
- **Parametric Bayesian centroids fail.** A free per-type variance let one component
  shrink-wrap a blob and take many cells at near-zero entropy: misspecification surfaced as
  false confidence.
- **High Leiden resolution is a bad operating point**: more, less stable communities and fewer
  cells labelled.
