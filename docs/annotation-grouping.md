# Why marker annotation pools cells into *coarse* clusters

The projection method (`lupin annotate --method projection`) calls each group of cells by its
most over-represented marker type. This note explains why those groups must stay coarse, and what
does and does not reduce false calls. The scoring code is `type_annotation::term_ora` in the
`legume-graph-embedding` crate.

## The constraint

Term-ORA calls a cluster by its **most over-represented** marker type. Over-representation is a
*discovery* statistic: it ranks terms by how **surprising** a count is, not by how **likely** the
term is. Those two rankings agree only when the cluster is large.

- In a cluster of several hundred cells you need many cells of a type before the count is
  surprising, so most-enriched and most-abundant coincide. The distinction is invisible.
- In a group of a few dozen cells it inverts. A type with only a handful of cells in the entire
  dataset has an expected count near zero, so two of them are more "enriched" than the many cells
  of the type that actually fills the group.

Anything that shrinks the groups walks into this: a high `--resolution`, or replacing the
partition with per-cell neighbourhoods.

## The negative control

A marker panel can carry types the tissue is known not to contain (for example mature or terminal
types in a progenitor-rich sample). Commenting them out is the normal setting; enabling them turns
every cell they win into a *countable false positive*, with no ground truth needed. We call the
fraction of cells they capture the **absent-type share**, and compare groupings by it.

On such a control we compared the untested per-cell `argmin`, Leiden communities at a low
resolution, per-cell kNN neighbourhoods (small k, and k grown to Leiden's own community size), and
each of those with the panel bootstrap (since removed; see `lupin docs annotation` §4). The findings:

- **Leiden at a low resolution cut the absent-type share to well under half of the untested
  `argmin`**, and assigned no cells at all to the types `argmin` used only sparsely.
- **The clustering-free variant was worse than the untested `argmin` it was meant to filter.**
  Small per-cell neighbourhoods raised the absent-type share and manufactured cells for types with
  almost none in the dataset. Growing the neighbourhood to Leiden's pooling size narrowed the gap
  but did not close it.
- Calling by **plurality** instead of by enrichment (the vote picks the label; the p-value only
  gates whether to call at all) removed the manufactured rare-type cells, but left the absent-type
  share high. Those cells were a symptom; the remaining false positives are coherent *blobs* in the
  embedding, and a small neighbourhood can sit entirely inside one. Smoothing cannot remove them.

## Why Milo's device does not transfer

Milo (Dann et al., *Nat. Biotechnol.* 2022) tests kNN neighbourhoods for **differential
abundance**, and it works because its labels (condition, sample) are **external to the
embedding**. Ours are `argmin` over marker centroids in the *same* space that defines the
neighbourhoods, so a per-cell neighbourhood test partly re-tests the geometry against itself. That
is the disanalogy; do not lean on the precedent.

## What actually moves the number

Not the grouping. On the same negative control:

- **The panel bootstrap** (resample the marker panel *and* re-derive the clustering, ship the
  consensus) reduced the absent-type share by an order of magnitude, at the cost of abstaining on
  a substantial fraction of cells. It has since been removed as too stringent and too slow
  (`lupin docs annotation` §4).
- **The embedding.** When most types have fewer than half their markers trained, no grouping and no
  statistic can rescue a centroid built from genes the model never saw. The loss is type-dependent:
  a highly-variable-gene filter rewards variance, and a rare population's markers are high-variance
  by construction. senna's gem now trains every gene by default (`--n-hvg 0`); if you set
  `--n-hvg > 0`, pass `senna gem --markers <panel>` so the panel is trained regardless.

None of the groupings recovered most of the types the tissue actually contains.
The grouping question is second-order until the embedding and panel coverage are fixed.
