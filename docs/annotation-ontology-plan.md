# Plan: reference-free cluster annotation on a gene-annotated Cell Ontology DAG

## Status

Partly implemented. What `lupin annotate` does today:

- **Ontology walk** (`--obo`, with `--label-cl` or an automatic label map): places each cluster
  on the Cell Ontology `is_a` tree from the marker-enrichment scores with TreeBH, at the deepest
  node the data supports, abstaining on sibling ties. Writes `{out}.ontology_assignment.tsv` and
  `{out}.ontology_node_mass.parquet`.
- **Coarse first round**: groups the panel's types under their nearest shared analysis class in
  the ontology (or by shared markers) and calls clusters at that level first.

What remains a plan is the mechanism below: placing clusters by **coverage of their own top genes
on a gene-annotated ontology**, rather than by the scores of a user's marker panel.

---

## 1. The problem

Per-cell and per-cluster annotation posteriors are **flat**: the top cell type often wins with a
small share of the probability mass over a panel of a few dozen types. Calls are confident only
for transcriptionally distinct types; confusable types within one lineage split the mass. This
affects both the enrichment and the projection scorers, so it sits upstream of either.

The root cause:

- The dilution is a **collinearity / overlap** problem. Confusable types share lineage markers,
  and any **symmetric** weighting (IDF `ln(C/df)`, empirical specificity) credits the shared genes
  equally to all of them, so it **structurally cannot break the tie**.
- For some pairs the **discriminating evidence is not in the RNA assay at all**: it lives in
  receptor sequence, surface protein or chromatin. No RNA-only method can recover it.

**So the goal is not a sharper number.** The number *should* be flat when the evidence is flat.
The goal is:

> **Report the finest level of the cell-type hierarchy at which the call is identifiable; abstain
> ("not enough resolution") below that; and surface what the data contains that the ontology
> cannot explain.**

---

## 2. Core idea

Treat annotation as **placing activity marks on a gene-annotated Cell Ontology (CL) DAG to explain
each cluster's top-k genes**, under a submodular, information-theoretic objective.

- **Explaining-away.** A coarse mark covers shared lineage genes once, so a leaf mark only pays off
  if it has *unique* present genes: submodular diminishing returns. This removes the dilution
  structurally.
- **Adaptive granularity.** Marks land only at the depth that earns marginal coverage; otherwise
  the walk stops and reports the parent.
- **Abstention** is the refusal to break a sibling tie.
- **Disagreement** is the residual top-k genes no node explains.

---

## 3. Objects

- **DAG** `G = (V, is_a)`: CL restricted to the tissue. It is a DAG, not a tree: many terms have
  more than one parent.
- **Gene annotations** with the **true-path rule**: a gene annotated to a node is implicitly
  annotated to all its `is_a` ancestors. `genes(t)` = genes at `t` or below.
- **Information content** `IC(t) = −log(|genes(t)| / |genes(root)|)`, from the ontology's own
  annotation statistics: the specificity weight IDF could never be, and parameter-free.
- **Gene → mark affinity** by semantic similarity (Resnik): `a(g, t) = IC(MICA(d(g), t))`, where
  `d(g)` is gene `g`'s deepest annotation(s) and MICA is the most informative common ancestor. A
  gene annotated in two lineages contributes to both, which correctly *creates* a cross-lineage tie
  for the tie rules to resolve.
- **Per cluster** `c`: top-k genes `G_c` with weights `w_{c,g}`. **k and the ranking are part of the
  model** (§7).

---

## 4. Algorithm

### Phase A: depth (resolution walk), which emits "not enough resolution"

Top-down from the root; descend only while a single child dominates the evidence.

```
walk(node v, genes G_c):
    children = is_a-children of v supported by G_c
    for each child u: mass(u) = Σ_{g∈G_c} w_g · a(g, u)
    normalise mass over children → P(child | v);  h = entropy(P(child | v))
    u* = argmax mass
    if P(u*|v) ≥ τ and h ≤ h_max:     # one child dominates
        return walk(u*, G_c)           # commit deeper
    else:
        return STOP(v, candidates = {u: mass(u) high}, residual_entropy = h)
```

The output is the deepest unambiguous node, plus an explicit "resolution-limited, unresolved among
{children}" when the walk stops early; `h` at the stop is the confidence of the abstention.

### Phase B: multiplicity and disagreement (submodular residual)

After the primary path, explain the *uncovered* top-k genes; whatever survives is the
disagreement output.

```
R = { g ∈ G_c : a(g, v*) low }
while R not empty and gain ≥ gmin:
    t* = argmax_t Σ_{g∈R} w_g · a(g, t)   # facility-location marginal gain
    place mark t*; remove explained genes from R
report the remaining R as data-only residual   # a program the ontology cannot explain
```

### Objective choices (all submodular)

- **Selection: information-weighted facility location** `f(S) = Σ_g w_g · max_{t∈S} a(g, t)`:
  soft, discrimination-aware, provably submodular, no independence assumption.
- **Depth: conditional entropy / mutual information**, submodular under naive Bayes (genes
  independent given type). Shared genes contribute almost no entropy reduction automatically.
- Coverage alone is **discrimination-blind** (covering a shared gene scores like covering a
  discriminative one); do not use it as the only objective.
- Constraint: a **DAG antichain matroid** (a cluster maps to an antichain; no ancestor and
  descendant marked together unless both earn gain), so greedy keeps its guarantees.

Propagation and IC are computed once; per cluster the work is cheap and parallel.

---

## 5. Tie-breaking: most ties are not broken

A tie is the signal, not a nuisance. By location in the DAG:

1. **Sibling tie (same parent): do not break; this is the abstention.** Report the parent
   ("subtype unresolved"). `τ` is the only knob.
2. **Ancestor vs descendant: Occam plus necessity, default shallower.** Prefer the more specific
   node only if its extra specificity is free:
   `score(t) = explained(t) − λ · (annotated_in(t) but absent_in G_c)`. Still tied: go shallower.
3. **Cross-lineage tie (different subtrees): back off to the LCA, or report multiplicity.** An
   informative LCA is reported; two lineages each independently supported by different genes is a
   doublet or mixed state (Phase B finds it); an uninformative LCA is flagged "cross-lineage
   ambiguous".
4. **Exact numerical ties: deterministic, for reproducibility only.** A fixed key (higher prior,
   higher IC, lower CL id), kept apart from the semantic rules so it never makes a biological
   decision.

**Unifying principle: bias every unresolved tie toward claiming less** (shallower, higher LCA, or
abstain).

---

## 6. Outputs

- **Assignment**: a CL node per cluster at adaptive depth, with its ancestor path.
- **Abstention**: an explicit "resolution-limited" flag, the candidate children and the residual
  entropy.
- **Disagreement**: top-k genes no node covers ("the ontology cannot explain"), and clusters that
  accept only a coarse mark ("definitionally distinct, empirically unresolved in this assay").
- **Provenance**: per assignment, the genes that grounded it (with IC), covered vs residual.
- **Annotation depth**: per node, how many genes support it; a trust cap (§10).

---

## 7. The k / gene-ranking interaction

The ranking that defines "top-k genes" decides which objective works:

- top-k by **expression magnitude** is dominated by shared lineage genes, so coverage points coarse
  and the discriminative signal never enters the budget;
- top-k by **specificity** (differential against the rest) lets discriminative genes in.

Take top-k by a specificity-aware score, let the entropy objective decide whether those genes
actually disambiguate, and abstain if not.

---

## 8. Decision support in a viewer

The `τ` and `k` knobs become interactive, with their consequences shown on the embedding:

1. Cells coloured by their node at the resolved depth; abstained cells in the parent's colour,
   desaturated, so low-resolution regions show as desaturated blobs.
2. A residual-entropy layer: where the calls are ambiguous.
3. A linked DAG panel: clicking a node highlights the cells at or below it; a `τ` slider re-walks
   the DAG and recolours.
4. A disagreement view: clusters with a large uncovered residual.
5. A per-cluster evidence card: top-k genes, the node each supports (with IC), covered vs residual,
   and the candidate siblings at the stop node.

---

## 9. Prior art

The pieces exist individually; the combination does not.

**Abstention at internal nodes (not new on its own):**
- CellO (Bernstein et al., *iScience* 2021): hierarchical classification on CL, places cells at
  internal nodes when uncertain; supervised, trained on a reference atlas.
- Uncertainty-aware annotation with a hierarchical reject option (2023/24): partial and full
  reject; supervised.
- GPTAnno (2025): ontology-guided, automatic resolution selection; LLM-based.
- Hierarchical cross-entropy loss (2025): atlas-scale, visualised on the CL DAG.

**Explaining-away on an ontology DAG (in gene-set enrichment, not cell typing):**
- MGSA (Bauer, Gagneur, Robinson, *NAR* 2010): a Bayesian network selects a minimal set of
  categories that explain a gene list, accounting for overlap.
- topGO elim/weight (Alexa et al. 2006): DAG-conditional enrichment; evoGO (2025): redundancy
  minimisation.

**Submodular selection in genomics:**
- Submodular assay-panel selection (Wei, Libbrecht, Bilmes, *Genome Biology* 2016): facility
  location.
- scGeneFit (*Nat. Commun.* 2021): label-aware marker selection with hierarchical labels.

**Marker-based hierarchical typing:**
- Garnett (Pliner et al., *Nat. Methods* 2019): a user-specified hierarchy and markers, with an
  "unknown" class.
- OnClass (Wang et al., *Nat. Commun.* 2021): embeds the CL graph; zero-shot to unseen terms.

**Where this plan is new:**
1. **Reference-free**: driven by the ontology's own gene annotations and IC, with no labelled
   reference.
2. **Submodular coverage / conditional entropy** as the explicit assignment objective on a
   gene-annotated CL DAG, with greedy guarantees and an antichain matroid.
3. **Disagreement as a first-class output**: what the data has that the ontology lacks.
4. **Abstention as refusing to break sibling ties**, with one "claim less when tied" rule rather
   than a threshold.

The abstention contribution alone is well trodden; lead with the reference-free mechanism and the
disagreement output.

---

## 10. Open questions and risks

**Check before building:**
- **Method vs information**: do the discriminating genes for confusable pairs exist in the gene →
  CL annotations and vary in the data? If not, the method correctly abstains everywhere it
  matters; know that up front.
- **Annotation coverage** of the gene → CL map. Native CL gene axioms are sparse and may need a GO
  bridge or a marker-database projection. **Annotation completeness, not the algorithm, caps
  resolution**; instrument per-node support from the start.

**Decisions:**
- One blended objective vs **two stages** (facility location selects the antichain, entropy walks
  depth); start two-stage, which is easier to inspect.
- The gene → node annotation source: native CL axioms, a GO bridge, or a marker database; likely a
  union with provenance per edge.
- Calibrating `τ` and `h_max` against types that should stay deep.
- Query the DAG via the lowest common ancestors of the label set rather than flattening it.

**Risks:**
- It creates no information: an uncovered residual stays uncovered (correct, but set
  expectations).
- `a(g, t)` quality is everything: weak annotations give weak affinities.
- Greedy is approximate; with the necessity penalty, use distorted greedy (Harshaw et al.).

---

## 11. Where it would live in lupin

- Inputs: a cluster × gene top-k table (from `{out}.cluster_expression.parquet`), the Cell
  Ontology (found as for `lupin annotate`) and gene → CL annotations.
- Core: DAG propagation and IC, then Phase A/B per cluster in parallel; outputs as parquet
  (assignment, path, entropy, covered, residual), recorded in the round's manifest like the other
  annotation outputs.
- It reuses the existing ontology walk's tree handling and the first round's ontology resolution.
