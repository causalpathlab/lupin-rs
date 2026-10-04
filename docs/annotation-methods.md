# Marker-based cell-type annotation by projection: methods

What `lupin annotate --method projection` does and why each part of it is there. Every default
quoted here is the shipped default. Where a design choice is justified by a measurement, the
measurement was made on a marker panel that includes types absent from the tissue (the negative
control of §7); results are stated qualitatively.

Code: `type_annotation::{term_ora, panel_null}` in the `legume-graph-embedding`
crate, driven by `lupin annotate`.

---

## 1. Inputs and the space the call is made in

`senna gem` embeds cells and gene features jointly. Two of its outputs matter here:

| output | rows | meaning |
|---|---|---|
| `{out}.cell_embedding.parquet` | cells × H | the latent cell coordinate **θ_c** |
| `{out}.feature_coembedding.parquet` | feature rows × H | the co-embedded gene vector **e_g** (ρ re-placed onto the cell manifold; the raw ρ is `{out}.feature_embedding.parquet`) |

Feature rows are keyed `{gene}/count/{spliced,unspliced}`; annotation selects the spliced rows and
re-keys them by gene (`spliced` is the base track every gene has; `unspliced`, when present, is its
ridge-shrunk offset track).

**Annotation does not use the `β_g` dictionary**, and this is load-bearing. A Euclidean
nearest-centroid call is only meaningful if genes and cells inhabit one metric space. gem couples
β and θ through an *inner product*, which fixes their relative directions but not their relative
scale, and the fitted model exploits that freedom: β rows are typically orders of magnitude
shorter than cell vectors. At that ratio

    ‖θ − c‖² = ‖θ‖² − 2⟨θ, c⟩ + ‖c‖²

loses the ‖c‖² term, and ‖θ‖² is constant across candidate types, so `argmin` collapses to
`argmax ⟨θ, c⟩`: an unnormalized inner product in which a centroid's **norm** decides the winner
largely irrespective of its direction. On an untrained panel, a type's centroid norm and the share
of cells it captured were almost perfectly rank-correlated.

> **Prerequisite: the panel must be on the trained feature axis.** A marker that never entered
> the fit is not merely down-weighted, it is *absent* from `{out}.feature_coembedding.parquet`, and
> a type that entered with many markers and scores on one still produces a confident-looking call.
> With a panel only partly trained, spurious assignment to types absent from the tissue was several
> times higher than with a fully trained one.
>
> At gem's default `--n-hvg 0` this is satisfied by construction: every gene is trained and the
> per-gene softmax gate does the selecting. The prerequisite bites only when you set
> `--n-hvg > 0`; then pass `senna gem --markers <panel>` (or `--must-train-features <panel>`) with
> *the same marker file* the annotation will use.
>
> Either way, report `n_live / n_markers` (in `{out}.panel_null.tsv`, or the "marker liveness" log
> line) as a QC figure. If it is not near 100%, nothing downstream is interpretable.

---

## 2. Type prototypes and the per-cell call

For a marker panel assigning genes to types, each type *T* gets an **IDF-weighted centroid**

    e_T = ( Σ_{g ∈ markers(T)} w_g · e_g ) / ( Σ_{g ∈ markers(T)} w_g )

over its **live** markers only. A gene with an all-zero embedding row contributes nothing and is
excluded from *both* numerator and denominator: counting it in the denominator would shrink the
centroid toward the origin in proportion to a type's dead-marker fraction, and a short centroid is
not a weak competitor but a *magnet*. `w_g` is the inverse-document-frequency weight
down-weighting markers shared across many types (`--no-idf` disables).

Cells are assigned by nearest centroid, `t(c) = argmin_T ‖θ_c − e_T‖₂`. Zero-norm centroids are
excluded from the competition: they sit at constant distance `‖θ_c‖` from every cell and would
otherwise capture every cell nearer the origin than to any real prototype.

**QC prune.** Within each type, cells whose distance to their assigned centroid is a high-side
robust outlier (`> median + k·MAD`, `k = --assign-mad`, default **2.5**) are set to `unassigned`:
they took a type by `argmin` but do not actually sit near it (ambient RNA, doublets).

---

## 3. Over-representation within cell groups

A single cell's nearest-centroid call is close to a coin flip; pooling makes it testable. Cells
are grouped by **Leiden** community detection on their cosine kNN graph (`--knn`, default 30;
`--resolution`, default 1.0; see §8 for why a *low* resolution is preferable).

For each (group *K*, type *T*) the count `a = #{c ∈ K : t(c) = T}` is tested against the
hypergeometric null with margins `(N, m_T, n_K)`, where **N and m_T are counted over the cells**
(each cell once), *not* summed from the contingency table. The statistic is `S = −ln P(X ≥ a)`.

`S` is then **calibrated against a permutation null**: the per-cell labels are shuffled with the
group memberships held fixed, and `S` is pooled across groups within a type (the statistic is
relabeling-invariant). `--num-perm` (default 1000) draws; the pool is `n_perm × n_groups`, capped
at 10⁵ per type. The permutation p is Benjamini–Hochberg-adjusted across the types within each
group. A group is called by its top over-represented type if `q < --fdr-alpha` (default **0.1**),
else left uncalled; its cells inherit the call.

**The permutation and the hypergeometric are the same test.** Shuffling the cell labels with the
group memberships held fixed makes the count in (K, T) *exactly* Hypergeometric(N, m_T, n_K)
conditional on the margins (Fisher's argument), which is what the analytic form already computes;
the two agree to numerical precision on every run. The permutation is therefore a
*self-consistency check*, not an independent test.

Calibration diagnostics (`{out}.null_calibration.tsv`): analytic-vs-permutation agreement, and a
genomic-inflation-style λ plus a Kolmogorov–Smirnov statistic on the permutation null itself. λ is
computed with **mid-p** tie-splitting; without it the discrete statistic's ties at the floor pin λ
at a constant, an artifact of the numerical clamp rather than inflation.

---

## 4. Removed: stability bootstrap and set-valued annotation

Earlier builds ran a marker-panel stability bootstrap by default (`--n-boot 200`): each replicate
resampled every type's panel and re-derived the grouping, a call shipped only if its type won at
least `--min-support` of the replicates, and near-ties were reported as a set (`label_set`,
`--set-coverage`, `--max-set-size`). On the negative control of §7 it cut the absent-type share
by an order of magnitude, at the cost of abstaining on a substantial fraction of cells.

It was removed for both costs. **Too stringent:** on sorted progenitor cells annotated with a
mature-cell panel, enrichment's bootstrap abstained on nearly every cell, where the single-pass
FDR call labelled almost all of them. **Too slow:** it took most of the time of an enrichment pass.

What replaces it: each call is the single-pass call (FDR-gated, with Q as the cluster's evidence
split over the types). A contested cluster shows up as a low top share, and is settled by hand in
the annotate TUI (`lupin annotate` without `-o`), where every decision is kept with its rationale in the round's history.

---

## 6. Marker-panel permutation null (bias guard)

*A library option (`panel_perm` in `TermOraConfig`); `lupin annotate` does not expose it yet and runs with it off.*

A single call is blind to bias. The panel null asks the complementary question: *is this answer
better than one a panel that means nothing would have given?*

For each type *T* and each of `P` draws: replace **only** *T*'s panel with `|live(T)|` genes drawn
at random from the pool of *live* marker genes, keeping *T*'s IDF weight multiset, and leave every
rival type's panel real. Rebuild `e_T`, re-run the assignment, and score

    bar[c][T] = min_{S ≠ T} ‖θ_c − e_S‖²          (the rivals: real, and fixed)
    cost(T | panel) = Σ_c min( ‖θ_c − e_T(panel)‖², bar[c][T] )
    p_T = P( cost(T | random genes) ≤ cost(T | T's own genes) )

**The null draw is matched on gene norm**, and this is not a refinement: it decides the answer. A
type's centroid is the mean of its markers' embeddings, so a type whose markers are *long* vectors
gets a long centroid; and because `‖cell‖ ≫ ‖centroid‖` the Euclidean rule degenerates to
`argmax ⟨x, c⟩`, where a longer centroid wins cells almost irrespective of direction (§1). Draw the
null genes *uniformly* and every null panel inherits the pool's mean norm, so a type above that
mean beats its null on norm alone, and one below it loses on norm alone, with no biology tested
either way.

This is GOseq's bias (Young et al., *Genome Biology* 2010) in a different coordinate. GOseq
stratifies on gene *length* because length is the observable proxy for what biases the test
(reads, hence power). Our covariate is the embedding norm itself, which we can measure directly,
so we stratify on it exactly, and skip the noncentral (Wallenius) approximation GOseq needs
because it cannot permute.

The effect is not subtle. Drawn uniformly, the p-values were almost a monotone function of each
type's mean gene norm; stratified on norm, that dependence disappeared. Some types' apparent
significance under the uniform draw was pure norm artifact, and other types were *masked* by it
and became significant once the draw was stratified.

Four further design points, each of which we found to be necessary:

- **One type at a time.** Randomising every panel at once collapses all `C` null centroids onto
  the marker-pool mean; they become mutually indistinguishable and the null fails for reasons
  unrelated to any particular type. Holding the rivals real keeps the competition intact.
- **Same size, so the winner's curse cancels.** A type with few live markers has a high-variance
  centroid, and a noisy prototype wins cells it should not (the maximum of a noisy score is biased
  upward). The null panel is drawn at the same size and is *equally* wobbly, so the advantage
  appears on both sides and divides out. A small panel is asked only whether *these* genes beat
  *any* genes.
- **Null genes are drawn from the live pool.** A random *untrained* gene carries no signal at all,
  so a null of dead genes would be trivially beatable and every type would look significant. This
  holds "is the gene trained?" fixed and isolates "are these the right genes?".
- **The statistic is assignment cost, not cell count.** Occupancy measures whether any rival is
  nearby, not whether the panel is right: on a cleanly separated synthetic panel a random draw
  captures *as many cells as the real one*, because once *T*'s real centroid leaves the competition
  its cells have no near rival and anything in the neighbourhood sweeps them up by elimination.
  Cost separates them; it also has no perverse optimum (a centroid capturing nothing pays the
  maximum `Σ bar`).

On a matched panel the null discriminates: the types that hold cells pass and those that do not,
fail. On a *mismatched* panel it finds nothing significant, which is itself the correct verdict.

---

## 7. Validation design

**Negative control without ground truth.** Use a marker panel that includes types the tissue is
known not to contain (for example mature or terminal types in a progenitor-rich sample). Every
cell they capture is a *countable false positive*, needing no labels. We report this as the
**absent-type share**.

**Headline.** Training the whole panel (rather than part of it) cut the absent-type share several
fold. (The removed bootstrap of §4 cut it by an order of magnitude, by abstaining.)

---

## 8. Negative results (do not re-derive)

**Clustering-free / per-cell-neighbourhood ORA fails.** Replacing the Leiden partition with one
kNN neighbourhood per cell (the Milo device) was **worse than the untested `argmin`** and
manufactured cells for types with almost none in the dataset. Over-representation ranks types by
how **surprising** a count is, not how **likely**: a discovery statistic, not a classifier. The two
rankings coincide only when the group is large; in a small neighbourhood they invert. The pooling
must stay coarse, and Milo's precedent does not transfer because *its* labels (condition, sample)
are external to the embedding whereas ours are `argmin` in the same space that defines the
neighbourhoods. Full write-up: `lupin docs grouping`.

**Parametric Bayesian centroids fail.** A conjugate Gibbs sampler with the centroid as a latent
(`μ_T ~ N(m_T, (ω_T²/κ_T)·I)`) passed every synthetic control and then made the real negative
control *worse* than the baseline: the free per-type variance let one component shrink-wrap a blob
and take a large share of cells at near-zero entropy. Misspecification surfaces as false
confidence, not as uncertainty.

**High Leiden resolution is a bad operating point.** Across replicates, raising `--resolution`
multiplied the number of communities, made the count far less stable between runs, and lowered
the fraction of cells labelled. A low resolution (around 0.5) did best.

---

## 9. Parameters and outputs

| flag | default | § |
|---|---|---|
| `--knn` | 30 | 3 |
| `--resolution` | 1.0 (0.5 recommended) | 3, 8 |
| `--assign-mad` | 2.5 | 2 |
| `--num-perm` | 1000 (pool capped at 10⁵) | 3 |
| `--fdr-alpha` | 0.1 | 3 |
| `--seed` | 42 | · |

| output | contents |
|---|---|
| `{out}.annot.parquet` | per cell: `community`, `coarse_label` with its `coarse_z`/`_p`/`_q`, `fine_label` with its `fine_z`/`_p`/`_q`, and `fine_margin` |
| `{out}.panel_null.tsv` | per type: `n_live`, `occupancy`, `cost`, `null_cost`, `p` |
| `{out}.null_calibration.tsv` | permutation-null diagnostics (λ, KS, analytic agreement) |
| `{out}.cluster_term_{p,q,softq}.parquet` | group × type test matrices |
