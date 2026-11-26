# Confirming a round by data fission

`lupin confirm -f ROUND` checks a round's calls on half of the embedding its
clusters were not made from, and with half of each type's markers held out.
It adds a line of evidence beside the round's others and changes none of
them; `lupin why` weighs it as its own line, **fission**.

## The double dip it removes

A round clusters the cells on their embedding `z` and then tests the panel
on the same `z`. Clustering finds the directions the cells split along, so
any marker set aligned with them looks up in the clusters, even where the
cells are one population split by noise. On a homogeneous population (no
cluster, no type), clustered and tested on one `z`, a type's markers were up
at p < 0.05 in 44% of (cluster, type) pairs (the test
`fission_keeps_false_calls_at_their_level_where_the_double_dip_does_not`).

## Fission

With noise `w ~ N(0, Σ)` drawn per cell and `τ > 0`:

```text
f = z + τw        (clustered: Leiden, from legume-numeric, as the run's own)
g = z − w/τ       (tested)
```

When `z_i ~ N(μ_i, Σ)`, `f` and `g` are independent (Leiner, Duan, Wasserman
& Ramdas, *Data fission: splitting a single data point*, JASA 2025,
arXiv:2112.11079), so the clusters of `f` are fixed in advance as far as `g`
is concerned. The embedding is not Gaussian. The split is independent only
as far as it is asymptotically normal, a central-limit argument and no more.

`τ = 1` splits the information evenly. A smaller `τ` keeps more for the
clusters and less for the test, which suits clusters far apart.

### Σ is fitted, not assumed

`Σ` must be the cells' noise, and `z` alone does not tell the noise from the
cells' spread of means. A wrong `Σ` leaves `Cov(f, g) = Σ_noise − Σ`:

- too small (the covariance within the round's clusters, where they
  over-split): `f` and `g` split together, and a spurious cluster confirms
  itself;
- too large (the cells' total covariance, where types differ): a type split
  in two by the noise has `g` shifted the other way in each half, and false
  calls follow. A planted run showed this: halves of type B called C and D.

So `Σ` is fitted to the condition that holds when it is right: `f` and `g`
do not co-vary within any cluster of `f`. The slope of `g` on `f` within the
clusters, `M`, estimates `D V⁻¹`, with `D = Σ_noise − Σ` and
`V = D + (1 + τ²) Σ`. The slope is used rather than the covariance because
selecting cells on `f` cuts `f`'s spread short, and the covariance with it,
but not the slope. Each round adds `D = (I − M)⁻¹ (1 + τ²) M Σ` to `Σ` and
splits again, until the cross-covariance is within 1.5 times its size by
chance, or after five rounds. The fit starts from the round's within-cluster
covariance.

`dependence` reports what is left: the cross-covariance in units of its
chance size. It is about 1 for independent halves, well above 1 when they
still depend.

`--sigma total` and `--sigma within` skip the fit. Each is right only in its
own case: total when nothing differs, within when the round's clusters are
the types.

## The ranking

A call needs only each cluster's ranking of the genes. Cluster `k`'s ranking
comes from `g`, through the decoder's gene table `ρ`:

```text
t_gk = ρ_g · (ḡ_k − m) / √(ρ_gᵀ Σ_g ρ_g / n_k)
```

Here `ḡ_k` is the cluster's mean of `g`, `m` the mean of the cluster means,
and `Σ_g` the covariance of `g` within the clusters. So `t_gk` is the
cluster's mean logit for the gene against the balanced centre, in standard
errors. Its sufficient statistics are the clusters' sizes, their means of
`g`, and one pooled covariance.

## Two tests, both needed

1. **Against random gene sets**: the enrichment pass's engine on the ranking
   (weighted KS, multilevel tail, BH within the cluster; no specificity
   transform, since the ranking is already a contrast). The ranking projects
   every gene's `ρ_g` on one direction, so markers with similar rows move
   together, which random sets do not. Alone, this test is
   anti-conservative.
2. **The markers up** (`up p`): the markers' weighted mean of `t`, which is
   `u_T · (ḡ_k − m)` with `u_T = Σ w_g ρ_g / sd_g / Σ w_g`, against its exact
   variance under Gaussian `g`, `u_Tᵀ Σ_g u_T · v_k`. The factor
   `v_k = (1 − 1/K)²/n_k + Σ_{j≠k} 1/(K² n_j)` is the variance of a cluster
   mean less the mean of all the cluster means. The markers' correlation is
   in `Σ_g`.

A type stands in a cluster when its q (1) and its up p (2) are both below
0.05. The call is the type with the largest share of the cluster's Q, as the
enrichment pass calls.

## Held-out markers

Each type's markers are split in two halves by a hash of the gene name. One
half picks the cluster's call on the ranking from `f`; the other half tests
that call on the ranking from `g`, both ways round. The pick is made on `f`
because a call picked on `g` with one half would carry over to the other
half through the markers' correlation: the double dip again, across genes.

## Calibration

| | markers up at p < 0.05 | stands on both tests |
|---|---|---|
| clustered and tested on one `z`, 1500 cells | 44% | 18% |
| fission, `Σ` fitted, 1500 cells | 6.1% | 5.5% |
| fission, exact `Σ`, 600 cells | 3.8% | 3.7% |

These are homogeneous populations in 6 dimensions, 20 seeds each, so every
call is false and the nominal rate is 5%. With `Σ` fitted, the rate sits a
little above nominal (about 7% at 600 cells, 6% at 1500) and falls toward it
as cells are added, as the fit's error shrinks; with `Σ` exact it is at or
below nominal, so the excess is the fit's. A planted run (five clusters, four
with a type in the panel) is called right in every typed cluster, with every
held-out half holding. Its typeless cluster gets no call, though the round's
score pass, an argmax with no abstention, had labelled it.

## What it does not remove

- **The counts are not split.** They made the embedding and the gene table,
  so the round's enrichment and SuSiE, which read the counts, stay
  post-selection. Fission removes the double dip of clustering and testing
  on one embedding, given that embedding.
- **The round still enters twice.** Its within-cluster covariance starts the
  fit (which the fit's condition then corrects). Its clusters are matched to
  those of `f` by shared cells, for the report and for the **fission** line
  of `lupin why`, and a call is said to agree with its label. These describe
  and change no statistic.
- **The genes are the decoder's.** The ranking reads the logits the model
  implies, not observed expression.

## Output

Beside the round: `{round}.confirm.json` (what `lupin why` reads),
`{round}.confirm.tsv` (one row per cluster of `f`) and
`{round}.confirm.clusters.parquet` (the clusters of `f`).
