# Trajectory stack: what was removed and what may be reused

`lupin lineage`, `pseudotime`, `dyn-assoc` and `lineage-plot` were removed
after 0.2.10 (commit `6c1c6b5`) to make room for a redesign. The code is
still in history; restore any file with

    git show 6c1c6b5:<path>        # e.g. src/lineage/assoc/gam.rs

Restored files need dependencies the removal dropped from `Cargo.toml`:
`rand` (`SmallRng` in `orient.rs` and the assoc tests), `statrs`
(distributions in `gam.rs`), and the `mcmc` feature of `legume-numeric`
(`EssSampler` in the Bayesian tests).

## Why it was removed

- **The velocity half never ran.** Edge orientation, the velocity-flux root,
  rewiring by max-weight branching, the velocity-aware layout and arrow grid
  all needed `{from}.velocity.parquet`, and no senna command writes one. In
  practice `lineage` was k-means → MST → user root → Slingshot.
- **The result was the user's prior.** With undirected edges, the root
  (`--root-type`, taken from the marker annotation) fixed pseudotime and branch
  order. Changing it gave nearly uncorrelated pseudotime and non-overlapping
  `dyn-assoc` hits on the same embedding, so the output mostly restated the
  annotation it was rooted on.
- **Redundancy.** Two trajectory fits (`lineage`, `pseudotime`) and two
  figure paths (`lineage-plot`, `plot --colour-by pseudotime`).

## Worth a second look

The statistics in `src/lineage/assoc/` are independent of how the trajectory
is built; they take per-cell pseudotime + branch and a site matrix.

| file | lines | what |
|---|---|---|
| `assoc/gam.rs` | 447 | binomial / quasi-binomial spline GAM of `logit(k/n)` on pseudotime |
| `assoc/bayes_common.rs` | 180 | shared result row and shrinkage prior for the Bayesian tests |
| `assoc/contrast_bayes.rs` | 265 | between-branch contrast at matched pseudotime |
| `assoc/trend_bayes.rs` | 133 | along-branch trend |
| `assoc/contrast.rs`, `assoc/trend.rs` | 135 | frequentist versions (tradeSeq `patternTest` / `associationTest`) |
| `assoc/test_util.rs`, `assoc/*/tests.rs` | ~900 | simulators and tests for the above |

Also reusable as ideas rather than code:

- `orient.rs`: per-edge direction test (bootstrap CI + sign-flip permutation
  + BH), should a velocity table ever exist.
- `traj_annotation.rs`: the marker ORA run on graph nodes instead of Leiden
  communities.
- `docs/lineage-rooting.md` (`git show 6c1c6b5:docs/lineage-rooting.md`): the root-agreement diagnostic and the
  root-invariant split (detection does not depend on the root; trend sign
  and fate polarity do).

The generic numerics (`principal_curve`, `principal_graph`, `branching`,
`hypothesis`) live in `legume-numeric` and are untouched; senna and pinto use
some of them.

## Lessons for the redesign

- Make prior knowledge an explicit input (a root, or a partial order of
  types), and report how far the data supports it instead of hiding it in a
  default.
- Report root-invariant quantities first; direction second.
- Validate against a dataset with a known trajectory before adding machinery.
