# Trajectory bench: the scanpy reference for `lupin trajectory`

Phase 0 of `docs/trajectory-plan.md`. lupin's diffusion pseudotime and type
connectivity replicate scanpy 1.10's `neighbors(method='gauss')`, `diffmap`,
`dpt` and `paga`; this directory produces the reference they are compared
against. The Python scripts are local TDD tools and are not committed: they
produce the golden fixtures the Rust tests check, and this README records how.
Outputs go to `out/` (git-ignored). Nothing is written next to the input data.

## Data

Ainciburu 2022, donor `young4`: 8,911 CD34+ bone-marrow cells (10x v2) with
published BoneMarrowMap labels (`CellType`), from
`~/work/paper-senna/data/BoneMarrowMap/Ainciburu2022_young4.{zarr.zip,cell_metadata.tsv.gz}`.

The 10x 10k BMMNC sample was tried first and dropped: it has about 26 CD34+
cells out of 8,313, so no haematopoietic trajectory to order.

## Commands

From this directory, with the local scripts (`run_scanpy_dpt.py`,
`baseline.py`; uv, dependencies pinned in their headers):

```
D=~/work/paper-senna/data/BoneMarrowMap
mkdir -p out
senna svd $D/Ainciburu2022_young4.zarr.zip --out out/young4_svd30 -t 30

uv run --managed-python --python 3.12 run_scanpy_dpt.py \
  --table out/young4_svd30.latent.parquet --space signed \
  --labels $D/Ainciburu2022_young4.cell_metadata.tsv.gz \
  --root-type HSC --golden 500 --out out/young4_svd30
uv run --managed-python --python 3.12 baseline.py --prefix out/young4_svd30

# refresh the committed fixtures
cp out/young4_svd30.golden.tsv.gz ../../src/tests/fixtures/dpt_young4_svd30_golden.tsv.gz
cp out/young4_svd30.golden_paga.tsv ../../src/tests/fixtures/paga_young4_svd30_golden.tsv
```

`k = 15` counting the cell itself, exact neighbours, 15 diffusion components,
root = the HSC medoid in diffusion distance. The geometry came from senna
0.20.8; the scripts pin their direct dependencies and resolve the rest as of
2026-10-03 (`exclude-newer`).

## Geometry

`senna svd` with 30 components, prepared as lupin prepares a signed latent
(per-column z-score, population standard deviation).

Rejected geometries, both for senna-side reasons:

| run | problem |
|---|---|
| `senna vae` (32-d, 100 epochs) | latent clamped at ±8: 627 cells hit the clamp in every dimension and collapse onto shared corners (many of them HSC) |
| `senna topic` (10 topics) | log θ floored at −17.4 with exact values such as log ¼: saturated proportions, 22 duplicate cells |

With that many coincident cells some cells' neighbours are all at distance 0,
their Gaussian bandwidth is 0, scanpy's kernel turns NaN and ARPACK fails
(`error -9999`).

## Baseline (scanpy 1.10.3, exact neighbours, 2026-10-03)

From `baseline.py`. The expected stage of each staged type (other types are
left out of ρ):

| stage | types |
|---|---|
| 0 | HSC |
| 1 | MPP-MkEry, MPP-MyLy |
| 2 | LMPP, MEP |
| 3 | MLP, Early GMP, BFU-E |
| 4 | CLP, GMP-Neut, GMP-Mono, CFU-E, Megakaryocyte Precursor, EoBasoMast Precursor, Pre-pDC, Pre-cDC |
| 5 | Pro-Erythroblast, Early ProMono, pDC, Pro-B VDJ |
| 6 | Basophilic Erythroblast |
| 7 | Polychromatic Erythroblast |

- **Order**: Spearman ρ = **0.739** between pseudotime from the HSC medoid and
  the expected stage, over 7,224 staged cells; no cell at infinite distance.
  From 30 random HSC instead: median 0.737, range 0.670–0.741.
- **Median pseudotime by type** (≥ 20 cells): HSC 0.048, MPP-MkEry 0.057,
  MPP-MyLy 0.068, MEP 0.080, Pro-B VDJ 0.098, CLP 0.145, MLP-II 0.150, LMPP
  0.159, MLP 0.227, BFU-E 0.245, Megakaryocyte Precursor 0.260, Early GMP
  0.263, Cycling Progenitor 0.302, Pre-pDC 0.303, CFU-E 0.314, Large Pre-B
  0.330, Pre-cDC 0.340, Early ProMono 0.352, EoBasoMast Precursor 0.424,
  GMP-Mono 0.429, GMP-Cycle 0.431, pDC 0.440, Polychromatic Erythroblast
  0.538, Pro-Erythroblast 0.541, GMP-Neut 0.551, Pre-pDC Cycling 0.552,
  Basophilic Erythroblast 0.753.
- **PAGA** (all labels as groups; pairs among types with ≥ 20 cells): 351
  pairs, 290 non-zero, 61 saturated at 1.0. HSC's strongest: MEP, LMPP,
  MPP-MyLy, MPP-MkEry at 1.0, Polychromatic Erythroblast 0.96 (22 cells), MLP
  0.74. Saturation is common, so connectivity alone separates neighbouring
  types poorly.
- **Time**: scanpy's own steps (neighbours, diffusion map, DPT, PAGA) 0.9 s
  warm; the first run in a fresh environment takes about 2.8 s, mostly numba
  compiling.

Earlier exploration, not reproduced by the scripts: with a 10-d SVD the medoid
root gave ρ 0.56, and the first draft's root rule (the HSC farthest from the
terminal types) gave 0.21 — it picks an outlying HSC. scanpy's default
approximate neighbours (used from 8,192 cells) gave ρ 0.742 but an asymmetric
kernel with self-loops and a top eigenvalue of 1.003, which lupin will not
reproduce; hence exact neighbours.

## Fixtures

In `src/tests/fixtures/`, for Rust tests without Python:

Both are inputs to compare against, not data to prepare again: the geometry
was z-scored over all 8,911 cells before the subset was drawn, so a test feeds
it to the kernel as it is; and PAGA grouped every raw label, small ones and
all, so a test groups the same way, without `--min-cells`.

- `dpt_young4_svd30_golden.tsv.gz`: 511 cells stratified by type. A `#`
  comment line (with the subset's 15 eigenvalues, so a test can check the
  eigensolve itself), then a header: `cell`, `label`, the prepared geometry `x0..x29`
  (float32 values written with 9 significant digits), `is_root`
  (`true`/`false`), and scanpy's pseudotime with the subset's own graph.
- `paga_young4_svd30_golden.tsv`: scanpy's PAGA on that subset: `a`, `b`,
  `connectivity`. (PAGA's tree is not kept: with many pairs saturated at 1 its
  edges follow the order of the type names.)

The gzip stream has no timestamp or file name, so regenerating identical content gives
identical bytes.
