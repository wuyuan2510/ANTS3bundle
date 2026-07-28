# Reproducing the DavisLUT results

This bundle reproduces the validation of the `DavisLUT` surface rule against the old
microfacet model on LYSO AFM data (5 / 14 / 28 um etched surfaces). See
`../../../DavisLUT_implementation_notes.md` for the design and the physics discussion.

## Contents

```
prep_afm.py            AFM *.spm.txt  ->  leveled height matrix (um)
pool_luts.py           sum several per-location .lut files into one ensemble LUT
run_reproduce.py       build config + run lsim over depths + make the comparison plots
luts/                  the pooled ensemble LUTs actually used
                         pooled_{5um,14um,28um}_{fwd,rev}.lut
                         REGENERATION.md records inputs, seeds, counters, and validation
config/                config-28Yuan-airGaps-04-LUT.json  (crystal + dual-ended detector setup)
old_reference/         old-microfacet SensorSignals per roughness/depth (for the comparison)
example_heightmap/     one leveled 28um heightmap, as a generateSurfaceLut input example
results/               the reference figures this bundle should reproduce
validation/            actual-rule validation, ideal dual-monitor geometry experiment,
                         and interactive GUI visualization; see validation/VALIDATION_EN.md
```

The raw AFM scans (`LYSO-*.spm.txt`, ~428 MB total) are **not** included — too large for git.
Download them here: https://utexas.box.com/s/3h6rj82d5831qjb9yrdo3rha2mqdwwh6
(the `LYSO-20260312-filtered/` and `LYSO-20260322-filtered/` folders). They are only needed to
regenerate the LUTs from scratch (Level B below); Level A does not need them.

## Prerequisites

Build ANTS3 (qmake, Qt 6.x, ROOT 6.28+); you need `ants3bundle/bin/lsim` (and `bin/ants3` for
LUT generation). Python 3 with numpy, matplotlib (and scipy for the histogram fits).

## Level A — reproduce the simulation results from the LUTs (fast, no AFM data)

```bash
python3 run_reproduce.py 1000 4    # events/depth, parallel workers; use 500 for a quicker check
```
This runs `lsim` for each mirror-periodic roughness over depths [-13,-10,-5,0], then writes to
`results_reproduced/`:
- `pooled_DOI_old_vs_LUT.png` — DOI asymmetry for old microfacet and mirror-periodic LUTs
- `pooled_Rtheta.png` — LYSO->air reflectance R(theta) for the three roughnesses
- `*_hist_fit.png` — DOI-asymmetry histograms and Gaussian fits at each source depth
and prints the DOI table. Compare against `results/`.

## Level B — regenerate the LUTs from AFM data (full pipeline)

First download the raw AFM scans from Box (link above) and unpack the
`LYSO-2026*-filtered/` folders. Then, for one location:
```bash
# 1. level the raw AFM scan (trace = column 1), 10x10 um scan, 512x512 points
python3 prep_afm.py  <LYSO-...-loc1-filtered.spm.txt>  loc1_um.txt

# 1a. diagnose the direct periodic seam and prepare a continuous mirror-tiled input
python3 diagnose_heightmap_seam.py loc1_um.txt --mirror-output loc1_mirror_um.txt

# 2. generate the LUT in the ANTS3 Script window (JavaScript), or headless: ants3 -j gen.js
#    gen.js:
#      rules.generateSurfaceLut("loc1_mirror_um.txt", "loc1.lut",
#          { n1:1.824, n2:1.0, wavelength:420,
#            pixelSizeX:0.0195694, pixelSizeY:0.0195694,   // = 10um / 511
#            format:"matrix", photonsPerBin:40000, alsoReverse:true });
```
Repeat for each location, then pool (ensemble average):
```bash
python3 pool_luts.py  pooled_28um_fwd.lut  loc1.lut loc2.lut ... locN.lut
python3 pool_luts.py  pooled_28um_rev.lut  loc1_reverse.lut loc2_reverse.lut ...
```
Drop the pooled LUTs into `luts/` and run Level A. n1 = crystal index (LSO 1.824), n2 = outer
medium (air 1.0); `alsoReverse:true` writes the swapped-index reverse LUT for the air->crystal
direction.

The generator rejects heightmaps whose opposite edges do not match. This is intentional: its
periodic DDA represents only the height-field triangles and does not add a vertical wall at a tile
boundary. A generic leveled AFM patch must therefore be mirror-tiled (or made C0-periodic by
another explicit construction) before generation.

## Level C — runtime and complete-geometry conformance

```bash
cd validation
../../../bin/ants3 -j validate_rule.js
python3 check_runtime.py
python3 validate_geom.py
python3 validate_rotated_geom.py
```

The geometry commands launch independent `lsim` processes and check R/T, theta marginals, joint
theta-phi distributions, and missing photons. `validate_rotated_geom.py` additionally rotates the
complete setup by nontrivial Euler angles and checks local-to-global covariance against an
unrotated same-seed baseline. See `validation/VALIDATION_EN.md`.

## Notes

- Ensembles: 5 um and 14 um use 3 locations (0312); 28 um uses 8 locations (0312 + 0322),
  matching the `old_reference/` datasets they are compared against.
- The bundled LUTs were regenerated on 2026-07-24 from C0-continuous mirror-tiled AFM maps.
  Every one of the 14 forward/reverse generation pairs reported zero inconsistent escapes and
  zero anomalies. The pooled 5/14/28 um LUTs contain 120000/120000/320000 launched photons per
  incidence bin.
- `results/pooled_Rtheta_3way.png` is retained only as a historical escape-classification
  diagnostic from before the seam repair. The production reference plots are
  `pooled_Rtheta.png` and `pooled_DOI_old_vs_LUT.png`.
- Single-wavelength LUTs (fixed n1/n2), as in Roncali & Cherry 2013.
- `run_reproduce.py` resolves all paths relative to itself, so it works from a fresh checkout
  once ANTS3 is built.
