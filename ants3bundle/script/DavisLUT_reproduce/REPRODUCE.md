# Reproducing the DavisLUT results

This bundle reproduces the validation of the `DavisLUT` surface rule against the old
microfacet model on LYSO AFM data (5 / 14 / 28 um etched surfaces). See
`../../../DavisLUT_implementation_notes.md` for the design and the physics discussion.

## Contents

```
prep_afm.py            AFM *.spm.txt  ->  leveled height matrix (um)
pool_luts.py           sum several per-location .lut files into one ensemble LUT
run_reproduce.py       build config + run lsim over depths + make the comparison plots
luts/                  the pooled ensemble LUTs actually used (medium-based / paper rule)
                         pooled_{5um,14um,28um}_{fwd,rev}.lut
config/                config-28Yuan-airGaps-04-LUT.json  (crystal + dual-ended detector setup)
old_reference/         old-microfacet SensorSignals per roughness/depth (for the comparison)
example_heightmap/     one leveled 28um heightmap, as a generateSurfaceLut input example
results/               the reference figures this bundle should reproduce
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
python3 run_reproduce.py 1000      # events/depth; use 500 for a quicker check
```
This runs `lsim` for each roughness over depths [-13,-10,-5,0], then writes to
`results_reproduced/`:
- `pooled_DOI_old_vs_LUT.png` — DOI asymmetry vs depth, old (dashed) vs LUT (solid), 5/14/28 um
- `pooled_Rtheta.png` — LYSO->air reflectance R(theta) for the three roughnesses
and prints the DOI table. Compare against `results/`.

DOI histograms with Gaussian fits (as in `results/*_hist_fit.png`) are produced by the project
script `../plot_signal_vs_depth_hist_fits.py` (point its dataset list at the `repro_output/*`
dirs, or at your own output dirs).

## Level B — regenerate the LUTs from AFM data (full pipeline)

First download the raw AFM scans from Box (link above) and unpack the
`LYSO-2026*-filtered/` folders. Then, for one location:
```bash
# 1. level the raw AFM scan (trace = column 1), 10x10 um scan, 512x512 points
python3 prep_afm.py  <LYSO-...-loc1-filtered.spm.txt>  loc1_um.txt

# 2. generate the LUT in the ANTS3 Script window (JavaScript), or headless: ants3 -j gen.js
#    gen.js:
#      rules.generateSurfaceLut("loc1_um.txt", "loc1.lut",
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

## Notes

- Ensembles: 5 um and 14 um use 3 locations (0312); 28 um uses 8 locations (0312 + 0322),
  matching the `old_reference/` datasets they are compared against.
- The LUTs here use the paper-faithful "medium-based" escape handling (see the notes). The
  90-degree spike in the reflected angular distribution is expected and explained there.
- Single-wavelength LUTs (fixed n1/n2), as in Roncali & Cherry 2013.
- `run_reproduce.py` resolves all paths relative to itself, so it works from a fresh checkout
  once ANTS3 is built.
