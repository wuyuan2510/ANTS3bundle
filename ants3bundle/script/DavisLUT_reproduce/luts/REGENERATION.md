# Mirror-periodic LUT regeneration record

The six pooled LUTs in this directory were regenerated on 2026-07-24 after the periodic-seam
diagnosis documented in `../../../../DavisLUT_implementation_notes.md`.

## Inputs and settings

- AFM trace column, 10 x 10 um scan, 512 x 512 measured nodes;
- best-fit plane removed by `../prep_afm.py`;
- x/y mirror tiling by `../diagnose_heightmap_seam.py`, producing 1023 x 1023 nodes with exactly
  matching opposite edges;
- `n1=1.824`, `n2=1.0`, wavelength 420 nm;
- 40 incidence-theta bins, 45 outgoing-theta bins, 36 outgoing-phi bins;
- 40000 photons per incidence bin per location and direction;
- maximum 100 interface bounces;
- three locations for 5 um, three for 14 um, and eight for 28 um.

`wraps` are normal crossings of a continuous periodic boundary. The last column is the sum of
upward-medium-1 and downward-medium-2 inconsistent escapes; it must be zero.

| AFM input tag | Seed | Mean bounces fwd | Mean bounces rev | Wraps fwd / rev | Inconsistent escapes fwd / rev |
|---|---:|---:|---:|---:|---:|
| 5um_0312_loc1 | 2026072401 | 1.106977 | 1.022321 | 648038 / 334138 | 0 / 0 |
| 5um_0312_loc2 | 2026072402 | 1.081973 | 1.019549 | 293760 / 155257 | 0 / 0 |
| 5um_0312_loc3 | 2026072403 | 1.091567 | 1.020868 | 354539 / 210232 | 0 / 0 |
| 14um_0312_loc1 | 2026072404 | 1.261401 | 1.056241 | 817406 / 666146 | 0 / 0 |
| 14um_0312_loc2 | 2026072405 | 1.185781 | 1.037152 | 545330 / 388881 | 0 / 0 |
| 14um_0312_loc3 | 2026072406 | 1.301383 | 1.056613 | 1007181 / 471785 | 0 / 0 |
| 28um_0312_loc1 | 2026072407 | 1.253260 | 1.069969 | 731286 / 565581 | 0 / 0 |
| 28um_0312_loc2 | 2026072408 | 1.296615 | 1.068258 | 953501 / 669200 | 0 / 0 |
| 28um_0312_loc3 | 2026072409 | 1.284333 | 1.070932 | 879387 / 501156 | 0 / 0 |
| 28um_0322_loc1 | 2026072410 | 1.422532 | 1.113426 | 1125961 / 1009306 | 0 / 0 |
| 28um_0322_loc2 | 2026072411 | 1.276682 | 1.053503 | 1136919 / 691537 | 0 / 0 |
| 28um_0322_loc3 | 2026072412 | 1.225558 | 1.047895 | 744922 / 666882 | 0 / 0 |
| 28um_0322_loc4 | 2026072413 | 1.364856 | 1.093235 | 812087 / 594552 | 0 / 0 |
| 28um_0322_loc5 | 2026072414 | 1.244861 | 1.046939 | 713316 / 486653 | 0 / 0 |

Every individual forward and reverse run also reported zero anomalies and zero bounce-limit
discards.

## Pooled outputs and validation

Integer counts and two-dimensional histograms were summed with `../pool_luts.py`. Exact
`Launched = Reflected + Transmitted + Absorbed` conservation and equality between each branch
count and its histogram integral were checked for every incidence bin.

| Roughness | Locations | Launched photons per incidence bin |
|---|---:|---:|
| 5 um | 3 | 120000 |
| 14 um | 3 | 120000 |
| 28 um | 8 | 320000 |

The production-rule validation passed in both directions at 0, 10, 32, 33.3, 34, 45, 80, and
89 degrees. The independent-process geometry validation passed all 16 direction/angle cases with
zero missing photons, including theta and applicable theta-phi TVD checks. The DOI regression
then completed for all three roughnesses at source depths -13, -10, -5, and 0 mm with 1000 events
per depth.
