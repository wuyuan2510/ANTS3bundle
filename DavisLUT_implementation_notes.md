# DavisLUT surface model — implementation record

Adds a look-up-table (LUT) based optical interface rule implementing the approach of
E. Roncali & S. R. Cherry, *Phys. Med. Biol.* **58** (2013) 2185 ("simulation of light
transport based on 3D characterization of crystal surfaces"), plus an offline generator
that builds the LUTs by ray tracing over a measured surface topography (e.g. AFM), plus
GUI and scripting integration.

Rule type string: **`DavisLUT`**, abbreviation **`LUT`**.

Unlike the existing `RoughSurface` rule (a single sampled microfacet), the LUT contains
the full angular reflectance/transmittance of the surface precomputed with multiple
micro-reflections, shadowing and masking, so it stays accurate at large incidence angles.

---

## Design decisions (answering the three original questions)

These were the questions to settle before coding the rule; the answers drive everything below.

**1. How the data is stored (LUT structure).** A custom JSON `.lut`. For each incidence-angle
bin theta_i (default 40 bins over [0, 90 deg)) it stores integer counts
(Launched / Reflected / Transmitted / Absorbed) and two flattened 2-D histograms
(theta_out x phi_out, default 45 x 36; phi_out measured relative to the incidence plane), one
for reflected and one for transmitted rays, plus metadata (n1, n2, wavelength, source heightmap,
grid, pixel size). Integer counts keep it compact (~0.5 MB/LUT), make count-conservation exact,
and let several scans be pooled by summing. A `FormatVersion` field leaves room for a wavelength
axis later; for now it is single-wavelength (fixed n1/n2), as in the paper. The LUT is embedded in
the config JSON so it reaches the `lsim` workers unchanged. (Details: `ALutSurfaceData` below.)

**2. How random angles are generated (depends on 1).** At each surface hit: compute theta_i from
photon.normal; pick the theta_i bin (stochastic interpolation between adjacent bin centers to
avoid a staircase near the critical angle); draw a uniform to decide reflect/transmit from the
stored per-bin R/T; then sample (theta_out, phi_out) from the matching 2-D histogram by
inverse-CDF (CDFs precomputed once before the run, then read-only -> thread-safe). The outgoing
direction is rebuilt in a frame tied to the incidence plane (specular azimuth at phi_rel = 0),
signed so reflected -> back, transmitted -> forward. (Details: `ALutInterfaceRule` below.)

**3. GUI and scripting.** New interface-rule type `DavisLUT` with a GUI editor (load LUT,
R/T-vs-angle and outgoing-distribution plots, n1/n2-mismatch warning; auto-listed in the
rule-type combo and the interface-rule tester), and a new `rules` script unit
(`generateSurfaceLut`, `getLutInfo`, `setLutMaterialRule`, ...). (Details: `ainterfacewidgetfactory`
and `AInterfaceRules_SI` below.)

**Populating the LUT — "construct the data ourselves".** The LUTs are built by a self-contained
offline heightmap ray tracer over the measured AFM height field Z(x,y) (local-normal Fresnel +
Snell/TIR, multi-bounce, continuous mirror tiling). This is essentially the paper's own method, so it
**avoids the ROOT-upgrade / tessellated-volume path entirely** — no infrastructure change was
needed for a working rule + generator + validation. The tessellated-volume route (updated ROOT +
tessellated object shapes + navigator validation) remains a worthwhile future item, but it is not
a prerequisite for this workflow.

---

## New files

### `src/ants3/photonSim/interfaceRules/alutsurfacedata.{h,cpp}` — class `ALutSurfaceData`
Container + (de)serialization + runtime sampling of one surface LUT.
- Schema: metadata (`n1`, `n2`, wavelength, source heightmap, grid, generation stats) +
  per-incidence-angle integer counts (`Launched/ReflectedCounts/TransmittedCounts/AbsorbedCounts`)
  + flattened 2D `(thetaOut x phiOut)` histograms for the reflected and transmitted rays.
  `FormatVersion` gates the reader so a wavelength axis can be added later.
- `writeToJson()` / `readFromJson()` — JSON I/O (custom `.lut` format).
- `check()` — validates binning, array sizes and per-bin count conservation
  (`Refl+Trans+Abs == Launched`).
- `buildRuntime()` — precomputes per-bin reflection/transmission probabilities and the
  inverse-CDF sampling tables (called once before simulation; read-only afterwards).
- `selectThetaBin(thetaDeg, rnd)` — incidence-angle bin, stochastically interpolated
  between adjacent bin centers.
- `sampleOutgoing(reflected, iBin, r1,r2,r3, thetaOut, phiOut)` — samples an outgoing
  direction from the LUT distribution (inverse CDF + within-bin smoothing).

### `src/ants3/photonSim/interfaceRules/alutinterfacerule.{h,cpp}` — class `ALutInterfaceRule : AInterfaceRule`
The runtime interface rule.
- `calculate(photon, normal)` — computes the incidence angle vs the global normal, picks
  the angle bin, decides reflect/transmit/absorb from the LUT, samples `(thetaOut,phiOut)`,
  and rebuilds the outgoing direction in the incidence-plane frame. Returns `Back`/`Forward`.
- `loadLUT(fileName)` — reads a `.lut` file into `Data`.
- `getMaterialConsistencyWarning()` — non-blocking check that the LUT `n1/n2` match the
  refractive indices of the assigned materials.
- `initializeWaveResolved()` / `doCheckOverrideData()` — prepare the runtime CDFs.
- `doWriteToJson()` / `doReadFromJson()` — embeds the full LUT in the config JSON (required
  so it reaches the `lsim` worker processes).
- `canHaveRoughSurface()` = false (the LUT *is* the surface model; surface stays `Polished`).
- `canBeSymmetric()` = false (the LUT bakes in the `n1 -> n2` direction).
- The output direction is explicitly renormalized before returning — see "Bug found" below.

### `src/ants3/photonSim/interfaceRules/alutsurfacegenerator.{h,cpp}` — class `ALutSurfaceGenerator`
Offline generator (built into `ants3` only, not `lsim`).
- `loadHeightmapMatrix(file, dx, dy)` — ASCII matrix of heights (rows<->y, cols<->x).
- `loadHeightmapXYZ(file)` — 3-column x/y/z on a regular grid.
- `generate(result)` — ray traces `photonsPerThetaBin` photons per incidence-angle bin over
  the tiled surface; fills an `ALutSurfaceData`. Config fields: `n1,n2,wavelength,
  thetaIncBins/thetaOutBins/phiOutBins,photonsPerThetaBin,phiSteps,maxBounces,seed,
  reverseGeometry`, plus a `progressCallback` for abort/progress.
- Internals: grid->triangles; 2D DDA traversal (`findIntersection`) with Möller-Trumbore
  (`intersectTriangle`); per-facet unpolarized Fresnel (`fresnelReflection`, incl. TIR) with
  specular reflection or Snell refraction about the *local* normal; multi-bounce loop
  (`tracePhoton`); lateral tiling. Direct periodic repetition is only valid when opposite
  heightmap edges match; production inputs are mirror-tiled to make the boundary C0-continuous.
  Reports mean bounces, anomalies, wrap count and inconsistent-escape diagnostics.

### `src/ants3/script/ScriptInterfaces/ainterfacerules_si.{h,cpp}` — class `AInterfaceRules_SI`
Script unit registered as **`rules`**.
- `generateSurfaceLut(heightmapFile, outLutFile, params)` — runs the generator (params:
  `n1,n2,wavelength,pixelSizeX/Y,format,photonsPerBin,thetaBins,thetaOutBins,phiOutBins,
  phiSteps,maxBounces,seed,comment,reverseGeometry,alsoReverse`); returns generation stats.
- `getLutInfo(lutFile)` — returns a LUT's metadata/binning.
- `validateSurfaceLut(lutFile, params)` — drives the actual
  `ALutInterfaceRule::calculate()` path and checks R/T/A, conditional 1D/2D angular
  distributions, output-vector/status invariants, rotated-normal covariance, and a synthetic
  absorption branch; returns machine-readable measurements and pass/fail limits.
- `setLutMaterialRule(matFrom, matTo, lutFile)` / `setLutVolumeRule(volFrom, volTo, lutFile)`
  — create a `DavisLUT` rule from a `.lut` file and assign it.
- `clearMaterialRule(...)` / `clearVolumeRule(...)`.

---

## Modified files

- `photonSim/interfaceRules/ainterfacerule.h` — added `virtual bool canBeSymmetric()`
  (default true; false for direction-specific rules such as DavisLUT).
- `photonSim/interfaceRules/ainterfacerule.cpp` — registered `"DavisLUT"` in
  `interfaceRuleFactory()` and `getAllInterfaceRuleTypes()`.
- `photonSim/interfaceRules/ainterfacerulehub.h` — added `announceRulesChanged()` (emits the
  existing `rulesLoaded` signal so the GUI refreshes after script-side edits).
- `photonSim/interfaceRules/ainterfacerulehub.cpp` — **bug fix** in `checkAll()`: the volume-
  rule loop reported per-rule errors only when a material-rule error already existed
  (`if(!err..)` -> `if(!es..)`). Pre-existing bug affecting all rule types, exposed while
  wiring DavisLUT as a volume rule.
- `gui/photsim/ainterfacewidgetfactory.{h,cpp}` — new editor widget `ALutInterfaceWidget`
  (Load-LUT, R/T-vs-angle plot, per-bin outgoing distribution, n1/n2-mismatch warning) plus
  its `dynamic_cast` branch in `createEditWidget()`.
- `gui/photsim/ainterfaceruledialog.cpp` — disable/uncheck the "Symmetric" box when the rule
  reports `!canBeSymmetric()`.
- `script/ascripthub.cpp` — register the new `rules` script unit.
- `src/ants3/ants3.pro`, `src/lsim/lsim.pro` — build entries. `ants3.pro` also gained a
  fallback for `python3-config` when the default is too old for `--embed` (was pointing at a
  3.6 build and breaking the link on this machine).

No changes were needed in the photon tracer, the interface-rule tester, or the dispatcher/
farm path.

---

## Workflow: from an AFM surface to a DOI simulation

End-to-end recipe (this is exactly the pipeline used for the 5/14/28 um LYSO validation;
helper scripts live in `ants3bundle/script/LUT_test_results/`).

**1. Preprocess the AFM scan into a heightmap the generator can read.**
A raw AFM export (`*.spm.txt`) is a flat, multi-channel column list, not a height grid, so it
must be converted first. Extract the height column (col 1 = trace, col 4 = retrace), reshape to
the scan grid (e.g. 512x512, row-major), convert to a single length unit, and subtract a best-fit
plane (levels out sample tilt). Write it as either:
- a **matrix** file (rows = y, columns = x, one height per cell) -> use `format:"matrix"` and
  give `pixelSizeX/Y`; pixel size = scan_size / (N-1) in the same units as the heights, or
- an **x y z** three-column file on a regular grid -> use `format:"xyz"`.
Helper: `prep_afm.py` (does exactly this; heights written in um, pixel = 10 um / 511).

The leveled AFM patch is generally not periodic. Before LUT generation, diagnose and mirror-tile
it so that repeated boundaries are height-continuous:

```bash
python3 diagnose_heightmap_seam.py lyso_leveled.txt --mirror-output lyso_mirror.txt
```

**2. Generate the LUT** from the mirror-tiled heightmap, in the ANTS3 Script window:
```js
rules.generateSurfaceLut("lyso_mirror.txt", "lyso.lut",
    { n1:1.824, n2:1.0, wavelength:420,
      pixelSizeX:0.0195694, pixelSizeY:0.0195694,   // um (10um / 511)
      format:"matrix", photonsPerBin:40000, alsoReverse:true });
```
`n1` = index of the medium the photons come from (crystal), `n2` = medium behind the surface
(air/grease). `alsoReverse:true` additionally writes `lyso_reverse.lut` with n1/n2 swapped, for
the return direction (air->crystal). The call returns generation stats (mean bounces, etc.).
Inspect any file later with `rules.getLutInfo("lyso.lut")`.

*Ensemble over several scans (optional but recommended).* The script unit generates one LUT from
one heightmap; to average several locations of the same crystal, call `generateSurfaceLut` once
per location, then **sum their integer count arrays** (`Launched/Reflected/Transmitted/
AbsorbedCounts` and the two histograms) into one pooled `.lut` — equivalent to ray-tracing over
all patches together. This pooling is a small pure-Python post-step: `pool_luts.py`
(`python3 pool_luts.py pooled.lut loc1.lut loc2.lut ...`). It could be folded into
`generateSurfaceLut` later.

**3. Assign the rule** to the crystal<->outside interface, either from a script:
```js
rules.setLutMaterialRule("LYSO", "air", "lyso.lut");          // crystal -> air
rules.setLutMaterialRule("air", "LYSO", "lyso_reverse.lut");  // air -> crystal (reverse LUT)
```
(or `setLutVolumeRule` for a volume-name pair), **or** from the GUI (next section). The LUT data
are embedded into the config, so the `.lut` file is not needed at run time.

**4. Run the photon simulation** as usual — `lsim.simulate()` from the Script window, or the
headless worker `lsim <workdir> <config.json> <id>` (when running the worker directly, set
`PhotonSim.Run.EventFrom=0`, `EventTo=<Flood.Number>`, and a non-zero `Seed`, which the GUI
otherwise fills in per worker).

**5. Analyze** `SensorSignals.txt` (per-event sensor signals). For the dual-ended DOI studies
here the observable is `(s0-s1)/(s0+s1)` vs source depth. The reproducibility bundle's
`run_reproduce.py` plots it versus depth and also produces per-depth histograms with Gaussian
fits.

---

## Using the DavisLUT rule in the GUI

**Generating a LUT** is script-only (no GUI dialog yet): open the **Script** window and call
`rules.generateSurfaceLut(...)` as above. Everything else can be done in the GUI.

**Assigning / editing the rule:**
1. Open the **interface-rule window** (photon-simulation GUI). It shows a material x material
   matrix and a volume-pair list.
2. **Double-click the cell** for the interface you want (e.g. row = LYSO, column = air).
   The interface-rule dialog opens.
3. In the **rule-type combo box**, choose **`DavisLUT`**. The DavisLUT editor panel appears.
4. Click **Load LUT** and pick the `.lut` file. The info line shows `n1 -> n2`, wavelength,
   binning and mean bounces. A **red warning** appears if the LUT's `n1/n2` disagree (>1%) with
   the refractive indices of the two materials this cell connects — a guard against assigning a
   LUT to the wrong material pair or direction.
5. Inspect the LUT: **R/T vs angle** plots the reflection and transmission probability against
   incidence angle; enter an incidence angle, pick **Reflected/Transmitted**, and
   **Show angular distribution** to view the 2D (theta_out, phi_out) outgoing map for that angle.
6. The **Symmetric** checkbox is **disabled** for DavisLUT — the LUT bakes in the n1->n2
   direction, so assign the reverse-direction LUT to the opposite cell (air -> LYSO) separately,
   using a LUT generated with swapped indices (`alsoReverse:true`).
7. **Accept**. The rule (with its embedded LUT) is now in the configuration; save the config as
   usual. The rule also appears in the material matrix with the abbreviation **`LUT`**.

**Checking it before a run:** the interface-rule **tester** (in the same GUI) shoots test photons
at the selected rule and draws the resulting reflected/transmitted directions — a quick way to
confirm the LUT behaves as expected before launching a full simulation.

**Running:** build the geometry and run the photon simulation as normal; the DavisLUT rule is
applied automatically at the assigned interface, in both the GUI and headless/farm runs.

---

## Bug found and fixed while testing

`ALutInterfaceRule::calculate()` returns the photon direction directly (the tracer takes it
as-is for non-rough rules, without renormalizing). Even though the direction is unit by
construction, a floating-point `|v| = 1+epsilon` made a downstream
`acos(sensorNormal . v)` in `ASensorModel::getAngularFactor()` return NaN at a near-normal
sensor hit, which then indexed the angular-response array out of bounds and crashed `lsim`.
Fix: explicitly renormalize the output direction at the end of `calculate()`. (The old
`RoughSurface` rule avoided this because its Fresnel path in the tracer renormalizes.)

---

## What was tested

Verification used small standalone harnesses linked against the compiled object files, plus a
full end-to-end run in `lsim`.

1. **Flat-plane analytic check** — a zero heightmap, n1=1.824, n2=1.0:
   - generated `R(theta)` matches the analytic unpolarized Fresnel curve to within 0.003,
     with the total-internal-reflection step at the 33.3 deg critical angle;
   - reflected rays land at `thetaOut = thetaInc`, transmitted at the Snell angle, both in the
     incidence plane (`phiRel ~ 0`) — confirms the generator<->rule frame convention;
   - per-bin count conservation holds; `buildRuntime/selectThetaBin/sampleOutgoing` round-trip.

2. **Rough (egg-carton) surface** — multi-bounce (>1 avg), strongly broadened angular
   distributions, count conservation, and JSON write/read round-trip all pass. Anomaly rate
   scales with surface steepness (grazing/multi-bounce rays that cannot be cleanly classified
   are excluded from the LUT and reported in its metadata).

3. **End-to-end on real data** — `script/LYSO-28um-...-loc1.spm.txt` (AFM, plane-leveled) ->
   forward+reverse LUTs -> config with LYSO<->air rules swapped to DavisLUT -> `lsim` depth
   scan vs the old `RoughSurface` rule. The LUT method produced a clean, monotonic DOI
   response with a stronger depth slope and lower/more depth-dependent light collection,
   consistent with the paper. Artifacts and a comparison plot are in
   `ants3bundle/script/LUT_test_results/` (see its README for the pipeline and caveats).

4. **Runtime-rule and ideal-geometry conformance** — the validation bundle now calls the real
   `ALutInterfaceRule::calculate()` for both directions at 0, 10, 32, 33.3, 34, 45, 80 and
   89 degrees, including full `(theta_out, phi_out)` TVD checks, status/direction invariants,
   rotated normals and synthetic absorption. A separate dual-monitor experiment then exercises
   the complete geometry -> tracer -> interface -> monitor path with bulk losses disabled and
   measures absolute R/T plus the conditional theta and joint `(theta,phi)` distributions.
   A rigidly rotated complete-geometry regression uses Euler angles `(37,29,23) deg` at 45-degree
   incidence in both directions. Its same-seed baseline and rotated monitor histograms are
   identical bin by bin, with zero missing photons, directly validating the LUT rule's
   local-to-global outgoing-direction conversion.
   Each geometry case is run as a fresh static `lsim` process; the earlier multi-run GUI
   dispatcher loop was found to reuse source state across nominally different angles and has been
   deprecated. Both suites have independent process exit codes. The loadable
   `ants3bundle/script/DavisLUT_reproduce/validation/validate_geom_gui.json` configuration together with the
   `validate_geom_gui.txt` GUI script provides the same experiment interactively, with the
   setup/tracks in the Geometry window, reflected/transmitted theta overlays, and direct monitor
   `(theta_out, phi_out)` histograms compared with the LUT using TVD. The photon monitor now
   serializes an `AnglePhi` histogram and exposes it as `lsim.getMonitorAnglePhi()`.
   `validate_rotated_geom_gui.txt` rebuilds that setup with editable Euler angles and provides
   the corresponding interactive local-to-global conversion test.

Both `ants3` and `lsim` build cleanly (qmake, Qt 6.5.3, ROOT 6.28/04 on this machine).

---

## Generator periodic-seam diagnosis (2026-07-24)

Earlier generator runs produced a sizeable population whose direction disagreed with its medium:
typically a photon was still in medium 1 after reflection but `dir.z > 0`. Medium-based fallback
classification conserved it as reflected, and clamping the incompatible reflected polar angle
put it into the final `theta_out = 90 deg` bin. Earlier revisions of this document described that
population as a normal height-field multi-bounce escape. That explanation was wrong.

For a continuous height field `z = Z(x,y)`, a ray that starts on the lower side and eventually
rises above `Zmax` must intersect the surface again: the continuous signed height
`g(t) = z_ray(t)-Z(x_ray(t),y_ray(t))` changes from negative to positive. Therefore
"medium 1 + upward + no next intersection" proves that the numerical surface is open or an
intersection was missed; lack of overhangs cannot explain it.

The actual cause is the **discontinuous periodic tile seam**. The generator repeats AFM indices
periodically, but opposite edges of a generic leveled AFM patch do not have matching heights and
no vertical wall joins them. For the representative 28 um location, ordinary neighbor-step RMS
is only about `0.010--0.013 um`, while the x/y periodic edge-jump RMS is
`0.616/0.697 um`. Across all eight 28 um scans used in the pooled LUT, edge jumps are about
`20--132` times their internal neighbor-step RMS.

Diagnostics added to `ALutSurfaceGenerator` record whether every inconsistent escape crossed a
periodic boundary anywhere in its history. A same-seed A/B test on the representative scan gave:

| test | photons | UpEscape | DownEscape | reflected last-theta-bin / R |
|---|---:|---:|---:|---:|
| original discontinuous periodic tiling, 40 incidence bins | 40,000 | 2,783 | 154 | 9.82% |
| continuous x/y mirror tiling, same settings | 40,000 | 0 | 0 | 0.17% |
| original tiling, 45 deg only | 100,000 | 8,509 | 486 | 11.19% |
| continuous mirror tiling, 45 deg only | 100,000 | 0 | 0 | 0.17% |

All `8,509/8,509` upward and `486/486` downward inconsistent escapes in the 45-degree original
run had crossed a seam. The continuous mirror case still recorded `32,852` normal periodic-wrap
events, so the cure is continuity, not avoiding boundary crossings. A previously tested
watertight triangle-intersection replacement had no effect because ordinary shared triangle
edges were not the opening.

**Impact:** the `theta_out = 90 deg` reflected spike is a seam artifact, not expected surface
physics. The old pooled LUTs and plots generated with discontinuous tiling are legacy diagnostics.
Medium-based classification remains a useful invariant and count-conserving fallback, but it
must not be used to hide a non-zero inconsistent-escape count. For a production LUT, both
`upEscapeReclassified` and `downEscapeReclassified` should be zero (apart from an explicitly
justified numerical tolerance).

The preferred repair is continuous mirror tiling (or another explicitly C0-periodic
preprocessing), followed by regeneration of every per-location forward/reverse LUT, pooling, and
rerunning the rule/geometry/DOI validation. Simply spreading the last-bin spike or changing it to
direction-based classification would mask the open seam without repairing the ray geometry.

That repair was completed on 2026-07-24:

- 14 AFM locations were plane-leveled and mirror-tiled to `1023 x 1023` nodes; all x/y seam
  differences were exactly zero;
- 14 forward/reverse pairs were generated with 40 incidence bins and 40000 photons per bin;
  every pair reported zero anomalies, zero bounce-limit discards, and zero inconsistent escapes;
- pooled 5/14/28 um LUTs contain 120000/120000/320000 photons per incidence bin with exact
  count and histogram conservation;
- direct runtime validation passed in both directions; the isolated-process geometry validation
  passed all 16 cases with zero missing photons; the rotated complete-geometry covariance test
  passed in both directions; the 5/14/28 um, four-depth DOI rerun completed with 1000 events per
  depth;
- the reflected last-theta-bin fraction near 46 deg fell from about
  `4.0/9.9/12.5%` for the old 5/14/28 um forward LUTs to
  `0.15/0.14/0.16%` after the repair.

The generator now rejects a heightmap at load time when opposite x/y edges do not match within
numerical tolerance, preventing accidental regeneration of the open-seam artifact.

---

## Known limitations (v1)

- Single-wavelength LUT (applied to all photons regardless of `waveIndex`); a wavelength axis
  is the intended `FormatVersion 2` extension.
- LUT generation is single-threaded and exposed via script only (no GUI "generate" dialog yet).
- Periodic DDA requires matching opposite heightmap edges and does not synthesize vertical seam
  walls. Non-matching inputs are rejected; use continuous mirror tiling (or another C0-periodic
  construction) and require zero inconsistent-escape counters.
