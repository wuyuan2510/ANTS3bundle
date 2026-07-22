# DavisLUT validation

This directory separates three questions that are easy to conflate:

1. Does the production runtime rule reproduce the R/T/A probabilities and outgoing-angle
   distributions stored in the LUT?
2. Does the same rule remain correct when exercised through the complete
   **geometry -> photon tracer -> interface -> monitor** chain?
3. Can the experimental setup, reflected/transmitted directions, and measured distributions be
   inspected interactively in the ANTS3 GUI?

These checks validate the software path that consumes a LUT. Independent validation of the LUT
generator physics—for example, Fresnel/Snell behavior on a flat height map, the total-internal-
reflection threshold, or comparisons with another ray tracer or measurements—is a separate layer.
See `../../../../DavisLUT_implementation_notes.md` for that work.

## 1. Production-rule validation: `validate_rule.js`

This is the primary numerical validation. The script calls `rules.validateSurfaceLut()`, which
drives the production `ALutInterfaceRule::calculate()` implementation directly. The test does not
copy the outgoing-direction construction algorithm into a second implementation and then use that
copy as its oracle.

For both `pooled_28um_fwd.lut` and `pooled_28um_rev.lut`, the test covers:

- incidence angles `0, 10, 32, 33.3, 34, 45, 80, 89 deg`, including normal incidence, the
  LYSO-to-air critical-angle region, and grazing incidence;
- measured versus stored R/T/A probabilities;
- the conditional reflected and transmitted `theta_out` marginals;
- the complete conditional `(theta_out, phi_out)` distributions;
- interpolation between adjacent incidence-angle bins, using

  ```text
  H(out | reflected) = [(1-f) R0 H0 + f R1 H1] / [(1-f) R0 + f R1]
  ```

- unit length of every returned direction vector, the outgoing hemisphere, and the returned rule
  status;
- rotational covariance with a non-axis-aligned interface normal and incidence plane;
- an otherwise uncommon absorption path using a small synthetic LUT with
  `R/T/A = 0.2/0.3/0.5`.

Monte Carlo acceptance limits are derived from each expected multinomial distribution. A failure
sets `passed:false` in the result and makes `validate_rule.js` exit with a non-zero status.

Run from this directory:

```bash
../../../bin/ants3 -j validate_rule.js
python3 check_runtime.py
python3 plot_runtime_validation.py
```

Principal outputs:

- `data/runtime/{fwd,rev}.json`: numerical measurements and limits; when `include2D:true`, it also
  contains flattened two-dimensional distributions;
- `figures/runtime_probabilities.png`: expected and sampled R/T/A;
- `figures/runtime_theta_{fwd,rev}.png`: reflected and transmitted theta marginals;
- `figures/runtime_2d_tvd.png`: TVD of each complete two-dimensional distribution together with
  its statistical acceptance limit.

### Low-level sampler diagnostic

`validate_lut.cpp` remains as a small diagnostic for `ALutSurfaceData`. It invokes
`selectThetaBin()` and `sampleOutgoing()` directly and is useful for isolating CDF construction or
sampling problems, but it is not a rule-level validation. Its conditional interpolation oracle
also includes the reflection probability of each neighboring incidence bin.

## 2. Automated dual-monitor geometry validation: `validate_geom.js`

The geometry is deliberately ideal and compact:

- air is above `z=0`, and LYSO is below it;
- a collimated pencil beam starts only `0.002 mm` from the interface;
- two photon monitors are located at `z=+/-0.004 mm` and stop photons when hit;
- one monitor receives reflected photons and the other receives transmitted photons;
- material names and refractive indices come from the real DOI configuration, but bulk absorption
  and scattering in both LYSO and air are disabled, and the wavelength is fixed at `420 nm`;
- placing the monitors close to the interface prevents the loss of grazing photons through the
  sides, which biased R/T in the older geometry validation.

The forward and reverse LUTs are tested at the same eight angles used by the production-rule test.
The experiment therefore measures both absolute R/T and the conditional theta distribution of
each branch.

```bash
../../../bin/ants3 -j validate_geom.js
python3 check_geometry.py
python3 plot_geom.py
```

`check_geometry.py` exits non-zero if the probability error, theta TVD, or uncollected fraction
exceeds its corresponding limit. Numerical results are written to `data/geometry/`; plots are
written as `validation_geom_probabilities.png` and `validation_geom_{fwd,rev}.png`.

## 3. GUI geometry experiment: JSON configuration and TXT script

The GUI workflow is split into two files with distinct responsibilities:

| File | Where to load it | Responsibility |
|---|---|---|
| `validate_geom_gui.json` | Main window: **Configuration -> Load** | Materials, geometry, two monitors, both directional DavisLUT rules, ideal transport settings, and the default photon source |
| `validate_geom_gui.txt` | **Scripting -> JavaScript -> Load script** | Select direction and angle, run the simulation, read monitors, calculate expected distributions, show the geometry, and draw comparisons |

The `.txt` file still contains ANTS3 JavaScript. The extension is used because it is accepted by
the GUI script picker. The JSON file is a complete, independently loadable configuration; the
script no longer creates the geometry or injects interface rules at runtime.

### 3.1 GUI procedure

1. Start ANTS3.
2. In the main window select **Configuration -> Load** and open `validate_geom_gui.json`.
3. Optionally select **View -> Geometry** to inspect the static setup; the tested interface is at
   `z=0`.
4. Select **Scripting -> JavaScript**, click **Load script**, and open
   `validate_geom_gui.txt`.
5. Adjust the five controls at the beginning of the script if needed, then click **Run**.

```js
var VALIDATION_DIR = "";
var DIRECTION = "fwd"; // or "rev"
var ANGLE = 45.0;
var PHOTONS = 300000;
var SEED = 20260714;
```

If ANTS3 was started from this `validation/` directory, leave `VALIDATION_DIR` empty. Otherwise,
set it to the absolute path of this directory. LUT oracle files and output paths are resolved from
that location. Before running, the script checks the configuration name, the required GUI APIs,
direction and angle values, photon count, LUT file, and monitor count, so a misloaded setup fails
with an actionable message.

### 3.2 Geometry stored in `validate_geom_gui.json`

The JSON retains only the two materials used by this experiment: `air` (material 0) and `LYSO`
(material 1). The DOI configuration's historical `airbubble`, `epoxy`, and `grease` materials and
their interface rules are intentionally absent. The Interface Rules window should therefore show
only the two relevant DavisLUT rules.

| Object | Parent / material | Size and position | Purpose |
|---|---|---|---|
| `AirBox` | `World` / air | `30 x 30 x 20 mm`, centered at `z=0` | Provides air above the interface and contains the LYSO block |
| `LYSO` | `AirBox` / LYSO | `20 x 20 x 10 mm`, centered at `z=-5 mm` | Its top face is exactly at `z=0`, forming the only tested interface |
| `LowerMonitor...` | `LYSO` / LYSO | `12 x 12 mm`, global center at `z=-0.004 mm` | Receives photons leaving the interface downward |
| `UpperMonitor...` | `AirBox` / air | `12 x 12 mm`, global center at `z=+0.004 mm` | Receives photons leaving the interface upward |

Both monitors use 45 bins over `0..90 deg` for theta and 36 bins over `0..360 deg` for phi,
exactly matching the `45 x 36` outgoing-angle grid of the 28 um LUTs. Both stop tracking on a hit.
Their distance from the interface is only `4 um`: enough to distinguish reflection from
transmission, but short enough that even near-grazing photons rarely escape around a monitor.

The lower monitor has local `z=4.996 mm` inside the LYSO object. Since the LYSO center is at global
`z=-5 mm`, its global monitor position is `-5+4.996=-0.004 mm`. This explains the apparently
asymmetric z coordinates shown in the JSON and Geometry editor.

The JSON also stores:

- a fixed wavelength of `420 nm`;
- disabled bulk absorption and Rayleigh scattering in air and LYSO, so only the interface rule
  controls the result;
- a default single-source run with `fwd`, `45 deg`, and `300000 photons`;
- `pooled_28um_fwd.lut` embedded in the `LYSO -> air` rule and
  `pooled_28um_rev.lut` embedded in the `air -> LYSO` rule;
- enabled photon-monitor output and disabled storage of a large set of real photon tracks.

Consequently, after loading only the JSON, the GUI can already show the materials, geometry,
source, and two rules. The TXT script reads the external LUT files only to construct an independent
"stored LUT" plotting oracle. Photon tracking uses the DavisLUT data already embedded in the JSON.

### 3.3 Script walkthrough

#### A. Controls and preflight checks

`DIRECTION="fwd"` launches photons from LYSO toward air; `"rev"` launches them from air toward
LYSO. `ANGLE` is the incidence angle relative to the interface normal and must be in `[0,90)`.
A fixed `SEED` makes a run reproducible, while `PHOTONS` controls statistical noise and runtime.

`EXPECTED_CONFIG_NAME` verifies that the GUI validation JSON is loaded.
`lsim.countMonitors()==2` then verifies that the geometry hub contains both photon monitors. These
checks prevent the script from accidentally running against a normal DOI configuration and
producing plausible-looking but semantically incorrect plots.

#### B. Source position and direction

The source is placed on the incident side, `depth=0.002 mm` from the interface. With
`a = ANGLE*pi/180`:

```text
fwd: source z = -depth, direction = (sin(a), 0, +cos(a))
rev: source z = +depth, direction = (sin(a), 0, -cos(a))
source x = -depth*tan(a)
```

The x displacement while propagating to `z=0` is exactly `depth*tan(a)`, so both directions strike
the interface at `(0,0,0)`. The source lies between the interface and the incident-side monitor;
therefore the incoming beam is not stopped before reaching the interface. Only photons processed
by the interface can reach either monitor.

This section uses `config.replace()` only for per-run values: photon count, output directory,
source position, and direction. `config.updateConfig()` then synchronizes these JSON changes with
the runtime hubs. The script does not recreate the geometry or reassign rules.

#### C. Three guide lines in the Geometry window

Before simulation, `geowin.redraw()` displays the actual configured geometry and adds a red
incident line. Afterward, the script adds blue and green lines using the measured mean theta from
the two monitor distributions:

- red: incident pencil beam;
- blue: representative direction for the measured reflected conditional distribution;
- green: representative direction for the measured transmitted conditional distribution.

These are representative guide tracks, not individual Monte Carlo tracks loaded from the output.
This avoids storing hundreds of thousands of tracks and avoids presenting one random photon as if
it represented the complete distribution.

#### D. Simulation and monitor role reversal

`lsim.simulate()` runs the production photon tracer. The script then calls
`lsim.loadMonitorData("geom_gui_out/sim/PhotonMonitors.txt")` to load the results into the monitor
API.

Geometry creation order fixes monitor 0 below and monitor 1 above the interface. Their physical
roles change with launch direction:

| Launch direction | Reflection monitor | Transmission monitor |
|---|---:|---:|
| fwd, LYSO -> air | 0 (lower) | 1 (upper) |
| rev, air -> LYSO | 1 (upper) | 0 (lower) |

The script therefore uses:

```js
var reflectedIndex   = reverse ? 1 : 0;
var transmittedIndex = reverse ? 0 : 1;
```

Absolute probabilities are monitor hits divided by `PHOTONS`. The console value
`A-or-uncollected = 1-R-T` is deliberately named as a residual. In this compact ideal geometry it
should be close to interface absorption A, but strictly it also contains any photon missed by both
monitors.

#### E. Monitor normalization and angular coordinates

`lsim.getMonitorAngle(index)` returns `[binCenter, count]`. `normalizeMonitor()` calculates

```text
p(theta_i | branch) = count_i / (sum_j count_j * binWidth)
```

The y-axis unit is therefore `1/deg`, and the theta integral is one. This is the angular shape
conditional on reflection or transmission; it does not include the branch's absolute probability,
which is reported separately from monitor hits.

`meanAngle()` computes a conditional weighted mean only for the blue and green Geometry-window
guide lines. It is not part of the numerical pass/fail decision.

When a photon enters a monitor, the tracer transforms its direction into the monitor's local frame
and calculates

```text
phi_monitor = atan2(v_local_y, v_local_x), wrapped to [0,360) deg.
```

`AnglePhi` is filled once per photon as a genuine two-dimensional ROOT histogram. It is not
reconstructed from two independent one-dimensional marginals. It is serialized in
`PhotonMonitors.txt` as:

```text
AnglePhi = {
  Xbins:45, Xfrom:0, Xto:90,
  Ybins:36, Yfrom:0, Yto:360,
  Data:[[thetaCenter, phiCenter, count], ...],
  Entries:numberOfMonitorHits
}
```

`lsim.getMonitorAnglePhi(index)` returns only normal bins, also as
`[thetaCenter, phiCenter, count]`. The JSON `Data` array includes ROOT underflow and overflow bins,
so it has `(45+2)*(36+2)=1786` entries; the script API returns the `45*36=1620` normal cells.

DavisLUT phi is not an arbitrary global azimuth. Its basis is defined by `ex`, the tangential
component of the incoming direction, and `ey=N x ex`; phi increases from `ex` toward `ey`. In this
experiment the incidence plane is global xz and the monitors are not rotated:

- for `fwd`, `N=+z` and `ey=+y`, so monitor-local phi can be compared directly with LUT phi;
- for `rev`, `N=-z` and `ey=-y`, so the script applies
  `phi_LUT=(360-phi_monitor) mod 360`, equivalently reversing the phi bins.

At exact normal incidence, the incidence plane has no unique direction. The two-dimensional phi
TVD is therefore reported as N/A, while R/T and the theta marginals remain valid checks.

#### F. Independent expectation calculated from the LUT

LUT incidence-angle bin centers are

```text
theta_k = (k+0.5) * 90/N_inc.
```

`incidenceMix()` calculates `x=ANGLE/(90/N_inc)-0.5`, lower index `k0=floor(x)`, upper index
`k1=k0+1`, and interpolation weight `f=x-k0`; values at either endpoint are clamped to the nearest
bin.

`expectedProbability()` obtains each branch probability from the launched and
reflected/transmitted counts in each incidence bin, then interpolates:

```text
P_branch = (1-f) * C(k0)/L(k0) + f * C(k1)/L(k1).
```

`lutMarginal()` sums the LUT's two-dimensional `(theta_out,phi_out)` histogram over phi. The
important detail is that it first interpolates unconditional probabilities using launched counts
and only then conditions on the interpolated branch probability:

```text
H(theta_i | branch) =
  sum_phi[(1-f) N(k0,i,phi)/L(k0) + f N(k1,i,phi)/L(k1)]
  / (P_branch * delta_theta_out).
```

It would be incorrect to normalize the two neighboring bins separately and blend their
conditional distributions with raw weights `(1-f,f)`. If the bins have different R or T, their
conditional contributions must be weighted accordingly. This distinction is particularly
important near the critical angle.

`lutJoint()` preserves every two-dimensional cell and applies the same interpolate-first,
condition-second rule:

```text
P(i,j | branch) =
  [(1-f) N(k0,i,j)/L(k0) + f N(k1,i,j)/L(k1)] / P_branch.
```

`monitorJoint()` divides measured counts by the total hits in that branch and performs the reverse
phi mapping described above. The script compares the two joint probability distributions using
total variation distance:

```text
TVD = 0.5 * sum_ij |P_monitor(i,j) - P_LUT(i,j)|.
```

The acceptance limit is three times an approximation of the mean statistical TVD for the expected
multinomial distribution at the observed branch hit count, with a minimum of `0.01`. The console
reports reflected and transmitted `TVD / limit / PASS|FAIL` independently. This prevents a
low-statistics transmission branch from being rejected merely for its expected sampling noise.

#### G. Graph-window output

The one-dimensional overlay contains four datasets:

- blue solid line: LUT-expected reflected theta density;
- blue points: geometry-monitor measured reflected theta density;
- green solid line: LUT-expected transmitted theta density;
- green points: geometry-monitor measured transmitted theta density.

The plot is added to the Graph basket. The script then draws four `colz` two-dimensional plots;
each expected/measured pair uses a common color scale:

- reflected LUT expected and monitor measured;
- transmitted LUT expected and monitor measured.

All four are added to the Graph basket and saved automatically as:

```text
geom_gui_out/joint_<direction>_<angle>_reflected_expected.png
geom_gui_out/joint_<direction>_<angle>_reflected_measured.png
geom_gui_out/joint_<direction>_<angle>_transmitted_expected.png
geom_gui_out/joint_<direction>_<angle>_transmitted_measured.png
```

Finite Monte Carlo fluctuations are expected in these images. Pass/fail is determined by the
joint-TVD limits printed in the console, not by requiring visual pixel-by-pixel agreement. The
direct production-rule checker in Section 1 remains the principal automated test across multiple
angles, rotated normals, and uncommon branches.

## 4. Status of the upstream `dev` monitor implementation

Before this work, `andrmor/ANTS3bundle` `dev` at commit `cc1aa4bc` (fetched 2026-07-14) was compared
with the local branch. Its base photon-monitor implementation (`amonitor*`, `amonitorhub*`, the
geometry monitor delegate, and the script API) was already the same as the implementation then
present here. Recent upstream additions were primarily the Sensor Interface/PDE tester. A wholesale
merge would conflict with local GUI changes and, because upstream does not contain DavisLUT, would
also remove DavisLUT files; no such merge was required.

The local branch has since extended the monitor with the `AnglePhi` joint histogram, JSON
serialization, and the `lsim.getMonitorAnglePhi()` API used here. This extension must not be
assumed to exist in the current upstream `dev` branch.

## File index

- `validate_rule.js`, `check_runtime.py`, `plot_runtime_validation.py`: production runtime-rule
  validation;
- `validate_geom.js`, `check_geometry.py`, `plot_geom.py`: automated complete-geometry validation;
- `validate_geom_gui.json`: complete ideal experiment configuration loadable by the GUI;
- `validate_geom_gui.txt`: GUI JavaScript runner, readout, and visualization;
- `validate_lut.cpp`, `plot_validation.py`, `plot_lego.py`: supplementary low-level sampler tests;
- `data/`: numerical results;
- `figures/`: generated or reference figures.
