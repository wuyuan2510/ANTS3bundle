# Validating the DavisLUT rule

Two complementary checks that the runtime rule reproduces the LUT it samples from:
1. a **literal geometry experiment** — fire collimated light at a LYSO->air interface at a known
   incidence angle, catch the reflected photons with a detector, and compare their measured
   polar-angle (θ_out) distribution to the LUT (validates the full tracer -> rule -> sensor chain);
2. a **rule-level Monte-Carlo check** — sample the rule directly to verify the full 2D
   (θ_out, φ_out) reflected distribution, which the flat geometry monitor cannot resolve.

## The experiment (`validate_geom.js`, run headless with `ants3 -j`)

**Geometry** (all in mm; surface normal is +z):
- An air world box; inside it a **LYSO block** `[100 x 100 x 20]` placed so its **top face sits
  at z = 0** (the block spans z = -20..0). The top face is the **LYSO -> air interface**, where
  the `DavisLUT` rule acts.
- Materials from the loaded config: LYSO = material 1 (n = 1.824), air = material 0 (n = 1.0).
- The forward LUT (`pooled_28um_fwd.lut`) is assigned to LYSO -> air and the reverse LUT
  (`pooled_28um_rev.lut`) to air -> LYSO.

**Source** — a single collimated pencil beam of 3e5 photons at `(0, 0, -0.05)`, i.e. 0.05 mm
below the interface, aimed **upward** with direction `(sin θ, 0, cos θ)` so it strikes the
interface at incidence angle **θ = 10, 45, 80 deg** (via `PhGenOverrides.Direction = Fixed`).
The beam lies in the x-z plane, so the **incidence plane is x-z**.

**Detector** — a thin **photon monitor** (48 mm half-size, horizontal) placed at **z = -0.1 mm**,
just below the interface and below the source, inside the LYSO. Because the beam is emitted
upward, only **reflected** photons (which travel back downward after the interface) reach the
monitor. The monitor records, per photon, the angle `acos(N . v)` between its face normal
(vertical) and the photon direction — i.e. the **polar angle θ_out from the surface normal** —
into a 1D histogram (90 bins over 0..90 deg). It stops tracking on hit, so each reflected photon
is counted once. Transmitted photons go up into the air and never reach the monitor.

So the experiment measures the **θ_out distribution of reflected photons** at each incidence
angle, through the full ANTS3 chain: tracer computes the interface normal -> invokes the
DavisLUT rule -> reflected direction -> sensor. It is the direct analogue of the paper's
"collimated beam + hemispherical detector" measurement (fig. 8).

**Run:**
```bash
ants3 -j validate_geom.js      # needs a built bin/ants3; writes data/geom_{10,45,80}.txt
python3 plot_geom.py           # figures/validation_geom_linear.png, validation_geom_polar.png
```

**Result.** The measured θ_out distribution overlays the LUT prediction across the whole range
at all three angles (mean θ_out within ~0.6 deg), including the strong grazing component near
90 deg. The one systematic difference is a small **deficit right at θ_out -> 90 deg**: those
reflected photons are nearly horizontal, travel a long lateral distance, and a few percent leave
the block through its sides before reaching the finite monitor. That is a real
detector-collection effect (grazing-reflected light skims along the surface), not a rule error.

## Rule-level Monte-Carlo check — the 2D (θ_out, φ_out) view (`validate_lut.cpp`)

The geometry monitor gives the θ_out marginal cleanly but **not** the joint (θ_out, φ_out)
distribution: a flat horizontal monitor records `acos(N·v)` with N vertical, i.e. only the polar
angle θ_out — it integrates over the azimuth φ_out (and ANTS3 has no hemispherical monitor
shape). To check the **full 2D** reflected distribution (azimuth included, and the θ_out = 90 deg
ridge), this second check drives the **real runtime sampling** directly:

`validate_lut.cpp` calls `ALutSurfaceData::selectThetaBin` + `sampleOutgoing` (exactly what
`ALutInterfaceRule::calculate()` invokes) for 2e6 photons per angle, reconstructs each reflected
photon's (θ_out, φ_out) with the same frame `calculate()` uses, and histograms them — a direct
per-interaction Monte-Carlo, independent of any geometry/monitor. Build against the compiled
object files, then run:
```bash
cd <repo>/ants3bundle/src/ants3            # object files (*.o) live here after building ants3
g++ -O2 -std=c++17 -I. -Itools -IphotonSim -IphotonSim/interfaceRules $(root-config --cflags) \
    -I<Qt>/include -I<Qt>/include/QtCore \
    <this>/validate_lut.cpp alutsurfacedata.o ajsontools.o aerrorhub.o \
    -L<Qt>/lib -lQt6Core -o validate_lut
./validate_lut <this>/../luts/pooled_28um_fwd.lut <this>/data   # writes data/valid_{10,45,80}.txt
```
Then plot:
```bash
python3 plot_lego.py           # figures/validation_heatmap.png  (3D lego, LUT surface vs measured wireframe, log z)
python3 plot_validation.py     # figures/validation_polar.png    (theta_out marginal, polar)
```
**Result.** Measured vs stored-LUT agree to Monte-Carlo noise (θ_out total-variation distance
~0.002), reflectance R matches to ~0.001, and the θ_out = 90 deg bin fraction matches to three
digits (0.124 / 0.125 / 0.149 for 10 / 45 / 80 deg) — confirming the rule reproduces the full 2D
distribution, and that the 90 deg feature is a faithful, deterministic result of the
medium-based escape handling (see the main notes). The bundled `data/valid_*.txt` let the two
plotters regenerate the figures without rebuilding/rerunning.

The two checks are complementary: the geometry experiment confirms the θ_out marginal through the
full tracer -> rule -> sensor chain, and this MC check confirms the full 2D angular distribution
against the LUT.

## Files
- `validate_geom.js`, `plot_geom.py` — geometry experiment and its θ_out plots (linear + polar)
- `validate_lut.cpp`, `plot_lego.py`, `plot_validation.py` — rule-level MC check: 2D lego overlay + polar
- `data/` — intermediate outputs, so the plotters regenerate the figures immediately
- `figures/` — the reference figures
