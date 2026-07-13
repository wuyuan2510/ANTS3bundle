# Validating the DavisLUT rule

The validation is a **literal geometry experiment**: fire collimated light at a LYSO->air
interface at a known incidence angle, catch the reflected photons with a detector, and compare
their measured angular distribution to what the LUT predicts.

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

## Supplementary 2D angular view (`validation_heatmap.png`)

The monitor gives the θ_out marginal cleanly but not the joint (θ_out, φ_out) map. To also check
the **full 2D** reflected distribution (including the azimuth and the θ_out = 90 deg ridge),
`plot_lego.py` draws a 3D lego (log z): the **LUT prediction as a blue surface** and, overlaid,
a **direct per-interaction sampling of the same rule** as a red wireframe (bundled in
`data/valid_{10,45,80}.txt`). The wireframe tracks the surface everywhere, confirming the rule
reproduces the stored 2D distribution. (This 2D cross-check uses direct sampling of the rule
rather than the geometry monitor, which only records θ_out.)
```bash
python3 plot_lego.py           # figures/validation_heatmap.png
```

## Files
- `validate_geom.js`, `plot_geom.py` — the geometry experiment and its θ_out plots (linear + polar)
- `plot_lego.py` — the supplementary 2D angular lego overlay
- `data/` — intermediate outputs, so the plotters regenerate the figures immediately
- `figures/` — the reference figures
