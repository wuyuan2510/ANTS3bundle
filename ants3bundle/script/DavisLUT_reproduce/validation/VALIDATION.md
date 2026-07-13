# Validating the DavisLUT rule

Two complementary checks that the runtime rule reproduces the LUT it was built from, and that
the reflected angular distribution (including the grazing / theta_out=90 deg feature) is
faithful. Both use the 28 um ensemble LUT in `../luts/`.

## 1. Literal geometry experiment (source -> interface -> detector)

`validate_geom.js` — a full headless ANTS3 photon simulation. A collimated source inside a
LYSO block fires at a fixed incidence angle (10, 45, 80 deg) at the LYSO->air interface (where
the DavisLUT rule acts); a thin photon monitor just below the interface catches the reflected
photons and records their polar angle theta_out. This exercises the whole chain
tracer-normal -> rule -> reflected direction -> sensor.

Run (needs a built `bin/ants3`), from this directory:
```bash
ants3 -j validate_geom.js         # writes data/geom_{10,45,80}.txt
python3 plot_geom.py              # figures/validation_geom_linear.png, validation_geom_polar.png
```
Result: the measured theta_out distribution overlays the LUT prediction across the whole range,
including the grazing spike; the only difference is a small deficit at theta_out -> 90 deg,
where near-horizontal reflected photons skim the surface and a few percent leave the block
before reaching the finite monitor (a real detector-collection effect).

## 2. Rule-level Monte-Carlo check (samples the exact runtime code)

`validate_lut.cpp` — drives the real runtime sampling (`ALutSurfaceData::selectThetaBin` +
`sampleOutgoing`, i.e. what `ALutInterfaceRule::calculate()` calls) for 2e6 photons per angle
and reconstructs (theta_out, phi_out) of the reflected photons. Compile against the built
object files, e.g.:
```bash
cd <repo>/ants3bundle/src/ants3
g++ -O2 -std=c++17 -I. -Itools -IphotonSim -IphotonSim/interfaceRules $(root-config --cflags) \
    -I<Qt>/include -I<Qt>/include/QtCore \
    <this>/validate_lut.cpp alutsurfacedata.o ajsontools.o aerrorhub.o \
    -L<Qt>/lib -lQt6Core -std=c++17 -o validate_lut
./validate_lut <this>/../luts/pooled_28um_fwd.lut <this>/data     # writes data/valid_{10,45,80}.txt
```
Then plot:
```bash
python3 plot_validation.py        # figures/validation_polar.png  (measured vs LUT, polar)
python3 plot_lego.py              # figures/validation_heatmap.png (3D lego overlay, log z)
```
Result: measured vs stored-LUT agree to Monte-Carlo noise (theta_out total-variation distance
~0.002), reflectance R matches to ~0.001, and the theta_out=90 deg bin fraction matches to
three digits (0.124 / 0.125 / 0.149 for 10 / 45 / 80 deg) — confirming the 90 deg feature is a
faithful, deterministic consequence of the medium-based escape handling (see the main notes).

## Files
- `validate_geom.js`, `plot_geom.py` — geometry experiment + plots
- `validate_lut.cpp`, `plot_validation.py`, `plot_lego.py` — rule-level check + plots
- `data/` — the intermediate outputs, so the plotters regenerate the figures immediately
- `figures/` — the reference figures
