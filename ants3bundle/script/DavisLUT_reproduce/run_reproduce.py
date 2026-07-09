#!/usr/bin/env python3
"""Reproduce the DavisLUT vs old-microfacet DOI comparison from the bundled pooled LUTs.

For each roughness (5 / 14 / 28 um) it builds a config with the DavisLUT rule on the
LYSO<->air interface, runs the headless photon worker `lsim` over four source depths, and
plots (a) DOI asymmetry vs depth (paper LUT vs the bundled old-method reference) and
(b) R(theta) for the three roughnesses. Everything is resolved relative to this file, so it
runs from a fresh checkout after ANTS3 is built.

    python3 run_reproduce.py [events_per_depth=500]

Needs: a built `lsim` at <repo>/ants3bundle/bin/lsim, numpy + matplotlib.
This reproduces the *simulation* results from the LUTs; regenerating the LUTs from AFM data
is a separate step documented in REPRODUCE.md.
"""
import json, os, sys, subprocess
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
LSIM = os.path.normpath(os.path.join(HERE, "..", "..", "bin", "lsim"))
BASE_CONFIG = os.path.join(HERE, "config", "config-28Yuan-airGaps-04-LUT.json")
DEPTHS = [-13, -10, -5, 0]
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 500     # paper runs used 1000
ROUGH = {  # tag -> (fwd lut, rev lut, colour)
    "5um":  ("pooled_5um_fwd.lut",  "pooled_5um_rev.lut",  "tab:blue"),
    "14um": ("pooled_14um_fwd.lut", "pooled_14um_rev.lut", "tab:green"),
    "28um": ("pooled_28um_fwd.lut", "pooled_28um_rev.lut", "tab:red"),
}

def doi(path):
    a = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"): continue
        p = line.split()
        if len(p) >= 2 and float(p[0]) + float(p[1]) > 0:
            a.append((float(p[0]) - float(p[1])) / (float(p[0]) + float(p[1])))
    a = np.array(a); return a.mean(), a.std()

def lut_rule(mf, mt, lut, fn):
    return {"MatFrom": mf, "MatTo": mt, "Model": "DavisLUT", "Symmetric": False,
            "SurfaceProperties": {"Model": "Polished", "KillPhotonsRefractedBackward": False,
                                  "OrientationProbabilityCorrection": False},
            "LutFile": fn, "LUT": lut}

if not os.path.isfile(LSIM):
    sys.exit(f"lsim not found at {LSIM} -- build ANTS3 first (see REPRODUCE.md).")

paper = {}
for tag, (fwdf, revf, _) in ROUGH.items():
    fwd = json.load(open(os.path.join(HERE, "luts", fwdf)))
    rev = json.load(open(os.path.join(HERE, "luts", revf)))
    paper[tag] = {}
    for z in DEPTHS:
        wd = os.path.join(HERE, "repro_output", tag, str(z)); os.makedirs(wd, exist_ok=True)
        c = json.load(open(BASE_CONFIG))
        c["InterfaceRules"]["MaterialRules"] = [
            lut_rule(1, 4, fwd, fwdf) if (r["MatFrom"], r["MatTo"]) == (1, 4) else
            lut_rule(4, 1, rev, revf) if (r["MatFrom"], r["MatTo"]) == (4, 1) else r
            for r in c["InterfaceRules"]["MaterialRules"]]
        c["PhotonSim"]["PhotonBombs"]["Flood"]["Zfixed"] = z
        c["PhotonSim"]["PhotonBombs"]["Flood"]["Number"] = NB
        c["PhotonSim"]["Run"]["EventFrom"] = 0
        c["PhotonSim"]["Run"]["EventTo"]   = NB
        c["PhotonSim"]["Run"]["Seed"]      = 3000 + abs(z)
        json.dump(c, open(os.path.join(wd, "conf.json"), "w"))
        subprocess.run([LSIM, wd, "conf.json", "0"], cwd=wd,
                       stdout=open(os.path.join(wd, "log"), "w"), stderr=subprocess.STDOUT)
        paper[tag][z] = doi(os.path.join(wd, "SensorSignals.txt"))
        print(f"{tag} z={z:>4}: paper DOI = {paper[tag][z][0]:+.4f}", flush=True)

old = {tag: {z: doi(os.path.join(HERE, "old_reference", tag, str(z), "SensorSignals.txt"))
             for z in DEPTHS} for tag in ROUGH}

resdir = os.path.join(HERE, "results_reproduced"); os.makedirs(resdir, exist_ok=True)

# --- DOI vs depth: old vs paper, three roughnesses ---
fig, ax = plt.subplots(figsize=(7.5, 5))
for tag, (_, _, col) in ROUGH.items():
    ax.errorbar(DEPTHS, [old[tag][z][0] for z in DEPTHS], yerr=[old[tag][z][1] for z in DEPTHS],
                marker='o', ls='--', color=col, capsize=3, alpha=0.7, label=f"{tag} old microfacet")
    ax.errorbar(DEPTHS, [paper[tag][z][0] for z in DEPTHS], yerr=[paper[tag][z][1] for z in DEPTHS],
                marker='s', ls='-', color=col, capsize=3, label=f"{tag} LUT (paper)")
ax.set_xlabel("source depth Zfixed (mm)"); ax.set_ylabel("DOI asymmetry (s0-s1)/(s0+s1)")
ax.set_title(f"DOI vs depth: old microfacet (dashed) vs LUT (solid), {NB} ev")
ax.grid(alpha=.3); ax.legend(fontsize=8, ncol=3); fig.tight_layout()
fig.savefig(os.path.join(resdir, "pooled_DOI_old_vs_LUT.png"), dpi=140)

# --- R(theta), three roughnesses ---
def Rof(fn):
    d = json.load(open(os.path.join(HERE, "luts", fn))); n = d["Binning"]["ThetaIncBins"]
    return (np.arange(n) + 0.5) * 90 / n, np.array(d["ReflectedCounts"], float) / np.array(d["Launched"], float)
def fres(th, n1, n2):
    th = np.radians(th); c = np.cos(th); s2 = np.sin(th)**2; nsq = (n2/n1)**2; o = np.ones_like(th)
    ok = s2 <= nsq; f2 = np.sqrt(np.clip(nsq - s2, 0, None))
    o[ok] = 0.5*(((nsq*c - f2)/(nsq*c + f2))**2 + ((c - f2)/(c + f2))**2)[ok]; return o
fig, ax = plt.subplots(figsize=(7, 4.8)); thf = np.linspace(0, 90, 400)
ax.plot(thf, fres(thf, 1.824, 1.0), 'k--', lw=1.4, label='flat Fresnel')
for tag, (fwdf, _, col) in ROUGH.items():
    th, r = Rof(fwdf); ax.plot(th, r, '-o', ms=3, color=col, label=f"{tag} LUT")
ax.axvline(np.degrees(np.arcsin(1/1.824)), color='gray', ls=':', lw=1)
ax.set_xlabel("incidence angle (deg)"); ax.set_ylabel("reflection probability R")
ax.set_title("LYSO->air reflectance (paper LUT): 5 / 14 / 28 um"); ax.set_ylim(0, 1.05)
ax.grid(alpha=.3); ax.legend(fontsize=9, loc='lower right'); fig.tight_layout()
fig.savefig(os.path.join(resdir, "pooled_Rtheta.png"), dpi=140)

print("\nwrote plots to", resdir)
print(f"{'depth':>6} | " + " | ".join(f"{t} old/LUT" for t in ROUGH))
for z in DEPTHS:
    print(f"{z:>6} | " + " | ".join(f"{old[t][z][0]:+.3f}/{paper[t][z][0]:+.3f}" for t in ROUGH))
