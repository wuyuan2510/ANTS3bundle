#!/usr/bin/env python3
"""Plot results produced by validate_rule.js (actual ALutInterfaceRule runtime path)."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "runtime")
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

results = {}
for tag in ("fwd", "rev"):
    with open(os.path.join(DATA, f"{tag}.json"), encoding="utf-8") as stream:
        results[tag] = json.load(stream)

# R/T/A probabilities and their Monte-Carlo measurements.
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
for ax, tag in zip(axes, ("fwd", "rev")):
    records = results[tag]["angles"]
    angle = np.array([r["angle"] for r in records])
    for outcome, color in (("R", "tab:blue"), ("T", "tab:orange"), ("A", "tab:gray")):
        expected = np.array([r[f"expected{outcome}"] for r in records])
        measured = np.array([r[f"measured{outcome}"] for r in records])
        ax.plot(angle, expected, color=color, lw=2, label=f"{outcome} LUT")
        ax.plot(angle, measured, "o", color=color, ms=4, fillstyle="none", label=f"{outcome} runtime")
    ax.set_title(f"{tag}: {'LYSO -> air' if tag == 'fwd' else 'air -> LYSO'}")
    ax.set_xlabel("incidence angle, deg")
    ax.grid(alpha=0.3)
    ax.set_ylim(-0.02, 1.02)
axes[0].set_ylabel("probability")
axes[0].legend(ncol=2, fontsize=8)
fig.suptitle("DavisLUT actual-rule validation: outcome probabilities", y=0.98)
fig.tight_layout(rect=(0, 0, 1, 0.90))
fig.savefig(os.path.join(FIG, "runtime_probabilities.png"), dpi=150)

# Conditional theta_out marginals for both outcomes and both directions.
selected = (10.0, 33.3, 45.0, 80.0)
for tag in ("fwd", "rev"):
    records = results[tag]["angles"]
    by_angle = {round(r["angle"], 6): r for r in records}
    n_theta = results[tag]["thetaOutBins"]
    theta = (np.arange(n_theta) + 0.5)*90.0/n_theta
    fig, axes = plt.subplots(2, len(selected), figsize=(15, 7.6), sharex=True)
    for col, angle in enumerate(selected):
        rec = by_angle[round(angle, 6)]
        for row, outcome in enumerate(("Reflected", "Transmitted")):
            ax = axes[row, col]
            expected = np.asarray(rec[f"expected{outcome}Theta"])
            measured = np.asarray(rec[f"measured{outcome}Theta"])
            tvd = rec[f"{outcome.lower()}ThetaTvd"]
            ax.plot(theta, expected, color="tab:blue", lw=2, label="LUT")
            ax.plot(theta, measured, "o", color="tab:red", ms=2.5, label="actual rule")
            ax.set_title(f"{angle:g} deg, TVD={tvd:.4f}")
            ax.grid(alpha=0.25)
            if col == 0:
                ax.set_ylabel(f"{outcome}\nprobability / bin")
            if row == 1:
                ax.set_xlabel("theta_out, deg")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(f"DavisLUT {tag}: actual ALutInterfaceRule theta_out marginals", y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.91), h_pad=2.0)
    fig.savefig(os.path.join(FIG, f"runtime_theta_{tag}.png"), dpi=150)

# A quantitative view of the full two-dimensional theta/phi agreement.
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
for ax, tag in zip(axes, ("fwd", "rev")):
    records = results[tag]["angles"]
    angle = np.array([r["angle"] for r in records])
    for outcome, color in (("reflected", "tab:blue"), ("transmitted", "tab:orange")):
        tvd = np.array([r[f"{outcome}2dTvd"] for r in records])
        limit = np.array([r[f"{outcome}2dTvdLimit"] for r in records])
        ax.plot(angle, tvd, "o-", color=color, label=f"{outcome} TVD")
        ax.plot(angle, limit, "--", color=color, alpha=0.55, label=f"{outcome} limit")
    ax.set_title(tag)
    ax.set_xlabel("incidence angle, deg")
    ax.grid(alpha=0.3)
axes[0].set_ylabel("full 2D total-variation distance")
axes[0].legend(fontsize=8)
fig.suptitle("DavisLUT actual-rule 2D validation and statistical limits", y=0.98)
fig.tight_layout(rect=(0, 0, 1, 0.90))
fig.savefig(os.path.join(FIG, "runtime_2d_tvd.png"), dpi=150)

print("wrote figures/runtime_probabilities.png, runtime_theta_{fwd,rev}.png, runtime_2d_tvd.png")
