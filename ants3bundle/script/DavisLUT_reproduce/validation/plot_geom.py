#!/usr/bin/env python3
"""Plot the ideal dual-monitor geometry experiment produced by validate_geom.js."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "geometry")
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
LUTS = {
    "fwd": os.path.join(HERE, "..", "luts", "pooled_28um_fwd.lut"),
    "rev": os.path.join(HERE, "..", "luts", "pooled_28um_rev.lut"),
}
SELECTED = (10.0, 33.3, 45.0, 80.0)


def angle_tag(angle):
    return f"{angle:g}".replace(".", "p")


def incidence_mix(data, angle):
    n_inc = data["Binning"]["ThetaIncBins"]
    x = angle/(90.0/n_inc)-0.5
    k0 = int(np.floor(x))
    frac = x-k0
    if k0 < 0:
        return 0, 0, 0.0
    if k0 >= n_inc-1:
        return n_inc-1, n_inc-1, 0.0
    return k0, k0+1 if frac > 0 else k0, frac


def expected(data, angle, reflected):
    n_theta = data["Binning"]["ThetaOutBins"]
    n_phi = data["Binning"]["PhiOutBins"]
    counts = data["ReflectedCounts" if reflected else "TransmittedCounts"]
    histograms = data["ReflectedHist" if reflected else "TransmittedHist"]
    k0, k1, frac = incidence_mix(data, angle)
    weights = ((k0, 1.0-frac),) if k0 == k1 else ((k0, 1.0-frac), (k1, frac))
    probability = sum(weight*counts[k]/data["Launched"][k] for k, weight in weights)
    marginal = np.zeros(n_theta)
    for k, weight in weights:
        hist = np.asarray(histograms[k], dtype=float).reshape(n_theta, n_phi).sum(axis=1)
        marginal += weight*hist/data["Launched"][k]
    if probability > 0:
        marginal /= probability
    return probability, marginal


def measured(direction, outcome, angle):
    filename = os.path.join(DATA, f"geom_{direction}_{outcome}_{angle_tag(angle)}.txt")
    values = np.loadtxt(filename)
    counts = values[:, 1]
    return values[:, 0], counts/counts.sum() if counts.sum() else counts


with open(os.path.join(DATA, "summary.json"), encoding="utf-8") as stream:
    summary = json.load(stream)
launched = summary["photonsPerRun"]

# Absolute probabilities.
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
for ax, direction in zip(axes, ("fwd", "rev")):
    with open(LUTS[direction], encoding="utf-8") as stream:
        lut = json.load(stream)
    records = [r for r in summary["runs"] if r["direction"] == direction]
    angles = np.array([r["angle"] for r in records])
    for outcome, reflected, color in (("R", True, "tab:blue"), ("T", False, "tab:orange")):
        prediction = np.array([expected(lut, angle, reflected)[0] for angle in angles])
        key = "reflected" if reflected else "transmitted"
        observation = np.array([r[key]/launched for r in records])
        ax.plot(angles, prediction, color=color, lw=2, label=f"{outcome} LUT")
        ax.plot(angles, observation, "o", color=color, fillstyle="none", label=f"{outcome} geometry")
    ax.set_title(direction)
    ax.set_xlabel("incidence angle, deg")
    ax.grid(alpha=0.3)
    ax.set_ylim(-0.02, 1.02)
axes[0].set_ylabel("probability")
axes[0].legend(fontsize=8)
fig.suptitle("DavisLUT ideal geometry experiment: absolute R/T", y=0.98)
fig.tight_layout(rect=(0, 0, 1, 0.90))
fig.savefig(os.path.join(FIG, "validation_geom_probabilities.png"), dpi=150)

# Conditional theta distributions.
for direction in ("fwd", "rev"):
    with open(LUTS[direction], encoding="utf-8") as stream:
        lut = json.load(stream)
    n_theta = lut["Binning"]["ThetaOutBins"]
    theta_lut = (np.arange(n_theta)+0.5)*90.0/n_theta
    theta_width = 90.0/n_theta
    fig, axes = plt.subplots(2, len(SELECTED), figsize=(15, 7.6), sharex=True)
    for col, angle in enumerate(SELECTED):
        for row, (outcome, reflected) in enumerate((("reflected", True), ("transmitted", False))):
            _, prediction = expected(lut, angle, reflected)
            theta_measured, observation = measured(direction, outcome, angle)
            measured_width = theta_measured[1]-theta_measured[0]
            ax = axes[row, col]
            ax.plot(theta_lut, prediction/theta_width, color="tab:blue", lw=2, label="LUT")
            ax.step(theta_measured, observation/measured_width, where="mid", color="tab:red", lw=1.2,
                    label="geometry monitor")
            ax.set_title(f"{angle:g} deg")
            ax.grid(alpha=0.25)
            if col == 0:
                ax.set_ylabel(f"{outcome}\nprobability density, 1/deg")
            if row == 1:
                ax.set_xlabel("theta_out, deg")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(f"DavisLUT ideal geometry experiment: {direction}", y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.91), h_pad=2.0)
    fig.savefig(os.path.join(FIG, f"validation_geom_{direction}.png"), dpi=150)

print("wrote figures/validation_geom_probabilities.png and validation_geom_{fwd,rev}.png")
