#!/usr/bin/env python3
"""Numerically check the dual-monitor geometry experiment; exits non-zero on failure."""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "geometry")
LUTS = {
    "fwd": os.path.join(HERE, "..", "luts", "pooled_28um_fwd.lut"),
    "rev": os.path.join(HERE, "..", "luts", "pooled_28um_rev.lut"),
}


def angle_tag(angle):
    return f"{angle:g}".replace(".", "p")


def expected(data, angle, reflected):
    n_inc = data["Binning"]["ThetaIncBins"]
    n_theta = data["Binning"]["ThetaOutBins"]
    n_phi = data["Binning"]["PhiOutBins"]
    x = angle/(90.0/n_inc)-0.5
    k0, frac = int(np.floor(x)), x-math.floor(x)
    if k0 < 0:
        k0, frac = 0, 0.0
    if k0 >= n_inc-1:
        k0, frac = n_inc-1, 0.0
    k1 = k0+1 if frac > 0 else k0
    weights = ((k0, 1.0-frac),) if k0 == k1 else ((k0, 1.0-frac), (k1, frac))
    counts = data["ReflectedCounts" if reflected else "TransmittedCounts"]
    histograms = data["ReflectedHist" if reflected else "TransmittedHist"]
    probability = sum(w*counts[k]/data["Launched"][k] for k, w in weights)
    marginal = np.zeros(n_theta)
    for k, weight in weights:
        marginal += weight*np.asarray(histograms[k], float).reshape(n_theta, n_phi).sum(axis=1)/data["Launched"][k]
    if probability:
        marginal /= probability
    return probability, marginal


def noise_tvd_limit(expected_distribution, samples):
    if samples < 1:
        return 1.0
    mean_noise = 0.5*sum(math.sqrt(max(0.0, 2*p*(1-p)/(math.pi*samples))) for p in expected_distribution)
    return max(0.008, 3.0*mean_noise)


with open(os.path.join(DATA, "summary.json"), encoding="utf-8") as stream:
    summary = json.load(stream)
launched = summary["photonsPerRun"]
all_passed = True
metrics = []

for direction in ("fwd", "rev"):
    with open(LUTS[direction], encoding="utf-8") as stream:
        lut = json.load(stream)
    for run in (r for r in summary["runs"] if r["direction"] == direction):
        angle = run["angle"]
        missing_fraction = abs(run["missing"])/launched
        run_passed = missing_fraction <= 0.001
        for outcome, reflected in (("reflected", True), ("transmitted", False)):
            probability, prediction = expected(lut, angle, reflected)
            hits = run[outcome]
            measured_probability = hits/launched
            probability_limit = max(0.001, 6*math.sqrt(max(probability*(1-probability), 1/launched)/launched))
            filename = os.path.join(DATA, f"geom_{direction}_{outcome}_{angle_tag(angle)}.txt")
            counts = np.loadtxt(filename)[:, 1]
            # The monitor has 90 one-degree bins; LUT has 45 two-degree bins.
            observation = counts.reshape(len(prediction), -1).sum(axis=1)
            if observation.sum():
                observation /= observation.sum()
            tvd = 0.5*np.abs(observation-prediction).sum() if hits else 0.0
            tvd_limit = noise_tvd_limit(prediction, hits)
            passed = bool(abs(measured_probability-probability) <= probability_limit and tvd <= tvd_limit)
            run_passed = run_passed and passed
            metrics.append({"direction": direction, "angle": angle, "outcome": outcome,
                            "measuredProbability": measured_probability, "expectedProbability": probability,
                            "probabilityLimit": probability_limit, "thetaTvd": tvd, "thetaTvdLimit": tvd_limit,
                            "missingFraction": missing_fraction, "passed": passed})
        all_passed = all_passed and run_passed
        print(f"{direction:3s} {angle:4.1f} deg: missing={missing_fraction:.5f}  {'PASS' if run_passed else 'FAIL'}")

with open(os.path.join(DATA, "metrics.json"), "w", encoding="utf-8") as stream:
    json.dump({"passed": all_passed, "metrics": metrics}, stream, indent=2)

print("Geometry validation", "PASSED" if all_passed else "FAILED")
sys.exit(0 if all_passed else 1)
