#!/usr/bin/env python3
"""Validate DavisLUT through complete geometry transport in fresh lsim processes.

Each direction/angle case receives its own static configuration and worker
process. This avoids state races caused by repeatedly mutating the GUI config
while dispatcher simulations are still being handled.

Usage:
    python3 validate_geom.py
    python3 validate_geom.py --photons 100000 --angles 0 45 80
    python3 validate_geom.py --fwd /path/fwd.lut --rev /path/rev.lut
"""

import argparse
import copy
import json
import math
import subprocess
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
BUNDLE = HERE.parents[2]
DEFAULT_LSIM = BUNDLE / "bin/lsim"
DEFAULT_CONFIG = HERE / "validate_geom_gui.json"
DEFAULT_FWD = HERE.parent / "luts/pooled_28um_fwd.lut"
DEFAULT_REV = HERE.parent / "luts/pooled_28um_rev.lut"


def interpolated_expected(lut, angle, reflected):
    bins = lut["Binning"]
    n_inc = bins["ThetaIncBins"]
    n_theta = bins["ThetaOutBins"]
    n_phi = bins["PhiOutBins"]
    coordinate = angle / (90.0 / n_inc) - 0.5
    low = math.floor(coordinate)
    fraction = coordinate - low
    if low < 0:
        low, fraction = 0, 0.0
    if low >= n_inc - 1:
        low, fraction = n_inc - 1, 0.0
    high = low + 1 if fraction > 0 else low
    weights = ((low, 1.0 - fraction),) if low == high else (
        (low, 1.0 - fraction), (high, fraction))

    count_key = "ReflectedCounts" if reflected else "TransmittedCounts"
    histogram_key = "ReflectedHist" if reflected else "TransmittedHist"
    probability = sum(
        weight * lut[count_key][index] / lut["Launched"][index]
        for index, weight in weights)
    distribution = np.zeros((n_theta, n_phi))
    for index, weight in weights:
        histogram = np.asarray(
            lut[histogram_key][index], dtype=float).reshape(n_theta, n_phi)
        distribution += weight * histogram / lut["Launched"][index]
    if probability:
        distribution /= probability
    return probability, distribution


def noise_tvd_limit(distribution, samples):
    if samples < 1:
        return 1.0
    flat = np.asarray(distribution).ravel()
    mean_noise = 0.5 * sum(
        math.sqrt(max(0.0, 2 * p * (1 - p) / (math.pi * samples)))
        for p in flat)
    return max(0.008, 3.0 * mean_noise)


def replace_luts(config, fwd, rev, fwd_path, rev_path):
    found = set()
    for rule in config["InterfaceRules"]["MaterialRules"]:
        pair = (rule["MatFrom"], rule["MatTo"])
        if pair == (1, 0):
            rule["LUT"] = copy.deepcopy(fwd)
            rule["LutFile"] = str(fwd_path)
            found.add(pair)
        elif pair == (0, 1):
            rule["LUT"] = copy.deepcopy(rev)
            rule["LutFile"] = str(rev_path)
            found.add(pair)
    if found != {(1, 0), (0, 1)}:
        raise RuntimeError(
            "validation config must contain LYSO->air and air->LYSO LUT rules")


def monitor_distributions(monitor):
    angle = monitor["Angle"]
    theta = np.asarray(
        [record[1] for record in angle["Data"][1:-1]], dtype=float)
    if len(theta) != angle["Bins"]:
        raise RuntimeError("unexpected theta histogram shape")

    joint_record = monitor.get("AnglePhi")
    if not joint_record:
        raise RuntimeError(
            "monitor has no AnglePhi data; ANTS3 with 2D monitor support is required")
    rows = [
        record for record in joint_record["Data"]
        if 0 < record[0] < 90 and 0 < record[1] < 360
    ]
    joint = np.asarray([record[2] for record in rows], dtype=float)
    if joint.size != 45 * 36:
        raise RuntimeError(
            "expected a 45 x 36 AnglePhi histogram, got {} cells".format(
                joint.size))
    joint = joint.reshape(45, 36)
    hits = int(angle["Entries"])
    if int(theta.sum()) != hits or int(joint.sum()) != hits:
        raise RuntimeError("monitor histogram entries do not agree")
    return hits, theta, joint


def write_theta(path, counts, direction, angle, outcome, hits):
    with path.open("w") as stream:
        stream.write(
            "# direction={} angle={} outcome={} total_hits={}\n".format(
                direction, angle, outcome, hits))
        for index, count in enumerate(counts):
            stream.write("{} {}\n".format(2 * index + 1, int(count)))


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fwd", type=Path, default=DEFAULT_FWD)
    parser.add_argument("--rev", type=Path, default=DEFAULT_REV)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--lsim", type=Path, default=DEFAULT_LSIM)
    parser.add_argument("--photons", type=int, default=300000)
    parser.add_argument(
        "--angles", type=float, nargs="+",
        default=[0, 10, 32, 33.3, 34, 45, 80, 89])
    parser.add_argument(
        "--output-root", type=Path, default=HERE / "geom_out/batch")
    parser.add_argument(
        "--data-dir", type=Path, default=HERE / "data/geometry")
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.photons < 1:
        raise SystemExit("--photons should be positive")
    if not args.lsim.is_file():
        raise SystemExit("lsim not found: {}".format(args.lsim))

    args.output_root.mkdir(parents=True, exist_ok=True)
    args.data_dir.mkdir(parents=True, exist_ok=True)
    base = json.load(args.config.open())
    fwd = json.load(args.fwd.open())
    rev = json.load(args.rev.open())
    replace_luts(base, fwd, rev, args.fwd, args.rev)
    base["PhotonSim"]["PhotonBombs"]["PhotonsPerBomb"]["FixedNumber"] = args.photons
    base["PhotonSim"]["Run"]["SaveMonitors"] = True
    base["PhotonSim"]["Run"]["SaveSensorSignals"] = False
    base["PhotonSim"]["Run"]["SavePhotonBombs"] = False
    base["PhotonSim"]["Run"]["SaveStatistics"] = False

    metrics = []
    summaries = []
    all_passed = True
    for direction, reverse in (("fwd", False), ("rev", True)):
        lut = fwd if direction == "fwd" else rev
        for angle_index, angle in enumerate(args.angles):
            angle_tag = "{:g}".format(angle).replace(".", "p")
            case_dir = args.output_root / "{}_{}".format(direction, angle_tag)
            case_dir.mkdir(parents=True, exist_ok=True)
            config = copy.deepcopy(base)
            radians = math.radians(angle)
            depth = 0.002
            config["PhotonSim"]["PhotonBombs"]["Single"]["Position"] = [
                -depth * math.tan(radians), 0, depth if reverse else -depth]
            config["PhotonSim"]["PhGenOverrides"]["Direction"] = {
                "ConeAngle": 10,
                "DirectionMode": "Fixed",
                "DirectionVector": [
                    math.sin(radians), 0,
                    -math.cos(radians) if reverse else math.cos(radians)]
            }
            config["PhotonSim"]["Run"]["OutputDirectory"] = str(case_dir)
            config["PhotonSim"]["Run"]["Seed"] = (
                20260724 + 1000 * int(reverse) + angle_index)
            config_path = case_dir / "config.json"
            json.dump(config, config_path.open("w"))

            with (case_dir / "log").open("w") as log:
                completed = subprocess.run(
                    [str(args.lsim), str(case_dir), config_path.name, "0"],
                    cwd=case_dir, stdout=log, stderr=subprocess.STDOUT)
            if completed.returncode:
                raise RuntimeError(
                    "{} {} deg: lsim exited with {}".format(
                        direction, angle, completed.returncode))

            monitor_file = case_dir / "PhotonMonitors.txt"
            monitors = json.load(monitor_file.open())["ParticleMonitorData"]
            if len(monitors) != 2:
                raise RuntimeError("expected exactly two photon monitors")
            reflected_index = 1 if reverse else 0
            transmitted_index = 0 if reverse else 1
            reflected = monitor_distributions(monitors[reflected_index])
            transmitted = monitor_distributions(monitors[transmitted_index])
            missing = args.photons - reflected[0] - transmitted[0]
            run_passed = abs(missing) / args.photons <= 0.001

            for outcome, is_reflected, measured in (
                    ("reflected", True, reflected),
                    ("transmitted", False, transmitted)):
                hits, theta_counts, joint_counts = measured
                probability, expected_joint = interpolated_expected(
                    lut, angle, is_reflected)
                measured_probability = hits / args.photons
                probability_limit = max(
                    0.001,
                    6 * math.sqrt(max(
                        probability * (1 - probability), 1 / args.photons)
                        / args.photons))

                expected_theta = expected_joint.sum(axis=1)
                observed_theta = theta_counts / hits if hits else theta_counts
                theta_tvd = (
                    0.5 * np.abs(observed_theta - expected_theta).sum()
                    if hits else 0.0)
                theta_limit = noise_tvd_limit(expected_theta, hits)

                # At exact normal incidence there is no unique incidence plane and therefore
                # no physically defined phi=0. Probability and theta remain testable.
                joint_applicable = abs(math.sin(radians)) >= 1.0e-6
                observed_joint = joint_counts / hits if hits else joint_counts
                # The reverse-facing monitor has the opposite azimuth handedness.
                if reverse:
                    observed_joint = observed_joint[:, ::-1]
                joint_tvd = (
                    0.5 * np.abs(observed_joint - expected_joint).sum()
                    if hits and joint_applicable else 0.0)
                joint_limit = (
                    noise_tvd_limit(expected_joint, hits)
                    if joint_applicable else 0.0)

                passed = bool(
                    abs(measured_probability - probability) <= probability_limit
                    and theta_tvd <= theta_limit
                    and (not joint_applicable or joint_tvd <= joint_limit))
                run_passed = run_passed and passed
                metrics.append({
                    "direction": direction,
                    "angle": angle,
                    "outcome": outcome,
                    "measuredProbability": measured_probability,
                    "expectedProbability": probability,
                    "probabilityLimit": probability_limit,
                    "thetaTvd": float(theta_tvd),
                    "thetaTvdLimit": theta_limit,
                    "jointTvd": float(joint_tvd),
                    "jointTvdLimit": joint_limit,
                    "jointApplicable": joint_applicable,
                    "missingFraction": abs(missing) / args.photons,
                    "passed": passed,
                })
                write_theta(
                    args.data_dir / "geom_{}_{}_{}.txt".format(
                        direction, outcome, angle_tag),
                    theta_counts, direction, angle, outcome, hits)

            summaries.append({
                "direction": direction,
                "angle": angle,
                "reflected": reflected[0],
                "transmitted": transmitted[0],
                "launched": args.photons,
                "missing": missing,
            })
            all_passed = all_passed and run_passed
            print("{} {:4.1f} deg missing={:4d} {}".format(
                direction, angle, missing,
                "PASS" if run_passed else "FAIL"), flush=True)

    json.dump(
        {"photonsPerRun": args.photons, "wavelength": 420, "runs": summaries},
        (args.data_dir / "summary.json").open("w"), indent=2)
    json.dump(
        {"passed": bool(all_passed), "metrics": metrics},
        (args.data_dir / "metrics.json").open("w"), indent=2)
    print("Geometry validation", "PASSED" if all_passed else "FAILED")
    raise SystemExit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
