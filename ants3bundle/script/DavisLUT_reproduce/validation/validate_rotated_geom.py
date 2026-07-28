#!/usr/bin/env python3
"""Validate DavisLUT local-to-global conversion with a rotated interface.

For each propagation direction, this script runs an unrotated baseline and a
rigidly rotated copy of the complete interface/monitor setup with the same
random seed.  It checks both cases against the LUT and checks that their
monitor-frame theta/phi distributions are rotation invariant.
"""

import argparse
import copy
import json
import math
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from validate_geom import (
    DEFAULT_CONFIG,
    DEFAULT_FWD,
    DEFAULT_LSIM,
    DEFAULT_REV,
    interpolated_expected,
    monitor_distributions,
    noise_tvd_limit,
    replace_luts,
)


HERE = Path(__file__).resolve().parent


def axis_rotation(axis, angle_degrees):
    axis = np.asarray(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    angle = math.radians(angle_degrees)
    cross = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    return (
        np.eye(3) * math.cos(angle)
        + (1.0 - math.cos(angle)) * np.outer(axis, axis)
        + math.sin(angle) * cross
    )


def ants_rotation(phi, theta, psi):
    """Return the active rotation used by ANTS3/TGeo for Phi,Theta,Psi."""
    rz = axis_rotation([0, 0, 1], phi)
    x_after_phi = rz @ np.array([1.0, 0.0, 0.0])
    rtheta = axis_rotation(x_after_phi, theta)
    z_after_theta = rtheta @ np.array([0.0, 0.0, 1.0])
    rpsi = axis_rotation(z_after_theta, psi)
    return rpsi @ rtheta @ rz


def geometry_object(config, name):
    matches = [
        obj for obj in config["Geometry"]["WorldTree"]
        if obj.get("Name") == name
    ]
    if len(matches) != 1:
        raise RuntimeError(
            "expected one geometry object named {}, found {}".format(
                name, len(matches)))
    return matches[0]


def rotate_setup(config, rotation_angles):
    phi, theta, psi = rotation_angles
    rotation = ants_rotation(phi, theta, psi)
    ex = rotation @ np.array([1.0, 0.0, 0.0])
    normal = rotation @ np.array([0.0, 0.0, 1.0])

    air_box = geometry_object(config, "AirBox")
    lyso = geometry_object(config, "LYSO")
    upper = geometry_object(
        config, "UpperMonitor_reflection_or_transmission")
    for obj in (lyso, upper):
        obj["Phi"], obj["Theta"], obj["Psi"] = rotation_angles

    # Keep the LYSO/air interface centered at the global origin.  The lower
    # monitor is a LYSO child and inherits this transform.  The upper monitor
    # is an AirBox child, so its transform is set explicitly.
    lyso_center = -5.0 * normal
    lyso["X"], lyso["Y"], lyso["Z"] = map(float, lyso_center)
    upper_center = 0.004 * normal
    upper["X"], upper["Y"], upper["Z"] = map(float, upper_center)
    # The original horizontal setup only needs dz=10 mm.  Increase the mother
    # volume in the rotated variant so every corner of the tilted LYSO remains
    # strictly contained and no unrelated parent-boundary overlap is introduced.
    air_box["ShapeSpecific"].update({"dx": 20.0, "dy": 20.0, "dz": 20.0})
    return rotation, ex, normal


def set_beam(config, angle, reverse, ex, normal):
    radians = math.radians(angle)
    depth = 0.002
    normal_sign = -1.0 if reverse else 1.0
    source = (
        -depth * math.tan(radians) * ex
        - normal_sign * depth * normal
    )
    direction = (
        math.sin(radians) * ex
        + normal_sign * math.cos(radians) * normal
    )
    config["PhotonSim"]["PhotonBombs"]["Single"]["Position"] = (
        source.astype(float).tolist())
    config["PhotonSim"]["PhGenOverrides"]["Direction"] = {
        "ConeAngle": 10,
        "DirectionMode": "Fixed",
        "DirectionVector": direction.astype(float).tolist(),
    }


def run_case(config, lsim, case_dir):
    case_dir.mkdir(parents=True, exist_ok=True)
    config["PhotonSim"]["Run"]["OutputDirectory"] = str(case_dir)
    config_path = case_dir / "config.json"
    json.dump(config, config_path.open("w"), indent=2)
    with (case_dir / "log").open("w") as log:
        completed = subprocess.run(
            [str(lsim), str(case_dir), config_path.name, "0"],
            cwd=case_dir, stdout=log, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError(
            "{}: lsim exited with {}".format(case_dir.name, completed.returncode))
    monitors = json.load(
        (case_dir / "PhotonMonitors.txt").open())["ParticleMonitorData"]
    if len(monitors) != 2:
        raise RuntimeError("expected exactly two photon monitors")
    return monitors


def normalized(measured, reverse):
    hits, theta, joint = measured
    theta = theta / hits if hits else theta
    joint = joint / hits if hits else joint
    if reverse:
        joint = joint[:, ::-1]
    return hits, theta, joint


def outcome_metrics(
        baseline, rotated, expected_probability, expected_joint, photons):
    baseline_hits, baseline_theta, baseline_joint = baseline
    rotated_hits, rotated_theta, rotated_joint = rotated
    expected_theta = expected_joint.sum(axis=1)

    probability_limit = max(
        0.001,
        6.0 * math.sqrt(max(
            expected_probability * (1.0 - expected_probability),
            1.0 / photons) / photons),
    )
    rotated_probability = rotated_hits / photons
    theta_tvd = 0.5 * np.abs(rotated_theta - expected_theta).sum()
    joint_tvd = 0.5 * np.abs(rotated_joint - expected_joint).sum()
    theta_limit = noise_tvd_limit(expected_theta, rotated_hits)
    joint_limit = noise_tvd_limit(expected_joint, rotated_hits)

    baseline_rotated_theta_tvd = (
        0.5 * np.abs(baseline_theta - rotated_theta).sum())
    baseline_rotated_joint_tvd = (
        0.5 * np.abs(baseline_joint - rotated_joint).sum())
    invariance_limit = 0.002

    passed = bool(
        abs(rotated_probability - expected_probability) <= probability_limit
        and theta_tvd <= theta_limit
        and joint_tvd <= joint_limit
        and abs(rotated_hits - baseline_hits) / photons <= invariance_limit
        and baseline_rotated_theta_tvd <= invariance_limit
        and baseline_rotated_joint_tvd <= invariance_limit
    )
    return {
        "baselineHits": baseline_hits,
        "rotatedHits": rotated_hits,
        "rotatedProbability": rotated_probability,
        "expectedProbability": expected_probability,
        "probabilityLimit": probability_limit,
        "thetaTvdToLut": float(theta_tvd),
        "thetaTvdLimit": theta_limit,
        "jointTvdToLut": float(joint_tvd),
        "jointTvdLimit": joint_limit,
        "baselineRotatedThetaTvd": float(baseline_rotated_theta_tvd),
        "baselineRotatedJointTvd": float(baseline_rotated_joint_tvd),
        "invarianceLimit": invariance_limit,
        "passed": passed,
    }


def plot_theta(results, expected, path, angle, rotation_angles):
    theta = np.arange(45) * 2.0 + 1.0
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for row, direction in enumerate(("fwd", "rev")):
        for col, outcome in enumerate(("reflected", "transmitted")):
            ax = axes[row, col]
            ax.step(
                theta, expected[(direction, outcome)].sum(axis=1),
                where="mid", color="black", linewidth=1.5, label="LUT")
            ax.step(
                theta, results[(direction, "baseline", outcome)][1],
                where="mid", color="tab:blue", alpha=0.8, label="baseline")
            ax.step(
                theta, results[(direction, "rotated", outcome)][1],
                where="mid", color="tab:orange", linestyle="--",
                label="rotated")
            ax.set_title("{} {}".format(direction, outcome))
            ax.set_ylabel("conditional probability")
            ax.grid(alpha=0.25)
    for ax in axes[-1]:
        ax.set_xlabel("theta_out in interface/monitor local frame (deg)")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(
        "Rotated DavisLUT geometry validation: incidence {} deg, "
        "Euler ({:g}, {:g}, {:g}) deg".format(angle, *rotation_angles))
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fwd", type=Path, default=DEFAULT_FWD)
    parser.add_argument("--rev", type=Path, default=DEFAULT_REV)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--lsim", type=Path, default=DEFAULT_LSIM)
    parser.add_argument("--photons", type=int, default=300000)
    parser.add_argument("--angle", type=float, default=45.0)
    parser.add_argument(
        "--rotation", type=float, nargs=3, metavar=("PHI", "THETA", "PSI"),
        default=[37.0, 29.0, 23.0])
    parser.add_argument(
        "--output-root", type=Path, default=HERE / "geom_out/rotated")
    parser.add_argument(
        "--metrics", type=Path,
        default=HERE / "data/geometry/rotated_metrics.json")
    parser.add_argument(
        "--figure", type=Path,
        default=HERE / "figures/validation_rotated.png")
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.photons < 1:
        raise SystemExit("--photons should be positive")
    if not 0 < args.angle < 90:
        raise SystemExit("--angle should be strictly between 0 and 90 degrees")
    if not args.lsim.is_file():
        raise SystemExit("lsim not found: {}".format(args.lsim))

    base = json.load(args.config.open())
    fwd = json.load(args.fwd.open())
    rev = json.load(args.rev.open())
    replace_luts(base, fwd, rev, args.fwd, args.rev)
    base["PhotonSim"]["PhotonBombs"]["PhotonsPerBomb"]["FixedNumber"] = (
        args.photons)
    base["PhotonSim"]["Run"]["SaveMonitors"] = True
    base["PhotonSim"]["Run"]["SaveSensorSignals"] = False
    base["PhotonSim"]["Run"]["SavePhotonBombs"] = False
    base["PhotonSim"]["Run"]["SaveStatistics"] = False

    identity = np.eye(3)
    rotated_config = copy.deepcopy(base)
    rotation, rotated_ex, rotated_normal = rotate_setup(
        rotated_config, args.rotation)
    frames = {
        "baseline": (base, identity[:, 0], identity[:, 2]),
        "rotated": (rotated_config, rotated_ex, rotated_normal),
    }

    results = {}
    expected = {}
    metrics = []
    summaries = []
    all_passed = True
    for direction, reverse in (("fwd", False), ("rev", True)):
        lut = fwd if direction == "fwd" else rev
        seed = 20260726 + 1000 * int(reverse)
        for variant, (template, ex, normal) in frames.items():
            config = copy.deepcopy(template)
            set_beam(config, args.angle, reverse, ex, normal)
            config["PhotonSim"]["Run"]["Seed"] = seed
            monitors = run_case(
                config, args.lsim,
                args.output_root / "{}_{}".format(variant, direction))
            reflected_index = 1 if reverse else 0
            transmitted_index = 0 if reverse else 1
            results[(direction, variant, "reflected")] = normalized(
                monitor_distributions(monitors[reflected_index]), reverse)
            results[(direction, variant, "transmitted")] = normalized(
                monitor_distributions(monitors[transmitted_index]), reverse)

        baseline_missing = args.photons - sum(
            results[(direction, "baseline", outcome)][0]
            for outcome in ("reflected", "transmitted"))
        rotated_missing = args.photons - sum(
            results[(direction, "rotated", outcome)][0]
            for outcome in ("reflected", "transmitted"))
        conservation_passed = bool(
            abs(baseline_missing) / args.photons <= 0.001
            and abs(rotated_missing) / args.photons <= 0.001)
        summaries.append({
            "direction": direction,
            "baselineMissing": baseline_missing,
            "rotatedMissing": rotated_missing,
            "passed": conservation_passed,
        })
        all_passed = all_passed and conservation_passed
        print(
            "{} conservation: baseline missing={}, rotated missing={} {}".format(
                direction, baseline_missing, rotated_missing,
                "PASS" if conservation_passed else "FAIL"),
            flush=True)

        for outcome, is_reflected in (
                ("reflected", True), ("transmitted", False)):
            probability, expected_joint = interpolated_expected(
                lut, args.angle, is_reflected)
            expected[(direction, outcome)] = expected_joint
            record = outcome_metrics(
                results[(direction, "baseline", outcome)],
                results[(direction, "rotated", outcome)],
                probability, expected_joint, args.photons)
            record.update({"direction": direction, "outcome": outcome})
            metrics.append(record)
            all_passed = all_passed and record["passed"]
            print(
                "{} {:11s}: P={:.6f}/{:.6f}, thetaTVD={:.5f}, "
                "jointTVD={:.5f}, rotation jointTVD={:.6f} {}".format(
                    direction, outcome,
                    record["rotatedProbability"],
                    record["expectedProbability"],
                    record["thetaTvdToLut"],
                    record["jointTvdToLut"],
                    record["baselineRotatedJointTvd"],
                    "PASS" if record["passed"] else "FAIL"),
                flush=True)

    args.metrics.parent.mkdir(parents=True, exist_ok=True)
    json.dump({
        "passed": bool(all_passed),
        "photonsPerRun": args.photons,
        "incidenceAngle": args.angle,
        "rotationAngles": args.rotation,
        "rotationMatrix": rotation.tolist(),
        "rotatedEx": rotated_ex.tolist(),
        "rotatedNormal": rotated_normal.tolist(),
        "summaries": summaries,
        "metrics": metrics,
    }, args.metrics.open("w"), indent=2)
    plot_theta(
        results, expected, args.figure, args.angle, args.rotation)
    print("Rotated geometry validation", "PASSED" if all_passed else "FAILED")
    raise SystemExit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
