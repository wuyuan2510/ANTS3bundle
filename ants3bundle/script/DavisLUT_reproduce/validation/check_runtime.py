#!/usr/bin/env python3
"""Report actual-rule validation status with a reliable process exit code."""
import json
import os
import sys


HERE = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(HERE, "data", "runtime")
all_passed = True

for direction in ("fwd", "rev"):
    path = os.path.join(RESULT_DIR, direction + ".json")
    try:
        with open(path, encoding="utf-8") as stream:
            result = json.load(stream)
    except (OSError, ValueError) as error:
        print("{}: missing or invalid result: {}".format(direction, error))
        all_passed = False
        continue

    passed = bool(result.get("passed", False))
    angles = result.get("angles", [])
    failed_angles = [record.get("angle") for record in angles if not record.get("passed", False)]
    covariance = result.get("rotatedNormalCovariance", {})
    absorption = result.get("syntheticAbsorption", {})
    passed = passed and bool(covariance.get("passed", False)) and bool(absorption.get("passed", False))
    all_passed = all_passed and passed
    print("{}: {}  angles={}  failed={}  norm_error={:.3g}  covariance_mismatches={}".format(
        direction, "PASS" if passed else "FAIL", len(angles), failed_angles,
        result.get("maximumNormError", float("nan")), covariance.get("mismatches", "?")))

print("Actual-rule validation", "PASSED" if all_passed else "FAILED")
sys.exit(0 if all_passed else 1)
