#!/usr/bin/env python3
"""Pool several surface .lut files (same binning, same n1/n2) into one ensemble LUT.

Summing the integer count arrays is equivalent to ray-tracing photons over all the surface
patches combined, so it gives a properly statistics-weighted ensemble-average LUT.

Usage:
    python3 pool_luts.py pooled.lut  loc1.lut loc2.lut loc3.lut ...

Each per-location .lut is produced with the ANTS3 script unit, e.g.
    rules.generateSurfaceLut("loc1_leveled.txt", "loc1.lut",
        {n1:1.824, n2:1.0, pixelSizeX:0.0195694, pixelSizeY:0.0195694, format:"matrix"});
(and likewise with swapped n1/n2 for the reverse-direction LUTs).
"""
import datetime
import json
import sys

import numpy as np


def pool(out_path, in_paths):
    with open(in_paths[0]) as stream:
        base = json.load(stream)
    count_keys = ["Launched", "ReflectedCounts", "TransmittedCounts", "AbsorbedCounts"]
    acc = {k: np.asarray(base[k], dtype=np.int64) for k in count_keys}
    RH = np.asarray(base["ReflectedHist"], dtype=np.int64)
    TH = np.asarray(base["TransmittedHist"], dtype=np.int64)
    launch_weights = [int(acc["Launched"].sum())]
    mean_bounces = [float(base["Meta"].get("MeanBounces", 0.0))]
    anomaly_count = int(base["Meta"].get("AnomalyCount", 0))
    compatibility = {
        "FormatVersion": base.get("FormatVersion"),
        "Binning": base["Binning"],
        "n1": base["Meta"].get("n1"),
        "n2": base["Meta"].get("n2"),
        "Wavelength": base["Meta"].get("Wavelength"),
        "PixelSize": base["Meta"].get("PixelSize"),
    }

    for p in in_paths[1:]:
        with open(p) as stream:
            d = json.load(stream)
        candidate = {
            "FormatVersion": d.get("FormatVersion"),
            "Binning": d["Binning"],
            "n1": d["Meta"].get("n1"),
            "n2": d["Meta"].get("n2"),
            "Wavelength": d["Meta"].get("Wavelength"),
            "PixelSize": d["Meta"].get("PixelSize"),
        }
        assert candidate == compatibility, "incompatible LUT metadata in {}".format(p)
        for k in count_keys:
            acc[k] += np.asarray(d[k], dtype=np.int64)
        RH += np.asarray(d["ReflectedHist"], dtype=np.int64)
        TH += np.asarray(d["TransmittedHist"], dtype=np.int64)
        launch_weights.append(int(np.asarray(d["Launched"], dtype=np.int64).sum()))
        mean_bounces.append(float(d["Meta"].get("MeanBounces", 0.0)))
        anomaly_count += int(d["Meta"].get("AnomalyCount", 0))

    out = dict(base)
    for k in count_keys:
        out[k] = acc[k].astype(int).tolist()
    out["ReflectedHist"]   = [[int(v) for v in row] for row in RH]
    out["TransmittedHist"] = [[int(v) for v in row] for row in TH]
    out["Meta"] = dict(base["Meta"])
    out["Meta"]["MeanBounces"] = float(np.average(mean_bounces, weights=launch_weights))
    out["Meta"]["AnomalyCount"] = anomaly_count
    if np.all(acc["Launched"] == acc["Launched"][0]):
        out["Meta"]["PhotonsPerThetaBin"] = int(acc["Launched"][0])
    out["Meta"]["SourceHeightmap"] = "mirror-periodic ensemble of {} AFM scans".format(
        len(in_paths))
    out["Meta"]["Comment"] = (
        "pooled from {} independent mirror-periodic LUTs (summed integer counts)"
        .format(len(in_paths)))
    out["Meta"]["GenerationDate"] = (
        datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z")
    out["Meta"]["Seed"] = 0  # the inputs use independent seeds
    with open(out_path, "w") as stream:
        json.dump(out, stream)

    R = acc["ReflectedCounts"] / acc["Launched"]
    n = base["Binning"]["ThetaIncBins"]
    print("pooled {} LUTs -> {}".format(len(in_paths), out_path))
    print("  R(0deg)={:.3f}  R(45deg)={:.3f}  R(grazing)={:.3f}"
          "  ~{} photons/bin".format(
              R[0], R[n//2], R[-1], int(acc["Launched"][0])))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    pool(sys.argv[1], sys.argv[2:])
