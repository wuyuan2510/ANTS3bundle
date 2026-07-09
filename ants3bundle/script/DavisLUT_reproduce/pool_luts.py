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
import json, sys
import numpy as np

def pool(out_path, in_paths):
    base = json.load(open(in_paths[0]))
    count_keys = ["Launched", "ReflectedCounts", "TransmittedCounts", "AbsorbedCounts"]
    acc = {k: np.array(base[k], float) for k in count_keys}
    RH = np.array(base["ReflectedHist"], float)
    TH = np.array(base["TransmittedHist"], float)
    mean_bounces = [base["Meta"].get("MeanBounces", 0.0)]

    for p in in_paths[1:]:
        d = json.load(open(p))
        assert d["Binning"] == base["Binning"], f"binning mismatch in {p}"
        for k in count_keys:
            acc[k] += np.array(d[k], float)
        RH += np.array(d["ReflectedHist"], float)
        TH += np.array(d["TransmittedHist"], float)
        mean_bounces.append(d["Meta"].get("MeanBounces", 0.0))

    out = dict(base)
    for k in count_keys:
        out[k] = [int(v) for v in acc[k]]
    out["ReflectedHist"]   = [[int(v) for v in row] for row in RH]
    out["TransmittedHist"] = [[int(v) for v in row] for row in TH]
    out["Meta"] = dict(base["Meta"])
    out["Meta"]["MeanBounces"] = float(np.mean(mean_bounces))
    out["Meta"]["SourceHeightmap"] = f"ENSEMBLE of {len(in_paths)} scans"
    out["Meta"]["Comment"] = f"pooled from {len(in_paths)} LUTs (summed counts)"
    json.dump(out, open(out_path, "w"))

    R = acc["ReflectedCounts"] / acc["Launched"]
    n = base["Binning"]["ThetaIncBins"]
    print(f"pooled {len(in_paths)} LUTs -> {out_path}")
    print(f"  R(0deg)={R[0]:.3f}  R(45deg)={R[n//2]:.3f}  R(grazing)={R[-1]:.3f}  ~{acc['Launched'][0]:.0f} photons/bin")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    pool(sys.argv[1], sys.argv[2:])
