#!/usr/bin/env python3
"""Measure periodic edge discontinuities in a leveled heightmap matrix.

Optionally write a C0-continuous mirror-tiled matrix for an A/B LUT-generation
diagnostic. The output retains the input pixel spacing and has shape
(2*Ny-1, 2*Nx-1).

Usage:
    python3 diagnose_heightmap_seam.py input_matrix.txt
    python3 diagnose_heightmap_seam.py input_matrix.txt --mirror-output mirrored.txt
"""

import argparse
import numpy as np


def rms(values):
    return float(np.sqrt(np.mean(np.asarray(values, dtype=float) ** 2)))


def report(z):
    dx = np.diff(z, axis=1)
    dy = np.diff(z, axis=0)
    seam_x = z[:, -1] - z[:, 0]
    seam_y = z[-1, :] - z[0, :]

    dx_rms = rms(dx)
    dy_rms = rms(dy)
    sx_rms = rms(seam_x)
    sy_rms = rms(seam_y)

    print("grid: {} x {}".format(z.shape[1], z.shape[0]))
    print("height range: {:.8g}".format(float(np.ptp(z))))
    print("internal x-step RMS: {:.8g}".format(dx_rms))
    print("periodic x-seam RMS: {:.8g}  ratio: {:.3f}".format(
        sx_rms, sx_rms / dx_rms if dx_rms else float("inf")))
    print("internal y-step RMS: {:.8g}".format(dy_rms))
    print("periodic y-seam RMS: {:.8g}  ratio: {:.3f}".format(
        sy_rms, sy_rms / dy_rms if dy_rms else float("inf")))


def mirror_tile(z):
    # Do not duplicate the outermost node: Nx nodes describe Nx-1 cells.
    mirrored_x = np.concatenate((z, z[:, -2::-1]), axis=1)
    return np.concatenate((mirrored_x, mirrored_x[-2::-1, :]), axis=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("heightmap", help="leveled height matrix")
    parser.add_argument("--mirror-output",
                        help="write a C0-continuous x/y mirror-tiled matrix")
    args = parser.parse_args()

    z = np.loadtxt(args.heightmap)
    if z.ndim != 2 or min(z.shape) < 2:
        raise SystemExit("heightmap should be a 2D matrix with at least 2 x 2 nodes")

    report(z)

    if args.mirror_output:
        mirrored = mirror_tile(z)
        np.savetxt(args.mirror_output, mirrored, fmt="%.8g")
        print("mirror output: {} ({} x {})".format(
            args.mirror_output, mirrored.shape[1], mirrored.shape[0]))
        print("mirror periodic seam max abs: x={:.3g}, y={:.3g}".format(
            float(np.max(np.abs(mirrored[:, -1] - mirrored[:, 0]))),
            float(np.max(np.abs(mirrored[-1, :] - mirrored[0, :])))))


if __name__ == "__main__":
    main()
