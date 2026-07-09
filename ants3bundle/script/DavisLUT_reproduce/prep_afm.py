#!/usr/bin/env python3
"""Convert a raw AFM export (*.spm.txt) into a leveled height matrix for the LUT generator.

Extracts one height column, reshapes to the scan grid (row-major), converts nm -> um, and
subtracts a best-fit plane (removes sample tilt) -- matching AFM_surface_microfacet_analysis.

Usage:
    python3 prep_afm.py <input.spm.txt> <output_leveled_um.txt> [column=1] [scan_um=10] [N=512]

column: 1-based; 1 = height trace, 4 = height retrace (for these Bruker exports).
The output feeds rules.generateSurfaceLut(..., {format:"matrix", pixelSizeX/Y: scan_um/(N-1)}).
"""
import sys, numpy as np

def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    infile, outfile = sys.argv[1], sys.argv[2]
    col     = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    scan_um = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0
    N       = int(sys.argv[5]) if len(sys.argv) > 5 else 512

    vals = []
    for line in open(infile, encoding="latin1"):
        tok = line.split()
        if len(tok) < col: continue
        try: vals.append(float(tok[col-1]))
        except ValueError: continue           # header / non-numeric line
    h = np.array(vals[:N*N], dtype=float)
    assert h.size >= N*N, f"got {h.size} values, need {N*N}"
    Z = h.reshape(N, N) / 1000.0               # nm -> um

    x = np.linspace(0, scan_um, N); X, Y = np.meshgrid(x, x)
    A = np.column_stack([X.ravel(), Y.ravel(), np.ones(X.size)])
    a, b, c = np.linalg.lstsq(A, Z.ravel(), rcond=None)[0]
    Zd = Z - (a*X + b*Y + c)

    np.savetxt(outfile, Zd, fmt="%.6e")
    Sq = np.sqrt(np.mean((Zd - Zd.mean())**2)) * 1000
    px = scan_um / (N - 1)
    print(f"wrote {outfile}: {N}x{N}, pixel = {px*1000:.3f} nm, Sq = {Sq:.1f} nm")
    print(f"  -> generateSurfaceLut(..., {{format:'matrix', pixelSizeX:{px:.7f}, pixelSizeY:{px:.7f}}})")

if __name__ == "__main__":
    main()
