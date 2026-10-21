"""xcheck_fit.py XCHECK.tsv top|back : fit of the mesh outer-surface height against the slice index, per plate, on blocks away from the centre line (|x| > 30 mm: no bridge, strings, fingerboard or tailpiece above the top plate),
with the plane-fit gradient and the outer residual, for the slope (voxel pitch) test. Read-only."""
import sys

import numpy as np

VOX = 0.03566
top = sys.argv[2] == "top"
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
hdr = rows[0]
a = np.array([[float(v) for v in r] for r in rows[1:]])
c = {h: a[:, i] for i, h in enumerate(hdr)}
k = c["k_start_med"] if top else c["k_end_med"]
yo = c["map_outer_y"]
m0 = (np.abs(c["x_mm"]) > 30) & (c["frac_ok"] > 0.8)
print("blocks with |x| > 30 mm and 80 %% columns agreeing: %d" % m0.sum())
sel = m0.copy()
for it in range(8):
    A = np.c_[np.ones(sel.sum()), k[sel]]
    (y0, s), *_ = np.linalg.lstsq(A, yo[sel], rcond=None)
    res = yo - (y0 + s * k)
    sd = 1.4826 * np.median(np.abs(res[sel] - np.median(res[sel])))
    sel = m0 & (np.abs(res) < max(3 * sd, 0.1))
print("fit y_map = %.4f + %.5f k on %d blocks; slope / %.5f = %.4f; robust scatter %.3f mm" % (y0, s, sel.sum(), VOX, s / VOX, sd))
print("k range of the fit %d..%d (%.1f mm); map outer y range %.2f..%.2f mm" % (k[sel].min(), k[sel].max(), (k[sel].max() - k[sel].min()) * VOX, yo[sel].min(), yo[sel].max()))
# slope dependence on the local gradient: split by the surface steepness
g = np.hypot(c["grad_x"], c["grad_z"])
for lo, hi in ((0, 0.1), (0.1, 0.25), (0.25, 0.5), (0.5, 5)):
    mm = sel & (g >= lo) & (g < hi)
    if mm.sum() > 20:
        print("  gradient %.2f..%.2f: %d blocks, median residual (map y - fitted from k with the voxel slope and the fitted offset) %.3f mm" % (lo, hi, mm.sum(), np.median((yo - (np.median((yo - VOX * k)[sel]) + VOX * k))[mm])))
# the fitted line against k in bins
edges = np.percentile(k[sel], np.linspace(0, 100, 7))
for i in range(6):
    mm = sel & (k >= edges[i]) & (k <= edges[i + 1])
    print("  k %5d..%5d: %d blocks, mean map y %.3f, mean y from k (voxel slope, offset %.3f) %.3f, difference %.3f mm" % (edges[i], edges[i + 1], mm.sum(), yo[mm].mean(), np.median((yo - VOX * k)[sel]), np.median((yo - VOX * k)[sel]) + VOX * k[mm].mean(), yo[mm].mean() - (np.median((yo - VOX * k)[sel]) + VOX * k[mm].mean())))
