"""xcheck_post.py VOL.npy SOUNDPOST_SECTIONS.tsv : independent check of the in-plane pixel pitch of the volume against the mesh: the soundpost cross-section in several frontal slices
(area = sum of grey / core grey over the disc, in pixels^2 x pitch^2) against the mesh soundpost section area at the same height (y = offset + 0.03566 k, offset from the plate fits). Read-only."""
import sys

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

VOX = 0.03566
S = 8 * VOX
MAP_X0, MAP_DX = -101.49264, 13.5
MAP_Z0, MAP_DZ = 88.69451, 824.5
vol = np.load(sys.argv[1], mmap_mode="r")
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[2]) if not l.startswith("frankensim.")]
tab = np.array([[float(v) for v in r] for r in rows[1:]])
OFF = -56.54  # mean of the two plate offsets from the lateral fit (-56.572 top, -56.513 back)
print("k\ty_mm\tvol_area_px\tvol_area_mm2_nominal\tmesh_area_mm2\tratio_mesh/volume\tvol_eq_diam_nominal\tmesh_eq_diam\tcentre_col\tcentre_row")
ratios = []
for k in range(2150, 2480, 20):
    a = np.asarray(vol[k]).astype(float)
    # the post lies near column 429, row 1537 (registered grid) and tilts little over this height: search a window of 40 pixels around it, locate the disc by the brightest 7 x 7 mean
    c0, r0 = 429, 1537
    win = a[r0 - 40:r0 + 40, c0 - 40:c0 + 40]
    loc = sliding_window_view(win, (7, 7)).mean(axis=(2, 3))
    iy, ix = np.unravel_index(loc.argmax(), loc.shape)
    cy, cx = iy + 3, ix + 3
    core = win[cy - 2:cy + 3, cx - 2:cx + 3]
    g0 = np.median(core)
    big = win[max(cy - 14, 0):cy + 15, max(cx - 14, 0):cx + 15]
    area_px = (big.clip(0, g0) / g0).sum()
    y = OFF + VOX * k
    ma = np.interp(y, tab[:, 0], tab[:, 1])
    va = area_px * S * S
    ratios.append(ma / va)
    print("%d\t%.2f\t%.1f\t%.2f\t%.2f\t%.4f\t%.3f\t%.3f\t%d\t%d" % (k, y, area_px, va, ma, ma / va, 2 * np.sqrt(va / np.pi), 2 * np.sqrt(ma / np.pi), c0 - 40 + cx, r0 - 40 + cy))
r = np.array(ratios)
print("mesh area / volume area (nominal pitch): median %.4f, mean %.4f, sd %.4f over %d slices -> implied linear ratio (mesh / volume nominal) %.4f" % (np.median(r), r.mean(), r.std(), len(r), np.sqrt(np.median(r))))
