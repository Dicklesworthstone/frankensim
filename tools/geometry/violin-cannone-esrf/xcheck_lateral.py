"""xcheck_lateral.py XCHECK_TOP.tsv XCHECK_BACK.tsv : lateral registration error of the volume against the mesh, from the outer-surface height residual on the arched plates.
If the volume's pixel (x, z) is mapped with a translation error (tx, tz) and a scale error (ex, ez) (x_true = x_map (1 + ex) + tx ...), the height residual d = y_map - (offset + 0.03566 k) of a block is
d = a - gx tx - gz tz - gx ex x - gz ez (z - zc) to first order, with (gx, gz) the surface gradient of the volume (plane fit in the block). Regress d on [1, gx, gz, gx x, gz (z - zc)], robustly, per plate and jointly.
Prints the lateral translation (mm) and scale errors (percent) with the residual scatter before and after. Read-only."""
import sys

import numpy as np

VOX = 0.03566
ZC = -87.0


def load(p, top):
    rows = [l.rstrip("\n").split("\t") for l in open(p)]
    hdr = rows[0]
    a = np.array([[float(v) for v in r] for r in rows[1:]])
    c = {h: a[:, i] for i, h in enumerate(hdr)}
    k = c["k_start_med"] if top else c["k_end_med"]
    m = (c["frac_ok"] > 0.8) & (np.hypot(c["grad_x"], c["grad_z"]) < 0.6)
    d = (c["map_outer_y"] - VOX * k)
    return [v[m] for v in (c["x_mm"], c["z_mm"], d, c["grad_x"], c["grad_z"])]


def rsd(v):
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def fit(x, z, d, gx, gz, label):
    sel = np.abs(d - np.median(d)) < 0.8
    for it in range(8):
        A = np.c_[np.ones(len(d)), gx, gz, gx * x, gz * (z - ZC)]
        co = np.linalg.lstsq(A[sel], d[sel], rcond=None)[0]
        res = d - A @ co
        s = rsd(res[sel])
        sel = np.abs(res) < max(3 * s, 0.1)
    print("%s: %d blocks (%d kept); scatter of d about its median %.3f mm -> %.3f mm after the lateral model" % (label, len(d), sel.sum(), rsd(d[np.abs(d - np.median(d)) < 0.8]), s))
    print("   offset %.3f mm; translation: x %.3f mm, z %.3f mm; scale error: x %.2f %%, z %.2f %%  (signs: x_true = x_map (1 + e) + t, with the volume mapped at the nominal 0.28528 mm per pixel)" % (co[0], -co[1], -co[2], -100 * co[3], -100 * co[4]))
    return co


top = load(sys.argv[1], True)
back = load(sys.argv[2], False)
fit(*top, "top plate")
fit(*back, "back plate")
