"""xcheck_tilt.py XCHECK_TOP.tsv XCHECK_BACK.tsv : is the offset between the mesh y and the slice index (at the voxel pitch) a plane in (x, z)? A small rotation between the mesh frame and the volume frame
would show as the same plane coefficients on both plates. Fits d = y_map - 0.03566 * k = a + b x + c z (outer surface), robustly, per plate and on both together with one tilt and one offset per plate.
Prints coefficients (mm per mm, i.e. radians), robust scatter before and after. Read-only."""
import sys

import numpy as np

VOX = 0.03566


def load(p, top):
    rows = [l.rstrip("\n").split("\t") for l in open(p)]
    hdr = rows[0]
    a = np.array([[float(v) for v in r] for r in rows[1:]])
    c = {h: a[:, i] for i, h in enumerate(hdr)}
    k = c["k_start_med"] if top else c["k_end_med"]
    m = c["frac_ok"] > 0.8
    return c["x_mm"][m], c["z_mm"][m], (c["map_outer_y"] - VOX * k)[m], c["len_med_slices"][m] * VOX, np.abs(c["map_inner_y"] - c["map_outer_y"])[m]


def rsd(v):
    return 1.4826 * np.median(np.abs(v - np.median(v)))


data = {}
for name, p, top in (("top", sys.argv[1], True), ("back", sys.argv[2], False)):
    x, z, d, vl, mv = load(p, top)
    sel = np.abs(d - np.median(d)) < 0.6
    for it in range(8):
        A = np.c_[np.ones(sel.sum()), x[sel], z[sel]]
        co = np.linalg.lstsq(A, d[sel], rcond=None)[0]
        res = d - np.c_[np.ones(len(d)), x, z] @ co
        s = rsd(res[sel])
        sel = np.abs(res) < max(3 * s, 0.1)
    print("%s plate: %d blocks; offset-only robust scatter %.3f mm; plane d = %.3f + %.5f x + %.5f z (mm, mm/mm); scatter after %.3f mm (%d blocks)" % (
        name, len(d), rsd(d[np.abs(d - np.median(d)) < 0.6]), co[0], co[1], co[2], s, sel.sum()))
    data[name] = (x, z, d, sel, co)
# joint fit: one tilt, two offsets
xs = np.concatenate([data["top"][0][data["top"][3]], data["back"][0][data["back"][3]]])
zs = np.concatenate([data["top"][1][data["top"][3]], data["back"][1][data["back"][3]]])
ds = np.concatenate([data["top"][2][data["top"][3]], data["back"][2][data["back"][3]]])
is_back = np.r_[np.zeros(data["top"][3].sum()), np.ones(data["back"][3].sum())]
A = np.c_[1 - is_back, is_back, xs, zs]
co = np.linalg.lstsq(A, ds, rcond=None)[0]
res = ds - A @ co
print("joint: offsets top %.3f back %.3f mm, common tilt %.5f per mm in x, %.5f per mm in z; scatter %.3f mm; the offset difference back - top is %.3f mm over the %.1f mm between the plates' mean slice indices" % (
    co[0], co[1], co[2], co[3], rsd(res), co[1] - co[0], 0.0))
