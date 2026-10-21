"""xcheck_stats.py XCHECK.tsv top|back : statistics of the volume-versus-mesh cross-check written by xcheck_plates.py.
(1) offset and slope between the slice index and the mesh y (outer surface), robustly fitted; the slope is a test of the voxel size (0.03566 mm) and of the mesh units;
(2) vertical thickness of the volume (median run length) against the maps-table inner_y - outer_y, on blocks whose run is not contaminated (bridge, bass bar, soundpost, linings) and whose columns agree;
(3) the same for the normal thickness with the plane-fit gradient of the volume and the table's own normal thickness; (4) the contaminated share and where. Read-only."""
import sys

import numpy as np

VOX = 0.03566
top = sys.argv[2] == "top"
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
hdr = rows[0]
a = np.array([[float(v) for v in r] for r in rows[1:]])
c = {h: a[:, i] for i, h in enumerate(hdr)}
x, z = c["x_mm"], c["z_mm"]
my_out, my_in, mt_norm = c["map_outer_y"], c["map_inner_y"], c["map_thick_normal"]
k_ref = c["k_start_med"] if top else c["k_end_med"]
k_oth = c["k_end_med"] if top else c["k_start_med"]
vlen = c["len_med_slices"] * VOX
vmap = np.abs(my_in - my_out)
print("%d blocks" % len(a))

# (1) outer surface: y_map = y0 + s * k
sel = np.ones(len(a), bool)
for it in range(6):
    A = np.c_[np.ones(sel.sum()), k_ref[sel]]
    (y0, s), *_ = np.linalg.lstsq(A, my_out[sel], rcond=None)
    res = my_out - (y0 + s * k_ref)
    sd = 1.4826 * np.median(np.abs(res[sel] - np.median(res[sel])))
    sel = np.abs(res) < max(3 * sd, 0.15)
print("outer surface: y_map = %.4f + %.5f * k over %d blocks (%.1f %%); slope / voxel pitch = %.4f; robust scatter %.3f mm; offset-removed residual median %.3f mm (p5 %.3f p95 %.3f)" % (
    y0, s, sel.sum(), 100 * sel.mean(), s / (VOX if top else VOX), sd, np.median(res[sel]), np.percentile(res[sel], 5), np.percentile(res[sel], 95)))
# fix the slope to the voxel pitch and refit the offset
y0f = np.median((my_out - VOX * k_ref)[sel])
res_f = my_out - (y0f + VOX * k_ref)
print("with the slope fixed to the voxel pitch: offset %.3f mm; residual median %.3f, p5 %.3f, p95 %.3f, robust sd %.3f mm" % (
    y0f, np.median(res_f[sel]), np.percentile(res_f[sel], 5), np.percentile(res_f[sel], 95), 1.4826 * np.median(np.abs(res_f[sel] - np.median(res_f[sel])))))

# (2) vertical thickness
dv = vlen - vmap
good = (np.abs(dv) < 1.0) & (c["frac_ok"] > 0.7) & sel
print("\nvertical thickness (volume run length minus table |inner_y - outer_y|):")
for name, m in (("all blocks", np.ones(len(a), bool)), ("agreeing (|d| < 1 mm, 70 %% of columns within 0.5 mm of the median, outer surface within the fit)", good)):
    d = dv[m]
    print("  %-90s n=%d: median %.3f mm, mean %.3f, sd %.3f, p5 %.3f, p95 %.3f; table median %.3f mm, volume median %.3f mm" % (name, m.sum(), np.median(d), d.mean(), d.std(), np.percentile(d, 5), np.percentile(d, 95), np.median(vmap[m]), np.median(vlen[m])))
print("  share of blocks where the volume run is more than 1 mm longer (bar, post, linings, bridge, repairs): %.1f %%; more than 1 mm shorter: %.1f %%; columns disagreeing: %.1f %% of blocks" % (
    100 * np.mean(dv > 1.0), 100 * np.mean(dv < -1.0), 100 * np.mean(c["frac_ok"] <= 0.7)))

# (3) normal thickness
gnorm = np.sqrt(1 + c["grad_x"] ** 2 + c["grad_z"] ** 2)
vnorm = vlen / gnorm
dn = vnorm - mt_norm
print("\nnormal thickness (volume run / sqrt(1 + gradient^2) minus table thickness_normal_mm):")
for name, m in (("all blocks", np.ones(len(a), bool)), ("agreeing blocks", good)):
    d = dn[m]
    print("  %-20s n=%d: median %.3f mm, mean %.3f, sd %.3f, p5 %.3f, p95 %.3f; table median %.3f, volume median %.3f" % (name, m.sum(), np.median(d), d.mean(), d.std(), np.percentile(d, 5), np.percentile(d, 95), np.median(mt_norm[m]), np.median(vnorm[m])))
rel = dn[good] / mt_norm[good]
print("  relative difference on agreeing blocks: median %.2f %%, p5 %.2f %%, p95 %.2f %%" % (100 * np.median(rel), 100 * np.percentile(rel, 5), 100 * np.percentile(rel, 95)))
print("  correlation of the block thickness maps (agreeing blocks): %.4f" % np.corrcoef(vnorm[good], mt_norm[good])[0, 1])

# (4) a few sections
def section(label, mask, key, fmt_key):
    idx = np.nonzero(mask)[0]
    idx = idx[np.argsort(key[idx])]
    print("\n%s (%d blocks; every 4th shown)" % (label, len(idx)))
    print("  %s\tmap_vert\tvol_vert\tmap_norm\tvol_norm\tmap_outer\tvol_outer(offset fitted)" % fmt_key)
    for j in idx[::4]:
        print("  %.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f%s" % (key[j], vmap[j], vlen[j], mt_norm[j], vnorm[j], my_out[j], y0f + VOX * k_ref[j], "" if good[j] else "   (not agreeing)"))
    m = idx
    print("  section median differences: vertical %.3f mm, normal %.3f mm (all shown blocks, including not-agreeing ones: %d)" % (np.median(dv[m]), np.median(dn[m]), (~good[m]).sum()))


def near(arr, t):
    u = np.unique(arr)
    return u[np.argmin(np.abs(u - t))]


zt = near(z, -100.0)
section("transverse section z = %.1f mm" % zt, np.abs(z - zt) < 0.1, x, "x_mm")
xt = near(x, 0.0)
section("longitudinal section x = %.1f mm (centre line)" % xt, np.abs(x - xt) < 0.1, z, "z_mm")
xt = near(x, -40.0 if top else 40.0)
section("longitudinal section x = %.1f mm" % xt, np.abs(x - xt) < 0.1, z, "z_mm")
