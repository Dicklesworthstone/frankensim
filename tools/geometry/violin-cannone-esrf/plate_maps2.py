"""plate_maps2.py GLB OUT.npz OUTER_SIGN : like plate_maps.py but rasterises the SURFACE, not only its vertices: every triangle is sampled on a barycentric lattice whose spacing
is at most STEP mm (so a triangle larger than a cell still reaches every (x, z) cell it covers), and the extreme Y per cell is kept. Frame, cell size and saved fields as in
plate_maps.py (outer = smallest Y for OUTER_SIGN -1, largest for +1). Adds `samples` (number of surface samples per cell)."""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cannone_io as io  # noqa: E402

CELL, STEP = 0.5, 0.2
path, out, sign = sys.argv[1], sys.argv[2], int(sys.argv[3])
t0 = time.time()
v, t, info = io.glb_triangles(path)
print("loaded", len(v), "vertices", len(t), "triangles in %.0f s" % (time.time() - t0), flush=True)
a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
lmax = np.maximum(np.maximum(np.linalg.norm(b - a, axis=1), np.linalg.norm(c - b, axis=1)), np.linalg.norm(a - c, axis=1))
n_of = np.maximum(1, np.ceil(lmax / STEP)).astype(np.int64)
print("triangles needing subdivision: %d of %d, largest edge %.2f mm, max n %d" % (int((n_of > 1).sum()), len(t), lmax.max(), n_of.max()), flush=True)
lo_all = np.minimum(np.minimum(a, b), c).min(0)
hi_all = np.maximum(np.maximum(a, b), c).max(0)
x0, z0 = np.floor(lo_all[0]) - 1.0, np.floor(lo_all[2]) - 1.0
nx = int((hi_all[0] - x0) / CELL) + 2
nz = int((hi_all[2] - z0) / CELL) + 2
lo = np.full(nx * nz, np.inf)
hi = np.full(nx * nz, -np.inf)
cnt = np.zeros(nx * nz, np.int64)


def bin_points(p):
    ix = ((p[:, 0] - x0) / CELL).astype(np.int64)
    iz = ((p[:, 2] - z0) / CELL).astype(np.int64)
    idx = ix * nz + iz
    np.minimum.at(lo, idx, p[:, 1])
    np.maximum.at(hi, idx, p[:, 1])
    cnt[:] += np.bincount(idx, minlength=nx * nz)


nsamp = 0
for n in np.unique(n_of):
    sel = np.nonzero(n_of == n)[0]
    if n == 1:
        for tri_part in (a[sel], b[sel], c[sel], (a[sel] + b[sel] + c[sel]) / 3.0):
            bin_points(tri_part)
            nsamp += len(sel)
        continue
    ii, jj = np.meshgrid(np.arange(n + 1), np.arange(n + 1), indexing="ij")
    keep = (ii + jj) <= n
    wi, wj = ii[keep] / n, jj[keep] / n
    for s in range(0, len(sel), 200000):
        part = sel[s:s + 200000]
        pa, pb, pc = a[part], b[part], c[part]
        pts = pa[:, None, :] + (pb - pa)[:, None, :] * wi[None, :, None] + (pc - pa)[:, None, :] * wj[None, :, None]
        pts = pts.reshape(-1, 3)
        bin_points(pts)
        nsamp += len(pts)
print("samples binned: %d in %.0f s" % (nsamp, time.time() - t0), flush=True)
empty = cnt == 0
lo[empty] = np.nan
hi[empty] = np.nan
outer, inner = (lo, hi) if sign < 0 else (hi, lo)
np.savez_compressed(out, x0=x0, z0=z0, cell=CELL, nx=nx, nz=nz, outer=outer.reshape(nx, nz), inner=inner.reshape(nx, nz), count=cnt.reshape(nx, nz),
                    bbox_lo=lo_all, bbox_hi=hi_all, source=os.path.basename(path), step=STEP)
print("saved", out, "grid", nx, "x", nz, "occupied cells", int((~empty).sum()), "in %.0f s" % (time.time() - t0), flush=True)
