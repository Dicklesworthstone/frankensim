"""cannone_sections.py : helpers shared by the il Cannone table builders (numpy only).
  sections(v, t, axis, step, lo, hi) -> rows (s, area, c_u0, c_u1, u0_min, u0_max, u1_min, u1_max) by planes perpendicular to `axis` (Green's theorem on the oriented segments)
  surface_samples(v, t, step) -> iterator of (k, 3) arrays of points covering every triangle with lattice spacing <= step
  silhouette_runs(points, axes, cell) -> {row index: [(start, end), ...]} of occupied cells, for the projection onto the two chosen axes
"""
import numpy as np


def sections(v, t, axis, step, lo=None, hi=None):
    others = [i for i in range(3) if i != axis]
    a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    lo_t = np.minimum(np.minimum(a[:, axis], b[:, axis]), c[:, axis])
    hi_t = np.maximum(np.maximum(a[:, axis], b[:, axis]), c[:, axis])
    lo = lo_t.min() + step / 2 if lo is None else lo
    hi = hi_t.max() if hi is None else hi
    rows = []
    for s0 in np.arange(lo, hi, step):
        s = s0 + 1.2345e-5  # keep every plane off the vertices of the mesh
        cand = np.nonzero((lo_t <= s) & (hi_t > s))[0]
        if len(cand) == 0:
            continue
        pa, pb, pc = a[cand], b[cand], c[cand]
        da, db, dc = pa[:, axis] - s, pb[:, axis] - s, pc[:, axis] - s
        P, Q = [], []
        for (p, q, r, dp, dq, dr) in ((pa, pb, pc, da, db, dc), (pb, pc, pa, db, dc, da), (pc, pa, pb, dc, da, db)):
            alone = ((dp > 0) != (dq > 0)) & ((dp > 0) != (dr > 0))
            if not alone.any():
                continue
            pp, qq, rr, d1, d2, d3 = p[alone], q[alone], r[alone], dp[alone], dq[alone], dr[alone]
            e1 = pp + (qq - pp) * (d1 / (d1 - d2))[:, None]
            e2 = pp + (rr - pp) * (d1 / (d1 - d3))[:, None]
            pos = (d1 > 0)[:, None]
            P.append(np.where(pos, e1, e2))
            Q.append(np.where(pos, e2, e1))
        if not P:
            continue
        P, Q = np.concatenate(P), np.concatenate(Q)
        u0, u1 = others
        cross = P[:, u0] * Q[:, u1] - Q[:, u0] * P[:, u1]
        area = 0.5 * cross.sum()
        if abs(area) < 1e-9:
            continue
        cx = ((P[:, u0] + Q[:, u0]) * cross).sum() / (6 * area)
        cy = ((P[:, u1] + Q[:, u1]) * cross).sum() / (6 * area)
        allp = np.concatenate([P, Q])
        rows.append((float(s0), abs(float(area)), float(cx), float(cy), float(allp[:, u0].min()), float(allp[:, u0].max()), float(allp[:, u1].min()), float(allp[:, u1].max())))
    return rows


def surface_samples(v, t, step):
    a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    lmax = np.maximum(np.maximum(np.linalg.norm(b - a, axis=1), np.linalg.norm(c - b, axis=1)), np.linalg.norm(a - c, axis=1))
    n_of = np.maximum(1, np.ceil(lmax / step)).astype(np.int64)
    for n in np.unique(n_of):
        sel = np.nonzero(n_of == n)[0]
        if n == 1:
            for part in (a[sel], b[sel], c[sel], (a[sel] + b[sel] + c[sel]) / 3.0):
                yield part
            continue
        ii, jj = np.meshgrid(np.arange(n + 1), np.arange(n + 1), indexing="ij")
        keep = (ii + jj) <= n
        wi, wj = ii[keep] / n, jj[keep] / n
        for s in range(0, len(sel), 100000):
            part = sel[s:s + 100000]
            pa, pb, pc = a[part], b[part], c[part]
            yield (pa[:, None, :] + (pb - pa)[:, None, :] * wi[None, :, None] + (pc - pa)[:, None, :] * wj[None, :, None]).reshape(-1, 3)


def silhouette_runs(v, t, axes, cell, step):
    """Occupancy of the projection of the surface onto axes (u0, u1): returns (origin_u0, origin_u1, boolean array [n_u0, n_u1])."""
    u0, u1 = axes
    lo = np.array([v[:, u0].min(), v[:, u1].min()]) - cell
    hi = np.array([v[:, u0].max(), v[:, u1].max()]) + cell
    n0 = int((hi[0] - lo[0]) / cell) + 1
    n1 = int((hi[1] - lo[1]) / cell) + 1
    occ = np.zeros((n0, n1), bool)
    for pts in surface_samples(v, t, step):
        i0 = ((pts[:, u0] - lo[0]) / cell).astype(np.int64)
        i1 = ((pts[:, u1] - lo[1]) / cell).astype(np.int64)
        occ[i0, i1] = True
    return lo[0], lo[1], occ


def runs_of(row):
    """start/end indices (inclusive) of True runs in a 1-D boolean array"""
    d = np.diff(np.concatenate([[0], row.astype(np.int8), [0]]))
    return list(zip(np.nonzero(d == 1)[0].tolist(), (np.nonzero(d == -1)[0] - 1).tolist()))


def plane_segments(a, b, c, axis, s):
    """oriented intersection segments (P, Q) of the triangles (a, b, c) with the plane coordinate[axis] = s; consistent orientation, for Green's theorem"""
    da, db, dc = a[:, axis] - s, b[:, axis] - s, c[:, axis] - s
    P, Q = [], []
    for (p, q, r, dp, dq, dr) in ((a, b, c, da, db, dc), (b, c, a, db, dc, da), (c, a, b, dc, da, db)):
        alone = ((dp > 0) != (dq > 0)) & ((dp > 0) != (dr > 0))
        if not alone.any():
            continue
        pp, qq, rr, d1, d2, d3 = p[alone], q[alone], r[alone], dp[alone], dq[alone], dr[alone]
        e1 = pp + (qq - pp) * (d1 / (d1 - d2))[:, None]
        e2 = pp + (rr - pp) * (d1 / (d1 - d3))[:, None]
        pos = (d1 > 0)[:, None]
        P.append(np.where(pos, e1, e2))
        Q.append(np.where(pos, e2, e1))
    if not P:
        return np.zeros((0, 3)), np.zeros((0, 3))
    return np.concatenate(P), np.concatenate(Q)


def green_area(P, Q, u0, u1):
    cross = P[:, u0] * Q[:, u1] - Q[:, u0] * P[:, u1]
    return 0.5 * float(cross.sum())
