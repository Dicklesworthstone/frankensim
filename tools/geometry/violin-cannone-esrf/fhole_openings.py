"""fhole_openings.py TABLE_ALL_GLB OUT_DIR : the il Cannone f-hole through-openings and walls from the top-plate mesh
(2416550405_table_all.glb of the ESRF data collection DOI 10.15151/ESRF-DC-2413577471, CC BY 4.0), as data tables of the
repo description (TSV, mm, the frame of the source files: x across the body, y vertical with the top plate's outer
surface at the smaller y, z along the body). Numpy only; reads nothing but the GLB (refused unless its SHA-256 starts with
the prefix specimen.tsv records).

Vertical rays (along y) through the centre of every cell of a plan (x, z) grid; the sorted crossings of a ray bound the
wood. A cell is OPEN THROUGH THE PLATE when its ray crosses no surface.
1. Every 0.5 mm cell of the plate's bounding box: the open cells not connected to the box edge are the openings inside the
   plate. Openings of 100 mm2 or more are the two f-holes (hole 1 at negative x, the bass-bar side); every smaller opening
   within 10 mm of one joins it (an eye whose neck the 0.5 mm rays miss), and the rest are listed, not used.
2. Every 0.1 mm cell of each hole's box (its 0.5 mm cells plus 3 mm): the hole is the union of the open components of this
   grid that meet its 0.5 mm cells and stay off the box edge.
3. The 0.5 mm raster (top-plate-fholes.tsv, the grid of the specimen's other rasters, cell edges on multiples of 0.5 mm):
   a cell is open when at least 13 of its 25 cells of 0.1 mm are open (at least half). The rule it replaces counted a
   0.5 mm cell as wood when any surface sample fell in it ("edge as wood"); that rule is evaluated here on the same 0.1 mm
   data (a cell open only when all 25 of its cells are) and its areas are written beside the new ones.
4. Walls: the outer and inner surfaces fitted by cubics in (x, z) on cells 1.5 to 5 mm from the hole that the ray crosses
   once through the wood (robust: three passes dropping residuals above 4 sigma); the opening area A(d) at depth fraction
   d of the thickness in two layer families, VERTICAL (planes of constant fitted depth along y) and NORMAL (surfaces
   parallel to the outer surface, offset along its normal); the through-resistance R = integral of t dd / A(d) of each
   family; wall sections at stations along the hole's principal axis. The family whose walls run straighter (the smaller
   mean drift of the opening's mid-point from d = 0.1 to d = 0.9 over the stations |s| <= 0.25 L) is the hole's wall
   direction, and its R is the hole's through-resistance.
   The hole's straight bore is t_k = R_k x A_k (A_k the 0.5 mm raster's area): a straight duct of the raster's outline
   and length t_k has the measured through-resistance. One bore for both holes keeps their parallel bore conductance:
   t_pair = (A_1 + A_2) / (1 / R_1 + 1 / R_2).
Writes OUT_DIR/top-plate-fholes.tsv, top-plate-fholes-0p1mm.tsv (the 0.1 mm runs it is coarsened from) and
top-plate-fhole-walls.tsv (one quantity per row), each under a comment line naming the source, as the lane's other
builders do (assemble.py adds the schema lines). Prints a summary.
Run on the fetched mesh: python3 -I fhole_openings.py RAW/stage2/2416550405_table_all.glb OUT_DIR (about 3 GB of memory)."""
import hashlib
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cannone_io as io  # noqa: E402

GLB_SHA_PREFIX = "bddb2d9dcbc44e8d"  # specimen.tsv source_files: 2416550405_table_all.glb
COARSE, FINE = 0.5, 0.1
SUB = int(round(COARSE / FINE))  # 5 x 5 fine cells per coarse cell
OPEN_MIN = (SUB * SUB + 1) // 2  # 13 of 25: at least half open
JIT = 1.2345e-6  # keeps ray positions off shared edges and vertices
HOLE_MIN_MM2, JOIN_MM, FINE_MARGIN = 100.0, 10.0, 3.0
STATIONS = (-0.40, -0.25, -0.10, 0.0, 0.10, 0.25, 0.40)
D = np.linspace(0.05, 0.95, 19)
HEAD = ("# il Cannone (Guarneri del Gesu 1743), ESRF BM18 synchrotron CT, data collection DOI 10.15151/ESRF-DC-2413577471 (CC BY 4.0), "
        "tensioned-string scan; top-plate mesh 2416550405_table_all.glb (SHA-256 %s...); fhole_openings.py; frame of the source files, mm\n")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def crossings(P, x0, z0, h, nx, nz):
    """Vertical-ray crossings of the triangles P (n, 3, 3) on the grid of cell centres x0 + (i + 1/2) h, z0 + (j + 1/2) h:
    the count per cell (nx, nz) and the crossing y per cell, sorted, as one array ordered by cell (CSR by the counts)."""
    txmin, txmax = P[:, :, 0].min(1), P[:, :, 0].max(1)
    tzmin, tzmax = P[:, :, 2].min(1), P[:, :, 2].max(1)
    sel = np.nonzero((txmax >= x0) & (txmin <= x0 + nx * h) & (tzmax >= z0) & (tzmin <= z0 + nz * h))[0]
    ia = np.clip(np.ceil((txmin[sel] - x0) / h - 0.5).astype(np.int64), 0, nx)
    ib = np.clip(np.floor((txmax[sel] - x0) / h - 0.5).astype(np.int64), -1, nx - 1)
    ja = np.clip(np.ceil((tzmin[sel] - z0) / h - 0.5).astype(np.int64), 0, nz)
    jb = np.clip(np.floor((tzmax[sel] - z0) / h - 0.5).astype(np.int64), -1, nz - 1)
    cx, cz = ib - ia + 1, jb - ja + 1
    ok = np.nonzero((cx > 0) & (cz > 0))[0]
    cells, ys = [], []
    for c0 in range(0, len(ok), 200000):
        g = ok[c0:c0 + 200000]
        nc = cx[g] * cz[g]
        rep = np.repeat(np.arange(len(g)), nc)
        off = np.arange(nc.sum()) - np.repeat(np.cumsum(nc) - nc, nc)
        ii = ia[g][rep] + off % cx[g][rep]
        jj = ja[g][rep] + off // cx[g][rep]
        px = x0 + (ii + 0.5) * h + JIT
        pz = z0 + (jj + 0.5) * h + JIT
        T = P[sel[g]][rep]
        ax, az, ay = T[:, 0, 0], T[:, 0, 2], T[:, 0, 1]
        bx, bz, by = T[:, 1, 0], T[:, 1, 2], T[:, 1, 1]
        qx, qz, qy = T[:, 2, 0], T[:, 2, 2], T[:, 2, 1]
        d = (bx - ax) * (qz - az) - (bz - az) * (qx - ax)
        nzd = np.abs(d) > 1e-14
        dd = np.where(nzd, d, 1.0)
        wb = np.where(nzd, ((px - ax) * (qz - az) - (pz - az) * (qx - ax)) / dd, -1.0)
        wc = np.where(nzd, ((bx - ax) * (pz - az) - (bz - az) * (px - ax)) / dd, -1.0)
        wa = 1.0 - wb - wc
        hit = nzd & (wa > 0) & (wb > 0) & (wc > 0)
        cells.append((ii * nz + jj)[hit])
        ys.append((wa * ay + wb * by + wc * qy)[hit])
    cells = np.concatenate(cells) if cells else np.zeros(0, np.int64)
    ys = np.concatenate(ys) if ys else np.zeros(0)
    o = np.lexsort((ys, cells))
    cnt = np.bincount(cells, minlength=nx * nz).reshape(nx, nz)
    return cnt, ys[o]


def label(mask):
    """4-connected components of a boolean grid: labels (0 = not in the mask) and the count."""
    lab = np.zeros(mask.shape, np.int32)
    n = 0
    nx, nz = mask.shape
    for i, j in zip(*np.nonzero(mask)):
        if lab[i, j]:
            continue
        n += 1
        lab[i, j] = n
        stack = [(i, j)]
        while stack:
            a, b = stack.pop()
            for p, q in ((a + 1, b), (a - 1, b), (a, b + 1), (a, b - 1)):
                if 0 <= p < nx and 0 <= q < nz and mask[p, q] and not lab[p, q]:
                    lab[p, q] = n
                    stack.append((p, q))
    return lab, n


def dilate(m, conn8=True):
    o = m.copy()
    o[1:, :] |= m[:-1, :]
    o[:-1, :] |= m[1:, :]
    o[:, 1:] |= m[:, :-1]
    o[:, :-1] |= m[:, 1:]
    if conn8:
        o[1:, 1:] |= m[:-1, :-1]
        o[:-1, :-1] |= m[1:, 1:]
        o[1:, :-1] |= m[:-1, 1:]
        o[:-1, 1:] |= m[1:, :-1]
    return o


def runs(mask, x0, z0, h, hid):
    """Row runs of a (nx, nz) mask: (hole, z centre, x start, x end) per maximal run of open cells along x."""
    out = []
    for j in range(mask.shape[1]):
        col = np.r_[0, mask[:, j].astype(np.int8), 0]
        d = np.diff(col)
        for a, b in zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]):
            out.append((hid, z0 + (j + 0.5) * h, x0 + a * h, x0 + b * h))
    return out


def basis(x, z, sx, sz):
    u, w = (x - sx[0]) / sx[1], (z - sz[0]) / sz[1]
    return np.stack([np.ones_like(u), u, w, u * u, u * w, w * w, u ** 3, u * u * w, u * w * w, w ** 3], -1)


def principal_axis(mask, x0, z0, h):
    """Area centroid, principal-axis angle from +z toward +x and extent along it of a plan mask."""
    ii, jj = np.nonzero(mask)
    X = np.c_[x0 + (ii + 0.5) * h, z0 + (jj + 0.5) * h]
    c = X.mean(0)
    ev, evec = np.linalg.eigh(np.cov((X - c).T))
    ax = evec[:, 1] if evec[1, 1] >= 0 else -evec[:, 1]
    proj = (X - c) @ ax
    return c, math.atan2(ax[0], ax[1]), float(proj.max() - proj.min())


def main():
    glb, out = sys.argv[1], sys.argv[2]
    digest = sha256(glb)
    if not digest.startswith(GLB_SHA_PREFIX):
        raise SystemExit("refusing %s: SHA-256 %s does not start with %s" % (glb, digest, GLB_SHA_PREFIX))
    os.makedirs(out, exist_ok=True)
    v, t, _ = io.glb_triangles(glb)
    P = v[t]
    print("mesh: %d vertices, %d triangles, SHA-256 %s" % (len(v), len(t), digest[:16]), flush=True)
    W = []  # walls table rows: hole, quantity, value, unit, note

    def put(hid, q, val, unit, note):
        W.append((str(hid), q, ("%.6g" % val) if isinstance(val, float) else str(val), unit, note))

    # 1. the plate at 0.5 mm: openings inside the plate
    lo, hi = P.reshape(-1, 3).min(0), P.reshape(-1, 3).max(0)
    cx0 = math.floor(lo[0] / COARSE) * COARSE - COARSE
    cz0 = math.floor(lo[2] / COARSE) * COARSE - COARSE
    cnx = int(math.ceil((hi[0] - cx0) / COARSE)) + 1
    cnz = int(math.ceil((hi[2] - cz0) / COARSE)) + 1
    ccnt, _ = crossings(P, cx0, cz0, COARSE, cnx, cnz)
    lab, n = label(ccnt == 0)
    edge = set(np.unique(np.r_[lab[0], lab[-1], lab[:, 0], lab[:, -1]]).tolist()) - {0}
    comps = []
    for c in range(1, n + 1):
        if c in edge:
            continue
        ii, jj = np.nonzero(lab == c)
        comps.append(dict(id=c, area=len(ii) * COARSE ** 2, x=cx0 + (ii.mean() + 0.5) * COARSE, z=cz0 + (jj.mean() + 0.5) * COARSE,
                          cells=(ii, jj)))
    big = sorted([c for c in comps if c["area"] >= HOLE_MIN_MM2], key=lambda c: c["x"])
    if len(big) != 2:
        raise SystemExit("expected two openings of %g mm2 or more inside the plate, found %d: %s" % (HOLE_MIN_MM2, len(big), [(c["area"], c["x"], c["z"]) for c in big]))
    groups = {1: [big[0]], 2: [big[1]]}
    left = []
    for c in comps:
        if c["area"] >= HOLE_MIN_MM2:
            continue
        ii, jj = c["cells"]
        best = None
        for hid in (1, 2):
            bi, bj = big[hid - 1]["cells"]
            d = np.sqrt(((ii[:, None] - bi[None, :]) * COARSE) ** 2 + ((jj[:, None] - bj[None, :]) * COARSE) ** 2).min() if len(ii) * len(bi) < 4e7 else 1e9
            if d <= JOIN_MM and (best is None or d < best[1]):
                best = (hid, d)
        if best:
            groups[best[0]].append(c)
        else:
            left.append(c)
    print("openings inside the plate at 0.5 mm: %d; f-holes %.2f and %.2f mm2; joined to them %s; others %d (%.2f mm2 in all)"
          % (len(comps), big[0]["area"], big[1]["area"], {h: [round(c["area"], 2) for c in g[1:]] for h, g in groups.items()},
             len(left), sum(c["area"] for c in left)), flush=True)
    put(0, "plate-openings-other-count", len(left), "1", "openings inside the plate at 0.5 mm that are neither f-hole nor within %g mm of one (not used)" % JOIN_MM)
    put(0, "plate-openings-other-area", float(sum(c["area"] for c in left)), "mm2", "their total plan area at 0.5 mm")
    for c in left:
        put(0, "plate-opening-other", "%.1f %.1f %.2f" % (c["x"], c["z"], c["area"]), "mm mm mm2", "centroid x, z and area at 0.5 mm of an opening inside the plate that is not used")

    rows05, rows01 = [], []
    res = {}
    for hid in (1, 2):
        ii = np.concatenate([c["cells"][0] for c in groups[hid]])
        jj = np.concatenate([c["cells"][1] for c in groups[hid]])
        # 2. the hole's box at 0.1 mm, edges on the 0.5 mm grid
        bx0 = cx0 + ii.min() * COARSE - FINE_MARGIN
        bz0 = cz0 + jj.min() * COARSE - FINE_MARGIN
        bx1 = cx0 + (ii.max() + 1) * COARSE + FINE_MARGIN
        bz1 = cz0 + (jj.max() + 1) * COARSE + FINE_MARGIN
        nx, nz = int(round((bx1 - bx0) / FINE)), int(round((bz1 - bz0) / FINE))
        cnt, ys = crossings(P, bx0, bz0, FINE, nx, nz)
        flab, fn = label(cnt == 0)
        fedge = set(np.unique(np.r_[flab[0], flab[-1], flab[:, 0], flab[:, -1]]).tolist()) - {0}
        seed = np.zeros((nx, nz), bool)
        for a, b in zip(ii, jj):
            fi = int(round((cx0 + a * COARSE - bx0) / FINE))
            fj = int(round((cz0 + b * COARSE - bz0) / FINE))
            seed[fi:fi + SUB, fj:fj + SUB] = True
        keep = sorted(set(np.unique(flab[seed & (flab > 0)]).tolist()) - fedge)
        hole = np.isin(flab, keep)
        area01 = float(hole.sum() * FINE * FINE)
        # 3. coarsening to 0.5 mm
        cnx_, cnz_ = nx // SUB, nz // SUB
        nopen = hole[:cnx_ * SUB, :cnz_ * SUB].reshape(cnx_, SUB, cnz_, SUB).sum((1, 3))
        m05 = nopen >= OPEN_MIN
        m05_edge = nopen == SUB * SUB
        area05, area05e = float(m05.sum() * COARSE ** 2), float(m05_edge.sum() * COARSE ** 2)
        # the rule's own error: the same coarsening on the 25 grids shifted by whole 0.1 mm cells against the 0.1 mm area
        shifted = []
        for si in range(SUB):
            for sj in range(SUB):
                hs = hole[si:si + (cnx_ - 1) * SUB, sj:sj + (cnz_ - 1) * SUB]
                ns = hs.reshape(cnx_ - 1, SUB, cnz_ - 1, SUB).sum((1, 3))
                shifted.append(float((ns >= OPEN_MIN).sum() * COARSE ** 2) - area01)
        shifted = np.array(shifted)
        l05, n05 = label(m05)
        rows05 += runs(m05, bx0, bz0, COARSE, hid)
        rows01 += runs(hole, bx0, bz0, FINE, hid)
        print("hole %d: 0.1 mm open %.2f mm2 in %d component(s) (%d open components in the box, %d at its edge); 0.5 mm raster %.2f mm2 "
              "(%d cells, %d 4-connected pieces); edge-as-wood rule %.2f mm2" % (hid, area01, len(keep), fn, len(fedge), area05, m05.sum(), n05, area05e), flush=True)
        put(hid, "box", "%.1f %.1f %.1f %.1f" % (bx0, bx1, bz0, bz1), "mm", "x from, x to, z from, z to of the 0.1 mm grid")
        put(hid, "area-0p1mm", area01, "mm2", "cells of 0.1 mm open through the plate (%d component(s) of the 0.1 mm grid)" % len(keep))
        put(hid, "area-0p5mm", area05, "mm2", "top-plate-fholes.tsv: 0.5 mm cells with at least %d of 25 cells of 0.1 mm open" % OPEN_MIN)
        put(hid, "pieces-0p5mm", n05, "1", "4-connected pieces of the 0.5 mm raster (more than 1: an eye whose neck is under half a cell wide)")
        put(hid, "area-0p5mm-edge-as-wood", area05e, "mm2", "the replaced rule on the same data: a 0.5 mm cell open only when all 25 of its cells are")
        put(hid, "area-0p5mm-rule-error", float(np.abs(shifted).max()), "mm2", "largest |0.5 mm area - 0.1 mm area| of the rule over the 25 grid offsets of whole 0.1 mm cells (mean %+.2f, RMS %.2f mm2); the specimen grid's own offset gives %+.2f"
            % (shifted.mean(), np.sqrt((shifted ** 2).mean()), area05 - area01))
        # 4. walls
        off = np.cumsum(cnt.ravel()) - cnt.ravel()
        Y = np.full((nx * nz, 6), np.nan)
        for m in range(6):
            s = cnt.ravel() > m
            Y[s, m] = ys[off[s] + m]
        Y = Y.reshape(nx, nz, 6)
        xc = bx0 + (np.arange(nx) + 0.5) * FINE
        zc = bz0 + (np.arange(nz) + 0.5) * FINE
        XX, ZZ = np.meshgrid(xc, zc, indexing="ij")
        dist = np.full((nx, nz), 999, np.int32)
        dist[hole] = 0
        front = hole.copy()
        for s in range(1, 61):
            front = dilate(front) | front
            dist[front & (dist == 999)] = s
        full = cnt == 2
        ylo, yhi = Y[:, :, 0], Y[:, :, 1]
        ann = full & (dist >= 15) & (dist <= 50)
        sx, sz = (float(XX[ann].mean()), 20.0), (float(ZZ[ann].mean()), 40.0)
        B = basis(XX[ann], ZZ[ann], sx, sz)

        def robust(y):
            k = np.ones(len(y), bool)
            for _ in range(3):
                c, *_ = np.linalg.lstsq(B[k], y[k], rcond=None)
                r = y - B @ c
                k = np.abs(r) <= 4 * r[k].std()
            return c, r[k], int((~k).sum())

        co, ro, dro = robust(ylo[ann])
        ci, ri, dri = robust(yhi[ann])
        BA = basis(XX, ZZ, sx, sz)
        Yo, Yi = BA @ co, BA @ ci
        eps = 1e-3
        gx = (basis(XX + eps, ZZ, sx, sz) @ co - basis(XX - eps, ZZ, sx, sz) @ co) / (2 * eps)
        gz = (basis(XX, ZZ + eps, sx, sz) @ co - basis(XX, ZZ - eps, sx, sz) @ co) / (2 * eps)
        nrm = np.sqrt(1 + gx * gx + gz * gz)
        cos_t = 1 / nrm
        tn = (Yi - Yo) * cos_t
        rim = full & (dist >= 8) & (dist <= 20)
        t_rim_v = yhi[rim] - ylo[rim]
        t_rim_n = t_rim_v * cos_t[rim]
        ring = (dist >= 1) & (dist <= 8)
        wgw = int((cnt[ring] == 4).sum())
        put(hid, "surface-fit-cells", int(ann.sum()), "1", "cells 1.5 to 5 mm from the hole crossed once through the wood; cubic fits in (x, z)")
        put(hid, "surface-fit-rms-outer", float(ro.std()), "mm", "residual RMS of the outer fit after dropping %d cells above 4 sigma" % dro)
        put(hid, "surface-fit-rms-inner", float(ri.std()), "mm", "residual RMS of the inner fit after dropping %d cells above 4 sigma" % dri)
        put(hid, "outer-slope-median", float(np.degrees(np.arccos(np.median(cos_t[hole])))), "deg", "fitted outer surface over the hole")
        put(hid, "rim-thickness-vertical", float(np.median(t_rim_v)), "mm", "median over cells 0.8 to 2 mm from the hole (p10 %.3f, p90 %.3f)" % (np.percentile(t_rim_v, 10), np.percentile(t_rim_v, 90)))
        put(hid, "rim-thickness-normal", float(np.median(t_rim_n)), "mm", "the same along the outer surface's normal (p10 %.3f, p90 %.3f)" % (np.percentile(t_rim_n, 10), np.percentile(t_rim_n, 90)))
        put(hid, "wood-gap-wood-cells", wgw, "1", "cells of 0.1 mm within 0.8 mm of the hole whose ray crosses wood, a gap and wood (unresolved: a repair, a crack or varnish)")

        def is_open(px, py, pz):
            i = np.floor((px - bx0) / FINE).astype(np.int64)
            j = np.floor((pz - bz0) / FINE).astype(np.int64)
            inside = (i >= 0) & (i < nx) & (j >= 0) & (j < nz)
            i, j = np.clip(i, 0, nx - 1), np.clip(j, 0, nz - 1)
            yy = Y[i, j]
            wood = np.zeros(np.shape(px), bool)
            for m in (0, 2, 4):
                wood |= (yy[..., m] <= py) & (py <= yy[..., m + 1])
            return inside & ~wood

        near = dist <= 30
        fam = {}
        for name in ("vertical", "normal"):
            A = []
            for d in D:
                if name == "vertical":
                    op = is_open(XX, Yo + d * (Yi - Yo), ZZ) & near
                    A.append(float(op.sum() * FINE * FINE))
                else:
                    op = is_open(XX + d * tn * (-gx / nrm), Yo + d * tn / nrm, ZZ + d * tn * (-gz / nrm)) & near
                    A.append(float(np.sum(op / cos_t) * FINE * FINE))
            A = np.array(A)
            t_ref = float(np.median(((Yi - Yo) if name == "vertical" else tn)[hole]))
            R = float(np.trapezoid(t_ref / np.r_[A[0], A, A[-1]], np.r_[0.0, D, 1.0]))
            fam[name] = dict(A=A, t=t_ref, R=R)
        c, ang, length = principal_axis(hole, bx0, bz0, FINE)
        u = np.array([math.sin(ang), math.cos(ang)])
        w = np.array([math.cos(ang), -math.sin(ang)])
        r = np.arange(-15, 15.0001, 0.05)
        mid = int(np.argmin(np.abs(r)))
        for name in ("vertical", "normal"):
            drifts, secs = [], []
            for f in STATIONS:
                s = f * length
                px0, pz0 = c[0] + s * u[0] + r * w[0], c[1] + s * u[1] + r * w[1]
                i = np.clip(np.floor((px0 - bx0) / FINE).astype(int), 0, nx - 1)
                j = np.clip(np.floor((pz0 - bz0) / FINE).astype(int), 0, nz - 1)
                edges = []
                for d in (0.1, 0.5, 0.9):
                    if name == "vertical":
                        op = is_open(px0, Yo[i, j] + d * (Yi[i, j] - Yo[i, j]), pz0)
                    else:
                        op = is_open(px0 + d * tn[i, j] * (-gx[i, j] / nrm[i, j]), Yo[i, j] + d * tn[i, j] / nrm[i, j],
                                     pz0 + d * tn[i, j] * (-gz[i, j] / nrm[i, j]))
                    idx = np.nonzero(op)[0]
                    if len(idx) == 0:
                        edges.append(None)
                        continue
                    c0 = mid if op[mid] else idx[np.argmin(np.abs(idx - mid))]
                    a, b = c0, c0
                    while a > 0 and op[a - 1]:
                        a -= 1
                    while b < len(r) - 1 and op[b + 1]:
                        b += 1
                    edges.append((r[a] - 0.025, r[b] + 0.025))
                secs.append((f, s, edges))
                if all(e is not None for e in edges) and abs(f) <= 0.25:
                    drifts.append(0.5 * (edges[2][0] + edges[2][1]) - 0.5 * (edges[0][0] + edges[0][1]))
            fam[name]["drift"] = float(np.mean(drifts)) if drifts else float("nan")
            fam[name]["sections"] = secs
        wall = min(("vertical", "normal"), key=lambda k: abs(fam[k]["drift"]))
        Rk = fam[wall]["R"]
        res[hid] = dict(R=Rk, A05=area05, A01=area01, wall=wall, fam=fam)
        for name in ("vertical", "normal"):
            F = fam[name]
            put(hid, "through-resistance-%s" % name, F["R"], "1/mm", "integral of t dd / A(d) over %s layers at t %.3f mm; A(d) at d 0.05/0.20/0.35/0.50/0.65/0.80/0.95: %s mm2"
                % (name, F["t"], " ".join("%.1f" % a for a in F["A"][::3])))
            put(hid, "wall-drift-%s" % name, F["drift"], "mm", "mean shift of the opening's mid-point from d 0.1 to d 0.9 across the axis, stations |s| <= 0.25 L, %s layers" % name)
            for f, s, edges in F["sections"]:
                if any(e is None for e in edges):
                    put(hid, "station%+.2fL-%s" % (f, name), "none", "mm", "no opening on the section at some depth")
                    continue
                put(hid, "station%+.2fL-%s" % (f, name), "%.2f %.2f %.2f %+.2f %+.2f" % (edges[0][1] - edges[0][0], edges[1][1] - edges[1][0], edges[2][1] - edges[2][0],
                                                                                         edges[0][0] - edges[2][0], edges[2][1] - edges[0][1]),
                    "mm", "s = %.1f mm along the axis: width at d 0.1, 0.5, 0.9 and the flare of each side (cavity side wider: positive)" % s)
        put(hid, "wall-direction", wall, "-", "the layer family with the smaller |drift|")
        put(hid, "through-resistance", Rk, "1/mm", "the %s family's" % wall)
        put(hid, "bore", Rk * area05, "mm", "straight bore of the 0.5 mm raster with the same through-resistance: R x %.2f mm2" % area05)
        put(hid, "bore-at-0p1mm-area", Rk * area01, "mm", "R x the 0.1 mm open area %.2f mm2" % area01)
        put(hid, "axis", "%.4f %.4f %.6f %.4f" % (c[0], c[1], ang, length), "mm mm rad mm", "area centroid x, z; principal-axis angle from +z toward +x; extent along it (0.1 mm cells)")
    t_pair = (res[1]["A05"] + res[2]["A05"]) / (1 / res[1]["R"] + 1 / res[2]["R"])
    put(0, "bore-pair", t_pair, "mm", "one bore for both holes with their parallel bore conductance: (A1 + A2) / (1/R1 + 1/R2), A the 0.5 mm raster areas")
    for name, rows, h in (("top-plate-fholes.tsv", rows05, COARSE), ("top-plate-fholes-0p1mm.tsv", rows01, FINE)):
        with open(os.path.join(out, name), "w") as f:
            f.write(HEAD % digest[:16])
            f.write("hole\tz_mm\tx_start_mm\tx_end_mm\n")
            for hid, z, a, b in rows:
                f.write("%d\t%.2f\t%.2f\t%.2f\n" % (hid, z, a, b))
        print("wrote", name, len(rows), "rows", flush=True)
    with open(os.path.join(out, "top-plate-fhole-walls.tsv"), "w") as f:
        f.write(HEAD % digest[:16])
        f.write("hole\tquantity\tvalue\tunit\tnote\n")
        for row in W:
            f.write("\t".join(row) + "\n")
    print("wrote top-plate-fhole-walls.tsv", len(W), "rows")
    for hid in (1, 2):
        q = res[hid]
        print("hole %d: walls %s (drift vertical %+.2f, normal %+.2f mm); R vertical %.5f, normal %.5f /mm; bore %.3f mm at %.2f mm2"
              % (hid, q["wall"], q["fam"]["vertical"]["drift"], q["fam"]["normal"]["drift"], q["fam"]["vertical"]["R"], q["fam"]["normal"]["R"], q["R"] * q["A05"], q["A05"]))
    print("pair bore %.4f mm" % t_pair)


if __name__ == "__main__":
    main()
