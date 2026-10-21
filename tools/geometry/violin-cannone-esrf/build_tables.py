"""build_tables.py OUTDIR : turns the il Cannone plate height maps and the bass bar, soundpost and bridge meshes into the data tables of the repo description (TSV, mm, the frame
of the source files: x across the body, y vertical (top plate at negative y, back at positive y), z along the body (tail end at negative z)). Reads only the files fetched from
the public ESRF data collection (DOI 10.15151/ESRF-DC-2413577471, CC BY 4.0). Nothing is fitted or smoothed except where a column header says so."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cannone_io as io  # noqa: E402
import cannone_sections as cs  # noqa: E402

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
DER = "/data/tmp/wa-inst/geometry/violin-cannone-esrf/derived"
RAW1 = "/data/tmp/wa-inst/geometry/violin-cannone-esrf/raw/stage1"
HEAD = "# il Cannone (Guarneri del Gesu 1743), ESRF BM18 synchrotron CT, data collection DOI 10.15151/ESRF-DC-2413577471 (CC BY 4.0), tensioned-string scan; frame of the source files, mm\n"


def write(name, header, rows):
    with open(os.path.join(out, name), "w") as f:
        f.write(HEAD)
        f.write(header + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print("wrote", name, len(rows), "rows")


def fill_outside(mask):
    free = ~mask
    o = np.zeros_like(mask)
    o[0, :], o[-1, :], o[:, 0], o[:, -1] = free[0, :], free[-1, :], free[:, 0], free[:, -1]
    while True:
        new = o.copy()
        new[1:, :] |= o[:-1, :] & free[1:, :]
        new[:-1, :] |= o[1:, :] & free[:-1, :]
        new[:, 1:] |= o[:, :-1] & free[:, 1:]
        new[:, :-1] |= o[:, 1:] & free[:, :-1]
        if (new == o).all():
            return o
        o = new


def erosion_distance(F):
    dist = np.zeros(F.shape, np.int32)
    cur = F.copy()
    s = 0
    while cur.any() and s < 600:
        s += 1
        dist[cur] = s
        er = cur.copy()
        er[1:, :] &= cur[:-1, :]
        er[:-1, :] &= cur[1:, :]
        er[:, 1:] &= cur[:, :-1]
        er[:, :-1] &= cur[:, 1:]
        cur = er
    return dist


def plate(name, npz, kind):
    d = np.load(npz, allow_pickle=True)
    x0, z0, h = float(d["x0"]), float(d["z0"]), float(d["cell"])
    outer, inner, cnt = d["outer"], d["inner"], d["count"]
    nx, nz = outer.shape
    M = cnt > 0
    F = ~fill_outside(M)
    dist = erosion_distance(F) * h  # mm from the footprint edge (cell centres)
    # outline: the first and last occupied cell of every 2nd row (1 mm in z)
    rows = []
    for j in range(0, nz, 2):
        ii = np.nonzero(F[:, j])[0]
        if len(ii):
            rows.append(("%.2f" % (z0 + (j + 0.5) * h), "%.2f" % (x0 + ii[0] * h), "%.2f" % (x0 + (ii[-1] + 1) * h)))
    write("%s-outline.tsv" % name, "z_mm\tx_left_mm\tx_right_mm", rows)
    # f-holes
    holes = F & ~M
    lab = np.zeros(holes.shape, np.int32)
    k = 0
    for i, j in zip(*np.nonzero(holes)):
        if lab[i, j]:
            continue
        k += 1
        stack = [(i, j)]
        lab[i, j] = k
        while stack:
            a, b = stack.pop()
            for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                p, q = a + da, b + db
                if 0 <= p < nx and 0 <= q < nz and holes[p, q] and not lab[p, q]:
                    lab[p, q] = k
                    stack.append((p, q))
    big = sorted([(int((lab == c).sum()), c) for c in range(1, k + 1)], reverse=True)
    big = [(s, c) for s, c in big if s * h * h >= 100.0]
    if big:
        frows, info = [], []
        for hid, (s, c) in enumerate(sorted(big, key=lambda sc: np.nonzero(lab == sc[1])[0].mean()), 1):
            ii, jj = np.nonzero(lab == c)
            info.append((hid, s * h * h, x0 + (ii.mean() + 0.5) * h, z0 + (jj.mean() + 0.5) * h))
            for j in np.unique(jj):
                col = np.nonzero(lab[:, j] == c)[0]
                for a, b in cs.runs_of(np.isin(np.arange(nx), col)):
                    frows.append((str(hid), "%.2f" % (z0 + (j + 0.5) * h), "%.2f" % (x0 + a * h), "%.2f" % (x0 + (b + 1) * h)))
        write("%s-fholes.tsv" % name, "hole\tz_mm\tx_start_mm\tx_end_mm", frows)
        print("   f-holes (id, area mm2, centroid x, centroid z):", [(i, round(a, 1), round(x, 1), round(z, 1)) for i, a, x, z in info])
    # maps: 3 mm blocks (6 x 6 cells); medians over the valid cells of the block
    tv = np.abs(outer - inner)
    tv[(cnt < 2) | ~M] = np.nan
    gx = (np.roll(outer, -2, 0) - np.roll(outer, 2, 0)) / (4 * h)
    gz = (np.roll(outer, -2, 1) - np.roll(outer, 2, 1)) / (4 * h)
    tn = tv / np.sqrt(1 + gx * gx + gz * gz)
    rows = []
    B = 6
    for bi in range(0, nx - B + 1, B):
        for bj in range(0, nz - B + 1, B):
            blk = (slice(bi, bi + B), slice(bj, bj + B))
            ok = ~np.isnan(outer[blk]) & ~np.isnan(inner[blk]) & ~np.isnan(tn[blk])
            if ok.sum() < 18:
                continue
            xo = x0 + (bi + B / 2) * h
            zo = z0 + (bj + B / 2) * h
            rows.append(("%.1f" % xo, "%.1f" % zo, "%.2f" % np.median(outer[blk][ok]), "%.2f" % np.median(inner[blk][ok]), "%.2f" % np.median(tn[blk][ok]),
                         "%.1f" % np.min(dist[blk]), "%d" % ok.sum()))
    write("%s-maps.tsv" % name, "x_mm\tz_mm\touter_y_mm\tinner_y_mm\tthickness_normal_mm\tdist_to_outline_mm\tvalid_cells_of_36", rows)
    return d, F, M, dist, tv, tn


for nm, npz, kind in (("top-plate", DER + "/table_maps2.npz", "TOP"), ("back-plate", DER + "/dos_maps2.npz", "BACK")):
    print("==", nm)
    plate(nm, npz, kind)

# bass bar, soundpost: sections
v, t, _ = io.glb_triangles(RAW1 + "/2416550443_barre_harmonie.glb")
rows = cs.sections(v, t, 2, 5.0)
write("bass-bar-sections.tsv", "z_mm\tarea_mm2\tcentroid_x_mm\tcentroid_y_mm\tx_min_mm\tx_max_mm\ty_min_mm\ty_max_mm", [tuple("%.2f" % q for q in r) for r in rows])
v, t, _ = io.glb_triangles(RAW1 + "/2416550418_ame.glb")
rows = cs.sections(v, t, 1, 2.0)
write("soundpost-sections.tsv", "y_mm\tarea_mm2\tcentroid_x_mm\tcentroid_z_mm\tx_min_mm\tx_max_mm\tz_min_mm\tz_max_mm", [tuple("%.2f" % q for q in r) for r in rows])
# bridge: sections by y, silhouette in the x-y plane
v, t, _ = io.glb_triangles(RAW1 + "/2416550436_chevalet.glb")
rows = cs.sections(v, t, 1, 1.0)
write("bridge-sections.tsv", "y_mm\tarea_mm2\tcentroid_x_mm\tcentroid_z_mm\tx_min_mm\tx_max_mm\tz_min_mm\tz_max_mm", [tuple("%.2f" % q for q in r) for r in rows])
o0, o1, occ = cs.silhouette_runs(v, t, (0, 1), 0.25, 0.1)
rows = []
for j in range(occ.shape[1]):
    for a, b in cs.runs_of(occ[:, j]):
        rows.append(("%.2f" % (o1 + (j + 0.5) * 0.25), "%.2f" % (o0 + a * 0.25), "%.2f" % (o0 + (b + 1) * 0.25)))
write("bridge-silhouette.tsv", "y_mm\tx_start_mm\tx_end_mm", rows)
