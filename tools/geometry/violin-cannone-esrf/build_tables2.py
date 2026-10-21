"""build_tables2.py OUTDIR : second set of il Cannone tables: rib and lining cross-sections around the body (planes through the vertical axis every 2 degrees), the four strings,
the neck, the fingerboard and the heel block (sections along z), and the volume of every part. Same frame, units and source as build_tables.py."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cannone_io as io  # noqa: E402
import cannone_sections as cs  # noqa: E402

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
RAW1 = "/data/tmp/wa-inst/geometry/violin-cannone-esrf/raw/stage1"
RAW2 = "/data/tmp/wa-inst/geometry/violin-cannone-esrf/raw/stage2"
HEAD = "# il Cannone (Guarneri del Gesu 1743), ESRF BM18 synchrotron CT, data collection DOI 10.15151/ESRF-DC-2413577471 (CC BY 4.0), tensioned-string scan; frame of the source files, mm\n"


def write(name, header, rows):
    with open(os.path.join(out, name), "w") as f:
        f.write(HEAD)
        f.write(header + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print("wrote", name, len(rows), "rows", flush=True)


def volume(v, t):
    a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


vols = []


def load(path, name):
    v, t, _ = io.glb_triangles(path)
    vols.append((name, volume(v, t), len(t), v.min(0), v.max(0)))
    return v, t


# --- rings: ribs and linings, sections by half-planes through the vertical axis (y), every 2 degrees
def ring(name, path, split_y=None):
    v, t = load(path, name)
    lo, hi = v.min(0), v.max(0)
    cx, cz = 0.5 * (lo[0] + hi[0]), 0.5 * (lo[2] + hi[2])
    cen = v[t].mean(1)
    th_c = np.degrees(np.arctan2(cen[:, 2] - cz, cen[:, 0] - cx)) % 360
    rows = []
    for ang in range(0, 360, 2):
        sel = np.nonzero(np.abs(((th_c - ang + 180) % 360) - 180) < 5.0)[0]
        if len(sel) == 0:
            continue
        w = v[t[sel]]  # (k, 3, 3) world
        c, s = np.cos(np.radians(ang)), np.sin(np.radians(ang))
        xr = (w[:, :, 0] - cx) * c + (w[:, :, 2] - cz) * s  # radial
        zr = -(w[:, :, 0] - cx) * s + (w[:, :, 2] - cz) * c  # tangential
        W = np.stack([xr, w[:, :, 1], zr], axis=2)
        P, Q = cs.plane_segments(W[:, 0], W[:, 1], W[:, 2], 2, 0.0)
        keep = (P[:, 0] > 5.0) & (Q[:, 0] > 5.0) if len(P) else np.zeros(0, bool)
        P, Q = P[keep], Q[keep]
        if len(P) == 0:
            continue
        # world position of the section plane point at radial r: (cx + r cos, cz + r sin)
        parts = [("", np.ones(len(P), bool))] if split_y is None else [("top", 0.5 * (P[:, 1] + Q[:, 1]) < split_y), ("bottom", 0.5 * (P[:, 1] + Q[:, 1]) >= split_y)]
        row = ["%d" % ang, "%.2f" % (cx + 0.0), "%.2f" % cz]
        for lab, m in parts:
            if m.sum() < 3:
                row += ["", "", "", "", "", ""]
                continue
            Pm, Qm = P[m], Q[m]
            area = abs(cs.green_area(Pm, Qm, 0, 1))
            allp = np.concatenate([Pm, Qm])
            row += ["%.2f" % area, "%.2f" % allp[:, 0].min(), "%.2f" % allp[:, 0].max(), "%.2f" % allp[:, 1].min(), "%.2f" % allp[:, 1].max(), "%.2f" % (area / max(allp[:, 1].max() - allp[:, 1].min(), 1e-9))]
        rows.append(tuple(row))
    cols = ["area_mm2", "r_min_mm", "r_max_mm", "y_min_mm", "y_max_mm", "mean_width_mm"]
    if split_y is None:
        header = "theta_deg\tcentre_x_mm\tcentre_z_mm\t" + "\t".join(cols)
    else:
        header = "theta_deg\tcentre_x_mm\tcentre_z_mm\t" + "\t".join("top_" + c for c in cols) + "\t" + "\t".join("bottom_" + c for c in cols)
    write(name + "-sections.tsv", header, rows)
    return cx, cz


ring("ribs", RAW2 + "/2416550419_eclisses_all.glb")
ring("linings", RAW2 + "/2416550428_contre_elisses.glb", split_y=8.5)

# --- strings: sections along z every 20 mm (area -> diameter, centre)
for nm, fid in (("string-G", "2416550410_corde_sol"), ("string-D", "2416550448_corde_re"), ("string-A", "2416550453_corde_la"), ("string-E", "2416550432_corde_mi")):
    v, t = load(RAW2 + "/%s.glb" % fid, nm)
    rows = cs.sections(v, t, 2, 20.0)
    write(nm + "-sections.tsv", "z_mm\tarea_mm2\tcentroid_x_mm\tcentroid_y_mm\tx_min_mm\tx_max_mm\ty_min_mm\ty_max_mm", [tuple("%.3f" % q for q in r) for r in rows])

# --- neck, fingerboard, heel block: sections along z
for nm, path, step in (("neck", RAW2 + "/2416550416_manche_all.glb", 10.0), ("fingerboard", RAW2 + "/2416550444_touche.glb", 10.0), ("heel-block", RAW2 + "/2416550411_talon_bloc_principal.glb", 5.0)):
    v, t = load(path, nm)
    rows = cs.sections(v, t, 2, step)
    write(nm + "-sections.tsv", "z_mm\tarea_mm2\tcentroid_x_mm\tcentroid_y_mm\tx_min_mm\tx_max_mm\ty_min_mm\ty_max_mm", [tuple("%.2f" % q for q in r) for r in rows])

# --- volumes of the parts already loaded here plus the ones from stage 1 and the plates
for nm, path in (("top-plate", RAW2 + "/2416550405_table_all.glb"), ("back-plate", RAW2 + "/2416550427_dos_all.glb"), ("bass-bar", RAW1 + "/2416550443_barre_harmonie.glb"),
                 ("soundpost", RAW1 + "/2416550418_ame.glb"), ("bridge", RAW1 + "/2416550436_chevalet.glb"), ("tailpiece", RAW1 + "/2416550421_cordier.glb")):
    load(path, nm)
rows = []
for nm, vol, n, lo, hi in vols:
    rows.append((nm, "%.1f" % vol, "%d" % n, "%.2f" % lo[0], "%.2f" % hi[0], "%.2f" % lo[1], "%.2f" % hi[1], "%.2f" % lo[2], "%.2f" % hi[2]))
write("part-volumes.tsv", "part\tvolume_mm3\ttriangles\tx_min_mm\tx_max_mm\ty_min_mm\ty_max_mm\tz_min_mm\tz_max_mm", rows)
