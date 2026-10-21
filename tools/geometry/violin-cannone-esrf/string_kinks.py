"""string_kinks.py : per string, section centroids every 1 mm along z; the bridge and nut contact points are where the centre line of a string bends: report the z of the
largest change of slope dy/dz (and dx/dz) in two windows (near the bridge, near the nut), and the straight-segment slopes either side. Read-only."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cannone_io as io  # noqa: E402
import cannone_sections as cs  # noqa: E402

RAW2 = "/data/tmp/wa-inst/geometry/violin-cannone-esrf/raw/stage2"
OUT = sys.argv[1] if len(sys.argv) > 1 else None
rows_out = []
for nm, fid in (("G", "2416550410_corde_sol"), ("D", "2416550448_corde_re"), ("A", "2416550453_corde_la"), ("E", "2416550432_corde_mi")):
    v, t, _ = io.glb_triangles(RAW2 + "/%s.glb" % fid)
    print("== string", nm, "z range %.2f .. %.2f, y range %.2f .. %.2f" % (v[:, 2].min(), v[:, 2].max(), v[:, 1].min(), v[:, 1].max()))
    rows = np.array(cs.sections(v, t, 2, 1.0))
    z, area, cx, cy = rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3]
    mask = area < 3 * np.median(area)  # drop stations that include the winding ends or the loops
    for lab, (a, b) in (("bridge window", (-125, -90)), ("nut window", (190, 240))):
        m = (z >= a) & (z <= b) & mask
        zz, yy, xx = z[m], cy[m], cx[m]
        # slope on each side by fitting 12-station windows; kink = z maximising |slope_right - slope_left|
        best = None
        for i in range(12, len(zz) - 12):
            sl = np.polyfit(zz[i - 12:i], yy[i - 12:i], 1)[0]
            sr = np.polyfit(zz[i:i + 12], yy[i:i + 12], 1)[0]
            d = abs(sr - sl)
            if best is None or d > best[0]:
                best = (d, zz[i], sl, sr, np.interp(zz[i], zz, xx), np.interp(zz[i], zz, yy))
        print("  %s: kink at z=%.1f  (slope dy/dz %.4f -> %.4f, change %.4f)  x=%.2f y=%.2f" % (lab, best[1], best[2], best[3], best[0], best[4], best[5]))
        rows_out.append((nm, "bridge" if lab.startswith("bridge") else "nut", best[1], best[4], best[5], best[2], best[3]))
    print("  median section area %.3f mm2 -> diameter %.3f mm; min area %.4f" % (np.median(area[mask]), 2 * np.sqrt(np.median(area[mask]) / np.pi), area.min()))

if OUT:
    with open(OUT, "w") as f:
        f.write("frankensim.specimen-string-contacts.v1\n")
        f.write("string\tcontact\tz_mm\tx_mm\ty_mm\tslope_dy_dz_before\tslope_dy_dz_after\n")
        for r in rows_out:
            f.write("%s\t%s\t%.1f\t%.2f\t%.2f\t%.4f\t%.4f\n" % r)
    print("wrote", OUT)
