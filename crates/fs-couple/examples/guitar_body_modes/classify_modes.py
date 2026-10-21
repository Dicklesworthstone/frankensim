"""Classify the air-coupled body modes written by guitar_body_modes --shapes: for the top and the back plate count the
antinodal regions along the long axis (strip |y| < 15 mm) and across the lower bout (strip |x - 0.37 m| < 15 mm), and give the
phase relation (sign of the mean vertical displacement) of top and back.
usage: python3 -I classify_modes.py SHAPES.tsv
"""
import sys

rows = [l.split("\t") for l in open(sys.argv[1], encoding="utf-8").read().splitlines()]
hdr = rows[0]
nm = len(hdr) - 4
data = {"top": [], "back": []}
for r in rows[1:]:
    data[r[0]].append((float(r[2]), float(r[3]), [float(v) for v in r[4:]]))


def antinodes(points, key, band_axis, band_centre, band_half, mode, level=0.10):
    sel = sorted(
        ((p[key], p[2][mode]) for p in points if abs(p[band_axis] - band_centre) < band_half),
        key=lambda t: t[0],
    )
    if len(sel) < 4:
        return 0
    peak = max(abs(v) for _, v in sel)
    count = 0
    sign = 0
    for _, v in sel:
        if abs(v) < level * peak:
            continue
        s = 1 if v > 0 else -1
        if s != sign:
            count += 1
            sign = s
    return count


def rms(points, mode):
    return (sum(p[2][mode] ** 2 for p in points) / len(points)) ** 0.5


def meanw(points, mode):
    return sum(p[2][mode] for p in points) / len(points)


print("mode  top_rms  back_rms  top_mean  back_mean  phase  top_long  top_cross  back_long  back_cross")
for m in range(nm):
    t, b = data["top"], data["back"]
    tm, bm = meanw(t, m), meanw(b, m)
    phase = "in" if tm * bm > 0 else "anti"
    print(
        "%4d %8.2e %9.2e %9.2e %10.2e  %-4s %8d %9d %10d %11d"
        % (
            m + 1,
            rms(t, m),
            rms(b, m),
            tm,
            bm,
            phase,
            antinodes(t, 0, 1, 0.0, 0.015, m),
            antinodes(t, 1, 0, 0.37, 0.015, m),
            antinodes(b, 0, 1, 0.0, 0.015, m),
            antinodes(b, 1, 0, 0.37, 0.015, m),
        )
    )
