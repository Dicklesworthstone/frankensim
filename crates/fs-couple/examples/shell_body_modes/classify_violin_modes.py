"""classify_violin_modes.py TABLE OUTPUT.txt SHAPES.tsv : apply the identification rules of body-mode-preregistration.tsv (row `identification`) to a shell_body_modes run.

TABLE is the shell-assembly table that was run (for nodal areas), OUTPUT.txt the text the example printed, SHAPES.tsv the shapes file it wrote (director displacement of the top and back
nodes for the first 36 air-coupled modes). For every air-coupled mode with the neck open it prints: frequency, neck share, net cavity volume change, mean top and back displacement along the outward
directors, the odd-in-x fraction of the top displacement (area-weighted on a 6 mm grid about the mid line, mirrored cells only), the share of the top displacement energy at the waist
(z -120 to -40 mm) and the sign relation of the plates; then the rules in the order written there. Numpy only.
"""
import sys

import numpy as np

table_path, out_path, shapes_path = sys.argv[1:4]

# nodal areas of every component from the table
areas = {}
lines = open(table_path, encoding="utf-8").read().splitlines()
i = 1
while i < len(lines):
    f = lines[i].split("\t")
    if f[0] == "component":
        name, nr, nn, nt = f[1], int(f[4]), int(f[5]), int(f[6])
        P = np.array([[float(x) for x in lines[i + 1 + nr + k].split("\t")[1:4]] for k in range(nn)])
        T = np.array([[int(x) for x in lines[i + 1 + nr + nn + k].split("\t")[1:4]] for k in range(nt)])
        a = 0.5 * np.linalg.norm(np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]]), axis=1)
        w = np.zeros(nn)
        for k in range(3):
            np.add.at(w, T[:, k], a / 3.0)
        areas[name] = w
        i += 1 + nr + nn + nt
    else:
        i += 1

# the coupled list with the neck open
text = open(out_path, encoding="utf-8").read().splitlines()
start = next(k for k, l in enumerate(text) if l.startswith("open sound hole"))
rows = []
for l in text[start + 2:]:
    t = l.split()
    if not t or not t[0].isdigit():
        break
    rows.append(dict(rank=int(t[0]), f=float(t[1]), neck=float(t[2]), dv=float(t[3]), tw=float(t[4]), bw=float(t[5]), q=float(t[6]), qrad=float(t[7]), qtot=float(t[8])))

# shapes
sh = [l.split("\t") for l in open(shapes_path, encoding="utf-8").read().splitlines()]
head, body = sh[0], sh[1:]
nm = len(head) - 5
part = [r[0] for r in body]
node = np.array([int(r[1]) for r in body])
xyz = np.array([[float(r[2]), float(r[3]), float(r[4])] for r in body])
W = np.array([[float(v) for v in r[5:]] for r in body])
area = np.array([areas[p][n] for p, n in zip(part, node)])
is_top = np.array([p.startswith("top-") for p in part])
x0 = 0.5 * (xyz[is_top, 0].min() + xyz[is_top, 0].max())
cell = 0.006


def odd_fraction(m):
    xs, zs = xyz[is_top, 0] - x0, xyz[is_top, 2]
    w, a = W[is_top, m], area[is_top]
    ix = np.floor((xs - xs.min()) / cell).astype(int)
    iz = np.floor((zs - zs.min()) / cell).astype(int)
    nx, nz = ix.max() + 1, iz.max() + 1
    num, den = np.zeros((nx, nz)), np.zeros((nx, nz))
    np.add.at(num, (ix, iz), a * w)
    np.add.at(den, (ix, iz), a)
    ok = den > 0
    mean = np.where(ok, num / np.where(ok, den, 1.0), 0.0)
    # mirror about the mid line: cell (i, j) <-> (nx_sym - 1 - i, j) with the grid centred on x = 0
    lo = xs.min()
    c = int(np.floor((0.0 - lo) / cell))  # the cell holding the mid line
    e_odd = e_even = 0.0
    for i in range(nx):
        i2 = 2 * c - i
        if not (0 <= i2 < nx) or i2 <= i:
            continue
        for j in range(nz):
            if ok[i, j] and ok[i2, j]:
                a_ = min(den[i, j], den[i2, j])
                e_odd += a_ * ((mean[i, j] - mean[i2, j]) / 2) ** 2
                e_even += a_ * ((mean[i, j] + mean[i2, j]) / 2) ** 2
    return e_odd / (e_odd + e_even) if e_odd + e_even > 0 else float("nan")


def waist_share(m):
    zs = xyz[is_top, 2]
    w2 = area[is_top] * W[is_top, m] ** 2
    return float(w2[(zs > -0.12) & (zs < -0.04)].sum() / w2.sum())


print("rank  f_Hz  neck_share  net_dV  top_w  back_w  odd_frac(top)  waist_share(top)  sign(top*back)  Q_solid(KV, not a result)")
info = []
for r in rows[:nm]:
    m = r["rank"] - 1
    of, ws = odd_fraction(m), waist_share(m)
    sgn = "same" if r["tw"] * r["bw"] > 0 else "opposite"
    info.append(dict(r, odd=of, waist=ws, sign=sgn))
    print("%4d %7.2f %9.3f %8.2e %8.2e %8.2e %9.3f %12.3f %12s %12.1f" % (r["rank"], r["f"], r["neck"], r["dv"], r["tw"], r["bw"], of, ws, sgn, r["q"]))

# the pre-registered rules
a0 = next((r for r in info if r["neck"] >= 0.2), None)
print()
print("a0  :", "rank %d at %.2f Hz (neck share %.3f)" % (a0["rank"], a0["f"], a0["neck"]) if a0 else "not found")
cbr = next((r for r in info if a0 and r["rank"] > a0["rank"] and r["odd"] >= 0.6 and r["waist"] >= 0.3), None)
print("cbr :", "rank %d at %.2f Hz (odd %.2f, waist share %.2f)" % (cbr["rank"], cbr["f"], cbr["odd"], cbr["waist"]) if cbr else "not found")
dv0 = abs(a0["dv"]) if a0 else float("nan")
b1m = next((r for r in info if cbr and r["rank"] > cbr["rank"] and r["odd"] < 0.4 and r["sign"] == "opposite" and abs(r["dv"]) < 0.2 * dv0), None)
print("b1- :", "rank %d at %.2f Hz (odd %.2f, dV %.2e against a0 %.2e)" % (b1m["rank"], b1m["f"], b1m["odd"], b1m["dv"], dv0) if b1m else "not found")
b1p = next((r for r in info if b1m and r["rank"] > b1m["rank"] and r["odd"] < 0.4 and r["sign"] == "same"), None)
print("b1+ :", "rank %d at %.2f Hz (odd %.2f, sign same)" % (b1p["rank"], b1p["f"], b1p["odd"]) if b1p else "not found")
bands = dict(a0=(270, 290), cbr=(380, 420), b1m=(430, 480), b1p=(510, 580))
for name, r in (("a0", a0), ("cbr", cbr), ("b1m", b1m), ("b1p", b1p)):
    lo, hi = bands[name]
    if r is None:
        print("%-4s NOT FOUND (band %d to %d Hz)" % (name, lo, hi))
    else:
        off = 0.0 if lo <= r["f"] <= hi else (r["f"] - lo if r["f"] < lo else r["f"] - hi)
        print("%-4s %.1f Hz: %s (band %d to %d Hz, signed offset %.1f Hz); Q_solid (Kelvin-Voigt, not a result) %.0f, Q_rad %.0f" % (name, r["f"], "MET" if off == 0 else "MISS", lo, hi, off, r["q"], r["qrad"]))
