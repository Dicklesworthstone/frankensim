"""Compact view of guitar_body_modes logs: the first rows of the open-hole and sealed-hole coupled tables per variant, with
top/back phase and Q, next to the specimen's measured modes and the pre-registered windows.
usage: python3 -I summarize_modes.py LOG [LOG ...]
"""
import re
import sys

MEAS = [("A0", 98, 0.15), ("T1", 200, 0.20), ("T2", 240, 0.20), ("T3", 277, 0.20), ("T4", 317, 0.20), ("A1", 408, 0.25)]


def parse(path):
    text = open(path, encoding="utf-8", errors="replace").read().splitlines()
    meta = {}
    sections = {}
    cur = None
    for line in text:
        if line.startswith("open sound hole"):
            cur = "open"
            sections[cur] = []
        elif line.startswith("sealed sound hole"):
            cur = "sealed"
            sections[cur] = []
        elif line.startswith("in-vacuo shell modes"):
            cur = "vacuo"
            sections[cur] = []
        elif line.startswith("cavity:"):
            meta["cavity"] = line.strip()
            cur = None if cur in ("open", "sealed") and False else cur
        elif line.startswith("reduce:"):
            meta["reduce"] = line.strip()
        elif cur and re.match(r"^\s+\d+\s", line):
            sections[cur].append(line.split("|")[0].split())
    return meta, sections


for path in sys.argv[1:]:
    meta, sec = parse(path)
    print("=== %s" % path)
    for k in ("reduce", "cavity"):
        if k in meta:
            print(meta[k][:200])
    if "vacuo" in sec:
        print("in-vacuo: " + " ".join("%.1f" % float(r[1]) for r in sec["vacuo"][:14]) + " ...")
    for name in ("open", "sealed"):
        if name not in sec:
            continue
        print("%s hole, first coupled modes (f Hz, neck share, top/back phase, Q_solid):" % name)
        for r in sec[name][:14]:
            f, neck, tw, bw, q = float(r[1]), float(r[2]), float(r[4]), float(r[5]), "%s (rad %s, total %s)" % (r[6], r[7], r[8])
            phase = "in-phase" if tw * bw > 0 else "antiphase"
            print("   %8.2f  neck %.2f  %-9s  Q %s" % (f, neck, phase, q))
    if "open" in sec:
        fs = [float(r[1]) for r in sec["open"]]
        print("windows (pre-registered, relative to the specimen's own measured value):")
        for name, f0, tol in MEAS:
            lo, hi = f0 * (1 - tol), f0 * (1 + tol)
            hit = [f for f in fs if lo <= f <= hi]
            print("   %-2s measured %3d Hz window %5.1f-%5.1f : %s" % (name, f0, lo, hi, ", ".join("%.1f" % f for f in hit) if hit else "NO COMPUTED MODE IN WINDOW"))
