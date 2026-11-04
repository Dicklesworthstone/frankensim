"""cello_fhole_ring_check.py SPEC_DIR...: clearance of each f-hole raster (top-plate-fholes.tsv run ends) from the rib ring
(ribs-sections.tsv: centre + mean of r_min and r_max at each angle, 3-point mean over theta, sign +1 for z, as the
generator builds its ring before its own inward pull). Positive = inside the ring. Read-only: a hole that crosses the ring
puts a wing over the rib glue line, which the shell generator does not detect."""
import math
import sys

import numpy as np


def rows(path):
    lines = open(path, encoding="utf-8").read().splitlines()
    return np.array([[float(x) for x in l.split("\t")] for l in lines[2:] if l.strip()])


def seg_dist(poly, pts):
    P = np.vstack([poly, poly[:1]])
    d = np.full(len(pts), np.inf)
    for a, b in zip(P[:-1], P[1:]):
        ab = b - a
        t = np.clip(((pts - a) @ ab) / float(ab @ ab), 0, 1)
        d = np.minimum(d, np.hypot(*(pts - (a + t[:, None] * ab)).T))
    return d


def inside(poly, pts):
    out = np.zeros(len(pts), bool)
    n = len(poly)
    for i in range(n):
        (x1, z1), (x2, z2) = poly[i], poly[(i + 1) % n]
        if z1 == z2:
            continue
        m = (z1 > pts[:, 1]) != (z2 > pts[:, 1])
        out ^= m & (pts[:, 0] < (x2 - x1) * (pts[:, 1] - z1) / (z2 - z1) + x1)
    return out


for d in sys.argv[1:]:
    rs = rows(d + "/ribs-sections.tsv")
    th = np.radians(rs[:, 0])
    rmid = 0.5 * (rs[:, 4] + rs[:, 5])
    rmid = (np.roll(rmid, 1) + rmid + np.roll(rmid, -1)) / 3.0
    ring = np.c_[rs[:, 1] + rmid * np.cos(th), rs[:, 2] + rmid * np.sin(th)]
    fh = rows(d + "/top-plate-fholes.tsv")
    for h in (1, 2):
        r = fh[fh[:, 0] == h]
        pts = np.vstack([np.c_[r[:, 2], r[:, 1]], np.c_[r[:, 3], r[:, 1]]])
        dist = seg_dist(ring, pts) * np.where(inside(ring, pts), 1.0, -1.0)
        k = int(np.argmin(dist))
        print("%s hole %d: %d run ends, %d outside the ring; min clearance %.2f mm at (%.2f, %.2f)" % (
            d.rstrip("/").rsplit("/", 1)[-1], h, len(pts), int((dist < 0).sum()), dist[k], *pts[k]))
