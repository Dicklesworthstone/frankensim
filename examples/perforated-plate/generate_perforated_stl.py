#!/usr/bin/env python3
"""Generate a perforated plate: ONE closed shell with an N x N grid of square holes.

The mesher's scaling corpus body (bead q61wp.46). Plate side 6 mm * (N + 1),
thickness 2 mm, square 3 mm holes on a 6 mm pitch, all in metres. Every face
is cut on one global tensor grid into diagonal-split rectangles (no
T-junctions, no fans), as in examples/plate-hole. Facet count grows as
~21.7 N^2, so N = 20 gives 8652 facets.

`--rotate` writes the same shell rotated 35 deg about z, then 21 deg about x,
and shifted by (0.1, 0.2, 0.05) m. Every facet is then oblique, and the
coordinates are written with round-trip-exact digits (Python repr). This is
the oblique-geometry stress case. Output is deterministic, so a project can
pin its source hash and a lane can generate it instead of tracking megabytes.

Usage: python3 examples/perforated-plate/generate_perforated_stl.py N OUT [--rotate]
"""
import math
import sys


def build(n):
    pitch, hole, z_top = 0.006, 0.003, 0.002
    side = pitch * (n + 1)
    xs = [0.0]
    for k in range(n):
        c = pitch * (k + 1)
        xs += [c - hole / 2, c + hole / 2]
    xs.append(side)
    holes = {(2 * k + 1, 2 * j + 1) for k in range(n) for j in range(n)}
    facets = []

    def quad(corners, outward):
        a, b, c, d = corners
        u = [b[i] - a[i] for i in range(3)]
        v = [c[i] - a[i] for i in range(3)]
        nrm = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
        if sum(nrm[i] * outward[i] for i in range(3)) < 0:
            a, b, c, d = a, d, c, b
        facets.append((tuple(outward), a, b, c))
        facets.append((tuple(outward), a, c, d))

    cells = len(xs) - 1
    for i in range(cells):
        for j in range(cells):
            if (i, j) in holes:
                continue
            x0, x1, y0, y1 = xs[i], xs[i + 1], xs[j], xs[j + 1]
            quad([(x0, y0, 0), (x1, y0, 0), (x1, y1, 0), (x0, y1, 0)], (0, 0, -1))
            quad([(x0, y0, z_top), (x1, y0, z_top), (x1, y1, z_top), (x0, y1, z_top)], (0, 0, 1))
    for i in range(cells):
        a0, a1 = xs[i], xs[i + 1]
        for y, ny in ((0.0, -1), (side, 1)):
            quad([(a0, y, 0), (a1, y, 0), (a1, y, z_top), (a0, y, z_top)], (0, ny, 0))
        for x, nx in ((0.0, -1), (side, 1)):
            quad([(x, a0, 0), (x, a1, 0), (x, a1, z_top), (x, a0, z_top)], (nx, 0, 0))
    for (i, j) in sorted(holes):
        hx0, hx1, hy0, hy1 = xs[i], xs[i + 1], xs[j], xs[j + 1]
        for x, s in ((hx0, 1), (hx1, -1)):
            quad([(x, hy0, 0), (x, hy1, 0), (x, hy1, z_top), (x, hy0, z_top)], (s, 0, 0))
        for y, s in ((hy0, 1), (hy1, -1)):
            quad([(hx0, y, 0), (hx1, y, 0), (hx1, y, z_top), (hx0, y, z_top)], (0, s, 0))

    edges = {}
    for _, a, b, c in facets:
        for e in ((a, b), (b, c), (c, a)):
            edges[e] = edges.get(e, 0) + 1
    assert all(k == 1 and (e[1], e[0]) in edges for e, k in edges.items()), "not closed"
    return facets


def rotate(facets):
    cz, sz = math.cos(math.radians(35)), math.sin(math.radians(35))
    cx, sx = math.cos(math.radians(21)), math.sin(math.radians(21))

    def turn(p):
        x, y, z = p
        x, y = cz * x - sz * y, sz * x + cz * y
        y, z = cx * y - sx * z, sx * y + cx * z
        return (x, y, z)

    def place(p):
        x, y, z = turn(p)
        return (x + 0.1, y + 0.2, z + 0.05)

    return [(turn(n), place(a), place(b), place(c)) for n, a, b, c in facets]


def main():
    args = [a for a in sys.argv[1:] if a != "--rotate"]
    oblique = "--rotate" in sys.argv[1:]
    n, out = int(args[0]), args[1]
    facets = build(n)
    if oblique:
        facets = rotate(facets)
    fmt = "%r" if oblique else "%f"
    with open(out, "w") as f:
        f.write("solid perforated\n")
        for nrm, a, b, c in facets:
            f.write(("  facet normal " + " ".join([fmt] * 3) + "\n    outer loop\n") % nrm)
            for v in (a, b, c):
                f.write(("      vertex " + " ".join([fmt] * 3) + "\n") % v)
            f.write("    endloop\n  endfacet\n")
        f.write("endsolid perforated\n")
    print(f"n={n} facets={len(facets)} rotate={oblique} -> {out}")


if __name__ == "__main__":
    main()
