#!/usr/bin/env python3
"""Generate plate-hole.stl: ONE closed genus-1 shell, a plate with a square through-hole.

Plate 60 x 40 x 4 mm at the origin with a 12 x 12 mm through-hole centred at
(30, 20) mm, all in metres. Every face is cut on one global tensor grid
(x = 0/24/36/60 mm, y = 0/14/26/40 mm) into axis-aligned rectangles split on
one diagonal, as in examples/heatsink-fan/generate_heatsink_stl.py: no
T-junctions, no fans, no slivers. The hole walls are four rectangles whose
edges are exactly the top and bottom cell edges around the hole, so the
shell is closed and the through-opening is a genuine handle (one inner
boundary loop on each face), which is what the corpus body exists to test.

The generator asserts every edge is used exactly once in each direction
(closed, consistently oriented) and the analytic volume
(60*40 - 12*12) * 4 mm^3 = 9.024 cm^3 and exterior area
2*(2400 - 144) + 4*(200 + 48) mm^2 = 5504 mm^2.

Usage: python3 examples/plate-hole/generate_plate_hole_stl.py examples/plate-hole/plate-hole.stl
"""
import sys

XS = [0.0, 0.024, 0.036, 0.060]
YS = [0.0, 0.014, 0.026, 0.040]
Z = 0.004
HOLE = (1, 1)  # grid cell (i, j) removed through the thickness


def build():
    facets = []

    def quad(corners, outward):
        a, b, c, d = corners
        u = [b[i] - a[i] for i in range(3)]
        v = [c[i] - a[i] for i in range(3)]
        n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
        dot = sum(n[i] * outward[i] for i in range(3))
        assert dot != 0.0, "degenerate quad"
        if dot < 0:
            a, b, c, d = a, d, c, b
        facets.append((tuple(outward), a, b, c))
        facets.append((tuple(outward), a, c, d))

    for i in range(3):
        for j in range(3):
            if (i, j) == HOLE:
                continue
            x0, x1, y0, y1 = XS[i], XS[i + 1], YS[j], YS[j + 1]
            quad([(x0, y0, 0.0), (x1, y0, 0.0), (x1, y1, 0.0), (x0, y1, 0.0)], (0, 0, -1))
            quad([(x0, y0, Z), (x1, y0, Z), (x1, y1, Z), (x0, y1, Z)], (0, 0, 1))
    # Outer walls, split at the grid lines so their edges meet the face cells.
    for i in range(3):
        x0, x1 = XS[i], XS[i + 1]
        for y, ny in ((YS[0], -1), (YS[-1], 1)):
            quad([(x0, y, 0.0), (x1, y, 0.0), (x1, y, Z), (x0, y, Z)], (0, ny, 0))
    for j in range(3):
        y0, y1 = YS[j], YS[j + 1]
        for x, nx in ((XS[0], -1), (XS[-1], 1)):
            quad([(x, y0, 0.0), (x, y1, 0.0), (x, y1, Z), (x, y0, Z)], (nx, 0, 0))
    # Hole walls: outward normals point from the solid into the opening.
    hx0, hx1 = XS[HOLE[0]], XS[HOLE[0] + 1]
    hy0, hy1 = YS[HOLE[1]], YS[HOLE[1] + 1]
    for x, nx in ((hx0, 1), (hx1, -1)):
        quad([(x, hy0, 0.0), (x, hy1, 0.0), (x, hy1, Z), (x, hy0, Z)], (nx, 0, 0))
    for y, ny in ((hy0, 1), (hy1, -1)):
        quad([(hx0, y, 0.0), (hx1, y, 0.0), (hx1, y, Z), (hx0, y, Z)], (0, ny, 0))

    edges = {}
    for _, a, b, c in facets:
        for e in ((a, b), (b, c), (c, a)):
            edges[e] = edges.get(e, 0) + 1
    bad = [e for e, k in edges.items() if k != 1 or (e[1], e[0]) not in edges]
    assert not bad, f"non-manifold edges: {bad[:4]}"
    vol = area = 0.0
    for _, a, b, c in facets:
        vol += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0
        u = [b[i] - a[i] for i in range(3)]
        v = [c[i] - a[i] for i in range(3)]
        n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
        area += 0.5 * sum(x * x for x in n) ** 0.5
    assert abs(vol - 9.024e-6) < 1e-15 and vol > 0, vol
    assert abs(area - 5.504e-3) < 1e-12, area
    # Euler characteristic of a closed genus-1 surface is 0.
    verts = {p for _, a, b, c in facets for p in (a, b, c)}
    assert len(verts) - len(edges) // 2 + len(facets) == 0, "not a single handle"
    return facets, vol, area


def main():
    out = sys.argv[1]
    facets, vol, area = build()
    with open(out, "w") as f:
        f.write("solid plate_hole\n")
        for nrm, a, b, c in facets:
            f.write("  facet normal %f %f %f\n    outer loop\n" % nrm)
            for v in (a, b, c):
                f.write("      vertex %f %f %f\n" % v)
            f.write("    endloop\n  endfacet\n")
        f.write("endsolid plate_hole\n")
    print(f"facets={len(facets)} volume_m3={vol:.6e} area_m2={area:.6e} -> {out}")


if __name__ == "__main__":
    main()
