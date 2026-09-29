#!/usr/bin/env python3
"""Generate heatsink.stl for the heatsink-fan worked example: ONE closed manifold shell.

Base plate 80 x 60 x 5 mm at the origin; NFINS (default 4) fins 6 x 60 x 20 mm
at x-offsets 8 + 18k mm, all in metres. The comb profile lives in the x-z
plane and is extruded along y. Every facet is an axis-aligned rectangle of the
breakpoint grid (x split at every fin edge) cut on one diagonal, so the shell
has no T-junctions and no slivers; the mesh side (crates/fs-mesh/tests/
comb_prism.rs) builds the identical triangulation in Rust. Facet vertices are
written with 6 decimals, exactly like the values used to build them.

Usage: python3 examples/heatsink-fan/generate_heatsink_stl.py examples/heatsink-fan/heatsink.stl [NFINS]
       python3 examples/heatsink-fan/generate_heatsink_stl.py OUT [NFINS] --chip X0 X1 Y0 Y1
       python3 examples/heatsink-fan/generate_heatsink_stl.py OUT [NFINS] --rotate ZDEG XDEG --shift X Y Z

`--rotate`/`--shift` write the SAME shell rotated (about z, then x) and
translated: the rotation-invariance twin for the mesh corpus (bead
q61wp.46). Every facet is then oblique and no coordinate is dyadic, which is
what real CAD looks like. Rotated vertices are written with round-trip-exact
digits (Python repr). The tracked twin `heatsink-rotated.stl` is
`--rotate 35 21 --shift 0.1 0.2 0.05`; MEASURED 2026-09-02 it meshes with
the exact volume (378 tets, 98 Steiner points; fs-mesh CONTRACT items 18-21
record the kernel and recovery defects that took). A copy rounded to nine
decimals does NOT mesh: every coplanar facet grid becomes a 1e-10-noisy
point cloud. The precision an STL carries is part of its meshability.

`--chip X0 X1 Y0 Y1` (metres) writes the same body with a chip FOOTPRINT on the
base bottom (`build_chip`): every face is cut on one global tensor grid (fin
breakpoints plus X0/X1, times 0/Y0/Y1/BASE_Y) into diagonal-split rectangles,
so the footprint is a union of whole facets that a `(box ...)` assignment
selector picks out exactly, with no T-junctions and no fans. Without `--chip` the output is byte-identical
to before. The tracked `heatsink-chip.stl` is `--chip 0.030 0.050 0.020 0.040`,
a 20 x 20 mm die under the middle fins.
"""
import math
import sys

BASE_X, BASE_Y, BASE_Z = 0.080, 0.060, 0.005
FIN_W, FIN_H = 0.006, 0.020


def build(nfins):
    fin_x = [(0.008 + 0.018 * k, 0.008 + 0.018 * k + FIN_W) for k in range(nfins)]
    xs = [0.0]
    for x0, x1 in fin_x:
        xs += [x0, x1]
    xs.append(BASE_X)
    top = BASE_Z + FIN_H
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

    for i in range(len(xs) - 1):
        x0, x1 = xs[i], xs[i + 1]
        is_fin = i % 2 == 1
        quad([(x0, 0.0, 0.0), (x1, 0.0, 0.0), (x1, BASE_Y, 0.0), (x0, BASE_Y, 0.0)], (0, 0, -1))
        z_top = top if is_fin else BASE_Z
        quad([(x0, 0.0, z_top), (x1, 0.0, z_top), (x1, BASE_Y, z_top), (x0, BASE_Y, z_top)], (0, 0, 1))
        for y, ny in ((0.0, -1), (BASE_Y, 1)):
            quad([(x0, y, 0.0), (x1, y, 0.0), (x1, y, BASE_Z), (x0, y, BASE_Z)], (0, ny, 0))
            if is_fin:
                quad([(x0, y, BASE_Z), (x1, y, BASE_Z), (x1, y, top), (x0, y, top)], (0, ny, 0))
        if is_fin:
            quad([(x0, 0.0, BASE_Z), (x0, BASE_Y, BASE_Z), (x0, BASE_Y, top), (x0, 0.0, top)], (-1, 0, 0))
            quad([(x1, 0.0, BASE_Z), (x1, BASE_Y, BASE_Z), (x1, BASE_Y, top), (x1, 0.0, top)], (1, 0, 0))
    quad([(0.0, 0.0, 0.0), (0.0, BASE_Y, 0.0), (0.0, BASE_Y, BASE_Z), (0.0, 0.0, BASE_Z)], (-1, 0, 0))
    quad([(BASE_X, 0.0, 0.0), (BASE_X, BASE_Y, 0.0), (BASE_X, BASE_Y, BASE_Z), (BASE_X, 0.0, BASE_Z)], (1, 0, 0))

    # Closed-manifold check: every directed edge once, its reverse once.
    edges = {}
    for _, a, b, c in facets:
        for e in ((a, b), (b, c), (c, a)):
            edges[e] = edges.get(e, 0) + 1
    bad = [e for e, k in edges.items() if k != 1 or (e[1], e[0]) not in edges]
    assert not bad, f"non-manifold edges: {bad[:4]}"
    # Outward normals: signed volume (divergence theorem) equals the analytic volume.
    vol = 0.0
    for _, a, b, c in facets:
        vol += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0
    expected = BASE_X * BASE_Y * BASE_Z + nfins * FIN_W * BASE_Y * FIN_H
    assert abs(vol - expected) < 1e-12 and vol > 0, (vol, expected)
    return facets, vol


def build_chip(nfins, chip):
    """`build` with a chip footprint on the base bottom (see module docs).

    Every face is cut on ONE global tensor grid, the fin breakpoints plus the
    chip x-cuts times 0/Y0/Y1/BASE_Y, and every cell is an axis-aligned
    rectangle split on one diagonal, as in `build`. Fanning a face through
    collinear edge vertices instead leaves facets fs-mesh cannot recover.
    """
    cx0, cx1, cy0, cy1 = chip
    fin_x = [(0.008 + 0.018 * k, 0.008 + 0.018 * k + FIN_W) for k in range(nfins)]
    xs = sorted({0.0, BASE_X, cx0, cx1} | {x for pair in fin_x for x in pair})
    ys = [0.0, cy0, cy1, BASE_Y]
    top = BASE_Z + FIN_H
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

    def is_fin(x0, x1):
        mid = 0.5 * (x0 + x1)
        return any(f0 <= mid <= f1 for f0, f1 in fin_x)

    cells = [(xs[i], xs[i + 1], is_fin(xs[i], xs[i + 1])) for i in range(len(xs) - 1)]
    for x0, x1, fin in cells:
        z_top = top if fin else BASE_Z
        for j in range(len(ys) - 1):
            y0, y1 = ys[j], ys[j + 1]
            quad([(x0, y0, 0.0), (x1, y0, 0.0), (x1, y1, 0.0), (x0, y1, 0.0)], (0, 0, -1))
            quad([(x0, y0, z_top), (x1, y0, z_top), (x1, y1, z_top), (x0, y1, z_top)], (0, 0, 1))
        for y, ny in ((0.0, -1), (BASE_Y, 1)):
            quad([(x0, y, 0.0), (x1, y, 0.0), (x1, y, BASE_Z), (x0, y, BASE_Z)], (0, ny, 0))
            if fin:
                quad([(x0, y, BASE_Z), (x1, y, BASE_Z), (x1, y, top), (x0, y, top)], (0, ny, 0))
    # Fin walls: x-planes where a fin cell meets a non-fin cell.
    for k, (x0, x1, fin) in enumerate(cells):
        if not fin:
            continue
        for x, nx, neighbour in ((x0, -1, k - 1), (x1, 1, k + 1)):
            if 0 <= neighbour < len(cells) and cells[neighbour][2]:
                continue
            for j in range(len(ys) - 1):
                y0, y1 = ys[j], ys[j + 1]
                quad([(x, y0, BASE_Z), (x, y1, BASE_Z), (x, y1, top), (x, y0, top)], (nx, 0, 0))
    for x, nx in ((0.0, -1), (BASE_X, 1)):
        for j in range(len(ys) - 1):
            y0, y1 = ys[j], ys[j + 1]
            quad([(x, y0, 0.0), (x, y1, 0.0), (x, y1, BASE_Z), (x, y0, BASE_Z)], (nx, 0, 0))

    edges = {}
    for _, a, b, c in facets:
        for e in ((a, b), (b, c), (c, a)):
            edges[e] = edges.get(e, 0) + 1
    bad = [e for e, k in edges.items() if k != 1 or (e[1], e[0]) not in edges]
    assert not bad, f"non-manifold edges: {bad[:4]}"
    vol = 0.0
    for _, a, b, c in facets:
        vol += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0
    expected = BASE_X * BASE_Y * BASE_Z + nfins * FIN_W * BASE_Y * FIN_H
    assert abs(vol - expected) < 1e-12 and vol > 0, (vol, expected)
    footprint = sum(
        1 for _, a, b, c in facets
        if all(q[2] == 0.0 and cx0 <= q[0] <= cx1 and cy0 <= q[1] <= cy1 for q in (a, b, c))
    )
    assert footprint > 0, "footprint has no whole facets"
    return facets, vol


def transform(rotate, shift):
    """Rotation about z by rotate[0] degrees, then about x by rotate[1], then shift."""
    if rotate is None and shift is None:
        return None
    zr = math.radians(rotate[0]) if rotate else 0.0
    xr = math.radians(rotate[1]) if rotate else 0.0
    sz, cz, sx, cx = math.sin(zr), math.cos(zr), math.sin(xr), math.cos(xr)
    dx, dy, dz = shift if shift else (0.0, 0.0, 0.0)

    def point(p):
        q = (cz * p[0] - sz * p[1], sz * p[0] + cz * p[1], p[2])
        r = (q[0], cx * q[1] - sx * q[2], sx * q[1] + cx * q[2])
        return (r[0] + dx, r[1] + dy, r[2] + dz)

    def vector(v):
        q = (cz * v[0] - sz * v[1], sz * v[0] + cz * v[1], v[2])
        return (q[0], cx * q[1] - sx * q[2], sx * q[1] + cx * q[2])

    return point, vector


def main():
    argv = sys.argv[1:]
    rotate = shift = None
    if "--rotate" in argv:
        i = argv.index("--rotate")
        rotate = (float(argv[i + 1]), float(argv[i + 2]))
        del argv[i:i + 3]
    if "--shift" in argv:
        i = argv.index("--shift")
        shift = (float(argv[i + 1]), float(argv[i + 2]), float(argv[i + 3]))
        del argv[i:i + 4]
    chip = None
    if "--chip" in argv:
        i = argv.index("--chip")
        chip = tuple(float(argv[i + k]) for k in range(1, 5))
        assert 0.0 < chip[0] < chip[1] < BASE_X and 0.0 < chip[2] < chip[3] < BASE_Y, chip
        del argv[i:i + 5]
    out = argv[0]
    nfins = int(argv[1]) if len(argv) > 1 else 4
    facets, vol = build(nfins) if chip is None else build_chip(nfins, chip)
    mapping = transform(rotate, shift)
    fmt = "%f" if mapping is None else "%r"
    with open(out, "w") as f:
        f.write("solid heatsink\n")
        for nrm, a, b, c in facets:
            if mapping is not None:
                point, vector = mapping
                nrm, a, b, c = vector(nrm), point(a), point(b), point(c)
            f.write(("  facet normal %s %s %s\n    outer loop\n" % ((fmt,) * 3)) % nrm)
            for v in (a, b, c):
                f.write(("      vertex %s %s %s\n" % ((fmt,) * 3)) % v)
            f.write("    endloop\n  endfacet\n")
        f.write("endsolid heatsink\n")
    print(f"facets={len(facets)} fins={nfins} volume_m3={vol:.6e} rotate={rotate} shift={shift} -> {out}")


if __name__ == "__main__":
    main()
