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
"""
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


def main():
    out = sys.argv[1]
    nfins = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    facets, vol = build(nfins)
    with open(out, "w") as f:
        f.write("solid heatsink\n")
        for nrm, a, b, c in facets:
            f.write("  facet normal %f %f %f\n    outer loop\n" % nrm)
            for v in (a, b, c):
                f.write("      vertex %f %f %f\n" % v)
            f.write("    endloop\n  endfacet\n")
        f.write("endsolid heatsink\n")
    print(f"facets={len(facets)} fins={nfins} volume_m3={vol:.6e} -> {out}")


if __name__ == "__main__":
    main()
