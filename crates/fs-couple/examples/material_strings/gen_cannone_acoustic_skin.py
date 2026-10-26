"""Write a closed acoustic skin of the il Cannone violin body for the one-way observer (pexfp.299), as the
`acoustic-skin-v1` section of its v4 body profile (README.md, "Finite-distance pressure through a declared skin").

usage (trj or Mac): uv run --python 3.12 --with triangle --with numpy python3 -I gen_cannone_acoustic_skin.py \
    TABLE PROFILE H_MM K_HOLE1_MM K_HOLE2_MM OUT_PROFILE

TABLE is the shell-assembly table (data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026/shell-assembly-h6.tsv)
and PROFILE the v4 profile `material_strings assemble-shell` wrote from it. OUT_PROFILE is PROFILE followed by the skin.
The skin reads only the table's own parts, so every skin vertex is a point of a structural triangle:
- the top and back plates inside the rib attachment ring: the ring's own plate nodes, every few, and the planform
  between them remeshed at H_MM (Shewchuk's Triangle) and lifted onto the plate triangle under each point;
- the ribs from plate midsurface to plate midsurface: the rib strip's own nodes, every few columns and layers, plus the
  ring's plate nodes above and below, which close the half-plate-thickness slot between the rib edge and each plate;
- the two f-holes closed by caps on every one of their boundary nodes; cap interior points lie on the least-squares plane of the
  hole's boundary and move with nothing. A cap carries its hole's share of the cavity's neck flow, K_i / (K_1 + K_2),
  from each hole's potential-flow conductance alone (fs-bem plate_aperture, b224d16a6: 27.022 and 22.330 mm at the
  table's 3.05 mm top); together the holes carry 0.918 of that sum, and how the coupling splits the flow is not known.
Left out, and stated: the plate overhang beyond the ring (2.5-8 mm fins whose two faces move together), the bridge, the
bass bar and soundpost (inside), and the parts the table already omits (neck, blocks, linings, fingerboard, tailpiece).
Checked before writing: every edge used twice in opposite directions, Euler characteristic 2, positive enclosed volume,
every anchored point equal to its vertex, and the table's frame equal to PROFILE's node rows.
"""
import math
import sys
from collections import defaultdict

import numpy as np
import triangle

TABLE, PROFILE, H, K1, K2, OUT = sys.argv[1], sys.argv[2], float(sys.argv[3]) * 1e-3, float(sys.argv[4]), float(sys.argv[5]), sys.argv[6]


def read_table(path):
    lines = iter(open(path, encoding="utf-8").read().splitlines())
    assert next(lines) == "frankensim.shell-assembly.v1"
    parts, ties = [], []
    for line in lines:
        f = line.split("\t")
        if f[0] == "component":
            regions, nodes, tris = int(f[4]), int(f[5]), int(f[6])
            for _ in range(regions):
                assert next(lines).startswith("region\t")
            P = np.array([[float(v) for v in next(lines).split("\t")[1:4]] for _ in range(nodes)])
            T = np.array([[int(v) for v in next(lines).split("\t")[1:4]] for _ in range(tris)])
            parts.append(dict(name=f[1], P=P, T=T))
        elif f[0] == "tie":
            master = [int(v) for v in next(lines).split("\t")[1].split(" ")]
            slave = [int(v) for v in next(lines).split("\t")[1].split(" ")]
            ties.append((f[1], f[2], master, slave))
    return parts, ties


parts, ties = read_table(TABLE)
index = {p["name"]: i for i, p in enumerate(parts)}

# The frame assemble-shell writes (assembly_table.rs): x/y fractions of the bounding box of every node, z in metres.
allP = np.vstack([p["P"] for p in parts])
x0, y0 = allP[:, 0].min(), allP[:, 1].min()
L, W = allP[:, 0].max() - x0, allP[:, 1].max() - y0
profile = open(PROFILE, encoding="utf-8").read()
rows = [r.split() for r in profile.splitlines() if r.startswith("node ")]
assert len(rows) == len(allP), "PROFILE is not this table's profile"
for r, p in zip(rows, allP):
    assert (float(r[1]), float(r[2]), float(r[3])) == ((p[0] - x0) / L, (p[1] - y0) / W, p[2]), "frame differs"
assert "acoustic-skin-v1" not in profile and "interior-fe-v1" not in profile


def corner_anchor(c, v):
    """An anchor exactly at node v of part c: a triangle of c that uses v, barycentric 1 at that corner."""
    T = parts[c]["T"]
    t = int(np.nonzero((T == v).any(axis=1))[0][0])
    b = [0.0, 0.0, 0.0]
    b[list(T[t]).index(v)] = 1.0
    return (c, t, b)


# --- the ring: the rib strip is column-major (node i*(K+1)+k at ring column i, layer k), tied to each plate's ring chain
rib = index["rib"]
ring = {}
for master, slave, mn, sn in ties:
    if slave == "rib":
        side = "top" if master.startswith("top") else "back"
        ring.setdefault(side, []).append((index[master], mn, sn))
n_rib = len(parts[rib]["P"])
cols = sum(len(sn) for _, _, sn in ring["top"])
K = n_rib // cols - 1
assert cols * (K + 1) == n_rib
plate_ring = {}
for side, row in (("top", 0), ("back", K)):
    chain = [(c, m, s) for c, mn, sn in ring[side] for m, s in zip(mn, sn)]
    assert [s for _, _, s in chain] == [i * (K + 1) + row for i in range(cols)], "rib strip layout differs"
    plate_ring[side] = [(c, m) for c, m, _ in chain]

ring_xz = np.array([parts[c]["P"][m][[0, 2]] for c, m in plate_ring["top"]])
spacing = np.median(np.linalg.norm(np.roll(ring_xz, -1, axis=0) - ring_xz, axis=1))
step = max(1, round(H / spacing))
sel_cols = list(range(0, cols, step))
heights = np.linalg.norm(parts[rib]["P"][[i * (K + 1) + K for i in range(cols)]] - parts[rib]["P"][[i * (K + 1) for i in range(cols)]], axis=1)
m_rows = max(1, math.ceil(heights.max() / H))
sel_rows = sorted({round(j * K / m_rows) for j in range(m_rows + 1)})

# --- f-holes: closed loops of free edges of the welded top plate (its zones are separate parts with coincident nodes)
top_parts = [i for i, p in enumerate(parts) if p["name"].startswith("top-")]
back_parts = [i for i, p in enumerate(parts) if p["name"].startswith("back-")]
weld, owner = {}, []
for c in top_parts:
    for v, p in enumerate(parts[c]["P"]):
        key = tuple(p)
        if key not in weld:
            weld[key] = len(owner)
            owner.append((c, v))
edge_use = defaultdict(int)
for c in top_parts:
    P, T = parts[c]["P"], parts[c]["T"]
    for t in T:
        w = [weld[tuple(P[v])] for v in t]
        for a, b in ((w[0], w[1]), (w[1], w[2]), (w[2], w[0])):
            edge_use[(min(a, b), max(a, b))] += 1
adj = defaultdict(list)
for (a, b), k in edge_use.items():
    if k == 1:
        adj[a].append(b)
        adj[b].append(a)
loops, seen = [], set()
for start in adj:
    if start in seen or len(adj[start]) != 2:
        continue
    loop, prev, cur = [start], None, start
    while True:
        seen.add(cur)
        nxt = [v for v in adj[cur] if v != prev]
        if len(adj[cur]) != 2 or not nxt:
            loop = None
            break
        prev, cur = cur, nxt[0]
        if cur == start:
            break
        loop.append(cur)
    if loop:
        loops.append(loop)


def plan(w):
    c, v = owner[w]
    return parts[c]["P"][v][[0, 2]]


def area(poly):
    x, z = poly[:, 0], poly[:, 1]
    return 0.5 * (x @ np.roll(z, -1) - z @ np.roll(x, -1))


def inside(poly, p):
    x, z = poly[:, 0], poly[:, 1]
    x2, z2 = np.roll(x, -1), np.roll(z, -1)
    hit = ((z > p[1]) != (z2 > p[1])) & (p[0] < x + (p[1] - z) * (x2 - x) / (z2 - z + 1e-300))
    return bool(hit.sum() % 2)


holes = [l for l in loops if abs(area(np.array([plan(w) for w in l]))) > 1e-4 and all(inside(ring_xz, plan(w)) for w in l)]
assert len(holes) == 2, f"expected two f-holes inside the ring, found {len(holes)}"
holes.sort(key=lambda l: np.mean([plan(w)[0] for w in l]))  # hole 1 is the bass side (negative x), as top-plate-fholes.tsv


hole_nodes = holes  # every boundary node of each hole: a chord of the outline would cut into the hole itself

# --- skin vertices: positions in the table frame [m] and anchors (part, triangle, barycentric) with weight 1
V, A = [], []


def vertex(p, anchor):
    V.append(np.asarray(p, float))
    A.append(anchor)
    return len(V) - 1


def plate_point(plate, xz):
    """The plate triangle under a plan point: its barycentric lift is the vertex."""
    best = None
    for c in plate:
        P, T = parts[c]["P"], parts[c]["T"]
        a, b, d = P[T[:, 0]][:, [0, 2]], P[T[:, 1]][:, [0, 2]], P[T[:, 2]][:, [0, 2]]
        det = (b[:, 0] - a[:, 0]) * (d[:, 1] - a[:, 1]) - (d[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
        l1 = ((xz[0] - a[:, 0]) * (d[:, 1] - a[:, 1]) - (d[:, 0] - a[:, 0]) * (xz[1] - a[:, 1])) / det
        l2 = ((b[:, 0] - a[:, 0]) * (xz[1] - a[:, 1]) - (xz[0] - a[:, 0]) * (b[:, 1] - a[:, 1])) / det
        lo = np.minimum(np.minimum(1 - l1 - l2, l1), l2)
        t = int(np.argmax(lo))
        if best is None or lo[t] > best[0]:
            best = (lo[t], c, t, l1[t], l2[t])
    lo, c, t, l1, l2 = best
    assert lo > -1e-9, "a skin point lies off its plate"
    bary = np.clip([1 - l1 - l2, l1, l2], 0.0, None)
    bary = bary / bary.sum()
    P = parts[c]["P"][parts[c]["T"][t]]
    return bary @ P, (c, t, [float(bary[0]), float(bary[1]), 1.0 - float(bary[0]) - float(bary[1])])


tris, cap_of = [], []


def plate_cap(side, plate, holes_here, flip):
    boundary = [vertex(parts[c]["P"][m], corner_anchor(c, m)) for c, m in (plate_ring[side][i] for i in sel_cols)]
    pts = [V[i][[0, 2]] for i in boundary]
    segs = [(i, (i + 1) % len(boundary)) for i in range(len(boundary))]
    ids = list(boundary)
    hole_ids, regions = [], [[*np.mean(pts, axis=0), 0, 0]]
    for h, loop in enumerate(holes_here):
        start = len(pts)
        hv = [vertex(parts[owner[w][0]]["P"][owner[w][1]], corner_anchor(*owner[w])) for w in loop]
        hole_ids.append(hv)
        ids += hv
        pts += [V[i][[0, 2]] for i in hv]
        segs += [(start + i, start + (i + 1) % len(loop)) for i in range(len(loop))]
        # a point inside the (S-shaped) hole: the centroid of one triangle of the hole polygon itself
        poly = np.array([V[i][[0, 2]] for i in hv])
        alone = triangle.triangulate(dict(vertices=poly, segments=np.array([(i, (i + 1) % len(poly)) for i in range(len(poly))])), "p")
        probe = alone["vertices"][alone["triangles"][0]].mean(axis=0)
        regions.append([probe[0], probe[1], h + 1, 0])
    out = triangle.triangulate(dict(vertices=np.array(pts), segments=np.array(segs), regions=np.array(regions)),
                               "pq28a%.15fAY" % (math.sqrt(3) / 4 * H * H))  # fixed point: Triangle reads no exponent; Y: no Steiner point on the ring
    assert np.array_equal(out["vertices"][: len(pts)], np.array(pts)), "Triangle moved an input vertex"
    planes = []
    for hv in hole_ids:
        Q = np.array([V[i] for i in hv])
        G = np.c_[np.ones(len(Q)), Q[:, 0], Q[:, 2]]
        planes.append(np.linalg.lstsq(G, Q[:, 1], rcond=None)[0])
    attr = out["triangle_attributes"][:, 0].astype(int)
    new = {}
    for j, xz in enumerate(out["vertices"]):
        if j < len(pts):
            new[j] = ids[j]
            continue
        owners = {attr[t] for t in np.nonzero((out["triangles"] == j).any(axis=1))[0]}
        if 0 in owners:  # on the plate, or on an f-hole outline between two of its kept nodes
            assert len(owners) <= 2
            p, anchor = plate_point(plate, xz)
            new[j] = vertex(p, anchor)
        else:
            assert len(owners) == 1
            a0, a1, a2 = planes[owners.pop() - 1]
            new[j] = vertex([xz[0], a0 + a1 * xz[0] + a2 * xz[1], xz[1]], None)
    for t, g in zip(out["triangles"], attr):
        tri = [new[int(v)] for v in t]
        tris.append(tri[::-1] if flip else tri)
        cap_of.append(g if g > 0 else 0)
    return boundary


# Triangle's counter-clockwise (x, z) triangles face -y, the top's outside (the top plate is at negative y).
top_ring = plate_cap("top", top_parts, hole_nodes, flip=False)
back_ring = plate_cap("back", back_parts, [], flip=True)
caps_end = len(tris)  # both plates and the f-hole caps face outward by construction
band = [top_ring]
for k in sel_rows:
    band.append([vertex(parts[rib]["P"][i * (K + 1) + k], corner_anchor(rib, i * (K + 1) + k)) for i in sel_cols])
band.append(back_ring)
for r in range(len(band) - 1):
    for i in range(len(sel_cols)):
        i2 = (i + 1) % len(sel_cols)
        a, b, c, d = band[r][i], band[r][i2], band[r + 1][i2], band[r + 1][i]
        tris += [[a, b, c], [a, c, d]]
        cap_of += [0, 0]

# --- orientation: propagate the caps' outward faces across shared edges, then check closure and volume
E = defaultdict(list)
for t, tri in enumerate(tris):
    for e in range(3):
        a, b = tri[e], tri[(e + 1) % 3]
        E[(min(a, b), max(a, b))].append(t)
assert all(len(ts) == 2 for ts in E.values()), "the skin is not closed"
queue = list(range(caps_end))
done = set(queue)
while queue:
    t = queue.pop()
    for e in range(3):
        a, b = tris[t][e], tris[t][(e + 1) % 3]
        for u in E[(min(a, b), max(a, b))]:
            if u in done:
                continue
            # a consistent neighbour traverses the shared edge b -> a
            if any(tris[u][f] == a and tris[u][(f + 1) % 3] == b for f in range(3)):
                tris[u] = tris[u][::-1]
            done.add(u)
            queue.append(u)
assert len(done) == len(tris)
for (a, b), ts in E.items():
    dirs = [any(tris[t][f] == a and tris[t][(f + 1) % 3] == b for f in range(3)) for t in ts]
    assert dirs[0] != dirs[1], "the skin's orientation is inconsistent"
Vn = np.array(V)
volume = sum(np.dot(Vn[a], np.cross(Vn[b], Vn[c])) for a, b, c in tris) / 6.0
assert volume > 0.0, "the skin faces inward"
assert len(Vn) - len(E) + len(tris) == 2
assert len(tris) <= 8192
assert min(np.linalg.norm(Vn[a] - Vn[b]) for a, b in E) > 1e-5, "a skin edge is shorter than 10 um"

# every anchored point is its vertex, exactly
for p, anchor in zip(V, A):
    if anchor is not None:
        c, t, b = anchor
        q = np.array(b) @ parts[c]["P"][parts[c]["T"][t]]
        assert np.linalg.norm(q - p) < 1e-12

share1 = K1 / (K1 + K2)
shares = [share1, 1.0 - share1]
out = [profile.rstrip("\n"), "acoustic-skin-v1", f"skin-vertices {len(V)}"]
for p, anchor in zip(V, A):
    out.append(f"skin-vertex {float((p[0] - x0) / L)!r} {float((p[1] - y0) / W)!r} {float(p[2])!r} {0 if anchor is None else 1}")
    if anchor is not None:
        c, t, b = anchor
        out.append(f"skin-anchor 1 {c} {t} {float(b[0])!r} {float(b[1])!r} {float(b[2])!r}")
out.append(f"skin-triangles {len(tris)}")
out += [f"skin-triangle {a} {b} {c}" for a, b, c in tris]
out.append("skin-apertures 2")
for h in (1, 2):
    caps = [t for t, g in enumerate(cap_of) if g == h]
    out.append(f"skin-aperture f-hole-{h} {float(shares[h - 1])!r} {len(caps)}")
    out += [f"skin-cap {t}" for t in caps]
with open(OUT, "x", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")

areas = [0.5 * np.linalg.norm(np.cross(Vn[b] - Vn[a], Vn[c] - Vn[a])) for a, b, c in tris]
edges = [np.linalg.norm(Vn[a] - Vn[b]) for a, b in E]
cap_area = [sum(areas[t] for t, g in enumerate(cap_of) if g == h) * 1e6 for h in (1, 2)]
print(f"skin H={H * 1e3:g} mm: {len(V)} vertices, {len(tris)} panels ({sum(1 for g in cap_of if g)} on the f-hole caps), "
      f"edges {min(edges) * 1e3:.2f}-{max(edges) * 1e3:.2f} mm (median {np.median(edges) * 1e3:.2f}), enclosed {volume * 1e3:.4f} L")
print(f"ring columns {len(sel_cols)} of {cols} (every {step}), rib layers {sel_rows} of 0..{K}, f-hole outlines {[len(h) for h in hole_nodes]} nodes")
print(f"f-hole caps {cap_area[0]:.1f} and {cap_area[1]:.1f} mm2, flow shares {shares[0]:.6f} {shares[1]:.6f}")
