"""Write a closed acoustic skin of the 1964 Esteso guitar body (pexfp.377) as the `acoustic-skin-v1` section of its v4 body
profile (material_strings README, "Finite-distance pressure through a declared skin"), for the two-way `body-radiation
skin` load with `cavity exterior-neck`. Adapted from gen_cannone_acoustic_skin.py (171da641f) to the guitar table's layout.

usage (trj): uv run --python 3.12 --with triangle --with numpy python3 -I gen_esteso_acoustic_skin.py \
    TABLE PROFILE H_MM OUT_PROFILE

TABLE is the shell-assembly table the renders use (the generator's bridge=1 table, h12-bridge.tsv) and PROFILE the v4
profile `material_strings assemble-shell` wrote from it for one string's drive point; OUT_PROFILE is PROFILE followed by the
skin. The skin is geometry only, so every string's profile takes the same block. It reads only the table's own parts, so
every skin vertex is a point of a structural triangle:
- the top and the back inside the rib attachment ring: the ring's own plate nodes, every few, and the planform between
  them remeshed at H_MM (Shewchuk's Triangle) and lifted onto the plate triangle under each point. The top is flat; the
  back is the table's cylinder over the rib-depth profile, single-valued over the plan;
- the ribs from plate midsurface to plate midsurface: the two rib strips' own nodes on the same columns and a few interior
  layers. The treble strip runs tail to neck and the bass strip neck to tail; their end columns coincide and are used once.
  The strips' edge layers (half a plate thickness from each midsurface, tied to the plate rings) are left out, so the
  plate ring itself closes each rib band;
- the sound hole closed by one cap on every one of its boundary nodes; cap interior points lie on the least-squares plane of
  the hole's boundary (the top's own plane) and move with nothing. The cap carries the whole of the cavity's neck flow.
Left out, and stated: the bridge (a wall standing 8.7 mm off the top's outer face; the skin runs over the top beneath it),
the braces and back bars (inside), and the parts the table already omits (neck and heel, blocks, linings, fingerboard).
Checked before writing: every edge used twice in opposite directions, Euler characteristic 2, positive enclosed volume,
every anchored point equal to its vertex, and the table's frame equal to PROFILE's node rows.
"""
import math
import sys
from collections import defaultdict

import numpy as np
import triangle

TABLE, PROFILE, H, OUT = sys.argv[1], sys.argv[2], float(sys.argv[3]) * 1e-3, sys.argv[4]


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
TOP, BACK, TREBLE, BASS = index["top"], index["back"], index["rib-treble"], index["rib-bass"]

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


# --- the ring: each rib strip is column-major (node i*(K+1)+k at column i, layer k; layer 0 under the top, K over the
# back), tied column by column to a chain of plate nodes. Treble columns 0..n-1 run tail to neck, bass neck to tail.
def chain(plate, rib):
    (m, s), = [(m, s) for a, b, m, s in ties if a == plate and b == rib]
    return m, s


K = None
ring = []  # (rib part, column, top node, back node), once around
for rib, name in ((TREBLE, "rib-treble"), (BASS, "rib-bass")):
    mt, st = chain("top", name)
    mb, sb = chain("back", name)
    k = len(parts[rib]["P"]) // len(st) - 1
    assert len(parts[rib]["P"]) == len(st) * (k + 1) and len(sb) == len(st)
    assert st == [i * (k + 1) for i in range(len(st))] and sb == [i * (k + 1) + k for i in range(len(st))], "rib layout differs"
    K = k if K is None else K
    assert k == K
    cols = list(zip(range(len(st)), mt, mb))
    if rib == BASS:
        cols = cols[1:-1]  # its end columns coincide with the treble strip's
    ring += [(rib, i, t, b) for i, t, b in cols]
mt_treble, _ = chain("top", "rib-treble")
mt_bass, _ = chain("top", "rib-bass")
assert mt_bass[0] == mt_treble[-1] and mt_bass[-1] == mt_treble[0], "the strips do not meet end to end"
assert len({t for _, _, t, _ in ring}) == len(ring), "a plate ring node repeats"
gap = max(np.linalg.norm(parts[rib]["P"][i * (K + 1) + (0 if plate == TOP else K)][:2] - parts[plate]["P"][node][:2])
          for rib, i, t, b in ring for node, plate in ((t, TOP), (b, BACK)))
assert gap < 1e-3, f"a rib column stands {gap * 1e3:.3f} mm off its plate node in plan"

ring_xy = np.array([parts[TOP]["P"][t][:2] for _, _, t, _ in ring])
spacing = np.median(np.linalg.norm(np.roll(ring_xy, -1, axis=0) - ring_xy, axis=1))
step = max(1, round(H / spacing))
sel = list(range(0, len(ring), step))
heights = np.array([parts[TOP]["P"][t][2] - parts[BACK]["P"][b][2] for _, _, t, b in ring])
bands = max(2, math.ceil(heights.max() / H))
sel_rows = sorted({round(j * K / bands) for j in range(1, bands)})  # interior layers only
assert sel_rows[0] > 0 and sel_rows[-1] < K

# --- the sound hole: the top's free-edge loop that does not run along the ring
use = defaultdict(int)
for t in parts[TOP]["T"]:
    for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
        use[(min(a, b), max(a, b))] += 1
adj = defaultdict(list)
for (a, b), k in use.items():
    if k == 1:
        adj[a].append(b)
        adj[b].append(a)
loops, seen = [], set()
for start in adj:
    if start in seen:
        continue
    loop, prev, cur = [start], None, start
    while True:
        seen.add(cur)
        assert len(adj[cur]) == 2, "the top's free edges branch"
        prev, cur = cur, [v for v in adj[cur] if v != prev][0]
        if cur == start:
            break
        loop.append(cur)
    loops.append(loop)
ring_top = {t for _, _, t, _ in ring}
holes = [lp for lp in loops if not ring_top & set(lp)]
outline = [lp for lp in loops if set(lp) == ring_top]
assert len(holes) == 1 and len(outline) == 1 and len(loops) == 2, f"expected the outline and one sound hole, found {len(loops)} loops"
hole = holes[0]

V, A = [], []


def vertex(p, anchor):
    V.append(np.asarray(p, float))
    A.append(anchor)
    return len(V) - 1


def plate_point(c, xy):
    """The plate triangle under a plan point: its barycentric lift is the vertex."""
    P, T = parts[c]["P"], parts[c]["T"]
    a, b, d = P[T[:, 0]][:, :2], P[T[:, 1]][:, :2], P[T[:, 2]][:, :2]
    det = (b[:, 0] - a[:, 0]) * (d[:, 1] - a[:, 1]) - (d[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
    l1 = ((xy[0] - a[:, 0]) * (d[:, 1] - a[:, 1]) - (d[:, 0] - a[:, 0]) * (xy[1] - a[:, 1])) / det
    l2 = ((b[:, 0] - a[:, 0]) * (xy[1] - a[:, 1]) - (xy[0] - a[:, 0]) * (b[:, 1] - a[:, 1])) / det
    lo = np.minimum(np.minimum(1 - l1 - l2, l1), l2)
    t = int(np.argmax(lo))
    assert lo[t] > -1e-9, "a skin point lies off its plate"
    bary = np.clip([1 - l1[t] - l2[t], l1[t], l2[t]], 0.0, None)
    bary = bary / bary.sum()
    return bary @ P[T[t]], (c, t, [float(bary[0]), float(bary[1]), 1.0 - float(bary[0]) - float(bary[1])])


tris, cap_of = [], []


def plate_cap(plate, ring_nodes, hole_loop, flip):
    boundary = [vertex(parts[plate]["P"][n], corner_anchor(plate, n)) for n in ring_nodes]
    pts = [V[i][:2] for i in boundary]
    segs = [(i, (i + 1) % len(boundary)) for i in range(len(boundary))]
    ids = list(boundary)
    regions = [[*np.mean(pts, axis=0), 0, 0]]
    plane = None
    if hole_loop:
        start = len(pts)
        hv = [vertex(parts[plate]["P"][n], corner_anchor(plate, n)) for n in hole_loop]
        ids += hv
        pts += [V[i][:2] for i in hv]
        segs += [(start + i, start + (i + 1) % len(hv)) for i in range(len(hv))]
        regions.append([*np.mean([V[i][:2] for i in hv], axis=0), 1, 0])  # the round hole holds its own centroid
        Q = np.array([V[i] for i in hv])
        plane = np.linalg.lstsq(np.c_[np.ones(len(Q)), Q[:, 0], Q[:, 1]], Q[:, 2], rcond=None)[0]
    out = triangle.triangulate(dict(vertices=np.array(pts), segments=np.array(segs), regions=np.array(regions)),
                               "pq28a%.15fAY" % (math.sqrt(3) / 4 * H * H))  # fixed point: Triangle reads no exponent; Y: no Steiner point on a boundary
    assert np.array_equal(out["vertices"][: len(pts)], np.array(pts)), "Triangle moved an input vertex"
    attr = out["triangle_attributes"][:, 0].astype(int)
    new = {}
    for j, xy in enumerate(out["vertices"]):
        if j < len(pts):
            new[j] = ids[j]
            continue
        owners = {attr[t] for t in np.nonzero((out["triangles"] == j).any(axis=1))[0]}
        if 0 in owners:
            p, anchor = plate_point(plate, xy)
            new[j] = vertex(p, anchor)
        else:
            a0, a1, a2 = plane
            new[j] = vertex([xy[0], xy[1], a0 + a1 * xy[0] + a2 * xy[1]], None)
    for t, g in zip(out["triangles"], attr):
        tri = [new[int(v)] for v in t]
        tris.append(tri[::-1] if flip else tri)
        cap_of.append(int(g))
    return boundary


# Triangle's counter-clockwise (x, y) triangles face +z, the top's outside; the back's outside faces -z.
top_ring = plate_cap(TOP, [ring[i][2] for i in sel], hole, flip=False)
back_ring = plate_cap(BACK, [ring[i][3] for i in sel], None, flip=True)
caps_end = len(tris)
band = [top_ring]
for k in sel_rows:
    band.append([vertex(parts[ring[i][0]]["P"][ring[i][1] * (K + 1) + k], corner_anchor(ring[i][0], ring[i][1] * (K + 1) + k))
                 for i in sel])
band.append(back_ring)
for r in range(len(band) - 1):
    for i in range(len(sel)):
        i2 = (i + 1) % len(sel)
        a, b, c, d = band[r][i], band[r][i2], band[r + 1][i2], band[r + 1][i]
        tris += [[a, b, c], [a, c, d]]
        cap_of += [0, 0]

# --- orientation: propagate the plates' outward faces across shared edges, then check closure and volume
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
for p, anchor in zip(V, A):
    if anchor is not None:
        c, t, b = anchor
        assert np.linalg.norm(np.array(b) @ parts[c]["P"][parts[c]["T"][t]] - p) < 1e-12

out = [profile.rstrip("\n"), "acoustic-skin-v1", f"skin-vertices {len(V)}"]
for p, anchor in zip(V, A):
    out.append(f"skin-vertex {float((p[0] - x0) / L)!r} {float((p[1] - y0) / W)!r} {float(p[2])!r} {0 if anchor is None else 1}")
    if anchor is not None:
        c, t, b = anchor
        out.append(f"skin-anchor 1 {c} {t} {float(b[0])!r} {float(b[1])!r} {float(b[2])!r}")
out.append(f"skin-triangles {len(tris)}")
out += [f"skin-triangle {a} {b} {c}" for a, b, c in tris]
caps = [t for t, g in enumerate(cap_of) if g == 1]
out += ["skin-apertures 1", f"skin-aperture sound-hole 1.0 {len(caps)}"] + [f"skin-cap {t}" for t in caps]
with open(OUT, "x", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")

areas = [0.5 * np.linalg.norm(np.cross(Vn[b] - Vn[a], Vn[c] - Vn[a])) for a, b, c in tris]
edges = [np.linalg.norm(Vn[a] - Vn[b]) for a, b in E]
box = Vn.max(axis=0) - Vn.min(axis=0)
print(f"skin H={H * 1e3:g} mm: {len(V)} vertices, {len(tris)} panels ({len(caps)} on the sound-hole cap, "
      f"{sum(areas[t] for t in caps) * 1e6:.0f} mm2), edges {min(edges) * 1e3:.2f}-{max(edges) * 1e3:.2f} mm "
      f"(median {np.median(edges) * 1e3:.2f}), largest panel {max(areas) * 1e6:.0f} mm2, enclosed {volume * 1e3:.4f} L")
print(f"rib columns stand at most {gap * 1e3:.3f} mm off their plate nodes in plan")
print(f"ring columns {len(sel)} of {len(ring)} (every {step}), rib layers {sel_rows} of 0..{K}, hole outline {len(hole)} nodes; "
      f"box {box[0] * 1e3:.1f} x {box[1] * 1e3:.1f} x {box[2] * 1e3:.1f} mm")
