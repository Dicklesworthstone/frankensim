"""Build the il Cannone violin body shell assembly (top with f-holes, back, ribs, bass bar, soundpost, bridge) from the landed specimen tables, as the text format read by
crates/fs-couple/examples/shell_body_modes.rs (this script lives in crates/fs-couple/examples/shell_body_modes/).

usage (trj or Mac): uv run --python 3.12 --with triangle --with numpy python3 -I gen_violin_shell_assembly.py SPECIMEN_DIR OUT.tsv H_MM [key=value ...]

Every choice that is not a number of the specimen tables is written into the output as a `meta` row, so it is carried with the data. Frame and units of the tables: millimetres,
x across (bass bar at negative x), y vertical (top plate at negative y, back at positive y), z along the body (tail at negative z). The output is in metres in the same frame.
Plate thickness comes from the normal-thickness maps in uniform-thickness zones (the shell owner needs one thickness per shared node, so a thickness step is a separate component
tied to its neighbour along the shared edge chain); the attachment bands of the ribs, bar, post and bridge are zones of their own so that every tie chain lies inside one component.
"""
import math
import sys
from collections import defaultdict, deque

import numpy as np
import triangle

D = sys.argv[1].rstrip("/") + "/"
OUT = sys.argv[2]
H = float(sys.argv[3])  # target element edge length [mm]
opts = dict(a.split("=", 1) for a in sys.argv[4:])
SC = float(opts.get("scale_xz", "1.0"))  # declared sensitivity (u-inplane-scale): plan-view scale of every part about the body centre; heights are not scaled
T_SHIFT = float(opts.get("thickness_shift_mm", "0.0"))  # declared sensitivity (u-segmentation-error): added to every plate zone thickness
PACKS = dict(top=opts.get("top_pack", "spruce"), back=opts.get("back_pack", "maple"), rib=opts.get("rib_pack", "maple"), bar=opts.get("bar_pack", "spruce"),
             post=opts.get("post_pack", "spruce"), bridge=opts.get("bridge_pack", "maple"))
BAR_ON = opts.get("bar", "1") == "1"
POST_ON = opts.get("post", "1") == "1"
BRIDGE_ON = opts.get("bridge", "1") == "1"
EDGE = float(opts.get("edge", "0.45"))  # spacing of the outline and attachment-ring vertices in units of H: a fine ring keeps the ring band narrow
UNIFORM_T = float(opts.get("uniform_t_mm", "0"))  # diagnostic: every plate zone gets this thickness, so that splitting a plate into zones must not change any frequency (a test of the welds)
CAVITY_EXACT = opts.get("cavity_volume", "facet") == "exact"  # "exact": the air volume enclosed by the inner skins (the divergence sum over open midsurfaces is origin dependent, P13 mail 5583); default keeps the old facet sum
ONLY = opts.get("only", "")  # diagnostic: "top" or "back" builds that plate alone (with bar=0 post=0 bridge=0), free in space, to read its modes without the box
FANFIX = opts.get("fanfix", "0") == "1"  # repairs for meshes that leave a vertex with several fans in one component (needed by some of the scaled sensitivity meshes; the base table does not need them)
MIN_PATCH = int(opts.get("min_patch", "24"))  # interior zone patches smaller than this many triangles are merged into their neighbour
ZSTEP = dict(top=float(opts.get("zone_step_top_mm", "0.2")), back=float(opts.get("zone_step_back_mm", "0.4")))
CENTRE = np.array([0.0, -87.0])  # plan-view centre of the body (x, z) for the declared scale

META = []


def meta(k, v):
    META.append((k, v))


def table(name):
    lines = open(D + name, encoding="utf-8").read().splitlines()
    return lines[1].split("\t"), [l.split("\t") for l in lines[2:] if l.strip()]


def ntable(name):
    head, rows = table(name)
    return head, np.array([[float(x) for x in r] for r in rows])


def dimension(part, feature, quantity):
    for r in table("dimensions.tsv")[1]:
        if r[0] == part and r[1] == feature and r[2] == quantity:
            return float(r[3]) * (1000.0 if r[4] == "m" else 1.0)
    raise KeyError((part, feature, quantity))


# ---------------------------------------------------------------- plan-view helpers (points are (x, z) in mm)
def to_plan(P):  # table coordinates -> model plan coordinates (declared scale)
    return CENTRE + (np.asarray(P, float) - CENTRE) * SC


def from_plan(P):  # model plan coordinates -> table coordinates, where the fields live
    return CENTRE + (np.asarray(P, float) - CENTRE) / SC


def resample_closed(P, h):
    Q = np.vstack([P, P[:1]])
    s = np.r_[0, np.cumsum(np.hypot(*(Q[1:] - Q[:-1]).T))]
    t = np.linspace(0, s[-1], max(8, int(round(s[-1] / h))), endpoint=False)
    return np.c_[np.interp(t, s, Q[:, 0]), np.interp(t, s, Q[:, 1])]


def resample_open(P, h):
    s = np.r_[0, np.cumsum(np.hypot(*(P[1:] - P[:-1]).T))]
    t = np.linspace(0, s[-1], max(2, int(round(s[-1] / h)) + 1))
    return np.c_[np.interp(t, s, P[:, 0]), np.interp(t, s, P[:, 1])]


def smooth_closed(P, k):
    for _ in range(k):
        P = 0.25 * np.roll(P, 1, 0) + 0.5 * P + 0.25 * np.roll(P, -1, 0)
    return P


def poly_area(P):
    return 0.5 * float(np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1]))


def inside(poly, p):
    c = False
    n = len(poly)
    for i in range(n):
        x1, z1 = poly[i]
        x2, z2 = poly[(i + 1) % n]
        if (z1 > p[1]) != (z2 > p[1]) and p[0] < (x2 - x1) * (p[1] - z1) / (z2 - z1) + x1:
            c = not c
    return c


def inside_many(poly, pts):
    out = np.zeros(len(pts), bool)
    xs, zs = pts[:, 0], pts[:, 1]
    n = len(poly)
    for i in range(n):
        x1, z1 = poly[i]
        x2, z2 = poly[(i + 1) % n]
        if z1 == z2:
            continue
        m = (z1 > zs) != (z2 > zs)
        xin = (x2 - x1) * (zs - z1) / (z2 - z1) + x1
        out ^= m & (xs < xin)
    return out


def seg_cross(a, b, c, d):
    r, s = b - a, d - c
    den = r[0] * s[1] - r[1] * s[0]
    if abs(den) < 1e-12:
        return None
    t = ((c[0] - a[0]) * s[1] - (c[1] - a[1]) * s[0]) / den
    u = ((c[0] - a[0]) * r[1] - (c[1] - a[1]) * r[0]) / den
    return (t, u) if 1e-6 < t < 1 - 1e-6 and 1e-6 < u < 1 - 1e-6 else None


def insert_crossings(A, B):
    """Insert the crossing points of two open polylines into both."""
    A, B = [np.array(p) for p in A], [np.array(p) for p in B]
    for i in range(len(A) - 1):
        for j in range(len(B) - 1):
            c = seg_cross(A[i], A[i + 1], B[j], B[j + 1])
            if c:
                p = A[i] + c[0] * (A[i + 1] - A[i])
                A.insert(i + 1, p)
                B.insert(j + 1, p.copy())
                return insert_crossings(A, B)
    return np.array(A), np.array(B)


def trace_loops(occ, cell, x0, y0):
    """Boundary loops of a boolean raster occ[row, col] with the cell (col, row) at x0 + col cell, y0 + row cell: closed polylines of the cell corners. Edges run counter-clockwise around
    the occupied cells; at a corner where two occupied cells touch diagonally the walk takes the leftmost turn, which keeps them apart (4-connectivity)."""
    ny, nx = occ.shape
    out = defaultdict(list)
    for j in range(ny):
        for i in range(nx):
            if not occ[j, i]:
                continue
            if j == 0 or not occ[j - 1, i]:
                out[(i, j)].append((i + 1, j))
            if i == nx - 1 or not occ[j, i + 1]:
                out[(i + 1, j)].append((i + 1, j + 1))
            if j == ny - 1 or not occ[j + 1, i]:
                out[(i + 1, j + 1)].append((i, j + 1))
            if i == 0 or not occ[j, i - 1]:
                out[(i, j + 1)].append((i, j))
    used, loops = set(), []
    for a in list(out):
        for b in out[a]:
            if (a, b) in used:
                continue
            loop, cur = [a], (a, b)
            while cur not in used:
                used.add(cur)
                p, q = cur
                loop.append(q)
                opts_ = [r for r in out[q] if (q, r) not in used] or [r for r in out[q] if (q, r) == (q, a) or True]
                if not [r for r in out[q] if (q, r) not in used]:
                    break
                d = (q[0] - p[0], q[1] - p[1])
                r = max([r for r in out[q] if (q, r) not in used], key=lambda r: d[0] * (r[1] - q[1]) - d[1] * (r[0] - q[0]))
                cur = (q, r)
            if len(loop) > 3 and loop[0] == loop[-1]:
                loops.append(np.array([[x0 + i * cell, y0 + j * cell] for i, j in loop[:-1]]))
    return loops


def ellipse_conductance(poly):
    """Acoustic conductance [mm] of a thin-plate aperture of the polygon's shape, from the ellipse with its area and its principal semi-axis from the second moments: K0 = pi a / K(k) = 2 a AGM(1, b/a)
    (the capacitance of an elliptic disc over 4 eps0; a circle gives 2 a). Returns (K0, area, a, b)."""
    P = np.asarray(poly, float)
    x, z = P[:, 0], P[:, 1]
    x1, z1 = np.roll(x, -1), np.roll(z, -1)
    cr = x * z1 - x1 * z
    A = 0.5 * cr.sum()
    cx = ((x + x1) * cr).sum() / (6 * A)
    cz = ((z + z1) * cr).sum() / (6 * A)
    ixx = ((x * x + x * x1 + x1 * x1) * cr).sum() / 12.0 - A * cx * cx
    izz = ((z * z + z * z1 + z1 * z1) * cr).sum() / 12.0 - A * cz * cz
    ixz = ((x * z1 + 2 * x * z + 2 * x1 * z1 + x1 * z) * cr).sum() / 24.0 - A * cx * cz
    cov = np.array([[ixx, ixz], [ixz, izz]]) / A
    lam = np.sort(np.linalg.eigvalsh(cov))[::-1]
    a = 2.0 * math.sqrt(lam[0])  # the semi-axis along the principal direction of the second moments; the other follows from the area: pi a b = A
    b = abs(A) / (math.pi * a)
    g1, g2 = 1.0, b / a
    for _ in range(30):
        g1, g2 = 0.5 * (g1 + g2), math.sqrt(g1 * g2)
    return 2.0 * a * g1, abs(A), a, b


def interior_point(poly):
    """A point well inside a closed polygon: the grid point with the largest distance to the vertices."""
    lo, hi = poly.min(axis=0), poly.max(axis=0)
    gx, gz = np.meshgrid(np.linspace(lo[0], hi[0], 41), np.linspace(lo[1], hi[1], 41))
    g = np.c_[gx.ravel(), gz.ravel()]
    g = g[inside_many(poly, g)]
    d = np.array([np.min(np.hypot(*(poly - p).T)) for p in g])
    return g[int(np.argmax(d))]


def douglas_peucker(P, tol, closed=True):
    if closed:
        k = int(np.argmax(np.hypot(*(P - P[0]).T)))
        a = douglas_peucker(P[:k + 1], tol, False)
        b = douglas_peucker(np.vstack([P[k:], P[:1]]), tol, False)
        return np.vstack([a[:-1], b[:-1]])
    if len(P) < 3:
        return P
    a, b = P[0], P[-1]
    ab = b - a
    L = np.hypot(*ab)
    d = np.hypot(*(P - a).T) if L < 1e-12 else np.abs(ab[0] * (P[:, 1] - a[1]) - ab[1] * (P[:, 0] - a[0])) / L
    k = int(np.argmax(d))
    if d[k] <= tol:
        return np.vstack([a, b])
    return np.vstack([douglas_peucker(P[:k + 1], tol, False)[:-1], douglas_peucker(P[k:], tol, False)])


# ---------------------------------------------------------------- plate fields from the maps tables (table coordinates inside; plan coordinates outside)
class PlateField:
    """3 mm block maps of one plate: midsurface height, outward normal and normal thickness at any plan (x, z)."""

    def __init__(self, prefix, outward):
        _, m = ntable(prefix + "-plate-maps.tsv")
        self.x0, self.z0 = float(m[:, 0].min()), float(m[:, 1].min())
        self.nx = int(round((m[:, 0].max() - self.x0) / 3.0)) + 1
        self.nz = int(round((m[:, 1].max() - self.z0) / 3.0)) + 1
        ix = np.rint((m[:, 0] - self.x0) / 3.0).astype(int)
        iz = np.rint((m[:, 1] - self.z0) / 3.0).astype(int)
        gs = []
        for col in (2, 3, 4):
            g = np.full((self.nx, self.nz), np.nan)
            g[ix, iz] = m[:, col]
            gs.append(self.fill(g))
        self.ymid = self.smooth(0.5 * (gs[0] + gs[1]), 2)
        self.thick = self.smooth(gs[2], 3)
        self.fx, self.fz = np.gradient(self.ymid, 3.0, 3.0)
        self.outward = outward  # -1: toward -y (top plate), +1: toward +y (back plate)

    @staticmethod
    def fill(g):
        g = g.copy()
        for _ in range(60):
            nan = np.isnan(g)
            if not nan.any():
                break
            pad = np.pad(g, 1, constant_values=np.nan)
            nb = np.stack([pad[:-2, 1:-1], pad[2:, 1:-1], pad[1:-1, :-2], pad[1:-1, 2:]])
            cnt = (~np.isnan(nb)).sum(0)
            avg = np.where(cnt > 0, np.nansum(nb, 0) / np.maximum(cnt, 1), np.nan)
            g[nan] = avg[nan]
        return g

    @staticmethod
    def smooth(g, k):
        for _ in range(k):
            pad = np.pad(g, 1, mode="edge")
            g = sum(pad[1 + a:1 + a + g.shape[0], 1 + b:1 + b + g.shape[1]] for a in (-1, 0, 1) for b in (-1, 0, 1)) / 9.0
        return g

    def bil(self, g, P):
        P = from_plan(np.atleast_2d(P))
        u = np.clip((P[:, 0] - self.x0) / 3.0, 0, self.nx - 1.001)
        v = np.clip((P[:, 1] - self.z0) / 3.0, 0, self.nz - 1.001)
        i, j = np.floor(u).astype(int), np.floor(v).astype(int)
        fu, fv = u - i, v - j
        return g[i, j] * (1 - fu) * (1 - fv) + g[i + 1, j] * fu * (1 - fv) + g[i, j + 1] * (1 - fu) * fv + g[i + 1, j + 1] * fu * fv

    def y(self, P):
        return self.bil(self.ymid, P)

    def t(self, P):
        return self.bil(self.thick, P)

    def normal(self, P):
        """Unit normal of y = f(x, z) on the outward side; the plan scale divides the slopes."""
        fx, fz = self.bil(self.fx, P) / SC, self.bil(self.fz, P) / SC
        n = np.stack([-fx * self.outward, np.full_like(fx, self.outward, dtype=float), -fz * self.outward], axis=1)
        return n / np.linalg.norm(n, axis=1, keepdims=True)

    def xyz(self, P):
        P = np.atleast_2d(P)
        return np.c_[P[:, 0], self.y(P), P[:, 1]]


# ---------------------------------------------------------------- components and ties
class Comp:
    def __init__(self, name, pack, t_mm, exterior, cavity):
        self.name, self.pack, self.t = name, pack, t_mm
        self.exterior, self.cavity = exterior, cavity
        self.nodes = []  # x, y, z mm, then the director
        self.tris = []  # i, j, k, region
        self.regions = []

    def add(self, p, d):
        self.nodes.append((float(p[0]), float(p[1]), float(p[2]), float(d[0]), float(d[1]), float(d[2])))
        return len(self.nodes) - 1


components, ties = [], []  # ties: (master name, slave name, side, master node ids, slave node ids)


def orient_tris(comp):
    """Flip each triangle so the facet normal agrees with the mean nodal director."""
    P = np.array([n[:3] for n in comp.nodes])
    Dd = np.array([n[3:] for n in comp.nodes])
    out = []
    for (i, j, k, g) in comp.tris:
        n = np.cross(P[j] - P[i], P[k] - P[i])
        if n @ (Dd[i] + Dd[j] + Dd[k]) < 0:
            j, k = k, j
        out.append((i, j, k, g))
    comp.tris = out


def quad_strip(comp, rows, closed_cols=False, region=0):
    """Triangulate a grid of node ids rows[k][i] (k along the height, i along the length)."""
    nk, ni = len(rows), len(rows[0])
    for k in range(nk - 1):
        for i in range(ni - (0 if closed_cols else 1)):
            i2 = (i + 1) % ni
            a, b, c, e = rows[k][i], rows[k][i2], rows[k + 1][i2], rows[k + 1][i]
            if (i + k) % 2 == 0:
                comp.tris += [(a, b, c, region), (a, c, e, region)]
            else:
                comp.tris += [(a, b, e, region), (b, c, e, region)]


def add_ring_ties(master, slave, side, master_nodes, slave_nodes):
    """A tie along a closed ring: the owner refuses a node twice in one tie and a slave node in two ties, so the ring is cut into two chains that share no node; the two segments between
    the chains are left to their tied end nodes (the linear fields of the two surfaces agree along them)."""
    assert len(master_nodes) == len(slave_nodes)
    if master_nodes[0] == master_nodes[-1] and slave_nodes[0] == slave_nodes[-1]:
        master_nodes, slave_nodes = master_nodes[:-1], slave_nodes[:-1]
    h = len(master_nodes) // 2
    for a_, b_ in ((0, h), (h, len(master_nodes))):
        ties.append((master, slave, side, list(master_nodes[a_:b_]), list(slave_nodes[a_:b_])))


def chains_from_edges(edges):
    """Ordered vertex chains of a set of undirected edges: paths end at degree != 2 vertices, closed loops repeat their first vertex."""
    adj = defaultdict(set)
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    used, out = set(), []

    def walk(a, b):
        chain = [a, b]
        used.add(frozenset((a, b)))
        while len(adj[chain[-1]]) == 2:
            nxt = [v for v in adj[chain[-1]] if frozenset((chain[-1], v)) not in used]
            if not nxt:
                break
            used.add(frozenset((chain[-1], nxt[0])))
            chain.append(nxt[0])
            if nxt[0] == chain[0]:
                break
        return chain

    for a in sorted(adj):
        if len(adj[a]) != 2:
            for b in sorted(adj[a]):
                if frozenset((a, b)) not in used:
                    out.append(walk(a, b))
    for a in sorted(adj):
        for b in sorted(adj[a]):
            if frozenset((a, b)) not in used:
                out.append(walk(a, b))
    return out


# ---------------------------------------------------------------- one plate: mesh, zones, components
def f_hole_rings(h_ring):
    """Closed plan-view rings of the two f-holes from the 0.5 mm raster runs, with one interior point each."""
    _, runs = ntable("top-plate-fholes.tsv")
    rings, pts = [], []
    for hole in (1, 2):
        r = runs[runs[:, 0] == hole]
        cell = 0.5
        x0, z0 = r[:, 2].min() - cell, r[:, 1].min() - cell - 0.25
        nx = int(math.ceil((r[:, 3].max() - x0) / cell)) + 2
        nz = int(round((r[:, 1].max() - z0) / cell)) + 2
        occ = np.zeros((nz, nx), bool)
        for _, z, xa, xb in r:
            j = int(round((z - z0) / cell - 0.5))
            occ[j, int(round((xa - x0) / cell)):int(round((xb - x0) / cell))] = True
        loops = trace_loops(occ, cell, x0, z0)
        loop = max(loops, key=lambda L: abs(poly_area(L)))
        ring = douglas_peucker(loop, 0.3)
        ring = resample_closed(smooth_closed(resample_closed(ring, 1.0), 1), h_ring)
        rings.append(ring)
        cj, ci = np.nonzero(occ)
        k = len(cj) // 2
        pts.append(np.array([x0 + (ci[k] + 0.5) * cell, z0 + (cj[k] + 0.5) * cell]))
    return rings, pts


def build_plate(prefix, outward):
    fld = PlateField(prefix, outward)
    # outline ring from the 1 mm rows (the end rows narrower than 10 mm are the notch and button features, left out and named)
    _, a = ntable(prefix + "-plate-outline.tsv")
    wide = a[(a[:, 2] - a[:, 1]) >= 10.0]
    outline = np.vstack([np.c_[wide[:, 2], wide[:, 0]], np.c_[wide[::-1, 1], wide[::-1, 0]]])
    outline = to_plan(smooth_closed(resample_closed(outline, 3.0), 2))
    outline = resample_closed(outline, EDGE * H)
    meta(prefix + "_outline", "%d vertices from the 1 mm outline rows (%d rows, %d rows narrower than 10 mm at the ends left out), smoothed twice, spaced %.1f mm; plan area %.0f mm2" % (len(outline), len(a), len(a) - len(wide), EDGE * H, abs(poly_area(outline))))
    return fld, outline


# the chains that constrain the mesh of a plate: name -> (points, closed)
def mesh_plate(prefix, fld, outline, ring, extra):
    verts, segs, key = [], set(), {}

    def vid(p):
        k = (round(float(p[0]), 3), round(float(p[1]), 3))
        if k not in key:
            key[k] = len(verts)
            verts.append((float(p[0]), float(p[1])))
        return key[k]

    chains = {}

    def add(name, pts, closed):
        ids = [vid(p) for p in pts]
        if closed:
            ids.append(ids[0])
        for u, v in zip(ids[:-1], ids[1:]):
            if u != v:
                segs.add((min(u, v), max(u, v)))
        chains[name] = ids

    add("outline", outline, True)
    add("ring", ring, True)
    holes = []
    for name, (pts, closed, hole_pt) in extra.items():
        add(name, pts, closed)
        if hole_pt is not None:
            holes.append(hole_pt)
    pslg = dict(vertices=np.array(verts), segments=np.array(sorted(segs)))
    if holes:
        pslg["holes"] = np.array(holes)
    res = triangle.triangulate(pslg, "pq25YYa%.4f" % (0.4330127 * H * H * 1.15))
    V, T = res["vertices"], res["triangles"]
    for i, v in enumerate(verts):
        assert abs(V[i][0] - v[0]) < 1e-9 and abs(V[i][1] - v[1]) < 1e-9, (prefix, i)
    P = V[T]
    ang = []
    for k in range(3):
        e1, e2 = P[:, (k + 1) % 3] - P[:, k], P[:, (k + 2) % 3] - P[:, k]
        ang.append(np.degrees(np.arccos(np.clip((e1 * e2).sum(1) / (np.hypot(*e1.T) * np.hypot(*e2.T) + 1e-300), -1, 1))))
    mn = np.min(ang, axis=0)
    meta("mesh_" + prefix, "%d vertices, %d triangles, minimum angle %.1f deg, mean of the per-triangle minimum angles %.1f deg" % (len(V), len(T), mn.min(), mn.mean()))
    return V, T, chains


def zone_components(prefix, fld, V, T, chains, ring, feature_chains, circles, pack, level_step):
    """Split the triangles of one plate into zones: the ring band (all triangles with a vertex on or outside the rib attachment ring), the feature bands (a vertex on a feature chain or
    inside a feature circle) and thickness-quantised interior zones; connected patches become components, each with one thickness, tied along their shared edges."""
    nT = len(T)
    ring_poly = ring
    on_ring = np.zeros(len(V), bool)
    on_ring[chains["ring"]] = True
    out_ring = ~inside_many(ring_poly, V) | on_ring
    zone = np.full(nT, -1, int)  # -2 ring, -3 feature, >= 0 level
    cen = V[T].mean(axis=1)
    tc = fld.t(cen)
    feat_v = np.zeros(len(V), bool)
    for name in feature_chains:
        feat_v[chains[name]] = True
    in_circle = np.zeros(nT, bool)
    for c, r in circles:
        in_circle |= np.hypot(cen[:, 0] - c[0], cen[:, 1] - c[1]) < r
    is_ring = out_ring[T].any(axis=1)
    is_feat = (feat_v[T].any(axis=1) | in_circle) & ~is_ring
    rest = ~(is_ring | is_feat)
    levels = np.floor((tc - tc[rest].min()) / level_step + 1e-9).astype(int)
    zone[rest] = levels[rest]
    zone[is_feat] = -3
    zone[is_ring] = -2
    # connected patches of equal zone label (edge adjacency)
    edge_tris = defaultdict(list)
    for t_i, (a, b, c) in enumerate(T):
        for u, v in ((a, b), (b, c), (c, a)):
            edge_tris[(min(u, v), max(u, v))].append(t_i)
    vert_tris = defaultdict(list)
    for t_i, tri in enumerate(T):
        for v in tri:
            vert_tris[int(v)].append(t_i)

    def edge_patches():
        patch_ = np.full(nT, -1, int)
        count = 0
        for s_ in range(nT):
            if patch_[s_] >= 0:
                continue
            patch_[s_] = count
            dq = deque([s_])
            while dq:
                t_i = dq.popleft()
                a, b, c = T[t_i]
                for u, v in ((a, b), (b, c), (c, a)):
                    for o in edge_tris[(min(u, v), max(u, v))]:
                        if patch_[o] < 0 and zone[o] == zone[s_]:
                            patch_[o] = count
                            dq.append(o)
            count += 1
        return patch_, count

    # a component must be edge-connected: where two edge-connected pieces of a ring or feature zone touch at one vertex (a pinch), every triangle round that vertex joins the zone
    absorbed = 0
    for _ in range(8):
        patch, npatch = edge_patches()
        changed = False
        for v, ts in vert_tris.items():
            for label in (-2, -3):
                mine = [t_i for t_i in ts if zone[t_i] == label]
                if len({patch[t_i] for t_i in mine}) > 1:
                    for t_i in ts:
                        if zone[t_i] == -5 - label or zone[t_i] in (-2, -3) and zone[t_i] != label:
                            raise AssertionError("pinch between the ring zone and a feature zone at vertex %d" % v)
                        if zone[t_i] != label:
                            zone[t_i] = label
                            absorbed += 1
                            changed = True
        if not changed:
            break
    patch, npatch = edge_patches()
    if absorbed:
        meta("pinch_fix_" + prefix, "%d interior triangles joined the ring or feature zones so that every component is edge-connected" % absorbed)
    # merge interior patches smaller than MIN_PATCH triangles into the neighbouring patch sharing the most edges
    for _ in range(6):
        sizes = np.bincount(patch, minlength=patch.max() + 1)
        changed = False
        for p in np.nonzero((sizes < MIN_PATCH) & (sizes > 0))[0]:
            members = np.nonzero(patch == p)[0]
            if zone[members[0]] < 0:
                continue
            votes = defaultdict(int)
            for t_i in members:
                a, b, c = T[t_i]
                for u, v in ((a, b), (b, c), (c, a)):
                    for o in edge_tris[(min(u, v), max(u, v))]:
                        if patch[o] != p and zone[o] >= 0:
                            votes[patch[o]] += 1
            if votes:
                q = max(votes, key=votes.get)
                patch[members] = q
                zone[members] = zone[np.nonzero(patch == q)[0][0]]
                changed = True
        if not changed:
            break
    if FANFIX:
        # a component's triangles round every vertex must form one fan: where a patch meets itself at a vertex (a loop of the patch round another zone touching itself there), the smaller
        # fans of an interior-zone patch go to the neighbouring patch with which they share the most edges, and the patches are recomputed (a detached piece becomes a patch of its own)
        reassigned = 0
        for _ in range(30):
            patch, npatch = edge_patches()
            moved_any = False
            for v, ts in vert_tris.items():
                by_patch = defaultdict(list)
                for t_i in ts:
                    by_patch[patch[t_i]].append(t_i)
                for p_, mine in by_patch.items():
                    if len(mine) < 2 or zone[mine[0]] < 0 and zone[mine[0]] in (-2, -3):
                        continue
                    par = {t_i: t_i for t_i in mine}

                    def fnd(a_):
                        while par[a_] != a_:
                            par[a_] = par[par[a_]]
                            a_ = par[a_]
                        return a_

                    for x_ in mine:
                        for y_ in mine:
                            if x_ < y_ and len(set(T[x_].tolist()) & set(T[y_].tolist())) >= 2:  # they share an edge, which contains v since both contain v
                                par[fnd(x_)] = fnd(y_)
                    groups = defaultdict(list)
                    for t_i in mine:
                        groups[fnd(t_i)].append(t_i)
                    if len(groups) > 1 and zone[mine[0]] >= 0:
                        keep = max(groups.values(), key=len)
                        for g in groups.values():
                            if g is keep:
                                continue
                            votes = defaultdict(int)
                            for t_i in g:
                                a, b, c = T[t_i]
                                for u, w_ in ((a, b), (b, c), (c, a)):
                                    for o in edge_tris[(min(u, w_), max(u, w_))]:
                                        if patch[o] != p_ and zone[o] >= 0:
                                            votes[o] += 1
                            if votes:
                                tgt = max(votes, key=votes.get)
                                for t_i in g:
                                    zone[t_i] = zone[tgt]
                                reassigned += len(g)
                                moved_any = True
            if not moved_any:
                break
        patch, npatch = edge_patches()
        if reassigned:
            meta("fan_fix_" + prefix, "%d triangles moved to a neighbouring interior zone so that every component has one fan at every vertex" % reassigned)
    if FANFIX:
        # ring and feature zones: wherever their triangles round a vertex form several fans, every triangle round that vertex joins the zone
        extra = 0
        for _ in range(8):
            changed = False
            for v, ts in vert_tris.items():
                for label in (-2, -3):
                    mine = [t_i for t_i in ts if zone[t_i] == label]
                    if len(mine) < 2:
                        continue
                    par = {t_i: t_i for t_i in mine}

                    def fnd2(a_):
                        while par[a_] != a_:
                            par[a_] = par[par[a_]]
                            a_ = par[a_]
                        return a_

                    for x_ in mine:
                        for y_ in mine:
                            if x_ < y_ and len(set(T[x_].tolist()) & set(T[y_].tolist())) >= 2:
                                par[fnd2(x_)] = fnd2(y_)
                    if len({fnd2(t_i) for t_i in mine}) > 1:
                        for t_i in ts:
                            if zone[t_i] != label:
                                if zone[t_i] in (-2, -3):
                                    raise AssertionError("fan fix would take a triangle of another special zone at vertex %d" % v)
                                zone[t_i] = label
                                extra += 1
                                changed = True
            if not changed:
                break
        patch, npatch = edge_patches()
        if extra:
            meta("fan_fix_special_" + prefix, "%d interior triangles joined the ring or feature zones so that their components have one fan at every vertex" % extra)
    ids = sorted(set(patch.tolist()))
    comps, local = [], {}
    P3 = fld.xyz(V)
    # nodal directors = area-weighted vertex averages of the facet normals of the whole plate mesh (so that every facet agrees with its three directors), signed to the outward side
    fn = np.cross(P3[T[:, 1]] - P3[T[:, 0]], P3[T[:, 2]] - P3[T[:, 0]])
    fn *= np.where((fn[:, 1] * fld.outward) < 0, -1.0, 1.0)[:, None]  # outward = +y for the back plate, -y for the top plate
    N3 = np.zeros_like(P3)
    for k in range(3):
        np.add.at(N3, T[:, k], fn)
    N3 /= np.linalg.norm(N3, axis=1, keepdims=True)
    tri_area = 0.5 * np.linalg.norm(np.cross(P3[T[:, 1]] - P3[T[:, 0]], P3[T[:, 2]] - P3[T[:, 0]]), axis=1)
    counters = defaultdict(int)
    for p in ids:
        tris = np.nonzero(patch == p)[0]
        z = zone[tris[0]]
        kind = "ring" if z == -2 else ("band" if z == -3 else "z")
        counters[kind] += 1
        name = "%s-%s%d" % (prefix, kind, counters[kind])
        w = tri_area[tris]
        t_eq = (float(np.cbrt(np.sum(w * tc[tris] ** 3) / w.sum())) if not UNIFORM_T else UNIFORM_T) + T_SHIFT
        comp = Comp(name, pack, t_eq, True, True)
        comp.regions = [(0.0, 0.0, 1.0)]
        comp.zone_kind = kind
        comp.area = float(w.sum())
        used = sorted(set(T[tris].ravel().tolist()))
        loc = {v: comp.add(P3[v], N3[v]) for v in used}
        for t_i in tris:
            a, b, c = T[t_i]
            comp.tris.append((loc[a], loc[b], loc[c], 0))
        orient_tris(comp)
        comps.append(comp)
        local[name] = loc
    # ties between zone components along their shared edge chains (same vertices on both sides, midsurface rigid)
    owner = {}
    for p, comp in zip(ids, comps):
        for t_i in np.nonzero(patch == p)[0]:
            owner[t_i] = comp.name
    pair_edges = defaultdict(list)
    for e, ts in edge_tris.items():
        if len(ts) == 2 and owner[ts[0]] != owner[ts[1]]:
            pair_edges[tuple(sorted((owner[ts[0]], owner[ts[1]])))].append(e)
    # a vertex shared by several components is welded once: the node of the first component (creation order) is the master and every other component's node there a slave (the owner refuses
    # a node that is a slave in one tie and a master in another). For each (master, slave) pair the vertices form chains along their interface edges where consecutive vertices are both
    # in the pair's vertex set; vertices left over (their neighbours have another master, or the two components meet there only) are chained to their nearest leftover vertex: the two
    # chords are the same segment on both sides, so the gap is zero, and every slave node is in exactly one tie
    order = {c.name: k for k, c in enumerate(comps)}
    members = defaultdict(list)
    for c in comps:
        for v in local[c.name]:
            members[v].append(c.name)
    pair_vertices = defaultdict(list)
    for v, ns in members.items():
        if len(ns) > 1:
            m_ = min(ns, key=order.get)
            for x_ in ns:
                if x_ != m_:
                    pair_vertices[(m_, x_)].append(v)
    unwelded = 0
    for (ma, sb), vs in sorted(pair_vertices.items(), key=lambda kv: (order[kv[0][0]], order[kv[0][1]])):
        vset = set(vs)
        es = [e for e in pair_edges.get(tuple(sorted((ma, sb))), []) if e[0] in vset and e[1] in vset]
        covered, chains_ = set(), []
        for ch in chains_from_edges(es):
            if len(ch) > 2 and ch[0] == ch[-1]:
                ch = ch[:-1]  # a closed loop: the owner refuses a node twice in one tie
            if len(ch) >= 2:
                chains_.append(list(ch))
                covered.update(ch)
        left = [v for v in vs if v not in covered]
        while len(left) >= 2:
            v0 = left.pop(0)
            j = min(range(len(left)), key=lambda k: float(np.hypot(*(V[left[k]] - V[v0]))))
            chains_.append([v0, left.pop(j)])
        if left:  # one vertex alone: join the nearest chain of this pair if there is one, else it stays unwelded
            if chains_:
                chains_[-1].append(left[0])
            else:
                unwelded += 1
        for ch in chains_:
            ties.append((ma, sb, 0, [local[ma][u] for u in ch], [local[sb][u] for u in ch]))
    if unwelded:
        meta("unwelded_vertices_" + prefix, "%d vertices where two components meet and no second vertex of the pair exists to make a chain: not welded there" % unwelded)
    meta("zones_" + prefix, "%d components (%s); thickness of each = cube root of the area-weighted mean cube of the smoothed normal-thickness map over its triangles (equal bending stiffness), %s; zone step %.2f mm; %d interface chains tied at the midsurface" % (
        len(comps), ", ".join("%s %.3f mm" % (c.name, c.t) for c in comps), "shifted by %.3f mm" % T_SHIFT if T_SHIFT else "no shift", level_step, sum(1 for t in ties if t[2] == 0 and t[0].startswith(prefix) and t[1].startswith(prefix))))
    return comps, local, patch, owner, P3, N3


# ---------------------------------------------------------------- the body
def main():
    # --- rib attachment ring (the rib midsurface seen from above) from the 2 degree half-plane sections
    _, rs = ntable("ribs-sections.tsv")
    th = np.radians(rs[:, 0])
    rmid = 0.5 * (rs[:, 4] + rs[:, 5])
    rmid = (np.roll(rmid, 1) + rmid + np.roll(rmid, -1)) / 3.0
    outl_top = None
    top_fld, top_outline = build_plate("top", -1)
    back_fld, back_outline = build_plate("back", +1)
    best = None
    for sgn in (1.0, -1.0):
        ring_raw = np.c_[rs[:, 1] + rmid * np.cos(th), rs[:, 2] + sgn * rmid * np.sin(th)]
        n_in = int(inside_many(top_outline, to_plan(ring_raw)).sum())
        if best is None or n_in > best[0]:
            best = (n_in, sgn, ring_raw)
    n_in, sgn, ring_raw = best
    meta("rib_ring", "rib midsurface radius = mean of the section's r_min and r_max (ribs-sections.tsv) at each 2 degrees (3-point mean over theta), angle measured with sign %+d for z; %d of %d points lie inside the top outline" % (sgn, n_in, len(ring_raw)))
    ring = resample_closed(to_plan(ring_raw), EDGE * H)
    # overhang: the plate outline beyond the rib midsurface
    fine_t, fine_b = resample_closed(top_outline, 0.5), resample_closed(back_outline, 0.5)

    def signed_gap(p):
        """(signed distance to the nearer outline, nearest outline point, inward direction): positive when p is inside that outline."""
        best = None
        for poly, fine in ((top_outline, fine_t), (back_outline, fine_b)):
            dd = np.hypot(*(fine - p).T)
            k = int(np.argmin(dd))
            ins = inside(poly, p)
            sd = dd[k] if ins else -dd[k]
            inward = (p - fine[k]) / max(dd[k], 1e-9) * (1.0 if ins else -1.0)
            if best is None or sd < best[0]:
                best = (sd, fine[k], inward)
        return best

    sd0 = np.array([signed_gap(p)[0] for p in ring])
    meta("plate_overhang", "signed distance from the rib attachment ring to the nearer of the top and back outlines (positive inside): mean %.2f mm, minimum %.2f, maximum %.2f; %d of %d ring points start outside an outline" % (sd0.mean(), sd0.min(), sd0.max(), int((sd0 < 0).sum()), len(ring)))
    MIN_BEYOND = 1.6  # mm of plate required beyond the ring (both plates) so that the ring is a mesh chain with well-shaped triangles on both sides
    moved, shift_max = 0, 0.0
    ring0 = ring.copy()
    for i in range(len(ring)):
        for _ in range(8):
            sd, q, inward = signed_gap(ring[i])
            if sd >= MIN_BEYOND - 1e-6:
                break
            ring[i] = ring[i] + inward * (MIN_BEYOND - sd + 0.05)
        shift = float(np.hypot(*(ring[i] - ring0[i])))
        if shift > 1e-9:
            moved += 1
            shift_max = max(shift_max, shift)
    ring = smooth_closed(ring, 1)
    ring = np.array([ring[i] if True else ring0[i] for i in range(len(ring))])
    sd1 = np.array([signed_gap(p)[0] for p in ring])
    meta("ring_moved_inward", "%d of %d ring points (the corner tips of the C-bouts and the ends, where an outline and the rib sections nearly meet) moved inward by up to %.2f mm and the ring smoothed once, so that at least %.2f mm of each plate lies beyond it (minimum now %.2f mm)" % (moved, len(ring), shift_max, MIN_BEYOND, sd1.min()))
    assert sd1.min() > 0.8, "the ring is still too close to an outline"

    # --- feature chains: f-holes, bass bar line, bridge feet (top plate); soundpost feet (both plates)
    extra_top = {}
    hole_rings, hole_pts = f_hole_rings(0.6 * H)
    for k, (hr, hp) in enumerate(zip(hole_rings, hole_pts)):
        extra_top["fhole%d" % (k + 1)] = (to_plan(hr), True, to_plan(hp))
    meta("f_holes", "plan outline of each hole traced from the 0.5 mm raster runs (top-plate-fholes.tsv), simplified to 0.3 mm, spaced %.1f mm; polygon areas %s mm2 (table areas 409 and 313)" % (0.6 * H, ", ".join("%.0f" % abs(poly_area(r)) for r in hole_rings)))
    feature_top, circles_top, circles_back, extra_back = [], [], [], {}
    bar = bridge = post = None
    bar_pts = None
    if BAR_ON:
        _, bs = ntable("bass-bar-sections.tsv")  # z, area, cx, cy, xmin, xmax, ymin, ymax
        hb = bs[:, 7] - bs[:, 6]
        wb = bs[:, 1] / hb
        keep = hb > 2.0
        bar = dict(z=bs[keep, 0], h=hb[keep], w=float(np.median(wb[keep])))
        bar_pts = resample_open(to_plan(np.c_[bs[keep, 2], bs[keep, 0]]), 0.8 * H)
        meta("bass_bar", "centre line from the section centroids (%d sections, z %.0f to %.0f mm, section height above 2 mm), wall thickness = median of area/height = %.2f mm, wall height from the section y extents (%.1f to %.1f mm), rectangular cross-section, grain along the line" % (keep.sum(), bar["z"].min(), bar["z"].max(), bar["w"], hb[keep].min(), hb[keep].max()))
    if BRIDGE_ON:
        _, sil = ntable("bridge-silhouette.tsv")  # y, x_start, x_end
        _, bsec = ntable("bridge-sections.tsv")
        _, pv = table("part-volumes.tsv")
        vol_bridge = [float(r[1]) for r in pv if r[0] == "bridge"][0]
        zc = float(np.average(bsec[:, 3], weights=bsec[:, 1]))
        ybot = sil[:, 0].max()
        bottom = sil[sil[:, 0] > ybot - 3.0]  # the two legs: the bottom 3 mm of rows
        # two feet: the x runs of the bottom rows, merged
        runs = sorted((float(r[1]), float(r[2])) for r in bottom)
        feet = []
        for a_, b_ in runs:
            if feet and a_ <= feet[-1][1] + 0.6:
                feet[-1][1] = max(feet[-1][1], b_)
            else:
                feet.append([a_, b_])
        feet = [f for f in feet if f[1] - f[0] > 3.0]
        assert len(feet) == 2, feet
        area_sil = float(np.sum(sil[:, 2] - sil[:, 1]) * 0.25)
        t_bridge = vol_bridge / area_sil
        bridge = dict(sil=sil, zc=zc, feet=feet, t=t_bridge, ybot=ybot)
        meta("bridge", "silhouette from the 0.25 mm raster rows (bridge-silhouette.tsv), plane z = %.2f mm (area-weighted mid plane of the sections), uniform thickness = mesh volume %.0f mm3 / silhouette area %.0f mm2 = %.2f mm, feet = the two legs in the bottom 3 mm of rows (x %.1f to %.1f and %.1f to %.1f mm), grain across x, string notches and the heart cut-out as in the silhouette" % (zc, vol_bridge, area_sil, t_bridge, feet[0][0], feet[0][1], feet[1][0], feet[1][1]))
        foot_pts = []
        for k, (xa, xb) in enumerate(feet):
            pts = to_plan(np.c_[np.linspace(xa, xb, max(3, int(round((xb - xa) / (0.6 * H))) + 1)), np.full(max(3, int(round((xb - xa) / (0.6 * H))) + 1), zc)])
            foot_pts.append(pts)
        if bar_pts is not None:
            bar_pts, foot_pts[0] = insert_crossings(bar_pts, foot_pts[0])
        extra_top["foot1"] = (foot_pts[0], False, None)
        extra_top["foot2"] = (foot_pts[1], False, None)
        feature_top += ["foot1", "foot2"]
    if bar_pts is not None:
        extra_top["bar"] = (bar_pts, False, None)
        feature_top.append("bar")
    if POST_ON:
        _, ps = ntable("soundpost-sections.tsv")  # y, area, cx, cz, ...
        ax = np.polyfit(ps[:, 0], ps[:, 2], 1)
        az = np.polyfit(ps[:, 0], ps[:, 3], 1)
        r_post = 0.5 * dimension("soundpost", "section", "equivalent-diameter-mean")
        t_post = 0.5 * r_post
        r_c = r_post - 0.5 * t_post  # the tube midsurface radius, so the outer radius is the post radius
        foot = {}
        for name, fld in (("top", top_fld), ("back", back_fld)):
            y = -19.0 if name == "top" else 33.0
            for _ in range(6):
                P = to_plan(np.array([[ax[0] * y + ax[1], az[0] * y + az[1]]]))[0]
                y = float(fld.y(P)[0] + (0.5 if name == "top" else -0.5) * float(fld.t(P)[0]))
            foot[name] = (P, y)
        post = dict(ax=ax, az=az, r=r_post, t=t_post, r_c=r_c, foot=foot)
        n_c = 16  # the owner refuses neighbouring facets whose normals differ by more than 25.8 degrees: an octagon (45 degrees) is refused, a 16-gon (22.5) is not
        ang = 2 * math.pi * np.arange(n_c) / n_c
        for name, extra, circles in (("top", extra_top, circles_top), ("back", extra_back, circles_back)):
            c = foot[name][0]
            extra["post"] = (c + r_c * np.c_[np.cos(ang), np.sin(ang)], True, None)
            circles.append((c, r_c))
        feature_top.append("post")
        meta("soundpost", "straight axis through the section centroids (linear fits x(y), z(y)), diameter %.2f mm, modelled as a tube shell of midsurface radius %.2f mm and thickness %.2f mm (outer radius = the post radius; its area is 75 percent of the solid post, which is still far stiffer axially than either plate), grain along the axis, feet = %d-gon rings tied to the inner skins of the plates at the foot positions top %s (y %.1f mm) and back %s (y %.1f mm)" % (2 * r_post, r_c, t_post, n_c, np.round(foot["top"][0], 1), foot["top"][1], np.round(foot["back"][0], 1), foot["back"][1]))

    # --- meshes and zone components
    meshes = {}
    for name, fld, outline, extra, circles, feats, step in (("top", top_fld, top_outline, extra_top, circles_top, feature_top, ZSTEP["top"]),
                                                            ("back", back_fld, back_outline, extra_back, circles_back, ["post"] if POST_ON else [], ZSTEP["back"])):
        if ONLY and name != ONLY:
            continue
        V, T, chains = mesh_plate(name, fld, outline, ring, extra)
        comps, local, patch, owner, P3, N3 = zone_components(name, fld, V, T, chains, ring, feats, circles, PACKS[name], step)
        components.extend(comps)
        meshes[name] = dict(V=V, T=T, chains=chains, comps=comps, local=local, owner=owner, P3=P3, N3=N3, fld=fld)

    def comp_of_chain(name, chain_name, kind=None):
        m = meshes[name]
        ids = m["chains"][chain_name]
        for c in m["comps"]:
            if all(v in m["local"][c.name] for v in ids) and (kind is None or c.zone_kind == kind):
                return c
        raise AssertionError("chain %s of %s is not inside one component: %s" % (chain_name, name, {c.name: (sum(1 for v in ids if v in m["local"][c.name]), c.zone_kind, round(c.area)) for c in m["comps"]}))

    if not ONLY:
        # --- ribs: one closed strip along the attachment ring, from the lower skin of the top to the inner skin of the back
        _, rr = ntable("ribs-sections.tsv")
        t_rib = float(np.percentile(rr[:, 8], 10))
        ctop, cback = comp_of_chain("top", "ring", "ring"), comp_of_chain("back", "ring", "ring")
        ids_t, ids_b = meshes["top"]["chains"]["ring"][:-1], meshes["back"]["chains"]["ring"][:-1]
        n = len(ids_t)
        assert n == len(ids_b)
        Pt, Nt, Pb, Nb = meshes["top"]["P3"], meshes["top"]["N3"], meshes["back"]["P3"], meshes["back"]["N3"]
        top_end = np.array([Pt[v] - 0.5 * ctop.t * Nt[v] for v in ids_t])
        bot_end = np.array([Pb[v] - 0.5 * cback.t * Nb[v] for v in ids_b])
        height = np.linalg.norm(bot_end - top_end, axis=1)
        K = max(2, int(math.ceil(height.max() / (0.9 * H))))
        rib = Comp("rib", PACKS["rib"], t_rib, True, True)
        ccw = poly_area(ring) > 0
        rows = [[] for _ in range(K + 1)]
        tang = []
        nsec = 16
        for i in range(n):
            tv = ring[(i + 1) % n] - ring[i - 1]
            tv = tv / np.hypot(*tv)
            tang.append(tv)
            out = np.array([tv[1], 0.0, -tv[0]]) if ccw else np.array([-tv[1], 0.0, tv[0]])
            for k in range(K + 1):
                rows[k].append(rib.add(top_end[i] + (bot_end[i] - top_end[i]) * k / K, out))
        rib.regions = []
        sec = lambda i: min(nsec - 1, i * nsec // n)
        for s_ in range(nsec):
            m = np.mean([tang[i] for i in range(n) if sec(i) == s_], axis=0)
            m /= np.hypot(*m)
            rib.regions.append((m[0], 0.0, m[1]))
        for k in range(K):
            for i in range(n):
                i2 = (i + 1) % n
                a_, b_, c_, e_ = rows[k][i], rows[k][i2], rows[k + 1][i2], rows[k + 1][i]
                g = sec(i)
                if (i + k) % 2 == 0:
                    rib.tris += [(a_, b_, c_, g), (a_, c_, e_, g)]
                else:
                    rib.tris += [(a_, b_, e_, g), (b_, c_, e_, g)]
        orient_tris(rib)
        components.append(rib)
        lt, lb = meshes["top"]["local"][ctop.name], meshes["back"]["local"][cback.name]
        add_ring_ties(ctop.name, "rib", -1, [lt[v] for v in ids_t], list(rows[0]))
        add_ring_ties(cback.name, "rib", -1, [lb[v] for v in ids_b], list(rows[K]))
        meta("ribs", "one closed strip along the attachment ring (%d vertices), %d layers, thickness %.2f mm = 10th percentile of the section mean widths (a radial section overstates the wall where it is oblique), grain along the local tangent in %d sectors, top edge on the lower skin of the top ring zone (%.2f mm) and bottom edge on the inner skin of the back ring zone (%.2f mm), height %.1f to %.1f mm; linings and corner blocks are not modelled" % (n, K, t_rib, nsec, ctop.t, cback.t, height.min(), height.max()))

    # --- bass bar wall on the lower skin of the top
    if bar is not None:
        m = meshes["top"]
        cb = comp_of_chain("top", "bar", "band")
        ids = m["chains"]["bar"]
        wall = Comp("bass-bar", PACKS["bar"], bar["w"], False, False)
        P = np.array([m["P3"][v] for v in ids])
        spacing = float(np.mean(np.hypot(*(P[1:, [0, 2]] - P[:-1, [0, 2]]).T)))
        hz = lambda z: float(np.interp(z, bar["z"], bar["h"]))
        hmax = max(hz(p[2]) for p in P)
        K = max(2, int(math.ceil(hmax / (0.8 * spacing))))
        wall_rows = [[] for _ in range(K + 1)]
        for j, v in enumerate(ids):
            nrm = m["N3"][v]
            base = m["P3"][v] - 0.5 * cb.t * nrm
            up = -nrm
            u = P[min(j + 1, len(ids) - 1)] - P[max(j - 1, 0)]
            u = u / np.linalg.norm(u)
            d = np.cross(u, up)
            d = d / np.linalg.norm(d)
            h_here = hz(from_plan(P[j][[0, 2]])[1])
            for k in range(K + 1):
                wall_rows[k].append(wall.add(base + up * h_here * k / K, d))
        wall.regions = [(float(P[-1][0] - P[0][0]), 0.0, float(P[-1][2] - P[0][2]))]
        wall.regions[0] = tuple(np.array(wall.regions[0]) / np.linalg.norm(wall.regions[0]))
        quad_strip(wall, wall_rows)
        orient_tris(wall)
        components.append(wall)
        ties.append((cb.name, wall.name, -1, [m["local"][cb.name][v] for v in ids], wall_rows[0]))
        meta("bass_bar_wall", "wall on the lower skin of the top band zone %s (%.3f mm), %d chain vertices, %d layers" % (cb.name, cb.t, len(ids), K))

    # --- soundpost tube between the inner skins
    if post is not None:
        tube = Comp("soundpost", PACKS["post"], post["t"], False, False)
        ends = {}
        for name in ("top", "back"):
            m = meshes[name]
            cp = comp_of_chain(name, "post", "band")
            ids = m["chains"]["post"][:-1]
            ends[name] = (cp, ids, np.array([m["P3"][v] - 0.5 * cp.t * m["N3"][v] for v in ids]))
        axis = np.mean(ends["back"][2], axis=0) - np.mean(ends["top"][2], axis=0)
        axis_u = axis / np.linalg.norm(axis)
        K = max(3, int(math.ceil(np.linalg.norm(axis) / (1.0 * H))))
        tube_rows = [[] for _ in range(K + 1)]
        nc = len(ends["top"][1])
        for i in range(nc):
            for k in range(K + 1):
                p = ends["top"][2][i] + (ends["back"][2][i] - ends["top"][2][i]) * k / K
                radial = p - (np.mean(ends["top"][2], axis=0) + axis * k / K)
                radial = radial - (radial @ axis_u) * axis_u
                tube_rows[k].append(tube.add(p, radial / np.linalg.norm(radial)))
        tube.regions = [tuple(axis_u)]
        quad_strip(tube, [[row[i] for i in range(nc)] for row in tube_rows], closed_cols=True)
        orient_tris(tube)
        components.append(tube)
        for name, row in (("top", tube_rows[0]), ("back", tube_rows[K])):
            cp, ids, _ = ends[name]
            add_ring_ties(cp.name, "soundpost", -1, [meshes[name]["local"][cp.name][v] for v in ids], list(row))
        meta("soundpost_tube", "length %.1f mm between the inner skins, %d layers, %d nodes per ring" % (np.linalg.norm(axis), K, nc))

    # --- bridge shell standing on the outer skin of the top at the two feet
    if bridge is not None:
        m = meshes["top"]
        sil = bridge["sil"]
        cell = 0.25
        xs0 = float(sil[:, 1].min()) - cell
        ys0 = float(sil[:, 0].min()) - 1.5 * cell  # rows are 0.25 mm cell centres; this origin makes the row index an exact integer
        nx = int(math.ceil((sil[:, 2].max() - xs0) / cell)) + 2
        ny = int(round((sil[:, 0].max() - ys0) / cell)) + 2
        occ = np.zeros((ny, nx), bool)
        for y, xa, xb in sil:
            occ[int(round((y - ys0) / cell - 0.5)), int(round((xa - xs0) / cell)):int(round((xb - xs0) / cell))] = True
        loops = trace_loops(occ, cell, xs0, ys0)
        loops = sorted(loops, key=lambda L: -abs(poly_area(L)))
        outer, holes_ = loops[0], [L for L in loops[1:] if abs(poly_area(L)) > 3.0]
        meta('bridge_holes', '%d cut-outs of the silhouette larger than 3 mm2 (areas %s mm2)' % (len(holes_), ', '.join('%.1f' % abs(poly_area(L)) for L in holes_)))
        hs = max(1.0, 0.45 * H)

        def subdivide(P, h):  # keep every corner, split the long edges into pieces no longer than h
            out_ = []
            for a_, b_ in zip(P, np.roll(P, -1, 0)):
                n_ = max(1, int(math.ceil(np.hypot(*(b_ - a_)) / h)))
                out_ += [a_ + (b_ - a_) * q / n_ for q in range(n_)]
            return np.array(out_)

        outer = subdivide(douglas_peucker(outer, 0.3), hs)
        # replace the bottom run of each foot by the exact slave positions (outer skin of the top at the foot chain)
        foot_nodes = []
        for k, name in enumerate(("foot1", "foot2")):
            cf = comp_of_chain("top", name, "band")
            ids = m["chains"][name]
            pos = np.array([m["P3"][v] + 0.5 * cf.t * m["N3"][v] for v in ids])
            foot_nodes.append((cf, ids, pos))
        bridge_poly = outer.copy()

        def leg_bottom(x):  # lowest silhouette point (largest y) of the legs at x
            ys = [r[0] for r in sil if r[0] > bridge["ybot"] - 3.0 and r[1] - 0.13 <= x <= r[2] + 0.13]
            return max(ys) if ys else None

        for (cf, ids, pos), (xa, xb) in zip(foot_nodes, bridge["feet"]):
            xa_, xb_ = (xa - 0.3) * SC, (xb + 0.3) * SC
            idx = []
            for i, p in enumerate(bridge_poly):
                yb = leg_bottom(p[0] / SC)
                if xa_ <= p[0] <= xb_ and yb is not None and abs(p[1] - yb) < 0.7:
                    idx.append(i)
            assert idx and idx == list(range(idx[0], idx[-1] + 1)), ("foot edge not one contiguous run", idx)
            lo, hi = idx[0], idx[-1]
            seq = pos[:, [0, 1]]
            if bridge_poly[lo][0] > bridge_poly[hi][0]:
                seq = seq[::-1]
            bridge_poly = np.vstack([bridge_poly[:lo], seq, bridge_poly[hi + 1:]])

        def simple(poly, name):
            n_ = len(poly)
            d_ = np.hypot(*(np.roll(poly, -1, 0) - poly).T)
            assert d_.min() > 1e-6, (name, "duplicate vertices", int(d_.argmin()))
            for a_ in range(n_):
                for b_ in range(a_ + 2, n_):
                    if a_ == 0 and b_ == n_ - 1:
                        continue
                    assert seg_cross(poly[a_], poly[(a_ + 1) % n_], poly[b_], poly[(b_ + 1) % n_]) is None, (name, "self-intersection", a_, b_)

        simple(bridge_poly, "bridge outline")
        # the plan of the bridge is the (x, y) silhouette; the plate chain nodes carry the exact feet
        verts, segs = [tuple(p) for p in bridge_poly], []
        nv = len(verts)
        segs = [(i, (i + 1) % nv) for i in range(nv)]
        hole_pts_b = []
        for L in holes_:
            Ls = subdivide(douglas_peucker(L, 0.3), hs)
            o = len(verts)
            verts += [tuple(p) for p in Ls]
            segs += [(o + i, o + (i + 1) % len(Ls)) for i in range(len(Ls))]
            hole_pts_b.append(interior_point(Ls))
        # subdivide the outer polygon edges so that no segment is longer than hs (Triangle without Steiner points on segments)
        res_v, res_s = [], []
        pslg = dict(vertices=np.array(verts), segments=np.array(segs))
        if hole_pts_b:
            pslg["holes"] = np.array(hole_pts_b)
        res = triangle.triangulate(pslg, "pq25YYa%.4f" % (0.4330127 * (1.2 * hs) ** 2))
        Vb, Tb = res["vertices"], res["triangles"]
        br = Comp("bridge", PACKS["bridge"], bridge["t"], False, False)
        ids_node = {}
        foot_lookup = {}
        for (cf, ids, pos) in foot_nodes:
            for v, p in zip(ids, pos):
                foot_lookup[(round(float(p[0]), 4), round(float(p[1]), 4))] = p
        for i, (x, y) in enumerate(Vb):
            p3 = foot_lookup.get((round(float(x), 4), round(float(y), 4)))
            ids_node[i] = br.add(p3 if p3 is not None else (x, y, bridge["zc"] if SC == 1.0 else CENTRE[1] + (bridge["zc"] - CENTRE[1]) * SC), (0.0, 0.0, 1.0))
        for (i, j, k) in Tb:
            br.tris.append((ids_node[int(i)], ids_node[int(j)], ids_node[int(k)], 0))
        br.regions = [(1.0, 0.0, 0.0)]
        orient_tris(br)
        components.append(br)
        for (cf, ids, pos) in foot_nodes:
            slave = []
            for p in pos:
                k = foot_lookup.get((round(float(p[0]), 4), round(float(p[1]), 4)))
                found = [i for i in range(len(Vb)) if round(float(Vb[i][0]), 4) == round(float(p[0]), 4) and round(float(Vb[i][1]), 4) == round(float(p[1]), 4)]
                assert found, "foot node missing from the bridge mesh"
                slave.append(ids_node[found[0]])
            ties.append((cf.name, "bridge", +1, [meshes["top"]["local"][cf.name][v] for v in ids], slave))
        meta("bridge_mesh", "%d vertices, %d triangles in the (x, y) plane; the feet nodes sit exactly on the outer skin of the top band zone at the foot chains" % (len(Vb), len(Tb)))

    # --- cavity: volume enclosed by the midsurfaces of the cavity components (divergence theorem), the two f-holes closed by planar patches, minus half the wall thickness
    vol_mid, area_by = 0.0, {}
    for comp in components:
        if not comp.cavity:
            continue
        Pn = np.array([nn[:3] for nn in comp.nodes])
        a_tot = 0.0
        for (i, j, k, g) in comp.tris:
            nrm = np.cross(Pn[j] - Pn[i], Pn[k] - Pn[i])
            a_tot += 0.5 * np.linalg.norm(nrm)
            vol_mid += ((Pn[i] + Pn[j] + Pn[k]) / 3.0) @ (nrm / 2.0) / 3.0
        area_by[comp.name] = a_tot
    hole_area = 0.0
    for hr in hole_rings:
        ah = abs(poly_area(to_plan(hr)))
        yh = float(np.mean(top_fld.y(to_plan(hr))))
        vol_mid += (-yh) * ah / 3.0
        hole_area += ah
    inner = sum(area_by[c.name] * c.t / 2.0 for c in components if c.cavity)
    vol_inner = vol_mid - inner
    if CAVITY_EXACT:
        # exact inner volume: integral over the planform of the rib's inner skin (the rib ring moved in by half the rib thickness) of the distance between the inner skins of the back and the top, on a 0.5 mm grid,
        # minus the shells of the bass bar and the soundpost (the model's solids inside the box); the f-holes are closed by the planar patches, so the volume under them counts
        t_rib_half = 0.5 * float(np.percentile(rr_t, 10)) if False else 0.5 * [c.t for c in components if c.name == "rib"][0]
        cen = ring.mean(axis=0)
        dirs = ring - cen
        ring_in = ring.copy()
        n_ = len(ring)
        for i in range(n_):
            tv = ring[(i + 1) % n_] - ring[i - 1]
            tv = tv / np.hypot(*tv)
            nv = np.array([tv[1], -tv[0]])
            if nv @ dirs[i] > 0:
                nv = -nv  # inward
            ring_in[i] = ring[i] + nv * t_rib_half
        lo, hi = ring_in.min(axis=0), ring_in.max(axis=0)
        gx, gz = np.meshgrid(np.arange(lo[0], hi[0], 0.5) + 0.25, np.arange(lo[1], hi[1], 0.5) + 0.25)
        G = np.c_[gx.ravel(), gz.ravel()]
        G = G[inside_many(ring_in, G)]
        h_in = (back_fld.y(G) - 0.5 * back_fld.t(G)) - (top_fld.y(G) + 0.5 * top_fld.t(G))
        vol_exact = float(np.sum(h_in)) * 0.25
        solids = 0.0
        for c in components:
            if c.name in ("bass-bar", "soundpost"):
                Pn = np.array([nn[:3] for nn in c.nodes])
                a_c = sum(0.5 * np.linalg.norm(np.cross(Pn[j] - Pn[i], Pn[k] - Pn[i])) for (i, j, k, g) in c.tris)
                solids += a_c * c.t
        meta("cavity_volume_facet_sum_L", "%.4f (the divergence-theorem sum minus half thicknesses; open surfaces make it origin dependent, so it is NOT the cavity record when cavity_volume=exact)" % (vol_inner / 1e6))
        meta("cavity_volume_exact", "inner-skin planform integral %.4f L (planform = the rib ring moved in by half the rib thickness %.3f mm, %d grid cells of 0.5 mm, height from the inner skins of back and top) minus the shells of the bass bar and soundpost %.4f L" % (vol_exact / 1e6, t_rib_half, len(G), solids / 1e6))
        vol_inner = vol_exact - solids
    t_holes = float(np.mean(top_fld.t(to_plan(np.vstack(hole_rings)))))
    # the two f-holes as one circular neck of the same conductance: each hole is a thin-plate aperture of its elliptic equivalent (K0 = 2 a AGM(1, b/a)) in series with its through-thickness
    # resistance t / A; the neck radius follows from K = pi r^2 / (t + 2 x 0.8216 r), the library's flanged end correction at both ends, at the top thickness at the holes
    k_parts = []
    for hr in hole_rings:
        k0, a_h, ea, eb = ellipse_conductance(to_plan(hr))
        k_parts.append((1.0 / (1.0 / k0 + t_holes / a_h), k0, a_h, ea, eb))
    k_tot = sum(k[0] for k in k_parts)
    r_eq = (1.6432 * k_tot + math.sqrt((1.6432 * k_tot) ** 2 + 4 * math.pi * k_tot * t_holes)) / (2 * math.pi)
    f_rigid = 343.0 / (2 * math.pi) * math.sqrt(k_tot * 1e-3 / (vol_inner * 1e-9))
    meta("cavity_volume_midsurface_L", "%.4f" % (vol_mid / 1e6))
    meta("cavity_volume_L", ("%.4f (the exact inner-skin volume minus the bass bar and soundpost shells, see cavity_volume_exact)" if CAVITY_EXACT else "%.4f (midsurface facets oriented outward, the f-holes closed by planar patches at their mean height, minus half the thickness of every wall)") % (vol_inner / 1e6))
    meta("neck", "one circular neck of the conductance of the two f-holes: per hole the equivalent ellipse (semi-axes %s mm) has thin-plate conductance K0 %s mm, in series with t/A at the top thickness %.2f mm gives %s mm; total %.2f mm; the neck radius that gives it with the flanged end correction 0.8216 r at both ends is %.2f mm (area %.0f mm2 against the hole area %.0f mm2), length = the top thickness; rigid-wall Helmholtz frequency of this neck and cavity %.0f Hz (c = 343 m/s); the equivalent-ellipse conductance of a curved f-hole is an approximation (UNKNOWN size of the error; named in the pre-registered scope)" % (
        " and ".join("%.1f x %.1f" % (k[3], k[4]) for k in k_parts), " and ".join("%.1f" % k[1] for k in k_parts), t_holes, " and ".join("%.1f" % k[0] for k in k_parts), k_tot, r_eq, math.pi * r_eq ** 2, hole_area, f_rigid))
    meta("scale_and_shift", "plan-view scale %.4f (declared sensitivity u-inplane-scale), thickness shift %.3f mm (declared sensitivity u-segmentation-error)" % (SC, T_SHIFT))
    meta("packs", ", ".join("%s=%s" % kv for kv in sorted(PACKS.items())))
    meta("omitted", "neck and upper block, tail block, tailpiece, strings, chinrest, linings and corner blocks, purfling, varnish, glue compliance, the button and tail-end notch features of the plates (end rows narrower than 10 mm)")
    write(vol_inner, r_eq, t_holes)


def write(vol_inner_mm3, r_eq_mm, neck_len_mm):
    total_nodes = sum(len(c.nodes) for c in components)
    total_tris = sum(len(c.tris) for c in components)
    meta("totals", "%d components, %d nodes, %d triangles, %d ties" % (len(components), total_nodes, total_tris, len(ties)))
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("frankensim.shell-assembly.v1\n")
        for k, v in META:
            f.write("meta\t%s\t%s\n" % (k, v))
        for c in components:
            f.write("component\t%s\t%s\t%.9f\t%d\t%d\t%d\t%d\t%d\n" % (c.name, c.pack, c.t / 1000.0, len(c.regions), len(c.nodes), len(c.tris), int(c.exterior), int(c.cavity)))
            for r in c.regions:
                f.write("region\t%.9f\t%.9f\t%.9f\n" % tuple(r))
            for nn in c.nodes:
                d = np.array(nn[3:])
                d = d / np.linalg.norm(d)
                f.write("node\t%.9f\t%.9f\t%.9f\t%.12f\t%.12f\t%.12f\n" % (nn[0] / 1000.0, nn[1] / 1000.0, nn[2] / 1000.0, d[0], d[1], d[2]))
            for t in c.tris:
                f.write("tri\t%d\t%d\t%d\t%d\n" % t)
        for (m, s_, side, mn, sn) in ties:
            f.write("tie\t%s\t%s\t%d\t%d\t%d\n" % (m, s_, side, len(mn), len(sn)))
            f.write("master\t%s\n" % " ".join(map(str, mn)))
            f.write("slave\t%s\n" % " ".join(map(str, sn)))
        f.write("cavity\t%.9f\t%.9f\t%.9f\n" % (vol_inner_mm3 / 1e9, r_eq_mm / 1000.0, neck_len_mm / 1000.0))
        if ONLY:
            first = next(c for c in components if c.name.startswith(ONLY + "-z"))
            f.write("drive\tdiagnostic centre, along the director\t%s\t0.0\t%.6f\t-0.1\t0\t1\t0\n" % (first.name, first.nodes[0][1]))
        if BRIDGE_ON:
            contact = {(r[0], r[1]): r for r in table("string-contacts.tsv")[1]}  # string, contact, z, x, y, slopes
            for label, key in (("bass string (G) at the bridge, along x", ("g", "bridge")), ("treble string (E) at the bridge, along x", ("e", "bridge"))):
                r = contact[key]
                p = to_plan(np.array([[float(r[3]), float(r[2])]]))[0]
                f.write("drive\t%s\tbridge\t%.6f\t%.6f\t%.6f\t1\t0\t0\n" % (label, p[0] / 1000.0, float(r[4]) / 1000.0, p[1] / 1000.0))
            gz = to_plan(np.array([[0.0, float(contact[("g", "bridge")][2])]]))[0]
            f.write("drive\tbridge centre, along y (into the top)\tbridge\t%.6f\t%.6f\t%.6f\t0\t1\t0\n" % (0.0, -0.0520, gz[1] / 1000.0))
    for k, v in META:
        print(k, "=", v)


if __name__ == "__main__":
    main()
