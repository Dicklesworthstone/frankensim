"""Build the Esteso 1964 closed-body shell assembly (top + brace walls + back + back bars + ribs) from the
specimen tables as the text format read by crates/fs-couple/examples/guitar_body_modes.rs.

usage (trj): uv run --python 3.12 --with triangle --with numpy python3 -I gen_shell_assembly.py SPECIMEN_DIR OUT.tsv H_MM [options]

Every choice that is not a printed or vector-read number is written into the output as a `meta` row so it is
carried with the data (UNKNOWN or authored choices are named as such).
"""
import math
import sys

import numpy as np
import triangle

D = sys.argv[1].rstrip("/") + "/"
OUT = sys.argv[2]
H = float(sys.argv[3])  # target element edge length [mm]
opts = dict(a.split("=", 1) for a in sys.argv[4:])
BACK_ARCH = opts.get("back_arch", "parabola")  # parabola through the three printed rib depths, or flat
BRACE_H_SCALE = float(opts.get("brace_h_scale", "1.0"))  # sensitivity only; 1.0 = printed heights
BACK_BAR_H_SCALE = float(opts.get("back_bar_h_scale", "1.0"))
TOP_T_SCALE = float(opts.get("top_t_scale", "1.0"))
# Constant-strain triangles lock badly in the in-plane bending of a deep wall: a brace or bar needs several element layers through its
# height whatever the plate mesh size is (wall_layers_min, default 1 = the first table's behaviour).
WALL_K_MIN = int(opts.get("wall_layers_min", "1"))
# table format: "shared" = frankensim.shell-assembly.v1 of the shell_body_modes example (with drive records); "guitar" = the first format
FORMAT = opts.get("format", "shared")
TOP_MAT = opts.get("top_pack", "spruce")
BACK_MAT = opts.get("back_pack", "rosewood")
BRACE_MAT = opts.get("brace_pack", "spruce")
BACKBAR_MAT = opts.get("backbar_pack", "cedar")
RIB_MAT = opts.get("rib_pack", "rosewood")

META = []


def meta(k, v):
    META.append((k, v))


def rows(name):
    return [l.split("\t") for l in open(D + name, encoding="utf-8").read().splitlines()[2:] if l.strip()]


def rows_with_header(name):
    lines = open(D + name, encoding="utf-8").read().splitlines()
    return lines[1].split("\t"), [l.split("\t") for l in lines[2:] if l.strip()]


# ---------------------------------------------------------------- thickness (declared equivalents)
def t_equiv(name):
    t = np.array([float(r[3]) for r in rows(name)])
    return float(np.cbrt(np.mean(t**3))), float(np.mean(t)), float(np.std(t)), len(t)


t_top, t_top_mean, t_top_sd, n_top = t_equiv("top-thickness.tsv")
t_back, t_back_mean, t_back_sd, n_back = t_equiv("back-thickness.tsv")
t_rib, t_rib_mean, t_rib_sd, n_rib = t_equiv("side-thickness.tsv")
t_top *= TOP_T_SCALE
meta("top_thickness_mm", "%.4f (cube-root of mean cube of %d printed points; mean %.3f, sd %.3f; scale %.3f)" % (t_top, n_top, t_top_mean, t_top_sd, TOP_T_SCALE))
meta("back_thickness_mm", "%.4f (cube-root of mean cube of %d printed points; mean %.3f, sd %.3f)" % (t_back, n_back, t_back_mean, t_back_sd))
meta("rib_thickness_mm", "%.4f (cube-root of mean cube of %d printed points; mean %.3f, sd %.3f)" % (t_rib, n_rib, t_rib_mean, t_rib_sd))
meta("thickness_policy", "uniform per component (the shell owner requires one thickness per component); the equivalent is the cube root of the mean cube of the printed points, so the plate bending stiffness of the points is kept on average; UNKNOWN between points")

# ---------------------------------------------------------------- planform
outline_main = np.array([[float(r[2]), float(r[3])] for r in rows("top-outline-vector.tsv") if r[0] == "outline"])
outline_second = np.array([[float(r[2]), float(r[3])] for r in rows("top-outline-vector.tsv") if r[0] == "second-thick-path"])
hole_pts = np.array([[float(r[2]), float(r[3])] for r in rows("top-outline-vector.tsv") if r[0] == "soundhole-circle"])
# `outline` runs from (281.70, 164.78) round the tail to the neck end (0.01, -30.00); `second-thick-path` runs from the other side
# of the neck end (-0.45, 30.19) to (282.00, 164.92), which is within 0.4 mm of the start of `outline`: it is the treble upper bout.
# The 60.2 mm gap at the neck end is the printed neck width at the body edge (60 mm): closed here by a straight segment.
gap_join = float(np.hypot(*(outline_main[0] - outline_second[-1])))
gap_neck = float(np.hypot(*(outline_main[-1] - outline_second[0])))
assert gap_join < 1.0, gap_join
outline = np.vstack([outline_main, outline_second])
meta("outline_closure", "outline path (486 points) + second-thick-path (209 points, the treble upper bout; joins the outline start within %.2f mm) + a straight neck-end segment of %.2f mm (the plan prints a neck width of 60 mm at the body edge)" % (gap_join, gap_neck))
# drop zero-length neighbours (the vector reads repeat points)
keep = [0]
for i in range(1, len(outline)):
    if np.hypot(*(outline[i] - outline[keep[-1]])) > 0.05:
        keep.append(i)
outline = outline[keep]


def area(p):
    x, y = p[:, 0], p[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


if area(outline) < 0:
    outline = outline[::-1]


MAX_TURN = float(opts.get("max_turn", "0.30"))  # radians of outline turning per ring vertex (0.30 rad = 17 deg)
MIN_RING_STEP = float(opts.get("min_ring_step", "3.0"))


def resample_closed(p, h):
    """Ring vertices with spacing min(h, MAX_TURN / curvature) along the outline (curvature from the turning of the dense
    vector polyline over a window of about 8 mm on each side)."""
    q = np.vstack([p, p[:1]])
    seg = np.hypot(*(q[1:] - q[:-1]).T)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = s[-1]
    m = len(p)
    # turning angle at each dense vertex
    ang = np.zeros(m)
    for i in range(m):
        a_ = p[i] - p[i - 1]
        b_ = p[(i + 1) % m] - p[i]
        ang[i] = abs(math.atan2(a_[0] * b_[1] - a_[1] * b_[0], a_ @ b_))
    mean_step = total / m
    win = max(1, int(round(8.0 / mean_step)))
    kappa = np.array([sum(ang[(i + k) % m] for k in range(-win, win + 1)) / ((2 * win + 1) * mean_step) for i in range(m)])
    h_local = np.minimum(h, np.maximum(MIN_RING_STEP, MAX_TURN / np.maximum(kappa, 1e-9)))
    # integrate ds / h_local over the dense polyline
    ds = seg
    mid_h = 0.5 * (h_local + np.roll(h_local, -1))
    cum = np.concatenate([[0.0], np.cumsum(ds / mid_h)])
    n = max(8, int(round(cum[-1])))
    targets = np.linspace(0.0, cum[-1], n, endpoint=False)
    u = np.interp(targets, cum, s)
    return np.column_stack([np.interp(u, s, q[:, 0]), np.interp(u, s, q[:, 1])])


# start the ring at the tail (largest x, nearest y = 0)
tail = int(np.argmax(outline[:, 0] - 1e-3 * np.abs(outline[:, 1])))
outline = np.roll(outline, -tail, axis=0)
ring = resample_closed(outline, H)
# keep the start vertex at the tail extreme
ring = np.roll(ring, 0, axis=0)
NR = len(ring)
_t = [abs(math.atan2((ring[i] - ring[i - 1])[0] * (ring[(i + 1) % NR] - ring[i])[1] - (ring[i] - ring[i - 1])[1] * (ring[(i + 1) % NR] - ring[i])[0], (ring[i] - ring[i - 1]) @ (ring[(i + 1) % NR] - ring[i]))) for i in range(NR)]
_len = [float(np.hypot(*(ring[(i + 1) % NR] - ring[i]))) for i in range(NR)]
meta("ring_vertices", "%d (outline resampled with spacing min(%.1f mm, %.2f rad / curvature), at least %.1f mm; max turning angle per vertex %.1f deg, edge lengths %.1f to %.1f mm)" % (NR, H, MAX_TURN, MIN_RING_STEP, math.degrees(max(_t)), min(_len), max(_len)))

# sound hole circle (fit)
hc = hole_pts.mean(axis=0)
for _ in range(20):
    r = np.hypot(*(hole_pts - hc).T)
    # algebraic circle fit
    A = np.column_stack([2 * hole_pts[:, 0], 2 * hole_pts[:, 1], np.ones(len(hole_pts))])
    b = (hole_pts**2).sum(axis=1)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    hc = sol[:2]
    R_hole = math.sqrt(sol[2] + hc @ hc)
    break
meta("hole_centre_mm", "%.3f,%.3f (fit to 33 vector points)" % (hc[0], hc[1]))
meta("hole_radius_mm", "%.3f (fit; printed R43,5)" % R_hole)
nh = max(12, int(round(2 * math.pi * R_hole / (0.8 * H))))
hole_ring = np.array([[hc[0] + R_hole * math.cos(2 * math.pi * k / nh), hc[1] + R_hole * math.sin(2 * math.pi * k / nh)] for k in range(nh)])


def inside(poly, pt):
    x, y = pt
    c = False
    n = len(poly)
    j = n - 1
    for i in range(n):
        xi, yi = poly[i]
        xj, yj = poly[j]
        if (yi > y) != (yj > y) and x < (xj - xi) * (y - yi) / (yj - yi + 1e-300) + xi:
            c = not c
        j = i
    return c


def clip_segment(a, b, poly, hole_c=None, hole_r=None):
    """Return the longest sub-segment of a-b inside poly (and outside the hole)."""
    n = max(2, int(np.hypot(*(b - a)) / 0.25))
    ts = np.linspace(0, 1, n + 1)
    ok = []
    for t in ts:
        p = a + t * (b - a)
        good = inside(poly, p)
        if good and hole_c is not None and np.hypot(*(p - hole_c)) < hole_r:
            good = False
        ok.append(good)
    best = None
    start = None
    for i, g in enumerate(ok + [False]):
        if g and start is None:
            start = i
        if not g and start is not None:
            if best is None or (i - start) > (best[1] - best[0]):
                best = (start, i - 1)
            start = None
    if best is None:
        return None
    return a + ts[best[0]] * (b - a), a + ts[best[1]] * (b - a)


def snap_to_ring(p, ring_pts, tol):
    d = np.hypot(*(ring_pts - p).T)
    k = int(np.argmin(d))
    return (ring_pts[k].copy(), k) if d[k] < tol else (p, None)


# ---------------------------------------------------------------- brace centre lines (top)
hdr, strips = rows_with_header("top-braces-vector.tsv")
col = {c: i for i, c in enumerate(hdr)}
top_lines = []  # dict(name, a, b, width_mm, height_mm)
by_feature = {}
for r in strips:
    by_feature.setdefault(r[col["feature"]], []).append(r)
for feat, rs in by_feature.items():
    if len(rs) == 2:  # S06 + S07: one slanted bar drawn in two dashed runs; joined (the undrawn gap is bridged)
        pts = []
        for r in rs:
            pts += [np.array([float(r[col["x0_mm"]]), float(r[col["y0_mm"]])]), np.array([float(r[col["x1_mm"]]), float(r[col["y1_mm"]])])]
        far = max(((np.hypot(*(p - q)), i, j) for i, p in enumerate(pts) for j, q in enumerate(pts)), key=lambda z: z[0])
        a, b = pts[far[1]], pts[far[2]]
        meta("joined_strips_" + feat, "the two dashed runs are bridged into one straight bar (the gap between them is not drawn dashed; UNKNOWN whether the bar is continuous there)")
    else:
        r = rs[0]
        a = np.array([float(r[col["x0_mm"]]), float(r[col["y0_mm"]])])
        b = np.array([float(r[col["x1_mm"]]), float(r[col["y1_mm"]])])
    r0 = rs[0]
    top_lines.append(dict(name=feat, a=a, b=b, w=float(r0[col["width_mm"]]), h=float(r0[col["annotated_height_mm"]]) * BRACE_H_SCALE))

BRIDGE = opts.get("bridge", "0") == "1"
BRIDGE_MAT = opts.get("bridge_pack", "rosewood")
# Saddle station: the plan prints no bridge position; the centre of the gap between the two axis strips is INFERRED. It places the bridge
# wall when bridge=1 and the three drive points (outer strings and the G string) of the table in every case.
_neck_strip = [l for l in top_lines if l["name"] == "axis-strip-neck-side"][0]
_tail_strip = [l for l in top_lines if l["name"] == "axis-strip-tail-side"][0]
XB_SADDLE = float(opts.get("bridge_x", 0.5 * (max(_neck_strip["a"][0], _neck_strip["b"][0]) + min(_tail_strip["a"][0], _tail_strip["b"][0]))))
if BRIDGE:
    _bd = {(r[1], r[2]): float(r[3]) * 1000.0 for r in rows_with_header("dimensions.tsv")[1] if r[0] == "bridge"}
    b_len, b_w, b_h = _bd[("outline", "length")], _bd[("outline", "width")], _bd[("outline", "height")]
    b_mid, b_scoop = _bd[("central-section", "length")], _bd[("end-scoop", "depth")]
    xb = XB_SADDLE

    def bridge_h(y, b_len=b_len, b_mid=b_mid, b_h=b_h, b_scoop=b_scoop):
        ay = abs(y)
        if ay <= 0.5 * b_mid:
            return b_h
        f = min(1.0, (ay - 0.5 * b_mid) / max(1e-9, 0.5 * (b_len - b_mid)))
        return b_h - b_scoop * f

    top_lines.append(dict(name="bridge", a=np.array([xb, -0.5 * b_len]), b=np.array([xb, 0.5 * b_len]), w=b_w, h=b_h, side=1, hfun=bridge_h, pack=BRIDGE_MAT))
    meta("bridge_model", "bridge as a vertical wall on the OUTER skin of the top at x = %.2f mm (centre of the gap between the two axis strips, INFERRED; the plan prints no bridge position), %.1f mm long, thickness %.1f mm, height %.1f mm over the central %.1f mm falling by the printed %.1f mm scoop to the ends; pack %s (the plan names no bridge wood); the wall is thicker than it is tall (outside thin-shell validity, kept for its mass and its beam stiffness); the patch under the bridge (204 x 30 mm) and the string holes are not modelled" % (xb, b_len, b_w, b_h, b_mid, b_scoop, BRIDGE_MAT))

# ---------------------------------------------------------------- depth profile and back surface
# read the three printed depths by their parameter names
dep = {}
for r in rows_with_header("dimensions.tsv")[1]:
    if r[0] == "ribs" and r[1] == "depth":
        dep[r[2]] = float(r[3]) * 1000.0
x_lb = 179.0  # station of the printed 101 mm read from the side view (dimensions.tsv note); mm from the upper edge
xmax = float(outline[:, 0].max())
xmin = float(outline[:, 0].min())
pts_d = np.array([[xmin, dep["at-upper-edge"]], [x_lb, dep["at-lower-bout-start"]], [xmax, dep["at-tail-end"]]])
coef = np.polyfit(pts_d[:, 0], pts_d[:, 1], 2)
meta("rib_depth_profile", "depth below the top plane line through the printed 94 (x=%.1f), 101 (x=179) and 99 mm (x=%.1f): quadratic %.6g x^2 + %.6g x + %.6g mm; the profile between the three printed values is authored" % (xmin, xmax, coef[0], coef[1], coef[2]))


def depth(x):
    if BACK_ARCH == "flat":
        return float(np.mean(pts_d[:, 1])) + 0 * x
    return np.polyval(coef, x)


def ddepth(x):
    if BACK_ARCH == "flat":
        return 0.0 * x
    return np.polyval(np.polyder(coef), x)


# top: midsurface z = 0, outer face +z/2. back outer face at depth d below the top plane (top outer face), back midsurface:
def z_back_mid(x):
    return t_top / 2.0 - depth(x) + t_back / 2.0


def z_back_inner(x):  # inner skin of the back (towards the cavity)
    return z_back_mid(x) + t_back / 2.0


z_top_inner = -t_top / 2.0
meta("vertical_frame", "z=0 top midsurface, +z toward the player (outer face of the top at +t/2); the back outer face lies at depth d(x) below the top outer face; ribs run from the top inner skin to the back inner skin")


# ---------------------------------------------------------------- geometry containers (output, SI metres)
class Comp:
    def __init__(self, name, pack, t_mm, exterior, cavity):
        self.name, self.pack, self.t = name, pack, t_mm
        self.exterior, self.cavity = exterior, cavity
        self.nodes = []  # x,y,z mm, dx,dy,dz
        self.tris = []  # i,j,k,region
        self.regions = []  # axis vectors

    def add_node(self, p, d):
        self.nodes.append((p[0], p[1], p[2], d[0], d[1], d[2]))
        return len(self.nodes) - 1


components = []
ties = []  # (master, slave, side, [m nodes], [s nodes])


def orient_tris(comp, want_dir_from_nodes=True):
    """Flip each triangle winding so the facet normal agrees with the mean nodal director."""
    P = np.array([n[:3] for n in comp.nodes])
    Dd = np.array([n[3:] for n in comp.nodes])
    out = []
    for (i, j, k, g) in comp.tris:
        n = np.cross(P[j] - P[i], P[k] - P[i])
        dm = Dd[i] + Dd[j] + Dd[k]
        if n @ dm < 0:
            j, k = k, j
        out.append((i, j, k, g))
    comp.tris = out


def tri_min_angle(P, T):
    mins = []
    for (i, j, k) in T:
        a, b, c = P[i], P[j], P[k]
        angs = []
        for u, v, w in ((a, b, c), (b, c, a), (c, a, b)):
            e1, e2 = v - u, w - u
            cs = e1 @ e2 / (np.linalg.norm(e1) * np.linalg.norm(e2) + 1e-300)
            angs.append(math.degrees(math.acos(max(-1, min(1, cs)))))
        mins.append(min(angs))
    return min(mins), float(np.mean(mins))


def plate_mesh(ring_pts, hole_pts_, lines, h, label):
    """Constrained Delaunay mesh (Shewchuk Triangle, no Steiner points on the boundary or on the given lines).
    lines: list of dict(a,b,...) 2-D centre lines in mm. Returns vertices (mm), triangles, per-line ordered vertex index lists."""
    verts = list(map(tuple, ring_pts))
    segs = [(i, (i + 1) % len(ring_pts)) for i in range(len(ring_pts))]
    off = len(verts)
    if hole_pts_ is not None:
        verts += list(map(tuple, hole_pts_))
        segs += [(off + i, off + (i + 1) % len(hole_pts_)) for i in range(len(hole_pts_))]
    key = {}

    def vid(p):
        k = (round(p[0], 2), round(p[1], 2))
        if k in key:
            return key[k]
        for i, q in enumerate(verts):  # existing vertex within 1e-2 mm
            if abs(q[0] - p[0]) < 0.01 and abs(q[1] - p[1]) < 0.01:
                key[k] = i
                return i
        verts.append((float(p[0]), float(p[1])))
        key[k] = len(verts) - 1
        return key[k]

    # clip each line, snap its ends to ring vertices, cut at mutual intersections
    prepared = []
    for ln in lines:
        c = clip_segment(ln["a"], ln["b"], ring_pts, None if hole_pts_ is None else hc, R_hole + 0.4 * h)
        if c is None:
            continue
        a, b = c
        a, ka = snap_to_ring(a, ring_pts, 0.9 * h)
        b, kb = snap_to_ring(b, ring_pts, 0.9 * h)
        prepared.append(dict(ln, a=a, b=b))
    # intersections
    cuts = [[] for _ in prepared]
    for i in range(len(prepared)):
        for j in range(i + 1, len(prepared)):
            A, B = prepared[i]["a"], prepared[i]["b"]
            C, Dp = prepared[j]["a"], prepared[j]["b"]
            r_, s_ = B - A, Dp - C
            den = r_[0] * s_[1] - r_[1] * s_[0]
            if abs(den) < 1e-9:
                continue
            t = ((C[0] - A[0]) * s_[1] - (C[1] - A[1]) * s_[0]) / den
            u = ((C[0] - A[0]) * r_[1] - (C[1] - A[1]) * r_[0]) / den
            if 0.02 < t < 0.98 and 0.02 < u < 0.98:
                P = A + t * r_
                cuts[i].append(t)
                cuts[j].append(u)
    line_vertices = []
    for ln, ct in zip(prepared, cuts):
        a, b = ln["a"], ln["b"]
        ts = sorted(set([0.0, 1.0] + ct))
        chain = []
        for t0, t1 in zip(ts[:-1], ts[1:]):
            p0, p1 = a + t0 * (b - a), a + t1 * (b - a)
            n = max(1, int(round(np.hypot(*(p1 - p0)) / h)))
            sub = [vid(p0 + (p1 - p0) * k / n) for k in range(n + 1)]
            for u_, v_ in zip(sub[:-1], sub[1:]):
                if u_ != v_:
                    segs.append((u_, v_))
            if chain and chain[-1] == sub[0]:
                sub = sub[1:]
            chain += sub
        line_vertices.append(chain)
    # merge duplicate segments
    segs = sorted(set(tuple(sorted(s)) for s in segs))
    pslg = dict(vertices=np.array(verts), segments=np.array(segs))
    if hole_pts_ is not None:
        pslg["holes"] = np.array([hc])
    res = triangle.triangulate(pslg, "pq25YYa%.4f" % (0.4330127 * h * h * 1.15))
    V = res["vertices"]
    T = res["triangles"]
    mn, mean = tri_min_angle(V, T)
    meta("mesh_%s" % label, "%d vertices, %d triangles, min angle %.1f deg, mean of per-triangle minimum angles %.1f deg" % (len(V), len(T), mn, mean))
    assert len(V) >= len(verts)
    # the first len(verts) vertices keep their input order (Triangle without Steiner points on segments)
    for i, v in enumerate(verts):
        assert abs(V[i][0] - v[0]) < 1e-9 and abs(V[i][1] - v[1]) < 1e-9, (label, i)
    # per-line ordered vertex lists
    lv = []
    for ln, chain in zip(prepared, line_vertices):
        lv.append((ln, chain))
    return V, T, lv


# ---------------------------------------------------------------- top plate
top = Comp("top", TOP_MAT, t_top, True, True)
V, T, lv = plate_mesh(ring, hole_ring, top_lines, H, "top")
zt = lambda p: 0.0
for p in V:
    top.add_node((p[0], p[1], 0.0), (0.0, 0.0, 1.0))
top.regions = [(1.0, 0.0, 0.0)]
for (i, j, k) in T:
    top.tris.append((int(i), int(j), int(k), 0))
orient_tris(top)
top_ring_idx = list(range(NR))  # first NR vertices are the ring in order
top_lines_mesh = lv
components.append(top)

# brace walls (and the bridge, which stands on the outer skin)
for ln, chain in top_lines_mesh:
    pts = np.array([V[c] for c in chain])
    d = ln["b"] - ln["a"]
    L = np.hypot(*d)
    u = d / L
    nrm = np.array([-u[1], u[0], 0.0])
    side = ln.get("side", -1)
    sgn = 1.0 if side == 1 else -1.0
    z0 = sgn * t_top / 2.0
    h_mm = ln["h"]
    w_mm = ln["w"]
    hfun = ln.get("hfun")
    spacing = np.hypot(*(pts[1:] - pts[:-1]).T).mean()
    K = max(WALL_K_MIN if h_mm >= 6.0 else 1, int(math.ceil(h_mm / (0.8 * spacing))))
    # walls are not in group "top" (a component's group is its name up to the first '-'): the plate means of the shared body module use the plates only
    wall = Comp(("brace-" + ln["name"]) if side == -1 else ln["name"], ln.get("pack", BRACE_MAT), w_mm, False, False)
    wall.regions = [(u[0], u[1], 0.0)]
    idx = []
    for k in range(K + 1):
        row = []
        for p in pts:
            h_i = hfun(p[1]) if hfun else h_mm
            row.append(wall.add_node((p[0], p[1], z0 + sgn * h_i * k / K), tuple(nrm)))
        idx.append(row)
    for k in range(K):
        for i in range(len(pts) - 1):
            a, b, c, e = idx[k][i], idx[k][i + 1], idx[k + 1][i + 1], idx[k + 1][i]
            if (i + k) % 2 == 0:
                wall.tris += [(a, b, c, 0), (a, c, e, 0)]
            else:
                wall.tris += [(a, b, e, 0), (b, c, e, 0)]
    orient_tris(wall)
    components.append(wall)
    ties.append(("top", wall.name, side, [int(c) for c in chain], idx[0][:]))
meta("brace_model", "each top brace is a vertical shell wall (thickness = printed width, height = printed/assigned h:, grain along the bar, rectangular cross-section, full length, no taper or scallop) tied rigidly to the lower skin of the top along its centre line; bar crossings are not notched; brace cross-section shape and ends are UNKNOWN")

# ---------------------------------------------------------------- back plate with ladder bars
bars = []
bar_x = None
dd = {}
for r in rows_with_header("dimensions.tsv")[1]:
    dd[(r[0], r[1], r[2])] = float(r[3]) * 1000.0
x1 = dd[("back", "first-brace", "distance-from-neck-end")] + xmin
sp12 = 0.5 * (dd[("back-braces", "brace-1-to-2", "spacing-treble-side")] + dd[("back-braces", "brace-1-to-2", "spacing-bass-side")])
sp23 = 0.5 * (dd[("back-braces", "brace-2-to-3", "spacing-treble-side")] + dd[("back-braces", "brace-2-to-3", "spacing-bass-side")])
sp34 = 0.5 * (dd[("back-braces", "brace-3-to-4", "spacing-treble-side")] + dd[("back-braces", "brace-3-to-4", "spacing-bass-side")])
bar_x = [x1, x1 + sp12, x1 + sp12 + sp23, x1 + sp12 + sp23 + sp34]
meta("back_bar_positions_mm", "%s (bar 1 at the printed 91 mm from the neck end, then the mean of the printed treble-side and bass-side spacings 106/107.5, 88/85.5, 85/83; the bars are drawn slightly slanted, modelled straight and transverse; the back view is a mirror of the front view whose handedness is UNKNOWN)" % ",".join("%.2f" % x for x in bar_x))
back_lines = []
for k, x in enumerate(bar_x):
    back_lines.append(dict(name="bar-%d" % (k + 1), a=np.array([x, -250.0]), b=np.array([x, 250.0]),
                           w=dd[("back-braces", "brace-%d" % (k + 1), "width")], h=dd[("back-braces", "brace-%d" % (k + 1), "annotated-height")] * BACK_BAR_H_SCALE))
back = Comp("back", BACK_MAT, t_back, True, True)
Vb, Tb, lvb = plate_mesh(ring, None, back_lines, H, "back")
for p in Vb:
    zm = float(z_back_mid(p[0]))
    # normal of the surface z = zm(x) pointing outward from the cavity (toward -z): (-dz/dx, 0, -1)/|.|
    dz = -float(ddepth(p[0]))  # d(z_back_mid)/dx = -d depth/dx
    n = np.array([dz, 0.0, -1.0])
    n /= np.linalg.norm(n)
    back.add_node((p[0], p[1], zm), tuple(n))
back.regions = [(1.0, 0.0, 0.0)]
for (i, j, k) in Tb:
    back.tris.append((int(i), int(j), int(k), 0))
orient_tris(back)
components.append(back)
for ln, chain in lvb:
    pts = np.array([Vb[c] for c in chain])
    d = ln["b"] - ln["a"]
    u = d / np.hypot(*d)
    h_mm, w_mm = ln["h"], ln["w"]
    spacing = np.hypot(*(pts[1:] - pts[:-1]).T).mean()
    K = max(WALL_K_MIN if h_mm >= 6.0 else 1, int(math.ceil(h_mm / (0.8 * spacing))))
    wall = Comp(ln["name"], BACKBAR_MAT, w_mm, False, False)
    wall.regions = [(u[0], u[1], 0.0)]
    idx = []
    # The bar stands on the inner skin along the back plate's NODAL DIRECTOR (the tie shifts the attachment by -t/2 times that
    # director), so its bottom nodes sit exactly on the tied skin curve and the wall rises along -director.
    for k in range(K + 1):
        row = []
        for c in chain:
            dn = np.array(back.nodes[c][3:6])
            dn = dn / np.linalg.norm(dn)
            base = np.array(back.nodes[c][:3]) - 0.5 * t_back * dn
            e_up = -dn
            wall_dir = np.cross(u.tolist() + [0.0], e_up)
            wall_dir = wall_dir / np.linalg.norm(wall_dir)
            row.append(wall.add_node(tuple(base + e_up * h_mm * k / K), tuple(wall_dir)))
        idx.append(row)
    for k in range(K):
        for i in range(len(pts) - 1):
            a, b, c, e = idx[k][i], idx[k][i + 1], idx[k + 1][i + 1], idx[k + 1][i]
            if (i + k) % 2 == 0:
                wall.tris += [(a, b, c, 0), (a, c, e, 0)]
            else:
                wall.tris += [(a, b, e, 0), (b, c, e, 0)]
    orient_tris(wall)
    components.append(wall)
    ties.append(("back", wall.name, -1, [int(c) for c in chain], idx[0][:]))
meta("back_model", "back planform = the top outline (the back view outline is 495 mm long against 489.3 mm, UNKNOWN reason); back surface = cylinder over the quadratic rib-depth profile (no cross arch; the printed arch 4 mm is not modelled); four transverse ladder bars as vertical walls of printed width and height 17 mm (shape, taper and arch cut of the bars UNKNOWN); the central strip (h 1.5 mm) and the 6 x 15 mm patches are not modelled")

# ---------------------------------------------------------------- ribs: two strips (treble side, bass side) meeting at the tail and the neck end
tang = []
for i in range(NR):
    t = ring[(i + 1) % NR] - ring[i - 1]
    tang.append(t / np.hypot(*t))
tang = np.array(tang)
m_neck = int(np.argmin(ring[:, 0]))
Kmax = max(2, int(math.ceil(max(z_top_inner - float(z_back_inner(x)) for x in ring[:, 0]) / (0.9 * H))))
SECT_PER = 12


def build_rib(name, c0, c1):
    comp = Comp(name, RIB_MAT, t_rib, True, True)
    cs = list(range(c0, c1 + 1))
    nsec = SECT_PER
    sect = lambda q: min(nsec - 1, (q * nsec) // (len(cs) - 1))
    axes = []
    for s_ in range(nsec):
        ids = [cs[q] % NR for q in range(len(cs)) if sect(q) == s_] or [cs[0] % NR]
        m = tang[ids].mean(axis=0)
        m /= np.hypot(*m)
        axes.append((m[0], m[1], 0.0))
    comp.regions = axes
    cols_ = []
    for q, c in enumerate(cs):
        i = c % NR
        x, y = ring[i]
        n = np.array([tang[i][1], -tang[i][0], 0.0])
        # top end on the lower skin of the top (z = -t/2 under the ring vertex); bottom end exactly on the back's inner skin, which
        # the tie places at the back vertex minus t/2 times its nodal director
        top_end = np.array([x, y, z_top_inner])
        dn = np.array(back.nodes[i][3:6])
        dn = dn / np.linalg.norm(dn)
        bottom_end = np.array(back.nodes[i][:3]) - 0.5 * t_back * dn
        cols_.append([comp.add_node(tuple(top_end + (bottom_end - top_end) * k / Kmax), tuple(n)) for k in range(Kmax + 1)])
    for q in range(len(cs) - 1):
        g = sect(q)
        for k in range(Kmax):
            a, b, cc, e = cols_[q][k], cols_[q + 1][k], cols_[q + 1][k + 1], cols_[q][k + 1]
            if (q + k) % 2 == 0:
                comp.tris += [(a, b, cc, g), (a, cc, e, g)]
            else:
                comp.tris += [(a, b, e, g), (b, cc, e, g)]
    orient_tris(comp)
    return comp, cols_, cs


rib_t, cols_t, cs_t = build_rib("rib-treble", 0, m_neck)
rib_b, cols_b, cs_b = build_rib("rib-bass", m_neck, NR)
components += [rib_t, rib_b]
meta("rib_model", "two strips (treble side from the tail to the neck end, bass side back to the tail), %d vertical layers, uniform thickness, grain along the local tangent in %d sectors each; the strips meet at the tail and at the neck end through vertical rigid ties; neck block, tail block and linings are not modelled" % (Kmax, SECT_PER))
# rib <-> top and rib <-> back along the outline (the master polyline is the plate boundary, repeating vertex 0 to close the ring)
ties.append(("top", "rib-treble", -1, [c % NR for c in cs_t], [col[0] for col in cols_t]))
ties.append(("top", "rib-bass", -1, [c % NR for c in cs_b], [col[0] for col in cols_b]))
ties.append(("back", "rib-treble", -1, [c % NR for c in cs_t], [col[Kmax] for col in cols_t]))
ties.append(("back", "rib-bass", -1, [c % NR for c in cs_b], [col[Kmax] for col in cols_b]))
# vertical seams between the strips (interior nodes only: the end nodes already follow the plates)
ties.append(("rib-treble", "rib-bass", 0, cols_t[-1][1:Kmax], cols_b[0][1:Kmax]))
ties.append(("rib-treble", "rib-bass", 0, cols_t[0][1:Kmax], cols_b[-1][1:Kmax]))

# ---------------------------------------------------------------- cavity volume (divergence theorem on the midsurface facets)
vol_mid = 0.0
area_by = {}
for comp in components:
    if not comp.cavity:
        continue
    P = np.array([n[:3] for n in comp.nodes])
    a_tot = 0.0
    for (i, j, k, g) in comp.tris:
        n = np.cross(P[j] - P[i], P[k] - P[i])
        A = 0.5 * np.linalg.norm(n)
        c = (P[i] + P[j] + P[k]) / 3.0
        vol_mid += (c @ (n / 2.0)) / 3.0  # (1/3) x . n dA with n dA = 0.5 * cross
        a_tot += A
    area_by[comp.name] = a_tot
# close the sound hole with a flat disc at z = 0 (outward normal +z): contributes (1/3) * 0 = 0 since z = 0 on the disc
inner_corr = sum(area_by[c.name] * c.t / 2.0 for c in components if c.cavity)  # mm^3
vol_inner_mm3 = vol_mid - inner_corr
meta("cavity_volume_midsurface_L", "%.4f" % (vol_mid / 1e6))
meta("cavity_volume_L", "%.4f (midsurface facets with outward-from-cavity orientation, hole closed by a disc in the top midsurface plane, minus half-thickness on each wall; the plan prints no volume)" % (vol_inner_mm3 / 1e6))
meta("cavity_areas_mm2", ",".join("%s=%.0f" % kv for kv in area_by.items()))
meta("neck", "radius %.4f mm (fit to the vector circle), length = top thickness %.4f mm, both ends flanged (authored choice)" % (R_hole, t_top))

# ---------------------------------------------------------------- write
total_nodes = sum(len(c.nodes) for c in components)
total_tris = sum(len(c.tris) for c in components)
meta("totals", "%d components, %d nodes, %d triangles" % (len(components), total_nodes, total_tris))
OMITTED = "bridge%s, neck and heel, tail and heel blocks, linings, rosette, fingerboard, strings, the central back strip and the back patches, glue compliance, finish" % (" (modelled: wall on the outer skin)" if BRIDGE else "")
if FORMAT == "shared":
    meta("omitted", OMITTED)
with open(OUT, "w", encoding="utf-8") as f:
    f.write(("frankensim.shell-assembly.v1" if FORMAT == "shared" else "frankensim.guitar-shell-assembly.v1") + chr(10))
    for k, v in META:
        f.write("meta\t%s\t%s\n" % (k, v))
    for c in components:
        f.write("component\t%s\t%s\t%.9f\t%d\t%d\t%d\t%d\t%d\n" % (c.name, c.pack, c.t / 1000.0, len(c.regions), len(c.nodes), len(c.tris), int(c.exterior), int(c.cavity)))
        for r in c.regions:
            f.write("region\t%.9f\t%.9f\t%.9f\n" % tuple(r))
        for n in c.nodes:
            d = np.array(n[3:])
            d = d / np.linalg.norm(d)
            f.write("node\t%.9f\t%.9f\t%.9f\t%.12f\t%.12f\t%.12f\n" % (n[0] / 1000.0, n[1] / 1000.0, n[2] / 1000.0, d[0], d[1], d[2]))
        for t in c.tris:
            f.write("tri\t%d\t%d\t%d\t%d\n" % t)
    for (m, s, side, mn, sn) in ties:
        f.write("tie\t%s\t%s\t%d\t%d\t%d\n" % (m, s, side, len(mn), len(sn)))
        f.write("master\t%s\n" % " ".join(map(str, mn)))
        f.write("slave\t%s\n" % " ".join(map(str, sn)))
    if FORMAT == "shared":
        # three saddle drive points on the top: the two outer string positions and the G string (string-hole span 58.5 mm, holes 11.7 mm apart;
        # which outer string is on which side is UNKNOWN), unit force along +z, at the inferred saddle station
        for label, y_mm in (("outer string, -y side", -29.25), ("G string position", 5.85), ("outer string, +y side", 29.25)):
            f.write("drive\t%s\ttop\t%.6f\t%.6f\t0.000000\t0\t0\t1\n" % (label, XB_SADDLE / 1000.0, y_mm / 1000.0))
    f.write("cavity\t%.9f\t%.9f\t%.9f\n" % (vol_inner_mm3 / 1e9, R_hole / 1000.0, t_top / 1000.0))
for k, v in META:
    print(k, "=", v)
