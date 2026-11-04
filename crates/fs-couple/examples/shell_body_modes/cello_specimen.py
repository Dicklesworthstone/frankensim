"""A CELLO specimen directory for gen_violin_shell_assembly.py (instrument=cello), made by mapping the il Cannone
specimen tables (ESRF CT, the h6 violin's inputs) to cello dimensions from published sources. Every mapped quantity is
written with its source or labelled ASSUMPTION in PROVENANCE.tsv beside the tables. Nothing here is tuned to any
recording or mode frequency. cello_bridge_19g.py then sets the bridge to its published mass; the two together write
data/instruments/geometry/cello-castelbarco-1697-shell (check_cello_shell.py reproduces it byte for byte).

usage: python3 -I cello_specimen.py CANNONE_DIR OUT_DIR   (OUT_DIR must not exist)

Sources (values in mm unless stated):
  LOC   = Library of Congress catalog record 2022560102, Stradivari cello "Castelbarco" 1697, via the repo table
          data/instruments/geometry/loc-stradivari-family-bouts/dimensions.tsv: length of back 749, bouts 346 / 234 / 440.
  JOH   = Hans Johannsson, cello specification table, hansjohannsson.com/resources/mescello.html, 1/1 column
          (copy md5 921b8776ddb22adc55c01e9e6f9e44cc): rib height average 120, rib thickness 1.7, arching back 28 / front 29,
          back thickness: centre at sound post 7, flanks 3.5, fluting 3.5, centre of bouts 3.6, edge 5;
          front thickness: centre at sound post 4.7, flanks 3.9, fluting 3.9, centre of bouts 4, edge 5.5;
          purfling distance from edge 5.5-6; lining width 16; soundhole length 130, between top holes 98, bottom hole to edge 28;
          bass bar length approx. 580, thickness approx. 11, height middle approx. 24; soundpost thickness 11;
          bridge feet width 90, thickness 11, height approx. 92; string spacing at bridge 47; stoplength 405.
  LOCM  = the measured f-hole drawing of the same "Castelbarco" cello in the Library of Congress record above
          (item ihas.200154818, image 0021: "Cello the Castelbarco Antonius Stradivarius 1697", John Montgomery,
          December 2005, "all measurements in millimeters with calipers except where noted"; copy md5 c7f17bbeb193):
          f-hole lengths 130.5 (bass) / 128.5 (treble); upper eyes 102.3 apart (inner edges; 103.0 by tape), diameters
          13.5 / 13.0; lower eyes 220.0 apart (inner edges) and 258.2 (outer edges), diameters 19.7 / 19.2; lower eye to
          the outside of the instrument 32.8 / 34; stop to the notches 419 / 420; notches 162.0 apart.
  CAN   = the il Cannone tables themselves (shape only: outline, arch shape, the back's thickness pattern, hole shapes,
          bar profile, bridge silhouette, positions), mapped as stated.
The f-holes are placed from LOCM, and the outputs that depend on them (such as the top's island thickness zone) follow
them; PROVENANCE.tsv has one row per printed f-hole value (f_hole_*).
"""
import math
import os
import sys

import numpy as np

SRC = sys.argv[1].rstrip("/") + "/"
OUT = sys.argv[2].rstrip("/") + "/"
if os.path.exists(OUT):
    raise SystemExit("refusing: %s exists" % OUT)
os.makedirs(OUT)
PROV = []


def prov(item, value, source, note):
    PROV.append((item, value, source, note))


def nt(name):
    lines = open(SRC + name, encoding="utf-8").read().splitlines()
    return lines[0], lines[1].split("\t"), np.array([[float(x) for x in l.split("\t")] for l in lines[2:] if l.strip()])


def write(name, magic, head, rows, fmt):
    with open(OUT + name, "x", encoding="utf-8") as f:
        f.write(magic + "\n" + "\t".join(head) + "\n")
        for r in rows:
            f.write("\t".join(fm % v if not isinstance(v, str) else v for fm, v in zip(fmt, r)) + "\n")


# ------------------------------------------------------------------ targets
L_BACK, B_UP, B_C, B_LO = 749.0, 346.0, 234.0, 440.0  # LOC
RIB_H, RIB_T = 120.0, 1.7  # JOH
ARCH_TOP, ARCH_BACK = 29.0, 28.0  # JOH
T_TOP = dict(post=4.7, flank=3.9, flute=3.9, bout=4.0, edge=5.5)  # JOH front
T_BACK = dict(post=7.0, flank=3.5, flute=3.5, bout=3.6, edge=5.0)  # JOH back
PURF, LINING_W = 6.0, 16.0  # JOH purfling distance 5.5-6 (upper end), lining width 16
# (the f-holes come from LOCM, the measured drawing of this cello, not from JOH 130 / 98 / 28; see the f-hole block)
BAR_L, BAR_W, BAR_H = 580.0, 11.0, 24.0  # JOH
POST_D = 11.0  # JOH
BR_FEET, BR_T_FOOT, BR_H = 90.0, 11.0, 92.0  # JOH
STR_SPAN = 47.0  # JOH

# ------------------------------------------------------------------ the Cannone references
_, _, top_ol = nt("top-plate-outline.tsv")
_, _, back_ol = nt("back-plate-outline.tsv")
wt = top_ol[:, 2] - top_ol[:, 1]
L_CAN = float(top_ol[wt >= 10, 0].max() - top_ol[wt >= 10, 0].min() + 1.0)  # 355.0: the body without the back's button
wb = back_ol[:, 2] - back_ol[:, 1]
zb = back_ol[:, 0]
kL = int(np.argmax(np.where(zb < -87, wb, -1)))
kU = int(np.argmax(np.where(zb >= -87, wb, -1)))
sel = (zb > zb[kL]) & (zb < zb[kU])
kC = int(np.argmin(np.where(sel, wb, 1e9)))
SZ = L_BACK / L_CAN
ST_Z = np.array([zb[kL], zb[kC], zb[kU]])
ST_S = np.array([B_LO / wb[kL], B_C / wb[kC], B_UP / wb[kU]])
prov("plan_length_scale", "%.6f" % SZ, "LOC / CAN", "z' = %.4f z: LOC length of back 749 mm over the Cannone body length %.1f mm (top outline rows at least 10 mm wide, i.e. without the back's button)" % (SZ, L_CAN))
prov("plan_width_scale", "lower %.4f at z %.2f, C %.4f at z %.2f, upper %.4f at z %.2f (Cannone z)" % (ST_S[0], ST_Z[0], ST_S[1], ST_Z[1], ST_S[2], ST_Z[2]), "LOC / CAN",
     "x' = s(z) x, s piecewise linear in the Cannone z through the three bout stations of the Cannone back (widths %.1f, %.1f, %.1f mm -> LOC 440, 234, 346), constant beyond; ASSUMPTION: the cello outline is the Cannone outline so stretched (no cello outline curve is published in reach)" % (wb[kL], wb[kC], wb[kU]))


def sx(zc):  # Cannone z -> width scale
    return np.interp(zc, ST_Z, ST_S)


def fwd(x, z):  # Cannone plan -> cello plan
    x, z = np.asarray(x, float), np.asarray(z, float)
    return sx(z) * x, SZ * z


def inv(xp, zp):  # cello plan -> Cannone plan
    z = np.asarray(zp, float) / SZ
    return np.asarray(xp, float) / sx(z), z


# ------------------------------------------------------------------ Cannone plate fields (the generator's own processing: 3 mm blocks, fill, smooth 2 / 3)
class Field:
    def __init__(self, pre):
        _, _, m = nt(pre + "-plate-maps.tsv")
        self.x0, self.z0 = m[:, 0].min(), m[:, 1].min()
        self.nx = int(round((m[:, 0].max() - self.x0) / 3)) + 1
        self.nz = int(round((m[:, 1].max() - self.z0) / 3)) + 1
        ix, iz = np.rint((m[:, 0] - self.x0) / 3).astype(int), np.rint((m[:, 1] - self.z0) / 3).astype(int)
        g = []
        for col in (2, 3, 4):
            a = np.full((self.nx, self.nz), np.nan)
            a[ix, iz] = m[:, col]
            g.append(self.fill(a))
        self.ymid = self.smooth(0.5 * (g[0] + g[1]), 2)
        self.thick = self.smooth(g[2], 3)

    @staticmethod
    def fill(g):
        g = g.copy()
        for _ in range(200):
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

    def bil(self, g, x, z):
        u = np.clip((np.asarray(x) - self.x0) / 3, 0, self.nx - 1.001)
        v = np.clip((np.asarray(z) - self.z0) / 3, 0, self.nz - 1.001)
        i, j = np.floor(u).astype(int), np.floor(v).astype(int)
        fu, fv = u - i, v - j
        return g[i, j] * (1 - fu) * (1 - fv) + g[i + 1, j] * fu * (1 - fv) + g[i, j + 1] * (1 - fu) * fv + g[i + 1, j + 1] * fu * fv


FT, FB = Field("top"), Field("back")
_, _, rs = nt("ribs-sections.tsv")
th = np.radians(rs[:, 0])
rmid = 0.5 * (rs[:, 4] + rs[:, 5])


def ring_points(cx, cz, r, sgn):
    return np.c_[cx + r * np.cos(th), cz + sgn * r * np.sin(th)]


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
        out ^= m & (xs < (x2 - x1) * (zs - z1) / (z2 - z1) + x1)
    return out


def outline_poly(a):
    w = a[:, 2] - a[:, 1]
    k = w >= 10.0
    return np.vstack([np.c_[a[k, 2], a[k, 0]], np.c_[a[k, 1][::-1], a[k, 0][::-1]]])


TOPPOLY = outline_poly(top_ol)
sgn = max((1.0, -1.0), key=lambda s: inside_many(TOPPOLY, ring_points(rs[:, 1], rs[:, 2], rmid, s)).sum())
RING_C = ring_points(rs[:, 1], rs[:, 2], rmid, sgn)  # Cannone ring (x, z)
# Cannone glue planes: the plates' inner skins at the ring (where the generator's rib strip ends)
inner_t_c = FT.bil(FT.ymid, RING_C[:, 0], RING_C[:, 1]) + 0.5 * FT.bil(FT.thick, RING_C[:, 0], RING_C[:, 1])
inner_b_c = FB.bil(FB.ymid, RING_C[:, 0], RING_C[:, 1]) - 0.5 * FB.bil(FB.thick, RING_C[:, 0], RING_C[:, 1])
E_T, E_B = float(FT.bil(FT.ymid, RING_C[:, 0], RING_C[:, 1]).mean()), float(FB.bil(FB.ymid, RING_C[:, 0], RING_C[:, 1]).mean())

# ------------------------------------------------------------------ cello plan outlines
def map_outline(a):
    w = a[:, 2] - a[:, 1]
    k = w >= 10.0  # the generator leaves out rows narrower than 10 mm (button and notch): decided in Cannone units
    xl, z2 = fwd(a[k, 1], a[k, 0])
    xr, _ = fwd(a[k, 2], a[k, 0])
    return np.c_[z2, xl, xr]


top_rows, back_rows = map_outline(top_ol), map_outline(back_ol)
for name, rows in (("top", top_rows), ("back", back_rows)):
    write(name + "-plate-outline.tsv", "frankensim.specimen-outline-rows.v1", ["z_mm", "x_left_mm", "x_right_mm"], rows, ["%.3f"] * 3)
cpoly_t = np.vstack([np.c_[top_rows[:, 2], top_rows[:, 0]], np.c_[top_rows[::-1, 1], top_rows[::-1, 0]]])
cpoly_b = np.vstack([np.c_[back_rows[:, 2], back_rows[:, 0]], np.c_[back_rows[::-1, 1], back_rows[::-1, 0]]])
wcb = back_rows[:, 2] - back_rows[:, 1]
prov("outline_check", "back length %.1f, bouts lower %.1f / C %.1f / upper %.1f mm" % (back_rows[:, 0].max() - back_rows[:, 0].min() + SZ, wcb[back_rows[:, 0] < -183].max(),
     wcb[(back_rows[:, 0] > ST_Z[0] * SZ) & (back_rows[:, 0] < ST_Z[2] * SZ)].min(), wcb[back_rows[:, 0] > -183].max()), "check", "the mapped back outline against LOC 749 / 440 / 234 / 346 (length includes the button rows that are at least 10 mm wide in Cannone units)")


def dist_to_poly(poly, pts):
    P = np.vstack([poly, poly[:1]])
    d = np.full(len(pts), np.inf)
    for a, b in zip(P[:-1], P[1:]):
        ab = b - a
        L2 = float(ab @ ab)
        if L2 == 0:
            continue
        t = np.clip(((pts - a) @ ab) / L2, 0, 1)
        q = a + t[:, None] * ab
        d = np.minimum(d, np.hypot(*(pts - q).T))
    return d


# ------------------------------------------------------------------ f-holes (needed for the top's island before the maps)
# Position and size from the printed measurements of THIS cello (LOCM, see the docstring); only the hole SHAPE stays the
# Cannone's. Placing the top eyes at JOH 98 mm and rotating each hole until its bottom-eye centroid is JOH 28 mm from the
# outline instead puts the bass hole's lower wing over the rib glue ring.
FH_LEN_P = {1: 130.5, 2: 128.5}  # LOCM f-hole length end to end (the drawn line from the top end to the bottom end), bass / treble
UP_INNER = 102.3  # LOCM between the upper eyes (inner edges), caliper (103.0 by tape)
UP_DIA = {1: 13.5, 2: 13.0}  # LOCM upper eye diameters
LO_INNER, LO_OUTER = 220.0, 258.2  # LOCM between the lower eyes' inner edges, and between their outer edges
LO_DIA = {1: 19.7, 2: 19.2}  # LOCM lower eye diameters (check only)
STOP_P = {1: 419.0, 2: 420.0}  # LOCM stop: top edge to the f-hole notches, bass / treble
NOTCH_GAP_P = 162.0  # LOCM between the notches (check only)
LO_EDGE_P = {1: 32.8, 2: 34.0}  # LOCM lower eye to the outside of the instrument, bass / treble (check only)
_, _, fh = nt("top-plate-fholes.tsv")
_, _, bsec_c = nt("bridge-sections.tsv")
Z_NOTCH_C = float(np.average(bsec_c[:, 3], weights=bsec_c[:, 1]))  # ASSUMPTION: the Cannone's notch line is its bridge plane
Z_TOP = SZ * (float(top_ol[:, 0].max()) + 0.5)  # ASSUMPTION: the top's edge on the centre line is half a 1 mm row beyond its last row


def eye_circle(r, zlo, zhi):
    """Largest inscribed circle of raster hole r (rows at k 0.5 + 0.25, runs on the 0.5 mm grid) centred in zlo < z < zhi."""
    cells = set()
    for _, z, xa, xb in r:
        for i in range(int(round((xb - xa) / 0.5))):
            cells.add((round((xa + 0.25 + 0.5 * i) * 4), round(z * 4)))
    out = {(cx + dx, cz + dz) for cx, cz in cells for dx, dz in ((2, 0), (-2, 0), (0, 2), (0, -2))} - cells
    O = np.array(sorted(out), float) / 4.0
    P = np.array(sorted(c for c in cells if zlo < c[1] / 4.0 < zhi), float) / 4.0
    d = np.min(np.hypot(P[:, None, 0] - O[None, :, 0], P[:, None, 1] - O[None, :, 1]), axis=1) - 0.25
    k = int(np.argmax(d))
    return P[k], float(d[k])


holes = {}
for h in (1, 2):
    r = fh[fh[:, 0] == h]
    zmin, zmax = float(r[:, 1].min()), float(r[:, 1].max())
    span = zmax - zmin
    eu, ru = eye_circle(r, zmax - 0.3 * span, zmax + 1.0)
    el, rl = eye_circle(r, zmin - 1.0, zmin + 0.3 * span)
    tt, tb = r[r[:, 1] == zmax], r[r[:, 1] == zmin]
    tip_t = np.array([0.5 * (tt[:, 2].min() + tt[:, 3].max()), zmax + 0.25])
    tip_b = np.array([0.5 * (tb[:, 2].min() + tb[:, 3].max()), zmin - 0.25])
    nrow = r[np.abs(r[:, 1] - Z_NOTCH_C) <= 0.25 + 1e-9]
    if len(nrow) != 1:
        raise SystemExit("refusing: hole %d has %d runs on the notch row (Cannone z %.2f), expected the stem alone" % (h, len(nrow), Z_NOTCH_C))
    nx = float(nrow[0, 3]) if h == 1 else float(nrow[0, 2])  # the stem's inner edge (towards the centre line)
    holes[h] = dict(r=r, eye=eu, r_up=ru, beye=el, r_lo=rl, tips=(tip_t, tip_b), notch=np.array([nx, float(nrow[0, 1])]),
                    len_c=float(np.hypot(*(tip_t - tip_b))))


def rot(phi):
    return np.array([[math.cos(phi), -math.sin(phi)], [math.sin(phi), math.cos(phi)]])


for h, side in ((1, -1.0), (2, 1.0)):
    H_ = holes[h]
    H_["s"] = FH_LEN_P[h] / H_["len_c"]  # uniform scale of the Cannone shape to the printed length
    xu = side * 0.5 * (UP_INNER + UP_DIA[h])  # upper eye centre: inner edge + radius; ASSUMPTION: symmetric about the centre line
    xl = side * 0.25 * (LO_INNER + LO_OUTER)  # lower eye centre: mean of the inner- and outer-edge half spacings; same assumption
    v = H_["s"] * (H_["beye"] - H_["eye"])
    Lv, dx = float(np.hypot(*v)), xl - xu
    if abs(dx) >= Lv:
        raise SystemExit("refusing: hole %d printed eye offset %.2f mm exceeds the scaled eye-to-eye distance %.2f mm" % (h, dx, Lv))
    # rotate about the upper eye so the lower eye lands at the printed x, below the upper eye
    phi = -math.acos(dx / Lv) - math.atan2(v[1], v[0])
    H_["phi"] = math.atan2(math.sin(phi), math.cos(phi))
    zn = float((rot(H_["phi"]) @ (H_["s"] * (H_["notch"] - H_["eye"])))[1])
    H_["eye2"] = np.array([xu, Z_TOP - STOP_P[h] - zn])  # the notch at the printed stop below the top edge
cell = 0.5
fh_rows = []
for h in (1, 2):
    H_ = holes[h]
    r = H_["r"]
    R_inv = rot(-H_["phi"])
    # corners of the hole's box in cello coordinates, then every 0.5 mm cell centre inside it mapped back into the Cannone raster
    box = np.array([[r[:, 2].min() - 1, r[:, 1].min() - 1], [r[:, 3].max() + 1, r[:, 1].min() - 1], [r[:, 3].max() + 1, r[:, 1].max() + 1], [r[:, 2].min() - 1, r[:, 1].max() + 1]])
    box2 = np.array([H_["eye2"] + H_["s"] * rot(H_["phi"]) @ (p - H_["eye"]) for p in box])
    xs_new = np.arange(math.floor(box2[:, 0].min() / cell) * cell, box2[:, 0].max() + cell, cell)
    zs_new = np.arange(math.floor(box2[:, 1].min() / cell) * cell + 0.25, box2[:, 1].max() + cell, cell)
    rows_by_z = {}
    for _, z, xa, xb in r:
        rows_by_z.setdefault(round(z * 4), []).append((xa, xb))
    for z2 in zs_new:
        P2 = np.c_[xs_new + 0.25, np.full(len(xs_new), z2)]
        Pc = H_["eye"] + ((P2 - H_["eye2"]) @ R_inv.T) / H_["s"]
        zrow = np.round((Pc[:, 1] - 0.25) / 0.5) * 0.5 + 0.25  # the Cannone raster row (rows at k 0.5 + 0.25)
        occ = np.zeros(len(xs_new), bool)
        for i_ in range(len(xs_new)):
            for xa, xb in rows_by_z.get(round(zrow[i_] * 4), ()):
                if xa <= Pc[i_, 0] < xb:
                    occ[i_] = True
                    break
        k = 0
        while k < len(occ):
            if occ[k]:
                j = k
                while j < len(occ) and occ[j]:
                    j += 1
                fh_rows.append((h, z2, xs_new[k], xs_new[k] + cell * (j - k)))
                k = j
            else:
                k += 1
write("top-plate-fholes.tsv", "frankensim.specimen-raster-runs.v1", ["hole", "z_mm", "x_start_mm", "x_end_mm"], fh_rows, ["%d", "%.2f", "%.2f", "%.2f"])
fhr = np.array(fh_rows)
for h in (1, 2):
    rr_ = fhr[fhr[:, 0] == h]
    holes[h]["len2"] = float(rr_[:, 1].max() - rr_[:, 1].min() + cell)
    holes[h]["area2"] = float(np.sum(rr_[:, 3] - rr_[:, 2]) * cell)
    fw_ = lambda p, H_=holes[h]: H_["eye2"] + H_["s"] * rot(H_["phi"]) @ (p - H_["eye"])
    holes[h]["b2"] = fw_(holes[h]["beye"])
    holes[h]["n2"] = fw_(holes[h]["notch"])
    holes[h]["len_new"] = float(np.hypot(*(fw_(holes[h]["tips"][0]) - fw_(holes[h]["tips"][1]))))
    holes[h]["lo_edge"] = float(dist_to_poly(cpoly_t, holes[h]["b2"][None, :])[0]) - holes[h]["s"] * holes[h]["r_lo"]
H1, H2 = holes[1], holes[2]
LOCM_REF = "LOC item ihas.200154818 image 0021 (Montgomery 2005, Castelbarco cello, caliper mm)"
prov("f_hole_length", "%.1f / %.1f mm (bass / treble), end to end" % (FH_LEN_P[1], FH_LEN_P[2]), "LOCM",
     "%s; USED: each Cannone hole shape is scaled uniformly by printed / Cannone end-to-end length (%.2f / %.2f mm), scale %.4f / %.4f; placed length %.2f / %.2f mm" % (
         LOCM_REF, H1["len_c"], H2["len_c"], H1["s"], H2["s"], H1["len_new"], H2["len_new"]))
prov("f_hole_upper_eye_spacing", "%.1f mm between the upper eyes' inner edges (103.0 by tape)" % UP_INNER, "LOCM",
     "%s; USED: upper eye centre x = -/+(spacing + diameter)/2 = %.2f / %.2f mm; ASSUMPTION: symmetric about the centre line" % (LOCM_REF, H1["eye2"][0], H2["eye2"][0]))
prov("f_hole_upper_eye_diameter", "%.1f / %.1f mm (bass / treble)" % (UP_DIA[1], UP_DIA[2]), "LOCM",
     "%s; USED only to place the upper eye centre (inner edge + radius); the eyes are not resized: the scaled Cannone upper eyes are %.2f / %.2f mm" % (
         LOCM_REF, 2 * H1["s"] * H1["r_up"], 2 * H2["s"] * H2["r_up"]))
prov("f_hole_lower_eye_spacing_inner", "%.1f mm between the lower eyes' inner edges" % LO_INNER, "LOCM",
     "%s; USED with the outer-edge spacing: lower eye centre x = -/+(inner + outer)/4 = -/+%.2f mm, reached by rotating each hole about its upper eye (%.2f / %.2f deg); ASSUMPTION: symmetric about the centre line" % (
         LOCM_REF, 0.25 * (LO_INNER + LO_OUTER), math.degrees(H1["phi"]), math.degrees(H2["phi"])))
prov("f_hole_lower_eye_spacing_outer", "%.1f mm between the lower eyes' outer edges" % LO_OUTER, "LOCM",
     "%s; USED: see f_hole_lower_eye_spacing_inner (placed lower eye centres x %.2f / %.2f mm)" % (LOCM_REF, H1["b2"][0], H2["b2"][0]))
prov("f_hole_stop", "%.0f / %.0f mm from the top edge to the notches (bass / treble)" % (STOP_P[1], STOP_P[2]), "LOCM",
     "%s ('bass ff hole higher than treble'); USED: each hole is moved along z until the stem's inner edge on the Cannone notch row lies this far below the top edge (z %.2f): notches at z %.2f / %.2f mm; ASSUMPTION: the Cannone notch line = its bridge plane (Cannone z %.2f)" % (
         LOCM_REF, Z_TOP, H1["n2"][1], H2["n2"][1], Z_NOTCH_C))
prov("f_hole_lower_eye_diameter", "%.1f / %.1f mm (bass / treble)" % (LO_DIA[1], LO_DIA[2]), "LOCM",
     "%s; CHECK only, not used: the scaled Cannone lower eyes are %.2f / %.2f mm" % (LOCM_REF, 2 * H1["s"] * H1["r_lo"], 2 * H2["s"] * H2["r_lo"]))
prov("f_hole_notch_gap", "%.1f mm between the notches" % NOTCH_GAP_P, "LOCM",
     "%s; CHECK only, not used: the placed notch inner edges are %.2f mm apart" % (LOCM_REF, H2["n2"][0] - H1["n2"][0]))
prov("f_hole_lower_eye_to_edge", "%.1f / %.1f mm from the lower eye to the outside of the instrument (bass / treble)" % (LO_EDGE_P[1], LO_EDGE_P[2]), "LOCM",
     "%s; CHECK only, not used: placed lower eye edge to the top outline %.2f / %.2f mm" % (LOCM_REF, H1["lo_edge"], H2["lo_edge"]))
prov("f_hole_shape", "the il Cannone hole rasters (top-plate-fholes.tsv of the Cannone specimen), scaled, rotated and moved as in the rows above", "CAN",
     "ASSUMPTIONS: the Cannone hole shape; eye centres = the largest inscribed circles of the Cannone raster in the top and bottom 30 %% of each hole (Cannone radii %.2f, %.2f / %.2f, %.2f mm); hole ends = the centres of the first and last raster rows; the top's edge on the centre line half a 1 mm row beyond the last outline row" % (
         H1["r_up"], H1["r_lo"], H2["r_up"], H2["r_lo"]))
prov("f_holes","scale %.4f / %.4f; rotation %.2f / %.2f deg about the upper eye; length %.2f / %.2f mm; z extents %.1f / %.1f mm, areas %.0f / %.0f mm2; upper eyes (%.2f, %.2f) / (%.2f, %.2f), lower eyes (%.2f, %.2f) / (%.2f, %.2f), notches (%.2f, %.2f) / (%.2f, %.2f)" % (
    H1["s"], H2["s"], math.degrees(H1["phi"]), math.degrees(H2["phi"]), H1["len_new"], H2["len_new"], H1["len2"], H2["len2"], H1["area2"], H2["area2"],
    *H1["eye2"], *H2["eye2"], *H1["b2"], *H2["b2"], *H1["n2"], *H2["n2"]),
     "LOCM / CAN", "LOCM (printed on this cello): lengths 130.5 / 128.5 end to end; upper eyes 102.3 apart at their inner edges, diameters 13.5 / 13.0; lower eyes 220.0 apart at their inner edges and 258.2 at their outer edges; stop to the notches 419 / 420 from the top edge. Each Cannone hole shape is scaled uniformly to its printed length (Cannone end-to-end %.2f / %.2f mm), its upper eye centre put at x = -/+(102.3 + d_upper)/2, then rotated about that eye until the lower eye centre is at x = -/+(220.0 + 258.2)/4 = -/+119.55 below it, and moved along z until the stem's inner edge on the Cannone notch row lies at the printed stop below the top edge (z %.2f). ASSUMPTIONS: the Cannone hole shape; eye centres = the largest inscribed circles of the Cannone raster in the top and bottom 30 %% of each hole (Cannone radii %.2f, %.2f / %.2f, %.2f mm); hole ends = the centres of the first and last raster rows; the Cannone notch line = its bridge plane (Cannone z %.2f; the treble hole's inner nick is visible on the raster at z -106.75, the bass hole's is not resolved at 0.5 mm); the top's edge on the centre line half a 1 mm row beyond the last outline row; each spacing split symmetrically about the centre line" % (
    H1["len_c"], H2["len_c"], Z_TOP, H1["r_up"], H1["r_lo"], H2["r_up"], H2["r_lo"], Z_NOTCH_C))
prov("f_hole_check", "eye diameters upper %.2f / %.2f (LOCM 13.5 / 13.0), lower %.2f / %.2f (LOCM 19.7 / 19.2) mm; notch inner edges %.2f apart (LOCM 162.0); lower eye edge to the top outline %.2f / %.2f mm (LOCM 32.8 / 34.0 to the outside of the instrument); bridge plane z %.2f = %.2f mm below the top edge, %.2f / %.2f mm above the printed notch line (the bridge is not moved)" % (
    2 * H1["s"] * H1["r_up"], 2 * H2["s"] * H2["r_up"], 2 * H1["s"] * H1["r_lo"], 2 * H2["s"] * H2["r_lo"], H2["n2"][0] - H1["n2"][0], H1["lo_edge"], H2["lo_edge"],
    SZ * Z_NOTCH_C, Z_TOP - SZ * Z_NOTCH_C, STOP_P[1] - (Z_TOP - SZ * Z_NOTCH_C), STOP_P[2] - (Z_TOP - SZ * Z_NOTCH_C)), "check", "computed from the placed holes; not fitted")

# island between the holes for the top's 'centre at sound post' thickness
isl = []
for z2 in np.unique(fhr[:, 1]):
    a1 = fhr[(fhr[:, 0] == 1) & (fhr[:, 1] == z2)]
    a2 = fhr[(fhr[:, 0] == 2) & (fhr[:, 1] == z2)]
    if len(a1) and len(a2):
        isl.append((z2, a1[:, 3].max(), a2[:, 2].min()))
isl = np.array(isl)
ISLAND = np.vstack([np.c_[isl[:, 1], isl[:, 0]], np.c_[isl[::-1, 2], isl[::-1, 0]]])

# ------------------------------------------------------------------ thickness fields of the cello
_, _, mb_ = nt("back-plate-maps.tsv")
_, _, ps = nt("soundpost-sections.tsv")
PX, PZ = float(ps[:, 2].mean()), float(ps[:, 3].mean())
tb_interior = mb_[mb_[:, 5] > 15, 4]
TB_FLANK = float(np.percentile(tb_interior, 10))
TB_POST = float(FB.bil(FB.thick, PX, PZ))


def edge_w(d):
    return np.clip((LINING_W - d) / (LINING_W - PURF), 0.0, 1.0)


def t_top(xp, zp, d):
    di = dist_to_poly(ISLAND, np.c_[xp, zp])
    ins = inside_many(ISLAND, np.c_[xp, zp])
    w_isl = np.where(ins, 1.0, np.clip(1.0 - di / 20.0, 0.0, 1.0))
    t_int = T_TOP["flank"] + (T_TOP["post"] - T_TOP["flank"]) * w_isl
    we = edge_w(d)
    return we * T_TOP["edge"] + (1 - we) * t_int


def t_back(xp, zp, d):
    xc, zc = inv(xp, zp)
    tc = FB.bil(FB.thick, xc, zc)
    t_int = T_BACK["flank"] + (T_BACK["post"] - T_BACK["flank"]) * np.clip((tc - TB_FLANK) / (TB_POST - TB_FLANK), 0.0, 1.0)
    we = edge_w(d)
    return we * T_BACK["edge"] + (1 - we) * t_int


prov("thickness_top", "edge %.1f to %.0f mm from the outline, linear to the interior by %.0f mm; interior flank %.1f; island between the f-holes %.1f, blended to the flank over 20 mm" % (T_TOP["edge"], PURF, LINING_W, T_TOP["flank"], T_TOP["post"]),
     "JOH front", "edge 5.5, flanks = fluting 3.9, centre at sound post 4.7; centre of bouts 4.0 is within 0.1 mm of the flank and below the generator's 0.2 mm top zone step, so not drawn; ASSUMPTION: the 4.7 region is the island between the holes (bounded by the holes' inner edges over the z span where both holes are open) because the Cannone top shows no thickening at its post (3.08 mm within 6 mm of the post vs interior median 3.16); edge band width = purfling distance 6 and ramp end = lining width 16 (JOH) are ASSUMPTIONS of where the edge value applies")
prov("thickness_back", "edge %.1f to %.0f mm, linear by %.0f mm; interior = %.1f + %.1f x clip((t_Cannone - %.2f) / (%.2f - %.2f)) at the mapped Cannone point" % (T_BACK["edge"], PURF, LINING_W, T_BACK["flank"], T_BACK["post"] - T_BACK["flank"], TB_FLANK, TB_POST, TB_FLANK),
     "JOH back / CAN", "edge 5, flanks = fluting 3.5, centre at sound post 7; centre of bouts 3.6 within 0.1 mm of the flank (below the 0.4 mm back zone step); the Cannone back's measured thickening (5.9 mm at its post, %.2f mm = its interior 10th percentile at the flanks) sets WHERE the cello back is thick; ASSUMPTION: the cello back's thickness pattern is the Cannone's, rescaled to JOH's levels (Sacconi's 7 -> 6 -> 5 -> 4 centre contours in Hargrave 1987 p. 37 have the same central-strip form)" % TB_FLANK)

# ------------------------------------------------------------------ arching and rib height: solve the vertical scales
cring = np.column_stack(fwd(RING_C[:, 0], RING_C[:, 1]))
d_ring_t = dist_to_poly(cpoly_t, cring)
d_ring_b = dist_to_poly(cpoly_b, cring)
t_ring_t = t_top(cring[:, 0], cring[:, 1], d_ring_t)
t_ring_b = t_back(cring[:, 0], cring[:, 1], d_ring_b)
arch_t_c = E_T - FT.bil(FT.ymid, RING_C[:, 0], RING_C[:, 1])  # arch values at the ring (about 0)
arch_b_c = FB.bil(FB.ymid, RING_C[:, 0], RING_C[:, 1]) - E_B

# grid of the cello maps
GX = np.arange(math.floor(min(cpoly_t[:, 0].min(), cpoly_b[:, 0].min()) / 3) * 3 - 6 + 0.5, max(cpoly_t[:, 0].max(), cpoly_b[:, 0].max()) + 6, 3.0)
GZ = np.arange(math.floor(min(cpoly_t[:, 1].min(), cpoly_b[:, 1].min()) / 3) * 3 - 6 + 0.5, max(cpoly_t[:, 1].max(), cpoly_b[:, 1].max()) + 6, 3.0)
gx, gz = np.meshgrid(GX, GZ, indexing="ij")
G = np.c_[gx.ravel(), gz.ravel()]


def plate_grid(poly):
    d = dist_to_poly(poly, G)
    ins = inside_many(poly, G)
    keep = ins | (d <= 4.5)
    return G[keep], np.where(ins[keep], d[keep], 0.0)


Gt, dt = plate_grid(cpoly_t)
Gb, db = plate_grid(cpoly_b)
xt_c, zt_c = inv(Gt[:, 0], Gt[:, 1])
xb_c, zb_c = inv(Gb[:, 0], Gb[:, 1])
arch_t = E_T - FT.bil(FT.ymid, xt_c, zt_c)
arch_b = FB.bil(FB.ymid, xb_c, zb_c) - E_B
tt = t_top(Gt[:, 0], Gt[:, 1], dt)
tb = t_back(Gb[:, 0], Gb[:, 1], db)
# top: mid = -k_t arch; glue plane = mean inner skin at the ring = mean(-k_t arch_ring + t_ring/2); arch height = glue - min(outer)
# (outer = mid - t/2) = 29  ->  linear in k_t
in_tri = inside_many(cpoly_t, Gt)
A1 = float(np.mean(-arch_t_c)) ; B1 = float(np.mean(0.5 * t_ring_t))


def arch_height_top(k):
    return (k * A1 + B1) - float(np.min(-k * arch_t[in_tri] - 0.5 * tt[in_tri]))


lo_k, hi_k = 0.5, 6.0
for _ in range(80):
    mid_k = 0.5 * (lo_k + hi_k)
    lo_k, hi_k = (mid_k, hi_k) if arch_height_top(mid_k) < ARCH_TOP else (lo_k, mid_k)
K_T = 0.5 * (lo_k + hi_k)
in_bk = inside_many(cpoly_b, Gb)


def arch_height_back(k):  # back: mid = E + k arch; glue = mean(E + k arch_ring - t_ring/2); height = max(outer) - glue, outer = mid + t/2
    return float(np.max(k * arch_b[in_bk] + 0.5 * tb[in_bk])) - float(np.mean(k * arch_b_c - 0.5 * t_ring_b))


lo_k, hi_k = 0.5, 6.0
for _ in range(80):
    mid_k = 0.5 * (lo_k + hi_k)
    lo_k, hi_k = (mid_k, hi_k) if arch_height_back(mid_k) < ARCH_BACK else (lo_k, mid_k)
K_B = 0.5 * (lo_k + hi_k)
# rib height: mean over the ring of (inner_back - inner_top) = 120  ->  E_back'
inner_t_ring = -K_T * arch_t_c + 0.5 * t_ring_t
E_BP = RIB_H - float(np.mean(K_B * arch_b_c - 0.5 * t_ring_b)) + float(np.mean(inner_t_ring))
mid_t = -K_T * arch_t
mid_b = E_BP + K_B * arch_b
inner_b_ring = E_BP + K_B * arch_b_c - 0.5 * t_ring_b
rib_h = inner_b_ring - inner_t_ring
prov("arching", "top midsurface = -%.4f x (Cannone ring level %.2f - Cannone midsurface), back = %.2f + %.4f x (Cannone midsurface - its ring level %.2f)" % (K_T, E_T, E_BP, K_B, E_B),
     "JOH / CAN", "scales solved so that the outer surface's highest point lies JOH 29 (front) / 28 (back) mm above the plate's glue plane (mean inner skin at the rib ring); ASSUMPTION: the Cannone's arch shape, and JOH 'arching height' read as that height")
prov("rib_height", "mean %.2f mm (min %.2f, max %.2f) between the inner skins at the ring" % (rib_h.mean(), rib_h.min(), rib_h.max()), "JOH", "rib height average 120 (the variation along the ring is the Cannone's edge undulation x the arch scales)")
for name, Gp, mid, t, d, outward in (("top", Gt, mid_t, tt, dt, -1), ("back", Gb, mid_b, tb, db, +1)):
    outer = mid + outward * 0.5 * t
    inner = mid - outward * 0.5 * t
    rows = [(x, z, o, i, tk, dd, 36) for (x, z), o, i, tk, dd in zip(Gp, outer, inner, t, d)]
    write(name + "-plate-maps.tsv", "frankensim.specimen-plate-maps.v1", ["x_mm", "z_mm", "outer_y_mm", "inner_y_mm", "thickness_normal_mm", "dist_to_outline_mm", "valid_cells_of_36"], rows,
          ["%.1f", "%.1f", "%.3f", "%.3f", "%.3f", "%.2f", "%d"])

# ------------------------------------------------------------------ ribs (ring radii about the mapped centre at the same angles)
cx0, cz0 = rs[0, 1], rs[0, 2]
CX2, CZ2 = [float(v) for v in fwd(cx0, cz0)]


def polar_resample(poly):  # radius of a closed cello-plan polygon at the table's angles about (CX2, CZ2)
    out = []
    P = np.vstack([poly, poly[:1]])
    for t_ in th:
        u = np.array([math.cos(t_), sgn * math.sin(t_)])
        best = None
        for a, b in zip(P[:-1], P[1:]):
            e = b - a
            M = np.array([[u[0], -e[0]], [u[1], -e[1]]])
            if abs(np.linalg.det(M)) < 1e-12:
                continue
            s, q = np.linalg.solve(M, a - np.array([CX2, CZ2]))
            if s > 0 and -1e-9 <= q <= 1 + 1e-9 and (best is None or s < best):
                best = s
        out.append(best)
    return np.array(out)


rmin2 = polar_resample(np.column_stack(fwd(*ring_points(cx0, cz0, rs[:, 4], sgn).T)))
rmax2 = polar_resample(np.column_stack(fwd(*ring_points(cx0, cz0, rs[:, 5], sgn).T)))
# The stretched Cannone ring leaves the outlines at the C-bout corners (the Cannone's own -1.6 mm there, x 2.1). The
# generator moves such points in to 1.6 mm and asserts 0.8 mm after one smoothing pass, which the doubled corner moves
# break. Move only the points closer than 2.0 mm (the generator's 1.6 + 0.4) along the nearer outline's inward normal,
# so the generator has nothing left to move; everywhere else the stretched Cannone overhang is kept (mean about 5.3 mm,
# against JOH's 3 mm plate overhang + 0.85 mm half rib = 3.85: not imposed, named).
OVERHANG = 2.0
r_mid2 = 0.5 * (rmin2 + rmax2)
r_mid_in = r_mid2.copy()
U = np.c_[np.cos(th), sgn * np.sin(th)]


def ring_xy(r):
    rr_ = (np.roll(r, 1) + r + np.roll(r, -1)) / 3.0
    return np.c_[CX2 + rr_ * U[:, 0], CZ2 + rr_ * U[:, 1]]


def overhang(r):
    P = ring_xy(r)
    dt_ = np.where(inside_many(cpoly_t, P), 1, -1) * dist_to_poly(cpoly_t, P)
    db_ = np.where(inside_many(cpoly_b, P), 1, -1) * dist_to_poly(cpoly_b, P)
    return np.minimum(dt_, db_)


oh0 = overhang(r_mid_in)


def nearest_inward(poly, p):
    P = np.vstack([poly, poly[:1]])
    best = None
    for a, b in zip(P[:-1], P[1:]):
        ab = b - a
        L2 = float(ab @ ab)
        if L2 == 0:
            continue
        t = min(1.0, max(0.0, float((p - a) @ ab) / L2))
        q = a + t * ab
        dd = float(np.hypot(*(p - q)))
        if best is None or dd < best[0]:
            best = (dd, q)
    dd, q = best
    ins = inside_many(poly, p[None, :])[0]
    v = (p - q) / max(dd, 1e-9)
    return (dd if ins else -dd), (v if ins else -v)


# move the (dense) stretched ring polygon along the nearer outline's inward normal, then read it back at the table's angles
poly = np.column_stack(fwd(*ring_points(cx0, cz0, 0.5 * (rs[:, 4] + rs[:, 5]), sgn).T))
for _ in range(12):
    changed = False
    for i in range(len(poly)):
        sd_t, in_t = nearest_inward(cpoly_t, poly[i])
        sd_b, in_b = nearest_inward(cpoly_b, poly[i])
        sd, inw = (sd_t, in_t) if sd_t < sd_b else (sd_b, in_b)
        if sd < OVERHANG:
            poly[i] = poly[i] + inw * (OVERHANG + 0.1 - sd)
            changed = True
    if not changed:
        break
r_mid_in = polar_resample(poly)
for _ in range(40):  # the generator's 3-angle smoothing can pull a corner back out: small radial touch-ups
    oh = overhang(r_mid_in)
    bad = oh < OVERHANG - 0.3
    if not bad.any():
        break
    r_mid_in[bad] -= (OVERHANG - 0.2 - oh[bad])
oh1 = overhang(r_mid_in)
moved = r_mid2 - r_mid_in
half = 0.5 * (rmax2 - rmin2)
rmin2, rmax2 = r_mid_in - half, r_mid_in + half
prov("rib_ring_overhang", "stretched Cannone ring: overhang to the nearer outline mean %.2f, min %.2f mm; after pulling in: min %.2f mm; %d of %d angles moved, by up to %.2f mm" % (oh0.mean(), oh0.min(), oh1.min(), int((moved > 1e-9).sum()), len(moved), moved.max()),
     "generator / CAN", "points closer than 2.0 mm to either outline (the generator's own 1.6 mm minimum + 0.4) moved along the nearer outline's inward normal, then re-read at the table's angles; elsewhere the stretched Cannone overhang is kept (JOH's plate overhang 3 mm + half rib 0.85 mm = 3.85 is NOT imposed: the stretched Cannone mean is larger)")
rows = [(rs[i, 0], CX2, CZ2, RIB_T * RIB_H, rmin2[i], rmax2[i], float(inner_t_ring.mean()), float(inner_b_ring.mean()), RIB_T) for i in range(len(rs))]
write("ribs-sections.tsv", "frankensim.specimen-ring-sections.v1", ["theta_deg", "centre_x_mm", "centre_z_mm", "area_mm2", "r_min_mm", "r_max_mm", "y_min_mm", "y_max_mm", "mean_width_mm"], rows,
      ["%.0f", "%.3f", "%.3f", "%.2f", "%.3f", "%.3f", "%.2f", "%.2f", "%.2f"])
prov("ribs", "ring = the Cannone ring mapped through the plan scales, radii re-measured about the mapped centre (%.2f, %.2f) at the table's 2 degree angles; wall %.1f mm" % (CX2, CZ2, RIB_T), "JOH / CAN", "rib thickness 1.7 (JOH; written as every section's mean width, so the generator's 10th percentile is 1.7); y columns and area are not read by the generator")

# ------------------------------------------------------------------ bass bar
_, _, bs = nt("bass-bar-sections.tsv")
hb_c = bs[:, 7] - bs[:, 6]
kh = BAR_H / hb_c.max()
xb2, zb2 = fwd(bs[:, 2], bs[:, 0])
zmid = 0.5 * (zb2.min() + zb2.max())
klen = BAR_L / (zb2.max() - zb2.min())
zb3 = zmid + (zb2 - zmid) * klen
xb3 = sx(zb3 / SZ) * (bs[:, 2])  # centre line follows the plan map at the stretched z
# keep the bar under the bass foot where the Cannone has it: at the bridge plane, the bar's x relative to the feet is the
# Cannone's relative x scaled by the bridge's own x scale (the bridge is scaled about its feet's middle)
_, _, sil0 = nt("bridge-silhouette.tsv")
_, _, bsec0 = nt("bridge-sections.tsv")
yb0 = sil0[:, 0].max()
bot0 = sil0[sil0[:, 0] > yb0 - 3.0]
xfeet0 = 0.5 * (bot0[:, 1].min() + bot0[:, 2].max())
kxb0 = BR_FEET / (bot0[:, 2].max() - bot0[:, 1].min())
zbr_c = float(np.average(bsec0[:, 3], weights=bsec0[:, 1]))
bar_x_at_bridge_c = float(np.interp(zbr_c, bs[:, 0], bs[:, 2]))
bar_x_target = (bar_x_at_bridge_c - xfeet0) * kxb0
BAR_SHIFT = bar_x_target - float(np.interp(SZ * zbr_c, zb3, xb3))
xb3 = xb3 + BAR_SHIFT
h3 = hb_c * kh
rows = [(z, BAR_W * h, x, 0.0, x - 0.5 * BAR_W, x + 0.5 * BAR_W, 0.0, h) for z, h, x in zip(zb3, h3, xb3)]
write("bass-bar-sections.tsv", "frankensim.specimen-sections.v1", ["z_mm", "area_mm2", "centroid_x_mm", "centroid_y_mm", "x_min_mm", "x_max_mm", "y_min_mm", "y_max_mm"], rows, ["%.3f"] * 8)
keep = h3 > 2.0
prov("bass_bar", "length %.1f mm of sections (z %.1f to %.1f), %.1f mm above 2 mm height; max height %.1f; wall %.1f mm; centre line x %.1f to %.1f" % (zb3.max() - zb3.min(), zb3.min(), zb3.max(), zb3[keep].max() - zb3[keep].min(), h3.max(), BAR_W, xb3.min(), xb3.max()),
     "JOH / CAN", "JOH bar length approx. 580 (the Cannone's whole section range stretched to it about its mapped middle), height middle approx. 24 (the Cannone height profile x %.3f), thickness approx. 11; ASSUMPTION: Cannone bar profile and line (centre line x = the plan width scale at each z, then shifted %.1f mm in x so that at the bridge plane it sits under the bass foot at the Cannone's position across the foot, x %.1f mm)" % (kh, BAR_SHIFT, bar_x_target))

# ------------------------------------------------------------------ soundpost
P2x, P2z = [float(v) for v in fwd(PX, PZ)]
y_top_post = float(-K_T * (E_T - FT.bil(FT.ymid, PX, PZ)))
y_back_post = float(E_BP + K_B * (FB.bil(FB.ymid, PX, PZ) - E_B))
ys = np.linspace(y_top_post + 6, y_back_post - 6, 28)
rows = [(y, math.pi * (0.5 * POST_D) ** 2, P2x, P2z, P2x - 0.5 * POST_D, P2x + 0.5 * POST_D, P2z - 0.5 * POST_D, P2z + 0.5 * POST_D) for y in ys]
write("soundpost-sections.tsv", "frankensim.specimen-sections.v1", ["y_mm", "area_mm2", "centroid_x_mm", "centroid_z_mm", "x_min_mm", "x_max_mm", "z_min_mm", "z_max_mm"], rows, ["%.3f"] * 8)
prov("soundpost", "diameter %.1f mm, vertical axis at (%.2f, %.2f)" % (POST_D, P2x, P2z), "JOH / CAN", "JOH soundpost thickness 11; ASSUMPTION: position = the Cannone post centroid (16.6, -114.3) through the plan map, axis vertical (the Cannone's tilt 0.034 / 0.008 not carried)")

# ------------------------------------------------------------------ bridge
_, _, sil = nt("bridge-silhouette.tsv")
_, _, bsec = nt("bridge-sections.tsv")
ybot, ytop = sil[:, 0].max(), sil[:, 0].min()
bottom = sil[sil[:, 0] > ybot - 3.0]
feet_c = bottom[:, 2].max() - bottom[:, 1].min()
kx_b = BR_FEET / feet_c
ky_b = BR_H / (ybot - ytop + 0.25)
zc_c = float(np.average(bsec[:, 3], weights=bsec[:, 1]))
ZC2 = float(SZ * zc_c)
foot_secs = bsec[(bsec[:, 0] > -25.5) & (bsec[:, 0] < -21.5)]
t_foot_c = float(np.mean(foot_secs[:, 7] - foot_secs[:, 6]))
area_c = float(np.sum(sil[:, 2] - sil[:, 1]) * 0.25)
vol_c = 3170.3
t_unif_c = vol_c / area_c
T_BR = t_unif_c * BR_T_FOOT / t_foot_c
xfeet = 0.5 * (bottom[:, 1].min() + bottom[:, 2].max())  # the Cannone feet's middle in x
# feet bottoms onto the cello top's outer skin: mean of the outer skin under the two feet at the bridge plane
xs_feet = np.linspace((bottom[:, 1].min() - xfeet) * kx_b, (bottom[:, 2].max() - xfeet) * kx_b, 41)
xs_c, zs_c = inv(xs_feet, np.full_like(xs_feet, ZC2))
outer_feet = -K_T * (E_T - FT.bil(FT.ymid, xs_c, zs_c)) - 0.5 * t_top(xs_feet, np.full_like(xs_feet, ZC2), dist_to_poly(cpoly_t, np.c_[xs_feet, np.full_like(xs_feet, ZC2)]))
YBOT2 = float(np.max(outer_feet))  # the lowest point (largest y) the feet reach
# Each Cannone foot is fitted to the Cannone arch (bass foot bottom about 1 mm lower than the treble foot); scaled 2.65x in
# height that offset would be about 3 mm and the generator (feet = runs in the bottom 3 mm of rows) would see one foot.
# So each foot's bottom is put on the cello's outer skin under that foot's centre: the inverse map is
#   x_c = xfeet + x' / kx,  y_c = yfoot(x_c) - (Yskin(x') - y') / ky,
# yfoot and Yskin linear in x between the two feet (a shear of the whole silhouette, a few tenths of a degree).
legs = sil[sil[:, 0] > ybot - 6.0]
foot_c = []
for lo_, hi_ in ((-1e9, xfeet), (xfeet, 1e9)):
    rr_ = legs[(0.5 * (legs[:, 1] + legs[:, 2]) > lo_) & (0.5 * (legs[:, 1] + legs[:, 2]) <= hi_)]
    yb_ = rr_[:, 0].max()
    at = rr_[rr_[:, 0] == yb_]
    foot_c.append((float(0.5 * (at[:, 1].min() + at[:, 2].max())), float(yb_) + 0.125))
(fxl, fyl), (fxr, fyr) = foot_c
fxl2, fxr2 = (fxl - xfeet) * kx_b, (fxr - xfeet) * kx_b
sk = []
for fx2 in (fxl2, fxr2):
    xc_, zc_ = inv(np.array([fx2]), np.array([ZC2]))
    sk.append(float(-K_T * (E_T - FT.bil(FT.ymid, xc_, zc_))[0] - 0.5 * t_top(np.array([fx2]), np.array([ZC2]), dist_to_poly(cpoly_t, np.array([[fx2, ZC2]])))[0]))
skl, skr = sk


def yfoot(xc):
    return fyl + (fyr - fyl) * (xc - fxl) / (fxr - fxl)


def yskin(x2):
    return skl + (skr - skl) * (x2 - fxl2) / (fxr2 - fxl2)


rows_b = []
xs2 = np.arange(math.floor((sil[:, 1].min() - xfeet) * kx_b / 0.25) * 0.25 - 1.0, (sil[:, 2].max() - xfeet) * kx_b + 1.0, 0.25)
xc_all = xfeet + (xs2 + 0.125) / kx_b
rows_by_y = {}
for y_, a_, b_ in sil:
    rows_by_y.setdefault(int(round(y_ * 4)), []).append((a_, b_))
ys_new = np.arange(math.floor((min(skl, skr) - BR_H - 2.0) / 0.25) * 0.25, max(skl, skr) + 1.0, 0.25)
for y2 in ys_new:
    yc = yfoot(xc_all) - (yskin(xs2 + 0.125) - (y2 + 0.125)) / ky_b
    occ = np.zeros(len(xs2), bool)
    for i_, (xc_, ycv) in enumerate(zip(xc_all, yc)):
        for a_, b_ in rows_by_y.get(int(round(math.floor(ycv * 4 + 0.5) )), ()):
            if a_ <= xc_ < b_:
                occ[i_] = True
                break
    k = 0
    while k < len(occ):
        if occ[k]:
            j = k
            while j < len(occ) and occ[j]:
                j += 1
            rows_b.append((y2 + 0.125, xs2[k], xs2[k] + 0.25 * (j - k)))
            k = j
        else:
            k += 1
YBOT2 = max(r[0] for r in rows_b)
write("bridge-silhouette.tsv", "frankensim.specimen-raster-runs.v1", ["y_mm", "x_start_mm", "x_end_mm"], rows_b, ["%.2f", "%.2f", "%.2f"])
rb = np.array(rows_b)
area2 = float(np.sum(rb[:, 2] - rb[:, 1]) * 0.25)
write("bridge-sections.tsv", "frankensim.specimen-sections.v1", ["y_mm", "area_mm2", "centroid_x_mm", "centroid_z_mm", "x_min_mm", "x_max_mm", "z_min_mm", "z_max_mm"],
      [(YBOT2 - 0.5 * BR_H, 100.0, 0.0, ZC2, -0.5 * BR_FEET, 0.5 * BR_FEET, ZC2 - 0.5 * T_BR, ZC2 + 0.5 * T_BR)], ["%.3f"] * 8)
write("part-volumes.tsv", "frankensim.specimen-part-volumes.v1", ["part", "volume_mm3", "triangles", "x_min_mm", "x_max_mm", "y_min_mm", "y_max_mm", "z_min_mm", "z_max_mm"],
      [("bridge", T_BR * area2, 0, rb[:, 1].min(), rb[:, 2].max(), rb[:, 0].min(), rb[:, 0].max(), ZC2 - 0.5 * T_BR, ZC2 + 0.5 * T_BR)], ["%s", "%.1f", "%d", "%.2f", "%.2f", "%.2f", "%.2f", "%.2f", "%.2f"])
prov("bridge", "silhouette x %.4f (feet %.1f -> 90 mm), y %.4f (height %.2f -> 92 mm), plane z %.2f, feet bottoms on the outer skin at y %.2f (bass) / %.2f (treble), Cannone feet bottoms %.2f / %.2f; uniform shell thickness %.2f mm; area %.0f mm2; volume %.0f mm3" % (kx_b, feet_c, ky_b, ybot - ytop + 0.25, ZC2, skl, skr, fyl, fyr, T_BR, area2, T_BR * area2),
     "JOH / CAN", "JOH bridge feet width 90, height approx. 92, thickness 11; thickness = the Cannone's uniform shell value (volume 3170 mm3 / silhouette %.0f mm2 = %.2f mm) x 11 / the Cannone's foot thickness %.2f mm (z extent of its sections 1-3 mm above the feet); ASSUMPTIONS: the Cannone silhouette's shape, JOH 'thickness 11' read as the foot thickness, the Cannone's mean-to-foot thickness ratio; position = the Cannone bridge plane through the length scale (stop %.1f mm from the top's neck-end edge vs JOH stoplength 405)" % (area_c, t_unif_c, t_foot_c, top_rows[:, 0].max() - ZC2))

# ------------------------------------------------------------------ strings at the bridge
xs_str = np.array([-1.5, -0.5, 0.5, 1.5]) * STR_SPAN / 3.0
rows_s = []
for name, x in zip("cgda", xs_str):
    occ = rb[(rb[:, 1] <= x) & (rb[:, 2] >= x)]
    y_top = float(occ[:, 0].min()) - 0.125
    rows_s.append((name, "bridge", ZC2, x, y_top, 0.0, 0.0))
write("string-contacts.tsv", "frankensim.specimen-string-contacts.v1", ["string", "contact", "z_mm", "x_mm", "y_mm", "slope_dy_dz_before", "slope_dy_dz_after"], rows_s,
      ["%s", "%s", "%.2f", "%.3f", "%.3f", "%.1f", "%.1f"])
prov("strings_at_bridge", " ".join("%s (%.2f, %.2f)" % (r[0], r[3], r[4]) for r in rows_s), "JOH / CAN", "JOH string spacing at bridge 47 (equal gaps, C on the bass-bar side, -x); y = the top of the scaled silhouette at that x (notch depth not modelled); slopes not used by the generator")

# ------------------------------------------------------------------ dimensions the generator reads
write("dimensions.tsv", "frankensim.specimen-dimensions.v1", ["part", "feature", "quantity", "value_si", "unit", "sigma_si", "route", "source_id", "note"],
      [("soundpost", "section", "equivalent-diameter-mean", "%.6f" % (POST_D / 1000), "m", "unstated", "published-maker-table", "johannsson-cello-1-1", "soundpost thickness 11 mm (hansjohannsson.com/resources/mescello.html, 1/1 column)")],
      ["%s"] * 9)
with open(OUT + "PROVENANCE.tsv", "x", encoding="utf-8") as f:
    f.write("item\tvalue\tsource\tnote\n")
    for r in PROV:
        f.write("\t".join(r) + "\n")
for r in PROV:
    print(r[0], "=", r[1])
