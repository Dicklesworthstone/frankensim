"""A DOUBLE BASS specimen directory for the shell-assembly generator gen_violin_shell_assembly.py, made by the same mapping as the cello
specimen script: the il Cannone specimen tables (ESRF CT, the h6 violin's inputs) mapped to bass dimensions from a
published maker table. This file is that script with the targets below and the changes listed in PROVENANCE.tsv (f-holes by total
area, thickness interior, rib height, bridge). Every mapped quantity is written with its source or labelled ASSUMPTION in
PROVENANCE.tsv beside the tables. Nothing here is tuned to any recording or mode frequency.

usage: python3 -I bass_specimen.py CANNONE_DIR OUT_DIR   (OUT_DIR must not exist)

Sources (values in mm unless stated):
  JOH   = Hans Johannsson, bass specification table, hansjohannsson.com/resources/mesbass.html, the single column "BASS mm."
          (the copy read has md5 e9fd0c688c24fc8525314876865df672, 15084 bytes; the page says "This table is being modified and revised"):
          body length 1130; widths upper ca. 530, middle ca. 380, lower ca. 650; rib height average 220-160, rib thickness 2.5;
          arching back 44 / front 48; back thickness: centre at sound post 12, flanks 5, fluting 5.5, centre of bouts 8, edge 7.5;
          front thickness: centre at sound post 8, flanks 5, fluting 6, centre of bouts 7, edge 7.5; overhang 3.5;
          purfling distance from edge 7.5; lining width 25; between top holes 170, bottom hole to edge 50 (soundhole LENGTH blank);
          bass bar length approx. 850, thickness approx. 20 (height 'approx. 450' not used); soundpost thickness 18;
          bridge feet width 155, thickness 20, height approx. 157; string spacing at bridge 85 (4-string); string length approx. 1070.
          'stoplength 6250' and 'necklength 4360' are not used (the stop is taken from the string length).
  MB    = the mapping rules fixed before any render: the Cannone outline scaled to Johannsson's 1130 length and
          530/380/650 bouts, labelled as such; f-holes scaled so the total area is BRN's 101 cm2; bar height 35 = the cello's
          24 x 850/580; ribs 220->160 if the generator takes a taper, else 190 = the midpoint (said); stop from the string length;
          bridge 155 x 157 with a 20 mm foot, a published bass bridge mass if found, else the uniform shell, labelled.
  BRN   = Brown 2004 (thesis, PDF md5 1578d330), pp. 56-57: total f-hole area 101 cm2 (Fischbach 3/4 test bass, a different bass).
  CAN   = the il Cannone tables themselves (shape only: outline, arch shape, the back's thickness pattern, hole shapes,
          bar profile, bridge silhouette, positions), mapped as stated.
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
L_BACK, B_UP, B_C, B_LO = 1130.0, 530.0, 380.0, 650.0  # JOH body length and the three widths ("ca." as printed); the Cannone outline is scaled to these
RIB_H, RIB_T = 190.0, 2.5  # JOH rib height average 220-160: the MIDPOINT 190 (the generator reads rib height only through the two plates' maps; a taper would be an unpublished tilt of the back); thickness 2.5
ARCH_TOP, ARCH_BACK = 48.0, 44.0  # JOH
T_TOP = dict(post=8.0, flank=5.0, flute=6.0, bout=7.0, edge=7.5)  # JOH front
T_BACK = dict(post=12.0, flank=5.0, flute=5.5, bout=8.0, edge=7.5)  # JOH back
# ASSUMPTION: the bass table's flanks / fluting / centre-of-bouts differ by 1-3 mm and the table does not say where they apply; the
# interior base value of each plate is the arithmetic mean of the three printed stations. The cello's three stations agreed within 0.1 mm.
T_TOP_INT = (T_TOP["flank"] + T_TOP["flute"] + T_TOP["bout"]) / 3.0
T_BACK_INT = (T_BACK["flank"] + T_BACK["flute"] + T_BACK["bout"]) / 3.0
PURF, LINING_W = 7.5, 25.0  # JOH purfling distance from edge 7.5, lining width 25
FH_AREA_TOTAL, FH_EYES = 10100.0, 170.0  # BRN total f-hole area 101 cm2; JOH between top holes 170. JOH f-hole LENGTH is blank.
BAR_L, BAR_W, BAR_H = 850.0, 20.0, 35.0  # JOH length approx. 850, thickness approx. 20; height 35 = the cello's 24 x 850/580 ('approx. 450' not used)
POST_D = 18.0  # JOH
BR_FEET, BR_T_FOOT, BR_H = 155.0, 20.0, 157.0  # JOH
STR_SPAN = 85.0  # JOH string spacing at the bridge, 4-string (outer string to outer string)

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
prov("plan_length_scale", "%.6f" % SZ, "JOH / MB 16036 / CAN", "z' = %.4f z: JOH body length 1130 mm over the Cannone body length %.1f mm (top outline rows at least 10 mm wide, i.e. without the back's button)" % (SZ, L_CAN))
prov("plan_width_scale", "lower %.4f at z %.2f, C %.4f at z %.2f, upper %.4f at z %.2f (Cannone z)" % (ST_S[0], ST_Z[0], ST_S[1], ST_Z[1], ST_S[2], ST_Z[2]), "JOH / MB 16036 / CAN",
     "x' = s(z) x, s piecewise linear in the Cannone z through the three bout stations of the Cannone back (widths %.1f, %.1f, %.1f mm -> JOH 650, 380, 530 'ca.'), constant beyond. LABEL (MB 16036): violin-form outline scaled to the bass table; no measured bass outline reachable (Brown 2004 prints back length, bout widths and rib depth only; PlumCanyon 16029); many basses have sloping viol shoulders, which this outline does not have" % (wb[kL], wb[kC], wb[kU]))


def sx(zc):  # Cannone z -> width scale
    return np.interp(zc, ST_Z, ST_S)


def fwd(x, z):  # Cannone plan -> bass plan
    x, z = np.asarray(x, float), np.asarray(z, float)
    return sx(z) * x, SZ * z


def inv(xp, zp):  # bass plan -> Cannone plan
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

# ------------------------------------------------------------------ bass plan outlines
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
ZMID_B = -87.0 * SZ  # the Cannone body centre (z -87, see CENTRE in the generator) in bass units
prov("outline_check", "back length %.1f, bouts lower %.1f / C %.1f / upper %.1f mm" % (back_rows[:, 0].max() - back_rows[:, 0].min() + SZ, wcb[back_rows[:, 0] < ZMID_B].max(),
     wcb[(back_rows[:, 0] > ST_Z[0] * SZ) & (back_rows[:, 0] < ST_Z[2] * SZ)].min(), wcb[back_rows[:, 0] > ZMID_B].max()), "check", "the mapped back outline against JOH 1130 / 650 / 380 / 530 (length includes the button rows that are at least 10 mm wide in Cannone units)")


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
_, _, fh = nt("top-plate-fholes.tsv")
holes = {}
for h in (1, 2):
    r = fh[fh[:, 0] == h]
    w = r[:, 3] - r[:, 2]
    cz = float(np.sum(r[:, 1] * w) / w.sum())
    cx = float(np.sum(0.5 * (r[:, 2] + r[:, 3]) * w) / w.sum())
    zmin, zmax = r[:, 1].min(), r[:, 1].max()
    top = r[r[:, 1] > zmax - 6]
    ex = float(np.sum(0.5 * (top[:, 2] + top[:, 3]) * (top[:, 3] - top[:, 2])) / np.sum(top[:, 3] - top[:, 2]))
    ez = float(np.sum(top[:, 1] * (top[:, 3] - top[:, 2])) / np.sum(top[:, 3] - top[:, 2]))
    bot = r[r[:, 1] < zmin + 6]
    bx = float(np.sum(0.5 * (bot[:, 2] + bot[:, 3]) * (bot[:, 3] - bot[:, 2])) / np.sum(bot[:, 3] - bot[:, 2]))
    bz = float(np.sum(bot[:, 1] * (bot[:, 3] - bot[:, 2])) / np.sum(bot[:, 3] - bot[:, 2]))
    holes[h] = dict(r=r, cx=cx, cz=cz, len=zmax - zmin + 0.5, eye=np.array([ex, ez]), beye=np.array([bx, bz]), area=float(np.sum(w) * 0.5))
# JOH's soundhole length is blank, so the uniform scale is set by the TOTAL area of both holes = BRN's 101 cm2 (not by a length)
AREA_C = holes[1]["area"] + holes[2]["area"]
SF = math.sqrt(FH_AREA_TOTAL / AREA_C)
FH_EDGE = 50.0  # JOH bottom hole to edge


def rot(phi):
    return np.array([[math.cos(phi), -math.sin(phi)], [math.sin(phi), math.cos(phi)]])


for h, side in ((1, -1.0), (2, 1.0)):
    H_ = holes[h]
    H_["eye2"] = np.array([side * 0.5 * FH_EYES, SZ * H_["eye"][1]])  # top eye at -85 / +85 mm, at the mapped Cannone eye z

    def bdist(phi, H_=H_):
        b2 = H_["eye2"] + SF * rot(phi) @ (H_["beye"] - H_["eye"])
        return float(dist_to_poly(cpoly_t, b2[None, :])[0]), b2
    # rotate the hole about its top eye (bottom moving outward) until the bottom eye is JOH 50 mm from the outline
    cands = [p for p in np.radians(np.linspace(-45, 45, 361))]
    best = min(cands, key=lambda p: (abs(bdist(p)[0] - FH_EDGE) if inside_many(cpoly_t, bdist(p)[1][None, :])[0] else 1e9, abs(p)))
    a_, b_ = best - math.radians(0.25), best + math.radians(0.25)
    for _ in range(40):
        m_ = 0.5 * (a_ + b_)
        if (bdist(a_)[0] - FH_EDGE) * (bdist(m_)[0] - FH_EDGE) <= 0:
            b_ = m_
        else:
            a_ = m_
    H_["phi"] = 0.5 * (a_ + b_)
cell = 0.5
fh_rows = []
for h in (1, 2):
    H_ = holes[h]
    r = H_["r"]
    R_inv = rot(-H_["phi"])
    # corners of the hole's box in bass-plan coordinates, then every 0.5 mm cell centre inside it mapped back into the Cannone raster
    box = np.array([[r[:, 2].min() - 1, r[:, 1].min() - 1], [r[:, 3].max() + 1, r[:, 1].min() - 1], [r[:, 3].max() + 1, r[:, 1].max() + 1], [r[:, 2].min() - 1, r[:, 1].max() + 1]])
    box2 = np.array([H_["eye2"] + SF * rot(H_["phi"]) @ (p - H_["eye"]) for p in box])
    xs_new = np.arange(math.floor(box2[:, 0].min() / cell) * cell, box2[:, 0].max() + cell, cell)
    zs_new = np.arange(math.floor(box2[:, 1].min() / cell) * cell + 0.25, box2[:, 1].max() + cell, cell)
    rows_by_z = {}
    for _, z, xa, xb in r:
        rows_by_z.setdefault(round(z * 4), []).append((xa, xb))
    for z2 in zs_new:
        P2 = np.c_[xs_new + 0.25, np.full(len(xs_new), z2)]
        Pc = H_["eye"] + ((P2 - H_["eye2"]) @ R_inv.T) / SF
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
    holes[h]["b2"] = holes[h]["eye2"] + SF * rot(holes[h]["phi"]) @ (holes[h]["beye"] - holes[h]["eye"])
    holes[h]["bottom_eye_to_edge"] = float(dist_to_poly(cpoly_t, holes[h]["b2"][None, :])[0])
prov("f_holes", "scale %.4f; rotation %.1f / %.1f deg about the top eye; z extents %.1f / %.1f mm, areas %.0f / %.0f mm2; top eyes (%.1f, %.1f) / (%.1f, %.1f), bottom eyes (%.1f, %.1f) / (%.1f, %.1f)" % (
    SF, math.degrees(holes[1]["phi"]), math.degrees(holes[2]["phi"]), holes[1]["len2"], holes[2]["len2"], holes[1]["area2"], holes[2]["area2"], *holes[1]["eye2"], *holes[2]["eye2"], *holes[1]["b2"], *holes[2]["b2"]),
     "BRN / JOH / MB 16036 / CAN", "each Cannone hole shape scaled uniformly (one factor %.4f for both) so that the TOTAL area of both holes is BRN's 101 cm2 = 10100 mm2 (Brown 2004 pp. 56-57, measured on the Fischbach 3/4 test bass, a different bass; JOH's soundhole length is blank; Cannone raster areas %.0f + %.0f mm2, lengths %.1f and %.1f mm); top eyes (centroid of the top 6 Cannone-mm) JOH 170 mm apart at the mapped Cannone eye z; each hole then rotated about its top eye until its bottom eye (centroid of the bottom 6 Cannone-mm) is JOH 50 mm from the top outline (nearest point); ASSUMPTIONS: Cannone hole shape, the eye centroids, rotation as the way to meet 170 and 50 together, the area of a different bass" % (SF, holes[1]["area"], holes[2]["area"], holes[1]["len"], holes[2]["len"]))
prov("f_hole_check", "bottom eye to the top outline %.1f / %.1f mm; raster areas %.0f + %.0f = %.0f mm2 (target 10100)" % (holes[1]["bottom_eye_to_edge"], holes[2]["bottom_eye_to_edge"], holes[1]["area2"], holes[2]["area2"], holes[1]["area2"] + holes[2]["area2"]), "check", "JOH bottom hole to edge 50; BRN total area 101 cm2")

# island between the holes for the top's 'centre at sound post' thickness
isl = []
for z2 in np.unique(fhr[:, 1]):
    a1 = fhr[(fhr[:, 0] == 1) & (fhr[:, 1] == z2)]
    a2 = fhr[(fhr[:, 0] == 2) & (fhr[:, 1] == z2)]
    if len(a1) and len(a2):
        isl.append((z2, a1[:, 3].max(), a2[:, 2].min()))
isl = np.array(isl)
ISLAND = np.vstack([np.c_[isl[:, 1], isl[:, 0]], np.c_[isl[::-1, 2], isl[::-1, 0]]])

# ------------------------------------------------------------------ thickness fields of the bass
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
    t_int = T_TOP_INT + (T_TOP["post"] - T_TOP_INT) * w_isl
    we = edge_w(d)
    return we * T_TOP["edge"] + (1 - we) * t_int


def t_back(xp, zp, d):
    xc, zc = inv(xp, zp)
    tc = FB.bil(FB.thick, xc, zc)
    t_int = T_BACK_INT + (T_BACK["post"] - T_BACK_INT) * np.clip((tc - TB_FLANK) / (TB_POST - TB_FLANK), 0.0, 1.0)
    we = edge_w(d)
    return we * T_BACK["edge"] + (1 - we) * t_int


prov("thickness_top", "edge %.1f to %.0f mm from the outline, linear to the interior by %.0f mm; interior base %.2f (mean of flanks %.1f, fluting %.1f, centre of bouts %.1f); island between the f-holes %.1f, blended to the interior over 20 mm" % (T_TOP["edge"], PURF, LINING_W, T_TOP_INT, T_TOP["flank"], T_TOP["flute"], T_TOP["bout"], T_TOP["post"]),
     "JOH front", "JOH front: edge 7.5, flanks 5, fluting 6, centre of bouts 7, centre at sound post 8. ASSUMPTIONS (the table gives values, not where they apply): (a) the interior base is the arithmetic mean of the three printed interior stations (the cello's three agreed within 0.1 mm; the bass's differ by up to 2 mm, so this is a guess of the area-weighted level, not a pattern: fluting and bout-centre are NOT drawn as separate zones); (b) the 8 region is the island between the holes (bounded by the holes' inner edges over the z span where both holes are open), as the cello; (c) edge band width = purfling distance 7.5 and ramp end = lining width 25 (JOH) are ASSUMPTIONS of where the edge value applies; (d) the 20 mm island blend is the cello script's, not scaled")
prov("thickness_back", "edge %.1f to %.0f mm, linear by %.0f mm; interior = %.2f + %.2f x clip((t_Cannone - %.2f) / (%.2f - %.2f)) at the mapped Cannone point" % (T_BACK["edge"], PURF, LINING_W, T_BACK_INT, T_BACK["post"] - T_BACK_INT, TB_FLANK, TB_POST, TB_FLANK),
     "JOH back / CAN", "JOH back: edge 7.5, flanks 5, fluting 5.5, centre of bouts 8, centre at sound post 12. ASSUMPTIONS: (a) interior base = mean of flanks, fluting and centre of bouts (%.2f), fluting and bout-centre not drawn separately; (b) the Cannone back's measured thickening (5.9 mm at its post, %.2f mm = its interior 10th percentile at the flanks) sets WHERE the bass back is thick, rescaled to JOH's levels up to the 12 mm at the post, as the cello (Sacconi's 7 -> 6 -> 5 -> 4 centre contours in Hargrave 1987 p. 37 have the same central-strip form for the cello; no bass contour map in reach)" % (T_BACK_INT, TB_FLANK))

# ------------------------------------------------------------------ arching and rib height: solve the vertical scales
cring = np.column_stack(fwd(RING_C[:, 0], RING_C[:, 1]))
d_ring_t = dist_to_poly(cpoly_t, cring)
d_ring_b = dist_to_poly(cpoly_b, cring)
t_ring_t = t_top(cring[:, 0], cring[:, 1], d_ring_t)
t_ring_b = t_back(cring[:, 0], cring[:, 1], d_ring_b)
arch_t_c = E_T - FT.bil(FT.ymid, RING_C[:, 0], RING_C[:, 1])  # arch values at the ring (about 0)
arch_b_c = FB.bil(FB.ymid, RING_C[:, 0], RING_C[:, 1]) - E_B

# grid of the bass maps
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
     "JOH / CAN", "scales solved so that the outer surface's highest point lies JOH 48 (front) / 44 (back) mm above the plate's glue plane (mean inner skin at the rib ring); ASSUMPTION: the Cannone's arch shape (a violin-form arch, with a flat or sloping bass back not represented), and JOH 'arching height' read as that height")
prov("rib_height", "mean %.2f mm (min %.2f, max %.2f) between the inner skins at the ring" % (rib_h.mean(), rib_h.min(), rib_h.max()), "JOH / MB 16036", "JOH rib height average '220-160' (a taper): the generator reads rib height only as the distance between the two plates' inner skins at the ring, so a taper would need a tilted back (an unpublished direction and slope); MB 16036: 'if not, use 190 and say it is the midpoint': the MIDPOINT 190 is used, same mean height, same enclosed volume to first order, no taper (the variation along the ring is the Cannone's edge undulation x the arch scales)")
for name, Gp, mid, t, d, outward in (("top", Gt, mid_t, tt, dt, -1), ("back", Gb, mid_b, tb, db, +1)):
    outer = mid + outward * 0.5 * t
    inner = mid - outward * 0.5 * t
    rows = [(x, z, o, i, tk, dd, 36) for (x, z), o, i, tk, dd in zip(Gp, outer, inner, t, d)]
    write(name + "-plate-maps.tsv", "frankensim.specimen-plate-maps.v1", ["x_mm", "z_mm", "outer_y_mm", "inner_y_mm", "thickness_normal_mm", "dist_to_outline_mm", "valid_cells_of_36"], rows,
          ["%.1f", "%.1f", "%.3f", "%.3f", "%.3f", "%.2f", "%d"])

# ------------------------------------------------------------------ ribs (ring radii about the mapped centre at the same angles)
cx0, cz0 = rs[0, 1], rs[0, 2]
CX2, CZ2 = [float(v) for v in fwd(cx0, cz0)]


def polar_resample(poly):  # radius of a closed bass-plan polygon at the table's angles about (CX2, CZ2)
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
prov("ribs", "ring = the Cannone ring mapped through the plan scales, radii re-measured about the mapped centre (%.2f, %.2f) at the table's 2 degree angles; wall %.1f mm" % (CX2, CZ2, RIB_T), "JOH / CAN", "rib thickness 2.5 (JOH; written as every section's mean width, so the generator's 10th percentile is 2.5); y columns and area are not read by the generator")

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
# Found by the generator's refusal (chain 'bar' of the top is not inside one component): with JOH's top eyes 170 apart (read centre to centre, as for the cello), JOH's feet 155 apart and the Cannone bar position under the bass foot, the 20 mm bar overlaps the bass-side f-hole by about 6 mm.
# The bar position is the only one of the four that is not a published number, so the whole bar line is moved toward the centre (+x) by the smallest
# constant shift that leaves BAR_CLEAR mm between the bar's outer edge and the f-hole's inner edge anywhere along the hole's z span (+-10 mm window).
BAR_CLEAR = 10.0
h1r = fhr[fhr[:, 0] == 1]
gaps = []
for z_, x_ in zip(zb3, xb3):
    near = h1r[np.abs(h1r[:, 1] - z_) <= 10.0]
    if len(near):
        gaps.append((x_ - 0.5 * BAR_W) - float(near[:, 3].max()))
BAR_HOLE_GAP0 = min(gaps)
BAR_HOLE_SHIFT = max(0.0, BAR_CLEAR - BAR_HOLE_GAP0)
xb3 = xb3 + BAR_HOLE_SHIFT
BAR_AT_BRIDGE = float(np.interp(SZ * zbr_c, zb3, xb3))
h3 = hb_c * kh
rows = [(z, BAR_W * h, x, 0.0, x - 0.5 * BAR_W, x + 0.5 * BAR_W, 0.0, h) for z, h, x in zip(zb3, h3, xb3)]
write("bass-bar-sections.tsv", "frankensim.specimen-sections.v1", ["z_mm", "area_mm2", "centroid_x_mm", "centroid_y_mm", "x_min_mm", "x_max_mm", "y_min_mm", "y_max_mm"], rows, ["%.3f"] * 8)
keep = h3 > 2.0
prov("bass_bar", "length %.1f mm of sections (z %.1f to %.1f), %.1f mm above 2 mm height; max height %.1f; wall %.1f mm; centre line x %.1f to %.1f" % (zb3.max() - zb3.min(), zb3.min(), zb3.max(), zb3[keep].max() - zb3[keep].min(), h3.max(), BAR_W, xb3.min(), xb3.max()),
     "JOH / CAN", "JOH bar length approx. 850 (the Cannone's whole section range stretched to it about its mapped middle), thickness approx. 20, height 35 by MB 16036 (= the cello's 24 x 850/580; JOH's printed 'height middle approx. 450' is implausible as mm and is not used; the Cannone height profile x %.3f); ASSUMPTION: Cannone bar profile and line (centre line x = the plan width scale at each z, then shifted %.1f mm in x so that at the bridge plane it sits under the bass foot at the Cannone's position across the foot, x %.1f mm). FURTHER SHIFT, a conflict between published numbers found by the generator's refusal: with JOH's top eyes 170 apart (centre to centre, the cello script's reading), feet 155 apart and the Cannone bar position the 20 mm bar overlapped the bass-side f-hole (gap %.1f mm, negative = overlap); the whole bar line is moved %.1f mm toward the centre so that its outer edge clears the hole's inner edge by %.1f mm; the bar centre at the bridge plane is then x %.1f mm against the bass foot at x -77.5 (the bar is no longer under the foot). The bar position is the only one of these not published." % (kh, BAR_SHIFT, bar_x_target, BAR_HOLE_GAP0, BAR_HOLE_SHIFT, BAR_CLEAR, BAR_AT_BRIDGE))

# ------------------------------------------------------------------ soundpost
P2x, P2z = [float(v) for v in fwd(PX, PZ)]
y_top_post = float(-K_T * (E_T - FT.bil(FT.ymid, PX, PZ)))
y_back_post = float(E_BP + K_B * (FB.bil(FB.ymid, PX, PZ) - E_B))
ys = np.linspace(y_top_post + 6, y_back_post - 6, 28)
rows = [(y, math.pi * (0.5 * POST_D) ** 2, P2x, P2z, P2x - 0.5 * POST_D, P2x + 0.5 * POST_D, P2z - 0.5 * POST_D, P2z + 0.5 * POST_D) for y in ys]
write("soundpost-sections.tsv", "frankensim.specimen-sections.v1", ["y_mm", "area_mm2", "centroid_x_mm", "centroid_z_mm", "x_min_mm", "x_max_mm", "z_min_mm", "z_max_mm"], rows, ["%.3f"] * 8)
prov("soundpost", "diameter %.1f mm, vertical axis at (%.2f, %.2f)" % (POST_D, P2x, P2z), "JOH / CAN", "JOH soundpost thickness 18 (read as a diameter); ASSUMPTION: position = the Cannone post centroid (16.6, -114.3) through the plan map, axis vertical (the Cannone's tilt 0.034 / 0.008 not carried)")

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
# feet bottoms onto the bass top's outer skin: mean of the outer skin under the two feet at the bridge plane
xs_feet = np.linspace((bottom[:, 1].min() - xfeet) * kx_b, (bottom[:, 2].max() - xfeet) * kx_b, 41)
xs_c, zs_c = inv(xs_feet, np.full_like(xs_feet, ZC2))
outer_feet = -K_T * (E_T - FT.bil(FT.ymid, xs_c, zs_c)) - 0.5 * t_top(xs_feet, np.full_like(xs_feet, ZC2), dist_to_poly(cpoly_t, np.c_[xs_feet, np.full_like(xs_feet, ZC2)]))
YBOT2 = float(np.max(outer_feet))  # the lowest point (largest y) the feet reach
# Each Cannone foot is fitted to the Cannone arch (bass foot bottom about 1 mm lower than the treble foot); scaled 2.65x in
# height that offset would be about 3 mm and the generator (feet = runs in the bottom 3 mm of rows) would see one foot.
# So each foot's bottom is put on the bass's outer skin under that foot's centre: the inverse map is
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
     "JOH / MB 16036 / CAN", "JOH bridge feet width 155, height approx. 157, thickness 20; thickness = the Cannone's uniform shell value (volume 3170 mm3 / silhouette %.0f mm2 = %.2f mm) x 20 / the Cannone's foot thickness %.2f mm (z extent of its sections 1-3 mm above the feet); ASSUMPTIONS: the Cannone silhouette's shape (a violin bridge scaled, not a bass bridge), JOH 'thickness 20' read as the foot thickness, the Cannone's mean-to-foot thickness ratio. MASS: NO PUBLISHED BASS BRIDGE MASS FOUND (searched: Brown 2004 text, Askenfelt 1982 text, web 11-03). Askenfelt 1982 Fig. 2 caption gives the magnet (2 g) + accelerometer (0.65 g) mass loading as 'about 3 %% of the weight of the bridge for the bass' (and 30 %% for the violin): 2.65 g / 0.03 = 88 g, but the same caption's violin value gives 8.8 g against about 2 g for the Cannone bridge, so that ratio is not used as an absolute mass; forum posts weigh bridge BLANKS at 116-128 g (unverified, not a finished bridge): not used. So the UNIFORM SHELL stays (MB 16036), labelled: a violin-form bridge scaled to 155 x 157 mm with a uniform shell thickness, not a measured bridge; position = the Cannone bridge plane through the length scale (stop %.1f mm from the top's neck-end edge; JOH stoplength 6250 not used, MB 16036)" % (area_c, t_unif_c, t_foot_c, top_rows[:, 0].max() - ZC2))

# ------------------------------------------------------------------ strings at the bridge
xs_str = np.array([-1.5, -0.5, 0.5, 1.5]) * STR_SPAN / 3.0
rows_s = []
for name, x in zip("cgda", xs_str):
    occ = rb[(rb[:, 1] <= x) & (rb[:, 2] >= x)]
    y_top = float(occ[:, 0].min()) - 0.125
    rows_s.append((name, "bridge", ZC2, x, y_top, 0.0, 0.0))
write("string-contacts.tsv", "frankensim.specimen-string-contacts.v1", ["string", "contact", "z_mm", "x_mm", "y_mm", "slope_dy_dz_before", "slope_dy_dz_after"], rows_s,
      ["%s", "%s", "%.2f", "%.3f", "%.3f", "%.1f", "%.1f"])
prov("strings_at_bridge", " ".join("%s (%.2f, %.2f)" % (r[0], r[3], r[4]) for r in rows_s), "JOH / CAN", "JOH string spacing at bridge 85 (4-string; equal gaps of 28.33 mm, E on the bass-bar side, -x; names c g d a below = E A D G); y = the top of the scaled silhouette at that x (notch depth not modelled); slopes not used by the generator")

# ------------------------------------------------------------------ dimensions the generator reads
write("dimensions.tsv", "frankensim.specimen-dimensions.v1", ["part", "feature", "quantity", "value_si", "unit", "sigma_si", "route", "source_id", "note"],
      [("soundpost", "section", "equivalent-diameter-mean", "%.6f" % (POST_D / 1000), "m", "unstated", "published-maker-table", "johannsson-bass", "soundpost thickness 18 mm (hansjohannsson.com/resources/mesbass.html, column BASS mm.)")],
      ["%s"] * 9)
with open(OUT + "PROVENANCE.tsv", "x", encoding="utf-8") as f:
    f.write("item\tvalue\tsource\tnote\n")
    for r in PROV:
        f.write("\t".join(r) + "\n")
for r in PROV:
    print(r[0], "=", r[1])
