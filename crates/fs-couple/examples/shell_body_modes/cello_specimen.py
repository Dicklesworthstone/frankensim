"""A CELLO specimen directory for gen_violin_shell_assembly.py (instrument=cello), made by mapping the il Cannone
specimen tables (ESRF CT, the h6 violin's inputs) to cello dimensions from published sources. Every mapped quantity is
written with its source or labelled ASSUMPTION in PROVENANCE.tsv beside the tables. Nothing here is tuned to any
recording or mode frequency. cello_bridge_19g.py then sets the bridge to its published mass; the two together write
data/instruments/geometry/cello-castelbarco-1697-shell (check_cello_shell.py reproduces it byte for byte).

usage: python3 -I cello_specimen.py CANNONE_DIR OUT_DIR [instrument=viola] [fholes=RASTER_RUNS.tsv] [bar_slide_mm=MM]   (OUT_DIR must not exist)
       python3 -I cello_specimen.py CANNONE_DIR OUT_DIR fholes=RASTER_RUNS.tsv outline=MEASURED_OUTLINE.tsv fhstop=print|photo [fhole_walls=WALLS.tsv fhole_conductance=TABLE.tsv]

fhole_walls=WALLS.tsv fhole_conductance=TABLE.tsv (optional, both together, outline= only; without them every output byte is as before): the f-holes'
neck for gen_violin_shell_assembly.py neck=bem. WALLS.tsv (frankensim.specimen-fhole-walls.v1, written by cello_fhole_walls.py from this specimen) holds the
DECLARED bores (no CT of this cello exists: each hole's bore is the mean normal thickness of top-plate-maps.tsv along its traced ring, the pair's by the
parallel straight-bore rule); TABLE.tsv (frankensim.specimen-fhole-conductance.v1) is crates/fs-bem/examples/fhole_conductance.rs's output on this
specimen's raster and those walls. Both are checked against this specimen's raster (cells, areas, bores) and copied in as top-plate-fhole-walls.tsv and
top-plate-fhole-conductance.tsv; dimensions.tsv gains the three bore rows the generator reads, and PROVENANCE.tsv the rows f_hole_bores and
f_hole_conductance.

outline=MEASURED_OUTLINE.tsv (optional, needs fholes= and fhstop=, refuses bar_slide_mm=; without it every output byte is as before): the plan outlines are the
Castelbarco's own, measured on the LOC photographs (outline-measured-loc-photographs.tsv: z above the bottom of the back, left / right x of the back and of the
top, derived numbers only). Every borrowed il Cannone field (rib ring, arching, thickness zones, bar, post, bridge) goes through ONE plan map instead of the
uniform stretch: a piecewise-linear z map through five landmarks (bottom edge, lower ledge, waist, upper ledge, top edge = the printed back length) found by one
rule on the Cannone back and on the measured back, and an x map by the per-row ratio of the measured to the Cannone half-width about least-squares centre lines.
What a printed or photographed distance fixes is placed, not warped: the f-hole trace (rigid, moved in z only by fhstop=), the bridge plane (the notch line),
the bar (its length 580, and its slide: rigid in z until its ends sit equidistant from the top plate's own ends, the rule of the posted pages' bar_slide_mm=12.27,
so bar_slide_mm= is refused with outline=). fhstop=print places the trace as stamped (bass notch 419 mm below the top edge); fhstop=photo places the notch line
at the height above the front's measured lower ledge that LOC photograph 0002 shows (the printed stop and the photographs differ by about 14 mm; the rows state both).

bar_slide_mm=MM (optional, cello only; without it every output byte is as before): after the bar is placed, slide it MM mm toward the neck
along z as one piece (length, height profile, width and tilt kept), then re-set its x so that its x at the bridge plane is unchanged, and
add the bass_bar_slide row to PROVENANCE.tsv. The live cello pages carry bar_slide_mm=12.27: the mapped bar's ends sit 71.2 / 95.7 mm from
the top's ends, and Goldblatt's cello specification sheet prints 80 to 85 mm (Green 80), so 12.27 puts them at 83.5 / 83.4 mm.

fholes=RASTER_RUNS.tsv (optional): the two hole rasters traced from the LOCM drawing itself (frankensim.specimen-raster-runs.v1,
already in this specimen's plan mm), used in place of the scaled Cannone holes; the top's island follows them, and the
f_hole_* rows of PROVENANCE.tsv then state the trace and its scored checks (TRACE_* below) instead of the placement.
The trace of image 0021 is data/instruments/geometry/cello-castelbarco-1697-shell/fholes-trace-image0021.tsv (kept as data;
the tools that traced it are not in the tree). No table of that directory is made with it.

With instrument=viola (off by default; without it every output byte is the cello's) the same mapping makes a VIOLA
specimen from the viola's printed values instead, with the same three kinds of source:
  LOC   = Library of Congress catalog record 2022560103, Stradivari viola "Cassavetti" 1727, via the same repo table:
          length of back 410.6, bouts 183 / 123 / 237.
  JOH   = Hans Johannsson, viola specification table, hansjohannsson.com/resources/mesviola.html, Medium column
          (copy md5 1276ed9e3df166102f6d63c558b1ebef; the column whose body length range 405-415 holds LOC 410.6):
          rib height average 38, rib thickness 1.2, arching back 16 / front 17; back thickness: centre at sound post 4.5,
          flanks 2.7, fluting 2.3, centre of bouts 2.6, edge 4; front thickness: centre at sound post 3.2, flanks 2.6,
          fluting 2.4, centre of bouts 2.6, edge 4; overhang 2.3; purfling distance from edge 4; lining width 8;
          soundhole length 78, between top holes 48; stoplength 225; bass bar thickness approx. 6, height middle approx. 14
          (its length cell is blank in this column); soundpost thickness 6.5; bridge feet width 46, thickness 5,
          height approx. 39; string spacing at bridge 36.
  LOCM  = the measured f-hole drawing of the same "Cassavetti" viola in that record's images (resource
          music.musihas-200154830, image 0015: "Viola the Cassevetti Antonius Stradivarius 1727", John Montgomery,
          February 2006, "all measurements in millimeters"; copy md5 805760becef25579715be299729849b1): f-hole lengths
          79.0 (bass) / 78.3 (treble); upper eyes 49.2 apart (inner edges), diameters 8.0 / 8.0; lower eyes 123.5 apart
          (inner edges) and 142.0 (outer edges), diameters 9.7 / 9.25; lower eye to the outside of the instrument
          17.7 / 15.7; stop to the notches 223.0 / 219.0 ("Treble ff hole higher than bass ff"); notches 85.7 apart.
The viola has no printed bridge mass and no printed bass bar length: the bridge keeps the Cannone-ratio volume from the
JOH sizes, and the bar keeps the Cannone's length through the plan length scale; both are ASSUMPTION rows.

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
TRACE, VIOLA, BAR_SLIDE, OUTLINE, FHSTOP, GIB = None, False, None, None, None, None
FH_WALLS, FH_COND = None, None
for arg in sys.argv[3:]:
    key, eq, value = arg.partition("=")
    if key == "fholes" and eq and value and TRACE is None:
        TRACE = value
    elif key == "fhole_walls" and eq and value and FH_WALLS is None:
        FH_WALLS = value
    elif key == "fhole_conductance" and eq and value and FH_COND is None:
        FH_COND = value
    elif key == "outline" and eq and value and OUTLINE is None:
        OUTLINE = value
    elif key == "fhstop" and value in ("print", "photo") and FHSTOP is None:
        FHSTOP = value
    elif key == "bar_slide_mm" and eq and value and BAR_SLIDE is None:
        try:
            BAR_SLIDE = float(value)
        except ValueError:
            raise SystemExit("refusing: bar_slide_mm=%r is not a number" % value)
        if not (0.0 < BAR_SLIDE < 50.0):
            raise SystemExit("refusing: bar_slide_mm=%r: the slide is toward the neck, more than 0 and less than 50 mm" % value)
    elif arg == "instrument=viola" and not VIOLA:
        VIOLA = True
    elif key == "gibson" and eq and value and GIB is None:
        GIB = value
    else:
        raise SystemExit("refusing: %r; the options are one fholes=RASTER_RUNS.tsv, one bar_slide_mm=MM, one outline=MEASURED_OUTLINE.tsv, one fhstop=print|photo, (outline= only) one fhole_walls=WALLS.tsv with one fhole_conductance=TABLE.tsv, instrument=viola and (viola only) gibson=POINTS.tsv" % arg)
if (FH_WALLS is None) != (FH_COND is None):
    raise SystemExit("refusing: fhole_walls= and fhole_conductance= go together (the conductance table was computed through those walls)")
if FH_WALLS and not OUTLINE:
    raise SystemExit("refusing: fhole_walls= and fhole_conductance= are outline= options (the walls and the conductance are the measured-outline specimen's)")
if GIB and not VIOLA:
    raise SystemExit("refusing: gibson=POINTS.tsv is the viola's (Stradivari 'Gibson' viola, Hargrave 1986)")
if TRACE and VIOLA:
    raise SystemExit("refusing: the fholes trace is the Castelbarco cello's (image 0021); the viola places its holes from image 0015")
if BAR_SLIDE is not None and VIOLA:
    raise SystemExit("refusing: bar_slide_mm= is the cello's (Goldblatt's cello specification sheet); the viola keeps its mapped bar")
if OUTLINE and (VIOLA or not TRACE):
    raise SystemExit("refusing: outline= is the Castelbarco cello's measured outline (LOC photographs) and goes with its f-hole trace: give fholes= and not instrument=viola")
if OUTLINE and FHSTOP is None:
    raise SystemExit("refusing: outline= needs fhstop=print (notch line at the printed stop below the top edge) or fhstop=photo (notch line at the photographs' height above the lower ledge)")
if FHSTOP and not OUTLINE:
    raise SystemExit("refusing: fhstop= is an outline= option")
if OUTLINE and BAR_SLIDE is not None:
    raise SystemExit("refusing: outline= places the bar's slide itself (equidistant ends from the top plate's ends); bar_slide_mm= would stack a second slide")
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
if VIOLA:  # the viola's printed values (docstring); the f-holes again come from LOCM, not from JOH 78 / 48 / 13-15
    L_BACK, B_UP, B_C, B_LO = 410.6, 183.0, 123.0, 237.0  # LOC
    RIB_H, RIB_T = 38.0, 1.2  # JOH Medium
    ARCH_TOP, ARCH_BACK = 17.0, 16.0  # JOH
    T_TOP = dict(post=3.2, flank=2.6, flute=2.4, bout=2.6, edge=4.0)  # JOH front
    T_BACK = dict(post=4.5, flank=2.7, flute=2.3, bout=2.6, edge=4.0)  # JOH back
    PURF, LINING_W = 4.0, 8.0  # JOH purfling distance from edge 4, lining width 8
    BAR_L, BAR_W, BAR_H = None, 6.0, 14.0  # JOH thickness and height; no printed length (None: the Cannone bar through the plan length scale)
    POST_D = 6.5  # JOH
    BR_FEET, BR_T_FOOT, BR_H = 46.0, 5.0, 39.0  # JOH
    STR_SPAN = 36.0  # JOH
GIB_PTS = None
if GIB:  # measured Stradivari viola (the 'Gibson', Hargrave, The Strad 1986-09 pp. 113-117) in place of Johannsson's values for the rib and the top interior; the back edge only
    RIB_H, RIB_T = 36.3125, 1.0  # Hargrave p. 116 rib heights 35.5 35.5 35.75 36.75 36.5 37 36.5 36.5 (mean 36.3125); p. 115 rib thickness about 1 mm
    T_BACK = dict(T_BACK, edge=4.75)  # Hargrave p. 116 edge thickness 'taken from the back only': upper and lower bouts 4.75
    _gl = [l.rstrip("\n").split("\t") for l in open(GIB, encoding="utf-8") if l.strip() and not l.startswith("#")]
    assert _gl[0][:3] == ["thickness_mm", "u", "v"], _gl[0]
    GIB_PTS = np.array([[float(x) for x in r[:3]] for r in _gl[1:]])
    assert GIB_PTS.shape == (47, 3), GIB_PTS.shape
NAME = "viola" if VIOLA else "cello"

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
LOC_TXT = dict(length="%g" % L_BACK, bouts="%g, %g, %g" % (B_LO, B_C, B_UP))  # as printed: 749 and 440, 234, 346 for the cello
if not OUTLINE:
    prov("plan_length_scale", "%.6f" % SZ, "LOC / CAN", "z' = %.4f z: LOC length of back %s mm over the Cannone body length %.1f mm (top outline rows at least 10 mm wide, i.e. without the back's button)" % (SZ, LOC_TXT["length"], L_CAN))
    prov("plan_width_scale", "lower %.4f at z %.2f, C %.4f at z %.2f, upper %.4f at z %.2f (Cannone z)" % (ST_S[0], ST_Z[0], ST_S[1], ST_Z[1], ST_S[2], ST_Z[2]), "LOC / CAN",
         "x' = s(z) x, s piecewise linear in the Cannone z through the three bout stations of the Cannone back (widths %.1f, %.1f, %.1f mm -> LOC %s), constant beyond; ASSUMPTION: the %s outline is the Cannone outline so stretched (no %s outline curve is published in reach)" % (wb[kL], wb[kC], wb[kU], LOC_TXT["bouts"], NAME, NAME))


def sx(zc):  # Cannone z -> width scale
    return np.interp(zc, ST_Z, ST_S)


def fwd(x, z):  # Cannone plan -> cello plan
    x, z = np.asarray(x, float), np.asarray(z, float)
    return sx(z) * x, SZ * z


def inv(xp, zp):  # cello plan -> Cannone plan
    z = np.asarray(zp, float) / SZ
    return np.asarray(xp, float) / sx(z), z


def ZF(zc):  # the plan map's z part alone: Cannone z -> cello plan z (SZ z without outline=)
    return SZ * zc


# ------------------------------------------------------------------ outline=: the MEASURED Castelbarco outline, and ONE landmark map for every borrowed violin field
# The measured outline is data/instruments/geometry/cello-castelbarco-1697-shell/outline-measured-loc-photographs.tsv (derived numbers only: z above the bottom of the
# back, left / right x per plate, from the LOC photographs p1 (back) and p2 (front) with the scale of the printed bouts 346 / 234 / 440 mm). With it the plan map is
# ONE piecewise-linear z map through five landmarks (bottom, lower ledge, waist, upper ledge, top: the same five found on the il Cannone back and on the measured back)
# and, for x, the per-row ratio of the measured to the Cannone half-width about each outline's own centre line (3-row mean). Outline, ribs, arching, thickness
# zones, bar, post and bridge all go through this one map (no blend); the f-holes, the bridge plane and the bar length are re-placed from printed distances below.
if OUTLINE:
    def _smooth(a, k):  # moving mean over k rows (odd), edge-padded
        a = np.asarray(a, float)
        if k <= 1:
            return a
        return np.convolve(np.pad(a, k // 2, mode="edge"), np.ones(k) / k, mode="valid")

    def _landmarks(z, hw, k):
        """(lower ledge, waist, upper ledge) z of one half-width profile: waist = minimum in 35-65 % of the span; a ledge = the fastest narrowing going up
        (smallest d hw / dz of the k-row mean) between 30 % of the span and the waist (lower), the waist and 80 % of the span (upper)."""
        hs = _smooth(hw, k)
        span = z[-1] - z[0]
        d = np.gradient(hs, z)
        iw = int(np.argmin(np.where((z > z[0] + 0.35 * span) & (z < z[0] + 0.65 * span), hs, np.inf)))
        sl = (z > z[0] + 0.30 * span) & (z < z[iw])
        su = (z > z[iw]) & (z < z[0] + 0.80 * span)
        return float(z[sl][np.argmin(d[sl])]), float(z[iw]), float(z[su][np.argmin(d[su])])

    def _waist_band(z, hw, k):
        """(first, last) z of the contiguous band within 0.5 mm of the half-width minimum: the minimum is soft, and the waist landmark is the row of the minimum itself.
        The smoothed profile is rounded to 1e-6 mm first: a row exactly 0.5 mm above the minimum is a tie that the last bit of the arithmetic decides, and the last bit differs between architectures."""
        hs = np.round(_smooth(hw, k), 6)
        span = z[-1] - z[0]
        iw = int(np.argmin(np.where((z > z[0] + 0.35 * span) & (z < z[0] + 0.65 * span), hs, np.inf)))
        lo_i, hi_i = iw, iw
        while lo_i > 0 and hs[lo_i - 1] <= hs[iw] + 0.5:
            lo_i -= 1
        while hi_i < len(z) - 1 and hs[hi_i + 1] <= hs[iw] + 0.5:
            hi_i += 1
        return float(z[lo_i]), float(z[hi_i])

    _ol = open(OUTLINE, encoding="utf-8").read().splitlines()
    if _ol[0] != "frankensim.measured-outline-rows.v1" or _ol[1].split("\t") != ["z_mm_above_bottom", "back_x_left_mm", "back_x_right_mm", "front_x_left_mm", "front_x_right_mm"]:
        raise SystemExit("refusing: %s is not the measured-outline table" % OUTLINE)
    _M = np.array([[float(v) if v != "" else np.nan for v in l.split("\t")] for l in _ol[2:] if l.strip()])

    def _meas(c0):  # (z, left, right) of the rows where the plate was measured
        ok = ~np.isnan(_M[:, c0]) & ~np.isnan(_M[:, c0 + 1])
        return _M[ok, 0], _M[ok, c0], _M[ok, c0 + 1]

    MB_Z, MB_L, MB_R = _meas(1)  # back (photograph p1)
    MT_Z, MT_L, MT_R = _meas(3)  # top / front (photograph p2)
    # frame: bottom edge of the plate at the old frame's, the top edge at the old frame's too (so the length is still the printed 749 and Z_TOP the trace's 191.47)
    L_TOP_C = float(top_ol[wt >= 10, 0].max()) + 0.5  # Cannone z of the top edge (half a 1 mm row beyond the last row at least 10 mm wide)
    Z_TOP_OLD = SZ * L_TOP_C
    Z_BOT_EDGE = Z_TOP_OLD - L_BACK
    # landmarks: the Cannone back rows between the plate's bottom and top edges (the button is above), the measured back rows
    _kb = (back_ol[:, 0] < L_TOP_C) & ((back_ol[:, 2] - back_ol[:, 1]) >= 10.0)
    cz, chw = back_ol[_kb, 0], 0.5 * (back_ol[_kb, 2] - back_ol[_kb, 1])
    C_LM = _landmarks(cz, chw, 3)
    M_LM = _landmarks(MB_Z, 0.5 * (MB_R - MB_L), 5)
    MF_LM = _landmarks(MT_Z, 0.5 * (MT_R - MT_L), 5)  # the front's own landmarks: the f-holes are cut in the top, so the notch height reads against the FRONT's lower ledge
    _spread = lambda z_, hw_, ks_=(3, 5, 7, 9): (np.array([_landmarks(z_, hw_, k_) for k_ in ks_]).min(0), np.array([_landmarks(z_, hw_, k_) for k_ in ks_]).max(0))
    SP_B, SP_F, SP_C = _spread(MB_Z, 0.5 * (MB_R - MB_L)), _spread(MT_Z, 0.5 * (MT_R - MT_L)), _spread(cz, chw, (1, 3, 5, 7))
    WB_B, WB_F, WB_C = _waist_band(MB_Z, 0.5 * (MB_R - MB_L), 5), _waist_band(MT_Z, 0.5 * (MT_R - MT_L), 5), _waist_band(cz, chw, 3)
    LM_C = np.array([L_TOP_C - L_CAN, *C_LM, L_TOP_C])
    LM_M = Z_BOT_EDGE + np.array([0.0, *M_LM, L_BACK])
    _s0, _s1 = (LM_M[1] - LM_M[0]) / (LM_C[1] - LM_C[0]), (LM_M[4] - LM_M[3]) / (LM_C[4] - LM_C[3])
    ZC_N = np.r_[LM_C[0] - 100.0, LM_C, LM_C[-1] + 100.0]
    ZM_N = np.r_[LM_M[0] - 100.0 * _s0, LM_M, LM_M[-1] + 100.0 * _s1]

    def ZF(zc):  # noqa: F811 - Cannone z -> cello plan z
        return np.interp(zc, ZC_N, ZM_N)

    def ZINV(zp):  # cello plan z -> Cannone z
        return np.interp(zp, ZM_N, ZC_N)

    def _ratio_rows(zc_rows, hwc_rows, mz, ml, mr):
        """Per-row ratio (measured / Cannone half-width) at the Cannone rows, held at the ends of the measured coverage, clipped to 0.4-3, 3-row mean."""
        zp = ZF(zc_rows) - Z_BOT_EDGE
        inside = (zp >= mz[0]) & (zp <= mz[-1]) & (hwc_rows >= 5.0)
        raw = np.interp(zp, mz, 0.5 * (mr - ml)) / np.maximum(hwc_rows, 5.0)
        w = np.interp(zc_rows, zc_rows[inside], raw[inside])
        return _smooth(np.clip(w, 0.4, 3.0), 3)

    # map A (every borrowed field): the BACK outline's rows; the plates' own outlines below use their own measured rows
    _zb, _lb, _rb = back_ol[:, 0], back_ol[:, 1], back_ol[:, 2]
    W_A = _ratio_rows(_zb, 0.5 * (_rb - _lb), MB_Z, MB_L, MB_R)

    def _centre_line(z, mid):  # least-squares line through the row midpoints: a corner's left and right edges are not at one z, so the midpoint of a corner row is not the centre
        b_, a_ = np.polyfit(z, mid, 1)
        return lambda zz: a_ + b_ * np.asarray(zz, float)

    CC_LIN = _centre_line(cz, 0.5 * (back_ol[_kb, 1] + back_ol[_kb, 2]))  # the il Cannone back's centre line (rows between its plate edges, at least 10 mm wide)
    CM_B = _centre_line(MB_Z, 0.5 * (MB_L + MB_R))  # the measured centre lines (about 0 by the photographs' own fit)
    CM_T = _centre_line(MT_Z, 0.5 * (MT_L + MT_R))

    def _cm_A(zp):
        return CM_B(zp - Z_BOT_EDGE)

    def fwd(x, z):  # noqa: F811 - Cannone plan -> cello plan (one map for every borrowed field)
        x, z = np.asarray(x, float), np.asarray(z, float)
        zp = ZF(z)
        return _cm_A(zp) + np.interp(z, _zb, W_A) * (x - CC_LIN(z)), zp

    def inv(xp, zp):  # noqa: F811 - cello plan -> Cannone plan
        zp = np.asarray(zp, float)
        z = ZINV(zp)
        return CC_LIN(z) + (np.asarray(xp, float) - _cm_A(zp)) / np.interp(z, _zb, W_A), z

    _xt, _zt = fwd(np.array([-30.0, 0.0, 25.0]), np.array([-200.0, -60.0, 40.0]))
    _xi, _zi = inv(_xt, _zt)
    assert abs(_xi - np.array([-30.0, 0.0, 25.0])).max() < 1e-9 and abs(_zi - np.array([-200.0, -60.0, 40.0])).max() < 1e-9, "fwd / inv are not inverse"
    # how well the map reproduces the measured back outline (the check that the fields and the outline are one geometry), z 40-700 above the bottom edge
    _zq = np.arange(40.0, 700.0, 1.0)
    _zcq = ZINV(Z_BOT_EDGE + _zq)
    _hwq = np.interp(_zcq, _zb, W_A) * np.interp(_zcq, _zb, 0.5 * (_rb - _lb))
    _hwm = np.interp(_zq, MB_Z, 0.5 * (MB_R - MB_L))
    MAP_ERR = (float(np.sqrt(np.mean((_hwq - _hwm) ** 2))), float(np.abs(_hwq - _hwm).max()), float(_zq[np.argmax(np.abs(_hwq - _hwm))]))
    prov("plan_map_landmarks", "z above the bottom edge: measured back lower ledge %.1f, waist %.1f, upper ledge %.1f, top %.1f; Cannone back (own z) %.2f, %.2f, %.2f, top edge %.2f; segment scales %s; front lower ledge %.1f, waist %.1f, upper ledge %.1f; spread over the 3-9 row means (il Cannone 1-7): back %s to %s, front %s to %s, il Cannone %s to %s; the waist minimum is soft: the band within 0.5 mm of it spans z %.1f-%.1f on the back, %.1f-%.1f on the front, %.1f-%.1f on the il Cannone back (the landmark is the minimum row itself)" % (
        M_LM[0], M_LM[1], M_LM[2], L_BACK, C_LM[0], C_LM[1], C_LM[2], L_TOP_C, " / ".join("%.3f" % ((LM_M[i + 1] - LM_M[i]) / (LM_C[i + 1] - LM_C[i])) for i in range(4)), MF_LM[0], MF_LM[1], MF_LM[2],
        "/".join("%.0f" % v for v in SP_B[0]), "/".join("%.0f" % v for v in SP_B[1]), "/".join("%.0f" % v for v in SP_F[0]), "/".join("%.0f" % v for v in SP_F[1]), "/".join("%.0f" % v for v in SP_C[0]), "/".join("%.0f" % v for v in SP_C[1]), *WB_B, *WB_F, *WB_C), "LOC photographs / CAN",
         "ONE plan map for every borrowed il Cannone field (outline rows of the rib ring, arching, thickness zones, bar, post, bridge): z' = piecewise linear through five landmarks (bottom edge, lower ledge, waist, upper ledge, top edge = the printed back length %g mm; ledge = the fastest narrowing of the half-width going up, waist = its minimum, found on the 3-row mean of the Cannone back rows and the 5-row mean of the measured back rows, same rule), linear beyond with the end segments' slopes; x' = measured centre line + (measured / Cannone half-width, 3-row mean, held beyond the measured rows, clipped 0.4-3) x (x - the Cannone back's centre line); the centre lines are least-squares lines through the row midpoints (a corner row's midpoint is not a centre: its left and right edges are at different z). No blend with the old uniform stretch. ASSUMPTION: the il Cannone violin's fields are this cello's fields up to this map (arching and thickness maps: the Cannone's, as before)" % L_BACK)
    _ktp = (top_ol[:, 2] - top_ol[:, 1]) >= 10.0
    _ztc = top_ol[_ktp, 0]
    _zpt = ZF(_ztc) - Z_BOT_EDGE
    _etp = np.interp(_ztc, _zb, W_A) * 0.5 * (top_ol[_ktp, 2] - top_ol[_ktp, 1]) - np.interp(_zpt, MT_Z, 0.5 * (MT_R - MT_L))
    _bt = (_zpt >= 40) & (_zpt <= 700)
    prov("plan_map_check_top", "map half-width (the back's ratio, as the top's fields use it) against the measured FRONT, z 40-700: rms %.2f mm, max %.2f mm at z %.0f" % (
        float(np.sqrt(np.mean(_etp[_bt] ** 2))), float(np.abs(_etp[_bt]).max()), float(_zpt[_bt][np.argmax(np.abs(_etp[_bt]))])), "check", "computed; not fitted. The top outline itself is the measured front's rows, so this is the size of the mismatch between the top's borrowed fields and its outline")
    prov("plan_map_check", "map half-width against the measured back, z 40-700: rms %.2f mm, max %.2f mm at z %.0f" % MAP_ERR, "check", "computed; not fitted. The plates' outlines themselves are the measured rows (below), so this is the size of the mismatch between the fields' geometry and the outline")


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


def map_outline_measured(a, mz, ml, mr, plate):
    """The plate's outline rows at the Cannone rows' mapped z: the MEASURED rows where the photograph has them, the Cannone rows through the plate's
    own held half-width ratio where it has not (the top above the fingerboard's cut, z 730). The back stops where its measured rows stop."""
    ac = a[(a[:, 2] - a[:, 1]) >= 10.0]  # as map_outline: the generator leaves out rows narrower than 10 mm (decided in Cannone units)
    zp = ZF(ac[:, 0])
    zm_ = zp - Z_BOT_EDGE
    Wt = _ratio_rows(ac[:, 0], 0.5 * (ac[:, 2] - ac[:, 1]), mz, ml, mr)
    cc = CC_LIN(ac[:, 0])
    cm = (CM_B if plate == "back" else CM_T)(zm_)
    inside = (zm_ >= mz[0]) & (zm_ <= mz[-1])
    xl = np.where(inside, np.interp(zm_, mz, ml), cm + Wt * (ac[:, 1] - cc))
    xr = np.where(inside, np.interp(zm_, mz, mr), cm + Wt * (ac[:, 2] - cc))
    keep = inside if plate == "back" else np.ones(len(ac), bool)
    return np.c_[zp[keep], xl[keep], xr[keep]]


if OUTLINE:
    top_rows = map_outline_measured(top_ol, MT_Z, MT_L, MT_R, "top")
    back_rows = map_outline_measured(back_ol, MB_Z, MB_L, MB_R, "back")
else:
    top_rows, back_rows = map_outline(top_ol), map_outline(back_ol)
for name, rows in (("top", top_rows), ("back", back_rows)):
    write(name + "-plate-outline.tsv", "frankensim.specimen-outline-rows.v1", ["z_mm", "x_left_mm", "x_right_mm"], rows, ["%.3f"] * 3)
cpoly_t = np.vstack([np.c_[top_rows[:, 2], top_rows[:, 0]], np.c_[top_rows[::-1, 1], top_rows[::-1, 0]]])
cpoly_b = np.vstack([np.c_[back_rows[:, 2], back_rows[:, 0]], np.c_[back_rows[::-1, 1], back_rows[::-1, 0]]])
wcb = back_rows[:, 2] - back_rows[:, 1]
Z_SPLIT = -87.0 * SZ if VIOLA else -183  # lower / upper bout split: the Cannone's z -87 mapped (the cello's -183 is that, rounded)
if OUTLINE:
    _wt = top_rows[:, 2] - top_rows[:, 1]
    _zw_ = Z_BOT_EDGE + M_LM[1]
    _sel_c = (back_rows[:, 0] > Z_BOT_EDGE + M_LM[0]) & (back_rows[:, 0] < Z_BOT_EDGE + M_LM[2])
    prov("outline_measured", "top %d rows (measured z %.1f-%.1f above the bottom edge, then the Cannone top rows through the held ratio to the top edge), back %d rows (measured z %.1f-%.1f, none beyond); frame: bottom edge z %.2f, top edge z %.2f (printed back length %g)" % (
        len(top_rows), MT_Z[0], MT_Z[-1], len(back_rows), MB_Z[0], MB_Z[-1], Z_BOT_EDGE, Z_TOP_OLD, L_BACK), "LOC photographs",
         "the outline rows are the measured ones (LOC item ihas.200154818 photographs p1 back and p2 front, 3.7056 px/mm from the printed bouts 346 / 234 / 440 mm, each edge by threshold, centre line fitted; derived numbers in outline-measured-loc-photographs.tsv): the plate's own rows at the Cannone rows' mapped z, 1 mm table interpolated; the photograph's front stops at z 730 where the fingerboard hides the plate edge, the back at z 748 where the neck heel hides it (the back therefore loses the button rows that the stretched il Cannone back carried: rows at least 10 mm wide up to z 770 above its bottom). ASSUMPTIONS: the photographs are plan views at one scale (the three bout scales differ by 0.7 %); the printed length 749 puts the top edge")
    prov("outline_check", "back widths: lower bout %.1f / C %.1f / upper bout %.1f mm; top: lower %.1f / C %.1f / upper %.1f mm" % (
        wcb[back_rows[:, 0] < _zw_].max(), wcb[_sel_c].min(), wcb[back_rows[:, 0] > _zw_].max(),
        _wt[top_rows[:, 0] < _zw_].max(), _wt[(top_rows[:, 0] > Z_BOT_EDGE + M_LM[0]) & (top_rows[:, 0] < Z_BOT_EDGE + M_LM[2])].min(), _wt[top_rows[:, 0] > _zw_].max()), "check",
         "against LOC printed 440 / 234 / 346 mm (the photographs' scale is their mean, so the three agree to the 0.7 % spread, not exactly)")
else:
    prov("outline_check", "back length %.1f, bouts lower %.1f / C %.1f / upper %.1f mm" % (back_rows[:, 0].max() - back_rows[:, 0].min() + SZ, wcb[back_rows[:, 0] < Z_SPLIT].max(),
         wcb[(back_rows[:, 0] > ST_Z[0] * SZ) & (back_rows[:, 0] < ST_Z[2] * SZ)].min(), wcb[back_rows[:, 0] > Z_SPLIT].max()), "check", "the mapped back outline against LOC %s / %g / %g / %g (length includes the button rows that are at least 10 mm wide in Cannone units)" % (LOC_TXT["length"], B_LO, B_C, B_UP))


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
FH_PHOTO_ABOVE_LEDGE, FH_PHOTO_STOP = 19.0, 403.6  # outline= only: LOC photograph p2 (3.7056 px/mm from the printed bouts 346 / 234 / 440): notch nicks (inner stem edges) at rows 1589 / 1588, lower ledge row 1659, top edge row 93: mm above the lower ledge, mm below the top edge
LO_EDGE_P = {1: 32.8, 2: 34.0}  # LOCM lower eye to the outside of the instrument, bass / treble (check only)
LOCM_REF = "LOC item ihas.200154818 image 0021 (Montgomery 2005, Castelbarco cello, caliper mm)"
LOCM_TXT = dict(tape=" (103.0 by tape)", stop="'bass ff hole higher than treble'", lo_dia="%.1f / %.1f")
if VIOLA:  # the Cassavetti sheet draws the length from the upper wing's tip to the lower wing's tip, the whole hole as here
    FH_LEN_P = {1: 79.0, 2: 78.3}
    UP_INNER = 49.2
    UP_DIA = {1: 8.0, 2: 8.0}
    LO_INNER, LO_OUTER = 123.5, 142.0
    LO_DIA = {1: 9.7, 2: 9.25}
    STOP_P = {1: 223.0, 2: 219.0}
    NOTCH_GAP_P = 85.7
    LO_EDGE_P = {1: 17.7, 2: 15.7}
    LOCM_REF = "LOC resource music.musihas-200154830 image 0015 (Montgomery 2006, Cassavetti viola, mm)"
    LOCM_TXT = dict(tape="", stop="'Treble ff hole higher than bass ff'", lo_dia="%.1f / %.2f")
_, _, fh = nt("top-plate-fholes.tsv")
_, _, bsec_c = nt("bridge-sections.tsv")
Z_NOTCH_C = float(np.average(bsec_c[:, 3], weights=bsec_c[:, 1]))  # ASSUMPTION: the Cannone's notch line is its bridge plane
Z_TOP = SZ * (float(top_ol[:, 0].max()) + 0.5)  # ASSUMPTION: the top's edge on the centre line is half a 1 mm row beyond its last row
if OUTLINE and abs(Z_TOP - Z_TOP_OLD) > 1e-6:
    raise SystemExit("refusing: the frame's top edge %.4f differs from the uniform stretch's %.4f" % (Z_TOP_OLD, Z_TOP))


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
if TRACE:
    tl = open(TRACE, encoding="utf-8").read().splitlines()
    if tl[:2] != ["frankensim.specimen-raster-runs.v1", "hole\tz_mm\tx_start_mm\tx_end_mm"]:
        raise SystemExit("refusing: %s is not a two-hole raster-runs table" % TRACE)
    fh_rows = [(int(f[0]), float(f[1]), float(f[2]), float(f[3])) for f in (l.split("\t") for l in tl[2:] if l.strip())]
    if sorted({r[0] for r in fh_rows}) != [1, 2]:
        raise SystemExit("refusing: %s must hold hole 1 (bass) and hole 2 (treble)" % TRACE)
    FH_DZ = 0.0  # z shift of the whole trace (rigid: the drawing's own distances are not touched); only with outline=
    if OUTLINE:
        # the trace's notches as traced and placed (bass notch at the printed 419 below the top edge, z 191.47 - 419 = -227.53; treble -224.83)
        _n_raw = 0.5 * (-227.53 + -224.83)
        if FHSTOP == "print":
            FH_DZ = 0.0  # the trace as stamped: the bass notch at the printed stop below the top edge
        else:
            # the photographs: the notch line is FH_PHOTO_ABOVE_LEDGE above the lower ledge (LOC photograph p2: nicks at rows 1589 / 1588, ledge row 1659, 3.7056 px/mm)
            FH_DZ = (Z_BOT_EDGE + MF_LM[0] + FH_PHOTO_ABOVE_LEDGE) - _n_raw
        fh_rows = [(h_, z_ + FH_DZ, xa_, xb_) for h_, z_, xa_, xb_ in fh_rows]
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
prov("f_hole_length", "%.1f / %.1f mm (bass / treble), end to end" % (FH_LEN_P[1], FH_LEN_P[2]), "LOCM",
     "%s; USED: each Cannone hole shape is scaled uniformly by printed / Cannone end-to-end length (%.2f / %.2f mm), scale %.4f / %.4f; placed length %.2f / %.2f mm" % (
         LOCM_REF, H1["len_c"], H2["len_c"], H1["s"], H2["s"], H1["len_new"], H2["len_new"]))
prov("f_hole_upper_eye_spacing", "%.1f mm between the upper eyes' inner edges%s" % (UP_INNER, LOCM_TXT["tape"]), "LOCM",
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
     "%s (%s); USED: each hole is moved along z until the stem's inner edge on the Cannone notch row lies this far below the top edge (z %.2f): notches at z %.2f / %.2f mm; ASSUMPTION: the Cannone notch line = its bridge plane (Cannone z %.2f)" % (
         LOCM_REF, LOCM_TXT["stop"], Z_TOP, H1["n2"][1], H2["n2"][1], Z_NOTCH_C))
prov("f_hole_lower_eye_diameter", (LOCM_TXT["lo_dia"] + " mm (bass / treble)") % (LO_DIA[1], LO_DIA[2]), "LOCM",
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
     "LOCM / CAN", "LOCM (printed on this %s): lengths %.1f / %.1f end to end; upper eyes %.1f apart at their inner edges, diameters %.1f / %.1f; lower eyes %.1f apart at their inner edges and %.1f at their outer edges; stop to the notches %.0f / %.0f from the top edge. Each Cannone hole shape is scaled uniformly to its printed length (Cannone end-to-end %.2f / %.2f mm), its upper eye centre put at x = -/+(%.1f + d_upper)/2, then rotated about that eye until the lower eye centre is at x = -/+(%.1f + %.1f)/4 = -/+%.2f below it, and moved along z until the stem's inner edge on the Cannone notch row lies at the printed stop below the top edge (z %.2f). ASSUMPTIONS: the Cannone hole shape; eye centres = the largest inscribed circles of the Cannone raster in the top and bottom 30 %% of each hole (Cannone radii %.2f, %.2f / %.2f, %.2f mm); hole ends = the centres of the first and last raster rows; the Cannone notch line = its bridge plane (Cannone z %.2f; the treble hole's inner nick is visible on the raster at z -106.75, the bass hole's is not resolved at 0.5 mm); the top's edge on the centre line half a 1 mm row beyond the last outline row; each spacing split symmetrically about the centre line" % (
    NAME, FH_LEN_P[1], FH_LEN_P[2], UP_INNER, UP_DIA[1], UP_DIA[2], LO_INNER, LO_OUTER, STOP_P[1], STOP_P[2],
    H1["len_c"], H2["len_c"], UP_INNER, LO_INNER, LO_OUTER, 0.25 * (LO_INNER + LO_OUTER), Z_TOP, H1["r_up"], H1["r_lo"], H2["r_up"], H2["r_lo"], Z_NOTCH_C))
prov("f_hole_check", ("eye diameters upper %.2f / %.2f (LOCM %.1f / %.1f), lower %.2f / %.2f (LOCM " + LOCM_TXT["lo_dia"] + ") mm; notch inner edges %.2f apart (LOCM %.1f); lower eye edge to the top outline %.2f / %.2f mm (LOCM %.1f / %.1f to the outside of the instrument); bridge plane z %.2f = %.2f mm below the top edge, %.2f / %.2f mm above the printed notch line (the bridge is not moved)") % (
    2 * H1["s"] * H1["r_up"], 2 * H2["s"] * H2["r_up"], UP_DIA[1], UP_DIA[2], 2 * H1["s"] * H1["r_lo"], 2 * H2["s"] * H2["r_lo"], LO_DIA[1], LO_DIA[2],
    H2["n2"][0] - H1["n2"][0], NOTCH_GAP_P, H1["lo_edge"], H2["lo_edge"], LO_EDGE_P[1], LO_EDGE_P[2],
    SZ * Z_NOTCH_C, Z_TOP - SZ * Z_NOTCH_C, STOP_P[1] - (Z_TOP - SZ * Z_NOTCH_C), STOP_P[2] - (Z_TOP - SZ * Z_NOTCH_C)), "check", "computed from the placed holes; not fitted")
if TRACE:
    # The f_hole_* rows state the trace instead (image 0021 trace, run8: classify.py, measure.py v2, not in the tree; the checks were
    # scored against predictions stamped before measuring). Points in this plan's mm: x = (px - 1504.0) / 5.9284,
    # z = -227.53 - (py - 1395) / 5.9284 (the bass notch row at the printed 419 below the top edge, z 191.47).
    if abs(Z_TOP - 191.47) > 0.005:
        raise SystemExit("refusing: the trace was placed with the top edge at z 191.47, this specimen has %.3f" % Z_TOP)
    T_EYE = {(1, "up"): (-58.13, -174.46, 6.61), (2, "up"): (57.54, -173.40, 6.01), (1, "lo"): (-119.70, -270.05, 9.60), (2, "lo"): (119.58, -266.07, 9.21)}
    T_NOTCH = {1: (-79.87, -227.53), 2: (80.54, -224.83)}
    if OUTLINE:  # the whole trace moved rigidly by FH_DZ
        T_EYE = {k_: (v_[0], v_[1] + FH_DZ, v_[2]) for k_, v_ in T_EYE.items()}
        T_NOTCH = {k_: (v_[0], v_[1] + FH_DZ) for k_, v_ in T_NOTCH.items()}
        BR_PLANE = 0.5 * (T_NOTCH[1][1] + T_NOTCH[2][1])  # the bridge stands on the notch line: the mean of the two traced notches' z
    lo_t = {h: float(dist_to_poly(cpoly_t, np.array([T_EYE[(h, "lo")][:2]]))[0]) - T_EYE[(h, "lo")][2] for h in (1, 2)}
    fh_new = [
        ("f_hole_length", "%.1f / %.1f mm (bass / treble), end to end" % (FH_LEN_P[1], FH_LEN_P[2]), "LOCM",
         "%s; CHECK only: the traced holes are 130.78 / 128.24 mm eye to eye over the arch (130.24 / 127.23 projected), so the printed lengths are eye to eye; the traced wing tips reach beyond the eyes" % LOCM_REF),
        ("f_hole_upper_eye_spacing", "%.1f mm between the upper eyes' inner edges (103.0 by tape)" % UP_INNER, "LOCM",
         "%s; USED for the trace scale, one of three chords (611.0 px on the drawing: 5.9723 px/mm)" % LOCM_REF),
        ("f_hole_upper_eye_diameter", "%.1f / %.1f mm (bass / treble)" % (UP_DIA[1], UP_DIA[2]), "LOCM",
         "%s; CHECK only: traced across-body widths over the arch 13.22 / 12.07 mm; treble -0.93 is a KNOWN MISS (T3, see f_hole_check)" % LOCM_REF),
        ("f_hole_lower_eye_spacing_inner", "%.1f mm between the lower eyes' inner edges" % LO_INNER, "LOCM",
         "%s; USED for the trace scale, one of three chords (1307.2 px: 5.9416 px/mm)" % LOCM_REF),
        ("f_hole_lower_eye_spacing_outer", "%.1f mm between the lower eyes' outer edges" % LO_OUTER, "LOCM",
         "%s; CHECK only: 258.12 mm on the trace at the mean scale" % LOCM_REF),
        ("f_hole_stop", "%.0f / %.0f mm from the top edge to the notches (bass / treble)" % (STOP_P[1], STOP_P[2]), "LOCM",
         "%s ('bass ff hole higher than treble'); USED for the bass hole only: the drawing is moved along z until the bass stem's inner notch lies %.1f mm below the top edge (z %.2f), notch z %.2f; the treble notch follows as drawn, %.2f mm below the top edge: KNOWN MISS (T9, see f_hole_check)" % (
             LOCM_REF, STOP_P[1], Z_TOP, T_NOTCH[1][1], Z_TOP - T_NOTCH[2][1])),
        ("f_hole_lower_eye_diameter", "%.1f / %.1f mm (bass / treble)" % (LO_DIA[1], LO_DIA[2]), "LOCM",
         "%s; CHECK only: traced across-body widths over the arch 19.55 / 18.10 mm; treble -1.10 is a KNOWN MISS (T3, see f_hole_check)" % LOCM_REF),
        ("f_hole_notch_gap", "%.1f mm between the notches" % NOTCH_GAP_P, "LOCM",
         "%s; USED for the trace scale, one of three chords (951.1 px: 5.8712 px/mm); at the mean scale the traced notches are 160.44 mm apart (-1.56)" % LOCM_REF),
        ("f_hole_lower_eye_to_edge", "%.1f / %.1f mm from the lower eye to the outside of the instrument (bass / treble)" % (LO_EDGE_P[1], LO_EDGE_P[2]), "LOCM",
         "%s; CHECK only, not used: traced lower eye circle edge to the mapped top outline %.2f / %.2f mm" % (LOCM_REF, lo_t[1], lo_t[2])),
        ("f_hole_shape", "traced from the LOCM drawing (image 0021) at 5.9284 px/mm: hole = the pixels darker than the half-intensity threshold between the hole and wood medians (53.1 between 20.1 and 86.1); the drawn ink lines are masked and each masked segment is decided from its two ends; the four necks are constructed as capsules at the printed widths 1.4 / 1.3 / 1.3 / 1.3 mm between the drawing's tip markers; cut into 0.5 mm raster rows", "LOCM",
         "ASSUMPTIONS: the drawing is a plan view at one scale (its caliper chords run over the arch, which adds 0.54 / 1.01 mm to the traced lengths); the centre line is the mean chord midpoint (x 1504.0 px); no rotation (the drawing's chords tilt treble side up by 0.53 to 0.96 deg); the lit wall band inside the treble eyes' rims (0.5 to 0.85 mm) is wood by the half-intensity rule; the necks are not measured"),
        ("f_holes", "traced: %d / %d raster runs, z extents %.1f / %.1f mm, areas %.0f / %.0f mm2; upper eyes (%.2f, %.2f) / (%.2f, %.2f), lower eyes (%.2f, %.2f) / (%.2f, %.2f), notches (%.2f, %.2f) / (%.2f, %.2f)" % (
            int((fhr[:, 0] == 1).sum()), int((fhr[:, 0] == 2).sum()), H1["len2"], H2["len2"], H1["area2"], H2["area2"], *T_EYE[(1, "up")][:2], *T_EYE[(2, "up")][:2],
            *T_EYE[(1, "lo")][:2], *T_EYE[(2, "lo")][:2], *T_NOTCH[1], *T_NOTCH[2]), "LOCM",
         "scale = mean of the three printed chords (upper eyes inner 102.3, notches 162.0, lower eyes inner 220.0; +0.74 / -0.96 / +0.22 %% from the mean); eye centres = circle fits to the traced eye boundaries (radii %.2f, %.2f / %.2f, %.2f mm)" % (
             T_EYE[(1, "up")][2], T_EYE[(1, "lo")][2], T_EYE[(2, "up")][2], T_EYE[(2, "lo")][2])),
        ("f_hole_check", "scored against predictions stamped before measuring: eye widths over the arch upper 13.22 / 12.07 (LOCM 13.5 / 13.0), lower 19.55 / 18.10 (19.7 / 19.2) mm; lengths eye to eye over the arch 130.78 / 128.24 (130.5 / 128.5); notch gap 160.44 (162.0); treble stop 416.30 (420); bridge plane z %.2f = %.2f mm below the top edge, %.2f / %.2f mm above the traced notches (the bridge is not moved)" % (
            SZ * Z_NOTCH_C, Z_TOP - SZ * Z_NOTCH_C, SZ * Z_NOTCH_C - T_NOTCH[1][1], SZ * Z_NOTCH_C - T_NOTCH[2][1]), "check",
         "KNOWN MISSES, kept as measured (the trace is placed as stamped; nothing is moved after measuring): T9 the treble stop is 3.70 mm short of the printed 420 (tolerance 1.0): the drawing's treble notch lies 2.70 mm above its bass notch and all three chords tilt treble side up (-0.53, -0.96, -0.96 deg); T3 the treble eyes are 0.93 / 1.10 mm narrower than printed (tolerance 0.5): the lit wall band inside their rims is wood by the half-intensity rule. In tolerance: the bass eyes (-0.28, -0.15) and both lengths (+0.28, -0.26)"),
    ]
    if OUTLINE:  # the stop and the check state the placement on the measured outline
        _stop_b, _stop_t = Z_TOP - T_NOTCH[1][1], Z_TOP - T_NOTCH[2][1]
        _led = Z_BOT_EDGE + MF_LM[0]
        _what = ("fhstop=print: the trace as stamped, the bass notch at the printed stop below the top edge" if FHSTOP == "print" else
                 "fhstop=photo: the trace moved up %.2f mm so that the mean notch height is %.1f mm above the front's measured lower ledge (z %.2f), the height read on LOC photograph p2" % (FH_DZ, FH_PHOTO_ABOVE_LEDGE, _led))
        for i_, row_ in enumerate(fh_new):
            if row_[0] == "f_hole_stop":
                fh_new[i_] = ("f_hole_stop", "%.0f / %.0f mm from the top edge to the notches (bass / treble)" % (STOP_P[1], STOP_P[2]), "LOCM",
                              "%s ('bass ff hole higher than treble'); USED on the measured outline as: %s; placed notches at z %.2f / %.2f = %.2f / %.2f mm below the top edge (z %.2f), %.2f / %.2f mm above the lower ledge. DISAGREEMENT, kept: LOC photograph p2 (3.7056 px/mm from the printed bouts) reads the notch nicks at %.1f mm below the top edge (printed %.0f / %.0f) and %.1f mm above the lower ledge; the printed stop and the photographs differ by about %.0f mm" % (
                                  LOCM_REF, _what, T_NOTCH[1][1], T_NOTCH[2][1], _stop_b, _stop_t, Z_TOP, T_NOTCH[1][1] - _led, T_NOTCH[2][1] - _led, FH_PHOTO_STOP, STOP_P[1], STOP_P[2], FH_PHOTO_ABOVE_LEDGE, 0.5 * (STOP_P[1] + STOP_P[2]) - FH_PHOTO_STOP))
            elif row_[0] == "f_hole_check":
                fh_new[i_] = ("f_hole_check", "scored against predictions stamped before measuring, as traced at the stamped placement (the over-the-arch values depend on the arch under the holes, which the move of the trace changes by up to 13.6 mm in z): eye widths over the arch upper 13.22 / 12.07 (LOCM 13.5 / 13.0), lower 19.55 / 18.10 (19.7 / 19.2) mm; lengths eye to eye over the arch 130.78 / 128.24 (130.5 / 128.5); notch gap 160.44 (162.0). On the measured outline: stop %.2f / %.2f mm (printed %.0f / %.0f, photograph %.1f), notch line %.2f mm above the lower ledge (photograph %.1f), lower eye to the top outline %.2f / %.2f mm (LOCM %.1f / %.1f); bridge plane z %.2f = the mean notch z, %.2f mm below the top edge" % (
                    _stop_b, _stop_t, STOP_P[1], STOP_P[2], FH_PHOTO_STOP, BR_PLANE - _led, FH_PHOTO_ABOVE_LEDGE, lo_t[1], lo_t[2], LO_EDGE_P[1], LO_EDGE_P[2], BR_PLANE, Z_TOP - BR_PLANE), "check",
                    "KNOWN MISSES of the trace, kept as measured: T9 the treble notch lies 2.70 mm above the bass notch on the drawing (all three chords tilt treble side up, -0.53, -0.96, -0.96 deg; the photograph reads 0.3 mm); T3 the treble eyes are 0.93 / 1.10 mm narrower than printed. Nothing is moved after measuring")
    at = [i for i, p in enumerate(PROV) if p[0].startswith("f_hole")]
    if at != list(range(at[0], at[0] + len(fh_new))):
        raise SystemExit("refusing: the f_hole rows are not the %d contiguous rows the trace replaces" % len(fh_new))
    PROV[at[0]:at[0] + len(fh_new)] = fh_new

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


def gib_top_interior(xp, zp):
    """The 47 printed belly thicknesses as a map: Delaunay triangulation (triangle, no quality or area switches) of the points' (u, v), barycentric linear inside, the nearest printed point's value outside the hull. (u, v) of a plate point: u = (z_top - z)/(z_top - z_bottom) of the mapped top outline, v = (x - outline centre)/(half the outline width) on its row."""
    import triangle
    zt, zb = top_rows[:, 0].max(), top_rows[:, 0].min()
    order = np.argsort(top_rows[:, 0])
    zs, xl, xr = top_rows[order, 0], top_rows[order, 1], top_rows[order, 2]
    xp, zp = np.asarray(xp, float), np.asarray(zp, float)
    u = (zt - zp) / (zt - zb)
    v = (xp - 0.5 * (np.interp(zp, zs, xl) + np.interp(zp, zs, xr))) / (0.5 * (np.interp(zp, zs, xr) - np.interp(zp, zs, xl)))
    P = GIB_PTS[:, [1, 2]]
    T = triangle.triangulate({"vertices": P}, "")["triangles"]
    out = np.full(len(u), np.nan)
    A = P[T[:, 0]]
    B = P[T[:, 1]] - A
    C = P[T[:, 2]] - A
    den = B[:, 0] * C[:, 1] - B[:, 1] * C[:, 0]
    for k in range(len(u)):
        q = np.array([u[k], v[k]]) - A
        l1 = (q[:, 0] * C[:, 1] - q[:, 1] * C[:, 0]) / den
        l2 = (B[:, 0] * q[:, 1] - B[:, 1] * q[:, 0]) / den
        l0 = 1.0 - l1 - l2
        ok = np.nonzero((l0 >= -1e-12) & (l1 >= -1e-12) & (l2 >= -1e-12))[0]
        if len(ok):
            i = ok[0]
            out[k] = l0[i] * GIB_PTS[T[i, 0], 0] + l1[i] * GIB_PTS[T[i, 1], 0] + l2[i] * GIB_PTS[T[i, 2], 0]
        else:
            out[k] = GIB_PTS[np.argmin(((P - np.array([u[k], v[k]])) ** 2).sum(axis=1)), 0]
    return out


def t_top(xp, zp, d):
    if GIB:
        we = edge_w(d)
        return we * T_TOP["edge"] + (1 - we) * gib_top_interior(xp, zp)
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


THICK_NOTES = (  # the JOH values each note reads, and what is not drawn
    "edge 5.5, flanks = fluting 3.9, centre at sound post 4.7; centre of bouts 4.0 is within 0.1 mm of the flank and below the generator's 0.2 mm top zone step, so not drawn",
    "edge 5, flanks = fluting 3.5, centre at sound post 7; centre of bouts 3.6 within 0.1 mm of the flank (below the 0.4 mm back zone step)")
if VIOLA:
    THICK_NOTES = (
        "edge 4, flanks 2.6, fluting 2.4, centre at sound post 3.2, centre of bouts 2.6 (= the flank); fluting 2.4 is 0.2 mm below the flank, one 0.2 mm top zone step, and is NOT drawn: JOH prints no fluting width, so the edge ramp runs straight to the flank (ASSUMPTION)",
        "edge 4, flanks 2.7, fluting 2.3, centre at sound post 4.5; centre of bouts 2.6 within 0.1 mm of the flank (below the 0.4 mm back zone step); fluting 2.3 is 0.4 mm below the flank, one back zone step, and is NOT drawn: JOH prints no fluting width, so the edge ramp runs straight to the flank (ASSUMPTION)")
prov("thickness_top", "edge %.1f to %.0f mm from the outline, linear to the interior by %.0f mm; interior flank %.1f; island between the f-holes %.1f, blended to the flank over 20 mm" % (T_TOP["edge"], PURF, LINING_W, T_TOP["flank"], T_TOP["post"]),
     "JOH front", "%s; ASSUMPTION: the %g region is the island between the holes (bounded by the holes' inner edges over the z span where both holes are open) because the Cannone top shows no thickening at its post (3.08 mm within 6 mm of the post vs interior median 3.16); edge band width = purfling distance %g and ramp end = lining width %g (JOH) are ASSUMPTIONS of where the edge value applies" % (THICK_NOTES[0], T_TOP["post"], PURF, LINING_W))
prov("thickness_back", "edge %.1f to %.0f mm, linear by %.0f mm; interior = %.1f + %.1f x clip((t_Cannone - %.2f) / (%.2f - %.2f)) at the mapped Cannone point" % (T_BACK["edge"], PURF, LINING_W, T_BACK["flank"], T_BACK["post"] - T_BACK["flank"], TB_FLANK, TB_POST, TB_FLANK),
     "JOH back / CAN", "%s; the Cannone back's measured thickening (5.9 mm at its post, %.2f mm = its interior 10th percentile at the flanks) sets WHERE the %s back is thick; ASSUMPTION: the %s back's thickness pattern is the Cannone's, rescaled to JOH's levels (Sacconi's 7 -> 6 -> 5 -> 4 centre contours in Hargrave 1987 p. 37 have the same central-strip form)" % (THICK_NOTES[1], TB_FLANK, NAME, NAME))

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
     "JOH / CAN", "scales solved so that the outer surface's highest point lies JOH %g (front) / %g (back) mm above the plate's glue plane (mean inner skin at the rib ring); ASSUMPTION: the Cannone's arch shape, and JOH 'arching height' read as that height" % (ARCH_TOP, ARCH_BACK))
prov("rib_height", "mean %.2f mm (min %.2f, max %.2f) between the inner skins at the ring" % (rib_h.mean(), rib_h.min(), rib_h.max()), "JOH", "rib height average %g (the variation along the ring is the Cannone's edge undulation x the arch scales)" % RIB_H)
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
     "generator / CAN", "points closer than 2.0 mm to either outline (the generator's own 1.6 mm minimum + 0.4) moved along the nearer outline's inward normal, then re-read at the table's angles; elsewhere the stretched Cannone overhang is kept (%s)" % (
         "JOH's plate overhang 2.3 mm + half rib 0.6 mm = 2.9 is NOT imposed: the stretched Cannone mean is %s" % ("larger" if oh0.mean() > 2.9 else "smaller") if VIOLA else
         "JOH's plate overhang 3 mm + half rib 0.85 mm = 3.85 is NOT imposed: the stretched Cannone mean is larger"))
rows = [(rs[i, 0], CX2, CZ2, RIB_T * RIB_H, rmin2[i], rmax2[i], float(inner_t_ring.mean()), float(inner_b_ring.mean()), RIB_T) for i in range(len(rs))]
write("ribs-sections.tsv", "frankensim.specimen-ring-sections.v1", ["theta_deg", "centre_x_mm", "centre_z_mm", "area_mm2", "r_min_mm", "r_max_mm", "y_min_mm", "y_max_mm", "mean_width_mm"], rows,
      ["%.0f", "%.3f", "%.3f", "%.2f", "%.3f", "%.3f", "%.2f", "%.2f", "%.2f"])
prov("ribs", "ring = the Cannone ring mapped through the plan scales, radii re-measured about the mapped centre (%.2f, %.2f) at the table's 2 degree angles; wall %.1f mm" % (CX2, CZ2, RIB_T), "JOH / CAN", "rib thickness %g (JOH; written as every section's mean width, so the generator's 10th percentile is %g); y columns and area are not read by the generator" % (RIB_T, RIB_T))

# ------------------------------------------------------------------ bass bar
_, _, bs = nt("bass-bar-sections.tsv")
hb_c = bs[:, 7] - bs[:, 6]
kh = BAR_H / hb_c.max()
xb2, zb2 = fwd(bs[:, 2], bs[:, 0])
zmid = 0.5 * (zb2.min() + zb2.max())
klen = 1.0 if BAR_L is None else BAR_L / (zb2.max() - zb2.min())  # None: no printed length, the plan length scale alone
zb3 = zmid + (zb2 - zmid) * klen
if OUTLINE:
    xb3 = fwd(bs[:, 2], ZINV(zb3))[0]  # centre line follows the one plan map at the stretched z (the bar is then slid rigidly in z below, as on the posted pages)
else:
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
BAR_SHIFT = bar_x_target - float(np.interp(BR_PLANE if OUTLINE else SZ * zbr_c, zb3, xb3))
xb3 = xb3 + BAR_SHIFT
if OUTLINE:
    # the rule of the posted pages (bar_slide_mm=12.27; Goldblatt prints the bar's ends 80-85 mm from the plate's ends, Green 80): the bar slides rigidly in z until its ends sit equidistant
    # from the TOP plate's own two ends (its first and last outline rows); its length, height profile, width and tilt stay; x is re-set so that its x at the bridge plane (the notch line) stays too
    _bar_lo0, _bar_hi0 = float(zb3.min()), float(zb3.max())
    OUTLINE_SLIDE = 0.5 * ((float(top_rows[:, 0].min()) + float(top_rows[:, 0].max())) - (_bar_lo0 + _bar_hi0))
    zb3 = zb3 + OUTLINE_SLIDE
    xb3 = xb3 + (bar_x_target - float(np.interp(BR_PLANE, zb3, xb3)))
    prov("bass_bar_slide", "bar slid %.2f mm along z (toward the %s), rigid; x at the bridge plane unchanged; ends %.2f / %.2f mm from the top plate's ends (the mapped bar sat %.2f / %.2f mm)" % (
        OUTLINE_SLIDE, "neck" if OUTLINE_SLIDE > 0 else "tail", zb3.min() - float(top_rows[:, 0].min()), float(top_rows[:, 0].max()) - zb3.max(),
        _bar_lo0 - float(top_rows[:, 0].min()), float(top_rows[:, 0].max()) - _bar_hi0), "makers' tables",
         "the rule of the posted pages, re-applied on the measured outline: Goldblatt prints the bar's ends 80-85 mm from the plate's ends, Green 80; the ends are measured from the top plate's own first and last outline rows, so equidistant ends are (top length - %.0f) / 2 from each. Makers' tables, borrowed and named, not measurements of this cello" % (zb3.max() - zb3.min()))
elif BAR_SLIDE is not None:  # bar_slide_mm=: rigid slide toward the neck, x re-set so that the bar's x at the bridge plane stays
    zb3 = zb3 + BAR_SLIDE
    xb3 = xb3 + (bar_x_target - float(np.interp(SZ * zbr_c, zb3, xb3)))
    prov("bass_bar_slide", "bar slid %.2f mm toward the neck along z, rigid; x at the bridge plane unchanged" % BAR_SLIDE, "makers' tables",
         "Goldblatt prints the bar's ends 80-85 mm from the plate's ends, Green 80; the mapped bar sat 71.2 / 95.7 mm. Makers' tables, borrowed and named, not measurements of this cello")
h3 = hb_c * kh
rows = [(z, BAR_W * h, x, 0.0, x - 0.5 * BAR_W, x + 0.5 * BAR_W, 0.0, h) for z, h, x in zip(zb3, h3, xb3)]
write("bass-bar-sections.tsv", "frankensim.specimen-sections.v1", ["z_mm", "area_mm2", "centroid_x_mm", "centroid_y_mm", "x_min_mm", "x_max_mm", "y_min_mm", "y_max_mm"], rows, ["%.3f"] * 8)
keep = h3 > 2.0
prov("bass_bar", "length %.1f mm of sections (z %.1f to %.1f), %.1f mm above 2 mm height; max height %.1f; wall %.1f mm; centre line x %.1f to %.1f" % (zb3.max() - zb3.min(), zb3.min(), zb3.max(), zb3[keep].max() - zb3[keep].min(), h3.max(), BAR_W, xb3.min(), xb3.max()),
     "JOH / CAN", "%s, height middle approx. %g (the Cannone height profile x %.3f), thickness approx. %g; ASSUMPTION: Cannone bar profile and line (centre line x = the plan width scale at each z, then shifted %.1f mm in x so that at the bridge plane it sits under the bass foot at the Cannone's position across the foot, x %.1f mm)" % (
         "JOH bar length approx. %g (the Cannone's whole section range stretched to it about its mapped middle)" % BAR_L if BAR_L is not None else
         "JOH prints no bar length in this column (blank; the Small and Large cells read '5 to 4' and '4 to 3'), so ASSUMPTION: the Cannone's whole section range through the plan length scale %.4f alone" % SZ,
         BAR_H, kh, BAR_W, BAR_SHIFT, bar_x_target))

# ------------------------------------------------------------------ soundpost
P2x, P2z = [float(v) for v in fwd(PX, PZ)]
y_top_post = float(-K_T * (E_T - FT.bil(FT.ymid, PX, PZ)))
y_back_post = float(E_BP + K_B * (FB.bil(FB.ymid, PX, PZ) - E_B))
ys = np.linspace(y_top_post + 6, y_back_post - 6, 28)
rows = [(y, math.pi * (0.5 * POST_D) ** 2, P2x, P2z, P2x - 0.5 * POST_D, P2x + 0.5 * POST_D, P2z - 0.5 * POST_D, P2z + 0.5 * POST_D) for y in ys]
write("soundpost-sections.tsv", "frankensim.specimen-sections.v1", ["y_mm", "area_mm2", "centroid_x_mm", "centroid_z_mm", "x_min_mm", "x_max_mm", "z_min_mm", "z_max_mm"], rows, ["%.3f"] * 8)
prov("soundpost", "diameter %.1f mm, vertical axis at (%.2f, %.2f)" % (POST_D, P2x, P2z), "JOH / CAN", "JOH soundpost thickness %g; ASSUMPTION:" % POST_D + " position = the Cannone post centroid (16.6, -114.3) through the plan map, axis vertical (the Cannone's tilt 0.034 / 0.008 not carried)")

# ------------------------------------------------------------------ bridge
_, _, sil = nt("bridge-silhouette.tsv")
_, _, bsec = nt("bridge-sections.tsv")
ybot, ytop = sil[:, 0].max(), sil[:, 0].min()
bottom = sil[sil[:, 0] > ybot - 3.0]
feet_c = bottom[:, 2].max() - bottom[:, 1].min()
kx_b = BR_FEET / feet_c
ky_b = BR_H / (ybot - ytop + 0.25)
zc_c = float(np.average(bsec[:, 3], weights=bsec[:, 1]))
ZC2 = BR_PLANE if OUTLINE else float(SZ * zc_c)
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
STOP_JOH = 225 if VIOLA else 405
prov("bridge", "silhouette x %.4f (feet %.1f -> %g mm), y %.4f (height %.2f -> %g mm), plane z %.2f, feet bottoms on the outer skin at y %.2f (bass) / %.2f (treble), Cannone feet bottoms %.2f / %.2f; uniform shell thickness %.2f mm; area %.0f mm2; volume %.0f mm3" % (kx_b, feet_c, BR_FEET, ky_b, ybot - ytop + 0.25, BR_H, ZC2, skl, skr, fyl, fyr, T_BR, area2, T_BR * area2),
     "JOH / CAN", "JOH bridge feet width %g, height approx. %g, thickness %g; thickness = the Cannone's uniform shell value (volume 3170 mm3 / silhouette %.0f mm2 = %.2f mm) x %g / the Cannone's foot thickness %.2f mm (z extent of its sections 1-3 mm above the feet); ASSUMPTIONS: the Cannone silhouette's shape, JOH 'thickness %g' read as the foot thickness, the Cannone's mean-to-foot thickness ratio; position = the Cannone bridge plane through the length scale (stop %.1f mm from the top's neck-end edge vs JOH stoplength %d)" % (
         BR_FEET, BR_H, BR_T_FOOT, area_c, t_unif_c, BR_T_FOOT, t_foot_c, BR_T_FOOT, top_rows[:, 0].max() - ZC2, STOP_JOH)
     + ("; MASS: no printed viola bridge mass was found, so the bridge keeps this volume: %.2f g at the maple pack's 626 kg/m3 (Sonderegger 2013), an ESTIMATE from the sizes above, not a measured or typical mass" % (T_BR * area2 * 626.0e-6) if VIOLA else ""))

# ------------------------------------------------------------------ strings at the bridge
xs_str = np.array([-1.5, -0.5, 0.5, 1.5]) * STR_SPAN / 3.0
rows_s = []
for name, x in zip("cgda", xs_str):
    occ = rb[(rb[:, 1] <= x) & (rb[:, 2] >= x)]
    y_top = float(occ[:, 0].min()) - 0.125
    rows_s.append((name, "bridge", ZC2, x, y_top, 0.0, 0.0))
write("string-contacts.tsv", "frankensim.specimen-string-contacts.v1", ["string", "contact", "z_mm", "x_mm", "y_mm", "slope_dy_dz_before", "slope_dy_dz_after"], rows_s,
      ["%s", "%s", "%.2f", "%.3f", "%.3f", "%.1f", "%.1f"])
prov("strings_at_bridge", " ".join("%s (%.2f, %.2f)" % (r[0], r[3], r[4]) for r in rows_s), "JOH / CAN", "JOH string spacing at bridge %g (equal gaps, C on the bass-bar side, -x); y = the top of the scaled silhouette at that x (notch depth not modelled); slopes not used by the generator" % STR_SPAN)

# ------------------------------------------------------------------ the f-holes' neck (fhole_walls= fhole_conductance=)
NECK_DIMS, NECK_PROV, NECK_FILES = [], [], []
if FH_WALLS:
    def _rows(path, magic):
        lines_ = open(path, encoding="utf-8").read().splitlines()
        if not lines_ or lines_[0] != magic:
            raise SystemExit("refusing: %s is not a %s table" % (path, magic))
        return lines_[1].split("\t"), [l.split("\t") for l in lines_[2:] if l.strip()]

    _wh, _wr = _rows(FH_WALLS, "frankensim.specimen-fhole-walls.v1")
    _ch, _cr = _rows(FH_COND, "frankensim.specimen-fhole-conductance.v1")
    if _wh != ["hole", "quantity", "value", "unit", "note"]:
        raise SystemExit("refusing: %s: unexpected columns %s" % (FH_WALLS, _wh))
    _w = {(r[0], r[1]): r[2] for r in _wr}
    _c = {r[0]: dict(zip(_ch, r)) for r in _cr}
    _configs = ("pair", "hole1-at-pair-bore", "hole2-at-pair-bore", "hole1-at-own-bore", "hole2-at-own-bore")
    if sorted(_c) != sorted(_configs):
        raise SystemExit("refusing: %s must hold the rows %s" % (FH_COND, ", ".join(_configs)))
    # the raster as written above: cells per hole and the rows' offset from the solver's grid (rows at k/2 + 0.25 mm)
    _runs = [(h_, float("%.2f" % z_), float("%.2f" % xa_), float("%.2f" % xb_)) for h_, z_, xa_, xb_ in fh_rows]
    _cells = {h_: int(round(sum(b_ - a_ for hh, _, a_, b_ in _runs if hh == h_) / 0.5)) for h_ in (1, 2)}
    _off = sorted({round((z_ - 0.25) % 0.5, 6) % 0.5 for _, z_, _, _ in _runs})
    if len(_off) != 1:
        raise SystemExit("refusing: the raster rows are not on one 0.5 mm grid (offsets %s)" % _off)
    for h_ in (1, 2):
        if abs(float(_w[(str(h_), "area-0p5mm")]) - 0.25 * _cells[h_]) > 1e-6:
            raise SystemExit("refusing: %s: hole %d area %s mm2, this raster's %d cells give %.2f" % (FH_WALLS, h_, _w[(str(h_), "area-0p5mm")], _cells[h_], 0.25 * _cells[h_]))
    _want = {"pair": (_cells[1] + _cells[2], ("0", "bore-pair")), "hole1-at-pair-bore": (_cells[1], ("0", "bore-pair")), "hole2-at-pair-bore": (_cells[2], ("0", "bore-pair")),
             "hole1-at-own-bore": (_cells[1], ("1", "bore")), "hole2-at-own-bore": (_cells[2], ("2", "bore"))}
    for name_, (cells_, bore_) in _want.items():
        if int(_c[name_]["cells"]) != cells_ or abs(float(_c[name_]["thickness_mm"]) - float(_w[bore_])) > 1e-9:
            raise SystemExit("refusing: %s row %s is %s cells through %s mm; this raster and the walls give %d cells through %s mm" % (
                FH_COND, name_, _c[name_]["cells"], _c[name_]["thickness_mm"], cells_, _w[bore_]))
    _bore_note = ("DECLARED, not measured (no CT of this cello exists): the mean normal thickness of this specimen's top-plate-maps.tsv "
                  "(thickness_top: JOH front) along the hole's traced ring, as the shell generator reads t_holes but per hole (top-plate-fhole-walls.tsv)")
    NECK_DIMS = [("top-plate", "f-hole-%d" % h_, "bore", "%.9g" % (float(_w[(str(h_), "bore")]) / 1000), "m", "unstated", "declared", "johannsson-cello-1-1", _bore_note)
                 for h_ in (1, 2)]
    NECK_DIMS.append(("top-plate", "f-holes", "bore-pair", "%.9g" % (float(_w[("0", "bore-pair")]) / 1000), "m", "unstated", "declared", "johannsson-cello-1-1",
                      "one straight bore for both holes that keeps their parallel bore conductance: (A1 + A2) / (A1/t1 + A2/t2), A the 0.5 mm raster areas; "
                      "from the two declared bores above"))
    _p = _c["pair"]
    _layers = max(1, int(round(float(_p["thickness_mm"]) / (0.5 / 8))))
    NECK_PROV = [
        ("f_hole_bores", "hole 1 %s mm, hole 2 %s mm, pair %s mm (the dimensions.tsv bore rows)" % (_w[("1", "bore")], _w[("2", "bore")], _w[("0", "bore-pair")]),
         "DECLARED / JOH front",
         "no CT of this cello exists: each hole's bore is the mean normal thickness of this specimen's top-plate-maps.tsv along the hole's traced ring, "
         "as the shell generator reads t_holes but per hole (top-plate-fhole-walls.tsv, written by cello_fhole_walls.py); the pair bore keeps the holes' "
         "parallel straight-bore conductance, (A1 + A2) / (A1/t1 + A2/t2), A the 0.5 mm raster areas %s and %s mm2" % (_w[("1", "area-0p5mm")], _w[("2", "area-0p5mm")])),
        ("f_hole_conductance", "pair %s mm, the Richardson limit of subdivisions 1, 2, 4 and 8 (%s, %s, %s, %s mm; orders %s and %s, step %s mm); each hole alone "
         "at the pair bore %s and %s mm, at its own bore %s and %s mm" % (
             _p["K_mm"], _p["K_m1_mm"], _p["K_m2_mm"], _p["K_m4_mm"], _p["K_m8_mm"], _p["order_first"], _p["order_last"], _p["step_mm"],
             _c["hole1-at-pair-bore"]["K_mm"], _c["hole2-at-pair-bore"]["K_mm"], _c["hole1-at-own-bore"]["K_mm"], _c["hole2-at-own-bore"]["K_mm"]),
         "COMPUTED (fs_bem::plate_aperture)",
         "top-plate-fhole-conductance.tsv, written by crates/fs-bem/examples/fhole_conductance.rs (99327c3ab) from this specimen's raster moved %.2f mm in z "
         "onto the solver's grid (rows at k/2 + 0.25 mm; a translation, and the solve depends on cell offsets only) and the bores above; the pair's %d cells "
         "are %d at subdivision 8 through %d layers; the example's budgets 2^22 cells, 2^26 padded FFT points and 2^26 duct voxels, tolerance 1e-10; "
         "the plate of the solve is flat, rigid and infinite; nothing is fitted to a sound, a frequency or a page" % (
             0.0 - _off[0], int(_p["cells"]), 64 * int(_p["cells"]), _layers))]
    NECK_FILES = [("top-plate-fhole-walls.tsv", open(FH_WALLS, "rb").read()), ("top-plate-fhole-conductance.tsv", open(FH_COND, "rb").read())]

# ------------------------------------------------------------------ dimensions the generator reads
POST_ROW =(("johannsson-viola-medium", "soundpost thickness 6.5 mm (hansjohannsson.com/resources/mesviola.html, Medium column)") if VIOLA else
            ("johannsson-cello-1-1", "soundpost thickness 11 mm (hansjohannsson.com/resources/mescello.html, 1/1 column)"))
write("dimensions.tsv", "frankensim.specimen-dimensions.v1", ["part", "feature", "quantity", "value_si", "unit", "sigma_si", "route", "source_id", "note"],
      [("soundpost", "section", "equivalent-diameter-mean", "%.6f" % (POST_D / 1000), "m", "unstated", "published-maker-table") + POST_ROW] + NECK_DIMS,
      ["%s"] * 9)
for name_, bytes_ in NECK_FILES:
    with open(OUT + name_, "xb") as f:
        f.write(bytes_)
if OUTLINE:
    def _rep(item, old, new):
        i_ = [p[0] for p in PROV].index(item)
        assert old in PROV[i_][3], (item, old)
        PROV[i_] = (PROV[i_][0], PROV[i_][1], PROV[i_][2], PROV[i_][3].replace(old, new))

    _rep("bridge", "position = the Cannone bridge plane through the length scale", "position = the notch line of the placed f-hole trace (mean of the two notches' z, %.2f), NOT the plan map" % ZC2)
    _rep("soundpost", "through the plan map", "through the one landmark plan map (the warped il Cannone centroid; it does not follow the bridge plane's move to the notch line)")
    _rep("bass_bar", "about its mapped middle)", "about its mapped middle, then slid rigidly in z to equidistant ends: row bass_bar_slide)")
    _rep("bass_bar", "the plan width scale at each z", "the one landmark plan map's x at each z")
    _keep_bar = h3 > 2.0
    _bar_lo, _bar_hi = float(zb3[_keep_bar].min()), float(zb3[_keep_bar].max())
    _zw_c = float(ZF(C_LM[1]))
    _post_old = (SZ * PZ, SZ * zbr_c)  # the old uniform map's post and bridge plane z
    prov("placement_check", "bridge plane z %.2f (%.2f mm below the top edge z %.2f, %.2f mm above the front's lower ledge z %.2f); bass bar z %.2f to %.2f (%.2f mm from the top plate's bottom end, %.2f mm from its top end; length %.1f mm); soundpost z %.2f (%.2f mm %s the bridge plane; the old uniform map had %.2f mm); the warped il Cannone bridge plane would have been z %.2f (%.2f mm from the notch line)" % (
        ZC2, Z_TOP - ZC2, Z_TOP, ZC2 - (Z_BOT_EDGE + MF_LM[0]), Z_BOT_EDGE + MF_LM[0], _bar_lo, _bar_hi, _bar_lo - float(top_rows[:, 0].min()), float(top_rows[:, 0].max()) - _bar_hi, _bar_hi - _bar_lo, P2z, abs(P2z - ZC2),
        "below" if P2z < ZC2 else "above", abs(_post_old[0] - _post_old[1]), float(ZF(zbr_c)), float(ZF(zbr_c)) - ZC2), "check",
         "expected against placed, computed; not fitted. Items placed from a printed or photographed distance (f-hole trace, bridge plane = the notch line, bar length 580 and its place relative to the bridge) are listed beside the one-map warped value, so none moves by accident")
    prov("ring_overhang_check", "stretched ring on the measured outlines: mean %.2f mm, min %.2f mm; after pulling in: min %.2f mm; %d of %d angles moved, by up to %.2f mm" % (oh0.mean(), oh0.min(), oh1.min(), int((moved > 1e-9).sum()), len(moved), moved.max()), "check", "the rib_ring_overhang row above, repeated here so the byte check can compare it with the uniform stretch's (mean 5.27, min -3.57 mm, 48 of 180 angles moved by up to 24.55 mm)")
if GIB:
    for _i, _r in enumerate(PROV):
        if _r[0] in ("thickness_top", "rib_height", "ribs", "thickness_back"):
            _v = _r[1]
            if _r[0] == "thickness_top":
                _v = "edge %.1f to %.0f mm from the outline, linear to the interior by %.0f mm; interior = the Gibson belly map of the gibson_top_interior row" % (T_TOP["edge"], PURF, LINING_W)
            elif _r[0] == "thickness_back":
                _v = _v.replace("edge %.1f to" % T_BACK["edge"], "edge %.2f to" % T_BACK["edge"])
            PROV[_i] = (_r[0], _v, _r[2] + " -> GIBSON where noted in the gibson_* rows", _r[3])
    PROV += [
        ("gibson_top_interior", "47 printed belly thicknesses 2.0-3.1 mm (mean %.3f), as a piecewise-linear map" % GIB_PTS[:, 0].mean(), "hargrave-strad-1986-09-gibson-viola p. 117",
         "REPLACES Johannsson's flank 2.6 / island 3.2 as the top interior; each point at its printed position on the Gibson belly outline, (u, v) = (distance from the top of the printed outline over its height, offset from the outline's centre over half its width on that row), bass side on the viewer's left; the map is the Delaunay triangulation (triangle, default switches) of the points' (u, v) with barycentric linear interpolation, the nearest printed point outside their hull; a plate point's (u, v) is read the same way on the mapped viola top outline; no smoothing, no weight, nothing fitted; the top edge stays 4.0 (Hargrave prints edge thickness from the back only)"),
        ("gibson_rib_height", "36.3125 mm", "hargrave-strad-1986-09-gibson-viola p. 116", "REPLACES Johannsson's 38.0: the mean of the eight printed rib heights (neck root 35.5/35.5, upper corner 35.75/36.75, lower corner 36.5/37, endpin 36.5/36.5)"),
        ("gibson_rib_wall", "1.0 mm", "hargrave-strad-1986-09-gibson-viola p. 115", "REPLACES Johannsson's 1.2: measured in the C bout through the f hole, 'approximately 1 mm'"),
        ("gibson_back_edge", "4.75 mm", "hargrave-strad-1986-09-gibson-viola p. 116", "REPLACES Johannsson's 4.0 as the back's edge thickness: the printed value at the upper and lower bouts (corners 4.75-5, centre bouts 4.8-5, button 4.5-5); the back interior stays Johannsson/Cannone")]
PROV += NECK_PROV
with open(OUT + "PROVENANCE.tsv", "x", encoding="utf-8") as f:
    f.write("item\tvalue\tsource\tnote\n")
    for r in PROV:
        f.write("\t".join(r) + "\n")
for r in PROV:
    print(r[0], "=", r[1])
