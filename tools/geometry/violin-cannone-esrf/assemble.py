"""assemble.py OUT_TABLES_DIR FINAL_DIR : builds the repo directory data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026/ from the tables written by build_tables.py,
build_tables2.py and string_kinks.py: the six standard specimen tables (specimen, frames, dimensions, parts, unknowns; fingerings do not apply) with the dimension values computed
from the data tables, and every data table under its own schema line. Everything numeric is computed here from the tables; nothing is typed in except the text of the identity,
parts and unknowns tables."""
import os
import sys

import numpy as np

src, dst = sys.argv[1], sys.argv[2]
os.makedirs(dst, exist_ok=True)
SPEC = "violin-cannone-guarneri-1743-esrf2026"
SRC = "esrf-il-cannone-ct-2026"
PROV = "ESRF data collection DOI 10.15151/ESRF-DC-2413577471 (CC BY 4.0), dataset 2416550332 (strings under tension)"


def rd(name):
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(src, name)) if not l.startswith("#")]
    if rows[0][0].startswith("frankensim."):
        rows = rows[1:]
    return rows[0], np.array([[float(x) if x != "" else np.nan for x in r] for r in rows[1:]])


def rd_text(name):
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(src, name)) if not l.startswith("#")]
    if rows[0][0].startswith("frankensim."):
        rows = rows[1:]
    return rows[0], rows[1:]


# --- data tables with their own schema line
SCHEMAS = {
    "top-plate-outline.tsv": "frankensim.specimen-outline-rows.v1", "back-plate-outline.tsv": "frankensim.specimen-outline-rows.v1",
    "top-plate-fholes.tsv": "frankensim.specimen-raster-runs.v1", "bridge-silhouette.tsv": "frankensim.specimen-raster-runs.v1",
    "top-plate-maps.tsv": "frankensim.specimen-plate-maps.v1", "back-plate-maps.tsv": "frankensim.specimen-plate-maps.v1",
    "ribs-sections.tsv": "frankensim.specimen-ring-sections.v1", "linings-sections.tsv": "frankensim.specimen-ring-sections.v1",
    "bass-bar-sections.tsv": "frankensim.specimen-sections.v1", "soundpost-sections.tsv": "frankensim.specimen-sections.v1", "bridge-sections.tsv": "frankensim.specimen-sections.v1",
    "neck-sections.tsv": "frankensim.specimen-sections.v1", "fingerboard-sections.tsv": "frankensim.specimen-sections.v1", "tail-block-sections.tsv": "frankensim.specimen-sections.v1",
    "string-G-sections.tsv": "frankensim.specimen-sections.v1", "string-D-sections.tsv": "frankensim.specimen-sections.v1", "string-A-sections.tsv": "frankensim.specimen-sections.v1",
    "string-E-sections.tsv": "frankensim.specimen-sections.v1", "string-contacts.tsv": "frankensim.specimen-string-contacts.v1", "part-volumes.tsv": "frankensim.specimen-part-volumes.v1",
}
RENAME = {"heel-block-sections.tsv": "tail-block-sections.tsv", "string-G-sections.tsv": "string-g-sections.tsv", "string-D-sections.tsv": "string-d-sections.tsv",
          "string-A-sections.tsv": "string-a-sections.tsv", "string-E-sections.tsv": "string-e-sections.tsv"}
written = []
for name in sorted(os.listdir(src)):
    if name in ("dimensions.tsv", "SHA256.txt") or not name.endswith(".tsv"):
        continue
    out_name = RENAME.get(name, name)
    schema = SCHEMAS.get(name) or SCHEMAS.get(out_name) or SCHEMAS.get(name.replace("string-g", "string-G"))
    if schema is None:
        raise SystemExit("no schema for " + name)
    head, rows = rd_text(name)
    if name == "part-volumes.tsv":  # part names to slugs
        rows = [[{"heel-block": "tail-block"}.get(r[0], r[0]).lower()] + r[1:] for r in rows]
    if name == "string-contacts.tsv":
        rows = [[r[0].lower()] + r[1:] for r in rows]
    with open(os.path.join(dst, out_name), "w") as f:
        f.write(schema + "\n" + "\t".join(head) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    written.append(out_name)

# --- dimensions
D = []


def add(part, feature, quantity, value, unit, sigma, note):
    D.append((part, feature, quantity, float(value), unit, float(sigma), note))


HALF = 0.25
mm, mm2, mm3 = 1e-3, 1e-6, 1e-9
for pl in ("top-plate", "back-plate"):
    h, o = rd("%s-outline.tsv" % pl)
    z, xl, xr = o[:, 0], o[:, 1], o[:, 2]
    w = xr - xl
    z0, z1 = z.min(), z.max()
    add(pl, "footprint", "length-z", (z1 - z0 + 1.0) * mm, "m", 0.5 * mm, "extent in z of the outline rows (1 mm rows), including the end features of the plate")
    add(pl, "footprint", "width-max", w.max() * mm, "m", 2 * HALF * mm, "widest outline row (right edge minus left edge of the occupied cells); row at z = %.2f mm" % z[w.argmax()])
    third = (z1 - z0) / 3
    for lab, (a, b), fn in (("lower-bout-width-max", (z0, z0 + third), np.argmax), ("waist-width-min", (z0 + third, z0 + 2 * third), np.argmin), ("upper-bout-width-max", (z0 + 2 * third, z1), np.argmax)):
        m = (z >= a) & (z <= b) & (w > 40)
        i = np.nonzero(m)[0][fn(w[m])]
        add(pl, "outline", lab, w[i] * mm, "m", 2 * HALF * mm, "extreme of the outline width within the tail-end, middle or neck-end third of the plate")
        add(pl, "outline", lab.replace("-width-max", "-z-position").replace("-width-min", "-z-position"), z[i] * mm, "m", 0.5 * mm, "z of the row of the preceding width")
    mid = 0.5 * (xl + xr)
    add(pl, "footprint", "midline-x-mean", mid[w > 100].mean() * mm, "m", 0.5 * mm, "mean of the outline midpoints over rows wider than 100 mm: offset of the plate's mid line from the x origin of the scan frame")
    h, m = rd("%s-maps.tsv" % pl)
    x, zz, yo, yi, tn, dist, nv = m.T
    inter = dist >= 8.0
    add(pl, "thickness", "normal-median-interior", np.median(tn[inter]) * mm, "m", 0.1 * mm, "median of the 3 mm block values farther than 8 mm from the outline; thickness measured normal to the outer surface")
    add(pl, "thickness", "normal-p05-interior", np.percentile(tn[inter], 5) * mm, "m", 0.1 * mm, "5th percentile of the same cells")
    add(pl, "thickness", "normal-p95-interior", np.percentile(tn[inter], 95) * mm, "m", 0.1 * mm, "95th percentile of the same cells")
    i = np.argmax(np.where(inter, tn, -1))
    add(pl, "thickness", "normal-max-interior", tn[i] * mm, "m", 0.1 * mm, "largest block value of the same cells, at x = %.1f mm, z = %.1f mm" % (x[i], zz[i]))
    edge = (dist >= 1.0) & (dist <= 3.0)
    ref = np.median(yo[edge])
    ic = np.argmin(yo) if pl == "top-plate" else np.argmax(yo)
    add(pl, "arching", "outer-crown-height-above-edge-plane", abs(yo[ic] - ref) * mm, "m", 0.5 * mm,
        "outer surface at the crown (x = %.1f mm, z = %.1f mm) against the median outer surface of the 1 to 3 mm band inside the outline; the plate thickness is NOT added; other authors may define arch height differently" % (x[ic], zz[ic]))
    add(pl, "arching", "edge-plane-y", ref * mm, "m", 0.5 * mm, "median outer-surface y of the 1 to 3 mm band inside the outline (scan frame)")

h, fh = rd("top-plate-fholes.tsv")
for hid in (1, 2):
    s = fh[fh[:, 0] == hid]
    zf, a, b = s[:, 1], s[:, 2], s[:, 3]
    area = ((b - a) * 0.5).sum()
    add("top-plate", "f-hole-%d" % hid, "area", area * mm2, "m2", 5 * mm2, "area of the unoccupied cells inside the outline (0.5 mm raster); hole 1 lies at negative x (the bass-bar side), hole 2 at positive x")
    add("top-plate", "f-hole-%d" % hid, "centroid-x", (((a + b) / 2) * (b - a)).sum() / (b - a).sum() * mm, "m", 0.5 * mm, "area centroid of the raster")
    add("top-plate", "f-hole-%d" % hid, "centroid-z", (zf * (b - a)).sum() / (b - a).sum() * mm, "m", 0.5 * mm, "area centroid of the raster")
    add("top-plate", "f-hole-%d" % hid, "z-extent", (zf.max() - zf.min() + 0.5) * mm, "m", 0.5 * mm, "z from %.2f to %.2f mm" % (zf.min() - 0.25, zf.max() + 0.25))
    add("top-plate", "f-hole-%d" % hid, "x-extent", (b.max() - a.min()) * mm, "m", 0.5 * mm, "x from %.2f to %.2f mm" % (a.min(), b.max()))
    X = np.concatenate([np.stack([np.linspace(aa, bb, max(int((bb - aa) / 0.25), 1)), np.full(max(int((bb - aa) / 0.25), 1), zz_)], 1) for zz_, aa, bb in zip(zf, a, b)])
    ev, evec = np.linalg.eigh(np.cov((X - X.mean(0)).T))
    ax = evec[:, 1] if evec[1, 1] >= 0 else -evec[:, 1]
    proj = (X - X.mean(0)) @ ax
    add("top-plate", "f-hole-%d" % hid, "length-along-principal-axis", (proj.max() - proj.min()) * mm, "m", 0.5 * mm, "extent along the principal axis of the raster")
    add("top-plate", "f-hole-%d" % hid, "principal-axis-angle-from-z", np.arctan2(ax[0], ax[1]), "rad", np.radians(1.0), "angle of the principal axis from +z toward +x")

h, b = rd("bass-bar-sections.tsv")
zb, ar, cx, cy, x0, x1, y0, y1 = b.T
pv = {{"heel-block": "tail-block"}.get(r[0], r[0]).lower(): r for r in rd_text("part-volumes.tsv")[1]}
def ext(name, ax):
    r = pv[name]
    return float(r[3 + 2 * ax + 1]) - float(r[3 + 2 * ax])
add("bass-bar", "extent", "length-z", ext("bass-bar", 2) * mm, "m", 0.1 * mm, "z extent of the mesh")
i = ar.argmax()
add("bass-bar", "section", "area-max", ar[i] * mm2, "m2", 1 * mm2, "largest section area (planes every 5 mm), at z = %.2f mm" % zb[i])
add("bass-bar", "section", "width-x-at-area-max", (x1[i] - x0[i]) * mm, "m", 0.2 * mm, "x extent of that section")
add("bass-bar", "section", "height-y-at-area-max", (y1[i] - y0[i]) * mm, "m", 0.2 * mm, "y extent of that section")
add("bass-bar", "position", "centroid-x-at-area-max", cx[i] * mm, "m", 0.2 * mm, "section centroid, scan frame")
add("bass-bar", "position", "centroid-x-at-bridge-z", np.interp(-106.0, zb, cx) * mm, "m", 0.3 * mm, "section centroid interpolated at the string-bridge contact z (-106.0 mm)")

h, p = rd("soundpost-sections.tsv")
yp, ar, cx, cz, x0, x1, z0_, z1_ = p.T
add("soundpost", "extent", "length-y", ext("soundpost", 1) * mm, "m", 0.1 * mm, "y extent of the mesh")
add("soundpost", "section", "equivalent-diameter-mean", 2 * np.sqrt(ar.mean() / np.pi) * mm, "m", 0.1 * mm, "2 sqrt(mean section area / pi); mean area %.2f mm2" % ar.mean())
add("soundpost", "axis", "tilt-x-per-y", np.polyfit(yp, cx, 1)[0], "1", 0.002, "dx/dy of the section centroids (dimensionless)")
add("soundpost", "axis", "tilt-z-per-y", np.polyfit(yp, cz, 1)[0], "1", 0.002, "dz/dy of the section centroids (dimensionless)")
add("soundpost", "position", "centroid-x-at-mid-length", np.interp(yp.mean(), yp, cx) * mm, "m", 0.2 * mm, "scan frame")
add("soundpost", "position", "centroid-z-at-mid-length", np.interp(yp.mean(), yp, cz) * mm, "m", 0.2 * mm, "scan frame")

h, g = rd("bridge-sections.tsv")
yb, ar, cx, cz, x0, x1, z0_, z1_ = g.T
add("bridge", "extent", "height-y", ext("bridge", 1) * mm, "m", 0.1 * mm, "y extent of the mesh")
add("bridge", "extent", "width-x-max", (x1 - x0).max() * mm, "m", 0.3 * mm, "largest x extent of any section (planes every 1 mm)")
add("bridge", "extent", "thickness-z-max", (z1_ - z0_).max() * mm, "m", 0.2 * mm, "largest z extent of any section")
add("bridge", "extent", "thickness-z-min-above-y-minus-52", (z1_ - z0_)[yb < -52].min() * mm, "m", 0.2 * mm, "smallest z extent among the sections above y = -52 mm")
add("bridge", "position", "z-mid-plane-mean", ((z0_ + z1_) / 2)[ar > 50].mean() * mm, "m", 0.2 * mm, "mean of the mid z over sections of area above 50 mm2")

# ribs and linings
h, rb = rd("ribs-sections.tsv")
th, cxr, czr, area, rmin, rmax, ymin, ymax, mw = rb.T
add("ribs", "section", "height-y-median", np.median(ymax - ymin) * mm, "m", 0.3 * mm, "median over the 180 half-plane sections of y_max - y_min")
add("ribs", "section", "radial-mean-width-p10", np.percentile(mw, 10) * mm, "m", 0.2 * mm, "10th percentile of section area / height; a RADIAL section overstates the wall thickness where the wall is oblique to the radial direction")
add("ribs", "section", "radial-mean-width-median", np.median(mw) * mm, "m", 0.2 * mm, "median of the same quantity")
add("ribs", "extent", "x-extent", ext("ribs", 0) * mm, "m", 0.1 * mm, "x extent of the mesh")
add("ribs", "extent", "z-extent", ext("ribs", 2) * mm, "m", 0.1 * mm, "z extent of the mesh")
add("ribs", "extent", "y-extent", ext("ribs", 1) * mm, "m", 0.1 * mm, "y extent of the mesh")
h, ln = rd("linings-sections.tsv")
t_area, t_rmin, t_rmax, t_ymin, t_ymax, t_mw, b_area, b_rmin, b_rmax, b_ymin, b_ymax, b_mw = ln[:, 3:].T
small = (t_area < 20) & (b_area < 25)
add("linings", "top-lining", "height-y-median", np.nanmedian((t_ymax - t_ymin)[small]) * mm, "m", 0.3 * mm, "median over the sections that are not corner blocks (top area < 20 mm2 and bottom area < 25 mm2)")
add("linings", "top-lining", "radial-mean-width-median", np.nanmedian(t_mw[small]) * mm, "m", 0.2 * mm, "section area / height, radial section")
add("linings", "bottom-lining", "height-y-median", np.nanmedian((b_ymax - b_ymin)[small]) * mm, "m", 0.3 * mm, "same selection")
add("linings", "bottom-lining", "radial-mean-width-median", np.nanmedian(b_mw[small]) * mm, "m", 0.2 * mm, "section area / height, radial section")
add("linings", "corner-blocks", "sections-merged-into-the-mesh", int((~small).sum()), "1", 0.0, "number of the 2-degree half-plane sections whose top or bottom area is that of a block (top >= 20 mm2 or bottom >= 25 mm2): corner and end blocks are merged into this mesh")

# strings
h, sc = rd_text("string-contacts.tsv")
cont = {(r[0], r[1]): r for r in sc}
for s in ("G", "D", "A", "E"):
    h2, st = rd("string-%s-sections.tsv" % s)
    mask = st[:, 1] < 3 * np.median(st[:, 1])
    dia = 2 * np.sqrt(np.median(st[mask, 1]) / np.pi)
    add("string-%s" % s.lower(), "section", "equivalent-diameter-median", dia * mm, "m", 0.02 * mm, "2 sqrt(median section area / pi) over the stations along the string (planes every 20 mm); envelope of the segmented cross-section, windings included where present")
    zb_ = float(cont[(s, "bridge")][2])
    add("string-%s" % s.lower(), "bridge-contact", "z", zb_ * mm, "m", 1.0 * mm, "z of the largest bend of the centre line near the bridge (1 mm sections, 12-station slope fits)")
    add("string-%s" % s.lower(), "bridge-contact", "x", float(cont[(s, "bridge")][3]) * mm, "m", 0.3 * mm, "x of the centre line at that z")
    add("string-%s" % s.lower(), "bridge-contact", "y", float(cont[(s, "bridge")][4]) * mm, "m", 0.3 * mm, "y of the centre line at that z")
    add("string-%s" % s.lower(), "nut-contact", "z", float(cont[(s, "nut")][2]) * mm, "m", 1.0 * mm, "z of the largest bend of the centre line near the nut" + (" (the G detection differs by 8 mm from the other three and is not used for the vibrating length)" if s == "G" else ""))
nut_z = np.mean([float(cont[(s, "nut")][2]) for s in "DAE"])
br_z = np.mean([float(cont[(s, "bridge")][2]) for s in "GDAE"])
add("fingerboard", "nut-to-bridge", "string-vibrating-length", (nut_z - br_z) * mm, "m", 1.5 * mm, "mean nut contact z of the D, A and E strings (%.1f mm) minus the mean bridge contact z of the four strings (%.1f mm); the fingerboard mesh ends at z = %.2f mm" % (nut_z, br_z, float(pv["fingerboard"][8])))

# fingerboard, neck, tail block, tailpiece, volumes
h, fb = rd("fingerboard-sections.tsv")
add("fingerboard", "extent", "length-z", ext("fingerboard", 2) * mm, "m", 0.1 * mm, "z extent of the mesh")
add("fingerboard", "section", "width-x-first-station", (fb[0, 5] - fb[0, 4]) * mm, "m", 0.3 * mm, "x extent at the station nearest the body (z = %.2f mm)" % fb[0, 0])
add("fingerboard", "section", "width-x-last-station", (fb[-1, 5] - fb[-1, 4]) * mm, "m", 0.3 * mm, "x extent at the station nearest the nut (z = %.2f mm)" % fb[-1, 0])
h, nk = rd("neck-sections.tsv")
add("neck", "extent", "z-extent", ext("neck", 2) * mm, "m", 0.1 * mm, "z extent of the mesh (it includes the neck foot, pegbox and scroll)")
add("tailpiece", "extent", "z-extent", ext("tailpiece", 2) * mm, "m", 0.1 * mm, "z extent of the mesh")
add("tail-block", "extent", "z-extent", ext("tail-block", 2) * mm, "m", 0.1 * mm, "z extent of the mesh (file talon_bloc_principal)")
add("tail-block", "extent", "x-extent", ext("tail-block", 0) * mm, "m", 0.1 * mm, "x extent of the mesh")
VOL_PARTS = ["top-plate", "back-plate", "ribs", "linings", "bass-bar", "soundpost", "bridge", "tailpiece", "neck", "fingerboard", "tail-block", "string-g", "string-d", "string-a", "string-e"]
for name in VOL_PARTS:
    r = pv[name]
    add(name, "mesh", "enclosed-volume", float(r[1]) * mm3, "m3", 0.01 * float(r[1]) * mm3, "signed volume of the closed part mesh (a few non-manifold edges in the tailpiece); 1 percent is a stated allowance, not a measured bound")

with open(os.path.join(dst, "dimensions.tsv"), "w") as f:
    f.write("frankensim.specimen-dimensions.v1\n")
    f.write("part\tfeature\tquantity\tvalue_si\tunit\tsigma_si\troute\tsource_id\tnote\n")
    for part, feat, q, v, unit, s, note in D:
        f.write("%s\t%s\t%s\t%.7g\t%s\t%.3g\tpublished-measured\t%s\t%s\n" % (part, feat, q, v, unit, s, SRC, note.replace("\t", " ")))
written.append("dimensions.tsv")

# --- specimen, frames, parts, unknowns
sha = {}
for line in open(os.path.join(src, "SHA256.txt")):
    h_, n_ = line.split()
    sha[n_] = h_
files_used = ";".join("%s=%s" % (n, h_[:16]) for n, h_ in sorted(sha.items()))
FRAME = ["the centre of the exported whole-violin surface of the ESRF segmentation (VGStudio export in mm; STL bounding box x -101.97..101.67, y -54.71..55.33, z -293.95..293.04 mm)",
         "+x across the body toward the treble side (the E string and the soundpost lie at positive x, the G string and the bass bar at negative x)",
         "+y from the strings side toward the back (the top plate lies at negative y, the back plate at positive y)",
         "+z along the body from the tail end toward the scroll"]
rows = [
    ("specimen_id", SPEC),
    ("maker", "Bartolomeo Giuseppe Guarneri 'del Gesu', Cremona (the producer field of the ESRF record)"),
    ("model", "'il Cannone', the violin of Niccolo Paganini (ESRF record title)"),
    ("serial", "unstated"),
    ("year", "1743 (ESRF record)"),
    ("system", "four-string violin scanned with its strings tensioned (the dataset the ESRF names accorded); repairs and later additions are present, see unknowns.tsv"),
    ("holder", "Palazzo Doria-Tursi, Genoa city hall (the repository institution of the ESRF record); not accessible to FrankenSim"),
    ("reference_state", "strung"),
    ("environment", "unstated; scanned at the ESRF beamline BM18, session 2024-03-09 to 2024-03-11, proposal HG-236"),
    ("source_ids", SRC),
    ("target_relation", "ONE REAL violin from a published CT study for pexfp.199 (mail 5226). A 1743 del Gesu, NOT the UIowa violin of the recordings, whose own dimensions stay UNKNOWN (u-not-iowa-violin); population-grade stand-in only"),
    ("doi", "10.15151/ESRF-DC-2413577471"),
    ("licence", "CC-BY-4.0 (ESRF DataCite record, issued 2026); attribution: ESRF; data authors Paul Tafforeau, Luigi Paolasini and Alessandro Mirone, segmentation by Paul Tafforeau (the scan metadata file); the record states no licence of the instrument's owner beyond the ESRF's"),
    ("scan", "propagation phase-contrast synchrotron tomography, helical, 12000 projections, 35.66 um voxel, volume 5947 x 3165 x 16933 voxels (the scan metadata file of dataset 2416550332)"),
    ("mesh_resolution", "measured on the fetched meshes with tools/geometry/violin-cannone-esrf/mesh_resolution.py: the part meshes are decimated copies of the CT surface, not voxel-sized: median edge 0.19 mm for the top plate (5.4 voxels, 90th percentile 0.46 mm), 0.16 mm for the back plate, 0.17 mm for the ribs; the whole-violin STL named 200K holds 3508724 triangles and the one named 10M holds 10027320, so the names are not triangle counts; their volumes differ by 0.06 percent (694157 and 694565 mm3), their areas by 0.4 percent and their sections at z = -100 and -190 mm, y = 5 mm and x = 30 mm by at most 0.09 percent; every part mesh checked is closed (each edge used twice, a few exceptions in the tailpiece) and the STL header states Unit:Millimeter"),
    ("source_files", "files of dataset 2416550332 under surface_models/ (and the STL of the whole violin, used only for the frame): name=first 16 hex of the SHA-256 of the file as fetched: " + files_used),
    ("units_and_frame", "every length in the data tables is in millimetres in the frame of the source files (frames.tsv); dimensions.tsv is in SI"),
    ("extraction", "this lane's numpy-only scripts (tools/geometry/violin-cannone-esrf/): plate outlines, holes and thickness from a 0.5 mm raster of the surface samples of each closed plate mesh (3 mm blocks in the map tables); sections of bar, post, bridge, neck, fingerboard, tail block and strings by planes; rib and lining sections by half-planes through the vertical axis every 2 degrees"),
    ("cross_checks", "the 4 string meshes give equivalent diameters 1.095, 1.094, 0.826 and 0.266 mm and the soundpost 6.35 mm (1/4 inch); the string bridge contacts agree to 0.2 mm in z and the nut contacts of three strings to 0.2 mm; the STL header states Unit:Millimeter"),
]
with open(os.path.join(dst, "specimen.tsv"), "w") as f:
    f.write("frankensim.specimen.v1\nkey\tvalue\n")
    for k, v in rows:
        f.write("%s\t%s\n" % (k, v))
written.append("specimen.tsv")

PARTS = [
    ("top-plate", "unstated in the ESRF record (spruce is usual for the type but not stated)", "assembled", "ribs;bass-bar;bridge;soundpost;linings", "file table_all.glb; closed shell of 7190594 triangles holding the outer and inner surfaces, the f-holes and the purfling channel; whether varnish is inside it is not established (u-varnish-and-purfling)"),
    ("back-plate", "unstated in the ESRF record (maple is usual for the type but not stated)", "assembled", "ribs;soundpost;linings;neck;tail-block", "file dos_all.glb; closed shell of 10389200 triangles including the neck button; one piece in this scan's maps (the record does not say whether the back is one or two pieces)"),
    ("ribs", "unstated in the ESRF record", "assembled", "top-plate;back-plate;linings;tail-block", "file eclisses_all.glb; the six rib strips are one mesh"),
    ("linings", "unstated in the ESRF record", "assembled", "ribs;top-plate;back-plate", "file contre_elisses.glb (the French name for the counter-ribs); top and bottom linings merged with the corner blocks (u-body-corner-blocks)"),
    ("bass-bar", "unstated in the ESRF record", "assembled", "top-plate", "file barre_harmonie.glb; a single closed mesh; whether it is the original bar is not stated (u-repairs-and-modifications)"),
    ("soundpost", "unstated in the ESRF record", "assembled", "top-plate;back-plate", "file ame.glb; a single closed mesh"),
    ("bridge", "unstated in the ESRF record", "assembled", "strings;top-plate", "file chevalet.glb; a single closed mesh"),
    ("tailpiece", "unstated in the ESRF record", "assembled", "strings", "file cordier.glb; closed mesh with 4 non-manifold edges out of 1.5 million"),
    ("neck", "unstated in the ESRF record", "assembled", "fingerboard;back-plate;tail-block", "file manche_all.glb: the neck with its foot, pegbox and scroll; separate meshes of later neck additions exist in the dataset and are not extracted (u-repairs-and-modifications)"),
    ("fingerboard", "unstated in the ESRF record", "assembled", "neck;strings", "file touche.glb"),
    ("tail-block", "unstated in the ESRF record", "assembled", "ribs;top-plate;back-plate", "file talon_bloc_principal.glb: located at the tail end (z -261.57..-250.52 mm) and as tall as the ribs, so read here as the tail (end) block; the file's own label is talon"),
    ("string-g", "unstated in the ESRF record", "strung", "bridge;neck", "file corde_sol.glb"),
    ("string-d", "unstated in the ESRF record", "strung", "bridge;neck", "file corde_re.glb"),
    ("string-a", "unstated in the ESRF record", "strung", "bridge;neck", "file corde_la.glb"),
    ("string-e", "unstated in the ESRF record", "strung", "bridge;neck", "file corde_mi.glb"),
]
with open(os.path.join(dst, "parts.tsv"), "w") as f:
    f.write("frankensim.specimen-parts.v1\npart\tmaterial\tstate\tjoins\tnote\n")
    for r in PARTS:
        f.write("\t".join(r) + "\n")
written.append("parts.tsv")
with open(os.path.join(dst, "frames.tsv"), "w") as f:
    f.write("frankensim.specimen-frames.v1\npart\torigin\taxis_x\taxis_y\taxis_z\tnote\n")
    for r in PARTS:
        f.write("%s\t%s\t%s\t%s\t%s\tone assembly frame for every part; the x sign convention is inferred from the string meshes\n" % ((r[0],) + tuple(FRAME)))
written.append("frames.tsv")
UNK = [
    ("u-segmentation-error", "the authors' threshold choice and mesh decimation of the exported surfaces: the record states no uncertainty for them; the sigma values in dimensions.tsv bound only this lane's extraction from the meshes; the meshes are decimated (median edge 4 to 5 voxels, mesh_resolution in specimen.tsv)", "unbounded beyond the 35.66 um voxel", "scan", "any claim finer than about 0.1 mm and thickness-sensitivity studies", "the authors' 8-bit reduced volume (dataset 2416550332, 7.06 GB) or the full 35.66 um volume to cross-check mesh thickness at sections"),
    ("u-tension-state", "the scan has the strings tensioned; the free geometry differs by the displacements the authors measure by digital volume correlation; the slack-string scan (dataset 2416532883) has volumes only and no segmentation", "unbounded", "state", "unloaded-reference claims, neck and plate deformation", "segment the slack-string volume (71.32 um JPEG-2000, 8 GB) or use the authors' displacement fields"),
    ("u-wood-material", "no wood species, density, grain direction or moisture is in the record and the grey values are not density-calibrated in the published files", "unbounded", "materials", "stiffness, mass and loss of every wooden part", "sourced spruce and maple tables with the model-form mismatch stated; a measured density of this violin is not available"),
    ("u-repairs-and-modifications", "the authors' part meshes show repairs, nails, later neck additions, a repaired fracture at an f-hole and the effects of age; none of these parts is extracted here and their effect on local stiffness is not separated", "unbounded", "history", "pristine-instrument claims and local plate stiffness near repairs", "the authors' part meshes (manche_rajout, clous_*, talon_reparation_base) and their papers"),
    ("u-varnish-and-purfling", "the plate meshes named _all may include the varnish layer (separate _vernis meshes exist for the top, back, ribs and neck) and the purfling inlay; whether the thickness here includes varnish is not established", "varnish thickness unknown", "surface", "plate thickness and mass to better than the varnish thickness", "compare table_all with table_vernis and table_mixed in the dataset"),
    ("u-body-corner-blocks", "the corner blocks and the lining strips are one mesh (contre_elisses); the upper block belongs to the neck mesh; neither is separated here", "unbounded", "assembly", "block mass and stiffness", "cut the combined meshes by the block positions"),
    ("u-bass-treble-side", "the sign convention of x (bass at negative x) is inferred from the string meshes (G at negative x, E at positive x), not stated by the record", "categorical", "frame", "a mirror error in any consumer of the maps", "none needed if the strings are read first"),
    ("u-not-iowa-violin", "this violin is not the UIowa MIS violin of the recordings; no same-specimen claim against them", "categorical", "identity", "same-specimen claims against the recordings", "owner measurement of the UIowa violin or population-level comparison only"),
    ("u-contacts", "bridge feet against the top plate, post ends against the plates, bar glue line and string grooves at the nut and bridge are separate meshes and the contact fits are not extracted", "unbounded", "contacts", "contact stiffness and gap claims", "section the paired meshes at the contact planes"),
    ("u-neck-pegbox-scroll", "pegbox, scroll, pegs, nails and the nut are not extracted; the nut position is taken from the string bends and the end of the fingerboard mesh", "unbounded", "neck", "neck-end radiation and nut-boundary claims", "the neck mesh and the peg meshes of the dataset"),
    ("u-tail-end-features", "the tail end of the top plate carries a rectangular notch and the back plate a button at the neck end; end pin, tailgut and chinrest are not in the extracted meshes", "unbounded", "tail", "tail-end boundary conditions", "the other meshes of the dataset"),
]
with open(os.path.join(dst, "unknowns.tsv"), "w") as f:
    f.write("frankensim.specimen-unknowns.v1\nid\twhat\tbound\tcorrelation_group\tblocks\troute\n")
    for r in UNK:
        f.write("\t".join(r) + "\n")
written.append("unknowns.tsv")
print("files:", sorted(written))
print(len(D), "dimension rows")
