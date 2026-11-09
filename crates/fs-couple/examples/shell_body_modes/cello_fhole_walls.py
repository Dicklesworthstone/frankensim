"""The f-hole wall table of a cello outline specimen, and the raster the conductance solve reads.

The cello has no CT, so its bores are DECLARED: each hole's bore is the mean normal thickness of the specimen's
top-plate-maps.tsv along that hole's traced ring, read exactly as gen_violin_shell_assembly.py reads t_holes
(PlateField('top').t over f_hole_rings(0.6 H), with the options of shell-assembly-h12.5-outline-r24.tsv), but per hole.
The pair bore keeps the holes' parallel straight-bore conductance (R_i = t_i / A_i): (A1 + A2) / (A1 / t1 + A2 / t2),
A the 0.5 mm raster areas. The generator's pooled t_holes is written beside them (the table's neck meta prints it).

crates/fs-bem/examples/fhole_conductance.rs reads raster rows at z = k/2 + 0.25 mm only. The outline specimen's rows sit
at another offset of the same 0.5 mm grid, so the solver's raster is the specimen raster with every z moved onto that
grid and nothing else changed; the solve (a flat, rigid, infinite plate) depends on cell offsets only, so the move does
not change any conductance.

usage: uv run --python 3.12 --with triangle --with numpy python3 -I cello_fhole_walls.py SPECIMEN_DIR OUT_DIR
writes OUT_DIR/top-plate-fhole-walls.tsv and OUT_DIR/top-plate-fholes-solver-grid.tsv (OUT_DIR must exist; nothing is
overwritten). check_cello_shell.py N24 compares the first with the committed top-plate-fhole-walls-outline-r24.tsv.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SPEC, OUT = sys.argv[1].rstrip("/"), Path(sys.argv[2])
H = 12.5
sys.argv = ["gen", SPEC, "/dev/null", str(H), "neck=ellipse", "foot_band_mm=8", "fanfix=move", "crossing_snap_mm=2.5",
            "bridge_dp_mm=0.5", "rib_layers=24", "instrument=cello"]
spec = importlib.util.spec_from_file_location("gen", HERE / "gen_violin_shell_assembly.py")
gen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gen)

fld = gen.PlateField("top", -1)
rings, _ = gen.f_hole_rings(0.6 * H)
t_pool = float(np.mean(fld.t(gen.to_plan(np.vstack(rings)))))
t = [float(np.mean(fld.t(gen.to_plan(r)))) for r in rings]
t_min = [float(np.min(fld.t(gen.to_plan(r)))) for r in rings]
t_max = [float(np.max(fld.t(gen.to_plan(r)))) for r in rings]

# the raster: areas, and the copy on the solver's grid
lines = open(SPEC + "/top-plate-fholes.tsv", encoding="utf-8").read().splitlines()
if lines[:2] != ["frankensim.specimen-raster-runs.v1", "hole\tz_mm\tx_start_mm\tx_end_mm"]:
    raise SystemExit("refusing: %s/top-plate-fholes.tsv is not a raster-runs table" % SPEC)
runs = [l.split("\t") for l in lines[2:] if l.strip()]
offsets = sorted({round((float(z) - 0.25) % 0.5, 6) % 0.5 for _, z, _, _ in runs})
if len(offsets) != 1:
    raise SystemExit("refusing: the raster rows are not on one 0.5 mm grid (offsets %s)" % offsets)
shift = 0.0 - offsets[0]
area = {1: 0.0, 2: 0.0}
out_rows = []
for h, z, a, b in runs:
    for v in (float(a), float(b)):
        if abs(round(v / 0.5) * 0.5 - v) > 1e-9:
            raise SystemExit("refusing: x %r is off the 0.5 mm grid" % v)
    area[int(h)] += float(b) - float(a)
    zs = float(z) + shift
    assert abs(round((zs - 0.25) / 0.5) * 0.5 + 0.25 - zs) < 1e-9, (z, zs)
    out_rows.append("%s\t%.2f\t%s\t%s" % (h, zs, a, b))
area = {h: v * 0.5 for h, v in area.items()}
pair = (area[1] + area[2]) / (area[1] / t[0] + area[2] / t[1])

with open(OUT / "top-plate-fholes-solver-grid.tsv", "x", encoding="utf-8") as f:
    f.write("frankensim.specimen-raster-runs.v1\nhole\tz_mm\tx_start_mm\tx_end_mm\n")
    f.write("\n".join(out_rows) + "\n")

note_bore = ("DECLARED, not measured: mean normal thickness of top-plate-maps.tsv along the hole's traced ring (the shell "
             "generator's PlateField t over f_hole_rings at 7.5 mm spacing, %d points; min %.3f, max %.3f mm)")
with open(OUT / "top-plate-fhole-walls.tsv", "x", encoding="utf-8") as f:
    f.write("frankensim.specimen-fhole-walls.v1\nhole\tquantity\tvalue\tunit\tnote\n")
    for h in (1, 2):
        f.write("%d\tarea-0p5mm\t%.2f\tmm2\tthe raster's 0.5 mm cells\n" % (h, area[h]))
        f.write("%d\tbore\t%.6g\tmm\t%s\n" % (h, t[h - 1], note_bore % (len(rings[h - 1]), t_min[h - 1], t_max[h - 1])))
    f.write("0\tbore-pair\t%.6g\tmm\tone bore for both holes with their parallel straight-bore conductance: (A1 + A2) / (A1/t1 + A2/t2), A the 0.5 mm raster areas\n" % pair)
    f.write("0\tthickness-at-holes-pooled\t%.6g\tmm\tthe generator's t_holes (both rings pooled), as in the table's neck meta\n" % t_pool)

print("t_pool %.6f (meta prints %.2f)" % (t_pool, t_pool))
for h in (1, 2):
    print("hole %d: ring points %d, bore %.6f (min %.3f max %.3f), area %.2f mm2" % (h, len(rings[h - 1]), t[h - 1], t_min[h - 1], t_max[h - 1], area[h]))
print("pair bore %.6f; raster runs %d moved %.2f mm in z onto the solver's grid" % (pair, len(out_rows), shift))
