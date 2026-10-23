"""make_results.py OUT_MODES OUT_OUTCOME RUN_ID=FILE ... : the air-coupled mode lists (neck open) of the named run outputs as one table (body-mode-results.tsv).
OUT_OUTCOME is accepted and not written: the outcome table (body-mode-outcome.tsv) is written by make_outcome.py from the classifier output. Standard library only."""
import sys

out_modes, out_outcome = sys.argv[1:3]
runs = [a.split("=", 1) for a in sys.argv[3:]]
SUP = "SUPERSEDED, a run of the cracked model (456 of 1239 zone-interface nodes were not welded; repaired in 9080c6a69, found by audit_welds.py): "
# Every run below used the maple pack as it was before 4f9abefe1 (density 563 kg/m3, the oven-dry value read by mistake; 626 kg/m3 from then on).
MAPLE = "; maple pack as before 4f9abefe1 (density 563 kg/m3, the oven-dry value; the pack carries 626 kg/m3 from then on)"
REP = "repaired model (shell-assembly-h8.tsv of 9080c6a69, 48 components, 7042 nodes, 119 ties; cavity record 1.9653 L from the facet sum unless the run name starts with vex; reproduce with gen_violin_shell_assembly.py SPECIMEN_DIR OUT.tsv 8 for the facet record and ... 8 cavity_volume=exact for the corrected one), binary of 8a0fb7ec1 on trj, packs spruce (guitar-spruce-sitka-fpl-gtr282-shell) and maple (violin-maple-sycamore-sonderegger2013-shell): "
NOTES = {
    "cracked-base": SUP + "shell-assembly-h8.tsv of 8a0dc81b3, 28 shell modes between 20 and 700 Hz",
    "cracked-t018": SUP + "every plate zone thickness 0.18 mm thinner, 26 shell modes",
    "cracked-sc102": SUP + "base table scaled by 1.02 in plan by scale_table.py, 26 shell modes",
    "cracked-sc098": SUP + "base table scaled by 0.98 in plan by scale_table.py, 26 shell modes",
    "base-n5": REP + "the registered run, 5 shell modes between 20 and 900 Hz (the in-vacuo list is identical for 3, 4 and 5 modes and for window ends 700 and 900 Hz)",
    "base-n4": REP + "4 shell modes, 20 to 900 Hz: truncation sensitivity of the coupled list",
    "base-n3": REP + "3 shell modes, 20 to 700 Hz: truncation sensitivity of the coupled list",
    "base-n6": REP + "6 shell modes, 20 to 900 Hz: truncation sensitivity of the coupled list",
    "wh10": "repaired model re-meshed with a target edge of 10 mm instead of 8 (38 components, 4924 nodes, 6988 triangles, 91 ties; same generator, same packs, same binary), 5 shell modes: mesh ladder",
    "wh12": "repaired model re-meshed with a target edge of 12 mm (32 components, 3644 nodes, 5185 triangles, 74 ties), 5 shell modes: mesh ladder",
    "alpha-0.1": "repaired model, drilling penalty 1e-4 (x0.1 of the default 1e-3) by --drilling-alpha, binary of b85323bb5, 5 shell modes: declared sensitivity of the fully tied T-junctions",
    "alpha-10": "repaired model, drilling penalty 1e-2 (x10 of the default) by --drilling-alpha, binary of b85323bb5, 5 shell modes: declared sensitivity of the fully tied T-junctions",
    "fplshear": "repaired model with the maple shear moduli as a declared sensitivity: G_LR 1210 and G_LT 687 MPa (Wood Handbook sugar maple ratios 0.111 and 0.063 times the same E_L) instead of the ultrasonic 1680 and 1070 MPa, G_RT unchanged (pack from gen_violin_maple_pack.py shear=fpl), 5 shell modes",
    "npost": "repaired model WITHOUT the soundpost (generator option post=0: 46 components, 6778 nodes, 9792 triangles, 113 ties), 5 shell modes, chosen AFTER the first result because a measured body-only violin (no neck, no soundpost; Woodhouse 2014 Fig. 18, Stoppani, A0 243, CBR 394, B1- 403, B1+ 562 Hz, a different violin) exists; not a registered test",
    "vex-n5": REP + "the same geometry with the cavity record corrected: 2.0610 L (exact inner-skin volume minus the bass bar and soundpost shells) instead of the origin-dependent facet sum 1.9653 L (P13 mail 5583; generator option cavity_volume=exact), 5 shell modes: the corrected A0 run",
    "vex-n6": REP + "as vex-n5 with 6 shell modes (truncation sensitivity)",
    "tw018": REP + "declared sensitivity: every plate zone thickness 0.18 mm thinner (the largest volume-against-mesh difference), 5 shell modes",
    "tw004": REP + "declared sensitivity: every plate zone thickness 0.04 mm thinner (the smallest volume-against-mesh difference), 5 shell modes",
    "scw098": REP + "declared sensitivity: plan view scaled by 0.98 about the body centre by scale_table.py (topology, heights, ties unchanged), 5 shell modes",
    "scw102": REP + "declared sensitivity: plan view scaled by 1.02 about the body centre by scale_table.py, 5 shell modes",
}


def coupled(path):
    t = open(path, encoding="utf-8").read().splitlines()
    s = next(i for i, l in enumerate(t) if l.startswith("open sound hole"))
    rows = []
    for l in t[s + 2:]:
        p = l.split()
        if not p or not p[0].isdigit():
            break
        rows.append(p)
    return rows


lines = ["frankensim.specimen-body-mode-results.v1", "run\trank\tf_hz\tneck_share\tnet_cavity_dv\ttop_mean_w\tback_mean_w\tq_solid_kv_not_a_result\tq_rad\tnote"]
for run, path in runs:
    for i, p in enumerate(coupled(path)):
        lines.append("\t".join([run, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], NOTES.get(run, "run " + run) + MAPLE if i == 0 else "same run"]))
open(out_modes, "w", encoding="utf-8").write("\n".join(lines) + "\n")
print("wrote", out_modes, len(lines) - 2, "rows")
