"""make_results.py OUT_MODES OUT_OUTCOME RUN_ID=FILE ... : the air-coupled mode lists (neck open) of the named run outputs as one table (body-mode-results.tsv), and the outcome of the
pre-registered identification of the base run (body-mode-outcome.tsv). Standard library only."""
import sys

out_modes, out_outcome = sys.argv[1:3]
runs = [a.split("=", 1) for a in sys.argv[3:]]
NOTES = {
    "base": "the registered model: shell-assembly-h8.tsv, 28 shell modes between 20 and 700 Hz, packs spruce (guitar-spruce-sitka-fpl-gtr282-shell) and maple (violin-maple-sycamore-sonderegger2013-shell)",
    "t018": "declared sensitivity: every plate zone thickness 0.18 mm thinner (the largest of the volume-against-mesh differences), 26 shell modes",
    "sc102": "declared sensitivity: the base table scaled by 1.02 in plan about the body centre by scale_table.py (topology, heights, ties unchanged), 26 shell modes",
    "sc098": "declared sensitivity: the base table scaled by 0.98 in plan about the body centre by scale_table.py (topology, heights, ties unchanged), 26 shell modes",
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
        lines.append("\t".join([run, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], NOTES.get(run, "run") if i == 0 else "same run"]))
open(out_modes, "w", encoding="utf-8").write("\n".join(lines) + "\n")
print("wrote", out_modes, len(lines) - 2, "rows")
