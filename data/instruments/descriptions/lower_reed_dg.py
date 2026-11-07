#!/usr/bin/env python3
"""Put the reed of a wind description on a DG bore.

usage: python3 -I lower_reed_dg.py DESCRIPTION.tsv TEMPLATE.performance OUT.performance

TEMPLATE is a frankensim-material-modal-dg-valve-v1 card. Its bore, gas, radiation, observation, limits and
schedule rows are kept byte for byte. This script replaces exactly three things:

  * the `valve` row with the description's rest opening;
  * the `contact` row with the lay that lower_reed_performance.py's reed-performance card gives a massive reed
    (render.rs: stiffness 1e7 x width, exponent 2, loss chi = c / (k_lay H^2) with c = 2 zeta sqrt(k m));
  * the plate block (`loss` or `solid` up to `compile_limits`) with ONE typed `reed` record whose six numbers
    are lower_reed_performance.reed_record(desc) and whose source token names the description rows.

The reed is therefore lowered in one place: the same six numbers, the same lay, the DG bore the only change.
`slit_length` is set to none, which the reed record requires (the live recipe has no viscous slit).

The description rows that sourced the reed are printed to stdout verbatim (key, value, unit, status, source), so a
page can carry the table's own basis words as they stand.
"""
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script directory
import lower_reed_performance as lrp  # noqa: E402

PLAIN_KEYS = ("reed.rest_opening_m", "reed.width_m", "reed.closing_pressure_pa",
              "reed.mass_kg", "reed.stiffness_n_m", "reed.damping_ratio")
PER_AREA_KEYS = ("reed.rest_opening_m", "reed.width_m", "reed.effective_area_m2",
                 "reed.stiffness_per_area_pa_m", "reed.mass_per_area_kg_m2", "reed.damping_rate_1_s")
MAX_TOKEN_BYTES = 256


def description_rows(path):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    return {row[0]: row for row in rows[2:] if row}


def main(argv):
    if len(argv) != 4:
        sys.exit(__doc__)
    desc_path, template_path, out_path = argv[1:]
    desc = lrp.read_description(desc_path)
    rows = description_rows(desc_path)
    fields = lrp.reed_record(desc).split()
    if len(fields) != 6:
        lrp.refuse(f"reed_record gave {len(fields)} fields, not the six of a reed row")
    height, width, _closing, mass, stiffness, zeta = fields
    h0, w, m, k, z = (float(x) for x in (height, width, mass, stiffness, zeta))
    k_lay = 1.0e7 * w
    chi = 2.0 * z * math.sqrt(max(k * m, 0.0)) / max(k_lay * h0 * h0, 1.0e-18)
    keys = PER_AREA_KEYS if "reed.stiffness_per_area_pa_m" in desc else PLAIN_KEYS
    if "reed.tip_opening_cap_m" in desc and "reed.stiffness_per_area_pa_m" in desc:
        keys += ("reed.tip_opening_cap_m",)
    for key in keys:
        if key not in rows:
            lrp.refuse(f"description has no {key} row")
    description_id = lrp.need(desc, "description_id")
    token = f"{os.path.relpath(desc_path)}#{description_id}/{','.join(keys)}"
    if len(token.encode()) > MAX_TOKEN_BYTES or not all(c.isascii() and c.isprintable() and c != " " for c in token):
        lrp.refuse(f"source token {token!r} is not a printable token of at most {MAX_TOKEN_BYTES} bytes")
    reed = f"reed {' '.join(fields)} {token}"
    contact = f"contact {k_lay!r} 2 {chi!r}"

    with open(template_path, encoding="utf-8") as handle:
        lines = handle.read().split("\n")
    if lines[0] != "frankensim-material-modal-dg-valve-v1":
        lrp.refuse(f"{template_path}: first line must be frankensim-material-modal-dg-valve-v1")

    def only(prefix):
        found = [i for i, line in enumerate(lines) if line == prefix or line.startswith(prefix + " ")]
        if len(found) != 1:
            lrp.refuse(f"{template_path}: expected exactly one `{prefix}` row, found {len(found)}")
        return found[0]

    valve, contact_at, limits_at = only("valve"), only("contact"), only("limits")
    slit, compile_at = only("slit_length"), only("compile_limits")
    if not (valve < limits_at < contact_at < compile_at):
        lrp.refuse(f"{template_path}: valve, limits, contact and compile_limits are out of order")
    plate = lines[contact_at + 1].split(" ")[0]
    if plate not in ("loss", "solid"):
        lrp.refuse(f"{template_path}: the plate block must start with `loss` or `solid`, not `{plate}`")
    lines[slit] = "slit_length none"
    lines[valve] = f"valve {height}"
    lines[contact_at] = contact
    out = lines[:contact_at + 1] + [reed] + lines[compile_at:]
    with open(out_path, "x", encoding="utf-8") as handle:  # never overwrite
        handle.write("\n".join(out))
    print(f"wrote {out_path}: {reed}\n{contact}")
    print("# description rows (key, value, unit, status, source), verbatim:")
    for key in keys:
        print("\t".join(rows[key]))


if __name__ == "__main__":
    main(sys.argv)
