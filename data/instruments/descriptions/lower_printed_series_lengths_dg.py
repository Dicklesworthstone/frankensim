#!/usr/bin/env python3
"""Put the printed series length t_a of Lefebvre and Scavone on a DG card as a passive bore shortening.

usage: python3 -I lower_printed_series_lengths_dg.py LIVE.performance DG.performance OUT.performance

LIVE is the reed-performance card of the note (its `hole RADIUS HEIGHT BORE_RADIUS OPEN` rows give each hole's radius b, chimney
height t, bore radius a at the hole and open/closed state). DG is a frankensim-material-modal-dg-valve-v1 card whose
`dg_side_chimney PORT ELEMENT X Y RADIUS LENGTH K P` rows are the same holes in the same order (PORT = hole number + 1).

THE PHYSICS. The printed thick-wall hole of Lefebvre and Scavone (J. Acoust. Soc. Am. 131(4), 2012) is a shunt (the chimney, whose
lengths lower_printed_hole_lengths_dg.py lowers) plus a SERIES impedance Z_a = -i k Z_0 t_a acting on the bore at the hole,
with the signed series length (fs-duct hole_lengths.rs, `series_m`, with hand values in its tests)
    open   (Eq. 33):  t_a = (-0.35 + 0.06 tanh(2.7 t / b)) b d^2
    closed (Eq. 34):  t_a = (-0.12 - 0.17 tanh(2.4 t / b)) b d^2        d = b / a
Both are NEGATIVE: an inertance -rho |t_a| / (pi a^2). A DG junction has no series element and a negative lumped inertance would
not be passive. The passive equivalent is a bore SHORTENING by |t_a| at the hole: it removes the same inertance at the hole's own
area. The compliance of the removed volume goes with it (A |t_a| / (rho c^2)); that is second order in k a and is not corrected.

WHAT THE SCRIPT DOES, per hole in order from the reed, on a card whose hole faces are exact `dg_segment` boundaries (the card writes the
segment list so; the face station is the bore length minus X_M):
  * the dg_segment that ENDS at the hole's face is shortened by |t_a| (a cone keeps its end radii);
  * every `gas_station` keeps its RELATIVE position inside its segment (the profile follows the geometry; the first station stays at 0
    and the last becomes the new bore length, which the engine requires);
  * X_M of each hole is the new bore length minus its new face station (X_M is measured from the bell);
  * ELEMENT, K, P, the radii, the chimney lengths and every other row are untouched.
It REFUSES by name unless: every hole face is a segment boundary within 1e-9 m, no gas station lands within 1e-9 m of a boundary
before or after (so no span appears or disappears), the DG card has exactly as many spans as elements (so the mesh allocation cannot
change), no segment is shortened by more than half of its length, b < a at every hole, and the new bore length is the old
one minus the sum of the shortenings.

NAMED, NOT HIDDEN: G(ka) of Eq. 32 is still frozen at 1 and the mutual length of open holes is not represented (as in
lower_printed_hole_lengths_dg.py). The location of the removed inertance is the hole's face, as the printed series element is.
Writing never overwrites an existing file.
"""
import math
import sys

TOL = 1e-9  # metres


def series_length_m(b, t, a, is_open):
    """Eq. 33 (open) or Eq. 34 (closed): the signed series length referred to the bore."""
    d = b / a
    if is_open:
        return (-0.35 + 0.06 * math.tanh(2.7 * t / b)) * b * d * d
    return (-0.12 - 0.17 * math.tanh(2.4 * t / b)) * b * d * d


def refuse(message):
    sys.exit(f"lower_printed_series_lengths_dg: {message}")


def live_holes(text):
    holes = []
    for line in text.split("\n"):
        fields = line.split()
        if fields and fields[0] == "hole":
            if len(fields) != 5 or fields[4] not in ("0", "1"):
                refuse(f"malformed live hole row: {line!r}")
            holes.append((float(fields[1]), float(fields[2]), float(fields[3]), fields[4] == "1"))
    return holes


def segment_length_index(fields):
    if fields[1] == "cylinder" and len(fields) == 4:
        return 3
    if fields[1] == "cone" and len(fields) == 5:
        return 4
    refuse(f"unsupported dg_segment row: {' '.join(fields)!r}")


def lower(live_text, dg_text):
    """Return (new DG text, table rows). Refuses any mismatch between the live holes and the DG chimneys or any broken invariant."""
    holes = live_holes(live_text)
    lines = dg_text.split("\n")
    if not lines or lines[0] != "frankensim-material-modal-dg-valve-v1":
        refuse("the DG card's first line must be frankensim-material-modal-dg-valve-v1")
    seg_at = [i for i, line in enumerate(lines) if line.startswith("dg_segment ")]
    chim_at = [i for i, line in enumerate(lines) if line.startswith("dg_side_chimney ")]
    gas_at = [i for i, line in enumerate(lines) if line.startswith("gas_station ")]
    bore_at = [i for i, line in enumerate(lines) if line.startswith("dg_bore ")]
    if len(bore_at) != 1 or not seg_at or not gas_at:
        refuse("the DG card needs one dg_bore row, dg_segment rows and gas_station rows")
    if not holes or len(holes) != len(chim_at):
        refuse(f"the live card has {len(holes)} holes and the DG card {len(chim_at)} chimneys")
    bore = lines[bore_at[0]].split()
    elements = int(bore[3])
    if int(bore[1]) != len(seg_at):
        refuse(f"dg_bore says {bore[1]} segments, the card has {len(seg_at)}")
    seg_fields = [lines[i].split() for i in seg_at]
    length_index = [segment_length_index(f) for f in seg_fields]
    seg_len = [float(f[k]) for f, k in zip(seg_fields, length_index)]
    bounds = [0.0]
    for length in seg_len:
        bounds.append(bounds[-1] + length)
    total = bounds[-1]
    gas_z = [float(lines[i].split()[1]) for i in gas_at]
    if abs(gas_z[0]) > TOL or abs(gas_z[-1] - total) > TOL:
        refuse("the gas stations must run from 0 to the full bore length")
    inside = [z for z in gas_z if TOL < z < total - TOL]
    for z in inside:
        if min(abs(z - b) for b in bounds) <= TOL:
            refuse(f"a gas station at {z} m sits on a segment boundary")
    spans = len(seg_len) + len(inside)
    if spans != elements:
        refuse(f"the card has {spans} spans for {elements} elements; the mesh allocation could change, this lowering needs one element per span")
    closed_ports = {int(line.split()[1]) for line in lines if line.startswith("dg_terminal_closed ")}
    shorten = [0.0] * len(seg_len)
    rows = []
    for k, ((b, t, a, is_open), i) in enumerate(zip(holes, chim_at)):
        fields = lines[i].split(" ")
        if len(fields) != 9:
            refuse(f"chimney row {i + 1} does not have nine fields: {lines[i]!r}")
        port, x_m, radius = int(fields[1]), float(fields[3]), float(fields[5])
        if port != k + 2 or not math.isclose(radius, b, rel_tol=1e-12):
            refuse(f"hole {k + 1}: live (b {b}) does not match DG chimney port {port} (radius {radius})")
        if (port not in closed_ports) != is_open:
            refuse(f"hole {k + 1}: live card says {'open' if is_open else 'closed'}, the DG card's dg_terminal_closed list says the opposite")
        if not 0.0 < b < a:
            refuse(f"hole {k + 1}: b {b} is not inside the bore radius {a}; the printed fits are for d = b / a <= 1")
        face = total - x_m
        ends = [s for s in range(len(seg_len)) if abs(bounds[s + 1] - face) <= TOL]
        if len(ends) != 1:
            refuse(f"hole {k + 1} (port {port}): its face at {face} m is not a segment boundary")
        t_a = series_length_m(b, t, a, is_open)
        if not t_a < 0.0:
            refuse(f"hole {k + 1}: t_a {t_a} is not negative")
        shorten[ends[0]] += -t_a
        rows.append([k + 1, port, is_open, b, a, t, t_a, ends[0], x_m])
    for s, cut in enumerate(shorten):
        if cut > 0.5 * seg_len[s]:
            refuse(f"segment {s + 1} would lose {cut} of {seg_len[s]} m (more than half)")
    new_len = [length - cut for length, cut in zip(seg_len, shorten)]
    new_bounds = [0.0]
    for length in new_len:
        new_bounds.append(new_bounds[-1] + length)
    new_total = new_bounds[-1]
    if abs(new_total - (total - sum(shorten))) > 1e-12:
        refuse("the new bore length is not the old one minus the shortenings")
    # segment rows
    for i, fields, k, length in zip(seg_at, seg_fields, length_index, new_len):
        fields = list(fields)
        fields[k] = repr(length)
        lines[i] = " ".join(fields)
    # gas stations keep their relative position inside their segment
    new_gas = []
    for i, z in zip(gas_at, gas_z):
        if abs(z) <= TOL:
            z_new = 0.0
        elif abs(z - total) <= TOL:
            z_new = new_total
        else:
            s = max(s for s in range(len(seg_len)) if bounds[s] < z)
            z_new = new_bounds[s] + (z - bounds[s]) / seg_len[s] * new_len[s]
        new_gas.append(z_new)
        fields = lines[i].split(" ")
        fields[1] = repr(z_new)
        lines[i] = " ".join(fields)
    for z in new_gas[1:-1]:
        if min(abs(z - b) for b in new_bounds) <= TOL:
            refuse(f"after the shortening a gas station at {z} m sits on a segment boundary")
    # X_M of each hole: the new bore length minus its new face station
    for (k, port, is_open, b, a, t, t_a, end, x_m), i in zip(rows, chim_at):
        fields = lines[i].split(" ")
        fields[3] = repr(new_total - new_bounds[end + 1])
        lines[i] = " ".join(fields)
    return "\n".join(lines), rows, total, new_total


def main(argv):
    if len(argv) != 4:
        sys.exit(__doc__)
    live_path, dg_path, out_path = argv[1:]
    with open(live_path, encoding="utf-8") as handle:
        live_text = handle.read()
    with open(dg_path, encoding="utf-8") as handle:
        dg_text = handle.read()
    text, rows, total, new_total = lower(live_text, dg_text)
    with open(out_path, "x", encoding="utf-8") as handle:  # never overwrite
        handle.write(text)
    print(f"wrote {out_path}: {len(rows)} holes shortened the bore from {total!r} to {new_total!r} m (sum {total - new_total!r}); G(ka) frozen at 1")
    print("# hole port state b_m a_m t_m t_a_m (the bore is shortened by -t_a at the hole's face)")
    for k, port, is_open, b, a, t, t_a, end, x_m in rows:
        print(f"{k} {port} {'open' if is_open else 'closed'} {b!r} {a!r} {t!r} {t_a!r}")


if __name__ == "__main__":
    main(sys.argv)
