#!/usr/bin/env python3
"""Put the printed thick-wall tone-hole lengths of Lefebvre and Scavone on a DG chimney card.

usage: python3 -I lower_printed_hole_lengths_dg.py LIVE.performance DG.performance OUT.performance

LIVE is the reed-performance card of the note (its `hole RADIUS HEIGHT BORE_RADIUS OPEN` rows give each hole's radius b, chimney
height t, bore radius a at the hole and open/closed state). DG is a frankensim-material-modal-dg-valve-v1 card whose
`dg_side_chimney PORT ELEMENT X Y RADIUS LENGTH K P` rows are the same holes in the same order (PORT = hole number + 1). The script
rewrites exactly the LENGTH field of each chimney and nothing else; the DG route has no `tonehole_lengths` row, so this is its
lowering of the same printed set that `fs-duct/src/hole_lengths.rs` (`HoleLengths::Printed`) gives the live route:

  * open hole: the chimney grows by  t_i + t_m + (t_r - 0.8216 b)
        t_i = b * (0.822 - 0.095 d - 1.566 d^2 + 2.138 d^3 - 1.640 d^4 + 0.502 d^5)          Eq. 31, d = b / a
        t_m = b d (1 + 0.207 d^3) / 8                                                        Eq. 8 (Nederveen et al. 1998)
        t_r = (0.822 - 0.47 [b / (a + t)]^0.8) b                                             Eq. 11 (Dalmont and Nederveen 2001)
    where 0.8216 b is the flanged-pipe end correction the radiation load already carries at low frequency (the DG route's
    computed baffled-piston radiation has it too), so only the difference is added.
  * closed hole: the sealed part of the chimney grows by t_m.

NAMED, NOT HIDDEN:
  * G(ka) of Eq. 32 is FROZEN AT 1 (its ka -> 0 value): a DG chimney is a duct of one constant length, and the printed inner
    length's frequency factor G = 1 + H(d) I(ka) cannot be a function of frequency there. On the live route G acts (ka <= 1);
    the difference is t_i (G - 1), a few percent of t_i at ka = 0.3 for a clarinet-like hole.
  * The signed SERIES length t_a (Eq. 33 open, Eq. 34 closed), an inertance acting on the BORE at the hole, is ABSENT: a DG
    junction has no series element, and this script does not fake one. The pitch share of t_a is therefore not in this card.
  * A hole with b >= a keeps the engine's terms on the live route (the fits are printed for d <= 1); here it is REFUSED by name.

Every number is printed in a table with its equation terms. Writing never overwrites an existing file.
"""
import math
import sys

FLANGED_END_CORRECTION = 0.8216  # fs-duct hole_lengths.rs FLANGED_END_CORRECTION: the flanged load keeps 0.8216 b


def inner_over_b(d):
    """Eq. 31: the inner length over b."""
    return 0.822 + d * (-0.095 + d * (-1.566 + d * (2.138 + d * (-1.640 + d * 0.502))))


def matching_m(b, d):
    """Eq. 8: the matching-volume length."""
    return b * d * (1.0 + 0.207 * d * d * d) / 8.0


def radiation_m(b, a, t):
    """Eq. 11: the radiation length of a hole drilled through a thick pipe of outer radius a + t."""
    return (0.822 - 0.47 * (b / (a + t)) ** 0.8) * b


def open_extra_m(b, a, t):
    """Length added to an open chimney's height, G(ka) frozen at 1."""
    d = b / a
    return inner_over_b(d) * b + matching_m(b, d) + (radiation_m(b, a, t) - FLANGED_END_CORRECTION * b)


def sealed_extra_m(b, a):
    """Length added to a closed chimney's height."""
    return matching_m(b, b / a)


def refuse(message):
    sys.exit(f"lower_printed_hole_lengths_dg: {message}")


def live_holes(text):
    holes = []
    for line in text.split("\n"):
        fields = line.split()
        if fields and fields[0] == "hole":
            if len(fields) != 5 or fields[4] not in ("0", "1"):
                refuse(f"malformed live hole row: {line!r}")
            holes.append((float(fields[1]), float(fields[2]), float(fields[3]), fields[4] == "1"))
    return holes


def lower(live_text, dg_text):
    """Return (new DG text, table rows). Refuses any mismatch between the live holes and the DG chimneys."""
    holes = live_holes(live_text)
    lines = dg_text.split("\n")
    if not lines or lines[0] != "frankensim-material-modal-dg-valve-v1":
        refuse("the DG card's first line must be frankensim-material-modal-dg-valve-v1")
    at = [i for i, line in enumerate(lines) if line.startswith("dg_side_chimney ")]
    if not holes or len(holes) != len(at):
        refuse(f"the live card has {len(holes)} holes and the DG card {len(at)} chimneys")
    closed_ports = {int(line.split()[1]) for line in lines if line.startswith("dg_terminal_closed ")}
    rows = []
    for k, ((b, t, a, is_open), i) in enumerate(zip(holes, at)):
        fields = lines[i].split(" ")
        if len(fields) != 9:
            refuse(f"chimney row {i + 1} does not have nine fields: {lines[i]!r}")
        port, radius, length = int(fields[1]), float(fields[5]), float(fields[6])
        if port != k + 2 or not math.isclose(radius, b, rel_tol=1e-12) or not math.isclose(length, t, rel_tol=1e-12):
            refuse(f"hole {k + 1}: live (b {b}, t {t}) does not match DG chimney port {port} (radius {radius}, length {length})")
        if (port not in closed_ports) != is_open:
            refuse(f"hole {k + 1}: live card says {'open' if is_open else 'closed'}, the DG card's dg_terminal_closed list says the opposite")
        if not 0.0 < b < a:
            refuse(f"hole {k + 1}: b {b} is not inside the bore radius {a}; the printed fits are for d = b / a <= 1")
        extra = open_extra_m(b, a, t) if is_open else sealed_extra_m(b, a)
        fields[6] = repr(t + extra)
        lines[i] = " ".join(fields)
        rows.append((k + 1, port, is_open, b, a, t, extra, t + extra))
    return "\n".join(lines), rows


def main(argv):
    if len(argv) != 4:
        sys.exit(__doc__)
    live_path, dg_path, out_path = argv[1:]
    with open(live_path, encoding="utf-8") as handle:
        live_text = handle.read()
    with open(dg_path, encoding="utf-8") as handle:
        dg_text = handle.read()
    text, rows = lower(live_text, dg_text)
    with open(out_path, "x", encoding="utf-8") as handle:  # never overwrite
        handle.write(text)
    print(f"wrote {out_path}: {len(rows)} chimney lengths; G(ka) frozen at 1; series length t_a absent")
    print("# hole port state b_m a_m t_m extra_m new_length_m  (extra = t_i + t_m + t_r - 0.8216 b when open, t_m when closed)")
    for k, port, is_open, b, a, t, extra, new in rows:
        print(f"{k} {port} {'open' if is_open else 'closed'} {b!r} {a!r} {t!r} {extra!r} {new!r}")


if __name__ == "__main__":
    main(sys.argv)
