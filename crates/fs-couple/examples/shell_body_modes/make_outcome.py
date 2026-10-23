"""make_outcome.py OUT.tsv RUN=CLASSIFY_OUTPUT.txt ... : the registered outcome of the identification rules per run (body-mode-outcome.tsv) from the printed output of classify_violin_modes.py.
A run name starting with `cracked` is labelled superseded. Standard library only."""
import re
import sys

out = sys.argv[1]
BANDS = {"a0": (270, 290), "cbr": (380, 420), "b1-": (430, 480), "b1+": (510, 580)}
KEY = {"a0": "a0", "cbr": "cbr", "b1-": "b1m", "b1+": "b1p"}
EXTRA = {
    ("base-n5", "cbr"): "the registered text says 'concentrated at the waist' without a number; the classifier used a waist share of 0.3; the outcome is the same for any threshold up to 0.58 (the mode has 0.585)",
    ("base-n5", "b1-"): "ranks 4 to 6 above the CBR mode (522.17, 621.58, 648.83 Hz) all have top and back displacement of the same sign along the outward directors, the rule needs opposite signs; the registered B1- definition (opposite signs along the outward directors, net volume change below 20 percent of A0, above the CBR mode) CONTRADICTS the page named as its source, which says (read after the run) that B1- and B1+ both have a significant volume change with the plates moving in opposite directions: the rule was mis-specified, the outcome stays as registered; by that description the only mode with a large volume change besides A0 is rank 4 at 522.17 Hz and none lies in the B1- band",
    ("base-n5", "b1+"): "the chain stops at b1-; for information only, outside the chain: rank 4 at 522.17 Hz is even in x (odd fraction 0.30), has top and back of the same sign and the largest net cavity volume change of the list, and lies inside the band",
}
EXTRA[("vex-n5", "a0")] = "the corrected cavity record (2.0610 L instead of the origin-dependent facet sum 1.9653 L): this A0 supersedes the base-n5 value, the geometry and every other record of the table are the same"
EXTRA[("vex-n5", "cbr")] = "unchanged by the cavity record (it does not couple to the air): 480.73 against 480.75 Hz"
for _i in ("a0", "cbr", "b1-", "b1+"):
    EXTRA[("npost", _i)] = "the registered rules applied to a configuration that was NOT registered (the body without the soundpost, chosen after the first result because a measured body-only violin exists); information only"
EXTRA[("corrected-n6", "b1-")] = ("the mode the registered rule finds carries 88 percent of its kinetic energy in the bridge (lumped nodal masses from the --node-shapes "
                               "output): a bridge mode, not a corpus mode; the rule reads plate displacement only; the registered outcome stands, information only")
EXTRA[("corrected-n6", "b1+")] = "no mode of the list above the b1- row (rank 6) qualifies; for information only, outside the chain: rank 3 at 530.72 Hz is the breathing mode, inside the band"
EXTRA[("alt-laemmlein-n6", "b1-")] = ("rank 6 (700.12 Hz) has the B1- rule's shape and 0.787 of its kinetic energy in the bridge: left out of the corpus chain as a "
                                   "bridge mode (e4b27a28e); the classifier before that commit would have named it B1-")
EXTRA[("c5-maple-shear-n6", "b1-")] = ("rank 6 (691.26 Hz) has the B1- rule's shape and 0.883 of its kinetic energy in the bridge: left out of the corpus chain as a "
                                    "bridge mode (e4b27a28e); the classifier before that commit names it B1- (MISS, +211.3 Hz)")
rows = ["frankensim.specimen-body-mode-outcome.v1", "run\tid\tband_lo_hz\tband_hi_hz\tregistered_rule_outcome\tcandidate_rank\tcandidate_f_hz\tsigned_offset_hz\tnote"]
for a in sys.argv[2:]:
    run, path = a.split("=", 1)
    text = open(path, encoding="utf-8").read().splitlines()
    found = {}
    for l in text:
        m = re.match(r"^(a0|cbr|b1-|b1\+)\s*:\s*(.*)$", l)
        if m:
            found[m.group(1)] = m.group(2)
    scored = {}
    for l in text:
        m = re.match(r"^(a0|cbr|b1m|b1p)\s+(.*)$", l)
        if m:
            scored[m.group(1)] = m.group(2)
    prefix = "SUPERSEDED (cracked model). " if run.startswith("cracked") else ""
    for ident, (lo, hi) in BANDS.items():
        f = found.get(ident, "not found")
        s = scored.get(KEY[ident], "")
        rk = re.search(r"rank (\d+) at ([0-9.]+) Hz", f)
        off = re.search(r"signed offset (-?[0-9.]+) Hz", s)
        if rk:
            outcome, rank, freq = "found, %s" % ("MET" if "MET" in s else "MISS"), rk.group(1), rk.group(2)
            offset = off.group(1) if off else "0.0"
        else:
            outcome, rank, freq, offset = "not found", "none", "none", "none"
        rows.append("\t".join([run, ident, str(lo), str(hi), outcome, rank, freq, offset, prefix + f + ("; " + EXTRA[(run, ident)] if (run, ident) in EXTRA else "")]))
open(out, "w", encoding="utf-8").write("\n".join(rows) + "\n")
print("wrote", out, len(rows) - 2, "rows")
