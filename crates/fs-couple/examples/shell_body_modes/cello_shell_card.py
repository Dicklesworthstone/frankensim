"""A cello shell card with ONLY its two body-derived records rebuilt from a table's `material_strings assemble-shell`
output: the 'body' record (box length, width, reference thickness) and the 'cavity physical-neck' record (cavity volume,
neck radius, neck length = the f-holes' equivalent neck). Every other byte is the input card's. Refuses to overwrite.
On a card whose records already are the table's, the output equals the input (check_cello_shell.py, check E).
usage: python3 -I cello_shell_card.py CARD ASSEMBLE_SHELL_OUTPUT OUT_CARD"""
import sys
from pathlib import Path

card, log, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
KEYS = ["body-length-m", "body-width-m", "reference-thickness-m", "cavity-volume-m3", "cavity-neck-radius-m", "cavity-neck-length-m"]
val = {}
for l in log.read_text().split("\n"):
    f = l.split(" ")
    if len(f) == 2 and f[0] in KEYS:
        assert f[0] not in val, ("repeated", f[0])
        float(f[1])
        val[f[0]] = f[1]
assert sorted(val) == sorted(KEYS), ("missing", sorted(set(KEYS) - set(val)))
lines = card.read_text().split("\n")
b = [i for i, l in enumerate(lines) if l.startswith("body ")]
c = [i for i, l in enumerate(lines) if l.startswith("cavity ")]
assert len(b) == 1 and len(c) == 1, (b, c)
fb, fc = lines[b[0]].split(" "), lines[c[0]].split(" ")
assert len(fb) == 7 and fb[1] == "violin-spruce-norway-viala2020-shell.fsmatpk" and fb[5:] == ["256", "1"], fb
assert len(fc) == 8 and fc[1] == "physical-neck" and fc[5:] == ["flanged", "flanged", "1"], fc
new = list(lines)
new[b[0]] = " ".join(fb[:2] + [val["body-length-m"], val["body-width-m"], val["reference-thickness-m"]] + fb[5:])
new[c[0]] = " ".join(fc[:2] + [val["cavity-volume-m3"], val["cavity-neck-radius-m"], val["cavity-neck-length-m"]] + fc[5:])
changed = [i for i, (x, y) in enumerate(zip(lines, new)) if x != y]
assert len(new) == len(lines) and set(changed) <= {b[0], c[0]}, changed
with open(out, "x") as f:
    f.write("\n".join(new))
print("card %s -> %s: changed lines %s" % (card.name, out.name, [i + 1 for i in changed]))
for i in changed:
    print("  old %d: %s\n  new %d: %s" % (i + 1, lines[i], i + 1, new[i]))
