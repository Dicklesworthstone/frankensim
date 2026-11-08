"""Write a NEW '<parent>-declared-lowband' matdb pack SOURCE directory (manifest.tsv + properties.tsv) from a parent shell pack source directory, changing ONLY
  - the lower edge of every `validity <claim> omega LO HI rad/s` line (to 2 pi x LO_HZ, printed with 9 significant digits like the parent's edges),
  - the observation notes (the declaration is appended to the notes field of every observation, i.e. to the source text of every claim),
  - the manifest's pack_id, citation and redistribution lines (the variant names its parent and says what differs).
Every scalar, uncertainty and temperature-validity line is byte-identical to the parent's; the script verifies that line by line and refuses otherwise.
It never overwrites (files open in mode x) and never edits the parent.

This is the generator behind adb4789a1, written for it on 11-03 and not committed then (bf0604080 first committed a maple-only script under this name, which
this file replaces in place). The two variants of the double-bass shell body (adb4789a1) were written by this script from the parents violin-spruce-norway-viala2020-shell (LO_HZ 90) and
violin-maple-sycamore-sonderegger2013-shell (LO_HZ 90): run on the parent seeds as committed before 825ec4b44 it rebuilds both committed variant seeds byte for byte.
Since 825ec4b44 the maple parent takes eta_R from Alkadri et al. 2018's plate C (radial first-mode scan, free-free bending at 80-500 Hz) in place of the eta_LR
surrogate; for a parent that carries it (`loss_factor_radial_plate_c` in its loss observation) the maple declaration names that scan beside Laemmlein's sweep, and
refuses a lower edge below 80 Hz, the bottom of the scan. For every other parent the declaration is the earlier text byte for byte.

usage: python3 -I gen_declared_lowband_pack.py PARENT_DIR OUT_DIR KIND LO_HZ        (KIND: spruce | maple)
"""
import math
import sys
from pathlib import Path

parent, out, kind, lo_hz = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], float(sys.argv[4])
SPEC = {
    "spruce": dict(
        parent_id="violin-spruce-norway-viala2020-shell", parent_lo_hz=100.0, parent_hi_hz=5000.0, parent_band="100-5000 Hz",
        stiffness=("stiffness rows are Viala's audio-band (100-5000 Hz) model-updating values with the Wood Handbook nu_LR, "
                   "used unchanged below 100 Hz"),
        loss="loss rows extrapolate the parent's constant loss factor below its source band (Viala's audio-band model updating, 100-5000 Hz)",
        edge="100 Hz"),
    "maple": dict(
        parent_id="violin-maple-sycamore-sonderegger2013-shell", parent_lo_hz=150.0, parent_hi_hz=10000.0, parent_band="150-10000 Hz",
        stiffness=("stiffness rows are static (E_L, E_R) and 1 MHz ultrasonic (G_LT, G_RT) measurements and Laemmlein's audio-band G_LR, "
                   "used unchanged below 150 Hz"),
        loss="loss rows extrapolate the parent's constant loss factor below its measured sweep (Laemmlein's loudspeaker sweep, 150-10000 Hz)",
        edge="150 Hz"),
}[kind]
PLATE_C_SCAN_LO_HZ = 80.0  # Alkadri et al. 2018, radial first-mode scans 80-500 Hz
new_id = SPEC["parent_id"] + "-declared-lowband"
assert out.name == new_id, f"OUT_DIR must be named {new_id}"
assert 0.0 < lo_hz < SPEC["parent_lo_hz"], "the new lower edge must lie below the parent's"
fmt = lambda x: f"{x:.9g}"
parent_lo, parent_hi = fmt(2 * math.pi * SPEC["parent_lo_hz"]), fmt(2 * math.pi * SPEC["parent_hi_hz"])
new_lo = fmt(2 * math.pi * lo_hz)

pm = (parent / "manifest.tsv").read_text().split("\n")
pp = (parent / "properties.tsv").read_text().split("\n")
if kind == "maple" and any(l.startswith("observation\t") and "loss_factor_radial_plate_c" in l for l in pp):
    assert lo_hz >= PLATE_C_SCAN_LO_HZ, "the parent's eta_R comes from a scan that starts at 80 Hz; a lower edge would be an extrapolation the text does not say"
    SPEC["loss"] = ("loss rows extrapolate the parent's constant loss factors below their measured sweeps (eta_L and eta_LR: Laemmlein's loudspeaker sweep, 150-10000 Hz; "
                    f"eta_R: the radial first-mode scan of Alkadri's plate C, free-free bending at 80-500 Hz, which holds the lower edge {lo_hz:g} Hz and everything above it)")
declaration = (
    f" DECLARED LOW BAND: the band below {SPEC['edge']} (this pack's lower edge {lo_hz:g} Hz = omega {new_lo} rad/s) is DECLARED BY THE PACK AUTHOR, "
    f"NOT MEASURED. Parent pack {SPEC['parent_id']} (band {SPEC['parent_band']}); every value is copied byte for byte from it and only the band's lower edge moved. "
    f"{SPEC['stiffness'][0].upper() + SPEC['stiffness'][1:]}; {SPEC['loss']}. "
    f"No measurement of this wood below the parent's band is claimed or used.")

nm, np_ = [], []
for line in pm:
    f = line.split("\t")
    if f[0] == "pack_id":
        assert f[1] == SPEC["parent_id"], f"parent pack_id is {f[1]}"
        f[1] = new_id
    elif f[0] == "citation":
        f[1] += (f"; variant {new_id} of {SPEC['parent_id']}: band's lower edge declared by the pack author, not measured")
    elif f[0] == "redistribution":
        assert len(f) == 3 and f[1] == "permitted", "redistribution grammar: permitted<TAB>text"
        f[2] += (f"; this variant changes only the lower edge of the angular-frequency band and the observation notes of {SPEC['parent_id']}, "
                 f"whose own terms apply")
    nm.append("\t".join(f))
edges = set()
for line in pp:
    f = line.split("\t")
    if f[0] == "validity" and len(f) == 6 and f[2] == "omega":
        assert f[5] == "rad/s", line
        edges.add((f[3], f[4]))
        f[3] = new_lo
    elif f[0] == "observation":
        assert len(f) == 5, f"observation with {len(f)} fields"
        f[4] += declaration
    np_.append("\t".join(f))
assert edges == {(parent_lo, parent_hi)}, f"the parent's omega edges are {edges}, expected {(parent_lo, parent_hi)}"

# Verify line by line against the parent: only the three allowed kinds of difference exist.
assert len(nm) == len(pm) and len(np_) == len(pp)
counts = {"manifest": 0, "omega-lower-edge": 0, "observation-notes": 0, "identical": 0}
for a, b in zip(pm, nm):
    if a != b:
        assert a.split("\t")[0] in ("pack_id", "citation", "redistribution"), a[:60]
        counts["manifest"] += 1
for a, b in zip(pp, np_):
    fa, fb = a.split("\t"), b.split("\t")
    if a == b:
        counts["identical"] += 1
    elif fa[0] == "validity":
        assert fa[:3] == fb[:3] and fa[4:] == fb[4:] and fa[2] == "omega" and fb[3] == new_lo, a
        counts["omega-lower-edge"] += 1
    elif fa[0] == "observation":
        assert fa[:4] == fb[:4] and fb[4] == fa[4] + declaration, a[:60]
        counts["observation-notes"] += 1
    else:
        raise AssertionError("an unexpected line differs: " + a[:80])
for a, b in zip(pp, np_):  # scalars and uncertainties are byte-identical
    if a.startswith(("scalar\t", "uncertainty\t")) or (a.startswith("validity\t") and a.split("\t")[2] == "T"):
        assert a == b
out.mkdir()
with open(out / "manifest.tsv", "x") as g:
    g.write("\n".join(nm))
with open(out / "properties.tsv", "x") as g:
    g.write("\n".join(np_))
print(f"wrote {out}: {counts}; new lower edge {lo_hz:g} Hz -> omega {new_lo} rad/s (parent {parent_lo}); {sum(1 for l in pp if l.startswith('scalar'))} scalars byte-identical")
