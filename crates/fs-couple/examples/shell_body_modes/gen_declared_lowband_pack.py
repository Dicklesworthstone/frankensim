"""Write the seed of the declared low-band sycamore maple shell pack of the double-bass body from the seed of its parent pack.

The double-bass shell body has four assembled modes (94.8, 113.9, 123.0 and 126.0 Hz in vacuo) below the 150 Hz lower edge of violin-maple-sycamore-sonderegger2013-shell, which
refuses them. The variant violin-maple-sycamore-sonderegger2013-shell-declared-lowband (adb4789a1) copies every value of the parent byte for byte and moves only the lower edge
of the angular-frequency band to 90 Hz (omega 565.486678 rad/s); that edge is DECLARED BY THE PACK AUTHOR, NOT MEASURED, and the two observation notes and the manifest say so.
The first seed was derived by hand. This script is that derivation, so that a change of the parent (825ec4b44 took Alkadri et al. 2018's plate C radial loss, eta_R 0.0240, in
place of the eta_LR surrogate 0.0162) reaches the variant by the same steps and not by a second hand edit.

  properties.tsv: every `validity ... omega 942.477796 62831.8531 rad/s` row takes the lower edge 565.486678 (all of them are asserted to be that row); the two observation
                  rows (sond-shell, laemmlein-maple-loss) gain the DECLARED LOW BAND sentence in their last field. Nothing else changes.
  manifest.tsv:   pack_id gains `-declared-lowband`; the redistribution text and the citation each gain one clause naming the parent.

The sentence names the loss rows' measured sweeps. When the parent's loss observation carries plate C (`loss_factor_radial_plate_c`), eta_R's sweep is Alkadri's radial
first-mode scan (free-free bending, 80-500 Hz), not Laemmlein's 150-10000 Hz, and the sentence says so; for the earlier parent it is the earlier text byte for byte.

usage: python3 -I gen_declared_lowband_pack.py PARENT_SEED_DIR OUT_SEED_DIR     (OUT_SEED_DIR must not exist; the parent's directory name must be the maple shell pack's)
Check used for the earlier seed: run on the parent seed as committed before 825ec4b44, the output equals the committed declared-lowband seed byte for byte (both files).
"""
import math
import os
import sys

PARENT_ID = "violin-maple-sycamore-sonderegger2013-shell"
LOW_ID = PARENT_ID + "-declared-lowband"
OMEGA_PARENT, OMEGA_LOW = "942.477796", "565.486678"
OMEGA_TOP = "62831.8531"
assert abs(float(OMEGA_PARENT) - 2.0 * math.pi * 150.0) < 1e-5 and abs(float(OMEGA_LOW) - 2.0 * math.pi * 90.0) < 1e-5

LOSS_SURROGATE = "loss rows extrapolate the parent's constant loss factor below its measured sweep (Laemmlein's loudspeaker sweep, 150-10000 Hz)."
LOSS_PLATE_C = ("loss rows extrapolate the parent's constant loss factors below their measured sweeps (eta_L and eta_LR: Laemmlein's loudspeaker sweep, 150-10000 Hz; "
                "eta_R: the radial first-mode scan of Alkadri's plate C, free-free bending at 80-500 Hz, which holds the lower edge 90 Hz and everything above it).")
NOTE = (" DECLARED LOW BAND: the band below 150 Hz (this pack's lower edge 90 Hz = omega " + OMEGA_LOW + " rad/s) is DECLARED BY THE PACK AUTHOR, NOT MEASURED. Parent pack "
        + PARENT_ID + " (band 150-10000 Hz); every value is copied byte for byte from it and only the band's lower edge moved. Stiffness rows are static (E_L, E_R) and 1 MHz "
        "ultrasonic (G_LT, G_RT) measurements and Laemmlein's audio-band G_LR, used unchanged below 150 Hz; @LOSS@ No measurement of this wood below the parent's band is claimed or used.")
REDIS_ADD = "; this variant changes only the lower edge of the angular-frequency band and the observation notes of " + PARENT_ID + ", whose own terms apply"
CITE_ADD = "; variant " + LOW_ID + " of " + PARENT_ID + ": band's lower edge declared by the pack author, not measured"


def properties(text):
    lines = text.split("\n")
    assert lines[-1] == "" and lines[0] == "frankensim.matdb-source.v1"
    loss = [l.split("\t") for l in lines if l.startswith("observation\tlaemmlein-maple-loss\t")]
    assert len(loss) == 1
    plate_c = "loss_factor_radial_plate_c" in loss[0][3]
    note = NOTE.replace("@LOSS@", LOSS_PLATE_C if plate_c else LOSS_SURROGATE)
    out, n_omega, n_scalar, n_obs = [], 0, 0, 0
    for l in lines:
        f = l.split("\t")
        if f[0] == "scalar":
            n_scalar += 1
        if f[0] == "validity" and f[2] == "omega":
            assert f[3:] == [OMEGA_PARENT, OMEGA_TOP, "rad/s"], l
            f[3] = OMEGA_LOW
            n_omega += 1
        elif f[0] == "observation" and f[1] in ("sond-shell", "laemmlein-maple-loss"):
            assert len(f) == 5 and not f[4].endswith(note), l[:80]
            f[4] += note
            n_obs += 1
        out.append("\t".join(f))
    assert n_omega == n_scalar and n_omega > 0 and n_obs == 2 and OMEGA_PARENT not in "\n".join(out), (n_omega, n_scalar, n_obs)
    return "\n".join(out), n_omega, plate_c


def manifest(text):
    lines = text.split("\n")
    assert lines[-1] == "" and lines[0] == "frankensim.matdb-manifest.v1"
    f1, f2, f3 = (lines[i].split("\t") for i in (1, 2, 3))
    assert f1 == ["pack_id", PARENT_ID] and f2[0] == "redistribution" and len(f2) == 3 and f3[0] == "citation" and len(f3) == 2
    f1[1] = LOW_ID
    f2[2] += REDIS_ADD
    f3[1] += CITE_ADD
    lines[1], lines[2], lines[3] = "\t".join(f1), "\t".join(f2), "\t".join(f3)
    return "\n".join(lines)


def main():
    parent, out = sys.argv[1].rstrip("/"), sys.argv[2].rstrip("/")
    assert os.path.basename(parent) == PARENT_ID, "the parent directory must be named " + PARENT_ID
    props, n, plate_c = properties(open(parent + "/properties.tsv", encoding="utf-8", newline="").read())
    man = manifest(open(parent + "/manifest.tsv", encoding="utf-8", newline="").read())
    os.makedirs(out)
    for name, body in (("properties.tsv", props), ("manifest.tsv", man)):
        with open(out + "/" + name, "x", encoding="utf-8", newline="") as fh:
            fh.write(body)
    print("%s: %d validity rows moved to %s rad/s; parent carries plate C: %s" % (out, n, OMEGA_LOW, plate_c))


if __name__ == "__main__":
    main()
