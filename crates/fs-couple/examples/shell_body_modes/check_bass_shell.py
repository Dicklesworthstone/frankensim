"""Reproduction checks for data/instruments/geometry/double-bass-johannsson-shell, the body of the posted double-bass notes E1, A1, D2 and G2.

usage: python3 -I check_bass_shell.py [--repo REPO_ROOT] [--bin PATH_TO_material_strings] [--generator]

  C  specimen (numpy only; any machine): bass_specimen.py is run on the committed il Cannone tables and must write the 14 files listed in specimen-files.md5, and its PROVENANCE.tsv
     must equal the committed one.
  T  table and cards (standard library only): the table's counts equal its 'totals' row and the committed md5 below; every card names the declared-lowband spruce pack and carries the
     table's cavity row and one common body record; every assembly names the declared-lowband maple pack and one bridge point.
  E  profiles (needs --bin, the R5C material_strings or later): `assemble-shell` on the table and each assembly file must write the profile whose md5 is in profiles.md5, and print the
     body length, width and cavity numbers the cards carry.
  B  table from the specimen (needs --generator; Linux x86-64 only, because the triangulation differs by architecture with the same numpy and triangle versions): the shell-assembly
     generator, run with the options of PROVENANCE-generator.tsv on the regenerated specimen, must write the committed table. Skipped, with the reason printed, when the machine is not
     Linux x86-64, when `triangle` is missing, or when the generator on this checkout does not have the options (instrument, hole_seed, foot_tol_mm, foot_band_mm, crossing_snap_mm, rib_layers, bridge_dp_mm).

Exit status 1 when any check that ran failed. Skips are printed, never silent.
"""
import argparse
import hashlib
import os
import platform
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TABLE_MD5 = "3c12b931275bce326db4e0e58efae4bb"
NOTES = (("E", "E1"), ("A", "A1"), ("D", "D2"), ("G", "G2"))
SPRUCE_LB = "violin-spruce-norway-viala2020-shell-declared-lowband.fsmatpk"
MAPLE_LB = "violin-maple-sycamore-sonderegger2013-shell-declared-lowband.fsmatpk"
GEN_OPTIONS = ["instrument=bass", "hole_seed=median-run", "neck=ellipse", "fanfix=move", "rib_layers=6", "foot_band_mm=13", "foot_tol_mm=1.5", "bridge_dp_mm=0.3", "crossing_snap_mm=3.8"]
FAILED = []


def md5(path):
    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def report(tag, ok, text):
    print("%s %s: %s" % ("PASS" if ok else "FAIL", tag, text))
    if not ok:
        FAILED.append(tag)


def skip(tag, why):
    print("SKIP %s: %s" % (tag, why))


def md5_list(path):
    out = {}
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if line:
            h, name = line.split(None, 1)
            out[name.strip().lstrip("*")] = h
    return out


def run(cmd, **kw):
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, **kw)


def check_c(geo, cannone, tmp):
    spec = os.path.join(tmp, "specimen")
    r = run([sys.executable, "-I", os.path.join(HERE, "bass_specimen.py"), cannone, spec])
    if r.returncode != 0:
        report("C", False, "bass_specimen.py exited %d: %s" % (r.returncode, r.stderr.strip()[-300:]))
        return None
    want = md5_list(os.path.join(geo, "specimen-files.md5"))
    got = {n: md5(os.path.join(spec, n)) for n in sorted(os.listdir(spec))}
    bad = sorted(n for n in set(want) | set(got) if want.get(n) != got.get(n))
    report("C", not bad, "%d specimen files regenerated, %d differ%s" % (len(got), len(bad), (": " + ", ".join(bad)) if bad else ""))
    same = md5(os.path.join(geo, "PROVENANCE.tsv")) == got.get("PROVENANCE.tsv")
    report("C", same, "the committed PROVENANCE.tsv equals the regenerated one (md5 %s)" % got.get("PROVENANCE.tsv", "?")[:12])
    return spec


def read_table(path):
    kinds = {}
    cavity = None
    totals = None
    for line in open(path, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        kinds[f[0]] = kinds.get(f[0], 0) + 1
        if f[0] == "cavity":
            cavity = f[1:4]
        if f[0] == "meta" and f[1] == "totals":
            totals = f[2]
    return kinds, cavity, totals


def check_t(geo):
    table = os.path.join(geo, "bass-h19-s2d.tsv")
    report("T", md5(table) == TABLE_MD5, "table md5 %s" % md5(table)[:12])
    kinds, cavity, totals = read_table(table)
    n = (kinds.get("component"), kinds.get("node"), kinds.get("tri"), kinds.get("tie"))
    want = "%d components, %d nodes, %d triangles, %d ties" % n
    report("T", totals == want, "counts %s; the totals row says %r" % (want, totals))
    body = set()
    for s, note in NOTES:
        card = open(os.path.join(geo, "bass-shell-lb-%s.material-strings" % note), encoding="utf-8").read().split("\n")
        b = [l for l in card if l.startswith("body ")]
        c = [l for l in card if l.startswith("cavity ")]
        ok = len(b) == 1 and len(c) == 1 and b[0].split(" ")[1] == SPRUCE_LB and c[0].split(" ")[:5] == ["cavity", "physical-neck"] + cavity[:3]
        report("T", ok, "card %s: spruce variant on the body record, cavity row = the table's %s" % (note, " ".join(cavity)))
        if b:
            body.add(" ".join(b[0].split(" ")[2:]))
        opt = open(os.path.join(geo, "bass-h19-%s-lb.assembly" % s), encoding="utf-8").read().split("\n")
        packs = [l for l in opt if l.startswith("pack ")]
        bridges = [l for l in opt if l.startswith("bridge ")]
        report("T", packs == ["pack spruce @body", "pack maple " + MAPLE_LB] and len(bridges) == 1, "assembly %s: maple variant, one bridge point (%s)" % (s, bridges[0] if bridges else "none"))
    report("T", len(body) == 1, "the four cards share one body record: %s" % (sorted(body)[0] if body else "none"))


def check_e(geo, binary, tmp):
    table = os.path.join(geo, "bass-h19-s2d.tsv")
    want = md5_list(os.path.join(geo, "profiles.md5"))
    _, cavity, _ = read_table(table)
    for s, note in NOTES:
        out = os.path.join(tmp, "bass-h19-%s-lb.profile" % s)
        r = run([binary, "assemble-shell", table, os.path.join(geo, "bass-h19-%s-lb.assembly" % s), out])
        if r.returncode != 0:
            report("E", False, "assemble-shell %s exited %d: %s" % (s, r.returncode, r.stderr.strip()[-200:]))
            continue
        printed = dict(l.split(None, 1) for l in r.stdout.splitlines() if " " in l)
        card = open(os.path.join(geo, "bass-shell-lb-%s.material-strings" % note), encoding="utf-8").read().split("\n")
        body = [l for l in card if l.startswith("body ")][0].split(" ")
        same = body[2] == printed.get("body-length-m", "").strip() and body[3] == printed.get("body-width-m", "").strip()
        cav = [printed.get("cavity-volume-m3", "").strip(), printed.get("cavity-neck-radius-m", "").strip(), printed.get("cavity-neck-length-m", "").strip()] == cavity
        name = "bass-h19-%s-lb.profile" % s
        report("E", md5(out) == want[name], "profile %s md5 %s (expected %s)" % (s, md5(out)[:12], want[name][:12]))
        report("E", same and cav, "profile %s: printed body length, width and cavity row equal the card's and the table's" % s)


def check_b(geo, generator, spec, tmp):
    if platform.system() != "Linux" or platform.machine() not in ("x86_64", "AMD64"):
        return skip("B", "the table reproduces on Linux x86-64 only (this is %s %s)" % (platform.system(), platform.machine()))
    text = open(generator, encoding="utf-8").read()
    missing = [o.split("=")[0] for o in GEN_OPTIONS if o.split("=")[0] not in text]
    if missing:
        return skip("B", "the generator on this checkout lacks the options %s" % ", ".join(missing))
    out = os.path.join(tmp, "bass-h19-s2d.tsv")
    r = run([sys.executable, "-I", generator, spec, out, "19"] + GEN_OPTIONS)
    if r.returncode != 0:
        return report("B", False, "the generator exited %d: %s" % (r.returncode, r.stderr.strip()[-300:]))
    report("B", md5(out) == TABLE_MD5, "generated table md5 %s (committed %s)" % (md5(out)[:12], TABLE_MD5[:12]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=os.path.abspath(os.path.join(HERE, "..", "..", "..", "..")))
    ap.add_argument("--bin")
    ap.add_argument("--generator", action="store_true")
    a = ap.parse_args()
    geo = os.path.join(a.repo, "data", "instruments", "geometry", "double-bass-johannsson-shell")
    cannone = os.path.join(a.repo, "data", "instruments", "geometry", "violin-cannone-guarneri-1743-esrf2026")
    with tempfile.TemporaryDirectory() as tmp:
        spec = check_c(geo, cannone, tmp)
        check_t(geo)
        if a.bin:
            check_e(geo, a.bin, tmp)
        else:
            skip("E", "no --bin given (needs material_strings from the R5C pin or later)")
        if a.generator and spec:
            check_b(geo, os.path.join(HERE, "gen_violin_shell_assembly.py"), spec, tmp)
        else:
            skip("B", "no --generator given")
    print("FAILED: " + ", ".join(sorted(set(FAILED))) if FAILED else "all checks that ran passed")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
