"""Reproduce the cello shell data of data/instruments/geometry/cello-castelbarco-1697-shell byte for byte from origin's tables:
  C  cello_specimen.py on the il Cannone tables, then cello_bridge_19g.py, give every specimen table of the directory;
  B  gen_violin_shell_assembly.py instrument=cello on that specimen gives shell-assembly-h12.5.tsv;
  E  cello_shell_card.py with the directory's assemble-shell output leaves cello-shell-D3.material-strings unchanged
     (its body and cavity records are the table's);
  A  gen_violin_shell_assembly.py with no cello option still gives il Cannone's shell-assembly-h6.tsv.
The two tables reproduce only on the architecture that made them: with the same numpy 2.5.4 and triangle 20250106 the
triangulation differs between x86-64 Linux and arm64 macOS. The cello table was made on x86-64 Linux, so B runs there;
il Cannone's h6 table was made on arm64 macOS, so A runs there; elsewhere each is reported as skipped. C and E run anywhere.
assemble-shell needs the material_strings binary; its output for the D assembly is kept as assemble-shell-D.txt and is
used here only as E's input.
usage: uv run --python 3.12 --with triangle --with numpy python3 -I check_cello_shell.py WORK_DIR   (WORK_DIR must not exist)"""
import platform
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CANNONE = ROOT / "data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026"
DATA = ROOT / "data/instruments/geometry/cello-castelbarco-1697-shell"
CELLO_FLAGS = ["12.5", "neck=ellipse", "foot_band_mm=8", "fanfix=move", "crossing_snap_mm=2.5", "bridge_dp_mm=0.5", "rib_layers=6", "instrument=cello"]
WORK = Path(sys.argv[1])
if WORK.exists():
    raise SystemExit("refusing: %s exists" % WORK)
WORK.mkdir(parents=True)
failed = []
MACHINE = (sys.platform, platform.machine())
LINUX_X86 = MACHINE[0].startswith("linux") and MACHINE[1] == "x86_64"
MAC_ARM = MACHINE[0] == "darwin" and MACHINE[1] == "arm64"


def run(*args):
    r = subprocess.run([sys.executable, "-I", *map(str, args)], capture_output=True, text=True)
    if r.returncode:
        raise SystemExit("%s failed (rc %d):\n%s" % (Path(str(args[0])).name, r.returncode, r.stderr[-2000:]))


def same(label, made, kept):
    ok = made.read_bytes() == kept.read_bytes()
    print("%s %s: %s" % ("PASS" if ok else "FAIL", label, kept.relative_to(ROOT)))
    if not ok:
        failed.append(label)


# C: the specimen
run(HERE / "cello_specimen.py", CANNONE, WORK / "specimen-raw")
run(HERE / "cello_bridge_19g.py", WORK / "specimen-raw", WORK / "specimen")
for made in sorted((WORK / "specimen").iterdir()):
    same("C " + made.name, made, DATA / made.name)
# B: the cello table (made on x86-64 Linux)
if LINUX_X86:
    run(HERE / "gen_violin_shell_assembly.py", WORK / "specimen", WORK / "shell-assembly-h12.5.tsv", *CELLO_FLAGS)
    same("B table", WORK / "shell-assembly-h12.5.tsv", DATA / "shell-assembly-h12.5.tsv")
else:
    print("SKIP B table: made on x86-64 Linux, this is %s %s" % MACHINE)
# E: the D3 card's table-derived records
run(HERE / "cello_shell_card.py", DATA / "cello-shell-D3.material-strings", DATA / "assemble-shell-D.txt", WORK / "cello-shell-D3.material-strings")
same("E card", WORK / "cello-shell-D3.material-strings", DATA / "cello-shell-D3.material-strings")
# A: the violin table is unchanged by the cello options (made on arm64 macOS)
if MAC_ARM:
    run(HERE / "gen_violin_shell_assembly.py", CANNONE, WORK / "shell-assembly-h6.tsv", "6", "fanfix=move")
    same("A violin table", WORK / "shell-assembly-h6.tsv", CANNONE / "shell-assembly-h6.tsv")
else:
    print("SKIP A violin table: made on arm64 macOS, this is %s %s" % MACHINE)
print("FAILED: " + ", ".join(failed) if failed else "all checks pass")
sys.exit(1 if failed else 0)
