"""The trumpet phrases play the graded holds' embouchure: each note's lip pre-stress and sustain pressure are the HOLDS
default of trumpet/render.sh. Run with
    python3 -I apps/instrument-explorer/tools/phrases/test_lip_holds.py [DESCRIPTIONS_DIR]
DESCRIPTIONS_DIR defaults to this checkout's data/instruments/descriptions."""
import json
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "../../../../data/instruments/descriptions")
NOTES = ["Bb3", "F4", "Bb4", "D5"]
ATTACK = 1.4e4


def holds(src):
    with open(os.path.join(src, "trumpet", "render.sh")) as f:
        (line,) = re.findall(r'^HOLDS=\$\{HOLDS-"([^"]*)"\}$', f.read(), re.M)
    return {n: (float(h.split(":")[0]), float(h.split(":")[2])) for n, h in zip(NOTES, line.split())}


def phrases(src, out):
    return subprocess.run([sys.executable, "-I", os.path.join(HERE, "phrases.py"), src, out],
                          capture_output=True, text=True, timeout=120)


def track(card, header):
    """A gesture track's values in time order: its initial value, then each event's."""
    with open(card) as f:
        lines = f.read().splitlines()
    i = lines.index(header)
    n = int(lines[i + 2].split("\t")[1])
    return [float(lines[i + 1].split("\t")[2])] + [float(x.split("\t")[4]) for x in lines[i + 3:i + 3 + n]]


def src_with_trumpet_render(text):
    """A descriptions dir whose families are links to SRC's, except trumpet/render.sh, which holds `text`."""
    tmp = tempfile.mkdtemp()
    for name in os.listdir(SRC):
        if name != "trumpet":
            os.symlink(os.path.abspath(os.path.join(SRC, name)), os.path.join(tmp, name))
    os.mkdir(os.path.join(tmp, "trumpet"))
    for name in os.listdir(os.path.join(SRC, "trumpet")):
        if name != "render.sh":
            os.symlink(os.path.abspath(os.path.join(SRC, "trumpet", name)), os.path.join(tmp, "trumpet", name))
    with open(os.path.join(tmp, "trumpet", "render.sh"), "w") as f:
        f.write(text)
    return tmp


def test_each_note_steps_to_its_holds_pre_stress_and_sustain():
    h = holds(SRC)
    out = tempfile.mkdtemp()
    run = phrases(SRC, out)
    assert run.returncode == 0, run.stderr
    with open(os.path.join(out, "plan.json")) as f:
        plan = json.load(f)
    for name in ("fanfare", "slur"):
        seq = [note for note, _ in plan[f"trumpet/{name}"]["seq"]]
        card = os.path.join(out, "trumpet", f"{name}.performance")
        got = track(card, "track\tlip\tpre-stress")
        assert got == [h[n][0] for n in seq], (name, got, h)
    # the tongued fanfare settles every note (all longer than 0.4 s) on its hold's sustain
    blow = track(os.path.join(out, "trumpet", "fanfare.performance"), "track\tblow\tblowing-pressure")
    assert {p for p in blow if p not in (0.0, ATTACK)} == {h[n][1] for n in NOTES}, (blow, h)


def test_a_render_sh_without_one_hold_per_note_is_refused():
    with open(os.path.join(SRC, "trumpet", "render.sh")) as f:
        text = f.read()
    line = re.search(r'^HOLDS=\$\{HOLDS-"([^"]*)"\}$', text, re.M)
    three = text.replace(line.group(1), " ".join(line.group(1).split()[:3]))
    for bad, why in ((text.replace(line.group(0), ""), "found 0"), (three, "is not one pre-stress:attack:sustain")):
        run = phrases(src_with_trumpet_render(bad), tempfile.mkdtemp())
        assert run.returncode != 0 and why in run.stderr, (why, run.returncode, run.stderr[-300:])


for name, fn in sorted(globals().copy().items()):
    if name.startswith("test_") and callable(fn):
        fn()
        print("ok", name)
print("all passed")
