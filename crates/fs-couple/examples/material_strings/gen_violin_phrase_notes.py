"""Write one material-strings note per violin phrase render on the il Cannone body (MagentaBadger 7268/7276).

usage: python3 -I gen_violin_phrase_notes.py REPO_ROOT BASE_NOTE OUT_DIR [--stop K C] [--note STRING SEMITONES SECONDS]

BASE_NOTE is the open-G3 note on the h6 body: the committed cannone-dominant-g3.material-strings with its body and
cavity records replaced by what `material_strings assemble-shell` prints for shell-assembly-h6.tsv. Every note keeps
BASE_NOTE's header unchanged (body, cavity, strings and their packs, gas, friction, prepared budget, audio clock, full
scale and callback) except the audio horizon, and gets its own bow schedule. Render one note with
  material_strings NOTE OUT.wav --body-profile cannone-h6.body

The phrases are those of the instrument-tests site (scale, Ode to Joy at 100 bpm, four double stops in G major). Each
note sounds for its duration and then rings for 0.45 s with the bow lifted. Strings 0-3 are G D A E. A note n semitones
above an open string is played on the string with the smallest n >= 0, and is named s{string}-o{n}-d{duration}.
Double stops are two bow tracks in one note, named chord{k}-d{duration}.

The bow, the same rule as the site's bowed cards (MagentaBadger 7276):
- set down at the speed 0 with the full force of data/instruments/descriptions/violin/arco-open-{string}.performance;
- the speed ramps to that card's speed in 80 ms at constant force;
- the force is released over the last 20 ms of the note;
- the station is 0.1 of the speaking length from the bridge. The bridge is each string row's own station, and a
  stopped string speaks over 2^(-n/12) of the open length.

The finger stop (MagentaBadger 7276 item 1), with --stop K C. Without it, stopped notes are listed, not written.
- Each stopped string gets two opposed stops along the bow tangent, at station bridge (1 - 2^(-n/12)) on a fixed
  receiver (the table has no fingerboard or neck; stated), with exponent 1. A stopped note NAME.material-strings
  renders with --stop-profile NAME.stops.
- The pair is a linear spring K [N/m] and dashpot C [N s/m] (the pad's measured shear damping). The stop's Hunt-Crossley
  loss chi K p dp/dt vanishes at zero penetration, so each side is preloaded into the string by p0 = 10 N / K with
  stiffness K / 2 and chi = C / (K p0). While both sides stay engaged (|x| < p0) their sum is exactly
  -K x - C dx/dt. The 5 N preloads cancel at rest.
- K is a STATED UNKNOWN: the ladder 1e4, 1e5, 1e6 N/m approaches the rigid limit (the string sticks on the ebony).
- The plane normal to the bow carries no stop: the bow drives only its own plane.
--note writes that one note (open or stopped) instead of the phrase set.
"""
import math
import os
import sys
from decimal import Decimal

REPO, BASE, OUT = sys.argv[1].rstrip("/"), sys.argv[2], sys.argv[3].rstrip("/")
options = sys.argv[4:]
STOP, ONLY = None, None
while options:
    if options[0] == "--stop" and len(options) >= 3:
        STOP, options = (float(options[1]), float(options[2])), options[3:]
    elif options[0] == "--note" and len(options) >= 4:
        ONLY, options = (int(options[1]), int(options[2]), Decimal(options[3])), options[4:]
    else:
        sys.exit(__doc__)
if STOP and not (STOP[0] > 0 and STOP[1] >= 0):
    sys.exit("--stop needs K > 0 and C >= 0")
PRELOAD_N = 5.0  # each side's preload force K/2 p0, so p0 = 10 N / K
CARDS = os.path.join(REPO, "data/instruments/descriptions/violin")
OPEN = ["G3", "D4", "A4", "E5"]
TAIL = Decimal("0.45")  # seconds rendered after the bow lifts
RAMP = Decimal("0.08")  # the speed reaches the card's value in 80 ms
RELEASE = Decimal("0.02")  # the force falls to zero over the note's last 20 ms
BETA = 0.1  # bow station, as a fraction of the speaking length from the bridge
SCALE = ["G3", "A3", "B3", "C4", "D4", "E4", "F#4", "G4"]
ODE = ["B4", "B4", "C5", "D5", "D5", "C5", "B4", "A4", "G4", "G4", "A4", "B4", "B4", "A4", "A4"]
ODE_BEATS = ["1"] * 12 + ["1.5", "0.5", "2"]
BEAT = Decimal("0.6")  # 100 bpm
STOPS = [[(0, 0), (1, 0)], [(1, 0), (2, 2)], [(1, 4), (2, 0)], [(1, 5), (2, 2)]]
STOP_DURATIONS = [Decimal("1.6")] * 3 + [Decimal("2.4")]
NOTE = {"C": 0, "C#": 1, "D": 2, "D#": 3, "E": 4, "F": 5, "F#": 6, "G": 7, "G#": 8, "A": 9, "A#": 10, "B": 11}


def midi(name):
    return 12 * (int(name[-1]) + 1) + NOTE[name[:-1]]


def pick_string(name):
    return min((midi(name) - midi(s), i) for i, s in enumerate(OPEN) if midi(name) >= midi(s))[::-1]


def rust_e(x):
    """Rust's `{:e}` for an f64: the shortest round-trip digits, which the canonical gesture bytes require."""
    x = float(x)
    if x == 0:
        return "0e0"
    sign, digits, exp = Decimal(repr(x)).as_tuple()
    digits = list(digits)
    while len(digits) > 1 and digits[-1] == 0:
        digits.pop()
        exp += 1
    m = str(digits[0]) + ("." + "".join(map(str, digits[1:])) if len(digits) > 1 else "")
    return ("-" if sign else "") + f"{m}e{exp + len(digits) - 1}"


def card_bow(string):
    """The site card's force (its initial row) and speed (its first event), as written."""
    rows = [line.split("\t") for line in open(os.path.join(CARDS, f"arco-open-{OPEN[string]}.performance"))
            if line.startswith(("initial\t", "event\t"))]
    initial = next(r for r in rows if r[0] == "initial")
    first = next(r for r in rows if r[0] == "event")
    if initial[1] != "bow" or first[3] != "bow":
        sys.exit(f"{OPEN[string]}: the card's bow rows are not bow values")
    return float(first[4]), float(initial[3])


header, schedule = open(BASE, encoding="utf-8").read().split("\nschedule\n")
header = header.split("\n")
audio = next(i for i, line in enumerate(header) if line.startswith("audio "))
_, rate, _, full_scale, block = header[audio].split(" ")
rate = int(rate)
rows = [line.split(" ") for line in header if line.startswith("string ")]
if [r[1].split("-")[3] for r in rows] != ["g", "d", "a", "e"]:
    sys.exit("BASE_NOTE must carry the four Dominant strings in the order G D A E")
bridge = [float(r[8]) for r in rows]
control = next(line for line in schedule.split("\n") if line.startswith("control_rate_hz\t"))
base_initial = next(line for line in schedule.split("\n") if line.startswith("initial\t")).split("\t")
base_track = next(line for line in schedule.split("\n") if line.startswith("track\t")).split("\t")


def station(string, semitones):
    return bridge[string] * (1.0 - BETA * 2.0 ** (-semitones / 12))


# The rule reproduces the base note's own open bow station bit for bit.
if rust_e(station(int(base_track[3]), 0)) != base_initial[4]:
    sys.exit(f"the 0.1-of-length rule gives {rust_e(station(int(base_track[3]), 0))}, the base note {base_initial[4]}")


def track(track_id, string, semitones, duration):
    speed, force = card_bow(string)
    at = rust_e(station(string, semitones))
    return [f"track\t{track_id}\tbow\t{string}",
            f"initial\tbow\t0e0\t{rust_e(force)}\t{at}",
            "events\t2",
            f"event\t0e0\t{rust_e(RAMP)}\tbow\t{rust_e(speed)}\t{rust_e(force)}\t{at}",
            f"event\t{rust_e(duration - RELEASE)}\t{rust_e(RELEASE)}\tbow\t{rust_e(speed)}\t0e0\t{at}"]


def stops(parts):
    """The finger-stop profile: an opposed preloaded pair on every stopped string."""
    k, c = STOP
    p0 = 2.0 * PRELOAD_N / k
    lines = ["frankensim-string-stops-v1"]
    for string, semitones in parts:
        if semitones == 0:
            continue
        at = bridge[string] * (1.0 - 2.0 ** (-semitones / 12))
        for side in ("positive", "negative"):
            # id string station side receiver clearance stiffness alpha loss max-reaction max-penetration
            lines.append(" ".join(["stop", f"finger{string}{side[0]}", str(string), repr(at), side, "fixed",
                                   repr(-p0), repr(0.5 * k), "1", repr(c / (k * p0)), repr(10.0 * PRELOAD_N),
                                   repr(3.0 * p0)]))
    return "\n".join(lines) + "\n"


written, pending = {}, set()


def note(name, parts, duration):
    """parts: (string, semitones) per bow track."""
    stopped = any(n > 0 for _, n in parts)
    if stopped and STOP is None:
        pending.add(name)
        return
    samples = (duration + TAIL) * rate
    if samples != samples.to_integral_value():
        sys.exit(f"{name}: {duration} s plus the tail is not a whole number of samples at {rate} Hz")
    lines = list(header)
    lines[audio] = f"audio {rate} {int(samples)} {full_scale} {block}"
    ids = ["bow"] if len(parts) == 1 else [f"bow{s}" for s, _ in parts]
    lines += ["schedule", "frankensim-gesture-schedule-v1", control, f"tracks\t{len(parts)}"]
    for track_id, (string, semitones) in zip(ids, parts):
        lines += track(track_id, string, semitones, duration)
    text = "\n".join(lines) + "\n"
    if written.setdefault(name, text) != text:
        sys.exit(f"{name}: two renders share one name")
    if stopped:
        written[name.removesuffix(".material-strings") + ".stops"] = stops(parts)


def single(name, duration):
    string, semitones = pick_string(name)
    note(f"violin/s{string}-o{semitones}-d{duration.normalize():f}.material-strings",
         [(string, semitones)], duration)


if ONLY:
    string, semitones, duration = ONLY
    note(f"violin/s{string}-o{semitones}-d{duration.normalize():f}.material-strings", [(string, semitones)], duration)
else:
    for name, duration in zip(SCALE + SCALE[-2::-1], [Decimal("0.5")] * 7 + [Decimal("0.8")]
                              + [Decimal("0.5")] * 6 + [Decimal("1.3")], strict=True):
        single(name, duration)
    for name, beats in zip(ODE, ODE_BEATS, strict=True):
        single(name, Decimal(beats) * BEAT)
    for k, (pair, duration) in enumerate(zip(STOPS, STOP_DURATIONS, strict=True)):
        note(f"violin/chord{k}-d{duration.normalize():f}.material-strings", pair, duration)

for name, text in sorted(written.items()):
    path = os.path.join(OUT, name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "x", encoding="utf-8") as f:
        f.write(text)
    print("wrote", name)
for name in sorted(pending):
    print("pending the stop law", name)
