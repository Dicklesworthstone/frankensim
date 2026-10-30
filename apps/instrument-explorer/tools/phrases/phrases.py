"""Write the phrase renders for the instrument-tests site and the plan that joins them.

usage: python3 phrases.py SRC_DESCRIPTIONS_DIR OUT_DIR

Every input is derived from a committed description at eaf6e6539; nothing is tuned to a recording.
  bowed   one card per note: the open-string card with the string shortened to 2^(-n/12) of its length
          (a rigid stop, no finger), mode count scaled so the top mode stays where the open string's was,
          its own duration and a bow lift at the end. Double stops are two cards in one `ensemble` run.
  wind    the fingering's own file with its duration and attack/release moved; nothing else changes.
  voice   the vowel/pitch yielding-wall file render.sh grades, with its duration and subglottal schedule moved.
  lip     ONE render per phrase: the head plus a gesture schedule that steps lip pre-stress and tongues
          the mouth pressure per note (the pre-stress and pressures of render.sh's held notes).
OUT_DIR/jobs.txt: `<mode> <out.wav> <args...>` per render; OUT_DIR/plan.json: how each phrase is joined.
"""
import json
import math
import os
import sys
from decimal import Decimal

SRC, OUT = sys.argv[1], sys.argv[2]
NOTE = {"C": 0, "C#": 1, "Db": 1, "D": 2, "D#": 3, "Eb": 3, "E": 4, "F": 5, "F#": 6, "Gb": 6, "G": 7,
        "G#": 8, "Ab": 8, "A": 9, "A#": 10, "Bb": 10, "B": 11}


def midi(name):
    return 12 * (int(name[-1]) + 1) + NOTE[name[:-1]]


def nm(m):
    return ["C", "C#", "D", "Eb", "E", "F", "F#", "G", "Ab", "A", "Bb", "B"][m % 12] + str(m // 12 - 1)


def e(x):
    """A gesture number in the canonical form the schedule parser demands (Rust's shortest `{:e}`)."""
    x = round(float(x), 6)
    if x == 0:
        return "0e0"
    sign, digits, exp = Decimal(repr(x)).normalize().as_tuple()
    m = str(digits[0]) + ("." + "".join(map(str, digits[1:])) if len(digits) > 1 else "")
    return ("-" if sign else "") + f"{m}e{exp + len(digits) - 1}"


def dec(x):
    """A plain decimal without trailing zeros, as the voice files write their subglottal schedule."""
    return repr(round(float(x), 6)).rstrip("0").rstrip(".")


def read(path):
    with open(path) as f:
        return f.read().splitlines()


def write(path, lines):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


jobs, plan = [], {}


def job(line):
    if line not in jobs:
        jobs.append(line)


# ---------------------------------------------------------------- bowed strings
BOWED = {
    "violin": dict(strings=["G3", "D4", "A4", "E5"], fs="0.01", key="G major",
                   scale=["G3", "A3", "B3", "C4", "D4", "E4", "F#4", "G4"],
                   ode=["B4", "B4", "C5", "D5", "D5", "C5", "B4", "A4", "G4", "G4", "A4", "B4", "B4", "A4", "A4"],
                   stops=[[(0, 0), (1, 0)], [(1, 0), (2, 2)], [(1, 4), (2, 0)], [(1, 5), (2, 2)]],
                   stop_names="G major: I (G-D), I6 (D-B), V (F#-A), I (G-B)"),
    "viola": dict(strings=["C3", "G3", "D4", "A4"], fs="0.02", key="C major",
                  scale=["C3", "D3", "E3", "F3", "G3", "A3", "B3", "C4"],
                  ode=["E4", "E4", "F4", "G4", "G4", "F4", "E4", "D4", "C4", "C4", "D4", "E4", "E4", "D4", "D4"],
                  stops=[[(0, 0), (1, 0)], [(1, 0), (2, 2)], [(1, 4), (2, 0)], [(1, 5), (2, 2)]],
                  stop_names="C major: I (C-G), I6 (G-E), V (B-D), I (C-E)"),
    "cello": dict(strings=["C2", "G2", "D3", "A3"], fs="0.04", key="C major",
                  scale=["C2", "D2", "E2", "F2", "G2", "A2", "B2", "C3"],
                  ode=["E3", "E3", "F3", "G3", "G3", "F3", "E3", "D3", "C3", "C3", "D3", "E3", "E3", "D3", "D3"],
                  stops=[[(0, 0), (1, 0)], [(1, 0), (2, 2)], [(1, 4), (2, 0)], [(1, 5), (2, 2)]],
                  stop_names="C major: I (C-G), I6 (G-E), V (B-D), I (C-E)"),
    "double-bass": dict(strings=["E1", "A1", "D2", "G2"], fs="0.1", key="G major",
                        scale=["G1", "A1", "B1", "C2", "D2", "E2", "F#2", "G2"],
                        ode=["B2", "B2", "C3", "D3", "D3", "C3", "B2", "A2", "G2", "G2", "A2", "B2", "B2", "A2", "A2"],
                        stops=[[(0, 3), (2, 0)], [(2, 0), (3, 4)], [(2, 4), (3, 2)], [(2, 5), (3, 4)]],
                        stop_names="G major: I (G-D), I6 (D-B), V (F#-A), I (G-B)"),
}
ODE_BEATS = [1] * 12 + [1.5, 0.5, 2]  # Beethoven, Symphony No. 9 (1824), public domain: first phrase
TAIL = 0.45  # seconds rendered after the bow lifts, so each note can ring into the next
RAMP = 0.08  # bow speed reaches its value in 80 ms (the committed cards use 150 ms for 3 s notes)
# The bow is set on the string with this share of the note's force, and the force grows with the speed over the ramp.
# The top of Schelleng's playable window scales with bow speed, so the full force at a near-zero speed sits above it:
# the open cello C2 then scrapes for 680 ms and settles twice as bright (centroid 1122 Hz, against 333 Hz once it
# speaks). Too little force for the acceleration falls under Guettler's lower limit instead (multiple slips). Swept on
# the renders' own lock time and steady regime, 0.02 to 1 of the force with the 80 ms ramp, one value per string
# (strings in BOWED order). The string has more than one steady regime (Helmholtz and brighter multiple-slip states)
# and the set-down picks one, not monotonically: each string takes the value that reaches its usual regime soonest.
# Cello C 0.1 (C2 in 20 ms, against 680 ms at full force), G 0.7 (G2 120, B2 20, C3 0 ms), D 0.35 (D3 0, E3 20 ms);
# bass E 0.2 (G1 480 ms, all slow), A 0.1 (A1 180 ms), D 0.5 (D2 200 ms), G 0.2 (A2 120, B2 160 ms). Real players:
# cello C2 120-280 ms, bass A1 40-120 ms. The violin and viola speak at once from full force, while a low set-down
# leaves the violin A string quiet and bright, sometimes on its second mode. 1 where no phrase bows the string.
SET_DOWN = {"violin": (1, 1, 1, 1), "viola": (1, 1, 1, 1), "cello": (0.1, 0.7, 0.35, 1), "double-bass": (0.2, 0.1, 0.5, 0.2)}
# Those lock times came from the audio, which cannot tell Helmholtz motion from a steady multiple slip. The solver's
# own slip traces (one slip per period, slip time beta of the period) show the cello D string at 0.35 in multiple
# slip for 0.9 s (D3 locks at 1.0 s) and the bass D and G strings out of Helmholtz motion for 0.3-0.5 s (bass F#2,
# G2, A2, B2). On these strings the card's own onset leaves less time out of it in sum over their double-stop notes
# (D3 0.25 s, F#2, G2, A2 0.15 s, B2 0.2 s; E3 0.5 s against 0.4 s; bass D2 0.3 s either way): the force constant from rest and the speed
# ramp of the card (1 m/s^2, its sources.tsv bow_acceleration from Guettler 2002), the gesture its force was derived
# for. Strings listed here keep the card's initial and first event as written; the others keep SET_DOWN and RAMP.
CARD_ONSET = {"cello": (2,), "double-bass": (2, 3)}


def bow_card(fam, si, off, dur):
    f = BOWED[fam]
    base = read(f"{SRC}/{fam}/arco-open-{f['strings'][si]}.performance")
    out, ratio = [], 2 ** (-off / 12)
    n_audio = round((dur + TAIL) * 48000)
    init = next(l for l in base if l.startswith("initial"))
    _, _, _, force, beta = init.split("\t")
    card_onset = si in CARD_ONSET.get(fam, ())
    ramp = next(l for l in base if l.startswith("event\t0e0\t")).split("\t")[2] if card_onset else e(RAMP)
    for l in base:
        if l.startswith("audio "):
            a = l.split()
            out.append(f"audio {a[1]} {n_audio} {a[3]}")
        elif l.startswith("string "):
            s = l.split()
            modes = max(8, math.floor(int(s[6]) * ratio))
            out.append(" ".join(["string", f"{float(s[1]) * ratio:.6f}"] + s[2:6] + [str(modes)]))
        elif l.startswith("damping "):
            d = l.split()
            # The committed cards damp every mode alike; a per-mode list would have to be resampled, not copied.
            assert len(set(d[1:])) == 1, f"{fam} {f['strings'][si]}: per-mode damping"
            out.append(" ".join(["damping"] + [d[1]] * modes))
        elif l.startswith("initial"):
            out.append(l if card_onset else f"initial\tbow\t0e0\t{e(SET_DOWN[fam][si] * float(force))}\t{beta}")
        elif l.startswith("events"):
            out.append("events\t2")
        elif l.startswith("event\t"):
            continue
        else:
            out.append(l)
    i = out.index(next(l for l in out if l.startswith("events")))
    out[i + 1:i + 1] = [f"event\t0e0\t{ramp}\tbow\t1.5e-1\t{force}\t{beta}",
                        f"event\t{e(dur - 0.02)}\t2e-2\tbow\t1.5e-1\t0e0\t{beta}"]
    name = f"{fam}/s{si}-o{off}-d{dur:g}.performance"
    write(f"{OUT}/{name}", out)
    return name


def pick_string(fam, m):
    opens = [midi(s) for s in BOWED[fam]["strings"]]
    best = min((m - o, i) for i, o in enumerate(opens) if m - o >= 0)
    return best[1], best[0]


for fam, f in BOWED.items():
    fs = f["fs"]
    # scale up and down, 0.5 s a note, the top and the last held
    seq = f["scale"] + f["scale"][-2::-1]
    durs = [0.5] * 7 + [0.8] + [0.5] * 6 + [1.3]
    notes, t = [], 0.0
    for name, d in zip(seq, durs, strict=True):
        si, off = pick_string(fam, midi(name))
        card = bow_card(fam, si, off, d)
        wav = card.replace(".performance", ".wav")
        job(f"ensemble {wav} --bow {card} --full-scale-pa {fs} --block 4096")
        notes.append(dict(wav=wav, at=round(t, 4), dur=d, note=name, string=f["strings"][si], semitones=off))
        t += d
    plan[f"{fam}/scale"] = dict(family=fam, kind="scale", title=f"{f['key']} scale, up and down",
                                notes=notes, join="detache")
    # Ode to Joy at 100 bpm
    notes, t = [], 0.0
    for name, b in zip(f["ode"], ODE_BEATS, strict=True):
        d = round(b * 0.6, 3)
        si, off = pick_string(fam, midi(name))
        card = bow_card(fam, si, off, d)
        wav = card.replace(".performance", ".wav")
        job(f"ensemble {wav} --bow {card} --full-scale-pa {fs} --block 4096")
        notes.append(dict(wav=wav, at=round(t, 4), dur=d, note=name, string=f["strings"][si], semitones=off))
        t += d
    plan[f"{fam}/melody"] = dict(family=fam, kind="melody", title="Ode to Joy (Beethoven), first phrase",
                                 notes=notes, join="detache")
    # double stops: two cards in one ensemble run each. The two parts' pressures add, so the PCM scale is twice the
    # one-string scale (at 0.04 Pa the cello C2+G2 render clipped 23 samples, peak 0.0443 Pa).
    notes, t = [], 0.0
    for k, pair in enumerate(f["stops"]):
        d = 1.6 if k < 3 else 2.4
        cards = [bow_card(fam, si, off, d) for si, off in pair]
        wav = f"{fam}/chord{k}-d{d:g}.wav"
        job(f"ensemble {wav} " + " ".join(f"--bow {c}" for c in cards) + f" --full-scale-pa {2 * float(fs):g} --block 4096")
        names = [nm(midi(f["strings"][si]) + off) for si, off in pair]
        notes.append(dict(wav=wav, at=round(t, 4), dur=d, note="+".join(names),
                          string="+".join(f["strings"][si] for si, _ in pair)))
        t += d
    plan[f"{fam}/chords"] = dict(family=fam, kind="chords", title=f"Double stops, {f['stop_names']}",
                                 notes=notes, join="detache")

# ---------------------------------------------------------------- flues (tin whistle, pan pipes)


def wind_note(fam, stem, dur, attack=0.03):
    base = read(f"{SRC}/{fam}/{stem}.performance")
    out = []
    n_audio = round((dur + 0.3) * 48000)
    for l in base:
        if l.startswith("output "):
            o = l.split()
            out.append(" ".join(["output", str(n_audio)] + o[2:]))
        elif l.startswith("event\t0e0\t") and "pressure-pa" in l:
            p = l.split("\t")
            out.append("\t".join([p[0], p[1], e(attack)] + p[3:]))
        elif l.startswith("event\t") and l.rstrip().endswith("0e0") and "pressure-pa" in l:
            p = l.split("\t")
            out.append("\t".join([p[0], e(dur)] + p[2:]))
        else:
            out.append(l)
    name = f"{fam}/{stem}-d{dur:g}.performance"
    write(f"{OUT}/{name}", out)
    return name


WHISTLE = {"D5": "d5", "E5": "e5", "F#5": "fs5", "G5": "g5", "A5": "a5", "B5": "b5", "C#6": "cs6", "D6": "d6"}
for kind, title, seq, durs in [
    ("scale", "D major scale, up and down", list(WHISTLE) + list(WHISTLE)[-2::-1], [0.4] * 7 + [0.7] + [0.4] * 6 + [1.1]),
    ("melody", "Ode to Joy (Beethoven), first phrase",
     ["F#5", "F#5", "G5", "A5", "A5", "G5", "F#5", "E5", "D5", "D5", "E5", "F#5", "F#5", "E5", "E5"],
     [round(b * 0.5, 3) for b in ODE_BEATS])]:
    notes, t = [], 0.0
    for name, d in zip(seq, durs, strict=True):
        card = wind_note("tin-whistle", f"tin-whistle-{WHISTLE[name]}", round(d - 0.04, 3))
        wav = card.replace(".performance", ".wav")
        job(f"wind {wav} {card}")
        notes.append(dict(wav=wav, at=round(t, 4), dur=d, note=name))
        t += d
    plan[f"tin-whistle/{kind}"] = dict(family="tin-whistle", kind=kind, title=title, notes=notes, join="tongued")

PIPES = ["D3", "G3", "D4", "G4", "B4"]
notes, t = [], 0.0
for name, d in zip(PIPES + PIPES[-2::-1], [0.6] * 4 + [0.9] + [0.6] * 3 + [1.3], strict=True):
    card = wind_note("pan-pipes", f"pan-pipes-{name.lower()}", round(d - 0.05, 3))
    wav = card.replace(".performance", ".wav")
    job(f"wind {wav} {card}")
    notes.append(dict(wav=wav, at=round(t, 4), dur=d, note=name))
    t += d
plan["pan-pipes/arpeggio"] = dict(family="pan-pipes", kind="arpeggio", title="G major arpeggio across five tubes, up and down",
                                  notes=notes, join="tongued")

# ---------------------------------------------------------------- voice


def voice_note(stem, dur):
    base = read(f"{SRC}/voice/{stem}.performance")
    out, done = [], False
    pressure = max((l.split()[2] for l in base if l.startswith("subglottal ")), key=float)
    for l in base:
        if l.startswith("audio "):
            a = l.split()
            out.append(f"audio {a[1]} {round((dur + 0.25) * 48000)} {a[3]}")
        elif l.startswith("subglottal "):
            if not done:
                out += ["subglottal 0 0", f"subglottal 0.06 {pressure}",
                        f"subglottal {dec(dur - 0.1)} {pressure}", f"subglottal {dec(dur)} 0"]
                done = True
        else:
            out.append(l)
    name = f"voice/{stem}-d{dur:g}.performance"
    write(f"{OUT}/{name}", out)
    return name


for pitch in ["C3", "C4"]:
    notes, t = [], 0.0
    for v in "aeiou":
        d = 0.9 if v != "u" else 1.4
        # The cards voice/render.sh grades: the cited yielding wall, at the 192 kHz clock except /a/ (render.sh header).
        card = voice_note(f"story2008-{v}-{pitch}-flanagan" + ("" if v == "a" else "-r4"), d + 0.15)
        wav = card.replace(".performance", ".wav")
        job(f"voice {wav} {card}")
        notes.append(dict(wav=wav, at=round(t, 4), dur=d, note=f"/{v}/ {pitch}"))
        t += d
    plan[f"voice/vowels-{pitch}"] = dict(family="voice", kind="vowels", title=f"Five vowels on {pitch}",
                                         notes=notes, join="legato")

# ---------------------------------------------------------------- lip valves: one render per phrase
LIP = {
    "trumpet": dict(head="trumpet-open.head", rate=32000, attack=1.4e4,
                    t={"Bb3": (0.93, 5.71e3), "F4": (1.42, 8.56e3), "Bb4": (1.9, 1.142e4), "D5": (2.42, 1.439e4)}),
}


def lip_phrase(fam, name, seq, tongued):
    L = LIP[fam]
    head = read(f"{SRC}/{fam}/{L['head']}")
    total = sum(d for _, d in seq) + 0.8
    pres, lip, t = [], [], 0.25
    first = L["t"][seq[0][0]][0]
    for k, (note, d) in enumerate(seq):
        tension, sustain = L["t"][note]
        if k:
            lip.append(f"event\t{e(t - 0.03)}\t3e-2\tmultiplier\t{e(tension)}")
        if tongued or k == 0:
            if k:
                pres.append(f"event\t{e(t - 0.05)}\t3e-2\tpressure-pa\t0e0")
            pres.append(f"event\t{e(t)}\t2e-2\tpressure-pa\t{e(L['attack'])}")
            if tongued and d > 0.4:
                pres.append(f"event\t{e(t + min(0.25, d / 3))}\t1.5e-1\tpressure-pa\t{e(sustain)}")
        t += d
    pres.append(f"event\t{e(t)}\t6e-2\tpressure-pa\t0e0")
    lines = ["frankensim-lip-performance-v1", f"audio {L['rate']} {round(L['rate'] * (total + 0.25))}"] + head + [
        "schedule", "frankensim-gesture-schedule-v1", "control_rate_hz\t700", "tracks\t2",
        "track\tblow\tblowing-pressure", "initial\tpressure-pa\t0e0", f"events\t{len(pres)}"] + pres + [
        "track\tlip\tpre-stress", f"initial\tmultiplier\t{e(first)}", f"events\t{len(lip)}"] + lip
    card = f"{fam}/{name}.performance"
    write(f"{OUT}/{card}", lines)
    wav = card.replace(".performance", ".wav")
    job(f"lip {wav} {card}")
    return wav


fan = [("F4", .45), ("Bb4", .45), ("D5", .45), ("Bb4", .45), ("F4", .45), ("Bb4", .45), ("D5", 1.2),
       ("D5", .45), ("Bb4", .45), ("F4", .45), ("Bb3", 1.6)]
plan["trumpet/fanfare"] = dict(family="trumpet", kind="melody", title="Bugle fanfare on the open notes, tongued",
                               native=lip_phrase("trumpet", "fanfare", fan, True), seq=fan)
slur = [("Bb3", .8), ("F4", .8), ("Bb4", .8), ("D5", 1.2), ("Bb4", .8), ("F4", .8), ("Bb3", 1.4)]
plan["trumpet/slur"] = dict(family="trumpet", kind="slur", title="Lip slur through the open notes, one breath",
                            native=lip_phrase("trumpet", "slur", slur, False), seq=slur)

# ---------------------------------------------------------------- rows cut from the playable-today holds
# (each family's committed descriptions rendered at one origin snapshot, copied into OUT/playable/<family>/: first
# ts1 mb-scales at eaf6e6539, now trj /data/tmp/mb-phrases/playable6 at d49a124c2, the guitar by its own render.sh)
CLAR = ["D3", "Eb3", "F3", "G3", "Ab3", "Bb3", "C4", "D4", "Eb4", "F4", "Gb4", "G4", "A4", "Bb4", "C5", "D5", "Eb5",
        "F5", "G5", "A5", "Bb5"]
ROWS = {
    "clarinet": ("Every fingering it has, from D3 up to Bb5", [(f"{i:02d}-{n}.wav", n) for i, n in enumerate(CLAR, 1)], 0.7),
    "flute": ("The three fingerings it has: D4, G4, A5", [("d4.wav", "D4"), ("g4.wav", "G4"), ("a5.wav", "A5")], 1.2),
    "recorder": ("The four fingerings it has: F4, A4, Bb4, B4", [("f1.wav", "F4"), ("a1.wav", "A4"), ("bb1.wav", "Bb4"),
                                                             ("b1.wav", "B4")], 1.0),
    "voice": ("/a/ at the six pitches it can sing", [(f"a-{n}.wav", n.replace("Ds", "D#")) for n in
                                                    ["C3", "Ds3", "A3", "C4", "D4", "F4"]], 1.0),
    # data/instruments/descriptions/guitar/render.sh: the plain nylon trebles only (the wound basses are blocked)
    "guitar": ("The five notes it has: three open strings (G3, B3, E4), then F#4 and F5 stopped on the B string",
               [("gt-s1-solo-g3.wav", "G3"), ("gt-s1-solo-b3.wav", "B3"), ("gt-s1-solo-e4.wav", "E4"),
                ("gt-s1-fret-gb4.wav", "F#4"), ("gt-s1-fret-f5.wav", "F5")], 1.6),
    "trombone": ("Six lip tensions in first position (the note is whatever the lip finds)",
                 [(f"hold-t{t}.wav", f"lip {v}") for t, v in [("5e-1", "0.5"), ("6e-1", "0.6"), ("7e-1", "0.7"),
                                                              ("8.5e-1", "0.85"), ("1e0", "1.0"), ("1.15e0", "1.15")]], 1.2),
}
for fam, (title, row, d) in ROWS.items():
    notes = [dict(wav=f"playable/{fam}/{w}", at=round(k * d, 3), dur=d, note=n) for k, (w, n) in enumerate(row)]
    plan[f"{fam}/row"] = dict(family=fam, kind="row", title=title, notes=notes, join="cut")

write(f"{OUT}/jobs.txt", jobs)
with open(f"{OUT}/plan.json", "w") as f:
    json.dump(plan, f, indent=1)
print(len(jobs), "renders,", len(plan), "phrases")
