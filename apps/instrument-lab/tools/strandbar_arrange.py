"""Strandbar groove (Todd Terje, 2013) for the instrument lab.

Piano rhythm and voicings: music_theory_data_tool's case study of the record (RH stabs on sixteenths
0 3 5 7 9 12 14, LH pickups on 2 and 11, kick on every beat, second chord on slot 7).
Tempo and pitch: the reference recording plays it 130.45 BPM, a semitone up (+1 here), and mostly
vamps on the first loop bar, so four vamp bars precede loop bars A-D.
Hi-hat accents: per-slot onset strength of the recording's 6-14 kHz band (60-120 s and 240-290 s).

usage: python3 tools/strandbar_arrange.py OUT_DIR
writes strandbar-arrangement.json (for tools/strandbar.py), strandbar-piano.csv (grand_piano
--performance), hits-kick.csv and hits-hihat-closed.csv (fs-vkplate perform) and duration.txt."""
import json, sys
from pathlib import Path
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
BPM, SHIFT, SR = 130.45, 1, 48000
PRE = 0.15  # seconds of silence before the downbeat of bar 1
NAMES = {"C": 0, "Db": 1, "D": 2, "Eb": 3, "E": 4, "F": 5, "Gb": 6, "G": 7, "Ab": 8, "A": 9, "Bb": 10, "B": 11}
def m(n):
    return 12 * (int(n[-1]) + 1) + NAMES[n[:-1]] + SHIFT
BARS = {
    "A": [(0, "Ab3 C4 Eb4 F4", "Db2"), (7, "G3 C4 Eb4 F4", "Ab2")],
    "B": [(0, "F3 Bb3 D4", "G2"), (7, "F3 A3 C4 F4", "D2")],
    "C": [(0, "Gb3 Bb3 F4 Db5", "Eb2"), (7, "F3 Ab3 C4 Eb4 F4", "F2")],
    "D": [(0, "Bb3 D4 F4", "G2"), (5, "G3 Bb3 Db4 Eb4", "Db2"), (7, "A3 D4 E4 F4 C5", "D2")],
}
FORM = ["A", "A", "A", "A", "A", "B", "C", "D"]
STABS, PICKUPS = [0, 3, 5, 7, 9, 12, 14], [2, 11]
HAT = [0.44, 0.28, 0.86, 0.36, 0.91, 0.58, 1.00, 0.63, 0.53, 0.65, 0.99, 0.38, 0.55, 0.81, 0.28, 0.68]
piano, drums = [], []
for b, name in enumerate(FORM):
    segs = BARS[name]
    seg_at = lambda slot: [s for s in segs if s[0] <= slot][-1]
    for slot in STABS:
        for n in seg_at(slot)[1].split():
            piano.append({"t": b * 4 + slot / 4, "dur": 0.3, "midi": m(n), "vel": 0.86, "part": "rh"})
    for i, (s0, _, bass) in enumerate(segs):
        s1 = segs[i + 1][0] if i + 1 < len(segs) else 16
        piano.append({"t": b * 4 + s0 / 4, "dur": (s1 - s0) / 4 - 0.12, "midi": m(bass), "vel": 0.8, "part": "bass"})
    for slot in PICKUPS:
        # the octave above the bass, unless the right hand strikes that key a sixteenth later or the
        # bass holds it (a real action cannot restrike a key that has not risen): then the fifth
        bass = m(seg_at(slot)[2])
        busy = {m(n) for n in seg_at(slot + 1)[1].split()} | {bass}
        pick = next(c for c in (bass + 12, bass + 7, bass + 19) if c not in busy)
        piano.append({"t": b * 4 + slot / 4, "dur": 0.2, "midi": pick, "vel": 0.7, "part": "lh"})
    for beat in range(4):
        drums.append({"t": b * 4 + beat, "piece": "kick", "vel": 0.9})
    for s in range(16):
        # a small, repeatable spread so sixteen identical strokes do not sound like a machine
        wob = ((b * 16 + s) * 2654435761 % 1000) / 1000 - 0.5
        drums.append({"t": b * 4 + s / 4, "piece": "hihat-closed", "vel": round(HAT[s] * (1 + 0.08 * wob), 3)})
piano.sort(key=lambda e: (e["t"], e["midi"]))
beat_s = 60 / BPM
bars = len(FORM)
loop_s = bars * 4 * beat_s
arr = {"title": "Strandbar", "after": "Todd Terje, Strandbar (2013)", "bpm": BPM, "beats_per_bar": 4, "bars": bars,
       "form": FORM, "transpose_semitones": SHIFT, "pre_s": PRE, "loop_s": loop_s, "piano": piano, "drums": drums}
(OUT / "strandbar-arrangement.json").write_text(json.dumps(arr, indent=1))
ev = []
for e in piano:
    t0 = PRE + e["t"] * beat_s
    ev.append((round(t0 * SR), 1, f"note_on,{e['midi']},{0.6 + 3.4 * e['vel']:.3f}"))
    ev.append((round((t0 + e["dur"] * beat_s) * SR), 0, f"note_off,{e['midi']},0"))
ev.sort()
(OUT / "strandbar-piano.csv").write_text("sample,event,key,value\n" + "".join(f"{s},{x}\n" for s, _, x in ev))
for piece, speed in (("kick", lambda v: 3.0 * v / 0.9), ("hihat-closed", lambda v: 1.2 + 1.6 * v)):
    rows = [f"{PRE + d['t'] * beat_s:.6f},{speed(d['vel']):.4f}" for d in drums if d["piece"] == piece]
    (OUT / f"hits-{piece}.csv").write_text("time_s,speed_m_s\n" + "\n".join(rows) + "\n")
total = PRE + loop_s + 2.5
(OUT / "duration.txt").write_text(f"{total:.3f}")
print(f"{len(piano)} piano notes, {len(drums)} drum hits, loop {loop_s:.3f} s, render {total:.2f} s")
