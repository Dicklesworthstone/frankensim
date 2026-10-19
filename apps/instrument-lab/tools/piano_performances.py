"""Write the grand_piano performance CSVs behind the page's piano clips.

usage: python tools/piano_performances.py OUT_DIR

Each file is `sample,event,key,value` on the 48 kHz output clock (a2-solo: 24 kHz). Render them with
the graded flags (README, "Grading against real recordings"):
  hero.csv            Satie, Gymnopedie No. 1 (1888, public domain): bars 1-2, then the melody of bars 5-9;
                      stereo (--microphone 0.2,0.8,1.0 --microphone-right 1.2,0.8,1.0), --duration from hero.dur
  a3-{pp,mf,ff}.csv   one A3 strike at 0.6 / 2.0 / 4.5 m/s; a0-mf.csv, c6-mf.csv at 2.0 m/s
  silent-held.csv     C4 E4 G4 held down silently, C3 struck staccato; silent-control.csv the same strike alone,
                      silent-pedal.csv with the sustain pedal down (--duration 4)
  octaves-high.csv    A4 + A5 struck and held, then A5 + A6 (--duration 5.4; again with --stretch-tuning 440)
  a2-solo.csv         one A2 strike for the soundboard mixer (--sample-rate 24000 --board-band-hz 1200
                      --modal-pressure-csv / --pressure-basis-json)
"""
import sys

SR = 48000


def write(path, events, sr=SR):
    rows = [(int(round(t * sr)), ev, key, val) for t, ev, key, val in events]
    order = {"sustain": 0, "sostenuto": 1, "una_corda": 2, "note_off": 3, "silent_key_down": 4, "note_on": 5}
    rows.sort(key=lambda r: (r[0], order.get(r[1], 9)))
    with open(path, "w") as f:
        f.write("sample,event,key,value\n")
        for s, ev, key, val in rows:
            f.write(f"{s},{ev},{key},{val:g}\n")


def note(ev, t, key, vel, dur):
    ev.append((t, "note_on", key, vel))
    ev.append((t + dur, "note_off", key, 0))


# Satie: a crotchet of 0.78 s, three to the bar; bass on beat 1, chord on beat 2, melody from bar 3.
Q = 0.78
BAR = 3 * Q
G2, D2 = 43, 38
CH1, CH2 = [59, 62, 66], [57, 61, 66]
MELODY = {2: [(1, 78, 1), (2, 81, 1)], 3: [(0, 79, 1), (1, 78, 1), (2, 73, 1)], 4: [(0, 71, 1), (1, 73, 1), (2, 74, 1)],
          5: [(0, 69, 3)], 6: [(0, 66, 3)]}


def hero(nbars):
    ev = []
    for b in range(nbars):
        t0 = 0.05 + b * BAR
        # Pedal change at each bar line: lift just before, catch just after the bass note.
        if b > 0:
            ev.append((t0 - 0.06, "sustain", 0, 0))
        ev.append((t0 + 0.09, "sustain", 0, 1))
        bass, chord = (G2, CH1) if b % 2 == 0 else (D2, CH2)
        note(ev, t0, bass, 1.25, BAR - 0.35)
        # A key the melody is holding cannot be struck again by the chord (bar 7: the melody's F#5 is down).
        held = {key for beat, key, dur in MELODY.get(b, []) if beat < 1 < beat + dur}
        for k in chord:
            if k not in held:
                note(ev, t0 + Q, k, 0.85 if k != 66 else 0.9, 2 * Q - 0.32)
        for beat, key, dur in MELODY.get(b, []):
            note(ev, t0 + beat * Q, key, 1.7, dur * Q - 0.12)
    end = 0.05 + nbars * BAR
    ev.append((end + 1.6, "sustain", 0, 0))
    return ev, end + 2.6


out = sys.argv[1]
ev, duration = hero(7)
write(f"{out}/hero.csv", ev)
with open(f"{out}/hero.dur", "w") as f:
    f.write(f"{duration:.2f}")
for name, key, vel in [("a3-pp", 57, 0.6), ("a3-mf", 57, 2.0), ("a3-ff", 57, 4.5), ("a0-mf", 21, 2.0), ("c6-mf", 84, 2.0)]:
    write(f"{out}/{name}.csv", [(0, "note_on", key, vel)])
write(f"{out}/a2-solo.csv", [(0, "note_on", 45, 2.5)], sr=24000)
strike = [(0.10, "note_on", 48, 4.0), (0.32, "note_off", 48, 0)]
write(f"{out}/silent-held.csv", [(0, "silent_key_down", k, 0) for k in (60, 64, 67)] + strike)
write(f"{out}/silent-control.csv", strike)
write(f"{out}/silent-pedal.csv", [(0, "sustain", 0, 1)] + strike)
write(f"{out}/octaves-high.csv", [(0.05, "note_on", 69, 2.2), (0.05, "note_on", 81, 2.2), (2.30, "note_off", 69, 0),
                                  (2.30, "note_off", 81, 0), (2.60, "note_on", 81, 2.2), (2.60, "note_on", 93, 2.2)])
print(f"wrote the page's piano performances to {out}; hero lasts {duration:.2f} s")
