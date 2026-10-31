"""phrasekit.whole_cycles on planted tones whose centre is known. Run with
uv run --no-project --with numpy --with pillow python -I test_whole_cycles.py <this directory>"""
import sys

import numpy as np

sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else ".")
from phrasekit import SR, hz, pitch, whole_cycles  # noqa: E402

RATE = 5.9
PLAIN = [1.0, 0.5, 0.3, 0.2, 0.15, 0.1, 0.08, 0.05]
SEVENTH = [10 ** (-11 / 20), 0.3, 0.3, 0.4, 0.5, 0.7, 1.0, 0.5, 0.2]  # the 7th harmonic 11 dB over the 1st, as sung /a/ C3 has


def tone(note, centre, swing, amps, low_louder_db=0.0, secs=1.2):
    t = np.arange(int(secs * SR)) / SR
    s = np.sin(2 * np.pi * RATE * t + 0.7)
    phase = 2 * np.pi * np.cumsum(hz(note) * 2 ** ((centre + swing * s) / 1200)) / SR
    return sum(a * np.sin((k + 1) * phase) for k, a in enumerate(amps)) * 10 ** (-low_louder_db * s / 20)


def read(x, note):
    m = pitch(x, 0.16, 0.86, note)
    return m, whole_cycles(x, 0.26, 0.86, RATE, m, note)


def test_a_swing_louder_low_reads_its_centre_and_a_single_reading_does_not():
    m, w = read(tone("C4", -20, 83, PLAIN, 6.0), "C4")
    assert abs(w["cents"] + 20) < 1.3, w  # the mean of frequency sits about 1 cent above the mean of cents here
    assert abs(m["cents"] + 20) > 30, m  # the planted failure: pitch() is what this replaces on these notes
    assert w["whole_cycles"] == 3


def test_the_reading_does_not_depend_on_the_level():
    flat = read(tone("C4", -20, 83, PLAIN), "C4")[1]["cents"]
    swung = read(tone("C4", -20, 83, PLAIN, 6.0), "C4")[1]["cents"]
    assert abs(flat - swung) <= 0.1, (flat, swung)


def test_a_weak_fundamental_under_a_strong_seventh_still_reads():
    w = read(tone("C3", 12, 75, SEVENTH, 4.0), "C3")[1]
    assert abs(w["cents"] - 12) < 1.3, w


def test_a_straight_tone_reads_as_pitch_reads_it():
    m, w = read(tone("C4", -20, 0, PLAIN), "C4")
    assert abs(w["cents"] - m["cents"]) <= 0.1, (m, w)


def test_a_swing_has_a_steady_pitch():
    m, w = read(tone("C4", -20, 83, SEVENTH, 6.0), "C4")
    assert w["periodicity"] >= 0.5, w  # under 0.5 the phrase page says the note has no steady pitch
    assert m["periodicity"] < w["periodicity"] - 0.05, (m, w)  # the planted failure: one span across the swing reads lower


def test_a_span_under_one_cycle_keeps_the_first_reading():
    m = pitch(tone("C4", -20, 83, PLAIN), 0.16, 0.86, "C4")
    assert whole_cycles(tone("C4", -20, 83, PLAIN), 0.26, 0.40, RATE, m, "C4") == m


def test_a_vibrato_note_keeps_its_reading_through_the_analyzer():
    argv, sys.argv = sys.argv, ["analyze.py", ".", ".", ".", "."]  # analyze.py reads its directories on import
    import analyze
    sys.argv = argv
    x = tone("C4", -20, 83, PLAIN, 6.0)
    frames = len(x) // analyze.HOP + 1
    right = hz("C4") * 2 ** (-19 / 1200)
    fr = dict(level=[0.0] * frames, centroid=[1000.0] * frames, hnr=[20.0] * frames, flat=[-30.0] * frames,
              rough=[0.1] * frames, f0=[hz("C4") * 2 ** (61 / 1200)] * frames, cents=[61.0] * frames)
    n = dict(at=0.0, dur=1.0, note="C4", target_hz=hz("C4"), f0_hz=right, cents=-19.0)
    got = analyze.note_measures(x, fr, n, "voice", True)
    assert got["cents"] == -19.0 and "speaks_ms" not in got, got
    plain = analyze.note_measures(x, fr, n, "voice")
    assert plain["cents"] == 61.0 and "speaks_ms" in plain, plain  # the planted failure: the mistracked frames win
    tracked = analyze.note_measures(x, dict(fr, f0=[right] * frames), n, "voice")
    assert got["partials_db"] == tracked["partials_db"], (got, tracked)


for name, fn in sorted(globals().copy().items()):
    if name.startswith("test_") and callable(fn):
        fn()
        print("ok", name)
print("all passed")
