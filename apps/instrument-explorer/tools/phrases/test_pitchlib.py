"""Planted negatives and regression guards for pitchlib.f0, on synthesized signals only (no audio is read or committed).

    uv run -q --with numpy python3 -I apps/instrument-explorer/tools/phrases/test_pitchlib.py

To run the same cases against another copy of the detector (to see a case fail on an older file), name it:

    PITCHLIB=/path/to/pitchlib.py uv run -q --with numpy python3 -I apps/instrument-explorer/tools/phrases/test_pitchlib.py

Every signal goes through f0(path) as a 16-bit mono WAV 2 s long, so the steady middle half is 1 s, as for a site render.
Phases and noise are seeded, so a run is repeatable. Each case states which kind it is:
  planted negative: reads wrong on the file before the own_ladder change (the bowed bass open strings, 275.8 Hz for a 55 Hz A1);
  guard: reads right before and after, and keeps the neighbouring reading rules from being relaxed to fix the negative.
"""
import importlib.util
import os
import tempfile
import unittest
import wave

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("pitchlib_under_test", os.environ.get("PITCHLIB") or os.path.join(HERE, "pitchlib.py"))
pl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pl)

SR = 48000

# Levels in dB of harmonics 1..12 of a bowed double bass open string, 1 s middle of the site's own render (fundamental
# 39 dB under one body-boosted harmonic, the lines beside it 11-19 dB under it, a line at every slot). Beyond h12 the
# lines fall 7 dB per octave of harmonic number, as in the render.
BASS_A1_DB = [-39.1, -32.0, -26.2, -19.1, 0.0, -15.8, -21.6, -24.8, -26.9, -28.6, -35.0, -31.1]  # 55 Hz, boosted h5
BASS_E1_DB = [-39.5, -33.0, -28.3, -24.0, -19.0, -11.2, 0.0, -13.4, -18.2, -21.1, -27.6, -24.6]  # 41.2 Hz, boosted h7


def ladder(f0, levels_db, top=45):
    """{frequency: dB}: harmonics of f0 at the given levels, then -7 dB per octave of harmonic number up to harmonic `top`."""
    lines = {f0 * (h + 1): db for h, db in enumerate(levels_db)}
    for h in range(len(levels_db) + 1, top + 1):
        lines[f0 * h] = levels_db[-1] - 7.0 * np.log2(h / len(levels_db))
    return lines


def synth(lines, seconds=2.0, noise_db=-75.0, seed=1):
    """Sum of sines {frequency: dB} with seeded phases and seeded white noise noise_db under the 0 dB line, peak 0.5."""
    rng = np.random.default_rng(seed)
    t = np.arange(int(seconds * SR)) / SR
    x = np.zeros_like(t)
    for f, db in lines.items():
        if f < 0.45 * SR:
            x += 10 ** (db / 20) * np.sin(2 * np.pi * f * t + rng.uniform(0, 2 * np.pi))
    x += 10 ** (noise_db / 20) * rng.standard_normal(len(t)) / 3
    return 0.5 * x / np.abs(x).max()


def f0_of(x):
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "synth.wav")
        with wave.open(path, "wb") as w:
            w.setnchannels(1)
            w.setsampwidth(2)
            w.setframerate(SR)
            w.writeframes(np.round(x * 32767).astype("<i2").tobytes())
        return pl.f0(path)[0]


class PitchReader(unittest.TestCase):
    def assertReads(self, x, hz, tol=0.01):
        got = f0_of(x)
        self.assertLess(abs(got / hz - 1), tol, f"read {got:.2f} Hz, expected {hz} Hz +-{tol:.0%}")

    def test_bass_a1_fundamental_39_db_under_the_boosted_5th_reads_55_not_275(self):
        """PLANTED NEGATIVE. 55 Hz with h1 at -39 dB and h5 (275 Hz) boosted: the autocorrelation peaks at 55 Hz, the
        55 Hz candidate has the best contrast, and subharmonic() used to reject it as a weak series under 275 Hz (the
        mean of its seven off-slot levels is 25 dB under h5) so the reader answered 275.86 Hz."""
        self.assertReads(synth(ladder(55.0, BASS_A1_DB)), 55.0)

    def test_bass_e1_fundamental_under_the_boosted_7th_reads_41_2(self):
        """Guard (reads right on the file before the change too; the 7th-harmonic case of the site's E1 render is the
        287.9 Hz misreading, which needs that render's own tail to reproduce)."""
        self.assertReads(synth(ladder(41.2, BASS_E1_DB)), 41.2)

    def test_series_of_even_multiples_of_116_94_reads_233_88_not_116_94(self):
        """Guard. Lines only at 233.88 x m (the even multiples of 116.94), falling 1/m: 116.94 is a candidate whose
        odd slots are all empty, and must stay rejected whatever is done to the weak-series test."""
        self.assertReads(synth({233.88 * m: -20 * np.log10(m) for m in range(1, 9)}), 233.88)

    def test_clean_463_8_hz_tone_reads_463_8(self):
        """Guard. 463.8 Hz at 48 kHz is a period of 103.5 samples, so the autocorrelation is higher at its octave below
        (207 whole samples) than at the true lag, and the spectrum has to decide."""
        self.assertReads(synth({463.8: 0.0}), 463.8)

    def test_a_real_weak_series_an_octave_below_a_tone_is_still_rejected(self):
        """Guard. A 123.47 Hz tone (lines falling 6 dB per octave) carrying its own series at the odd multiples of
        61.735 Hz, 26 dB then 39 dB down: every slot of that series is a real line, but the tone has a comb of its own
        (each of its lines 20 dB or more over the two beside it), so it is the weak series, not the ladder, and the
        reader must stay on 123.47."""
        for top in (-26.0, -32.0):
            lines = {123.47 * m: -6.0 * np.log2(m) for m in range(1, 9)}
            for i, h in enumerate((1, 3, 5, 7, 9, 11)):
                lines[61.735 * h] = top - 13.0 * i / 5
            with self.subTest(sub_series_top_db=top):
                self.assertReads(synth(lines), 123.47)

    def test_one_stray_line_among_a_weak_series_does_not_make_it_a_ladder(self):
        """Guard. The same 123.47 Hz tone with its sub-series 30 dB down (falling 2 dB per line) and ONE stray line at
        -12 dB on one of its odd slots, on the grid and 0.4 % off it (a breath or noise peak): the stray stands over its
        neighbours but the tone's strongest line still has one slot beside it 20 dB or more under it, so the reader
        stays on 123.47. A statistic that took the strongest off slot (mean(on) - max(off)) fails these."""
        for slot in (1, 3, 5, 7, 9):
            for off_pct in (0.0, 0.4):
                lines = {123.47 * m: -6.0 * np.log2(m) for m in range(1, 9)}
                for i, h in enumerate((1, 3, 5, 7, 9, 11)):
                    lines[61.735 * h] = -30.0 - 2.0 * i
                lines[61.735 * slot * (1 + off_pct / 100)] = -12.0
                with self.subTest(stray_slot=slot, stray_offset_pct=off_pct):
                    self.assertReads(synth(lines), 123.47)


if __name__ == "__main__":
    unittest.main(verbosity=2)
