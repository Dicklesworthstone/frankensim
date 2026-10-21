"""formants.py on synthetic vowels of known formants.

uv run --with numpy python3 -I test_formants.py -v
The sources are radiated glottal-like pulse trains (harmonics falling 6 dB/octave) through three two-pole resonators at the
formant frequencies and bandwidths below, so every expected value is set by the filter, not by a recording.
"""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import formants as fm
from formants import SR


def vowel(f0, formants_hz, bws=(70.0, 90.0, 140.0), dur=2.4):
    n = int(dur * SR)
    t = np.arange(n) / SR
    x = np.zeros(n)
    for k in range(1, int(8000 // f0) + 1):
        f = k * f0
        h = 1.0 / k  # glottal source -12 dB/octave times lip radiation +6 dB/octave: the radiated pressure falls 6 dB/octave
        for fc, bw in zip(formants_hz, bws):  # |H| of a two-pole resonator at the harmonic, unit gain at DC removed
            h *= 1.0 / np.sqrt(((1 - (f / fc) ** 2) ** 2) + (f * bw / fc ** 2) ** 2)
        x += h * np.sin(2 * np.pi * f * t + 0.7 * k)
    return 0.3 * x / np.abs(x).max()


class FormantTests(unittest.TestCase):
    def test_recovers_known_formants_across_the_pitch_range(self):
        for f0 in (130.8, 155.6, 220.0, 293.7, 349.2):
            for target in ((700.0, 1200.0, 2600.0), (300.0, 2200.0, 3000.0), (500.0, 900.0, 2400.0)):
                fr, bw = fm.formants(vowel(f0, target), f0)
                err = np.abs(fr - np.array(target))
                # half the harmonic spacing is the most a harmonic-sampled estimate can promise; at 349 Hz the harmonics
                # (698, 1047, 1396) are too sparse to separate resonances 500 Hz apart, so there the bound is 0.75 spacings
                tol = max(70.0, (0.5 if f0 < 300.0 else 0.75) * f0)
                self.assertTrue(np.all(err < tol), f"f0 {f0} target {target} got {fr.round()} err {err.round()}")

    def test_a_harmonic_is_not_mistaken_for_a_resonance(self):
        # 155.5 Hz, /a/-like F2 at 1180 Hz: the 6th harmonic is at 933 Hz and must not be reported as a formant
        fr, _ = fm.formants(vowel(155.5, (722.0, 1180.0, 2690.0)), 155.5)
        self.assertGreater(fr[1], 1050.0)
        self.assertLess(fr[1], 1300.0)

    def test_story_table_parses_the_measured_means(self):
        p = Path(tempfile.mkdtemp()) / "t.tsv"
        p.write_text("# comment\n\ti\tɑ\tu\nFo\t159\t154\t155\nF1N Sd\t295±17\t715±9\t430±5\nF2N Sd\t1923±28\t1180±31\t917±14\n"
                     "F3N Sd\t2797±60\t2691±21\t2087±14\nF1C\t325\t694\t314\n")
        t = fm.story_table(p)
        self.assertEqual(t["i"], (295.0, 1923.0, 2797.0, 159.0))
        self.assertEqual(t["a"], (715.0, 1180.0, 2691.0, 154.0))
        self.assertNotIn("e", t)

    def test_harmonic_summary_tells_a_one_peak_spectrum_from_a_three_formant_one(self):
        three = fm.harmonic_summary(vowel(155.5, (700.0, 2000.0, 3000.0)), 155.5)
        one = fm.harmonic_summary(vowel(155.5, (700.0, 7000.0, 9000.0)), 155.5)  # F2 and F3 far above the band
        self.assertAlmostEqual(three["peak_hz"], one["peak_hz"], delta=160.0)
        self.assertGreater(three["mid_db"], one["mid_db"] + 8.0)
        self.assertGreater(three["hi_db"], one["hi_db"] + 8.0)

    def test_harmonic_summary_reads_the_strongest_low_harmonic_and_the_band_contrast(self):
        x = vowel(155.5, (722.0, 1180.0, 2690.0))
        s = fm.harmonic_summary(x, 155.5)
        self.assertTrue(600.0 < s["peak_hz"] < 1300.0)
        self.assertLess(s["hi_db"], s["mid_db"] + 1.0)  # the higher band is no stronger than the middle one

    def test_too_little_signal_refuses(self):
        with self.assertRaises(ValueError):
            fm.formants(np.zeros(SR // 10), 155.0)


if __name__ == "__main__":
    unittest.main()
