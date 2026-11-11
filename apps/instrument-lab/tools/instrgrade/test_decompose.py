"""decompose.py on synthetic notes with known answers.

uv run --with numpy --with matplotlib python3 -I test_decompose.py -v
Every expected value is fixed by how the signal is built (a delay, a detuning, a missing partial), not by a recording.
"""
from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import decompose as dc
from decompose import SR, Note

F0 = 220.0


def tone(f0=F0, dur=2.6, attack=0.04, lead=0.3, delay=0.0, odd_only=False, cents=0.0, noise_db=None, tilt=1.0, seed=1):
    """Sustained harmonic tone with a raised-cosine attack, `lead` seconds of leading silence, optional noise floor."""
    n = int(dur * SR)
    t = np.arange(n) / SR
    f = f0 * 2 ** (cents / 1200)
    x = np.zeros(n)
    for k in range(1, 25):
        if odd_only and k % 2 == 0:
            continue
        if k * f > 15000:
            break
        x += np.sin(2 * np.pi * k * f * t + 0.3 * k) / k ** tilt
    env = np.clip((t - lead - delay) / attack, 0, 1)
    x *= 0.5 - 0.5 * np.cos(np.pi * env)
    x *= 0.2 / np.abs(x).max()
    if noise_db is not None:
        x += np.random.default_rng(seed).standard_normal(n) * 10 ** (noise_db / 20)
    return x


def note_of(x, f0=F0, kind="sustained", max_cents=50.0):
    floor = dc.render_floor(dc.highpass(x, f0), f0)
    y, info = dc.align(x, f0, floor)
    return Note(y, f0, kind, False, "synthetic", info, max_cents=max_cents)


class AlignmentTests(unittest.TestCase):
    def test_delayed_copy_aligns_to_the_millisecond(self):
        a, b = note_of(tone()), note_of(tone(delay=0.137, noise_db=-80))
        self.assertLess(abs(a.info["cut_onset_s"] - b.info["cut_onset_s"]) * 1000, dc.ALIGN_TOL_MS)
        self.assertAlmostEqual(a.info["cut_onset_s"], dc.preroll_s(F0), delta=0.0015)
        self.assertGreater(abs(a.info["onset_s"] - b.info["onset_s"]), 0.1)  # the cut really moved them

    def test_legato_reattack_is_found_through_the_previous_notes_tail(self):
        n = int(2.0 * SR)
        t = np.arange(n) / SR
        first = np.where(t < 0.5, np.exp(-(t - 0.0) * 6), 0) * np.sin(2 * np.pi * 233 * t)  # decaying tail, still audible at 0.5 s
        env = np.clip((t - 0.5) / 0.03, 0, 1)
        second = 3.0 * env * np.sin(2 * np.pi * 247 * t)
        x = 0.1 * (first + second) / 3.0
        # the search starts 20 ms before the new note, inside the tail
        ta, tb = dc.onset_times(dc.highpass(x, 247.0), -100.0, anchor_s=0.52, f0_hz=247.0)
        self.assertAlmostEqual(max(ta, tb), 0.5, delta=0.01)

    def test_slow_ramp_from_rest_with_no_quiet_frames_aligns(self):
        # a bowed render: starts from silence at sample 0, ramps over 0.4 s, never has a quiet 20 ms frame afterwards
        n = note_of(tone(attack=0.4, lead=0.0, dur=3.0))
        self.assertAlmostEqual(n.info["cut_onset_s"], dc.preroll_s(F0), delta=0.0015)
        again = note_of(n.y)  # cutting at the onset and detecting again finds the same crossing
        self.assertAlmostEqual(again.info["cut_onset_s"], n.info["cut_onset_s"], delta=0.0015)

    def test_note_that_begins_on_a_loud_tail_is_refused_not_guessed(self):
        n = int(2.0 * SR)
        t = np.arange(n) / SR
        x = 0.1 * (np.sin(2 * np.pi * 233 * t) * np.where(t < 0.5, 1.0, 0) + 1.0 * np.clip((t - 0.5) / 0.03, 0, 1) * np.sin(2 * np.pi * 247 * t))
        with self.assertRaisesRegex(dc.OnsetMisaligned, "NO_QUIET_BEFORE_ONSET"):
            dc.onset_times(dc.highpass(x, 247.0), -100.0, anchor_s=0.52, f0_hz=247.0)

    def test_anchor_inside_a_slow_rise_finds_the_start_of_the_rise(self):
        x = tone(attack=0.2, lead=0.5)
        y = dc.highpass(x, F0)
        ta, tb = dc.onset_times(y, -100.0, anchor_s=0.9, f0_hz=F0)  # the anchor lags well inside the ramp
        first = dc.onset_times(y, -100.0, f0_hz=F0)
        self.assertAlmostEqual(max(ta, tb), max(first), delta=0.002)

    def test_a_recording_that_begins_at_its_own_onset_aligns_only_from_the_file_start(self):
        n = int(2.0 * SR)
        tt = np.arange(n) / SR
        x = 0.2 * np.exp(-tt * 2.0) * np.sin(2 * np.pi * 196 * tt)  # a plucked string: loudest at sample 0
        floor = -90.0
        cut, info = dc.align(x, 196.0, floor, anchor_s=0.02, at_start=True)
        self.assertAlmostEqual(info["cut_onset_s"], dc.preroll_s(196.0), delta=0.0015)
        self.assertLess(abs(info["onset_s"]), 0.01)  # in the time of the file: its first sample
        with self.assertRaisesRegex(dc.OnsetMisaligned, "NO_QUIET_BEFORE_ONSET"):
            dc.align(x, 196.0, floor, anchor_s=0.02, at_start=False)  # a window cut from the middle of a file: legato
        render, rinfo = dc.align(x, 196.0, floor)  # a render that starts at sample 0 is the same note
        self.assertAlmostEqual(rinfo["cut_onset_s"], dc.preroll_s(196.0), delta=0.0015)

    def test_silence_refuses(self):
        with self.assertRaises(dc.OnsetMisaligned):
            dc.onset_times(np.zeros(SR), -120.0)


class GapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = note_of(tone())

    def gaps(self, x):
        return dc.compare(note_of(x), self.ref)

    def test_identical_note_has_zero_gaps(self):
        g = dc.compare(self.ref, self.ref)
        for m in ("spec_db", "stft_lsd_db", "cwt_lsd_db", "cepstral_env_db", "partial_freq_cents", "f0_cents", "attack_log2",
                  "centroid_st", "flux_db", "odd_even_db"):
            self.assertAlmostEqual(g[m][0], 0.0, places=6, msg=m)

    def test_level_and_delay_do_not_register(self):
        g = self.gaps(0.25 * tone(delay=0.2, noise_db=-90))
        self.assertLess(g["stft_lsd_db"][0], 2.0)
        self.assertLess(g["attack_log2"][0], 0.1)
        self.assertLess(g["spec_db"][0], 1.0)

    def test_detuning_is_read_in_cents(self):
        g = self.gaps(tone(cents=25.0))
        self.assertAlmostEqual(g["f0_cents"][1], 25.0, delta=1.0)
        self.assertLess(g["partial_freq_cents"][0], 3.0)  # a pure detune moves f0, not the partial ratios

    def test_attack_ratio_is_log2(self):
        g = self.gaps(tone(attack=0.16))  # 4x slower
        self.assertAlmostEqual(g["attack_log2"][0], 2.0, delta=0.35)

    def test_odd_only_spectrum_is_far_in_odd_even_and_harmonic_levels(self):
        g = self.gaps(tone(odd_only=True))
        self.assertGreater(g["odd_even_db"][0], 6.0)
        self.assertGreater(g["harm_db"][0], 3.0)

    def test_darker_tone_has_lower_centroid(self):
        g = self.gaps(tone(tilt=2.0))
        self.assertLess(g["centroid_st"][1], -1.0)
        self.assertGreater(g["cepstral_env_db"][0], 1.0)

    def test_stretched_partials_show_in_partial_frequency(self):
        n = int(2.6 * SR)
        t = np.arange(n) / SR
        B = 4e-4
        x = sum(np.sin(2 * np.pi * k * F0 * math.sqrt(1 + B * k * k) * t) / k for k in range(1, 20))
        env = np.clip((t - 0.3) / 0.04, 0, 1)
        x = 0.2 * x * (0.5 - 0.5 * np.cos(np.pi * env)) / np.abs(x).max()
        g = self.gaps(x)
        self.assertGreater(g["inharmonicity_B"][0], 1e-4)
        self.assertGreater(g["partial_freq_cents"][0], 3.0)


class PitchTests(unittest.TestCase):
    def weak_fundamental(self, f0, cents=0.0):
        n = int(2.6 * SR)
        t = np.arange(n) / SR
        f = f0 * 2 ** (cents / 1200)
        amps = {1: 0.35, 2: 1.0, 3: 0.2, 4: 0.12, 5: 0.06}  # the second partial dominates, as in a bowed bass
        x = sum(a * np.sin(2 * np.pi * k * f * t + k) for k, a in amps.items())
        return 0.2 * x / np.abs(x).max() * (0.5 - 0.5 * np.cos(np.pi * np.clip((t - 0.3) / 0.05, 0, 1)))

    def test_comb_finds_the_fundamental_when_the_second_partial_dominates(self):
        for cents in (0.0, 4.0, -30.0, 120.0):
            x = self.weak_fundamental(55.0, cents)
            y = dc.highpass(x, 55.0)
            cut, _ = dc.align(y, 55.0, dc.render_floor(y, 55.0))
            f = dc.f0_comb(cut, 55.0)
            self.assertAlmostEqual(1200 * math.log2(f / (55.0 * 2 ** (cents / 1200))), 0.0, delta=1.0, msg=str(cents))

    def test_note_pitch_gap_is_the_true_detune_not_an_estimator_artifact(self):
        ref = note_of(self.weak_fundamental(55.0), 55.0)
        for cents in (4.0, 150.0):
            ours = note_of(self.weak_fundamental(55.0, cents), 55.0, max_cents=400.0)
            self.assertAlmostEqual(dc.compare(ours, ref)["f0_cents"][1], cents, delta=2.0)

    def test_a_tone_with_nothing_at_its_fundamental_is_refused(self):
        n = int(2.6 * SR)
        tt = np.arange(n) / SR
        env = 0.5 - 0.5 * np.cos(np.pi * np.clip((tt - 0.3) / 0.05, 0, 1))
        for f in (1558.2, 2047.3):  # the flute lock-on tones against a D4 / G4 name
            x = 0.2 * env * (np.sin(2 * np.pi * f * tt) + 0.2 * np.sin(2 * np.pi * 3 * f * tt))
            with self.assertRaisesRegex(dc.GraderRefusalError, "NOT_THE_NAMED_NOTE|PITCH_MISMATCH"):
                note_of(x, 293.66 if f < 2000 else 392.0, max_cents=400.0)

    def test_the_refusal_says_the_named_fundamental_and_the_graders_own_estimate_apart(self):
        # P8 5066: "nothing at the fundamental near 208.5 Hz" for an F3 (174.6 Hz) read as the note's pitch; 208.5 Hz was the estimate
        n = int(2.6 * SR)
        tt = np.arange(n) / SR
        env = 0.5 - 0.5 * np.cos(np.pi * np.clip((tt - 0.3) / 0.05, 0, 1))
        x = 0.2 * env * np.sin(2 * np.pi * 1200.0 * tt)  # a parasitic line, nothing near 174.6 or 208.5 Hz
        with self.assertRaises(dc.GraderRefusalError) as named:
            dc.check_fundamental(x, 208.5, "sustained", "ours", named_hz=174.6)
        msg = str(named.exception)
        self.assertIn("IRPDECOMP_NOT_THE_NAMED_NOTE", msg)
        self.assertIn("named note's fundamental is 174.6 Hz", msg)
        self.assertIn("pitch estimate for it is 208.5 Hz", msg)
        self.assertIn("strongest spectral line is 1200 Hz", msg)
        with self.assertRaises(dc.GraderRefusalError) as plain:  # the form without a name is unchanged
            dc.check_fundamental(x, 208.5, "sustained", "ours")
        self.assertIn("nothing at the fundamental near 208.5 Hz", str(plain.exception))

    def test_a_stretched_harmonic_note_passes_the_fundamental_check(self):
        n = int(2.6 * SR)
        tt = np.arange(n) / SR
        x = sum(np.sin(2 * np.pi * k * 110.0 * math.sqrt(1 + 6e-4 * k * k) * tt) / k for k in range(1, 16))
        x = 0.2 * x * (0.5 - 0.5 * np.cos(np.pi * np.clip((tt - 0.3) / 0.05, 0, 1))) / np.abs(x).max()
        note_of(x, 110.0)

    def test_a_reference_more_than_fifty_cents_off_its_name_is_refused(self):
        with self.assertRaisesRegex(dc.GraderRefusalError, "PITCH_MISMATCH"):
            note_of(self.weak_fundamental(55.0, 90.0), 55.0)


# the keys the page generator reads from a references.json entry: id, licence (a CREDITS key) and the url of the recording
REFERENCE_ENTRIES = {i: dict(id=i, licence="uiowa-mis", url=f"https://theremin.music.uiowa.edu/{i}.aiff") for i in ("syn-mf", "syn-ff")}


def cwt_db_full_rate(y, per_oct=24, f_lo=60.0, f_hi=dc.F_HI, w0=20.0):
    """The scalogram as first written: one full-length inverse FFT per scale (the reference the fast one must match)."""
    import math
    n = len(y)
    nfft = 1 << int(math.ceil(math.log2(n + SR // 4)))
    X = np.fft.rfft(y, nfft)
    f = np.fft.rfftfreq(nfft, 1 / SR)
    fc = f_lo * 2 ** (np.arange(int(per_oct * math.log2(f_hi / f_lo)) + 1) / per_oct)
    frames = n // dc.HOP
    out = np.empty((len(fc), frames), np.float32)
    for i, c in enumerate(fc):
        lo, hi = np.searchsorted(f, c * (1 - 6 / w0)), np.searchsorted(f, c * (1 + 6 / w0))
        Y = np.zeros(nfft, complex)
        Y[lo:hi] = X[lo:hi] * np.exp(-0.5 * ((f[lo:hi] - c) / (c / w0)) ** 2)
        z = np.fft.ifft(Y)[:n]
        out[i] = ((z.real ** 2 + z.imag ** 2)[: frames * dc.HOP].reshape(frames, dc.HOP)).mean(1)
    return (10 * np.log10(out + 1e-20)).astype(np.float32), fc


class ScalogramTests(unittest.TestCase):
    def test_fast_scalogram_matches_the_full_rate_one_where_it_matters(self):
        import time
        x = tone(f0=196.0, dur=2.2, lead=0.1, attack=0.05, noise_db=-70)[: int(2.0 * SR)]
        t = time.time()
        ref, fc = cwt_db_full_rate(x)
        t_ref = time.time() - t
        t = time.time()
        new, fc2 = dc.cwt_db(x)
        t_new = time.time() - t
        self.assertEqual(ref.shape, new.shape)
        live = ref > ref.max() - 60.0  # cells within 60 dB of the peak: the ones every metric and figure uses
        err = (new - ref)[live]
        self.assertLess(float(np.sqrt(np.mean(err ** 2))), 0.6, "rms dB difference over live cells")
        self.assertLess(float(np.abs(err).max()), 3.0)
        self.assertLess(t_new, t_ref)


class OnsetRuleTests(unittest.TestCase):
    """The onset is where the envelope first reaches 20 dB under the signal's OWN steady level, for any signal."""

    def step_then_growth(self, f0=233.08, noise_db=None, lead=2.0, seed=2):
        # digital silence, a priming step 26 dB under the final level held 70 ms, then growth over 0.25 s to a steady note
        n = int((lead + 3.0) * SR)
        tt = np.arange(n) / SR
        env = np.zeros(n)
        t0 = lead
        env[(tt >= t0) & (tt < t0 + 0.07)] = 0.05
        g = (tt >= t0 + 0.07) & (tt < t0 + 0.32)
        env[g] = 0.05 + 0.75 * ((tt[g] - t0 - 0.07) / 0.25) ** 2
        env[tt >= t0 + 0.32] = 0.8
        x = env * (np.sin(2 * np.pi * f0 * tt) + 0.4 * np.sin(2 * np.pi * 2 * f0 * tt))
        if noise_db is not None:
            x = x + np.random.default_rng(seed).standard_normal(n) * 10 ** (noise_db / 20)
        return x, t0

    def test_a_low_priming_step_does_not_count_as_the_onset(self):
        f0 = 233.08
        x, t0 = self.step_then_growth(f0)
        y = dc.highpass(x, f0)
        e, tms = dc.env_rms_db(y, win=dc.env_window(f0))
        old_rule_onset = float(tms[np.flatnonzero(e > e.max() - 30.0)[0]])  # the previous rule: peak - 30 dB
        self.assertLess(old_rule_onset, t0 + 0.07, "the old rule latches onto the priming step: the test discriminates")
        cut, info = dc.align(x, f0, dc.render_floor(y, f0))  # raises if the 1.5 ms re-detection fails
        self.assertGreater(info["onset_s"], t0 + 0.07)  # after the step, in the growth
        self.assertLess(info["onset_s"], t0 + 0.30)
        self.assertAlmostEqual(info["cut_onset_s"], dc.preroll_s(f0), delta=0.0015)

    def test_a_recording_of_the_same_note_gets_the_same_onset_through_the_reference_path(self):
        f0 = 233.08
        render, t0 = self.step_then_growth(f0)
        recording, _ = self.step_then_growth(f0, noise_db=-62.0)  # a room: noise floor 50+ dB under the steady level
        _, ri = dc.align(render, f0, dc.render_floor(dc.highpass(render, f0), f0))
        floor = dc.floor_of_file(recording, f0)
        _, qi = dc.align(recording, f0, floor, anchor_s=t0 + 0.3, at_start=False)  # as the library gives it: an anchor inside the note
        self.assertLess(abs(ri["onset_s"] - qi["onset_s"]), 0.003)  # the same instant, found by the same rule on both

    def pluck_with_early_touch(self, f0=233.08, noise_db=None, lead=2.0, seed=3):
        # a plucked note (a step to full level, then a 0.4 s decay) preceded by a 20 ms touch of the string 24 dB under it and a 20 ms gap:
        # the envelope crosses steady - 20 dB at the touch, drops back under it, and crosses again at the pluck
        n = int((lead + 3.0) * SR)
        tt = np.arange(n) / SR
        env = np.zeros(n)
        env[(tt >= lead - 0.040) & (tt < lead - 0.020)] = 0.05
        k = tt >= lead
        env[k] = 0.8 * np.exp(-(tt[k] - lead) / 0.4)
        x = env * np.sin(2 * np.pi * f0 * tt)
        if noise_db is not None:
            x = x + np.random.default_rng(seed).standard_normal(n) * 10 ** (noise_db / 20)
        return x, lead

    def test_an_early_touch_under_the_steady_level_does_not_make_a_reference_unalignable(self):
        # the failure of the regrade: guitar Gb4/F5 and violin B5 references were refused (cut onset 1.4-2.7 ms against 10 ms) because the
        # anchored detection took the last dip under the threshold before the note and the re-detection on the cut took the FIRST crossing
        f0 = 233.08
        recording, t0 = self.pluck_with_early_touch(f0, noise_db=-62.0)
        y = dc.highpass(recording, f0)
        e, tms = dc.env_rms_db(y, win=dc.env_window(f0))
        thr = dc.steady_level_db(e[int((t0 - 0.1) * 1000): int((t0 + 3.0) * 1000)]) - dc.STEADY_DROP_DB
        self.assertLess(float(tms[np.flatnonzero(e > thr)[0]]), t0 - 0.015, "the setup has an early crossing: the test discriminates")
        _, qi = dc.align(recording, f0, dc.floor_of_file(recording, f0), anchor_s=t0 + 0.3, at_start=False)  # raised ONSET_MISALIGNED (4.28 ms)
        self.assertAlmostEqual(qi["cut_onset_s"], dc.preroll_s(f0), delta=0.0015)
        # the pluck, not the touch: the touch ended 20 ms before the pluck and the envelope is back under the threshold from 15.7 ms before it;
        # the zero-phase high-pass rings ahead of a step edge, so the final rise starts a few ms before t0 (12.7 ms in the noise-free signal)
        self.assertGreater(qi["onset_s"], t0 - 0.0157)
        self.assertLess(qi["onset_s"], t0 + 0.005)

    def test_a_steep_onset_re_detects_after_the_cut_even_with_an_overshoot_and_noise(self):
        # A guard for the class of failure that refused a UIowa clarinet note (a steep 40 dB rise, a wandering steady level, noise
        # on the peak). This synthetic signal passes under the old rule too: the real recording (BbClar.mf.C4B4 @ Bb4, which is not
        # in the repository) is what reproduced it, refused with a cut onset of 8.14 ms before and aligned after.
        f0 = 466.16
        n = int(5.0 * SR)
        tt = np.arange(n) / SR
        env = np.clip((tt - 2.0) / 0.05, 0, 1) ** 2 * (1.0 + 0.35 * np.exp(-(tt - 2.0) / 0.08)) * (1.0 + 0.2 * np.sin(2 * np.pi * 0.7 * tt))  # overshoot, steady, a slow +-1.6 dB wander
        x = env * np.sin(2 * np.pi * f0 * tt)
        x = x + np.random.default_rng(5).standard_normal(n) * 10 ** (-55 / 20)
        floor = dc.floor_of_file(x, f0)
        for shift in (0.0, 0.0007, 0.0013, 0.0021):  # the window slides by a few ms, which moved a max-based threshold
            _, info = dc.align(x[int(shift * SR):], f0, floor, anchor_s=2.03 - shift, at_start=False)
            self.assertAlmostEqual(info["cut_onset_s"], dc.preroll_s(f0), delta=0.0015)


class FluxTests(unittest.TestCase):
    def test_flux_curve_is_what_the_old_inline_formula_gave_a_note(self):
        n = note_of(tone(noise_db=-55))  # noise puts live cells between the partials, where the flux lives
        S = n.S[:, : n.frames].astype(np.float64)
        m = (n.sf > dc.F_LO) & (n.sf < dc.F_HI)
        live = np.maximum(S[m][:, 1:], S[m][:, :-1]) > S[m].max() - dc.SHOW_DB
        pos = np.where(live, np.maximum(np.diff(S[m], axis=1), 0.0), 0.0)
        old = np.concatenate([[0.0], np.sqrt(np.sum(pos ** 2, axis=0) / np.maximum(live.sum(axis=0), 1))])
        np.testing.assert_array_equal(n.flux_db, old)
        np.testing.assert_array_equal(dc.flux_curve(S, n.sf), old)
        self.assertGreater(n.flux_db[n.sound].mean(), 0.5)  # not a vacuous equality of zeros


class FloorTests(unittest.TestCase):
    def test_block_wise_floor_matches_the_whole_file_floor_and_is_faster(self):
        import time
        rng = np.random.default_rng(3)
        parts = []
        for k in range(14):  # a range file: notes of 2.2 s with 0.6 s gaps of room noise
            parts += [tone(f0=F0 * 2 ** (k / 12), dur=2.2, lead=0.0, attack=0.05, noise_db=-70, seed=k), 10 ** (-72 / 20) * rng.standard_normal(int(0.6 * SR))]
        x = np.concatenate(parts)
        t = time.time()
        whole = dc.floor_of(dc.highpass(x, F0))
        t_whole = time.time() - t
        t = time.time()
        blocks = dc.floor_of_file(x, F0, block=1 << 17)
        t_blocks = time.time() - t
        self.assertLess(abs(whole - blocks), 2.0, f"{whole} vs {blocks}")
        self.assertLess(t_blocks, t_whole * 1.5)  # small files are not a speed test; the point is it never costs more

    @staticmethod
    def struck(f0=1046.5, dur=2.0, lead=0.0, dust_db=None, seed=5):
        """A struck note that starts at its first sample (a render) or after `lead` s of dust, with a piano's two-stage decay: 60 dB/s
        for 0.5 s (the prompt sound), then 6 dB/s. Most of its frames sit 30-40 dB under its peak, and its quietest 20 ms are its own
        decay, about 40 dB down: the shape of the 2 s C6 render (peak envelope -10.8 dB, steady -41.9 dB, quietest -53.3 dB)."""
        n = int((lead + dur) * SR)
        t = np.arange(n) / SR
        tt = np.maximum(t - lead, 0.0)
        level_db = np.where(tt < 0.5, -60.0 * tt, -30.0 - 6.0 * (tt - 0.5))
        x = sum(np.sin(2 * np.pi * k * f0 * tt) / k for k in range(1, 8)) * 10 ** (level_db / 20) * np.clip(tt / 0.002, 0, 1)
        x = np.where(t >= lead, 0.3 * x / np.abs(x).max(), 0.0)
        if dust_db is not None:
            x = x + np.random.default_rng(seed).standard_normal(n) * 10 ** (dust_db / 20)
        return x

    def test_a_decaying_render_takes_its_floor_from_its_pre_roll_not_from_its_own_decay(self):
        f0 = 1046.5  # C6, 2 s: the shape of the C6 piano render refused as IRPDECOMP_NO_SIGNAL in the 2026-10-21 reality check
        x = self.struck(f0)
        y = dc.highpass(x, f0)
        e, _ = dc.env_rms_db(y, win=dc.env_window(f0))
        steady, quietest = dc.steady_level_db(e), dc.floor_of(y)
        old = min(quietest, float(e.max()) - 50.0)
        print(f"\n  C6-like render: steady {steady:.1f} dB, quietest 20 ms {quietest:.1f} dB, old floor {old:.1f} dB, "
              f"new floor {dc.render_floor(y, f0):.1f} dB")
        self.assertLess(steady, old + 20.0, "the old rule refuses this render: the test discriminates")
        self.assertEqual(dc.render_floor(y, f0), -120.0)  # it begins at its own onset: its pre-roll is the silence align() adds
        cut, info = dc.align(x, f0, dc.render_floor(y, f0))
        self.assertAlmostEqual(info["onset_s"], 0.0, delta=0.002)
        self.assertAlmostEqual(info["cut_onset_s"], dc.preroll_s(f0), delta=0.0015)

    def test_a_render_with_a_pre_roll_takes_that_pre_roll_as_its_floor(self):
        f0 = 1046.5
        x = self.struck(f0, lead=0.5, dust_db=-90.0)  # 0.5 s of numerical dust at -90 dB before the strike
        fl = dc.render_floor(dc.highpass(x, f0), f0)
        print(f"\n  render with 0.5 s of -90 dB dust: floor {fl:.1f} dB")
        self.assertAlmostEqual(fl, -90.0, delta=3.0)  # the dust's level (high-passed white noise loses a little), not the decay's
        _, info = dc.align(x, f0, fl)
        self.assertAlmostEqual(info["onset_s"], 0.5, delta=0.002)


class PartialGateTests(unittest.TestCase):
    def test_noise_peaks_far_below_the_real_partials_do_not_make_inharmonicity_or_frequency_gaps(self):
        def two_partials(junk):
            n = int(2.6 * SR)
            tt = np.arange(n) / SR
            x = np.sin(2 * np.pi * 220 * tt) + 0.5 * np.sin(2 * np.pi * 440 * tt) + 0.3 * np.sin(2 * np.pi * 660 * tt)
            if junk:  # spectral lines nowhere near the harmonic series, 70 dB down: what a render's numerical floor looks like
                for f in (890.0, 1118.0, 1345.0, 1565.0, 1790.0, 2010.0):  # just off harmonics 4-9: inside the tracker's search windows
                    x = x + 3e-4 * np.sin(2 * np.pi * f * tt)
            x = 0.2 * x / np.abs(x).max() * (0.5 - 0.5 * np.cos(np.pi * np.clip((tt - 0.3) / 0.05, 0, 1)))
            return x
        clean, noisy = note_of(two_partials(False)), note_of(two_partials(True))
        old_gate, dc.STRONG_DB = dc.STRONG_DB, 1000.0  # the gate switched off: every tracked peak counts
        try:
            ungated = dc.compare(note_of(two_partials(True)), note_of(two_partials(False)))
            self.assertGreater(ungated["partial_freq_cents"][0], 8.0)  # the artifact is real without the gate
        finally:
            dc.STRONG_DB = old_gate
        g = dc.compare(noisy, clean)
        self.assertLess(g["partial_freq_cents"][0], 5.0)
        self.assertLess(g["inharmonicity_B"][0], 2e-5)
        self.assertLess(noisy.B, 2e-5)
        self.assertEqual(int(noisy.strong.sum()), 3)  # exactly the three real partials count


class GlideTests(unittest.TestCase):
    def test_a_note_whose_partials_are_all_offset_together_has_no_stiffness(self):
        n = int(2.6 * SR)
        tt = np.arange(n) / SR
        # f0 estimated 0.6% below where every partial really is (an unsettled pitch): the series is still exactly harmonic
        x = sum(a * np.sin(2 * np.pi * k * 221.3 * tt) for k, a in ((1, 1.0), (2, 0.6), (3, 0.3)))
        x = 0.2 * x / np.abs(x).max() * (0.5 - 0.5 * np.cos(np.pi * np.clip((tt - 0.3) / 0.05, 0, 1)))
        from unittest import mock
        with mock.patch.object(dc, "fit_stiff_string_law", wraps=dc.fit_stiff_string_law) as fit:
            note = note_of(x, 220.0)
        # the B fit is the one Note makes itself, on strong partials, with f0 free (describe()'s own calls pass f0_fixed)
        own = [c for c in fit.call_args_list if "f0_fixed" not in c.kwargs]
        self.assertTrue(own, "Note fits B with f0 free")
        self.assertLess(note.B, 1e-5)


class ListeningOnlyTests(unittest.TestCase):
    """The listening-only pair obeys every rule P8's generator enforces for a render with no reference."""

    def make(self, root, **kw):
        from instrgrade import write_wav
        wav = Path(tempfile.mkdtemp()) / "r.wav"
        write_wav(wav, 0.3 * tone(f0=523.25, dur=2.4, lead=0.2, noise_db=-80))
        c = dict(family="pan-pipes", note="C5", dyn="mf", wav=str(wav), git_sha="abc", observer="receiver 1 m", recipe="r.sh",
                 reason="No licensed recording of this instrument exists in the library.", listening_only=True)
        c.update(kw)
        dc._REFS = None
        return dc._batch_job((c, root)), c

    def test_the_pair_has_exactly_the_shape_p8_accepts(self):
        import json
        from instrgrade import load
        root = Path(tempfile.mkdtemp())
        res, c = self.make(root)
        self.assertTrue(res["ok"], res)
        d = root / "pan-pipes" / "C5"
        pj = json.loads((d / "pair.json").read_text())
        self.assertIsNone(pj["real"])
        self.assertTrue(pj["no_reference"]["reason"].strip())
        self.assertIsNone(pj["alignment"])
        self.assertIsNone(pj["partials"])
        self.assertEqual(pj["metrics"], [])
        self.assertIsNone(pj["audio"]["real"])
        self.assertEqual(pj["audio"]["ours"], "ours.wav")
        for kind in ("spectrogram", "scalogram"):
            f = pj["figures"][kind]
            self.assertIsNone(f["real"])
            self.assertIsNone(f["diff"])
            self.assertTrue((d / f["ours"]).stat().st_size > 1000)
            self.assertEqual((f["f_lo_hz"], f["f_hi_hz"], f["span_db"]), ((40.0 if kind == "spectrogram" else 60.0), dc.F_HI, 80.0))
        x = load(d / "ours.wav")
        self.assertAlmostEqual(len(x) / SR, pj["duration_s"], delta=0.05)
        self.assertLess(float(np.abs(x).max()), 0.999)
        self.assertGreater(float(np.abs(x).max()), 0.5)
        self.assertEqual(pj["ours"]["git_sha"], "abc")
        self.assertEqual(sorted(p.name for p in d.iterdir()), sorted(["pair.json", "ours.wav", "spectrogram-ours.webp", "scalogram-ours.webp"]))

    def test_calibrated_listening_retains_the_material_pressure_ratio(self):
        from instrgrade import load
        root = Path(tempfile.mkdtemp())
        signal = tone(f0=523.25, dur=1.2, lead=0.05)
        meta = dict(wav="native.wav", git_sha="fixture", observer="receiver 1 m", recipe="material_strings", sample_rate=SR, note="C5")
        for label, gain in (("reference", 1.0), ("changed", 0.25)):
            pj = dc.write_listening(root, "violin", label, signal * gain, {"onset_s": 0.05}, meta,
                                    "Manufactured material comparison; no reference recording.", full_scale_pa=2.0)
            self.assertEqual(pj["ours"]["full_scale_pa"], 2.0)
            self.assertEqual(pj["audio"]["gain_db"], 0.0)
            self.assertIn("without peak normalization", pj["audio"]["level_match"])
        a = load(root / "violin" / "reference" / "ours.wav")
        b = load(root / "violin" / "changed" / "ours.wav")
        np.testing.assert_allclose(b, 0.25 * a, rtol=1e-7, atol=1e-9)
        self.assertAlmostEqual(float(np.sqrt(np.mean(b * b) / np.mean(a * a))), 0.25, places=7)

    def test_invalid_listening_pressure_scales_refuse_before_staging(self):
        root = Path(tempfile.mkdtemp())
        for scale in (0.0, -1.0, float("nan"), float("inf")):
            with self.subTest(scale=scale), self.assertRaisesRegex(ValueError, "LISTEN_SCALE"):
                dc.write_listening(root, "violin", "G3", np.zeros(100), {}, {}, "No paired recording", full_scale_pa=scale)
        self.assertEqual(list(root.iterdir()), [])

    def test_a_second_render_moves_the_first_aside_and_a_missing_reason_is_a_refusal_record(self):
        import json
        root = Path(tempfile.mkdtemp())
        self.assertTrue(self.make(root)[0]["ok"])
        self.assertTrue(self.make(root, git_sha="def")[0]["ok"])
        d = root / "pan-pipes" / "C5"
        self.assertEqual(len(list((d / "superseded").iterdir())), 1)
        self.assertEqual(json.loads((d / "pair.json").read_text())["ours"]["git_sha"], "def")
        res, _ = self.make(root, reason="  ", git_sha="ghi")
        self.assertFalse(res["ok"])
        self.assertIn("LISTEN_REASON", res["error"])
        self.assertIn("LISTEN_REASON", json.loads((d / "pair.json").read_text())["refused"]["reason"])

    def test_a_band_cap_is_stated_and_cuts_the_figures(self):
        import json
        root = Path(tempfile.mkdtemp())
        res, _ = self.make(root, band_hz=3900.0)
        self.assertTrue(res["ok"], res)
        pj = json.loads((root / "pan-pipes" / "C5" / "pair.json").read_text())
        self.assertEqual(pj["figures"]["spectrogram"]["f_hi_hz"], 3900.0)
        self.assertIn("low-passed at 3900 Hz", pj["audio"]["level_match"])


class StrictRulerTests(unittest.TestCase):
    def rows(self, gap_spec, ruler):
        real, ours = note_of(tone()), note_of(tone(cents=5))
        gaps = {m: (1.0, None) for m in dc.METRICS}
        gaps["spec_db"] = (gap_spec, None)
        gaps["f0_cents"] = (gap_spec, None)
        full = {c: [] for c in ("other_dynamic", "other_take", "neighbour", "other_articulation", "other_instrument")}
        full.update(ruler)
        return {m["id"]: m for m in dc.metric_rows(gaps, real, ours, full)}

    def test_a_neighbouring_note_widens_only_the_loose_ruler(self):
        ruler = dict(other_dynamic=[dict(label="ff", with_="x", gaps=dict(spec_db=2.0, f0_cents=2.0))],
                     neighbour=[dict(label="next note", with_="y", gaps=dict(spec_db=10.0, f0_cents=10.0))])
        r = self.rows(5.0, ruler)
        s = r["spec_db"]
        self.assertEqual([x["value"] for x in s["ruler"]], [2.0, 10.0])
        self.assertEqual([x["value"] for x in s["ruler_strict"]], [2.0])
        self.assertTrue(s["inside_loose"])  # 5 <= 10: flattering, because a neighbouring note is in the ruler
        self.assertFalse(s["inside_strict"])  # 5 > 2: further from the recording than the same note's other dynamic is from it
        f = r["f0_cents"]  # pitch-absolute measures never had neighbours in the ruler: strict equals loose
        self.assertEqual(f["ruler"], f["ruler_strict"])
        self.assertEqual(f["inside_loose"], f["inside_strict"])

    def test_no_same_note_recording_means_no_strict_claim_not_a_failure(self):
        ruler = dict(neighbour=[dict(label="next note", with_="y", gaps=dict(spec_db=10.0))])
        s = self.rows(5.0, ruler)["spec_db"]
        self.assertEqual(s["ruler_strict"], [])
        self.assertTrue(s["inside_loose"])
        self.assertIsNone(s["inside_strict"])

    def test_a_reference_only_pair_has_the_strict_ruler_and_no_verdicts(self):
        real = note_of(tone())
        gaps = None
        ruler = dict(other_dynamic=[dict(label="ff", with_="x", gaps={m: 1.0 for m in dc.METRICS})], neighbour=[])
        full = {c: [] for c in ("other_dynamic", "other_take", "neighbour", "other_articulation", "other_instrument")}
        full.update(ruler)
        rows = dc.metric_rows(gaps, real, None, full)
        self.assertTrue(all(m["inside_loose"] is None and m["inside_strict"] is None for m in rows))
        self.assertTrue(any(m["ruler_strict"] for m in rows))

    def test_verdict_counts_add_up(self):
        ms = [dict(gap=1.0, inside_loose=True, inside_strict=True), dict(gap=2.0, inside_loose=True, inside_strict=False),
              dict(gap=3.0, inside_loose=False, inside_strict=False), dict(gap=4.0, inside_loose=True, inside_strict=None),
              dict(gap=None, inside_loose=None, inside_strict=None)]
        v = dc.verdict_counts(ms)
        self.assertEqual(v, {"measures": 4, "loose": {"inside": 3, "outside": 1, "no_ruler": 0}, "strict": {"inside": 1, "outside": 2, "no_ruler": 1}})


class ContractTests(unittest.TestCase):
    """write_pair against P8's validator (tools/build_instruments.py) on a synthetic library."""

    @classmethod
    def setUpClass(cls):
        import json
        cls.tmp = Path(tempfile.mkdtemp())
        wav = cls.tmp / "wav"
        wav.mkdir()
        # a UIowa-style file: three notes, each preceded by silence, plus a quieter second dynamic of the middle note
        notes = [(57, F0), (58, F0 * 2 ** (1 / 12)), (59, F0 * 2 ** (2 / 12))]

        def scale(scale_gain):
            out = [np.zeros(SR)]
            for _, f in notes:
                out += [tone(f0=f, dur=2.0, lead=0.0, attack=0.04) * scale_gain, np.zeros(SR)]
            return np.concatenate(out)
        from instrgrade import write_wav
        write_wav(wav / "Syn.mf.A3B3.wav", scale(1.0))
        write_wav(wav / "Syn.ff.A3B3.wav", 2.0 * scale(0.5))
        lib = []
        for fname, dyn, g in (("Syn.mf.A3B3.wav", "mf", 1.0), ("Syn.ff.A3B3.wav", "ff", 1.0)):
            for k, (m, f) in enumerate(notes):
                t0 = 1.0 + 3.0 * k
                lib.append(dict(file=fname, start=t0 - 0.005, end=t0 + 2.2, note=m, name=dc.name_of(m), kind="sustained", inst="syn", art="",
                                dyn=dyn, sul="", mallet="", range="A3B3", desc={}))
        (cls.tmp / "lib.json").write_text(json.dumps(lib))
        (cls.tmp / "refs.json").write_text(json.dumps({"files": [
            dict(path="x/Syn.mf.A3B3.aiff", page="p.html", chain="c", id="syn-mf"), dict(path="x/Syn.ff.A3B3.aiff", page="p.html", chain="c", id="syn-ff")]}))
        cls.refs = dc.Refs(cls.tmp / "lib.json", wav, cls.tmp / "refs.json")
        cls.render = cls.tmp / "render.wav"
        write_wav(cls.render, 0.3 * tone(f0=F0 * 2 ** (1 / 12), dur=2.4, lead=0.2, noise_db=-85, cents=8))

    def validate(self, root, family, note):
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
        import build_instruments as bi
        import json
        base = Path(root) / family / note
        pj = json.loads((base / "pair.json").read_text())
        problems = bi.validate_pair(pj, family, note, base, REFERENCE_ENTRIES) + bi.check_audio(base, pj)
        self.assertEqual(problems, [])
        return pj

    def test_reference_only_pair_meets_the_contract(self):
        root = self.tmp / "ref"
        ref_e = self.refs.find("syn", 58, "mf")[0]
        y, info = self.refs.aligned(ref_e)
        pj = dc.write_pair(root, "syn", "Bb3", ref_e, self.refs, y, info)
        self.validate(root, "syn", "Bb3")
        self.assertIsNone(pj["ours"])
        self.assertTrue(any(m["ruler"] for m in pj["metrics"]))
        self.assertTrue(all(m["gap"] is None for m in pj["metrics"]))

    def test_render_pair_meets_the_contract_and_reads_the_detune(self):
        root = self.tmp / "pair"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="abc", observer="test", recipe="r.sh")
        pj = dc.run_pair(self.refs, c, root)
        self.validate(root, "syn", "Bb3")
        by = {m["id"]: m for m in pj["metrics"]}
        self.assertAlmostEqual(by["f0_cents"]["gap"], 8.0, delta=1.5)
        self.assertAlmostEqual(pj["alignment"]["level_offset_db"], 20 * math.log10(1 / 0.3), delta=1.0)
        self.assertEqual(by["f0_cents"]["gap_unit"], "cents")
        self.assertNotIn("neighbour", " ".join(r["label"] for r in by["f0_cents"]["ruler"]))  # pitch-absolute: same note only
        self.assertTrue(any("neighbouring" in r["label"] for r in by["attack_log2"]["ruler"]))

    def test_a_render_through_its_sessions_chain_says_so_in_the_pair(self):
        """With a chain (chain.py) the render is graded as the session would have recorded it, and pair.json records which
        chain, from which floor (and the hashes of its files), with which flags, at which level match. P8 still accepts it."""
        import chain as chn
        floor = np.random.default_rng(9).standard_normal(20 * chn.FS_REC) * 10 ** (-60 / 20)
        ch = chn.Chain(id="syn-chain@p.html", page="p.html", record={}, stretches=[chn.Stretch("amb", "a" * 64, 0, len(floor), "ambient")],
                       flags=[dict(flag="capsule_not_applied")], bits=16, sources=[dict(path="amb", sha256="a" * 64, role="ambient")],
                       _audio={"amb": floor})
        root = self.tmp / "chained"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="abc", observer="test", recipe="r.sh", chain=ch)
        pj = dc.run_pair(self.refs, c, root)
        self.validate(root, "syn", "Bb3")
        rc = pj["chain"]
        print(f"\n  chain {rc['chain_id']}: floor {rc['floor_source']}, ambient {rc['ambient_sha256']}, {rc['bits']} bit at {rc['fs_hz']} Hz, "
              f"level match {rc['level_match_db']:+.2f} dB, flags {[f['flag'] for f in rc['flags']]}")
        self.assertEqual((rc["chain_id"], rc["floor_source"], rc["ambient_sha256"]), ("syn-chain@p.html", "ambient", ["a" * 64]))
        self.assertAlmostEqual(rc["level_match_db"], 20 * math.log10(1 / 0.3), delta=1.0)  # the render was 0.3 of the library's level
        plain = dc.run_pair(self.refs, dict(c, chain=None), self.tmp / "plain")
        self.assertNotIn("chain", plain)  # no chain unless one is asked for

    def test_a_chain_for_a_session_with_no_record_refuses_the_pair(self):
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="abc", observer="test", recipe="r.sh",
                 chain_refdir=str(self.tmp))
        with self.assertRaises(ValueError) as e:
            dc.run_pair(self.refs, c, self.tmp / "nochain")
        print(f"\n  {e.exception}")
        self.assertIn("IRPCHAIN_UNKNOWN_SESSION", str(e.exception))
        self.assertIsInstance(e.exception, dc.REFUSALS)  # recorded as a refusal of the pair, never graded without its chain

    def test_render_pair_refuses_misaligned_onsets(self):
        ref_e = self.refs.find("syn", 58, "mf")[0]
        y, info = self.refs.aligned(ref_e)
        raw = tone(f0=F0 * 2 ** (1 / 12), dur=2.4, lead=0.2)
        cy, ci = dc.align(raw, F0 * 2 ** (1 / 12), dc.render_floor(raw, F0 * 2 ** (1 / 12)))
        # the render's onset sits 10 ms late in its own cut: planted on the grid onset and on the sub-frame onset the pair check compares
        ci = dict(ci, cut_onset_s=ci["cut_onset_s"] + 0.01, cut_onset_subframe_s=ci["cut_onset_subframe_s"] + 0.01)
        meta = dict(wav="r.wav", git_sha="x", observer="o", sample_rate=48000, recipe="r.sh")
        with self.assertRaises(dc.OnsetMisaligned):
            dc.write_pair(self.tmp / "bad", "syn", "Bb3", ref_e, self.refs, y, info, cy, ci, meta)
        # an info made by an earlier grader has no sub-frame onset: the grid onset decides on that side, as before
        old = {k: v for k, v in ci.items() if k != "cut_onset_subframe_s"}
        with self.assertRaises(dc.OnsetMisaligned):
            dc.write_pair(self.tmp / "bad-old-info", "syn", "Bb3", ref_e, self.refs, y, info, cy, old, meta)

    def test_reference_only_run_never_replaces_a_graded_pair(self):
        root = self.tmp / "keep"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="abc", observer="test", recipe="r.sh")
        dc.run_pair(self.refs, c, root)
        self.assertTrue(dc.has_render(root, "syn", "Bb3"))
        self.assertFalse(dc.has_render(root, "syn", "A3"))
        dc._REFS = self.refs
        res = dc._ref_job((dict(family="syn", art="", dyn="mf", notes=["Bb3"], review=False), root))
        self.assertIn("untouched", res["kept"])
        self.assertIsNotNone(__import__("json").loads((root / "syn" / "Bb3" / "pair.json").read_text())["ours"])

    def test_reference_only_run_leaves_a_refusal_record_of_a_render_alone(self):
        root = self.tmp / "keeprefusal"
        dc.refuse_pair(root, "syn", "Bb3", "IRPDECOMP_NOT_THE_NAMED_NOTE: test", render_sha="zz", render_wav="r.wav", dyn="mf")
        self.assertTrue(dc.has_render(root, "syn", "Bb3"))
        dc.refuse_pair(root, "syn", "A3", "stale pair superseded")  # no render_sha: it names no render
        self.assertFalse(dc.has_render(root, "syn", "A3"))

    def test_a_restarted_reference_run_skips_what_it_already_regenerated(self):
        import time
        root = self.tmp / "resume"
        dc._REFS = self.refs
        ref_e = self.refs.find("syn", 58, "mf")[0]
        y, info = self.refs.aligned(ref_e)
        dc.write_pair(root, "syn", "Bb3", ref_e, self.refs, y, info)
        job = dict(family="syn", art="", dyn="mf", notes=["Bb3"], review=False)
        self.assertEqual(dc._ref_job((dict(job, since=time.time() - 60), root))["kept"], "already regenerated since --since")
        res = dc._ref_job((dict(job, since=time.time() + 60), root))  # an older pair than --since is regenerated
        self.assertNotIn("kept", res)
        self.assertEqual(len(list((root / "syn" / "Bb3" / "superseded").iterdir())), 1)

    def test_reference_only_run_falls_back_to_the_next_verifying_note(self):
        root = self.tmp / "fallback"
        dc._REFS = self.refs
        res = dc._ref_job((dict(family="syn", art="", dyn="mf", notes=["C9", "Bb3"], review=False), root))
        self.assertTrue(res["ok"])
        self.assertEqual(res["note"], "Bb3")
        self.assertEqual([r["note"] for r in res["refused"]], ["C9"])

    def test_refusal_supersedes_by_moving_and_never_deletes(self):
        import json
        root = self.tmp / "refuse"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="abc", observer="test", recipe="r.sh")
        dc.run_pair(self.refs, c, root)
        before = {p.name for p in (root / "syn" / "Bb3").iterdir()}
        self.assertIn("ours.wav", before)
        dc.refuse_pair(root, "syn", "Bb3", "IRPDECOMP_PITCH_MISMATCH: test", "abc", self.render, "mf")
        d = root / "syn" / "Bb3"
        rec = json.loads((d / "pair.json").read_text())
        self.assertIsNone(rec["ours"])
        self.assertEqual(rec["refused"]["render_sha"], "abc")
        self.assertIn("PITCH_MISMATCH", rec["refused"]["reason"])
        stamps = list((d / "superseded").iterdir())
        self.assertEqual(len(stamps), 1)
        self.assertEqual({p.name for p in stamps[0].iterdir()}, before)  # everything moved, including the old pair.json
        self.assertEqual({p.name for p in d.iterdir()}, {"pair.json", "superseded"})
        dc.refuse_pair(root, "syn", "Bb3", "again")  # a second refusal moves only the first record, keeping the first stamp
        self.assertEqual(len(list((d / "superseded").iterdir())), 2)

    def test_a_regrade_moves_the_earlier_files_aside_instead_of_overwriting(self):
        root = self.tmp / "regrade"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="abc", observer="test", recipe="r.sh")
        dc.run_pair(self.refs, c, root)
        first = (root / "syn" / "Bb3" / "ours.wav").read_bytes()
        dc.run_pair(self.refs, dict(c, git_sha="def"), root)
        d = root / "syn" / "Bb3"
        stamps = list((d / "superseded").iterdir())
        self.assertEqual(len(stamps), 1)
        self.assertEqual((stamps[0] / "ours.wav").read_bytes(), first)
        self.assertTrue((stamps[0] / "pair.json").exists())
        self.assertEqual(__import__("json").loads((d / "pair.json").read_text())["ours"]["git_sha"], "def")

    def test_a_recording_outside_references_json_carries_its_own_source(self):
        import json
        refs = dc.Refs(self.tmp / "lib.json", self.tmp / "wav", self.tmp / "refs.json")
        refs.external = dict(name="Set", url="https://example.org/set", licence="CC BY 4.0", licence_url="https://example.org/by", credit="Someone")
        root = self.tmp / "external"
        ref_e = refs.find("syn", 58, "mf")[0]
        y, info = refs.aligned(ref_e)
        pj = dc.write_pair(root, "voice", "Bb3", ref_e, refs, y, info)
        self.assertEqual(pj["real"]["source"]["credit"], "Someone")
        self.assertEqual(pj["real"]["source"]["id"], "syn-mf")
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
        import build_instruments as bi
        base = root / "voice" / "Bb3"
        self.assertEqual(bi.validate_pair(json.loads((base / "pair.json").read_text()), "voice", "Bb3", base, {}), [])

    def test_a_render_that_locks_onto_another_tone_is_refused_for_that_not_for_its_onset(self):
        import json
        root = self.tmp / "lock"
        n = int(2.4 * SR)
        tt = np.arange(n) / SR
        x = 0.2 * np.sin(2 * np.pi * 1915.8 * tt) * np.clip(tt / 0.01, 0, 1)  # sounds at once, on the wrong tone
        from instrgrade import write_wav
        wav = self.tmp / "lock.wav"
        write_wav(wav, x)
        res = dc._batch_job((dict(family="syn", note="Bb3", dyn="mf", wav=str(wav), git_sha="q", observer="o", recipe="r.sh"), root))
        self.assertFalse(res["ok"])
        self.assertIn("NOT_THE_NAMED_NOTE", res["error"])
        self.assertIn("NOT_THE_NAMED_NOTE", json.loads((root / "syn" / "Bb3" / "pair.json").read_text())["refused"]["reason"])

    def test_diff_reports_what_moved_between_two_gradings_of_a_note(self):
        from instrgrade import write_wav
        root = self.tmp / "diffroot"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="old1", observer="o", recipe="r.sh")
        dc.run_pair(self.refs, c, root)
        far = self.tmp / "far.wav"
        write_wav(far, 0.3 * tone(f0=F0 * 2 ** (1 / 12), dur=2.4, lead=0.2, noise_db=-85, cents=60))
        dc.run_pair(self.refs, dict(c, wav=str(far), git_sha="new2"), root)
        d = root / "syn" / "Bb3"
        old = __import__("json").loads(next((d / "superseded").glob("*/pair.json")).read_text())
        new = __import__("json").loads((d / "pair.json").read_text())
        rows, tot = dc.diff_pairs(old, new)
        f0 = next(r for r in rows if r["id"] == "f0_cents")
        self.assertAlmostEqual(f0["gap_old"], 8.0, delta=1.5)
        self.assertAlmostEqual(f0["gap_new"], 60.0, delta=2.0)
        self.assertGreater(f0["gap_new"], f0["gap_old"])
        self.assertNotEqual(old["ours"]["git_sha"], new["ours"]["git_sha"])
        import argparse, io, contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            dc.cmd_diff(argparse.Namespace(out=str(root), family=None, against="old1"))
        self.assertIn("old1 -> new2", buf.getvalue())
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            dc.cmd_diff(argparse.Namespace(out=str(root), family=None, against="nomatch"))
        self.assertEqual(buf.getvalue(), "")
        self.assertEqual(sum(tot["old"]), sum(tot["new"]))  # the same measures are judged both times

    def test_a_band_limited_render_is_graded_only_inside_its_band(self):
        import json
        from instrgrade import write_wav
        root = self.tmp / "band"
        f0 = F0 * 2 ** (1 / 12)
        full = 0.3 * tone(f0=f0, dur=2.4, lead=0.2, noise_db=-85)
        limited = self.tmp / "limited.wav"
        write_wav(limited, dc.lowpass(full, 3500.0))  # a render that models nothing above 3.5 kHz
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(limited), git_sha="b", observer="o", recipe="r.sh")
        wide = dc.run_pair(self.refs, c, root / "wide")
        narrow = dc.run_pair(self.refs, dict(c, band_hz=3900.0), root / "narrow")
        gap = lambda pj, m: next(x["gap"] for x in pj["metrics"] if x["id"] == m)  # noqa: E731
        self.assertLess(gap(narrow, "stft_lsd_db"), 0.6 * gap(wide, "stft_lsd_db"))
        self.assertEqual(narrow["figures"]["spectrogram"]["f_hi_hz"], 3900.0)
        self.assertEqual(narrow["figures"]["scalogram"]["f_hi_hz"], 3900.0)
        self.assertEqual(narrow["alignment"]["band_hz"], 3900.0)
        self.assertIn("low-passed at 3900 Hz", narrow["audio"]["level_match"])
        self.assertIsNone(wide["alignment"]["band_hz"])

    def test_a_failed_regrade_leaves_the_notes_earlier_pair_in_place(self):
        import json
        from unittest import mock
        root = self.tmp / "atomic"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="first", observer="o", recipe="r.sh")
        dc.run_pair(self.refs, c, root)
        d = root / "syn" / "Bb3"
        before = (d / "pair.json").read_bytes()
        with mock.patch.object(self.refs, "ruler", side_effect=RuntimeError("boom")):
            with self.assertRaises(RuntimeError):
                dc.run_pair(self.refs, dict(c, git_sha="second"), root)
        self.assertEqual((d / "pair.json").read_bytes(), before)  # the old grading is still there, untouched, mid-failure
        aborted = list((d / "superseded").glob("aborted-*"))  # the failed attempt's partial files are kept aside, not deleted
        self.assertEqual(len(aborted), 1)
        self.assertTrue(any(aborted[0].iterdir()))
        dc.run_pair(self.refs, dict(c, git_sha="third"), root)
        self.assertEqual(json.loads((d / "pair.json").read_text())["ours"]["git_sha"], "third")
        self.assertFalse([p for p in d.iterdir() if p.name.startswith(".staging")])  # nothing left staged

    def test_an_io_error_in_a_reference_run_is_not_written_into_the_data_as_a_refusal(self):
        from unittest import mock
        root = self.tmp / "ioerr"
        dc._REFS = self.refs
        ref_e = self.refs.find("syn", 58, "mf")[0]
        y, info = self.refs.aligned(ref_e)
        dc.write_pair(root, "syn", "Bb3", ref_e, self.refs, y, info)
        before = (root / "syn" / "Bb3" / "pair.json").read_bytes()
        with mock.patch.object(dc, "run_pair", side_effect=OSError("disk went away")):
            res = dc._ref_job((dict(family="syn", art="", dyn="mf", notes=["Bb3"], review=False), root))
        self.assertFalse(res["ok"])
        self.assertIn("I/O error", res["refused"][0]["reason"])
        self.assertEqual((root / "syn" / "Bb3" / "pair.json").read_bytes(), before)  # the earlier pair is untouched

    def test_pair_json_carries_both_verdicts_and_the_candidate_record(self):
        import json
        root = self.tmp / "strict"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="s", observer="o", recipe="r.sh", band_hz=None)
        pj = dc.run_pair(self.refs, c, root)
        self.validate(root, "syn", "Bb3")  # P8's validator still accepts it
        for m in pj["metrics"]:
            self.assertIn("ruler_strict", m)
            self.assertIn(m["inside_strict"], (True, False, None))
        v = pj["verdicts"]
        self.assertEqual(v["measures"], sum(1 for m in pj["metrics"] if m["gap"] is not None))
        for kind in ("loose", "strict"):
            self.assertEqual(sum(v[kind].values()), v["measures"])
        self.assertLessEqual(v["strict"]["inside"] + v["strict"]["outside"], v["loose"]["inside"] + v["loose"]["outside"])  # strict judges no more measures
        cand = pj["decompose"]["candidate"]
        self.assertEqual((cand["note"], cand["dyn"], cand["art"]), ("Bb3", "mf", ""))
        self.assertEqual(cand["ref_contains"], "Syn.mf.A3B3")  # enough to rebuild the same candidate

    def test_a_slot_gives_one_note_several_pairs(self):
        import json
        root = self.tmp / "slots"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="s", observer="o", recipe="r.sh")
        a = dc.run_pair(self.refs, dict(c, slot="Bb3-a"), root)
        b = dc.run_pair(self.refs, dict(c, slot="Bb3-e"), root)
        self.assertEqual((a["note"], b["note"]), ("Bb3-a", "Bb3-e"))
        self.assertEqual(json.loads((root / "syn" / "Bb3-e" / "pair.json").read_text())["note"], "Bb3-e")
        self.assertTrue((root / "syn" / "Bb3-a" / "ours.wav").exists())

    def test_an_io_error_in_a_render_job_is_not_recorded_as_a_refusal_and_leaves_the_pair(self):
        import json
        from unittest import mock
        root = self.tmp / "ioerr2"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="ok", observer="o", recipe="r.sh")
        dc._REFS = self.refs
        self.assertTrue(dc._batch_job((c, root))["ok"])
        before = (root / "syn" / "Bb3" / "pair.json").read_bytes()
        with mock.patch.object(dc, "_write_pair_staged", side_effect=OSError("staging directory vanished")):
            res = dc._batch_job((dict(c, git_sha="bad"), root))
        self.assertFalse(res["ok"])
        self.assertIn("not a refusal", res["error"])
        d = root / "syn" / "Bb3"
        self.assertEqual((d / "pair.json").read_bytes(), before)
        self.assertIsNone(json.loads((d / "pair.json").read_text()).get("refused"))

    def test_a_manifest_of_another_recording_set_must_carry_its_attribution(self):
        import json
        man = self.tmp / "vocalset.json"
        entry = dict(path="VocalSet/x.wav", page="p", chain="c", id="vs-1", licence="vocalset-cc-by-4.0")
        man.write_text(json.dumps({"files": [entry]}))
        with self.assertRaisesRegex(ValueError, "SOURCE_MISSING"):
            dc.Refs(self.tmp / "lib.json", self.tmp / "wav", man)
        man.write_text(json.dumps({"files": [entry], "external": {"name": "VocalSet", "url": "u", "licence": "CC BY 4.0"}}))
        with self.assertRaisesRegex(ValueError, "credit"):  # CC BY needs the credit line
            dc.Refs(self.tmp / "lib.json", self.tmp / "wav", man)
        man.write_text(json.dumps({"files": [entry], "external": {"name": "VocalSet", "url": "u", "licence": "CC BY 4.0",
                                                                  "licence_url": "lu", "credit": "Wilkins et al."}}))
        self.assertEqual(dc.Refs(self.tmp / "lib.json", self.tmp / "wav", man).external["credit"], "Wilkins et al.")
        man.write_text(json.dumps({"files": [dict(entry, licence="uiowa-mis")]}))  # the UIowa set needs none
        self.assertIsNone(dc.Refs(self.tmp / "lib.json", self.tmp / "wav", man).external)

    def test_refs_can_be_asked_for_exact_notes(self):
        import argparse, io, contextlib, json
        root = self.tmp / "exact"
        args = argparse.Namespace(lib=str(self.tmp / "lib.json"), wav_dir=str(self.tmp / "wav"), refs_json=str(self.tmp / "refs.json"), out=str(root),
                                  family=["syn"], all=False, notes="Bb3", jobs=1, review=False, since=None)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), self.assertRaises(SystemExit) as cm:
            dc.cmd_refs(args)
        self.assertEqual(cm.exception.code, 0)
        rows = [json.loads(l) for l in buf.getvalue().splitlines()]
        self.assertEqual([r["note"] for r in rows], ["Bb3"])  # exactly the note asked for
        self.assertTrue((root / "syn" / "Bb3" / "pair.json").exists())

    def test_a_refused_render_batch_job_writes_the_record(self):
        import json
        root = self.tmp / "refbatch"
        dc._REFS = self.refs
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.tmp / "nonexistent.wav"), git_sha="zz", observer="o", recipe="r.sh")
        res = dc._batch_job((c, root))
        self.assertFalse(res["ok"])
        self.assertEqual(json.loads((root / "syn" / "Bb3" / "pair.json").read_text())["refused"]["render_sha"], "zz")

    def test_json_has_no_nan(self):
        import json
        s = json.dumps(dc.clean({"a": float("nan"), "b": np.float32(1.5), "c": [np.inf]}), allow_nan=False)
        self.assertEqual(json.loads(s), {"a": None, "b": 1.5, "c": [None]})


class SubframeOnsetTests(unittest.TestCase):
    """The pair check's onset gap is measured between the 1 ms analysis frames (bead frankensim-grader-onset-subframe-zlv3e).

    The envelope frames sit on a 1 ms grid, so a crossing taken at a frame centre moves in 1 ms steps with the phase of the signal against that grid.
    The sub-frame crossing slides the same window one sample at a time over the whole envelope. Every expected value here is fixed by how the signal is
    built: an edge moved by a number of samples, a delay planted on an aligned cut, a swell too steep for align()'s iteration to settle. decompose-2 places
    the cut at the sample-resolution onset, so the onset inside a cut is the pre-roll to within a sample or two at every phase; the planted-delay test
    states what a delay does to the pair check, and the swell test reaches the per-side refusal."""

    SHIFTS = range(48)  # one analysis hop
    ONE_SAMPLE_MS = 1000.0 / SR

    @staticmethod
    def edge(shift=0, attack=0.002):
        """A tone with a 2 ms raised-cosine attack after 0.3 s of silence, the WHOLE waveform moved later by `shift` samples (zeros in front). Moving only the
        edge would change the carrier phase at the edge, which moves an energy threshold crossing by itself."""
        return np.concatenate([np.zeros(shift), tone(attack=attack, lead=0.3, dur=1.2)])

    def test_onset_times_is_the_grid_part_of_the_crossings(self):
        y = dc.highpass(self.edge(), F0)
        ta, tb, sa, sb = dc._onset_crossings(y, -100.0, None, F0)
        self.assertEqual((ta, tb), dc.onset_times(y, -100.0, None, F0))
        self.assertEqual((sa, sb), dc.onset_times_subframe(y, -100.0, None, F0))
        self.assertLessEqual(sb, tb)  # between the last frame under the threshold and the first over it
        self.assertGreater(sb, tb - 0.001 - 1e-9)

    def test_subframe_crossing_follows_the_edge_sample_for_sample_and_the_grid_one_does_not(self):
        sub, grid = [], []
        for sh in self.SHIFTS:
            y = dc.highpass(self.edge(sh), F0)
            ta, tb, sa, sb = dc._onset_crossings(y, -100.0, None, F0)
            sub.append(sb - sh / SR)
            grid.append(tb - sh / SR)
        # moving the edge by sh samples moves the sub-frame crossing by sh samples (to within one sample) ...
        self.assertLessEqual((max(sub) - min(sub)) * 1000, self.ONE_SAMPLE_MS + 1e-9)
        # ... while the frame-centre crossing is off by up to a whole millisecond from the edge over one hop
        self.assertGreater((max(grid) - min(grid)) * 1000, 0.5)

    def test_align_keeps_the_grid_onset_and_adds_the_subframe_one(self):
        for sh in (0, 11, 23, 47):
            x = self.edge(sh)
            floor = dc.render_floor(dc.highpass(x, F0), F0)
            y, info = dc.align(x, F0, floor)
            ra, rb = dc.onset_times(y, floor, None, F0)
            self.assertEqual(info["cut_onset_s"], max(ra, rb))  # unchanged: re-detection on the cut, on the grid
            self.assertLessEqual(info["cut_onset_subframe_s"], info["cut_onset_s"])
            self.assertGreater(info["cut_onset_subframe_s"], info["cut_onset_s"] - 0.001 - 1e-9)

    def test_gap_is_the_subframe_gap_and_the_grid_gap_is_kept_for_the_record(self):
        ours = dict(cut_onset_s=0.0103, cut_onset_subframe_s=0.0091)
        real = dict(cut_onset_s=0.0113, cut_onset_subframe_s=0.0104)
        gap, grid = dc.onset_gap_ms(ours, real)
        self.assertAlmostEqual(gap, 1.3, places=9)
        self.assertAlmostEqual(grid, 1.0, places=9)
        # an info from an earlier grader has no sub-frame onset: that side falls back to its grid onset
        gap_old, _ = dc.onset_gap_ms({"cut_onset_s": 0.0103}, real)
        self.assertAlmostEqual(gap_old, abs(0.0103 - 0.0104) * 1000, places=9)

    def test_the_onset_inside_the_cut_is_the_preroll_at_every_phase(self):
        # decompose-2 places the cut at the sample-resolution onset: the onset inside the cut is `preroll_s` to within a sample or two at every phase of
        # the signal against the 1 ms grid. Placed on the grid (decompose-1) it was anywhere in the hop before it.
        pre_ms = dc.preroll_s(F0) * 1000
        subs = []
        for sh in self.SHIFTS:
            x = self.edge(sh)
            _, info = dc.align(x, F0, dc.render_floor(dc.highpass(x, F0), F0))
            subs.append(info["cut_onset_subframe_s"] * 1000)
        self.assertLessEqual(max(abs(v - pre_ms) for v in subs), 2 * self.ONE_SAMPLE_MS + 1e-9, f"{min(subs):.3f}..{max(subs):.3f} ms")

    def test_step_edge_crossing_is_recovered_within_a_sample_of_the_power_linear_window(self):
        # A square wave has a constant instantaneous power, so the RMS window over a step edge is EXACTLY a linear ramp of the window start: with the edge at
        # e0 (samples, fractional: the first sample carries the fraction a of the amplitude, so e0 = n0 + 1 - a^2) the window power is A^2 (s + W - e0) / W,
        # and criterion B (steady - 20 dB, q = 1 %) is first over its threshold at the first window start s > e0 - W (1 - q). The centre of that window is
        # the crossing. The bound is one sample: the start is an integer.
        W = dc.env_window(F0)
        q = 10.0 ** (-dc.STEADY_DROP_DB / 10.0)
        A, n0 = 0.1, 14400
        for a in (0.05, 0.3, 0.5, 0.7, 0.95):
            for sh in (0, 7, 19, 31, 47):
                n = int(1.0 * SR)
                y = np.zeros(n)
                half = 109  # square wave of period 218 samples
                sq = A * np.where((np.arange(n - (n0 + sh)) // half) % 2 == 0, 1.0, -1.0)
                y[n0 + sh:] = sq
                y[n0 + sh] *= a
                e0 = n0 + sh + 1 - a * a
                s_min = math.floor(e0 - W * (1 - q)) + 1
                want = (s_min + W / 2) / SR
                _, _, _, sb = dc._onset_crossings(y, -120.0, None, F0)
                self.assertAlmostEqual(sb, want, delta=1.0 / SR + 1e-9, msg=f"a={a} shift={sh}")
                self.assertAlmostEqual(sb, (e0 - W / 2 + q * W) / SR, delta=1.0 / SR + 1e-9)  # and the closed form

    def test_planted_delay_on_an_aligned_cut_is_judged_the_same_at_every_phase(self):
        # align() cuts every signal at its own onset, so a delay of the RAW render is removed by design and cannot refuse. The phase of the raw render
        # against the 1 ms grid is set first (0..47 samples), the render is aligned, and the delay is planted on the aligned cut.
        floor = dc.render_floor(dc.highpass(self.edge(), F0), F0)
        real_cut, real_info = dc.align(self.edge(), F0, floor)

        def gaps(delay_ms):
            out = []
            for sh in self.SHIFTS:
                y, _ = dc.align(self.edge(sh), F0, floor)
                y = np.concatenate([np.zeros(int(round(delay_ms * SR / 1000))), y])
                ra, rb, sa, sb = dc._onset_crossings(y, floor, None, F0)
                out.append(dc.onset_gap_ms(dict(cut_onset_s=max(ra, rb), cut_onset_subframe_s=max(sa, sb)), real_info)[0])
            return out

        for delay, refused in ((0.0, False), (1.0, False), (1.7, True), (2.5, True)):
            g = gaps(delay)
            self.assertEqual({v > dc.ALIGN_TOL_MS for v in g}, {refused}, f"{delay} ms delay: gaps {min(g):.3f}..{max(g):.3f}")
            self.assertTrue(all(abs(v - delay) <= 3 * self.ONE_SAMPLE_MS for v in g), f"{delay} ms delay: gaps {min(g):.3f}..{max(g):.3f}")


    @staticmethod
    def swell(slope_db_s, dur=6.0, t0=0.3):
        """A tone whose level rises linearly in dB from t0 to the end of the signal (an exponential swell), 220 Hz with 11 harmonics."""
        n = int(dur * SR)
        t = np.arange(n) / SR
        x = np.zeros(n)
        on = t >= t0
        env = 10.0 ** (slope_db_s * (t[on] - t0) / 20.0)
        x[on] = env * sum(np.sin(2 * np.pi * k * F0 * t[on] + 0.3 * k) / k for k in range(1, 12))
        return 0.2 * x / np.abs(x).max()

    def test_a_swell_the_iteration_cannot_settle_on_is_refused_at_the_per_side_reading(self):
        # The per-side assertion: after the cut is placed at the sample-resolution onset, the onset re-detected on the cut (and, for a recording, re-detected again by
        # walking back from the anchor) must sit within ALIGN_TOL_MS of the pre-roll. Nothing is mocked: these are real notes through the real align().
        # A level rising linearly in dB for longer than the 3 s cut has no steady level of its own: the median of the sounding frames follows where the cut ends, so
        # every pass of align()'s four-pass iteration moves the cut by a fixed step, and a swell steeper than about 13 dB/s does not settle. A 14 dB/s swell lasting
        # 6 s is refused (about 137 ms off the pre-roll, with or without an anchor, so the walk-back re-detection of a recording does not rescue it either); the
        # same note at 10 dB/s converges and lands on the pre-roll to within a sample.
        pre_ms = dc.preroll_s(F0) * 1000
        for anchor in (None, 0.8):  # a render, then a recording with an anchor inside the rise (the walk-back re-detection runs)
            steep, gentle = self.swell(14.0), self.swell(10.0)
            with self.assertRaisesRegex(dc.OnsetMisaligned, r"IRPDECOMP_ONSET_MISALIGNED: cut onset -?[0-9.]+ ms, expected 10\.0 ms") as r:
                dc.align(steep, F0, dc.render_floor(dc.highpass(steep, F0), F0), anchor_s=anchor)
            off = float(str(r.exception).split("cut onset ")[1].split(" ms")[0])
            self.assertGreater(abs(off - pre_ms), dc.ALIGN_TOL_MS, f"anchor {anchor}: refused at {off} ms")
            _, info = dc.align(gentle, F0, dc.render_floor(dc.highpass(gentle, F0), F0), anchor_s=anchor)  # the control: within tolerance, so it is graded
            self.assertLessEqual(abs(info["cut_onset_subframe_s"] * 1000 - pre_ms), dc.ALIGN_TOL_MS, f"anchor {anchor}")


class BandDeadZoneTests(unittest.TestCase):
    """A render valid only below a band is graded only up to it (decompose-3; bead frankensim-grader-onset-subframe-zlv3e, found on trombone Bb2 at 3.9 kHz).

    Above the band the signals hold only a floor (the filter's numerical floor, the splatter of the cut), whose level follows the attack and the phase of the frames against it, so a
    measure that reads above the band compares two floors. On the 9 band-limited pairs of the site, with the dead band made 20 dB louder, the gaps of cepstral_env_db, stft_lsd_db,
    cwt_lsd_db, centroid_st, flux_db and attack_log2 moved by 0.3 to 4.2 (their own units), and hnr_db, partial_traj_db, spec_db and env_db by up to 0.04. A Note with a declared band is
    therefore analysed as the signal band-limited to it (band_limit), and the cepstral gap is taken over the grid points below the band. Every expected value here is fixed by how the
    signals are built: the same low-passed tone, one with a different floor planted above the band after the low-pass, or one with a 40 dB notch."""

    BAND = 3900.0

    @staticmethod
    def planted_floor(x, rel_db, lo=4200.0):
        """`x` plus white noise confined above `lo` Hz, at `rel_db` re the RMS of `x`: the spectrum below `lo` is `x`'s, the floor above it is new."""
        spec = np.fft.rfft(np.random.default_rng(3).standard_normal(len(x)))
        spec[np.fft.rfftfreq(len(x), 1 / SR) < lo] = 0.0
        z = np.fft.irfft(spec, len(x))
        return x + np.sqrt(np.mean(x * x)) * 10 ** (rel_db / 20) * z / np.sqrt(np.mean(z * z))

    @staticmethod
    def notch_db(f, f_c, depth_db=40.0, octaves=0.5):
        """The dB gain of a notch `depth_db` deep at `f_c` and `octaves` wide: a raised cosine in log2 frequency, 0 dB at f_c * 2**(+-octaves/2)."""
        u = np.log2(np.maximum(f, 1e-9) / f_c) / octaves
        return -depth_db * np.where(np.abs(u) < 0.5, 0.5 + 0.5 * np.cos(2 * np.pi * u), 0.0)

    def notched(self, x, f_c):
        return np.fft.irfft(np.fft.rfft(x) * 10 ** (self.notch_db(np.fft.rfftfreq(len(x), 1 / SR), f_c) / 20), len(x))

    def note(self, x, band):
        """x is already what the pipeline would hand over (low-passed, or not)."""
        floor = dc.render_floor(dc.highpass(x, F0), F0)
        y, info = dc.align(x, F0, floor)
        return Note(y, F0, "sustained", False, "synthetic", info, band_hz=band)

    def test_a_different_out_of_band_floor_moves_no_measure(self):
        """Every key compare() returns, floors planted above the band, banded against unbanded, at four phases of the frames against the attack.

        -47 dB re the tone's RMS is the loudest dead band measured on the 9 band-limited pairs (render, recording and ruler notes: -47 to -119 dB re the RMS below the band, median -77); -30 dB
        is 17 dB louder. The tolerance of each is what is left once the floor is outside the analysis: its leakage into the band across the finite cut, which scales with the floor. Measured on
        this tone: at most 0.004 (cepstral_env_db) at -47 and 0.22 (stft_lsd_db) at -30, in each measure's own unit (dB, cents, semitones, octaves of time); the tolerances are 0.02 and 0.5. A
        measure that reads above the band moves by much more than that: unbanded, at -47 hnr_db moves 100 dB, cepstral_env_db 38 and cwt_lsd_db 27; at -30 stft_lsd_db 54, cepstral_env_db
        46, cwt_lsd_db 34, spec_db and env_db 14, flux_db 1.9. The keys that cannot see a floor of this size even unbanded (f0, inharmonicity, partial frequencies, odd/even, harmonic levels,
        centroid, attack, partial trajectories) are held to the same tolerance, which is all a test on a clean tone can say about them; the table of the site's regrade says what they do on
        the real pairs."""
        base = dc.lowpass(tone(attack=0.04, lead=0.3, dur=2.6), self.BAND)
        saw = {"banded": set(), "control": set()}
        for phase in (0, 11, 23, 37):
            pad = np.zeros(phase)
            ref_b, ref_u = self.note(np.concatenate([pad, base]), self.BAND), self.note(np.concatenate([pad, base]), None)
            for rel_db, tol, control in ((-47.0, 0.02, ("hnr_db", "cepstral_env_db", "cwt_lsd_db")), (-30.0, 0.5, ("stft_lsd_db", "spec_db", "env_db", "flux_db"))):
                moved = np.concatenate([pad, self.planted_floor(base, rel_db)])
                banded = dc.compare(self.note(moved, self.BAND), ref_b)
                unbanded = dc.compare(self.note(moved, None), ref_u)
                for k, (gap, _) in banded.items():
                    if gap is None:
                        self.assertEqual(k, "decay_db_s", f"phase {phase}, {rel_db} dB: {k} has no value")  # a sustained note has no decay
                        continue
                    self.assertLess(gap, tol, f"phase {phase}, floor {rel_db} dB: {k} reads above the band ({gap:.4f})")
                for k in control:  # the control: without the band the same floor is seen by these measures, so the test could fail
                    self.assertGreater(unbanded[k][0], 1.0, f"phase {phase}, floor {rel_db} dB: the unbanded control must see the floor in {k}")

    def test_a_notch_inside_the_band_is_a_timbre_gap_and_the_same_notch_above_it_is_not(self):
        """Two tones low-passed at 3.9 kHz, with noise between their partials and a -47 dB floor above 4.2 kHz; one has a 40 dB notch, 0.5 octave wide.

        The notch is a raised cosine in log2 frequency (notch_db): at 2 kHz it spans 1.68 to 2.38 kHz and covers 3 of the 32 grid points of the envelope below the band (1/6 octave apart), at
        dips of 40, 12 and 8 dB; the RMS of those dips over the 32 points is 7.5 dB (computed below from the grid). The lifter (0.6/f0 = 2.7 ms for 220 Hz) resolves about 370 Hz and the notch is
        700 Hz wide, so the smoothed gap is a fraction of that: measured 0.56 of it (4.2 dB); the threshold is 0.4 of it, 3.0 dB. (Without the noise between the partials the log spectrum there sits at the
        numerical floor, the smoothing is dominated by it, and the same notch reads 2.8 dB.) The same notch at 5 kHz (4.2 to 5.9 kHz) is above the band, where
        there is nothing to compare but the floor: its gap must be under 0.5 dB, and without the band (the control) it is a gap of the same size as the one inside."""
        floor = self.planted_floor(dc.lowpass(tone(attack=0.04, lead=0.3, dur=2.6, noise_db=-60), self.BAND), -47.0)  # -60 dB FS of noise between the partials, as a recording has
        ref = self.note(floor, self.BAND)
        grid = ref.env_grid[ref.env_grid <= self.BAND]  # the points of the envelope that compare() reads for a banded pair
        self.assertEqual(len(grid), 32)
        for f_c, inside in ((2000.0, True), (5000.0, False)):
            expected = float(np.sqrt(np.mean(self.notch_db(grid, f_c) ** 2)))  # the unsmoothed RMS dip over the grid points below the band
            banded = dc.compare(self.note(self.notched(floor, f_c), self.BAND), ref)["cepstral_env_db"][0]
            if inside:
                self.assertAlmostEqual(expected, 7.5, delta=0.1)
                self.assertGreater(banded, 0.4 * expected, f"a 40 dB notch at {f_c:.0f} Hz, inside the band: {banded:.3f} dB against {expected:.2f} unsmoothed")
            else:
                self.assertEqual(expected, 0.0)
                self.assertLess(banded, 0.5, f"a 40 dB notch at {f_c:.0f} Hz, above the band: {banded:.3f} dB")
                unbanded = dc.compare(self.note(self.notched(floor, f_c), None), self.note(floor, None))["cepstral_env_db"][0]
                self.assertGreater(unbanded, 3.0, "the control: without the band the notch in the floor is a gap")

    def test_the_cepstral_gap_is_taken_over_the_grid_points_below_the_band_only(self):
        """Band-limited, the dead band of a Note is still not exactly alike on both sides (the leakage of the in-band lines into it is not), and the smoothing is of dB. Whatever the envelope
        holds above the band must not be in the gap: 20 dB added to every grid point above it changes nothing for a banded pair, and is a large gap for an unbanded one."""
        a0, b0 = tone(attack=0.02, dur=2.6, tilt=1.3), tone(attack=0.02, dur=2.6)
        a, b = self.note(dc.lowpass(a0, self.BAND), self.BAND), self.note(dc.lowpass(b0, self.BAND), self.BAND)
        clean = dc.compare(a, b)["cepstral_env_db"]
        self.assertGreater(clean[0], 0.5)  # a real in-band gap, from the tilt
        a.env_db = a.env_db + 20.0 * (a.env_grid > self.BAND)
        self.assertEqual(dc.compare(a, b)["cepstral_env_db"], clean)
        a, b = self.note(a0, None), self.note(b0, None)
        clean = dc.compare(a, b)["cepstral_env_db"]
        a.env_db = a.env_db + 20.0 * (a.env_grid > self.BAND)
        self.assertGreater(dc.compare(a, b)["cepstral_env_db"][0], clean[0] + 5.0)  # the control: with no band the same change is a gap

    def test_the_gap_of_a_banded_pair_is_the_same_at_every_phase(self):
        base = dc.lowpass(tone(attack=0.01, lead=0.3, dur=2.6), self.BAND)
        other = dc.lowpass(tone(attack=0.01, lead=0.3, dur=2.6, tilt=1.4), self.BAND)  # a different spectral tilt: a real in-band gap
        gaps = []
        for phase in range(0, 48, 6):
            a = self.note(np.concatenate([np.zeros(phase), other]), self.BAND)
            b = self.note(base, self.BAND)
            gaps.append(dc.compare(a, b)["cepstral_env_db"][0])
        self.assertGreater(min(gaps), 0.5)  # the gap is real
        self.assertLess(max(gaps) - min(gaps), 0.05, f"{min(gaps):.3f}..{max(gaps):.3f}")

    def test_a_note_with_no_band_is_graded_over_the_whole_grid_as_before(self):
        a, b = self.note(tone(attack=0.02, dur=2.6, tilt=1.3), None), self.note(tone(attack=0.02, dur=2.6), None)
        self.assertEqual(len(a.env_grid), 42)  # 100 Hz to 12 kHz, 1/6 octave apart
        off = b.level_db - a.level_db
        de = a.env_db + off - b.env_db
        self.assertEqual(dc.compare(a, b)["cepstral_env_db"], (float(np.sqrt(np.mean(de * de))), float(np.mean(de))))

    def test_band_limit_removes_what_is_above_the_band_and_nothing_else(self):
        x = tone(attack=0.02, dur=1.0, noise_db=-50)  # broadband
        for n in (len(x), len(x) - 1):  # even and odd lengths
            y = dc.band_limit(x[:n], self.BAND)
            X, Y = np.fft.rfft(x[:n]), np.fft.rfft(y)
            f = np.fft.rfftfreq(n, 1 / SR)
            self.assertLess(np.abs(Y[f > self.BAND]).max(), 1e-9)
            self.assertLess(np.abs(Y[f <= self.BAND] - X[f <= self.BAND]).max(), 1e-9 * np.abs(X).max())
            self.assertLess(np.abs(dc.band_limit(y, self.BAND) - y).max(), 1e-12)  # idempotent


class FileFloorOrderTests(unittest.TestCase):
    """The noise floor of a recording file does not depend on which of its notes is asked for first (bead frankensim-tgc2z).

    Refs.aligned() used to cache floor_of_file(x, f0) under (file, band) with the f0 of the first note asked, so a later note was aligned with another note's floor. floor_of_file()
    reads the note only through the corner of the note-aware high-pass, min(50 Hz, f0 / 2): every note of 100 Hz and up has the same floor, a lower note has its own. The file here
    holds A3 (220 Hz), A4 (440 Hz) and C2 (65.4 Hz) on a bed of noise and a 40 Hz hum: the hum passes the C2 corner (32.7 Hz) and is halved by the 50 Hz one, so the two floors differ by
    about 6 dB, and a hum is what a real recording has below its notes."""

    NOTES = (("hi", 57, 220.0), ("hi2", 69, 440.0), ("lo", 36, 65.406))

    @classmethod
    def setUpClass(cls):
        import json
        cls.tmp = Path(tempfile.mkdtemp())
        wav = cls.tmp / "wav"
        wav.mkdir()
        from instrgrade import write_wav
        parts = [np.zeros(SR)]
        for _, _, f in cls.NOTES:
            parts += [tone(f0=f, dur=2.0, lead=0.0, attack=0.04), np.zeros(SR)]
        x = np.concatenate(parts)
        t = np.arange(len(x)) / SR
        x = x + np.random.default_rng(11).standard_normal(len(x)) * 10 ** (-70 / 20) + 10 ** (-45 / 20) * np.sin(2 * np.pi * 40.0 * t)
        write_wav(wav / "Syn.mf.C2A4.wav", x)
        lib = []
        for k, (_, m, _) in enumerate(cls.NOTES):
            t0 = 1.0 + 3.0 * k
            lib.append(dict(file="Syn.mf.C2A4.wav", start=t0 - 0.005, end=t0 + 2.2, note=m, name=dc.name_of(m), kind="sustained", inst="syn", art="", dyn="mf", sul="", mallet="",
                            range="C2A4", desc={}))
        (cls.tmp / "lib.json").write_text(json.dumps(lib))
        (cls.tmp / "refs.json").write_text(json.dumps({"files": [dict(path="x/Syn.mf.C2A4.aiff", page="p.html", chain="c", id="syn-mf")]}))
        cls.wav = wav

    def refs(self):
        return dc.Refs(self.tmp / "lib.json", self.wav, self.tmp / "refs.json")

    def graded(self, order):
        """A fresh Refs (a fresh worker), the notes asked in this order: {name: (y, info)}."""
        refs = self.refs()
        by = {n: refs.find("syn", m, "mf")[0] for n, m, _ in self.NOTES}
        return {n: refs.aligned(by[n]) for n in order}

    def test_a_note_is_aligned_the_same_whichever_note_of_its_file_was_asked_first(self):
        first_hi, first_lo = self.graded(["hi", "hi2", "lo"]), self.graded(["lo", "hi", "hi2"])
        for n, _, _ in self.NOTES:
            (ya, ia), (yb, ib) = first_hi[n], first_lo[n]
            self.assertEqual(ia, ib, f"{n}: alignment info depends on the order")
            self.assertTrue(np.array_equal(ya, yb), f"{n}: the cut depends on the order")
        # the control: the floors of the two corners really differ, so the test can fail
        self.assertGreater(abs(first_hi["hi"][1]["floor_db"] - first_hi["lo"][1]["floor_db"]), 3.0)

    def test_notes_with_one_corner_share_one_floor_and_a_lower_note_has_its_own(self):
        from unittest import mock
        calls = []
        real = dc.floor_of_file

        def counting(x, f0, **kw):
            calls.append(f0)
            return real(x, f0, **kw)
        refs = self.refs()
        by = {n: refs.find("syn", m, "mf")[0] for n, m, _ in self.NOTES}
        with mock.patch.object(dc, "floor_of_file", counting):
            refs.aligned(by["hi"])
            refs.aligned(by["hi2"])
            self.assertEqual(len(calls), 1, "two notes of 100 Hz and up share one floor computation")
            refs.aligned(by["lo"])
            self.assertEqual(len(calls), 2, "a note below 100 Hz has its own")
            refs.aligned(by["lo"])
            refs.aligned(by["hi"])
            self.assertEqual(len(calls), 2, "and neither is computed twice")

    def test_the_corner_is_the_one_the_high_pass_uses(self):
        x = np.random.default_rng(5).standard_normal(SR)
        ref = dc.highpass(x, 220.0)
        for f0 in (100.0, 150.0, 440.0, 1000.0):
            self.assertEqual(dc.floor_corner_hz(f0), 50.0)
            self.assertTrue(np.array_equal(dc.highpass(x, f0), ref), f"highpass differs at {f0} Hz")
        for f0 in (98.0, 65.406, 41.2):
            self.assertEqual(dc.floor_corner_hz(f0), 0.5 * f0)
            self.assertFalse(np.array_equal(dc.highpass(x, f0), ref), f"highpass is the same at {f0} Hz")


def shaped_noise(dur, peak_hz=6000.0, seed=1, rms=0.05):
    """Stationary noise with one broad spectral peak (a fricative's shape, roughly): white noise weighted by a Gaussian in log2 f."""
    n = int(round(dur * SR))
    X = np.fft.rfft(np.random.default_rng(seed).standard_normal(n))
    f = np.fft.rfftfreq(n, 1 / SR)
    X *= np.exp(-0.5 * (np.log2(np.maximum(f, 1.0) / peak_hz) / 0.5) ** 2) + 0.01
    y = np.fft.irfft(X, n)
    return y * rms / np.sqrt(np.mean(y * y))


UNVOICED_COMPUTED = ("env_db", "stft_lsd_db", "cwt_lsd_db", "cepstral_env_db", "centroid_st", "flux_db")


class UnvoicedTests(unittest.TestCase):
    """pair --unvoiced: a recording with no pitch (a fricative, a flow rig's noise), chosen by its file, graded only on the measures
    that need no pitch, its pitch measures written not_applicable and counted nowhere; and the voiced path beside it unchanged."""

    @classmethod
    def setUpClass(cls):
        import json
        from instrgrade import write_wav
        cls.tmp = Path(tempfile.mkdtemp())
        wav = cls.tmp / "wav"
        wav.mkdir()
        write_wav(wav / "noise-a.wav", shaped_noise(4.0, seed=1))
        write_wav(wav / "noise-b.wav", shaped_noise(4.0, seed=2))  # another take of the same sound: the ruler
        cls.seg = shaped_noise(2.0, seed=4)
        write_wav(wav / "noise-seg.wav", cls.seg)  # a file that is exactly its segment, for the planted level
        write_wav(wav / "Syn.mf.A3.wav", np.concatenate([np.zeros(SR), tone(f0=F0, dur=2.0, lead=0.0, attack=0.04), np.zeros(SR)]))
        entry = dict(name="s", kind="sustained", inst="rig", art="", dyn="mf", sul="", mallet="", range="", desc={}, note=None, unvoiced=True, full_scale_pa=12.62)
        lib = [dict(entry, file="noise-a.wav", start=1.0, end=3.0), dict(entry, file="noise-b.wav", start=1.0, end=3.0),
               dict(entry, file="noise-seg.wav", name="z", start=0.0, end=2.0),
               dict(file="Syn.mf.A3.wav", start=0.995, end=3.2, note=57, name="A3", kind="sustained", inst="syn", art="", dyn="mf", sul="", mallet="", range="A3", desc={})]
        (cls.tmp / "lib.json").write_text(json.dumps(lib))
        (cls.tmp / "refs.json").write_text(json.dumps({
            "external": dict(name="A synthetic rig", url="https://example.org/rig", licence="CC0 1.0", licence_url="https://creativecommons.org/publicdomain/zero/1.0/",
                             credit="synthetic noise for the tests"),
            "files": [dict(path=f"x/{s}.wav", page="p.html", chain="c", id=s, licence="test-cc0") for s in ("noise-a", "noise-b", "noise-seg")]
                     + [dict(path="x/Syn.mf.A3.aiff", page="p.html", chain="c", id="syn-mf", licence="uiowa-mis")]}))
        cls.refs = dc.Refs(cls.tmp / "lib.json", wav, cls.tmp / "refs.json")
        cls.render = cls.tmp / "render.wav"
        write_wav(cls.render, shaped_noise(3.0, seed=3, rms=0.02))
        cls.shifted = cls.tmp / "shifted.wav"
        write_wav(cls.shifted, shaped_noise(3.0, peak_hz=2500.0, seed=3, rms=0.02))

    def c(self, **kw):
        c = dict(family="rig", note="s", slot="s-test", ref_file="noise-a.wav", wav=str(self.render), git_sha="abc", observer="rig mic", recipe="r.sh", full_scale_pa=1.0)
        c.update(kw)
        return c

    def validate(self, root, family, note):
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
        import build_instruments as bi
        import json
        base = Path(root) / family / note
        pj = json.loads((base / "pair.json").read_text())
        self.assertEqual(bi.validate_pair(pj, family, note, base, REFERENCE_ENTRIES) + bi.check_audio(base, pj), [])
        return pj

    def test_planted_noise_pitch_measures_are_not_applicable_and_the_others_finite(self):
        root = self.tmp / "u1"
        dc.run_unvoiced(self.refs, self.c(), root)
        pj = self.validate(root, "rig", "s-test")
        by = {m["id"]: m for m in pj["metrics"]}
        for mid, why in dc.UNVOICED_NA.items():
            m = by[mid]
            self.assertEqual(m["not_applicable"], why)
            self.assertEqual((m["gap"], m["real"], m["ours"], m["ruler"], m["ruler_strict"], m["inside_loose"], m["inside_strict"]), (None, None, None, [], [], None, None))
        self.assertEqual({"spec_db", "harm_db", "partial_traj_db", "partial_freq_cents", "f0_cents", "inharmonicity_B", "odd_even_db", "hnr_db"},
                         {k for k, v in dc.UNVOICED_NA.items() if v == "unvoiced"})
        for mid in UNVOICED_COMPUTED + ("level_db",):
            m = by[mid]
            self.assertNotIn("not_applicable", m)
            self.assertTrue(math.isfinite(m["gap"]) and m["gap"] >= 0, (mid, m["gap"]))
            self.assertEqual(len(m["ruler"]), 1, mid)  # noise-b, the other take; noise-seg is another sound ("z")
            self.assertIn("another recording", m["ruler"][0]["label"])
            self.assertEqual(m["ruler"], m["ruler_strict"])  # no neighbouring note: the loose ruler is the strict one
        self.assertNotIn("decay_db_s", by)
        self.assertEqual(set(by), set(dc.UNVOICED_NA) | set(UNVOICED_COMPUTED) | {"level_db"})
        self.assertEqual(pj["verdicts"]["measures"], 7)
        self.assertEqual(sum(pj["verdicts"]["loose"].values()), 7)
        self.assertIsNone(pj["partials"])
        self.assertEqual(pj["decompose"]["mode"], "unvoiced")
        self.assertEqual(pj["dynamic"], "mf")
        self.assertEqual(pj["real"]["segment_s"][0], 1.0)
        self.assertEqual(pj["alignment"]["segment_ours_s"][0], 1.0)  # the 3 s render's last 2 s
        self.assertEqual(pj["alignment"]["highpass_hz"], 50.0)
        self.assertEqual(pj["real"]["source"]["id"], "noise-a")
        self.assertEqual(pj["decompose"]["candidate"]["ref_file"], "noise-a.wav")

    def test_a_different_shape_is_further_than_another_take_of_the_same_one(self):
        root = self.tmp / "u2"
        pj = dc.run_unvoiced(self.refs, self.c(wav=str(self.shifted)), root)
        env = next(m for m in pj["metrics"] if m["id"] == "env_db")
        cen = next(m for m in pj["metrics"] if m["id"] == "centroid_st")
        self.assertGreater(env["gap"], 3 * env["ruler"][0]["value"])
        self.assertFalse(env["inside_loose"])
        self.assertLess(cen["signed_gap"], -3.0)  # the peak moved down by more than an octave: the centroid falls

    def test_the_level_row_reads_the_planted_pressure_and_nothing_else_moves(self):
        from instrgrade import write_wav
        root = self.tmp / "u3"
        wav = self.tmp / "half.wav"
        write_wav(wav, self.seg * (12.62 / 2.0) * 0.5)  # the recording's own samples at half its pressure, on a 2 Pa full scale
        pj = dc.run_unvoiced(self.refs, self.c(ref_file="noise-seg", wav=str(wav), full_scale_pa=2.0, slot="z-half"), root)
        by = {m["id"]: m for m in pj["metrics"]}
        self.assertAlmostEqual(by["level_db"]["signed_gap"], -20 * math.log10(2.0), delta=1e-3)
        self.assertAlmostEqual(by["level_db"]["real"] - by["level_db"]["ours"], 20 * math.log10(2.0), delta=2e-3)
        for mid in UNVOICED_COMPUTED:
            self.assertLess(by[mid]["gap"], 0.05, mid)
        self.assertEqual(by["level_db"]["ruler"], [])  # noise-seg's sound "z" has no other take
        self.assertIsNone(by["level_db"]["inside_loose"])

    def test_without_a_declared_scale_there_is_no_level_row_and_the_pair_says_why(self):
        root = self.tmp / "u4"
        pj = dc.run_unvoiced(self.refs, self.c(full_scale_pa=None), root)
        self.assertNotIn("level_db", {m["id"] for m in pj["metrics"]})
        self.assertIn("No level row: the render no --full-scale-pa.", pj["decompose"]["method"])
        self.assertNotIn("full_scale_pa", pj["ours"])
        self.assertEqual(pj["verdicts"]["measures"], 6)

    def test_a_slot_name_is_a_typed_refusal_not_a_crash(self):
        import json
        with self.assertRaisesRegex(dc.GraderRefusalError, "IRPDECOMP_NO_PITCH: 's1-11-s-rig' is not a note name"):
            dc.midi_of("s1-11-s-rig")
        self.assertTrue(issubclass(dc.GraderRefusalError, dc.REFUSALS))
        root = self.tmp / "u5"
        dc._REFS = self.refs
        res = dc._batch_job((dict(family="rig", note="s1-11-s-rig", dyn="mf", wav=str(self.render), git_sha="abc", observer="o", recipe="r"), root))
        dc._REFS = None
        self.assertFalse(res["ok"])
        self.assertIn("IRPDECOMP_NO_PITCH", json.loads((root / "rig" / "s1-11-s-rig" / "pair.json").read_text())["refused"]["reason"])

    def test_a_voiced_tone_grades_through_the_old_path_beside_unvoiced_entries(self):
        from instrgrade import write_wav
        self.assertEqual([e["file"] for e in self.refs.lib], ["Syn.mf.A3.wav"])
        self.assertEqual(len(self.refs.unvoiced), 3)
        render = self.tmp / "tone.wav"
        write_wav(render, 0.3 * tone(f0=F0, dur=2.4, lead=0.2, noise_db=-85, cents=8))
        root = self.tmp / "u6"
        pj = dc.run_pair(self.refs, dict(family="syn", note="A3", dyn="mf", wav=str(render), git_sha="abc", observer="test", recipe="r.sh"), root)
        self.validate(root, "syn", "A3")
        by = {m["id"]: m for m in pj["metrics"]}
        self.assertAlmostEqual(by["f0_cents"]["gap"], 8.0, delta=1.5)
        self.assertFalse(any("not_applicable" in m for m in pj["metrics"]))
        self.assertFalse(any("noise" in w for ws in pj["ruler_with"].values() for w in ws))
        self.assertNotIn("mode", pj["decompose"])

    def test_the_library_and_the_command_refuse_what_they_cannot_grade(self):
        import json
        bad = self.tmp / "bad-lib.json"
        bad.write_text(json.dumps([dict(file="noise-a.wav", start=1.0, end=3.0, note=60, name="s", kind="sustained", inst="rig", art="", dyn="mf", unvoiced=True)]))
        refs = dc.Refs(bad, self.tmp / "wav")  # the entry is set aside with its reason; the library itself still serves
        self.assertEqual((refs.lib, refs.unvoiced, [f for f, _ in refs.unvoiced_refused]), ([], [], ["noise-a.wav"]))
        self.assertIn("IRPDECOMP_UNVOICED_NOTE", refs.unvoiced_refused[0][1])
        with self.assertRaisesRegex(dc.GraderRefusalError, "IRPDECOMP_UNVOICED_NOTE"):
            dc.run_unvoiced(refs, self.c(), self.tmp / "u7")
        with self.assertRaisesRegex(SystemExit, "IRPDECOMP_REFERENCE_MISSING"):
            dc.run_unvoiced(self.refs, self.c(ref_file="no-such.wav"), self.tmp / "u7")
        from instrgrade import write_wav
        short = self.tmp / "short.wav"
        write_wav(short, shaped_noise(1.0, seed=5))
        with self.assertRaisesRegex(dc.GraderRefusalError, "IRPDECOMP_RENDER_SHORT"):
            dc.run_unvoiced(self.refs, self.c(wav=str(short)), self.tmp / "u7")
        with self.assertRaisesRegex(dc.GraderRefusalError, "IRPDECOMP_UNVOICED_SCALE"):
            dc.run_unvoiced(self.refs, self.c(full_scale_pa=-1.0), self.tmp / "u7")

    def test_an_unusable_entry_is_set_aside_with_its_typed_reason_and_the_rest_serves(self):
        import json
        import re
        s = dict(name="s", kind="sustained", inst="rig", art="", dyn="mf", unvoiced=True)
        lib = [dict(s, file="noise-a.wav", start=1.0, end=3.0, full_scale_pa=12.62),
               dict(s, file="start-text.wav", start="1.0", end=3.0), dict(s, file="start-missing.wav", end=3.0),
               dict(s, file="start-bool.wav", start=True, end=3.0), dict(s, file="start-nan.wav", start=float("nan"), end=3.0),
               dict(s, file="end-first.wav", start=2.0, end=1.0), dict(s, file="name-number.wav", name=5, start=1.0, end=3.0),
               dict(s, start=1.0, end=3.0)]  # the last has no file: set aside, and no pair can name it
        p = self.tmp / "set-aside-lib.json"
        p.write_text(json.dumps(lib))
        refs = dc.Refs(p, self.tmp / "wav", self.tmp / "refs.json")
        self.assertEqual([e["file"] for e in refs.unvoiced], ["noise-a.wav"])
        self.assertEqual([f for f, _ in refs.unvoiced_refused], [e.get("file") for e in lib[1:]])
        self.assertTrue(all(why.startswith("IRPDECOMP_UNVOICED_ENTRY: ") for _, why in refs.unvoiced_refused))
        for e in lib[1:-1]:
            with self.assertRaisesRegex(dc.GraderRefusalError, "IRPDECOMP_UNVOICED_ENTRY: .*" + re.escape(e["file"])):
                dc.run_unvoiced(refs, self.c(ref_file=e["file"]), self.tmp / "u10")
        for rf in (None, "", 5):  # a candidate that names no recording (a batch entry without ref_file)
            with self.assertRaisesRegex(dc.GraderRefusalError, "IRPDECOMP_USAGE: an unvoiced pair names its recording by ref_file"):
                dc.run_unvoiced(refs, self.c(ref_file=rf), self.tmp / "u10")
        self.assertEqual(dc.run_unvoiced(refs, self.c(), self.tmp / "u10")["verdicts"]["measures"], 7)  # the good entry still grades

    def test_a_scale_that_is_not_a_real_number_of_pascals_is_refused_on_either_side(self):
        import json
        for v in (True, False, "12.62", "abc", float("nan"), float("inf"), 0.0, -1.0, [1.0]):
            with self.assertRaisesRegex(dc.GraderRefusalError, "IRPDECOMP_UNVOICED_SCALE: the render"):
                dc.run_unvoiced(self.refs, self.c(full_scale_pa=v), self.tmp / "u11")
        e = dict(name="s", kind="sustained", inst="rig", art="", dyn="mf", unvoiced=True, file="noise-a.wav", start=1.0, end=3.0)
        p = self.tmp / "scale-lib.json"
        for v in (True, "12.62", "abc"):
            p.write_text(json.dumps([dict(e, full_scale_pa=v)]))
            with self.assertRaisesRegex(dc.GraderRefusalError, "IRPDECOMP_UNVOICED_SCALE: real"):
                dc.run_unvoiced(dc.Refs(p, self.tmp / "wav", self.tmp / "refs.json"), self.c(), self.tmp / "u11")

    def test_the_commands_turn_a_bad_library_into_typed_refusals_with_no_traceback(self):
        import json
        import subprocess
        from instrgrade import write_wav
        root = self.tmp / "u12"

        def run(cmd, lib, *a):
            return subprocess.run([sys.executable, "-I", str(Path(dc.__file__)), cmd, "--lib", str(lib), "--wav-dir", str(self.tmp / "wav"),
                                   "--refs-json", str(self.tmp / "refs.json"), "--out", str(root), *a], capture_output=True, text=True)

        def typed(p, code):
            self.assertEqual(p.returncode, 1, p.stderr)
            self.assertNotIn("Traceback", p.stderr)
            self.assertIn(code, p.stderr)

        lib = json.loads((self.tmp / "lib.json").read_text())
        planted = self.tmp / "planted-lib.json"  # noise-a now carries a MIDI note; noise-b, noise-seg and the voiced A3 are as before
        planted.write_text(json.dumps([dict(lib[0], note=60)] + lib[1:]))
        pair = ("--render", str(self.render), "--family", "rig", "--note", "s", "--git-sha", "abc", "--observer", "o", "--recipe", "r")
        p = run("pair", planted, "--unvoiced", "--ref-file", "noise-a.wav", "--slot", "s-planted", *pair)
        typed(p, "IRPDECOMP_UNVOICED_NOTE")
        self.assertIn("set aside: IRPDECOMP_UNVOICED_NOTE", p.stderr)
        self.assertIn("IRPDECOMP_UNVOICED_NOTE", json.loads((root / "rig" / "s-planted" / "pair.json").read_text())["refused"]["reason"])
        p = run("pair", planted, "--unvoiced", "--ref-file", "noise-b.wav", "--slot", "s-other", "--full-scale-pa", "1.0", *pair)
        self.assertEqual(p.returncode, 0, p.stderr)  # the library's other entries still serve, the voiced one too
        tone_wav = self.tmp / "tone-cli.wav"
        write_wav(tone_wav, 0.3 * tone(f0=F0, dur=2.4, lead=0.2, noise_db=-85, cents=8))
        p = run("pair", planted, "--render", str(tone_wav), "--family", "syn", "--note", "A3", "--git-sha", "abc", "--observer", "o", "--recipe", "r")
        self.assertEqual(p.returncode, 0, p.stderr)
        w = dict(family="rig", note="s", unvoiced=True, wav=str(self.render), git_sha="abc", observer="o", recipe="r", full_scale_pa=1.0)
        cands = self.tmp / "cands.json"
        cands.write_text(json.dumps([dict(w, slot="b-planted", ref_file="noise-a.wav"), dict(w, slot="b-noref"), dict(w, slot="b-ok", ref_file="noise-b.wav")]))
        p = run("batch", planted, str(cands))
        self.assertEqual(p.returncode, 1, p.stderr)
        self.assertNotIn("Traceback", p.stderr)
        res = [json.loads(line) for line in p.stdout.splitlines()]
        self.assertEqual([r["ok"] for r in res], [False, False, True])
        self.assertIn("IRPDECOMP_UNVOICED_NOTE", res[0]["error"])
        self.assertIn("IRPDECOMP_USAGE", res[1]["error"])
        for cmd, args in (("pair", ("--unvoiced", "--ref-file", "noise-b.wav", "--slot", "s-nolib", *pair)), ("batch", (str(cands),)), ("refs", ())):
            typed(run(cmd, self.tmp / "no-such-lib.json", *args), "IRPDECOMP_LIBRARY_UNREADABLE")
        self.assertFalse((root / "rig" / "s-nolib").exists())  # the run's own trouble: nothing is recorded
        (self.tmp / "shape-lib.json").write_text(json.dumps({"entries": []}))
        (self.tmp / "text-lib.json").write_text("not json")
        for name in ("shape-lib.json", "text-lib.json"):
            for cmd, args in (("pair", ("--unvoiced", "--ref-file", "noise-b.wav", "--slot", "s-shape", *pair)), ("batch", (str(cands),)), ("refs", ())):
                p = run(cmd, self.tmp / name, *args)
                typed(p, "IRPDECOMP_LIBRARY_SHAPE")
                self.assertIn("nothing was graded or recorded", p.stderr)
        self.assertFalse((root / "rig" / "s-shape").exists())  # a library fault is not the note's: nothing is recorded

    def test_a_library_fault_leaves_the_live_pair_untouched(self):
        """MB 18275: a good pair, rerun against a library that cannot be used, stops typed with no Traceback; its pair.json stays
        byte-identical and no superseded/ is written. Voiced and unvoiced alike, for text that is not JSON, JSON that is not a
        list, and a library file that is missing."""
        import json
        import subprocess
        from instrgrade import write_wav
        root = self.tmp / "u14"

        def run(lib, *a):
            return subprocess.run([sys.executable, "-I", str(Path(dc.__file__)), "pair", "--lib", str(lib), "--wav-dir", str(self.tmp / "wav"),
                                   "--refs-json", str(self.tmp / "refs.json"), "--out", str(root), "--git-sha", "abc", "--observer", "o",
                                   "--recipe", "r", *a], capture_output=True, text=True)

        tone_wav = self.tmp / "tone-fault.wav"
        write_wav(tone_wav, 0.3 * tone(f0=F0, dur=2.4, lead=0.2, noise_db=-85, cents=8))
        pairs = {("rig", "s-fault"): ("--unvoiced", "--ref-file", "noise-b.wav", "--slot", "s-fault", "--full-scale-pa", "1.0", "--render",
                                      str(self.render), "--family", "rig", "--note", "s"),
                 ("syn", "A3"): ("--render", str(tone_wav), "--family", "syn", "--note", "A3")}
        bad = {"not JSON": self.tmp / "fault-text-lib.json", "not a list": self.tmp / "fault-shape-lib.json", "missing": self.tmp / "no-such-fault-lib.json"}
        bad["not JSON"].write_text("not json")
        bad["not a list"].write_text(json.dumps({"entries": []}))
        for (family, note), args in pairs.items():
            p = run(self.tmp / "lib.json", *args)
            self.assertEqual(p.returncode, 0, p.stderr)
            d = root / family / note
            live = (d / "pair.json").read_bytes()
            self.assertFalse((d / "superseded").exists())
            for why, lib in bad.items():
                p = run(lib, *args)
                self.assertEqual(p.returncode, 1, (why, p.stderr))
                self.assertNotIn("Traceback", p.stderr)
                self.assertIn("IRPDECOMP_LIBRARY_UNREADABLE" if why == "missing" else "IRPDECOMP_LIBRARY_SHAPE", p.stderr)
                self.assertEqual((d / "pair.json").read_bytes(), live, f"{family}/{note}, library {why}: the live pair is untouched")
                self.assertFalse((d / "superseded").exists(), f"{family}/{note}, library {why}: nothing is superseded")

    def test_the_report_of_an_unvoiced_pair_says_onset_n_a_and_a_voiced_report_is_unchanged(self):
        import subprocess
        from instrgrade import write_wav
        root = self.tmp / "u13"
        pj = dc.run_unvoiced(self.refs, self.c(slot="s-report"), root)
        lines = dc.pair_report(pj).splitlines()
        self.assertTrue(lines[0].endswith(f'(ours set {pj["alignment"]["level_offset_db"]:+.1f} dB; onset n/a)'), lines[0])
        self.assertEqual(sum(line.endswith(" n/a") for line in lines[2:]), len(dc.UNVOICED_NA))
        self.assertEqual(len(lines) - 2, len(pj["metrics"]))
        no_level = dict(pj, alignment=dict(pj["alignment"], level_offset_db=None))
        self.assertTrue(dc.pair_report(no_level).splitlines()[0].endswith("(level n/a; onset n/a)"))
        render = self.tmp / "tone-report.wav"
        write_wav(render, 0.3 * tone(f0=F0, dur=2.4, lead=0.2, noise_db=-85, cents=8))
        vj = dc.run_pair(self.refs, dict(family="syn", note="A3", dyn="mf", wav=str(render), git_sha="abc", observer="test", recipe="r.sh"), root)
        al = vj["alignment"]
        head = dc.pair_report(vj).splitlines()[0]
        self.assertTrue(head.endswith(f'(ours set {al["level_offset_db"]:+.1f} dB; onset {al["onset_ours_s"]:.3f} s / {al["onset_real_s"]:.3f} s)'), head)
        self.assertNotIn(" n/a", dc.pair_report(vj))
        p = subprocess.run([sys.executable, "-I", str(Path(dc.__file__)), "report", "--out", str(root)], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("onset n/a)", p.stdout)

    def test_the_command_line_writes_the_pair_and_prints_the_not_applicable_rows(self):
        import json
        import subprocess
        root = self.tmp / "u8"
        common = ["--lib", str(self.tmp / "lib.json"), "--wav-dir", str(self.tmp / "wav"), "--refs-json", str(self.tmp / "refs.json"), "--out", str(root)]
        run = lambda *a: subprocess.run([sys.executable, "-I", str(Path(dc.__file__)), "pair", *common, *a], capture_output=True, text=True)  # noqa: E731
        p = run("--unvoiced", "--ref-file", "noise-a.wav", "--render", str(self.render), "--family", "rig", "--note", "s", "--slot", "s-cli",
                "--git-sha", "abc", "--observer", "o", "--recipe", "r", "--full-scale-pa", "1.0")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("not applicable: unvoiced", p.stdout)
        self.assertIn("not applicable: no onset", p.stdout)
        self.validate(root, "rig", "s-cli")
        p = run("--unvoiced", "--render", str(self.render), "--family", "rig", "--note", "s", "--git-sha", "abc", "--observer", "o", "--recipe", "r")
        self.assertIn("IRPDECOMP_USAGE", p.stderr)
        p = run("--ref-file", "noise-a.wav", "--render", str(self.render), "--family", "rig", "--note", "A3", "--git-sha", "abc", "--observer", "o", "--recipe", "r")
        self.assertIn("IRPDECOMP_USAGE: --ref-file and --full-scale-pa belong to --unvoiced", p.stderr)
        p = run("--render", str(self.render), "--family", "rig", "--note", "s1-11-s-rig", "--git-sha", "abc", "--observer", "o", "--recipe", "r")
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("IRPDECOMP_NO_PITCH", p.stderr)
        self.assertNotIn("Traceback", p.stderr)
        self.assertIn("IRPDECOMP_NO_PITCH", json.loads((root / "rig" / "s1-11-s-rig" / "pair.json").read_text())["refused"]["reason"])

    def test_the_unvoiced_high_pass_is_the_corner_of_every_note_of_100_hz_and_up(self):
        x = np.random.default_rng(6).standard_normal(SR)
        self.assertTrue(np.array_equal(dc.unvoiced_highpass(x), dc.highpass(x, 220.0)))

    def test_the_explorer_shows_not_applicable_rows_and_counts_none_of_them(self):
        import importlib.util
        path = Path(__file__).resolve().parents[3] / "instrument-explorer" / "tools" / "build.py"
        spec = importlib.util.spec_from_file_location("explorer_build", path)
        eb = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(eb)
        root = self.tmp / "u9"
        pj = dc.run_unvoiced(self.refs, self.c(), root)
        self.assertEqual(eb.pair_kind(pj), "graded")
        views = [eb.metric_view(m, True) for m in pj["metrics"]]
        na = [v for v in views if v.get("not_applicable")]
        self.assertEqual({v["id"] for v in na}, set(dc.UNVOICED_NA))
        self.assertTrue(all(v["gap"] is None and v["inside"] is None for v in na))
        self.assertEqual(eb.tally(views)["measures"], 7)
        self.assertNotIn("not_applicable", eb.metric_view(dict(id="f0_cents", label="f0", gap=1.0, ruler=[]), True))  # a voiced row is unchanged


if __name__ == "__main__":
    unittest.main()
