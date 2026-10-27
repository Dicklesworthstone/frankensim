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
        ci = dict(ci, cut_onset_s=ci["cut_onset_s"] + 0.01)  # the render's onset sits 10 ms late in its own cut
        meta = dict(wav="r.wav", git_sha="x", observer="o", sample_rate=48000, recipe="r.sh")
        with self.assertRaises(dc.OnsetMisaligned):
            dc.write_pair(self.tmp / "bad", "syn", "Bb3", ref_e, self.refs, y, info, cy, ci, meta)

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


if __name__ == "__main__":
    unittest.main()
