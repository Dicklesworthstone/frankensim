"""G0/G1/G2/G3 unit and synthetic-oracle tests for M7.1-W2 grader blind spots and Q3 descriptors.

Run with:
  uv run --with numpy --with matplotlib python -m unittest discover -s apps/instrument-lab/tools/instrgrade -p 'test_*.py' -v
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

TOOLS = Path(__file__).resolve().parent
ROOT = TOOLS.parents[3]
sys.path.insert(0, str(TOOLS))
import instrgrade as ig  # noqa: E402

try:
    import descriptors as desc_mod  # noqa: E402
except ImportError:
    desc_mod = None


def required_env(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise AssertionError(f"IRPGRADE_INPUT_MISSING: set {name}")
    return value


def required_file(path: Path) -> Path:
    path = Path(path)
    if not path.is_file():
        raise AssertionError(f"IRPGRADE_INPUT_MISSING: {path}")
    return path


class BlindSpotsAndDescriptors(unittest.TestCase):
    """M7.1-W2 §7.2 blind-spot fixes, Q3 bound descriptors, and §7.2 checks."""

    def test_attack_window_sees_attack_only_change(self):
        """G1: 40 ms onset noise burst changes attack_erb_distance_db >= 3 dB while sustain <= 0.05 dB."""
        rng = np.random.default_rng(20261016)
        t = np.arange(2 * ig.SR, dtype=np.float64) / ig.SR
        base = sum(np.sin(2.0 * np.pi * 220.0 * k * t) / k for k in range(1, 11))
        burst = np.zeros_like(t)
        n_burst = int(0.040 * ig.SR)
        burst[:n_burst] = rng.normal(0.0, 0.8, n_burst) * np.hanning(n_burst)
        x_clean = base.copy()
        x_burst = base + burst

        da = ig.describe(x_clean, kind="sustained", f0_hint=220.0, noise_floor_db=-80.0) if "noise_floor_db" in ig.describe.__code__.co_varnames else ig.describe(x_clean, kind="sustained", f0_hint=220.0)
        db = ig.describe(x_burst, kind="sustained", f0_hint=220.0, noise_floor_db=-80.0) if "noise_floor_db" in ig.describe.__code__.co_varnames else ig.describe(x_burst, kind="sustained", f0_hint=220.0)
        dist = ig.distance(da, db)

        attack_dist = dist.get("attack_erb_distance_db", dist.get("spec_db", 0.0))
        sustain_dist = dist.get("sustain_erb_distance_db", dist.get("spec_db", 0.0))
        self.assertIn("attack_erb_distance_db", dist, "W1 grader lacks attack_erb_distance_db")
        self.assertGreaterEqual(attack_dist, 3.0)
        self.assertLessEqual(sustain_dist, 0.05)

    def test_release_window_sees_release_only_change(self):
        """G1: identical attack/sustain with release tau 0.1 s vs 0.4 s gives release >= 3 dB, sustain <= 0.05 dB."""
        t = np.arange(int(2.4 * ig.SR), dtype=np.float64) / ig.SR
        carrier = sum(np.sin(2.0 * np.pi * 220.0 * k * t) / k for k in range(1, 11))
        t_rel = 1.0
        env_fast = np.where(t <= t_rel, 1.0, np.exp(-(t - t_rel) / 0.10))
        env_slow = np.where(t <= t_rel, 1.0, np.exp(-(t - t_rel) / 0.40))
        attack_ramp = np.minimum(t / 0.03, 1.0)
        x_fast = carrier * attack_ramp * env_fast
        x_slow = carrier * attack_ramp * env_slow

        if "noise_floor_db" in ig.describe.__code__.co_varnames:
            da = ig.describe(x_fast, kind="sustained", f0_hint=220.0, noise_floor_db=-80.0)
            db = ig.describe(x_slow, kind="sustained", f0_hint=220.0, noise_floor_db=-80.0)
        else:
            da = ig.describe(x_fast, kind="sustained", f0_hint=220.0)
            db = ig.describe(x_slow, kind="sustained", f0_hint=220.0)
        dist = ig.distance(da, db)

        self.assertIn("release_erb_distance_db", dist, "W1 grader lacks release_erb_distance_db")
        self.assertGreaterEqual(dist["release_erb_distance_db"], 3.0)
        self.assertLessEqual(dist["sustain_erb_distance_db"], 0.05)

    def test_two_slope_decay_recovered(self):
        """G1: level = 10*log10(0.9*e^(-t/0.25) + 0.1*e^(-t/2.5)) on 330 Hz tone, 6 s, floor -90 dB."""
        rng = np.random.default_rng(33006)
        t = np.arange(6 * ig.SR, dtype=np.float64) / ig.SR
        power_env = 0.9 * np.exp(-t / 0.25) + 0.1 * np.exp(-t / 2.5)
        x = np.sqrt(2.0 * power_env) * np.sin(2.0 * np.pi * 330.0 * t)
        x += rng.normal(0.0, 10.0 ** (-90.0 / 20.0), len(t))

        s1_true = -10.0 * math.log10(math.e) / 0.25
        s2_true = -10.0 * math.log10(math.e) / 2.5
        knee_true = math.log(9.0) / (1.0 / 0.25 - 1.0 / 2.5)

        if "noise_floor_db" in ig.describe.__code__.co_varnames:
            d = ig.describe(x, kind="decaying", f0_hint=330.0, noise_floor_db=-90.0)
        else:
            d = ig.describe(x, kind="decaying", f0_hint=330.0)

        self.assertIn("slope1_db_s", d, f"W1 decay_profile only fits single slope decay_db_s={d.get('decay_db_s')}")
        self.assertLessEqual(abs(d["slope1_db_s"] - s1_true) / abs(s1_true), 0.05)
        self.assertLessEqual(abs(d["slope2_db_s"] - s2_true) / abs(s2_true), 0.05)
        frame_hop_s = d.get("decay_frame_hop_s", 1024.0 / ig.SR)
        self.assertLessEqual(abs(d["knee_s"] - knee_true), frame_hop_s)

    def test_beat_rate_and_depth_recovered(self):
        """G1: partials at 440 and 441.3 Hz, amplitude ratio 0.5 -> rate 1.30 +/- 0.03 Hz, depth 9.54 +/- 0.5 dB."""
        t = np.arange(4 * ig.SR, dtype=np.float64) / ig.SR
        x = np.sin(2.0 * np.pi * 440.0 * t) + 0.5 * np.sin(2.0 * np.pi * 441.3 * t)
        if "noise_floor_db" in ig.describe.__code__.co_varnames:
            d = ig.describe(x, kind="decaying", f0_hint=440.0, noise_floor_db=-80.0, allow_unended_decay=True)
        else:
            d = ig.describe(x, kind="decaying", f0_hint=440.0)

        self.assertIn("beat_rate_hz", d, "W1 decay_profile does not recover beat_rate_hz or beat_depth_db")
        self.assertAlmostEqual(d["beat_rate_hz"], 1.30, delta=0.03)
        expected_depth_db = 20.0 * math.log10(1.5 / 0.5)
        self.assertAlmostEqual(d["beat_depth_db"], expected_depth_db, delta=0.5)

    def test_stiff_partials_are_not_noise(self):
        """G1: f0 = 110 Hz, B = 4e-4, 30 partials, no noise -> band HNR >= 60 dB and fitted B within 2%."""
        t = np.arange(2 * ig.SR, dtype=np.float64) / ig.SR
        f0_true = 110.0
        b_true = 4.0e-4
        x = np.zeros_like(t)
        for k in range(1, 31):
            fk = k * f0_true * math.sqrt(1.0 + b_true * k * k)
            x += np.sin(2.0 * np.pi * fk * t) / k

        d = ig.describe(x, kind="sustained", f0_hint=f0_true)
        self.assertGreaterEqual(d["hnr_db"], 60.0, f"W1 comb mask misclassified stiff partials as noise: hnr_db={d['hnr_db']:.2f}")
        self.assertIn("inharmonicity_B", d, "W1 describe does not fit stiff-string inharmonicity B")
        self.assertLessEqual(abs(d["inharmonicity_B"] - b_true) / b_true, 0.02)

    def test_32hz_fundamental_survives(self):
        """G1: 32.70 Hz (note C1) 1/n tone -> fundamental level after note-aware high-pass within 0.5 dB."""
        t = np.arange(2 * ig.SR, dtype=np.float64) / ig.SR
        f0 = 32.70
        x = sum(np.sin(2.0 * np.pi * k * f0 * t) / k for k in range(1, 16))
        y = ig.highpass(x, f0) if "f0_expected" in ig.highpass.__code__.co_varnames else ig.highpass(x)

        f_pre, s_pre = ig.spectrum(x, 0.2, 1.8)
        f_post, s_post = ig.spectrum(y, 0.2, 1.8)
        band = (f_pre >= f0 * 0.95) & (f_pre <= f0 * 1.05)
        lvl_pre_db = 10.0 * math.log10(float(s_pre[band].sum()) + 1e-30)
        lvl_post_db = 10.0 * math.log10(float(s_post[band].sum()) + 1e-30)
        self.assertLessEqual(abs(lvl_post_db - lvl_pre_db), 0.5,
                             f"W1 fixed 30-50 Hz highpass attenuated 32.70 Hz fundamental by {lvl_pre_db - lvl_post_db:.2f} dB")

    def test_pitch_synchronous_vs_stationary_noise(self):
        """G1: equal-power noise (20% duty at 196 Hz vs stationary), same band HNR within 1 dB -> cyclo >= 0.5 vs <= 0.1."""
        rng = np.random.default_rng(19620)
        t = np.arange(2 * ig.SR, dtype=np.float64) / ig.SR
        f0 = 196.0
        tone = sum(np.sin(2.0 * np.pi * k * f0 * t) / k for k in range(1, 31))
        noise = rng.normal(0.0, 0.05, len(t))
        gate = (((t * f0) % 1.0) < 0.20).astype(np.float64) * math.sqrt(5.0)
        x_sync = tone + noise * gate
        x_stat = tone + noise

        d_sync = ig.describe(x_sync, kind="sustained", f0_hint=f0)
        d_stat = ig.describe(x_stat, kind="sustained", f0_hint=f0)
        self.assertLessEqual(abs(d_sync["hnr_db"] - d_stat["hnr_db"]), 1.0)
        self.assertIn("cyclo_coherence_kf0", d_sync, "W1 describe lacks cyclo_coherence_kf0")
        self.assertGreaterEqual(d_sync["cyclo_coherence_kf0"], 0.5)
        self.assertLessEqual(d_stat["cyclo_coherence_kf0"], 0.1)

    def test_decay_analysed_to_noise_floor(self):
        """G1: 6 s decaying 82.4 Hz (E2) note with tau = 1.2 s over -80 dB floor -> end >= 5 s within 0.1 s of crossing; 3 s refuses AnalysisTooShort."""
        rng = np.random.default_rng(82406)
        f0 = 82.4
        tau = 1.2
        # Signal crosses -80 dB floor at t_cross = 5.10 s (>= 5 s) so signal+noise crosses -77 dB at t_cross
        t_cross_true = 5.10
        alpha = (80.0 / 20.0) * math.log(10.0) / t_cross_true
        t6 = np.arange(6 * ig.SR, dtype=np.float64) / ig.SR
        noise6 = rng.normal(0.0, 10.0 ** (-80.0 / 20.0), len(t6))
        x6 = math.sqrt(2.0) * np.exp(-alpha * t6) * np.sin(2.0 * np.pi * f0 * t6) + noise6
        x3 = x6[: 3 * ig.SR]

        self.assertIsNotNone(desc_mod, "W1 lacks descriptors module and AnalysisTooShort refusal")
        win6 = desc_mod.windows(x6, "decaying", noise_floor_db=-80.0)
        t_end = win6["release"][1]
        self.assertGreaterEqual(t_end, 5.0)
        self.assertLessEqual(abs(t_end - t_cross_true), 0.10)
        with self.assertRaises(desc_mod.AnalysisTooShort):
            desc_mod.windows(x3, "decaying", noise_floor_db=-80.0)
        with self.assertRaises(desc_mod.AnalysisTooShort):
            ig.describe(x3, kind="decaying", f0_hint=f0, noise_floor_db=-80.0)

    def test_pitch_shift_invariance(self):
        """G3: a note resampled +-20 cents moves spec_db and the attack, sustain and release ERB distances < 1 dB.

        The notes have what real ones have: a 30 ms attack, 0.1 dB rms slow amplitude ripple, a release at 1.6 s and a
        -80 dB noise floor. Triangular ERB bands read a partial just outside their edge through window leakage (one
        band moved 26 dB), bands far below f0 read only a partial's filter skirt, and a knee picked within 0.02 dB of
        the plateau jittered by tens of ms: at 330 and 587 Hz these distances moved 3.4-7.4 dB while spec_db stayed
        below 1 dB."""
        self.assertIsNotNone(desc_mod, "descriptors module required")
        for f0 in (147.0, 330.0, 587.0):
            rng = np.random.default_rng(20261028)
            t = np.arange(int(2.4 * ig.SR), dtype=np.float64) / ig.SR
            tone = sum(np.sin(2.0 * np.pi * k * f0 * t + rng.uniform(0.0, 2.0 * np.pi)) / (k ** 1.2) for k in range(1, 20))
            slow = np.convolve(rng.normal(size=len(t)), np.hanning(int(0.1 * ig.SR)), "same")
            env = np.minimum(t / 0.03, 1.0) * np.where(t <= 1.6, 1.0, np.exp(-(t - 1.6) / 0.15))
            x = tone * env * 10.0 ** (0.1 * slow / slow.std() / 20.0)
            x = 0.5 * x / np.abs(x).max() + rng.normal(0.0, 0.5e-4, len(t))
            da = ig.describe(x, kind="sustained", f0_hint=f0)
            for cents in (20.0, -20.0):
                ratio = 2.0 ** (cents / 1200.0)
                n_shift = int(round(len(x) / ratio))
                X = np.fft.rfft(x)
                Y = np.zeros(n_shift // 2 + 1, dtype=complex)
                m = min(len(X), len(Y))
                Y[:m] = X[:m]
                x_shifted = np.fft.irfft(Y, n_shift) * (n_shift / len(x))
                d = ig.distance(da, ig.describe(x_shifted, kind="sustained", f0_hint=f0 * ratio))
                for key in ("spec_db", "attack_erb_distance_db", "sustain_erb_distance_db", "release_erb_distance_db"):
                    self.assertLess(d[key], 1.0, f"{f0:.0f} Hz {cents:+.0f} c: {key}")

    def test_sustained_windows_follow_the_note(self):
        """G1: a sustained note released at 1.6 s keeps its sustain to the release whether it rises 1.5 dB over its
        first 0.7 s or dips 2 dB for 60 ms mid-note (the first 1 dB drop after 0.20 s ended it at 0.20 s and 0.745 s;
        80 % of UIowa sustained notes got a sustain of 0.10 s or less). A release tail that touches the floor in brief
        dips ends where it stays there, not at the first dip."""
        self.assertIsNotNone(desc_mod, "descriptors module required")
        t = np.arange(int(2.4 * ig.SR), dtype=np.float64) / ig.SR
        tone = sum(np.sin(2.0 * np.pi * 220.0 * k * t) / k for k in range(1, 11))
        gate = np.minimum(t / 0.03, 1.0) * np.where(t <= 1.6, 1.0, np.exp(-(t - 1.6) / 0.15))
        rise = 10.0 ** (np.clip((t - 0.15) / 0.55, 0.0, 1.0) * 1.5 / 20.0)
        dip = 10.0 ** (-2.0 * ((t > 0.8) & (t < 0.86)) / 20.0)
        for name, env in (("rise", rise), ("dip", dip)):
            w = desc_mod.windows(tone * gate * env, "sustained", noise_floor_db=-80.0)
            self.assertAlmostEqual(w["sustain"][1], 1.6 - 0.05, delta=0.05, msg=name)
        # A tail 3 dB over the floor threshold (-37 dB rms against -40 dB) from 1.6 s, with three 40 ms dips under it at
        # 1.8, 2.0 and 2.2 s, and at the floor from 2.6 s.
        t3 = np.arange(int(3.0 * ig.SR), dtype=np.float64) / ig.SR
        level_db = np.where(t3 <= 1.6, 0.0, np.where(t3 <= 2.6, -34.0, -60.0))
        for t_dip in (1.8, 2.0, 2.2):
            level_db[(t3 > t_dip) & (t3 < t_dip + 0.04)] = -60.0
        tail = np.sin(2.0 * np.pi * 220.0 * t3) * 10.0 ** (level_db / 20.0)
        w = desc_mod.windows(tail, "sustained", noise_floor_db=-43.0)
        self.assertAlmostEqual(w["release"][1], 2.6, delta=0.05)

    def test_identity_and_symmetry(self):
        """G0: self-distance 0; symmetric within 1e-12."""
        self.assertIsNotNone(desc_mod, "descriptors module required")
        t = np.arange(2 * ig.SR, dtype=np.float64) / ig.SR
        x = sum(np.sin(2.0 * np.pi * 220.0 * k * t) / k for k in range(1, 13))
        y = sum(np.sin(2.0 * np.pi * 220.0 * k * t) / (k * k) for k in range(1, 13))
        da = ig.describe(x, kind="sustained", f0_hint=220.0, noise_floor_db=-80.0)
        db = ig.describe(y, kind="sustained", f0_hint=220.0, noise_floor_db=-80.0)
        d_aa = ig.distance(da, da)
        self.assertEqual(d_aa["spec_db"], 0.0)
        self.assertEqual(d_aa["attack_erb_distance_db"], 0.0)
        self.assertEqual(d_aa["sustain_erb_distance_db"], 0.0)
        self.assertEqual(d_aa["release_erb_distance_db"], 0.0)
        d_ab = ig.distance(da, db)
        d_ba = ig.distance(db, da)
        for key in ("spec_db", "env_db", "harm_db", "attack_erb_distance_db", "sustain_erb_distance_db", "release_erb_distance_db"):
            self.assertAlmostEqual(d_ab[key], d_ba[key], delta=1e-12)

    def test_partial_onset_order(self):
        """G1: synthetic partial onsets staggered 0/5/10/15 ms recover order tau = 1."""
        self.assertIsNotNone(desc_mod, "descriptors module required")
        t = np.arange(int(1.0 * ig.SR), dtype=np.float64) / ig.SR
        f0 = 220.0
        staggers_s = [0.000, 0.005, 0.010, 0.015]
        x = np.zeros_like(t)
        for idx, onset_s in enumerate(staggers_s):
            k = idx + 1
            gate = np.clip((t - (0.020 + onset_s)) / 0.002, 0.0, 1.0)
            x += gate * np.sin(2.0 * np.pi * k * f0 * t) / k
        partials = desc_mod.track_partials(x, f0, kmax=4)
        self.assertEqual(len(partials), 4)
        tau = desc_mod.partial_onset_order([p.onset_s for p in partials], [0.020 + s for s in staggers_s])
        self.assertAlmostEqual(tau, 1.0, delta=1e-9)

    def test_intonation_fingerprint(self):
        """G1: partials detuned +3/-2/+5 cents recover the fingerprint within 0.3 cents."""
        self.assertIsNotNone(desc_mod, "descriptors module required")
        t = np.arange(2 * ig.SR, dtype=np.float64) / ig.SR
        f0 = 220.0
        detune_cents = [0.0, 3.0, -2.0, 5.0]
        x = np.zeros_like(t)
        for idx, c in enumerate(detune_cents):
            k = idx + 1
            fk = k * f0 * (2.0 ** (c / 1200.0))
            x += np.sin(2.0 * np.pi * fk * t) / k
        partials = desc_mod.track_partials(x[int(0.15 * ig.SR):int(1.85 * ig.SR)], f0, kmax=4, fit_inharmonicity=False)
        fp = desc_mod.intonation_fingerprint(partials, f0)
        self.assertEqual(len(fp), 4)
        for got, exp in zip(fp[1:], detune_cents[1:]):
            self.assertAlmostEqual(got, exp, delta=0.3)

    def test_other_descriptors_and_refusals(self):
        """G0/G1: log-attack time, spectral flux, AM/FM modulation, jitter, shimmer, band split, and typed refusals."""
        self.assertIsNotNone(desc_mod, "descriptors module required")
        t = np.arange(2 * ig.SR, dtype=np.float64) / ig.SR
        f0 = 220.0
        am = 1.0 + 0.25 * np.sin(2.0 * np.pi * 5.0 * t)
        phase = 2.0 * np.pi * f0 * t + (0.5 / 6.0) * np.sin(2.0 * np.pi * 6.0 * t)
        ramp = np.clip(t / 0.05, 0.0, 1.0)
        x = ramp * am * (np.sin(phase) + 0.5 * np.sin(2.0 * phase) + 0.25 * np.sin(3.0 * phase))
        d = ig.describe(x, kind="sustained", f0_hint=f0, noise_floor_db=-80.0, band_split_hz=1000.0)
        self.assertIn("q3", d)
        self.assertIn("descriptors", d)
        self.assertIn("log_attack_time_log10_s", d["descriptors"])
        self.assertIn("spectral_flux_db", d["descriptors"])
        self.assertIn("am_fm_modulation", d["descriptors"])
        self.assertIn("jitter_rel", d["descriptors"])
        self.assertIn("shimmer_db", d["descriptors"])
        self.assertIn("band_split", d["descriptors"])
        self.assertAlmostEqual(d["descriptors"]["am_fm_modulation"]["am_rate_hz"]["value"], 5.0, delta=0.2)

        with self.assertRaises(desc_mod.MissingNoiseFloorEstimate):
            desc_mod.windows(x, "sustained", noise_floor_db=None)
        with self.assertRaises(desc_mod.MissingBandSplitFrequency):
            desc_mod.band_split_statistics(x, [], None)

    def test_baseline_ordering(self):
        """G2: for every graded reference note, real-vs-real dynamics ruler < every null baseline distance."""
        self.assertIsNotNone(desc_mod, "descriptors module required")
        scratch = Path(required_env("IRP_GRADE_SCRATCH"))
        libpath = required_file(scratch / "grade" / os.environ.get("LIBFILE", "library.json"))
        library = json.loads(libpath.read_text())
        for note_name in ("D3", "E3"):
            targets = ig.find(library, "clarinet", ig.midi_of(note_name), dyn="mf")
            self.assertTrue(targets, f"missing clarinet {note_name} mf in library")
            target = targets[0]
            rl = ig.rulers(library, target)
            dyn_ruler = min(rl["other_dynamic"])
            baselines = desc_mod.null_baselines(target["desc"], dyn_ruler_db=dyn_ruler)
            for name in ("sine", "1/n", "odd_1/n", "1/n2"):
                self.assertIn(name, baselines["distances_db"])
                self.assertGreater(baselines["distances_db"][name], dyn_ruler)
                self.assertGreater(baselines["margins_vs_dynamic_ruler_db"][name], 0.0)
            ordered = [baselines["distances_db"][k] for k in baselines["ordering"]]
            self.assertEqual(ordered, sorted(ordered))

    def test_segmentation_counts(self):
        """G2: every Iowa scale file used in IRP_UIOWA_WAV48 yields its expected note count,
        each span sounding its named note; a refusal must be a measured header defect."""
        self.assertIsNotNone(desc_mod, "descriptors module required")
        wavdir = Path(required_env("IRP_UIOWA_WAV48"))
        wavs = sorted(wavdir.glob("*.wav"))
        self.assertTrue(wavs, "no WAV files found in IRP_UIOWA_WAV48")
        wrong, refused = [], []
        for wav in wavs:
            meta = ig.parse_ref_name(wav.stem)
            if not meta["range"]:
                continue
            exp = ig.expected_notes(meta["range"])
            x = ig.highpass(ig.load(wav), ig.hz_of(exp[0]))
            try:
                spans = ig.segment_notes(x, exp)
            except ig.ReferenceMismatch as refusal:
                if refusal.kind == "note-absent":
                    # The other named notes verified at pitch; these never sound.
                    refused.append((wav.name, refusal.kind, refusal.detail))
                elif refusal.kind == "no-named-note":
                    # Harmonic-comb pitch of each silence span; the header factor is
                    # the one that maps the majority of spans onto named notes.
                    pitches = []
                    for a, b in ig.segment_scale(x):
                        f, s = ig.spectrum(x[a:b], 0.0, (b - a) / ig.SR)
                        pitches.append(ig.f0_spectral(f, s, 25.0, 2000.0))
                    refused.append((wav.name, refusal.kind, pitches))
                else:
                    wrong.append(f"{wav.name}: {refusal}")
                continue
            if len(spans) != len(exp):
                wrong.append(f"{wav.name}: got {len(spans)} spans, expected {len(exp)}")
        self.assertEqual(wrong, [])
        # A file sounding none of its notes must be a sample-rate header mislabel
        # (96 or 48 kHz data declared 44.1 kHz), measured, not assumed: undoing
        # that factor puts a majority of its spans on named notes.
        mislabels = (44_100 / 96_000, 44_100 / 48_000)
        report = []
        for name, kind, detail in refused:
            if kind != "no-named-note":
                report.append(dict(file=name, kind=kind, absent=detail))
                continue
            exp = ig.expected_notes(ig.parse_ref_name(Path(name).stem)["range"])
            hits = {r: sum(any(abs(12 * math.log2(p / r / ig.hz_of(m))) <= 0.5 for m in exp) for p in detail)
                    for r in (1.0, *mislabels)}
            factor = max(mislabels, key=hits.get)
            report.append(dict(file=name, kind=kind, spans=len(detail), named_hits=hits[1.0],
                               factor=factor, named_hits_after_factor=hits[factor]))
            self.assertTrue(2 * hits[factor] > len(detail), f"{name}: {hits} over {len(detail)} spans")
        print(json.dumps(dict(test="segmentation_counts", files=len(wavs), refused=report)))


if __name__ == "__main__":
    unittest.main()
