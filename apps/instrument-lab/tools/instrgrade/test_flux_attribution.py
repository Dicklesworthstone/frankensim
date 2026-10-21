from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import decompose as dc  # noqa: E402
import flux_attribution as fa  # noqa: E402
from decompose import SR  # noqa: E402

F0 = 200.0


def voice_like(vib_cents=0.0, vib_hz=4.0, shimmer_db=0.0, shim_hz=2.5, noise_db=None, dur=2.5, seed=1, partial_am_db=0.0):
    """Fourteen harmonics at 1/k, with a sinusoidal pitch movement, a sinusoidal level movement (all partials together), a level movement
    of each partial on its own phase (`partial_am_db`, 2.5 Hz) and broadband noise (dB re the signal's RMS)."""
    n = int(dur * SR)
    t = np.arange(n) / SR
    cents = vib_cents * np.sin(2 * np.pi * vib_hz * t)
    phase = 2 * np.pi * np.cumsum(F0 * 2 ** (cents / 1200)) / SR
    x = sum(np.sin(k * phase + 0.3 * k) / k * 10 ** (partial_am_db * np.sin(2 * np.pi * 2.5 * t + 2.1 * k) / 20) for k in range(1, 15))
    x = x * 10 ** (shimmer_db * np.sin(2 * np.pi * shim_hz * t) / 20) * np.clip(t / 0.03, 0, 1)
    x = 0.2 * x / np.abs(x).max()
    if noise_db is not None:
        x = x + np.random.default_rng(seed).standard_normal(n) * np.sqrt(np.mean(x * x)) * 10 ** (noise_db / 20)
    return x


class FluxAttributionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.steady = voice_like()
        cls.vibrato = voice_like(vib_cents=60.0)
        cls.small_vibrato = voice_like(vib_cents=12.0, vib_hz=5.0)
        cls.shimmer = voice_like(shimmer_db=3.0)
        cls.noisy = voice_like(noise_db=-32.0)

    def share(self, res, key, steady=False):
        v = res[key]["steady_frames"] if steady else res[key]
        return v["gap_explained_share"]

    def test_the_f0_follows_a_vibrato_to_a_cent_and_a_steady_tone_reads_steady(self):
        look = fa.Look(self.vibrato, F0)
        t = np.arange(len(self.vibrato)) / SR
        truth = 60.0 * np.sin(2 * np.pi * 4.0 * t)
        c = 1200 * np.log2(look.f0_inst / np.median(look.f0_inst[look.sound_samples]))
        ok = (t > 0.2) & (t < 2.2)
        self.assertLess(np.sqrt(np.mean((c[ok] - truth[ok]) ** 2)), 3.0)  # cents, against a 42 cent rms movement
        s = fa.track_stats(look)
        self.assertAlmostEqual(s["vibrato_rate_hz"], 4.0, delta=0.5)
        self.assertGreater(s["vibrato_extent_cents_pp"], 80.0)  # 120 pp, less what an 85 ms window averages away
        calm = fa.track_stats(fa.Look(self.steady, F0))  # the measurement floor on a tone that does not move
        self.assertLess(calm["f0_jitter_cents_rms"], 2.0)
        self.assertLess(calm["f0_movement_cents_rms"], 4.0)
        self.assertLess(calm["level_shimmer_db_rms"], 0.5)  # the floor: the onset ramp pulls the smoothed baseline a third of a dB

    def test_two_identical_signals_have_no_gap_and_nothing_to_explain(self):
        res = fa.attribute(voice_like(), self.steady, F0)
        self.assertEqual(res["baseline"]["gap"], 0.0)
        self.assertIsNone(res["f0_movement_removed_from_both"]["gap_explained_share"])

    def test_a_pitch_movement_gap_is_carried_by_the_pitch_movement(self):
        for sig, flatten, impose in ((self.vibrato, 0.8, 0.8), (self.small_vibrato, 0.7, 0.5)):
            res = fa.attribute(sig, self.steady, F0)
            self.assertGreater(res["baseline"]["real"], res["baseline"]["ours"] * 5)  # the moving one has the larger flux
            self.assertGreater(self.share(res, "f0_movement_removed_from_both"), flatten)
            self.assertGreater(self.share(res, "real_f0_movement_put_into_ours"), impose)
            self.assertLess(self.share(res, "level_movement_removed_from_both"), 0.1)
            self.assertLess(self.share(res, "real_level_movement_put_into_ours"), 0.1)

    def test_a_level_movement_gap_is_carried_by_the_level_movement_in_the_body_of_the_note(self):
        res = fa.attribute(self.shimmer, self.steady, F0)
        self.assertGreater(res["baseline"]["gap"], 0.2)
        # the body of the note: the onset and release are its shape and are left alone, so over all frames the share has a ceiling
        self.assertGreater(self.share(res, "level_movement_removed_from_both", steady=True), 0.7)
        self.assertGreater(self.share(res, "real_level_movement_put_into_ours", steady=True), 0.7)
        self.assertLess(self.share(res, "f0_movement_removed_from_both"), 0.1)
        self.assertLess(self.share(res, "real_f0_movement_put_into_ours"), 0.1)

    def test_a_noise_gap_sits_between_the_harmonics_and_neither_movement_carries_it(self):
        res = fa.attribute(self.noisy, self.steady, F0)
        self.assertGreater(res["baseline"]["gap"], 0.5)
        self.assertLess(self.share(res, "f0_movement_removed_from_both"), 0.1)
        self.assertLess(self.share(res, "level_movement_removed_from_both"), 0.1)
        self.assertGreater(res["restricted_to_between_cells"]["real"], 1.3 * res["restricted_to_harmonic_cells"]["real"])
        self.assertLess(res["restricted_to_between_cells"]["ours"], 0.1)  # the clean tone has nothing between its harmonics
        shares = res["cells"]["real"]["all"]
        self.assertAlmostEqual(shares["harmonic"]["live_cell_share"] + shares["between"]["live_cell_share"], 1.0, places=9)
        self.assertGreater(shares["between"]["live_cell_share"], 2 * res["cells"]["ours"]["all"]["between"]["live_cell_share"])  # noise makes cells live between the partials
        for b in ("40-1000", "1000-2000", "2000-4000"):  # the bands the fourteen harmonics of 200 Hz reach; above 2.8 kHz neither signal has any
            self.assertLess(res["real"]["bands"][b]["harmonic_over_between_db_median"], res["ours"]["bands"][b]["harmonic_over_between_db_median"], b)

    def test_the_harmonic_split_is_exact_and_finds_the_residual_where_the_noise_is(self):
        for sig in (self.steady, self.noisy, self.vibrato):
            look = fa.Look(sig, F0)
            harm, resid = fa.split_harmonic(sig, look.f0_inst)
            np.testing.assert_allclose(harm + resid, sig, atol=1e-12, rtol=0)
        # the residual is defined as y - harmonic, so that proves nothing: with a pitch so low that every bin lies on a harmonic the comb
        # passes everything but the DC bin, so the harmonic part must be y, the first samples included (four overlapping windows must
        # cover them: with three they came out a quarter low). The last window of a signal that stops abruptly holds a net DC, which the
        # comb excludes, so the end is not asserted.
        low = np.full(len(self.steady), 20.0)
        h_all, _ = fa.split_harmonic(self.steady, low)
        err = h_all - self.steady
        rms = np.sqrt(np.mean(self.steady ** 2))
        self.assertLess(np.sqrt(np.mean(err[:2000] ** 2)) / rms, 1e-4)
        self.assertLess(np.sqrt(np.mean(err[8000:-8000] ** 2)) / rms, 1e-4)
        db = lambda x, y: 10 * np.log10(np.mean(x * x) / np.mean(y * y))  # noqa: E731
        h, r = fa.split_harmonic(self.steady, fa.Look(self.steady, F0).f0_inst)
        self.assertLess(db(r, h), -30.0)  # a clean tone has little outside its comb: its onset ramp and its abrupt end
        h, r = fa.split_harmonic(self.noisy, fa.Look(self.noisy, F0).f0_inst)
        self.assertTrue(-39.0 < db(r, h) < -30.0, db(r, h))  # -32 dB white noise, a little over half of it outside a comb of this width
        h, r = fa.split_harmonic(self.vibrato, fa.Look(self.vibrato, F0).f0_inst)
        self.assertLess(db(r, h), -30.0)  # the comb follows a +-60 cent vibrato about its 85 ms mean

    def test_the_flux_rises_until_cells_between_the_partials_are_live_and_then_reads_the_level_of_noise(self):
        # A cell that is noise (not a stable partial's skirt) changes by about 3 dB from frame to frame whatever its level, and the flux
        # counts every cell within 60 dB of the peak. A clean tone reads next to nothing; noise 45 dB under it is below the window and
        # changes little; from 36 dB under it the noise cells are live and the flux sits at 2.3-3.1 dB however much more noise there
        # is, where white noise itself reads about 3 dB. So two signals with the same partials can differ by 2 dB of flux with the
        # same spectrum to the eye, and removing one cause of an unstable cell from a signal that has another leaves the flux where it was.
        def mean_flux(sig):
            look = fa.Look(sig, F0)
            return float(np.mean(look.flux[look.sound][20:-20]))

        self.assertLess(mean_flux(self.steady), 0.2)
        self.assertLess(mean_flux(voice_like(noise_db=-45.0)), 1.0)
        for noise_db in (-36.0, -24.0, -8.0):
            self.assertTrue(2.0 < mean_flux(voice_like(noise_db=noise_db)) < 3.3, noise_db)
        white = fa.Look(np.random.default_rng(5).standard_normal(int(2.5 * SR)), F0)
        self.assertTrue(2.8 < float(np.mean(white.flux[white.sound])) < 3.3)

    def test_a_residual_that_is_noise_closes_a_noise_gap_in_the_body_of_the_note(self):
        res = fa.attribute(self.noisy, self.steady, F0)
        self.assertGreater(self.share(res, "real_residual_put_into_ours", steady=True), 0.5)
        self.assertAlmostEqual(res["residual_to_harmonic_power_db"]["real"], -34.0, delta=3.0)

    def test_the_flux_by_cell_level_puts_noise_near_the_floor_and_movement_on_the_partials(self):
        noisy = fa.attribute(self.noisy, self.steady, F0)["level_classes"]["real"]
        moving = fa.attribute(self.small_vibrato, self.steady, F0)["level_classes"]["real"]
        self.assertGreater(noisy["40-60"]["squared_step_share"], 0.9)  # noise cells live near the floor of the window
        self.assertGreater(noisy["40-60"]["positive_step_rms_db"], 2.0)
        self.assertLess(moving["40-60"]["squared_step_share"], noisy["40-60"]["squared_step_share"])
        self.assertGreater(moving["0-20"]["squared_step_share"] + moving["20-40"]["squared_step_share"], 0.3)  # a moving partial: its peak and skirts
        for c in (noisy, moving):
            self.assertAlmostEqual(sum(v["squared_step_share"] for v in c.values()), 1.0, places=9)

    def test_partial_flutter_sees_a_level_movement_and_not_a_small_pitch_movement(self):
        calm = fa.partial_flutter(fa.Look(self.steady, F0))["median_step_db_rms"]
        shim = fa.partial_flutter(fa.Look(self.shimmer, F0))["median_step_db_rms"]
        vib = fa.partial_flutter(fa.Look(self.small_vibrato, F0))["median_step_db_rms"]
        self.assertLess(calm, 0.05)
        self.assertGreater(shim, 0.2)  # 3 dB at 2.5 Hz, as an 85 ms window sees it
        self.assertLess(vib, 0.1)  # +-12 cents stays inside a +-2 bin sum, although its cells step by a couple of dB

    def test_flattening_then_imposing_the_same_movement_leaves_the_flux_alone(self):
        res = fa.attribute(self.vibrato, self.steady, F0)
        c = res["roundtrip_control_real"]
        self.assertLess(abs(c["mean_flux_after"] - c["mean_flux_before"]) / c["mean_flux_before"], 0.05)

    def test_zero_gains_leave_the_signal_alone_and_a_gain_and_its_inverse_cancel(self):
        look = fa.Look(self.steady, F0)
        zero = np.zeros((12, look.frames))
        np.testing.assert_allclose(fa.modulate_partials(self.steady, look.f0_inst, zero, +1), self.steady, atol=1e-12, rtol=0)
        d = np.tile(3.0 * np.sin(np.arange(look.frames) / 9.0), (12, 1))
        back = fa.modulate_partials(fa.modulate_partials(self.steady, look.f0_inst, d, +1), look.f0_inst, d, -1)
        self.assertLess(np.sqrt(np.mean((back - self.steady) ** 2)) / np.sqrt(np.mean(self.steady ** 2)), 0.02)

    def test_a_partial_level_gap_shows_in_the_partial_levels_and_the_control_prices_the_operation(self):
        res = fa.attribute(voice_like(partial_am_db=3.0), self.steady, F0)
        self.assertGreater(res["baseline"]["steady_frames"]["gap"], 0.2)
        self.assertGreater(res["partial_flutter"]["real"]["median_step_db_rms"], 10 * res["partial_flutter"]["ours"]["median_step_db_rms"])
        self.assertLess(self.share(res, "level_movement_removed_from_both", steady=True), 0.3)  # the partials' movements do not add up to a broadband one
        self.assertLess(self.share(res, "f0_movement_removed_from_both", steady=True), 0.1)
        # the operation that puts a partial level movement into the render has a flux of its own, which the constant-gain control shows:
        # the movement counts only for what it adds to the control
        moved = res["real_partial_level_movement_put_into_ours"]["steady_frames"]["ours"]
        control = res["control_constant_partial_gains_put_into_ours"]["steady_frames"]["ours"]
        self.assertGreater(control, 0.1)  # the operation is not free
        self.assertLess(control, 0.6)  # ...but it is small beside the 2-3 dB of noise or pitch movement
        self.assertGreater(moved, control + 0.3)  # the imposed movement shows above it

    def test_the_baseline_is_the_graded_flux(self):
        look = fa.Look(self.noisy, F0)
        n = dc.Note(self.noisy, F0, "sustained", False, "synthetic", None)
        self.assertAlmostEqual(float(np.mean(look.flux[look.sound])), n.scalars()["mean_flux_db"], places=9)


if __name__ == "__main__":
    unittest.main()
