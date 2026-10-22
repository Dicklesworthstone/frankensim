"""freeze.py on synthetic notes whose movement, noise and partials are known by construction.

uv run --with numpy --with matplotlib python3 -I test_freeze.py -v
Each test prints what it built and what it measured, so a failure log says why.
"""
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import decompose as dc  # noqa: E402
import freeze  # noqa: E402
from flux_attribution import f0_instant  # noqa: E402
from decompose import SR  # noqa: E402

F0 = 220.0


def living(f0=F0, dur=2.5, vib_cents=15.0, vib_hz=5.0, shimmer_db=1.5, noise_db=-30.0, seed=7, B=0.0, decay_db_s=0.0, kmax=30):
    """A note that moves: 1/k partials (stretched by B), a sinusoidal pitch movement, each partial's level wandering on its own
    phase at 3-7 Hz (`shimmer_db` peak), a 30 ms attack, an optional decay, and white noise `noise_db` under the note's RMS.
    Returns (signal, harmonic part, noise)."""
    n = int(dur * SR)
    t = np.arange(n) / SR
    rng = np.random.default_rng(seed)
    cents = vib_cents * np.sin(2 * np.pi * vib_hz * t)
    phase = 2 * np.pi * np.cumsum(f0 * 2 ** (cents / 1200)) / SR
    env = np.clip(t / 0.03, 0, 1) * 10 ** (-decay_db_s * t / 20)
    h = np.zeros(n)
    for k in range(1, kmax + 1):
        stretch = math.sqrt(1 + B * k * k)
        if k * f0 * stretch > 15000:
            break
        wobble = 10 ** (shimmer_db * np.sin(2 * np.pi * rng.uniform(3, 7) * t + rng.uniform(0, 2 * np.pi)) / 20)
        h += wobble * np.sin(k * stretch * phase + 0.7 * k) / k
    h *= env
    h *= 0.2 / np.abs(h).max()
    noise = rng.standard_normal(n) * np.sqrt(np.mean(h * h)) * 10 ** (noise_db / 20)
    return h + noise, h, noise


def mean_flux(y, f0=F0, skip_s=0.1):
    """The graded flux (decompose.flux_curve) averaged over the sounding frames past the first `skip_s`."""
    S, sf = dc.stft_db(y, 4096)
    frames = max(1, (len(y) - 2048) // dc.HOP)
    S = S[:, :frames].astype(np.float64)
    fl = dc.flux_curve(S, sf)
    rms = 10 * np.log10(np.sum(10 ** (S / 10), axis=0) + 1e-20)
    t = np.arange(frames) * dc.HOP / SR
    use = (rms > rms.max() - 40) & (t > skip_s)
    return float(np.mean(fl[use]))


def level_db(y):
    return 10 * math.log10(float(np.mean(y * y)) + 1e-30)


class FreezeTests(unittest.TestCase):
    def test_the_split_is_exact_and_the_residual_is_the_noise(self):
        x, h, noise = living()
        s = freeze.split(x, F0)
        err = float(np.max(np.abs(s.harmonic + s.residual - x)))
        st = s.stats(x)
        print(f"\n  built: 1/k partials, 15 c vibrato at 5 Hz, 1.5 dB shimmer, noise {level_db(noise) - level_db(x):.1f} dB re note")
        print(f"  measured: residual {st['residual_db']:.1f} dB re note, |h + r - x| max {err:.1e}, {st['partials']} partials, f0 sd {st['f0_sd_cents']:.2f} c")
        self.assertLess(err, 1e-12)
        sl = slice(int(0.2 * SR), int(2.3 * SR))  # the brick-wall band edges ring at the ends
        self.assertAlmostEqual(level_db(s.residual[sl]), level_db(noise[sl]), delta=1.0)  # the residual is the noise, to a dB
        self.assertLess(level_db(s.harmonic[sl] - h[sl]) - level_db(h[sl]), -25.0)  # and the harmonic part the partials

    def test_the_f0_track_follows_the_vibrato(self):
        x, _, _ = living(noise_db=-40.0)
        s = freeze.split(x, F0)
        t = np.arange(len(x)) / SR
        truth = 15.0 * np.sin(2 * np.pi * 5.0 * t)
        got = 1200 * np.log2(s.f0_hz / F0)
        ok = (t > 0.2) & (t < 2.3)
        rms = float(np.sqrt(np.mean((got[ok] - truth[ok]) ** 2)))
        print(f"\n  vibrato 15 c at 5 Hz: tracked within {rms:.2f} c rms; median f0 {s.f0_median_hz:.3f} Hz")
        self.assertLess(rms, 2.0)
        self.assertAlmostEqual(s.f0_median_hz, F0, delta=F0 * (2 ** (2 / 1200) - 1))

    def test_a_frozen_note_keeps_its_spectrum_and_loses_its_movement(self):
        x, _, _ = living()
        s = freeze.split(x, F0)
        fx, fh, ff = mean_flux(x), mean_flux(s.harmonic), mean_flux(s.frozen)
        fin = f0_instant(s.frozen, s.f0_median_hz)
        sl = slice(int(0.3 * SR), int(2.2 * SR))
        sd = float(np.std(1200 * np.log2(fin[sl] / s.f0_median_hz)))
        Sx, sf = dc.stft_db(s.harmonic[sl], 8192)
        Sz, _ = dc.stft_db(s.frozen[sl], 8192)
        px, pz = np.mean(10 ** (Sx.astype(float) / 10), axis=1), np.mean(10 ** (Sz.astype(float) / 10), axis=1)
        # each partial's power over a +-0.4 f0 band: the vibrato moves partial 12 by +-23 Hz, across several bins
        bands = [np.abs(sf - k * s.f0_median_hz) < 0.4 * s.f0_median_hz for k in range(1, 13)]
        dk = np.abs(np.array([10 * np.log10(px[b].sum() / pz[b].sum()) for b in bands]))
        print(f"\n  flux (dB per frame): original {fx:.3f}, harmonic part {fh:.3f}, frozen {ff:.3f}; frozen f0 sd {sd:.3f} c; "
              f"partial 1-12 mean level change max {dk.max():.2f} dB")
        self.assertGreater(fx, 0.3)
        self.assertLess(ff, 0.1)  # WhiteGull's T-1 on recordings: 0.00-0.08 frozen against 0.45-1.63 real
        self.assertLess(sd, 0.5)
        self.assertLess(float(dk.max()), 1.0)  # the same partial levels: only the movement went

    def test_stiff_partials_need_their_inharmonicity(self):
        B = 4e-4  # a piano string's order: partial 20 is 7.8 % sharp, far outside a 30 Hz heterodyne band
        x, h, _ = living(f0=110.0, B=B, vib_cents=0.0, shimmer_db=0.0, noise_db=-60.0, decay_db_s=10.0, kmax=60)
        with_b, without = freeze.split(x, 110.0, B=B), freeze.split(x, 110.0, B=0.0)
        rb, r0 = with_b.stats(x)["residual_db"], without.stats(x)["residual_db"]
        print(f"\n  stiff string B {B:g}: residual {rb:.1f} dB re note with B, {r0:.1f} dB with B = 0")
        self.assertLess(rb, -40.0)
        self.assertGreater(r0, rb + 15.0)

    def test_the_attack_and_the_decay_survive_freezing(self):
        x, _, _ = living(decay_db_s=8.0, noise_db=-50.0)
        s = freeze.split(x, F0)

        def attack_ms(y):
            e = np.sqrt(np.convolve(y * y, np.ones(240) / 240, "same"))
            top = e[: int(0.3 * SR)].max()
            return 1000 * (np.argmax(e > 0.9 * top) - np.argmax(e > 0.1 * top)) / SR

        def decay_db_s(y):
            w = int(0.05 * SR)
            n = len(y) // w
            lv = 10 * np.log10(np.mean(y[: n * w].reshape(n, w) ** 2, axis=1) + 1e-30)
            t = (np.arange(n) + 0.5) * w / SR
            sel = (t > 0.5) & (t < 2.3)
            return float(np.polyfit(t[sel], lv[sel], 1)[0])

        ax, af = attack_ms(x), attack_ms(s.frozen)
        dx, df = decay_db_s(x), decay_db_s(s.frozen)
        print(f"\n  attack 10-90 %: original {ax:.1f} ms, frozen {af:.1f} ms; decay: original {dx:.2f} dB/s, frozen {df:.2f} dB/s")
        self.assertLess(abs(af - ax), max(3.0, 0.2 * ax))
        self.assertLess(abs(df - dx), 1.0)


    def test_the_four_tone_nulls_carry_the_notes_pitch_and_shape_and_nothing_else(self):
        x, _, _ = living(decay_db_s=8.0, noise_db=-40.0)
        nulls = freeze.null_signals(x, F0)
        s = freeze.split(x, F0)
        self.assertEqual(tuple(nulls), freeze.NULLS)
        self.assertTrue(np.array_equal(nulls["frozen"], s.frozen))
        w = int(0.05 * SR)
        n = len(x) // w

        def lv(y):
            return 10 * np.log10(np.mean(y[: n * w].reshape(n, w) ** 2, axis=1) + 1e-30)

        t = (np.arange(n) + 0.5) * w / SR
        sel = (t > 0.3) & (t < 2.2)
        nfft = 1 << 19
        fr = np.fft.rfftfreq(nfft, 1 / SR)
        for name in ("sine", "1/n", "odd_1/n", "1/n2"):
            z = nulls[name]
            shape = float(np.max(np.abs(lv(z)[sel] - lv(s.frozen)[sel])))
            A = np.abs(np.fft.rfft(z[sel.argmax() * w:] * np.hanning(len(z) - sel.argmax() * w), nfft))
            peak = float(fr[np.argmax(A)])
            print(f"\n  null {name}: level curve within {shape:.2f} dB of the frozen note's, strongest line {peak:.2f} Hz "
                  f"(note median {s.f0_median_hz:.2f} Hz), flux after the attack {mean_flux(z, skip_s=0.25):.3f} dB")
            self.assertLess(shape, 1.0)
            self.assertAlmostEqual(peak, s.f0_median_hz, delta=0.5)
            self.assertLess(mean_flux(z, skip_s=0.25), 0.1)  # the kept attack is the note's; after it, nothing moves


if __name__ == "__main__":
    unittest.main()
