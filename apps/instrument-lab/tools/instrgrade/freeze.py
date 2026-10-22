#!/usr/bin/env python3
"""freeze: a real recording of one note with its life taken out, the frozen-real null baseline.

A recording is split into its partials and everything else, then the partials are put back without their fast movement:

  partials    a_k(t) = the complex envelope of partial k: the recording heterodyned at f_k(t) = k f0(t) sqrt(1 + B k^2) and
              low-passed at bw Hz; h = sum_k 2 Re(a_k(t) e^{i phi_k(t)}) is the harmonic part, r = x - h the residual (bow, breath
              and reed noise, the hammer's knock, the room's tail between partials);
  frozen      hf = sum_k 2 A_k(t) cos(2 pi f_k t + c_k): f0 held at its median over the note, one constant phase c_k per partial,
              A_k = |a_k| smoothed over `smooth_s` (the note's slow shape: swell, decay, slow beating) except over the attack,
              where |a_k| is kept as it is and handed over to the smoothed curve along a raised cosine;
  frozen      ... and hf + r ("frozen + noise") for the comparison of the two causes.

A frozen note keeps the recording's spectrum, its envelope and its attack, so a statistic built on the spectrum alone cannot tell
it from the recording, yet its flux is that of a periodic render (WhiteGull's T-1, 2026-10-21: 0.00-0.08 dB per frame against
0.45-1.63 for the untouched recordings). A valid realism gate must fail it. This is a measurement on recordings: nothing here is
an instrument model, and nothing is fitted to a render.

usage: freeze.py IN.wav F0_HZ OUT_PREFIX [--B B] [--smooth-s S]
       writes OUT_PREFIX-{harmonic,residual,frozen,frozen+noise}.wav and prints the split's statistics.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from flux_attribution import f0_instant  # noqa: E402
from instrgrade import SR, load, write_wav  # noqa: E402

SMOOTH_S = 0.25  # half-width of the Hann that keeps a partial's slow shape: movement faster than ~2 Hz is the life removed
ATTACK_S = 0.08  # the attack is kept as recorded for this long after the onset, then handed over to the smoothed envelope
F_TOP = 16000.0  # partials are taken up to here (the grader's maps stop at 16 kHz)
LOUD_DB = 20.0  # the note's median f0 and each partial's constant phase are taken where the note is within this of its peak


@dataclass
class Split:
    harmonic: np.ndarray  # h: the partials as recorded
    residual: np.ndarray  # r = x - h
    frozen: np.ndarray  # hf: the partials with their fast movement removed
    f0_hz: np.ndarray  # f0(t) at every sample
    f0_median_hz: float
    partials: int
    bw_hz: float
    onset_s: float
    loud: np.ndarray  # samples where the note is within LOUD_DB of its peak
    envelope: np.ndarray  # the frozen note's RMS amplitude at every sample: sqrt(sum_k 2 A_k^2)

    def stats(self, x):
        p = float(np.mean(x * x)) + 1e-30
        cents = 1200 * np.log2(self.f0_hz[self.loud] / self.f0_median_hz)
        return dict(partials=self.partials, bw_hz=self.bw_hz, f0_median_hz=self.f0_median_hz, onset_s=self.onset_s,
                    residual_db=10 * math.log10(float(np.mean(self.residual ** 2)) / p + 1e-30),
                    frozen_vs_harmonic_db=10 * math.log10((float(np.mean(self.frozen ** 2)) + 1e-30) / (float(np.mean(self.harmonic ** 2)) + 1e-30)),
                    f0_sd_cents=float(np.std(cents)) if len(cents) else float("nan"))


def lowpass_complex(z, bw):
    """Zero-phase brick-wall low-pass of a complex signal at `bw` Hz, padded so the transform does not wrap the end onto the start."""
    pad = int(2.0 * SR / bw)
    n = 1 << int(math.ceil(math.log2(len(z) + 2 * pad)))
    Z = np.fft.fft(np.concatenate([np.zeros(pad, complex), z]), n)
    Z[np.abs(np.fft.fftfreq(n, 1 / SR)) > bw] = 0
    return np.fft.ifft(Z)[pad:pad + len(z)]


def smooth(v, seconds):
    """Hann-weighted running mean of half-width `seconds`, normalised at the ends (so a curve is not pulled to zero there).
    By FFT: the kernel is 24001 taps at 0.25 s, too long for a direct convolution on every partial."""
    k = max(1, int(seconds * SR))
    v = np.asarray(v, dtype=np.float64)
    n = len(v)
    w = np.hanning(2 * k + 1)
    L = 1 << int(math.ceil(math.log2(n + 2 * k + 1)))
    W = np.fft.rfft(w, L)
    num = np.fft.irfft(np.fft.rfft(v, L) * W, L)[k:k + n]
    den = np.fft.irfft(np.fft.rfft(np.ones(n), L) * W, L)[k:k + n]
    return num / np.maximum(den, 1e-12 * w.sum())


def stretch_bias(x, f0_hz, B, kmax=6):
    """What f0_instant reads above f0 when the partials are stiff: it averages f_k / k over the first six partials with weights
    k^2 times their power, and f_k / k = f0 sqrt(1 + B k^2). The same average of sqrt(1 + B k^2), with each partial's power over
    the whole note (+-0.45 f0 about it, as f0_instant's bands), is the factor to divide by: 1 for B = 0. At B = 4e-4 it is about
    1.003, which would put partial 60 a full 30 Hz band off."""
    if B <= 0:
        return 1.0
    X = np.abs(np.fft.rfft(x)) ** 2
    fr = np.fft.rfftfreq(len(x), 1 / SR)
    num = den = 0.0
    for k in range(1, min(kmax, max(1, int(3000 // f0_hz))) + 1):
        s = math.sqrt(1.0 + B * k * k)
        w = k * k * float(X[np.abs(fr - k * f0_hz * s) < 0.45 * f0_hz].sum())
        num, den = num + w * s, den + w
    return num / den if den > 0 else 1.0


def onset_of(x):
    """First sample where the 5 ms RMS reaches 30 dB under its maximum (the note is cut close to its onset by every caller)."""
    w = int(0.005 * SR)
    e = np.convolve(x * x, np.ones(w) / w, mode="same")
    return int(np.argmax(e > e.max() * 1e-3))


def split(x, f0_hz, B=0.0, smooth_s=SMOOTH_S, attack_s=ATTACK_S, bw_hz=None):
    """Take one note apart (see the module docstring). `f0_hz` is the note's nominal fundamental, `B` the stiff-string
    inharmonicity of its partials (0 for a wind, a bowed string or a voice, whose partials are locked harmonic)."""
    x = np.asarray(x, dtype=np.float64)
    n = len(x)
    t = np.arange(n) / SR
    bw = float(bw_hz) if bw_hz else min(30.0, 0.35 * f0_hz)
    f0 = f0_instant(x, f0_hz) / stretch_bias(x, f0_hz, B)  # Hz at every sample, from the phase of the first six harmonics
    w = int(0.02 * SR)
    lvl = np.convolve(x * x, np.ones(w) / w, mode="same")
    loud = lvl > lvl.max() * 10 ** (-LOUD_DB / 10)
    f0_med = float(np.median(f0[loud]))
    f0 = np.where(np.isfinite(f0) & (np.abs(np.log2(f0 / f0_med)) < 0.25), f0, f0_med)  # where no harmonic is heard, hold the median
    o = onset_of(x)
    ramp = np.clip((t - o / SR - attack_s) / attack_s, 0.0, 1.0)
    keep = 0.5 + 0.5 * np.cos(np.pi * ramp)  # 1 through the attack, falling to 0 over the next attack_s
    h, hf, power = np.zeros(n), np.zeros(n), np.zeros(n)
    K = 0
    for k in range(1, int(F_TOP // f0_med) + 1):
        stretch = math.sqrt(1.0 + B * k * k)
        if k * f0_med * stretch >= F_TOP:
            break
        phi = 2 * np.pi * np.cumsum(k * f0 * stretch) / SR
        a = lowpass_complex(x * np.exp(-1j * phi), bw)
        h += 2 * np.real(a * np.exp(1j * phi))
        mag = np.abs(a)
        env = keep * mag + (1 - keep) * smooth(mag, smooth_s)
        c = float(np.angle(np.sum(a[loud] * np.exp(1j * (phi[loud] - 2 * np.pi * k * f0_med * stretch * t[loud])))))
        hf += 2 * env * np.cos(2 * np.pi * k * f0_med * stretch * t + c)
        power += 2 * env * env
        K = k
    return Split(harmonic=h, residual=x - h, frozen=hf, f0_hz=f0, f0_median_hz=f0_med, partials=K, bw_hz=bw, onset_s=o / SR, loud=loud,
                 envelope=np.sqrt(power))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("wav")
    ap.add_argument("f0_hz", type=float)
    ap.add_argument("out_prefix")
    ap.add_argument("--B", type=float, default=0.0, help="stiff-string inharmonicity of the partials (0: harmonic)")
    ap.add_argument("--smooth-s", type=float, default=SMOOTH_S)
    a = ap.parse_args(argv)
    x = load(Path(a.wav))
    s = split(x, a.f0_hz, a.B, a.smooth_s)
    g = 0.7 / max(float(np.abs(x).max()), 1e-12)  # one gain for every part, so their levels compare
    for name, y in (("harmonic", s.harmonic), ("residual", s.residual), ("frozen", s.frozen), ("frozen+noise", s.frozen + s.residual)):
        write_wav(Path(f"{a.out_prefix}-{name}.wav"), y * g)
    print(json.dumps(dict(wav=a.wav, **s.stats(x))))


if __name__ == "__main__":
    main()
