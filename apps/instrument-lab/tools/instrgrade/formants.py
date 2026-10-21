#!/usr/bin/env python3
"""formants: F1-F3 of a sustained vowel-like sound, measured the same way on a render and on recordings.

A voiced sound at f0 samples its vocal tract only at the harmonics k f0, so an LPC fit to the raw waveform
(order ~14) chases a harmonic, not a resonance, when f0 is not small against the formant spacing (a 155 Hz
voice's 6th harmonic sits at 933 Hz). This estimator reads the level of every harmonic, joins them by linear
interpolation in dB into a smooth power spectrum up to 5.5 kHz, and fits an all-pole model (autocorrelation
LPC, order 8) to THAT: the model sees the envelope between harmonics, not their comb.

  uv run --with numpy python3 -I formants.py compare --render R.wav --vowel a --f0 155.56 \\
      --lib voice/library.json --wav-dir refs-wav48-voice --story story2008-table-II-formants.tsv [--t0 0.8 --t1 2.0]

compare prints the render's F1-F3, the same measure on every library note of that vowel (the spread of real
singers: min / median / max over takes, and over the takes within +-3 semitones of the render's f0), and Story's
own measured F1N-F3N for the same speaker where the table has that vowel. The estimator can place a formant no better than about
half the harmonic spacing (f0 / 2) up to f0 ~ 300 Hz, and above that closely spaced F1/F2 (700 and 1200 Hz at 349 Hz)
are not separable from harmonics 2 apart: test_formants.py checks both on synthetic vowels of known formants.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from instrgrade import SR, hz_of, load  # noqa: E402

TOP_HZ = 5500.0
ORDER = 8  # 4 resonances over 0-5.5 kHz: F1-F3 and one for what the envelope cannot place (order 10-12 splits one formant into two poles)
PRE_EMPHASIS = 0.97
MIN_FORMANT_HZ = 250.0
VOWEL_IPA = {"i": "i", "e": "e", "a": "ɑ", "o": "o", "u": "u"}  # Story's table columns for the five VocalSet vowels


def harmonic_levels(x, f0, t0, t1, kmax=None):
    """[(frequency, level dB)] of each harmonic k f0 below TOP_HZ in x[t0:t1], and the shared measuring window."""
    seg = x[int(t0 * SR):int(t1 * SR)]
    if len(seg) < SR // 5:
        raise ValueError(f"formants: only {len(seg) / SR:.2f} s of signal in [{t0}, {t1}]")
    n = 1 << 18
    A = np.abs(np.fft.rfft(seg * np.hanning(len(seg)), n))
    fr = np.fft.rfftfreq(n, 1 / SR)
    out = []
    for k in range(1, int(TOP_HZ // f0) + 1 if kmax is None else kmax + 1):
        m = (fr > k * f0 * (1 - 0.04 / k)) & (fr < k * f0 * (1 + 0.04 / k))  # +-4% of f0 around the line, tighter in cents for large k
        i = int(np.argmax(A * m))
        out.append((float(fr[i]), 20 * math.log10(float(A[i]) + 1e-12)))
    return out


def lpc_from_power(P, order):
    """Autocorrelation-method LPC coefficients a (a[0] = 1) of a one-sided power spectrum P on a uniform grid 0..fs/2."""
    r = np.fft.irfft(P)[: order + 1]
    a = np.zeros(order + 1)
    a[0] = 1.0
    e = r[0]
    for i in range(1, order + 1):
        k = -(r[i] + a[1:i] @ r[i - 1:0:-1]) / e
        a_new = a.copy()
        a_new[1:i] = a[1:i] + k * a[i - 1:0:-1]
        a_new[i] = k
        a, e = a_new, e * (1 - k * k)
    return a


def harmonic_summary(x, f0, t0=0.8, t1=2.0):
    """Estimator-free description of the radiated vowel spectrum, from the harmonic levels alone.

    peak_hz: frequency of the strongest harmonic below 1.2 kHz (where F1 sits, to within one harmonic spacing).
    mid_db, hi_db: mean power of the harmonics in 1.5-2.5 kHz and 2.5-4 kHz relative to that strongest harmonic [dB]:
    how much second- and third-formant structure the sound carries above its first.
    """
    pts = harmonic_levels(x, f0, t0, t1)
    f = np.array([p[0] for p in pts])
    L = np.array([p[1] for p in pts])
    low = f < 1200.0
    i = int(np.argmax(np.where(low, L, -1e9)))

    def band(lo, hi):
        m = (f >= lo) & (f < hi)
        return float(10 * math.log10(np.mean(10 ** ((L[m] - L[i]) / 10)))) if m.any() else float("nan")
    return dict(peak_hz=float(f[i]), mid_db=band(1500.0, 2500.0), hi_db=band(2500.0, 4000.0))


def formants(x, f0, t0=0.8, t1=2.0, n_formants=3):
    """(F1..F3 [Hz], bandwidths [Hz]) of x by the harmonic-envelope all-pole fit (poles under 250 Hz or wider than
    500 Hz dropped). Read it WITH harmonic_summary: an all-pole model puts poles on a smooth spectrum that has only
    one real peak, so a third "formant" is not evidence of a third resonance unless the band levels support it."""
    pts = harmonic_levels(x, f0, t0, t1)
    f = np.array([p[0] for p in pts])
    L = np.array([p[1] for p in pts])
    grid = np.linspace(0.0, TOP_HZ, 512)
    env_db = np.interp(grid, np.concatenate([[0.0], f]), np.concatenate([[L[0]], L]))
    # first-order pre-emphasis, as in speech LPC: +6 dB/octave takes the tilt of source and radiation out of the model
    w = np.pi * grid / TOP_HZ
    P = 10 ** (env_db / 10) * np.abs(1 - PRE_EMPHASIS * np.exp(-1j * w)) ** 2
    a = lpc_from_power(P, ORDER)
    fs = 2 * TOP_HZ
    z = np.roots(a)
    z = z[np.imag(z) > 1e-6]
    freq = np.angle(z) * fs / (2 * math.pi)
    bw = -np.log(np.abs(z)) * fs / math.pi
    keep = (freq > MIN_FORMANT_HZ) & (bw < 500.0)
    order = np.argsort(freq[keep])
    fr, bws = freq[keep][order], bw[keep][order]
    return fr[:n_formants], bws[:n_formants]


def story_table(path):
    """{vowel letter: (F1N, F2N, F3N, F0)} measured means from Story (2008) Table II, parsed from the repo's TSV."""
    lines = [ln.rstrip("\n").split("\t") for ln in Path(path).read_text().splitlines() if ln and not ln.startswith("#")]
    head = lines[0][1:]
    rows = {ln[0].replace(" Sd", ""): ln[1:] for ln in lines[1:]}
    out = {}
    for letter, ipa in VOWEL_IPA.items():
        if ipa in head:
            j = head.index(ipa)
            out[letter] = tuple(float(rows[k][j].split("±")[0]) for k in ("F1N", "F2N", "F3N")) + (float(rows["Fo"][j]),)
    return out


def fmt(fr):
    return "  ".join(f"F{k + 1} {v:5.0f}" for k, v in enumerate(fr)) + ("   (fewer than three poles survive the filter)" if len(fr) < 3 else "")


def cmd_compare(a):
    x = load(Path(a.render))
    fr, bw = formants(x, a.f0, a.t0, a.t1)
    hs = harmonic_summary(x, a.f0, a.t0, a.t1)
    print(f"render {Path(a.render).name}  vowel /{a.vowel}/  f0 {a.f0:.2f} Hz  window {a.t0}-{a.t1} s")
    print(f"  ours:           {fmt(fr)} Hz")
    print(f"                  strongest harmonic below 1.2 kHz {hs['peak_hz']:.0f} Hz; harmonics at 1.5-2.5 kHz {hs['mid_db']:+.0f} dB and 2.5-4 kHz {hs['hi_db']:+.0f} dB re it")
    lib = [e for e in json.loads(Path(a.lib).read_text()) if e["art"] == a.vowel]
    rows = []
    for e in lib:
        y = load(Path(a.wav_dir) / e["file"])
        t0 = e["start"] + 0.3
        t1 = min(e["end"], t0 + 1.2)
        try:
            ff, _ = formants(y, e["f0_hz"], t0, t1)
            h = harmonic_summary(y, e["f0_hz"], t0, t1)
        except ValueError:
            continue
        rows.append(dict(f0=e["f0_hz"], ff=ff, **h))
    if not rows:
        print("  no library notes of this vowel")
        return
    near = [r for r in rows if abs(12 * math.log2(r["f0"] / a.f0)) <= 3.0]

    def rng(v):
        v = np.array([u for u in v if np.isfinite(u)])
        return "none" if not len(v) else f"{v.min():5.0f} / {np.median(v):5.0f} / {v.max():5.0f}"

    def line(label, rs):
        if not rs:
            return f"  {label}: none"
        return (f"  {label} (n={len(rs)}): F1 {rng([r['ff'][0] for r in rs if len(r['ff']) > 0])}   F2 {rng([r['ff'][1] for r in rs if len(r['ff']) > 1])}"
                f"   F3 {rng([r['ff'][2] for r in rs if len(r['ff']) > 2])}   strongest harmonic {rng([r['peak_hz'] for r in rs])}"
                f"   1.5-2.5 kHz {rng([r['mid_db'] for r in rs])} dB   2.5-4 kHz {rng([r['hi_db'] for r in rs])} dB   (min / median / max)")
    print(line("VocalSet, every take", rows))
    print(line("VocalSet, f0 within 3 semitones", near))
    if a.story:
        st = story_table(a.story).get(a.vowel)
        if st:
            print(f"  Story 2008 measured (the speaker the tract came from, f0 {st[3]:.0f} Hz): F1N {st[0]:.0f}  F2N {st[1]:.0f}  F3N {st[2]:.0f} Hz")
            print(f"  ours minus Story: " + "  ".join(f"F{k + 1} {fr[k] - st[k]:+.0f}" for k in range(min(3, len(fr)))) + " Hz")
    ref = near if len(near) >= 3 else rows
    for key, label, v in (("mid_db", "harmonics at 1.5-2.5 kHz", hs["mid_db"]), ("hi_db", "harmonics at 2.5-4 kHz", hs["hi_db"])):
        vals = np.array([r[key] for r in ref if np.isfinite(r[key])])
        if len(vals):
            where = "inside" if vals.min() <= v <= vals.max() else f"OUTSIDE by {min(abs(v - vals.min()), abs(v - vals.max())):.0f} dB"
            print(f"  {label}: ours {v:+.0f} dB re the strongest low harmonic against real singers {vals.min():+.0f} to {vals.max():+.0f}: {where}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("compare")
    p.add_argument("--render", required=True)
    p.add_argument("--vowel", required=True)
    p.add_argument("--f0", type=float, required=True, help="the render's fundamental, Hz (measure it first; do not trust the label)")
    p.add_argument("--lib", required=True)
    p.add_argument("--wav-dir", required=True)
    p.add_argument("--story")
    p.add_argument("--t0", type=float, default=0.8)
    p.add_argument("--t1", type=float, default=2.0)
    p.set_defaults(fn=cmd_compare)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
