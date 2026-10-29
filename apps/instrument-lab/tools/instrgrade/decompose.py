#!/usr/bin/env python3
"""decompose: take a render and a real recording of the same note apart, side by side, and read every gap
against how far REAL recordings differ from one another.

  refs     reference-only pairs (no render yet): a real note, its spectrogram, scalogram, partial tracks and metrics,
           each metric with its real-vs-real ruler (other dynamics, other takes, neighbouring notes)
  pair     one render against the library's recording of the same note: the same, plus ours and the difference
  batch    pair for each candidate in a CANDIDATES.json (list of the pair arguments)
  listen   a render with no recording to compare it with: audio, maps and the reason, no gaps
  report   print every graded pair as a table, the worst gap against its real-vs-real ruler first
  diff     for every re-graded note, what moved against the ruler since the grading kept in superseded/

Output is P8's site contract: <out>/<family>/<note>/pair.json (frankensim.instrument-pair.v1) with the aligned,
loudness-matched real.wav / ours.wav and bare spectrogram / scalogram maps. See tools/build_instruments.py.

On trj (defaults point at the library built there; numpy + pillow + matplotlib):
  uv run --with numpy --with pillow --with matplotlib python3 -I decompose.py refs [--family clarinet] [--all] --jobs 16
  uv run --with numpy --with pillow --with matplotlib python3 -I decompose.py pair --render R.wav --family clarinet \
      --note D3 --git-sha SHA --observer "baffled_outlet 0 0 0.5" --recipe data/instruments/descriptions/clarinet/render.sh
The library (--lib) comes from buildlib.py, the 48 kHz mono WAVs (--wav-dir) from convert_refs.py.

Conventions (verisim.py, instrgrade.py): 48 kHz mono, the SAME note-aware high-pass on both signals, onset found
against the recording's noise floor and checked against a second, noise-independent criterion, loudness matched over
the same window. A gap is only ever a number computed here; nothing is typed.
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing
import os
import sys
import time
from datetime import datetime, timezone
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from descriptors import GraderRefusalError, fit_stiff_string_law, track_partials  # noqa: E402
from instrgrade import GRADER_VERSION, SR, describe, distance, highpass, hz_of, load, midi_of, name_of, write_wav  # noqa: E402

DECOMPOSE_VERSION = "decompose-1"
PREROLL_S = 0.010  # both signals are cut this long before the detected onset
MAX_S = 3.0
ALIGN_TOL_MS = 1.5
HOP = SR // 100  # 10 ms: STFT hop and the grid of every curve
F_LO, F_HI = 40.0, 16000.0
SHOW_DB = 60.0  # cells within this many dB of the reference peak count in a spectrogram distance
N_PARTIALS = 12
STRONG_DB = 40.0  # a partial this far under the strongest is not a measurement of its own frequency


class OnsetMisaligned(ValueError):
    """IRPDECOMP_ONSET_MISALIGNED: the two onsets do not coincide after alignment."""


# ------------------------------------------------------------------------------------------ alignment
def env_rms_db(y, win=96, hop=48):
    """RMS envelope over `win` samples on a 1 ms grid, in dB re full scale (power)."""
    c = np.concatenate([[0.0], np.cumsum(y * y)])
    n = max(0, (len(y) - win) // hop + 1)
    idx = np.arange(n) * hop
    ms = (c[idx + win] - c[idx]) / win
    return 10 * np.log10(ms + 1e-18), (idx + win / 2) / SR


def env_window(f0_hz):
    """Envelope window: 2 ms, or two periods of the note if that is longer (a shorter window ripples with the tone)."""
    return 96 if not f0_hz else int(min(0.05 * SR, max(96, 2.0 / f0_hz * SR)))


def preroll_s(f0_hz):
    """Lead-in kept before the onset: 10 ms, or enough that the envelope window centred on the onset fits inside it."""
    return max(PREROLL_S, env_window(f0_hz) / 2 / SR + 0.005)


def floor_of(y):
    """Noise floor [dB re FS power]: 5th percentile of the 20 ms RMS over the whole recording."""
    hop = int(0.02 * SR)
    n = len(y) // hop
    if n < 2:
        return -120.0
    ms = np.mean(y[: n * hop].reshape(n, hop) ** 2, axis=1)
    return float(max(-120.0, 10 * np.log10(np.percentile(ms, 5) + 1e-18)))


def floor_of_file(x, f0_hz, block=1 << 18):
    """The noise floor of a long recording: floor_of() of its high-passed signal, high-passed in blocks.

    One FFT high-pass over a 100 MB range file costs a minute and a gigabyte; the floor is a 5th percentile of
    20 ms frame levels, which a block-wise high-pass (a ringing frame at a block edge is loud, never quiet) gives
    to within a dB or so for a hundredth of the work.
    """
    if len(x) <= 2 * block:
        return floor_of(highpass(x, f0_hz))
    return floor_of(np.concatenate([highpass(x[i:i + block], f0_hz) for i in range(0, len(x), block)]))


STEADY_DROP_DB = 20.0  # the onset is where the envelope first reaches this far under the note's own steady level
STEADY_FRAMES_DB = 40.0  # "sounding" frames for the steady level: within this of the loudest


def steady_level_db(e):
    """The note's own steady level [dB]: the MEDIAN of the sounding frames (those within STEADY_FRAMES_DB of the loudest).

    A median, not a maximum: a maximum moves by a dB or two when the window slides by a few milliseconds (a noise peak,
    an overshoot), and on a steep onset that is milliseconds of onset. An overshoot or a priming step lasting a small
    part of the note does not move the median; a decaying note's median sits part-way down its decay, which only moves
    the threshold lower on an onset that is steep anyway.
    """
    e = np.asarray(e)
    top = float(e.max())
    return float(np.median(e[e > top - STEADY_FRAMES_DB]))


def onset_times(y, floor_db, anchor_s=None, f0_hz=None, steady_db=None):
    """(noise-floor crossing, steady-relative crossing) in seconds from the start of y; the onset is the later of the two.

    A: the envelope crosses floor + 12 dB.  B: it crosses (the signal's own steady level) - 20 dB, so a low priming
    step, a click or the skirt of a slow attack does not count and the same rule applies to a render and to a
    recording. Both are absolute levels, so cutting the signal at the onset and detecting again finds the same
    crossing. A alone fires on dust in a digital-silence render; B alone on a recording whose steady level is
    barely above its noise.

    With no anchor the first crossing is taken (a render). With an anchor (a time inside the note, from
    the library's segmentation, or the expected onset in a cut of such a recording) the crossing is the one that
    precedes it: walk back from the anchor to where the level was last under the threshold. If B never was under it
    within the window, the note begins legato on the previous one and has no level onset: refused, never guessed
    (align() hands a recording that begins at its own onset zeros to start from, so a plucked string is not mistaken
    for that).
    """
    e, t = env_rms_db(y, win=env_window(f0_hz))
    n = len(e)
    if not n:
        raise OnsetMisaligned("IRPDECOMP_NO_SIGNAL: empty signal")
    j0 = 0 if anchor_s is None else min(n - 1, int(anchor_s * 1000))
    pk = steady_level_db(e[j0: j0 + int((MAX_S + 1.0) * 1000)]) if steady_db is None else float(steady_db)
    if pk < floor_db + 20:
        raise OnsetMisaligned(f"IRPDECOMP_NO_SIGNAL: steady level {pk:.1f} dB is under floor {floor_db:.1f} dB + 20")

    def crossing(th, strict):
        above = e > th
        if anchor_s is None:
            idx = np.flatnonzero(above)
            return float(t[idx[0]]) if len(idx) else None
        j = min(n - 1, j0 + 25)
        if not above[j]:
            idx = np.flatnonzero(above[j:])
            return float(t[j + idx[0]]) if len(idx) else None
        i = j
        while i > 0 and above[i - 1]:
            i -= 1
        if i == 0:
            if strict:
                raise OnsetMisaligned(f"IRPDECOMP_NO_QUIET_BEFORE_ONSET: the level never falls under {th:.1f} dB before the note "
                                      "(it begins legato on the previous note)")
            return None
        return float(t[i])

    tb = crossing(pk - STEADY_DROP_DB, True)
    ta = crossing(floor_db + 12.0, False)
    if tb is None:
        raise OnsetMisaligned("IRPDECOMP_NO_SIGNAL: the note never reaches its own steady window")
    return (tb if ta is None else ta), tb


def render_floor(y, f0_hz=None):
    """A render's noise floor is that of its pre-roll: the 20 ms frames before the note first reaches STEADY_DROP_DB under its own
    steady level, less 5 ms. A render that begins at its own onset has no pre-roll but the silence align() puts in front of it,
    so its floor is digital silence (floor_of's -120 dB) and the onset is criterion B's. The quietest 20 ms of the whole render is
    not a floor: on a decaying note they lie inside the note's own decay (a 2 s C6 piano render, steady level -41.9 dB, was refused
    as IRPDECOMP_NO_SIGNAL against a floor of -60.8 dB taken there)."""
    e, _ = env_rms_db(y, win=env_window(f0_hz))
    if not len(e):
        return -120.0
    first = int(np.flatnonzero(e > steady_level_db(e) - STEADY_DROP_DB)[0])
    return floor_of(y[: max(0, first * 48 - int(0.005 * SR))])  # env_rms_db's frame i starts at sample 48 i


def align(raw, f0_hz, floor_db, anchor_s=None, at_start=False):
    """High-pass with the note's corner, find the onset, cut `preroll_s` before it. Returns (y, info).

    A signal that begins at its own onset (a render, or a reference window that starts at the file's first sample:
    `at_start`) is given silence to start from before detection, so that a plucked note's step edge is found by the
    same detector, with the same window, on the signal and on its cut. A window that starts in the middle of a
    file is not padded: a note that begins on the tail of the previous one is refused by onset_times.
    The steady level that B is measured from is that of the cut itself (found by iterating), so slow ramps stay
    idempotent; re-detecting on the cut at the lead-in is the alignment assertion. Times in `info` are in the
    time of `raw`.
    """
    y = highpass(raw, f0_hz)
    pre = preroll_s(f0_hz)
    win = env_window(f0_hz)
    pad = int((pre + win / SR + 0.005) * SR) if (anchor_s is None or at_start) else 0
    if pad:
        y = np.concatenate([np.zeros(pad), y])
        anchor_s = None if anchor_s is None else anchor_s + pad / SR
    pk = None
    for _ in range(4):
        ta, tb = onset_times(y, floor_db, anchor_s, f0_hz, pk)
        t_on = max(ta, tb)
        a = a_in = int(round((t_on - pre) * SR))
        ypad = y
        if a < 0:
            ypad, a = np.concatenate([np.zeros(-a), y]), 0
        cut = ypad[a:a + int(MAX_S * SR)]
        pk_cut = steady_level_db(env_rms_db(cut, win=win)[0][: int((MAX_S + 1.0) * 1000)])
        if pk is not None and abs(pk_cut - pk) < 0.01:
            break
        pk = pk_cut
    ra, rb = onset_times(cut, floor_db, None, f0_hz)
    got = max(ra, rb)
    if anchor_s is not None and abs(got - pre) * 1000 > ALIGN_TOL_MS:
        # A recording's onset was found by walking back from the note to the last dip under the threshold. The first crossing of its cut
        # can be earlier: a touch, a click or a fret noise that reached the threshold, dropped back under it and was followed by the note
        # lies inside the lead-in (guitar Gb4/F5 and violin B5 were refused for it). Re-detect the way the onset was found, on the cut
        # given silence to start from (as above), and accept it only if that lands on the lead-in to the same tolerance.
        padc = int((pre + win / SR + 0.005) * SR)
        try:
            ra, rb = onset_times(np.concatenate([np.zeros(padc), cut]), floor_db, padc / SR + pre, f0_hz, pk_cut)
            got = max(ra, rb) - padc / SR
        except OnsetMisaligned:
            pass  # no better: the first-crossing miss below stands
    if abs(got - pre) * 1000 > ALIGN_TOL_MS:
        raise OnsetMisaligned(f"IRPDECOMP_ONSET_MISALIGNED: cut onset {got * 1000:.2f} ms, expected {pre * 1000:.1f} ms")
    off = pad / SR
    return cut, dict(onset_s=float(t_on) - off, noise_criterion_s=float(ta) - off, peak_criterion_s=float(tb) - off,
                     criteria_spread_ms=float(abs(ta - tb) * 1000), floor_db=float(floor_db), cut_onset_s=float(got),
                     cut_start_s=a_in / SR - off, preroll_s=pre)


# ------------------------------------------------------------------------------- time-frequency maps
def stft_db(y, n):
    """Hann STFT power [dB re 1; a full-scale sine reads -6 dB], centred frames, 10 ms hop."""
    w = np.hanning(n)
    pad = np.concatenate([np.zeros(n // 2), y, np.zeros(n)])
    cols = min(max(1, (len(pad) - n) // HOP), int(len(y) / HOP) + 1)
    idx = np.arange(cols)[:, None] * HOP + np.arange(n)[None, :]
    mag = np.abs(np.fft.rfft(pad[idx] * w, axis=1)) ** 2 / (w.sum() ** 2)
    return (10 * np.log10(mag.T + 1e-20)).astype(np.float32), np.fft.rfftfreq(n, 1 / SR)


def cwt_db(y, per_oct=24, f_lo=60.0, f_hi=F_HI, w0=20.0):
    """Analytic Morlet scalogram [dB, a full-scale sine reads -6 dB]: 24 scales per octave, Gaussian bandwidth f/20
    (about a semitone, so the low partials of a pitched note stay resolved), 10 ms frames.

    Each scale occupies a narrow band of the spectrum and only its envelope power per frame is wanted, so the
    band is inverse-transformed on its own at a reduced sample rate (at least ~500 Hz, enough for 10 ms frames)
    instead of at the full rate: a magnitude is unchanged by moving the band to baseband, and the 10 ms frame mean
    is taken from the integral of that envelope power. About 30 times less work than one full-length inverse FFT
    per scale, which cost 50 ms each on a loaded trj.
    """
    n = len(y)
    nfft = 1 << int(math.ceil(math.log2(n + SR // 4)))
    X = np.fft.rfft(y, nfft)
    f = np.fft.rfftfreq(nfft, 1 / SR)
    fc = f_lo * 2 ** (np.arange(int(per_oct * math.log2(f_hi / f_lo)) + 1) / per_oct)
    frames = n // HOP
    edges = np.arange(frames + 1) * HOP / SR
    out = np.empty((len(fc), frames), np.float32)
    for i, c in enumerate(fc):
        lo, hi = np.searchsorted(f, c * (1 - 6 / w0)), np.searchsorted(f, c * (1 + 6 / w0))
        band = X[lo:hi] * np.exp(-0.5 * ((f[lo:hi] - c) / (c / w0)) ** 2)
        m = max(4096, 1 << int(math.ceil(math.log2(max(hi - lo, 2)))))
        w = np.fft.ifft(band, m)  # numpy normalises by 1/m; the full-rate envelope is (m / nfft) times this
        p = (m / nfft) ** 2 * (w.real ** 2 + w.imag ** 2)
        dt = nfft / (m * SR)  # seconds per decimated sample; sample l sits at l * dt
        cum = np.concatenate([[0.0], np.cumsum(p)]) * dt
        e = np.interp(edges, (np.arange(m + 1) - 0.5) * dt, cum)
        out[i] = (e[1:] - e[:-1]) / (HOP / SR)
    return (10 * np.log10(out + 1e-20)).astype(np.float32), fc


def f0_comb(y, hint_hz, kind="sustained", max_cents=400.0, nh=8):
    """Fundamental of a harmonic note near `hint_hz` [Hz] by a harmonic sum of the AMPLITUDE spectrum.

    Candidates lie on a 1-cent grid within +-max_cents of the hint, so an octave or twelfth error is out of
    range by construction; each scores the sum of the spectral amplitudes at k f, k = 1..nh. Amplitude (not log)
    sums are dominated by the real peaks, which is what a bowed bass with a weak fundamental needs. The winner
    is refined by an amplitude-weighted least-squares fit of f_k / k over the peaks found within 1.2% of k f.
    """
    t0, t1 = (0.05, 0.6) if kind == "decaying" else (0.15, 1.5)
    seg = y[int(t0 * SR):int(t1 * SR)]
    if len(seg) < SR // 10:
        seg = y
    n = 1 << 19
    A = np.abs(np.fft.rfft(seg * np.hanning(len(seg)), n))
    fr = np.fft.rfftfreq(n, 1 / SR)
    grid = hint_hz * 2 ** (np.arange(-int(max_cents), int(max_cents) + 1) / 1200)
    score = sum(np.interp(k * grid, fr, A) for k in range(1, nh + 1) if k * grid[0] < 0.45 * SR)
    f = float(grid[int(np.argmax(score))])
    num = den = 0.0
    for k in range(1, nh + 1):
        if k * f > 0.45 * SR:
            break
        m = (fr > k * f * 0.988) & (fr < k * f * 1.012)
        i = int(np.argmax(A * m))
        if m[i] and A[i] > 0.02 * A.max():
            w = A[i]
            num += w * k * fr[i]
            den += w * k * k
    return num / den if den else f


FUNDAMENTAL_FLOOR_DB = -50.0


def check_fundamental(y, f0, kind, label, named_hz=None):
    """Refuse a note that has no energy at its own measured fundamental. `f0` is the grader's pitch estimate for the signal; when
    `named_hz` (the named note's fundamental) is given the refusal says both, since the estimate of a signal that has locked
    onto a parasitic line is not the note's pitch and a reader takes "the fundamental near X Hz" for the note's.

    A pitched note carries its fundamental within 50 dB of its strongest spectral line (a bowed bass or a
    clarinet's even harmonics are weaker, not absent). A render that locks onto one high-order mode (the flute
    sounding 1558 Hz against a D4 name) has a strongest line 5 or more harmonic numbers above anything at the
    fundamental, so there is no note to compare, whatever the harmonic-sum estimator was able to place.
    """
    t0, t1 = (0.05, 0.6) if kind == "decaying" else (0.15, 1.5)
    seg = y[int(t0 * SR):int(t1 * SR)]
    if len(seg) < SR // 10:
        seg = y
    n = 1 << 18
    A = np.abs(np.fft.rfft(seg * np.hanning(len(seg)), n))
    fr = np.fft.rfftfreq(n, 1 / SR)
    band = (fr > 50.0) & (fr < 10000.0)
    top = float(A[band].max())
    near = (fr > f0 * 0.98) & (fr < f0 * 1.02)
    lvl = 20 * math.log10(float(A[near].max()) / top + 1e-12)
    if lvl < FUNDAMENTAL_FLOOR_DB:
        f_top = float(fr[band][int(np.argmax(A[band]))])
        if named_hz is None:
            raise GraderRefusalError(f"IRPDECOMP_NOT_THE_NAMED_NOTE: {label or 'signal'}: the strongest spectral line is {f_top:.0f} Hz "
                                     f"and there is nothing at the fundamental near {f0:.1f} Hz ({lvl:.0f} dB re that line)")
        raise GraderRefusalError(f"IRPDECOMP_NOT_THE_NAMED_NOTE: {label or 'signal'}: the named note's fundamental is {named_hz:.1f} Hz; the "
                                 f"strongest spectral line is {f_top:.0f} Hz, the grader's pitch estimate for it is {f0:.1f} Hz, and there is "
                                 f"nothing at that estimate ({lvl:.0f} dB re the strongest line)")


def lowpass(x, cap_hz):
    """Raised-cosine FFT low-pass: unity below 0.8 cap, zero above cap. Used when a render is valid only up to a band
    (a receiver admitted to ka <= 20, say): the same filter goes on the render, the recording and every ruler note,
    so the grade says nothing about what the render does not model."""
    x = np.asarray(x, dtype=np.float64)
    n = 1 << int(np.ceil(np.log2(len(x) + SR // 2)))
    X = np.fft.rfft(x, n)
    f = np.fft.rfftfreq(n, 1 / SR)
    g = np.clip((cap_hz - f) / (0.2 * cap_hz), 0.0, 1.0)
    return np.fft.irfft(X * (0.5 - 0.5 * np.cos(np.pi * g)), n)[: len(x)]


def cepstral_envelope(P_db, freqs, f0, grid):
    """Smooth spectral envelope [dB] of the mean sounding spectrum, liftered below 0.6/f0, read on `grid` [Hz]."""
    n = (len(freqs) - 1) * 2
    c = np.fft.irfft(P_db, n)
    mc = max(4, int(0.6 / max(f0, 20.0) * SR))
    w = np.zeros(n)
    w[:mc] = 0.5 + 0.5 * np.cos(np.pi * np.arange(mc) / mc)
    w[n - mc + 1:] = w[1:mc][::-1]
    env = np.fft.rfft(c * w).real
    return np.interp(grid, freqs, env)


def modal_partials(y, f0, P, f, n_modes=N_PARTIALS):
    """Strongest spectral peaks of the first 0.6 s, for bars and plates whose modes are not integer multiples."""
    seg = y[: int(0.6 * SR)]
    nfft = 1 << 17
    S = np.abs(np.fft.rfft(seg * np.hanning(len(seg)), nfft)) ** 2
    fr = np.fft.rfftfreq(nfft, 1 / SR)
    pk = np.flatnonzero((S[1:-1] > S[:-2]) & (S[1:-1] >= S[2:])) + 1
    pk = pk[(fr[pk] > 0.5 * f0) & (fr[pk] < F_HI) & (S[pk] > S.max() * 1e-6)]
    picked = []
    for i in pk[np.argsort(S[pk])[::-1]]:
        if all(abs(fr[i] - q) > 0.03 * q for q in picked):
            picked.append(float(fr[i]))
        if len(picked) >= n_modes:
            break
    return sorted(picked)


# ------------------------------------------------------------------------------------------- one note
def flux_cells(S, sf):
    """The cells the spectral flux is made of, for an STFT map `S` [dB] with bin frequencies `sf`: (positive dB step between
    consecutive frames, whether the step is live, the rows of `S` that are in F_LO..F_HI). A step is live when the louder of
    its two cells is within SHOW_DB of the map's peak."""
    m = (sf > F_LO) & (sf < F_HI)
    X = S[m]
    live = np.maximum(X[:, 1:], X[:, :-1]) > X.max() - SHOW_DB
    return np.where(live, np.maximum(np.diff(X, axis=1), 0.0), 0.0), live, m


def flux_curve(S, sf):
    """Spectral flux [dB] per frame (0 for the first): the RMS over the live cells of the positive frame-to-frame dB step."""
    pos, live, _ = flux_cells(S, sf)
    return np.concatenate([[0.0], np.sqrt(np.sum(pos ** 2, axis=0) / np.maximum(live.sum(axis=0), 1))])


def match_level_db(y, kind):
    """The level a pair matches loudness on [dB re full scale]: mean power over the kind's window of an aligned signal
    (0-0.3 s of a decaying note, 0.1-0.8 s of a sustained one)."""
    w0, w1 = (0.0, 0.3) if kind == "decaying" else (0.1, 0.8)
    seg = y[int(w0 * SR):int(w1 * SR)]
    return float(10 * np.log10(np.mean(seg * seg) + 1e-18)) if len(seg) else float("nan")


class Note:
    """Everything the decomposition needs from one aligned signal, computed once."""

    def __init__(self, y, f0_hz, kind="sustained", inharmonic=False, label="", info=None, max_cents=50.0):
        self.y, self.kind, self.inharmonic, self.label, self.info = y, kind, inharmonic, label, info or {}
        self.f0_hint = f0_hz
        self.dur = len(y) / SR
        w0, w1 = (0.0, 0.3) if kind == "decaying" else (0.1, 0.8)
        self.match_window = [w0, min(w1, self.dur)]
        self.level_db = match_level_db(y, kind)
        f0_m = None if inharmonic else f0_comb(y, f0_hz, kind, 400.0)
        self.desc = describe(y, kind, f0_hz, inharmonic=inharmonic, _skip_q3_heavy=True, f0_override=f0_m)
        if self.desc is None:
            raise GraderRefusalError(f"IRPDECOMP_NO_PITCH: {label or 'signal'} has no pitch near {f0_hz:.1f} Hz")
        self.f0 = float(self.desc["f0"])
        off_c = 1200 * math.log2(self.f0 / f0_hz)
        if abs(off_c) > max_cents:
            raise GraderRefusalError(f"IRPDECOMP_PITCH_MISMATCH: {label or 'signal'} measures {off_c:+.0f} cents from the named note")
        if not inharmonic:
            check_fundamental(y, self.f0, kind, label, f0_hz)
        n = 8192 if self.f0 < 80 else 4096
        self.S, self.sf = stft_db(y, n)
        self.C, self.cf = cwt_db(y)
        # only frames whose whole window lies inside the signal: the zero-padded tail splatters
        self.frames = min(self.S.shape[1], self.C.shape[1], max(1, (len(y) - n // 2) // HOP))
        self.S, self.C = self.S[:, : self.frames], self.C[:, : self.frames]
        self.t = np.arange(self.frames) * HOP / SR
        self._curves()
        self._partials()
        grid = 100 * 2 ** (np.arange(0, int(6 * math.log2(12000 / 100)) + 1) / 6)
        ok = self.rms_db[: self.S.shape[1]] > self.rms_db.max() - 30
        mean_db = 10 * np.log10(np.mean(10 ** (self.S[:, : len(ok)][:, ok] / 10), axis=1) + 1e-20) if ok.any() else self.S[:, 0]
        self.env_grid, self.env_db = grid, cepstral_envelope(mean_db, self.sf, self.f0, grid)
        self.mean_db = mean_db

    def _partials(self):
        d = self.desc
        if self.inharmonic:
            modes = modal_partials(self.y, self.f0, self.S, self.sf)
            self.ratios = np.array([m / self.f0 for m in modes])
            self.pf = np.array(modes)
            tr = []
            for m in modes:
                sel = np.abs(self.sf - m) <= max(0.02 * m, 3 * (self.sf[1] - self.sf[0]))
                tr.append(10 * np.log10(np.sum(10 ** (self.S[sel][:, : self.frames] / 10), axis=0) + 1e-20))
            self.pl = np.array(tr).reshape(len(modes), self.frames)
            self.pk = np.arange(1, len(modes) + 1)
            tot = 10 * np.log10(np.mean(10 ** (self.rms_db[self.sound] / 10)) + 1e-20)
            self.pmean = np.array([10 * np.log10(np.mean(10 ** (p[self.sound] / 10)) + 1e-20) for p in self.pl]) - tot
            self.pdev = np.zeros(len(modes))
            self.pfreq_t = [np.full(self.frames, m) for m in modes]
            self.decay = np.array([_decay_db_s(self.t, p) for p in self.pl]) if self.kind == "decaying" else None
            self.strong = np.ones(len(modes), bool)
            self.B = float(d["inharmonicity_B"])
            return
        t0, t1 = d["windows"]["sustain"]
        ps = track_partials(self.y, self.f0, kmax=min(40, int(15000 // self.f0)), fit_inharmonicity=True, t_window=(t0, t1), f0_fixed=self.f0)
        self.pk = np.array([p.k for p in ps], int)
        self.pf = np.array([p.mean_f_hz for p in ps])
        self.pmean = np.array([p.mean_level_db for p in ps])
        self.ratios = self.pf / self.f0
        self.pl = np.array([np.interp(self.t, p.t_s, p.level_db) for p in ps]).reshape(len(ps), self.frames)
        self.pdev = 1200 * np.log2(self.pf / (self.pk * self.f0)) if len(ps) else np.array([])
        self.pfreq_t = [np.interp(self.t, p.t_s, p.f_hz) for p in ps]
        self.decay = np.array([_decay_db_s(self.t, self.pl[i]) for i in range(len(ps))]) if self.kind == "decaying" else None
        # A partial's frequency is a measurement only if the partial is there: within STRONG_DB of the strongest. Below
        # that the tracker is reporting whichever noise peak lies in its search window (a render with two real
        # partials has "partials" 4-12 at -60 to -90 dB), and a stiff-string fit or a frequency gap built on them is
        # an artifact. B is fitted on the strong partials only, with f0 free: a note whose pitch has not settled has
        # every partial offset together from the measured f0, and a fit that holds f0 would call that stiffness.
        # It is 0 with fewer than three strong partials.
        self.strong = self.pmean > self.pmean.max() - STRONG_DB if len(self.pmean) else np.zeros(0, bool)
        if self.strong.sum() >= 3:
            w = 10 ** ((self.pmean[self.strong] - self.pmean.max()) / 20.0)
            self.B = float(fit_stiff_string_law(self.pk[self.strong].astype(float), self.pf[self.strong], w)[1])
        else:
            self.B = 0.0

    def _curves(self):
        S = self.S[:, : self.frames].astype(np.float64)
        band = (self.sf > 80) & (self.sf < 10000)
        P = 10 ** (S[band] / 10)
        fb = self.sf[band][:, None]
        tot = P.sum(axis=0)
        self.rms_db = 10 * np.log10(np.sum(10 ** (S / 10), axis=0) + 1e-20)
        sound = self.rms_db > self.rms_db.max() - 40
        self.sound = sound
        self.centroid_hz = np.where(sound, (fb * P).sum(axis=0) / (tot + 1e-30), np.nan)
        self.flux_db = flux_curve(S, self.sf)

    def scalars(self):
        d = self.desc
        out = dict(f0_hz=self.f0, inharmonicity_B=self.B, odd_even_db=d["odd_even_db"], hnr_db=d["hnr_db"],
                   attack_ms=d["attack_ms"], centroid_rel=d["centroid_rel"], hi_share_db=d["hi_share_db"],
                   f0_sd_cents=d["f0_sd_cents"], level_db=self.level_db, duration_s=self.dur,
                   mean_centroid_hz=float(np.nanmean(self.centroid_hz)) if self.sound.any() else None,
                   mean_flux_db=float(np.mean(self.flux_db[self.sound])) if self.sound.any() else None)
        for k in ("decay_db_s", "slope1_db_s", "slope2_db_s", "knee_s"):
            if k in d:
                out[k] = d[k]
        dec = [v for v in (self.decay[:6] if self.decay is not None else []) if np.isfinite(v)]
        out["decay_med_db_s"] = float(np.median(dec)) if dec else None
        return out


def _decay_db_s(t, lvl):
    i = int(np.argmax(lvl))
    sel = np.arange(i, len(lvl))
    sel = sel[lvl[sel] > lvl[i] - 50]
    return float(np.polyfit(t[sel], lvl[sel], 1)[0]) if len(sel) >= 6 else float("nan")


# --------------------------------------------------------------------------------------- the gaps
ALL_CATS = ("other_dynamic", "other_take", "neighbour")
B_RESOLUTION = 0.2  # B x 1e-4: smaller differences are inside the stiff-string fit's own scatter
SAME_NOTE = ("other_dynamic", "other_take")  # a neighbouring note is not a spread of THIS note's pitch-absolute measures


def _m(label, unit, gap_unit, note, cats=ALL_CATS, value=None, scale=1.0):
    return dict(label=label, unit=unit, gap_unit=gap_unit, note=note, cats=cats, value=value, scale=scale)


METRICS = OrderedDict([
    ("spec_db", _m("Spectral distance (grader)", "dB", "dB", "instrgrade spec_db: the RMS of the third-octave envelope gap and the harmonic-level gap.")),
    ("env_db", _m("Third-octave band envelope", "dB", "dB", "RMS gap of the mean spectrum in third-octave bands, relative to the total.")),
    ("harm_db", _m("Harmonic level pattern", "dB", "dB", "RMS gap of the tracked partials' mean levels, relative to the total.")),
    ("stft_lsd_db", _m("Spectrogram difference", "dB", "dB", "RMS dB difference of the two STFT maps over cells within 60 dB of the reference peak; onsets aligned, loudness matched.", SAME_NOTE)),
    ("cwt_lsd_db", _m("Scalogram difference", "dB", "dB", "RMS dB difference of the two Morlet scalograms over cells within 60 dB of the reference peak.", SAME_NOTE)),
    ("cepstral_env_db", _m("Spectral envelope (cepstral)", "dB", "dB", "RMS dB gap of the cepstrally smoothed mean spectra, liftered below 0.6/f0, 100 Hz to 12 kHz.")),
    ("partial_traj_db", _m("Partial level trajectories", "dB", "dB", "Median over the first eight partials of the RMS difference of their level-versus-time curves with the mean removed.")),
    ("partial_freq_cents", _m("Partial frequencies", "cents", "cents", "RMS over partials of the difference in frequency deviation from the harmonic series (stretch, inharmonicity).")),
    ("f0_cents", _m("Fundamental frequency", "Hz", "cents", "Distance in cents between the two fitted fundamentals.", SAME_NOTE, "f0_hz")),
    ("inharmonicity_B", _m("Inharmonicity B", "B x 1e-4", "B x 1e-4", "Stiff-string law f_k = k f0 sqrt(1 + B k^2) fitted to the tracked partials; shown in units of 1e-4.", ALL_CATS, "inharmonicity_B", 1e4)),
    ("odd_even_db", _m("Odd/even balance", "dB", "dB", "Level of the odd partials (3, 5, 7, 9) against the even ones (2-10).", ALL_CATS, "odd_even_db")),
    ("hnr_db", _m("Harmonics-to-noise", "dB", "dB", "Energy within 14 Hz of the tracked partials against the energy between them, 80 Hz to 10 kHz, with the pitch held steady first, so vibrato and drift count as the tone.", ALL_CATS, "hnr_db")),
    ("attack_log2", _m("Attack time 10-90%", "ms", "octaves of time", "Time for the 5 ms RMS envelope to rise from 10% to 90% of its steady (or peak) level; the gap is log2 of the ratio.", ALL_CATS, "attack_ms")),
    ("decay_db_s", _m("Partial decay rate", "dB/s", "dB/s", "Median over the first partials of the level slope after each partial's peak (decaying notes only).", ALL_CATS, "decay_med_db_s")),
    ("centroid_st", _m("Spectral centroid over time", "Hz", "semitones", "RMS over time of 12 log2(centroid ratio), sounding frames, 80 Hz to 10 kHz.", ALL_CATS, "mean_centroid_hz")),
    ("flux_db", _m("Spectral flux over time", "dB", "dB", "RMS over time of the difference in positive frame-to-frame spectral change.", ALL_CATS, "mean_flux_db")),
])


def _map_gap(A, B, fa, fb, off):
    """RMS dB difference over cells that count in either map, plus the signed mean (A - B)."""
    A, fa, B, fb = _unify(A, fa, B, fb)
    nf = min(A.shape[1], B.shape[1])
    A, B = A[:, :nf].astype(np.float64) + off, B[:, :nf].astype(np.float64)
    m = (fa > F_LO) & (fa < F_HI)
    A, B = A[m], B[m]
    live = np.maximum(A, B) > B.max() - SHOW_DB
    if not live.any():
        return float("nan"), float("nan")
    d = np.clip(A - B, -SHOW_DB, SHOW_DB)[live]
    return float(np.sqrt(np.mean(d * d))), float(np.mean(d))


def _halve(M):
    """Average a spectrogram's bin pairs in power: an 8192-point map onto the 4096-point grid."""
    n = (M.shape[0] - 1) // 2
    p = 10 ** (M.astype(np.float64) / 10)
    return (10 * np.log10(np.concatenate([(p[0:2 * n:2] + p[1:2 * n:2]) / 2, p[-1:]]) + 1e-20)).astype(np.float32)


def _unify(ma, fa, mb, fb):
    """Put two maps on one frequency grid (the coarser), and return them with that grid."""
    if ma.shape[0] > mb.shape[0]:
        ma, fa = _halve(ma), fb
    elif mb.shape[0] > ma.shape[0]:
        mb, fb = _halve(mb), fa
    return ma, fa, mb, fb


def _partial_match(a, b):
    """Index pairs (ia, ib) of partials present in both: same integer k, or the nearest mode ratio within 4%."""
    if a.inharmonic or b.inharmonic:
        pairs = []
        for i, r in enumerate(a.ratios):
            if not len(b.ratios):
                break
            j = int(np.argmin(np.abs(b.ratios - r)))
            if abs(b.ratios[j] / r - 1) < 0.04:
                pairs.append((i, j))
        return pairs
    kb = {int(k): j for j, k in enumerate(b.pk)}
    return [(i, kb[int(k)]) for i, k in enumerate(a.pk) if int(k) in kb]


def compare(a: Note, b: Note):
    """Gaps between candidate `a` and reference `b` for every metric: {id: (gap >= 0 | None, signed a - b | None)}.

    Loudness: a's absolute levels are moved to b's by the difference of their match-window levels.
    """
    off = b.level_db - a.level_db
    g = {}
    dist = distance(a.desc, b.desc)
    for k in ("spec_db", "env_db", "harm_db"):
        g[k] = (dist[k], None)
    g["stft_lsd_db"] = _map_gap(a.S, b.S, a.sf, b.sf, off)
    g["cwt_lsd_db"] = _map_gap(a.C, b.C, a.cf, b.cf, off)
    de = a.env_db + off - b.env_db
    g["cepstral_env_db"] = (float(np.sqrt(np.mean(de * de))), float(np.mean(de)))
    pairs = _partial_match(a, b)
    nf = min(a.frames, b.frames)
    traj, cents = [], []
    for ia, ib in pairs[:N_PARTIALS]:
        la, lb = a.pl[ia, :nf] + off, b.pl[ib, :nf]
        live = np.maximum(la, lb) > max(la.max(), lb.max()) - 40
        if ia < 8 and live.sum() >= 5 and b.pmean[ib] > b.pmean.max() - 50:
            d = (la[live] - la[live].mean()) - (lb[live] - lb[live].mean())
            traj.append(float(np.sqrt(np.mean(d * d))))
        if a.strong[ia] and b.strong[ib]:
            cents.append(1200 * math.log2(a.ratios[ia] / b.ratios[ib]) if (a.inharmonic or b.inharmonic) else float(a.pdev[ia] - b.pdev[ib]))
    g["partial_traj_db"] = (float(np.median(traj)) if traj else None, None)
    g["partial_freq_cents"] = (float(np.sqrt(np.mean(np.square(cents)))) if cents else None, float(np.mean(cents)) if cents else None)
    sc = 1200 * math.log2(a.f0 / b.f0)
    g["f0_cents"] = (abs(sc), sc)
    g["inharmonicity_B"] = (abs(a.B - b.B), a.B - b.B)
    for k in ("odd_even_db", "hnr_db"):
        g[k] = (abs(a.desc[k] - b.desc[k]), a.desc[k] - b.desc[k])
    r = math.log2(max(a.desc["attack_ms"], 0.02) / max(b.desc["attack_ms"], 0.02))
    g["attack_log2"] = (abs(r), r)
    dd = ([a.decay[i] - b.decay[j] for i, j in pairs[:6] if np.isfinite(a.decay[i]) and np.isfinite(b.decay[j])]
          if a.decay is not None and b.decay is not None else [])
    g["decay_db_s"] = (float(np.median(np.abs(dd))), float(np.median(dd))) if dd else (None, None)
    both = a.sound[:nf] & b.sound[:nf]
    if both.sum() >= 5:
        cs = 12 * np.log2(a.centroid_hz[:nf][both] / b.centroid_hz[:nf][both])
        fl = a.flux_db[:nf][both] - b.flux_db[:nf][both]
        g["centroid_st"] = (float(np.sqrt(np.mean(cs * cs))), float(np.mean(cs)))
        g["flux_db"] = (float(np.sqrt(np.mean(fl * fl))), float(np.mean(fl)))
    else:
        g["centroid_st"] = g["flux_db"] = (None, None)
    return {k: ((v[0] if v[0] is not None and math.isfinite(v[0]) else None), v[1]) for k, v in g.items()}


# ------------------------------------------------------------------------------- reference access
FAMILY = {"clarinet": "clarinet", "altosax": "saxophone", "tenortrombone": "trombone", "doublebass": "double-bass"}
INST = {v: k for k, v in FAMILY.items()}
CANON_ART = ("", "novib", "arco")


def family_of(inst):
    return FAMILY.get(inst, inst)


def canon_art(arts):
    for a in CANON_ART:
        if a in arts:
            return a
    return sorted(arts)[0]


class Refs:
    """Library entries and their audio: raw files cached, notes analysed on demand."""

    def __init__(self, lib_path, wav_dir, refs_json=None):
        self.lib = json.loads(Path(lib_path).read_text())
        self.wav_dir = Path(wav_dir)
        self.refs_json = refs_json  # the manifest: also where a recording chain finds its session (session_chain)
        self.meta, self.external = {}, None
        if refs_json:
            manifest = json.loads(Path(refs_json).read_text())
            self.external = manifest.get("external")  # a recording set that is not in references.json (VocalSet)
            foreign = sorted({e["licence"] for e in manifest["files"] if e.get("licence")} - {"uiowa-mis"})
            if foreign:  # a set under another licence cannot be graded without the attribution its licence requires
                missing = [k for k in ("name", "url", "licence", "licence_url", "credit") if not (self.external or {}).get(k)]
                if missing:
                    raise ValueError(f"IRPDECOMP_SOURCE_MISSING: {refs_json} lists recordings under {foreign} but its `external` block "
                                     f"lacks {missing}: a pair would carry no attribution, and the page refuses it")
            for e in manifest["files"]:
                self.meta[Path(e["path"]).stem.replace(" ", "")] = dict(page=e["page"], chain=e["chain"], id=e["id"], name=Path(e["path"]).name)
        self._raw, self._floor, self._notes = OrderedDict(), {}, OrderedDict()
        self.band_hz = None  # set per pair by run_pair: every signal read through this Refs is low-passed there

    def info_of(self, e):
        return self.meta.get(Path(e["file"]).stem, {})

    def page(self, e):
        return self.info_of(e).get("page")

    def raw(self, fname):
        key = (fname, self.band_hz)
        if key not in self._raw:
            x = load(self.wav_dir / fname)
            self._raw[key] = lowpass(x, self.band_hz) if self.band_hz else x
            while len(self._raw) > 8:
                self._raw.popitem(last=False)
        return self._raw[key]

    def aligned(self, e):
        """(aligned samples, alignment info) of library entry `e`, in the file's own time."""
        f0 = hz_of(e["note"])
        x = self.raw(e["file"])
        fkey = (e["file"], self.band_hz)
        if fkey not in self._floor:
            self._floor[fkey] = floor_of_file(x, f0)
        a = max(0, int((e["start"] - 2.0) * SR))
        b = min(len(x), int(e["end"] * SR) + int(0.1 * SR))
        y, info = align(x[a:b], f0, self._floor[fkey], anchor_s=e["start"] - a / SR, at_start=(a == 0))
        info = dict(info, source_onset_s=a / SR + info["onset_s"], source_cut_start_s=a / SR + info["cut_start_s"])
        return y, info

    def note(self, e) -> Note:
        key = (e["file"], e["start"], self.band_hz)
        if key not in self._notes:
            y, info = self.aligned(e)
            self._notes[key] = Note(y, hz_of(e["note"]), e["kind"], e["inst"] in ("xylophone", "marimba"), f'{e["file"]}@{e["name"]}.{e["dyn"]}', info)
            while len(self._notes) > 160:
                self._notes.popitem(last=False)
        return self._notes[key]

    def find(self, inst, midi, dyn=None, art=None):
        return [e for e in self.lib if e["inst"] == inst and e["note"] == midi and (dyn is None or e["dyn"] == dyn)
                and (art is None or e["art"] == art)]

    def ruler_set(self, t, max_each=4, max_instruments=8):
        """Real recordings to measure against target entry `t`, by kind of difference."""
        same = [e for e in self.lib if e["inst"] == t["inst"] and e["art"] == t["art"] and e is not t]

        def prefer(es):  # the same string and session first
            return sorted(es, key=lambda e: (e["sul"] != t["sul"], self.page(e) != self.page(t), e["file"] != t["file"], e["file"]))
        cats = OrderedDict()
        cats["other_dynamic"] = prefer([e for e in same if e["note"] == t["note"] and e["dyn"] != t["dyn"]])[:max_each]
        cats["other_take"] = prefer([e for e in same if e["note"] == t["note"] and e["dyn"] == t["dyn"] and e["file"] != t["file"]])[:max_each]
        cats["neighbour"] = prefer([e for e in same if e["dyn"] == t["dyn"] and abs(e["note"] - t["note"]) in (1, 2)])[:max_each]
        cats["other_articulation"] = [e for e in self.lib if e["inst"] == t["inst"] and e["art"] != t["art"] and e["kind"] == t["kind"]
                                      and e["note"] == t["note"] and e["dyn"] == t["dyn"]][:2]
        best = {}
        for e in self.lib:
            if e["inst"] != t["inst"] and e["kind"] == t["kind"] and e["dyn"] == t["dyn"] and abs(e["note"] - t["note"]) <= 1:
                k = e["inst"] + e["art"]
                if k not in best or abs(e["note"] - t["note"]) < abs(best[k]["note"] - t["note"]):
                    best[k] = e
        cats["other_instrument"] = list(best.values())[:max_instruments]
        return cats

    @staticmethod
    def label(cat, e, t):
        sul = f' ({e["sul"]})' if e["sul"] and e["sul"] != t["sul"] else ""
        if cat == "other_dynamic":
            return f'{e["dyn"]} against {t["dyn"]}, same note{sul}'
        if cat == "other_take":
            return f'another recording of {t["name"]} {t["dyn"]}{sul or " (" + Path(e["file"]).stem + ")"}'
        if cat == "neighbour":
            return f'neighbouring note {e["name"]}, {e["dyn"]}'
        if cat == "other_articulation":
            return f'{e["art"] or "plain"} against {t["art"] or "plain"}'
        return f'{e["inst"]}{" " + e["art"] if e["art"] else ""} {e["name"]}'

    def ruler(self, t):
        """{category: [{label, with_, gaps: {metric: gap}}]} for target entry `t`: real recordings against real."""
        tn = self.note(t)
        out = OrderedDict()
        for cat, es in self.ruler_set(t).items():
            rows = []
            for e in es:
                try:
                    gp = compare(self.note(e), tn)
                except (OnsetMisaligned, GraderRefusalError, ValueError, IndexError, KeyError):
                    continue
                rows.append(dict(label=self.label(cat, e, t), with_=f'{e["file"]}@{e["name"]}.{e["dyn"]}', gaps={k: v[0] for k, v in gp.items()}))
            out[cat] = rows
        return out


# ------------------------------------------------------------------------------------------- images
def _lut(name):
    from matplotlib import colormaps
    return (colormaps[name](np.linspace(0, 1, 256))[:, :3] * 255).astype(np.uint8)


def log_rows(M, f, f_lo, f_hi, rows=360):
    """Resample a (bins x frames) dB map onto `rows` log-spaced frequencies, high frequency at row 0."""
    fr = f_lo * (f_hi / f_lo) ** (1 - np.arange(rows) / (rows - 1))
    i = np.clip(np.searchsorted(f, fr) - 1, 0, len(f) - 2)
    w = np.clip((fr - f[i]) / (f[i + 1] - f[i]), 0, 1)[:, None]
    return M[i] * (1 - w) + M[i + 1] * w


def save_map(rgb, path, width):
    from PIL import Image
    im = Image.fromarray(rgb).resize((width, rgb.shape[0]), Image.BILINEAR)
    im.save(path, "WEBP", quality=86, method=6)


def map_images(kind, real: Note, ours: Note | None, off, out: Path, f_lo, f_hi, span_db=80.0, diff_db=30.0):
    """Bare spectrogram/scalogram images: real, ours, diff (ours minus real). x = time, y = log frequency."""
    Mr, fr = (real.S, real.sf) if kind == "spectrogram" else (real.C, real.cf)
    maps = {"real": (Mr, fr)}
    if ours is not None:
        maps["ours"] = (ours.S, ours.sf) if kind == "spectrogram" else (ours.C, ours.cf)
    nf = min(m.shape[1] for m, _ in maps.values())
    if ours is not None:
        Mr, fr, Mo, fo = _unify(Mr, fr, maps["ours"][0], maps["ours"][1])
        maps = {"real": (Mr, fr), "ours": (Mo, fo)}
    ref_pk = maps["real"][0][:, :nf][(maps["real"][1] > f_lo) & (maps["real"][1] < f_hi)].max()
    names = {}
    mag = _lut("magma")
    stack = {}
    for side, (M, f) in maps.items():
        o = off if side == "ours" else 0.0
        stack[side] = log_rows(M[:, :nf].astype(np.float64) + o - ref_pk, f, f_lo, f_hi)
        v = np.clip((stack[side] + span_db) / span_db, 0, 1)
        names[side] = f"{kind}-{side}.webp"
        save_map(mag[(v * 255).astype(int)], out / names[side], nf * 4)
    names["ours"] = names.get("ours")
    names["diff"] = None
    if ours is not None:
        D = np.clip(stack["ours"] - stack["real"], -diff_db, diff_db)
        D[np.maximum(stack["ours"], stack["real"]) <= -span_db] = 0.0
        names["diff"] = f"{kind}-diff.webp"
        save_map(_lut("RdBu_r")[((D / diff_db + 1) / 2 * 255).astype(int)], out / names["diff"], nf * 4)
    return dict(real=names["real"], ours=names["ours"], diff=names["diff"], f_lo_hz=f_lo, f_hi_hz=f_hi, span_db=span_db, diff_range_db=diff_db)


# ------------------------------------------------------------------------------------ review figures
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def save_fig(fig, path):
    fig.savefig(path, dpi=100)
    _plt().close(fig)


def fig_review(a: Note, b: Note, off, path):
    """One axes-and-titles sheet for people reading the numbers: maps, partials, envelope, centroid/flux, level."""
    plt = _plt()
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(4, 3, height_ratios=[1.4, 1.4, 1, 1])
    for row, (M, F) in enumerate((("S", "sf"), ("C", "cf"))):
        ma, fa, mb, fb = _unify(getattr(a, M), getattr(a, F), getattr(b, M), getattr(b, F))
        nf = min(ma.shape[1], mb.shape[1])
        keep = fa >= 0.8 * F_LO
        ma, mb, fk = ma[keep, :nf].astype(float), mb[keep, :nf].astype(float), fa[keep]
        pk = mb[(fk > F_LO) & (fk < F_HI)].max()
        t = np.arange(nf) * HOP / SR
        for k, (X, ttl) in enumerate(((mb, "real"), (ma + off, "ours"))):
            ax = fig.add_subplot(gs[row, k])
            ax.pcolormesh(t, fk, X - pk, vmin=-80, vmax=0, cmap="magma", shading="auto")
            ax.set_yscale("log"); ax.set_ylim(F_LO, F_HI); ax.set_title(f'{ttl}: {"STFT" if row == 0 else "Morlet scalogram"}', fontsize=9)
        ax = fig.add_subplot(gs[row, 2])
        D = np.clip(ma + off - mb, -30, 30)
        D[np.maximum(ma + off, mb) <= pk - SHOW_DB] = np.nan
        im = ax.pcolormesh(t, fk, D, vmin=-30, vmax=30, cmap="RdBu_r", shading="auto")
        ax.set_yscale("log"); ax.set_ylim(F_LO, F_HI); ax.set_title("ours - real (dB)", fontsize=9)
        fig.colorbar(im, ax=ax, pad=0.01)
    ax = fig.add_subplot(gs[2, 0])
    inh = a.inharmonic or b.inharmonic
    top = b.pmean.max() if len(b.pmean) else 0.0
    for n_, c_, lab, dx in ((b, "#222", "real", -0.12), (a, "#d55e00", "ours", 0.12)):
        x = (n_.ratios if inh else n_.pk) + dx
        ax.vlines(x, -80, n_.pmean - top, color=c_, lw=1.5)
        ax.plot(x, n_.pmean - top, "o", color=c_, ms=3, label=lab)
    ax.set_ylim(-70, 6); ax.set_title("partial levels (dB re strongest real partial)", fontsize=9); ax.legend(fontsize=7)
    ax = fig.add_subplot(gs[2, 1])
    nf = min(a.frames, b.frames)
    for n_, ls, o in ((b, "-", 0.0), (a, "--", off)):
        for rank, (ia, ib) in enumerate(_partial_match(a, b)[:6]):
            i = ib if n_ is b else ia
            ax.plot(n_.t[:nf], n_.pl[i, :nf] + o, ls, color=f"C{rank}", lw=1)
    if len(b.pl):
        ax.set_ylim(b.pl[:, :nf].max() - 70, b.pl[:, :nf].max() + 3)
    ax.set_title("partial level trajectories (solid real, dashed ours)", fontsize=9)
    ax = fig.add_subplot(gs[2, 2])
    for n_, c_, o in ((b, "#222", 0.0), (a, "#d55e00", off)):
        ax.semilogx(n_.sf[1:], n_.mean_db[1:] + o, color=c_, alpha=0.25, lw=0.6)
        ax.semilogx(n_.env_grid, n_.env_db + o, color=c_, lw=1.6)
    ax.set_xlim(80, 14000); ax.set_ylim(b.env_db.max() - 80, b.env_db.max() + 8); ax.set_title("cepstral envelope over the mean spectrum", fontsize=9)
    for k, (attr, ttl) in enumerate((("centroid_hz", "spectral centroid (Hz)"), ("flux_db", "spectral flux (dB / 10 ms)"), ("rms_db", "RMS level (dB)"))):
        ax = fig.add_subplot(gs[3, k])
        for n_, c_, o in ((b, "#222", 0.0), (a, "#d55e00", off if attr == "rms_db" else 0.0)):
            ax.plot(n_.t, getattr(n_, attr)[: len(n_.t)] + o, color=c_, lw=1.1)
        ax.set_title(ttl, fontsize=9)
    fig.tight_layout()
    save_fig(fig, path)


# ------------------------------------------------------------------------------------------ outputs
def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, (np.floating, float)):
        return float(f"{o:.5g}") if math.isfinite(o) else None
    if isinstance(o, np.integer):
        return int(o)
    return o


def wav_sample_rate(path):
    b = Path(path).read_bytes()[:4096]
    i = b.index(b"fmt ") + 8
    return int.from_bytes(b[i + 4:i + 8], "little")


def tracks(n: Note, nf):
    return [dict(n=int(n.pk[i]), f_hz=np.round(n.pfreq_t[i][:nf], 3), level_db=np.round(n.pl[i, :nf], 2)) for i in range(min(len(n.pk), N_PARTIALS))]


def metric_rows(gaps, real: Note, ours: Note | None, ruler):
    rs, rv = real.scalars(), (ours.scalars() if ours else {})
    rows = []
    for mid, spec in METRICS.items():
        k = spec["scale"]
        sc = lambda v: None if v is None else v * k  # noqa: E731
        gp = gaps[mid][0] if gaps else None
        key = spec["value"]
        have = [r["gaps"].get(mid) for c in spec["cats"] for r in ruler.get(c, [])]
        if key is None and not any(v is not None for v in have) and gp is None:
            continue
        if mid == "decay_db_s" and rs.get(key) is None:
            continue
        ruler_rows = [dict(label=r["label"], value=r["gaps"][mid] * k) for c in spec["cats"] for r in ruler.get(c, [])
                      if r["gaps"].get(mid) is not None]
        context = [dict(label=r["label"], value=r["gaps"][mid] * k) for c in ("other_articulation", "other_instrument") for r in ruler.get(c, [])
                   if r["gaps"].get(mid) is not None]
        # the STRICT ruler: real recordings of the SAME note only (other dynamics, other takes); neighbouring notes excluded
        strict_rows = [dict(label=r["label"], value=r["gaps"][mid] * k) for c in spec["cats"] if c in SAME_NOTE for r in ruler.get(c, [])
                       if r["gaps"].get(mid) is not None]
        note = spec["note"]
        if mid == "inharmonicity_B" and ruler_rows and max(r["value"] for r in ruler_rows) < B_RESOLUTION:
            ruler_rows, note = [], note + " Real recordings of this note differ by less than the fit resolves, so no ruler is offered."
        if mid == "inharmonicity_B" and strict_rows and max(r["value"] for r in strict_rows) < B_RESOLUTION:
            strict_rows = []
        gap = sc(gp) if ours else None
        rows.append(dict(id=mid, label=spec["label"], unit=spec["unit"], gap_unit=spec["gap_unit"], real=sc(rs.get(key)) if key else None,
                         ours=sc(rv.get(key)) if key and ours else None, gap=gap,
                         signed_gap=sc(gaps[mid][1]) if gaps else None, ruler=ruler_rows, ruler_strict=strict_rows,
                         inside_loose=_inside(gap, ruler_rows), inside_strict=_inside(gap, strict_rows),
                         other_instruments_and_articulations=context, note=note))
    return rows


def _inside(gap, rows):
    """True/False: is the gap no larger than the largest ruler value? None: no gap (no render) or no ruler to judge by."""
    if gap is None or not rows:
        return None
    return bool(gap <= max(r["value"] for r in rows))


def verdict_counts(metrics):
    """Per-pair counts of the two verdicts: against the loose ruler (this note, other dynamics and takes, neighbouring notes)
    and the strict one (this note only). no_ruler = measures with a gap but nothing to judge it by."""
    out = {"measures": sum(1 for m in metrics if m["gap"] is not None)}
    for kind in ("loose", "strict"):
        v = [m[f"inside_{kind}"] for m in metrics if m["gap"] is not None]
        out[kind] = {"inside": sum(1 for x in v if x is True), "outside": sum(1 for x in v if x is False), "no_ruler": sum(1 for x in v if x is None)}
    return out


def staged(final, build):
    """Run build(stage) with a fresh staging directory inside note directory `final`; when it returns, swap the staged
    files in whole (the earlier files move to superseded/, nothing is deleted or overwritten). If it raises, the
    partial files move to superseded/aborted-<UTC>-<pid>-<n>/ and the earlier grading stays in place."""
    final.mkdir(parents=True, exist_ok=True)
    stage = final / f".staging-{os.getpid()}-{time.time_ns()}"
    stage.mkdir()
    try:
        result = build(stage)
        supersede_existing(final)
        for p in sorted(stage.iterdir()):
            p.rename(final / p.name)
        stage.rmdir()
        return result
    except BaseException:
        if stage.exists():
            if any(stage.iterdir()):
                aside = final / "superseded" / f"aborted-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{os.getpid()}-{time.time_ns() % 10**6}"
                aside.parent.mkdir(parents=True, exist_ok=True)
                stage.rename(aside)
            else:
                stage.rmdir()
        raise


def write_pair(root, family, note_name, ref_e, refs: Refs, ref_y, ref_info, cand_y=None, cand_info=None, ours_meta=None, review=False, cand_meta=None,
               chain=None):
    """Write <root>/<family>/<note>/pair.json and the files it names (the P8 site contract).

    Everything is built in a staging directory inside the note's directory and swapped in only when complete, so a
    note never sits without its pair while the slow ruler is computed, and a failure leaves the earlier grading in
    place. A failed attempt's partial files are moved to superseded/aborted-<UTC>-<pid>/, never deleted. `chain` is the
    receipt of the recording chain the render went through (chain id, floor source, ambient and gap hashes, flags), if any.
    """
    final = Path(root) / family / note_name
    return staged(final, lambda stage: _write_pair_staged(stage, final, family, note_name, ref_e, refs, ref_y, ref_info,
                                                          cand_y, cand_info, ours_meta, review, cand_meta, chain))


def _write_pair_staged(stage, final, family, note_name, ref_e, refs, ref_y, ref_info, cand_y, cand_info, ours_meta, review, cand_meta=None,
                       chain=None):
    f0 = hz_of(ref_e["note"])
    inh = ref_e["inst"] in ("xylophone", "marimba")
    n = len(ref_y) if cand_y is None else min(len(ref_y), len(cand_y))
    real = Note(ref_y[:n], f0, ref_e["kind"], inh, "real", ref_info)
    ours = Note(cand_y[:n], f0, ref_e["kind"], inh, "ours", cand_info, max_cents=400.0) if cand_y is not None else None
    if ours is not None:
        d_on = abs(ours.info["cut_onset_s"] - real.info["cut_onset_s"]) * 1000
        if d_on > ALIGN_TOL_MS:
            raise OnsetMisaligned(f"IRPDECOMP_ONSET_MISALIGNED: ours and real onsets are {d_on:.2f} ms apart after alignment")
    nf = min(real.frames, ours.frames) if ours else real.frames
    dur = nf * HOP / SR
    off = real.level_db - ours.level_db if ours else 0.0
    out = stage
    ns = int(round(dur * SR))
    yr, yo = ref_y[:ns], (cand_y[:ns] * 10 ** (off / 20) if ours else None)
    peak = max(np.abs(yr).max(), np.abs(yo).max() if yo is not None else 0.0)
    g = 0.7 / peak
    write_wav(out / "real.wav", yr * g)
    if yo is not None:
        write_wav(out / "ours.wav", yo * g)
    gaps = compare(ours, real) if ours else None
    ruler = refs.ruler(ref_e)
    metrics = metric_rows(gaps, real, ours, ruler)
    top_hz = float(refs.band_hz) if refs.band_hz else F_HI
    figs = dict(spectrogram=map_images("spectrogram", real, ours, off, out, F_LO, top_hz),
                scalogram=map_images("scalogram", real, ours, off, out, 60.0, top_hz))
    pm = _partial_match(ours, real) if ours else []
    parts = None
    if len(real.pk):
        ot = [dict(n=int(real.pk[ib]), f_hz=np.round(ours.pfreq_t[ia][:nf], 3), level_db=np.round(ours.pl[ia, :nf] + off, 2)) for ia, ib in pm
              if ib < N_PARTIALS] if ours else None
        parts = dict(t_s=np.round(real.t[:nf], 4), f0_hz=dict(real=real.f0, ours=ours.f0 if ours else None), real=tracks(real, nf), ours=ot or None)
    hp = min(50.0, 0.5 * f0)
    pj = dict(
        schema="frankensim.instrument-pair.v1", family=family, note=note_name, dynamic=ref_e["dyn"], duration_s=round(dur, 4),
        real=dict(reference_id=refs.info_of(ref_e).get("id"), file=refs.info_of(ref_e).get("name"),
                  segment_s=[round(ref_info["source_cut_start_s"], 4), round(ref_info["source_cut_start_s"] + dur, 4)],
                  **({"source": dict(refs.external, id=refs.info_of(ref_e).get("id"))} if refs.external else {})),
        ours=dict(wav=str(ours_meta["wav"]), git_sha=ours_meta["git_sha"], observer=ours_meta["observer"], sample_rate=ours_meta["sample_rate"],
                  recipe=ours_meta["recipe"]) if ours else None,
        audio=dict(real="real.wav", ours="ours.wav" if ours else None,
                   level_match=((f"Both clips are low-passed at {refs.band_hz:.0f} Hz (the render is valid only below it). " if refs.band_hz else "")
                                + f"One shared gain of {20 * math.log10(g):+.1f} dB puts the louder clip at 0.70 of full scale; "
                                + (f"ours was set {off:+.1f} dB beforehand so its RMS over {real.match_window[0]:.1f}-{real.match_window[1]:.1f} s "
                                   "equals the real recording's." if ours else "the real recording is otherwise untouched."))),
        alignment=dict(onset_real_s=ref_info["source_onset_s"], onset_ours_s=cand_info["onset_s"] if cand_info else None, highpass_hz=hp,
                       level_offset_db=off if ours else None, band_hz=refs.band_hz, preroll_ms=ref_info["preroll_s"] * 1000, tolerance_ms=ALIGN_TOL_MS,
                       onset_criteria_spread_ms=dict(real=ref_info["criteria_spread_ms"], ours=cand_info["criteria_spread_ms"] if cand_info else None)),
        figures=figs, partials=parts, metrics=metrics,
        verdicts=verdict_counts(metrics) if ours else None,
        ruler_with={c: [r["with_"] for r in rs] for c, rs in ruler.items()},
        decompose=dict(version=DECOMPOSE_VERSION, grader=GRADER_VERSION, library_note=f'{ref_e["file"]}@{ref_e["name"]}.{ref_e["dyn"]}',
                       candidate=cand_meta),
        **({"chain": chain} if chain else {}),
    )
    (out / "pair.json").write_text(json.dumps(clean(pj), indent=1, allow_nan=False))
    if review and ours:
        fig_review(ours, real, off, out / "review.png")
    return pj


def single_map_image(kind, M, f, nf, out, f_lo, f_hi, span_db=80.0):
    """A bare spectrogram/scalogram of one signal, same axes, colours and span as map_images, normalised to its own peak."""
    M = M[:, :nf].astype(np.float64)
    pk = M[(f > f_lo) & (f < f_hi)].max()
    stack = log_rows(M - pk, f, f_lo, f_hi)
    v = np.clip((stack + span_db) / span_db, 0, 1)
    name = f"{kind}-ours.webp"
    save_map(_lut("magma")[(v * 255).astype(int)], out / name, nf * 4)
    return dict(real=None, ours=name, diff=None, f_lo_hz=f_lo, f_hi_hz=f_hi, span_db=span_db)


def write_listening(root, family, slot, y, info, meta, reason, dyn="mf", band_hz=None, full_scale_pa=None):
    """A render with no recording to compare it with (an instrument with no licensed reference): audio, its two maps,
    and the reason there is no reference, in the contract P8 accepts (real null, no_reference.reason, alignment null,
    partials null, metrics []). No gap, no ruler, no claim."""
    if not reason or not str(reason).strip():
        raise ValueError("IRPDECOMP_LISTEN_REASON: a listening-only pair states why there is no reference")
    if full_scale_pa is not None and (not math.isfinite(full_scale_pa) or full_scale_pa <= 0):
        raise ValueError("IRPDECOMP_LISTEN_SCALE: full_scale_pa must be finite positive pascals")
    final = Path(root) / family / slot

    def build(stage):
        n = 8192 if hz_of(midi_of(meta["note"])) < 80 else 4096
        S, sf = stft_db(y, n)
        C, cf = cwt_db(y)
        nf = min(S.shape[1], C.shape[1], max(1, (len(y) - n // 2) // HOP))
        top = float(band_hz) if band_hz else F_HI
        figs = dict(spectrogram=single_map_image("spectrogram", S, sf, nf, stage, F_LO, top),
                    scalogram=single_map_image("scalogram", C, cf, nf, stage, 60.0, top))
        dur = nf * HOP / SR
        yo = y[: int(round(dur * SR))]
        # A material/atmosphere experiment must retain its pressure ratios.
        # The declared native Pa-to-PCM scale is a calibration, never a gain to
        # fit this particular note. The default historical listening view is unchanged.
        g = 1.0 if full_scale_pa is not None else 0.7 / max(float(np.abs(yo).max()), 1e-12)
        write_wav(stage / "ours.wav", yo * g)
        level = (f"The render's level is retained at {full_scale_pa:g} Pa per full scale, without peak normalization. "
                 "Use the same full-scale pressure when comparing material or atmosphere variants."
                 if full_scale_pa is not None else
                 f"The peak is set to 0.70 of full scale ({20 * math.log10(g):+.1f} dB); there is no recording to match it to.")
        pj = dict(schema="frankensim.instrument-pair.v1", family=family, note=slot, dynamic=dyn, duration_s=round(dur, 4), real=None,
                  no_reference=dict(reason=str(reason).strip()),
                  ours=dict(wav=str(meta["wav"]), git_sha=meta["git_sha"], observer=meta["observer"], sample_rate=meta["sample_rate"], recipe=meta["recipe"]),
                  audio=dict(real=None, ours="ours.wav",
                             level_match=((f"The clip is low-passed at {band_hz:.0f} Hz (the render is valid only below it). " if band_hz else "")
                                          + level)),
                  alignment=None, figures=figs, partials=None, metrics=[],
                  decompose=dict(version=DECOMPOSE_VERSION, grader=GRADER_VERSION, mode="listening-only", onset_s=info["onset_s"]))
        if full_scale_pa is not None:
            pj["ours"]["full_scale_pa"] = full_scale_pa
            pj["audio"]["gain_db"] = 0.0
        (stage / "pair.json").write_text(json.dumps(clean(pj), indent=1, allow_nan=False))
        return pj

    return staged(final, build)


def run_listening(c, root):
    """The listening-only twin of run_pair: no reference note is looked up and nothing is compared."""
    f0 = hz_of(midi_of(c["note"]))  # only the high-pass corner and the analysis window length use the nominal pitch
    try:
        raw = load(Path(c["wav"]))
    except OSError as e:
        raise GraderRefusalError(f"IRPDECOMP_RENDER_UNREADABLE: {c['wav']}: {e}") from None
    if c.get("band_hz"):
        raw = lowpass(raw, c["band_hz"])
    y, info = _labelled("render", align, raw, f0, render_floor(highpass(raw, f0), f0))
    meta = dict(wav=c["wav"], git_sha=c["git_sha"], observer=c["observer"], recipe=c["recipe"], sample_rate=wav_sample_rate(c["wav"]), note=c["note"])
    return write_listening(root, c["family"], c.get("slot") or c["note"], y, info, meta, c.get("reason"), c.get("dyn", "mf"), c.get("band_hz"), c.get("full_scale_pa"))


DEFAULTS = dict(lib="/data/tmp/wa-inst/p2-grade/grade/library-v6-light.json", wav_dir="/data/tmp/wa-inst/refs-wav48",
                refs_json="/data/tmp/wa-inst/snap/a64e0c7b37d1/frankensim/data/instruments/references.json",
                out="/data/tmp/wa-inst/analysis")


def pick_entry(refs, family, note, dyn, art, contains=None):
    inst = INST.get(family, family)
    cands = refs.find(inst, midi_of(note), dyn, art)
    if art is None:
        arts = {e["art"] for e in refs.lib if e["inst"] == inst and e["note"] == midi_of(note) and e["dyn"] == dyn}
        if not arts:
            raise SystemExit(f"IRPDECOMP_REFERENCE_MISSING: {family} {note} {dyn}")
        cands = refs.find(inst, midi_of(note), dyn, canon_art(arts))
    if contains:
        exact = [e for e in cands if Path(e["file"]).stem.lower() == contains.lower()]
        cands = exact or [e for e in cands if contains.lower() in e["file"].lower()]
    if not cands:
        raise SystemExit(f"IRPDECOMP_REFERENCE_MISSING: {family} {note} {dyn} {art or ''} {contains or ''}")
    return sorted(cands, key=lambda e: e["file"])[0]


def _labelled(signal, fn, *a, **k):
    """Run an alignment and say WHICH signal it was when the onset check refuses."""
    try:
        return fn(*a, **k)
    except OnsetMisaligned as e:
        raise OnsetMisaligned(f"{e} [{signal}]") from None


_CHAINS = {}


def session_chain(refs, ref_e, refdir):
    """The recording chain (chain.py) of the session reference entry `ref_e` was recorded in, built once per process and session
    from the manifest this Refs was given and the fetched files in `refdir`. A session that has no chain record, or whose floor
    cannot be built, is a ChainError: the pair is refused, never graded without the chain it asked for."""
    import chain as chn

    if not refs.refs_json:
        raise chn.ChainError("IRPCHAIN_UNKNOWN_SESSION: a chain needs --refs-json, the manifest that records the sessions")
    stem = Path(ref_e["file"]).stem
    f = next((f for f in json.loads(Path(refs.refs_json).read_text()).get("files", []) if Path(f["path"]).stem.replace(" ", "") == stem), None)
    if f is None:
        raise chn.ChainError(f"IRPCHAIN_UNKNOWN_SESSION: {refs.refs_json} has no session record for {stem}")
    if f["page"] not in _CHAINS:
        _CHAINS[f["page"]] = chn.Chain.from_references(Path(refs.refs_json), f["page"], Path(refdir), refs.lib)
    return _CHAINS[f["page"]]


def run_pair(refs, c, root):
    refs.band_hz = c.get("band_hz")  # a render valid only up to a band is graded only up to it, signal and ruler alike
    ref_e = c.get("ref_entry") or pick_entry(refs, c["family"], c["note"], c.get("dyn", "mf"), c.get("art"), c.get("ref_contains"))
    cand_y = cand_info = ours_meta = None
    if c.get("wav"):  # the render first: it is cheap, and a render that is not the note should be refused for that
        f0 = hz_of(ref_e["note"])
        try:
            raw = load(Path(c["wav"]))
        except OSError as e:
            raise GraderRefusalError(f"IRPDECOMP_RENDER_UNREADABLE: {c['wav']}: {e}") from None
        if refs.band_hz:
            raw = lowpass(raw, refs.band_hz)
        hp = highpass(raw, f0)
        if ref_e["inst"] not in ("xylophone", "marimba"):
            check_fundamental(hp, f0_comb(hp, f0, ref_e["kind"], 400.0), ref_e["kind"], "ours", f0)
        cand_y, cand_info = _labelled("render", align, raw, f0, render_floor(hp, f0))
        ours_meta = dict(wav=c["wav"], git_sha=c["git_sha"], observer=c["observer"], recipe=c["recipe"], sample_rate=wav_sample_rate(c["wav"]))
    ref_y, ref_info = _labelled(f"reference {ref_e['file']}@{ref_e['name']}", refs.aligned, ref_e)
    receipt = None
    if c.get("wav") and (c.get("chain") is not None or c.get("chain_refdir")):
        # the render as its reference's session would have recorded it: level-matched to the reference as the pair matches it,
        # then the session's measured floor and its files' grid (chain.py); no chain is ever applied by default
        ch = c["chain"] if c.get("chain") is not None else session_chain(refs, ref_e, c["chain_refdir"])
        gain_db = match_level_db(ref_y, ref_e["kind"]) - match_level_db(cand_y, ref_e["kind"])
        y, receipt = ch.apply(raw * 10 ** (gain_db / 20), f'{c["family"]} {c.get("slot") or c["note"]} {c["git_sha"]}')
        receipt = dict({k: v for k, v in receipt.items() if k != "blocks"}, level_match_db=gain_db, blocks=receipt["blocks"])
        raw = lowpass(y, refs.band_hz) if refs.band_hz else y
        cand_y, cand_info = _labelled("render through its session's chain", align, raw, f0, render_floor(highpass(raw, f0), f0))
    # `slot` names the pair's directory and note id when one musical note has several pairs (a vowel at C3: "C3-a")
    cand_meta = dict(note=c["note"], slot=c.get("slot"), art=ref_e["art"], dyn=ref_e["dyn"], band_hz=c.get("band_hz"),
                     ref_contains=Path(ref_e["file"]).stem) if c.get("wav") else None  # what a regrade needs to rebuild the same candidate
    return write_pair(root, c["family"], c.get("slot") or c["note"], ref_e, refs, ref_y, ref_info, cand_y, cand_info, ours_meta, c.get("review", False),
                      cand_meta, receipt)


_REFS = None


def _init(lib, wav, refs_json):
    global _REFS
    _REFS = Refs(lib, wav, refs_json)


# What is a fact about a note. An OSError is not (a moved directory, a full disk): it is the run's trouble and is never
# written into the data as a refusal. A render file that cannot be read IS a fact about the candidate: run_pair turns
# that one OSError into IRPDECOMP_RENDER_UNREADABLE.
REFUSALS = (OnsetMisaligned, GraderRefusalError, SystemExit, ValueError, IndexError, KeyError)
NOTE_REFUSALS = REFUSALS


def supersede_existing(d):
    """Move everything in note directory `d` (bar a superseded/ folder) into superseded/<UTC>/; nothing is deleted or
    overwritten. Returns the folder it moved into, or None if there was nothing to move."""
    old = [p for p in Path(d).iterdir() if p.name != "superseded" and not p.name.startswith(".staging")] if Path(d).is_dir() else []
    if not old:
        return None
    stamp, k = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"), 1
    keep = Path(d) / "superseded" / stamp
    while keep.exists():  # never rename over an earlier superseded file
        k += 1
        keep = Path(d) / "superseded" / f"{stamp}-{k}"
    keep.mkdir(parents=True)
    for p in old:
        p.rename(keep / p.name)
    return keep


def refuse_pair(root, family, note, reason, render_sha=None, render_wav=None, dyn=None):
    """Record that a note could not be graded, without deleting anything.

    Everything already in <root>/<family>/<note>/ (the old pair.json, audio, maps) is MOVED into
    superseded/<UTC>/, and pair.json is rewritten as a refusal record in the same contract: ours null,
    and a top-level `refused` {reason, render_sha, render_wav, time}. The page then says the note was
    refused and why, and leaves it out of the badge.
    """
    d = Path(root) / family / note
    d.mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc)
    supersede_existing(d)
    rec = dict(schema="frankensim.instrument-pair.v1", family=family, note=note, dynamic=dyn, ours=None, partials=None, metrics=[],
               refused=dict(reason=str(reason), render_sha=render_sha, render_wav=str(render_wav) if render_wav else None,
                            time=now.strftime("%Y-%m-%dT%H:%M:%SZ"), decompose_version=DECOMPOSE_VERSION))
    (d / "pair.json").write_text(json.dumps(clean(rec), indent=1, allow_nan=False))
    return rec


def _written_since(root, family, note, epoch):
    """True if the note's pair.json is a reference-only pair written at or after `epoch` (seconds): a restarted run skips it."""
    p = Path(root) / family / note / "pair.json"
    try:
        return p.stat().st_mtime >= float(epoch) and json.loads(p.read_text()).get("ours") is None and not json.loads(p.read_text()).get("refused")
    except (OSError, ValueError):
        return False


def has_render(root, family, note):
    """True if <root>/<family>/<note>/pair.json already speaks for a render: ours filled in, or a refusal record that
    names the render it refused. A reference-only run leaves both alone."""
    try:
        j = json.loads((Path(root) / family / note / "pair.json").read_text())
    except (OSError, ValueError):
        return False
    return j.get("ours") is not None or bool((j.get("refused") or {}).get("render_sha"))


def _ref_job(args):
    """Reference-only pair for the first note of `notes` whose onset and pitch verify; every refusal is reported."""
    c, root = args
    refused = []
    for name in c["notes"]:
        if c.get("since") and _written_since(root, c["family"], name, c["since"]):
            return dict(family=c["family"], art=c["art"], note=name, ok=True, kept="already regenerated since --since", refused=refused)
        if has_render(root, c["family"], name):
            return dict(family=c["family"], art=c["art"], note=name, ok=True, kept="a graded or refused render pair is left untouched", refused=refused)
        try:
            pj = run_pair(_REFS, dict(c, note=name), root)
            return dict(family=c["family"], art=c["art"], note=name, ok=True, metrics=len(pj["metrics"]), refused=refused)
        except NOTE_REFUSALS as e:
            refused.append(dict(note=name, reason=str(e)))
            if (Path(root) / c["family"] / name / "pair.json").exists():  # an older pair of this note is stale: supersede it
                refuse_pair(root, c["family"], name, e, dyn=c.get("dyn"))
        except OSError as e:  # the run's own trouble (disk, a moved directory): reported, never written into the data as a refusal
            refused.append(dict(note=name, reason=f"I/O error, not a refusal: {e}"))
    return dict(family=c["family"], art=c["art"], ok=False, refused=refused)


def cmd_refs(a):
    """Reference-only pairs (ours null): four spread notes per family, mf, canonical articulation. A note whose
    onset or pitch does not verify emits nothing and the next-nearest note stands in; the refusal is printed."""
    refs = Refs(a.lib, a.wav_dir, a.refs_json)
    jobs = []
    for inst in sorted({e["inst"] for e in refs.lib}):
        fam = family_of(inst)
        if a.family and fam not in a.family:
            continue
        es = [e for e in refs.lib if e["inst"] == inst and e["dyn"] == "mf"]
        art = canon_art({e["art"] for e in es})
        notes = sorted({e["note"] for e in es if e["art"] == art})
        if a.notes:  # exactly these notes (names, comma-separated), no fallback to a neighbour
            want = {midi_of(n) for n in a.notes.split(",")}
            centres = [k for k, m in enumerate(notes) if m in want]
        else:
            centres = list(range(len(notes))) if a.all else sorted({min(len(notes) - 1, int(q * len(notes))) for q in (0.1, 0.35, 0.6, 0.85)})
        reach = 0 if (a.all or a.notes) else min(4, max(1, len(notes) // 8))
        for ci in centres:
            order = [ci] + [ci + d * s_ for d in range(1, reach + 1) for s_ in (1, -1)]
            order = [k for k in order if 0 <= k < len(notes)]
            jobs.append((dict(family=fam, art=art, dyn="mf", notes=[name_of(notes[k]) for k in order], review=a.review, since=a.since), a.out))
    if a.jobs > 1:
        with ProcessPoolExecutor(min(a.jobs, 16), mp_context=multiprocessing.get_context("fork"), initializer=_init,
                                 initargs=(a.lib, a.wav_dir, a.refs_json)) as pool:
            # one instrument's notes go to one worker in turn: they share its range files and its analysed ruler notes
            res = list(pool.map(_ref_job, jobs, chunksize=4))
    else:
        _init(a.lib, a.wav_dir, a.refs_json)
        res = [_ref_job(j) for j in jobs]
    for r in res:
        print(json.dumps(r))
    raise SystemExit(0 if any(r["ok"] for r in res) else 1)


def cmd_pair(a):
    refs = Refs(a.lib, a.wav_dir, a.refs_json)
    c = dict(family=a.family, note=a.note, dyn=a.dyn, art=a.art, ref_contains=a.ref_contains, band_hz=a.band_hz, wav=a.render, git_sha=a.git_sha, observer=a.observer, recipe=a.recipe, review=a.review,
             chain_refdir=a.chain_refdir)
    if not all((a.render, a.git_sha, a.observer, a.recipe)):
        raise SystemExit("IRPDECOMP_USAGE: a render needs --git-sha, --observer and --recipe (the page shows them)")
    try:
        pj = run_pair(refs, c, a.out)
    except REFUSALS as e:
        refuse_pair(a.out, a.family, a.note, e, a.git_sha, a.render, a.dyn)
        raise SystemExit(f"{e}\n(recorded as a refusal in {a.out}/{a.family}/{a.note}/pair.json)")
    for m in pj["metrics"]:
        top = max((r["value"] for r in m["ruler"]), default=None)
        verdict = "no ruler" if top is None else ("inside" if m["gap"] is not None and m["gap"] <= top else "outside")
        print(f'{m["id"]:20s} gap {m["gap"] if m["gap"] is None else round(m["gap"], 3)!s:>9} {m["gap_unit"]:16s} real-vs-real max {top if top is None else round(top, 3)!s:>9}  {verdict}')


def _batch_job(args):
    c, root = args
    try:
        pj = run_listening(c, root) if c.get("listening_only") else run_pair(_REFS, c, root)
        return dict(family=c["family"], note=c.get("slot") or c["note"], wav=c["wav"], ok=True, metrics=len(pj["metrics"]))
    except REFUSALS as e:
        refuse_pair(root, c["family"], c.get("slot") or c["note"], e, c.get("git_sha"), c.get("wav"), c.get("dyn"))
        return dict(family=c["family"], note=c.get("slot") or c["note"], wav=c["wav"], ok=False, error=str(e))
    except OSError as e:  # the run's own trouble: reported here, never recorded as a refusal, the note's pair untouched
        return dict(family=c["family"], note=c.get("slot") or c["note"], wav=c["wav"], ok=False, error=f"I/O error, not a refusal: {e}")


def cmd_listen(a):
    c = dict(family=a.family, note=a.note, slot=a.slot, dyn=a.dyn, wav=a.render, git_sha=a.git_sha, observer=a.observer, recipe=a.recipe,
             reason=a.reason, band_hz=a.band_hz, full_scale_pa=a.full_scale_pa)
    try:
        pj = run_listening({k: v for k, v in c.items() if v is not None}, a.out)
    except REFUSALS as e:
        refuse_pair(a.out, a.family, a.slot or a.note, e, a.git_sha, a.render, a.dyn)
        raise SystemExit(f"{e}\n(recorded as a refusal in {a.out}/{a.family}/{a.slot or a.note}/pair.json)")
    print(f'listening-only pair written: {a.family}/{pj["note"]}  {pj["duration_s"]} s  ({pj["no_reference"]["reason"]})')


def cmd_batch(a):
    """Pair for every candidate in a JSON list of {family, note, wav, git_sha, observer, recipe, [dyn, art, ref_contains]}."""
    jobs = [(dict(c, dyn=c.get("dyn", "mf"), **({"chain_refdir": a.chain_refdir} if a.chain_refdir else {})), a.out)
            for c in json.loads(Path(a.candidates).read_text())]
    if a.jobs > 1:
        with ProcessPoolExecutor(min(a.jobs, 16), mp_context=multiprocessing.get_context("fork"), initializer=_init,
                                 initargs=(a.lib, a.wav_dir, a.refs_json)) as pool:
            res = list(pool.map(_batch_job, jobs, chunksize=1))
    else:
        _init(a.lib, a.wav_dir, a.refs_json)
        res = [_batch_job(j) for j in jobs]
    for r in res:
        print(json.dumps(r))
    raise SystemExit(0 if all(r["ok"] for r in res) else 1)


def pair_report(pj):
    """Text table of one pair with a render: each metric's gap against the largest real-vs-real distance, worst first."""
    rows = []
    for m in pj["metrics"]:
        top = max((r["value"] for r in m["ruler"]), default=None)
        ratio = None if m["gap"] is None or not top else m["gap"] / top
        rows.append((ratio, m, top))
    rows.sort(key=lambda r: (r[0] is None, -(r[0] or 0)))
    o = pj["ours"]
    out = [f'{pj["family"]} {pj["note"]} {pj["dynamic"]}  ours {Path(o["wav"]).name} @{o["git_sha"]}  vs {pj["real"]["file"]}  '
           f'(ours set {pj["alignment"]["level_offset_db"]:+.1f} dB; onset {pj["alignment"]["onset_ours_s"]:.3f} s / {pj["alignment"]["onset_real_s"]:.3f} s)',
           f'  {"metric":34s} {"real":>10s} {"ours":>10s} {"gap":>9s}  {"ruler max":>9s} {"gap/max":>7s}']
    for ratio, m, top in rows:
        f = lambda v: "-" if v is None else f"{v:.4g}"  # noqa: E731
        out.append(f'  {m["label"] + " [" + m["gap_unit"] + "]":34.34s} {f(m["real"]):>10s} {f(m["ours"]):>10s} {f(m["gap"]):>9s}  {f(top):>9s} '
                   f'{("x" + format(ratio, ".2f")) if ratio is not None else "no ruler":>7s}')
    return "\n".join(out)


def pair_view(pj):
    """{metric id: (gap, largest real-vs-real distance | None)} of a pair with a render."""
    return {m["id"]: (m["gap"], max((r["value"] for r in m["ruler"]), default=None)) for m in pj["metrics"] if m["gap"] is not None}


def diff_pairs(old, new):
    """What moved between two gradings of the same note: per metric the old and new gap, the old and new gap / ruler
    ratio (None where there is no ruler) and whether the verdict flipped, plus the inside/outside totals of each."""
    vo, vn = pair_view(old), pair_view(new)
    rows, tot = [], {"old": [0, 0], "new": [0, 0]}
    for mid in vn:
        if mid not in vo:
            continue
        (go, to), (gn, tn) = vo[mid], vn[mid]
        ro = None if not to else go / to
        rn = None if not tn else gn / tn
        for key, g, t in (("old", go, to), ("new", gn, tn)):
            if t is not None:
                tot[key][0 if g <= t else 1] += 1
        rows.append(dict(id=mid, label=next(m["label"] for m in new["metrics"] if m["id"] == mid), unit=next(m["gap_unit"] for m in new["metrics"] if m["id"] == mid),
                         gap_old=go, gap_new=gn, ratio_old=ro, ratio_new=rn,
                         flipped=(ro is not None and rn is not None and (ro <= 1.0) != (rn <= 1.0))))
    return rows, tot


def cmd_diff(a):
    """For every graded note that has an earlier grading in superseded/, print what moved against the real-vs-real ruler."""
    root = Path(a.out)
    for p in sorted(root.glob("*/*/pair.json")):
        pj = json.loads(p.read_text())
        if not pj.get("ours") or (a.family and pj["family"] not in a.family):
            continue
        prev = sorted((p.parent / "superseded").glob("*/pair.json")) if (p.parent / "superseded").is_dir() else []
        prev = [q for q in prev if (json.loads(q.read_text()).get("ours") or {}).get("git_sha") != pj["ours"]["git_sha"] and json.loads(q.read_text()).get("ours")]
        if a.against:  # compare with the newest earlier grading of a named render sha, not merely the latest different one
            prev = [q for q in prev if json.loads(q.read_text())["ours"]["git_sha"].startswith(a.against)]
        if not prev:
            continue
        old = json.loads(prev[-1].read_text())
        rows, tot = diff_pairs(old, pj)
        print(f'{pj["family"]} {pj["note"]}: {old["ours"]["git_sha"][:9]} -> {pj["ours"]["git_sha"][:9]}   inside {tot["old"][0]} -> {tot["new"][0]}, outside {tot["old"][1]} -> {tot["new"][1]}')
        for r in sorted(rows, key=lambda r: -abs((r["ratio_new"] or 0) - (r["ratio_old"] or 0))):
            if r["flipped"] or (r["ratio_old"] is not None and r["ratio_new"] is not None and abs(r["ratio_new"] - r["ratio_old"]) >= 0.5):
                arrow = "  FLIPPED" if r["flipped"] else ""
                print(f'    {r["label"]} [{r["unit"]}]: gap {r["gap_old"]:.4g} -> {r["gap_new"]:.4g}   x{r["ratio_old"]:.2f} -> x{r["ratio_new"]:.2f}{arrow}')


def cmd_report(a):
    root = Path(a.out)
    for p in sorted(root.glob("*/*/pair.json")):
        pj = json.loads(p.read_text())
        if pj.get("ours") and (not a.family or pj["family"] in a.family):
            print(pair_report(pj) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--lib", default=DEFAULTS["lib"])
        p.add_argument("--wav-dir", default=DEFAULTS["wav_dir"])
        p.add_argument("--refs-json", default=DEFAULTS["refs_json"])
        p.add_argument("--out", default=DEFAULTS["out"])
        p.add_argument("--review", action="store_true", help="also write review.png, an annotated sheet of every decomposition")
    p = sub.add_parser("pair")
    common(p)
    p.add_argument("--render", required=True)
    p.add_argument("--family", required=True)
    p.add_argument("--note", required=True)
    p.add_argument("--dyn", default="mf")
    p.add_argument("--art")
    p.add_argument("--ref-contains", help="substring of the library file name, to pick the string/session of the reference")
    p.add_argument("--band-hz", type=float, help="grade only up to this frequency: low-pass render, recording and ruler there (a render valid only below it)")
    p.add_argument("--git-sha")
    p.add_argument("--observer")
    p.add_argument("--recipe")
    p.add_argument("--chain-refdir", help="put the render through its reference session's recording chain (chain.py), reading the session's "
                                           "fetched files here (FRANKENSIM_REFERENCE_DIR); there is no default chain")
    p.set_defaults(fn=cmd_pair)
    p = sub.add_parser("listen", help="a render with no recording to compare it with: audio, maps and the reason, no gaps")
    p.add_argument("--render", required=True)
    p.add_argument("--family", required=True)
    p.add_argument("--note", required=True, help="nominal note name (sets the high-pass corner only); nothing is compared with it")
    p.add_argument("--slot", help="directory and note id of the pair when it differs from --note")
    p.add_argument("--dyn", default="mf")
    p.add_argument("--git-sha", required=True)
    p.add_argument("--observer", required=True)
    p.add_argument("--recipe", required=True)
    p.add_argument("--reason", required=True, help="one sentence: why there is no recording to compare with")
    p.add_argument("--band-hz", type=float)
    p.add_argument("--full-scale-pa", type=float, help="retain the native render's amplitude at this declared Pa-to-PCM scale instead of peak-normalizing it")
    p.add_argument("--out", default=DEFAULTS["out"])
    p.set_defaults(fn=cmd_listen)
    p = sub.add_parser("batch")
    common(p)
    p.add_argument("candidates")
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--chain-refdir", help="put the render through its reference session's recording chain (chain.py), reading the session's "
                                           "fetched files here (FRANKENSIM_REFERENCE_DIR); there is no default chain")
    p.set_defaults(fn=cmd_batch)
    p = sub.add_parser("refs")
    common(p)
    p.add_argument("--family", action="append")
    p.add_argument("--all", action="store_true", help="every mf note, not four spread ones")
    p.add_argument("--notes", help="only these notes, e.g. G5,C4 (with --family); no fallback to a neighbouring note")
    p.add_argument("--since", type=float, help="epoch seconds: skip notes whose reference-only pair.json is at least this new (resume a run)")
    p.add_argument("--jobs", type=int, default=1)
    p.set_defaults(fn=cmd_refs)
    p = sub.add_parser("diff", help="for every re-graded note, what moved against the real-vs-real ruler since the grading in superseded/")
    p.add_argument("--out", default=DEFAULTS["out"])
    p.add_argument("--family", action="append")
    p.add_argument("--against", help="sha (prefix) of the earlier grading to compare with")
    p.set_defaults(fn=cmd_diff)
    p = sub.add_parser("report", help="print each graded pair as a table, worst gap against its ruler first")
    p.add_argument("--out", default=DEFAULTS["out"])
    p.add_argument("--family", action="append")
    p.set_defaults(fn=cmd_report)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
