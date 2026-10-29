#!/usr/bin/env python3
"""instrgrade: grade simulated instrument notes against real recordings by spectrogram descriptors.

The question it answers is the owner's acceptance test: is a simulated note's spectrogram
indistinguishable from a well-recorded real instrument playing the same note? Every distance is
reported beside rulers measured only between REAL recordings:
  - same instrument, same note, another dynamic (how far a player's own variation moves it);
  - same instrument, a neighbouring note at the same dynamic;
  - a different instrument at the same note (how far "another instrument" sits).
A candidate that sits farther from its target than another instrument does is not recognisable
as that instrument, whatever its pitch.

Conventions follow apps/instrument-lab/tools/verisim.py: 48 kHz mono float, one channel,
identical high-pass on both signals of a pair, onset relative to the noise floor.
"""
from __future__ import annotations

import json
import math
import re
import struct
import sys
from pathlib import Path

import numpy as np

from descriptors import (  # noqa: E402
    AnalysisTooShort,
    GraderRefusalError,
    MissingBandSplitFrequency,
    MissingNoiseFloorEstimate,
    Partial,
    am_fm_modulation_spectrum,
    band_split_statistics,
    beats,
    blackman_harris,
    compute_window_erb_profiles,
    cyclo_coherence,
    erb_matrix_distance_db,
    estimate_noise_floor_db,
    fit_stiff_string_law,
    hnr_band,
    intonation_fingerprint,
    jitter_and_shimmer,
    null_baselines,
    partial_envelopes,
    partial_onset_order,
    spectral_flux_db,
    track_partials,
    two_slope_decay,
    windows,
)

GRADER_VERSION = "instrgrade-3"
SR = 48000
NOTE = {"C": 0, "Db": 1, "D": 2, "Eb": 3, "E": 4, "F": 5, "Gb": 6, "G": 7, "Ab": 8, "A": 9, "Bb": 10, "B": 11}
NAMES = ["C", "Db", "D", "Eb", "E", "F", "Gb", "G", "Ab", "A", "Bb", "B"]


def midi_of(name: str) -> int:
    m = re.fullmatch(r"([A-G]b?)(-?\d)", name)
    return 12 * (int(m.group(2)) + 1) + NOTE[m.group(1)]


def name_of(midi: int) -> str:
    return f"{NAMES[midi % 12]}{midi // 12 - 1}"


def hz_of(midi: float) -> float:
    return 440.0 * 2 ** ((midi - 69) / 12)


def cents(f: float, ref: float) -> float:
    return 1200 * math.log2(f / ref)


# ----------------------------------------------------------------------------- io
def load(path: Path, channel: int = 0) -> np.ndarray:
    b = Path(path).read_bytes()
    assert b[:4] == b"RIFF" and b[8:12] == b"WAVE", path
    fmt = data = None
    i = 12
    while i + 8 <= len(b):
        cid, size = b[i:i + 4], struct.unpack("<I", b[i + 4:i + 8])[0]
        if cid == b"fmt ":
            fmt = b[i + 8:i + 8 + size]
        elif cid == b"data":
            data = b[i + 8:i + 8 + size]
        i += 8 + size + (size & 1)
    tag, ch, sr, _, _, bits = struct.unpack("<HHIIHH", fmt[:16])
    if tag == 0xFFFE:
        tag = struct.unpack("<H", fmt[24:26])[0]
    if (tag, bits) == (1, 16):
        x = np.frombuffer(data[: len(data) // 2 * 2], "<i2").astype(np.float64) / 32768
    elif (tag, bits) == (1, 24):
        b3 = np.frombuffer(data[: len(data) // 3 * 3], np.uint8).reshape(-1, 3).astype(np.int32)
        v = b3[:, 0] | (b3[:, 1] << 8) | (b3[:, 2] << 16)
        x = np.where(v >= 1 << 23, v - (1 << 24), v).astype(np.float64) / (1 << 23)
    elif (tag, bits) == (3, 32):
        x = np.frombuffer(data[: len(data) // 4 * 4], "<f4").astype(np.float64)
    else:
        raise SystemExit(f"{path}: unsupported WAV tag={tag} bits={bits}")
    x = x[: len(x) // ch * ch].reshape(-1, ch)[:, channel]
    if sr != SR:  # band-limited resample by FFT (renders at 12/24 kHz)
        n_out = int(round(len(x) * SR / sr))
        X = np.fft.rfft(x)
        Y = np.zeros(n_out // 2 + 1, complex)
        m = min(len(X), len(Y))
        Y[:m] = X[:m]
        x = np.fft.irfft(Y, n_out) * (n_out / len(x))
    return x


def write_wav(path: Path, x: np.ndarray, sr: int = SR):
    y = np.asarray(x, "<f4")
    data = y.tobytes()
    Path(path).write_bytes(b"RIFF" + struct.pack("<I", 36 + len(data)) + b"WAVEfmt " +
                           struct.pack("<IHHIIHH", 16, 3, 1, sr, sr * 4, 4, 32) + b"data" + struct.pack("<I", len(data)) + data)


def highpass(x, f0_expected=None, declared_family_corner=50.0, *, f0=None, f1=None):
    """Note-aware raised-cosine FFT high-pass: f1 = min(declared_family_corner, 0.5 * f0_expected)."""
    x = np.asarray(x, dtype=np.float64)
    n = 1 << int(np.ceil(np.log2(len(x) + SR // 2)))
    X = np.fft.rfft(x, n)
    f = np.fft.rfftfreq(n, 1 / SR)
    if f0 is not None and f1 is not None:
        f0_use, f1_use = float(f0), float(f1)
    else:
        if isinstance(f0_expected, str):
            f_note = hz_of(midi_of(f0_expected))
        elif f0_expected is not None and math.isfinite(float(f0_expected)) and float(f0_expected) > 0:
            f_note = float(f0_expected)
        else:
            # Detect low-frequency fundamental (< 100 Hz) if no explicit note hint was supplied
            low_mask = (f >= 20.0) & (f <= 100.0)
            mid_mask = (f > 100.0) & (f <= 2000.0)
            p_spec = np.abs(X) ** 2
            low_pk = float(p_spec[low_mask].max()) if low_mask.any() else 0.0
            mid_pk = float(p_spec[mid_mask].max()) if mid_mask.any() else 0.0
            if low_pk > 0.05 * max(mid_pk, 1e-30):
                f_note = float(f[low_mask][int(np.argmax(p_spec[low_mask]))])
            else:
                f_note = 2.0 * float(declared_family_corner)
        f1_use = min(float(declared_family_corner), 0.5 * f_note)
        f0_use = 0.6 * f1_use
    g = np.clip((f - f0_use) / max(f1_use - f0_use, 1e-6), 0, 1)
    X *= 0.5 - 0.5 * np.cos(np.pi * g)
    return np.fft.irfft(X, n)[: len(x)]


def env_db(x, win=480):
    e = np.sqrt(np.convolve(x * x, np.ones(win) / win, "same"))
    return 20 * np.log10(e + 1e-12)


def onset(x):
    env = np.convolve(np.abs(x), np.ones(96) / 96, "same")
    head = env[: int(0.05 * SR)]
    lo, hi = np.percentile(head, [10, 90])
    noise = float(hi) if env.argmax() > len(head) and hi <= 10 * lo else 0.0
    thr = max(env.max() * 10 ** (-30 / 20), noise * 10 ** (15 / 20))
    return max(0, int(np.argmax(env > thr)) - int(0.005 * SR))


# ------------------------------------------------------------- scale segmentation
def segment_scale(x, min_gap_s=0.18):
    """Split a UIowa chromatic scale (notes separated by silence) into note spans."""
    hop = 480
    n = len(x) // hop
    rms = np.sqrt(np.mean(x[: n * hop].reshape(n, hop) ** 2, axis=1))
    floor = np.percentile(rms, 10)
    thr = max(1e-5, 0.0125 * rms.max(), 4 * floor)
    act = rms > thr
    spans, i = [], 0
    while i < n:
        if not act[i]:
            i += 1
            continue
        j = i
        while j < n:
            if act[j]:
                j += 1
                continue
            k = j
            while k < n and not act[k]:
                k += 1
            if k - j >= int(min_gap_s * SR / hop) or k >= n:
                break
            j = k
        if j - i >= 3:
            spans.append((max(0, i * hop - int(0.02 * SR)), min(len(x), j * hop + int(0.08 * SR))))
        i = j + 1
    return spans


def segment_by_pitch(x, exp):
    """Legato scales without silences: split where the YIN pitch settles on a new note."""
    tr = yin_track(x, fmin=max(25.0, hz_of(exp[0]) * 0.7), fmax=hz_of(exp[-1]) * 1.4)
    if not len(tr):
        return []
    e = env_db(x)
    thr = e.max() - 45
    spans, cur, start = [], None, 0
    for t, f, p in tr:
        lvl = e[min(len(e) - 1, int(t * SR) + 1024)]
        m = int(round(69 + 12 * math.log2(f / 440))) if (np.isfinite(f) and p > 0.6 and lvl > thr) else None
        if m != cur:
            if cur is not None and t - start >= 0.35 and exp[0] <= cur <= exp[-1]:
                spans.append((int(start * SR), int(t * SR) + 2048))
            cur, start = m, t
    if cur is not None and exp[0] <= cur <= exp[-1]:
        spans.append((int(start * SR), len(x)))
    return spans


def segment_monotone(x, exp, min_note_s=0.35):
    """Align the YIN track to the file's known ascending notes `exp` (Viterbi).

    States: leading silence, then note k and the gap after it, visited in order;
    each file plays each note once, so no state can be re-entered. A voiced
    frame costs its squared semitone distance to note k, capped at 1.5 so a
    YIN harmonic confusion (oboe 3f/2, 4f/3) stays inside some note and the
    order decides which; voiced frames cost 2 in a gap, unvoiced frames cost 1
    inside a note and 0 in a gap. No energy threshold enters, so legato bow
    changes, decay re-crossings and gapless loud files segment alike. Notes
    shorter than `min_note_s` are dropped, never invented; `segment_notes`
    verifies every span's pitch.
    """
    tr = yin_track(x, fmin=max(25.0, hz_of(exp[0]) * 0.7), fmax=hz_of(exp[-1]) * 1.4)
    if not len(tr):
        return []
    e = env_db(x)
    lvl = e[np.minimum(len(e) - 1, (tr[:, 0] * SR).astype(int) + 1024)]
    f = tr[:, 1]
    voiced = np.isfinite(f) & (f > 0) & (tr[:, 2] > 0.6) & (lvl > e.max() - 45)
    semis = np.where(voiced, 69 + 12 * np.log2(np.where(voiced, f, 440.0) / 440), np.nan)
    k = len(exp)
    states = 2 * k + 1  # 0 = leading silence; 2j+1 = note j; 2j+2 = gap after note j
    cost = np.empty((len(tr), states))
    cost[:, 0] = np.where(voiced, 2.0, 0.0)
    for j, m in enumerate(exp):
        d = np.minimum((semis - m) ** 2, 1.5)
        cost[:, 2 * j + 1] = np.where(voiced, d, 1.0)
        cost[:, 2 * j + 2] = np.where(voiced, 2.0, 0.0)
    # Predecessors: stay, or advance from the previous state; a note may also
    # follow the previous note directly (legato, no gap).
    acc = np.full(states, np.inf)
    acc[0] = cost[0, 0]
    acc[1] = cost[0, 1]
    back = np.zeros((len(tr), states), dtype=np.int64)
    for t in range(1, len(tr)):
        nxt = np.empty(states)
        for s in range(states):
            best, arg = acc[s], s
            if s >= 1 and acc[s - 1] < best:
                best, arg = acc[s - 1], s - 1
            if s >= 3 and s % 2 == 1 and acc[s - 2] < best:
                best, arg = acc[s - 2], s - 2
            nxt[s] = best + cost[t, s]
            back[t, s] = arg
        acc = nxt
    s = 2 * k if acc[2 * k] <= acc[2 * k - 1] else 2 * k - 1
    path = np.empty(len(tr), dtype=np.int64)
    for t in range(len(tr) - 1, -1, -1):
        path[t] = s
        s = back[t, s]
    notes = []
    for j in range(k):
        idx = np.flatnonzero(path == 2 * j + 1)
        if len(idx) and tr[idx[-1], 0] - tr[idx[0], 0] >= min_note_s:
            a = int(tr[idx[0], 0] * SR)
            b = int(tr[idx[-1], 0] * SR) + 2048
            notes.append((j, (max(0, a - int(0.02 * SR)), min(len(x), b + int(0.08 * SR)))))
    return notes


class ReferenceMismatch(ValueError):
    """IRPGRADE_REFERENCE_MISMATCH: a reference whose audio contradicts its file name.

    `kind` is "pitch" (a span sounds off its note), "note-absent" (the other
    named notes verify, these never sound) or "no-named-note" (none sounds).
    """

    def __init__(self, kind, detail, verified=()):
        super().__init__(f"IRPGRADE_REFERENCE_MISMATCH[{kind}]: {detail}")
        self.kind = kind
        self.detail = detail
        # (midi, span) pairs that did verify at pitch ("note-absent" only).
        self.verified = list(verified)


def span_semitones(x, span, midi):
    """Harmonic-comb f0 of one span relative to `midi` [semitones].

    Subharmonic summation over [0.7, 1.4] x the named pitch, so YIN harmonic
    confusions (an oboe F4 tracked as Bb4) cannot fail a correct span.
    """
    a, b = span
    f, s = spectrum(x[a:b], 0.0, (b - a) / SR)
    return 12 * math.log2(f0_spectral(f, s, hz_of(midi) * 0.7, hz_of(midi) * 1.4) / hz_of(midi))


def sounding_seconds(x, midi, tol=0.5):
    """Total voiced time anywhere in `x` within `tol` semitones of `midi` [s]."""
    tr = yin_track(x, fmin=hz_of(midi) * 0.7, fmax=hz_of(midi) * 1.4)
    if not len(tr):
        return 0.0
    f = tr[(tr[:, 2] > 0.6) & np.isfinite(tr[:, 1]), 1]
    return float(np.sum(np.abs(12 * np.log2(f / hz_of(midi))) <= tol)) * 480 / SR


def segment_notes(x, exp, max_semitones=0.5, min_note_s=0.35):
    """Note spans of a UIowa file whose ascending notes `exp` are known from its name.

    The silence segmenter is kept wherever it finds exactly the expected count;
    otherwise the monotone pitch alignment decides. A full count is accepted
    only if every span sounds its named note within `max_semitones`. A short
    count refuses when measured: named notes that never sound anywhere
    ("note-absent", the rest verified) or a file that sounds none of them
    ("no-named-note", e.g. 96 kHz data under a 44.1 kHz header plays 0.459x
    low). Anything else keeps the historical routing; the mismatch stays visible.
    """
    spans = segment_scale(x)
    if len(spans) != len(exp):
        aligned = segment_monotone(x, exp, min_note_s)
        if len(aligned) == len(exp):
            spans = [s for _, s in aligned]
        else:
            # Harmonic-comb pitch of each silence-segmented span over a wide range.
            # YIN octave errors on a mislabelled file land on named notes, and a
            # comb octave error can too, so a minority of named spans refuses.
            sounding = [69 + 12 * math.log2(f0_spectral(*spectrum(x[a:b], 0.0, (b - a) / SR), 25.0, 2000.0) / 440)
                        for a, b in spans]
            named = sum(any(abs(p - m) <= max_semitones for m in exp) for p in sounding)
            if spans and 2 * named < len(spans):
                raise ReferenceMismatch("no-named-note", [name_of(m) for m in exp])
            present = {j for j, _ in aligned}
            missing = [m for j, m in enumerate(exp) if j not in present]
            if (aligned and all(sounding_seconds(x, m) < min_note_s for m in missing)
                    and all(abs(span_semitones(x, s, exp[j])) <= max_semitones for j, s in aligned)):
                raise ReferenceMismatch("note-absent", [name_of(m) for m in missing],
                                        [(exp[j], s) for j, s in aligned])
            if len(exp) > 2 and len(spans) < len(exp) / 2:
                return segment_by_pitch(x, exp)
            return spans
    deviations = [span_semitones(x, s, m) for s, m in zip(spans, exp)]
    if not all(abs(d) <= max_semitones for d in deviations):
        raise ReferenceMismatch("pitch", deviations)
    return spans


def expected_notes(range_token: str):
    toks = re.findall(r"[A-G]b?\d", range_token)
    if len(toks) == 1:
        return [midi_of(toks[0])]
    a, b = midi_of(toks[0]), midi_of(toks[1])
    return list(range(a, b + 1))


# --------------------------------------------------------------------- pitch (YIN)
def yin_track(x, fmin=40.0, fmax=2000.0, frame=2048, hop=480, thresh=0.15):
    taus_max = int(SR / fmin)
    taus_min = int(SR / fmax)
    frame = max(frame, 2 * taus_max + 2)
    out = []
    for s in range(0, len(x) - frame, hop):
        w = x[s:s + frame]
        if np.sqrt(np.mean(w * w)) < 1e-6:
            out.append((s / SR, np.nan, 0.0))
            continue
        W = frame // 2
        # difference function via autocorrelation identity
        a = w[:W + taus_max]
        e = np.concatenate([[0.0], np.cumsum(a * a)])
        n2 = 1 << int(np.ceil(np.log2(2 * len(a))))
        A = np.fft.rfft(a, n2)
        B = np.fft.rfft(w[:W][::-1], n2)
        corr = np.fft.irfft(A * B, n2)[W - 1:W - 1 + taus_max + 1]
        tau = np.arange(taus_max + 1)
        d = e[W] + (e[W + tau] - e[tau]) - 2 * corr
        d[0] = 0
        cm = d[1:] * np.arange(1, taus_max + 1) / np.maximum(np.cumsum(d[1:]), 1e-20)
        cm = np.concatenate([[1.0], cm])
        cand = np.where(cm[taus_min:] < thresh)[0]
        if len(cand):
            t = cand[0] + taus_min
            while t + 1 < len(cm) and cm[t + 1] < cm[t]:
                t += 1
        else:
            t = taus_min + int(np.argmin(cm[taus_min:]))
        if 1 <= t < len(cm) - 1:
            y0, y1, y2 = cm[t - 1], cm[t], cm[t + 1]
            den = y0 - 2 * y1 + y2
            tf = t + (0.5 * (y0 - y2) / den if abs(den) > 1e-12 else 0.0)
        else:
            tf = float(t)
        out.append((s / SR, SR / tf, 1 - float(cm[t])))
    return np.array(out) if out else np.zeros((0, 3))


# ---------------------------------------------------------------------- descriptors
THIRD_OCT = 1000 * 2 ** (np.arange(-14, 15) / 3)  # 62.5 Hz .. 16 kHz centres


def steady_window(x, kind, noise_floor_db=None, *, allow_unended_decay=True):
    """Analysis window relative to the (already onset-aligned) note."""
    nf = noise_floor_db if noise_floor_db is not None else estimate_noise_floor_db(x)
    w = windows(x, kind, nf, sr=SR, allow_unended_decay=allow_unended_decay)
    return w["sustain"]


def spectrum(x, t0, t1, nfft=1 << 18):
    seg = x[int(t0 * SR):int(t1 * SR)]
    if len(seg) == 0:
        seg = x
    seg = seg * np.hanning(len(seg))
    s = np.abs(np.fft.rfft(seg, nfft)) ** 2
    return np.fft.rfftfreq(nfft, 1 / SR), s


def f0_spectral(f, s, lo, hi, nh=8):
    """Subharmonic summation on the log power spectrum, 1-cent grid in [lo, hi]."""
    lp = np.log10(s + 1e-30)
    grid = lo * 2 ** (np.arange(0, int(1200 * math.log2(hi / lo)) + 1) / 1200)
    score = np.zeros_like(grid)
    for k in range(1, nh + 1):
        score += np.interp(k * grid, f, lp) / k
    # prefer the lowest of near-equal candidates (octave errors upward are the common failure)
    best = grid[int(np.argmax(score))]
    return float(best)


def partial_ratios(f, s, fmin=100.0, fmax=12000.0, n=5):
    """Frequencies of the n strongest spectral peaks (>= 3% apart), ascending, and ratios to the lowest."""
    sel = (f > fmin) & (f < fmax)
    ff, ss = f[sel], s[sel]
    order = np.argsort(ss)[::-1]
    picked = []
    for i in order:
        if all(abs(ff[i] - q) > 0.03 * q for q in picked):
            picked.append(ff[i])
        if len(picked) >= n:
            break
    picked.sort()
    return [float(q) for q in picked], [float(q / picked[0]) for q in picked] if picked else []


def describe(
    x,
    kind="sustained",
    f0_hint=None,
    chord=False,
    inharmonic=False,
    *,
    noise_floor_db=None,
    band_split_hz=None,
    allow_unended_decay=False,
    _skip_q3_heavy=False,
    f0_override=None,
):
    """x is onset-aligned. Returns a dict of descriptors (all dB values relative, not absolute)."""
    x = np.asarray(x, dtype=np.float64)
    nf_db = float(noise_floor_db) if noise_floor_db is not None else estimate_noise_floor_db(x)
    allow_unended = allow_unended_decay or (noise_floor_db is None)
    win_dict = windows(x, kind, nf_db, sr=SR, allow_unended_decay=allow_unended)
    t0, t1 = win_dict["sustain"]
    tr = yin_track(x[int(t0 * SR):int(t1 * SR)], fmin=30.0, fmax=2400.0)
    good = tr[(tr[:, 2] > 0.7) & np.isfinite(tr[:, 1])] if len(tr) else tr
    f0 = float(np.median(good[:, 1])) if len(good) else (f0_hint or float("nan"))
    if f0_override:  # a caller that measured f0 more carefully (decompose.f0_comb) says so; nothing is estimated here
        f0 = float(f0_override)
    elif f0_hint and not chord:  # trust the spectrum within +-7 semitones of the target note
        fs_, ss_ = spectrum(x, t0, t1)
        f0 = f0_spectral(fs_, ss_, f0_hint * 2 ** (-7 / 12), f0_hint * 2 ** (7 / 12))
    if f0 is None or not np.isfinite(f0):
        f0 = f0_hint
    if f0 is None or not np.isfinite(f0) or f0 <= 0:
        return None
    periodicity = float(np.median(tr[:, 2])) if len(tr) else 0.0
    f, s = spectrum(x, t0, t1)
    tot = s[(f > 25) & (f < 16000)].sum() + 1e-30

    kmax = int(min(30, 15000 // f0)) if f0 and f0 > 0 else 0
    partials = track_partials(
        x,
        f0,
        kmax=kmax,
        sr=SR,
        fit_inharmonicity=not (chord or inharmonic),
        t_window=(t0, t1),
        f0_fixed=f0_override,
    )
    if partials and not (chord or inharmonic):
        ks_arr = np.array([p.k for p in partials], dtype=np.float64)
        fs_arr = np.array([p.mean_f_hz for p in partials], dtype=np.float64)
        lv_arr = np.array([p.mean_level_db for p in partials], dtype=np.float64)
        hi_snr = lv_arr >= (float(np.max(lv_arr)) - 55.0)
        if int(hi_snr.sum()) >= 3:
            w_arr = 10.0 ** ((lv_arr[hi_snr] - float(np.max(lv_arr))) / 20.0)
            f0r, b_fit = fit_stiff_string_law(ks_arr[hi_snr], fs_arr[hi_snr], w_arr, f0_fixed=f0_override)
        else:
            f0r, b_fit = float(f0), 0.0
        harm = np.array([p.mean_level_db for p in partials], dtype=np.float64)
    else:
        f0r, b_fit = float(f0), 0.0
        harm = np.array([], dtype=np.float64)

    # Smooth third-octave envelope: Gaussian weights in log2(f) (FWHM = 1/3 octave)
    bands = []
    lf = np.log2(np.maximum(f, 1.0))
    sig = (1 / 3) / 2.3548
    for fc in THIRD_OCT:
        w = np.exp(-0.5 * ((lf - math.log2(fc)) / sig) ** 2)
        bands.append(10 * np.log10((s * w).sum() / (tot * 1.0645 * sig * 3) + 1e-30))
    bands = np.array(bands)

    # Tracked-partial HNR (±3 bin widths around tracked partials, including stiff partials)
    seg_sustain = x[int(t0 * SR):int(t1 * SR)]
    if len(seg_sustain) < 64:
        seg_sustain = x
    if partials:
        hnr = hnr_band(seg_sustain, partials, (80.0, 10000.0), sr=SR)
    else:
        m_fallback = (f > 80) & (f < 10000)
        hm = np.zeros_like(f, bool)
        for k_idx in range(1, int(10000 // f0r) + 1):
            hm |= np.abs(f - k_idx * f0r) < 0.03 * f0r
        hnr = float(10 * np.log10((s[m_fallback & hm].sum() + 1e-30) / (s[m_fallback & ~hm].sum() + 1e-30)))

    m = (f > 80) & (f < 10000)
    centroid = float((f[m] * s[m]).sum() / (s[m].sum() + 1e-30))
    hi_share = 10 * np.log10(s[(f > 5000) & (f < 16000)].sum() / tot + 1e-30)
    odd = sum(10 ** (h / 10) for i, h in enumerate(harm[:10]) if (i + 1) % 2 == 1 and i > 0)
    even = sum(10 ** (h / 10) for i, h in enumerate(harm[:10]) if (i + 1) % 2 == 0)
    oe = 10 * np.log10((odd + 1e-30) / (even + 1e-30))
    if len(good) > 3:
        c = 1200 * np.log2(good[:, 1] / np.median(good[:, 1]))
        f0_sd = float(np.std(c))
    else:
        f0_sd = float("nan")

    # Attack: onset -> 90% of steady RMS (sustained) or peak (decaying)
    e = np.sqrt(np.convolve(x * x, np.ones(240) / 240, "same"))
    if kind == "sustained":
        target = np.median(e[int(t0 * SR):int(t1 * SR)]) if int(t1 * SR) > int(t0 * SR) else e.max()
    else:
        target = e[: int(0.3 * SR)].max()
    a10 = float(np.argmax(e > 0.1 * target) / SR)
    a90 = float(np.argmax(e > 0.9 * target) / SR)
    attack_dur_s = max(1.0 / SR, a90 - a10)
    log_attack_time = float(math.log10(attack_dur_s))

    cyclo_k1 = cyclo_coherence(seg_sustain, f0r, k=1, partials=partials, sr=SR)
    windows_serializable = {k_w: [float(v_w[0]), float(v_w[1])] for k_w, v_w in win_dict.items()}
    out = dict(
        f0=f0r,
        inharmonicity_B=float(b_fit),
        periodicity=periodicity,
        f0_sd_cents=f0_sd,
        harm_db=harm.tolist(),
        bands_db=bands.tolist(),
        hnr_db=float(hnr),
        cyclo_coherence_kf0=float(cyclo_k1),
        centroid_hz=centroid,
        centroid_rel=centroid / f0r,
        hi_share_db=float(hi_share),
        odd_even_db=float(oe),
        attack_ms=1000.0 * (a90 - a10),
        log_attack_time_log10_s=log_attack_time,
        window=[float(t0), float(t1)],
        windows=windows_serializable,
        noise_floor_db=float(nf_db),
    )
    pf, pr = partial_ratios(f, s)
    out.update(peak_hz=pf, peak_ratios=pr)
    sub = (np.abs(f - 0.5 * f0r) < 0.03 * f0r) | (np.abs(f - 1.5 * f0r) < 0.03 * f0r)
    fam = (np.abs(f - f0r) < 0.03 * f0r) | (np.abs(f - 2 * f0r) < 0.03 * f0r)
    out["subharmonic_db"] = float(10 * np.log10((s[sub].sum() + 1e-30) / (s[fam].sum() + 1e-30)))
    out.update(amp_fluct(x, f0r, t0, t1))
    if kind == "decaying":
        out.update(decay_profile(x, f0r, t1=win_dict["release"][1], partials=partials))

    if not _skip_q3_heavy:
        erb_profiles = compute_window_erb_profiles(x, win_dict, sr=SR)
        out["erb_profiles"] = erb_profiles
        int_fp = intonation_fingerprint(partials, f0r, B=b_fit)
        p_envs = partial_envelopes(partials)
        p_onsets = [float(p.onset_s) for p in partials]
        onset_tau = partial_onset_order(p_onsets)
        cyclo_multi = [
            float(cyclo_k1),
            float(cyclo_coherence(seg_sustain, f0r, k=2, partials=partials, sr=SR)),
            float(cyclo_coherence(seg_sustain, f0r, k=3, partials=partials, sr=SR)),
        ]
        out["q3"] = dict(
            partial_onset_order=dict(value=float(onset_tau), onsets_s=p_onsets, unit="kendall_tau"),
            partial_envelopes=dict(value=p_envs, unit="db_per_s_and_hz"),
            partial_freqs_and_B=dict(
                f0_hz=float(f0r),
                B=float(b_fit),
                partial_freqs_hz=[float(p.mean_f_hz) for p in partials],
                unit="Hz",
            ),
            cyclo_coherence_kf0=dict(value=float(cyclo_k1), harmonics=cyclo_multi, unit="coherence"),
            intonation_fingerprint=dict(value=int_fp, unit="cents"),
        )
        flux_db = spectral_flux_db(x, sr=SR)
        am_fm = am_fm_modulation_spectrum(partials, f0r)
        jit, shim = jitter_and_shimmer(partials, good[:, 1] if len(good) else None)
        desc_block = dict(
            log_attack_time_log10_s=dict(value=float(log_attack_time), unit="log10_s"),
            spectral_flux_db=dict(value=float(flux_db), unit="dB"),
            am_fm_modulation=am_fm,
            jitter_rel=dict(value=float(jit), unit="ratio"),
            shimmer_db=dict(value=float(shim), unit="dB"),
            band_hnr_db=dict(value=float(hnr), band_hz=[80.0, 10000.0], unit="dB"),
        )
        if band_split_hz is not None:
            desc_block["band_split"] = band_split_statistics(seg_sustain, partials, band_split_hz, sr=SR)
            out["band_split_hz"] = float(band_split_hz)
        out["descriptors"] = desc_block
    return out


def amp_fluct(x, f0, t0, t1, nh=6, n=4096, hop=480):
    """Short-time level fluctuation of the first harmonics after removing a linear trend (dB std, median
    over harmonics). Real players and turbulent sources wobble; a frozen simulated partial reads ~0."""
    win = np.hanning(n)
    ff = np.fft.rfftfreq(n, 1 / SR)
    starts = range(int(t0 * SR), max(int(t0 * SR) + 1, int(t1 * SR) - n), hop)
    rows = [np.abs(np.fft.rfft(x[s0:s0 + n] * win)) ** 2 for s0 in starts if s0 + n <= len(x)]
    if len(rows) < 8:
        return dict(amp_fluct_db=float("nan"))
    R = np.array(rows)
    tt = np.arange(len(rows), dtype=float)
    sds = []
    for k in range(1, nh + 1):
        sel = np.abs(ff - k * f0) < max(0.03 * k * f0, 12)
        if not sel.any():
            break
        lv = 10 * np.log10(R[:, sel].sum(axis=1) + 1e-30)
        resid = lv - np.polyval(np.polyfit(tt, lv, 1), tt)
        sds.append(float(np.std(resid)))
    return dict(amp_fluct_db=float(np.median(sds)) if sds else float("nan"))


def decay_profile(x, f0, t0=0.04, t1=None, nh=8, partials=None):
    """Two-slope decay and beat analysis down to the noise-floor horizon t1."""
    n, hop = 2048, 480
    dur = len(x) / SR
    t1_use = min(dur - 0.02, float(t1) if t1 is not None else dur - 0.02)
    if t1_use <= t0 + n / SR:
        t0 = 0.0
        t1_use = dur
    frames = range(int(t0 * SR), max(int(t0 * SR) + 1, int(t1_use * SR) - n), hop)
    win = np.hanning(n)
    ff = np.fft.rfftfreq(n, 1 / SR)
    rows, tt = [], []
    for s0 in frames:
        if s0 + n > len(x):
            break
        X = np.abs(np.fft.rfft(x[s0:s0 + n] * win)) ** 2
        rows.append(X)
        tt.append((s0 + n / 2) / SR)
    if len(rows) < 4:
        return dict(
            decay_db_s=float("nan"),
            slope1_db_s=float("nan"),
            slope2_db_s=float("nan"),
            knee_s=float("nan"),
            decay_residual_db=float("nan"),
            beat_rate_hz=0.0,
            beat_depth_db=0.0,
            decay_frame_hop_s=hop / SR,
            harm_decay_db_s=[],
        )
    R = np.array(rows)
    tt = np.array(tt, dtype=np.float64)
    bb = 10 * np.log10(R.sum(axis=1) + 1e-30)
    slope = float(np.polyfit(tt, bb, 1)[0])
    s1, s2, knee, resid = two_slope_decay(bb, tt)
    beat_rate, beat_depth = beats(bb, tt)
    hd = []
    for k in range(1, nh + 1):
        fk = float(partials[k - 1].mean_f_hz) if (partials and k - 1 < len(partials)) else k * f0
        sel = np.abs(ff - fk) < max(0.025 * fk, 15.0)
        if not sel.any():
            break
        lv = 10 * np.log10(R[:, sel].sum(axis=1) + 1e-30)
        ok = lv > lv.max() - 50
        hd.append(float(np.polyfit(tt[ok], lv[ok], 1)[0]) if ok.sum() > 3 else float("nan"))
    return dict(
        decay_db_s=slope,
        slope1_db_s=float(s1),
        slope2_db_s=float(s2),
        knee_s=float(knee),
        decay_residual_db=float(resid),
        beat_rate_hz=float(beat_rate),
        beat_depth_db=float(beat_depth),
        decay_frame_hop_s=float(hop / SR),
        harm_decay_db_s=hd,
    )


# ------------------------------------------------------------------------ distances
def rms_gap(a, b, floor=-60.0, gate=-50.0):
    a, b = np.asarray(a, float), np.asarray(b, float)
    n = min(len(a), len(b))
    if n == 0:
        return float("nan")
    a, b = np.maximum(a[:n], floor), np.maximum(b[:n], floor)
    keep = np.maximum(a, b) > gate
    return float(np.sqrt(np.mean((a[keep] - b[keep]) ** 2))) if keep.any() else float("nan")


def distance(da, db, nharm=20):
    env = rms_gap(da["bands_db"], db["bands_db"])
    harm = rms_gap(da["harm_db"][:nharm], db["harm_db"][:nharm])
    spec = float(math.hypot(env, harm)) if np.isfinite(harm) else float(env * math.sqrt(2))
    out = dict(
        env_db=env,
        harm_db=harm,
        spec_db=spec,
        hnr_gap=da["hnr_db"] - db["hnr_db"],
        centroid_ratio=da["centroid_rel"] / db["centroid_rel"],
        hi_share_gap=da["hi_share_db"] - db["hi_share_db"],
    )
    if "erb_profiles" in da and "erb_profiles" in db:
        att_d = erb_matrix_distance_db(da["erb_profiles"]["attack"], db["erb_profiles"]["attack"])
        sus_d = erb_matrix_distance_db(da["erb_profiles"]["sustain"], db["erb_profiles"]["sustain"])
        rel_d = erb_matrix_distance_db(da["erb_profiles"]["release"], db["erb_profiles"]["release"])
        out.update(
            attack_erb_distance_db=att_d,
            sustain_erb_distance_db=sus_d,
            release_erb_distance_db=rel_d,
        )
    if "q3" in da and "q3" in db:
        on_a = da["q3"]["partial_onset_order"]["onsets_s"]
        on_b = db["q3"]["partial_onset_order"]["onsets_s"]
        out["partial_onset_order_tau"] = partial_onset_order(on_a, on_b)
    return out


# ---------------------------------------------------------------- reference library
DYN = ("pp", "mf", "ff")


def parse_ref_name(stem: str):
    """Instrument, articulation, dynamic, string and note range from a UIowa file stem."""
    parts = stem.replace(" ", "").split(".")
    dyn = next((p for p in parts if p in DYN), None)
    rng = next((p for p in reversed(parts) if re.fullmatch(r"([A-G]b?\d){1,2}", p)), None)
    inst = parts[0].lower()
    if inst == "bass":
        inst = "doublebass"
    if inst == "bbclar":
        inst = "clarinet"
    art = "pizz" if "pizz" in parts else ("arco" if "arco" in parts else
                                          ("vib" if ("vib" in parts or "Vib" in parts) else ("novib" if any(p.lower() in ("novib", "nonvib") for p in parts) else "")))
    sul = next((p for p in parts if p.startswith("sul")), "")
    mallet = next((p for p in parts if p in ("rosewood", "hardrubber")), "")
    return dict(inst=inst, dyn=dyn, range=rng, art=art, sul=sul, mallet=mallet)


def library_entries(p: Path, heavy: bool = True):
    """Segment one UIowa file into verified notes and describe each; returns library entries."""
    p = Path(p)
    meta = parse_ref_name(p.stem)
    if not meta["range"] or not meta["dyn"]:
        return []
    exp = expected_notes(meta["range"])
    x = highpass(load(p), hz_of(exp[0]) if exp else None)
    try:
        spans = segment_notes(x, exp)
        labels = exp if len(spans) == len(exp) else None
    except ReferenceMismatch as refusal:
        print(f"{p.name}: refused, {refusal}", file=sys.stderr)
        if refusal.kind != "note-absent":
            return []
        labels = [m for m, _ in refusal.verified]
        spans = [s for _, s in refusal.verified]
    kind = "decaying" if (meta["art"] == "pizz" or meta["inst"] in ("guitar", "xylophone", "marimba", "piano")) else "sustained"
    entries = []
    for k, (a, b) in enumerate(spans):
        y = x[a:b]
        o = onset(y)
        y = y[o:]
        if len(y) < int(0.35 * SR):
            continue
        hint = hz_of(labels[k]) if labels else None
        d = describe(y, kind, hint, inharmonic=meta["inst"] in ("xylophone", "marimba"), _skip_q3_heavy=not heavy)
        if d is None:
            continue
        midi = int(round(69 + 12 * math.log2(d["f0"] / 440))) if d["f0"] and np.isfinite(d["f0"]) else None
        if hint is None and midi is not None and exp and not (exp[0] - 1 <= midi <= exp[-1] + 1):
            continue
        note = labels[k] if hint is not None else midi
        if note is None:
            continue
        entries.append(dict(file=p.name, start=(a + o) / SR, end=b / SR, note=note, name=name_of(note), kind=kind,
                            cents_vs_et=cents(d["f0"], hz_of(note)), **meta, desc=d))
    print(f"{p.name}: {len(spans)} spans, expected {len(exp)}, kept {len(entries)}", file=sys.stderr)
    return entries


def build_library(wav_dir: Path, cache: Path, jobs: int = 1, heavy: bool = True):
    """Per-file library of described reference notes; `jobs` > 1 fans the files over processes."""
    if cache.exists():
        return json.loads(cache.read_text())
    paths = sorted(Path(wav_dir).glob("*.wav"))
    if jobs > 1:
        import multiprocessing
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(jobs, mp_context=multiprocessing.get_context("fork")) as pool:
            parts = list(pool.map(library_entries, paths, [heavy] * len(paths), chunksize=1))
    else:
        parts = [library_entries(p, heavy) for p in paths]
    lib = [e for part in parts for e in part]
    cache.write_text(json.dumps(lib))
    return lib


def ref_audio(entry, wav_dir: Path, dur=None):
    f0_exp = hz_of(entry["note"]) if "note" in entry and entry["note"] is not None else None
    x = highpass(load(Path(wav_dir) / entry["file"]), f0_exp)
    a = int(entry["start"] * SR)
    b = int(entry["end"] * SR) if dur is None else a + int(dur * SR)
    return x[a:b]


def find(lib, inst, note, dyn=None, art=None, sul=None, mallet=None, tol=0):
    out = []
    for e in lib:
        if e["inst"] != inst or abs(e["note"] - note) > tol:
            continue
        if dyn and e["dyn"] != dyn:
            continue
        if art and e["art"] != art:
            continue
        if sul and e["sul"] != sul:
            continue
        if mallet and e["mallet"] != mallet:
            continue
        out.append(e)
    return out


def rulers(lib, target):
    """Real-vs-real distances around one reference note."""
    inst, note, dyn = target["inst"], target["note"], target["dyn"]
    same = [e for e in lib if e["inst"] == inst and e["note"] == note and e["art"] == target["art"] and e is not target]
    other_dyn = [distance(target["desc"], e["desc"])["spec_db"] for e in same if e["dyn"] != dyn]
    neigh = [distance(target["desc"], e["desc"])["spec_db"] for e in lib
             if e["inst"] == inst and e["art"] == target["art"] and e["dyn"] == dyn and abs(e["note"] - note) in (1, 2)]
    others = {}
    for e in lib:
        if e["inst"] != inst and abs(e["note"] - note) <= 1 and e["dyn"] == dyn and e["kind"] == target["kind"]:
            d = distance(target["desc"], e["desc"])["spec_db"]
            others[e["inst"]] = min(others.get(e["inst"], 1e9), d)
    return dict(other_dynamic=other_dyn, neighbour=neigh, other_instruments=others)


def nearest_instruments(lib, cand_desc, f0, kind, tol=2):
    """Which real instrument does the candidate sit nearest to, at nearby pitch, any dynamic?"""
    midi = 69 + 12 * math.log2(f0 / 440)
    best = {}
    for e in lib:
        if abs(e["note"] - midi) > tol + 0.5 or e["kind"] != kind:
            continue
        d = distance(cand_desc, e["desc"])["spec_db"]
        key = e["inst"] + ("." + e["art"] if e["art"] in ("pizz",) else "")
        if d < best.get(key, (1e9, None))[0]:
            best[key] = (d, f"{e['file']}@{e['name']}")
    return sorted(((k, v[0], v[1]) for k, v in best.items()), key=lambda t: t[1])
