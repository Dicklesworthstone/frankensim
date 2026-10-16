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

GRADER_VERSION = "instrgrade-1"
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


def highpass(x, f0=30.0, f1=50.0):
    n = 1 << int(np.ceil(np.log2(len(x) + SR // 2)))
    X = np.fft.rfft(x, n)
    f = np.fft.rfftfreq(n, 1 / SR)
    g = np.clip((f - f0) / (f1 - f0), 0, 1)
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


def steady_window(x, kind):
    """Analysis window relative to the (already onset-aligned) note."""
    dur = len(x) / SR
    if kind == "sustained":
        t0 = 0.30
        t1 = min(1.10, dur - 0.10)
        if t1 - t0 < 0.25:
            t0, t1 = 0.12 * dur, 0.9 * dur
    else:  # decaying: early sound carries the timbre
        t0, t1 = 0.02, min(0.52, dur)
    return t0, t1


def spectrum(x, t0, t1, nfft=1 << 18):
    seg = x[int(t0 * SR):int(t1 * SR)]
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
    return [float(q) for q in picked], [float(q / picked[0]) for q in picked]


def describe(x, kind="sustained", f0_hint=None, chord=False, inharmonic=False):
    """x is onset-aligned. Returns a dict of descriptors (all dB values relative, not absolute)."""
    t0, t1 = steady_window(x, kind)
    tr = yin_track(x[int(t0 * SR):int(t1 * SR)], fmin=30.0, fmax=2400.0)
    good = tr[(tr[:, 2] > 0.7) & np.isfinite(tr[:, 1])] if len(tr) else tr
    f0 = float(np.median(good[:, 1])) if len(good) else (f0_hint or float("nan"))
    if f0_hint and not chord:  # trust the spectrum within +-7 semitones of the target note
        fs_, ss_ = spectrum(x, t0, t1)
        f0 = f0_spectral(fs_, ss_, f0_hint * 2 ** (-7 / 12), f0_hint * 2 ** (7 / 12))
    if f0 is None or not np.isfinite(f0):
        f0 = f0_hint
    if f0 is None or not np.isfinite(f0) or f0 <= 0:
        return None
    periodicity = float(np.median(tr[:, 2])) if len(tr) else 0.0
    f, s = spectrum(x, t0, t1)
    tot = s[(f > 40) & (f < 16000)].sum() + 1e-30
    # refine f0 from harmonic peaks (least squares on k*f0)
    harm = []
    kmax = int(min(30, 15000 // f0)) if f0 and f0 > 0 else 0
    pk = []
    for k in range(1, kmax + 1):
        lo, hi = k * f0 * 0.975, k * f0 * 1.025
        sel = (f >= lo) & (f <= hi)
        if not sel.any():
            break
        j = np.argmax(s[sel])
        fk = f[sel][j]
        # integrate +-1.5% band so window leakage does not bias levels
        band = (f >= fk * 0.99) & (f <= fk * 1.01)
        harm.append(10 * np.log10(s[band].sum() / tot + 1e-30))
        pk.append((k, fk))
    harm = np.array(harm)
    if chord or inharmonic:
        harm = np.array([])
    if len(pk) >= 3 and not (chord or inharmonic):
        ks = np.array([p[0] for p in pk], float)
        fs = np.array([p[1] for p in pk], float)
        f0r = float(np.sum(ks * fs) / np.sum(ks * ks))
    else:
        f0r = f0
    # third-octave envelope (dB re total)
    # Smooth third-octave envelope: Gaussian weights in log2(f) (FWHM = 1/3 octave), so a small pitch
    # shift moves energy continuously between neighbouring bands instead of jumping at a hard edge.
    bands = []
    lf = np.log2(np.maximum(f, 1.0))
    sig = (1 / 3) / 2.3548
    for fc in THIRD_OCT:
        w = np.exp(-0.5 * ((lf - math.log2(fc)) / sig) ** 2)
        bands.append(10 * np.log10((s * w).sum() / (tot * 1.0645 * sig * 3) + 1e-30))
    bands = np.array(bands)
    # harmonic-to-noise: harmonic bins (+-3% of f0 around k f0) vs rest, 80 Hz..10 kHz
    m = (f > 80) & (f < 10000)
    hm = np.zeros_like(f, bool)
    for k in range(1, int(10000 // f0r) + 1):
        hm |= np.abs(f - k * f0r) < 0.03 * f0r
    hnr = 10 * np.log10((s[m & hm].sum() + 1e-30) / (s[m & ~hm].sum() + 1e-30))
    centroid = float((f[m] * s[m]).sum() / (s[m].sum() + 1e-30))
    hi_share = 10 * np.log10(s[(f > 5000) & (f < 16000)].sum() / tot + 1e-30)
    odd = sum(10 ** (h / 10) for i, h in enumerate(harm[:10]) if (i + 1) % 2 == 1 and i > 0)
    even = sum(10 ** (h / 10) for i, h in enumerate(harm[:10]) if (i + 1) % 2 == 0)
    oe = 10 * np.log10((odd + 1e-30) / (even + 1e-30))
    # f0 stability
    if len(good) > 3:
        c = 1200 * np.log2(good[:, 1] / np.median(good[:, 1]))
        f0_sd = float(np.std(c))
    else:
        f0_sd = float("nan")
    # attack: onset -> 90% of steady RMS (sustained) or peak (decaying)
    e = np.sqrt(np.convolve(x * x, np.ones(240) / 240, "same"))
    if kind == "sustained":
        target = np.median(e[int(t0 * SR):int(t1 * SR)])
    else:
        target = e[: int(0.3 * SR)].max()
    a10 = np.argmax(e > 0.1 * target) / SR
    a90 = np.argmax(e > 0.9 * target) / SR
    out = dict(f0=f0r, periodicity=periodicity, f0_sd_cents=f0_sd, harm_db=harm.tolist(), bands_db=bands.tolist(),
               hnr_db=float(hnr), centroid_hz=centroid, centroid_rel=centroid / f0r, hi_share_db=float(hi_share),
               odd_even_db=float(oe), attack_ms=1000 * (a90 - a10), window=[t0, t1])
    pf, pr = partial_ratios(f, s)
    out.update(peak_hz=pf, peak_ratios=pr)
    # Subharmonic (period-doubling) energy relative to the fundamental's harmonic family.
    sub = (np.abs(f - 0.5 * f0r) < 0.03 * f0r) | (np.abs(f - 1.5 * f0r) < 0.03 * f0r)
    fam = (np.abs(f - f0r) < 0.03 * f0r) | (np.abs(f - 2 * f0r) < 0.03 * f0r)
    out["subharmonic_db"] = float(10 * np.log10((s[sub].sum() + 1e-30) / (s[fam].sum() + 1e-30)))
    out.update(amp_fluct(x, f0r, t0, t1))
    if kind == "decaying":
        out.update(decay_profile(x, f0r))
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


def decay_profile(x, f0, t0=0.05, t1=1.5, nh=8):
    """dB/s slopes of broadband RMS and of the first harmonics (STFT, 4096/1024)."""
    n, hop = 4096, 1024
    t1 = min(t1, len(x) / SR - 0.1)
    frames = range(int(t0 * SR), int(t1 * SR) - n, hop)
    win = np.hanning(n)
    ff = np.fft.rfftfreq(n, 1 / SR)
    rows, tt = [], []
    for s0 in frames:
        X = np.abs(np.fft.rfft(x[s0:s0 + n] * win)) ** 2
        rows.append(X)
        tt.append((s0 + n / 2) / SR)
    if len(rows) < 4:
        return dict(decay_db_s=float("nan"), harm_decay_db_s=[])
    R = np.array(rows)
    tt = np.array(tt)
    bb = 10 * np.log10(R.sum(axis=1) + 1e-30)
    slope = float(np.polyfit(tt, bb, 1)[0])
    hd = []
    for k in range(1, nh + 1):
        sel = np.abs(ff - k * f0) < max(0.02 * k * f0, 15)
        if not sel.any():
            break
        lv = 10 * np.log10(R[:, sel].sum(axis=1) + 1e-30)
        ok = lv > lv.max() - 50
        hd.append(float(np.polyfit(tt[ok], lv[ok], 1)[0]) if ok.sum() > 3 else float("nan"))
    return dict(decay_db_s=slope, harm_decay_db_s=hd)


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
    return dict(env_db=env, harm_db=harm, spec_db=spec,
                hnr_gap=da["hnr_db"] - db["hnr_db"], centroid_ratio=da["centroid_rel"] / db["centroid_rel"],
                hi_share_gap=da["hi_share_db"] - db["hi_share_db"])


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


def build_library(wav_dir: Path, cache: Path):
    if cache.exists():
        return json.loads(cache.read_text())
    lib = []
    for p in sorted(Path(wav_dir).glob("*.wav")):
        meta = parse_ref_name(p.stem)
        if not meta["range"] or not meta["dyn"]:
            continue
        x = highpass(load(p))
        spans = segment_scale(x)
        exp = expected_notes(meta["range"])
        if len(exp) > 2 and len(spans) < len(exp) / 2:
            spans = segment_by_pitch(x, exp)
        kind = "decaying" if (meta["art"] == "pizz" or meta["inst"] in ("guitar", "xylophone", "piano")) else "sustained"
        assigned = []
        for k, (a, b) in enumerate(spans):
            y = x[a:b]
            o = onset(y)
            y = y[o:]
            if len(y) < int(0.35 * SR):
                continue
            hint = hz_of(exp[k]) if len(spans) == len(exp) else None
            d = describe(y, kind, hint, inharmonic=meta["inst"] in ("xylophone", "marimba"))
            if d is None:
                continue
            midi = int(round(69 + 12 * math.log2(d["f0"] / 440))) if d["f0"] and np.isfinite(d["f0"]) else None
            if hint is None and midi is not None and exp and not (exp[0] - 1 <= midi <= exp[-1] + 1):
                continue
            note = exp[k] if hint is not None else midi
            if note is None:
                continue
            lib.append(dict(file=p.name, start=(a + o) / SR, end=b / SR, note=note, name=name_of(note), kind=kind,
                            cents_vs_et=cents(d["f0"], hz_of(note)), **meta, desc=d))
            assigned.append(note)
        print(f"{p.name}: {len(spans)} spans, expected {len(exp)}, kept {len(assigned)}", file=sys.stderr)
    cache.write_text(json.dumps(lib))
    return lib


def ref_audio(entry, wav_dir: Path, dur=None):
    x = highpass(load(Path(wav_dir) / entry["file"]))
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
