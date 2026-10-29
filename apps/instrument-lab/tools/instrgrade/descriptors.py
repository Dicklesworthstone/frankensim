#!/usr/bin/env python3
"""M7.1-W2: Grader blind-spot fixes, Q3 bound descriptors, and §7.2 descriptors.

Numpy-only implementation (no scipy). All refusal conditions raise typed
GraderRefusalError subclasses rather than silently defaulting.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np

SR = 48000
N_ERB = 32
ERB_FMIN_HZ = 30.0
ERB_FMAX_HZ = 16000.0
ERB_RANGE_DB = 50.0  # levels this far under a window's loudest band are clamped (erb_spectrogram_db)


class GraderRefusalError(ValueError):
    """Base typed error for grader refusals."""


class AnalysisTooShort(GraderRefusalError):
    """Raised when a recording ends before decaying to the noise floor or minimum window."""


class MissingBandSplitFrequency(GraderRefusalError):
    """Raised when a band-split calculation is requested without an explicit split frequency."""


class MissingNoiseFloorEstimate(GraderRefusalError):
    """Raised when windowing requires a noise-floor estimate and none is provided or measurable."""


@dataclass
class Partial:
    k: int
    f_hz: np.ndarray
    level_db: np.ndarray
    onset_s: float
    t_s: np.ndarray
    mean_f_hz: float
    mean_level_db: float

    def to_dict(self) -> dict:
        return dict(
            k=int(self.k),
            mean_f_hz=float(self.mean_f_hz),
            mean_level_db=float(self.mean_level_db),
            onset_s=float(self.onset_s),
        )


def blackman_harris(n: int) -> np.ndarray:
    if n <= 1:
        return np.ones(max(n, 1), dtype=np.float64)
    idx = np.arange(n, dtype=np.float64)
    p = 2.0 * math.pi * idx / (n - 1)
    return 0.35875 - 0.48829 * np.cos(p) + 0.14128 * np.cos(2.0 * p) - 0.01168 * np.cos(3.0 * p)


def estimate_noise_floor_db(x: np.ndarray, sr: int = SR) -> float:
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 1 or len(x) < int(0.05 * sr) or not np.all(np.isfinite(x)):
        raise MissingNoiseFloorEstimate("IRPGRADE_NOISE_FLOOR_MISSING: invalid or too-short signal")
    hop = max(1, int(0.02 * sr))
    n_frames = len(x) // hop
    if n_frames < 2:
        raise MissingNoiseFloorEstimate("IRPGRADE_NOISE_FLOOR_MISSING: insufficient frames for floor estimate")
    frames = x[: n_frames * hop].reshape(n_frames, hop)
    rms = np.sqrt(np.mean(frames * frames, axis=1))
    floor_rms = float(np.percentile(rms, 5))
    return float(20.0 * math.log10(max(floor_rms, 1e-6)))


def _rms_50ms_db(x: np.ndarray, sr: int = SR) -> tuple[np.ndarray, np.ndarray]:
    win = max(1, int(0.05 * sr))
    if len(x) < win:
        rms = math.sqrt(float(np.mean(x * x))) if len(x) else 0.0
        return np.array([0.5 * len(x) / sr], dtype=np.float64), np.array([20.0 * math.log10(max(rms, 1e-15))], dtype=np.float64)
    csum = np.concatenate([[0.0], np.cumsum(x * x, dtype=np.float64)])
    step = max(1, int(0.005 * sr))
    starts = np.arange(0, len(x) - win + 1, step, dtype=np.int64)
    mean_sq = (csum[starts + win] - csum[starts]) / float(win)
    t_rms = (starts + 0.5 * win) / float(sr)
    rms_db = 10.0 * np.log10(np.maximum(mean_sq, 1e-30))
    return t_rms, rms_db


def windows(
    x: np.ndarray,
    kind: str,
    noise_floor_db: float | None,
    *,
    sr: int = SR,
    allow_unended_decay: bool = False,
) -> dict[str, tuple[float, float]]:
    if noise_floor_db is None or not isinstance(noise_floor_db, (int, float, np.floating)) or not math.isfinite(float(noise_floor_db)):
        raise MissingNoiseFloorEstimate("IRPGRADE_NOISE_FLOOR_MISSING: noise_floor_db must be a finite dB value")
    noise_floor_db = float(noise_floor_db)
    x = np.asarray(x, dtype=np.float64)
    dur = len(x) / float(sr)
    if dur <= 0.16:
        raise AnalysisTooShort(f"IRPGRADE_ANALYSIS_TOO_SHORT: duration {dur:.3f}s <= 0.16s")

    t_rms, rms_db = _rms_50ms_db(x, sr=sr)
    thr_db = noise_floor_db + 3.0
    i_peak = int(np.argmax(rms_db))

    if kind == "decaying":
        post_peak = rms_db[i_peak:]
        below = np.where(post_peak <= thr_db)[0]
        if len(below) == 0:
            if not allow_unended_decay and dur < 5.0:
                raise AnalysisTooShort(
                    f"IRPGRADE_ANALYSIS_TOO_SHORT: decaying note of {dur:.2f}s ended at {float(rms_db[-1]):.1f} dB "
                    f"above noise floor threshold {thr_db:.1f} dB"
                )
            t_floor = dur
        else:
            t_floor = float(t_rms[i_peak + int(below[0])])
            t_floor = min(dur, max(0.20, t_floor))
        t0 = 0.15
        t1 = min(0.55, max(0.20, 0.15 + 0.20 * (t_floor - 0.15)))
        return dict(attack=(0.0, 0.15), sustain=(t0, t1), release=(t1, t_floor))

    # Sustained note: attack (0, 0.15), sustain (0.15, t1), release (t1, t_floor)
    t0 = 0.15
    sus_mask = (t_rms >= 0.15) & (t_rms <= min(dur, max(0.45, 0.65 * dur)))
    if not sus_mask.any():
        sus_mask = t_rms >= 0.15
    plateau_db = float(np.percentile(rms_db[sus_mask], 80)) if sus_mask.any() else float(np.max(rms_db))
    # upper half only: on a short note the mask reaches into the release, which must not count as ripple
    ripple_db = plateau_db - float(np.percentile(rms_db[sus_mask], 50)) if sus_mask.any() else 0.0
    post_t0 = np.where(t_rms >= 0.20)[0]
    t1 = min(1.10, max(0.25, dur - 0.10))
    # The release begins at the last frame where the note is still on its plateau, within tol, twice the plateau's own
    # ripple (80th minus 50th percentile; 0.02 dB for a clean tone, at most 1 dB), minus 30 ms of the 50 ms RMS window.
    # It is a release only if at least five frames (25 ms) follow below it; otherwise the note sustains to the end of
    # the signal. The first 1 dB drop after 0.20 s is not the release: a note still rising to its plateau, or dipping
    # mid-note, would end its sustain there. A fixed 0.02 dB tolerance on a note with 0.1 dB ripple moves the knee by
    # tens of ms. The crossing is interpolated between frames: snapped to the 5 ms grid, it jumps a whole frame when a
    # pitch shift moves the level near the knee by a hundredth of a dB.
    knee_db = plateau_db - min(1.0, max(0.02, 2.0 * ripple_db))
    at_plateau = post_t0[rms_db[post_t0] >= knee_db]
    if len(at_plateau) and int(at_plateau[-1]) + 5 < len(rms_db):
        j = int(at_plateau[-1])
        frac = (rms_db[j] - knee_db) / (rms_db[j] - rms_db[j + 1])
        t1 = max(0.20, float(t_rms[j] + frac * (t_rms[j + 1] - t_rms[j])) - 0.030)
    # The note has ended where its level stays at the floor (+3 dB) for 0.10 s, or to the end of the signal. A
    # reverberant tail touches the floor in brief dips; ending at the first one moved a UIowa violin A3's release end
    # from 3.15 s to 3.63 s when the floor estimate moved by 2 dB.
    below = rms_db <= thr_db
    hold = 20  # 0.10 s of 5 ms frames
    ended = [k for k in np.where((t_rms >= t1) & below)[0] if below[k : k + hold].all()]
    t_floor = float(t_rms[ended[0]]) if ended else dur
    if t_floor <= t1:
        t_floor = dur
    return dict(attack=(0.0, 0.15), sustain=(t0, t1), release=(t1, t_floor))


def _erb_filterbank(freqs: np.ndarray, n_bands: int = N_ERB, fmin: float = ERB_FMIN_HZ, fmax: float = ERB_FMAX_HZ) -> tuple[np.ndarray, np.ndarray]:
    def hz_to_erb(f):
        return 21.4 * np.log10(1.0 + 0.00437 * np.asarray(f, dtype=np.float64))

    def erb_to_hz(e):
        return (10.0 ** (np.asarray(e, dtype=np.float64) / 21.4) - 1.0) / 0.00437

    e_edges = np.linspace(hz_to_erb(fmin), hz_to_erb(fmax), n_bands + 2)
    centers = erb_to_hz(e_edges)[1:-1]
    # Rounded-exponential auditory filters, W(g) = (1 + p g) exp(-p g) with g = |f - fc| / fc and p = 4 fc / ERB(fc)
    # (Glasberg & Moore 1990). Their skirts are smooth, so a partial moving by 20 cents changes a band it reaches only
    # through a skirt by about 1 dB. A finite-support triangle on FFT bins does not: a partial just outside its edge is
    # read through the analysis window's leakage, and moving it 25 Hz changed one band by 26 dB (2026-10-28).
    freqs = np.asarray(freqs, dtype=np.float64)
    erb_hz = 24.7 * (4.37 * centers / 1000.0 + 1.0)
    pg = (4.0 * centers / erb_hz)[:, None] * np.abs(freqs[None, :] - centers[:, None]) / centers[:, None]
    weights = (1.0 + pg) * np.exp(-pg)
    # Each band is the weighted mean of its bins, as before, so levels keep their scale against the -80 dB floor.
    weights /= weights.sum(axis=1, keepdims=True)
    return centers, weights


def erb_spectrogram_db(
    x: np.ndarray,
    t_start: float,
    t_end: float,
    *,
    ref_power: float | None = None,
    nfft: int = 1024,
    hop: int = 240,
    sr: int = SR,
    n_frames_out: int = 24,
) -> np.ndarray:
    a = max(0, int(round(t_start * sr)))
    b = min(len(x), int(round(t_end * sr)))
    seg = np.asarray(x[a:b], dtype=np.float64)
    if len(seg) < 64:
        return np.full((N_ERB, n_frames_out), -80.0, dtype=np.float64)
    n_use = min(nfft, max(64, 1 << int(math.floor(math.log2(len(seg))))))
    hop_use = max(16, min(hop, max(16, (len(seg) - n_use) // max(1, n_frames_out))))
    win = np.hanning(n_use)
    win_norm = float(np.sum(win * win)) + 1e-30
    starts = list(range(0, max(1, len(seg) - n_use + 1), hop_use))
    frames = []
    for s0 in starts:
        chunk = seg[s0 : s0 + n_use]
        if len(chunk) < n_use:
            chunk = np.pad(chunk, (0, n_use - len(chunk)))
        frames.append(np.abs(np.fft.rfft(chunk * win, nfft)) ** 2 / win_norm)
    power = np.stack(frames, axis=1)  # (n_freq, n_time)
    freqs = np.fft.rfftfreq(nfft, 1.0 / sr)
    _, fb = _erb_filterbank(freqs)
    erb_pow = fb @ power  # (N_ERB, n_time)
    if ref_power is None or ref_power <= 0:
        ref_power = float(np.max(erb_pow)) + 1e-30
    # Average power into n_frames_out equal time cells, each frame weighted by its overlap with the cell, so windows of
    # slightly different durations compare cleanly. Sampling the frames at n_frames_out points instead aliases: a 60 ms
    # event lands on a sample point or between two, and a 1.2 % time scale (a 20 c transposition) decides which.
    cells = np.arange(erb_pow.shape[1], dtype=np.float64)[:, None]
    edges = np.linspace(0.0, erb_pow.shape[1], n_frames_out + 1)[None, :]
    overlap = np.clip(np.minimum(cells + 1.0, edges[:, 1:]) - np.maximum(cells, edges[:, :-1]), 0.0, None)
    erb_pow = erb_pow @ (overlap / overlap.sum(axis=0, keepdims=True))
    erb_db = 10.0 * np.log10(erb_pow / ref_power + 1e-30)
    # Levels more than ERB_RANGE_DB under the window's loudest band are clamped there: such a band reads only the far
    # skirts of louder partials (masked, and a restatement of their frequency that a 20 c transposition moves by 2-5 dB
    # per band at f0 = 880 Hz) or noise. The floor is one per window: a per-frame floor follows a decaying tail's
    # loudest band, and every clamped band then inherits that band's change (a violin B3 release moved 0.85 -> 1.67 dB).
    return np.maximum(erb_db, max(float(erb_db.max()) - ERB_RANGE_DB, -80.0))


def compute_window_erb_profiles(x: np.ndarray, win_dict: dict[str, tuple[float, float]], sr: int = SR) -> dict[str, list[list[float]]]:
    s0, s1 = win_dict["sustain"]
    sus_slice = x[max(0, int(s0 * sr)) : min(len(x), int(s1 * sr))]
    ref_power = float(np.mean(sus_slice * sus_slice)) if len(sus_slice) else float(np.mean(x * x))
    ref_power = max(ref_power, 1e-24)
    att = erb_spectrogram_db(x, win_dict["attack"][0], win_dict["attack"][1], ref_power=ref_power, nfft=1024, hop=120, sr=sr, n_frames_out=20)
    sus = erb_spectrogram_db(x, win_dict["sustain"][0], win_dict["sustain"][1], ref_power=ref_power, nfft=2048, hop=480, sr=sr, n_frames_out=20)
    # Align release duration on a fixed physical horizon (0.80 s after sustain end, clipped to t_floor)
    rel_start = win_dict["release"][0]
    rel_end = min(win_dict["release"][1], rel_start + 0.80)
    if rel_end <= rel_start + 0.02:
        rel_end = win_dict["release"][1]
    rel = erb_spectrogram_db(x, rel_start, rel_end, ref_power=ref_power, nfft=1024, hop=240, sr=sr, n_frames_out=20)
    return dict(
        attack=att.tolist(),
        sustain=sus.tolist(),
        release=rel.tolist(),
    )


def erb_matrix_distance_db(a: Iterable, b: Iterable) -> float:
    aa = np.asarray(a, dtype=np.float64)
    bb = np.asarray(b, dtype=np.float64)
    if aa.shape != bb.shape or aa.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean((aa - bb) ** 2)))


def _parabolic_peak(freqs: np.ndarray, log_power: np.ndarray, idx: int) -> float:
    if 1 <= idx < len(freqs) - 1:
        y0, y1, y2 = float(log_power[idx - 1]), float(log_power[idx]), float(log_power[idx + 1])
        den = y0 - 2.0 * y1 + y2
        if abs(den) > 1e-12:
            delta = float(np.clip(0.5 * (y0 - y2) / den, -0.5, 0.5))
            df = float(freqs[1] - freqs[0])
            return float(freqs[idx]) + delta * df
    return float(freqs[idx])


def fit_stiff_string_law(ks: np.ndarray, fks: np.ndarray, weights: np.ndarray | None = None, f0_fixed: float | None = None) -> tuple[float, float]:
    ks = np.asarray(ks, dtype=np.float64)
    fks = np.asarray(fks, dtype=np.float64)
    if len(ks) == 0:
        return (float(f0_fixed or 0.0), 0.0)
    if len(ks) == 1:
        return (float(f0_fixed if f0_fixed is not None else fks[0] / ks[0]), 0.0)
    w = np.ones_like(ks) if weights is None else np.asarray(weights, dtype=np.float64)
    w = np.maximum(w, 1e-12)
    y = (fks / ks) ** 2
    x2 = ks * ks
    if f0_fixed is not None:
        f0 = float(f0_fixed)
        b_num = float(np.sum(w * x2 * (y / (f0 * f0) - 1.0)))
        b_den = float(np.sum(w * x2 * x2)) + 1e-30
        b_est = max(0.0, b_num / b_den)
        return f0, b_est

    # Weighted linear regression on y = a + b * k^2 with a = f0^2, b = f0^2 * B
    W = np.diag(w)
    A = np.column_stack([np.ones_like(x2), x2])
    coef, *_ = np.linalg.lstsq(np.sqrt(w)[:, None] * A, np.sqrt(w) * y, rcond=None)
    a, b = float(coef[0]), float(coef[1])
    if a <= 0 or b < 0:
        f0 = float(np.sum(w * ks * fks) / np.sum(w * ks * ks))
        return f0, 0.0
    f0 = math.sqrt(a)
    b_est = max(0.0, b / a)
    # 4 Gauss-Newton steps on fk = k * f0 * sqrt(1 + B * k^2)
    for _ in range(4):
        rad = np.sqrt(1.0 + b_est * x2)
        pred = ks * f0 * rad
        r = fks - pred
        j_f0 = ks * rad
        j_b = 0.5 * ks * f0 * x2 / np.maximum(rad, 1e-12)
        J = np.column_stack([j_f0, j_b])
        step, *_ = np.linalg.lstsq(np.sqrt(w)[:, None] * J, np.sqrt(w) * r, rcond=None)
        f0_new = f0 + float(step[0])
        b_new = max(0.0, b_est + float(step[1]))
        if f0_new > 0:
            f0, b_est = f0_new, b_new
    return float(f0), float(b_est)


def track_partials(
    x: np.ndarray,
    f0: float,
    kmax: int = 30,
    *,
    sr: int = SR,
    fit_inharmonicity: bool = True,
    t_window: tuple[float, float] | None = None,
    f0_fixed: float | None = None,
) -> list[Partial]:
    x = np.asarray(x, dtype=np.float64)
    if f0 is None or not math.isfinite(f0) or f0 <= 0 or len(x) < 64:
        return []
    if t_window is not None:
        a = max(0, int(round(t_window[0] * sr)))
        b = min(len(x), int(round(t_window[1] * sr)))
        seg = x[a:b] if b > a + 64 else x
    else:
        seg = x

    nfft = max(1 << 18, 1 << int(math.ceil(math.log2(len(seg) * 2))))
    win = blackman_harris(len(seg))
    win_norm = float(np.sum(win * win)) + 1e-30
    spec = (np.abs(np.fft.rfft(seg * win, nfft)) ** 2) / win_norm
    freqs = np.fft.rfftfreq(nfft, 1.0 / sr)
    log_spec = np.log(spec + 1e-30)
    tot = float(spec[(freqs > 25.0) & (freqs < 16000.0)].sum()) + 1e-30
    median_floor = float(np.median(spec[(freqs > 50.0) & (freqs < 15000.0)])) + 1e-30

    k_limit = int(min(kmax, (0.48 * sr) // f0))
    f0_curr = float(f0)
    b_curr = 0.0
    found_k: list[int] = []
    found_f: list[float] = []
    found_lvl: list[float] = []
    found_snr: list[float] = []

    for k in range(1, k_limit + 1):
        f_pred = k * f0_curr * math.sqrt(1.0 + b_curr * k * k) if fit_inharmonicity else k * f0_curr
        if f_pred >= 0.48 * sr or f_pred >= 15500.0:
            break
        tol_rel = 0.035 if len(found_k) < 4 else 0.022
        lo = max(15.0, f_pred * (1.0 - tol_rel))
        hi = min(0.49 * sr, f_pred * (1.0 + tol_rel))
        sel = np.where((freqs >= lo) & (freqs <= hi))[0]
        if len(sel) == 0:
            break
        local_idx = int(sel[int(np.argmax(spec[sel]))])
        fk = _parabolic_peak(freqs, log_spec, local_idx)
        df = float(freqs[1] - freqs[0])
        band = (freqs >= fk - max(3.5 * df, 0.008 * fk)) & (freqs <= fk + max(3.5 * df, 0.008 * fk))
        p_band = float(spec[band].sum())
        lvl_db = 10.0 * math.log10(p_band / tot + 1e-30)
        local_peak = float(spec[local_idx])
        snr = local_peak / median_floor
        found_k.append(k)
        found_f.append(fk)
        found_lvl.append(lvl_db)
        found_snr.append(snr)

        # Update stiff-string law parameters only using high-SNR partials
        good = [i for i, s_val in enumerate(found_snr) if s_val > 25.0]
        if fit_inharmonicity and len(good) >= 3:
            ks_g = np.array([found_k[i] for i in good], dtype=np.float64)
            fs_g = np.array([found_f[i] for i in good], dtype=np.float64)
            ws_g = np.array([math.log10(max(found_snr[i], 10.0)) for i in good], dtype=np.float64)
            f0_curr, b_curr = fit_stiff_string_law(ks_g, fs_g, ws_g, f0_fixed=f0_fixed)

    # Short-time tracking for f_hz(t), level_db(t), and narrowband onset_s
    n_stft = min(4096, max(256, 1 << int(math.floor(math.log2(max(256, len(x) // 4))))))
    hop_stft = max(128, min(512, n_stft // 4))
    w_stft = blackman_harris(n_stft)
    w_stft_norm = float(np.sum(w_stft * w_stft)) + 1e-30
    starts = list(range(0, max(1, len(x) - n_stft + 1), hop_stft))
    t_s = np.array([(s0 + 0.5 * n_stft) / float(sr) for s0 in starts], dtype=np.float64)
    nfft_stft = max(4096, n_stft * 2)
    f_stft = np.fft.rfftfreq(nfft_stft, 1.0 / sr)
    df_stft = float(f_stft[1] - f_stft[0])
    stft_pow = np.zeros((len(starts), len(f_stft)), dtype=np.float64)
    for i_st, s0 in enumerate(starts):
        chunk = x[s0 : s0 + n_stft]
        if len(chunk) < n_stft:
            chunk = np.pad(chunk, (0, n_stft - len(chunk)))
        stft_pow[i_st] = (np.abs(np.fft.rfft(chunk * w_stft, nfft_stft)) ** 2) / w_stft_norm
    stft_log = np.log(stft_pow + 1e-30)

    # Narrowband heterodyne onset detector on the first 0.5 s for sub-ms onset resolution
    n_onset = min(len(x), int(0.50 * sr))
    t_onset = np.arange(n_onset, dtype=np.float64) / float(sr)
    ker_len = max(32, int(round(max(0.0025, min(0.006, 1.5 / max(f0, 50.0))) * sr)))
    if ker_len % 2 == 0:
        ker_len += 1
    ker = np.hanning(ker_len)
    ker /= np.sum(ker)

    partials: list[Partial] = []
    for k, fk, mean_lvl in zip(found_k, found_f, found_lvl):
        bw = max(3.5 * df_stft, 0.015 * fk)
        sel = np.where((f_stft >= fk - bw) & (f_stft <= fk + bw))[0]
        if len(sel) == 0:
            f_traj = np.full(len(starts), fk, dtype=np.float64)
            lvl_traj = np.full(len(starts), mean_lvl, dtype=np.float64)
        else:
            f_traj = np.zeros(len(starts), dtype=np.float64)
            lvl_traj = np.zeros(len(starts), dtype=np.float64)
            for i_st in range(len(starts)):
                row = stft_pow[i_st]
                loc = int(sel[int(np.argmax(row[sel]))])
                f_traj[i_st] = _parabolic_peak(f_stft, stft_log[i_st], loc)
                lvl_traj[i_st] = 10.0 * math.log10(float(row[sel].sum()) + 1e-30)

        # Heterodyne demodulation at fk for accurate onset_s
        demod = x[:n_onset] * np.exp(-2j * math.pi * fk * t_onset)
        env_k = np.abs(np.convolve(demod, ker, mode="same"))
        pk_env = float(np.max(env_k)) if len(env_k) else 0.0
        if pk_env > 1e-12:
            idx_on = int(np.argmax(env_k >= 0.20 * pk_env))
            onset_s = float(idx_on / float(sr))
        else:
            onset_s = 0.0

        partials.append(
            Partial(
                k=int(k),
                f_hz=f_traj,
                level_db=lvl_traj,
                onset_s=onset_s,
                t_s=t_s,
                mean_f_hz=float(fk),
                mean_level_db=float(mean_lvl),
            )
        )
    return partials


def hnr_band(
    x: np.ndarray,
    partials: list[Partial] | list[tuple[int, float]],
    band: tuple[float, float] | float | None = (80.0, 10000.0),
    *,
    sr: int = SR,
) -> float:
    if band is None:
        raise MissingBandSplitFrequency("IRPGRADE_BAND_SPLIT_MISSING: band cannot be None")
    if isinstance(band, (int, float, np.floating)):
        if not math.isfinite(float(band)) or float(band) <= 0:
            raise MissingBandSplitFrequency(f"IRPGRADE_BAND_SPLIT_INVALID: {band}")
        f_lo, f_hi = float(band), 10000.0
    else:
        if len(band) != 2 or not (math.isfinite(float(band[0])) and math.isfinite(float(band[1]))) or float(band[0]) <= 0 or float(band[1]) <= float(band[0]):
            raise MissingBandSplitFrequency(f"IRPGRADE_BAND_SPLIT_INVALID: {band}")
        f_lo, f_hi = float(band[0]), float(band[1])

    x = np.asarray(x, dtype=np.float64)
    if len(x) < 64:
        return float("nan")
    win = blackman_harris(len(x))
    spec = np.abs(np.fft.rfft(x * win)) ** 2
    freqs = np.fft.rfftfreq(len(x), 1.0 / sr)
    df = float(freqs[1] - freqs[0])
    hm = np.zeros_like(freqs, dtype=bool)
    for p in partials:
        fk = float(p.mean_f_hz if isinstance(p, Partial) else p[1])
        hm |= np.abs(freqs - fk) <= 3.5 * df
    mask = (freqs >= f_lo) & (freqs <= f_hi)
    h_pow = float(spec[mask & hm].sum()) + 1e-30
    n_pow = float(spec[mask & ~hm].sum()) + 1e-30
    return float(10.0 * math.log10(h_pow / n_pow))


def two_slope_decay(level_db: np.ndarray, t: np.ndarray) -> tuple[float, float, float, float]:
    level_db = np.asarray(level_db, dtype=np.float64)
    t = np.asarray(t, dtype=np.float64)
    valid = np.isfinite(level_db) & np.isfinite(t)
    level_db, t = level_db[valid], t[valid]
    if len(t) < 8:
        slope = float(np.polyfit(t, level_db, 1)[0]) if len(t) >= 2 else float("nan")
        return slope, slope, float(t[len(t) // 2]) if len(t) else float("nan"), 0.0

    p = 10.0 ** ((level_db - np.max(level_db)) / 10.0)
    lv_Shift = level_db - np.max(level_db)
    best: tuple[float, float, float, float] | None = None
    margin = max(3, min(8, len(t) // 8))

    for i in range(margin, len(t) - margin):
        tk = float(t[i])
        # Candidate 1: Continuous piecewise-linear hinge on frame grid
        X = np.column_stack([np.ones_like(t), np.minimum(t - tk, 0.0), np.maximum(t - tk, 0.0)])
        coef, *_ = np.linalg.lstsq(X, level_db, rcond=None)
        pred_hinge = X @ coef
        res_hinge = float(np.sqrt(np.mean((level_db - pred_hinge) ** 2)))
        cand = (res_hinge, float(coef[1]), float(coef[2]), tk)

        # Candidate 2: Two-mode exponential power sum (log-sum-exp with mode equality at knee tk)
        tail = t >= (tk + 0.35 * (t[-1] - tk))
        if int(tail.sum()) >= 4:
            s2, a2 = np.polyfit(t[tail], lv_Shift[tail], 1)
            head = t <= tk
            p1 = p[head] - 10.0 ** ((a2 + s2 * t[head]) / 10.0)
            if len(p1) >= 4 and np.all(p1 > 1e-15):
                lv1 = 10.0 * np.log10(p1)
                s1, a1 = np.polyfit(t[head], lv1, 1)
                c = 0.5 * ((a1 + s1 * tk) + (a2 + s2 * tk))
                l1 = np.clip((c + s1 * (t - tk)) / 10.0, -30.0, 30.0)
                l2 = np.clip((c + s2 * (t - tk)) / 10.0, -30.0, 30.0)
                pred_lse = np.max(level_db) + 10.0 * np.log10(10.0 ** l1 + 10.0 ** l2)
                res_lse = float(np.sqrt(np.mean((level_db - pred_lse) ** 2)))
                if res_lse < cand[0]:
                    cand = (res_lse, float(s1), float(s2), tk)

        if best is None or cand[0] < best[0]:
            best = cand

    assert best is not None
    res_db, s1_out, s2_out, knee_out = best
    return float(s1_out), float(s2_out), float(knee_out), float(res_db)


def beats(level_db: np.ndarray, t: np.ndarray) -> tuple[float, float]:
    level_db = np.asarray(level_db, dtype=np.float64)
    t = np.asarray(t, dtype=np.float64)
    valid = np.isfinite(level_db) & np.isfinite(t)
    level_db, t = level_db[valid], t[valid]
    if len(t) < 8:
        return 0.0, 0.0

    dt = float(np.median(np.diff(t)))
    if dt <= 0:
        return 0.0, 0.0
    slope, intercept = np.polyfit(t, level_db, 1)
    detrended_db = level_db - (intercept + slope * t)
    p = 10.0 ** (detrended_db / 10.0)
    p_zero = p - np.mean(p)
    win = np.hanning(len(t))
    nfft = max(65536, 1 << int(math.ceil(math.log2(len(t) * 32))))
    spec = np.abs(np.fft.rfft(p_zero * win, nfft))
    freqs = np.fft.rfftfreq(nfft, dt)
    f_hi = min(25.0, 0.45 / dt)
    mask = (freqs >= 0.15) & (freqs <= f_hi)
    if not mask.any():
        return 0.0, 0.0
    sel = np.where(mask)[0]
    idx = int(sel[int(np.argmax(spec[sel]))])
    log_spec = np.log(spec + 1e-30)
    rate_hz = _parabolic_peak(freqs, log_spec, idx)

    # Exact least-squares demodulation at rate_hz
    tn = t - float(np.mean(t))
    A = np.column_stack([np.ones_like(t), tn, np.cos(2.0 * math.pi * rate_hz * t), np.sin(2.0 * math.pi * rate_hz * t)])
    coef, *_ = np.linalg.lstsq(A, p, rcond=None)
    p0 = max(float(coef[0]), 1e-30)
    amp = float(math.hypot(float(coef[2]), float(coef[3])))
    m_p = min(0.999999, amp / p0)
    depth_db = float(10.0 * math.log10((1.0 + m_p) / (1.0 - m_p)))
    return float(rate_hz), float(depth_db)


def cyclo_coherence(
    x: np.ndarray,
    f0: float,
    k: int = 1,
    *,
    partials: list[Partial] | None = None,
    sr: int = SR,
) -> float:
    x = np.asarray(x, dtype=np.float64)
    if len(x) < 256 or f0 is None or not math.isfinite(f0) or f0 <= 0 or k < 1:
        return 0.0
    period_samples = float(sr) / float(f0)
    n_periods = int(math.floor(len(x) / period_samples))
    n = int(round(n_periods * period_samples)) if n_periods >= 4 else len(x)
    n = max(256, min(len(x), n))
    seg = x[:n]
    freqs = np.fft.rfftfreq(n, 1.0 / sr)
    df = float(freqs[1] - freqs[0])
    X = np.fft.rfft(seg)
    mask = (freqs >= max(80.0, 1.5 * f0)) & (freqs <= min(12000.0, 0.45 * sr))
    if partials:
        for p in partials:
            fk = float(p.mean_f_hz)
            mask &= np.abs(freqs - fk) > max(4.0 * df, 0.025 * f0)
    for m in range(1, int(12000.0 // f0) + 1):
        mask &= np.abs(freqs - m * f0) > max(4.0 * df, 0.025 * f0)
    X_res = np.where(mask, X, 0.0)
    r = np.fft.irfft(X_res, n)
    power = r * r
    t_arr = np.arange(n, dtype=np.float64) / float(sr)
    alpha = float(k) * float(f0)
    num = float(np.abs(np.mean(power * np.exp(-2j * math.pi * alpha * t_arr))))
    den = float(np.mean(power)) + 1e-30
    return float(np.clip(num / den, 0.0, 1.0))


def partial_onset_order(onsets_a: Iterable[float], onsets_b: Iterable[float] | None = None) -> float:
    a = np.asarray(list(onsets_a), dtype=np.float64)
    b = np.arange(len(a), dtype=np.float64) if onsets_b is None else np.asarray(list(onsets_b), dtype=np.float64)
    n = min(len(a), len(b))
    if n < 2:
        return 0.0
    a, b = a[:n], b[:n]
    concordant = 0
    discordant = 0
    ties_a = 0
    ties_b = 0
    for i in range(n):
        for j in range(i + 1, n):
            da = a[j] - a[i]
            db = b[j] - b[i]
            if abs(da) <= 1e-12 and abs(db) <= 1e-12:
                ties_a += 1
                ties_b += 1
            elif abs(da) <= 1e-12:
                ties_a += 1
            elif abs(db) <= 1e-12:
                ties_b += 1
            elif da * db > 0:
                concordant += 1
            else:
                discordant += 1
    n0 = n * (n - 1) // 2
    den = math.sqrt(max(1, (n0 - ties_a) * (n0 - ties_b)))
    return float((concordant - discordant) / den)


def intonation_fingerprint(partials: list[Partial], f0: float, B: float = 0.0) -> list[float]:
    if f0 is None or not math.isfinite(f0) or f0 <= 0:
        return []
    out: list[float] = []
    for p in partials:
        f_ref = p.k * f0 * math.sqrt(max(1.0, 1.0 + B * p.k * p.k))
        if p.mean_f_hz > 0 and f_ref > 0:
            out.append(float(1200.0 * math.log2(p.mean_f_hz / f_ref)))
        else:
            out.append(float("nan"))
    return out


def partial_envelopes(partials: list[Partial]) -> list[dict]:
    envs: list[dict] = []
    for p in partials[:12]:
        if len(p.t_s) >= 8:
            s1, s2, knee, res = two_slope_decay(p.level_db, p.t_s)
            br, bd = beats(p.level_db, p.t_s)
        else:
            s1 = s2 = knee = res = float("nan")
            br = bd = 0.0
        envs.append(
            dict(
                k=int(p.k),
                onset_s=float(p.onset_s),
                slope1_db_s=float(s1),
                slope2_db_s=float(s2),
                knee_s=float(knee),
                residual_db=float(res),
                beat_rate_hz=float(br),
                beat_depth_db=float(bd),
            )
        )
    return envs


def spectral_flux_db(x: np.ndarray, sr: int = SR, nfft: int = 2048, hop: int = 480) -> float:
    x = np.asarray(x, dtype=np.float64)
    if len(x) < nfft + hop:
        return 0.0
    win = np.hanning(nfft)
    starts = range(0, len(x) - nfft + 1, hop)
    frames = [10.0 * np.log10(np.abs(np.fft.rfft(x[s0 : s0 + nfft] * win)) ** 2 + 1e-20) for s0 in starts]
    if len(frames) < 2:
        return 0.0
    mat = np.stack(frames, axis=0)
    mat = np.maximum(mat, np.max(mat) - 80.0)
    diffs = np.diff(mat, axis=0)
    return float(np.mean(np.sqrt(np.mean(diffs * diffs, axis=1))))


def am_fm_modulation_spectrum(partials: list[Partial], f0: float) -> dict:
    if not partials or len(partials[0].t_s) < 8:
        return dict(
            am_rate_hz=dict(value=0.0, unit="Hz"),
            am_depth_db=dict(value=0.0, unit="dB"),
            fm_rate_hz=dict(value=0.0, unit="Hz"),
            fm_depth_cents=dict(value=0.0, unit="cents"),
        )
    p0 = partials[0]
    am_rate, am_depth = beats(p0.level_db, p0.t_s)
    t = p0.t_s
    dt = float(np.median(np.diff(t)))
    f_traj = np.maximum(p0.f_hz, 1e-6)
    cents_traj = 1200.0 * np.log2(f_traj / float(np.median(f_traj)))
    cents_zero = cents_traj - np.mean(cents_traj)
    if dt > 0 and len(t) >= 8:
        w = np.hanning(len(t))
        nfft = max(16384, 1 << int(math.ceil(math.log2(len(t) * 16))))
        spec = np.abs(np.fft.rfft(cents_zero * w, nfft))
        freqs = np.fft.rfftfreq(nfft, dt)
        mask = (freqs >= 0.2) & (freqs <= min(25.0, 0.45 / dt))
        if mask.any():
            sel = np.where(mask)[0]
            idx = int(sel[int(np.argmax(spec[sel]))])
            fm_rate = _parabolic_peak(freqs, np.log(spec + 1e-30), idx)
            A = np.column_stack([np.ones_like(t), np.cos(2.0 * math.pi * fm_rate * t), np.sin(2.0 * math.pi * fm_rate * t)])
            coef, *_ = np.linalg.lstsq(A, cents_zero, rcond=None)
            fm_depth = float(math.hypot(float(coef[1]), float(coef[2])))
        else:
            fm_rate, fm_depth = 0.0, 0.0
    else:
        fm_rate, fm_depth = 0.0, 0.0
    return dict(
        am_rate_hz=dict(value=float(am_rate), unit="Hz"),
        am_depth_db=dict(value=float(am_depth), unit="dB"),
        fm_rate_hz=dict(value=float(fm_rate), unit="Hz"),
        fm_depth_cents=dict(value=float(fm_depth), unit="cents"),
    )


def jitter_and_shimmer(partials: list[Partial], f0_track: np.ndarray | None = None) -> tuple[float, float]:
    jitter = float("nan")
    if f0_track is not None and len(f0_track) >= 4:
        good = f0_track[np.isfinite(f0_track) & (f0_track > 0)]
        if len(good) >= 4:
            periods = 1.0 / good
            jitter = float(np.mean(np.abs(np.diff(periods))) / np.mean(periods))
    elif partials and len(partials[0].f_hz) >= 4:
        good = partials[0].f_hz[np.isfinite(partials[0].f_hz) & (partials[0].f_hz > 0)]
        if len(good) >= 4:
            periods = 1.0 / good
            jitter = float(np.mean(np.abs(np.diff(periods))) / np.mean(periods))

    shimmer = float("nan")
    if partials and len(partials[0].level_db) >= 4:
        lv = partials[0].level_db[np.isfinite(partials[0].level_db)]
        if len(lv) >= 4:
            shimmer = float(np.mean(np.abs(np.diff(lv))))
    return jitter, shimmer


def band_split_statistics(
    x: np.ndarray,
    partials: list[Partial],
    band_split_hz: float | None,
    *,
    sr: int = SR,
) -> dict:
    if band_split_hz is None or not isinstance(band_split_hz, (int, float, np.floating)) or not math.isfinite(float(band_split_hz)) or float(band_split_hz) <= 0:
        raise MissingBandSplitFrequency("IRPGRADE_BAND_SPLIT_MISSING: band_split_hz must be a positive finite frequency in Hz")
    f_split = float(band_split_hz)
    x = np.asarray(x, dtype=np.float64)
    win = blackman_harris(len(x))
    spec = np.abs(np.fft.rfft(x * win)) ** 2
    freqs = np.fft.rfftfreq(len(x), 1.0 / sr)
    tot = float(spec[(freqs >= 30.0) & (freqs <= 16000.0)].sum()) + 1e-30
    above = (freqs >= f_split) & (freqs <= 16000.0)
    above_pow = float(spec[above].sum()) + 1e-30
    share_above_db = float(10.0 * math.log10(above_pow / tot))
    centroid_above_hz = float((freqs[above] * spec[above]).sum() / above_pow) if above.any() else float("nan")
    hnr_above_db = hnr_band(x, partials, (f_split, min(16000.0, max(f_split + 500.0, 10000.0))), sr=sr)
    partials_above = [p for p in partials if p.mean_f_hz >= f_split]
    return dict(
        band_split_hz=dict(value=f_split, unit="Hz"),
        hnr_above_db=dict(value=float(hnr_above_db), unit="dB"),
        energy_share_above_db=dict(value=float(share_above_db), unit="dB"),
        centroid_above_hz=dict(value=float(centroid_above_hz), unit="Hz"),
        partial_count_above=dict(value=int(len(partials_above)), unit="count"),
    )


def null_baselines(
    target_desc: dict,
    *,
    cand_spec_db: float | None = None,
    dyn_ruler_db: float | None = None,
    sr: int = SR,
) -> dict:
    import instrgrade as ig

    f0 = float(target_desc["f0"])
    t = np.arange(int(1.6 * sr), dtype=np.float64) / float(sr)
    kmax = int(min(30, 15000 // max(f0, 20.0)))
    waveforms = {
        "sine": np.sin(2.0 * math.pi * f0 * t),
        "1/n": sum(np.sin(2.0 * math.pi * k * f0 * t) / k for k in range(1, kmax + 1)),
        "odd_1/n": sum(np.sin(2.0 * math.pi * k * f0 * t) / k for k in range(1, kmax + 1, 2)),
        "1/n2": sum(np.sin(2.0 * math.pi * k * f0 * t) / (k * k) for k in range(1, kmax + 1)),
    }
    distances_db: dict[str, float] = {}
    for name, wav in waveforms.items():
        d_null = ig.describe(wav, kind="sustained", f0_hint=f0, noise_floor_db=-90.0, _skip_q3_heavy=True)
        distances_db[name] = float(ig.distance(d_null, target_desc)["spec_db"])
    ordering = sorted(distances_db.keys(), key=lambda k: distances_db[k])
    margins_vs_cand = {
        name: (float(dist - cand_spec_db) if cand_spec_db is not None and math.isfinite(cand_spec_db) else None)
        for name, dist in distances_db.items()
    }
    margins_vs_dyn = {
        name: (float(dist - dyn_ruler_db) if dyn_ruler_db is not None and math.isfinite(dyn_ruler_db) else None)
        for name, dist in distances_db.items()
    }
    return dict(
        unit="dB",
        distances_db=distances_db,
        ordering=ordering,
        margins_vs_candidate_db=margins_vs_cand,
        margins_vs_dynamic_ruler_db=margins_vs_dyn,
    )
