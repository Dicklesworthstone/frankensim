"""Where a spectral-flux gap between a render and a recording lives.

The graded flux (decompose.flux_curve) is the RMS over live STFT cells of the positive frame-to-frame dB step, on a 4096-point Hann
window with a 10 ms hop. Nothing here models the voice: every number is that same flux, recomputed on the REAL signal and on the
render with one thing changed, so the attribution is measured.

  cells      the flux restricted to the cells at the harmonics of the frame's f0 and to the cells between them
  f0         the pitch movement (jitter, vibrato, drift) taken out of both signals by time-warping to a constant f0, and the
             real signal's movement put into the render the same way
  level      the same for the level movement (shimmer, tremor): the level curve's departure from its running median
  roundtrip  flatten then impose on one signal: the resampling's own effect on the flux (it must be small)

and the measured statistics the flux responds to: f0 jitter per frame, vibrato rate and extent, level shimmer, and the ratio of
power at the harmonics to power between them, per band.

A share is (gap - gap after) / gap, with the gap the graded one: the RMS of the frame-by-frame flux difference over the frames
both signals sound in. The first class (cells) is a restriction of the metric, the others are ablations of the signal.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import decompose as dc  # noqa: E402
from decompose import HOP, SR, flux_cells, stft_db  # noqa: E402

N_FFT = 4096  # decompose.Note's window for every note above 80 Hz: the flux is a property of this window
BANDS = (40.0, 1000.0, 2000.0, 4000.0, 8000.0, 16000.0)
VIB_BAND = (3.5, 9.0)  # Hz: where singing vibrato lives
LEVEL_BAND = (80.0, 5000.0)
MEDIAN_FRAMES = 101  # the level curve's running median: movement faster than ~1 Hz is shimmer, slower is the note's own shape
SHIMMER_CLIP_DB = 6.0
EDGE_FRAMES = 20  # 0.2 s: the onset and release of a note are its shape, not shimmer


class Look:
    """One signal seen the way the grader sees it: the STFT map, its sounding frames, its f0 track and its flux cells."""

    def __init__(self, y, f0_hz):
        S, self.sf = stft_db(y, N_FFT)
        self.frames = max(1, (len(y) - N_FFT // 2) // HOP)
        self.S = S[:, : self.frames].astype(np.float64)
        self.rms_db = 10 * np.log10(np.sum(10 ** (self.S / 10), axis=0) + 1e-20)
        self.sound = self.rms_db > self.rms_db.max() - 40
        self.flux = dc.flux_curve(self.S, self.sf)
        self.pos, self.live, self.rows = flux_cells(self.S, self.sf)
        comb = f0_track(self.S, self.sf, f0_hz, self.sound)
        f0_med = float(np.nanmedian(comb)) if np.isfinite(comb).any() else f0_hz
        self.f0_inst = f0_instant(y, f0_med)
        self.f0_t = np.where(self.sound, frame_mean(self.f0_inst, self.frames), np.nan)
        self.sound_samples = np.interp(np.arange(len(y)), np.arange(self.frames) * HOP, self.sound.astype(float)) > 0.5
        lv = (self.sf >= LEVEL_BAND[0]) & (self.sf <= LEVEL_BAND[1])
        self.level_db = 10 * np.log10(np.sum(10 ** (self.S[lv] / 10), axis=0) + 1e-20)


def f0_track(S, sf, f0_hz, sound):
    """f0 [Hz] per frame, NaN where the frame is not sounding: the frequency whose harmonics carry the most amplitude, found on a
    coarse grid (+-100 cents around the note) and then a fine one (+-30 cents around the median), with a parabolic peak."""
    A = 10 ** (S / 20)
    df = float(sf[1] - sf[0])
    ks = np.arange(1, int(min(24, max(8, 4000.0 // f0_hz))) + 1)

    def score(freqs):
        pos = (freqs[:, None] * ks[None, :] - sf[0]) / df
        i0 = np.clip(np.floor(pos).astype(int), 0, A.shape[0] - 2)
        w = np.clip(pos - i0, 0.0, 1.0)[..., None]
        return (A[i0] * (1 - w) + A[i0 + 1] * w).sum(axis=1)

    coarse = f0_hz * 2 ** (np.arange(-100.0, 100.01, 2.0) / 1200)
    f_med = float(np.median(coarse[np.argmax(score(coarse), axis=0)][sound])) if sound.any() else f0_hz
    cents = np.arange(-30.0, 30.01, 0.5)
    sc = score(f_med * 2 ** (cents / 1200))
    idx = np.argmax(sc, axis=0)
    cols = np.arange(sc.shape[1])
    lo, hi = np.clip(idx - 1, 0, len(cents) - 1), np.clip(idx + 1, 0, len(cents) - 1)
    a, b, c = sc[lo, cols], sc[idx, cols], sc[hi, cols]
    den = a - 2 * b + c
    delta = np.where((np.abs(den) > 1e-30) & (idx > 0) & (idx < len(cents) - 1), 0.5 * (a - c) / np.where(den == 0, 1, den), 0.0)
    out = f_med * 2 ** ((cents[idx] + np.clip(delta, -1, 1) * 0.5) / 1200)
    return np.where(sound, out, np.nan)


def f0_instant(y, f0_med, kmax=6, smooth_ms=3.0):
    """f0 [Hz] at every sample: the phase derivative of each of the first harmonics (a raised-cosine band of +-0.45 f0 about the
    harmonic, so a pitch movement of several tens of cents stays inside it), divided by its number, weighted by k^2 times the band's
    power (the phase noise of a harmonic falls as 1/k and as its amplitude rises), then smoothed over `smooth_ms`."""
    n = len(y)
    Y = np.fft.fft(y)
    fr = np.fft.fftfreq(n, 1 / SR)
    win = max(3, int(smooth_ms * 1e-3 * SR))
    box = np.ones(win) / win
    num, den = np.zeros(n), np.zeros(n)
    for k in range(1, min(kmax, max(1, int(3000 // f0_med))) + 1):
        d = np.abs(fr - k * f0_med)
        half = 0.45 * f0_med
        band = np.where((d < half) & (fr > 0), 0.5 * (1 + np.cos(np.pi * d / half)), 0.0)
        a = np.fft.ifft(2 * Y * band)
        inst = np.gradient(np.unwrap(np.angle(a))) * SR / (2 * np.pi) / k
        w = k * k * np.convolve(np.abs(a) ** 2, box, mode="same")
        num += w * inst
        den += w
    f = num / np.maximum(den, 1e-300)
    return np.where(den > 0, np.convolve(f, box, mode="same"), f0_med)


def frame_mean(x, frames):
    """x averaged with the analysis window (Hann, N_FFT) at the centre of each STFT frame: what the flux 'sees' of a sample-rate curve."""
    w = np.hanning(N_FFT)
    w /= w.sum()
    sm = np.convolve(x, w, mode="same") / np.convolve(np.ones(len(x)), w, mode="same")  # normalised, so the ends are not pulled to zero
    return sm[np.minimum(np.arange(frames) * HOP, len(x) - 1)]


COMB_F0_SECONDS = 0.3  # the comb follows the pitch smoothed over this long: its reach (3 bins + 0.5 %) is wide enough for the movement left over,
# and a track that follows the estimate's few-cent wobble moves the comb's edge bins from frame to frame


def comb_track(f0_inst):
    """f0 [Hz] per sample as the comb follows it: `f0_inst` smoothed with a Hann window of COMB_F0_SECONDS (normalised at the ends)."""
    n = len(f0_inst)
    hw = np.hanning(int(COMB_F0_SECONDS * SR) | 1)
    hw /= hw.sum()
    return np.convolve(f0_inst, hw, mode="same") / np.convolve(np.ones(n), hw, mode="same")


def comb_weight(fr, f0, df):
    """Weight (0..1) of each bin of frequency `fr` in the comb about the multiples of `f0`: 1 within 3 bins plus 0.5 % of its frequency
    of one, a raised-cosine fall over the next two bins, 0 beyond. The edge is soft because a hard one flips bins in and out as the
    f0 estimate wobbles by a few cents, and a gain on those bins then steps the skirt cells from frame to frame."""
    h = np.rint(fr / f0)
    over = np.abs(fr - h * f0) - (3 * df + 0.005 * fr)
    return np.where(h >= 1, np.where(over <= 0, 1.0, np.where(over < 2 * df, 0.5 * (1 + np.cos(np.pi * over / (2 * df))), 0.0)), 0.0)


def _comb_filter(y, f0_inst, gain):
    """Filter `y` with a gain that follows its pitch: Hann-windowed frames at 75 % overlap (the windows sum to 2), each zero-padded to
    twice its length before the FFT so that the gain's impulse response does not wrap around the frame, `gain(fr, f0, c)` giving the
    per-bin gain at the frame centred on sample c whose comb track reads f0."""
    n, N = len(y), N_FFT
    H, pad = N // 4, N  # with a pad of N, four windows cover every sample of y, including the first and the last
    yp = np.concatenate([np.zeros(pad), y, np.zeros(pad + 2 * N)])
    w = np.hanning(N + 1)[:-1]  # periodic: the windows then sum to exactly 2 at this hop
    fr = np.fft.rfftfreq(2 * N, 1 / SR)
    f0_win = comb_track(f0_inst)
    out = np.zeros(len(yp) + N)
    for s in range(0, len(yp) - N + 1, H):
        c = int(np.clip(s + N // 2 - pad, 0, n - 1))
        out[s: s + 2 * N] += np.fft.irfft(np.fft.rfft(yp[s: s + N] * w, 2 * N) * gain(fr, f0_win[c], c), 2 * N)
    return out[pad: pad + n] / 2.0


def split_harmonic(y, f0_inst):
    """y = harmonic + residual, exactly: the part of y inside a comb about the harmonics of its own f0 track (`comb_weight`, the comb
    following the pitch smoothed over COMB_F0_SECONDS), and the rest. Returns (harmonic, residual). The residual is the breath, the room and
    whatever else is not at a multiple of the pitch; the harmonic part keeps the noise that falls inside the comb. The comb excludes DC,
    so a signal cut off abruptly has the mean of its last window in the residual."""
    harm = _comb_filter(y, f0_inst, lambda fr, f0, c: comb_weight(fr, f0, SR / N_FFT))
    return harm, y - harm


def harmonic_cells(look):
    """(cells, steps) True where the cell is at a harmonic of the frame's f0: within 3 bins plus 0.5 % of its own frequency of a
    multiple, which is the reach of a partial under this window and a pitch that moves a few cents. The frame is the later one of
    the step."""
    f = look.sf[look.rows][:, None]
    with np.errstate(invalid="ignore"):
        f0 = look.f0_t[None, 1:]
        h = np.rint(f / f0)
        return (h >= 1) & (np.abs(f - h * f0) <= 3 * float(look.sf[1] - look.sf[0]) + 0.005 * f)


def class_flux(look, cls):
    """The flux curve of the cells in `cls` only (RMS over the live cells of that class)."""
    num = np.sum(np.where(cls, look.pos, 0.0) ** 2, axis=0)
    return np.concatenate([[0.0], np.sqrt(num / np.maximum((look.live & cls).sum(axis=0), 1))])


def flux_gap(fa, sa, fb, sb, within=None):
    """Mean flux of each signal over its sounding frames, and the graded gap: the RMS of the frame-by-frame difference over the
    frames both sound in (decompose's definition). `within` (a frame mask) restricts all three to those frames."""
    nf = min(len(fa), len(fb))
    if within is not None:
        w = np.zeros(max(len(fa), len(fb)), bool)
        w[: len(within)] = within[: len(w)]
        sa, sb = sa & w[: len(sa)], sb & w[: len(sb)]
    both = sa[:nf] & sb[:nf]
    d = fa[:nf][both] - fb[:nf][both]
    return dict(real=float(np.mean(fa[sa])), ours=float(np.mean(fb[sb])), gap=float(np.sqrt(np.mean(d * d))), bias=float(np.mean(d)))


def steady_frames(sound):
    """The sounding frames with EDGE_FRAMES shaved off each end of the sounding span: the note's body, without its onset and release."""
    snd = np.flatnonzero(sound)
    m = np.zeros(len(sound), bool)
    if len(snd) > 2 * EDGE_FRAMES:
        m[snd[0] + EDGE_FRAMES: snd[-1] - EDGE_FRAMES + 1] = True
    return m & sound


def _fractional(y, t):
    """y at fractional sample positions t: FFT-upsampled four times, then linear (a plain linear interpolation of y itself would
    lose the top octave by a position-dependent amount, which is amplitude modulation there and would show up as flux)."""
    n, up = len(y), 4
    Y = np.fft.rfft(y)
    Y4 = np.zeros(up * n // 2 + 1, dtype=complex)
    Y4[: len(Y)] = Y[: len(Y4)]
    y4 = np.fft.irfft(Y4, n=up * n) * up
    return np.interp(np.asarray(t) * up, np.arange(up * n), y4)


def retime(y, ratio_s, mode):
    """Time-warp `y` by a pitch-ratio track (f0 / its median, one value per sample, 1 where the signal is not sounding). 'flatten'
    removes that movement from y; 'impose' gives y that movement (a track shorter than y holds its last value)."""
    n = len(y)
    r = np.interp(np.arange(n), np.arange(len(ratio_s)), ratio_s)
    theta = np.concatenate([[0.0], np.cumsum(r[:-1])])  # samples of reference time elapsed at each output sample
    s = np.arange(n, dtype=float)
    t = np.interp(s, theta, s) if mode == "flatten" else theta
    return _fractional(y, np.clip(t, 0, n - 1))


def ratio_track(look):
    """f0 / its median over the sounding samples, per sample, 1 outside them."""
    s = look.sound_samples
    if not s.any():
        return np.ones(len(look.f0_inst))
    return np.where(s, look.f0_inst / np.median(look.f0_inst[s]), 1.0)


def _departure(level_db, sound):
    """A level curve's departure from its smoothed running median [dB, clipped], zero outside the sounding frames."""
    k = MEDIAN_FRAMES
    pad = np.pad(level_db, (k // 2, k // 2), mode="edge")
    med = np.median(np.lib.stride_tricks.sliding_window_view(pad, k), axis=1)
    # the running median of a swinging level steps; the baseline is its moving average, so that removing the departure puts no step into the gain
    base = np.convolve(np.pad(med, (k // 2, k // 2), mode="edge"), np.ones(k) / k, mode="valid")
    d = np.clip(level_db - base, -SHIMMER_CLIP_DB, SHIMMER_CLIP_DB)
    # the note's own onset and release are its shape, not shimmer: the curve is zero outside the sounding frames and fades in and out
    # over EDGE_FRAMES frames inside them
    snd = np.flatnonzero(sound)
    if not len(snd):
        return np.zeros(len(d))
    idx = np.arange(len(d))
    edge = np.clip(np.minimum(idx - snd[0], snd[-1] - idx) / EDGE_FRAMES - 1.0, 0.0, 1.0)
    return np.where(sound, d * 0.5 * (1 - np.cos(np.pi * edge)), 0.0)


def shimmer_curve(look):
    """The broadband level's departure from its smoothed running median [dB]."""
    return _departure(look.level_db, look.sound)


def partial_levels(look, kmax=12):
    """Each of the first `kmax` partials' own level per frame [dB]: the power within +-2 bins of k f0 (NaN f0 frames read the floor)."""
    df = float(look.sf[1] - look.sf[0])
    P = 10 ** (look.S / 10)
    lv = []
    for k in range(1, kmax + 1):
        with np.errstate(invalid="ignore"):
            near = np.abs(look.sf[:, None] - k * look.f0_t[None, :]) <= 2 * df
        lv.append(10 * np.log10((P * near).sum(axis=0) + 1e-20))
    return np.array(lv)


def partial_deltas(look, kmax=12):
    """(kmax, frames) departure of each partial's own level from its smoothed running median [dB]; zero for a partial more than 35 dB
    under the strongest, whose 'level' is the noise under it."""
    lv = partial_levels(look, kmax)
    body = steady_frames(look.sound) & np.isfinite(look.f0_t)
    mean_db = np.array([10 * np.log10(np.mean(10 ** (l[body] / 10)) + 1e-20) if body.any() else -200 for l in lv])
    d = np.array([_departure(l, look.sound) for l in lv])
    d[mean_db < mean_db.max() - 35] = 0.0
    return d


def modulate_partials(y, f0_inst, deltas_db, sign):
    """y with every harmonic's level moved by sign * deltas_db[k - 1, frame] dB (the harmonic comb of `split_harmonic`, a gain on each
    partial's bins that fades out over the comb's edge; a harmonic above the last row moves like the last): sign -1 takes a partial level
    movement out, +1 puts one in."""
    K, last = deltas_db.shape

    def gain(fr, f0, c):
        h = np.rint(fr / f0)
        frame = min(int(round(c / HOP)), last - 1)
        return 10 ** (comb_weight(fr, f0, SR / N_FFT) * sign * deltas_db[np.clip(h.astype(int), 1, K) - 1, frame] / 20)

    return _comb_filter(y, f0_inst, gain)


def regain(y, d_db, sign):
    g = np.interp(np.arange(len(y)), np.arange(len(d_db)) * HOP, d_db)
    return y * 10 ** (sign * g / 20)


def _seen(y, f0_hz):
    S, sf = stft_db(y, N_FFT)
    frames = max(1, (len(y) - N_FFT // 2) // HOP)
    S = S[:, :frames].astype(np.float64)
    rms = 10 * np.log10(np.sum(10 ** (S / 10), axis=0) + 1e-20)
    return dc.flux_curve(S, sf), rms > rms.max() - 40


def track_stats(look):
    """Measured per-signal statistics: f0 jitter (RMS frame-to-frame step, cents), the slow f0 movement and its vibrato-band part
    (rate Hz, peak-to-peak extent cents), the level shimmer (RMS frame-to-frame step and RMS departure from the running median, dB),
    and per band the median over sounding frames of 10 log10 (power at the harmonics / power between them)."""
    ok = np.flatnonzero(np.isfinite(look.f0_t))
    out = dict(sounding_frames=int(look.sound.sum()))
    if len(ok) < 8:
        return out
    c = 1200 * np.log2(look.f0_t / np.median(look.f0_t[ok]))
    run = ok[0], ok[-1]
    cc = np.interp(np.arange(run[0], run[1] + 1), ok, c[ok])
    consecutive = np.diff(ok) == 1
    step = np.diff(c[ok])[consecutive]
    out["f0_jitter_cents_rms"] = float(np.sqrt(np.mean(step ** 2))) if len(step) else None
    cc = cc - np.polyval(np.polyfit(np.arange(len(cc)), cc, 1), np.arange(len(cc)))
    out["f0_movement_cents_rms"] = float(np.std(cc))
    w = np.hanning(len(cc))
    X = np.abs(np.fft.rfft(cc * w))
    fr = np.fft.rfftfreq(len(cc), HOP / SR)
    band = (fr >= VIB_BAND[0]) & (fr <= VIB_BAND[1])
    if band.any() and len(cc) >= 32:
        k = np.argmax(np.where(band, X, 0))
        out["vibrato_rate_hz"] = float(fr[k])
        out["vibrato_extent_cents_pp"] = float(2 * 2 * X[k] / w.sum())
        out["vibrato_band_cents_rms"] = float(np.sqrt(np.sum((2 * X[band] / w.sum()) ** 2) / 2))
    body = steady_frames(look.sound) & np.isfinite(look.f0_t)
    ls = np.diff(look.level_db)[body[1:] & body[:-1]]
    out["level_step_db_rms"] = float(np.sqrt(np.mean(ls ** 2))) if len(ls) else None
    out["level_shimmer_db_rms"] = float(np.sqrt(np.mean(shimmer_curve(look)[look.sound] ** 2)))
    P = 10 ** (look.S / 10)
    f = look.sf[:, None]
    with np.errstate(invalid="ignore"):
        h = np.rint(f / look.f0_t[None, :])
        hm = (h >= 1) & (np.abs(f - h * look.f0_t[None, :]) <= 3 * float(look.sf[1] - look.sf[0]) + 0.005 * f)
    bands = {}
    for lo, hi in zip(BANDS[:-1], BANDS[1:]):
        r = (look.sf >= lo) & (look.sf < hi)
        ph, pb = (P[r] * hm[r]).sum(axis=0), (P[r] * ~hm[r]).sum(axis=0)
        v = (body & (pb > 0) & (ph > 0))
        ratio = 10 * np.log10(ph[v] / pb[v]) if v.any() else np.array([np.nan])
        bands[f"{lo:.0f}-{hi:.0f}"] = dict(harmonic_over_between_db_median=float(np.median(ratio)), over_time_db_sd=float(np.std(ratio)))
    out["bands"] = bands
    return out


def cell_shares(look):
    """Over the steps in the body of the note (`steady_frames`): the share of the live cells and of the summed squared step that sits at
    harmonics and between them, in total and per band."""
    hc = harmonic_cells(look)
    snd = steady_frames(look.sound)[1:]
    sq = look.pos ** 2
    fr = look.sf[look.rows]
    res = {}
    for name, sel in (("all", np.ones(len(fr), bool)),) + tuple((f"{lo:.0f}-{hi:.0f}", (fr >= lo) & (fr < hi)) for lo, hi in zip(BANDS[:-1], BANDS[1:])):
        tot_sq = sq[sel][:, snd].sum()
        tot_live = look.live[sel][:, snd].sum()
        row = {}
        for cname, cls in (("harmonic", hc), ("between", ~hc)):
            row[cname] = dict(live_cell_share=float((look.live & cls)[sel][:, snd].sum() / max(tot_live, 1)),
                              squared_step_share=float((sq * cls)[sel][:, snd].sum() / max(tot_sq, 1e-30)))
        row["live_cells_per_step"] = float(tot_live / max(snd.sum(), 1))
        res[name] = row
    return res


def level_shares(look):
    """The flux cells by how far below the map's peak the louder cell of the step is: partial peaks (0-20 dB), their skirts and
    strong noise (20-40), and cells near the live floor (40-60). Per class, over the steps in the body of the note (`steady_frames`): live
    cells per step, the share of the summed squared step, and the RMS positive step in that class (0 where the cell falls)."""
    X = look.S[look.rows]
    below = X.max() - np.maximum(X[:, 1:], X[:, :-1])
    snd = steady_frames(look.sound)[1:]
    sq = (look.pos ** 2)[:, snd]
    out = {}
    for lo, hi in ((0, 20), (20, 40), (40, 60)):
        cls = (look.live & (below >= lo) & (below < hi))[:, snd]
        n = int(cls.sum())
        out[f"{lo}-{hi}"] = dict(live_cells_per_step=float(n / max(snd.sum(), 1)), squared_step_share=float((sq * cls).sum() / max(sq.sum(), 1e-30)),
                                 positive_step_rms_db=float(np.sqrt((sq * cls).sum() / max(n, 1))))
    return out


def partial_flutter(look, kmax=8):
    """Frame-to-frame flutter of each partial's own level [dB, RMS step over consecutive sounding frames]: the power within +-2 bins of
    k f0 (so a pitch movement of a few cents does not read as level), for the partials within 35 dB of the strongest. Returns the
    per-partial list and its median; the broadband level flutter (`level_step_db_rms`) is the same measure on all partials at once."""
    lv = partial_levels(look, kmax)
    ok = steady_frames(look.sound) & np.isfinite(look.f0_t)  # the body of the note: the onset ramp is not flutter
    run = ok[1:] & ok[:-1]
    mean_db = np.array([10 * np.log10(np.mean(10 ** (l[ok] / 10)) + 1e-20) if ok.any() else -200 for l in lv])
    use = [k for k in range(kmax) if mean_db[k] > mean_db.max() - 35]
    steps = [float(np.sqrt(np.mean(np.diff(lv[k])[run] ** 2))) if run.any() else None for k in use]
    return dict(partials=[k + 1 for k in use], step_db_rms=steps, median_step_db_rms=float(np.median(steps)) if steps and None not in steps else None)


def attribute(y_real, y_ours, f0_hz):
    """The whole measured attribution of one pair's flux gap (see the module docstring). `y_real` and `y_ours` are the aligned signals
    the grader compared; like the grader, this cuts both to the length of the shorter."""
    n = min(len(y_real), len(y_ours))
    y_real, y_ours = y_real[:n], y_ours[:n]
    R, O = Look(y_real, f0_hz), Look(y_ours, f0_hz)
    base = flux_gap(R.flux, R.sound, O.flux, O.sound)
    nf = min(R.frames, O.frames)
    steady = steady_frames(R.sound[:nf] & O.sound[:nf])
    base_steady = flux_gap(R.flux, R.sound, O.flux, O.sound, steady)
    base["steady_frames"] = base_steady
    out = dict(baseline=base, real=track_stats(R), ours=track_stats(O), cells=dict(real=cell_shares(R), ours=cell_shares(O)),
               level_classes=dict(real=level_shares(R), ours=level_shares(O)), partial_flutter=dict(real=partial_flutter(R), ours=partial_flutter(O)))

    def variant(fr, sr, fo, so):
        """The graded numbers of the variant, and the same over the steady frames only (`steady_frames`), each with the share of
        the baseline gap it explains."""
        v = flux_gap(fr, sr, fo, so)
        v["gap_explained_share"] = float(1 - v["gap"] / base["gap"]) if base["gap"] > 0 else None
        s = flux_gap(fr, sr, fo, so, steady)
        s["gap_explained_share"] = float(1 - s["gap"] / base_steady["gap"]) if base_steady["gap"] > 0 else None
        v["steady_frames"] = s
        return v

    hr, ho = harmonic_cells(R), harmonic_cells(O)
    out["restricted_to_harmonic_cells"] = variant(class_flux(R, hr), R.sound, class_flux(O, ho), O.sound)
    out["restricted_to_between_cells"] = variant(class_flux(R, ~hr), R.sound, class_flux(O, ~ho), O.sound)
    rr, ro = ratio_track(R), ratio_track(O)
    sh_r, sh_o = shimmer_curve(R), shimmer_curve(O)
    fr_flat, sr_flat = _seen(retime(y_real, rr, "flatten"), f0_hz)
    fo_flat, so_flat = _seen(retime(y_ours, ro, "flatten"), f0_hz)
    out["f0_movement_removed_from_both"] = variant(fr_flat, sr_flat, fo_flat, so_flat)
    f_imp, s_imp = _seen(retime(y_ours, rr, "impose"), f0_hz)
    out["real_f0_movement_put_into_ours"] = variant(R.flux, R.sound, f_imp, s_imp)
    fr_g, sr_g = _seen(regain(y_real, sh_r, -1), f0_hz)
    fo_g, so_g = _seen(regain(y_ours, sh_o, -1), f0_hz)
    out["level_movement_removed_from_both"] = variant(fr_g, sr_g, fo_g, so_g)
    f_gi, s_gi = _seen(regain(y_ours, sh_r, +1), f0_hz)
    out["real_level_movement_put_into_ours"] = variant(R.flux, R.sound, f_gi, s_gi)
    f_b, s_b = _seen(regain(retime(y_real, rr, "flatten"), sh_r, -1), f0_hz)
    f_c, s_c = _seen(regain(retime(y_ours, ro, "flatten"), sh_o, -1), f0_hz)
    out["f0_and_level_movement_removed_from_both"] = variant(f_b, s_b, f_c, s_c)
    f_d, s_d = _seen(regain(retime(y_ours, rr, "impose"), sh_r, +1), f0_hz)
    out["real_f0_and_level_movement_put_into_ours"] = variant(R.flux, R.sound, f_d, s_d)
    # the residual (what is not at a multiple of the pitch): taken out of both, and the real's put into the render at the real's
    # harmonic-to-residual ratio; then the three measured components together
    hr_, nr_ = split_harmonic(y_real, R.f0_inst)
    ho_, no_ = split_harmonic(y_ours, O.f0_inst)
    f_h, s_h = _seen(hr_, f0_hz)
    f_ho, s_ho = _seen(ho_, f0_hz)
    out["harmonic_part_only_both"] = variant(f_h, s_h, f_ho, s_ho)
    g = float(np.sqrt(np.mean(ho_ ** 2) / max(np.mean(hr_ ** 2), 1e-30)))
    f_n, s_n = _seen(y_ours + g * nr_, f0_hz)
    out["real_residual_put_into_ours"] = variant(R.flux, R.sound, f_n, s_n)
    f_hh, s_hh = _seen(regain(retime(hr_, rr, "flatten"), sh_r, -1), f0_hz)
    f_hho, s_hho = _seen(regain(retime(ho_, ro, "flatten"), sh_o, -1), f0_hz)
    out["harmonic_part_with_f0_and_level_movement_removed_from_both"] = variant(f_hh, s_hh, f_hho, s_hho)
    f_all, s_all = _seen(regain(retime(y_ours + g * nr_, rr, "impose"), sh_r, +1), f0_hz)
    out["real_f0_level_and_residual_put_into_ours"] = variant(R.flux, R.sound, f_all, s_all)
    out["residual_to_harmonic_power_db"] = dict(
        real=float(10 * np.log10(np.mean(nr_ ** 2) / max(np.mean(hr_ ** 2), 1e-30))),
        ours=float(10 * np.log10(np.mean(no_ ** 2) / max(np.mean(ho_ ** 2), 1e-30))))
    # each partial's own level movement (not the broadband one): taken out of both, and the real's put into the render
    pd_r, pd_o = partial_deltas(R), partial_deltas(O)
    f_pr, s_pr = _seen(modulate_partials(y_real, R.f0_inst, pd_r, -1), f0_hz)
    f_po, s_po = _seen(modulate_partials(y_ours, O.f0_inst, pd_o, -1), f0_hz)
    out["partial_level_movement_removed_from_both"] = variant(f_pr, s_pr, f_po, s_po)
    f_pi, s_pi = _seen(modulate_partials(y_ours, O.f0_inst, pd_r, +1), f0_hz)
    out["real_partial_level_movement_put_into_ours"] = variant(R.flux, R.sound, f_pi, s_pi)
    # The gain operation adds flux of its own (a comb filter whose frames sit differently against the analysis frames: ~0.3 dB on a clean
    # tone for gains of +-3 dB, all of it in the near-floor cells). The control gives the render the same gain MAGNITUDES, constant in time:
    # whatever the control reads is the operation, not the movement.
    rms = np.sqrt(np.mean(pd_r[:, steady_frames(R.sound)] ** 2, axis=1, keepdims=True)) * np.ones_like(pd_r)
    f_pc, s_pc = _seen(modulate_partials(y_ours, O.f0_inst, rms, +1), f0_hz)
    out["control_constant_partial_gains_put_into_ours"] = variant(R.flux, R.sound, f_pc, s_pc)
    f_a4, s_a4 = _seen(modulate_partials(retime(y_ours + g * nr_, rr, "impose"), O.f0_inst * rr, pd_r, +1), f0_hz)
    out["real_f0_partial_level_and_residual_put_into_ours"] = variant(R.flux, R.sound, f_a4, s_a4)
    f_rt, s_rt = _seen(retime(retime(y_real, rr, "flatten"), rr, "impose"), f0_hz)
    out["roundtrip_control_real"] = dict(mean_flux_before=base["real"], mean_flux_after=float(np.mean(f_rt[s_rt])))
    return out


def _main(argv):
    """python3 flux_attribution.py FAMILY OUT_DIR [slot ...]: the attribution of every graded pair of a family in the analysis tree."""
    fam, out_dir = argv[1], Path(argv[2])
    only = set(argv[3:])
    out_dir.mkdir(parents=True, exist_ok=True)
    root = Path(dc.DEFAULTS["out"])
    if fam == "voice":
        lib, wav_dir, refs_json = ("/data/tmp/wa-inst/p2-grade/voice/library.json", "/data/tmp/wa-inst/refs-wav48-voice",
                                   "/data/tmp/wa-inst/p2-grade/voice/references-vocalset.json")
    else:
        lib, wav_dir, refs_json = dc.DEFAULTS["lib"], dc.DEFAULTS["wav_dir"], dc.DEFAULTS["refs_json"]
    refs = dc.Refs(lib, wav_dir, refs_json)
    for p in sorted((root / fam).glob("*/pair.json")):
        pj = json.loads(p.read_text())
        slot = p.parent.name
        if pj.get("refused") or not pj.get("ours") or pj.get("real") is None or (only and slot not in only):
            continue
        fname, rest = pj["decompose"]["library_note"].split("@", 1)
        name, dyn = rest.rsplit(".", 1)
        e = [x for x in refs.lib if x["file"] == fname and x["name"] == name and x["dyn"] == dyn][0]
        refs.band_hz = pj["alignment"].get("band_hz")
        f0 = dc.hz_of(e["note"])
        y_real, _ = refs.aligned(e)
        raw = dc.load(Path(pj["ours"]["wav"]))
        if refs.band_hz:
            raw = dc.lowpass(raw, refs.band_hz)
        y_ours, _ = dc.align(raw, f0, dc.render_floor(dc.highpass(raw, f0), f0))
        res = attribute(y_real, y_ours, f0)
        m = {x["id"]: x for x in pj["metrics"]}["flux_db"]
        res["pair"] = dict(family=fam, slot=slot, render_sha=pj["ours"]["git_sha"], reference=pj["real"]["file"],
                           graded=dict(real=m["real"], ours=m["ours"], gap=m["gap"]))
        (out_dir / f"{fam}-{slot}.json").write_text(json.dumps(res, indent=1))
        b = res["baseline"]
        print(f"{slot:11s} recomputed real {b['real']:.3f} ours {b['ours']:.3f} gap {b['gap']:.3f}   graded real {m['real']:.3f} ours {m['ours']:.3f} gap {m['gap']:.3f}", flush=True)


if __name__ == "__main__":
    _main(sys.argv)
