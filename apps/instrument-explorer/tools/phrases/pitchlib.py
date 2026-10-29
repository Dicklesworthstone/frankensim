"""The pitch detector of f0check.py as a library: f0(path) on the steady middle half, f0_samples(x, sr) on any span."""
import math
import re
import struct
import subprocess

import numpy as np


def resample(x, sr, to):
    """x at rate sr, at rate `to`, by ffmpeg's windowed-sinc resampler (filter_size 64): the same length and timing,
    flat to 15 kHz, and what lies above the new Nyquist removed (-103 dB at 30 kHz from 96 or 144 kHz). Plain linear
    interpolation, used before 2026-10-26, cost 3 dB at 15 kHz from 44.1 kHz and folded a 96 kHz clip's content above
    24 kHz back into the band at full level."""
    if sr == to:
        return x
    out = subprocess.run(["ffmpeg", "-nostdin", "-v", "error", "-f", "f32le", "-ar", str(sr), "-ac", "1", "-i", "-",
                          "-af", f"aresample={to}:filter_size=64", "-f", "f32le", "-ac", "1", "-"],
                         input=np.asarray(x, "<f4").tobytes(), capture_output=True, check=True).stdout
    return np.frombuffer(out, "<f4").astype(float)

def read_wav(path):
    """Mono samples at full scale +-1 and the rate, from 16/24/32-bit PCM or 32/64-bit float WAV, plain or extensible
    (the 48 kHz Iowa set is float, which the wave module refuses)."""
    with open(path, "rb") as f:
        b = f.read()
    if b[:4] != b"RIFF" or b[8:12] != b"WAVE":
        raise ValueError(f"{path}: not a WAV file")
    fmt = data = None
    i = 12
    while i + 8 <= len(b):
        cid, size = b[i:i + 4], int.from_bytes(b[i + 4:i + 8], "little")
        if cid == b"fmt ":
            fmt = b[i + 8:i + 8 + size]
        elif cid == b"data":
            data = b[i + 8:i + 8 + size]
        i += 8 + size + (size & 1)
    if fmt is None or data is None:
        raise ValueError(f"{path}: no fmt or data chunk")
    tag, ch, sr = struct.unpack("<HHI", fmt[:8])
    bits = struct.unpack("<H", fmt[14:16])[0]
    if tag == 0xFFFE:  # extensible: the real format is the sub-format GUID's first two bytes
        tag = struct.unpack("<H", fmt[24:26])[0]
    if tag == 3 and bits in (32, 64):
        x = np.frombuffer(data[: len(data) // (bits // 8) * (bits // 8)], "<f4" if bits == 32 else "<f8").astype(float)
    elif tag == 1 and bits == 16:
        x = np.frombuffer(data[: len(data) // 2 * 2], "<i2") / 32768.0
    elif tag == 1 and bits == 32:
        x = np.frombuffer(data[: len(data) // 4 * 4], "<i4") / 2147483648.0
    elif tag == 1 and bits == 24:
        u = np.frombuffer(data[: len(data) // 3 * 3], np.uint8).reshape(-1, 3).astype(np.int32)
        v = u[:, 0] | (u[:, 1] << 8) | (u[:, 2] << 16)
        x = np.where(v >= 1 << 23, v - (1 << 24), v) / 8388608.0
    else:
        raise ValueError(f"{path}: WAV format {tag}, {bits}-bit")
    n = len(x) // ch * ch
    return x[:n].reshape(-1, ch).mean(axis=1), sr


def highpass(x, sr, fc=25.0, top=None):
    """x without its DC offset and rumble below fc (a raised-cosine ramp from fc/2 to fc), by FFT. Iowa takes carry
    both, up to 30 dB above a low note, and a pitch search keyed to the autocorrelation's first zero crossing needs
    them gone. Given top, also without what lies above it (a raised-cosine roll-off from top to 1.5 top)."""
    if not len(x):
        return x
    n = 1 << (len(x) + sr // 10).bit_length()
    f = np.fft.rfftfreq(n, 1 / sr)
    g = np.clip((f - fc / 2) / (fc / 2), 0.0, 1.0)
    gain = 0.5 - 0.5 * np.cos(np.pi * g)
    if top:
        gain *= 0.5 + 0.5 * np.cos(np.pi * np.clip((f - top) / (0.5 * top), 0.0, 1.0))
    return np.fft.irfft(np.fft.rfft(x, n) * gain, n)[: len(x)]


def f0(path):
    x, sr = read_wav(path)
    return f0_samples(x[len(x) // 4: 3 * len(x) // 4], sr)


def f0_samples(x, sr):
    if len(x) > sr:  # the middle second is plenty for a pitch, and a long span makes every FFT below slow
        x = x[(len(x) - sr) // 2:][:sr]
    x = highpass(x - x.mean(), sr)
    if not np.any(x):
        return None, 0.0
    # Candidate periods come from the autocorrelation: its highest peak, the exact divisors of that lag and its
    # multiples, each kept when its own peak is within 10 % of the highest. Time-domain periodicity alone is fooled
    # both ways: a clean tone's autocorrelation peaks as high at 2T as at T, and a low bowed note whose waveform is
    # dominated by a body resonance near its 3rd to 7th harmonic repeats almost as well at that harmonic's period.
    # So the spectrum decides, by harmonic contrast: how far the spectrum at h*f (h = 1..8) stands above the spectrum
    # midway between them. The answer is the lowest candidate with at least 70 % of the best contrast. A subharmonic
    # f/k of a real tone has peaks at only every k-th slot, so its contrast is about 1/k of the true one; a low note
    # with a weak fundamental (a bowed bass 26-39 dB down on its body-boosted 5th harmonic) keeps its full contrast.
    # Before scoring, any candidate whose lines all sit at multiples of k is lifted to k times itself (lift below):
    # when every candidate's contrast is small, 1/k of it is not small enough to lose.
    nfft = 1 << (2 * len(x) - 1).bit_length()  # a power of two: an odd length can cost 50 times as much
    spec = np.fft.rfft(x, nfft)
    ac = np.fft.irfft(spec * np.conj(spec), nfft)[: len(x)]
    ac /= ac[0]
    lo, hi = int(sr / 2500), int(sr / 30)
    # Search past the first zero crossing. The autocorrelation is still falling from its zero-lag peak at the shortest
    # lag searched, and for a low note it is highest there (a breathy cello E2 read 2595 Hz); a zero-mean periodic
    # signal's autocorrelation always crosses zero within one period, so nothing real is lost.
    neg = np.flatnonzero(ac[1:hi] <= 0)
    if len(neg):
        lo = max(lo, int(neg[0]) + 1)
    big = lo + int(np.argmax(ac[lo:hi]))
    best = float(ac[big])

    def peak_near(c):
        w = max(2, int(0.03 * c))
        i0, i1 = int(round(c)) - w, int(round(c)) + w + 1
        if i0 < lo or i1 >= min(hi * 4, len(ac)):
            return None
        j = i0 + int(np.argmax(ac[i0:i1]))
        return j if i0 < j < i1 - 1 else None

    cands = {big}
    for k in range(2, 9):
        for c in (big / k, big * k):
            j = peak_near(c)
            # A longer period only has to be clearly periodic (the spectral test below rejects subharmonics).
            if j is not None and ac[j] >= (0.9 if c < big else 0.5) * best:
                cands.add(j)
    n = 1 << 18
    seg = x[:n]
    mag = np.abs(np.fft.rfft(seg * np.hanning(len(seg)), n))
    freqs = np.fft.rfftfreq(n, 1 / sr)
    db = 20 * np.log10(mag / mag.max() + 1e-12)

    db = np.maximum(db, -80.0)
    df = float(freqs[1])

    def window(fc, rel):
        """db over [fc (1 - rel), fc (1 + rel)], by index (the bins are uniform)."""
        i0 = max(0, int(math.ceil(fc * (1 - rel) / df)))
        i1 = min(len(db) - 1, int(fc * (1 + rel) / df))
        return db[i0:i1 + 1]

    def level_at(f, rel=0.02):
        w = window(f, rel)
        return float(w.max()) if len(w) else -120.0

    def valley(fc):
        w = window(fc, 0.02)
        return float(np.median(w)) if len(w) else -80.0

    def contrast(f):
        """Mean level at the harmonics h*f (h = 1..8) above the mean level midway between them."""
        peaks, valleys = [], []
        for hnum in range(1, 9):
            if (hnum + 0.5) * f >= sr / 2:
                break
            peaks.append(level_at(hnum * f))
            valleys.append(valley((hnum + 0.5) * f))
        return float(np.mean(peaks) - np.mean(valleys)) if peaks else 0.0

    def peaked(f):
        """The harmonic slots h = 1..12 of f that carry a real line: within 45 dB of the strongest, and 15 dB or more
        above the spectrum midway to both neighbouring slots (the maximum of plain noise over a +-2 % window stands
        about 10 dB above its median, so a lower bar counts noise as harmonics)."""
        out = []
        for hnum in range(1, 13):
            if (hnum + 0.5) * f >= sr / 2:
                break
            lv = level_at(hnum * f, 0.01)  # +-1 %: +-2 % lets h >= 5 match by chance
            if lv >= -45 and lv - max(valley((hnum - 0.5) * f), valley((hnum + 0.5) * f)) >= 15:
                out.append(hnum)
        return out

    def has_fundamental(j):
        # A bowed bass's fundamental can sit 26-39 dB under its body-boosted 5th harmonic and is still there; a
        # subharmonic has nothing at all at its own frequency.
        return level_at(sr / j) >= -50

    def backed(j):
        """The spectrum supports this period: a line at its own frequency, or lines that are not all multiples of
        some k >= 2 (a stopped pan pipe's D3 is 5, 7, 9, 11 x 145 Hz with nothing at 145 Hz, and 145 Hz is its
        period)."""
        if has_fundamental(j):
            return True
        s = peaked(sr / j)
        return bool(s) and math.gcd(*s) == 1

    def lift(j):
        """A period with no line at its own frequency whose lines all sit at multiples of k is k times too long (a
        recorder note with lines at 670 and 2012 Hz reads as 134 Hz, 670/5, when only the contrast decides)."""
        for _ in range(3):
            if has_fundamental(j):
                break
            s = peaked(sr / j)
            g = math.gcd(*s) if s else 1
            if g <= 1:
                break
            jj = peak_near(j / g)
            j = jj if jj is not None else max(1, int(round(j / g)))
        return j

    real = {j for j in map(lift, cands) if 1 <= j < len(ac) - 1 and backed(j)}
    if not real:
        # Nothing the autocorrelation offers is backed (a squeal whose waveform repeats only over many cycles):
        # take the strongest line.
        band = (freqs >= 30) & (freqs <= 2500)
        if band.any():
            fsx = float(freqs[band][int(np.argmax(db[band]))])
            j = peak_near(sr / fsx)
            real = {j if j is not None else int(round(sr / fsx))}
        real = {j for j in real if 1 <= j < len(ac) - 1} or {big}
    scored = {j: contrast(sr / j) for j in real if sr / j >= 30}
    top = max(scored.values()) if scored else 0.0

    def subharmonic(j):
        """j is k times the period of another candidate and its slots off the multiples of k sit 20 dB or more under
        the slots on them: a weak series at f/k (a real cello B2 carried one 26-39 dB down, and its contrast alone
        was 79 % of the true pitch's). A low note's weak fundamental does not trip this: its other odd slots, the
        body-boosted 3rd and 5th among them, are strong."""
        f = sr / j
        for j2 in scored:
            k = round(j / j2)
            if k < 2 or abs(j / j2 - k) > 0.03 * k:
                continue
            on = [level_at(h * f, 0.01) for h in range(1, 9) if h % k == 0 and (h + 0.5) * f < sr / 2]
            off = [level_at(h * f, 0.01) for h in range(1, 9) if h % k and (h + 0.5) * f < sr / 2]
            if on and off and np.mean(on) - np.mean(off) >= 20:
                return True
        return False

    lag = big
    for j in sorted(scored, reverse=True):  # longest period (lowest pitch) first
        if scored[j] >= 0.7 * top and not subharmonic(j):
            lag = j
            break
    a, b, cc = ac[lag - 1], ac[lag], ac[lag + 1]
    # parabolic refinement; a flat or inverted neighbourhood can throw the vertex anywhere, so keep it within the bin
    shift = 0.5 * (a - cc) / (a - 2 * b + cc) if (a - 2 * b + cc) < 0 else 0.0
    shift = max(-0.5, min(0.5, shift))
    return sr / (lag + shift), float(b)


def midi_of(note):
    """'F#5' / 'Bb4' / 'D3' -> MIDI number, or None for a non-note target."""
    m = re.fullmatch(r"([A-G])([#b]?)(-?\d)", note or "")
    if not m:
        return None
    base = {"C": 0, "D": 2, "E": 4, "F": 5, "G": 7, "A": 9, "B": 11}[m.group(1)]
    return base + {"#": 1, "b": -1, "": 0}[m.group(2)] + 12 * (int(m.group(3)) + 1)


