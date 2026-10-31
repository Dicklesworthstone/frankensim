"""Measure every phrase and a matching real-recording phrase, the same way (run on trj after phrasekit.py).

usage: python3 analyze.py PHRASES_DIR REFS_WAV_DIR SCREENED_CLIPS_DIR SCREENED.json [BANK.json]
  PHRASES_DIR      phrasekit.py's output (phrases.json, <id>.wav, <id>.webp); written to in place
  REFS_WAV_DIR     University of Iowa MIS chromatic runs as 48 kHz WAV (one file per run, notes separated by silence)
  SCREENED_*       the screened reference clips (references-screened.json and clips/<family>/<id>.wav), used for
                   families the Iowa library lacks (recorder, pan pipes, tin whistle, voice)
  BANK.json        cache of the segmented runs (built on first use)

For each phrase whose notes are known, a reference phrase is assembled from real single-note recordings of the
same notes: each note's recording is cut from its onset to the simulation note's slot and placed at the same time,
joined the same way (a note on the same string stops the one before). It is a splice of real notes, not a
performance of the phrase, and the page says so. Phrases that already carry a real recording (the kit groove) use it.

Both tracks are then measured by the same code:
  frames (every 10 ms)  level (A-weighted RMS, dB), brightness (spectral centroid), harmonicity (HNR from the
                        normalised autocorrelation), noisiness (spectral flatness), pitch (and cents from the note
                        asked for), roughness (share of energy away from the harmonics of the frame's pitch), onsets
                        (spectral flux)
  notes                 the same, summarised over each note's steady part, plus attack time, onset roughness (first
                        120 ms), the first 16 harmonic levels, partial stretch, the noise floor between harmonics,
                        vibrato/jitter, decay rate, and a timbre distance between the two tracks (third-octave bands)
  phrase                long-term average spectrum (third-octave), clipping, crest factor, and a summary of the
                        per-note differences
Writes PHRASES_DIR/<id>.analysis.json, <id>-ref.wav and <id>-ref.webp, and adds `analysis` and `ref` to phrases.json.
"""
import json
import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor
import wave
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from phrasekit import spectrogram  # noqa: E402
from pitchlib import f0_samples, highpass, midi_of, read_wav, resample  # noqa: E402

SR = 48000
HOP = 480  # 10 ms
WIN = 2048
PH, REFS, CLIPS, SCREENED = (Path(a) for a in sys.argv[1:5])
BANK = Path(sys.argv[5]) if len(sys.argv) > 5 else PH / "refbank.json"
NAMES = ["C", "Db", "D", "Eb", "E", "F", "Gb", "G", "Ab", "A", "Bb", "B"]
IOWA = {"violin": "Violin", "viola": "Viola", "cello": "Cello", "double-bass": "Bass", "flute": "Flute",
        "clarinet": "BbClar", "trumpet": "Trumpet", "trombone": "TenorTrombone", "oboe": "Oboe",
        "bassoon": "Bassoon", "horn": "Horn", "tuba": "Tuba", "saxophone": "AltoSax", "guitar": "Guitar"}
IOWA_CREDIT = "University of Iowa Electronic Music Studios, Musical Instrument Samples (free to use without restriction)"
# The lowest note each instrument can sound (open string, lowest key, lowest sung note). A real recording is cleaned
# of what lies below 3/4 of it: the Iowa takes carry DC offsets and room rumble up to 30 dB above a low note, which
# the splice turns into thumps at every join and which the pitch search reads as the note.
FLOOR_HZ = {"violin": 196.0, "viola": 130.81, "cello": 65.41, "double-bass": 41.2, "flute": 261.63,
            "clarinet": 146.83, "trumpet": 164.81, "trombone": 82.41, "oboe": 233.08, "bassoon": 58.27,
            "horn": 61.74, "tuba": 36.71, "saxophone": 138.59, "guitar": 82.41, "recorder": 261.63,
            "tin-whistle": 587.33, "pan-pipes": 98.0, "voice": 82.41}


def hz_of(midi):
    return 440.0 * 2 ** ((midi - 69) / 12)


def name_of(midi):
    return f"{NAMES[midi % 12]}{midi // 12 - 1}"


def wav_mono(path):
    x, sr = read_wav(path)
    return resample(x, sr, SR)  # the Iowa clips are 44.1, 32 and 96 kHz, a material-strings render 144 kHz


def write_wav(path, x):
    peak = np.max(np.abs(x)) or 1.0
    with wave.open(str(path), "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(SR)
        w.writeframes((x / peak * 10 ** (-1 / 20) * 32767).astype("<i2").tobytes())


def rms_env(x, hop=HOP):
    n = len(x) // hop
    return np.sqrt(np.mean(x[: n * hop].reshape(n, hop) ** 2, axis=1) + 1e-20)


# ------------------------------------------------------------------ the reference bank
def run_range(name):
    """The MIDI span an Iowa run's name declares ('Cello.arco.mf.sulC.C2Ab2.wav' -> 36..44), or None."""
    for tok in name.split("."):
        m = re.fullmatch(r"([A-G]b?-?\d)([A-G]b?-?\d)?", tok)
        if m:
            lo = midi_of(m.group(1))
            return lo, midi_of(m.group(2)) if m.group(2) else lo
    return None


def segment_run(path):
    """Notes of one Iowa chromatic run: stretches above -40 dB of the run's peak, apart by >= 150 ms of quiet.

    A run plays every semitone of the span its name declares, in order, so when the count matches, the k-th note is
    the span's k-th semitone and the measured pitch only has to agree with it (to within half a semitone, in any
    octave: the detector can still mistake the octave of a low note, which labelled a D3 as D2 and a C4 as C3).
    Otherwise each note's measured pitch is moved by octaves into the span, and dropped if it cannot be."""
    fam = next((k for k, v in IOWA.items() if v == path.name.split(".")[0]), None)
    x = wav_mono(path)
    if fam in FLOOR_HZ:
        x = highpass(x, SR, 0.75 * FLOOR_HZ[fam])
    env = rms_env(x)
    env = 20 * np.log10(env / (env.max() or 1))
    loud = env > -40
    segs, i, n = [], 0, len(loud)
    while i < n:
        if not loud[i]:
            i += 1
            continue
        j = i
        quiet = 0
        while j < n and quiet < 15:
            quiet = quiet + 1 if not loud[j] else 0
            j += 1
        a, b = i, j - quiet
        if b - a >= 25:
            segs.append((a * HOP / SR, b * HOP / SR))
        i = j
    span = run_range(path.name)
    by_order = span is not None and len(segs) == span[1] - span[0] + 1
    out = []
    for k, (a, b) in enumerate(segs):
        mid = x[int((a + 0.25 * (b - a)) * SR): int((a + 0.75 * (b - a)) * SR)]
        f, clarity = f0_samples(mid, SR)
        if not f or f <= 0 or clarity < 0.5:
            continue
        m = 69 + 12 * math.log2(f / 440)
        if by_order:
            want = span[0] + k
        elif span is not None:
            fits = [round(m) + 12 * j for j in range(-3, 4) if span[0] <= round(m) + 12 * j <= span[1]]
            if not fits:
                continue
            want = min(fits, key=lambda w: abs(w - m))
        else:
            want = round(m)
        d = (m - want + 6) % 12 - 6  # the measured pitch's distance from the label, ignoring octaves
        if abs(d) > 0.4:
            continue
        out.append(dict(start=round(a, 3), end=round(b, 3), midi=int(want), cents=round(100 * d, 1)))
    return out


def build_bank():
    if BANK.exists():
        return json.loads(BANK.read_text())
    bank = []
    runs = [p for p in sorted(REFS.glob("*.wav")) if p.name.split(".")[0] in IOWA.values()]
    # one process per run file: segmenting 600-odd runs one after another takes most of an hour on a busy host
    with ProcessPoolExecutor(max_workers=min(12, os.cpu_count() or 4)) as pool:
        for p, segs in zip(runs, pool.map(segment_run, runs)):
            fam = next(k for k, v in IOWA.items() if v == p.name.split(".")[0])
            low = p.name.lower()
            dyn = next((d for d in ("pp", "mf", "ff") if f".{d}." in low), "")
            m = re.search(r"\.sul([a-g])\.", low)
            tech = "pizz" if ".pizz." in low else ("arco" if ".arco." in low else ("vib" if ".vib." in low else "plain"))
            for s in segs:
                bank.append(dict(family=fam, file=p.name, dyn=dyn, string=m.group(1).upper() if m else None, tech=tech,
                                 mono=".mono." in low, **s))
    BANK.write_text(json.dumps(bank))
    return bank


def screened_clips():
    """Valid screened clips by family and MIDI note, with their licence and author (publish.py drops the audio of
    any whose licence forbids redistribution)."""
    d = json.loads(SCREENED.read_text())
    out = {}
    for f in d["files"]:
        for clips in [f.get("clips") or []] + [m.get("clips") or [] for m in f.get("members") or []]:
            for c in clips:
                m = midi_of(c.get("note"))
                p = CLIPS / c["family"] / f"{c['id']}.wav"
                # only recordings the public site may play: never the all-rights-reserved YouTube licence
                if c.get("status") == "VALID" and m is not None and p.exists() and "standard" not in f["licence"]:
                    out.setdefault((c["family"], m), []).append(dict(path=str(p), licence=f["licence"],
                                                                     author=f.get("author") or "", id=c["id"]))
    return out


def pick(bank, screened, fam, midi, string):
    """A recording of exactly this note: an Iowa arco/plain mf segment, on the same string when the library has one."""
    tech_ok = ("arco",) if fam in ("violin", "viola", "cello", "double-bass") else ("plain", "vib", "arco")
    cands = [b for b in bank if b["family"] == fam and b["midi"] == midi and b["tech"] in tech_ok]
    if fam == "flute":  # the Iowa flute has vibrato and non-vibrato runs; the simulation has no vibrato
        cands.sort(key=lambda b: b["tech"] == "vib")
    if cands:
        def score(b):
            return ((string is not None and b["string"] != string) * 4 + (b["dyn"] != "mf") * 2
                    + (not b["mono"]) * 1 + (b["tech"] == "vib") * 1 + abs(b["cents"]) / 100)
        b = min(cands, key=score)
        x = wav_mono(REFS / b["file"])[int(b["start"] * SR): int(b["end"] * SR)]
        return x, dict(note=name_of(midi), source=f"Iowa MIS {b['file'].rsplit('.', 1)[0]}", licence="UIowa",
                       credit=IOWA_CREDIT)
    for c in screened.get((fam, midi), []):
        x = wav_mono(c["path"])
        return x, dict(note=name_of(midi), source=c["id"], licence=c["licence"], credit=c["author"])
    return None, None


def onset(x):
    """First sample within 30 dB of the clip's loudest 10 ms, less 5 ms."""
    env = rms_env(x)
    k = int(np.argmax(env >= env.max() * 10 ** (-30 / 20)))
    return max(0, k * HOP - SR // 200)


def reference_phrase(p, bank, screened):
    """The phrase rebuilt from real single-note recordings, joined like the simulation, and the indices of the notes no
    recording covers (left silent: the open-licence whistle, pan pipe and voice clips cover only some notes). None if
    a note is not a note name or fewer than two notes have a recording."""
    if p.get("real"):
        return None
    notes = p.get("notes") or []
    if not notes or p["family"] not in set(IOWA) | {"recorder", "pan-pipes", "tin-whistle", "voice"}:
        return None
    tracks, used, missing = [], [], []
    for i, n in enumerate(notes):
        parts = n["note"].split("+") if "+" in n["note"] else [n["note"].split()[-1]]
        strings = (n.get("string") or "").split("+")
        if any(midi_of(nm) is None for nm in parts):
            return None
        picks = []
        for k, nm in enumerate(parts):
            s = strings[k] if k < len(strings) and strings[k] else None
            picks.append(pick(bank, screened, p["family"], midi_of(nm), s[0] if s else None))
        if any(x is None for x, _ in picks):  # a double stop needs both of its notes
            missing.append(i)
            continue
        xs = []
        for x, info in picks:
            used.append(info)
            x = highpass(x, SR, 0.75 * FLOOR_HZ.get(p["family"], 40.0))
            x = x[onset(x):][: int((n["dur"] + 0.04) * SR)].copy()
            f = min(len(x), int(0.04 * SR))
            x[-f:] *= np.linspace(1, 0, f)
            xs.append(x / (np.max(np.abs(x)) or 1))
        y = np.zeros(max(len(x) for x in xs))
        for x in xs:
            y[: len(x)] += x / len(xs)
        tracks.append((set(strings) if n.get("string") else {None}, int(round(n["at"] * SR)), y))
    if len(tracks) < 2:
        return None
    total = max(i + len(y) for _, i, y in tracks)
    out = np.zeros(total + int(0.3 * SR))
    for k, (held, i, y) in enumerate(tracks):
        # a later note on the same string (or, for winds and voice, any later note) ends this one, as in phrasekit
        end = len(y)
        for held2, i2, _ in tracks[k + 1:]:
            if held & held2 and i2 > i:
                end = min(end, i2 - i)
                break
        seg = y[:end].copy()
        if end < len(y):
            f = min(len(seg), int(0.04 * SR))
            seg[-f:] *= np.linspace(1, 0, f)
        out[i: i + len(seg)] += seg
    return out, used, missing


# ------------------------------------------------------------------ measurement
A_W = None


def a_weight(freqs):
    f2 = freqs ** 2
    ra = (12194 ** 2 * f2 ** 2) / ((f2 + 20.6 ** 2) * np.sqrt((f2 + 107.7 ** 2) * (f2 + 737.9 ** 2)) * (f2 + 12194 ** 2))
    return ra / ra[np.argmin(np.abs(freqs - 1000))]


def frames(x, targets, f_lo, anchors, chords=lambda t: None):
    """Per-frame measures. targets(t) -> wanted Hz or None; anchors(t) -> the note's own measured pitch (from the long
    window, which can tell a subharmonic from the real tone), used to settle a frame's octave; chords(t) -> the
    fundamentals of every note of a chord sounding at t, or None.

    Roughness is the share of power off the harmonic comb: the comb of every note of a chord (the second note's
    partials are off the first one's comb, so a single comb read the other note as noise), and on a frame with no
    reliable pitch (a scraping onset) the comb of the note's own steady pitch rather than of a lag that does not
    repeat, so the frames where a start is rough are measured against the note it settles into.

    Pitch and harmonicity come from the autocorrelation of a longer frame holding three periods of the lowest pitch
    (85 ms at 50 Hz), zero-padded so it does not wrap, and divided by the autocorrelation of its own window (Boersma
    1993, as Praat does): a 43 ms Hann frame keeps only half of a 65 Hz tone's correlation at one period, and its
    autocorrelation is still falling from the zero-lag peak at the shortest lag searched, so the highest value was the
    shortest lag (2526 Hz) on every low note. The search starts past the first zero crossing."""
    global A_W
    win = np.hanning(WIN)
    pad = np.concatenate([np.zeros(WIN // 2), x, np.zeros(WIN)])
    nfr = 1 + (len(pad) - WIN) // HOP
    freqs = np.fft.rfftfreq(WIN, 1 / SR)
    if A_W is None:
        A_W = a_weight(np.maximum(freqs, 1.0))
    band = (freqs >= 100) & (freqs <= 10000)
    lag_lo, lag_hi = int(SR / 2500), int(SR / max(f_lo, 30))
    pwl = 1 << math.ceil(math.log2(max(WIN, 3 * lag_hi)))
    pwin = np.hanning(pwl)
    rw = np.fft.irfft(np.abs(np.fft.rfft(pwin, 2 * pwl)) ** 2)[:pwl]
    rw = np.maximum(rw / rw[0], 1e-3)
    lags = np.maximum(np.arange(pwl), 1).astype(float)
    # Pitch and harmonicity are read below 6 kHz. A flue render whose strongest lines sat at 19-22 kHz (its 31st
    # harmonic) correlated better at five periods than at one, because a whole-sample lag misses a 20 kHz line's peak
    # by up to a quarter of its cycle: a recorder F4 sounding 669 Hz read as 134 Hz, and one sounding 1418 Hz (its third
    # mode) read 34 cents from the Bb4 asked for.
    ppad = np.concatenate([np.zeros(pwl // 2), highpass(x, SR, top=6000), np.zeros(pwl)])
    out = {k: [] for k in ("level", "centroid", "hnr", "flat", "f0", "cents", "rough", "flux")}
    prev = None
    peak = np.max(np.abs(x)) or 1.0
    for k in range(nfr):
        seg = pad[k * HOP: k * HOP + WIN]
        spec = np.fft.rfft(seg * win)
        pw = np.abs(spec) ** 2
        lvl = 10 * np.log10(np.sum(pw * A_W ** 2) / (np.sum(win ** 2) * peak ** 2 * WIN / 2) + 1e-14)
        out["level"].append(round(max(lvl, -100.0), 1))
        tot = pw[freqs >= 50].sum()
        quiet = lvl < -60
        out["centroid"].append(None if quiet or tot <= 0 else round(float((freqs * pw)[freqs >= 50].sum() / tot)))
        b = pw[band] + 1e-20
        out["flat"].append(None if quiet else round(float(10 * np.log10(np.exp(np.mean(np.log(b))) / np.mean(b))), 1))
        hnr = f0 = cents = rough = None
        if not quiet:
            ps = ppad[k * HOP: k * HOP + pwl]
            ac = np.fft.irfft(np.abs(np.fft.rfft((ps - ps.mean()) * pwin, 2 * pwl)) ** 2)[:pwl]
            ac = ac / (ac[0] or 1) / rw
            neg = np.flatnonzero(ac[1:lag_hi] <= 0)
            start = max(lag_lo, int(neg[0]) + 1) if len(neg) else lag_hi
        if not quiet and start < lag_hi - 1:
            # Corrected for the window, a steady tone correlates about equally at T, 2T, 3T, 4T (0.9992 against
            # 0.9999 on a trumpet D5), so the plain maximum is often a multiple. Praat's octave cost (0.01 per
            # octave of lag) makes the shortest well-correlated period win.
            j = start + int(np.argmax(ac[start:lag_hi] - 0.01 * np.log2(lags[start:lag_hi])))
            r = float(np.clip(ac[j], 1e-4, 0.9999))
            hnr = round(10 * math.log10(r / (1 - r)), 1)
            want = targets(k * HOP / SR)
            anchor = anchors(k * HOP / SR)
            if r > 0.5:
                # The octave cost can still lose to a period double (or win against a real one). Settle the octave
                # toward the note's own measured pitch, never toward the note asked for: that would pull a wrong note
                # an octave nearer the target and hide the error.
                if anchor:
                    j0, best = j, abs(math.log2(SR / j / anchor))
                    for mult in (2, 3, 4, 1 / 2, 1 / 3, 1 / 4):
                        c = int(round(j0 * mult))
                        if not start + 2 <= c < lag_hi - 2:
                            continue
                        jj = c - 2 + int(np.argmax(ac[c - 2: c + 3]))  # the peak itself, not the rounded lag
                        dist = abs(math.log2(SR / jj / anchor))
                        if ac[jj] > 0.9 * r and dist < best:
                            j, best = jj, dist
                jj = min(max(j, 1), pwl - 2)
                a_, b_, c_ = ac[jj - 1], ac[jj], ac[jj + 1]
                # parabolic refinement (one lag is 36 cents at 1 kHz), only where the neighbourhood is a real peak
                shift = 0.5 * (a_ - c_) / (a_ - 2 * b_ + c_) if (a_ - 2 * b_ + c_) < 0 else 0.0
                f0v = SR / (jj + max(-0.5, min(0.5, shift)))
                f0 = round(f0v, 1)
                cents = round(1200 * math.log2(f0v / want), 1) if want else None
                harm = np.zeros(len(freqs), bool)
                for g in chords(k * HOP / SR) or [f0v if r > 0.5 or not anchor else anchor]:
                    hn = np.maximum(np.round(freqs / g), 1)
                    harm |= np.abs(freqs - hn * g) <= np.maximum(0.03 * hn * g, 1.5 * SR / WIN)
                # Below about 70 Hz (and on a low chord) the comb's tolerance, one and a half bins of this window,
                # closes every gap between partials up to 5 kHz: there is nothing off the comb to measure, not zero.
                if harm[(freqs >= 50) & (freqs <= 5000)].mean() < 0.9:
                    rough = round(float(pw[(freqs >= 50) & ~harm].sum() / (tot or 1)), 3)
        out["hnr"].append(hnr)
        out["f0"].append(f0)
        out["cents"].append(cents)
        out["rough"].append(rough)
        lm = np.log10(np.abs(spec) + 1e-9)
        out["flux"].append(0.0 if prev is None else round(float(np.sum(np.maximum(lm - prev, 0)) / len(lm) * 100), 2))
        prev = lm
    return out


def third_octave(x, lo=31.5, hi=16000):
    centres = [1000 * 2 ** (k / 3) for k in range(-15, 13) if lo <= 1000 * 2 ** (k / 3) <= hi]
    n = 1 << 15
    pw = np.zeros(n // 2 + 1)
    for i in range(0, max(1, len(x) - n), n // 2):
        seg = x[i: i + n]
        if len(seg) < n:
            seg = np.pad(seg, (0, n - len(seg)))
        pw += np.abs(np.fft.rfft(seg * np.hanning(n))) ** 2
    freqs = np.fft.rfftfreq(n, 1 / SR)
    lv = []
    for c in centres:
        sel = (freqs >= c / 2 ** (1 / 6)) & (freqs < c * 2 ** (1 / 6))
        lv.append(10 * math.log10(pw[sel].sum() + 1e-20))
    top = max(lv)
    return [round(c, 1) for c in centres], [round(v - top, 1) for v in lv]


def partials(seg, f0, count=16):
    """Levels (dB re the strongest) and frequencies of the first `count` harmonics of seg, and the noise floor
    between them (dB re the strongest)."""
    n = 1 << 16
    s = seg[:n]
    mag = np.abs(np.fft.rfft(s * np.hanning(len(s)), n))
    freqs = np.fft.rfftfreq(n, 1 / SR)
    lv, fr, valleys = [], [], []
    for hn in range(1, count + 1):
        c = hn * f0
        if c > 0.45 * SR:
            break
        sel = (freqs >= c * 0.97) & (freqs <= c * 1.03)
        i = np.flatnonzero(sel)
        j = i[np.argmax(mag[i])]
        lv.append(float(mag[j]))
        fr.append(float(freqs[j]))
        v = (freqs >= (hn + 0.35) * f0) & (freqs <= (hn + 0.65) * f0)
        if v.any():
            valleys.append(float(np.median(mag[v])))
    top = max(lv) if lv else 1.0
    db = [round(20 * math.log10(v / top + 1e-9), 1) for v in lv]
    floor = round(20 * math.log10(np.median(valleys) / top + 1e-9), 1) if valleys else None
    stretch = None  # mean cents of partials 6..12 from exact multiples (stiffness or dispersion makes it positive)
    hi = [(k + 1, f) for k, f in enumerate(fr) if 6 <= k + 1 <= 12 and db[k] > -50]
    if len(hi) >= 3:
        stretch = round(float(np.mean([1200 * math.log2(f / (k * f0)) for k, f in hi])), 1)
    return db, floor, stretch


def mean(v):
    v = [a for a in v if a is not None]
    return round(float(np.mean(v)), 2) if v else None


def common_period(fs):
    """A chord's period: the shortest multiple of its lowest note's period (up to 6) in which every other note completes
    a whole number of cycles to within 6 % of one (equal-tempered fifths, fourths and thirds all do by 5 periods)."""
    lo = min(fs)
    for m in range(1, 7):
        T = m / lo
        if all(abs(T * f - round(T * f)) < 0.06 for f in fs):
            return T
    return None


def speaks_ms(x, a, d, period):
    """How long a note takes to speak: from its start until the sound repeats at its period (the normalised correlation
    of a window with the next one a period on, the period searched within 2 %, above 0.9) and keeps repeating for
    50 ms. A scraping or multiply-slipping bowed start, a breathy flue transient and a reed that has not locked are all
    aperiodic: this is what a spectrogram shows as smear before the lines. Returns (ms, settled); a note that never
    settles returns its own length and False."""
    P = period * SR
    win = int(max(0.06, 2 * period) * SR)
    i0, i1 = int(a * SR), int(min(len(x), (a + d) * SR))
    run = 0
    lags = range(int(0.98 * P), int(1.02 * P) + 2)
    for i in range(i0, i1 - win - lags[-1], SR // 100):
        u = x[i:i + win] - x[i:i + win].mean()
        nu = np.linalg.norm(u)
        best = 0.0
        if nu > 0:
            for L in lags:
                v = x[i + L:i + L + win]
                v = v - v.mean()
                best = max(best, float(np.dot(u, v)) / (nu * (np.linalg.norm(v) or 1)))
        run = run + 1 if best > 0.9 else 0
        if run == 5:
            return round((i - i0) / SR * 1000 - 40), True
    return round(d * 1000), False


def note_measures(x, fr, n, family):
    a, d = n["at"], n["dur"]
    i0, i1 = int(a * SR), int(min(len(x), (a + d) * SR))
    k0, k1 = int(a * SR / HOP), int((a + d) * SR / HOP)
    lead = min(0.12, 0.3 * d)
    s0, s1 = int((a + lead) / (HOP / SR)), int(max(a + lead + 0.05, a + d - 0.04) / (HOP / SR))
    steady = slice(s0, max(s0 + 1, s1))
    env = np.array(fr["level"][k0:max(k0 + 1, k1)])
    out = {}
    if len(env):
        top = env.max()
        k90 = int(np.argmax(env >= top - 1))  # within 1 dB of the note's peak (90 % of amplitude)
        k10 = int(np.argmax(env >= top - 20))  # 10 % of amplitude
        out["attack_ms"] = round(max(0, k90 - k10) * HOP / SR * 1000)
        out["level_db"] = round(float(np.mean(fr["level"][steady])), 1) if s1 > s0 else round(float(top), 1)
        if family in ("piano", "kit", "guitar", "guitar-steel") or n.get("decays"):
            after = env[int(np.argmax(env)):]
            if len(after) > 10:
                t = np.arange(len(after)) * HOP / SR
                out["decay_db_s"] = round(float(np.polyfit(t, after, 1)[0]), 1)
    for key, name in (("centroid", "centroid_hz"), ("hnr", "hnr_db"), ("flat", "flatness_db"), ("rough", "roughness")):
        out[name] = mean(fr[key][steady])
    on = slice(k0, min(k1, k0 + 12))
    out["onset_roughness"] = mean(fr["rough"][on])
    c = [v for v in fr["cents"][steady] if v is not None]
    if len(c) >= 5:
        med = float(np.median(c))
        out["cents"] = round(med, 1)
        # vibrato plus jitter, after removing a slow drift, over the frames within 300 cents of the note's median: a
        # vibrato is never wider, and one mistracked onset frame (a real cello C2's first frame read 100 Hz) would
        # otherwise read as 131 cents of wobble. A note with fewer than half its frames there has no steady pitch.
        k = np.array([abs(v - med) <= 300 for v in c])
        if k.sum() >= max(5, len(c) / 2):
            t = np.arange(len(c))[k]
            cc = np.array(c)[k]
            resid = cc - np.polyval(np.polyfit(t, cc, 1), t)
            out["wobble_cents"] = round(float(np.std(resid)), 1)
    f0s = [v for v in fr["f0"][steady] if v is not None]
    f0 = float(np.median(f0s)) if f0s else (n.get("target_hz") or None)
    if family not in ("piano", "kit") and not n.get("decays"):
        names = (n.get("note") or "").split("+")
        period = (common_period([hz_of(midi_of(p)) for p in names]) if len(names) > 1 and all(midi_of(p) is not None for p in names)
                  else 1 / f0 if f0 else None)
        if period:
            out["speaks_ms"], out["speaks"] = speaks_ms(x, a, d, period)
    if f0 and i1 - i0 > SR // 20:
        seg = x[int((a + lead) * SR): i1]
        if len(seg) > SR // 20:
            db, floor, stretch = partials(seg, f0)
            out["partials_db"] = db
            out["noise_floor_db"] = floor
            if stretch is not None:
                out["stretch_cents"] = stretch
    bands, lv = third_octave(x[i0:i1]) if i1 - i0 > SR // 10 else (None, None)
    if lv:
        out["bands_db"] = lv
    return out


def timbre_distance(a, b):
    """Mean absolute difference (dB) of two third-octave profiles, after aligning their overall level."""
    if not a or not b:
        return None
    a, b = np.array(a), np.array(b)
    ok = (a > -60) | (b > -60)
    if not ok.any():
        return None
    d = (a - b)[ok]
    return round(float(np.mean(np.abs(d - np.median(d)))), 1)


def envelope_correlation(fa, fb):
    a, b = np.array(fa["level"]), np.array(fb["level"])
    n = min(len(a), len(b))
    if n < 10:
        return None
    a, b = np.maximum(a[:n], -60), np.maximum(b[:n], -60)
    if a.std() == 0 or b.std() == 0:
        return None
    return round(float(np.corrcoef(a, b)[0, 1]), 3)


def summarise(notes):
    def med(key, f=lambda s, r: s - r):
        v = [f(n["sim"][key], n["ref"][key]) for n in notes
             if n.get("ref") and n["sim"].get(key) is not None and n["ref"].get(key) is not None]
        return round(float(np.median(v)), 2) if v else None
    out = dict(
        attack_ratio=med("attack_ms", lambda s, r: (s + 5) / (r + 5)),
        brightness_ratio=med("centroid_hz", lambda s, r: s / r if r else None),
        harmonicity_diff_db=med("hnr_db"),
        noise_floor_diff_db=med("noise_floor_db"),
        roughness_diff=med("roughness"),
        onset_roughness_diff=med("onset_roughness"),
        speaks_ratio=med("speaks_ms", lambda s, r: (s + 20) / (r + 20)),
        wobble_diff_cents=med("wobble_cents"),
        level_diff_db=med("level_db"),
    )
    td = [n["timbre_distance_db"] for n in notes if n.get("timbre_distance_db") is not None]
    out["timbre_distance_db"] = round(float(np.median(td)), 1) if td else None
    pd = []
    for n in notes:
        s, r = (n.get("sim") or {}).get("partials_db"), (n.get("ref") or {}).get("partials_db")
        if s and r:
            k = min(len(s), len(r), 12)
            pd.append(float(np.mean(np.abs(np.array(s[:k]) - np.array(r[:k])))))
    out["partials_diff_db"] = round(float(np.median(pd)), 1) if pd else None
    return {k: v for k, v in out.items() if v is not None}


def targets_fn(notes, key="target_hz"):
    spans = [(n["at"], n["at"] + n["dur"], n.get(key)) for n in notes]

    def f(t):
        for a, b, hz in spans:
            if a <= t < b:
                return hz
        return None
    return f


def chords_fn(notes):
    """t -> the fundamentals of a chord's notes (equal temperament from the names, 'C2+G2') while it sounds, else None."""
    spans = [(n["at"], n["at"] + n["dur"], [440 * 2 ** ((midi_of(p) - 69) / 12) for p in n["note"].split("+")])
             for n in notes if "+" in (n.get("note") or "") and all(midi_of(p) is not None for p in n["note"].split("+"))]

    def f(t):
        for a, b, hz in spans:
            if a <= t < b:
                return hz
        return None
    return f


def own_pitch(x, notes):
    """Each note's pitch measured on this track's steady part with the long-window detector (the simulation's notes
    carry theirs from phrasekit; the reference track is measured here the same way)."""
    out = []
    for n in notes:
        a, d = n["at"], n["dur"]
        seg = x[int((a + min(0.12, 0.3 * d)) * SR): int(max(a + 0.2, a + d - 0.04) * SR)]
        f = f0_samples(seg, SR)[0] if n.get("target_hz") and len(seg) > SR // 20 else None
        out.append(dict(at=a, dur=d, f0_hz=f))
    return out


def quality(x):
    peak = float(np.max(np.abs(x)) or 1.0)
    clipped = int(np.sum(np.abs(x) >= 0.999 * peak)) if peak >= 0.999 else 0
    rms = float(np.sqrt(np.mean(x ** 2)) or 1e-12)
    return dict(clipped_samples=clipped, crest_db=round(20 * math.log10(peak / rms), 1))


def main():
    doc = json.loads((PH / "phrases.json").read_text())
    bank = build_bank()
    screened = screened_clips()
    for p in doc["phrases"]:
        x = wav_mono(PH / p["audio"])
        notes = p.get("notes") or []
        tf = targets_fn(notes)
        ref_x, used, ref_how, missing = None, None, None, []
        if p.get("real") and p["real"].get("audio"):
            ref_x = wav_mono(PH / p["real"]["audio"])
            ref_how = "the real recording shown with the phrase"
        else:
            got = reference_phrase(p, bank, screened)
            if got:
                ref_x, used, missing = got
                # the same length as the simulation, so the two maps and every lane share one time axis
                ref_x = np.pad(ref_x, (0, max(0, len(x) - len(ref_x))))[: len(x)]
                write_wav(PH / f"{p['id']}-ref.wav", ref_x)
                spectrogram(ref_x, p["f_lo_hz"], p["f_hi_hz"]).save(PH / f"{p['id']}-ref.webp", quality=88, method=6)
                libs = sorted({u["licence"] for u in used})
                bowed = p["family"] in ("violin", "viola", "cello", "double-bass")
                ref_how = ("Each note is a real recording of that same note"
                           + (" (University of Iowa MIS, mezzo-forte" + (", on the same string where the library has it"
                                                                         if bowed else "") + ")"
                              if libs == ["UIowa"] else "")
                           + ", cut from its onset to the simulation note's slot and placed at the same time, "
                           "joined the same way. It is a splice of single notes, not a performance of the phrase: "
                           "listen to the notes, not the joins.")
                if any("+" in n["note"] for n in notes):
                    ref_how += " A double stop is two single-note recordings played together."
                if missing:
                    gone = sorted({notes[i]["note"] for i in missing})
                    ref_how += (f" Only {len(notes) - len(missing)} of the {len(notes)} notes have an open-licence "
                                f"recording; the others ({', '.join(gone)}) are silent in this track and are not "
                                "compared.")
                p["ref"] = dict(audio=f"{p['id']}-ref.wav", map=f"{p['id']}-ref.webp",
                                duration_s=round(len(ref_x) / SR, 3), how=ref_how,
                                missing=[notes[i]["note"] for i in missing],
                                sources=[dict(note=u["note"], source=u["source"], licence=u["licence"],
                                              credit=u["credit"]) for u in used])
        f_lo = p.get("f_lo_hz", 50.0)
        cf = chords_fn(notes)
        fs = frames(x, tf, f_lo, targets_fn(notes, "f0_hz"), cf)
        fr_ref = frames(ref_x, tf, f_lo, targets_fn(own_pitch(ref_x, notes), "f0_hz"), cf) if ref_x is not None else None
        out_notes = []
        refused = [k for k, n in enumerate(notes) if n.get("refused")]
        for k, n in enumerate(notes):
            if k in refused:  # the solver refused this note and its slot is silent: nothing to measure
                out_notes.append(dict(i=k, refused=n["refused"], sim={}))
                continue
            row = dict(i=k, sim=note_measures(x, fs, n, p["family"]))
            if n.get("whole_cycles"):  # a median of frames across a vibrato swing is not its centre; the map's reading is
                row["sim"]["cents"] = n["cents"]
            if fr_ref is not None and k not in missing and n["at"] + 0.05 < len(ref_x) / SR:
                row["ref"] = note_measures(ref_x, fr_ref, n, p["family"])
                row["timbre_distance_db"] = timbre_distance(row["sim"].get("bands_db"), row["ref"].get("bands_db"))
            out_notes.append(row)
        bands, lt_sim = third_octave(x)
        an = dict(id=p["id"], hop_s=HOP / SR, frames=dict(sim=fs), notes=out_notes,
                  ltas=dict(bands_hz=bands, sim=lt_sim), quality=dict(sim=quality(x)))
        if fr_ref is not None:
            an["frames"]["ref"] = fr_ref
            an["ltas"]["ref"] = third_octave(ref_x)[1]
            an["quality"]["ref"] = quality(ref_x)
            if not missing and not refused:  # silent gaps in either track would only measure the gaps
                an["envelope_correlation"] = envelope_correlation(fs, fr_ref)
            an["summary"] = summarise([r for r in out_notes if "refused" not in r]) if notes else {}
            if ref_how:
                an["ref_how"] = ref_how
        (PH / f"{p['id']}.analysis.json").write_text(json.dumps(an, separators=(",", ":")))
        p["analysis"] = f"{p['id']}.analysis.json"
        print(p["id"], "ref" if fr_ref is not None else "no ref", len(out_notes), "notes", flush=True)
    (PH / "phrases.json").write_text(json.dumps(doc, indent=1))


if __name__ == "__main__":
    main()
