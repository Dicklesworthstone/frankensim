"""Verisimilitude scorecard: FrankenSim renders against real recordings.

Compares offline FrankenSim clips with University of Iowa Musical Instrument
Samples (Steinway & Sons Model B piano; anechoic splash and crash cymbals; the
university states they "may be downloaded and used for any projects without
restrictions"). Every descriptor is computed identically on both signals, and
each gap is placed against a RULER MADE OF REAL RECORDINGS (real-vs-real
differences), never against an invented tolerance.

usage: uv run --with numpy --with pillow python tools/verisim.py CAMPAIGN_DIR REF_WAV48_DIR PUBLIC_DIR

REF_WAV48_DIR holds the Iowa files converted to 48 kHz float WAV
(ffmpeg -i X.aiff -ar 48000 -c:a pcm_f32le X.wav); originals are identified by
URL and SHA-256 in the output, not redistributed beyond short level-matched
excerpts for listening.
"""
from __future__ import annotations
import hashlib, json, math, re, struct, subprocess, sys
from pathlib import Path
import numpy as np
from PIL import Image

CAMPAIGN = REFS = OUT = Path(".")  # set by configure(); tools/diagnose.py imports this module
SR = 48000
IOWA_PIANO = "https://theremin.music.uiowa.edu/sound%20files/MIS/Piano_Other/piano/"
IOWA_CYM = "https://theremin.music.uiowa.edu/sound%20files/MIS/Percussion/cymbals/suspended%20cymbals/stick/"
RIGHTS = "University of Iowa Electronic Music Studios, Musical Instrument Samples: \"may be downloaded and used for any projects without restrictions\" (theremin.music.uiowa.edu/MIS.html)"


def r(v, d=5):
    return None if v is None or not np.isfinite(v) else float(f"{v:.{d}g}")


def load(path: Path, channel=0):
    """RIFF/WAVE reader for the encodings used here (16- and 24-bit PCM, 32-bit float); anything else refuses."""
    b = path.read_bytes()
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
    assert fmt is not None and data is not None, path
    tag, ch, sr, _, _, bits = struct.unpack("<HHIIHH", fmt[:16])
    if tag == 0xFFFE:  # WAVE_FORMAT_EXTENSIBLE: the subformat GUID starts with the plain tag
        tag = struct.unpack("<H", fmt[24:26])[0]
    if (tag, bits) == (1, 16):
        x = np.frombuffer(data, "<i2").astype(np.float64) / 32768
    elif (tag, bits) == (1, 24):
        b3 = np.frombuffer(data[: len(data) // 3 * 3], np.uint8).reshape(-1, 3).astype(np.int32)
        v = b3[:, 0] | (b3[:, 1] << 8) | (b3[:, 2] << 16)
        x = np.where(v >= 1 << 23, v - (1 << 24), v).astype(np.float64) / (1 << 23)
    elif (tag, bits) == (3, 32):
        x = np.frombuffer(data, "<f4").astype(np.float64)
    else:
        raise SystemExit(f"{path}: unsupported WAV encoding tag={tag} bits={bits}")
    assert sr == SR, path
    # One channel only: averaging spaced microphones can cancel partials.
    return x[: len(x) // ch * ch].reshape(-1, ch)[:, channel]


def highpass(x, f0=40.0, f1=70.0):
    """Zero-phase raised-cosine high-pass (0 below f0, 1 above f1), applied identically to both
    signals of every comparison. Room rumble below 70 Hz precedes and underlies the Iowa pp
    piano notes. The piano notes compared here start at 220 Hz; the simulated splash keeps modes
    down to about 50 Hz, so its 63 Hz octave band is measured after this filter on both sides.
    Half a second of zero padding keeps the filter's two-sided response from wrapping around."""
    n = 1 << int(np.ceil(np.log2(len(x) + SR // 2)))
    X = np.fft.rfft(x, n)
    f = np.fft.rfftfreq(n, 1 / SR)
    g = np.clip((f - f0) / (f1 - f0), 0, 1)
    X *= 0.5 - 0.5 * np.cos(np.pi * g)
    return np.fft.irfft(X, n)[: len(x)]


def onset(x):
    # Relative to BOTH the peak and the recording's own noise floor: a pp note in a room
    # sits less than 30 dB above the hiss before it, and the hiss is not the note.
    # The first 50 ms count as a noise floor only when steady (10th to 90th percentile within
    # 20 dB): a render whose note starts 10 ms in holds silence and an onset there, not hiss.
    env = np.convolve(np.abs(x), np.ones(96) / 96, "same")
    head = env[: max(int(0.05 * SR), 1)]
    lo, hi = np.percentile(head, [10, 90])
    noise = float(hi) if env.argmax() > len(head) and hi <= 10 * lo else 0.0
    thr = max(env.max() * 10 ** (-30 / 20), noise * 10 ** (15 / 20))
    return max(0, int(np.argmax(env > thr)) - int(0.005 * SR))


def segment(x, dur):
    o = onset(x)
    y = x[o:o + int(dur * SR)]
    y = np.pad(y, (0, int(dur * SR) - len(y)))
    env = np.convolve(np.abs(y), np.ones(96) / 96, "same")
    first = int(np.argmax(env > env.max() * 10 ** (-20 / 20))) / SR
    assert first < 0.012, f"onset not aligned: the note starts {first*1000:.1f} ms into the window"
    return y


def early_rms(x, t=0.3):
    return float(np.sqrt(np.mean(x[: int(t * SR)] ** 2)))


# ------------------------------------------------------------------ descriptors
def spectrum(x, t0, t1, nfft=1 << 20):
    seg = x[int(t0 * SR):int(t1 * SR)]
    s = np.abs(np.fft.rfft(seg * np.hanning(len(seg)), nfft))
    return np.fft.rfftfreq(nfft, 1 / SR), s


def partials(x, f_guess, count=16):
    f, s = spectrum(x, 0.05, 1.05)
    floor = np.median(s[(f > 50) & (f < 12000)])
    found, B, misses = [], 0.0, 0
    f0 = f_guess
    for n in range(1, count + 1):
        pred = n * f0 * math.sqrt(1 + B * n * n)
        half = 0.035 / max(1, n ** 0.3)
        sel = (f > pred * (1 - half)) & (f < pred * (1 + half))
        if not sel.any():
            break
        k = np.argmax(np.where(sel, s, 0))
        if s[k] < 8 * floor:  # below the noise floor: not a partial, never fit it
            misses += 1
            if misses >= 3:
                break
            continue
        misses = 0
        if 0 < k < len(s) - 1:
            a, b, c = np.log(s[k - 1:k + 2] + 1e-30)
            fk = f[k] + 0.5 * (a - c) / (a - 2 * b + c) * (f[1] - f[0])
        else:
            fk = f[k]
        found.append((n, fk, s[k]))
        if len(found) >= 2:  # refit f0 and B from (f_n/n)^2 = f0^2 + f0^2 B n^2
            ns = np.array([q[0] for q in found], float)
            ys = (np.array([q[1] for q in found]) / ns) ** 2
            c1, c0 = np.polyfit(ns ** 2, ys, 1)
            f0, B = math.sqrt(max(c0, 1e-9)), max(0.0, c1 / max(c0, 1e-9))
    # Every level is relative to partial 1, so a fit that missed it is refused, not relabelled.
    assert found and found[0][0] == 1, f"partial 1 not found near {f_guess} Hz"
    return {"f1": found[0][1], "B": B, "f0": f0, "count": len(found),
            "partials": [{"n": n, "f": r(fk, 7), "db": r(20 * math.log10(a / found[0][2]), 4)} for n, fk, a in found]}


def partial_decay(x, found, t0=0.2, t1=2.6):
    """Linear dB/s fit per partial, keyed by the partial's NUMBER (a skipped partial stays skipped)."""
    n, hop = 8192, 1200
    win = np.hanning(n)
    t1 = min(t1, len(x) / SR)
    frames = [x[i:i + n] * win for i in range(int(t0 * SR), int(t1 * SR) - n, hop)]
    if len(frames) < 6:
        return {}, t1
    S = np.abs(np.fft.rfft(np.array(frames), axis=1))
    bins = np.fft.rfftfreq(n, 1 / SR)
    t = np.array([t0 + (k * hop + n / 2) / SR for k in range(len(frames))])
    out = {}
    for p in found:
        k = int(np.argmin(np.abs(bins - p["f"])))
        a = 20 * np.log10(S[:, max(0, k - 2):k + 3].max(axis=1) + 1e-12)
        out[p["n"]] = r(float(np.polyfit(t, a, 1)[0]), 4)
    return out, t1


def centroid_track(x, dur):
    n, hop = 2048, 480
    win = np.hanning(n)
    f = np.fft.rfftfreq(n, 1 / SR)
    tt, cc = [], []
    for i in range(0, int(dur * SR) - n, hop):
        S = np.abs(np.fft.rfft(x[i:i + n] * win)) ** 2
        tot = S.sum()
        if tot <= 0:
            continue
        tt.append(r((i + n / 2) / SR, 4))
        cc.append(r(float((f * S).sum() / tot), 4))
    return {"t": tt[::2], "hz": cc[::2]}


def low_knock(x, f1, f_guess):
    """Strongest component between 80 Hz and 0.77*f_guess (the same band on both signals),
    relative to this signal's partial 1, early and late."""
    out = {}
    for tag, t0, t1 in (("early", 0.05, 0.5), ("late", 1.5, 2.9)):
        if t1 * SR > len(x):
            continue
        f, s = spectrum(x, t0, t1, 1 << 19)
        lo = (f > 80) & (f < 0.77 * f_guess)
        p1 = (f > f1 * 0.96) & (f < f1 * 1.04)
        if not lo.any() or not p1.any():
            continue
        k, k1 = np.argmax(np.where(lo, s, 0)), np.argmax(np.where(p1, s, 0))
        out[tag] = {"hz": r(f[k], 5), "db": r(20 * math.log10(s[k] / s[k1]), 3)}
    return out


def band_ratio(x, t0=0.0, t1=1.0):
    f, s = spectrum(x, t0, t1, 1 << 16)
    p = s ** 2
    hi = p[(f >= 2000) & (f < 8000)].sum()
    lo = p[(f >= 200) & (f < 2000)].sum()
    return 10 * math.log10(hi / lo) if lo > 0 and hi > 0 else None


def rise_ms(x):
    # The hammer/stick transient only: a room's slow build-up is not the attack.
    env = np.convolve(np.abs(x), np.ones(48) / 48, "same")[: int(0.06 * SR)]
    pk = env.max()
    a = int(np.argmax(env >= 0.1 * pk))
    b = int(np.argmax(env >= 0.9 * pk))
    return (b - a) / SR * 1000


OCT = [63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000]


def octave_levels(x, t0, t1):
    f, s = spectrum(x, t0, t1, 1 << 16)
    p = s ** 2
    out = []
    for fc in OCT:
        sel = (f >= fc / math.sqrt(2)) & (f < fc * math.sqrt(2))
        out.append(10 * math.log10(p[sel].sum() + 1e-30))
    tot = 10 * math.log10(sum(10 ** (v / 10) for v in out))
    return [v - tot for v in out]


def octave_decay(x, t0=0.05, t1=1.5):
    n, hop = 4096, 960
    win = np.hanning(n)
    frames = np.array([x[i:i + n] * win for i in range(int(t0 * SR), int(t1 * SR) - n, hop)])
    S = np.abs(np.fft.rfft(frames, axis=1)) ** 2
    f = np.fft.rfftfreq(n, 1 / SR)
    t = np.array([t0 + (k * hop + n / 2) / SR for k in range(len(frames))])
    out = []
    for fc in OCT:
        sel = (f >= fc / math.sqrt(2)) & (f < fc * math.sqrt(2))
        e = 10 * np.log10(S[:, sel].sum(axis=1) + 1e-30)
        out.append(float(np.polyfit(t, e, 1)[0]) if np.isfinite(e).all() else float("nan"))
    return out


def rms_db(a, b):
    a, b = np.array(a, float), np.array(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.sqrt(np.mean((a[ok] - b[ok]) ** 2))) if ok.any() else None


# --------------------------------------------------------------------- images
STOPS = [(0.00, (8, 9, 14)), (0.18, (32, 22, 48)), (0.36, (92, 30, 58)), (0.55, (170, 64, 38)),
         (0.72, (222, 128, 42)), (0.88, (244, 200, 112)), (1.00, (255, 246, 222))]
DIV = [(0.0, (58, 120, 190)), (0.42, (24, 40, 60)), (0.5, (14, 15, 18)), (0.58, (70, 34, 20)), (1.0, (232, 106, 60))]


def lut(stops):
    t = np.zeros((256, 3), np.uint8)
    for k in range(256):
        u = k / 255
        for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
            if t0 <= u <= t1:
                w = (u - t0) / (t1 - t0)
                t[k] = [round(a + (b - a) * w) for a, b in zip(c0, c1)]
                break
    return t


LUT, DLUT = lut(STOPS), lut(DIV)
F_LO, F_HI, ROWS = 40.0, 16000.0, 260


def spec_db(x, cols=700):
    n = 4096
    hop = max(1, (len(x) - n) // (cols - 1))
    win = np.hanning(n)
    pad = np.concatenate([np.zeros(n // 2), x, np.zeros(n)])
    frames = np.stack([pad[i * hop:i * hop + n] * win for i in range(cols)])
    mag = np.abs(np.fft.rfft(frames, axis=1))
    f = np.fft.rfftfreq(n, 1 / SR)
    edges = np.geomspace(F_LO, F_HI, ROWS + 1)
    out = np.empty((ROWS, cols))
    for i in range(ROWS):
        sel = (f >= edges[i]) & (f < edges[i + 1])
        out[i] = mag[:, sel].max(axis=1) if sel.any() else np.array([np.interp(math.sqrt(edges[i] * edges[i + 1]), f, m) for m in mag])
    return 20 * np.log10(out + 1e-12)


def save(img, name):
    Image.fromarray(img, "RGB").save(OUT / "img" / f"{name}.webp", "WEBP", quality=88, method=6)
    return f"img/{name}.webp"


def images(real, ours, pid):
    a, b = spec_db(real), spec_db(ours)
    ref = max(a.max(), b.max())
    span = 80.0
    to_img = lambda d: LUT[(np.clip((d - (ref - span)) / span, 0, 1)[::-1] * 255).astype(np.uint8)]
    # Difference only where either signal is above the shared floor: silence vs silence is not a match.
    floor = ref - 60
    diff = np.clip(b - a, -30, 30)
    mask = (np.maximum(a, b) > floor)
    d = np.where(mask, (diff + 30) / 60, 0.5)
    return {"real": save(to_img(a), f"vs-{pid}-real"), "ours": save(to_img(b), f"vs-{pid}-ours"),
            "diff": save(DLUT[(d[::-1] * 255).astype(np.uint8)], f"vs-{pid}-diff"),
            "f_lo": F_LO, "f_hi": F_HI, "span_db": span, "diff_range_db": 30}


def encode(x, name):
    tmp = OUT / "audio" / f".{name}.wav"
    y = x.astype("<f4")
    data = y.tobytes()
    tmp.write_bytes(b"RIFF" + struct.pack("<I", 36 + len(data)) + b"WAVEfmt " +
                    struct.pack("<IHHIIHH", 16, 3, 1, SR, SR * 4, 4, 32) + b"data" + struct.pack("<I", len(data)) + data)
    dst = OUT / "audio" / f"{name}.m4a"
    for codec in (["-c:a", "aac_at", "-b:a", "224k"], ["-c:a", "aac", "-b:a", "224k"]):
        if subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", str(tmp), *codec, "-movflags", "+faststart", str(dst)]).returncode == 0:
            break
    else:
        raise SystemExit(f"ffmpeg could not encode {dst}")
    tmp.unlink()
    return f"audio/{name}.m4a"


def sha(p: Path):
    return hashlib.sha256(p.read_bytes()).hexdigest()


# ---------------------------------------------------------------------- pairs
def ours_clip(kind, cid):
    if kind == "piano":
        return load(CAMPAIGN / "piano" / f"{cid}.wav")
    return load(CAMPAIGN / "perc" / f"{cid}.wav")


def piano_label(cid):
    """What the render simulated, read from its own command line."""
    log = CAMPAIGN / "piano" / f"{cid}.log"
    text = log.read_text(errors="replace") if log.exists() else ""
    cmd = next((l for l in text.splitlines() if l.startswith("cmd:")), "")
    extra = []
    # A stretching file can select linear strings for every key; the run reports how many channels it stretched.
    stretched = re.search(r"String extension: .*?; (\d+) geometric speaking/duplex channels; (\d+) retained longitudinal", text)
    if "--string-stretching" in cmd and stretched and int(stretched.group(1)) > 0:
        extra.append("longitudinal string waves whose axial force drives the bridge" if int(stretched.group(2)) > 0
                     else "string stretching (tension modulation)")
    if "--string-polarization" in cmd:
        extra.append("both transverse string polarizations")
    return ("FrankenSim offline render (RT-0425 Model D scale; " + "".join(e + "; " for e in extra)
            + "one pressure point 1 m above the board, infinite baffle, no room)")


def piano_pair(pid, title, ours_id, real_name, f_guess, dur):
    real_raw = highpass(load(REFS / f"{real_name}.wav"))
    ours_raw = highpass(ours_clip("piano", ours_id))
    real, ours = segment(real_raw, dur), segment(ours_raw, dur)
    g_r, g_o = 1 / early_rms(real), 1 / early_rms(ours)
    real_n, ours_n = real * g_r, ours * g_o
    pk = max(np.abs(real_n).max(), np.abs(ours_n).max())
    play = 10 ** (-1 / 20) / pk
    pr, po = partials(real, f_guess), partials(ours, f_guess)
    kr, ko = low_knock(real, pr["f1"], f_guess), low_knock(ours, po["f1"], f_guess)
    dr, dt1 = partial_decay(real, [p for p in pr["partials"] if p["n"] <= 6])
    do, _ = partial_decay(ours, [p for p in po["partials"] if p["n"] <= 6])
    decay_n = list(range(1, 7))
    env_r = {p["n"]: p["db"] for p in pr["partials"]}
    env_o = {p["n"]: p["db"] for p in po["partials"]}
    common = sorted(set(env_r) & set(env_o) - {1})
    return {
        "id": pid, "kind": "piano", "title": title, "duration_s": dur,
        "real": {"label": f"Steinway & Sons Model B · {real_name}", "url": IOWA_PIANO + real_name + ".aiff",
                 "sha256_wav48": sha(REFS / f"{real_name}.wav"), "rights": RIGHTS,
                 "recording": "two Neumann KM 84s 8 in above the bass and treble strings (the bass-side channel is used), 16-bit 44.1 kHz, in a normal room"},
        "ours": {"clip": ours_id, "label": piano_label(ours_id)},
        "audio": {"real": encode(real_n * play, f"vs-{pid}-real"), "ours": encode(ours_n * play, f"vs-{pid}-ours"),
                  "level_match": "both scaled to equal RMS over the first 300 ms after onset; absolute levels are not comparable (different instruments, microphones and distances)"},
        "images": images(real_n, ours_n, pid),
        "descriptors": {
            "f1_hz": {"real": r(pr["f1"], 7), "ours": r(po["f1"], 7), "note": "first partial; the Iowa piano's tuning is not stated"},
            "inharmonicity_B": {"real": r(pr["B"], 3), "ours": r(po["B"], 3), "partials_used": [pr["count"], po["count"]], "note": "fitted from (f_n/n)^2 = f0^2(1+Bn^2) over partials above the noise floor"},
            "partial_envelope_db": {"real": [env_r[n] for n in common], "ours": [env_o[n] for n in common], "n": common,
                                    "rms_gap_db": r(rms_db([env_r[n] for n in common], [env_o[n] for n in common]), 3),
                                    "note": "partial level relative to partial 1, 0.05-1.05 s after onset"},
            "partial_decay_db_s": {"n": decay_n, "real": [dr.get(n) for n in decay_n], "ours": [do.get(n) for n in decay_n],
                                   "note": f"linear fit of each partial's level 0.2-{dt1:.1f} s after onset; a dash is a partial below the noise floor"},
            "high_band_db": {"real": r(band_ratio(real), 4), "ours": r(band_ratio(ours), 4), "note": "energy 2-8 kHz relative to 200-2000 Hz, first second"},
            "attack_rise_ms": {"real": r(rise_ms(real), 3), "ours": r(rise_ms(ours), 3), "note": "10% to 90% of the envelope peak"},
            **({"low_knock_early_db": {"real": kr["early"]["db"], "ours": ko["early"]["db"],
                                       "hz": {"real": kr["early"]["hz"], "ours": ko["early"]["hz"]},
                                       "note": f"strongest component 80 Hz-{0.77 * f_guess:.0f} Hz in the first 0.05-0.5 s, relative to partial 1 (real at {kr['early']['hz']} Hz, ours at {ko['early']['hz']} Hz)"}}
               if "early" in kr and "early" in ko else {}),
            **({"low_knock_late_db": {"real": kr["late"]["db"], "ours": ko["late"]["db"],
                                      "hz": {"real": kr["late"]["hz"], "ours": ko["late"]["hz"]},
                                      "note": f"the same, 1.5-2.9 s after onset (real at {kr['late']['hz']} Hz, ours at {ko['late']['hz']} Hz)"}}
               if "late" in kr and "late" in ko else {}),
        },
        "centroid": {"real": centroid_track(real_n, dur), "ours": centroid_track(ours_n, dur)},
        "_env": (env_r, env_o),
        "_level": (float(np.abs(real).max()), float(np.abs(ours).max())),
    }


def cymbal_desc(x):
    return {"early": octave_levels(x, 0.0, 0.25), "decay": octave_decay(x)}


def kit_grid_mm(cid):
    """The grid spacing a kit render reports in its own log ("splash 8: N nodes at X mm")."""
    log = CAMPAIGN / "perc" / f"{cid}.log"
    m = re.search(r"nodes at ([\d.]+) mm", log.read_text(errors="replace")) if log.exists() else None
    if not m:
        raise SystemExit(f"{log}: no grid spacing; render the splash with the kit example and keep its log")
    return m.group(1)


def cymbal_pair(pid, title, ours_id, real_name, dur):
    real_raw = highpass(load(REFS / f"{real_name}.wav"))
    ours_raw = highpass(ours_clip("perc", ours_id))
    real, ours = segment(real_raw, dur), segment(ours_raw, dur)
    real_n, ours_n = real / early_rms(real), ours / early_rms(ours)
    pk = max(np.abs(real_n).max(), np.abs(ours_n).max())
    play = 10 ** (-1 / 20) / pk
    cr, co = cymbal_desc(real), cymbal_desc(ours)
    hi = lambda x: band_ratio(x, 0.0, 0.5)
    return {
        "id": pid, "kind": "cymbal", "title": title, "duration_s": dur,
        "real": {"label": f"Iowa splash cymbal · {real_name}", "url": IOWA_CYM + real_name + ".aif",
                 "sha256_wav48": sha(REFS / f"{real_name}.wav"), "rights": RIGHTS,
                 "recording": "anechoic chamber, Earthworks QTC40 stereo at 5 ft (first channel used), 24-bit 44.1 kHz"},
        "ours": {"clip": ours_id, "label": f"FrankenSim offline render (fs-vkplate, the drum kit's solver: an 8-inch splash as a von Kármán shell on a {kit_grid_mm(ours_id)} mm grid; one point about 1 m above it, no room)"},
        "audio": {"real": encode(real_n * play, f"vs-{pid}-real"), "ours": encode(ours_n * play, f"vs-{pid}-ours"),
                  "level_match": "both scaled to equal RMS over the first 300 ms after onset"},
        "images": images(real_n, ours_n, pid),
        "descriptors": {
            "octave_early_db": {"bands_hz": OCT, "real": [r(v, 4) for v in cr["early"]], "ours": [r(v, 4) for v in co["early"]],
                                "rms_gap_db": r(rms_db(cr["early"], co["early"]), 3), "note": "octave-band level relative to total, first 250 ms"},
            "octave_decay_db_s": {"bands_hz": OCT, "real": [r(v, 4) for v in cr["decay"]], "ours": [r(v, 4) for v in co["decay"]],
                                  "note": "linear fit of each octave band 0.05-1.5 s after onset"},
            "high_band_db": {"real": r(hi(real), 4), "ours": r(hi(ours), 4), "note": "energy 2-8 kHz relative to 200-2000 Hz, first 0.5 s"},
            "attack_rise_ms": {"real": r(rise_ms(real), 3), "ours": r(rise_ms(ours), 3), "note": "10% to 90% of the envelope peak"},
        },
        "centroid": {"real": centroid_track(real_n, dur), "ours": centroid_track(ours_n, dur)},
        "_early": (cr["early"], co["early"]),
    }


def configure(campaign, refs, out):
    global CAMPAIGN, REFS, OUT
    CAMPAIGN, REFS, OUT = Path(campaign), Path(refs), Path(out)


def main():
    (OUT / "img").mkdir(parents=True, exist_ok=True)
    (OUT / "audio").mkdir(parents=True, exist_ok=True)
    pairs = []
    ladder = {}
    for dyn, speed in [("pp", "0.6"), ("mf", "2.0"), ("ff", "4.5")]:
        p = piano_pair(f"a3-{dyn}", f"A3 {dyn} · {speed} m/s: real Steinway vs FrankenSim at {speed} m/s", f"a3-{dyn}", f"Piano.{dyn}.A3", 220.0, 3.0)
        pairs.append(p)
        ladder[dyn] = p
    if (CAMPAIGN / "piano" / "a3-fff.wav").exists():
        p8 = piano_pair("a3-ff8", "A3 ff · 8 m/s: real Steinway ff vs FrankenSim at 8 m/s, the hardest stroke rendered", "a3-fff", "Piano.ff.A3", 220.0, 3.0)
        pairs.append(p8)
        ladder["ff8"] = p8
    pairs.append(piano_pair("c6-mf", "C6 mf · 2.0 m/s: real Steinway vs FrankenSim at 2.0 m/s", "c6-mf", "Piano.mf.C6", 1046.5, 2.0))

    # Piano ruler: how far apart are two DYNAMIC LAYERS of the same real piano?
    def env_gap(a, b):
        common = sorted(set(a) & set(b) - {1})
        return rms_db([a[n] for n in common], [b[n] for n in common])
    piano_ruler = {
        "real_pp_vs_real_mf_db": r(env_gap(ladder["pp"]["_env"][0], ladder["mf"]["_env"][0]), 3),
        "real_mf_vs_real_ff_db": r(env_gap(ladder["mf"]["_env"][0], ladder["ff"]["_env"][0]), 3),
        "note": "RMS difference of partial envelopes between two dynamic layers of the SAME real piano: what a change of dynamic alone does",
    }
    dyn_levels = {
        "real_db": {d: r(20 * math.log10(ladder[d]["_level"][0] / ladder["mf"]["_level"][0]), 3) for d in ladder},
        "ours_db": {d: r(20 * math.log10(ladder[d]["_level"][1] / ladder["mf"]["_level"][1]), 3) for d in ladder},
        "real_high_band_db": {d: ladder[d]["descriptors"]["high_band_db"]["real"] for d in ladder},
        "ours_high_band_db": {d: ladder[d]["descriptors"]["high_band_db"]["ours"] for d in ladder},
        "note": "peak level relative to mf, and 2-8 kHz brightness, within each source. Assumes one recording gain across the Iowa dynamic layers (their files are not level-normalized).",
    }

    cym = []
    for dyn, ours_id, speed in [("pp", "vk-splash-soft", "0.6"), ("mf", "vk-splash", "2")]:
        c = cymbal_pair(f"splash-{dyn}", f"Splash {dyn} · {speed} m/s: real Iowa splash vs FrankenSim with a {speed} m/s stick", ours_id, f"splash.stick.normal.{dyn}", 2.5)
        pairs.append(c)
        cym.append(c)

    # Cymbal ruler: distance of each REAL crash, and of FrankenSim, from the real mf splash.
    splash_mf = segment(highpass(load(REFS / "splash.stick.normal.mf.wav")), 2.5)
    ref_early = octave_levels(splash_mf, 0.0, 0.25)
    ruler = []
    for name, label in [("13crash.stick.normal.mf", "real 13-inch crash"), ("17crash.stick.normal.mf", "real 17-inch crash"),
                        ("18crash.stick.normal.mf", "real 18-inch crash"), ("splash.stick.normal.ff", "same real splash, ff"),
                        ("splash.stick.normal.pp", "same real splash, pp")]:
        x = segment(highpass(load(REFS / f"{name}.wav")), 2.5)
        ruler.append({"label": label, "rms_gap_db": r(rms_db(octave_levels(x, 0.0, 0.25), ref_early), 3), "real": True})
    for cid, label in [("vk-splash", "FrankenSim splash, 2 m/s"), ("vk-splash-soft", "FrankenSim splash, 0.6 m/s")]:
        x = segment(highpass(ours_clip("perc", cid)), 2.5)
        ruler.append({"label": label, "rms_gap_db": r(rms_db(octave_levels(x, 0.0, 0.25), ref_early), 3), "real": False})
    # The same 2 m/s stroke on other grids (vk-splash-grid-<spacing>mm.wav), when the campaign has them.
    for path in sorted((CAMPAIGN / "perc").glob("vk-splash-grid-*mm.wav"), key=lambda q: -float(q.stem.split("-")[-1][:-2])):
        x = segment(highpass(load(path)), 2.5)
        ruler.append({"label": f"FrankenSim splash, 2 m/s, {path.stem.split('-')[-1][:-2]} mm grid", "grid_mm": float(path.stem.split("-")[-1][:-2]),
                      "rms_gap_db": r(rms_db(octave_levels(x, 0.0, 0.25), ref_early), 3), "real": False})

    for p in pairs:
        for k in [k for k in p if k.startswith("_")]:
            del p[k]
    (OUT / "data" / "verisim.json").write_text(json.dumps({
        "pairs": pairs, "piano_ruler": piano_ruler, "piano_dynamics": dyn_levels,
        "cymbal_ruler": {"reference": "real Iowa splash, stick, mf", "rows": ruler,
                         "note": "RMS difference of octave-band levels (63 Hz-16 kHz, relative to total, first 250 ms) from the real mf splash"},
        "sources": {"iowa_piano": "https://theremin.music.uiowa.edu/MISpiano.html", "iowa_cymbals": "https://theremin.music.uiowa.edu/MIScymbals.html", "rights": RIGHTS},
    }, separators=(",", ":")))
    print(json.dumps({"piano_ruler": piano_ruler, "dyn": dyn_levels, "cym_ruler": ruler}, indent=1))
    for p in pairs:
        d = p["descriptors"]
        print(p["id"], {k: (v.get("real"), v.get("ours"), v.get("rms_gap_db")) if isinstance(v, dict) and not isinstance(v.get("real"), list) else ("gap", v.get("rms_gap_db")) for k, v in d.items()})


if __name__ == "__main__":
    configure(*sys.argv[1:4])
    main()
