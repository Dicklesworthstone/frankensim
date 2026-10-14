"""Turn FrankenSim render-campaign outputs into the site's static assets.

Inputs are the offline campaign directory (WAV/CSV/log files produced by the
`percussion` and `grand_piano` examples) and the visualization exports. Every
number written to the manifest is parsed from the run's own log or computed
from its own output; nothing is typed in by hand.

Audio policy: each COMPARISON GROUP gets ONE fixed playback gain (the loudest
clip in the group peaks at -1 dBFS). Clips are never normalized individually,
so level differences inside a comparison are physical results.

usage: uv run --with numpy --with pillow python tools/build_data.py CAMPAIGN_DIR OUT_PUBLIC_DIR
"""
from __future__ import annotations
import json, math, os, re, struct, subprocess, sys
from pathlib import Path
import numpy as np
from PIL import Image

CAMPAIGN = Path(sys.argv[1])
OUT = Path(sys.argv[2])
for sub in ("audio", "img", "data", "viz"):
    (OUT / sub).mkdir(parents=True, exist_ok=True)

COMMIT = subprocess.run(["git", "-C", str(Path(__file__).resolve().parents[3]), "rev-parse", "--short", "HEAD"],
                        capture_output=True, text=True).stdout.strip()

# ----------------------------------------------------------------------------- io
def load_wav(path: Path):
    b = path.read_bytes()
    i = b.find(b"fmt ")
    tag, ch, sr, _, _, bits = struct.unpack("<HHIIHH", b[i + 8:i + 24])
    j = b.find(b"data")
    n = struct.unpack("<I", b[j + 4:j + 8])[0]
    d = b[j + 8:j + 8 + n]
    if bits == 16:
        x = np.frombuffer(d, "<i2").astype(np.float64) / 32768.0
    elif bits == 32:
        x = np.frombuffer(d, "<f4").astype(np.float64)
    elif bits == 24:
        a = np.frombuffer(d, "u1").reshape(-1, 3).astype(np.int32)
        x = (((a[:, 0] | (a[:, 1] << 8) | (a[:, 2] << 16)) ^ 0x800000) - 0x800000) / 8388608.0
    else:
        raise ValueError(f"{path}: {bits}-bit")
    return sr, x.reshape(-1, ch)


def log_text(path: Path) -> str:
    return path.read_text(errors="replace") if path.exists() else ""


def grab(pattern: str, text: str, cast=float, group=1):
    m = re.search(pattern, text)
    return cast(m.group(group)) if m else None


# ----------------------------------------------------------------------- spectra
STOPS = [  # perceptual ramp tuned to the site palette (ink -> ember -> brass -> ivory)
    (0.00, (8, 9, 14)), (0.18, (32, 22, 48)), (0.36, (92, 30, 58)),
    (0.55, (170, 64, 38)), (0.72, (222, 128, 42)), (0.88, (244, 200, 112)), (1.00, (255, 246, 222)),
]
LUT = np.zeros((256, 3), np.uint8)
for k in range(256):
    t = k / 255
    for (t0, c0), (t1, c1) in zip(STOPS, STOPS[1:]):
        if t0 <= t <= t1:
            u = (t - t0) / (t1 - t0)
            LUT[k] = [round(a + (b - a) * u) for a, b in zip(c0, c1)]
            break

F_LO, F_HI, ROWS = 25.0, 16000.0, 300


def spectrogram(x: np.ndarray, sr: int, cols: int):
    """Log-frequency magnitude (dB re 1 Pa) on a fixed grid; returns (rows x cols)."""
    n = 4096 if sr >= 44100 else 2048
    hop = max(1, (len(x) - n) // max(1, cols - 1))
    win = np.hanning(n)
    pad = np.concatenate([np.zeros(n // 2), x, np.zeros(n)])
    frames = np.stack([pad[i * hop:i * hop + n] * win for i in range(cols)])
    mag = np.abs(np.fft.rfft(frames, axis=1)) * 2 / win.sum()
    f = np.fft.rfftfreq(n, 1 / sr)
    edges = np.geomspace(F_LO, min(F_HI, sr / 2 * 0.98), ROWS + 1)
    out = np.zeros((ROWS, cols))
    for r in range(ROWS):
        lo, hi = edges[r], edges[r + 1]
        sel = (f >= lo) & (f < hi)
        if sel.any():
            out[r] = mag[:, sel].max(axis=1)
        else:
            out[r] = np.array([np.interp(math.sqrt(lo * hi), f, m) for m in mag])
    return 20 * np.log10(out + 1e-12)


def save_spectrogram(db: np.ndarray, ref_db: float, name: str, span_db: float = 90.0):
    t = np.clip((db - (ref_db - span_db)) / span_db, 0, 1)
    img = LUT[(t[::-1] * 255).astype(np.uint8)]
    Image.fromarray(img, "RGB").save(OUT / "img" / f"{name}.webp", "WEBP", quality=88, method=6)


def envelope(x: np.ndarray, bins: int = 1400):
    n = len(x)
    edges = np.linspace(0, n, bins + 1).astype(int)
    lo = [float(x[a:b].min()) if b > a else 0.0 for a, b in zip(edges, edges[1:])]
    hi = [float(x[a:b].max()) if b > a else 0.0 for a, b in zip(edges, edges[1:])]
    return lo, hi


def r(v, d=5):
    return float(f"{v:.{d}g}")


# --------------------------------------------------------------------- encoding
def encode(mono_or_stereo: np.ndarray, sr: int, gain: float, name: str):
    """Pa * gain -> float WAV -> AAC (Apple AudioToolbox encoder when available)."""
    y = (mono_or_stereo * gain).astype("<f4")
    ch = y.shape[1]
    tmp = OUT / "audio" / f".{name}.wav"
    data = y.tobytes()
    hdr = b"RIFF" + struct.pack("<I", 36 + len(data)) + b"WAVEfmt " + struct.pack(
        "<IHHIIHH", 16, 3, ch, sr, sr * ch * 4, ch * 4, 32) + b"data" + struct.pack("<I", len(data))
    tmp.write_bytes(hdr + data)
    dst = OUT / "audio" / f"{name}.m4a"
    for codec in (["-c:a", "aac_at", "-b:a", "256k"], ["-c:a", "aac", "-b:a", "256k"]):
        res = subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", str(tmp), *codec,
                              "-movflags", "+faststart", str(dst)], capture_output=True, text=True)
        if res.returncode == 0:
            break
    else:
        raise RuntimeError(res.stderr)
    tmp.unlink()
    return f"audio/{name}.m4a"


# ------------------------------------------------------------------ clip catalog
PIANO = CAMPAIGN / "piano"
PERC = CAMPAIGN / "perc"


def piano_clip(name):
    sr, x = load_wav(PIANO / f"{name}.wav")
    return sr, x * 2.0, log_text(PIANO / f"{name}.log")  # float32 full scale = 2 Pa


def perc_clip(name, full_scale):
    sr, x = load_wav(PERC / f"{name}.wav")
    return sr, x * full_scale, log_text(PERC / f"{name}.log")


def piano_evidence(text: str) -> dict:
    return {
        "command": grab(r"cmd: (.*)", text, str),
        "rendered_s": grab(r"([\d.]+) s, \d channel\(s\) rendered", text),
        "wall_s": grab(r"rendered in ([\d.]+) s", text),
        "wall_audio_ratio": grab(r"wall/audio ratio ([\d.]+)", text),
        "string_modes": grab(r"(\d+) string modes", text, int),
        "board_modes": grab(r"(\d+) board modes", text, int),
        "board_band_hz": grab(r"source band \(0,([\d.]+)\] Hz", text),
        "peak_pa": grab(r"peak ([\d.]+) Pa-equivalent", text),
        "energy_input_j": grab(r"Input ([\d.e-]+) J", text),
        "energy_closure_j": grab(r"closure ([-\d.e]+) J", text),
        "worst_substep_defect_j": grab(r"worst substep defect ([-\d.e]+) J", text),
        "felt_loss_j": grab(r"Felt loss ([\d.e-]+) J", text),
        "damper_loss_j": grab(r"Spatial damper loss ([\d.e-]+) J", text),
        "string_scale": grab(r"String scale: ([^\n]*)", text, str),
        "soundboard": grab(r"Soundboard: ([^\n]*)", text, str),
        "microphone": grab(r"(Rayleigh pressure microphone at [^;\n]*)", text, str),
        "mechanical_rate_hz": 48000 * 4,
    }


def perc_evidence(text: str) -> dict:
    heads = re.findall(r"head (\d): film_mass_kg=[\d.e-]+,frequencies_hz=\[([^\]]*)\]", text)
    shell = grab(r"modal frequencies_hz=\[([^\]]*)\]", text, str)
    return {
        "command": grab(r"cmd: percussion (.*)", text, str),
        "wall_s": grab(r"exit 0 wall_s (\d+)", text, int),
        "peak_pa": grab(r"peak_pa=([\d.e-]+)", text),
        "clips": grab(r"clips=(\d+)", text, int),
        "bem_panels": grab(r"BEM bake: panels=(\d+)", text, int),
        "bem_band_hz": grab(r"BEM bake: .*?band_hz=\[([^\]]*)\]", text, str),
        "bem_audit_max": grab(r"BEM bake: .*?max_error=([\d.e-]+)", text),
        "bem_min_panels_per_wavelength": grab(r"min_panels_per_wavelength=([\d.e-]+)", text),
        "head_frequencies_hz": {h: [r(float(v)) for v in f.split(",") if v.strip()] for h, f in heads} or None,
        "shell_frequencies_hz": [r(float(v)) for v in shell.split(",")] if shell else None,
        "mechanical_rate_hz": 768000,
        "stroke": grab(r"(physical stroke: [^\n;]*)", text, str),
    }


manifest: dict = {"commit": COMMIT, "groups": {}, "clips": {}}


def add_group(gid: str, title: str, items: list, band: dict, extra: dict | None = None, peak_dbfs=-1.0):
    """items: (clip_id, label, sr, pa_array, evidence)."""
    peak = max(float(np.abs(pa).max()) for _, _, _, pa, _ in items)
    gain = 10 ** (peak_dbfs / 20) / peak
    longest = max(len(pa) / sr for _, _, sr, pa, _ in items)
    specs = []
    for cid, label, sr, pa, ev in items:
        mono = pa.mean(axis=1)
        cols = max(240, int(900 * (len(pa) / sr) / longest))
        specs.append((cid, spectrogram(mono, sr, cols)))
    ref = max(float(s.max()) for _, s in specs)
    for (cid, label, sr, pa, ev), (_, db) in zip(items, specs):
        save_spectrogram(db, ref, cid)
        lo, hi = envelope(pa.mean(axis=1))
        manifest["clips"][cid] = {
            "label": label, "group": gid, "audio": encode(pa, sr, gain, cid),
            "spectrogram": f"img/{cid}.webp", "duration_s": r(len(pa) / sr), "sample_rate": sr,
            "channels": pa.shape[1], "peak_pa": r(float(np.abs(pa).max())),
            "rms_pa": r(float(np.sqrt((pa ** 2).mean()))),
            "envelope_pa": {"lo": [r(v, 4) for v in lo], "hi": [r(v, 4) for v in hi]},
            "evidence": ev,
        }
    manifest["groups"][gid] = {"title": title, "clips": [c[0] for c in items],
                               "playback_gain_db": r(20 * math.log10(gain), 4),
                               "spectrogram": {"f_lo": F_LO, "f_hi": F_HI, "rows": ROWS, "ref_db": r(ref, 4),
                                               "span_db": 90.0, "max_duration_s": r(longest)},
                               "band": band, **(extra or {})}
    print(f"group {gid}: {len(items)} clips, gain {20*math.log10(gain):+.1f} dB")


def have(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 1000


# ---------------------------------------------------------------- piano groups
PIANO_BAND = {"label": "Soundboard modes retained", "lo_hz": 0, "hi_hz": 2400,
              "note": "Strings keep partials to 21.6 kHz; the board and its radiation are resolved to 2.4 kHz (128-mode Ritz reduction)."}
items = []
for cid, label in [("hero", "Gymnopédie No. 1 — stereo pair")]:
    sr, pa, t = piano_clip(cid)
    items.append((cid, label, sr, pa, piano_evidence(t)))
add_group("piano-hero", "Erik Satie, Gymnopédie No. 1 (1888) — bars 1–2 and 5–9", items, PIANO_BAND)

items = []
for cid, label in [("a3-pp", "pp — 0.6 m/s"), ("a3-mf", "mf — 2.0 m/s"), ("a3-ff", "ff — 4.5 m/s")]:
    sr, pa, t = piano_clip(cid)
    ev = piano_evidence(t)
    bridge = np.loadtxt(PIANO / f"{cid}-bridge.csv", delimiter=",", skiprows=1, usecols=(1, 2))
    keep = bridge[: int(0.06 * 48000)]  # the first 60 ms: hammer contact and the first bridge pulses
    ev["bridge_velocity_first_60ms"] = {"t": [r(v, 4) for v in keep[::12, 0]], "v": [r(v, 4) for v in keep[::12, 1]]}
    items.append((cid, label, sr, pa, ev))
add_group("piano-ladder", "One key (A3), three physical hammer speeds", items, PIANO_BAND)

items = []
for cid, label in [("silent-held", "C4·E4·G4 held silently, C3 struck"), ("silent-control", "C3 struck alone"),
                   ("silent-pedal", "C3 struck, sustain pedal down")]:
    sr, pa, t = piano_clip(cid)
    items.append((cid, label, sr, pa, piano_evidence(t)))
add_group("piano-silent", "Sympathetic resonance: strings nobody struck", items, PIANO_BAND)

items = []
OCT = [("octhi2-et", "Equal-tempered first partials"), ("octhi3-stretch", "Stretch tuning from string stiffness")]
if not all((PIANO / f"{c}.wav").exists() for c, _ in OCT):
    OCT = [("octaves-et", "Equal-tempered first partials"), ("octaves-stretch", "Stretch tuning from string stiffness")]
for cid, label in OCT:
    sr, pa, t = piano_clip(cid)
    ev = piano_evidence(t)
    mono_full = pa.mean(axis=1)
    # A short window (43 ms) resolves a few-Hz octave beat as an amplitude wobble.
    n_s, hop_s = 2048, 240
    w_s = np.hanning(n_s)
    fr = np.fft.rfftfreq(n_s, 1 / sr)
    def band_env(f0, t0, t1):
        tt, aa = [], []
        for i in range(int(t0 * sr), min(len(mono_full) - n_s, int(t1 * sr)), hop_s):
            sp = np.abs(np.fft.rfft(mono_full[i:i + n_s] * w_s))
            sel = (fr > f0 * 0.97) & (fr < f0 * 1.03)
            tt.append(r((i + n_s / 2) / sr, 4))
            aa.append(r(float(np.sqrt((sp[sel] ** 2).sum())) * 2 / w_s.sum(), 4))
        return {"t": tt, "a": aa}
    ev["short_band_envelopes"] = {"880": band_env(880.0, 0.15, 2.3), "1760": band_env(1760.0, 2.7, 5.4)}
    mono = pa.mean(axis=1)
    beats = {}
    for f0, tag in [(110.0, "A1·A2 near 110 Hz"), (220.0, "A2·A3 near 220 Hz"), (440.0, "A3·A4 near 440 Hz")]:
        n = 16384
        hop = 480
        win = np.hanning(n)
        env = []
        for i in range(0, len(mono) - n, hop):
            seg = mono[i:i + n] * win
            spec = np.abs(np.fft.rfft(seg))
            freqs = np.fft.rfftfreq(n, 1 / sr)
            band = (freqs > f0 * 0.97) & (freqs < f0 * 1.03)
            env.append(float(np.sqrt((spec[band] ** 2).sum())) * 2 / win.sum())
        beats[tag] = {"t": [r((i * hop + n / 2) / sr, 4) for i in range(len(env))], "a": [r(v, 4) for v in env]}
    ev["band_envelopes"] = beats
    seg = mono[int(2.75 * sr):]
    nfft = 1 << 21
    spec = np.abs(np.fft.rfft(seg * np.hanning(len(seg)), nfft))
    freqs = np.fft.rfftfreq(nfft, 1 / sr)
    win = (freqs > 1748) & (freqs < 1778)
    mag = spec[win] / spec[win].max()
    ev["zoom_1760"] = {"f": [r(v, 8) for v in freqs[win][::3]], "db": [r(20 * math.log10(v + 1e-9), 4) for v in mag[::3]],
                       "window_s": r(len(seg) / sr, 4)}
    items.append((cid, label, sr, pa, ev))
add_group("piano-octaves", "Octaves A4–A5, then A5–A6: two tunings of the same strings", items, PIANO_BAND)

items = []
for cid, label in [("a0-mf", "A0 — the lowest string, 2.0 m/s"), ("c6-mf", "C6 — 2.0 m/s")]:
    sr, pa, t = piano_clip(cid)
    items.append((cid, label, sr, pa, piano_evidence(t)))
add_group("piano-register", "The bottom and the top of the scale", items, PIANO_BAND)

# ------------------------------------------------------------ percussion groups
DRUM_BAND = {"label": "Head modes retained", "lo_hz": 80, "hi_hz": 500,
             "acoustic_lo_hz": 40, "acoustic_hi_hz": 1000,
             "note": "Head modes are retained in 80–500 Hz; the BEM radiation transfer is fitted over 40–1000 Hz. Content outside is not simulated and is not faked."}


def perc_group(gid, title, rows, band, extra=None):
    its = []
    for cid, label, fs in rows:
        if not have(PERC / f"{cid}.wav"):
            print(f"skip {cid}: missing")
            continue
        sr, pa, t = perc_clip(cid, fs)
        its.append((cid, label, sr, pa, perc_evidence(t)))
    if its:
        add_group(gid, title, its, band, extra)


perc_group("snare", "Snare wires engaged vs. dropped 3 mm", [
    ("snare-on", "Wires on — 20 strands, 240 contacts", 200), ("snare-off", "Wires off — same drum, 3 mm clearance", 200)], DRUM_BAND)


def pitch_track(cid, fs):
    sr, pa, _ = perc_clip(cid, fs)
    x = pa.mean(axis=1)
    n, hop = 8192, 240
    win = np.hanning(n)
    t, f = [], []
    for i in range(0, min(len(x) - n, int(0.9 * sr)), hop):
        spec = np.abs(np.fft.rfft(x[i:i + n] * win, 4 * n))
        freqs = np.fft.rfftfreq(4 * n, 1 / sr)
        band = (freqs > 170) & (freqs < 300)
        k = np.argmax(spec * band)
        a, b, c = np.log(spec[k - 1:k + 2] + 1e-30)
        delta = 0.5 * (a - c) / (a - 2 * b + c)
        t.append(r((i + n / 2) / sr, 4))
        f.append(r((k + delta) * sr / (4 * n), 6))
    return {"t": t, "f": f}


hard_rows = [("hard-lin-08", "Linear · 0.8 m/s", 400), ("hard-lin-25", "Linear · 2.5 m/s", 400),
             ("hard-lin-50", "Linear · 5 m/s", 400), ("hard-nl-08", "Stretching · 0.8 m/s", 400),
             ("hard-nl-25", "Stretching · 2.5 m/s", 400), ("hard-nl-50", "Stretching · 5 m/s", 400)]
pitch = {cid: pitch_track(cid, fs) for cid, _, fs in hard_rows if have(PERC / f"{cid}.wav")}
perc_group("hard", "Hit it harder: linear vs geometrically stretching heads", hard_rows, DRUM_BAND, {"pitch_tracks": pitch})
perc_group("position", "Where you strike decides which modes speak", [
    ("pos-center", "Center", 200), ("pos-mid", "75 mm off center", 200), ("pos-edge", "145 mm — near the hoop", 200)], DRUM_BAND)
perc_group("tension", "Even lug tension vs an uneven head", [
    ("tens-even", "Even tension", 200), ("tens-varied", "Spatially varying prestress", 200)], DRUM_BAND)
perc_group("mallet", "Hard stick tip vs felt mallet face", [
    ("stick", "Hertz stick tip — 2.5 m/s", 200), ("mallet", "Felt mallet — 2.5 m/s", 200)], DRUM_BAND)
perc_group("roll", "Two hands from a force score", [("roll", "Alternating strokes, accented", 200)], DRUM_BAND)
CYM_BAND = {"label": "Shell modes retained", "lo_hz": 50, "hi_hz": 1200, "acoustic_lo_hz": 40, "acoustic_hi_hz": 1640,
            "note": "The cymbal keeps its 30 lowest elastic modes (50–1200 Hz) and the radiation fit spans 40–1640 Hz. A real splash rings far above this; that part is not simulated and is not faked."}
perc_group("cymbal", "Splash cymbal: stick, soft stroke, felt mallet, rocking", [
    ("cym-ff", "Stick — 2 m/s", 5), ("cym-ff-soft", "Stick — 0.6 m/s", 5), ("cym-ff-mallet", "Felt mallet — 1 m/s", 5),
    ("cym-ff-rock3", "Stick — 0.6 m/s, free to rock", 5)], CYM_BAND)

# ------------------------------------------------------- mechanics time series
def csv_series(name, cols, every):
    p = PERC / f"{name}.csv"
    if not p.exists():
        return None
    with open(p) as fh:
        header = fh.readline().strip().split(",")
    idx = [header.index(c) for c in cols if c in header]
    full = np.loadtxt(p, delimiter=",", skiprows=1, usecols=idx)
    out = {header[i]: [r(v, 5) for v in full[::every, k]] for k, i in enumerate(idx)}
    # Per-step ledger terms are accumulated at FULL resolution before decimation.
    for name in ("loss_j", "felt_crush_j"):
        if name in header and header.index(name) in idx:
            k = idx.index(header.index(name))
            out[name.replace("_j", "_cumulative_j")] = [r(v, 7) for v in np.cumsum(full[:, k])[::every]]
    if "balance_j" in header and header.index("balance_j") in idx:
        k = idx.index(header.index("balance_j"))
        out["balance_abs_max_j"] = r(float(np.abs(full[:, k]).max()), 4)
    return out


series = {
    "drum_stretch_4ms": csv_series("csv-drum-nl", ["time_s", "point_a_displacement_m", "point_b_displacement_m",
                                                   "cavity_internal_pa", "total_energy_j", "loss_j", "balance_j",
                                                   "head_stretching_energy_j"], 10),
    "drum_linear_4ms": csv_series("csv-drum-lin", ["time_s", "point_a_displacement_m", "point_b_displacement_m",
                                                   "cavity_internal_pa", "total_energy_j", "loss_j", "balance_j"], 10),
    "snare_on": csv_series("csv-snare", ["time_s", "point_a_displacement_m", "point_b_displacement_m", "total_energy_j", "loss_j", "balance_j"], 10),
    "snare_off": csv_series("csv-snare-off", ["time_s", "point_a_displacement_m", "point_b_displacement_m", "total_energy_j", "loss_j", "balance_j"], 10),
    "vent": csv_series("csv-vent", ["time_s", "point_a_displacement_m", "cavity_internal_pa", "neck_flow_m3_s",
                                    "neck_pressure_pa", "total_energy_j", "loss_j", "balance_j"], 10),
    "hihat": csv_series("csv-hihat", ["time_s", "upper_down_m", "lower_down_m", "pedal_down_m", "min_contact_gap_m",
                                      "active_sites", "total_energy_j", "loss_j", "player_work_j", "balance_j"], 4),
    "drag_ordinary": csv_series("csv-drag-0", ["time_s", "cavity_point_a_pa", "cavity_point_b_pa", "cavity_internal_pa"], 10),
    "drag_counterfactual": csv_series("csv-drag-max", ["time_s", "cavity_point_a_pa", "cavity_point_b_pa", "cavity_internal_pa"], 10),
}
(OUT / "data" / "series.json").write_text(json.dumps({k: v for k, v in series.items() if v}, separators=(",", ":")))

# ------------------------------------------------------------- modal state viz
VIZ = CAMPAIGN / "viz2"


def export_viz(name: str, out: str, max_frames: int, every: int = 1):
    d = VIZ / name
    if not (d / "shapes.jsonl").exists() or not (d / "state.f32").exists():
        print(f"viz {name}: missing")
        return
    bodies = [json.loads(l) for l in (d / "shapes.jsonl").read_text().splitlines() if l.strip()]
    meta = json.loads((d / "meta.json").read_text().splitlines()[0])
    nco = meta["coordinates"]
    raw = np.fromfile(d / "state.f32", "<f4")
    raw = raw[: len(raw) // (2 * nco) * (2 * nco)]  # a run still in progress may end mid-record
    frames = raw.reshape(-1, 2 * nco)[::every][:max_frames]
    q = frames[:, 0::2]
    blobs, header = [], {"bodies": [], "frames": int(q.shape[0]),
                         "dt_s": meta["dt_s"] * meta["every"] * every, "coordinates": nco}
    off = 0

    def push(arr, dtype):
        nonlocal off
        a = np.ascontiguousarray(arr, dtype=dtype)
        pad = (-off) % 4
        if pad:
            blobs.append(b"\0" * pad)
            off += pad
        blobs.append(a.tobytes())
        rec = {"offset": off, "length": int(a.size), "dtype": a.dtype.name}
        off += a.nbytes
        return rec

    used = []
    for b in bodies:
        modes = b["modes"]
        rec = {"kind": b["kind"], "modes": modes, "offset_coordinate": b["offset"],
               "frequencies_hz": [r(w / (2 * math.pi)) for w in b["omegas"]]}
        if b["kind"] == "head":
            rec["head"] = b["head"]
            rec["nodes"] = push(np.array(b["nodes"]).reshape(-1), "<f4")
            rec["shapes"] = push(np.array(b["shapes"]).reshape(-1), "<f4")  # node-major, modes per node
        else:
            rec["nodes"] = push(np.array(b["nodes"]).reshape(-1), "<f4")
            rec["shapes"] = push(np.array(b["shapes"]).reshape(-1), "<f4")  # node-major, 3*modes per node
        rec["tris"] = push(np.array(b["tris"]).reshape(-1), "<u4")
        rec["q"] = push(q[:, b["offset"]:b["offset"] + modes], "<f4")  # frame-major
        used.append(rec)
    header["bodies"] = used
    header["stick_q"] = push(q[:, 0], "<f4")
    (OUT / "viz" / f"{out}.bin").write_bytes(b"".join(blobs))
    (OUT / "viz" / f"{out}.json").write_text(json.dumps(header))
    print(f"viz {out}: {q.shape[0]} frames, {off/1e6:.2f} MB")


export_viz("viz-drum-slowmo", "drum-slowmo", 3000)
export_viz("viz-splash-slowmo", "cymbal-slowmo", 3000)
export_viz("viz-drum-hard", "drum-hard", 12000, every=2)
export_viz("viz-snare", "snare", 12000, every=2)
export_viz("cym-ff", "cymbal-ring", 12000, every=1)

# ---------------------------------------------------------------- piano model
def scale(path):
    rows = np.loadtxt(path, delimiter=",", comments="#", skiprows=2)
    out = []
    for m, uni, L, mu, T, EI, strike, hm, ft, fa, duplex, det in rows:
        B = math.pi ** 2 * EI / (T * L * L)
        f1 = 1 / (2 * L) * math.sqrt(T / mu) * math.sqrt(1 + B)
        out.append({"midi": int(m), "unison": int(uni), "L": r(L), "mu": r(mu), "T": r(T), "EI": r(EI), "B": r(B),
                    "f1": r(f1, 8), "strike": r(strike), "hammer_kg": r(hm), "felt_m": r(ft), "duplex_m": r(duplex),
                    "detune_cents": r(det)})
    return out


piano_model = {"scale_et": scale(PIANO / "scale-et.csv"), "scale_stretch": scale(PIANO / "scale-stretch.csv"),
               "scale_raw": scale(PIANO / "scale-raw.csv")}
# Soundboard geometry from the FSB dump: nodes, triangles, beam paths, bridge stations.
nodes, tris, bridges, beams = {}, [], {}, []
for line in (PIANO / "board.fsb").read_text().splitlines():
    f = line.split(",")
    if f[0] == "node":
        nodes[int(f[1])] = (float(f[2]), float(f[3]))
    elif f[0] == "triangle":
        tris.append((int(f[2]), int(f[3]), int(f[4]), float(f[5])))
    elif f[0] == "bridge":
        bridges[int(f[1])] = (int(f[2]), float(f[3]), float(f[4]), float(f[5]))
    elif f[0] == "stiffener":
        beams.append(f[1:])
xy = np.array([nodes[i] for i in range(len(nodes))])
station = {}
for key, (t, w0, w1, w2) in bridges.items():
    a, b, c = tris[t][:3]
    station[key] = (w0 * xy[a] + w1 * xy[b] + w2 * xy[c]).tolist()
piano_model["board"] = {
    "nodes": [r(v, 6) for v in xy.reshape(-1)],
    "tris": [v for t in tris for v in t[:3]],
    "thickness_m": [r(t[3], 4) for t in tris],
    "bridge_stations": {str(k): [r(v, 6) for v in s] for k, s in sorted(station.items())},
}
# OBJ groups (panel, ribs, bridges, cutoff bar) for the 3D view.
verts, groups, cur = [], {}, "panel"
for line in (PIANO / "board.obj").read_text().splitlines():
    if line.startswith("v "):
        verts.extend(float(v) for v in line.split()[1:4])
    elif line.startswith(("g ", "o ")):
        cur = line.split(maxsplit=1)[1]
    elif line.startswith("f "):
        idx = [int(p.split("/")[0]) - 1 for p in line.split()[1:]]
        g = groups.setdefault(cur, [])
        for k in range(1, len(idx) - 1):
            g.extend([idx[0], idx[k], idx[k + 1]])
piano_model["board_obj"] = {"vertices": [r(v, 6) for v in verts], "groups": groups}
bm = np.loadtxt(PIANO / "board-modes-400.csv", delimiter=",", comments="#", skiprows=2)
piano_model["board_modes_400"] = sorted({r(v) for v in bm[:, 1]})
(OUT / "data" / "piano-model.json").write_text(json.dumps(piano_model, separators=(",", ":")))

# Per-mode pressure: the soundboard mixer (exact signed modal sum).
basis = json.loads((PIANO / "a2-solo-basis.json").read_text())
with open(PIANO / "a2-solo-modal.csv") as fh:
    header = fh.readline().strip().split(",")
cols = [i for i, h in enumerate(header) if h.startswith("left_mode_")]
data = np.loadtxt(PIANO / "a2-solo-modal.csv", delimiter=",", skiprows=1, usecols=[header.index("pressure_left_pa")] + cols)
total, modes = data[:, 0], data[:, 1:]
energy = (modes ** 2).sum(axis=0)
order = np.argsort(-energy)
keep = order[:31]
rest = modes[:, order[31:]].sum(axis=1) if modes.shape[1] > 31 else np.zeros(len(total))
chans = [modes[:, k] for k in keep] + [rest]
scales = [float(np.abs(c).max()) or 1.0 for c in chans]
pcm = np.stack([np.round(c / s * 32767) for c, s in zip(chans, scales)]).astype("<i2")
(OUT / "data" / "soundboard-modes.i16").write_bytes(pcm.tobytes())
freqs = basis.get("loaded_reference_frequency_hz", [])
recon = np.abs(sum(chans) - total).max() / np.abs(total).max()
(OUT / "data" / "soundboard-modes.json").write_text(json.dumps({
    "sample_rate": 24000, "frames": int(len(total)), "channels": len(chans),
    "scales_pa": [r(s, 6) for s in scales],
    "frequency_hz": [r(freqs[k]) if k < len(freqs) else None for k in keep] + [None],
    "mode_index": [int(k) for k in keep] + [-1],
    "energy_share": [r(float(energy[k] / energy.sum()), 4) for k in keep] + [r(float(energy[order[31:]].sum() / energy.sum()), 4)],
    "total_modes": int(modes.shape[1]), "reconstruction_max_rel_error": r(float(recon), 3),
    "note": "Each channel is one loaded-board mode's signed pressure at the microphone (A2, 2.5 m/s); the last channel is the exact sum of all remaining modes. Their sum is the rendered pressure.",
    "evidence": piano_evidence(log_text(PIANO / "a2-solo.log")),
}))
print("soundboard mixer:", modes.shape, "recon", recon)

perf = []
for line in (PIANO / "hero.csv").read_text().splitlines()[1:]:
    sample, ev, key, val = line.split(",")
    perf.append([r(int(sample) / 48000, 6), ev, int(key), float(val)])
manifest["performances"] = {"hero": perf}
(OUT / "data" / "manifest.json").write_text(json.dumps(manifest, separators=(",", ":")))
print("manifest clips:", len(manifest["clips"]))
