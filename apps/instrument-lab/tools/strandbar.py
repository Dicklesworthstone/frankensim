"""Build the Strandbar section's data: the studio loop and the live engines' score and levels.

usage: uv run --with numpy --with soundfile python tools/strandbar.py ARRANGEMENT.json STEMS_DIR RECORDING.wav PUBLIC_DIR

ARRANGEMENT.json is written by the arrangement script (bars, voicings, strokes; see the section's
evidence). STEMS_DIR holds the studio renders, all starting `pre_s` before the first downbeat:
  piano.wav   grand_piano --performance strandbar-piano.csv (log: piano.log)
  kick.wav    fs-vkplate `perform kick hits-kick.csv`        (log: kick.log)
  hihat.wav   fs-vkplate `perform hihat-closed ...`          (log: hihat.log, or hihat-*.log per chunk)
RECORDING.wav is the reference recording; only its band balance is used.

The three stems are mixed with one fixed gain each, chosen so the mix's low (30-90 Hz), middle
(250-2500 Hz) and high (6-14 kHz) bands stand in the same proportion as the recording's. The
loop's tail is folded onto its start so it repeats seamlessly. The live engines get the same
strokes and keys, and fixed gains that give each live part the RMS of its studio part
(measured by running the live engines' own WebAssembly in Node: tools/strandbar_live_rms.mjs).
One master level serves both, the highest at which the studio loop peaks at -1 dBFS and the
live parts' sum at -2 dBFS (a dB to spare for the board and grids a device picks).
"""
from __future__ import annotations
import json, math, re, subprocess, sys
from pathlib import Path
import numpy as np
import soundfile as sf

SR = 48000
BANDS = [(30, 90), (250, 2500), (6000, 14000)]
HAT_STRIKE_X = 0.85 * 0.1778  # the offline hi-hat stroke lands at 0.85 of the 14-inch radius
KICK_XY = (0.0, 0.02)


def load(p: Path) -> np.ndarray:
    x, sr = sf.read(p, dtype="float64", always_2d=True)
    assert sr == SR, f"{p}: {sr} Hz"
    return x[:, 0]


def band_energy(x: np.ndarray) -> np.ndarray:
    n = 1 << 15
    hop = n // 2
    acc = np.zeros(len(BANDS))
    f = np.fft.rfftfreq(n, 1 / SR)
    w = np.hanning(n)
    for s in range(0, max(1, len(x) - n), hop):
        S = np.abs(np.fft.rfft(x[s:s + n] * w)) ** 2
        acc += [S[(f >= a) & (f < b)].sum() for a, b in BANDS]
    return acc


def balance(stems: dict[str, np.ndarray], target: np.ndarray) -> dict[str, float]:
    """Gains (amplitude) so that the mix's band energies are proportional to `target`."""
    E = {k: band_energy(v) + 1e-24 for k, v in stems.items()}
    owner = {"kick": 0, "piano": 1, "hihat": 2}
    g = {k: 1.0 for k in stems}
    for _ in range(200):
        mix = sum(g[k] ** 2 * E[k] for k in stems)
        err = np.log(mix / mix.sum()) - np.log(target / target.sum())
        for k, j in owner.items():
            g[k] = min(1e3, max(1e-3, g[k] * math.exp(-0.25 * err[j])))
    mix = sum(g[k] ** 2 * E[k] for k in stems)
    print("band shares  mix", np.round(mix / mix.sum(), 4), "recording", np.round(target / target.sum(), 4))
    return g


def wall_seconds(stems: Path, name: str) -> float:
    total = 0.0
    for log in sorted(stems.glob(f"{name}*.log")):
        t = log.read_text(errors="replace")
        m = re.search(r"wall_s (\d+)", t) or re.search(r"rendered in ([0-9.]+) s", t.replace("\r", "\n"))
        if m:
            total += float(m.group(1))
    return total


def main(arr_path: Path, stems: Path, recording: Path, pub: Path):
    arr = json.loads(arr_path.read_text())
    beat = 60 / arr["bpm"]
    pre = round(arr["pre_s"] * SR)
    loop = round(arr["loop_s"] * SR)
    raw = {k: load(stems / f"{k}.wav") for k in ("piano", "kick", "hihat")}
    n = min(len(v) for v in raw.values())
    raw = {k: v[:n] for k, v in raw.items()}
    rec = load(recording)
    target = band_energy(rec[60 * SR:290 * SR])
    g = balance(raw, target)
    mix = sum(g[k] * raw[k] for k in raw)
    G = 10 ** (-1 / 20) / np.abs(mix).max()
    body = mix[pre:pre + loop].copy()
    tail = mix[pre + loop:]
    body[:len(tail)] += tail  # the next pass rings on under the start of this one
    peak = np.abs(body).max() * G
    G *= min(1.0, 10 ** (-1 / 20) / peak)

    # the same score for the live engines, in frames from the first downbeat
    fr = lambda beats: round(beats * beat * SR)
    piano = []
    for e in arr["piano"]:
        piano.append({"offset": fr(e["t"]), "on": e["midi"], "speed": round(0.6 + 3.4 * e["vel"], 4)})
        piano.append({"offset": fr(e["t"] + e["dur"]), "off": e["midi"]})
    kick = [{"offset": fr(d["t"]), "x": KICK_XY[0], "y": KICK_XY[1], "speed": round(3.0 * d["vel"] / 0.9, 4)} for d in arr["drums"] if d["piece"] == "kick"]
    hihat = [{"offset": fr(d["t"]), "x": HAT_STRIKE_X, "y": 0.0, "speed": round(1.2 + 1.6 * d["vel"], 4)} for d in arr["drums"] if d["piece"] == "hihat-closed"]
    piano.sort(key=lambda e: (e["offset"], "on" in e))
    events = {"loopFrames": loop, "tailFrames": len(tail), "piano": piano, "kick": kick, "hihat": hihat}
    ev_path = stems / "live-events.json"
    ev_path.write_text(json.dumps(events))
    out = subprocess.run(["node", str(Path(__file__).parent / "strandbar_live_rms.mjs"), str(ev_path), str(pub.resolve()), "5", str(stems / "live")], capture_output=True, text=True, check=True).stdout
    live = json.loads(out.strip().splitlines()[-1])
    studio_rms = {k: float(np.sqrt(np.mean((g[k] * raw[k][pre:]) ** 2))) for k in raw}
    live_gain = {k: G * studio_rms[k] / live[k] for k in ("piano", "kick", "hihat")}
    # The live parts are summed with nothing limiting them, so their sum has to fit as well, and the
    # studio loop comes down with it to keep the two at the same level. Measured on Node's render (board
    # 32, 5x grid) to -2 dBFS: the boards and grids a device picks peak up to about 1 dB higher.
    lanes = {k: np.fromfile(stems / f"live-{k}.f32", dtype="<f4").astype(np.float64) for k in live_gain}
    live_peak = np.abs(sum(v * lanes[k] for k, v in live_gain.items())).max()
    fit = min(1.0, 10 ** (-2 / 20) / live_peak)
    G *= fit
    live_gain = {k: v * fit for k, v in live_gain.items()}
    print("gains", {k: round(v, 5) for k, v in g.items()}, "G", round(G, 5), "live", {k: round(v, 5) for k, v in live_gain.items()}, "live peak", round(live_peak * fit, 4), "piano live realtime (node)", round(live["pianoRealtime"], 2))

    body *= G
    (pub / "audio").mkdir(parents=True, exist_ok=True)
    tmp = pub / "audio" / ".strandbar-studio.wav"
    sf.write(tmp, body.astype(np.float32), SR, subtype="FLOAT")
    dst = pub / "audio" / "strandbar-studio.flac"
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", str(tmp), "-c:a", "flac", "-sample_fmt", "s32", "-compression_level", "8", str(dst)], check=True)
    tmp.unlink()

    shift = arr["transpose_semitones"]
    names = {"A": [(0, "D♭maj9"), (7, "A♭maj7(13)")], "B": [(0, "Gm7"), (7, "Dm7")], "C": [(0, "E♭m9"), (7, "Fm7")], "D": [(0, "Gm7"), (5, "E♭7/D♭"), (7, "Dm7(9)")]}
    up = {"D♭": "D", "A♭": "A", "G": "A♭", "D": "E♭", "E♭": "E", "F": "F♯"}
    def transpose(ch: str) -> str:
        if shift != 1:
            return ch
        out = ch
        for root in sorted(up, key=len, reverse=True):
            if ch.startswith(root):
                out = up[root] + ch[len(root):]
                break
        return re.sub(r"/(D♭)$", "/D", out)
    data = {
        "title": "Strandbar",
        "credit": "After Todd Terje, Strandbar (Olsen Records, 2013). The piano riff's rhythm and voicings follow a published transcription study of the record; here it is a semitone higher and at 130.45 BPM, like the version we matched. The drum part is four-on-the-floor with closed hi-hat sixteenths, accented as in that version.",
        "bpm": arr["bpm"], "beatsPerBar": arr["beats_per_bar"], "bars": arr["bars"], "loopFrames": loop,
        "form": [{"name": f, "chords": [{"slot": s, "name": transpose(c)} for s, c in names[f]]} for f in arr["form"]],
        "audio": "audio/strandbar-studio.flac",
        "piano": piano, "kick": kick, "hihat": hihat,
        "live": {k: {"gain": round(v, 8)} for k, v in live_gain.items()},
        "studio": [
            {"part": "Piano", "detail": "grand_piano, Steinway D preset: all 88 strings and their dampers on one FE soundboard (128 retained board coordinates up to 2.4 kHz), felt hammers driven by the score, a Rayleigh microphone 1 m above the board.", "computeS": wall_seconds(stems, "piano")},
            {"part": "Kick", "detail": "fs-vkplate: two Mylar heads on a 4.5 mm grid with the air between them, a convex felt beater crown, every stroke into the same ringing drum.", "computeS": wall_seconds(stems, "kick")},
            {"part": "Hi-hat", "detail": "fs-vkplate: two bronze shells on a 5.5 mm grid (von Kármán shells), clamped by the pedal at their rims, stepped at 288 kHz.", "computeS": wall_seconds(stems, "hihat")},
        ],
        "liveParts": [
            {"part": "Piano", "detail": f"engine-piano: the same hammer, string and damper equations for the keys being played, on {live['pianoInfo'].get('boardCoordinates', '?')} soundboard coordinates."},
            {"part": "Kick and hi-hat", "detail": "fs-vkplate itself, compiled to WebAssembly, on a grid several times coarser (picked for this device when you press play): fewer, larger cells, so less detail above a few kHz."},
        ],
        "balance": "Levels: one fixed gain per part. In the studio render they set the low (30–90 Hz), middle (250–2,500 Hz) and high (6–14 kHz) bands in the same proportion as the recording we matched; each live part then gets the gain that gives it the RMS of its studio part, corrected for the drum grid your device runs. Both versions share one master level, low enough that neither clips. Nothing is compressed, limited or equalised.",
    }
    (pub / "data").mkdir(parents=True, exist_ok=True)
    (pub / "data" / "strandbar.json").write_text(json.dumps(data, separators=(",", ":")))
    print("wrote", dst, f"{dst.stat().st_size / 1e6:.2f} MB", "and data/strandbar.json")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]))
