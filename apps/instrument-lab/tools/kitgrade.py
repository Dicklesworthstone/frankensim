"""Grade the fs-vkplate drum kit against a real kit, piece by piece.

usage: uv run --with numpy --with pillow python tools/kitgrade.py KIT_DIR SALAMANDER_OH_DIR PUBLIC_DIR

KIT_DIR holds fs-vkplate renders (crates/fs-vkplate, `kit PIECE OUT.wav`): kick.wav, snare.wav,
snare-off.wav, tom-12.wav, tom-14.wav, hihat-closed.wav, hihat-open.wav, crash.wav, ride.wav.
SALAMANDER_OH_DIR holds the overhead files of the Salamander Drumkit (Alexander Holm, CC BY-SA 3.0,
archive.org/details/SalamanderDrumkit), whose sizes the kit copies.

Every descriptor is computed identically on both signals and each gap is placed against a ruler
made of REAL recordings: two strokes of the same drum at the same dynamic (what a drummer cannot
repeat), and the same drum one dynamic apart.
Writes PUBLIC_DIR/data/kit.json, level-matched excerpts in audio/ and spectrograms in img/.
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import verisim as V  # noqa: E402  (same reader, onset, spectrogram and encoder)

SR = V.SR
SOURCE = "https://archive.org/details/SalamanderDrumkit"
RIGHTS = "Salamander Drumkit by Alexander Holm, CC BY-SA 3.0 (creativecommons.org/licenses/by-sa/3.0); excerpts are trimmed and level-matched"

# id, title, ours, real stroke, a second real stroke of the same dynamic, the same drum one dynamic apart, seconds
PIECES = [
    ("kick", "Kick 18x18, felt beater", "kick", "kick_OH_F_1", "kick_OH_F_5", "kick_OH_P_3", 1.0),
    ("snare", "Snare 14x6, wires on", "snare", "snare2_OH_F_4", "snare2_OH_F_9", "snare2_OH_MP_5", 1.0),
    ("snare-off", "Snare 14x6, wires off", "snare-off", "snare2OFF_OH_F_1", "snare2OFF_OH_F_6", "snare2OFF_OH_P_2", 1.0),
    ("tom-12", "Rack tom 12x7", "tom-12", "hiTom_OH_F_4", "hiTom_OH_F_8", "hiTom_OH_P_3", 2.0),
    ("tom-14", "Floor tom 14x14", "tom-14", "loTom_OH_FF_5", "loTom_OH_FF_2", "loTom_OH_MP_3", 2.0),
    ("hihat-closed", "Hi-hat 14, closed", "hihat-closed", "hihatClosed_OH_F_13", "hihatClosed_OH_F_4", "hihatClosed_OH_P_8", 0.6),
    ("hihat-open", "Hi-hat 14, open", "hihat-open", "hihatOpen_OH_F_1", "hihatOpen_OH_F_3", "hihatOpen_OH_P_3", 2.5),
    ("crash", "Crash 18", "crash", "crash1_OH_FF_6", "crash1_OH_FF_2", "crash1_OH_P_2", 3.0),
    ("ride", "Ride 20", "ride", "ride1_OH_MP_1", "ride1_OH_MP_4", "ride1_OH_FF_2", 3.0),
]
OCT = [63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000]


def highpass(x):
    """20-30 Hz zero-phase high-pass: the kick's 44 Hz fundamental must survive (the piano grader cuts at 70)."""
    return V.highpass(x, 20.0, 30.0)


def clip(path: Path, dur: float):
    return V.segment(highpass(V.load(path)), dur)


def fundamental(x):
    f, s = V.spectrum(x, 0.0, 0.25, 1 << 18)
    sel = (f > 30) & (f < 1200)
    return float(f[np.argmax(np.where(sel, s, 0))])


def shares(x):
    f, s = V.spectrum(x, 0.0, 0.25, 1 << 17)
    p = s ** 2
    tot = p[f > 30].sum()
    return float(100 * p[f > 1000].sum() / tot), float(100 * p[f > 4000].sum() / tot), float((f * p).sum() / p.sum())


def decay_ms(x):
    """Per octave band: time from the band's peak to 20 dB below it (None if it never falls that far)."""
    n = int(0.01 * SR)
    fr = np.array([x[i:i + n] for i in range(0, len(x) - n, n)])
    S = np.abs(np.fft.rfft(fr * np.hanning(n), axis=1)) ** 2
    f = np.fft.rfftfreq(n, 1 / SR)
    out = []
    for fc in OCT[1:8]:
        e = 10 * np.log10(S[:, (f >= fc / math.sqrt(2)) & (f < fc * math.sqrt(2))].sum(1) + 1e-30)
        k = int(np.argmax(e))
        d = np.where(e[k:] < e[k] - 20)[0]
        out.append(int(d[0] * 10) if len(d) else None)
    return out


def balance(x):
    return V.octave_levels(x, 0.0, 0.25)


def log_decay_gap(a, b):
    """RMS of log2 decay-time ratios over bands where both fall 20 dB, in octaves of time."""
    v = [math.log2(p / q) for p, q in zip(a, b) if p and q]
    return float(np.sqrt(np.mean(np.square(v)))) if v else None


def main(kit_dir: Path, oh: Path):
    (V.OUT / "img").mkdir(parents=True, exist_ok=True)
    (V.OUT / "audio").mkdir(parents=True, exist_ok=True)
    pairs = []
    for pid, title, ours_id, real_id, again_id, other_id, dur in PIECES:
        p = kit_dir / f"{ours_id}.wav"
        if not p.exists():
            print("skip", pid)
            continue
        real, ours = clip(oh / f"{real_id}.wav", dur), clip(p, dur)
        again, other = clip(oh / f"{again_id}.wav", dur), clip(oh / f"{other_id}.wav", dur)
        rn, on = real / V.early_rms(real), ours / V.early_rms(ours)
        play = 10 ** (-1 / 20) / max(np.abs(rn).max(), np.abs(on).max())
        d = {}
        sr_, so = shares(real), shares(ours)
        d["fundamental_hz"] = {"real": V.r(fundamental(real), 4), "ours": V.r(fundamental(ours), 4), "note": "strongest component 30-1200 Hz, first 250 ms"}
        d["centroid_hz"] = {"real": V.r(sr_[2], 4), "ours": V.r(so[2], 4), "note": "spectral centroid, first 250 ms"}
        d["above_1k_pct"] = {"real": V.r(sr_[0], 3), "ours": V.r(so[0], 3), "note": "share of energy above 1 kHz, first 250 ms"}
        d["above_4k_pct"] = {"real": V.r(sr_[1], 3), "ours": V.r(so[1], 3), "note": "share of energy above 4 kHz, first 250 ms"}
        dr, do = decay_ms(real), decay_ms(ours)
        d["decay_ms"] = {"bands_hz": OCT[1:8], "real": dr, "ours": do, "note": "time for each octave band to fall 20 dB from its peak; a dash never falls that far in the excerpt"}
        br, bo, ba, bx = balance(real), balance(ours), balance(again), balance(other)
        gap = V.rms_db(br, bo)
        d["balance_gap_db"] = {"real": [V.r(v, 3) for v in br], "ours": [V.r(v, 3) for v in bo], "bands_hz": OCT, "rms_gap_db": V.r(gap, 3),
                               "note": "octave-band levels relative to total, first 250 ms"}
        pairs.append({
            "id": pid, "title": title, "duration_s": dur,
            "real": {"file": real_id, "url": SOURCE, "rights": RIGHTS, "sha256_wav": V.sha(oh / f"{real_id}.wav"),
                     "recording": "overhead microphones in a room (first channel used), 24-bit 48 kHz"},
            "ours": {"clip": ours_id, "label": "fs-vkplate render: one finite-difference plate solver for every piece; overhead point about 1 m away"},
            "audio": {"real": V.encode(rn * play, f"kit-{pid}-real"), "ours": V.encode(on * play, f"kit-{pid}-ours"),
                      "level_match": "both scaled to equal RMS over the first 300 ms after onset"},
            "images": V.images(rn, on, f"kit-{pid}"),
            "descriptors": d,
            "ruler": {"real_same_stroke_db": V.r(V.rms_db(br, ba), 3), "real_other_dynamic_db": V.r(V.rms_db(br, bx), 3), "ours_db": V.r(gap, 3),
                      "decay_same_stroke_oct": V.r(log_decay_gap(dr, decay_ms(again)), 3), "decay_ours_oct": V.r(log_decay_gap(dr, do), 3),
                      "strokes": [real_id, again_id, other_id]},
        })
        print(pid, {k: (v.get("real"), v.get("ours")) for k, v in d.items() if not isinstance(v.get("real"), list)}, "gap", round(gap, 1),
              "ruler", pairs[-1]["ruler"]["real_same_stroke_db"], pairs[-1]["ruler"]["real_other_dynamic_db"])
    groove = build_groove(kit_dir, oh)
    (V.OUT / "data" / "kit.json").write_text(json.dumps({"pairs": pairs, "source": SOURCE, "rights": RIGHTS, **({"groove": groove} if groove else {})}, separators=(",", ":")))


# Two bars at 100 bpm on a 16th grid: crash and kick on the downbeat, closed hi-hat eighths, snare on 2
# and 4, the ride in bar 2, a tom fill on beat 4, and the crash again on the last downbeat.
BPM = 100
GROOVE = (
    [("crash", 0), ("kick", 0), ("kick", 8), ("kick", 10), ("kick", 16), ("kick", 24), ("kick", 26), ("crash", 32), ("kick", 32)]
    + [("snare", s) for s in (4, 12, 20)]
    + [("hihat-closed", s) for s in range(2, 16, 2)]
    + [("ride", s) for s in range(16, 28, 2)]
    + [("tom-12", 28), ("tom-12", 29), ("tom-14", 30), ("tom-14", 31)]
)
REAL_OF = {pid: real_id for pid, _, _, real_id, *_ in PIECES}


def build_groove(kit_dir: Path, oh: Path):
    need = sorted({p for p, _ in GROOVE})
    if not all((kit_dir / f"{p}.wav").exists() for p in need):
        print("groove skipped: missing renders")
        return None
    step = 60 / BPM / 4
    total = int((33 * step + 3.0) * SR)
    mix = {"real": np.zeros(total), "ours": np.zeros(total)}
    for piece in need:
        real = V.segment(highpass(V.load(oh / f"{REAL_OF[piece]}.wav")), 3.0)
        ours = V.segment(highpass(V.load(kit_dir / f"{piece}.wav")), 3.0)
        # FrankenSim's level for each piece is matched to the real kit's own balance, piece by piece
        ours = ours * V.early_rms(real) / V.early_rms(ours)
        for p, s in GROOVE:
            if p != piece:
                continue
            i = int(round(s * step * SR))
            for key, x in (("real", real), ("ours", ours)):
                n = min(len(x), total - i)
                mix[key][i:i + n] += x[:n]
    peak = max(np.abs(mix["real"]).max(), np.abs(mix["ours"]).max())
    g = 10 ** (-1 / 20) / peak
    return {"real": V.encode(mix["real"] * g, "kit-groove-real"), "ours": V.encode(mix["ours"] * g, "kit-groove-ours"),
            "pattern": "Two bars at 100 bpm: crash and kick on the downbeat, closed hi-hat in eighths, snare on 2 and 4, the ride in bar 2, a tom fill, and the crash again. Each kit repeats one stroke per piece, so the two differ only in how each piece sounds.",
            "level_match": "FrankenSim's level for each piece is matched to the real kit's recorded balance (equal RMS over each stroke's first 300 ms); both grooves are then scaled together."}


if __name__ == "__main__":
    V.configure(".", sys.argv[2], sys.argv[3])
    main(Path(sys.argv[1]), Path(sys.argv[2]))
