#!/usr/bin/env python3
"""Per-octave-band decays of card-built drum renders against the .365 references.

Bead pexfp.370, DONE-WHEN 3: a MEASURED finding, printed. Nothing is tuned to it. The renders
come from fs-vkplate's `kit` example, whose drums are built from the sourced material cards
by default:

    kit kick dw3-kick.wav 0 2.5; kit snare dw3-snare.wav 0 1.5
    kit tom-12 dw3-tom12.wav 0 3.0; kit tom-14 dw3-tom14.wav 0 3.0

The references are the VALID clips of data/instruments/references-screened.json. Each is pinned
by its clip_sha256; on trj they live under /data/tmp/wa-inst/refs/u4/clips/kit.

Estimator (declared before any run, applied identically to renders and references):
- The channels are averaged to one signal.
- Each octave band, fc in 63..8000 Hz, is cut from the whole signal by a zero-phase
  Butterworth band-pass magnitude (third-order prototype, edges fc/sqrt2 and fc*sqrt2,
  as in ANSI S1.11), applied to the signal's FFT.
- The band's energy is summed in 10 ms frames. The analysis ends 100 ms before the clip's
  end, so that the band filter's response to the cut is excluded. A render still ringing at
  its last sample ends in a step, which spreads into the final frames of every band.
  (Added after the first run, where that step produced spurious end-of-clip peaks; the
  cells it did not touch are unchanged.)
- decay = the time from the band's peak frame to the first later frame 20 dB below that peak.
  This is kitgrade.decay_ms's definition; only the band filter differs. kitgrade's 10 ms
  frame FFT has 100 Hz bins and cannot resolve the 63-250 Hz bands.
- A recording's band whose floor (the median frame energy of the clip's last 10 %) lies
  within 26 dB of its peak is NOT MEASURABLE ("nm"). Renders have no noise, so they skip
  this test.
- A band that never falls 20 dB inside the clip reports ">" its remaining length, a lower
  bound.
- Each band's resolution floor is the decay the estimator reports for a unit impulse. A
  measured decay under 1.5 times that floor is marked "~" (at resolution).
- If a band decays exponentially, T60 = 3 x the 20 dB time.

usage: uv run --with numpy python3 -I kitdecay.py MANIFEST.json CLIPS_DIR RENDERS_DIR OUT.json
"""
import hashlib
import json
import math
import struct
import sys
from pathlib import Path

import numpy as np

OCT = [63, 125, 250, 500, 1000, 2000, 4000, 8000]
FRAME_S = 0.010
EDGE_S = 0.100
DROP_DB = 20.0
FLOOR_MARGIN_DB = 26.0

# Render file -> the references it is printed against (VALID clips of references-screened.json).
# The sizes of the reference drums are not stated in the manifest; DRSKit tom1 is taken as
# the rack tom and tom2/tom3 as the larger toms. Rimshots are a different stroke and are left out.
PIECES = [
    ("kick", "dw3-kick.wav", "kick 18x18, 2-ply batter, pillow; preset felt beater",
     ["kit-drskit-kdrum_with_contact-19", "kit-drskit-kdrum_without_contact-29"]),
    ("snare", "dw3-snare.wav", "snare 14x6, wires on; hickory stick card",
     ["kit-vcsl-snare3m-hitsn-v3-rr2-mid", "kit-vcsl-snare4-hitsn-v5-rr1-mid"]),
    ("tom-12", "dw3-tom12.wav", "rack tom 12x7; hickory stick card",
     ["kit-drskit-tom1-27", "kit-vcsl-tomh-hits-v3-rr2-mid", "kit-vcsl-tomh-hitm-v2-rr1-mid",
      "kit-vcsl-tomh-hitm-v2-rr2-mid", "kit-vcsl-tomh-hitm-v3-rr2-mid", "kit-vcsl-tomh-hitm-v4-rr2-mid"]),
    ("tom-14", "dw3-tom14.wav", "floor tom 14x14; hickory stick card",
     ["kit-drskit-tom2-10", "kit-drskit-tom2-11", "kit-drskit-tom2-20", "kit-drskit-tom2-21",
      "kit-drskit-tom2-30", "kit-drskit-tom3-11", "kit-drskit-tom3-20", "kit-drskit-tom3-21",
      "kit-drskit-tom3-30", "kit-vcsl-toml-hits-v3-rr2-mid", "kit-vcsl-toml-hitm-v2-rr1-mid",
      "kit-vcsl-toml-hitm-v3-rr1-mid", "kit-vcsl-toml-hitm-v3-rr2-mid", "kit-vcsl-toml-hitm-v4-rr2-mid"]),
]


def read_wav(path):
    data = Path(path).read_bytes()
    if data[:4] != b"RIFF" or data[8:12] != b"WAVE":
        raise ValueError(f"{path}: not a RIFF/WAVE file")
    pos, fmt, samples = 12, None, None
    while pos + 8 <= len(data):
        cid, size = data[pos:pos + 4], struct.unpack("<I", data[pos + 4:pos + 8])[0]
        body = data[pos + 8:pos + 8 + size]
        if cid == b"fmt ":
            tag, ch, sr = struct.unpack("<HHI", body[:8])
            bits = struct.unpack("<H", body[14:16])[0]
            if tag == 0xFFFE:
                tag = struct.unpack("<H", body[24:26])[0]
            fmt = (tag, ch, sr, bits)
        elif cid == b"data":
            samples = body
        pos += 8 + size + (size & 1)
    if fmt is None or samples is None:
        raise ValueError(f"{path}: no fmt or data chunk")
    tag, ch, sr, bits = fmt
    if tag == 3 and bits == 32:
        x = np.frombuffer(samples, "<f4").astype(np.float64)
    elif tag == 3 and bits == 64:
        x = np.frombuffer(samples, "<f8").astype(np.float64)
    elif tag == 1 and bits == 16:
        x = np.frombuffer(samples, "<i2") / 32768.0
    elif tag == 1 and bits == 24:
        b = np.frombuffer(samples, np.uint8).reshape(-1, 3).astype(np.int32)
        v = b[:, 0] | (b[:, 1] << 8) | (b[:, 2] << 16)
        x = np.where(v >= 1 << 23, v - (1 << 24), v) / float(1 << 23)
    else:
        raise ValueError(f"{path}: unsupported format tag {tag}, {bits} bits")
    x = x[: len(x) // ch * ch].reshape(-1, ch)
    return x.mean(axis=1), sr


def band_gain(f, fc, order=3):
    """Zero-phase magnitude of a Butterworth band-pass, edges fc/sqrt2 and fc*sqrt2."""
    lo, hi = fc / math.sqrt(2), fc * math.sqrt(2)
    w0, bw = math.sqrt(lo * hi), hi - lo
    with np.errstate(divide="ignore"):
        q = np.where(f > 0, (f * f - w0 * w0) / (f * bw), np.inf)
    return 1.0 / np.sqrt(1.0 + q ** (2 * order))


def band_energy_db(x, sr):
    n = 1 << int(math.ceil(math.log2(len(x) + sr)))  # pad by 1 s so the filter does not wrap
    spec = np.fft.rfft(x, n)
    f = np.fft.rfftfreq(n, 1.0 / sr)
    hop = int(round(FRAME_S * sr))
    out = []
    for fc in OCT:
        y = np.fft.irfft(spec * band_gain(f, fc), n)[: len(x)]
        frames = len(y) // hop
        e = (y[: frames * hop] ** 2).reshape(frames, hop).sum(axis=1)
        out.append(10.0 * np.log10(e + 1e-300))
    return out


def analysed_band_energy_db(x, sr):
    """band_energy_db without its last EDGE_S, where the band filter's response to the cut lies."""
    edge = int(round(EDGE_S / FRAME_S))
    return [e[: len(e) - edge] for e in band_energy_db(x, sr)]


def decays(x, sr, floors=None, recording=False):
    rows = []
    for k, e in enumerate(analysed_band_energy_db(x, sr)):
        p = int(np.argmax(e))
        tail = e[int(0.9 * len(e)):]
        floor = float(np.median(tail)) if len(tail) else -300.0
        below = np.where(e[p:] < e[p] - DROP_DB)[0]
        row = {"fc_hz": OCT[k], "peak_s": round(p * FRAME_S, 3), "peak_db": round(float(e[p]), 2),
               "floor_db": round(floor, 2)}
        if recording and e[p] - floor < FLOOR_MARGIN_DB:
            row["decay_ms"], row["mark"] = None, "nm"
        elif len(below):
            ms = float(below[0] * FRAME_S * 1000.0)
            row["decay_ms"] = round(ms, 1)
            row["mark"] = "~" if floors and ms < 1.5 * floors[k] else ""
        else:
            row["decay_ms"], row["mark"] = round((len(e) - p) * FRAME_S * 1000.0, 1), ">"
        rows.append(row)
    return rows


def show(rows):
    cells = []
    for r in rows:
        if r["mark"] == "nm":
            cells.append("     nm")
        else:
            cells.append(f"{r['mark']:>1}{r['decay_ms']:6.0f}")
    return " ".join(cells)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(manifest, clips_dir, renders_dir, out_path):
    m = json.loads(Path(manifest).read_text())
    clips = {}
    for f in m["files"]:
        for item in (f["members"] if f.get("kind") == "archive" else [f]):
            for c in item.get("clips", []):
                clips[c["id"]] = c
    sr0 = 48000
    impulse = np.zeros(3 * sr0)
    impulse[sr0 // 2] = 1.0
    floors = [r["decay_ms"] for r in decays(impulse, sr0)]
    report = {"bead": "frankensim-irp-v5-instrument-realism-pexfp.370", "done_when": 3,
              "estimator": __doc__.split("Estimator (")[1].split("usage:")[0].strip(),
              "bands_hz": OCT, "resolution_floor_ms": floors, "pieces": []}
    print("bands (Hz)            " + " ".join(f"{b:>7}" for b in OCT))
    print("resolution floor (ms) " + " ".join(f"{v:7.0f}" for v in floors))
    for name, render, what, refs in PIECES:
        x, sr = read_wav(f"{renders_dir}/{render}")
        r_rows = decays(x, sr, floors)
        piece = {"piece": name, "render": render, "render_sha256": sha(f"{renders_dir}/{render}"),
                 "render_seconds": round(len(x) / sr, 3), "what": what, "render_decays": r_rows,
                 "references": []}
        print(f"\n{name}: {what}")
        print(f"  render {render:<34}{show(r_rows)}")
        for cid in refs:
            c = clips[cid]
            if c["status"] != "VALID":
                raise SystemExit(f"{cid} is not VALID in the manifest")
            path = f"{clips_dir}/{cid}.wav"
            digest = sha(path)
            if digest != c["clip_sha256"]:
                raise SystemExit(f"{cid}: sha256 {digest} does not match the manifest")
            y, ysr = read_wav(path)
            rows = decays(y, ysr, floors, recording=True)
            piece["references"].append({"id": cid, "clip_sha256": digest, "dynamic": c.get("dynamic"),
                                        "articulation": c.get("articulation"), "seconds": round(len(y) / ysr, 3),
                                        "decays": rows})
            print(f"  ref {cid:<37}{show(rows)}")
        report["pieces"].append(piece)
    Path(out_path).write_text(json.dumps(report, indent=1))
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit(__doc__.split("usage: ")[1])
    main(*sys.argv[1:])
