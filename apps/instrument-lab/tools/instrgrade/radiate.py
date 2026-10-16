#!/usr/bin/env python3
"""DIAGNOSTIC-ONLY: never a graded render.

Quickie test: does moving the observer from the mouthpiece to the far field close the gap?

Independent textbook TMM (not fs-duct): lossy plane-wave cylinders, spherical-wave cone
segments, Levine-Schwinger unflanged or flanged-piston radiation load. The transfer
H(f) = p_far / p_mouthpiece with p_far ~ j w rho U_end / (4 pi r) (compact monopole, valid while
ka < ~1; above that the far field of a flaring bell is direction-dependent: this is a screen, not a
radiation model). The filter is applied zero-phase to the simulated mouthpiece pressure, then the
result is a diagnostic screen. grade_run.py refuses its output.
"""
from __future__ import annotations

import json
import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from instrgrade import GRADER_VERSION, SR, describe, distance, highpass, hz_of, load, midi_of, nearest_instruments, onset, write_wav  # noqa: E402

RHO, C = 1.204, 343.2


def k_lossy(f, r):
    w = 2 * math.pi * f
    alpha = 3.0e-5 * np.sqrt(f) / r  # Np/m, wide-tube viscothermal (Benade)
    return w / C * (1 + alpha * C / w) - 1j * alpha  # dispersion and attenuation of the same order (Benade)


def cyl(f, r, L):
    k = k_lossy(f, r)
    Zc = RHO * C / (math.pi * r * r)
    return np.array([[np.cos(k * L), 1j * Zc * np.sin(k * L)], [1j * np.sin(k * L) / Zc, np.cos(k * L)]])


def cone(f, r1, r2, L, n=48):
    """Staircase of n short cylinders at midpoint radii (adequate for a screen)."""
    M = np.eye(2, dtype=complex)
    for i in range(n):
        r = r1 + (r2 - r1) * (i + 0.5) / n
        M = M @ cyl(f, r, L / n)
    return M


def zrad(f, r, flanged):
    k = 2 * math.pi * f / C
    Zc = RHO * C / (math.pi * r * r)
    if flanged:
        return Zc * (0.5 * (k * r) ** 2 + 1j * 0.8216 * k * r)
    return Zc * (0.25 * (k * r) ** 2 + 1j * 0.6133 * k * r)


def transfer(fs, segs, flanged):
    H = np.zeros(len(fs), complex)
    for i, f in enumerate(fs):
        if f <= 0:
            continue
        M = np.eye(2, dtype=complex)
        for s in segs:
            M = M @ (cyl(f, s[1], s[2]) if s[0] == "cyl" else cone(f, s[1], s[2], s[3]))
        r_end = segs[-1][1] if segs[-1][0] == "cyl" else segs[-1][2]
        ZL = zrad(f, r_end, flanged)
        # [p_in; U_in] = M [p_end; U_end], U_end = p_end / ZL
        p_in = M[0, 0] * ZL + M[0, 1]  # per unit U_end
        w = 2 * math.pi * f
        p_far = 1j * w * RHO / (4 * math.pi * 1.0)  # per unit U_end at 1 m
        H[i] = p_far / p_in
    return H


def apply(x, segs, flanged):
    n = 1 << int(np.ceil(np.log2(len(x) + SR)))
    X = np.fft.rfft(x, n)
    f = np.fft.rfftfreq(n, 1 / SR)
    fg = np.geomspace(20, 20000, 600)
    Hg = np.abs(transfer(fg, segs, flanged))
    H = np.interp(f, fg, Hg, left=0, right=Hg[-1])
    return np.fft.irfft(X * H, n)[: len(x)], fg, Hg


if __name__ == "__main__":
    print("DIAGNOSTIC-ONLY: never a graded render", flush=True)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("library", type=Path)
    parser.add_argument("wav", type=Path)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--instrument", required=True)
    parser.add_argument("--note", required=True)
    parser.add_argument("--segments", required=True, help='JSON, e.g. [["cyl", 0.0073, 0.57]]')
    parser.add_argument("--flanged", action="store_true")
    parser.add_argument("--window", type=float, nargs=2, metavar=("START", "END"))
    args = parser.parse_args()
    lib = json.loads(args.library.read_text())
    x = load(args.wav)
    if args.window:
        x = x[int(args.window[0] * SR):int(args.window[1] * SR)]
    x = highpass(x)
    x = x[onset(x):]
    y, fg, Hg = apply(x, json.loads(args.segments), args.flanged)
    hint = hz_of(midi_of(args.note))
    d0, d1 = describe(x, "sustained", hint), describe(y, "sustained", hint)
    midi = int(round(69 + 12 * math.log2(d1["f0"] / 440)))
    es = [e for e in lib if e["inst"] == args.instrument and abs(e["note"] - midi) <= 1 and e["dyn"] == "mf"]
    if not es:
        raise SystemExit("IRPGRADE_REFERENCE_MISSING: diagnostic needs an mf reference within one semitone")
    best0 = min(distance(d0, e["desc"])["spec_db"] for e in es)
    best1 = min(distance(d1, e["desc"])["spec_db"] for e in es)
    result = dict(grader_version=GRADER_VERSION, radiated=True, diagnostic_only=True,
                  label="DIAGNOSTIC-ONLY: never a graded render", mouthpiece_spec_db=best0,
                  diagnostic_spec_db=best1, reference_files=[e["file"] for e in es],
                  nearest=nearest_instruments(lib, d1, d1["f0"], "sustained")[:3])
    args.outdir.mkdir(parents=True, exist_ok=True)
    stem = args.wav.stem + "-radiated"
    write_wav(args.outdir / (stem + ".wav"), y / max(1e-12, np.abs(y).max()) * 0.8)
    (args.outdir / (stem + ".json")).write_text(json.dumps(result, indent=2))
    print(json.dumps(result))
