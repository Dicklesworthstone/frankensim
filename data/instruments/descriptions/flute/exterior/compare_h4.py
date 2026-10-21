"""H4 test 1: old observer vs BEM re-observation vs the measured capsule pattern, per partial.

Usage: python compare_h4.py OUTDIR MEASURED.tsv BODY.tsv RUNG_A RUNG_B REPORT.json
OUTDIR holds Flute.nonvib.mf.<N>.wav (the existing observer, as rendered) and harmonics-<rung>-<note>.tsv
written by `reobserve_apertures harmonics` for rung A (reported) and rung B (its ladder partner).

Per note and partial k (inside the harmonics band):
- the PATTERN of each observer at the 32 capsules: dB relative to the capsule power mean, so the source strength
  cancels: BEM rung A and B at r206 (SOFA distance) and r105 (paper radius); measured (SOFA facts); and the OLD
  observer's own model (the same tapped flow phasors as point sources in one infinite baffle through the flute
  axis, openings facing up; undefined below the baffle, so compared on the 16 capsules above it);
- pattern error against the measured pattern: RMS over capsules [dB]; error bar: rung A vs rung B pattern RMS;
- the level change from the old observer (as rendered, with its band) to the BEM capsule-mean level, per partial.
Partial-power spectral centroids of the old WAV, the capsule-mean BEM spectrum and the capsule-mean measured
spectrum are reported as facts; they are not the verdict. Nothing here is fitted or tuned.
"""
import json
import struct
import sys

import numpy as np

out, measured_path, body_path, rung_a, rung_b, report_path = sys.argv[1:7]
NOTES = ("d4", "g4", "a5")
WINDOW_S = (1.0, 2.2)


def read_wav(path):
    b = open(path, "rb").read()
    assert b[:4] == b"RIFF" and b[8:12] == b"WAVE"
    i, fmt, rate, data, bits = 12, None, None, None, None
    while i < len(b):
        cid, size = b[i:i + 4], struct.unpack("<I", b[i + 4:i + 8])[0]
        body = b[i + 8:i + 8 + size]
        if cid == b"fmt ":
            fmt, _, rate = struct.unpack("<HHI", body[:8])
            bits = struct.unpack("<H", body[14:16])[0]
        elif cid == b"data":
            data = body
        i += 8 + size + (size & 1)
    if fmt == 3:
        return np.frombuffer(data, "<f4").astype(float), rate
    assert fmt == 1 and bits == 16
    return np.frombuffer(data, "<i2").astype(float), rate


def partial_rms(x, rate, f):
    n0, n1 = int(WINDOW_S[0] * rate), int(WINDOW_S[1] * rate)
    seg, w = x[n0:n1], np.hanning(n1 - n0)
    t = np.arange(n1 - n0) / rate
    return abs(np.sum(seg * w * np.exp(-2j * np.pi * f * t))) * 2 / w.sum() / np.sqrt(2)


def read_harmonics(path):
    f0, flows, pressure, freq = None, {}, {}, {}
    for line in open(path):
        if line.startswith("f0_hz\t"):
            f0 = float(line.split("\t")[1])
        elif line.startswith(("flow\t", "pressure\t")):
            kind, k, f, name, mag, ph = line.rstrip("\n").split("\t")
            k = int(k)
            freq[k] = float(f)
            (flows if kind == "flow" else pressure).setdefault(k, {})[name] = (float(mag), float(ph))
    return f0, freq, flows, pressure


def pattern(levels_pa):
    p = np.asarray(levels_pa, float)
    return 20 * np.log10(np.maximum(p, 1e-15)) - 10 * np.log10(np.mean(p ** 2))


def rms(a):
    return float(np.sqrt(np.mean(np.square(a))))


measured = {}
for line in open(measured_path):
    if line.startswith(("frankensim", "#", "note")):
        continue
    note, k, f, cap, az, el, dist, value = line.rstrip("\n").split("\t")
    measured.setdefault(note.lower(), {}).setdefault(int(k), {})[int(cap)] = (float(f), float(value), float(az), float(el))

positions, medium, profile_x = {}, None, []
for line in open(body_path):
    f = line.rstrip("\n").split("\t")
    if f[0] == "medium":
        medium = (float(f[1]), float(f[2]))
    elif f[0] == "wall":
        positions[f[1]] = float(f[2])
    elif f[0] == "profile":
        profile_x.append(float(f[1]))
HEAD, FOOT = min(profile_x), max(profile_x)
positions["foot"] = FOOT
RHO, C = medium
MID = 0.5 * (HEAD + FOOT)  # the declared pose of receivers.tsv: tube midpoint at the array centre

report = {"window_s": WINDOW_S, "rungs": [rung_a, rung_b], "notes": {}}
for note in NOTES:
    old, rate = read_wav(f"{out}/Flute.nonvib.mf.{note.upper()}.wav")
    f0, freq, flows, pa = read_harmonics(f"{out}/harmonics-{rung_a}-{note}.tsv")
    _, _, _, pb = read_harmonics(f"{out}/harmonics-{rung_b}-{note}.tsv")
    up = [c for c in range(32) if measured[note][1][c][3] > 0]
    rows = []
    for k in sorted(freq):
        if k not in measured[note]:
            break
        f = freq[k]
        meas = np.array([measured[note][k][c][1] for c in range(32)])
        lev = {(r, tag): np.array([p[k][f"{tag}-{c:02d}"][0] for c in range(32)])
               for r, p in ((rung_a, pa), (rung_b, pb)) for tag in ("r206", "r105")}
        kw = 2 * np.pi * f / C
        old_lev = []
        for c in up:
            az, el = np.radians(measured[note][k][c][2]), np.radians(measured[note][k][c][3])
            xa = 2.06 * np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
            p = 0j
            for name, (mag, ph) in flows[k].items():
                d = np.linalg.norm(xa - np.array([0.0, MID - positions[name], 0.0]))
                p += mag * np.exp(1j * ph) * np.exp(1j * kw * d) / d
            old_lev.append(abs(2 * np.pi * f * RHO / (2 * np.pi) * p) / np.sqrt(2))
        pm_all, pm_up = pattern(meas), pattern(meas[up])
        row = {"partial": k, "f_hz": f, "measured_f_hz": measured[note][k][0][0],
               "measured_spread_db": float(pm_all.max() - pm_all.min())}
        for (r, tag), l in lev.items():
            row[f"err_{r}_{tag}_all_db"] = rms(pattern(l) - pm_all)
            row[f"err_{r}_{tag}_up_db"] = rms(pattern(l[up]) - pm_up)
        row["err_old_model_up_db"] = rms(pattern(old_lev) - pm_up)
        for tag in ("r206", "r105"):
            row[f"ladder_{tag}_db"] = rms(pattern(lev[(rung_a, tag)]) - pattern(lev[(rung_b, tag)]))
        row["old_wav_rms"] = float(partial_rms(old, rate, f))
        row["bem_mean_rms_pa"] = float(np.sqrt(np.mean(lev[(rung_a, "r105")] ** 2)))
        row["measured_mean_rms_pa"] = float(np.sqrt(np.mean(meas ** 2)))
        rows.append(row)
    fk = np.array([r["f_hz"] for r in rows])
    cent = {}
    for key in ("old_wav_rms", "bem_mean_rms_pa", "measured_mean_rms_pa"):
        pw = np.array([r[key] for r in rows]) ** 2
        cent[key] = float(np.sum(fk * pw) / np.sum(pw))
    # Level change old -> BEM per partial, relative to the fundamental (the old WAV is PCM16 at its own scale).
    ref_old, ref_bem = rows[0]["old_wav_rms"], rows[0]["bem_mean_rms_pa"]
    for r in rows:
        r["bem_minus_old_rel_db"] = float(20 * np.log10(r["bem_mean_rms_pa"] / ref_bem) - 20 * np.log10(max(r["old_wav_rms"], 1e-12) / ref_old))
    report["notes"][note] = {"f0_hz": f0, "partials": rows, "centroid_hz": cent, "capsules_above_baffle": up}
    print(f"== {note} f0 {f0:.3f} Hz; centroid old {cent['old_wav_rms']:.0f} BEM {cent['bem_mean_rms_pa']:.0f} measured {cent['measured_mean_rms_pa']:.0f} Hz")
    print(" k      f  err_old_up  err_A_up(206)  err_A_all(206)  err_A_all(105)  ladder206  ladder105  spread  BEM-old(rel dB)")
    for r in rows:
        print(f"{r['partial']:2d} {r['f_hz']:7.0f}  {r['err_old_model_up_db']:9.1f}  {r[f'err_{rung_a}_r206_up_db']:12.1f}  "
              f"{r[f'err_{rung_a}_r206_all_db']:13.1f}  {r[f'err_{rung_a}_r105_all_db']:13.1f}  {r['ladder_r206_db']:8.2f}  "
              f"{r['ladder_r105_db']:8.2f}  {r['measured_spread_db']:6.1f}  {r['bem_minus_old_rel_db']:8.1f}")
json.dump(report, open(report_path, "w"), indent=1)
