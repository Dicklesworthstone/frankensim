#!/usr/bin/env python3
"""Grade a list of candidate renders against the UIowa library and write figures, pairs and JSON.

usage: grade_run.py SCRATCH CANDIDATES.json OUTDIR
Run with: uv run --with numpy --with matplotlib python3 -I grade_run.py ...
Each candidate: {"id", "label", "wav", "inst", "note" (name), "dyn", "art", "kind",
                 optional "t0"/"t1" (seconds: excerpt of the wav), "chord": [names], "stagger_s"}
Optional reproduction fields: expected_db, tolerance_db, conversion_cmd,
render_git_sha, render_command, performance_sha256.
Receipts report descriptor reproduction only, never instrument acceptance.
"""
from __future__ import annotations

import json
import hashlib
import math
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from instrgrade import (GRADER_VERSION, SR, cents, describe, distance, find, highpass, hz_of, load, midi_of, nearest_instruments,
                        null_baselines, onset, ref_audio, rulers, write_wav)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

if len(sys.argv) != 4:
    raise SystemExit("IRPGRADE_USAGE: grade_run.py SCRATCH CANDIDATES.json OUTDIR")
STARTED = time.monotonic()
S = Path(sys.argv[1])
CANDS = json.loads(Path(sys.argv[2]).read_text())
for candidate in CANDS:
    if candidate.get("radiated") is True or Path(candidate["wav"]).stem.lower().endswith("-radiated"):
        print(json.dumps(dict(error="IRPGRADE_RADIATED_CANDIDATE_REFUSED", id=candidate["id"],
                              message="DIAGNOSTIC-ONLY: never a graded render")), file=sys.stderr)
        raise SystemExit(2)
OUT = Path(sys.argv[3])
(OUT / "fig").mkdir(parents=True, exist_ok=True)
(OUT / "pairs").mkdir(parents=True, exist_ok=True)
WAVDIR = Path(os.environ["IRP_UIOWA_WAV48"]) if "IRP_UIOWA_WAV48" in os.environ else S / "real" / "wav48"
LIBPATH = S / "grade" / os.environ.get("LIBFILE", "library.json")
LIB = json.loads(LIBPATH.read_text())
TOOLDIR = Path(__file__).resolve().parent
TOOL_SHA256 = {name: hashlib.sha256((TOOLDIR / name).read_bytes()).hexdigest()
               for name in ("instrgrade.py", "descriptors.py", "grade_run.py", "radiate.py", "summarize.py", "buildlib.py", "uiowa-links.json")}
GIT = subprocess.run(["git", "-C", str(TOOLDIR), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=10)
CONTEXT = dict(grader_version=GRADER_VERSION, library_sha256=hashlib.sha256(LIBPATH.read_bytes()).hexdigest(),
               tool_sha256=TOOL_SHA256, host=platform.node(), git_sha=GIT.stdout.strip() or None,
               command=[sys.executable, "-I", *sys.argv], utc=datetime.now(timezone.utc).isoformat())
RECEIPTS = []


def receipt(candidate, result, started, target=None):
    dist = result.get("distance", {})
    expected, tolerance = candidate.get("expected_db"), candidate.get("tolerance_db")
    observed = dist.get("spec_db")
    if "error" in result:
        verdict = "error"
    elif expected is None:
        verdict = "measured"
    elif (tolerance is not None and math.isfinite(expected) and math.isfinite(tolerance)
          and tolerance >= 0 and math.isfinite(observed) and abs(observed - expected) <= tolerance):
        verdict = "pass"
    else:
        verdict = "fail"
    path = Path(candidate["wav"])
    return dict(CONTEXT, event="candidate", id=candidate["id"], verdict=verdict,
                candidate_wav_sha256=hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None,
                reference_file=target["file"] if target else None,
                reference_note=target["name"] if target else candidate.get("chord"),
                reference_dynamic=target["dyn"] if target else candidate.get("dyn"),
                reference_wav_sha256=hashlib.sha256((WAVDIR / target["file"]).read_bytes()).hexdigest() if target else None,
                windows=result.get("windows"),
                band_split_hz=result.get("band_split_hz"),
                spec_db=observed, env_db=dist.get("env_db"), harm_db=dist.get("harm_db"),
                attack_erb_distance_db=dist.get("attack_erb_distance_db"),
                sustain_erb_distance_db=dist.get("sustain_erb_distance_db"),
                release_erb_distance_db=dist.get("release_erb_distance_db"),
                q3=result.get("q3"),
                descriptors=result.get("descriptors"),
                baselines=result.get("baselines"),
                rulers=result.get("rulers"), expected_db=expected, tolerance_db=tolerance,
                conversion_cmd=candidate.get("conversion_cmd", os.environ.get("IRP_CONVERSION_CMD")),
                render_git_sha=candidate.get("render_git_sha"), render_command=candidate.get("render_command"),
                performance_sha256=candidate.get("performance_sha256"), error=result.get("error"),
                elapsed_s=time.monotonic() - started)


def stft_db(x, n=4096, hop=256):
    win = np.hanning(n)
    pad = np.concatenate([np.zeros(n // 2), x, np.zeros(n)])
    cols = max(1, (len(pad) - n) // hop)
    frames = np.stack([pad[i * hop:i * hop + n] * win for i in range(cols)])
    mag = np.abs(np.fft.rfft(frames, axis=1)) ** 2
    return 10 * np.log10(mag.T + 1e-20), np.fft.rfftfreq(n, 1 / SR), np.arange(cols) * hop / SR


def rms(x):
    return float(np.sqrt(np.mean(x * x))) if len(x) else 0.0


def prep_candidate(c):
    f0_exp = hz_of(midi_of(c["note"])) if "note" in c else (hz_of(midi_of(c["chord"][0])) if c.get("chord") else None)
    x = highpass(load(Path(c["wav"])), f0_exp)
    if "t0" in c:
        x = x[int(c["t0"] * SR):int(c["t1"] * SR)]
    o = onset(x)
    return x[o:]


def real_chord(names, dyn, inst, stagger):
    parts = []
    for k, nm in enumerate(names):
        es = find(LIB, inst, midi_of(nm), dyn=dyn)
        if not es:
            return None
        y = ref_audio(es[0], WAVDIR)
        y = np.concatenate([np.zeros(int(k * stagger * SR)), y])
        parts.append(y)
    n = max(len(p) for p in parts)
    return sum(np.pad(p, (0, n - len(p))) for p in parts)


def figure(cid, title, panels, desc_pairs, path):
    fig = plt.figure(figsize=(14, 8.6), dpi=110)
    gs = fig.add_gridspec(2, len(panels), height_ratios=[2.2, 1])
    vmax = None
    specs = []
    for lab, x in panels:
        d, f, t = stft_db(x[: int(3.0 * SR)])
        specs.append((lab, d, f, t))
        m = d[(f > 40) & (f < 16000)].max()
        vmax = m if vmax is None else max(vmax, m)
    for i, (lab, d, f, t) in enumerate(specs):
        ax = fig.add_subplot(gs[0, i])
        # each panel normalised to its own peak: we compare SHAPE, the level is matched separately
        dd = d - d[(f > 40) & (f < 16000)].max()
        ax.pcolormesh(t, f, dd, vmin=-90, vmax=0, cmap="magma", shading="auto")
        ax.set_yscale("log")
        ax.set_ylim(40, 16000)
        ax.set_xlim(0, min(3.0, t[-1] if len(t) else 3.0))
        ax.set_title(lab, fontsize=10)
        ax.set_xlabel("s")
        if i == 0:
            ax.set_ylabel("Hz (log)")
    ax = fig.add_subplot(gs[1, : max(1, len(panels) // 2)])
    for lab, d in desc_pairs:
        h = np.array(d["harm_db"][:24])
        if not len(h):
            continue
        ax.plot(np.arange(1, len(h) + 1), h - h.max(), marker="o", ms=3, label=lab)
    ax.set_xlabel("harmonic number")
    ax.set_ylabel("dB re strongest")
    ax.set_ylim(-80, 3)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax2 = fig.add_subplot(gs[1, max(1, len(panels) // 2):])
    fc = 1000 * 2 ** (np.arange(-14, 15) / 3)
    for lab, d in desc_pairs:
        b = np.array(d["bands_db"])
        ax2.semilogx(fc, b, marker=".", label=lab)
    ax2.set_xlabel("third-octave band centre (Hz)")
    ax2.set_ylabel("dB re total")
    ax2.set_ylim(-70, 0)
    ax2.grid(alpha=0.3)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


results = []
for c in CANDS:
    candidate_started = time.monotonic()
    try:
        ours = prep_candidate(c)
    except Exception as e:  # noqa: BLE001
        result = dict(id=c["id"], grader_version=GRADER_VERSION, error=str(e))
        results.append(result)
        RECEIPTS.append(receipt(c, result, candidate_started))
        continue
    kind = c.get("kind", "sustained")
    band_split_hz = c.get("band_split_hz")
    if c.get("chord"):
        real = real_chord(c["chord"], c["dyn"], c["inst"], c.get("stagger_s", 0.08))
        f0_hint = hz_of(midi_of(c["chord"][0]))
        dr = describe(highpass(real, f0_hint), kind, f0_hint, chord=True, band_split_hz=band_split_hz)
        do = describe(ours, kind, f0_hint, chord=True, band_split_hz=band_split_hz)
        target = None
        rl = {"other_dynamic": [], "neighbour": [], "other_instruments": {}}
        for shift in (1, 2, -1):
            nb = [__import__("instrgrade").name_of(midi_of(n) + shift) for n in c["chord"]]
            rc = real_chord(nb, c["dyn"], c["inst"], c.get("stagger_s", 0.08))
            if rc is not None:
                rl["neighbour"].append(distance(dr, describe(highpass(rc, f0_hint), kind, f0_hint, chord=True))["spec_db"])
    else:
        nm = midi_of(c["note"])
        cand_t = find(LIB, c["inst"], nm, dyn=c.get("dyn"), art=c.get("art"))
        if not cand_t:
            cand_t = find(LIB, c["inst"], nm, art=c.get("art"), tol=1)
        if not cand_t:
            result = dict(id=c["id"], grader_version=GRADER_VERSION,
                          error=f"IRPGRADE_REFERENCE_MISSING: {c['inst']} {c['note']}")
            results.append(result)
            RECEIPTS.append(receipt(c, result, candidate_started))
            continue
        target = cand_t[0]
        real = ref_audio(target, WAVDIR)
        inh = c["inst"] in ("xylophone", "marimba")
        dr = describe(real, kind, hz_of(nm), inharmonic=inh, band_split_hz=band_split_hz)
        do = describe(ours, kind, hz_of(nm), inharmonic=inh, band_split_hz=band_split_hz)
        rl = rulers(LIB, dict(target, desc=dr))
    if do is None:
        result = dict(id=c["id"], grader_version=GRADER_VERSION, error="no pitch found in candidate")
        results.append(result)
        RECEIPTS.append(receipt(c, result, candidate_started, target))
        continue
    dist = distance(do, dr)
    near = nearest_instruments(LIB, do, do["f0"], kind) if not c.get("chord") else []
    same_ruler = min(rl.get("other_dynamic") or [float("nan")]) if rl else float("nan")
    baselines = null_baselines(dr, cand_spec_db=dist["spec_db"], dyn_ruler_db=same_ruler)
    q3_summary = dict(
        attack_erb_distance_db=dict(value=dist.get("attack_erb_distance_db"), unit="dB"),
        partial_onset_order_vs_reference=dict(value=dist.get("partial_onset_order_tau"), unit="kendall_tau"),
        ours=do.get("q3"),
        real=dr.get("q3"),
    )
    desc_summary = dict(
        erb_distance_db=dict(
            attack_db=dist.get("attack_erb_distance_db"),
            sustain_db=dist.get("sustain_erb_distance_db"),
            release_db=dist.get("release_erb_distance_db"),
            unit="dB",
        ),
        ours=do.get("descriptors"),
        real=dr.get("descriptors"),
    )
    res = dict(
        id=c["id"],
        grader_version=GRADER_VERSION,
        label=c["label"],
        inst=c["inst"],
        target=(target["file"] + "@" + target["name"] + "." + target["dyn"]) if target else "chord",
        windows=do.get("windows"),
        band_split_hz=band_split_hz,
        ours_f0=do["f0"],
        ours_cents_vs_et=cents(do["f0"], hz_of(round(69 + 12 * math.log2(do["f0"] / 440)))),
        ours=do,
        real=dr,
        distance=dist,
        q3=q3_summary,
        descriptors=desc_summary,
        baselines=baselines,
        rulers=rl,
        nearest=near[:6],
    )
    if rl and rl.get("neighbour"):
        res["ratio_to_neighbour_ruler"] = dist["spec_db"] / max(1e-9, float(np.median(rl["neighbour"])))
    if rl and rl.get("other_dynamic"):
        res["ratio_to_dynamic_ruler"] = dist["spec_db"] / max(1e-9, float(np.median(rl["other_dynamic"])))
    if rl and rl.get("other_instruments"):
        res["closest_other_instrument_db"] = min(rl["other_instruments"].values())
    results.append(res)
    # level-matched listening pair (equal RMS over the first second, shared -1 dBFS peak)
    a, b = real[: int(4 * SR)], ours[: int(4 * SR)]
    ga, gb = 1 / max(rms(a[:SR]), 1e-12), 1 / max(rms(b[:SR]), 1e-12)
    a, b = a * ga, b * gb
    pk = max(np.abs(a).max(), np.abs(b).max())
    a, b = a / pk * 0.89, b / pk * 0.89
    write_wav(OUT / "pairs" / f"{c['id']}--real.wav", a)
    write_wav(OUT / "pairs" / f"{c['id']}--frankensim.wav", b)
    tlab = res["target"]
    figure(c["id"], f"{c['label']}  vs  real {c['inst']} ({tlab})   spec distance {dist['spec_db']:.1f} dB",
           [(f"REAL {c['inst']} {tlab}", real), (f"FRANKENSIM: {c['label']}", ours)],
           [("real", dr), ("frankensim", do)], OUT / "fig" / f"{c['id']}.png")
    RECEIPTS.append(receipt(c, res, candidate_started, target))
    print(json.dumps(dict(id=c["id"], grader_version=GRADER_VERSION, f0=round(do["f0"], 2), dist=round(dist["spec_db"], 1), env=round(dist["env_db"], 1),
                          harm=round(dist["harm_db"], 1), ruler_dyn=[round(v, 1) for v in rl.get("other_dynamic", [])] if rl else None, ruler_nb=[round(v, 1) for v in rl.get("neighbour", [])][:4] if rl else None, peaks=[round(v, 2) for v in do.get("peak_ratios", [])], peaks_real=[round(v, 2) for v in dr.get("peak_ratios", [])],
                          others={k: round(v, 1) for k, v in (rl.get("other_instruments") or {}).items()} if rl else None,
                          nearest=[(n[0], round(n[1], 1)) for n in near[:4]],
                          hnr=(round(do["hnr_db"], 1), round(dr["hnr_db"], 1)),
                          cent_rel=(round(do["centroid_rel"], 2), round(dr["centroid_rel"], 2)),
                          attack_ms=(round(do["attack_ms"]), round(dr["attack_ms"])),
                          f0_sd=(round(do["f0_sd_cents"], 1) if do["f0_sd_cents"] == do["f0_sd_cents"] else None,
                                 round(dr["f0_sd_cents"], 1) if dr["f0_sd_cents"] == dr["f0_sd_cents"] else None))))
(OUT / "results.json").write_text(json.dumps(results, indent=1, default=float))
failed = not RECEIPTS or any(r["verdict"] in ("fail", "error") for r in RECEIPTS)
RECEIPTS.append(dict(CONTEXT, event="terminal", verdict="fail" if failed else "pass", candidates=len(results),
                     elapsed_s=time.monotonic() - STARTED))
(OUT / "receipt.jsonl").write_text("".join(json.dumps(r, default=float) + "\n" for r in RECEIPTS))
if failed:
    raise SystemExit(1)
