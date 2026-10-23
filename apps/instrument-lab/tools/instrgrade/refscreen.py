#!/usr/bin/env python3
"""refscreen: cut a real recording to one note with its natural decay, screen it, and pin it by hash (pexfp.365).

  cut      SRC OUT.wav [--note A4] [--kind decaying|sustained] [--hint S] [--until S] [--channels 0,1]
           find the note's onset (the grader's detector), keep LEAD_S before it and everything up to where the note
           has decayed to its noise floor, write those source samples unchanged (source rate and channels, float32),
           then screen the clip; prints the row as JSON
  screen   CLIP.wav [--note A4] [--kind ...]       the six checks below on a clip already cut
  verify   MANIFEST CLIP_ROOT                      every clip of a references manifest (data/instruments/references-screened.json):
                                                   its sha256 and its screen verdict, recomputed from CLIP_ROOT/<family>/<id>.wav
  coverage MANIFEST                                per family: VALID pitches x dynamics, and what failed and why

The screen is a SANITY check that a clip is one clean note of the pitch its source names. It is not a realism grader
and nothing is tuned to these recordings. Its thresholds were fixed here before any recording was screened:
  onset   decompose.align finds the note's onset in the clip (the grader's two-criterion detector, quiet lead-in)
  snr     the note's loudest 20 ms is at least SNR_MIN_DB over the clip's noise, measured where there is no note (noise_floor):
          its lead-in, else a stationary tail, else NOT MEASURABLE and the check fails (a library sample cut at its own attack
          that fades out has no silence to read the noise in; a 5th percentile of its frames lies inside the note). The peak,
          not the grader's steady level: the median of a decaying note's frames sits part-way down its decay, so a
          clean struck note 60 dB over its noise read 33 dB that way (test_refscreen, before any recording was screened)
  pitch   pitched notes: the harmonic-sum f0 (decompose.f0_comb) and the median YIN f0 both lie within
          PITCH_TOL_CENTS of the named note. Unpitched notes (drums, cymbals) name no pitch: not applicable
  single  no second attack: from 50 ms after the onset to the end, no 10 ms frame is REATTACK_DB over the loudest frame of
          the 20-300 ms before it. A strike jumps over the decaying envelope; a beating dip between unison strings (a UIowa piano
          ff Ab6 notch 20 dB deep in 50 ms) only climbs back to it, which a rise measured from the dip took for a second note
  clip    no run of 3 or more samples at full scale in any channel of the source's own samples
  decay   the clip ends in its own decay: its last 20 ms are under max(noise + 10 dB, peak - 40 dB)
VALID = all six. Audio never enters the repo (it is PUBLIC); a manifest row carries the source URL and sha256, the
licence, the note, the dynamic, the trim times in the source and the clip's sha256.

On trj: uv run --with numpy python3 -I refscreen.py ... (ffmpeg decodes every source format at its own rate).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import struct
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import decompose as dc  # noqa: E402
from instrgrade import SR, cents, highpass, hz_of, midi_of, yin_track  # noqa: E402

SCREEN_VERSION = "refscreen-4"
SNR_MIN_DB = 40.0
PITCH_TOL_CENTS = 50.0
REATTACK_DB = 10.0
FULL_SCALE = 0.9999
LEAD_S = 0.10  # kept before the onset
LEAD_MIN_S = 0.04  # lead-in, ending LEAD_GUARD_S before the onset, that measures the floor
LEAD_GUARD_S = 0.02  # the FFT resample rings ahead of a sharp attack for some ms: that is not the recording's noise
LEAD_EDGE_S = 0.02  # nor is the clip's first 20 ms, where the resample meets the step from its padding
TAIL_FRAMES = 10  # 200 ms of 20 ms frames at the clip's end ...
TAIL_SPREAD_DB = 6.0  # ... read as noise only if stationary: their 90th - 10th percentile level within this
TAIL_UNDER_PEAK_DB = 20.0  # ... and this far under the note's peak: a loop sample cut mid-sustain ends ON the note
SILENT_DB = -150.0  # a frame this quiet is digital silence (a fade's end, zero padding), not noise
TAIL_OVER_FLOOR_DB = 6.0  # the note has decayed when its 20 ms level reaches this far over the floor
TAIL_KEEP_S = 0.05
CHECKS = ("onset", "snr", "pitch", "single", "clip", "decay")


# ------------------------------------------------------------------------------------------------ audio
def decode(path):
    """(samples [n, channels] float64, rate) of any source, at its own rate, by ffmpeg."""
    probe = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "a:0", "-show_entries", "stream=sample_rate,channels",
                            "-of", "json", str(path)], capture_output=True, text=True, check=True, stdin=subprocess.DEVNULL)
    st = json.loads(probe.stdout)["streams"][0]
    sr, ch = int(st["sample_rate"]), int(st["channels"])
    raw = subprocess.run(["ffmpeg", "-nostdin", "-v", "error", "-i", str(path), "-map", "0:a:0", "-f", "f32le", "-c:a", "pcm_f32le", "-"],
                         capture_output=True, check=True, stdin=subprocess.DEVNULL).stdout
    x = np.frombuffer(raw[: len(raw) // (4 * ch) * 4 * ch], "<f4").astype(np.float64).reshape(-1, ch)
    return x, sr


def to_sr(x, sr):
    """One channel at the grader's rate, by the band-limited FFT resample instrgrade.load uses, on the signal zero-padded to
    q 2^k samples (q = sr / gcd(sr, SR)) so that both transforms have fast lengths: unpadded, a 30 s 44.1 kHz file's length
    has large prime factors and one resample took 34 s."""
    if sr == SR:
        return np.asarray(x, np.float64)
    q = sr // math.gcd(sr, SR)
    n = q * (1 << max(0, math.ceil(math.log2(max(1, -(-len(x) // q))))))
    n_out = n * SR // sr
    X = np.fft.rfft(x, n)
    Y = np.zeros(n_out // 2 + 1, complex)
    m = min(len(X), len(Y))
    Y[:m] = X[:m]
    return np.fft.irfft(Y, n_out)[: int(round(len(x) * SR / sr))] * (n_out / n)


def write_wav_f32(path, x, sr):
    """[n, channels] float32 WAV (format 3), so 16- and 24-bit source samples are kept exactly."""
    y = np.ascontiguousarray(np.asarray(x, "<f4"))
    ch = y.shape[1]
    data = y.tobytes()
    Path(path).write_bytes(b"RIFF" + struct.pack("<I", 36 + len(data)) + b"WAVEfmt " +
                           struct.pack("<IHHIIHH", 16, 3, ch, sr, sr * 4 * ch, 4 * ch, 32) + b"data" + struct.pack("<I", len(data)) + data)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def level_db(y, hop):
    """Mean-square level of consecutive `hop`-sample frames [dB re full scale]."""
    n = len(y) // hop
    if n < 1:
        return np.zeros(0)
    return 10 * np.log10(np.mean(y[: n * hop].reshape(n, hop) ** 2, axis=1) + 1e-18)


def clipped_runs(x, run=3):
    """Runs of `run` or more consecutive samples at full scale, over all channels."""
    total = 0
    for c in range(x.shape[1]):
        hot = np.abs(x[:, c]) >= FULL_SCALE
        if not hot.any():
            continue
        edges = np.diff(np.concatenate([[0], hot.astype(np.int8), [0]]))
        starts, ends = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)
        total += int(np.sum(ends - starts >= run))
    return total


# ----------------------------------------------------------------------------------------------- screen
def noise_floor(y, hp, named=None):
    """(noise [dB re full scale] or None, method) of a clip (y at the grader's rate, hp its high-passed signal), measured
    where there is no note.

    lead-in: the samples from LEAD_EDGE_S to LEAD_GUARD_S before the note's first 10 ms frame within STEADY_DROP_DB of its steady
    level (the grader's criterion B, which needs no floor; on 10 ms frames, where noise varies by tenths of a dB, not on the
    grader's 2 ms envelope, whose noise frames 30 dB under a note reach 20 dB under it), when there are LEAD_MIN_S of them,
    high-passed ON THEIR OWN: the
    zero-phase FFT high-pass of the whole clip rings 40-60 ms ahead of a sharp attack (a burst 70 dB over its noise read 64).
    tail: else the median of the last TAIL_FRAMES 20 ms frames above digital silence, when they are stationary (90th - 10th
    percentile within TAIL_SPREAD_DB) and TAIL_UNDER_PEAK_DB under the clip's loudest frame: noise plus whatever is left of the
    note, so an upper bound on the noise (a FreePats loop sample, cut mid-sustain, ended stationary 0.6 dB under its own peak,
    which is the note, not noise). none: neither."""
    e = level_db(hp, SR // 100)
    if not len(e):
        return None, "none"
    above = np.flatnonzero(e > dc.steady_level_db(e) - dc.STEADY_DROP_DB)
    lead = y[int(LEAD_EDGE_S * SR): max(0, int(round((above[0] * 0.01 - LEAD_GUARD_S) * SR)))] if len(above) else y[:0]
    if len(lead) >= LEAD_MIN_S * SR:
        return float(10 * np.log10(np.mean(highpass(lead, named) ** 2) + 1e-18)), "lead-in"
    live = level_db(hp, int(0.02 * SR))
    live = live[live > SILENT_DB][-TAIL_FRAMES:]
    top = float(level_db(hp, int(0.02 * SR)).max()) if len(hp) >= int(0.02 * SR) else -120.0
    if (len(live) == TAIL_FRAMES and np.percentile(live, 90) - np.percentile(live, 10) <= TAIL_SPREAD_DB
            and np.median(live) <= top - TAIL_UNDER_PEAK_DB):
        return float(np.median(live)), "tail"
    return None, "none"


def screen_signal(x_native, sr, note=None, kind="decaying"):
    """The six checks on one clip: (measured values, failed check names with their reasons)."""
    named = hz_of(midi_of(note)) if note else None
    y = to_sr(x_native[:, 0], sr)
    hp = highpass(y, named)
    noise, how = noise_floor(y, hp, named)
    # with no measurable noise the onset is criterion B's alone: the grader's convention for a signal that begins at its onset
    floor = noise if noise is not None else -120.0
    m = dict(screen=SCREEN_VERSION, noise_db=None if noise is None else round(noise, 2), noise_method=how,
             floor_p5_db=round(dc.floor_of(hp), 2), clipped_runs=clipped_runs(x_native))
    fails = []
    pre = dc.preroll_s(named)
    try:
        cut, info = dc.align(y, named, floor)
        m["onset_s"] = round(info["onset_s"], 4)
        body = hp[max(0, int(round((info["onset_s"] - pre) * SR))):]  # the whole note from the cut's start, not the grader's 3 s
    except dc.OnsetMisaligned as e:
        fails.append(("onset", str(e)))
        cut = body = hp[int(LEAD_S * SR):]
    steady = dc.steady_level_db(dc.env_rms_db(cut, win=dc.env_window(named))[0])
    peak = float(level_db(body, int(0.02 * SR)).max()) if len(body) >= int(0.02 * SR) else -120.0
    m.update(steady_db=round(steady, 2), peak_db=round(peak, 2), snr_db=None if noise is None else round(peak - noise, 2))
    if noise is None:
        fails.append(("snr", "the noise is not measurable: no lead-in before the onset and no stationary tail "
                             f"(the clip begins at its attack and its last {TAIL_FRAMES} frames above digital silence are a decay or a fade)"))
    elif peak - noise < SNR_MIN_DB:
        fails.append(("snr", f"peak level {peak:.1f} dB is {peak - noise:.1f} dB over the noise {noise:.1f} dB ({how}; < {SNR_MIN_DB:.0f})"))
    if named:
        f_comb = dc.f0_comb(cut, named, kind)
        t0, t1 = (0.05, 0.6) if kind == "decaying" else (0.15, 1.5)
        seg = cut[int(t0 * SR):int(t1 * SR)]
        tr = yin_track(seg, fmin=max(20.0, named / 2.2), fmax=min(5000.0, named * 2.2))
        good = tr[(tr[:, 2] > 0.8) & np.isfinite(tr[:, 1])] if len(tr) else tr
        f_yin = float(np.median(good[:, 1])) if len(good) else float("nan")
        c_comb, c_yin = cents(f_comb, named), (cents(f_yin, named) if math.isfinite(f_yin) else float("nan"))
        m.update(f0_comb_hz=round(f_comb, 3), f0_yin_hz=round(f_yin, 3) if math.isfinite(f_yin) else None,
                 cents_comb=round(c_comb, 1), cents_yin=round(c_yin, 1) if math.isfinite(c_yin) else None, yin_frames=int(len(good)))
        if not (abs(c_comb) <= PITCH_TOL_CENTS and math.isfinite(c_yin) and abs(c_yin) <= PITCH_TOL_CENTS):
            fails.append(("pitch", f"named {note} {named:.2f} Hz; harmonic-sum f0 {f_comb:.2f} Hz ({c_comb:+.0f} c), "
                                   f"YIN {f_yin:.2f} Hz ({c_yin:+.0f} c) on {len(good)} periodic frames; tolerance {PITCH_TOL_CENTS:.0f} c"))
    e = level_db(body, SR // 100)
    on = int(round((pre + 0.05) * 100))
    rise, at = -120.0, None
    for j in range(on + 3, len(e)):
        if e[j] > floor + 10.0:
            r = float(e[j] - e[max(on, j - 30):j - 2].max())
            if r > rise:
                rise, at = r, j
    m["reattack_db"] = round(max(rise, 0.0), 2)
    if rise >= REATTACK_DB:
        fails.append(("single", f"{(at - on) / 100 + 0.05:.2f} s after the onset a 10 ms frame is {rise:.1f} dB over the loudest of the "
                                f"300 ms before it: a second note or strike"))
    if m["clipped_runs"]:
        fails.append(("clip", f"{m['clipped_runs']} runs of 3 or more samples at full scale"))
    tail = float(level_db(hp[-int(0.02 * SR):], int(0.02 * SR))[0]) if len(hp) >= int(0.02 * SR) else 0.0
    m["tail_db"] = round(tail, 2)
    if tail > max(floor + 10.0, peak - 40.0):
        fails.append(("decay", f"the clip's last 20 ms are at {tail:.1f} dB, over max(floor + 10, peak - 40) = "
                               f"{max(floor + 10.0, peak - 40.0):.1f} dB: it stops before the note has decayed"))
    return m, fails


def verdict(m, fails):
    return dict(m, status="INVALID" if fails else "VALID", failed=[f"{k}: {why}" for k, why in fails])


# ------------------------------------------------------------------------------------------------- cut
def find_cut(x, sr, note=None, hint=None, until=None):
    """(start, end) in seconds of the source: LEAD_S before the note's onset, to TAIL_KEEP_S after the first 20 ms frame
    past the note's peak at TAIL_OVER_FLOOR_DB over the noise floor, never past `until`. With `hint` (a time inside the
    note) the onset is the last rise before it, as the grader finds a note inside a longer recording."""
    named = hz_of(midi_of(note)) if note else None
    t_end = len(x) / sr if until is None else min(until, len(x) / sr)
    w0 = 0.0 if hint is None else max(0.0, hint - 2.0)
    y = to_sr(x[int(w0 * sr):int(t_end * sr), 0], sr)
    hp = highpass(y, named)
    floor = dc.floor_of_file(y, named) if hint is None else dc.floor_of(hp)
    try:
        ta, tb = dc.onset_times(hp, floor, None if hint is None else hint - w0, named)
    except dc.OnsetMisaligned as e:  # no onset to cut at: keep the whole window, and the screen says why
        return w0, t_end, dict(floor_db=round(floor, 2), onset_error=str(e))
    t_on = max(ta, tb)
    hop = int(0.02 * SR)
    e = level_db(hp[int(t_on * SR):], hop)
    pk = int(np.argmax(e[: int(3.0 / 0.02)])) if len(e) else 0
    quiet = np.flatnonzero(e[pk:] <= floor + TAIL_OVER_FLOOR_DB)
    end = t_on + ((pk + quiet[0] + 1) * hop / SR + TAIL_KEEP_S if len(quiet) else len(e) * hop / SR)
    return max(0.0, w0 + t_on - LEAD_S), min(t_end, w0 + end), dict(floor_db=round(floor, 2), onset_in_source_s=round(w0 + t_on, 4))


def cut(src, out, note=None, kind="decaying", hint=None, until=None, channels=None):
    x, sr = decode(src)
    if channels is not None:
        x = x[:, channels]
    t0, t1, info = find_cut(x, sr, note, hint, until)
    i0, i1 = int(round(t0 * sr)), int(round(t1 * sr))
    clip = x[i0:i1]
    write_wav_f32(out, clip, sr)
    m, fails = screen_signal(clip, sr, note, kind)
    return dict(note=note, kind=kind, sample_rate_hz=sr, channels=int(clip.shape[1]), trim_start_s=round(i0 / sr, 6),
                trim_end_s=round(i1 / sr, 6), clip_frames=int(i1 - i0), clip_sha256=sha256(out),
                source_sha256=sha256(src), source_bytes=Path(src).stat().st_size, **{"cut": info}, **verdict(m, fails))


# -------------------------------------------------------------------------------------------- manifest
def coverage(manifest):
    """Per family: VALID pitches (or pieces, for unpitched rows) by dynamic, and the failure reasons counted."""
    fam = defaultdict(lambda: dict(valid=defaultdict(set), rows=0, n_valid=0, failed=defaultdict(int)))
    for r in manifest["rows"]:
        f = fam[r["family"]]
        f["rows"] += 1
        what = r["note"] or r.get("piece") or "unpitched"
        if r["status"] == "VALID":
            f["n_valid"] += 1
            f["valid"][r["dynamic"]].add(what)
        for why in r["failed"]:
            f["failed"][why.split(":")[0]] += 1
    out = {}
    for name, f in sorted(fam.items()):
        dyns = {d: sorted(v, key=lambda s: (midi_of(s) if s[0] in "ABCDEFG" and s[-1].isdigit() else 999, s)) for d, v in sorted(f["valid"].items())}
        best = max((len(v) for v in dyns.values()), default=0)
        out[name] = dict(rows=f["rows"], valid=f["n_valid"], dynamics_with_valid=len(dyns), most_pitches_in_one_dynamic=best,
                         distinct_valid=len(set().union(*dyns.values())) if dyns else 0, done_when=best >= 3,
                         by_dynamic=dyns, failed_checks=dict(f["failed"]))
    return out


def cmd_cut(a):
    ch = [int(c) for c in a.channels.split(",")] if a.channels else None
    print(json.dumps(cut(a.src, a.out, a.note, a.kind, a.hint, a.until, ch)))


def cmd_screen(a):
    x, sr = decode(a.clip)
    print(json.dumps(verdict(*screen_signal(x, sr, a.note, a.kind))))


def manifest_clips(man):
    """Every clip of a references manifest: those of each file and those of each archive member."""
    for f in man["files"]:
        yield from f.get("clips", [])
        for m in f.get("members", []):
            yield from m.get("clips", [])


def verify_clip(job):
    """(clip id, problems) for one manifest clip: its file under root, its sha256, its verdict recomputed."""
    root, r = job
    p = Path(root) / r["family"] / f"{r['id']}.wav"
    problems = []
    if not p.exists():
        problems.append("missing")
    elif sha256(p) != r["clip_sha256"]:
        problems.append("sha256 differs")
    else:
        x, sr = decode(p)
        v = verdict(*screen_signal(x, sr, r.get("note"), r["kind"]))  # an unpitched clip names no note
        if v["status"] != r["status"]:
            problems.append(f"status {v['status']} != manifest {r['status']}")
    return r["id"], problems


def cmd_verify(a):
    from concurrent.futures import ProcessPoolExecutor
    man = json.loads(Path(a.manifest).read_text())
    clips = list(manifest_clips(man))
    jobs = [(a.root, r) for r in clips if "clip_sha256" in r]  # a source that could not be cut claims nothing
    jobs_n = getattr(a, "jobs", 1) or 1
    if jobs_n > 1:
        with ProcessPoolExecutor(jobs_n) as pool:
            results = list(pool.map(verify_clip, jobs, chunksize=4))
    else:
        results = [verify_clip(j) for j in jobs]
    bad = 0
    for cid, problems in results:
        bad += bool(problems)
        print(cid, "ok" if not problems else "; ".join(problems))
    print(json.dumps(dict(clips=len(clips), checked=len(jobs), mismatched=bad)))
    return 1 if bad else 0


def cmd_coverage(a):
    man = json.loads(Path(a.manifest).read_text())
    print(json.dumps(coverage(dict(rows=[dict(c, note=c.get("note"), failed=c.get("failed", [])) for c in manifest_clips(man)])), indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("cut")
    s.add_argument("src"), s.add_argument("out")
    s.add_argument("--note"), s.add_argument("--kind", default="decaying", choices=("decaying", "sustained"))
    s.add_argument("--hint", type=float), s.add_argument("--until", type=float), s.add_argument("--channels")
    s.set_defaults(f=cmd_cut)
    s = sp.add_parser("screen")
    s.add_argument("clip"), s.add_argument("--note"), s.add_argument("--kind", default="decaying", choices=("decaying", "sustained"))
    s.set_defaults(f=cmd_screen)
    s = sp.add_parser("verify")
    s.add_argument("manifest"), s.add_argument("root"), s.add_argument("--jobs", type=int, default=1)
    s.set_defaults(f=cmd_verify)
    s = sp.add_parser("coverage")
    s.add_argument("manifest")
    s.set_defaults(f=cmd_coverage)
    a = p.parse_args()
    sys.exit(a.f(a) or 0)


if __name__ == "__main__":
    main()
