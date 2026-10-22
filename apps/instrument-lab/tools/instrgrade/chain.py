#!/usr/bin/env python3
"""chain: put a render through the recording chain of a reference session, built from measured chain facts only (M2.6-W2).

The UIowa recordings were made through a chain: microphone, room, preamp/mixer, recorder, then 16-bit/44.1 kHz files, which
convert_refs.py turns into the grader's 48 kHz floats. A render that skips the chain differs from them in ways no instrument
physics can fix. This module gives a render what the chain measurably adds, and flags everything it cannot measure.

  floor       ACTUAL recorded floor samples of the session: its ambient file(s) where the page lists one, plus the gaps
              between the notes of its graded files, with breath, key and decay excluded by GAP_RULE (declared below,
              before any grading; never tuned). Blocks of BLOCK_S are drawn from the clean stretches by a counter stream
              keyed by (description, chain id, block), and joined by equal-power crossfades. Nothing is synthesised.
  grid        the render (already level-matched to its reference) is moved onto the recording's 44.1 kHz grid, the floor is
              added there, the sum is quantized to the graded files' bit depth (as probed: round half to even, no dither), back to 48 kHz
              by the same ffmpeg aresample that convert_refs.py applies to every reference.
  flags       every stage the session pages do not document is a FLAG in the receipt, never an invented stage: the capsule's
              response and pattern (M2.6-W1, not applied here), the A/D anti-alias filter, dither, mixer EQ, the 2011+
              downsampler, and the chamber's low cutoff (descriptors below LOW_FLAG_HZ stay flagged until it is sourced).

The facts come from data/instruments/references.json (the chain and page records, the files and their sha256). Audio is never
committed: the files are read from FRANKENSIM_REFERENCE_DIR (fetch-references.sh), and a missing or altered file is a refusal.

NOT VALIDATED: the gap floor. Only one session documents an ambient file (B-flat clarinet, 2.24 s), and there the gap floor does
not match it (bead pexfp.117, 2026-10-22; 248 s of gaps against the ambient):
- 2-4 dB lower below 130 Hz, where the ambient holds more rumble;
- 3-11 dB higher at 1-5 kHz, broadband with no lines (the player breathing between notes?), which a level rule cannot separate
  from the floor;
- 17 of 31 third-octave bands differ by more than 1 dB beyond their 95 % interval.
Every chain whose floor uses gaps therefore carries the flag gap_floor_unvalidated.

  chain.py floor PAGE --manifest references.json --refdir DIR --lib LIBRARY.json [--out floor.wav]
  chain.py apply IN.wav OUT.wav PAGE --description ID --manifest ... --refdir ... --lib ...
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from instrgrade import SR, THIRD_OCT, load, segment_scale, write_wav  # noqa: E402

CHAIN_VERSION = "chain-1"
FS_REC = 44100  # the usual rate of the graded files; each session's own rate is probed (the 2011 guitar .mono files are 96 kHz)

# GAP_RULE, declared before any grading and never tuned against graded notes. A gap is a stretch of a graded file that no
# sounding span covers: neither a library note nor any span instrgrade.segment_scale finds in the file (a note the library did
# not verify is still a note), with LEAD_S before each span (breath intake, key preparation) and TAIL_S after it (its decay and
# the chamber's tail) given to the span. Its 50 ms frames are judged against the session's own floor: a frame louder than the
# median of all gap and ambient frames by EVENT_MADS robust standard deviations (1.4826 MAD, in dB) is an event, and GUARD_S on
# either side of it goes with it. What is left in runs of at least MIN_RUN_S is floor.
# Version 2 (2026-10-22): version 1 took library notes only. Against the B-flat clarinet session's ambient file, notes the
# library had not verified leaked into its gaps (a 660.8 Hz line 35 dB over the ambient, E5). Found on the measured ambient
# file, not on any graded note.
GAP_RULE = dict(version=2, lead_s=0.15, tail_s=0.5, frame_s=0.05, event_mads=4.0, guard_s=0.05, min_run_s=0.25)
BLOCK_S = 0.5  # floor is laid down in blocks this long, each drawn by the counter stream
XFADE_S = 0.01  # equal-power (sqrt-Hann) crossfade between blocks: uncorrelated noise keeps its level across a join
LOW_FLAG_HZ = 100.0
FFMPEG = ["ffmpeg", "-v", "error", "-nostdin"]


class ChainError(ValueError):
    """IRPCHAIN_*: a chain that cannot be built or applied from measured facts. Never defaulted around."""


@dataclass
class Stretch:
    file: str  # path within the reference dir, as references.json names it
    sha256: str
    start: int  # samples at the session's rate
    end: int
    source: str  # "ambient" or "gap"


@dataclass
class Chain:
    id: str
    page: str
    record: dict
    stretches: list[Stretch]
    flags: list[dict]
    bits: int  # the graded files' bits per sample, as probed
    fs: int = FS_REC  # their sample rate, as probed (after convert_refs.py's sample-rate-header repair)
    sources: list[dict] = field(default_factory=list)  # every file read: path, sha256, role
    _audio: dict = field(default_factory=dict, repr=False)

    # ------------------------------------------------------------------ building
    @classmethod
    def from_references(cls, manifest: Path, page: str, refdir: Path, library: list[dict]):
        """The chain of session `page` (an MIS page in references.json): its chain record, its ambient file(s), and the gaps of its
        graded files located with the library's note spans (decompose's library: file stem.wav, start, end in seconds)."""
        m = json.loads(Path(manifest).read_text())
        prec = next((p for p in m["pages"] if p["page"] == page), None)
        if prec is None or not prec.get("chain"):
            raise ChainError(f"IRPCHAIN_UNKNOWN_SESSION: {page} has no chain record in {manifest}")
        rec = next(c for c in m["chains"] if c["id"] == prec["chain"])
        files = {f["id"]: f for f in m["files"]}
        flags = chain_flags(rec, prec)
        audio, sources, stretches, fmt = {}, [], [], set()
        amb = prec.get("ambient")
        if amb and amb != "absent-on-page":
            f = files.get(amb)
            if f is None:
                raise ChainError(f"IRPCHAIN_AMBIENT_MISSING: {page} names ambient {amb}, which references.json does not list")
            x, bits, rate = read_pcm(refdir, f)
            fmt.add((bits, rate))
            audio[f["path"]] = x
            sources.append(dict(path=f["path"], sha256=f["sha256"], role="ambient"))
            stretches.append(Stretch(f["path"], f["sha256"], 0, len(x), "ambient"))
        by_stem = {Path(f["path"]).stem.replace(" ", ""): f for f in m["files"] if f["page"] == page and f["kind"] == "chromatic-scale" and f["gradable"]}
        spans = {}
        for e in library:
            stem = Path(e["file"]).stem
            if stem in by_stem:
                spans.setdefault(stem, []).append((e["start"], e["end"]))
        for stem in sorted(spans):
            f = by_stem[stem]
            x, bits, rate = read_pcm(refdir, f)
            fmt.add((bits, rate))
            audio[f["path"]] = x
            sources.append(dict(path=f["path"], sha256=f["sha256"], role="gaps"))
            sounding = [(a / rate, b / rate) for a, b in segment_scale(x)]  # every sounding span, verified as a note or not
            stretches += [Stretch(f["path"], f["sha256"], a, b, "gap") for a, b in gaps(len(x), spans[stem] + sounding, fs=rate)]
        if not stretches:
            raise ChainError(f"IRPCHAIN_NO_FLOOR: {page} has no ambient file and no graded file with library notes to take gaps from")
        if len(fmt) != 1:
            raise ChainError(f"IRPCHAIN_FORMAT: {page}'s files disagree on (bits, rate): {sorted(fmt)}")
        bits, fs = fmt.pop()
        clean = exclude_events(stretches, audio, fs=fs)
        if not any(s.end - s.start >= int((BLOCK_S + XFADE_S) * fs) for s in clean):
            raise ChainError(f"IRPCHAIN_FLOOR_TOO_SHORT: no clean floor stretch of {BLOCK_S + XFADE_S:.2f} s in {page}")
        if not amb or amb == "absent-on-page":
            flags.append(dict(flag="ambient_absent_on_page", note="the floor comes from inter-note gaps only"))
        if any(s.source == "gap" for s in clean):
            flags.append(dict(flag="gap_floor_unvalidated", note="on the one session with an ambient file the gap floor does not match it "
                              "(+3..+11 dB at 1-5 kHz, -2..-4 dB below 130 Hz; pexfp.117)"))
        if fs != FS_REC:
            flags.append(dict(flag="rate_differs_from_record", note=f"the graded files are {bits}-bit at {fs} Hz; the chain record says "
                              f"{rec.get('format')}"))
        return cls(id=f'{rec["id"]}@{page}', page=page, record=rec, stretches=clean, flags=flags, bits=bits, fs=fs, sources=sources,
                   _audio=audio)

    # ------------------------------------------------------------------ floor
    def floor(self, n: int, description: str) -> tuple[np.ndarray, list[dict]]:
        """`n` samples at the session's rate of recorded floor for `description`: blocks drawn by the counter stream, crossfaded."""
        blk, xf = int(BLOCK_S * self.fs), int(XFADE_S * self.fs)
        usable = [(s, s.end - s.start - blk - xf) for s in self.stretches if s.end - s.start - blk - xf >= 0]
        total = sum(u + 1 for _, u in usable)
        nb = max(1, int(math.ceil(n / blk)))
        out = np.zeros(nb * blk + xf)
        w = np.sin(0.5 * np.pi * (np.arange(xf) + 0.5) / xf)  # sqrt-Hann rise; its mirror falls; sin^2 + cos^2 = 1
        picks = []
        for b in range(nb):
            u = stream(description, self.id, b) % total
            for s, room in usable:
                if u <= room:
                    break
                u -= room + 1
            a = s.start + u
            seg = self._audio[s.file][a:a + blk + xf].copy()
            seg[:xf] *= w if b else 1.0
            seg[blk:] *= w[::-1]
            out[b * blk:b * blk + blk + xf] += seg
            picks.append(dict(block=b, file=s.file, start=int(a), source=s.source))
        return out[:n], picks

    # ------------------------------------------------------------------ apply
    def apply(self, y48: np.ndarray, description: str) -> tuple[np.ndarray, dict]:
        """A level-matched render at the grader's 48 kHz, as the session would have recorded it: on the 44.1 kHz grid, plus the
        session's floor, quantized to its files' bit depth, and back to 48 kHz by convert_refs.py's resampler. Returns (signal, receipt)."""
        x = resample(np.asarray(y48, dtype=np.float64), SR, self.fs)
        fl, picks = self.floor(len(x), description)
        q, clips = quantize(x + fl, self.bits)
        if clips:
            raise ChainError(f"IRPCHAIN_CLIPS: {clips} samples exceed 16-bit full scale after the floor is added: the render is "
                             "louder than the session could record at its reference's level")
        y = resample_pcm(q, self.bits, self.fs, SR)
        return y, self.receipt(description, fl, picks)

    def receipt(self, description, fl, picks):
        return dict(stage="chain", version=CHAIN_VERSION, chain_id=self.id, description=description,
                    ambient_sha256=[s["sha256"] for s in self.sources if s["role"] == "ambient"],
                    gap_sources=[dict(path=s["path"], sha256=s["sha256"]) for s in self.sources if s["role"] == "gaps"],
                    gap_rule=GAP_RULE, block_s=BLOCK_S, xfade_s=XFADE_S, blocks=picks,
                    floor_db_per_band=[round(v, 2) for v in third_octave_db(fl, self.fs)], bands_hz=[round(f, 1) for f in THIRD_OCT],
                    bits=self.bits, fs_hz=self.fs, quantize="round half to even, no dither", resampler=resampler_id(), flags=self.flags,
                    low_flag_hz=LOW_FLAG_HZ)


# ---------------------------------------------------------------------- the pieces
def chain_flags(rec: dict, page: dict) -> list[dict]:
    """Every chain stage the session pages leave undocumented, as a flag (never an invented stage)."""
    flags = [dict(flag="capsule_not_applied", note=f'{rec["mic"]} response and pattern are M2.6-W1; this chain adds no capsule'),
             dict(flag="antialias_filter_undocumented", note="the recorder's own filter is not documented; ffmpeg aresample puts the render on the 44.1 kHz grid"),
             dict(flag="dither_undocumented", note="no dither is applied")]
    if rec.get("low_cutoff_status") != "sourced":
        flags.append(dict(flag="low_cutoff_undocumented", note=f"descriptors below {LOW_FLAG_HZ:.0f} Hz are flagged until the chamber cutoff is sourced"))
    if any("mixer" in s.lower() for s in rec.get("signal_chain", [])):
        flags.append(dict(flag="mixer_eq_undocumented", note="the mixer's EQ setting is not documented; none is applied"))
    if "96" in rec.get("format", "") or rec.get("era") == "2011+":
        flags.append(dict(flag="downsampler_undocumented", note="the session's 96 to 44.1 kHz downsampler (and its dither) is not documented"))
    return flags


def read_pcm(refdir: Path, f: dict) -> tuple[np.ndarray, int, int]:
    """(samples as floats on the file's integer grid, its bits per sample, its sample rate) of a references.json file, after
    checking its sha256: its own rate, never resampled, with channel 0 and the sample-rate-header repair exactly as convert_refs.py
    does. Refuses a missing or altered file, and one that is not 16- or 24-bit PCM."""
    p = Path(refdir).joinpath(*f["path"].split("/"))
    if not p.is_file():
        raise ChainError(f"IRPCHAIN_AMBIENT_MISSING: {p} is not there (fetch it with data/instruments/fetch-references.sh)")
    h = hashlib.sha256(p.read_bytes()).hexdigest()
    if h != f["sha256"]:
        raise ChainError(f"IRPCHAIN_HASH_MISMATCH: {p} has sha256 {h}, references.json pins {f['sha256']}")
    r = subprocess.run(["ffprobe", "-v", "error", "-show_entries", "stream=sample_rate,bits_per_sample", "-of", "json", str(p)],
                       capture_output=True, text=True)
    st = json.loads(r.stdout or "{}").get("streams") or [{}]
    rate, bits = int(st[0].get("sample_rate", 0)), int(st[0].get("bits_per_sample", 0))
    defect = any(d["kind"] == "sample-rate-header" for d in f.get("defects", []))
    if bits not in (16, 24) or rate <= 0:
        raise ChainError(f"IRPCHAIN_FORMAT: {p} is {bits}-bit at {rate} Hz; the chain knows 16- and 24-bit PCM")
    true_rate = 96000 if defect else rate  # convert_refs.py: 96 kHz samples under a 44.1 kHz header
    filters = (["pan=mono|c0=c0"] if f.get("channels") == "stereo" else []) + ([f"asetrate={true_rate}"] if defect else [])
    r = subprocess.run(FFMPEG + ["-i", str(p)] + (["-af", ",".join(filters)] if filters else []) + ["-ac", "1", "-f", "f64le", "pipe:1"],
                       capture_output=True)
    if r.returncode:
        raise ChainError(f"IRPCHAIN_DECODE: {p}: {r.stderr.decode()[:200]}")
    return np.frombuffer(r.stdout, "<f8").astype(np.float64), bits, true_rate


def gaps(n: int, spans: list[tuple[float, float]], rule=GAP_RULE, fs: int = FS_REC) -> list[tuple[int, int]]:
    """Sample ranges of a file of `n` samples at `fs` that no note covers, each note given LEAD_S before and TAIL_S after."""
    taken = np.zeros(n, bool)
    for a, b in spans:
        taken[max(0, int((a - rule["lead_s"]) * fs)):min(n, int((b + rule["tail_s"]) * fs))] = True
    edges = np.flatnonzero(np.diff(np.concatenate([[1], taken.astype(np.int8), [1]])))
    return [(int(a), int(b)) for a, b in zip(edges[::2], edges[1::2]) if b > a]


def frame_levels(x: np.ndarray, frame: int) -> np.ndarray:
    k = len(x) // frame
    return 10 * np.log10(np.mean(x[: k * frame].reshape(k, frame) ** 2, axis=1) + 1e-30)


def exclude_events(stretches: list[Stretch], audio: dict, rule=GAP_RULE, fs: int = FS_REC) -> list[Stretch]:
    """GAP_RULE's event exclusion over all stretches of a session together: what is left in runs of MIN_RUN_S or longer."""
    fr = int(rule["frame_s"] * fs)
    levels = [frame_levels(audio[s.file][s.start:s.end], fr) for s in stretches]
    allv = np.concatenate([v for v in levels if len(v)]) if any(len(v) for v in levels) else np.zeros(0)
    if not len(allv):
        return []
    med = float(np.median(allv))
    sd = 1.4826 * float(np.median(np.abs(allv - med)))
    th = med + rule["event_mads"] * max(sd, 1e-6)
    guard = int(math.ceil(rule["guard_s"] / rule["frame_s"]))
    out = []
    for s, v in zip(stretches, levels):
        bad = v > th
        if guard and bad.any():
            bad = np.convolve(bad.astype(int), np.ones(2 * guard + 1, int), mode="same") > 0
        run = 0
        for i, b in enumerate(np.concatenate([bad, [True]])):
            if not b:
                run += 1
                continue
            if run * fr >= rule["min_run_s"] * fs:
                out.append(Stretch(s.file, s.sha256, s.start + (i - run) * fr, s.start + i * fr, s.source))
            run = 0
    return out


def stream(description: str, chain_id: str, block: int) -> int:
    """The counter stream: a 64-bit draw that depends only on (description, chain id, block), never on a thread or the clock."""
    return int.from_bytes(hashlib.sha256(f"{CHAIN_VERSION}|{description}|{chain_id}|{block}".encode()).digest()[:8], "little")


def quantize(x: np.ndarray, bits: int = 16) -> tuple[np.ndarray, int]:
    """Quantization to `bits`: round half to even on the 2^-(bits-1) grid of full scale, no dither. Returns (integer samples,
    samples that clipped)."""
    full = float(1 << (bits - 1))
    v = np.rint(np.asarray(x, dtype=np.float64) * full)
    clips = int(np.sum((v > full - 1) | (v < -full)))
    return np.clip(v, -full, full - 1).astype(np.int32), clips


def _ffmpeg(args: list[str], data: bytes) -> bytes:
    if shutil.which("ffmpeg") is None:
        raise ChainError("IRPCHAIN_NO_RESAMPLER: ffmpeg (convert_refs.py's resampler) is not installed")
    r = subprocess.run(FFMPEG + args, input=data, capture_output=True)
    if r.returncode:
        raise ChainError(f"IRPCHAIN_RESAMPLE: {r.stderr.decode()[:200]}")
    return r.stdout


def resample(x: np.ndarray, fs_in: int, fs_out: int) -> np.ndarray:
    """ffmpeg aresample (swr defaults: the resampler convert_refs.py uses) on a float signal."""
    out = _ffmpeg(["-f", "f64le", "-ar", str(fs_in), "-ac", "1", "-i", "pipe:0", "-af", f"aresample={fs_out}", "-f", "f64le", "pipe:1"],
                  np.asarray(x, dtype="<f8").tobytes())
    return np.frombuffer(out, "<f8").astype(np.float64)


def resample_pcm(q: np.ndarray, bits: int, fs_in: int, fs_out: int) -> np.ndarray:
    """Integer samples of `bits` to float at fs_out, as convert_refs.py converts a reference: the sample format an AIFF of that
    depth decodes to (s16, or s32 holding 24 bits), aresample, pcm_f32le."""
    q = np.asarray(q, dtype=np.int64)
    fmt, data = ("s16le", q.astype("<i2").tobytes()) if bits == 16 else ("s32le", (q << 8).astype("<i4").tobytes())
    out = _ffmpeg(["-f", fmt, "-ar", str(fs_in), "-ac", "1", "-i", "pipe:0", "-af", f"aresample={fs_out}", "-ac", "1", "-c:a", "pcm_f32le",
                   "-f", "f32le", "pipe:1"], data)
    return np.frombuffer(out, "<f4").astype(np.float64)


def resampler_id() -> str:
    try:
        v = subprocess.run(["ffmpeg", "-version"], capture_output=True, text=True).stdout.splitlines()[0]
    except (OSError, IndexError):
        v = "ffmpeg missing"
    return f"{v}; aresample (swr defaults), as convert_refs.py"


def third_octave_db(x: np.ndarray, fs: int) -> list[float]:
    """Power in each THIRD_OCT band [dB re full-scale power per band], Welch-averaged over 8192-point Hann frames."""
    n = 8192
    if len(x) < n:
        raise ChainError(f"IRPCHAIN_TOO_SHORT: {len(x)} samples, a third-octave PSD needs {n}")
    w = np.hanning(n)
    frames = np.lib.stride_tricks.sliding_window_view(x, n)[:: n // 2]
    P = np.mean(np.abs(np.fft.rfft(frames * w, axis=1)) ** 2, axis=0) / (w ** 2).sum()
    f = np.fft.rfftfreq(n, 1 / fs)
    out = []
    for fc in THIRD_OCT:
        sel = (f >= fc * 2 ** (-1 / 6)) & (f < fc * 2 ** (1 / 6))
        out.append(float(10 * np.log10(P[sel].sum() * 2 / n + 1e-30)) if sel.any() else float("nan"))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=("floor", "apply"))
    ap.add_argument("args", nargs="+")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--refdir", required=True)
    ap.add_argument("--lib", required=True)
    ap.add_argument("--description", default="chain.py floor")
    ap.add_argument("--seconds", type=float, default=10.0)
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    lib = json.loads(Path(a.lib).read_text())
    if a.cmd == "floor":
        ch = Chain.from_references(Path(a.manifest), a.args[0], Path(a.refdir), lib)
        fl, picks = ch.floor(int(a.seconds * ch.fs), a.description)
        if a.out:
            write_wav(Path(a.out), resample_pcm(quantize(fl, ch.bits)[0], ch.bits, ch.fs, SR))
        print(json.dumps(ch.receipt(a.description, fl, picks)))
        return
    src, dst, page = a.args
    ch = Chain.from_references(Path(a.manifest), page, Path(a.refdir), lib)
    y, rcpt = ch.apply(load(Path(src)), a.description)
    write_wav(Path(dst), y)
    print(json.dumps(rcpt))


if __name__ == "__main__":
    main()
