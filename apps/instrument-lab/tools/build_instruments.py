#!/usr/bin/env python3
"""Build the data behind the #instruments section: one card per instrument family.

    uv run --with numpy --with pillow python tools/build_instruments.py [--analysis DIR] [--pull]

Nothing on the page is typed by hand. Every number comes from three places:

  1. P2's analysis, one directory per graded note:   <analysis>/<family>/<note>/pair.json  (+ the files it names)
  2. the family lane's status file, prose only:       data/instruments/descriptions/<family>/status.json
  3. the reference library's provenance:              data/instruments/references.json

and the output (all under public/data/instruments/, which is gitignored like the rest of public/data) is

     index.json                      every family with its badge and counts (the overview)
     <family>/card.json              one family's card
     <family>/<note>/{real,ours}.m4a                       the A/B clips (AAC, like every other clip on the page);
                                     .mp4 with --audio-ext mp4 (same AAC, for a host that serves no .m4a)
     <family>/<note>/{spectrogram,scalogram}-{real,ours,diff}.webp

The default analysis directory is ~/.local/share/wa-inst/analysis (a mirror of trj:/data/tmp/wa-inst/analysis);
--pull refreshes it with rsync first.

----------------------------------------------------------------------------------------------------------------
pair.json (frankensim.instrument-pair.v1), written by P2's decompose.py. Paths are relative to its own directory.
A reference-only pair (no render yet) has "ours": null and no ours/diff files.

  schema          "frankensim.instrument-pair.v1"
  family, note    must equal the two directory names (family id as in FAMILIES below; note like "D3")
  dynamic         "pp" | "mf" | "ff"
  duration_s      length of BOTH clips in audio, from the aligned onset
  real            {reference_id: an id in data/instruments/references.json files[],
                   file:        the library file name,
                   segment_s:   [t0, t1] inside that file for this note, or null,
                   source:      only when the recording is NOT in references.json (e.g. VocalSet):
                                {id, name, credit, licence, url, licence_url?}, every string non-empty.
                                credit is the attribution text shown beside the clip (CC BY requires it).}
                  A library recording needs no source: its credit and licence come from references.json
                  and CREDITS below.
  ours            null, or {wav: the render's path on trj, git_sha, observer: "baffled_outlet 0 0 0.5" | ...,
                            sample_rate, recipe: repo path of the render.sh that made it}
  audio           {real: "real.wav", ours: "ours.wav" | null, level_match: one sentence saying how loudness was matched}
                  Aligned at the onset, trimmed to duration_s, and already loudness-matched: the page applies no
                  gain. The generator refuses a clip whose peak reaches full scale, or whose length is not
                  duration_s.
  alignment       {onset_real_s, onset_ours_s | null, highpass_hz, level_offset_db | null}   (shown, not used)
  figures         {spectrogram: F, scalogram: F} with F = {real, ours | null, diff | null,   (png | webp | jpg)
                                                          f_lo_hz, f_hi_hz, span_db, diff_range_db}
                  Bare maps: no axes, titles or colour bar. x = time 0..duration_s linear, y = frequency
                  f_lo_hz..f_hi_hz LOG with the low end at the bottom, the same geometry for all three images.
                  diff = ours minus real in dB (positive: ours louder), +-diff_range_db.
  partials        null, or {t_s: [frame times], f0_hz: {real, ours | null},
                            real: [{n, f_hz: [...], level_db: [...]}, ...],      one entry per tracked partial
                            ours: same shape | null}                              (null where untrackable;
                  ours null with a render present means none of its partials matched the recording's)
                  level_db is dB re the same reference for both sides (the clips are loudness-matched).
  metrics         [{id, label, unit, real, ours | null, gap | null, gap_unit (default: unit),
                    ruler: [{label, value}, ...], note}]
                  real/ours: a number (or null) shown side by side. gap: the distance between ours and real in
                  that metric. ruler: the same distance between REAL recordings (other dynamics, neighbouring
                  notes, other takes), value in gap_unit. An empty ruler means "no real-vs-real ruler for this
                  measure" and the page then makes no inside/outside claim for it.
                  A second, STRICT ruler: ruler_strict: [{label, value}, ...] holds only real recordings of the SAME note
                  (other dynamics, other takes); [] when there is none. inside_loose / inside_strict: true | false | null
                  (gap <= the largest ruler value; null = no gap or no ruler, no claim). verdicts (per pair with a render):
                  {measures, loose: {inside, outside, no_ruler}, strict: {...}}. The generator recomputes every one of
                  these from the gap and the rulers and refuses a pair whose own values disagree. The page shows both
                  counts; the badge is the weaker of the two. A pair without ruler_strict shows the loose verdict only.

A note the grader cannot grade honestly (reference off its name, onset check failed) is not deleted: decompose.py
overwrites its pair.json with {"refused": {"reason", "render_sha", "time"}, "ours": null} (family and note optional,
checked when present) and moves the old figures and audio into <note>/superseded/<UTC>/. The card lists the note as
"refused by the grader: <reason>", no badge counts it, and the build still exits 0 (the data is saying what is true).
A family whose only renders were refused reads "Rendered, not graded".

A render of an instrument with NO reference recording (pan pipes, tin whistle) is a listening-only note: "real": null and
"no_reference": {"reason": one sentence on why there is no recording}. It has ours and audio.ours (plus audio.level_match),
figures {spectrogram: F, scalogram: F} with F = {real: null, ours, diff: null, f_lo_hz, f_hi_hz, span_db}, and
alignment null, partials null, metrics []. Anything else on such a pair is refused. The card shows our audio and our
spectrogram and scalogram, no metrics and no recording credit, and a family with only such notes reads "Rendered, no
reference to grade against". Next to graded notes they do not change the family's badge.

--only FAMILY=NOTE,NOTE shows just those notes of a family (the review build's file limit). The others get no files, are
listed on the card as card.withheld, and are in no badge: the badge covers exactly the notes the page shows.
--recordings-only-if-rendered does the same for a family with no render at all: its recording-only notes (a real note
and its pictures, no FrankenSim) get no files and are listed as withheld, so a file-limited build spends its files on
what has been rendered.

status.json (frankensim.instrument-status.v1), written by the family lane, prose only (no measured numbers):
  family, title, reference_instrument (a references.json files[].instrument, or null), summary,
  physics_present: [str], physics_missing: [str], blocked: null | str, beads: [str]

The badge is computed here from the metrics: a measure is "within" when its gap is no larger than the largest
real-vs-real distance in its ruler, "outside" otherwise.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
APP = HERE.parent
REPO = APP.parents[1]
OUT = APP / "public" / "data" / "instruments"
REFERENCES = REPO / "data" / "instruments" / "references.json"
DESCRIPTIONS = REPO / "data" / "instruments" / "descriptions"
DEFAULT_ANALYSIS = Path.home() / ".local" / "share" / "wa-inst" / "analysis"
TRJ_ANALYSIS = "trj:/data/tmp/wa-inst/analysis/"

# Attribution line shown beside every clip, per references.json licence id (owner ruling 2026-10-20: the UIowa MIS
# recordings are free to use without restriction; credit the Electronic Music Studios).
CREDITS = {"uiowa-mis": "University of Iowa Electronic Music Studios"}
LIBRARY_NAME = {"uiowa-mis": "University of Iowa Musical Instrument Samples"}

PAIR_SCHEMA = "frankensim.instrument-pair.v1"
STATUS_SCHEMA = "frankensim.instrument-status.v1"

# id, title, references.json instrument (None: the library has no recording of this instrument)
FAMILIES: list[tuple[str, str, str | None]] = [
    ("clarinet", "Clarinet", "Bb clarinet"),
    ("saxophone", "Alto saxophone", "Alto saxophone"),
    ("oboe", "Oboe", "Oboe"),
    ("bassoon", "Bassoon", "Bassoon"),
    ("flute", "Flute", "Flute"),
    ("trumpet", "Trumpet", "Bb trumpet"),
    ("trombone", "Tenor trombone", "Tenor trombone"),
    ("horn", "French horn", "French horn"),
    ("tuba", "Tuba", "Tuba"),
    ("violin", "Violin", "Violin"),
    ("viola", "Viola", "Viola"),
    ("cello", "Cello", "Cello"),
    ("double-bass", "Double bass", "Double bass"),
    ("guitar", "Classical guitar", "Guitar"),
    ("marimba", "Marimba", "Marimba"),
    ("xylophone", "Xylophone", "Xylophone"),
    ("vibraphone", "Vibraphone", None),
    ("recorder", "Recorder", None),
    ("tin-whistle", "Tin whistle", None),
    ("pan-pipes", "Pan pipes", None),
    ("duduk", "Duduk", None),
]

MAX_PARTIALS = 8
MAX_FRAMES = 2000
BADGE_LABEL = {
    "none": "No render yet",
    "blocked": "Blocked",
    "ungraded": "Rendered, not graded",
    "listening": "Rendered, no reference to grade against",
    "close": "Inside the real spread",
    "mixed": "Partly inside the real spread",
    "off": "Outside the real spread",
}


class Refused(Exception):
    """A pair that cannot be shown honestly; carries every problem found."""

    def __init__(self, problems: list[str]):
        super().__init__("; ".join(problems))
        self.problems = problems


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def r(v, d=3):
    return None if v is None else round(float(v), d)


def sig(v, digits=6):
    """Round to significant digits, never to a fixed decimal place: an inharmonicity ruler of 3e-6 must not become 0."""
    return None if v is None else float(f"{float(v):.{digits}g}")


def finite(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def num_or_none(v):
    return float(v) if finite(v) else None


# ------------------------------------------------------------------ audio
def decode(path: Path) -> np.ndarray:
    """Mono float32 samples at 48 kHz, whatever the file's format."""
    res = subprocess.run(["ffmpeg", "-v", "error", "-i", str(path), "-ac", "1", "-ar", "48000", "-f", "f32le", "-"], capture_output=True)
    if res.returncode != 0:
        raise Refused([f"{path.name}: ffmpeg cannot decode it: {res.stderr.decode(errors='replace').strip()[:200]}"])
    return np.frombuffer(res.stdout, dtype="<f4")


def encode(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and dst.stat().st_mtime >= src.stat().st_mtime:
        return
    tmp = dst.with_suffix(f".tmp{dst.suffix}")  # ffmpeg picks the muxer from the extension, so it must match the final name
    for codec in (["-c:a", "aac_at", "-b:a", "256k"], ["-c:a", "aac", "-b:a", "256k"]):
        res = subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", str(src), "-ac", "1", *codec, str(tmp)], capture_output=True)
        if res.returncode == 0:
            tmp.replace(dst)
            return
    raise RuntimeError(f"ffmpeg could not encode {src}: {res.stderr.decode(errors='replace')[:300]}")


def convert_image(src: Path, dst: Path) -> tuple[int, int]:
    dst.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(src) as im:
        size = im.size
        if dst.exists() and dst.stat().st_mtime >= src.stat().st_mtime:
            return size
        im.convert("RGB").save(dst, "WEBP", quality=88, method=6)
    return size


# ------------------------------------------------------------------ pair validation
def check_figure_set(base: Path, name: str, fig, has_ours: bool, problems: list[str]) -> None:
    if not isinstance(fig, dict):
        problems.append(f"figures.{name} is missing")
        return
    for key in ("f_lo_hz", "f_hi_hz", "span_db", "diff_range_db"):
        if not finite(fig.get(key)) or fig[key] <= 0:
            problems.append(f"figures.{name}.{key} must be a positive number")
    if finite(fig.get("f_lo_hz")) and finite(fig.get("f_hi_hz")) and fig["f_hi_hz"] <= fig["f_lo_hz"]:
        problems.append(f"figures.{name}: f_hi_hz must exceed f_lo_hz")
    for side in ("real", "ours", "diff"):
        want = side == "real" or has_ours
        v = fig.get(side)
        if want and not (isinstance(v, str) and (base / v).is_file()):
            problems.append(f"figures.{name}.{side} names no file ({v!r})")
        if not want and v is not None:
            problems.append(f"figures.{name}.{side} must be null for a pair with no render")


def check_tracks(tracks, t_len: int, side: str, problems: list[str]) -> None:
    if not isinstance(tracks, list) or not tracks:
        problems.append(f"partials.{side} must be a non-empty list")
        return
    for tr in tracks:
        for key in ("f_hz", "level_db"):
            a = tr.get(key) if isinstance(tr, dict) else None
            if not isinstance(a, list) or len(a) != t_len:
                problems.append(f"partials.{side}[n={tr.get('n') if isinstance(tr, dict) else '?'}].{key} must have {t_len} entries (one per frame in t_s)")
                return


def is_listening(pj) -> bool:
    """A render with no reference recording: real is null and no_reference says why."""
    return isinstance(pj, dict) and pj.get("real") is None and pj.get("no_reference") is not None


def validate_listening(pj: dict, base: Path) -> list[str]:
    """The parts of a listening-only note: our render, our figures, and nothing that needs a recording."""
    problems: list[str] = []
    nr = pj["no_reference"]
    if not isinstance(nr, dict) or not isinstance(nr.get("reason"), str) or not nr["reason"].strip():
        problems.append("no_reference must be {reason} with a non-empty reason")
    ours = pj.get("ours")
    if not isinstance(ours, dict):
        problems.append("a note with no reference needs a render: ours is missing")
    else:
        for key in ("wav", "git_sha", "observer", "recipe"):
            if not isinstance(ours.get(key), str) or not ours[key]:
                problems.append(f"ours.{key} is required")
    audio = pj.get("audio")
    if not isinstance(audio, dict):
        problems.append("audio is missing")
    else:
        if audio.get("real") is not None:
            problems.append("audio.real must be null for a note with no reference")
        v = audio.get("ours")
        if not (isinstance(v, str) and (base / v).is_file()):
            problems.append(f"audio.ours names no file ({v!r})")
        if not isinstance(audio.get("level_match"), str) or not audio["level_match"]:
            problems.append("audio.level_match must say how loudness was set")
    figs = pj.get("figures")
    if not isinstance(figs, dict):
        problems.append("figures is missing")
    else:
        for name in ("spectrogram", "scalogram"):
            fig = figs.get(name)
            if not isinstance(fig, dict):
                problems.append(f"figures.{name} is missing")
                continue
            for key in ("f_lo_hz", "f_hi_hz", "span_db"):
                if not finite(fig.get(key)) or fig[key] <= 0:
                    problems.append(f"figures.{name}.{key} must be a positive number")
            if finite(fig.get("f_lo_hz")) and finite(fig.get("f_hi_hz")) and fig["f_hi_hz"] <= fig["f_lo_hz"]:
                problems.append(f"figures.{name}: f_hi_hz must exceed f_lo_hz")
            if not (isinstance(fig.get("ours"), str) and (base / fig["ours"]).is_file()):
                problems.append(f"figures.{name}.ours names no file ({fig.get('ours')!r})")
            for side in ("real", "diff"):
                if fig.get(side) is not None:
                    problems.append(f"figures.{name}.{side} must be null for a note with no reference")
    for key in ("alignment", "partials"):
        if pj.get(key) is not None:
            problems.append(f"{key} must be null for a note with no reference")
    if pj.get("metrics") != []:
        problems.append("metrics must be [] for a note with no reference: there is nothing to measure it against")
    return problems


def validate_pair(pj: dict, family: str, note: str, base: Path, refs_by_id: dict) -> list[str]:
    problems: list[str] = []
    if pj.get("schema") != PAIR_SCHEMA:
        problems.append(f"schema must be {PAIR_SCHEMA!r}, found {pj.get('schema')!r}")
    if pj.get("family") != family:
        problems.append(f"family {pj.get('family')!r} is not its directory {family!r}")
    if pj.get("note") != note:
        problems.append(f"note {pj.get('note')!r} is not its directory {note!r}")
    if pj.get("dynamic") not in ("pp", "mf", "ff"):
        problems.append(f"dynamic must be pp, mf or ff, found {pj.get('dynamic')!r}")
    if not finite(pj.get("duration_s")) or pj["duration_s"] <= 0:
        problems.append("duration_s must be a positive number")
    real = pj.get("real")
    if is_listening(pj):
        return problems + validate_listening(pj, base)
    if not isinstance(real, dict):
        problems.append("real is missing (a render with no reference recording sets real null and no_reference.reason)")
    elif real.get("reference_id") in refs_by_id:
        entry = refs_by_id[real["reference_id"]]
        lic = entry.get("licence") if isinstance(entry, dict) else None
        if lic not in CREDITS:
            problems.append(f"reference {real['reference_id']!r} has licence {lic!r}, which has no credit line in CREDITS (a library entry needs a licence id that CREDITS knows)")
    else:
        src = real.get("source")
        if not isinstance(src, dict):
            problems.append(f"real.reference_id {real.get('reference_id')!r} is not in data/instruments/references.json and real.source {{id, name, credit, licence, url}} is missing")
        else:
            for key in ("id", "name", "credit", "licence", "url"):
                if not isinstance(src.get(key), str) or not src[key].strip():
                    problems.append(f"real.source.{key} must be a non-empty string")
            if src.get("licence_url") is not None and not isinstance(src["licence_url"], str):
                problems.append("real.source.licence_url must be a string")
        if not isinstance(real.get("file"), str) or not real["file"]:
            problems.append("real.file is required for a recording that is not in references.json")
    has_ours = pj.get("ours") is not None
    if has_ours:
        ours = pj["ours"]
        for key in ("wav", "git_sha", "observer", "recipe"):
            if not isinstance(ours.get(key), str) or not ours[key]:
                problems.append(f"ours.{key} is required for a pair with a render")
    audio = pj.get("audio")
    if not isinstance(audio, dict):
        problems.append("audio is missing")
    else:
        for side in ("real", "ours"):
            want = side == "real" or has_ours
            v = audio.get(side)
            if want and not (isinstance(v, str) and (base / v).is_file()):
                problems.append(f"audio.{side} names no file ({v!r})")
            if not want and v is not None:
                problems.append("audio.ours must be null for a pair with no render")
        if not isinstance(audio.get("level_match"), str) or not audio["level_match"]:
            problems.append("audio.level_match must say how loudness was matched")
    figs = pj.get("figures")
    if not isinstance(figs, dict):
        problems.append("figures is missing")
    else:
        check_figure_set(base, "spectrogram", figs.get("spectrogram"), has_ours, problems)
        check_figure_set(base, "scalogram", figs.get("scalogram"), has_ours, problems)
    parts = pj.get("partials")
    if parts is not None:
        t = parts.get("t_s") if isinstance(parts, dict) else None
        if not isinstance(t, list) or not 2 <= len(t) <= MAX_FRAMES:
            problems.append(f"partials.t_s must hold 2..{MAX_FRAMES} frames")
        else:
            check_tracks(parts.get("real"), len(t), "real", problems)
            if has_ours:
                if parts.get("ours") is not None:  # null: the render exists but none of its partials matched the recording's
                    check_tracks(parts.get("ours"), len(t), "ours", problems)
            elif parts.get("ours") is not None:
                problems.append("partials.ours must be null for a pair with no render")
    metrics = pj.get("metrics")
    if not isinstance(metrics, list):
        problems.append("metrics must be a list")
    else:
        for i, m in enumerate(metrics):
            where = f"metrics[{i}]({m.get('id') if isinstance(m, dict) else '?'})"
            if not isinstance(m, dict) or not all(isinstance(m.get(k), str) and m[k] for k in ("id", "label")):
                problems.append(f"{where}: id and label are required")
                continue
            if not isinstance(m.get("ruler"), list):
                problems.append(f"{where}: ruler must be a list (empty when there is none)")
            else:
                for e in m["ruler"]:
                    if not (isinstance(e, dict) and isinstance(e.get("label"), str) and finite(e.get("value")) and e["value"] >= 0):
                        problems.append(f"{where}: every ruler entry needs a label and a non-negative value")
                        break
            if m.get("gap") is not None and (not finite(m["gap"]) or m["gap"] < 0):
                problems.append(f"{where}: gap must be a non-negative number or null")
            if m.get("gap") is not None and not has_ours:
                problems.append(f"{where}: a pair with no render has no gap")
            rs = m.get("ruler_strict")
            if "ruler_strict" in m and not (isinstance(rs, list) and all(isinstance(e, dict) and isinstance(e.get("label"), str) and finite(e.get("value")) and e["value"] >= 0 for e in rs)):
                problems.append(f"{where}: ruler_strict must be a list of entries with a label and a non-negative value (empty when the note has no other recording of the same note)")
            # decompose.py's own verdicts must be the ones the numbers give: two truths on one page are worse than a refused pair.
            for key, rows in (("inside_loose", m.get("ruler")), ("inside_strict", rs)):
                if key not in m or not isinstance(rows, list) or not all(isinstance(e, dict) and finite(e.get("value")) for e in rows) or not (m.get("gap") is None or finite(m["gap"])):
                    continue
                want = inside(m.get("gap"), [e["value"] for e in rows])
                if m[key] is not want:
                    problems.append(f"{where}: {key} is {m[key]!r} but gap {m.get('gap')!r} against its ruler gives {want!r}")
        verdicts = pj.get("verdicts")
        if verdicts is not None and not problems:
            if not isinstance(verdicts, dict):
                problems.append("verdicts must be an object or null")
            else:
                want = pair_verdicts([metric_view(m, has_ours) for m in metrics])
                for key in ("measures", "loose", "strict"):
                    if want[key] is not None and verdicts.get(key) != want[key]:
                        problems.append(f"verdicts.{key} is {verdicts.get(key)!r} but the measures give {want[key]!r}")
    return problems


def grader_refusal(pj, family: str, note: str) -> dict | None:
    """The grader's own refusal of a note, or None when pj is an ordinary pair.

    decompose.py refuses a note it cannot grade honestly by overwriting that note's pair.json with
    {"refused": {"reason", "render_sha", "time"}} (ours null) and moving the old figures and audio into superseded/.
    The note is listed on the card with its reason and counts toward no badge.
    """
    if not isinstance(pj, dict) or pj.get("refused") is None:
        return None
    rec, problems = pj["refused"], []
    if not isinstance(rec, dict) or not isinstance(rec.get("reason"), str) or not rec["reason"].strip():
        problems.append("refused must be {reason, render_sha, time} with a non-empty reason")
    for key in ("render_sha", "time"):
        if isinstance(rec, dict) and rec.get(key) is not None and not isinstance(rec[key], str):
            problems.append(f"refused.{key} must be a string")
    for key, want in (("family", family), ("note", note)):
        if pj.get(key) is not None and pj[key] != want:
            problems.append(f"{key} {pj[key]!r} is not its directory {want!r}")
    if pj.get("ours") is not None:
        problems.append("a refused note must have ours null")
    if problems:
        raise Refused(problems)
    return {"note": note, "kind": "grader", "reason": rec["reason"].strip(), "render_sha": rec.get("render_sha"), "time": rec.get("time"), "problems": [rec["reason"].strip()]}


def holds_a_render(pj_path: Path) -> bool:
    """True when this pair.json carries a render (a graded pair or a listening-only note), not a grader refusal."""
    try:
        pj = json.loads(pj_path.read_text())
    except json.JSONDecodeError:
        return False
    return isinstance(pj, dict) and pj.get("ours") is not None and pj.get("refused") is None


def check_audio(base: Path, pj: dict) -> list[str]:
    problems = []
    for side in ("real", "ours"):
        if pj["audio"].get(side) is None:
            continue
        x = decode(base / pj["audio"][side])
        dur = len(x) / 48000
        if abs(dur - pj["duration_s"]) > 0.05:
            problems.append(f"audio.{side} lasts {dur:.3f} s but duration_s is {pj['duration_s']}")
        peak = float(np.max(np.abs(x))) if len(x) else 0.0
        if peak >= 0.999:
            problems.append(f"audio.{side} peaks at full scale ({peak:.4f}): the A/B would clip; match loudness with headroom")
        if peak == 0.0:
            problems.append(f"audio.{side} is silent")
    return problems


# ------------------------------------------------------------------ card pieces
def inside(gap, values) -> bool | None:
    """True/False: the gap is no larger than the largest real-vs-real distance. None: no gap, or no ruler to judge it by (no claim)."""
    if gap is None or not values:
        return None
    return bool(gap <= max(values))


def metric_view(m: dict, has_ours: bool) -> dict:
    """One measure. ruler is the loose real-vs-real ruler (this and neighbouring notes, other dynamics and takes); ruler_strict
    the strict one (the same note only), None when the pair predates it. Both verdicts are recomputed here from the raw gap and
    ruler values, the way decompose.py does, so the page never shows a verdict the numbers do not give."""
    gap = num_or_none(m.get("gap")) if has_ours else None
    loose = [e["value"] for e in m["ruler"]]
    strict = None if m.get("ruler_strict") is None else [e["value"] for e in m["ruler_strict"]]
    top = max(loose, default=None)
    s_top = max(strict, default=None) if strict else None
    verdict = "no-gap" if gap is None else "no-ruler" if top is None else "within" if gap <= top else "outside"
    verdict_strict = None if strict is None else "no-gap" if gap is None else "no-ruler" if s_top is None else "within" if gap <= s_top else "outside"
    shown = lambda v: [num_or_none(x) for x in v] if isinstance(v, list) else num_or_none(v)  # noqa: E731
    rows = lambda rs: [{"label": e["label"], "value": sig(e["value"])} for e in rs]  # noqa: E731
    return {
        "id": m["id"], "label": m["label"], "unit": m.get("unit", ""), "gap_unit": m.get("gap_unit", m.get("unit", "")),
        "real": shown(m.get("real")), "ours": shown(m.get("ours")) if has_ours else None,
        "gap": sig(gap), "ruler": rows(m["ruler"]), "ruler_max": sig(top), "ruler_min": sig(min(loose, default=None)),
        "verdict": verdict, "ratio": sig(None if gap is None or not top else gap / top, 4), "inside_loose": inside(gap, loose),
        "ruler_strict": None if strict is None else rows(m["ruler_strict"]), "ruler_strict_max": sig(s_top), "ruler_strict_min": sig(min(strict, default=None)) if strict else None,
        "verdict_strict": verdict_strict, "ratio_strict": sig(None if gap is None or not s_top else gap / s_top, 4), "inside_strict": None if strict is None else inside(gap, strict),
        "note": m.get("note", ""),
    }


def pair_verdicts(views: list[dict]) -> dict:
    """The per-pair counts of both verdicts, in pair.json's own shape: the loose ruler and the strict one (None when a measure has
    no strict ruler field). no_ruler: a measure that has a gap and nothing to judge it by."""
    gapped = [v for v in views if v["gap"] is not None]
    tally = lambda key: {"inside": sum(x[key] is True for x in gapped), "outside": sum(x[key] is False for x in gapped), "no_ruler": sum(x[key] is None for x in gapped)}  # noqa: E731
    return {"measures": len(gapped), "loose": tally("inside_loose"), "strict": tally("inside_strict") if all(v["ruler_strict"] is not None for v in views) else None}


def partials_view(parts: dict | None) -> dict | None:
    """First MAX_PARTIALS tracked partials: level, and frequency as cents from the real recording's median."""
    if parts is None:
        return None
    by_n = lambda tracks: {int(t["n"]): t for t in (tracks or [])}  # noqa: E731
    real, ours = by_n(parts["real"]), by_n(parts.get("ours"))
    tracks = []
    for n in sorted(real)[:MAX_PARTIALS]:
        fr = [num_or_none(v) for v in real[n]["f_hz"]]
        ok = [v for v in fr if v is not None and v > 0]
        ref = statistics.median(ok) if ok else None

        def side(tr):
            if tr is None:
                return None
            f = [num_or_none(v) for v in tr["f_hz"]]
            cents = [None if ref is None or v is None or v <= 0 else r(1200 * math.log2(v / ref), 2) for v in f]
            return {"level_db": [r(num_or_none(v), 2) for v in tr["level_db"]], "cents": cents}

        tracks.append({"n": n, "ref_hz": r(ref, 3), "real": side(real[n]), "ours": side(ours.get(n))})
    return {"t_s": [r(v, 4) for v in parts["t_s"]], "f0_hz": {k: r(num_or_none(v), 4) for k, v in (parts.get("f0_hz") or {}).items()}, "tracks": tracks}


def real_source(real: dict, ref: dict | None, licences: dict) -> dict:
    """Who made the recording and on what terms, for the credit line: from references.json, else from the pair's own source."""
    if ref is not None:
        lid = ref.get("licence")
        lic = licences.get(lid, {})
        return {"id": ref.get("id") or real.get("reference_id"), "name": LIBRARY_NAME.get(lid, str(lid)), "credit": CREDITS.get(lid, ""), "licence": lic.get("text", str(lid)),
                "licence_url": lic.get("url"), "url": ref.get("url"), "external": False}
    s = real["source"]
    return {"id": s["id"].strip(), "name": s["name"].strip(), "credit": s["credit"].strip(), "licence": s["licence"].strip(), "licence_url": s.get("licence_url"), "url": s["url"].strip(), "external": True}


def build_pair(pj: dict, base: Path, out_dir: Path, refs_by_id: dict, chains: dict, licences: dict, audio_ext: str = "m4a") -> dict:
    has_ours = pj.get("ours") is not None
    note = pj["note"]
    rel = f"{out_dir.parent.name}/{out_dir.name}"
    listening = is_listening(pj)
    ref = None if listening else refs_by_id.get(pj["real"]["reference_id"])
    source = None if listening else real_source(pj["real"], ref, licences)
    audio = {}
    for side in ("real", "ours"):
        if pj["audio"].get(side):
            encode(base / pj["audio"][side], out_dir / f"{side}.{audio_ext}")
            audio[side] = f"{rel}/{side}.{audio_ext}"
        else:
            audio[side] = None
    figures = {}
    for kind in ("spectrogram", "scalogram"):
        fig = pj["figures"][kind]
        entry = {k: fig.get(k) for k in ("f_lo_hz", "f_hi_hz", "span_db", "diff_range_db")}  # a note with no reference has no diff_range_db
        for side in ("real", "ours", "diff"):
            if fig.get(side):
                w, h = convert_image(base / fig[side], out_dir / f"{kind}-{side}.webp")
                entry[side] = f"{rel}/{kind}-{side}.webp"
                entry["size"] = [w, h]
            else:
                entry[side] = None
        figures[kind] = entry
    chain = chains.get(ref.get("chain"), {}) if ref else {}
    metrics = [metric_view(m, has_ours) for m in pj["metrics"]]
    return {
        "note": note, "dynamic": pj["dynamic"], "duration_s": r(pj["duration_s"], 3),
        "real": None if listening else {
            "reference_id": source["id"], "file": pj["real"].get("file") or (Path(ref.get("path") or "").name if ref else ""),
            "segment_s": pj["real"].get("segment_s"), "url": source["url"], "sha256": ref.get("sha256") if ref else None,
            "licence": ref.get("licence") if ref else source["licence"], "source": source,
            "mic": chain.get("mic"), "pattern": chain.get("pattern"), "distance_m": chain.get("distance_m"), "distance_status": chain.get("distance_status"),
            "room": chain.get("room"), "format": chain.get("format"),
        },
        "no_reference": pj["no_reference"]["reason"].strip() if listening else None,
        "ours": pj.get("ours"),
        "audio": audio, "level_match": pj["audio"]["level_match"],
        "alignment": pj.get("alignment"),
        "figures": figures,
        "partials": partials_view(pj.get("partials")),
        "metrics": metrics, "verdicts": pair_verdicts(metrics) if has_ours and not listening else None,
        "source_sha256": None,  # filled by the caller
    }


def judge(status: dict | None, pairs: list[dict], grader_refused: list[dict] | None = None) -> dict:
    """The badge, and the sentence under it. Both are computed from the metrics; a note the grader refused counts in none."""
    listening = [p for p in pairs if p["real"] is None]  # renders with no reference recording: nothing to grade them against
    pairs = [p for p in pairs if p["real"] is not None]
    rendered = [p for p in pairs if p["ours"]]
    blocked = (status or {}).get("blocked")
    refused_renders = [x for x in (grader_refused or []) if x.get("render_sha")]
    counts = {"notes_with_render": len(rendered) + len(listening), "notes_reference_only": len(pairs) - len(rendered), "notes_no_reference": len(listening),
              "notes_refused_by_grader": len(grader_refused or []), "within": 0, "outside": 0, "unruled": 0}
    if not rendered:
        if listening:
            n = len(listening)
            return {"badge": "listening", "label": BADGE_LABEL["listening"], "counts": counts,
                    "line": f"{n} note{'s' if n != 1 else ''} rendered. There is no reference recording for this instrument, so nothing is compared or graded: this is audio to listen to, with its spectrogram and scalogram."}
        if refused_renders:
            why = "; ".join(f"{x['note']}: {x['reason']}" for x in refused_renders)
            return {"badge": "ungraded", "label": BADGE_LABEL["ungraded"], "counts": counts,
                    "line": f"{len(refused_renders)} note{'s' if len(refused_renders) != 1 else ''} rendered, but the grader refused to compare {'them' if len(refused_renders) != 1 else 'it'} with the recording, so nothing is graded. {why}"}
        if blocked:
            return {"badge": "blocked", "label": BADGE_LABEL["blocked"], "line": blocked, "counts": counts}
        extra = f" The reference library has real recordings for {len(pairs)} note{'s' if len(pairs) != 1 else ''} below." if pairs else ""
        return {"badge": "none", "label": BADGE_LABEL["none"], "line": "No render has been made for this family yet." + extra, "counts": counts}
    counts.update(strict_within=0, strict_outside=0, strict_unruled=0, strict_notes=0)
    worst = worst_strict = None
    for p in rendered:
        if p["metrics"] and all(m["ruler_strict"] is not None for m in p["metrics"]):
            counts["strict_notes"] += 1  # this note was graded against the same-note ruler too
        for m in p["metrics"]:
            if m["verdict"] == "within":
                counts["within"] += 1
            elif m["verdict"] == "outside":
                counts["outside"] += 1
                if worst is None or (m["ratio"] or 0) > (worst[0]["ratio"] or 0):
                    worst = (m, p["note"])
            elif m["verdict"] == "no-ruler":
                counts["unruled"] += 1
            if m["verdict_strict"] == "within":
                counts["strict_within"] += 1
            elif m["verdict_strict"] == "outside":
                counts["strict_outside"] += 1
                if worst_strict is None or (m["ratio_strict"] or 0) > (worst_strict[0]["ratio_strict"] or 0):
                    worst_strict = (m, p["note"])
            elif m["verdict_strict"] == "no-ruler":
                counts["strict_unruled"] += 1
    ruled = counts["within"] + counts["outside"]
    n = len(rendered)
    notes = f"{n} note{'s' if n != 1 else ''}"
    if ruled == 0:
        return {"badge": "ungraded", "label": BADGE_LABEL["ungraded"], "line": f"{notes} rendered; none has a measure with a real-vs-real ruler yet.", "counts": counts}
    s_ruled = counts["strict_within"] + counts["strict_outside"]
    level = lambda inside_n, outside_n: "close" if outside_n == 0 else "mixed" if inside_n / (inside_n + outside_n) >= 0.5 else "off"  # noqa: E731
    levels = [level(counts["within"], counts["outside"])] + ([level(counts["strict_within"], counts["strict_outside"])] if s_ruled else [])
    badge = min(levels, key=("off", "mixed", "close").index)  # the weaker of the two verdicts
    line = f"Over {notes}: inside the spread of real recordings of this and neighbouring notes, {counts['within']} of {ruled} measures. "
    if s_ruled:
        none = f" ({counts['strict_unruled']} {'has' if counts['strict_unruled'] == 1 else 'have'} no same-note comparison)" if counts["strict_unruled"] else ""
        part = f"; computed for {counts['strict_notes']} of {n} notes" if counts["strict_notes"] < n else ""
        line += f"Inside the spread of real recordings of the same note, {counts['strict_within']} of {s_ruled} measures{none}{part}. The badge is the weaker of the two."
    elif counts["strict_notes"]:
        line += "No same-note comparison could be made: there is no other recording of the same note to judge against."
    else:
        line += "The same-note comparison has not been computed for these notes."
    if s_ruled and worst_strict:
        m, note = worst_strict
        line += f" Largest gap against the same-note spread: {m['label']}, {m['ratio_strict'] or 0:.1f}× it ({note})."
    elif worst:
        m, note = worst
        line += f" Largest gap against the neighbouring-note spread: {m['label']}, {m['ratio'] or 0:.1f}× it ({note})."
    return {"badge": badge, "label": BADGE_LABEL[badge], "line": line, "counts": counts}


def load_status(family: str, root: Path = DESCRIPTIONS) -> dict | None:
    path = root / family / "status.json"
    if not path.is_file():
        return None
    s = json.loads(path.read_text())
    problems = []
    if s.get("schema") != STATUS_SCHEMA:
        problems.append(f"schema must be {STATUS_SCHEMA!r}")
    for key in ("physics_present", "physics_missing", "beads"):
        if not isinstance(s.get(key, []), list) or not all(isinstance(x, str) for x in s.get(key, [])):
            problems.append(f"{key} must be a list of strings")
    if s.get("blocked") is not None and not isinstance(s["blocked"], str):
        problems.append("blocked must be null or a string")
    if problems:
        raise SystemExit(f"{path}: {'; '.join(problems)}")
    return s


def library_summary(refs: dict, instrument: str) -> dict:
    files = [f for f in refs["files"] if f["instrument"] == instrument and f["kind"] != "ambient-silence"]
    by_dyn: dict[str, int] = {}
    for f in files:
        by_dyn[f["dynamic"] or "unlabelled"] = by_dyn.get(f["dynamic"] or "unlabelled", 0) + 1
    chains = sorted({f["chain"] for f in files})
    licence_id = files[0]["licence"] if files else None
    lic = {l["id"]: l for l in refs["licences"]}.get(licence_id, {})
    return {"instrument": instrument, "files": len(files), "by_dynamic": by_dyn, "chains": chains,
            "licence": lic.get("text", licence_id), "licence_url": lic.get("url"), "credit": CREDITS.get(licence_id)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--analysis", type=Path, default=DEFAULT_ANALYSIS, help=f"P2's analysis tree (default {DEFAULT_ANALYSIS})")
    ap.add_argument("--pull", action="store_true", help=f"rsync {TRJ_ANALYSIS} into --analysis first")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--descriptions", type=Path, default=DESCRIPTIONS, help="where <family>/status.json live")
    ap.add_argument("--audio-ext", choices=("m4a", "mp4"), default="m4a", help="extension of the A/B clips; the codec is AAC either way (mp4 for a host that serves no .m4a)")
    ap.add_argument("--only", action="append", default=[], metavar="FAMILY=NOTE,NOTE", help="show only these notes of a family; the rest are listed on its card as withheld and get no files")
    ap.add_argument("--recordings-only-if-rendered", action="store_true", help="a family with no render shows none of its recording-only notes; they are listed as withheld and get no files")
    args = ap.parse_args()
    only: dict[str, set[str]] = {}
    for spec in args.only:
        fam, sep, notes = spec.partition("=")
        if not sep or not fam or not notes:
            ap.error(f"--only wants FAMILY=NOTE,NOTE (got {spec!r})")
        only.setdefault(fam, set()).update(notes.split(","))

    if args.pull:
        args.analysis.mkdir(parents=True, exist_ok=True)
        subprocess.run(["rsync", "-a", "--delete", TRJ_ANALYSIS, f"{args.analysis}/"], check=True)

    refs = json.loads(REFERENCES.read_text())
    refs_by_id = {f["id"]: f for f in refs["files"]}
    chains = {c["id"]: c for c in refs["chains"]}
    licences = {l["id"]: l for l in refs["licences"]}

    known = {fid for fid, _, _ in FAMILIES}
    extra = sorted(p.parent.name for p in args.descriptions.glob("*/status.json") if p.parent.name not in known) if args.descriptions.is_dir() else []
    families = list(FAMILIES) + [(fid, (load_status(fid, args.descriptions) or {}).get("title", fid), (load_status(fid, args.descriptions) or {}).get("reference_instrument")) for fid in extra]

    refused_total = 0
    total_grader_refused = 0
    index = []
    args.out.mkdir(parents=True, exist_ok=True)
    for fid, title, ref_instrument in families:
        status = load_status(fid, args.descriptions)
        pairs, refused, withheld = [], [], []
        fam_dir = args.analysis / fid
        pair_files = sorted(fam_dir.glob("*/pair.json")) if fam_dir.is_dir() else []
        has_render = any(holds_a_render(p) for p in pair_files)
        for pj_path in pair_files:
            note = pj_path.parent.name
            base = pj_path.parent
            try:
                pj = json.loads(pj_path.read_text())
                refusal = grader_refusal(pj, fid, note)
                if refusal:
                    refused.append(refusal)
                    continue
                if (fid in only and note not in only[fid]) or (args.recordings_only_if_rendered and pj.get("ours") is None and not has_render):
                    withheld.append(note)
                    continue
                problems = validate_pair(pj, fid, note, base, refs_by_id)
                if not problems:
                    problems = check_audio(base, pj)
                if problems:
                    raise Refused(problems)
                pair = build_pair(pj, base, args.out / fid / note, refs_by_id, chains, licences, args.audio_ext)
                pair["source_sha256"] = sha256(pj_path)
                pairs.append(pair)
            except (Refused, json.JSONDecodeError) as e:
                refused.append({"note": note, "kind": "site", "problems": e.problems if isinstance(e, Refused) else [f"pair.json is not JSON: {e}"]})
                print(f"REFUSED {fid}/{note}: {'; '.join(refused[-1]['problems'])}", file=sys.stderr)
        seen = {p["note"] for p in pairs} | set(withheld) | {x["note"] for x in refused}
        for missing in sorted(only.get(fid, set()) - seen):
            print(f"WARNING --only {fid}={missing}: no such note in the analysis", file=sys.stderr)
        refused_total += sum(1 for x in refused if x["kind"] == "site")
        grader_refused =[x for x in refused if x["kind"] == "grader"]
        pairs.sort(key=lambda p: (p["ours"] is None, p["dynamic"], p["note"]))
        verdict = judge(status, pairs, grader_refused)
        card = {
            "schema": "frankensim.instrument-card.v1", "family": fid, "title": (status or {}).get("title", title),
            "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "reference": library_summary(refs, ref_instrument) if ref_instrument in {f["instrument"] for f in refs["files"]} else None,
            "status": {
                **verdict, "summary": (status or {}).get("summary", ""), "blocked": (status or {}).get("blocked"),
                "physics_present": (status or {}).get("physics_present", []), "physics_missing": (status or {}).get("physics_missing", []),
                "beads": (status or {}).get("beads", []), "described": status is not None,
            },
            "pairs": pairs, "refused": refused, "withheld": sorted(withheld),
        }
        (args.out / fid).mkdir(parents=True, exist_ok=True)
        (args.out / fid / "card.json").write_text(json.dumps(card, ensure_ascii=False, separators=(",", ":")))
        index.append({
            "id": fid, "title": card["title"], "reference": card["reference"]["instrument"] if card["reference"] else None, "badge": verdict["badge"], "label": verdict["label"],
            "line": verdict["line"], "counts": verdict["counts"], "notes": [{"note": p["note"], "dynamic": p["dynamic"], "render": p["ours"] is not None} for p in pairs],
            "refused": len(refused),
        })
        total_grader_refused += len(grader_refused)
    (args.out / "index.json").write_text(json.dumps({"schema": "frankensim.instrument-index.v1", "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "families": index}, ensure_ascii=False, separators=(",", ":")))
    graded = sum(1 for f in index if f["badge"] in ("close", "mixed", "off"))
    print(f"{len(index)} families, {sum(len(f['notes']) for f in index)} notes shown, {graded} graded, {total_grader_refused} notes refused by the grader, {refused_total} pairs refused by this generator -> {args.out}")
    return 2 if refused_total else 0


if __name__ == "__main__":
    raise SystemExit(main())
