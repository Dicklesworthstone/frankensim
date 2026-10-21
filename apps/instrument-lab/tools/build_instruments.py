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

A note the grader cannot grade honestly (reference off its name, onset check failed) is not deleted: decompose.py
overwrites its pair.json with {"refused": {"reason", "render_sha", "time"}, "ours": null} (family and note optional,
checked when present) and moves the old figures and audio into <note>/superseded/<UTC>/. The card lists the note as
"refused by the grader: <reason>", no badge counts it, and the build still exits 0 (the data is saying what is true).
A family whose only renders were refused reads "Rendered, not graded".

--only FAMILY=NOTE,NOTE shows just those notes of a family (the review build's file limit). The others get no files, are
listed on the card as card.withheld, and are in no badge: the badge covers exactly the notes the page shows.

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
    if not isinstance(real, dict):
        problems.append("real is missing")
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
def metric_view(m: dict, has_ours: bool) -> dict:
    gap = num_or_none(m.get("gap")) if has_ours else None
    ruler = [{"label": e["label"], "value": sig(e["value"])} for e in m["ruler"]]
    top = max((e["value"] for e in ruler), default=None)
    verdict = "no-gap" if gap is None else "no-ruler" if top is None else "within" if gap <= top else "outside"
    ratio = None if gap is None or not top else gap / top
    shown = lambda v: [num_or_none(x) for x in v] if isinstance(v, list) else num_or_none(v)  # noqa: E731
    return {
        "id": m["id"], "label": m["label"], "unit": m.get("unit", ""), "gap_unit": m.get("gap_unit", m.get("unit", "")),
        "real": shown(m.get("real")), "ours": shown(m.get("ours")) if has_ours else None,
        "gap": sig(gap), "ruler": ruler, "ruler_max": sig(top), "ruler_min": sig(min((e["value"] for e in ruler), default=None)),
        "verdict": verdict, "ratio": sig(ratio, 4), "note": m.get("note", ""),
    }


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
    ref = refs_by_id.get(pj["real"]["reference_id"])
    source = real_source(pj["real"], ref, licences)
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
        entry = {k: fig[k] for k in ("f_lo_hz", "f_hi_hz", "span_db", "diff_range_db")}
        for side in ("real", "ours", "diff"):
            if fig.get(side):
                w, h = convert_image(base / fig[side], out_dir / f"{kind}-{side}.webp")
                entry[side] = f"{rel}/{kind}-{side}.webp"
                entry["size"] = [w, h]
            else:
                entry[side] = None
        figures[kind] = entry
    chain = chains.get(ref.get("chain"), {}) if ref else {}
    return {
        "note": note, "dynamic": pj["dynamic"], "duration_s": r(pj["duration_s"], 3),
        "real": {
            "reference_id": source["id"], "file": pj["real"].get("file") or (Path(ref.get("path") or "").name if ref else ""),
            "segment_s": pj["real"].get("segment_s"), "url": source["url"], "sha256": ref.get("sha256") if ref else None,
            "licence": ref.get("licence") if ref else source["licence"], "source": source,
            "mic": chain.get("mic"), "pattern": chain.get("pattern"), "distance_m": chain.get("distance_m"), "distance_status": chain.get("distance_status"),
            "room": chain.get("room"), "format": chain.get("format"),
        },
        "ours": pj.get("ours"),
        "audio": audio, "level_match": pj["audio"]["level_match"],
        "alignment": pj.get("alignment"),
        "figures": figures,
        "partials": partials_view(pj.get("partials")),
        "metrics": [metric_view(m, has_ours) for m in pj["metrics"]],
        "source_sha256": None,  # filled by the caller
    }


def judge(status: dict | None, pairs: list[dict], grader_refused: list[dict] | None = None) -> dict:
    """The badge, and the sentence under it. Both are computed from the metrics; a note the grader refused counts in none."""
    rendered = [p for p in pairs if p["ours"]]
    blocked = (status or {}).get("blocked")
    refused_renders = [x for x in (grader_refused or []) if x.get("render_sha")]
    counts = {"notes_with_render": len(rendered), "notes_reference_only": len(pairs) - len(rendered), "notes_refused_by_grader": len(grader_refused or []), "within": 0, "outside": 0, "unruled": 0}
    if not rendered:
        if refused_renders:
            why = "; ".join(f"{x['note']}: {x['reason']}" for x in refused_renders)
            return {"badge": "ungraded", "label": BADGE_LABEL["ungraded"], "counts": counts,
                    "line": f"{len(refused_renders)} note{'s' if len(refused_renders) != 1 else ''} rendered, but the grader refused to compare {'them' if len(refused_renders) != 1 else 'it'} with the recording, so nothing is graded. {why}"}
        if blocked:
            return {"badge": "blocked", "label": BADGE_LABEL["blocked"], "line": blocked, "counts": counts}
        extra = f" The reference library has real recordings for {len(pairs)} note{'s' if len(pairs) != 1 else ''} below." if pairs else ""
        return {"badge": "none", "label": BADGE_LABEL["none"], "line": "No render has been made for this family yet." + extra, "counts": counts}
    worst = None
    for p in rendered:
        for m in p["metrics"]:
            if m["verdict"] == "within":
                counts["within"] += 1
            elif m["verdict"] == "outside":
                counts["outside"] += 1
                if worst is None or m["ratio"] > worst[0]["ratio"]:
                    worst = (m, p["note"])
            elif m["verdict"] == "no-ruler":
                counts["unruled"] += 1
    ruled = counts["within"] + counts["outside"]
    if ruled == 0:
        return {"badge": "ungraded", "label": BADGE_LABEL["ungraded"], "line": f"{len(rendered)} note{'s' if len(rendered) != 1 else ''} rendered; none has a measure with a real-vs-real ruler yet.", "counts": counts}
    frac = counts["within"] / ruled
    badge = "close" if counts["outside"] == 0 else "mixed" if frac >= 0.5 else "off"
    line = f"{counts['within']} of {ruled} graded measures, over {len(rendered)} note{'s' if len(rendered) != 1 else ''}, are no further from the real recording than real recordings are from each other."
    if worst:
        m, note = worst
        line += f" Largest gap against that spread: {m['label']}, {m['ratio']:.1f}× it ({note})."
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
        for pj_path in sorted(fam_dir.glob("*/pair.json")) if fam_dir.is_dir() else []:
            note = pj_path.parent.name
            base = pj_path.parent
            try:
                pj = json.loads(pj_path.read_text())
                refusal = grader_refusal(pj, fid, note)
                if refusal:
                    refused.append(refusal)
                    continue
                if fid in only and note not in only[fid]:
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
