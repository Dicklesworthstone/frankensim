#!/usr/bin/env python3
"""Build the data behind the "Reference candidates" section of the review page: real recordings found for every instrument family,
for the owner to mark VALID or INVALID. Nothing here is used for grading until he confirms it.

    uv run --with numpy python tools/build_candidates.py --candidates DIR [--out DIR]

Input: a mirror of trj:/data/tmp/wa-inst/refs-candidates (never in the repo), laid out as its README says:

    <DIR>/families.json                         per family the notes our renders are graded on (optional; shown on the family's page)
    <DIR>/<family>/candidates.json              the merged, validated file for the family (P2)
    <DIR>/<family>/candidates-library*.json     until it exists: open libraries we hold (UIowa, TinySOL, VocalSet, ...), read together with
    <DIR>/<family>/candidates-external*.json    ... everything else (other libraries, Freesound, Commons, YouTube), one file per hunter
    <DIR>/<family>/<id>/clip.wav, facts.json    the clip a record names and its measured facts

Each file is {schema: frankensim.reference-candidates.v1, family, made_by, made_utc, looked: [{source, query_or_url, result}],
none_found: null | {reason, looked_summary}, candidates: [...]}; a candidate has id, family, instrument, notes, dynamic, articulation,
source_type (open-library | youtube | freesound | other), url, title, channel_or_author, licence, licence_evidence, start_s, end_s,
conditions (what the source says and what was measured, and that it was not listened to), sample_rate_hz, reason, files {clip, facts},
sha256_clip, and optionally usable_for_render_notes (a list, a note -> category map or a sentence: read by coverage()). A record that fails a check is refused and listed, never shown as if it were fine.

Output (public/data/candidates/): index.json, and <family>/<id>.mp4 (AAC; every clip at the same loudness so no candidate is louder
than another, the gain applied is in its record and the measured levels of the original are kept as measured).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
APP = HERE.parent
OUT = APP / "public" / "data" / "candidates"
SCHEMA = "frankensim.reference-candidates.v1"
INDEX_SCHEMA = "frankensim.reference-candidates-index.v1"
SOURCE_TYPES = ("open-library", "youtube", "freesound", "other")
TARGET_RMS_DBFS = -26.0   # every clip is played at this RMS, so the owner is not judging loudness
PEAK_LIMIT_DBFS = -1.0    # unless that would put its peak above this
DURATION_SLACK_S = 0.3    # a record's end_s - start_s against the clip it names (an mp3 or a cut on a frame edge is not exact)
SAFE_ID = re.compile(r"^[A-Za-z0-9._-]+$")

# The families the owner named, in his order: id, title.
FAMILIES = [
    ("violin", "Violin"), ("viola", "Viola"), ("cello", "Cello"), ("double-bass", "Double bass"), ("guitar", "Classical guitar"),
    ("trumpet", "Trumpet"), ("trombone", "Trombone"), ("horn", "French horn"), ("tuba", "Tuba"),
    ("clarinet", "Clarinet"), ("oboe", "Oboe"), ("bassoon", "Bassoon"), ("saxophone", "Saxophone (alto, tenor)"),
    ("flute", "Flute"), ("recorder", "Alto recorder"), ("tin-whistle", "Tin whistle"), ("pan-pipes", "Pan pipes"), ("duduk", "Duduk"),
    ("xylophone", "Xylophone"), ("marimba", "Marimba"), ("vibraphone", "Vibraphone"), ("voice", "Voice (sustained vowels)"),
]


class Refused(Exception):
    def __init__(self, problems: list[str]):
        super().__init__("; ".join(problems))
        self.problems = problems


def text(v) -> bool:
    return isinstance(v, str) and bool(v.strip())


def number(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def decode(path: Path) -> tuple[np.ndarray, int]:
    """Mono float samples at the file's own rate."""
    probe = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "a:0", "-show_entries", "stream=sample_rate", "-of", "json", str(path)], capture_output=True, text=True)
    if probe.returncode != 0 or not json.loads(probe.stdout or "{}").get("streams"):
        raise Refused([f"{path.name}: ffprobe cannot read it: {probe.stderr.strip()[:200]}"])
    sr = int(json.loads(probe.stdout)["streams"][0]["sample_rate"])
    res = subprocess.run(["ffmpeg", "-v", "error", "-i", str(path), "-ac", "1", "-f", "f32le", "-"], capture_output=True)
    if res.returncode != 0:
        raise Refused([f"{path.name}: ffmpeg cannot decode it: {res.stderr.decode(errors='replace').strip()[:200]}"])
    return np.frombuffer(res.stdout, dtype="<f4").astype(np.float64), sr


def playback_gain_db(x: np.ndarray) -> float:
    """The gain that brings the clip to TARGET_RMS_DBFS, held back so its peak stays under PEAK_LIMIT_DBFS."""
    rms = math.sqrt(float(np.mean(x ** 2))) if len(x) else 0.0
    peak = float(np.max(np.abs(x))) if len(x) else 0.0
    if rms <= 0 or peak <= 0:
        raise Refused(["the clip is silent"])
    return round(min(TARGET_RMS_DBFS - 20 * math.log10(rms), PEAK_LIMIT_DBFS - 20 * math.log10(peak)), 2)


def encode(src: Path, dst: Path, gain_db: float) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(f".tmp{dst.suffix}")  # ffmpeg picks the muxer from the extension, so it must match the final name
    res = None
    for codec in (["-c:a", "aac_at", "-b:a", "96k"], ["-c:a", "aac", "-b:a", "96k"]):
        res = subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", str(src), "-ac", "1", "-af", f"volume={gain_db}dB", *codec, str(tmp)], capture_output=True)
        if res.returncode == 0:
            tmp.replace(dst)
            return
    raise RuntimeError(f"ffmpeg could not encode {src}: {res.stderr.decode(errors='replace')[:300]}")


def is_cc(licence: str) -> bool:
    """Not a search: YouTube's standard-licence wording can say 'not Creative Commons'."""
    return licence.strip().lower().startswith(("creative commons", "cc ", "cc-", "cc0"))


NOTE_NAME = re.compile(r"\b[A-G](?:b|#)?\d\b")


def coverage(v) -> dict | None:
    """Which of the notes our renders are graded on a clip covers, whatever shape the hunter wrote (the field is optional and loosely
    described): a list of notes, a note -> category map (families.json's categories), or a sentence such as "graded notes covered: A3, D4".
    Returns {graded, reference_only, other, text}, or None for a shape that is none of these."""
    out = {"graded": [], "reference_only": [], "other": [], "text": None}
    if v is None:
        return out
    if isinstance(v, list) and all(isinstance(n, str) for n in v):
        out["graded"] = [n.strip() for n in v if n.strip()]
    elif isinstance(v, dict) and all(isinstance(k, str) and isinstance(x, str) for k, x in v.items()):
        for note, kind in v.items():
            out["graded" if kind == "rendered_and_graded" else "reference_only" if kind == "reference_only" else "other"].append(note)
    elif isinstance(v, str):
        m = re.match(r"\s*graded notes covered:\s*([^;(]*)", v, re.I)
        out["graded"] = NOTE_NAME.findall(m.group(1)) if m else []
        out["text"] = v.strip()
    else:
        return None
    return out


def validate(rec, family: str, base: Path) -> tuple[list[str], dict | None]:
    """Every check a record must pass before the owner is asked about it. Returns the problems and the clip's facts.json (or None)."""
    if not isinstance(rec, dict):
        return ["the candidate is not an object"], None
    p: list[str] = []
    if not (text(rec.get("id")) and SAFE_ID.match(rec["id"])):
        p.append("id must be a filesystem-safe string (letters, digits, . _ -)")
    if rec.get("family") != family:
        p.append(f"family {rec.get('family')!r} is not {family!r} (its folder)")
    for key in ("instrument", "reason", "dynamic", "articulation", "title", "channel_or_author"):
        if not text(rec.get(key)):
            p.append(f"{key} must be a non-empty string (say 'not stated' when the source does not)")
    if not isinstance(rec.get("notes"), list) or not rec["notes"] or not all(text(n) for n in rec["notes"]):
        p.append("notes must be a non-empty list of note names (or ['various'] for a scale)")
    if rec.get("source_type") not in SOURCE_TYPES:
        p.append(f"source_type must be one of {', '.join(SOURCE_TYPES)}")
    if not (text(rec.get("url")) and rec["url"].startswith(("http://", "https://"))):
        p.append("url must be an http(s) address")
    for key in ("licence", "licence_evidence"):
        if not text(rec.get(key)):
            p.append(f"{key} must be stated, with where it was read (a video with no licence field is 'YouTube standard licence (all rights reserved)')")
    if not (number(rec.get("start_s")) and number(rec.get("end_s")) and 0 <= rec["start_s"] < rec["end_s"]):
        p.append("start_s < end_s, in seconds of the source")
    if not (text(rec.get("conditions")) and "not listened to" in rec["conditions"].lower()):
        p.append("conditions must say what the source states, what was measured, and that it was NOT LISTENED TO (nobody who found a candidate could hear it)")
    if not (isinstance(rec.get("sample_rate_hz"), int) and rec["sample_rate_hz"] > 0):
        p.append("sample_rate_hz must be a whole number of hertz")
    if coverage(rec.get("usable_for_render_notes")) is None:
        p.append("usable_for_render_notes must be a list of note names, a note -> category map, or a sentence")
    files = rec.get("files")
    facts = None
    if not (isinstance(files, dict) and text(files.get("clip"))):
        return p + ["files.clip must name the clip"], None
    clip = base / files["clip"]
    if not clip.is_file():
        return p + [f"clip file {files['clip']!r} is missing"], None
    if hashlib.sha256(clip.read_bytes()).hexdigest() != rec.get("sha256_clip"):
        p.append("sha256_clip does not match the clip file")
    try:
        x, sr = decode(clip)
    except Refused as e:
        return p + e.problems, None
    if sr != rec.get("sample_rate_hz"):
        p.append(f"sample_rate_hz says {rec.get('sample_rate_hz')} but the clip is {sr} Hz")
    if number(rec.get("start_s")) and number(rec.get("end_s")) and abs(len(x) / sr - (rec["end_s"] - rec["start_s"])) > DURATION_SLACK_S:
        p.append(f"the clip lasts {len(x) / sr:.2f} s but start_s..end_s is {rec['end_s'] - rec['start_s']:.2f} s")
    if text(files.get("facts")):
        fp = base / files["facts"]
        try:
            facts = json.loads(fp.read_text())
            if not isinstance(facts, dict):
                raise ValueError("not an object")
        except (OSError, ValueError) as e:
            p.append(f"files.facts {files['facts']!r} cannot be read: {e}")
            facts = None
    return p, facts


def read_family_files(fam_dir: Path, fid: str) -> tuple[list[dict], list[dict], list[dict], list[str]]:
    """The records of a family: (candidates, looked, none_found notes, made_by), from candidates.json or else every library and external file."""
    paths = [fam_dir / "candidates.json"] if (fam_dir / "candidates.json").is_file() else sorted([*fam_dir.glob("candidates-library*.json"), *fam_dir.glob("candidates-external*.json")])
    cands, looked, nones, made = [], [], [], []
    for path in paths:
        try:
            doc = json.loads(path.read_text())
        except json.JSONDecodeError as e:
            raise Refused([f"{path.name} is not JSON: {e}"])
        p = []
        if not isinstance(doc, dict) or doc.get("schema") != SCHEMA:
            p.append(f"{path.name}: schema must be {SCHEMA!r}")
        elif doc.get("family") != fid:
            p.append(f"{path.name}: family {doc.get('family')!r} is not {fid!r}")
        elif not (isinstance(doc.get("looked"), list) and doc["looked"] and all(isinstance(x, dict) and all(text(x.get(k)) for k in ("source", "query_or_url", "result")) for x in doc["looked"])):
            p.append(f"{path.name}: looked must list every source and query tried, each with source, query_or_url and result")
        elif not isinstance(doc.get("candidates"), list):
            p.append(f"{path.name}: candidates must be a list")
        elif not doc["candidates"] and not (isinstance(doc.get("none_found"), dict) and text(doc["none_found"].get("reason"))):
            p.append(f"{path.name}: no candidates and no none_found.reason")
        if p:
            raise Refused(p)
        cands += [(c, path.name) for c in doc["candidates"]]
        looked += [{**x, "file": path.name} for x in doc["looked"]]
        if isinstance(doc.get("none_found"), dict):
            nones.append({"reason": doc["none_found"]["reason"].strip(), "looked_summary": str(doc["none_found"].get("looked_summary", "")).strip(), "file": path.name})
        made.append(str(doc.get("made_by", "unknown")))
    return cands, looked, nones, made


def rank(rec: dict, facts: dict | None) -> float:
    """Which candidates to keep when the host's file limit forces a choice: coverage of the notes our renders are graded on counts most, then
    an open library (known provenance), a Creative Commons licence, how steady the measured pitch is and how much steady sound there is.
    All of it is measured or stated; none of it says anything about timbre."""
    steady = ((facts or {}).get("pitch") or {}).get("fraction_within_30_cents_of_median")
    sound = sum(n.get("clip_t1", 0) - n.get("clip_t0", 0) for n in (facts or {}).get("steady_notes_by_pitch_tracker", []) if isinstance(n, dict))
    return float(3 * min(len((coverage(rec.get("usable_for_render_notes")) or {"graded": []})["graded"]), 4) + (1 if rec["source_type"] == "open-library" else 0) + (1 if is_cc(rec["licence"]) else 0)
                 + 2 * (steady if number(steady) else 0.0) + min(sound, 12.0) / 12.0)


def build_candidate(rec: dict, facts: dict | None, base: Path, out_dir: Path) -> dict:
    x, sr = decode(base / rec["files"]["clip"])
    gain = playback_gain_db(x)
    encode(base / rec["files"]["clip"], out_dir / f"{rec['id']}.mp4", gain)
    keep = ("id", "instrument", "notes", "dynamic", "articulation", "source_type", "url", "title", "channel_or_author", "licence", "licence_evidence", "start_s", "end_s",
            "conditions", "sample_rate_hz", "reason", "sha256_clip")
    return {**{k: rec.get(k) for k in keep}, "render_coverage": coverage(rec.get("usable_for_render_notes")), "creative_commons": is_cc(rec["licence"]), "audio": f"{out_dir.name}/{rec['id']}.mp4", "clip_duration_s": round(len(x) / sr, 3),
            "playback_gain_db": gain, "facts": facts}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--candidates", type=Path, required=True, help="a mirror of trj:/data/tmp/wa-inst/refs-candidates")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--max-per-family", type=int, default=5, help="show at most this many candidates of a family (the best by rank()); the rest are listed as omitted and get no files")
    args = ap.parse_args()
    if not args.candidates.is_dir():
        raise SystemExit(f"no candidates directory {args.candidates}")
    render_notes = {}
    if (args.candidates / "families.json").is_file():
        render_notes = json.loads((args.candidates / "families.json").read_text()).get("families", {})
    known = {fid for fid, _ in FAMILIES}
    families, refused_total = [], 0
    for fid, title in FAMILIES + [(p.name, p.name) for p in sorted(args.candidates.iterdir()) if p.is_dir() and p.name not in known and not p.name.startswith(("_", "."))]:
        fam_dir = args.candidates / fid
        out_dir = args.out / fid
        valid, refused, looked, nones, made = [], [], [], [], []
        if fam_dir.is_dir():
            try:
                raw, looked, nones, made = read_family_files(fam_dir, fid)
            except Refused as e:
                raw = []
                refused.append({"id": "(file)", "problems": e.problems})
                print(f"REFUSED {fid}: {'; '.join(e.problems)}", file=sys.stderr)
            seen: set[str] = set()
            for rec, origin in raw:
                cid = rec.get("id") if isinstance(rec, dict) else None
                problems, facts = validate(rec, fid, fam_dir)
                if isinstance(cid, str) and cid in seen:
                    problems.append("id appears twice in the family")
                seen.add(str(cid))
                if problems:
                    refused.append({"id": str(cid), "file": origin, "problems": problems})
                    print(f"REFUSED {fid}/{cid} ({origin}): {'; '.join(problems)}", file=sys.stderr)
                else:
                    valid.append((rec, facts))
        # The page host takes a limited number of files: show the best few of a family and list the rest as omitted.
        valid.sort(key=lambda rf: (-rank(*rf), rf[0]["id"]))
        shown, left = valid[:args.max_per_family], valid[args.max_per_family:]
        cands = [build_candidate(rec, facts, fam_dir, out_dir) for rec, facts in sorted(shown, key=lambda rf: rf[0]["id"])]
        none = None if cands else (nones[0] if nones else None)
        refused_total += len(refused)
        families.append({"id": fid, "title": title, "candidates": cands, "omitted": sorted(rec["id"] for rec, _ in left), "looked": looked, "none_found": none, "made_by": sorted(set(made)),
                         "render_notes": render_notes.get(fid), "refused": refused})
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "index.json").write_text(json.dumps({"schema": INDEX_SCHEMA, "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "families": families}, ensure_ascii=False, separators=(",", ":")))
    n = sum(len(f["candidates"]) for f in families)
    omitted = sum(len(f["omitted"]) for f in families)
    empty = [f["id"] for f in families if not f["candidates"] and not f["none_found"]]
    print(f"{len(families)} families, {n} candidates shown ({omitted} more omitted), {sum(1 for f in families if f['none_found'])} none-found, {len(empty)} with nothing recorded ({', '.join(empty) or 'none'}), {refused_total} refused -> {args.out}")
    return 2 if refused_total else 0


if __name__ == "__main__":
    raise SystemExit(main())
