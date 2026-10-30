#!/usr/bin/env python3
"""Build the data of the local Instrument Explorer from a mirror of the instrument program's renders and grades.

Usage:
  python3 apps/instrument-explorer/tools/build.py [--mirror DIR] [--out DIR] [--ref origin/main] [--jobs N]

Inputs (nothing is typed by hand except the FAMILIES crosswalk below):
  --mirror  a fresh mirror made by tools/sync.sh (default: the newest ~/.local/share/instrument-explorer/mirror-*):
            analysis/<family>/<note>/{pair.json, ours.wav, real.wav, *.webp, superseded/<UTC>/pair.json}
            (frankensim.instrument-pair.v1, written by the grader decompose.py), piano-pairs/out-*/piano/<note>/ (same
            schema), u4/clips/<family>/<clip id>.wav (screened reference recordings), kit/ (the kit's renders, with the
            commit of each piece in source.json, graded piece by piece by kitgrade.py) and renders/ (extra renders with
            their logs or provenance sidecars).
  the repo  committed files read with `git show <ref>:<path>` (never the shared working tree, which holds other agents'
            unfinished work): data/instruments/descriptions/<family>/status.json, references-screened.json,
            grades/, geometry/violin-cannone-guarneri-1743-esrf2026/body-mode-*.tsv, instrument-claims.json,
            data/listening/*. Beads come from .beads/issues.jsonl (the tracker's own export; read only) and the history
            from `git log <ref>`.

Output (--out, default ~/.local/share/instrument-explorer/site): index.json, process.json, family/<id>.json, audio/**
(FLAC, lossless) and img/** (the grader's spectrogram and scalogram images, copied). The audio is never committed: several
reference recordings are licensed for local use only (see each clip's licence).

Verdicts are recomputed from the raw numbers the way the grader states them: a measure is "inside" when its gap is no
larger than the largest real-vs-real distance in its ruler (loose: this and neighbouring notes, other dynamics and
takes; strict: the same note only).
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import re
import shutil
import struct
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
HOME_DATA = Path.home() / ".local" / "share" / "instrument-explorer"
GIT = "/Users/jemanuel/bin/git" if Path("/Users/jemanuel/bin/git").exists() else "git"
PROGRAM_START = "2026-10-15"

# The one hand-typed table: how each family is named in the beads, the claims register and commit subjects.
# id, title, section, bead labels, claim fillings, subject keywords (whole words, case-insensitive)
FAMILIES = [
    ("clarinet", "Clarinet", "Reeds", ["irp-clarinet"], ["wind-reed"], ["clarinet"]),
    ("saxophone", "Alto saxophone", "Reeds", ["irp-saxophone"], ["wind-reed"], ["saxophone"]),
    ("oboe", "Oboe", "Reeds", ["irp-double-reeds"], ["wind-reed"], ["oboe", "double reeds?"]),
    ("bassoon", "Bassoon", "Reeds", ["irp-double-reeds"], ["wind-reed"], ["bassoon", "double reeds?"]),
    ("duduk", "Duduk", "Reeds", ["irp-duduk"], ["wind-reed"], ["duduk"]),
    ("flute", "Flute", "Flutes and pipes", ["irp-flute"], ["flue", "jet"], ["flute"]),
    ("recorder", "Recorder", "Flutes and pipes", ["irp-recorder"], ["flue", "jet"], ["recorder"]),
    ("tin-whistle", "Tin whistle", "Flutes and pipes", ["irp-tin-whistle"], ["flue", "jet"], ["tin[- ]whistle"]),
    ("pan-pipes", "Pan pipes", "Flutes and pipes", ["irp-pan-pipes"], ["flue", "jet"], ["pan[- ]?pipes?"]),
    ("organ", "Organ pipe", "Flutes and pipes", [], ["flue"], ["organ"]),
    ("trumpet", "Trumpet", "Brass", ["irp-brass"], ["brass", "valve"], ["trumpet", "brass"]),
    ("trombone", "Tenor trombone", "Brass", ["irp-brass"], ["brass"], ["trombone", "brass"]),
    ("horn", "French horn", "Brass", ["irp-horn"], ["brass"], ["french horn", "horn"]),
    ("tuba", "Tuba", "Brass", ["irp-tuba"], ["brass"], ["tuba"]),
    ("violin", "Violin", "Bowed strings", ["irp-violin"], ["bowed", "string-plate-cavity"], ["violin", "cannone", "bowed"]),
    ("viola", "Viola", "Bowed strings", ["irp-viola"], ["bowed"], ["viola", "bowed"]),
    ("cello", "Cello", "Bowed strings", ["irp-cello"], ["bowed"], ["cello", "bowed"]),
    ("double-bass", "Double bass", "Bowed strings", ["irp-double-bass"], ["bowed"], ["double[- ]bass", "bowed"]),
    ("guitar", "Classical guitar", "Plucked strings", ["irp-guitar"], ["string", "string-plate-cavity", "plate"], ["guitar", "esteso"]),
    ("guitar-steel", "Steel-string guitar", "Plucked strings", ["irp-steel-guitar"], ["electric-string"], ["steel[- ]string", "electric guitar"]),
    ("piano", "Piano", "Keyboard", [], ["hammer-felt", "piano-lining", "string-plate-nl"], ["piano", "soundboard", "hammer"]),
    ("xylophone", "Xylophone", "Mallets", ["irp-bars"], ["bar", "structural-modal"], ["xylophone", "bars?-s1"]),
    ("marimba", "Marimba", "Mallets", ["irp-bars"], ["bar"], ["marimba"]),
    ("vibraphone", "Vibraphone", "Mallets", ["irp-vibraphone"], ["bar"], ["vibraphone"]),
    ("timpani", "Timpani", "Drums", [], [], ["timpani", "kettle ?drum"]),
    ("kit", "Drum kit", "Drums", [], [], ["drum", "kit", "snare", "cymbals?", "kick", "toms?"]),
    ("voice", "Voice", "Voice", ["irp-voice"], ["voice-glottis", "voice-tract"], ["voice", "vowels?", "glott\\w*", "vocal tract"]),
]
SECTIONS = ["Reeds", "Flutes and pipes", "Brass", "Bowed strings", "Plucked strings", "Keyboard", "Mallets", "Drums", "Voice"]

# Committed listening fixtures from the August program (before the 2026-10-21 first-principles pivot), by family.
LISTENING_FIXTURES = {
    "brass-emergent": "trumpet", "guitar-strum": "guitar", "piano-ladder": "piano", "reed-sustained": "clarinet",
    "vowel-ladder": "voice", "xylophone-strike": "xylophone", "electric-clean-vs-driven": "guitar-steel",
}

# Every note each family can play today (sync.sh PLAYABLE), in playing order: (title, file stems).
# An empty stem list means "every hold-t<tension> render, by tension": brass pitch emerges from lip tension.
# The titles say only what is played. How well each note lands is measured at build time from the mirror's pitch.json
# (pitch_summary); a hand-typed verdict here went stale as soon as the renders changed.
PLAYABLE = {
    "tin-whistle": ("D-major fingerings, D5 to D6", ["d5", "e5", "fs5", "g5", "a5", "b5", "cs6", "d6"]),
    "pan-pipes": ("Five tubes, D3 to B4", ["d3", "g3", "d4", "g4", "b4"]),
    "recorder": ("Four first-register fingerings, F4 to B4", ["f1", "a1", "bb1", "b1"]),
    "flute": ("Three fingerings: D4, G4, A5", ["d4", "g4", "a5"]),
    "clarinet": ("Chromatic fingerings, concert D3 to Bb5",
                 ["01-D3", "02-Eb3", "03-F3", "04-G3", "05-Ab3", "06-Bb3", "07-C4", "08-D4", "09-Eb4", "10-F4", "11-Gb4",
                  "12-G4", "13-A4", "14-Bb4", "15-C5", "16-D5", "17-Eb5", "18-F5", "19-G5", "20-A5", "21-Bb5"]),
    "voice": ("Vowel /a/ at six pitches, C3 to F4", ["a-C3", "a-Ds3", "a-A3", "a-C4", "a-D4", "a-F4"]),
    "violin": ("Four open strings, G3 D4 A4 E5", [f"violin-arco-mf-open-{n}" for n in ("G3", "D4", "A4", "E5")]),
    "viola": ("Four open strings, C3 G3 D4 A4", [f"viola-arco-mf-open-{n}" for n in ("C3", "G3", "D4", "A4")]),
    "cello": ("Four open strings, C2 G2 D3 A3", [f"cello-arco-mf-open-{n}" for n in ("C2", "G2", "D3", "A3")]),
    "double-bass": ("Four open strings, E1 A1 D2 G2", [f"double-bass-arco-mf-open-{n}" for n in ("E1", "A1", "D2", "G2")]),
    "trumpet": ("Valves up: four lip tensions aimed at Bb3, F4, Bb4 and D5", []),
    "trombone": ("Slide in first position: the partials the lip finds at six lip tensions (no target pitch)", []),
}

UIOWA_CREDIT = "University of Iowa Electronic Music Studios, Musical Instrument Samples"

BADGE_LABEL = {
    "none": "No render yet",
    "listening": "Rendered; no recording to grade against",
    "ungraded": "Rendered, not graded",
    "close": "Inside the real spread",
    "mixed": "Partly inside the real spread",
    "off": "Outside the real spread",
}


# ------------------------------------------------------------------ small helpers
def log(*a):
    print(*a, file=sys.stderr, flush=True)


def git(*args: str) -> bytes:
    res = subprocess.run([GIT, "-C", str(REPO), *args], capture_output=True)
    if res.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)}: {res.stderr.decode(errors='replace')[:300]}")
    return res.stdout


def git_show(ref: str, path: str) -> bytes | None:
    try:
        return git("show", f"{ref}:{path}")
    except RuntimeError:
        return None


def git_ls(ref: str, prefix: str) -> list[str]:
    return [p for p in git("ls-tree", "-r", "--name-only", ref, "--", prefix).decode().splitlines() if p]


def r(v, d=4):
    if isinstance(v, (int, float)) and v == v and abs(v) != float("inf"):
        return round(float(v), d)
    return None


def tsv_rows(text: str) -> list[dict]:
    lines = [l for l in text.splitlines() if l.strip()]
    start = 1 if lines and "\t" not in lines[0] else 0  # first line is a schema id
    head = lines[start].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[start + 1:]]


def wav_info(path: Path) -> dict:
    """Channels, rate and duration from a RIFF/WAVE header (PCM or IEEE float)."""
    try:
        with open(path, "rb") as f:
            if f.read(4) != b"RIFF":
                return {}
            f.read(4)
            if f.read(4) != b"WAVE":
                return {}
            ch = rate = bits = None
            while True:
                hdr = f.read(8)
                if len(hdr) < 8:
                    return {}
                cid, size = hdr[:4], struct.unpack("<I", hdr[4:])[0]
                if cid == b"fmt ":
                    fmt = f.read(size)
                    ch, rate = struct.unpack("<HI", fmt[2:8])
                    bits = struct.unpack("<H", fmt[14:16])[0]
                elif cid == b"data":
                    if not (ch and rate and bits):
                        return {}
                    return {"channels": ch, "rate": rate, "duration_s": round(size / (ch * rate * bits / 8), 3)}
                else:
                    f.seek(size + (size & 1), 1)
    except OSError:
        return {}


class Media:
    """Collects (source, destination) pairs and makes them in parallel: WAV to lossless FLAC, images copied.

    The paths handed back carry ?v=<hash of the source's bytes>: a re-render keeps its path, and audio and images are
    served with a day's max-age, so without the version a visitor would hear the old render for up to a day. The hash is
    of the content, not the mtime, because a phrase join rewrites every phrase, and an unchanged one must not be fetched
    again. Hashes are kept in _cache/versions.json against each source's size and mtime, so a build hashes only what
    changed."""

    def __init__(self, out: Path):
        self.out = out
        self.jobs: dict[Path, Path] = {}
        self.versions_path = out / "_cache" / "versions.json"
        try:
            self.versions = json.loads(self.versions_path.read_text())
        except (OSError, ValueError):
            self.versions = {}
        self.used: set[str] = set()

    def versioned(self, src: Path, rel: str) -> str:
        st = src.stat()
        key = f"{src}|{st.st_size}|{st.st_mtime_ns}"
        if key not in self.versions:
            with open(src, "rb") as f:
                self.versions[key] = hashlib.file_digest(f, "sha1").hexdigest()[:10]
        self.used.add(key)
        return f"{rel}?v={self.versions[key]}"

    def save_versions(self) -> None:
        self.versions_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.versions_path.with_suffix(".tmp")
        tmp.write_text(json.dumps({k: self.versions[k] for k in self.used}, sort_keys=True))
        tmp.replace(self.versions_path)

    def audio(self, src: Path, rel: str) -> str | None:
        if not src.is_file() or src.stat().st_size < 1000:
            return None
        rel = rel if rel.endswith(".flac") else rel + ".flac"
        self.jobs[self.out / rel] = src
        return self.versioned(src, rel)

    def image(self, src: Path, rel: str) -> str | None:
        if not src.is_file():
            return None
        self.jobs[self.out / rel] = src
        return self.versioned(src, rel)

    def run(self, workers: int) -> list[str]:
        failures = []

        def make(item):
            # Each made file takes its source's mtime, and is remade whenever the two differ: a replacement source can be
            # older than the file made from its predecessor (a re-graded directory swapped in), so "newer" is not enough.
            dst, src = item
            st = src.stat()
            if dst.exists() and dst.stat().st_mtime_ns == st.st_mtime_ns:
                return None
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.suffix != ".flac":
                shutil.copy2(src, dst)
                return None
            tmp = dst.with_name(dst.stem + ".tmp.flac")
            res = subprocess.run(["ffmpeg", "-nostdin", "-y", "-v", "error", "-i", str(src), "-c:a", "flac",
                                  "-sample_fmt", "s32", "-compression_level", "5", str(tmp)], capture_output=True)
            if res.returncode != 0:
                return f"{src}: {res.stderr.decode(errors='replace')[:200]}"
            os.utime(tmp, ns=(st.st_atime_ns, st.st_mtime_ns))
            tmp.replace(dst)
            return None

        todo = list(self.jobs.items())
        log(f"media: {len(todo)} files ({sum(1 for d, _ in todo if d.suffix == '.flac')} audio)")
        with cf.ThreadPoolExecutor(max_workers=workers) as ex:
            for i, err in enumerate(ex.map(make, todo)):
                if err:
                    failures.append(err)
                if (i + 1) % 200 == 0:
                    log(f"  media {i + 1}/{len(todo)}")
        return failures


# ------------------------------------------------------------------ grading pairs
def inside(gap, values):
    if gap is None or not values:
        return None
    return bool(gap <= max(values))


def metric_view(m: dict, has_ours: bool) -> dict:
    gap = r(m.get("gap")) if has_ours else None
    raw = m["gap"] if gap is not None else None  # the verdict compares unrounded values, as the grader does
    loose = [e["value"] for e in (m.get("ruler") or []) if r(e.get("value")) is not None]
    strict = None if m.get("ruler_strict") is None else [e["value"] for e in m["ruler_strict"] if r(e.get("value")) is not None]
    return {
        "id": m["id"], "label": m.get("label", m["id"]), "unit": m.get("gap_unit") or m.get("unit", ""),
        "gap": gap, "note": m.get("note", ""),
        "ruler": [{"label": e["label"], "value": r(e["value"])} for e in (m.get("ruler") or []) if r(e.get("value")) is not None],
        "ruler_max": r(max(loose)) if loose else None,
        "ruler_strict_max": r(max(strict)) if strict else None,
        "inside": inside(raw, loose),
        "inside_strict": None if strict is None else inside(raw, strict),
    }


def tally(views: list[dict]) -> dict:
    g = [v for v in views if v["gap"] is not None]
    t = lambda k: {"inside": sum(v[k] is True for v in g), "outside": sum(v[k] is False for v in g), "no_ruler": sum(v[k] is None for v in g)}  # noqa: E731
    return {"measures": len(g), "loose": t("inside"), "strict": t("inside_strict")}


def pair_kind(pj: dict) -> str:
    if pj.get("ours") is None:
        return "refused" if pj.get("refused") else "reference-only"
    return "listening" if pj.get("real") is None else "graded"


def stamp_of(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).strftime("%Y-%m-%dT%H:%MZ")


def sounds_text(s: str) -> str:
    """The pitch check's "G3-1c" as "G3 −1¢"."""
    m = re.fullmatch(r"([A-G][#b]?\d)([+-]\d+)c", s)
    if not m:
        return s
    c = int(m[2])
    return f"{m[1]} {'±' if c == 0 else '+' if c > 0 else '−'}{abs(c)}¢"


def pitch_summary(entries: list[dict]) -> str:
    """How a family's notes land, from the pitch check: "all 21 sound their note, median 9¢ off" or "none sounds its note".
    A note sounds its note when its pitch is within 50 cents of the target; notes without a target are not counted."""
    graded = [e for e in entries if e.get("target")]
    if not graded:
        return ""
    hit = [abs(e["cents_off"]) for e in graded if e.get("cents_off") is not None and abs(e["cents_off"]) <= 50]
    if not hit:
        return "none sounds its note" if len(graded) > 1 else "it does not sound its note"
    hit.sort()
    mid = round(hit[len(hit) // 2] if len(hit) % 2 else (hit[len(hit) // 2 - 1] + hit[len(hit) // 2]) / 2)
    if len(hit) == len(graded):
        return f"all {len(hit)} sound their note, median {mid}¢ off" if len(hit) > 1 else f"it sounds its note, {mid}¢ off"
    return f"{len(hit)} of {len(graded)} sound their note, those a median {mid}¢ off"


def join_wavs(wavs: list[Path], dst: Path) -> bool:
    """Join renders end to end, as 48 kHz mono, into one WAV; remade only when a source is newer."""
    if dst.is_file() and all(dst.stat().st_mtime >= w.stat().st_mtime for w in wavs):
        return True
    dst.parent.mkdir(parents=True, exist_ok=True)
    inputs = [a for w in wavs for a in ("-i", str(w))]
    chain = "".join(f"[{i}:a]aresample=48000,aformat=channel_layouts=mono[a{i}];" for i in range(len(wavs)))
    chain += "".join(f"[a{i}]" for i in range(len(wavs))) + f"concat=n={len(wavs)}:v=0:a=1[joined]"
    res = subprocess.run(["ffmpeg", "-nostdin", "-y", "-v", "error", *inputs, "-filter_complex", chain,
                          "-map", "[joined]", str(dst)], capture_output=True)
    return res.returncode == 0


def history_of(note_dir: Path, current: dict) -> list[dict]:
    snaps = []
    sup = note_dir / "superseded"
    entries = []
    if sup.is_dir():
        for d in sorted(sup.iterdir()):
            pj = d / "pair.json"
            if pj.is_file():
                # The directory is named for when the grading was replaced; the file keeps the time it was made.
                entries.append((stamp_of(pj), pj, False))
    entries.sort(key=lambda e: e[0])
    entries.append((stamp_of(note_dir / "pair.json"), note_dir / "pair.json", True))
    for when, pj_path, is_current in entries:
        try:
            pj = json.loads(pj_path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        kind = pair_kind(pj)
        views = [metric_view(m, kind == "graded") for m in (pj.get("metrics") or [])] if kind == "graded" else []
        snaps.append({
            "when": when, "current": is_current, "kind": kind,
            "git_sha": (pj.get("ours") or {}).get("git_sha") or (pj.get("refused") or {}).get("render_sha"),
            "refused": (pj.get("refused") or {}).get("reason"),
            "tally": tally(views) if views else None,
            "gaps": {v["id"]: v["gap"] for v in views},
            "ratios": {v["id"]: (r(v["gap"] / v["ruler_max"], 3) if v["gap"] is not None and v["ruler_max"] else None) for v in views},
        })
    return snaps


def score_history(items: list[tuple[str, dict]]) -> list[dict]:
    """Share of measures inside the real spread after each grading, every note counted at its latest grading so far.
    A refused grading takes its note out until it is graded again."""
    events = sorted(((s["when"], f"{fid}/{p['id']}", s.get("tally")) for fid, p in items for s in p["history"]),
                    key=lambda e: (e[0], e[1]))
    state: dict[str, tuple[int, int]] = {}
    points = []
    for when, key, t in events:
        if t and t["loose"]["inside"] + t["loose"]["outside"]:
            state[key] = (t["loose"]["inside"], t["loose"]["outside"])
        else:
            state.pop(key, None)
        if state:
            i = sum(a for a, _ in state.values())
            points.append({"when": when, "inside": i, "measured": i + sum(b for _, b in state.values()), "notes": len(state)})
    # One point per hour (its last state): a regrading batch lands many snapshots within minutes.
    by_hour = {p["when"][:13]: p for p in points}
    return list(by_hour.values())


def measure_table(family_docs: dict) -> list[dict]:
    """Every grader measure over every current graded note: how often it is inside, how far off it is (gap over the
    largest real-vs-real distance), and which families miss it. Hardest first."""
    acc: dict[str, dict] = {}
    for doc in family_docs.values():
        for p in doc["pairs"]:
            for m in p["metrics"] if p["kind"] == "graded" else []:
                if m["gap"] is None or m["inside"] is None:
                    continue
                e = acc.setdefault(m["id"], {"id": m["id"], "label": m["label"], "unit": m["unit"], "graded": 0, "inside": 0,
                                             "ratios": [], "missed_by": Counter()})
                e["graded"] += 1
                e["inside"] += m["inside"]
                if m["ruler_max"]:
                    e["ratios"].append(m["gap"] / m["ruler_max"])
                if not m["inside"]:
                    e["missed_by"][doc["title"]] += 1
    rows = []
    for e in acc.values():
        rs = sorted(e.pop("ratios"))
        e["median_ratio"] = r(rs[len(rs) // 2], 3) if rs else None
        e["missed_by"] = dict(e["missed_by"].most_common())
        rows.append(e)
    return sorted(rows, key=lambda e: e["inside"] / e["graded"])


def build_pair(fid: str, note_id: str, note_dir: Path, media: Media, origin: str) -> dict | None:
    pj_path = note_dir / "pair.json"
    try:
        pj = json.loads(pj_path.read_text())
    except (OSError, json.JSONDecodeError) as e:
        log(f"skip {note_dir}: {e}")
        return None
    kind = pair_kind(pj)
    base = f"{origin}/{fid}/{note_id}"
    audio = {}
    for side in ("ours", "real"):
        name = (pj.get("audio") or {}).get(side)
        if name:
            audio[side] = media.audio(note_dir / name, f"audio/{base}/{side}")
    figures = {}
    for kind_name, fig in (pj.get("figures") or {}).items():
        if not isinstance(fig, dict):
            continue
        imgs = {}
        for side in ("real", "ours", "diff"):
            name = fig.get(side)
            if isinstance(name, str):
                imgs[side] = media.image(note_dir / name, f"img/{base}/{name}")
        figures[kind_name] = {"img": imgs, "f_lo_hz": fig.get("f_lo_hz"), "f_hi_hz": fig.get("f_hi_hz"),
                              "span_db": fig.get("span_db"), "diff_range_db": fig.get("diff_range_db")}
    views = [metric_view(m, kind == "graded") for m in (pj.get("metrics") or [])] if kind == "graded" else []
    real = pj.get("real") or None
    ours = pj.get("ours") or None
    src = (real or {}).get("source") or {}
    # The grader's decompose library is the University of Iowa set; its pairs carry no source block.
    library_note = (pj.get("decompose") or {}).get("library_note")
    if real is not None and not src and (str(real.get("reference_id") or "").startswith("uiowa") or library_note):
        src = {"credit": UIOWA_CREDIT, "licence": "UIowa MIS: free for this project's grading use; not redistributed",
               "url": "https://theremin.music.uiowa.edu/MIS.html"}
        if not real.get("file") and library_note:
            real = {**real, "file": library_note.split("@")[0]}
    parts = pj.get("partials") or {}
    f0 = parts.get("f0_hz") if isinstance(parts, dict) else None
    return {
        "id": note_id, "note": pj.get("note", note_id), "dynamic": pj.get("dynamic"), "kind": kind,
        "duration_s": r(pj.get("duration_s"), 3), "origin": origin,
        "audio": audio, "level_match": (pj.get("audio") or {}).get("level_match"),
        "figures": figures,
        "ours": None if ours is None else {"git_sha": ours.get("git_sha"), "observer": ours.get("observer"),
                                            "recipe": ours.get("recipe"), "sample_rate": ours.get("sample_rate")},
        "real": None if real is None else {"reference_id": real.get("reference_id"), "file": real.get("file"),
                                            "segment_s": real.get("segment_s"), "credit": src.get("credit"),
                                            "licence": src.get("licence"), "url": src.get("url")},
        "no_reference": (pj.get("no_reference") or {}).get("reason"),
        "refused": pj.get("refused"),
        "f0_hz": f0 if isinstance(f0, dict) else None,
        "metrics": views,
        "tally": tally(views) if views else None,
        "history": history_of(note_dir, pj),
    }


SALAMANDER = {"credit": "Salamander Drumkit by Alexander Holm", "licence": "CC BY-SA 3.0",
              "url": "https://archive.org/details/SalamanderDrumkit"}
# 96ccf8f1b (2026-10-14) adopted the kit edits that graded closer to these Salamander strokes; the card route still uses them
# (ruling 7486 replaces them). Shown with every kit comparison until it lands, because it compares against the same strokes.
KIT_FITTED_CAVEAT = ("Several of this kit's settings were chosen on 2026-10-14 (commit 96ccf8f1b) because they made it sound "
                     "closer to these same Salamander recordings: the head tensions, strike points, the stick's contact loss, "
                     "the kick beater's felt, the snare wires' tension, preload and contact, the ride's shell and the hi-hat's "
                     "rim chatter. That is fitting, which this project does not allow, so closeness here is not evidence. "
                     "They are being replaced with measured values or declared estimates.")
# verisim.py's diverging map (its DIV stops), which differs from the grader's: the key must show the colours drawn
VERISIM_DIFF_KEY = ("linear-gradient(90deg, rgb(58,120,190), rgb(24,40,60) 42%, rgb(14,15,18) 50%, "
                    "rgb(70,34,20) 58%, rgb(232,106,60))")


def kit_commit(src: dict, pid: str) -> str:
    """The commit a kit piece was rendered at: source.json's "commit", unless "pieces" names a later one for that piece
    (a piece re-rendered after its own physics landed). The groove carries every piece, so it names them all."""
    pieces = src.get("pieces") or {}
    if pid == "groove" and pieces:
        return f"{src['commit']} ({', '.join(f'{p} at {c}' for p, c in pieces.items())})"
    return pieces.get(pid, src["commit"])


def kit_pairs(kit: Path, media: Media) -> list[dict]:
    """The kit grade (kitgrade.py on renders at known commits): each piece beside one Salamander stroke, level-matched,
    with the grader's descriptors and its real-vs-real ruler. Compared, not graded: the kit has no grade gate."""
    try:
        grade = json.loads((kit / "grade" / "data" / "kit.json").read_text())
        src = json.loads((kit / "source.json").read_text())
    except (OSError, json.JSONDecodeError) as e:
        log(f"kit grade: {e}")
        return []
    out = []
    for g in grade["pairs"]:
        pid, d, im = g["id"], g["descriptors"], g["images"]
        base = f"kit9/{pid}"
        audio = {side: media.audio(kit / "grade" / Path(g["audio"][side]).with_suffix(".wav"), f"audio/{base}/{side}") for side in ("ours", "real")}
        imgs = {side: media.image(kit / "grade" / im[side], f"img/{base}/{side}.webp") for side in ("real", "ours", "diff")}
        out.append({
            # "14x6" -> "14×6": Inter's contextual alternates already draw the list's x as a times sign, Fraunces does not
            "id": f"kit9-{pid}", "note": re.sub(r"(\d)x(\d)", "\\1×\\2", g["title"]), "dynamic": None, "kind": "compared",
            "duration_s": g["duration_s"],
            "origin": "kit9", "audio": {k: v for k, v in audio.items() if v}, "level_match": g["audio"]["level_match"],
            "figures": {"spectrogram": {"img": {k: v for k, v in imgs.items() if v}, "f_lo_hz": im["f_lo"], "f_hi_hz": im["f_hi"],
                                        "span_db": im["span_db"], "diff_range_db": im["diff_range_db"], "diff_key": VERISIM_DIFF_KEY}},
            "ours": {"git_sha": kit_commit(src, pid), "observer": g["ours"]["label"], "recipe": f"kit {pid} (fs-vkplate)", "sample_rate": 48000},
            "real": {"reference_id": g["real"]["file"], "file": f"{g['real']['file']}.wav", "segment_s": None, **SALAMANDER},
            "no_reference": None, "refused": None,
            "f0_hz": {"ours": d["fundamental_hz"]["ours"], "real": d["fundamental_hz"]["real"]},
            "metrics": [], "tally": None, "history": [],
            "compare": {"descriptors": {k: {"ours": v["ours"], "real": v["real"], "note": v["note"]} for k, v in d.items()
                                        if k in ("centroid_hz", "above_1k_pct", "above_4k_pct")},
                        "balance_gap_db": d["balance_gap_db"]["rms_gap_db"], "ruler": g["ruler"],
                        "caveat": src.get("caveat", KIT_FITTED_CAVEAT)},
        })
    return out


def kit_decay_evidence(kit: Path) -> list[dict]:
    """Per-octave 20 dB decays of the kit renders against their Salamander strokes (kitdecays.py), drums and cymbals
    apart so each chart keeps one colour per piece."""
    try:
        d = json.loads((kit / "grade" / "data" / "decays.json").read_text())
        src = json.loads((kit / "source.json").read_text())
    except (OSError, json.JSONDecodeError):
        return []
    cymbal = lambda p: p["id"].startswith(("hihat", "crash", "ride"))  # noqa: E731
    how = ("Time for each octave band to fall 20 dB after its peak. Grey dots are two strokes of the real one at the same "
           "dynamic. A triangle means the render was still ringing when the excerpt ended, so its value is a lower bound: "
           "the time from the band's peak to the end of the excerpt. ")
    charts = (
        # a kit whose renders carry a sourced loss term says so in source.json's "drum_loss"
        ("Drum decay per octave band: renders against recordings", lambda p: not cymbal(p),
         how + src.get("drum_loss", "The heads' film and edge loss enter these renders as a stated unknown (zero). No "
                       "measurement of a drumhead's edge, hoop or shell loss has been found, and fitting one to these "
                       "recordings is not allowed, so the drums ring longer than the real ones.")),
        ("Cymbal decay per octave band: renders against recordings", cymbal,
         how + "The cymbals' loss is still each preset's declared decay time at 500 Hz and 10 kHz, an estimate rather than "
               "material data; lane .372 is replacing it with the felt clamp and air loss from data."),
    )
    out = []
    for title, keep, caption in charts:
        pieces = [{"piece": p["id"], "what": f"{p['title']}; recording {p['real']['id']}, excerpt {p['duration_s']} s",
                   "render_decays": p["ours"], "references": [p["real"], p["real_again"]]} for p in d["pieces"] if keep(p)]
        if pieces:
            out.append({"kind": "band-decays", "title": title, "bands_hz": d["bands_hz"], "pieces": pieces, "caption": caption,
                        "source": f"kitdecays.py on renders at commit {kit_commit(src, 'groove')} and the Salamander Drumkit overheads ({d['method']})"})
    return out


def judge(pairs: list[dict], status: dict | None) -> dict:
    graded = [p for p in pairs if p["kind"] == "graded"]
    listening = [p for p in pairs if p["kind"] == "listening"]
    refused = [p for p in pairs if p["kind"] == "refused"]
    counts = {"graded": len(graded), "listening": len(listening), "refused": len(refused),
              "compared": sum(p["kind"] == "compared" for p in pairs),
              "reference_only": sum(p["kind"] == "reference-only" for p in pairs)}
    li = sum(p["tally"]["loose"]["inside"] for p in graded)
    lo = sum(p["tally"]["loose"]["outside"] for p in graded)
    si = sum(p["tally"]["strict"]["inside"] for p in graded)
    so = sum(p["tally"]["strict"]["outside"] for p in graded)
    counts.update(inside=li, outside=lo, strict_inside=si, strict_outside=so)
    if not graded:
        if listening:
            badge = "listening"
        elif refused or counts["compared"]:
            badge = "ungraded"
        else:
            badge = "none"
    elif li + lo == 0:
        badge = "ungraded"
    else:
        level = lambda i, o: "close" if o == 0 else "mixed" if i / (i + o) >= 0.5 else "off"  # noqa: E731
        levels = [level(li, lo)] + ([level(si, so)] if si + so else [])
        badge = min(levels, key=("off", "mixed", "close").index)
    return {"badge": badge, "label": BADGE_LABEL[badge], "counts": counts, "blocked": (status or {}).get("blocked")}


def derived_summary(counts: dict, n_extras: int, n_fixtures: int, n_refs: int) -> str:
    """A factual stand-in for families the program's status file does not describe: what is on record, nothing more."""
    plural = lambda n, word: f"{n} {word}{'' if n == 1 else 's'}"  # noqa: E731
    facts = []
    if counts["graded"]:
        measured = counts["inside"] + counts["outside"]
        facts.append(f"{plural(counts['graded'], 'note')} graded against a recording ({counts['inside']} of {measured} measures inside the real spread)")
    if counts["listening"]:
        facts.append(f"{plural(counts['listening'], 'listening pair')} not graded")
    if counts["refused"]:
        facts.append(f"{plural(counts['refused'], 'render')} the grader refused")
    if counts.get("compared"):
        facts.append(f"{plural(counts['compared'], 'piece')} compared with a recording, not graded")
    if n_extras:
        facts.append(f"{plural(n_extras, 'other render')} not graded")
    if n_fixtures:
        facts.append(f"{plural(n_fixtures, 'listening fixture')} from before the 2026-10-21 pivot")
    if n_refs:
        facts.append(f"{plural(n_refs, 'screened reference recording')}")
    on_record = "; ".join(facts) if facts else "nothing yet"
    return f"The program's status file does not describe this family yet. On record: {on_record}."


# ------------------------------------------------------------------ references
def build_references(refs: dict, mirror: Path, media: Media) -> dict[str, list[dict]]:
    lic = {l["id"]: l for l in refs.get("licences", [])}
    out: dict[str, list[dict]] = defaultdict(list)
    for f in refs["files"]:
        # Archive members carry no licence or author of their own: every clip is credited from its file `f`.
        for clips in [f.get("clips") or []] + [m.get("clips") or [] for m in f.get("members") or []]:
            for c in clips:
                if c.get("status") != "VALID":
                    continue
                src = mirror / "u4" / "clips" / c["family"] / f"{c['id']}.wav"
                rel = media.audio(src, f"audio/refs/{c['family']}/{c['id']}")
                if rel is None:
                    continue
                L = lic.get(f.get("licence"), {})
                info = wav_info(src)
                out[c["family"]].append({
                    "id": c["id"], "note": c.get("note"), "piece": c.get("piece"), "dynamic": c.get("dynamic"),
                    "kind": c.get("kind"), "audio": rel, "duration_s": info.get("duration_s"),
                    "rate": info.get("rate"), "channels": info.get("channels"),
                    "instrument": f.get("instrument"), "source": f.get("id"),
                    "url": f.get("source_page") or f.get("url"), "author": f.get("author"),
                    "licence": L.get("name", f.get("licence")), "licence_id": f.get("licence"),
                    "licence_url": L.get("url"), "licence_notes": L.get("notes"),
                    "credit": L.get("credit") or (UIOWA_CREDIT if f.get("licence") == "uiowa-mis" else None),
                })
    for fam in out:
        out[fam].sort(key=lambda c: (note_key(c.get("note")), c.get("dynamic") or "", c["id"]))
    return out


NOTE_RE = re.compile(r"^([A-G])([#b]?)(-?\d)")


def note_key(n):
    m = NOTE_RE.match(n or "")
    if not m:
        return (99, n or "")
    semis = {"C": 0, "D": 2, "E": 4, "F": 5, "G": 7, "A": 9, "B": 11}[m.group(1)] + {"#": 1, "b": -1, "": 0}[m.group(2)]
    return (int(m.group(3)) * 12 + semis, n)


# ------------------------------------------------------------------ beads and history
def load_beads() -> list[dict]:
    path = REPO / ".beads" / "issues.jsonl"
    beads = []
    with open(path) as fh:
        for line in fh:
            try:
                d = json.loads(line)
            except json.JSONDecodeError:
                continue
            if "instrument-realism" in d.get("id", "") or "irp" in (d.get("labels") or []):
                if d.get("status") == "tombstone":
                    continue
                beads.append(d)
    return beads


def word_re(words: list[str]) -> re.Pattern:
    return re.compile(r"\b(" + "|".join(words) + r")\b", re.I)


def blockers(beads_by_id: dict, bid: str) -> set[str]:
    """Transitive prerequisites: the beads that block bid, and theirs."""
    seen, stack = set(), [bid]
    while stack:
        b = beads_by_id.get(stack.pop())
        for dep in (b or {}).get("dependencies") or []:
            if dep.get("type") == "blocks" and dep.get("issue_id") == b["id"]:
                x = dep.get("depends_on_id")
                if x and x not in seen and x in beads_by_id:
                    seen.add(x)
                    stack.append(x)
    return seen


MILESTONE_RE = re.compile(r"\bM-([A-Z]+)(?:-(QP\d|Q\d)|(0))\b")
GATE_TEXT = {"QP1": "Passive chain + speaks", "Q2": "Recognisable", "Q3": "Inside the ruler", "Q4": "Blind tests",
             "Q5": "Phrases", "0": "First audible baseline"}


def milestones_for(labels: list[str], beads: list[dict], beads_by_id: dict) -> list[dict]:
    out = []
    for b in beads:
        if not set(labels) & set(b.get("labels") or []):
            continue
        m = MILESTONE_RE.search(b.get("title", ""))
        if not m:
            continue
        gate = m.group(2) or m.group(3)
        pre = blockers(beads_by_id, b["id"])
        st = Counter(beads_by_id[x]["status"] for x in pre)
        out.append({"id": b["id"], "code": f"M-{m.group(1)}-{gate}" if m.group(2) else f"M-{m.group(1)}0", "gate": gate,
                    "gate_text": GATE_TEXT.get(gate, gate),
                    "title": re.sub(r"^\[IRP\]\s*", "", b["title"]), "status": b["status"], "closed_at": b.get("closed_at"),
                    "prereq": {"total": len(pre), "closed": st.get("closed", 0), "in_progress": st.get("in_progress", 0),
                               "open": st.get("open", 0) + st.get("blocked", 0)}})
    order = ["0", "QP1", "Q2", "Q3", "Q4", "Q5"]
    out.sort(key=lambda x: (x["code"].split("-")[1] if x["code"].count("-") > 1 else x["code"], order.index(x["gate"]) if x["gate"] in order else 9))
    return out


def load_commits(ref: str) -> list[dict]:
    # a bare date in --since means that day at the current time of day, so the boundary moved with each build's hour
    raw = git("log", ref, f"--since={PROGRAM_START} 00:00:00 +0000", "--format=%h%x1f%at%x1f%s%x1e").decode(errors="replace")
    out = []
    for rec in raw.split("\x1e"):
        rec = rec.strip()
        if not rec:
            continue
        sha, stamp, subject = (rec.split("\x1f") + ["", ""])[:3]
        when = datetime.fromtimestamp(int(stamp), timezone.utc).strftime("%Y-%m-%d %H:%M")
        tags = re.findall(r"\[(frankensim-[\w.-]+)\]", subject)
        kind = re.match(r"^(\w+)(?:\(([^)]*)\))?!?:", subject)
        out.append({"sha": sha, "when": when, "date": when[:10], "subject": re.sub(r"\s*\[frankensim-[\w.-]+\]", "", subject).strip(),
                    "beads": tags, "type": kind.group(1) if kind else "other", "scope": (kind.group(2) or "") if kind else ""})
    return out


# ------------------------------------------------------------------ qualification evidence (committed tables)
def violin_body_modes(ref: str) -> dict | None:
    base = "data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026/"
    pre, outc = git_show(ref, base + "body-mode-preregistration.tsv"), git_show(ref, base + "body-mode-outcome.tsv")
    if not pre or not outc:
        return None
    bands = [row for row in tsv_rows(pre.decode()) if row.get("quantity") == "frequency"]
    rows = tsv_rows(outc.decode())
    runs: dict[str, dict] = {}
    for row in rows:
        run = runs.setdefault(row["run"], {"run": row["run"], "modes": {}})
        f = row.get("candidate_f_hz")
        run["modes"][row["id"]] = {"f_hz": r(float(f), 2) if f not in (None, "", "none") else None,
                                   "outcome": row.get("registered_rule_outcome"), "note": row.get("note", "")}
    return {"kind": "body-modes", "title": "il Cannone body modes against the pre-registered bands",
            "source": base + "body-mode-outcome.tsv",
            "caption": "Each run of the assembled shell body (plates, ribs, bass bar, soundpost, enclosed air) is a row; each "
                       "registered mode is placed where the identification rules found it. The shaded bands are the "
                       "pre-registered ranges for finished violins (with neck and strings); runs whose name starts with "
                       "'cracked' used a superseded table.",
            "bands": [{"id": b["id"], "lo": float(b["lo"]), "hi": float(b["hi"])} for b in bands],
            "runs": list(runs.values())}


def bowed_sweeps(ref: str, fid: str) -> list[dict]:
    out = []
    for path in git_ls(ref, f"data/instruments/grades/{fid}/"):
        if "bowed-card-sweep" not in path or not path.endswith(".jsonl"):
            continue
        rows = []
        for line in (git_show(ref, path) or b"").decode().splitlines():
            try:
                d = json.loads(line)
            except json.JSONDecodeError:
                continue
            rows.append({"beta": r(d.get("beta")), "force_n": r(d.get("force_n"), 5), "speed": r(d.get("speed_m_s")),
                         "slip": r(d.get("slip_fraction")), "cls": d.get("class"), "peak_hz": r(d.get("peak_hz"), 2),
                         "f1_hz": r(d.get("f1_hz"), 2), "sweep": d.get("sweep")})
        if rows:
            out.append({"kind": "bowed-sweep", "title": f"Bow force sweep ({Path(path).stem.replace('bowed-card-sweep', '').strip('-') or 'main'})",
                        "source": path, "caption": "Each point is one bowed render of the family's performance card at a bow "
                        "force and bow-bridge distance; colour is the regime the render fell into (Helmholtz motion is the "
                        "playable one).", "rows": rows})
    return out


def bars_before_after(ref: str) -> dict | None:
    rows = []
    for name in ("bar-s1-before", "bar-s1-after"):
        raw = git_show(ref, f"data/instruments/grades/bars/{name}.json")
        if not raw:
            continue
        d = json.loads(raw.decode().replace("NaN", "null"))
        dist = d.get("distance") or {}
        rows.append({"name": name, "label": d.get("label"), "target": d.get("target"),
                     "ours_f0": r(d.get("ours_f0"), 2), "real_f0": r((d.get("real") or {}).get("f0"), 2),
                     "distance": {k: r(v, 3) for k, v in dist.items() if isinstance(v, (int, float))}})
    if not rows:
        return None
    return {"kind": "before-after", "title": "Xylophone bar: before and after the 3-D designed undercut",
            "source": "data/instruments/grades/bars/", "caption": "Grader distances of one bar note (Db5) against the "
            "University of Iowa recording, before and after the bar was modelled as a 3-D designed undercut with material "
            "contact. Smaller is closer.", "rows": rows}


def claims_for(claims: list[dict], fillings: list[str]) -> list[dict]:
    out = []
    for c in claims:
        if c.get("filling") in fillings:
            out.append({"filling": c.get("filling"), "image": c.get("image"), "qoi": c.get("qoi"), "gate": c.get("gate"),
                        "exactness": c.get("exactness"), "live_default": c.get("live_default"),
                        "evidence": [e.get("ref") for e in c.get("evidence") or []][:6]})
    return out


# ------------------------------------------------------------------ main
def newest_mirror() -> Path:
    mirrors = sorted(HOME_DATA.glob("mirror-*"))
    if not mirrors:
        raise SystemExit(f"no mirror in {HOME_DATA}; run tools/sync.sh first")
    return mirrors[-1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--mirror", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=HOME_DATA / "site")
    ap.add_argument("--ref", default="origin/main", help="git ref for committed tables and history (default origin/main)")
    ap.add_argument("--jobs", type=int, default=max(2, (os.cpu_count() or 4) - 2))
    ap.add_argument("--phrases", type=Path, default=None,
                    help="a phrasekit output (phrases.json + WAV + webp); default the newest phrases-<UTC> next to the mirrors")
    args = ap.parse_args()
    mirror = args.mirror or newest_mirror()
    # Every input below is optional when absent, so a mistyped path would build a site with no grades at all
    if not (mirror / "analysis").is_dir():
        raise SystemExit(f"{mirror} is not a mirror (no analysis/ in it); run tools/sync.sh or pass --mirror")
    if args.phrases is not None and not (args.phrases / "phrases.json").is_file():
        raise SystemExit(f"{args.phrases} has no phrases.json")
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    media = Media(out)
    ref_sha = git("rev-parse", "--short=9", args.ref).decode().strip()
    log(f"mirror {mirror}\nrepo {REPO} at {args.ref} ({ref_sha})\nout {out}")

    refs = json.loads((git_show(args.ref, "data/instruments/references-screened.json") or b"{}").decode())
    claims_doc = json.loads((git_show(args.ref, "instrument-claims.json") or b"{}").decode() or "{}")
    claims = claims_doc.get("rows") or []
    baseline_path = mirror / "analysis" / "cross-library-baseline.json"
    baseline = json.loads(baseline_path.read_text()).get("families", {}) if baseline_path.is_file() else {}
    beads = load_beads()
    beads_by_id = {b["id"]: b for b in beads}
    commits = load_commits(args.ref)
    references = build_references(refs, mirror, media)
    coverage = refs.get("coverage") or {}

    # Committed listening fixtures (August program), by family.
    fixtures: dict[str, list[dict]] = defaultdict(list)
    for stem, fid in LISTENING_FIXTURES.items():
        wav = git_show(args.ref, f"data/listening/{stem}.wav")
        if not wav:
            continue
        cache = out / "_cache" / f"{stem}.wav"
        cache.parent.mkdir(parents=True, exist_ok=True)
        if not cache.exists() or cache.read_bytes() != wav:
            cache.write_bytes(wav)
        receipt = (git_show(args.ref, f"data/listening/{stem}.listening-receipt") or b"").decode()
        prov = json.loads((git_show(args.ref, f"data/listening/{stem}.provenance.json") or b"{}").decode() or "{}")
        q = dict(l.split("\t", 1) for l in receipt.splitlines() if "\t" in l)
        last = git("log", "-1", "--format=%as", args.ref, "--", f"data/listening/{stem}.wav").decode().strip()
        fixtures[fid].append({"id": stem, "title": stem.replace("-", " "), "audio": media.audio(cache, f"audio/fixtures/{stem}"),
                              "question": q.get("question"), "verdict": q.get("verdict"), "date": last,
                              "era": "August 2026 listening fixture (before the 2026-10-21 first-principles pivot)",
                              "provenance": {k: prov.get(k) for k in ("fixture", "sample_rate", "git_sha", "command") if k in prov}})

    # Extra renders: ungraded violin, flute and timpani renders.
    extras: dict[str, list[dict]] = defaultdict(list)
    va = mirror / "renders" / "violin-assembled"
    for wav in sorted(va.glob("*.wav")) if va.is_dir() else []:
        extras["violin"].append({"id": wav.stem, "title": "Assembled violin, G3 (2 s)", "audio": media.audio(wav, f"audio/extra/violin/{wav.stem}"),
                                 "date": stamp_of(wav), "era": "Speed-lane render (.358) of the assembled instrument: strings, il Cannone body "
                                                               "and enclosed air; not graded"})
    fl = mirror / "renders" / "flute-exterior"
    if fl.is_dir():
        for wav in sorted(fl.glob("*/*.wav")):
            extras["flute"].append({"id": wav.stem, "title": wav.stem, "audio": media.audio(wav, f"audio/extra/flute/{wav.stem}"),
                                    "date": stamp_of(wav), "era": f"Ungraded render with the flute's exterior body (commit {wav.parent.name[:9]})"})
    tp = mirror / "renders" / "timpani"
    for wav in sorted(tp.glob("*/*.wav")) if tp.is_dir() else []:
        extras["timpani"].append({"id": wav.stem, "title": "Kettle B, one strike (0.5 s)", "audio": media.audio(wav, f"audio/extra/timpani/{wav.stem}"),
                                  "date": stamp_of(wav), "era": f"timpani/render.sh at commit {wav.parent.name[:9]}: partial ratios checked "
                                                                "against the model's own coupled modes, not graded against a recording"})

    # Every note a family can play today, joined into one phrase, then each note on its own.
    playable = mirror / "renders" / "playable"
    pitch_doc = json.loads((playable / "pitch.json").read_text()) if (playable / "pitch.json").is_file() else {}
    pitch = {(n["family"], n["stem"]): n for n in pitch_doc.get("notes") or []}
    pitch_table: dict[str, list[dict]] = {}
    pitch_sha: dict[str, str] = {}
    for fid, (title, stems) in PLAYABLE.items():
        d = playable / fid
        if not d.is_dir():
            continue
        if stems:
            notes = [d / f"{s}.wav" for s in stems if (d / f"{s}.wav").is_file()]
        else:
            notes = sorted(d.glob("hold-t*.wav"), key=lambda w: float(w.stem.removeprefix("hold-t").split("-p")[0]))
        if not notes:
            continue
        # each family's notes may come from a different commit: its own source.json first, then the set's
        src = d / "source.json" if (d / "source.json").is_file() else playable / "source.json"
        source = json.loads(src.read_text()) if src.is_file() else {}
        # identical_at: a later origin binary re-rendered one of these inputs and matched it byte for byte;
        # engine: what rendered them when it is not music_render (the assembled violin is the material_strings example)
        same = source.get("identical_at") or {}
        era = "Listening render, not graded: everything this family can play today" + (
            f" ({source.get('engine', 'music_render')} at commit {source['sha']}" + (
                f"; a note re-rendered at {same['sha']} on {same['date']} came out byte-identical" if same.get("sha") else "")
            + ")" if source.get("sha") else "")
        if source.get("sha"):
            pitch_sha[fid] = source["sha"]
        joined = out / "_cache" / "playable" / f"{fid}.wav"
        if join_wavs(notes, joined):
            missing = f", {len(stems) - len(notes)} not rendered" if stems and len(notes) < len(stems) else ""
            landed = pitch_summary([pitch.get((fid, w.stem)) or {} for w in notes])
            extras[fid].append({"id": "playable", "title": f"{title} ({len(notes)} notes joined{missing}{f'; {landed}' if landed else ''})",
                                "audio": media.audio(joined, f"audio/extra/{fid}/playable"), "date": stamp_of(notes[-1]), "era": era})
        for wav in notes:
            pm = pitch.get((fid, wav.stem)) or {}
            label = pm.get("target") or wav.stem
            heard = f": sounds {sounds_text(pm['sounds'])}" if pm.get("sounds") else ""
            # the two pitch trackers disagree and the reading of record misses the asked note: no single clear pitch. (When the
            # reading of record IS the asked note, the other tracker has taken a harmonic, as on the bowed bass's open strings.)
            alt = pm.get("sounds_alt") if not (pm.get("cents_off") is not None and abs(pm["cents_off"]) <= 50) else None
            if heard and alt:
                heard += f" (a second pitch tracker reads {sounds_text(alt)})"
            extras[fid].append({"id": f"playable-{wav.stem}", "title": f"Note {label}{heard}", "pitch": pm or None,
                                "audio": media.audio(wav, f"audio/extra/{fid}/playable-{wav.stem}"), "date": stamp_of(wav), "era": era})
            if pm:
                pitch_table.setdefault(fid, []).append({"label": label, "sounds_alt": alt,
                                                        **{k: pm.get(k) for k in ("target", "sounds", "f0_hz", "cents_off")}})

    # Phrases: scales, double stops, melodies, grooves and rows of notes, joined by phrasekit with a spectrogram and
    # every note's measured pitch.
    phrase_dir = args.phrases or max(HOME_DATA.glob("phrases-*"), default=None)
    phrases: dict[str, list[dict]] = defaultdict(list)
    if phrase_dir and (phrase_dir / "phrases.json").is_file():
        for p in json.loads((phrase_dir / "phrases.json").read_text())["phrases"]:
            p["audio"] = media.audio(phrase_dir / p["audio"], f"audio/phrases/{p['id']}")
            p["map"] = media.image(phrase_dir / p["map"], f"img/phrases/{p['id']}.webp")
            if not (p["audio"] and p["map"]):
                log(f"phrase {p['id']}: no audio or map in {phrase_dir.name}, left out")
                continue
            if p.get("real"):
                p["real"]["audio"] = media.audio(phrase_dir / p["real"]["audio"], f"audio/phrases/{p['id']}-real")
                p["real"]["map"] = media.image(phrase_dir / p["real"]["map"], f"img/phrases/{p['id']}-real.webp")
                if not (p["real"]["audio"] and p["real"]["map"]):
                    p.pop("real")
            if p.get("ref"):  # the phrase rebuilt from real single-note recordings (analyze.py)
                p["ref"]["audio"] = media.audio(phrase_dir / p["ref"]["audio"], f"audio/phrases/{p['id']}-ref")
                p["ref"]["map"] = media.image(phrase_dir / p["ref"]["map"], f"img/phrases/{p['id']}-ref.webp")
                if not (p["ref"]["audio"] and p["ref"]["map"]):
                    p.pop("ref")
            an = phrase_dir / (p.pop("analysis", None) or f"{p['id']}.analysis.json")
            if an.is_file():
                (out / "phrase").mkdir(parents=True, exist_ok=True)
                (out / "phrase" / f"{p['id']}.json").write_text(an.read_text())
                p["analysis"] = f"phrase/{p['id']}.json"
            phrases[p["family"]].append(p)
        log(f"phrases {phrase_dir.name}: {sum(map(len, phrases.values()))}")

    # Grading pairs: the analysis tree plus the piano pairs.
    pair_dirs: dict[str, list[tuple[str, Path, str]]] = defaultdict(list)
    adir = mirror / "analysis"
    for fam_dir in sorted(p for p in adir.iterdir() if p.is_dir()) if adir.is_dir() else []:
        for nd in sorted(p for p in fam_dir.iterdir() if p.is_dir() and (p / "pair.json").is_file()):
            pair_dirs[fam_dir.name].append((nd.name, nd, "analysis"))
    for out_dir in sorted((mirror / "piano-pairs").glob("out-*")):
        dyn = out_dir.name.split("-")[-1]
        for nd in sorted((out_dir / "piano").glob("*")):
            if (nd / "pair.json").is_file():
                pair_dirs["piano"].append((f"{nd.name}-{dyn}", nd, "piano-pairs"))

    families_index, family_docs = [], {}
    fam_beads_all: dict[str, set] = {}
    for fid, title, section, labels, fillings, words in FAMILIES:
        status_raw = git_show(args.ref, f"data/instruments/descriptions/{fid}/status.json")
        status = json.loads(status_raw.decode()) if status_raw else None
        draft = HERE.parent / "status-drafts" / f"{fid}.json"  # a card awaiting its commit; it says so on the page
        if draft.is_file():
            status = json.loads(draft.read_text())
            if not status.get("draft"):
                sys.exit(f"{draft}: a status draft must say why it is not the committed card (its 'draft' field)")
        pairs = [p for p in (build_pair(fid, nid, nd, media, origin) for nid, nd, origin in pair_dirs.get(fid, [])) if p]
        pairs.sort(key=lambda p: ({"graded": 0, "listening": 1, "refused": 2, "reference-only": 3}[p["kind"]], note_key(p["note"]), p["id"]))
        kit_dir = mirror / "kit" if (mirror / "kit").is_dir() else mirror / "kit9"  # mirrors before 2026-10-26T03:30Z
        if fid == "kit":
            pairs += kit_pairs(kit_dir, media)
        verdict = judge(pairs, status)
        if verdict["badge"] == "none" and (extras.get(fid) or fixtures.get(fid)):
            verdict.update(badge="ungraded", label=BADGE_LABEL["ungraded"])
        rx = word_re(words)
        fam_beads = {b["id"] for b in beads if set(labels) & set(b.get("labels") or []) or rx.search(b.get("title", ""))}
        fam_beads_all[fid] = fam_beads
        timeline = [c for c in commits if fam_beads & set(c["beads"]) or rx.search(c["subject"])]
        closes = [{"id": b["id"], "title": re.sub(r"^\[IRP\]\s*", "", b["title"]), "closed_at": b.get("closed_at"),
                   "reason": (b.get("close_reason") or "")[:400]}
                  for b in beads if b["id"] in fam_beads and b.get("status") == "closed" and b.get("closed_at")]
        closes.sort(key=lambda x: x["closed_at"], reverse=True)
        bead_status = Counter(beads_by_id[x]["status"] for x in fam_beads)
        ms = milestones_for(labels, beads, beads_by_id)
        evidence = []
        if fid == "violin":
            v = violin_body_modes(args.ref)
            if v:
                evidence.append(v)
        if section == "Bowed strings":
            evidence += bowed_sweeps(args.ref, fid)
        if fid == "xylophone":
            b = bars_before_after(args.ref)
            if b:
                evidence.append(b)
        if fid == "kit":
            evidence += kit_decay_evidence(kit_dir)
        base = baseline.get(fid)
        cov = coverage.get(fid) or {}
        doc = {
            "id": fid, "title": title, "section": section, "status": status, "verdict": verdict,
            "pairs": pairs, "references": references.get(fid, []),
            "reference_coverage": {"rows": cov.get("rows"), "valid": cov.get("valid"), "by_dynamic": cov.get("by_dynamic"),
                                   "failed_checks": cov.get("failed_checks")},
            "baseline": None if not base else {
                "source": (base.get("source") or {}).get("name"), "compared_with": base.get("compared_with"),
                "recordings": base.get("recordings"), "real": (base.get("real_recordings") or {}).get("loose"),
                "real_strict": (base.get("real_recordings") or {}).get("strict"),
                "ours": (base.get("our_renders_same_notes") or {}).get("loose"),
                "ours_strict": (base.get("our_renders_same_notes") or {}).get("strict")},
            "fixtures": fixtures.get(fid, []), "extras": extras.get(fid, []), "phrases": phrases.get(fid, []),
            "milestones": ms, "beads": dict(sorted(bead_status.items())), "closes": closes[:60],
            "timeline": timeline[:400], "timeline_total": len(timeline),
            "claims": claims_for(claims, fillings), "evidence": evidence,
        }
        doc["summary"] = (status or {}).get("summary")
        if not doc["summary"] and fid == "kit":
            # the kit has no status file; its render set describes itself in source.json's "summary"
            try:
                doc["summary"] = json.loads((kit_dir / "source.json").read_text()).get("summary")
            except (OSError, json.JSONDecodeError):
                pass
        doc["summary_derived"] = not doc["summary"]
        if doc["summary_derived"]:
            doc["summary"] = derived_summary(verdict["counts"], len(extras.get(fid, [])), len(fixtures.get(fid, [])),
                                             len(references.get(fid, [])))
        family_docs[fid] = doc
        last = max([c["date"] for c in timeline[:1]] + [p["history"][-1]["when"][:10] for p in pairs if p["history"]] or [None], key=lambda x: x or "")
        families_index.append({
            "id": fid, "title": title, "section": section, "badge": verdict["badge"], "label": verdict["label"],
            "counts": verdict["counts"], "blocked": verdict["blocked"],
            "summary": doc["summary"], "summary_derived": doc["summary_derived"],
            "renders": sum(p["kind"] in ("graded", "listening", "refused", "compared") for p in pairs) + len(extras.get(fid, [])) + len(fixtures.get(fid, [])),
            "references": len(references.get(fid, [])), "phrases": len(phrases.get(fid, [])),
            "milestones": [{"code": m["code"], "gate": m["gate"], "status": m["status"], "closed": m["prereq"]["closed"], "total": m["prereq"]["total"]} for m in ms],
            "commits": len(timeline), "last_activity": last,
            "baseline_real_inside": (doc["baseline"] or {}).get("real", {}).get("inside_fraction") if doc["baseline"] else None,
        })

    # Program-wide process view.
    by_day = defaultdict(Counter)
    sec_of = {f[0]: f[2] for f in FAMILIES}
    fam_of_commit = defaultdict(set)
    for fid, doc in family_docs.items():
        for c in doc["timeline"]:
            fam_of_commit[c["sha"]].add(sec_of[fid])
    irp = [c for c in commits if any("instrument-realism" in t or t == "frankensim-o7qbp" for t in c["beads"])]
    for c in commits:
        secs = fam_of_commit.get(c["sha"])
        key = sorted(secs)[0] if secs and len(secs) == 1 else ("Several families" if secs else "Shared physics and tools")
        by_day[c["date"]][key] += 1
    closes_by_day = Counter((b.get("closed_at") or "")[:10] for b in beads if b.get("status") == "closed" and b.get("closed_at"))
    status_all = Counter(b["status"] for b in beads)
    process = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"), "ref": args.ref, "ref_sha": ref_sha,
        "since": PROGRAM_START, "commits_total": len(commits), "commits_irp": len(irp),
        "commits_by_day": [{"date": d, "by": dict(by_day[d])} for d in sorted(by_day)],
        "closes_by_day": [{"date": d, "n": closes_by_day[d]} for d in sorted(closes_by_day) if d >= PROGRAM_START],
        "beads": dict(sorted(status_all.items())), "beads_total": len(beads),
        "labels": dict(Counter(l for b in beads for l in b.get("labels") or [] if l in ("batch-pending", "rework", "owner-decision", "irp-milestone", "irp-data", "irp-experiment"))),
        "feed": [c for c in commits if c["type"] in ("feat", "fix", "data", "test", "perf")][:150],
        "sections": SECTIONS,
        # The grading method page: the score over time, the measures, the rulers' own check, and the pitch of every
        # note each family can play today.
        "program_score": score_history([(fid, p) for fid, doc in family_docs.items() for p in doc["pairs"]]),
        "scores": {fid: pts for fid, doc in family_docs.items() if (pts := score_history([(fid, p) for p in doc["pairs"]]))},
        "titles": {fid: doc["title"] for fid, doc in family_docs.items()},
        "measures": measure_table(family_docs),
        "calibration": [{"id": fid, "title": doc["title"], "source": doc["baseline"]["source"],
                         "recordings": doc["baseline"]["recordings"], "real": doc["baseline"]["real"]["inside_fraction"],
                         "ours": (doc["baseline"]["ours"] or {}).get("inside_fraction")}
                        for fid, doc in family_docs.items() if doc["baseline"] and doc["baseline"].get("real")],
        "pitch": [{"id": fid, "title": family_docs[fid]["title"], "sha": pitch_sha.get(fid), "notes": rows}
                  for fid, rows in pitch_table.items() if fid in family_docs],
        "pitch_method": pitch_doc.get("method"),
    }
    # The decision log: the hand-written record next to this tool, plus the owner-decision beads as the tracker has them.
    log_path = HERE.parent / "decisions.json"
    decisions = json.loads(log_path.read_text()) if log_path.is_file() else {"entries": []}
    decisions["owner_beads"] = sorted(({"id": b["id"], "title": re.sub(r"^\[IRP\]\s*", "", b["title"]), "status": b["status"],
                                        "created": (b.get("created_at") or "")[:10], "closed": (b.get("closed_at") or "")[:10] or None,
                                        "reason": (b.get("close_reason") or "")[:600]}
                                       for b in beads if "owner-decision" in (b.get("labels") or [])),
                                      key=lambda b: (b["status"] == "closed", b["created"]))
    index = {
        "built_utc": process["built_utc"], "ref": args.ref, "ref_sha": ref_sha, "mirror": mirror.name,
        "sections": SECTIONS, "families": families_index,
        "totals": {
            "families": len(FAMILIES),
            "graded_notes": sum(f["counts"]["graded"] for f in families_index),
            "listening_notes": sum(f["counts"]["listening"] for f in families_index),
            "measures_inside": sum(f["counts"]["inside"] for f in families_index),
            "measures_outside": sum(f["counts"]["outside"] for f in families_index),
            "references": sum(f["references"] for f in families_index),
            "milestones_total": sum(len(f["milestones"]) for f in families_index),
            "milestones_closed": sum(m["status"] == "closed" for f in families_index for m in f["milestones"]),
            "beads_closed": status_all.get("closed", 0), "beads_total": len(beads),
            "commits": len(commits),
        },
    }

    failures = media.run(args.jobs)
    media.save_versions()
    (out / "family").mkdir(exist_ok=True)
    for fid, doc in family_docs.items():
        (out / "family" / f"{fid}.json").write_text(json.dumps(doc, separators=(",", ":")))
    (out / "process.json").write_text(json.dumps(process, separators=(",", ":")))
    (out / "decisions.json").write_text(json.dumps(decisions, separators=(",", ":")))
    (out / "phrases.json").write_text(json.dumps({"phrases": [p for fid, *_ in FAMILIES for p in phrases.get(fid, [])],
                                                  "titles": {fid: title for fid, title, *_ in FAMILIES}}, separators=(",", ":")))
    (out / "index.json").write_text(json.dumps(index, indent=1))
    log(f"built {len(family_docs)} families, {index['totals']['graded_notes']} graded notes, "
        f"{index['totals']['references']} reference clips; {len(failures)} media failures")
    for f in failures[:20]:
        log("  FAIL", f)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
