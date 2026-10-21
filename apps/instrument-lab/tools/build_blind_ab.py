"""Blind A/B pairs for the review page: two renders that differ in ONE thing, and the real reference, ready for the owner to judge blind.

The first-principles analysis decides which physics matters by one-change-at-a-time render pairs that the owner judges without knowing
which is which. A lane (P9 bars, P10 flute, P7 guitar, ...) hands over one folder per pair:

    <pairs>/<pair-id>/pair.json      (+ the three audio files it names, any format ffmpeg reads)

    {"schema": "frankensim.blind-ab-pair.v1",
     "id": "<pair-id>",                      the folder's name; lower case, digits, hyphens
     "set": "P9 bars",                       the lane's set; groups the pairs
     "title": "Xylophone C5, pair 1",        SHOWN BEFORE THE ANSWER: must not say what differs or which clip is which
     "question": "Which is closer to the real xylophone C5?",   SHOWN BEFORE THE ANSWER, same rule
     "change": "silicone vs boxwood mallet, everything else identical",   shown only AFTER the answer is recorded
     "made_by": "P9",
     "level_match": "rms",                   "rms" (default): A, B and the reference are played at one RMS; "none": as rendered
     "a": {"label": "silicone mallet", "file": "a.wav", "commit": "<git sha, 7 to 40 hex>", "recipe": "<command or script>"},
     "b": {"label": "boxwood mallet",  "file": "b.wav", "commit": "...", "recipe": "..."},
     "reference": {"label": "real xylophone C5", "file": "real.wav", "source": "<where the recording is from>", "licence": "<optional>"}}

tools/build_blind_ab.py --pairs <pairs> --out <dir> [--placeholder-from <analysis-dir>]

Writes <out>/index.json and, per pair, three AAC .mp4 clips under opaque names (a hash of the encoded file, never a or b), so nothing in a
file name says which clip is which. The page shows the clips as X and Y in a random order the browser draws, plays the reference on its own
button, and reveals the labels, commits and `change` only after the owner's answer is recorded. A pair that fails a check is listed in
index.json (and on the page) and the build goes on: exit code 2, as for the candidates.

--placeholder-from builds a placeholder pair from two renders already in an analysis tree (clarinet A4 and Bb4, against the real A4), flagged
placeholder in the data and on the page, so the section can be tried before the first real pair arrives. It is not a test of anything.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_candidates import MAX_PLAY_S, PEAK_LIMIT_DBFS, TARGET_RMS_DBFS, Refused, decode, encode, text  # noqa: E402

SCHEMA = "frankensim.blind-ab-pair.v1"
INDEX_SCHEMA = "frankensim.blind-ab-index.v1"
SLUG = re.compile(r"^[a-z0-9][a-z0-9-]*$")
COMMIT = re.compile(r"^[0-9a-f]{7,40}$")
MIN_S = 0.2
PLACEHOLDER_ID = "placeholder-clarinet-a4-vs-bb4"


def rel_path(v) -> bool:
    return text(v) and not Path(v).is_absolute() and ".." not in Path(v).parts


def validate(rec, folder: Path) -> list[str]:
    """Everything about a pair.json that can be wrong without decoding any audio."""
    if not isinstance(rec, dict):
        return ["pair.json is not an object"]
    p: list[str] = []
    if rec.get("schema") != SCHEMA:
        p.append(f"schema must be {SCHEMA}")
    pid = rec.get("id")
    if not (isinstance(pid, str) and SLUG.match(pid)):
        p.append("id must be lower case letters, digits and hyphens")
    elif pid != folder.name:
        p.append(f"id {pid!r} is not the folder's name {folder.name!r}")
    for k in ("set", "title", "question", "change", "made_by"):
        if not text(rec.get(k)):
            p.append(f"{k} must be a non-empty string")
    if rec.get("level_match", "rms") not in ("rms", "none"):
        p.append("level_match must be 'rms' or 'none'")
    if not isinstance(rec.get("placeholder", False), bool):
        p.append("placeholder must be true or false")
    labels = []
    for slot in ("a", "b"):
        c = rec.get(slot)
        if not isinstance(c, dict):
            p.append(f"{slot} must be an object")
            continue
        for k in ("label", "recipe"):
            if not text(c.get(k)):
                p.append(f"{slot}.{k} must be a non-empty string")
        if not rel_path(c.get("file")):
            p.append(f"{slot}.file must be a relative path inside the pair's folder")
        elif not (folder / c["file"]).is_file():
            p.append(f"{slot}.file {c['file']} is missing")
        if not (isinstance(c.get("commit"), str) and COMMIT.match(c["commit"])):
            p.append(f"{slot}.commit must be a git sha of 7 to 40 hex digits")
        if text(c.get("label")):
            labels.append(c["label"].strip().lower())
    if len(labels) == 2 and labels[0] == labels[1]:
        p.append("a.label and b.label are the same: the owner could not be told afterwards which was which")
    shown = f"{rec.get('title', '')} {rec.get('question', '')}".lower()
    for lab in labels:
        if lab and lab in shown:
            p.append(f"title or question contains the label {lab!r}: both are shown BEFORE the answer and would give it away")
    for slot in ("a", "b"):
        sha = (rec.get(slot) or {}).get("commit") if isinstance(rec.get(slot), dict) else None
        if isinstance(sha, str) and sha and sha.lower() in shown:
            p.append(f"title or question contains the commit of {slot}")
    r = rec.get("reference")
    if not isinstance(r, dict):
        p.append("reference must be an object")
    else:
        for k in ("label", "source"):
            if not text(r.get(k)):
                p.append(f"reference.{k} must be a non-empty string")
        if "licence" in r and not text(r["licence"]):
            p.append("reference.licence, when given, must be a non-empty string")
        if not rel_path(r.get("file")):
            p.append("reference.file must be a relative path inside the pair's folder")
        elif not (folder / r["file"]).is_file():
            p.append(f"reference.file {r['file']} is missing")
    return p


def rms_db(x: np.ndarray) -> float:
    return 20 * math.log10(math.sqrt(float(np.mean(x ** 2))))


def peak_db(x: np.ndarray) -> float:
    return 20 * math.log10(float(np.max(np.abs(x))))


def load(path: Path) -> dict:
    x, sr = decode(path)
    if len(x) / sr < MIN_S:
        raise Refused([f"{path.name} lasts {len(x) / sr:.2f} s, under {MIN_S} s"])
    if float(np.max(np.abs(x))) <= 0:
        raise Refused([f"{path.name} is silent"])
    return {"path": path, "sha": hashlib.sha256(path.read_bytes()).hexdigest(), "seconds": len(x) / sr, "rms": rms_db(x), "peak": peak_db(x), "rate": sr}


def build_pair(rec: dict, folder: Path, out: Path) -> dict:
    """Decode, level-match, encode. Raises Refused for a clip that cannot be used."""
    clips = {"a": load(folder / rec["a"]["file"]), "b": load(folder / rec["b"]["file"]), "reference": load(folder / rec["reference"]["file"])}
    if clips["a"]["sha"] == clips["b"]["sha"]:
        raise Refused(["a.file and b.file are the same bytes: there is nothing to tell apart"])
    match = rec.get("level_match", "rms")
    # One common target for all three, so loudness is not a cue: the page-wide target, lowered only as far as the peakiest clip needs.
    target = min([TARGET_RMS_DBFS] + [c["rms"] + (PEAK_LIMIT_DBFS - c["peak"]) for c in clips.values()]) if match == "rms" else None
    entries = {}
    folder_out = out / rec["id"]
    folder_out.mkdir(parents=True, exist_ok=True)
    for slot, c in clips.items():
        gain = round(target - c["rms"], 2) if target is not None else 0.0
        seconds = c["seconds"] if c["seconds"] <= MAX_PLAY_S else MAX_PLAY_S
        tmp = folder_out / f".enc-{slot}.mp4"
        encode(c["path"], tmp, gain, None if c["seconds"] <= MAX_PLAY_S else MAX_PLAY_S)
        name = f"{hashlib.sha256(tmp.read_bytes()).hexdigest()[:10]}.mp4"
        tmp.replace(folder_out / name)
        entries[slot] = {"audio": f"{rec['id']}/{name}", "gain_db": gain, "rms_dbfs": round(c["rms"], 2), "peak_dbfs": round(c["peak"], 2),
                         "source_duration_s": round(c["seconds"], 3), "playback_s": round(seconds, 3)}
    return {
        "id": rec["id"], "set": rec["set"], "title": rec["title"], "question": rec["question"], "change": rec["change"], "made_by": rec["made_by"],
        "placeholder": bool(rec.get("placeholder", False)), "level_match": match, "target_rms_dbfs": None if target is None else round(target, 2),
        "a": {**entries["a"], "label": rec["a"]["label"], "commit": rec["a"]["commit"], "recipe": rec["a"]["recipe"]},
        "b": {**entries["b"], "label": rec["b"]["label"], "commit": rec["b"]["commit"], "recipe": rec["b"]["recipe"]},
        "reference": {**entries["reference"], "label": rec["reference"]["label"], "source": rec["reference"]["source"], **({"licence": rec["reference"]["licence"]} if "licence" in rec["reference"] else {})},
    }


def placeholder_pair(analysis: Path, into: Path) -> Path | None:
    """A pair folder made of two existing renders (clarinet A4 and Bb4, against the real A4), or None when the analysis tree lacks them."""
    try:
        a = json.loads((analysis / "clarinet/A4/pair.json").read_text())
        b = json.loads((analysis / "clarinet/Bb4/pair.json").read_text())
        files = [analysis / "clarinet/A4/ours.wav", analysis / "clarinet/Bb4/ours.wav", analysis / "clarinet/A4/real.wav"]
        if not all(f.is_file() for f in files) or not a.get("ours") or not b.get("ours"):
            return None
    except (OSError, ValueError):
        return None
    folder = into / PLACEHOLDER_ID
    folder.mkdir(parents=True)
    for f, name in zip(files, ("a.wav", "b.wav", "real.wav")):
        shutil.copyfile(f, folder / name)
    ours = lambda p: p["ours"]  # noqa: E731
    rec = {
        "schema": SCHEMA, "id": PLACEHOLDER_ID, "set": "Page trial", "placeholder": True, "made_by": "P8",
        "title": "Placeholder pair: clarinet", "question": "Which of these two is closer to the real clarinet recording? (Placeholder: answer anything.)",
        "change": "NOT A REAL TEST. One clip is our clarinet A4 render and the other our clarinet Bb4 render (a different note): two existing renders used only so the page can be tried before the first real pair arrives. The recording is the real A4, so the A4 render is the expected answer.",
        "a": {"label": "our clarinet A4 render", "file": "a.wav", "commit": ours(a)["git_sha"], "recipe": f"{ours(a)['recipe']} (observer {ours(a).get('observer', 'unstated')})"},
        "b": {"label": "our clarinet Bb4 render", "file": "b.wav", "commit": ours(b)["git_sha"], "recipe": f"{ours(b)['recipe']} (observer {ours(b).get('observer', 'unstated')})"},
        "reference": {"label": "real clarinet A4", "file": "real.wav", "source": "University of Iowa Musical Instrument Samples, the recording the grader uses for the A4 pair"},
    }
    (folder / "pair.json").write_text(json.dumps(rec, indent=1))
    return folder


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--pairs", type=Path, help="a folder of <pair-id>/pair.json folders")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--placeholder-from", type=Path, help="an analysis tree: add a placeholder pair made of two of its existing clarinet renders")
    args = ap.parse_args()
    if args.pairs is None and args.placeholder_from is None:
        raise SystemExit("name --pairs and/or --placeholder-from")
    if args.pairs is not None and not args.pairs.is_dir():
        raise SystemExit(f"no pairs directory {args.pairs}")
    folders = sorted(p for p in args.pairs.iterdir() if p.is_dir() and not p.name.startswith(("_", "."))) if args.pairs else []
    pairs, refused, seen = [], [], set()
    with tempfile.TemporaryDirectory() as tmp:
        if args.placeholder_from is not None:
            ph = placeholder_pair(args.placeholder_from, Path(tmp))
            if ph is None:
                print(f"WARNING no placeholder: {args.placeholder_from} lacks clarinet A4/Bb4 renders", file=sys.stderr)
            else:
                folders.append(ph)
        for folder in folders:
            f = folder / "pair.json"
            try:
                rec = json.loads(f.read_text())
            except (OSError, ValueError) as e:
                refused.append({"id": folder.name, "problems": [f"pair.json cannot be read: {e}"]})
                continue
            problems = validate(rec, folder)
            if not problems and rec["id"] in seen:
                problems.append("id appears twice")
            if not problems:
                try:
                    pairs.append(build_pair(rec, folder, args.out))
                    seen.add(rec["id"])
                except Refused as e:
                    problems = e.problems
            if problems:
                refused.append({"id": folder.name, "problems": problems})
                print(f"REFUSED {folder.name}: {'; '.join(problems)}", file=sys.stderr)
    args.out.mkdir(parents=True, exist_ok=True)
    pairs.sort(key=lambda p: (p["placeholder"], p["set"], p["id"]))
    (args.out / "index.json").write_text(json.dumps({"schema": INDEX_SCHEMA, "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "pairs": pairs, "refused": refused}, ensure_ascii=False, separators=(",", ":")))
    real = sum(1 for p in pairs if not p["placeholder"])
    print(f"{len(pairs)} blind A/B pairs ({real} real, {len(pairs) - real} placeholder), {len(refused)} refused -> {args.out}")
    return 2 if refused else 0


if __name__ == "__main__":
    raise SystemExit(main())
