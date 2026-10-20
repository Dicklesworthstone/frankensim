#!/usr/bin/env python3
"""Convert the fetched UIowa references to the grader's 48 kHz mono float WAVs.

usage: convert_refs.py REFERENCES.json REFDIR OUTDIR [JOBS]
REFDIR is the FRANKENSIM_REFERENCE_DIR that data/instruments/fetch-references.sh filled.
Only entries with gradable=true are converted (the .stereo partners of chromatic scales are not),
plus entries whose measured defect is a sample-rate header: those are re-read at the true rate.
Stereo-only sessions keep channel 0 (references.json channel_policy). Existing outputs are kept.
The WAV is named <file stem>.wav, which is what instrgrade.parse_ref_name reads.
"""
from __future__ import annotations

import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

SR = 48000
TRUE_RATE_HZ = 96000  # viola sulC/sulG 2012 mono files: 96 kHz samples under a 44.1 kHz header


def command(entry, src: Path, dst: Path) -> list[str]:
    filters = []
    if any(d["kind"] == "sample-rate-header" for d in entry.get("defects", [])):
        filters.append(f"asetrate={TRUE_RATE_HZ}")
    filters.append(f"aresample={SR}")
    if entry["channels"] == "stereo":
        filters.insert(0, "pan=mono|c0=c0")
    return ["ffmpeg", "-v", "error", "-n", "-i", str(src), "-af", ",".join(filters), "-ac", "1", "-c:a", "pcm_f32le", str(dst)]


def main():
    if len(sys.argv) not in (4, 5):
        raise SystemExit(__doc__)
    manifest, refdir, outdir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    jobs = int(sys.argv[4]) if len(sys.argv) == 5 else 8
    outdir.mkdir(parents=True, exist_ok=True)
    todo, missing = [], 0
    for e in json.loads(manifest.read_text())["files"]:
        defect_rate = any(d["kind"] == "sample-rate-header" for d in e.get("defects", []))
        if e["kind"] == "ambient-silence" or not (e["gradable"] or defect_rate):
            continue
        src = refdir.joinpath(*e["path"].split("/"))
        dst = outdir / (Path(e["path"]).stem.replace(" ", "") + ".wav")
        if not src.is_file():
            missing += 1
        elif not dst.exists():
            todo.append((e, src, dst))

    def run(job):
        e, src, dst = job
        part = dst.with_suffix(".part.wav")
        r = subprocess.run(command(e, src, part), capture_output=True, text=True)
        if r.returncode:
            return f"FAILED {src.name}: {r.stderr.strip()[:200]}"
        part.rename(dst)
        return None

    with ThreadPoolExecutor(jobs) as pool:
        errors = [m for m in pool.map(run, todo) if m]
    print(json.dumps(dict(converted=len(todo) - len(errors), failed=len(errors), source_missing=missing, out=str(outdir))))
    for m in errors:
        print(m, file=sys.stderr)
    raise SystemExit(1 if errors else 0)


if __name__ == "__main__":
    main()
