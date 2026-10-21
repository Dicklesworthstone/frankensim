#!/usr/bin/env python3
"""Library and manifest for the VocalSet male long tones, in the shape decompose.py reads for the UIowa files.

usage: build_vocalset_lib.py VOCALSET_DIR WAV48_DIR LIBRARY.json MANIFEST.json [JOBS]
VOCALSET_DIR holds FULL/male<N>/long_tones/{straight,forte,pp}/m<N>_long_<technique>_<vowel>.wav
(Wilkins, Seetharaman, Wahl, Pardo 2018, Zenodo 1193957, CC BY 4.0). Each file is one singer holding one vowel on a
pitch the singer chose; the files state no pitch, so the note is the MEASURED f0 rounded to the nearest semitone and
cents_vs_et records how far from it. Each file is converted to 48 kHz mono float WAV in WAV48_DIR and becomes one
library entry whose `art` is the vowel and whose `dyn` is straight -> mf, forte -> ff, pp -> pp, so that
decompose.Refs.ruler reads "another singer on this note", "forte against straight", "a neighbouring note" from it.
Takes whose f0 is not steady for a second are left out and listed under `refused` in the manifest.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from instrgrade import SR, cents, hz_of, load, name_of, yin_track  # noqa: E402

DYN = {"straight": "mf", "forte": "ff", "pp": "pp"}
SOURCE = "https://zenodo.org/records/1193957"
LICENCE = "CC BY 4.0 (Wilkins, Seetharaman, Wahl, Pardo 2018, VocalSet: a singing voice dataset, ISMIR)"


def convert(src: Path, dst: Path) -> None:
    if dst.exists():
        return
    part = dst.with_suffix(".part.wav")
    r = subprocess.run(["ffmpeg", "-v", "error", "-n", "-i", str(src), "-ac", "1", "-ar", str(SR), "-c:a", "pcm_f32le", str(part)],
                       capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"{src.name}: {r.stderr.strip()[:200]}")
    part.rename(dst)


def analyse(job):
    src, wav48 = job
    src, wav48 = Path(src), Path(wav48)
    stem = src.stem  # m1_long_straight_a
    singer, _, technique, vowel = stem.split("_")
    dst = wav48 / f"{stem}.wav"
    convert(src, dst)
    x = load(dst)
    tr = yin_track(x, fmin=70.0, fmax=700.0)
    if not len(tr):
        return dict(stem=stem, refused="no pitch track")
    hop = 480
    lvl = np.array([20 * math.log10(np.sqrt(np.mean(x[int(t * SR):int(t * SR) + 2048] ** 2)) + 1e-9) for t in tr[:, 0]])
    voiced = np.isfinite(tr[:, 1]) & (tr[:, 2] > 0.8) & (lvl > lvl.max() - 20.0)
    if voiced.sum() * hop / SR < 1.0:
        return dict(stem=stem, refused=f"only {voiced.sum() * hop / SR:.2f} s voiced")
    f = tr[voiced, 1]
    f0 = float(np.median(f))
    spread = float(np.std(1200 * np.log2(f / f0)))
    if spread > 40.0:
        return dict(stem=stem, refused=f"f0 wanders {spread:.0f} cents (s.d.): not a steady long tone")
    midi = int(round(69 + 12 * math.log2(f0 / 440.0)))
    t_on, t_off = float(tr[voiced, 0][0]), float(tr[voiced, 0][-1]) + 2048 / SR
    return dict(stem=stem, singer=singer, technique=technique, vowel=vowel, f0_hz=f0, f0_sd_cents=spread, note=midi, name=name_of(midi),
                cents_vs_et=cents(f0, hz_of(midi)), start=t_on, end=min(len(x) / SR, t_off), src=str(src),
                sha256=hashlib.sha256(src.read_bytes()).hexdigest(), bytes=src.stat().st_size)


def main():
    if len(sys.argv) not in (5, 6):
        raise SystemExit(__doc__)
    root, wav48, lib_path, man_path = (Path(a) for a in sys.argv[1:5])
    jobs = int(sys.argv[5]) if len(sys.argv) == 6 else 8
    wav48.mkdir(parents=True, exist_ok=True)
    files = sorted(p for p in root.glob("FULL/male*/long_tones/*/*.wav") if p.parent.name in DYN)
    with ProcessPoolExecutor(min(jobs, 16)) as pool:
        res = list(pool.map(analyse, [(str(p), str(wav48)) for p in files], chunksize=2))
    lib, entries, refused = [], [], []
    for r in res:
        if "refused" in r:
            refused.append(dict(file=r["stem"], reason=r["refused"]))
            continue
        lib.append(dict(file=f'{r["stem"]}.wav', start=r["start"], end=r["end"], note=r["note"], name=r["name"], kind="sustained",
                        cents_vs_et=r["cents_vs_et"], inst="voice", art=r["vowel"], dyn=DYN[r["technique"]], sul=r["singer"], mallet="",
                        range=r["name"], f0_hz=r["f0_hz"], f0_sd_cents=r["f0_sd_cents"]))
        entries.append(dict(id=f'vocalset-{r["stem"].replace("_", "-")}', path=f'VocalSet/FULL/{r["singer"].replace("m", "male")}/long_tones/'
                                                                              f'{r["technique"]}/{r["stem"]}.wav',
                            url=SOURCE, sha256=r["sha256"], bytes=r["bytes"], kind="long-tone", instrument="Voice (male singer)",
                            articulation=r["technique"], dynamic=DYN[r["technique"]], range=r["name"], page=SOURCE,
                            chain=f'VocalSet 44.1 kHz mono, singer {r["singer"]}', licence="vocalset-cc-by-4.0", channels="mono", gradable=True,
                            note_source="measured f0 rounded to the nearest semitone; the file states no pitch"))
    lib_path.write_text(json.dumps(lib))
    man_path.write_text(json.dumps(dict(schema="frankensim.instrument-references.v1", source=SOURCE,
                                        licences=[dict(id="vocalset-cc-by-4.0", text=LICENCE, url="https://creativecommons.org/licenses/by/4.0/")],
                                        files=entries, refused=refused), indent=1))
    print(json.dumps(dict(files=len(files), entries=len(lib), refused=len(refused))))
    for r in refused:
        print("refused", r["file"], r["reason"], file=sys.stderr)


if __name__ == "__main__":
    main()
