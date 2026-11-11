#!/usr/bin/env python3
"""Library and manifest for the MusicSphere Ilmenau organ pipe tones (Part 3), in the shape decompose.py reads.

usage: build_musicsphere_lib.py MEMBERS_DIR WAV48_DIR LIBRARY.json MANIFEST.json KEY [KEY ...]
MEMBERS_DIR holds members of the Zenodo 19497225 zips at their paths in the zip, "Tone <T>/Mic <m> Register 1 Tone <T>.wav"
and "Tone <T>/MIDI <T>.mid" (Stolz, Treybig, McLeod, Klein, Werner 2026, MusicSphere Part 3, CC BY 4.0). Register 1 is the
Principal 16'. KEY is a MIDI key name such as B5. Each WAV is one microphone hearing the register while an Orgamat plays the
MIDI file: every key of tone T alone and in chords, held 4 s with 4 s of silence after, the whole sequence three times.
The members are range-read, never the whole zip (MB 27665): https://zenodo.org/api/records/19497225/files/Tone%20B.zip/content
(4281976773 bytes, Zenodo md5 aad934f7f334b5db6a1674d9b52f6639) and .../Tone%20C.zip/content (8762222921 bytes, md5
9daa6c312d7d295089ea77a41b51230d); the manifest pins each WAV member's sha256 and byte count.

Only a key played ALONE becomes a library entry. Each of its presses on each microphone is one take, cut to its own 48 kHz
mono float WAV in WAV48_DIR from 2.5 s before the key (or the file's first sample, if the key comes sooner) to 3.5 s after
its release (the next key is 4 s after the release), so that decompose.Refs.ruler reads another press as another take
and another microphone's file as another take too.
The files name the key, not the pitch the pipe sounds: the note is the MEASURED f0 of the held tone rounded to the nearest
semitone and cents_vs_et records how far from it (the organ sounds about 39 c under A440, and a 16' pipe an octave under
its key). The WAV's clock is placed against the MIDI clock by the first onset in the file; a press whose own onset lies
more than 0.25 s from where the MIDI puts it, or whose held f0 wanders more than 10 c (s.d.), is refused and listed under
`refused` in the manifest. A library entry's `start` is 0.3 s after the measured onset (decompose walks back from it to
the onset) and its `end` is the release.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import struct
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from instrgrade import SR, cents, hz_of, load, midi_of, name_of, write_wav, yin_track  # noqa: E402

SOURCE = "https://doi.org/10.5281/zenodo.19497225"
LICENCE = ("CC BY 4.0 (Stolz, Treybig, McLeod, Klein, Werner 2026, MusicSphere: Room and Organ Acoustic Measurements of "
           "St. Jakobus Church Ilmenau, Part 3: Recordings of Individual Pipe Tones)")
LEAD_S, TAIL_S, ONSET_TOL_S, WANDER_C = 2.5, 3.5, 0.25, 10.0
MEMBER = re.compile(r"Mic (\d+) Register 1 Tone (\w+)\.wav")


def midi_notes(path: Path):
    """[(on_s, off_s, key)] of a Standard MIDI File, the tempo map applied; a key held in a chord is marked by `alone`."""
    b = path.read_bytes()
    if b[:4] != b"MThd":
        raise ValueError(f"{path} is not a Standard MIDI File")
    ntrk, div = struct.unpack(">HH", b[10:14])
    i, ev, tempos = 14, [], [(0, 500000)]

    def vlq(j):
        v = 0
        while True:
            c = b[j]
            j += 1
            v = (v << 7) | (c & 0x7F)
            if not c & 0x80:
                return v, j
    for _ in range(ntrk):
        n = struct.unpack(">I", b[i + 4:i + 8])[0]
        j, end, tick, status = i + 8, i + 8 + n, 0, 0
        while j < end:
            d, j = vlq(j)
            tick += d
            c = b[j]
            if c == 0xFF:
                ln, k = vlq(j + 2)
                if b[j + 1] == 0x51:
                    tempos.append((tick, int.from_bytes(b[k:k + 3], "big")))
                j = k + ln
                continue
            if c in (0xF0, 0xF7):
                ln, k = vlq(j + 1)
                j = k + ln
                continue
            if c & 0x80:
                status, j = c, j + 1
            hi = status & 0xF0
            if hi in (0x80, 0x90, 0xA0, 0xB0, 0xE0):
                if hi == 0x90 and b[j + 1] > 0:
                    ev.append((tick, 1, b[j]))
                elif hi in (0x80, 0x90):
                    ev.append((tick, 0, b[j]))
                j += 2
            else:
                j += 1
        i = end
    tempos.sort()

    def secs(t):
        s, last, us = 0.0, 0, 500000
        for tt, tu in tempos:
            if tt > t:
                break
            s, last, us = s + (tt - last) * us / 1e6 / div, tt, tu
        return s + (t - last) * us / 1e6 / div
    held, notes = {}, []
    for tick, on, key in sorted(ev):
        if on:
            held[key] = secs(tick)
        elif key in held:
            notes.append((held.pop(key), secs(tick), key))
    starts = {}
    for on, _, _ in notes:
        starts[on] = starts.get(on, 0) + 1
    return [dict(on=on, off=off, key=key, alone=starts[on] == 1) for on, off, key in sorted(notes)]


def level_db(x, hop):
    n = len(x) // hop
    return 10 * np.log10(np.mean(x[: n * hop].reshape(n, hop) ** 2, axis=1) + 1e-18)


def first_onset_s(x, hop=SR // 100):
    """The first 10 ms frame 30 dB over the file's quietest second."""
    lv = level_db(x, hop)
    quiet = float(np.min(np.convolve(lv, np.ones(100) / 100, mode="valid")))
    return float(np.flatnonzero(lv > quiet + 30.0)[0] * hop / SR), quiet


def onset_near_s(x, t, quiet, hop=SR // 1000):
    """The onset nearest MIDI time t: walking back from t + 0.3 s to the last 1 ms frame under quiet + 30 dB."""
    a = max(0, int((t - 1.0) * SR))
    lv = level_db(x[a: int((t + 0.3) * SR)], hop)
    under = np.flatnonzero(lv <= quiet + 30.0)
    return a / SR + (float(under[-1] + 1) * hop / SR if len(under) else 0.0)


def convert(src: Path, dst: Path) -> None:
    part = dst.with_suffix(".part.wav")
    r = subprocess.run(["ffmpeg", "-v", "error", "-n", "-i", str(src), "-ac", "1", "-ar", str(SR), "-c:a", "pcm_f32le", str(part)],
                       capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"{src.name}: {r.stderr.strip()[:200]}")
    part.rename(dst)


def main():
    if len(sys.argv) < 6:
        raise SystemExit(__doc__)
    members, wav48, lib_path, man_path = (Path(a) for a in sys.argv[1:5])
    keys = [midi_of(k) for k in sys.argv[5:]]
    for p in (lib_path, man_path):
        if p.exists():
            raise SystemExit(f"{p} exists: refusing to overwrite")
    wav48.mkdir(parents=True, exist_ok=True)
    lib, files, refused = [], [], []
    for src in sorted(members.glob("Tone */Mic * Register 1 Tone *.wav")):
        m = MEMBER.fullmatch(src.name)
        mic, tone = m.group(1), m.group(2)
        mid = src.parent / f"MIDI {tone}.mid"
        notes = midi_notes(mid)
        whole = wav48 / f"musicsphere-p3-tone{tone}-mic{mic}-reg1-whole.wav"
        if not whole.exists():
            convert(src, whole)
        x = load(whole)
        t0, quiet = first_onset_s(x)
        offset = notes[0]["on"] - t0  # MIDI time - WAV time
        sha = hashlib.sha256(src.read_bytes()).hexdigest()
        press = {}
        for nt in notes:
            if not nt["alone"] or nt["key"] not in keys:
                continue
            press[nt["key"]] = press.get(nt["key"], 0) + 1
            take = press[nt["key"]]
            kname = name_of(nt["key"])
            on_w, off_w = nt["on"] - offset, nt["off"] - offset
            got = onset_near_s(x, on_w, quiet)
            tag = f"tone{tone} mic {mic} key {kname} take {take}"
            if abs(got - on_w) > ONSET_TOL_S:
                refused.append(dict(file=src.name, take=take, key=kname, reason=f"onset at {got:.3f} s, the MIDI puts it at {on_w:.3f} s"))
                continue
            lead = min(LEAD_S, got)  # the first press of a file begins a fraction of a second into it: its cut starts at the file's first sample
            a, b = int(round((got - lead) * SR)), int(round((off_w + TAIL_S) * SR))
            if b > len(x):
                refused.append(dict(file=src.name, take=take, key=kname, reason=f"{tag}: the cut {a / SR:.2f}-{b / SR:.2f} s leaves the file"))
                continue
            seg = x[a:b]
            tr = yin_track(seg[int((lead + 1.0) * SR): int((off_w - got + lead - 0.2) * SR)], fmin=30.0, fmax=2000.0)
            ok = np.isfinite(tr[:, 1]) & (tr[:, 2] > 0.8) if len(tr) else np.zeros(0, bool)
            if not ok.any():
                refused.append(dict(file=src.name, take=take, key=kname, reason=f"{tag}: no pitch track on the hold"))
                continue
            f = tr[ok, 1]
            f0 = float(np.median(f))
            sd = float(np.std(1200 * np.log2(f / f0)))
            if sd > WANDER_C:
                refused.append(dict(file=src.name, take=take, key=kname, reason=f"{tag}: the hold wanders {sd:.1f} c (s.d.)"))
                continue
            midi = int(round(69 + 12 * math.log2(f0 / 440.0)))
            stem = f"musicsphere-p3-mic{mic}-reg1-key{kname}-take{take}"
            dst = wav48 / f"{stem}.wav"
            if dst.exists():
                raise SystemExit(f"{dst} exists: refusing to overwrite")
            write_wav(dst, seg)
            lib.append(dict(file=dst.name, start=lead + 0.3, end=off_w - got + lead, note=midi, name=name_of(midi), kind="sustained",
                            cents_vs_et=cents(f0, hz_of(midi)), inst="organ", art="", dyn="mf", sul=f"Mic {mic}", mallet="",
                            range=name_of(midi), f0_hz=f0, f0_sd_cents=sd, key=kname))
            files.append(dict(id=stem, path=f"MusicSphere-Part3/{stem}.wav", url=SOURCE,
                              member=f"Tone {tone}/{src.name}", sha256=sha, bytes=src.stat().st_size, kind="pipe-tone",
                              instrument="Organ, St. Jakobus Ilmenau, Principal 16'", articulation="held", dynamic="mf",
                              range=name_of(midi), page=SOURCE, chain=f"MusicSphere Part 3, 96 kHz 24-bit mono, Mic {mic}",
                              licence="musicsphere-cc-by-4.0", channels="mono", gradable=True, key=kname, take=take,
                              cut_s=[round(a / SR, 6), round(b / SR, 6)], midi_minus_wav_s=round(offset, 6),
                              onset_minus_midi_s=round(got - on_w, 4),
                              note_source=f"key {kname} on the Principal 16'; the note is the measured f0 {f0:.2f} Hz rounded to the nearest semitone"))
    external = dict(name="MusicSphere: Room and Organ Acoustic Measurements of St. Jakobus Church Ilmenau - Part 3", url=SOURCE,
                    licence="CC BY 4.0", licence_url="https://creativecommons.org/licenses/by/4.0/",
                    credit=("Stolz, Treybig, McLeod, Klein, Werner, MusicSphere: Room and Organ Acoustic Measurements of St. Jakobus Church "
                            "Ilmenau, Part 3: Recordings of Individual Pipe Tones, Zenodo, 2026, doi 10.5281/zenodo.19497225. Cut from the "
                            "longer take and resampled; otherwise unchanged."))
    lib_path.write_text(json.dumps(lib, indent=1))
    man_path.write_text(json.dumps(dict(schema="frankensim.instrument-references.v1", source=SOURCE, external=external,
                                        licences=[dict(id="musicsphere-cc-by-4.0", text=LICENCE, url="https://creativecommons.org/licenses/by/4.0/")],
                                        files=files, refused=refused), indent=1))
    print(json.dumps(dict(notes=len(lib), refused=len(refused))))
    for e in lib:
        print(e["file"], e["name"], f'{e["f0_hz"]:.3f} Hz', f'{e["cents_vs_et"]:+.1f} c', f'sd {e["f0_sd_cents"]:.2f} c', f'end {e["end"]:.3f}')
    for r in refused:
        print("refused", r)


if __name__ == "__main__":
    main()
