#!/usr/bin/env python3
"""Library and manifest, in the shape decompose.py reads, for one recorded instrument's screened clips.

usage: build_screened_lib.py SCREENED.json SOURCE_ROOT SET WAV_DIR LIBRARY.json MANIFEST.json
SCREENED.json is a refscreen manifest (data/instruments/references-screened.json, or a hunt's fragment: a list of its
file rows). SET names one recorded instrument in SETS below: its family, which file rows are that instrument, how an
entry's note, dynamic and `sul` are read, and the attribution its licence asks for, the manifest's `external` block
(decompose refuses a non-Iowa licence without one). decompose.Refs.ruler reads other takes, other dynamics and
neighbouring notes only within one instrument, so one SET is one instrument (one pack, one set of drums), never a mix.

A library file is the whole SOURCE, not the clip: decompose reads a file's noise floor as the 5th percentile of its
20 ms frames, and a clip cut 0.1 s before a held note's onset has too little quiet for that (a pan pipe's 8 s clip read
its floor inside the note and was refused, IRPDECOMP_NO_SIGNAL). Each source is found under SOURCE_ROOT by its URL's
file name and its row's sha256, decoded by refscreen.decode (ffmpeg, its own rate and channels) and written whole to
WAV_DIR/<row id>.wav. Each clip is rebuilt from the decoded source between its trim times the way refscreen cut wrote
it (WAV_DIR/clips/<clip id>.wav) and must reproduce the row's clip_sha256, so the entry's samples are the screened ones.
The grader reads the first channel and resamples it to 48 kHz.

Only VALID clips become library entries. A clip whose only failure is pitch, re-screened VALID at the nearest note
(`rescreen_as_nearest`), is filed at that note: the source named the pitch class and the instrument sounds the next one.
A set whose rows name no pitch (a struck drum: `note_from` "principal") files each clip at its measured principal: the
strongest line inside the set's principal band over instrgrade describe's sustain window gives a first note, and the
entry is filed at the note where decompose's own reading of the library file (Refs.note: align, then Note, the strongest
line within +-7 semitones of the note) lies within the grader's 50 c, moving to the note the grader measures until it
holds. A clip with no such note is left out with the reason. An entry's `start` is the screen's onset and its `end` the
clip's end, both in the source's time.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path
from urllib.parse import unquote, urlparse

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(1, str(Path(__file__).parent.parent))  # tools/: verisim, whose high-pass the kit excerpts carry
import decompose as dc  # noqa: E402
import refscreen  # noqa: E402
import verisim  # noqa: E402
from instrgrade import (SR, cents, describe, estimate_noise_floor_db, highpass, hz_of, load, midi_of, name_of, spectrum,  # noqa: E402
                        strongest_line, windows)

LOSSY = "The reference is a lossy MP3 preview (about 130-180 kbit/s), the file Freesound serves openly"
VCSL = dict(url="https://github.com/sgossner/VCSL", licence="CC0 1.0", licence_url="https://creativecommons.org/publicdomain/zero/1.0/")
SETS = {
    "pan-pipes-timohanes": dict(
        family="pan-pipes", inst="pan-pipes", rows=lambda f: f.get("author") == "timohanes", note_from="named", dyn_from="mf", sul_from="",
        external=dict(name="Freesound pack 'Pan_pipe' by timohanes", url="https://freesound.org/people/timohanes/packs/4159/",
                      licence="CC0 1.0", licence_url="http://creativecommons.org/publicdomain/zero/1.0/",
                      credit=f"timohanes, Freesound pack 4159 'Pan_pipe' (2008), CC0 1.0; maker, material and tuning not stated. {LOSSY}; "
                             "cut to one note and resampled, otherwise unchanged.")),
    "tin-whistle-waltons": dict(
        family="tin-whistle", inst="tin-whistle", rows=lambda f: f.get("author") == "cernoel", note_from="named", dyn_from="mf", sul_from="",
        external=dict(name="Freesound pack 'Waltons Irish ThinWhistle in D' by cernoel", url="https://freesound.org/people/cernoel/packs/22189/",
                      licence="CC0 1.0", licence_url="http://creativecommons.org/publicdomain/zero/1.0/",
                      credit=f"cernoel, Freesound pack 22189 'Waltons Irish ThinWhistle in D' (2017, recorded with an SM57), CC0 1.0. {LOSSY}; "
                             "the short note of each file, cut to one note and resampled, otherwise unchanged.")),
    # The timpani card tunes its head to VCSL Timpani 2B's principal, so 2B is never its reference (MB 27665, Q1 ruling C):
    # the drum is graded only against kettles it was not tuned to.
    "timpani-vcsl2-not-2b": dict(
        family="timpani", inst="timpani", rows=lambda f: f["id"].startswith("vcsl2-timpani") and not f["id"].startswith("vcsl2-timpani2b_"),
        excluded=dict(prefix="vcsl2-timpani2b_", reason="the timpani card tunes its head to this kettle's principal (MB 27665, Q1 ruling C)"),
        note_from="principal", principal_band_hz=(50.0, 400.0), dyn_from="clip", sul_from="piece",
        external=dict(name="Versilian Community Sample Library (VCSL), Timpani 2", **VCSL,
                      credit="Versilian Studios, Versilian Community Sample Library (VCSL), Membranophones / Timpani 2 hits, CC0 1.0 "
                             "(github.com/sgossner/VCSL at c1ea7bcc3c73). Cut to one strike and resampled, otherwise unchanged; VCSL names "
                             "each kettle and velocity layer, not its pitch, so the note is the strike's measured principal.")),
    # The kit page grades kitgrade.py's excerpts of its renders (the raw renders are not copied off their host, MB 27684). Those
    # excerpts were high-passed by verisim's 20-30 Hz raised cosine, so the same filter is applied to each VCSL source before the
    # grader cuts it, and the low band is compared like for like. Each VCSL drum and stroke is its own instrument (`art`): one drum
    # is never another drum's other take. A drum is a pitched decaying note at its measured principal; a cymbal or hi-hat stroke is
    # unvoiced (pair --unvoiced), its segment 5 ms before the screen's onset for as long as the excerpt of the render piece it is
    # paired with, since an unvoiced pair compares the render's last span of the segment's length. VCSL's concert bass drums are
    # not a kit kick, and the kick stays ungraded (MB 27665). The kit's crash and ride are struck with a stick's tip; of VCSL's
    # two orchestral suspended cymbals only the strokes named hit_stick are stick strokes, one VALID take of each.
    "kit-vcsl": dict(
        family="kit", inst="kit", rows=lambda f: f["id"].startswith("vcsl-") and "bass drum" not in (f.get("instrument") or ""),
        note_from="principal", principal_band_hz=(50.0, 400.0), dyn_from="clip", sul_from="", art_from="instrument+articulation",
        prefilter_hz=(20.0, 30.0), unvoiced_s={"hi-hat cymbal, HitC": 0.6, "hi-hat cymbal, HitO": 2.5,
                                              "suspended cymbal 1 (orchestral), hit_stick": 3.0,
                                              "suspended cymbal 2 (orchestral), hit_stick": 3.0}, unvoiced_lead_s=0.005,
        pitched_arts=("snare drum", "tom "),
        external=dict(name="Versilian Community Sample Library (VCSL), drums and cymbals", **VCSL,
                      credit="Versilian Studios, Versilian Community Sample Library (VCSL), Membranophones and Idiophones (snare drums, toms, "
                             "hi-hat, suspended cymbals), CC0 1.0 (github.com/sgossner/VCSL at c1ea7bcc3c73). High-passed 20-30 Hz (raised "
                             "cosine) as the render's excerpt was, cut to one stroke and resampled, otherwise unchanged.")),
}


def principal(clip, sr, onset_s, band):
    """(f_hz, why_not): the clip's principal as decompose reads an inharmonic decaying note, or None and the reason."""
    y = refscreen.to_sr(clip[:, 0], sr)
    y = highpass(y, None)
    y = y[max(0, int((onset_s - 0.01) * SR)):]
    t0, t1 = windows(y, "decaying", estimate_noise_floor_db(y), sr=SR, allow_unended_decay=True)["sustain"]
    f, s = spectrum(y, t0, t1)
    wide = strongest_line(f, s, *band)
    n = int(round(69 + 12 * math.log2(wide / 440.0)))
    d = describe(y, "decaying", hz_of(n), inharmonic=True, _skip_q3_heavy=True)
    if d is None or abs(1200 * math.log2(float(d["f0"]) / wide)) > 1.0:
        return None, f"no clear principal: the strongest line in {band[0]:g}-{band[1]:g} Hz is {wide:.2f} Hz, the strongest within 7 semitones of {name_of(n)} is {d and d['f0']}"
    return wide, None


def grader_note(path, start, end, kind, inst, n0):
    """(MIDI note, f0) at which decompose itself reads library file `path` [start, end] within its 50 c check, starting from n0 and
    moving to the note the grader measures until it holds: the alignment and Note of decompose.Refs.note, read from the very file
    the grader reads. None when no note holds (the grader's reading moves with the note it is told), with the reason."""
    x = load(path)
    seen, n = set(), n0
    while n not in seen:
        seen.add(n)
        f0 = hz_of(n)
        a = max(0, int((start - 2.0) * SR))
        b = min(len(x), int(end * SR) + int(0.1 * SR))
        try:
            y, info = dc.align(x[a:b], f0, dc.floor_of_file(x, f0), anchor_s=start - a / SR, at_start=(a == 0))
            nt = dc.Note(y, f0, kind, inst in dc.INHARMONIC, path.name, info, max_cents=1200.0)
        except (dc.OnsetMisaligned, dc.GraderRefusalError, ValueError) as e:
            return None, None, f"decompose refuses it at {name_of(n)}: {e}"
        off = 1200 * math.log2(nt.f0 / f0)
        if abs(off) <= 50.0:
            return n, nt.f0, None
        n = int(round(n + off / 100.0))
    return None, None, f"decompose's reading moves with the note it is told (tried {', '.join(name_of(m) for m in sorted(seen))})"


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    if len(sys.argv) != 7:
        raise SystemExit(__doc__)
    screened, src_root, name = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
    wav_dir, lib_path, man_path = Path(sys.argv[4]), Path(sys.argv[5]), Path(sys.argv[6])
    if name not in SETS:
        raise SystemExit(f"unknown set {name!r}: one of {sorted(SETS)}")
    for p in (lib_path, man_path):
        if p.exists():
            raise SystemExit(f"{p} exists: refusing to overwrite")
    st = SETS[name]
    man = json.loads(screened.read_text())
    all_rows = [f for f in (man["files"] if isinstance(man, dict) else man) if f.get("family") == st["family"]]
    rows = [f for f in all_rows if st["rows"](f)]
    (wav_dir / "clips").mkdir(parents=True, exist_ok=True)
    lib, files, left_out = [], [], []
    ex = st.get("excluded")
    if ex:
        left_out += [dict(row=f["id"], reason=ex["reason"]) for f in all_rows if f["id"].startswith(ex["prefix"])]
    for f in rows:
        chosen = []
        for c in f["clips"]:
            note, scr, status = c.get("note"), c.get("screen") or {}, c["status"]
            rs = c.get("rescreen_as_nearest")
            if status != "VALID" and c.get("failed") == ["pitch"] and rs and rs.get("status") == "VALID":
                note, scr, status = rs["note"], rs, "VALID"
            if status != "VALID":
                left_out.append(dict(clip=c["id"], status=status, failed=c.get("failed")))
            else:
                chosen.append((c, note, scr))
        if not chosen:
            continue
        base = unquote(Path(urlparse(f["url"]).path).name)
        found = [p for p in src_root.rglob(base) if p.is_file() and sha256(p) == f["sha256"]]
        if not found:
            raise SystemExit(f"{f['id']}: no file named {base} with sha256 {f['sha256']} under {src_root}")
        x, sr = refscreen.decode(found[0])
        dst = wav_dir / f"{f['id']}.wav"
        if not dst.exists():
            if st.get("prefilter_hz"):  # the render's excerpt filter, on the grader's own channel and rate (instrgrade.load reads channel 0 at 48 kHz)
                y = verisim.highpass(refscreen.to_sr(x[:, 0], sr), *st["prefilter_hz"])
                refscreen.write_wav_f32(dst, y[:, None], SR)
            else:
                refscreen.write_wav_f32(dst, x, sr)
        kept = []
        for c, note, scr in chosen:
            i0, i1 = int(round(c["trim_start_s"] * sr)), int(round(c["trim_end_s"] * sr))
            clip = x[i0:i1]
            cp = wav_dir / "clips" / f"{c['id']}.wav"
            if not cp.exists():
                refscreen.write_wav_f32(cp, clip, sr)
            if sha256(cp) != c["clip_sha256"]:
                raise SystemExit(f"{c['id']}: the clip rebuilt from the source is not the row's clip_sha256 {c['clip_sha256']}")
            f0 = scr.get("f0_comb_hz")
            art = f"{f.get('instrument')}, {c.get('articulation')}" if st.get("art_from") == "instrument+articulation" else ""
            dyn = c["dynamic"] if st["dyn_from"] == "clip" else "mf"
            if art in st.get("unvoiced_s", {}):  # a cymbal or hi-hat stroke: unvoiced, no note
                a = c["trim_start_s"] + max(0.0, float(scr["onset_s"]) - st["unvoiced_lead_s"])
                b = a + st["unvoiced_s"][art]
                if b > c["trim_end_s"]:
                    left_out.append(dict(clip=c["id"], status="VALID", reason=f"the clip ends {c['trim_end_s'] - a:.3f} s after the segment "
                                                                                f"starts, shorter than the {st['unvoiced_s'][art]} s excerpt it pairs with"))
                    continue
                lib.append(dict(file=dst.name, start=a, end=b, note=None, unvoiced=True, name=art, kind="unvoiced", inst=st["inst"], art=art,
                                dyn=dyn, sul="", mallet="", clip=c["id"]))
                kept.append(dict(id=c["id"], note=None, dynamic=dyn, trim_s=[c["trim_start_s"], c["trim_end_s"]], clip_sha256=c["clip_sha256"],
                                 segment_s=[round(a, 6), round(b, 6)], articulation=art,
                                 segment_source=f"{st['unvoiced_lead_s'] * 1000:g} ms before the screen's onset, as long as the render excerpt"))
                continue
            if st.get("pitched_arts") and not str(f.get("instrument")).startswith(st["pitched_arts"]):
                left_out.append(dict(clip=c["id"], status="VALID", reason=f"{art}: no render piece of this family is paired with this stroke"))
                continue
            if st["note_from"] == "principal":
                f0, why = principal(clip, sr, float(scr["onset_s"]), st["principal_band_hz"])
                if f0 is not None:
                    midi, f0, why = grader_note(dst, c["trim_start_s"] + float(scr["onset_s"]), c["trim_end_s"], c["kind"], st["inst"],
                                                int(round(69 + 12 * math.log2(f0 / 440.0))))
                if f0 is None:
                    left_out.append(dict(clip=c["id"], status="VALID", reason=why))
                    continue
            else:
                midi = midi_of(note)
            sul = (c.get("piece") or "") if st["sul_from"] == "piece" else ""
            lib.append(dict(file=dst.name, start=c["trim_start_s"] + float(scr["onset_s"]), end=c["trim_end_s"], note=midi, name=name_of(midi),
                            kind=c["kind"], cents_vs_et=cents(f0, hz_of(midi)) if f0 else None, inst=st["inst"], art=art, dyn=dyn, sul=sul,
                            mallet="", range=name_of(midi), f0_hz=f0, clip=c["id"]))
            kept.append(dict(id=c["id"], note=name_of(midi), dynamic=dyn, piece=c.get("piece"), trim_s=[c["trim_start_s"], c["trim_end_s"]],
                             clip_sha256=c["clip_sha256"], note_source=c.get("pitch_from") or (
                                 f"the measured principal {f0:.2f} Hz rounded to the nearest semitone" if st["note_from"] == "principal" else ""),
                             articulation=c.get("articulation", "")))
        if kept:
            files.append(dict(id=f["id"], path=f"{st['family']}/{dst.name}", url=f["url"], sha256=f["sha256"], bytes=f.get("bytes"),
                              page=f.get("source_page") or st["external"]["url"], chain=f.get("notes") or f.get("format") or f.get("instrument") or "",
                              licence=f["licence"], clips=kept, gradable=True,
                              dynamic_source=("the source's own velocity layer" if st["dyn_from"] == "clip" else
                                              "the source states no dynamic; filed as mf, the dynamic a render is paired at")))
    if not lib:
        raise SystemExit(f"set {name}: no VALID clip")
    lib_path.write_text(json.dumps(lib, indent=1))
    licences = sorted({(f["licence"], f.get("licence_name", ""), f.get("licence_url", "")) for f in rows})
    man_path.write_text(json.dumps(dict(schema="frankensim.instrument-references.v1", source=st["external"]["url"], external=st["external"],
                                        licences=[dict(id=i, text=t, url=u) for i, t, u in licences], files=files, left_out=left_out), indent=1))
    print(json.dumps(dict(set=name, rows=len(rows), files=len(files), entries=len(lib), left_out=len(left_out))))
    for e in lib:
        print(e["file"], e["clip"], e["name"], e["dyn"], e["art"] or e["sul"], f'{e["f0_hz"]:.2f} Hz' if e.get("f0_hz") else "-",
              f'{e["cents_vs_et"]:+.1f} c' if e.get("cents_vs_et") is not None else "", f'start {e["start"]:.4f} end {e["end"]:.4f}')
    for r in left_out:
        if r.get("reason"):
            print("left out", r)


if __name__ == "__main__":
    main()
