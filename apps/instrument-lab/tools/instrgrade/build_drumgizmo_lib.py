#!/usr/bin/env python3
"""Library and manifest, in the shape decompose.py reads, for the graded DRSKit 2.1 (DrumGizmo) pieces: one overhead per piece.

usage: build_drumgizmo_lib.py REFERENCES_SCREENED.json UNZIP_ROOT WAV_DIR LIBRARY.json MANIFEST.json
REFERENCES_SCREENED.json is data/instruments/references-screened.json; UNZIP_ROOT holds DRSKit/, the unzipped
https://drumgizmo.org/kits/DRSKit/DRSKit2_1.zip (the drskit-2.1-archive entry's sha256), whose members the entry names.
Each hit is one 13-channel 44.1 kHz float WAV, one channel per microphone; the instrument XMLs give each hit's `power`.

The rows: the entry's clips are refscreen cuts of a member's two overheads (Rode NT5), the first channel the one the clip
was cut and screened on. A piece is heard by ONE overhead (MB 29109 ruling 1): the right one (OHR) for the toms and the
right crash, the left one (OHL) for the ride. A piece's clips are the VALID ones whose first source channel is its
overhead; each is cut again from its member and must give the row's clip_sha256. The grader's file is that channel of
the whole member, resampled to 48 kHz as refscreen.to_sr does and high-passed by verisim's 20-30 Hz raised cosine, as the
kit render's excerpt was (build_screened_lib's kit-vcsl rule); the library times are the member's.

A tom is a pitched decaying note at its measured principal (build_screened_lib.principal, then grader_note: the note at
which decompose's own reading holds within its 50 c); a piece is graded in its largest named-note group (MB 29109 ruling
2). DrumGizmo's editor fades every hit to zero at its last sample, one gain on all 13 channels; its length is not printed,
so it is fitted (fade_fit). A tom entry ends 0.1 s before the fitted fade starts, or at the screen's cut if that is earlier
(MB 29298 ruling 3), so neither the pair nor the ruler (which reads 0.1 s past an entry's end) reads a fade frame. A cymbal
stroke is unvoiced, its segment 5 ms before the screen's onset and 3.0 s long, the render excerpt's length, inside the cut;
its fade fit is recorded and not used. An entry's dynamic is the hit's power, so decompose.Refs.ruler reads the set's other
hits (MB 29109 ruling 2); the reference is the median-power hit, the lower middle for an even count (MB 29217 ruling 2).

The strokes follow the cards' strike positions (MB 29109 ruling 3): the crash card strikes at 0.95 of its radius, at the
edge, so the stick's shank; the ride card at 0.6, on the bow, so the tip (fs-vkplate examples/kit.rs at b08487277 and
d58291406, the commits the crash and ride renders were built from).
"""
from __future__ import annotations

import hashlib
import json
import statistics
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(1, str(Path(__file__).parent.parent))  # tools/: verisim, whose high-pass the kit excerpts carry
import build_screened_lib as bsl  # noqa: E402
import refscreen  # noqa: E402
import verisim  # noqa: E402
from instrgrade import SR, cents, hz_of, name_of  # noqa: E402

ARCHIVE_ID = "drskit-2.1-archive"
LICENCE_ID = "drskit-cc-by-4"
CC_BY = "https://creativecommons.org/licenses/by/4.0/"
EXTERNAL_CREDIT = ("DRSKit 2.1 (August 2018), a collaboration between the DrumGizmo team and Jes Eiler of DRSDrums; kit setup, "
                   "microphone setup and sampling by Deva and Lars Muldjord, edited by Lars Muldjord; CC BY 4.0. One overhead "
                   "microphone per piece, high-passed 20-30 Hz (raised cosine) as the render's excerpt was, cut and level-matched, "
                   "otherwise unchanged.")
PREFILTER_HZ = (20.0, 30.0)
PRINCIPAL_BAND_HZ = (50.0, 400.0)  # build_screened_lib's kit-vcsl band
UNVOICED_LEAD_S = 0.005
FADE_GUARD_S = 0.1  # decompose.Refs.aligned reads 0.1 s past an entry's end
FADE_GRID_S = np.round(np.arange(0.02, 1.5001, 0.01), 3)
FADE_READ_S = 2.0
C_AIR = 343.0
FILECH = {"OHL": 6, "OHR": 7}
OVERHEAD = dict(OHR="the right overhead (Rode NT5, facing the right crash and the ride)",
                OHL="the left overhead (Rode NT5, facing the left crash and the hi-hat)")
# (DrumGizmo instrument, overhead, close microphone its distance is measured behind or None, art, unvoiced segment [s] or None)
PIECES = [
    ("Tom1", "OHR", "Tom1", "DRSKit Tom1 (the hanging tom), stick", None),
    ("Tom3", "OHR", "Tom3", "DRSKit Tom3 (the lower floor tom), stick", None),
    ("Crash_right_shank", "OHR", None, "DRSKit right crash (Paiste Formula 602), stick shank at the edge", 3.0),
    ("Ride_tip", "OHL", "Ride", "DRSKit ride (Paiste Formula 602 thin crash), stick tip on the bow", 3.0),
]
LEFT_OUT = [
    dict(piece="Tom2", reason="the 14-inch tom is graded against Tom3; at Tom2's louder hits its render reads 399.9 c from the named "
                              "note, inside the grader's 400 c by 0.1 c (MB 29106, 29109)"),
    dict(piece="Snare", reason="one hit (Snare-8) passes the screen on the right overhead, and the render reads -402 c from its note, "
                               "outside the grader's 400 c (MB 29217 ruling 1)"),
    dict(piece="Kdrum_with_contact, Kdrum_without_contact", reason="no principal on the right overhead: most hits are under the grader's "
                              "floor and the rest scatter over 70-192 Hz; the close microphones are another receiver than the card's "
                              "overhead point (MB 29109 ruling 4)"),
    dict(piece="snare without wires", reason="DRSKit has no stroke with the wires off"),
]


def sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def kit_hits(unzip: Path, instr: str):
    """{hit number: (power, {microphone: 0-based file channel}, member path)} from the instrument's XML."""
    ip = unzip / "DRSKit" / instr / f"{instr}.xml"
    out = {}
    for s in ET.parse(ip).getroot().iter("sample"):
        afs = list(s.iter("audiofile"))
        out[int(s.get("name").rsplit("-", 1)[1])] = (float(s.get("power")), {a.get("channel"): int(a.get("filechannel")) - 1 for a in afs},
                                                     (ip.parent / afs[0].get("file")).resolve())
    return out


def onset_s(x, sr):
    a = np.abs(x)
    return float(np.flatnonzero(a >= 0.1 * a.max())[0]) / sr


def fade_fit(x, sr):
    """The editor's fade of one hit, fitted on all its channels (native rate, unfiltered): the log RMS of 10 ms frames counted
    back from the last non-zero sample, at distance d before it, is a_c + b_c d + log(min(1, d / L)) per channel c (its own
    level and exponential envelope), L common to the channels, a linear gain ramp ending at that sample. a_c, b_c by least
    squares for each L of FADE_GRID_S, L the least summed squared residual, on frames with d in [0.01, D], D = min(2.0, 0.7 of
    the file). Returns dict(end_s, length_s, start_s, at_grid_edge, residual_rms, residual_rms_no_fade, read_s)."""
    n_end = int(np.flatnonzero(np.any(x != 0, axis=1))[-1]) + 1
    hop = int(round(0.01 * sr))
    big_d = min(FADE_READ_S, 0.7 * n_end / sr)
    k = int(big_d * sr) // hop
    rms = np.sqrt(np.mean(x[n_end - k * hop: n_end].reshape(k, hop, x.shape[1]) ** 2, axis=1))
    d = (k - np.arange(k) - 0.5) * hop / sr
    keep = (d >= 0.01) & np.all(rms > 0, axis=1)
    y, dd = np.log(rms[keep]), d[keep]
    a = np.stack([np.ones_like(dd), dd], axis=1)

    def rss(logg):
        r = y - logg[:, None]
        coef, *_ = np.linalg.lstsq(a, r, rcond=None)
        return float(np.sum((r - a @ coef) ** 2))

    sums = np.array([rss(np.log(np.minimum(1.0, dd / length))) for length in FADE_GRID_S])
    j = int(np.argmin(sums))
    length, end = float(FADE_GRID_S[j]), n_end / sr
    return dict(end_s=round(end, 6), length_s=length, start_s=round(end - length, 6), at_grid_edge=j in (0, len(FADE_GRID_S) - 1),
                residual_rms=round(float(np.sqrt(sums[j] / y.size)), 4), residual_rms_no_fade=round(float(np.sqrt(rss(np.zeros_like(dd)) / y.size)), 4),
                read_s=round(big_d, 3))


def main():
    if len(sys.argv) != 6:
        raise SystemExit(__doc__)
    refs_p, unzip, wav_dir, lib_path, man_path = (Path(sys.argv[1]), Path(sys.argv[2]).resolve(), Path(sys.argv[3]), Path(sys.argv[4]),
                                                  Path(sys.argv[5]))
    for p in (lib_path, man_path):
        if p.exists():
            raise SystemExit(f"{p} exists: refusing to overwrite")
    refs = json.loads(refs_p.read_text())
    arc = next(f for f in refs["files"] if f["id"] == ARCHIVE_ID)
    lic = next(x for x in refs["licences"] if x["id"] == LICENCE_ID)
    (wav_dir / "clips").mkdir(parents=True, exist_ok=True)
    lib, files, left_out, pieces = [], [], list(LEFT_OUT), []
    for ins, oh, close, art, seg_s in PIECES:
        hits = kit_hits(unzip, ins)
        delays = []
        for n, (_, chans, member) in hits.items():
            if close:
                x, sr = refscreen.decode(member)
                delays.append(onset_s(x[:, chans[oh]], sr) - onset_s(x[:, chans[close]], sr))
        dist = statistics.median(delays) * C_AIR if delays else None
        where = (f"{OVERHEAD[oh]}, {dist:.2f} m from the piece by its arrival after the {close} close microphone (median of "
                 f"{len(delays)} hits, at {C_AIR:g} m/s)" if dist is not None else f"{OVERHEAD[oh]}; no close microphone on this piece, "
                 "so its distance is not measured")
        entries = []
        for m in arc["members"]:
            if m["archive_member"].split("/")[1] != ins:
                continue
            n = int(m["archive_member"].rsplit("/", 1)[1].split("-")[0])
            power, chans, member = hits[n]
            for c in m["clips"]:
                if c["status"] != "VALID" or (c.get("source_channels") or [None])[0] != FILECH[oh]:
                    continue
                x, sr = refscreen.decode(member)
                if sha256(member) != m["sha256"]:
                    raise SystemExit(f"{member}: sha256 is not the row's {m['sha256']}")
                i0, i1 = int(round(c["trim_start_s"] * sr)), int(round(c["trim_end_s"] * sr))
                cp = wav_dir / "clips" / f"{c['id']}.wav"
                if not cp.exists():
                    refscreen.write_wav_f32(cp, x[i0:i1][:, [ch - 1 for ch in c["source_channels"]]], sr)
                if sha256(cp) != c["clip_sha256"]:
                    raise SystemExit(f"{c['id']}: the clip cut again from its member is not the row's clip_sha256 {c['clip_sha256']}")
                dst = wav_dir / f"{c['id']}.wav"
                if not dst.exists():
                    refscreen.write_wav_f32(dst, verisim.highpass(refscreen.to_sr(x[:, c["source_channels"][0] - 1], sr), *PREFILTER_HZ)[:, None], SR)
                fade = fade_fit(x, sr)
                on, dyn = c["trim_start_s"] + float(c["screen"]["onset_s"]), f"power {power:.6g}"
                rec = dict(id=c["id"], power=power, dynamic=dyn, trim_s=[c["trim_start_s"], c["trim_end_s"]], clip_sha256=c["clip_sha256"],
                           fade=fade, row_cut_at=c.get("cut_at", "the row as committed in 4778d99eb2"))
                if seg_s:  # a cymbal stroke: unvoiced, no note
                    a = c["trim_start_s"] + max(0.0, float(c["screen"]["onset_s"]) - UNVOICED_LEAD_S)
                    b = a + seg_s
                    if b > c["trim_end_s"]:
                        left_out.append(dict(clip=c["id"], status="VALID", reason=f"the cut ends {c['trim_end_s'] - a:.3f} s after the "
                                                                                    f"segment starts, shorter than the {seg_s} s excerpt"))
                        continue
                    e = dict(file=dst.name, start=a, end=b, note=None, unvoiced=True, name=art, kind="unvoiced", inst="kit", art=art, dyn=dyn,
                             sul="", mallet="", clip=c["id"])
                    read_from = fade["end_s"] - fade["read_s"]
                    if b >= read_from:
                        raise SystemExit(f"{c['id']}: the {seg_s} s segment ends at {b:.3f} s, inside the last {fade['read_s']} s the fade fit reads")
                    rec.update(note=None, segment_s=[round(a, 6), round(b, 6)],
                               fade_use=(f"not used: the segment ends {read_from - b:.3f} s before the last {fade['read_s']:g} s the fit reads"
                                         + ("; fade not resolved within the last 2 s read" if fade["at_grid_edge"] else "")))
                else:
                    if fade["at_grid_edge"]:
                        raise SystemExit(f"{c['id']}: the fade fit stops at its grid's edge ({fade['length_s']} s): no end rule")
                    end = min(c["trim_end_s"], fade["start_s"]) - FADE_GUARD_S
                    f0, why = bsl.principal(x[i0:i1][:, [ch - 1 for ch in c["source_channels"]]], sr, float(c["screen"]["onset_s"]),
                                            PRINCIPAL_BAND_HZ)
                    if f0 is not None:
                        midi, f0, why = bsl.grader_note(dst, on, end, "decaying", "kit", int(round(69 + 12 * np.log2(f0 / 440.0))))
                    if f0 is None:
                        left_out.append(dict(clip=c["id"], status="VALID", reason=why))
                        continue
                    e = dict(file=dst.name, start=on, end=end, note=midi, name=name_of(midi), kind="decaying", cents_vs_et=cents(f0, hz_of(midi)),
                             inst="kit", art=art, dyn=dyn, sul="", mallet="", range=name_of(midi), f0_hz=f0, clip=c["id"])
                    rec.update(note=name_of(midi), f0_hz=round(f0, 3), entry_s=[round(on, 6), round(end, 6)],
                               end_rule=("0.1 s before the fitted fade start" if fade["start_s"] < c["trim_end_s"] else "0.1 s before the screen's cut"),
                               note_source=f"the measured principal {f0:.2f} Hz rounded to the nearest semitone")
                entries.append((e, rec, dict(id=c["id"], path=f"kit/{dst.name}", url=arc["url"], member=m["archive_member"],
                                             file_channel=c["source_channels"][0], sha256=m["sha256"], bytes=m["bytes"], page=arc["source_page"],
                                             chain=f"one 13-microphone take; this file is {where}", licence=LICENCE_ID, instrument=art,
                                             gradable=True, dynamic_source="DRSKit's power for the hit, the kit's own velocity measure")))
        if not seg_s:  # a drum is graded in its largest named-note group
            groups = {}
            for t in entries:
                groups.setdefault(t[0]["name"], []).append(t)
            sizes = sorted((len(v), k) for k, v in groups.items())
            if len(sizes) > 1 and sizes[-1][0] == sizes[-2][0]:
                raise SystemExit(f"{ins}: two largest note groups {sizes[-2:]}")
            big = sizes[-1][1]
            for k, v in groups.items():
                if k != big:
                    left_out += [dict(clip=t[1]["id"], status="VALID", reason=f"filed at {k}, not the piece's largest note group ({big} x{len(groups[big])})")
                                 for t in v]
            entries = groups[big]
        entries.sort(key=lambda t: t[1]["power"])
        ref = entries[(len(entries) - 1) // 2][1]["id"]
        for e, rec, f in entries:
            lib.append(e)
            files.append(dict(f, clips=[dict(rec, reference=rec["id"] == ref)]))
        pieces.append(dict(instrument=ins, overhead=oh, close_microphone=close, distance_m=None if dist is None else round(dist, 3),
                           hits=len(hits), graded=len(entries), note=None if seg_s else entries[0][0]["name"], art=art, reference=ref))
    lib_path.write_text(json.dumps(lib, indent=1))
    man_path.write_text(json.dumps(dict(schema="frankensim.instrument-references.v1", source=arc["source_page"],
                                        external=dict(name="DRSKit 2.1 (DrumGizmo)", url=arc["source_page"], licence=lic["name"], licence_url=CC_BY,
                                                      credit=EXTERNAL_CREDIT),
                                        licences=[dict(id=LICENCE_ID, text=lic["name"], url=CC_BY)],
                                        archive=dict(id=ARCHIVE_ID, url=arc["url"], sha256=arc["sha256"], bytes=arc["bytes"]),
                                        pieces=pieces, files=files, left_out=left_out), indent=1))
    print(json.dumps(dict(entries=len(lib), files=len(files), left_out=len(left_out))))
    for p in pieces:
        print("piece", json.dumps(p))
    for f in files:
        c = f["clips"][0]
        print(f["id"], c["dynamic"], c.get("note") or "unvoiced", c.get("entry_s") or c.get("segment_s"), c.get("end_rule") or c.get("fade_use"),
              f"fade L {c['fade']['length_s']} start {c['fade']['start_s']}{' (grid edge)' if c['fade']['at_grid_edge'] else ''}",
              "REFERENCE" if c["reference"] else "")
    for r in left_out:
        print("left out", json.dumps(r)[:300])


if __name__ == "__main__":
    main()
