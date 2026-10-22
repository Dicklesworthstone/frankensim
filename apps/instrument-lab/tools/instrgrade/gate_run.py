#!/usr/bin/env python3
"""gate_run: one run of the Q2 gate as the beads define it today (identity.bead_q2) over its controls and the renders, as one table.

  real controls   real recordings of the instrument from ANOTHER library, graded as if they were renders (P2's certify candidates:
                  TinySOL, CC BY 4.0; VocalSet's held-out singer m6, CC BY 4.0, against a voice library without him). Required:
                  >= 80 % pass.
  piano           the POD-round piano renders, which the owner accepted. Required: A3 mf, A3 ff and C6 mf pass.
  nulls           the five null baselines (freeze.null_signals) of every graded note's reference recording, judged against that
                  reference. Required: none passes.
  pitch shift     every real control and piano render shifted +20 cents (band-limited). Required: no margin moves by 1 dB or more.
  renders         every live render pair under ANALYSIS: reported, never required. Beside each family, Q3: the share of measures
                  inside the loose real-vs-real ruler (decompose's pair.json verdicts) for its renders and for its real controls.
  candidates      with --eval, identity_eval.py's numbers for the UNREGISTERED candidates A, B, C, measured on library recordings
                  only: printed under the table, never applied to a render here.

A note passes the current definition when d_o - d_t >= identity.BEAD_MARGIN_DB, d_t the spectral distance (spec_db) to the pair's
reference recording, d_o to the nearest note of another instrument of the same kind within 2 semitones, any dynamic or
articulation (UIowa; for the voice too). Every candidate goes the render route (decompose.render_floor, align, Note), every library
note the reference route (decompose.Refs.note). Writes OUT/gate.json (every margin) and OUT/gate.txt (the table).
It is a final validation: no number here may be used to choose, tune or fit physics.

usage (trj; numpy, pillow, matplotlib):
  gate_run.py --out DIR --piano PIANO.json [--eval DIR/identity-eval.json] [--analysis ...] [--certify CANDS.json ...] [--jobs 12]
PIANO.json: [{"note": "A3", "dyn": "mf", "wav": ..., "git_sha": ...}, ...]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import sys
import time
from collections import OrderedDict, defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import decompose as dc  # noqa: E402
import freeze  # noqa: E402
import identity as idn  # noqa: E402
import identity_eval  # noqa: E402
from instrgrade import hz_of, load, midi_of  # noqa: E402

LIBS = {
    "iowa": dict(lib="/data/tmp/stormyhill/piano/lib-light-plus-piano.json", wav="/data/tmp/stormyhill/wav48",
                 refs=dc.DEFAULTS["refs_json"]),
    "voice": dict(lib="/data/tmp/wa-inst/p2/certify/voice-lib-no-m6.json", wav="/data/tmp/wa-inst/refs-wav48-voice",
                  refs="/data/tmp/wa-inst/p2-grade/voice/references-vocalset.json"),
}
CERTIFY = ["/data/tmp/wa-inst/p2/certify/cands-tinysol-mf.json", "/data/tmp/wa-inst/p2/certify/cands-tinysol-mf-cello-viola-clarinet-flute.json",
           "/data/tmp/wa-inst/p2/certify/cands-voice-m6.json"]
OTHERS_SEMITONES = 2
REQUIRED_REAL = 0.80
PIANO_REQUIRED = (("A3", "mf"), ("A3", "ff"), ("C6", "mf"))
REFUSALS = dc.REFUSALS + (FloatingPointError, ZeroDivisionError)


def ekey(e, band=None):
    return f'{e["file"]}@{e["start"]:.3f}' + (f"|{band:.0f}" if band else "")


_R = {}


def _init(libs):
    for k, s in libs.items():
        _R[k] = dc.Refs(s["lib"], s["wav"], s["refs"])


def member_job(args):
    which, e, band = args
    refs = _R[which]
    refs.band_hz = band
    try:
        return ekey(e, band), idn.features(refs.note(e)), None
    except REFUSALS as err:
        return ekey(e, band), None, f"{type(err).__name__}: {err}"


def as_candidate(raw, f0, kind, inh, label):
    hp = dc.highpass(raw, f0)
    if not inh:
        dc.check_fundamental(hp, dc.f0_comb(hp, f0, kind, 400.0), kind, label, f0)
    y, info = dc.align(raw, f0, dc.render_floor(hp, f0))
    return dc.Note(y, f0, kind, inh, label, info, max_cents=400.0)


def cand_job(c):
    """Features of one candidate (a real control, a piano render or a render) and of its +20 c copy, or of a reference's nulls."""
    out = dict(c)
    try:
        f0 = hz_of(midi_of(c["note"]))
        inh = c["inst"] in ("xylophone", "marimba")
        if c["source"] == "null":
            refs = _R[c["which"]]
            refs.band_hz = c.get("band_hz")
            e = c["ref_entry"]
            y, _ = refs.aligned(e)
            B = refs.note(e).B if c["kind"] == "decaying" and not inh else 0.0
            out["nulls"] = {}
            for name, z in freeze.null_signals(y, f0, B).items():
                if name == "frozen" and inh:
                    out["nulls"][name] = dict(refused="a bar's partials are not k f0: there is no frozen null of a bar")
                    continue
                try:
                    out["nulls"][name] = idn.features(as_candidate(z, f0, c["kind"], inh, f"null {name}"))
                except REFUSALS as err:
                    out["nulls"][name] = dict(refused=f"{type(err).__name__}: {err}")
            return out
        raw = load(Path(c["wav"]))
        if c.get("band_hz"):
            raw = dc.lowpass(raw, c["band_hz"])
        out["sha256"] = hashlib.sha256(Path(c["wav"]).read_bytes()).hexdigest()
        out["features"] = idn.features(as_candidate(raw, f0, c["kind"], inh, c["label"]))
        if c["source"] in ("real", "piano"):
            try:
                out["shifted"] = idn.features(as_candidate(identity_eval.shift(raw, identity_eval.SHIFT_CENTS), f0, c["kind"], inh, c["label"] + " +20c"))
            except REFUSALS as err:
                out["shifted"] = dict(refused=f"{type(err).__name__}: {err}")
    except REFUSALS as err:
        out["refused"] = f"{type(err).__name__}: {err}"
    except OSError as err:  # a render file that is gone: reported in its row, never a verdict
        out["refused"] = f"unreadable: {err}"
    return out


def build_candidates(a, lib):
    cands = []

    def ref_entry(which, inst, note, dyn, art, contains):
        es = [e for e in lib[which] if e["inst"] == inst and e["note"] == midi_of(note) and e["dyn"] == dyn and e["art"] == art]
        if contains:
            es = [e for e in es if Path(e["file"]).stem.lower() == contains.lower()] or [e for e in es if contains.lower() in e["file"].lower()]
        return sorted(es, key=lambda e: e["file"])[0] if es else None

    for path in a.certify:
        for c in json.loads(Path(path).read_text()):
            which = "voice" if c["family"] == "voice" else "iowa"
            inst = dc.INST.get(c["family"], c["family"])
            e = ref_entry(which, inst, c["note"], c.get("dyn", "mf"), c["art"], c.get("ref_contains"))
            cands.append(dict(source="real", which=which, family=c["family"], inst=inst, art=c["art"], note=c["note"], slot=c.get("slot") or c["note"],
                              kind=e["kind"] if e else "sustained", wav=c["wav"], label=f'{c["family"]} {c.get("slot") or c["note"]}',
                              band_hz=c.get("band_hz"), git_sha=c.get("git_sha"), ref_entry=e))
    for p in json.loads(Path(a.piano).read_text()):
        e = ref_entry("iowa", "piano", p["note"], p["dyn"], "", None)
        cands.append(dict(source="piano", which="iowa", family="piano", inst="piano", art="", note=p["note"], slot=f'{p["note"]}-{p["dyn"]}',
                          dyn=p["dyn"], kind="decaying", wav=p["wav"], label=f'piano {p["note"]} {p["dyn"]}', band_hz=None, git_sha=p["git_sha"], ref_entry=e))
    for pj in sorted(Path(a.analysis).glob("*/*/pair.json")):
        j = json.loads(pj.read_text())
        cm = (j.get("decompose") or {}).get("candidate")
        if not j.get("ours") or not cm or not (j.get("decompose") or {}).get("library_note"):
            continue
        which = "voice" if j["family"] == "voice" else "iowa"
        inst = dc.INST.get(j["family"], j["family"])
        e = ref_entry(which, inst, cm["note"], cm.get("dyn", "mf"), cm.get("art", ""), cm.get("ref_contains"))
        cands.append(dict(source="render", which=which, family=j["family"], inst=inst, art=cm.get("art", ""), note=cm["note"], slot=pj.parent.name,
                          kind=e["kind"] if e else "sustained", wav=j["ours"]["wav"], label=f'{j["family"]} {pj.parent.name}',
                          band_hz=cm.get("band_hz"), git_sha=j["ours"].get("git_sha"), ref_entry=e))
    seen = OrderedDict()
    for c in cands:
        if c["ref_entry"] is not None:
            seen.setdefault((c["which"], ekey(c["ref_entry"], c.get("band_hz"))), c)
    nulls = [dict(source="null", which=c["which"], family=c["family"], inst=c["inst"], art=c["art"], note=c["ref_entry"]["name"], slot=c["slot"],
                  kind=c["kind"], label=f'nulls of {c["ref_entry"]["file"]}@{c["ref_entry"]["name"]}', band_hz=c.get("band_hz"),
                  ref_entry=c["ref_entry"]) for c in seen.values()]
    return cands, nulls


def others_of(c, lib):
    """The UIowa notes d_o is taken over: another instrument, the same kind, within OTHERS_SEMITONES, any dynamic or articulation."""
    m = midi_of(c["note"])
    return [e for e in lib["iowa"] if e["inst"] != c["inst"] and e["kind"] == c["kind"] and abs(e["note"] - m) <= OTHERS_SEMITONES]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--piano", required=True)
    ap.add_argument("--eval", help="identity_eval.py's identity-eval.json: the unregistered candidates, library recordings only")
    ap.add_argument("--analysis", default="/data/tmp/wa-inst/analysis")
    ap.add_argument("--certify", nargs="*", default=CERTIFY)
    ap.add_argument("--certify-analysis", default="/data/tmp/wa-inst/p2/certify/analysis", help="P2's graded certify pairs (for Q3)")
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--rows", help="an earlier run's gate.json: rebuild the table from its rows (and --eval) without recomputing")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.rows:
        prev = json.loads(Path(a.rows).read_text())
        ev = json.loads(Path(a.eval).read_text()) if a.eval else None
        table = report(prev["rows"], q3(Path(a.analysis), Path(a.certify_analysis)), ev)
        (out / "gate.json").write_text(json.dumps(dict(prev, table=table, eval=a.eval, candidates_on_library_recordings=ev, rebuilt_from=a.rows,
                                                       rebuilt_utc=datetime.now(timezone.utc).isoformat()), indent=1, default=float))
        (out / "gate.txt").write_text(table + "\n")
        print(table)
        return
    t_start = time.time()
    lib = {w: json.loads(Path(s["lib"]).read_text()) for w, s in LIBS.items()}
    cands, nulls = build_candidates(a, lib)
    want = OrderedDict()  # every reference and every other-instrument note some candidate needs
    for c in cands + nulls:
        if c["ref_entry"] is not None:
            want.setdefault((c["which"], ekey(c["ref_entry"], c.get("band_hz"))), (c["which"], c["ref_entry"], c.get("band_hz")))
        for e in others_of(c, lib):
            want.setdefault(("iowa", ekey(e, c.get("band_hz"))), ("iowa", e, c.get("band_hz")))
    cache_p = out / "features-cache.json"
    cache = json.loads(cache_p.read_text()) if cache_p.exists() else {}
    cache = {k: v for k, v in cache.items() if (v.get("features") or {}).get("version") == idn.IDENTITY_VERSION or "refused" in v}
    # by band and file: a worker's Refs keeps the last few low-passed files, and a band-limited range file costs seconds to re-read
    todo = sorted((v for k, v in want.items() if f"{k[0]}:{k[1]}" not in cache), key=lambda v: (v[0], v[2] or 0.0, v[1]["file"], v[1]["start"]))
    print(f"candidates {len(cands)}, null sources {len(nulls)}, library notes {len(want)} ({len(todo)} to compute)", flush=True)
    with ProcessPoolExecutor(min(a.jobs, 16), mp_context=multiprocessing.get_context("fork"), initializer=_init, initargs=(LIBS,)) as pool:
        for (k, f, err), (which, _, _) in zip(pool.map(member_job, todo, chunksize=8), todo):
            cache[f"{which}:{k}"] = dict(features=f) if f else dict(refused=err)
        cache_p.write_text(json.dumps(cache))
        done = list(pool.map(cand_job, cands + nulls, chunksize=1))
    feats = {k: v["features"] for k, v in cache.items() if v.get("features")}
    rows = []
    for c in done:
        r = {k: v for k, v in c.items() if k not in ("ref_entry", "features", "shifted", "nulls")}
        ref = feats.get(f"{c['which']}:{ekey(c['ref_entry'], c.get('band_hz'))}") if c.get("ref_entry") else None
        r["reference"] = ekey(c["ref_entry"]) if c.get("ref_entry") else None
        others = [(ekey(e), feats[k]) for e in others_of(c, lib) if (k := f"iowa:{ekey(e, c.get('band_hz'))}") in feats]
        if ref is None and "refused" not in r:
            r["refused"] = "the reference recording has no features (refused or missing)"
        if c["source"] == "null":
            r["nulls"] = {n: (dict(recognisable=None, refused=f["refused"]) if "refused" in f else idn.bead_q2(f, ref, others))
                          for n, f in c["nulls"].items()} if ref is not None else {}
        elif "features" in c and ref is not None:
            r["q2"] = idn.bead_q2(c["features"], ref, others)
            if "shifted" in c and "refused" not in c["shifted"]:
                qs = idn.bead_q2(c["shifted"], ref, others)
                r["shift"] = dict(q2=qs, margin_change_db=abs(qs["margin_db"] - r["q2"]["margin_db"]))
        rows.append(r)
    ev = json.loads(Path(a.eval).read_text()) if a.eval else None
    table = report(rows, q3(Path(a.analysis), Path(a.certify_analysis)), ev)
    meta = dict(schema="frankensim.q2-gate-run.v1", definition="bead (identity.bead_q2)", margin_db=idn.BEAD_MARGIN_DB,
                identity_version=idn.IDENTITY_VERSION, utc=datetime.now(timezone.utc).isoformat(), elapsed_s=round(time.time() - t_start, 1),
                tools_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(Path(__file__).resolve().parent.glob("*.py"))},
                libraries=LIBS, analysis=a.analysis, certify=a.certify, piano=a.piano, eval=a.eval)
    (out / "gate.json").write_text(json.dumps(dict(meta, rows=rows, table=table, candidates_on_library_recordings=ev), indent=1, default=float))
    (out / "gate.txt").write_text(table + "\n")
    print(table)


def q3(analysis, certify):
    """{family: {renders: [inside, measures], real: [inside, measures]}} from the loose verdicts of the graded pairs."""
    out = defaultdict(lambda: dict(renders=[0, 0], real=[0, 0]))
    for root, k in ((analysis, "renders"), (certify, "real")):
        for pj in sorted(root.glob("*/*/pair.json")):
            v = (json.loads(pj.read_text()).get("verdicts") or {}).get("loose")
            if v:
                out[pj.parent.parent.name][k][0] += v["inside"]
                out[pj.parent.parent.name][k][1] += v["inside"] + v["outside"]
    return out


def rate(vs):
    n = [x for x in vs if x is not None]
    return (sum(n), len(n), len(vs) - len(n))


def fmt(t):
    k, n, none = t
    return (f"{k:3d}/{n:<3d} {100 * k / n:5.1f} %" if n else "  -/-        ") + (f" ({none} refused)" if none else "")


def passes(r):
    return (r.get("q2") or {}).get("recognisable")


def report(rows, q3s, ev):
    L = [f"Q2 gate run, the CURRENT definition (identity.bead_q2: spec_db margin to the nearest other instrument >= {idn.BEAD_MARGIN_DB:g} dB)", ""]
    L.append(f'{"group":44s} {"passes":>24s}   {"required":20s} verdict')
    real = [r for r in rows if r["source"] == "real"]
    k, n, none = rate([passes(r) for r in real])
    ok_real = n > 0 and k / n >= REQUIRED_REAL and none == 0
    L.append(f'{"real recordings, another library (all)":44s} {fmt((k, n, none)):>24s}   {">= 80 %, all graded":20s} {"PASS" if ok_real else "FAIL"}')
    for f in sorted({r["family"] for r in real}):
        rf = [r for r in real if r["family"] == f]
        med = np.median([r["q2"]["margin_db"] for r in rf if r.get("q2")]) if any(r.get("q2") for r in rf) else float("nan")
        L.append(f'{"  " + f:44s} {fmt(rate([passes(r) for r in rf])):>24s}   median margin {med:+.1f} dB')
    ok_piano = True
    for r in [r for r in rows if r["source"] == "piano"]:
        req = (r["note"], r.get("dyn")) in PIANO_REQUIRED
        got = passes(r)
        ok_piano &= (got is True) or not req
        mg = (r.get("q2") or {}).get("margin_db")
        txt = (f"{got} ({mg:+.1f} dB)" if mg is not None else str(r.get("refused", got)))[:24]
        L.append(f'{"piano render " + r["slot"]:44s} {txt:>24s}   {"pass" if req else "reported":20s} {("PASS" if got else "FAIL") if req else ""}')
    ok_null = True
    for name in freeze.NULLS:
        vs = [r["nulls"][name].get("recognisable") for r in rows if r["source"] == "null" and name in r.get("nulls", {})]
        k, n, none = rate(vs)
        ok_null &= k == 0
        L.append(f'{"null " + name:44s} {fmt((k, n, none)):>24s}   {"none":20s} {"PASS" if k == 0 else "FAIL"}')
    sh = [r["shift"]["margin_change_db"] for r in rows if r.get("shift")]
    flips = sum(1 for r in rows if r.get("shift") and r["shift"]["q2"]["recognisable"] != r["q2"]["recognisable"])
    ok_shift = bool(sh) and max(sh) < 1.0
    L.append(f'{"+20 c shift: largest margin change":44s} {(f"{max(sh):.3f} dB" if sh else "-"):>24s}   {"< 1 dB":20s} {"PASS" if ok_shift else "FAIL"}'
             f"   ({len(sh)} notes, median {np.median(sh) if sh else float('nan'):.3f} dB, {flips} verdict flips)")
    rend = [r for r in rows if r["source"] == "render"]
    for f in sorted({r["family"] for r in rend}):
        q = q3s.get(f, {})
        q3txt = "  Q3 inside: " + ", ".join(f"{k} {q[k][0]}/{q[k][1]} ({100 * q[k][0] / q[k][1]:.0f} %)" for k in ("renders", "real") if q[k][1]) if q else ""
        realf = [passes(x) for x in real if x["family"] == f]
        L.append(f'{"renders " + f:44s} {fmt(rate([passes(r) for r in rend if r["family"] == f])):>24s}   reported'
                 + (f"  (real controls {fmt(rate(realf)).strip()})" if realf else "") + q3txt)
    L += ["", f'GATE {"VALID" if (ok_real and ok_piano and ok_null and ok_shift) else "NOT VALID"}: real controls {ok_real}, piano {ok_piano}, '
               f"nulls {ok_null}, pitch shift {ok_shift}"]
    if ev:
        L += ["", "UNREGISTERED candidates, library recordings only (identity_eval.py; never applied to a render here):",
              identity_eval.report(ev)]
    return "\n".join(L)


if __name__ == "__main__":
    main()
