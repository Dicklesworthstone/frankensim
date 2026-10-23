#!/usr/bin/env python3
"""gate_run: one run of Q2 Alive (identity.q2_alive, registered on pexfp.357) over its controls and the renders, as one table.

  real controls   real recordings of the instrument from ANOTHER library, graded as if they were renders (P2's certify candidates:
                  TinySOL, CC BY 4.0; VocalSet's held-out singer m6, CC BY 4.0, against a voice library without him). Required:
                  >= 80 % pass, every one graded.
  nulls           the five null baselines (freeze.null_signals) of every sustained graded note's reference recording, judged
                  against that note's set. Required: none passes.
  pitch shift     every real control shifted +20 cents (band-limited). Required: no Q2 Alive margin moves by 1 dB or more.
  piano           the POD-round piano renders: struck, so Q2 Alive makes no claim (owner, 2026-10-22); listed, never required.
  renders         every live render pair under ANALYSIS: reported, never required; struck and plucked families get no claim.
                  Beside each family, Q3: the share of measures inside the loose real-vs-real ruler (decompose's pair.json
                  verdicts) for its renders and for its real controls.
  beside          never gating: identity (identity.nearest_set: is the note's nearest real recording on the tilt-corrected
                  cepstral envelope of its own instrument?), confusion (does a real control pass Q2 Alive against ANOTHER
                  instrument's set?), and the retired definition identity.bead_q2, so the 2026-10-22 table can be reproduced.
  candidates      with --eval, identity_eval.py's numbers for the candidates Q2 Alive was chosen from, on library recordings only.

Every candidate goes the render route (decompose.render_floor, align, Note), every library note the reference route
(decompose.Refs.note). With --chain-refdir a render first goes through its reference session's recording chain (chain.py).
Writes OUT/gate.json (every margin) and OUT/gate.txt (the table). It is a final validation: no number here may be used to
choose, tune or fit physics.

usage (trj; numpy, pillow, matplotlib):
  gate_run.py --out DIR --piano PIANO.json [--eval DIR/identity-eval.json] [--analysis ...] [--certify CANDS.json ...] [--jobs 12]
PIANO.json: [{"note": "A3", "dyn": "mf", "wav": ..., "git_sha": ...}, ...]
"""
from __future__ import annotations

import argparse
import hashlib
import math
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
# Sustained families whose real recordings fail Q2 Alive: the family's realism claim is NOT valid yet, so its renders are not graded
# (treated like struck/plucked). clarinet: its steady reed tone has genuinely low post-attack flux, near the frozen null; whether that
# is a reference-set gap or a real limit of post-attack flux is pexfp.360 / the clarinet-calibration follow-up (verifier, 2026-10-22).
BLOCKED_FAMILIES = {"clarinet"}
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


_CHAINS, _MANIFESTS = {}, {}


def chain_for(c):
    """(the recording chain of the session c's reference comes from, None) or (None, why), built once per worker and session."""
    import chain as chn

    manifest = Path(LIBS[c["which"]]["refs"])
    m = _MANIFESTS.setdefault(str(manifest), json.loads(manifest.read_text()))
    stem = Path(c["ref_entry"]["file"]).stem
    f = next((f for f in m.get("files", []) if Path(f["path"]).stem.replace(" ", "") == stem), None)
    if f is None or "pages" not in m:
        return None, f"IRPCHAIN_UNKNOWN_SESSION: {manifest.name} has no session record for {stem}"
    if f["page"] not in _CHAINS:
        try:
            _CHAINS[f["page"]] = (chn.Chain.from_references(manifest, f["page"], Path(c["chain_refdir"]), _R[c["which"]].lib), None)
        except chn.ChainError as err:
            _CHAINS[f["page"]] = (None, str(err))
    return _CHAINS[f["page"]]


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
        note = as_candidate(raw, f0, c["kind"], inh, c["label"])
        if c.get("chain_refdir") and c["source"] == "render":
            ch, why = chain_for(c)
            if ch is None:
                out["chain"] = dict(refused=why)
            else:
                refs = _R[c["which"]]
                refs.band_hz = c.get("band_hz")
                gain = 10 ** ((refs.note(c["ref_entry"]).level_db - note.level_db) / 20)  # decompose's level match, on the raw render
                y, rc = ch.apply(raw * gain, f'{c["label"]} {out["sha256"][:12]}')
                out["chain"] = dict(rc, level_match_db=20 * math.log10(gain))
                note = as_candidate(y, f0, c["kind"], inh, c["label"] + " through its session's chain")
        out["features"] = idn.features(note)
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
                  ref_entry=c["ref_entry"]) for c in seen.values() if c["kind"] == "sustained"]  # Q2 Alive judges sustained notes only
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
    ap.add_argument("--chain-refdir", help="put every render through its reference session's recording chain (chain.py), reading the "
                                           "session's files from this FRANKENSIM_REFERENCE_DIR; there is no default chain")
    ap.add_argument("--sources", default="real,piano,render,null", help="which candidates to grade (comma-separated)")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.rows:
        prev = json.loads(Path(a.rows).read_text())
        if prev.get("statistic") != idn.STATISTIC:
            raise SystemExit(f"IRPGATE_STATISTIC: {a.rows} was graded with {prev.get('statistic') or prev.get('definition')!r}, not {idn.STATISTIC}")
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
    keep = set(a.sources.split(","))
    cands = [dict(c, chain_refdir=a.chain_refdir) for c in cands if c["source"] in keep]
    nulls = nulls if "null" in keep else []
    want = OrderedDict()  # every reference and every other-instrument note some candidate needs
    for c in cands + nulls:
        if c["ref_entry"] is not None:
            want.setdefault((c["which"], ekey(c["ref_entry"], c.get("band_hz"))), (c["which"], c["ref_entry"], c.get("band_hz")))
        for e in others_of(c, lib):
            want.setdefault(("iowa", ekey(e, c.get("band_hz"))), ("iowa", e, c.get("band_hz")))
        for e in lib[c["which"]]:  # the candidate's own set, which Q2 Alive judges it against
            if idn.in_set(e, c["inst"], c["art"], c["kind"], midi_of(c["note"])):
                want.setdefault((c["which"], ekey(e, c.get("band_hz"))), (c["which"], e, c.get("band_hz")))
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
    canon = {}  # each instrument's canonical articulation (decompose's choice): the set another instrument's note is judged against
    for e in lib["iowa"]:
        canon.setdefault(e["inst"], set()).add(e["art"])
    canon = {inst: dc.canon_art(arts) for inst, arts in canon.items()}
    rulers = {}

    def alive(f, S, kind):
        k = (tuple(la for la, _ in S), kind)
        if k not in rulers and kind == "sustained" and len(S) >= idn.MIN_SET:
            rulers[k] = {m: idn.ruler(S, m) for m in idn.LIFE_MEASURES}
        return idn.q2_alive(f, S, kind, rulers.get(k))

    rows = []
    for c in done:
        r = {k: v for k, v in c.items() if k not in ("ref_entry", "features", "shifted", "nulls")}
        band = c.get("band_hz")
        ref = feats.get(f"{c['which']}:{ekey(c['ref_entry'], band)}") if c.get("ref_entry") else None
        r["reference"] = ekey(c["ref_entry"]) if c.get("ref_entry") else None
        r["mean_flux_db"] = (c.get("features") or {}).get("mean_flux_db")
        m = midi_of(c["note"])
        S = [(ekey(e), feats[k]) for e in lib[c["which"]] if idn.in_set(e, c["inst"], c["art"], c["kind"], m)
             and (k := f"{c['which']}:{ekey(e, band)}") in feats]
        others = [(ekey(e), feats[k]) for e in others_of(c, lib) if (k := f"iowa:{ekey(e, band)}") in feats]
        by_inst = defaultdict(list)  # other instruments' sets at this pitch, canonical articulation
        for e in others_of(c, lib):
            if e["art"] == canon.get(e["inst"]) and (k := f"iowa:{ekey(e, band)}") in feats:
                by_inst[e["inst"]].append((ekey(e), feats[k]))
        r["set_size"] = len(S)
        if c["source"] == "null":
            r["nulls"] = {n: (dict(recognisable=None, refused=f["refused"]) if "refused" in f else alive(f, S, c["kind"]))
                          for n, f in c["nulls"].items()}
        elif "features" in c:
            r["q2"] = alive(c["features"], S, c["kind"])
            r["identity"] = idn.nearest_set(c["features"], S, dict(by_inst)) if S and by_inst else None
            r["retired_bead_q2"] = idn.bead_q2(c["features"], ref, others) if ref is not None else None
            if "shifted" in c and "refused" not in c["shifted"] and r["q2"]["margin_db"] is not None:
                qs = alive(c["shifted"], S, c["kind"])
                r["shift"] = dict(q2=qs, margin_change_db=abs(qs["margin_db"] - r["q2"]["margin_db"]) if qs["margin_db"] is not None else None)
            if c["source"] == "real":
                r["confusion"] = {inst: alive(c["features"], So, c["kind"])["recognisable"] for inst, So in by_inst.items() if len(So) >= idn.MIN_SET}
        rows.append(r)
    ev = json.loads(Path(a.eval).read_text()) if a.eval else None
    table = report(rows, q3(Path(a.analysis), Path(a.certify_analysis)), ev)
    meta = dict(schema="frankensim.q2-gate-run.v2", statistic=idn.STATISTIC, definition="Q2 Alive (identity.q2_alive, pexfp.357)",
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
    return (f"{k:3d}/{n:<3d} {100 * k / n:5.1f} %" if n else "  -/-        ") + (f" ({none} no claim)" if none else "")


def passes(r):
    return (r.get("q2") or {}).get("recognisable")


def ident(rs):
    """Identity beside Q2 Alive: how many have their own instrument's recording nearest (identity.nearest_set)."""
    return fmt(rate([(r.get("identity") or {}).get("own_nearest") for r in rs]))


def retired(rs):
    return fmt(rate([(r.get("retired_bead_q2") or {}).get("recognisable") for r in rs]))


def report(rows, q3s, ev):
    L = [f"Q2 Alive gate run ({idn.STATISTIC}: mean spectral flux after the attack, as a level, inside the spread of the instrument's "
         f"real recordings within {idn.SEMITONES} semitones; sustained families only)", ""]
    L.append(f'{"group":44s} {"Q2 Alive passes":>24s}   {"required":20s} {"verdict":8s} {"identity (beside)":>22s} {"retired bead_q2":>22s}')
    real = [r for r in rows if r["source"] == "real"]
    k, n, none = rate([passes(r) for r in real])
    ok_real = n > 0 and k / n >= REQUIRED_REAL  # refusals ("could not grade") are not failures: n is the graded count, refusals enumerated below (pexfp.357, option a)
    L.append(f'{"real recordings, another library (all)":44s} {fmt((k, n, none)):>24s}   {">= 80 % of graded":20s} {"PASS" if ok_real else "FAIL":8s} '
             f"{ident(real):>22s} {retired(real):>22s}")
    for f in sorted({r["family"] for r in real}):
        rf = [r for r in real if r["family"] == f]
        blocked = "   BLOCKED: family Q2 Alive not valid" if f in BLOCKED_FAMILIES else ""
        L.append(f'{"  " + f:44s} {fmt(rate([passes(r) for r in rf])):>24s}   {"":20s} {"":8s} {ident(rf):>22s} {retired(rf):>22s}{blocked}')
    refusals = [r for r in real if passes(r) is None]
    if refusals:
        L.append(f'refusals (could not grade; not failures, excluded from the {n} graded above): {len(refusals)}')
        for r in sorted(refusals, key=lambda r: r.get("label", "")):
            why = r.get("refused") or (r.get("q2") or {}).get("reason") or "no claim"
            L.append(f'  {r.get("label") or r.get("note") or "?":40s} {why}')
    nl = [r for r in rows if r["source"] == "null"]
    ok_null = bool(nl)  # nulls not graded in this run: not a pass
    for name in freeze.NULLS:
        vs = [r["nulls"][name].get("recognisable") for r in nl if name in r.get("nulls", {})]
        k, n, none = rate(vs)
        ok_null &= k == 0
        L.append(f'{"null " + name + " (sustained references)":44s} {fmt((k, n, none)):>24s}   {"none":20s} '
                 f'{("PASS" if k == 0 else "FAIL") if n else "not graded"}')
    sh = [r["shift"]["margin_change_db"] for r in real if r.get("shift") and r["shift"]["margin_change_db"] is not None]
    flips = sum(1 for r in real if r.get("shift") and r["shift"]["q2"]["recognisable"] != r["q2"]["recognisable"])
    ok_shift = bool(sh) and max(sh) < 1.0
    L.append(f'{"+20 c shift: largest Q2 Alive margin change":44s} {(f"{max(sh):.3f} dB" if sh else "-"):>24s}   {"< 1 dB":20s} '
             f'{"PASS" if ok_shift else "FAIL":8s} ({len(sh)} real controls, median {np.median(sh) if sh else float("nan"):.3f} dB, {flips} verdict flips)')
    for r in [r for r in rows if r["source"] == "piano"]:
        L.append(f'{"piano render " + r["slot"]:44s} {"no claim (struck)":>24s}   {"owner 2026-10-22":20s}')
    rend = [r for r in rows if r["source"] == "render"]
    for f in sorted({r["family"] for r in rend}):
        rr = [r for r in rend if r["family"] == f]
        q = q3s.get(f, {})
        q3txt = "  Q3 inside: " + ", ".join(f"{k} {q[k][0]}/{q[k][1]} ({100 * q[k][0] / q[k][1]:.0f} %)" for k in ("renders", "real") if q[k][1]) if q else ""
        realf = [passes(x) for x in real if x["family"] == f]
        if f in BLOCKED_FAMILIES:
            verdict = "no claim (family not valid)"
        elif any(r["kind"] == "sustained" for r in rr):
            verdict = fmt(rate([passes(r) for r in rr]))
        else:
            verdict = "no claim (struck/plucked)"
        L.append(f'{"renders " + f:44s} {verdict:>24s}   {"reported":20s} {"":8s} {ident(rr):>22s} {retired(rr):>22s}'
                 + (f"  real controls {fmt(rate(realf)).strip()}" if realf else "") + q3txt)
    conf = [x for r in real for x in (r.get("confusion") or {}).values()]
    L.append(f'{"confusion: real control passes as ANOTHER":44s} {fmt(rate(conf)):>24s}   {"reported":20s}')
    ch = [r for r in rows if r.get("chain")]
    if ch:
        through = [r for r in ch if "refused" not in r["chain"]]
        why = sorted({r["chain"]["refused"].split(":")[0] + " (" + r["family"] + ")" for r in ch if "refused" in r["chain"]})
        flags = sorted({f["flag"] for r in through for f in r["chain"]["flags"]})
        L.append(f"chain: {len(through)} renders graded through their reference session's recording chain, {len(ch) - len(through)} "
                 f"without one {why}; flags {flags}")
    blocked = ", ".join(sorted(BLOCKED_FAMILIES))
    L += ["", f'GATE {"VALID" if (ok_real and ok_null and ok_shift) else "NOT VALID"}: real controls {ok_real} (>= 80 % of graded; '
               f"refusals enumerated, not failures), nulls {ok_null}, pitch shift {ok_shift}",
          f"  no claim / renders not graded: the piano and every struck or plucked family (owner 2026-10-22)"
          + (f"; and sustained families whose real recordings fail Q2 Alive: {blocked} (not valid yet; pexfp.360 / clarinet calibration)" if blocked else "")]
    if ev:
        L += ["", "The candidates Q2 Alive was chosen from, library recordings only (identity_eval.py):", identity_eval.report(ev)]
    return "\n".join(L)

if __name__ == "__main__":
    main()
