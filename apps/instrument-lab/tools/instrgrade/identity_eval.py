#!/usr/bin/env python3
"""identity_eval: the UNREGISTERED Q2 candidates of identity.py measured on library recordings only, for the owner's choice.

No render and no recording from a second library is read here (those are gate_run.py's). Each library note sampled below is a
candidate against its own instrument's set with the note held out:
  leave-era-out     UIowa recorded the violin, viola, cello and bass twice, on different instruments in different sessions:
                    the set is the OTHER instrument's notes (the closest a single library comes to a second library);
  leave-singer-out  VocalSet: the set has no note of the candidate's singer;
  leave-file-out    every other instrument: the set has no note of the candidate's file (one run of a scale is near-copies).
For each candidate: candidates A, B and C of identity.py, the five nulls of freeze.null_signals judged against the same set, and a
+20 c copy (band-limited resampling). Also the confusion of A and B: the candidate judged against every OTHER instrument's set.

  phase 1   features of every library note        -> DIR/lib-{iowa,voice}.json
  phase 2   the sampled notes, their nulls and +20 c copies -> DIR/cands.json
  phase 3   the tables                             -> DIR/identity-eval.json, printed

usage (trj; numpy, pillow, matplotlib): identity_eval.py DIR PHASE [--jobs 12]
"""
from __future__ import annotations

import argparse
import json
import multiprocessing
import random
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import decompose as dc  # noqa: E402
import freeze  # noqa: E402
import identity as idn  # noqa: E402
from instrgrade import hz_of  # noqa: E402

SETS = {
    "iowa": dict(lib="/data/tmp/stormyhill/piano/lib-light-plus-piano.json", wav="/data/tmp/stormyhill/wav48",
                 refs=dc.DEFAULTS["refs_json"]),
    "voice": dict(lib="/data/tmp/wa-inst/p2-grade/voice/library.json", wav="/data/tmp/wa-inst/refs-wav48-voice",
                  refs="/data/tmp/wa-inst/p2-grade/voice/references-vocalset.json"),
}
CANON = {"altosax": "novib", "flute": "novib", "trumpet": "novib"}  # winds recorded with and without vibrato: the plain one
PER_GROUP = 30  # notes sampled per instrument, articulation and kind, spread over pitch and dynamics
SEED = 357
SHIFT_CENTS = 20.0
CANDIDATES = ("A", "B", "C")
REFUSALS = dc.REFUSALS + (FloatingPointError, ZeroDivisionError)


def key(e):
    return f'{e["file"]}@{e["start"]:.3f}'


def shift(y, cents):
    """`y` played `cents` higher by band-limited (FFT) resampling."""
    n = len(y)
    m = int(round(n / 2 ** (cents / 1200)))
    return np.fft.irfft(np.fft.rfft(y)[: m // 2 + 1], m) * (m / n)


_R = {}


def _init():
    for k, s in SETS.items():
        _R[k] = dc.Refs(s["lib"], s["wav"], s["refs"])


def lib_job(args):
    which, e = args
    try:
        return key(e), idn.features(_R[which].note(e)), None
    except REFUSALS as err:
        return key(e), None, f"{type(err).__name__}: {err}"


def as_candidate(raw, f0, kind, inh, label):
    """A signal graded as a render is: the render floor, the alignment, the grader's Note."""
    y, info = dc.align(raw, f0, dc.render_floor(dc.highpass(raw, f0), f0))
    return dc.Note(y, f0, kind, inh, label, info, max_cents=400.0)


def cand_job(args):
    which, e = args
    out = dict(key=key(e), which=which)
    try:
        refs = _R[which]
        y, _ = refs.aligned(e)
        n = refs.note(e)
        f0 = hz_of(e["note"])
        inh = e["inst"] in ("xylophone", "marimba")
        out["real"] = idn.features(n)
        nulls = freeze.null_signals(y, f0, n.B if e["kind"] == "decaying" and not inh else 0.0)
        if inh:
            nulls.pop("frozen")  # a bar's partials are not k f0: there is no frozen null of a bar
        out["nulls"] = {}
        for name, z in nulls.items():
            try:
                out["nulls"][name] = idn.features(as_candidate(z, f0, e["kind"], inh, name))
            except REFUSALS as err:
                out["nulls"][name] = dict(refused=f"{type(err).__name__}: {err}")
        try:
            out["shift"] = idn.features(as_candidate(shift(y, SHIFT_CENTS), f0, e["kind"], inh, "+20c"))
        except REFUSALS as err:
            out["shift"] = dict(refused=f"{type(err).__name__}: {err}")
    except REFUSALS as err:
        out["error"] = f"{type(err).__name__}: {err}"
    return out


def load_lib(which):
    return json.loads(Path(SETS[which]["lib"]).read_text())


def canon(e):
    return CANON.get(e["inst"], e["art"])


def sample(d):
    rng = random.Random(SEED)
    out = []
    for which in SETS:
        feats = json.loads((d / f"lib-{which}.json").read_text())["features"]
        groups = defaultdict(list)
        for e in load_lib(which):
            if key(e) in feats and e["art"] == canon(e):
                groups[(e["inst"], e["art"], e["kind"])].append(e)
        for _, es in sorted(groups.items()):
            es = sorted(es, key=lambda e: (e["note"], e["dyn"], e["file"]))
            if len(es) > PER_GROUP:
                step = len(es) / PER_GROUP
                es = [es[int(i * step + rng.random() * step)] for i in range(PER_GROUP)]
            out += [(which, e) for e in es]
    return out


def phase1(d, jobs):
    for which in SETS:
        lib = sorted(load_lib(which), key=lambda e: (e["file"], e["start"]))
        with ProcessPoolExecutor(jobs, mp_context=multiprocessing.get_context("fork"), initializer=_init) as pool:
            res = list(pool.map(lib_job, [(which, e) for e in lib], chunksize=16))
        feats = {k: f for k, f, _ in res if f is not None}
        (d / f"lib-{which}.json").write_text(json.dumps(dict(features=feats, refused={k: r for k, _, r in res if r is not None})))
        print(which, "features", len(feats), "refused", len(res) - len(feats), flush=True)


def phase2(d, jobs):
    cands = sample(d)
    print("sampled notes", len(cands), flush=True)
    with ProcessPoolExecutor(jobs, mp_context=multiprocessing.get_context("fork"), initializer=_init) as pool:
        res = list(pool.map(cand_job, cands, chunksize=1))
    (d / "cands.json").write_text(json.dumps(res))
    print("done", len(res), "errors", sum(1 for r in res if "error" in r), flush=True)


_RULERS = {}


def judge(c, members, ms):
    k = (tuple(la for la, _ in members), ms)
    if k not in _RULERS:
        _RULERS[k] = {m: idn.ruler(members, m) for m in ms}
    return idn.judge(c, members, ms, _RULERS[k])


def verdicts(c, S, others):
    """{A, B, C: True | False | None} of features c against own set S and the other instruments' sets `others`."""
    a, b = judge(c, S, idn.SET_MEASURES)["recognisable"], judge(c, S, idn.LIFE_MEASURES)["recognisable"]
    ns = idn.nearest_set(c, S, others) if len(S) >= idn.MIN_SET else None
    return dict(A=a, B=b, C=None if b is None or ns is None else bool(b and ns["own_nearest"]))


def phase3(d):
    libs = {w: json.loads((d / f"lib-{w}.json").read_text())["features"] for w in SETS}
    entries = {w: {key(e): e for e in load_lib(w) if key(e) in libs[w]} for w in SETS}
    pages = {}
    for w in SETS:
        for f in json.loads(Path(SETS[w]["refs"]).read_text())["files"]:
            pages[Path(f["path"]).stem.replace(" ", "")] = f["page"]
    page = lambda e: pages.get(Path(e["file"]).stem, "?")  # noqa: E731
    eras = defaultdict(set)
    for e in entries["iowa"].values():
        eras[e["inst"]].add(page(e))
    insts = {w: sorted({(x["inst"], canon(x)) for x in entries[w].values()}) for w in SETS}

    def members(w, e, protocol, inst=None, art=None):
        inst, art = inst or e["inst"], art or e["art"]
        out = []
        for k, x in entries[w].items():
            if k == key(e) or not idn.in_set(x, inst, art, e["kind"], e["note"]):
                continue
            if inst == e["inst"] and ((protocol == "leave-file-out" and x["file"] == e["file"])
                                      or (protocol == "leave-era-out" and page(x) == page(e))
                                      or (protocol == "leave-singer-out" and x["sul"] == e["sul"])):
                continue
            out.append((k, libs[w][k]))
        return out

    def protocol(w, e):  # the most distant hold-out the library allows for this note
        if w == "voice":
            return "leave-singer-out"
        return "leave-era-out" if len(eras[e["inst"]]) > 1 and e["kind"] == "sustained" else "leave-file-out"

    real = defaultdict(lambda: defaultdict(list))  # (group, protocol) -> candidate -> verdicts
    nulls = defaultdict(lambda: defaultdict(list))  # null -> candidate -> verdicts
    shifts = defaultdict(list)  # candidate -> |margin change| dB
    flips = defaultdict(int)
    confusion = defaultdict(lambda: defaultdict(list))
    by_kind = defaultdict(lambda: defaultdict(list))  # (kind, "real" | null name | "shift") -> candidate -> verdicts or dB
    res = json.loads((d / "cands.json").read_text())
    for r in res:
        if "error" in r:
            continue
        w = r["which"]
        e = entries[w][r["key"]]
        g = f'{e["inst"]}.{e["art"] or "-"}.{e["kind"]}'
        p = protocol(w, e)
        S = members(w, e, p)
        others = {inst: So for inst, art in insts[w] if inst != e["inst"] and (So := members(w, e, "other", inst, art))}
        v = verdicts(r["real"], S, others)
        for c in CANDIDATES:
            real[(g, p)][c].append(v[c])
            by_kind[(e["kind"], "real")][c].append(v[c])
        for name, nf in r["nulls"].items():
            vn = dict.fromkeys(CANDIDATES) if "refused" in nf else verdicts(nf, S, others)
            for c in CANDIDATES:
                nulls[name][c].append(vn[c])
                by_kind[(e["kind"], name)][c].append(vn[c])
        if "refused" not in r["shift"] and len(S) >= idn.MIN_SET:
            vs = verdicts(r["shift"], S, others)
            moved = {}
            for c, ms in (("A", idn.SET_MEASURES), ("B", idn.LIFE_MEASURES)):
                ja, jb = judge(r["real"], S, ms), judge(r["shift"], S, ms)
                if ja["margin_db"] is not None and jb["margin_db"] is not None:
                    moved[c] = max(abs(ja["measures"][m]["margin_db"] - jb["measures"][m]["margin_db"]) for m in ms)
                    shifts[c].append(moved[c])
                    by_kind[(e["kind"], "shift")][c].append(moved[c])
            na, nb = idn.nearest_set(r["real"], S, others), idn.nearest_set(r["shift"], S, others)
            if na and nb and "B" in moved:  # C moves by the larger of its two parts
                shifts["C"].append(max(moved["B"], abs(na["margin_db"] - nb["margin_db"])))
                by_kind[(e["kind"], "shift")]["C"].append(shifts["C"][-1])
            for c in CANDIDATES:
                flips[c] += int(v[c] is not None and vs[c] is not None and v[c] != vs[c])
        for inst, So in others.items():
            if len(So) >= idn.MIN_SET:
                for c, ms in (("A", idn.SET_MEASURES), ("B", idn.LIFE_MEASURES)):
                    confusion[g][c].append(judge(r["real"], So, ms)["recognisable"])

    def rate(vs):
        n = [x for x in vs if x is not None]
        return [sum(n), len(n)]

    out = dict(schema="frankensim.identity-eval.v1", identity_version=idn.IDENTITY_VERSION, utc=datetime.now(timezone.utc).isoformat(),
               per_group=PER_GROUP, seed=SEED, sets=SETS,
               real={f"{g}|{p}": {c: rate(vv[c]) for c in CANDIDATES} for (g, p), vv in sorted(real.items())},
               real_pooled={c: rate([x for vv in real.values() for x in vv[c]]) for c in CANDIDATES},
               real_cross_specimen={c: rate([x for (g, p), vv in real.items() if p != "leave-file-out" for x in vv[c]]) for c in CANDIDATES},
               nulls={name: {c: rate(vv[c]) for c in CANDIDATES} for name, vv in sorted(nulls.items())},
               shift={c: dict(median_db=float(np.median(v)), max_db=float(np.max(v)), n=len(v), share_over_1db=float(np.mean(np.array(v) >= 1.0)),
                              verdict_flips=flips[c]) for c, v in shifts.items() if v},
               confusion={g: {c: rate(vv[c]) for c in ("A", "B")} for g, vv in sorted(confusion.items())},
               confusion_pooled={c: rate([x for vv in confusion.values() for x in vv[c]]) for c in ("A", "B")},
               by_kind={kind: dict(real={c: rate(by_kind[(kind, "real")][c]) for c in CANDIDATES},
                                   nulls={c: rate([x for name in freeze.NULLS for x in by_kind[(kind, name)][c]]) for c in CANDIDATES},
                                   shift_over_1db={c: [int(np.sum(np.array(by_kind[(kind, "shift")][c]) >= 1.0)), len(by_kind[(kind, "shift")][c])]
                                                   for c in CANDIDATES},
                                   shift_max_db={c: float(np.max(by_kind[(kind, "shift")][c])) if by_kind[(kind, "shift")][c] else None
                                                 for c in CANDIDATES})
                        for kind in ("sustained", "decaying")})
    (d / "identity-eval.json").write_text(json.dumps(out, indent=1))
    print(report(out))


def pct(kn):
    k, n = kn
    return f"{k:3d}/{n:<3d} {100 * k / n:5.1f} %" if n else "  -/-        "


def report(o):
    L = [f"UNREGISTERED Q2 candidates on library recordings only ({o['identity_version']}; A set membership, B life in the set, "
         "C life + nearest-set identity)", "", f'{"real notes, held out":42s}' + "".join(f"{c:>18s}" for c in CANDIDATES)]
    for gp, cc in o["real"].items():
        L.append(f"  {gp.replace('|', '  '):40s}" + "".join(f"{pct(cc[c]):>18s}" for c in CANDIDATES))
    L.append(f'{"  ALL":42s}' + "".join(f"{pct(o['real_pooled'][c]):>18s}" for c in CANDIDATES))
    L.append(f'{"  across specimens/singers only":42s}' + "".join(f"{pct(o['real_cross_specimen'][c]):>18s}" for c in CANDIDATES))
    L += ["", f'{"nulls passing (must be 0)":42s}' + "".join(f"{c:>18s}" for c in CANDIDATES)]
    for name, cc in o["nulls"].items():
        L.append(f"  {name:40s}" + "".join(f"{pct(cc[c]):>18s}" for c in CANDIDATES))
    L += ["", "+20 c: largest margin change per note (dB): " + "   ".join(
        f"{c} median {s['median_db']:.3f} max {s['max_db']:.3f} ({100 * s['share_over_1db']:.1f} % >= 1 dB, {s['verdict_flips']} flips)"
        for c, s in o["shift"].items())]
    L.append("confusion, a real note passing as ANOTHER instrument: " + "   ".join(f"{c} {pct(o['confusion_pooled'][c])}" for c in ("A", "B")))
    for kind, k in o.get("by_kind", {}).items():
        L.append(f"{kind:9s} real " + "  ".join(f"{c} {pct(k['real'][c])}" for c in CANDIDATES)
                 + " | nulls passing " + "  ".join(f"{c} {pct(k['nulls'][c])}" for c in CANDIDATES)
                 + " | +20 c >= 1 dB " + "  ".join(f"{c} {pct(k['shift_over_1db'][c])} (max {k['shift_max_db'][c] or 0:.1f})" for c in CANDIDATES))
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir")
    ap.add_argument("phase", choices=("1", "2", "3"))
    ap.add_argument("--jobs", type=int, default=12)
    a = ap.parse_args()
    d = Path(a.dir)
    d.mkdir(parents=True, exist_ok=True)
    {"1": lambda: phase1(d, min(a.jobs, 16)), "2": lambda: phase2(d, min(a.jobs, 16)), "3": lambda: phase3(d)}[a.phase]()


if __name__ == "__main__":
    main()
