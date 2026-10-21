"""Would a real recording from another library pass the grader? Read off the pair files of /data/tmp/wa-inst/p2/certify/analysis."""
import collections
import glob
import json
import statistics as st
import sys

ROOT = "/data/tmp/wa-inst/p2/certify/analysis"
LIVE = "/data/tmp/wa-inst/analysis"
only_dyn = sys.argv[1] if len(sys.argv) > 1 else None  # e.g. "mf": restrict to slots with that dynamic in the name

for fam in ("violin", "trumpet", "voice", "cello", "viola", "clarinet", "flute"):
    rows = []
    for p in sorted(glob.glob(f"{ROOT}/{fam}/*/pair.json")):
        j = json.load(open(p))
        slot = p.split("/")[-2]
        if only_dyn and fam != "voice" and f"-{only_dyn}" not in slot:
            continue
        if j.get("refused"):
            rows.append((slot, None, j["refused"]["reason"]))
            continue
        rows.append((slot, j, None))
    if not rows:
        continue
    print(f"\n==== {fam}: {len(rows)} real recordings from another source graded against the UIowa/VocalSet reference of the same note")
    agg_l, agg_s = collections.Counter(), collections.Counter()
    nl, ns = collections.Counter(), collections.Counter()
    loose_in, strict_in, render_loose, render_strict, measures = [], [], [], [], []
    pool = collections.Counter()
    nrow_l, nrow_s = [], []
    for slot, j, why in rows:
        if j is None:
            print(f"  {slot:12s} REFUSED: {why[:150]}")
            continue
        note = slot.split("-")[0] if fam != "voice" else slot
        live = json.load(open(f"{LIVE}/{fam}/{note}/pair.json"))
        v, lv = j["verdicts"], live["verdicts"]
        out_l = [m["id"] for m in j["metrics"] if m["inside_loose"] is False]
        out_s = [m["id"] for m in j["metrics"] if m["inside_strict"] is False]
        for m in j["metrics"]:
            if m["gap"] is not None:
                nrow_l.append(len(m["ruler"])); nrow_s.append(len(m["ruler_strict"]))
            if m["inside_loose"] is not None:
                nl[m["id"]] += 1
                agg_l[m["id"]] += m["inside_loose"] is False
            if m["inside_strict"] is not None:
                ns[m["id"]] += 1
                agg_s[m["id"]] += m["inside_strict"] is False
        frac = lambda c: None if c["inside"] + c["outside"] == 0 else c["inside"] / (c["inside"] + c["outside"])  # noqa: E731  (a measure with no ruler gives no verdict)
        for kind in ("loose", "strict"):
            pool["real_" + kind + "_in"] += v[kind]["inside"]; pool["real_" + kind + "_out"] += v[kind]["outside"]
            pool["render_" + kind + "_in"] += lv[kind]["inside"]; pool["render_" + kind + "_out"] += lv[kind]["outside"]
        pool["notes_real"] += 1
        pool.setdefault("notes", set()) if False else None
        loose_in.append(frac(v["loose"])); strict_in.append(frac(v["strict"])); measures.append(v["measures"])
        render_loose.append(frac(lv["loose"])); render_strict.append(frac(lv["strict"]))
        f0 = {m["id"]: m for m in j["metrics"]}["f0_cents"]
        cs = lambda c: f"{c['inside']:2d}/{c['inside'] + c['outside']:<2d}" + (f"(+{c['no_ruler']} no ruler)" if c["no_ruler"] else "")  # noqa: E731
        print(f"  {slot:12s} loose {cs(v['loose'])} strict {cs(v['strict'])}   our render ({live['ours']['git_sha'][:9]}): loose {cs(lv['loose'])} strict {cs(lv['strict'])}   f0 gap {f0['gap']:5.1f} c (ruler max {max((r['value'] for r in f0['ruler']), default=float('nan')):5.1f})   OUTSIDE loose: {','.join(out_l)}")
    if not loose_in:
        continue
    mdn = lambda xs: "-" if not [x for x in xs if x is not None] else f"{st.median([x for x in xs if x is not None]):.2f} (range {min(x for x in xs if x is not None):.2f}-{max(x for x in xs if x is not None):.2f}, n={len([x for x in xs if x is not None])})"  # noqa: E731
    print(f"  -- FRACTION of the measures that have a ruler that a real recording from another source lands inside: loose median {mdn(loose_in)}; strict median {mdn(strict_in)}")
    print(f"  -- the same fraction for OUR RENDERS of the same notes:                                             loose median {mdn(render_loose)}; strict median {mdn(render_strict)}")
    level = lambda i, o: "Inside the real spread" if o == 0 else ("Partly inside the real spread" if i / (i + o) >= 0.5 else "Outside the real spread")  # noqa: E731  (P8's rule)
    weaker = lambda a_, b_: min((a_, b_), key=("Outside the real spread", "Partly inside the real spread", "Inside the real spread").index)  # noqa: E731
    for who in ("real", "render"):
        li, lo_, si, so = pool[who + "_loose_in"], pool[who + "_loose_out"], pool[who + "_strict_in"], pool[who + "_strict_out"]
        print(f"  -- POOLED over these notes ({'all ' + str(pool['notes_real']) + ' real recordings' if who == 'real' else 'our ' + str(pool['notes_real']) + ' render pairs, one per recording above, repeats counted'}): {who} inside loose {li} of {li + lo_}, strict {si} of {si + so}  ->  the page's rule gives \"{weaker(level(li, lo_), level(si, so))}\"")
    print(f"  -- the rulers these are judged by are built from a median of {st.median(nrow_l)} real-vs-real distances per measure (loose; range {min(nrow_l)}-{max(nrow_l)}) and {st.median(nrow_s)} (strict; range {min(nrow_s)}-{max(nrow_s)}); the verdict is against the LARGEST")
    print("  -- measure by measure: real recordings OUTSIDE the loose ruler / outside the strict ruler, of n with a ruler")
    for k in sorted(nl, key=lambda k: -agg_l[k] / nl[k]):
        print(f"       {k:18s} loose {agg_l[k]:2d}/{nl[k]:2d}   strict {agg_s[k]:2d}/{ns[k]:2d}")

# ---- all families together: which measures does a real recording from another source miss?
tot_out, tot_n, tot_out_s, tot_n_s = collections.Counter(), collections.Counter(), collections.Counter(), collections.Counter()
nrec, per_fam = 0, {}
for fam in ("violin", "trumpet", "voice", "cello", "viola", "clarinet", "flute"):
    for p in glob.glob(f"{ROOT}/{fam}/*/pair.json"):
        j = json.load(open(p))
        if j.get("refused"):
            continue
        nrec += 1
        for m in j["metrics"]:
            if m["inside_loose"] is not None:
                tot_n[m["id"]] += 1; tot_out[m["id"]] += m["inside_loose"] is False
            if m["inside_strict"] is not None:
                tot_n_s[m["id"]] += 1; tot_out_s[m["id"]] += m["inside_strict"] is False
print(f"\n==== ALL {nrec} real recordings from another source, 7 families: share OUTSIDE the ruler, per measure (loose | strict)")
for k in sorted(tot_n, key=lambda k: -tot_out[k] / tot_n[k]):
    print(f"   {k:18s} {tot_out[k]:2d}/{tot_n[k]:2d} = {100 * tot_out[k] / tot_n[k]:3.0f}%  |  {tot_out_s[k]:2d}/{tot_n_s[k]:2d} = {100 * tot_out_s[k] / max(tot_n_s[k], 1):3.0f}%")
