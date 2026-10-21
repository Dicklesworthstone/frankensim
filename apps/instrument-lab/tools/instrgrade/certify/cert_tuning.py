"""Tuning of each real candidate and of the UIowa/VocalSet reference it is graded against, in cents re equal temperament (A4 = 440)."""
import glob
import json
import math
import statistics as st

NAMES = {"C": 0, "D": 2, "E": 4, "F": 5, "G": 7, "A": 9, "B": 11}


def et(note):
    n = NAMES[note[0]] + {"b": -1, "#": 1}.get(note[1], 0) if len(note) > 2 and note[1] in "b#" else NAMES[note[0]]
    octv = int(note[-1])
    return 440.0 * 2 ** ((12 * (octv + 1) + n - 69) / 12)


for fam in ("violin", "trumpet", "voice", "cello", "viola", "clarinet", "flute"):
    cand, ref, gaps, rulers = [], [], [], []
    for p in sorted(glob.glob(f"/data/tmp/wa-inst/p2/certify/analysis/{fam}/*/pair.json")):
        j = json.load(open(p))
        if j.get("refused"):
            continue
        slot = p.split("/")[-2]
        note = slot.split("-")[0]
        m = {x["id"]: x for x in j["metrics"]}["f0_cents"]
        f = et(note)
        cand.append(1200 * math.log2(m["ours"] / f)); ref.append(1200 * math.log2(m["real"] / f)); gaps.append(m["gap"])
        rulers.append(max((r["value"] for r in m["ruler"]), default=float("nan")))
    if not cand:
        continue
    print(f"{fam:8s} n={len(cand):2d}  candidate vs ET: median |{st.median([abs(c) for c in cand]):5.1f}| c (range {min(cand):+.1f}..{max(cand):+.1f})   reference vs ET: median |{st.median([abs(r) for r in ref]):5.1f}| c (range {min(ref):+.1f}..{max(ref):+.1f})   "
          f"f0 gap candidate-reference: median {st.median(gaps):5.1f} c; f0 ruler max median {st.median([r for r in rulers if r == r] or [float('nan')]):5.1f} c")
