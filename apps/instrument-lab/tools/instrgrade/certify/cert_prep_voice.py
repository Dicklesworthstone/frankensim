"""Voice leg of the certification: one VocalSet singer (m6) held OUT of the library, graded as if he were a render against the
reference singer of each graded slot, with the ruler built from the other singers only.

Writes /data/tmp/wa-inst/p2/certify/voice-lib-no-m6.json, the clips under voice-cands-wav/ and cands-voice-m6.json."""
import json
import sys
from pathlib import Path

sys.path.insert(0, "/data/tmp/wa-inst/p2/tools")
import decompose as D  # noqa: E402
import numpy as np  # noqa: E402

OUT = Path("/data/tmp/wa-inst/p2/certify")
LIB = json.load(open("/data/tmp/wa-inst/p2-grade/voice/library.json"))
HOLD = "m6_"
WAV = Path("/data/tmp/wa-inst/refs-wav48-voice")
PAD_BEFORE, PAD_AFTER, MIN_GAP = 0.4, 0.1, 0.45

json.dump([e for e in LIB if not e["file"].startswith(HOLD)], open(OUT / "voice-lib-no-m6.json", "w"))
(OUT / "voice-cands-wav").mkdir(exist_ok=True)
cands, skipped = [], []
for p in sorted(Path("/data/tmp/wa-inst/analysis/voice").glob("*/pair.json")):
    slot = p.parent.name
    if slot == "C4-u-rigid":
        continue
    j = json.load(open(p))
    fname, rest = j["decompose"]["library_note"].split("@", 1)
    name, dyn = rest.rsplit(".", 1)
    ref = [e for e in LIB if e["file"] == fname and e["name"] == name and e["dyn"] == dyn][0]
    # the held-out singer's recording of the same vowel (and the same dynamic) at the same pitch
    hits = [e for e in LIB if e["file"].startswith(HOLD) and e["art"] == ref["art"] and e["name"] == name and e["dyn"] == dyn]
    if not hits:
        skipped.append((slot, f"{HOLD[:-1]} has no {ref['art']} at {name} {dyn}"))
        continue
    c = hits[0]
    same = sorted((e for e in LIB if e["file"] == c["file"]), key=lambda e: e["start"])
    k = same.index(c)
    if k > 0 and c["start"] - same[k - 1]["end"] < MIN_GAP:
        skipped.append((slot, f"previous note in {c['file']} ends {c['start'] - same[k - 1]['end']:.2f} s before"))
        continue
    x = D.load(WAV / c["file"])
    a, b = max(0, int((c["start"] - PAD_BEFORE) * D.SR)), min(len(x), int((c["end"] + PAD_AFTER) * D.SR))
    clip = OUT / "voice-cands-wav" / f"{slot}-{Path(c['file']).stem}.wav"
    D.write_wav(clip, x[a:b] / max(1e-9, float(np.abs(x[a:b]).max())) * 0.7)
    cands.append(dict(family="voice", note=name, slot=slot, dyn=dyn, art=ref["art"], ref_contains=Path(fname).stem, band_hz=None, wav=str(clip),
                      git_sha="vocalset-" + c["file"][:2], observer=f"a real recording, not a render: VocalSet (CC BY 4.0) {c['file']}, held out of the ruler",
                      recipe="none: recorded"))
    print(f"{slot:6s} reference {fname:26s} candidate {c['file']:26s} {c['start']:.2f}-{c['end']:.2f} s")
json.dump(cands, open(OUT / "cands-voice-m6.json", "w"), indent=1)
print(len(cands), "candidates;", len(skipped), "skipped:")
for s in skipped:
    print("  skipped", s)
