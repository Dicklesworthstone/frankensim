"""What a REAL instrument scores on the page's own rulers: the cross-library baseline, as one JSON file for the page.

The rulers (decompose.py) are built from one library: one room, one microphone, one player per instrument. So a real instrument from
somewhere else need not land inside them. This measures how often it does. Each real recording is graded as if it were a render,
with the same tool, against the same reference of the same note, with the unmodified rulers (pair.json under `certify/analysis`,
made by `certify/cert_prep_tinysol.py` / `cert_prep_voice.py` and `decompose.py batch --out ...`). This script only counts those
pair files; it grades nothing.

  inside fraction   inside / (inside + outside), pooled over every graded recording of the family and over every measure that has a
                    ruler (a measure without one gives no verdict and is counted apart as no_ruler). Pooled, not a mean of per-recording
                    fractions: the counts are in the file so the page can show "n of m".
  loose / strict    the two rulers of the page: all real-vs-real rows, and the same-note rows only (decompose.SAME_NOTE)
  badge             the page's badge rule applied to the pooled counts (inside only when nothing is outside; partly when at least half
                    of the ruled measures are inside; else outside; the weaker of loose and strict)
  our renders       the live pair of the same note as each recording, pooled the same way; one render pair counts once per recording
                    graded against its note, so it is repeated when several recordings share a note

Families with no recording here have no baseline; `families_without_baseline` lists those that have pairs on the page.
"""
from __future__ import annotations

import argparse
import collections
import datetime
import glob
import json
import statistics
import sys
from pathlib import Path

SCHEMA = "frankensim.instrument-cross-library-baseline.v1"
CERTIFY = Path("/data/tmp/wa-inst/p2/certify/analysis")
LIVE = Path("/data/tmp/wa-inst/analysis")
OUT = LIVE / "cross-library-baseline.json"

SOURCES = {
    "tinysol": {
        "prefix": "a real recording, not a render: TinySOL (IRCAM Studio On Line, CC BY 4.0)",
        "name": "TinySOL (IRCAM Studio On Line), v6.0",
        "url": "https://zenodo.org/records/3685367",
        "licence": "CC BY 4.0",
        "licence_url": "https://creativecommons.org/licenses/by/4.0/",
        "credit": "Cella, Ghisi, Lostanlen, Lévy, Fineberg and Maresz, IRCAM",
        "independence": "other_library",
        "versus": "UIowa Musical Instrument Samples (a different room, microphone and player)",
    },
    "vocalset": {
        "prefix": "a real recording, not a render: VocalSet (CC BY 4.0)",
        "name": "VocalSet: a singing voice dataset, singer m6, held out of the voice rulers",
        "url": "https://zenodo.org/records/1193957",
        "licence": "CC BY 4.0",
        "licence_url": "https://creativecommons.org/licenses/by/4.0/",
        "credit": "Wilkins, Seetharaman, Wahl and Pardo, 2018",
        "independence": "other_singer_same_library",
        "versus": "the other VocalSet singers (the same recording set, room and microphone; only the singer differs)",
    },
}
FAMILY_SOURCE = {"violin": "tinysol", "cello": "tinysol", "viola": "tinysol", "trumpet": "tinysol", "clarinet": "tinysol", "flute": "tinysol", "voice": "vocalset"}
RETUNED = "TinySOL's own authors digitally retuned the sample (the last field of its file name is not N)"

BADGES = ("Outside the real spread", "Partly inside the real spread", "Inside the real spread")


def badge(inside: int, outside: int) -> str | None:
    """The page's rule for one ruler: inside only when no measure is outside, partly when at least half of the ruled ones are inside."""
    if inside + outside == 0:
        return None
    return BADGES[2] if outside == 0 else BADGES[1] if inside / (inside + outside) >= 0.5 else BADGES[0]


def weaker(a: str | None, b: str | None) -> str | None:
    """The weaker of the loose and the strict badge (a ruler that gave no verdict does not speak)."""
    got = [x for x in (a, b) if x is not None]
    return min(got, key=BADGES.index) if got else None


def fraction(inside: int, outside: int) -> float | None:
    return None if inside + outside == 0 else round(inside / (inside + outside), 4)


def pooled(counts: dict) -> dict:
    i, o = counts["inside"], counts["outside"]
    return {"inside": i, "outside": o, "no_ruler": counts["no_ruler"], "inside_fraction": fraction(i, o)}


def source_of(pair: dict) -> str:
    """Which source table row a recording belongs to, from its own observer text: the table cannot silently disagree with the data."""
    obs = pair["ours"]["observer"]
    for key, s in SOURCES.items():
        if obs.startswith(s["prefix"]):
            return key
    raise SystemExit(f"{pair['family']}/{pair['note']}: observer text matches no known source: {obs[:120]}")


def real_note(slot: str, family: str) -> str:
    """The note a certification slot is a recording of: 'A4-mf-s2' -> 'A4'; a voice slot is its own page slot ('C3-a')."""
    return slot if family == "voice" else slot.split("-")[0]


def read_family(root: Path, live: Path, family: str) -> dict | None:
    files = sorted(glob.glob(str(root / family / "*" / "pair.json")))
    if not files:
        return None
    rows, refused, sources, versions = [], [], set(), set()
    real_c = {k: collections.Counter() for k in ("loose", "strict")}
    ours_c = {k: collections.Counter() for k in ("loose", "strict")}
    miss = {k: collections.Counter() for k in ("loose", "strict")}
    seen = {k: collections.Counter() for k in ("loose", "strict")}
    ruler_rows = {"loose": [], "strict": []}
    earns = collections.Counter()
    missing_live = 0
    for f in files:
        slot = Path(f).parent.name
        pair = json.loads(Path(f).read_text())
        if pair.get("refused"):
            refused.append({"slot": slot, "reason": pair["refused"]["reason"]})
            continue
        sources.add(source_of(pair))
        versions.add((pair["decompose"]["version"], pair["decompose"]["grader"]))
        verdicts = pair["verdicts"]
        for kind in ("loose", "strict"):
            for c in ("inside", "outside", "no_ruler"):
                real_c[kind][c] += verdicts[kind][c]
        for m in pair["metrics"]:
            if m["gap"] is None:
                continue
            ruler_rows["loose"].append(len(m["ruler"]))
            ruler_rows["strict"].append(len(m["ruler_strict"]))
            for kind, key in (("loose", "inside_loose"), ("strict", "inside_strict")):
                if m[key] is not None:
                    seen[kind][m["id"]] += 1
                    miss[kind][m["id"]] += m[key] is False
        both_clean = True
        for kind in ("loose", "strict"):
            clean = verdicts[kind]["outside"] == 0 and verdicts[kind]["inside"] > 0
            earns[kind] += clean
            both_clean &= clean
        earns["both"] += both_clean
        live_file = live / family / real_note(slot, family) / "pair.json"
        ours = json.loads(live_file.read_text()) if live_file.exists() else None
        if ours is None or ours.get("refused"):
            missing_live += 1
        else:
            for kind in ("loose", "strict"):
                for c in ("inside", "outside", "no_ruler"):
                    ours_c[kind][c] += ours["verdicts"][kind][c]
        stem = Path(pair["ours"]["wav"]).stem
        rows.append({
            "slot": slot, "note": real_note(slot, family), "recording": stem,
            "reference_file": pair["real"]["file"],
            "retuned_by_authors": family != "voice" and not stem.endswith("-N"),
            "measures": verdicts["measures"],
            "loose": verdicts["loose"], "strict": verdicts["strict"],
            "outside_loose": [m["id"] for m in pair["metrics"] if m["inside_loose"] is False],
            "outside_strict": [m["id"] for m in pair["metrics"] if m["inside_strict"] is False],
        })
    if not rows:
        return {"source": None, "recordings": 0, "refused": refused}
    assert len(sources) == 1, f"{family}: recordings from more than one source {sources}"
    source = SOURCES[next(iter(sources))]
    by_measure = {k: {"loose": {"outside": int(miss["loose"][k]), "ruled": seen["loose"][k]},
                      "strict": {"outside": int(miss["strict"][k]), "ruled": seen["strict"][k]}} for k in sorted(seen["loose"] | seen["strict"])}
    real = {k: pooled(real_c[k]) for k in ("loose", "strict")}
    ours = {k: pooled(ours_c[k]) for k in ("loose", "strict")}
    for who in (real, ours):
        who["badge_by_page_rule"] = weaker(badge(who["loose"]["inside"], who["loose"]["outside"]), badge(who["strict"]["inside"], who["strict"]["outside"]))
    return {
        "source": {k: source[k] for k in ("name", "url", "licence", "licence_url", "credit")},
        "independence": source["independence"], "compared_with": source["versus"],
        "grader": sorted(f"{v}/{g}" for v, g in versions),
        "recordings": len(rows), "refused": refused, "retuned_by_authors": sum(r["retuned_by_authors"] for r in rows),
        "real_recordings": real,
        "recordings_with_nothing_outside": {k: int(earns[k]) for k in ("loose", "strict", "both")},
        "our_renders_same_notes": {**ours, "render_pairs": len(rows) - missing_live, "recordings_whose_note_has_no_live_render": missing_live},
        "ruler_rows_per_measure": {k: {"median": statistics.median(v), "min": min(v), "max": max(v)} for k, v in ruler_rows.items()},
        "by_measure": by_measure, "recording_list": rows,
    }


CAVEATS = [
    "A cross-library baseline is a handful of recordings per family: read the counts, not only the percentages.",
    "Each real recording is graded against ONE reference of its note, and the rulers are built from a median of only a few real-vs-real "
    "distances per measure (ruler_rows_per_measure), so a ruler is a weak bound and a real recording can land outside it by sampling alone.",
    "Only the mf dynamic was graded; pp and ff were not.",
    "Some TinySOL samples were digitally retuned by their authors (retuned_by_authors counts them), and the TinySOL notes sit closer to "
    "equal temperament than the UIowa references do, so f0_cents partly tests a tuning standard rather than an instrument.",
    "The voice baseline is a different SINGER from the same library (VocalSet), not a different library: same room, same microphone. "
    "The held-out singer is not a row of any voice ruler. It is a weaker test of 'a different library' than the TinySOL families.",
    "The counts are of measures inside a ruler, not of audible similarity; a recording can sit inside every ruler and still sound different.",
]


def build(root: Path = CERTIFY, live: Path = LIVE) -> dict:
    families = {}
    for family in sorted(p.name for p in root.iterdir() if p.is_dir()):
        got = read_family(root, live, family)
        if got is not None:
            families[family] = got
    page_families = sorted({p.parent.parent.name for p in live.glob("*/*/pair.json")})
    all_c = collections.Counter()
    out_c = collections.Counter()
    for fam in families.values():
        for k, v in fam.get("by_measure", {}).items():
            for kind in ("loose", "strict"):
                all_c[(kind, k)] += v[kind]["ruled"]
                out_c[(kind, k)] += v[kind]["outside"]
    all_measures = {k: {kind: {"outside": out_c[(kind, k)], "ruled": all_c[(kind, k)]} for kind in ("loose", "strict")} for k in sorted({k for _, k in all_c})}
    return {
        "schema": SCHEMA,
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "question": "How much of the page's rulers does a REAL instrument from somewhere else land inside?",
        "method": __doc__.strip().split("\n\n")[0:2],
        "page_rule": "Inside the real spread only when no ruled measure is outside; Partly inside when at least half of the ruled measures are inside; otherwise Outside; the weaker of the loose and the strict ruler.",
        "families": families,
        "families_without_baseline": [f for f in page_families if f not in families],
        "all_recordings_by_measure": all_measures,
        "caveats": CAVEATS,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--certify", type=Path, default=CERTIFY, help="the certification pair tree (default %(default)s)")
    ap.add_argument("--live", type=Path, default=LIVE, help="the page's analysis tree, for our renders of the same notes (default %(default)s)")
    ap.add_argument("--out", type=Path, default=OUT, help="the file to write (default %(default)s)")
    args = ap.parse_args(argv)
    doc = build(args.certify, args.live)
    tmp = args.out.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + "\n")
    tmp.replace(args.out)
    for fam, f in doc["families"].items():
        r, o = f["real_recordings"], f["our_renders_same_notes"]
        print(f"{fam:9s} {f['recordings']:2d} recordings  real loose {r['loose']['inside']}/{r['loose']['inside'] + r['loose']['outside']}  strict {r['strict']['inside']}/{r['strict']['inside'] + r['strict']['outside']}"
              f"   ours loose {o['loose']['inside']}/{o['loose']['inside'] + o['loose']['outside']}  strict {o['strict']['inside']}/{o['strict']['inside'] + o['strict']['outside']}")
    print("wrote", args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
