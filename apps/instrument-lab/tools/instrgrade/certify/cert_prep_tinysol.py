"""TinySOL leg of the certification: real notes from another library (IRCAM Studio On Line, CC BY 4.0), graded as if they were renders
against the UIowa reference of each note our renders are graded against. Families: violin, trumpet (default), viola, cello, clarinet, flute.

usage: cert_prep_tinysol.py DYNS [FAMS]   (e.g. mf  or  pp,ff;  FAMS comma-separated)   -> extracts, converts, writes
cands-tinysol-<DYNS>[-<FAMS>].json for `decompose.py batch --out certify/analysis`; cross_library_baseline.py then counts the pairs."""
import csv
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, "/data/tmp/wa-inst/p2/tools")
OUT = Path("/data/tmp/wa-inst/p2/certify")
DL = OUT / "tinysol-dl"
dyns = sys.argv[1].split(",")
fams = sys.argv[2].split(",") if len(sys.argv) > 2 else ["violin", "trumpet"]
FLAT = {"Db": "C#", "Eb": "D#", "Gb": "F#", "Ab": "G#", "Bb": "A#"}


def sharp(n):  # TinySOL names the black keys with sharps
    return FLAT.get(n[:2], n[:2]) + n[2:] if n[:2] in FLAT else n


INST = {"violin": "Vn", "trumpet": "TpC", "viola": "Va", "cello": "Vc", "clarinet": "ClBb", "flute": "Fl"}
TARGETS = {f: sorted(p.split("/")[-2] for p in __import__("glob").glob(f"/data/tmp/wa-inst/analysis/{f}/*/pair.json")
                     if json.load(open(p)).get("ours") and not json.load(open(p)).get("refused") and json.load(open(p)).get("real")) for f in fams}
LIB = json.load(open("/data/tmp/wa-inst/p2-grade/grade/library-v6-light.json"))
rows = list(csv.DictReader(open(DL / "TinySOL_metadata.csv")))

want = []
for fam, notes in TARGETS.items():
    for note in notes:
        pair = json.load(open(f"/data/tmp/wa-inst/analysis/{fam}/{note}/pair.json"))
        fname, rest = pair["decompose"]["library_note"].split("@", 1)
        name, dyn = rest.rsplit(".", 1)
        ref = [e for e in LIB if e["file"] == fname and e["name"] == name and e["dyn"] == dyn][0]
        for r in rows:
            if r["Instrument (abbr.)"] == INST[fam] and r["Pitch"] == sharp(note) and r["Dynamics"] in dyns:
                want.append((fam, note, ref, pair["alignment"]["band_hz"], r))

(OUT / "tinysol").mkdir(exist_ok=True)
(OUT / "tinysol-wav48").mkdir(exist_ok=True)
names = sorted({Path(w[4]["Path"]).name for w in want})
todo = [n for n in names if not list((OUT / "tinysol").rglob(n))]
if todo:  # one pass over the archive for everything missing
    subprocess.run(["tar", "-xzf", str(DL / "TinySOL.tar.gz"), "-C", str(OUT / "tinysol"), "--wildcards"] + ["*/" + n for n in todo], check=True)
cands = []
for fam, note, ref, band, r in want:
    src = next((OUT / "tinysol").rglob(Path(r["Path"]).name))
    stem = src.stem
    dst = OUT / "tinysol-wav48" / (stem + ".wav")
    if not dst.exists():
        subprocess.run(["ffmpeg", "-nostdin", "-v", "error", "-y", "-i", str(src), "-ac", "1", "-ar", "48000", "-c:a", "pcm_f32le", str(dst)], check=True)
    string = r["String ID (if applicable)"]
    slot = f"{note}-{r['Dynamics']}" + (f"-s{int(float(string))}" if string else "")
    retuned = r["Needed digital retuning"] == "TRUE"
    cands.append(dict(family=fam, note=note, slot=slot, dyn=ref["dyn"], art=ref["art"], ref_contains=Path(ref["file"]).stem, band_hz=band, wav=str(dst),
                      git_sha="tinysol", observer=f"a real recording, not a render: TinySOL (IRCAM Studio On Line, CC BY 4.0) {r['Path']}"
                      + (", digitally retuned by its authors" if retuned else ""), recipe="none: recorded"))
    print(f"{fam:8s} {slot:12s} {r['Path'].split('/')[-1]:34s} reference {ref['file']:36s} band {band} {'RETUNED' if retuned else ''}")
json.dump(cands, open(OUT / f"cands-tinysol-{'-'.join(dyns)}-{'-'.join(fams)}.json", "w"), indent=1)
print(len(cands), "candidates")
