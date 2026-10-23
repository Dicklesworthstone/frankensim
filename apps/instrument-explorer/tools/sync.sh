#!/usr/bin/env bash
# Mirror the instrument program's renders, grades and screened reference recordings into a fresh
# ~/.local/share/instrument-explorer/mirror-<UTC>/ for tools/build.py. Read-only on the hosts.
#
# A fresh directory per mirror, because an rsync into an old mirror keeps files the hosts have since moved away.
# Files unchanged since the previous mirror are hard-linked from it (rsync --link-dest), so a re-sync
# transfers only what changed. A step that fails is reported and the others still run.
#
# Hosts: IE_TRJ (default trj-ts) holds the grading tree and the recordings; IE_TS1 (default ts1) holds lane renders.
set -u
ROOT="$HOME/.local/share/instrument-explorer"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
GIT=/Users/jemanuel/bin/git
[ -x "$GIT" ] || GIT=git
TRJ="${IE_TRJ:-trj-ts}"
TS1="${IE_TS1:-ts1}"
SSH="ssh -o ConnectTimeout=15"

# Work-in-progress renders by lanes that have not committed or graded them yet. Edit when the lanes move them.
KIT_DECAYS="/data/tmp/proudbison-370/dw3/"                                          # on TRJ: dw3-decays.json
KIT_RENDERS="/home/ubuntu/irp-verify/out-ProudBison/dw3-*.wav"                      # on TS1
VIOLIN_ASSEMBLED="/home/ubuntu/irp-verify/out-ProudBison/violin-g3-rung-v2/violin-g3-rung-2s.wav"  # on TS1
FLUTE_EXTERIOR="/data/tmp/wa-inst/renders/flute-exterior/25a5514944ea/"             # on TRJ

PREV="$(ls -d "$ROOT"/mirror-* 2>/dev/null | tail -1)"
M="$ROOT/mirror-$(date -u +%Y%m%dT%H%MZ)"
mkdir -p "$M/analysis" "$M/piano-pairs" "$M/u4" "$M/kit-dw3/renders" "$M/renders/flute-exterior/25a5514944ea" "$M/renders/violin-assembled"
FAILED=0
step() { echo "== $(date -u +%H:%M:%SZ) $*"; }
done_step() { local rc=$1; echo "   exit $rc"; [ "$rc" -eq 0 ] || FAILED=$((FAILED + 1)); }
linkdest() { if [ -n "$PREV" ] && [ -d "$PREV/$1" ]; then echo "--link-dest=$PREV/$1"; fi; }

step "reference list (references-screened.json at origin/main)"
"$GIT" -C "$REPO" fetch -q origin
"$GIT" -C "$REPO" show origin/main:data/instruments/references-screened.json > "$M/references-screened.json"
python3 -I -c '
import json, sys
d = json.load(open(sys.argv[1]))
for f in d["files"]:
    for owner in [f] + (f.get("members") or []):
        for c in owner.get("clips") or []:
            if c.get("status") == "VALID":
                print("clips/" + c["family"] + "/" + c["id"] + ".wav")
' "$M/references-screened.json" > "$M/u4-valid-files.txt"
done_step $?
echo "   $(wc -l < "$M/u4-valid-files.txt") valid clips listed"

step "grading tree (current files, and only pair.json of superseded gradings)"
rsync -a --prune-empty-dirs -e "$SSH" $(linkdest analysis) \
  --include='*/' --include='superseded/**/pair.json' --exclude='superseded/**' \
  --include='*.json' --include='*.wav' --include='*.webp' --include='*.txt' --exclude='*' \
  "$TRJ:/data/tmp/wa-inst/analysis/" "$M/analysis/"
done_step $?

step "piano pairs"
rsync -a -e "$SSH" $(linkdest piano-pairs) --exclude='*.f64le' "$TRJ:/data/tmp/stormyhill/piano-pairs/" "$M/piano-pairs/"
done_step $?

step "screened reference recordings"
rsync -a -e "$SSH" $(linkdest u4) --files-from="$M/u4-valid-files.txt" "$TRJ:/data/tmp/wa-inst/refs/u4/" "$M/u4/"
done_step $?

step "kit lane: band decays and renders (work in progress)"
rsync -a -e "$SSH" --include='*.json' --include='*.log' --exclude='*' "$TRJ:$KIT_DECAYS" "$M/kit-dw3/"
done_step $?
rsync -a -e "$SSH" "$TS1:$KIT_RENDERS" "$M/kit-dw3/renders/"
done_step $?

step "assembled violin render (work in progress)"
rsync -a -e "$SSH" "$TS1:$VIOLIN_ASSEMBLED" "$M/renders/violin-assembled/"
done_step $?

step "flute renders with the exterior body (ungraded)"
rsync -a -e "$SSH" --include='*.wav' --include='*.provenance.json' --exclude='*' "$TRJ:$FLUTE_EXTERIOR" "$M/renders/flute-exterior/25a5514944ea/"
done_step $?

step "done: $M ($FAILED failed steps)"
echo "next: python3 $HERE/build.py --mirror $M"
[ "$FAILED" -eq 0 ]
