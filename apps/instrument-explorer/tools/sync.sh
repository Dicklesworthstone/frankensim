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
KIT="/data/tmp/mb-phrases/kit16/"  # on TRJ: source.json (render commits), renders/, grade/ (kitgrade.py: data/kit.json, data/decays.json, audio/*.wav, img/*.webp)
VIOLIN_ASSEMBLED="/home/ubuntu/irp-verify/out-ProudBison/violin-g3-rung-v2/violin-g3-rung-2s.wav"  # on TS1
FLUTE_EXTERIOR="/data/tmp/wa-inst/renders/flute-exterior/25a5514944ea/"             # on TRJ
TIMPANI="/home/ubuntu/irp-verify/sb-timpani-full/timpani-2b-4s.wav"  # on TS1: kettle B at a95d54fa4, timpani/render.sh's card and command, 192000 frames (pexfp.376)
PLAYABLE="/data/tmp/mb-phrases/scales18/"   # on TRJ: every note each family can play today, per family dir with its own source.json (the commit it was rendered at), and pitch.json (measured f0). scales18 = scales17 with the three notes that changed at 30fa619b4816 (trumpet pramp-t1.25e0, voice a-A3, viola open G3)
# On TRJ: regradings made outside the shared tree with the committed grader. Each <family>/<note>/ here replaces that note
# of the grading tree (its superseded/ carries the earlier gradings forward). Every graded pair the site shows was regraded
# on 2026-10-29 by instrgrade-4 (the HNR fix, f95e8c6c6; ig4b/run-ig4b.sh), the reference-only notes' rulers included.
# The renders are those of grade-overrides-1028b: the violin's four open strings on the assembled violin (material_strings
# at c6854036b), the flute's D4, G4 and A5 with the player's face and one-sided lip (music_render at 29b7a3088), the guitar's
# five notes from its committed recipe (render.sh at 2cc845201, pexfp.377; grade-g377land), the clarinet's 21 notes with the
# measured orifice law (music_render at e60c12cf3), the voice's F4-e and C4-u-rigid (music_render at e2414832f), and the
# grading tree's own. 1028b's refused pairs and the recorder's and tin whistle's listening pages (no reference, so no
# verdicts) are carried forward unregraded.
GRADE_OVERRIDES="/data/tmp/mb-phrases/grade-overrides-ig4g/"

PREV="$(ls -d "$ROOT"/mirror-* 2>/dev/null | tail -1)"
M="$ROOT/mirror-$(date -u +%Y%m%dT%H%MZ)"
mkdir -p "$M/analysis" "$M/piano-pairs" "$M/u4" "$M/kit" "$M/renders/flute-exterior/25a5514944ea" "$M/renders/violin-assembled" "$M/renders/timpani/a95d54fa4" "$M/renders/playable"
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

step "regradings that replace notes of the grading tree ($GRADE_OVERRIDES)"
for note in $($SSH "$TRJ" "cd $GRADE_OVERRIDES && ls -d */*/ 2>/dev/null" | tr -d '\r'); do
  # A fresh directory per note, so no file of the shared tree's grading survives beside the regrading's.
  mv "$M/analysis/$note" "$M/analysis-replaced-$(echo "$note" | tr / _)" 2>/dev/null
  mkdir -p "$M/analysis/$note"
  rsync -a --prune-empty-dirs -e "$SSH" \
    --include='*/' --include='superseded/**/pair.json' --exclude='superseded/**' \
    --include='*.json' --include='*.wav' --include='*.webp' --include='*.txt' --exclude='*' \
    "$TRJ:$GRADE_OVERRIDES$note" "$M/analysis/$note" || FAILED=$((FAILED + 1))
  echo "   $note"
done

step "piano pairs"
# The graded set re-rendered 2026-10-28 with the radiation-fit storage damping (acb563d4c with fit.rs as landed in
# 88b352677, bead i0cam) and graded by instrgrade-3 (eec0dd6e9). Before: /data/tmp/stormyhill/piano-pairs.
PIANO_PAIRS=/data/tmp/mb-phrases/piano-pairs-ig3
KEPT=""
for d in $(ls -d "$ROOT"/mirror-* 2>/dev/null | grep -v "^$M\$"); do [ -n "$(ls -A "$d/piano-pairs" 2>/dev/null)" ] && KEPT=$d; done
if $SSH "$TRJ" "test -d $PIANO_PAIRS"; then
  rsync -a -e "$SSH" $(linkdest piano-pairs) --exclude='*.f64le' "$TRJ:$PIANO_PAIRS/" "$M/piano-pairs/"
  done_step $?
elif [ -n "$KEPT" ]; then
  # The lane's directory left the host (gone from trj between 07:10Z and 15:13Z on 2026-10-26); the pairs are finished
  # gradings, so the newest mirrored copy stands, carried forward hard-linked and named here.
  echo "   $TRJ:$PIANO_PAIRS is gone; keeping the copy from $(basename "$KEPT") ($(find "$KEPT/piano-pairs" -type f | wc -l | tr -d ' ') files)"
  rsync -a --link-dest="$KEPT/piano-pairs" "$KEPT/piano-pairs/" "$M/piano-pairs/"
  done_step $?
else
  echo "   $TRJ:$PIANO_PAIRS is gone and no earlier mirror holds the pairs"
  done_step 1
fi

step "screened reference recordings"
rsync -a -e "$SSH" $(linkdest u4) --files-from="$M/u4-valid-files.txt" "$TRJ:/data/tmp/wa-inst/refs/u4/" "$M/u4/"
done_step $?

step "kit: renders at known commits, graded piece by piece against the Salamander Drumkit"
rsync -a -e "$SSH" $(linkdest kit) --include='/source.json' --include='/grade/***' --include='/renders/' --include='/renders/run*.log' --exclude='*' "$TRJ:$KIT" "$M/kit/"
done_step $?

step "assembled violin render (work in progress)"
rsync -a -e "$SSH" "$TS1:$VIOLIN_ASSEMBLED" "$M/renders/violin-assembled/"
done_step $?

step "flute renders with the exterior body (ungraded)"
rsync -a -e "$SSH" --include='*.wav' --include='*.provenance.json' --exclude='*' "$TRJ:$FLUTE_EXTERIOR" "$M/renders/flute-exterior/25a5514944ea/"
done_step $?

step "timpani render at its landing commit (not graded)"
rsync -a -e "$SSH" "$TS1:$TIMPANI" "$M/renders/timpani/a95d54fa4/timpani-2b.wav"
done_step $?

step "playable notes per family (listening renders, not graded)"
rsync -a -e "$SSH" --include='*/' --include='*.wav' --include='*.provenance.json' --include='source.json' --include='pitch.json' --exclude='*' "$TRJ:$PLAYABLE" "$M/renders/playable/"
done_step $?

step "done: $M ($FAILED failed steps)"
echo "next: python3 $HERE/build.py --mirror $M"
[ "$FAILED" -eq 0 ]
