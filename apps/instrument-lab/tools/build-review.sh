#!/usr/bin/env bash
# build-review.sh <analysis-dir> <out-dir> <commit>
#
# The instruments-only review build: instruments.html (served as index.html too) and only the data, figures and audio
# tools/build_instruments.py writes for it. A note the grader refused has no files, so none ships. Output goes to a NEW
# directory (an existing one is never touched); the build fails rather than exceed what the private page host takes.
# The host serves no .m4a under any content type, so the clips are written as .mp4 (the same AAC audio; only the file
# name differs) and the cards, the content-hash map and the files all say so.
#
#   cd apps/instrument-lab && bash tools/build-review.sh ~/.local/share/wa-inst/analysis ~/.local/share/wa-inst/review-<sha> <sha>
set -euo pipefail
cd "$(dirname "$0")/.."
ANALYSIS="${1:?usage: build-review.sh <analysis-dir> <out-dir> <commit>}"
OUT="${2:?usage: build-review.sh <analysis-dir> <out-dir> <commit>}"
SHA="${3:?usage: build-review.sh <analysis-dir> <out-dir> <commit>}"
MAX_FILES="${REVIEW_MAX_FILES:-500}"   # the host's cap is 511
MAX_MB="${REVIEW_MAX_MB:-60}"
# Voice has 13 graded notes at 8 files each, more than the host's file limit leaves room for: show one per vowel at C3
# (a rendered note, C3 before C4). The card says how many notes it left out. REVIEW_ONLY="family=note,note family=note"
# changes the choice and REVIEW_ONLY="" shows every note.
# A family with no render shows none of its recording-only notes (the card lists them as left out), so the files go to what
# has been rendered; the flag follows the data, so a family that gets a render shows its recordings again.
ONLY="${REVIEW_ONLY-voice=C3-a,C3-e,C3-i,C3-o,C3-u}"
only_args=()
for spec in $ONLY; do only_args+=(--only "$spec"); done

[ ! -e "$OUT" ] || { echo "refusing: $OUT already exists (a review build is never overwritten; give it a new name)" >&2; exit 64; }
[ -d "$ANALYSIS" ] || { echo "no analysis directory $ANALYSIS" >&2; exit 64; }

STAGE="$(mktemp -d "${TMPDIR:-/tmp}/review-public.XXXXXX")"
echo "staging the data in $STAGE"
status=0
uv run --with numpy --with pillow python tools/build_instruments.py --analysis "$ANALYSIS" --out "$STAGE/data/instruments" --audio-ext mp4 --recordings-only-if-rendered ${only_args[@]+"${only_args[@]}"} || status=$?
# 2 = the generator itself rejected a pair; the card lists it, the review still builds. Anything else is a real failure.
[ "$status" -eq 0 ] || [ "$status" -eq 2 ] || { echo "build_instruments.py failed ($status)" >&2; exit "$status"; }
# The reference-candidates section, when a candidates tree is named (REVIEW_CANDIDATES=<mirror of trj:/data/tmp/wa-inst/refs-candidates>).
# A record that fails a check is listed on its family's page and the build goes on (exit 2), as for a pair.
if [ -n "${REVIEW_CANDIDATES:-}" ]; then
  [ -d "$REVIEW_CANDIDATES" ] || { echo "no candidates directory $REVIEW_CANDIDATES" >&2; exit 64; }
  cstatus=0
  # REVIEW_CAND_MAX: how many candidates of a family to show (the best ranked; the rest are listed as omitted), to fit the host's file limit.
  uv run --with numpy python tools/build_candidates.py --candidates "$REVIEW_CANDIDATES" --out "$STAGE/data/candidates" --max-per-family "${REVIEW_CAND_MAX:-5}" || cstatus=$?
  [ "$cstatus" -eq 0 ] || [ "$cstatus" -eq 2 ] || { echo "build_candidates.py failed ($cstatus)" >&2; exit "$cstatus"; }
  [ "$status" -ne 0 ] || status="$cstatus"
fi

./node_modules/.bin/tsc --noEmit
REVIEW_PUBLIC="$STAGE" REVIEW_OUT="$OUT" VITE_REVIEW_SHA="$SHA" ./node_modules/.bin/vite build --config vite.instruments.config.ts

# What shipped: only the page, its bundle and its own data.
stray="$(cd "$OUT" && find . -type f ! -name instruments.html ! -name index.html ! -path './assets/*' ! -path './data/instruments/*' ! -path './data/candidates/*' | head -5)"
[ -z "$stray" ] || { echo "unexpected files in the review build:" >&2; echo "$stray" >&2; exit 1; }
m4a="$(cd "$OUT" && find . -type f -name '*.m4a' | head -3)"
[ -z "$m4a" ] || { echo "the host serves no .m4a, but the review build holds:" >&2; echo "$m4a" >&2; exit 1; }
FILES="$(find "$OUT" -type f | wc -l | tr -d ' ')"
BYTES="$(find "$OUT" -type f -exec stat -f %z {} + | awk '{s+=$1} END {print s+0}')"
MB="$(awk -v b="$BYTES" 'BEGIN {printf "%.1f", b/1e6}')"
echo "REVIEW $OUT  files=$FILES  size_mb=$MB  commit=$SHA  generator_exit=$status"
if [ "$FILES" -gt "$MAX_FILES" ] || awk -v b="$BYTES" -v m="$MAX_MB" 'BEGIN {exit !(b/1e6 > m)}'; then
  echo "TOO BIG for the host: $FILES files (cap $MAX_FILES), $MB MB (cap $MAX_MB)" >&2
  exit 1
fi
