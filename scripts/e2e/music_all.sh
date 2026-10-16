#!/usr/bin/env bash
# Program-level music e2e umbrella (music bead frankensim-music-v8-root-3ez8g.16):
# chains every AVAILABLE track lane and emits one program receipt. Tracks
# whose fixtures are not yet gated are reported BY NAME as not-yet-gated —
# the NO-DATA-vs-skipped discipline: an ungated track is a named absence,
# never a silent omission that reads as coverage.
#
# usage: scripts/e2e/music_all.sh --out-dir PATH [--binary PATH] [--enforce-baseline]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"
cd "$REPO_ROOT"

OUT_DIR=""
EXTRA_ARGS=()

die() {
  printf '{"suite":"music-e2e","lane":"all","verdict":"refused","reason":"%s"}\n' "$*" >&2
  exit 2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) [[ $# -ge 2 ]] || die '--out-dir requires a path'; OUT_DIR="$2"; shift 2 ;;
    --binary|--seconds) [[ $# -ge 2 ]] || die "$1 requires a value"; EXTRA_ARGS+=("$1" "$2"); shift 2 ;;
    --enforce-baseline) EXTRA_ARGS+=("$1"); shift ;;
    -h|--help)
      printf 'usage: scripts/e2e/music_all.sh --out-dir PATH [--binary PATH] [--seconds S] [--enforce-baseline]\n'
      exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done
[[ -n "$OUT_DIR" ]] || die '--out-dir is required'
[[ ! -e "$OUT_DIR" ]] || die "refusing to overwrite existing out dir: $OUT_DIR"
mkdir -p "$OUT_DIR"

# Available lanes: track name -> lane script. Extend as tracks gate.
readonly -a AVAILABLE=("wind:music_wind.sh" "string:music_string.sh")
# Not-yet-gated tracks, BY NAME (from the v8 program graph). Each entry
# names the bead that gates it, so the absence is actionable.
readonly -a UNGATED=(
  "brass:3ez8g.4.4"
  "piano:3ez8g.5.3"
  "bowed:3ez8g.7.5"
  "voice:3ez8g.8.3"
  "electric:3ez8g.9.5"
  "jet:3ez8g.10.4"
  "bar:3ez8g.12.1"
  "guitar-body:3ez8g.7.4"
)

PROGRAM_RECEIPT="$OUT_DIR/program-receipt.jsonl"
: > "$PROGRAM_RECEIPT"
GREEN=0
RED=0
for entry in "${AVAILABLE[@]}"; do
  track="${entry%%:*}"
  lane="${entry#*:}"
  if "$HERE/$lane" --out-dir "$OUT_DIR/$track" "${EXTRA_ARGS[@]}" > "$OUT_DIR/$track-lane.stdout" 2>&1; then
    verdict="green"; GREEN=$((GREEN + 1))
  else
    verdict="red"; RED=$((RED + 1))
  fi
  printf '{"track":"%s","verdict":"%s","receipts":"%s/receipts.jsonl"}\n' \
    "$track" "$verdict" "$OUT_DIR/$track" | tee -a "$PROGRAM_RECEIPT"
done
for entry in "${UNGATED[@]}"; do
  track="${entry%%:*}"
  bead="${entry#*:}"
  printf '{"track":"%s","verdict":"not-yet-gated","gating_bead":"%s"}\n' \
    "$track" "$bead" | tee -a "$PROGRAM_RECEIPT"
done
printf '{"suite":"music-e2e","lane":"all","green":%d,"red":%d,"not_yet_gated":%d,"verdict":"%s"}\n' \
  "$GREEN" "$RED" "${#UNGATED[@]}" "$([[ $RED -eq 0 ]] && echo green || echo red)" \
  | tee -a "$PROGRAM_RECEIPT"
[[ $RED -eq 0 ]] || exit 1
