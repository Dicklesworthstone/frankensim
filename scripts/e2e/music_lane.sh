#!/usr/bin/env bash
# Music-program e2e lane engine (music bead frankensim-music-v8-root-3ez8g.16).
#
# Renders ONE pinned fixture through the committed music_render binary,
# proves same-host replay (two renders, byte-identical WAV + sidecar),
# checks the never-overwrite refusal arm, and emits one JSON-lines
# receipt per stage plus a terminal receipt chaining them. Artifacts
# (WAV + provenance sidecar) are the retained evidence.
#
# This validates software determinism and encoding law only. It is not
# a physical-validation, realism, or listening claim; those live in the
# claims registry's gates and listening receipts.
#
# Baseline discipline: WAV hashes are ONE-HOST bit-deterministic (the
# music stack claims no cross-ISA replay; zero cross-ISA goldens exist).
# The committed baseline in scripts/e2e/music-baselines.json is therefore
# compared as a NOTE by default (verdict matched|drifted|no-baseline) and
# enforced only under --enforce-baseline (the recording host / CI lane).
# A drift under enforcement is a golden event: explain the semantic cause
# and update the baseline in the same commit, never loosen the check.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

readonly BASELINES="scripts/e2e/music-baselines.json"

FIXTURE=""
OUT_DIR=""
BINARY=""
SECONDS_ARG="0.25"
ENFORCE_BASELINE=0
DRY_RUN=0

usage() {
  printf '%s\n' \
    'usage: scripts/e2e/music_lane.sh --fixture <reed|string> --out-dir PATH [options]' \
    '' \
    'Options:' \
    '  --binary PATH        music_render executable (default: auto-discover under' \
    '                       $CARGO_TARGET_DIR or target/, debug then release).' \
    '  --seconds S          Render length in seconds (default 0.25).' \
    '  --enforce-baseline   Fail (not note) on baseline hash drift for this host lane.' \
    '  --dry-run            Validate inputs, print the plan, create nothing.' \
    '  -h|--help            This help.' \
    '' \
    'The out dir must not exist; the lane creates it and refuses to overwrite' \
    'evidence. Receipts go to <out-dir>/receipts.jsonl.'
}

die() {
  printf '{"suite":"music-e2e","lane":"%s","verdict":"refused","elapsed_s":%s,"reason":"%s"}\n' \
    "${FIXTURE:-unset}" "$SECONDS" "$*" >&2
  exit 2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fixture) [[ $# -ge 2 ]] || die '--fixture requires reed|string'; FIXTURE="$2"; shift 2 ;;
    --out-dir) [[ $# -ge 2 ]] || die '--out-dir requires a path'; OUT_DIR="$2"; shift 2 ;;
    --binary)  [[ $# -ge 2 ]] || die '--binary requires a path'; BINARY="$2"; shift 2 ;;
    --seconds) [[ $# -ge 2 ]] || die '--seconds requires a number'; SECONDS_ARG="$2"; shift 2 ;;
    --enforce-baseline) ENFORCE_BASELINE=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[[ "$FIXTURE" == "reed" || "$FIXTURE" == "string" ]] || die '--fixture must be reed or string'
[[ -n "$OUT_DIR" ]] || die '--out-dir is required'
[[ ! -e "$OUT_DIR" ]] || die "refusing to overwrite existing out dir: $OUT_DIR"

if [[ -z "$BINARY" ]]; then
  for candidate in \
    "${CARGO_TARGET_DIR:-target}/debug/music_render" \
    "${CARGO_TARGET_DIR:-target}/release/music_render" \
    "target/debug/music_render" "target/release/music_render"; do
    if [[ -x "$candidate" ]]; then BINARY="$candidate"; break; fi
  done
fi
[[ -n "$BINARY" && -x "$BINARY" ]] || \
  die 'music_render binary not found; build it (cargo build -p fs-couple --bin music_render) or pass --binary'

if [[ $DRY_RUN -eq 1 ]]; then
  printf '{"suite":"music-e2e","lane":"%s","verdict":"dry-run","binary":"%s","seconds":%s,"out_dir":"%s"}\n' \
    "$FIXTURE" "$BINARY" "$SECONDS_ARG" "$OUT_DIR"
  exit 0
fi

mkdir -p "$OUT_DIR"
RECEIPTS="$OUT_DIR/receipts.jsonl"
: > "$RECEIPTS"

receipt() { printf '%s\n' "$1" | tee -a "$RECEIPTS"; }

hash_field() { # extract wav_blake3 from a provenance sidecar
  sed -n 's/.*"wav_blake3":"\([0-9a-f]*\)".*/\1/p' "$1"
}

# Stage 1+2: two independent renders (replay proof).
for run in a b; do
  WAV="$OUT_DIR/$FIXTURE-$run.wav"
  if ! "$BINARY" "$FIXTURE" "$WAV" --seconds "$SECONDS_ARG" > "$OUT_DIR/render-$run.stdout"; then
    die "render $run refused; see $OUT_DIR/render-$run.stdout"
  fi
  [[ -s "$WAV" ]] || die "render $run produced no WAV"
  [[ -s "$OUT_DIR/$FIXTURE-$run.provenance.json" ]] || die "render $run produced no sidecar"
done
HASH_A="$(hash_field "$OUT_DIR/$FIXTURE-a.provenance.json")"
HASH_B="$(hash_field "$OUT_DIR/$FIXTURE-b.provenance.json")"
[[ -n "$HASH_A" && "$HASH_A" == "$HASH_B" ]] || \
  die "replay hash mismatch: $HASH_A vs $HASH_B (same-host determinism regression)"
cmp -s "$OUT_DIR/$FIXTURE-a.wav" "$OUT_DIR/$FIXTURE-b.wav" || \
  die 'replay byte mismatch: WAV files differ despite matching sidecars'
receipt "{\"stage\":\"replay\",\"lane\":\"$FIXTURE\",\"verdict\":\"bit-identical\",\"wav_blake3\":\"$HASH_A\",\"seconds\":$SECONDS_ARG}"

# Stage 3: the never-overwrite refusal arm must fire.
if "$BINARY" "$FIXTURE" "$OUT_DIR/$FIXTURE-a.wav" --seconds "$SECONDS_ARG" \
    > "$OUT_DIR/refusal.stdout" 2>&1; then
  die 'overwrite refusal arm did not fire; the lane cannot trust the artifacts'
fi
grep -q 'refuses to overwrite evidence' "$OUT_DIR/refusal.stdout" || \
  die 'overwrite refusal fired with the wrong reason'
receipt "{\"stage\":\"refusal-arm\",\"lane\":\"$FIXTURE\",\"verdict\":\"refused-as-required\"}"

# Stage 4: baseline comparison (note by default; gate under enforcement).
BASELINE_VERDICT="no-baseline"
if [[ -f "$BASELINES" ]]; then
  RECORDED="$(sed -n "s/.*\"$FIXTURE\": *\"\([0-9a-f]*\)\".*/\1/p" "$BASELINES" | head -1)"
  if [[ -n "$RECORDED" ]]; then
    if [[ "$RECORDED" == "$HASH_A" ]]; then BASELINE_VERDICT="matched"; else BASELINE_VERDICT="drifted"; fi
  fi
fi
receipt "{\"stage\":\"baseline\",\"lane\":\"$FIXTURE\",\"verdict\":\"$BASELINE_VERDICT\",\"recorded\":\"${RECORDED:-}\",\"observed\":\"$HASH_A\",\"determinism_class\":\"one-host\"}"
if [[ $ENFORCE_BASELINE -eq 1 && "$BASELINE_VERDICT" == "drifted" ]]; then
  die "baseline drift under enforcement: recorded ${RECORDED:-} vs observed $HASH_A (golden event: explain and update in the same commit)"
fi

# Terminal receipt: chain the stages.
receipt "{\"stage\":\"terminal\",\"lane\":\"$FIXTURE\",\"verdict\":\"green\",\"stages\":[\"replay\",\"refusal-arm\",\"baseline\"],\"wav_blake3\":\"$HASH_A\",\"baseline\":\"$BASELINE_VERDICT\",\"elapsed_s\":$SECONDS}"
