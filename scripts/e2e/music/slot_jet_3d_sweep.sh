#!/usr/bin/env bash
# T-Jet 3-D slot-jet Reynolds sweep lane (music bead
# frankensim-music-v8-root-3ez8g.10.1): the recorded heavy-run recipe.
#
# Runs the committed slot_jet_3d_sweep binary through RCH in release
# mode (RCH_MIN_LOCAL_TIME_MS pins execution away from a quick local
# fallback; CARGO_TARGET_DIR lands under RCH_TARGET_BASE, resolving on
# the execution host). One geometry per invocation, a comma list of
# second-order rates as the Re ladder; every rung settles/records/
# classifies independently with typed refusals.
#
# Receipts: OUT/run.jsonl (binary-owned, fail-closed against reruns;
# the same lines stream to stdout) plus the stderr log. The binary
# ALSO receives --out so its internal refusal guard and the script's
# guard protect the SAME path.
#
# Box-sensitivity discipline (one spanwise octave) is a SECOND
# invocation at --nz 64 restricted to the rung(s) under scrutiny;
# compare the two terminal verdicts, not a single number.
#
# No-claim boundary: lattice measurements only — no experimental,
# video-backed, or absolute-level flue-noise claim survives this lane.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO_ROOT"

usage() {
  printf '%s\n' \
    'usage: scripts/e2e/music/slot_jet_3d_sweep.sh --out DIR [--budget N] [--max-chunks M] [--r2 1.95,1.97,...] [sweep flags]' >&2
}

OUT=""
BUDGET=0
MAX_CHUNKS=40
ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out)
      [[ $# -ge 2 ]] || { usage; exit 2; }
      OUT="$2"; shift 2 ;;
    --budget)
      [[ $# -ge 2 ]] || { usage; exit 2; }
      BUDGET="$2"; shift 2 ;;
    --max-chunks)
      [[ $# -ge 2 ]] || { usage; exit 2; }
      MAX_CHUNKS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) ARGS+=("$1"); shift ;;
  esac
done
[[ -n "$OUT" ]] || { usage; exit 2; }
if [[ -e "$OUT" ]]; then
  echo "refusing existing output path $OUT (fail-closed rerun)" >&2
  exit 3
fi
mkdir -p "$OUT"

if ((BUDGET > 0)); then
  # Budget (checkpoint/resume) lane: repeatedly invoke the binary
  # with a per-invocation step budget until the terminal receipt
  # appears or the chunk cap fires. The driver owns run.jsonl; every
  # invocation's lines append verbatim and its full stdout is kept
  # as a chunk log.
  HEADER="{\"schema\":\"fs-aeroac.slot-jet-3d.sweep/v1\",\"mode\":\"budget\",\"steps_budget\":${BUDGET},\"args\":\"$(printf '%s ' "${ARGS[@]}" | sed 's/"/\\"/g')\"}"
  printf '%s\n' "$HEADER" >>"$OUT/run.jsonl"
  for ((chunk = 1; chunk <= MAX_CHUNKS; chunk++)); do
    CHUNK_LOG="$OUT/chunk-$(printf '%04d' "$chunk").log"
    export RCH_MIN_LOCAL_TIME_MS=9999999
    if ! rch exec -- cargo run --release -p fs-aeroac --bin slot_jet_3d_sweep -- \
      --out "$OUT" --steps-budget "$BUDGET" "${ARGS[@]}" | tee "$CHUNK_LOG"; then
      echo "chunk $chunk invocation failed" >&2
      exit 4
    fi
    grep -E '^\{' "$CHUNK_LOG" >>"$OUT/run.jsonl"
    if grep -q 'slot-jet-3d.terminal/v1' "$CHUNK_LOG"; then
      echo "terminal receipt at chunk $chunk"
      break
    fi
  done
  echo "receipts: $OUT/run.jsonl"
  grep -c . "$OUT/run.jsonl"
  exit 0
fi


ARG_STRING="--out $(printf '%q' "$OUT")"
if ((${#ARGS[@]} > 0)); then
  ARG_STRING+=" $(printf '%q ' "${ARGS[@]}")"
fi

export RCH_MIN_LOCAL_TIME_MS=9999999
# Single-quoted remote body: RCH_TARGET_BASE resolves on the
# execution host, never on this one.
rch exec -- sh -c '
  set -eu
  export CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/cargo-slot-jet-3d"
  mkdir -p "$CARGO_TARGET_DIR"
  cargo run --release -p fs-aeroac --bin slot_jet_3d_sweep -- '"$ARG_STRING"'
' | tee "$OUT/stdout.log"

echo "receipts: $OUT/run.jsonl"
grep -c . "$OUT/run.jsonl"
