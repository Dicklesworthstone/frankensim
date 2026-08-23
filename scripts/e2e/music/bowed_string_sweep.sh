#!/usr/bin/env bash
# Bowed-string playable-region sweep + listening receipt
# (bead frankensim-music-v8-root-3ez8g.7.5).
#
# Runs the Schelleng-style (F_n, v_bow) regime sweep through the release
# build and renders the open-string stroke listening receipt. Both artifacts
# are deterministic on one host; the sweep JSONL is the retained receipt.
#
# Usage: scripts/e2e/music/bowed_string_sweep.sh [OUT_DIR]
set -u -o pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
OUT_DIR="${1:-$REPO_ROOT/.e2e-out/bowed-sweep}"
mkdir -p "$OUT_DIR"

cargo run -r -q -p fs-couple --bin bowed_sweep -- \
  --out "$OUT_DIR/sweep.jsonl" || exit $?

cargo run -r -q -p fs-couple --bin bowed_listen -- \
  --wav "$OUT_DIR/open_string_stroke.wav" \
  --receipt "$OUT_DIR/listening-receipt.jsonl" || exit $?

echo "artifacts: $OUT_DIR/sweep.jsonl $OUT_DIR/open_string_stroke.wav $OUT_DIR/listening-receipt.jsonl"
