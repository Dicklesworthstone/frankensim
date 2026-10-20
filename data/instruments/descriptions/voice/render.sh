#!/usr/bin/env bash
# Renders the voice description set with the generic `music_render voice` command.
# Usage (from the repository root): render.sh <path-to-music_render> <output-directory>
# Each file declares 3 s at 48 kHz mono; outputs are physically scaled, never normalized.
# story2008-u-Ds3 (rigid wall) is deliberately absent: its mouth-flow ladder does not converge and the
# command refuses it. The -wall files carry an AUTHORED ESTIMATE wall (see their headers); a-Ds3-wall
# exists to measure what the wall alone changes against a-Ds3.
set -euo pipefail
BIN=${1:?path to the music_render binary}
OUT=${2:?output directory}
mkdir -p "$OUT"
for d in a-Ds3 e-Ds3 i-Ds3 o-Ds3 a-A3 a-D4 u-Ds3-wall a-Ds3-wall; do
  "$BIN" voice "data/instruments/descriptions/voice/story2008-$d.performance" "$OUT/story2008-$d.wav"
done
