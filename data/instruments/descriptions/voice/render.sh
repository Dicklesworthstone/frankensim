#!/usr/bin/env bash
# Renders the voice description set with the generic `music_render voice` command.
# Usage (from the repository root): render.sh <path-to-music_render> <output-directory>
# Each file declares 3 s at 48 kHz mono; outputs are physically scaled, never normalized.
# story2008-u-Ds3 and story2008-u-Ds3-wall (96 kHz clock) are deliberately absent: their mouth-flow ladder does
# not converge within 65536 taps and the command refuses them. The -r1 files use a 48 kHz physical clock, so the
# same cap spans twice the response; a-Ds3-r1 against a-Ds3 is a first mechanics-rate comparison. The -wall file
# carries an AUTHORED ESTIMATE wall (see its header) to measure what the wall alone changes against a-Ds3.
set -euo pipefail
BIN=${1:?path to the music_render binary}
OUT=${2:?output directory}
mkdir -p "$OUT"
for d in a-Ds3 e-Ds3 i-Ds3 o-Ds3 a-A3 a-D4 a-Ds3-wall u-Ds3-r1 a-Ds3-r1; do
  "$BIN" voice "data/instruments/descriptions/voice/story2008-$d.performance" "$OUT/story2008-$d.wav"
done
