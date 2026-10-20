#!/usr/bin/env bash
# Violin open-string arco notes: 3.0 s, 48 kHz, mono, one WAV per string, bowed at constant
# speed and force with a constant-acceleration onset (inputs and sources in sources.tsv,
# citations.tsv; what is missing in unknowns.tsv).
# usage: render.sh MUSIC_RENDER_BIN OUT_DIR
# Matching UIowa MIS reference (the first note of each file is the open string):
#   G3  Violin.arco.mf.sulG.G3B3
#   D4  Violin.arco.mf.sulD.D4A4
#   A4  Violin.arco.mf.sulA.A4B4
#   E5  Violin.arco.mf.sulE.E5B5
set -euo pipefail
bin="${1:?usage: render.sh MUSIC_RENDER_BIN OUT_DIR}"
out="${2:?usage: render.sh MUSIC_RENDER_BIN OUT_DIR}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
mkdir -p "$out"
for note in G3 D4 A4 E5; do
    "$bin" ensemble "$out/violin-arco-mf-open-$note.wav" \
        --bow "$here/arco-open-$note.performance" --full-scale-pa 1 --block 4096
done
