#!/usr/bin/env bash
# Cello open-string arco notes: 3.0 s, 48 kHz, mono, one WAV per string, bowed at constant
# speed and force with a constant-acceleration onset (inputs and sources in sources.tsv,
# citations.tsv; what is missing in unknowns.tsv). A labelled diagnostic baseline: one-string
# point bow with a one-way one-oscillator body, not the instrument.
# PCM full scale 0.04 Pa is chosen only to use the 16-bit range (probe peaks were 1e-3 to 5e-2 Pa).
# usage: render.sh MUSIC_RENDER_BIN OUT_DIR
# Matching UIowa MIS reference (the first note of each file is the open string):
#   C2  Cello.arco.mf.sulC.C2Ab2
#   G2  Cello.arco.mf.sulG.G2Bb2
#   D3  Cello.arco.mf.sulD.D3B3
#   A3  Cello.arco.mf.sulA.A3B3
set -euo pipefail
bin="${1:?usage: render.sh MUSIC_RENDER_BIN OUT_DIR}"
out="${2:?usage: render.sh MUSIC_RENDER_BIN OUT_DIR}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
mkdir -p "$out"
for note in C2 G2 D3 A3; do
    "$bin" ensemble "$out/cello-arco-mf-open-$note.wav" \
        --bow "$here/arco-open-$note.performance" --full-scale-pa 0.04 --block 4096
done
