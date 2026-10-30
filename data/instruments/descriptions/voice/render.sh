#!/usr/bin/env bash
# Renders the voice description set with the generic `music_render voice` command.
# Usage (from the repository root): render.sh <path-to-music_render> <output-directory>
# Each file declares 3 s at 48 kHz mono; outputs are physically scaled, never normalized.
# story2008-u-Ds3, -u-C3, -u-C4 and their -wall/-flanagan variants run at the 96 kHz clock; they were refused before the DC-bin
# fix 2e629d5f3 (a defect in the mouth-flow line), which is why the -r1 files (48 kHz clock) exist. Use binary >= da31fbd67ee9.
# The -wall file carries an AUTHORED ESTIMATE wall (see its header). The -flanagan files carry the CITED yielding wall
# (Flanagan et al. via Birkholz, Haesner and Kuerbis, ICASSP 2022); each differs from its rigid twin only by that wall and its
# pitch factor q. C3, C4 and F4 are the notes VocalSet's singers sang (the only notes with a same-note reference).
# F4 exists for /a/ /e/ /o/ rigid and /a/ /o/ with the wall only (no stable F4 was found for the others).
# The -flanagan-r4 files are the same descriptions at a 192 kHz mechanics clock (ratio 4): the 96 kHz clock leaves the fold contact
# under-resolved (with an ideal-open mouth the /a/ C4 non-harmonic floor is converged by 192 kHz). The graded harmonics-to-noise barely
# moves (mean +0.6 dB over nine slots), so the clock is chosen for the converged fold contact, not for the grade. /a/ has none: the
# flanged mouth bound was lifted (ka 20 to 24) and /a/ was measured at 192 kHz, but C3 and C4 do not change (21.50 to 21.52 and
# 26.59 to 26.42 dB) and F4 gains 5 dB with its f0 outside the 3-cent trim, so /a/ stays at 96 kHz (status.json).
set -euo pipefail
BIN=${1:?path to the music_render binary}
OUT=${2:?output directory}
mkdir -p "$OUT"
for d in a-Ds3 e-Ds3 i-Ds3 o-Ds3 a-A3 a-D4 a-Ds3-wall u-Ds3 u-Ds3-r1 a-Ds3-r1 \
         a-C3 e-C3 i-C3 o-C3 u-C3 u-C3-r1 a-C4 e-C4 i-C4 o-C4 u-C4 u-C4-r1 a-F4 e-F4 o-F4 \
         a-C3-flanagan e-C3-flanagan i-C3-flanagan o-C3-flanagan u-C3-flanagan \
         a-C4-flanagan e-C4-flanagan i-C4-flanagan o-C4-flanagan u-C4-flanagan a-F4-flanagan o-F4-flanagan \
         e-C3-flanagan-r4 i-C3-flanagan-r4 o-C3-flanagan-r4 u-C3-flanagan-r4 \
         e-C4-flanagan-r4 i-C4-flanagan-r4 o-C4-flanagan-r4 u-C4-flanagan-r4 o-F4-flanagan-r4; do
  "$BIN" voice "data/instruments/descriptions/voice/story2008-$d.performance" "$OUT/story2008-$d.wav"
done
# /a/ at C3 with the spoken tract of the trained singer of Tom, Titze, Hoffman and Story (1999) (bead pexfp.325.2): only the
# tract differs from story2008-a-C3-flanagan. His sung (falsetto) and spoken tracts at C4 and F4 graded lower and are not kept.
"$BIN" voice data/instruments/descriptions/voice/tom1999-a-C3-speech-comf-flanagan.performance "$OUT/tom1999-a-C3-speech-comf-flanagan.wav"
