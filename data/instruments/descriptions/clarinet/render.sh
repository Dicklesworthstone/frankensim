#!/usr/bin/env bash
# Native renders of the committed Buffet-Prestige/Taillard reed performances (pexfp.228).
# Run on trj against a pushed SHA whose release music_render exists:
#   SHA12=<sha12> data/instruments/descriptions/clarinet/render.sh
# Each note is mf at 48 kHz mono, 2.5 s, observed at the physical baffled_outlet receiver.
# References (grading only, never an input): UIowa MIS BbClar.mf.D3B3 / C4B4 / C5B5 aiff, concert pitch.
set -euo pipefail

SHA12=${SHA12:?set SHA12 to the pushed sha12 that was built}
BIN=${BIN:-/data/tmp/wa-inst/bin/$SHA12/music_render}
SNAP=${SNAP:-/data/tmp/wa-inst/snap/$SHA12/frankensim}
OUT=${OUT:-/data/tmp/wa-inst/renders/clarinet/$SHA12}
NOTES=${NOTES:-"D3 G3 Bb3 D4 G4 D5 Bb5"}
PERF=data/instruments/descriptions/clarinet/performances

mkdir -p "$OUT"
cd "$SNAP"
for note in $NOTES; do
  nohup nice -n 10 "$BIN" wind "$PERF/BbClar.mf.$note.performance" "$OUT/BbClar.mf.$note.wav" \
    >"$OUT/BbClar.mf.$note.log" 2>&1 &
done
wait
ls -l "$OUT"/*.wav
