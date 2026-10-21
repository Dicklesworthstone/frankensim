#!/usr/bin/env bash
# H4 test 1 (pexfp.104): re-observe the committed flute D4, G4 and A5 performances through a closed exterior BEM
# body of the published-Boehm flute, every opening radiating, with NO change to the source. The renders use the
# same command as ../render.sh plus --aperture-flows; tests/flue_aperture_flow_tap.rs proves the tapped WAV is
# bit-identical to the untapped one. The WAVs written here are therefore the existing (baffled Rayleigh) observer.
#
# Inputs in this directory:
#   body-coarse.tsv, body-fine.tsv  two panel-ladder rungs of the same body (4340 and 7840 panels); their
#                                   difference is the error bar on every comparison
#   receivers.tsv                   the 32 TU Berlin capsule directions at the SOFA distance 2.06 m and at the
#                                   paper's 1.05 m radius, through a declared pose of the flute
#   bind-<note>.tsv                 tapped network node -> body aperture (MOUTH is filled from the tap header)
# Band: 48 kHz kernels, pass band 10 kHz, C4 taper to 12 kHz, 64-sample latency, DFT ladder 256..2048.
#
# Run on trj:  SHA12=<sha12 holding these files> BIN_SHA12=<sha12 with music_render and reobserve_apertures> bash \
#                data/instruments/descriptions/flute/exterior/reobserve.sh
set -euo pipefail
SHA12=${SHA12:?set SHA12 to the pushed sha12 whose snapshot holds these files}
BIN_SHA12=${BIN_SHA12:-$SHA12}
BIN=/data/tmp/wa-inst/bin/$BIN_SHA12
SNAP=/data/tmp/wa-inst/snap/$SHA12/frankensim
OUT=${OUT:-/data/tmp/wa-inst/renders/flute-exterior/$SHA12}
THREADS=${THREADS:-8}
EXT=data/instruments/descriptions/flute/exterior

cd "$SNAP"
mkdir -p "$OUT"
for note in d4 g4 a5; do
  upper=$(printf '%s' "$note" | tr 'a-z' 'A-Z')
  if [ ! -e "$OUT/flute-$note.flows" ]; then
    nice -n 10 "$BIN/music_render" wind "data/instruments/descriptions/flute/flute-$note.performance" \
      "$OUT/Flute.nonvib.mf.$upper.wav" --aperture-flows "$OUT/flute-$note.flows" >"$OUT/render-$note.log" 2>&1 &
  fi
done
wait

BAND=(--fs 48000 --pass 10000 --stop 12000 --latency 64 --fft 256,512,1024,2048 --clearance 0.05)
for rung in coarse fine; do
  if [ ! -e "$OUT/$rung.kernels" ]; then
    nice -n 10 "$BIN/reobserve_apertures" bake "$EXT/body-$rung.tsv" "$EXT/receivers.tsv" "$OUT/$rung" \
      "${BAND[@]}" --threads "$THREADS" --samples "$OUT/$rung.samples" >"$OUT/bake-$rung.log" 2>&1
  fi
done

for note in d4 g4 a5; do
  mouth=$(head -c 4096 "$OUT/flute-$note.flows" | grep -a -m1 '^aperture mouth node' | cut -d' ' -f4)
  sed "s/^bind\tMOUTH\t/bind\t$mouth\t/" "$EXT/bind-$note.tsv" >"$OUT/bind-$note.tsv"
  for rung in coarse fine; do
    if [ ! -d "$OUT/$rung-$note" ]; then
      "$BIN/reobserve_apertures" replay "$OUT/$rung.kernels" "$OUT/flute-$note.flows" "$OUT/bind-$note.tsv" \
        "$OUT/$rung-$note" >"$OUT/replay-$rung-$note.log" 2>&1
    fi
  done
done
ls -l "$OUT"
