#!/usr/bin/env bash
# H4 test 1 (pexfp.104): re-observe the committed flute D4, G4 and A5 performances through a closed exterior BEM
# body of the published-Boehm flute, every opening radiating, with NO change to the source. The renders use the
# same command as ../render.sh plus --aperture-flows; tests/flue_aperture_flow_tap.rs proves the tapped WAV is
# bit-identical to the untapped one. The WAVs written here are therefore the existing (baffled Rayleigh) observer.
#
# Inputs in this directory:
#   body-rung-a.tsv, body-rung-b.tsv  the panel ladder used (3024 and 2140 panels; the dense solver measured 585 s
#                                     per frequency at 2140 panels on trj); their difference is the error bar
#   body-coarse.tsv, body-fine.tsv    4340 and 7840 panels, kept for a faster solver route (FIR=1 below)
#   receivers.tsv                     the 32 TU Berlin capsule directions at the SOFA distance 2.06 m and at the
#                                     paper's 1.05 m radius, through a declared pose of the flute
#   bind-<note>.tsv                   tapped network node -> body aperture (MOUTH is filled from the tap header)
#   measured-partials.tsv             the measured per-capsule partial levels (facts only)
#   compare_h4.py                     patterns of old observer, BEM rungs and measurement; error bars
# Steady window 1.0-2.2 s; partials up to 10 kHz.
#
# Run on trj:  SHA12=<sha12 holding these files> BIN_SHA12=<sha12 with music_render and reobserve_apertures> bash \
#                data/instruments/descriptions/flute/exterior/reobserve.sh
set -euo pipefail
SHA12=${SHA12:?set SHA12 to the pushed sha12 whose snapshot holds these files}
BIN_SHA12=${BIN_SHA12:-$SHA12}
BIN=/data/tmp/wa-inst/bin/$BIN_SHA12
SNAP=/data/tmp/wa-inst/snap/$SHA12/frankensim
OUT=${OUT:-/data/tmp/wa-inst/renders/flute-exterior/$SHA12}
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
for note in d4 g4 a5; do
  mouth=$(head -c 4096 "$OUT/flute-$note.flows" | grep -a -m1 '^aperture mouth node' | cut -d' ' -f4)
  sed "s/^bind\tMOUTH\t/bind\t$mouth\t/" "$EXT/bind-$note.tsv" >"$OUT/bind-$note.tsv"
done

# Steady-state partials through the body, both rungs; about 24 solver threads in total.
declare -A F0=([d4]=293.665 [g4]=391.995 [a5]=880.0)
declare -A THREADS=([d4]=8 [g4]=6 [a5]=4)
harmonics() { # rung note threads
  local tsv="$OUT/harmonics-$1-$2.tsv"
  [ -e "$tsv" ] || nice -n 10 "$BIN/reobserve_apertures" harmonics "$EXT/body-$1.tsv" "$EXT/receivers.tsv" \
    "$OUT/flute-$2.flows" "$OUT/bind-$2.tsv" "$tsv" --f0-nominal "${F0[$2]}" --band 10000 --window 1.0,2.2 \
    --clearance 0.05 --threads "$3" >"$OUT/harmonics-$1-$2.log" 2>&1
}
for note in d4 g4 a5; do harmonics rung-a "$note" "${THREADS[$note]}" & done
(for note in d4 g4 a5; do harmonics rung-b "$note" 6; done) &
wait

uv run --python 3.12 --with numpy python "$EXT/compare_h4.py" "$OUT" "$EXT/measured-partials.tsv" \
  "$EXT/body-rung-a.tsv" rung-a rung-b "$OUT/h4-report.json" | tee "$OUT/h4-report.txt"

# Optional time-domain re-observation (listening pairs): FIR bake and replay; days on the dense owner today.
if [ "${FIR:-0}" = 1 ]; then
  BAND=(--fs 48000 --pass 10000 --stop 12000 --latency 64 --fft 256,512,1024,2048 --clearance 0.05)
  BAND+=(--only embouchure,foot,hole-07,hole-08,hole-09,hole-10,hole-11,hole-12,hole-13,hole-14,hole-15,hole-16)
  for rung in coarse fine; do
    [ -e "$OUT/$rung.kernels" ] || nice -n 10 "$BIN/reobserve_apertures" bake "$EXT/body-$rung.tsv" \
      "$EXT/receivers.tsv" "$OUT/$rung" "${BAND[@]}" --threads 8 --samples "$OUT/$rung.samples" >"$OUT/bake-$rung.log" 2>&1
    for note in d4 g4 a5; do
      [ -d "$OUT/$rung-$note" ] || "$BIN/reobserve_apertures" replay "$OUT/$rung.kernels" "$OUT/flute-$note.flows" \
        "$OUT/bind-$note.tsv" "$OUT/$rung-$note" >"$OUT/replay-$rung-$note.log" 2>&1
    done
  done
fi
ls -l "$OUT"
