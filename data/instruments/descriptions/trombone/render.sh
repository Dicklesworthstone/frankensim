#!/usr/bin/env bash
# Native lip-valve renders of the tenor trombone description (pexfp.254, lane P5).
# Run on trj. The binary must contain `music_render lip` (the brass_loop::performance owner):
#   SHA12=<sha12 whose snapshot holds these files> BIN_SHA12=<sha12 of a music_render with `lip`> \
#     data/instruments/descriptions/trombone/render.sh
#
# Composition: tenor-trombone.head (every bore, lip, cup and microphone number, each with its source)
# plus a canonical gesture schedule written below: the player's mouth pressure and lip pre-stress.
# Slide in first position, mf. There is no pitch, note or frequency anywhere in an input: the note is
# the lock the loop finds between the lip and the bore.
#
#  * sweep-pNNN: lip pre-stress ramps 0.5 -> 2.0 at a held mouth pressure. A diagnostic, not a note:
#    it records which regime the passive bore + lip speaks as the lip tightens (the one-dimensional
#    cut through the blowing-pressure x lip-tension threshold map).
#  * hold-tNNN-pNNN: the same lip held at a fixed pre-stress for 3 s. The player's embouchure setting
#    is the free player input; which setting yields which note is read from the render afterwards.
#
# Mouth pressure 1.5 kPa: mf on trombone in the artificial-mouth and player literature (Boutin et al.
# 2015, JASA 137:1195, B-flat2 dc mouth pressure 1.4 kPa +-10 %; Bromage 2007 PhD Edinburgh Table 4.3,
# E3 mf 14-16 mbar). Per-note and per-dynamic tables for trombone are UNKNOWN.
#
# Clock: 19.2 kHz mechanics = output, not 48 kHz. The fs-duct flanged radiation load is admitted only
# to ka = 20, and the 0.1086 m bell reaches ka = 20 at 10.05 kHz, so a 19.2 kHz clock (Nyquist 9.6 kHz,
# ka 19.1) is the admitted choice; 48 kHz mechanics would be refused. The radiation receiver's resolved
# band is 3.9 kHz (see the microphone record). The WAV is written at 19200 Hz and never resampled;
# grade up to 3.9 kHz.
#
# PIECES=<n> re-renders with n cones per Bessel element (default 16): the convergence ladder on the
# horn expansion, written to $OUT/pieces-<n>.
#
# References (grading only, never an input): UIowa MIS TenorTrombone.mf.* range files, ids
# uiowa-tenortrombone-tenortrombone-mf-e2b2 / -c3b3 / -c4b4 / -c5 in references.json.
set -euo pipefail

SHA12=${SHA12:?set SHA12 to the pushed sha12 whose snapshot holds these files}
BIN_SHA12=${BIN_SHA12:-$SHA12}
BIN=${BIN:-/data/tmp/wa-inst/bin/$BIN_SHA12/music_render}
SNAP=${SNAP:-/data/tmp/wa-inst/snap/$SHA12/frankensim}
OUT=${OUT:-/data/tmp/wa-inst/renders/trombone/$SHA12}
PIECES=${PIECES:-16}
PRESSURE=${PRESSURE:-1.5e3}
DIR=data/instruments/descriptions/trombone
RATE=19200
HOLD_TENSIONS=${HOLD_TENSIONS:-"6e-1 8e-1 1e0 1.2e0 1.4e0 1.6e0 1.8e0"}

cd "$SNAP"
if [ "$PIECES" != 16 ]; then
  OUT="$OUT/pieces-$PIECES"
fi
mkdir -p "$OUT"

# compose <file> <seconds> <mouth-pressure-pa> <initial-pre-stress> [<final-pre-stress>]
compose() {
  local file=$1 seconds=$2 pressure=$3 initial=$4 final=${5:-}
  {
    printf 'frankensim-lip-performance-v1\naudio %d %d\n' "$RATE" $((RATE * seconds))
    sed -E "s/^(bessel .* )16\$/\1$PIECES/" "$DIR/tenor-trombone.head"
    printf 'schedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t2\n'
    printf 'track\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t1\n'
    printf 'event\t0e0\t5e-2\tpressure-pa\t%s\n' "$pressure"
    printf 'track\tlip\tpre-stress\ninitial\tmultiplier\t%s\n' "$initial"
    if [ -n "$final" ]; then
      printf 'events\t1\nevent\t0e0\t%se0\tmultiplier\t%s\n' "$seconds" "$final"
    else
      printf 'events\t0\n'
    fi
  } > "$file"
}

# compose_pressure_ramp <file> <seconds-up> <peak-pressure-pa> <pre-stress>: mouth pressure ramps
# 0 -> peak over <seconds-up>, then back to 0 over the same time, at a fixed lip pre-stress. The
# pressure at which the note starts on the way up and stops on the way down is the oscillation
# threshold and the extinction pressure; the gap between them is the hysteresis.
compose_pressure_ramp() {
  local file=$1 up=$2 peak=$3 tension=$4
  {
    printf 'frankensim-lip-performance-v1\naudio %d %d\n' "$RATE" $((RATE * up * 2))
    sed -E "s/^(bessel .* )16\$/\1$PIECES/" "$DIR/tenor-trombone.head"
    printf 'schedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t2\n'
    printf 'track\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t2\n'
    printf 'event\t0e0\t%se0\tpressure-pa\t%s\n' "$up" "$peak"
    printf 'event\t%se0\t%se0\tpressure-pa\t0e0\n' "$up" "$up"
    printf 'track\tlip\tpre-stress\ninitial\tmultiplier\t%s\nevents\t0\n' "$tension"
  } > "$file"
}

# At most 8 renders at once (lane budget on the shared trj host).
MAX_JOBS=8
render() {
  local name=$1
  while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do wait -n || true; done
  { "$BIN" lip "$OUT/$name.performance" "$OUT/$name.wav" > "$OUT/$name.log" 2>&1 \
      && echo "rendered $name" || echo "FAILED $name (see $OUT/$name.log)"; } &
}

compose "$OUT/sweep-p1500.performance" 8 "$PRESSURE" 5e-1 2e0
render sweep-p1500
# Threshold and extinction pressure at the mid pre-stress (up to 4 kPa, 8 s each way).
compose_pressure_ramp "$OUT/pramp-t1e0.performance" 8 4e3 1e0
render pramp-t1e0
for tension in $HOLD_TENSIONS; do
  name="hold-t${tension}-p1500"
  compose "$OUT/$name.performance" 3 "$PRESSURE" "$tension"
  render "$name"
done
wait
ls -l "$OUT"/*.wav
