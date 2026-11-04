#!/usr/bin/env bash
# Native lip-valve renders of the tenor trombone description (pexfp.254, lane P5).
# Run on ts1. The binary must contain `music_render lip` (the brass_loop::performance owner):
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
#  * pramp-tNNN: mouth pressure ramps 0 -> 12 kPa -> 0 at a fixed lip pre-stress. A diagnostic: at this rate
#    the pressures where the note starts and stops are mostly lag around one threshold, the time a small
#    signal takes to grow and a note takes to decay, so a ramp can't place a sustain (see the rule below).
#  * hold-tNNN: the same lip held at a fixed pre-stress with an attack (see ATTACK_PRESSURE below). The
#    player's embouchure setting is the free player input; which setting yields which note is read from
#    the render afterwards, never set from a pitch.
#
# Mouth pressure. The published mouth pressures are 1.4 kPa +-10 % for one trombonist playing mp on B-flat2 (Boutin
# et al. 2015, JASA 137:1195, Sec. II B 2: the subjects were asked to play the notes mp) and 14-16 mbar for an artificial mouth on E3 (Bromage 2007 PhD Edinburgh
# Table 4.3); per-note and per-dynamic tables for trombone are UNKNOWN. The first native renders of THIS
# lip (Cullen's artificial lips, Table A.2 embouchure 7) at 1.5 kPa, and a ramp to 4 kPa, never started:
# the flow stayed steady and the differentiating receiver radiated exactly zero. This lip needs several
# kPa, so the ramps measure the threshold of the model directly instead of assuming a mouth pressure.
#
# Needs a music_render with the DC-bin fix of the mouth-flow line (bin da31fbd67ee9 or later): without it the
# default mouth-flow ladder cannot converge on this bore (see the head).
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
DIR=data/instruments/descriptions/trombone
RATE=19200
# Each hold is "pre-stress:sustain-Pa" (see the sustain rule below).
HOLDS=${HOLDS-"5e-1:6e3 6e-1:6e3 7e-1:6e3 8.5e-1:6e3 1e0:7.132e3 1.15e0:6.22e3"}
# The first native renders (1.5 kPa holds at pre-stress 0.6 to 1.8, and a ramp to 4 kPa at 1.0) never
# started the lip: the mouth flow stayed steady and the receiver, which differentiates flow, radiated
# exactly zero. These ramps find where THIS model's lip first speaks.
RAMP_TENSIONS=${RAMP_TENSIONS-"7e-1 1e0 1.3e0 1.6e0 2e0 2.5e0"}
RAMP_PEAK=${RAMP_PEAK:-1.2e4}
# Held notes use an attack: mouth pressure rises to ATTACK_PRESSURE in 50 ms, stays there for ATTACK_SECONDS,
# then relaxes over 0.5 s to the hold's sustain.
#
# The sustain rule. Every hold stood at 6 kPa until the reflection moved to the junction's own sample, which
# raised the small-signal threshold of two of them: pre-stress 1.0 from 5664 to 6733 Pa, and 1.15 from 8161 to
# 8460 Pa. A hold left under its threshold keeps the ratio of sustain to threshold it had: 6000 / 5664 x 6733 =
# 7132 Pa, and 6000 / 8161 x 8460 = 6220 Pa. Every other hold is still above its threshold and stays at 6 kPa.
# A threshold is the root of this model's linearized voice on its own bore response, not a recording and not a
# published mouth pressure (those belong to other lips; see the note above). Pre-stress 1.15 sat under its
# threshold before the move too, so the rule keeps a ratio below one and that note dies away.
ATTACK_PRESSURE=${ATTACK_PRESSURE:-1.2e4}
ATTACK_SECONDS=${ATTACK_SECONDS:-8e-1}
HOLD_SECONDS=${HOLD_SECONDS:-4}
SWEEP_PRESSURE=${SWEEP_PRESSURE:-8e3}

cd "$SNAP"
if [ "$PIECES" != 16 ]; then
  OUT="$OUT/pieces-$PIECES"
fi
if [ -n "${CUP_VOLUME_M3:-}" ]; then
  OUT="$OUT/cup-$CUP_VOLUME_M3"
fi
mkdir -p "$OUT"

# head_text: the description head with the Bessel piece count (PIECES) and, when CUP_VOLUME_M3 is set, the
# lumped cup volume replaced. The head's own 12.7 ml is Cullen's rig mouthpiece; Braden's mouthpiece-and-leadpipe
# element starts at a 1.286 cm radius (a cup rim), so the cup may be counted twice (once lumped, once inside
# that element, whose interior is UNKNOWN). CUP_VOLUME_M3=1e-9 removes the lumped cup: a sensitivity bracket
# over an unknown, not a choice.
head_text() {
  sed -E "s/^(bessel .* )16\$/\1$PIECES/" "$DIR/tenor-trombone.head" |
    if [ -n "${CUP_VOLUME_M3:-}" ]; then sed -E "s/^cup_volume_m3 .*/cup_volume_m3 $CUP_VOLUME_M3/"; else cat; fi
}

# compose <file> <seconds> <mouth-pressure-pa> <initial-pre-stress> [<final-pre-stress>]
compose() {
  local file=$1 seconds=$2 pressure=$3 initial=$4 final=${5:-}
  {
    printf 'frankensim-lip-performance-v1\naudio %d %d\n' "$RATE" $((RATE * seconds))
    head_text
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
# 0 -> peak over <seconds-up>, then back to 0 over the same time, at a fixed lip pre-stress. See
# pramp-tNNN at the top for what its start and stop pressures do and don't measure.
compose_pressure_ramp() {
  local file=$1 up=$2 peak=$3 tension=$4
  {
    printf 'frankensim-lip-performance-v1\naudio %d %d\n' "$RATE" $((RATE * up * 2))
    head_text
    printf 'schedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t2\n'
    printf 'track\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t2\n'
    printf 'event\t0e0\t%se0\tpressure-pa\t%s\n' "$up" "$peak"
    printf 'event\t%se0\t%se0\tpressure-pa\t0e0\n' "$up" "$up"
    printf 'track\tlip\tpre-stress\ninitial\tmultiplier\t%s\nevents\t0\n' "$tension"
  } > "$file"
}

# compose_attack <file> <seconds> <pre-stress> <sustain-pa>: the held-note gesture described above.
compose_attack() {
  local file=$1 seconds=$2 tension=$3 sustain=$4
  {
    printf 'frankensim-lip-performance-v1\naudio %d %d\n' "$RATE" $((RATE * seconds))
    head_text
    printf 'schedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t2\n'
    printf 'track\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t2\n'
    printf 'event\t0e0\t5e-2\tpressure-pa\t%s\n' "$ATTACK_PRESSURE"
    printf 'event\t%s\t5e-1\tpressure-pa\t%s\n' "$ATTACK_SECONDS" "$sustain"
    printf 'track\tlip\tpre-stress\ninitial\tmultiplier\t%s\nevents\t0\n' "$tension"
  } > "$file"
}

# At most 8 renders at once (lane budget on the shared ts1 host).
MAX_JOBS=8
render() {
  local name=$1
  while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do wait -n || true; done
  { "$BIN" lip "$OUT/$name.performance" "$OUT/$name.wav" > "$OUT/$name.log" 2>&1 \
      && echo "rendered $name" || echo "FAILED $name (see $OUT/$name.log)"; } &
}

compose "$OUT/sweep-p$SWEEP_PRESSURE.performance" 8 "$SWEEP_PRESSURE" 5e-1 2e0
render "sweep-p$SWEEP_PRESSURE"
# Threshold and extinction pressure at each pre-stress (8 s up, 8 s down).
for tension in $RAMP_TENSIONS; do
  compose_pressure_ramp "$OUT/pramp-t$tension.performance" 8 "$RAMP_PEAK" "$tension"
  render "pramp-t$tension"
done
for hold in $HOLDS; do
  IFS=: read -r tension sustain <<< "$hold"
  name="hold-t${tension}"
  compose_attack "$OUT/$name.performance" "$HOLD_SECONDS" "$tension" "$sustain"
  render "$name"
done
wait
ls -l "$OUT"/*.wav
