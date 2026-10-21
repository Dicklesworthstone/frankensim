#!/usr/bin/env bash
# Native lip-valve renders of the Bb trumpet description, valves up (pexfp.254, lane P5).
# Run on trj. The binary must contain `music_render lip` (the brass_loop::performance owner):
#   SHA12=<sha12 whose snapshot holds these files> BIN_SHA12=<sha12 of a music_render with `lip`> \
#     data/instruments/descriptions/trumpet/render.sh
#
# Composition: trumpet-open.head (bore, lip, microphone, each number with its source) plus a canonical
# gesture schedule written below: the player's mouth pressure and lip pre-stress. VALVES UP ONLY: the
# bore file is straightened, so only the open harmonics are valid. There is no pitch, note or frequency
# anywhere in an input: the note is the lock the loop finds between the lip and the bore.
#
#  * sweep-pNNN: lip pre-stress ramps 0.7 -> 3.5 at a held mouth pressure (default 8 kPa). A diagnostic, not a note: it records
#    which regime the bore + lip speaks as the lip tightens.
#  * pramp-tNNN: mouth pressure ramps 0 -> 12 kPa -> 0 at a fixed lip pre-stress. The pressure at which
#    the note starts on the way up and stops on the way down is the oscillation threshold and the
#    extinction pressure; their gap is the hysteresis.
#  * hold-tNNN-pNNN: the same lip held at a fixed pre-stress and mouth pressure for 3 s. The embouchure
#    setting is the free player input; which setting yields which note is READ FROM THE RENDER afterwards.
#
# Mouth pressure: threshold p0 = K f with K of 4 to 8 Pa/Hz for three trumpeters (Fletcher and
# Tarnopolsky 1999, JASA 105:874, p.877), and 4 to 13 kPa pp to ff on high notes (Bouhuys, as quoted
# there); per-note mf tables are UNKNOWN. The first trombone renders with this same lip never started at
# 1.5 kPa (nor to 4 kPa), so the ramps go to 12 kPa and measure the threshold of THIS model directly;
# holds are placed above the measured threshold.
#
# Clock: 32 kHz mechanics = output, not 48 kHz. The fs-duct flanged radiation load is admitted only to
# ka = 20, and the 0.061 m bell reaches ka = 20 at 17.9 kHz, so Nyquist 16 kHz (ka 17.9) is admitted and
# 48 kHz (Nyquist 24 kHz, ka 26.8) would be refused. The radiation receiver's resolved band is 6.9 kHz
# (see the microphone record). The WAV is written at 32000 Hz and never resampled; grade up to 6.9 kHz.
#
# References (grading only, never an input): UIowa MIS Trumpet.novib.mf.* range files (the recording
# uses a vibrato set too), ids uiowa-bbtrumpet-trumpet-novib-mf-e3b3 / -c4b4 / -c5b5 in references.json;
# the notes of the open series sound at concert Bb3, F4, Bb4, D5.
set -euo pipefail

SHA12=${SHA12:?set SHA12 to the pushed sha12 whose snapshot holds these files}
BIN_SHA12=${BIN_SHA12:-$SHA12}
BIN=${BIN:-/data/tmp/wa-inst/bin/$BIN_SHA12/music_render}
SNAP=${SNAP:-/data/tmp/wa-inst/snap/$SHA12/frankensim}
OUT=${OUT:-/data/tmp/wa-inst/renders/trumpet/$SHA12}
DIR=data/instruments/descriptions/trumpet
RATE=32000
HOLD_TENSIONS=${HOLD_TENSIONS-"1e0 1.5e0 2.2e0 2.9e0"}
RAMP_TENSIONS=${RAMP_TENSIONS-"1e0 1.5e0 2.2e0 2.9e0"}
RAMP_PEAK=${RAMP_PEAK:-1.2e4}
# Held notes use an attack: mouth pressure rises to ATTACK_PRESSURE in 50 ms (above the onset the pressure
# ramps measured), stays there for ATTACK_SECONDS, then relaxes over 0.5 s to SUSTAIN_PRESSURE (inside the
# hysteresis band: above the extinction pressure, below the onset). Both pressures come from the model's own
# measured threshold map, not from a recording and not from a published mf pressure (the published mouth
# pressures belong to other lips; see the note above).
ATTACK_PRESSURE=${ATTACK_PRESSURE:-1.2e4}
ATTACK_SECONDS=${ATTACK_SECONDS:-8e-1}
SUSTAIN_PRESSURE=${SUSTAIN_PRESSURE:-6e3}
HOLD_SECONDS=${HOLD_SECONDS:-4}
SWEEP_PRESSURE=${SWEEP_PRESSURE:-8e3}

cd "$SNAP"
mkdir -p "$OUT"

# compose <file> <seconds> <mouth-pressure-pa> <initial-pre-stress> [<final-pre-stress>]
compose() {
  local file=$1 seconds=$2 pressure=$3 initial=$4 final=${5:-}
  {
    printf 'frankensim-lip-performance-v1\naudio %d %d\n' "$RATE" $((RATE * seconds))
    cat "$DIR/trumpet-open.head"
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

# compose_pressure_ramp <file> <seconds-up> <peak-pressure-pa> <pre-stress>
compose_pressure_ramp() {
  local file=$1 up=$2 peak=$3 tension=$4
  {
    printf 'frankensim-lip-performance-v1\naudio %d %d\n' "$RATE" $((RATE * up * 2))
    cat "$DIR/trumpet-open.head"
    printf 'schedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t2\n'
    printf 'track\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t2\n'
    printf 'event\t0e0\t%se0\tpressure-pa\t%s\n' "$up" "$peak"
    printf 'event\t%se0\t%se0\tpressure-pa\t0e0\n' "$up" "$up"
    printf 'track\tlip\tpre-stress\ninitial\tmultiplier\t%s\nevents\t0\n' "$tension"
  } > "$file"
}

# compose_attack <file> <seconds> <pre-stress>: the held-note gesture described above.
compose_attack() {
  local file=$1 seconds=$2 tension=$3
  {
    printf 'frankensim-lip-performance-v1\naudio %d %d\n' "$RATE" $((RATE * seconds))
    cat "$DIR/trumpet-open.head"
    printf 'schedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t2\n'
    printf 'track\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t2\n'
    printf 'event\t0e0\t5e-2\tpressure-pa\t%s\n' "$ATTACK_PRESSURE"
    printf 'event\t%s\t5e-1\tpressure-pa\t%s\n' "$ATTACK_SECONDS" "$SUSTAIN_PRESSURE"
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

compose "$OUT/sweep-p$SWEEP_PRESSURE.performance" 8 "$SWEEP_PRESSURE" 7e-1 3.5e0
render "sweep-p$SWEEP_PRESSURE"
for tension in $RAMP_TENSIONS; do
  compose_pressure_ramp "$OUT/pramp-t$tension.performance" 8 "$RAMP_PEAK" "$tension"
  render "pramp-t$tension"
done
for tension in $HOLD_TENSIONS; do
  name="hold-t${tension}"
  compose_attack "$OUT/$name.performance" "$HOLD_SECONDS" "$tension"
  render "$name"
done
wait
ls -l "$OUT"/*.wav
