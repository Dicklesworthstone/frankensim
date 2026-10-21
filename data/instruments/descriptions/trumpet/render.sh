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
#  * sweep-pNNN: lip pre-stress ramps 0.7 -> 3.5 at a held mouth pressure (default 14 kPa). A diagnostic, not a note: it records
#    which regime the bore + lip speaks as the lip tightens.
#  * pramp-tNNN: mouth pressure ramps 0 -> 25 kPa -> 0 at a fixed lip pre-stress. The pressure at which
#    the note starts on the way up and stops on the way down is the oscillation threshold and the
#    extinction pressure; their gap is the hysteresis.
#  * hold-tNNN-pNNN: the same lip held at a fixed pre-stress and mouth pressure for 3 s. The embouchure
#    setting is the free player input; which setting yields which note is READ FROM THE RENDER afterwards.
#
# Mouth pressure is a CITED player input. Fletcher and Tarnopolsky (1999, JASA 105:874) measured the threshold
# blowing pressure of trumpet players as p0 = K f (their Eq. 1 and Fig 4, p.877: K = 7 Pa/Hz for the professional
# GC, 8 and 4 for two amateurs) and found loud playing at 3.3 to 3.9 times threshold (their Fig 6: C4 2.0 / 3.9 /
# 6.6 kPa, C5 3.3 / 6.3 / 13.0, C6 6.7 / 11.0 / 24.0 kPa). The held notes use 3.5 x 7 Pa/Hz x f: Bb3 (233.08 Hz)
# 5.71 kPa, F4 (349.23 Hz) 8.56, Bb4 (466.16 Hz) 11.42, D5 (587.33 Hz) 14.39 kPa. The 25 kPa ramp peak is the
# highest blowing pressure they measured (Fig 1, p.875); the ramps measure the threshold of THIS lip directly.
# Earlier revisions held all four notes at 14 kPa, which gave radiated levels near the measured ones but a
# shank fundamental 11 dB too strong; at the cited pressures this model's notes are 10 to 21 dB(A) quieter
# and duller than F&T measured (the status card has both results), and that is the setting kept here. The
# lip is Cullen (2000) embouchure 19 (see the head for why, and for the disclosure that playability entered
# the choice): the first renders used embouchure 7, whose upper open notes started only at 19 to 27.5 kPa and
# sat on the extinction edge at 21 kPa (the D5 hold fell silent).
#
# Needs a music_render with the DC-bin fix of the mouth-flow line (bin da31fbd67ee9 or later): without it the
# default mouth-flow ladder cannot converge on this bore (see the head).
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
RAMP_TENSIONS=${RAMP_TENSIONS-"1e0 1.25e0 1.5e0 1.75e0 2e0 2.5e0"}
RAMP_PEAK=${RAMP_PEAK:-2.5e4}
# Held notes use an attack: mouth pressure rises to an ATTACK pressure in 50 ms (above the onset the pressure
# ramps measured at that pre-stress), stays there for ATTACK_SECONDS, then relaxes over 0.5 s to a SUSTAIN
# pressure (inside the hysteresis band: above the extinction pressure, below the onset). Each hold is
# "pre-stress:attack-Pa:sustain-Pa". The SUSTAIN is the cited Fletcher and Tarnopolsky pressure above. This lip
# needs more than its sustain pressure to start: a slow 0 -> 16 -> 0 kPa ramp at the landed pre-stress starts
# Bb3, F4, Bb4 and D5 at 5.7, 11.7, 9.7 and 13.0 kPa and stops them at 2.3, 3.0, 3.7 and 7.0 kPa (a diagnostic;
# their players start a note with an attack and settle at the pressure they hold, which the hysteresis
# reproduces), so each note starts with a 14 kPa attack for ATTACK_SECONDS (a player's tongued attack, a gesture
# and not a measured value) and relaxes over 0.5 s to its sustain; D5's sustain is already the attack. The
# pre-stress of each hold was set from the emergent pitch of a first
# render against equal temperament (a player tuning an embouchure; pitch is never an input): first guesses
# 0.9, 1.38, 1.86, 2.38 gave Bb3 33 cents flat, F4 26, Bb4 21, D5 21 cents flat; 0.95, 1.42, 1.9, 2.42 gave
# Bb3 +20, F4 +5, Bb4 +2, D5 -1 cents; Bb3 at 0.92 gave -16 and at 0.93 -5 cents (steady-state, 2 to 4 s).
# The emergent pitches are the bore's: its passive impedance peaks (music_render lip --impedance) lie 41 to
# 62 cents flat of the equal-tempered open series for partials 3 to 6, which the lip pulls toward the note.
ATTACK_SECONDS=${ATTACK_SECONDS:-8e-1}
HOLD_SECONDS=${HOLD_SECONDS:-4}
HOLDS=${HOLDS-"9.3e-1:1.4e4:5.71e3 1.42e0:1.4e4:8.56e3 1.9e0:1.4e4:1.142e4 2.42e0:1.439e4:1.439e4"}
SWEEP_PRESSURE=${SWEEP_PRESSURE:-1.4e4}

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

# compose_attack <file> <seconds> <pre-stress> <attack-pa> <sustain-pa>: the held-note gesture described above.
compose_attack() {
  local file=$1 seconds=$2 tension=$3 attack=$4 sustain=$5
  {
    printf 'frankensim-lip-performance-v1\naudio %d %d\n' "$RATE" $((RATE * seconds))
    cat "$DIR/trumpet-open.head"
    printf 'schedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t2\n'
    printf 'track\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t2\n'
    printf 'event\t0e0\t5e-2\tpressure-pa\t%s\n' "$attack"
    printf 'event\t%s\t5e-1\tpressure-pa\t%s\n' "$ATTACK_SECONDS" "$sustain"
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
for hold in $HOLDS; do
  IFS=: read -r tension attack sustain <<< "$hold"
  name="hold-t${tension}"
  compose_attack "$OUT/$name.performance" "$HOLD_SECONDS" "$tension" "$attack" "$sustain"
  render "$name"
done
wait
ls -l "$OUT"/*.wav
