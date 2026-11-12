#!/usr/bin/env bash
# Native lip-valve renders of the Bb horn description, valves open (lane: horn, MB 27394 and 27456).
# Run on ts1 (or any host with the binary below). The binary must contain `music_render lip` and every row of the head.
#   BIN=<music_render> SNAP=<checkout holding these files> PACKS=<gas pack directory> OUT=<new directory> data/instruments/descriptions/horn/render.sh
#   DRY=1 composes the cards only (no render): the way to check a card against the one that was rendered.
#
# Composition: horn-bb-open.head (mouthpiece, bore, lip, microphone, each number labelled) plus, per note, the lip-area exponent row and a canonical flat-hold gesture schedule written below: the mouth pressure and the lip
# pre-stress. There is no pitch, note or frequency in any input: the note is the lock the loop finds between the lip and the bore. The open Bb series only (the valve loops are not measured).
#
# Rows byte for byte the live horn cards' (the D4 card of 11-09): the Cullen (2000) embouchure 19 lip, the bore gas pack, the flanged modal radiation load, the constriction noise, the observer gas, the declared
# mouth-flow tolerance 4e-6, 19.2 kHz clock, 4 s. The blow of the stages discovery and final is the 6.22 kPa flat blow rising in 0.05 s: a hold BORROWED from the old trombone recipe (trombone/render.sh, "The sustain rule": 6 kPa chosen for an
# artificial lip that did not start at 1.5 kPa, rescaled to 6.22 kPa when that model's threshold moved; not measured by anyone, and not a horn's); those two stages are the ladder that picked the lip settings and are kept as they were run.
#
# THE POSTED NOTES (stage gesture) blow a horn player's MEASURED gesture instead: Stevenson (2009, Edinburgh PhD, Figures 6.7 and 6.8): F3 on a horn, one player (JC), no tonguing, a transparent mouthpiece, dynamic not stated; read from the PDF's vector data
# (Figure 6.7, no valves): the mouth pressure is 1.56 kPa at t = 0, rises at 7.84 kPa/s to the 2.0 kPa at which the oscillation begins (50 ms) and averages 2.32 kPa over 150 to 200 ms. The thesis text prints the rise as 10 Pa/s; the unit read from the figure's axes is
# kPa/s (that reading is this lane's). The card blows 1.56 kPa at t = 0, a rise of 0.76 kPa in 9.7e-2 s (7.84 kPa/s), then 2.32 kPa held. It is MEASURED at F3; the same gesture on Bb3, D4 and F4 is a labelled EXTRAPOLATION of one player's F3 measurement. The lip rows and
# pre-stress of each note are the live ones (picked by the ladder at 6.22 kPa); a card of this stage differs from the live card of its note in the two lines of the blow track only, and the held pitch moved by 1.7 cents or less.
#
# TWO STAGES, because the pre-stress of the lip is the free player input and which partial a setting plays is read from the render:
#  * discovery (d-mNNN): 22 renders of 1 s at the pre-stress settings listed in DISCOVERY below (0.30 to 0.46 and 0.50 to 0.62 in steps of 0.04, 0.70 to 0.86, 0.90, 0.94, 1.00, 1.04, 1.10, 1.14, 1.18, 1.22), each with the per-note lip exponent of the partial it was expected to reach. Read: the fundamental of 0.5 to 1.0 s
#    by a harmonic comb search around the engine's passive impedance maxima (music_render lip CARD OUT.wav --impedance OUT.tsv: 42.63, 114.62, 169.97, 229.55, 286.90, 343.47, 400.58, 456.19 Hz for this head).
#    Result at 0.5 to 1.0 s: Bb2 at 0.30 to 0.42 (115.19 to 116.39 Hz); F3 at 0.50 to 0.62 (171.81 to 172.36 Hz); Bb3 at 0.74 to 0.90 (230.90 to 234.96 Hz); D4 at 1.00 to 1.14 (289.29 to 293.81 Hz; a test render at 0.97 gave 288.9 Hz); F4 at 1.18 and 1.22
#    (345.66 and 346.08 Hz); 0.46, 0.70 and 0.94 are not steady (the first half is another partial).
#  * final (f-NOTE-mNNN): four renders per note of 4 s at settings that bracket the discovery result. The setting of a posted note is the steady one whose pitch is nearest the equal-tempered note: the embouchure a player
#    would pick; the pitch is never an input of the card.
#  * gesture (ges-NOTE; STAGE=gesture, not part of STAGE=both): the four posted notes, 4 s each, at the lip settings of the posted notes (0.68, 0.86, 1.14, 1.34) and the measured gesture above.
#
# The lip-area exponent q is the trombonist's, per note, from Bromage, Campbell and Gilbert (2010), Table 2 (1.2 at Bb1, 1.3 at Bb2, 1.6 at F3, 1.9 at F4), linear in MIDI between the printed notes, none above F4:
# F3 1.6, Bb3 1.725 (53 to 65 at 58), D4 1.825, F4 1.9; Bb2 1.3. Borrowed, not a horn player's.
#
# References (grading only, never an input): UIowa MIS horn files, ids in references-screened.json (refscreen-4). Their screen status there: horn-uiowa-horn-ff-c3b3-F3-05 and -Bb3-10 (the only F3 and Bb3 clips, ff) are INVALID
# (failed snr and pitch); D4 and F4 have ff and mf clips INVALID (snr and pitch) and pp clips VALID (horn-uiowa-horn-pp-c4b4-D4-02 and -F4-05). A clip that is INVALID claims nothing in the manifest.
set -euo pipefail

BIN=${BIN:?set BIN to a music_render with the lip owner}
SNAP=${SNAP:?set SNAP to the checkout whose data/instruments/descriptions/horn holds the head}
PACKS=${PACKS:?set PACKS to the directory of the N2, O2, Ar and CO2 species, nasa9 and transport packs}
OUT=${OUT:?set OUT to a new directory}
DRY=${DRY:-0}
STAGE=${STAGE:-both}
DIR=data/instruments/descriptions/horn
RATE=${RATE:-19200}
MAX_JOBS=${MAX_JOBS:-24}
DISCOVERY="d-m030:1.3e0:3e-1 d-m034:1.3e0:3.4e-1 d-m038:1.3e0:3.8e-1 d-m042:1.3e0:4.2e-1 d-m046:1.3e0:4.6e-1 d-m050:1.6e0:5e-1 d-m054:1.6e0:5.4e-1 d-m058:1.6e0:5.8e-1 d-m062:1.6e0:6.2e-1 d-m070:1.725e0:7e-1 d-m074:1.725e0:7.4e-1 d-m078:1.725e0:7.8e-1 d-m082:1.725e0:8.2e-1 d-m086:1.725e0:8.6e-1 d-m090:1.825e0:9e-1 d-m094:1.825e0:9.4e-1 d-m100:1.825e0:1e0 d-m104:1.825e0:1.04e0 d-m110:1.9e0:1.1e0 d-m114:1.9e0:1.14e0 d-m118:1.9e0:1.18e0 d-m122:1.9e0:1.22e0"
FINAL="f-F3-m058:1.6e0:5.8e-1 f-F3-m062:1.6e0:6.2e-1 f-F3-m065:1.6e0:6.5e-1 f-F3-m068:1.6e0:6.8e-1 f-Bb3-m082:1.725e0:8.2e-1 f-Bb3-m086:1.725e0:8.6e-1 f-Bb3-m088:1.725e0:8.8e-1 f-Bb3-m090:1.725e0:9e-1 f-D4-m104:1.825e0:1.04e0 f-D4-m110:1.825e0:1.1e0 f-D4-m112:1.825e0:1.12e0 f-D4-m114:1.825e0:1.14e0 f-F4-m122:1.9e0:1.22e0 f-F4-m126:1.9e0:1.26e0 f-F4-m130:1.9e0:1.3e0 f-F4-m134:1.9e0:1.34e0"

GESTURE="ges-F3:1.6e0:6.8e-1 ges-Bb3:1.725e0:8.6e-1 ges-D4:1.825e0:1.14e0 ges-F4:1.9e0:1.34e0"
# Stevenson (2009) Figure 6.7, F3 on a horn: initial mouth pressure, the time of the rise and the held pressure, in Pa and s (see the header).
GESTURE_BLOW="1.56e3 9.7e-2 2.32e3"

cd "$SNAP"
mkdir -p "$OUT"

# compose <file> <seconds> <lip-area-exponent> <pre-stress> [<initial-Pa> <rise-s> <hold-Pa>]: gesture numbers are in canonical shortest e-notation, as the schedule reader requires.
# Without the blow arguments the blow is the flat 6.22 kPa rising in 0.05 s of the stages discovery and final.
compose() {
  local file=$1 seconds=$2 q=$3 tension=$4 p0=${5:-0e0} rise=${6:-5e-2} hold=${7:-6.22e3}
  {
    printf 'frankensim-lip-performance-v1\naudio %d %d\n' "$RATE" $((RATE * seconds))
    cat "$DIR/horn-bb-open.head"
    printf 'lip_area_exponent %s\n' "$q"
    printf 'schedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t2\n'
    printf 'track\tblow\tblowing-pressure\ninitial\tpressure-pa\t%s\nevents\t1\nevent\t0e0\t%s\tpressure-pa\t%s\n' "$p0" "$rise" "$hold"
    printf 'track\tlip\tpre-stress\ninitial\tmultiplier\t%s\nevents\t0\n' "$tension"
  } > "$file"
}

GAS=()
for sp in N2 O2 Ar CO2; do GAS+=(--gas-component "$PACKS/$sp.species.pack" "$PACKS/$sp.nasa9.pack" "$PACKS/$sp.transport.pack"); done
GAS+=(--gas-water none --gas-phase frozen)

render() {
  local name=$1
  [ "$DRY" = 1 ] && return 0
  while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do wait -n || true; done
  { nice -n 10 "$BIN" lip "$OUT/$name.performance" "$OUT/$name.wav" "${GAS[@]}" > "$OUT/$name.log" 2>&1 \
      && echo "rendered $name" || echo "FAILED $name (see $OUT/$name.log)"; } &
}

if [ "$STAGE" = both ] || [ "$STAGE" = discovery ]; then
  for entry in $DISCOVERY; do
    IFS=: read -r name q tension <<< "$entry"
    compose "$OUT/$name.performance" 1 "$q" "$tension"
    render "$name"
  done
fi
if [ "$STAGE" = both ] || [ "$STAGE" = final ]; then
  for entry in $FINAL; do
    IFS=: read -r name q tension <<< "$entry"
    compose "$OUT/$name.performance" 4 "$q" "$tension"
    render "$name"
  done
fi
if [ "$STAGE" = gesture ]; then
  read -r G_P0 G_RISE G_HOLD <<< "$GESTURE_BLOW"
  for entry in $GESTURE; do
    IFS=: read -r name q tension <<< "$entry"
    compose "$OUT/$name.performance" 4 "$q" "$tension" "$G_P0" "$G_RISE" "$G_HOLD"
    render "$name"
  done
fi
wait
ls -l "$OUT"/*.performance | wc -l
