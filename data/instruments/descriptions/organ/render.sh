#!/usr/bin/env bash
# Native renders of the committed organ-pipe performances (pexfp.381, lane P4). Checked against the pipe's own measurements, not graded: the
# U4 organ references are registrations at Eb3, A2 and Gb6, none within an octave of this pipe. Run on a Linux host. The performance files
# must exist in the snapshot of a pushed SHA; the music_render binary may come from an earlier green SHA (BIN_SHA12):
#   SHA12=<sha12 containing these files> BIN_SHA12=<sha12 with a release music_render> data/instruments/descriptions/organ/render.sh
# Each file is 48 kHz mono, observed as the summed pressure of the open end and the mouth at a physical receiver 1 m above the mouth.
#
# Geometry: data/instruments/geometry/organpipe-tin-principal-castellengo2004 (an open tin flue pipe measured by M. Castellengo, ISMA 2004,
# Nara, Sec. 2; every number is a dimensions.tsv row). The network is the tin whistle's without finger holes: the mouth channel through the
# wall (a riser of equal area to the 21 x 7 mm mouth), a junction, and the 27 mm bore 308.5 mm from the mouth centre to the open end. The
# mouth's mass is computed, not a constant: `window wall-cut 12 cheeks parallel` solves the 3-D potential flow through the 21 mm wide by 7 mm
# cut in the wall of the 13.5 mm bore, side walls parallel, closed flush at the mouth's lower edge (the languid) and free space outside
# (fs_bem::wall_window, a staircase of cut cells at 12, 24 and 48 axial cells, extrapolated at the order the three rungs show), in place of
# the 0.61 a_eq mouth inertance; `mouth_reaction window-cylinder` radiates the cut's outer face from the rigid outer cylinder.
# The wall is not stated (u-wall-thickness): the cards carry 0.5 mm, the low end of the bracket that the same laboratory's two metal
# Principals give (0.5 and 0.7 mm; Acustica 85, 1999, Table I). WALLS="0.5 0.7" also renders the high end into $OUT/w07: its clock is 576 kHz
# and its channel 0.7 mm, nothing else changes. The window is 31.78 mm of bore at 0.5 mm (rungs 32.57, 32.10, 31.91 mm, order 1.29) and
# 32.91 mm at 0.7 mm (order 1.50).
# Wind: the measured foot pressure of the same pipe on its organ chest (Sec. 3): 360 Pa steady. The short attack passes through the two
# stated points of the paper's short transient, 125 Pa at 9 ms and 360 Pa at 45 ms, piecewise linearly (u-attack-profile); the rise cards
# take the paper's four-second slow rise from 0 to 360 Pa (its end pressure is not stated; 360 Pa is the pipe's steady pressure). Dry
# room air at 20 C: organ wind is room air, and the room's temperature and humidity are not stated (u-environment). The foot's own acoustics
# are not on the route (u-foot). No windway (`windway none`, a uniform exit at the Bernoulli speed): the flue's depth is not stated
# (u-flue-depth). Labium offset and jet angle zero (UNKNOWN). Flue height 0.25 mm, the centre of the paper's 0.2 to 0.3.
# Jet, two arms, both rendered by default (ARMS="authored forced"):
#   authored: the repository interim slot-jet card with the authored displacement gain 90, as on every other flue card;
#   forced:   the same card with `jet_response forced-bickley 0.1 180`, the forced Bickley jet (rulings 7376, 7518; uniform exit:
#             b = 0.75 h, U0 = U). 0.1 is the grid spacing [b] and 180 the recorded window [b/U0], numerical controls, not physics. The
#             production row is 0.1 90; the window is 180 here because this jet crosses 37 b to the labium against about 9 on the
#             whistle. NUMERICS=1 also renders the forced short card at spacing 0.05 and 0.2 (window 180) and window 90 and 360
#             (spacing 0.1) into $OUT/numerics. Spacing 0.2 moves the 4.68 kHz line 0.8 cents. The window does NOT converge on this
#             pipe: at 90, 180 and 360 the dominant line is 6.37, 4.68 and 7.44 kHz. So the forced arm's sounding mode is not a
#             result, only that it is a high one (status.json).
#
# Numerical controls forced by the generic route (none is a physical parameter, none was tuned): mechanics clock 720 kHz, because the mouth
# channel (the 0.5 mm wall) must span one transit sample (at least 686 kHz); jet history 262144 samples (65536 is refused at 720 kHz for a
# 7 mm cut-up at the 0.2 m/s minimum jet speed); window cells 12 (8 is refused: its rungs round the arc differently; 10 gives 0.04 mm less
# bore); compact baffled radiation at the open end, which holds for ka <= 0.5 and so caps the receiver band at 2.0 kHz on the 13.5 mm end;
# the pressure release lasts 10.5 ms (a release that ends exactly on a control tick can leave a rounding-sized pressure that the jet domain
# check refuses).
set -euo pipefail

SHA12=${SHA12:?set SHA12 to the pushed sha12 whose snapshot holds these files}
BIN_SHA12=${BIN_SHA12:-$SHA12}
BIN=${BIN:-/data/tmp/wa-inst/bin/$BIN_SHA12/music_render}
SNAP=${SNAP:-/data/tmp/wa-inst/snap/$SHA12/frankensim}
OUT=${OUT:-/data/tmp/wa-inst/renders/organ/$SHA12}
CARDS=${CARDS:-"short rise"}
ARMS=${ARMS:-"authored forced"}
WALLS=${WALLS:-"0.5"}
NUMERICS=${NUMERICS:-0}
DIR=data/instruments/descriptions/organ

cd "$SNAP"
render() { # card file, output directory
    local name
    name=$(basename "$1" .performance)
    mkdir -p "$2"
    "$BIN" wind "$1" "$2/$name.wav" | tee "$2/$name.json"
}
for arm in $ARMS; do
    for card in $CARDS; do
        file="$DIR/organ-principal-$card.performance"
        [ "$arm" = forced ] && file="$DIR/organ-principal-forced-$card.performance"
        for wall in $WALLS; do
            case "$wall" in
            0.5) render "$file" "$OUT/w05" ;;
            0.7)
                mkdir -p "$OUT/w07"
                variant="$OUT/w07/$(basename "$file")"
                sed -e 's/^clock 720000 48000$/clock 576000 48000/' \
                    -e 's/^duct_section 0 1 0\.0005000 /duct_section 0 1 0.0007000 /' "$file" > "$variant"
                render "$variant" "$OUT/w07"
                ;;
            *) echo "WALLS takes 0.5 and 0.7" >&2; exit 2 ;;
            esac
        done
    done
done
if [ "$NUMERICS" = 1 ]; then
    mkdir -p "$OUT/numerics"
    for row in "0.05 180" "0.2 180" "0.1 90" "0.1 360"; do
        variant="$OUT/numerics/organ-principal-forced-short-${row// /-}.performance"
        sed -e "s/^jet_response forced-bickley 0\.1 180$/jet_response forced-bickley $row/" \
            "$DIR/organ-principal-forced-short.performance" > "$variant"
        render "$variant" "$OUT/numerics"
    done
fi
