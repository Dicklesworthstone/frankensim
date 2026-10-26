#!/usr/bin/env bash
# Native renders of the committed alto-recorder performances (pexfp.315, lane P4). Listening only: there is no licensed reference recording yet
# (bead pexfp.219), so these are not graded. Run on trj. The performance files must exist in the snapshot of a pushed SHA; the music_render
# binary may come from an earlier green SHA (BIN_SHA12), it does not read these files at build time:
#   SHA12=<sha12 containing these files> BIN_SHA12=<sha12 with a release music_render> \
#     data/instruments/descriptions/recorder/render.sh
# Each file is 48 kHz mono, 2.5 s (the four held notes) or 7 s (the four crescendos), observed as the summed pressure of every open aperture
# (open holes, the open end and the window) at a physical receiver 1 m above the window.
#
# Geometry: data/instruments/geometry/recorder-alto-bressan-brueggen-x-bouterse2010, the alto recorder 'd' by Peter Bressan of the Frans Brueggen
# collection (drawing No. X by Fred Morgan, Zen-On 1981) as tabulated by Jan Bouterse in FoMRHI Comms. 1898, 1910 and 1929 (2009 to 2010): head,
# middle joint and foot bore, thumb hole and six finger holes (positions, sizes, outer diameters), window, labium, step, windway length. Every
# geometry number is a dimensions.tsv row. The network is: the window (wall thickness at the window, DERIVED from the mean outer diameter of the
# other four Bressan altos, as a riser of equal area to the 12.0 x 4.6 mm opening) on a junction with the closed bore up to the block line
# behind it, then the bore as equal-volume cylinders of at most 40 mm with a chimney section at every finger hole (chimney = (outer diameter -
# local bore) / 2), to the open end of the foot. Hole 7 (5.25 mm, chimney (32.0 - 14.02) / 2 mm) opens 24 mm below the foot's upper end:
# its position is not published, the dimensions row is DERIVED from the 16.0 mm tenon socket and the 32.0 mm bulb, with the bracket 18 to
# 30 mm in unknowns.tsv (u-hole-7-position); it is open in a1, b-flat1 and b1 and closed in f1.
# The window's mass is computed, not a constant: each card's row `window wall-cut 4 cheeks parallel ramp 0.00051 18.6210` solves the 3-D
# potential flow through the cut in the 6.575 mm wall of the 9.83 mm bore: parallel cheeks 12.0 mm apart, 4.6 mm long, its downstream face
# the labium ramp, straight up to the edge 0.51 mm above the bore (windway exit 1.06 - step 0.8 + edge thickness 0.25 mm, so the edge
# inherits u-windway-height, bounds 0 to 1.25 mm) and then rising 6.065 mm over the side walls' 18.0 mm beyond the window (22.6 - 4.6 mm,
# 18.62 degrees), with the closed bore behind it and free space outside (fs_bem::wall_window, a staircase of cut cells at 4, 8 and 16
# axial cells, extrapolated at the order the three rungs show: 1.61; rungs 8, 16 and 32 give 52.567 mm at order 1.62, too costly for
# a render). That mass, 52.56 mm of bore, replaces the 0.61 a_eq mouth inertance; the jet's source coefficient keeps its own law. The network
# keeps the 6.575 mm channel only as a one-sample inlet line (1.20 mm at 288 kHz), the computed mass carrying the rest: the channel's full
# length as a 1-D line rang against the window mass at 18 to 21 kHz, above the channel's first cross mode (13.5 kHz across 12.9 mm).
# The windway HEIGHT is not published for this alto. The recipe uses 1.8 mm / 1.7 = 1.06 mm (the opening height of the Berlin Bressan alto e
# divided by the opening-to-exit area ratio of Bouterse's copy), with the bounds 0.5 to 1.8 mm in unknowns.tsv; the flue width is the window width.
# The viscous windway law uses that exit height over the 61.5 mm windway length (an upper bound on the resistance of a tapering windway). The jet
# card is the repository interim slot-jet card, chosen BEFORE any render because a windway jet is a slot jet (as for the tin whistle, whose
# measured register ladder it reproduces).
# Player: Terrien, Blandin, Vergez and Fabre, Acta Acustica united with Acustica 101 (2015) 300-316, artificial mouth on an alto Bressan Zen-On
# recorder, F4, A4, Bb4 and B4 fingerings: the pressure at which the first register jumps to the second on an increasing ramp, read by eye from
# their Fig. 2 (F4 230, A4 305, Bb4 505, B4 490 Pa, error about 20 Pa), minus the hysteresis of their Fig. 3 (95, 95, 90, 60 Pa): below that pressure
# only the first register exists. The held notes (recorder-<fingering>) are blown at that pressure: F4 135, A4 210, Bb4 415, B4 430 Pa. The
# crescendos (recorder-<fingering>-crescendo) are blown with NO pressure chosen: a linear ramp from 3 Pa to 600 Pa at 100 Pa/s (their ramps
# ran at up to 500 Pa/s), the protocol of that paper. The fingerings are Bouterse's f1 (0 1 2 3 4 5 6 7), a1 (0 1 2 3 4 5), b-flat1 (0 1 2 3 4 . 6),
# b1 (0 1 2 3 . 5 6). The gesture of the held notes is 50 ms attack, hold, 10.5 ms release.
# Labium offset +0.1455 mm on the pipe side, DERIVED: the edge underside sits step - h/2 = 0.8 - 0.53 = 0.27 mm on the bore side of the jet
# centreline, and the jet splits at the middle of the 0.25 mm edge face, where the dividing streamline of a jet meeting a flat face ends; the
# face ends bracket it (+0.02 to +0.27 mm) and the windway-height bounds move it from -0.225 to +0.425 mm. Jet angle zero: the 0.3 mm rise of
# the windway roof moves the jet 0.01 to 0.02 mm at the edge, inside that bracket.
# Unknown or authored: the windway height and shape, onset ramp (50 ms), mouth loss (0.02), jet
# displacement gain (90, as in the flue gates), hole undercutting, hole 7, the thumb-hole venting that gives the second and third registers (not
# modelled, so only fingerings of the first register are rendered), hole 7's position, the labium edge height. The binary must contain the
# jet-island phase fix (fs-couple 1718d79fb or later), the computed window (the commit that adds this row) and the one-sample channel line
# (fs-couple 988c2eafd or later).
# Gas: the woodwind breath law, exhaled breath saturated with water at 27 C as the mean temperature of the bore, with its published 4.5 % CO2
# and 16.5 % O2 (dry): the shared woodwind pack, row `gas pack breath-woodwind-exhaled-coltman1966-nederveen1969` (ledger
# coltman-nederveen-fletcher-woodwind-breath; it renders bit-identical to the explicit mixture this recipe carried before, 300.15 K,
# 101325 Pa, mole fractions N2 0.7534 O2 0.1591 Ar 0.0090 CO2 0.0434 H2O 0.0351) (Coltman 1966 p. 101,
# Nederveen 1969 sec. 12 p. 28, Fletcher and Rossing sec. 8.1 p. 191; the clarinet recipes' mixture, pexfp.234), because these data were
# played by mouth. Its sound speed is 345.36 m/s; against saturated air the CO2 and the oxygen depletion lower every passive mode by 19.1 to
# 19.3 cents (u-breath-co2). The window and holes radiate into the same gas (exterior same). Not in it: the temperature gradient along the
# bore; Bouterse's rooms, 15 to 27 C, bracket the temperature (u-temperature: about 40 cents). Bouterse puts this alto at about a1 = 409 Hz.
#
# Numerical controls forced by the generic route (none is a physical parameter, none was tuned): mechanics clock 288 kHz because a viscothermal loss
# cell needs four transit samples and under 0.5 rad at its band edge (3 kHz); loss only on bore sections of 17 mm or more; compact baffled radiation
# loads, which hold for ka <= 0.5 and so cap the receiver band at 2.8 kHz; minimum jet speed 0.2 m/s for the held notes and 0.1 m/s for the
# crescendos (a viscous windway jet slower than the minimum at a control tick is refused; the crescendo starts at 3 Pa for that reason); 8192
# control ticks for the 7 s files; the pressure release lasts 10.5 ms because a release that ends exactly on a control tick can leave a
# rounding-sized pressure that the jet domain check refuses.
#
# SEED=<n> re-renders with another jet seed into $OUT/seed-<n> (seed-independence check).
set -euo pipefail

SHA12=${SHA12:?set SHA12 to the pushed sha12 whose snapshot holds these files}
BIN_SHA12=${BIN_SHA12:-$SHA12}
BIN=${BIN:-/data/tmp/wa-inst/bin/$BIN_SHA12/music_render}
SNAP=${SNAP:-/data/tmp/wa-inst/snap/$SHA12/frankensim}
OUT=${OUT:-/data/tmp/wa-inst/renders/recorder/$SHA12}
NOTES=${NOTES:-"f1 a1 bb1 b1 f1-crescendo a1-crescendo bb1-crescendo b1-crescendo"}
DIR=data/instruments/descriptions/recorder

cd "$SNAP"
if [ -n "${SEED:-}" ]; then
  OUT="$OUT/seed-$SEED"
fi
mkdir -p "$OUT"
for note in $NOTES; do
  src="$DIR/recorder-$note.performance"
  if [ -n "${SEED:-}" ]; then
    src="$OUT/recorder-$note.seed.performance"
    sed -E "s/^jet ([^ ]+) ([^ ]+) [0-9]+ ([0-9]+)\$/jet \\1 \\2 $SEED \\3/" "$DIR/recorder-$note.performance" >"$src"
  fi
  nohup nice -n 10 "$BIN" wind "$src" "$OUT/recorder-$note.wav" \
    >"$OUT/recorder-$note.log" 2>&1 &
done
wait
ls -l "$OUT"/*.wav
