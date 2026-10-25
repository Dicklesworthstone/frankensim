#!/usr/bin/env bash
# Native renders of the committed tin-whistle performances (pexfp.338, lane P4). Listening only: there is no licensed reference recording yet
# (bead pexfp.219), so these are not graded. Run on trj. The performance files must exist in the snapshot of a pushed SHA; the music_render
# binary may come from an earlier green SHA (BIN_SHA12), it does not read these files at build time:
#   SHA12=<sha12 containing these files> BIN_SHA12=<sha12 with a release music_render> \
#     data/instruments/descriptions/tin-whistle/render.sh
# Each note is 48 kHz mono, 2.5 s, observed as the summed pressure of every open aperture (open holes, open end and the window) at a physical
# receiver 1 m above the window.
#
# Geometry: data/instruments/geometry/tinwhistle-clare-d-davies2006 (a Clare D whistle measured by P. O. A. L. Davies, J. Pinho and E. J. English,
# Proc. Inst. Acoust. 28(1), 2006, pp. 597-606; fingerings from Timoney et al., DAFx'04 Fig. 2 and the paper's Sec. 3.2). Every geometry number
# is a dimensions.tsv row. The network is: the window (2.4 mm wall as a riser of equal area to the 8 x 5.3 mm opening) on a junction with a
# closed pocket behind it (8 mm beyond the flue plus half the window, a cylinder of the bore radius: an upper bound, u-pocket-shape), the bore to
# hole 6 and on to hole 1 and the open end, each finger hole a chimney section ending open or closed. The window's mass is computed, not a
# constant: each card's row `window wall-cut 16` solves the 3-D potential flow through the 8 mm arc by 5.3 mm cut in the 2.4 mm wall of the
# 6.1 mm bore, with the pocket behind it and free space outside (fs_bem::wall_window at 16 and 32 axial cells, extrapolated at the re-entrant
# edges' order 4/3). That mass, 15.71 mm of bore, replaces the 0.61 a_eq mouth inertance; the jet's source coefficient keeps its own law.
# Player: the air flow the same paper measured for each note of the major scale at its tempered frequency (Table 1: l/min D5 6.2, E5 6.9,
# F#5 9.0, G5 8.2, A5 10.0, B5 9.55, C#6 11.3, D6 with only hole 6 open 11.9; and the paper's overblown series with every hole closed, Fig. 2 and
# Table 1: D6 12.0, A6 17.9, D7 23.8, files d6c, a6c, d7c, a test of the overblowing ladder whose answer is a result), turned into mouth pressure by the island's own viscous windway law
# p = 12 mu L v / h^2 + rho v^2 / 2 with L = 24.55 mm, h = 1.48 mm (exit depth over the whole length: an upper bound, u-flue-taper), flue width 8 mm,
# v = flow / (8 x 1.48 mm^2), dry air at 20 C. The gesture is 50 ms attack, hold, 10.5 ms release.
# Dry air at 20 C is the measurement's own gas, not an omission: the paper blew the whistle from a regulated, metered air supply in an anechoic
# room at 20 C (specimen.tsv, environment), not by mouth. Do not change it to breath; breath gas belongs to recipes whose data were played.
# Unknown or authored: labium offset and jet angle (zero, UNKNOWN), onset ramp (50 ms), mouth loss (0.02), jet displacement gain (90, as in the
# flue gates), the jet card (the repository interim slot-jet card, chosen BEFORE looking at any render because a windway jet is a slot jet; it is
# outside its Reynolds band, u-jet-profile). The binary must contain the jet-island phase fix (fs-couple 1718d79fb or later).
#
# Numerical controls forced by the generic route (none is a physical parameter, none was tuned): mechanics clock 288 kHz because a viscothermal loss
# cell needs four transit samples and under 0.5 rad at its band edge (3 kHz); loss only on bore sections of 17 mm or more; finger holes are chimney
# sections of 1.5 mm (a section under one propagation step is refused, the wall is 0.3 mm: u-hole-chimney); compact baffled radiation loads, which
# hold for ka <= 0.5 and so cap the receiver band at 2.8 kHz; minimum jet speed 0.2 m/s (the viscous windway gives a jet slower than 1 m/s at the
# first control tick of the attack, which the default limit refuses); the pressure release lasts 10.5 ms because a release that ends exactly on a
# control tick can leave a rounding-sized pressure that the jet domain check refuses.
#
# SEED=<n> re-renders with another jet seed into $OUT/seed-<n> (seed-independence check).
set -euo pipefail

SHA12=${SHA12:?set SHA12 to the pushed sha12 whose snapshot holds these files}
BIN_SHA12=${BIN_SHA12:-$SHA12}
BIN=${BIN:-/data/tmp/wa-inst/bin/$BIN_SHA12/music_render}
SNAP=${SNAP:-/data/tmp/wa-inst/snap/$SHA12/frankensim}
OUT=${OUT:-/data/tmp/wa-inst/renders/tin-whistle/$SHA12}
NOTES=${NOTES:-"d5 e5 fs5 g5 a5 b5 cs6 d6 d6c a6c d7c"}
DIR=data/instruments/descriptions/tin-whistle

cd "$SNAP"
if [ -n "${SEED:-}" ]; then
  OUT="$OUT/seed-$SEED"
fi
mkdir -p "$OUT"
for note in $NOTES; do
  src="$DIR/tin-whistle-$note.performance"
  if [ -n "${SEED:-}" ]; then
    src="$OUT/tin-whistle-$note.seed.performance"
    sed -E "s/^jet ([^ ]+) ([^ ]+) [0-9]+ ([0-9]+)\$/jet \\1 \\2 $SEED \\3/" "$DIR/tin-whistle-$note.performance" >"$src"
  fi
  nohup nice -n 10 "$BIN" wind "$src" "$OUT/tin-whistle-$note.wav" \
    >"$OUT/tin-whistle-$note.log" 2>&1 &
done
wait
ls -l "$OUT"/*.wav
