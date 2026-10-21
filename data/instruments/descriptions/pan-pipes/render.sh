#!/usr/bin/env bash
# Native renders of the committed Andean pan-pipe performances (pexfp.345, lane P4). Listening only: there is no licensed reference recording yet
# (bead pexfp.219), so these are not graded. Run on trj. The performance files must exist in the snapshot of a pushed SHA; the music_render
# binary may come from an earlier green SHA (BIN_SHA12), it does not read these files at build time:
#   SHA12=<sha12 containing these files> BIN_SHA12=<sha12 with a release music_render> \
#     data/instruments/descriptions/pan-pipes/render.sh
# Each note is 48 kHz mono, 2.5 s, one stopped tube blown at a time, observed as the radiated pressure of the open top (the window) at a
# physical receiver 1 m above it.
#
# Geometry: data/instruments/geometry/panpipes-zampona-gmajor-fletcher2005. There is NO measured instrument: no published table of tube lengths of a
# siku or zampona was reachable. The note set is the theoretical table of the zampona of Sanchez Barrero and Ghezzi, Rev. Boliviana de Fisica 17
# (2010) pp. 47-52 (G major, 13 tubes in two rows; the five notes D3, G3, D4, G4, B4 are rendered), the bore of each tube follows N. H. Fletcher,
# J. Acoust. Soc. Am. 117(1), 2005, pp. 370-374, Eq. 1 (d = 1.32 L^0.43, mm, Eastern European panpipes), and the length of each tube is computed
# from the stopped-pipe formula L = c / (4 f) - 0.6133 a (open-end correction of a thin-walled pipe, Davies, Pinho and English 2006, Eq. 5), with
# dry air at 20 C. Those derived rows are labelled DERIVED in fingerings.tsv and unknowns.tsv (u-tube-lengths).
# Player: Fletcher 2005 Sec. IV and VI: jet speed 20 to 40 m/s (mouth pressure = rho V0^2 / 2) and a higher pressure for high notes; initial jet
# half-width 1 mm to 0.5 mm (lip height 2 b0 reduced for high notes); jet length 10 mm to 5 mm (shorter for high notes); the three controls are
# linear in log2(f) between the lowest and highest tube of the 13-note set. The labium offset is the cited stopped-pipe value of 0.5 jet half-widths
# (Eq. 6, alpha about 2, gamma about 0.4): the steady entrained flow is blocked by the closed end, so the jet midplane passes outside the lip; it is
# applied to the island's constant half-width, because the island does not model the jet spreading (u-jet-spreading). The lip aperture is as wide as the
# tube (the jet cannot be wider; u-lip-aperture-width). The gesture is 50 ms attack, hold, 10.5 ms release.
# Unknown or authored: onset ramp (50 ms), mouth loss (0.02), jet displacement gain (90, as in the flue gates), the jet card (the inline X-Est lip-jet card
# of the flute, chosen BEFORE looking at any render because the pan-pipe jet comes from the lips). The jet length is only 5 times its height, just
# under the 5.6 the island's early-tap correction needs for that card, so the loop delay is about 12 percent too long. The binary must contain
# the jet-island phase fix (fs-couple 1718d79fb or later).
#
# Numerical controls forced by the generic route (none is a physical parameter, none was tuned): mechanics clock 288 kHz because a viscothermal loss
# cell needs four transit samples and under 0.5 rad at its band edge (3 kHz); the compact baffled load holds for ka <= 0.5 and caps the receiver band
# at 2.5 kHz for the widest tube; the input requires one radiating outlet and a closed end listed as the outlet becomes an open radiating end, so the
# network is the open top (inlet), a junction, the closed tube and a 0.1 mm radius, 2 mm long open vent (area 1e-4 of the bore) as the required
# outlet, with the 2 mm riser at the top and the 2 mm vent each at least one propagation step of the clock; the pressure release lasts 10.5 ms because a release that ends exactly
# on a control tick can leave a rounding-sized pressure that the jet domain check refuses.
#
# SEED=<n> re-renders with another jet seed into $OUT/seed-<n> (seed-independence check).
set -euo pipefail

SHA12=${SHA12:?set SHA12 to the pushed sha12 whose snapshot holds these files}
BIN_SHA12=${BIN_SHA12:-$SHA12}
BIN=${BIN:-/data/tmp/wa-inst/bin/$BIN_SHA12/music_render}
SNAP=${SNAP:-/data/tmp/wa-inst/snap/$SHA12/frankensim}
OUT=${OUT:-/data/tmp/wa-inst/renders/pan-pipes/$SHA12}
NOTES=${NOTES:-"d3 g3 d4 g4 b4"}
DIR=data/instruments/descriptions/pan-pipes

cd "$SNAP"
if [ -n "${SEED:-}" ]; then
  OUT="$OUT/seed-$SEED"
fi
mkdir -p "$OUT"
for note in $NOTES; do
  src="$DIR/pan-pipes-$note.performance"
  if [ -n "${SEED:-}" ]; then
    src="$OUT/pan-pipes-$note.seed.performance"
    sed -E "s/^jet ([^ ]+) ([^ ]+) [0-9]+ ([0-9]+)\$/jet \\1 \\2 $SEED \\3/" "$DIR/pan-pipes-$note.performance" >"$src"
  fi
  nohup nice -n 10 "$BIN" wind "$src" "$OUT/pan-pipes-$note.wav" \
    >"$OUT/pan-pipes-$note.log" 2>&1 &
done
wait
ls -l "$OUT"/*.wav
