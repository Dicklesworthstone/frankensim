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
# Gas: the tubes hold the gas Coltman 1968 measured in a played flute's bore (the cards' explicit `gas mixture` row, 26.15 C saturated, 2.5 percent dry
# CO2; ledger coltman-1968-flute-bore-gas), applied to stopped pipes by analogy; it sounds every tube 16.6 to 16.7 cents above the dry air the lengths
# were cut for (unknowns u-bore-gas, whose bracket runs down to that dry air). The observer and radiation exterior stay dry 20 C air.
# Player: no measurement of a pan-pipe player's pressure or jet length is in the library. Blowing pressure and jet length at each tube's equal-tempered
# pitch are Fletcher 1975's median relations measured on four flute players (J. Acoust. Soc. Am. 57(1), pp. 233-237, Eq. 1 and Eq. 2): p [mbar] =
# 0.008 nu [Hz] and l [mm] = 1.8 + 100 nu^(-1/2), i.e. 117.5 Pa and 10.05 mm at D3 to 395.1 Pa and 6.30 mm at B4, the card pressure printed to five
# significant figures; Fletcher 2005 (Sec. IV and VI) points to that paper for the jet speed "as in the flute" and says the jet is shorter and the
# pressure higher for high notes, so the step from the flute to a stopped pipe is an analogy, and D3 and G3 lie at or below the lowest pitches the
# source reports (unknowns u-player-controls, ledger fletcher-1975-flute-blowing-pressure-and-jet-length). Initial jet half-width 1 mm to 0.5 mm
# (lip height 2 b0 reduced for high notes, Fletcher 2005 Sec. VI) is still the published range ends paired linearly in log2(f) between the lowest
# and highest tube of the 13-note set. The cards used before took the range ends of Fletcher 2005 for all three controls (V0 20 to 40 m/s, jet length
# 10 mm to 5 mm), paired the same way. The labium offset is the cited stopped-pipe value of 0.5 jet half-widths
# (Eq. 6, alpha about 2, gamma about 0.4): the steady entrained flow is blocked by the closed end, so the jet midplane passes outside the lip; it is
# applied to the island's constant half-width, because the island does not model the jet spreading (u-jet-spreading). The lip aperture is as wide as the
# tube (the jet cannot be wider; u-lip-aperture-width). The gesture is 50 ms attack, hold, 10.5 ms release.
# Unknown or authored: onset ramp (50 ms), mouth loss (0.02), jet displacement gain (90, as in the flue gates), the jet card (the inline X-Est lip-jet card
# of the flute, chosen BEFORE looking at any render because the pan-pipe jet comes from the lips). The jet length is 5.03 (D3), 5.08 (G3), 5.34 (D4), 5.76 (G4)
# and 6.30 (B4) times its height; the island's early-tap correction needs 5.6 for that card, so for D3, G3 and D4 the loop delay is too long (about 12 percent at a ratio of 5, not re-derived at 5.08 and 5.34). The binary must contain
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
