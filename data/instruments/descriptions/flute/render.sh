#!/usr/bin/env bash
# Native renders of the committed published-Boehm flute performances (pexfp.320, lane P4).
# Run on trj. The performance files must exist in the snapshot of a pushed SHA; the music_render
# binary may come from an earlier green SHA (BIN_SHA12), it does not read these files at build time:
#   SHA12=<sha12 containing these files> BIN_SHA12=<sha12 with a release music_render> \
#     data/instruments/descriptions/flute/render.sh
# Each note is mf at 48 kHz mono, 2.5 s, observed as the summed pressure of every open aperture
# (open holes, open end and the embouchure) at a physical receiver 1 m above the embouchure.
#
# Numerical controls forced by the generic route (none is a physical parameter, none was tuned):
# mechanics clock 288 kHz because a viscothermal loss cell needs four transit samples and under 0.5 rad
# at its band edge (here 3 kHz); loss only on bore sections of 17 mm or more (the 3.3 mm chimneys and the
# 5 mm riser are lossless); compact baffled loads for every opening, which hold only for ka <= 0.5 and so
# cap the receiver band at 2.8 kHz (the computed mutual operator is capped at eight openings and fails
# its quadrature at this clock); the pressure release lasts 10.5 ms because a release that ends exactly
# on a control tick can leave a rounding-sized positive pressure that the jet domain check refuses.
#
# Geometry: data/instruments/geometry/flute-boehm-silver-bfoot-boehm1871 (Boehm 1871 as recomputed by
# Miller 1922, Table II, rescaled from A=435 to A=440 by Boehm's published Schema rule; riser of the
# Pearl PF-661 from Wolfe, Smith and Green 2003). Every geometry number is a dimensions.tsv row.
# Player: experienced flautist, mf, de la Cuadra, Fabre, Montgermont and Chafe, Acta Acustica 94(5),
# 2008, pp.740-749. Mean mouth pressure 300 Pa (D4) to 600 Pa (A5), jet length W 8 mm to 4.5 mm
# (Sec. 5.1); lip-hole area ~6e-6 to ~3.5e-6 m2 (Fig. 4) and height ~0.9 mm (Fig. 6) are read off
# the plots; values between the two published endpoints are linear in log2(f). Never read from a
# recording. Labium offset +0.29 mm on the pipe side: one player's measured reference jet (Onogi,
# Yokoyama and Iida, Acta Acustica 5 (2021) 11, Table 1: y/h = 0.24 +- 0.25 at h = 1.2 mm), so the
# bracket is -0.01 to +0.59 mm; jet angle 0, because the measured 39 degrees are inside that offset.
# Onset ramp (50 ms), mouth loss (0.02), the jet displacement gain (90) are UNKNOWN or authored: see
# unknowns.tsv and status.json. The jet card is an inline X-Est
# lip-jet card, not the interim 2-D slot-jet card: theta/b 1.75 places the island's band peak (0.55 omega_c,
# omega_c = U / 4 theta) at the published Strouhal number of maximum jet amplification, Str_h = 0.025 (de la
# Cuadra 2008 Sec. 7), and the convection ratio is 0.4 (Auvray 2012; Terrien 2014; de la Cuadra: 0.3 to 0.5).
# The binary must contain the jet-island phase fix (fs-couple 1718d79fb or later).
#
# References (grading only, never an input): UIowa MIS Flute.nonvib.mf.B3B4 (D4, G4) and
# Flute.nonvib.mf.C5B5 (A5), ids uiowa-flute-flute-nonvib-mf-b3b4 / -c5b5 in references.json.
#
# SEED=<n> re-renders with another jet seed into $OUT/seed-<n> (seed-independence check).
set -euo pipefail

SHA12=${SHA12:?set SHA12 to the pushed sha12 whose snapshot holds these files}
BIN_SHA12=${BIN_SHA12:-$SHA12}
BIN=${BIN:-/data/tmp/wa-inst/bin/$BIN_SHA12/music_render}
SNAP=${SNAP:-/data/tmp/wa-inst/snap/$SHA12/frankensim}
OUT=${OUT:-/data/tmp/wa-inst/renders/flute/$SHA12}
NOTES=${NOTES:-"d4 g4 a5"}
DIR=data/instruments/descriptions/flute

cd "$SNAP"
if [ -n "${SEED:-}" ]; then
  OUT="$OUT/seed-$SEED"
fi
mkdir -p "$OUT"
for note in $NOTES; do
  upper=$(printf '%s' "$note" | tr 'a-z' 'A-Z')
  src="$DIR/flute-$note.performance"
  if [ -n "${SEED:-}" ]; then
    src="$OUT/flute-$note.seed.performance"
    sed -E "s/^jet ([^ ]+) ([^ ]+) [0-9]+ ([0-9]+)\$/jet \\1 \\2 $SEED \\3/" "$DIR/flute-$note.performance" >"$src"
  fi
  nohup nice -n 10 "$BIN" wind "$src" "$OUT/Flute.nonvib.mf.$upper.wav" \
    >"$OUT/Flute.nonvib.mf.$upper.log" 2>&1 &
done
wait
ls -l "$OUT"/*.wav
