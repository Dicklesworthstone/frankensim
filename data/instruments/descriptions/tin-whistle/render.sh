#!/usr/bin/env bash
# Native renders of the committed tin-whistle performances (pexfp.338, lane P4). Listening only: there is no licensed reference recording yet
# (bead pexfp.219), so these are not graded. Run on trj. The performance files must exist in the snapshot of a pushed SHA; the music_render
# binary may come from an earlier green SHA (BIN_SHA12), it does not read these files at build time:
#   SHA12=<sha12 containing these files> BIN_SHA12=<sha12 with a release music_render> \
#     data/instruments/descriptions/tin-whistle/render.sh
# Each note is 48 kHz mono, 2.5 s, observed as the summed pressure of every open aperture (open holes, open end and the window) at a physical
# receiver 1 m above the window. Every note is rendered at both ends of the labium edge bracket, into $OUT/prism and $OUT/wedge
# (ENDS=prism or ENDS=wedge renders one).
#
# Geometry: data/instruments/geometry/tinwhistle-clare-d-davies2006 (a Clare D whistle measured by P. O. A. L. Davies, J. Pinho and E. J. English,
# Proc. Inst. Acoust. 28(1), 2006, pp. 597-606; fingerings from Timoney et al., DAFx'04 Fig. 2 and the paper's Sec. 3.2). Every geometry number
# is a dimensions.tsv row. The network is: the window (2.4 mm wall as a riser of equal area to the 8 x 5.3 mm opening) on a junction with a
# closed pocket behind it (8 mm beyond the flue plus half the window, a cylinder of the bore radius: an upper bound, u-pocket-shape), the bore to
# hole 6 and on to hole 1 and the open end, each finger hole a chimney section ending open or closed. The window's mass is computed, not a
# constant: each card's row `window wall-cut 8 cheeks parallel` solves the 3-D potential flow through the 8 mm wide by 5.3 mm cut in the
# 2.4 mm wall of the 6.1 mm bore, its side walls parallel and 8 mm apart like the windway's, with the pocket behind it and free space
# outside (fs_bem::wall_window, a staircase of cut cells at 8, 16 and 32 axial cells, extrapolated at the order the three rungs show).
# Where the wedge's tip sits inside the wall is not published (u-labium-edge-height): the paper gives its 25 degree angle and 11 mm upper
# surface only. So the edge height is a stated bracket and both ends are rendered, neither preferred: prism, the card's face straight
# through the wall (17.84 mm of bore, order 1.20), and wedge, the 25 degree face rising from a tip 0.5 mm above the
# bore (`ramp 0.0005 25` added to the row: 16.31 mm, order 1.51). That mass replaces the 0.61 a_eq mouth inertance; the
# jet's source coefficient keeps its own law. The binary must contain the parallel-cheek window (the commit that adds this row).
# The window radiates; its loss is not authored. Each card's row `mouth_reaction window-cylinder` takes the flux that the same solve puts
# through the cut's outer face (8.5 mm from the axis, on the head's rigid outer cylinder, the solve's own exterior) and loads the mouth with
# that face's radiation: resistance from the cylindrical-harmonic series of Junger and Feit (propagating axial wavenumbers only), static mass
# from it by Kramers-Kronig, realized as a passive chain of parallel R-L-C cells fitted from 20 Hz to 24 kHz (prism: 12 cells; error 1.4 % in
# that band, 3.6 % after the 288 kHz clock's warp, 4.8 % up to its 144 kHz Nyquist; wedge: 11 cells, 1.7 %, 3.3 % and 5.6 %; kA 0.003 to 3.7
# on both). The radiation's static mass (prism 65.9, wedge 55.4 Pa s^2/m^3) comes out of the computed window mass, so the outer mass is not
# counted twice. It replaces mouth_loss_rel 0.02: the flue
# row's last field is 0. The uniform cut-up piston on the same cylinder (`mouth_reaction window-cylinder-piston`, static mass 70.3) is the
# shape bracket, reported with the landing renders, not rendered here. The binary must contain the window radiation (the commit that adds
# `mouth_reaction window-cylinder`).
# Player: the air flow the same paper measured for each note of the major scale at its tempered frequency (Table 1: l/min D5 6.2, E5 6.9,
# F#5 9.0, G5 8.2, A5 10.0, B5 9.55, C#6 11.3, D6 with only hole 6 open 11.9; and the paper's overblown series with every hole closed, Fig. 2 and
# Table 1: D6 12.0, A6 17.9, D7 23.8, files d6c, a6c, d7c, a test of the overblowing ladder whose answer is a result), turned into mouth pressure by
# the windway law the cards name, `windway slit-developing 0.02455` (pexfp.135.4; the binary must contain e2414832f): Shah's (1978) laminar
# developing-flow law for a straight slit fed through a uniform inlet, p = rho v^2/2 + mu L v (f_app Re)/(2 h^2) with
# f_app Re = 3.44/sqrt(x+) + (24 + 0.674/(4 x+) - 3.44/sqrt(x+))/(1 + 0.000029/x+^2) and x+ = mu L/(4 rho v h^2) (form and constants as McCreery
# et al., NURETH-15 2013, eq. 2, state them), L = 24.55 mm, h = 1.48 mm (exit depth over the whole length: the straight stand-in most likely
# overstates the loss, u-flue-taper), flue width 8 mm, v = flow / (8 x 1.48 mm^2), dry air at 20 C (rho 1.2041 kg/m^3, mu 1.8134e-5 Pa s, the
# gas model's own). Pressures D5 to D7c in card order d5 e5 fs5 g5 a5 b5 cs6 d6 d6c a6c d7c: 93.447, 112.43, 178.71, 151.84, 215.06, 198.32,
# 266.83, 292.44, 296.81, 606.64, 1016.5 Pa, 1.34 to 1.40 times those of the fully developed law 12 mu L v/h^2 + rho v^2/2 the cards used
# before; the voice holds the measured flow either way. Each held note is inside the law's cited validity, which the binary checks and the
# sidecar records (windway_held): x+ 1.26e-3 to 4.84e-3 (floor 4e-4), Re_Dh 1715 to 6585, Re_L 2.0e4 to 6.7e4, the core still inviscid at the
# exit with Re_delta* 183 to 373 (bound 520), Chen's low-Re inlet term moving the drop 0.6 to 1.9 % (bound 2.4 %). The attack and release
# use the same quasi-steady law and are not covered. Cross-code check of the law against Segoufin's laminar RNS-p straight channel
# (thesis Fig. 5.14; laminar vs laminar, it verifies the law in x+ only, laminarity is bounded separately): the thesis states two Reynolds
# length scales, so the status is split, MET under Re_Dh = 2 Re (2/2) and MISSED under Re_Dh = 4 Re (1/3). Reading II misses by 7 to 9 %
# at x+ 1.6e-3 to 3.3e-3, which is a6c's x+ (1.68e-3); d7c (x+ 1.26e-3) lies below the graded range. Chen's 0.64 + 38/Re and Shah and
# London's agreement statement behind the x+ floor are secondary, search snippet only; primary not read. The four *-nolle gate files keep
# `windway slit` and its pressures: the Nolle/Price transfer refuses this windway. The gesture is 50 ms attack, hold, 10.5 ms release.
# Dry air at 20 C is the measurement's own gas, not an omission: the paper blew the whistle from a regulated, metered air supply in an anechoic
# room at 20 C (specimen.tsv, environment), not by mouth. Do not change it to breath; breath gas belongs to recipes whose data were played.
# Unknown or authored: labium offset and jet angle (zero, UNKNOWN), onset ramp (50 ms), jet displacement gain (90, as in the
# flue gates), the jet card (the repository interim slot-jet card, chosen BEFORE looking at any render because a windway jet is a slot jet; it is
# outside its Reynolds band, u-jet-profile). The binary must contain the jet-island phase fix (fs-couple 1718d79fb or later).
#
# Numerical controls forced by the generic route (none is a physical parameter, none was tuned): mechanics clock 288 kHz because a viscothermal loss
# cell needs four transit samples and under 0.5 rad at its band edge (3 kHz); loss only on bore sections of 17 mm or more; finger holes are chimney
# sections of 1.5 mm (a section under one propagation step is refused, the wall is 0.3 mm: u-hole-chimney); compact baffled radiation loads, which
# hold for ka <= 0.5 and so cap the receiver band at 2.8 kHz; minimum jet speed 0.2 m/s (the windway gives a jet slower than 1 m/s at the
# first control tick of the attack, which the default limit refuses); the pressure release lasts 10.5 ms because a release that ends exactly on a
# control tick can leave a rounding-sized pressure that the jet domain check refuses.
#
# SEED=<n> re-renders with another jet seed into $OUT/seed-<n>/<end> (seed-independence check).
set -euo pipefail

SHA12=${SHA12:?set SHA12 to the pushed sha12 whose snapshot holds these files}
BIN_SHA12=${BIN_SHA12:-$SHA12}
BIN=${BIN:-/data/tmp/wa-inst/bin/$BIN_SHA12/music_render}
SNAP=${SNAP:-/data/tmp/wa-inst/snap/$SHA12/frankensim}
OUT=${OUT:-/data/tmp/wa-inst/renders/tin-whistle/$SHA12}
NOTES=${NOTES:-"d5 e5 fs5 g5 a5 b5 cs6 d6 d6c a6c d7c"}
ENDS=${ENDS:-"prism wedge"}
DIR=data/instruments/descriptions/tin-whistle

ROW="window wall-cut 8 cheeks parallel"
cd "$SNAP"
seed=()
if [ -n "${SEED:-}" ]; then
  OUT="$OUT/seed-$SEED"
  seed=(-e "s/^jet ([^ ]+) ([^ ]+) [0-9]+ ([0-9]+)\$/jet \\1 \\2 $SEED \\3/")
fi
for end in $ENDS; do
  case "$end" in
    prism) row="$ROW" ;;
    wedge) row="$ROW ramp 0.0005 25" ;;
    *) echo "ENDS takes prism and wedge, not $end" >&2; exit 64 ;;
  esac
  mkdir -p "$OUT/$end"
  for note in $NOTES; do
    src="$OUT/$end/tin-whistle-$note.performance"
    sed -E -e "s/^$ROW\$/$row/" ${seed[@]+"${seed[@]}"} "$DIR/tin-whistle-$note.performance" >"$src"
    grep -qx "$row" "$src" || { echo "$src lacks the row: $row" >&2; exit 65; }
    nohup nice -n 10 "$BIN" wind "$src" "$OUT/$end/tin-whistle-$note.wav" \
      >"$OUT/$end/tin-whistle-$note.log" 2>&1 &
  done
done
wait
ls -l "$OUT"/*/*.wav
