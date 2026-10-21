#!/usr/bin/env bash
# Gt-S1 classical-guitar open-string renders (trj). Plain nylon trebles only; the wound basses are
# blocked on UNKNOWN core/winding geometry (see sources.tsv).
#   render.sh BIN_SHA12 [DATA_SHA12 [OUT_DIR]]
# BIN_SHA12 names a wa.sh build that produced both `xtask` and the `material_strings` example
# (data is read from the snapshot DATA_SHA12, default BIN_SHA12; pass a later snapshot when the data landed after the build):
#   ssh trj /data/tmp/wa-inst/wa.sh build <sha> --release -p fs-couple -p xtask --example material_strings --bin xtask
# Reference notes in the Iowa library (48 kHz mono, mf):
#   gt-s1-solo-g3 -> Guitar.mf.sulG.G3B3 @ G3   gt-s1-solo-b3 -> Guitar.mf.sulB.B3 @ B3   gt-s1-solo-e4 -> Guitar.mf.sul_E.E4B4 @ E4
# The solo files (one plucked string, 12 modes per plane) fit the ordinary 42-mode connected exact block. The three-string
# gt-s1-{g3,b3,e4}.material-strings files (84 string modes) need the prepared propagator at bank construction, which the
# generic bank does not yet offer: the bank is built with the ordinary owner before prepared-linear is applied.
set -euo pipefail
SHA12="${1:?usage: render.sh BIN_SHA12 [DATA_SHA12 [OUT_DIR]]}"
DATA_SHA12="${2:-${SHA12}}"
BIN="/data/tmp/wa-inst/bin/${SHA12}"
XTASK="${XTASK:-${BIN}/xtask}"   # packs only need any xtask that predates a pack-format change
SNAP="/data/tmp/wa-inst/snap/${DATA_SHA12}/frankensim"
OUT="${3:-/data/tmp/wa-inst/renders/guitar/${DATA_SHA12}-${SHA12}}"
D="${SNAP}/data/instruments/descriptions/guitar"
mkdir -p "${OUT}"
cd "${OUT}"
for pack in nylon-g3 nylon-b3 nylon-e4 cedar-top; do
  cp "${D}/${pack}.tsv" "${D}/${pack}-manifest.tsv" .
  (cd "${SNAP}" && "${XTASK}" matdb-pack --manifest "${OUT}/${pack}-manifest.tsv" --out "${OUT}/${pack}.fsmatpk")
done
cp "${D}"/gt-s1-solo-*.material-strings "${D}"/gt-s1-solo-*.body "${D}"/gt-s1-fret-*.material-strings .
# Independent single-core renders (about 5 min each). Open strings: g3 b3 e4. Stopped notes on the B string
# (gb4 = fret 7, f5 = fret 18) reuse the B-string pack, tension and body profile unchanged.
# NOTES="gb4 f5" renders only some of them.
for note in ${NOTES:-g3 b3 e4 gb4 f5}; do
  case "${note}" in
    gb4|f5) stem="gt-s1-fret-${note}"; body="gt-s1-solo-b3.body" ;;
    *) stem="gt-s1-solo-${note}"; body="gt-s1-solo-${note}.body" ;;
  esac
  nice -n 10 "${BIN}/material_strings" "${stem}.material-strings" "${OUT}/${stem}.wav" \
    --body-profile "${body}" > "${stem}.log" 2>&1 &
done
wait
tail -n 2 gt-s1-*.log
