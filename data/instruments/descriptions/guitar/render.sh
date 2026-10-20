#!/usr/bin/env bash
# Gt-S1 classical-guitar open-string renders (trj). Plain nylon trebles only; the wound basses are
# blocked on UNKNOWN core/winding geometry (see sources.tsv).
#   render.sh SHA12 [OUT_DIR]
# SHA12 names a wa.sh build that produced both `xtask` and the `material_strings` example:
#   ssh trj /data/tmp/wa-inst/wa.sh build <sha> --release -p fs-couple -p xtask --example material_strings --bin xtask
# Reference notes in the Iowa library (48 kHz mono, mf):
#   gt-s1-g3 -> Guitar.mf.sulG.G3B3 @ G3     gt-s1-b3 -> Guitar.mf.sulB.B3 @ B3     gt-s1-e4 -> Guitar.mf.sul_E.E4B4 @ E4
set -euo pipefail
SHA12="${1:?usage: render.sh SHA12 [OUT_DIR]}"
BIN="/data/tmp/wa-inst/bin/${SHA12}"
SNAP="/data/tmp/wa-inst/snap/${SHA12}/frankensim"
OUT="${2:-/data/tmp/wa-inst/renders/guitar/${SHA12}}"
D="${SNAP}/data/instruments/descriptions/guitar"
mkdir -p "${OUT}"
cd "${OUT}"
for pack in nylon-g3 nylon-b3 nylon-e4 cedar-top; do
  cp "${D}/${pack}.tsv" "${D}/${pack}-manifest.tsv" .
  (cd "${SNAP}" && "${BIN}/xtask" matdb-pack --manifest "${OUT}/${pack}-manifest.tsv" --out "${OUT}/${pack}.fsmatpk")
done
cp "${D}"/gt-s1-*.material-strings "${D}/gt-s1-top.body" .
for note in g3 b3 e4; do
  nice -n 10 "${BIN}/material_strings" "gt-s1-${note}.material-strings" "${OUT}/gt-s1-${note}.wav" \
    --body-profile gt-s1-top.body 2>&1 | tee "gt-s1-${note}.log"
done
