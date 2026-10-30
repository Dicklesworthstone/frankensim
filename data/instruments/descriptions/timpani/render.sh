#!/usr/bin/env bash
# Tp-S1 timpani render (pexfp.376) on trj: Fleischer's kettle B with a 0.2 mm PET head, struck by
# Rhaouti, Chaigne and Joly's measured mallet, tuned to the principal of VCSL Timpani 2B.
#   render.sh BIN_SHA12 [DATA_SHA12 [OUT_DIR]]
# BIN_SHA12 names a wa.sh build that produced the `kettle` example:
#   ssh trj /data/tmp/wa-inst/wa.sh build <sha> --release -p fs-couple --example kettle
# A SHORT STRIKE: 0.5 s at 48 kHz mono (24000 frames), 16-bit PCM with 20 Pa full scale (a clip
# refuses the WAV), observed 1 m off the axis and 1 m above the head plane in free field. Pitch and
# partial ratios are graded; the decay is not (ruling 8371). The exterior load's 317 poles make each
# 288 kHz step cost about 50x the mechanics (pexfp.376.2), so 0.5 s takes most of a day
# single-threaded on loaded trj. The bake's receipt (reciprocity defect, power margin, load fit
# errors) is the head of timpani-2b.log.
# Reference (grading only, never an input): VCSL Timpani 2B one-shots, principal 112.1 Hz, ids
# timpani-vcsl-timpani2b-hit-* in data/instruments/references-screened.json.
set -euo pipefail
SHA12="${1:?usage: render.sh BIN_SHA12 [DATA_SHA12 [OUT_DIR]]}"
DATA_SHA12="${2:-${SHA12}}"
BIN="${BIN:-/data/tmp/wa-inst/bin/${SHA12}/kettle}"
SNAP="/data/tmp/wa-inst/snap/${DATA_SHA12}/frankensim"
OUT="${3:-/data/tmp/wa-inst/renders/timpani/${DATA_SHA12}-${SHA12}}"
D="${SNAP}/data/instruments/descriptions/timpani"
mkdir -p "${OUT}"
cd "${OUT}"
cp "${D}/timpani-2b.fsk" .
# --wav FULL_SCALE_PA MIN_HZ MAX_HZ INTERVALS ORDER X_M Y_M Z_M (receiver in the kettle's frame,
# z up from the kettle bottom; the head plane is z = 0.53 m). The 40..400 Hz band covers every
# retained head mode's coupled frequency (the highest, (5,1), near 330 Hz).
# The WAV is published only whole: a refusal leaves timpani-2b.wav.part and prints the log's end.
if ! nice -n 10 "${BIN}" timpani-2b.fsk 24000 --wav 20 40 400 20 10 1.0 0 1.53 \
  > timpani-2b.wav.part 2> timpani-2b.log; then
  tail -n 4 timpani-2b.log >&2
  echo "render.sh: the kettle refused; see ${OUT}/timpani-2b.log" >&2
  exit 1
fi
mv timpani-2b.wav.part timpani-2b.wav
tail -n 4 timpani-2b.log
