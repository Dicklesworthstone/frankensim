#!/bin/zsh
# Build the live piano for the browser and install it under
# apps/instrument-lab/public/live/piano/ as engine.js + engine_bg.wasm.
# The crate links fs-couple and its neighbours by path, so the wasm holds whatever those files say:
# build from an export of a commit (git archive, siblings beside it) when the shared tree carries other
# work in progress, and the wasm is that commit's.
set -eu
HERE=${0:A:h}
TARGET=${FSIM_LIVE_PIANO_TARGET:-/var/folders/vt/n2xyn_s51b97_j3yh2qbqcnc0000gn/T/rch_target_frankensim_live_piano_wasm}
OUT=$HERE/../public/live/piano
cd $HERE
# simd128 (wasm SIMD: Chrome 91, Firefox 89, Safari 16.4): the kernel's mode loops are written as f64x2.
SIMD="-C target-feature=+simd128"
run() { env RCH_SHIM_LOCAL_IDE=1 CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-8} CARGO_TARGET_DIR=$TARGET CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_RUSTFLAGS=$SIMD cargo "$@"; }
STAGE=$TARGET/bindgen-stage
mkdir -p $STAGE $OUT
run build --release --lib --target wasm32-unknown-unknown
wasm-bindgen --target web --out-dir $STAGE --out-name engine $TARGET/wasm32-unknown-unknown/release/fs_live_piano.wasm
wasm-opt -O3 --enable-simd --enable-bulk-memory --enable-nontrapping-float-to-int --enable-sign-ext --enable-mutable-globals \
  $STAGE/engine_bg.wasm -o $STAGE/engine_bg.opt.wasm
# RUSTFLAGS in the environment replaces the target flags above and would drop simd128 without a word
if command -v wasm-objdump > /dev/null && ! wasm-objdump -d $STAGE/engine_bg.opt.wasm | grep -q 'f64x2\.'; then
  echo "error: no f64x2 instructions in the piano wasm; was +simd128 dropped (RUSTFLAGS set?)" >&2
  exit 1
fi
cp $STAGE/engine.js $OUT/engine.js
cp $STAGE/engine.d.ts $OUT/engine.d.ts
cp $STAGE/engine_bg.wasm.d.ts $OUT/engine_bg.wasm.d.ts
cp $STAGE/engine_bg.opt.wasm $OUT/engine_bg.wasm
echo "engine_bg.wasm raw=$(wc -c < $OUT/engine_bg.wasm) gzip=$(gzip -9 -c $OUT/engine_bg.wasm | wc -c)"
