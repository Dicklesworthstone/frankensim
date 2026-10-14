#!/bin/zsh
# Build the live drum kit for the browser and install it under
# apps/instrument-lab/public/live/kit/ as engine.js + engine_bg.wasm.
# fs-vkplate is std only, so this builds straight from the working tree.
set -eu
HERE=${0:A:h}
TARGET=${FSIM_LIVE_KIT_TARGET:-/var/folders/vt/n2xyn_s51b97_j3yh2qbqcnc0000gn/T/rch_target_frankensim_live_kit}
OUT=$HERE/../public/live/kit
cd $HERE
run() { env RCH_SHIM_LOCAL_IDE=1 CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-8} CARGO_TARGET_DIR=$TARGET cargo "$@"; }
STAGE=$TARGET/bindgen-stage
mkdir -p $STAGE $OUT
run build --release --lib --target wasm32-unknown-unknown
wasm-bindgen --target web --out-dir $STAGE --out-name engine $TARGET/wasm32-unknown-unknown/release/fs_live_kit.wasm
wasm-opt -O3 --enable-bulk-memory --enable-nontrapping-float-to-int --enable-sign-ext --enable-mutable-globals \
  --enable-reference-types --enable-multivalue $STAGE/engine_bg.wasm -o $STAGE/engine_bg.opt.wasm
cp $STAGE/engine.js $OUT/engine.js
cp $STAGE/engine.d.ts $OUT/engine.d.ts
cp $STAGE/engine_bg.opt.wasm $OUT/engine_bg.wasm
echo "engine_bg.wasm raw=$(wc -c < $OUT/engine_bg.wasm) gzip=$(gzip -9 -c $OUT/engine_bg.wasm | wc -c)"
