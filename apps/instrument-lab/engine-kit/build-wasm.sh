#!/bin/zsh
# Build the live drum kit for the browser and install it under
# apps/instrument-lab/public/live/kit/ as engine.js + engine_bg.wasm.
# fs-vkplate uses the shared material owner; its path dependencies and sibling
# constellation checkouts must be available, as for the live drum and piano.
set -eu
HERE=${0:A:h}
TARGET=${FSIM_LIVE_KIT_TARGET:-/var/folders/vt/n2xyn_s51b97_j3yh2qbqcnc0000gn/T/rch_target_frankensim_live_kit}
OUT=$HERE/../public/live/kit
cd $HERE
# simd128 (wasm SIMD: Chrome 91, Firefox 89, Safari 16.4; the live piano already requires it). LLVM turns
# fs-vkplate's fixed-order loops (Band::solve's 8-accumulator dot and its axpy, the elementwise force and
# predict loops) into f64x2 without reordering any sum, so the samples are bit-identical to a scalar build.
SIMD="-C target-feature=+simd128"
run() { env RCH_SHIM_LOCAL_IDE=1 CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-8} CARGO_TARGET_DIR=$TARGET CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_RUSTFLAGS=$SIMD cargo "$@"; }
STAGE=$TARGET/bindgen-stage
mkdir -p $STAGE $OUT
run build --release --lib --target wasm32-unknown-unknown
wasm-bindgen --target web --out-dir $STAGE --out-name engine $TARGET/wasm32-unknown-unknown/release/fs_live_kit.wasm
wasm-opt -O3 --enable-simd --enable-bulk-memory --enable-nontrapping-float-to-int --enable-sign-ext --enable-mutable-globals \
  --enable-reference-types --enable-multivalue $STAGE/engine_bg.wasm -o $STAGE/engine_bg.opt.wasm
# RUSTFLAGS in the environment replaces the target flags above and would drop simd128 without a word
if command -v wasm-objdump > /dev/null && ! wasm-objdump -d $STAGE/engine_bg.opt.wasm | grep -q 'f64x2\.'; then
  echo "error: no f64x2 instructions in the kit wasm; was +simd128 dropped (RUSTFLAGS set?)" >&2
  exit 1
fi
cp $STAGE/engine.js $OUT/engine.js
cp $STAGE/engine.d.ts $OUT/engine.d.ts
cp $STAGE/engine_bg.opt.wasm $OUT/engine_bg.wasm
# each piece's stroke loudness per grid, for the Worker's grid-independent level (levels.js); the
# table names this wasm's SHA-256 and is ignored beside any other build
node $HERE/../tools/kit_levels.mjs $OUT > $STAGE/engine_levels.json
cp $STAGE/engine_levels.json $OUT/engine_levels.json
echo "engine_bg.wasm raw=$(wc -c < $OUT/engine_bg.wasm) gzip=$(gzip -9 -c $OUT/engine_bg.wasm | wc -c)"
