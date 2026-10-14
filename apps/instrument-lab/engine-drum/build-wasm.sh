#!/bin/zsh
# Build the live drum engine for the browser and install it under
# apps/instrument-lab/public/live/drum/ as engine.js + engine_bg.wasm.
#
# The real checkout cannot resolve its sibling repos (they drifted off their
# pins), so the crate is mirrored into the build-only export and built there.
# Nothing under the export's crates/ tree is touched.
#
#   ./build-wasm.sh            wasm release build + wasm-bindgen + wasm-opt
#   ./build-wasm.sh native     native release build of the bench example + tests
set -eu
HERE=${0:A:h}
B=${FSIM_BUILD_ROOT:-/var/folders/vt/n2xyn_s51b97_j3yh2qbqcnc0000gn/T/fsim-build/frankensim}
TARGET=${FSIM_LIVE_DRUM_TARGET:-/var/folders/vt/n2xyn_s51b97_j3yh2qbqcnc0000gn/T/rch_target_frankensim_live_drum}
OUT=$HERE/../public/live/drum
MIRROR=$B/apps/instrument-lab/engine-drum

mkdir -p $MIRROR
# The committed Cargo.lock (if any) wins; otherwise keep the mirror's own.
if [[ -f $HERE/Cargo.lock ]]; then
  rsync -a --delete --exclude target $HERE/ $MIRROR/
else
  rsync -a --delete --exclude Cargo.lock --exclude target $HERE/ $MIRROR/
fi
cd $MIRROR
run() { env RCH_SHIM_LOCAL_IDE=1 CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-8} CARGO_TARGET_DIR=$TARGET cargo "$@"; }

if [[ "${1:-wasm}" == "native" ]]; then
  run build --release --example bench
  run test --release
  echo "bench: $TARGET/release/examples/bench"
  exit 0
fi

STAGE=$TARGET/bindgen-stage
mkdir -p $STAGE $OUT
RAW=$TARGET/wasm32-unknown-unknown/release/fs_live_drum.wasm
WASM_OPT_FLAGS=(-O3 --enable-bulk-memory --enable-nontrapping-float-to-int --enable-sign-ext
  --enable-mutable-globals --enable-reference-types --enable-multivalue)

# 1. Baseline: runs everywhere; scalar FMA is compiler-builtins software.
run build --release --lib --target wasm32-unknown-unknown
wasm-bindgen --target web --out-dir $STAGE --out-name engine $RAW
wasm-opt $WASM_OPT_FLAGS $STAGE/engine_bg.wasm -o $STAGE/engine_bg.opt.wasm
cp $STAGE/engine.js $OUT/engine.js
cp $STAGE/engine.d.ts $OUT/engine.d.ts
cp $STAGE/engine_bg.opt.wasm $OUT/engine_bg.wasm

# 2. Hardware FMA through relaxed-simd (see fma-shim/SAFETY.md): bit-identical
#    results, several times faster where the CPU fuses. The worker falls back
#    to the baseline when the browser cannot compile relaxed-simd.
run build --release --lib --target wasm32-unknown-unknown --features hardware-fma
wasm-bindgen --target web --out-dir $STAGE --out-name engine_fma $RAW
wasm-opt $WASM_OPT_FLAGS --enable-simd --enable-relaxed-simd $STAGE/engine_fma_bg.wasm -o $STAGE/engine_fma_bg.opt.wasm
cp $STAGE/engine_fma.js $OUT/engine_fma.js
cp $STAGE/engine_fma_bg.opt.wasm $OUT/engine_fma_bg.wasm

# Keep the resolved lockfile beside the source for reproducible rebuilds.
cp $MIRROR/Cargo.lock $HERE/Cargo.lock
for f in engine_bg.wasm engine_fma_bg.wasm; do
  echo "$f raw=$(wc -c < $OUT/$f) gzip=$(gzip -9 -c $OUT/$f | wc -c)"
done
