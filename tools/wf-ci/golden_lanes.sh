#!/bin/bash
# E6.2 golden-lane CI (bead wf-root-guzez.7.2). Runs every lane this
# host can execute against its per-lane golden and emits JSONL.
# Six-lane matrix status:
#   aarch64 x {debug, release}     — runs here (native goldens)
#   wasm-in-node x {debug, release}— runs here (wasm golden)
#   x86 x {debug, release}         — needs an x86 CI host (recorded)
# Cross-lane identity (native vs wasm) is REQUIRED-IDENTICAL (bead
# guzez.7.2.1 closed 2026-08-21: the divergence was platform libm —
# acos in the fs-airscrew Prandtl factors, tanh in fs-flyer contact
# friction — routed through det::, NOT FP contraction). ONE golden
# now pins all lanes and a cross-lane mismatch FAILS the run.
# Repro: tools/wf-ci/golden_lanes.sh
set -u
cd "$(dirname "$0")/../.."
GOLDEN="2c50f8a672617cd3e872dfbfb706d4ff26d5828b25d710c910d8e8d632c50714"
NATIVE_GOLDEN="$GOLDEN"
WASM_GOLDEN="$GOLDEN"
PKG="${WF_PKG_DIR:-/tmp/wf-ci-pkg}"
fail=0
lane() { # name digest golden
  if [ "$2" = "$3" ]; then
    echo "{\"suite\":\"wf-golden-ci\",\"lane\":\"$1\",\"pass\":true,\"digest\":\"$2\"}"
  else
    echo "{\"suite\":\"wf-golden-ci\",\"lane\":\"$1\",\"pass\":false,\"digest\":\"$2\",\"golden\":\"$3\"}"
    fail=1
  fi
}
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-${RCH_TARGET_BASE:-/tmp}/rch_target_frankensim_visco}"
REL=$(cargo run -q -p fs-flyer --release --bin canonical_digest 2>/dev/null | python3 -c 'import json,sys;print(json.load(sys.stdin)["digest"])')
lane "aarch64-release" "$REL" "$NATIVE_GOLDEN"
DBG=$(cargo run -q -p fs-flyer --bin canonical_digest 2>/dev/null | python3 -c 'import json,sys;print(json.load(sys.stdin)["digest"])')
lane "aarch64-debug" "$DBG" "$NATIVE_GOLDEN"
( cd crates/fs-flyer-wasm && CARGO_TARGET_DIR="${RCH_TARGET_BASE:-/tmp}/rch_target_fsflyerwasm" \
  wasm-pack build --target nodejs --release --out-dir "$PKG-rel" >/dev/null 2>&1 )
WREL=$(node -e 'const w=require(process.argv[1]);const r=JSON.parse(w.flyer_selftest());console.log(r.ok?r.ok.digest:"SELFTEST-REFUSED:"+r.refusal.code);' "$PKG-rel/fs_flyer_wasm.js")
lane "wasm-node-release" "$WREL" "$WASM_GOLDEN"
( cd crates/fs-flyer-wasm && CARGO_TARGET_DIR="${RCH_TARGET_BASE:-/tmp}/rch_target_fsflyerwasm" \
  wasm-pack build --target nodejs --dev --out-dir "$PKG-dbg" >/dev/null 2>&1 )
WDBG=$(node -e 'const w=require(process.argv[1]);const r=JSON.parse(w.flyer_selftest());console.log(r.ok?r.ok.digest:"SELFTEST-REFUSED:"+r.refusal.code);' "$PKG-dbg/fs_flyer_wasm.js")
lane "wasm-node-debug" "$WDBG" "$WASM_GOLDEN"
echo "{\"suite\":\"wf-golden-ci\",\"lane\":\"x86-release\",\"pass\":null,\"status\":\"NO-DATA: needs an x86 CI host\"}"
echo "{\"suite\":\"wf-golden-ci\",\"lane\":\"x86-debug\",\"pass\":null,\"status\":\"NO-DATA: needs an x86 CI host\"}"
if [ "$REL" = "$WREL" ]; then
  echo "{\"suite\":\"wf-golden-ci\",\"check\":\"cross-lane-identity\",\"pass\":true,\"digest\":\"$REL\"}"
else
  echo "{\"suite\":\"wf-golden-ci\",\"check\":\"cross-lane-identity\",\"pass\":false,\"required_identical\":true,\"closed\":\"guzez.7.2.1\",\"native\":\"$REL\",\"wasm\":\"$WREL\"}"
  fail=1
fi
exit $fail
