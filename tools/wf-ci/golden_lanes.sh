#!/bin/bash
# E6.2 golden-lane CI (bead wf-root-guzez.7.2). Runs every lane this
# host can execute against its per-lane golden and emits JSONL.
# Six-lane matrix status:
#   aarch64 x {debug, release}     — runs here (native goldens)
#   wasm-in-node x {debug, release}— runs here (wasm golden)
#   x86 x {debug, release}         — needs an x86 CI host (recorded)
# Cross-lane identity (native vs wasm) is EXPECTED-DIVERGENT until the
# FMA-contraction remediation lands (bead guzez.7.2.1) — the pair is
# compared and reported LOUDLY either way.
# Repro: tools/wf-ci/golden_lanes.sh
set -u
cd "$(dirname "$0")/../.."
NATIVE_GOLDEN="76119b72d2a3296a5b1f41ff9f21291ac4b54172e47d24116c5543017f07ee39"
WASM_GOLDEN="57e7f8769ab0e6c6252aa716681b34e65760649ac4f8080343287d4f3eb2a5ef"
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
  echo "{\"suite\":\"wf-golden-ci\",\"check\":\"cross-lane-identity\",\"pass\":true,\"note\":\"REMEDIATED? flip guzez.7.2.1 and unify the goldens\"}"
else
  echo "{\"suite\":\"wf-golden-ci\",\"check\":\"cross-lane-identity\",\"pass\":false,\"expected_divergent\":true,\"tracked\":\"guzez.7.2.1\",\"native\":\"$REL\",\"wasm\":\"$WREL\"}"
fi
exit $fail
