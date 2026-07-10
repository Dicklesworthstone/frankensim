#!/usr/bin/env bash
# quality_lanes.sh (bead huq.18) — the DSR quality lanes the root
# default-feature gate does not cover:
#
#   1. FEATURE MATRIX: every [[test]] target with required-features,
#      DERIVED from `cargo metadata` (never a hand-kept list): adding or
#      removing a gated target changes the lane set automatically.
#   2. FS-WASM STANDALONE: the nested fs-wasm workspace (native unit
#      tests; the browser build itself stays a wasm-pack lane).
#
# Every lane writes its COMPLETE log to an artifact path (no truncation)
# and a JSONL verdict row carrying the tested HEAD, tree-dirty flag, and
# log path, so a red gate is diagnosable and cannot be misattributed
# across concurrent commits. A lane that cannot run is SKIPPED WITH A
# NAMED REASON, never silently passed. Exit is non-zero if any lane
# fails. An inventory.json (counts + source hash, derived from the
# manifests) is emitted alongside the logs.
set -euo pipefail

cd "$(dirname "$0")/../.."

HEAD_SHA=$(git rev-parse HEAD)
DIRTY=$(test -n "$(git status --porcelain 2>/dev/null)" && echo true || echo false)
LOG_DIR="${FSIM_QUALITY_LOG_DIR:-$PWD/target/quality-lanes/${HEAD_SHA:0:12}}"
mkdir -p "$LOG_DIR"

FAILURES=0
row() { # lane status detail log
  printf '{"check":"quality-lane","lane":"%s","status":"%s","head":"%s","dirty":%s,"detail":"%s","log":"%s"}\n' \
    "$1" "$2" "$HEAD_SHA" "$DIRTY" "$3" "$4"
}

# ---- Lane set derived from cargo metadata (required-features tests) ----
LANES_FILE="$LOG_DIR/gated_lanes.tsv"
cargo metadata --format-version 1 --no-deps 2>/dev/null | python3 -c '
import json, sys
meta = json.load(sys.stdin)
rows = []
for pkg in meta["packages"]:
    for t in pkg["targets"]:
        if "test" in t["kind"] and t.get("required-features"):
            rows.append((pkg["name"], t["name"], ",".join(t["required-features"])))
for r in sorted(rows):
    print("\t".join(r))
' > "$LANES_FILE"

LANE_COUNT=$(wc -l < "$LANES_FILE" | tr -d " ")
if [ "$LANE_COUNT" -eq 0 ]; then
  row "feature-matrix" "fail" "derived zero gated targets — metadata derivation broke (there are known gated targets)" "$LANES_FILE"
  FAILURES=$((FAILURES + 1))
fi

while IFS=$'\t' read -r pkg target feats; do
  lane="gated:${pkg}:${target}"
  log="$LOG_DIR/${pkg}__${target}.log"
  if cargo test -p "$pkg" --features "$feats" --test "$target" >"$log" 2>&1; then
    row "$lane" "pass" "features=$feats" "$log"
  else
    row "$lane" "fail" "features=$feats — see full log" "$log"
    FAILURES=$((FAILURES + 1))
  fi
done < "$LANES_FILE"

# ---- fs-wasm standalone workspace (native tests) ----
WASM_MANIFEST="crates/fs-wasm/Cargo.toml"
WASM_LOG="$LOG_DIR/fs-wasm-native.log"
if [ -f "$WASM_MANIFEST" ]; then
  if cargo test --manifest-path "$WASM_MANIFEST" >"$WASM_LOG" 2>&1; then
    row "fs-wasm-native" "pass" "standalone workspace native tests" "$WASM_LOG"
  else
    row "fs-wasm-native" "fail" "standalone workspace native tests — see full log" "$WASM_LOG"
    FAILURES=$((FAILURES + 1))
  fi
else
  row "fs-wasm-native" "skip" "crates/fs-wasm/Cargo.toml missing" "-"
fi
# The wasm32 browser BUILD needs wasm-pack + target; skip with the reason
# when absent rather than pretending coverage.
if command -v wasm-pack >/dev/null 2>&1 && rustup target list --installed 2>/dev/null | grep -q wasm32-unknown-unknown; then
  WASM_BUILD_LOG="$LOG_DIR/fs-wasm-build.log"
  if (cd crates/fs-wasm && wasm-pack build --dev --target web >"$WASM_BUILD_LOG" 2>&1); then
    row "fs-wasm-build" "pass" "wasm-pack dev build" "$WASM_BUILD_LOG"
  else
    row "fs-wasm-build" "fail" "wasm-pack dev build — see full log" "$WASM_BUILD_LOG"
    FAILURES=$((FAILURES + 1))
  fi
else
  row "fs-wasm-build" "skip" "wasm-pack or wasm32-unknown-unknown target not installed on this host" "-"
fi

# ---- Derived inventory (counts come from the tree, never by hand) ----
python3 - "$LOG_DIR" "$HEAD_SHA" "$DIRTY" <<'PY'
import hashlib, json, pathlib, sys
log_dir, head, dirty = sys.argv[1], sys.argv[2], sys.argv[3]
crates = sorted(p.parent.name for p in pathlib.Path("crates").glob("*/Cargo.toml"))
contracts = sorted(p.parent.name for p in pathlib.Path("crates").glob("*/CONTRACT.md"))
test_files = sorted(str(p) for p in pathlib.Path("crates").glob("*/tests/*.rs"))
gated = sorted(pathlib.Path(log_dir, "gated_lanes.tsv").read_text().splitlines())
src = "\n".join(crates + contracts + test_files + gated)
inv = {
    "head": head,
    "dirty": dirty == "true",
    "crates": len(crates),
    "contracts": len(contracts),
    "test_files": len(test_files),
    "gated_test_targets": len(gated),
    "source_hash": hashlib.sha256(src.encode()).hexdigest()[:16],
}
out = pathlib.Path(log_dir, "inventory.json")
out.write_text(json.dumps(inv, indent=1) + "\n")
print(json.dumps({"check": "inventory", **inv}))
PY

echo "{\"check\":\"quality-lanes-summary\",\"head\":\"$HEAD_SHA\",\"dirty\":$DIRTY,\"lanes\":$((LANE_COUNT + 2)),\"failures\":$FAILURES,\"log_dir\":\"$LOG_DIR\"}"
exit "$FAILURES"
