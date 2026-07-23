#!/usr/bin/env bash
# quality_lanes.sh (bead huq.18) — the DSR quality lanes the root
# default-feature gate does not cover:
#
#   1. FEATURE MATRIX: every test/bench/bin/example with required-features,
#      DERIVED from `cargo metadata` (never a hand-kept list): adding or
#      removing a gated target changes the lane set automatically.
#   2. HIGH-PRECISION ORACLE: the isolated tools/oracle workspace checks
#      fs-math ULP budgets and fs-ivl enclosures against MPFR at >=200 bits.
#   3. FS-WASM STANDALONE: the nested fs-wasm workspace (native unit
#      tests; the browser build itself stays a wasm-pack lane).
#
# Every lane writes its COMPLETE log and a provisional JSONL verdict to a
# unique run directory. The terminal seal binds all preceding rows to root +
# clean pinned-constellation snapshots before and after execution, so a stable
# dirty root is attributable while a moving root or sibling is a named proof
# failure. DSR requires fs-wasm native
# tests, a locked wasm32 Cargo check, and the wasm-pack browser build: missing
# manifest/tooling is a failure, never a skip.
set -euo pipefail

cd "$(dirname "$0")/../.."

root_dirty() {
  local tracked untracked
  tracked=$(git -c core.fileMode=true -c core.excludesFile=/dev/null \
    status --porcelain --untracked-files=all) || return 1
  untracked=$(git -c core.excludesFile=/dev/null \
    ls-files --others --exclude-per-directory=.gitignore) || return 1
  if [[ -n "$tracked" || -n "$untracked" ]]; then
    printf '%s\n' true
  else
    printf '%s\n' false
  fi
}

absolute_path() {
  case "$1" in
    /*) printf '%s\n' "$1" ;;
    *) printf '%s\n' "$PWD/$1" ;;
  esac
}

START_HEAD=$(git rev-parse HEAD)
HEAD_SHA="$START_HEAD"
DIRTY=$(root_dirty)
RUSTC_VERSION=$(rustc --version)
CARGO_VERSION=$(cargo --version)
LOG_ROOT=$(absolute_path "${FSIM_QUALITY_LOG_DIR:-target/quality-lanes}")
mkdir -p "$LOG_ROOT"
LOG_DIR=$(mktemp -d "$LOG_ROOT/${START_HEAD:0:12}-XXXXXXXX")
RUN_KEY=$(basename "$LOG_DIR")
QUALITY_TARGET_ROOT=$(absolute_path "${FSIM_QUALITY_TARGET_DIR:-target/quality-lanes-builds}")
QUALITY_TARGET_DIR="$QUALITY_TARGET_ROOT/$RUN_KEY"
mkdir -p "$QUALITY_TARGET_DIR"
VERDICTS="$LOG_DIR/verdicts.jsonl"
SNAPSHOT_LOG="$LOG_DIR/snapshot.log"
SNAPSHOT_BEFORE="<unavailable>"
SNAPSHOT_AFTER="<unavailable>"
CONSTELLATION_LOCK_HASH="<unavailable>"
PROVENANCE_STATE="unsealed"
PROOF_STATE="incomplete"
OVERALL_STATUS="fail"
SEAL_WRITTEN=false
: >"$VERDICTS"

FAILURES=0
row() { # lane status detail log
  python3 - "$1" "$2" "$HEAD_SHA" "$DIRTY" "$SNAPSHOT_BEFORE" \
    "$SNAPSHOT_AFTER" "$RUSTC_VERSION" "$CARGO_VERSION" \
    "$CONSTELLATION_LOCK_HASH" "$PROVENANCE_STATE" "$3" "$4" <<'PY' |
import json
import sys

lane, status, head, dirty, before, after, rustc, cargo, lock_hash, provenance, detail, log = sys.argv[1:]
optional = lambda value: None if value == "<unavailable>" else value
print(json.dumps({
    "check": "quality-lane",
    "lane": lane,
    "status": status,
    "head": head,
    "dirty": dirty == "true",
    "provenance_state": provenance,
    "snapshot_before": optional(before),
    "snapshot_after": None,
    "rustc": rustc,
    "cargo": cargo,
    "constellation_lock_sha256": optional(lock_hash),
    "detail": detail,
    "log": log,
}, separators=(",", ":")))
PY
    tee -a "$VERDICTS"
}

emit_proof_seal() { # provenance status
  local provenance="$1" status="$2" prefix_hash seal_row
  if ! prefix_hash=$(shasum -a 256 "$VERDICTS" | awk '{print $1}'); then
    prefix_hash="<unavailable>"
    provenance="incomplete"
    status="fail"
  fi
  if ! seal_row=$(python3 - "$HEAD_SHA" "$DIRTY" "$SNAPSHOT_BEFORE" \
      "$SNAPSHOT_AFTER" "$CONSTELLATION_LOCK_HASH" "$prefix_hash" \
      "$provenance" "$status" "$LOG_DIR" <<'PY'
import json
import sys

head, dirty, before, after, lock_hash, prefix_hash, provenance, status, log_dir = sys.argv[1:]
optional = lambda value: None if value == "<unavailable>" else value
print(json.dumps({
    "identity_domain": "org.frankensim.ci.quality-proof-record.v2",
    "identity_version": 2,
    "check": "quality-proof-seal",
    "provenance_state": provenance,
    "status": status,
    "head": head,
    "dirty": dirty == "true",
    "snapshot_before": optional(before),
    "snapshot_after": optional(after),
    "constellation_lock_sha256": optional(lock_hash),
    "verdicts_prefix_sha256": optional(prefix_hash),
    "log_dir": log_dir,
}, separators=(",", ":")))
PY
  ); then
    return 1
  fi
  if ! printf '%s\n' "$seal_row" | tee -a "$VERDICTS"; then
    return 1
  fi
  SEAL_WRITTEN=true
}

seal_on_exit() { # original-status
  local status="$1"
  trap - EXIT
  if [[ "$SEAL_WRITTEN" != true ]]; then
    # A zero status without the normal terminal seal is itself a failed proof.
    if [[ "$status" -eq 0 ]]; then
      status=1
    fi
    set +e
    emit_proof_seal "incomplete" "fail"
  fi
  exit "$status"
}
trap 'seal_on_exit "$?"' EXIT

if [[ ! -f constellation.lock ]]; then
  printf '%s\n' 'required constellation.lock is missing' >"$SNAPSHOT_LOG"
  row "proof-snapshot" "fail" "required constellation.lock is missing" "$SNAPSHOT_LOG"
  exit 1
fi
if INITIAL_SNAPSHOT="$(scripts/ci/checkout_constellation.sh --snapshot 2>"$SNAPSHOT_LOG")"; then
  SNAPSHOT_BEFORE="$INITIAL_SNAPSHOT"
else
  row "proof-snapshot" "fail" \
    "could not establish initial root + clean pinned-constellation snapshot" "$SNAPSHOT_LOG"
  exit 1
fi
HEAD_SHA=$(git rev-parse HEAD)
DIRTY=$(root_dirty)
if ! CONSTELLATION_LOCK_HASH=$(shasum -a 256 constellation.lock | awk '{print $1}'); then
  row "proof-snapshot" "fail" "could not hash admitted constellation.lock" "$SNAPSHOT_LOG"
  exit 1
fi
if CONFIRMED_SNAPSHOT="$(scripts/ci/checkout_constellation.sh --snapshot 2>>"$SNAPSHOT_LOG")"; then
  if [[ "$CONFIRMED_SNAPSHOT" != "$SNAPSHOT_BEFORE" ]]; then
    row "proof-snapshot" "fail" "root or constellation moved during initial admission" "$SNAPSHOT_LOG"
    exit 1
  fi
else
  row "proof-snapshot" "fail" "could not confirm initial proof snapshot" "$SNAPSHOT_LOG"
  exit 1
fi
CONFIRMED_HEAD=$(git rev-parse HEAD)
CONFIRMED_DIRTY=$(root_dirty)
if ! CONFIRMED_LOCK_HASH=$(shasum -a 256 constellation.lock | awk '{print $1}'); then
  row "proof-snapshot" "fail" "could not confirm admitted lock identity" "$SNAPSHOT_LOG"
  exit 1
fi
if [[ "$CONFIRMED_HEAD" != "$HEAD_SHA" \
    || "$CONFIRMED_DIRTY" != "$DIRTY" \
    || "$CONFIRMED_LOCK_HASH" != "$CONSTELLATION_LOCK_HASH" ]]; then
  row "proof-snapshot" "fail" \
    "HEAD, dirty state, or lock identity moved during initial admission" "$SNAPSHOT_LOG"
  exit 1
fi
PROVENANCE_STATE="provisional"

# ---- Lane set derived from cargo metadata (all required-feature targets) ----
LANES_FILE="$LOG_DIR/gated_lanes.tsv"
METADATA_LOG="$LOG_DIR/metadata.log"
METADATA_OK=true
if ! cargo metadata --locked --format-version 1 --no-deps 2>>"$METADATA_LOG" | python3 -c '
import json, sys
meta = json.load(sys.stdin)
selectors = ("test", "bench", "bin", "example")
rows = []
for pkg in meta["packages"]:
    for t in pkg["targets"]:
        features = t.get("required-features") or []
        kind = next((item for item in selectors if item in t["kind"]), None)
        if features and kind is not None:
            rows.append((pkg["name"], kind, t["name"], ",".join(features)))
for r in sorted(rows):
    print("\t".join(r))
' 2>>"$METADATA_LOG" >"$LANES_FILE"; then
  METADATA_OK=false
  FAILURES=$((FAILURES + 1))
  row "feature-matrix" "fail" "locked Cargo metadata derivation failed" "$METADATA_LOG"
fi

LANE_COUNT=$(wc -l < "$LANES_FILE" | tr -d " ")
if [[ "$METADATA_OK" == true && "$LANE_COUNT" -eq 0 ]]; then
  row "feature-matrix" "fail" "derived zero gated targets — metadata derivation broke (there are known gated targets)" "$LANES_FILE"
  FAILURES=$((FAILURES + 1))
fi

while IFS=$'\t' read -r pkg kind target feats; do
  lane="gated:${pkg}:${kind}:${target}"
  log="$LOG_DIR/${pkg}__${kind}__${target}.log"
  if [[ "$kind" == "test" ]]; then
    cargo_cmd=(cargo test --locked -p "$pkg" --features "$feats" --test "$target")
  else
    cargo_cmd=(cargo check --locked -p "$pkg" --features "$feats" "--$kind" "$target")
  fi
  if env CARGO_TARGET_DIR="$QUALITY_TARGET_DIR" "${cargo_cmd[@]}" >"$log" 2>&1; then
    row "$lane" "pass" "kind=$kind features=$feats" "$log"
  else
    row "$lane" "fail" "kind=$kind features=$feats — see full log" "$log"
    FAILURES=$((FAILURES + 1))
  fi
done < "$LANES_FILE"

# ---- External high-precision oracle (isolated development workspace) ----
ORACLE_MANIFEST="tools/oracle/Cargo.toml"
ORACLE_LOCK="tools/oracle/Cargo.lock"
ORACLE_TARGET_DIR="$QUALITY_TARGET_DIR/high-precision-oracle"
ORACLE_SAMPLES=4096
ORACLE_PRECISION_BITS=256
ORACLE_FMT_LOG="$LOG_DIR/high-precision-oracle-fmt.log"
ORACLE_CLIPPY_LOG="$LOG_DIR/high-precision-oracle-clippy.log"
ORACLE_TEST_LOG="$LOG_DIR/high-precision-oracle-tests.log"
ORACLE_AUDIT_LOG="$LOG_DIR/high-precision-oracle-audit.jsonl"

if [[ -f "$ORACLE_MANIFEST" && -f "$ORACLE_LOCK" ]]; then
  if cargo fmt --manifest-path "$ORACLE_MANIFEST" --check >"$ORACLE_FMT_LOG" 2>&1; then
    row "high-precision-oracle-fmt" "pass" "isolated oracle formatting" "$ORACLE_FMT_LOG"
  else
    row "high-precision-oracle-fmt" "fail" \
      "isolated oracle formatting failed" "$ORACLE_FMT_LOG"
    FAILURES=$((FAILURES + 1))
  fi

  if env CARGO_TARGET_DIR="$ORACLE_TARGET_DIR" \
      cargo clippy --locked --manifest-path "$ORACLE_MANIFEST" --all-targets \
        -- -D warnings >"$ORACLE_CLIPPY_LOG" 2>&1; then
    row "high-precision-oracle-clippy" "pass" \
      "isolated oracle strict clippy" "$ORACLE_CLIPPY_LOG"
  else
    row "high-precision-oracle-clippy" "fail" \
      "isolated oracle strict clippy failed" "$ORACLE_CLIPPY_LOG"
    FAILURES=$((FAILURES + 1))
  fi

  if env CARGO_TARGET_DIR="$ORACLE_TARGET_DIR" \
      cargo test --locked --manifest-path "$ORACLE_MANIFEST" \
        >"$ORACLE_TEST_LOG" 2>&1; then
    row "high-precision-oracle-tests" "pass" \
      "isolated oracle comparison-plumbing tests" "$ORACLE_TEST_LOG"
  else
    row "high-precision-oracle-tests" "fail" \
      "isolated oracle comparison-plumbing tests failed" "$ORACLE_TEST_LOG"
    FAILURES=$((FAILURES + 1))
  fi

  if env CARGO_TARGET_DIR="$ORACLE_TARGET_DIR" \
      cargo run --quiet --locked --manifest-path "$ORACLE_MANIFEST" -- \
        --samples "$ORACLE_SAMPLES" --precision-bits "$ORACLE_PRECISION_BITS" \
        >"$ORACLE_AUDIT_LOG" 2>&1; then
    row "high-precision-oracle-audit" "pass" \
      "MPFR precision_bits=$ORACLE_PRECISION_BITS samples_per_family=$ORACLE_SAMPLES" \
      "$ORACLE_AUDIT_LOG"
  else
    row "high-precision-oracle-audit" "fail" \
      "MPFR comparison found a budget or enclosure violation" "$ORACLE_AUDIT_LOG"
    FAILURES=$((FAILURES + 1))
  fi
else
  ORACLE_MISSING="required isolated oracle manifest or reviewed lock is missing"
  for oracle_lane in fmt clippy tests audit; do
    oracle_log="$LOG_DIR/high-precision-oracle-${oracle_lane}.log"
    printf '%s\n' "$ORACLE_MISSING" >"$oracle_log"
    row "high-precision-oracle-${oracle_lane}" "fail" "$ORACLE_MISSING" "$oracle_log"
    FAILURES=$((FAILURES + 1))
  done
fi

# ---- fs-wasm standalone workspace (native tests) ----
WASM_MANIFEST="crates/fs-wasm/Cargo.toml"
WASM_LOCK="crates/fs-wasm/Cargo.lock"
WASM_LOG="$LOG_DIR/fs-wasm-native.log"
WASM_CHECK_LOG="$LOG_DIR/fs-wasm-browser-locked-check.log"
WASM_BUILD_LOG="$LOG_DIR/fs-wasm-build.log"
WASM_OUT_DIR="$LOG_DIR/fs-wasm-pkg"
WASM_OUT_REL=$(python3 - "$WASM_OUT_DIR" "$PWD/crates/fs-wasm" <<'PY'
import os
import sys

print(os.path.relpath(sys.argv[1], sys.argv[2]))
PY
)
if [[ -f "$WASM_MANIFEST" ]]; then
  if env CARGO_TARGET_DIR="$QUALITY_TARGET_DIR" \
      cargo test --locked --manifest-path "$WASM_MANIFEST" >"$WASM_LOG" 2>&1; then
    row "fs-wasm-native" "pass" "standalone workspace native tests" "$WASM_LOG"
  else
    row "fs-wasm-native" "fail" "standalone workspace native tests — see full log" "$WASM_LOG"
    FAILURES=$((FAILURES + 1))
  fi
else
  printf '%s\n' "required manifest $WASM_MANIFEST is missing" >"$WASM_LOG"
  row "fs-wasm-native" "fail" "required standalone workspace manifest is missing" "$WASM_LOG"
  FAILURES=$((FAILURES + 1))
fi

if [[ -f "$WASM_MANIFEST" ]]; then
  WASM_TOOL_FAILURE=""
  if [[ ! -f "$WASM_LOCK" ]]; then
    WASM_TOOL_FAILURE="required nested lock $WASM_LOCK is missing"
  elif ! command -v wasm-pack >/dev/null 2>&1; then
    WASM_TOOL_FAILURE="required command wasm-pack is not installed"
  elif ! rustup target list --installed 2>/dev/null | grep -q '^wasm32-unknown-unknown$'; then
    WASM_TOOL_FAILURE="required Rust target wasm32-unknown-unknown is not installed"
  fi
  if [[ -n "$WASM_TOOL_FAILURE" ]]; then
    printf '%s\n' "$WASM_TOOL_FAILURE" >"$WASM_BUILD_LOG"
    row "fs-wasm-tooling" "fail" "$WASM_TOOL_FAILURE" "$WASM_BUILD_LOG"
    FAILURES=$((FAILURES + 1))
  else
    if env CARGO_TARGET_DIR="$QUALITY_TARGET_DIR" \
        cargo check --locked --manifest-path "$WASM_MANIFEST" \
          --target wasm32-unknown-unknown >"$WASM_CHECK_LOG" 2>&1; then
      row "fs-wasm-browser-locked-check" "pass" \
        "locked wasm32 Cargo graph checks clean" "$WASM_CHECK_LOG"
    else
      row "fs-wasm-browser-locked-check" "fail" \
        "locked wasm32 Cargo check failed" "$WASM_CHECK_LOG"
      FAILURES=$((FAILURES + 1))
    fi

    if ! WASM_LOCK_BEFORE=$(shasum -a 256 "$WASM_LOCK" | awk '{print $1}'); then
      printf '%s\n' "could not hash required nested lock $WASM_LOCK" >"$WASM_BUILD_LOG"
      row "fs-wasm-build" "fail" "could not hash reviewed nested Cargo.lock" "$WASM_BUILD_LOG"
      FAILURES=$((FAILURES + 1))
    else
      WASM_BUILD_OK=true
      if ! (cd crates/fs-wasm && env CARGO_TARGET_DIR="$QUALITY_TARGET_DIR" \
          wasm-pack build --dev --target web --out-dir "$WASM_OUT_REL" \
            -- --locked >"$WASM_BUILD_LOG" 2>&1); then
        WASM_BUILD_OK=false
      fi
      if ! WASM_LOCK_AFTER=$(shasum -a 256 "$WASM_LOCK" | awk '{print $1}'); then
        row "fs-wasm-build" "fail" \
          "could not re-hash nested Cargo.lock after wasm-pack" "$WASM_BUILD_LOG"
        FAILURES=$((FAILURES + 1))
      elif [[ "$WASM_BUILD_OK" == true && "$WASM_LOCK_AFTER" == "$WASM_LOCK_BEFORE" ]]; then
        row "fs-wasm-build" "pass" "required wasm-pack dev browser build" "$WASM_BUILD_LOG"
      elif [[ "$WASM_LOCK_AFTER" != "$WASM_LOCK_BEFORE" ]]; then
        row "fs-wasm-build" "fail" \
          "wasm-pack modified the reviewed nested Cargo.lock" "$WASM_BUILD_LOG"
        FAILURES=$((FAILURES + 1))
      else
        row "fs-wasm-build" "fail" "wasm-pack dev browser build failed" "$WASM_BUILD_LOG"
        FAILURES=$((FAILURES + 1))
      fi
    fi
  fi
fi

# ---- Content-derived inventory (counts and hashes come from actual bytes) ----
INVENTORY="$LOG_DIR/inventory.json"
if INVENTORY_ROW=$(python3 - "$LOG_DIR" "$HEAD_SHA" "$DIRTY" "$SNAPSHOT_BEFORE" <<'PY'
import hashlib
import json
import pathlib
import sys

log_dir, head, dirty, snapshot = sys.argv[1:]
root = pathlib.Path(".")
patterns = (
    ".cargo/config.toml", "Cargo.toml", "Cargo.lock", "rust-toolchain.toml",
    "constellation.lock", "scripts/ci/*.sh", "xtask/src/**/*.rs",
    "tools/oracle/Cargo.toml", "tools/oracle/Cargo.lock",
    "tools/oracle/README.md", "tools/oracle/src/**/*.rs",
    "crates/*/Cargo.toml", "crates/*/Cargo.lock", "crates/*/.cargo/*.toml",
    "crates/*/CONTRACT.md", "crates/*/build.rs",
    "crates/*/src/**/*.rs", "crates/*/tests/**/*.rs", "crates/*/benches/**/*.rs",
    "crates/*/examples/**/*.rs",
)
paths = sorted({path for pattern in patterns for path in root.glob(pattern) if path.is_file()})
hasher = hashlib.sha256()

def add(label, data):
    encoded = label.encode()
    hasher.update(len(encoded).to_bytes(8, "big"))
    hasher.update(encoded)
    hasher.update(len(data).to_bytes(8, "big"))
    hasher.update(data)

for path in paths:
    add(path.as_posix(), path.read_bytes())
gated_path = pathlib.Path(log_dir, "gated_lanes.tsv")
add("derived:gated_lanes.tsv", gated_path.read_bytes())
crates = sorted(path.parent.name for path in root.glob("crates/*/Cargo.toml"))
contracts = sorted(path.parent.name for path in root.glob("crates/*/CONTRACT.md"))
test_files = sorted(str(path) for path in root.glob("crates/*/tests/*.rs"))
gated = gated_path.read_text().splitlines()
gated_tests = sum(line.split("\t")[1] == "test" for line in gated if "\t" in line)
inventory = {
    "check": "inventory",
    "head": head,
    "dirty": dirty == "true",
    "snapshot": snapshot,
    "crates": len(crates),
    "contracts": len(contracts),
    "test_files": len(test_files),
    "gated_test_targets": gated_tests,
    "gated_required_feature_targets": len(gated),
    "hashed_inputs": len(paths) + 1,
    "source_hash": hasher.hexdigest(),
}
pathlib.Path(log_dir, "inventory.json").write_text(json.dumps(inventory, indent=1) + "\n")
print(json.dumps(inventory, separators=(",", ":")))
PY
); then
  printf '%s\n' "$INVENTORY_ROW" | tee -a "$VERDICTS"
else
  row "inventory" "fail" "could not hash canonical quality-lane inputs" "$INVENTORY"
  FAILURES=$((FAILURES + 1))
fi

if ! LOCK_HASH_AFTER=$(shasum -a 256 constellation.lock | awk '{print $1}'); then
  LOCK_HASH_AFTER="<unavailable>"
  row "proof-snapshot" "fail" "could not hash final constellation.lock" "$SNAPSHOT_LOG"
  FAILURES=$((FAILURES + 1))
elif [[ "$LOCK_HASH_AFTER" != "$CONSTELLATION_LOCK_HASH" ]]; then
  row "proof-snapshot" "fail" "constellation.lock moved during quality proof" "$SNAPSHOT_LOG"
  FAILURES=$((FAILURES + 1))
fi

# This is the single terminal source observation.  Do not perform another
# fallible HEAD/dirty probe after it and then leave that late failure out of the
# exit status: snapshot v2 already binds HEAD, index/worktree bytes, the exact
# lock, and every clean pinned sibling.
if ! SNAPSHOT_AFTER="$(scripts/ci/checkout_constellation.sh --snapshot 2>>"$SNAPSHOT_LOG")"; then
  SNAPSHOT_AFTER="<unavailable>"
  row "proof-snapshot" "fail" \
    "could not establish final root + clean pinned-constellation snapshot" "$SNAPSHOT_LOG"
  FAILURES=$((FAILURES + 1))
elif [[ "$SNAPSHOT_AFTER" != "$SNAPSHOT_BEFORE" ]]; then
  row "proof-snapshot" "fail" \
    "root content, HEAD, dirty state, or constellation moved during quality proof" "$SNAPSHOT_LOG"
  FAILURES=$((FAILURES + 1))
fi

if [[ "$SNAPSHOT_AFTER" == "$SNAPSHOT_BEFORE" \
    && "$LOCK_HASH_AFTER" == "$CONSTELLATION_LOCK_HASH" ]]; then
  PROOF_STATE="sealed"
fi

python3 - "$HEAD_SHA" "$DIRTY" "$SNAPSHOT_BEFORE" "$SNAPSHOT_AFTER" \
  "$CONSTELLATION_LOCK_HASH" "$((LANE_COUNT + 7))" "$FAILURES" "$LOG_DIR" <<'PY' | tee -a "$VERDICTS"
import json
import sys

head, dirty, before, after, lock_hash, lanes, failures, log_dir = sys.argv[1:]
optional = lambda value: None if value == "<unavailable>" else value
print(json.dumps({
    "check": "quality-lanes-summary",
    "provenance_state": "provisional",
    "status": "pass" if int(failures) == 0 else "fail",
    "head": head,
    "dirty": dirty == "true",
    "snapshot_before": optional(before),
    "snapshot_after": optional(after),
    "constellation_lock_sha256": optional(lock_hash),
    "lanes": int(lanes),
    "failures": int(failures),
    "log_dir": log_dir,
}, separators=(",", ":")))
PY

if [[ "$FAILURES" -eq 0 && "$PROOF_STATE" == "sealed" ]]; then
  OVERALL_STATUS="pass"
fi
emit_proof_seal "$PROOF_STATE" "$OVERALL_STATUS"

if [[ "$OVERALL_STATUS" != "pass" ]]; then
  exit 1
fi
