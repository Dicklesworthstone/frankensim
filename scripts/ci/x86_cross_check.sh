#!/usr/bin/env bash
# x86_cross_check.sh (bead ebro) - the ARM-INVISIBLE BREAKAGE gate.
#
# The fleet develops on aarch64. This lane cross-checks the workspace's
# default target surface, then derives every required-feature test, bench,
# binary, and example from locked Cargo metadata. `cargo check --target` does
# not link, so runtime-only x86 behavior still belongs to Threadripper lanes.
#
# The retained JSONL verdict binds HEAD, constellation.lock, the exact root
# tracked/untracked content, and every clean pinned sibling tree. A unique run
# directory retains the full log. The composite snapshot is recomputed after
# Cargo exits, so root or constellation movement is a named proof failure.
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
LOG_ROOT=$(absolute_path "${FSIM_X86_CROSS_LOG_DIR:-target/x86-cross-logs}")
mkdir -p "$LOG_ROOT"
LOG_DIR=$(mktemp -d "$LOG_ROOT/${START_HEAD:0:12}-XXXXXXXX")
RUN_KEY=$(basename "$LOG_DIR")
CROSS_TARGET_ROOT=$(absolute_path "${FSIM_X86_CROSS_TARGET_DIR:-target/x86-cross-builds}")
CROSS_TARGET_DIR="$CROSS_TARGET_ROOT/$RUN_KEY"
mkdir -p "$CROSS_TARGET_DIR"
LOG="$LOG_DIR/check.log"
GATED_TARGETS="$LOG_DIR/required_feature_targets.tsv"
VERDICTS="$LOG_DIR/verdicts.jsonl"
SNAPSHOT_LOG="$LOG_DIR/snapshot.log"
SNAPSHOT_BEFORE="<unavailable>"
SNAPSHOT_AFTER="<unavailable>"
CONSTELLATION_LOCK_HASH="<unavailable>"
PROVENANCE_STATE="unsealed"

row() { # status detail [log]
  local row_log="${3:-$LOG}"
  python3 - "$1" "$HEAD_SHA" "$DIRTY" "$SNAPSHOT_BEFORE" "$SNAPSHOT_AFTER" \
    "$RUSTC_VERSION" "$CARGO_VERSION" "$CONSTELLATION_LOCK_HASH" \
    "$PROVENANCE_STATE" "$2" "$row_log" <<'PY' |
import json
import sys

status, head, dirty, snapshot_before, snapshot_after, rustc, cargo, lock_hash, provenance, detail, log = sys.argv[1:]
optional = lambda value: None if value in ("<unavailable>", "<snapshot-failed>") else value
print(json.dumps({
    "identity_domain": "org.frankensim.ci.x86-cross-verdict.v1",
    "identity_version": 1,
    "check": "quality-lane",
    "lane": "x86-cross-all-targets",
    "status": status,
    "head": head,
    "dirty": dirty == "true",
    "provenance_state": provenance,
    "snapshot_before": optional(snapshot_before),
    "snapshot_after": optional(snapshot_after),
    "rustc": rustc,
    "cargo": cargo,
    "constellation_lock_sha256": optional(lock_hash),
    "detail": detail,
    "log": log,
}, separators=(",", ":")))
PY
    tee -a "$VERDICTS"
}

if [[ ! -f constellation.lock ]]; then
  printf '%s\n' 'required constellation.lock is missing' >"$LOG"
  row "fail" "required constellation.lock is missing"
  exit 1
fi
if INITIAL_SNAPSHOT="$(scripts/ci/checkout_constellation.sh --snapshot 2>"$SNAPSHOT_LOG")"; then
  SNAPSHOT_BEFORE="$INITIAL_SNAPSHOT"
else
  printf '%s\n' 'could not establish initial root + clean pinned-constellation snapshot' >"$LOG"
  row "fail" "could not establish initial root + clean pinned-constellation snapshot" "$SNAPSHOT_LOG"
  exit 1
fi
HEAD_SHA=$(git rev-parse HEAD)
DIRTY=$(root_dirty)
if ! CONSTELLATION_LOCK_HASH=$(shasum -a 256 constellation.lock | awk '{print $1}'); then
  printf '%s\n' 'could not hash admitted constellation.lock' >"$LOG"
  row "fail" "could not hash admitted constellation.lock" "$SNAPSHOT_LOG"
  exit 1
fi
if CONFIRMED_SNAPSHOT="$(scripts/ci/checkout_constellation.sh --snapshot 2>>"$SNAPSHOT_LOG")"; then
  if [[ "$CONFIRMED_SNAPSHOT" != "$SNAPSHOT_BEFORE" ]]; then
    printf '%s\n' 'root or constellation moved during initial admission' >"$LOG"
    row "fail" "root or constellation moved during initial admission" "$SNAPSHOT_LOG"
    exit 1
  fi
else
  printf '%s\n' 'could not confirm initial proof snapshot' >"$LOG"
  row "fail" "could not confirm initial proof snapshot" "$SNAPSHOT_LOG"
  exit 1
fi
CONFIRMED_HEAD=$(git rev-parse HEAD)
CONFIRMED_DIRTY=$(root_dirty)
if ! CONFIRMED_LOCK_HASH=$(shasum -a 256 constellation.lock | awk '{print $1}'); then
  printf '%s\n' 'could not confirm admitted lock identity' >"$LOG"
  row "fail" "could not confirm admitted lock identity" "$SNAPSHOT_LOG"
  exit 1
fi
if [[ "$CONFIRMED_HEAD" != "$HEAD_SHA" \
    || "$CONFIRMED_DIRTY" != "$DIRTY" \
    || "$CONFIRMED_LOCK_HASH" != "$CONSTELLATION_LOCK_HASH" ]]; then
  printf '%s\n' 'HEAD, dirty state, or lock identity moved during initial admission' >"$LOG"
  row "fail" "HEAD, dirty state, or lock identity moved during initial admission" "$SNAPSHOT_LOG"
  exit 1
fi
PROVENANCE_STATE="provisional"

FAILURES=0
TARGET_AVAILABLE=true
if ! rustup target list --installed 2>/dev/null | grep -q '^x86_64-unknown-linux-gnu$'; then
  TARGET_AVAILABLE=false
  printf '%s\n' \
    'required rustup target x86_64-unknown-linux-gnu is not installed; run: rustup target add x86_64-unknown-linux-gnu' >"$LOG"
  FAILURES=$((FAILURES + 1))
fi

if [[ "$TARGET_AVAILABLE" == true ]]; then
  {
    printf 'rustc: %s\ncargo: %s\n' "$RUSTC_VERSION" "$CARGO_VERSION"
    printf 'command: cargo check --locked --workspace --all-targets --target x86_64-unknown-linux-gnu\n'
    if ! env CARGO_TARGET_DIR="$CROSS_TARGET_DIR" \
        cargo check --locked --workspace --all-targets --target x86_64-unknown-linux-gnu; then
      FAILURES=$((FAILURES + 1))
    fi
  } >"$LOG" 2>&1

  if ! cargo metadata --locked --format-version 1 --no-deps 2>>"$LOG" | python3 -c '
import json, sys
meta = json.load(sys.stdin)
selectors = ("test", "bench", "bin", "example")
rows = []
for pkg in meta["packages"]:
    for target in pkg["targets"]:
        features = target.get("required-features") or []
        kind = next((item for item in selectors if item in target["kind"]), None)
        if features and kind is not None:
            rows.append((pkg["name"], kind, target["name"], ",".join(features)))
for item in sorted(rows):
    print("\t".join(item))
' 2>>"$LOG" >"$GATED_TARGETS"; then
    printf '%s\n' 'failed to derive required-feature targets from locked Cargo metadata' >>"$LOG"
    FAILURES=$((FAILURES + 1))
  else
    if [ ! -s "$GATED_TARGETS" ]; then
      printf '%s\n' 'derived zero required-feature targets; refusing incomplete x86 coverage' >>"$LOG"
      FAILURES=$((FAILURES + 1))
    fi
    while IFS=$'\t' read -r pkg kind target features; do
      [ -n "$pkg" ] || continue
      printf '\ncommand: cargo check --locked -p %q --features %q --%s %q --target x86_64-unknown-linux-gnu\n' \
        "$pkg" "$features" "$kind" "$target" >>"$LOG"
      if ! env CARGO_TARGET_DIR="$CROSS_TARGET_DIR" \
          cargo check --locked -p "$pkg" --features "$features" "--$kind" "$target" \
            --target x86_64-unknown-linux-gnu >>"$LOG" 2>&1; then
        FAILURES=$((FAILURES + 1))
      fi
    done <"$GATED_TARGETS"
  fi
else
  : >"$GATED_TARGETS"
fi

if ! LOCK_HASH_AFTER=$(shasum -a 256 constellation.lock | awk '{print $1}'); then
  LOCK_HASH_AFTER="<unavailable>"
  printf '%s\n' 'could not hash final constellation.lock' >>"$LOG"
  FAILURES=$((FAILURES + 1))
elif [ "$LOCK_HASH_AFTER" != "$CONSTELLATION_LOCK_HASH" ]; then
  printf '%s\n' 'constellation.lock moved during x86 proof' >>"$LOG"
  FAILURES=$((FAILURES + 1))
fi

# One terminal source observation drives both the seal and process verdict.
# Snapshot v2 already binds HEAD, the complete root source state, the exact
# lock, and every clean pinned sibling; a later auxiliary probe must not create
# an `incomplete` seal while leaving the command green.
if ! SNAPSHOT_AFTER=$(scripts/ci/checkout_constellation.sh --snapshot 2>>"$SNAPSHOT_LOG"); then
  SNAPSHOT_AFTER="<snapshot-failed>"
  printf '%s\n' 'failed to hash the final root + constellation snapshot' >>"$LOG"
  FAILURES=$((FAILURES + 1))
elif [ "$SNAPSHOT_AFTER" != "$SNAPSHOT_BEFORE" ]; then
  printf 'root or constellation moved during proof: before=%s after=%s\n' \
    "$SNAPSHOT_BEFORE" "$SNAPSHOT_AFTER" >>"$LOG"
  FAILURES=$((FAILURES + 1))
fi

PROVENANCE_STATE="incomplete"
if [[ "$SNAPSHOT_AFTER" == "$SNAPSHOT_BEFORE" \
    && "$LOCK_HASH_AFTER" == "$CONSTELLATION_LOCK_HASH" ]]; then
  PROVENANCE_STATE="sealed"
fi

if [ -f "$GATED_TARGETS" ]; then
  GATED_COUNT=$(wc -l <"$GATED_TARGETS" | tr -d ' ')
else
  GATED_COUNT=0
fi
if [[ "$FAILURES" -eq 0 && "$PROVENANCE_STATE" == "sealed" ]]; then
  row "pass" "locked default surface plus ${GATED_COUNT} required-feature targets cross-check clean for x86_64-unknown-linux-gnu"
  exit 0
fi

ERRORS=$(grep -cE '^error(\[|:)' "$LOG" || true)
if [[ "$TARGET_AVAILABLE" == false ]]; then
  row "fail" "required rustup target x86_64-unknown-linux-gnu is not installed; run: rustup target add x86_64-unknown-linux-gnu"
else
  row "fail" "x86 cross-check had ${FAILURES} failing lane or provenance condition(s), with ${ERRORS} compiler error(s)"
fi
exit 1
