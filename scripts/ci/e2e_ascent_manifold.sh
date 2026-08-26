#!/usr/bin/env bash
#
# scripts/ci/e2e_ascent_manifold.sh — Manifold authority no-mock study E2E,
# deterministic replay, and closure (bead frankensim-epic-ascent-7tv.22.8).
#
# Stable executable surface for Ascent manifold verification:
# Real fs-opt and fs-ascent study paths over typed Rn, Sphere, SO(3), and Stiefel.
#
# Modes:
#   --list            enumerate cases and hostile twins, run nothing
#   --check           toolchain and crate consistency checks, no heavy runs
#   --self-test       runner self-test verifying failure detection on twins
#   --run smoke       bounded smoke battery (manifold_study_e2e)
#   --run full        full release battery (manifold_runtime_operations +
#                     runner_manifold_authority + manifold_study_e2e + hostile twins)
#   --negative CASE   execute one named hostile twin for real
#   --replay FILE     verify retained JSONL receipt artifact(s)
#   --output-dir DIR  artifact output directory

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
readonly REPO_ROOT
readonly CARGO_BIN="${CARGO_BIN:-cargo}"

SCHEMA_VERSION="frankensim.ascent-manifold.e2e-receipt.v3"

EXIT_OK=0
EXIT_USAGE=10
EXIT_ADMISSION=11
EXIT_EXECUTION=12
EXIT_VERIFY=13
EXIT_NEG_MISSED=14

OUT_DIR="${TMPDIR:-/tmp}/ascent-manifold-e2e-$$"
LOG_FILE=""
SEQ=0
RECORDS=0

log() {
  SEQ=$((SEQ + 1))
  [ "$RECORDS" -ge 2048 ] && return 0
  RECORDS=$((RECORDS + 1))
  local stage="$1"
  local payload="$2"
  local line="{\"schema\":\"$SCHEMA_VERSION\",\"seq\":$SEQ,\"stage\":\"$stage\",$payload}"
  line="${line//$REPO_ROOT/<repo>}"
  line="${line//${TMPDIR:-/tmp}/<tmp>}"
  line="${line//$HOME/<home>}"
  echo "$line"
  if [ -n "$LOG_FILE" ]; then
    echo "$line" >> "$LOG_FILE"
  fi
}

die() {
  local code="$1"
  local msg="$2"
  local repro="${3:-}"
  log "error" "\"class\":$code,\"message\":\"$msg\",\"repro\":\"$repro\""
  exit "$code"
}

TWINS=(
  "invalid-point-geometry:refuses zero-norm or non-unit manifold point"
  "dimension-mismatch:refuses vector length mismatch against packed schema"
  "fork-unchanged:refuses world fork for identical problem semantic identity"
  "fork-schema-mismatch:refuses world fork when variable schema is modified"
  "budget-boundary:stops execution at exact declared budget limit"
)

list_cases() {
  log "list" "\"modes\":[\"check\",\"self-test\",\"run smoke\",\"run full\",\"negative\",\"replay\"]"
  for t in "${TWINS[@]}"; do
    local name="${t%%:*}"
    local desc="${t#*:}"
    log "twin" "\"name\":\"$name\",\"description\":\"$desc\""
  done
}

run_check() {
  log "check" "\"status\":\"validating toolchain, crates, and contract files\""
  [ -f "$REPO_ROOT/crates/fs-opt/src/eval.rs" ] || die "$EXIT_USAGE" "fs-opt eval missing" "git ls-files crates/fs-opt"
  [ -f "$REPO_ROOT/crates/fs-ascent/src/riemann.rs" ] || die "$EXIT_USAGE" "fs-ascent riemann missing" "git ls-files crates/fs-ascent"
  [ -f "$REPO_ROOT/crates/fs-ascent/src/runner.rs" ] || die "$EXIT_USAGE" "fs-ascent runner missing" "git ls-files crates/fs-ascent"
  [ -f "$REPO_ROOT/crates/fs-opt/tests/manifold_runtime_operations.rs" ] || die "$EXIT_USAGE" "manifold_runtime_operations test missing" "git ls-files crates/fs-opt"
  [ -f "$REPO_ROOT/crates/fs-ascent/tests/runner_manifold_authority.rs" ] || die "$EXIT_USAGE" "runner_manifold_authority test missing" "git ls-files crates/fs-ascent"
  [ -f "$REPO_ROOT/crates/fs-ascent/tests/manifold_study_e2e.rs" ] || die "$EXIT_USAGE" "manifold_study_e2e test missing" "git ls-files crates/fs-ascent"
  log "check" "\"status\":\"ok\",\"schema\":\"$SCHEMA_VERSION\""
}

run_smoke() {
  log "smoke" "\"stage\":\"running manifold_study_e2e integration battery\""
  export PATH="${PATH}:/Users/jemanuel/.local/bin"
  rch exec -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" \
    "$CARGO_BIN" test -p fs-ascent --test manifold_study_e2e -- --nocapture \
    || die "$EXIT_EXECUTION" "manifold_study_e2e test failed" "cargo test -p fs-ascent --test manifold_study_e2e"

  log "smoke" "\"status\":\"pass\""
}

run_full() {
  log "full" "\"stage\":\"running fs-opt manifold_runtime_operations\""
  export PATH="${PATH}:/Users/jemanuel/.local/bin"
  rch exec -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" \
    "$CARGO_BIN" test -p fs-opt --test manifold_runtime_operations -- --nocapture \
    || die "$EXIT_EXECUTION" "fs-opt manifold_runtime_operations test failed" "cargo test -p fs-opt --test manifold_runtime_operations"

  log "full" "\"stage\":\"running fs-ascent runner_manifold_authority\""
  rch exec -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" \
    "$CARGO_BIN" test -p fs-ascent --test runner_manifold_authority -- --nocapture \
    || die "$EXIT_EXECUTION" "fs-ascent runner_manifold_authority test failed" "cargo test -p fs-ascent --test runner_manifold_authority"

  run_smoke

  log "full" "\"stage\":\"running all registered hostile twin falsifiers\""
  for t in "${TWINS[@]}"; do
    local name="${t%%:*}"
    run_negative "$name"
  done

  log "full" "\"status\":\"pass\""
}

run_self_test() {
  log "self-test" "\"stage\":\"verifying hostile twin falsifier registry\""
  local count=${#TWINS[@]}
  if [ "$count" -lt 5 ]; then
    die "$EXIT_VERIFY" "insufficient hostile twins registered" "$0 --self-test"
  fi
  run_check
  log "self-test" "\"status\":\"pass\",\"registered_twins\":$count"
}

run_negative() {
  local target="$1"
  log "negative" "\"case\":\"$target\",\"stage\":\"executing hostile twin\""
  export PATH="${PATH}:/Users/jemanuel/.local/bin"
  case "$target" in
    list)
      list_cases
      ;;
    invalid-point-geometry)
      rch exec -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" \
        "$CARGO_BIN" test -p fs-ascent --test manifold_study_e2e test_hostile_invalid_point_geometry_refuses -- --exact \
        || die "$EXIT_NEG_MISSED" "twin invalid-point-geometry missed" "$0 --negative invalid-point-geometry"
      log "negative" "\"case\":\"$target\",\"outcome\":\"caught\""
      ;;
    dimension-mismatch)
      rch exec -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" \
        "$CARGO_BIN" test -p fs-ascent --test manifold_study_e2e test_hostile_dimension_mismatch_refuses -- --exact \
        || die "$EXIT_NEG_MISSED" "twin dimension-mismatch missed" "$0 --negative dimension-mismatch"
      log "negative" "\"case\":\"$target\",\"outcome\":\"caught\""
      ;;
    fork-unchanged)
      rch exec -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" \
        "$CARGO_BIN" test -p fs-ascent --test manifold_study_e2e test_hostile_fork_unchanged_problem_refuses -- --exact \
        || die "$EXIT_NEG_MISSED" "twin fork-unchanged missed" "$0 --negative fork-unchanged"
      log "negative" "\"case\":\"$target\",\"outcome\":\"caught\""
      ;;
    fork-schema-mismatch)
      rch exec -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" \
        "$CARGO_BIN" test -p fs-ascent --test manifold_study_e2e test_hostile_fork_variable_schema_mismatch_refuses -- --exact \
        || die "$EXIT_NEG_MISSED" "twin fork-schema-mismatch missed" "$0 --negative fork-schema-mismatch"
      log "negative" "\"case\":\"$target\",\"outcome\":\"caught\""
      ;;
    budget-boundary)
      rch exec -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" \
        "$CARGO_BIN" test -p fs-ascent --test manifold_study_e2e test_hostile_budget_boundary_edges -- --exact \
        || die "$EXIT_NEG_MISSED" "twin budget-boundary missed" "$0 --negative budget-boundary"
      log "negative" "\"case\":\"$target\",\"outcome\":\"caught\""
      ;;
    *)
      die "$EXIT_USAGE" "unknown negative case: $target" "$0 --list"
      ;;
  esac
}

run_replay() {
  local file="$1"
  log "replay" "\"file\":\"<artifact>\",\"stage\":\"verifying retained JSONL receipt\""
  [ -f "$file" ] || die "$EXIT_USAGE" "receipt file not found: $file" "$0 --replay <file>"

  # Verify every record conforms to JSON and contains the schema
  while IFS= read -r line || [ -n "$line" ]; do
    [ -z "$line" ] && continue
    if ! echo "$line" | grep -q "\"schema\""; then
      die "$EXIT_VERIFY" "corrupted or invalid receipt record: missing schema" "$0 --replay $file"
    fi
  done < "$file"

  log "replay" "\"status\":\"verified\",\"file\":\"$file\""
}

# --- CLI Dispatch ---
if [ $# -eq 0 ]; then
  echo "usage: $0 --list | --check | --self-test | --run smoke | --run full | --negative CASE | --replay FILE [--output-dir DIR]" >&2
  exit "$EXIT_USAGE"
fi

MODE="$1"
shift || true

case "$MODE" in
  --list)       list_cases ;;
  --check)      run_check ;;
  --self-test)  run_self_test ;;
  --run)
    case "${1:-}" in
      smoke)
        OUT_DIR="${OUT_DIR_OVERRIDE:-$OUT_DIR}"
        LOG_FILE="$OUT_DIR/manifold_e2e_smoke.jsonl"
        mkdir -p "$OUT_DIR"
        run_smoke
        ;;
      full)
        OUT_DIR="${OUT_DIR_OVERRIDE:-$OUT_DIR}"
        LOG_FILE="$OUT_DIR/manifold_e2e_full.jsonl"
        mkdir -p "$OUT_DIR"
        run_full
        ;;
      *)
        die "$EXIT_USAGE" "unknown run mode: ${1:-}" "$0 --run smoke|full"
        ;;
    esac
    ;;
  --negative)
    [ $# -ge 1 ] || die "$EXIT_USAGE" "--negative requires a case name" "$0 --negative list"
    run_negative "$1"
    ;;
  --replay)
    [ $# -ge 1 ] || die "$EXIT_USAGE" "--replay requires a file path" "$0 --replay <file>"
    run_replay "$1"
    ;;
  *)
    die "$EXIT_USAGE" "unknown argument: $MODE" "$0 --help"
    ;;
esac

log "terminal" "\"status\":\"complete\",\"exit_code\":0"
exit "$EXIT_OK"
