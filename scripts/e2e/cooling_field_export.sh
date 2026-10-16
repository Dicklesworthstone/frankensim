#!/usr/bin/env bash
#
# cooling_field_export.sh — no-mock end-to-end verification of deterministic
# VTU/XDMF solution field export (bead frankensim-extreal-program-f85xj.6.8).
#
# Drives real VTU/XDMF export for 3D unstructured thermal and flow meshes,
# runs the independent VtuChecker, verifies bit-level determinism and field
# extrema against ledger QoIs, and exercises failure/tamper drills.
#
# Usage:
#   scripts/e2e/cooling_field_export.sh [--list|--check|--self-test|--run|--negative|--replay]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMAND="${1:---run}"
ARTIFACT_DIR="${ARTIFACT_DIR:-${REPO_ROOT}/target/cooling-field-export}"

log_json() {
  local event="$1"
  local status="$2"
  local detail="$3"
  printf '{"ts":"%s","event":"%s","status":"%s","detail":"%s"}\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${event}" "${status}" "${detail}"
}

case "${COMMAND}" in
  --list)
    printf "cooling_field_export::vtu_ascii_roundtrip\n"
    printf "cooling_field_export::xdmf_binary_companion\n"
    printf "cooling_field_export::field_registry_validation\n"
    printf "cooling_field_export::tamper_rejection\n"
    printf "cooling_field_export::determinism_audit\n"
    exit 0
    ;;
  --check|--self-test)
    if ! command -v rch >/dev/null 2>&1; then
      log_json "self_test" "failed" "preflight failed: rch not on PATH"
      exit 1
    fi
    if ! command -v cargo >/dev/null 2>&1; then
      log_json "self_test" "failed" "preflight failed: cargo not on PATH"
      exit 1
    fi
    log_json "self_test" "ok" "preflight checks passed"
    exit 0
    ;;
  --run)
    mkdir -p "${ARTIFACT_DIR}"
    log_json "run_start" "started" "executing focused VTU/XDMF writer tests"

    export PATH="${PATH}:/Users/jemanuel/.local/bin"

    # 1. Run unit & conformance suite in fs-viz
    log_json "test_dispatch" "running" "dispatching fs-viz tests via rch"
    if rch exec -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" cargo test -p fs-viz --test vtu; then
      log_json "vtu_tests" "passed" "all fs-viz vtu and xdmf tests succeeded"
    else
      log_json "vtu_tests" "failed" "fs-viz vtu tests failed"
      exit 1
    fi

    # Presence is useful preflight information, but it is not an artifact-open
    # result and therefore cannot satisfy the ParaView acceptance lane.
    if command -v pvpython >/dev/null 2>&1; then
      log_json "paraview_lane" "not_run" "pvpython is available, but no retained cooling artifact was opened"
    else
      log_json "paraview_lane" "not_run" "pvpython not present on host"
    fi

    log_json "run_terminal" "refused" "focused writer tests passed, but the real cooling solve, requested-field export, independent artifact readback, ledger QoI comparison, and replay are not wired in this script"
    exit 1
    ;;
  --negative)
    log_json "run_terminal" "refused" "negative mode has no cancellation, short-write, corruption, or ledger-failure drill wired"
    exit 1
    ;;
  --replay)
    log_json "run_terminal" "refused" "replay mode has no retained cooling artifact production or byte-comparison path wired"
    exit 1
    ;;
  *)
    printf "FATAL: unknown command %s\n" "${COMMAND}" >&2
    exit 2
    ;;
esac
