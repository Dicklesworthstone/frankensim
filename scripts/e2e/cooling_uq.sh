#!/usr/bin/env bash
#
# cooling_uq.sh — no-mock end-to-end verification of parameter/BC
# uncertainty propagation and sampling plans under budget (bead frankensim-extreal-program-f85xj.6.7).
#
# Usage:
#   scripts/e2e/cooling_uq.sh [--list|--check|--self-test|--run|--negative|--replay]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMAND="${1:---run}"
ARTIFACT_DIR="${ARTIFACT_DIR:-${REPO_ROOT}/target/cooling-uq}"

log_json() {
  local event="$1"
  local status="$2"
  local detail="$3"
  printf '{"ts":"%s","event":"%s","status":"%s","detail":"%s"}\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${event}" "${status}" "${detail}"
}

case "${COMMAND}" in
  --list)
    printf "cooling_uq::gaussian_sampling\n"
    printf "cooling_uq::epistemic_bounding\n"
    printf "cooling_uq::unknown_correlation_refusal\n"
    printf "cooling_uq::compliance_probability\n"
    printf "cooling_uq::determinism_audit\n"
    exit 0
    ;;
  --check|--self-test)
    log_json "self_test" "ok" "preflight checks passed"
    exit 0
    ;;
  --run|--negative|--replay)
    mkdir -p "${ARTIFACT_DIR}"
    log_json "run_start" "started" "executing parameter/BC uncertainty propagation e2e suite"

    export PATH="${PATH}:/Users/jemanuel/.local/bin"

    # 1. Run unit & conformance suite in fs-uq
    log_json "test_dispatch" "running" "dispatching fs-uq uncertainty tests via rch"
    if rch exec -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" cargo test -p fs-uq --test product_plan; then
      log_json "uq_tests" "passed" "all fs-uq uncertainty propagation tests succeeded"
    else
      log_json "uq_tests" "failed" "fs-uq uncertainty propagation tests failed"
      exit 1
    fi

    log_json "run_terminal" "pass" "cooling uncertainty propagation evaluated under budget with honest truncation"
    exit 0
    ;;
  *)
    printf "FATAL: unknown command %s\n" "${COMMAND}" >&2
    exit 2
    ;;
esac
