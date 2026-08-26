#!/usr/bin/env bash
#
# cooling_convergence.sh — no-mock end-to-end verification of automated
# mesh convergence ladder evaluation with Richardson extrapolation
# (bead frankensim-extreal-program-f85xj.6.6).
#
# Usage:
#   scripts/e2e/cooling_convergence.sh [--list|--check|--self-test|--run|--negative|--replay]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMAND="${1:---run}"
ARTIFACT_DIR="${ARTIFACT_DIR:-${REPO_ROOT}/target/cooling-convergence}"

log_json() {
  local event="$1"
  local status="$2"
  local detail="$3"
  printf '{"ts":"%s","event":"%s","status":"%s","detail":"%s"}\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${event}" "${status}" "${detail}"
}

case "${COMMAND}" in
  --list)
    printf "cooling_convergence::order_2_asymptotic\n"
    printf "cooling_convergence::insufficient_rungs_refusal\n"
    printf "cooling_convergence::oscillatory_refusal\n"
    printf "cooling_convergence::richardson_extrapolation\n"
    printf "cooling_convergence::gci_evaluation\n"
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
  --run|--negative|--replay)
    mkdir -p "${ARTIFACT_DIR}"
    log_json "run_start" "started" "executing mesh convergence ladder e2e suite"

    export PATH="${PATH}:/Users/jemanuel/.local/bin"

    # 1. Run unit & conformance suite in fs-ladder
    log_json "test_dispatch" "running" "dispatching fs-ladder convergence tests via rch"
    if rch exec -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" cargo test -p fs-ladder --test convergence; then
      log_json "convergence_tests" "passed" "all fs-ladder convergence and extrapolation tests succeeded"
    else
      log_json "convergence_tests" "failed" "fs-ladder convergence tests failed"
      exit 1
    fi

    log_json "run_terminal" "pass" "cooling convergence ladder evaluated with verified order and GCI"
    exit 0
    ;;
  *)
    printf "FATAL: unknown command %s\n" "${COMMAND}" >&2
    exit 2
    ;;
esac
