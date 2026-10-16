#!/usr/bin/env bash
#
# cooling_01.sh — Journey A acceptance lane (beads frankensim-rc-root-q61wp.13
# and frankensim-extreal-program-f85xj.6.11): the tracked reference cooling
# project goes through the REAL `frankensim` binary — validate, import, the
# seven-stage solve (import-verify, assign, material-resolve, flow-network,
# conduction, qoi, report), the report/package export verbs, the
# determinism repeat and the resume drills — and the run retains the spine
# receipt that `xtask check-maturity` reads for the L3 claim.
#
# This file used to assert the OLD truth (solve exits 5 at the conduction
# gap; report/package refuse). That truth is gone: every stage executes.
# Rather than keep two lanes that can drift, this lane delegates to the
# no-mock producers harness and adds only what Journey A needs on top:
# the retained receipt path and the exit-code contract callers rely on.
#
# Usage:
#   scripts/e2e/cooling_01.sh [--list|--check|--self-test|--run|--retain]
#
#   --run     full profile, receipt written to the artifact dir only
#   --retain  full profile AND overwrite the tracked spine-e2e-summary.json
#             (commit it in the same change as anything it must attest)
#
# Environment: FRANKENSIM_BIN (optional prebuilt binary), ARTIFACT_DIR.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMAND="${1:---run}"
ARTIFACT_DIR="${ARTIFACT_DIR:-${REPO_ROOT}/target/cooling-01}"
LANE="${REPO_ROOT}/scripts/ci/solve_stage_producers_e2e.sh"
TRACKED_RECEIPT="${REPO_ROOT}/spine-e2e-summary.json"

BINARY="${FRANKENSIM_BIN:-}"
if [ -z "${BINARY}" ]; then
  for candidate in \
    "${CARGO_TARGET_DIR:-${REPO_ROOT}/target}/release/frankensim" \
    "${CARGO_TARGET_DIR:-${REPO_ROOT}/target}/debug/frankensim" \
    "${REPO_ROOT}/target/release/frankensim" \
    "${REPO_ROOT}/target/debug/frankensim"; do
    if [ -x "${candidate}" ]; then BINARY="${candidate}"; break; fi
  done
  if [ -z "${BINARY}" ] && command -v frankensim >/dev/null 2>&1; then
    BINARY="$(command -v frankensim)"
  fi
fi

case "${COMMAND}" in
  --list)
    printf "cooling_01::seven_stage_solve_on_reference_project\n"
    printf "cooling_01::report_and_package_exports_reproduce_retained_bytes\n"
    printf "cooling_01::determinism_repeat_and_resume_drills\n"
    printf "cooling_01::refusal_drills\n"
    printf "cooling_01::spine_receipt_retained\n"
    exit 0
    ;;
  --check|--self-test)
    if [ -n "${BINARY}" ] && [ -x "${BINARY}" ] && [ -x "${LANE}" ]; then
      printf 'ok: binary %s, lane %s\n' "${BINARY}" "${LANE}"
      exit 0
    fi
    printf 'failed: need an executable frankensim binary (FRANKENSIM_BIN) and %s\n' "${LANE}" >&2
    exit 1
    ;;
  --run|--retain)
    mkdir -p "${ARTIFACT_DIR}"
    RECEIPT="${ARTIFACT_DIR}/spine-e2e-summary.json"
    if [ "${COMMAND}" = "--retain" ]; then RECEIPT="${TRACKED_RECEIPT}"; fi
    ARGS=(--profile full --through report --artifact-dir "${ARTIFACT_DIR}" --retain-receipt "${RECEIPT}")
    if [ -n "${BINARY}" ]; then ARGS+=(--binary "${BINARY}"); fi
    "${LANE}" "${ARGS[@]}"
    test -s "${RECEIPT}"
    grep -q '"stage": "conduction", "capability": "thermal.conduction-solve", "status": "executed"' "${RECEIPT}"
    printf 'cooling_01: seven stages executed; receipt at %s\n' "${RECEIPT}"
    exit 0
    ;;
  *)
    printf "FATAL: unknown command %s\n" "${COMMAND}" >&2
    exit 2
    ;;
esac
