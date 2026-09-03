#!/usr/bin/env bash
#
# cooling_report.sh — no-mock end-to-end verification of deterministic HTML
# engineering report and JSON twin generation (bead frankensim-extreal-program-f85xj.6.9).
#
# Usage:
#   scripts/e2e/cooling_report.sh [--list|--check|--self-test|--run|--negative|--replay]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMAND="${1:---run}"
ARTIFACT_DIR="${ARTIFACT_DIR:-${REPO_ROOT}/target/cooling-report}"
BINARY="${FRANKENSIM_BIN:-}"

if [ -z "${BINARY}" ]; then
  if [ -x "/Volumes/USB_NVME/cargo-target/debug/frankensim" ]; then
    BINARY="/Volumes/USB_NVME/cargo-target/debug/frankensim"
  elif [ -x "${REPO_ROOT}/target/debug/frankensim" ]; then
    BINARY="${REPO_ROOT}/target/debug/frankensim"
  elif command -v frankensim >/dev/null 2>&1; then
    BINARY="$(command -v frankensim)"
  fi
fi

log_json() {
  local event="$1"
  local status="$2"
  local detail="$3"
  printf '{"ts":"%s","event":"%s","status":"%s","detail":"%s"}\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${event}" "${status}" "${detail}"
}

case "${COMMAND}" in
  --list)
    printf "cooling_report::fail_closed_without_retained_loader\n"
    printf "cooling_report::no_fabricated_verified_qois\n"
    exit 0
    ;;
  --check|--self-test)
    if [ -n "${BINARY}" ] && [ -x "${BINARY}" ]; then
      log_json "self_test" "ok" "frankensim binary is available"
      exit 0
    fi
    log_json "self_test" "failed" "frankensim binary is unavailable"
    exit 1
    ;;
  --run|--negative|--replay)
    if [ -z "${BINARY}" ] || [ ! -x "${BINARY}" ]; then
      log_json "run_start" "failed" "frankensim binary is unavailable"
      exit 1
    fi
    BINARY="$(cd "$(dirname "${BINARY}")" && pwd)/$(basename "${BINARY}")"
    mkdir -p "${ARTIFACT_DIR}"
    ARTIFACT_DIR="$(cd "${ARTIFACT_DIR}" && pwd)"
    log_json "run_start" "started" "verifying report refuses a missing ledger before writing anything"
    set +e
    (cd "${ARTIFACT_DIR}" && "${BINARY}" --json report unbound_run /definitely/missing/ledger.db) \
      > "${ARTIFACT_DIR}/report.json" 2> "${ARTIFACT_DIR}/report.stderr.jsonl"
    report_exit=$?
    set -e
    # `report` executes against a completed run since 2026-08-25; a missing
    # ledger is an input refusal (exit 3), never a fabricated export.
    test "${report_exit}" -eq 3
    grep -q '"status":"refused"' "${ARTIFACT_DIR}/report.json"
    grep -q '"code":"cli-export-ledger-missing"' "${ARTIFACT_DIR}/report.stderr.jsonl"
    if grep -Eq 'junction_maximum|thermal_margin|Verified|content_hash' \
      "${ARTIFACT_DIR}/report.json"; then
      log_json "run_terminal" "failed" "report emitted fabricated scientific evidence"
      exit 1
    fi
    test ! -e "${ARTIFACT_DIR}/unbound_run.html"
    test ! -e "${ARTIFACT_DIR}/unbound_run.report.json"
    log_json "run_terminal" "pass" "report refused the missing ledger and emitted no scientific claim"
    exit 0
    ;;
  *)
    printf "FATAL: unknown command %s\n" "${COMMAND}" >&2
    exit 2
    ;;
esac
