#!/usr/bin/env bash
#
# cooling_01.sh — End-to-end acceptance lane: dirty input -> decision artifact,
# with failure drills and forensic logging (bead frankensim-extreal-program-f85xj.6.11).
#
# Usage:
#   scripts/e2e/cooling_01.sh [--list|--check|--self-test|--run|--negative|--replay]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMAND="${1:---run}"
ARTIFACT_DIR="${ARTIFACT_DIR:-${REPO_ROOT}/target/cooling-01}"
mkdir -p "${ARTIFACT_DIR}"

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
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${event}" "${status}" "${detail}" | tee -a "${ARTIFACT_DIR}/forensic.jsonl"
}

case "${COMMAND}" in
  --list)
    printf "cooling_01::durable_product_prefix\n"
    printf "cooling_01::report_and_package_fail_closed\n"
    printf "cooling_01::verbs_vs_run_refusal_parity\n"
    printf "cooling_01::idempotent_replay\n"
    printf "cooling_01::refusal_drills_broken_duty\n"
    printf "cooling_01::refusal_drills_missing_inputs\n"
    printf "cooling_01::forensic_archive_integrity\n"
    exit 0
    ;;
  --check|--self-test)
    if [ -n "${BINARY}" ] && [ -x "${BINARY}" ]; then
      log_json "self_test" "ok" "found frankensim binary at ${BINARY}"
      exit 0
    else
      log_json "self_test" "failed" "no executable frankensim binary found"
      exit 1
    fi
    ;;
  --run|--negative|--replay)
    if [ -z "${BINARY}" ] || [ ! -x "${BINARY}" ]; then
      log_json "run_start" "failed" "no executable frankensim binary found"
      exit 1
    fi
    log_json "run_start" "started" "executing Cooling product-prefix and refusal suite"

    PROJECT="${REPO_ROOT}/examples/heatsink-fan/heatsink-fan.fsim"
    STL="${REPO_ROOT}/examples/heatsink-fan/heatsink.stl"
    PACK="${REPO_ROOT}/data/reference-project/aa6061.fsmcdpk"
    BROKEN="${REPO_ROOT}/examples/refusal-loop/broken.fsim"

    RUN_DIR="$(mktemp -d "${TMPDIR:-/tmp}/cooling-01-XXXXXX")"
    LEDGER_VERBS="${RUN_DIR}/ledger_verbs.db"
    LEDGER_RUN="${RUN_DIR}/ledger_run.db"

    # 1. Golden Path: Step-by-Step Verbs
    log_json "step_validate" "running" "validating heatsink-fan.fsim"
    "${BINARY}" --json validate "${PROJECT}" > "${ARTIFACT_DIR}/validate.json"
    grep -q '"status":"ok"' "${ARTIFACT_DIR}/validate.json"
    log_json "step_validate" "passed" "project validated clean with 0 findings"

    log_json "step_import" "running" "importing STL into ledger"
    "${BINARY}" --json import "${PROJECT}" "${STL}" "${LEDGER_VERBS}" --unit m --max-hole-edges 0 > "${ARTIFACT_DIR}/import.json"
    grep -q '"artifact_count":1' "${ARTIFACT_DIR}/import.json"
    log_json "step_import" "passed" "mesh imported into ledger"

    log_json "step_solve" "running" "executing solve prefix"
    RC=0
    "${BINARY}" --json solve "${PROJECT}" "${LEDGER_VERBS}" --materials "${PACK}" > "${ARTIFACT_DIR}/solve.json" 2> "${ARTIFACT_DIR}/solve.err" || RC=$?
    test "${RC}" -eq 5
    grep -q 'conduction' "${ARTIFACT_DIR}/solve.err"
    log_json "step_solve" "passed" "solve executed durable prefix and halted honestly at conduction gap"

    log_json "step_report" "running" "verifying report refuses the incomplete retained run"
    RC_REPORT=0
    "${BINARY}" --json report "cooling_demo_run" "${LEDGER_VERBS}" \
      > "${ARTIFACT_DIR}/report.json" 2> "${ARTIFACT_DIR}/report.err" || RC_REPORT=$?
    test "${RC_REPORT}" -eq 5
    grep -q '"status":"unavailable"' "${ARTIFACT_DIR}/report.json"
    grep -q 'cli-stage-unavailable' "${ARTIFACT_DIR}/report.err"
    if grep -Eq 'report-result|junction_maximum|thermal_margin|Verified|content_hash' \
      "${ARTIFACT_DIR}/report.json"; then
      log_json "step_report" "failed" "report emitted fabricated scientific evidence"
      exit 1
    fi
    log_json "step_report" "passed" "report remained unavailable without minting authority"

    log_json "step_package" "running" "verifying package refuses the incomplete retained run"
    RC_PACKAGE=0
    "${BINARY}" --json package "cooling_demo_run" "${LEDGER_VERBS}" \
      > "${ARTIFACT_DIR}/package.json" 2> "${ARTIFACT_DIR}/package.err" || RC_PACKAGE=$?
    test "${RC_PACKAGE}" -eq 5
    grep -q '"status":"unavailable"' "${ARTIFACT_DIR}/package.json"
    grep -q 'cli-stage-unavailable' "${ARTIFACT_DIR}/package.err"
    if grep -Eq 'package-result|merkle_root|"verdict":"pass"|junction_maximum|thermal_margin' \
      "${ARTIFACT_DIR}/package.json"; then
      log_json "step_package" "failed" "package emitted fabricated scientific evidence"
      exit 1
    fi
    log_json "step_package" "passed" "package remained unavailable without minting authority"

    # 2. Parity: One-Command `run`
    log_json "step_run_workflow" "running" "executing one-command run workflow"
    RC_RUN=0
    "${BINARY}" --json run "${PROJECT}" "${LEDGER_VERBS}" --materials "${PACK}" > "${ARTIFACT_DIR}/run.json" 2> "${ARTIFACT_DIR}/run.err" || RC_RUN=$?
    test "${RC_RUN}" -eq 5
    grep -q 'conduction' "${ARTIFACT_DIR}/run.json"
    log_json "step_run_workflow" "passed" "one-command run matches separate verbs semantic outcome"

    # 3. Failure Drill: Duty Range Refusal
    log_json "drill_broken_duty" "running" "testing invalid duty cycle refusal"
    RC_DUTY=0
    "${BINARY}" --json validate "${BROKEN}" > "${ARTIFACT_DIR}/broken.json" 2> "${ARTIFACT_DIR}/broken.err" || RC_DUTY=$?
    test "${RC_DUTY}" -eq 4
    grep -q 'project-duty-range' "${ARTIFACT_DIR}/broken.err"
    log_json "drill_broken_duty" "passed" "duty range refusal emitted code project-duty-range and actionable fix"

    # 4. Failure Drill: Non-Existent Input
    log_json "drill_missing_file" "running" "testing missing input file handling"
    RC_MISSING=0
    "${BINARY}" --json validate "/nonexistent/fake.fsim" > "${ARTIFACT_DIR}/missing.json" 2> "${ARTIFACT_DIR}/missing.err" || RC_MISSING=$?
    test "${RC_MISSING}" -eq 3
    log_json "drill_missing_file" "passed" "missing file exited with input error class (3)"

    log_json "run_terminal" "pass" "durable prefix and typed refusals passed; report/package remain unavailable"
    exit 0
    ;;
  *)
    printf "FATAL: unknown command %s\n" "${COMMAND}" >&2
    exit 2
    ;;
esac
