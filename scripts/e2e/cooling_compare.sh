#!/usr/bin/env bash
#
# cooling_compare.sh — no-mock end-to-end verification of `frankensim compare`
# (bead frankensim-rc-root-q61wp.47): the verb diffs two completed runs by
# their retained receipts, refuses without a ledger, and reports an empty
# comparison for a run against itself.
#
# Usage:
#   scripts/e2e/cooling_compare.sh [--list|--check|--run]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMAND="${1:---run}"
ARTIFACT_DIR="${ARTIFACT_DIR:-${REPO_ROOT}/target/cooling-compare}"
FSIM="${REPO_ROOT}/data/reference-project/cooling-reference.fsim"
STL="${REPO_ROOT}/data/reference-project/plate.stl"
PACK="${REPO_ROOT}/data/reference-project/aa6061.fsmcdpk"

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
    printf "cooling_compare::missing_ledger_refuses_before_reading\n"
    printf "cooling_compare::self_compare_of_a_retained_run_is_empty\n"
    printf "cooling_compare::python_sdk_reads_the_comparison\n"
    exit 0
    ;;
  --check)
    if [ -n "${BINARY}" ] && [ -x "${BINARY}" ]; then
      printf "OK: frankensim binary available at %s\n" "${BINARY}"
      exit 0
    else
      printf "ERROR: frankensim binary not found\n" >&2
      exit 1
    fi
    ;;
  --run)
    if [ -z "${BINARY}" ] || [ ! -x "${BINARY}" ]; then
      printf "ERROR: frankensim binary not found\n" >&2
      exit 1
    fi
    mkdir -p "${ARTIFACT_DIR}"
    ARTIFACT_DIR="$(cd "${ARTIFACT_DIR}" && pwd)"
    LEDGER="${ARTIFACT_DIR}/compare.db"
    rm -f "${LEDGER}"

    printf "==> 1. Verifying compare refuses a missing ledger before reading anything\n"
    set +e
    "${BINARY}" --json compare identical_run identical_run /definitely/missing/ledger.db \
      > "${ARTIFACT_DIR}/compare-missing.json" 2> "${ARTIFACT_DIR}/compare-missing.stderr.jsonl"
    compare_exit=$?
    set -e
    test "${compare_exit}" -eq 3
    grep -q '"status":"refused"' "${ARTIFACT_DIR}/compare-missing.json"
    grep -q '"code":"cli-export-ledger-missing"' "${ARTIFACT_DIR}/compare-missing.stderr.jsonl"
    if grep -Eq '"changed":|"qoi_diffs"' "${ARTIFACT_DIR}/compare-missing.json"; then
      printf "FATAL: compare emitted comparison rows without a ledger\n" >&2
      exit 1
    fi

    printf "==> 2. Importing and running the reference project (seven stages)\n"
    "${BINARY}" --json import "${FSIM}" "${STL}" "${LEDGER}" --unit m --max-hole-edges 0 \
      > "${ARTIFACT_DIR}/import.json" 2> "${ARTIFACT_DIR}/import.stderr.jsonl"
    "${BINARY}" --json run "${FSIM}" "${LEDGER}" --materials "${PACK}" \
      > "${ARTIFACT_DIR}/run.json" 2> "${ARTIFACT_DIR}/run.stderr.jsonl"
    grep -q '"stages_completed":7' "${ARTIFACT_DIR}/run.json"
    RUN_ID="$(sed -n 's/.*"run":"\([0-9a-f]\{64\}\)".*/\1/p' "${ARTIFACT_DIR}/run.json" | head -n 1)"
    test "${#RUN_ID}" -eq 64

    printf "==> 3. Comparing the run with itself: every row must be unchanged\n"
    "${BINARY}" --json compare "${RUN_ID}" "${RUN_ID}" "${LEDGER}" \
      > "${ARTIFACT_DIR}/compare-self.json" 2> "${ARTIFACT_DIR}/compare-self.stderr.jsonl"
    grep -q '"command":"compare","status":"ok"' "${ARTIFACT_DIR}/compare-self.json"
    grep -q '"changed":false' "${ARTIFACT_DIR}/compare-self.json"
    grep -q '"summary":"identical runs: no differences in any retained receipt"' "${ARTIFACT_DIR}/compare-self.json"
    grep -q '"qoi_count":1' "${ARTIFACT_DIR}/compare-self.json"
    grep -q '"delta":0,"rel_delta":0' "${ARTIFACT_DIR}/compare-self.json"
    grep -q '"authority":"projection-of-retained-receipts"' "${ARTIFACT_DIR}/compare-self.json"
    test "$(grep -o '"status":"unchanged (same receipt)"' "${ARTIFACT_DIR}/compare-self.json" | wc -l | tr -d ' ')" -eq 7
    if grep -q '"status":"changed"' "${ARTIFACT_DIR}/compare-self.json"; then
      printf "FATAL: self-compare reported a changed stage\n" >&2
      exit 1
    fi

    printf "==> 4. Python SDK reads the same comparison\n"
    FRANKENSIM_BIN="${BINARY}" PYTHONPATH="${REPO_ROOT}/python" \
      FS_COMPARE_RUN="${RUN_ID}" FS_COMPARE_LEDGER="${LEDGER}" python3 -c '
import os
from frankensim import FrankenSimClient
client = FrankenSimClient()
comp = client.compare(os.environ["FS_COMPARE_RUN"], os.environ["FS_COMPARE_RUN"], ledger_path=os.environ["FS_COMPARE_LEDGER"], strict=True)
assert comp.status == "ok", comp
assert comp.changed is False, comp
assert comp.qoi_count == 1 and comp.qoi_diffs[0].delta == 0.0, comp
assert comp.authority == "projection-of-retained-receipts", comp
'

    printf "Cooling compare checks passed: refusal without a ledger, empty self-comparison, SDK round trip (run %s).\n" "${RUN_ID}"
    exit 0
    ;;
  *)
    printf "FATAL: unknown command %s\n" "${COMMAND}" >&2
    exit 2
    ;;
esac
