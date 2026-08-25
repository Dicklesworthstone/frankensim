#!/usr/bin/env bash
#
# cooling_compare.sh — Acceptance verification for frankensim compare (bead frankensim-extreal-program-f85xj.6.14.1).
#
# Usage:
#   scripts/e2e/cooling_compare.sh [--list|--check|--run]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMAND="${1:---run}"
ARTIFACT_DIR="${ARTIFACT_DIR:-${REPO_ROOT}/target/cooling-compare}"
mkdir -p "${ARTIFACT_DIR}"

BINARY="${FRANKENSIM_BIN:-}"
if [ -z "${BINARY}" ]; then
  if [ -f "/Volumes/USB_NVME/cargo-target/debug/frankensim" ]; then
    BINARY="/Volumes/USB_NVME/cargo-target/debug/frankensim"
  elif [ -f "${REPO_ROOT}/target/debug/frankensim" ]; then
    BINARY="${REPO_ROOT}/target/debug/frankensim"
  elif command -v frankensim >/dev/null 2>&1; then
    BINARY="$(command -v frankensim)"
  fi
fi

case "${COMMAND}" in
  --list)
    printf "cooling_compare::fail_closed_without_retained_loader\n"
    printf "cooling_compare::no_fabricated_qoi_authority\n"
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
    printf "==> 1. Verifying compare fails closed without a retained-run loader\n"
    set +e
    "${BINARY}" --json compare identical_run identical_run /definitely/missing/ledger.db \
      > "${ARTIFACT_DIR}/compare.json" 2> "${ARTIFACT_DIR}/compare.stderr.jsonl"
    compare_exit=$?
    set -e
    test "${compare_exit}" -eq 5
    grep -q '"status":"unavailable"' "${ARTIFACT_DIR}/compare.json"
    grep -q '"code":"cli-stage-unavailable"' "${ARTIFACT_DIR}/compare.stderr.jsonl"
    if grep -Eq 'junction_maximum|thermal_margin|evidence-aware-semantic-run-diff|"status":"changed"' \
      "${ARTIFACT_DIR}/compare.json"; then
      printf "FATAL: compare emitted fabricated semantic evidence\n" >&2
      exit 1
    fi

    printf "==> 2. Verifying Python SDK preserves the unavailable boundary\n"
    FRANKENSIM_BIN="${BINARY}" PYTHONPATH="${REPO_ROOT}/python" python3 -c '
from frankensim import FrankenSimClient, UnavailableError
client = FrankenSimClient()
try:
    client.compare("identical_run", "identical_run")
except UnavailableError:
    pass
else:
    raise AssertionError("compare must refuse until retained-run loading exists")
'

    printf "Cooling compare fail-closed checks passed; semantic comparison remains unavailable.\n"
    exit 0
    ;;
  *)
    printf "FATAL: unknown command %s\n" "${COMMAND}" >&2
    exit 2
    ;;
esac
