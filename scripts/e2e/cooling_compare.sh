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
    printf "cooling_compare::text_diff\n"
    printf "cooling_compare::json_diff\n"
    printf "cooling_compare::qoi_semantics\n"
    printf "cooling_compare::evidence_regime\n"
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
    printf "==> 1. Testing text compare output\n"
    "${BINARY}" compare baseline_run candidate_run > "${ARTIFACT_DIR}/compare.txt"
    grep -q "FrankenSim Semantic Run Comparison" "${ARTIFACT_DIR}/compare.txt"
    grep -q "junction_maximum" "${ARTIFACT_DIR}/compare.txt"
    grep -q "thermal_margin" "${ARTIFACT_DIR}/compare.txt"

    printf "==> 2. Testing JSON structured compare output\n"
    "${BINARY}" --json compare baseline_run candidate_run > "${ARTIFACT_DIR}/compare.json"
    grep -q '"schema": "frankensim.cli.compare-result.v1"' "${ARTIFACT_DIR}/compare.json"
    grep -q '"status": "changed"' "${ARTIFACT_DIR}/compare.json"
    grep -q '"name": "junction_maximum"' "${ARTIFACT_DIR}/compare.json"
    grep -q '"name": "thermal_margin"' "${ARTIFACT_DIR}/compare.json"
    grep -q '"authority": "evidence-aware-semantic-run-diff"' "${ARTIFACT_DIR}/compare.json"

    printf "==> 3. Testing Python SDK compare integration\n"
    PYTHONPATH="${REPO_ROOT}/python" python3 -c '
from frankensim import FrankenSimClient
client = FrankenSimClient()
res = client.compare("baseline_run", "candidate_run")
assert res.status == "changed"
assert len(res.qoi_diffs) >= 2
print(f"Verified {len(res.qoi_diffs)} QoI diffs via Python SDK")
'

    printf "All cooling_compare acceptance checks passed!\n"
    exit 0
    ;;
  *)
    printf "FATAL: unknown command %s\n" "${COMMAND}" >&2
    exit 2
    ;;
esac
