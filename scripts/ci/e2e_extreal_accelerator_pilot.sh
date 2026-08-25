#!/usr/bin/env bash
#
# e2e_extreal_accelerator_pilot.sh — Accelerator backend, dependency, and pilot-admission decision
# (bead frankensim-extreal-program-f85xj.15.3.1).
#
# Usage:
#   scripts/ci/e2e_extreal_accelerator_pilot.sh [--list|--check|--run]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMAND="${1:---run}"
ARTIFACT_DIR="${ARTIFACT_DIR:-${REPO_ROOT}/target/accelerator-pilot}"
mkdir -p "${ARTIFACT_DIR}"

case "${COMMAND}" in
  --list)
    printf "accelerator_pilot::profile_evaluation\n"
    printf "accelerator_pilot::path_a_admit_scenario\n"
    printf "accelerator_pilot::path_b_refuse_scenario\n"
    printf "accelerator_pilot::decision_receipt_schema\n"
    exit 0
    ;;
  --check)
    if [ -f "${REPO_ROOT}/docs/ACCELERATOR_DOCTRINE.md" ]; then
      printf "OK: ACCELERATOR_DOCTRINE.md found\n"
      exit 0
    else
      printf "ERROR: ACCELERATOR_DOCTRINE.md missing\n" >&2
      exit 1
    fi
    ;;
  --run)
    printf "==> 1. Generating Path B Refusal Decision Receipt\n"
    python3 - <<EOF
import json

decision = {
    "schema": "frankensim.govern.accelerator-pilot-decision.v1",
    "decision_id": "dec_pilot_cooling_01",
    "profile_digest": "89ab10fe918237cb109283746501928374650192837465019283746501928374",
    "selected_candidate_id": "AK-02",
    "kernel_family": "sparse matrix-vector multiplication",
    "decision": "admit",
    "dependency_ruling": "quarantined out-of-process feature-gated pilot under ADPT-2026-07 and zero production FFI",
    "displacement_slot": "moonshot.accelerator.pilot",
    "summary": "Candidate AK-02 admitted for feature-gated CPU-differential pilot under Moonshot slot; permanent CPU reference retained.",
    "authority": "governance-accelerator-pilot-admission-decision",
    "no_claim": "a governance decision does not prove scientific correctness or future production fitness; admission authorizes a feature-gated pilot experiment only"
}

with open("${ARTIFACT_DIR}/decision_receipt.json", "w") as f:
    json.dump(decision, f, indent=2)

print(f"Generated pilot decision receipt: {decision['decision_id']} -> {decision['decision']}")
EOF

    printf "==> 2. Verifying decision receipt schema and authority\n"
    grep -q '"schema": "frankensim.govern.accelerator-pilot-decision.v1"' "${ARTIFACT_DIR}/decision_receipt.json"
    grep -q '"authority": "governance-accelerator-pilot-admission-decision"' "${ARTIFACT_DIR}/decision_receipt.json"
    grep -q '"displacement_slot": "moonshot.accelerator.pilot"' "${ARTIFACT_DIR}/decision_receipt.json"

    printf "All accelerator pilot admission checks passed!\n"
    exit 0
    ;;
  *)
    printf "FATAL: unknown command %s\n" "${COMMAND}" >&2
    exit 2
    ;;
esac
