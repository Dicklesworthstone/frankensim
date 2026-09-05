#!/usr/bin/env bash
#
# cooling_uq.sh — admission for the unfinished Cooling UQ product lane
# (bead frankensim-extreal-program-f85xj.6.7). The fs-uq product_plan library
# tests do not exercise project inputs, Cooling solves, or product outputs.
#
# Usage:
#   scripts/e2e/cooling_uq.sh [--list|--check|--self-test|--run|--negative|--replay]

set -euo pipefail

COMMAND="${1:---run}"

case "${COMMAND}" in
  --list)
    printf '%s\n' '# Required cases; the product UQ driver is not implemented.'
    printf "cooling_uq::gaussian_sampling\n"
    printf "cooling_uq::epistemic_bounding\n"
    printf "cooling_uq::unknown_correlation_refusal\n"
    printf "cooling_uq::compliance_probability\n"
    printf "cooling_uq::determinism_audit\n"
    exit 0
    ;;
  --check|--self-test|--run|--negative|--replay)
    printf '%s\n' '{"event":"run_terminal","status":"unavailable","code":"cooling-uq-product-driver-unavailable","detail":"No project-to-Cooling-QoI uncertainty driver is wired. Library product_plan tests are focused proof only; no product sampling, probability, budget truncation, negative case, or replay ran.","fix":"Implement frankensim-extreal-program-f85xj.6.7 and invoke the real binary, solver, ledger, report/package and checker before admitting this lane."}'
    exit 5
    ;;
  *)
    printf "FATAL: unknown command %s\n" "${COMMAND}" >&2
    exit 2
    ;;
esac
