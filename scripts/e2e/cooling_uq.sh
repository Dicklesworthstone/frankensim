#!/usr/bin/env bash
#
# cooling_uq.sh — product admission for empirical Cooling UQ.
#
# Usage:
#   scripts/e2e/cooling_uq.sh [--list|--check|--self-test|--run|--negative|--replay]

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"

COMMAND="${1:---run}"
BASE="examples/cooling-network/fan-correlated-hotspot.json"
PLAN="examples/cooling-network/uq-fan-hotspot.json"
NEGATIVE="examples/cooling-network/uq-unknown-dependence.json"

run_uq() {
  cargo run --quiet -p fs-cli --bin frankensim -- --json cooling-network-uq "${BASE}" "${PLAN}"
}

case "${COMMAND}" in
  --list)
    printf "cooling_uq::real_coupled_sampling\n"
    printf "cooling_uq::explicit_dependence\n"
    printf "cooling_uq::unknown_correlation_refusal\n"
    printf "cooling_uq::empirical_compliance_probability\n"
    printf "cooling_uq::deterministic_replay\n"
    exit 0
    ;;
  --check)
    cargo check -p fs-cli --bin frankensim
    ;;
  --self-test)
    cargo test -p fs-cli --bin frankensim uq_command
    cargo test -p fs-cli --test cooling_network_uq
    ;;
  --run)
    run_uq
    ;;
  --negative)
    if output="$(cargo run --quiet -p fs-cli --bin frankensim -- --json cooling-network-uq "${BASE}" "${NEGATIVE}" 2>&1)"; then
      printf '%s\n' "${output}" >&2
      printf '%s\n' 'FATAL: unknown multivariate dependence unexpectedly passed' >&2
      exit 1
    fi
    case "${output}" in
      *"explicit dependence"*|*"joint probability"*) printf '%s\n' "${output}" ;;
      *)
        printf '%s\n' "${output}" >&2
        printf '%s\n' 'FATAL: negative UQ case failed for the wrong reason' >&2
        exit 1
        ;;
    esac
    ;;
  --replay)
    first="$(run_uq)"
    second="$(run_uq)"
    if [[ "${first}" != "${second}" ]]; then
      printf '%s\n' 'FATAL: deterministic Cooling UQ replay diverged' >&2
      exit 1
    fi
    printf '%s\n' "${first}"
    ;;
  *)
    printf "FATAL: unknown command %s\n" "${COMMAND}" >&2
    exit 2
    ;;
esac
