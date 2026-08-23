#!/usr/bin/env bash
#
# examples_freshness_e2e.sh — every worked example keeps running, or this
# lane breaks (bead frankensim-extreal-program-f85xj.6.12).
#
# Examples that rot are worse than none. This harness executes each worked
# example's documented commands against the REAL `frankensim` binary:
#
#   1. examples/heated-plate   — minimal schema walkthrough validates clean.
#   2. data/reference-project  — the cooling reference fixture validates
#                                clean (the enclosure example's subject).
#   3. examples/refusal-loop   — broken.fsim keeps refusing with exactly
#                                code `project-duty-range`, and its one-
#                                token repair stays byte-equal to the
#                                tracked reference project.
#
# FROZEN BYTES: the canonical project hashes are frozen as literals in the
# G0 battery (`crates/fs-cli/tests/cli.rs`,
# g0_the_worked_example_fixtures_stay_fresh_through_the_real_cli_verb),
# which runs wherever `cargo test` runs and fails on any fixture drift.
# This lane is the human-runnable wrapper; it needs a NATIVE frankensim
# binary. Under the RCH offload regime, set FRANKENSIM_BIN (or --binary)
# explicitly, or rely on the G0 battery, which needs no local binary.
#
# Usage:
#   scripts/ci/examples_freshness_e2e.sh [--binary PATH]
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BINARY="${FRANKENSIM_BIN:-}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --binary) BINARY="${2:-}"; shift 2 ;;
    -h|--help) sed -n '3,26p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) printf 'FATAL: unknown argument: %s\n' "$1" >&2; exit 2 ;;
  esac
done

HEATED="${REPO_ROOT}/examples/heated-plate/heated-plate.fsim"
REFERENCE="${REPO_ROOT}/data/reference-project/cooling-reference.fsim"
BROKEN="${REPO_ROOT}/examples/refusal-loop/broken.fsim"

for f in "${HEATED}" "${REFERENCE}" "${BROKEN}"; do
  [[ -f "${f}" ]] || { printf 'FATAL: missing %s\n' "${f}" >&2; exit 2; }
done

if [[ -z "${BINARY}" ]]; then
  printf 'FATAL: no frankensim binary. Set FRANKENSIM_BIN or pass --binary PATH.\n' >&2
  printf 'The frozen-hash freshness assertions live in the fs-cli G0 battery,\n' >&2
  printf 'which runs under plain `cargo test -p fs-cli --test cli` anywhere.\n' >&2
  exit 2
fi
[[ -x "${BINARY}" ]] || { printf 'FATAL: not executable: %s\n' "${BINARY}" >&2; exit 2; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/examples-freshness-XXXXXX")"
trap 'rm -rf "${WORK}"' EXIT
LOG="${WORK}/run.ndjson"
: > "${LOG}"

FAILURES=0
CHECKS=0

log() {
  local kind="$1"; shift
  local msg="$1"; shift
  printf '{"schema":"frankensim.ci.examples-freshness.v1","kind":"%s","message":%s}\n' \
    "${kind}" "$(printf '%s' "${msg}" | python3 -c 'import json,sys; print(json.dumps(sys.stdin.read()))')" >> "${LOG}"
  printf '[%-6s] %s\n' "${kind}" "${msg}" >&2
}

check() {
  local desc="$1"; shift
  CHECKS=$((CHECKS + 1))
  if "$@"; then
    log check "PASS ${desc}"
  else
    FAILURES=$((FAILURES + 1))
    log check "FAIL ${desc}"
  fi
}

validate_ok() {
  "${BINARY}" --json validate "$1" > "${WORK}/v.json" 2> "${WORK}/v.err"
}

# ---- 1. heated plate: minimal example validates clean ----------------------
check "heated-plate validates ok" validate_ok "${HEATED}"
check "heated-plate reports zero findings" \
  grep -q '"finding_count":0' "${WORK}/v.json"

# ---- 2. cooling reference: enclosure example's subject stays valid --------
check "cooling-reference validates ok" validate_ok "${REFERENCE}"
check "cooling-reference reports zero findings" \
  grep -q '"finding_count":0' "${WORK}/v.json"

# ---- 3. refusal loop: broken fixture refuses with the documented code -----
RC=0
"${BINARY}" --json validate "${BROKEN}" > "${WORK}/b.json" 2> "${WORK}/b.err" || RC=$?
check "broken.fsim exits nonzero (observed rc=${RC})" test "${RC}" -ne 0
check "refusal names project-duty-range" grep -q 'project-duty-range' "${WORK}/b.err"
check "refusal states the duty fix" grep -q 'duty must lie in 0.0..=1.0' "${WORK}/b.err"

# ---- 4. one-token repair stays byte-equal to the reference -----------------
sed 's/:duty 2\.0/:duty 1.0/' "${BROKEN}" > "${WORK}/repaired.fsim"
check "one-token repair reproduces the tracked reference bytes" \
  cmp -s "${WORK}/repaired.fsim" "${REFERENCE}"

# ------------------------------------------------------------------- verdict
log summary "checks=${CHECKS} failures=${FAILURES}"
if [[ "${FAILURES}" -gt 0 ]]; then
  printf 'FAILED: %d of %d freshness checks failed; full NDJSON log: %s\n' \
    "${FAILURES}" "${CHECKS}" "${LOG}" >&2
  exit 1
fi
printf 'OK: all %d examples-freshness checks passed\n' "${CHECKS}"
