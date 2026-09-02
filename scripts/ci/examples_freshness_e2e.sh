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
#   4. examples/heatsink-fan   — the finned heatsink validates, imports,
#                                SOLVES all seven stages (conduction with
#                                a derived airflow-convection law), and
#                                the report/package verbs export exactly
#                                the retained bytes of that run.
#   5. heatsink-fan-ladder      — the same project with solver fidelity
#                                "ladder": three uniform 1->8 rungs, the
#                                QoI budget's discretization term measured
#                                (interval), the report's convergence
#                                section present. Minutes in a debug build.
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

# ---- 5. heatsink-fan: full cooling contract validates clean -----------------
HEATSINK="${REPO_ROOT}/examples/heatsink-fan/heatsink-fan.fsim"
STL="${REPO_ROOT}/examples/heatsink-fan/heatsink.stl"
PACK="${REPO_ROOT}/data/reference-project/aa6061.fsmcdpk"

check "heatsink-fan validates ok" validate_ok "${HEATSINK}"
check "heatsink-fan reports zero findings" \
  grep -q '"finding_count":0' "${WORK}/v.json"

# ---- 6. import and solve orchestration --------------------------------------
import_ok() {
  "${BINARY}" --json import "${HEATSINK}" "${STL}" "${WORK}/ledger.db" --unit m --max-hole-edges 0 > "${WORK}/imp.json" 2> "${WORK}/imp.err"
}
check "heatsink import into ledger ok" import_ok
check "import reports artifact count >= 1" grep -q '"artifact_count":1' "${WORK}/imp.json"

solve_completes() {
  "${BINARY}" --json solve "${HEATSINK}" "${WORK}/ledger.db" --materials "${PACK}" > "${WORK}/s.json" 2> "${WORK}/s.err"
}
check "solve completes every stage on the finned heatsink (exit 0)" solve_completes
check "solve reports seven completed stages" grep -q '"stages_completed":7' "${WORK}/s.json"
check "conduction stage executed (derived airflow convection, not a typed gap)" \
  grep -q '"stage":"conduction","ordinal":4,"status":"ok"' "${WORK}/s.err"

# ---- 7. report and package are projections of the retained run --------------
RUN_ID="$(grep -oE '"run":"[0-9a-f]{64}"' "${WORK}/s.json" | head -1 | cut -d'"' -f4)"
check "solve result names its 64-hex run id" test "${#RUN_ID}" -eq 64

report_ok() {
  (cd "${WORK}" && "${BINARY}" --json report "${RUN_ID}" "${WORK}/ledger.db" > "${WORK}/rep.json" 2> "${WORK}/rep.err")
}
check "report exports the retained HTML and JSON twin" report_ok
check "report export names the retained content hash" grep -q '"content_hash":"' "${WORK}/rep.json"
check "report verdict is the retained Estimated/indeterminate one" grep -q '"verdict":"indeterminate"' "${WORK}/rep.json"

package_ok() {
  (cd "${WORK}" && "${BINARY}" --json package "${RUN_ID}" "${WORK}/ledger.db" > "${WORK}/pkg.json" 2> "${WORK}/pkg.err")
}
check "package exports the retained evidence package" package_ok
check "package passes the checker" grep -q '"checker":"pass"' "${WORK}/pkg.json"

unknown_run_refuses() {
  RC=0
  "${BINARY}" --json report "0000000000000000000000000000000000000000000000000000000000000000" "${WORK}/ledger.db" > "${WORK}/unk.json" 2> "${WORK}/unk.err" || RC=$?
  test "${RC}" -eq 4 && grep -q 'cli-solve-unknown-run' "${WORK}/unk.err"
}
check "report of an unknown run refuses with cli-solve-unknown-run (exit 4)" unknown_run_refuses

# ---- 8. the ladder variant: three uniform rungs and a measured discretization term
LADDER="${REPO_ROOT}/examples/heatsink-fan/heatsink-fan-ladder.fsim"
check "heatsink-fan-ladder validates ok" validate_ok "${LADDER}"
ladder_import_ok() {
  "${BINARY}" --json import "${LADDER}" "${STL}" "${WORK}/ladder.db" --unit m --max-hole-edges 0 > "${WORK}/limp.json" 2> "${WORK}/limp.err"
}
check "ladder variant imports into its own ledger" ladder_import_ok
ladder_solve_completes() {
  "${BINARY}" --json solve "${LADDER}" "${WORK}/ladder.db" --materials "${PACK}" > "${WORK}/ls.json" 2> "${WORK}/ls.err"
}
check "ladder solve completes every stage (three uniform rungs; minutes in a debug build)" ladder_solve_completes
check "ladder solve reports seven completed stages" grep -q '"stages_completed":7' "${WORK}/ls.json"
check "QoI stage measured exactly one budget term" grep -q '"budget_terms_measured":1' "${WORK}/ls.err"
check "the weakest term is now seven-no-data, not all-eight" grep -q '"weakest_term":"seven-no-data"' "${WORK}/ls.err"
LADDER_RUN="$(grep -oE '"run":"[0-9a-f]{64}"' "${WORK}/ls.json" | head -1 | cut -d'"' -f4)"
ladder_report_ok() {
  (cd "${WORK}" && "${BINARY}" --json report "${LADDER_RUN}" "${WORK}/ladder.db" > "${WORK}/lrep.json" 2> "${WORK}/lrep.err")
}
check "ladder report exports" ladder_report_ok
check "report JSON twin carries the interval discretization term" \
  grep -q '"state": "interval"' "${WORK}/${LADDER_RUN}.report.json"
check "report JSON twin names the grid-refinement method" \
  grep -Eq 'richardson-gci|eca-hoekstra-data-range|bitwise-agreement' "${WORK}/${LADDER_RUN}.report.json"
check "report JSON twin has a convergence section" grep -q '"convergence"' "${WORK}/${LADDER_RUN}.report.json"
check "ladder verdict stays the honest Estimated/indeterminate one" grep -q '"verdict":"indeterminate"' "${WORK}/lrep.json"

# ------------------------------------------------------------------- verdict
log summary "checks=${CHECKS} failures=${FAILURES}"
if [[ "${FAILURES}" -gt 0 ]]; then
  printf 'FAILED: %d of %d freshness checks failed; full NDJSON log: %s\n' \
    "${FAILURES}" "${CHECKS}" "${LOG}" >&2
  exit 1
fi
printf 'OK: all %d examples-freshness checks passed\n' "${CHECKS}"
