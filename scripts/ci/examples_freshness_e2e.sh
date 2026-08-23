#!/usr/bin/env bash
#
# examples_freshness_e2e.sh — every worked example keeps running, or this
# lane breaks (bead frankensim-extreal-program-f85xj.6.12).
#
# Examples that rot are worse than none. This harness executes each worked
# example's documented commands against the REAL `frankensim` binary and
# compares the observed results against frozen expectations:
#
#   1. examples/heated-plate   — minimal schema walkthrough validates clean,
#                                byte-for-byte the frozen canonical hash.
#   2. data/reference-project  — the cooling reference fixture validates
#                                clean (the enclosure example's subject).
#   3. examples/refusal-loop   — broken.fsim must keep refusing with exactly
#                                code `project-duty-range`, and its one-token
#                                repair must stay byte-equal to the tracked
#                                reference project.
#
# THE LOAD-BEARING RULE: an expectation that stops matching is a FAILURE,
# never a silent pass. If you intentionally change a fixture, regenerate its
# frozen hash with --update and review the diff; the commit that changes a
# fixture and the commit that changes its frozen hash must be the same.
#
# Usage:
#   scripts/ci/examples_freshness_e2e.sh [--binary PATH] [--update]
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXPECTED="${REPO_ROOT}/scripts/ci/examples-freshness-expected.json"
BINARY="${FRANKENSIM_BIN:-}"
UPDATE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --binary) BINARY="${2:-}"; shift 2 ;;
    --update) UPDATE=1; shift ;;
    -h|--help) sed -n '3,26p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) printf 'FATAL: unknown argument: %s\n' "$1" >&2; exit 2 ;;
  esac
done

HEATED="${REPO_ROOT}/examples/heated-plate/heated-plate.fsim"
REFERENCE="${REPO_ROOT}/data/reference-project/cooling-reference.fsim"
BROKEN="${REPO_ROOT}/examples/refusal-loop/broken.fsim"

for f in "${HEATED}" "${REFERENCE}" "${BROKEN}" "${EXPECTED}"; do
  [[ -f "${f}" ]] || { printf 'FATAL: missing %s\n' "${f}" >&2; exit 2; }
done

if [[ -z "${BINARY}" ]]; then
  cargo build -q -p fs-cli --bin frankensim
  BINARY="$(find "${CARGO_TARGET_DIR:-${REPO_ROOT}/target}" -name frankensim -type f -perm -u+x 2>/dev/null | head -1)"
fi
[[ -x "${BINARY}" ]] || { printf 'FATAL: frankensim binary not found: %s\n' "${BINARY}" >&2; exit 2; }

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

frozen_hash() {
  python3 - "$1" "$2" <<'PY'
import json, sys
path, field = sys.argv[1], sys.argv[2]
data = json.load(open(path))
print(data["fixtures"][field]["project_hash"])
PY
}

observed_hash() {
  python3 -c 'import json,sys; data=json.load(sys.stdin); print(data["project_hash"])'
}

validate_ok() {
  "${BINARY}" --json validate "$1" > "${WORK}/v.json" 2> "${WORK}/v.err"
}

# ---- 1. heated plate: minimal example validates clean at the frozen hash ---
check "heated-plate validates ok" validate_ok "${HEATED}"
check "heated-plate reports zero findings" \
  grep -q '"finding_count":0' "${WORK}/v.json"
OBS_HEATED="$(observed_hash < "${WORK}/v.json")"
if [[ "${UPDATE}" != 1 ]]; then
  FROZEN="$(frozen_hash "${EXPECTED}" "heated_plate")"
  check "heated-plate canonical hash matches frozen expectation (observed=${OBS_HEATED} frozen=${FROZEN})" \
    test "${OBS_HEATED}" = "${FROZEN}"
fi

# ---- 2. cooling reference: enclosure example's subject stays valid --------
check "cooling-reference validates ok" validate_ok "${REFERENCE}"
check "cooling-reference reports zero findings" \
  grep -q '"finding_count":0' "${WORK}/v.json"
OBS_REFERENCE="$(observed_hash < "${WORK}/v.json")"
if [[ "${UPDATE}" != 1 ]]; then
  FROZEN="$(frozen_hash "${EXPECTED}" "cooling_reference")"
  check "cooling-reference canonical hash matches frozen expectation (observed=${OBS_REFERENCE} frozen=${FROZEN})" \
    test "${OBS_REFERENCE}" = "${FROZEN}"
fi

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

if [[ "${UPDATE}" == 1 ]]; then
  python3 - "$EXPECTED" "${OBS_HEATED}" "${OBS_REFERENCE}" <<'PY'
import json, sys
path, heated, reference = sys.argv[1], sys.argv[2], sys.argv[3]
data = {"fixtures": {
    "heated_plate": {"project_hash": heated},
    "cooling_reference": {"project_hash": reference}}}
with open(path, "w") as f:
    json.dump(data, f, indent=2, sort_keys=True)
    f.write("\n")
print(f"updated {path}")
PY
fi


# ------------------------------------------------------------------- verdict
log summary "checks=${CHECKS} failures=${FAILURES}"
if [[ "${FAILURES}" -gt 0 ]]; then
  printf 'FAILED: %d of %d freshness checks failed; full NDJSON log: %s\n' \
    "${FAILURES}" "${CHECKS}" "${LOG}" >&2
  exit 1
fi
printf 'OK: all %d examples-freshness checks passed\n' "${CHECKS}"
