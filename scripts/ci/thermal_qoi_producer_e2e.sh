#!/usr/bin/env bash
#
# thermal_qoi_producer_e2e.sh — no-mock end-to-end proof for the canonical
# thermal QoI producer (bead frankensim-extreal-program-f85xj.5.10).
#
# This drives the REAL production path in fs-airflow: a solved fan/system
# operating point and a REAL fs-conduction FEM solve feed
# `extract_thermal_qois` and the operating-envelope audit. It never
# fabricates an evidence object: every asserted value is either produced by
# the public producer API or independently recomputed by the python oracle
# from the raw inputs (fan table, resistances, nodal temperatures) emitted
# alongside the results.
#
# THE LOAD-BEARING RULES:
#   1. The producer's own event stream is hash-chained (BLAKE3). A byte flip
#      anywhere breaks verification — proven live by the tamper drill, which
#      MUST fail through the real verifier.
#   2. The oracle recomputes junction maximum, thermal margin, spread bounds,
#      and the fan/network operating point from RAW inputs, not from the
#      producer's outputs, so a producer bug cannot grade its own homework.
#   3. All cargo work runs through RCH offloading; this script never builds
#      locally. Lane artifacts travel to/from the worker through stdout base64
#      (producer) and the repo-relative staged copy (verifier).
#
# Usage:
#   scripts/ci/thermal_qoi_producer_e2e.sh
#       [--profile pr|full|recovery]
#       [--artifact-dir PATH]
#       [--skip-tamper-drill]
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PROFILE="pr"
ARTIFACT_DIR=""
SKIP_TAMPER=0

die() { printf 'FATAL: %s\n' "$*" >&2; exit 2; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --profile) PROFILE="${2:?}"; shift 2 ;;
    --artifact-dir) ARTIFACT_DIR="${2:?}"; shift 2 ;;
    --skip-tamper-drill) SKIP_TAMPER=1; shift ;;
    *) die "unknown argument: $1" ;;
  esac
done
case "${PROFILE}" in pr|full|recovery) ;; *) die "unknown profile: ${PROFILE}" ;; esac

if [[ -z "${ARTIFACT_DIR}" ]]; then
  ARTIFACT_DIR="${REPO_ROOT}/target/tqoi-e2e-${PROFILE}"
fi
mkdir -p "${ARTIFACT_DIR}"

FAILURES=0
CHECKS=0

log() {
  printf '{"ts":"%s","lane":"thermal-qoi-producer-e2e","stage":"%s","detail":"%s"}\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" "$2" >&2
}

check() {
  CHECKS=$((CHECKS + 1))
  if [[ "$1" != "0" ]]; then
    FAILURES=$((FAILURES + 1))
    log check "FAIL: $2"
  else
    log check "pass: $2"
  fi
}

REVISION="$(git -C "${REPO_ROOT}" rev-parse HEAD 2>/dev/null || echo unknown)"
log setup "profile=${PROFILE} artifacts=${ARTIFACT_DIR} revision=${REVISION}"

# Lane events are staged INSIDE the repo tree (never target/, which RCH
# excludes from sync) so the verifier-mode run on the worker sees them.
LANE_DIR_REL=".tqoi-lane-events-${PROFILE}"
TAMPER_DIR_REL=".tqoi-tamper-events-${PROFILE}"
rm -rf "${REPO_ROOT:?}/${LANE_DIR_REL}" "${REPO_ROOT:?}/${TAMPER_DIR_REL}"
mkdir -p "${REPO_ROOT}/${LANE_DIR_REL}" "${REPO_ROOT}/${TAMPER_DIR_REL}"
trap 'rm -rf "${REPO_ROOT:?}/${LANE_DIR_REL}" "${REPO_ROOT:?}/${TAMPER_DIR_REL}"' EXIT

# TQOI_* variables cross to the worker ONLY through the explicit allowlist;
# the verifier reads its events from the staged repo-relative directory.
ALLOWLIST="RCH_ENV_ALLOWLIST=TQOI_MODE,TQOI_PROFILE,TQOI_SOURCE_REVISION,TQOI_ARTIFACT_DIR,TQOI_EMIT_STDOUT"

run_lane() {
  local mode="$1"
  local profile="$2"
  local artifact_rel="$3"
  local emit_stdout="$4"
  (
    cd "${REPO_ROOT}" &&
      env \
        "${ALLOWLIST}" \
        "TQOI_MODE=${mode}" \
        "TQOI_PROFILE=${profile}" \
        "TQOI_SOURCE_REVISION=${REVISION}" \
        "TQOI_ARTIFACT_DIR=${artifact_rel}" \
        "TQOI_EMIT_STDOUT=${emit_stdout}" \
        rch exec -- cargo test -p fs-airflow --test qoi_producer_e2e -- --nocapture --test-threads=1
  )
}

decode_stream() {
  local payload="$1"
  local destination="$2"
  python3 - "${payload}" "${destination}" <<'PY'
import base64, sys
with open(sys.argv[2], "wb") as handle:
    handle.write(base64.b64decode(sys.argv[1]))
PY
}

# ---------------------------------------------------------- phase 1: produce
log produce "running the real producer battery via RCH"
PRODUCE_OUT="${ARTIFACT_DIR}/produce.stdout"
PRODUCE_OK=0
if run_lane produce "${PROFILE}" "/tmp/tqoi-remote-produce-${PROFILE}" 1 \
    >"${PRODUCE_OUT}" 2>"${ARTIFACT_DIR}/produce.stderr"; then
  PRODUCE_OK=1
fi
check "${PRODUCE_OK}" "producer battery completed green"

B64_LINE=""
if [[ "${PRODUCE_OK}" == "1" ]]; then
  B64_LINE="$(grep -o 'TQOI_EVENTS_B64:[A-Za-z0-9+/=]*' "${PRODUCE_OUT}" | tail -1 || true)"
fi
[[ -n "${B64_LINE}" ]] || { log produce "no event stream on stdout; tail of stderr follows"; tail -40 "${ARTIFACT_DIR}/produce.stderr" >&2 || true; exit 1; }
decode_stream "${B64_LINE#TQOI_EVENTS_B64:}" "${ARTIFACT_DIR}/events.ndjson"

STREAM_NONEMPTY=0
[[ -s "${ARTIFACT_DIR}/events.ndjson" ]] && STREAM_NONEMPTY=1
check "${STREAM_NONEMPTY}" "decoded non-empty events.ndjson"
cp "${ARTIFACT_DIR}/events.ndjson" "${REPO_ROOT}/${LANE_DIR_REL}/events.ndjson"

# ----------------------------------------- phase 2: verifier, pristine stream
log verify "solver-free chain/coverage verification via RCH"
VERIFY_OK=0
if run_lane verify "${PROFILE}" "${LANE_DIR_REL}" 0 \
    >"${ARTIFACT_DIR}/verify.stdout" 2>"${ARTIFACT_DIR}/verify.stderr"; then
  VERIFY_OK=1
fi
check "${VERIFY_OK}" "pristine stream passes the real verifier"
VERIFY_CLEAN=0
if ! grep -qE "digest mismatch|prev chain mismatch|absent from stream|no successful run_end" \
  "${ARTIFACT_DIR}/verify.stderr" 2>/dev/null; then
  VERIFY_CLEAN=1
fi
check "${VERIFY_CLEAN}" "green verifier did not name any violation"

# --------------------------------------------- phase 3: python physics oracle
log oracle "independent recomputation from raw inputs"
ORACLE_OK=0
python3 - "${ARTIFACT_DIR}/events.ndjson" >"${ARTIFACT_DIR}/oracle.json" <<'PY' || ORACLE_OK=1
import json, struct, sys

def f64(text):
    return struct.unpack(">d", bytes.fromhex(text))[0]

events = [json.loads(line) for line in open(sys.argv[1]) if line.strip()]
inputs = next(event for event in events if event["event"] == "input_artifacts")
records = {event["qoi"]: event for event in events if event["event"] == "qoi_record"}
problems, notes = [], []

temps = [f64(t.strip('"')) for t in inputs["temperature_bits"]]
junction = json.loads(inputs["junction_vertices"])
raw_max = max(temps[index] for index in junction)
limit_k = f64(inputs["limit_k_bits"].strip('"'))
value_of = lambda name: f64(records[name]["value_bits"].strip('"'))

if value_of("junction_maximum").to_bytes(8, "big") != raw_max.to_bytes(8, "big"):
    problems.append(f"junction_maximum {value_of('junction_maximum')!r} != raw region max {raw_max!r}")
else:
    notes.append("junction_maximum equals raw region maximum bitwise")

margin_expected = limit_k - raw_max
if value_of("thermal_margin").to_bytes(8, "big") != margin_expected.to_bytes(8, "big"):
    problems.append(f"thermal_margin {value_of('thermal_margin')!r} != limit-max {margin_expected!r}")
else:
    notes.append("thermal_margin equals limit minus raw maximum bitwise")

field_min, field_max = min(temps), max(temps)
spread = value_of("uniformity_spread")
if not (0.0 <= spread <= field_max - field_min + 1e-12):
    problems.append(f"uniformity_spread {spread!r} outside [0, field range]")
else:
    notes.append("uniformity_spread within honest field-range bound")

mean = value_of("uniformity_mean")
if not (field_min - 1e-9 <= mean <= field_max + 1e-9):
    problems.append(f"uniformity_mean {mean!r} outside field range")
else:
    notes.append("uniformity_mean within [min(T), max(T)]")

resistances = [f64(r.strip('"')) for r in inputs["resistances"]]
r_series, r_leak = sum(resistances[:3]), resistances[3]
raw_points = json.loads(inputs["fan_points"])
fan_points = [(f64(p["q"].strip('"')), f64(p["p"].strip('"'))) for p in raw_points]

def fan_pressure(q):
    for (qa, pa), (qb, pb) in zip(fan_points, fan_points[1:]):
        if qa <= q <= qb:
            return pa + (pb - pa) * (q - qa) / (qb - qa)
    return None

def split_and_drop(q_total):
    lo, hi = 0.0, max(q_total, 0.0)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if r_series * mid * mid > r_leak * (q_total - mid) ** 2:
            hi = mid
        else:
            lo = mid
    branch = 0.5 * (lo + hi)
    return branch, r_series * branch * branch

flow_lo, flow_hi = fan_points[0][0], fan_points[-1][0]
oracle_flow = None
for _ in range(200):
    mid = 0.5 * (flow_lo + flow_hi)
    p_fan = fan_pressure(mid)
    if p_fan is None:
        break
    _, dp_sys = split_and_drop(mid)
    if p_fan > dp_sys:
        flow_lo = mid
    else:
        flow_hi = mid
oracle_flow = 0.5 * (flow_lo + flow_hi)
_, oracle_dp = split_and_drop(oracle_flow)

emitted_flow = f64(inputs["operating_flow_bits"].strip('"'))
emitted_dp = f64(inputs["operating_pressure_bits"].strip('"'))
if abs(oracle_flow - emitted_flow) > 0.01:
    problems.append(f"operating flow {emitted_flow!r} vs bisection oracle {oracle_flow!r} beyond tolerance")
else:
    notes.append(f"operating flow within declared tolerance of oracle ({oracle_flow:.6f} m^3/s)")
if abs(oracle_dp - emitted_dp) > max(1.0, 0.05 * abs(oracle_dp)):
    problems.append(f"pressure drop {emitted_dp!r} vs network oracle {oracle_dp!r} beyond 5%")
else:
    notes.append(f"pressure drop within 5% of network oracle ({oracle_dp:.3f} Pa)")

efficiency = f64(inputs["efficiency_bits"].strip('"'))
power_expected = emitted_dp * emitted_flow / efficiency
power_emitted = value_of("fan_power")
if abs(power_expected - power_emitted) > max(1e-6, 1e-6 * abs(power_expected)):
    problems.append(f"fan_power {power_emitted!r} != emitted dp*q/eta {power_expected!r}")
else:
    notes.append("fan_power equals emitted operating point times efficiency")

for name, record in records.items():
    terms = record["terms"] if isinstance(record["terms"], list) else json.loads(record["terms"])
    if len(terms) != 8:
        problems.append(f"{name}: {len(terms)} budget terms, expected exactly 8")

print(json.dumps({"problems": problems, "notes": notes}, indent=2))
sys.exit(1 if problems else 0)
PY
check "${ORACLE_OK}" "python oracle reproduces the producer from raw inputs"

# ------------------------------------------------------- phase 4: tamper drill
if [[ "${SKIP_TAMPER}" != "1" ]]; then
  log tamper "corrupting one digest byte; the verifier MUST refuse"
  python3 - "${REPO_ROOT}/${LANE_DIR_REL}/events.ndjson" \
    "${REPO_ROOT}/${TAMPER_DIR_REL}/events.ndjson" <<'PY'
import re, sys
lines = open(sys.argv[1]).read().splitlines()
victim = len(lines) // 2
match = re.search(r'"digest":"([0-9a-f])', lines[victim])
assert match, "no digest found to corrupt"
flipped = "0" if match.group(1) != "0" else "1"
lines[victim] = lines[victim][: match.start(1)] + flipped + lines[victim][match.end(1):]
open(sys.argv[2], "w").write("\n".join(lines) + "\n")
PY
  TAMPER_REFUSED=0
  if run_lane verify "${PROFILE}" "${TAMPER_DIR_REL}" 0 \
      >"${ARTIFACT_DIR}/tamper.stdout" 2>"${ARTIFACT_DIR}/tamper.stderr"; then
    TAMPER_REFUSED=0
  else
    TAMPER_REFUSED=1
  fi
  check "${TAMPER_REFUSED}" "tampered stream refused by the real verifier"
  NAMED_CHAIN=0
  if grep -q "digest mismatch" "${ARTIFACT_DIR}/tamper.stderr" 2>/dev/null; then
    NAMED_CHAIN=1
  fi
  check "${NAMED_CHAIN}" "refusal names the broken hash chain specifically"
else
  log tamper "skipped by flag"
fi

# ------------------------------------------- phase 5: deterministic repeat
if [[ "${PROFILE}" == "full" ]]; then
  log repeat "second full produce; artifact bytes must be identical"
  REPEAT_OUT="${ARTIFACT_DIR}/produce_repeat.stdout"
  REPEAT_OK=0
  if run_lane produce "${PROFILE}" "/tmp/tqoi-remote-repeat-${PROFILE}" 1 \
      >"${REPEAT_OUT}" 2>"${ARTIFACT_DIR}/produce_repeat.stderr"; then
    REPEAT_OK=1
  fi
  check "${REPEAT_OK}" "repeat battery completed green"
  if [[ "${REPEAT_OK}" == "1" ]]; then
    B64_REPEAT="$(grep -o 'TQOI_EVENTS_B64:[A-Za-z0-9+/=]*' "${REPEAT_OUT}" | tail -1 || true)"
    decode_stream "${B64_REPEAT#TQOI_EVENTS_B64:}" "${ARTIFACT_DIR}/events_repeat.ndjson"
    IDENTICAL=0
    cmp -s "${ARTIFACT_DIR}/events.ndjson" "${ARTIFACT_DIR}/events_repeat.ndjson" && IDENTICAL=1
    check "${IDENTICAL}" "repeat produce is byte-identical across runs"
  fi
fi

# ------------------------------------------------------------------ summary
cat > "${ARTIFACT_DIR}/summary.json" <<JSON
{
  "bead": "frankensim-extreal-program-f85xj.5.10",
  "profile": "${PROFILE}",
  "revision": "${REVISION}",
  "checks": ${CHECKS},
  "failures": ${FAILURES},
  "events": "${ARTIFACT_DIR}/events.ndjson",
  "oracle": "raw-input recomputation (max/margin/spread-bounds/network-bisection/power)"
}
JSON

log summary "checks=${CHECKS} failures=${FAILURES}"
log summary "artifacts written to ${ARTIFACT_DIR}"
if [[ "${FAILURES}" -ne 0 ]]; then
  printf 'FAIL: %d of %d checks failed (profile %s)\n' "${FAILURES}" "${CHECKS}" "${PROFILE}"
  exit 1
fi
printf 'OK: %d checks passed; thermal QoI producer lane green (profile %s)\n' "${CHECKS}" "${PROFILE}"
