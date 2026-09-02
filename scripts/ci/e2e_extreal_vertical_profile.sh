#!/usr/bin/env bash
#
# e2e_extreal_vertical_profile.sh — End-to-end vertical profiling: wall-clock of
# the CLI phases (validate, import, solve, report, package) on one host
# (bead frankensim-extreal-program-f85xj.15.2). Kernel attribution, energy and
# memory are NOT measured by this lane and are recorded as absent, never
# estimated (bridge plan A9 / q61wp.60); the accelerator falsifier therefore
# refuses with evidence until a named profiler supplies measured kernel rows.
#
# Usage:
#   scripts/ci/e2e_extreal_vertical_profile.sh [--list|--check|--run]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMAND="${1:---run}"
ARTIFACT_DIR="${ARTIFACT_DIR:-${REPO_ROOT}/target/vertical-profile}"
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
    printf "extreal_profile::stage_attribution\n"
    printf "extreal_profile::kernel_roofline_positions\n"
    printf "extreal_profile::accelerator_falsifier_check\n"
    printf "extreal_profile::decision_receipt_schema\n"
    exit 0
    ;;
  --check)
    if [ -n "${BINARY}" ] && [ -x "${BINARY}" ]; then
      printf "OK: frankensim binary found at %s\n" "${BINARY}"
      exit 0
    else
      printf "ERROR: frankensim binary not found\n" >&2
      exit 1
    fi
    ;;
  --run)
    printf "==> 1. Running end-to-end pipeline profiling campaign\n"
    PROJECT="${REPO_ROOT}/examples/heatsink-fan/heatsink-fan.fsim"
    STL="${REPO_ROOT}/examples/heatsink-fan/heatsink.stl"
    PACK="${REPO_ROOT}/data/reference-project/aa6061.fsmcdpk"

    RUN_DIR="$(mktemp -d "${TMPDIR:-/tmp}/profile-run-XXXXXX")"
    LEDGER="${RUN_DIR}/profile_ledger.db"

    # Step 1: Validate
    T0=$(python3 -c 'import time; print(time.perf_counter())')
    "${BINARY}" --json validate "${PROJECT}" > "${ARTIFACT_DIR}/validate.json"
    T1=$(python3 -c 'import time; print(time.perf_counter())')

    # Step 2: Import
    "${BINARY}" --json import "${PROJECT}" "${STL}" "${LEDGER}" --unit m --max-hole-edges 0 > "${ARTIFACT_DIR}/import.json"
    T2=$(python3 -c 'import time; print(time.perf_counter())')

    # Step 3: Solve prefix
    RC=0
    "${BINARY}" --json solve "${PROJECT}" "${LEDGER}" --materials "${PACK}" > "${ARTIFACT_DIR}/solve.json" 2> "${ARTIFACT_DIR}/solve.err" || RC=$?
    T3=$(python3 -c 'import time; print(time.perf_counter())')

    # Step 4: Report — the export verbs take the REAL run id the solve
    # retained; an export of an unknown run refuses and writes nothing.
    RUN_ID="$(grep -oE '"run":"[0-9a-f]{64}"' "${ARTIFACT_DIR}/solve.json" | head -1 | cut -d'"' -f4)"
    if [ "${#RUN_ID}" -ne 64 ]; then
      printf "ERROR: solve retained no run id (exit %s); see %s\n" "${RC}" "${ARTIFACT_DIR}/solve.err" >&2
      exit 1
    fi
    (cd "${RUN_DIR}" && "${BINARY}" --json report "${RUN_ID}" "${LEDGER}" > "${ARTIFACT_DIR}/report.json")
    T4=$(python3 -c 'import time; print(time.perf_counter())')

    # Step 5: Package
    (cd "${RUN_DIR}" && "${BINARY}" --json package "${RUN_ID}" "${LEDGER}" > "${ARTIFACT_DIR}/package.json")
    T5=$(python3 -c 'import time; print(time.perf_counter())')

    printf "==> 2. Synthesizing PipelineAttributionReceipt\n"
    python3 - <<EOF
import json
import platform

t_validate = max(0.001, $T1 - $T0)
t_import = max(0.001, $T2 - $T1)
t_solve = max(0.001, $T3 - $T2)
t_report = max(0.001, $T4 - $T3)
t_package = max(0.001, $T5 - $T4)
t_total = t_validate + t_import + t_solve + t_report + t_package

def to_bps(part, total):
    return int(round((part / total) * 10000))

phases = [
    {
        "phase": "validate",
        "wall_s": round(t_validate, 6),
        "wall_share_bps": to_bps(t_validate, t_total),
        "energy_j": None,
        "memory_bytes": None,
        "is_accelerator_addressable": False
    },
    {
        "phase": "import",
        "wall_s": round(t_import, 6),
        "wall_share_bps": to_bps(t_import, t_total),
        "energy_j": None,
        "memory_bytes": None,
        "is_accelerator_addressable": False
    },
    {
        "phase": "solve_prefix",
        "wall_s": round(t_solve, 6),
        "wall_share_bps": to_bps(t_solve, t_total),
        "energy_j": None,
        "memory_bytes": None,
        "is_accelerator_addressable": True
    },
    {
        "phase": "report",
        "wall_s": round(t_report, 6),
        "wall_share_bps": to_bps(t_report, t_total),
        "energy_j": None,
        "memory_bytes": None,
        "is_accelerator_addressable": False
    },
    {
        "phase": "package",
        "wall_s": round(t_package, 6),
        "wall_share_bps": to_bps(t_package, t_total),
        "energy_j": None,
        "memory_bytes": None,
        "is_accelerator_addressable": False
    }
]

# No profiler ran in this lane. Kernel attribution is therefore NOT measured
# and is not synthesized: an earlier revision wrote literal fractions of the
# solve wall-clock as kernel rows (0.45 / 0.35 / 0.10 with invented
# arithmetic intensities) under the "measured" authority string and let the
# accelerator falsifier admit a candidate from them. That was a fabricated
# receipt. Until a named profiler produces kernel rows, the list is empty,
# the falsifier refuses with evidence, and the receipt says why.
top_three_kernels = []
top_three_bps = 0
lead_bps = 0
top_three_meets = False
lead_meets = False
decision = "refused-with-evidence"
reason = (
    "no measured kernel attribution: this lane times the CLI phases only "
    "(validate, import, solve, report, package) and ran no profiler; the "
    "accelerator doctrine falsifier cannot evaluate top-three or lead-kernel "
    "shares from stage wall-clock alone"
)

receipt = {
    "schema": "frankensim.roofline.pipeline-attribution.v1",
    "workflow_id": "cooling_01_profile_campaign",
    "machine_fingerprint": f"{platform.system().lower()}-{platform.machine()}",
    "isa_family": platform.machine(),
    "total_wall_s": round(t_total, 6),
    "total_energy_j": None,
    "phases": phases,
    "top_three_kernels": top_three_kernels,
    "falsifier": {
        "top_three_wall_share_bps": top_three_bps,
        "top_three_meets_gate": top_three_meets,
        "selected_kernel_wall_share_bps": lead_bps,
        "selected_kernel_meets_gate": lead_meets,
        "decision": decision,
        "reason": reason
    },
    "authority": "measured-stage-wall-clock-only",
    "kernel_attribution": {"measured": False, "profiler": None, "reason": "no profiler ran; kernel rows are absent, not estimated"},
    "measured": {"phase_wall_s": True, "energy_j": False, "memory_bytes": False, "kernels": False},
    "no_claim": "this lane measures CLI phase wall-clock on one host; it attributes nothing to kernels, measures no energy or memory, and authorizes no device execution or speedup claim; the accelerator falsifier is refused by construction until a named profiler supplies measured kernel rows"
}

with open("${ARTIFACT_DIR}/pipeline_attribution_receipt.json", "w") as f:
    json.dump(receipt, f, indent=2)

print(f"Attribution summary: total wall {t_total:.4f}s across {len(phases)} phases")
print(f"Top 3 kernel share: {top_three_bps/100:.1f}% (lead: {lead_bps/100:.1f}%) -> Verdict: {decision}")
EOF

    printf "==> 3. Verifying receipt artifact integrity\n"
    grep -q '"schema": "frankensim.roofline.pipeline-attribution.v1"' "${ARTIFACT_DIR}/pipeline_attribution_receipt.json"
    grep -q '"authority": "measured-stage-wall-clock-only"' "${ARTIFACT_DIR}/pipeline_attribution_receipt.json"
    grep -q '"decision": "refused-with-evidence"' "${ARTIFACT_DIR}/pipeline_attribution_receipt.json"
    grep -q '"top_three_kernels": \[\]' "${ARTIFACT_DIR}/pipeline_attribution_receipt.json"

    printf "All vertical profiling checks passed!\n"
    exit 0
    ;;
  *)
    printf "FATAL: unknown command %s\n" "${COMMAND}" >&2
    exit 2
    ;;
esac
