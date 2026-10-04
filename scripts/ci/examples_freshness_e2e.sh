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
#   6. heatsink-fan-rotated     — the same project on the shell rotated
#                                35/21 deg (every facet oblique, full-precision
#                                ASCII import): imports, solves all seven
#                                stages, and its T_max matches the
#                                axis-aligned body's within 1 mK.
#   7. heatsink-fan-chip        — the same heatsink with its 3 W entering
#                                through a declared 20 x 20 mm surface
#                                (fsim v8): seven stages, and hotter than
#                                the volumetric source.
#   8. examples/plate-hole      — a genus-1 body (square through-hole),
#                                no air network: seven stages, and the
#                                exact energy-balance bracket
#                                T_min <= T_ref + P/(hA) <= T_max with the
#                                hole walls counted in A.
#   9. examples/contact-pair    — two sources in one import, a declared
#                                card-backed contact joint: seven stages
#                                (the exact -5 W / -0.5 K oracles are the
#                                G1 test's).
#  10. examples/perforated-plate — the mesher complexity gate: an 8652-facet
#                                OBLIQUE plate generated in the work dir,
#                                imported, solved (7 stages) in < 120 s.
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
    -h|--help) sed -n '3,44p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
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
check "report exports the published field as VTU (mesh + nodal temperature)" \
  grep -q '<DataArray[^>]*Name="temperature"' "${WORK}/${RUN_ID}.field.vtu"
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
# All eight terms measured: the ladder's discretization half-width, the
# declared-input propagation's boundary-condition, model-form, solver,
# parameters (declared AA6061 conductivity tolerance) and geometry (declared
# +/-0.05 mm surface offset) terms, the componentwise adjoint roundoff bound,
# and a negligible measurement term. The conservative sum (about 21.4 K)
# clears the 5 K margin, so the verdict is an Estimated `satisfied`.
# The solver term is a verified coupled enclosure; on the finest rung the solid
# inverse comes from bounded interval LDL elimination (driver 30), inside the
# ladder project's declared 256 MiB. Missing evidence would be NO-DATA.
check "QoI stage measured all eight budget terms" grep -q '"budget_terms_measured":8' "${WORK}/ls.err"
check "no budget term remains NO-DATA" grep -q '"weakest_term":"none-no-data"' "${WORK}/ls.err"
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
check "ladder verdict is the Estimated satisfied decision of the complete budget" grep -q '"verdict":"satisfied"' "${WORK}/lrep.json"

# ---- 9. the rotated twin: rotation invariance of the whole product path -----
# The same shell rotated 35 deg about z and 21 deg about x and translated
# (generate_heatsink_stl.py --rotate 35 21 --shift 0.1 0.2 0.05): every facet
# oblique, imported at the file's full ASCII precision. It must import, mesh
# (fs-mesh CONTRACT items 18-23), solve every stage, and reproduce the
# axis-aligned body's maximum temperature on a different mesh within the
# discretization scale. TOLERANCE 1 mK: MEASURED 2026-09-30 after the f64
# ASCII import, axis 301.99571 K vs rotated 301.99605 K, |dT_max| = 0.34 mK
# (2026-09-03 under the f32 import: 713 vs 704 tets, 0.157 mK) on a rise of
# 8.85 K, ~6x headroom, and the ladder's data-range discretization bound on
# this body is 1.1 mK — two meshes of one body agreeing better than that
# bound is the honest expectation.
ROTATED="${REPO_ROOT}/examples/heatsink-fan/heatsink-fan-rotated.fsim"
ROTATED_STL="${REPO_ROOT}/examples/heatsink-fan/heatsink-rotated.stl"
check "heatsink-fan-rotated validates ok" validate_ok "${ROTATED}"
rotated_import_ok() {
  "${BINARY}" --json import "${ROTATED}" "${ROTATED_STL}" "${WORK}/rotated.db" --unit m --max-hole-edges 0 > "${WORK}/rimp.json" 2> "${WORK}/rimp.err"
}
check "rotated twin imports into its own ledger" rotated_import_ok
rotated_solve_completes() {
  "${BINARY}" --json solve "${ROTATED}" "${WORK}/rotated.db" --materials "${PACK}" > "${WORK}/rs.json" 2> "${WORK}/rs.err"
}
check "rotated twin solves every stage (oblique facets, full-precision import)" rotated_solve_completes
check "rotated solve reports seven completed stages" grep -q '"stages_completed":7' "${WORK}/rs.json"
ROTATED_RUN="$(grep -oE '"run":"[0-9a-f]{64}"' "${WORK}/rs.json" | head -1 | cut -d'"' -f4)"
rotated_report_ok() {
  (cd "${WORK}" && "${BINARY}" --json report "${ROTATED_RUN}" "${WORK}/rotated.db" > "${WORK}/rrep.json" 2> "${WORK}/rrep.err")
}
check "rotated report exports" rotated_report_ok
t_max_of() {
  grep -oE '"temperature_max": ?"[0-9.]+"' "$1" | head -1 | grep -oE '[0-9]+\.[0-9]+'
}
AXIS_TMAX="$(t_max_of "${WORK}/${RUN_ID}.report.json")"
ROTATED_TMAX="$(t_max_of "${WORK}/${ROTATED_RUN}.report.json")"
rotation_invariant() {
  test -n "${AXIS_TMAX}" && test -n "${ROTATED_TMAX}" \
    && awk -v a="${AXIS_TMAX}" -v r="${ROTATED_TMAX}" 'BEGIN { d = a - r; if (d < 0) d = -d; exit !(d <= 0.001) }'
}
check "rotated T_max matches the axis-aligned body within 1 mK (axis ${AXIS_TMAX:-?} K, rotated ${ROTATED_TMAX:-?} K)" rotation_invariant

# ---- 10. chip footprint: heat entering through a declared surface -----------
# The heatsink with its 3 W entering through a 20 x 20 mm die contact (an
# fsim v8 `(surface ...)` entity lowered to an inward Neumann flux). The
# concentrated source must run hotter than the volumetric one: MEASURED
# 2026-09-29 302.693 K vs 301.996 K.
CHIP="${REPO_ROOT}/examples/heatsink-fan/heatsink-fan-chip.fsim"
CHIP_STL="${REPO_ROOT}/examples/heatsink-fan/heatsink-chip.stl"
check "heatsink-fan-chip validates ok" validate_ok "${CHIP}"
chip_solve_completes() {
  "${BINARY}" --json import "${CHIP}" "${CHIP_STL}" "${WORK}/chip.db" --unit m --max-hole-edges 0 > "${WORK}/cimp.json" 2> "${WORK}/cimp.err" \
    && "${BINARY}" --json solve "${CHIP}" "${WORK}/chip.db" --materials "${PACK}" > "${WORK}/cs.json" 2> "${WORK}/cs.err"
}
check "chip-footprint heatsink imports and solves" chip_solve_completes
check "chip solve reports seven completed stages" grep -q '"stages_completed":7' "${WORK}/cs.json"
CHIP_RUN="$(grep -oE '"run":"[0-9a-f]{64}"' "${WORK}/cs.json" | head -1 | cut -d'"' -f4)"
chip_report_ok() {
  (cd "${WORK}" && "${BINARY}" --json report "${CHIP_RUN}" "${WORK}/chip.db" > "${WORK}/crep.json" 2> "${WORK}/crep.err")
}
check "chip report exports" chip_report_ok
CHIP_TMAX="$(t_max_of "${WORK}/${CHIP_RUN}.report.json")"
chip_hotter() {
  test -n "${AXIS_TMAX}" && test -n "${CHIP_TMAX}" \
    && awk -v a="${AXIS_TMAX}" -v c="${CHIP_TMAX}" 'BEGIN { exit !(c > a) }'
}
check "chip footprint runs hotter than the volumetric source (volumetric ${AXIS_TMAX:-?} K, chip ${CHIP_TMAX:-?} K)" chip_hotter

# ---- 11. plate with a through-hole: a genus-1 body, no air network ---------
# generate_plate_hole_stl.py asserts the analytic volume (9.024 cm^3), area
# (5504 mm^2, hole walls included) and Euler characteristic 0. With no fan,
# vent or leakage the flow-network stage retains a not-declared receipt. At
# steady state all 2 W leave by h = 10 W/m^2/K convection, so the
# area-weighted mean surface temperature is EXACTLY 293.15 + 2/(10 A) =
# 329.4872 K and must lie between T_min and T_max. A body whose hole was lost
# (A = 5600 mm^2) has mean 328.864 K and a T_max below the bracket. MEASURED
# 2026-09-29: T_min 329.4643, T_max 329.4992 K on 180 tets.
PLATE="${REPO_ROOT}/examples/plate-hole/plate-hole.fsim"
PLATE_STL="${REPO_ROOT}/examples/plate-hole/plate-hole.stl"
check "plate-hole validates ok" validate_ok "${PLATE}"
plate_solve_completes() {
  "${BINARY}" --json import "${PLATE}" "${PLATE_STL}" "${WORK}/plate.db" --unit m --max-hole-edges 0 > "${WORK}/pimp.json" 2> "${WORK}/pimp.err" \
    && "${BINARY}" --json solve "${PLATE}" "${WORK}/plate.db" --materials "${PACK}" > "${WORK}/ps.json" 2> "${WORK}/ps.err"
}
check "plate-hole imports and solves" plate_solve_completes
check "plate-hole solve reports seven completed stages" grep -q '"stages_completed":7' "${WORK}/ps.json"
PLATE_RUN="$(grep -oE '"run":"[0-9a-f]{64}"' "${WORK}/ps.json" | head -1 | cut -d'"' -f4)"
plate_report_ok() {
  (cd "${WORK}" && "${BINARY}" --json report "${PLATE_RUN}" "${WORK}/plate.db" > "${WORK}/prep.json" 2> "${WORK}/prep.err")
}
check "plate-hole report exports" plate_report_ok
check "plate-hole flow stage retains the declared absence of an air network" \
  grep -q 'declared-absence-of-an-air-network' "${WORK}/${PLATE_RUN}.report.json"
PLATE_TMAX="$(t_max_of "${WORK}/${PLATE_RUN}.report.json")"
PLATE_TMIN="$(grep -oE '"temperature_min": ?"[0-9.]+"' "${WORK}/${PLATE_RUN}.report.json" | head -1 | grep -oE '[0-9]+\.[0-9]+')"
plate_bracket() {
  test -n "${PLATE_TMIN}" && test -n "${PLATE_TMAX}" \
    && awk -v lo="${PLATE_TMIN}" -v hi="${PLATE_TMAX}" 'BEGIN { m = 293.15 + 2.0 / (10.0 * 5.504e-3); exit !(lo <= m && m <= hi) }'
}
check "plate-hole brackets the exact mean surface temperature 329.4872 K (T_min ${PLATE_TMIN:-?}, T_max ${PLATE_TMAX:-?})" plate_bracket

# ---- 12. contact pair: two sources, one declared contact joint -------------
PAIR_DIR="${REPO_ROOT}/examples/contact-pair"
PAIR="${PAIR_DIR}/contact-pair.fsim"
check "contact-pair validates ok" validate_ok "${PAIR}"
pair_solve_completes() {
  "${BINARY}" --json import "${PAIR}" "${PAIR_DIR}/cold-body.stl" "${PAIR_DIR}/hot-body.stl" "${WORK}/pair.db" --unit m --max-hole-edges 0 > "${WORK}/cpimp.json" 2> "${WORK}/cpimp.err" \
    && "${BINARY}" --json solve "${PAIR}" "${WORK}/pair.db" --materials "${PACK}" --interfaces "${PAIR_DIR}/cold-hot.fsintpk" > "${WORK}/cps.json" 2> "${WORK}/cps.err"
}
check "contact pair imports both sources and solves" pair_solve_completes
check "contact pair import retained two artifacts" grep -q '"artifact_count":2' "${WORK}/cpimp.json"
check "contact pair solve reports seven completed stages" grep -q '"stages_completed":7' "${WORK}/cps.json"
pair_swapped_refuses() {
  ! "${BINARY}" --json import "${PAIR}" "${PAIR_DIR}/hot-body.stl" "${PAIR_DIR}/cold-body.stl" "${WORK}/swapped.db" --unit m --max-hole-edges 0 > /dev/null 2> "${WORK}/cpswap.err" \
    && grep -q 'cli-import-source-hash-mismatch' "${WORK}/cpswap.err"
}
check "contact pair refuses swapped sources on the pinned hash" pair_swapped_refuses

# ---- 13. complexity gate: 8652-facet oblique perforated plate ---------------
# generate_perforated_stl.py 20 --rotate is deterministic, so the project
# pins its source hash and the lane generates the ~3 MB STL instead of
# tracking it. TOLERANCE 120 s wall for the solve: MEASURED 2026-09-30 23.8 s
# conduction (~5x headroom on a shared host). Before that day's mesher fixes
# this body refused outright (f32-import slivers), and a 5644-facet
# axis-aligned plate took 199 s.
PERF_DIR="${REPO_ROOT}/examples/perforated-plate"
PERF="${PERF_DIR}/perforated-plate-rotated.fsim"
check "perforated-plate-rotated validates ok" validate_ok "${PERF}"
perf_generate() {
  python3 "${PERF_DIR}/generate_perforated_stl.py" 20 "${WORK}/perforated.stl" --rotate > /dev/null
}
check "perforated plate generator writes the pinned 8652-facet body" perf_generate
perf_solve_in_budget() {
  "${BINARY}" --json import "${PERF}" "${WORK}/perforated.stl" "${WORK}/perf.db" --unit m --max-hole-edges 0 > "${WORK}/pfimp.json" 2> "${WORK}/pfimp.err" || return 1
  local start end
  start=$(date +%s)
  "${BINARY}" --json solve "${PERF}" "${WORK}/perf.db" --materials "${PACK}" > "${WORK}/pfs.json" 2> "${WORK}/pfs.err" || return 1
  end=$(date +%s)
  log perf "{\"facets\":8652,\"solve_wall_s\":$((end - start))}"
  [[ $((end - start)) -lt 120 ]]
}
check "8652-facet oblique plate imports and solves within 120 s" perf_solve_in_budget
check "oblique plate solve reports seven completed stages" grep -q '"stages_completed":7' "${WORK}/pfs.json"

# ---- 14. passive heatsink: natural convection by fixed point (fsim v10) ----
# Same body and 3 W, no fan: the Churchill-Chu card's coefficient is
# iterated against the solved wall temperature. MEASURED 2026-10-01: 16
# iterations, h 5.90 W/m2K, T_max 316.74 K. The passive part must run hotter
# than the fan-cooled one, and the receipt must record the converged law.
NATURAL="${REPO_ROOT}/examples/heatsink-fan/heatsink-natural.fsim"
check "heatsink-natural validates ok" validate_ok "${NATURAL}"
natural_solve_completes() {
  "${BINARY}" --json import "${NATURAL}" "${REPO_ROOT}/examples/heatsink-fan/heatsink.stl" "${WORK}/natural.db" --unit m --max-hole-edges 0 > "${WORK}/nimp.json" 2> "${WORK}/nimp.err" \
    && "${BINARY}" --json solve "${NATURAL}" "${WORK}/natural.db" --materials "${PACK}" > "${WORK}/ns.json" 2> "${WORK}/ns.err"
}
check "passive heatsink imports and solves" natural_solve_completes
check "passive heatsink solve reports seven completed stages" grep -q '"stages_completed":7' "${WORK}/ns.json"
NATURAL_RUN="$(grep -oE '"run":"[0-9a-f]{64}"' "${WORK}/ns.json" | head -1 | cut -d'"' -f4)"
natural_report_ok() {
  (cd "${WORK}" && "${BINARY}" --json report "${NATURAL_RUN}" "${WORK}/natural.db" > "${WORK}/nrep.json" 2> "${WORK}/nrep.err")
}
check "passive heatsink report exports" natural_report_ok
NATURAL_TMAX="$(t_max_of "${WORK}/${NATURAL_RUN}.report.json")"
natural_hotter() {
  test -n "${AXIS_TMAX}" && test -n "${NATURAL_TMAX}" \
    && awk -v a="${AXIS_TMAX}" -v n="${NATURAL_TMAX}" 'BEGIN { exit !(n > a + 5.0) }'
}
check "the passive heatsink runs more than 5 K hotter than the fan-cooled one (fan ${AXIS_TMAX:-?} K, passive ${NATURAL_TMAX:-?} K)" natural_hotter

# ------------------------------------------------------------------- verdict
log summary "checks=${CHECKS} failures=${FAILURES}"
if [[ "${FAILURES}" -gt 0 ]]; then
  printf 'FAILED: %d of %d freshness checks failed; full NDJSON log: %s\n' \
    "${FAILURES}" "${CHECKS}" "${LOG}" >&2
  exit 1
fi
printf 'OK: all %d examples-freshness checks passed\n' "${CHECKS}"
