#!/usr/bin/env bash
# ASCENT conformance profile runner.
#
# The PR profile is a deliberately small, family-representative selector. The
# nightly profile is the complete fs-ascent package with Cargo fail-fast
# disabled. Both profiles are lockfile-closed and enforce an aggregate wall
# budget that intentionally includes compilation. This is an execution guard,
# not a machine-independent performance claim.

set -euo pipefail

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
readonly REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
readonly CARGO_BIN="${CARGO_BIN:-cargo}"
readonly PR_BUDGET_SECONDS="${FS_ASCENT_PR_BUDGET_SECONDS:-900}"
readonly NIGHTLY_BUDGET_SECONDS="${FS_ASCENT_NIGHTLY_BUDGET_SECONDS:-7200}"

readonly -a PR_TARGETS=(
  ascent_battery
  constrained_battery
  pareto_battery
  runner_battery
  budget_trend_manifest
  bbob_budget_ledger
)

usage() {
  cat >&2 <<'EOF'
usage:
  scripts/ci/ascent_conformance_profile.sh pr
  scripts/ci/ascent_conformance_profile.sh nightly
  scripts/ci/ascent_conformance_profile.sh --list pr
  scripts/ci/ascent_conformance_profile.sh --list nightly

environment:
  CARGO_BIN                         Cargo executable path (default: cargo)
  FS_ASCENT_PR_BUDGET_SECONDS       Aggregate PR wall budget (default: 900)
  FS_ASCENT_NIGHTLY_BUDGET_SECONDS  Aggregate nightly wall budget (default: 7200)

The aggregate wall budget includes compilation and test execution. Callers may
override it for a declared host/target policy; the emitted receipt retains the
effective value.
EOF
}

require_positive_integer() {
  local label="$1"
  local value="$2"
  if [[ ! "${value}" =~ ^[1-9][0-9]*$ ]]; then
    printf 'invalid %s: expected a positive integer, got %q\n' "${label}" "${value}" >&2
    exit 2
  fi
}

profile_budget() {
  local profile="$1"
  case "${profile}" in
    pr) printf '%s\n' "${PR_BUDGET_SECONDS}" ;;
    nightly) printf '%s\n' "${NIGHTLY_BUDGET_SECONDS}" ;;
    *) return 2 ;;
  esac
}

list_profile() {
  local profile="$1"
  local budget
  budget="$(profile_budget "${profile}")" || {
    usage
    exit 2
  }
  printf '{"schema":"frankensim-ascent-conformance-profile-v1","profile":"%s","budget_seconds":%s,"build_time_included":true}\n' \
    "${profile}" "${budget}"
  if [[ "${profile}" == "pr" ]]; then
    local target
    for target in "${PR_TARGETS[@]}"; do
      printf '{"profile":"pr","package":"fs-ascent","target":"%s","selector":"cargo test --locked -p fs-ascent --test %s -- --nocapture"}\n' \
        "${target}" "${target}"
    done
  else
    printf '%s\n' '{"profile":"nightly","package":"fs-ascent","target":"all","selector":"cargo test --locked -p fs-ascent --no-fail-fast -- --nocapture"}'
  fi
}

emit_target_result() {
  local profile="$1"
  local target="$2"
  local status="$3"
  local elapsed_seconds="$4"
  printf '{"schema":"frankensim-ascent-conformance-target-v1","profile":"%s","package":"fs-ascent","target":"%s","status":"%s","elapsed_seconds":%s}\n' \
    "${profile}" "${target}" "${status}" "${elapsed_seconds}"
}

run_pr_target() {
  local target="$1"
  local started_at finished_at exit_code status
  started_at="$(date +%s)"
  set +e
  "${CARGO_BIN}" test --locked -p fs-ascent --test "${target}" -- --nocapture
  exit_code=$?
  set -e
  finished_at="$(date +%s)"
  status="pass"
  if (( exit_code != 0 )); then
    status="fail"
  fi
  emit_target_result "pr" "${target}" "${status}" "$((finished_at - started_at))"
  return "${exit_code}"
}

run_nightly() {
  local started_at finished_at exit_code status
  started_at="$(date +%s)"
  set +e
  "${CARGO_BIN}" test --locked -p fs-ascent --no-fail-fast -- --nocapture
  exit_code=$?
  set -e
  finished_at="$(date +%s)"
  status="pass"
  if (( exit_code != 0 )); then
    status="fail"
  fi
  emit_target_result "nightly" "all" "${status}" "$((finished_at - started_at))"
  return "${exit_code}"
}

if (( $# == 2 )) && [[ "$1" == "--list" ]]; then
  require_positive_integer "PR budget" "${PR_BUDGET_SECONDS}"
  require_positive_integer "nightly budget" "${NIGHTLY_BUDGET_SECONDS}"
  list_profile "$2"
  exit 0
fi

if (( $# != 1 )) || [[ "$1" != "pr" && "$1" != "nightly" ]]; then
  usage
  exit 2
fi

readonly PROFILE="$1"
readonly BUDGET_SECONDS="$(profile_budget "${PROFILE}")"
require_positive_integer "${PROFILE} budget" "${BUDGET_SECONDS}"

cd "${REPO_ROOT}"
readonly STARTED_AT="$(date +%s)"
printf '{"schema":"frankensim-ascent-conformance-run-v1","event":"start","profile":"%s","budget_seconds":%s,"build_time_included":true}\n' \
  "${PROFILE}" "${BUDGET_SECONDS}"

overall_status=0
if [[ "${PROFILE}" == "pr" ]]; then
  for target in "${PR_TARGETS[@]}"; do
    if ! run_pr_target "${target}"; then
      overall_status=1
    fi
  done
elif ! run_nightly; then
  overall_status=1
fi

readonly FINISHED_AT="$(date +%s)"
readonly ELAPSED_SECONDS="$((FINISHED_AT - STARTED_AT))"
budget_status="within"
if (( ELAPSED_SECONDS > BUDGET_SECONDS )); then
  budget_status="exceeded"
  overall_status=1
fi
run_status="pass"
if (( overall_status != 0 )); then
  run_status="fail"
fi
printf '{"schema":"frankensim-ascent-conformance-run-v1","event":"finish","profile":"%s","status":"%s","budget_status":"%s","budget_seconds":%s,"elapsed_seconds":%s,"build_time_included":true}\n' \
  "${PROFILE}" "${run_status}" "${budget_status}" "${BUDGET_SECONDS}" "${ELAPSED_SECONDS}"

exit "${overall_status}"
