#!/usr/bin/env bash
# IRP S1 structural subset, following material_sources.sh receipt conventions.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
readonly REPO_ROOT
readonly -a EXPECTED_TESTS=(
  g1_mms_orders_isotropic_and_rotated_orthotropic
  g1_curved_boundary_projection_reduces_geometry_error
  g2_cooks_membrane_plane_strain_equivalent
  g2_macneal_harder_straight_cantilever
  g1_free_free_bar_eigen_ladder
)
readonly REDUCE_TEST='reduce::reduce_tests::rl_005_tet10_cantilever_card'
readonly -a SHELL_TESTS=(
  g2_scordelis_lo_roof
  g2_pinched_cylinder_with_diaphragms
  g2_pinched_hemisphere_18deg_hole
  g2_pinched_hemisphere_full
  g2_obstacle_course_distorted_meshes
  g0_drilling_alpha_insensitivity
  report_dkt_rung_on_obstacle_course
)

list() {
  local test
  for test in "${EXPECTED_TESTS[@]}"; do printf 'fs-solid:tet10::%s\n' "$test"; done
  printf 'fs-solid:lib::%s\n' "$REDUCE_TEST"
  for test in "${SHELL_TESTS[@]}"; do printf 'fs-plate:shell_obstacle_course::%s\n' "$test"; done
}

check() {
  local test
  for test in "${EXPECTED_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-solid/tests/tet10.rs" || {
      printf 'missing test: %s\n' "$test" >&2; return 1;
    }
  done
  grep -Fq 'fn rl_005_tet10_cantilever_card(' "${REPO_ROOT}/crates/fs-solid/src/reduce/tests.rs"
  for test in "${SHELL_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-plate/tests/shell_obstacle_course.rs" || {
      printf 'missing shell test: %s\n' "$test" >&2; return 1;
    }
  done
  printf 'IRP S1: five solid G1/G2 cases, DEF-08 and seven shell cases present\n'
}

run() {
  check
  # This swarm uses central Linux verification. Never start a Cargo workload
  # on the disk-constrained Mac, including an accidental RCH local fallback.
  if [[ "$(uname -s)" == Darwin ]]; then
    printf 'IRP S1 --run belongs to the central Linux verification host\n' >&2
    return 2
  fi
  local output_dir started git_head host producer_exit=0 status=0 test missing=0
  output_dir="$(mktemp -d "${IRP_S1_RECEIPT_BASE:-${TMPDIR:-/tmp}}/frankensim-irp-s1.XXXXXX")"
  started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  git_head="$(git -C "$REPO_ROOT" rev-parse HEAD)"
  host="$(uname -n)"
  printf 'IRP S1 live receipt: %s\n' "$output_dir"
  local -a tet_command lib_command shell_command
  if [[ -n "${IRP_S1_TET10_BIN:-}" && -n "${IRP_S1_LIB_BIN:-}" && -n "${IRP_S1_SHELL_BIN:-}" ]]; then
    tet_command=("$IRP_S1_TET10_BIN" --exact --test-threads=1 --show-output "${EXPECTED_TESTS[@]}")
    lib_command=("$IRP_S1_LIB_BIN" --exact --test-threads=1 --show-output "$REDUCE_TEST")
    shell_command=("$IRP_S1_SHELL_BIN" --exact --test-threads=1 --show-output "${SHELL_TESTS[@]}")
  else
    tet_command=(cargo test --locked -j 12 -p fs-solid --test tet10 -- --exact --test-threads=1 --show-output "${EXPECTED_TESTS[@]}")
    lib_command=(cargo test --locked -j 12 -p fs-solid --lib -- --exact --test-threads=1 --show-output "$REDUCE_TEST")
    shell_command=(cargo test --locked -j 12 -p fs-plate --test shell_obstacle_course -- --exact --test-threads=1 --show-output "${SHELL_TESTS[@]}")
  fi
  {
    printf 'host=%s\ngit_head=%s\nstarted_utc=%s\n' "$host" "$git_head" "$started"
    printf '%q ' "${tet_command[@]}"; printf '\n'
    printf '%q ' "${lib_command[@]}"; printf '\n'
    printf '%q ' "${shell_command[@]}"; printf '\n'
  } >"${output_dir}/command.txt"
  export CARGO_BUILD_JOBS=12
  export CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test"
  set +e
  (cd "$REPO_ROOT" && "${tet_command[@]}") >"${output_dir}/stdout.log" 2>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  (cd "$REPO_ROOT" && "${lib_command[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  (cd "$REPO_ROOT" && "${shell_command[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  set -e
  for test in "${EXPECTED_TESTS[@]}" "$REDUCE_TEST" "${SHELL_TESTS[@]}"; do
    if ! grep -Fqx "test ${test} ... ok" "${output_dir}/stdout.log" "${output_dir}/stderr.log"; then
      printf 'missing or failed: %s\n' "$test" >&2; missing=1
    fi
  done
  if ! grep -Eq '^test result: ok\. [1-9][0-9]* passed' "${output_dir}/stdout.log" "${output_dir}/stderr.log"; then missing=1; fi
  python3 -I - "$output_dir" "$host" "$git_head" "$started" "$producer_exit" "$missing" <<'PY'
import datetime
import json
import pathlib
import sys
directory = pathlib.Path(sys.argv[1])
host, head, started = sys.argv[2:5]
producer_exit, missing = map(int, sys.argv[5:7])
command = (directory / "command.txt").read_text()
events = []
for name in ("stdout.log", "stderr.log"):
    for line in (directory / name).read_text(errors="replace").splitlines():
        try:
            row = json.loads(line)
        except ValueError:
            continue
        if not isinstance(row, dict) or row.get("suite") != "irp-s1":
            continue
        row.update(host=host, git_head=head, command=command,
                   utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        events.append(row)
with (directory / "events.jsonl").open("w") as output:
    for event in events:
        output.write(json.dumps(event, allow_nan=False) + "\n")
passed = producer_exit == 0 and missing == 0 and bool(events)
summary = dict(schema="frankensim.irp-s1-structures.e2e.v1", host=host,
               git_head=head, started_utc=started, producer_exit=producer_exit,
               missing_or_failed=missing, events=len(events),
               command_file="command.txt", stdout_log="stdout.log",
               stderr_log="stderr.log", events_log="events.jsonl",
               verdict="pass" if passed else "failed")
(directory / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
sys.exit(0 if passed else 1)
PY
}

case "${1:---check}" in
  --list) list ;;
  --check) check ;;
  --run) run ;;
  *) printf 'usage: %s [--list|--check|--run]\n' "$0" >&2; exit 2 ;;
esac
