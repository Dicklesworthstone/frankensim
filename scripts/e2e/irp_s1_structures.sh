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
# M1.1-W1a (modal ZOH core moved to fs-modal::time): the cap test and the bitwise witnesses it must not move.
readonly MODAL_TIME_TEST='g0_mode_cap_4096_admitted_4097_refused'
# M1.3-W2a linear orthotropic viscoelastic TET10 elements: the G1 Prony
# ladder and complex-stiffness checks, and the G3 rotation metamorphic.
readonly -a VISCO_TESTS=(
  g1_prony_ladder_modal_loss_certified
  g1_time_domain_element_force_matches_complex_stiffness
  g3_rotation_metamorphic
)
readonly -a ND_TESTS=(
  g0_nd_is_a_permutation_and_separators_separate
  g0_nd_deterministic_repeat_bitwise
  g0_nd_leaf_size_bounds
  g0_factor_budget_cap_and_cap_plus_one
  g1_nd_inertia_matches_analytic_3d_grid
)
readonly -a WINDOW_TESTS=(
  g1_parallel_windows_certified_and_thread_invariant
  g1_one_symbolic_analysis_serves_all_shifts
  g0_window_count_cap_and_cap_plus_one
  g4_cancellation_drains_and_returns_cancelled
)
readonly -a HEX_TESTS=(
  g1_roller_block_exact_p_modes_order
  g2_slender_square_bar_free_free_vs_timoshenko
  g3_hex27_vs_tet10_same_bar
)
readonly -a OUTLINE_TESTS=(
  g0_cdt_conforms_to_outline_and_holes
  g0_quality_bound
  g1_curved_outline_chord_ladder
)
readonly -a ARCHING_TESTS=(
  g0_arching_surface_from_parameters
  g1_plate_with_hole_modes_converge
  g1_annular_tensioned_film_matches_bessel_roots
)
# M1.3-W4: generated and imported solid parts reach audited, curved TET10 meshes.
readonly -a SOLID_PARTS_TESTS=(
  g1_undercut_bar_with_nodal_holes_from_parameters
  g1_bridge_outline_with_holes_extruded
  g0_profile_generator_is_deterministic_and_closed
  g0_refuses_self_intersecting_outline_and_missing_reference_state
  g0_policy_cap_and_cap_plus_one
)
readonly -a TET_IMPROVE_TESTS=(
  g0_flip23_valid_and_nonconvex_refused
  g0_flip32_removes_a_degree_three_edge
  g0_inverted_interior_vertex_is_untangled_and_boundary_frozen
  g0_refusals
)
readonly -a GMSH_TESTS=(
  g0_tet10_hex20_hex27_round_trip
  g0_midside_nodes_at_midpoints_pin_ordering
  g0_unknown_type_refuses
  g0_element_cap_and_cap_plus_one
  g1_imported_mesh_promotes_and_assembles
  g0_promotion_refuses_undeclared_and_unsupported
)
readonly STEP_TEST='g1_step_tessellation_reaches_volume_mesh'
# M1.2-W1c: MITC3+ and DKT against a TET10 reference on real spruce.
readonly -a CROSSCHECK_TESTS=(
  g0_symmetry_class_merge_recovers_full_model_spectrum
  g0_mac_pairing_refuses_ambiguous_pairs
  g1_3d_reference_ladder_converges
  g2_shear_shell_matches_3d_reference_2_to_8_khz
  g2_dkt_deviation_table_and_refusal_band
  g1_isotropic_thin_plate_three_way_agreement
)
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
  for test in "${VISCO_TESTS[@]}"; do printf 'fs-solid:visco_elements::%s\n' "$test"; done
  printf 'fs-solid:lib::%s\n' "$REDUCE_TEST"
  for test in "${HEX_TESTS[@]}"; do printf 'fs-solid:hex::%s\n' "$test"; done
  for test in "${SHELL_TESTS[@]}"; do printf 'fs-plate:shell_obstacle_course::%s\n' "$test"; done
  for test in "${OUTLINE_TESTS[@]}"; do printf 'fs-mesh:planar_outline::%s\n' "$test"; done
  for test in "${SOLID_PARTS_TESTS[@]}"; do printf 'fs-mesh:solid_parts::%s\n' "$test"; done
  for test in "${TET_IMPROVE_TESTS[@]}"; do printf 'fs-mesh:tet_improve::%s\n' "$test"; done
  for test in "${GMSH_TESTS[@]}"; do printf 'fs-io:gmsh_high_order::%s\n' "$test"; done
  printf 'fs-io:step_import::%s\n' "$STEP_TEST"
  for test in "${CROSSCHECK_TESTS[@]}"; do printf 'fs-couple:shell_solid_crosscheck::%s\n' "$test"; done
  for test in "${ARCHING_TESTS[@]}"; do printf 'fs-plate:arched_outline_shell::%s\n' "$test"; done
  for test in "${ND_TESTS[@]}" g1_nd_fill_beats_amd_on_3d_grids; do printf 'fs-sparse:nested_dissection::%s\n' "$test"; done
  for test in "${WINDOW_TESTS[@]}" g2_million_dof_window_timing_receipt report_elastic_p1_tet10_fill; do printf 'fs-modal:parallel_windows::%s\n' "$test"; done
  printf 'fs-modal:lib::time::tests::%s\n' "$MODAL_TIME_TEST"
  printf 'fs-couple:bakeoff_string::string_bakeoff_receipt_is_reproduced_and_golden_pinned\n'
  printf 'fs-couple:bakeoff_dispersion::committed_dispersion_receipt_decides_and_modal_rows_are_untouched\n'
  printf 'engine-piano:clone::a_copy_renders_what_the_original_would\n'
  printf 'engine-piano:audio_rate::factor_one_retains_the_original_pressure_and_float_bits\n'
}

check() {
  local test
  for test in "${ND_TESTS[@]}" g1_nd_fill_beats_amd_on_3d_grids; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-sparse/tests/nested_dissection.rs" || return 1
  done
  for test in "${WINDOW_TESTS[@]}" g2_million_dof_window_timing_receipt report_elastic_p1_tet10_fill; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-modal/tests/parallel_windows.rs" || return 1
  done
  for test in "${VISCO_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-solid/tests/visco_elements.rs" || {
      printf 'missing test: %s\n' "$test" >&2; return 1;
    }
  done
  for test in "${EXPECTED_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-solid/tests/tet10.rs" || {
      printf 'missing test: %s\n' "$test" >&2; return 1;
    }
  done
  grep -Fq 'fn rl_005_tet10_cantilever_card(' "${REPO_ROOT}/crates/fs-solid/src/reduce/tests.rs"
  for test in "${HEX_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-solid/tests/hex.rs" || return 1
  done
  for test in "${SHELL_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-plate/tests/shell_obstacle_course.rs" || {
      printf 'missing shell test: %s\n' "$test" >&2; return 1;
    }
  done
  for test in "${OUTLINE_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-mesh/tests/planar_outline.rs" || return 1
  done
  for test in "${ARCHING_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-plate/tests/arched_outline_shell.rs" || return 1
  done
  for test in "${SOLID_PARTS_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-mesh/tests/solid_parts.rs" || { printf 'missing test: %s\n' "$test" >&2; return 1; }
  done
  for test in "${TET_IMPROVE_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-mesh/tests/tet_improve.rs" || { printf 'missing test: %s\n' "$test" >&2; return 1; }
  done
  for test in "${GMSH_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-io/tests/gmsh_high_order.rs" || { printf 'missing test: %s\n' "$test" >&2; return 1; }
  done
  grep -Fq "fn ${STEP_TEST}(" "${REPO_ROOT}/crates/fs-io/tests/step_import.rs" || { printf 'missing test: %s\n' "$STEP_TEST" >&2; return 1; }
  for test in "${CROSSCHECK_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-couple/tests/shell_solid_crosscheck.rs" || { printf 'missing test: %s\n' "$test" >&2; return 1; }
  done
  grep -Fq "fn ${MODAL_TIME_TEST}(" "${REPO_ROOT}/crates/fs-modal/src/time.rs" || { printf 'missing test: %s\n' "$MODAL_TIME_TEST" >&2; return 1; }
  grep -Fq 'fn string_bakeoff_receipt_is_reproduced_and_golden_pinned(' "${REPO_ROOT}/crates/fs-couple/tests/bakeoff_string.rs" || return 1
  grep -Fq 'fn committed_dispersion_receipt_decides_and_modal_rows_are_untouched(' "${REPO_ROOT}/crates/fs-couple/tests/bakeoff_dispersion.rs" || return 1
  grep -Fq 'fn a_copy_renders_what_the_original_would(' "${REPO_ROOT}/apps/instrument-lab/engine-piano/tests/clone.rs" || return 1
  grep -Fq 'fn factor_one_retains_the_original_pressure_and_float_bits(' "${REPO_ROOT}/apps/instrument-lab/engine-piano/tests/audio_rate.rs" || return 1
  printf 'IRP S1: solid, shell, outline, arching and solid-part cases present\n'
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
  local -a tet_command lib_command shell_command outline_command arching_command hex_command
  local -a visco_command
  if [[ -n "${IRP_S1_TET10_BIN:-}" && -n "${IRP_S1_LIB_BIN:-}" && -n "${IRP_S1_SHELL_BIN:-}" ]]; then
    tet_command=("$IRP_S1_TET10_BIN" --exact --test-threads=1 --show-output "${EXPECTED_TESTS[@]}")
    lib_command=("$IRP_S1_LIB_BIN" --exact --test-threads=1 --show-output "$REDUCE_TEST")
    shell_command=("$IRP_S1_SHELL_BIN" --exact --test-threads=1 --show-output "${SHELL_TESTS[@]}")
  else
    tet_command=(cargo test --locked -j 12 -p fs-solid --test tet10 -- --exact --test-threads=1 --show-output "${EXPECTED_TESTS[@]}")
    lib_command=(cargo test --locked -j 12 -p fs-solid --lib -- --exact --test-threads=1 --show-output "$REDUCE_TEST")
    shell_command=(cargo test --locked -j 12 -p fs-plate --test shell_obstacle_course -- --exact --test-threads=1 --show-output "${SHELL_TESTS[@]}")
  fi
  if [[ -n "${IRP_S1_OUTLINE_BIN:-}" ]]; then
    outline_command=("$IRP_S1_OUTLINE_BIN" --exact --test-threads=1 --show-output "${OUTLINE_TESTS[@]}")
  else
    outline_command=(cargo test --locked -j 12 -p fs-mesh --test planar_outline -- --exact --test-threads=1 --show-output "${OUTLINE_TESTS[@]}")
  fi
  if [[ -n "${IRP_S1_ARCHING_BIN:-}" ]]; then
    arching_command=("$IRP_S1_ARCHING_BIN" --exact --test-threads=1 --show-output "${ARCHING_TESTS[@]}")
  else
    arching_command=(cargo test --locked -j 12 -p fs-plate --test arched_outline_shell -- --exact --test-threads=1 --show-output "${ARCHING_TESTS[@]}")
  fi
  if [[ -n "${IRP_S1_HEX_BIN:-}" ]]; then
    hex_command=("$IRP_S1_HEX_BIN" --exact --test-threads=1 --show-output "${HEX_TESTS[@]}")
  else
    hex_command=(cargo test --locked -j 12 -p fs-solid --test hex -- --exact --test-threads=1 --show-output "${HEX_TESTS[@]}")
  fi
  if [[ -n "${IRP_S1_VISCO_BIN:-}" ]]; then
    visco_command=("$IRP_S1_VISCO_BIN" --exact --test-threads=1 --show-output "${VISCO_TESTS[@]}")
  else
    visco_command=(cargo test --release --locked -j 12 -p fs-solid --test visco_elements -- --exact --test-threads=1 --show-output "${VISCO_TESTS[@]}")
  fi
  {
    printf 'host=%s\ngit_head=%s\nstarted_utc=%s\n' "$host" "$git_head" "$started"
    printf '%q ' "${tet_command[@]}"; printf '\n'
    printf '%q ' "${lib_command[@]}"; printf '\n'
    printf '%q ' "${shell_command[@]}"; printf '\n'
    printf '%q ' "${outline_command[@]}"; printf '\n'
    printf '%q ' "${arching_command[@]}"; printf '\n'
    printf '%q ' "${hex_command[@]}"; printf '\n'
    printf '%q ' "${visco_command[@]}"; printf '\n'
    printf '%q ' cargo test --locked -j 12 -p fs-sparse --test nested_dissection -- --exact --test-threads=1 --show-output "${ND_TESTS[@]}"; printf '\n'
    printf '%q ' cargo test --locked -j 12 -p fs-modal --test parallel_windows -- --exact --test-threads=1 --show-output "${WINDOW_TESTS[@]}"; printf '\n'
    printf '%q ' cargo test --release --locked -j 12 -p fs-mesh --test solid_parts --test tet_improve -- --exact --test-threads=1 --show-output "${SOLID_PARTS_TESTS[@]}" "${TET_IMPROVE_TESTS[@]}"; printf '\n'
    printf '%q ' cargo test --locked -j 12 -p fs-io --test gmsh_high_order --test step_import -- --exact --test-threads=1 --show-output "${GMSH_TESTS[@]}" "$STEP_TEST"; printf '\n'
    printf '%q ' cargo test --release --locked -j 12 -p fs-couple --test shell_solid_crosscheck -- --exact --test-threads=1 --show-output "${CROSSCHECK_TESTS[@]}"; printf '\n'
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
  (cd "$REPO_ROOT" && "${outline_command[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  (cd "$REPO_ROOT" && "${arching_command[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  (cd "$REPO_ROOT" && "${hex_command[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  (cd "$REPO_ROOT" && "${visco_command[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  set -e
  set +e
  (cd "$REPO_ROOT" && cargo test --locked -j 12 -p fs-sparse --test nested_dissection -- --exact --test-threads=1 --show-output "${ND_TESTS[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  (cd "$REPO_ROOT" && cargo test --locked -j 12 -p fs-modal --test parallel_windows -- --exact --test-threads=1 --show-output "${WINDOW_TESTS[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  # The undercut-bar ladder reaches 1.7 M tets on its third rung: release only.
  (cd "$REPO_ROOT" && cargo test --release --locked -j 12 -p fs-mesh --test solid_parts --test tet_improve -- --exact --test-threads=1 --show-output "${SOLID_PARTS_TESTS[@]}" "${TET_IMPROVE_TESTS[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  (cd "$REPO_ROOT" && cargo test --locked -j 12 -p fs-io --test gmsh_high_order --test step_import -- --exact --test-threads=1 --show-output "${GMSH_TESTS[@]}" "$STEP_TEST") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  # The spruce TET10 ladder reaches ~7.9e4 DOF per mirror class: release only.
  (cd "$REPO_ROOT" && cargo test --release --locked -j 12 -p fs-couple --test shell_solid_crosscheck -- --exact --test-threads=1 --show-output "${CROSSCHECK_TESTS[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  set -e
  for test in "${ND_TESTS[@]}" "${WINDOW_TESTS[@]}" "${SOLID_PARTS_TESTS[@]}" "${TET_IMPROVE_TESTS[@]}" "${GMSH_TESTS[@]}" "$STEP_TEST" "${CROSSCHECK_TESTS[@]}"; do
    if ! grep -Fqx "test ${test} ... ok" "${output_dir}/stdout.log" "${output_dir}/stderr.log"; then missing=1; fi
  done
  for test in "${EXPECTED_TESTS[@]}" "${VISCO_TESTS[@]}" "$REDUCE_TEST" "${SHELL_TESTS[@]}" "${OUTLINE_TESTS[@]}" "${ARCHING_TESTS[@]}" "${HEX_TESTS[@]}"; do
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

# M1.1-W1a lift-parity stage. Runs each declared render command from binaries built at ONE commit and prints one JSON line per
# fixture. Run it at the parent export and at the change commit on the same host, then compare with --compare.
#   IRP_S1_MUSIC_RENDER_BIN  music_render        IRP_S1_GRAND_PIANO_BIN  grand_piano example
#   IRP_S1_PIANO_RENDER_BIN  engine-piano render example
renders() {
  local out="${1:?usage: $0 --renders OUT.jsonl}"
  [[ "$(uname -s)" != Darwin ]] || { echo 'IRP S1 --renders belongs to the central Linux verification host' >&2; return 2; }
  local music="${IRP_S1_MUSIC_RENDER_BIN:?set IRP_S1_MUSIC_RENDER_BIN}" piano="${IRP_S1_GRAND_PIANO_BIN:?set IRP_S1_GRAND_PIANO_BIN}"
  local live="${IRP_S1_PIANO_RENDER_BIN:?set IRP_S1_PIANO_RENDER_BIN}" work head host
  work="$(mktemp -d "${IRP_S1_RECEIPT_BASE:-${TMPDIR:-/tmp}}/frankensim-irp-s1-renders.XXXXXX")"
  head="$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
  host="$(uname -n)"
  : >"$out"
  printf '0,note_on,60,2.5\n192000,note_off,60,0\n' >"$work/c4_held.csv"
  emit() {
    local fixture="$1" wav="$2" started="$3"; shift 3
    python3 -I - "$fixture" "$wav" "$started" "$head" "$host" "$*" >>"$out" <<'PY'
import hashlib, json, os, sys, time
fixture, wav, started, head, host, command = sys.argv[1:7]
data = open(wav, 'rb').read() if os.path.exists(wav) else b''
row = dict(suite='irp-s1', bead='M1.1-W1a', fixture=fixture, command=command, wav_sha256=hashlib.sha256(data).hexdigest(),
           bytes=len(data), seconds=round(time.time() - float(started), 3), host=host, git_head=head)
side = os.path.splitext(wav)[0] + '.provenance.json'
if os.path.exists(side):
    row['wav_blake3'] = json.load(open(side)).get('wav_blake3')
print(json.dumps(row, allow_nan=False))
PY
  }
  local t
  t="$(date +%s.%N)"; "$music" string "$work/string.wav" --seconds 0.25 >/dev/null; emit music-string "$work/string.wav" "$t" "$music" string --seconds 0.25
  t="$(date +%s.%N)"; "$piano" --note 69 --duration 0.1 --modes 12 --render "$work/piano-smoke.wav" >/dev/null; emit piano-smoke "$work/piano-smoke.wav" "$t" "$piano" --note 69 --duration 0.1 --modes 12
  t="$(date +%s.%N)"; "$live" "$work/c4_held.csv" 5 "$work/c4_held.wav" >/dev/null; emit c4_held-r32 "$work/c4_held.wav" "$t" "$live" c4_held.csv 5
}

# Compares two --renders outputs fixture by fixture: identical or drifted. Any drift (or a missing fixture) fails; it is never re-baselined.
compare() {
  python3 -I - "${1:?usage: $0 --compare BEFORE.jsonl AFTER.jsonl}" "${2:?usage: $0 --compare BEFORE.jsonl AFTER.jsonl}" <<'PY'
import json, sys
def load(path):
    rows = {}
    for line in open(path):
        line = line.strip()
        if line:
            row = json.loads(line)
            rows[row['fixture']] = row
    return rows
before, after = load(sys.argv[1]), load(sys.argv[2])
bad = 0
for fixture in sorted(set(before) | set(after)):
    a, b = before.get(fixture), after.get(fixture)
    same = bool(a and b and a['wav_sha256'] == b['wav_sha256'] and a.get('wav_blake3') == b.get('wav_blake3') and a['bytes'] == b['bytes'] and a['bytes'] > 0)
    bad += not same
    print(json.dumps(dict(fixture=fixture, before=a and a['wav_sha256'], after=b and b['wav_sha256'], verdict='identical' if same else 'drifted')))
sys.exit(1 if bad else 0)
PY
}

run_eigen_scale() {
  [[ "$(uname -s)" != Darwin ]] || { echo 'IRP .67 scale belongs to ts1, no Mac builds' >&2; return 2; }
  cd "$REPO_ROOT"
  local output_dir
  output_dir="$(mktemp -d "${IRP_S1_RECEIPT_BASE:-${TMPDIR:-/tmp}}/frankensim-irp-s1-eigen.XXXXXX")"
  printf 'IRP .67 release receipts: %s\n' "$output_dir"
  # No ignored tests. Scale functions exist only in release builds.
  env CARGO_BUILD_JOBS=12 cargo test --release --locked -p fs-sparse --test nested_dissection g1_nd_fill_beats_amd_on_3d_grids -- --exact --nocapture --test-threads=1 2>&1 | tee "$output_dir/nd.log"
  env CARGO_BUILD_JOBS=12 cargo test --release --locked -p fs-modal --test parallel_windows -- --exact --nocapture --test-threads=1 g2_million_dof_window_timing_receipt report_elastic_p1_tet10_fill 2>&1 | tee "$output_dir/windows.log"
  python3 -I - "$output_dir" "$REPO_ROOT" <<'PY'
import json,pathlib,platform,subprocess,sys
out,root=map(pathlib.Path,sys.argv[1:])
sha=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()
rows=[]
for name in ('nd.log','windows.log'):
    for line in (out/name).read_text().splitlines():
        at=line.find('{"suite":"irp-s1"')
        if at<0:continue
        row=json.loads(line[at:]);row.update(git_sha=sha,host=platform.node(),command='scripts/e2e/irp_s1_structures.sh --run-eigen-scale',profile='release',native_exit_code=0)
        rows.append(row)
assert len([r for r in rows if r.get('test')=='g1_nd_fill_beats_amd_on_3d_grids'])==3
assert len([r for r in rows if r.get('test')=='g2_million_dof_window_timing_receipt'])==2
with (out/'events.jsonl').open('x') as f:
    for row in rows:f.write(json.dumps(row,allow_nan=False)+'\n')
PY
}

case "${1:---check}" in
  --list) list ;;
  --check) check ;;
  --run) run ;;
  --run-eigen-scale) run_eigen_scale ;;
  --renders) renders "${2:-}" ;;
  --compare) compare "${2:-}" "${3:-}" ;;
  *) printf 'usage: %s [--list|--check|--run|--run-eigen-scale|--renders OUT.jsonl|--compare BEFORE.jsonl AFTER.jsonl]\n' "$0" >&2; exit 2 ;;
esac
