#!/usr/bin/env bash
# IRP S0 materials & media subset, following material_sources.sh receipt conventions.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
readonly REPO_ROOT
readonly -a PRONY_TENSOR_TESTS=(
  g0_admission_requires_pd_equilibrium_psd_branches_positive_tau
  g0_random_multiaxial_histories_never_create_energy
  g1_isotropic_proportional_law_reduces_to_scalar_maxwell_bitwise
  g1_kelvin_constructor_matches_general_form
  g3_rotation_preserves_spectrum_and_dissipation
  g1_kramers_kronig_consistency
  g1_relaxation_function_fourier_matches_complex_modulus
  g0_plane_stress_schur_complement_is_passive_and_matches_elastic_limit
  g0_out_of_band_refuses
)
readonly -a PRONY_TENSOR_LOWERING_TESTS=(
  g0_scalar_target_reproduces_lower_to_prony_bitwise
  g1_kelvin_fractional_target_lowers_per_mode_within_certificate
  g1_complex_engineering_constants_target_lowers_with_psd_branches
  g1_plane_stress_schur_target_lowers
  g1_certificate_encloses_dense_independent_check
  g3_dropping_any_branch_breaks_the_certificate
  g0_term_count_and_band_refusals
  g0_lowered_tensor_is_bitwise_deterministic
)
readonly -a HYPERELASTIC_TESTS=(
  g2_treloar_ogden1972_reference_curves
  g1_hgo_fibre_uniaxial_analytic
  g0_reference_state_and_coalescent_stretches
)
# M0.3-W3: prescribed active fibre stress (the published-curve gate joins this list when its source is ledgered).
readonly -a ACTIVE_TISSUE_TESTS=(
  g1_isometric_active_stress_analytic
  g1_free_contraction_equilibrium
)
# M0.2-W5 certified modal damping: fs-modal owner plus the plate's own matrix.
readonly -a MODAL_DAMPING_TESTS=(
  g0_certified_radius_encloses_dense_truth
  g0_rouche_radius_is_second_order
  g1_clustered_pair_triggers_full_block_fallback
  g0_cluster_cap_and_cap_plus_one
  g0_budget_inputs_are_required_and_finite
  g0_bessel_far_coupling_bound_dominates_true_norm
  g0_windowed_certificate_encloses_dense_truth
  g0_window_neighbours_are_exactly_the_entries_read
  g0_far_row_bound_path_encloses_dense_truth_and_refuses_bad_bounds
  g0_residual_far_row_bound_encloses_exact_far_norm
)
readonly PLATE_DAMPING_MODULE=thin_plate::material_loss::ports::tensor_tests
readonly -a PLATE_DAMPING_TESTS=(
  g1_proportional_plate_is_all_diagonal_and_matches_energy_quotient
  g1_clustered_pair_triggers_full_block_fallback
  g0_bessel_far_coupling_bound_dominates_true_norm
)

list() {
  local test
  for test in "${PRONY_TENSOR_TESTS[@]}"; do
    printf 'fs-material:prony_tensor::%s\n' "$test"
  done
  for test in "${PRONY_TENSOR_LOWERING_TESTS[@]}"; do
    printf 'fs-material:prony_tensor_lowering::%s\n' "$test"
  done
  for test in "${HYPERELASTIC_TESTS[@]}"; do
    printf 'fs-material:hyperelastic_ogden_hgo::%s\n' "$test"
  done
  for test in "${ACTIVE_TISSUE_TESTS[@]}"; do
    printf 'fs-material:active_tissue::%s\n' "$test"
  done
  for test in "${MODAL_DAMPING_TESTS[@]}"; do
    printf 'fs-modal:certified_modal_damping::%s\n' "$test"
  done
  for test in "${PLATE_DAMPING_TESTS[@]}"; do
    printf 'fs-couple:lib:%s::%s\n' "$PLATE_DAMPING_MODULE" "$test"
  done
}

check() {
  local test
  for test in "${PRONY_TENSOR_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-material/tests/prony_tensor.rs" || {
      printf 'missing test: %s\n' "$test" >&2
      return 1
    }
  done
  for test in "${PRONY_TENSOR_LOWERING_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-material/tests/prony_tensor_lowering.rs" || {
      printf 'missing test: %s\n' "$test" >&2
      return 1
    }
  done
  for test in "${HYPERELASTIC_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-material/tests/hyperelastic_ogden_hgo.rs" || {
      printf 'missing test: %s\n' "$test" >&2
      return 1
    }
  done
  for test in "${ACTIVE_TISSUE_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-material/tests/active_tissue.rs" || {
      printf 'missing test: %s\n' "$test" >&2
      return 1
    }
  done
  for test in "${MODAL_DAMPING_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-modal/tests/certified_modal_damping.rs" || {
      printf 'missing test: %s\n' "$test" >&2
      return 1
    }
  done
  for test in "${PLATE_DAMPING_TESTS[@]}"; do
    grep -Fq "fn ${test}(" "${REPO_ROOT}/crates/fs-couple/src/thin_plate/material_loss/ports/tensor_tests.rs" || {
      printf 'missing test: %s\n' "$test" >&2
      return 1
    }
  done
  printf 'IRP S0: %s prony_tensor + %s prony_tensor_lowering + %s hyperelastic + %s active tissue + %s modal damping + %s plate damping tests present\n' \
    "${#PRONY_TENSOR_TESTS[@]}" "${#PRONY_TENSOR_LOWERING_TESTS[@]}" "${#HYPERELASTIC_TESTS[@]}" "${#ACTIVE_TISSUE_TESTS[@]}" \
    "${#MODAL_DAMPING_TESTS[@]}" "${#PLATE_DAMPING_TESTS[@]}"
}

run() {
  check
  local output_dir started git_head host producer_exit=0 status=0 test missing=0
  output_dir="$(mktemp -d "${IRP_S0_RECEIPT_BASE:-${TMPDIR:-/tmp}}/frankensim-irp-s0.XXXXXX")"
  started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  git_head="$(git -C "$REPO_ROOT" rev-parse HEAD)"
  host="$(uname -n)"
  printf 'IRP S0 live receipt: %s\n' "$output_dir"

  local cargo_bin="${CARGO_BIN:-cargo}"
  if [[ -x "/Users/jemanuel/.rustup/toolchains/nightly-2026-07-06-aarch64-apple-darwin/bin/cargo" ]]; then
    cargo_bin="/Users/jemanuel/.rustup/toolchains/nightly-2026-07-06-aarch64-apple-darwin/bin/cargo"
  fi
  export RCH_CARGO_WRAPPER_BYPASS=1
  export CARGO_BUILD_JOBS=12
  export CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test"

  local -a prony_command lowering_command hyperelastic_command active_command modal_command plate_command
  local -a plate_filters=()
  for test in "${PLATE_DAMPING_TESTS[@]}"; do
    plate_filters+=("${PLATE_DAMPING_MODULE}::${test}")
  done
  if [[ -n "${IRP_S0_PRONY_TENSOR_BIN:-}" ]]; then
    prony_command=("$IRP_S0_PRONY_TENSOR_BIN" --exact --test-threads=1 --show-output "${PRONY_TENSOR_TESTS[@]}")
  else
    prony_command=("$cargo_bin" test --locked -j 12 -p fs-material --test prony_tensor -- --exact --test-threads=1 --show-output "${PRONY_TENSOR_TESTS[@]}")
  fi
  if [[ -n "${IRP_S0_PRONY_TENSOR_LOWERING_BIN:-}" ]]; then
    lowering_command=("$IRP_S0_PRONY_TENSOR_LOWERING_BIN" --exact --test-threads=1 --show-output "${PRONY_TENSOR_LOWERING_TESTS[@]}")
  else
    lowering_command=("$cargo_bin" test --locked -j 12 -p fs-material --test prony_tensor_lowering -- --exact --test-threads=1 --show-output "${PRONY_TENSOR_LOWERING_TESTS[@]}")
  fi
  if [[ -n "${IRP_S0_HYPERELASTIC_BIN:-}" ]]; then
    hyperelastic_command=("$IRP_S0_HYPERELASTIC_BIN" --exact --test-threads=1 --show-output "${HYPERELASTIC_TESTS[@]}")
  else
    hyperelastic_command=("$cargo_bin" test --locked -j 12 -p fs-material --test hyperelastic_ogden_hgo -- --exact --test-threads=1 --show-output "${HYPERELASTIC_TESTS[@]}")
  fi
  if [[ -n "${IRP_S0_ACTIVE_TISSUE_BIN:-}" ]]; then
    active_command=("$IRP_S0_ACTIVE_TISSUE_BIN" --exact --test-threads=1 --show-output "${ACTIVE_TISSUE_TESTS[@]}")
  else
    active_command=("$cargo_bin" test --locked -j 12 -p fs-material --test active_tissue -- --exact --test-threads=1 --show-output "${ACTIVE_TISSUE_TESTS[@]}")
  fi
  if [[ -n "${IRP_S0_MODAL_DAMPING_BIN:-}" ]]; then
    modal_command=("$IRP_S0_MODAL_DAMPING_BIN" --exact --test-threads=1 --show-output "${MODAL_DAMPING_TESTS[@]}")
  else
    modal_command=("$cargo_bin" test --locked -j 12 -p fs-modal --test certified_modal_damping -- --exact --test-threads=1 --show-output "${MODAL_DAMPING_TESTS[@]}")
  fi
  if [[ -n "${IRP_S0_PLATE_DAMPING_BIN:-}" ]]; then
    plate_command=("$IRP_S0_PLATE_DAMPING_BIN" --exact --test-threads=1 --show-output "${plate_filters[@]}")
  else
    plate_command=("$cargo_bin" test --locked -j 12 -p fs-couple --lib -- --exact --test-threads=1 --show-output "${plate_filters[@]}")
  fi

  {
    printf 'host=%s\ngit_head=%s\nstarted_utc=%s\n' "$host" "$git_head" "$started"
    printf '%q ' "${prony_command[@]}"
    printf '\n'
    printf '%q ' "${lowering_command[@]}"
    printf '\n'
    printf '%q ' "${hyperelastic_command[@]}"
    printf '\n'
    printf '%q ' "${active_command[@]}"
    printf '\n'
    printf '%q ' "${modal_command[@]}"
    printf '\n'
    printf '%q ' "${plate_command[@]}"
    printf '\n'
  } >"${output_dir}/command.txt"

  set +e
  (cd "$REPO_ROOT" && "${prony_command[@]}") >"${output_dir}/stdout.log" 2>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  (cd "$REPO_ROOT" && "${lowering_command[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  (cd "$REPO_ROOT" && "${hyperelastic_command[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  (cd "$REPO_ROOT" && "${active_command[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  (cd "$REPO_ROOT" && "${modal_command[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  (cd "$REPO_ROOT" && "${plate_command[@]}") >>"${output_dir}/stdout.log" 2>>"${output_dir}/stderr.log"
  status=$?
  if [[ "$status" -ne 0 ]]; then producer_exit="$status"; fi
  set -e

  for test in "${PRONY_TENSOR_TESTS[@]}" "${PRONY_TENSOR_LOWERING_TESTS[@]}" "${HYPERELASTIC_TESTS[@]}" \
    "${ACTIVE_TISSUE_TESTS[@]}" "${MODAL_DAMPING_TESTS[@]}" "${plate_filters[@]}"; do
    if ! grep -Fqx "test ${test} ... ok" "${output_dir}/stdout.log" "${output_dir}/stderr.log"; then
      printf 'missing or failed: %s\n' "$test" >&2
      missing=1
    fi
  done
  if ! grep -Eq '^test result: ok\. [1-9][0-9]* passed' "${output_dir}/stdout.log" "${output_dir}/stderr.log"; then
    missing=1
  fi

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
        if not isinstance(row, dict) or row.get("suite") not in ("irp-s0", "irp-s0-materials"):
            continue
        row.update(
            host=host,
            git_head=head,
            command=command,
            utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        )
        events.append(row)
with (directory / "events.jsonl").open("w") as output:
    for event in events:
        output.write(json.dumps(event, allow_nan=False) + "\n")
passed = producer_exit == 0 and missing == 0 and bool(events)
summary = dict(
    schema="frankensim.irp-s0-materials.e2e.v1",
    host=host,
    git_head=head,
    started_utc=started,
    producer_exit=producer_exit,
    missing_or_failed=missing,
    events=len(events),
    command_file="command.txt",
    stdout_log="stdout.log",
    stderr_log="stderr.log",
    events_log="events.jsonl",
    verdict="pass" if passed else "failed",
)
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
