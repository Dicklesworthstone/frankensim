#!/usr/bin/env bash
# Focused source-pack compiler/store/consumer E2E receipt.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
readonly REPO_ROOT
readonly TEST_FILE="${REPO_ROOT}/xtask/tests/matdb_pack_cli.rs"
readonly -a EXPECTED_TESTS=(
  "common_material_acquisition::g0_g3_cryogenic_aluminum_copper_curves_and_discovery"
  "common_material_acquisition::g0_g3_iapws_liquid_water_curves_and_discovery"
  "common_material_acquisition::g1_g3_sourced_water_conduction_at_two_states"
  "common_material_acquisition::g1_g3_sourced_liquid_water_enthalpy_heating"
  "common_material_acquisition::g1_g3_sourced_stainless_store_to_conduction"
  "common_material_acquisition::g1_g3_sourced_warm_316_reaches_conduction"
  "common_material_acquisition::g0_g3_stainless_thermomechanical_envelope"
  "common_material_acquisition::g1_g3_sourced_silicon_tensor_reaches_oriented_solid"
  "common_material_acquisition::g1_g3_sourced_dry_air_reaches_acoustic_loss"
  "common_material_acquisition::g1_g3_sourced_metal_profiles_reach_thermoelastic_plate"
  "common_material_acquisition::g1_g3_sourced_conductors_reach_circuit_dissipation"
  "common_material_acquisition::g1_g3_sourced_humid_air_reaches_acoustic_transport"
  "common_material_acquisition::g1_g3_sourced_glycols_reach_heat_and_flow"
  "common_material_acquisition::g1_g3_sourced_mineral_oil_reaches_heat_and_flow"
  "common_material_acquisition::g0_g3_sourced_ensinger_tecafine_pe300_natural_2017_observations"
  "glass_reference::g1_g3_glass_reference_profiles_reach_heat"
  "common_material_acquisition::g1_g3_source_card_slab_proteus_pp"
  "common_material_acquisition::g1_g3_polymer_reference_comparison_reaches_heat"
  "common_material_acquisition::g1_g3_housing_polymers_source_card_slab"
)

usage() {
  printf 'usage: %s [--list|--check|--run]\n' "$0" >&2
}

list_cases() {
  printf '%s\n' "${EXPECTED_TESTS[@]}"
}

check() {
  [[ -f "$TEST_FILE" ]] || { printf 'missing %s\n' "$TEST_FILE" >&2; return 1; }
  if [[ -n "${MATDB_PACK_TEST_BIN:-}" ]]; then
    [[ -x "$MATDB_PACK_TEST_BIN" ]] || { printf 'MATDB_PACK_TEST_BIN is not executable\n' >&2; return 1; }
  else
    command -v rch >/dev/null || { printf 'rch is required for --run\n' >&2; return 1; }
    command -v cargo >/dev/null || { printf 'cargo is required for --run\n' >&2; return 1; }
    rch exec --help 2>&1 | grep -Fq -- '--no-self-healing' || {
      printf 'rch exec lacks --no-self-healing\n' >&2; return 1;
    }
  fi
  local test
  for test in "${EXPECTED_TESTS[@]}"; do
    grep -Fq "fn ${test##*::}" "$TEST_FILE" || {
      printf 'expected test missing: %s\n' "$test" >&2
      return 1
    }
  done
  printf 'material_sources check: %s expected tests present\n' "${#EXPECTED_TESTS[@]}"
}

write_summary() {
  local output_dir="$1" producer_exit="$2" verdict="$3"
  {
    printf '{\n  "schema":"frankensim.material-sources.e2e.v1",\n'
    printf '  "started_utc":"%s",\n' "$STARTED_UTC"
    printf '  "producer_exit":%s,\n' "$producer_exit"
    printf '  "verdict":"%s",\n' "$verdict"
    printf '  "command_file":"command.txt",\n  "stdout_log":"stdout.log",\n'
    printf '  "stderr_log":"stderr.log",\n  "events_log":"events.jsonl",\n'
    printf '  "expected_tests":[\n'
    local index
    for index in "${!EXPECTED_TESTS[@]}"; do
      printf '    "%s"%s\n' "${EXPECTED_TESTS[$index]}" \
        "$([[ "$index" -eq $((${#EXPECTED_TESTS[@]} - 1)) ]] && printf '' || printf ',')"
    done
    printf '  ]\n}\n'
  } >"${output_dir}/summary.json"
}

run() {
  check
  local output_dir started command_source producer_exit missing=0 test
  output_dir="$(mktemp -d "${TMPDIR:-/tmp}/frankensim-material-sources.XXXXXX")"
  started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  STARTED_UTC="$started"
  export STARTED_UTC
  printf '{"utc":"%s","event":"start","expected_test_count":%s}\n' \
    "$started" "${#EXPECTED_TESTS[@]}" >"${output_dir}/events.jsonl"

  local -a command
  if [[ -n "${MATDB_PACK_TEST_BIN:-}" ]]; then
    command=("$MATDB_PACK_TEST_BIN" --exact --test-threads=1 --show-output "${EXPECTED_TESTS[@]}")
    command_source="supplied-native-test-binary"
  else
    export RCH_REQUIRE_REMOTE=1
    command=(rch exec --no-self-healing -- env \
      "CARGO_TARGET_DIR=${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch-greenosprey-study" \
      cargo +nightly-2026-07-06 test --locked -j 4 -p xtask --test matdb_pack_cli -- \
      --exact --test-threads=1 --show-output "${EXPECTED_TESTS[@]}")
    command_source="rch-required-remote-cargo"
  fi
  {
    printf 'utc=%s\nsource=%s\nRCH_REQUIRE_REMOTE=%s\ncommand=' \
      "$started" "$command_source" "${RCH_REQUIRE_REMOTE:-unset}"
    printf '%q ' "${command[@]}"
    printf '\n'
  } >"${output_dir}/command.txt"
  printf 'material source E2E receipt (live logs): %s\n' "$output_dir"

  set +e
  (cd "$REPO_ROOT" && "${command[@]}") >"${output_dir}/stdout.log" 2>"${output_dir}/stderr.log"
  producer_exit=$?
  set -e
  printf '{"utc":"%s","event":"producer-exit","status":%s}\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$producer_exit" >>"${output_dir}/events.jsonl"

  for test in "${EXPECTED_TESTS[@]}"; do
    # RCH can relay both remote streams on its stderr; native libtest uses
    # stdout. Require the complete result line in either retained stream.
    if ! grep -Fqx "test ${test} ... ok" "${output_dir}/stdout.log" "${output_dir}/stderr.log"; then
      printf '{"utc":"%s","event":"expected-test-missing-or-failed","test":"%s"}\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$test" >>"${output_dir}/events.jsonl"
      missing=1
    fi
  done
  if ! grep -Eq '^test result: ok\. [1-9][0-9]* passed' "${output_dir}/stdout.log" "${output_dir}/stderr.log"; then
    missing=1
  fi

  if [[ "$producer_exit" -eq 0 && "$missing" -eq 0 ]]; then
    write_summary "$output_dir" "$producer_exit" "pass"
    printf 'material source E2E passed; receipt: %s\n' "$output_dir"
    return 0
  fi
  write_summary "$output_dir" "$producer_exit" "failed"
  printf 'material source E2E failed; retained receipt: %s\n' "$output_dir" >&2
  tail -n 80 "${output_dir}/stdout.log"
  tail -n 80 "${output_dir}/stderr.log" >&2
  return 1
}

case "${1:---check}" in
  --list) list_cases ;;
  --check) check ;;
  --run) run ;;
  *) usage; exit 2 ;;
esac
