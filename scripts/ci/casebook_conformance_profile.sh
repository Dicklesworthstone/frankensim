#!/usr/bin/env bash
# Casebook conformance profile runner.
#
# The PR profile is an explicit cheap selector. The nightly-full profile is
# the complete source-discovered inventory of ordinary (non-ignored) Cargo
# integration targets that contain an fs_casebook source token, protected by a
# reviewed minimum baseline. Locked Cargo metadata is used only for target
# discovery; no filename convention is trusted.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
readonly SCRIPT_DIR
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
readonly REPO_ROOT
readonly CARGO_BIN="${CARGO_BIN:-cargo}"
readonly PR_BUDGET_SECONDS="${FS_CASEBOOK_PR_BUDGET_SECONDS:-900}"
readonly FULL_BUDGET_SECONDS="${FS_CASEBOOK_FULL_BUDGET_SECONDS:-7200}"
readonly TERMINATION_GRACE_SECONDS="${FS_CASEBOOK_TERMINATION_GRACE_SECONDS:-5}"
readonly KILL_DRAIN_SECONDS="${FS_CASEBOOK_KILL_DRAIN_SECONDS:-5}"

# Deliberately cheap, family-representative merge selectors. Every entry must
# also be in REQUIRED_FULL_TARGETS and in live metadata/source-token discovery.
readonly -a PR_TARGETS=(
  "fs-ad:conformance"
  "fs-casebook:casebook"
  "fs-cheb:conformance"
  "fs-fft:conformance"
  "fs-ivl:structured_conformance_casebook"
  "fs-la:conformance"
  "fs-math:conformance"
  "fs-rand:conformance"
  "fs-simd:conformance"
  "fs-sparse:conformance"
)

# Reviewed minimum Casebook inventory. The full profile is source-derived and
# auto-adopts additional discovered targets; this baseline prevents a scanner
# regression or accidental target removal from silently shrinking coverage.
readonly -a REQUIRED_FULL_TARGETS=(
  "fs-ad:conformance"
  "fs-ad:la_dual_bridge_casebook"
  "fs-archive:conformance"
  "fs-ascent:frankenscipy_optimizer_oracle_casebook"
  "fs-bo:bo_study_replay"
  "fs-bo:mf_study_replay"
  "fs-bo:turbo_study_replay"
  "fs-casebook:casebook"
  "fs-cheb:conformance"
  "fs-cheb:dct_bridge_casebook"
  "fs-cheb:frankenscipy_integrate_oracle_casebook"
  "fs-fft:conformance"
  "fs-fft:frankenscipy_fft_oracle_casebook"
  "fs-ir:conformance_ir"
  "fs-ivl:eft_interval_bridge_casebook"
  "fs-ivl:structured_conformance_casebook"
  "fs-la:conformance"
  "fs-la:eigen_replay_casebook"
  "fs-la:frankenscipy_linalg_oracle_casebook"
  "fs-la:rand_gemm_replay_casebook"
  "fs-la:rand_nla_casebook"
  "fs-math:conformance"
  "fs-math:frankenscipy_special_oracle_casebook"
  "fs-rand:conformance"
  "fs-rand:qmc_replay_casebook"
  "fs-simd:conformance"
  "fs-sparse:conformance"
  "fs-sparse:frankenscipy_oracle_casebook"
  "fs-sparse:preconditioner_casebook"
  "fs-time:frankenscipy_ode_oracle_casebook"
)

usage() {
  cat >&2 <<'EOF'
usage:
  bash scripts/ci/casebook_conformance_profile.sh --check
  bash scripts/ci/casebook_conformance_profile.sh --list pr
  bash scripts/ci/casebook_conformance_profile.sh --list nightly-full
  bash scripts/ci/casebook_conformance_profile.sh pr
  bash scripts/ci/casebook_conformance_profile.sh nightly-full
  bash scripts/ci/casebook_conformance_profile.sh --self-test

environment:
  CARGO_BIN                              Cargo executable (default: cargo)
  FS_CASEBOOK_PR_BUDGET_SECONDS          Aggregate PR wall budget (default: 900)
  FS_CASEBOOK_FULL_BUDGET_SECONDS        Aggregate full wall budget (default: 7200)
  FS_CASEBOOK_TERMINATION_GRACE_SECONDS  TERM grace period (default: 5)
  FS_CASEBOOK_KILL_DRAIN_SECONDS          KILL drain period (default: 5)

--check and --list run only `cargo metadata --locked --no-deps` plus source
inspection; they do not build or execute tests. --self-test uses synthetic
inventories and a fake expired deadline and never invokes Cargo.
EOF
}

require_positive_integer() {
  local label="$1"
  local value="$2"
  if [[ ! "${value}" =~ ^[1-9][0-9]*$ ]]; then
    printf 'invalid %s: expected a positive integer, got %q\n' "${label}" "${value}" >&2
    return 2
  fi
}

profile_budget() {
  case "$1" in
    pr) printf '%s\n' "${PR_BUDGET_SECONDS}" ;;
    nightly-full) printf '%s\n' "${FULL_BUDGET_SECONDS}" ;;
    *) return 2 ;;
  esac
}

emit_registry_rows() {
  local target
  for target in "${PR_TARGETS[@]}"; do
    printf 'pr\t%s\n' "${target}"
  done
  for target in "${REQUIRED_FULL_TARGETS[@]}"; do
    printf 'required\t%s\n' "${target}"
  done
}

discover_casebook_targets() {
  "${CARGO_BIN}" metadata --locked --format-version 1 --no-deps | python3 -c '
import json
import pathlib
import re
import sys

meta = json.load(sys.stdin)
root = pathlib.Path(meta["workspace_root"]).resolve()
token = re.compile(rb"(?:\buse\s+fs_casebook\b|\bfs_casebook\s*::)")
ignored = re.compile(rb"^\s*#\s*\[\s*ignore(?:\s|=|\])", re.MULTILINE)
rows = []
errors = []

for package in meta["packages"]:
    package_name = package["name"]
    for target in package["targets"]:
        if "test" not in target["kind"]:
            continue
        target_name = target["name"]
        source = pathlib.Path(target["src_path"]).resolve()
        try:
            source.relative_to(root)
        except ValueError:
            errors.append(f"target source escapes workspace: {package_name}:{target_name}={source}")
            continue
        data = source.read_bytes()
        if not token.search(data):
            continue
        identity = f"{package_name}:{target_name}"
        if ignored.search(data):
            errors.append(
                f"{identity} mixes #[ignore] with Casebook coverage; classify the ignored lane separately"
            )
            continue
        features = ",".join(sorted(target.get("required-features") or [])) or "-"
        rows.append((identity, features, source.relative_to(root).as_posix()))

if errors:
    for error in sorted(errors):
        print(error, file=sys.stderr)
    raise SystemExit(1)

for identity, features, source in sorted(rows):
    print(f"discovered\t{identity}\t{features}\t{source}")
'
}

# Input rows are tagged as pr/required/discovered. The validator is deliberately
# reusable by --self-test so the negative fixtures exercise the live policy.
validate_inventory_payload() {
  python3 -c '
import collections
import json
import re
import sys

valid_identity = re.compile(r"^[A-Za-z0-9_-]+:[A-Za-z0-9_-]+$")
groups = {"pr": [], "required": [], "discovered": []}
errors = []

for number, raw in enumerate(sys.stdin, 1):
    line = raw.rstrip("\n")
    if not line:
        continue
    fields = line.split("\t")
    kind = fields[0]
    if kind not in groups:
        errors.append({"code": "malformed_row", "detail": f"line {number}: {line!r}"})
        continue
    expected_fields = 4 if kind == "discovered" else 2
    if len(fields) != expected_fields:
        errors.append({"code": f"malformed_{kind}", "detail": f"line {number}: {line!r}"})
        continue
    identity = fields[1]
    if not valid_identity.fullmatch(identity):
        errors.append({"code": "invalid_identity", "detail": identity})
        continue
    groups[kind].append(identity)
    if kind == "discovered":
        features, source = fields[2:]
        if not features or not source:
            errors.append({"code": "malformed_discovered", "detail": identity})

for kind, values in groups.items():
    for identity, count in sorted(collections.Counter(values).items()):
        if count > 1:
            errors.append({"code": f"duplicate_{kind}", "detail": identity})

pr = set(groups["pr"])
required = set(groups["required"])
discovered = set(groups["discovered"])
for identity in sorted(pr - discovered):
    errors.append({"code": "missing_pr", "detail": identity})
for identity in sorted(pr - required):
    errors.append({"code": "pr_not_required", "detail": identity})
for identity in sorted(required - discovered):
    errors.append({"code": "stale_required", "detail": identity})
if not pr:
    errors.append({"code": "empty_pr", "detail": "PR profile has zero targets"})
if not required:
    errors.append({"code": "empty_required", "detail": "required baseline has zero targets"})
if not discovered:
    errors.append({"code": "empty_full", "detail": "full discovery has zero targets"})

receipt = {
    "schema": "frankensim-casebook-profile-inventory-v1",
    "status": "fail" if errors else "pass",
    "pr_targets": len(pr),
    "required_full_targets": len(required),
    "full_targets": len(discovered),
    "discovered_targets": len(discovered),
    "errors": sorted(errors, key=lambda item: (item["code"], item["detail"])),
}
print(json.dumps(receipt, sort_keys=True, separators=(",", ":")))
raise SystemExit(1 if errors else 0)
'
}

DISCOVERY_ROWS=""

refresh_and_validate_inventory() {
  local receipt
  DISCOVERY_ROWS="$(discover_casebook_targets)" || return $?
  if receipt="$({ emit_registry_rows; printf '%s\n' "${DISCOVERY_ROWS}"; } | validate_inventory_payload)"; then
    printf '%s\n' "${receipt}"
  else
    local status=$?
    printf '%s\n' "${receipt}" >&2
    return "${status}"
  fi
}

lookup_discovery() {
  local wanted="$1"
  local kind identity features source
  LOOKUP_FEATURES=""
  LOOKUP_SOURCE=""
  while IFS=$'\t' read -r kind identity features source; do
    if [[ "${kind}" == "discovered" && "${identity}" == "${wanted}" ]]; then
      if [[ "${features}" == "-" ]]; then
        LOOKUP_FEATURES=""
      else
        LOOKUP_FEATURES="${features}"
      fi
      LOOKUP_SOURCE="${source}"
      return 0
    fi
  done <<< "${DISCOVERY_ROWS}"
  return 1
}

profile_targets() {
  local profile="$1"
  local target
  if [[ "${profile}" == "pr" ]]; then
    for target in "${PR_TARGETS[@]}"; do
      printf '%s\n' "${target}"
    done
  else
    local kind identity features source
    while IFS=$'\t' read -r kind identity features source; do
      if [[ "${kind}" == "discovered" ]]; then
        printf '%s\n' "${identity}"
      fi
    done <<< "${DISCOVERY_ROWS}"
  fi
}

emit_selector() {
  local profile="$1"
  local identity="$2"
  local package="${identity%%:*}"
  local target="${identity#*:}"
  lookup_discovery "${identity}"
  python3 - "${profile}" "${package}" "${target}" "${LOOKUP_FEATURES}" \
    "${LOOKUP_SOURCE}" "${CARGO_BIN}" <<'PY'
import json
import sys

profile, package, target, features, source, cargo_bin = sys.argv[1:]
command = [cargo_bin, "test", "--locked", "-p", package, "--test", target]
if features:
    command.extend(["--features", features])
command.extend(["--", "--nocapture"])
print(json.dumps({
    "schema": "frankensim-casebook-profile-selector-v1",
    "profile": profile,
    "package": package,
    "target": target,
    "required_features": features.split(",") if features else [],
    "source": source,
    "command": command,
}, sort_keys=True, separators=(",", ":")))
PY
}

list_profile() {
  local profile="$1"
  local budget target
  budget="$(profile_budget "${profile}")" || return 2
  printf '{"schema":"frankensim-casebook-profile-v1","profile":"%s","budget_seconds":%s,"build_time_included":true,"deadline_enforced":true,"ignored_tests_included":false}\n' \
    "${profile}" "${budget}"
  while IFS= read -r target; do
    emit_selector "${profile}" "${target}"
  done < <(profile_targets "${profile}")
}

deadline_refusal_receipt() {
  local profile="$1"
  local identity="$2"
  local now="$3"
  local deadline="$4"
  local package="${identity%%:*}"
  local target="${identity#*:}"
  if (( now < deadline )); then
    return 0
  fi
  printf '{"schema":"frankensim-casebook-profile-target-v1","profile":"%s","package":"%s","target":"%s","status":"budget_exceeded","launched":false,"exit_code":null,"budget_status":"exceeded","drain_status":"not_applicable","drained_process_group":false}\n' \
    "${profile}" "${package}" "${target}"
  return 124
}

# Execute one selector in its own session/process group. Python enforces the
# absolute aggregate deadline, TERM -> bounded wait -> KILL -> bounded drain,
# and emits the complete target receipt. It never targets unrelated processes.
run_target_until_deadline() {
  local deadline_epoch="$1"
  local profile="$2"
  local identity="$3"
  shift 3
  python3 - "${deadline_epoch}" "${TERMINATION_GRACE_SECONDS}" \
    "${KILL_DRAIN_SECONDS}" "${profile}" "${identity}" "$@" <<'PY'
import json
import os
import signal
import subprocess
import sys
import time

deadline_epoch = int(sys.argv[1])
term_grace = int(sys.argv[2])
kill_grace = int(sys.argv[3])
profile = sys.argv[4]
identity = sys.argv[5]
command = sys.argv[6:]
package, target = identity.split(":", 1)
started = time.monotonic()

def emit_receipt(**fields):
    receipt = {
        "schema": "frankensim-casebook-profile-target-v1",
        "profile": profile,
        "package": package,
        "target": target,
        "command": command,
    }
    receipt.update(fields)
    print(json.dumps(receipt, sort_keys=True, separators=(",", ":")))

if time.time() >= deadline_epoch:
    emit_receipt(
        status="budget_exceeded",
        launched=False,
        exit_code=None,
        leader_exit_code=None,
        elapsed_seconds=0,
        budget_status="exceeded",
        drain_status="not_applicable",
        drain_trigger="deadline_before_spawn",
        drained_process_group=False,
    )
    raise SystemExit(124)

try:
    process = subprocess.Popen(command, start_new_session=True)
except OSError as error:
    emit_receipt(
        status="fail",
        launched=False,
        exit_code=126,
        leader_exit_code=None,
        elapsed_seconds=max(0, int(time.monotonic() - started)),
        budget_status="within",
        drain_status="not_applicable",
        drain_trigger="spawn_failure",
        drained_process_group=False,
        spawn_error=f"{type(error).__name__}: {error}",
    )
    raise SystemExit(126)

def group_is_running():
    try:
        os.killpg(process.pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True

def wait_for_group(seconds):
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        process.poll()
        if not group_is_running():
            return True
        time.sleep(0.05)
    process.poll()
    return not group_is_running()

def drain_owned_group():
    if not group_is_running():
        return "complete"
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    if not wait_for_group(term_grace):
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        wait_for_group(kill_grace)
    try:
        process.wait(timeout=0.1)
    except subprocess.TimeoutExpired:
        pass
    return (
        "complete"
        if process.poll() is not None and not group_is_running()
        else "incomplete"
    )

remaining_seconds = max(0.0, deadline_epoch - time.time())
timed_out = False
try:
    leader_exit_code = process.wait(timeout=remaining_seconds)
except subprocess.TimeoutExpired:
    timed_out = True

if timed_out:
    drain_status = drain_owned_group()
    leader_exit_code = process.poll()
    exit_code = 124 if drain_status == "complete" else 125
    status = "budget_exceeded"
    budget_status = "exceeded"
    drain_trigger = "deadline"
    wrapper_exit_code = exit_code
elif group_is_running():
    drain_status = drain_owned_group()
    exit_code = 1 if drain_status == "complete" else 127
    status = "fail"
    budget_status = "within"
    drain_trigger = "leader_exit_with_live_group"
    wrapper_exit_code = exit_code
else:
    drain_status = "not_needed"
    exit_code = leader_exit_code
    status = "pass" if leader_exit_code == 0 else "fail"
    budget_status = "within"
    drain_trigger = "none"
    wrapper_exit_code = 0 if leader_exit_code == 0 else 1

emit_receipt(
    status=status,
    launched=True,
    exit_code=exit_code,
    leader_exit_code=leader_exit_code,
    elapsed_seconds=max(0, int(time.monotonic() - started)),
    budget_status=budget_status,
    drain_status=drain_status,
    drain_trigger=drain_trigger,
    drained_process_group=drain_status == "complete",
)
raise SystemExit(wrapper_exit_code)
PY
}

run_profile() {
  local profile="$1"
  local budget="$2"
  local started deadline now identity package target
  local status overall_status=0 budget_status="within" run_status="pass"
  local completed=0 failed=0
  local -a command=()
  started="$(date +%s)"
  deadline="$((started + budget))"
  printf '{"schema":"frankensim-casebook-profile-run-v1","event":"start","profile":"%s","budget_seconds":%s,"deadline_epoch_seconds":%s,"build_time_included":true,"deadline_enforced":true,"ignored_tests_included":false}\n' \
    "${profile}" "${budget}" "${deadline}"

  while IFS= read -r identity; do
    now="$(date +%s)"
    if deadline_refusal_receipt "${profile}" "${identity}" "${now}" "${deadline}"; then
      :
    else
      status=$?
      if (( status == 124 )); then
        overall_status=1
        budget_status="exceeded"
        failed=$((failed + 1))
        break
      fi
      return "${status}"
    fi
    package="${identity%%:*}"
    target="${identity#*:}"
    lookup_discovery "${identity}"
    command=("${CARGO_BIN}" test --locked -p "${package}" --test "${target}")
    if [[ -n "${LOOKUP_FEATURES}" ]]; then
      command+=(--features "${LOOKUP_FEATURES}")
    fi
    command+=(-- --nocapture)
    if run_target_until_deadline "${deadline}" "${profile}" "${identity}" "${command[@]}"; then
      status=0
      completed=$((completed + 1))
    else
      status=$?
      failed=$((failed + 1))
      overall_status=1
      if (( status == 124 || status == 125 )); then
        budget_status="exceeded"
        break
      elif (( status == 127 )); then
        break
      fi
    fi
  done < <(profile_targets "${profile}")

  now="$(date +%s)"
  if [[ "${budget_status}" == "exceeded" ]]; then
    run_status="budget_exceeded"
  elif (( overall_status != 0 )); then
    run_status="fail"
  fi
  printf '{"schema":"frankensim-casebook-profile-run-v1","event":"finish","profile":"%s","status":"%s","budget_status":"%s","budget_seconds":%s,"elapsed_seconds":%s,"completed_targets":%s,"failed_targets":%s}\n' \
    "${profile}" "${run_status}" "${budget_status}" "${budget}" \
    "$((now - started))" "${completed}" "${failed}"
  return "${overall_status}"
}

expect_validation_failure() {
  local expected_code="$1"
  local payload="$2"
  local output status
  if output="$(printf '%s\n' "${payload}" | validate_inventory_payload)"; then
    printf 'self-test expected %s failure, validator passed: %s\n' \
      "${expected_code}" "${output}" >&2
    return 1
  else
    status=$?
  fi
  if (( status == 0 )) || [[ "${output}" != *"\"code\":\"${expected_code}\""* ]]; then
    printf 'self-test expected %s, observed: %s\n' "${expected_code}" "${output}" >&2
    return 1
  fi
}

run_self_tests() {
  local base missing stale duplicate deadline output status
  base=$'pr\tfs-a:casebook\nrequired\tfs-a:casebook\nrequired\tfs-b:conformance\ndiscovered\tfs-a:casebook\t-\tcrates/fs-a/tests/casebook.rs\ndiscovered\tfs-b:conformance\t-\tcrates/fs-b/tests/conformance.rs'
  missing="${base}"$'\npr\tfs-c:missing_case\nrequired\tfs-c:missing_case'
  stale="${base}"$'\nrequired\tfs-c:removed_case'
  duplicate="${base}"$'\nrequired\tfs-b:conformance'
  expect_validation_failure "missing_pr" "${missing}"
  expect_validation_failure "stale_required" "${stale}"
  expect_validation_failure "duplicate_required" "${duplicate}"

  deadline=100
  if output="$(deadline_refusal_receipt pr fs-a:casebook 100 "${deadline}")"; then
    printf 'self-test deadline refusal unexpectedly allowed launch\n' >&2
    return 1
  else
    status=$?
  fi
  if (( status != 124 )) || [[ "${output}" != *'"launched":false'* ]] \
      || [[ "${output}" != *'"drain_status":"not_applicable"'* ]]; then
    printf 'self-test deadline refusal receipt mismatch: %s\n' "${output}" >&2
    return 1
  fi
  printf '%s\n' '{"schema":"frankensim-casebook-profile-self-test-v1","status":"pass","cases":4,"cargo_invocations":0}'
}

if (( $# == 1 )) && [[ "$1" == "--self-test" ]]; then
  run_self_tests
  exit $?
fi

require_positive_integer "PR budget" "${PR_BUDGET_SECONDS}"
require_positive_integer "full budget" "${FULL_BUDGET_SECONDS}"
require_positive_integer "termination grace" "${TERMINATION_GRACE_SECONDS}"
require_positive_integer "kill drain" "${KILL_DRAIN_SECONDS}"
cd "${REPO_ROOT}"

if (( $# == 1 )) && [[ "$1" == "--check" ]]; then
  refresh_and_validate_inventory
  exit $?
fi

if (( $# == 2 )) && [[ "$1" == "--list" ]]; then
  profile_budget "$2" >/dev/null || {
    usage
    exit 2
  }
  refresh_and_validate_inventory
  list_profile "$2"
  exit $?
fi

if (( $# != 1 )) || [[ "$1" != "pr" && "$1" != "nightly-full" ]]; then
  usage
  exit 2
fi

readonly PROFILE="$1"
BUDGET_SECONDS="$(profile_budget "${PROFILE}")"
readonly BUDGET_SECONDS
refresh_and_validate_inventory
run_profile "${PROFILE}" "${BUDGET_SECONDS}"
exit $?
