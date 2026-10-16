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
readonly PYTHON_BIN="${PYTHON_BIN:-python3}"
readonly PR_BUDGET_SECONDS="${FS_CASEBOOK_PR_BUDGET_SECONDS:-900}"
readonly FULL_BUDGET_SECONDS="${FS_CASEBOOK_FULL_BUDGET_SECONDS:-7200}"
readonly TERMINATION_GRACE_SECONDS="${FS_CASEBOOK_TERMINATION_GRACE_SECONDS:-5}"
readonly KILL_DRAIN_SECONDS="${FS_CASEBOOK_KILL_DRAIN_SECONDS:-5}"
readonly MAX_PROFILE_BUDGET_SECONDS=604800
readonly MAX_DRAIN_SECONDS=300

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
  "fs-bo:hetero_gp_study_replay"
  "fs-bo:mf_study_replay"
  "fs-bo:qei_study_replay"
  "fs-bo:sparse_gp_study_replay"
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
  bash scripts/ci/casebook_conformance_profile.sh --run pr
  bash scripts/ci/casebook_conformance_profile.sh --run nightly-full
  bash scripts/ci/casebook_conformance_profile.sh --replay <receipt.jsonl>
  bash scripts/ci/casebook_conformance_profile.sh --self-test

environment:
  CARGO_BIN                              Cargo executable (default: cargo)
  PYTHON_BIN                             Python 3 executable (default: python3)
  FS_CASEBOOK_PR_BUDGET_SECONDS          Aggregate PR wall budget (default: 900)
  FS_CASEBOOK_FULL_BUDGET_SECONDS        Aggregate full wall budget (default: 7200)
  FS_CASEBOOK_TERMINATION_GRACE_SECONDS  TERM grace period (default: 5)
  FS_CASEBOOK_KILL_DRAIN_SECONDS          KILL drain period (default: 5)

--check and --list run only `cargo metadata --locked --no-deps` plus source
inspection; they do not build or execute tests. --replay validates a recorded
profile run JSONL stream against the live profile inventory without executing
tests. --self-test uses synthetic inventories and no-Cargo child processes to
exercise discovery, budgets, cleanup, signal forwarding, receipt isolation, and
replay verification; it never invokes Cargo.
EOF
}

require_bounded_positive_integer() {
  local label="$1"
  local value="$2"
  local maximum="$3"
  if [[ ! "${value}" =~ ^[1-9][0-9]*$ ]]; then
    printf 'invalid %s: expected a positive integer, got %q\n' "${label}" "${value}" >&2
    return 2
  fi
  if (( ${#value} > ${#maximum} )) \
      || { (( ${#value} == ${#maximum} )) && (( 10#${value} > 10#${maximum} )); }; then
    printf 'invalid %s: expected at most %s, got %q\n' \
      "${label}" "${maximum}" "${value}" >&2
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

casebook_discovery_engine() {
  local mode="$1"
  "${PYTHON_BIN}" - "${mode}" 3<&0 <<'PY'
import json
import os
import pathlib
import re
import sys

mode = sys.argv[1]


class ScanError(ValueError):
    pass


def blank_non_newlines(output, start, end):
    for offset in range(start, end):
        if output[offset] not in (10, 13):
            output[offset] = 32


def raw_string_end(data, start):
    for prefix in (b"br", b"cr", b"r"):
        if not data.startswith(prefix, start):
            continue
        cursor = start + len(prefix)
        while cursor < len(data) and data[cursor] == 35:
            cursor += 1
        if cursor >= len(data) or data[cursor] != 34:
            continue
        marker = b"\"" + (b"#" * (cursor - start - len(prefix)))
        end = data.find(marker, cursor + 1)
        if end < 0:
            raise ScanError(f"unterminated raw string at byte {start}")
        return end + len(marker)
    return None


def quoted_end(data, quote_at, quote):
    cursor = quote_at + 1
    while cursor < len(data):
        byte = data[cursor]
        if byte == 92:
            cursor += 2
            continue
        if byte == quote:
            return cursor + 1
        cursor += 1
    raise ScanError(f"unterminated quoted literal at byte {quote_at}")


def char_literal_end(data, start):
    # Lifetimes and loop labels also begin with an apostrophe. A real character
    # literal with a one-byte identifier character closes immediately; otherwise
    # an ASCII identifier start denotes a lifetime rather than a literal.
    if start + 1 >= len(data):
        return None
    first = data[start + 1]
    if first == 95 or 65 <= first <= 90 or 97 <= first <= 122:
        if start + 2 < len(data) and data[start + 2] == 39:
            return start + 3
        return None
    # Escapes and UTF-8 scalar values are bounded conservatively here.
    cursor = start + 1
    limit = min(len(data), start + 32)
    while cursor < limit and data[cursor] not in (10, 13):
        if data[cursor] == 92:
            cursor += 2
            continue
        if data[cursor] == 39:
            return cursor + 1
        cursor += 1
    return None


def code_only(data):
    output = bytearray(data)
    cursor = 0
    while cursor < len(data):
        if data.startswith(b"//", cursor):
            end = data.find(b"\n", cursor + 2)
            end = len(data) if end < 0 else end
            blank_non_newlines(output, cursor, end)
            cursor = end
            continue
        if data.startswith(b"/*", cursor):
            start = cursor
            cursor += 2
            depth = 1
            while cursor < len(data) and depth:
                if data.startswith(b"/*", cursor):
                    depth += 1
                    cursor += 2
                elif data.startswith(b"*/", cursor):
                    depth -= 1
                    cursor += 2
                else:
                    cursor += 1
            if depth:
                raise ScanError(f"unterminated block comment at byte {start}")
            blank_non_newlines(output, start, cursor)
            continue
        raw_end = raw_string_end(data, cursor)
        if raw_end is not None:
            blank_non_newlines(output, cursor, raw_end)
            cursor = raw_end
            continue
        quote_at = None
        if data[cursor] == 34:
            quote_at = cursor
        elif data[cursor:cursor + 2] in (b"b\"", b"c\""):
            quote_at = cursor + 1
        if quote_at is not None:
            end = quoted_end(data, quote_at, 34)
            blank_non_newlines(output, cursor, end)
            cursor = end
            continue
        apostrophe_at = None
        if data[cursor] == 39:
            apostrophe_at = cursor
        elif (
            data[cursor] == 98
            and cursor + 1 < len(data)
            and data[cursor + 1] == 39
        ):
            apostrophe_at = cursor + 1
        if apostrophe_at is not None:
            end = char_literal_end(data, apostrophe_at)
            if end is not None:
                blank_non_newlines(output, cursor, end)
                cursor = end
                continue
        cursor += 1
    return bytes(output)


casebook_token = re.compile(
    rb"(?:"
    rb"\b(?:pub(?:\s*\([^)]*\))?\s+)?use\s+(?:::)?\s*fs_casebook\b|"
    rb"\buse\b[^;]*\bfs_casebook\b|"
    rb"\bextern\s+crate\s+fs_casebook\b|"
    rb"\bfs_casebook\s*::"
    rb")",
    re.DOTALL,
)


def matching_delimiter(data, opening_at):
    pairs = {40: 41, 91: 93, 123: 125}
    opening = data[opening_at]
    if opening not in pairs:
        raise ScanError(f"expected opening delimiter at byte {opening_at}")
    stack = [pairs[opening]]
    cursor = opening_at + 1
    while cursor < len(data):
        byte = data[cursor]
        if byte in pairs:
            stack.append(pairs[byte])
        elif byte in (41, 93, 125):
            if byte != stack[-1]:
                raise ScanError(f"mismatched delimiter at byte {cursor}")
            stack.pop()
            if not stack:
                return cursor
        cursor += 1
    raise ScanError(f"unterminated delimiter at byte {opening_at}")


def split_top_level_commas(data):
    parts = []
    start = 0
    cursor = 0
    while cursor < len(data):
        if data[cursor] in (40, 91, 123):
            cursor = matching_delimiter(data, cursor) + 1
            continue
        if data[cursor] == 44:
            parts.append(data[start:cursor])
            start = cursor + 1
        cursor += 1
    parts.append(data[start:])
    return parts


def attribute_path_and_tail(body):
    stripped = body.lstrip()
    match = re.match(rb"([A-Za-z_][A-Za-z0-9_]*)", stripped)
    if match is None:
        return None, b""
    return match.group(1), stripped[match.end():].lstrip()


def attribute_expands_to_ignore(body):
    path, tail = attribute_path_and_tail(body)
    if path == b"ignore":
        return not tail or tail[:1] in (b"=", b"(")
    if path != b"cfg_attr" or not tail.startswith(b"("):
        return False
    closing_at = matching_delimiter(tail, 0)
    if tail[closing_at + 1:].strip():
        return False
    arguments = split_top_level_commas(tail[1:closing_at])
    # cfg_attr's first argument is the cfg predicate. Only subsequent arguments
    # are attributes; a predicate named `ignore` must not create a false hit.
    return any(attribute_expands_to_ignore(argument) for argument in arguments[1:])


def contains_ignore_attribute(code):
    cursor = 0
    while cursor < len(code):
        marker = code.find(b"#", cursor)
        if marker < 0:
            return False
        opening_at = marker + 1
        while opening_at < len(code) and code[opening_at] in b" \t\r\n":
            opening_at += 1
        if opening_at < len(code) and code[opening_at] == 33:
            opening_at += 1
            while opening_at < len(code) and code[opening_at] in b" \t\r\n":
                opening_at += 1
        if opening_at >= len(code) or code[opening_at] != 91:
            cursor = marker + 1
            continue
        closing_at = matching_delimiter(code, opening_at)
        if attribute_expands_to_ignore(code[opening_at + 1:closing_at]):
            return True
        cursor = closing_at + 1
    return False


def classify(data):
    code = code_only(data)
    return (
        casebook_token.search(code) is not None,
        contains_ignore_attribute(code),
    )


if mode == "self-test":
    positive = [
        b"use fs_casebook::Casebook;",
        b"pub(crate) use ::fs_casebook as casebook;",
        b"use {fs_casebook as casebook};",
        b"extern crate fs_casebook;",
        b"fn f() { fs_casebook :: replay(); }",
        b"fn f<\x27a>(x: &\x27a str) { fs_casebook::replay(x); }",
    ]
    negative = [
        b"// use fs_casebook::Casebook;",
        b"/* fs_casebook::replay(); /* nested */ */ fn f() {}",
        b"const S: &str = \"fs_casebook::replay()\";",
        b"const S: &str = r###\"use fs_casebook::Casebook;\"###;",
        b"let my_fs_casebook = 1;",
        b"fn f<\x27fs_casebook>() {}",
    ]
    ignored = [
        b"use fs_casebook::Casebook; #[ignore] #[test] fn slow() {}",
        b"use fs_casebook::Casebook; #[cfg_attr(feature = \"slow\", ignore)] fn slow() {}",
        b"use fs_casebook::Casebook; #[cfg_attr(any(foo, bar), ignore = \"why\")] fn slow() {}",
        b"use fs_casebook::Casebook; #[cfg_attr(foo, cfg_attr(bar, ignore))] fn slow() {}",
    ]
    not_ignored = [
        b"use fs_casebook::Casebook; #[cfg_attr(ignore, test)] fn ordinary() {}",
        b"use fs_casebook::Casebook; #[cfg_attr(any(ignore, slow), allow(dead_code))] fn ordinary() {}",
    ]
    failures = []
    for index, data in enumerate(positive):
        if classify(data) != (True, False):
            failures.append(f"positive[{index}]")
    for index, data in enumerate(negative):
        if classify(data) != (False, False):
            failures.append(f"negative[{index}]")
    for index, data in enumerate(ignored):
        if classify(data) != (True, True):
            failures.append(f"ignored[{index}]")
    for index, data in enumerate(not_ignored):
        if classify(data) != (True, False):
            failures.append(f"not_ignored[{index}]")
    try:
        classify(b"use fs_casebook::Casebook; /*")
    except ScanError:
        pass
    else:
        failures.append("unterminated_comment")
    receipt = {
        "schema": "frankensim-casebook-source-scanner-self-test-v1",
        "status": "fail" if failures else "pass",
        "cases": len(positive) + len(negative) + len(ignored) + len(not_ignored) + 1,
        "failures": failures,
    }
    print(json.dumps(receipt, sort_keys=True, separators=(",", ":")))
    raise SystemExit(1 if failures else 0)

if mode != "discover":
    raise SystemExit(f"unsupported discovery-engine mode: {mode}")

try:
    meta = json.load(os.fdopen(3))
    root = pathlib.Path(meta["workspace_root"]).resolve(strict=True)
except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
    print(f"invalid locked Cargo metadata: {type(error).__name__}: {error}", file=sys.stderr)
    raise SystemExit(1)

rows = []
errors = []

for package in meta.get("packages", []):
    package_name = package.get("name", "")
    for target in package.get("targets", []):
        if "test" not in target.get("kind", []):
            continue
        target_name = target.get("name", "")
        try:
            source = pathlib.Path(target["src_path"]).resolve(strict=True)
        except (KeyError, OSError) as error:
            errors.append(
                f"unreadable target source: {package_name}:{target_name}: "
                f"{type(error).__name__}: {error}"
            )
            continue
        try:
            relative_source = source.relative_to(root).as_posix()
        except ValueError:
            errors.append(f"target source escapes workspace: {package_name}:{target_name}={source}")
            continue
        if any(character in relative_source for character in "\t\r\n"):
            errors.append(f"target source contains a row delimiter: {package_name}:{target_name}")
            continue
        try:
            has_casebook, has_ignore = classify(source.read_bytes())
        except (OSError, ScanError) as error:
            errors.append(
                f"cannot classify target source: {package_name}:{target_name}: "
                f"{type(error).__name__}: {error}"
            )
            continue
        if not has_casebook:
            continue
        identity = f"{package_name}:{target_name}"
        if has_ignore:
            errors.append(
                f"{identity} mixes #[ignore] with Casebook coverage; classify the ignored lane separately"
            )
            continue
        required_features = sorted(target.get("required-features") or [])
        if any(any(character in feature for character in "\t\r\n,") for feature in required_features):
            errors.append(f"{identity} has an unrepresentable required feature")
            continue
        features = ",".join(required_features) or "-"
        rows.append((identity, features, relative_source))

if errors:
    for error in sorted(errors):
        print(error, file=sys.stderr)
    raise SystemExit(1)

for identity, features, source in sorted(rows):
    print(f"discovered\t{identity}\t{features}\t{source}")
PY
}

discover_casebook_targets() {
  "${CARGO_BIN}" metadata --locked --format-version 1 --no-deps \
    | casebook_discovery_engine discover
}

# Input rows are tagged as pr/required/discovered. The validator is deliberately
# reusable by --self-test so the negative fixtures exercise the live policy.
validate_inventory_payload() {
  "${PYTHON_BIN}" -c '
import collections
import json
import pathlib
import re
import sys

valid_identity = re.compile(r"^[A-Za-z0-9_-]+:[A-Za-z0-9_-]+$")
valid_features = re.compile(r"^(?:-|[A-Za-z0-9_+./:-]+(?:,[A-Za-z0-9_+./:-]+)*)$")
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
        source_path = pathlib.PurePosixPath(source)
        if not valid_features.fullmatch(features):
            errors.append({"code": "invalid_features", "detail": identity})
        elif features != "-" and features.split(",") != sorted(set(features.split(","))):
            errors.append({"code": "noncanonical_features", "detail": identity})
        if (
            not source
            or source_path.is_absolute()
            or ".." in source_path.parts
            or any(character in source for character in "\r\n")
        ):
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
    "provenance_status": "unsealed",
    "provenance_scope": "locked-metadata-entrypoint-inventory-only",
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
  "${PYTHON_BIN}" - "${profile}" "${package}" "${target}" "${LOOKUP_FEATURES}" \
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
    "provenance_status": "unsealed",
    "provenance_scope": "selector-only",
}, sort_keys=True, separators=(",", ":")))
PY
}

list_profile() {
  local profile="$1"
  local budget target
  budget="$(profile_budget "${profile}")" || return 2
  printf '{"schema":"frankensim-casebook-profile-v1","profile":"%s","budget_seconds":%s,"build_time_included":true,"deadline_enforced":true,"deadline_clock":"monotonic","ignored_tests_included":false,"provenance_status":"unsealed","provenance_scope":"selector-only"}\n' \
    "${profile}" "${budget}"
  while IFS= read -r target; do
    emit_selector "${profile}" "${target}"
  done < <(profile_targets "${profile}")
}

monotonic_window() {
  "${PYTHON_BIN}" - "$1" <<'PY'
import sys
import time

started_ns = time.monotonic_ns()
budget_seconds = int(sys.argv[1])
print(f"{started_ns}\t{started_ns + budget_seconds * 1_000_000_000}")
PY
}

monotonic_elapsed_seconds() {
  "${PYTHON_BIN}" - "$1" <<'PY'
import sys
import time

started_ns = int(sys.argv[1])
print(max(0, (time.monotonic_ns() - started_ns) // 1_000_000_000))
PY
}

ACTIVE_WRAPPER_JOB=0
REQUESTED_SIGNAL=""

# shellcheck disable=SC2329 # Invoked indirectly by the signal traps below.
forward_active_signal() {
  local signal_name="$1"
  REQUESTED_SIGNAL="${signal_name}"
  if (( ACTIVE_WRAPPER_JOB != 0 )); then
    # Address the Bash job, not its numeric PID. Once wait has reaped the job,
    # `%+` ceases to resolve, so this cannot signal an unrelated PID that the OS
    # reused between wait returning and the active marker being cleared.
    kill -s "${signal_name}" %+ 2>/dev/null || true
  fi
}

trap 'forward_active_signal HUP' HUP
trap 'forward_active_signal INT' INT
trap 'forward_active_signal TERM' TERM

# Execute one selector in its own session/process group. Python enforces the
# monotonic aggregate deadline, TERM -> bounded wait -> KILL -> bounded drain,
# and emits the complete target receipt. The unreaped session leader pins the
# PGID while /bin/ps verifies live members, preventing PID reuse during drain.
run_target_until_deadline_with_grace() {
  local deadline_monotonic_ns="$1"
  local profile="$2"
  local identity="$3"
  local term_grace="$4"
  local kill_grace="$5"
  local wrapper_pid status
  shift 5
  "${PYTHON_BIN}" - "${deadline_monotonic_ns}" "${term_grace}" \
    "${kill_grace}" "${profile}" "${identity}" "$@" <<'PY' &
import json
import os
import signal
import subprocess
import sys
import time

deadline_monotonic_ns = int(sys.argv[1])
term_grace = float(sys.argv[2])
kill_grace = float(sys.argv[3])
profile = sys.argv[4]
identity = sys.argv[5]
command = sys.argv[6:]
package, target = identity.split(":", 1)
started_ns = time.monotonic_ns()
process = None
requested_signal = None
inspection_errors = []
observed_descendants = set()
term_signal_sent = False
kill_signal_sent = False

def emit_receipt(**fields):
    receipt = {
        "schema": "frankensim-casebook-profile-target-v1",
        "profile": profile,
        "package": package,
        "target": target,
        "command": command,
        "deadline_clock": "monotonic",
        "containment_scope": "new-session-process-group",
        "escaped_session_descendants_claimed": False,
        "child_output_channel": "stderr",
        "receipt_channel": "stdout-jsonl",
        "provenance_status": "unsealed",
        "provenance_scope": "selector-and-runtime-outcome-only",
    }
    receipt.update(fields)
    print(json.dumps(receipt, sort_keys=True, separators=(",", ":")), flush=True)

class RunnerInterrupted(Exception):
    def __init__(self, signum):
        super().__init__(signum)
        self.signum = signum

def remember_signal(signum, _frame):
    global requested_signal
    if requested_signal is None:
        requested_signal = signum

def check_interrupted():
    if requested_signal is not None:
        raise RunnerInterrupted(requested_signal)

for caught_signal in (signal.SIGHUP, signal.SIGINT, signal.SIGTERM):
    signal.signal(caught_signal, remember_signal)

if time.monotonic_ns() >= deadline_monotonic_ns:
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
        process_group_identity_pinned_until_drain=False,
        live_descendants_observed=0,
        inspection_errors=[],
        term_signal_sent=False,
        kill_signal_sent=False,
    )
    raise SystemExit(123)

try:
    if requested_signal is not None:
        signal_name = signal.Signals(requested_signal).name
        emit_receipt(
            status="interrupted",
            launched=False,
            exit_code=128 + requested_signal,
            leader_exit_code=None,
            elapsed_seconds=0,
            budget_status="within",
            drain_status="not_applicable",
            drain_trigger=f"signal_{signal_name}_before_spawn",
            drained_process_group=False,
            process_group_identity_pinned_until_drain=False,
            live_descendants_observed=0,
            inspection_errors=[],
            term_signal_sent=False,
            kill_signal_sent=False,
        )
        raise SystemExit(128 + requested_signal)
    process = subprocess.Popen(
        command,
        start_new_session=True,
        stdout=sys.stderr,
        stderr=sys.stderr,
    )
except OSError as error:
    emit_receipt(
        status="fail",
        launched=False,
        exit_code=126,
        leader_exit_code=None,
        elapsed_seconds=max(0, (time.monotonic_ns() - started_ns) // 1_000_000_000),
        budget_status="within",
        drain_status="not_applicable",
        drain_trigger="spawn_failure",
        drained_process_group=False,
        process_group_identity_pinned_until_drain=False,
        live_descendants_observed=0,
        inspection_errors=[],
        term_signal_sent=False,
        kill_signal_sent=False,
        spawn_error=f"{type(error).__name__}: {error}",
    )
    raise SystemExit(126)

def leader_status_without_reaping():
    info = os.waitid(
        os.P_PID,
        process.pid,
        os.WEXITED | os.WNOHANG | os.WNOWAIT,
    )
    if info is None:
        return None
    if info.si_code == os.CLD_EXITED:
        return int(info.si_status)
    return -int(info.si_status)

def live_group_descendants():
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    try:
        result = subprocess.run(
            ["/bin/ps", "-axo", "pid=,pgid=,stat="],
            check=False,
            capture_output=True,
            text=True,
            timeout=2,
            env=environment,
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise RuntimeError(f"process inspection failed: {type(error).__name__}: {error}")
    if result.returncode != 0:
        raise RuntimeError(
            f"process inspection exited {result.returncode}: {result.stderr.strip()}"
        )
    live = set()
    for line in result.stdout.splitlines():
        fields = line.split(None, 2)
        if len(fields) != 3:
            continue
        pid_text, pgid_text, state = fields
        try:
            pid = int(pid_text)
            pgid = int(pgid_text)
        except ValueError:
            continue
        try:
            sid = os.getsid(pid)
        except ProcessLookupError:
            continue
        except PermissionError as error:
            raise RuntimeError(
                f"cannot inspect session for pid {pid}: {type(error).__name__}: {error}"
            )
        if (
            pgid == process.pid
            and sid == process.pid
            and pid != process.pid
            and not state.startswith("Z")
        ):
            live.add(pid)
    observed_descendants.update(live)
    return live

def snapshot_live_descendants():
    try:
        return live_group_descendants()
    except RuntimeError as error:
        inspection_errors.append(str(error))
        return None

def wait_for_drain(deadline_ns):
    leader_exit_code = None
    while True:
        check_interrupted()
        leader_exit_code = leader_status_without_reaping()
        live = snapshot_live_descendants()
        if leader_exit_code is not None and live == set():
            return True, leader_exit_code
        now_ns = time.monotonic_ns()
        if now_ns >= deadline_ns:
            return False, leader_exit_code
        time.sleep(min(0.05, (deadline_ns - now_ns) / 1_000_000_000))

def signal_owned_group(signum):
    try:
        os.killpg(process.pid, signum)
    except ProcessLookupError:
        return False
    except PermissionError as error:
        inspection_errors.append(f"group signal denied: {type(error).__name__}: {error}")
        return False
    return True

def reap_if_exited(leader_exit_code):
    if leader_exit_code is None:
        return
    try:
        process.wait(timeout=0.2)
    except subprocess.TimeoutExpired:
        inspection_errors.append("leader was observable as exited but could not be reaped")

def drain_owned_group():
    global term_signal_sent
    global kill_signal_sent
    leader_exit_code = leader_status_without_reaping()
    live = snapshot_live_descendants()
    if leader_exit_code is not None and live == set():
        reap_if_exited(leader_exit_code)
        return "not_needed", leader_exit_code

    term_signal_sent = signal_owned_group(signal.SIGTERM) or term_signal_sent
    term_deadline_ns = time.monotonic_ns() + int(term_grace * 1_000_000_000)
    complete, observed_exit_code = wait_for_drain(term_deadline_ns)
    if not complete:
        kill_signal_sent = signal_owned_group(signal.SIGKILL) or kill_signal_sent
        kill_deadline_ns = time.monotonic_ns() + int(kill_grace * 1_000_000_000)
        complete, observed_exit_code = wait_for_drain(kill_deadline_ns)
    reap_if_exited(observed_exit_code)
    if inspection_errors:
        complete = False
    return ("complete" if complete else "incomplete"), observed_exit_code

def wait_for_leader_until(deadline_ns):
    while True:
        check_interrupted()
        leader_exit_code = leader_status_without_reaping()
        if leader_exit_code is not None:
            return leader_exit_code
        now_ns = time.monotonic_ns()
        if now_ns >= deadline_ns:
            return None
        time.sleep(min(0.05, (deadline_ns - now_ns) / 1_000_000_000))

try:
    leader_exit_code = wait_for_leader_until(deadline_monotonic_ns)
    if leader_exit_code is None:
        drain_status, leader_exit_code = drain_owned_group()
        exit_code = 124 if drain_status in ("complete", "not_needed") else 125
        status = "budget_exceeded"
        budget_status = "exceeded"
        drain_trigger = "deadline"
        wrapper_exit_code = exit_code
    else:
        drain_status, leader_exit_code = drain_owned_group()
        if drain_status == "not_needed":
            exit_code = leader_exit_code
            status = "pass" if leader_exit_code == 0 else "fail"
            budget_status = "within"
            drain_trigger = "none"
            wrapper_exit_code = 0 if leader_exit_code == 0 else 1
        else:
            exit_code = 1 if drain_status == "complete" else 127
            status = "fail"
            budget_status = "within"
            drain_trigger = "leader_exit_with_live_group"
            wrapper_exit_code = exit_code
except RunnerInterrupted as interruption:
    for caught_signal in (signal.SIGHUP, signal.SIGINT, signal.SIGTERM):
        signal.signal(caught_signal, signal.SIG_IGN)
    requested_signal = None
    drain_status, leader_exit_code = drain_owned_group()
    signal_name = signal.Signals(interruption.signum).name
    exit_code = 128 + interruption.signum
    status = "interrupted"
    budget_status = "within"
    drain_trigger = f"signal_{signal_name}"
    # Preserve the conventional signal status only after bounded cleanup
    # succeeds. A failed drain is a containment failure, not a clean interrupt.
    wrapper_exit_code = exit_code if drain_status in ("complete", "not_needed") else 127

emit_receipt(
    status=status,
    launched=True,
    exit_code=exit_code,
    leader_exit_code=leader_exit_code,
    elapsed_seconds=max(0, (time.monotonic_ns() - started_ns) // 1_000_000_000),
    budget_status=budget_status,
    drain_status=drain_status,
    drain_trigger=drain_trigger,
    drained_process_group=drain_status == "complete",
    process_group_identity_pinned_until_drain=True,
    live_descendants_observed=len(observed_descendants),
    inspection_errors=sorted(set(inspection_errors)),
    term_signal_sent=term_signal_sent,
    kill_signal_sent=kill_signal_sent,
)
raise SystemExit(wrapper_exit_code)
PY
  wrapper_pid=$!
  ACTIVE_WRAPPER_JOB=1
  if [[ -n "${REQUESTED_SIGNAL}" ]]; then
    kill -s "${REQUESTED_SIGNAL}" %+ 2>/dev/null || true
  fi
  while :; do
    if wait "${wrapper_pid}"; then
      status=0
      break
    else
      status=$?
    fi
    # A trapped shell signal interrupts wait before the child is reaped. Only
    # retry in that case and only while Bash still reports this exact child as a
    # running job. Never probe the numeric PID after a nonzero child exit: it may
    # already have been reaped and reused.
    if [[ -n "${REQUESTED_SIGNAL}" ]] \
        && [[ $'\n'"$(jobs -pr)"$'\n' == *$'\n'"${wrapper_pid}"$'\n'* ]]; then
      continue
    fi
    break
  done
  ACTIVE_WRAPPER_JOB=0
  return "${status}"
}

run_target_until_deadline() {
  local deadline_monotonic_ns="$1"
  local profile="$2"
  local identity="$3"
  shift 3
  run_target_until_deadline_with_grace \
    "${deadline_monotonic_ns}" "${profile}" "${identity}" \
    "${TERMINATION_GRACE_SECONDS}" "${KILL_DRAIN_SECONDS}" "$@"
}

run_profile() {
  local profile="$1"
  local budget="$2"
  local window started_ns deadline_ns identity package target targets elapsed_seconds
  local status overall_status=0 budget_status="within" run_status="pass"
  local interrupted_status=0
  local terminal_exit_code=0
  local total=0 attempted=0 passed=0 failed=0 budget_exceeded=0 unattempted=0
  local -a command=()
  targets="$(profile_targets "${profile}")"
  while IFS= read -r identity; do
    [[ -n "${identity}" ]] && total=$((total + 1))
  done <<< "${targets}"
  window="$(monotonic_window "${budget}")"
  IFS=$'\t' read -r started_ns deadline_ns <<< "${window}"
  printf '{"schema":"frankensim-casebook-profile-run-v1","event":"start","profile":"%s","budget_seconds":%s,"deadline_clock":"monotonic","build_time_included":true,"deadline_enforced":true,"ignored_tests_included":false,"total_targets":%s,"child_output_channel":"stderr","receipt_channel":"stdout-jsonl","containment_scope":"new-session-process-group","escaped_session_descendants_claimed":false,"provenance_status":"unsealed","provenance_scope":"selector-and-runtime-outcome-only"}\n' \
    "${profile}" "${budget}" "${total}"

  while IFS= read -r identity; do
    [[ -z "${identity}" ]] && continue
    if [[ -n "${REQUESTED_SIGNAL}" ]]; then
      case "${REQUESTED_SIGNAL}" in
        HUP) interrupted_status=129 ;;
        INT) interrupted_status=130 ;;
        TERM) interrupted_status=143 ;;
        *) interrupted_status=1 ;;
      esac
      overall_status=1
      run_status="interrupted"
      break
    fi
    package="${identity%%:*}"
    target="${identity#*:}"
    lookup_discovery "${identity}"
    command=("${CARGO_BIN}" test --locked -p "${package}" --test "${target}")
    if [[ -n "${LOOKUP_FEATURES}" ]]; then
      command+=(--features "${LOOKUP_FEATURES}")
    fi
    command+=(-- --nocapture)
    attempted=$((attempted + 1))
    if run_target_until_deadline "${deadline_ns}" "${profile}" "${identity}" "${command[@]}"; then
      status=0
      passed=$((passed + 1))
    else
      status=$?
      overall_status=1
      case "${status}" in
        1)
          failed=$((failed + 1))
          ;;
        123|124|125)
          budget_exceeded=$((budget_exceeded + 1))
          budget_status="exceeded"
          break
          ;;
        126|127)
          failed=$((failed + 1))
          break
          ;;
        129|130|143)
          failed=$((failed + 1))
          interrupted_status="${status}"
          run_status="interrupted"
          break
          ;;
        *)
          failed=$((failed + 1))
          break
          ;;
      esac
    fi
  done <<< "${targets}"

  if [[ -n "${REQUESTED_SIGNAL}" ]] && (( interrupted_status == 0 )); then
    case "${REQUESTED_SIGNAL}" in
      HUP) interrupted_status=129 ;;
      INT) interrupted_status=130 ;;
      TERM) interrupted_status=143 ;;
      *) interrupted_status=1 ;;
    esac
    overall_status=1
  fi
  elapsed_seconds="$(monotonic_elapsed_seconds "${started_ns}")"
  unattempted=$((total - attempted))
  if (( interrupted_status != 0 )); then
    run_status="interrupted"
  elif [[ "${budget_status}" == "exceeded" ]]; then
    run_status="budget_exceeded"
  elif (( overall_status != 0 )); then
    run_status="fail"
  fi
  if (( interrupted_status != 0 )); then
    terminal_exit_code="${interrupted_status}"
  elif (( overall_status != 0 )); then
    terminal_exit_code=1
  fi
  printf '{"schema":"frankensim-casebook-profile-run-v1","event":"finish","profile":"%s","status":"%s","terminal_exit_code":%s,"budget_status":"%s","budget_seconds":%s,"elapsed_seconds":%s,"total_targets":%s,"attempted_targets":%s,"expected_target_receipts":%s,"passed_targets":%s,"failed_targets":%s,"budget_exceeded_targets":%s,"unattempted_targets":%s,"deadline_clock":"monotonic","child_output_channel":"stderr","receipt_channel":"stdout-jsonl","provenance_status":"unsealed","provenance_scope":"selector-and-runtime-outcome-only"}\n' \
    "${profile}" "${run_status}" "${terminal_exit_code}" "${budget_status}" "${budget}" \
    "${elapsed_seconds}" "${total}" "${attempted}" "${attempted}" "${passed}" \
    "${failed}" "${budget_exceeded}" "${unattempted}"
  if (( interrupted_status != 0 )); then
    return "${interrupted_status}"
  fi
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

expect_validation_success() {
  local payload="$1"
  local output
  if ! output="$(printf '%s\n' "${payload}" | validate_inventory_payload)"; then
    printf 'self-test expected valid inventory, observed: %s\n' "${output}" >&2
    return 1
  fi
  if [[ "${output}" != *'"status":"pass"'* ]]; then
    printf 'self-test valid inventory receipt mismatch: %s\n' "${output}" >&2
    return 1
  fi
}

self_test_require_contains() {
  local output="$1"
  local expected="$2"
  local label="$3"
  if [[ "${output}" != *"${expected}"* ]]; then
    printf 'self-test %s missing %s: %s\n' "${label}" "${expected}" "${output}" >&2
    return 1
  fi
}

self_test_require_single_json_line() {
  local output="$1"
  local label="$2"
  if [[ -z "${output}" || "${output}" == *$'\n'* \
      || "${output}" != \{* || "${output}" != *\} ]]; then
    printf 'self-test %s did not emit exactly one JSONL receipt: %s\n' \
      "${label}" "${output}" >&2
    return 1
  fi
  if ! printf '%s\n' "${output}" | "${PYTHON_BIN}" -c \
      'import json, sys; json.load(sys.stdin)'; then
    printf 'self-test %s emitted invalid JSON: %s\n' "${label}" "${output}" >&2
    return 1
  fi
}

self_test_deadline_ns() {
  "${PYTHON_BIN}" - "$1" <<'PY'
import sys
import time

print(time.monotonic_ns() + int(sys.argv[1]) * 1_000_000)
PY
}

replay_profile_receipt() {
  local receipt_path="$1"
  local pr_list full_list
  pr_list="$(profile_targets pr)"
  full_list="$(profile_targets nightly-full)"
  "${PYTHON_BIN}" - "${receipt_path}" "${pr_list}" "${full_list}" 3<&0 <<'PY'
import json
import os
import sys

receipt_path, pr_raw, full_raw = sys.argv[1:4]
expected_by_profile = {
    "pr": [line for line in pr_raw.splitlines() if line],
    "nightly-full": [line for line in full_raw.splitlines() if line],
}
errors = []

try:
    if receipt_path == "-":
        with os.fdopen(3, "r", encoding="utf-8") as handle:
            raw_lines = handle.read().splitlines()
    else:
        with open(receipt_path, "r", encoding="utf-8") as handle:
            raw_lines = handle.read().splitlines()
except OSError as error:
    out = {
        "schema": "frankensim-casebook-profile-replay-v1",
        "status": "fail",
        "profile": None,
        "verified_targets": 0,
        "errors": [{"code": "unreadable_receipt", "detail": f"{type(error).__name__}: {error}"}],
    }
    print(json.dumps(out, sort_keys=True, separators=(",", ":")))
    raise SystemExit(1)

records = []
for index, line in enumerate(raw_lines, 1):
    if not line.strip():
        continue
    try:
        obj = json.loads(line)
    except json.JSONDecodeError as error:
        errors.append({"code": "invalid_json", "detail": f"line {index}: {error}"})
        continue
    if not isinstance(obj, dict):
        errors.append({"code": "non_object_row", "detail": f"line {index}"})
        continue
    records.append(obj)

if records and records[0].get("schema") == "frankensim-casebook-profile-inventory-v1":
    inv = records.pop(0)
    if inv.get("status") != "pass":
        errors.append({"code": "inventory_failed", "detail": str(inv.get("status"))})

start_row = None
finish_row = None
target_rows = []
for index, row in enumerate(records):
    schema = row.get("schema")
    if schema == "frankensim-casebook-profile-run-v1" and row.get("event") == "start":
        if start_row is not None or index != 0:
            errors.append({"code": "misplaced_start", "detail": f"record {index}"})
        start_row = row
    elif schema == "frankensim-casebook-profile-run-v1" and row.get("event") == "finish":
        if finish_row is not None or index != len(records) - 1:
            errors.append({"code": "misplaced_finish", "detail": f"record {index}"})
        finish_row = row
    elif schema == "frankensim-casebook-profile-target-v1":
        target_rows.append(row)
    else:
        errors.append({"code": "unexpected_schema", "detail": str(schema)})

profile = start_row.get("profile") if start_row else None
if start_row is None:
    errors.append({"code": "missing_start", "detail": "no start event found"})
if finish_row is None:
    errors.append({"code": "missing_finish", "detail": "no finish event found"})
if profile not in expected_by_profile:
    errors.append({"code": "invalid_profile", "detail": str(profile)})
else:
    expected_targets = expected_by_profile[profile]
    observed_targets = [
        f"{row.get('package', '')}:{row.get('target', '')}" for row in target_rows
    ]
    if observed_targets != expected_targets:
        errors.append({
            "code": "target_sequence_mismatch",
            "detail": f"expected {len(expected_targets)} targets, observed {len(observed_targets)}",
        })
    for row in target_rows:
        identity = f"{row.get('package', '')}:{row.get('target', '')}"
        if (
            row.get("profile") != profile
            or row.get("status") != "pass"
            or row.get("launched") is not True
            or row.get("exit_code") != 0
            or row.get("budget_status") != "within"
            or row.get("drain_status") != "not_needed"
            or row.get("escaped_session_descendants_claimed") is not False
            or row.get("inspection_errors") != []
        ):
            errors.append({"code": "target_failed", "detail": identity})
    if finish_row is not None:
        total_expected = len(expected_targets)
        if (
            finish_row.get("profile") != profile
            or finish_row.get("status") != "pass"
            or finish_row.get("terminal_exit_code") != 0
            or finish_row.get("budget_status") != "within"
            or finish_row.get("total_targets") != total_expected
            or finish_row.get("attempted_targets") != total_expected
            or finish_row.get("passed_targets") != total_expected
            or finish_row.get("failed_targets") != 0
            or finish_row.get("budget_exceeded_targets") != 0
            or finish_row.get("unattempted_targets") != 0
        ):
            errors.append({"code": "finish_mismatch", "detail": str(finish_row.get("status"))})

receipt = {
    "schema": "frankensim-casebook-profile-replay-v1",
    "status": "fail" if errors else "pass",
    "profile": profile,
    "verified_targets": len(target_rows) if not errors else 0,
    "errors": errors,
}
print(json.dumps(receipt, sort_keys=True, separators=(",", ":")))
raise SystemExit(1 if errors else 0)
PY
}

run_self_tests() {
  local base missing stale duplicate_required duplicate_pr duplicate_discovered
  local pr_not_required malformed scanner selector saved_discovery
  local deadline output status
  base=$'pr\tfs-a:casebook\nrequired\tfs-a:casebook\nrequired\tfs-b:conformance\ndiscovered\tfs-a:casebook\t-\tcrates/fs-a/tests/casebook.rs\ndiscovered\tfs-b:conformance\t-\tcrates/fs-b/tests/conformance.rs'
  missing="${base}"$'\npr\tfs-c:missing_case\nrequired\tfs-c:missing_case'
  stale="${base}"$'\nrequired\tfs-c:removed_case'
  duplicate_required="${base}"$'\nrequired\tfs-b:conformance'
  duplicate_pr="${base}"$'\npr\tfs-a:casebook'
  duplicate_discovered="${base}"$'\ndiscovered\tfs-b:conformance\t-\tcrates/fs-b/tests/conformance.rs'
  pr_not_required="${base}"$'\npr\tfs-c:casebook\ndiscovered\tfs-c:casebook\t-\tcrates/fs-c/tests/casebook.rs'
  malformed="${base}"$'\ndiscovered\tfs-c:casebook\tfeature-b,feature-a\t../escape.rs'
  expect_validation_success "${base}"
  expect_validation_failure "missing_pr" "${missing}"
  expect_validation_failure "stale_required" "${stale}"
  expect_validation_failure "duplicate_required" "${duplicate_required}"
  expect_validation_failure "duplicate_pr" "${duplicate_pr}"
  expect_validation_failure "duplicate_discovered" "${duplicate_discovered}"
  expect_validation_failure "pr_not_required" "${pr_not_required}"
  expect_validation_failure "noncanonical_features" "${malformed}"
  expect_validation_failure "malformed_discovered" "${malformed}"

  scanner="$(casebook_discovery_engine self-test)"
  self_test_require_single_json_line "${scanner}" scanner
  self_test_require_contains "${scanner}" '"status":"pass"' scanner

  saved_discovery="${DISCOVERY_ROWS}"
  DISCOVERY_ROWS=$'discovered\tfs-a:casebook\tfeature-a,feature-b\tcrates/fs-a/tests/casebook.rs'
  selector="$(emit_selector pr fs-a:casebook)"
  DISCOVERY_ROWS="${saved_discovery}"
  self_test_require_single_json_line "${selector}" selector
  self_test_require_contains "${selector}" \
    '"required_features":["feature-a","feature-b"]' selector
  self_test_require_contains "${selector}" \
    '"command":[' selector

  deadline="$(self_test_deadline_ns -1)"
  if output="$(run_target_until_deadline_with_grace \
      "${deadline}" pr fs-a:casebook 0.05 0.05 \
      /definitely-not-a-casebook-self-test-command 2>/dev/null)"; then
    printf 'self-test deadline refusal unexpectedly allowed launch\n' >&2
    return 1
  else
    status=$?
  fi
  if (( status != 123 )); then
    printf 'self-test deadline refusal returned %s: %s\n' "${status}" "${output}" >&2
    return 1
  fi
  self_test_require_single_json_line "${output}" deadline-refusal
  self_test_require_contains "${output}" '"status":"budget_exceeded"' deadline-refusal
  self_test_require_contains "${output}" '"launched":false' deadline-refusal
  self_test_require_contains "${output}" \
    '"drain_status":"not_applicable"' deadline-refusal

  deadline="$(self_test_deadline_ns 5000)"
  if output="$(run_target_until_deadline_with_grace \
      "${deadline}" pr fs-a:casebook 0.05 0.05 \
      /definitely-not-a-casebook-self-test-command 2>/dev/null)"; then
    printf 'self-test spawn failure unexpectedly succeeded\n' >&2
    return 1
  else
    status=$?
  fi
  if (( status != 126 )); then
    printf 'self-test spawn failure returned %s: %s\n' "${status}" "${output}" >&2
    return 1
  fi
  self_test_require_single_json_line "${output}" spawn-failure
  self_test_require_contains "${output}" '"drain_trigger":"spawn_failure"' spawn-failure

  deadline="$(self_test_deadline_ns 5000)"
  if output="$(run_target_until_deadline_with_grace \
      "${deadline}" pr fs-a:casebook 0.05 0.05 \
      "${PYTHON_BIN}" -c 'print("CHILD_STDOUT_MUST_NOT_ENTER_RECEIPTS")' \
      2>/dev/null)"; then
    status=0
  else
    status=$?
  fi
  if (( status != 0 )); then
    printf 'self-test stdout isolation failed (%s): %s\n' "${status}" "${output}" >&2
    return 1
  fi
  self_test_require_single_json_line "${output}" stdout-isolation
  self_test_require_contains "${output}" '"status":"pass"' stdout-isolation

  deadline="$(self_test_deadline_ns 5000)"
  if output="$(run_target_until_deadline_with_grace \
      "${deadline}" pr fs-a:casebook 0.05 0.05 \
      "${PYTHON_BIN}" -c 'raise SystemExit(7)' 2>/dev/null)"; then
    status=0
  else
    status=$?
  fi
  if (( status != 1 )); then
    printf 'self-test ordinary failure returned %s: %s\n' "${status}" "${output}" >&2
    return 1
  fi
  self_test_require_single_json_line "${output}" ordinary-failure
  self_test_require_contains "${output}" '"leader_exit_code":7' ordinary-failure
  self_test_require_contains "${output}" '"drain_status":"not_needed"' ordinary-failure

  deadline="$(self_test_deadline_ns 5000)"
  if output="$(run_target_until_deadline_with_grace \
      "${deadline}" pr fs-a:casebook 0.10 0.50 \
      "${PYTHON_BIN}" -c \
      'import subprocess,sys,time; subprocess.Popen([sys.executable,"-c","import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(30)"]); time.sleep(0.2)' \
      2>/dev/null)"; then
    status=0
  else
    status=$?
  fi
  if (( status != 1 )); then
    printf 'self-test lingering descendant returned %s: %s\n' "${status}" "${output}" >&2
    return 1
  fi
  self_test_require_single_json_line "${output}" descendant-drain
  self_test_require_contains "${output}" \
    '"drain_trigger":"leader_exit_with_live_group"' descendant-drain
  self_test_require_contains "${output}" '"drain_status":"complete"' descendant-drain
  self_test_require_contains "${output}" '"term_signal_sent":true' descendant-drain
  self_test_require_contains "${output}" '"kill_signal_sent":true' descendant-drain

  deadline="$(self_test_deadline_ns 1000)"
  if output="$(run_target_until_deadline_with_grace \
      "${deadline}" pr fs-a:casebook 0.10 0.50 \
      "${PYTHON_BIN}" -c \
      'import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(30)' \
      2>/dev/null)"; then
    status=0
  else
    status=$?
  fi
  if (( status != 124 )); then
    printf 'self-test timeout returned %s: %s\n' "${status}" "${output}" >&2
    return 1
  fi
  self_test_require_single_json_line "${output}" timeout-drain
  self_test_require_contains "${output}" '"status":"budget_exceeded"' timeout-drain
  self_test_require_contains "${output}" '"drain_status":"complete"' timeout-drain
  self_test_require_contains "${output}" '"term_signal_sent":true' timeout-drain
  self_test_require_contains "${output}" '"kill_signal_sent":true' timeout-drain

  deadline="$(self_test_deadline_ns 5000)"
  if output="$(run_target_until_deadline_with_grace \
      "${deadline}" pr fs-a:casebook 0.10 0.50 \
      "${PYTHON_BIN}" -c \
      'import os,signal,time; os.kill(os.getppid(), signal.SIGTERM); time.sleep(30)' \
      2>/dev/null)"; then
    status=0
  else
    status=$?
  fi
  if (( status != 143 )); then
    printf 'self-test interrupt cleanup returned %s: %s\n' "${status}" "${output}" >&2
    return 1
  fi
  self_test_require_single_json_line "${output}" interrupt-drain
  self_test_require_contains "${output}" '"status":"interrupted"' interrupt-drain
  self_test_require_contains "${output}" '"drain_status":"complete"' interrupt-drain

  if [[ -n "${REQUESTED_SIGNAL}" ]]; then
    printf 'self-test unexpectedly signalled the outer shell: %s\n' \
      "${REQUESTED_SIGNAL}" >&2
    return 1
  fi

  if [[ "${output}" != *'"escaped_session_descendants_claimed":false'* ]]; then
    printf 'self-test no-claim boundary missing from interrupt receipt: %s\n' \
      "${output}" >&2
    return 1
  fi

  # Exercise the shell trap itself, not only the wrapper's Python signal
  # handler. The synthetic child finds the wrapper's parent shell and signals
  # that shell; the trap must address the active Bash job, the wrapper must
  # drain its owned child group, and the shell must survive with status 143.
  REQUESTED_SIGNAL=""
  deadline="$(self_test_deadline_ns 5000)"
  if run_target_until_deadline_with_grace \
      "${deadline}" pr fs-a:casebook 0.10 0.50 \
      "${PYTHON_BIN}" -c \
      'import os,signal,subprocess,time; wrapper=os.getppid(); shell=int(subprocess.check_output(["/bin/ps","-o","ppid=","-p",str(wrapper)]).strip()); os.kill(shell, signal.SIGTERM); time.sleep(30)' \
      >/dev/null 2>/dev/null; then
    status=0
  else
    status=$?
  fi
  if (( status != 143 )) || [[ "${REQUESTED_SIGNAL}" != "TERM" ]]; then
    printf 'self-test shell signal forwarding returned %s with request %q\n' \
      "${status}" "${REQUESTED_SIGNAL}" >&2
    return 1
  fi
  REQUESTED_SIGNAL=""

  # Exercise deterministic replay verification over synthetic run streams.
  saved_discovery="${DISCOVERY_ROWS}"
  DISCOVERY_ROWS="$(
    for target in "${PR_TARGETS[@]}"; do
      printf 'discovered\t%s\t-\tcrates/%s/tests/%s.rs\n' \
        "${target}" "${target%%:*}" "${target#*:}"
    done
  )"
  local replay_stream
  replay_stream="$(
    printf '{"schema":"frankensim-casebook-profile-run-v1","event":"start","profile":"pr","total_targets":%s}\n' "${#PR_TARGETS[@]}"
    for target in "${PR_TARGETS[@]}"; do
      printf '{"schema":"frankensim-casebook-profile-target-v1","profile":"pr","package":"%s","target":"%s","status":"pass","launched":true,"exit_code":0,"budget_status":"within","drain_status":"not_needed","escaped_session_descendants_claimed":false,"inspection_errors":[]}\n' \
        "${target%%:*}" "${target#*:}"
    done
    printf '{"schema":"frankensim-casebook-profile-run-v1","event":"finish","profile":"pr","status":"pass","terminal_exit_code":0,"budget_status":"within","total_targets":%s,"attempted_targets":%s,"passed_targets":%s,"failed_targets":0,"budget_exceeded_targets":0,"unattempted_targets":0}\n' \
      "${#PR_TARGETS[@]}" "${#PR_TARGETS[@]}" "${#PR_TARGETS[@]}"
  )"
  output="$(printf '%s\n' "${replay_stream}" | replay_profile_receipt -)"
  self_test_require_single_json_line "${output}" replay-pass
  self_test_require_contains "${output}" '"status":"pass"' replay-pass
  if output="$(printf '%s\n' "${replay_stream/\"status\":\"pass\"/\"status\":\"fail\"}" | replay_profile_receipt -)"; then
    printf 'self-test tampered replay unexpectedly passed: %s\n' "${output}" >&2
    DISCOVERY_ROWS="${saved_discovery}"
    return 1
  fi
  self_test_require_single_json_line "${output}" replay-tampered
  self_test_require_contains "${output}" '"status":"fail"' replay-tampered
  DISCOVERY_ROWS="${saved_discovery}"

  printf '%s\n' '{"schema":"frankensim-casebook-profile-self-test-v2","status":"pass","cases":21,"scanner_fixture_cases":19,"cargo_invocations":0,"temporary_files":0}'
}

if (( $# == 1 )) && [[ "$1" == "--self-test" ]]; then
  run_self_tests
  exit $?
fi

require_bounded_positive_integer \
  "PR budget" "${PR_BUDGET_SECONDS}" "${MAX_PROFILE_BUDGET_SECONDS}"
require_bounded_positive_integer \
  "full budget" "${FULL_BUDGET_SECONDS}" "${MAX_PROFILE_BUDGET_SECONDS}"
require_bounded_positive_integer \
  "termination grace" "${TERMINATION_GRACE_SECONDS}" "${MAX_DRAIN_SECONDS}"
require_bounded_positive_integer \
  "kill drain" "${KILL_DRAIN_SECONDS}" "${MAX_DRAIN_SECONDS}"
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

if (( $# == 2 )) && [[ "$1" == "--replay" ]]; then
  refresh_and_validate_inventory >/dev/null
  replay_profile_receipt "$2"
  exit $?
fi

PROFILE_ARG=""
if (( $# == 1 )) && [[ "$1" == "pr" || "$1" == "nightly-full" ]]; then
  PROFILE_ARG="$1"
elif (( $# == 2 )) && [[ "$1" == "--run" ]] && [[ "$2" == "pr" || "$2" == "nightly-full" ]]; then
  PROFILE_ARG="$2"
else
  usage
  exit 2
fi

readonly PROFILE="${PROFILE_ARG}"
BUDGET_SECONDS="$(profile_budget "${PROFILE}")"
readonly BUDGET_SECONDS
refresh_and_validate_inventory
run_profile "${PROFILE}" "${BUDGET_SECONDS}"
exit $?
