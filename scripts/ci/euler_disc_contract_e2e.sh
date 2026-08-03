#!/usr/bin/env bash
# Focused and clean-HEAD-bookended closure-candidate proof for the Euler-disc
# scientific-contract leaf. This runner validates software structure only: it
# never upgrades a synthetic protocol fixture into physical validation,
# mechanism evidence, or a prediction about the Steve Mould observations.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

PROFILE="focused"
VERIFY_BUNDLE=""
VERIFY_BUNDLE_SET=0
SELF_TEST=0
RETAINED_LOG_CHECKER_SMOKE_COMMAND='cargo test --locked -p fs-euler-disc-e2e --test scientific_contract -- g0_check_receipts_logs_and_domain_separation_are_exact_and_bounded --exact --test-threads=1'

usage() {
  printf '%s\n' \
    'usage: scripts/ci/euler_disc_contract_e2e.sh [--profile focused|closure]' \
    '       scripts/ci/euler_disc_contract_e2e.sh --verify-bundle DIR' \
    '       scripts/ci/euler_disc_contract_e2e.sh --self-test' \
    '' \
    'focused  Runs the crate and structural registry gates. Because this may' \
    '         run on untracked or dirty work, source-manifest authority is' \
    '         recorded as NO_DATA rather than inferred.' \
    'closure  Additionally requires the complete repository tree to be tracked' \
    '         and clean, runs docs/source-manifest gates, and independently' \
    '         verifies source-manifest membership. This is a clean-HEAD-observed' \
    '         candidate for DSR, not proof against transient concurrent edits.' \
    '' \
    '--verify-bundle verifies the exact sealed inventory, ordered terminal' \
    'trace, verdict-prefix bytes/hash, and every lane-log size/hash without' \
    'running Cargo.' \
    '--self-test exercises supervisor, drain, cap, sentinel, and hostile bundle' \
    'verification failure paths without running Cargo or deleting evidence.' \
    '' \
    'Execution environment:' \
    '  FSIM_EULER_DISC_E2E_EXECUTOR=local|rch|dsr (caller declaration only;' \
    '    default: local; this script does not itself route through RCH or DSR)' \
    '  FSIM_EULER_DISC_E2E_ALLOW_LOCAL=1 (required for local Cargo)' \
    '  CARGO_TARGET_DIR=... (required when executor is declared rch)' \
    '  FSIM_EULER_DISC_E2E_CARGO=/path/to/cargo' \
    '  FSIM_EULER_DISC_E2E_LOG_DIR=target/euler-disc-contract-e2e' \
    '  FSIM_EULER_DISC_E2E_LANE_TIMEOUT_SECONDS=3600' \
    '  FSIM_EULER_DISC_E2E_RUN_TIMEOUT_SECONDS=14400' \
    '  FSIM_EULER_DISC_E2E_LANE_LOG_MAX_BYTES=16777216'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --profile)
      if [[ $# -lt 2 ]]; then
        printf '%s\n' '--profile requires focused or closure' >&2
        usage >&2
        exit 2
      fi
      PROFILE="${2:-}"
      shift 2
      ;;
    --verify-bundle)
      if [[ $# -lt 2 ]]; then
        printf '%s\n' '--verify-bundle requires a proof-bundle directory' >&2
        usage >&2
        exit 2
      fi
      VERIFY_BUNDLE="${2:-}"
      VERIFY_BUNDLE_SET=1
      if [[ -z "$VERIFY_BUNDLE" ]]; then
        printf '%s\n' '--verify-bundle requires a nonempty proof-bundle directory' >&2
        usage >&2
        exit 2
      fi
      shift 2
      ;;
    --self-test)
      SELF_TEST=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'unknown argument: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

case "$PROFILE" in
  focused|closure) ;;
  *)
    printf 'unknown profile: %s\n' "$PROFILE" >&2
    exit 2
    ;;
esac
if [[ "$SELF_TEST" == 1 && "$PROFILE" != "focused" ]]; then
  printf '%s\n' '--self-test requires the focused profile' >&2
  exit 2
fi

# Producer and verifier resource ceilings are one contract. Keep these values
# in shell so normal execution, no-Cargo self-tests, and offline verification
# cannot silently drift onto different budgets.
PROOF_MAX_SUMMARY_BYTES=$((64 * 1024))
PROOF_MAX_VERDICTS_BYTES=$((2 * 1024 * 1024))
PROOF_MAX_RECORD_BYTES=$((64 * 1024))
PROOF_MAX_RECORDS=96
PROOF_MAX_LANE_LOG_BYTES=$((64 * 1024 * 1024))
PROOF_MAX_TOTAL_LOG_BYTES=$((512 * 1024 * 1024))
PROOF_MAX_SNAPSHOT_BYTES=$((16 * 1024 * 1024))
PROOF_MAX_BUNDLE_ENTRIES=$((PROOF_MAX_RECORDS * 2 + 16))
PROOF_LOG_METADATA_RESERVE_BYTES=$((32 * 1024))
PROOF_TERMINAL_LOG_RESERVE_BYTES=$((256 * 1024))
SUCCESS_FINALIZATION_LOG_MAX_BYTES=$((1024 * 1024))
SUCCESS_FINALIZATION_TOTAL_LOG_MAX_BYTES=$((
  SUCCESS_FINALIZATION_LOG_MAX_BYTES + PROOF_LOG_METADATA_RESERVE_BYTES + 64 * 1024
))
PROOF_OPERATIONAL_LOG_BUDGET=$((
  PROOF_MAX_TOTAL_LOG_BYTES - PROOF_TERMINAL_LOG_RESERVE_BYTES
))
PROOF_MAX_TIMEOUT_SECONDS=$((7 * 24 * 60 * 60))
# The hostile normal-harness fixtures launch a complete focused run inside an
# already supervised self-test lane. Give that nested process tree a distinct,
# finite containment budget so transient shared-host load cannot masquerade as
# the semantic failure each fixture is meant to observe. The recursive
# candidate-corruption case executes all four fixtures and therefore retains a
# separately bounded aggregate allowance with one fixture's worth of headroom.
NORMAL_HARNESS_HOSTILE_SELF_TEST_TIMEOUT_SECONDS=120
NORMAL_HARNESS_NESTED_LANE_TIMEOUT_SECONDS=20
NORMAL_HARNESS_NESTED_RUN_TIMEOUT_SECONDS=90
RECURSIVE_SELF_TEST_TIMEOUT_SECONDS=$((
  NORMAL_HARNESS_HOSTILE_SELF_TEST_TIMEOUT_SECONDS * 5
))
FINALIZER_HANDSHAKE_TIMEOUT_SECONDS=30
FINALIZER_ACK_TIMEOUT_SECONDS=10
WRAPPER_SIGNAL_SELF_TEST_TIMEOUT_SECONDS=120
WRAPPER_SIGNAL_READINESS_TIMEOUT_SECONDS=30
WRAPPER_SIGNAL_TERMINATION_TIMEOUT_SECONDS=60

verify_proof_bundle() { # bundle-directory
  python3 - "$1" \
    "$PROOF_MAX_SUMMARY_BYTES" "$PROOF_MAX_VERDICTS_BYTES" \
    "$PROOF_MAX_RECORD_BYTES" "$PROOF_MAX_RECORDS" \
    "$PROOF_MAX_LANE_LOG_BYTES" "$PROOF_MAX_TOTAL_LOG_BYTES" \
    "$PROOF_MAX_SNAPSHOT_BYTES" "$PROOF_MAX_BUNDLE_ENTRIES" \
    "$PROOF_MAX_TIMEOUT_SECONDS" "$PROOF_LOG_METADATA_RESERVE_BYTES" <<'PY'
import hashlib
import json
import os
import pathlib
import re
import stat
import sys
import time

(
    root_text,
    max_summary_text,
    max_verdicts_text,
    max_record_text,
    max_records_text,
    max_log_text,
    max_total_log_text,
    max_snapshot_text,
    max_entries_text,
    max_timeout_text,
    supervisor_metadata_reserve_text,
) = sys.argv[1:]
MAX_SUMMARY_BYTES = int(max_summary_text)
MAX_VERDICTS_BYTES = int(max_verdicts_text)
MAX_RECORD_BYTES = int(max_record_text)
MAX_RECORDS = int(max_records_text)
MAX_LOG_BYTES = int(max_log_text)
MAX_TOTAL_LOG_BYTES = int(max_total_log_text)
MAX_SNAPSHOT_BYTES = int(max_snapshot_text)
MAX_BUNDLE_ENTRIES = int(max_entries_text)
MAX_TIMEOUT_SECONDS = int(max_timeout_text)
SUPERVISOR_METADATA_RESERVE_BYTES = int(supervisor_metadata_reserve_text)
EXPECTED_NO_CLAIM = (
    "Software/protocol structure only; no physical validation or emergent "
    "Euler-disc prediction."
)
EXPECTED_CONCURRENCY_NO_CLAIM = (
    "Clean HEAD is observed at bookends; transient concurrent edits between "
    "observations are not excluded."
)
EXPECTED_PROOF_BOUNDARY_DETAIL = (
    "software/protocol structure only; no physical validation or emergent prediction"
)
EXPECTED_PROOF_BOUNDARY_LOG = (
    "software/protocol structure only\n"
    "no physical validation\n"
    "no mechanism attribution\n"
    "no emergent Euler-disc prediction\n"
    "retained_log_checker_smoke_command=cargo test --locked -p fs-euler-disc-e2e "
    "--test scientific_contract -- "
    "g0_check_receipts_logs_and_domain_separation_are_exact_and_bounded --exact "
    "--test-threads=1\n"
    "retained command is not packet/case replay and resolves no artifacts\n"
).encode("ascii")
EXPECTED_FOCUSED_SOURCE_MANIFEST_DETAIL = (
    "focused profile does not claim untracked/dirty Euler-disc paths are covered"
)
EXPECTED_FOCUSED_SOURCE_MANIFEST_LOG = (
    b"NO_DATA: focused profile does not cover untracked or dirty candidate paths\n"
)
EXPECTED_CLOSURE_SOURCE_MANIFEST_DETAIL = (
    "docs/source-manifest gates skipped because closure preconditions failed"
)
EXPECTED_CLOSURE_SOURCE_MANIFEST_LOG = (
    b"NO_DATA: closure refused before source-manifest evaluation\n"
)
HEX_256 = re.compile(r"[0-9a-f]{64}\Z")

root_argument = pathlib.Path(root_text)

def require(condition, detail):
    if not condition:
        raise SystemExit(f"proof-bundle verification failed: {detail}")

directory_flags = (
    os.O_RDONLY
    | getattr(os, "O_CLOEXEC", 0)
    | getattr(os, "O_DIRECTORY", 0)
    | getattr(os, "O_NOFOLLOW", 0)
)
file_flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
try:
    root_lstat = root_argument.lstat()
except OSError as error:
    raise SystemExit(
        f"proof-bundle verification failed: cannot inspect bundle root: {error}"
    ) from error
require(stat.S_ISDIR(root_lstat.st_mode), "bundle root is not a real directory")
try:
    root_descriptor = os.open(root_argument, directory_flags)
except OSError as error:
    raise SystemExit(
        f"proof-bundle verification failed: cannot open bundle root without following links: {error}"
    ) from error
root_opened = os.fstat(root_descriptor)
require(
    (root_opened.st_dev, root_opened.st_ino) == (root_lstat.st_dev, root_lstat.st_ino),
    "bundle root changed identity while opening",
)
root = root_argument.absolute()
verified_file_digests = {}

def metadata_tuple(value):
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )

def open_relative_regular(relative_name, context):
    relative = pathlib.PurePosixPath(relative_name)
    require(
        relative.parts
        and not relative.is_absolute()
        and ".." not in relative.parts
        and "\\" not in relative_name,
        f"unsafe path for {context}",
    )
    parent = os.dup(root_descriptor)
    try:
        for part in relative.parts[:-1]:
            child = os.open(part, directory_flags, dir_fd=parent)
            os.close(parent)
            parent = child
        descriptor = os.open(relative.parts[-1], file_flags, dir_fd=parent)
    except OSError as error:
        raise SystemExit(
            f"proof-bundle verification failed: cannot open {context} without following links: {error}"
        ) from error
    finally:
        os.close(parent)
    opened = os.fstat(descriptor)
    require(stat.S_ISREG(opened.st_mode), f"{context} is not a regular file")
    require(opened.st_nlink == 1, f"{context} is multiply linked")
    return descriptor, opened

def bounded_regular_bytes(relative_name, limit, context, *, remember=True):
    descriptor, opened = open_relative_regular(relative_name, context)
    try:
        require(opened.st_size <= limit, f"{context} exceeds the {limit}-byte bound")
        chunks = []
        remaining = limit + 1
        while remaining > 0:
            chunk = os.read(descriptor, min(65536, remaining))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        data = b"".join(chunks)
        require(len(data) <= limit, f"{context} exceeds the {limit}-byte bound")
        after = os.fstat(descriptor)
        require(metadata_tuple(after) == metadata_tuple(opened), f"{context} changed while reading")
        require(len(data) == opened.st_size, f"{context} was not read completely")
        if remember:
            verified_file_digests[relative_name] = hashlib.sha256(data).digest()
        return data
    finally:
        os.close(descriptor)

def inventory_bundle():
    observed = {}
    pending = [(os.dup(root_descriptor), "")]
    observed_entries = 0
    while pending:
        directory, prefix = pending.pop()
        try:
            entry_names = []
            with os.scandir(directory) as iterator:
                for entry in iterator:
                    observed_entries += 1
                    require(
                        observed_entries <= MAX_BUNDLE_ENTRIES,
                        f"bundle exceeds the {MAX_BUNDLE_ENTRIES}-entry inventory bound",
                    )
                    entry_names.append(entry.name)
            entry_names.sort()
            for entry_name in entry_names:
                relative_name = f"{prefix}/{entry_name}" if prefix else entry_name
                try:
                    entry_stat = os.stat(
                        entry_name,
                        dir_fd=directory,
                        follow_symlinks=False,
                    )
                except OSError as error:
                    raise SystemExit(
                        "proof-bundle verification failed: cannot inspect bundle entry "
                        f"{relative_name}: {error}"
                    ) from error
                require(
                    not stat.S_ISLNK(entry_stat.st_mode),
                    f"bundle inventory contains symlink {relative_name}",
                )
                if stat.S_ISDIR(entry_stat.st_mode):
                    try:
                        child = os.open(entry_name, directory_flags, dir_fd=directory)
                    except OSError as error:
                        raise SystemExit(
                            "proof-bundle verification failed: cannot open bundle directory "
                            f"{relative_name}: {error}"
                        ) from error
                    child_stat = os.fstat(child)
                    require(
                        metadata_tuple(child_stat) == metadata_tuple(entry_stat),
                        f"bundle directory {relative_name} changed while opening",
                    )
                    observed[relative_name] = ("directory", metadata_tuple(entry_stat))
                    pending.append((child, relative_name))
                elif stat.S_ISREG(entry_stat.st_mode):
                    require(
                        entry_stat.st_nlink == 1,
                        f"bundle inventory contains multiply linked file {relative_name}",
                    )
                    observed[relative_name] = ("file", metadata_tuple(entry_stat))
                else:
                    require(False, f"bundle inventory contains non-regular entry {relative_name}")
        finally:
            os.close(directory)
    return observed

initial_root_metadata = metadata_tuple(root_opened)
initial_inventory = inventory_bundle()

def reject_duplicate_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result

def parse_bounded_integer(text):
    digits = text[1:] if text.startswith("-") else text
    if len(digits) > 20:
        raise ValueError("JSON integer exceeds the 20-digit bound")
    return int(text)

def reject_float(text):
    raise ValueError(f"JSON floating constant {text!r} is forbidden")

def reject_non_finite(text):
    raise ValueError(f"non-finite JSON constant {text!r} is forbidden")

def strict_json_object(raw, context):
    require(len(raw) <= MAX_RECORD_BYTES, f"{context} exceeds the record bound")
    require(raw.endswith(b"\n"), f"{context} is not newline terminated")
    try:
        text = raw.decode("utf-8", "strict")
        value = json.loads(
            text,
            object_pairs_hook=reject_duplicate_keys,
            parse_int=parse_bounded_integer,
            parse_float=reject_float,
            parse_constant=reject_non_finite,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise SystemExit(
            f"proof-bundle verification failed: {context} is not strict JSON: {error}"
        ) from error
    require(type(value) is dict, f"{context} is not a JSON object")
    encoded = (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    require(encoded == raw, f"{context} is not canonical JSON")
    return value

def exact_keys(value, expected, context):
    require(set(value) == expected, f"{context} does not have the exact field set")

def exact_string(value, field, context, *, empty=False, maximum=4096):
    require(type(value) is str, f"{context} field {field} is not a string")
    require(empty or bool(value), f"{context} field {field} is empty")
    require(len(value.encode("utf-8")) <= maximum, f"{context} field {field} is oversized")

def exact_integer(value, field, context, *, minimum=0, maximum=2**31 - 1):
    require(type(value) is int, f"{context} field {field} is not an integer")
    require(minimum <= value <= maximum, f"{context} field {field} is out of range")

def exact_boolean(value, field, context):
    require(type(value) is bool, f"{context} field {field} is not a boolean")

summary_bytes = bounded_regular_bytes("summary.json", MAX_SUMMARY_BYTES, "summary.json")
verdict_bytes = bounded_regular_bytes("verdicts.jsonl", MAX_VERDICTS_BYTES, "verdicts.jsonl")
prefix_bytes = bounded_regular_bytes(
    "verdicts-prefix.jsonl", MAX_VERDICTS_BYTES, "verdicts-prefix.jsonl"
)
lines = verdict_bytes.splitlines(keepends=True)
require(lines, "verdicts.jsonl is empty")
require(all(line.endswith(b"\n") for line in lines), "every verdict must be newline terminated")
require(len(lines) <= MAX_RECORDS, f"bundle exceeds the {MAX_RECORDS}-record bound")
seal_line = lines[-1]
prefix = b"".join(lines[:-1])
require(prefix_bytes == prefix, "verdicts-prefix.jsonl is not the exact verdict prefix")
all_records = [
    strict_json_object(line, f"verdict record {index + 1}")
    for index, line in enumerate(lines)
]
seal = all_records[-1]
summary = strict_json_object(summary_bytes, "summary.json")
lane_keys = {
    "schema", "lane", "status", "authority", "detail", "log", "log_bytes",
    "log_sha256", "head", "host_isa", "profile", "executor_declaration",
    "executor_attestation", "provenance_state", "terminal",
}
seal_keys = {
    "schema", "record_type", "status", "checks", "failures", "head",
    "host_isa", "profile", "executor_declaration", "executor_attestation",
    "source_manifest_coverage", "candidate_ready_for_dsr", "dsr_status",
    "proof_scope", "provenance_state", "snapshot_before_sha256",
    "snapshot_after_sha256", "verdicts_prefix_sha256", "verdicts",
    "proof_seal_locator", "terminal_exit_code", "terminal", "no_claim",
    "concurrency_no_claim",
}
exact_keys(seal, seal_keys, "proof seal")
exact_keys(summary, seal_keys, "summary.json")
require(sum(record.get("record_type") == "proof-seal" for record in all_records) == 1,
        "bundle must contain exactly one proof seal")
require(seal.get("record_type") == "proof-seal" and seal.get("terminal") is True,
        "last verdict is not a terminal proof seal")
require(seal.get("schema") == "frankensim.euler-disc-contract-e2e.proof-seal.v2",
        "unknown proof-seal schema")
require(summary.get("record_type") == "proof-seal" and summary.get("terminal") is True,
        "summary.json is not the terminal proof seal")
require(summary_bytes == seal_line, "summary.json is not byte-identical to the final proof seal")
require(summary.get("proof_seal_locator") == "verdicts.jsonl#last-line",
        "summary proof-seal locator is not canonical")
require(summary.get("verdicts") == "verdicts.jsonl",
        "summary verdict locator is not canonical")
require(summary.get("executor_attestation") == "caller-declared-unverified",
        "summary executor attestation is not the exact v2 declaration")
require(summary.get("dsr_status") == "not-run-by-this-harness",
        "summary DSR status is not the exact v2 no-claim")
require(summary.get("no_claim") == EXPECTED_NO_CLAIM,
        "summary software/physical no-claim changed")
require(summary.get("concurrency_no_claim") == EXPECTED_CONCURRENCY_NO_CLAIM,
        "summary concurrency no-claim changed")
require(hashlib.sha256(prefix).hexdigest() == seal.get("verdicts_prefix_sha256"),
        "verdict prefix hash mismatch")
require(summary.get("verdicts_prefix_sha256") == seal.get("verdicts_prefix_sha256"),
        "summary and seal disagree about verdict-prefix identity")
for field in (
    "status", "head", "profile", "executor_declaration", "proof_scope",
    "provenance_state", "snapshot_before_sha256", "snapshot_after_sha256",
        "terminal_exit_code",
):
    require(summary.get(field) == seal.get(field), f"summary/seal mismatch for {field}")

for field in (
    "schema", "record_type", "status", "head", "host_isa", "profile",
    "executor_declaration", "executor_attestation", "source_manifest_coverage",
    "dsr_status", "proof_scope", "provenance_state", "snapshot_before_sha256",
    "snapshot_after_sha256", "verdicts_prefix_sha256", "verdicts",
    "proof_seal_locator", "no_claim", "concurrency_no_claim",
):
    exact_string(summary[field], field, "summary.json", empty=field.startswith("snapshot_"))
for field in ("checks", "failures"):
    exact_integer(summary[field], field, "summary.json", maximum=MAX_RECORDS)
exact_integer(summary["terminal_exit_code"], "terminal_exit_code", "summary.json", maximum=255)
exact_boolean(summary["candidate_ready_for_dsr"], "candidate_ready_for_dsr", "summary.json")
exact_boolean(summary["terminal"], "terminal", "summary.json")
require(HEX_256.fullmatch(summary["verdicts_prefix_sha256"]) is not None,
        "verdict prefix hash is not lowercase SHA-256 hex")
for field in ("snapshot_before_sha256", "snapshot_after_sha256"):
    require(not summary[field] or HEX_256.fullmatch(summary[field]) is not None,
            f"{field} is not empty or lowercase SHA-256 hex")

lane_names = set()
lane_order = []
lane_statuses = {}
lane_records = {}
lane_log_data = {}
lane_logs = set()
failed_lanes = 0
total_log_bytes = 0
for record in all_records[:-1]:
    exact_keys(record, lane_keys, f"lane record {record.get('lane', '<missing>')}")
    require(record.get("schema") == "frankensim.euler-disc-contract-e2e.verdict.v1",
            "unknown nonterminal verdict schema")
    if "lane" not in record:
        raise SystemExit("proof-bundle verification failed: non-seal record lacks lane")
    require(record["lane"] not in lane_names, f"duplicate lane {record['lane']}")
    lane_names.add(record["lane"])
    lane_order.append(record["lane"])
    require(record.get("status") in {"PASS", "FAIL", "NO_DATA"},
            f"unknown lane status for {record['lane']}")
    lane_statuses[record["lane"]] = record["status"]
    lane_records[record["lane"]] = record
    for field in (
        "schema", "lane", "status", "authority", "detail", "log", "log_sha256",
        "head", "host_isa", "profile", "executor_declaration",
        "executor_attestation", "provenance_state",
    ):
        exact_string(record[field], field, f"lane {record['lane']}")
    exact_integer(
        record["log_bytes"], "log_bytes", f"lane {record['lane']}", maximum=MAX_LOG_BYTES
    )
    exact_boolean(record["terminal"], "terminal", f"lane {record['lane']}")
    require(HEX_256.fullmatch(record["log_sha256"]) is not None,
            f"lane {record['lane']} has a malformed log hash")
    for field in ("head", "host_isa", "profile", "executor_declaration"):
        require(record.get(field) == summary.get(field),
                f"lane {record['lane']} is not bound to summary field {field}")
    require(record.get("executor_attestation") == "caller-declared-unverified",
            f"lane {record['lane']} has an unknown executor attestation")
    require(record.get("provenance_state") == "provisional" and
            record.get("terminal") is False,
            f"lane {record['lane']} has invalid nonterminal provenance state")
    failed_lanes += record.get("status") == "FAIL"
    raw_relative = record["log"]
    relative = pathlib.PurePosixPath(raw_relative)
    require(not relative.is_absolute() and ".." not in relative.parts and "\\" not in raw_relative,
            f"unsafe lane log path for {record['lane']}")
    require(relative.parts and relative.parts[0] == "logs" and relative.as_posix() == raw_relative,
            f"noncanonical lane log path for {record['lane']}")
    require(raw_relative not in lane_logs, f"duplicate lane log path {raw_relative}")
    lane_logs.add(raw_relative)
    total_log_bytes += record["log_bytes"]
    require(total_log_bytes <= MAX_TOTAL_LOG_BYTES, "aggregate retained logs exceed the bound")
    data = bounded_regular_bytes(raw_relative, MAX_LOG_BYTES, f"lane log {raw_relative}")
    lane_log_data[record["lane"]] = data
    require(len(data) == record["log_bytes"], f"lane log size mismatch for {record['lane']}")
    require(hashlib.sha256(data).hexdigest() == record["log_sha256"],
            f"lane log hash mismatch for {record['lane']}")

require(summary.get("checks") == len(lane_names), "summary check count is not derived from lanes")
require(summary.get("failures") == failed_lanes, "summary failure count is not derived from lanes")
for field, filename in (
    ("snapshot_before_sha256", "snapshot-before.txt"),
    ("snapshot_after_sha256", "snapshot-after.txt"),
):
    digest = summary.get(field)
    if digest:
        data = bounded_regular_bytes(filename, MAX_SNAPSHOT_BYTES, filename)
        require(hashlib.sha256(data).hexdigest() == digest, f"{filename} hash mismatch")

expected_files = {
    "summary.json",
    "verdicts.jsonl",
    "verdicts-prefix.jsonl",
    *lane_logs,
}
if summary.get("snapshot_before_sha256"):
    expected_files.add("snapshot-before.txt")
if summary.get("snapshot_after_sha256"):
    expected_files.add("snapshot-after.txt")
expected_directories = set()
for relative_name in expected_files:
    parent = pathlib.PurePosixPath(relative_name).parent
    while parent != pathlib.PurePosixPath("."):
        expected_directories.add(parent.as_posix())
        parent = parent.parent

observed_files = {
    relative_name
    for relative_name, (kind, _) in initial_inventory.items()
    if kind == "file"
}
observed_directories = {
    relative_name
    for relative_name, (kind, _) in initial_inventory.items()
    if kind == "directory"
}

require(
    observed_files == expected_files,
    "bundle file inventory differs from the seal-derived exact inventory: "
    f"missing={sorted(expected_files - observed_files)!r} "
    f"unexpected={sorted(observed_files - expected_files)!r}",
)
require(
    observed_directories == expected_directories,
    "bundle directory inventory differs from the seal-derived exact inventory: "
    f"missing={sorted(expected_directories - observed_directories)!r} "
    f"unexpected={sorted(observed_directories - expected_directories)!r}",
)
status = summary.get("status")
exit_code = summary.get("terminal_exit_code")
ready = summary.get("candidate_ready_for_dsr")
allowed_terminal_statuses = {
    "READY_FOR_DSR", "FOCUSED_PASS", "SELF_TEST_PASS", "NO_DATA", "FAIL",
    "SELF_TEST_FAIL", "INTERRUPTED", "INCOMPLETE",
}
require(status in allowed_terminal_statuses, f"unknown terminal status {status!r}")
require(ready == (status == "READY_FOR_DSR"),
        "candidate_ready_for_dsr does not exactly match READY_FOR_DSR status")
core_sequence = (
    "crate-fmt",
    "crate-check",
    "retained-log-checker-smoke",
    "retained-log-checker-smoke-sentinel",
    "crate-unit-integration",
    "crate-doctest-hostile-boundary",
    "crate-clippy",
    "xtask-check-layers",
    "xtask-check-deps",
    "xtask-check-contracts",
    "xtask-check-schemas",
    "xtask-check-consolidation",
    "xtask-check-identities",
    "xtask-check-goldens",
)
focused_sequence = (
    "proof-boundary",
    "source-manifest",
    "constellation-verify",
    "constellation-snapshot-before",
    *core_sequence,
    "constellation-snapshot-after",
    "source-stability",
)
closure_sequence = (
    "proof-boundary",
    "closure-root-preflight",
    "constellation-verify",
    "constellation-snapshot-before",
    *core_sequence,
    "xtask-check-docs",
    "xtask-check-source-manifest",
    "source-manifest-membership",
    "constellation-snapshot-after",
    "closure-root-bookend",
    "source-stability",
)
closure_refusal_sequence = (
    "proof-boundary",
    "closure-root-preflight",
    "source-manifest",
)
self_test_sequence = (
    "self-test-ordinary-pass",
    "self-test-exact-nonzero",
    "self-test-reserved-child-exit",
    "self-test-child-signal-exit",
    "self-test-spawn-failure",
    "self-test-timeout-stubborn-group",
    "self-test-leader-exit-live-group",
    "self-test-setpgid-escape-drained",
    "self-test-output-truncation-drained",
    "self-test-cap-plus-timeout",
    "self-test-lane-log-cap-boundary",
    "self-test-numeric-config-boundaries",
    "self-test-completion-classification-signal",
    "self-test-snapshot-byte-boundaries",
    "self-test-snapshot-hash-failure",
    "self-test-timeout-hanging-helper",
    "self-test-supervisor-exception",
    "self-test-special-index-flag-refusal",
    "self-test-skip-worktree-flag-refusal",
    "self-test-fsmonitor-valid-flag-refusal",
    "self-test-real-repository-read-failure",
    "self-test-zero-smoke-refusal",
    "self-test-invalid-consolidation-disposition-refusal",
    "self-test-consolidation-scope-mutation-refusal",
    "self-test-success-deadline-publication-refusal",
    "self-test-postpublication-deadline-retraction",
    "self-test-valid-bundle",
    "self-test-valid-success-bundle",
    "self-test-valid-closure-ready-bundle",
    "self-test-closure-ready-missing-lane-refusal",
    "self-test-closure-ready-snapshot-mismatch-refusal",
    "self-test-valid-self-test-pass-bundle",
    "self-test-self-test-pass-authority-refusal",
    "self-test-self-test-pass-detail-refusal",
    "self-test-self-test-pass-log-refusal",
    "self-test-wrong-authority-refusal",
    "self-test-wrong-command-refusal",
    "self-test-supervisor-state-contradiction-refusal",
    "self-test-supervisor-containment-contradiction-refusal",
    "self-test-supervisor-flag-contradiction-refusal",
    "self-test-supervisor-deadline-contradiction-refusal",
    "self-test-proof-boundary-body-refusal",
    "self-test-source-manifest-body-refusal",
    "self-test-failed-control-continuation-refusal",
    "self-test-premature-snapshot-refusal",
    "self-test-toctou-mutation-refusal",
    "self-test-publication-gap-mutation-refusal",
    "self-test-publication-destination-race-refusal",
    "self-test-mutated-bundle-refusal",
    "self-test-duplicate-seal-refusal",
    "self-test-truncated-seal-refusal",
    "self-test-summary-mismatch-refusal",
    "self-test-unsafe-path-refusal",
    "self-test-duplicate-json-key-refusal",
    "self-test-nonfinite-json-refusal",
    "self-test-oversized-json-refusal",
    "self-test-record-count-refusal",
    "self-test-readiness-mismatch-refusal",
    "self-test-no-claim-mutation-refusal",
    "self-test-unknown-terminal-refusal",
    "self-test-inventory-entry-cap-refusal",
    "self-test-unreferenced-file-refusal",
    "self-test-prefix-file-mismatch-refusal",
    "self-test-interrupted-status-matrix-refusal",
    "self-test-incomplete-status-matrix-refusal",
    "self-test-closure-refusal-writer",
    "self-test-closure-infrastructure-failure",
    "self-test-aggregate-budget-exhaustion",
    "self-test-invalid-candidate-not-published",
    "self-test-incomplete-seal",
    "self-test-publication-fsync-failure-retains-bundle",
    "self-test-prepublication-signal",
    "self-test-publication-window-hup",
    "self-test-publication-window-int",
    "self-test-publication-window-term",
    "self-test-wrapper-arbitrary-term",
    "self-test-wrapper-hup",
    "self-test-wrapper-int",
    "self-test-wrapper-term",
)
focused_required = set(focused_sequence)
closure_required = set(closure_sequence)
self_test_required = set(self_test_sequence)
focused_statuses = {lane: "PASS" for lane in focused_required}
focused_statuses["source-manifest"] = "NO_DATA"
closure_statuses = {lane: "PASS" for lane in closure_required}
self_test_statuses = {lane: "PASS" for lane in self_test_required}
normal_executors = {"local", "rch", "dsr"}
profile = summary.get("profile")
executor = summary.get("executor_declaration")
coverage = summary.get("source_manifest_coverage")
provenance = summary.get("provenance_state")
proof_scope = summary.get("proof_scope")
before_snapshot = summary.get("snapshot_before_sha256")
after_snapshot = summary.get("snapshot_after_sha256")

def require_named_lane_status(lane, expected_status):
    require(lane in lane_names, f"terminal state is missing required lane {lane}")
    require(
        lane_statuses[lane] == expected_status,
        f"terminal lane {lane} must be {expected_status}",
    )

def require_wrapper_signal_record():
    signal_number = exit_code - 128
    record = lane_records["wrapper-signal"]
    require(record["authority"] == "bounded-interrupt-cleanup",
            "wrapper-signal has the wrong authority")
    require(
        record["detail"] ==
        f"wrapper received signal {signal_number} and stopped scheduling",
        "wrapper-signal has the wrong detail",
    )
    require(record["log"] == "logs/wrapper-signal.log",
            "wrapper-signal has the wrong log path")
    expected = (
        f"wrapper-signal={signal_number}\n"
        "no-later-lanes-launched=true\n"
    ).encode()
    require(lane_log_data["wrapper-signal"] == expected,
            "wrapper-signal log does not match the terminal signal")

def require_internal_incomplete_record(expected_exit_code=None):
    record = lane_records["internal-incomplete"]
    require(record["authority"] == "harness-integrity",
            "internal-incomplete has the wrong authority")
    require(record["detail"] == "EXIT trap sealed an otherwise incomplete run",
            "internal-incomplete has the wrong detail")
    require(record["log"] == "logs/internal-incomplete.log",
            "internal-incomplete has the wrong log path")
    data = lane_log_data["internal-incomplete"]
    if expected_exit_code is None:
        match = re.fullmatch(
            rb"unexpected-exit-code=([1-9][0-9]{0,2})\n"
            rb"terminal-seal-was-missing=true\n",
            data,
        )
        require(match is not None and int(match.group(1)) <= 255,
                "internal-incomplete log has an invalid pre-transition exit")
    else:
        expected = (
            f"unexpected-exit-code={expected_exit_code}\n"
            "terminal-seal-was-missing=true\n"
        ).encode()
        require(data == expected,
                "internal-incomplete log does not match the terminal exit")

def require_profile_scope():
    if profile == "focused":
        require(proof_scope == "focused-software-only", "focused state has wrong scope")
    elif profile == "closure":
        require(
            proof_scope == "head-bookended-closure-candidate",
            "closure state has wrong scope",
        )
    else:
        require(False, "terminal state has an unknown profile")

scope_paths = (
    "Cargo.toml", "Cargo.lock", "README.md", "consolidation-review.json",
    "docs/CONSOLIDATION_REVIEW.md", "docs/CONVENTIONS.md", "docs/SCHEMA_POLICY.md",
    "doc-facts-inventory.json", "schema-policy.json", "identity-authorities.json",
    "identity-schemas.json", "golden-couplings.json", "xtask/src/identities.rs",
    "scripts/ci/euler_disc_contract_e2e.sh", "crates/fs-euler-disc-e2e/Cargo.toml",
    "crates/fs-euler-disc-e2e/CONTRACT.md", "crates/fs-euler-disc-e2e/src/lib.rs",
    "crates/fs-euler-disc-e2e/src/contract.rs", "crates/fs-euler-disc-e2e/src/protocol.rs",
    "crates/fs-euler-disc-e2e/tests/scientific_contract.rs",
)
normal_authorities = {
    "proof-boundary": "declaration-only",
    "closure-root-preflight": "full-root-clean-head-observation",
    "constellation-verify": "constellation-source-preflight",
    "constellation-snapshot-before": "source-provenance-snapshot",
    "crate-fmt": "static-hygiene-only",
    "crate-check": "focused-software-evidence",
    "retained-log-checker-smoke": "focused-software-evidence",
    "retained-log-checker-smoke-sentinel": "non-vacuity-evidence",
    "crate-unit-integration": "focused-software-evidence",
    "crate-doctest-hostile-boundary": "compile-fail-api-evidence",
    "crate-clippy": "static-hygiene-only",
    "xtask-check-layers": "workspace-structural-gate",
    "xtask-check-deps": "workspace-structural-gate",
    "xtask-check-contracts": "workspace-structural-gate",
    "xtask-check-schemas": "workspace-structural-gate",
    "xtask-check-consolidation": "workspace-structural-gate",
    "xtask-check-identities": "workspace-structural-gate",
    "xtask-check-goldens": "workspace-structural-gate",
    "xtask-check-docs": "workspace-documentation-gate",
    "xtask-check-source-manifest": "head-bound-source-inventory",
    "source-manifest-membership": "independent-path-membership",
    "constellation-snapshot-after": "source-provenance-snapshot",
    "closure-root-bookend": "full-root-clean-head-observation",
    "source-stability": "source-provenance-snapshot",
}
command_lanes = set(normal_authorities) - {"proof-boundary"}
command_failure_details = {
    "aggregate run timeout exhausted before command launch",
    "aggregate retained-log budget exhausted before command launch",
    "command output was truncated after its process group drained",
    "command leader exited before same-session descendants; the supervisor drained the group",
    "command exceeded its monotonic deadline",
    "process-group drain or supervisor metadata integrity could not be established",
    "command could not be launched",
    "command executable was not found",
    "wrapper signal interrupted the command after bounded process-group cleanup",
}
command_run_metadata = None
supervisor_result_keys = {
    "schema",
    "configured_lane_timeout_seconds",
    "run_deadline_monotonic_ns",
    "supervisor_started_monotonic_ns",
    "lane_deadline_monotonic_ns",
    "effective_deadline_monotonic_ns",
    "deadline_kind",
    "configured_output_cap_bytes",
    "shell_effective_output_cap_bytes",
    "retained_output_cap_bytes",
    "total_log_cap_bytes",
    "initial_log_bytes",
    "leader_exit_code",
    "wrapper_exit_code",
    "shutdown_reason",
    "interrupted_signal",
    "process_group_drained",
    "process_session_drained",
    "output_pipe_eof",
    "output_truncated",
    "shell_cap_reduced",
    "python_cap_reduced",
    "metadata_complete",
    "inspection_error_count",
    "term_sent",
    "kill_sent",
}

def expected_authority(lane):
    if lane == "source-manifest":
        return (
            "manifest-not-evaluated"
            if profile == "focused"
            else "closure-refused-before-manifest-evaluation"
        )
    return normal_authorities.get(lane)

def require_command_argv(lane, argv):
    cargo_tails = {
        "crate-fmt": ["fmt", "-p", "fs-euler-disc-e2e", "--check"],
        "crate-check": ["check", "--locked", "-p", "fs-euler-disc-e2e", "--all-targets"],
        "retained-log-checker-smoke": [
            "test", "--locked", "-p", "fs-euler-disc-e2e", "--test",
            "scientific_contract", "--",
            "g0_check_receipts_logs_and_domain_separation_are_exact_and_bounded",
            "--exact", "--test-threads=1",
        ],
        "crate-unit-integration": [
            "test", "--locked", "--no-fail-fast", "-p", "fs-euler-disc-e2e",
            "--lib", "--tests", "--", "--test-threads=1",
        ],
        "crate-doctest-hostile-boundary": [
            "test", "--locked", "-p", "fs-euler-disc-e2e", "--doc",
        ],
        "crate-clippy": [
            "clippy", "--locked", "-p", "fs-euler-disc-e2e", "--all-targets",
            "--no-deps", "--", "-D", "warnings",
        ],
    }
    for gate in (
        "check-layers", "check-deps", "check-contracts", "check-schemas",
        "check-consolidation", "check-identities", "check-goldens", "check-docs",
        "check-source-manifest",
    ):
        cargo_tails[f"xtask-{gate}"] = ["run", "--locked", "-q", "-p", "xtask", "--", gate]
    if lane in cargo_tails:
        require(len(argv) >= 2 and argv[0], f"{lane} command lacks an executable")
        require(argv[1:] == cargo_tails[lane], f"{lane} command arguments are not exact")
    elif lane == "constellation-verify":
        require(argv == ["scripts/ci/checkout_constellation.sh", "--verify-only"],
                "constellation-verify command arguments are not exact")
    elif lane in {"closure-root-preflight", "closure-root-bookend"}:
        require(
            argv == ["bash", "-c", 'closure_root_preflight_command "$@"', "_", summary["head"], *scope_paths],
            f"{lane} command arguments are not exact",
        )
    elif lane in {"constellation-snapshot-before", "constellation-snapshot-after"}:
        expected_name = "snapshot-before.txt" if lane.endswith("before") else "snapshot-after.txt"
        require(
            len(argv) == 7 + len(scope_paths)
            and argv[:4] == ["bash", "-c", 'capture_snapshot_command "$@"', "_"]
            and pathlib.PurePosixPath(argv[4]).name == expected_name
            and pathlib.PurePosixPath(argv[5]).name.startswith(f".{lane}-candidate-")
            and argv[6] == summary["head"]
            and tuple(argv[7:]) == scope_paths,
            f"{lane} command arguments are not exact",
        )
    elif lane == "retained-log-checker-smoke-sentinel":
        require(
            len(argv) == 6
            and argv[:4] == ["bash", "-c", 'checker_smoke_sentinel_command "$@"', "_"]
            and pathlib.PurePosixPath(argv[4]).name == "retained-log-checker-smoke.log"
            and argv[5] == "g0_check_receipts_logs_and_domain_separation_are_exact_and_bounded",
            "retained-log checker sentinel command arguments are not exact",
        )
    elif lane == "source-manifest-membership":
        require(
            len(argv) == 5 + len(scope_paths)
            and argv[:4] == ["bash", "-c", 'source_manifest_membership_command "$@"', "_"]
            and pathlib.PurePosixPath(argv[4]).name == "frankensim-source-manifest.json"
            and tuple(argv[5:]) == scope_paths,
            "source-manifest-membership command arguments are not exact",
        )
    elif lane == "source-stability":
        require(
            len(argv) == 4 and argv[:2] == ["cmp", "-s"]
            and pathlib.PurePosixPath(argv[2]).name == "snapshot-before.txt"
            and pathlib.PurePosixPath(argv[3]).name == "snapshot-after.txt",
            "source-stability command arguments are not exact",
        )
    else:
        require(False, f"no command specification exists for lane {lane}")

def require_command_metadata(lane):
    global command_run_metadata
    data = lane_log_data[lane]
    parts = data.split(b"\n", 14)
    require(len(parts) == 15 and parts[14].startswith(b"--- command output ---\n"),
            f"{lane} log lacks the exact command.v1 header")
    decoded = []
    for index, raw in enumerate(parts[:14]):
        try:
            decoded.append(raw.decode("ascii"))
        except UnicodeDecodeError as error:
            require(False, f"{lane} command header field {index} is not ASCII: {error}")
    expected_pairs = (
        ("schema", "frankensim.euler-disc-contract-e2e.command.v1"),
        ("lane", lane),
        ("head", summary["head"]),
        ("host_isa", summary["host_isa"]),
        ("profile", profile),
        ("executor_declaration", executor),
        ("executor_attestation", "caller-declared-unverified"),
    )
    for index, (key, value) in enumerate(expected_pairs):
        require(decoded[index] == f"{key}={value}", f"{lane} command header has wrong {key}")
    numeric = {}
    for index, key in enumerate((
        "lane_timeout_seconds", "run_timeout_seconds", "run_started_monotonic_ns",
        "run_deadline_monotonic_ns", "lane_log_max_bytes",
    ), start=7):
        prefix = f"{key}="
        require(decoded[index].startswith(prefix), f"{lane} command header lacks {key}")
        value = decoded[index][len(prefix):]
        require(re.fullmatch(r"[1-9][0-9]*", value) is not None,
                f"{lane} command header {key} is not canonical")
        numeric[key] = int(value)
    require(numeric["lane_timeout_seconds"] <= MAX_TIMEOUT_SECONDS,
            f"{lane} command lane timeout exceeds the verifier bound")
    require(numeric["run_timeout_seconds"] <= MAX_TIMEOUT_SECONDS,
            f"{lane} command run timeout exceeds the verifier bound")
    require(numeric["lane_log_max_bytes"] <= MAX_LOG_BYTES,
            f"{lane} command log cap exceeds the verifier bound")
    require(
        numeric["run_deadline_monotonic_ns"] - numeric["run_started_monotonic_ns"]
        == numeric["run_timeout_seconds"] * 1_000_000_000,
        f"{lane} command monotonic deadline is inconsistent",
    )
    shared = tuple(numeric.values())
    if command_run_metadata is None:
        command_run_metadata = shared
    else:
        require(command_run_metadata == shared, f"{lane} command run metadata drifted")
    require(decoded[12] == f"authority={lane_records[lane]['authority']}",
            f"{lane} command authority header disagrees with its verdict")
    require(decoded[13].startswith("argv_json="), f"{lane} command header lacks argv_json")
    raw_argv = decoded[13][len("argv_json="):]
    try:
        argv = json.loads(raw_argv)
    except (json.JSONDecodeError, ValueError) as error:
        require(False, f"{lane} argv_json is invalid: {error}")
    require(type(argv) is list and argv and all(type(item) is str for item in argv),
            f"{lane} argv_json is not a nonempty string array")
    require(json.dumps(argv, ensure_ascii=True, separators=(",", ":")) == raw_argv,
            f"{lane} argv_json is not canonical")
    require_command_argv(lane, argv)
    result_prefix = b"supervisor_result_json="
    result_lines = [
        line for line in data.splitlines(keepends=True)
        if line.startswith(result_prefix)
    ]
    require(len(result_lines) == 1, f"{lane} log must contain exactly one supervisor result")
    require(data.endswith(result_lines[0]), f"{lane} supervisor result is not the final log line")
    result = strict_json_object(
        result_lines[0][len(result_prefix):],
        f"{lane} supervisor result",
    )
    exact_keys(result, supervisor_result_keys, f"{lane} supervisor result")
    require(
        result["schema"] == "frankensim.euler-disc-contract-e2e.supervisor-result.v1",
        f"{lane} supervisor result has an unknown schema",
    )
    integer_fields = (
        "configured_lane_timeout_seconds",
        "run_deadline_monotonic_ns",
        "supervisor_started_monotonic_ns",
        "lane_deadline_monotonic_ns",
        "effective_deadline_monotonic_ns",
        "configured_output_cap_bytes",
        "shell_effective_output_cap_bytes",
        "retained_output_cap_bytes",
        "total_log_cap_bytes",
        "initial_log_bytes",
        "wrapper_exit_code",
        "inspection_error_count",
    )
    for field in integer_fields:
        exact_integer(
            result[field], field, f"{lane} supervisor result",
            maximum=2**63 - 1,
        )
    for field in (
        "process_group_drained", "process_session_drained", "output_pipe_eof",
        "output_truncated", "shell_cap_reduced", "python_cap_reduced",
        "metadata_complete", "term_sent", "kill_sent",
    ):
        exact_boolean(result[field], field, f"{lane} supervisor result")
    for field in ("deadline_kind", "shutdown_reason"):
        exact_string(result[field], field, f"{lane} supervisor result")
    for field in ("leader_exit_code", "interrupted_signal"):
        require(
            result[field] is None or type(result[field]) is int,
            f"{lane} supervisor result field {field} is not null or an integer",
        )
    if result["leader_exit_code"] is not None:
        require(
            -255 <= result["leader_exit_code"] <= 255,
            f"{lane} supervisor leader exit code is out of range",
        )
    if result["interrupted_signal"] is not None:
        require(
            1 <= result["interrupted_signal"] <= 127,
            f"{lane} supervisor interrupt signal is out of range",
        )
    header_prefix = b"\n".join(parts[:14]) + b"\n--- command output ---\n"
    require(
        result["initial_log_bytes"] == len(header_prefix),
        f"{lane} supervisor initial log size disagrees with the exact command header",
    )
    require(
        result["configured_lane_timeout_seconds"] == numeric["lane_timeout_seconds"],
        f"{lane} supervisor configured timeout disagrees with its command header",
    )
    require(
        result["configured_output_cap_bytes"] == numeric["lane_log_max_bytes"],
        f"{lane} supervisor configured output cap disagrees with its command header",
    )
    require(
        result["run_deadline_monotonic_ns"] == numeric["run_deadline_monotonic_ns"],
        f"{lane} supervisor aggregate deadline disagrees with its command header",
    )
    require(
        result["supervisor_started_monotonic_ns"] >= numeric["run_started_monotonic_ns"],
        f"{lane} supervisor predates the declared run start",
    )
    require(
        result["lane_deadline_monotonic_ns"]
        == result["supervisor_started_monotonic_ns"]
        + result["configured_lane_timeout_seconds"] * 1_000_000_000,
        f"{lane} supervisor lane deadline is inconsistent",
    )
    require(
        result["effective_deadline_monotonic_ns"]
        == min(result["lane_deadline_monotonic_ns"], result["run_deadline_monotonic_ns"]),
        f"{lane} supervisor effective deadline is not the exact minimum",
    )
    expected_deadline_kind = (
        "aggregate"
        if result["run_deadline_monotonic_ns"] <= result["lane_deadline_monotonic_ns"]
        else "lane"
    )
    require(
        result["deadline_kind"] == expected_deadline_kind,
        f"{lane} supervisor deadline kind is inconsistent",
    )
    require(
        0 <= result["retained_output_cap_bytes"]
        <= result["shell_effective_output_cap_bytes"]
        <= result["configured_output_cap_bytes"]
        <= MAX_LOG_BYTES,
        f"{lane} supervisor output caps are not monotonically bounded",
    )
    require(
        result["initial_log_bytes"] <= result["total_log_cap_bytes"] <= MAX_LOG_BYTES,
        f"{lane} supervisor total log cap is invalid",
    )
    require(
        len(data) <= result["total_log_cap_bytes"],
        f"{lane} retained log exceeds its supervisor total cap",
    )
    expected_retained_output_cap = max(
        0,
        min(
            result["shell_effective_output_cap_bytes"],
            result["total_log_cap_bytes"]
            - result["initial_log_bytes"]
            - SUPERVISOR_METADATA_RESERVE_BYTES,
        ),
    )
    require(
        result["retained_output_cap_bytes"] == expected_retained_output_cap,
        f"{lane} supervisor retained-output cap is not the exact bounded remainder",
    )
    require(
        result["shell_cap_reduced"]
        == (result["shell_effective_output_cap_bytes"] < result["configured_output_cap_bytes"]),
        f"{lane} supervisor shell-cap reduction flag is inconsistent",
    )
    require(
        result["python_cap_reduced"]
        == (result["retained_output_cap_bytes"] < result["shell_effective_output_cap_bytes"]),
        f"{lane} supervisor Python-cap reduction flag is inconsistent",
    )
    rc = result["wrapper_exit_code"]
    reason = result["shutdown_reason"]
    if reason == "none":
        expected_detail = "command completed"
        require(rc == 0 and result["leader_exit_code"] == 0,
                f"{lane} successful supervisor result is contradictory")
        require(
            result["process_group_drained"]
            and result["process_session_drained"]
            and result["output_pipe_eof"]
            and result["metadata_complete"]
            and not result["output_truncated"]
            and not result["python_cap_reduced"]
            and result["inspection_error_count"] == 0,
            f"{lane} successful supervisor result lacks complete containment metadata",
        )
    elif reason == "retained-log-budget-before-launch":
        expected_detail = "aggregate retained-log budget exhausted before command launch"
        require(rc == 120,
                f"{lane} retained-log refusal has the wrong supervisor reason")
    elif reason == "aggregate-deadline-before-launch":
        expected_detail = "aggregate run timeout exhausted before command launch"
        require(rc == 121,
                f"{lane} aggregate deadline refusal has the wrong supervisor reason")
    elif reason == "output-truncated":
        expected_detail = "command output was truncated after its process group drained"
        require(rc == 122 and result["output_truncated"],
                f"{lane} output-truncation result is contradictory")
    elif reason == "leader-exit-with-live-group":
        expected_detail = "command leader exited before same-session descendants; the supervisor drained the group"
        require(rc == 123,
                f"{lane} live-descendant result has the wrong reason")
    elif reason == "timeout":
        expected_detail = "command exceeded its monotonic deadline"
        require(rc == 124, f"{lane} timeout result has the wrong reason")
    elif reason in {"supervisor-integrity-failure", "inspection-failure"}:
        expected_detail = "process-group drain or supervisor metadata integrity could not be established"
        require(rc == 125,
                f"{lane} integrity failure has the wrong supervisor reason")
    elif reason == "launch-error":
        expected_detail = "command could not be launched"
        require(rc == 126, f"{lane} launch failure has the wrong reason")
    elif reason == "launch-not-found":
        expected_detail = "command executable was not found"
        require(rc == 127, f"{lane} missing executable has the wrong reason")
    elif reason == "interrupt":
        expected_detail = "wrapper signal interrupted the command after bounded process-group cleanup"
        require(rc in {129, 130, 143},
                f"{lane} interrupt result has an unsupported wrapper exit code")
        require(result["interrupted_signal"] == rc - 128,
                f"{lane} interrupt signal disagrees with the wrapper exit code")
    elif reason == "leader-exit":
        expected_detail = f"command exited {rc}"
        require(rc != 0 and result["leader_exit_code"] == rc,
                f"{lane} leader exit disagrees with the wrapper exit code")
    elif reason == "leader-signal":
        expected_detail = f"command exited {rc}"
        require(129 <= rc <= 255 and result["leader_exit_code"] == -(rc - 128),
                f"{lane} leader signal disagrees with the wrapper exit code")
    else:
        raise SystemExit(
            "proof-bundle verification failed: "
            f"{lane} command exit has unknown supervisor reason {reason!r}"
        )
    if reason not in {"supervisor-integrity-failure", "inspection-failure"}:
        require(
            result["process_group_drained"]
            and result["process_session_drained"]
            and result["output_pipe_eof"]
            and result["metadata_complete"]
            and not result["python_cap_reduced"]
            and result["inspection_error_count"] == 0,
            f"{lane} admitted supervisor disposition lacks complete containment metadata",
        )
    if reason in {
        "none",
        "output-truncated",
        "leader-exit-with-live-group",
        "timeout",
        "interrupt",
        "leader-exit",
        "leader-signal",
    }:
        require(type(result["leader_exit_code"]) is int,
                f"{lane} launched supervisor disposition lacks a leader exit code")
    if reason in {
        "none",
        "retained-log-budget-before-launch",
        "aggregate-deadline-before-launch",
        "launch-error",
        "launch-not-found",
        "leader-exit",
        "leader-signal",
    }:
        require(not result["output_truncated"],
                f"{lane} supervisor disposition spuriously claims output truncation")
    if reason != "interrupt":
        require(result["interrupted_signal"] is None,
                f"{lane} non-interrupt disposition carries an interrupt signal")
    supervisor_started_ns = result["supervisor_started_monotonic_ns"]
    run_deadline_ns = result["run_deadline_monotonic_ns"]
    if reason == "aggregate-deadline-before-launch":
        require(supervisor_started_ns >= run_deadline_ns,
                f"{lane} aggregate prelaunch refusal precedes its run deadline")
    elif reason not in {
        "retained-log-budget-before-launch",
        "supervisor-integrity-failure",
        "inspection-failure",
    }:
        require(supervisor_started_ns < run_deadline_ns,
                f"{lane} launched disposition starts at or after its run deadline")
    if reason in {
        "retained-log-budget-before-launch",
        "aggregate-deadline-before-launch",
        "launch-error",
        "launch-not-found",
    }:
        require(
            result["leader_exit_code"] is None
            and result["process_group_drained"]
            and result["process_session_drained"]
            and result["output_pipe_eof"]
            and result["metadata_complete"]
            and result["interrupted_signal"] is None
            and not result["term_sent"]
            and not result["kill_sent"],
            f"{lane} prelaunch supervisor result is contradictory",
        )
    require(
        lane_records[lane]["detail"] == expected_detail,
        f"{lane} verdict detail contradicts its canonical supervisor result",
    )
    if lane_records[lane]["status"] == "PASS":
        require(rc == 0, f"{lane} PASS verdict contradicts a nonzero supervisor exit")
    else:
        require(rc != 0, f"{lane} non-PASS verdict contradicts a zero supervisor exit")

def require_normal_lane_metadata(base_lanes):
    for lane in base_lanes:
        record = lane_records[lane]
        authority = expected_authority(lane)
        require(authority is not None, f"normal lane {lane} has no authority specification")
        require(record["authority"] == authority, f"normal lane {lane} has the wrong authority")
        require(record["log"] == f"logs/{lane}.log", f"normal lane {lane} has the wrong log path")
        if lane == "proof-boundary":
            require(
                record["detail"] == EXPECTED_PROOF_BOUNDARY_DETAIL,
                "proof-boundary has the wrong exact no-claim detail",
            )
            require(
                lane_log_data[lane] == EXPECTED_PROOF_BOUNDARY_LOG,
                "proof-boundary log does not match the exact software-only declaration",
            )
        elif lane == "source-manifest":
            if profile == "focused":
                expected_detail = EXPECTED_FOCUSED_SOURCE_MANIFEST_DETAIL
                expected_log = EXPECTED_FOCUSED_SOURCE_MANIFEST_LOG
            else:
                expected_detail = EXPECTED_CLOSURE_SOURCE_MANIFEST_DETAIL
                expected_log = EXPECTED_CLOSURE_SOURCE_MANIFEST_LOG
            require(record["detail"] == expected_detail,
                    f"{profile} source-manifest has the wrong exact detail")
            require(lane_log_data[lane] == expected_log,
                    f"{profile} source-manifest has the wrong exact body")
        elif lane in command_lanes:
            if record["status"] == "PASS":
                require(record["detail"] == "command completed", f"PASS lane {lane} has wrong detail")
            else:
                require(
                    record["detail"] in command_failure_details
                    or re.fullmatch(r"command exited (?:[1-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])", record["detail"]),
                    f"non-PASS lane {lane} has an unknown command disposition",
                )
            require_command_metadata(lane)

def require_normal_producer_state(base_lanes):
    require(executor in normal_executors, "normal proof state has a non-normal executor")
    require_profile_scope()
    allowed = set(focused_sequence if profile == "focused" else closure_sequence)
    if profile == "closure":
        allowed.add("source-manifest")
    require(set(base_lanes) <= allowed, "terminal state contains lanes outside its profile's closed universe")
    producer_sequence = focused_sequence if profile == "focused" else closure_sequence
    if profile == "closure" and "closure-root-preflight" in base_lanes:
        preflight_status = lane_statuses["closure-root-preflight"]
        require(preflight_status in {"PASS", "NO_DATA", "FAIL"}, "closure-root-preflight has an impossible status")
        if preflight_status == "NO_DATA":
            producer_sequence = closure_refusal_sequence
    require_trace_prefix(base_lanes, producer_sequence, "normal producer trace")
    for lane in base_lanes:
        if lane == "proof-boundary":
            require(lane_statuses[lane] == "PASS", "proof-boundary must be PASS")
        elif lane == "source-manifest":
            require(lane_statuses[lane] == "NO_DATA", "source-manifest must be NO_DATA")
        elif lane == "closure-root-preflight":
            require(lane_statuses[lane] in {"PASS", "NO_DATA", "FAIL"},
                    "closure-root-preflight must be PASS, NO_DATA, or FAIL")
        else:
            require(lane_statuses[lane] in {"PASS", "FAIL"},
                    f"normal producer lane {lane} must be PASS or FAIL")
    for stop_lane in ("constellation-verify", "constellation-snapshot-before", "constellation-snapshot-after"):
        if stop_lane in base_lanes and lane_statuses[stop_lane] == "FAIL":
            require(base_lanes[-1] == stop_lane,
                    f"normal producer trace continues after failed control lane {stop_lane}")
    if profile == "closure" and "closure-root-preflight" in base_lanes \
            and lane_statuses["closure-root-preflight"] == "NO_DATA":
        require(lane_records["closure-root-preflight"]["detail"] == "command exited 1",
                "closure refusal requires the exact policy-refusal exit")
        require(base_lanes[-1] in {"closure-root-preflight", "source-manifest"},
                "closure refusal trace continues into the normal closure branch")
    if profile == "closure" and "closure-root-preflight" in base_lanes \
            and lane_statuses["closure-root-preflight"] == "FAIL":
        require(base_lanes[-1] == "closure-root-preflight",
                "closure trace continues after an infrastructure preflight failure")
    before_lane = "constellation-snapshot-before"
    after_lane = "constellation-snapshot-after"
    if before_lane not in base_lanes:
        require(not before_snapshot and not after_snapshot,
                "source snapshots exist before snapshot-before was attempted")
    else:
        if lane_statuses[before_lane] == "PASS":
            require(bool(before_snapshot), "successful snapshot-before lacks retained bytes")
        else:
            require(not before_snapshot, "failed snapshot-before retained authoritative bytes")
        if after_lane not in base_lanes:
            require(not after_snapshot, "snapshot-after exists before that lane was attempted")
    if after_lane in base_lanes and lane_statuses[after_lane] == "PASS":
        require(bool(after_snapshot), "successful snapshot-after lacks retained bytes")
    if after_lane in base_lanes and lane_statuses[after_lane] != "PASS":
        require(not after_snapshot, "failed snapshot-after retained authoritative bytes")
    require(not after_snapshot or bool(before_snapshot),
            "after snapshot cannot exist without a before snapshot")
    require_normal_lane_metadata(base_lanes)

def require_trace_prefix(trace, sequence, context):
    require(
        tuple(trace) == tuple(sequence[:len(trace)]),
        f"{context} is not an ordered producer-trace prefix",
    )

if ready:
    require(status == "READY_FOR_DSR" and exit_code == 0 and failed_lanes == 0,
            "DSR-ready state is inconsistent")
    require(summary.get("source_manifest_coverage") == "full-root-clean-head-bookended",
            "DSR-ready state lacks full-root bookended coverage")
    require(lane_names == closure_required,
            "DSR-ready bundle does not contain the exact closure lane set")
    require(lane_statuses == closure_statuses,
            "DSR-ready bundle does not contain the exact closure status matrix")
    require(tuple(lane_order) == closure_sequence,
            "DSR-ready bundle does not contain the exact ordered closure trace")
    require(summary.get("profile") == "closure" and
            executor in normal_executors and
            summary.get("proof_scope") == "head-bookended-closure-candidate" and
            summary.get("provenance_state") == "stable",
            "DSR-ready bundle has inconsistent profile, proof scope, or provenance")
    require(summary.get("snapshot_before_sha256") and
            summary.get("snapshot_before_sha256") == summary.get("snapshot_after_sha256"),
            "DSR-ready bundle lacks equal nonempty source snapshots")
    require_normal_producer_state(lane_order)
if status == "FOCUSED_PASS":
    require(exit_code == 0 and failed_lanes == 0 and not ready,
            "focused-pass state is inconsistent")
    require(lane_names == focused_required,
            "focused-pass bundle does not contain the exact focused lane set")
    require(lane_statuses == focused_statuses,
            "focused-pass bundle does not contain the exact focused status matrix")
    require(tuple(lane_order) == focused_sequence,
            "focused-pass bundle does not contain the exact ordered focused trace")
    require(summary.get("profile") == "focused" and
            executor in normal_executors and
            summary.get("proof_scope") == "focused-software-only" and
            summary.get("provenance_state") == "stable" and
            summary.get("source_manifest_coverage") == "not-checked",
            "focused-pass bundle has inconsistent profile, scope, provenance, or coverage")
    require(summary.get("snapshot_before_sha256") and
            summary.get("snapshot_before_sha256") == summary.get("snapshot_after_sha256"),
            "focused-pass bundle lacks equal nonempty source snapshots")
    require_normal_producer_state(lane_order)
self_test_log_label_overrides = {
    "self-test-special-index-flag-refusal": "special-index-flag",
    "self-test-skip-worktree-flag-refusal": "skip-worktree-flag",
    "self-test-fsmonitor-valid-flag-refusal": "fsmonitor-valid-flag",
    "self-test-zero-smoke-refusal": "zero-smoke-sentinel",
    "self-test-wrong-authority-refusal": "wrong-authority",
    "self-test-wrong-command-refusal": "wrong-command",
    "self-test-supervisor-state-contradiction-refusal": "supervisor-state-contradiction-refusal",
    "self-test-proof-boundary-body-refusal": "proof-boundary-body-refusal",
    "self-test-source-manifest-body-refusal": "source-manifest-body-refusal",
    "self-test-failed-control-continuation-refusal": "failed-control-continuation",
    "self-test-premature-snapshot-refusal": "premature-snapshot",
    "self-test-toctou-mutation-refusal": "toctou-mutation",
    "self-test-publication-destination-race-refusal": "publication-destination-race-refusal",
    "self-test-mutated-bundle-refusal": "mutated-bundle",
    "self-test-duplicate-seal-refusal": "duplicate-seal",
    "self-test-truncated-seal-refusal": "truncated-seal",
    "self-test-summary-mismatch-refusal": "summary-mismatch",
    "self-test-unsafe-path-refusal": "unsafe-path",
    "self-test-duplicate-json-key-refusal": "duplicate-json-key",
    "self-test-nonfinite-json-refusal": "nonfinite-json",
    "self-test-oversized-json-refusal": "oversized-json",
    "self-test-record-count-refusal": "record-count",
    "self-test-readiness-mismatch-refusal": "readiness-mismatch",
    "self-test-no-claim-mutation-refusal": "no-claim-mutation",
    "self-test-unknown-terminal-refusal": "unknown-terminal",
    "self-test-unreferenced-file-refusal": "unexpected-file",
    "self-test-prefix-file-mismatch-refusal": "prefix-mismatch",
    "self-test-interrupted-status-matrix-refusal": "invalid-interrupted-matrix",
    "self-test-incomplete-status-matrix-refusal": "invalid-incomplete-matrix",
}
if status == "SELF_TEST_PASS":
    require(exit_code == 0 and failed_lanes == 0 and not ready,
            "self-test-pass state is inconsistent")
    require(lane_names == self_test_required,
            "self-test-pass bundle does not contain the exact self-test lane set")
    require(lane_statuses == self_test_statuses,
            "self-test-pass bundle does not contain the exact self-test status matrix")
    require(tuple(lane_order) == self_test_sequence,
            "self-test-pass bundle does not contain the exact ordered self-test trace")
    require(summary.get("profile") == "focused" and
            summary.get("executor_declaration") == "self-test-no-cargo" and
            summary.get("proof_scope") == "harness-self-test-no-cargo" and
            summary.get("provenance_state") == "stable" and
            summary.get("source_manifest_coverage") == "not-applicable" and
            not summary.get("snapshot_before_sha256") and
            not summary.get("snapshot_after_sha256"),
            "self-test-pass bundle has inconsistent scope or provenance")
    for lane in self_test_sequence:
        record = lane_records[lane]
        require(record["authority"] == "harness-self-test",
                f"SELF_TEST_PASS lane {lane} has the wrong authority")
        require(
            record["detail"] == "harness self-test assertion matched expected disposition",
            f"SELF_TEST_PASS lane {lane} has the wrong detail",
        )
        expected_label = self_test_log_label_overrides.get(
            lane, lane.removeprefix("self-test-")
        )
        require(record["log"] == f"logs/self-test-{expected_label}.log",
                f"SELF_TEST_PASS lane {lane} has the wrong log locator")
        require(
            lane_log_data[lane].startswith(
                f"self-test-label={expected_label}\n".encode()
            ),
            f"SELF_TEST_PASS lane {lane} log is not bound to its lane label",
        )
        assertion_results = re.findall(
            rb"(?m)^self-test-assertion-result=(pass|fail)$",
            lane_log_data[lane],
        )
        require(
            assertion_results == [b"pass"]
            and lane_log_data[lane].endswith(b"self-test-assertion-result=pass\n"),
            f"SELF_TEST_PASS lane {lane} lacks one final passing assertion result",
        )
if status == "NO_DATA" and summary.get("profile") == "closure":
    require(exit_code == 4 and not ready and
            executor in normal_executors and
            lane_names == {"proof-boundary", "closure-root-preflight", "source-manifest"},
            "closure NO_DATA bundle is inconsistent")
    require(
        tuple(lane_order) == (
            "proof-boundary", "closure-root-preflight", "source-manifest"
        ),
        "closure NO_DATA bundle has the wrong ordered producer trace",
    )
    require(lane_statuses == {
                "proof-boundary": "PASS",
                "closure-root-preflight": "NO_DATA",
                "source-manifest": "NO_DATA",
            } and
            summary.get("proof_scope") == "head-bookended-closure-candidate" and
            summary.get("provenance_state") == "incomplete" and
            summary.get("source_manifest_coverage") == "not-checked" and
            not summary.get("snapshot_before_sha256") and
            not summary.get("snapshot_after_sha256"),
            "closure NO_DATA bundle has inconsistent statuses or provenance")
    require_normal_producer_state(lane_order)
if status == "NO_DATA":
    require(summary.get("profile") == "closure",
            "NO_DATA is only a closed terminal state for closure preflight refusal")
if status == "FAIL":
    require(exit_code in {1, 7} and failed_lanes > 0 and not ready,
            "FAIL terminal state is inconsistent")
    require_normal_producer_state(lane_order)
    producer_sequence = focused_sequence if profile == "focused" else closure_sequence
    require_trace_prefix(lane_order, producer_sequence, "FAIL trace")
    require_named_lane_status("proof-boundary", "PASS")
    if profile == "focused":
        require_named_lane_status("source-manifest", "NO_DATA")
    else:
        require(
            lane_statuses["closure-root-preflight"] in {"PASS", "FAIL"},
            "closure FAIL trace has an impossible preflight status",
        )
    for lane in lane_order:
        if lane not in {"proof-boundary", "source-manifest", "closure-root-preflight"}:
            require(
                lane_statuses[lane] in {"PASS", "FAIL"},
                f"FAIL trace lane {lane} has an impossible status",
            )
    early_traces = {
        tuple(producer_sequence[:3]),
        tuple(producer_sequence[:4]),
    }
    if profile == "closure":
        early_traces.add(tuple(producer_sequence[:2]))
    snapshot_after_index = producer_sequence.index("constellation-snapshot-after")
    snapshot_after_trace = tuple(producer_sequence[:snapshot_after_index + 1])
    full_trace = tuple(producer_sequence)
    observed_trace = tuple(lane_order)
    if exit_code == 7:
        require(observed_trace in early_traces, "exit-7 FAIL has an impossible trace")
        require(lane_statuses[lane_order[-1]] == "FAIL",
                "exit-7 FAIL trace does not end in failure")
        if profile == "closure" and len(lane_order) > 2:
            require_named_lane_status("closure-root-preflight", "PASS")
        if len(lane_order) == 4:
            require_named_lane_status("constellation-verify", "PASS")
        require(coverage == "not-checked" and provenance == "incomplete",
                "exit-7 FAIL has impossible coverage or provenance")
    else:
        require(
            observed_trace in {snapshot_after_trace, full_trace},
            "exit-1 FAIL is neither the snapshot-failure nor full producer trace",
        )
        require_named_lane_status("constellation-verify", "PASS")
        require_named_lane_status("constellation-snapshot-before", "PASS")
        if observed_trace == snapshot_after_trace:
            require_named_lane_status("constellation-snapshot-after", "FAIL")
            require(provenance == "incomplete",
                    "snapshot-after FAIL must have incomplete provenance")
        else:
            require_named_lane_status("constellation-snapshot-after", "PASS")
            moved = lane_statuses["source-stability"] == "FAIL"
            if profile == "closure":
                moved = moved or lane_statuses["closure-root-bookend"] == "FAIL"
            require(
                provenance == ("moved" if moved else "stable"),
                "full FAIL trace has provenance inconsistent with its bookends",
            )
        require(
            coverage == ("not-checked" if profile == "focused" else "unestablished"),
            "exit-1 FAIL has impossible manifest coverage",
        )
    if provenance == "stable":
        require(
            before_snapshot and before_snapshot == after_snapshot,
            "stable FAIL lacks equal nonempty source snapshots",
        )
    elif provenance == "moved":
        require(
            before_snapshot and after_snapshot,
            "moved FAIL lacks both source snapshots",
        )
if status == "SELF_TEST_FAIL":
    standard_self_test_failure = (
        tuple(lane_order) == self_test_sequence and
        set(lane_statuses.values()) <= {"PASS", "FAIL"}
    )
    signal_target_failure = (
        tuple(lane_order) == ("signal-target-unexpected-completion",) and
        lane_statuses["signal-target-unexpected-completion"] == "FAIL"
    )
    require(exit_code == 1 and failed_lanes > 0 and not ready and
            summary.get("profile") == "focused" and
            summary.get("executor_declaration") == "self-test-no-cargo" and
            summary.get("proof_scope") == "harness-self-test-no-cargo" and
            summary.get("source_manifest_coverage") == "not-applicable" and
            provenance == "incomplete" and
            not before_snapshot and not after_snapshot and
            (standard_self_test_failure or signal_target_failure),
            "SELF_TEST_FAIL terminal state is inconsistent")
    if standard_self_test_failure:
        for lane in self_test_sequence:
            record = lane_records[lane]
            expected_label = self_test_log_label_overrides.get(
                lane, lane.removeprefix("self-test-")
            )
            require(record["authority"] == "harness-self-test",
                    f"SELF_TEST_FAIL lane {lane} has the wrong authority")
            require(record["log"] == f"logs/self-test-{expected_label}.log",
                    f"SELF_TEST_FAIL lane {lane} has the wrong log locator")
            require(
                lane_log_data[lane].startswith(
                    f"self-test-label={expected_label}\n".encode()
                ),
                f"SELF_TEST_FAIL lane {lane} log is not bound to its lane label",
            )
            assertion_results = re.findall(
                rb"(?m)^self-test-assertion-result=(pass|fail)$",
                lane_log_data[lane],
            )
            expected_result = b"pass" if record["status"] == "PASS" else b"fail"
            expected_detail = (
                "harness self-test assertion matched expected disposition"
                if record["status"] == "PASS"
                else "harness self-test assertion did not match expected disposition"
            )
            require(record["detail"] == expected_detail,
                    f"SELF_TEST_FAIL lane {lane} has a contradictory detail")
            require(
                assertion_results == [expected_result]
                and lane_log_data[lane].endswith(
                    b"self-test-assertion-result=" + expected_result + b"\n"
                ),
                f"SELF_TEST_FAIL lane {lane} has a contradictory assertion result",
            )
if status == "INTERRUPTED":
    require(exit_code in {129, 130, 143} and failed_lanes > 0 and not ready and
            coverage == "not-checked" and provenance == "incomplete",
            "INTERRUPTED terminal state is inconsistent")
    require_named_lane_status("wrapper-signal", "FAIL")
    require(lane_order[-1] == "wrapper-signal",
            "wrapper-signal must be the final nonterminal lane")
    require_wrapper_signal_record()
    interrupted_prefix = lane_order[:-1]
    if "internal-incomplete" in lane_names:
        require(len(interrupted_prefix) > 0 and
                interrupted_prefix[-1] == "internal-incomplete",
                "internal-incomplete must immediately precede wrapper-signal")
        require_named_lane_status("internal-incomplete", "FAIL")
        require_internal_incomplete_record()
        interrupted_prefix = interrupted_prefix[:-1]
    if executor == "self-test-no-cargo":
        require(profile == "focused" and proof_scope == "focused-software-only",
                "self-test interruption has inconsistent profile or scope")
        require(
            lane_names <= self_test_required | {
                "self-test-signal-target-wait",
                "internal-incomplete",
                "wrapper-signal",
            },
            "self-test interruption contains an impossible lane",
        )
        require(not before_snapshot and not after_snapshot,
                "self-test interruption cannot carry source snapshots")
        if interrupted_prefix == ["self-test-signal-target-wait"]:
            require_named_lane_status("self-test-signal-target-wait", "FAIL")
        else:
            require_trace_prefix(
                interrupted_prefix, self_test_sequence, "self-test interruption trace"
            )
            for lane in interrupted_prefix:
                require(lane_statuses[lane] in {"PASS", "FAIL"},
                        f"self-test interruption lane {lane} must be PASS or FAIL")
    else:
        require_normal_producer_state(interrupted_prefix)
if status == "INCOMPLETE":
    require(exit_code != 0 and failed_lanes > 0 and not ready and
            coverage == "not-checked" and provenance == "incomplete",
            "INCOMPLETE terminal state is inconsistent")
    require_named_lane_status("internal-incomplete", "FAIL")
    require(lane_order[-1] == "internal-incomplete",
            "internal-incomplete must be the final nonterminal lane")
    require_internal_incomplete_record(exit_code)
    incomplete_prefix = lane_order[:-1]
    if executor == "self-test-no-cargo":
        require(profile == "focused" and proof_scope == "focused-software-only",
                "self-test incomplete state has inconsistent profile or scope")
        require(
            lane_names <= self_test_required | {"internal-incomplete"},
            "self-test incomplete state contains an impossible lane",
        )
        require(not before_snapshot and not after_snapshot,
                "self-test incomplete state cannot carry source snapshots")
        require_trace_prefix(
            incomplete_prefix, self_test_sequence, "self-test incomplete trace"
        )
        for lane in incomplete_prefix:
            require(lane_statuses[lane] in {"PASS", "FAIL"},
                    f"self-test incomplete lane {lane} must be PASS or FAIL")
    else:
        require_normal_producer_state(incomplete_prefix)
if exit_code == 0:
    require(status in {"FOCUSED_PASS", "READY_FOR_DSR", "SELF_TEST_PASS"},
            "zero exit has non-success status")

# This marker and delay create a deterministic observation window for an
# external hostile self-test. The verifier itself remains read-only with
# respect to the bundle under verification.
if os.environ.get("FSIM_EULER_DISC_E2E_SELF_TEST_PAUSE_DURING_VERIFY") == "1":
    marker_text = os.environ.get("FSIM_EULER_DISC_E2E_SELF_TEST_VERIFY_READY_MARKER", "")
    require(bool(marker_text), "verify-pause self-test lacks a readiness marker")
    marker = pathlib.Path(marker_text)
    marker.parent.mkdir(parents=True, exist_ok=True)
    try:
        with marker.open("x", encoding="utf-8") as handle:
            handle.write("semantic-verification-complete=true\n")
            handle.flush()
            os.fsync(handle.fileno())
    except OSError as error:
        require(False, f"cannot publish verify-pause readiness marker: {error}")
    time.sleep(1)

middle_inventory = inventory_bundle()
require(middle_inventory == initial_inventory, "bundle inventory or entry metadata changed during verification")
require(
    metadata_tuple(os.fstat(root_descriptor)) == initial_root_metadata,
    "bundle root metadata changed during verification",
)
for relative_name in sorted(expected_files):
    if relative_name == "summary.json":
        limit = MAX_SUMMARY_BYTES
    elif relative_name in {"verdicts.jsonl", "verdicts-prefix.jsonl"}:
        limit = MAX_VERDICTS_BYTES
    elif relative_name in {"snapshot-before.txt", "snapshot-after.txt"}:
        limit = MAX_SNAPSHOT_BYTES
    else:
        limit = MAX_LOG_BYTES
    final_data = bounded_regular_bytes(
        relative_name,
        limit,
        f"final stable read of {relative_name}",
        remember=False,
    )
    require(
        hashlib.sha256(final_data).digest() == verified_file_digests[relative_name],
        f"{relative_name} changed after its semantic verification",
    )
final_inventory = inventory_bundle()
require(final_inventory == initial_inventory, "bundle inventory or entry metadata changed during final verification")
require(
    metadata_tuple(os.fstat(root_descriptor)) == initial_root_metadata,
    "bundle root metadata changed during final verification",
)
try:
    final_root_lstat = root_argument.lstat()
except OSError as error:
    raise SystemExit(
        f"proof-bundle verification failed: bundle root path disappeared: {error}"
    ) from error
require(
    stat.S_ISDIR(final_root_lstat.st_mode)
    and (final_root_lstat.st_dev, final_root_lstat.st_ino)
    == (root_opened.st_dev, root_opened.st_ino),
    "bundle root path changed identity during verification",
)

binding_rows = [
    ["directory", relative_name]
    for relative_name in sorted(expected_directories)
]
binding_rows.extend(
    [
        "file",
        relative_name,
        initial_inventory[relative_name][1][4],
        verified_file_digests[relative_name].hex(),
    ]
    for relative_name in sorted(expected_files)
)
binding_rows.sort(key=lambda row: (row[1], row[0]))
binding_bytes = json.dumps(
    binding_rows,
    ensure_ascii=True,
    separators=(",", ":"),
    allow_nan=False,
).encode("ascii")
binding_sha256 = hashlib.sha256(binding_bytes).hexdigest()
print(f"proof-bundle verified: {root}")
print(
    "proof-bundle-binding: "
    f"dev={root_opened.st_dev} ino={root_opened.st_ino} "
    f"commitment_sha256={binding_sha256}"
)
PY
}

if [[ "$VERIFY_BUNDLE_SET" == 1 && "$SELF_TEST" == 1 ]]; then
  printf '%s\n' '--verify-bundle and --self-test are mutually exclusive' >&2
  exit 2
fi

if [[ "$VERIFY_BUNDLE_SET" == 1 ]]; then
  verify_proof_bundle "$VERIFY_BUNDLE"
  exit $?
fi

if [[ "$SELF_TEST" == 1 ]]; then
  EXECUTOR_DECLARATION="self-test-no-cargo"
else
  EXECUTOR_DECLARATION="${FSIM_EULER_DISC_E2E_EXECUTOR:-local}"
fi
case "$EXECUTOR_DECLARATION" in
  rch)
    if [[ -z "${CARGO_TARGET_DIR:-}" ]]; then
      printf '%s\n' \
        'rch executor declaration requires an explicit CARGO_TARGET_DIR' >&2
      exit 2
    fi
    ;;
  dsr|self-test-no-cargo) ;;
  local)
    if [[ "${FSIM_EULER_DISC_E2E_ALLOW_LOCAL:-0}" != 1 ]]; then
      printf '%s\n' \
        'local Cargo execution requires FSIM_EULER_DISC_E2E_ALLOW_LOCAL=1' >&2
      exit 2
    fi
    ;;
  *)
    printf 'invalid executor declaration: %s\n' "$EXECUTOR_DECLARATION" >&2
    exit 2
    ;;
esac

LANE_TIMEOUT_SECONDS="${FSIM_EULER_DISC_E2E_LANE_TIMEOUT_SECONDS:-3600}"
RUN_TIMEOUT_SECONDS="${FSIM_EULER_DISC_E2E_RUN_TIMEOUT_SECONDS:-14400}"
LANE_LOG_MAX_BYTES="${FSIM_EULER_DISC_E2E_LANE_LOG_MAX_BYTES:-16777216}"
HARD_MAX_RETAINED_LOG_BYTES="$PROOF_MAX_LANE_LOG_BYTES"
RETAINED_LOG_NONCHILD_RESERVE_BYTES=$((64 * 1024))
HARD_MAX_LANE_LOG_BYTES=$((HARD_MAX_RETAINED_LOG_BYTES - RETAINED_LOG_NONCHILD_RESERVE_BYTES))

# shellcheck disable=SC2329 # Exported for the bounded self-test child.
lane_log_cap_is_valid() { # child-output-byte-cap
  [[ "$1" =~ ^[1-9][0-9]{0,7}$ ]] && ((10#$1 <= HARD_MAX_LANE_LOG_BYTES))
}
if ! parsed_numeric_settings="$(
  python3 - \
    "$LANE_TIMEOUT_SECONDS" "$RUN_TIMEOUT_SECONDS" "$LANE_LOG_MAX_BYTES" \
    "$PROOF_MAX_TIMEOUT_SECONDS" "$HARD_MAX_LANE_LOG_BYTES" <<'PY'
import re
import sys

lane_text, run_text, log_text, timeout_max_text, log_max_text = sys.argv[1:]
timeout_max = int(timeout_max_text)
log_max = int(log_max_text)

def bounded_positive(text, maximum, label):
    maximum_text = str(maximum)
    if (
        re.fullmatch(r"[1-9][0-9]*", text) is None
        or len(text) > len(maximum_text)
        or (len(text) == len(maximum_text) and text > maximum_text)
    ):
        raise SystemExit(f"{label} must be a canonical positive integer <= {maximum}")
    return int(text)

lane = bounded_positive(lane_text, timeout_max, "lane timeout")
run = bounded_positive(run_text, timeout_max, "run timeout")
log = bounded_positive(log_text, log_max, "lane child-output cap")
print(lane, run, log)
PY
)"; then
  printf '%s\n' 'invalid bounded Euler-disc harness numeric configuration' >&2
  exit 2
fi
read -r LANE_TIMEOUT_SECONDS RUN_TIMEOUT_SECONDS LANE_LOG_MAX_BYTES \
  <<<"$parsed_numeric_settings"
read -r RUN_STARTED_MONOTONIC_NS RUN_DEADLINE_MONOTONIC_NS < <(
  python3 - "$RUN_TIMEOUT_SECONDS" <<'PY'
import sys
import time

started = time.monotonic_ns()
print(started, started + int(sys.argv[1]) * 1_000_000_000)
PY
)

absolute_path() {
  case "$1" in
    /*) printf '%s\n' "$1" ;;
    *) printf '%s\n' "$REPO_ROOT/$1" ;;
  esac
}

HEAD_SHA="$(git rev-parse HEAD)"
HOST_ISA="$(uname -m)"
if [[ "$SELF_TEST" == 1 ]]; then
  # The no-Cargo self-test must remain runnable on a host with no Cargo binary.
  CARGO_BIN="not-used-by-harness-self-test"
else
  CARGO_BIN="${FSIM_EULER_DISC_E2E_CARGO:-$(command -v cargo || true)}"
  if [[ -z "$CARGO_BIN" ]]; then
    printf '%s\n' 'Cargo is required for focused or closure execution' >&2
    exit 2
  fi
fi
LOG_ROOT="$(absolute_path "${FSIM_EULER_DISC_E2E_LOG_DIR:-target/euler-disc-contract-e2e}")"
mkdir -p "$LOG_ROOT"
LOG_DIR="$(mktemp -d "$LOG_ROOT/.candidate-${HEAD_SHA:0:12}-${PROFILE}-XXXXXXXX")"
PUBLISHED_DIR="$LOG_ROOT/${LOG_DIR##*/.candidate-}"
mkdir -p "$LOG_DIR/logs"
VERDICTS="$LOG_DIR/verdicts-prefix.jsonl"
SUMMARY="$LOG_DIR/summary.json"
: >"$VERDICTS"

FAILURES=0
CHECKS=0
LAST_LANE_RC=0
LAST_LANE_STATUS=PASS
SEALED=0
SEALING=0
LANE_COMMIT_CRITICAL=0
SNAPSHOT_BEFORE=""
SNAPSHOT_AFTER=""
ACTIVE_SUPERVISOR_PID=""
WRAPPER_SIGNAL=0
CONSUMED_WRAPPER_SIGNAL=0
SELF_TEST_ACTIVE_SUPERVISOR_MARKER="${FSIM_EULER_DISC_E2E_SELF_TEST_ACTIVE_SUPERVISOR_MARKER:-}"
PENDING_SELF_TEST_ACTIVE_SUPERVISOR_MARKER="$SELF_TEST_ACTIVE_SUPERVISOR_MARKER"
SELF_TEST_FINALIZER_READY_MARKER=""
SELF_TEST_LAST_OUTPUT_OFFSET=0
unset FSIM_EULER_DISC_E2E_SELF_TEST_ACTIVE_SUPERVISOR_MARKER
RETAINED_LOG_BYTES_TOTAL=0
PREPUBLICATION_SIGNAL_INJECTED=0
POSTPUBLICATION_SIGNAL_INJECTED=0
FINAL_COMMIT_SIGNAL_INJECTED=0
CANDIDATE_CORRUPTION_INJECTED=0
PREPUBLICATION_DEADLINE_INJECTED=0
POSTPUBLICATION_DEADLINE_INJECTED=0
POSTPUBLICATION_VERIFY_READY_MARKER=""
CANDIDATE_REJECTED=0
RETRACTION_LAST_LOG=""
RETRACTION_INTEGRITY_RECOVERY=0
RETRACTION_INTEGRITY_CONTEXT=""
FINALIZATION_REPORT_CURSOR=0
declare -a FINALIZATION_REPORT_LABELS=()
declare -a FINALIZATION_REPORT_RCS=()
declare -a FINALIZATION_REPORT_LOGS=()

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

record() { # lane status authority detail run-relative-log
  local lane="$1" status="$2" authority="$3" detail="$4" log_rel="$5"
  local log_path="$LOG_DIR/$log_rel"
  local log_bytes=0 log_sha=""
  if [[ -f "$log_path" ]]; then
    log_bytes="$(wc -c <"$log_path" | tr -d ' ')"
    log_sha="$(sha256_file "$log_path")"
  fi
  if ((log_bytes > PROOF_MAX_LANE_LOG_BYTES)); then
    printf 'refusing oversized retained lane log: lane=%s bytes=%s limit=%s\n' \
      "$lane" "$log_bytes" "$PROOF_MAX_LANE_LOG_BYTES" >&2
    return 125
  fi
  if ((RETAINED_LOG_BYTES_TOTAL + log_bytes > PROOF_MAX_TOTAL_LOG_BYTES)); then
    printf 'refusing aggregate retained-log overflow: lane=%s prior=%s bytes=%s limit=%s\n' \
      "$lane" "$RETAINED_LOG_BYTES_TOTAL" "$log_bytes" \
      "$PROOF_MAX_TOTAL_LOG_BYTES" >&2
    return 125
  fi
  python3 - "$lane" "$status" "$authority" "$detail" "$log_rel" \
    "$log_bytes" "$log_sha" "$HEAD_SHA" "$HOST_ISA" "$PROFILE" \
    "$EXECUTOR_DECLARATION" >>"$VERDICTS" <<'PY'
import json
import sys

(
    lane,
    status,
    authority,
    detail,
    log_path,
    log_bytes,
    log_sha,
    head,
    isa,
    profile,
    executor_declaration,
) = sys.argv[1:]
print(json.dumps({
    "schema": "frankensim.euler-disc-contract-e2e.verdict.v1",
    "lane": lane,
    "status": status,
    "authority": authority,
    "detail": detail,
    "log": log_path,
    "log_bytes": int(log_bytes),
    "log_sha256": log_sha,
    "head": head,
    "host_isa": isa,
    "profile": profile,
    "executor_declaration": executor_declaration,
    "executor_attestation": "caller-declared-unverified",
    "provenance_state": "provisional",
    "terminal": False,
}, sort_keys=True, separators=(",", ":"), allow_nan=False))
PY
  RETAINED_LOG_BYTES_TOTAL=$((RETAINED_LOG_BYTES_TOTAL + log_bytes))
  CHECKS=$((CHECKS + 1))
  if [[ "$status" == "FAIL" ]]; then
    FAILURES=$((FAILURES + 1))
  fi
  printf '[%-7s] %-28s %s\n' "$status" "$lane" "$detail" >&2
}

run_bounded_command() { # log lane-timeout run-deadline configured-cap shell-cap total-cap prelaunch-refusal command...
  local log_path="$1" timeout_seconds="$2" run_deadline_ns="$3"
  local configured_output_cap="$4" shell_output_cap="$5" max_total_bytes="$6"
  local prelaunch_refusal="$7"
  local supervisor_pid status finalizer_readiness_failed=0 forwarding_exercised=0
  local active_readiness_marker="$PENDING_SELF_TEST_ACTIVE_SUPERVISOR_MARKER"
  local finalizer_marker="$SELF_TEST_FINALIZER_READY_MARKER"
  local finalizer_release_marker=""
  local finalizer_ack_marker=""
  PENDING_SELF_TEST_ACTIVE_SUPERVISOR_MARKER=""
  SELF_TEST_FINALIZER_READY_MARKER=""
  shift 7
  # Python deliberately owns and appends the same retained log as its stderr.
  # shellcheck disable=SC2094
  python3 - "$log_path" "$timeout_seconds" "$run_deadline_ns" \
    "$configured_output_cap" "$shell_output_cap" "$max_total_bytes" \
    "$prelaunch_refusal" "$active_readiness_marker" "$finalizer_marker" \
    "$FINALIZER_HANDSHAKE_TIMEOUT_SECONDS" "$CONSUMED_WRAPPER_SIGNAL" "$@" \
    2>>"$log_path" <<'PY' &
import atexit
import json
import os
import pathlib
import selectors
import signal
import stat
import subprocess
import sys
import tempfile
import time

(
    log_path,
    timeout_text,
    run_deadline_text,
    configured_cap_text,
    shell_cap_text,
    total_cap_text,
    prelaunch_refusal,
    active_readiness_marker_text,
    finalizer_readiness_marker_text,
    finalizer_handshake_timeout_text,
    finalizer_signal_text,
    *command,
) = sys.argv[1:]
timeout_seconds = int(timeout_text)
run_deadline_ns = int(run_deadline_text)
configured_max_output_bytes = int(configured_cap_text)
shell_max_output_bytes = int(shell_cap_text)
max_total_bytes = int(total_cap_text)
finalizer_handshake_timeout_seconds = int(finalizer_handshake_timeout_text)
finalizer_signal_number = int(finalizer_signal_text)
metadata_reserve_bytes = 32 * 1024
initial_log_bytes = os.path.getsize(log_path)
max_output_bytes = max(
    0,
    min(
        shell_max_output_bytes,
        max_total_bytes - initial_log_bytes - metadata_reserve_bytes,
    ),
)
shell_cap_reduced = shell_max_output_bytes < configured_max_output_bytes
python_cap_reduced = max_output_bytes < shell_max_output_bytes
supervisor_started_ns = time.monotonic_ns()
lane_deadline_ns = supervisor_started_ns + timeout_seconds * 1_000_000_000
effective_deadline_ns = min(lane_deadline_ns, run_deadline_ns)
deadline_kind = "aggregate" if run_deadline_ns <= lane_deadline_ns else "lane"
written = 0
truncated = False
requested_signal = None
inspection_errors = []
observed_descendants = set()
shutdown_reason = None
interrupted_signal = None
term_sent_at_ns = None
drain_deadline_ns = None
term_signal_sent = False
kill_signal_sent = False
complete = False
eof = False
completion_delay_injected = False
completion_signal_injected = False

RESULT_SCHEMA = "frankensim.euler-disc-contract-e2e.supervisor-result.v1"

def canonical_result(
    *,
    wrapper_exit_code,
    reason,
    leader_exit_code=None,
    process_group_drained=False,
    process_session_drained=False,
    output_pipe_eof=False,
    metadata_complete=True,
):
    return {
        "configured_lane_timeout_seconds": timeout_seconds,
        "configured_output_cap_bytes": configured_max_output_bytes,
        "deadline_kind": deadline_kind,
        "effective_deadline_monotonic_ns": effective_deadline_ns,
        "initial_log_bytes": initial_log_bytes,
        "inspection_error_count": len(inspection_errors),
        "interrupted_signal": interrupted_signal,
        "kill_sent": kill_signal_sent,
        "lane_deadline_monotonic_ns": lane_deadline_ns,
        "leader_exit_code": leader_exit_code,
        "metadata_complete": metadata_complete,
        "output_pipe_eof": output_pipe_eof,
        "output_truncated": truncated,
        "process_group_drained": process_group_drained,
        "process_session_drained": process_session_drained,
        "python_cap_reduced": python_cap_reduced,
        "retained_output_cap_bytes": max_output_bytes,
        "run_deadline_monotonic_ns": run_deadline_ns,
        "schema": RESULT_SCHEMA,
        "shell_cap_reduced": shell_cap_reduced,
        "shell_effective_output_cap_bytes": shell_max_output_bytes,
        "shutdown_reason": reason,
        "supervisor_started_monotonic_ns": supervisor_started_ns,
        "term_sent": term_signal_sent,
        "total_log_cap_bytes": max_total_bytes,
        "wrapper_exit_code": wrapper_exit_code,
    }

def encode_result(**kwargs):
    payload = json.dumps(
        canonical_result(**kwargs),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )
    return f"supervisor_result_json={payload}\n".encode("ascii")

def append_prelaunch_result(log, *, wrapper_exit_code, reason, launch_detail=None):
    if launch_detail is not None:
        log.write(f"{launch_detail}\n".encode("utf-8", "backslashreplace"))
    result = encode_result(
        wrapper_exit_code=wrapper_exit_code,
        reason=reason,
        process_group_drained=True,
        process_session_drained=True,
        output_pipe_eof=True,
    )
    if log.tell() + len(result) > max_total_bytes:
        raise SystemExit(125)
    log.write(result)
    raise SystemExit(wrapper_exit_code)

def append_bounded(target, detail):
    bounded = str(detail)[:256]
    if len(target) < 16:
        target.append(bounded)
    elif len(target) == 16:
        target.append("additional-details-omitted")

def request_shutdown(signum, _frame):
    global requested_signal
    if requested_signal is None:
        requested_signal = signum

for handled_signal in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
    signal.signal(handled_signal, request_shutdown)

# Readiness paths arrive as supervisor-only argv and the legacy environment
# names are removed defensively before Popen. The finalizer marker participates
# in a ready/release/acknowledgement handshake: an armed supervisor cannot
# finish before the wrapper has exercised (or deliberately skipped) the
# consumed-signal forwarding latch, and the wrapper proves that the supervisor
# consumed its exact release payload.
os.environ.pop("FSIM_EULER_DISC_E2E_SELF_TEST_ACTIVE_SUPERVISOR_MARKER", None)
os.environ.pop("FSIM_EULER_DISC_E2E_SELF_TEST_FINALIZER_READY_MARKER", None)
# The arm is wrapper-local test control, not part of the supervised command's
# environment. In particular, a verifier recursively executing this harness
# must not inherit the consumed-signal mutant-killing hook.
os.environ.pop(
    "FSIM_EULER_DISC_E2E_SELF_TEST_ARM_CONSUMED_SIGNAL_FINALIZER", None
)
self_test_delay_completion_classification = (
    os.environ.pop(
        "FSIM_EULER_DISC_E2E_SELF_TEST_DELAY_COMPLETION_CLASSIFICATION", None
    )
    == "1"
)
self_test_signal_completion_classification = (
    os.environ.pop(
        "FSIM_EULER_DISC_E2E_SELF_TEST_SIGNAL_COMPLETION_CLASSIFICATION", None
    )
    == "1"
)

with open(log_path, "ab", buffering=0) as log:
    if prelaunch_refusal == "retained-log-budget-before-launch":
        append_prelaunch_result(
            log,
            wrapper_exit_code=120,
            reason="retained-log-budget-before-launch",
            launch_detail="aggregate retained-log budget exhausted before launch",
        )
    if prelaunch_refusal != "none":
        append_prelaunch_result(
            log,
            wrapper_exit_code=125,
            reason="supervisor-integrity-failure",
            launch_detail=f"unknown prelaunch refusal mode: {prelaunch_refusal}",
        )
    if supervisor_started_ns >= run_deadline_ns:
        append_prelaunch_result(
            log,
            wrapper_exit_code=121,
            reason="aggregate-deadline-before-launch",
            launch_detail="aggregate monotonic run budget exhausted before launch",
        )
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            bufsize=0,
        )
    except FileNotFoundError as error:
        append_prelaunch_result(
            log,
            wrapper_exit_code=127,
            reason="launch-not-found",
            launch_detail=f"launch failure: {error}",
        )
    except OSError as error:
        append_prelaunch_result(
            log,
            wrapper_exit_code=126,
            reason="launch-error",
            launch_detail=f"launch failure: {error}",
        )

    emergency_state = {"contained": False, "leader_reaped": False}

    def inspect_session_descendants():
        """Return live non-zombie members of the owned session, by PID/PGID.

        A descendant may call setpgid(2), so ownership is the session created by
        start_new_session=True, not merely the leader's original process group.
        """
        environment = os.environ.copy()
        environment["LC_ALL"] = "C"
        result = subprocess.run(
            ["/bin/ps", "-axo", "pid=,pgid=,stat="],
            check=False,
            capture_output=True,
            text=True,
            timeout=2,
            env=environment,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"process inspection exited {result.returncode}: {result.stderr.strip()}"
            )
        live = {}
        for line in result.stdout.splitlines():
            fields = line.split(None, 2)
            if len(fields) != 3:
                continue
            pid_text, pgid_text, state = fields
            try:
                pid = int(pid_text)
                pgid = int(pgid_text)
                sid = os.getsid(pid)
            except (ValueError, ProcessLookupError):
                continue
            if (
                sid == process.pid
                and pid != process.pid
                and not state.startswith("Z")
            ):
                live[pid] = pgid
        observed_descendants.update(live)
        return live

    def signal_owned_session(signum, action_log=None):
        """Signal the leader group and every currently observed session group/PID."""
        sent = False
        try:
            members = inspect_session_descendants()
        except (OSError, subprocess.SubprocessError, RuntimeError) as error:
            append_bounded(
                inspection_errors,
                f"{type(error).__name__}: {error}",
            )
            members = {}
        groups = set(members.values())
        for pgid in sorted(groups):
            try:
                os.killpg(pgid, signum)
                sent = True
                if action_log is not None:
                    append_bounded(action_log, f"signal-{signum}-pgid-{pgid}")
            except ProcessLookupError:
                if action_log is not None:
                    append_bounded(action_log, f"signal-{signum}-pgid-{pgid}-absent")
            except PermissionError as error:
                append_bounded(
                    inspection_errors,
                    f"group signal denied for {pgid}: {type(error).__name__}: {error}",
                )
        # A member can change groups between inspection and killpg. The direct
        # PID pass closes that migration window; subsequent scans repeat until
        # the session is empty or the drain deadline is exhausted.
        for pid in sorted(members):
            try:
                os.kill(pid, signum)
                sent = True
                if action_log is not None:
                    append_bounded(action_log, f"signal-{signum}-pid-{pid}")
            except ProcessLookupError:
                pass
            except PermissionError as error:
                append_bounded(
                    inspection_errors,
                    f"pid signal denied for {pid}: {type(error).__name__}: {error}",
                )
        try:
            os.kill(process.pid, signum)
            sent = True
            if action_log is not None:
                append_bounded(action_log, f"signal-{signum}-leader-{process.pid}")
        except ProcessLookupError:
            if action_log is not None:
                append_bounded(action_log, f"signal-{signum}-leader-{process.pid}-absent")
        except PermissionError as error:
            append_bounded(
                inspection_errors,
                f"leader signal denied: {type(error).__name__}: {error}",
            )
        return sent

    def emergency_cleanup():
        if emergency_state["contained"]:
            return
        actions = []
        drained = False
        try:
            signal_owned_session(signal.SIGTERM, actions)
            time.sleep(0.2)
            signal_owned_session(signal.SIGKILL, actions)
            if not emergency_state["leader_reaped"]:
                try:
                    process.wait(timeout=2)
                    emergency_state["leader_reaped"] = True
                    append_bounded(actions, "leader-reaped")
                except (subprocess.TimeoutExpired, ChildProcessError) as error:
                    append_bounded(actions, f"leader-reap-unestablished:{type(error).__name__}")
            try:
                if process.stdout is not None:
                    process.stdout.close()
            except OSError as error:
                append_bounded(actions, f"stdout-close-error:{error}")
            empty_deadline = time.monotonic() + 2.0
            while True:
                try:
                    live = inspect_session_descendants()
                except (OSError, subprocess.SubprocessError, RuntimeError) as error:
                    append_bounded(actions, f"session-probe-error:{type(error).__name__}:{error}")
                    live = None
                if live == {} and emergency_state["leader_reaped"]:
                    drained = True
                    break
                if time.monotonic() >= empty_deadline:
                    break
                time.sleep(0.05)
        except BaseException as error:
            append_bounded(actions, f"cleanup-error:{type(error).__name__}:{error}")
        try:
            with open(log_path, "ab", buffering=0) as emergency_log:
                payload = (
                        "\n--- emergency-supervisor-cleanup "
                        f"pgid={process.pid} sid={process.pid} "
                        f"leader_reaped={str(emergency_state['leader_reaped']).lower()} "
                        f"process_group_drained={str(drained).lower()} "
                        f"process_session_drained={str(drained).lower()} "
                        f"observed_descendants={len(observed_descendants)} "
                        f"actions={actions!r} "
                        "escaped_session_descendants_claimed=false ---\n"
                    ).encode("utf-8", "backslashreplace")
                remaining = max(0, max_total_bytes - emergency_log.tell())
                emergency_log.write(payload[:remaining])
        except BaseException:
            pass

    atexit.register(emergency_cleanup)
    def publish_regular_marker(marker, payload):
        """Expose a completely written marker without replacing any name."""
        if not hasattr(os, "O_NOFOLLOW"):
            raise RuntimeError("platform lacks no-follow marker publication")
        marker.parent.mkdir(parents=True, exist_ok=True)
        prepared = marker.with_name(f".{marker.name}.prepared-{os.getpid()}")
        flags = (
            os.O_WRONLY
            | os.O_CREAT
            | os.O_EXCL
            | os.O_NOFOLLOW
            | getattr(os, "O_CLOEXEC", 0)
        )
        descriptor = os.open(prepared, flags, 0o600)
        try:
            remaining = memoryview(payload)
            while remaining:
                written = os.write(descriptor, remaining)
                if written <= 0:
                    raise RuntimeError("prepared marker write made no progress")
                remaining = remaining[written:]
            os.fsync(descriptor)
            prepared_status = os.fstat(descriptor)
            if (
                not stat.S_ISREG(prepared_status.st_mode)
                or prepared_status.st_size != len(payload)
                or prepared_status.st_nlink != 1
            ):
                raise RuntimeError("prepared marker has invalid metadata")
        finally:
            os.close(descriptor)
        os.link(prepared, marker, follow_symlinks=False)
        published_status = marker.lstat()
        if (
            not stat.S_ISREG(published_status.st_mode)
            or (published_status.st_dev, published_status.st_ino)
            != (prepared_status.st_dev, prepared_status.st_ino)
            or published_status.st_size != len(payload)
            or published_status.st_nlink != 2
        ):
            raise RuntimeError("published marker identity is invalid")

    def publish_readiness_marker(readiness_marker_text):
        if not readiness_marker_text:
            return
        readiness_marker = pathlib.Path(readiness_marker_text)
        publish_regular_marker(
            readiness_marker,
            f"supervisor_pid={os.getpid()} child_pid={process.pid}\n".encode("ascii"),
        )

    release_marker = None
    acknowledgement_marker = None
    if finalizer_readiness_marker_text:
        if finalizer_signal_number <= 0:
            raise RuntimeError("finalizer handshake lacks a consumed signal number")
        release_marker = pathlib.Path(f"{finalizer_readiness_marker_text}.release")
        acknowledgement_marker = pathlib.Path(
            f"{finalizer_readiness_marker_text}.release-ack"
        )
        for unexpected in (release_marker, acknowledgement_marker):
            try:
                unexpected.lstat()
            except FileNotFoundError:
                continue
            raise RuntimeError(f"finalizer handshake marker already exists: {unexpected}")

    publish_readiness_marker(active_readiness_marker_text)
    publish_readiness_marker(finalizer_readiness_marker_text)
    if release_marker is not None and acknowledgement_marker is not None:
        handshake_deadline = min(
            time.monotonic() + finalizer_handshake_timeout_seconds,
            effective_deadline_ns / 1_000_000_000,
        )
        while True:
            try:
                flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
                if not hasattr(os, "O_NOFOLLOW"):
                    raise RuntimeError("platform lacks no-follow marker reads")
                descriptor = os.open(release_marker, flags | os.O_NOFOLLOW)
            except FileNotFoundError:
                descriptor = None
            if descriptor is not None:
                try:
                    release_status = os.fstat(descriptor)
                    if not stat.S_ISREG(release_status.st_mode):
                        raise RuntimeError("finalizer release marker is not a regular file")
                    if release_status.st_size > 128:
                        raise RuntimeError("finalizer release marker is oversized")
                    release_bytes = os.read(descriptor, 129)
                finally:
                    os.close(descriptor)
                expected_release = (
                    f"wrapper_release_signal={finalizer_signal_number}\n".encode("ascii")
                )
                if release_bytes != expected_release:
                    raise RuntimeError("finalizer release marker has invalid contents")
                break
            if time.monotonic() >= handshake_deadline:
                raise RuntimeError("finalizer release handshake timed out")
            time.sleep(0.01)
        acknowledgement = (
            f"supervisor_pid={os.getpid()} "
            f"wrapper_release_signal={finalizer_signal_number}\n"
        ).encode("ascii")
        publish_regular_marker(acknowledgement_marker, acknowledgement)
    assert process.stdout is not None
    descriptor = process.stdout.fileno()
    os.set_blocking(descriptor, False)
    selector = selectors.DefaultSelector()
    selector.register(descriptor, selectors.EVENT_READ)

    def leader_status_without_reaping():
        try:
            info = os.waitid(
                os.P_PID,
                process.pid,
                os.WEXITED | os.WNOHANG | os.WNOWAIT,
            )
        except ChildProcessError:
            return process.returncode
        if info is None:
            return None
        if info.si_code == os.CLD_EXITED:
            return int(info.si_status)
        return -int(info.si_status)

    def snapshot_descendants():
        try:
            return set(inspect_session_descendants())
        except (OSError, subprocess.SubprocessError, RuntimeError) as error:
            append_bounded(inspection_errors, f"{type(error).__name__}: {error}")
            return None

    if os.environ.get("FSIM_EULER_DISC_E2E_SELF_TEST_SUPERVISOR_EXCEPTION") == "1":
        observation_deadline = time.monotonic() + 1.0
        while time.monotonic() < observation_deadline:
            observed = snapshot_descendants()
            if observed:
                break
            time.sleep(0.02)
        raise RuntimeError("injected supervisor exception for no-Cargo self-test")

    def retain(data):
        global written, truncated
        available = max_output_bytes - written
        if available > 0:
            chunk = data[:available]
            log.write(chunk)
            written += len(chunk)
        if len(data) > max(available, 0):
            truncated = True

    def read_ready_output(wait_seconds):
        global eof
        for key, _ in selector.select(wait_seconds):
            try:
                data = os.read(key.fd, 65536)
            except BlockingIOError:
                continue
            if data:
                retain(data)
            else:
                eof = True
                selector.unregister(key.fd)

    leader_exit_code = None
    live_descendants = None
    while True:
        now_ns = time.monotonic_ns()
        if requested_signal is not None and shutdown_reason != "interrupt":
            shutdown_reason = "interrupt"
            interrupted_signal = requested_signal
            term_sent_at_ns = now_ns
            drain_deadline_ns = now_ns + 10_000_000_000
            snapshot_descendants()
            term_signal_sent = (
                signal_owned_session(interrupted_signal) or term_signal_sent
            )
        elif shutdown_reason is None and now_ns >= effective_deadline_ns:
            shutdown_reason = "timeout"
            term_sent_at_ns = now_ns
            drain_deadline_ns = now_ns + 10_000_000_000
            snapshot_descendants()
            term_signal_sent = signal_owned_session(signal.SIGTERM) or term_signal_sent

        leader_exit_code = leader_status_without_reaping()
        if leader_exit_code is not None:
            live_descendants = snapshot_descendants()
            if live_descendants is None and shutdown_reason is None:
                shutdown_reason = "inspection-failure"
                term_sent_at_ns = now_ns
                drain_deadline_ns = now_ns + 10_000_000_000
                term_signal_sent = signal_owned_session(signal.SIGTERM) or term_signal_sent
            elif live_descendants and shutdown_reason is None:
                shutdown_reason = "leader-exit-with-live-group"
                term_sent_at_ns = now_ns
                drain_deadline_ns = now_ns + 10_000_000_000
                term_signal_sent = signal_owned_session(signal.SIGTERM) or term_signal_sent

        # Inspection above can consume a material part of the remaining lane
        # budget. Re-sample before blocking, cap the poll to that remainder,
        # and fail closed if inspection itself crossed the deadline.
        now_ns = time.monotonic_ns()
        if shutdown_reason is None and now_ns >= effective_deadline_ns:
            shutdown_reason = "timeout"
            term_sent_at_ns = now_ns
            drain_deadline_ns = now_ns + 10_000_000_000
            snapshot_descendants()
            term_signal_sent = signal_owned_session(signal.SIGTERM) or term_signal_sent
        if (
            shutdown_reason is not None
            and not kill_signal_sent
            and now_ns >= term_sent_at_ns + 2_000_000_000
        ):
            kill_signal_sent = signal_owned_session(signal.SIGKILL) or kill_signal_sent

        wait_seconds = 0.05
        if shutdown_reason is None:
            wait_seconds = min(
                wait_seconds,
                max(0.0, (effective_deadline_ns - now_ns) / 1_000_000_000),
            )
        read_ready_output(wait_seconds)
        if leader_exit_code is not None:
            if live_descendants is None:
                live_descendants = snapshot_descendants()
        candidate_complete = (
            leader_exit_code is not None
            and live_descendants == set()
            and eof
        )
        if (
            candidate_complete
            and self_test_delay_completion_classification
            and not completion_delay_injected
        ):
            # Deterministic no-Cargo mutant killer: a stale loop-top timestamp
            # must never admit completion after the effective deadline.
            completion_delay_injected = True
            time.sleep(1.1)
        if (
            candidate_complete
            and self_test_signal_completion_classification
            and not completion_signal_injected
        ):
            # Deterministic no-Cargo mutant killer: a signal latched after the
            # loop-top check must outrank ordinary completion.
            completion_signal_injected = True
            os.kill(os.getpid(), signal.SIGTERM)

        now_ns = time.monotonic_ns()
        if requested_signal is not None and shutdown_reason != "interrupt":
            shutdown_reason = "interrupt"
            interrupted_signal = requested_signal
            term_sent_at_ns = now_ns
            drain_deadline_ns = now_ns + 10_000_000_000
            snapshot_descendants()
            term_signal_sent = signal_owned_session(interrupted_signal) or term_signal_sent
        elif shutdown_reason is None and now_ns >= effective_deadline_ns:
            shutdown_reason = "timeout"
            term_sent_at_ns = now_ns
            drain_deadline_ns = now_ns + 10_000_000_000
            snapshot_descendants()
            term_signal_sent = signal_owned_session(signal.SIGTERM) or term_signal_sent
        if candidate_complete:
            complete = not inspection_errors
            break

        if shutdown_reason is not None and now_ns >= drain_deadline_ns:
            break

    if not complete:
        kill_signal_sent = signal_owned_session(signal.SIGKILL) or kill_signal_sent
        final_deadline_ns = time.monotonic_ns() + 2_000_000_000
        while time.monotonic_ns() < final_deadline_ns:
            read_ready_output(0.05)
            leader_exit_code = leader_status_without_reaping()
            if leader_exit_code is not None:
                live_descendants = snapshot_descendants()
                if live_descendants == set() and eof and not inspection_errors:
                    complete = True
                    break

    if leader_exit_code is not None:
        try:
            process.wait(timeout=0.2)
            emergency_state["leader_reaped"] = True
        except subprocess.TimeoutExpired:
            append_bounded(inspection_errors, "exited leader could not be reaped")
            complete = False
    else:
        append_bounded(
            inspection_errors,
            "session leader did not exit within the drain bound",
        )
        complete = False
    process.stdout.close()
    selector.close()
    emergency_state["contained"] = complete

    metadata_truncated = False

    def write_metadata(payload):
        global metadata_truncated
        remaining = max(0, max_total_bytes - log.tell())
        if len(payload) > remaining:
            metadata_truncated = True
        log.write(payload[:remaining])

    if python_cap_reduced:
        write_metadata(
            f"\n--- shell-effective child-output cap {shell_max_output_bytes} reduced to "
            f"{max_output_bytes} to preserve the total retained-log bound ---\n".encode()
        )
    if truncated:
        write_metadata(
            f"\n--- output cap reached at {max_output_bytes} bytes; later observed bytes were consumed but not retained ---\n".encode()
        )
    if shutdown_reason == "timeout":
        write_metadata(
            f"\n--- monotonic timeout at {deadline_kind} deadline "
            f"{effective_deadline_ns} ---\n".encode()
        )
    elif shutdown_reason == "interrupt":
        write_metadata(
            f"\n--- supervisor received signal {interrupted_signal}; bounded shutdown requested ---\n".encode()
        )
    elif shutdown_reason == "leader-exit-with-live-group":
        write_metadata(b"\n--- leader exited while live same-session descendants remained ---\n")
    elif shutdown_reason == "inspection-failure":
        write_metadata(b"\n--- process-group inspection failed closed ---\n")
    write_metadata(
        (
            f"\n--- supervisor-state supervisor_sid={process.pid} "
            f"process_group_drained={str(complete).lower()} "
            f"process_session_drained={str(complete).lower()} "
            f"output_pipe_eof={str(eof).lower()} "
            f"observed_descendants={len(observed_descendants)} "
            f"completion_delay_hook_fired={str(completion_delay_injected).lower()} "
            f"completion_signal_hook_fired={str(completion_signal_injected).lower()} "
            f"term_sent={str(term_signal_sent).lower()} "
            f"kill_sent={str(kill_signal_sent).lower()} "
            f"retained_log_hard_limit={max_total_bytes} "
            f"inspection_errors={inspection_errors!r} ---\n"
        ).encode("utf-8", "backslashreplace")
    )

    if not complete or metadata_truncated or python_cap_reduced:
        wrapper_exit_code = 125
        result_reason = "supervisor-integrity-failure"
    elif interrupted_signal is not None:
        wrapper_exit_code = 128 + interrupted_signal
        result_reason = "interrupt"
    elif shutdown_reason == "timeout":
        wrapper_exit_code = 124
        result_reason = "timeout"
    elif shutdown_reason == "inspection-failure":
        wrapper_exit_code = 125
        result_reason = "inspection-failure"
    elif shutdown_reason == "leader-exit-with-live-group":
        wrapper_exit_code = 123
        result_reason = "leader-exit-with-live-group"
    elif truncated:
        wrapper_exit_code = 122
        result_reason = "output-truncated"
    elif leader_exit_code is None:
        wrapper_exit_code = 125
        result_reason = "supervisor-integrity-failure"
    elif leader_exit_code < 0:
        wrapper_exit_code = 128 + (-leader_exit_code)
        result_reason = "leader-signal"
    else:
        wrapper_exit_code = leader_exit_code
        result_reason = "none" if leader_exit_code == 0 else "leader-exit"

    result_payload = encode_result(
        wrapper_exit_code=wrapper_exit_code,
        reason=result_reason,
        leader_exit_code=leader_exit_code,
        process_group_drained=complete,
        process_session_drained=complete,
        output_pipe_eof=eof,
        metadata_complete=not metadata_truncated,
    )
    if log.tell() + len(result_payload) > max_total_bytes:
        raise SystemExit(125)
    log.write(result_payload)

raise SystemExit(wrapper_exit_code)
PY
  supervisor_pid=$!
  ACTIVE_SUPERVISOR_PID="$supervisor_pid"
  if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_ARM_CONSUMED_SIGNAL_FINALIZER:-0}" == 1 \
      && "$CONSUMED_WRAPPER_SIGNAL" != 0 ]]; then
    # Exactly one post-consumption finalizer must prove that its handlers are
    # armed and remain live until this wrapper has exercised the old forwarding
    # latch. Consuming the arm here prevents later verifier/publication helpers
    # from demanding a second handshake that was intentionally one-shot.
    unset FSIM_EULER_DISC_E2E_SELF_TEST_ARM_CONSUMED_SIGNAL_FINALIZER
    finalizer_release_marker="${finalizer_marker}.release"
    finalizer_ack_marker="${finalizer_marker}.release-ack"
    if [[ -z "$finalizer_marker" ]] || ! python3 - \
        "$finalizer_marker" "$finalizer_release_marker" "$finalizer_ack_marker" \
        "$supervisor_pid" "$CONSUMED_WRAPPER_SIGNAL" \
        "$FINALIZER_HANDSHAKE_TIMEOUT_SECONDS" <<'PY'
import os
import pathlib
import re
import stat
import sys
import time

marker = pathlib.Path(sys.argv[1])
release = pathlib.Path(sys.argv[2])
acknowledgement = pathlib.Path(sys.argv[3])
supervisor_pid = int(sys.argv[4])
signal_number = int(sys.argv[5])
deadline = time.monotonic() + int(sys.argv[6])

def read_regular_nofollow(path, maximum):
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    if not hasattr(os, "O_NOFOLLOW"):
        raise SystemExit("platform lacks no-follow marker reads")
    try:
        descriptor = os.open(path, flags | os.O_NOFOLLOW)
    except FileNotFoundError:
        return None
    try:
        status = os.fstat(descriptor)
        if not stat.S_ISREG(status.st_mode) or status.st_size > maximum:
            raise SystemExit(f"invalid finalizer marker metadata: {path}")
        return os.read(descriptor, maximum + 1)
    finally:
        os.close(descriptor)

for unexpected in (release, acknowledgement):
    try:
        unexpected.lstat()
    except FileNotFoundError:
        continue
    raise SystemExit(f"finalizer handshake marker existed before readiness: {unexpected}")

while True:
    ready = read_regular_nofollow(marker, 128)
    if ready is not None:
        match = re.fullmatch(rb"supervisor_pid=([0-9]+) child_pid=([0-9]+)\n", ready)
        if match is None or int(match.group(1)) != supervisor_pid:
            raise SystemExit("finalizer readiness marker is not bound to its supervisor")
        try:
            os.kill(supervisor_pid, 0)
        except ProcessLookupError as error:
            raise SystemExit("finalizer supervisor exited after readiness") from error
        break
    try:
        os.kill(supervisor_pid, 0)
    except ProcessLookupError as error:
        raise SystemExit("finalizer supervisor exited before handler readiness") from error
    if time.monotonic() >= deadline:
        raise SystemExit("finalizer supervisor handler readiness timed out")
    time.sleep(0.01)
print(f"consumed-wrapper-signal-finalizer-ready={signal_number}")
PY
    then
      finalizer_readiness_failed=1
      printf '%s\n' 'consumed wrapper signal finalizer readiness failed' >&2
    else
      if [[ "$WRAPPER_SIGNAL" != 0 ]]; then
        if kill -"$WRAPPER_SIGNAL" "$supervisor_pid" 2>/dev/null; then
          forwarding_exercised=1
        else
          finalizer_readiness_failed=1
          printf '%s\n' \
            'consumed wrapper signal forwarding exercise failed' >&2
        fi
      fi
      if ! python3 - "$finalizer_release_marker" "$finalizer_ack_marker" \
          "$supervisor_pid" "$CONSUMED_WRAPPER_SIGNAL" \
          "$FINALIZER_ACK_TIMEOUT_SECONDS" <<'PY'
import os
import pathlib
import stat
import sys
import time

release = pathlib.Path(sys.argv[1])
acknowledgement = pathlib.Path(sys.argv[2])
supervisor_pid = int(sys.argv[3])
signal_number = int(sys.argv[4])
deadline = time.monotonic() + int(sys.argv[5])

def publish_regular_marker(marker, payload):
    if not hasattr(os, "O_NOFOLLOW"):
        raise SystemExit("platform lacks no-follow marker publication")
    marker.parent.mkdir(parents=True, exist_ok=True)
    prepared = marker.with_name(f".{marker.name}.prepared-{os.getpid()}")
    flags = (
        os.O_WRONLY
        | os.O_CREAT
        | os.O_EXCL
        | os.O_NOFOLLOW
        | getattr(os, "O_CLOEXEC", 0)
    )
    descriptor = os.open(prepared, flags, 0o600)
    try:
        remaining = memoryview(payload)
        while remaining:
            written = os.write(descriptor, remaining)
            if written <= 0:
                raise SystemExit("prepared marker write made no progress")
            remaining = remaining[written:]
        os.fsync(descriptor)
        prepared_status = os.fstat(descriptor)
        if (
            not stat.S_ISREG(prepared_status.st_mode)
            or prepared_status.st_size != len(payload)
            or prepared_status.st_nlink != 1
        ):
            raise SystemExit("prepared marker has invalid metadata")
    finally:
        os.close(descriptor)
    os.link(prepared, marker, follow_symlinks=False)
    published_status = marker.lstat()
    if (
        not stat.S_ISREG(published_status.st_mode)
        or (published_status.st_dev, published_status.st_ino)
        != (prepared_status.st_dev, prepared_status.st_ino)
        or published_status.st_size != len(payload)
        or published_status.st_nlink != 2
    ):
        raise SystemExit("published marker identity is invalid")

publish_regular_marker(
    release,
    f"wrapper_release_signal={signal_number}\n".encode("ascii"),
)
print(f"consumed-wrapper-signal-finalizer-released={signal_number}")

read_flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
if not hasattr(os, "O_NOFOLLOW"):
    raise SystemExit("platform lacks no-follow marker reads")
while True:
    try:
        descriptor = os.open(acknowledgement, read_flags | os.O_NOFOLLOW)
    except FileNotFoundError:
        descriptor = None
    if descriptor is not None:
        try:
            status = os.fstat(descriptor)
            if not stat.S_ISREG(status.st_mode) or status.st_size > 128:
                raise SystemExit("invalid finalizer acknowledgement metadata")
            payload = os.read(descriptor, 129)
        finally:
            os.close(descriptor)
        expected = (
            f"supervisor_pid={supervisor_pid} "
            f"wrapper_release_signal={signal_number}\n"
        ).encode("ascii")
        if payload != expected:
            raise SystemExit("finalizer acknowledgement has invalid contents")
        break
    try:
        os.kill(supervisor_pid, 0)
    except ProcessLookupError as error:
        raise SystemExit("finalizer supervisor exited before acknowledgement") from error
    if time.monotonic() >= deadline:
        raise SystemExit("finalizer acknowledgement timed out")
    time.sleep(0.01)
print(f"consumed-wrapper-signal-finalizer-release-acknowledged={signal_number}")
PY
      then
        finalizer_readiness_failed=1
        printf '%s\n' \
          'consumed wrapper signal finalizer release failed' >&2
      fi
    fi
  fi
  if [[ "$WRAPPER_SIGNAL" != 0 && "$forwarding_exercised" == 0 ]]; then
    kill -"$WRAPPER_SIGNAL" "$supervisor_pid" 2>/dev/null || true
  fi
  while :; do
    if wait "$supervisor_pid"; then
      status=0
      break
    else
      status=$?
    fi
    if [[ "$WRAPPER_SIGNAL" != 0 ]] \
        && [[ $'\n'"$(jobs -pr)"$'\n' == *$'\n'"$supervisor_pid"$'\n'* ]]; then
      continue
    fi
    break
  done
  ACTIVE_SUPERVISOR_PID=""
  if [[ "$finalizer_readiness_failed" == 1 ]]; then
    return 125
  fi
  return "$status"
}

supervisor_disposition_for_log() { # log-path observed-wrapper-exit-code
  python3 - "$1" "$2" "$PROOF_LOG_METADATA_RESERVE_BYTES" <<'PY'
import json
import os
import pathlib
import sys

log_path = pathlib.Path(sys.argv[1])
try:
    observed_rc = int(sys.argv[2])
except ValueError as error:
    raise SystemExit(f"invalid observed wrapper exit code: {error}") from error
try:
    metadata_reserve_bytes = int(sys.argv[3])
except ValueError as error:
    raise SystemExit(f"invalid supervisor metadata reserve: {error}") from error

# The canonical supervisor result is the final line and is deliberately small.
# Read a fixed-size tail instead of loading a potentially 64 MiB lane log merely
# to classify its disposition.
tail_limit = 64 * 1024
with log_path.open("rb") as handle:
    size = os.fstat(handle.fileno()).st_size
    handle.seek(max(0, size - tail_limit))
    tail = handle.read(tail_limit + 1)
if len(tail) > tail_limit:
    raise SystemExit("supervisor result is not contained in the bounded log tail")

prefix = b"supervisor_result_json="
marker = b"\n" + prefix
offset = tail.rfind(marker)
if offset >= 0:
    encoded = tail[offset + len(marker):]
elif tail.startswith(prefix):
    encoded = tail[len(prefix):]
else:
    raise SystemExit("canonical supervisor result is missing from the lane log")
if not encoded.endswith(b"\n") or b"\n" in encoded[:-1]:
    raise SystemExit("canonical supervisor result is not the final bounded log line")
try:
    result = json.loads(encoded[:-1])
except (UnicodeDecodeError, json.JSONDecodeError) as error:
    raise SystemExit(f"canonical supervisor result is malformed: {error}") from error
if not isinstance(result, dict):
    raise SystemExit("canonical supervisor result is not an object")
if result.get("schema") != "frankensim.euler-disc-contract-e2e.supervisor-result.v1":
    raise SystemExit("canonical supervisor result has the wrong schema")
wrapper_rc = result.get("wrapper_exit_code")
if type(wrapper_rc) is not int or wrapper_rc != observed_rc:
    raise SystemExit("canonical supervisor result disagrees with the observed wrapper exit")
reason = result.get("shutdown_reason")
if not isinstance(reason, str):
    raise SystemExit("canonical supervisor result has no string shutdown reason")
initial_log_bytes = result.get("initial_log_bytes")
total_log_cap_bytes = result.get("total_log_cap_bytes")
shell_output_cap_bytes = result.get("shell_effective_output_cap_bytes")
retained_output_cap_bytes = result.get("retained_output_cap_bytes")
if not all(
    type(value) is int
    for value in (
        initial_log_bytes,
        total_log_cap_bytes,
        shell_output_cap_bytes,
        retained_output_cap_bytes,
    )
):
    raise SystemExit("canonical supervisor result has non-integer log caps")
if size > total_log_cap_bytes:
    raise SystemExit("retained log exceeds its supervisor total cap")
expected_retained_output_cap = max(
    0,
    min(
        shell_output_cap_bytes,
        total_log_cap_bytes - initial_log_bytes - metadata_reserve_bytes,
    ),
)
if retained_output_cap_bytes != expected_retained_output_cap:
    raise SystemExit("supervisor retained-output cap is not the exact bounded remainder")
if reason not in {"supervisor-integrity-failure", "inspection-failure"}:
    if not (
        result.get("process_group_drained") is True
        and result.get("process_session_drained") is True
        and result.get("output_pipe_eof") is True
        and result.get("metadata_complete") is True
        and result.get("python_cap_reduced") is False
        and result.get("inspection_error_count") == 0
    ):
        raise SystemExit("admitted supervisor disposition lacks complete containment metadata")
if reason in {
    "none",
    "output-truncated",
    "leader-exit-with-live-group",
    "timeout",
    "interrupt",
    "leader-exit",
    "leader-signal",
} and type(result.get("leader_exit_code")) is not int:
    raise SystemExit("launched supervisor disposition lacks a leader exit code")
if reason in {
    "none",
    "retained-log-budget-before-launch",
    "aggregate-deadline-before-launch",
    "launch-error",
    "launch-not-found",
    "leader-exit",
    "leader-signal",
} and result.get("output_truncated") is not False:
    raise SystemExit("supervisor disposition spuriously claims output truncation")
if reason != "interrupt" and result.get("interrupted_signal") is not None:
    raise SystemExit("non-interrupt supervisor disposition carries an interrupt signal")
supervisor_started_ns = result.get("supervisor_started_monotonic_ns")
run_deadline_ns = result.get("run_deadline_monotonic_ns")
if type(supervisor_started_ns) is not int or type(run_deadline_ns) is not int:
    raise SystemExit("supervisor disposition lacks monotonic start/deadline integers")
if reason == "aggregate-deadline-before-launch":
    if supervisor_started_ns < run_deadline_ns:
        raise SystemExit("aggregate prelaunch refusal precedes its run deadline")
elif reason not in {
    "retained-log-budget-before-launch",
    "supervisor-integrity-failure",
    "inspection-failure",
} and supervisor_started_ns >= run_deadline_ns:
    raise SystemExit("launched supervisor disposition starts at or after its run deadline")

fixed = {
    "retained-log-budget-before-launch": (
        120,
        "aggregate retained-log budget exhausted before command launch",
        False,
    ),
    "aggregate-deadline-before-launch": (
        121,
        "aggregate run timeout exhausted before command launch",
        False,
    ),
    "output-truncated": (
        122,
        "command output was truncated after its process group drained",
        False,
    ),
    "leader-exit-with-live-group": (
        123,
        "command leader exited before same-session descendants; the supervisor drained the group",
        False,
    ),
    "timeout": (124, "command exceeded its monotonic deadline", False),
    "supervisor-integrity-failure": (
        125,
        "process-group drain or supervisor metadata integrity could not be established",
        True,
    ),
    "inspection-failure": (
        125,
        "process-group drain or supervisor metadata integrity could not be established",
        True,
    ),
    "launch-error": (126, "command could not be launched", False),
    "launch-not-found": (127, "command executable was not found", False),
}
if reason in fixed:
    required_rc, detail, abort = fixed[reason]
    if wrapper_rc != required_rc:
        raise SystemExit("reserved supervisor reason has the wrong wrapper exit code")
elif reason == "interrupt":
    if wrapper_rc not in {129, 130, 143}:
        raise SystemExit("wrapper interrupt has an unsupported exit code")
    if result.get("interrupted_signal") != wrapper_rc - 128:
        raise SystemExit("wrapper interrupt signal disagrees with its exit code")
    detail = "wrapper signal interrupted the command after bounded process-group cleanup"
    abort = False
elif reason == "leader-exit":
    if wrapper_rc == 0 or result.get("leader_exit_code") != wrapper_rc:
        raise SystemExit("child exit disagrees with its wrapper exit code")
    detail = f"command exited {wrapper_rc}"
    abort = False
elif reason == "leader-signal":
    if not 129 <= wrapper_rc <= 255:
        raise SystemExit("child signal produced an invalid wrapper exit code")
    if result.get("leader_exit_code") != -(wrapper_rc - 128):
        raise SystemExit("child signal disagrees with its wrapper exit code")
    detail = f"command exited {wrapper_rc}"
    abort = False
else:
    raise SystemExit(f"unsupported supervisor shutdown reason: {reason!r}")

print(f"{reason}\t{detail}\t{int(abort)}")
PY
}

run_lane_with_status() { # lane authority nonzero-status command...
  local lane="$1" authority="$2" nonzero_status="$3"
  shift 3
  local owns_commit_boundary=0
  if [[ "$LANE_COMMIT_CRITICAL" == 0 ]]; then
    LANE_COMMIT_CRITICAL=1
    owns_commit_boundary=1
  fi
  local log_rel="logs/${lane}.log"
  local log_path="$LOG_DIR/$log_rel"
  {
    printf 'schema=frankensim.euler-disc-contract-e2e.command.v1\n'
    printf 'lane=%s\n' "$lane"
    printf 'head=%s\n' "$HEAD_SHA"
    printf 'host_isa=%s\n' "$HOST_ISA"
    printf 'profile=%s\n' "$PROFILE"
    printf 'executor_declaration=%s\n' "$EXECUTOR_DECLARATION"
    printf 'executor_attestation=caller-declared-unverified\n'
    printf 'lane_timeout_seconds=%s\n' "$LANE_TIMEOUT_SECONDS"
    printf 'run_timeout_seconds=%s\n' "$RUN_TIMEOUT_SECONDS"
    printf 'run_started_monotonic_ns=%s\n' "$RUN_STARTED_MONOTONIC_NS"
    printf 'run_deadline_monotonic_ns=%s\n' "$RUN_DEADLINE_MONOTONIC_NS"
    printf 'lane_log_max_bytes=%s\n' "$LANE_LOG_MAX_BYTES"
    printf 'authority=%s\n' "$authority"
    python3 - "$@" <<'PY'
import json
import sys

print("argv_json=" + json.dumps(sys.argv[1:], ensure_ascii=True, separators=(",", ":")))
PY
    printf '%s\n' '--- command output ---'
  } >"$log_path"
  local retained_budget_remaining total_budget_remaining per_log_total_cap
  local header_bytes effective_output_cap prelaunch_refusal="none"
  header_bytes="$(wc -c <"$log_path" | tr -d ' ')"
  retained_budget_remaining=$((
    PROOF_OPERATIONAL_LOG_BUDGET - RETAINED_LOG_BYTES_TOTAL
  ))
  total_budget_remaining=$((PROOF_MAX_TOTAL_LOG_BYTES - RETAINED_LOG_BYTES_TOTAL))
  if ((retained_budget_remaining <= header_bytes + PROOF_LOG_METADATA_RESERVE_BYTES)); then
    prelaunch_refusal="retained-log-budget-before-launch"
    effective_output_cap=0
    per_log_total_cap=$((header_bytes + PROOF_LOG_METADATA_RESERVE_BYTES))
    if ((per_log_total_cap > total_budget_remaining)); then
      per_log_total_cap="$total_budget_remaining"
    fi
  else
    per_log_total_cap="$retained_budget_remaining"
    if ((per_log_total_cap > PROOF_MAX_LANE_LOG_BYTES)); then
      per_log_total_cap="$PROOF_MAX_LANE_LOG_BYTES"
    fi
    effective_output_cap=$((
      per_log_total_cap - header_bytes - PROOF_LOG_METADATA_RESERVE_BYTES
    ))
    if ((effective_output_cap > LANE_LOG_MAX_BYTES)); then
      effective_output_cap="$LANE_LOG_MAX_BYTES"
    fi
  fi
  if run_bounded_command \
    "$log_path" "$LANE_TIMEOUT_SECONDS" "$RUN_DEADLINE_MONOTONIC_NS" \
    "$LANE_LOG_MAX_BYTES" "$effective_output_cap" \
    "$per_log_total_cap" "$prelaunch_refusal" "$@"; then
      LAST_LANE_RC=0
      LAST_LANE_STATUS=PASS
      record "$lane" PASS "$authority" "command completed" "$log_rel"
  else
    local rc=$?
    local detail disposition disposition_abort supervisor_reason
    local status="$nonzero_status"
    if disposition="$(supervisor_disposition_for_log "$log_path" "$rc")"; then
      IFS=$'\t' read -r supervisor_reason detail disposition_abort <<<"$disposition"
      if [[ -z "$supervisor_reason" || -z "$detail" \
          || ! "$disposition_abort" =~ ^[01]$ ]]; then
        disposition_abort=1
      fi
    else
      disposition_abort=1
    fi
    if [[ "$disposition_abort" == 1 ]]; then
      rc=125
      supervisor_reason="supervisor-integrity-failure"
      detail="process-group drain or supervisor metadata integrity could not be established"
      status=FAIL
    fi
    LAST_LANE_RC="$rc"
    if [[ "$nonzero_status" == NO_DATA && "$rc" != 1 ]]; then
      status=FAIL
    fi
    LAST_LANE_STATUS="$status"
    record "$lane" "$status" "$authority" "$detail" "$log_rel"
    if [[ "$supervisor_reason" == "supervisor-integrity-failure" \
        || "$supervisor_reason" == "inspection-failure" ]]; then
      printf '%s\n' \
        'supervisor containment is unestablished; aborting without scheduling another lane' \
        >&2
      exit 125
    fi
  fi
  if [[ "$owns_commit_boundary" == 1 ]]; then
    LANE_COMMIT_CRITICAL=0
    if [[ "$WRAPPER_SIGNAL" != 0 ]]; then
      exit_for_wrapper_signal "$WRAPPER_SIGNAL"
    fi
  fi
}

run_lane() { # lane authority command...
  local lane="$1" authority="$2"
  shift 2
  run_lane_with_status "$lane" "$authority" FAIL "$@"
}

run_refusal_lane() { # lane authority command...
  local lane="$1" authority="$2"
  shift 2
  run_lane_with_status "$lane" "$authority" NO_DATA "$@"
}

SCOPE_PATHS=(
  Cargo.toml
  Cargo.lock
  README.md
  consolidation-review.json
  docs/CONSOLIDATION_REVIEW.md
  docs/CONVENTIONS.md
  docs/SCHEMA_POLICY.md
  doc-facts-inventory.json
  schema-policy.json
  identity-authorities.json
  identity-schemas.json
  golden-couplings.json
  xtask/src/identities.rs
  scripts/ci/euler_disc_contract_e2e.sh
  crates/fs-euler-disc-e2e/Cargo.toml
  crates/fs-euler-disc-e2e/CONTRACT.md
  crates/fs-euler-disc-e2e/src/lib.rs
  crates/fs-euler-disc-e2e/src/contract.rs
  crates/fs-euler-disc-e2e/src/protocol.rs
  crates/fs-euler-disc-e2e/tests/scientific_contract.rs
)

# shellcheck disable=SC2329 # Exported into bounded `bash -c` lanes below.
closure_root_preflight_command() { # expected-head scope-path...
  local expected_head="$1" actual_head root_status special_flags fsmonitor_flags path
  local expected_blob actual_blob command_output fsmonitor_output command_rc
  shift
  if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_PREFLIGHT_INFRA_FAILURE:-0}" == 1 ]]; then
    printf '%s\n' 'self-test injected closure preflight infrastructure failure'
    return 2
  fi
  if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_PREFLIGHT_POLICY_REFUSAL:-0}" == 1 ]]; then
    printf '%s\n' 'self-test injected closure preflight policy refusal'
    return 1
  fi
  if actual_head="$(git rev-parse HEAD)"; then
    :
  else
    command_rc=$?
    printf 'infrastructure-command=git-rev-parse-head exit=%s\n' "$command_rc"
    return "$command_rc"
  fi
  printf 'expected-head=%s\nactual-head=%s\n' "$expected_head" "$actual_head"
  [[ "$actual_head" == "$expected_head" ]] || return 1
  # Read index flags before `git status`: status refresh may clear valid
  # fsmonitor bits and would otherwise make the closure check self-erasing.
  if command_output="$(git ls-files -v)"; then
    special_flags="$(printf '%s\n' "$command_output" | awk 'substr($0, 1, 1) ~ /[a-zS]/')"
  else
    command_rc=$?
    printf 'infrastructure-command=git-ls-files-special-flags exit=%s\n' "$command_rc"
    return "$command_rc"
  fi
  if fsmonitor_output="$(git ls-files -f)"; then
    fsmonitor_flags="$(printf '%s\n' "$fsmonitor_output" | awk 'substr($0, 1, 1) ~ /[a-z]/')"
  else
    command_rc=$?
    printf 'infrastructure-command=git-ls-files-fsmonitor-flags exit=%s\n' "$command_rc"
    return "$command_rc"
  fi
  if root_status="$(git status --porcelain=v1 --untracked-files=all)"; then
    :
  else
    command_rc=$?
    printf 'infrastructure-command=git-status exit=%s\n' "$command_rc"
    return "$command_rc"
  fi
  if [[ -n "$root_status" ]]; then
    printf '%s\n' 'repository-status=dirty' "$root_status"
    return 1
  fi
  printf '%s\n' 'repository-status=clean'
  if [[ -n "$special_flags" ]]; then
    printf '%s\n' 'special-index-flags=present' "$special_flags"
    return 1
  fi
  printf '%s\n' 'special-index-flags=absent'
  if [[ -n "$fsmonitor_flags" ]]; then
    printf '%s\n' 'fsmonitor-valid-flags=present' "$fsmonitor_flags"
    return 1
  fi
  printf '%s\n' 'fsmonitor-valid-flags=absent'
  for path in "$@"; do
    if command_output="$(git ls-files --error-unmatch -- "$path")"; then
      printf '%s\n' "$command_output"
    else
      command_rc=$?
      printf 'scope-index-lookup-failed=%s exit=%s\n' "$path" "$command_rc"
      return "$command_rc"
    fi
    if [[ ! -f "$path" || -L "$path" ]]; then
      printf 'scope-not-regular=%s\n' "$path"
      return 1
    fi
    if expected_blob="$(git rev-parse "${expected_head}:${path}")"; then
      :
    else
      command_rc=$?
      printf 'infrastructure-command=git-rev-parse-scope path=%s exit=%s\n' \
        "$path" "$command_rc"
      return "$command_rc"
    fi
    if actual_blob="$(git hash-object -- "$path")"; then
      :
    else
      command_rc=$?
      printf 'infrastructure-command=git-hash-object path=%s exit=%s\n' \
        "$path" "$command_rc"
      return "$command_rc"
    fi
    printf 'scope-head-blob=%s scope-worktree-blob=%s path=%s\n' \
      "$expected_blob" "$actual_blob" "$path"
    [[ "$actual_blob" == "$expected_blob" ]] || return 1
  done
}

# shellcheck disable=SC2329 # Exported into bounded `bash -c` lanes below.
bounded_snapshot_append() { # candidate maximum-bytes
  python3 /dev/fd/3 "$1" "$2" 3<<'PY'
import os
import pathlib
import stat
import sys

candidate = pathlib.Path(sys.argv[1])
maximum = int(sys.argv[2])
flags = os.O_WRONLY | os.O_APPEND | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
descriptor = os.open(candidate, flags)
initial = os.fstat(descriptor)
if not stat.S_ISREG(initial.st_mode) or initial.st_nlink != 1 or initial.st_size > maximum:
    os.close(descriptor)
    raise SystemExit("snapshot candidate violates its retained-byte contract")
remaining = maximum - initial.st_size
overflow = False
with os.fdopen(descriptor, "ab", buffering=0) as destination:
    while True:
        chunk = sys.stdin.buffer.read(65536)
        if not chunk:
            break
        retained = chunk[:remaining]
        if retained:
            destination.write(retained)
            sys.stdout.buffer.write(retained)
            remaining -= len(retained)
        if len(retained) != len(chunk):
            overflow = True
    destination.flush()
    os.fsync(destination.fileno())
    final = os.fstat(destination.fileno())
    if (
        final.st_dev != initial.st_dev
        or final.st_ino != initial.st_ino
        or final.st_nlink != 1
        or final.st_size > maximum
    ):
        raise SystemExit("snapshot candidate changed identity or exceeded its bound")
if overflow:
    print(
        f"snapshot retained-byte bound exceeded: maximum={maximum}",
        file=sys.stderr,
    )
    raise SystemExit(122)
PY
}

# shellcheck disable=SC2329 # Exported into bounded `bash -c` lanes below.
append_snapshot_scope_hash() { # candidate scope-path maximum-bytes
  local candidate="$1" path="$2" maximum_bytes="$3" digest command_rc
  if digest="$(sha256_file "$path")"; then
    :
  else
    command_rc=$?
    printf 'snapshot-scope-hash-failed=%s exit=%s\n' "$path" "$command_rc" >&2
    return "$command_rc"
  fi
  printf 'scope-sha256=%s path=%s\n' "$digest" "$path" | \
    bounded_snapshot_append "$candidate" "$maximum_bytes"
}

# shellcheck disable=SC2329 # Exported into bounded `bash -c` lanes below.
capture_snapshot_command() { # destination candidate expected-head scope-path...
  local destination="$1" candidate="$2" expected_head="$3"
  local actual_head command_rc path line
  shift 3
  set -o pipefail
  scripts/ci/checkout_constellation.sh --snapshot | \
    bounded_snapshot_append "$candidate" "$PROOF_MAX_SNAPSHOT_BYTES" || return $?
  if actual_head="$(git rev-parse HEAD)"; then
    :
  else
    command_rc=$?
    printf 'snapshot-head-read-failed=%s\n' "$command_rc" >&2
    return "$command_rc"
  fi
  line="expected-head=$expected_head actual-head=$actual_head"
  printf '%s\n' "$line" | \
    bounded_snapshot_append "$candidate" "$PROOF_MAX_SNAPSHOT_BYTES" || return $?
  [[ "$actual_head" == "$expected_head" ]] || return 1
  git status --porcelain=v1 --untracked-files=all | \
    sed 's/^/root-status=/' | \
    bounded_snapshot_append "$candidate" "$PROOF_MAX_SNAPSHOT_BYTES" || return $?
  for path in "$@"; do
    if [[ ! -f "$path" ]]; then
      printf 'scope-missing=%s\n' "$path" | \
        bounded_snapshot_append "$candidate" "$PROOF_MAX_SNAPSHOT_BYTES" || return $?
      return 1
    fi
    append_snapshot_scope_hash \
      "$candidate" "$path" "$PROOF_MAX_SNAPSHOT_BYTES" || return $?
  done
  mv "$candidate" "$destination"
}

# shellcheck disable=SC2329 # Exported into bounded `bash -c` lanes below.
source_manifest_membership_command() { # manifest scope-path...
  python3 - "$@" <<'PY'
import collections
import json
import sys

manifest_path, *expected = sys.argv[1:]
with open(manifest_path, encoding="utf-8") as handle:
    manifest = json.load(handle)
counts = collections.Counter(row["path"] for row in manifest["frankensim"]["files"])
invalid = sorted(path for path in expected if counts[path] != 1)
for path in sorted(expected):
    print(f"count={counts[path]} path={path}")
if invalid:
    raise SystemExit(1)
PY
}

# shellcheck disable=SC2329 # Exported into bounded `bash -c` lanes below.
checker_smoke_sentinel_command() { # lane-log exact-test-name
  python3 - "$@" <<'PY'
import re
import sys

log_path, test_name = sys.argv[1:]
with open(log_path, encoding="utf-8", errors="replace") as handle:
    text = handle.read()
expected = f"test {test_name} ... ok"
if text.count(expected) != 1:
    raise SystemExit(f"expected exactly one passing sentinel {expected!r}")
if not re.search(r"(?m)^running 1 test$", text):
    raise SystemExit("checker smoke did not report exactly one selected test")
if not re.search(r"test result: ok\. 1 passed; 0 failed; 0 ignored; 0 measured;", text):
    raise SystemExit("checker smoke did not report one non-vacuous passing result")
print(f"verified-one-test-sentinel={test_name}")
PY
}

# shellcheck disable=SC2329 # Exported into bounded `bash -c` lanes below.
wrapper_signal_self_test_command() { # script-path signal-name [target-mode]
  python3 - "$@" <<'PY'
import os
import pathlib
import signal
import subprocess
import sys
import tempfile
import time

script, signal_name, *mode = sys.argv[1:]
signum = getattr(signal, f"SIG{signal_name}")
environment = os.environ.copy()
publication_marker = None
readiness_marker = None
if mode == ["ordinary-lane"]:
    environment["FSIM_EULER_DISC_E2E_SELF_TEST_ORDINARY_SIGNAL_WINDOW"] = "1"
    environment[
        "FSIM_EULER_DISC_E2E_SELF_TEST_ARM_CONSUMED_SIGNAL_FINALIZER"
    ] = "1"
elif not mode or mode == ["dedicated-target"]:
    environment["FSIM_EULER_DISC_E2E_SELF_TEST_SIGNAL_TARGET"] = "1"
    environment[
        "FSIM_EULER_DISC_E2E_SELF_TEST_ARM_CONSUMED_SIGNAL_FINALIZER"
    ] = "1"
elif len(mode) == 2 and mode[0] == "publication-window":
    auxiliary = pathlib.Path(mode[1])
    auxiliary.mkdir(parents=True, exist_ok=True)
    publication_marker = auxiliary / "publication-ready.txt"
    environment["FSIM_EULER_DISC_E2E_SELF_TEST_PUBLICATION_SIGNAL_TARGET"] = "1"
    environment["FSIM_EULER_DISC_E2E_LOG_DIR"] = str(auxiliary / "bundles")
    environment["FSIM_EULER_DISC_E2E_SELF_TEST_PUBLICATION_GAP_MARKER"] = str(
        publication_marker
    )
else:
    raise SystemExit(f"unknown wrapper-signal self-test mode: {mode!r}")
if publication_marker is None:
    readiness_root = pathlib.Path(tempfile.mkdtemp(prefix="fsim-euler-supervisor-ready-"))
    readiness_marker = readiness_root / "active-supervisor.txt"
    environment["FSIM_EULER_DISC_E2E_SELF_TEST_ACTIVE_SUPERVISOR_MARKER"] = str(
        readiness_marker
    )
else:
    readiness_marker = publication_marker
process = subprocess.Popen([script, "--self-test"], env=environment)
readiness_timeout = int(os.environ["WRAPPER_SIGNAL_READINESS_TIMEOUT_SECONDS"])
termination_timeout = int(os.environ["WRAPPER_SIGNAL_TERMINATION_TIMEOUT_SECONDS"])
deadline = time.monotonic() + readiness_timeout
while not readiness_marker.exists():
    if process.poll() is not None:
        raise SystemExit(
            f"signal target exited before readiness marker: {process.returncode}"
        )
    if time.monotonic() >= deadline:
        process.kill()
        process.wait(timeout=5)
        raise SystemExit("signal target readiness marker was not created")
    time.sleep(0.02)
if publication_marker is not None:
    print(f"publication-window-marker={publication_marker}", flush=True)
else:
    print(f"active-supervisor-marker={readiness_marker}", flush=True)
os.kill(process.pid, signum)
try:
    status = process.wait(timeout=termination_timeout)
except subprocess.TimeoutExpired:
    process.kill()
    process.wait(timeout=5)
    raise SystemExit("signalled wrapper did not terminate within the self-test bound")
raise SystemExit(status)
PY
}

# shellcheck disable=SC2329 # Exported into one consolidated bounded lane below.
boundary_signal_injection_self_test_command() { # script-path retained-root
  python3 - "$@" <<'PY'
import os
import json
import pathlib
import subprocess
import sys

script, retained_root_text = sys.argv[1:]
retained_root = pathlib.Path(retained_root_text)
retained_root.mkdir(parents=True, exist_ok=True)
injection_names = (
    "FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_PREPUBLICATION_SIGNAL",
    "FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_POSTPUBLICATION_SIGNAL",
    "FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_FINAL_COMMIT_SIGNAL",
    "FSIM_EULER_DISC_E2E_SELF_TEST_REPLACE_PUBLIC_TOMBSTONE",
    "FSIM_EULER_DISC_E2E_SELF_TEST_PUBLICATION_SIGNAL_TARGET",
)
failed = False

def run_case(phase, expected_rc, log_root, updates):
    global failed
    environment = os.environ.copy()
    for injection_name in injection_names:
        environment.pop(injection_name, None)
    environment.update({
        "FSIM_EULER_DISC_E2E_LOG_DIR": str(log_root),
        "FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_INCOMPLETE": "1",
    })
    environment.update(updates)
    try:
        result = subprocess.run(
            [script, "--self-test"],
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=10,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        print(f"boundary-signal-mode={phase} timeout=true", flush=True)
        if error.stdout:
            print(error.stdout, end="", flush=True)
        failed = True
        return None
    print(f"boundary-signal-mode={phase} rc={result.returncode}", flush=True)
    print(result.stdout, end="", flush=True)
    if result.returncode != expected_rc:
        failed = True
    return result

nofollow_real = retained_root / "log-root-nofollow-real"
nofollow_link = retained_root / "log-root-nofollow-link"
nofollow_real.mkdir()
nofollow_link.symlink_to(nofollow_real.resolve(), target_is_directory=True)
nofollow_result = run_case("log-root-nofollow", 37, nofollow_link, {})
if nofollow_result is not None:
    publication_logs = sorted(
        nofollow_real.glob(".success-finalization-publication.log-*")
    )
    publication_text = "\n".join(
        path.read_text(encoding="utf-8", errors="replace")
        for path in publication_logs
    )
    public_successes = 0
    for candidate in nofollow_real.iterdir():
        if candidate.name.startswith(".") or not candidate.is_dir():
            continue
        summary = candidate / "summary.json"
        if not summary.is_file():
            continue
        status = json.loads(summary.read_text(encoding="utf-8"))["status"]
        if status in {"FOCUSED_PASS", "READY_FOR_DSR", "SELF_TEST_PASS"}:
            public_successes += 1
    nofollow_refused = (
        "require_atomic_retraction_support" in publication_text
        and ("NotADirectoryError" in publication_text or "Too many levels" in publication_text)
        and public_successes == 0
    )
    print(
        f"log-root-nofollow-refusal={str(nofollow_refused).lower()} "
        f"published-success-count={public_successes}",
        flush=True,
    )
    if not nofollow_refused:
        failed = True

phases = (
    (
        "tombstone-mismatch",
        125,
        {
            "FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_POSTPUBLICATION_SIGNAL": "1",
            "FSIM_EULER_DISC_E2E_SELF_TEST_REPLACE_PUBLIC_TOMBSTONE": "1",
            "FSIM_EULER_DISC_E2E_SELF_TEST_PUBLICATION_SIGNAL_TARGET": "1",
        },
    ),
    (
        "prepublication",
        143,
        {"FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_PREPUBLICATION_SIGNAL": "1"},
    ),
    (
        "final-commit",
        143,
        {"FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_FINAL_COMMIT_SIGNAL": "1"},
    ),
)
for phase, expected_rc, updates in phases:
    result = run_case(phase, expected_rc, retained_root / phase, updates)
    if phase == "tombstone-mismatch" and result is not None:
        mismatch_root = retained_root / phase
        public_successes = 0
        for candidate in mismatch_root.iterdir():
            if candidate.name.startswith(".") or not candidate.is_dir():
                continue
            summary = candidate / "summary.json"
            if not summary.is_file():
                continue
            status = json.loads(summary.read_text(encoding="utf-8"))["status"]
            if status in {"FOCUSED_PASS", "READY_FOR_DSR", "SELF_TEST_PASS"}:
                public_successes += 1
        print(
            f"tombstone-mismatch-published-success-count={public_successes}",
            flush=True,
        )
        if public_successes != 0:
            failed = True
raise SystemExit(1 if failed else 143)
PY
}

# shellcheck disable=SC2329 # Exported into a bounded self-test child.
publication_gap_mutation_self_test_command() { # script log-root marker
  python3 - "$@" <<'PY'
import os
import pathlib
import subprocess
import sys
import time

script, log_root, marker_text = sys.argv[1:]
marker = pathlib.Path(marker_text)
environment = os.environ.copy()
environment.update({
    "FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_INCOMPLETE": "1",
    "FSIM_EULER_DISC_E2E_LOG_DIR": log_root,
    "FSIM_EULER_DISC_E2E_SELF_TEST_PUBLICATION_GAP_MARKER": str(marker),
})
process = subprocess.Popen([script, "--self-test"], env=environment)
deadline = time.monotonic() + 15
while not marker.exists():
    if process.poll() is not None:
        raise SystemExit(
            f"publication-gap target exited before marker: {process.returncode}"
        )
    if time.monotonic() >= deadline:
        process.kill()
        process.wait(timeout=5)
        raise SystemExit("publication-gap marker was not created")
    time.sleep(0.02)
line = marker.read_text(encoding="utf-8").strip()
if not line.startswith("candidate="):
    process.kill()
    process.wait(timeout=5)
    raise SystemExit(f"publication-gap marker is malformed: {line!r}")
candidate = pathlib.Path(line.removeprefix("candidate="))
with (candidate / "summary.json").open("ab", buffering=0) as handle:
    handle.write(b"external-publication-gap-mutation\n")
print(f"external-publication-gap-mutation={candidate / 'summary.json'}", flush=True)
try:
    status = process.wait(timeout=20)
except subprocess.TimeoutExpired:
    process.kill()
    process.wait(timeout=5)
    raise SystemExit("publication-gap target did not terminate")
raise SystemExit(status)
PY
}

# shellcheck disable=SC2329 # Exported into a bounded self-test child.
verifier_toctou_mutation_self_test_command() { # script bundle marker
  python3 - "$@" <<'PY'
import os
import pathlib
import subprocess
import sys
import time

script, bundle_text, marker_text = sys.argv[1:]
bundle = pathlib.Path(bundle_text)
marker = pathlib.Path(marker_text)
environment = os.environ.copy()
environment.update({
    "FSIM_EULER_DISC_E2E_SELF_TEST_PAUSE_DURING_VERIFY": "1",
    "FSIM_EULER_DISC_E2E_SELF_TEST_VERIFY_READY_MARKER": str(marker),
})
process = subprocess.Popen([script, "--verify-bundle", str(bundle)], env=environment)
deadline = time.monotonic() + 10
while not marker.exists():
    if process.poll() is not None:
        raise SystemExit(
            f"TOCTOU verifier exited before readiness marker: {process.returncode}"
        )
    if time.monotonic() >= deadline:
        process.kill()
        process.wait(timeout=5)
        raise SystemExit("TOCTOU verifier readiness marker was not created")
    time.sleep(0.02)
with (bundle / "logs" / "crate-fmt.log").open("ab", buffering=0) as handle:
    handle.write(b"external concurrent mutation\n")
print(f"toctou-readiness-marker={marker}", flush=True)
try:
    status = process.wait(timeout=10)
except subprocess.TimeoutExpired:
    process.kill()
    process.wait(timeout=5)
    raise SystemExit("TOCTOU verifier did not terminate")
raise SystemExit(status)
PY
}

# shellcheck disable=SC2329 # Exported into a bounded self-test child.
publication_destination_collision_self_test_command() { # script log-root marker
  python3 - "$@" <<'PY'
import os
import pathlib
import subprocess
import sys
import time

script, log_root, marker_text = sys.argv[1:]
marker = pathlib.Path(marker_text)
environment = os.environ.copy()
environment.update({
    "FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_INCOMPLETE": "1",
    "FSIM_EULER_DISC_E2E_LOG_DIR": log_root,
    "FSIM_EULER_DISC_E2E_SELF_TEST_PUBLICATION_COLLISION_MARKER": str(marker),
})
process = subprocess.Popen([script, "--self-test"], env=environment)
deadline = time.monotonic() + 15
while not marker.exists():
    if process.poll() is not None:
        raise SystemExit(
            f"publication-collision target exited before marker: {process.returncode}"
        )
    if time.monotonic() >= deadline:
        process.kill()
        process.wait(timeout=5)
        raise SystemExit("publication-collision marker was not created")
    time.sleep(0.02)
values = {}
for line in marker.read_text(encoding="utf-8").splitlines():
    key, separator, value = line.partition("=")
    if not separator or key in values:
        process.kill()
        process.wait(timeout=5)
        raise SystemExit(f"publication-collision marker is malformed: {line!r}")
    values[key] = value
if set(values) != {"candidate", "published"}:
    process.kill()
    process.wait(timeout=5)
    raise SystemExit(f"publication-collision marker has wrong fields: {values!r}")
published = pathlib.Path(values["published"])
published.mkdir()
before = published.stat()
print(f"publication-collision-destination={published}", flush=True)
try:
    status = process.wait(timeout=25)
except subprocess.TimeoutExpired:
    process.kill()
    process.wait(timeout=5)
    raise SystemExit("publication-collision target did not terminate")
after = published.stat()
if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
    raise SystemExit("publication collision replaced the destination identity")
if any(published.iterdir()):
    raise SystemExit("publication collision mutated the destination directory")
print("publication-collision-destination-preserved=true", flush=True)
raise SystemExit(status)
PY
}

# shellcheck disable=SC2329 # Exported into bounded no-Cargo hostile lanes.
normal_harness_hostile_self_test_command() { # source-script fixture-root log-root mode
  python3 - "$@" <<'PY'
import json
import os
import pathlib
import shutil
import subprocess
import sys

source_script_text, fixture_text, log_root_text, mode = sys.argv[1:]
if mode not in {
    "invalid-consolidation",
    "scope-mutation",
    "deadline-expiry",
    "postpublication-deadline-expiry",
}:
    raise SystemExit(f"unknown normal-harness hostile mode: {mode!r}")
source_script = pathlib.Path(source_script_text).resolve()
source_root = source_script.parents[2]
fixture = pathlib.Path(fixture_text)
log_root = pathlib.Path(log_root_text)
fixture.mkdir(parents=True)
log_root.mkdir(parents=True)
bash_executable = shutil.which("bash")
if bash_executable is None:
    system_bash = pathlib.Path("/bin/bash")
    if system_bash.is_file() and os.access(system_bash, os.X_OK):
        bash_executable = str(system_bash)
if bash_executable is None:
    raise SystemExit("normal-harness hostile fixture requires Bash")
python_executable = sys.executable
if not python_executable or not os.access(python_executable, os.X_OK):
    python_executable = shutil.which("python3")
if python_executable is None:
    raise SystemExit("normal-harness hostile fixture requires Python 3")

scope_paths = (
    "Cargo.toml", "Cargo.lock", "README.md", "consolidation-review.json",
    "docs/CONSOLIDATION_REVIEW.md", "docs/CONVENTIONS.md", "docs/SCHEMA_POLICY.md",
    "doc-facts-inventory.json", "schema-policy.json",
    "identity-authorities.json", "identity-schemas.json", "golden-couplings.json",
    "xtask/src/identities.rs", "scripts/ci/euler_disc_contract_e2e.sh",
    "crates/fs-euler-disc-e2e/Cargo.toml", "crates/fs-euler-disc-e2e/CONTRACT.md",
    "crates/fs-euler-disc-e2e/src/lib.rs", "crates/fs-euler-disc-e2e/src/contract.rs",
    "crates/fs-euler-disc-e2e/src/protocol.rs",
    "crates/fs-euler-disc-e2e/tests/scientific_contract.rs",
)
for relative in scope_paths:
    destination = fixture / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    if relative == "scripts/ci/euler_disc_contract_e2e.sh":
        shutil.copy2(source_script, destination)
    elif relative == "consolidation-review.json":
        shutil.copy2(source_root / relative, destination)
    else:
        destination.write_text(f"hostile-fixture-path={relative}\n", encoding="utf-8")

checkout = fixture / "scripts/ci/checkout_constellation.sh"
true_executable = shutil.which("true")
if true_executable is None:
    system_true = pathlib.Path("/usr/bin/true")
    if system_true.is_file() and os.access(system_true, os.X_OK):
        true_executable = str(system_true)
if true_executable is None:
    raise SystemExit("normal-harness hostile fixture requires a true executable")
checkout.symlink_to(true_executable)

cargo_program = """import json
import os
import pathlib
import sys

repo = pathlib.Path(os.environ["FSIM_EULER_DISC_E2E_FAKE_REPO"])
mode = os.environ["FSIM_EULER_DISC_E2E_FAKE_MODE"]
gate = sys.argv[-1] if sys.argv[1:] else ""
review_path = repo / "consolidation-review.json"
if gate == "check-consolidation":
    review = json.loads(review_path.read_text(encoding="utf-8"))
    latest = review["reviews"][-1]
    root_count = latest["roots"].count("fs-euler-disc-e2e")
    if root_count != 1:
        print(f"fake-check-consolidation=invalid-euler-root-count count={root_count}")
        raise SystemExit(29)
    if mode == "scope-mutation":
        before = " [hostile-scope-before]"
        after = " [hostile-scope-after]"
        if before not in latest["method"]:
            raise SystemExit("scope-mutation fixture lacks its before marker")
        latest["method"] = latest["method"].replace(before, after, 1)
        review_path.write_text(
            json.dumps(review, sort_keys=True, separators=(",", ":")) + "\\n",
            encoding="utf-8",
        )
        print("fake-consolidation-mutation=before-to-after")
print("running 1 test")
print("test g0_check_receipts_logs_and_domain_separation_are_exact_and_bounded ... ok")
print("test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured;")
"""
for subcommand in ("fmt", "check", "test", "clippy", "run"):
    (fixture / subcommand).write_text(cargo_program, encoding="utf-8")

def git(*arguments, capture=False):
    return subprocess.run(
        ["git", "-C", str(fixture), *arguments],
        check=True,
        capture_output=capture,
        text=True,
    )

git("init", "-q", "-b", "main")
git("add", ".")
subprocess.run(
    [
        "git", "-C", str(fixture), "-c", "core.hooksPath=/dev/null",
        "-c", "commit.gpgSign=false", "-c", "user.name=FrankenSim Harness",
        "-c", "user.email=harness@invalid.example", "commit", "-q", "-m",
        "hostile normal-harness fixture",
    ],
    check=True,
)

review_path = fixture / "consolidation-review.json"
review = json.loads(review_path.read_text(encoding="utf-8"))
latest = review["reviews"][-1]
if latest["roots"].count("fs-euler-disc-e2e") != 1:
    raise SystemExit("source consolidation review lacks exactly one Euler workflow root")
if mode == "invalid-consolidation":
    latest["roots"] = [root for root in latest["roots"] if root != "fs-euler-disc-e2e"]
elif mode == "scope-mutation":
    latest["method"] += " [hostile-scope-before]"
if mode not in {"deadline-expiry", "postpublication-deadline-expiry"}:
    review_path.write_text(
        json.dumps(review, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )

status_before = git(
    "status", "--porcelain=v1", "--", "consolidation-review.json", capture=True,
).stdout.rstrip("\n")
expected_dirty = " M consolidation-review.json"
if mode in {"deadline-expiry", "postpublication-deadline-expiry"}:
    if status_before:
        raise SystemExit(f"deadline fixture unexpectedly dirty: {status_before!r}")
else:
    if status_before != expected_dirty:
        raise SystemExit(f"hostile fixture lacks exact M status: {status_before!r}")
print(f"fixture-status-before={status_before}", flush=True)

environment = os.environ.copy()
environment.pop("FSIM_EULER_DISC_E2E_SELF_TEST_PREPUBLICATION_DEADLINE_DELAY", None)
environment.pop("FSIM_EULER_DISC_E2E_SELF_TEST_POSTPUBLICATION_DEADLINE_DELAY", None)
environment.update({
    "FSIM_EULER_DISC_E2E_EXECUTOR": "local",
    "FSIM_EULER_DISC_E2E_ALLOW_LOCAL": "1",
    # Invoke the already installed interpreter directly. Its first Cargo-like
    # argument selects the same-named fixture program in this working tree,
    # avoiding platform-dependent execution of freshly written shebang files.
    "FSIM_EULER_DISC_E2E_CARGO": str(python_executable),
    "FSIM_EULER_DISC_E2E_LOG_DIR": str(log_root),
    "FSIM_EULER_DISC_E2E_FAKE_REPO": str(fixture),
    "FSIM_EULER_DISC_E2E_FAKE_MODE": mode,
    "FSIM_EULER_DISC_E2E_LANE_TIMEOUT_SECONDS": os.environ[
        "NORMAL_HARNESS_NESTED_LANE_TIMEOUT_SECONDS"
    ],
    "FSIM_EULER_DISC_E2E_RUN_TIMEOUT_SECONDS": os.environ[
        "NORMAL_HARNESS_NESTED_RUN_TIMEOUT_SECONDS"
    ],
})
if mode == "deadline-expiry":
    environment["FSIM_EULER_DISC_E2E_SELF_TEST_PREPUBLICATION_DEADLINE_DELAY"] = "1"
elif mode == "postpublication-deadline-expiry":
    environment["FSIM_EULER_DISC_E2E_SELF_TEST_POSTPUBLICATION_DEADLINE_DELAY"] = "1"
result = subprocess.run(
    [
        bash_executable,
        str(fixture / "scripts/ci/euler_disc_contract_e2e.sh"),
        "--profile",
        "focused",
    ],
    cwd=fixture,
    env=environment,
    check=False,
)

status_after = git(
    "status", "--porcelain=v1", "--", "consolidation-review.json", capture=True,
).stdout.rstrip("\n")
if mode == "scope-mutation" and status_after != expected_dirty:
    raise SystemExit(f"scope mutation did not retain exact M status: {status_after!r}")
print(f"fixture-status-after={status_after}", flush=True)

public_statuses = []
published_directories = []
for candidate in sorted(log_root.iterdir()):
    if candidate.name.startswith(".") or not candidate.is_dir():
        continue
    summary = candidate / "summary.json"
    if summary.is_file():
        public_statuses.append(json.loads(summary.read_text(encoding="utf-8"))["status"])
        published_directories.append(candidate)
success_count = sum(
    status in {"FOCUSED_PASS", "READY_FOR_DSR", "SELF_TEST_PASS"}
    for status in public_statuses
)
print(f"published-success-count={success_count}", flush=True)
print(f"published-statuses={','.join(public_statuses)}", flush=True)
if len(published_directories) == 1:
    # A deadline-refused success retains its original candidate under a hidden
    # name and publishes a fresh INCOMPLETE bundle. Inspect the unique retained
    # candidate that actually contains the normal harness lanes, rather than
    # assuming that the terminal public directory owns those earlier logs.
    representative_lane_logs = [
        candidate / "logs/constellation-verify.log"
        for candidate in log_root.iterdir()
        if candidate.is_dir()
        and (candidate / "logs/constellation-verify.log").is_file()
    ]
    if len(representative_lane_logs) != 1:
        raise SystemExit(
            "hostile fixture must retain exactly one representative supervised lane log"
        )
    representative_lane_log = representative_lane_logs[0]
    supervisor_rows = [
        line.removeprefix("supervisor_result_json=")
        for line in representative_lane_log.read_text(
            encoding="utf-8", errors="replace"
        ).splitlines()
        if line.startswith("supervisor_result_json=")
    ]
    if len(supervisor_rows) != 1:
        raise SystemExit(
            "representative supervised lane must retain exactly one result row"
        )
    representative_supervisor = json.loads(supervisor_rows[0])
    print(
        "nested-configured-lane-timeout-seconds="
        f"{representative_supervisor['configured_lane_timeout_seconds']}",
        flush=True,
    )
    consolidation_log = published_directories[0] / "logs/xtask-check-consolidation.log"
    if consolidation_log.is_file():
        for line in consolidation_log.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("fake-check-consolidation=") or line.startswith(
                "fake-consolidation-mutation="
            ):
                print(line, flush=True)

if mode in {"deadline-expiry", "postpublication-deadline-expiry"}:
    expected_rc = 124
    expected_statuses = ["INCOMPLETE"]
else:
    expected_rc = 1
    expected_statuses = ["FAIL"]
if result.returncode != expected_rc:
    raise SystemExit(
        f"hostile normal harness returned {result.returncode}, expected {expected_rc}"
    )
if public_statuses != expected_statuses or success_count != 0:
    raise SystemExit(
        "hostile normal harness published an unexpected terminal set: "
        f"{public_statuses!r}"
    )
raise SystemExit(result.returncode)
PY
}

export -f sha256_file closure_root_preflight_command bounded_snapshot_append
export -f append_snapshot_scope_hash capture_snapshot_command
export -f source_manifest_membership_command checker_smoke_sentinel_command
export -f wrapper_signal_self_test_command
export -f boundary_signal_injection_self_test_command
export -f publication_gap_mutation_self_test_command
export -f verifier_toctou_mutation_self_test_command
export -f publication_destination_collision_self_test_command
export -f normal_harness_hostile_self_test_command
export -f lane_log_cap_is_valid
export HARD_MAX_LANE_LOG_BYTES
export PROOF_MAX_SNAPSHOT_BYTES
export WRAPPER_SIGNAL_READINESS_TIMEOUT_SECONDS
export WRAPPER_SIGNAL_TERMINATION_TIMEOUT_SECONDS
export NORMAL_HARNESS_NESTED_LANE_TIMEOUT_SECONDS
export NORMAL_HARNESS_NESTED_RUN_TIMEOUT_SECONDS

capture_constellation_snapshot() { # lane destination
  local lane="$1" destination="$2" candidate failed_publication
  candidate="$(mktemp "$LOG_ROOT/.${lane}-candidate-XXXXXXXX")"
  LANE_COMMIT_CRITICAL=1
  run_lane "$lane" "source-provenance-snapshot" \
    bash -c 'capture_snapshot_command "$@"' _ \
    "$destination" "$candidate" "$HEAD_SHA" "${SCOPE_PATHS[@]}"
  if [[ "$LAST_LANE_STATUS" != PASS ]]; then
    if [[ -f "$destination" ]]; then
      failed_publication="${candidate}.published-before-failure"
      mv "$destination" "$failed_publication"
      candidate="$failed_publication"
    fi
    printf 'failed snapshot candidate retained outside proof bundle: %s\n' \
      "$candidate" >&2
  fi
  LANE_COMMIT_CRITICAL=0
  if [[ "$WRAPPER_SIGNAL" != 0 ]]; then
    exit_for_wrapper_signal "$WRAPPER_SIGNAL"
  fi
}

parse_verified_binding() { # verifier-output
  local output="$1" line count=0 dev="" ino="" commitment=""
  local binding_pattern='^proof-bundle-binding: dev=([0-9]+) ino=([0-9]+) commitment_sha256=([0-9a-f]{64})$'
  while IFS= read -r line || [[ -n "$line" ]]; do
    if [[ "$line" =~ $binding_pattern ]]; then
      count=$((count + 1))
      dev="${BASH_REMATCH[1]}"
      ino="${BASH_REMATCH[2]}"
      commitment="${BASH_REMATCH[3]}"
    fi
  done <<<"$output"
  if [[ "$count" != 1 ]]; then
    printf '%s\n' \
      'verified bundle output lacks exactly one publication binding' >&2
    return 1
  fi
  printf '%s %s %s\n' "$dev" "$ino" "$commitment"
}

# shellcheck disable=SC2329 # Exported into the bounded publication finalizer.
publish_bound_bundle() { # candidate published expected-dev expected-ino expected-commitment success-deadline
  python3 - "$1" "$2" "$3" "$4" "$5" "$6" \
    "$PROOF_MAX_BUNDLE_ENTRIES" "$PROOF_MAX_LANE_LOG_BYTES" \
    "$PROOF_MAX_TOTAL_LOG_BYTES" "$PROOF_MAX_VERDICTS_BYTES" \
    "$PROOF_MAX_SNAPSHOT_BYTES" "$PROOF_MAX_SUMMARY_BYTES" \
    "${FSIM_EULER_DISC_E2E_SELF_TEST_FAIL_PARENT_FSYNC:-0}" \
    "${FSIM_EULER_DISC_E2E_SELF_TEST_PUBLICATION_GAP_MARKER:-}" \
    "${FSIM_EULER_DISC_E2E_SELF_TEST_PUBLICATION_COLLISION_MARKER:-}" <<'PY'
import ctypes
import errno
import hashlib
import json
import os
import pathlib
import secrets
import stat
import sys
import time

(
    candidate_text,
    published_text,
    expected_dev_text,
    expected_ino_text,
    expected_commitment,
    success_deadline_text,
    max_entries_text,
    max_file_text,
    max_logs_text,
    max_verdicts_text,
    max_snapshots_text,
    max_summary_text,
    fail_parent_fsync_text,
    gap_marker_text,
    collision_marker_text,
) = sys.argv[1:]
candidate = pathlib.Path(candidate_text)
published = pathlib.Path(published_text)
expected_identity = (int(expected_dev_text), int(expected_ino_text))
success_deadline_ns = int(success_deadline_text)
max_entries = int(max_entries_text)
max_file = int(max_file_text)
max_total = (
    int(max_logs_text)
    + 2 * int(max_verdicts_text)
    + 2 * int(max_snapshots_text)
    + int(max_summary_text)
)
fail_parent_fsync = fail_parent_fsync_text == "1"
directory_flags = (
    os.O_RDONLY
    | getattr(os, "O_CLOEXEC", 0)
    | getattr(os, "O_DIRECTORY", 0)
    | getattr(os, "O_NOFOLLOW", 0)
)
file_flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)

def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)

def publish_regular_marker(marker, payload):
    """Atomically expose a fully written, identity-stable test marker."""
    require(hasattr(os, "O_NOFOLLOW"),
            "platform lacks no-follow marker publication")
    marker.parent.mkdir(parents=True, exist_ok=True)
    prepared = marker.with_name(f".{marker.name}.prepared-{os.getpid()}")
    flags = (
        os.O_WRONLY
        | os.O_CREAT
        | os.O_EXCL
        | os.O_NOFOLLOW
        | getattr(os, "O_CLOEXEC", 0)
    )
    descriptor = os.open(prepared, flags, 0o600)
    try:
        remaining = memoryview(payload)
        while remaining:
            written = os.write(descriptor, remaining)
            require(written > 0, "prepared marker write made no progress")
            remaining = remaining[written:]
        os.fsync(descriptor)
        prepared_status = os.fstat(descriptor)
        require(stat.S_ISREG(prepared_status.st_mode),
                "prepared marker is not a regular file")
        require(prepared_status.st_size == len(payload),
                "prepared marker has an incomplete payload")
        require(prepared_status.st_nlink == 1,
                "prepared marker acquired an unexpected link")
    finally:
        os.close(descriptor)
    os.link(prepared, marker, follow_symlinks=False)
    published_status = marker.lstat()
    require(stat.S_ISREG(published_status.st_mode),
            "published marker is not a regular file")
    require((published_status.st_dev, published_status.st_ino)
            == (prepared_status.st_dev, prepared_status.st_ino),
            "published marker does not bind the prepared inode")
    require(published_status.st_size == len(payload),
            "published marker has an incomplete payload")
    require(published_status.st_nlink == 2,
            "published marker has an unexpected link count")

def require_before_success_deadline(stage):
    if success_deadline_ns == 0:
        return
    observed = time.monotonic_ns()
    if observed >= success_deadline_ns:
        print(
            "aggregate success deadline expired "
            f"at stage={stage} deadline={success_deadline_ns} observed={observed}",
            file=sys.stderr,
        )
        raise SystemExit(124)

def rename_noreplace(source, destination):
    """Atomically publish without replacing any destination entry."""
    libc = ctypes.CDLL(None, use_errno=True)
    source_bytes = os.fsencode(source)
    destination_bytes = os.fsencode(destination)
    if sys.platform == "darwin":
        try:
            renamex_np = libc.renamex_np
        except AttributeError as error:
            raise RuntimeError(
                "atomic no-clobber publication is unavailable: renamex_np missing"
            ) from error
        renamex_np.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint]
        renamex_np.restype = ctypes.c_int
        result = renamex_np(source_bytes, destination_bytes, 0x00000004)
    elif sys.platform.startswith("linux"):
        try:
            renameat2 = libc.renameat2
        except AttributeError as error:
            raise RuntimeError(
                "atomic no-clobber publication is unavailable: renameat2 missing"
            ) from error
        renameat2.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        ]
        renameat2.restype = ctypes.c_int
        result = renameat2(-100, source_bytes, -100, destination_bytes, 0x00000001)
    else:
        raise RuntimeError(
            f"atomic no-clobber publication is unsupported on {sys.platform}"
        )
    if result == 0:
        return
    error_number = ctypes.get_errno()
    if error_number in {errno.EEXIST, errno.ENOTEMPTY}:
        raise RuntimeError(
            "atomic no-clobber publication refused existing destination"
        )
    raise OSError(
        error_number,
        f"atomic no-clobber publication failed: {os.strerror(error_number)}",
    )

def rename_exchange_at(source_directory, source_name, destination_directory, destination_name):
    """Atomically exchange two names beneath already-open directories."""
    libc = ctypes.CDLL(None, use_errno=True)
    source_bytes = os.fsencode(source_name)
    destination_bytes = os.fsencode(destination_name)
    if sys.platform == "darwin":
        try:
            renameatx_np = libc.renameatx_np
        except AttributeError as error:
            raise RuntimeError(
                "atomic identity-bound retraction is unavailable: renameatx_np missing"
            ) from error
        renameatx_np.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        ]
        renameatx_np.restype = ctypes.c_int
        result = renameatx_np(
            source_directory,
            source_bytes,
            destination_directory,
            destination_bytes,
            0x00000002,
        )
    elif sys.platform.startswith("linux"):
        try:
            renameat2 = libc.renameat2
        except AttributeError as error:
            raise RuntimeError(
                "atomic identity-bound retraction is unavailable: renameat2 missing"
            ) from error
        renameat2.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        ]
        renameat2.restype = ctypes.c_int
        result = renameat2(
            source_directory,
            source_bytes,
            destination_directory,
            destination_bytes,
            0x00000002,
        )
    else:
        raise RuntimeError(
            f"atomic identity-bound retraction is unsupported on {sys.platform}"
        )
    if result != 0:
        error_number = ctypes.get_errno()
        raise OSError(
            error_number,
            f"atomic identity-bound retraction failed: {os.strerror(error_number)}",
        )

def require_atomic_retraction_support(parent):
    """Exercise the exact same-filesystem exchange primitive before publication."""
    parent_descriptor = os.open(parent, directory_flags)
    preflight_descriptor = None
    left_descriptor = None
    right_descriptor = None
    try:
        for _ in range(128):
            preflight_name = (
                f".retraction-exchange-preflight-{secrets.token_hex(8)}"
            )
            try:
                os.mkdir(preflight_name, mode=0o700, dir_fd=parent_descriptor)
                break
            except FileExistsError:
                continue
        else:
            raise RuntimeError(
                "could not allocate a unique retraction-exchange preflight directory"
            )
        preflight = parent / preflight_name
        print(
            f"retraction-exchange-preflight-retained: {preflight}",
            file=sys.stderr,
        )
        preflight_descriptor = os.open(
            preflight_name,
            directory_flags,
            dir_fd=parent_descriptor,
        )
        os.mkdir("left", mode=0o700, dir_fd=preflight_descriptor)
        os.mkdir("right", mode=0o700, dir_fd=preflight_descriptor)
        left_descriptor = os.open(
            "left", directory_flags, dir_fd=preflight_descriptor
        )
        right_descriptor = os.open(
            "right", directory_flags, dir_fd=preflight_descriptor
        )
        left_identity = (
            os.fstat(left_descriptor).st_dev,
            os.fstat(left_descriptor).st_ino,
        )
        right_identity = (
            os.fstat(right_descriptor).st_dev,
            os.fstat(right_descriptor).st_ino,
        )
        rename_exchange_at(preflight_descriptor, "left", preflight_descriptor, "right")
        exchanged_left = os.stat(
            "left", dir_fd=preflight_descriptor, follow_symlinks=False
        )
        exchanged_right = os.stat(
            "right", dir_fd=preflight_descriptor, follow_symlinks=False
        )
        require(
            (exchanged_left.st_dev, exchanged_left.st_ino) == right_identity
            and (exchanged_right.st_dev, exchanged_right.st_ino) == left_identity,
            "atomic retraction exchange preflight returned inconsistent identities",
        )
    finally:
        for descriptor in (
            right_descriptor,
            left_descriptor,
            preflight_descriptor,
            parent_descriptor,
        ):
            if descriptor is not None:
                os.close(descriptor)

def metadata(value):
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )

def commitment(root_descriptor):
    rows = []
    pending = [(os.dup(root_descriptor), "")]
    entry_count = 0
    total_bytes = 0
    while pending:
        directory, prefix = pending.pop()
        try:
            entry_names = []
            with os.scandir(directory) as iterator:
                for entry in iterator:
                    entry_count += 1
                    require(
                        entry_count <= max_entries,
                        "publication inventory bound exceeded",
                    )
                    entry_names.append(entry.name)
            entry_names.sort()
            for entry_name in entry_names:
                relative = f"{prefix}/{entry_name}" if prefix else entry_name
                entry_stat = os.stat(
                    entry_name,
                    dir_fd=directory,
                    follow_symlinks=False,
                )
                require(not stat.S_ISLNK(entry_stat.st_mode), f"publication saw symlink {relative}")
                if stat.S_ISDIR(entry_stat.st_mode):
                    child = os.open(entry_name, directory_flags, dir_fd=directory)
                    child_stat = os.fstat(child)
                    require(metadata(child_stat) == metadata(entry_stat),
                            f"directory changed while binding {relative}")
                    rows.append(["directory", relative])
                    pending.append((child, relative))
                elif stat.S_ISREG(entry_stat.st_mode):
                    require(entry_stat.st_nlink == 1,
                            f"publication saw multiply linked file {relative}")
                    require(entry_stat.st_size <= max_file,
                            f"publication file bound exceeded for {relative}")
                    total_bytes += entry_stat.st_size
                    require(total_bytes <= max_total, "publication aggregate byte bound exceeded")
                    descriptor = os.open(entry_name, file_flags, dir_fd=directory)
                    try:
                        opened = os.fstat(descriptor)
                        require(metadata(opened) == metadata(entry_stat),
                                f"file changed while opening {relative}")
                        digest = hashlib.sha256()
                        observed = 0
                        while True:
                            chunk = os.read(descriptor, 65536)
                            if not chunk:
                                break
                            observed += len(chunk)
                            require(observed <= max_file,
                                    f"publication file grew beyond bound for {relative}")
                            digest.update(chunk)
                        require(observed == opened.st_size,
                                f"publication did not read all bytes for {relative}")
                        require(metadata(os.fstat(descriptor)) == metadata(opened),
                                f"file changed while binding {relative}")
                    finally:
                        os.close(descriptor)
                    rows.append(["file", relative, observed, digest.hexdigest()])
                else:
                    raise RuntimeError(f"publication saw non-regular entry {relative}")
        finally:
            os.close(directory)
    encoded = json.dumps(
        sorted(rows, key=lambda row: (row[1], row[0])),
        ensure_ascii=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()

require(candidate.parent.absolute() == published.parent.absolute(),
        "publication must remain in one parent directory")
require(not candidate.is_symlink() and candidate.is_dir(),
        "candidate proof bundle is not a real directory")
require(not os.path.lexists(published), "refusing to replace an existing proof-bundle path")
root_descriptor = os.open(candidate, directory_flags)
renamed = False
try:
    opened = os.fstat(root_descriptor)
    require((opened.st_dev, opened.st_ino) == expected_identity,
            "candidate directory identity differs from verified identity")
    require_before_success_deadline("before-publication-binding")
    require(commitment(root_descriptor) == expected_commitment,
            "candidate bytes differ from verified publication commitment")
    require_before_success_deadline("before-retraction-exchange-preflight")
    require_atomic_retraction_support(published.parent)
    require_before_success_deadline("after-retraction-exchange-preflight")

    if gap_marker_text and not os.path.lexists(gap_marker_text):
        marker = pathlib.Path(gap_marker_text)
        publish_regular_marker(
            marker,
            f"candidate={candidate}\n".encode("utf-8"),
        )
        print(f"publication-gap-marker-ready: {marker}", file=sys.stderr, flush=True)
        time.sleep(3)

    if collision_marker_text and not os.path.lexists(collision_marker_text):
        collision_marker = pathlib.Path(collision_marker_text)
        publish_regular_marker(
            collision_marker,
            f"candidate={candidate}\npublished={published}\n".encode("utf-8"),
        )
        print(
            f"publication-collision-marker-ready: {collision_marker}",
            file=sys.stderr,
            flush=True,
        )
        time.sleep(3)

    require(commitment(root_descriptor) == expected_commitment,
            "candidate bytes changed in the verifier-to-publication gap")
    require((os.fstat(root_descriptor).st_dev, os.fstat(root_descriptor).st_ino)
            == expected_identity,
            "candidate directory identity changed before publication")
    require_before_success_deadline("immediately-before-publication-rename")
    rename_noreplace(candidate, published)
    renamed = True
    published_stat = published.lstat()
    require(stat.S_ISDIR(published_stat.st_mode)
            and (published_stat.st_dev, published_stat.st_ino) == expected_identity,
            "published path does not name the verified directory")
    require(commitment(root_descriptor) == expected_commitment,
            "published bytes differ from the verified publication commitment")

    try:
        parent_descriptor = os.open(published.parent, os.O_RDONLY)
        try:
            if fail_parent_fsync:
                raise OSError(errno.EIO, "self-test injected parent-directory fsync failure")
            os.fsync(parent_descriptor)
        finally:
            os.close(parent_descriptor)
    except OSError as error:
        print(
            "proof-bundle parent-directory fsync failed after atomic publication: "
            f"{published}: {error}",
            file=sys.stderr,
        )
    require_before_success_deadline("after-publication-fsync")
except BaseException:
    if renamed:
        print(
            "publication failed after atomic rename; the supervising parent "
            "must perform identity-bound retraction",
            file=sys.stderr,
        )
    raise
finally:
    os.close(root_descriptor)
PY
}

# Publication runs in a bounded child so an over-deadline commitment, rename,
# or parent-directory sync cannot pin the producer inside success finalization.
# The parent remains responsible for retracting a path if the child was stopped
# after its atomic rename.
export -f publish_bound_bundle
export PROOF_MAX_BUNDLE_ENTRIES PROOF_MAX_LANE_LOG_BYTES
export PROOF_MAX_TOTAL_LOG_BYTES PROOF_MAX_VERDICTS_BYTES
export PROOF_MAX_SNAPSHOT_BYTES PROOF_MAX_SUMMARY_BYTES

retract_published_bundle() { # reason expected-dev expected-ino
  local reason="$1" expected_dev="$2" expected_ino="$3"
  local retracted_path="" retraction_rc retraction_log
  retraction_log="$(
    mktemp "$LOG_ROOT/.identity-bound-retraction.log-XXXXXXXX"
  )"
  RETRACTION_LAST_LOG="$retraction_log"
  if retracted_path="$(
    python3 - "$PUBLISHED_DIR" "$LOG_ROOT" "$reason" \
      "$expected_dev" "$expected_ino" 2>>"$retraction_log" <<'PY'
import ctypes
import errno
import os
import pathlib
import secrets
import stat
import sys

published_text, log_root_text, reason, expected_dev_text, expected_ino_text = sys.argv[1:]
published = pathlib.Path(published_text)
log_root = pathlib.Path(log_root_text)
expected_identity = (int(expected_dev_text), int(expected_ino_text))
directory_flags = (
    os.O_RDONLY
    | getattr(os, "O_CLOEXEC", 0)
    | getattr(os, "O_DIRECTORY", 0)
    | getattr(os, "O_NOFOLLOW", 0)
)

def identity(metadata):
    return metadata.st_dev, metadata.st_ino

def fail(detail):
    print(
        f"published candidate retraction failed: reason={reason} "
        f"path={published}: {detail}",
        file=sys.stderr,
    )
    raise SystemExit(1)

def exchange(source_directory, source_name, destination_directory, destination_name):
    libc = ctypes.CDLL(None, use_errno=True)
    source_bytes = os.fsencode(source_name)
    destination_bytes = os.fsencode(destination_name)
    if sys.platform == "darwin":
        try:
            renameatx_np = libc.renameatx_np
        except AttributeError as error:
            raise RuntimeError("renameatx_np is unavailable") from error
        renameatx_np.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        ]
        renameatx_np.restype = ctypes.c_int
        result = renameatx_np(
            source_directory,
            source_bytes,
            destination_directory,
            destination_bytes,
            0x00000002,
        )
    elif sys.platform.startswith("linux"):
        try:
            renameat2 = libc.renameat2
        except AttributeError as error:
            raise RuntimeError("renameat2 is unavailable") from error
        renameat2.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        ]
        renameat2.restype = ctypes.c_int
        result = renameat2(
            source_directory,
            source_bytes,
            destination_directory,
            destination_bytes,
            0x00000002,
        )
    else:
        raise RuntimeError(f"unsupported platform {sys.platform}")
    if result != 0:
        error_number = ctypes.get_errno()
        raise OSError(error_number, os.strerror(error_number))

if os.path.abspath(published.parent) != os.path.abspath(log_root):
    fail("publication parent and configured log root differ")

parent_descriptor = None
expected_descriptor = None
placeholder_descriptor = None
try:
    parent_descriptor = os.open(log_root, directory_flags)
    try:
        expected_descriptor = os.open(
            published.name,
            directory_flags,
            dir_fd=parent_descriptor,
        )
    except FileNotFoundError:
        raise SystemExit(0)
    except OSError as error:
        if error.errno in {errno.ELOOP, errno.ENOTDIR}:
            raise SystemExit(2) from error
        raise

    opened = os.fstat(expected_descriptor)
    if not stat.S_ISDIR(opened.st_mode) or identity(opened) != expected_identity:
        raise SystemExit(2)

    for _ in range(128):
        rejected_name = f".publication-retracted-{secrets.token_hex(8)}"
        try:
            os.mkdir(rejected_name, mode=0o700, dir_fd=parent_descriptor)
            break
        except FileExistsError:
            continue
    else:
        fail("could not allocate a unique retained directory")

    retained = log_root / rejected_name
    placeholder_descriptor = os.open(
        rejected_name,
        directory_flags,
        dir_fd=parent_descriptor,
    )
    placeholder_identity = identity(os.fstat(placeholder_descriptor))

    exchange(parent_descriptor, published.name, parent_descriptor, rejected_name)
    if os.environ.get(
        "FSIM_EULER_DISC_E2E_SELF_TEST_REPLACE_PUBLIC_TOMBSTONE"
    ) == "1":
        displaced_name = (
            f".self-test-displaced-public-tombstone-{secrets.token_hex(8)}"
        )
        os.rename(
            published.name,
            displaced_name,
            src_dir_fd=parent_descriptor,
            dst_dir_fd=parent_descriptor,
        )
        os.mkdir(published.name, mode=0o700, dir_fd=parent_descriptor)
        print(
            "self-test replaced public tombstone after exchange: "
            f"displaced={log_root / displaced_name} replacement={published}",
            file=sys.stderr,
        )
    hidden = os.stat(
        rejected_name,
        dir_fd=parent_descriptor,
        follow_symlinks=False,
    )
    public = os.stat(
        published.name,
        dir_fd=parent_descriptor,
        follow_symlinks=False,
    )
    if identity(hidden) != expected_identity:
        # Restore only if both names still designate the exact pair observed
        # after the exchange. Under continuing hostile namespace mutation there
        # is no safe path-only rollback; that case fails closed in place.
        rollback_hidden = os.stat(
            rejected_name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        rollback_public = os.stat(
            published.name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        if (
            identity(rollback_hidden) == identity(hidden)
            and identity(rollback_public) == identity(public)
            and identity(public) == placeholder_identity
        ):
            exchange(
                parent_descriptor,
                published.name,
                parent_descriptor,
                rejected_name,
            )
            restored = os.stat(
                published.name,
                dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
            if identity(restored) != identity(hidden):
                fail("identity-mismatch rollback did not restore the exchanged entry")
            raise SystemExit(2)
        fail("namespace changed during identity-bound exchange")

    if identity(public) != placeholder_identity:
        fail(
            "expected public tombstone identity was replaced after exchange; "
            f"namespace custody is not claimed; retained={retained}"
        )
    retained_stat = os.stat(
        rejected_name,
        dir_fd=parent_descriptor,
        follow_symlinks=False,
    )
    if identity(retained_stat) != expected_identity:
        fail("retained bundle identity changed after exchange")
    print(retained)
except SystemExit:
    raise
except BaseException as error:
    fail(f"{type(error).__name__}: {error}")
finally:
    for descriptor in (
        placeholder_descriptor,
        expected_descriptor,
        parent_descriptor,
    ):
        if descriptor is not None:
            os.close(descriptor)
PY
  )"; then
    if [[ -n "$retracted_path" ]]; then
      local retracted_root_name retracted_token
      LOG_DIR="$retracted_path"
      VERDICTS="$LOG_DIR/verdicts-prefix.jsonl"
      SUMMARY="$LOG_DIR/summary.json"
      retracted_root_name="${retracted_path##*/}"
      retracted_token="${retracted_root_name#.publication-retracted-}"
      PUBLISHED_DIR="$LOG_ROOT/retracted-${retracted_token}-terminal"
      printf 'published candidate retracted before terminal transition: reason=%s path=%s\n' \
        "$reason" "$LOG_DIR" >>"$retraction_log"
      queue_success_finalizer_report \
        "identity-bound-retraction" 0 "$retraction_log"
    fi
    return 0
  else
    retraction_rc=$?
  fi
  if [[ "$retraction_rc" == 2 ]]; then
    printf 'published candidate retraction refused an identity mismatch: reason=%s path=%s expected-dev=%s expected-ino=%s\n' \
      "$reason" "$PUBLISHED_DIR" "$expected_dev" "$expected_ino" \
      >>"$retraction_log"
    queue_success_finalizer_report \
      "identity-bound-retraction" 2 "$retraction_log"
    return 2
  fi
  printf 'published candidate retraction could not establish an identity-bound exchange: reason=%s path=%s exit=%s\n' \
    "$reason" "$PUBLISHED_DIR" "$retraction_rc" >>"$retraction_log"
  queue_success_finalizer_report \
    "identity-bound-retraction" "$retraction_rc" "$retraction_log"
  return 1
}

mark_retraction_integrity_failure() { # transition-context
  RETRACTION_INTEGRITY_CONTEXT="$1"
  # Establish the recovery guard before releasing SEALING. A signal delivered
  # anywhere after this point may still be latched for the retained diagnostic,
  # but the trap cannot replace authoritative integrity exit 125 with an
  # INTERRUPTED transition. The EXIT trap snapshots and clears that latch only
  # after it has ignored further wrapper signals.
  ((
    RETRACTION_INTEGRITY_RECOVERY = 1,
    CANDIDATE_REJECTED = 1,
    SEALING = 0,
    1
  ))
}

# shellcheck disable=SC2329 # Invoked from the EXIT-trap sealing path.
start_fresh_incomplete_candidate() {
  local rejected="$LOG_DIR"
  LOG_DIR="$(mktemp -d "$LOG_ROOT/.candidate-${HEAD_SHA:0:12}-${PROFILE}-incomplete-XXXXXXXX")"
  PUBLISHED_DIR="$LOG_ROOT/${LOG_DIR##*/.candidate-}"
  mkdir -p "$LOG_DIR/logs"
  VERDICTS="$LOG_DIR/verdicts-prefix.jsonl"
  SUMMARY="$LOG_DIR/summary.json"
  : >"$VERDICTS"
  CHECKS=0
  FAILURES=0
  RETAINED_LOG_BYTES_TOTAL=0
  SNAPSHOT_BEFORE=""
  SNAPSHOT_AFTER=""
  CANDIDATE_REJECTED=0
  printf 'rejected staged bundle retained outside terminal proof: %s\n' \
    "$rejected" >&2
  printf 'fresh incomplete proof candidate: %s\n' "$LOG_DIR" >&2
}

is_success_terminal_status() { # status
  case "$1" in
    FOCUSED_PASS|READY_FOR_DSR|SELF_TEST_PASS) return 0 ;;
    *) return 1 ;;
  esac
}

require_success_deadline() { # deadline-monotonic-ns stage
  local deadline_detail deadline_rc deadline_log
  if deadline_detail="$(python3 - "$1" "$2" 2>&1 <<'PY'
import sys
import time

deadline = int(sys.argv[1])
stage = sys.argv[2]
observed = time.monotonic_ns()
if observed >= deadline:
    print(
        "aggregate success deadline expired "
        f"at stage={stage} deadline={deadline} observed={observed}",
        file=sys.stderr,
    )
    raise SystemExit(124)
PY
  )"; then
    return 0
  else
    deadline_rc=$?
  fi
  deadline_log="$(
    mktemp "$LOG_ROOT/.success-deadline-${2}.log-XXXXXXXX"
  )"
  printf '%s\n' "$deadline_detail" >"$deadline_log"
  queue_success_finalizer_report \
    "aggregate-success-deadline-${2}" "$deadline_rc" "$deadline_log"
  return "$deadline_rc"
}

# This fail-only environment hook is used by the no-Cargo self-test. It moves
# the remaining aggregate success deadline close to the current instant, then
# waits beyond it so the real prepublication guard must refuse the candidate.
inject_prepublication_deadline_delay() { # current-deadline-monotonic-ns
  local injected_deadline
  read -r injected_deadline < <(
    python3 - <<'PY'
import time

print(time.monotonic_ns() + 100_000_000)
PY
  )
  RUN_DEADLINE_MONOTONIC_NS="$injected_deadline"
  printf 'self-test prepublication deadline delay started: prior=%s injected=%s\n' \
    "$1" "$RUN_DEADLINE_MONOTONIC_NS" >&2
  python3 - <<'PY'
import time

time.sleep(0.25)
PY
  printf '%s\n' 'self-test prepublication deadline delay complete' >&2
}

FINALIZATION_LAST_LOG=""

queue_success_finalizer_report() { # label exit-code retained-log
  FINALIZATION_REPORT_LABELS+=("$1")
  FINALIZATION_REPORT_RCS+=("$2")
  FINALIZATION_REPORT_LOGS+=("$3")
}

# Human-facing diagnostic output is deliberately postponed until a terminal
# bundle is already sealed. A blocked or broken stderr can therefore neither
# hold a renamed success candidate public nor prevent failure sealing.
flush_success_finalizer_reports() {
  local label finalization_rc finalization_log finalization_line
  if [[ "$SEALED" != 1 ]]; then
    return 0
  fi
  while ((FINALIZATION_REPORT_CURSOR < ${#FINALIZATION_REPORT_LABELS[@]})); do
    label="${FINALIZATION_REPORT_LABELS[$FINALIZATION_REPORT_CURSOR]}"
    finalization_rc="${FINALIZATION_REPORT_RCS[$FINALIZATION_REPORT_CURSOR]}"
    finalization_log="${FINALIZATION_REPORT_LOGS[$FINALIZATION_REPORT_CURSOR]}"
    FINALIZATION_REPORT_CURSOR=$((FINALIZATION_REPORT_CURSOR + 1))
    printf 'success finalization log retained: label=%s exit=%s path=%s\n' \
      "$label" "$finalization_rc" "$finalization_log" >&2 || true
    if [[ "$label" == publication && "$finalization_rc" == 0 \
        && -r "$finalization_log" ]]; then
      while IFS= read -r finalization_line || [[ -n "$finalization_line" ]]; do
        case "$finalization_line" in
          'proof-bundle parent-directory fsync failed after atomic publication:'*)
            printf '%s\n' "$finalization_line" >&2 || true
            ;;
          'retraction-exchange-preflight-retained:'*)
            printf '%s\n' "$finalization_line" >&2 || true
            ;;
        esac
      done <"$finalization_log" || true
    fi
    if [[ "$label" == identity-bound-retraction \
        || "$label" == postpublication-self-test-arm \
        || "$label" == postpublication-self-test-delay ]]; then
      if [[ -r "$finalization_log" ]]; then
        while IFS= read -r finalization_line || [[ -n "$finalization_line" ]]; do
          printf '%s\n' "$finalization_line" >&2 || true
        done <"$finalization_log" || true
      fi
    elif [[ "$finalization_rc" != 0 ]]; then
      printf '%s\n' \
        "--- failed success-finalization log replay: $finalization_log ---" \
        >&2 || true
      if [[ -r "$finalization_log" ]]; then
        while IFS= read -r finalization_line || [[ -n "$finalization_line" ]]; do
          printf '%s\n' "$finalization_line" >&2 || true
        done <"$finalization_log" || true
      fi
      printf '%s\n' '--- end failed success-finalization log replay ---' \
        >&2 || true
    fi
  done
  return 0
}

run_success_finalizer() { # label success-deadline-monotonic-ns command...
  local label="$1" success_deadline_ns="$2" supervisor_deadline_ns finalization_rc
  shift 2
  FINALIZATION_LAST_LOG="$(
    mktemp "$LOG_ROOT/.success-finalization-${label}.log-XXXXXXXX"
  )"
  supervisor_deadline_ns="$success_deadline_ns"
  if [[ "$supervisor_deadline_ns" == 0 ]]; then
    read -r supervisor_deadline_ns < <(
      python3 - "$LANE_TIMEOUT_SECONDS" <<'PY'
import sys
import time

print(time.monotonic_ns() + int(sys.argv[1]) * 1_000_000_000)
PY
    )
  fi
  {
    printf '%s\n' 'schema=frankensim.euler-disc-contract-e2e.success-finalization.v1'
    printf 'label=%s\nsuccess_deadline_monotonic_ns=%s\nargv=' \
      "$label" "$success_deadline_ns"
    printf '%q ' "$@"
    printf '\n'
  } >"$FINALIZATION_LAST_LOG"
  if run_bounded_command \
    "$FINALIZATION_LAST_LOG" "$LANE_TIMEOUT_SECONDS" "$supervisor_deadline_ns" \
    "$SUCCESS_FINALIZATION_LOG_MAX_BYTES" \
    "$SUCCESS_FINALIZATION_LOG_MAX_BYTES" \
    "$SUCCESS_FINALIZATION_TOTAL_LOG_MAX_BYTES" none "$@"; then
    finalization_rc=0
  else
    finalization_rc=$?
  fi
  # A command can exit just before the supervisor observes the aggregate
  # deadline. Admission still fails if the completed wrapper returns too late.
  if [[ "$success_deadline_ns" != 0 ]]; then
    if require_success_deadline "$success_deadline_ns" "after-${label}"; then
      :
    else
      finalization_rc=$?
    fi
  fi
  queue_success_finalizer_report \
    "$label" "$finalization_rc" "$FINALIZATION_LAST_LOG"
  return "$finalization_rc"
}

# This second fail-only hook arms expiry while the already-renamed success
# candidate is undergoing its independently supervised verifier pass. The
# marker lives outside the proof bundle and is retained for diagnosis.
arm_postpublication_deadline_delay() {
  local marker_root diagnostic_log
  marker_root="$(mktemp -d "$LOG_ROOT/.postpublication-deadline-XXXXXXXX")"
  POSTPUBLICATION_VERIFY_READY_MARKER="$marker_root/verifier-ready"
  read -r RUN_DEADLINE_MONOTONIC_NS < <(
    python3 - <<'PY'
import time

print(time.monotonic_ns() + 500_000_000)
PY
  )
  diagnostic_log="$(
    mktemp "$LOG_ROOT/.postpublication-self-test-arm.log-XXXXXXXX"
  )"
  printf 'self-test postpublication deadline armed: deadline=%s marker=%s\n' \
    "$RUN_DEADLINE_MONOTONIC_NS" "$POSTPUBLICATION_VERIFY_READY_MARKER" \
    >"$diagnostic_log"
  queue_success_finalizer_report \
    "postpublication-self-test-arm" 0 "$diagnostic_log"
}

write_summary_and_seal() { # status coverage ready proof-scope provenance-state exit-code
  local status="$1" coverage="$2" ready="$3" proof_scope="$4" provenance="$5"
  local terminal_exit_code="$6" success_deadline_ns=0 deadline_rc=0
  if [[ "$SEALED" == 1 ]]; then
    return
  fi
  if is_success_terminal_status "$status"; then
    success_deadline_ns="$RUN_DEADLINE_MONOTONIC_NS"
    if require_success_deadline \
      "$success_deadline_ns" "before-success-verification"; then
      :
    else
      deadline_rc=$?
      CANDIDATE_REJECTED=1
      printf 'success candidate rejected before verification: status=%s exit=%s\n' \
        "$status" "$deadline_rc" >&2
      return "$deadline_rc"
    fi
  fi
  SEALING=1
  local before_sha="" after_sha=""
  if [[ -n "$SNAPSHOT_BEFORE" && -f "$SNAPSHOT_BEFORE" ]]; then
    before_sha="$(sha256_file "$SNAPSHOT_BEFORE")"
  fi
  if [[ -n "$SNAPSHOT_AFTER" && -f "$SNAPSHOT_AFTER" ]]; then
    after_sha="$(sha256_file "$SNAPSHOT_AFTER")"
  fi
  local prefix_sha
  prefix_sha="$(sha256_file "$VERDICTS")"
  local candidate_summary="$LOG_DIR/summary.json"
  local candidate_verdicts="$LOG_DIR/verdicts.jsonl"
  if ! python3 - "$candidate_summary" "$candidate_verdicts" "$VERDICTS" \
    "$status" "$HEAD_SHA" "$HOST_ISA" \
    "$PROFILE" "$EXECUTOR_DECLARATION" "$CHECKS" "$FAILURES" "$coverage" \
    "$ready" "$proof_scope" "$provenance" "$before_sha" "$after_sha" \
    "$prefix_sha" "$terminal_exit_code" <<'PY'
import json
import os
import sys

(
    summary_path,
    verdicts,
    prefix_path,
    status,
    head,
    isa,
    profile,
    executor_declaration,
    checks,
    failures,
    coverage,
    ready,
    proof_scope,
    provenance,
    before_sha,
    after_sha,
    prefix_sha,
    terminal_exit_code,
) = sys.argv[1:]
seal = {
    "schema": "frankensim.euler-disc-contract-e2e.proof-seal.v2",
    "record_type": "proof-seal",
    "status": status,
    "checks": int(checks),
    "failures": int(failures),
    "head": head,
    "host_isa": isa,
    "profile": profile,
    "executor_declaration": executor_declaration,
    "executor_attestation": "caller-declared-unverified",
    "source_manifest_coverage": coverage,
    "candidate_ready_for_dsr": ready == "true",
    "dsr_status": "not-run-by-this-harness",
    "proof_scope": proof_scope,
    "provenance_state": provenance,
    "snapshot_before_sha256": before_sha,
    "snapshot_after_sha256": after_sha,
    "verdicts_prefix_sha256": prefix_sha,
    "verdicts": "verdicts.jsonl",
    "proof_seal_locator": "verdicts.jsonl#last-line",
    "terminal_exit_code": int(terminal_exit_code),
    "terminal": True,
    "no_claim": "Software/protocol structure only; no physical validation or emergent Euler-disc prediction.",
    "concurrency_no_claim": "Clean HEAD is observed at bookends; transient concurrent edits between observations are not excluded.",
}
encoded = (
    json.dumps(seal, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"
).encode("utf-8")
with open(prefix_path, "rb") as handle:
    prefix = handle.read()
with open(verdicts, "wb") as handle:
    handle.write(prefix)
    handle.write(encoded)
    handle.flush()
    os.fsync(handle.fileno())
with open(summary_path, "wb") as handle:
    handle.write(encoded)
    handle.flush()
    os.fsync(handle.fileno())
PY
  then
    CANDIDATE_REJECTED=1
    SEALING=0
    return 1
  fi
  if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_CORRUPT_CANDIDATE:-0}" == 1 \
      && "$status" == "SELF_TEST_PASS" \
      && "$CANDIDATE_CORRUPTION_INJECTED" == 0 ]]; then
    CANDIDATE_CORRUPTION_INJECTED=1
    printf '%s\n' ' ' >>"$candidate_summary"
    printf '%s\n' 'candidate corruption injected before verification' >&2
  fi
  local verification_output binding_dev binding_ino binding_commitment
  local verification_rc publication_rc retraction_rc post_verification_output
  local post_dev post_ino post_commitment postpublication_marker_log
  local -a post_verifier_command
  if run_success_finalizer \
    "prepublication-verification" "$success_deadline_ns" \
    "$BASH" \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$LOG_DIR"; then
    verification_output="$(<"$FINALIZATION_LAST_LOG")"
    if read -r binding_dev binding_ino binding_commitment \
      < <(parse_verified_binding "$verification_output"); then
      :
    else
      printf '%s\n' 'verified candidate did not yield a usable publication binding' >&2
      CANDIDATE_REJECTED=1
      SEALING=0
      return 1
    fi
  else
    verification_rc=$?
    if [[ "$WRAPPER_SIGNAL" != 0 && "$status" != "INTERRUPTED" ]]; then
      SEALING=0
      exit_for_wrapper_signal "$WRAPPER_SIGNAL"
    fi
    CANDIDATE_REJECTED=1
    SEALING=0
    return "$verification_rc"
  fi
  if [[ "$success_deadline_ns" != 0 ]]; then
    if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_PREPUBLICATION_DEADLINE_DELAY:-0}" == 1 \
        && "$PREPUBLICATION_DEADLINE_INJECTED" == 0 ]]; then
      PREPUBLICATION_DEADLINE_INJECTED=1
      inject_prepublication_deadline_delay "$success_deadline_ns"
      success_deadline_ns="$RUN_DEADLINE_MONOTONIC_NS"
    fi
    if require_success_deadline \
      "$success_deadline_ns" "before-success-publication"; then
      :
    else
      deadline_rc=$?
      CANDIDATE_REJECTED=1
      SEALING=0
      printf 'verified success candidate rejected before publication: status=%s exit=%s\n' \
        "$status" "$deadline_rc" >&2
      return "$deadline_rc"
    fi
  fi
  if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_PREPUBLICATION_SIGNAL:-0}" == 1 \
      && "$PREPUBLICATION_SIGNAL_INJECTED" == 0 \
      && "$status" != "INTERRUPTED" ]]; then
    PREPUBLICATION_SIGNAL_INJECTED=1
    WRAPPER_SIGNAL=15
  fi

  if [[ "$WRAPPER_SIGNAL" != 0 && "$status" != "INTERRUPTED" ]]; then
    SEALING=0
    exit_for_wrapper_signal "$WRAPPER_SIGNAL"
  fi
  if run_success_finalizer \
    "publication" "$success_deadline_ns" \
    bash -c 'publish_bound_bundle "$@"' _ \
    "$LOG_DIR" "$PUBLISHED_DIR" "$binding_dev" "$binding_ino" \
    "$binding_commitment" "$success_deadline_ns"; then
    LOG_DIR="$PUBLISHED_DIR"
    VERDICTS="$LOG_DIR/verdicts-prefix.jsonl"
    SUMMARY="$LOG_DIR/summary.json"
  else
    publication_rc=$?
    if retract_published_bundle \
      "bounded-publication-finalizer-exit-$publication_rc" \
      "$binding_dev" "$binding_ino"; then
      :
    else
      retraction_rc=$?
      # A different identity is the deliberately preserved destination-race
      # case, not our candidate and therefore not ours to move.
      if [[ "$retraction_rc" != 2 ]]; then
        mark_retraction_integrity_failure \
          "bounded-publication-finalizer-exit-$publication_rc"
        return 125
      fi
    fi
    if [[ "$WRAPPER_SIGNAL" != 0 && "$status" != "INTERRUPTED" ]]; then
      printf 'publication finalizer signal prevented or revoked candidate admission: signal=%s path=%s\n' \
        "$WRAPPER_SIGNAL" "$PUBLISHED_DIR" >&2
      SEALING=0
      exit_for_wrapper_signal "$WRAPPER_SIGNAL"
    fi
    CANDIDATE_REJECTED=1
    SEALING=0
    printf 'verified terminal bundle was not admitted at the publication path: %s\n' \
      "$LOG_DIR" >&2
    return "$publication_rc"
  fi

  if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_POSTPUBLICATION_SIGNAL:-0}" == 1 \
      && "$POSTPUBLICATION_SIGNAL_INJECTED" == 0 \
      && "$status" != "INTERRUPTED" ]]; then
    POSTPUBLICATION_SIGNAL_INJECTED=1
    WRAPPER_SIGNAL=15
  fi

  # Bash queues a handled signal while waiting for the foreground publication
  # subprocess. Revoke the just-published success before creating the exact
  # INTERRUPTED terminal state; publication-window HUP/INT/TERM is never ignored.
  if [[ "$WRAPPER_SIGNAL" != 0 && "$status" != "INTERRUPTED" ]]; then
    if ! retract_published_bundle \
      "wrapper-signal-$WRAPPER_SIGNAL" "$binding_dev" "$binding_ino"; then
      mark_retraction_integrity_failure "wrapper-signal-after-publication"
      return 125
    fi
    SEALING=0
    exit_for_wrapper_signal "$WRAPPER_SIGNAL"
  fi

  if [[ "$success_deadline_ns" != 0 ]]; then
    if require_success_deadline \
      "$success_deadline_ns" "before-postpublication-verification"; then
      :
    else
      deadline_rc=$?
      if ! retract_published_bundle \
        "aggregate-deadline-before-postverification" \
        "$binding_dev" "$binding_ino"; then
        mark_retraction_integrity_failure \
          "aggregate-deadline-before-postverification"
        return 125
      fi
      CANDIDATE_REJECTED=1
      SEALING=0
      printf 'published success candidate retracted before postverification: status=%s exit=%s\n' \
        "$status" "$deadline_rc" >&2
      return "$deadline_rc"
    fi
  fi

  post_verifier_command=(
    "$BASH"
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh"
    --verify-bundle "$LOG_DIR"
  )
  if [[ "$success_deadline_ns" != 0 \
      && "${FSIM_EULER_DISC_E2E_SELF_TEST_POSTPUBLICATION_DEADLINE_DELAY:-0}" == 1 \
      && "$POSTPUBLICATION_DEADLINE_INJECTED" == 0 ]]; then
    POSTPUBLICATION_DEADLINE_INJECTED=1
    arm_postpublication_deadline_delay
    success_deadline_ns="$RUN_DEADLINE_MONOTONIC_NS"
    post_verifier_command=(
      env
      FSIM_EULER_DISC_E2E_SELF_TEST_PAUSE_DURING_VERIFY=1
      "FSIM_EULER_DISC_E2E_SELF_TEST_VERIFY_READY_MARKER=$POSTPUBLICATION_VERIFY_READY_MARKER"
      "$BASH"
      "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh"
      --verify-bundle "$LOG_DIR"
    )
  fi

  if run_success_finalizer \
      "postpublication-verification" "$success_deadline_ns" \
      "${post_verifier_command[@]}"; then
    post_verification_output="$(<"$FINALIZATION_LAST_LOG")"
    if read -r post_dev post_ino post_commitment \
        < <(parse_verified_binding "$post_verification_output") \
        && [[ "$post_dev" == "$binding_dev" ]] \
        && [[ "$post_ino" == "$binding_ino" ]] \
        && [[ "$post_commitment" == "$binding_commitment" ]]; then
      :
    else
      publication_rc=1
    fi
  else
    publication_rc=$?
  fi
  if [[ -n "$POSTPUBLICATION_VERIFY_READY_MARKER" \
      && -f "$POSTPUBLICATION_VERIFY_READY_MARKER" ]]; then
    postpublication_marker_log="$(
      mktemp "$LOG_ROOT/.postpublication-self-test-delay.log-XXXXXXXX"
    )"
    printf 'self-test postpublication verifier delay observed: marker=%s\n' \
      "$POSTPUBLICATION_VERIFY_READY_MARKER" >"$postpublication_marker_log"
    queue_success_finalizer_report \
      "postpublication-self-test-delay" 0 "$postpublication_marker_log"
  fi
  if [[ "${publication_rc:-0}" == 0 ]]; then
    :
  else
    if ! retract_published_bundle \
      "post-publication-reverification-failed" \
      "$binding_dev" "$binding_ino"; then
      mark_retraction_integrity_failure \
        "post-publication-reverification-failed"
      return 125
    fi
    if [[ "$WRAPPER_SIGNAL" != 0 && "$status" != "INTERRUPTED" ]]; then
      SEALING=0
      exit_for_wrapper_signal "$WRAPPER_SIGNAL"
    fi
    CANDIDATE_REJECTED=1
    SEALING=0
    printf '%s\n' \
      'post-publication identity/byte/semantic binding could not be re-established' >&2
    return "$publication_rc"
  fi
  if [[ "$success_deadline_ns" != 0 ]]; then
    if require_success_deadline \
      "$success_deadline_ns" "after-postpublication-verification"; then
      :
    else
      deadline_rc=$?
      if ! retract_published_bundle \
        "aggregate-deadline-after-postverification" \
        "$binding_dev" "$binding_ino"; then
        mark_retraction_integrity_failure \
          "aggregate-deadline-after-postverification"
        return 125
      fi
      CANDIDATE_REJECTED=1
      SEALING=0
      printf 'published success candidate retracted after postverification: status=%s exit=%s\n' \
        "$status" "$deadline_rc" >&2
      return "$deadline_rc"
    fi
  fi
  if [[ "$WRAPPER_SIGNAL" != 0 && "$status" != "INTERRUPTED" ]]; then
    if ! retract_published_bundle \
      "wrapper-signal-$WRAPPER_SIGNAL-after-reverification" \
      "$binding_dev" "$binding_ino"; then
      mark_retraction_integrity_failure \
        "wrapper-signal-after-reverification"
      return 125
    fi
    SEALING=0
    exit_for_wrapper_signal "$WRAPPER_SIGNAL"
  fi
  if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_FINAL_COMMIT_SIGNAL:-0}" == 1 \
      && "$FINAL_COMMIT_SIGNAL_INJECTED" == 0 \
      && "$status" != "INTERRUPTED" ]]; then
    FINAL_COMMIT_SIGNAL_INJECTED=1
    WRAPPER_SIGNAL=15
  fi
  # This single arithmetic builtin is the success-commit point with respect to
  # Bash trap dispatch. A signal processed before it sees SEALING=1 and latches;
  # a signal processed after it sees the already-sealed success. The check just
  # below revokes any signal that was latched on the pre-commit side.
  ((SEALED = 1, SEALING = 0, 1))
  if [[ "$WRAPPER_SIGNAL" != 0 && "$status" != "INTERRUPTED" ]]; then
    SEALED=0
    SEALING=1
    if ! retract_published_bundle \
      "wrapper-signal-$WRAPPER_SIGNAL-at-success-commit" \
      "$binding_dev" "$binding_ino"; then
      mark_retraction_integrity_failure "wrapper-signal-at-success-commit"
      return 125
    fi
    CANDIDATE_REJECTED=1
    SEALING=0
    exit_for_wrapper_signal "$WRAPPER_SIGNAL"
  fi
  flush_success_finalizer_reports || true
  printf 'proof-bundle published: %s\n' "$LOG_DIR" >&2 || true
}

exit_for_wrapper_signal() { # numeric-signal
  local signal_number="$1"
  local rc=$((128 + signal_number))
  local proof_scope="focused-software-only"
  trap '' HUP INT TERM
  # The terminal transition has consumed this signal. Keep its number in a
  # diagnostic-only latch, but clear the forwarding latch before any sealing
  # helper is spawned; otherwise a newly armed supervisor can be signalled a
  # second time and make INTERRUPTED publication scheduler-dependent.
  WRAPPER_SIGNAL=0
  CONSUMED_WRAPPER_SIGNAL="$signal_number"
  # This no-Cargo self-test marker is a one-shot readiness witness. Replace its
  # already-created path with a fresh finalizer-handler witness before clearing
  # the pending supervisor binding. The first finalizer consumes the fresh
  # shell-local binding and completes a ready/release/acknowledgement handshake.
  if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_ARM_CONSUMED_SIGNAL_FINALIZER:-0}" == 1 \
      && -n "$SELF_TEST_ACTIVE_SUPERVISOR_MARKER" ]]; then
    SELF_TEST_FINALIZER_READY_MARKER="${SELF_TEST_ACTIVE_SUPERVISOR_MARKER}.finalizer-${signal_number}"
  fi
  if [[ "$PROFILE" == "closure" ]]; then
    proof_scope="head-bookended-closure-candidate"
  fi
  if [[ "$SEALED" == 0 ]]; then
    printf 'wrapper-signal=%s\nno-later-lanes-launched=true\n' "$signal_number" \
      >"$LOG_DIR/logs/wrapper-signal.log"
    printf 'wrapper-signal=%s\nno-later-lanes-launched=true\n' "$signal_number" \
      >&2
    record \
      "wrapper-signal" FAIL "bounded-interrupt-cleanup" \
      "wrapper received signal $signal_number and stopped scheduling" \
      "logs/wrapper-signal.log"
    write_summary_and_seal \
      INTERRUPTED not-checked false "$proof_scope" incomplete "$rc"
    printf 'Euler-disc harness interrupted; retained evidence=%s\n' \
      "$LOG_DIR" >&2
  fi
  exit "$rc"
}

# shellcheck disable=SC2329 # Invoked indirectly by the signal traps below.
forward_wrapper_signal() { # numeric-signal signal-name
  local signal_number="$1" signal_name="$2"
  if [[ "$WRAPPER_SIGNAL" == 0 ]]; then
    WRAPPER_SIGNAL="$signal_number"
  fi
  if [[ "$RETRACTION_INTEGRITY_RECOVERY" == 1 ]]; then
    return
  fi
  if [[ -n "$ACTIVE_SUPERVISOR_PID" ]]; then
    # Signal the active Bash job rather than probing a possibly reused PID.
    kill -s "$signal_name" %+ 2>/dev/null || true
    return
  fi
  if [[ "$SEALING" == 1 || "$LANE_COMMIT_CRITICAL" == 1 ]]; then
    return
  fi
  exit_for_wrapper_signal "$WRAPPER_SIGNAL"
}

# Invoked indirectly by the EXIT trap installed immediately below.
# shellcheck disable=SC2329
seal_on_exit() {
  local rc=$?
  local proof_scope="focused-software-only"
  trap - EXIT
  trap '' HUP INT TERM
  if [[ "$RETRACTION_INTEGRITY_RECOVERY" == 1 ]]; then
    if [[ "$WRAPPER_SIGNAL" != 0 && -n "$RETRACTION_LAST_LOG" ]]; then
      printf 'retraction-integrity-concurrent-wrapper-signal=%s context=%s\n' \
        "$WRAPPER_SIGNAL" "$RETRACTION_INTEGRITY_CONTEXT" \
        >>"$RETRACTION_LAST_LOG" || true
    fi
    WRAPPER_SIGNAL=0
  fi
  if [[ "$PROFILE" == "closure" ]]; then
    proof_scope="head-bookended-closure-candidate"
  fi
  if [[ "$SEALED" == 0 ]]; then
    if [[ "$CANDIDATE_REJECTED" == 1 ]]; then
      start_fresh_incomplete_candidate
    fi
    if [[ "$rc" == 0 ]]; then
      rc=7
    fi
    printf 'unexpected-exit-code=%s\nterminal-seal-was-missing=true\n' "$rc" \
      >"$LOG_DIR/logs/internal-incomplete.log"
    record \
      "internal-incomplete" FAIL "harness-integrity" \
      "EXIT trap sealed an otherwise incomplete run" \
      "logs/internal-incomplete.log"
    if write_summary_and_seal \
      INCOMPLETE not-checked false "$proof_scope" incomplete "$rc"; then
      :
    elif [[ "$CANDIDATE_REJECTED" == 1 ]]; then
      start_fresh_incomplete_candidate
      printf 'unexpected-exit-code=%s\nterminal-seal-was-missing=true\n' "$rc" \
        >"$LOG_DIR/logs/internal-incomplete.log"
      record \
        "internal-incomplete" FAIL "harness-integrity" \
        "EXIT trap sealed an otherwise incomplete run" \
        "logs/internal-incomplete.log"
      write_summary_and_seal \
        INCOMPLETE not-checked false "$proof_scope" incomplete "$rc" || true
    fi
    printf 'Euler-disc harness incomplete; retained evidence=%s\n' \
      "$LOG_DIR" >&2
  fi
  exit "$rc"
}

trap seal_on_exit EXIT
trap 'forward_wrapper_signal 1 HUP' HUP
trap 'forward_wrapper_signal 2 INT' INT
trap 'forward_wrapper_signal 15 TERM' TERM

self_test_capture() { # label timeout cap command...
  local label="$1" timeout_seconds="$2" cap_bytes="$3"
  shift 3
  LANE_COMMIT_CRITICAL=1
  SELF_TEST_LAST_LOG_REL="logs/self-test-${label}.log"
  SELF_TEST_LAST_LOG="$LOG_DIR/$SELF_TEST_LAST_LOG_REL"
  {
    printf 'self-test-label=%s\n' "$label"
    printf 'argv='
    printf '%q ' "$@"
    printf '\n'
    printf '%s\n' '--- supervised execution follows ---'
  } >"$SELF_TEST_LAST_LOG"
  SELF_TEST_LAST_OUTPUT_OFFSET="$(wc -c <"$SELF_TEST_LAST_LOG" | tr -d ' ')"
  if run_bounded_command \
    "$SELF_TEST_LAST_LOG" "$timeout_seconds" "$RUN_DEADLINE_MONOTONIC_NS" \
    "$cap_bytes" "$cap_bytes" "$HARD_MAX_RETAINED_LOG_BYTES" none "$@"; then
    SELF_TEST_LAST_RC=0
  else
    SELF_TEST_LAST_RC=$?
  fi
  if [[ "$WRAPPER_SIGNAL" != 0 ]]; then
    local interrupted_aux_root interrupted_aux_log
    interrupted_aux_root="$(mktemp -d "$LOG_ROOT/.interrupted-command-XXXXXXXX")"
    interrupted_aux_log="$interrupted_aux_root/${SELF_TEST_LAST_LOG_REL##*/}"
    mv "$SELF_TEST_LAST_LOG" "$interrupted_aux_log"
    printf 'interrupted command log retained outside proof bundle: %s\n' \
      "$interrupted_aux_log" >&2
    LANE_COMMIT_CRITICAL=0
    exit_for_wrapper_signal "$WRAPPER_SIGNAL"
  fi
}

self_test_supervised_log_content() {
  local size
  if [[ ! "$SELF_TEST_LAST_OUTPUT_OFFSET" =~ ^(0|[1-9][0-9]*)$ ]]; then
    return 1
  fi
  size="$(wc -c <"$SELF_TEST_LAST_LOG" | tr -d ' ')"
  if [[ ! "$size" =~ ^(0|[1-9][0-9]*)$ ]] \
      || ((10#$size < 10#$SELF_TEST_LAST_OUTPUT_OFFSET)); then
    return 1
  fi
  tail -c "+$((10#$SELF_TEST_LAST_OUTPUT_OFFSET + 1))" "$SELF_TEST_LAST_LOG"
}

self_test_append_producer_disposition() {
  local disposition detail disposition_abort supervisor_reason
  if disposition="$(
    supervisor_disposition_for_log "$SELF_TEST_LAST_LOG" "$SELF_TEST_LAST_RC"
  )"; then
    IFS=$'\t' read -r supervisor_reason detail disposition_abort <<<"$disposition"
    {
      printf 'producer-shutdown-reason=%s\n' "$supervisor_reason"
      printf 'producer-detail=%s\n' "$detail"
      printf 'producer-abort=%s\n' "$disposition_abort"
    } >>"$SELF_TEST_LAST_LOG"
  else
    printf '%s\n' 'producer-disposition-classification=failed' \
      >>"$SELF_TEST_LAST_LOG"
    SELF_TEST_LAST_RC=1
  fi
}

self_test_assert() { # label expected-rc detail marker...
  local label="$1" expected_rc="$2" detail="$3" marker content="" passed=1
  shift 3
  if ! content="$(self_test_supervised_log_content)"; then
    passed=0
    printf '%s\n' 'supervised-output-boundary-invalid=true' \
      >>"$SELF_TEST_LAST_LOG"
  fi
  if [[ "$SELF_TEST_LAST_RC" != "$expected_rc" ]]; then
    passed=0
    printf 'expected-rc=%s observed-rc=%s\n' "$expected_rc" "$SELF_TEST_LAST_RC" \
      >>"$SELF_TEST_LAST_LOG"
  fi
  for marker in "$@"; do
    if [[ "$marker" == absent:* ]]; then
      if [[ "$content" == *"${marker#absent:}"* ]]; then
        passed=0
        printf 'unexpected-marker=%s\n' "${marker#absent:}" >>"$SELF_TEST_LAST_LOG"
      fi
    elif [[ "$marker" == regex:* ]]; then
      if [[ ! "$content" =~ ${marker#regex:} ]]; then
        passed=0
        printf 'missing-regex=%s\n' "${marker#regex:}" >>"$SELF_TEST_LAST_LOG"
      fi
    elif [[ "$content" != *"$marker"* ]]; then
      passed=0
      printf 'missing-marker=%s\n' "$marker" >>"$SELF_TEST_LAST_LOG"
    fi
  done
  {
    printf 'self-test-expected-rc=%s\n' "$expected_rc"
    printf 'self-test-assertion-detail=%s\n' "$detail"
    if [[ "$passed" == 1 ]]; then
      printf '%s\n' 'self-test-assertion-result=pass'
    else
      printf '%s\n' 'self-test-assertion-result=fail'
    fi
  } >>"$SELF_TEST_LAST_LOG"
  if [[ "$passed" == 1 ]]; then
    record "$label" PASS "harness-self-test" \
      "harness self-test assertion matched expected disposition" \
      "$SELF_TEST_LAST_LOG_REL"
  else
    record "$label" FAIL "harness-self-test" \
      "harness self-test assertion did not match expected disposition" \
      "$SELF_TEST_LAST_LOG_REL"
  fi
  LANE_COMMIT_CRITICAL=0
  if [[ "$WRAPPER_SIGNAL" != 0 ]]; then
    exit_for_wrapper_signal "$WRAPPER_SIGNAL"
  fi
}

self_test_assert_external_session_empty() {
  local sid
  if ! sid="$(
    self_test_supervised_log_content | sed -n \
      -e 's/.*supervisor_sid=\([0-9][0-9]*\).*/\1/p' \
      -e 's/.*emergency-supervisor-cleanup .*sid=\([0-9][0-9]*\).*/\1/p' \
      | tail -n 1
  )"; then
    sid=""
  fi
  if [[ -z "$sid" ]]; then
    printf '%s\n' 'external_process_session_drained=false missing_sid=true' \
      >>"$SELF_TEST_LAST_LOG"
    SELF_TEST_LAST_RC=1
    return
  fi
  if python3 - "$sid" >>"$SELF_TEST_LAST_LOG" <<'PY'
import os
import subprocess
import sys

sid = int(sys.argv[1])
environment = os.environ.copy()
environment["LC_ALL"] = "C"
result = subprocess.run(
    ["/bin/ps", "-axo", "pid=,pgid=,stat="],
    check=False,
    capture_output=True,
    text=True,
    timeout=2,
    env=environment,
)
if result.returncode != 0:
    raise SystemExit(f"external ps failed with {result.returncode}: {result.stderr.strip()}")
live = []
for line in result.stdout.splitlines():
    fields = line.split(None, 2)
    if len(fields) != 3:
        continue
    pid_text, _pgid_text, state = fields
    try:
        pid = int(pid_text)
        observed_sid = os.getsid(pid)
    except (ValueError, ProcessLookupError):
        continue
    if observed_sid == sid and not state.startswith("Z"):
        live.append(pid)
print(f"external_process_session_drained={str(not live).lower()} live={live!r}")
raise SystemExit(0 if not live else 1)
PY
  then
    :
  else
    SELF_TEST_LAST_RC=1
  fi
}

self_test_verify_nested_bundle() { # retained-evidence-line-prefix expected-log-root [allow-other-roots]
  local prefix="$1" expected_root="${2%/}" allow_other_roots="${3:-0}"
  local original_rc="$SELF_TEST_LAST_RC" nested_dir="" line candidate relative
  local locator_count=0 matching_locator_count=0
  while IFS= read -r line; do
    case "$line" in
      "$prefix"*)
        locator_count=$((locator_count + 1))
        candidate="${line#"$prefix"}"
        case "$candidate" in
          "$expected_root"/*)
            relative="${candidate#"$expected_root"/}"
            if [[ -n "$relative" && "$relative" != .* && "$relative" != */* ]]; then
              matching_locator_count=$((matching_locator_count + 1))
              nested_dir="$candidate"
            fi
            ;;
        esac
        ;;
    esac
  done < <(self_test_supervised_log_content)
  if [[ "$matching_locator_count" != 1 ]] \
      || { [[ "$allow_other_roots" != 1 ]] && [[ "$locator_count" != 1 ]]; }; then
    printf 'nested-bundle-locator-invalid prefix=%s expected-root=%s locators=%s matching-direct-children=%s allow-other-roots=%s\n' \
      "$prefix" "$expected_root" "$locator_count" "$matching_locator_count" \
      "$allow_other_roots" \
      >>"$SELF_TEST_LAST_LOG"
    SELF_TEST_LAST_RC=1
    return
  fi
  printf 'nested-bundle-locator=%s\n' "$nested_dir" >>"$SELF_TEST_LAST_LOG"
  if run_bounded_command \
    "$SELF_TEST_LAST_LOG" 5 "$RUN_DEADLINE_MONOTONIC_NS" \
    1048576 1048576 "$HARD_MAX_RETAINED_LOG_BYTES" none \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$nested_dir"; then
    SELF_TEST_LAST_RC="$original_rc"
  else
    printf '%s\n' 'nested proof-bundle verification failed' \
      >>"$SELF_TEST_LAST_LOG"
    SELF_TEST_LAST_RC=1
  fi
}

create_verifier_self_test_fixtures() {
  python3 - "$PROOF_MAX_BUNDLE_ENTRIES" "$@" <<'PY'
import hashlib
import json
import pathlib
import re
import shutil
import sys

max_bundle_entries = int(sys.argv[1])
(
    valid,
    mutated,
    duplicate,
    truncated,
    summary_mismatch,
    unsafe_path,
    duplicate_key,
    nonfinite,
    oversized,
    entry_cap,
    too_many,
    valid_success,
    readiness_mismatch,
    no_claim_mutation,
    unknown_terminal,
    unexpected_file,
    prefix_mismatch,
    invalid_interrupted,
    invalid_incomplete,
    wrong_authority,
    wrong_command,
    control_continuation,
    premature_snapshot,
    toctou_mutation,
    proof_body_mutation,
    focused_source_body_mutation,
    closure_source_body_mutation,
    supervisor_contradiction,
    supervisor_containment_contradiction,
    supervisor_flag_contradiction,
    supervisor_log_cap_contradiction,
    supervisor_retained_cap_contradiction,
    supervisor_deadline_contradiction,
    valid_closure,
    closure_missing_lane,
    closure_snapshot_mismatch,
    valid_self_test,
    self_test_wrong_authority,
    self_test_wrong_detail,
    self_test_wrong_log,
) = map(pathlib.Path, sys.argv[2:])

NO_CLAIM = "Software/protocol structure only; no physical validation or emergent Euler-disc prediction."
CONCURRENCY_NO_CLAIM = "Clean HEAD is observed at bookends; transient concurrent edits between observations are not excluded."
PROOF_LOG = (
    "software/protocol structure only\n"
    "no physical validation\n"
    "no mechanism attribution\n"
    "no emergent Euler-disc prediction\n"
    "retained_log_checker_smoke_command=cargo test --locked -p fs-euler-disc-e2e "
    "--test scientific_contract -- "
    "g0_check_receipts_logs_and_domain_separation_are_exact_and_bounded --exact "
    "--test-threads=1\n"
    "retained command is not packet/case replay and resolves no artifacts\n"
).encode()
FOCUSED_MANIFEST_LOG = (
    b"NO_DATA: focused profile does not cover untracked or dirty candidate paths\n"
)
HEAD = "self-test-head"
ISA = "self-test-isa"
EXECUTOR = "local"
SCOPE_PATHS = (
    "Cargo.toml", "Cargo.lock", "README.md", "consolidation-review.json",
    "docs/CONSOLIDATION_REVIEW.md", "docs/CONVENTIONS.md", "docs/SCHEMA_POLICY.md",
    "doc-facts-inventory.json", "schema-policy.json", "identity-authorities.json",
    "identity-schemas.json", "golden-couplings.json", "xtask/src/identities.rs",
    "scripts/ci/euler_disc_contract_e2e.sh", "crates/fs-euler-disc-e2e/Cargo.toml",
    "crates/fs-euler-disc-e2e/CONTRACT.md", "crates/fs-euler-disc-e2e/src/lib.rs",
    "crates/fs-euler-disc-e2e/src/contract.rs", "crates/fs-euler-disc-e2e/src/protocol.rs",
    "crates/fs-euler-disc-e2e/tests/scientific_contract.rs",
)
AUTHORITIES = {
    "proof-boundary": "declaration-only",
    "source-manifest": "manifest-not-evaluated",
    "closure-root-preflight": "full-root-clean-head-observation",
    "constellation-verify": "constellation-source-preflight",
    "constellation-snapshot-before": "source-provenance-snapshot",
    "crate-fmt": "static-hygiene-only",
    "crate-check": "focused-software-evidence",
    "retained-log-checker-smoke": "focused-software-evidence",
    "retained-log-checker-smoke-sentinel": "non-vacuity-evidence",
    "crate-unit-integration": "focused-software-evidence",
    "crate-doctest-hostile-boundary": "compile-fail-api-evidence",
    "crate-clippy": "static-hygiene-only",
    "xtask-check-layers": "workspace-structural-gate",
    "xtask-check-deps": "workspace-structural-gate",
    "xtask-check-contracts": "workspace-structural-gate",
    "xtask-check-schemas": "workspace-structural-gate",
    "xtask-check-consolidation": "workspace-structural-gate",
    "xtask-check-identities": "workspace-structural-gate",
    "xtask-check-goldens": "workspace-structural-gate",
    "xtask-check-docs": "workspace-documentation-gate",
    "xtask-check-source-manifest": "head-bound-source-inventory",
    "source-manifest-membership": "independent-path-membership",
    "constellation-snapshot-after": "source-provenance-snapshot",
    "closure-root-bookend": "full-root-clean-head-observation",
    "source-stability": "source-provenance-snapshot",
}
CARGO_TAILS = {
    "crate-fmt": ["fmt", "-p", "fs-euler-disc-e2e", "--check"],
    "crate-check": ["check", "--locked", "-p", "fs-euler-disc-e2e", "--all-targets"],
    "retained-log-checker-smoke": [
        "test", "--locked", "-p", "fs-euler-disc-e2e", "--test",
        "scientific_contract", "--",
        "g0_check_receipts_logs_and_domain_separation_are_exact_and_bounded",
        "--exact", "--test-threads=1",
    ],
    "crate-unit-integration": [
        "test", "--locked", "--no-fail-fast", "-p", "fs-euler-disc-e2e",
        "--lib", "--tests", "--", "--test-threads=1",
    ],
    "crate-doctest-hostile-boundary": [
        "test", "--locked", "-p", "fs-euler-disc-e2e", "--doc",
    ],
    "crate-clippy": [
        "clippy", "--locked", "-p", "fs-euler-disc-e2e", "--all-targets",
        "--no-deps", "--", "-D", "warnings",
    ],
}
for gate in (
    "check-layers", "check-deps", "check-contracts", "check-schemas",
    "check-consolidation", "check-identities", "check-goldens", "check-docs",
    "check-source-manifest",
):
    CARGO_TAILS[f"xtask-{gate}"] = [
        "run", "--locked", "-q", "-p", "xtask", "--", gate,
    ]

def canonical(value):
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"
    ).encode()

def command_argv(lane):
    if lane in CARGO_TAILS:
        return ["/synthetic/cargo", *CARGO_TAILS[lane]]
    if lane == "constellation-verify":
        return ["scripts/ci/checkout_constellation.sh", "--verify-only"]
    if lane in {"closure-root-preflight", "closure-root-bookend"}:
        return [
            "bash", "-c", 'closure_root_preflight_command "$@"', "_",
            HEAD, *SCOPE_PATHS,
        ]
    if lane in {"constellation-snapshot-before", "constellation-snapshot-after"}:
        suffix = "before" if lane.endswith("before") else "after"
        return [
            "bash", "-c", 'capture_snapshot_command "$@"', "_",
            f"/synthetic/snapshot-{suffix}.txt",
            f"/synthetic/.{lane}-candidate-ABCDEFGH",
            HEAD,
            *SCOPE_PATHS,
        ]
    if lane == "retained-log-checker-smoke-sentinel":
        return [
            "bash", "-c", 'checker_smoke_sentinel_command "$@"', "_",
            "/synthetic/logs/retained-log-checker-smoke.log",
            "g0_check_receipts_logs_and_domain_separation_are_exact_and_bounded",
        ]
    if lane == "source-stability":
        return [
            "cmp", "-s", "/synthetic/snapshot-before.txt",
            "/synthetic/snapshot-after.txt",
        ]
    if lane == "source-manifest-membership":
        return [
            "bash", "-c", 'source_manifest_membership_command "$@"', "_",
            "/synthetic/frankensim-source-manifest.json", *SCOPE_PATHS,
        ]
    raise ValueError(f"no synthetic command for {lane}")

def command_log(
    lane,
    argv=None,
    *,
    profile="focused",
    executor=EXECUTOR,
    rc=0,
):
    if argv is None:
        argv = command_argv(lane)
    header_fields = [
        "schema=frankensim.euler-disc-contract-e2e.command.v1",
        f"lane={lane}",
        f"head={HEAD}",
        f"host_isa={ISA}",
        f"profile={profile}",
        f"executor_declaration={executor}",
        "executor_attestation=caller-declared-unverified",
        "lane_timeout_seconds=5",
        "run_timeout_seconds=10",
        "run_started_monotonic_ns=1000000000",
        "run_deadline_monotonic_ns=11000000000",
        "lane_log_max_bytes=1048576",
        f"authority={AUTHORITIES[lane]}",
        "argv_json=" + json.dumps(argv, ensure_ascii=True, separators=(",", ":")),
        "--- command output ---",
    ]
    header = ("\n".join(header_fields) + "\n").encode()
    if rc == 0:
        reason = "none"
        leader_exit_code = 0
    else:
        reason = "leader-exit"
        leader_exit_code = rc
    result = {
        "configured_lane_timeout_seconds": 5,
        "configured_output_cap_bytes": 1048576,
        "deadline_kind": "lane",
        "effective_deadline_monotonic_ns": 7000000000,
        "initial_log_bytes": len(header),
        "inspection_error_count": 0,
        "interrupted_signal": None,
        "kill_sent": False,
        "lane_deadline_monotonic_ns": 7000000000,
        "leader_exit_code": leader_exit_code,
        "metadata_complete": True,
        "output_pipe_eof": True,
        "output_truncated": False,
        "process_group_drained": True,
        "process_session_drained": True,
        "python_cap_reduced": False,
        "retained_output_cap_bytes": 1048576,
        "run_deadline_monotonic_ns": 11000000000,
        "schema": "frankensim.euler-disc-contract-e2e.supervisor-result.v1",
        "shell_cap_reduced": False,
        "shell_effective_output_cap_bytes": 1048576,
        "shutdown_reason": reason,
        "supervisor_started_monotonic_ns": 2000000000,
        "term_sent": False,
        "total_log_cap_bytes": 67108864,
        "wrapper_exit_code": rc,
    }
    result_line = b"supervisor_result_json=" + canonical(result)
    return header + f"synthetic-command-lane={lane}\n".encode() + result_line

def make_record(
    destination,
    lane,
    status,
    detail=None,
    *,
    data=None,
    authority=None,
    profile="focused",
    executor=EXECUTOR,
    log_name=None,
):
    relative = f"logs/{log_name or lane}.log"
    if detail is None:
        detail = "command completed" if status == "PASS" else "command exited 2"
    if data is None:
        if lane == "proof-boundary":
            data = PROOF_LOG
        elif lane == "source-manifest":
            data = FOCUSED_MANIFEST_LOG
        else:
            match = re.fullmatch(r"command exited ([0-9]+)", detail)
            rc = 0 if status == "PASS" else int(match.group(1)) if match else 2
            data = command_log(lane, profile=profile, executor=executor, rc=rc)
    (destination / relative).write_bytes(data)
    return {
        "schema": "frankensim.euler-disc-contract-e2e.verdict.v1",
        "lane": lane,
        "status": status,
        "authority": authority or AUTHORITIES[lane],
        "detail": detail,
        "log": relative,
        "log_bytes": len(data),
        "log_sha256": hashlib.sha256(data).hexdigest(),
        "head": HEAD,
        "host_isa": ISA,
        "profile": profile,
        "executor_declaration": executor,
        "executor_attestation": "caller-declared-unverified",
        "provenance_state": "provisional",
        "terminal": False,
    }

def base_seal(
    status,
    checks,
    failures,
    exit_code,
    *,
    provenance="incomplete",
    profile="focused",
    executor=EXECUTOR,
    coverage="not-checked",
    ready=False,
    proof_scope="focused-software-only",
):
    return {
        "schema": "frankensim.euler-disc-contract-e2e.proof-seal.v2",
        "record_type": "proof-seal",
        "status": status,
        "checks": checks,
        "failures": failures,
        "head": HEAD,
        "host_isa": ISA,
        "profile": profile,
        "executor_declaration": executor,
        "executor_attestation": "caller-declared-unverified",
        "source_manifest_coverage": coverage,
        "candidate_ready_for_dsr": ready,
        "dsr_status": "not-run-by-this-harness",
        "proof_scope": proof_scope,
        "provenance_state": provenance,
        "snapshot_before_sha256": "",
        "snapshot_after_sha256": "",
        "verdicts_prefix_sha256": "",
        "verdicts": "verdicts.jsonl",
        "proof_seal_locator": "verdicts.jsonl#last-line",
        "terminal_exit_code": exit_code,
        "terminal": True,
        "no_claim": NO_CLAIM,
        "concurrency_no_claim": CONCURRENCY_NO_CLAIM,
    }

def write_bundle(destination, records, seal):
    destination.mkdir(parents=True, exist_ok=True)
    prefix = b"".join(canonical(record) for record in records)
    terminal = dict(seal)
    terminal["verdicts_prefix_sha256"] = hashlib.sha256(prefix).hexdigest()
    seal_line = canonical(terminal)
    (destination / "verdicts-prefix.jsonl").write_bytes(prefix)
    (destination / "verdicts.jsonl").write_bytes(prefix + seal_line)
    (destination / "summary.json").write_bytes(seal_line)
    return terminal

def read_bundle(destination):
    lines = (destination / "verdicts.jsonl").read_bytes().splitlines()
    return [json.loads(line) for line in lines[:-1]], json.loads(lines[-1])

def rewrite_records(destination, records, terminal):
    terminal = dict(terminal)
    terminal["checks"] = len(records)
    terminal["failures"] = sum(record["status"] == "FAIL" for record in records)
    write_bundle(destination, records, terminal)

def rewrite_terminal(destination, **changes):
    lines = (destination / "verdicts.jsonl").read_bytes().splitlines(keepends=True)
    prefix = b"".join(lines[:-1])
    terminal = json.loads(lines[-1])
    terminal.update(changes)
    terminal["verdicts_prefix_sha256"] = hashlib.sha256(prefix).hexdigest()
    seal_line = canonical(terminal)
    (destination / "verdicts.jsonl").write_bytes(prefix + seal_line)
    (destination / "summary.json").write_bytes(seal_line)

def rewrite_lane_log(destination, lane, data):
    records, terminal = read_bundle(destination)
    matches = [record for record in records if record["lane"] == lane]
    if len(matches) != 1:
        raise ValueError(f"fixture lacks exactly one lane {lane}")
    record = matches[0]
    (destination / record["log"]).write_bytes(data)
    record["log_bytes"] = len(data)
    record["log_sha256"] = hashlib.sha256(data).hexdigest()
    rewrite_records(destination, records, terminal)

(valid / "logs").mkdir(parents=True)
proof_lane = make_record(
    valid, "proof-boundary", "PASS",
    "software/protocol structure only; no physical validation or emergent prediction",
)
manifest_lane = make_record(
    valid, "source-manifest", "NO_DATA",
    "focused profile does not claim untracked/dirty Euler-disc paths are covered",
)
failed_lane = make_record(valid, "constellation-verify", "FAIL", "command exited 2")
seal = write_bundle(
    valid,
    [proof_lane, manifest_lane, failed_lane],
    base_seal("FAIL", 3, 1, 7),
)
for destination in (
    mutated,
    duplicate,
    truncated,
    summary_mismatch,
    unsafe_path,
    duplicate_key,
    nonfinite,
    unknown_terminal,
    unexpected_file,
    prefix_mismatch,
    wrong_authority,
    wrong_command,
    premature_snapshot,
):
    shutil.copytree(valid, destination)

with (mutated / "logs" / "constellation-verify.log").open("ab") as handle:
    handle.write(b"tamper\n")

duplicate_prefix = (duplicate / "verdicts.jsonl").read_bytes()
with (duplicate / "verdicts.jsonl").open("ab") as handle:
    handle.write(canonical(seal))
(duplicate / "verdicts-prefix.jsonl").write_bytes(duplicate_prefix)

truncated_bytes = (truncated / "verdicts.jsonl").read_bytes()
(truncated / "verdicts.jsonl").write_bytes(truncated_bytes[:-1])

summary_only = json.loads((summary_mismatch / "summary.json").read_bytes())
summary_only["no_claim"] = "mutated canonical summary only"
(summary_mismatch / "summary.json").write_bytes(canonical(summary_only))

unsafe_records, unsafe_seal = read_bundle(unsafe_path)
unsafe_records[-1]["log"] = "../escape.log"
rewrite_records(unsafe_path, unsafe_records, unsafe_seal)

duplicate_lane = canonical(failed_lane)[:-2] + b',"status":"PASS"}\n'
duplicate_key_seal = dict(seal)
duplicate_key_seal["verdicts_prefix_sha256"] = hashlib.sha256(duplicate_lane).hexdigest()
duplicate_key_seal_line = canonical(duplicate_key_seal)
(duplicate_key / "verdicts.jsonl").write_bytes(duplicate_lane + duplicate_key_seal_line)
(duplicate_key / "summary.json").write_bytes(duplicate_key_seal_line)
(duplicate_key / "verdicts-prefix.jsonl").write_bytes(duplicate_lane)

finite_lane = canonical(failed_lane)
nonfinite_lane = finite_lane.replace(
    f'"log_bytes":{failed_lane["log_bytes"]}'.encode(), b'"log_bytes":NaN', 1
)
nonfinite_seal = dict(seal)
nonfinite_seal["verdicts_prefix_sha256"] = hashlib.sha256(nonfinite_lane).hexdigest()
nonfinite_seal_line = canonical(nonfinite_seal)
(nonfinite / "verdicts.jsonl").write_bytes(nonfinite_lane + nonfinite_seal_line)
(nonfinite / "summary.json").write_bytes(nonfinite_seal_line)
(nonfinite / "verdicts-prefix.jsonl").write_bytes(nonfinite_lane)

oversized.mkdir(parents=True)
(oversized / "verdicts.jsonl").write_bytes(b"x" * (2 * 1024 * 1024 + 1))
(oversized / "summary.json").write_bytes(canonical(seal))
(oversized / "verdicts-prefix.jsonl").write_bytes(b"")

entry_cap.mkdir(parents=True)
for index in range(max_bundle_entries + 1):
    (entry_cap / f"entry-{index:03d}").write_bytes(b"")

too_many_records = []
too_many.mkdir(parents=True)
for index in range(96):
    record = dict(failed_lane)
    record.update({
        "lane": f"synthetic-lane-{index}",
        "log": f"logs/lane-{index}.log",
    })
    too_many_records.append(record)
too_many_seal = dict(seal)
too_many_seal.update({"checks": 96, "failures": 96})
write_bundle(too_many, too_many_records, too_many_seal)

focused_sequence = (
    "proof-boundary", "source-manifest", "constellation-verify",
    "constellation-snapshot-before", "crate-fmt", "crate-check",
    "retained-log-checker-smoke", "retained-log-checker-smoke-sentinel",
    "crate-unit-integration", "crate-doctest-hostile-boundary", "crate-clippy",
    "xtask-check-layers", "xtask-check-deps", "xtask-check-contracts",
    "xtask-check-schemas", "xtask-check-consolidation", "xtask-check-identities",
    "xtask-check-goldens",
    "constellation-snapshot-after", "source-stability",
)
(valid_success / "logs").mkdir(parents=True)
success_records = []
for lane_name in focused_sequence:
    if lane_name == "proof-boundary":
        success_records.append(make_record(
            valid_success, lane_name, "PASS",
            "software/protocol structure only; no physical validation or emergent prediction",
        ))
    elif lane_name == "source-manifest":
        success_records.append(make_record(
            valid_success, lane_name, "NO_DATA",
            "focused profile does not claim untracked/dirty Euler-disc paths are covered",
        ))
    else:
        success_records.append(make_record(valid_success, lane_name, "PASS"))
snapshot = b"synthetic equal source snapshot\n"
(valid_success / "snapshot-before.txt").write_bytes(snapshot)
(valid_success / "snapshot-after.txt").write_bytes(snapshot)
success_seal = base_seal(
    "FOCUSED_PASS", len(success_records), 0, 0, provenance="stable",
)
success_seal.update({
    "snapshot_before_sha256": hashlib.sha256(snapshot).hexdigest(),
    "snapshot_after_sha256": hashlib.sha256(snapshot).hexdigest(),
})
write_bundle(valid_success, success_records, success_seal)

shutil.copytree(valid_success, proof_body_mutation)
rewrite_lane_log(
    proof_body_mutation,
    "proof-boundary",
    PROOF_LOG.replace(b"no mechanism attribution\n", b"mechanism attribution allowed\n"),
)
shutil.copytree(valid_success, focused_source_body_mutation)
rewrite_lane_log(
    focused_source_body_mutation,
    "source-manifest",
    b"NO_DATA: focused source-manifest body was consistently rehashed\n",
)

shutil.copytree(valid, supervisor_contradiction)
supervisor_records, supervisor_terminal = read_bundle(supervisor_contradiction)
supervisor_record = supervisor_records[-1]
supervisor_log_path = supervisor_contradiction / supervisor_record["log"]
supervisor_lines = supervisor_log_path.read_bytes().splitlines(keepends=True)
supervisor_prefix = b"supervisor_result_json="
supervisor_result = json.loads(supervisor_lines[-1][len(supervisor_prefix):])
supervisor_result.update({
    "leader_exit_code": 0,
    "shutdown_reason": "none",
    "wrapper_exit_code": 0,
})
supervisor_lines[-1] = supervisor_prefix + canonical(supervisor_result)
supervisor_log = b"".join(supervisor_lines)
supervisor_log_path.write_bytes(supervisor_log)
supervisor_record["log_bytes"] = len(supervisor_log)
supervisor_record["log_sha256"] = hashlib.sha256(supervisor_log).hexdigest()
rewrite_records(supervisor_contradiction, supervisor_records, supervisor_terminal)

shutil.copytree(valid, supervisor_containment_contradiction)
containment_records, containment_terminal = read_bundle(
    supervisor_containment_contradiction
)
containment_record = containment_records[-1]
containment_log_path = supervisor_containment_contradiction / containment_record["log"]
containment_lines = containment_log_path.read_bytes().splitlines(keepends=True)
containment_result = json.loads(containment_lines[-1][len(supervisor_prefix):])
containment_result.update({
    "output_pipe_eof": False,
    "process_group_drained": False,
    "process_session_drained": False,
})
containment_lines[-1] = supervisor_prefix + canonical(containment_result)
containment_log = b"".join(containment_lines)
containment_log_path.write_bytes(containment_log)
containment_record["log_bytes"] = len(containment_log)
containment_record["log_sha256"] = hashlib.sha256(containment_log).hexdigest()
rewrite_records(
    supervisor_containment_contradiction,
    containment_records,
    containment_terminal,
)

shutil.copytree(valid, supervisor_flag_contradiction)
flag_records, flag_terminal = read_bundle(supervisor_flag_contradiction)
flag_record = flag_records[-1]
flag_log_path = supervisor_flag_contradiction / flag_record["log"]
flag_lines = flag_log_path.read_bytes().splitlines(keepends=True)
flag_result = json.loads(flag_lines[-1][len(supervisor_prefix):])
flag_result["output_truncated"] = True
flag_lines[-1] = supervisor_prefix + canonical(flag_result)
flag_log = b"".join(flag_lines)
flag_log_path.write_bytes(flag_log)
flag_record["log_bytes"] = len(flag_log)
flag_record["log_sha256"] = hashlib.sha256(flag_log).hexdigest()
rewrite_records(supervisor_flag_contradiction, flag_records, flag_terminal)

shutil.copytree(valid, supervisor_log_cap_contradiction)
log_cap_records, log_cap_terminal = read_bundle(supervisor_log_cap_contradiction)
log_cap_record = log_cap_records[-1]
log_cap_path = supervisor_log_cap_contradiction / log_cap_record["log"]
log_cap_lines = log_cap_path.read_bytes().splitlines(keepends=True)
log_cap_result = json.loads(log_cap_lines[-1][len(supervisor_prefix):])
log_cap_result.update({
    "python_cap_reduced": True,
    "retained_output_cap_bytes": 0,
    "total_log_cap_bytes": log_cap_result["initial_log_bytes"],
})
log_cap_lines[-1] = supervisor_prefix + canonical(log_cap_result)
log_cap_data = b"".join(log_cap_lines)
log_cap_path.write_bytes(log_cap_data)
log_cap_record["log_bytes"] = len(log_cap_data)
log_cap_record["log_sha256"] = hashlib.sha256(log_cap_data).hexdigest()
rewrite_records(
    supervisor_log_cap_contradiction,
    log_cap_records,
    log_cap_terminal,
)

shutil.copytree(valid, supervisor_retained_cap_contradiction)
retained_cap_records, retained_cap_terminal = read_bundle(
    supervisor_retained_cap_contradiction
)
retained_cap_record = retained_cap_records[-1]
retained_cap_path = supervisor_retained_cap_contradiction / retained_cap_record["log"]
retained_cap_lines = retained_cap_path.read_bytes().splitlines(keepends=True)
retained_cap_result = json.loads(retained_cap_lines[-1][len(supervisor_prefix):])
retained_cap_result.update({
    "python_cap_reduced": True,
    "retained_output_cap_bytes": 0,
})
retained_cap_lines[-1] = supervisor_prefix + canonical(retained_cap_result)
retained_cap_data = b"".join(retained_cap_lines)
retained_cap_path.write_bytes(retained_cap_data)
retained_cap_record["log_bytes"] = len(retained_cap_data)
retained_cap_record["log_sha256"] = hashlib.sha256(retained_cap_data).hexdigest()
rewrite_records(
    supervisor_retained_cap_contradiction,
    retained_cap_records,
    retained_cap_terminal,
)

shutil.copytree(valid, supervisor_deadline_contradiction)
deadline_records, deadline_terminal = read_bundle(supervisor_deadline_contradiction)
deadline_record = deadline_records[-1]
deadline_log_path = supervisor_deadline_contradiction / deadline_record["log"]
deadline_lines = deadline_log_path.read_bytes().splitlines(keepends=True)
deadline_result = json.loads(deadline_lines[-1][len(supervisor_prefix):])
deadline_result.update({
    "deadline_kind": "aggregate",
    "effective_deadline_monotonic_ns": deadline_result["run_deadline_monotonic_ns"],
    "lane_deadline_monotonic_ns": (
        deadline_result["run_deadline_monotonic_ns"]
        + deadline_result["configured_lane_timeout_seconds"] * 1_000_000_000
    ),
    "supervisor_started_monotonic_ns": deadline_result["run_deadline_monotonic_ns"],
})
deadline_lines[-1] = supervisor_prefix + canonical(deadline_result)
deadline_log = b"".join(deadline_lines)
deadline_log_path.write_bytes(deadline_log)
deadline_record["log_bytes"] = len(deadline_log)
deadline_record["log_sha256"] = hashlib.sha256(deadline_log).hexdigest()
rewrite_records(
    supervisor_deadline_contradiction,
    deadline_records,
    deadline_terminal,
)

closure_sequence = (
    "proof-boundary", "closure-root-preflight", "constellation-verify",
    "constellation-snapshot-before", "crate-fmt", "crate-check",
    "retained-log-checker-smoke", "retained-log-checker-smoke-sentinel",
    "crate-unit-integration", "crate-doctest-hostile-boundary", "crate-clippy",
    "xtask-check-layers", "xtask-check-deps", "xtask-check-contracts",
    "xtask-check-schemas", "xtask-check-consolidation", "xtask-check-identities",
    "xtask-check-goldens",
    "xtask-check-docs", "xtask-check-source-manifest",
    "source-manifest-membership", "constellation-snapshot-after",
    "closure-root-bookend", "source-stability",
)

def build_closure_bundle(destination, sequence, before, after):
    (destination / "logs").mkdir(parents=True)
    records = []
    for lane_name in sequence:
        if lane_name == "proof-boundary":
            records.append(make_record(
                destination,
                lane_name,
                "PASS",
                "software/protocol structure only; no physical validation or emergent prediction",
                profile="closure",
            ))
        else:
            records.append(make_record(
                destination, lane_name, "PASS", profile="closure",
            ))
    (destination / "snapshot-before.txt").write_bytes(before)
    (destination / "snapshot-after.txt").write_bytes(after)
    closure_seal = base_seal(
        "READY_FOR_DSR",
        len(records),
        0,
        0,
        provenance="stable",
        profile="closure",
        coverage="full-root-clean-head-bookended",
        ready=True,
        proof_scope="head-bookended-closure-candidate",
    )
    closure_seal.update({
        "snapshot_before_sha256": hashlib.sha256(before).hexdigest(),
        "snapshot_after_sha256": hashlib.sha256(after).hexdigest(),
    })
    write_bundle(destination, records, closure_seal)

closure_snapshot = b"synthetic equal closure snapshot\n"
build_closure_bundle(
    valid_closure, closure_sequence, closure_snapshot, closure_snapshot,
)
build_closure_bundle(
    closure_missing_lane,
    tuple(lane for lane in closure_sequence if lane != "xtask-check-docs"),
    closure_snapshot,
    closure_snapshot,
)
build_closure_bundle(
    closure_snapshot_mismatch,
    closure_sequence,
    closure_snapshot,
    b"synthetic moved closure snapshot\n",
)

(closure_source_body_mutation / "logs").mkdir(parents=True)
closure_refusal_records = [
    make_record(
        closure_source_body_mutation,
        "proof-boundary",
        "PASS",
        "software/protocol structure only; no physical validation or emergent prediction",
        profile="closure",
    ),
    make_record(
        closure_source_body_mutation,
        "closure-root-preflight",
        "NO_DATA",
        "command exited 1",
        profile="closure",
    ),
    make_record(
        closure_source_body_mutation,
        "source-manifest",
        "NO_DATA",
        "docs/source-manifest gates skipped because closure preconditions failed",
        data=b"NO_DATA: closure source-manifest body was consistently rehashed\n",
        authority="closure-refused-before-manifest-evaluation",
        profile="closure",
    ),
]
closure_refusal_seal = base_seal(
    "NO_DATA",
    len(closure_refusal_records),
    0,
    4,
    provenance="incomplete",
    profile="closure",
    coverage="not-checked",
    proof_scope="head-bookended-closure-candidate",
)
write_bundle(
    closure_source_body_mutation,
    closure_refusal_records,
    closure_refusal_seal,
)

self_test_sequence = (
    "self-test-ordinary-pass", "self-test-exact-nonzero",
    "self-test-reserved-child-exit", "self-test-child-signal-exit",
    "self-test-spawn-failure", "self-test-timeout-stubborn-group",
    "self-test-leader-exit-live-group", "self-test-setpgid-escape-drained",
    "self-test-output-truncation-drained", "self-test-cap-plus-timeout",
    "self-test-lane-log-cap-boundary", "self-test-numeric-config-boundaries",
    "self-test-completion-classification-signal",
    "self-test-snapshot-byte-boundaries", "self-test-snapshot-hash-failure",
    "self-test-timeout-hanging-helper",
    "self-test-supervisor-exception", "self-test-special-index-flag-refusal",
    "self-test-skip-worktree-flag-refusal",
    "self-test-fsmonitor-valid-flag-refusal",
    "self-test-real-repository-read-failure",
    "self-test-zero-smoke-refusal",
    "self-test-invalid-consolidation-disposition-refusal",
    "self-test-consolidation-scope-mutation-refusal",
    "self-test-success-deadline-publication-refusal",
    "self-test-postpublication-deadline-retraction",
    "self-test-valid-bundle",
    "self-test-valid-success-bundle", "self-test-valid-closure-ready-bundle",
    "self-test-closure-ready-missing-lane-refusal",
    "self-test-closure-ready-snapshot-mismatch-refusal",
    "self-test-valid-self-test-pass-bundle",
    "self-test-self-test-pass-authority-refusal",
    "self-test-self-test-pass-detail-refusal",
    "self-test-self-test-pass-log-refusal", "self-test-wrong-authority-refusal",
    "self-test-wrong-command-refusal",
    "self-test-supervisor-state-contradiction-refusal",
    "self-test-supervisor-containment-contradiction-refusal",
    "self-test-supervisor-flag-contradiction-refusal",
    "self-test-supervisor-deadline-contradiction-refusal",
    "self-test-proof-boundary-body-refusal",
    "self-test-source-manifest-body-refusal",
    "self-test-failed-control-continuation-refusal",
    "self-test-premature-snapshot-refusal", "self-test-toctou-mutation-refusal",
    "self-test-publication-gap-mutation-refusal",
    "self-test-publication-destination-race-refusal",
    "self-test-mutated-bundle-refusal", "self-test-duplicate-seal-refusal",
    "self-test-truncated-seal-refusal", "self-test-summary-mismatch-refusal",
    "self-test-unsafe-path-refusal", "self-test-duplicate-json-key-refusal",
    "self-test-nonfinite-json-refusal", "self-test-oversized-json-refusal",
    "self-test-record-count-refusal", "self-test-readiness-mismatch-refusal",
    "self-test-no-claim-mutation-refusal", "self-test-unknown-terminal-refusal",
    "self-test-inventory-entry-cap-refusal",
    "self-test-unreferenced-file-refusal",
    "self-test-prefix-file-mismatch-refusal",
    "self-test-interrupted-status-matrix-refusal",
    "self-test-incomplete-status-matrix-refusal",
    "self-test-closure-refusal-writer",
    "self-test-closure-infrastructure-failure",
    "self-test-aggregate-budget-exhaustion",
    "self-test-invalid-candidate-not-published", "self-test-incomplete-seal",
    "self-test-publication-fsync-failure-retains-bundle",
    "self-test-prepublication-signal", "self-test-publication-window-hup",
    "self-test-publication-window-int", "self-test-publication-window-term",
    "self-test-wrapper-arbitrary-term", "self-test-wrapper-hup",
    "self-test-wrapper-int", "self-test-wrapper-term",
)
self_test_log_label_overrides = {
    "self-test-special-index-flag-refusal": "special-index-flag",
    "self-test-skip-worktree-flag-refusal": "skip-worktree-flag",
    "self-test-fsmonitor-valid-flag-refusal": "fsmonitor-valid-flag",
    "self-test-zero-smoke-refusal": "zero-smoke-sentinel",
    "self-test-wrong-authority-refusal": "wrong-authority",
    "self-test-wrong-command-refusal": "wrong-command",
    "self-test-supervisor-state-contradiction-refusal": "supervisor-state-contradiction-refusal",
    "self-test-proof-boundary-body-refusal": "proof-boundary-body-refusal",
    "self-test-source-manifest-body-refusal": "source-manifest-body-refusal",
    "self-test-failed-control-continuation-refusal": "failed-control-continuation",
    "self-test-premature-snapshot-refusal": "premature-snapshot",
    "self-test-toctou-mutation-refusal": "toctou-mutation",
    "self-test-publication-destination-race-refusal": "publication-destination-race-refusal",
    "self-test-mutated-bundle-refusal": "mutated-bundle",
    "self-test-duplicate-seal-refusal": "duplicate-seal",
    "self-test-truncated-seal-refusal": "truncated-seal",
    "self-test-summary-mismatch-refusal": "summary-mismatch",
    "self-test-unsafe-path-refusal": "unsafe-path",
    "self-test-duplicate-json-key-refusal": "duplicate-json-key",
    "self-test-nonfinite-json-refusal": "nonfinite-json",
    "self-test-oversized-json-refusal": "oversized-json",
    "self-test-record-count-refusal": "record-count",
    "self-test-readiness-mismatch-refusal": "readiness-mismatch",
    "self-test-no-claim-mutation-refusal": "no-claim-mutation",
    "self-test-unknown-terminal-refusal": "unknown-terminal",
    "self-test-unreferenced-file-refusal": "unexpected-file",
    "self-test-prefix-file-mismatch-refusal": "prefix-mismatch",
    "self-test-interrupted-status-matrix-refusal": "invalid-interrupted-matrix",
    "self-test-incomplete-status-matrix-refusal": "invalid-incomplete-matrix",
}

def build_self_test_bundle(destination):
    (destination / "logs").mkdir(parents=True)
    records = []
    for lane_name in self_test_sequence:
        log_label = self_test_log_label_overrides.get(
            lane_name, lane_name.removeprefix("self-test-")
        )
        data = (
            f"self-test-label={log_label}\n"
            "self-test-assertion-result=pass\n"
        ).encode()
        records.append(make_record(
            destination,
            lane_name,
            "PASS",
            "harness self-test assertion matched expected disposition",
            data=data,
            authority="harness-self-test",
            executor="self-test-no-cargo",
            log_name=f"self-test-{log_label}",
        ))
    self_test_seal = base_seal(
        "SELF_TEST_PASS",
        len(records),
        0,
        0,
        provenance="stable",
        executor="self-test-no-cargo",
        coverage="not-applicable",
        proof_scope="harness-self-test-no-cargo",
    )
    write_bundle(destination, records, self_test_seal)

build_self_test_bundle(valid_self_test)
for destination in (
    self_test_wrong_authority,
    self_test_wrong_detail,
    self_test_wrong_log,
):
    shutil.copytree(valid_self_test, destination)

authority_records, authority_seal = read_bundle(self_test_wrong_authority)
authority_records[0]["authority"] = "consistent-rehash-overclaim"
rewrite_records(self_test_wrong_authority, authority_records, authority_seal)

detail_records, detail_seal = read_bundle(self_test_wrong_detail)
detail_records[0]["detail"] = "consistently rehashed but semantically false"
rewrite_records(self_test_wrong_detail, detail_records, detail_seal)

log_records, log_seal = read_bundle(self_test_wrong_log)
old_log = self_test_wrong_log / log_records[0]["log"]
log_records[0]["log"] = "logs/consistently-rehashed-wrong-locator.log"
old_log.rename(self_test_wrong_log / log_records[0]["log"])
rewrite_records(self_test_wrong_log, log_records, log_seal)

for destination in (readiness_mismatch, no_claim_mutation):
    shutil.copytree(valid_success, destination)
rewrite_terminal(readiness_mismatch, status="READY_FOR_DSR", candidate_ready_for_dsr=False)
rewrite_terminal(no_claim_mutation, no_claim="physical authority silently widened")
rewrite_terminal(unknown_terminal, status="UNKNOWN_TERMINAL_STATE")
(unexpected_file / "unexpected.txt").write_text("unsealed extra evidence\n")
with (prefix_mismatch / "verdicts-prefix.jsonl").open("ab") as handle:
    handle.write(b"prefix mismatch\n")

wrong_authority_records, wrong_authority_seal = read_bundle(wrong_authority)
wrong_authority_records[-1]["authority"] = "self-test-overclaim"
rewrite_records(wrong_authority, wrong_authority_records, wrong_authority_seal)

wrong_command_records, wrong_command_seal = read_bundle(wrong_command)
wrong_command_log = command_log(
    "constellation-verify",
    ["scripts/ci/checkout_constellation.sh", "--snapshot"],
    rc=2,
)
(wrong_command / "logs" / "constellation-verify.log").write_bytes(wrong_command_log)
wrong_command_records[-1]["log_bytes"] = len(wrong_command_log)
wrong_command_records[-1]["log_sha256"] = hashlib.sha256(wrong_command_log).hexdigest()
rewrite_records(wrong_command, wrong_command_records, wrong_command_seal)

shutil.copytree(valid_success, control_continuation)
continuation_records, continuation_seal = read_bundle(control_continuation)
continuation_records[2]["status"] = "FAIL"
continuation_records[2]["detail"] = "command exited 2"
continuation_log = command_log("constellation-verify", rc=2)
(control_continuation / "logs" / "constellation-verify.log").write_bytes(
    continuation_log
)
continuation_records[2]["log_bytes"] = len(continuation_log)
continuation_records[2]["log_sha256"] = hashlib.sha256(continuation_log).hexdigest()
continuation_seal.update({
    "status": "FAIL",
    "terminal_exit_code": 7,
    "provenance_state": "incomplete",
})
rewrite_records(control_continuation, continuation_records, continuation_seal)

premature_snapshot_bytes = b"snapshot published before its lane passed\n"
(premature_snapshot / "snapshot-before.txt").write_bytes(premature_snapshot_bytes)
rewrite_terminal(
    premature_snapshot,
    snapshot_before_sha256=hashlib.sha256(premature_snapshot_bytes).hexdigest(),
)
shutil.copytree(valid_success, toctou_mutation)

(invalid_interrupted / "logs").mkdir(parents=True)
active_log = b"synthetic interrupted self-test command\n"
(invalid_interrupted / "logs" / "self-test-signal-target-wait.log").write_bytes(active_log)
active_lane = dict(failed_lane)
active_lane.update({
    "lane": "self-test-signal-target-wait",
    "status": "FAIL",
    "authority": "harness-self-test-interrupted",
    "detail": "synthetic interrupted command",
    "log": "logs/self-test-signal-target-wait.log",
    "log_bytes": len(active_log),
    "log_sha256": hashlib.sha256(active_log).hexdigest(),
    "executor_declaration": "self-test-no-cargo",
})
wrapper_log = b"wrapper-signal=15\nno-later-lanes-launched=true\n"
(invalid_interrupted / "logs" / "wrapper-signal.log").write_bytes(wrapper_log)
wrapper_lane = dict(active_lane)
wrapper_lane.update({
    "lane": "wrapper-signal",
    "status": "PASS",
    "authority": "bounded-interrupt-cleanup",
    "detail": "wrapper received signal 15 and stopped scheduling",
    "log": "logs/wrapper-signal.log",
    "log_bytes": len(wrapper_log),
    "log_sha256": hashlib.sha256(wrapper_log).hexdigest(),
})
interrupted_seal = dict(seal)
interrupted_seal.update({
    "status": "INTERRUPTED",
    "checks": 2,
    "failures": 1,
    "executor_declaration": "self-test-no-cargo",
    "provenance_state": "incomplete",
    "terminal_exit_code": 143,
})
write_bundle(invalid_interrupted, [active_lane, wrapper_lane], interrupted_seal)

(invalid_incomplete / "logs").mkdir(parents=True)
for source in ("proof-boundary.log", "source-manifest.log", "constellation-verify.log"):
    shutil.copy2(valid / "logs" / source, invalid_incomplete / "logs" / source)
internal_log = b"unexpected-exit-code=125\nterminal-seal-was-missing=true\n"
(invalid_incomplete / "logs" / "internal-incomplete.log").write_bytes(internal_log)
internal_lane = dict(failed_lane)
internal_lane.update({
    "lane": "internal-incomplete",
    "status": "PASS",
    "authority": "harness-integrity",
    "detail": "EXIT trap sealed an otherwise incomplete run",
    "log": "logs/internal-incomplete.log",
    "log_bytes": len(internal_log),
    "log_sha256": hashlib.sha256(internal_log).hexdigest(),
})
incomplete_seal = dict(seal)
incomplete_seal.update({
    "status": "INCOMPLETE",
    "checks": 4,
    "failures": 1,
    "provenance_state": "incomplete",
    "terminal_exit_code": 125,
})
write_bundle(
    invalid_incomplete,
    [proof_lane, manifest_lane, failed_lane, internal_lane],
    incomplete_seal,
)
PY
}

run_harness_self_test() {
  if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_PUBLICATION_SIGNAL_TARGET:-0}" == 1 ]]; then
    printf '%s\n' \
      'unexpected-exit-code=37' \
      'terminal-seal-was-missing=true' \
      >"$LOG_DIR/logs/internal-incomplete.log"
    record \
      "internal-incomplete" FAIL "harness-integrity" \
      "EXIT trap sealed an otherwise incomplete run" \
      "logs/internal-incomplete.log"
    write_summary_and_seal \
      INCOMPLETE not-checked false focused-software-only incomplete 37
    printf 'publication signal target unexpectedly completed; retained evidence=%s\n' \
      "$LOG_DIR" >&2
    exit 37
  fi

  if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_INCOMPLETE:-0}" == 1 ]]; then
    printf '%s\n' 'self-test injected an incomplete harness exit' >&2
    exit 37
  fi

  if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_SIGNAL_TARGET:-0}" == 1 ]]; then
    self_test_capture \
      "signal-target-wait" 60 1048576 \
      bash -c 'trap "" HUP INT TERM; while :; do sleep 1; done'
    if [[ "$WRAPPER_SIGNAL" != 0 ]]; then
      exit_for_wrapper_signal "$WRAPPER_SIGNAL"
    fi
    printf '%s\n' 'signal-target unexpectedly completed' \
      >>"$SELF_TEST_LAST_LOG"
    record \
      "signal-target-unexpected-completion" FAIL "harness-self-test" \
      "signal target completed without receiving a wrapper signal" \
      "$SELF_TEST_LAST_LOG_REL"
    LANE_COMMIT_CRITICAL=0
    if [[ "$WRAPPER_SIGNAL" != 0 ]]; then
      exit_for_wrapper_signal "$WRAPPER_SIGNAL"
    fi
    write_summary_and_seal \
      SELF_TEST_FAIL not-applicable false harness-self-test-no-cargo incomplete 1
    printf 'Euler-disc harness signal-target self-test FAILED; evidence=%s\n' \
      "$LOG_DIR" >&2
    exit 1
  fi

  if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_ORDINARY_SIGNAL_WINDOW:-0}" == 1 ]]; then
    self_test_capture \
      "ordinary-pass" 60 1048576 \
      bash -c 'printf "ordinary-pass-signal-window\n"; while :; do sleep 1; done'
  else
    self_test_capture "ordinary-pass" 5 1048576 \
      bash -c 'argv_only_marker="must-not-match-command-header"; printf "ordinary-pass\n"'
  fi
  self_test_assert \
    "self-test-ordinary-pass" 0 \
    "ordinary zero exit is retained and assertions ignore escaped argv text" \
    'ordinary-pass' 'process_group_drained=true' \
    'absent:must-not-match-command-header'

  self_test_capture "exact-nonzero" 5 1048576 bash -c 'printf "exact-nonzero\n"; exit 23'
  self_test_assert \
    "self-test-exact-nonzero" 23 \
    "exact nonzero exit is preserved" \
    'exact-nonzero' 'process_group_drained=true'

  self_test_capture \
    "reserved-child-exit" 5 1048576 \
    bash -c 'printf "reserved-child-exit=122\n"; exit 122'
  self_test_append_producer_disposition
  self_test_assert \
    "self-test-reserved-child-exit" 122 \
    "a child exit that collides with a reserved wrapper code remains an ordinary child exit" \
    'reserved-child-exit=122' '"shutdown_reason":"leader-exit"' \
    'producer-shutdown-reason=leader-exit' \
    'producer-detail=command exited 122' 'producer-abort=0'

  self_test_capture \
    "child-signal-exit" 5 1048576 \
    bash -c 'printf "child-signal=TERM\n"; kill -TERM "$$"'
  self_test_append_producer_disposition
  self_test_assert \
    "self-test-child-signal-exit" 143 \
    "a signal delivered to the child is not misclassified as a wrapper interrupt" \
    'child-signal=TERM' '"shutdown_reason":"leader-signal"' \
    'producer-shutdown-reason=leader-signal' \
    'producer-detail=command exited 143' 'producer-abort=0'

  self_test_capture "spawn-failure" 5 1048576 \
    /definitely-not-an-euler-disc-self-test-command
  self_test_assert \
    "self-test-spawn-failure" 127 \
    "missing executable is classified without Cargo" \
    'launch failure:'

  self_test_capture \
    "timeout-stubborn-group" 1 1048576 \
    bash -c 'trap "" TERM; (trap "" TERM; while :; do sleep 1; done) & while :; do sleep 1; done'
  self_test_assert \
    "self-test-timeout-stubborn-group" 124 \
    "timeout drains a stubborn leader and grandchild" \
    'regex:monotonic timeout at (lane|aggregate) deadline [0-9]+' 'process_group_drained=true' \
    'regex:observed_descendants=[1-9][0-9]*'

  self_test_capture \
    "leader-exit-live-group" 5 1048576 \
    bash -c '(trap "" TERM; exec >/dev/null 2>&1; while :; do sleep 1; done) & exit 0'
  self_test_assert \
    "self-test-leader-exit-live-group" 123 \
    "leader exit with a stdio-closed descendant is never PASS" \
    'leader exited while live same-session descendants remained' \
    'process_group_drained=true'

  self_test_capture \
    "setpgid-escape-drained" 5 1048576 \
    python3 -c $'import os\nimport signal\nimport time\nchild = os.fork()\nif child:\n    os._exit(0)\nos.setpgid(0, 0)\nsignal.signal(signal.SIGTERM, signal.SIG_IGN)\nprint(f"escaped-pgid={os.getpgrp()} retained-sid={os.getsid(0)}", flush=True)\nwhile True:\n    time.sleep(1)'
  self_test_assert_external_session_empty
  self_test_assert \
    "self-test-setpgid-escape-drained" 123 \
    "a same-session descendant cannot escape cleanup by changing process group" \
    'escaped-pgid=' 'retained-sid=' \
    'process_session_drained=true' \
    'external_process_session_drained=true'

  self_test_capture \
    "output-truncation-drained" 5 64 \
    bash -c 'for _ in {1..200}; do printf 1234567890; done'
  self_test_assert \
    "self-test-output-truncation-drained" 122 \
    "bounded output loss is distinct from unestablished process-group drain" \
    'output cap reached at 64 bytes' 'process_group_drained=true' \
    "retained_log_hard_limit=$HARD_MAX_RETAINED_LOG_BYTES"

  self_test_capture \
    "cap-plus-timeout" 1 64 \
    bash -c 'for _ in {1..200}; do printf 1234567890; done; trap "" TERM; while :; do sleep 1; done'
  self_test_assert \
    "self-test-cap-plus-timeout" 124 \
    "timeout and truncation remain separately observable" \
    'output cap reached at 64 bytes' \
    'regex:monotonic timeout at (lane|aggregate) deadline [0-9]+' \
    'process_group_drained=true'

  # The nested bash, not this wrapper, evaluates the helper expression.
  # shellcheck disable=SC2016
  self_test_capture \
    "lane-log-cap-boundary" 5 1048576 \
    bash -c \
    'max="$1"; if lane_log_cap_is_valid "$max" && ! lane_log_cap_is_valid "$((max + 1))"; then printf "max-accepted=%s max-plus-one-refused=%s\n" "$max" "$((max + 1))"; else exit 1; fi' \
    _ "$HARD_MAX_LANE_LOG_BYTES"
  self_test_assert \
    "self-test-lane-log-cap-boundary" 0 \
    "the exact safe child-output maximum is accepted and maximum-plus-one is refused" \
    "max-accepted=$HARD_MAX_LANE_LOG_BYTES" \
    "max-plus-one-refused=$((HARD_MAX_LANE_LOG_BYTES + 1))" \
    'process_group_drained=true'

  # Positional parameters and arithmetic are intentionally evaluated by the
  # bounded child, not by this wrapper.
  # shellcheck disable=SC2016
  self_test_capture \
    "numeric-config-boundaries" 20 1048576 \
    bash -c '
probe() {
  expected="$1" label="$2"
  shift 2
  "$@"
  observed=$?
  printf "numeric-probe=%s expected=%s observed=%s\n" "$label" "$expected" "$observed"
  [[ "$observed" == "$expected" ]]
}
script="$1" timeout_max="$2" log_max="$3" huge="$4"
probe 37 exact-max env \
  FSIM_EULER_DISC_E2E_LANE_TIMEOUT_SECONDS="$timeout_max" \
  FSIM_EULER_DISC_E2E_RUN_TIMEOUT_SECONDS="$timeout_max" \
  FSIM_EULER_DISC_E2E_LANE_LOG_MAX_BYTES="$log_max" \
  FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_INCOMPLETE=1 "$script" --self-test || exit 1
probe 2 lane-max-plus-one env \
  FSIM_EULER_DISC_E2E_LANE_TIMEOUT_SECONDS="$((timeout_max + 1))" "$script" --self-test || exit 1
probe 2 run-max-plus-one env \
  FSIM_EULER_DISC_E2E_RUN_TIMEOUT_SECONDS="$((timeout_max + 1))" "$script" --self-test || exit 1
probe 2 log-max-plus-one env \
  FSIM_EULER_DISC_E2E_LANE_LOG_MAX_BYTES="$((log_max + 1))" "$script" --self-test || exit 1
probe 2 lane-huge env FSIM_EULER_DISC_E2E_LANE_TIMEOUT_SECONDS="$huge" "$script" --self-test || exit 1
probe 2 run-huge env FSIM_EULER_DISC_E2E_RUN_TIMEOUT_SECONDS="$huge" "$script" --self-test || exit 1
probe 2 log-huge env FSIM_EULER_DISC_E2E_LANE_LOG_MAX_BYTES="$huge" "$script" --self-test || exit 1
probe 2 empty-verify "$script" --verify-bundle "" || exit 1
' _ "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    "$PROOF_MAX_TIMEOUT_SECONDS" "$HARD_MAX_LANE_LOG_BYTES" \
    9999999999999999999999999999999999999999
  local expired_deadline_ns near_deadline_ns deadline_before_launch_rc deadline_during_run_rc
  local completion_after_deadline_rc
  expired_deadline_ns="$(python3 - <<'PY'
import time
print(time.monotonic_ns() - 1)
PY
)"
  if run_bounded_command \
    "$SELF_TEST_LAST_LOG" 5 "$expired_deadline_ns" \
    1048576 1048576 "$HARD_MAX_RETAINED_LOG_BYTES" none \
    /definitely-not-launched-after-aggregate-deadline; then
    deadline_before_launch_rc=0
  else
    deadline_before_launch_rc=$?
  fi
  near_deadline_ns="$(python3 - <<'PY'
import time
print(time.monotonic_ns() + 200_000_000)
PY
)"
  if run_bounded_command \
    "$SELF_TEST_LAST_LOG" 5 "$near_deadline_ns" \
    1048576 1048576 "$HARD_MAX_RETAINED_LOG_BYTES" none \
    bash -c 'trap "" TERM; while :; do sleep 1; done'; then
    deadline_during_run_rc=0
  else
    deadline_during_run_rc=$?
  fi
  export FSIM_EULER_DISC_E2E_SELF_TEST_DELAY_COMPLETION_CLASSIFICATION=1
  if run_bounded_command \
    "$SELF_TEST_LAST_LOG" 1 "$RUN_DEADLINE_MONOTONIC_NS" \
    1048576 1048576 "$HARD_MAX_RETAINED_LOG_BYTES" none \
    bash -c 'printf "completion-deadline-candidate\n"'; then
    completion_after_deadline_rc=0
  else
    completion_after_deadline_rc=$?
  fi
  unset FSIM_EULER_DISC_E2E_SELF_TEST_DELAY_COMPLETION_CLASSIFICATION
  printf '%s\n' \
    "aggregate-deadline-before-launch-rc=$deadline_before_launch_rc" \
    "aggregate-deadline-during-run-rc=$deadline_during_run_rc" \
    "completion-after-deadline-rc=$completion_after_deadline_rc" \
    >>"$SELF_TEST_LAST_LOG"
  if [[ "$deadline_before_launch_rc" != 121 \
      || "$deadline_during_run_rc" != 124 \
      || "$completion_after_deadline_rc" != 124 ]]; then
    SELF_TEST_LAST_RC=1
  fi
  self_test_assert \
    "self-test-numeric-config-boundaries" 0 \
    "all numeric settings accept their exact maximum and refuse maximum-plus-one and huge text before Bash arithmetic" \
    'numeric-probe=exact-max expected=37 observed=37' \
    'numeric-probe=lane-max-plus-one expected=2 observed=2' \
    'numeric-probe=run-max-plus-one expected=2 observed=2' \
    'numeric-probe=log-max-plus-one expected=2 observed=2' \
    'numeric-probe=lane-huge expected=2 observed=2' \
    'numeric-probe=run-huge expected=2 observed=2' \
    'numeric-probe=log-huge expected=2 observed=2' \
    'numeric-probe=empty-verify expected=2 observed=2' \
    '--verify-bundle requires a nonempty proof-bundle directory' \
    'aggregate-deadline-before-launch-rc=121' \
    'aggregate-deadline-during-run-rc=124' \
    'completion-deadline-candidate' \
    'completion-after-deadline-rc=124' \
    'completion_delay_hook_fired=true' \
    '"deadline_kind":"aggregate"' \
    '"shutdown_reason":"aggregate-deadline-before-launch"' \
    'process_group_drained=true'

  export FSIM_EULER_DISC_E2E_SELF_TEST_SIGNAL_COMPLETION_CLASSIFICATION=1
  self_test_capture \
    "completion-classification-signal" 5 1048576 \
    bash -c 'printf "completion-signal-candidate\n"'
  unset FSIM_EULER_DISC_E2E_SELF_TEST_SIGNAL_COMPLETION_CLASSIFICATION
  self_test_assert \
    "self-test-completion-classification-signal" 143 \
    "a signal latched during final classification outranks ordinary completion" \
    'completion-signal-candidate' 'completion_signal_hook_fired=true' \
    '"interrupted_signal":15' '"shutdown_reason":"interrupt"' \
    'process_group_drained=true'

  local exact_snapshot_candidate="$LOG_ROOT/.snapshot-exact-${HEAD_SHA:0:12}-XXXXXXXX"
  local overflow_snapshot_candidate="$LOG_ROOT/.snapshot-overflow-${HEAD_SHA:0:12}-XXXXXXXX"
  exact_snapshot_candidate="$(mktemp "$exact_snapshot_candidate")"
  overflow_snapshot_candidate="$(mktemp "$overflow_snapshot_candidate")"
  # Positional parameters are intentionally evaluated by the bounded child.
  # shellcheck disable=SC2016
  self_test_capture \
    "snapshot-byte-boundaries" 5 1048576 \
    bash -c '
set -o pipefail
python3 -c '\''import sys; sys.stdout.buffer.write(b"x" * 64)'\'' | \
  bounded_snapshot_append "$1" 64 || exit 1
exact_bytes="$(wc -c <"$1" | tr -d " ")"
python3 -c '\''import sys; sys.stdout.buffer.write(b"y" * 65)'\'' | \
  bounded_snapshot_append "$2" 64
overflow_rc=$?
overflow_bytes="$(wc -c <"$2" | tr -d " ")"
printf "snapshot-exact-bytes=%s snapshot-overflow-rc=%s snapshot-overflow-bytes=%s\n" \
  "$exact_bytes" "$overflow_rc" "$overflow_bytes"
[[ "$exact_bytes" == 64 && "$overflow_rc" == 122 && "$overflow_bytes" == 64 ]]
' _ "$exact_snapshot_candidate" "$overflow_snapshot_candidate"
  self_test_assert \
    "self-test-snapshot-byte-boundaries" 0 \
    "snapshot streaming accepts the exact cap and retains a bounded nonauthoritative candidate on maximum-plus-one" \
    'snapshot-exact-bytes=64 snapshot-overflow-rc=122 snapshot-overflow-bytes=64' \
    'snapshot retained-byte bound exceeded: maximum=64' \
    'process_group_drained=true'

  local snapshot_hash_failure_candidate
  snapshot_hash_failure_candidate="$(
    mktemp "$LOG_ROOT/.snapshot-hash-failure-${HEAD_SHA:0:12}-XXXXXXXX"
  )"
  printf 'snapshot hash-failure candidate retained outside proof bundle: %s\n' \
    "$snapshot_hash_failure_candidate" >&2
  # The bounded child intentionally replaces only the exported hash helper.
  # shellcheck disable=SC2016
  self_test_capture \
    "snapshot-hash-failure" 5 1048576 \
    bash -c '
sha256_file() {
  printf "injected-snapshot-hash-read-failure\n" >&2
  return 74
}
append_snapshot_scope_hash "$1" "$2" "$3"
' _ "$snapshot_hash_failure_candidate" \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    "$PROOF_MAX_SNAPSHOT_BYTES"
  self_test_assert \
    "self-test-snapshot-hash-failure" 74 \
    "snapshot provenance refuses and preserves the real hash-read failure exit" \
    'injected-snapshot-hash-read-failure' \
    'snapshot-scope-hash-failed=' 'exit=74' \
    '"shutdown_reason":"leader-exit"' 'process_group_drained=true'

  self_test_capture \
    "timeout-hanging-helper" 1 1048576 \
    bash -c 'printf "helper-kind=snapshot-or-membership\n"; trap "" TERM; while :; do sleep 1; done'
  self_test_assert \
    "self-test-timeout-hanging-helper" 124 \
    "provenance helpers share the bounded supervisor path" \
    'helper-kind=snapshot-or-membership' \
    'regex:monotonic timeout at (lane|aggregate) deadline [0-9]+' \
    'process_group_drained=true'

  export FSIM_EULER_DISC_E2E_SELF_TEST_SUPERVISOR_EXCEPTION=1
  self_test_capture \
    "supervisor-exception" 5 1048576 \
    bash -c 'trap "" TERM; (trap "" TERM; while :; do sleep 1; done) & while :; do sleep 1; done'
  unset FSIM_EULER_DISC_E2E_SELF_TEST_SUPERVISOR_EXCEPTION
  self_test_assert_external_session_empty
  self_test_assert \
    "self-test-supervisor-exception" 1 \
    "an unexpected supervisor exception kills, reaps, and independently drains its owned session" \
    'injected supervisor exception for no-Cargo self-test' \
    'leader_reaped=true' 'process_group_drained=true' \
    'process_session_drained=true' 'external_process_session_drained=true' \
    'regex:observed_descendants=[1-9][0-9]*'

  local self_test_aux_root
  self_test_aux_root="$(
    mktemp -d "$LOG_ROOT/.self-test-aux-${HEAD_SHA:0:12}-XXXXXXXX"
  )"
  printf 'self-test auxiliary fixtures retained outside proof bundle: %s\n' \
    "$self_test_aux_root" >&2

  local special_index_repo="$self_test_aux_root/assume-unchanged-fixture"
  mkdir -p "$special_index_repo"
  git -c init.templateDir= -C "$special_index_repo" init -q -b main
  printf '%s\n' 'committed bytes' >"$special_index_repo/tracked.txt"
  git -C "$special_index_repo" add tracked.txt
  git -C "$special_index_repo" \
    -c core.hooksPath=/dev/null \
    -c commit.gpgSign=false \
    -c user.name='FrankenSim Harness' \
    -c user.email='harness@invalid.example' \
    commit -q -m 'self-test fixture'
  git -C "$special_index_repo" update-index --assume-unchanged tracked.txt
  printf '%s\n' 'concealed worktree mutation' >"$special_index_repo/tracked.txt"
  local special_index_head
  special_index_head="$(git -C "$special_index_repo" rev-parse HEAD)"
  # Positional parameters are intentionally expanded by the bounded child.
  # shellcheck disable=SC2016
  self_test_capture \
    "special-index-flag" 5 1048576 \
    bash -c 'cd "$1"; closure_root_preflight_command "$2" tracked.txt' _ \
    "$special_index_repo" "$special_index_head"
  self_test_assert \
    "self-test-special-index-flag-refusal" 1 \
    "closure refuses assume-unchanged concealment" \
    'repository-status=clean' 'special-index-flags=present' 'h tracked.txt'

  local skip_worktree_repo="$self_test_aux_root/skip-worktree-fixture"
  mkdir -p "$skip_worktree_repo"
  git -c init.templateDir= -C "$skip_worktree_repo" init -q -b main
  printf '%s\n' 'committed bytes' >"$skip_worktree_repo/tracked.txt"
  git -C "$skip_worktree_repo" add tracked.txt
  git -C "$skip_worktree_repo" \
    -c core.hooksPath=/dev/null \
    -c commit.gpgSign=false \
    -c user.name='FrankenSim Harness' \
    -c user.email='harness@invalid.example' \
    commit -q -m 'self-test fixture'
  git -C "$skip_worktree_repo" update-index --skip-worktree tracked.txt
  printf '%s\n' 'concealed worktree mutation' >"$skip_worktree_repo/tracked.txt"
  local skip_worktree_head
  skip_worktree_head="$(git -C "$skip_worktree_repo" rev-parse HEAD)"
  # Positional parameters are intentionally expanded by the bounded child.
  # shellcheck disable=SC2016
  self_test_capture \
    "skip-worktree-flag" 5 1048576 \
    bash -c 'cd "$1"; closure_root_preflight_command "$2" tracked.txt' _ \
    "$skip_worktree_repo" "$skip_worktree_head"
  self_test_assert \
    "self-test-skip-worktree-flag-refusal" 1 \
    "closure refuses skip-worktree concealment" \
    'repository-status=clean' 'special-index-flags=present' 'S tracked.txt'

  local fsmonitor_valid_repo="$self_test_aux_root/fsmonitor-valid-fixture"
  mkdir -p "$fsmonitor_valid_repo"
  git -c init.templateDir= -C "$fsmonitor_valid_repo" init -q -b main
  printf '%s\n' 'committed bytes' >"$fsmonitor_valid_repo/tracked.txt"
  git -C "$fsmonitor_valid_repo" add tracked.txt
  git -C "$fsmonitor_valid_repo" \
    -c core.hooksPath=/dev/null \
    -c commit.gpgSign=false \
    -c user.name='FrankenSim Harness' \
    -c user.email='harness@invalid.example' \
    commit -q -m 'self-test fixture'
  local fsmonitor_backend=hook
  if [[ "$(uname -s)" == Darwin ]]; then
    # Apple Git can suspend freshly created hook interpreters during its
    # fsmonitor probe. Use the supported built-in daemon, then stop it before
    # leaving the fixture so it cannot contaminate session-drain evidence.
    fsmonitor_backend=builtin
    git -C "$fsmonitor_valid_repo" config core.fsmonitor true
  else
    local fsmonitor_hook="$fsmonitor_valid_repo/.git/fsmonitor-self-test"
    printf '%s\n' \
      '#!/usr/bin/env python3' \
      'import sys' \
      'sys.stdout.buffer.write(b"self-test-token\x00")' \
      >"$fsmonitor_hook"
    chmod 700 "$fsmonitor_hook"
    git -C "$fsmonitor_valid_repo" config core.fsmonitor "$fsmonitor_hook"
    git -C "$fsmonitor_valid_repo" config core.fsmonitorHookVersion 2
  fi
  git -C "$fsmonitor_valid_repo" update-index --fsmonitor-valid tracked.txt
  local fsmonitor_valid_head
  fsmonitor_valid_head="$(git -C "$fsmonitor_valid_repo" rev-parse HEAD)"
  # Positional parameters are intentionally expanded by the bounded child.
  # shellcheck disable=SC2016
  self_test_capture \
    "fsmonitor-valid-flag" 5 1048576 \
    bash -c 'cd "$1"; closure_root_preflight_command "$2" tracked.txt' _ \
    "$fsmonitor_valid_repo" "$fsmonitor_valid_head"
  local fsmonitor_cleanup_status=not-applicable
  local fsmonitor_cleanup_expected='fsmonitor-backend=hook cleanup=not-applicable'
  if [[ "$fsmonitor_backend" == builtin ]]; then
    if git -C "$fsmonitor_valid_repo" fsmonitor--daemon stop \
      >>"$SELF_TEST_LAST_LOG" 2>&1; then
      fsmonitor_cleanup_status=complete
    else
      fsmonitor_cleanup_status=failed
      SELF_TEST_LAST_RC=1
    fi
    fsmonitor_cleanup_expected='fsmonitor-backend=builtin cleanup=complete'
  fi
  printf 'fsmonitor-backend=%s cleanup=%s\n' \
    "$fsmonitor_backend" "$fsmonitor_cleanup_status" >>"$SELF_TEST_LAST_LOG"
  self_test_assert \
    "self-test-fsmonitor-valid-flag-refusal" 1 \
    "closure refuses fsmonitor-valid index concealment" \
    'repository-status=clean' 'fsmonitor-valid-flags=present' 'h tracked.txt' \
    "$fsmonitor_cleanup_expected"

  local non_repository_fixture="$self_test_aux_root/non-repository-fixture"
  mkdir -p "$non_repository_fixture"
  # Positional parameters are intentionally expanded by the bounded child.
  # shellcheck disable=SC2016
  self_test_capture \
    "real-repository-read-failure" 5 1048576 \
    bash -c 'export GIT_CEILING_DIRECTORIES="$2"; cd "$1"; closure_root_preflight_command synthetic-head tracked.txt' _ \
    "$non_repository_fixture" "$self_test_aux_root"
  self_test_assert \
    "self-test-real-repository-read-failure" 128 \
    "a real Git repository-read failure preserves its non-policy exit code" \
    'fatal: not a git repository' \
    'infrastructure-command=git-rev-parse-head exit=128'

  local zero_smoke_input="$self_test_aux_root/zero-smoke-input.log"
  printf '%s\n' 'running 0 tests' >"$zero_smoke_input"
  self_test_capture \
    "zero-smoke-sentinel" 5 1048576 \
    bash -c 'checker_smoke_sentinel_command "$@"' _ \
    "$zero_smoke_input" exact_test_name
  self_test_assert \
    "self-test-zero-smoke-refusal" 1 \
    "zero-match exact test filters are refused" \
    'expected exactly one passing sentinel'

  local invalid_consolidation_repo="$self_test_aux_root/invalid-consolidation-repo"
  local invalid_consolidation_bundles="$self_test_aux_root/invalid-consolidation-bundles"
  self_test_capture \
    "invalid-consolidation-disposition-refusal" \
    "$NORMAL_HARNESS_HOSTILE_SELF_TEST_TIMEOUT_SECONDS" 1048576 \
    bash -c 'normal_harness_hostile_self_test_command "$@"' _ \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    "$invalid_consolidation_repo" "$invalid_consolidation_bundles" \
    invalid-consolidation
  self_test_verify_nested_bundle \
    'proof-bundle published: ' "$invalid_consolidation_bundles"
  self_test_assert \
    "self-test-invalid-consolidation-disposition-refusal" 1 \
    "an invalid Euler consolidation disposition fails its required structural lane and cannot publish focused success" \
    'fixture-status-before= M consolidation-review.json' \
    '[FAIL   ] xtask-check-consolidation' \
    'fake-check-consolidation=invalid-euler-root-count count=0' \
    'published-success-count=0' 'published-statuses=FAIL' \
    "nested-configured-lane-timeout-seconds=$NORMAL_HARNESS_NESTED_LANE_TIMEOUT_SECONDS" \
    'proof-bundle verified:'

  local consolidation_mutation_repo="$self_test_aux_root/consolidation-mutation-repo"
  local consolidation_mutation_bundles="$self_test_aux_root/consolidation-mutation-bundles"
  self_test_capture \
    "consolidation-scope-mutation-refusal" \
    "$NORMAL_HARNESS_HOSTILE_SELF_TEST_TIMEOUT_SECONDS" 1048576 \
    bash -c 'normal_harness_hostile_self_test_command "$@"' _ \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    "$consolidation_mutation_repo" "$consolidation_mutation_bundles" \
    scope-mutation
  self_test_verify_nested_bundle \
    'proof-bundle published: ' "$consolidation_mutation_bundles"
  self_test_assert \
    "self-test-consolidation-scope-mutation-refusal" 1 \
    "consolidation bytes changed between snapshots are refused even while Git status remains M at both bookends" \
    'fixture-status-before= M consolidation-review.json' \
    'fixture-status-after= M consolidation-review.json' \
    'fake-consolidation-mutation=before-to-after' \
    '[FAIL   ] source-stability' \
    'published-success-count=0' 'published-statuses=FAIL' \
    "nested-configured-lane-timeout-seconds=$NORMAL_HARNESS_NESTED_LANE_TIMEOUT_SECONDS" \
    'proof-bundle verified:'

  local success_deadline_repo="$self_test_aux_root/success-deadline-repo"
  local success_deadline_bundles="$self_test_aux_root/success-deadline-bundles"
  self_test_capture \
    "success-deadline-publication-refusal" \
    "$NORMAL_HARNESS_HOSTILE_SELF_TEST_TIMEOUT_SECONDS" 1048576 \
    bash -c 'normal_harness_hostile_self_test_command "$@"' _ \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    "$success_deadline_repo" "$success_deadline_bundles" deadline-expiry
  self_test_verify_nested_bundle \
    'Euler-disc harness incomplete; retained evidence=' \
    "$success_deadline_bundles"
  self_test_assert \
    "self-test-success-deadline-publication-refusal" 124 \
    "expiry during a deterministic prepublication delay refuses success while later incomplete sealing remains available" \
    'self-test prepublication deadline delay complete' \
    'aggregate success deadline expired at stage=before-success-publication' \
    'verified success candidate rejected before publication:' \
    'published-success-count=0' 'published-statuses=INCOMPLETE' \
    "nested-configured-lane-timeout-seconds=$NORMAL_HARNESS_NESTED_LANE_TIMEOUT_SECONDS" \
    'proof-bundle verified:'

  local post_deadline_repo="$self_test_aux_root/postpublication-deadline-repo"
  local post_deadline_bundles="$self_test_aux_root/postpublication-deadline-bundles"
  self_test_capture \
    "postpublication-deadline-retraction" \
    "$NORMAL_HARNESS_HOSTILE_SELF_TEST_TIMEOUT_SECONDS" 1048576 \
    bash -c 'normal_harness_hostile_self_test_command "$@"' _ \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    "$post_deadline_repo" "$post_deadline_bundles" \
    postpublication-deadline-expiry
  self_test_verify_nested_bundle \
    'Euler-disc harness incomplete; retained evidence=' \
    "$post_deadline_bundles"
  self_test_assert \
    "self-test-postpublication-deadline-retraction" 124 \
    "expiry inside supervised postpublication verification retracts public success before incomplete sealing" \
    'self-test postpublication deadline armed:' \
    'self-test postpublication verifier delay observed:' \
    'success finalization log retained: label=postpublication-verification exit=124' \
    'published candidate retracted before terminal transition:' \
    'published-success-count=0' 'published-statuses=INCOMPLETE' \
    "nested-configured-lane-timeout-seconds=$NORMAL_HARNESS_NESTED_LANE_TIMEOUT_SECONDS" \
    'proof-bundle verified:'

  local valid_fixture="$self_test_aux_root/verifier-valid"
  local mutated_fixture="$self_test_aux_root/verifier-mutated-log"
  local duplicate_fixture="$self_test_aux_root/verifier-duplicate-seal"
  local truncated_fixture="$self_test_aux_root/verifier-truncated-seal"
  local summary_mismatch_fixture="$self_test_aux_root/verifier-summary-mismatch"
  local unsafe_path_fixture="$self_test_aux_root/verifier-unsafe-path"
  local duplicate_key_fixture="$self_test_aux_root/verifier-duplicate-json-key"
  local nonfinite_fixture="$self_test_aux_root/verifier-nonfinite-json"
  local oversized_fixture="$self_test_aux_root/verifier-oversized-json"
  local entry_cap_fixture="$self_test_aux_root/verifier-entry-cap"
  local too_many_fixture="$self_test_aux_root/verifier-too-many-records"
  local valid_success_fixture="$self_test_aux_root/verifier-valid-focused-success"
  local readiness_mismatch_fixture="$self_test_aux_root/verifier-readiness-mismatch"
  local no_claim_mutation_fixture="$self_test_aux_root/verifier-no-claim-mutation"
  local unknown_terminal_fixture="$self_test_aux_root/verifier-unknown-terminal"
  local unexpected_file_fixture="$self_test_aux_root/verifier-unexpected-file"
  local prefix_mismatch_fixture="$self_test_aux_root/verifier-prefix-mismatch"
  local invalid_interrupted_fixture="$self_test_aux_root/verifier-invalid-interrupted"
  local invalid_incomplete_fixture="$self_test_aux_root/verifier-invalid-incomplete"
  local wrong_authority_fixture="$self_test_aux_root/verifier-wrong-authority"
  local wrong_command_fixture="$self_test_aux_root/verifier-wrong-command"
  local control_continuation_fixture="$self_test_aux_root/verifier-control-continuation"
  local premature_snapshot_fixture="$self_test_aux_root/verifier-premature-snapshot"
  local toctou_mutation_fixture="$self_test_aux_root/verifier-toctou-mutation"
  local proof_body_fixture="$self_test_aux_root/verifier-proof-boundary-body"
  local focused_source_body_fixture="$self_test_aux_root/verifier-focused-source-manifest-body"
  local closure_source_body_fixture="$self_test_aux_root/verifier-closure-source-manifest-body"
  local supervisor_contradiction_fixture="$self_test_aux_root/verifier-supervisor-contradiction"
  local supervisor_containment_fixture="$self_test_aux_root/verifier-supervisor-containment"
  local supervisor_flag_fixture="$self_test_aux_root/verifier-supervisor-flags"
  local supervisor_log_cap_fixture="$self_test_aux_root/verifier-supervisor-log-cap"
  local supervisor_retained_cap_fixture="$self_test_aux_root/verifier-supervisor-retained-cap"
  local supervisor_deadline_fixture="$self_test_aux_root/verifier-supervisor-deadline"
  local valid_closure_fixture="$self_test_aux_root/verifier-valid-closure-ready"
  local closure_missing_lane_fixture="$self_test_aux_root/verifier-closure-missing-lane"
  local closure_snapshot_mismatch_fixture="$self_test_aux_root/verifier-closure-snapshot-mismatch"
  local valid_self_test_fixture="$self_test_aux_root/verifier-valid-self-test-pass"
  local self_test_wrong_authority_fixture="$self_test_aux_root/verifier-self-test-wrong-authority"
  local self_test_wrong_detail_fixture="$self_test_aux_root/verifier-self-test-wrong-detail"
  local self_test_wrong_log_fixture="$self_test_aux_root/verifier-self-test-wrong-log"
  create_verifier_self_test_fixtures \
    "$valid_fixture" "$mutated_fixture" "$duplicate_fixture" \
    "$truncated_fixture" "$summary_mismatch_fixture" "$unsafe_path_fixture" \
    "$duplicate_key_fixture" "$nonfinite_fixture" "$oversized_fixture" \
    "$entry_cap_fixture" "$too_many_fixture" "$valid_success_fixture" \
    "$readiness_mismatch_fixture" "$no_claim_mutation_fixture" \
    "$unknown_terminal_fixture" "$unexpected_file_fixture" \
    "$prefix_mismatch_fixture" "$invalid_interrupted_fixture" \
    "$invalid_incomplete_fixture" "$wrong_authority_fixture" \
    "$wrong_command_fixture" "$control_continuation_fixture" \
    "$premature_snapshot_fixture" "$toctou_mutation_fixture" \
    "$proof_body_fixture" "$focused_source_body_fixture" \
    "$closure_source_body_fixture" "$supervisor_contradiction_fixture" \
    "$supervisor_containment_fixture" "$supervisor_flag_fixture" \
    "$supervisor_log_cap_fixture" "$supervisor_retained_cap_fixture" \
    "$supervisor_deadline_fixture" \
    "$valid_closure_fixture" "$closure_missing_lane_fixture" \
    "$closure_snapshot_mismatch_fixture" "$valid_self_test_fixture" \
    "$self_test_wrong_authority_fixture" "$self_test_wrong_detail_fixture" \
    "$self_test_wrong_log_fixture"
  self_test_capture \
    "valid-bundle" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$valid_fixture"
  self_test_assert \
    "self-test-valid-bundle" 0 \
    "valid terminal bundle verifies without Cargo" \
    'proof-bundle verified:'
  self_test_capture \
    "valid-success-bundle" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$valid_success_fixture"
  self_test_assert \
    "self-test-valid-success-bundle" 0 \
    "an exact focused success state verifies without Cargo" \
    'proof-bundle verified:'
  self_test_capture \
    "valid-closure-ready-bundle" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$valid_closure_fixture"
  self_test_assert \
    "self-test-valid-closure-ready-bundle" 0 \
    "an exact synthetic closure READY_FOR_DSR state verifies without claiming DSR execution" \
    'proof-bundle verified:'
  self_test_capture \
    "closure-ready-missing-lane-refusal" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$closure_missing_lane_fixture"
  self_test_assert \
    "self-test-closure-ready-missing-lane-refusal" 1 \
    "READY_FOR_DSR refuses a consistently sealed closure trace missing one required lane" \
    'DSR-ready bundle does not contain the exact closure lane set'
  self_test_capture \
    "closure-ready-snapshot-mismatch-refusal" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$closure_snapshot_mismatch_fixture"
  self_test_assert \
    "self-test-closure-ready-snapshot-mismatch-refusal" 1 \
    "READY_FOR_DSR refuses independently hashed but unequal bookend snapshots" \
    'DSR-ready bundle lacks equal nonempty source snapshots'
  self_test_capture \
    "valid-self-test-pass-bundle" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$valid_self_test_fixture"
  self_test_assert \
    "self-test-valid-self-test-pass-bundle" 0 \
    "the exact synthetic SELF_TEST_PASS authority/detail/log matrix verifies" \
    'proof-bundle verified:'
  self_test_capture \
    "self-test-pass-authority-refusal" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$self_test_wrong_authority_fixture"
  self_test_assert \
    "self-test-self-test-pass-authority-refusal" 1 \
    "SELF_TEST_PASS refuses a consistently rehashed lane-authority mutation" \
    'SELF_TEST_PASS lane self-test-ordinary-pass has the wrong authority'
  self_test_capture \
    "self-test-pass-detail-refusal" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$self_test_wrong_detail_fixture"
  self_test_assert \
    "self-test-self-test-pass-detail-refusal" 1 \
    "SELF_TEST_PASS refuses a consistently rehashed lane-detail mutation" \
    'SELF_TEST_PASS lane self-test-ordinary-pass has the wrong detail'
  self_test_capture \
    "self-test-pass-log-refusal" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$self_test_wrong_log_fixture"
  self_test_assert \
    "self-test-self-test-pass-log-refusal" 1 \
    "SELF_TEST_PASS refuses a consistently rehashed lane-log-locator mutation" \
    'SELF_TEST_PASS lane self-test-ordinary-pass has the wrong log locator'
  self_test_capture \
    "wrong-authority" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$wrong_authority_fixture"
  self_test_assert \
    "self-test-wrong-authority-refusal" 1 \
    "ordinary lane authority is an exact verifier contract" \
    'normal lane constellation-verify has the wrong authority'

  self_test_capture \
    "wrong-command" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$wrong_command_fixture"
  self_test_assert \
    "self-test-wrong-command-refusal" 1 \
    "canonical argv metadata binds every ordinary lane to its command" \
    'constellation-verify command arguments are not exact'

  self_test_capture \
    "supervisor-state-contradiction-refusal" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$supervisor_contradiction_fixture"
  self_test_assert \
    "self-test-supervisor-state-contradiction-refusal" 1 \
    "canonical supervisor results cannot contradict the lane verdict" \
    'verdict detail contradicts its canonical supervisor result'

  self_test_capture \
    "supervisor-containment-contradiction-refusal" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$supervisor_containment_fixture"
  self_test_assert \
    "self-test-supervisor-containment-contradiction-refusal" 1 \
    "an admitted supervisor disposition requires complete drain, EOF, and metadata" \
    'admitted supervisor disposition lacks complete containment metadata'

  # Positional parameters are intentionally evaluated by the bounded child.
  # shellcheck disable=SC2016
  self_test_capture \
    "supervisor-flag-contradiction-refusal" 5 1048576 \
    bash -c '
set +e
script="$1"
shift
observed=""
for fixture in "$@"; do
  "$script" --verify-bundle "$fixture"
  rc=$?
  printf "supervisor-hostile-fixture=%s exit=%s\n" "$fixture" "$rc"
  [[ "$rc" == 1 ]] || exit 99
  observed="${observed}x"
done
[[ "$observed" == xxx ]] || exit 98
exit 1
' _ "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    "$supervisor_flag_fixture" "$supervisor_log_cap_fixture" \
    "$supervisor_retained_cap_fixture"
  self_test_assert \
    "self-test-supervisor-flag-contradiction-refusal" 1 \
    "writer-impossible supervisor flags and cap arithmetic are refused after consistent rehashing" \
    'supervisor disposition spuriously claims output truncation' \
    'retained log exceeds its supervisor total cap' \
    'supervisor retained-output cap is not the exact bounded remainder' \
    'regex:supervisor-hostile-fixture=.* exit=1'

  self_test_capture \
    "supervisor-deadline-contradiction-refusal" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$supervisor_deadline_fixture"
  self_test_assert \
    "self-test-supervisor-deadline-contradiction-refusal" 1 \
    "a command cannot launch at or after the aggregate run deadline" \
    'launched disposition starts at or after its run deadline'

  self_test_capture \
    "proof-boundary-body-refusal" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$proof_body_fixture"
  self_test_assert \
    "self-test-proof-boundary-body-refusal" 1 \
    "proof-boundary authority requires the exact retained software-only declaration" \
    'proof-boundary log does not match the exact software-only declaration'

  # Both profiles must bind source-manifest NO_DATA to their exact body.
  # shellcheck disable=SC2016
  self_test_capture \
    "source-manifest-body-refusal" 10 1048576 \
    bash -c '
status=0
for bundle in "$2" "$3"; do
  if output="$("$1" --verify-bundle "$bundle" 2>&1)"; then
    printf "unexpected-source-body-verification=%s\n" "$bundle"
    status=1
  else
    rc=$?
    printf "%s\nsource-body-refusal-rc=%s bundle=%s\n" "$output" "$rc" "$bundle"
    [[ "$rc" == 1 ]] || status=1
  fi
done
exit "$status"
' _ "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    "$focused_source_body_fixture" "$closure_source_body_fixture"
  self_test_assert \
    "self-test-source-manifest-body-refusal" 0 \
    "focused and closure source-manifest NO_DATA authority require exact retained bodies" \
    'focused source-manifest has the wrong exact body' \
    'closure source-manifest has the wrong exact body' \
    'source-body-refusal-rc=1'

  self_test_capture \
    "failed-control-continuation" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$control_continuation_fixture"
  self_test_assert \
    "self-test-failed-control-continuation-refusal" 1 \
    "a failed control lane must stop the producer trace" \
    'normal producer trace continues after failed control lane constellation-verify'

  self_test_capture \
    "premature-snapshot" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$premature_snapshot_fixture"
  self_test_assert \
    "self-test-premature-snapshot-refusal" 1 \
    "snapshot authority is absent until its capture lane passes" \
    'source snapshots exist before snapshot-before was attempted'

  local toctou_ready_marker="$self_test_aux_root/verifier-toctou-ready.txt"
  self_test_capture \
    "toctou-mutation" 5 1048576 \
    bash -c 'verifier_toctou_mutation_self_test_command "$@"' _ \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    "$toctou_mutation_fixture" "$toctou_ready_marker"
  self_test_assert \
    "self-test-toctou-mutation-refusal" 1 \
    "root-bound final inventory and rehash refuse concurrent bundle mutation" \
    'bundle inventory or entry metadata changed during verification'

  local publication_gap_root="$self_test_aux_root/publication-gap"
  local publication_gap_marker="$self_test_aux_root/publication-gap-ready.txt"
  mkdir -p "$publication_gap_root"
  self_test_capture \
    "publication-gap-mutation-refusal" 30 1048576 \
    bash -c 'publication_gap_mutation_self_test_command "$@"' _ \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    "$publication_gap_root" "$publication_gap_marker"
  self_test_verify_nested_bundle \
    'Euler-disc harness incomplete; retained evidence=' \
    "$publication_gap_root"
  self_test_assert \
    "self-test-publication-gap-mutation-refusal" 37 \
    "an external verifier-to-publication mutation is refused before terminal publication" \
    'external-publication-gap-mutation=' \
    'candidate bytes changed in the verifier-to-publication gap' \
    'rejected staged bundle retained outside terminal proof:' \
    'proof-bundle verified:'

  local publication_collision_root="$self_test_aux_root/publication-collision"
  local publication_collision_marker="$self_test_aux_root/publication-collision-ready.txt"
  mkdir -p "$publication_collision_root"
  self_test_capture \
    "publication-destination-race-refusal" 30 1048576 \
    bash -c 'publication_destination_collision_self_test_command "$@"' _ \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    "$publication_collision_root" "$publication_collision_marker"
  self_test_verify_nested_bundle \
    'Euler-disc harness incomplete; retained evidence=' \
    "$publication_collision_root"
  self_test_assert \
    "self-test-publication-destination-race-refusal" 37 \
    "atomic no-clobber publication preserves a destination created after the early check" \
    'publication-collision-destination-preserved=true' \
    'atomic no-clobber publication refused existing destination' \
    'proof-bundle verified:'

  self_test_capture \
    "mutated-bundle" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$mutated_fixture"
  self_test_assert \
    "self-test-mutated-bundle-refusal" 1 \
    "mutated retained log invalidates its bundle" \
    'lane log size mismatch'

  self_test_capture \
    "duplicate-seal" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$duplicate_fixture"
  self_test_assert \
    "self-test-duplicate-seal-refusal" 1 \
    "a duplicate terminal seal invalidates its bundle" \
    'bundle must contain exactly one proof seal'

  self_test_capture \
    "truncated-seal" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$truncated_fixture"
  self_test_assert \
    "self-test-truncated-seal-refusal" 1 \
    "a truncated terminal seal invalidates its bundle" \
    'every verdict must be newline terminated'

  self_test_capture \
    "summary-mismatch" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$summary_mismatch_fixture"
  self_test_assert \
    "self-test-summary-mismatch-refusal" 1 \
    "summary bytes must exactly equal the terminal seal" \
    'summary.json is not byte-identical to the final proof seal'

  self_test_capture \
    "unsafe-path" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$unsafe_path_fixture"
  self_test_assert \
    "self-test-unsafe-path-refusal" 1 \
    "a traversal path cannot escape the retained bundle" \
    'unsafe lane log path'

  self_test_capture \
    "duplicate-json-key" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$duplicate_key_fixture"
  self_test_assert \
    "self-test-duplicate-json-key-refusal" 1 \
    "duplicate JSON keys are refused before semantic checks" \
    'duplicate JSON key'

  self_test_capture \
    "nonfinite-json" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$nonfinite_fixture"
  self_test_assert \
    "self-test-nonfinite-json-refusal" 1 \
    "non-finite JSON constants are refused" \
    'non-finite JSON constant'

  self_test_capture \
    "oversized-json" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$oversized_fixture"
  self_test_assert \
    "self-test-oversized-json-refusal" 1 \
    "oversized aggregate JSONL is refused before parsing" \
    'verdicts.jsonl exceeds the 2097152-byte bound'

  self_test_capture \
    "record-count" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$too_many_fixture"
  self_test_assert \
    "self-test-record-count-refusal" 1 \
    "maximum-plus-one proof records are refused" \
    'bundle exceeds the 96-record bound'

  self_test_capture \
    "readiness-mismatch" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$readiness_mismatch_fixture"
  self_test_assert \
    "self-test-readiness-mismatch-refusal" 1 \
    "READY status and readiness boolean cannot disagree" \
    'candidate_ready_for_dsr does not exactly match READY_FOR_DSR status'

  self_test_capture \
    "no-claim-mutation" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$no_claim_mutation_fixture"
  self_test_assert \
    "self-test-no-claim-mutation-refusal" 1 \
    "the proof seal cannot widen its diagnostic no-claim" \
    'summary software/physical no-claim changed'

  self_test_capture \
    "unknown-terminal" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$unknown_terminal_fixture"
  self_test_assert \
    "self-test-unknown-terminal-refusal" 1 \
    "unknown nonzero terminal states are refused" \
    'unknown terminal status'

  self_test_capture \
    "inventory-entry-cap-refusal" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$entry_cap_fixture"
  self_test_assert \
    "self-test-inventory-entry-cap-refusal" 1 \
    "maximum-plus-one bundle entries are refused before unbounded collection" \
    "bundle exceeds the ${PROOF_MAX_BUNDLE_ENTRIES}-entry inventory bound" \
    'process_group_drained=true'

  self_test_capture \
    "unexpected-file" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$unexpected_file_fixture"
  self_test_assert \
    "self-test-unreferenced-file-refusal" 1 \
    "an unreferenced file cannot enter the sealed proof inventory" \
    'bundle file inventory differs from the seal-derived exact inventory' \
    "unexpected=['unexpected.txt']"

  self_test_capture \
    "prefix-mismatch" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$prefix_mismatch_fixture"
  self_test_assert \
    "self-test-prefix-file-mismatch-refusal" 1 \
    "the retained prefix file must exactly equal the nonterminal verdict stream" \
    'verdicts-prefix.jsonl is not the exact verdict prefix'

  self_test_capture \
    "invalid-interrupted-matrix" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$invalid_interrupted_fixture"
  self_test_assert \
    "self-test-interrupted-status-matrix-refusal" 1 \
    "INTERRUPTED requires its exact terminal lane to fail" \
    'terminal lane wrapper-signal must be FAIL'

  self_test_capture \
    "invalid-incomplete-matrix" 5 1048576 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    --verify-bundle "$invalid_incomplete_fixture"
  self_test_assert \
    "self-test-incomplete-status-matrix-refusal" 1 \
    "INCOMPLETE requires its exact terminal lane to fail" \
    'terminal lane internal-incomplete must be FAIL'

  self_test_capture \
    "closure-refusal-writer" 15 1048576 \
    env FSIM_EULER_DISC_E2E_EXECUTOR=local \
    FSIM_EULER_DISC_E2E_ALLOW_LOCAL=1 \
    FSIM_EULER_DISC_E2E_CARGO=/definitely-not-used-by-closure-refusal \
    FSIM_EULER_DISC_E2E_SELF_TEST_PREFLIGHT_POLICY_REFUSAL=1 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" --profile closure
  self_test_verify_nested_bundle 'proof-bundle published: ' "$LOG_ROOT"
  self_test_assert \
    "self-test-closure-refusal-writer" 4 \
    "closure policy refusal emits the exact NO_DATA preflight/source-manifest trace" \
    '[NO_DATA] closure-root-preflight' 'command exited 1' \
    '[NO_DATA] source-manifest' 'proof-bundle verified:'

  self_test_capture \
    "closure-infrastructure-failure" 15 1048576 \
    env FSIM_EULER_DISC_E2E_EXECUTOR=local \
    FSIM_EULER_DISC_E2E_ALLOW_LOCAL=1 \
    FSIM_EULER_DISC_E2E_CARGO=/definitely-not-used-by-closure-failure \
    FSIM_EULER_DISC_E2E_SELF_TEST_PREFLIGHT_INFRA_FAILURE=1 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" --profile closure
  self_test_verify_nested_bundle 'proof-bundle published: ' "$LOG_ROOT"
  self_test_assert \
    "self-test-closure-infrastructure-failure" 7 \
    "closure infrastructure failure remains FAIL and does not mint a source-manifest NO_DATA lane" \
    '[FAIL   ] closure-root-preflight' 'command exited 2' \
    'proof-bundle verified:'

  self_test_capture \
    "aggregate-budget-exhaustion" 15 1048576 \
    env FSIM_EULER_DISC_E2E_EXECUTOR=local \
    FSIM_EULER_DISC_E2E_ALLOW_LOCAL=1 \
    FSIM_EULER_DISC_E2E_CARGO=/definitely-not-launched-by-budget-test \
    FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_AGGREGATE_BUDGET=1 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" --profile focused
  self_test_verify_nested_bundle 'proof-bundle published: ' "$LOG_ROOT"
  self_test_assert \
    "self-test-aggregate-budget-exhaustion" 7 \
    "incremental aggregate retained-log accounting refuses launch and seals a verifiable FAIL" \
    'aggregate retained-log budget exhausted before command launch' \
    '[FAIL   ] constellation-verify' 'proof-bundle verified:'

  if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_CORRUPT_CANDIDATE:-0}" == 1 ]]; then
    LANE_COMMIT_CRITICAL=1
    SELF_TEST_LAST_LOG_REL="logs/self-test-candidate-publication-parent.log"
    SELF_TEST_LAST_LOG="$LOG_DIR/$SELF_TEST_LAST_LOG_REL"
    printf '%s\n' \
      'recursive candidate-publication injection suppressed; parent self-test owns this case' \
      >"$SELF_TEST_LAST_LOG"
    record \
      "self-test-invalid-candidate-not-published" PASS "harness-self-test" \
      "parent run covers strict candidate verification before publication" \
      "$SELF_TEST_LAST_LOG_REL"
    LANE_COMMIT_CRITICAL=0
    if [[ "$WRAPPER_SIGNAL" != 0 ]]; then
      exit_for_wrapper_signal "$WRAPPER_SIGNAL"
    fi
  else
    # This child intentionally executes the complete 79-lane matrix so that the
    # injected corruption starts from a semantically valid SELF_TEST_PASS
    # candidate. Keep the recursive proof bounded while deriving enough
    # headroom for every nested normal-harness fixture on a loaded shared host.
    self_test_capture \
      "invalid-candidate-not-published" \
      "$RECURSIVE_SELF_TEST_TIMEOUT_SECONDS" 1048576 \
      env FSIM_EULER_DISC_E2E_SELF_TEST_CORRUPT_CANDIDATE=1 \
      "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" --self-test
    self_test_verify_nested_bundle \
      'Euler-disc harness incomplete; retained evidence=' "$LOG_ROOT"
    self_test_assert \
      "self-test-invalid-candidate-not-published" 1 \
      "an invalid staged success is refused and replaced by one verified incomplete terminal bundle" \
      'candidate corruption injected before verification' \
      'EXIT trap sealed an otherwise incomplete run' \
      'Euler-disc harness incomplete; retained evidence=' \
      "\"configured_lane_timeout_seconds\":$RECURSIVE_SELF_TEST_TIMEOUT_SECONDS" \
      'proof-bundle verified:'
  fi

  self_test_capture \
    "incomplete-seal" 10 1048576 \
    env FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_INCOMPLETE=1 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" --self-test
  self_test_verify_nested_bundle \
    'Euler-disc harness incomplete; retained evidence=' "$LOG_ROOT"
  self_test_assert \
    "self-test-incomplete-seal" 37 \
    "an unexpected harness exit retains one verifiable incomplete seal" \
    'EXIT trap sealed an otherwise incomplete run' \
    'Euler-disc harness incomplete; retained evidence=' \
    'proof-bundle verified:'

  self_test_capture \
    "publication-fsync-failure-retains-bundle" 10 1048576 \
    env FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_INCOMPLETE=1 \
    FSIM_EULER_DISC_E2E_SELF_TEST_FAIL_PARENT_FSYNC=1 \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" --self-test
  self_test_verify_nested_bundle \
    'Euler-disc harness incomplete; retained evidence=' "$LOG_ROOT"
  self_test_assert \
    "self-test-publication-fsync-failure-retains-bundle" 37 \
    "a post-rename parent fsync failure reports and verifies the published path" \
    'parent-directory fsync failed after atomic publication:' \
    'retraction-exchange-preflight-retained:' \
    'proof-bundle published:' \
    'Euler-disc harness incomplete; retained evidence=' \
    'proof-bundle verified:'

  self_test_capture \
    "prepublication-signal" 50 1048576 \
    bash -c 'boundary_signal_injection_self_test_command "$@"' _ \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
    "$self_test_aux_root/boundary-signal-injection"
  self_test_verify_nested_bundle \
    'Euler-disc harness interrupted; retained evidence=' \
    "$self_test_aux_root/boundary-signal-injection/final-commit" 1
  self_test_assert \
    "self-test-prepublication-signal" 143 \
    "boundary injection refuses an unusable log root and tombstone substitution while prepublication and final-commit signals reseal as interrupted" \
    'boundary-signal-mode=log-root-nofollow rc=37' \
    'log-root-nofollow-refusal=true published-success-count=0' \
    'boundary-signal-mode=tombstone-mismatch rc=125' \
    'self-test replaced public tombstone after exchange:' \
    'expected public tombstone identity was replaced after exchange' \
    'published candidate retraction failed: reason=wrapper-signal-15' \
    'retraction-integrity-concurrent-wrapper-signal=15 context=wrapper-signal-after-publication' \
    'tombstone-mismatch-published-success-count=0' \
    'boundary-signal-mode=prepublication rc=143' \
    'boundary-signal-mode=final-commit rc=143' \
    'reason=wrapper-signal-15-at-success-commit' \
    'wrapper received signal 15 and stopped scheduling' \
    'Euler-disc harness interrupted; retained evidence=' \
    'proof-bundle published:' \
    'proof-bundle verified:'

  local publication_signal_name publication_signal_number publication_signal_rc
  local publication_signal_label publication_signal_root
  while read -r publication_signal_name publication_signal_number \
      publication_signal_rc publication_signal_label; do
    publication_signal_root="$self_test_aux_root/publication-signal-${publication_signal_label}"
    self_test_capture \
      "publication-window-${publication_signal_label}" \
      "$WRAPPER_SIGNAL_SELF_TEST_TIMEOUT_SECONDS" 1048576 \
      bash -c 'wrapper_signal_self_test_command "$@"' _ \
      "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" \
      "$publication_signal_name" publication-window "$publication_signal_root"
    self_test_verify_nested_bundle \
      'Euler-disc harness interrupted; retained evidence=' \
      "$publication_signal_root/bundles"
    self_test_assert \
      "self-test-publication-window-${publication_signal_label}" \
      "$publication_signal_rc" \
      "$publication_signal_name during publication is retained as an interrupted terminal outcome" \
      'publication-window-marker=' \
      "wrapper received signal $publication_signal_number and stopped scheduling" \
      'publication finalizer signal prevented or revoked candidate admission:' \
      'proof-bundle verified:'
  done <<'EOF'
HUP 1 129 hup
INT 2 130 int
TERM 15 143 term
EOF

  self_test_capture \
    "wrapper-arbitrary-term" "$WRAPPER_SIGNAL_SELF_TEST_TIMEOUT_SECONDS" 1048576 \
    bash -c 'wrapper_signal_self_test_command "$@"' _ \
    "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" TERM ordinary-lane
  self_test_verify_nested_bundle \
    'Euler-disc harness interrupted; retained evidence=' "$LOG_ROOT"
  self_test_assert \
    "self-test-wrapper-arbitrary-term" 143 \
    "TERM during an ordinary self-test lane stops at its commit boundary without inventing that lane" \
    'interrupted command log retained outside proof bundle:' \
    'self-test-ordinary-pass.log' \
    'wrapper received signal 15 and stopped scheduling' \
    'consumed-wrapper-signal-finalizer-ready=15' \
    'consumed-wrapper-signal-finalizer-released=15' \
    'consumed-wrapper-signal-finalizer-release-acknowledged=15' \
    'no-later-lanes-launched=true' 'process_group_drained=true' \
    'absent:rejected staged bundle' 'absent:internal-incomplete' \
    'proof-bundle verified:'

  local signal_name signal_number expected_rc signal_label
  while read -r signal_name signal_number expected_rc signal_label; do
    self_test_capture \
      "wrapper-${signal_label}" "$WRAPPER_SIGNAL_SELF_TEST_TIMEOUT_SECONDS" 1048576 \
      bash -c 'wrapper_signal_self_test_command "$@"' _ \
      "$REPO_ROOT/scripts/ci/euler_disc_contract_e2e.sh" "$signal_name"
    self_test_verify_nested_bundle \
      'Euler-disc harness interrupted; retained evidence=' "$LOG_ROOT"
    self_test_assert \
      "self-test-wrapper-${signal_label}" "$expected_rc" \
      "$signal_name sent to the wrapper stops scheduling after bounded drain" \
      "wrapper received signal $signal_number and stopped scheduling" \
      "consumed-wrapper-signal-finalizer-ready=$signal_number" \
      "consumed-wrapper-signal-finalizer-released=$signal_number" \
      "consumed-wrapper-signal-finalizer-release-acknowledged=$signal_number" \
      'no-later-lanes-launched=true' 'process_group_drained=true' \
      'absent:rejected staged bundle' 'absent:internal-incomplete' \
      'proof-bundle verified:'
  done <<'EOF'
HUP 1 129 hup
INT 2 130 int
TERM 15 143 term
EOF

  if [[ "$FAILURES" != 0 ]]; then
    write_summary_and_seal \
      SELF_TEST_FAIL not-applicable false harness-self-test-no-cargo incomplete 1
    printf 'Euler-disc harness self-test FAILED; evidence=%s\n' "$LOG_DIR" >&2
    exit 1
  fi
  write_summary_and_seal \
    SELF_TEST_PASS not-applicable false harness-self-test-no-cargo stable 0
  printf 'Euler-disc harness self-test passed; retained evidence=%s\n' "$LOG_DIR" >&2
  exit 0
}

if [[ "$SELF_TEST" == 1 ]]; then
  run_harness_self_test
fi

printf 'Euler-disc contract e2e staging directory: %s\n' "$LOG_DIR" >&2
PROOF_BOUNDARY_LOG="$LOG_DIR/logs/proof-boundary.log"
LANE_COMMIT_CRITICAL=1
printf '%s\n' \
  'software/protocol structure only' \
  'no physical validation' \
  'no mechanism attribution' \
  'no emergent Euler-disc prediction' \
  "retained_log_checker_smoke_command=$RETAINED_LOG_CHECKER_SMOKE_COMMAND" \
  'retained command is not packet/case replay and resolves no artifacts' \
  >"$PROOF_BOUNDARY_LOG"
record \
  "proof-boundary" \
  PASS \
  "declaration-only" \
  "software/protocol structure only; no physical validation or emergent prediction" \
  "logs/proof-boundary.log"
LANE_COMMIT_CRITICAL=0
if [[ "$WRAPPER_SIGNAL" != 0 ]]; then
  exit_for_wrapper_signal "$WRAPPER_SIGNAL"
fi

if [[ "$PROFILE" == "closure" ]]; then
  run_refusal_lane \
    "closure-root-preflight" \
    "full-root-clean-head-observation" \
    bash -c 'closure_root_preflight_command "$@"' _ \
    "$HEAD_SHA" "${SCOPE_PATHS[@]}"
  if [[ "$LAST_LANE_RC" != 0 ]]; then
    if [[ "$LAST_LANE_STATUS" == NO_DATA ]]; then
      SOURCE_MANIFEST_LOG="$LOG_DIR/logs/source-manifest.log"
      LANE_COMMIT_CRITICAL=1
      printf '%s\n' \
        'NO_DATA: closure refused before source-manifest evaluation' \
        >"$SOURCE_MANIFEST_LOG"
      record \
        "source-manifest" \
        NO_DATA \
        "closure-refused-before-manifest-evaluation" \
        "docs/source-manifest gates skipped because closure preconditions failed" \
        "logs/source-manifest.log"
      LANE_COMMIT_CRITICAL=0
      if [[ "$WRAPPER_SIGNAL" != 0 ]]; then
        exit_for_wrapper_signal "$WRAPPER_SIGNAL"
      fi
      write_summary_and_seal \
        NO_DATA not-checked false head-bookended-closure-candidate incomplete 4
      exit 4
    fi
    write_summary_and_seal \
      FAIL not-checked false head-bookended-closure-candidate incomplete 7
    exit 7
  fi
else
  SOURCE_MANIFEST_LOG="$LOG_DIR/logs/source-manifest.log"
  LANE_COMMIT_CRITICAL=1
  printf '%s\n' \
    'NO_DATA: focused profile does not cover untracked or dirty candidate paths' \
    >"$SOURCE_MANIFEST_LOG"
  record \
    "source-manifest" \
    NO_DATA \
    "manifest-not-evaluated" \
    "focused profile does not claim untracked/dirty Euler-disc paths are covered" \
    "logs/source-manifest.log"
  LANE_COMMIT_CRITICAL=0
  if [[ "$WRAPPER_SIGNAL" != 0 ]]; then
    exit_for_wrapper_signal "$WRAPPER_SIGNAL"
  fi
fi

if [[ "${FSIM_EULER_DISC_E2E_SELF_TEST_INJECT_AGGREGATE_BUDGET:-0}" == 1 ]]; then
  RETAINED_LOG_BYTES_TOTAL="$PROOF_OPERATIONAL_LOG_BUDGET"
fi

run_lane "constellation-verify" "constellation-source-preflight" \
  scripts/ci/checkout_constellation.sh --verify-only
if [[ "$LAST_LANE_RC" != 0 ]]; then
  if [[ "$PROFILE" == "closure" ]]; then
    write_summary_and_seal \
      FAIL not-checked false head-bookended-closure-candidate incomplete 7
  else
    write_summary_and_seal \
      FAIL not-checked false focused-software-only incomplete 7
  fi
  exit 7
fi

SNAPSHOT_BEFORE="$LOG_DIR/snapshot-before.txt"
capture_constellation_snapshot "constellation-snapshot-before" "$SNAPSHOT_BEFORE"
if [[ "$LAST_LANE_RC" != 0 ]]; then
  if [[ "$PROFILE" == "closure" ]]; then
    write_summary_and_seal \
      FAIL not-checked false head-bookended-closure-candidate incomplete 7
  else
    write_summary_and_seal \
      FAIL not-checked false focused-software-only incomplete 7
  fi
  exit 7
fi

run_lane "crate-fmt" "static-hygiene-only" \
  "$CARGO_BIN" fmt -p fs-euler-disc-e2e --check
run_lane "crate-check" "focused-software-evidence" \
  "$CARGO_BIN" check --locked -p fs-euler-disc-e2e --all-targets
run_lane "retained-log-checker-smoke" "focused-software-evidence" \
  "$CARGO_BIN" test --locked -p fs-euler-disc-e2e \
  --test scientific_contract -- \
  g0_check_receipts_logs_and_domain_separation_are_exact_and_bounded \
  --exact --test-threads=1
run_lane "retained-log-checker-smoke-sentinel" "non-vacuity-evidence" \
  bash -c 'checker_smoke_sentinel_command "$@"' _ \
  "$LOG_DIR/logs/retained-log-checker-smoke.log" \
  g0_check_receipts_logs_and_domain_separation_are_exact_and_bounded
run_lane "crate-unit-integration" "focused-software-evidence" \
  "$CARGO_BIN" test --locked --no-fail-fast \
  -p fs-euler-disc-e2e --lib --tests -- --test-threads=1
run_lane "crate-doctest-hostile-boundary" "compile-fail-api-evidence" \
  "$CARGO_BIN" test --locked -p fs-euler-disc-e2e --doc
run_lane "crate-clippy" "static-hygiene-only" \
  "$CARGO_BIN" clippy --locked -p fs-euler-disc-e2e --all-targets --no-deps -- -D warnings

for gate in \
  check-layers check-deps check-contracts check-schemas check-consolidation \
  check-identities check-goldens; do
  run_lane "xtask-${gate}" "workspace-structural-gate" \
    "$CARGO_BIN" run --locked -q -p xtask -- "$gate"
done

if [[ "$PROFILE" == "closure" ]]; then
  run_lane "xtask-check-docs" "workspace-documentation-gate" \
    "$CARGO_BIN" run --locked -q -p xtask -- check-docs
  run_lane "xtask-check-source-manifest" "head-bound-source-inventory" \
    "$CARGO_BIN" run --locked -q -p xtask -- check-source-manifest

  run_lane \
    "source-manifest-membership" \
    "independent-path-membership" \
    bash -c 'source_manifest_membership_command "$@"' _ \
    "$REPO_ROOT/frankensim-source-manifest.json" "${SCOPE_PATHS[@]}"
fi

SNAPSHOT_AFTER="$LOG_DIR/snapshot-after.txt"
capture_constellation_snapshot "constellation-snapshot-after" "$SNAPSHOT_AFTER"
PROVENANCE_STATE=stable
if [[ "$LAST_LANE_RC" != 0 ]]; then
  PROVENANCE_STATE=incomplete
else
  if [[ "$PROFILE" == "closure" ]]; then
    run_lane \
      "closure-root-bookend" \
      "full-root-clean-head-observation" \
      bash -c 'closure_root_preflight_command "$@"' _ \
      "$HEAD_SHA" "${SCOPE_PATHS[@]}"
    if [[ "$LAST_LANE_RC" != 0 ]]; then
      PROVENANCE_STATE=moved
    fi
  fi
  run_lane \
    "source-stability" \
    "source-provenance-snapshot" \
    cmp -s "$SNAPSHOT_BEFORE" "$SNAPSHOT_AFTER"
  if [[ "$LAST_LANE_RC" != 0 ]]; then
    PROVENANCE_STATE=moved
  fi
fi

if [[ "$FAILURES" -ne 0 ]]; then
  if [[ "$PROFILE" == "closure" ]]; then
    write_summary_and_seal \
      FAIL unestablished false head-bookended-closure-candidate "$PROVENANCE_STATE" 1
  else
    write_summary_and_seal \
      FAIL not-checked false focused-software-only "$PROVENANCE_STATE" 1
  fi
  printf 'Euler-disc contract e2e FAILED: %s failure(s); summary=%s\n' \
    "$FAILURES" "$SUMMARY" >&2 || true
  exit 1
fi

if [[ "$PROFILE" == "closure" ]]; then
  write_summary_and_seal \
    READY_FOR_DSR full-root-clean-head-bookended true head-bookended-closure-candidate "$PROVENANCE_STATE" 0
else
  write_summary_and_seal \
    FOCUSED_PASS not-checked false focused-software-only "$PROVENANCE_STATE" 0
fi
printf 'Euler-disc contract e2e sealed successful run; summary=%s\n' \
  "$SUMMARY" >&2 || true
exit 0
