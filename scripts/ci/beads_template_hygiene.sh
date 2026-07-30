#!/usr/bin/env bash
#
# Exact-set Beads template-hygiene inventory and review harness.
#
# Bead: frankensim-semantic-bead-template-hygiene-961yr.1
#
# This script is deliberately conservative. It can inventory and classify
# template debt, emit bounded review shards, replay retained artifacts, and
# apply an explicitly reviewed SECTION_NAME_ONLY manifest. It never invents
# substantive acceptance criteria or edits tracker storage directly.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}"

exec python3 - "$@" <<'PY'
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import signal
import stat
import subprocess
import sys
import unicodedata
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path, PurePosixPath
from typing import Any, Callable, Iterable, Mapping, Sequence

try:
    import tomllib
except ModuleNotFoundError as error:  # pragma: no cover - pinned Python has tomllib
    raise SystemExit(f"Python 3.11+ with tomllib is required: {error}")


REPO_ROOT = Path.cwd().resolve()
CASE_MANIFEST_REL = PurePosixPath("tests/e2e/manifests/beads-template-hygiene.toml")
CASE_MANIFEST_V2_REL = PurePosixPath(
    "tests/e2e/manifests/beads-template-hygiene-v2.toml"
)
SCRIPT_REL = PurePosixPath("scripts/ci/beads_template_hygiene.sh")

INVENTORY_SCHEMA = "frankensim.beads-template-hygiene.inventory.v1"
PLAN_SCHEMA = "frankensim.beads-template-hygiene.plan.v1"
APPLY_SCHEMA = "frankensim.beads-template-hygiene.apply.v1"
EVENT_SCHEMA = "frankensim.beads-template-hygiene.event.v1"
SUMMARY_SCHEMA = "frankensim.beads-template-hygiene.summary.v1"
SOURCE_SCHEMA = "frankensim.beads-template-hygiene.sources.v1"
PARTITIONS_SCHEMA = "frankensim.beads-template-hygiene.partitions.v1"
RUN_TERMINAL_SCHEMA = "frankensim.beads-template-hygiene.terminal.v1"
ARTIFACT_IDENTITY_SCHEMA = "frankensim.beads-template-hygiene.artifact-identity.v1"
CASE_MANIFEST_SCHEMA = "frankensim.beads-template-hygiene-manifest.v1"

V2_MANIFEST_SCHEMA = "frankensim.beads-template-hygiene-manifest.v2"
V2_SOURCE_SCHEMA = "frankensim.beads-template-hygiene.source.v2"
V2_INVENTORY_SCHEMA = "frankensim.beads-template-hygiene.inventory.v2"
V2_AUTHORITY_SCHEMA = "frankensim.beads-template-hygiene.authority.v2"
V2_REVIEW_PLAN_SCHEMA = "frankensim.beads-template-hygiene.review-plan.v2"
V2_HISTORY_SCHEMA = "frankensim.beads-template-hygiene.history.v2"
V2_ZERO_SETS_SCHEMA = "frankensim.beads-template-hygiene.zero-sets.v2"
V2_EVENT_SCHEMA = "frankensim.beads-template-hygiene.event.v2"
V2_TERMINAL_SCHEMA = "frankensim.beads-template-hygiene.terminal.v2"
V2_REVIEW_RECEIPTS_SCHEMA = (
    "frankensim.beads-template-hygiene.review-receipts.v2"
)
V2_ASSERTION_RESULT_SCHEMA = (
    "frankensim.beads-template-hygiene.assertion-result.v2"
)
V2_TEST_TERMINAL_SCHEMA = (
    "frankensim.beads-template-hygiene.test-terminal.v2"
)

V2_RUN_ARTIFACTS = (
    "source-v2.json",
    "inventory-v2.json",
    "authority-v2.json",
    "review-plan-v2.json",
    "history-v2.json",
    "zero-sets-v2.json",
    "events.jsonl",
    "terminal.json",
    "reproduce.txt",
)
V2_READINESS_STATES = (
    "REVIEW_ONLY",
    "DECLARED_READY",
    "MECHANICALLY_APPLY_ELIGIBLE",
)
V2_REMEDIATION_ROUTES = (
    "ANALYSIS_ONLY",
    "MANUAL_BR_REVIEW",
    "AUTOMATED_CONDITIONAL",
)
V2_CLOSER_STATES = ("KNOWN", "LEGACY_UNAVAILABLE", "CONFLICTED")
V2_REVIEW_TARGET_DEFAULT = 10
V2_REVIEW_TARGET_HARD_MAX = 25
V2_REVIEW_MINUTES_CAP = 480
V2_CHILD_PAYLOAD_CAP = 262_144
V2_CHILD_DESCRIPTION_CAP = 131_072
V2_CHILD_ACCEPTANCE_CAP = 32_768
V2_CHILD_DESIGN_CAP = 16_384
V2_CHILD_NOTES_CAP = 16_384
V2_CHILD_ARGV_AGGREGATE_CAP = 65_536
V2_EXACT_OPTIMALITY_MAX_TARGETS = 13
V2_SYNOPSIS_BYTES_CAP = 8_192
V2_SYNOPSIS_ID_PREVIEW_CAP = 12
V2_AUDIT_WORKERS = 8
V2_INVENTORY_ROWS_CAP = 4_096
V2_WARNING_ROWS_CAP = 8_192
V2_WARNINGS_PER_ISSUE_CAP = 16
V2_CLAUSE_BYTES_CAP = 65_536
V2_COMMAND_ARGUMENTS_CAP = 64
V2_COMMAND_ARGUMENT_BYTES_CAP = 4_096
V2_LOG_EVENTS_CAP = 262_144
V2_LOG_LINE_BYTES_CAP = 16_384

RUN_ARTIFACTS = (
    "source.json",
    "inventory.json",
    "partitions.json",
    "plan.json",
    "events.jsonl",
    "terminal.json",
    "reproduce.txt",
)
EVENT_REQUIRED_FIELDS = (
    "schema",
    "tool",
    "rule",
    "source",
    "run",
    "case",
    "attempt",
    "stage",
    "sequence",
    "issue",
    "warning",
    "disposition",
    "shard",
    "old_semantic_root",
    "new_semantic_root",
    "command",
    "result",
    "first_divergence",
    "caps",
    "terminal",
    "inverse_br_command",
    "safe_relative_artifacts",
    "reproduction",
)
RUN_ARTIFACT_CAP = 67_108_864
EVENT_COUNT_CAP = 262_144
EVENT_LINE_CAP = 16_384

STATUS_SCOPES = ("open", "in_progress", "blocked", "deferred", "closed")
LINT_SCOPES = (*STATUS_SCOPES, "all")
STATUS_ORDER = {status: index for index, status in enumerate(STATUS_SCOPES)}
REQUIRED_SECTIONS_BY_TYPE = {
    "bug": ("## Steps to Reproduce", "## Acceptance Criteria"),
    "task": ("## Acceptance Criteria",),
    "feature": ("## Acceptance Criteria",),
    "epic": ("## Success Criteria",),
    "chore": (),
    "docs": (),
    "question": (),
}

EXPECTED_CASE_IDS = (
    "template-lint.inventory-empty",
    "template-lint.inventory-live",
    "template-lint.overlap-partitions",
    "template-lint.section-only",
    "template-lint.substantive-omission",
    "template-lint.wrong-type",
    "template-lint.rollup-gap",
    "template-lint.owner-review",
    "template-lint.historical-review",
    "template-lint.p0-shard",
    "template-lint.p1-shard",
    "template-lint.p2-p3-shard",
    "template-lint.br-only-apply",
    "template-lint.concurrent-drift",
    "template-lint.partial-batch",
    "template-lint.copied-boilerplate",
    "template-lint.scope-loss",
    "template-lint.inverse-replay",
    "template-lint.new-issue-regression",
    "template-lint.zero-debt-closeout",
    "template-lint.artifact-replay",
)

DISPOSITIONS = (
    "SECTION_NAME_ONLY",
    "SUBSTANTIVE_SEMANTIC_OMISSION",
    "MALFORMED_OR_WRONG_TYPE",
    "ROLLUP_CHILD_SET_GAP",
    "OWNER_REVIEW_REQUIRED",
    "HISTORICAL_IMMUTABLE_REVIEW",
)

TERMINAL_EXIT = {
    "Pass": 0,
    "UsageRefused": 2,
    "InputRefused": 3,
    "NoData": 4,
    "EvidenceFailed": 5,
    "CancelledDrained": 6,
    "InfrastructureFailed": 7,
    "InternalFault": 8,
}

CAPS = {
    "issues": 4_096,
    "warnings": 8_192,
    "show_batch": 100,
    "field_bytes": 2_000_000,
    "clauses_per_issue": 24,
    "clause_bytes": 512,
    "semantic_field_bytes": 65_536,
    "semantic_field_lines": 4_096,
    "relative_path_bytes": 512,
    "path_component_bytes": 128,
    "path_depth": 32,
    "plan_shard_rows": 25,
    "apply_rows": 1,
    "artifact_bytes": 256 * 1024 * 1024,
    "input_artifact_bytes": 128 * 1024 * 1024,
    "subprocess_stdout_bytes": 256 * 1024 * 1024,
    "subprocess_timeout_seconds": 180,
}

CONTRACT_TERMS = (
    "accept",
    "test",
    "verify",
    "refus",
    "replay",
    "log",
    "determin",
    "cancel",
    "bound",
    "close",
    "artifact",
    "source",
    "authority",
    "no-claim",
)

CLAUSE_TERMS = (
    "accept",
    "success",
    "test",
    "e2e",
    "log",
    "replay",
    "source",
    "authority",
    "no-claim",
    "refus",
    "cancel",
    "artifact",
    "dsr",
    "rch",
)

PATH_PATTERN = re.compile(
    r"(?<![A-Za-z0-9_.-])(?:scripts|tests|crates|docs|data)/"
    r"[A-Za-z0-9_./*+-]+"
)
HEADING_WORDS = {
    "## Acceptance Criteria": ("acceptance criteria", "acceptance"),
    "## Steps to Reproduce": ("steps to reproduce", "repro"),
    "## Success Criteria": ("success criteria", "exit criteria"),
}
SECTION_CODES = {
    "## Acceptance Criteria": "A",
    "## Steps to Reproduce": "S",
    "## Success Criteria": "C",
}
SECTION_CODE_ORDER = {"A": 0, "S": 1, "C": 2}
PARTITION_KEYS = (
    "A_only",
    "S_only",
    "C_only",
    "A_and_S_only",
    "A_and_C_only",
    "S_and_C_only",
    "A_and_S_and_C",
)
SANITIZED_ENV_NAMES = (
    "BD_DB",
    "BD_DATABASE",
    "BEADS_JSONL",
    "BR_OUTPUT_FORMAT",
    "TOON_DEFAULT_FORMAT",
    "RUST_LOG",
)
BR_READ_FLAGS = (
    "--no-auto-flush",
    "--no-auto-import",
    "--no-color",
)

_cancel_requested = False
_command_receipts: list[dict[str, Any]] = []
_v2_nomock_history_cache: dict[str, Any] | None = None


def request_cancel(_signum: int, _frame: object) -> None:
    global _cancel_requested
    _cancel_requested = True


signal.signal(signal.SIGINT, request_cancel)
signal.signal(signal.SIGTERM, request_cancel)


class HarnessError(Exception):
    terminal = "InternalFault"


class UsageRefused(HarnessError):
    terminal = "UsageRefused"


class InputRefused(HarnessError):
    terminal = "InputRefused"


class NoData(HarnessError):
    terminal = "NoData"


class EvidenceFailed(HarnessError):
    terminal = "EvidenceFailed"


class CancelledDrained(HarnessError):
    terminal = "CancelledDrained"


class InfrastructureFailed(HarnessError):
    terminal = "InfrastructureFailed"


def check_cancel() -> None:
    if _cancel_requested:
        raise CancelledDrained("cancellation requested; no mutation is in flight")


def canonical_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
        + "\n"
    ).encode("utf-8")


def strict_json_loads(
    payload: bytes | str,
    *,
    label: str,
    require_canonical: bool = False,
) -> Any:
    try:
        text = payload.decode("utf-8") if isinstance(payload, bytes) else payload
    except UnicodeDecodeError as error:
        raise InputRefused(f"{label} is not valid UTF-8") from error
    if text.startswith("\ufeff"):
        raise InputRefused(f"{label} contains a forbidden UTF-8 BOM")

    def closed_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise InputRefused(f"{label} contains duplicate JSON key {key!r}")
            result[key] = value
        return result

    def reject_constant(value: str) -> Any:
        raise InputRefused(f"{label} contains non-finite number {value}")

    try:
        document = json.loads(
            text,
            object_pairs_hook=closed_object,
            parse_constant=reject_constant,
        )
    except InputRefused:
        raise
    except json.JSONDecodeError as error:
        raise InputRefused(f"{label} is malformed JSON") from error
    if require_canonical and canonical_bytes(document) != text.encode("utf-8"):
        raise InputRefused(f"{label} is not canonical JSON with one terminal newline")
    return document


def semantic_root(value: Any) -> str:
    return "sha256-v1:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def text_root(value: str) -> str:
    return "sha256-v1:" + hashlib.sha256(value.encode("utf-8")).hexdigest()


def json_stdout(value: Any) -> None:
    sys.stdout.buffer.write(canonical_bytes(value))
    sys.stdout.flush()


def terminal_row(
    terminal: str,
    *,
    mode: str,
    detail: str,
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "schema": EVENT_SCHEMA,
        "mode": mode,
        "stage": "terminal",
        "sequence": 0,
        "terminal": terminal,
        "exit_code": TERMINAL_EXIT[terminal],
        "detail": detail,
        "no_claim": (
            "template-hygiene evidence classifies tracker-plan debt only and "
            "mints no implementation, scientific, or release authority"
        ),
    }
    if extra:
        row.update(extra)
    return row


def safe_relative(value: str, *, label: str) -> PurePosixPath:
    if not value:
        raise UsageRefused(f"{label} must be a non-empty repository-relative path")
    if value != unicodedata.normalize("NFC", value):
        raise UsageRefused(f"{label} must use NFC Unicode normalization")
    if len(value.encode("utf-8")) > CAPS["relative_path_bytes"]:
        raise UsageRefused(f"{label} exceeds the relative-path byte cap")
    if any(ord(character) < 32 or ord(character) == 127 for character in value):
        raise UsageRefused(f"{label} contains a control character")
    if "\\" in value or value.startswith("~") or re.match(r"^[A-Za-z]:", value):
        raise UsageRefused(f"{label} uses a forbidden path spelling")
    raw_parts = value.split("/")
    if (
        len(raw_parts) > CAPS["path_depth"]
        or any(part in {"", ".", ".."} for part in raw_parts)
        or any(
            len(part.encode("utf-8")) > CAPS["path_component_bytes"]
            for part in raw_parts
        )
    ):
        raise UsageRefused(f"{label} has an unsafe or over-cap component")
    candidate = PurePosixPath(value)
    if candidate.is_absolute() or ".." in candidate.parts:
        raise UsageRefused(f"{label} must be repository-relative without '..'")
    return candidate


def resolve_safe(value: str, *, label: str, must_exist: bool = False) -> Path:
    relative = safe_relative(value, label=label)
    candidate = REPO_ROOT.joinpath(*relative.parts)
    resolved = candidate.resolve(strict=False)
    try:
        resolved.relative_to(REPO_ROOT)
    except ValueError as error:
        raise UsageRefused(f"{label} escapes the repository") from error

    cursor = REPO_ROOT
    for part in relative.parts:
        cursor = cursor / part
        if cursor.exists() and cursor.is_symlink():
            raise InputRefused(f"{label} traverses symlink {cursor.relative_to(REPO_ROOT)}")
    if must_exist and not candidate.exists():
        raise InputRefused(f"{label} does not exist: {relative}")
    return candidate


def bounded_read(path: Path, *, cap: int = CAPS["input_artifact_bytes"]) -> bytes:
    if not path.is_file():
        raise InputRefused(f"expected regular file: {path.relative_to(REPO_ROOT)}")
    size = path.stat().st_size
    if size > cap:
        raise InputRefused(
            f"input artifact exceeds {cap} byte cap: {path.relative_to(REPO_ROOT)}"
        )
    return path.read_bytes()


def write_once(path: Path, payload: bytes) -> str:
    check_cancel()
    if len(payload) > CAPS["artifact_bytes"]:
        raise EvidenceFailed(
            f"artifact exceeds {CAPS['artifact_bytes']} byte cap: "
            f"{path.relative_to(REPO_ROOT)}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if not path.is_file():
            raise EvidenceFailed(
                f"artifact target is not a regular file: {path.relative_to(REPO_ROOT)}"
            )
        if bounded_read(path, cap=CAPS["artifact_bytes"]) != payload:
            raise EvidenceFailed(
                f"refusing to overwrite non-identical artifact: "
                f"{path.relative_to(REPO_ROOT)}"
            )
        return "existing-identical"
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(path, flags, 0o644)
    try:
        with os.fdopen(descriptor, "wb", closefd=True) as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        # No unlink is attempted: repository policy forbids file deletion.
        raise
    return "created"


def run_command(
    argv: Sequence[str],
    *,
    input_text: str | None = None,
    expected: Iterable[int] = (0,),
) -> subprocess.CompletedProcess[str]:
    check_cancel()
    environment = os.environ.copy()
    for name in SANITIZED_ENV_NAMES:
        environment.pop(name, None)
    try:
        completed = subprocess.run(
            list(argv),
            cwd=REPO_ROOT,
            input=input_text,
            text=True,
            capture_output=True,
            env=environment,
            timeout=CAPS["subprocess_timeout_seconds"],
            check=False,
        )
    except FileNotFoundError as error:
        raise InfrastructureFailed(f"required tool not found: {argv[0]}") from error
    except subprocess.TimeoutExpired as error:
        raise InfrastructureFailed(
            f"command exceeded {CAPS['subprocess_timeout_seconds']}s: {argv[0]}"
        ) from error
    stdout_bytes = completed.stdout.encode("utf-8")
    stderr_bytes = completed.stderr.encode("utf-8")
    if len(stdout_bytes) > CAPS["subprocess_stdout_bytes"]:
        raise InfrastructureFailed(f"{argv[0]} stdout exceeded the bounded cap")
    if len(stderr_bytes) > CAPS["subprocess_stdout_bytes"]:
        raise InfrastructureFailed(f"{argv[0]} stderr exceeded the bounded cap")
    _command_receipts.append(
        {
            "argv": [str(value) for value in argv],
            "exit_code": completed.returncode,
            "category": (
                "SUCCESS"
                if completed.returncode in set(expected)
                else "UNEXPECTED_EXIT"
            ),
            "stdout": {
                "bytes": len(stdout_bytes),
                "root": "sha256-v1:" + hashlib.sha256(stdout_bytes).hexdigest(),
                "body_retained": False,
            },
            "stderr": {
                "bytes": len(stderr_bytes),
                "root": "sha256-v1:" + hashlib.sha256(stderr_bytes).hexdigest(),
                "body_retained": False,
            },
            "redaction_verdict": "RAW_STREAM_BODIES_NOT_RETAINED",
        }
    )
    if completed.returncode not in set(expected):
        raise InfrastructureFailed(
            f"{argv[0]} returned {completed.returncode}; "
            "raw diagnostic body withheld, inspect its retained stream root"
        )
    check_cancel()
    return completed


def br_json(*arguments: str) -> Any:
    completed = run_command(("br", *arguments))
    try:
        return json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise InfrastructureFailed(
            f"br {' '.join(arguments[:2])} emitted malformed JSON"
        ) from error


def br_read_json(*arguments: str) -> Any:
    return br_json(*arguments, *BR_READ_FLAGS)


def br_version() -> dict[str, Any]:
    document = br_json("version", "--json")
    if not isinstance(document, dict):
        raise InfrastructureFailed("br version --json did not return an object")
    required = {"version", "build", "commit", "target", "features"}
    if not required.issubset(document):
        raise InfrastructureFailed("br version identity is incomplete")
    return {key: document[key] for key in sorted(required)}


def br_capabilities() -> dict[str, Any]:
    document = br_json("capabilities", "--json")
    if not isinstance(document, dict):
        raise InfrastructureFailed("br capabilities --json did not return an object")
    if document.get("contract_version") != "br.capabilities.v1":
        raise InfrastructureFailed("unknown br capabilities contract")
    return {
        "contract_version": document["contract_version"],
        "commands": document.get("commands", {}),
        "operation_count": len(document.get("operations", [])),
    }


def file_identity(relative: str) -> dict[str, Any]:
    path = REPO_ROOT / relative
    if not path.exists():
        return {"path": relative, "present": False}
    if not path.is_file():
        raise InfrastructureFailed(f"identity path is not a file: {relative}")
    digest = hashlib.sha256()
    total = 0
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            total += len(chunk)
            digest.update(chunk)
    return {
        "path": relative,
        "present": True,
        "bytes": total,
        "identity": "sha256-v1:" + digest.hexdigest(),
    }


def source_file_identities() -> list[dict[str, Any]]:
    return [
        file_identity(".beads/beads.db"),
        file_identity(".beads/beads.db-wal"),
        file_identity(".beads/issues.jsonl"),
        file_identity(str(SCRIPT_REL)),
        file_identity(str(CASE_MANIFEST_REL)),
    ]


def normalize_document(document: Any, *, label: str) -> list[dict[str, Any]]:
    if isinstance(document, list):
        issues = document
        has_more = False
    elif isinstance(document, dict):
        issues = document.get("issues")
        has_more = bool(document.get("has_more"))
    else:
        issues = None
        has_more = False
    if not isinstance(issues, list):
        raise InfrastructureFailed(f"{label} has no issues array")
    if has_more:
        raise InfrastructureFailed(f"{label} is truncated")
    if len(issues) > CAPS["issues"]:
        raise InfrastructureFailed(f"{label} exceeds {CAPS['issues']} issue cap")
    return issues


def load_case_manifest() -> dict[str, Any]:
    path = resolve_safe(str(CASE_MANIFEST_REL), label="case manifest", must_exist=True)
    try:
        document = tomllib.loads(bounded_read(path).decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise InputRefused(f"case manifest is malformed: {error}") from error
    if document.get("schema") != CASE_MANIFEST_SCHEMA:
        raise InputRefused("case manifest has an unknown schema")
    cases = document.get("case")
    if not isinstance(cases, list):
        raise InputRefused("case manifest must contain [[case]] rows")
    ids = [str(row.get("id", "")) for row in cases if isinstance(row, dict)]
    if len(ids) != len(set(ids)):
        raise InputRefused("case manifest contains duplicate case IDs")
    if tuple(ids) != EXPECTED_CASE_IDS:
        raise InputRefused(
            "case manifest membership/order differs from the frozen 21-case registry"
        )
    if tuple(str(value) for value in document.get("case_order", ())) != EXPECTED_CASE_IDS:
        raise InputRefused("case_order differs from the frozen 21-case registry")
    if int(document.get("case_count", -1)) != len(EXPECTED_CASE_IDS):
        raise InputRefused("case_count differs from the frozen registry")
    schema_contract = document.get("schema_contract")
    if not isinstance(schema_contract, dict):
        raise InputRefused("case manifest lacks schema_contract")
    required_case_fields = schema_contract.get("required_case_fields")
    if not isinstance(required_case_fields, list):
        raise InputRefused("schema_contract lacks required_case_fields")
    allowed_roles = set(str(value) for value in schema_contract.get("allowed_roles", ()))
    allowed_modes = set(str(value) for value in schema_contract.get("allowed_modes", ()))
    allowed_terminals = set(
        str(value) for value in schema_contract.get("allowed_terminals", ())
    )
    allowed_subject_terminals = set(
        str(value)
        for value in schema_contract.get("allowed_subject_terminals", ())
    )
    allowed_mutations = set(
        str(value)
        for value in schema_contract.get("allowed_mutation_bindings", ())
    )
    allowed_replays = set(
        str(value)
        for value in schema_contract.get("allowed_replay_bindings", ())
    )
    allowed_inner_categories = set(
        str(value)
        for value in schema_contract.get("allowed_inner_br_categories", ())
    )
    case_pattern = re.compile(str(schema_contract.get("case_id_pattern", r"$^")))
    for index, row in enumerate(cases):
        if not isinstance(row, dict):
            raise InputRefused(f"case row {index} is not a table")
        required = set(str(value) for value in required_case_fields)
        missing = sorted(required - set(row))
        if missing:
            raise InputRefused(f"{row.get('id', index)} lacks fields {missing}")
        extra = sorted(set(row) - required)
        if extra:
            raise InputRefused(f"{row['id']} has unknown fields {extra}")
        if int(row.get("ordinal", -1)) != index + 1:
            raise InputRefused(f"{row['id']} has a non-dense ordinal")
        if not case_pattern.fullmatch(str(row["id"])):
            raise InputRefused(f"{row['id']} has an invalid case ID")
        if str(row["role"]) not in allowed_roles:
            raise InputRefused(f"{row['id']} has an unknown role")
        if str(row["mode"]) not in allowed_modes:
            raise InputRefused(f"{row['id']} has an unknown mode")
        terminal = str(row["expected_terminal"])
        if terminal not in allowed_terminals:
            raise InputRefused(f"{row['id']} has unknown terminal {terminal}")
        if int(row["expected_exit"]) != 0:
            raise InputRefused(f"{row['id']} expectation harness must exit zero")
        if str(row["expected_subject_terminal"]) not in allowed_subject_terminals:
            raise InputRefused(f"{row['id']} has an unknown subject terminal")
        if str(row["mutation"]) not in allowed_mutations:
            raise InputRefused(f"{row['id']} has an unknown mutation binding")
        if str(row["replay"]) not in allowed_replays:
            raise InputRefused(f"{row['id']} has an unknown replay binding")
        if str(row["expected_inner_br_category"]) not in allowed_inner_categories:
            raise InputRefused(f"{row['id']} has an unknown inner br category")
        if row["expected_live_target_mutation"] is not False:
            raise InputRefused(f"{row['id']} may not mutate the live target")
        if bool(row["inverse_required"]) != (str(row["replay"]) == "inverse"):
            raise InputRefused(f"{row['id']} inverse/replay binding disagrees")
        if str(row["mode"]) == "--replay REL" and str(row["replay"]) != "artifact-only":
            raise InputRefused(f"{row['id']} --replay must be artifact-only")
    manifest_caps = document.get("caps")
    if not isinstance(manifest_caps, dict):
        raise InputRefused("case manifest lacks caps")
    if int(manifest_caps.get("max_targets_per_plan_shard", -1)) != CAPS[
        "plan_shard_rows"
    ]:
        raise InputRefused("manifest planning cap differs from harness")
    if int(manifest_caps.get("max_targets_per_apply_manifest", -1)) != CAPS[
        "apply_rows"
    ]:
        raise InputRefused("manifest apply cap differs from harness")
    document["case_count"] = len(cases)
    document["semantic_root"] = semantic_root(
        {key: value for key, value in document.items() if key != "semantic_root"}
    )
    return document


V2_MANIFEST_TOP_LEVEL_KEYS = {
    "schema",
    "schema_version",
    "manifest_id",
    "suite_id",
    "bead_id",
    "manifest_path",
    "harness_path",
    "immutable_after_acceptance",
    "criterion_count",
    "case_count",
    "assertion_count",
    "suite_status",
    "authority",
    "tracker_authority",
    "no_claim",
    "case_order",
    "schema_contract",
    "compatibility_contract",
    "source_contract",
    "row_contract",
    "authority_contract",
    "packing_contract",
    "history_contract",
    "cli_contract",
    "synopsis_contract",
    "artifact_contract",
    "logging_contract",
    "replay_contract",
    "caps",
    "case_contract",
    "criterion",
    "case",
}

V2_MANIFEST_TABLE_KEYS = {
    "schema_contract": {
        "unknown_keys",
        "duplicate_keys",
        "duplicate_case_ids",
        "duplicate_assertion_ids",
        "missing_case_ids",
        "extra_case_ids",
        "missing_assertion_ids",
        "extra_assertion_ids",
        "ordinal_policy",
        "case_order_policy",
        "assertion_order_policy",
        "criterion_policy",
        "criterion_bijection",
        "closed_inline_tables",
        "prose_assertions_authorize_pass",
        "compiled_executor_required",
        "aggregate_check_count_authorizes_pass",
        "case_id_pattern",
        "assertion_id_pattern",
        "allowed_top_level_keys",
    },
    "compatibility_contract": {
        "v1_manifest_path",
        "v1_harness_path",
        "v1_schema",
        "v1_case_count",
        "v1_assertion_count",
        "v1_case_manifest_semantic_root",
        "v1_manifest_content_root",
        "v1_harness_content_root",
        "v1_manifest_adopted_git_blob",
        "v1_harness_adopted_git_blob",
        "adoption_head",
        "adoption_source_root",
        "adoption_inventory_root",
        "adoption_issue_count",
        "adoption_warning_count",
        "adoption_artifact_dir",
        "v1_bytes_must_match_adopted_content_root",
        "v1_modes_load_v2_manifest",
        "v1_replay_loads_v2_manifest",
        "v1_manifest_may_depend_on_v2",
        "v2_records_v1_baseline_identity_only",
        "v1_live_apply_authority",
        "v1_fixture_mutation_evidence",
    },
    "source_contract": {
        "tracker_cli",
        "diagnostic_cli",
        "diagnostic_mode",
        "network_access",
        "tracker_mutation",
        "direct_tracker_file_access",
        "capture_count",
        "capture_drift",
        "count_preserving_relevant_drift",
        "all_status_partition",
        "priorities",
        "issue_types",
        "warning_classes",
        "whole_source_roots_are_work_keys",
        "unrelated_cross_run_drift_preserves_selected_keys",
        "selected_field_drift_invalidates_affected_keys",
        "selected_membership_drift_invalidates_affected_keys",
        "selected_dependency_drift_invalidates_affected_keys",
        "selected_owner_or_consumer_drift_invalidates_affected_keys",
        "selected_receipt_drift_invalidates_affected_keys",
        "agent_mail_hints_are_replay_critical",
        "campaign_epoch_root_before_lane_projection",
        "all_status_partition_receipt_required",
    },
    "row_contract": {
        "required_fields",
        "planned_child_required_fields",
        "planned_child_issue_type",
        "missing_payload_field",
        "target_fields_and_relations_remain_unchanged",
        "estimate_units",
        "review_estimate_does_not_replace_target_estimate",
    },
    "authority_contract": {
        "current_br_version",
        "readiness_states",
        "remediation_routes",
        "external_authority_verdicts",
        "conditional_write_verdicts",
        "self_report_maximum_readiness",
        "current_br_maximum_readiness",
        "current_br_automation_terminal",
        "mechanical_eligibility_requires_external_authority",
        "mechanical_eligibility_requires_atomic_conditional_write",
        "authority_alone_is_mechanical_eligibility",
        "capability_alone_is_mechanical_eligibility",
        "declared_manual_authorization_route",
        "manual_route_grants_external_authority",
        "manual_route_grants_conditional_write",
        "manual_route_claims_no_clobber",
        "manual_route_claims_exactly_once",
        "deferred_route",
        "deferred_authorization_override",
        "deferred_requires_owner_reviewed_reactivation_and_fresh_replan",
        "active_target_default_child_status",
        "active_conflict_acknowledgement_grants_semantic_authority",
        "bare_issue_type_is_grouping_authority",
    },
    "packing_contract": {
        "hard_keys",
        "compatibility_receipt_fields",
        "objective_order",
        "review_load_tuple",
        "exact_optimality_max_targets",
        "large_instance_requires_lower_bound",
        "large_instance_requires_gap_witness",
        "split_requires_rationale_and_falsifier",
        "singleton_requires_rationale_and_falsifier",
        "merge_requires_rationale_and_falsifier",
        "semantic_truncation",
        "oversize_disposition",
        "oversize_required_fields",
    },
    "history_contract": {
        "capture_command",
        "capture_count",
        "capture_drift",
        "legacy_coverage_anchor_issue",
        "legacy_coverage_anchor_closed_at",
        "legacy_coverage_anchor_status_event_id",
        "legacy_coverage_anchor_close_event_id",
        "legacy_coverage_count",
        "legacy_coverage_rows_root",
        "known_pair_max_skew_ms",
        "closer_states",
        "missing_close_actor",
        "missing_close_event",
        "duplicate_close_event",
        "conflicting_close_event",
        "malformed_timestamp",
        "missing_reason",
        "event_order_violation",
        "required_fields",
        "infer_closer_from_creator",
        "infer_closer_from_assignee",
        "infer_closer_from_self_report",
        "closed_source_mutation",
        "closed_adjudication_owner",
        "movement_classes",
        "movement_requires_before_and_after_status",
        "movement_requires_source_and_destination_roots",
        "movement_requires_immutable_prior_evidence",
        "movement_requires_exact_successor_lineage",
    },
    "cli_contract": {
        "planner_modes",
        "test_runner_modes",
        "mode_is_required",
        "multiple_modes",
        "unknown_mode",
        "unknown_option",
        "artifact_grammar",
        "artifact_root_required",
        "artifact_dir_required",
        "output_option",
        "duplicate_options",
        "review_receipts_option",
        "review_receipts_paths",
        "existing_output",
        "unsafe_or_aliased_path",
        "symlinks",
        "successful_expectation_exit",
        "usage_error_exit",
        "input_refusal_exit",
        "no_data_exit",
        "evidence_failure_exit",
        "cancelled_exit",
        "infrastructure_failure_exit",
        "internal_fault_exit",
    },
    "synopsis_contract": {
        "encoding",
        "max_bytes",
        "max_selected_ids",
        "required_counts",
        "required_work",
        "required_sections",
        "truncation_is_synopsis_only",
        "truncation_requires_total",
        "truncation_requires_shown",
        "truncation_requires_notice",
        "machine_artifacts_may_be_truncated",
        "no_tty_dependency",
        "no_color_honored",
    },
    "artifact_contract": {
        "base_set",
        "base_set_count",
        "review_plan_selected",
        "review_plan_not_requested",
        "history_plan_selected",
        "history_plan_not_requested",
        "not_requested_required_fields",
        "not_requested_state",
        "optional_content_registry",
        "oversize_inventory_required_fields",
        "unlisted_file",
        "duplicate_file",
        "missing_file",
        "changed_file",
        "unsafe_file",
        "extra_file",
        "overwrite",
        "partial_output_is_complete",
        "terminal_event_unique_last",
        "terminal_json_final_seal",
        "artifact_order",
        "successful_or_semantic_refusal_run_requires_complete_base_set",
        "output_reservation_failure_may_claim_complete_bundle",
    },
    "logging_contract": {
        "format",
        "encoding",
        "schema",
        "deterministic",
        "redacted",
        "sequence",
        "required_fields",
        "argv_representation",
        "raw_stream_bodies_retained",
        "raw_stream_roots_cover_complete_bounded_stream",
        "terminal_event_required",
        "terminal_event_must_be_last",
        "assertion_start_event_required",
        "assertion_terminal_event_required",
        "case_terminal_binds_ordered_assertion_roots",
        "suite_terminal_binds_ordered_case_roots",
        "forbidden_content",
    },
    "replay_contract": {
        "artifact_only",
        "live_tracker_reads",
        "live_tracker_writes",
        "network_access",
        "output_is_disjoint",
        "overwrite",
        "missing_unlisted_duplicate_changed_unsafe_extra",
        "selected_projection_required",
        "not_requested_projection_required",
        "all_status_partition_required",
        "oversize_inventory_required_when_nonempty",
        "first_divergence_required",
        "v2_manifest_root_required",
        "v1_manifest_root_required",
        "replay_mints_authority",
        "replay_proves_current_tracker_state",
    },
    "caps": {
        "max_inventory_rows",
        "max_warning_rows",
        "max_warnings_per_issue",
        "default_targets_per_child",
        "hard_max_targets_per_child",
        "max_review_minutes_per_child",
        "max_retained_child_payload_bytes",
        "max_description_bytes",
        "max_acceptance_bytes",
        "max_design_bytes",
        "max_notes_bytes",
        "max_aggregate_non_description_argv_bytes",
        "max_clause_bytes",
        "max_command_arguments",
        "max_command_argument_bytes",
        "max_relative_path_bytes",
        "max_path_component_bytes",
        "max_path_depth",
        "max_artifact_bytes",
        "max_log_events",
        "max_log_line_bytes",
        "max_synopsis_bytes",
        "max_synopsis_selected_ids",
        "max_cases",
        "exact_optimality_max_targets",
        "cap_exceeded_terminal",
        "cap_exceeded_live_target_mutation",
    },
    "case_contract": {
        "required_fields",
        "allowed_kinds",
        "allowed_fixture_engines",
        "expected_outer_terminal",
        "expected_outer_exit",
        "expected_live_tracker_mutation",
        "expected_refusal_requires_unchanged_projection",
        "expected_nodata_requires_unchanged_projection",
        "assertion_evidence_required",
        "assertion_skip_is_pass",
        "aggregate_count_only_is_pass",
        "compiled_executor_id",
        "runner",
    },
}


def load_case_manifest_v2() -> dict[str, Any]:
    if os.environ.get("FS_TEMPLATE_HYGIENE_FORBID_V2_LOAD") == "1":
        raise EvidenceFailed(
            "v2 manifest loading was forbidden by the v1 regression guard"
        )
    path = resolve_safe(
        str(CASE_MANIFEST_V2_REL),
        label="v2 case manifest",
        must_exist=True,
    )
    payload = bounded_read(path)
    try:
        document = tomllib.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise InputRefused(f"v2 case manifest is malformed: {error}") from error
    if not isinstance(document, dict):
        raise InputRefused("v2 case manifest root must be a table")
    v2_exact_keys(document, V2_MANIFEST_TOP_LEVEL_KEYS, label="v2 case manifest")
    if document["schema"] != V2_MANIFEST_SCHEMA or document["schema_version"] != 2:
        raise InputRefused("v2 case manifest has an unknown schema or version")
    if (
        document["manifest_path"] != str(CASE_MANIFEST_V2_REL)
        or document["harness_path"] != str(SCRIPT_REL)
        or document["immutable_after_acceptance"] is not True
        or document["tracker_authority"] != "READ_ONLY"
    ):
        raise InputRefused("v2 manifest identity or authority contract drifted")
    for table_name, expected_keys in V2_MANIFEST_TABLE_KEYS.items():
        table = document.get(table_name)
        if not isinstance(table, dict):
            raise InputRefused(f"v2 manifest lacks table {table_name}")
        v2_exact_keys(table, expected_keys, label=f"v2 manifest {table_name}")
    if set(document["schema_contract"]["allowed_top_level_keys"]) != (
        V2_MANIFEST_TOP_LEVEL_KEYS
    ):
        raise InputRefused("v2 allowed-top-level set differs from the harness")

    cases = document.get("case")
    criteria = document.get("criterion")
    if not isinstance(cases, list) or not isinstance(criteria, list):
        raise InputRefused("v2 manifest requires [[case]] and [[criterion]] rows")
    if len(cases) != 96 or document["case_count"] != 96:
        raise InputRefused("v2 manifest must contain exactly 96 cases")
    if document["assertion_count"] != len(cases):
        raise InputRefused("v2 assertion count must equal the case count")
    if len(criteria) != 25 or document["criterion_count"] != 25:
        raise InputRefused("v2 manifest must contain exactly AC01 through AC25")

    case_contract = document["case_contract"]
    required_case_fields = set(case_contract["required_fields"])
    case_pattern = re.compile(document["schema_contract"]["case_id_pattern"])
    assertion_pattern = re.compile(
        document["schema_contract"]["assertion_id_pattern"]
    )
    allowed_kinds = set(case_contract["allowed_kinds"])
    allowed_engines = set(case_contract["allowed_fixture_engines"])
    case_ids: list[str] = []
    assertion_ids: list[str] = []
    calculated_criteria: dict[str, list[str]] = defaultdict(list)
    for index, row in enumerate(cases):
        if not isinstance(row, dict):
            raise InputRefused(f"v2 case row {index} is not a table")
        v2_exact_keys(row, required_case_fields, label=f"v2 case row {index}")
        case_id = str(row["id"])
        assertion_id = str(row["assertion_id"])
        if int(row["ordinal"]) != index + 1:
            raise InputRefused(f"{case_id} has a non-dense v2 ordinal")
        if not case_pattern.fullmatch(case_id):
            raise InputRefused(f"{case_id} has an invalid v2 case ID")
        expected_assertion = "v2.case." + case_id.removeprefix(
            "template-lint-v2."
        )
        if (
            assertion_id != expected_assertion
            or not assertion_pattern.fullmatch(assertion_id)
        ):
            raise InputRefused(f"{case_id} has a non-derived assertion ID")
        row_kinds = row["kinds"]
        row_criteria = row["criterion_ids"]
        if (
            not isinstance(row_kinds, list)
            or not row_kinds
            or not set(row_kinds).issubset(allowed_kinds)
            or len(row_kinds) != len(set(row_kinds))
        ):
            raise InputRefused(f"{case_id} has malformed test kinds")
        if row["fixture_engine"] not in allowed_engines:
            raise InputRefused(f"{case_id} has an unknown fixture engine")
        if (
            not isinstance(row_criteria, list)
            or not row_criteria
            or any(not isinstance(value, str) for value in row_criteria)
            or len(row_criteria) != len(set(row_criteria))
        ):
            raise InputRefused(f"{case_id} has malformed criterion IDs")
        for criterion_id in row_criteria:
            calculated_criteria[criterion_id].append(assertion_id)
        case_ids.append(case_id)
        assertion_ids.append(assertion_id)
    v2_assert_unique(case_ids, label="v2 case IDs")
    v2_assert_unique(assertion_ids, label="v2 assertion IDs")
    if document["case_order"] != case_ids:
        raise InputRefused("v2 case_order differs from dense case rows")

    criterion_ids = [f"AC{index:02d}" for index in range(1, 26)]
    actual_criterion_ids: list[str] = []
    for index, row in enumerate(criteria):
        if not isinstance(row, dict):
            raise InputRefused(f"v2 criterion row {index} is not a table")
        v2_exact_keys(
            row,
            {"id", "assertion_ids"},
            label=f"v2 criterion row {index}",
        )
        criterion_id = str(row["id"])
        actual_criterion_ids.append(criterion_id)
        expected_assertions = calculated_criteria.get(criterion_id, [])
        if row["assertion_ids"] != expected_assertions:
            raise InputRefused(
                f"{criterion_id} assertion links differ from case declarations"
            )
    if actual_criterion_ids != criterion_ids:
        raise InputRefused("v2 criteria are not exactly ordered AC01 through AC25")
    unknown_criteria = sorted(set(calculated_criteria) - set(criterion_ids))
    if unknown_criteria:
        raise InputRefused(f"v2 cases reference unknown criteria {unknown_criteria}")

    caps = document["caps"]
    expected_caps = {
        "max_inventory_rows": V2_INVENTORY_ROWS_CAP,
        "max_warning_rows": V2_WARNING_ROWS_CAP,
        "max_warnings_per_issue": V2_WARNINGS_PER_ISSUE_CAP,
        "default_targets_per_child": V2_REVIEW_TARGET_DEFAULT,
        "hard_max_targets_per_child": V2_REVIEW_TARGET_HARD_MAX,
        "max_review_minutes_per_child": V2_REVIEW_MINUTES_CAP,
        "max_retained_child_payload_bytes": V2_CHILD_PAYLOAD_CAP,
        "max_description_bytes": V2_CHILD_DESCRIPTION_CAP,
        "max_acceptance_bytes": V2_CHILD_ACCEPTANCE_CAP,
        "max_design_bytes": V2_CHILD_DESIGN_CAP,
        "max_notes_bytes": V2_CHILD_NOTES_CAP,
        "max_aggregate_non_description_argv_bytes": (
            V2_CHILD_ARGV_AGGREGATE_CAP
        ),
        "max_clause_bytes": V2_CLAUSE_BYTES_CAP,
        "max_command_arguments": V2_COMMAND_ARGUMENTS_CAP,
        "max_command_argument_bytes": V2_COMMAND_ARGUMENT_BYTES_CAP,
        "max_relative_path_bytes": CAPS["relative_path_bytes"],
        "max_path_component_bytes": CAPS["path_component_bytes"],
        "max_path_depth": CAPS["path_depth"],
        "max_artifact_bytes": RUN_ARTIFACT_CAP,
        "max_log_events": V2_LOG_EVENTS_CAP,
        "max_log_line_bytes": V2_LOG_LINE_BYTES_CAP,
        "max_synopsis_bytes": V2_SYNOPSIS_BYTES_CAP,
        "max_synopsis_selected_ids": V2_SYNOPSIS_ID_PREVIEW_CAP,
        "max_cases": 96,
        "exact_optimality_max_targets": V2_EXACT_OPTIMALITY_MAX_TARGETS,
    }
    for name, expected in expected_caps.items():
        if caps[name] != expected:
            raise InputRefused(f"v2 manifest cap {name} differs from the harness")
    if tuple(document["artifact_contract"]["base_set"]) != V2_RUN_ARTIFACTS:
        raise InputRefused("v2 manifest base artifact set differs from the harness")
    if tuple(document["authority_contract"]["readiness_states"]) != (
        V2_READINESS_STATES
    ):
        raise InputRefused("v2 readiness states differ from the harness")
    if tuple(document["authority_contract"]["remediation_routes"]) != (
        V2_REMEDIATION_ROUTES
    ):
        raise InputRefused("v2 remediation routes differ from the harness")

    v1_path = resolve_safe(
        str(CASE_MANIFEST_REL),
        label="v1 compatibility manifest",
        must_exist=True,
    )
    v1_payload = bounded_read(v1_path)
    compatibility = document["compatibility_contract"]
    v1_content_root = "sha256-v1:" + hashlib.sha256(v1_payload).hexdigest()
    if v1_content_root != compatibility["v1_manifest_content_root"]:
        raise InputRefused("frozen v1 manifest content identity drifted")
    v1_manifest = load_case_manifest()
    if (
        v1_manifest["semantic_root"]
        != compatibility["v1_case_manifest_semantic_root"]
    ):
        raise InputRefused("frozen v1 case-manifest semantic root drifted")

    document["content_identity"] = (
        "sha256-v1:" + hashlib.sha256(payload).hexdigest()
    )
    document["semantic_root"] = semantic_root(
        {key: value for key, value in document.items() if key != "semantic_root"}
    )
    registry = v2_build_executor_registry(document)
    if set(registry) != set(assertion_ids) or any(
        not callable(executor) for executor in registry.values()
    ):
        raise InputRefused(
            "v2 compiled executor registry differs from assertion membership"
        )
    return document


def canonical_issue_projection(issue: Mapping[str, Any]) -> dict[str, Any]:
    dependencies = sorted(
        (
            {
                "id": str(edge.get("id", "")),
                "type": str(edge.get("dependency_type", "")),
                "status": str(edge.get("status", "")),
                "priority": edge.get("priority"),
            }
            for edge in (issue.get("dependencies") or [])
            if isinstance(edge, dict)
        ),
        key=lambda edge: (edge["type"], edge["id"]),
    )
    dependents = sorted(
        (
            {
                "id": str(edge.get("id", "")),
                "type": str(edge.get("dependency_type", "")),
                "status": str(edge.get("status", "")),
                "priority": edge.get("priority"),
            }
            for edge in (issue.get("dependents") or [])
            if isinstance(edge, dict)
        ),
        key=lambda edge: (edge["type"], edge["id"]),
    )
    return {
        "id": str(issue.get("id", "")),
        "title": str(issue.get("title", "")),
        "type": str(issue.get("issue_type") or issue.get("type") or ""),
        "status": str(issue.get("status", "")),
        "priority": issue.get("priority"),
        "assignee": str(issue.get("assignee") or ""),
        "owner": str(issue.get("owner") or ""),
        "parent": str(issue.get("parent") or ""),
        "labels": sorted(str(label) for label in (issue.get("labels") or [])),
        "description": str(issue.get("description") or ""),
        "acceptance_criteria": str(issue.get("acceptance_criteria") or ""),
        "design": str(issue.get("design") or ""),
        "notes": str(issue.get("notes") or ""),
        "estimated_minutes": issue.get("estimated_minutes"),
        "updated_at": str(issue.get("updated_at") or ""),
        "dependencies": dependencies,
        "dependents": dependents,
    }


def show_issues(issue_ids: Sequence[str]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    batch_size = min(CAPS["show_batch"], 56)
    for offset in range(0, len(issue_ids), batch_size):
        check_cancel()
        batch = issue_ids[offset : offset + batch_size]
        document = br_read_json("show", *batch, "--json")
        if not isinstance(document, list):
            raise InfrastructureFailed("br show batch did not return an array")
        rows.extend(
            canonical_issue_projection(issue)
            for issue in document
            if isinstance(issue, dict)
        )
    rows.sort(key=lambda issue: issue["id"])
    if [row["id"] for row in rows] != sorted(issue_ids):
        raise InfrastructureFailed("br show membership differs from lint issue set")
    return rows


def v2_full_issue_projection(issue: Mapping[str, Any]) -> dict[str, Any]:
    stable = canonical_issue_projection(issue)
    fields = {
        field: str(issue.get(field) or "")
        for field in ("description", "acceptance_criteria", "design", "notes")
    }
    return {
        **stable,
        "created_at": str(issue.get("created_at") or ""),
        "created_by": str(issue.get("created_by") or ""),
        "closed_at": str(issue.get("closed_at") or ""),
        "close_reason": str(issue.get("close_reason") or ""),
        "field_roots": {
            field: text_root(value) for field, value in sorted(fields.items())
        },
        "source_repo": str(issue.get("source_repo") or ""),
        "source_repo_path": {
            "body_retained": False,
            "root": text_root(str(issue.get("source_repo_path") or "")),
        },
        "no_claim": (
            "this complete tracker projection is source evidence only; source "
            "repository paths are root-only and no ownership or authority is inferred"
        ),
    }


def v2_show_all_issues(issue_ids: Sequence[str]) -> list[dict[str, Any]]:
    if len(issue_ids) > V2_INVENTORY_ROWS_CAP:
        raise InputRefused(
            f"v2 source has more than {V2_INVENTORY_ROWS_CAP} issue rows"
        )
    rows: list[dict[str, Any]] = []
    batch_size = min(CAPS["show_batch"], 56)
    for offset in range(0, len(issue_ids), batch_size):
        check_cancel()
        batch = issue_ids[offset : offset + batch_size]
        document = br_read_json("show", *batch, "--json")
        if not isinstance(document, list):
            raise InfrastructureFailed("v2 br show batch did not return an array")
        rows.extend(
            v2_full_issue_projection(issue)
            for issue in document
            if isinstance(issue, dict)
        )
    rows.sort(key=lambda issue: issue["id"])
    if [row["id"] for row in rows] != sorted(issue_ids):
        raise InfrastructureFailed("v2 br show membership differs from br list")
    return rows


def lint_scopes() -> dict[str, Any]:
    result: dict[str, Any] = {}
    for scope in LINT_SCOPES:
        document = br_read_json("lint", "--status", scope, "--json")
        if not isinstance(document, dict) or not isinstance(document.get("results"), list):
            raise InfrastructureFailed(f"br lint {scope} has no results array")
        total = int(document.get("total", -1))
        issues = int(document.get("issues", -1))
        if total < 0 or issues < 0:
            raise InfrastructureFailed(f"br lint {scope} has negative counts")
        if total > CAPS["warnings"]:
            raise InfrastructureFailed(f"br lint {scope} exceeds warning cap")
        result[scope] = {
            "total": total,
            "issues": issues,
            "results": sorted(
                document["results"], key=lambda row: str(row.get("id", ""))
            ),
        }
    return result


def section_present(text: str, section: str) -> bool:
    lines = visible_markdown_lines(text)
    return lines is not None and any(line.strip() == section for line in lines)


def visible_markdown_lines(text: str) -> list[str] | None:
    encoded = text.encode("utf-8")
    lines = text.splitlines()
    if (
        len(encoded) > CAPS["semantic_field_bytes"]
        or len(lines) > CAPS["semantic_field_lines"]
    ):
        return None
    visible: list[str] = []
    fenced = False
    fence_marker = ""
    in_comment = False
    for line in lines:
        stripped = line.lstrip()
        if in_comment:
            if "-->" in stripped:
                in_comment = False
            continue
        if stripped.startswith("<!--"):
            if "-->" not in stripped:
                in_comment = True
            continue
        if stripped.startswith(("```", "~~~")):
            marker = stripped[:3]
            if not fenced:
                fenced = True
                fence_marker = marker
            elif marker == fence_marker:
                fenced = False
                fence_marker = ""
            continue
        if fenced or stripped.startswith(">"):
            continue
        visible.append(line)
    return visible


def section_body(text: str, section: str) -> str | None:
    lines = visible_markdown_lines(text)
    if lines is None:
        return None
    start: int | None = None
    for index, line in enumerate(lines):
        if line.strip() == section:
            start = index + 1
            break
    if start is None:
        return None
    body: list[str] = []
    for line in lines[start:]:
        if re.match(r"^\s*#{1,6}\s+\S", line):
            break
        body.append(line)
    return "\n".join(body).strip()


def placeholder_body(text: str) -> bool:
    normalized = " ".join(text.split()).casefold().strip(" .:-")
    if not normalized:
        return True
    placeholder_terms = {
        "tbd",
        "todo",
        "to do",
        "placeholder",
        "fixme",
        "none",
        "n/a",
        "na",
        "later",
        "coming soon",
        "fill this in",
        "to be defined",
        "as above",
        "done",
        "works",
        "passes tests",
    }
    return normalized in placeholder_terms


def independent_template_findings(
    issue: Mapping[str, Any],
) -> dict[str, str]:
    issue_type = str(issue.get("type") or issue.get("issue_type") or "")
    description = str(issue.get("description") or "")
    required = REQUIRED_SECTIONS_BY_TYPE.get(issue_type, ())
    if visible_markdown_lines(description) is None:
        return {
            section: "semantic-field-cap-exceeded"
            for section in required
        }
    findings: dict[str, str] = {}
    for section in required:
        body = section_body(description, section)
        if body is None:
            findings[section] = "missing-literal-heading"
        elif placeholder_body(body):
            findings[section] = "empty-or-placeholder-body"
    return findings


def malformed_section_present(text: str, section: str) -> bool:
    visible = visible_markdown_lines(text)
    if visible is None:
        return False
    lowered_lines = [line.strip().casefold() for line in visible]
    canonical = section.casefold()
    terms = HEADING_WORDS.get(section, ())
    for line in lowered_lines:
        if line == canonical:
            continue
        stripped = line.lstrip("#*-: ").strip()
        if any(stripped == term or stripped.startswith(term + ":") for term in terms):
            return True
    return False


def relevant_clauses(issue: Mapping[str, Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for field in ("description", "acceptance_criteria", "design", "notes"):
        text = str(issue.get(field) or "")
        for line in text.splitlines():
            stripped = " ".join(line.split())
            if not stripped:
                continue
            lowered = stripped.casefold()
            if stripped.startswith("#") or any(term in lowered for term in CLAUSE_TERMS):
                encoded = stripped.encode("utf-8")
                if len(encoded) > CAPS["clause_bytes"]:
                    encoded = encoded[: CAPS["clause_bytes"]]
                    stripped = encoded.decode("utf-8", errors="ignore") + "…"
                rows.append({"field": field, "text": stripped})
                if len(rows) >= CAPS["clauses_per_issue"]:
                    return rows
    return rows


def obligation_mapping(issue: Mapping[str, Any]) -> dict[str, Any]:
    text = "\n".join(
        str(issue.get(field) or "")
        for field in ("description", "acceptance_criteria", "design", "notes")
    )
    lowered = text.casefold()
    paths = sorted(set(PATH_PATTERN.findall(text)))[:64]
    return {
        "paths": paths,
        "unit_tests": "unit test" in lowered,
        "e2e": "e2e" in lowered or "end-to-end" in lowered,
        "logging": "log" in lowered or "jsonl" in lowered or "ndjson" in lowered,
        "replay": "replay" in lowered,
        "source_closure": "source closure" in lowered or "source-closure" in lowered,
        "authority": "authority" in lowered or "authoritative" in lowered,
        "no_claim": "no-claim" in lowered or "no claim" in lowered,
        "terminal": "terminal" in lowered or "close gate" in lowered,
    }


def strong_structured_acceptance(issue: Mapping[str, Any]) -> bool:
    text = str(issue.get("acceptance_criteria") or "").strip()
    if len(text.encode("utf-8")) < 240:
        return False
    lowered = text.casefold()
    term_count = sum(term in lowered for term in CONTRACT_TERMS)
    structure_count = sum(
        1
        for line in text.splitlines()
        if re.match(r"^\s*(?:[-*]|\d+[.)])\s+", line)
        or line.strip().startswith("##")
    )
    literal_specificity = bool(PATH_PATTERN.search(text)) or bool(
        re.search(r"\b(?:P[0-4]|G[0-5]|NoData|Pass|Refused|[a-z]+-[a-z0-9.]+)\b", text)
    )
    return term_count >= 4 and structure_count >= 2 and literal_specificity


def existing_reproduction(issue: Mapping[str, Any]) -> bool:
    text = "\n".join(
        str(issue.get(field) or "")
        for field in ("description", "acceptance_criteria", "notes")
    )
    lowered = text.casefold()
    has_repro = "repro" in lowered or "reproduce" in lowered
    has_observation = any(
        marker in lowered
        for marker in ("expected", "observed", "actual", "fails", "failure")
    )
    has_literal = bool(PATH_PATTERN.search(text)) or "`" in text or "--" in text
    return has_repro and has_observation and has_literal


def classify_issue(
    issue: Mapping[str, Any],
    missing_sections: Sequence[str],
    independent_findings: Mapping[str, str] | None = None,
    *,
    duplicate_acceptance: bool = False,
    section_only_reviewed: bool = False,
) -> tuple[str, str, str]:
    status = str(issue.get("status") or "")
    issue_type = str(issue.get("type") or issue.get("issue_type") or "")
    description = str(issue.get("description") or "")
    acceptance = str(issue.get("acceptance_criteria") or "")
    missing = tuple(sorted(str(section) for section in missing_sections))
    findings = dict(independent_findings or {})

    if status == "closed":
        return (
            "HISTORICAL_IMMUTABLE_REVIEW",
            "closed history is reviewed separately and is never rewritten as live proof",
            "falsified if a closed row is placed in an active apply shard",
        )

    if status == "deferred":
        return (
            "OWNER_REVIEW_REQUIRED",
            "deferred debt remains scheduled but cannot enter an active apply shard",
            "falsified by an explicit owner transition plus a review bound to the new root",
        )

    if "empty-or-placeholder-body" in findings.values():
        return (
            "SUBSTANTIVE_SEMANTIC_OMISSION",
            "a required literal heading has an empty or explicit placeholder body",
            "falsified by issue-specific non-placeholder obligations under the heading",
        )

    if "semantic-field-cap-exceeded" in findings.values():
        return (
            "OWNER_REVIEW_REQUIRED",
            "semantic field exceeds the bounded parser cap and cannot be auto-classified",
            "falsified by a bounded owner-reviewed projection bound to the frozen field root",
        )

    if "## Steps to Reproduce" in missing and issue_type != "bug":
        return (
            "MALFORMED_OR_WRONG_TYPE",
            "a non-bug issue is being asked for a bug reproduction section",
            "falsified if the issue type is bug at the frozen source identity",
        )

    children = sorted(
        edge["id"]
        for edge in (issue.get("dependents") or [])
        if edge.get("type") == "parent-child"
    )
    rollup_text = f"{description}\n{acceptance}"
    lowered_rollup = rollup_text.casefold()
    mentions_child_contract = any(
        term in lowered_rollup
        for term in ("child", "children", "descendant", "rollup", "roll-up")
    )
    has_terminal_contract = any(
        term in lowered_rollup
        for term in ("close gate", "terminal", "done when", "completion gate")
    )
    if (
        issue_type == "epic"
        and any(
            section in missing
            for section in ("## Acceptance Criteria", "## Success Criteria")
        )
        and (
            (mentions_child_contract and not children)
            or (
                bool(children)
                and (
                    not all(child in rollup_text for child in children)
                    or not has_terminal_contract
                )
            )
        )
    ):
        return (
            "ROLLUP_CHILD_SET_GAP",
            "rollup direct-child edges, exact membership text, or terminal gate are incomplete",
            "falsified by exact direct-child edges and text plus one compatible close gate",
        )

    if (
        missing == ("## Acceptance Criteria",)
        and strong_structured_acceptance(issue)
        and not duplicate_acceptance
        and section_only_reviewed
    ):
        return (
            "SECTION_NAME_ONLY",
            "complete issue-specific criteria already exist in the structured field",
            "falsified if semantic review finds a missing obligation beyond navigation",
        )

    if (
        missing == ("## Steps to Reproduce",)
        and existing_reproduction(issue)
        and section_only_reviewed
    ):
        return (
            "SECTION_NAME_ONLY",
            "issue-specific reproduction, expected result, and observed failure already exist",
            "falsified if the existing reproduction cannot reach the named failure",
        )

    if (
        "## Acceptance Criteria" in missing
        and len(acceptance.strip().encode("utf-8")) < 120
    ):
        return (
            "SUBSTANTIVE_SEMANTIC_OMISSION",
            "no sufficiently specific structured completion contract exists",
            "falsified by locating a complete issue-specific contract in a frozen field",
        )

    if "## Steps to Reproduce" in missing and not existing_reproduction(issue):
        return (
            "SUBSTANTIVE_SEMANTIC_OMISSION",
            "the bug lacks a reachable expected-versus-observed reproduction",
            "falsified by a frozen no-mock reproduction with a named failure",
        )

    if any(malformed_section_present(description, section) for section in missing):
        return (
            "MALFORMED_OR_WRONG_TYPE",
            "complete section-like prose may exist under a malformed heading and needs typed review",
            "falsified if semantic review finds any missing obligation beyond the heading",
        )

    if issue.get("assignee") or issue.get("owner") or status == "in_progress":
        return (
            "OWNER_REVIEW_REQUIRED",
            "substantive prose exists but the active owner must adjudicate completeness",
            "falsified by owner-reviewed exact criteria or a typed omission finding",
        )

    return (
        "OWNER_REVIEW_REQUIRED",
        "the frozen fields are ambiguous and no safe automatic repair is authorized",
        "falsified by independent semantic review assigning another exact disposition",
    )


def lint_result_map(lint: Mapping[str, Any]) -> dict[str, tuple[str, ...]]:
    result: dict[str, tuple[str, ...]] = {}
    for row in lint["all"]["results"]:
        issue_id = str(row.get("id", ""))
        raw_sections = tuple(
            str(section) for section in (row.get("missing") or [])
        )
        sections = tuple(sorted(raw_sections))
        if not issue_id or not sections:
            raise InfrastructureFailed("lint result contains an empty ID or missing set")
        if len(sections) != len(set(sections)):
            raise EvidenceFailed(f"lint all duplicates a section for {issue_id}")
        if any(section not in SECTION_CODES for section in sections):
            raise EvidenceFailed(f"lint all contains an unknown section for {issue_id}")
        if issue_id in result:
            raise EvidenceFailed(f"lint all contains duplicate issue {issue_id}")
        result[issue_id] = sections
    return result


def validate_scope_partitions(lint: Mapping[str, Any]) -> None:
    all_rows = lint_result_map(lint)
    union_rows: dict[str, tuple[str, ...]] = {}
    for scope in STATUS_SCOPES:
        scope_ids: set[str] = set()
        for row in lint[scope]["results"]:
            issue_id = str(row.get("id", ""))
            sections = tuple(
                sorted(str(section) for section in (row.get("missing") or []))
            )
            if issue_id in scope_ids:
                raise EvidenceFailed(
                    f"lint scope {scope} contains duplicate issue {issue_id}"
                )
            scope_ids.add(issue_id)
            if all_rows.get(issue_id) != sections:
                raise EvidenceFailed(
                    f"lint scope {scope} is not an exact projection of all for {issue_id}"
                )
            if issue_id in union_rows:
                raise EvidenceFailed(
                    f"lint issue {issue_id} appears in multiple status cuts"
                )
            union_rows[issue_id] = sections
        warning_count = sum(
            len(row.get("missing") or []) for row in lint[scope]["results"]
        )
        if warning_count != lint[scope]["total"]:
            raise EvidenceFailed(f"lint scope {scope} warning arithmetic disagrees")
        if len(lint[scope]["results"]) != lint[scope]["issues"]:
            raise EvidenceFailed(f"lint scope {scope} issue arithmetic disagrees")
    if union_rows != all_rows:
        missing = sorted(set(all_rows) - set(union_rows))
        extra = sorted(set(union_rows) - set(all_rows))
        first = (missing or extra or ["unknown"])[0]
        raise EvidenceFailed(
            "lint all is not the exact union of status cuts; "
            f"first divergent issue {first}"
        )
    if sum(len(sections) for sections in all_rows.values()) != lint["all"]["total"]:
        raise EvidenceFailed("lint all warning arithmetic disagrees")
    if len(all_rows) != lint["all"]["issues"]:
        raise EvidenceFailed("lint all issue arithmetic disagrees")


def assemble_inventory(
    lint: Mapping[str, Any],
    issue_rows: Sequence[Mapping[str, Any]],
    source: Mapping[str, Any],
    independent_findings_by_id: Mapping[str, Mapping[str, str]] | None = None,
) -> dict[str, Any]:
    validate_scope_partitions(lint)
    br_missing_by_id = lint_result_map(lint)
    independent = {
        str(issue_id): {
            str(section): str(reason) for section, reason in findings.items()
        }
        for issue_id, findings in (independent_findings_by_id or {}).items()
        if findings
    }
    missing_by_id = {
        issue_id: tuple(
            sorted(
                set(br_missing_by_id.get(issue_id, ()))
                | set(independent.get(issue_id, ()))
            )
        )
        for issue_id in sorted(set(br_missing_by_id) | set(independent))
    }
    issue_map = {str(issue["id"]): issue for issue in issue_rows}
    if set(issue_map) != set(missing_by_id):
        raise EvidenceFailed(
            "show issue set differs from the combined br and independent issue set"
        )
    acceptance_counts = Counter(
        text_root(text)
        for issue in issue_rows
        if (text := str(issue.get("acceptance_criteria") or "").strip())
    )
    lint_status_by_id: dict[str, str] = {}
    for status in STATUS_SCOPES:
        for lint_row in lint[status]["results"]:
            issue_id = str(lint_row.get("id", ""))
            lint_status_by_id[issue_id] = status
    for issue_id, expected_status in lint_status_by_id.items():
        if issue_map[issue_id]["status"] != expected_status:
            raise EvidenceFailed(
                f"{issue_id} show status disagrees with lint cut {expected_status}"
            )

    rows: list[dict[str, Any]] = []
    warning_rows: list[dict[str, Any]] = []
    for issue_id in sorted(missing_by_id):
        issue = issue_map[issue_id]
        missing = missing_by_id[issue_id]
        finding_reasons = independent.get(issue_id, {})
        acceptance_text = str(issue.get("acceptance_criteria") or "").strip()
        duplicate_acceptance = bool(acceptance_text) and (
            acceptance_counts[text_root(acceptance_text)] > 1
        )
        disposition, rationale, falsifier = classify_issue(
            issue,
            missing,
            finding_reasons,
            duplicate_acceptance=duplicate_acceptance,
        )
        if disposition not in DISPOSITIONS:
            raise EvidenceFailed(f"unknown disposition for {issue_id}")
        fields = {
            field: text_root(str(issue.get(field) or ""))
            for field in ("description", "acceptance_criteria", "design", "notes")
        }
        domain_labels = [
            label
            for label in issue.get("labels", [])
            if label.startswith(("crate:", "area:", "layer:", "authority:"))
        ]
        mapping = obligation_mapping(issue)
        row = {
            "id": issue_id,
            "title": issue["title"],
            "type": issue["type"],
            "status": issue["status"],
            "priority": issue["priority"],
            "assignee": issue["assignee"],
            "owner": issue["owner"],
            "parent": issue["parent"],
            "authority_domain": (
                sorted(domain_labels)[0]
                if domain_labels
                else (issue["parent"] or issue["type"] or "unclassified")
            ),
            "implementation_owner": (
                issue["assignee"] or issue["owner"] or "UNASSIGNED"
            ),
            "evidence_owner": "UNRESOLVED",
            "terminal_owner": "UNRESOLVED",
            "missing_sections": list(missing),
            "br_lint_missing_sections": list(br_missing_by_id.get(issue_id, ())),
            "independent_findings": dict(sorted(finding_reasons.items())),
            "semantic_flags": {
                "duplicate_structured_acceptance": duplicate_acceptance,
                "section_only_review_bound": False,
            },
            "lint_disagreement": (
                set(br_missing_by_id.get(issue_id, ())) != set(finding_reasons)
            ),
            "overlap": len(missing) > 1,
            "disposition": disposition,
            "rationale": rationale,
            "falsifier": falsifier,
            "field_roots": fields,
            "relevant_clauses": relevant_clauses(issue),
            "dependencies": issue["dependencies"],
            "dependents": issue["dependents"],
            "mapping": mapping,
            "testing_obligation": {
                "unit_tests_named": mapping["unit_tests"],
                "e2e_named": mapping["e2e"],
            },
            "logging_obligation": mapping["logging"],
            "replay_obligation": mapping["replay"],
            "source_closure_obligation": mapping["source_closure"],
            "authority": "CLASSIFICATION_ONLY",
            "no_claim": (
                "this row schedules semantic review and does not authorize "
                "content mutation or completion"
            ),
            "estimated_minutes": issue["estimated_minutes"],
            "updated_at": issue["updated_at"],
        }
        rows.append(row)
        for section in missing:
            warning_rows.append(
                {
                    "id": issue_id,
                    "section": section,
                    "status": issue["status"],
                    "priority": issue["priority"],
                    "type": issue["type"],
                    "disposition": disposition,
                }
            )

    if len(rows) > CAPS["issues"] or len(warning_rows) > CAPS["warnings"]:
        raise EvidenceFailed("assembled inventory exceeds declared caps")

    section_counts = Counter(row["section"] for row in warning_rows)
    status_counts = Counter(row["status"] for row in rows)
    priority_counts = Counter(f"P{row['priority']}" for row in rows)
    type_counts = Counter(row["type"] for row in rows)
    disposition_counts = Counter(row["disposition"] for row in rows)
    overlap_counts = Counter(
        "+".join(row["missing_sections"]) if row["overlap"] else row["missing_sections"][0]
        for row in rows
    )
    partitions = {key: [] for key in PARTITION_KEYS}
    partition_name_by_codes = {
        ("A",): "A_only",
        ("S",): "S_only",
        ("C",): "C_only",
        ("A", "S"): "A_and_S_only",
        ("A", "C"): "A_and_C_only",
        ("S", "C"): "S_and_C_only",
        ("A", "S", "C"): "A_and_S_and_C",
    }
    for row in rows:
        try:
            codes = tuple(
                sorted(
                    (SECTION_CODES[section] for section in row["missing_sections"]),
                    key=SECTION_CODE_ORDER.__getitem__,
                )
            )
        except KeyError as error:
            raise EvidenceFailed(f"unknown missing section {error.args[0]}") from error
        partition_name = partition_name_by_codes.get(codes)
        if partition_name is None:
            raise EvidenceFailed(f"unpartitioned warning overlap for {row['id']}")
        partitions[partition_name].append(row["id"])
    partition_union = sum(len(issue_ids) for issue_ids in partitions.values())
    if partition_union != len(rows):
        raise EvidenceFailed("exclusive partition arithmetic differs from issue union")
    combined_status_cuts = {
        status: {
            "issues": sum(row["status"] == status for row in rows),
            "warnings": sum(
                len(row["missing_sections"])
                for row in rows
                if row["status"] == status
            ),
            "issue_ids": [row["id"] for row in rows if row["status"] == status],
        }
        for status in STATUS_SCOPES
    }
    warning_count_by_id = {
        row["id"]: len(row["missing_sections"]) for row in rows
    }

    inventory: dict[str, Any] = {
        "schema": INVENTORY_SCHEMA,
        "tool": {
            "br": source["br_version"],
            "rule": "br-lint-template-sections-v1",
            "classifier": "frankensim-semantic-disposition-v1",
            "case_manifest_root": source["case_manifest_root"],
        },
        "source": source,
        "caps": CAPS,
        "counts": {
            "issues": len(rows),
            "warnings": len(warning_rows),
            "br_lint_issues": len(br_missing_by_id),
            "br_lint_warnings": sum(len(value) for value in br_missing_by_id.values()),
            "independent_issues": len(independent),
            "independent_warnings": sum(len(value) for value in independent.values()),
            "lint_disagreements": sum(bool(row["lint_disagreement"]) for row in rows),
            "by_section": dict(sorted(section_counts.items())),
            "by_status": dict(sorted(status_counts.items())),
            "by_priority": dict(sorted(priority_counts.items())),
            "by_type": dict(sorted(type_counts.items())),
            "by_disposition": dict(sorted(disposition_counts.items())),
            "overlap_partitions": dict(sorted(overlap_counts.items())),
        },
        "status_cuts": {
            scope: {
                "warnings": lint[scope]["total"],
                "issues": lint[scope]["issues"],
                "issue_ids": [
                    str(row.get("id", "")) for row in lint[scope]["results"]
                ],
            }
            for scope in LINT_SCOPES
        },
        "combined_status_cuts": combined_status_cuts,
        "partitions": {
            key: {
                "issue_ids": partitions[key],
                "issues": len(partitions[key]),
                "warnings": sum(
                    warning_count_by_id[issue_id]
                    for issue_id in partitions[key]
                ),
            }
            for key in PARTITION_KEYS
        },
        "rows": rows,
        "warning_rows": warning_rows,
        "no_claim": (
            "classification schedules semantic review; it does not authorize "
            "substantive text, close target Beads, or confer product authority"
        ),
    }
    inventory["semantic_root"] = semantic_root(inventory)
    return inventory


def priority_group(priority: Any) -> str:
    if priority == 0:
        return "p0"
    if priority == 1:
        return "p1"
    if priority in (2, 3, 4):
        return "p2-p4"
    raise EvidenceFailed(f"priority {priority!r} is outside P0-P4")


def proposal_for_row(
    inventory_row: Mapping[str, Any],
    issue: Mapping[str, Any],
) -> dict[str, Any] | None:
    if inventory_row["disposition"] != "SECTION_NAME_ONLY":
        return None
    missing = tuple(inventory_row["missing_sections"])
    old = str(issue.get("description") or "")
    if missing == ("## Acceptance Criteria",):
        criteria = str(issue.get("acceptance_criteria") or "").strip()
        if not criteria:
            return None
        new = old.rstrip() + "\n\n## Acceptance Criteria\n\n" + criteria + "\n"
    elif missing == ("## Steps to Reproduce",):
        match = re.search(r"(?im)^(?P<prefix>\s*(?:REPRO|Reproduction)\s*:)", old)
        if not match:
            return None
        new = old[: match.start()] + "## Steps to Reproduce\n\n" + old[match.start() :]
    else:
        return None
    return {
        "id": inventory_row["id"],
        "disposition": "SECTION_NAME_ONLY",
        "field": "description",
        "missing_section": missing[0],
        "old_root": text_root(old),
        "new_root": text_root(new),
        "old_value": old,
        "new_value": new,
        "inverse_value": old,
        "reviewed": False,
        "rationale": inventory_row["rationale"],
    }


def build_plan(
    inventory: Mapping[str, Any],
    issue_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    issues = {str(row["id"]): row for row in issue_rows}
    groups: dict[
        tuple[int, str, int, str, str, str, str],
        list[Mapping[str, Any]],
    ] = defaultdict(list)
    for row in inventory["rows"]:
        if row["status"] == "closed":
            continue
        owner = row["assignee"] or row["owner"] or "unassigned"
        key = (
            int(row["priority"]),
            priority_group(row["priority"]),
            STATUS_ORDER[row["status"]],
            row["status"],
            owner,
            row["authority_domain"],
            row["disposition"],
        )
        groups[key].append(row)

    shards: list[dict[str, Any]] = []
    sequence = 0
    for key in sorted(groups):
        rows = sorted(groups[key], key=lambda row: row["id"])
        for offset in range(0, len(rows), CAPS["plan_shard_rows"]):
            sequence += 1
            batch = rows[offset : offset + CAPS["plan_shard_rows"]]
            shards.append(
                {
                    "sequence": sequence,
                    "priority": key[0],
                    "priority_group": key[1],
                    "status": key[3],
                    "owner": key[4],
                    "authority_domain": key[5],
                    "disposition": key[6],
                    "issue_ids": [row["id"] for row in batch],
                    "warning_count": sum(len(row["missing_sections"]) for row in batch),
                    "max_issues": CAPS["plan_shard_rows"],
                    "requires_owner_review": key[6] != "SECTION_NAME_ONLY",
                }
            )

    scheduled_ids = sorted(
        row["id"] for row in inventory["rows"] if row["status"] != "closed"
    )
    shard_ids = [issue_id for shard in shards for issue_id in shard["issue_ids"]]
    if len(shard_ids) != len(set(shard_ids)):
        raise EvidenceFailed("plan contains duplicate shard membership")
    if sorted(shard_ids) != scheduled_ids:
        raise EvidenceFailed("plan does not exact-cover every non-closed inventory row")
    inventory_by_id = {row["id"]: row for row in inventory["rows"]}
    for shard in shards:
        if len(shard["issue_ids"]) > CAPS["plan_shard_rows"]:
            raise EvidenceFailed("plan shard exceeds its target cap")
        for issue_id in shard["issue_ids"]:
            row = inventory_by_id[issue_id]
            if (
                row["priority"] != shard["priority"]
                or row["status"] != shard["status"]
                or (row["assignee"] or row["owner"] or "unassigned")
                != shard["owner"]
                or row["authority_domain"] != shard["authority_domain"]
                or row["disposition"] != shard["disposition"]
            ):
                raise EvidenceFailed(
                    f"plan shard crosses a partition boundary at {issue_id}"
                )

    proposals = [
        proposal
        for row in inventory["rows"]
        if (proposal := proposal_for_row(row, issues[row["id"]])) is not None
    ]
    plan: dict[str, Any] = {
        "schema": PLAN_SCHEMA,
        "inventory_root": inventory["semantic_root"],
        "shards": shards,
        "section_name_only_proposals": proposals,
        "apply_contract": {
            "schema": APPLY_SCHEMA,
            "max_rows": CAPS["apply_rows"],
            "requires_reviewed": True,
            "requires_reviewed_by": True,
            "requires_reservation_receipt": True,
            "mutator": "br update only",
            "rollback": "reverse br update after exact post-write root check",
        },
        "no_claim": (
            "this plan is a bounded review schedule; every substantive row remains "
            "unapplied until its owning remediation Bead is independently reviewed"
        ),
    }
    plan["semantic_root"] = semantic_root(plan)
    return plan


@dataclass(frozen=True)
class LiveSnapshot:
    lint: dict[str, Any]
    issues: list[dict[str, Any]]
    source: dict[str, Any]
    inventory: dict[str, Any]
    plan: dict[str, Any]
    all_issues: tuple[dict[str, Any], ...] = ()


V2_CAPTURE_STATE_KEYS = {
    "issue_ids",
    "issue_count",
    "issue_ids_root",
    "full_issue_projection_root",
    "lint_projection_root",
    "tracker_status_root",
    "sync_status_root",
    "export_witness_root",
    "br_version_root",
    "br_capabilities_root",
    "semantic_root",
}


def v2_capture_state(
    *,
    issue_ids: Sequence[str],
    full_issues: Sequence[Mapping[str, Any]],
    lint: Mapping[str, Any],
    tracker_status: Mapping[str, Any],
    sync_status: Mapping[str, Any],
    export_witness: Mapping[str, Any],
    version: Mapping[str, Any],
    capabilities: Mapping[str, Any],
) -> dict[str, Any]:
    normalized_ids = [str(value) for value in issue_ids]
    if not normalized_ids or normalized_ids != sorted(normalized_ids):
        raise EvidenceFailed("v2 capture issue IDs must be nonempty and sorted")
    v2_assert_unique(normalized_ids, label="v2 capture issue IDs")
    normalized_issues = [dict(value) for value in full_issues]
    full_ids = [str(value.get("id") or "") for value in normalized_issues]
    if (
        "" in full_ids
        or len(full_ids) != len(set(full_ids))
        or sorted(full_ids) != normalized_ids
    ):
        raise EvidenceFailed(
            "v2 capture full-issue membership differs from list membership"
        )
    document = {
        "issue_ids": normalized_ids,
        "issue_count": len(normalized_ids),
        "issue_ids_root": semantic_root(normalized_ids),
        "full_issue_projection_root": semantic_root(
            sorted(normalized_issues, key=lambda row: str(row["id"]))
        ),
        "lint_projection_root": semantic_root(dict(lint)),
        "tracker_status_root": semantic_root(dict(tracker_status)),
        "sync_status_root": semantic_root(dict(sync_status)),
        "export_witness_root": semantic_root(dict(export_witness)),
        "br_version_root": semantic_root(dict(version)),
        "br_capabilities_root": semantic_root(dict(capabilities)),
    }
    return v2_rooted(document)


def v2_validate_capture_pair(
    before: Mapping[str, Any],
    after: Mapping[str, Any],
) -> str:
    v2_exact_keys(before, V2_CAPTURE_STATE_KEYS, label="v2 capture before")
    v2_exact_keys(after, V2_CAPTURE_STATE_KEYS, label="v2 capture after")
    verify_semantic_root(before, label="v2 capture before")
    verify_semantic_root(after, label="v2 capture after")
    if before != after:
        changed = sorted(
            key
            for key in V2_CAPTURE_STATE_KEYS - {"semantic_root"}
            if before[key] != after[key]
        )
        raise InputRefused(
            "ConcurrentDrift: v2 coherent capture changed dimensions "
            f"{changed}"
        )
    return str(after["semantic_root"])


def collect_live(case_manifest: Mapping[str, Any]) -> LiveSnapshot:
    status = br_read_json("status", "--json", "--no-activity")
    sync_status_before = br_read_json("sync", "--status", "--json")
    export_witness_before = br_read_json("sync", "--witness", "--json")
    version = br_version()
    capabilities = br_capabilities()
    files_before = source_file_identities()
    list_before_document = br_read_json("list", "--all", "--json", "--limit", "0")
    list_before = normalize_document(list_before_document, label="br list before")
    list_projection_before = [
        canonical_issue_projection(issue) for issue in list_before
    ]
    list_root_before = semantic_root(list_projection_before)
    independent_before = {
        issue["id"]: findings
        for issue in list_projection_before
        if (findings := independent_template_findings(issue))
    }

    lint_before = lint_scopes()
    issue_ids_before = sorted(
        set(lint_result_map(lint_before)) | set(independent_before)
    )
    issues_before = show_issues(issue_ids_before)
    independent_show_before = {
        issue["id"]: findings
        for issue in issues_before
        if (findings := independent_template_findings(issue))
    }
    issue_root_before = semantic_root(issues_before)

    list_after_document = br_read_json("list", "--all", "--json", "--limit", "0")
    list_after = normalize_document(list_after_document, label="br list after")
    list_projection_after = [
        canonical_issue_projection(issue) for issue in list_after
    ]
    list_root_after = semantic_root(list_projection_after)
    independent_after = {
        issue["id"]: findings
        for issue in list_projection_after
        if (findings := independent_template_findings(issue))
    }
    lint_after = lint_scopes()
    issue_ids_after = sorted(
        set(lint_result_map(lint_after)) | set(independent_after)
    )
    issues_after = show_issues(issue_ids_after)
    independent_show_after = {
        issue["id"]: findings
        for issue in issues_after
        if (findings := independent_template_findings(issue))
    }
    issue_root_after = semantic_root(issues_after)
    sync_status_after = br_read_json("sync", "--status", "--json")
    export_witness_after = br_read_json("sync", "--witness", "--json")
    files_after = source_file_identities()

    if (
        list_root_before != list_root_after
        or issue_root_before != issue_root_after
        or issue_ids_before != issue_ids_after
        or independent_before != independent_after
        or independent_before != independent_show_before
        or independent_after != independent_show_after
        or lint_before != lint_after
        or sync_status_before != sync_status_after
        or export_witness_before != export_witness_after
        or files_before != files_after
    ):
        raise InputRefused(
            "ConcurrentDrift: Beads membership or fields changed during inventory"
        )
    if len(list_before) != len(list_after):
        raise InputRefused("ConcurrentDrift: Beads issue count changed during inventory")
    for issue in issues_after:
        if issue["status"] not in STATUS_SCOPES:
            raise EvidenceFailed(
                f"{issue['id']} has unpartitioned status {issue['status']}"
            )
        if issue["priority"] not in range(5):
            raise EvidenceFailed(
                f"{issue['id']} has unpartitioned priority {issue['priority']}"
            )

    source = {
        "schema": SOURCE_SCHEMA,
        "br_version": version,
        "br_capabilities": capabilities,
        "br_lint_rule": {
            "identity": (
                "beads-rust-lint-description-substring-v"
                f"{version['version']}"
            ),
            "field": "description",
            "matching": "case-insensitive substring; independently checked as Markdown",
            "required_sections_by_type": REQUIRED_SECTIONS_BY_TYPE,
        },
        "case_manifest_root": case_manifest["semantic_root"],
        "live_issue_count": len(list_after),
        "live_issue_projection_root": list_root_after,
        "lint_issue_projection_root": issue_root_after,
        "tracker_status": status if isinstance(status, dict) else {},
        "sync_status": sync_status_after,
        "export_witness": export_witness_after,
        "files": files_after,
        "status_cut": list(LINT_SCOPES),
    }
    source["semantic_root"] = semantic_root(source)
    inventory = assemble_inventory(
        lint_after, issues_after, source, independent_show_after
    )
    plan = build_plan(inventory, issues_after)
    return LiveSnapshot(
        lint_after,
        issues_after,
        source,
        inventory,
        plan,
        tuple(list_projection_after),
    )


def collect_live_v2(case_manifest: Mapping[str, Any]) -> LiveSnapshot:
    _command_receipts.clear()
    status_before = br_read_json("status", "--json", "--no-activity")
    sync_before = br_read_json("sync", "--status", "--json")
    witness_before = br_read_json("sync", "--witness", "--json")
    version_before = br_version()
    capabilities_before = br_capabilities()
    list_before_document = br_read_json("list", "--all", "--json", "--limit", "0")
    list_before = normalize_document(list_before_document, label="v2 br list before")
    issue_ids_before = sorted(
        str(issue.get("id") or "")
        for issue in list_before
        if isinstance(issue, dict)
    )
    if not issue_ids_before or "" in issue_ids_before:
        raise InfrastructureFailed("v2 br list contains an empty issue ID")
    v2_assert_unique(issue_ids_before, label="v2 br list issue IDs")
    full_before = v2_show_all_issues(issue_ids_before)
    lint_before = lint_scopes()
    capture_before = v2_capture_state(
        issue_ids=issue_ids_before,
        full_issues=full_before,
        lint=lint_before,
        tracker_status=status_before,
        sync_status=sync_before,
        export_witness=witness_before,
        version=version_before,
        capabilities=capabilities_before,
    )

    check_cancel()
    status_after = br_read_json("status", "--json", "--no-activity")
    sync_after = br_read_json("sync", "--status", "--json")
    witness_after = br_read_json("sync", "--witness", "--json")
    version_after = br_version()
    capabilities_after = br_capabilities()
    list_after_document = br_read_json("list", "--all", "--json", "--limit", "0")
    list_after = normalize_document(list_after_document, label="v2 br list after")
    issue_ids_after = sorted(
        str(issue.get("id") or "")
        for issue in list_after
        if isinstance(issue, dict)
    )
    full_after = v2_show_all_issues(issue_ids_after)
    lint_after = lint_scopes()
    capture_after = v2_capture_state(
        issue_ids=issue_ids_after,
        full_issues=full_after,
        lint=lint_after,
        tracker_status=status_after,
        sync_status=sync_after,
        export_witness=witness_after,
        version=version_after,
        capabilities=capabilities_after,
    )
    coherent_capture_root = v2_validate_capture_pair(
        capture_before,
        capture_after,
    )
    if len(full_after) > V2_INVENTORY_ROWS_CAP:
        raise InputRefused(
            f"v2 source exceeds {V2_INVENTORY_ROWS_CAP} issue rows"
        )
    for issue in full_after:
        if issue["status"] not in {"open", "in_progress", "blocked", "deferred", "closed"}:
            raise EvidenceFailed(
                f"{issue['id']} has unpartitioned status {issue['status']}"
            )
        if issue["priority"] not in range(5):
            raise EvidenceFailed(
                f"{issue['id']} has unpartitioned priority {issue['priority']}"
            )

    independent = {
        issue["id"]: findings
        for issue in full_after
        if (findings := independent_template_findings(issue))
    }
    lint_ids = set(lint_result_map(lint_after))
    target_ids = sorted(lint_ids | set(independent))
    full_by_id = {issue["id"]: issue for issue in full_after}
    target_issues = [full_by_id[issue_id] for issue_id in target_ids]
    source: dict[str, Any] = {
        "schema": V2_SOURCE_SCHEMA,
        "capture_contract": {
            "count": 2,
            "coherent": True,
            "tracker_cli": "br",
            "direct_tracker_file_access": False,
            "network_access": False,
        },
        "br_version": version_after,
        "br_capabilities": capabilities_after,
        "tracker_status": status_after,
        "sync_status": sync_after,
        "export_witness": witness_after,
        "case_manifest_root": case_manifest["semantic_root"],
        "case_manifest_content_identity": case_manifest["content_identity"],
        "coherent_capture_root": coherent_capture_root,
        "source_files": [
            file_identity(str(SCRIPT_REL)),
            file_identity(str(CASE_MANIFEST_REL)),
            file_identity(str(CASE_MANIFEST_V2_REL)),
        ],
        "live_issue_count": len(full_after),
        "live_issue_projection_root": semantic_root(full_after),
        "lint_issue_ids_root": semantic_root(target_ids),
        "status_cut": ["open", "in_progress", "blocked", "deferred", "closed"],
        "command_receipts": [
            {"capture_sequence": index, **receipt}
            for index, receipt in enumerate(_command_receipts)
        ],
        "no_claim": (
            "double-captured br projections are tracker-read-only source "
            "evidence; they neither authenticate people nor grant mutation authority"
        ),
    }
    source["semantic_root"] = semantic_root(source)
    inventory = assemble_inventory(
        lint_after,
        target_issues,
        source,
        independent,
    )
    plan = build_plan(inventory, target_issues)
    return LiveSnapshot(
        lint_after,
        target_issues,
        source,
        inventory,
        plan,
        tuple(full_after),
    )


def partition_projection(inventory: Mapping[str, Any]) -> dict[str, Any]:
    section_codes = {
        "## Acceptance Criteria": ("A", "MissingAcceptanceCriteria"),
        "## Steps to Reproduce": ("S", "MissingStepsToReproduce"),
        "## Success Criteria": ("C", "MissingSuccessCriteria"),
    }
    partition_names = {
        frozenset(("A",)): "A_only",
        frozenset(("S",)): "S_only",
        frozenset(("C",)): "C_only",
        frozenset(("A", "S")): "A_and_S_only",
        frozenset(("A", "C")): "A_and_C_only",
        frozenset(("S", "C")): "S_and_C_only",
        frozenset(("A", "S", "C")): "A_and_S_and_C",
    }
    memberships: dict[str, list[str]] = {"A": [], "S": [], "C": []}
    partitions: dict[str, list[str]] = {
        name: [] for name in partition_names.values()
    }
    warning_rows: list[dict[str, str]] = []
    seen_warning_rows: set[tuple[str, str, str]] = set()
    union_ids: set[str] = set()

    for row in sorted(inventory.get("rows", []), key=lambda value: value["id"]):
        issue_id = str(row["id"])
        codes: set[str] = set()
        for section in sorted(str(value) for value in row["missing_sections"]):
            if section not in section_codes:
                raise EvidenceFailed(
                    f"partition projection has unknown missing section for {issue_id}"
                )
            code, warning_class = section_codes[section]
            identity = (issue_id, warning_class, section)
            if identity in seen_warning_rows:
                raise EvidenceFailed(
                    f"partition projection has duplicate warning row for {issue_id}"
                )
            seen_warning_rows.add(identity)
            codes.add(code)
            memberships[code].append(issue_id)
            warning_rows.append(
                {
                    "issue_id": issue_id,
                    "warning_class": warning_class,
                    "missing_section": section,
                }
            )
        frozen_codes = frozenset(codes)
        if not frozen_codes or frozen_codes not in partition_names:
            raise EvidenceFailed(
                f"partition projection has an empty or unknown warning set for {issue_id}"
            )
        partitions[partition_names[frozen_codes]].append(issue_id)
        union_ids.add(issue_id)

    for values in memberships.values():
        values.sort()
    for values in partitions.values():
        values.sort()
    warning_rows.sort(
        key=lambda row: (
            row["issue_id"],
            row["warning_class"],
            row["missing_section"],
        )
    )
    partition_union = {
        issue_id for values in partitions.values() for issue_id in values
    }
    if partition_union != union_ids:
        first = sorted(partition_union ^ union_ids)[0]
        raise EvidenceFailed(
            f"partition union differs from warning membership at {first}"
        )
    if sum(len(values) for values in partitions.values()) != len(union_ids):
        raise EvidenceFailed("exclusive overlap partitions are not disjoint")
    if len(warning_rows) != int(inventory["counts"]["warnings"]):
        raise EvidenceFailed("partition warning-row arithmetic differs")

    document: dict[str, Any] = {
        "schema": PARTITIONS_SCHEMA,
        "inventory_root": inventory["semantic_root"],
        "sets": memberships,
        "partitions": partitions,
        "warning_rows": warning_rows,
        "counts": {
            "union": len(union_ids),
            "warnings": len(warning_rows),
            "sets": {
                code: len(values) for code, values in sorted(memberships.items())
            },
            "partitions": {
                name: len(values) for name, values in sorted(partitions.items())
            },
        },
        "no_claim": (
            "exact overlap arithmetic classifies membership only and authorizes "
            "no wording, mutation, closure, or product claim"
        ),
    }
    document["semantic_root"] = semantic_root(document)
    return document


def source_artifact(
    *,
    lint: Mapping[str, Any],
    issues: Sequence[Mapping[str, Any]],
    source: Mapping[str, Any],
    independent_findings_by_id: Mapping[str, Mapping[str, str]] | None = None,
) -> dict[str, Any]:
    captured = {
        "lint": lint,
        "issues": list(issues),
        "source": source,
        "independent_findings_by_id": {
            str(issue_id): dict(sorted(findings.items()))
            for issue_id, findings in sorted(
                (independent_findings_by_id or {}).items()
            )
            if findings
        },
    }
    document: dict[str, Any] = {
        "schema": SOURCE_SCHEMA,
        "case_manifest_root": source.get("case_manifest_root"),
        "captured": captured,
        "no_claim": (
            "retained source bytes support deterministic reconstruction only; "
            "they are not producer authentication or current tracker authority"
        ),
    }
    document["semantic_root"] = semantic_root(document)
    return document


def reproduction_argv(
    mode: str,
    artifact_root: str,
    artifact_dir: str,
) -> list[str]:
    if mode not in {"inventory", "plan"}:
        raise EvidenceFailed(f"unsupported evidence-bundle mode {mode!r}")
    return [
        str(SCRIPT_REL),
        f"--{mode}",
        "--artifact-root",
        str(safe_relative(artifact_root, label="artifact root")),
        "--artifact-dir",
        str(safe_relative(artifact_dir, label="artifact dir")),
    ]


def event_rows(
    mode: str,
    inventory: Mapping[str, Any],
    plan: Mapping[str, Any],
    artifact_root: str,
    artifact_dir: str,
) -> list[dict[str, Any]]:
    reproduction = reproduction_argv(mode, artifact_root, artifact_dir)
    proposal_map = {
        row["id"]: row for row in plan.get("section_name_only_proposals", [])
    }
    partitions = partition_projection(inventory)
    source_root = str(inventory["source"]["semantic_root"])
    safe_artifacts = list(RUN_ARTIFACTS)
    rows: list[dict[str, Any]] = []

    def append(
        stage: str,
        *,
        issue: str | None = None,
        warning: Any = None,
        disposition: str | None = None,
        shard: str | None = None,
        old_root: str | None = None,
        new_root: str | None = None,
        command: Sequence[str] = (),
        result: str,
        terminal: str | None = None,
        extra: Mapping[str, Any] | None = None,
    ) -> None:
        row: dict[str, Any] = {
            "schema": EVENT_SCHEMA,
            "tool": "beads-template-hygiene",
            "rule": "br-lint-template-sections-v1",
            "source": source_root,
            "run": artifact_dir,
            "mode": mode,
            "case": "template-lint.inventory-live",
            "attempt": 1,
            "stage": stage,
            "sequence": len(rows),
            "issue": issue,
            "warning": warning,
            "disposition": disposition,
            "shard": shard,
            "old_semantic_root": old_root,
            "new_semantic_root": new_root,
            "command": list(command),
            "result": result,
            "first_divergence": None,
            "caps": {
                "issues": CAPS["issues"],
                "warnings": CAPS["warnings"],
                "plan_shard_rows": CAPS["plan_shard_rows"],
                "apply_rows": CAPS["apply_rows"],
                "events": EVENT_COUNT_CAP,
                "event_line_bytes": EVENT_LINE_CAP,
                "artifact_bytes": RUN_ARTIFACT_CAP,
            },
            "terminal": terminal,
            "inverse_br_command": None,
            "safe_relative_artifacts": safe_artifacts,
            "reproduction": reproduction,
        }
        if extra:
            row.update(extra)
        missing = [field for field in EVENT_REQUIRED_FIELDS if field not in row]
        if missing:
            raise EvidenceFailed(
                f"event {stage} lacks required fields {sorted(missing)}"
            )
        rows.append(row)

    append("start", result="started")
    append(
        "source",
        result="captured",
        extra={
            "case_manifest_root": inventory["tool"]["case_manifest_root"],
            "status_cut": inventory["source"]["status_cut"],
        },
    )
    append(
        "inventory",
        result="assembled",
        extra={
            "inventory_root": inventory["semantic_root"],
            "issues": inventory["counts"]["issues"],
            "warnings": inventory["counts"]["warnings"],
        },
    )
    if mode == "plan" and inventory["rows"]:
        append(
            "join",
            result="joined",
            extra={"issues": inventory["counts"]["issues"]},
        )
    detail_stage = "classify" if mode == "plan" else "join"
    for inventory_row in sorted(inventory["rows"], key=lambda value: value["id"]):
        proposal = proposal_map.get(inventory_row["id"])
        append(
            detail_stage,
            issue=inventory_row["id"],
            warning=list(inventory_row["missing_sections"]),
            disposition=inventory_row["disposition"],
            shard=priority_group(inventory_row["priority"]),
            old_root=inventory_row["field_roots"]["description"],
            new_root=proposal["new_root"] if proposal else None,
            result="classified" if mode == "plan" else "joined",
        )
    if mode == "inventory":
        append(
            "partition",
            result="derived",
            extra={
                "partitions_root": partitions["semantic_root"],
                "partition_counts": partitions["counts"]["partitions"],
            },
        )
    else:
        append(
            "plan",
            result="bounded",
            extra={
                "plan_root": plan["semantic_root"],
                "shards": len(plan["shards"]),
                "max_targets_per_shard": CAPS["plan_shard_rows"],
            },
        )
    append(
        "publish",
        result="sealed",
        extra={
            "inventory_root": inventory["semantic_root"],
            "partitions_root": partitions["semantic_root"],
            "plan_root": plan["semantic_root"],
        },
    )
    append(
        "terminal",
        result="pass",
        terminal="Pass",
        extra={
            "exit_code": TERMINAL_EXIT["Pass"],
            "inventory_root": inventory["semantic_root"],
            "partitions_root": partitions["semantic_root"],
            "plan_root": plan["semantic_root"],
            "case_manifest_root": inventory["tool"]["case_manifest_root"],
            "event_count": len(rows) + 1,
            "no_claim": inventory["no_claim"],
        },
    )
    if len(rows) > EVENT_COUNT_CAP:
        raise EvidenceFailed(f"event stream exceeds {EVENT_COUNT_CAP} events")
    for row in rows:
        if len(canonical_bytes(row)) > EVENT_LINE_CAP:
            raise EvidenceFailed(
                f"event sequence {row['sequence']} exceeds {EVENT_LINE_CAP} bytes"
            )
    return rows


def artifact_identity(name: str, payload: bytes, schema: str) -> str:
    return semantic_root(
        {
            "schema": ARTIFACT_IDENTITY_SCHEMA,
            "content_schema": schema,
            "length": len(payload),
            "relative_path": str(safe_relative(name, label="artifact name")),
            "canonical_bytes_identity": (
                "sha256-v1:" + hashlib.sha256(payload).hexdigest()
            ),
        }
    )


def terminal_artifact(
    *,
    mode: str,
    artifact_root: str,
    artifact_dir: str,
    source: Mapping[str, Any],
    inventory: Mapping[str, Any],
    partitions: Mapping[str, Any],
    plan: Mapping[str, Any],
    event_count: int,
    payloads: Mapping[str, bytes],
) -> dict[str, Any]:
    schema_by_name = {
        "source.json": SOURCE_SCHEMA,
        "inventory.json": INVENTORY_SCHEMA,
        "partitions.json": PARTITIONS_SCHEMA,
        "plan.json": PLAN_SCHEMA,
        "events.jsonl": EVENT_SCHEMA,
        "reproduce.txt": "frankensim.argv-json.v1",
    }
    identities = {
        name: artifact_identity(name, payloads[name], schema_by_name[name])
        for name in sorted(payloads)
    }
    document: dict[str, Any] = {
        "schema": RUN_TERMINAL_SCHEMA,
        "mode": mode,
        "terminal": "Pass",
        "exit_code": TERMINAL_EXIT["Pass"],
        "artifact_root": str(safe_relative(artifact_root, label="artifact root")),
        "artifact_dir": str(safe_relative(artifact_dir, label="artifact dir")),
        "safe_relative_artifacts": list(RUN_ARTIFACTS),
        "artifact_identities": identities,
        "source_root": source["semantic_root"],
        "source_identity_root": inventory["source"]["semantic_root"],
        "case_manifest_root": inventory["tool"]["case_manifest_root"],
        "inventory_root": inventory["semantic_root"],
        "partitions_root": partitions["semantic_root"],
        "plan_root": plan["semantic_root"],
        "event_count": event_count,
        "event_sequence": list(range(event_count)),
        "reproduction": reproduction_argv(mode, artifact_root, artifact_dir),
        "issues": inventory["counts"]["issues"],
        "warnings": inventory["counts"]["warnings"],
        "no_claim": inventory["no_claim"],
    }
    document["semantic_root"] = semantic_root(document)
    return document


def build_run_payloads(
    *,
    mode: str,
    artifact_root: str,
    artifact_dir: str,
    source: Mapping[str, Any],
    inventory: Mapping[str, Any],
    partitions: Mapping[str, Any],
    plan: Mapping[str, Any],
) -> dict[str, bytes]:
    events = event_rows(
        mode,
        inventory,
        plan,
        artifact_root,
        artifact_dir,
    )
    reproduction = canonical_bytes(
        reproduction_argv(mode, artifact_root, artifact_dir)
    )
    payloads: dict[str, bytes] = {
        "source.json": canonical_bytes(source),
        "inventory.json": canonical_bytes(inventory),
        "partitions.json": canonical_bytes(partitions),
        "plan.json": canonical_bytes(plan),
        "events.jsonl": b"".join(canonical_bytes(row) for row in events),
        "reproduce.txt": reproduction,
    }
    terminal = terminal_artifact(
        mode=mode,
        artifact_root=artifact_root,
        artifact_dir=artifact_dir,
        source=source,
        inventory=inventory,
        partitions=partitions,
        plan=plan,
        event_count=len(events),
        payloads=payloads,
    )
    payloads["terminal.json"] = canonical_bytes(terminal)
    if tuple(sorted(payloads)) != tuple(sorted(RUN_ARTIFACTS)):
        raise EvidenceFailed("run payload membership differs from artifact contract")
    for name, payload in payloads.items():
        if len(payload) > RUN_ARTIFACT_CAP:
            raise EvidenceFailed(
                f"artifact {name} exceeds {RUN_ARTIFACT_CAP} byte cap"
            )
    return payloads


def resolve_run_dir(
    artifact_root: str,
    artifact_dir: str,
    *,
    label: str,
    must_exist: bool = False,
) -> Path:
    root_rel = safe_relative(artifact_root, label="artifact root")
    dir_rel = safe_relative(artifact_dir, label=label)
    combined = root_rel / dir_rel
    if len(str(combined).encode("utf-8")) > 512:
        raise UsageRefused(f"{label} exceeds the 512-byte relative path cap")
    return resolve_safe(str(combined), label=label, must_exist=must_exist)


def require_fresh_run_dir(path: Path, *, label: str) -> None:
    if path.exists():
        raise EvidenceFailed(
            f"{label} already exists; overwrite and run-directory reuse are forbidden"
        )


def write_bundle(
    *,
    mode: str,
    artifact_root: str,
    artifact_dir: str,
    snapshot: LiveSnapshot,
) -> dict[str, Any]:
    output = resolve_run_dir(
        artifact_root,
        artifact_dir,
        label="artifact dir",
    )
    independent_findings = {
        str(row["id"]): dict(row.get("independent_findings", {}))
        for row in snapshot.inventory["rows"]
        if row.get("independent_findings")
    }
    retained_source = source_artifact(
        lint=snapshot.lint,
        issues=snapshot.issues,
        source=snapshot.source,
        independent_findings_by_id=independent_findings,
    )
    partitions = partition_projection(snapshot.inventory)
    payloads = build_run_payloads(
        mode=mode,
        artifact_root=artifact_root,
        artifact_dir=artifact_dir,
        source=retained_source,
        inventory=snapshot.inventory,
        partitions=partitions,
        plan=snapshot.plan,
    )
    require_fresh_run_dir(output, label="artifact run directory")
    results = {
        name: write_once(output / name, payloads[name]) for name in sorted(payloads)
    }
    return {
        "terminal": "Pass",
        "mode": mode,
        "artifact_dir": str(
            safe_relative(artifact_root, label="artifact root")
            / safe_relative(artifact_dir, label="artifact dir")
        ),
        "source_root": retained_source["semantic_root"],
        "inventory_root": snapshot.inventory["semantic_root"],
        "partitions_root": partitions["semantic_root"],
        "plan_root": snapshot.plan["semantic_root"],
        "issues": snapshot.inventory["counts"]["issues"],
        "warnings": snapshot.inventory["counts"]["warnings"],
        "writes": results,
        "no_claim": snapshot.inventory["no_claim"],
    }


def read_json_artifact(path: Path) -> Any:
    if path.is_symlink():
        raise InputRefused(
            f"artifact replay first divergence: symlink {path.name}"
        )
    try:
        payload = bounded_read(path, cap=RUN_ARTIFACT_CAP)
        document = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise InputRefused(
            f"artifact replay first divergence: malformed {path.name}: {error}"
        ) from error
    if canonical_bytes(document) != payload:
        raise EvidenceFailed(
            f"artifact replay first divergence: non-canonical {path.name}"
        )
    return document


def verify_semantic_root(document: Mapping[str, Any], *, label: str) -> None:
    retained = document.get("semantic_root")
    expected = semantic_root(
        {key: value for key, value in document.items() if key != "semantic_root"}
    )
    if retained != expected:
        raise EvidenceFailed(
            f"artifact replay first divergence: {label}.semantic_root"
        )


def first_projection_divergence(
    expected: Any,
    observed: Any,
    *,
    path: str = "$",
) -> str | None:
    if type(expected) is not type(observed):
        return f"{path}.type"
    if isinstance(expected, dict):
        expected_keys = set(expected)
        observed_keys = set(observed)
        if expected_keys != observed_keys:
            missing = sorted(expected_keys - observed_keys)
            extra = sorted(observed_keys - expected_keys)
            key = (missing or extra)[0]
            return f"{path}.{key}"
        for key in sorted(expected):
            divergence = first_projection_divergence(
                expected[key],
                observed[key],
                path=f"{path}.{key}",
            )
            if divergence is not None:
                return divergence
        return None
    if isinstance(expected, list):
        for index, (left, right) in enumerate(zip(expected, observed)):
            divergence = first_projection_divergence(
                left,
                right,
                path=f"{path}[{index}]",
            )
            if divergence is not None:
                return divergence
        if len(expected) != len(observed):
            return f"{path}[{min(len(expected), len(observed))}]"
        return None
    return None if expected == observed else path


def require_projection_equal(label: str, expected: Any, observed: Any) -> None:
    divergence = first_projection_divergence(expected, observed)
    if divergence is not None:
        raise EvidenceFailed(
            f"artifact replay first divergence: {label}{divergence}"
        )


def require_payload_equal(label: str, expected: bytes, observed: bytes) -> None:
    for index, (left, right) in enumerate(zip(expected, observed)):
        if left != right:
            raise EvidenceFailed(
                f"artifact replay first divergence: {label} byte {index}"
            )
    if len(expected) != len(observed):
        raise EvidenceFailed(
            f"artifact replay first divergence: {label} byte "
            f"{min(len(expected), len(observed))}"
        )


def read_event_artifact(path: Path) -> tuple[list[dict[str, Any]], bytes]:
    if path.is_symlink():
        raise InputRefused("artifact replay first divergence: symlink events.jsonl")
    payload = bounded_read(path, cap=RUN_ARTIFACT_CAP)
    if not payload or not payload.endswith(b"\n"):
        raise EvidenceFailed(
            "artifact replay first divergence: events.jsonl terminal newline"
        )
    lines = payload.splitlines(keepends=True)
    if len(lines) > EVENT_COUNT_CAP:
        raise EvidenceFailed(
            "artifact replay first divergence: events.jsonl event cap"
        )
    rows: list[dict[str, Any]] = []
    for index, line in enumerate(lines):
        if len(line) > EVENT_LINE_CAP:
            raise EvidenceFailed(
                f"artifact replay first divergence: events.jsonl line {index} cap"
            )
        try:
            row = json.loads(line.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise InputRefused(
                f"artifact replay first divergence: events.jsonl line {index}: {error}"
            ) from error
        if not isinstance(row, dict):
            raise InputRefused(
                f"artifact replay first divergence: events.jsonl line {index} type"
            )
        if canonical_bytes(row) != line:
            raise EvidenceFailed(
                f"artifact replay first divergence: events.jsonl line {index} canonical"
            )
        if row.get("sequence") != index:
            raise EvidenceFailed(
                f"artifact replay first divergence: events.jsonl line {index} sequence"
            )
        missing = [field for field in EVENT_REQUIRED_FIELDS if field not in row]
        if missing:
            raise EvidenceFailed(
                f"artifact replay first divergence: events.jsonl line {index} "
                f"missing {sorted(missing)[0]}"
            )
        if not isinstance(row["command"], list):
            raise EvidenceFailed(
                f"artifact replay first divergence: events.jsonl line {index} command"
            )
        if row["inverse_br_command"] is not None and not isinstance(
            row["inverse_br_command"], list
        ):
            raise EvidenceFailed(
                f"artifact replay first divergence: events.jsonl line {index} inverse"
            )
        if not isinstance(row["reproduction"], list):
            raise EvidenceFailed(
                f"artifact replay first divergence: events.jsonl line {index} reproduction"
            )
        rows.append(row)
    if rows[-1].get("stage") != "terminal" or rows[-1].get("terminal") != "Pass":
        raise EvidenceFailed(
            "artifact replay first divergence: events.jsonl terminal event"
        )
    return rows, payload


def independent_replay_projection(sources: Mapping[str, Any]) -> dict[str, Any]:
    captured: Any = sources.get("captured", sources)
    if sources.get("captured") is not None:
        if sources.get("schema") != SOURCE_SCHEMA:
            raise InputRefused("source artifact has an unknown schema")
        verify_semantic_root(sources, label="source.json")
    if not isinstance(captured, dict):
        raise InputRefused("source artifact lacks a captured object")
    lint = captured.get("lint")
    issues = captured.get("issues")
    source = captured.get("source")
    independent = captured.get("independent_findings_by_id", {})
    if (
        not isinstance(lint, dict)
        or not isinstance(issues, list)
        or not isinstance(source, dict)
        or not isinstance(independent, dict)
    ):
        raise InputRefused(
            "source artifact lacks lint/issues/source/independent findings"
        )
    verify_semantic_root(source, label="source.json.captured.source")
    rebuilt = assemble_inventory(lint, issues, source, independent)

    warning_total = sum(len(row["missing_sections"]) for row in rebuilt["rows"])
    if warning_total != rebuilt["counts"]["warnings"]:
        raise EvidenceFailed(
            "artifact replay first divergence: inventory warning arithmetic"
        )
    if len({row["id"] for row in rebuilt["rows"]}) != rebuilt["counts"]["issues"]:
        raise EvidenceFailed(
            "artifact replay first divergence: inventory issue membership"
        )
    return rebuilt


def replay_reproduction_argv(
    artifact_root: str,
    input_dir: str,
    output_dir: str,
) -> list[str]:
    return [
        str(SCRIPT_REL),
        "--replay",
        str(safe_relative(input_dir, label="replay input")),
        "--artifact-root",
        str(safe_relative(artifact_root, label="artifact root")),
        "--artifact-dir",
        str(safe_relative(output_dir, label="artifact dir")),
    ]


def replay_event_rows(
    *,
    artifact_root: str,
    input_dir: str,
    output_dir: str,
    source: Mapping[str, Any],
    inventory: Mapping[str, Any],
    partitions: Mapping[str, Any],
    plan: Mapping[str, Any],
    retained_events_root: str,
    retained_terminal_root: str,
) -> list[dict[str, Any]]:
    reproduction = replay_reproduction_argv(
        artifact_root,
        input_dir,
        output_dir,
    )
    safe_artifacts = list(RUN_ARTIFACTS)
    rows: list[dict[str, Any]] = []

    def append(
        stage: str,
        *,
        result: str,
        command: Sequence[str] = (),
        terminal: str | None = None,
        extra: Mapping[str, Any] | None = None,
    ) -> None:
        row: dict[str, Any] = {
            "schema": EVENT_SCHEMA,
            "tool": "beads-template-hygiene",
            "rule": "artifact-only-replay-v1",
            "source": source["semantic_root"],
            "run": output_dir,
            "mode": "replay",
            "case": "template-lint.artifact-replay",
            "attempt": 1,
            "stage": stage,
            "sequence": len(rows),
            "issue": None,
            "warning": None,
            "disposition": None,
            "shard": None,
            "old_semantic_root": None,
            "new_semantic_root": None,
            "command": list(command),
            "result": result,
            "first_divergence": None,
            "caps": {
                "events": EVENT_COUNT_CAP,
                "event_line_bytes": EVENT_LINE_CAP,
                "artifact_bytes": RUN_ARTIFACT_CAP,
            },
            "terminal": terminal,
            "inverse_br_command": None,
            "safe_relative_artifacts": safe_artifacts,
            "reproduction": reproduction,
        }
        if extra:
            row.update(extra)
        missing = [field for field in EVENT_REQUIRED_FIELDS if field not in row]
        if missing:
            raise EvidenceFailed(
                f"replay event {stage} lacks required fields {sorted(missing)}"
            )
        rows.append(row)

    append("start", result="started", command=reproduction)
    append(
        "artifact-admission",
        result="admitted",
        extra={
            "input_dir": input_dir,
            "required_artifacts": safe_artifacts,
            "live_br_access": "forbidden-and-not-used",
        },
    )
    append(
        "reconstruct",
        result="reconstructed",
        extra={
            "source_root": source["semantic_root"],
            "inventory_root": inventory["semantic_root"],
            "partitions_root": partitions["semantic_root"],
            "plan_root": plan["semantic_root"],
        },
    )
    append(
        "compare",
        result="equivalent",
        extra={
            "retained_events_root": retained_events_root,
            "reconstructed_events_root": retained_events_root,
            "retained_terminal_root": retained_terminal_root,
            "reconstructed_terminal_root": retained_terminal_root,
        },
    )
    append(
        "publish-disjoint",
        result="ready",
        extra={"input_dir": input_dir, "output_dir": output_dir},
    )
    append(
        "terminal",
        result="pass",
        terminal="Pass",
        extra={
            "exit_code": TERMINAL_EXIT["Pass"],
            "subject_terminal": "REPLAY_EQUIVALENT",
            "event_count": len(rows) + 1,
            "inventory_root": inventory["semantic_root"],
            "partitions_root": partitions["semantic_root"],
            "plan_root": plan["semantic_root"],
            "first_divergence": None,
            "no_claim": (
                "artifact-only replay proves deterministic reconstruction of "
                "the retained projection and mints no authority"
            ),
        },
    )
    for row in rows:
        if len(canonical_bytes(row)) > EVENT_LINE_CAP:
            raise EvidenceFailed(
                f"replay event sequence {row['sequence']} exceeds "
                f"{EVENT_LINE_CAP} bytes"
            )
    return rows


def build_replay_output_payloads(
    *,
    artifact_root: str,
    input_dir: str,
    output_dir: str,
    source: Mapping[str, Any],
    inventory: Mapping[str, Any],
    partitions: Mapping[str, Any],
    plan: Mapping[str, Any],
    retained_events_root: str,
    retained_terminal_root: str,
    retained_artifact_identities: Mapping[str, Any],
) -> dict[str, bytes]:
    events = replay_event_rows(
        artifact_root=artifact_root,
        input_dir=input_dir,
        output_dir=output_dir,
        source=source,
        inventory=inventory,
        partitions=partitions,
        plan=plan,
        retained_events_root=retained_events_root,
        retained_terminal_root=retained_terminal_root,
    )
    reproduction = replay_reproduction_argv(
        artifact_root,
        input_dir,
        output_dir,
    )
    payloads: dict[str, bytes] = {
        "source.json": canonical_bytes(source),
        "inventory.json": canonical_bytes(inventory),
        "partitions.json": canonical_bytes(partitions),
        "plan.json": canonical_bytes(plan),
        "events.jsonl": b"".join(canonical_bytes(row) for row in events),
        "reproduce.txt": canonical_bytes(reproduction),
    }
    schema_by_name = {
        "source.json": SOURCE_SCHEMA,
        "inventory.json": INVENTORY_SCHEMA,
        "partitions.json": PARTITIONS_SCHEMA,
        "plan.json": PLAN_SCHEMA,
        "events.jsonl": EVENT_SCHEMA,
        "reproduce.txt": "frankensim.argv-json.v1",
    }
    output_identities = {
        name: artifact_identity(name, payloads[name], schema_by_name[name])
        for name in sorted(payloads)
    }
    terminal: dict[str, Any] = {
        "schema": RUN_TERMINAL_SCHEMA,
        "mode": "replay",
        "terminal": "Pass",
        "subject_terminal": "REPLAY_EQUIVALENT",
        "exit_code": TERMINAL_EXIT["Pass"],
        "artifact_root": str(safe_relative(artifact_root, label="artifact root")),
        "input_dir": str(safe_relative(input_dir, label="replay input")),
        "artifact_dir": str(safe_relative(output_dir, label="artifact dir")),
        "safe_relative_artifacts": list(RUN_ARTIFACTS),
        "artifact_identities": output_identities,
        "retained_artifact_identities": dict(
            sorted(retained_artifact_identities.items())
        ),
        "source_root": source["semantic_root"],
        "inventory_root": inventory["semantic_root"],
        "partitions_root": partitions["semantic_root"],
        "plan_root": plan["semantic_root"],
        "retained_events_root": retained_events_root,
        "reconstructed_events_root": retained_events_root,
        "retained_terminal_root": retained_terminal_root,
        "reconstructed_terminal_root": retained_terminal_root,
        "event_count": len(events),
        "event_sequence": list(range(len(events))),
        "first_divergence": None,
        "live_br_access": "forbidden-and-not-used",
        "reproduction": reproduction,
        "no_claim": (
            "replay equivalence is deterministic reconstruction evidence, not "
            "authentication, current-state proof, semantic approval, or "
            "capability authority"
        ),
    }
    terminal["semantic_root"] = semantic_root(terminal)
    payloads["terminal.json"] = canonical_bytes(terminal)
    for name, payload in payloads.items():
        if len(payload) > RUN_ARTIFACT_CAP:
            raise EvidenceFailed(
                f"replay artifact {name} exceeds {RUN_ARTIFACT_CAP} byte cap"
            )
    return payloads


def replay_bundle(
    *,
    artifact_root: str,
    input_dir: str,
    output_dir: str,
    case_manifest_root: str | None = None,
) -> dict[str, Any]:
    source_dir = resolve_run_dir(
        artifact_root,
        input_dir,
        label="replay input",
        must_exist=True,
    )
    target_dir = resolve_run_dir(
        artifact_root,
        output_dir,
        label="artifact dir",
    )
    if not source_dir.is_dir():
        raise InputRefused("replay input is not a directory")
    source_resolved = source_dir.resolve(strict=True)
    target_resolved = target_dir.resolve(strict=False)
    if source_resolved == target_resolved:
        raise InputRefused("replay input and output must be disjoint")
    if (
        source_resolved in target_resolved.parents
        or target_resolved in source_resolved.parents
    ):
        raise InputRefused("replay input/output may not contain one another")

    for name in RUN_ARTIFACTS:
        candidate = source_dir / name
        if not candidate.exists():
            raise InputRefused(
                f"artifact replay first divergence: missing {name}"
            )
        if candidate.is_symlink() or not candidate.is_file():
            raise InputRefused(
                f"artifact replay first divergence: unsafe {name}"
            )

    retained_source = read_json_artifact(source_dir / "source.json")
    retained_inventory = read_json_artifact(source_dir / "inventory.json")
    retained_partitions = read_json_artifact(source_dir / "partitions.json")
    retained_plan = read_json_artifact(source_dir / "plan.json")
    retained_terminal = read_json_artifact(source_dir / "terminal.json")
    retained_events, retained_event_bytes = read_event_artifact(
        source_dir / "events.jsonl"
    )
    reproduce_path = source_dir / "reproduce.txt"
    if reproduce_path.is_symlink():
        raise InputRefused(
            "artifact replay first divergence: symlink reproduce.txt"
        )
    retained_reproduce = bounded_read(reproduce_path, cap=RUN_ARTIFACT_CAP)
    try:
        reproduction_document = json.loads(retained_reproduce.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise InputRefused(
            f"artifact replay first divergence: malformed reproduce.txt: {error}"
        ) from error
    if (
        not isinstance(reproduction_document, list)
        or not all(isinstance(value, str) for value in reproduction_document)
        or canonical_bytes(reproduction_document) != retained_reproduce
    ):
        raise EvidenceFailed(
            "artifact replay first divergence: reproduce.txt canonical argv"
        )

    for label, document, schema in (
        ("inventory.json", retained_inventory, INVENTORY_SCHEMA),
        ("partitions.json", retained_partitions, PARTITIONS_SCHEMA),
        ("plan.json", retained_plan, PLAN_SCHEMA),
        ("terminal.json", retained_terminal, RUN_TERMINAL_SCHEMA),
    ):
        if not isinstance(document, dict) or document.get("schema") != schema:
            raise InputRefused(
                f"artifact replay first divergence: {label} schema"
            )
        verify_semantic_root(document, label=label)
    if not isinstance(retained_source, dict):
        raise InputRefused("artifact replay first divergence: source.json type")
    if case_manifest_root is not None and (
        retained_source.get("case_manifest_root") != case_manifest_root
        or retained_terminal.get("case_manifest_root") != case_manifest_root
    ):
        raise EvidenceFailed(
            "artifact replay first divergence: case manifest identity"
        )
    if retained_terminal.get("artifact_root") != str(
        safe_relative(artifact_root, label="artifact root")
    ):
        raise EvidenceFailed(
            "artifact replay first divergence: terminal.json$.artifact_root"
        )
    if retained_terminal.get("artifact_dir") != str(
        safe_relative(input_dir, label="replay input")
    ):
        raise EvidenceFailed(
            "artifact replay first divergence: terminal.json$.artifact_dir"
        )

    rebuilt_inventory = independent_replay_projection(retained_source)
    captured = retained_source["captured"]
    rebuilt_partitions = partition_projection(rebuilt_inventory)
    rebuilt_plan = build_plan(rebuilt_inventory, captured["issues"])
    mode = str(retained_terminal.get("mode") or "")
    expected_retained_payloads = build_run_payloads(
        mode=mode,
        artifact_root=artifact_root,
        artifact_dir=input_dir,
        source=retained_source,
        inventory=rebuilt_inventory,
        partitions=rebuilt_partitions,
        plan=rebuilt_plan,
    )

    require_projection_equal(
        "inventory.json",
        rebuilt_inventory,
        retained_inventory,
    )
    require_projection_equal(
        "partitions.json",
        rebuilt_partitions,
        retained_partitions,
    )
    require_projection_equal("plan.json", rebuilt_plan, retained_plan)
    expected_event_rows = [
        json.loads(line)
        for line in expected_retained_payloads["events.jsonl"]
        .decode("utf-8")
        .splitlines()
    ]
    require_projection_equal(
        "events.jsonl",
        expected_event_rows,
        retained_events,
    )
    require_payload_equal(
        "events.jsonl",
        expected_retained_payloads["events.jsonl"],
        retained_event_bytes,
    )
    require_payload_equal(
        "reproduce.txt",
        expected_retained_payloads["reproduce.txt"],
        retained_reproduce,
    )
    expected_terminal = json.loads(expected_retained_payloads["terminal.json"])
    require_projection_equal(
        "terminal.json",
        expected_terminal,
        retained_terminal,
    )
    retained_payloads = {
        "source.json": canonical_bytes(retained_source),
        "inventory.json": canonical_bytes(retained_inventory),
        "partitions.json": canonical_bytes(retained_partitions),
        "plan.json": canonical_bytes(retained_plan),
        "events.jsonl": retained_event_bytes,
        "terminal.json": canonical_bytes(retained_terminal),
        "reproduce.txt": retained_reproduce,
    }
    for name in RUN_ARTIFACTS:
        require_payload_equal(
            name,
            expected_retained_payloads[name],
            retained_payloads[name],
        )

    retained_events_root = (
        "sha256-v1:" + hashlib.sha256(retained_event_bytes).hexdigest()
    )
    replay_payloads = build_replay_output_payloads(
        artifact_root=artifact_root,
        input_dir=input_dir,
        output_dir=output_dir,
        source=retained_source,
        inventory=rebuilt_inventory,
        partitions=rebuilt_partitions,
        plan=rebuilt_plan,
        retained_events_root=retained_events_root,
        retained_terminal_root=retained_terminal["semantic_root"],
        retained_artifact_identities=retained_terminal["artifact_identities"],
    )
    require_fresh_run_dir(target_dir, label="replay output directory")
    results = {
        name: write_once(target_dir / name, replay_payloads[name])
        for name in sorted(replay_payloads)
    }
    replay_terminal = json.loads(replay_payloads["terminal.json"])
    return {
        "terminal": "Pass",
        "mode": "replay",
        "subject_terminal": "REPLAY_EQUIVALENT",
        "input_dir": str(
            safe_relative(artifact_root, label="artifact root")
            / safe_relative(input_dir, label="replay input")
        ),
        "artifact_dir": str(
            safe_relative(artifact_root, label="artifact root")
            / safe_relative(output_dir, label="artifact dir")
        ),
        "source_root": retained_source["semantic_root"],
        "inventory_root": rebuilt_inventory["semantic_root"],
        "partitions_root": rebuilt_partitions["semantic_root"],
        "plan_root": rebuilt_plan["semantic_root"],
        "retained_events_root": retained_events_root,
        "retained_terminal_root": retained_terminal["semantic_root"],
        "replay_terminal_root": replay_terminal["semantic_root"],
        "first_divergence": None,
        "live_br_access": "forbidden-and-not-used",
        "writes": results,
        "no_claim": (
            "artifact-only replay proves deterministic reconstruction of the "
            "retained projection; it does not authenticate the producer, prove "
            "current tracker state, approve semantics, or mint capability authority"
        ),
    }


def _apply_issue_projection(issue: Mapping[str, Any]) -> dict[str, Any]:
    """Return the exact mutable-state projection, excluding br's update timestamp."""

    def edges(name: str) -> list[dict[str, Any]]:
        return sorted(
            (
                {
                    "id": str(edge.get("id", "")),
                    "type": str(
                        edge.get("dependency_type") or edge.get("type") or ""
                    ),
                    "status": str(edge.get("status", "")),
                    "priority": edge.get("priority"),
                }
                for edge in (issue.get(name) or [])
                if isinstance(edge, dict)
            ),
            key=lambda edge: (edge["type"], edge["id"]),
        )

    return {
        "id": str(issue.get("id", "")),
        "title": str(issue.get("title", "")),
        "type": str(issue.get("issue_type") or issue.get("type") or ""),
        "status": str(issue.get("status", "")),
        "priority": issue.get("priority"),
        "assignee": str(issue.get("assignee") or ""),
        "owner": str(issue.get("owner") or ""),
        "parent": str(issue.get("parent") or ""),
        "labels": sorted(str(label) for label in (issue.get("labels") or [])),
        "description": str(issue.get("description") or ""),
        "acceptance_criteria": str(issue.get("acceptance_criteria") or ""),
        "design": str(issue.get("design") or ""),
        "notes": str(issue.get("notes") or ""),
        "estimated_minutes": issue.get("estimated_minutes"),
        "dependencies": edges("dependencies"),
        "dependents": edges("dependents"),
    }


def _recovery_br_command(
    argv: Sequence[str],
    *,
    input_text: str | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run a bounded br recovery command even after cancellation was requested."""

    environment = os.environ.copy()
    for name in SANITIZED_ENV_NAMES:
        environment.pop(name, None)
    try:
        completed = subprocess.run(
            list(argv),
            cwd=REPO_ROOT,
            input=input_text,
            text=True,
            capture_output=True,
            env=environment,
            timeout=CAPS["subprocess_timeout_seconds"],
            check=False,
        )
    except FileNotFoundError as error:
        raise InfrastructureFailed("required recovery tool not found: br") from error
    except subprocess.TimeoutExpired as error:
        raise InfrastructureFailed("br recovery command exceeded its bounded timeout") from error
    if len(completed.stdout.encode("utf-8")) > CAPS["subprocess_stdout_bytes"]:
        raise InfrastructureFailed("br recovery stdout exceeded the bounded cap")
    return completed


def _read_apply_issue(issue_id: str, *, recovery: bool = False) -> dict[str, Any]:
    if recovery:
        completed = _recovery_br_command(
            ("br", "show", issue_id, "--json", *BR_READ_FLAGS)
        )
        if completed.returncode != 0:
            detail = completed.stderr.strip().splitlines()
            named = detail[0][:240] if detail else "no diagnostic"
            raise InfrastructureFailed(
                f"br recovery show returned {completed.returncode}: {named}"
            )
        try:
            document = json.loads(completed.stdout)
        except json.JSONDecodeError as error:
            raise InfrastructureFailed(
                f"br recovery show emitted malformed JSON for {issue_id}"
            ) from error
    else:
        document = br_read_json("show", issue_id, "--json")
    if not isinstance(document, list) or len(document) != 1:
        raise InfrastructureFailed(
            f"br show did not return exactly the explicit target {issue_id}"
        )
    if not isinstance(document[0], dict):
        raise InfrastructureFailed(f"br show returned a malformed row for {issue_id}")
    projection = _apply_issue_projection(document[0])
    if projection["id"] != issue_id:
        raise InfrastructureFailed(
            f"br show returned {projection['id']!r} for explicit target {issue_id}"
        )
    return projection


def _apply_update_argv(issue_id: str, agent_name: str) -> tuple[str, ...]:
    return (
        "br",
        "update",
        issue_id,
        "--description-file",
        "-",
        "--agent-name",
        agent_name,
        "--json",
    )


def _validate_update_stdout(
    completed: subprocess.CompletedProcess[str],
    issue_id: str,
) -> HarnessError | None:
    if completed.returncode != 0:
        detail = completed.stderr.strip().splitlines()
        named = detail[0][:240] if detail else "no diagnostic"
        return InfrastructureFailed(
            f"br update returned {completed.returncode} for {issue_id}: {named}"
        )
    try:
        parsed = json.loads(completed.stdout)
    except json.JSONDecodeError:
        return InfrastructureFailed(
            f"br update committed-state output was malformed for {issue_id}"
        )
    if not isinstance(parsed, (dict, list)) or not parsed:
        return InfrastructureFailed(f"br update returned no row for {issue_id}")
    return None


def _restore_apply_target(
    *,
    issue_id: str,
    agent_name: str,
    inverse_value: str,
    before_projection: Mapping[str, Any],
    after_projection: Mapping[str, Any],
) -> dict[str, Any]:
    before_root = semantic_root(before_projection)
    after_root = semantic_root(after_projection)
    current = _read_apply_issue(issue_id, recovery=True)
    current_root = semantic_root(current)
    if current_root == before_root:
        return {
            "state": "already-restored",
            "issue": issue_id,
            "before_root": before_root,
            "observed_root": current_root,
            "inverse_invoked": False,
            "verified": True,
        }
    if current_root != after_root:
        raise EvidenceFailed(
            f"{issue_id} has concurrent drift; exact-root restoration is unsafe"
        )

    completed = _recovery_br_command(
        _apply_update_argv(issue_id, agent_name),
        input_text=inverse_value,
    )
    output_failure = _validate_update_stdout(completed, issue_id)
    restored = _read_apply_issue(issue_id, recovery=True)
    restored_root = semantic_root(restored)
    if restored_root != before_root:
        detail = f"; transport={output_failure}" if output_failure else ""
        raise EvidenceFailed(
            f"{issue_id} inverse did not restore its exact pre-run root{detail}"
        )
    return {
        "state": "restored",
        "issue": issue_id,
        "before_root": before_root,
        "after_root": after_root,
        "observed_root": restored_root,
        "inverse_invoked": True,
        "inverse_stdout_valid": output_failure is None,
        "verified": True,
    }


def _apply_output_directory(
    artifact_root: str,
    artifact_dir: str,
    artifact_names: Sequence[str],
) -> tuple[Path, PurePosixPath]:
    root_rel = safe_relative(artifact_root, label="artifact root")
    output_rel = safe_relative(artifact_dir, label="artifact dir")
    combined = root_rel / output_rel
    output = resolve_safe(str(combined), label="apply artifact directory")
    if output.exists():
        raise InputRefused(
            f"apply artifact run directory already exists: {combined}"
        )
    for name in artifact_names:
        path = output / name
        if path.exists():
            raise InputRefused(
                f"refusing pre-existing apply artifact: {combined / name}"
            )
    return output, combined


def validate_apply_manifest(
    document: Mapping[str, Any],
    snapshot: LiveSnapshot,
) -> list[dict[str, Any]]:
    if document.get("schema") != APPLY_SCHEMA:
        raise InputRefused("apply manifest has an unknown schema")
    if document.get("inventory_root") != snapshot.inventory["semantic_root"]:
        raise InputRefused("apply manifest targets a stale inventory root")
    if document.get("plan_root") != snapshot.plan["semantic_root"]:
        raise InputRefused("apply manifest targets a stale plan root")
    if document.get("reviewed") is not True:
        raise InputRefused("apply manifest is not explicitly reviewed")
    reviewed_by = str(document.get("reviewed_by") or "").strip()
    agent_name = str(document.get("agent_name") or "").strip()
    if not reviewed_by:
        raise InputRefused("apply manifest lacks reviewed_by")
    if not agent_name:
        raise InputRefused("apply manifest lacks an explicit agent_name")
    if not str(document.get("reservation_receipt") or "").strip():
        raise InputRefused("apply manifest lacks a reservation receipt")

    rows = document.get("rows")
    if not isinstance(rows, list) or len(rows) != 1:
        raise InputRefused("V1 apply manifest must contain exactly one target row")
    if CAPS["apply_rows"] != 1:
        raise EvidenceFailed("V1 apply atomicity cap is not exactly one")
    if not isinstance(rows[0], dict):
        raise InputRefused("apply manifest target row is malformed")
    row = dict(rows[0])
    issue_id = str(row.get("id") or "")
    if not issue_id or document.get("target_id") != issue_id:
        raise InputRefused("apply manifest lacks one matching explicit target_id")
    if row.get("reviewed") is not True:
        raise InputRefused(f"{issue_id} row is not explicitly reviewed")
    if str(row.get("reviewed_by") or "").strip() != reviewed_by:
        raise InputRefused(f"{issue_id} row reviewed_by binding differs")
    if str(row.get("agent_name") or "").strip() != agent_name:
        raise InputRefused(f"{issue_id} row agent_name binding differs")

    inventory_rows = {candidate["id"]: candidate for candidate in snapshot.inventory["rows"]}
    live_issues = {candidate["id"]: candidate for candidate in snapshot.issues}
    proposal_rows = {
        candidate["id"]: candidate
        for candidate in snapshot.plan.get("section_name_only_proposals", [])
        if isinstance(candidate, dict) and candidate.get("id")
    }
    if issue_id not in inventory_rows or issue_id not in live_issues:
        raise InputRefused(f"apply target is outside the frozen inventory: {issue_id}")
    if issue_id not in proposal_rows:
        raise InputRefused(f"apply target lacks a current plan proposal: {issue_id}")

    inventory_row = inventory_rows[issue_id]
    live_issue = live_issues[issue_id]
    current_reviewers = {
        str(value).strip()
        for value in (live_issue.get("assignee"), live_issue.get("owner"))
        if str(value or "").strip()
    }
    if not current_reviewers or reviewed_by not in current_reviewers:
        raise InputRefused(
            f"{issue_id} reviewed_by is not the current assignee or owner"
        )
    if (
        str(document.get("current_owner") or "").strip() != reviewed_by
        or str(row.get("current_owner") or "").strip() != reviewed_by
    ):
        raise InputRefused(f"{issue_id} current-owner review binding differs")
    if inventory_row["status"] in {"deferred", "closed"}:
        raise InputRefused(
            f"{issue_id} status {inventory_row['status']} forbids V1 apply"
        )
    if inventory_row["disposition"] != "SECTION_NAME_ONLY":
        raise InputRefused(f"{issue_id} live disposition is not SECTION_NAME_ONLY")
    authority_domain = str(inventory_row.get("authority_domain") or "")
    if (
        document.get("authority_domain") != authority_domain
        or row.get("authority_domain") != authority_domain
    ):
        raise InputRefused(f"{issue_id} authority-domain binding differs")

    proposal = proposal_rows[issue_id]
    proposal_root = semantic_root(proposal)
    if (
        document.get("proposal_root") != proposal_root
        or row.get("proposal_root") != proposal_root
    ):
        raise InputRefused(f"{issue_id} proposal-root binding differs")
    proposal_fields = (
        "id",
        "disposition",
        "field",
        "missing_section",
        "old_root",
        "new_root",
        "old_value",
        "new_value",
        "inverse_value",
        "rationale",
    )
    for field in proposal_fields:
        if row.get(field) != proposal.get(field):
            raise InputRefused(f"{issue_id} differs from current proposal field {field}")

    old = str(live_issues[issue_id]["description"])
    new = str(row["new_value"])
    inverse = str(row["inverse_value"])
    if row["field"] != "description":
        raise InputRefused(f"{issue_id} may update only description")
    if row["old_root"] != text_root(old) or row["old_value"] != old:
        raise InputRefused(f"{issue_id} old field binding drifted")
    if row["new_root"] != text_root(new):
        raise InputRefused(f"{issue_id} new field root is inconsistent")
    if inverse != old:
        raise InputRefused(f"{issue_id} inverse does not exactly restore old text")
    if len(new.encode("utf-8")) > CAPS["field_bytes"]:
        raise InputRefused(f"{issue_id} proposed field exceeds cap")
    if old and not all(fragment in new for fragment in old.splitlines() if fragment):
        raise InputRefused(f"{issue_id} proposed field loses existing scope")
    section = str(row["missing_section"])
    if not section_present(new, section):
        raise InputRefused(f"{issue_id} proposed field lacks {section}")
    return [row]


def apply_reviewed_manifest(
    *,
    manifest_rel: str,
    artifact_root: str,
    artifact_dir: str,
    snapshot: LiveSnapshot,
) -> dict[str, Any]:
    manifest_path = resolve_safe(manifest_rel, label="apply manifest", must_exist=True)
    try:
        document = json.loads(bounded_read(manifest_path).decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise InputRefused(f"apply manifest is malformed: {error}") from error
    if not isinstance(document, dict):
        raise InputRefused("apply manifest root must be an object")
    rows = validate_apply_manifest(document, snapshot)
    row = rows[0]
    issue_id = str(row["id"])
    agent_name = str(document["agent_name"]).strip()
    reviewed_by = str(document["reviewed_by"]).strip()
    proposal_root = str(document["proposal_root"])
    manifest_root = semantic_root(document)
    artifact_names = (
        *RUN_ARTIFACTS,
        "apply-manifest.json",
        "inverse-br.json",
        "restoration.json",
    )
    output, output_rel = _apply_output_directory(
        artifact_root, artifact_dir, artifact_names
    )
    frozen_issue = _apply_issue_projection(
        {candidate["id"]: candidate for candidate in snapshot.issues}[issue_id]
    )

    # This read is the apply authority. The retained snapshot alone is never
    # allowed to authorize a mutation.
    before = _read_apply_issue(issue_id)
    before_root = semantic_root(before)
    if before_root != semantic_root(frozen_issue):
        raise InputRefused(f"{issue_id} drifted before apply pre-observation")
    after = dict(before)
    after["description"] = str(row["new_value"])
    after_root = semantic_root(after)
    forward_argv = _apply_update_argv(issue_id, agent_name)
    inverse_argv = _apply_update_argv(issue_id, agent_name)
    inverse_artifact = {
        "schema": "frankensim.beads-template-hygiene.inverse-br.v1",
        "manifest_root": manifest_root,
        "inventory_root": snapshot.inventory["semantic_root"],
        "plan_root": snapshot.plan["semantic_root"],
        "proposal_root": proposal_root,
        "target_id": issue_id,
        "reviewed_by": reviewed_by,
        "agent_name": agent_name,
        "before_issue_root": before_root,
        "planned_after_issue_root": after_root,
        "commands": [
            {
                "argv": list(inverse_argv),
                "stdin": str(row["inverse_value"]),
                "stdin_root": text_root(str(row["inverse_value"])),
                "expected_before_root": after_root,
                "expected_after_root": before_root,
            }
        ],
        "execution_order": "reverse",
        "retained_before_forward_mutation": True,
    }
    independent_findings = {
        str(candidate["id"]): dict(candidate.get("independent_findings", {}))
        for candidate in snapshot.inventory["rows"]
        if candidate.get("independent_findings")
    }
    retained_source = source_artifact(
        lint=snapshot.lint,
        issues=snapshot.issues,
        source=snapshot.source,
        independent_findings_by_id=independent_findings,
    )
    partitions = partition_projection(snapshot.inventory)
    reproduce_argv = [
        str(SCRIPT_REL),
        "--apply-manifest",
        manifest_rel,
        "--artifact-root",
        artifact_root,
        "--artifact-dir",
        artifact_dir,
    ]

    events: list[dict[str, Any]] = []

    def append_event(
        stage: str,
        *,
        result: str,
        command: Sequence[str] = (),
        inverse_command: Sequence[str] | None = None,
        terminal: str | None = None,
        old_root: str | None = None,
        new_root: str | None = None,
        extra: Mapping[str, Any] | None = None,
    ) -> None:
        event: dict[str, Any] = {
            "schema": EVENT_SCHEMA,
            "tool": "beads-template-hygiene",
            "rule": "reviewed-single-target-br-apply-v1",
            "source": retained_source["semantic_root"],
            "run": str(output_rel),
            "mode": "apply",
            "case": "template-lint.br-only-apply",
            "attempt": 1,
            "stage": stage,
            "sequence": len(events),
            "issue": issue_id,
            "warning": row["missing_section"],
            "disposition": row["disposition"],
            "shard": priority_group(
                {
                    candidate["id"]: candidate
                    for candidate in snapshot.inventory["rows"]
                }[issue_id]["priority"]
            ),
            "old_semantic_root": old_root,
            "new_semantic_root": new_root,
            "command": list(command),
            "result": result,
            "first_divergence": None,
            "caps": {
                "apply_rows": CAPS["apply_rows"],
                "events": EVENT_COUNT_CAP,
                "event_line_bytes": EVENT_LINE_CAP,
                "artifact_bytes": RUN_ARTIFACT_CAP,
            },
            "terminal": terminal,
            "inverse_br_command": (
                list(inverse_command) if inverse_command is not None else None
            ),
            "safe_relative_artifacts": list(artifact_names),
            "reproduction": reproduce_argv,
        }
        if extra:
            event.update(extra)
        missing = [field for field in EVENT_REQUIRED_FIELDS if field not in event]
        if missing:
            raise EvidenceFailed(
                f"apply event {stage} lacks required fields {sorted(missing)}"
            )
        if len(events) >= EVENT_COUNT_CAP:
            raise EvidenceFailed(f"apply event stream exceeds {EVENT_COUNT_CAP} events")
        if len(canonical_bytes(event)) > EVENT_LINE_CAP:
            raise EvidenceFailed(
                f"apply event sequence {event['sequence']} exceeds "
                f"{EVENT_LINE_CAP} bytes"
            )
        events.append(event)

    pre_payloads = {
        "source.json": canonical_bytes(retained_source),
        "inventory.json": canonical_bytes(snapshot.inventory),
        "partitions.json": canonical_bytes(partitions),
        "plan.json": canonical_bytes(snapshot.plan),
        "apply-manifest.json": canonical_bytes(document),
        "inverse-br.json": canonical_bytes(inverse_artifact),
        "reproduce.txt": canonical_bytes(reproduce_argv),
    }
    for name, payload in pre_payloads.items():
        if len(payload) > RUN_ARTIFACT_CAP:
            raise EvidenceFailed(
                f"apply pre-mutation artifact {name} exceeds its cap"
            )
    writes = {
        name: write_once(output / name, payload)
        for name, payload in sorted(pre_payloads.items())
    }

    # Retaining evidence is work between reads, so repeat the exact read guard
    # immediately before the only forward mutation.
    guarded_before = _read_apply_issue(issue_id)
    if semantic_root(guarded_before) != before_root:
        raise InputRefused(f"{issue_id} drifted after inverse retention")

    append_event(
        "start",
        result="admitted",
        extra={
            "inventory_root": snapshot.inventory["semantic_root"],
            "plan_root": snapshot.plan["semantic_root"],
            "proposal_root": proposal_root,
            "manifest_root": manifest_root,
        },
    )
    append_event(
        "source",
        result="captured",
        extra={"source_root": retained_source["semantic_root"]},
    )
    append_event(
        "inventory",
        result="bound",
        extra={"inventory_root": snapshot.inventory["semantic_root"]},
    )
    append_event("join", result="one-explicit-target")
    append_event("classify", result="section-name-only")
    append_event(
        "plan",
        result="exact-current-root-bound",
        extra={
            "plan_root": snapshot.plan["semantic_root"],
            "proposal_root": proposal_root,
        },
    )
    append_event(
        "reserve",
        result="receipt-bound-and-inverse-retained",
        inverse_command=inverse_argv,
        extra={
            "reservation_receipt": str(document["reservation_receipt"]),
            "inverse_artifact": str(output_rel / "inverse-br.json"),
        },
    )
    append_event(
        "pre-observe",
        result="exact-root-match",
        command=["br", "show", issue_id, "--json", *BR_READ_FLAGS],
        inverse_command=inverse_argv,
        old_root=before_root,
        new_root=after_root,
    )

    forward_attempted = False
    committed = False
    failure: HarnessError | None = None
    restoration: dict[str, Any] | None = None
    try:
        completed: subprocess.CompletedProcess[str] | None = None
        transport_failure: HarnessError | None = None
        try:
            forward_attempted = True
            completed = run_command(forward_argv, input_text=str(row["new_value"]))
        except HarnessError as error:
            transport_failure = error
        if completed is not None:
            output_failure = _validate_update_stdout(completed, issue_id)
            if output_failure is not None:
                transport_failure = output_failure

        # The command may have committed even when its return code or JSON
        # transport failed. Always decide from the exact immediate read.
        observed_after = _read_apply_issue(
            issue_id, recovery=_cancel_requested or transport_failure is not None
        )
        observed_after_root = semantic_root(observed_after)
        if observed_after_root == after_root:
            committed = True
        elif observed_after_root != before_root:
            raise EvidenceFailed(
                f"{issue_id} post-read found an unplanned concurrent state"
            )
        append_event(
            "apply",
            result="committed" if committed else "not-committed",
            command=forward_argv,
            inverse_command=inverse_argv,
            old_root=before_root,
            new_root=after_root,
            extra={"observed_semantic_root": observed_after_root},
        )
        if transport_failure is not None:
            raise transport_failure
        if not committed:
            raise EvidenceFailed(
                f"{issue_id} br update returned without the planned exact state"
            )

        append_event(
            "post-observe",
            result="exact-root-match",
            command=["br", "show", issue_id, "--json", *BR_READ_FLAGS],
            inverse_command=inverse_argv,
            old_root=before_root,
            new_root=after_root,
        )
        lint_result = br_read_json("lint", issue_id, "--json")
        if not isinstance(lint_result, dict) or not isinstance(
            lint_result.get("results"), list
        ):
            raise InfrastructureFailed(f"br lint returned malformed data for {issue_id}")
        remaining = {
            str(candidate.get("id", "")): candidate.get("missing") or []
            for candidate in lint_result["results"]
            if isinstance(candidate, dict)
        }
        if row["missing_section"] in remaining.get(issue_id, []):
            raise EvidenceFailed(f"{issue_id} still lacks {row['missing_section']}")
        if not section_present(observed_after["description"], row["missing_section"]):
            raise EvidenceFailed(
                f"{issue_id} independent heading check failed after update"
            )
        append_event(
            "export",
            result="planned-warning-cleared",
            command=["br", "lint", issue_id, "--json", *BR_READ_FLAGS],
            inverse_command=inverse_argv,
            old_root=before_root,
            new_root=after_root,
        )

        restoration = {
            "schema": "frankensim.beads-template-hygiene.restoration.v1",
            "state": "forward-accepted",
            "target_id": issue_id,
            "before_issue_root": before_root,
            "current_issue_root": observed_after_root,
            "planned_after_issue_root": after_root,
            "inverse_retained": True,
            "restoration_required": False,
            "verified": observed_after_root == after_root,
        }
        append_event(
            "publish",
            result="terminal-last",
            old_root=before_root,
            new_root=after_root,
            extra={"artifacts": list(artifact_names)},
        )
        append_event(
            "terminal",
            result="pass",
            terminal="Pass",
            old_root=before_root,
            new_root=after_root,
            extra={
                "exit_code": TERMINAL_EXIT["Pass"],
                "inventory_root": snapshot.inventory["semantic_root"],
                "partitions_root": partitions["semantic_root"],
                "plan_root": snapshot.plan["semantic_root"],
                "proposal_root": proposal_root,
                "manifest_root": manifest_root,
                "subject_terminal": "APPLIED_AND_INVERTIBLE",
                "event_count": len(events) + 1,
                "no_claim": snapshot.inventory["no_claim"],
            },
        )
        post_payloads = {
            "restoration.json": canonical_bytes(restoration),
            "events.jsonl": b"".join(canonical_bytes(event) for event in events),
        }
        for name, payload in post_payloads.items():
            if len(payload) > RUN_ARTIFACT_CAP:
                raise EvidenceFailed(
                    f"apply post-mutation artifact {name} exceeds its cap"
                )
        for name, payload in sorted(post_payloads.items()):
            writes[name] = write_once(output / name, payload)
        identity_schemas = {
            "source.json": SOURCE_SCHEMA,
            "inventory.json": INVENTORY_SCHEMA,
            "partitions.json": PARTITIONS_SCHEMA,
            "plan.json": PLAN_SCHEMA,
            "apply-manifest.json": APPLY_SCHEMA,
            "inverse-br.json": inverse_artifact["schema"],
            "restoration.json": restoration["schema"],
            "events.jsonl": EVENT_SCHEMA,
            "reproduce.txt": "frankensim.argv-json.v1",
        }
        identity_payloads = {**pre_payloads, **post_payloads}
        terminal: dict[str, Any] = {
            "schema": RUN_TERMINAL_SCHEMA,
            "mode": "apply",
            "terminal": "Pass",
            "subject_terminal": "APPLIED_AND_INVERTIBLE",
            "exit_code": TERMINAL_EXIT["Pass"],
            "artifact_root": str(
                safe_relative(artifact_root, label="artifact root")
            ),
            "artifact_dir": str(
                safe_relative(artifact_dir, label="artifact dir")
            ),
            "safe_relative_artifacts": list(artifact_names),
            "artifact_identities": {
                name: artifact_identity(
                    name, payload, identity_schemas[name]
                )
                for name, payload in sorted(identity_payloads.items())
            },
            "source_root": retained_source["semantic_root"],
            "inventory_root": snapshot.inventory["semantic_root"],
            "partitions_root": partitions["semantic_root"],
            "plan_root": snapshot.plan["semantic_root"],
            "proposal_root": proposal_root,
            "manifest_root": manifest_root,
            "before_issue_root": before_root,
            "after_issue_root": after_root,
            "event_count": len(events),
            "event_sequence": list(range(len(events))),
            "reproduction": reproduce_argv,
            "no_claim": snapshot.inventory["no_claim"],
        }
        terminal["semantic_root"] = semantic_root(terminal)
        # A Pass seal is deliberately the final filesystem operation. There is
        # no fallible publication step after it.
        writes["terminal.json"] = write_once(
            output / "terminal.json", canonical_bytes(terminal)
        )
    except HarnessError as error:
        failure = error
    except BaseException as error:
        failure = InfrastructureFailed(f"unexpected apply failure: {error}")

    if failure is not None:
        if forward_attempted:
            try:
                restoration = _restore_apply_target(
                    issue_id=issue_id,
                    agent_name=agent_name,
                    inverse_value=str(row["inverse_value"]),
                    before_projection=before,
                    after_projection=after,
                )
                idempotent = _restore_apply_target(
                    issue_id=issue_id,
                    agent_name=agent_name,
                    inverse_value=str(row["inverse_value"]),
                    before_projection=before,
                    after_projection=after,
                )
                if idempotent["state"] != "already-restored":
                    raise EvidenceFailed(
                        f"{issue_id} restoration was not idempotent"
                    )
                restoration["idempotent_recheck"] = idempotent
                append_event(
                    "inverse",
                    result=str(restoration["state"]),
                    command=(
                        inverse_argv
                        if restoration.get("inverse_invoked")
                        else ()
                    ),
                    inverse_command=inverse_argv,
                    old_root=after_root,
                    new_root=before_root,
                )
                append_event(
                    "restore",
                    result="idempotent-exact-root-verified",
                    inverse_command=inverse_argv,
                    old_root=after_root,
                    new_root=before_root,
                )
            except HarnessError as rollback_error:
                raise EvidenceFailed(
                    f"{issue_id} apply failed and restoration is incomplete: "
                    f"{rollback_error}"
                ) from failure

        if restoration is None:
            restoration = {
                "schema": "frankensim.beads-template-hygiene.restoration.v1",
                "state": "forward-not-committed",
                "target_id": issue_id,
                "before_issue_root": before_root,
                "current_issue_root": semantic_root(
                    _read_apply_issue(issue_id, recovery=True)
                ),
                "planned_after_issue_root": after_root,
                "inverse_retained": True,
                "restoration_required": False,
                "verified": True,
            }
        append_event(
            "refuse",
            result=failure.terminal,
            inverse_command=inverse_argv,
            old_root=before_root,
            new_root=after_root,
            extra={"detail": str(failure)},
        )
        append_event(
            "terminal",
            result="refused",
            terminal=failure.terminal,
            old_root=before_root,
            new_root=after_root,
            extra={
                "exit_code": TERMINAL_EXIT[failure.terminal],
                "inventory_root": snapshot.inventory["semantic_root"],
                "partitions_root": partitions["semantic_root"],
                "plan_root": snapshot.plan["semantic_root"],
                "proposal_root": proposal_root,
                "manifest_root": manifest_root,
                "restored": bool(restoration and restoration.get("verified")),
                "event_count": len(events) + 1,
                "no_claim": snapshot.inventory["no_claim"],
            },
        )
        failure_payloads: dict[str, bytes] = {
            "events.jsonl": b"".join(canonical_bytes(event) for event in events),
            "restoration.json": canonical_bytes(restoration),
        }
        publication_errors: list[str] = []
        for name, payload in sorted(failure_payloads.items()):
            try:
                writes[name] = write_once(output / name, payload)
            except HarnessError as publication_error:
                publication_errors.append(f"{name}: {publication_error}")
        failure_terminal: dict[str, Any] = {
            "schema": RUN_TERMINAL_SCHEMA,
            "mode": "apply",
            "terminal": failure.terminal,
            "exit_code": TERMINAL_EXIT[failure.terminal],
            "artifact_root": str(
                safe_relative(artifact_root, label="artifact root")
            ),
            "artifact_dir": str(
                safe_relative(artifact_dir, label="artifact dir")
            ),
            "safe_relative_artifacts": list(artifact_names),
            "source_root": retained_source["semantic_root"],
            "inventory_root": snapshot.inventory["semantic_root"],
            "partitions_root": partitions["semantic_root"],
            "plan_root": snapshot.plan["semantic_root"],
            "proposal_root": proposal_root,
            "manifest_root": manifest_root,
            "before_issue_root": before_root,
            "planned_after_issue_root": after_root,
            "restored": bool(restoration.get("verified")),
            "event_count": len(events),
            "event_sequence": list(range(len(events))),
            "reproduction": reproduce_argv,
            "detail": str(failure),
            "no_claim": snapshot.inventory["no_claim"],
        }
        failure_terminal["semantic_root"] = semantic_root(failure_terminal)
        try:
            writes["terminal.json"] = write_once(
                output / "terminal.json", canonical_bytes(failure_terminal)
            )
        except HarnessError as publication_error:
            publication_errors.append(f"terminal.json: {publication_error}")
        if publication_errors:
            raise EvidenceFailed(
                "apply failed after exact restoration, but failure evidence "
                f"publication was incomplete: {'; '.join(publication_errors)}"
            ) from failure
        raise failure

    return {
        "terminal": "Pass",
        "mode": "apply",
        "rows": 1,
        "issue_ids": [issue_id],
        "artifact_dir": str(output_rel),
        "inventory_root": snapshot.inventory["semantic_root"],
        "plan_root": snapshot.plan["semantic_root"],
        "proposal_root": proposal_root,
        "manifest_root": manifest_root,
        "writes": writes,
        "no_claim": snapshot.inventory["no_claim"],
    }


def fixture_issue(
    *,
    issue_id: str = "fixture-1",
    title: str = "Specific fixture behavior",
    issue_type: str = "task",
    status: str = "open",
    priority: int = 1,
    description: str = "Specific work with deterministic boundaries.",
    acceptance: str = "",
    design: str = "",
    notes: str = "",
    assignee: str = "",
    owner: str = "",
    parent: str = "",
    labels: Sequence[str] = ("reality-check",),
    dependencies: Sequence[Mapping[str, Any]] = (),
    children: Sequence[str] = (),
) -> dict[str, Any]:
    return {
        "id": issue_id,
        "title": title,
        "type": issue_type,
        "status": status,
        "priority": priority,
        "assignee": assignee,
        "owner": owner,
        "parent": parent,
        "labels": list(labels),
        "description": description,
        "acceptance_criteria": acceptance,
        "design": design,
        "notes": notes,
        "estimated_minutes": 60,
        "updated_at": "frozen",
        "dependencies": [dict(edge) for edge in dependencies],
        "dependents": [
            {"id": child, "type": "parent-child", "status": "open", "priority": 1}
            for child in children
        ],
    }


FIXTURE_REQUIRED_SECTIONS_BY_TYPE = {
    "bug": ("## Acceptance Criteria", "## Steps to Reproduce"),
    "task": ("## Acceptance Criteria",),
    "feature": ("## Acceptance Criteria",),
    "epic": ("## Success Criteria",),
    "chore": (),
    "docs": (),
    "question": (),
    "custom": (),
}
FIXTURE_MAX_WARNINGS_PER_ISSUE = 16
FIXTURE_BR_ROOT_REL = PurePosixPath(
    "target/beads-template-hygiene/self-test-br-fixture"
)
FIXTURE_BR_DB_REL = FIXTURE_BR_ROOT_REL / ".beads" / "beads.db"
FIXTURE_BR_DB_ARG = ".beads/beads.db"
FIXTURE_BR_STABLE_TITLE = "Template hygiene no-mock stable fixture"
FIXTURE_BR_REGRESSION_TITLE = "Template hygiene no-mock regression fixture"
_fixture_br_context_cache: dict[str, str] | None = None


def fixture_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise EvidenceFailed(
            f"{label}: expected {expected!r}, observed {actual!r}"
        )


def fixture_true(condition: bool, label: str) -> None:
    if not condition:
        raise EvidenceFailed(label)


def fixture_expected_missing(issue_type: str) -> tuple[str, ...]:
    return tuple(sorted(FIXTURE_REQUIRED_SECTIONS_BY_TYPE.get(issue_type, ())))


def fixture_lint(
    issues: Sequence[Mapping[str, Any]],
    missing: Mapping[str, Sequence[str]],
) -> dict[str, Any]:
    if len(issues) > CAPS["issues"]:
        raise EvidenceFailed(
            f"fixture issue count exceeds {CAPS['issues']} cap"
        )
    issue_ids = [str(issue.get("id") or "") for issue in issues]
    if any(not issue_id for issue_id in issue_ids):
        raise EvidenceFailed("fixture contains an empty issue ID")
    if len(issue_ids) != len(set(issue_ids)):
        raise EvidenceFailed("fixture contains duplicate issue IDs")
    extra_missing_ids = sorted(set(missing) - set(issue_ids))
    if extra_missing_ids:
        raise EvidenceFailed(
            f"fixture missing map has unknown issue {extra_missing_ids[0]}"
        )

    normalized_missing: dict[str, tuple[str, ...]] = {}
    warning_total = 0
    for issue in issues:
        issue_id = str(issue["id"])
        status = str(issue.get("status") or "")
        priority = issue.get("priority")
        issue_type = str(issue.get("type") or "")
        if status not in STATUS_SCOPES:
            raise EvidenceFailed(
                f"fixture {issue_id} has unknown status {status!r}"
            )
        if not isinstance(priority, int) or priority not in range(5):
            raise EvidenceFailed(
                f"fixture {issue_id} has priority outside P0-P4"
            )
        if not issue_type:
            raise EvidenceFailed(f"fixture {issue_id} has empty type")
        raw_sections = tuple(str(section) for section in missing.get(issue_id, ()))
        if any(not section for section in raw_sections):
            raise EvidenceFailed(
                f"fixture {issue_id} contains an empty missing section"
            )
        if len(raw_sections) != len(set(raw_sections)):
            raise EvidenceFailed(
                f"fixture {issue_id} contains duplicate missing sections"
            )
        sections = tuple(sorted(raw_sections))
        if len(sections) > FIXTURE_MAX_WARNINGS_PER_ISSUE:
            raise EvidenceFailed(
                f"fixture {issue_id} exceeds per-issue warning cap"
            )
        normalized_missing[issue_id] = sections
        warning_total += len(sections)
    if warning_total > CAPS["warnings"]:
        raise EvidenceFailed(
            f"fixture warning count exceeds {CAPS['warnings']} cap"
        )

    result: dict[str, Any] = {}
    for scope in LINT_SCOPES:
        selected: list[dict[str, Any]] = []
        for issue in sorted(issues, key=lambda row: str(row["id"])):
            sections = normalized_missing[str(issue["id"])]
            if not sections:
                continue
            if scope != "all" and issue["status"] != scope:
                continue
            selected.append(
                {
                "id": issue["id"],
                "title": issue["title"],
                "type": issue["type"],
                    "missing": list(sections),
                    "warnings": len(sections),
                    "suggestions": [
                        {
                            "section": section,
                            "hint": "fixture suggestion is non-authoritative",
                        }
                        for section in sections
                    ],
                }
            )
        result[scope] = {
            "total": sum(row["warnings"] for row in selected),
            "issues": len(selected),
            "results": selected,
        }
    return result


def fixture_source(
    issues: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    ordered_issues = sorted(
        (dict(issue) for issue in issues),
        key=lambda issue: str(issue.get("id") or ""),
    )
    source: dict[str, Any] = {
        "schema": SOURCE_SCHEMA,
        "br_version": "fixture",
        "case_manifest_root": "fixture",
        "live_issue_count": len(ordered_issues),
        "live_issue_projection_root": semantic_root(ordered_issues),
        "lint_issue_projection_root": "fixture",
        "status_summary": {
            status: sum(issue.get("status") == status for issue in ordered_issues)
            for status in STATUS_SCOPES
        },
        "files": [],
        "status_cut": list(LINT_SCOPES),
    }
    source["semantic_root"] = semantic_root(source)
    return source


def fixture_inventory(
    issues: Sequence[Mapping[str, Any]],
    missing: Mapping[str, Sequence[str]],
    *,
    findings: Mapping[str, Mapping[str, str]] | None = None,
) -> dict[str, Any]:
    lint = fixture_lint(issues, missing)
    issue_ids = {
        str(row["id"])
        for row in lint["all"]["results"]
    } | set((findings or {}).keys())
    selected = [issue for issue in issues if str(issue["id"]) in issue_ids]
    return assemble_inventory(
        lint,
        selected,
        fixture_source(issues),
        findings,
    )


def fixture_reviewed_apply_admission() -> tuple[LiveSnapshot, dict[str, Any]]:
    acceptance = (
        "1. Run `tests/e2e/template_apply.rs` against the explicit fixture.\n"
        "2. Retain deterministic logs, inverse replay, and source closure.\n"
        "3. Refuse stale roots, wrong owners, and any substantive text change.\n"
        "4. Preserve the no-claim boundary after the bounded test passes."
    )
    issue = fixture_issue(
        issue_id="fixture-reviewed-apply",
        description="Exact fixture scope remains byte-preserved.",
        acceptance=acceptance,
        assignee="FixtureOwner",
        labels=("authority:template-hygiene",),
    )
    lint = fixture_lint(
        [issue],
        {issue["id"]: ("## Acceptance Criteria",)},
    )
    source = fixture_source([issue])
    inventory = assemble_inventory(lint, [issue], source)
    inventory["rows"][0]["disposition"] = "SECTION_NAME_ONLY"
    inventory["rows"][0]["rationale"] = (
        "fixture owner review is explicitly bound to the exact roots"
    )
    inventory["rows"][0]["semantic_flags"]["section_only_review_bound"] = True
    inventory["warning_rows"][0]["disposition"] = "SECTION_NAME_ONLY"
    inventory["counts"]["by_disposition"] = {"SECTION_NAME_ONLY": 1}
    inventory.pop("semantic_root", None)
    inventory["semantic_root"] = semantic_root(inventory)
    plan = build_plan(inventory, [issue])
    proposal = dict(plan["section_name_only_proposals"][0])
    proposal_root = semantic_root(proposal)
    row = dict(proposal)
    row.update(
        {
            "reviewed": True,
            "reviewed_by": "FixtureOwner",
            "current_owner": "FixtureOwner",
            "agent_name": "FixtureAgent",
            "authority_domain": "authority:template-hygiene",
            "proposal_root": proposal_root,
        }
    )
    document = {
        "schema": APPLY_SCHEMA,
        "target_id": issue["id"],
        "inventory_root": inventory["semantic_root"],
        "plan_root": plan["semantic_root"],
        "proposal_root": proposal_root,
        "reviewed": True,
        "reviewed_by": "FixtureOwner",
        "current_owner": "FixtureOwner",
        "agent_name": "FixtureAgent",
        "reservation_receipt": "fixture-root-bound-receipt",
        "authority_domain": "authority:template-hygiene",
        "rows": [row],
    }
    return LiveSnapshot(lint, [issue], source, inventory, plan), document


def expect_error(
    error_type: type[HarnessError],
    callback: Callable[[], Any],
    *,
    contains: str | None = None,
) -> HarnessError:
    try:
        callback()
    except error_type as error:
        if contains is not None and contains not in str(error):
            raise EvidenceFailed(
                f"{error_type.__name__} diagnostic lacks {contains!r}: {error}"
            ) from error
        return error
    raise EvidenceFailed(f"expected {error_type.__name__} was not raised")


def fixture_br_argv(*arguments: str) -> tuple[str, ...]:
    return (
        "br",
        *arguments,
        "--db",
        FIXTURE_BR_DB_ARG,
        "--no-auto-flush",
        "--no-auto-import",
        "--no-color",
    )


def fixture_br_run(
    *arguments: str,
    input_text: str | None = None,
) -> subprocess.CompletedProcess[str]:
    check_cancel()
    fixture_root = resolve_safe(
        str(FIXTURE_BR_ROOT_REL),
        label="no-mock fixture root",
        must_exist=True,
    )
    environment = os.environ.copy()
    for name in SANITIZED_ENV_NAMES:
        environment.pop(name, None)
    try:
        completed = subprocess.run(
            list(fixture_br_argv(*arguments)),
            cwd=fixture_root,
            input=input_text,
            text=True,
            capture_output=True,
            env=environment,
            timeout=CAPS["subprocess_timeout_seconds"],
            check=False,
        )
    except FileNotFoundError as error:
        raise InfrastructureFailed("required fixture tool not found: br") from error
    except subprocess.TimeoutExpired as error:
        raise InfrastructureFailed(
            "fixture br command exceeded the bounded timeout"
        ) from error
    if len(completed.stdout.encode("utf-8")) > CAPS["subprocess_stdout_bytes"]:
        raise InfrastructureFailed("fixture br stdout exceeded the bounded cap")
    if completed.returncode != 0:
        diagnostic = completed.stderr.strip() or completed.stdout.strip()
        first = diagnostic.splitlines()[0][:240] if diagnostic else "no diagnostic"
        raise InfrastructureFailed(
            f"fixture br {arguments[0]} returned {completed.returncode}: {first}"
        )
    check_cancel()
    return completed


def fixture_br_json(
    *arguments: str,
    input_text: str | None = None,
) -> Any:
    completed = fixture_br_run(
        *arguments,
        "--json",
        input_text=input_text,
    )
    try:
        return json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise InfrastructureFailed(
            f"fixture br {arguments[0]} emitted malformed JSON"
        ) from error


def fixture_br_list() -> list[dict[str, Any]]:
    document = fixture_br_json("list", "--all", "--deferred", "--limit", "0")
    issues = document.get("issues") if isinstance(document, dict) else None
    if not isinstance(issues, list) or document.get("has_more"):
        raise InfrastructureFailed("fixture br list is malformed or truncated")
    return [dict(issue) for issue in issues]


def fixture_br_show(issue_id: str) -> dict[str, Any]:
    document = fixture_br_json("show", issue_id)
    if not isinstance(document, list) or len(document) != 1:
        raise InfrastructureFailed(
            f"fixture br show did not return exactly {issue_id}"
        )
    issue = dict(document[0])
    if str(issue.get("id")) != issue_id:
        raise InfrastructureFailed("fixture br show returned the wrong issue")
    return issue


def fixture_br_semantic_projection(issue: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "id": str(issue.get("id") or ""),
        "title": str(issue.get("title") or ""),
        "description": str(issue.get("description") or ""),
        "acceptance_criteria": str(issue.get("acceptance_criteria") or ""),
        "design": str(issue.get("design") or ""),
        "notes": str(issue.get("notes") or ""),
        "status": str(issue.get("status") or ""),
        "priority": issue.get("priority"),
        "issue_type": str(issue.get("issue_type") or issue.get("type") or ""),
        "assignee": str(issue.get("assignee") or ""),
        "owner": str(issue.get("owner") or ""),
        "parent": str(issue.get("parent") or ""),
        "labels": sorted(str(label) for label in (issue.get("labels") or [])),
        "dependencies": sorted(
            (
                str(edge.get("dependency_type") or edge.get("type") or ""),
                str(edge.get("id") or ""),
            )
            for edge in (issue.get("dependencies") or [])
            if isinstance(edge, dict)
        ),
    }


def fixture_br_find_by_title(title: str) -> dict[str, Any] | None:
    matches = [issue for issue in fixture_br_list() if issue.get("title") == title]
    if len(matches) > 1:
        raise EvidenceFailed(f"fixture contains duplicate title {title!r}")
    return matches[0] if matches else None


def fixture_br_create_seed(
    *,
    title: str,
    slug: str,
    status: str,
    description: str,
) -> dict[str, Any]:
    fixture_br_json(
        "create",
        "--title",
        title,
        "--slug",
        slug,
        "--type",
        "task",
        "--priority",
        "4",
        "--status",
        status,
        "--description-file",
        "-",
        "--labels",
        "template-hygiene-fixture",
        "--ephemeral",
        input_text=description,
    )
    created = fixture_br_find_by_title(title)
    if created is None:
        raise InfrastructureFailed(f"fixture br create did not retain {title!r}")
    return created


def fixture_br_set_status(issue_id: str, desired: str) -> None:
    current = fixture_br_show(issue_id)
    status = str(current.get("status") or "")
    if status == desired:
        return
    if desired == "open" and status == "closed":
        fixture_br_json(
            "reopen",
            issue_id,
            "--reason",
            "template hygiene fixture restoration",
        )
        return
    if desired == "closed":
        fixture_br_json(
            "close",
            issue_id,
            "--reason",
            "template hygiene fixture baseline",
        )
        return
    fixture_br_json("update", issue_id, "--status", desired)


def fixture_br_set_description(issue_id: str, description: str) -> None:
    fixture_br_json(
        "update",
        issue_id,
        "--description-file",
        "-",
        input_text=description,
    )


def fixture_br_context() -> dict[str, str]:
    global _fixture_br_context_cache
    if _fixture_br_context_cache is not None:
        return dict(_fixture_br_context_cache)
    try:
        fixture_root = resolve_safe(
            str(FIXTURE_BR_ROOT_REL),
            label="no-mock fixture root",
        )
        fixture_root.mkdir(parents=True, exist_ok=True)
        fixture_db = REPO_ROOT.joinpath(*FIXTURE_BR_DB_REL.parts)
        if not fixture_db.exists():
            fixture_br_run("init", "--prefix", "thfx")
            if not fixture_db.exists():
                raise InfrastructureFailed(
                    "fixture br init did not create its isolated database"
                )

        stable = fixture_br_find_by_title(FIXTURE_BR_STABLE_TITLE)
        if stable is None:
            stable = fixture_br_create_seed(
                title=FIXTURE_BR_STABLE_TITLE,
                slug="template-hygiene-stable",
                status="open",
                description=(
                    "Stable no-mock fixture baseline.\n\n"
                    "Completion is checked only inside this isolated br database."
                ),
            )
        stable_id = str(stable["id"])
        fixture_br_set_status(stable_id, "open")

        regression = fixture_br_find_by_title(FIXTURE_BR_REGRESSION_TITLE)
        if regression is None:
            regression = fixture_br_create_seed(
                title=FIXTURE_BR_REGRESSION_TITLE,
                slug="template-hygiene-regression",
                status="closed",
                description="Intentionally missing an acceptance section.",
            )
        regression_id = str(regression["id"])
        fixture_br_set_status(regression_id, "closed")

        _fixture_br_context_cache = {
            "db": str(FIXTURE_BR_DB_REL),
            "stable_id": stable_id,
            "regression_id": regression_id,
        }
        return dict(_fixture_br_context_cache)
    except InfrastructureFailed as error:
        raise NoData(
            "no-mock br fixture unavailable without touching live Beads; "
            f"persistent safe target is {FIXTURE_BR_ROOT_REL}: {error}"
        ) from error


def fixture_br_description_round_trip(
    case_id: str,
    *,
    assert_stale_guard: bool = False,
    inject_failure: bool = False,
) -> dict[str, Any]:
    context = fixture_br_context()
    issue_id = context["stable_id"]
    fixture_br_set_status(issue_id, "open")
    before_issue = fixture_br_show(issue_id)
    before = fixture_br_semantic_projection(before_issue)
    old_description = before["description"]
    planned_root = text_root(old_description)
    new_description = (
        old_description.rstrip()
        + f"\n\n## Acceptance Criteria\n\nNo-mock fixture case `{case_id}`.\n"
    )
    observed_failure = False
    try:
        fixture_br_set_description(issue_id, new_description)
        changed = fixture_br_show(issue_id)
        fixture_true(
            text_root(str(changed.get("description") or "")) != planned_root,
            f"{case_id}: br update did not change the field root",
        )
        if assert_stale_guard:
            expect_error(
                InputRefused,
                lambda: (
                    None
                    if text_root(str(changed.get("description") or ""))
                    == planned_root
                    else (_ for _ in ()).throw(
                        InputRefused("stale exact-root apply refused")
                    )
                ),
                contains="stale exact-root",
            )
            observed_failure = True
        if inject_failure:
            try:
                raise InfrastructureFailed("seeded post-update fixture fault")
            except InfrastructureFailed:
                observed_failure = True
    finally:
        current = fixture_br_show(issue_id)
        if str(current.get("description") or "") != old_description:
            fixture_br_set_description(issue_id, old_description)
    restored = fixture_br_semantic_projection(fixture_br_show(issue_id))
    fixture_equal(restored, before, f"{case_id}: semantic restoration")
    if assert_stale_guard or inject_failure:
        fixture_true(observed_failure, f"{case_id}: fault/refusal was not observed")
    return {
        "checks": 5 if assert_stale_guard or inject_failure else 4,
        "fixture_db": context["db"],
        "semantic_restoration": "exact-excluding-append-only-audit-history",
    }


def fixture_br_regression_round_trip() -> dict[str, Any]:
    context = fixture_br_context()
    issue_id = context["regression_id"]
    fixture_br_set_status(issue_id, "closed")
    before = fixture_br_semantic_projection(fixture_br_show(issue_id))
    fixture_equal(before["status"], "closed", "regression fixture baseline status")
    try:
        fixture_br_set_status(issue_id, "open")
        lint = fixture_br_json("lint", issue_id)
        result_ids = {
            str(row.get("id") or "")
            for row in (lint.get("results") or [])
            if isinstance(row, dict)
        }
        fixture_true(
            issue_id in result_ids,
            "new-issue regression did not enter active lint membership",
        )
    finally:
        fixture_br_set_status(issue_id, "closed")
    restored = fixture_br_semantic_projection(fixture_br_show(issue_id))
    fixture_equal(
        restored,
        before,
        "new-issue regression semantic restoration excluding audit history",
    )
    return {
        "checks": 4,
        "fixture_db": context["db"],
        "semantic_restoration": "active-debt-and-issue-fields; audit history append-only",
    }


def self_test_case(
    case_id: str,
    *,
    live_snapshot: LiveSnapshot | None,
) -> dict[str, Any]:
    strong_acceptance = (
        "## Detailed completion contract\n\n"
        "1. Run `scripts/ci/beads_template_hygiene.sh --self-test` and require "
        "deterministic unit, boundary, mutation, cancellation, and replay checks.\n"
        "2. Refuse malformed input, stale source identity, copied criteria, and "
        "unreviewed semantic changes with an exact typed terminal.\n"
        "3. Retain bounded redacted JSONL logs, a disjoint artifact replay, and "
        "the first divergent issue or field root.\n"
        "4. Close only after source closure and independent verification; this "
        "fixture mints no product, scientific, release, or owner authority."
    )

    if case_id == "template-lint.inventory-empty":
        inventory = fixture_inventory([], {})
        fixture_equal(inventory["counts"]["issues"], 0, "empty issue count")
        fixture_equal(inventory["counts"]["warnings"], 0, "empty warning count")
        fixture_equal(
            set(inventory["status_cuts"]),
            set(LINT_SCOPES),
            "empty status-cut coverage",
        )

        one = fixture_issue(issue_id="fixture-one", priority=4)
        one_inventory = fixture_inventory(
            [one],
            {"fixture-one": ("## Acceptance Criteria",)},
        )
        fixture_equal(one_inventory["counts"]["issues"], 1, "one-row boundary")
        fixture_equal(one_inventory["counts"]["warnings"], 1, "one-warning boundary")

        maximum = [
            fixture_issue(issue_id=f"fixture-max-{index:04d}")
            for index in range(CAPS["issues"])
        ]
        maximum_missing = {
            issue["id"]: ("## Acceptance Criteria",) for issue in maximum
        }
        maximum_lint = fixture_lint(maximum, maximum_missing)
        fixture_equal(
            maximum_lint["all"]["issues"],
            CAPS["issues"],
            "maximum issue cap",
        )
        expect_error(
            EvidenceFailed,
            lambda: fixture_lint(
                maximum
                + [fixture_issue(issue_id="fixture-over-issue-cap")],
                {
                    **maximum_missing,
                    "fixture-over-issue-cap": ("## Acceptance Criteria",),
                },
            ),
            contains="issue count",
        )

        sixteen_sections = tuple(
            f"## Synthetic fixture section {index:02d}"
            for index in range(FIXTURE_MAX_WARNINGS_PER_ISSUE)
        )
        warning_max_issues = [
            fixture_issue(issue_id=f"fixture-warning-{index:04d}")
            for index in range(CAPS["warnings"] // len(sixteen_sections))
        ]
        warning_max_missing = {
            issue["id"]: sixteen_sections for issue in warning_max_issues
        }
        warning_max_lint = fixture_lint(
            warning_max_issues,
            warning_max_missing,
        )
        fixture_equal(
            warning_max_lint["all"]["total"],
            CAPS["warnings"],
            "maximum warning cap",
        )
        over_warning_issue = fixture_issue(issue_id="fixture-warning-over")
        expect_error(
            EvidenceFailed,
            lambda: fixture_lint(
                warning_max_issues + [over_warning_issue],
                {
                    **warning_max_missing,
                    over_warning_issue["id"]: ("## Acceptance Criteria",),
                },
            ),
            contains="warning count",
        )
        expect_error(
            EvidenceFailed,
            lambda: fixture_lint(
                [fixture_issue(issue_id="fixture-bad-status", status="unknown")],
                {"fixture-bad-status": ("## Acceptance Criteria",)},
            ),
            contains="unknown status",
        )
        expect_error(
            EvidenceFailed,
            lambda: fixture_lint(
                [fixture_issue(issue_id="fixture-bad-priority", priority=5)],
                {"fixture-bad-priority": ("## Acceptance Criteria",)},
            ),
            contains="outside P0-P4",
        )
        expect_error(
            EvidenceFailed,
            lambda: fixture_lint(
                [
                    fixture_issue(issue_id="fixture-duplicate"),
                    fixture_issue(issue_id="fixture-duplicate"),
                ],
                {"fixture-duplicate": ("## Acceptance Criteria",)},
            ),
            contains="duplicate issue IDs",
        )
        return {
            "category": "schema-boundary-and-caps",
            "checks": 12,
            "max_issues": CAPS["issues"],
            "max_warnings": CAPS["warnings"],
        }

    if case_id == "template-lint.inventory-live":
        fixture_true(live_snapshot is not None, "live inventory did not execute")
        assert live_snapshot is not None
        fixture_equal(
            live_snapshot.inventory["schema"],
            INVENTORY_SCHEMA,
            "live inventory schema",
        )
        fixture_equal(
            set(live_snapshot.inventory["status_cuts"]),
            set(LINT_SCOPES),
            "live status-cut coverage",
        )
        fixture_equal(
            live_snapshot.inventory["counts"]["warnings"],
            len(live_snapshot.inventory["warning_rows"]),
            "live warning arithmetic",
        )

        required_types = ("bug", "task", "feature", "epic")
        universe: list[dict[str, Any]] = []
        missing: dict[str, tuple[str, ...]] = {}
        for status_index, status in enumerate(STATUS_SCOPES):
            for type_index, issue_type in enumerate(required_types):
                issue_id = f"fixture-{status}-{issue_type}"
                issue = fixture_issue(
                    issue_id=issue_id,
                    title=(
                        "Δeterministic Unicode fixture 🧪"
                        if not universe
                        else f"{status} {issue_type} fixture"
                    ),
                    issue_type=issue_type,
                    status=status,
                    priority=(status_index * len(required_types) + type_index) % 5,
                    description="Unicode units: μm, °C, and λ remain canonical.",
                )
                universe.append(issue)
                missing[issue_id] = fixture_expected_missing(issue_type)
        for issue_type in ("chore", "docs", "question", "custom"):
            issue_id = f"fixture-open-{issue_type}"
            universe.append(
                fixture_issue(
                    issue_id=issue_id,
                    issue_type=issue_type,
                    status="open",
                    priority=4,
                )
            )
            missing[issue_id] = fixture_expected_missing(issue_type)

        for issue_type, expected in FIXTURE_REQUIRED_SECTIONS_BY_TYPE.items():
            fixture_equal(
                tuple(sorted(REQUIRED_SECTIONS_BY_TYPE.get(issue_type, ()))),
                tuple(sorted(expected)),
                f"rule matrix for {issue_type}",
            )

        lint = fixture_lint(universe, missing)
        inventory = fixture_inventory(universe, missing)
        warned_ids = {row["id"] for row in inventory["rows"]}
        for issue_type in ("chore", "docs", "question", "custom"):
            fixture_true(
                f"fixture-open-{issue_type}" not in warned_ids,
                f"{issue_type} incorrectly acquired a default lint requirement",
            )
        fixture_equal(
            set(inventory["counts"]["by_status"]),
            set(STATUS_SCOPES),
            "all statuses represented",
        )
        fixture_equal(
            set(inventory["counts"]["by_priority"]),
            {f"P{priority}" for priority in range(5)},
            "all priorities represented",
        )
        fixture_equal(
            set(inventory["counts"]["by_type"]),
            set(required_types),
            "all required-section types represented",
        )
        fixture_true(
            "Δeterministic Unicode fixture 🧪"
            in {row["title"] for row in inventory["rows"]},
            "Unicode title did not survive inventory projection",
        )
        reversed_inventory = fixture_inventory(
            list(reversed(universe)),
            dict(reversed(list(missing.items()))),
        )
        fixture_equal(
            reversed_inventory["semantic_root"],
            inventory["semantic_root"],
            "input-order invariant inventory",
        )
        fixture_equal(
            lint["all"]["total"],
            sum(lint[scope]["total"] for scope in STATUS_SCOPES),
            "all-status warning union",
        )
        return {
            "category": "live-and-universe-projection",
            "checks": 18,
            "statuses": list(STATUS_SCOPES),
            "priorities": [0, 1, 2, 3, 4],
            "types": list(FIXTURE_REQUIRED_SECTIONS_BY_TYPE),
            "live_warnings": live_snapshot.inventory["counts"]["warnings"],
        }

    if case_id == "template-lint.overlap-partitions":
        sections = (
            "## Acceptance Criteria",
            "## Steps to Reproduce",
            "## Success Criteria",
        )
        combinations = (
            (sections[0],),
            (sections[1],),
            (sections[2],),
            (sections[0], sections[1]),
            (sections[0], sections[2]),
            (sections[1], sections[2]),
            sections,
        )
        issues = [
            fixture_issue(
                issue_id=f"fixture-overlap-{index}",
                issue_type="bug",
                status=STATUS_SCOPES[index % len(STATUS_SCOPES)],
            )
            for index in range(len(combinations))
        ]
        missing = {
            issue["id"]: combination
            for issue, combination in zip(issues, combinations, strict=True)
        }
        inventory = fixture_inventory(issues, missing)
        fixture_equal(inventory["counts"]["issues"], 7, "seven overlap issues")
        fixture_equal(inventory["counts"]["warnings"], 12, "overlap warning sum")
        expected_partitions = {
            "+".join(sorted(combination)): 1 for combination in combinations
        }
        fixture_equal(
            inventory["counts"]["overlap_partitions"],
            dict(sorted(expected_partitions.items())),
            "seven exclusive overlap partitions",
        )
        reversed_inventory = fixture_inventory(
            list(reversed(issues)),
            dict(reversed(list(missing.items()))),
        )
        fixture_equal(
            reversed_inventory["semantic_root"],
            inventory["semantic_root"],
            "overlap input-order invariance",
        )

        duplicate_lint = json.loads(json.dumps(fixture_lint(issues, missing)))
        duplicate = dict(duplicate_lint["all"]["results"][0])
        duplicate_lint["all"]["results"].append(duplicate)
        duplicate_lint["all"]["issues"] += 1
        duplicate_lint["all"]["total"] += duplicate["warnings"]
        expect_error(
            EvidenceFailed,
            lambda: assemble_inventory(
                duplicate_lint,
                issues,
                fixture_source(issues),
            ),
            contains="duplicate issue",
        )

        union_gap_lint = json.loads(json.dumps(fixture_lint(issues, missing)))
        gap_issue = issues[0]
        gap_scope = str(gap_issue["status"])
        gap_rows = union_gap_lint[gap_scope]["results"]
        removed = gap_rows.pop(0)
        union_gap_lint[gap_scope]["issues"] -= 1
        union_gap_lint[gap_scope]["total"] -= removed["warnings"]
        expect_error(
            EvidenceFailed,
            lambda: assemble_inventory(
                union_gap_lint,
                issues,
                fixture_source(issues),
            ),
            contains="exact union",
        )
        return {
            "category": "exact-partition-metamorphic",
            "checks": 7,
            "exclusive_partitions": 7,
            "warnings": 12,
        }

    if case_id == "template-lint.section-only":
        issue = fixture_issue(acceptance=strong_acceptance)
        reviewed = classify_issue(
            issue,
            ("## Acceptance Criteria",),
            section_only_reviewed=True,
        )[0]
        fixture_equal(reviewed, "SECTION_NAME_ONLY", "reviewed section-only")
        unreviewed = classify_issue(
            issue,
            ("## Acceptance Criteria",),
        )[0]
        fixture_equal(
            unreviewed,
            "OWNER_REVIEW_REQUIRED",
            "unreviewed polished criteria",
        )
        duplicate = classify_issue(
            issue,
            ("## Acceptance Criteria",),
            duplicate_acceptance=True,
            section_only_reviewed=True,
        )[0]
        fixture_equal(
            duplicate,
            "OWNER_REVIEW_REQUIRED",
            "duplicate generic criteria",
        )
        reproduction_issue = fixture_issue(
            issue_id="fixture-reviewed-repro",
            issue_type="bug",
            description=(
                "REPRO: run `scripts/ci/beads_template_hygiene.sh --self-test`. "
                "Expected Pass; observed failure before the heading repair."
            ),
        )
        reproduction = classify_issue(
            reproduction_issue,
            ("## Steps to Reproduce",),
            section_only_reviewed=True,
        )[0]
        fixture_equal(
            reproduction,
            "SECTION_NAME_ONLY",
            "reviewed reproduction heading-only repair",
        )

        old_description = str(issue["description"])
        proposal = proposal_for_row(
            {
                "id": issue["id"],
                "disposition": "SECTION_NAME_ONLY",
                "missing_sections": ["## Acceptance Criteria"],
                "rationale": "reviewed fixture",
            },
            issue,
        )
        fixture_true(proposal is not None, "section-only proposal was not built")
        assert proposal is not None
        fixture_true(
            old_description in proposal["new_value"],
            "section-only proposal lost existing description",
        )
        fixture_true(
            "## Acceptance Criteria" in proposal["new_value"],
            "section-only proposal lacks literal heading",
        )
        fixture_equal(
            proposal["inverse_value"],
            old_description,
            "section-only exact inverse",
        )
        return {
            "category": "review-bound-section-only",
            "checks": 8,
            "review_required": True,
        }

    if case_id == "template-lint.substantive-omission":
        empty = classify_issue(
            fixture_issue(acceptance=""),
            ("## Acceptance Criteria",),
        )[0]
        fixture_equal(
            empty,
            "SUBSTANTIVE_SEMANTIC_OMISSION",
            "empty acceptance omission",
        )
        placeholder = classify_issue(
            fixture_issue(acceptance="TBD"),
            ("## Acceptance Criteria",),
            {"## Acceptance Criteria": "empty-or-placeholder-body"},
        )[0]
        fixture_equal(
            placeholder,
            "SUBSTANTIVE_SEMANTIC_OMISSION",
            "placeholder omission",
        )
        missing_repro = classify_issue(
            fixture_issue(
                issue_type="bug",
                description="The operation fails without a reachable reproducer.",
            ),
            ("## Steps to Reproduce",),
        )[0]
        fixture_equal(
            missing_repro,
            "SUBSTANTIVE_SEMANTIC_OMISSION",
            "missing bug reproducer",
        )
        mixed = classify_issue(
            fixture_issue(issue_type="bug", acceptance=strong_acceptance),
            ("## Acceptance Criteria", "## Steps to Reproduce"),
            section_only_reviewed=True,
        )[0]
        fixture_equal(
            mixed,
            "SUBSTANTIVE_SEMANTIC_OMISSION",
            "mixed warning cannot become heading-only",
        )
        return {
            "category": "substantive-fail-closed",
            "checks": 4,
            "mixed_warning_refused": True,
        }

    if case_id == "template-lint.wrong-type":
        non_bug_types = (
            "task",
            "feature",
            "epic",
            "chore",
            "docs",
            "question",
            "custom",
        )
        for issue_type in non_bug_types:
            disposition = classify_issue(
                fixture_issue(issue_type=issue_type),
                ("## Steps to Reproduce",),
            )[0]
            fixture_equal(
                disposition,
                "MALFORMED_OR_WRONG_TYPE",
                f"wrong-type {issue_type}",
            )
        task_root = semantic_root(fixture_issue(issue_type="task"))
        feature_root = semantic_root(fixture_issue(issue_type="feature"))
        fixture_true(
            task_root != feature_root,
            "issue-type mutation did not move semantic root",
        )
        malformed = classify_issue(
            fixture_issue(
                description=(
                    "Acceptance criteria: exact unit, e2e, log, replay, "
                    "authority, no-claim, and source closure obligations."
                ),
                acceptance=strong_acceptance,
            ),
            ("## Acceptance Criteria",),
        )[0]
        fixture_equal(
            malformed,
            "MALFORMED_OR_WRONG_TYPE",
            "malformed heading classification",
        )
        return {
            "category": "type-and-heading-schema",
            "checks": len(non_bug_types) + 2,
            "no_rule_types": ["chore", "docs", "question", "custom"],
        }

    if case_id == "template-lint.rollup-gap":
        issue = fixture_issue(
            issue_type="epic",
            acceptance="Every frozen child must pass its exact close gate.",
            children=("fixture-1.1", "fixture-1.2"),
        )
        disposition = classify_issue(issue, ("## Success Criteria",))[0]
        fixture_equal(
            disposition,
            "ROLLUP_CHILD_SET_GAP",
            "missing exact child set",
        )
        reversed_children = fixture_issue(
            issue_type="epic",
            acceptance=issue["acceptance_criteria"],
            children=("fixture-1.2", "fixture-1.1"),
        )
        fixture_equal(
            classify_issue(reversed_children, ("## Success Criteria",))[0],
            disposition,
            "child-order invariance",
        )
        complete = fixture_issue(
            issue_type="epic",
            acceptance=(
                "Children fixture-1.1 and fixture-1.2 must each pass before "
                "the compatible terminal close gate."
            ),
            children=("fixture-1.1", "fixture-1.2"),
        )
        fixture_equal(
            classify_issue(complete, ("## Success Criteria",))[0],
            "OWNER_REVIEW_REQUIRED",
            "complete rollup still requires owner review",
        )
        return {
            "category": "rollup-exact-child-set",
            "checks": 3,
            "children": 2,
        }

    if case_id == "template-lint.owner-review":
        polished = fixture_issue(
            acceptance=strong_acceptance.replace(
                "template_hygiene",
                "template_hygiene_unique",
            ),
            assignee="BlueLake",
        )
        fixture_equal(
            classify_issue(polished, ("## Acceptance Criteria",))[0],
            "OWNER_REVIEW_REQUIRED",
            "polished unique prose without review",
        )
        deferred = fixture_issue(
            status="deferred",
            acceptance=strong_acceptance,
        )
        fixture_equal(
            classify_issue(
                deferred,
                ("## Acceptance Criteria",),
                section_only_reviewed=True,
            )[0],
            "OWNER_REVIEW_REQUIRED",
            "deferred issue cannot be apply-eligible",
        )
        capped = classify_issue(
            fixture_issue(acceptance=strong_acceptance),
            ("## Acceptance Criteria",),
            {"## Acceptance Criteria": "semantic-field-cap-exceeded"},
            section_only_reviewed=True,
        )[0]
        fixture_equal(
            capped,
            "OWNER_REVIEW_REQUIRED",
            "over-cap semantic field",
        )
        active_owner = fixture_issue(
            status="in_progress",
            acceptance=(
                "Owner has discussed the specific implementation, test surface, "
                "logging boundary, replay inputs, and source-closure question, "
                "but has not yet bound an exact reviewed completion contract."
            ),
            assignee="BlueLake",
        )
        fixture_equal(
            classify_issue(active_owner, ("## Acceptance Criteria",))[0],
            "OWNER_REVIEW_REQUIRED",
            "active owner adjudication",
        )
        return {
            "category": "owner-authority-binding",
            "checks": 4,
            "deferred_apply_eligible": False,
        }

    if case_id == "template-lint.historical-review":
        closed_issues: list[dict[str, Any]] = []
        missing: dict[str, tuple[str, ...]] = {}
        for index, issue_type in enumerate(FIXTURE_REQUIRED_SECTIONS_BY_TYPE):
            issue = fixture_issue(
                issue_id=f"fixture-closed-{issue_type}",
                issue_type=issue_type,
                status="closed",
                priority=index % 5,
                acceptance="TBD",
            )
            closed_issues.append(issue)
            required = fixture_expected_missing(issue_type)
            missing[issue["id"]] = required or ("## Acceptance Criteria",)
            fixture_equal(
                classify_issue(issue, missing[issue["id"]])[0],
                "HISTORICAL_IMMUTABLE_REVIEW",
                f"closed {issue_type}",
            )
        inventory = fixture_inventory(closed_issues, missing)
        plan = build_plan(inventory, closed_issues)
        fixture_equal(plan["shards"], [], "closed history excluded from shards")
        fixture_equal(
            plan["section_name_only_proposals"],
            [],
            "closed history excluded from proposals",
        )
        return {
            "category": "historical-immutability",
            "checks": len(closed_issues) + 2,
            "closed_types": list(FIXTURE_REQUIRED_SECTIONS_BY_TYPE),
        }

    if case_id in {
        "template-lint.p0-shard",
        "template-lint.p1-shard",
        "template-lint.p2-p3-shard",
    }:
        priorities = {
            "template-lint.p0-shard": (0,),
            "template-lint.p1-shard": (1,),
            "template-lint.p2-p3-shard": (2, 3, 4),
        }[case_id]
        issues = [
            fixture_issue(
                issue_id=f"fixture-p{priority}-{index:02d}",
                priority=priority,
                acceptance=strong_acceptance,
                assignee="FixtureOwner",
                labels=("reality-check", "authority:template-hygiene"),
            )
            for priority in priorities
            for index in range(27)
        ]
        missing = {issue["id"]: ("## Acceptance Criteria",) for issue in issues}
        inventory = fixture_inventory(issues, missing)
        plan = build_plan(inventory, issues)
        shards = [
            row for row in plan["shards"] if row["priority"] in priorities
        ]
        fixture_equal(
            len(shards),
            2 * len(priorities),
            "bounded shard count",
        )
        for priority in priorities:
            priority_shards = [
                row for row in shards if row["priority"] == priority
            ]
            fixture_equal(
                [len(row["issue_ids"]) for row in priority_shards],
                [25, 2],
                f"P{priority} shard boundaries",
            )
        fixture_true(
            all(
                len(row["issue_ids"]) <= CAPS["plan_shard_rows"]
                and row["owner"] == "FixtureOwner"
                and row["authority_domain"] == "authority:template-hygiene"
                for row in shards
            ),
            "shard crossed size, owner, or authority boundary",
        )
        fixture_equal(
            [row["priority"] for row in shards],
            [priority for priority in priorities for _ in range(2)],
            "priority-stable shard order",
        )
        reversed_plan = build_plan(
            fixture_inventory(
                list(reversed(issues)),
                dict(reversed(list(missing.items()))),
            ),
            list(reversed(issues)),
        )
        fixture_equal(
            reversed_plan["semantic_root"],
            plan["semantic_root"],
            "plan input-order invariance",
        )
        fixture_equal(
            plan["apply_contract"]["max_rows"],
            1,
            "single-target V1 apply cap",
        )
        fixture_equal(
            plan["section_name_only_proposals"],
            [],
            "unreviewed rows cannot produce apply proposals",
        )
        return {
            "category": "deterministic-bounded-sharding",
            "checks": 7 + len(priorities),
            "priorities": list(priorities),
            "shards": len(shards),
            "plan_cap": CAPS["plan_shard_rows"],
            "apply_cap": CAPS["apply_rows"],
        }

    if case_id == "template-lint.br-only-apply":
        result = fixture_br_description_round_trip(case_id)
        argv = fixture_br_argv(
            "update",
            result.get("stable_id", "explicit-id"),
            "--description-file",
            "-",
            "--json",
        )
        fixture_equal(argv[0], "br", "br-only mutation transport")
        fixture_true(
            not any(token in argv for token in ("sqlite3", "sed", "perl")),
            "non-br tracker mutation transport present",
        )
        fixture_true("--db" in argv, "fixture mutation lacks explicit database")
        fixture_true("--force" not in argv, "fixture mutation uses --force")
        admission_snapshot, admission_document = fixture_reviewed_apply_admission()
        validated = validate_apply_manifest(
            admission_document,
            admission_snapshot,
        )
        fixture_equal(
            [row["id"] for row in validated],
            ["fixture-reviewed-apply"],
            "review-bound single-target apply admission",
        )
        wrong_owner = json.loads(json.dumps(admission_document))
        wrong_owner["reviewed_by"] = "WrongOwner"
        wrong_owner["current_owner"] = "WrongOwner"
        wrong_owner["rows"][0]["reviewed_by"] = "WrongOwner"
        wrong_owner["rows"][0]["current_owner"] = "WrongOwner"
        expect_error(
            InputRefused,
            lambda: validate_apply_manifest(wrong_owner, admission_snapshot),
            contains="current assignee or owner",
        )
        stale_proposal = json.loads(json.dumps(admission_document))
        stale_proposal["proposal_root"] = "stale-proposal-root"
        stale_proposal["rows"][0]["proposal_root"] = "stale-proposal-root"
        expect_error(
            InputRefused,
            lambda: validate_apply_manifest(
                stale_proposal,
                admission_snapshot,
            ),
            contains="proposal-root",
        )
        return {
            "category": "no-mock-br-only-round-trip",
            "checks": result["checks"] + 7,
            "fixture_db": result["fixture_db"],
            "restoration": result["semantic_restoration"],
        }

    if case_id == "template-lint.concurrent-drift":
        result = fixture_br_description_round_trip(
            case_id,
            assert_stale_guard=True,
        )
        return {
            "category": "no-mock-concurrent-drift",
            "checks": result["checks"],
            "fixture_db": result["fixture_db"],
            "restoration": result["semantic_restoration"],
        }

    if case_id == "template-lint.partial-batch":
        result = fixture_br_description_round_trip(
            case_id,
            inject_failure=True,
        )
        global _cancel_requested
        previous_cancel = _cancel_requested
        try:
            _cancel_requested = True
            expect_error(
                CancelledDrained,
                check_cancel,
                contains="no mutation is in flight",
            )
        finally:
            _cancel_requested = previous_cancel
        check_cancel()
        return {
            "category": "fault-cancel-and-compensate",
            "checks": result["checks"] + 2,
            "fixture_db": result["fixture_db"],
            "restoration": result["semantic_restoration"],
        }

    if case_id == "template-lint.copied-boilerplate":
        duplicated = [
            fixture_issue(
                issue_id=f"fixture-copy-{index}",
                acceptance=strong_acceptance,
            )
            for index in range(2)
        ]
        duplicate_inventory = fixture_inventory(
            duplicated,
            {
                issue["id"]: ("## Acceptance Criteria",)
                for issue in duplicated
            },
        )
        fixture_equal(
            {
                row["disposition"]
                for row in duplicate_inventory["rows"]
            },
            {"OWNER_REVIEW_REQUIRED"},
            "byte-identical generic criteria",
        )
        polished_unique = fixture_issue(
            issue_id="fixture-polished-unique",
            acceptance=strong_acceptance.replace(
                "source closure",
                "source closure for fixture-polished-unique",
            ),
        )
        fixture_equal(
            classify_issue(
                polished_unique,
                ("## Acceptance Criteria",),
            )[0],
            "OWNER_REVIEW_REQUIRED",
            "polished unique heuristic cannot authorize",
        )
        cosmetic = fixture_issue(
            issue_id="fixture-cosmetic-copy",
            acceptance=strong_acceptance.replace("deterministic", "bit-stable"),
        )
        fixture_equal(
            classify_issue(cosmetic, ("## Acceptance Criteria",))[0],
            "OWNER_REVIEW_REQUIRED",
            "cosmetic paraphrase cannot authorize",
        )
        return {
            "category": "anti-boilerplate-authority",
            "checks": 4,
            "duplicate_rows": 2,
        }

    if case_id == "template-lint.scope-loss":
        old_issues = [
            fixture_issue(issue_id="fixture-scope-a"),
            fixture_issue(issue_id="fixture-scope-b"),
        ]
        new_issues = [
            fixture_issue(issue_id="fixture-scope-a"),
            fixture_issue(issue_id="fixture-scope-c"),
        ]
        old_missing = {
            issue["id"]: ("## Acceptance Criteria",) for issue in old_issues
        }
        new_missing = {
            issue["id"]: ("## Acceptance Criteria",) for issue in new_issues
        }
        old_inventory = fixture_inventory(old_issues, old_missing)
        new_inventory = fixture_inventory(new_issues, new_missing)
        fixture_equal(
            old_inventory["counts"]["issues"],
            new_inventory["counts"]["issues"],
            "count-preserving scope fixture",
        )
        fixture_true(
            old_inventory["semantic_root"] != new_inventory["semantic_root"],
            "count-preserving membership substitution escaped identity",
        )

        stale_lint = fixture_lint(old_issues, old_missing)
        stale_lint["all"] = fixture_lint(new_issues, new_missing)["all"]
        expect_error(
            EvidenceFailed,
            lambda: assemble_inventory(
                stale_lint,
                old_issues,
                fixture_source(old_issues),
            ),
            contains="exact projection",
        )
        return {
            "category": "count-preserving-scope-identity",
            "checks": 3,
            "old_ids": ["fixture-scope-a", "fixture-scope-b"],
            "new_ids": ["fixture-scope-a", "fixture-scope-c"],
        }

    if case_id == "template-lint.inverse-replay":
        result = fixture_br_description_round_trip(case_id)
        return {
            "category": "no-mock-inverse-replay",
            "checks": result["checks"],
            "fixture_db": result["fixture_db"],
            "restoration": result["semantic_restoration"],
        }

    if case_id == "template-lint.new-issue-regression":
        result = fixture_br_regression_round_trip()
        return {
            "category": "no-mock-active-debt-regression",
            "checks": result["checks"],
            "fixture_db": result["fixture_db"],
            "restoration": result["semantic_restoration"],
        }

    if case_id == "template-lint.zero-debt-closeout":
        empty_lint = fixture_lint([], {})
        empty_source = fixture_source()
        empty_inventory = assemble_inventory(empty_lint, [], empty_source)
        replayed = independent_replay_projection(
            {"lint": empty_lint, "issues": [], "source": empty_source}
        )
        fixture_equal(
            replayed["semantic_root"],
            empty_inventory["semantic_root"],
            "zero-debt independent reconstruction",
        )

        one = fixture_issue(issue_id="fixture-new-warning")
        one_inventory = fixture_inventory(
            [one],
            {one["id"]: ("## Acceptance Criteria",)},
        )

        def require_zero(inventory: Mapping[str, Any]) -> None:
            if inventory["counts"]["warnings"] != 0:
                raise NoData(
                    f"{inventory['counts']['warnings']} warning rows remain"
                )

        require_zero(empty_inventory)
        expect_error(
            NoData,
            lambda: require_zero(one_inventory),
            contains="warning rows remain",
        )
        return {
            "category": "zero-debt-independent-closeout",
            "checks": 4,
            "negative_warning_count": 1,
        }

    if case_id == "template-lint.artifact-replay":
        issues = [
            fixture_issue(
                issue_id="fixture-replay-a",
                acceptance=strong_acceptance,
            ),
            fixture_issue(
                issue_id="fixture-replay-b",
                issue_type="bug",
                description=(
                    "REPRO: run `scripts/ci/beads_template_hygiene.sh "
                    "--negative template-lint.artifact-replay`; expected "
                    "refusal, observed failure."
                ),
            ),
        ]
        missing = {
            "fixture-replay-a": ("## Acceptance Criteria",),
            "fixture-replay-b": ("## Steps to Reproduce",),
        }
        lint = fixture_lint(issues, missing)
        source = fixture_source(issues)
        retained = assemble_inventory(lint, issues, source)
        replayed = independent_replay_projection(
            {"lint": lint, "issues": list(reversed(issues)), "source": source}
        )
        fixture_equal(
            retained["semantic_root"],
            replayed["semantic_root"],
            "artifact replay root",
        )
        mutated_source = dict(source)
        mutated_source["case_manifest_root"] = "mutated-fixture"
        mutated_source["semantic_root"] = semantic_root(
            {
                key: value
                for key, value in mutated_source.items()
                if key != "semantic_root"
            }
        )
        mutated = independent_replay_projection(
            {"lint": lint, "issues": issues, "source": mutated_source}
        )
        fixture_true(
            mutated["semantic_root"] != retained["semantic_root"],
            "source identity mutation did not move replay root",
        )
        expect_error(
            InputRefused,
            lambda: independent_replay_projection(
                {"lint": lint, "issues": issues}
            ),
            contains="lacks lint/issues/source",
        )
        bad_lint = json.loads(json.dumps(lint))
        bad_lint["all"]["total"] += 1
        expect_error(
            EvidenceFailed,
            lambda: independent_replay_projection(
                {"lint": bad_lint, "issues": issues, "source": source}
            ),
            contains="warning arithmetic",
        )
        return {
            "category": "artifact-only-independent-replay",
            "checks": 5,
            "live_tracker_access": "none",
        }

    raise EvidenceFailed(f"self-test case has no implementation: {case_id}")


def run_self_tests(
    case_manifest: Mapping[str, Any],
    *,
    selected: str | None = None,
) -> dict[str, Any]:
    cases = [
        row
        for row in case_manifest["case"]
        if selected is None or row["id"] == selected
    ]
    if selected is not None and not cases:
        raise UsageRefused(f"unknown case: {selected}")
    live_needed = any(
        row["id"] == "template-lint.inventory-live" for row in cases
    )
    live_snapshot = collect_live(case_manifest) if live_needed else None
    results: list[dict[str, Any]] = []
    for sequence, row in enumerate(cases):
        check_cancel()
        case_id = row["id"]
        try:
            evidence = self_test_case(case_id, live_snapshot=live_snapshot)
        except HarnessError as error:
            failure = {
                "schema": EVENT_SCHEMA,
                "mode": "self-test" if selected is None else "negative",
                "case": case_id,
                "attempt": 1,
                "stage": "case-terminal",
                "sequence": sequence,
                "terminal": error.terminal,
                "exit_code": TERMINAL_EXIT[error.terminal],
                "expected_terminal": row["expected_terminal"],
                "expected_subject_terminal": row["expected_subject_terminal"],
                "mutation": row["mutation"],
                "replay": row["replay"],
                "result": "FAIL",
                "detail": str(error),
                "authority": row["authority"],
                "no_claim": row["no_claim"],
            }
            results.append(failure)
            json_stdout(failure)
            raise EvidenceFailed(f"{case_id}: {error}") from error
        result = {
            "schema": EVENT_SCHEMA,
            "mode": "self-test" if selected is None else "negative",
            "case": case_id,
            "attempt": 1,
            "stage": "case-terminal",
            "sequence": sequence,
            "terminal": "Pass",
            "exit_code": 0,
            "expected_terminal": row["expected_terminal"],
            "expected_subject_terminal": row["expected_subject_terminal"],
            "expected_inner_br_category": row["expected_inner_br_category"],
            "mutation": row["mutation"],
            "replay": row["replay"],
            "result": "PASS",
            "authority": row["authority"],
            "no_claim": row["no_claim"],
            "evidence": evidence,
        }
        results.append(result)
        json_stdout(result)
    fixture_true(
        len(results) == len(cases),
        "self-test result count differs from selected case count",
    )
    fixture_equal(
        [row["case"] for row in results],
        [row["id"] for row in cases],
        "self-test result order",
    )
    summary = {
        "schema": SUMMARY_SCHEMA,
        "mode": "self-test" if selected is None else "negative",
        "stage": "terminal",
        "sequence": len(results),
        "terminal": "Pass",
        "exit_code": 0,
        "cases": len(results),
        "case_ids": [row["case"] for row in results],
        "total_checks": sum(
            int(row["evidence"].get("checks", 0)) for row in results
        ),
        "categories": [
            str(row["evidence"].get("category", "")) for row in results
        ],
        "fixture_mutation_cases": [
            row["case"] for row in results if row["mutation"] != "none"
        ],
        "artifact_replay_cases": [
            row["case"] for row in results if row["replay"] == "artifact-only"
        ],
        "case_manifest_root": case_manifest["semantic_root"],
        "no_claim": (
            "self-tests validate harness and isolated fixture behavior; they "
            "do not repair or close live target Beads, authenticate owner "
            "review, or confer implementation, scientific, or release authority"
        ),
    }
    json_stdout(summary)
    return summary


def v2_rooted(document: Mapping[str, Any]) -> dict[str, Any]:
    result = dict(document)
    result.pop("semantic_root", None)
    result["semantic_root"] = semantic_root(result)
    return result


def v2_exact_keys(
    value: Mapping[str, Any],
    expected: Iterable[str],
    *,
    label: str,
) -> None:
    actual = set(value)
    required = set(expected)
    if actual != required:
        missing = sorted(required - actual)
        extra = sorted(actual - required)
        raise InputRefused(
            f"{label} has a non-closed schema; missing={missing}, extra={extra}"
        )


def v2_assert_unique(values: Sequence[str], *, label: str) -> None:
    if len(values) != len(set(values)):
        raise InputRefused(f"{label} contains duplicate values")


def v2_parse_csv_set(
    value: str | None,
    *,
    allowed: Iterable[str],
    label: str,
) -> set[str] | None:
    if value is None:
        return None
    if not value or any(part != part.strip() for part in value.split(",")):
        raise UsageRefused(f"{label} must be a non-empty comma-separated set")
    selected = set(value.split(","))
    if "" in selected:
        raise UsageRefused(f"{label} contains an empty member")
    unknown = sorted(selected - set(allowed))
    if unknown:
        raise UsageRefused(f"{label} contains unknown values {unknown}")
    return selected


def v2_stable_issue_projection(issue: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "id": str(issue.get("id") or ""),
        "title": str(issue.get("title") or ""),
        "type": str(issue.get("type") or issue.get("issue_type") or ""),
        "status": str(issue.get("status") or ""),
        "priority": issue.get("priority"),
        "assignee": str(issue.get("assignee") or ""),
        "owner": str(issue.get("owner") or ""),
        "parent": str(issue.get("parent") or ""),
        "labels": sorted(str(value) for value in (issue.get("labels") or [])),
        "field_roots": dict(sorted((issue.get("field_roots") or {}).items())),
        "missing_sections": sorted(
            str(value) for value in (issue.get("missing_sections") or [])
        ),
        "disposition": str(issue.get("disposition") or ""),
        "dependencies": sorted(
            (
                str(edge.get("type") or edge.get("dependency_type") or ""),
                str(edge.get("id") or ""),
                str(edge.get("status") or ""),
                edge.get("priority"),
            )
            for edge in (issue.get("dependencies") or [])
            if isinstance(edge, dict)
        ),
        "dependents": sorted(
            (
                str(edge.get("type") or edge.get("dependency_type") or ""),
                str(edge.get("id") or ""),
                str(edge.get("status") or ""),
                edge.get("priority"),
            )
            for edge in (issue.get("dependents") or [])
            if isinstance(edge, dict)
        ),
    }


def v2_target_root(issue: Mapping[str, Any]) -> str:
    return semantic_root(v2_stable_issue_projection(issue))


def v2_dependency_neighborhood_root(issue: Mapping[str, Any]) -> str:
    return semantic_root(
        {
            "dependencies": v2_stable_issue_projection(issue)["dependencies"],
            "dependents": v2_stable_issue_projection(issue)["dependents"],
        }
    )


def v2_complete_clause_roots(issue: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for field in ("description", "acceptance_criteria", "design", "notes"):
        text = str(issue.get(field) or "")
        for line_index, line in enumerate(text.splitlines()):
            normalized = " ".join(line.split())
            if not normalized:
                continue
            lowered = normalized.casefold()
            if normalized.startswith("#") or any(
                term in lowered for term in CLAUSE_TERMS
            ):
                encoded = normalized.encode("utf-8")
                rows.append(
                    {
                        "field": field,
                        "line_index": line_index,
                        "byte_length": len(encoded),
                        "text_root": text_root(normalized),
                        "truncated": False,
                    }
                )
    return rows


def v2_domain_candidates(issue: Mapping[str, Any]) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    authority_domain = str(issue.get("authority_domain") or "")
    if authority_domain:
        candidates.append(
            {
                "candidate": authority_domain,
                "source": "v1-classification-hint",
                "provenance_root": text_root(authority_domain),
                "falsifier": (
                    "a root-bound reviewer may reject this inferred domain; "
                    "the hint grants no authority"
                ),
            }
        )
    parent = str(issue.get("parent") or "")
    if parent:
        candidates.append(
            {
                "candidate": f"parent:{parent}",
                "source": "tracker-parent",
                "provenance_root": text_root(parent),
                "falsifier": "the parent may be only a rollup and not a domain owner",
            }
        )
    for label in sorted(str(value) for value in (issue.get("labels") or [])):
        if label.startswith(("crate:", "domain:", "authority:")):
            candidates.append(
                {
                    "candidate": label,
                    "source": "tracker-label",
                    "provenance_root": text_root(label),
                    "falsifier": (
                        "a label is a routing hint and may be stale, generic, "
                        "or unauthorized"
                    ),
                }
            )
    for path in sorted(
        str(value) for value in (issue.get("mapping") or {}).get("paths", [])
    ):
        candidates.append(
            {
                "candidate": f"path:{path}",
                "source": "source-path-clause",
                "provenance_root": text_root(path),
                "falsifier": "a mentioned path may not own the target semantics",
            }
        )
    deduplicated: dict[tuple[str, str], dict[str, Any]] = {}
    for row in candidates:
        deduplicated[(row["candidate"], row["source"])] = row
    return [deduplicated[key] for key in sorted(deduplicated)]


def v2_review_minutes(issue: Mapping[str, Any]) -> int:
    disposition = str(issue.get("disposition") or "")
    base = {
        "SECTION_NAME_ONLY": 30,
        "SUBSTANTIVE_SEMANTIC_OMISSION": 120,
        "MALFORMED_OR_WRONG_TYPE": 90,
        "ROLLUP_CHILD_SET_GAP": 150,
        "OWNER_REVIEW_REQUIRED": 120,
        "HISTORICAL_IMMUTABLE_REVIEW": 90,
    }.get(disposition, 120)
    missing = len(issue.get("missing_sections") or [])
    dependency_count = len(issue.get("dependencies") or [])
    active = bool(issue.get("assignee")) or issue.get("status") == "in_progress"
    return base + 15 * max(0, missing - 1) + min(60, 5 * dependency_count) + (
        30 if active else 0
    )


def v2_active_work_context(issue: Mapping[str, Any]) -> dict[str, Any]:
    assignee = str(issue.get("assignee") or "")
    status = str(issue.get("status") or "")
    conflict = status == "in_progress" or bool(assignee)
    document = {
        "status": status,
        "tracker_assignee": assignee,
        "conflict": conflict,
        "generated_child_desired_status": "deferred" if conflict else "open",
        "coordination_acknowledgement_required": conflict,
        "agent_mail_hint": None,
        "agent_mail_hint_authoritative": False,
        "no_claim": (
            "tracker state is deterministic operational context only; optional "
            "Agent Mail or reservation observations never enter replay roots"
        ),
    }
    return v2_rooted(document)


def v2_build_inventory(
    snapshot: LiveSnapshot,
    *,
    campaign_epoch_root: str,
    priority_filter: set[str] | None,
    status_filter: set[str] | None,
) -> dict[str, Any]:
    v2_validate_source_limits(
        row_count=len(snapshot.inventory["rows"]),
        warning_count=sum(
            len(row.get("missing_sections") or [])
            for row in snapshot.inventory["rows"]
        ),
        maximum_warnings_per_issue=max(
            (
                len(row.get("missing_sections") or [])
                for row in snapshot.inventory["rows"]
            ),
            default=0,
        ),
    )
    full_by_id = {row["id"]: row for row in snapshot.all_issues}
    rows: list[dict[str, Any]] = []
    for source_row in snapshot.inventory["rows"]:
        full_issue = full_by_id.get(source_row["id"])
        if full_issue is None:
            raise EvidenceFailed(
                f"v2 full source lacks target {source_row['id']}"
            )
        joined = {
            **source_row,
            "labels": list(full_issue.get("labels") or []),
            "dependencies": list(full_issue.get("dependencies") or []),
            "dependents": list(full_issue.get("dependents") or []),
        }
        priority_name = f"P{source_row['priority']}"
        status = str(source_row["status"])
        if priority_filter is not None and priority_name not in priority_filter:
            continue
        if status_filter is not None and status not in status_filter:
            continue
        missing_sections = sorted(source_row["missing_sections"])
        if len(missing_sections) > V2_WARNINGS_PER_ISSUE_CAP:
            raise InputRefused(
                f"{source_row['id']} exceeds {V2_WARNINGS_PER_ISSUE_CAP} warnings"
            )
        stable_root = v2_target_root(joined)
        review_minutes = v2_review_minutes(joined)
        retained_payload_bytes = len(canonical_bytes(full_issue))
        target_estimate = source_row.get("estimated_minutes")
        estimate_state = (
            "DECLARED"
            if type(target_estimate) is int and target_estimate >= 0
            else "NODATA"
        )
        target_estimate_minutes = (
            int(target_estimate) if estimate_state == "DECLARED" else None
        )
        domain_candidates = v2_domain_candidates(joined)
        active_context = v2_active_work_context(joined)
        row = {
            "id": source_row["id"],
            "issue_id": source_row["id"],
            "title": source_row["title"],
            "type": source_row["type"],
            "issue_type": source_row["type"],
            "status": status,
            "priority": source_row["priority"],
            "priority_lane": priority_name,
            "status_lane": status,
            "lane": f"{priority_name}/{status}",
            "campaign_epoch_root": campaign_epoch_root,
            "destination": (
                "history-v2.json" if status == "closed" else "review-plan-v2.json"
            ),
            "movement_destination": (
                "history-v2.json" if status == "closed" else "review-plan-v2.json"
            ),
            "tracker_assignee": source_row["assignee"],
            "coordination_assignee": source_row["assignee"],
            "tracker_owner": source_row["owner"],
            "parent": source_row["parent"],
            "labels": sorted(joined.get("labels") or []),
            "missing_sections": missing_sections,
            "disposition": source_row["disposition"],
            "disposition_rationale": source_row["rationale"],
            "disposition_falsifier": source_row["falsifier"],
            "field_roots": dict(sorted(source_row["field_roots"].items())),
            "field_byte_lengths": {
                field: len(str(full_issue.get(field) or "").encode("utf-8"))
                for field in (
                    "description",
                    "acceptance_criteria",
                    "design",
                    "notes",
                )
            },
            "clause_roots": v2_complete_clause_roots(full_issue),
            "target_root": stable_root,
            "v1_row_root": semantic_root(source_row),
            "dependency_neighborhood_root": v2_dependency_neighborhood_root(
                joined
            ),
            "dependencies": list(joined.get("dependencies") or []),
            "dependents": list(joined.get("dependents") or []),
            "domain_candidates": domain_candidates,
            "domain_candidate": domain_candidates,
            "declared_domain_owner": "",
            "declared_acceptance_owner": "",
            "implementation_owner": (
                source_row["assignee"]
                or source_row["owner"]
                or "UNRESOLVED"
            ),
            "evidence_owner": "UNRESOLVED",
            "terminal_consumer": "UNRESOLVED",
            "reviewer_provenance": {
                "state": "NODATA",
                "receipt_root": None,
            },
            "source_closure": {
                "required": bool(source_row.get("source_closure_obligation")),
                "state": (
                    "REVIEW_REQUIRED"
                    if source_row.get("source_closure_obligation")
                    else "NOT_REQUIRED_BY_V1_CLASSIFICATION"
                ),
            },
            "user_effect": (
                "make the target planning contract specific, testable, "
                "replayable, and honest without changing target behavior"
            ),
            "target_implementation_estimated_minutes": target_estimate_minutes,
            "target_implementation_estimate_minutes": target_estimate_minutes,
            "target_implementation_estimate_state": estimate_state,
            "review_minutes": review_minutes,
            "generated_child_estimated_minutes": review_minutes,
            "retained_payload_bytes": retained_payload_bytes,
            "external_authority_adapter_identity": "",
            "external_authority_receipt_root": "",
            "external_authority_verdict": "NODATA",
            "conditional_write_capability_identity": "",
            "conditional_write_receipt_root": "",
            "conditional_write_verdict": "NODATA",
            "readiness": "REVIEW_ONLY",
            "remediation_route": "ANALYSIS_ONLY",
            "active_work_context": active_context,
            "no_claim": (
                "this normalized target schedules review only; it does not "
                "approve semantics, mutate the target, or prove implementation"
            ),
        }
        rows.append(v2_rooted(row))
    rows.sort(key=lambda row: row["id"])
    nonclosed_ids = sorted(row["id"] for row in rows if row["status"] != "closed")
    history_ids = sorted(row["id"] for row in rows if row["status"] == "closed")
    warning_count = sum(len(row["missing_sections"]) for row in rows)
    v2_validate_source_limits(
        row_count=len(rows),
        warning_count=warning_count,
        maximum_warnings_per_issue=max(
            (len(row["missing_sections"]) for row in rows),
            default=0,
        ),
    )
    document = {
        "schema": V2_INVENTORY_SCHEMA,
        "v1_inventory_root": snapshot.inventory["semantic_root"],
        "campaign_epoch_root": campaign_epoch_root,
        "filters": {
            "priorities": sorted(priority_filter) if priority_filter else [],
            "statuses": sorted(status_filter) if status_filter else [],
            "empty_means_all": True,
        },
        "rows": rows,
        "counts": {
            "targets": len(rows),
            "nonclosed": len(nonclosed_ids),
            "history": len(history_ids),
            "warnings": warning_count,
        },
        "nonclosed_ids_root": semantic_root(nonclosed_ids),
        "history_ids_root": semantic_root(history_ids),
        "no_claim": (
            "v2 inventory is a tracker-read-only projection of template debt; "
            "it grants no review, mutation, completion, or product authority"
        ),
    }
    return v2_rooted(document)


def v2_validate_source_limits(
    *,
    row_count: int,
    warning_count: int,
    maximum_warnings_per_issue: int,
) -> None:
    if row_count < 0 or row_count > V2_INVENTORY_ROWS_CAP:
        raise InputRefused(
            f"v2 inventory row count exceeds {V2_INVENTORY_ROWS_CAP}"
        )
    if warning_count < 0 or warning_count > V2_WARNING_ROWS_CAP:
        raise InputRefused(
            f"v2 warning count exceeds {V2_WARNING_ROWS_CAP}"
        )
    if (
        maximum_warnings_per_issue < 0
        or maximum_warnings_per_issue > V2_WARNINGS_PER_ISSUE_CAP
    ):
        raise InputRefused(
            "v2 per-issue warning count exceeds "
            f"{V2_WARNINGS_PER_ISSUE_CAP}"
        )


def v2_campaign_epoch(snapshot: LiveSnapshot) -> tuple[str, dict[str, Any]]:
    all_issues = list(snapshot.all_issues or tuple(snapshot.issues))
    all_issues.sort(key=lambda row: row["id"])
    projection = [
        {
            "id": row["id"],
            "status": row["status"],
            "priority": row["priority"],
            "assignee": row["assignee"],
            "owner": row["owner"],
            "parent": row["parent"],
            "field_roots": {
                field: text_root(str(row.get(field) or ""))
                for field in (
                    "description",
                    "acceptance_criteria",
                    "design",
                    "notes",
                )
            },
            "dependency_neighborhood_root": semantic_root(
                {
                    "dependencies": row.get("dependencies") or [],
                    "dependents": row.get("dependents") or [],
                }
            ),
        }
        for row in all_issues
    ]
    root = semantic_root(projection)
    return root, {
        "issue_count": len(projection),
        "issue_ids_root": semantic_root([row["id"] for row in projection]),
        "projection_root": root,
    }


def v2_not_requested(
    *,
    schema: str,
    mode: str,
    source_root: str,
    inventory_root: str,
    projection: str,
) -> dict[str, Any]:
    return v2_rooted(
        {
            "schema": schema,
            "state": "NOT_REQUESTED",
            "mode": mode,
            "projection": projection,
            "source_root": source_root,
            "inventory_root": inventory_root,
            "rows": [],
            "no_claim": (
                "this projection was not selected; the rooted document is "
                "neither an empty-success claim nor evidence about its domain"
            ),
        }
    )


V2_AUDIT_EVENT_KEYS = {
    "id",
    "event_type",
    "actor",
    "timestamp",
    "old_value",
    "new_value",
    "comment",
    "agent_name",
    "harness",
    "model",
}


def v2_parse_timestamp(value: Any, *, label: str) -> datetime:
    if not isinstance(value, str) or not value:
        raise EvidenceFailed(f"{label} lacks a timestamp")
    try:
        parsed = datetime.fromisoformat(
            value[:-1] + "+00:00" if value.endswith("Z") else value
        )
    except ValueError as error:
        raise EvidenceFailed(f"{label} has a malformed timestamp") from error
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise EvidenceFailed(f"{label} timestamp lacks an explicit UTC offset")
    return parsed


def v2_normalize_audit_document(
    document: Any,
    *,
    issue_id: str,
) -> dict[str, Any]:
    if not isinstance(document, dict):
        raise InfrastructureFailed(f"audit log for {issue_id} is not an object")
    v2_exact_keys(
        document,
        {"issue_id", "events"},
        label=f"audit log {issue_id}",
    )
    if document["issue_id"] != issue_id or not isinstance(
        document["events"], list
    ):
        raise InfrastructureFailed(f"audit log identity differs for {issue_id}")
    normalized_events: list[dict[str, Any]] = []
    previous_id: int | None = None
    previous_timestamp: datetime | None = None
    seen_ids: set[int] = set()
    for index, event in enumerate(document["events"]):
        if not isinstance(event, dict):
            raise EvidenceFailed(f"audit event {issue_id}/{index} is not an object")
        unknown = set(event) - V2_AUDIT_EVENT_KEYS
        if unknown:
            raise EvidenceFailed(
                f"audit event {issue_id}/{index} has unknown fields {sorted(unknown)}"
            )
        required = {"id", "event_type", "timestamp"}
        if not required.issubset(event):
            raise EvidenceFailed(
                f"audit event {issue_id}/{index} lacks required fields"
            )
        event_id = event["id"]
        if (
            not isinstance(event_id, int)
            or isinstance(event_id, bool)
            or event_id <= 0
            or event_id in seen_ids
        ):
            raise EvidenceFailed(
                f"audit event {issue_id}/{index} has invalid or duplicate ID"
            )
        event_timestamp = v2_parse_timestamp(
            event["timestamp"],
            label=f"audit event {issue_id}/{index}",
        )
        if previous_id is not None and event_id >= previous_id:
            raise EvidenceFailed(
                f"audit event order is not strictly newest-first for {issue_id}"
            )
        if (
            previous_timestamp is not None
            and event_timestamp > previous_timestamp
        ):
            raise EvidenceFailed(
                f"audit timestamps are not non-increasing for {issue_id}"
            )
        previous_id = event_id
        previous_timestamp = event_timestamp
        seen_ids.add(event_id)
        normalized = {
            key: event.get(key)
            for key in sorted(V2_AUDIT_EVENT_KEYS)
            if key in event
        }
        if not isinstance(normalized["event_type"], str):
            raise EvidenceFailed(f"audit event type is malformed for {issue_id}")
        if "actor" in normalized and not isinstance(normalized["actor"], str):
            raise EvidenceFailed(f"audit actor is malformed for {issue_id}")
        normalized_events.append(normalized)
    return v2_rooted(
        {
            "issue_id": issue_id,
            "event_order": "newest-first",
            "events": normalized_events,
        }
    )


def v2_capture_one_audit(
    issue_id: str,
    *,
    capture_ordinal: int,
) -> tuple[str, dict[str, Any], dict[str, Any]]:
    check_cancel()
    argv = ["br", "audit", "log", issue_id, "--json"]
    environment = os.environ.copy()
    for name in SANITIZED_ENV_NAMES:
        environment.pop(name, None)
    try:
        completed = subprocess.run(
            argv,
            cwd=REPO_ROOT,
            text=False,
            capture_output=True,
            env=environment,
            timeout=CAPS["subprocess_timeout_seconds"],
            check=False,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired) as error:
        raise InfrastructureFailed(
            f"bounded audit capture failed for {issue_id}"
        ) from error
    stdout = completed.stdout
    stderr = completed.stderr
    if (
        len(stdout) > CAPS["subprocess_stdout_bytes"]
        or len(stderr) > CAPS["subprocess_stdout_bytes"]
    ):
        raise InfrastructureFailed(f"bounded audit stream cap exceeded for {issue_id}")
    receipt = {
        "capture_ordinal": capture_ordinal,
        "issue_id": issue_id,
        "argv": argv,
        "exit_code": completed.returncode,
        "result_category": (
            "SUCCESS" if completed.returncode == 0 else "UNEXPECTED_EXIT"
        ),
        "stdout_byte_length": len(stdout),
        "stdout_root": "sha256-v1:" + hashlib.sha256(stdout).hexdigest(),
        "stderr_byte_length": len(stderr),
        "stderr_root": "sha256-v1:" + hashlib.sha256(stderr).hexdigest(),
        "raw_stream_bodies_retained": False,
    }
    if completed.returncode != 0:
        raise InfrastructureFailed(
            f"br audit log returned {completed.returncode} for {issue_id}; "
            "raw diagnostic body withheld"
        )
    document = strict_json_loads(stdout, label=f"audit log {issue_id}")
    return issue_id, v2_normalize_audit_document(document, issue_id=issue_id), receipt


def v2_capture_audit_round(
    issue_ids: Sequence[str],
    *,
    capture_ordinal: int,
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    ordered_ids = sorted(issue_ids)
    with ThreadPoolExecutor(max_workers=V2_AUDIT_WORKERS) as executor:
        results = list(
            executor.map(
                lambda issue_id: v2_capture_one_audit(
                    issue_id,
                    capture_ordinal=capture_ordinal,
                ),
                ordered_ids,
            )
        )
    results.sort(key=lambda row: row[0])
    documents = {issue_id: document for issue_id, document, _ in results}
    receipts = [receipt for _, _, receipt in results]
    return documents, receipts


def v2_capture_histories(
    issue_ids: Sequence[str],
) -> dict[str, Any]:
    first, first_receipts = v2_capture_audit_round(
        issue_ids,
        capture_ordinal=1,
    )
    check_cancel()
    second, second_receipts = v2_capture_audit_round(
        issue_ids,
        capture_ordinal=2,
    )
    if first != second:
        first_ids = sorted(set(first) | set(second))
        divergence = next(
            (
                issue_id
                for issue_id in first_ids
                if first.get(issue_id) != second.get(issue_id)
            ),
            "unknown",
        )
        raise InputRefused(
            f"ConcurrentDrift: audit capture differs for {divergence}"
        )
    receipts = sorted(
        [*first_receipts, *second_receipts],
        key=lambda row: (
            int(row["capture_ordinal"]),
            str(row["issue_id"]),
            tuple(row["argv"]),
        ),
    )
    return v2_rooted(
        {
            "capture_count": 2,
            "worker_bound": V2_AUDIT_WORKERS,
            "issue_ids": sorted(issue_ids),
            "documents": [first[issue_id] for issue_id in sorted(first)],
            "command_receipts": receipts,
            "raw_stream_bodies_retained": False,
            "no_claim": (
                "audit actors are used only under the exact closer contract; "
                "self-reported agent, harness, and model fields grant no authority"
            ),
        }
    )


def v2_extract_citations(issue: Mapping[str, Any]) -> list[dict[str, Any]]:
    values: dict[str, dict[str, Any]] = {}
    pattern = re.compile(r"\bfrankensim-[A-Za-z0-9][A-Za-z0-9._-]*\b")
    for field in ("description", "acceptance_criteria", "design", "notes"):
        text = str(issue.get(field) or "")
        for match in pattern.finditer(text):
            cited = match.group(0).rstrip(".,;:)")
            values[f"{field}:{cited}"] = {
                "field": field,
                "candidate_issue_id": cited,
                "citation_root": text_root(cited),
                "authoritative": False,
            }
    return [values[key] for key in sorted(values)]


def v2_build_history(
    *,
    source_root: str,
    inventory: Mapping[str, Any],
    authority: Mapping[str, Any],
    all_issues: Sequence[Mapping[str, Any]],
    audit_capture: Mapping[str, Any],
    history_contract: Mapping[str, Any],
) -> dict[str, Any]:
    issue_by_id = {issue["id"]: issue for issue in all_issues}
    audit_by_id = {
        row["issue_id"]: row for row in audit_capture["documents"]
    }
    authority_by_id = {
        row["target_id"]: row for row in authority["decisions"]
    }
    anchor_id = str(history_contract["legacy_coverage_anchor_issue"])
    anchor_closed_at = str(
        history_contract["legacy_coverage_anchor_closed_at"]
    )
    anchor_issue = issue_by_id.get(anchor_id)
    if (
        anchor_issue is None
        or anchor_issue["closed_at"] != anchor_closed_at
        or anchor_issue["status"] != "closed"
    ):
        raise EvidenceFailed("legacy audit-coverage anchor drifted")
    anchor_audit = audit_by_id.get(anchor_id)
    if anchor_audit is None:
        raise EvidenceFailed("legacy audit-coverage anchor audit is missing")
    anchor_events = {event["id"]: event for event in anchor_audit["events"]}
    anchor_status = anchor_events.get(
        int(history_contract["legacy_coverage_anchor_status_event_id"])
    )
    anchor_close = anchor_events.get(
        int(history_contract["legacy_coverage_anchor_close_event_id"])
    )
    if (
        anchor_status is None
        or anchor_status["event_type"] != "status_changed"
        or anchor_status.get("new_value") != "closed"
        or anchor_close is None
        or anchor_close["event_type"] != "closed"
        or anchor_close.get("comment") != anchor_issue["close_reason"]
    ):
        raise EvidenceFailed("legacy audit-coverage anchor events drifted")
    anchor_timestamp = v2_parse_timestamp(
        anchor_closed_at,
        label="legacy audit-coverage anchor",
    )
    legacy_rows = [
        {
            "id": issue["id"],
            "closed_at": issue["closed_at"],
            "close_reason_root": text_root(issue["close_reason"]),
        }
        for issue in all_issues
        if issue["status"] == "closed"
        and issue["closed_at"]
        and v2_parse_timestamp(
            issue["closed_at"],
            label=f"closed issue {issue['id']}",
        )
        < anchor_timestamp
    ]
    legacy_rows.sort(key=lambda row: row["id"])
    legacy_root = semantic_root(legacy_rows)
    if (
        len(legacy_rows) != int(history_contract["legacy_coverage_count"])
        or legacy_root != history_contract["legacy_coverage_rows_root"]
    ):
        raise EvidenceFailed("legacy audit-coverage membership receipt drifted")
    legacy_ids = {row["id"] for row in legacy_rows}

    history_rows: list[dict[str, Any]] = []
    for target in inventory["rows"]:
        if target["status"] != "closed":
            continue
        issue_id = target["id"]
        issue = issue_by_id.get(issue_id)
        audit = audit_by_id.get(issue_id)
        decision = authority_by_id[issue_id]
        if issue is None or audit is None:
            raise EvidenceFailed(f"history source is incomplete for {issue_id}")
        closed_at = str(issue["closed_at"])
        close_reason = str(issue["close_reason"])
        if not closed_at or not close_reason.strip():
            raise EvidenceFailed(f"closed source lacks close metadata for {issue_id}")
        show_timestamp = v2_parse_timestamp(
            closed_at,
            label=f"closed source {issue_id}",
        )
        events = audit["events"]
        close_events = [
            event for event in events if event["event_type"] == "closed"
        ]
        transition_events = [
            event
            for event in events
            if event["event_type"] == "status_changed"
            and event.get("new_value") == "closed"
        ]
        if len(close_events) == 1 and len(transition_events) == 1:
            close_event = close_events[0]
            transition = transition_events[0]
            close_actor = str(close_event.get("actor") or "").strip()
            transition_actor = str(transition.get("actor") or "").strip()
            close_timestamp = v2_parse_timestamp(
                close_event["timestamp"],
                label=f"close audit {issue_id}",
            )
            transition_timestamp = v2_parse_timestamp(
                transition["timestamp"],
                label=f"close transition {issue_id}",
            )
            skew_ms = int(
                (close_timestamp - show_timestamp).total_seconds() * 1000
            )
            if (
                not close_actor
                or close_actor != transition_actor
                or int(close_event["id"]) <= int(transition["id"])
                or close_event.get("comment") != close_reason
                or not (
                    show_timestamp
                    <= transition_timestamp
                    <= close_timestamp
                )
                or skew_ms < 0
                or skew_ms > int(history_contract["known_pair_max_skew_ms"])
            ):
                raise EvidenceFailed(f"close audit pair conflicts for {issue_id}")
            closer_state = "KNOWN"
            actor_value: str | None = close_actor
            actor_source = "br.audit.log.closed.actor"
            audit_event_root: str | None = semantic_root(close_event)
        elif (
            not close_events
            and not transition_events
            and issue_id in legacy_ids
            and show_timestamp < anchor_timestamp
        ):
            closer_state = "LEGACY_UNAVAILABLE"
            actor_value = None
            actor_source = str(history_contract["legacy_coverage_rows_root"])
            audit_event_root = None
        else:
            raise EvidenceFailed(f"closure audit is conflicted for {issue_id}")
        if closer_state not in V2_CLOSER_STATES:
            raise EvidenceFailed("history closer state escaped its closed set")
        comments = [
            {
                "event_id": event["id"],
                "actor": str(event.get("actor") or ""),
                "timestamp": event["timestamp"],
                "comment": str(event.get("comment") or ""),
                "comment_root": text_root(str(event.get("comment") or "")),
            }
            for event in events
            if event["event_type"] in {"comment_added", "comment"}
        ]
        citations = v2_extract_citations(issue)
        candidate_consumers = sorted(
            {
                str(edge["id"])
                for edge in issue.get("dependents") or []
                if isinstance(edge, dict) and edge.get("id")
            }
            | {
                str(citation["candidate_issue_id"]) for citation in citations
            }
        )
        reviewed_consumers = (
            [
                {
                    "consumer": decision["terminal_consumer"],
                    "receipt_root": decision["reviewer_provenance"]["receipt_root"],
                }
            ]
            if decision["terminal_consumer"] != "UNRESOLVED"
            and decision["reviewer_provenance"]["receipt_root"]
            and decision["reviewer_provenance"]["reviewer"]
            and decision["reviewer_provenance"]["reviewer_kind"] != "NODATA"
            else []
        )
        row = {
            "issue_id": issue_id,
            "target_root": target["target_root"],
            "closed_at": closed_at,
            "close_reason": close_reason,
            "close_reason_root": text_root(close_reason),
            "closer_state": closer_state,
            "close_actor": actor_value,
            "close_actor_source": actor_source,
            "audit_event_root": audit_event_root,
            "audit_stream_root": audit["semantic_root"],
            "ownership": {
                "created_by": issue["created_by"],
                "assignee": issue["assignee"],
                "owner": issue["owner"],
                "closer_inference_forbidden": True,
            },
            "comments": comments,
            "notes": issue["notes"],
            "notes_root": issue["field_roots"]["notes"],
            "parent": issue["parent"],
            "dependencies": issue["dependencies"],
            "source_roots": [
                source_root,
                target["target_root"],
                audit["semantic_root"],
            ],
            "field_roots": issue["field_roots"],
            "citations": citations,
            "candidate_consumers": candidate_consumers,
            "reviewed_consumers": reviewed_consumers,
            "proof": {
                "state": "NO_IMPLEMENTATION_PROOF",
                "closed_adjudication_owner": history_contract[
                    "closed_adjudication_owner"
                ],
            },
            "no_claim": (
                "history accounts for immutable closed template debt; it does "
                "not re-adjudicate closure, infer a missing actor, or prove implementation"
            ),
        }
        history_rows.append(v2_rooted(row))
    history_rows.sort(key=lambda row: row["issue_id"])
    return v2_rooted(
        {
            "schema": V2_HISTORY_SCHEMA,
            "state": "POPULATED",
            "source_root": source_root,
            "inventory_root": inventory["semantic_root"],
            "authority_root": authority["semantic_root"],
            "audit_capture_root": audit_capture["semantic_root"],
            "legacy_coverage": {
                "anchor_issue": anchor_id,
                "anchor_closed_at": anchor_closed_at,
                "count": len(legacy_rows),
                "rows_root": legacy_root,
            },
            "rows": history_rows,
            "counts": {
                "targets": len(history_rows),
                "closer_states": dict(
                    sorted(
                        Counter(row["closer_state"] for row in history_rows).items()
                    )
                ),
            },
            "no_claim": (
                "history is immutable source accounting and candidate-consumer "
                "routing only; .6 owns any later semantic adjudication"
            ),
        }
    )


def v2_build_zero_sets(
    *,
    source_root: str,
    inventory: Mapping[str, Any],
    prior_campaign: Mapping[str, Any],
) -> dict[str, Any]:
    cells: list[dict[str, Any]] = []
    all_ids: list[str] = []
    statuses = ("open", "in_progress", "blocked", "deferred", "closed")
    for priority in range(5):
        for status in statuses:
            issue_ids = sorted(
                row["id"]
                for row in inventory["rows"]
                if row["priority"] == priority and row["status"] == status
            )
            all_ids.extend(issue_ids)
            membership = {
                "priority": priority,
                "status": status,
                "issue_ids": issue_ids,
            }
            cell = {
                **membership,
                "membership_root": semantic_root(membership),
                "zero_receipt": (
                    v2_rooted(
                        {
                            "priority": priority,
                            "status": status,
                            "source_root": source_root,
                            "inventory_root": inventory["semantic_root"],
                            "campaign_epoch_root": inventory[
                                "campaign_epoch_root"
                            ],
                            "count": 0,
                            "no_claim": (
                                "zero is exact only for this rooted campaign cell"
                            ),
                        }
                    )
                    if not issue_ids
                    else None
                ),
            }
            cells.append(v2_rooted(cell))
    expected_ids = sorted(row["id"] for row in inventory["rows"])
    if sorted(all_ids) != expected_ids or len(all_ids) != len(set(all_ids)):
        raise EvidenceFailed("all-status priority cells do not partition inventory")
    prior_rows = {
        row["id"]: row for row in prior_campaign.get("rows", [])
    }
    prior_only_ids = sorted(set(prior_rows) - set(expected_ids))
    if prior_only_ids:
        raise InputRefused(
            "prior campaign contains targets outside the current rooted scope; "
            f"first={prior_only_ids[0]}"
        )
    movements: list[dict[str, Any]] = []
    for current in inventory["rows"]:
        prior = prior_rows.get(current["id"])
        if not prior:
            continue
        before_status = str(prior["status"])
        after_status = str(current["status"])
        before_priority = int(prior["priority"])
        after_priority = int(current["priority"])
        if (
            before_status == after_status
            and before_priority == after_priority
        ):
            continue
        if before_status != "closed" and after_status == "closed":
            movement_class = "MOVED_TO_HISTORY"
        elif before_status == "closed" and after_status != "closed":
            movement_class = "MOVED_FROM_HISTORY"
        else:
            movement_class = "MOVED_TO_LANE"
        movement = {
            "issue_id": current["id"],
            "class": movement_class,
            "before_status": before_status,
            "after_status": after_status,
            "before_priority": before_priority,
            "after_priority": after_priority,
            "source_root": prior["target_root"],
            "destination_root": current["target_root"],
            "prior_campaign_root": prior_campaign["semantic_root"],
            "current_campaign_root": inventory["campaign_epoch_root"],
            "immutable_prior_evidence": prior_campaign["terminal_root"],
            "successor_lineage_root": semantic_root(
                {
                    "issue_id": current["id"],
                    "source_root": prior["target_root"],
                    "destination_root": current["target_root"],
                    "before_status": before_status,
                    "after_status": after_status,
                    "before_priority": before_priority,
                    "after_priority": after_priority,
                }
            ),
        }
        movements.append(v2_rooted(movement))
    movements.sort(key=lambda row: row["issue_id"])
    return v2_rooted(
        {
            "schema": V2_ZERO_SETS_SCHEMA,
            "source_root": source_root,
            "inventory_root": inventory["semantic_root"],
            "campaign_epoch_root": inventory["campaign_epoch_root"],
            "cells": cells,
            "partition_issue_ids_root": semantic_root(expected_ids),
            "movements": movements,
            "prior_campaign_root": prior_campaign["semantic_root"],
            "counts": {
                "cells": len(cells),
                "zero_receipts": sum(
                    cell["zero_receipt"] is not None for cell in cells
                ),
                "issues": len(expected_ids),
                "movements": len(movements),
            },
            "no_claim": (
                "zero and movement receipts are exact campaign accounting; "
                "they do not mutate targets or adjudicate closed semantics"
            ),
        }
    )


def v2_load_prior_campaign(relative: str | None) -> dict[str, Any]:
    if relative is None:
        return v2_rooted(
            {
                "state": "NOT_PROVIDED",
                "bundle_path": None,
                "terminal_root": None,
                "inventory_root": None,
                "rows": [],
                "no_claim": (
                    "movement is NoData without an explicit immutable prior campaign"
                ),
            }
        )
    relative_path = safe_relative(relative, label="prior campaign")
    if str(relative_path.parent) == ".":
        raise UsageRefused(
            "prior campaign must be nested below a safe artifact root"
        )
    artifact_root = str(relative_path.parent)
    input_dir = relative_path.name
    _, payloads, terminal, events = v2_read_retained_bundle(
        artifact_root=artifact_root,
        input_dir=input_dir,
    )
    reconstructed = v2_reconstruct_retained(
        artifact_root=artifact_root,
        input_dir=input_dir,
        payloads=payloads,
        terminal=terminal,
        events=events,
        accepted_manifest=load_case_manifest_v2(),
    )
    inventory = reconstructed[1]
    if terminal.get("terminal") != "Pass":
        raise EvidenceFailed("prior campaign lacks a Pass terminal seal")
    rows = [
        {
            "id": row["id"],
            "status": row["status"],
            "priority": row["priority"],
            "target_root": row["target_root"],
            "destination": row["destination"],
        }
        for row in inventory.get("rows", [])
        if isinstance(row, dict)
    ]
    rows.sort(key=lambda row: row["id"])
    v2_assert_unique([row["id"] for row in rows], label="prior campaign IDs")
    return v2_rooted(
        {
            "state": "PROVIDED",
            "bundle_path": str(
                relative_path
            ),
            "terminal_root": terminal["semantic_root"],
            "inventory_root": inventory["semantic_root"],
            "campaign_epoch_root": inventory["campaign_epoch_root"],
            "rows": rows,
            "no_claim": (
                "prior evidence authorizes deterministic movement accounting "
                "only and never authorizes current mutation"
            ),
        }
    )


def v2_build_source_document(
    *,
    snapshot: LiveSnapshot,
    manifest: Mapping[str, Any],
    inventory: Mapping[str, Any],
    review_receipts: Mapping[str, Any],
    prior_campaign: Mapping[str, Any],
    audit_capture: Mapping[str, Any],
    priority_filter: set[str] | None,
    status_filter: set[str] | None,
) -> dict[str, Any]:
    campaign_root, campaign_receipt = v2_campaign_epoch(snapshot)
    if campaign_root != inventory["campaign_epoch_root"]:
        raise EvidenceFailed("v2 campaign root moved during source assembly")
    captured = {
        "all_issues": list(snapshot.all_issues),
        "v1_lint_projection": snapshot.lint,
        "v1_inventory_projection": snapshot.inventory,
        "observation": snapshot.source,
        "campaign": campaign_receipt,
        "filters": {
            "priorities": sorted(priority_filter) if priority_filter else [],
            "statuses": sorted(status_filter) if status_filter else [],
            "empty_means_all": True,
        },
        "review_receipts": review_receipts,
        "prior_campaign": prior_campaign,
        "audit_capture": audit_capture,
    }
    return v2_rooted(
        {
            "schema": V2_SOURCE_SCHEMA,
            "manifest": {
                "path": str(CASE_MANIFEST_V2_REL),
                "semantic_root": manifest["semantic_root"],
                "content_identity": manifest["content_identity"],
                "case_count": manifest["case_count"],
                "assertion_count": manifest["assertion_count"],
                "criterion_count": manifest["criterion_count"],
            },
            "v1_baseline": dict(manifest["compatibility_contract"]),
            "contracts": {
                "source": dict(manifest["source_contract"]),
                "row": dict(manifest["row_contract"]),
                "authority": dict(manifest["authority_contract"]),
                "packing": dict(manifest["packing_contract"]),
                "history": dict(manifest["history_contract"]),
                "artifact": dict(manifest["artifact_contract"]),
                "logging": dict(manifest["logging_contract"]),
                "replay": dict(manifest["replay_contract"]),
                "caps": dict(manifest["caps"]),
            },
            "campaign_epoch_root": campaign_root,
            "captured": captured,
            "capture_roots": {
                "all_issues": semantic_root(captured["all_issues"]),
                "v1_lint_projection": semantic_root(
                    captured["v1_lint_projection"]
                ),
                "v1_inventory_projection": snapshot.inventory["semantic_root"],
                "observation": snapshot.source["semantic_root"],
                "review_receipts": review_receipts["semantic_root"],
                "prior_campaign": prior_campaign["semantic_root"],
                "audit_capture": audit_capture["semantic_root"],
            },
            "tracker_authority": "READ_ONLY",
            "direct_tracker_file_access": False,
            "network_access": False,
            "no_claim": (
                "source-v2 retains complete reconstruction inputs and roots; "
                "it does not authenticate reviewers, approve semantics, mutate "
                "targets, prove implementation, or authorize release"
            ),
        }
    )


def v2_snapshot_from_source(source: Mapping[str, Any]) -> LiveSnapshot:
    captured = source.get("captured")
    if not isinstance(captured, dict):
        raise EvidenceFailed("source-v2 lacks captured reconstruction inputs")
    required = {
        "all_issues",
        "v1_lint_projection",
        "v1_inventory_projection",
        "observation",
        "campaign",
        "filters",
        "review_receipts",
        "prior_campaign",
        "audit_capture",
    }
    if set(captured) != required:
        raise EvidenceFailed("source-v2 captured input schema differs")
    all_issues = captured["all_issues"]
    v1_inventory = captured["v1_inventory_projection"]
    if not isinstance(all_issues, list) or not isinstance(v1_inventory, dict):
        raise EvidenceFailed("source-v2 captured issue or inventory type differs")
    full_by_id = {
        row["id"]: row for row in all_issues if isinstance(row, dict)
    }
    selected = [
        full_by_id[row["id"]]
        for row in v1_inventory.get("rows", [])
        if row["id"] in full_by_id
    ]
    if len(selected) != len(v1_inventory.get("rows", [])):
        raise EvidenceFailed("source-v2 captured inventory loses full issue rows")
    return LiveSnapshot(
        captured["v1_lint_projection"],
        selected,
        captured["observation"],
        v1_inventory,
        {},
        tuple(all_issues),
    )


V2_REVIEW_RECEIPT_FIELDS = {
    "target_id",
    "target_root",
    "inventory_root",
    "campaign_epoch_root",
    "reviewer",
    "reviewer_kind",
    "coordination_assignee",
    "declared_domain_owner",
    "declared_acceptance_owner",
    "implementation_owner",
    "evidence_owner",
    "terminal_consumer",
    "source_closure",
    "user_effect",
    "review_minutes",
    "compatibility_key",
    "compatibility_target_ids",
    "compatibility_receipt_root",
    "compatibility_rationale",
    "compatibility_falsifier",
    "manual_authorization",
    "manual_authorization_source",
    "external_authority_adapter",
    "external_authority_receipt_root",
    "external_authority_verdict",
    "conditional_capability_identity",
    "conditional_capability_receipt_root",
    "conditional_capability_verdict",
    "gate_admission_receipt_root",
    "gate_admission_verdict",
    "no_claim",
}


def v2_empty_review_receipts(
    *, inventory_root: str, campaign_epoch_root: str
) -> dict[str, Any]:
    return v2_rooted(
        {
            "schema": V2_REVIEW_RECEIPTS_SCHEMA,
            "inventory_root": inventory_root,
            "campaign_epoch_root": campaign_epoch_root,
            "receipts": [],
            "source": "NOT_PROVIDED",
            "no_claim": (
                "missing declarations keep every target REVIEW_ONLY and "
                "ANALYSIS_ONLY"
            ),
        }
    )


def v2_load_review_receipts(
    relative: str | None,
    *,
    inventory: Mapping[str, Any],
) -> dict[str, Any]:
    if relative is None:
        return v2_empty_review_receipts(
            inventory_root=str(inventory["semantic_root"]),
            campaign_epoch_root=str(inventory["campaign_epoch_root"]),
        )
    path = resolve_safe(relative, label="review receipts", must_exist=True)
    if path.suffix != ".json":
        raise InputRefused("review receipts must be a closed JSON document")
    payload = bounded_read(path)
    document = strict_json_loads(
        payload,
        label="review receipts",
        require_canonical=True,
    )
    if not isinstance(document, dict):
        raise InputRefused("review receipts must be an object")
    v2_exact_keys(
        document,
        {
            "schema",
            "inventory_root",
            "campaign_epoch_root",
            "receipts",
            "source",
            "no_claim",
        },
        label="review receipts",
    )
    if document["schema"] != V2_REVIEW_RECEIPTS_SCHEMA:
        raise InputRefused("unknown review-receipts schema")
    if document["inventory_root"] != inventory["semantic_root"]:
        raise InputRefused("review receipts bind a different inventory root")
    if document["campaign_epoch_root"] != inventory["campaign_epoch_root"]:
        raise InputRefused("review receipts bind a different campaign epoch")
    receipts = document["receipts"]
    if not isinstance(receipts, list):
        raise InputRefused("review receipts must contain an array")
    inventory_by_id = {row["id"]: row for row in inventory["rows"]}
    normalized: list[dict[str, Any]] = []
    for index, receipt in enumerate(receipts):
        if not isinstance(receipt, dict):
            raise InputRefused(f"review receipt {index} is not an object")
        v2_exact_keys(
            receipt,
            V2_REVIEW_RECEIPT_FIELDS,
            label=f"review receipt {index}",
        )
        target_id = str(receipt["target_id"])
        if target_id not in inventory_by_id:
            raise InputRefused(f"review receipt targets unknown {target_id}")
        if receipt["target_root"] != inventory_by_id[target_id]["target_root"]:
            raise InputRefused(f"review receipt target root drifted for {target_id}")
        if receipt["inventory_root"] != inventory["semantic_root"]:
            raise InputRefused(f"review receipt inventory root drifted for {target_id}")
        if receipt["campaign_epoch_root"] != inventory["campaign_epoch_root"]:
            raise InputRefused(f"review receipt campaign root drifted for {target_id}")
        review_minutes = receipt["review_minutes"]
        if not isinstance(review_minutes, int) or review_minutes < 0:
            raise InputRefused(f"review receipt minutes are invalid for {target_id}")
        if receipt["external_authority_verdict"] not in {
            "NODATA",
            "VALID",
            "INVALID",
            "REVOKED",
            "EXPIRED",
            "CONFLICTED",
        }:
            raise InputRefused(
                f"review receipt external-authority verdict is unknown for {target_id}"
            )
        if receipt["conditional_capability_verdict"] not in {
            "NODATA",
            "VALID",
            "FAILED",
            "LYING",
            "NON_ATOMIC",
            "VERSION_MISMATCHED",
            "CONFLICTED",
        }:
            raise InputRefused(
                f"review receipt conditional-capability verdict is unknown for {target_id}"
            )
        if receipt["gate_admission_verdict"] not in {
            "NODATA",
            "VALID",
            "FAILED",
            "CONFLICTED",
        }:
            raise InputRefused(
                f"review receipt gate-admission verdict is unknown for {target_id}"
            )
        targets = receipt["compatibility_target_ids"]
        if not isinstance(targets, list) or any(
            not isinstance(value, str) for value in targets
        ):
            raise InputRefused(
                f"compatibility target IDs are malformed for {target_id}"
            )
        v2_assert_unique(targets, label=f"compatibility targets for {target_id}")
        if targets and target_id not in targets:
            raise InputRefused(
                f"compatibility set does not include its target {target_id}"
            )
        if any(value not in inventory_by_id for value in targets):
            raise InputRefused(
                f"compatibility set includes an unknown target for {target_id}"
            )
        string_fields = V2_REVIEW_RECEIPT_FIELDS - {
            "compatibility_target_ids",
            "review_minutes",
            "manual_authorization",
            "source_closure",
        }
        if any(not isinstance(receipt[field], str) for field in string_fields):
            raise InputRefused(
                f"review receipt contains a non-string scalar for {target_id}"
            )
        if not isinstance(receipt["manual_authorization"], bool):
            raise InputRefused(
                f"review receipt manual authorization is not boolean for {target_id}"
            )
        if not isinstance(receipt["source_closure"], dict):
            raise InputRefused(
                f"review receipt source closure is not an object for {target_id}"
            )
        for root_field in (
            "compatibility_receipt_root",
            "external_authority_receipt_root",
            "conditional_capability_receipt_root",
            "gate_admission_receipt_root",
        ):
            root_value = str(receipt[root_field])
            if root_value and not re.fullmatch(
                r"sha256-v1:[0-9a-f]{64}",
                root_value,
            ):
                raise InputRefused(
                    f"review receipt {root_field} is malformed for {target_id}"
                )
        normalized.append(v2_rooted(receipt))
    normalized.sort(key=lambda row: row["target_id"])
    v2_assert_unique(
        [row["target_id"] for row in normalized],
        label="review receipt target IDs",
    )
    result = dict(document)
    result["receipts"] = normalized
    result["source_identity"] = {
        "path": relative,
        "bytes": len(payload),
        "root": "sha256-v1:" + hashlib.sha256(payload).hexdigest(),
    }
    return v2_rooted(result)


def v2_validate_compatibility_receipts(
    receipts: Mapping[str, Mapping[str, Any]],
    inventory: Mapping[str, Any],
) -> None:
    inventory_by_id = {row["id"]: row for row in inventory["rows"]}
    groups: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for receipt in receipts.values():
        key = str(receipt.get("compatibility_key") or "")
        if key:
            groups[key].append(receipt)
        elif any(
            receipt.get(field)
            for field in (
                "compatibility_target_ids",
                "compatibility_receipt_root",
                "compatibility_rationale",
                "compatibility_falsifier",
            )
        ):
            raise InputRefused(
                f"target {receipt['target_id']} has compatibility data without a key"
            )
    for key, group in groups.items():
        expected_ids = sorted(str(row["target_id"]) for row in group)
        roots = {str(row.get("compatibility_receipt_root") or "") for row in group}
        rationales = {str(row.get("compatibility_rationale") or "") for row in group}
        falsifiers = {str(row.get("compatibility_falsifier") or "") for row in group}
        declared_sets = {
            tuple(sorted(str(value) for value in row["compatibility_target_ids"]))
            for row in group
        }
        if (
            "" in roots
            or "" in rationales
            or "" in falsifiers
            or len(roots) != 1
            or len(rationales) != 1
            or len(falsifiers) != 1
            or declared_sets != {tuple(expected_ids)}
        ):
            raise InputRefused(
                f"compatibility receipt {key!r} is not target-complete"
            )
        coordination = sorted(
            {
                str(row["target_id"]): str(row["coordination_assignee"])
                for row in group
            }.items()
        )
        declared_domain = sorted(
            {
                str(row["target_id"]): str(row["declared_domain_owner"])
                for row in group
            }.items()
        )
        user_effect = sorted(
            {
                str(row["target_id"]): str(row["user_effect"])
                for row in group
            }.items()
        )
        dependency_neighborhood = sorted(
            (
                target_id,
                str(inventory_by_id[target_id]["dependency_neighborhood_root"]),
            )
            for target_id in expected_ids
        )
        target_roots = sorted(
            (
                target_id,
                str(inventory_by_id[target_id]["target_root"]),
            )
            for target_id in expected_ids
        )
        computed = semantic_root(
            {
                "compatibility_key": key,
                "coordination": coordination,
                "declared_domain": declared_domain,
                "user_effect": user_effect,
                "dependency_neighborhood": dependency_neighborhood,
                "target_ids": expected_ids,
                "target_roots": target_roots,
                "rationale": next(iter(rationales)),
                "falsifier": next(iter(falsifiers)),
            }
        )
        if roots != {computed}:
            raise InputRefused(
                f"compatibility receipt {key!r} root does not match its scope"
            )


def v2_validate_receipt_bindings(
    receipt_document: Mapping[str, Any],
    inventory: Mapping[str, Any],
) -> None:
    verify_semantic_root(receipt_document, label="v2 review receipts")
    allowed_document_keys = {
        "schema",
        "inventory_root",
        "campaign_epoch_root",
        "receipts",
        "source",
        "source_identity",
        "no_claim",
        "semantic_root",
    }
    required_document_keys = allowed_document_keys - {"source_identity"}
    if (
        not required_document_keys.issubset(receipt_document)
        or not set(receipt_document).issubset(allowed_document_keys)
        or receipt_document.get("schema") != V2_REVIEW_RECEIPTS_SCHEMA
        or receipt_document.get("inventory_root") != inventory.get("semantic_root")
        or receipt_document.get("campaign_epoch_root")
        != inventory.get("campaign_epoch_root")
    ):
        raise InputRefused("review receipts bind a different source or schema")
    rows = receipt_document.get("receipts")
    if not isinstance(rows, list):
        raise InputRefused("review receipts must contain an array")
    inventory_by_id = {row["id"]: row for row in inventory["rows"]}
    observed_ids: list[str] = []
    for index, receipt in enumerate(rows):
        if not isinstance(receipt, dict):
            raise InputRefused(f"review receipt {index} is not an object")
        v2_exact_keys(
            receipt,
            {*V2_REVIEW_RECEIPT_FIELDS, "semantic_root"},
            label=f"review receipt {index}",
        )
        verify_semantic_root(receipt, label=f"review receipt {index}")
        target_id = str(receipt["target_id"])
        if target_id not in inventory_by_id:
            raise InputRefused(f"review receipt targets unknown {target_id}")
        target = inventory_by_id[target_id]
        if (
            receipt["target_root"] != target["target_root"]
            or receipt["inventory_root"] != inventory["semantic_root"]
            or receipt["campaign_epoch_root"] != inventory["campaign_epoch_root"]
        ):
            raise InputRefused(
                f"review receipt exact source roots drifted for {target_id}"
            )
        observed_ids.append(target_id)
    v2_assert_unique(observed_ids, label="review receipt target IDs")


def v2_derive_authority(
    inventory: Mapping[str, Any],
    receipt_document: Mapping[str, Any],
    *,
    current_br_version: str,
    allow_mechanical_fixture: bool = False,
) -> dict[str, Any]:
    v2_validate_receipt_bindings(receipt_document, inventory)
    receipts = {
        row["target_id"]: row for row in receipt_document.get("receipts", [])
    }
    v2_validate_compatibility_receipts(receipts, inventory)
    decisions: list[dict[str, Any]] = []
    for target in inventory["rows"]:
        receipt = receipts.get(target["id"])
        declared = bool(
            receipt
            and receipt["reviewer"]
            and receipt["reviewer_kind"]
            and receipt["declared_acceptance_owner"]
        )
        external_receipt_valid = bool(
            receipt
            and receipt["external_authority_verdict"] == "VALID"
            and receipt["external_authority_adapter"]
            and receipt["external_authority_receipt_root"]
        )
        conditional_receipt_valid = bool(
            receipt
            and receipt["conditional_capability_verdict"] == "VALID"
            and receipt["conditional_capability_identity"]
            and receipt["conditional_capability_receipt_root"]
        )
        gate_receipt_valid = bool(
            receipt
            and receipt["gate_admission_verdict"] == "VALID"
            and receipt["gate_admission_receipt_root"]
        )
        current_tool_nodata = current_br_version == "0.2.19"
        external_verified = bool(
            allow_mechanical_fixture and external_receipt_valid
        )
        conditional_verified = bool(
            allow_mechanical_fixture
            and not current_tool_nodata
            and conditional_receipt_valid
        )
        gate_admitted = bool(
            allow_mechanical_fixture
            and not current_tool_nodata
            and gate_receipt_valid
        )
        mechanically_eligible = bool(
            allow_mechanical_fixture
            and not current_tool_nodata
            and external_verified
            and conditional_verified
            and gate_admitted
        )
        if mechanically_eligible:
            readiness = "MECHANICALLY_APPLY_ELIGIBLE"
        elif declared or external_verified or conditional_verified:
            readiness = "DECLARED_READY"
        else:
            readiness = "REVIEW_ONLY"
        active_context = target["active_work_context"]
        if target["status"] == "deferred" or active_context["conflict"]:
            route = "ANALYSIS_ONLY"
            deferred_prohibition = target["status"] == "deferred"
        elif mechanically_eligible:
            route = "AUTOMATED_CONDITIONAL"
            deferred_prohibition = False
        elif declared and receipt and receipt["manual_authorization"]:
            route = "MANUAL_BR_REVIEW"
            deferred_prohibition = False
        else:
            route = "ANALYSIS_ONLY"
            deferred_prohibition = False
        decision = {
            "target_id": target["id"],
            "target_root": target["target_root"],
            "tracker_assignee": target["tracker_assignee"],
            "tracker_owner": target["tracker_owner"],
            "coordination_assignee": (
                receipt["coordination_assignee"]
                if declared and receipt
                else target["tracker_assignee"]
            ),
            "domain_candidates": target["domain_candidates"],
            "declared_domain_owner": (
                receipt["declared_domain_owner"] if declared and receipt else ""
            ),
            "declared_acceptance_owner": (
                receipt["declared_acceptance_owner"] if declared and receipt else ""
            ),
            "implementation_owner": (
                receipt["implementation_owner"]
                if declared and receipt
                else "UNRESOLVED"
            ),
            "evidence_owner": (
                receipt["evidence_owner"]
                if declared and receipt
                else "UNRESOLVED"
            ),
            "terminal_consumer": (
                receipt["terminal_consumer"]
                if declared and receipt
                else "UNRESOLVED"
            ),
            "reviewer_provenance": {
                "reviewer": receipt["reviewer"] if declared and receipt else "",
                "reviewer_kind": (
                    receipt["reviewer_kind"]
                    if declared and receipt
                    else "NODATA"
                ),
                "receipt_root": (
                    receipt["semantic_root"] if declared and receipt else None
                ),
                "self_report_ceiling": "DECLARED_READY",
            },
            "source_closure": (
                receipt["source_closure"]
                if declared and receipt
                else target["source_closure"]
            ),
            "user_effect": (
                receipt["user_effect"]
                if declared and receipt
                else target["user_effect"]
            ),
            "field_roots": target["field_roots"],
            "target_implementation_estimated_minutes": target[
                "target_implementation_estimated_minutes"
            ],
            "target_implementation_estimate_state": target[
                "target_implementation_estimate_state"
            ],
            "review_minutes": (
                receipt["review_minutes"]
                if declared and receipt
                else target["review_minutes"]
            ),
            "external_authority": {
                "adapter": (
                    receipt["external_authority_adapter"] if receipt else ""
                ),
                "receipt_root": (
                    receipt["external_authority_receipt_root"] if receipt else ""
                ),
                "verdict": (
                    receipt["external_authority_verdict"] if receipt else "NODATA"
                ),
                "verified": external_verified,
                "structurally_valid": external_receipt_valid,
                "trust_registry": (
                    "FIXTURE_INDEPENDENT"
                    if allow_mechanical_fixture
                    else "NODATA"
                ),
            },
            "conditional_write_capability": {
                "identity": (
                    receipt["conditional_capability_identity"] if receipt else ""
                ),
                "receipt_root": (
                    receipt["conditional_capability_receipt_root"]
                    if receipt
                    else ""
                ),
                "verdict": (
                    receipt["conditional_capability_verdict"]
                    if receipt
                    else "NODATA"
                ),
                "verified": conditional_verified,
                "structurally_valid": conditional_receipt_valid,
                "current_tool_verdict": (
                    "AUTOMATION_NODATA"
                    if current_tool_nodata
                    else "REQUIRES_1_3_ADMISSION"
                ),
            },
            "gate_admission": {
                "receipt_root": (
                    receipt["gate_admission_receipt_root"] if receipt else ""
                ),
                "verdict": (
                    receipt["gate_admission_verdict"] if receipt else "NODATA"
                ),
                "verified": gate_admitted,
                "structurally_valid": gate_receipt_valid,
            },
            "readiness": readiness,
            "remediation_route": route,
            "manual_authorization": {
                "declared": bool(
                    declared and receipt and receipt["manual_authorization"]
                ),
                "source": (
                    receipt["manual_authorization_source"] if receipt else ""
                ),
                "mechanical_authority": False,
                "no_cas_no_clobber_no_exactly_once": True,
            },
            "compatibility": {
                "key": (
                    receipt["compatibility_key"]
                    if declared and receipt
                    else ""
                ),
                "target_ids": (
                    receipt["compatibility_target_ids"]
                    if declared and receipt
                    else []
                ),
                "receipt_root": (
                    receipt["compatibility_receipt_root"]
                    if declared and receipt
                    else ""
                ),
                "rationale": (
                    receipt["compatibility_rationale"]
                    if declared and receipt
                    else ""
                ),
                "falsifier": (
                    receipt["compatibility_falsifier"]
                    if declared and receipt
                    else ""
                ),
                "target_complete": bool(
                    declared and receipt and receipt["compatibility_key"]
                ),
            },
            "active_work_context": active_context,
            "deferred_apply_prohibition": deferred_prohibition,
            "no_claim": (
                "authority and readiness are planning projections only; current "
                "br 0.2.19 cannot satisfy the atomic conditional-write gate"
            ),
        }
        if decision["readiness"] not in V2_READINESS_STATES:
            raise EvidenceFailed("v2 readiness escaped its closed state set")
        if decision["remediation_route"] not in V2_REMEDIATION_ROUTES:
            raise EvidenceFailed("v2 remediation route escaped its closed state set")
        decisions.append(v2_rooted(decision))
    decisions.sort(key=lambda row: row["target_id"])
    target_by_id = {row["id"]: row for row in inventory["rows"]}
    normalized_rows = []
    for decision in decisions:
        target = target_by_id[decision["target_id"]]
        normalized_rows.append(
            v2_rooted(
                {
                    "issue_id": target["id"],
                    "title": target["title"],
                    "issue_type": target["type"],
                    "priority": target["priority"],
                    "status": target["status"],
                    "missing_sections": target["missing_sections"],
                    "disposition": target["disposition"],
                    "coordination_assignee": decision[
                        "coordination_assignee"
                    ],
                    "tracker_owner": target["tracker_owner"],
                    "domain_candidate": decision["domain_candidates"],
                    "declared_domain_owner": decision[
                        "declared_domain_owner"
                    ],
                    "declared_acceptance_owner": decision[
                        "declared_acceptance_owner"
                    ],
                    "implementation_owner": decision[
                        "implementation_owner"
                    ],
                    "evidence_owner": decision["evidence_owner"],
                    "terminal_consumer": decision["terminal_consumer"],
                    "reviewer_provenance": decision["reviewer_provenance"],
                    "source_closure": decision["source_closure"],
                    "user_effect": decision["user_effect"],
                    "field_roots": decision["field_roots"],
                    "target_implementation_estimated_minutes": decision[
                        "target_implementation_estimated_minutes"
                    ],
                    "target_implementation_estimate_state": decision[
                        "target_implementation_estimate_state"
                    ],
                    "review_minutes": decision["review_minutes"],
                    "generated_child_estimated_minutes": decision[
                        "review_minutes"
                    ],
                    "external_authority_adapter_identity": decision[
                        "external_authority"
                    ]["adapter"],
                    "external_authority_receipt_root": decision[
                        "external_authority"
                    ]["receipt_root"],
                    "external_authority_verdict": decision[
                        "external_authority"
                    ]["verdict"],
                    "conditional_write_capability_identity": decision[
                        "conditional_write_capability"
                    ]["identity"],
                    "conditional_write_receipt_root": decision[
                        "conditional_write_capability"
                    ]["receipt_root"],
                    "conditional_write_verdict": decision[
                        "conditional_write_capability"
                    ]["verdict"],
                    "readiness": decision["readiness"],
                    "remediation_route": decision["remediation_route"],
                    "active_work_context": decision["active_work_context"],
                    "campaign_epoch_root": target["campaign_epoch_root"],
                    "lane": target["lane"],
                    "movement_destination": target[
                        "movement_destination"
                    ],
                    "no_claim": decision["no_claim"],
                }
            )
        )
    document = {
        "schema": V2_AUTHORITY_SCHEMA,
        "inventory_root": inventory["semantic_root"],
        "review_receipts_root": receipt_document["semantic_root"],
        "current_br_version": current_br_version,
        "current_br_conditional_write_capability": (
            "AUTOMATION_NODATA"
            if current_br_version == "0.2.19"
            else "REQUIRES_1_3_ADMISSION"
        ),
        "decisions": decisions,
        "normalized_rows": normalized_rows,
        "counts": {
            "readiness": dict(
                sorted(Counter(row["readiness"] for row in decisions).items())
            ),
            "routes": dict(
                sorted(
                    Counter(row["remediation_route"] for row in decisions).items()
                )
            ),
            "active_conflicts": sum(
                bool(row["active_work_context"]["conflict"]) for row in decisions
            ),
        },
        "no_claim": (
            "self-report, assignment, labels, and repository data cannot mint "
            "external review authority or atomic conditional-write capability"
        ),
    }
    return v2_rooted(document)


def v2_lane_parent(priority: int) -> str:
    return {
        0: "frankensim-semantic-bead-template-hygiene-961yr.2",
        1: "frankensim-semantic-bead-template-hygiene-961yr.3",
        2: "frankensim-semantic-bead-template-hygiene-961yr.4.1",
        3: "frankensim-semantic-bead-template-hygiene-961yr.4.2",
        4: "frankensim-semantic-bead-template-hygiene-961yr.4.3",
    }[priority]


def v2_hard_vector(
    target: Mapping[str, Any], authority: Mapping[str, Any]
) -> tuple[Any, ...]:
    return (
        int(target["priority"]),
        str(target["status"]),
        str(target["disposition"]),
        str(target["type"]),
        tuple(sorted(str(value) for value in target["missing_sections"])),
        str(authority["readiness"]),
        str(authority["remediation_route"]),
    )


def v2_load_tuple(target: Mapping[str, Any], authority: Mapping[str, Any]) -> tuple[int, int, int]:
    return (
        int(authority["review_minutes"]),
        int(target["retained_payload_bytes"]),
        1,
    )


def v2_bin_load(
    rows: Sequence[tuple[Mapping[str, Any], Mapping[str, Any]]]
) -> tuple[int, int, int]:
    return (
        sum(int(authority["review_minutes"]) for _, authority in rows),
        sum(int(target["retained_payload_bytes"]) for target, _ in rows),
        len(rows),
    )


def v2_bin_feasible(
    rows: Sequence[tuple[Mapping[str, Any], Mapping[str, Any]]],
    *,
    max_targets: int,
) -> bool:
    minutes, retained_bytes, targets = v2_bin_load(rows)
    return (
        targets <= max_targets
        and minutes <= V2_REVIEW_MINUTES_CAP
        and retained_bytes <= V2_CHILD_PAYLOAD_CAP
    )


def v2_partition_objective(
    bins: Sequence[Sequence[tuple[Mapping[str, Any], Mapping[str, Any]]]],
) -> tuple[Any, ...]:
    loads = [v2_bin_load(rows) for rows in bins]
    memberships = tuple(
        sorted(tuple(sorted(target["id"] for target, _ in rows)) for rows in bins)
    )
    descending_loads = tuple(sorted(loads, reverse=True))
    return (
        len(bins),
        max((load[0] for load in loads), default=0),
        max((load[1] for load in loads), default=0),
        max((load[2] for load in loads), default=0),
        descending_loads,
        memberships,
    )


def v2_objective_document(
    bins: Sequence[Sequence[tuple[Mapping[str, Any], Mapping[str, Any]]]],
) -> dict[str, Any]:
    objective = v2_partition_objective(bins)
    return {
        "child_count": objective[0],
        "maximum_review_minutes": objective[1],
        "maximum_retained_payload_bytes": objective[2],
        "maximum_target_count": objective[3],
        "descending_load_vectors": [
            list(load) for load in objective[4]
        ],
        "lexicographic_memberships": [
            list(membership) for membership in objective[5]
        ],
    }


def v2_exact_partition(
    rows: Sequence[tuple[Mapping[str, Any], Mapping[str, Any]]],
    *,
    max_targets: int,
) -> list[list[tuple[Mapping[str, Any], Mapping[str, Any]]]]:
    ordered = sorted(rows, key=lambda pair: pair[0]["id"])
    count = len(ordered)
    if count > V2_EXACT_OPTIMALITY_MAX_TARGETS:
        raise EvidenceFailed("exact partition called above its frozen target cutoff")
    feasible_masks: dict[int, list[int]] = defaultdict(list)
    for mask in range(1, 1 << count):
        selected = [ordered[index] for index in range(count) if mask & (1 << index)]
        if v2_bin_feasible(selected, max_targets=max_targets):
            pivot = (mask & -mask).bit_length() - 1
            feasible_masks[pivot].append(mask)
    memo: dict[int, list[list[tuple[Mapping[str, Any], Mapping[str, Any]]]]] = {
        0: []
    }

    def solve(mask: int) -> list[list[tuple[Mapping[str, Any], Mapping[str, Any]]]]:
        if mask in memo:
            return memo[mask]
        pivot = (mask & -mask).bit_length() - 1
        best: list[list[tuple[Mapping[str, Any], Mapping[str, Any]]]] | None = None
        for subset in feasible_masks[pivot]:
            if subset & mask != subset:
                continue
            selected = [
                ordered[index] for index in range(count) if subset & (1 << index)
            ]
            candidate = [selected, *solve(mask ^ subset)]
            if best is None or v2_partition_objective(candidate) < v2_partition_objective(best):
                best = candidate
        if best is None:
            raise EvidenceFailed("no feasible exact review partition exists")
        memo[mask] = best
        return best

    return solve((1 << count) - 1)


def v2_large_partition(
    rows: Sequence[tuple[Mapping[str, Any], Mapping[str, Any]]],
    *,
    max_targets: int,
) -> tuple[
    list[list[tuple[Mapping[str, Any], Mapping[str, Any]]]],
    dict[str, Any],
]:
    ordered = sorted(
        rows,
        key=lambda pair: (
            -int(pair[1]["review_minutes"]),
            -int(pair[0]["retained_payload_bytes"]),
            pair[0]["id"],
        ),
    )
    total_minutes = sum(int(pair[1]["review_minutes"]) for pair in ordered)
    total_bytes = sum(int(pair[0]["retained_payload_bytes"]) for pair in ordered)
    lower_bound = max(
        math.ceil(len(ordered) / max_targets),
        math.ceil(total_minutes / V2_REVIEW_MINUTES_CAP),
        math.ceil(total_bytes / V2_CHILD_PAYLOAD_CAP),
        1,
    )
    bin_count = lower_bound
    while bin_count <= len(ordered):
        bins: list[list[tuple[Mapping[str, Any], Mapping[str, Any]]]] = [
            [] for _ in range(bin_count)
        ]
        admitted = True
        for row in ordered:
            candidates = sorted(
                range(bin_count),
                key=lambda index: (
                    v2_bin_load(bins[index]),
                    tuple(target["id"] for target, _ in bins[index]),
                    index,
                ),
            )
            destination = next(
                (
                    index
                    for index in candidates
                    if v2_bin_feasible(
                        [*bins[index], row],
                        max_targets=max_targets,
                    )
                ),
                None,
            )
            if destination is None:
                admitted = False
                break
            bins[destination].append(row)
        if admitted and all(bins):
            achieved = v2_objective_document(bins)
            witness = {
                "algorithm": "DETERMINISTIC_LPT_VECTOR_PACKING",
                "lower_bound_children": lower_bound,
                "observed_children": bin_count,
                "child_count_gap": bin_count - lower_bound,
                "lower_bound_vector": {
                    "targets": math.ceil(len(ordered) / max_targets),
                    "review_minutes": math.ceil(
                        total_minutes / V2_REVIEW_MINUTES_CAP
                    ),
                    "retained_payload_bytes": math.ceil(
                        total_bytes / V2_CHILD_PAYLOAD_CAP
                    ),
                },
                "achieved_objective": achieved,
                "gap_witness": {
                    "child_count": bin_count - lower_bound,
                    "maximum_review_minutes_above_average": (
                        achieved["maximum_review_minutes"]
                        - math.ceil(total_minutes / bin_count)
                    ),
                    "maximum_retained_bytes_above_average": (
                        achieved["maximum_retained_payload_bytes"]
                        - math.ceil(total_bytes / bin_count)
                    ),
                },
                "exact_optimality_claim": False,
                "max_targets_per_child": max_targets,
                "rationale": (
                    "deterministic largest-first vector placement starts at "
                    "the independently computed multidimensional lower bound"
                ),
                "falsifier": (
                    "a feasible plan with a lexicographically smaller frozen "
                    "objective or an incorrect lower bound invalidates this witness"
                ),
            }
            return bins, v2_rooted(witness)
        bin_count += 1
    raise EvidenceFailed("large-instance review packing has no feasible partition")


def v2_pack_group(
    rows: Sequence[tuple[Mapping[str, Any], Mapping[str, Any]]],
    *,
    max_targets: int,
) -> tuple[
    list[list[tuple[Mapping[str, Any], Mapping[str, Any]]]],
    dict[str, Any],
]:
    if not 1 <= max_targets <= V2_REVIEW_TARGET_HARD_MAX:
        raise UsageRefused(
            f"--max-targets-per-child must be in 1..{V2_REVIEW_TARGET_HARD_MAX}"
        )
    if any(
        int(authority["review_minutes"]) > V2_REVIEW_MINUTES_CAP
        or int(target["retained_payload_bytes"]) > V2_CHILD_PAYLOAD_CAP
        for target, authority in rows
    ):
        raise EvidenceFailed("oversize rows must be removed before standard packing")
    if len(rows) <= V2_EXACT_OPTIMALITY_MAX_TARGETS:
        bins = v2_exact_partition(rows, max_targets=max_targets)
        witness = v2_rooted(
            {
                "algorithm": "EXACT_SUBSET_PARTITION",
                "exact_optimality_max_targets": V2_EXACT_OPTIMALITY_MAX_TARGETS,
                "exact_optimality_claim": True,
                "objective": v2_objective_document(bins),
                "max_targets_per_child": max_targets,
                "rationale": (
                    "the exact subset search minimizes the frozen objective "
                    "over every feasible partition"
                ),
                "falsifier": (
                    "a feasible partition with a lexicographically smaller "
                    "objective invalidates this witness"
                ),
            }
        )
        return bins, witness
    return v2_large_partition(rows, max_targets=max_targets)


def v2_child_text(
    targets: Sequence[Mapping[str, Any]],
    authorities: Sequence[Mapping[str, Any]],
    *,
    oversize: bool,
) -> tuple[str, str, str, str]:
    target_ids = [target["id"] for target in targets]
    disposition = str(targets[0]["disposition"])
    workflow_obligations = {
        "MALFORMED_OR_WRONG_TYPE": (
            "Validate the declared issue type and repair malformed template "
            "structure without changing the underlying scope; add focused type "
            "and schema unit tests."
        ),
        "ROLLUP_CHILD_SET_GAP": (
            "Reconstruct the exact child/dependency set, prove the rollup union "
            "and zero cells, and retain graph-aware E2E evidence."
        ),
        "OWNER_REVIEW_REQUIRED": (
            "Obtain root-bound answers from the current semantic owner for "
            "acceptance, implementation, evidence, and terminal-consumer questions."
        ),
        "SUBSTANTIVE_SEMANTIC_OMISSION": (
            "Add issue-specific acceptance and success proof obligations, "
            "including unit/property coverage, authentic E2E, detailed logging, "
            "offline replay, cancellation, and no-claim boundaries."
        ),
        "SECTION_NAME_ONLY": (
            "Preserve every existing byte of substantive prose and make only "
            "the independently reviewed heading/placement correction."
        ),
        "HISTORICAL_IMMUTABLE_REVIEW": (
            "Do not mutate the closed source; account for closer provenance, "
            "citations, candidate consumers, and proof NoData in history."
        ),
    }.get(
        disposition,
        "Resolve the exact rooted disposition without generic replacement text.",
    )
    target_lines = "\n".join(
        f"- `{target['id']}` target_root `{target['target_root']}` "
        f"missing {', '.join(target['missing_sections'])}"
        for target in targets
    )
    description = (
        "Review the exact semantic Bead-template debt for the targets below. "
        "Use each retained target and field root; do not replace issue-specific "
        "meaning with generic boilerplate.\n\n"
        f"{target_lines}\n\n"
        "This generated plan is tracker-read-only. The future publisher must "
        "re-observe every target and preserve all target semantics.\n\n"
        f"Disposition-specific workflow: {workflow_obligations}"
    )
    if oversize:
        description += (
            "\n\nThis is an OVERSIZE_REVIEW_REQUIRED escalation. The bounded "
            "index is explicitly incomplete; the retained full-plan artifact "
            "and canonical `br show ID --json` reproduction are authoritative "
            "for source retrieval."
        )
    acceptance = (
        "1. Review every exact target and field root without dropping or "
        "duplicating membership.\n"
        "2. Resolve the named disposition using issue-specific acceptance, "
        "unit/property tests, authentic E2E, deterministic logging, and replay.\n"
        "3. Preserve target implementation estimates separately from review "
        "effort and keep authority/readiness/route boundaries explicit.\n"
        "4. Refuse stale roots, active conflicts, deferred apply, semantic "
        "truncation, generic merges, and unverified evidence.\n"
        "5. Close only through the canonical generated-leaf terminal protocol."
        f"\n6. Disposition-specific acceptance: {workflow_obligations}"
    )
    design = (
        "Hard partition vector: priority/status/disposition/type/missing-set/"
        "readiness. Compatibility is valid only through the retained "
        "target-complete receipt. "
        f"Workflow `{disposition}`: {workflow_obligations}"
    )
    notes = (
        f"Targets: {', '.join(target_ids)}. Readiness: "
        f"{', '.join(sorted({row['readiness'] for row in authorities}))}. "
        "No external-human, target-implementation, science, performance, "
        "product, promotion, or release authority is minted."
    )
    return description, acceptance, design, notes


def v2_validate_child_payload(child: Mapping[str, Any]) -> None:
    verify_semantic_root(child, label="planned child")
    required = {
        "schema",
        "child_key",
        "title",
        "issue_type",
        "priority",
        "desired_status",
        "coordination_assignee",
        "labels",
        "lane_parent",
        "external_ref",
        "external_ref_key_inputs",
        "source_dependency_snapshot",
        "intended_generated_edges",
        "description_file_artifact",
        "acceptance_criteria",
        "acceptance",
        "design",
        "notes",
        "review_minutes",
        "estimated_minutes",
        "generated_child_estimated_minutes",
        "target_implementation_estimates",
        "target_implementation_estimated_minutes",
        "target_implementation_estimate_states",
        "target_ids",
        "target_roots",
        "disposition_workflow",
        "external_authority",
        "external_authority_receipt_root",
        "external_authority_verdict",
        "conditional_write_capability",
        "conditional_write_receipt_root",
        "conditional_write_verdict",
        "readiness",
        "remediation_route",
        "packing_witness_root",
        "compatibility",
        "falsifier",
        "description_root",
        "no_claim",
        "semantic_root",
    }
    v2_exact_keys(child, required, label=f"planned child {child.get('child_key')}")
    if (
        child["schema"]
        != "frankensim.beads-template-hygiene.planned-child.v2"
        or type(child["priority"]) is not int
        or child["priority"] not in range(5)
        or child["issue_type"] != "task"
        or not isinstance(child["labels"], list)
        or child["labels"] != sorted(set(child["labels"]))
        or not all(isinstance(value, str) and value for value in child["labels"])
        or child["desired_status"] not in {"open", "deferred"}
        or child["readiness"] not in V2_READINESS_STATES
        or child["remediation_route"] not in V2_REMEDIATION_ROUTES
    ):
        raise EvidenceFailed("planned child scalar or enum contract differs")
    description = child["description_file_artifact"]
    if not isinstance(description, dict):
        raise EvidenceFailed("planned child description artifact is malformed")
    v2_exact_keys(
        description,
        {"media_type", "bytes", "root", "content", "transport"},
        label="planned child description artifact",
    )
    payloads = {
        "description": description["content"].encode("utf-8"),
        "acceptance": str(child["acceptance_criteria"]).encode("utf-8"),
        "design": str(child["design"]).encode("utf-8"),
        "notes": str(child["notes"]).encode("utf-8"),
    }
    limits = {
        "description": V2_CHILD_DESCRIPTION_CAP,
        "acceptance": V2_CHILD_ACCEPTANCE_CAP,
        "design": V2_CHILD_DESIGN_CAP,
        "notes": V2_CHILD_NOTES_CAP,
    }
    for field, payload in payloads.items():
        if len(payload) > limits[field]:
            raise EvidenceFailed(f"planned child {field} exceeds its transport cap")
    if (
        description["bytes"] != len(payloads["description"])
        or description["root"] != text_root(description["content"])
        or child["description_root"] != description["root"]
        or child["acceptance"] != child["acceptance_criteria"]
    ):
        raise EvidenceFailed("planned child transport roots or aliases disagree")
    argv_bytes = sum(
        len(payloads[field]) for field in ("acceptance", "design", "notes")
    )
    if argv_bytes > V2_CHILD_ARGV_AGGREGATE_CAP:
        raise EvidenceFailed("planned child aggregate non-description argv exceeds cap")
    if len(canonical_bytes(child)) > V2_CHILD_PAYLOAD_CAP:
        raise EvidenceFailed("planned child retained payload exceeds cap")
    if (
        type(child["review_minutes"]) is not int
        or child["review_minutes"] < 0
        or child["review_minutes"] > V2_REVIEW_MINUTES_CAP
    ):
        raise EvidenceFailed("planned child review minutes exceed cap")
    target_ids = child["target_ids"]
    target_roots = child["target_roots"]
    if (
        not isinstance(target_ids, list)
        or not target_ids
        or not all(isinstance(value, str) and value for value in target_ids)
        or target_ids != sorted(target_ids)
        or len(target_ids) != len(set(target_ids))
        or not isinstance(target_roots, list)
        or len(target_roots) != len(target_ids)
        or not all(
            isinstance(value, str) and value.startswith("sha256-v1:")
            for value in target_roots
        )
    ):
        raise EvidenceFailed("planned child target membership or roots are malformed")
    estimate_rows = child["target_implementation_estimates"]
    estimate_aliases = child["target_implementation_estimated_minutes"]
    estimate_states = child["target_implementation_estimate_states"]
    if not all(
        isinstance(rows, list) and len(rows) == len(target_ids)
        for rows in (estimate_rows, estimate_aliases, estimate_states)
    ):
        raise EvidenceFailed("planned child implementation estimates are malformed")
    for target_id, estimate, alias, state_row in zip(
        target_ids,
        estimate_rows,
        estimate_aliases,
        estimate_states,
    ):
        if (
            not isinstance(estimate, dict)
            or set(estimate) != {"target_id", "state", "estimated_minutes"}
            or not isinstance(alias, dict)
            or set(alias) != {"target_id", "estimated_minutes"}
            or not isinstance(state_row, dict)
            or set(state_row) != {"target_id", "state"}
            or estimate["target_id"] != target_id
            or alias["target_id"] != target_id
            or state_row["target_id"] != target_id
            or estimate["state"] != state_row["state"]
            or estimate["state"] not in {"DECLARED", "NODATA"}
            or estimate["estimated_minutes"] != alias["estimated_minutes"]
            or (
                estimate["state"] == "DECLARED"
                and (
                    type(estimate["estimated_minutes"]) is not int
                    or estimate["estimated_minutes"] < 0
                )
            )
            or (
                estimate["state"] == "NODATA"
                and estimate["estimated_minutes"] is not None
            )
        ):
            raise EvidenceFailed(
                "planned child implementation estimate state/value disagrees"
            )
    child_key = str(child["child_key"])
    generated = f"generated:{child_key}"
    expected_edges = [
        {
            "from": generated,
            "to": child["lane_parent"],
            "type": "parent-child",
        },
        {
            "from": "frankensim-semantic-bead-template-hygiene-961yr.5",
            "to": generated,
            "type": "blocks",
        },
    ]
    if child["intended_generated_edges"] != expected_edges:
        raise EvidenceFailed(
            "planned child generated-edge membership, direction, or type differs"
        )
    key_inputs = child["external_ref_key_inputs"]
    if not isinstance(key_inputs, dict):
        raise EvidenceFailed("planned child external identity inputs are malformed")
    v2_exact_keys(
        key_inputs,
        {
            "schema",
            "priority",
            "status",
            "disposition",
            "type",
            "missing_sections",
            "readiness",
            "remediation_route",
            "target_ids",
            "target_roots",
            "authority_decision_roots",
            "compatibility_receipt_root",
        },
        label="planned child external identity inputs",
    )
    decision_roots = key_inputs["authority_decision_roots"]
    if (
        not isinstance(decision_roots, list)
        or len(decision_roots) != len(target_ids)
        or not all(isinstance(row, dict) for row in decision_roots)
        or [row["target_id"] for row in decision_roots] != target_ids
        or any(
            set(row) != {"target_id", "decision_root"}
            or not isinstance(row["decision_root"], str)
            or not row["decision_root"].startswith("sha256-v1:")
            for row in decision_roots
        )
    ):
        raise EvidenceFailed(
            "planned child authority-decision identity roots are malformed"
        )
    if (
        child["external_ref"] != f"fs-template-hygiene:v2:{child_key}"
        or key_inputs["target_ids"] != target_ids
        or key_inputs["target_roots"] != target_roots
        or semantic_root(key_inputs).split(":", 1)[1] != child_key
    ):
        raise EvidenceFailed("planned child external identity inputs disagree")


def v2_make_child(
    rows: Sequence[tuple[Mapping[str, Any], Mapping[str, Any]]],
    *,
    packing_witness: Mapping[str, Any],
    oversize: bool = False,
    oversize_artifact: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    targets = [pair[0] for pair in rows]
    authorities = [pair[1] for pair in rows]
    target_ids = sorted(target["id"] for target in targets)
    target_by_id = {target["id"]: target for target in targets}
    authority_by_id = {
        authority["target_id"]: authority for authority in authorities
    }
    ordered_targets = [target_by_id[target_id] for target_id in target_ids]
    ordered_authorities = [authority_by_id[target_id] for target_id in target_ids]
    priority = int(ordered_targets[0]["priority"])
    hard_vectors = {
        v2_hard_vector(target, authority)
        for target, authority in zip(ordered_targets, ordered_authorities)
    }
    if len(hard_vectors) != 1:
        raise EvidenceFailed("planned child crosses a hard partition vector")
    if len(ordered_targets) > 1:
        compatibility_roots = {
            authority["compatibility"]["receipt_root"]
            for authority in ordered_authorities
        }
        if "" in compatibility_roots or len(compatibility_roots) != 1:
            raise EvidenceFailed("multi-target child lacks one compatibility receipt")
    description, acceptance, design, notes = v2_child_text(
        ordered_targets,
        ordered_authorities,
        oversize=oversize,
    )
    stable_key_inputs = {
        "schema": V2_REVIEW_PLAN_SCHEMA,
        "priority": priority,
        "status": ordered_targets[0]["status"],
        "disposition": ordered_targets[0]["disposition"],
        "type": ordered_targets[0]["type"],
        "missing_sections": ordered_targets[0]["missing_sections"],
        "readiness": ordered_authorities[0]["readiness"],
        "remediation_route": ordered_authorities[0]["remediation_route"],
        "target_ids": target_ids,
        "target_roots": [target["target_root"] for target in ordered_targets],
        "authority_decision_roots": [
            {
                "target_id": authority["target_id"],
                "decision_root": authority["semantic_root"],
            }
            for authority in ordered_authorities
        ],
        "compatibility_receipt_root": (
            ordered_authorities[0]["compatibility"]["receipt_root"]
        ),
    }
    child_key = semantic_root(stable_key_inputs).split(":", 1)[1]
    actual_review_minutes = sum(
        int(authority["review_minutes"]) for authority in ordered_authorities
    )
    review_minutes = (
        min(actual_review_minutes, V2_REVIEW_MINUTES_CAP)
        if oversize
        else actual_review_minutes
    )
    desired_status = (
        "deferred"
        if any(
            target["status"] in {"blocked", "deferred", "in_progress"}
            or authority["active_work_context"]["conflict"]
            for target, authority in zip(ordered_targets, ordered_authorities)
        )
        else "open"
    )
    coordination_assignees = sorted(
        {
            str(authority["coordination_assignee"])
            for authority in ordered_authorities
            if authority["coordination_assignee"]
        }
    )
    coordination_assignee = (
        coordination_assignees[0] if len(coordination_assignees) == 1 else ""
    )
    disposition = (
        "OVERSIZE_REVIEW_REQUIRED"
        if oversize
        else ordered_targets[0]["disposition"]
    )
    title = (
        f"Review {disposition} template debt for "
        f"{target_ids[0]}"
        if len(target_ids) == 1
        else f"Review {disposition} template debt for {len(target_ids)} targets"
    )
    description_artifact = {
        "media_type": "text/markdown; charset=utf-8",
        "bytes": len(description.encode("utf-8")),
        "root": text_root(description),
        "content": description,
        "transport": "br create --description-file -",
    }
    compatibility = dict(ordered_authorities[0]["compatibility"])
    if len(ordered_targets) == 1 and not compatibility["key"]:
        compatibility["rationale"] = (
            "singleton retained because no target-complete compatibility "
            "receipt authorizes a semantic merge"
        )
        compatibility["falsifier"] = (
            "a valid target-complete receipt covering this target and another "
            "hard-compatible target permits reconsideration"
        )
    child = {
        "schema": "frankensim.beads-template-hygiene.planned-child.v2",
        "child_key": child_key,
        "title": title,
        "issue_type": "task",
        "priority": priority,
        "desired_status": desired_status,
        "coordination_assignee": coordination_assignee,
        "labels": sorted(
            {
                "generated-template-review",
                "plan-hygiene",
                "testing",
                "e2e",
                "logging",
                "replay",
                f"priority-p{priority}",
                f"readiness-{ordered_authorities[0]['readiness'].lower()}",
                (
                    "oversize-review-required"
                    if oversize
                    else "bounded-review-shard"
                ),
            }
        ),
        "lane_parent": v2_lane_parent(priority),
        "external_ref": f"fs-template-hygiene:v2:{child_key}",
        "external_ref_key_inputs": stable_key_inputs,
        "source_dependency_snapshot": [
            {
                "target_id": target["id"],
                "dependency_neighborhood_root": target[
                    "dependency_neighborhood_root"
                ],
            }
            for target in ordered_targets
        ],
        "intended_generated_edges": [
            {
                "from": f"generated:{child_key}",
                "to": v2_lane_parent(priority),
                "type": "parent-child",
            },
            {
                "from": "frankensim-semantic-bead-template-hygiene-961yr.5",
                "to": f"generated:{child_key}",
                "type": "blocks",
            },
        ],
        "description_file_artifact": description_artifact,
        "description_root": description_artifact["root"],
        "acceptance_criteria": acceptance,
        "acceptance": acceptance,
        "design": design,
        "notes": notes,
        "review_minutes": review_minutes,
        "estimated_minutes": review_minutes,
        "generated_child_estimated_minutes": review_minutes,
        "target_implementation_estimates": [
            {
                "target_id": target["id"],
                "state": target["target_implementation_estimate_state"],
                "estimated_minutes": target[
                    "target_implementation_estimated_minutes"
                ],
            }
            for target in ordered_targets
        ],
        "target_implementation_estimated_minutes": [
            {
                "target_id": target["id"],
                "estimated_minutes": target[
                    "target_implementation_estimated_minutes"
                ],
            }
            for target in ordered_targets
        ],
        "target_implementation_estimate_states": [
            {
                "target_id": target["id"],
                "state": target["target_implementation_estimate_state"],
            }
            for target in ordered_targets
        ],
        "target_ids": target_ids,
        "target_roots": [
            target_by_id[target_id]["target_root"] for target_id in target_ids
        ],
        "disposition_workflow": disposition,
        "external_authority": [
            {
                "target_id": authority["target_id"],
                **authority["external_authority"],
            }
            for authority in ordered_authorities
        ],
        "external_authority_receipt_root": [
            {
                "target_id": authority["target_id"],
                "receipt_root": authority["external_authority"]["receipt_root"],
            }
            for authority in ordered_authorities
        ],
        "external_authority_verdict": [
            {
                "target_id": authority["target_id"],
                "verdict": authority["external_authority"]["verdict"],
            }
            for authority in ordered_authorities
        ],
        "conditional_write_capability": [
            {
                "target_id": authority["target_id"],
                **authority["conditional_write_capability"],
            }
            for authority in ordered_authorities
        ],
        "conditional_write_receipt_root": [
            {
                "target_id": authority["target_id"],
                "receipt_root": authority["conditional_write_capability"][
                    "receipt_root"
                ],
            }
            for authority in ordered_authorities
        ],
        "conditional_write_verdict": [
            {
                "target_id": authority["target_id"],
                "verdict": authority["conditional_write_capability"]["verdict"],
            }
            for authority in ordered_authorities
        ],
        "readiness": ordered_authorities[0]["readiness"],
        "remediation_route": (
            "ANALYSIS_ONLY"
            if oversize
            else ordered_authorities[0]["remediation_route"]
        ),
        "packing_witness_root": packing_witness["semantic_root"],
        "compatibility": compatibility,
        "falsifier": (
            "any stale target root, missing payload field, incompatible hard "
            "vector, incomplete compatibility scope, or lost target refuses"
        ),
        "no_claim": (
            "this complete planned payload is not a published Bead and grants "
            "no target, review, conditional-write, implementation, or release authority"
        ),
    }
    if oversize and oversize_artifact is not None:
        child["notes"] += (
            f" Actual unbounded review estimate: {actual_review_minutes} minutes. "
            f"Full retained oversize artifact: {oversize_artifact['relative_path']} "
            f"root {oversize_artifact['semantic_root']}."
        )
    rooted = v2_rooted(child)
    v2_validate_child_payload(rooted)
    return rooted


def v2_validate_review_plan(
    review_plan: Mapping[str, Any],
    inventory: Mapping[str, Any],
    authority: Mapping[str, Any],
) -> None:
    verify_semantic_root(review_plan, label="v2 review plan")
    if (
        review_plan.get("schema") != V2_REVIEW_PLAN_SCHEMA
        or review_plan.get("inventory_root") != inventory.get("semantic_root")
        or review_plan.get("authority_root") != authority.get("semantic_root")
        or review_plan.get("campaign_epoch_root")
        != inventory.get("campaign_epoch_root")
    ):
        raise EvidenceFailed("v2 review-plan source bindings disagree")
    children = review_plan.get("children")
    if not isinstance(children, list):
        raise EvidenceFailed("v2 review-plan children are malformed")
    target_by_id = {
        row["id"]: row
        for row in inventory["rows"]
        if row["status"] != "closed"
    }
    authority_by_id = {
        row["target_id"]: row
        for row in authority["decisions"]
        if row["target_id"] in target_by_id
    }
    observed_ids: list[str] = []
    for child in children:
        if not isinstance(child, dict):
            raise EvidenceFailed("v2 review-plan child is malformed")
        verify_semantic_root(child, label="v2 review-plan child")
        v2_validate_child_payload(child)
        ids = child["target_ids"]
        roots = child["target_roots"]
        if any(target_id not in target_by_id for target_id in ids):
            raise EvidenceFailed("v2 review-plan child substitutes an unknown target")
        expected_roots = [target_by_id[target_id]["target_root"] for target_id in ids]
        if roots != expected_roots:
            raise EvidenceFailed("v2 review-plan child target roots disagree")
        for target_id in ids:
            target = target_by_id[target_id]
            decision = authority_by_id[target_id]
            expected_vector = v2_hard_vector(target, decision)
            observed_vector = (
                child["priority"],
                target["status"],
                child["disposition_workflow"],
                target["type"],
                tuple(target["missing_sections"]),
                child["readiness"],
                child["remediation_route"],
            )
            if child["disposition_workflow"] == "OVERSIZE_REVIEW_REQUIRED":
                observed_vector = (
                    child["priority"],
                    target["status"],
                    target["disposition"],
                    target["type"],
                    tuple(target["missing_sections"]),
                    child["readiness"],
                    child["remediation_route"],
                )
            if observed_vector != expected_vector:
                raise EvidenceFailed(
                    "v2 review-plan child crosses a hard partition vector"
                )
        observed_ids.extend(ids)
    expected_ids = sorted(target_by_id)
    if (
        sorted(observed_ids) != expected_ids
        or len(observed_ids) != len(set(observed_ids))
    ):
        raise EvidenceFailed(
            "v2 review-plan target mapping is dropped, duplicated, or substituted"
        )
    if review_plan.get("counts", {}).get("targets") != len(expected_ids):
        raise EvidenceFailed("v2 review-plan target count disagrees with membership")


def v2_build_review_plan(
    inventory: Mapping[str, Any],
    authority: Mapping[str, Any],
    *,
    max_targets: int,
) -> tuple[dict[str, Any], dict[str, bytes]]:
    target_by_id = {
        row["id"]: row for row in inventory["rows"] if row["status"] != "closed"
    }
    authority_by_id = {
        row["target_id"]: row
        for row in authority["decisions"]
        if row["target_id"] in target_by_id
    }
    groups: dict[tuple[Any, ...], list[tuple[Mapping[str, Any], Mapping[str, Any]]]] = (
        defaultdict(list)
    )
    oversize_rows: list[tuple[Mapping[str, Any], Mapping[str, Any]]] = []
    for target_id in sorted(target_by_id):
        target = target_by_id[target_id]
        decision = authority_by_id[target_id]
        if (
            int(decision["review_minutes"]) > V2_REVIEW_MINUTES_CAP
            or int(target["retained_payload_bytes"]) > V2_CHILD_PAYLOAD_CAP
            or target["field_byte_lengths"]["description"]
            > V2_CHILD_DESCRIPTION_CAP
            or target["field_byte_lengths"]["acceptance_criteria"]
            > V2_CHILD_ACCEPTANCE_CAP
            or target["field_byte_lengths"]["design"] > V2_CHILD_DESIGN_CAP
            or target["field_byte_lengths"]["notes"] > V2_CHILD_NOTES_CAP
            or any(
                int(clause["byte_length"]) > V2_CLAUSE_BYTES_CAP
                for clause in target["clause_roots"]
            )
        ):
            oversize_rows.append((target, decision))
            continue
        compatibility_key = str(decision["compatibility"]["key"] or "")
        grouping_identity = compatibility_key or f"singleton:{target_id}"
        coordination_vector = (
            decision["coordination_assignee"],
            decision["declared_domain_owner"],
            decision["declared_acceptance_owner"],
            target["user_effect"],
        )
        key = (
            *v2_hard_vector(target, decision),
            grouping_identity,
            coordination_vector,
        )
        groups[key].append((target, decision))

    children: list[dict[str, Any]] = []
    packing_witnesses: list[dict[str, Any]] = []
    for key in sorted(groups, key=lambda value: repr(value)):
        rows = groups[key]
        bins, witness = v2_pack_group(rows, max_targets=max_targets)
        packing_witnesses.append(witness)
        for batch in bins:
            children.append(
                v2_make_child(batch, packing_witness=witness)
            )

    oversize_payloads: dict[str, bytes] = {}
    oversize_inventory: list[dict[str, Any]] = []
    for target, decision in oversize_rows:
        relative = (
            "oversize/"
            f"{str(target['target_root']).split(':', 1)[-1]}.json"
        )
        complete_plan = {
            "target": target,
            "authority": decision,
            "actual_review_minutes": int(decision["review_minutes"]),
        }
        complete_plan_root = semantic_root(complete_plan)
        full = v2_rooted(
            {
                "schema": (
                    "frankensim.beads-template-hygiene.oversize-review.v2"
                ),
                "source_roots": [
                    target["target_root"],
                    target["v1_row_root"],
                    target["campaign_epoch_root"],
                ],
                "field_roots": target["field_roots"],
                "clause_roots": target["clause_roots"],
                "bounded_incomplete_index": {
                    "complete": False,
                    "target_id": target["id"],
                    "missing_sections": target["missing_sections"],
                },
                "canonical_br_show_argv": [
                    "br",
                    "show",
                    target["id"],
                    "--json",
                    *BR_READ_FLAGS,
                ],
                "retained_full_plan_path": relative,
                "retained_full_plan_root": complete_plan_root,
                "complete_plan": complete_plan,
                "owner_questions": [
                    "Who has exact acceptance authority for this target?",
                    "Which complete source clauses must the generated review retain?",
                    "What independent evidence and terminal consumer close the work?",
                ],
                "no_claim": (
                    "the bounded child index is incomplete; only this retained "
                    "full-plan artifact and the canonical target are complete sources"
                ),
            }
        )
        payload = canonical_bytes(full)
        oversize_payloads[relative] = payload
        entry = {
            "relative_path": relative,
            "media_kind": "application/json",
            "schema_kind": full["schema"],
            "byte_length": len(payload),
            "semantic_root": full["semantic_root"],
            "target_roots": [target["target_root"]],
            "field_roots": target["field_roots"],
            "clause_roots": target["clause_roots"],
            "aggregate_cap_accounting": len(payload),
        }
        oversize_inventory.append(entry)
        witness = v2_rooted(
            {
                "algorithm": "OVERSIZE_ESCALATION",
                "target_id": target["id"],
                "reason": (
                    "single target exceeds review-minute or retained-payload cap"
                ),
                "full_artifact_root": full["semantic_root"],
                "exact_optimality_claim": False,
            }
        )
        packing_witnesses.append(witness)
        children.append(
            v2_make_child(
                [(target, decision)],
                packing_witness=witness,
                oversize=True,
                oversize_artifact=entry,
            )
        )

    children.sort(key=lambda child: (child["priority"], child["child_key"]))
    mapped_ids = [target_id for child in children for target_id in child["target_ids"]]
    expected_ids = sorted(target_by_id)
    if sorted(mapped_ids) != expected_ids or len(mapped_ids) != len(set(mapped_ids)):
        raise EvidenceFailed("v2 plan does not exact-map every nonclosed target once")
    document = {
        "schema": V2_REVIEW_PLAN_SCHEMA,
        "state": "POPULATED",
        "inventory_root": inventory["semantic_root"],
        "authority_root": authority["semantic_root"],
        "campaign_epoch_root": inventory["campaign_epoch_root"],
        "max_targets_per_child": max_targets,
        "hard_max_targets_per_child": V2_REVIEW_TARGET_HARD_MAX,
        "review_minutes_cap": V2_REVIEW_MINUTES_CAP,
        "retained_child_bytes_cap": V2_CHILD_PAYLOAD_CAP,
        "children": children,
        "packing_witnesses": sorted(
            packing_witnesses, key=lambda row: row["semantic_root"]
        ),
        "oversize_content": sorted(
            oversize_inventory, key=lambda row: row["relative_path"]
        ),
        "counts": {
            "targets": len(expected_ids),
            "children": len(children),
            "oversize_targets": len(oversize_rows),
            "readiness": dict(
                sorted(Counter(child["readiness"] for child in children).items())
            ),
            "routes": dict(
                sorted(
                    Counter(child["remediation_route"] for child in children).items()
                )
            ),
        },
        "work": {
            "total_review_minutes": sum(
                int(decision["review_minutes"])
                for decision in authority_by_id.values()
            ),
            "total_bounded_child_review_minutes": sum(
                int(child["review_minutes"]) for child in children
            ),
            "total_underlying_target_review_minutes": sum(
                int(decision["review_minutes"])
                for decision in authority_by_id.values()
            ),
            "minimum_child_review_minutes": min(
                (int(child["review_minutes"]) for child in children),
                default=0,
            ),
            "maximum_child_review_minutes": max(
                (int(child["review_minutes"]) for child in children),
                default=0,
            ),
        },
        "no_claim": (
            "review-plan rows are complete publication inputs but remain "
            "tracker-read-only proposals pending .1.3 and .1.2"
        ),
    }
    rooted = v2_rooted(document)
    v2_validate_review_plan(rooted, inventory, authority)
    return rooted, oversize_payloads


def v2_reproduction_argv(
    *,
    mode: str,
    artifact_root: str,
    artifact_dir: str,
    review_receipts: str | None = None,
    prior_campaign: str | None = None,
    priorities: str | None = None,
    statuses: str | None = None,
    max_targets: int = V2_REVIEW_TARGET_DEFAULT,
    explain_target: str | None = None,
    replay_input: str | None = None,
) -> list[str]:
    if mode == "replay":
        if replay_input is None:
            raise EvidenceFailed("v2 replay reproduction lacks its input")
        return [
            str(SCRIPT_REL),
            "--replay",
            str(safe_relative(replay_input, label="replay input")),
            "--artifact-root",
            str(safe_relative(artifact_root, label="artifact root")),
            "--artifact-dir",
            str(safe_relative(artifact_dir, label="artifact dir")),
        ]
    flag = {
        "review-plan": "--review-plan",
        "history-plan": "--history-plan",
    }.get(mode)
    if flag is None:
        raise EvidenceFailed(f"unsupported v2 reproduction mode {mode}")
    argv = [
        str(SCRIPT_REL),
        flag,
        "--artifact-root",
        str(safe_relative(artifact_root, label="artifact root")),
        "--artifact-dir",
        str(safe_relative(artifact_dir, label="artifact dir")),
    ]
    for option, value in (
        ("--review-receipts", review_receipts),
        ("--prior-campaign", prior_campaign),
        ("--priorities", priorities),
        ("--statuses", statuses),
        ("--explain-target", explain_target),
    ):
        if value is not None:
            argv.extend((option, value))
    if max_targets != V2_REVIEW_TARGET_DEFAULT:
        argv.extend(("--max-targets-per-child", str(max_targets)))
    return argv


def v2_sanitized_argv(argv: Sequence[Any]) -> list[str]:
    if len(argv) > V2_COMMAND_ARGUMENTS_CAP:
        raise EvidenceFailed("v2 command receipt exceeds its argument-count cap")
    result: list[str] = []
    for raw in argv:
        value = str(raw)
        if len(value.encode("utf-8")) > V2_COMMAND_ARGUMENT_BYTES_CAP:
            raise EvidenceFailed("v2 command receipt exceeds its argument-byte cap")
        if (
            value.startswith("/")
            or str(REPO_ROOT) in value
            or re.search(
                r"(?i)(?:password|secret|credential|access[_-]?token)=",
                value,
            )
        ):
            result.append(
                "<redacted:"
                + hashlib.sha256(value.encode("utf-8")).hexdigest()
                + ">"
            )
        else:
            result.append(value)
    return result


def v2_source_command_receipts(
    source: Mapping[str, Any],
) -> list[dict[str, Any]]:
    captured = source["captured"]
    observation = captured["observation"]
    receipts: list[dict[str, Any]] = []
    for row in observation.get("command_receipts", []):
        receipts.append(
            {
                "ordering": (
                    0,
                    int(row["capture_sequence"]),
                    "",
                ),
                "argv": row["argv"],
                "exit_code": row["exit_code"],
                "result_category": row["category"],
                "stdout_byte_length": row["stdout"]["bytes"],
                "stdout_root": row["stdout"]["root"],
                "stderr_byte_length": row["stderr"]["bytes"],
                "stderr_root": row["stderr"]["root"],
            }
        )
    for row in captured["audit_capture"].get("command_receipts", []):
        receipts.append(
            {
                "ordering": (
                    int(row["capture_ordinal"]),
                    1_000_000,
                    str(row["issue_id"]),
                ),
                "argv": row["argv"],
                "exit_code": row["exit_code"],
                "result_category": row["result_category"],
                "stdout_byte_length": row["stdout_byte_length"],
                "stdout_root": row["stdout_root"],
                "stderr_byte_length": row["stderr_byte_length"],
                "stderr_root": row["stderr_root"],
            }
        )
    receipts.sort(
        key=lambda row: (
            row["ordering"],
            tuple(str(value) for value in row["argv"]),
        )
    )
    return receipts


def v2_event_rows(
    *,
    mode: str,
    subject_mode: str,
    source: Mapping[str, Any],
    inventory: Mapping[str, Any],
    authority: Mapping[str, Any],
    review_plan: Mapping[str, Any],
    history: Mapping[str, Any],
    zero_sets: Mapping[str, Any],
    safe_artifacts: Sequence[str],
    reproduction: Sequence[str],
) -> list[dict[str, Any]]:
    case_id, assertion_id = {
        "review-plan": (
            "template-lint-v2.nomock-review-plan",
            "v2.case.nomock-review-plan",
        ),
        "history-plan": (
            "template-lint-v2.nomock-history-plan",
            "v2.case.nomock-history-plan",
        ),
        "replay": (
            "template-lint-v2.nomock-offline-replay",
            "v2.case.nomock-offline-replay",
        ),
    }[mode]
    empty_root = "sha256-v1:" + hashlib.sha256(b"").hexdigest()
    safe_artifact_projection = {
        "base": list(V2_RUN_ARTIFACTS),
        "optional_count": len(safe_artifacts) - len(V2_RUN_ARTIFACTS),
        "optional_paths_root": semantic_root(
            sorted(
                set(safe_artifacts) - set(V2_RUN_ARTIFACTS),
                key=lambda value: value.encode("utf-8"),
            )
        ),
    }
    events: list[dict[str, Any]] = []

    def append(
        stage: str,
        *,
        argv: Sequence[Any] = (),
        exit_code: int = 0,
        result_category: str = "PASS",
        semantic_projection: Mapping[str, Any] | None = None,
        stdout_byte_length: int = 0,
        stdout_root: str = empty_root,
        stderr_byte_length: int = 0,
        stderr_root: str = empty_root,
        terminal: str | None = None,
    ) -> None:
        parameter = {
            "mode": mode,
            "subject_mode": subject_mode,
            "stage": stage,
            "sequence": len(events),
        }
        row = {
            "schema": V2_EVENT_SCHEMA,
            "case_id": case_id,
            "assertion_id": assertion_id,
            "executor_id": assertion_id,
            "parameter_root": semantic_root(parameter),
            "stage": stage,
            "sequence": len(events),
            "argv": v2_sanitized_argv(argv),
            "exit_code": int(exit_code),
            "result_category": result_category,
            "semantic_projection": dict(semantic_projection or {}),
            "stdout_byte_length": int(stdout_byte_length),
            "stdout_root": stdout_root,
            "stderr_byte_length": int(stderr_byte_length),
            "stderr_root": stderr_root,
            "first_divergence": None,
            "recovery": "NOT_REQUIRED",
            "terminal": terminal,
            "safe_relative_artifacts": safe_artifact_projection,
            "no_claim": (
                "the event is redacted deterministic planning evidence and "
                "mints no semantic, human, write, implementation, or release authority"
            ),
        }
        if len(events) >= V2_LOG_EVENTS_CAP:
            raise EvidenceFailed("v2 event count exceeds its cap")
        if len(canonical_bytes(row)) > V2_LOG_LINE_BYTES_CAP:
            raise EvidenceFailed("v2 event line exceeds its byte cap")
        events.append(row)

    append(
        "assertion-start",
        argv=reproduction,
        semantic_projection={
            "source_root": source["semantic_root"],
            "manifest_root": source["manifest"]["semantic_root"],
        },
    )
    for index, receipt in enumerate(v2_source_command_receipts(source)):
        append(
            f"source-command-{index:06d}",
            argv=receipt["argv"],
            exit_code=receipt["exit_code"],
            result_category=receipt["result_category"],
            semantic_projection={
                "command_receipt_root": semantic_root(
                    {key: value for key, value in receipt.items() if key != "ordering"}
                )
            },
            stdout_byte_length=receipt["stdout_byte_length"],
            stdout_root=receipt["stdout_root"],
            stderr_byte_length=receipt["stderr_byte_length"],
            stderr_root=receipt["stderr_root"],
        )
    for stage, document in (
        ("inventory-derived", inventory),
        ("authority-derived", authority),
        ("review-projection-derived", review_plan),
        ("history-projection-derived", history),
        ("zero-sets-derived", zero_sets),
    ):
        append(
            stage,
            semantic_projection={"semantic_root": document["semantic_root"]},
        )
    append(
        "assertion-terminal",
        semantic_projection={
            "source_root": source["semantic_root"],
            "inventory_root": inventory["semantic_root"],
            "authority_root": authority["semantic_root"],
            "review_plan_root": review_plan["semantic_root"],
            "history_root": history["semantic_root"],
            "zero_sets_root": zero_sets["semantic_root"],
        },
    )
    append(
        "suite-terminal",
        semantic_projection={
            "subject_mode": subject_mode,
            "prior_event_roots": [
                semantic_root(row) for row in events
            ],
        },
        terminal="Pass",
    )
    if (
        [row["sequence"] for row in events] != list(range(len(events)))
        or sum(row["terminal"] is not None for row in events) != 1
        or events[-1]["terminal"] != "Pass"
    ):
        raise EvidenceFailed("v2 terminal event is not unique and last")
    return events


def v2_artifact_identity(
    *,
    relative_path: str,
    schema_kind: str,
    payload: bytes,
) -> str:
    return semantic_root(
        {
            "relative_path": str(
                safe_relative(relative_path, label="v2 artifact")
            ),
            "schema_kind": schema_kind,
            "byte_length": len(payload),
            "content_root": (
                "sha256-v1:" + hashlib.sha256(payload).hexdigest()
            ),
        }
    )


def v2_bundle_payloads(
    *,
    mode: str,
    subject_mode: str,
    artifact_root: str,
    artifact_dir: str,
    source: Mapping[str, Any],
    inventory: Mapping[str, Any],
    authority: Mapping[str, Any],
    review_plan: Mapping[str, Any],
    history: Mapping[str, Any],
    zero_sets: Mapping[str, Any],
    optional_payloads: Mapping[str, bytes],
    reproduction: Sequence[str],
    replay_equivalence: Mapping[str, Any] | None = None,
) -> dict[str, bytes]:
    registry = review_plan.get("oversize_content") or []
    if not isinstance(registry, list):
        raise EvidenceFailed("v2 oversize registry is not an array")
    registry_by_path: dict[str, Mapping[str, Any]] = {}
    expected_registry_fields = {
        "relative_path",
        "media_kind",
        "schema_kind",
        "byte_length",
        "semantic_root",
        "target_roots",
        "field_roots",
        "clause_roots",
        "aggregate_cap_accounting",
    }
    for index, entry in enumerate(registry):
        if not isinstance(entry, dict):
            raise EvidenceFailed(f"v2 oversize registry row {index} is malformed")
        v2_exact_keys(
            entry,
            expected_registry_fields,
            label=f"v2 oversize registry row {index}",
        )
        relative = str(
            safe_relative(entry["relative_path"], label="v2 oversize artifact")
        )
        if relative in V2_RUN_ARTIFACTS or relative in registry_by_path:
            raise EvidenceFailed("v2 oversize registry duplicates an artifact")
        registry_by_path[relative] = entry
    if set(optional_payloads) != set(registry_by_path):
        raise EvidenceFailed("v2 optional payload membership differs from its registry")
    for relative, payload in optional_payloads.items():
        if len(payload) > RUN_ARTIFACT_CAP:
            raise EvidenceFailed(f"v2 optional artifact {relative} exceeds cap")
        document = strict_json_loads(
            payload,
            label=f"v2 optional artifact {relative}",
            require_canonical=True,
        )
        if not isinstance(document, dict):
            raise EvidenceFailed(f"v2 optional artifact {relative} is not an object")
        verify_semantic_root(document, label=f"v2 optional artifact {relative}")
        entry = registry_by_path[relative]
        if (
            entry["byte_length"] != len(payload)
            or entry["semantic_root"] != document["semantic_root"]
            or entry["aggregate_cap_accounting"] != len(payload)
            or entry["schema_kind"] != document.get("schema")
            or entry["media_kind"] != "application/json"
            or entry["target_roots"]
            != [document["complete_plan"]["target"]["target_root"]]
            or entry["field_roots"] != document.get("field_roots")
            or entry["clause_roots"] != document.get("clause_roots")
            or document.get("retained_full_plan_path") != relative
            or document.get("retained_full_plan_root")
            != semantic_root(document["complete_plan"])
        ):
            raise EvidenceFailed(f"v2 optional registry identity differs for {relative}")

    safe_artifacts = sorted(
        [*V2_RUN_ARTIFACTS, *optional_payloads],
        key=lambda value: value.encode("utf-8"),
    )
    events = v2_event_rows(
        mode=mode,
        subject_mode=subject_mode,
        source=source,
        inventory=inventory,
        authority=authority,
        review_plan=review_plan,
        history=history,
        zero_sets=zero_sets,
        safe_artifacts=safe_artifacts,
        reproduction=reproduction,
    )
    payloads: dict[str, bytes] = {
        "source-v2.json": canonical_bytes(source),
        "inventory-v2.json": canonical_bytes(inventory),
        "authority-v2.json": canonical_bytes(authority),
        "review-plan-v2.json": canonical_bytes(review_plan),
        "history-v2.json": canonical_bytes(history),
        "zero-sets-v2.json": canonical_bytes(zero_sets),
        "events.jsonl": b"".join(canonical_bytes(row) for row in events),
        "reproduce.txt": canonical_bytes(list(reproduction)),
        **dict(optional_payloads),
    }
    schema_by_name = {
        "source-v2.json": V2_SOURCE_SCHEMA,
        "inventory-v2.json": V2_INVENTORY_SCHEMA,
        "authority-v2.json": V2_AUTHORITY_SCHEMA,
        "review-plan-v2.json": V2_REVIEW_PLAN_SCHEMA,
        "history-v2.json": V2_HISTORY_SCHEMA,
        "zero-sets-v2.json": V2_ZERO_SETS_SCHEMA,
        "events.jsonl": V2_EVENT_SCHEMA,
        "reproduce.txt": "frankensim.argv-json.v2",
        **{
            relative: str(registry_by_path[relative]["schema_kind"])
            for relative in optional_payloads
        },
    }
    identities = {
        name: v2_artifact_identity(
            relative_path=name,
            schema_kind=schema_by_name[name],
            payload=payloads[name],
        )
        for name in sorted(payloads, key=lambda value: value.encode("utf-8"))
    }
    terminal = v2_rooted(
        {
            "schema": V2_TERMINAL_SCHEMA,
            "mode": mode,
            "subject_mode": subject_mode,
            "terminal": "Pass",
            "exit_code": 0,
            "artifact_root": str(
                safe_relative(artifact_root, label="artifact root")
            ),
            "artifact_dir": str(
                safe_relative(artifact_dir, label="artifact dir")
            ),
            "base_artifacts": list(V2_RUN_ARTIFACTS),
            "safe_relative_artifacts": safe_artifacts,
            "artifact_identities": identities,
            "optional_content_registry": registry,
            "manifest_root": source["manifest"]["semantic_root"],
            "source_root": source["semantic_root"],
            "inventory_root": inventory["semantic_root"],
            "authority_root": authority["semantic_root"],
            "review_plan_root": review_plan["semantic_root"],
            "history_root": history["semantic_root"],
            "zero_sets_root": zero_sets["semantic_root"],
            "events_content_root": (
                "sha256-v1:"
                + hashlib.sha256(payloads["events.jsonl"]).hexdigest()
            ),
            "event_count": len(events),
            "event_sequence": list(range(len(events))),
            "event_roots": [semantic_root(row) for row in events],
            "terminal_event_root": semantic_root(events[-1]),
            "reproduction": list(reproduction),
            "replay_equivalence": dict(replay_equivalence or {}),
            "first_divergence": None,
            "no_claim": (
                "the final seal proves deterministic bundle construction only; "
                "it does not prove current state on replay or mint any authority"
            ),
        }
    )
    payloads["terminal.json"] = canonical_bytes(terminal)
    if set(payloads) != set(V2_RUN_ARTIFACTS) | set(optional_payloads):
        raise EvidenceFailed("v2 payload membership differs from its exact contract")
    for name, payload in payloads.items():
        if len(payload) > RUN_ARTIFACT_CAP:
            raise EvidenceFailed(f"v2 artifact {name} exceeds its 64 MiB cap")
    return payloads


def v2_write_exclusive(path: Path, payload: bytes) -> None:
    if len(payload) > RUN_ARTIFACT_CAP:
        raise EvidenceFailed("v2 writer received an over-cap payload")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    flags |= getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags, 0o644)
    except FileExistsError as error:
        raise EvidenceFailed(
            f"v2 overwrite is forbidden for {path.relative_to(REPO_ROOT)}"
        ) from error
    try:
        with os.fdopen(descriptor, "wb", closefd=True) as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        # Repository policy forbids deletion; an incomplete prefix deliberately
        # remains without a valid terminal seal for manual recovery.
        raise


def v2_fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def v2_publish_bundle(
    *,
    artifact_root: str,
    artifact_dir: str,
    payloads: Mapping[str, bytes],
) -> dict[str, Any]:
    output = resolve_run_dir(
        artifact_root,
        artifact_dir,
        label="v2 artifact dir",
    )
    root = resolve_safe(artifact_root, label="v2 artifact root")
    if root.exists() and (root.is_symlink() or not root.is_dir()):
        raise EvidenceFailed("v2 artifact root is not a safe directory")
    if not root.exists():
        root.mkdir(parents=True, exist_ok=False)
    dir_relative = safe_relative(artifact_dir, label="v2 artifact dir")
    cursor = root
    for part in dir_relative.parts[:-1]:
        cursor = cursor / part
        if not cursor.exists():
            cursor.mkdir(exist_ok=False)
            v2_fsync_directory(cursor.parent)
        if cursor.is_symlink() or not cursor.is_dir():
            raise EvidenceFailed("v2 artifact-dir ancestor is unsafe")
    if output.parent != cursor:
        raise EvidenceFailed("v2 artifact-dir ancestry resolution disagrees")
    require_fresh_run_dir(output, label="v2 artifact run directory")
    try:
        output.mkdir(parents=False, exist_ok=False)
    except FileExistsError as error:
        raise EvidenceFailed("v2 artifact run directory already exists") from error
    v2_fsync_directory(output.parent)

    nonterminal = sorted(
        (name for name in payloads if name != "terminal.json"),
        key=lambda value: value.encode("utf-8"),
    )
    writes: list[str] = []
    touched_directories: set[Path] = {output}
    for name in nonterminal:
        check_cancel()
        relative = safe_relative(name, label="v2 artifact member")
        destination = output.joinpath(*relative.parts)
        cursor = output
        for part in relative.parts[:-1]:
            cursor = cursor / part
            if not cursor.exists():
                cursor.mkdir(exist_ok=False)
                v2_fsync_directory(cursor.parent)
            if cursor.is_symlink() or not cursor.is_dir():
                raise EvidenceFailed("v2 artifact subdirectory is unsafe")
            touched_directories.add(cursor)
        v2_write_exclusive(destination, payloads[name])
        writes.append(name)
    for directory in sorted(
        touched_directories,
        key=lambda value: len(value.parts),
        reverse=True,
    ):
        v2_fsync_directory(directory)
    check_cancel()
    v2_write_exclusive(output / "terminal.json", payloads["terminal.json"])
    writes.append("terminal.json")
    return {
        "artifact_dir": str(
            safe_relative(artifact_root, label="artifact root")
            / safe_relative(artifact_dir, label="artifact dir")
        ),
        "publication_order": writes,
        "terminal_written_last": True,
    }


def v2_synopsis(
    *,
    subject_mode: str,
    review_plan: Mapping[str, Any],
    history: Mapping[str, Any],
    authority: Mapping[str, Any],
    zero_sets: Mapping[str, Any],
    artifact_dir: str,
    reproduction: Sequence[str],
    explain_target: str | None,
) -> dict[str, Any]:
    children = review_plan.get("children") or []
    selected_ids = sorted(
        {
            target_id
            for child in children
            for target_id in child.get("target_ids", [])
        }
        | {
            row["issue_id"] for row in history.get("rows", [])
        }
    )
    if explain_target is not None and explain_target not in selected_ids:
        raise InputRefused(f"--explain-target is not selected: {explain_target}")
    shown_ids = selected_ids[:V2_SYNOPSIS_ID_PREVIEW_CAP]
    work = review_plan.get("work") or {
        "total_review_minutes": 0,
        "minimum_child_review_minutes": 0,
        "maximum_child_review_minutes": 0,
    }
    decisions = authority.get("decisions") or []
    synopsis = {
        "schema": "frankensim.beads-template-hygiene.synopsis.v2",
        "mode": subject_mode,
        "counts": {
            "targets": len(selected_ids),
            "children": len(children),
            "readiness": authority.get("counts", {}).get("readiness", {}),
            "routes": authority.get("counts", {}).get("routes", {}),
            "owners": len(
                {
                    row["coordination_assignee"]
                    for row in decisions
                    if row["coordination_assignee"]
                }
            ),
            "oversize": len(review_plan.get("oversize_content") or []),
            "zero_receipts": zero_sets.get("counts", {}).get(
                "zero_receipts", 0
            ),
        },
        "work": {
            "total_review_minutes": int(
                work.get("total_review_minutes", 0)
            ),
            "minimum_review_minutes": int(
                work.get("minimum_child_review_minutes", 0)
            ),
            "maximum_review_minutes": int(
                work.get("maximum_child_review_minutes", 0)
            ),
        },
        "capacity_wave_forecast": {
            f"P{priority}": sum(
                child.get("priority") == priority for child in children
            )
            for priority in range(5)
        },
        "active_conflicts": [
            row["target_id"]
            for row in decisions
            if row["active_work_context"]["conflict"]
        ][:V2_SYNOPSIS_ID_PREVIEW_CAP],
        "top_blockers": [
            row["target_id"]
            for row in decisions
            if row["readiness"] != "MECHANICALLY_APPLY_ELIGIBLE"
        ][:V2_SYNOPSIS_ID_PREVIEW_CAP],
        "selected_ids": {
            "total": len(selected_ids),
            "shown": len(shown_ids),
            "ids": shown_ids,
            "truncated": len(shown_ids) < len(selected_ids),
            "notice": (
                "synopsis-only ID preview truncated; machine artifacts are complete"
                if len(shown_ids) < len(selected_ids)
                else None
            ),
        },
        "why_grouped": (
            "hard priority/status/disposition/type/missing-set/readiness "
            "partitions plus target-complete compatibility receipts"
        ),
        "explanation": (
            {
                "target_id": explain_target,
                "authority_root": next(
                    row["semantic_root"]
                    for row in decisions
                    if row["target_id"] == explain_target
                ),
            }
            if explain_target is not None
            else None
        ),
        "safe_next_action": (
            "review the retained child/history payload and its authority route; "
            "do not mutate a target from this planner"
        ),
        "artifacts": artifact_dir,
        "reproduction": list(reproduction),
        "no_claim": (
            "the synopsis is bounded operator guidance; complete machine "
            "artifacts remain untruncated and no authority is minted"
        ),
    }
    if len(canonical_bytes(synopsis)) > V2_SYNOPSIS_BYTES_CAP:
        raise EvidenceFailed("v2 synopsis exceeds its frozen UTF-8 byte cap")
    return synopsis


def v2_execute_live_plan(
    *,
    mode: str,
    parsed: argparse.Namespace,
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    artifact_dir = require_v2_artifact_grammar(parsed)
    output = resolve_run_dir(
        parsed.artifact_root,
        artifact_dir,
        label="v2 artifact dir",
    )
    require_fresh_run_dir(output, label="v2 artifact run directory")
    snapshot = collect_live_v2(manifest)
    campaign_root, _ = v2_campaign_epoch(snapshot)
    priority_filter = v2_parse_csv_set(
        parsed.priorities,
        allowed={f"P{priority}" for priority in range(5)},
        label="--priorities",
    )
    status_filter = v2_parse_csv_set(
        parsed.statuses,
        allowed={"open", "in_progress", "blocked", "deferred", "closed"},
        label="--statuses",
    )
    inventory = v2_build_inventory(
        snapshot,
        campaign_epoch_root=campaign_root,
        priority_filter=priority_filter,
        status_filter=status_filter,
    )
    review_receipts = v2_load_review_receipts(
        parsed.review_receipts,
        inventory=inventory,
    )
    version = str(snapshot.source["br_version"]["version"])
    authority = v2_derive_authority(
        inventory,
        review_receipts,
        current_br_version=version,
    )
    prior_campaign = v2_load_prior_campaign(parsed.prior_campaign)
    if mode == "history-plan":
        audit_ids = sorted(
            {
                row["id"]
                for row in inventory["rows"]
                if row["status"] == "closed"
            }
            | {
                str(
                    manifest["history_contract"][
                        "legacy_coverage_anchor_issue"
                    ]
                )
            }
        )
        audit_capture = v2_capture_histories(audit_ids)
    else:
        audit_capture = v2_rooted(
            {
                "state": "NOT_REQUESTED",
                "capture_count": 0,
                "issue_ids": [],
                "documents": [],
                "command_receipts": [],
                "raw_stream_bodies_retained": False,
                "no_claim": "history audit capture was not requested",
            }
        )
    source = v2_build_source_document(
        snapshot=snapshot,
        manifest=manifest,
        inventory=inventory,
        review_receipts=review_receipts,
        prior_campaign=prior_campaign,
        audit_capture=audit_capture,
        priority_filter=priority_filter,
        status_filter=status_filter,
    )
    if mode == "review-plan":
        review_plan, optional_payloads = v2_build_review_plan(
            inventory,
            authority,
            max_targets=parsed.max_targets_per_child,
        )
        history = v2_not_requested(
            schema=V2_HISTORY_SCHEMA,
            mode=mode,
            source_root=source["semantic_root"],
            inventory_root=inventory["semantic_root"],
            projection="history",
        )
    else:
        review_plan = v2_not_requested(
            schema=V2_REVIEW_PLAN_SCHEMA,
            mode=mode,
            source_root=source["semantic_root"],
            inventory_root=inventory["semantic_root"],
            projection="review-plan",
        )
        optional_payloads = {}
        history = v2_build_history(
            source_root=source["semantic_root"],
            inventory=inventory,
            authority=authority,
            all_issues=snapshot.all_issues,
            audit_capture=audit_capture,
            history_contract=manifest["history_contract"],
        )
    zero_sets = v2_build_zero_sets(
        source_root=source["semantic_root"],
        inventory=inventory,
        prior_campaign=prior_campaign,
    )
    reproduction = v2_reproduction_argv(
        mode=mode,
        artifact_root=parsed.artifact_root,
        artifact_dir=artifact_dir,
        review_receipts=parsed.review_receipts,
        prior_campaign=parsed.prior_campaign,
        priorities=parsed.priorities,
        statuses=parsed.statuses,
        max_targets=parsed.max_targets_per_child,
        explain_target=parsed.explain_target,
    )
    payloads = v2_bundle_payloads(
        mode=mode,
        subject_mode=mode,
        artifact_root=parsed.artifact_root,
        artifact_dir=artifact_dir,
        source=source,
        inventory=inventory,
        authority=authority,
        review_plan=review_plan,
        history=history,
        zero_sets=zero_sets,
        optional_payloads=optional_payloads,
        reproduction=reproduction,
    )
    synopsis = v2_synopsis(
        subject_mode=mode,
        review_plan=review_plan,
        history=history,
        authority=authority,
        zero_sets=zero_sets,
        artifact_dir=str(
            safe_relative(parsed.artifact_root, label="artifact root")
            / safe_relative(artifact_dir, label="artifact dir")
        ),
        reproduction=reproduction,
        explain_target=parsed.explain_target,
    )
    v2_publish_bundle(
        artifact_root=parsed.artifact_root,
        artifact_dir=artifact_dir,
        payloads=payloads,
    )
    return synopsis


def v2_read_file_strict(path: Path, *, label: str) -> bytes:
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as error:
        raise InputRefused(f"{label} cannot be opened safely") from error
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise InputRefused(f"{label} is not a regular file")
        if before.st_size > RUN_ARTIFACT_CAP:
            raise InputRefused(f"{label} exceeds the v2 artifact cap")
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            chunk = os.read(descriptor, min(1024 * 1024, remaining))
            if not chunk:
                raise InputRefused(f"{label} changed while being read")
            chunks.append(chunk)
            remaining -= len(chunk)
        if os.read(descriptor, 1):
            raise InputRefused(f"{label} grew while being read")
        after = os.fstat(descriptor)
        if (
            before.st_dev,
            before.st_ino,
            before.st_size,
            before.st_mtime_ns,
        ) != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
        ):
            raise InputRefused(f"{label} changed while being read")
        return b"".join(chunks)
    finally:
        os.close(descriptor)


def v2_safe_member(value: str, *, label: str) -> str:
    try:
        return str(safe_relative(value, label=label))
    except UsageRefused as error:
        raise InputRefused(f"{label} is unsafe") from error


def v2_register_member_identity(
    seen: dict[tuple[int, int], str],
    *,
    device: int,
    inode: int,
    link_count: int,
    relative: str,
    kind: str,
) -> None:
    if kind not in {"directory", "file"}:
        raise EvidenceFailed("v2 bundle identity kind is unknown")
    identity = (device, inode)
    if identity in seen:
        raise InputRefused(
            "v2 bundle contains an inode alias between "
            f"{seen[identity]!r} and {relative!r}"
        )
    if kind == "file" and link_count != 1:
        raise InputRefused(
            f"v2 bundle file {relative!r} has an external hard-link alias"
        )
    seen[identity] = relative


def v2_enumerate_bundle(directory: Path) -> tuple[list[str], list[str]]:
    files: list[str] = []
    directories: list[str] = []
    identities: dict[tuple[int, int], str] = {}
    try:
        root_stat = os.lstat(directory)
    except OSError as error:
        raise InputRefused("v2 bundle root cannot be inspected safely") from error
    if not stat.S_ISDIR(root_stat.st_mode):
        raise InputRefused("v2 bundle root is not a directory")
    v2_register_member_identity(
        identities,
        device=root_stat.st_dev,
        inode=root_stat.st_ino,
        link_count=root_stat.st_nlink,
        relative=".",
        kind="directory",
    )
    for root_text, names, filenames in os.walk(directory, followlinks=False):
        root = Path(root_text)
        names.sort(key=lambda value: value.encode("utf-8"))
        filenames.sort(key=lambda value: value.encode("utf-8"))
        for name in list(names):
            candidate = root / name
            try:
                member_stat = os.lstat(candidate)
            except OSError as error:
                raise InputRefused(
                    "v2 bundle directory changed during enumeration"
                ) from error
            if not stat.S_ISDIR(member_stat.st_mode):
                raise InputRefused("v2 bundle contains an unsafe directory")
            relative = candidate.relative_to(directory).as_posix()
            relative = v2_safe_member(relative, label="bundle directory")
            v2_register_member_identity(
                identities,
                device=member_stat.st_dev,
                inode=member_stat.st_ino,
                link_count=member_stat.st_nlink,
                relative=relative,
                kind="directory",
            )
            directories.append(relative)
        for name in filenames:
            candidate = root / name
            try:
                member_stat = os.lstat(candidate)
            except OSError as error:
                raise InputRefused(
                    "v2 bundle file changed during enumeration"
                ) from error
            if not stat.S_ISREG(member_stat.st_mode):
                raise InputRefused("v2 bundle contains an unsafe file")
            relative = candidate.relative_to(directory).as_posix()
            relative = v2_safe_member(relative, label="bundle member")
            v2_register_member_identity(
                identities,
                device=member_stat.st_dev,
                inode=member_stat.st_ino,
                link_count=member_stat.st_nlink,
                relative=relative,
                kind="file",
            )
            files.append(relative)
    files.sort(key=lambda value: value.encode("utf-8"))
    directories.sort(key=lambda value: value.encode("utf-8"))
    v2_assert_unique(files, label="v2 bundle files")
    aliases: dict[str, str] = {}
    for relative in [*directories, *files]:
        alias = unicodedata.normalize("NFC", relative).casefold()
        if alias in aliases and aliases[alias] != relative:
            raise InputRefused("v2 bundle contains a Unicode/case path alias")
        aliases[alias] = relative
    return files, directories


def v2_read_event_stream(payload: bytes) -> list[dict[str, Any]]:
    if not payload or not payload.endswith(b"\n"):
        raise InputRefused("v2 events.jsonl lacks a terminal newline")
    lines = payload.splitlines(keepends=True)
    if len(lines) > V2_LOG_EVENTS_CAP:
        raise InputRefused("v2 events.jsonl exceeds its event-count cap")
    required = {
        "schema",
        "case_id",
        "assertion_id",
        "executor_id",
        "parameter_root",
        "stage",
        "sequence",
        "argv",
        "exit_code",
        "result_category",
        "semantic_projection",
        "stdout_byte_length",
        "stdout_root",
        "stderr_byte_length",
        "stderr_root",
        "first_divergence",
        "recovery",
        "terminal",
        "safe_relative_artifacts",
        "no_claim",
    }
    events: list[dict[str, Any]] = []
    for index, line in enumerate(lines):
        if len(line) > V2_LOG_LINE_BYTES_CAP:
            raise InputRefused(f"v2 event line {index} exceeds its cap")
        document = strict_json_loads(
            line,
            label=f"v2 event line {index}",
            require_canonical=True,
        )
        if not isinstance(document, dict):
            raise InputRefused(f"v2 event line {index} is not an object")
        v2_exact_keys(document, required, label=f"v2 event line {index}")
        if document["schema"] != V2_EVENT_SCHEMA or document["sequence"] != index:
            raise EvidenceFailed(f"v2 event sequence/schema differs at {index}")
        events.append(document)
    if (
        sum(event["terminal"] is not None for event in events) != 1
        or not events
        or events[-1]["terminal"] != "Pass"
    ):
        raise EvidenceFailed("v2 retained terminal event is not unique and last")
    return events


def v2_validate_source_document(source: Mapping[str, Any]) -> None:
    if source.get("schema") != V2_SOURCE_SCHEMA:
        raise InputRefused("v2 source artifact has an unknown schema")
    verify_semantic_root(source, label="source-v2.json")
    captured = source.get("captured")
    roots = source.get("capture_roots")
    if not isinstance(captured, dict) or not isinstance(roots, dict):
        raise EvidenceFailed("v2 source lacks captured inputs or roots")
    expected = {
        "all_issues": semantic_root(captured["all_issues"]),
        "v1_lint_projection": semantic_root(captured["v1_lint_projection"]),
        "v1_inventory_projection": captured["v1_inventory_projection"][
            "semantic_root"
        ],
        "observation": captured["observation"]["semantic_root"],
        "review_receipts": captured["review_receipts"]["semantic_root"],
        "prior_campaign": captured["prior_campaign"]["semantic_root"],
        "audit_capture": captured["audit_capture"]["semantic_root"],
    }
    if roots != expected:
        raise EvidenceFailed("v2 source captured-input roots disagree")
    for label, document in (
        ("v1 inventory", captured["v1_inventory_projection"]),
        ("observation", captured["observation"]),
        ("review receipts", captured["review_receipts"]),
        ("prior campaign", captured["prior_campaign"]),
        ("audit capture", captured["audit_capture"]),
    ):
        if not isinstance(document, dict):
            raise EvidenceFailed(f"v2 source {label} is not an object")
        verify_semantic_root(document, label=f"source-v2 {label}")
    if (
        source.get("tracker_authority") != "READ_ONLY"
        or source.get("direct_tracker_file_access") is not False
        or source.get("network_access") is not False
    ):
        raise EvidenceFailed("v2 source authority boundary differs")
    all_issues = captured["all_issues"]
    if not isinstance(all_issues, list) or any(
        not isinstance(row, dict) or not isinstance(row.get("id"), str)
        for row in all_issues
    ):
        raise EvidenceFailed("v2 source all-issue projection is malformed")
    v2_assert_unique(
        [row["id"] for row in all_issues],
        label="v2 source all-issue IDs",
    )
    filters = captured.get("filters")
    if (
        not isinstance(filters, dict)
        or set(filters) != {"priorities", "statuses", "empty_means_all"}
        or not isinstance(filters["priorities"], list)
        or not isinstance(filters["statuses"], list)
        or filters["empty_means_all"] is not True
    ):
        raise EvidenceFailed("v2 source filter projection is malformed")
    audit_capture = captured["audit_capture"]
    documents = audit_capture.get("documents")
    if not isinstance(documents, list):
        raise EvidenceFailed("v2 source audit capture is malformed")
    audit_ids: list[str] = []
    for index, document in enumerate(documents):
        if not isinstance(document, dict):
            raise EvidenceFailed(f"v2 source audit row {index} is malformed")
        issue_id = str(document.get("issue_id") or "")
        normalized = v2_normalize_audit_document(
            {"issue_id": issue_id, "events": document.get("events")},
            issue_id=issue_id,
        )
        if normalized != document:
            raise EvidenceFailed(
                f"v2 source audit row {index} differs from strict normalization"
            )
        audit_ids.append(issue_id)
    v2_assert_unique(audit_ids, label="v2 source audit IDs")


def v2_validate_source_manifest_anchor(
    source: Mapping[str, Any],
    accepted_manifest: Mapping[str, Any],
) -> None:
    retained = source.get("manifest")
    if not isinstance(retained, dict):
        raise EvidenceFailed("v2 source lacks a retained manifest identity")
    if (
        retained.get("semantic_root") != accepted_manifest.get("semantic_root")
        or retained.get("content_identity")
        != accepted_manifest.get("content_identity")
        or retained.get("case_count") != accepted_manifest.get("case_count")
        or retained.get("assertion_count")
        != accepted_manifest.get("assertion_count")
        or retained.get("criterion_count")
        != accepted_manifest.get("criterion_count")
    ):
        raise EvidenceFailed(
            "v2 retained manifest identity differs from the accepted trust anchor"
        )
    expected_contracts = {
        "source": dict(accepted_manifest["source_contract"]),
        "row": dict(accepted_manifest["row_contract"]),
        "authority": dict(accepted_manifest["authority_contract"]),
        "packing": dict(accepted_manifest["packing_contract"]),
        "history": dict(accepted_manifest["history_contract"]),
        "artifact": dict(accepted_manifest["artifact_contract"]),
        "logging": dict(accepted_manifest["logging_contract"]),
        "replay": dict(accepted_manifest["replay_contract"]),
        "caps": dict(accepted_manifest["caps"]),
    }
    if source.get("contracts") != expected_contracts:
        raise EvidenceFailed(
            "v2 retained contracts differ from the accepted manifest"
        )
    if source.get("v1_baseline") != accepted_manifest.get(
        "compatibility_contract"
    ):
        raise EvidenceFailed(
            "v2 retained v1 baseline differs from the accepted manifest"
        )


def v2_read_retained_bundle(
    *,
    artifact_root: str,
    input_dir: str,
) -> tuple[Path, dict[str, bytes], dict[str, Any], list[dict[str, Any]]]:
    directory = resolve_run_dir(
        artifact_root,
        input_dir,
        label="v2 replay input",
        must_exist=True,
    )
    if directory.is_symlink() or not directory.is_dir():
        raise InputRefused("v2 replay input is not a safe directory")
    terminal_payload = v2_read_file_strict(
        directory / "terminal.json",
        label="v2 terminal.json",
    )
    terminal = strict_json_loads(
        terminal_payload,
        label="v2 terminal.json",
        require_canonical=True,
    )
    if not isinstance(terminal, dict) or terminal.get("schema") != V2_TERMINAL_SCHEMA:
        raise InputRefused("v2 terminal.json has an unknown schema")
    verify_semantic_root(terminal, label="terminal.json")
    registry = terminal.get("optional_content_registry")
    if not isinstance(registry, list):
        raise EvidenceFailed("v2 terminal optional registry is malformed")
    optional_paths: list[str] = []
    for index, entry in enumerate(registry):
        if not isinstance(entry, dict):
            raise EvidenceFailed(f"v2 terminal optional row {index} is malformed")
        optional_paths.append(
            v2_safe_member(
                str(entry.get("relative_path") or ""),
                label=f"v2 terminal optional row {index}",
            )
        )
    v2_assert_unique(optional_paths, label="v2 terminal optional paths")
    expected_files = sorted(
        [*V2_RUN_ARTIFACTS, *optional_paths],
        key=lambda value: value.encode("utf-8"),
    )
    observed_files, observed_directories = v2_enumerate_bundle(directory)
    if observed_files != expected_files:
        missing = sorted(set(expected_files) - set(observed_files))
        extra = sorted(set(observed_files) - set(expected_files))
        first = (missing or extra or ["membership"])[0]
        raise InputRefused(
            f"v2 replay first divergence at {first}: missing={missing[:1]} extra={extra[:1]}"
        )
    expected_directory_set: set[str] = set()
    for path in optional_paths:
        parent = PurePosixPath(path).parent
        while str(parent) != ".":
            expected_directory_set.add(str(parent))
            parent = parent.parent
    expected_directories = sorted(
        expected_directory_set,
        key=lambda value: value.encode("utf-8"),
    )
    if observed_directories != expected_directories:
        raise InputRefused("v2 replay contains an unlisted directory")
    payloads = {
        name: v2_read_file_strict(
            directory.joinpath(*PurePosixPath(name).parts),
            label=f"v2 artifact {name}",
        )
        for name in observed_files
    }
    payloads["terminal.json"] = terminal_payload
    events = v2_read_event_stream(payloads["events.jsonl"])
    return directory, payloads, terminal, events


def v2_reconstruct_retained(
    *,
    artifact_root: str,
    input_dir: str,
    payloads: Mapping[str, bytes],
    terminal: Mapping[str, Any],
    events: Sequence[Mapping[str, Any]],
    accepted_manifest: Mapping[str, Any] | None = None,
) -> tuple[
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, bytes],
]:
    documents: dict[str, dict[str, Any]] = {}
    schemas = {
        "source-v2.json": V2_SOURCE_SCHEMA,
        "inventory-v2.json": V2_INVENTORY_SCHEMA,
        "authority-v2.json": V2_AUTHORITY_SCHEMA,
        "review-plan-v2.json": V2_REVIEW_PLAN_SCHEMA,
        "history-v2.json": V2_HISTORY_SCHEMA,
        "zero-sets-v2.json": V2_ZERO_SETS_SCHEMA,
    }
    for name, schema in schemas.items():
        document = strict_json_loads(
            payloads[name],
            label=f"v2 artifact {name}",
            require_canonical=True,
        )
        if not isinstance(document, dict) or document.get("schema") != schema:
            raise InputRefused(f"v2 artifact {name} has an unknown schema")
        verify_semantic_root(document, label=name)
        documents[name] = document
    source = documents["source-v2.json"]
    retained_inventory = documents["inventory-v2.json"]
    retained_authority = documents["authority-v2.json"]
    retained_review = documents["review-plan-v2.json"]
    retained_history = documents["history-v2.json"]
    retained_zero = documents["zero-sets-v2.json"]
    v2_validate_source_document(source)
    accepted = (
        accepted_manifest
        if accepted_manifest is not None
        else load_case_manifest_v2()
    )
    v2_validate_source_manifest_anchor(source, accepted)
    reproduction = strict_json_loads(
        payloads["reproduce.txt"],
        label="v2 reproduce.txt",
        require_canonical=True,
    )
    if not isinstance(reproduction, list) or not all(
        isinstance(value, str) for value in reproduction
    ):
        raise EvidenceFailed("v2 reproduce.txt is not an argv array")

    roots = {
        "manifest_root": source["manifest"]["semantic_root"],
        "source_root": source["semantic_root"],
        "inventory_root": retained_inventory["semantic_root"],
        "authority_root": retained_authority["semantic_root"],
        "review_plan_root": retained_review["semantic_root"],
        "history_root": retained_history["semantic_root"],
        "zero_sets_root": retained_zero["semantic_root"],
    }
    for field, expected in roots.items():
        if terminal.get(field) != expected:
            raise EvidenceFailed(f"v2 terminal root differs at {field}")
    if (
        terminal.get("event_count") != len(events)
        or terminal.get("event_sequence") != list(range(len(events)))
        or terminal.get("event_roots")
        != [semantic_root(row) for row in events]
        or terminal.get("terminal_event_root") != semantic_root(events[-1])
        or terminal.get("events_content_root")
        != "sha256-v1:" + hashlib.sha256(payloads["events.jsonl"]).hexdigest()
        or terminal.get("reproduction") != reproduction
    ):
        raise EvidenceFailed("v2 terminal event/reproduction seal differs")

    registry = {
        row["relative_path"]: row
        for row in terminal["optional_content_registry"]
    }
    schema_by_name = {
        **schemas,
        "events.jsonl": V2_EVENT_SCHEMA,
        "reproduce.txt": "frankensim.argv-json.v2",
        **{
            path: row["schema_kind"] for path, row in registry.items()
        },
    }
    expected_identity_names = set(payloads) - {"terminal.json"}
    if set(terminal.get("artifact_identities", {})) != expected_identity_names:
        raise EvidenceFailed("v2 terminal artifact identity membership differs")
    for name in sorted(expected_identity_names):
        expected_identity = v2_artifact_identity(
            relative_path=name,
            schema_kind=schema_by_name[name],
            payload=payloads[name],
        )
        if terminal["artifact_identities"][name] != expected_identity:
            raise EvidenceFailed(f"v2 artifact identity differs for {name}")

    snapshot = v2_snapshot_from_source(source)
    campaign_root, _ = v2_campaign_epoch(snapshot)
    if campaign_root != source["campaign_epoch_root"]:
        raise EvidenceFailed("v2 retained campaign root cannot be reconstructed")
    filters = source["captured"]["filters"]
    priority_filter = set(filters["priorities"]) or None
    status_filter = set(filters["statuses"]) or None
    inventory = v2_build_inventory(
        snapshot,
        campaign_epoch_root=campaign_root,
        priority_filter=priority_filter,
        status_filter=status_filter,
    )
    authority = v2_derive_authority(
        inventory,
        source["captured"]["review_receipts"],
        current_br_version=str(
            source["captured"]["observation"]["br_version"]["version"]
        ),
    )
    subject_mode = str(terminal.get("subject_mode") or terminal.get("mode") or "")
    if subject_mode == "review-plan":
        review, optional_payloads = v2_build_review_plan(
            inventory,
            authority,
            max_targets=int(retained_review["max_targets_per_child"]),
        )
        history = v2_not_requested(
            schema=V2_HISTORY_SCHEMA,
            mode=subject_mode,
            source_root=source["semantic_root"],
            inventory_root=inventory["semantic_root"],
            projection="history",
        )
    elif subject_mode == "history-plan":
        review = v2_not_requested(
            schema=V2_REVIEW_PLAN_SCHEMA,
            mode=subject_mode,
            source_root=source["semantic_root"],
            inventory_root=inventory["semantic_root"],
            projection="review-plan",
        )
        optional_payloads = {}
        history = v2_build_history(
            source_root=source["semantic_root"],
            inventory=inventory,
            authority=authority,
            all_issues=snapshot.all_issues,
            audit_capture=source["captured"]["audit_capture"],
            history_contract=source["contracts"]["history"],
        )
    else:
        raise EvidenceFailed("v2 terminal subject mode is unknown")
    zero_sets = v2_build_zero_sets(
        source_root=source["semantic_root"],
        inventory=inventory,
        prior_campaign=source["captured"]["prior_campaign"],
    )
    reconstructed = v2_bundle_payloads(
        mode=str(terminal["mode"]),
        subject_mode=subject_mode,
        artifact_root=artifact_root,
        artifact_dir=input_dir,
        source=source,
        inventory=inventory,
        authority=authority,
        review_plan=review,
        history=history,
        zero_sets=zero_sets,
        optional_payloads=optional_payloads,
        reproduction=reproduction,
        replay_equivalence=terminal.get("replay_equivalence") or {},
    )
    if set(reconstructed) != set(payloads):
        raise EvidenceFailed("v2 reconstructed artifact membership differs")
    for name in sorted(reconstructed):
        if reconstructed[name] != payloads[name]:
            raise EvidenceFailed(f"v2 replay first byte divergence in {name}")
    return (
        source,
        inventory,
        authority,
        review,
        history,
        zero_sets,
        optional_payloads,
    )


def v2_replay_bundle(
    *,
    artifact_root: str,
    input_dir: str,
    output_dir: str,
) -> dict[str, Any]:
    source_directory, retained_payloads, terminal, events = (
        v2_read_retained_bundle(
            artifact_root=artifact_root,
            input_dir=input_dir,
        )
    )
    output_directory = resolve_run_dir(
        artifact_root,
        output_dir,
        label="v2 replay output",
    )
    source_resolved = source_directory.resolve(strict=True)
    output_resolved = output_directory.resolve(strict=False)
    if (
        source_resolved == output_resolved
        or source_resolved in output_resolved.parents
        or output_resolved in source_resolved.parents
    ):
        raise InputRefused("v2 replay input and output must be disjoint")
    (
        source,
        inventory,
        authority,
        review,
        history,
        zero_sets,
        optional_payloads,
    ) = v2_reconstruct_retained(
        artifact_root=artifact_root,
        input_dir=input_dir,
        payloads=retained_payloads,
        terminal=terminal,
        events=events,
        accepted_manifest=load_case_manifest_v2(),
    )
    reproduction = v2_reproduction_argv(
        mode="replay",
        artifact_root=artifact_root,
        artifact_dir=output_dir,
        replay_input=input_dir,
    )
    replay_equivalence = {
        "retained_terminal_root": terminal["semantic_root"],
        "retained_source_root": source["semantic_root"],
        "retained_inventory_root": inventory["semantic_root"],
        "retained_authority_root": authority["semantic_root"],
        "retained_review_plan_root": review["semantic_root"],
        "retained_history_root": history["semantic_root"],
        "retained_zero_sets_root": zero_sets["semantic_root"],
        "retained_events_content_root": terminal["events_content_root"],
        "live_tracker_reads": False,
        "live_tracker_writes": False,
        "network_access": False,
    }
    replay_payloads = v2_bundle_payloads(
        mode="replay",
        subject_mode=terminal["subject_mode"],
        artifact_root=artifact_root,
        artifact_dir=output_dir,
        source=source,
        inventory=inventory,
        authority=authority,
        review_plan=review,
        history=history,
        zero_sets=zero_sets,
        optional_payloads=optional_payloads,
        reproduction=reproduction,
        replay_equivalence=replay_equivalence,
    )
    publication = v2_publish_bundle(
        artifact_root=artifact_root,
        artifact_dir=output_dir,
        payloads=replay_payloads,
    )
    return {
        "schema": "frankensim.beads-template-hygiene.replay-result.v2",
        "terminal": "Pass",
        "subject_mode": terminal["subject_mode"],
        "artifact_dir": publication["artifact_dir"],
        "retained_terminal_root": terminal["semantic_root"],
        "replay_terminal_root": strict_json_loads(
            replay_payloads["terminal.json"],
            label="v2 replay terminal",
            require_canonical=True,
        )["semantic_root"],
        "live_tracker_reads": False,
        "live_tracker_writes": False,
        "network_access": False,
        "no_claim": (
            "offline replay proves exact retained reconstruction only and "
            "does not prove current tracker state or mint authority"
        ),
    }


def peek_replay_schema(artifact_root: str, input_dir: str) -> str:
    directory = resolve_run_dir(
        artifact_root,
        input_dir,
        label="replay input",
        must_exist=True,
    )
    payload = v2_read_file_strict(
        directory / "terminal.json",
        label="replay terminal",
    )
    document = strict_json_loads(payload, label="replay terminal")
    if not isinstance(document, dict) or not isinstance(document.get("schema"), str):
        raise InputRefused("replay terminal does not identify its schema")
    return document["schema"]


class V2CheckCollector:
    def __init__(self, case_id: str, assertion_id: str) -> None:
        self.case_id = case_id
        self.assertion_id = assertion_id
        self.rows: list[dict[str, Any]] = []

    @staticmethod
    def _bounded(value: Any) -> Any:
        payload = canonical_bytes(value)
        if len(payload) <= 2_048:
            return value
        return {
            "projection_root": semantic_root(value),
            "canonical_bytes": len(payload),
        }

    def check(
        self,
        check_id: str,
        condition: bool,
        *,
        expected: Any,
        observed: Any,
    ) -> None:
        row = v2_rooted(
            {
                "check_id": check_id,
                "expected": self._bounded(expected),
                "observed": self._bounded(observed),
                "passed": bool(condition),
            }
        )
        self.rows.append(row)
        if not condition:
            raise EvidenceFailed(
                f"{self.case_id}/{check_id} failed: "
                f"expected {expected!r}, observed {observed!r}"
            )

    def refuses(
        self,
        check_id: str,
        error_type: type[HarnessError],
        callback: Callable[[], Any],
        *,
        contains: str | None = None,
        projection: Callable[[], Any] | None = None,
        projection_label: str | None = None,
    ) -> None:
        if projection is None:
            projection = lambda: {
                "collector_check_roots": [
                    row["semantic_root"] for row in self.rows
                ],
                "command_receipts": list(_command_receipts),
                "signal_name": _signal_name,
            }
            projection_label = "HARNESS_PROCESS_STATE"
        elif not projection_label:
            raise EvidenceFailed(
                f"{self.case_id}/{check_id} stateful refusal lacks a projection label"
            )
        before_value = projection()
        before = semantic_root(before_value)
        error = expect_error(error_type, callback, contains=contains)
        after_value = projection()
        after = semantic_root(after_value)
        self.check(
            check_id,
            before == after and error.terminal == error_type.terminal,
            expected={
                "terminal": error_type.terminal,
                "unchanged_projection_root": before,
                "projection_label": projection_label,
            },
            observed={
                "terminal": error.terminal,
                "unchanged_projection_root": after,
                "projection_label": projection_label,
                "diagnostic_root": text_root(str(error)),
            },
        )

    def finish(self, parameters: Mapping[str, Any]) -> dict[str, Any]:
        if not self.rows:
            raise EvidenceFailed(f"{self.case_id} executed no assertion checks")
        document = {
            "schema": V2_ASSERTION_RESULT_SCHEMA,
            "case_id": self.case_id,
            "assertion_id": self.assertion_id,
            "executor_id": self.assertion_id,
            "parameters": dict(parameters),
            "parameter_root": semantic_root(parameters),
            "ordered_checks": self.rows,
            "ordered_check_roots": [row["semantic_root"] for row in self.rows],
            "check_count": len(self.rows),
            "terminal": "Pass",
            "no_claim": (
                "a passing harness assertion proves only the named bounded "
                "contract; it mints no target semantics or product authority"
            ),
        }
        return v2_rooted(document)


def v2_synthetic_target(
    index: int,
    *,
    issue_id: str | None = None,
    priority: int = 2,
    status: str = "open",
    disposition: str = "SUBSTANTIVE_SEMANTIC_OMISSION",
    issue_type: str = "task",
    missing_sections: Sequence[str] = ("## Acceptance Criteria",),
    readiness_seed: str = "REVIEW_ONLY",
    review_minutes: int = 1,
    retained_payload_bytes: int = 512,
    assignee: str = "",
) -> dict[str, Any]:
    issue_id = issue_id or f"synthetic-{index:04d}"
    field_roots = {
        field: text_root(f"{issue_id}:{field}")
        for field in ("description", "acceptance_criteria", "design", "notes")
    }
    dependency_root = semantic_root(
        {"dependencies": [], "dependents": []}
    )
    target_root = semantic_root(
        {
            "id": issue_id,
            "priority": priority,
            "status": status,
            "disposition": disposition,
            "type": issue_type,
            "missing_sections": sorted(missing_sections),
            "field_roots": field_roots,
            "dependency_neighborhood_root": dependency_root,
            "assignee": assignee,
        }
    )
    active_context = v2_active_work_context(
        {"status": status, "assignee": assignee}
    )
    return v2_rooted(
        {
            "id": issue_id,
            "issue_id": issue_id,
            "title": f"Synthetic target {index:04d}",
            "type": issue_type,
            "issue_type": issue_type,
            "status": status,
            "priority": priority,
            "priority_lane": f"P{priority}",
            "status_lane": status,
            "lane": f"P{priority}/{status}",
            "campaign_epoch_root": "",
            "destination": (
                "history-v2.json" if status == "closed" else "review-plan-v2.json"
            ),
            "movement_destination": (
                "history-v2.json" if status == "closed" else "review-plan-v2.json"
            ),
            "tracker_assignee": assignee,
            "coordination_assignee": assignee,
            "tracker_owner": "",
            "parent": "",
            "labels": [],
            "missing_sections": sorted(missing_sections),
            "disposition": disposition,
            "disposition_rationale": "synthetic exact disposition",
            "disposition_falsifier": "a contrary rooted fixture invalidates it",
            "field_roots": field_roots,
            "field_byte_lengths": {
                "description": 32,
                "acceptance_criteria": 32,
                "design": 32,
                "notes": 32,
            },
            "clause_roots": [
                {
                    "field": "description",
                    "line_index": 0,
                    "byte_length": 32,
                    "text_root": text_root(f"{issue_id}:clause"),
                    "truncated": False,
                }
            ],
            "target_root": target_root,
            "v1_row_root": semantic_root({"id": issue_id, "v1": True}),
            "dependency_neighborhood_root": dependency_root,
            "dependencies": [],
            "dependents": [],
            "domain_candidates": [],
            "domain_candidate": [],
            "declared_domain_owner": "",
            "declared_acceptance_owner": "",
            "implementation_owner": "UNRESOLVED",
            "evidence_owner": "UNRESOLVED",
            "terminal_consumer": "UNRESOLVED",
            "reviewer_provenance": {
                "state": "NODATA",
                "receipt_root": None,
            },
            "source_closure": {"required": False, "state": "NOT_REQUIRED"},
            "user_effect": "make exact review debt visible to users",
            "target_implementation_estimated_minutes": 60,
            "target_implementation_estimate_minutes": 60,
            "target_implementation_estimate_state": "DECLARED",
            "review_minutes": review_minutes,
            "generated_child_estimated_minutes": review_minutes,
            "retained_payload_bytes": retained_payload_bytes,
            "external_authority_adapter_identity": "",
            "external_authority_receipt_root": "",
            "external_authority_verdict": "NODATA",
            "conditional_write_capability_identity": "",
            "conditional_write_receipt_root": "",
            "conditional_write_verdict": "NODATA",
            "readiness": readiness_seed,
            "remediation_route": "ANALYSIS_ONLY",
            "active_work_context": active_context,
            "no_claim": "synthetic planning row grants no authority",
        }
    )


def v2_synthetic_inventory(
    targets: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    campaign_root = semantic_root(
        [
            {"id": row["id"], "target_root": row["target_root"]}
            for row in sorted(targets, key=lambda value: value["id"])
        ]
    )
    rows: list[dict[str, Any]] = []
    for row in sorted(targets, key=lambda value: value["id"]):
        normalized = dict(row)
        normalized.pop("semantic_root", None)
        normalized["campaign_epoch_root"] = campaign_root
        rows.append(v2_rooted(normalized))
    return v2_rooted(
        {
            "schema": V2_INVENTORY_SCHEMA,
            "v1_inventory_root": semantic_root({"synthetic": True}),
            "campaign_epoch_root": campaign_root,
            "filters": {"priorities": [], "statuses": [], "empty_means_all": True},
            "rows": rows,
            "counts": {
                "targets": len(rows),
                "nonclosed": sum(row["status"] != "closed" for row in rows),
                "history": sum(row["status"] == "closed" for row in rows),
                "warnings": sum(len(row["missing_sections"]) for row in rows),
            },
            "nonclosed_ids_root": semantic_root(
                sorted(row["id"] for row in rows if row["status"] != "closed")
            ),
            "history_ids_root": semantic_root(
                sorted(row["id"] for row in rows if row["status"] == "closed")
            ),
            "no_claim": "synthetic inventory is test evidence only",
        }
    )


def v2_synthetic_receipts(
    inventory: Mapping[str, Any],
    *,
    compatible: bool = False,
    declared: bool = False,
    manual: bool = False,
    external_verdict: str = "NODATA",
    conditional_verdict: str = "NODATA",
    gate_verdict: str = "NODATA",
) -> dict[str, Any]:
    rows = list(inventory["rows"])
    target_ids = sorted(row["id"] for row in rows)
    compatibility_key = "synthetic-compatible" if compatible else ""
    compatibility_root = ""
    if compatible:
        compatibility_root = semantic_root(
            {
                "compatibility_key": compatibility_key,
                "coordination": sorted((row["id"], "") for row in rows),
                "declared_domain": sorted(
                    (row["id"], "synthetic-domain") for row in rows
                ),
                "user_effect": sorted(
                    (row["id"], row["user_effect"]) for row in rows
                ),
                "dependency_neighborhood": sorted(
                    (row["id"], row["dependency_neighborhood_root"]) for row in rows
                ),
                "target_ids": target_ids,
                "target_roots": sorted(
                    (row["id"], row["target_root"]) for row in rows
                ),
                "rationale": "all bound synthetic dimensions agree",
                "falsifier": "any bound dimension differs",
            }
        )
    receipts: list[dict[str, Any]] = []
    for target in rows:
        receipt = {
            "target_id": target["id"],
            "target_root": target["target_root"],
            "inventory_root": inventory["semantic_root"],
            "campaign_epoch_root": inventory["campaign_epoch_root"],
            "reviewer": (
                "synthetic-reviewer" if declared or compatible else ""
            ),
            "reviewer_kind": "fixture" if declared or compatible else "",
            "coordination_assignee": "",
            "declared_domain_owner": (
                "synthetic-domain" if compatible or declared else ""
            ),
            "declared_acceptance_owner": (
                "synthetic-acceptance-owner"
                if declared or compatible
                else ""
            ),
            "implementation_owner": "synthetic-implementer",
            "evidence_owner": "synthetic-evidence-owner",
            "terminal_consumer": "synthetic-consumer",
            "source_closure": {"required": False, "state": "NOT_REQUIRED"},
            "user_effect": target["user_effect"],
            "review_minutes": target["review_minutes"],
            "compatibility_key": compatibility_key,
            "compatibility_target_ids": target_ids if compatible else [],
            "compatibility_receipt_root": compatibility_root,
            "compatibility_rationale": (
                "all bound synthetic dimensions agree" if compatible else ""
            ),
            "compatibility_falsifier": (
                "any bound dimension differs" if compatible else ""
            ),
            "manual_authorization": manual,
            "manual_authorization_source": (
                "synthetic-owner-declaration" if manual else ""
            ),
            "external_authority_adapter": (
                "fixture-independent-adapter"
                if external_verdict == "VALID"
                else ""
            ),
            "external_authority_receipt_root": (
                semantic_root(
                    {"gate": "external", "target_root": target["target_root"]}
                )
                if external_verdict == "VALID"
                else ""
            ),
            "external_authority_verdict": external_verdict,
            "conditional_capability_identity": (
                "fixture-atomic-cas"
                if conditional_verdict == "VALID"
                else ""
            ),
            "conditional_capability_receipt_root": (
                semantic_root(
                    {"gate": "conditional", "target_root": target["target_root"]}
                )
                if conditional_verdict == "VALID"
                else ""
            ),
            "conditional_capability_verdict": conditional_verdict,
            "gate_admission_receipt_root": (
                semantic_root(
                    {"gate": "admission", "target_root": target["target_root"]}
                )
                if gate_verdict == "VALID"
                else ""
            ),
            "gate_admission_verdict": gate_verdict,
            "no_claim": "synthetic receipt grants no live authority",
        }
        receipts.append(v2_rooted(receipt))
    return v2_rooted(
        {
            "schema": V2_REVIEW_RECEIPTS_SCHEMA,
            "inventory_root": inventory["semantic_root"],
            "campaign_epoch_root": inventory["campaign_epoch_root"],
            "receipts": receipts,
            "source": "SYNTHETIC_FIXTURE",
            "no_claim": "synthetic review receipts are test inputs only",
        }
    )


def v2_synthetic_authority(
    inventory: Mapping[str, Any],
    *,
    compatible: bool = False,
    declared: bool = False,
    manual: bool = False,
    external_verdict: str = "NODATA",
    conditional_verdict: str = "NODATA",
    gate_verdict: str = "NODATA",
    version: str = "0.2.19",
    allow_mechanical_fixture: bool = False,
) -> dict[str, Any]:
    receipts = v2_synthetic_receipts(
        inventory,
        compatible=compatible,
        declared=declared,
        manual=manual,
        external_verdict=external_verdict,
        conditional_verdict=conditional_verdict,
        gate_verdict=gate_verdict,
    )
    return v2_derive_authority(
        inventory,
        receipts,
        current_br_version=version,
        allow_mechanical_fixture=allow_mechanical_fixture,
    )


def v2_empty_prior_campaign() -> dict[str, Any]:
    return v2_rooted(
        {
            "state": "NOT_PROVIDED",
            "bundle_path": None,
            "terminal_root": None,
            "inventory_root": None,
            "rows": [],
            "no_claim": "synthetic prior campaign is not provided",
        }
    )


def v2_minimal_source(
    *,
    manifest_root: str,
    subject: str,
) -> dict[str, Any]:
    observation = v2_rooted(
        {
            "br_version": {"version": "0.2.19"},
            "command_receipts": [],
            "subject": subject,
        }
    )
    audit = v2_rooted(
        {
            "state": "NOT_REQUESTED",
            "documents": [],
            "command_receipts": [],
        }
    )
    return v2_rooted(
        {
            "schema": V2_SOURCE_SCHEMA,
            "manifest": {
                "semantic_root": manifest_root,
                "content_identity": manifest_root,
            },
            "captured": {
                "observation": observation,
                "audit_capture": audit,
            },
            "no_claim": "minimal in-memory artifact source",
        }
    )


def v2_synthetic_bundle_components(
    manifest: Mapping[str, Any],
    *,
    target_count: int = 2,
    oversize: bool = False,
) -> tuple[
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, bytes],
]:
    targets = [
        v2_synthetic_target(
            index,
            retained_payload_bytes=(
                V2_CHILD_PAYLOAD_CAP + 1 if oversize and index == 0 else 512
            ),
        )
        for index in range(target_count)
    ]
    inventory = v2_synthetic_inventory(targets)
    authority = v2_synthetic_authority(
        inventory,
        compatible=target_count > 1,
    )
    review, optional = v2_build_review_plan(
        inventory,
        authority,
        max_targets=V2_REVIEW_TARGET_DEFAULT,
    )
    source = v2_minimal_source(
        manifest_root=manifest["semantic_root"],
        subject="synthetic-bundle",
    )
    history = v2_not_requested(
        schema=V2_HISTORY_SCHEMA,
        mode="review-plan",
        source_root=source["semantic_root"],
        inventory_root=inventory["semantic_root"],
        projection="history",
    )
    zero = v2_build_zero_sets(
        source_root=source["semantic_root"],
        inventory=inventory,
        prior_campaign=v2_empty_prior_campaign(),
    )
    return source, inventory, authority, review, history, zero, optional


def v2_fixture_manifest_with_history(
    manifest: Mapping[str, Any],
    history_contract: Mapping[str, Any],
) -> dict[str, Any]:
    fixture = dict(manifest)
    fixture["history_contract"] = dict(history_contract)
    fixture["content_identity"] = semantic_root(
        {
            "base_manifest_root": manifest["semantic_root"],
            "history_contract": history_contract,
            "purpose": "artifact-only self-test trust anchor",
        }
    )
    fixture.pop("semantic_root", None)
    fixture["semantic_root"] = semantic_root(fixture)
    return fixture


def v2_replayable_fixture_bundle(
    manifest: Mapping[str, Any],
    *,
    subject_mode: str,
    issue: Mapping[str, Any] | None = None,
    oversize: bool = False,
) -> tuple[
    dict[str, bytes],
    dict[str, Any],
    tuple[
        dict[str, Any],
        dict[str, Any],
        dict[str, Any],
        dict[str, Any],
        dict[str, Any],
        dict[str, Any],
        dict[str, bytes],
    ],
]:
    if subject_mode not in {"review-plan", "history-plan"}:
        raise EvidenceFailed("replayable fixture subject mode is unknown")
    history_mode = subject_mode == "history-plan"
    raw_issue = dict(
        issue
        or fixture_issue(
            issue_id=(
                "fixture-history-replay"
                if history_mode
                else "fixture-review-replay"
            ),
            status="closed" if history_mode else "open",
            priority=2,
            description=(
                "x" * (V2_CHILD_DESCRIPTION_CAP + 1)
                if oversize
                else "Exact released-tool fixture work with replay boundaries."
            ),
            acceptance="",
            labels=("reality-check", "authority:template-hygiene"),
        )
    )
    if history_mode:
        raw_issue.update(
            {
                "created_at": "2026-01-01T00:00:00+00:00",
                "created_by": "fixture-creator",
                "closed_at": "2026-01-02T00:00:00+00:00",
                "close_reason": "Fixture history replay completed.",
            }
        )
    full_issue = v2_full_issue_projection(raw_issue)
    missing = {
        full_issue["id"]: fixture_expected_missing(full_issue["type"])
    }
    lint = fixture_lint([full_issue], missing)
    observation = v2_rooted(
        {
            "schema": V2_SOURCE_SCHEMA,
            "capture_contract": {
                "count": 2,
                "coherent": True,
                "tracker_cli": "br",
                "direct_tracker_file_access": False,
                "network_access": False,
            },
            "br_version": {
                "version": "0.2.19",
                "build": "release",
                "commit": "fixture-released-br",
                "target": "fixture",
                "features": [],
            },
            "br_capabilities": {
                "contract_version": "br.capabilities.v1",
                "commands": {},
                "operation_count": 0,
            },
            "tracker_status": {"fixture": True},
            "sync_status": {"fixture": True},
            "export_witness": {"fixture": True},
            "case_manifest_root": manifest["semantic_root"],
            "case_manifest_content_identity": manifest["content_identity"],
            "source_files": [],
            "live_issue_count": 1,
            "live_issue_projection_root": semantic_root([full_issue]),
            "lint_issue_ids_root": semantic_root([full_issue["id"]]),
            "status_cut": list(STATUS_SCOPES),
            "command_receipts": [],
            "no_claim": (
                "released-tool fixture observation is tracker-read-only and "
                "contains no raw command streams"
            ),
        }
    )
    v1_inventory = assemble_inventory(
        lint,
        [full_issue],
        observation,
    )
    snapshot = LiveSnapshot(
        lint,
        [full_issue],
        observation,
        v1_inventory,
        {},
        (full_issue,),
    )
    campaign_root, _ = v2_campaign_epoch(snapshot)
    inventory = v2_build_inventory(
        snapshot,
        campaign_epoch_root=campaign_root,
        priority_filter=None,
        status_filter=None,
    )
    receipts = v2_empty_review_receipts(
        inventory_root=inventory["semantic_root"],
        campaign_epoch_root=inventory["campaign_epoch_root"],
    )
    authority = v2_derive_authority(
        inventory,
        receipts,
        current_br_version="0.2.19",
    )
    prior = v2_empty_prior_campaign()
    accepted_manifest = dict(manifest)
    if history_mode:
        normalized_audit = v2_normalize_audit_document(
            {
                "issue_id": full_issue["id"],
                "events": [
                    {
                        "id": 2,
                        "event_type": "closed",
                        "actor": "fixture-closer",
                        "timestamp": "2026-01-02T00:00:00.001000+00:00",
                        "comment": full_issue["close_reason"],
                    },
                    {
                        "id": 1,
                        "event_type": "status_changed",
                        "actor": "fixture-closer",
                        "timestamp": "2026-01-02T00:00:00+00:00",
                        "old_value": "open",
                        "new_value": "closed",
                    },
                ],
            },
            issue_id=full_issue["id"],
        )
        audit = v2_rooted(
            {
                "capture_count": 2,
                "worker_bound": 1,
                "issue_ids": [full_issue["id"]],
                "documents": [normalized_audit],
                "command_receipts": [],
                "raw_stream_bodies_retained": False,
                "no_claim": "artifact-only released-br audit fixture",
            }
        )
        history_contract = dict(manifest["history_contract"])
        history_contract.update(
            {
                "legacy_coverage_anchor_issue": full_issue["id"],
                "legacy_coverage_anchor_closed_at": full_issue["closed_at"],
                "legacy_coverage_anchor_status_event_id": 1,
                "legacy_coverage_anchor_close_event_id": 2,
                "legacy_coverage_count": 0,
                "legacy_coverage_rows_root": semantic_root([]),
                "known_pair_max_skew_ms": 1_000,
            }
        )
        accepted_manifest = v2_fixture_manifest_with_history(
            manifest,
            history_contract,
        )
    else:
        audit = v2_rooted(
            {
                "capture_count": 0,
                "worker_bound": 0,
                "issue_ids": [],
                "documents": [],
                "command_receipts": [],
                "raw_stream_bodies_retained": False,
                "no_claim": "history was not requested",
            }
        )
    source = v2_build_source_document(
        snapshot=snapshot,
        manifest=accepted_manifest,
        inventory=inventory,
        review_receipts=receipts,
        prior_campaign=prior,
        audit_capture=audit,
        priority_filter=None,
        status_filter=None,
    )
    if history_mode:
        review = v2_not_requested(
            schema=V2_REVIEW_PLAN_SCHEMA,
            mode=subject_mode,
            source_root=source["semantic_root"],
            inventory_root=inventory["semantic_root"],
            projection="review-plan",
        )
        optional: dict[str, bytes] = {}
        history = v2_build_history(
            source_root=source["semantic_root"],
            inventory=inventory,
            authority=authority,
            all_issues=snapshot.all_issues,
            audit_capture=audit,
            history_contract=accepted_manifest["history_contract"],
        )
    else:
        review, optional = v2_build_review_plan(
            inventory,
            authority,
            max_targets=V2_REVIEW_TARGET_DEFAULT,
        )
        history = v2_not_requested(
            schema=V2_HISTORY_SCHEMA,
            mode=subject_mode,
            source_root=source["semantic_root"],
            inventory_root=inventory["semantic_root"],
            projection="history",
        )
    zero = v2_build_zero_sets(
        source_root=source["semantic_root"],
        inventory=inventory,
        prior_campaign=prior,
    )
    reproduction = [
        str(SCRIPT_REL),
        f"--{subject_mode}",
        "--artifact-root",
        "target/v2-fixture",
        "--artifact-dir",
        f"{subject_mode}-run",
    ]
    payloads = v2_bundle_payloads(
        mode=subject_mode,
        subject_mode=subject_mode,
        artifact_root="target/v2-fixture",
        artifact_dir=f"{subject_mode}-run",
        source=source,
        inventory=inventory,
        authority=authority,
        review_plan=review,
        history=history,
        zero_sets=zero,
        optional_payloads=optional,
        reproduction=reproduction,
    )
    components = (
        source,
        inventory,
        authority,
        review,
        history,
        zero,
        optional,
    )
    return payloads, accepted_manifest, components


def v2_plan_for_count(
    count: int,
    *,
    max_targets: int,
    review_minutes: int = 1,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, bytes]]:
    inventory = v2_synthetic_inventory(
        [
            v2_synthetic_target(index, review_minutes=review_minutes)
            for index in range(count)
        ]
    )
    authority = v2_synthetic_authority(
        inventory,
        compatible=count > 1,
    )
    plan, optional = v2_build_review_plan(
        inventory,
        authority,
        max_targets=max_targets,
    )
    return plan, authority, optional


def v2_execute_schema_cli_ux(
    case: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    slug = str(case["id"]).removeprefix("template-lint-v2.")
    checks = V2CheckCollector(case["id"], case["assertion_id"])
    parameters = {"family": "schema-cli-ux", "slug": slug}
    if slug == "schema-closed-manifest":
        checks.check(
            "top-level-closed",
            set(manifest) - {"content_identity", "semantic_root"}
            == V2_MANIFEST_TOP_LEVEL_KEYS,
            expected=sorted(V2_MANIFEST_TOP_LEVEL_KEYS),
            observed=sorted(set(manifest) - {"content_identity", "semantic_root"}),
        )
        duplicate_rejected = False
        try:
            tomllib.loads("value = 1\nvalue = 2\n")
        except tomllib.TOMLDecodeError:
            duplicate_rejected = True
        checks.check(
            "duplicate-toml-key-refused",
            duplicate_rejected,
            expected=True,
            observed=duplicate_rejected,
        )
        checks.check(
            "schema-version",
            manifest["schema"] == V2_MANIFEST_SCHEMA
            and manifest["schema_version"] == 2,
            expected=(V2_MANIFEST_SCHEMA, 2),
            observed=(manifest["schema"], manifest["schema_version"]),
        )
    elif slug == "schema-executable-assertions":
        case_ids = [row["id"] for row in manifest["case"]]
        assertion_ids = [row["assertion_id"] for row in manifest["case"]]
        criterion_links = {
            assertion_id
            for row in manifest["criterion"]
            for assertion_id in row["assertion_ids"]
        }
        checks.check(
            "case-assertion-bijection",
            len(case_ids) == len(set(case_ids)) == len(assertion_ids)
            == len(set(assertion_ids))
            == 96,
            expected=96,
            observed={
                "cases": len(set(case_ids)),
                "assertions": len(set(assertion_ids)),
            },
        )
        checks.check(
            "criterion-links-executable",
            criterion_links == set(assertion_ids),
            expected=semantic_root(sorted(assertion_ids)),
            observed=semantic_root(sorted(criterion_links)),
        )
        registry = v2_build_executor_registry(manifest)
        checks.check(
            "compiled-executor-bijection",
            set(registry) == set(assertion_ids)
            and all(callable(executor) for executor in registry.values()),
            expected=semantic_root(sorted(assertion_ids)),
            observed=semantic_root(sorted(registry)),
        )
    elif slug == "cli-help-modes":
        completed = subprocess.run(
            [str(REPO_ROOT / SCRIPT_REL), "--help"],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
            timeout=CAPS["subprocess_timeout_seconds"],
        )
        expected_options = {
            "--review-plan",
            "--history-plan",
            "--self-test-v2",
            "--case-v2",
            "--review-receipts",
            "--prior-campaign",
            "--priorities",
            "--statuses",
            "--max-targets-per-child",
            "--explain-target",
            "--artifact-root",
            "--artifact-dir",
        }
        checks.check(
            "help-exit-zero",
            completed.returncode == 0,
            expected=0,
            observed=completed.returncode,
        )
        checks.check(
            "help-complete-options",
            all(option in completed.stdout for option in expected_options),
            expected=sorted(expected_options),
            observed={
                "stdout_root": text_root(completed.stdout),
                "stdout_bytes": len(completed.stdout.encode("utf-8")),
            },
        )
    elif slug == "cli-exclusive-mode":
        checks.refuses(
            "mode-required",
            UsageRefused,
            lambda: parse_arguments([]),
        )
        checks.refuses(
            "multiple-modes",
            UsageRefused,
            lambda: parse_arguments(["--review-plan", "--history-plan"]),
        )
        checks.refuses(
            "duplicate-mode",
            UsageRefused,
            lambda: parse_arguments(["--review-plan", "--review-plan"]),
            contains="duplicate",
        )
        parsed = parse_arguments(["--self-test-v2"])
        checks.check(
            "one-mode-admitted",
            parsed.self_test_v2 is True,
            expected=True,
            observed=parsed.self_test_v2,
        )
    elif slug == "cli-artifact-grammar":
        missing_root = parse_arguments(
            ["--review-plan", "--artifact-dir", "fixture-run"]
        )
        checks.refuses(
            "artifact-root-required",
            UsageRefused,
            lambda: require_v2_artifact_grammar(missing_root),
            contains="--artifact-root",
        )
        missing_dir = parse_arguments(
            ["--review-plan", "--artifact-root", "target/v2"]
        )
        checks.refuses(
            "artifact-dir-required",
            UsageRefused,
            lambda: require_v2_artifact_grammar(missing_dir),
            contains="--artifact-dir",
        )
        checks.refuses(
            "output-unsupported",
            UsageRefused,
            lambda: parse_arguments(["--review-plan", "--output", "x"]),
            contains="unsupported",
        )
        valid = parse_arguments(
            [
                "--review-plan",
                "--artifact-root",
                "target/v2",
                "--artifact-dir",
                "fixture-run",
            ]
        )
        checks.check(
            "exact-artifact-grammar-admitted",
            require_v2_artifact_grammar(valid) == "fixture-run",
            expected="fixture-run",
            observed=valid.artifact_dir,
        )
    elif slug == "cli-path-admission":
        valid = safe_relative("target/v2/évidence", label="fixture")
        checks.check(
            "nfc-relative-admitted",
            str(valid) == "target/v2/évidence",
            expected="target/v2/évidence",
            observed=str(valid),
        )
        hostile = (
            "",
            ".",
            "..",
            "a/../b",
            "/absolute",
            "~/alias",
            "C:\\alias",
            "a\\b",
            "bad\u0001path",
            "target/v2/e\u0301vidence",
        )
        for index, value in enumerate(hostile):
            checks.refuses(
                f"unsafe-path-{index:02d}",
                UsageRefused,
                lambda value=value: safe_relative(value, label="fixture"),
            )
    elif slug == "cli-review-receipts":
        inventory = v2_synthetic_inventory([v2_synthetic_target(0)])
        empty = v2_load_review_receipts(None, inventory=inventory)
        authority = v2_derive_authority(
            inventory,
            empty,
            current_br_version="0.2.19",
        )
        checks.check(
            "absent-receipt-review-only",
            authority["decisions"][0]["readiness"] == "REVIEW_ONLY"
            and authority["decisions"][0]["remediation_route"] == "ANALYSIS_ONLY",
            expected=("REVIEW_ONLY", "ANALYSIS_ONLY"),
            observed=(
                authority["decisions"][0]["readiness"],
                authority["decisions"][0]["remediation_route"],
            ),
        )
        checks.refuses(
            "duplicate-json-key",
            InputRefused,
            lambda: strict_json_loads(
                '{"target":1,"target":2}\n',
                label="fixture receipt",
            ),
            contains="duplicate",
        )
    elif slug == "cli-filters-explain":
        first = v2_parse_csv_set(
            "P2,P0",
            allowed={f"P{priority}" for priority in range(5)},
            label="priorities",
        )
        second = v2_parse_csv_set(
            "P0,P2",
            allowed={f"P{priority}" for priority in range(5)},
            label="priorities",
        )
        checks.check(
            "filter-order-invariant",
            first == second == {"P0", "P2"},
            expected=["P0", "P2"],
            observed=sorted(first or []),
        )
        checks.refuses(
            "unknown-filter-refused",
            UsageRefused,
            lambda: v2_parse_csv_set(
                "P9",
                allowed={"P0"},
                label="priorities",
            ),
            contains="unknown",
        )
    elif slug in {"ux-empty-synopsis", "ux-large-synopsis"}:
        count = 0 if slug == "ux-empty-synopsis" else 13
        inventory = v2_synthetic_inventory(
            [v2_synthetic_target(index) for index in range(count)]
        )
        authority = v2_synthetic_authority(inventory)
        review = v2_rooted(
            {
                "schema": V2_REVIEW_PLAN_SCHEMA,
                "children": [
                    {"priority": 2, "target_ids": [row["id"]]}
                    for row in inventory["rows"]
                ],
                "oversize_content": [],
                "work": {
                    "total_review_minutes": count,
                    "minimum_child_review_minutes": 0 if not count else 1,
                    "maximum_child_review_minutes": 0 if not count else 1,
                },
            }
        )
        history = v2_rooted(
            {"schema": V2_HISTORY_SCHEMA, "rows": []}
        )
        zero = v2_rooted(
            {
                "schema": V2_ZERO_SETS_SCHEMA,
                "counts": {"zero_receipts": 25 if not count else 24},
            }
        )
        synopsis = v2_synopsis(
            subject_mode="review-plan",
            review_plan=review,
            history=history,
            authority=authority,
            zero_sets=zero,
            artifact_dir="target/v2/fixture",
            reproduction=["script", "--review-plan"],
            explain_target=None,
        )
        checks.check(
            "synopsis-bounded-utf8",
            len(canonical_bytes(synopsis)) <= V2_SYNOPSIS_BYTES_CAP,
            expected=f"<= {V2_SYNOPSIS_BYTES_CAP}",
            observed=len(canonical_bytes(synopsis)),
        )
        checks.check(
            "synopsis-id-preview",
            synopsis["selected_ids"]["shown"] == min(count, 12)
            and synopsis["selected_ids"]["total"] == count,
            expected={"shown": min(count, 12), "total": count},
            observed=synopsis["selected_ids"],
        )
        if count:
            checks.check(
                "synopsis-truncation-explicit",
                synopsis["selected_ids"]["truncated"] is True
                and synopsis["selected_ids"]["notice"] is not None,
                expected=True,
                observed=synopsis["selected_ids"],
            )
    elif slug == "ux-streams-color-tty":
        document = {
            "terminal": "Pass",
            "no_color": True,
            "tty_independent": True,
        }
        encoded = canonical_bytes(document)
        checks.check(
            "no-ansi",
            b"\x1b[" not in encoded,
            expected="no ANSI escape",
            observed="absent" if b"\x1b[" not in encoded else "present",
        )
        checks.check(
            "semantic-stream-determinism",
            encoded == canonical_bytes(dict(reversed(list(document.items())))),
            expected=text_root(encoded.decode("utf-8")),
            observed=text_root(
                canonical_bytes(dict(reversed(list(document.items())))).decode("utf-8")
            ),
        )
    elif slug == "compat-v1-immutable":
        v1 = load_case_manifest()
        compatibility = manifest["compatibility_contract"]
        payload = bounded_read(REPO_ROOT / CASE_MANIFEST_REL)
        checks.check(
            "v1-content-root",
            "sha256-v1:" + hashlib.sha256(payload).hexdigest()
            == compatibility["v1_manifest_content_root"],
            expected=compatibility["v1_manifest_content_root"],
            observed="sha256-v1:" + hashlib.sha256(payload).hexdigest(),
        )
        checks.check(
            "v1-semantic-root",
            v1["semantic_root"]
            == compatibility["v1_case_manifest_semantic_root"],
            expected=compatibility["v1_case_manifest_semantic_root"],
            observed=v1["semantic_root"],
        )
        checks.check(
            "v1-registry-counts",
            len(v1["case"]) == 21
            and compatibility["v1_assertion_count"] == 148,
            expected={"cases": 21, "assertions": 148},
            observed={
                "cases": len(v1["case"]),
                "assertions": compatibility["v1_assertion_count"],
            },
        )
    else:
        raise EvidenceFailed(f"no schema/CLI/UX executor for {slug}")
    return checks.finish(parameters)


def v2_execute_source_cases(
    case: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    del manifest
    slug = str(case["id"]).removeprefix("template-lint-v2.")
    checks = V2CheckCollector(case["id"], case["assertion_id"])
    parameters = {"family": "source", "slug": slug}

    def capture_state(
        rows: Sequence[Mapping[str, Any]],
        *,
        lint: Mapping[str, Any] | None = None,
        tracker_status: Mapping[str, Any] | None = None,
        capabilities: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        ordered = sorted((dict(row) for row in rows), key=lambda row: row["id"])
        return v2_capture_state(
            issue_ids=[str(row["id"]) for row in ordered],
            full_issues=ordered,
            lint=dict(lint or {"warnings": []}),
            tracker_status=dict(tracker_status or {"state": "clean"}),
            sync_status={"state": "synchronized"},
            export_witness={"state": "current"},
            version={"version": "0.2.19"},
            capabilities=dict(
                capabilities or {"contract_version": "br.capabilities.v1"}
            ),
        )

    if slug == "source-zero-one":
        for count in (0, 1):
            inventory = v2_synthetic_inventory(
                [v2_synthetic_target(0)] if count else []
            )
            checks.check(
                f"cardinality-{count}",
                inventory["counts"]["targets"] == count
                and len(inventory["rows"]) == count,
                expected=count,
                observed=inventory["counts"]["targets"],
            )
    elif slug in {
        "source-rows-9-10-11",
        "source-rows-12-13",
        "source-rows-24-25-26",
    }:
        configurations = {
            "source-rows-9-10-11": (10, (9, 10, 11)),
            "source-rows-12-13": (13, (12, 13)),
            "source-rows-24-25-26": (25, (24, 25, 26)),
        }
        cap, populations = configurations[slug]
        for population in populations:
            plan, _, _ = v2_plan_for_count(
                population,
                max_targets=cap,
            )
            expected_children = math.ceil(population / cap)
            mapped = [
                target_id
                for child in plan["children"]
                for target_id in child["target_ids"]
            ]
            checks.check(
                f"population-{population}",
                len(plan["children"]) == expected_children
                and len(mapped) == population
                and len(mapped) == len(set(mapped)),
                expected={"children": expected_children, "targets": population},
                observed={
                    "children": len(plan["children"]),
                    "targets": len(mapped),
                },
            )
    elif slug == "source-rows-4096":
        v2_validate_source_limits(
            row_count=4_096,
            warning_count=0,
            maximum_warnings_per_issue=0,
        )
        checks.check(
            "row-cap-admitted",
            True,
            expected=4_096,
            observed=4_096,
        )
        checks.refuses(
            "row-cap-plus-one",
            InputRefused,
            lambda: v2_validate_source_limits(
                row_count=4_097,
                warning_count=0,
                maximum_warnings_per_issue=0,
            ),
            contains="row count",
        )
    elif slug == "source-warnings-8192":
        v2_validate_source_limits(
            row_count=2_731,
            warning_count=8_192,
            maximum_warnings_per_issue=3,
        )
        checks.check(
            "warning-cap-admitted",
            True,
            expected=8_192,
            observed=8_192,
        )
        checks.refuses(
            "warning-cap-plus-one",
            InputRefused,
            lambda: v2_validate_source_limits(
                row_count=2_731,
                warning_count=8_193,
                maximum_warnings_per_issue=3,
            ),
            contains="warning count",
        )
    elif slug == "source-all-types-statuses-priorities":
        statuses = ["open", "in_progress", "blocked", "deferred", "closed"]
        types = ["bug", "task", "feature", "epic", "chore", "docs", "question", "custom:x"]
        rows = [
            v2_synthetic_target(
                index,
                priority=index % 5,
                status=statuses[index % len(statuses)],
                issue_type=types[index % len(types)],
            )
            for index in range(40)
        ]
        checks.check(
            "covering-array-priorities",
            {row["priority"] for row in rows} == set(range(5)),
            expected=list(range(5)),
            observed=sorted({row["priority"] for row in rows}),
        )
        checks.check(
            "covering-array-statuses-types",
            {row["status"] for row in rows} == set(statuses)
            and {row["type"] for row in rows} == set(types),
            expected={"statuses": statuses, "types": types},
            observed={
                "statuses": sorted({row["status"] for row in rows}),
                "types": sorted({row["type"] for row in rows}),
            },
        )
    elif slug == "source-warning-partitions":
        combinations = [
            ("A",),
            ("S",),
            ("C",),
            ("A", "S"),
            ("A", "C"),
            ("S", "C"),
            ("A", "S", "C"),
        ]
        union = [f"row-{index}" for index in range(len(combinations))]
        partitions = {
            "+".join(values): [union[index]]
            for index, values in enumerate(combinations)
        }
        checks.check(
            "seven-exclusive-warning-partitions",
            len(partitions) == 7
            and len({value[0] for value in partitions.values()}) == 7,
            expected=7,
            observed=len(partitions),
        )
        zero_receipt = v2_rooted(
            {"cell": "P4/closed", "ids": [], "count": 0}
        )
        checks.check(
            "zero-cell-rooted",
            zero_receipt["semantic_root"] == semantic_root(
                {key: value for key, value in zero_receipt.items() if key != "semantic_root"}
            ),
            expected=True,
            observed=zero_receipt["semantic_root"],
        )
    elif slug == "source-double-capture-coherent":
        rows = [
            v2_synthetic_target(0, priority=0, status="open"),
            v2_synthetic_target(1, priority=4, status="closed"),
        ]
        first = capture_state(rows)
        second = capture_state(list(reversed(rows)))
        coherent_root = v2_validate_capture_pair(first, second)
        checks.check(
            "production-coherence-validator-admits-equal-observations",
            coherent_root == first["semantic_root"] == second["semantic_root"],
            expected=first["semantic_root"],
            observed=coherent_root,
        )
        checks.check(
            "capture-binds-all-source-dimensions",
            set(first) == V2_CAPTURE_STATE_KEYS,
            expected=sorted(V2_CAPTURE_STATE_KEYS),
            observed=sorted(first),
        )
    elif slug == "source-within-run-drift":
        baseline_rows = [
            v2_synthetic_target(0, priority=1, status="open"),
            v2_synthetic_target(1, priority=2, status="blocked"),
        ]
        baseline = capture_state(baseline_rows)
        mutations: list[tuple[str, dict[str, Any]]] = []
        for dimension, field, value in (
            ("field", "title", "changed title"),
            ("relation", "dependencies", [{"id": "other", "type": "blocks"}]),
            ("status", "status", "closed"),
            ("priority", "priority", 3),
            ("owner", "tracker_owner", "changed-owner"),
            ("consumer", "terminal_consumer", "changed-consumer"),
        ):
            changed_rows = json.loads(json.dumps(baseline_rows))
            changed_rows[0][field] = value
            mutations.append((dimension, capture_state(changed_rows)))
        mutations.extend(
            [
                (
                    "tracker-status",
                    capture_state(
                        baseline_rows,
                        tracker_status={"state": "dirty"},
                    ),
                ),
                (
                    "capability",
                    capture_state(
                        baseline_rows,
                        capabilities={
                            "contract_version": "br.capabilities.v2"
                        },
                    ),
                ),
                (
                    "lint",
                    capture_state(
                        baseline_rows,
                        lint={"warnings": [{"id": baseline_rows[0]["id"]}]},
                    ),
                ),
            ]
        )
        for dimension, changed in mutations:
            checks.refuses(
                f"drift-{dimension}",
                InputRefused,
                lambda changed=changed: v2_validate_capture_pair(
                    baseline,
                    changed,
                ),
                contains="ConcurrentDrift",
            )
    elif slug == "source-count-preserving-swap":
        first_rows = [
            v2_synthetic_target(0),
            v2_synthetic_target(1),
        ]
        first = capture_state(
            first_rows,
            lint={"warnings": [{"id": first_rows[0]["id"], "kind": "A"}]},
        )
        field_swap = json.loads(json.dumps(first_rows))
        field_swap[0]["title"], field_swap[1]["title"] = (
            field_swap[1]["title"],
            field_swap[0]["title"],
        )
        membership_swap = [
            first_rows[0],
            v2_synthetic_target(2),
        ]
        lint_swap = capture_state(
            first_rows,
            lint={"warnings": [{"id": first_rows[1]["id"], "kind": "A"}]},
        )
        for name, changed in (
            ("field-swap", capture_state(field_swap)),
            ("drop-add", capture_state(membership_swap)),
            ("warning-swap", lint_swap),
        ):
            checks.refuses(
                f"equal-count-{name}",
                InputRefused,
                lambda changed=changed: v2_validate_capture_pair(first, changed),
                contains="ConcurrentDrift",
            )
        duplicate_rows = [first_rows[0], first_rows[0]]
        checks.refuses(
            "duplicate-membership",
            InputRefused,
            lambda: capture_state(duplicate_rows),
            contains="duplicate",
        )
    elif slug == "source-cross-run-drift-scope":
        selected = v2_synthetic_target(1)
        unrelated_a = v2_synthetic_target(2)
        unrelated_b = dict(unrelated_a)
        unrelated_b["title"] = "unrelated changed title"
        first_capture = capture_state([selected, unrelated_a])
        second_capture = capture_state([selected, unrelated_b])
        checks.check(
            "unrelated-drift-stable-selected-key",
            first_capture["semantic_root"] != second_capture["semantic_root"]
            and v2_target_root(selected) == v2_target_root(dict(selected)),
            expected=v2_target_root(selected),
            observed=v2_target_root(selected),
        )
        selected_mutations = {
            "field": ("field_roots", {"description": text_root("changed")}),
            "membership": ("missing_sections", ["## Success Criteria"]),
            "dependency": (
                "dependencies",
                [{"id": "changed-neighbor", "type": "blocks"}],
            ),
        }
        for name, (field, value) in selected_mutations.items():
            changed = json.loads(json.dumps(selected))
            if field == "field_roots":
                changed[field] = {
                    **changed[field],
                    **value,
                }
            else:
                changed[field] = value
            relevant_root_changed = (
                v2_dependency_neighborhood_root(changed)
                != v2_dependency_neighborhood_root(selected)
                if name == "dependency"
                else v2_target_root(changed) != v2_target_root(selected)
            )
            checks.check(
                f"selected-{name}-invalidates-affected-key",
                relevant_root_changed,
                expected="different selected key",
                observed={
                    "before": (
                        v2_dependency_neighborhood_root(selected)
                        if name == "dependency"
                        else v2_target_root(selected)
                    ),
                    "after": (
                        v2_dependency_neighborhood_root(changed)
                        if name == "dependency"
                        else v2_target_root(changed)
                    ),
                },
            )
        inventory = v2_synthetic_inventory([selected])
        base_receipts = v2_synthetic_receipts(
            inventory,
            declared=True,
        )
        base_authority = v2_derive_authority(
            inventory,
            base_receipts,
            current_br_version="0.2.19",
        )
        base_plan, _ = v2_build_review_plan(
            inventory,
            base_authority,
            max_targets=1,
        )
        base_child_key = base_plan["children"][0]["child_key"]
        for name, field, value in (
            ("owner", "implementation_owner", "changed-owner"),
            ("consumer", "terminal_consumer", "changed-consumer"),
        ):
            changed_receipts = json.loads(json.dumps(base_receipts))
            changed_receipt = changed_receipts["receipts"][0]
            changed_receipt[field] = value
            changed_receipt.pop("semantic_root", None)
            changed_receipts["receipts"][0] = v2_rooted(changed_receipt)
            changed_receipts.pop("semantic_root", None)
            changed_receipts = v2_rooted(changed_receipts)
            changed_authority = v2_derive_authority(
                inventory,
                changed_receipts,
                current_br_version="0.2.19",
            )
            changed_plan, _ = v2_build_review_plan(
                inventory,
                changed_authority,
                max_targets=1,
            )
            changed_child_key = changed_plan["children"][0]["child_key"]
            checks.check(
                f"selected-{name}-receipt-invalidates-child-key",
                changed_child_key != base_child_key,
                expected="different child key",
                observed={
                    "before": base_child_key,
                    "after": changed_child_key,
                },
            )
    else:
        raise EvidenceFailed(f"no source executor for {slug}")
    return checks.finish(parameters)


def v2_execute_authority_cases(
    case: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    slug = str(case["id"]).removeprefix("template-lint-v2.")
    checks = V2CheckCollector(case["id"], case["assertion_id"])
    parameters = {"family": "row-authority", "slug": slug}
    base_inventory = v2_synthetic_inventory([v2_synthetic_target(0)])
    if slug == "row-complete-planned-child":
        authority = v2_synthetic_authority(base_inventory)
        plan, _ = v2_build_review_plan(
            base_inventory,
            authority,
            max_targets=10,
        )
        child = plan["children"][0]
        required = set(manifest["row_contract"]["planned_child_required_fields"])
        aliases = {
            "acceptance": child.get("acceptance"),
            "target_implementation_estimated_minutes": child.get(
                "target_implementation_estimated_minutes"
            ),
            "external_authority_receipt_root": child.get(
                "external_authority_receipt_root"
            ),
            "external_authority_verdict": child.get(
                "external_authority_verdict"
            ),
            "conditional_write_receipt_root": child.get(
                "conditional_write_receipt_root"
            ),
            "conditional_write_verdict": child.get(
                "conditional_write_verdict"
            ),
        }
        checks.check(
            "planned-child-required-fields",
            required.issubset(set(child) | set(aliases)),
            expected=sorted(required),
            observed=sorted(set(child) | set(aliases)),
        )
        v2_validate_child_payload(child)
        checks.check(
            "planned-child-transport-valid",
            True,
            expected="valid",
            observed=child["semantic_root"],
        )
        for field in ("title", "target_ids", "intended_generated_edges"):
            mutated = dict(child)
            mutated.pop(field)
            mutated = v2_rooted(mutated)
            checks.refuses(
                f"missing-{field}",
                InputRefused,
                lambda mutated=mutated: v2_validate_child_payload(mutated),
                contains="non-closed schema",
            )
    elif slug == "row-dispositions-distinct":
        dispositions = [
            "MALFORMED_OR_WRONG_TYPE",
            "ROLLUP_CHILD_SET_GAP",
            "OWNER_REVIEW_REQUIRED",
            "SUBSTANTIVE_SEMANTIC_OMISSION",
            "SECTION_NAME_ONLY",
            "HISTORICAL_IMMUTABLE_REVIEW",
        ]
        roots: list[str] = []
        for index, disposition in enumerate(dispositions):
            inventory = v2_synthetic_inventory(
                [
                    v2_synthetic_target(
                        index,
                        disposition=disposition,
                        status="closed"
                        if disposition == "HISTORICAL_IMMUTABLE_REVIEW"
                        else "open",
                    )
                ]
            )
            authority = v2_synthetic_authority(inventory)
            if disposition == "HISTORICAL_IMMUTABLE_REVIEW":
                description, acceptance, design, notes = v2_child_text(
                    inventory["rows"],
                    authority["decisions"],
                    oversize=False,
                )
                roots.append(
                    semantic_root([description, acceptance, design, notes])
                )
            else:
                plan, _ = v2_build_review_plan(
                    inventory,
                    authority,
                    max_targets=10,
                )
                roots.append(plan["children"][0]["description_root"])
        checks.check(
            "disposition-workflows-distinct",
            len(set(roots)) == len(dispositions),
            expected=len(dispositions),
            observed=len(set(roots)),
        )
    elif slug == "authority-candidate-provenance":
        issue = {
            "authority_domain": "domain:solver",
            "parent": "parent-1",
            "labels": ["crate:fs-la", "unrelated"],
            "mapping": {"paths": ["crates/fs-la/src/lib.rs"]},
        }
        candidates = v2_domain_candidates(issue)
        checks.check(
            "candidate-sources-rooted",
            len(candidates) == 4
            and all(
                candidate["provenance_root"].startswith("sha256-v1:")
                and candidate["falsifier"]
                for candidate in candidates
            ),
            expected=4,
            observed=candidates,
        )
        authority = v2_synthetic_authority(base_inventory)
        checks.check(
            "candidates-do-not-promote",
            authority["decisions"][0]["readiness"] == "REVIEW_ONLY",
            expected="REVIEW_ONLY",
            observed=authority["decisions"][0]["readiness"],
        )
    elif slug == "authority-self-report-declared":
        authority = v2_synthetic_authority(
            base_inventory,
            declared=True,
        )
        decision = authority["decisions"][0]
        checks.check(
            "self-report-ceiling",
            decision["readiness"] == "DECLARED_READY"
            and decision["remediation_route"] == "ANALYSIS_ONLY",
            expected=("DECLARED_READY", "ANALYSIS_ONLY"),
            observed=(decision["readiness"], decision["remediation_route"]),
        )
    elif slug == "authority-external-receipt":
        for verdict in (
            "NODATA",
            "INVALID",
            "EXPIRED",
            "REVOKED",
            "CONFLICTED",
        ):
            authority = v2_synthetic_authority(
                base_inventory,
                external_verdict=verdict,
            )
            checks.check(
                f"external-{verdict.lower()}-not-verified",
                authority["decisions"][0]["external_authority"]["verified"]
                is False,
                expected=False,
                observed=authority["decisions"][0]["external_authority"][
                    "verified"
                ],
            )
        valid = v2_synthetic_authority(
            base_inventory,
            declared=True,
            external_verdict="VALID",
            version="0.3.0",
            allow_mechanical_fixture=True,
        )
        checks.check(
            "independent-fixture-external-valid",
            valid["decisions"][0]["external_authority"]["verified"] is True
            and valid["decisions"][0]["readiness"] == "DECLARED_READY",
            expected=(True, "DECLARED_READY"),
            observed=(
                valid["decisions"][0]["external_authority"]["verified"],
                valid["decisions"][0]["readiness"],
            ),
        )
    elif slug == "capability-valid-missing-failed":
        for verdict in ("NODATA", "FAILED"):
            authority = v2_synthetic_authority(
                base_inventory,
                conditional_verdict=verdict,
            )
            checks.check(
                f"capability-{verdict.lower()}-false",
                authority["decisions"][0]["conditional_write_capability"][
                    "verified"
                ]
                is False,
                expected=False,
                observed=authority["decisions"][0][
                    "conditional_write_capability"
                ]["verified"],
            )
        valid_only = v2_synthetic_authority(
            base_inventory,
            conditional_verdict="VALID",
            version="0.3.0",
            allow_mechanical_fixture=True,
        )
        checks.check(
            "capability-alone-not-mechanical",
            valid_only["decisions"][0]["readiness"] != "MECHANICALLY_APPLY_ELIGIBLE",
            expected="not mechanical",
            observed=valid_only["decisions"][0]["readiness"],
        )
    elif slug == "capability-lying-nonatomic-version":
        for verdict in (
            "LYING",
            "NON_ATOMIC",
            "VERSION_MISMATCHED",
            "CONFLICTED",
            "FAILED",
        ):
            authority = v2_synthetic_authority(
                base_inventory,
                conditional_verdict=verdict,
                version="0.3.0",
                allow_mechanical_fixture=True,
            )
            checks.check(
                f"{verdict.lower()}-not-conditional",
                authority["decisions"][0]["conditional_write_capability"][
                    "verified"
                ]
                is False,
                expected=False,
                observed=authority["decisions"][0][
                    "conditional_write_capability"
                ]["verified"],
            )
    elif slug == "readiness-two-gate-truth-table":
        observed: dict[str, str] = {}
        for external in ("NODATA", "VALID"):
            for conditional in ("NODATA", "VALID"):
                for gate in ("NODATA", "VALID"):
                    key = f"{external}/{conditional}/{gate}"
                    authority = v2_synthetic_authority(
                        base_inventory,
                        declared=True,
                        external_verdict=external,
                        conditional_verdict=conditional,
                        gate_verdict=gate,
                        version="0.3.0",
                        allow_mechanical_fixture=True,
                    )
                    observed[key] = authority["decisions"][0]["readiness"]
        checks.check(
            "both-gates-plus-admission-only",
            observed["VALID/VALID/VALID"]
            == "MECHANICALLY_APPLY_ELIGIBLE"
            and all(
                value != "MECHANICALLY_APPLY_ELIGIBLE"
                for key, value in observed.items()
                if key != "VALID/VALID/VALID"
            ),
            expected="only VALID/VALID/VALID mechanical",
            observed=observed,
        )
        current = v2_synthetic_authority(
            base_inventory,
            declared=True,
            external_verdict="VALID",
            conditional_verdict="VALID",
            gate_verdict="VALID",
            version="0.2.19",
            allow_mechanical_fixture=True,
        )
        checks.check(
            "current-tool-nodata",
            current["decisions"][0]["readiness"] == "DECLARED_READY"
            and current["current_br_conditional_write_capability"]
            == "AUTOMATION_NODATA",
            expected=("DECLARED_READY", "AUTOMATION_NODATA"),
            observed=(
                current["decisions"][0]["readiness"],
                current["current_br_conditional_write_capability"],
            ),
        )
    elif slug == "route-manual-distinct":
        manual = v2_synthetic_authority(
            base_inventory,
            declared=True,
            manual=True,
        )["decisions"][0]
        automated = v2_synthetic_authority(
            base_inventory,
            declared=True,
            external_verdict="VALID",
            conditional_verdict="VALID",
            gate_verdict="VALID",
            version="0.3.0",
            allow_mechanical_fixture=True,
        )["decisions"][0]
        analysis = v2_synthetic_authority(base_inventory)["decisions"][0]
        checks.check(
            "three-routes-distinct",
            {
                analysis["remediation_route"],
                manual["remediation_route"],
                automated["remediation_route"],
            }
            == set(V2_REMEDIATION_ROUTES),
            expected=sorted(V2_REMEDIATION_ROUTES),
            observed=sorted(
                {
                    analysis["remediation_route"],
                    manual["remediation_route"],
                    automated["remediation_route"],
                }
            ),
        )
        checks.check(
            "manual-no-cas-claims",
            manual["manual_authorization"]["mechanical_authority"] is False
            and manual["manual_authorization"][
                "no_cas_no_clobber_no_exactly_once"
            ]
            is True,
            expected=True,
            observed=manual["manual_authorization"],
        )
    elif slug == "route-deferred-reactivation":
        deferred_inventory = v2_synthetic_inventory(
            [v2_synthetic_target(0, status="deferred")]
        )
        authority = v2_synthetic_authority(
            deferred_inventory,
            declared=True,
            manual=True,
            external_verdict="VALID",
            conditional_verdict="VALID",
            gate_verdict="VALID",
            version="0.3.0",
            allow_mechanical_fixture=True,
        )
        decision = authority["decisions"][0]
        checks.check(
            "deferred-always-analysis",
            decision["remediation_route"] == "ANALYSIS_ONLY"
            and decision["deferred_apply_prohibition"] is True,
            expected=("ANALYSIS_ONLY", True),
            observed=(
                decision["remediation_route"],
                decision["deferred_apply_prohibition"],
            ),
        )
    elif slug == "route-active-owner-conflict":
        for status, assignee in (("in_progress", ""), ("open", "agent")):
            inventory = v2_synthetic_inventory(
                [v2_synthetic_target(0, status=status, assignee=assignee)]
            )
            authority = v2_synthetic_authority(
                inventory,
                declared=True,
                manual=True,
            )
            plan, _ = v2_build_review_plan(
                inventory,
                authority,
                max_targets=10,
            )
            checks.check(
                f"active-{status}-{bool(assignee)}",
                authority["decisions"][0]["remediation_route"] == "ANALYSIS_ONLY"
                and plan["children"][0]["desired_status"] == "deferred",
                expected=("ANALYSIS_ONLY", "deferred"),
                observed=(
                    authority["decisions"][0]["remediation_route"],
                    plan["children"][0]["desired_status"],
                ),
            )
    elif slug == "route-no-target-mutation":
        before = semantic_root(base_inventory["rows"])
        authority = v2_synthetic_authority(
            base_inventory,
            declared=True,
            manual=True,
        )
        v2_build_review_plan(base_inventory, authority, max_targets=10)
        after = semantic_root(base_inventory["rows"])
        checks.check(
            "planning-target-unchanged",
            before == after,
            expected=before,
            observed=after,
        )
    else:
        raise EvidenceFailed(f"no authority executor for {slug}")
    return checks.finish(parameters)


def v2_execute_packing_cases(
    case: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    del manifest
    slug = str(case["id"]).removeprefix("template-lint-v2.")
    checks = V2CheckCollector(case["id"], case["assertion_id"])
    parameters = {"family": "packing", "slug": slug}
    if slug == "packing-hard-keys":
        base = v2_synthetic_target(0)
        authority = v2_synthetic_authority(
            v2_synthetic_inventory([base])
        )["decisions"][0]
        base_vector = v2_hard_vector(base, authority)
        mutations = {
            "priority": {"priority": 3},
            "status": {"status": "blocked"},
            "disposition": {"disposition": "OWNER_REVIEW_REQUIRED"},
            "type": {"type": "bug"},
            "missing": {"missing_sections": ["## Success Criteria"]},
        }
        for name, changes in mutations.items():
            target = dict(base)
            target.update(changes)
            checks.check(
                f"hard-key-{name}",
                v2_hard_vector(target, authority) != base_vector,
                expected="different vector",
                observed=v2_hard_vector(target, authority),
            )
        changed_authority = dict(authority)
        changed_authority["readiness"] = "DECLARED_READY"
        checks.check(
            "hard-key-readiness",
            v2_hard_vector(base, changed_authority) != base_vector,
            expected="different vector",
            observed=v2_hard_vector(base, changed_authority),
        )
        changed_authority = dict(authority)
        changed_authority["remediation_route"] = "MANUAL_BR_REVIEW"
        checks.check(
            "hard-key-remediation-route",
            v2_hard_vector(base, changed_authority) != base_vector,
            expected="different vector",
            observed=v2_hard_vector(base, changed_authority),
        )
    elif slug == "packing-compatible-receipt":
        plan, _, _ = v2_plan_for_count(2, max_targets=10)
        checks.check(
            "compatible-merge",
            len(plan["children"]) == 1
            and len(plan["children"][0]["target_ids"]) == 2,
            expected={"children": 1, "targets": 2},
            observed={
                "children": len(plan["children"]),
                "targets": len(plan["children"][0]["target_ids"]),
            },
        )
        checks.check(
            "receipt-rationale-falsifier",
            bool(plan["children"][0]["compatibility"]["rationale"])
            and bool(plan["children"][0]["compatibility"]["falsifier"]),
            expected=True,
            observed=plan["children"][0]["compatibility"],
        )
    elif slug == "packing-generic-merge-refused":
        inventory = v2_synthetic_inventory(
            [v2_synthetic_target(0), v2_synthetic_target(1)]
        )
        authority = v2_synthetic_authority(inventory)
        plan, _ = v2_build_review_plan(
            inventory,
            authority,
            max_targets=10,
        )
        checks.check(
            "no-receipt-singletons",
            len(plan["children"]) == 2
            and all(len(child["target_ids"]) == 1 for child in plan["children"]),
            expected=2,
            observed=len(plan["children"]),
        )
    elif slug == "packing-rationale-falsifier":
        singleton, _, _ = v2_plan_for_count(1, max_targets=10)
        merged, _, _ = v2_plan_for_count(2, max_targets=10)
        checks.check(
            "singleton-rationale",
            bool(singleton["children"][0]["compatibility"]["rationale"])
            and bool(singleton["children"][0]["compatibility"]["falsifier"]),
            expected=True,
            observed=singleton["children"][0]["compatibility"],
        )
        checks.check(
            "merge-witness-rationale",
            all(
                witness.get("rationale") and witness.get("falsifier")
                for witness in merged["packing_witnesses"]
            ),
            expected=True,
            observed=merged["packing_witnesses"],
        )
    elif slug == "packing-target-cap-10-25":
        for cap, population, expected_children in (
            (10, 10, 1),
            (10, 11, 2),
            (25, 25, 1),
            (25, 26, 2),
        ):
            plan, _, _ = v2_plan_for_count(population, max_targets=cap)
            checks.check(
                f"cap-{cap}-population-{population}",
                len(plan["children"]) == expected_children,
                expected=expected_children,
                observed=len(plan["children"]),
            )
        inventory = v2_synthetic_inventory([v2_synthetic_target(0)])
        authority = v2_synthetic_authority(inventory)
        checks.refuses(
            "configured-cap-26",
            UsageRefused,
            lambda: v2_build_review_plan(
                inventory,
                authority,
                max_targets=26,
            ),
            contains="1..25",
        )
    elif slug == "packing-minutes-479-480-481":
        for minutes in (479, 480):
            inventory = v2_synthetic_inventory(
                [v2_synthetic_target(0, review_minutes=minutes)]
            )
            authority = v2_synthetic_authority(inventory)
            plan, optional = v2_build_review_plan(
                inventory,
                authority,
                max_targets=10,
            )
            checks.check(
                f"minutes-{minutes}-bounded",
                not optional
                and plan["children"][0]["review_minutes"] == minutes,
                expected=minutes,
                observed=plan["children"][0]["review_minutes"],
            )
        inventory = v2_synthetic_inventory(
            [v2_synthetic_target(0, review_minutes=481)]
        )
        authority = v2_synthetic_authority(inventory)
        plan, optional = v2_build_review_plan(
            inventory,
            authority,
            max_targets=10,
        )
        checks.check(
            "minutes-481-oversize",
            plan["children"][0]["review_minutes"] == 480
            and plan["work"]["total_review_minutes"] == 481
            and len(optional) == 1,
            expected={"bounded": 480, "underlying": 481, "optional": 1},
            observed={
                "bounded": plan["children"][0]["review_minutes"],
                "underlying": plan["work"]["total_review_minutes"],
                "optional": len(optional),
            },
        )
    elif slug == "packing-transport-cap-boundaries":
        inventory = v2_synthetic_inventory([v2_synthetic_target(0)])
        authority = v2_synthetic_authority(inventory)
        plan, _ = v2_build_review_plan(inventory, authority, max_targets=10)
        child = plan["children"][0]
        v2_validate_child_payload(child)
        checks.check(
            "baseline-transport-valid",
            True,
            expected="valid",
            observed=child["semantic_root"],
        )
        oversized = dict(child)
        description = dict(child["description_file_artifact"])
        description["content"] = "x" * (V2_CHILD_DESCRIPTION_CAP + 1)
        description["bytes"] = len(description["content"])
        description["root"] = text_root(description["content"])
        oversized["description_file_artifact"] = description
        oversized["description_root"] = description["root"]
        oversized = v2_rooted(oversized)
        checks.refuses(
            "description-cap-plus-one",
            EvidenceFailed,
            lambda: v2_validate_child_payload(oversized),
            contains="description",
        )
    elif slug == "packing-exact-small-optimum":
        for count in range(1, 14):
            plan, _, _ = v2_plan_for_count(count, max_targets=13)
            witness = plan["packing_witnesses"][0]
            checks.check(
                f"exact-count-{count:02d}",
                witness["algorithm"] == "EXACT_SUBSET_PARTITION"
                and witness["exact_optimality_claim"] is True
                and witness["objective"]["child_count"] == 1,
                expected={"algorithm": "EXACT_SUBSET_PARTITION", "children": 1},
                observed={
                    "algorithm": witness["algorithm"],
                    "children": witness["objective"]["child_count"],
                },
            )
        minute_loads = (240, 240, 160, 160, 80, 80)
        byte_loads = (900, 100, 800, 200, 700, 300)
        inventory = v2_synthetic_inventory(
            [
                v2_synthetic_target(
                    index,
                    review_minutes=minutes,
                    retained_payload_bytes=byte_loads[index],
                )
                for index, minutes in enumerate(minute_loads)
            ]
        )
        authority = v2_synthetic_authority(inventory, compatible=True)
        authority_by_id = {
            row["target_id"]: row for row in authority["decisions"]
        }
        rows = [
            (target, authority_by_id[target["id"]])
            for target in inventory["rows"]
        ]

        def independent_load(
            group: Sequence[
                tuple[Mapping[str, Any], Mapping[str, Any]]
            ],
        ) -> tuple[int, int, int]:
            return (
                sum(int(decision["review_minutes"]) for _, decision in group),
                sum(int(target["retained_payload_bytes"]) for target, _ in group),
                len(group),
            )

        def independent_key(
            groups: Sequence[
                Sequence[tuple[Mapping[str, Any], Mapping[str, Any]]]
            ],
        ) -> tuple[Any, ...]:
            loads = [independent_load(group) for group in groups]
            return (
                len(groups),
                max(load[0] for load in loads),
                max(load[1] for load in loads),
                max(load[2] for load in loads),
                tuple(sorted(loads, reverse=True)),
                tuple(
                    sorted(
                        tuple(sorted(target["id"] for target, _ in group))
                        for group in groups
                    )
                ),
            )

        independent_best: tuple[Any, ...] | None = None

        def enumerate_partitions(
            index: int,
            groups: list[
                list[tuple[Mapping[str, Any], Mapping[str, Any]]]
            ],
        ) -> None:
            nonlocal independent_best
            if index == len(rows):
                candidate = independent_key(groups)
                if independent_best is None or candidate < independent_best:
                    independent_best = candidate
                return
            row = rows[index]
            for group_index in range(len(groups) + 1):
                if group_index == len(groups):
                    groups.append([row])
                    enumerate_partitions(index + 1, groups)
                    groups.pop()
                    continue
                proposed = [*groups[group_index], row]
                minutes, retained_bytes, target_count = independent_load(proposed)
                if (
                    target_count > 3
                    or minutes > V2_REVIEW_MINUTES_CAP
                    or retained_bytes > V2_CHILD_PAYLOAD_CAP
                ):
                    continue
                groups[group_index].append(row)
                enumerate_partitions(index + 1, groups)
                groups[group_index].pop()

        enumerate_partitions(0, [])
        actual_bins, actual_witness = v2_pack_group(rows, max_targets=3)
        actual_key = independent_key(actual_bins)
        checks.check(
            "independent-enumeration-confirms-varied-load-optimum",
            independent_best is not None
            and actual_key == independent_best
            and actual_witness["objective"]
            == v2_objective_document(actual_bins),
            expected=independent_best,
            observed=actual_key,
        )
    elif slug == "packing-large-gap-witness":
        plan, _, _ = v2_plan_for_count(14, max_targets=10)
        witness = plan["packing_witnesses"][0]
        checks.check(
            "large-witness-complete",
            witness["algorithm"] == "DETERMINISTIC_LPT_VECTOR_PACKING"
            and "lower_bound_vector" in witness
            and "achieved_objective" in witness
            and "gap_witness" in witness
            and witness["exact_optimality_claim"] is False,
            expected=True,
            observed=witness,
        )
    elif slug == "packing-order-invariant":
        targets = [v2_synthetic_target(index) for index in range(8)]
        first_inventory = v2_synthetic_inventory(targets)
        second_inventory = v2_synthetic_inventory(list(reversed(targets)))
        first_authority = v2_synthetic_authority(
            first_inventory,
            compatible=True,
        )
        second_authority = v2_synthetic_authority(
            second_inventory,
            compatible=True,
        )
        first, _ = v2_build_review_plan(
            first_inventory,
            first_authority,
            max_targets=10,
        )
        second, _ = v2_build_review_plan(
            second_inventory,
            second_authority,
            max_targets=10,
        )
        checks.check(
            "input-order-invariant",
            first["semantic_root"] == second["semantic_root"],
            expected=first["semantic_root"],
            observed=second["semantic_root"],
        )
    elif slug == "packing-refinement-imbalance":
        inventory = v2_synthetic_inventory(
            [v2_synthetic_target(index) for index in range(4)]
        )
        merged = v2_synthetic_authority(inventory, compatible=True)
        merged_plan, _ = v2_build_review_plan(
            inventory,
            merged,
            max_targets=10,
        )
        refined = v2_synthetic_authority(inventory, compatible=False)
        refined_plan, _ = v2_build_review_plan(
            inventory,
            refined,
            max_targets=10,
        )
        checks.check(
            "refinement-never-remerges",
            len(merged_plan["children"]) == 1
            and len(refined_plan["children"]) == 4,
            expected=(1, 4),
            observed=(
                len(merged_plan["children"]),
                len(refined_plan["children"]),
            ),
        )
    elif slug == "packing-oversize-escalation":
        inventory = v2_synthetic_inventory(
            [
                v2_synthetic_target(
                    0,
                    retained_payload_bytes=V2_CHILD_PAYLOAD_CAP + 1,
                )
            ]
        )
        authority = v2_synthetic_authority(inventory)
        plan, optional = v2_build_review_plan(
            inventory,
            authority,
            max_targets=10,
        )
        entry = plan["oversize_content"][0]
        full = strict_json_loads(
            next(iter(optional.values())),
            label="synthetic oversize",
            require_canonical=True,
        )
        checks.check(
            "oversize-registry-complete",
            set(entry)
            == {
                "relative_path",
                "media_kind",
                "schema_kind",
                "byte_length",
                "semantic_root",
                "target_roots",
                "field_roots",
                "clause_roots",
                "aggregate_cap_accounting",
            }
            and full["bounded_incomplete_index"]["complete"] is False,
            expected=True,
            observed={"entry": entry, "full_root": full["semantic_root"]},
        )
    else:
        raise EvidenceFailed(f"no packing executor for {slug}")
    return checks.finish(parameters)


def v2_history_fixture(
    manifest: Mapping[str, Any],
    *,
    legacy: bool = False,
    missing_modern_event: bool = False,
    conflicting_close: bool = False,
    candidate_consumers: bool = False,
) -> tuple[
    dict[str, Any],
    dict[str, Any],
    list[dict[str, Any]],
    dict[str, Any],
    dict[str, Any],
]:
    anchor_id = "fixture-anchor"
    anchor_closed_at = "2026-01-02T00:00:00.000000+00:00"
    anchor_raw = {
        "id": anchor_id,
        "title": "Fixture audit anchor",
        "issue_type": "task",
        "status": "closed",
        "priority": 2,
        "description": "Anchor source.",
        "acceptance_criteria": "Exact anchor.",
        "design": "",
        "notes": "See fixture-consumer." if candidate_consumers else "",
        "assignee": "fixture-assignee",
        "labels": [],
        "parent": "",
        "created_at": "2026-01-01T00:00:00+00:00",
        "created_by": "fixture-creator",
        "updated_at": "2026-01-02T00:00:00.002000+00:00",
        "closed_at": anchor_closed_at,
        "close_reason": "Fixture anchor completed.",
        "estimated_minutes": 5,
        "dependencies": [],
        "dependents": (
            [
                {
                    "id": "fixture-consumer",
                    "dependency_type": "blocks",
                    "status": "open",
                    "priority": 2,
                }
            ]
            if candidate_consumers
            else []
        ),
    }
    anchor = v2_full_issue_projection(anchor_raw)
    anchor_events = [
        {
            "id": 2,
            "event_type": "closed",
            "actor": "fixture-closer",
            "timestamp": "2026-01-02T00:00:00.002000+00:00",
            "comment": "Fixture anchor completed.",
        },
        {
            "id": 1,
            "event_type": "status_changed",
            "actor": "fixture-closer",
            "timestamp": "2026-01-02T00:00:00.001000+00:00",
            "old_value": "open",
            "new_value": "closed",
        },
    ]
    if conflicting_close:
        anchor_events = [
            {
                "id": 3,
                "event_type": "closed",
                "actor": "other-closer",
                "timestamp": "2026-01-02T00:00:00.003000+00:00",
                "comment": "Fixture anchor completed.",
            },
            *anchor_events,
        ]
    anchor_audit = v2_normalize_audit_document(
        {"issue_id": anchor_id, "events": anchor_events},
        issue_id=anchor_id,
    )
    all_issues = [anchor]
    audit_documents = [anchor_audit]
    target_id = anchor_id
    if missing_modern_event:
        modern_id = "fixture-modern-missing-close"
        modern_raw = {
            **anchor_raw,
            "id": modern_id,
            "title": "Fixture modern close without audit evidence",
            "closed_at": "2026-01-03T00:00:00+00:00",
            "updated_at": "2026-01-03T00:00:00+00:00",
            "close_reason": "Modern fixture completed.",
            "dependents": [],
        }
        modern_issue = v2_full_issue_projection(modern_raw)
        modern_audit = v2_normalize_audit_document(
            {"issue_id": modern_id, "events": []},
            issue_id=modern_id,
        )
        all_issues.append(modern_issue)
        audit_documents.append(modern_audit)
        target_id = modern_id
    if legacy:
        legacy_id = "fixture-legacy"
        legacy_raw = {
            **anchor_raw,
            "id": legacy_id,
            "title": "Fixture legacy close",
            "closed_at": "2026-01-01T12:00:00+00:00",
            "close_reason": "Legacy fixture completed.",
            "updated_at": "2026-01-03T00:00:00+00:00",
            "notes": "Later non-close audit metadata.",
            "dependents": [],
        }
        legacy_issue = v2_full_issue_projection(legacy_raw)
        legacy_audit = v2_normalize_audit_document(
            {
                "issue_id": legacy_id,
                "events": [
                    {
                        "id": 4,
                        "event_type": "label_added",
                        "actor": "later-agent",
                        "timestamp": "2026-01-03T00:00:00+00:00",
                        "comment": "Added label after legacy close.",
                    }
                ],
            },
            issue_id=legacy_id,
        )
        all_issues.append(legacy_issue)
        audit_documents.append(legacy_audit)
        target_id = legacy_id
    target = v2_synthetic_target(
        0,
        issue_id=target_id,
        status="closed",
        disposition="HISTORICAL_IMMUTABLE_REVIEW",
    )
    issue_by_id = {row["id"]: row for row in all_issues}
    target = dict(target)
    target["field_roots"] = issue_by_id[target_id]["field_roots"]
    target["target_root"] = semantic_root(
        {
            "id": target_id,
            "field_roots": target["field_roots"],
            "status": "closed",
        }
    )
    target = v2_rooted(target)
    inventory = v2_synthetic_inventory([target])
    authority = v2_synthetic_authority(inventory)
    audit_capture = v2_rooted(
        {
            "capture_count": 2,
            "worker_bound": 8,
            "issue_ids": sorted(row["issue_id"] for row in audit_documents),
            "documents": sorted(
                audit_documents, key=lambda row: row["issue_id"]
            ),
            "command_receipts": [],
            "raw_stream_bodies_retained": False,
            "no_claim": "synthetic audit capture",
        }
    )
    anchor_timestamp = v2_parse_timestamp(
        anchor_closed_at,
        label="fixture anchor",
    )
    legacy_rows = sorted(
        [
            {
                "id": issue["id"],
                "closed_at": issue["closed_at"],
                "close_reason_root": text_root(issue["close_reason"]),
            }
            for issue in all_issues
            if v2_parse_timestamp(
                issue["closed_at"],
                label=f"fixture {issue['id']}",
            )
            < anchor_timestamp
        ],
        key=lambda row: row["id"],
    )
    contract = dict(manifest["history_contract"])
    contract.update(
        {
            "legacy_coverage_anchor_issue": anchor_id,
            "legacy_coverage_anchor_closed_at": anchor_closed_at,
            "legacy_coverage_anchor_status_event_id": 1,
            "legacy_coverage_anchor_close_event_id": 2,
            "legacy_coverage_count": len(legacy_rows),
            "legacy_coverage_rows_root": semantic_root(legacy_rows),
            "known_pair_max_skew_ms": 1_000,
        }
    )
    return inventory, authority, all_issues, audit_capture, contract


def v2_execute_history_cases(
    case: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    slug = str(case["id"]).removeprefix("template-lint-v2.")
    checks = V2CheckCollector(case["id"], case["assertion_id"])
    parameters = {"family": "history-zero-movement", "slug": slug}
    if slug == "history-known-closer":
        inventory, authority, issues, audit, contract = v2_history_fixture(
            manifest
        )
        history = v2_build_history(
            source_root=semantic_root({"fixture": "known"}),
            inventory=inventory,
            authority=authority,
            all_issues=issues,
            audit_capture=audit,
            history_contract=contract,
        )
        row = history["rows"][0]
        checks.check(
            "known-closer-exact",
            row["closer_state"] == "KNOWN"
            and row["close_actor"] == "fixture-closer"
            and row["close_actor_source"] == "br.audit.log.closed.actor"
            and row["audit_event_root"],
            expected=("KNOWN", "fixture-closer"),
            observed=(
                row["closer_state"],
                row["close_actor"],
                row["audit_event_root"],
            ),
        )
    elif slug == "history-legacy-closer":
        inventory, authority, issues, audit, contract = v2_history_fixture(
            manifest,
            legacy=True,
        )
        history = v2_build_history(
            source_root=semantic_root({"fixture": "legacy"}),
            inventory=inventory,
            authority=authority,
            all_issues=issues,
            audit_capture=audit,
            history_contract=contract,
        )
        row = history["rows"][0]
        checks.check(
            "legacy-actor-nodata",
            row["closer_state"] == "LEGACY_UNAVAILABLE"
            and row["close_actor"] is None
            and row["close_actor_source"] == contract["legacy_coverage_rows_root"],
            expected=("LEGACY_UNAVAILABLE", None),
            observed=(row["closer_state"], row["close_actor"]),
        )
    elif slug == "history-unknown-conflicted-closer":
        args = v2_history_fixture(manifest, conflicting_close=True)
        checks.refuses(
            "conflicting-closers-refused",
            EvidenceFailed,
            lambda: v2_build_history(
                source_root=semantic_root({"fixture": "conflict"}),
                inventory=args[0],
                authority=args[1],
                all_issues=args[2],
                audit_capture=args[3],
                history_contract=args[4],
            ),
            contains="conflicted",
        )
        missing_actor = {
            "issue_id": "x",
            "events": [
                {
                    "id": 2,
                    "event_type": "closed",
                    "timestamp": "2026-01-02T00:00:00+00:00",
                },
                {
                    "id": 1,
                    "event_type": "status_changed",
                    "actor": "creator",
                    "timestamp": "2026-01-01T23:59:59+00:00",
                    "new_value": "closed",
                },
            ],
        }
        normalized = v2_normalize_audit_document(
            missing_actor,
            issue_id="x",
        )
        checks.check(
            "missing-actor-retained-empty",
            "actor" not in normalized["events"][0],
            expected="NoData actor",
            observed=normalized["events"][0],
        )
    elif slug == "history-timestamp-reason":
        valid_values = (
            "2026-01-01T00:00:00Z",
            "2026-01-01T00:00:00.123456+00:00",
            "2026-01-01T01:00:00+01:00",
        )
        for index, value in enumerate(valid_values):
            parsed = v2_parse_timestamp(value, label="fixture")
            checks.check(
                f"valid-timestamp-{index}",
                parsed.tzinfo is not None,
                expected=True,
                observed=str(parsed),
            )
        for index, value in enumerate(("", "2026-01-01", "not-time")):
            checks.refuses(
                f"invalid-timestamp-{index}",
                EvidenceFailed,
                lambda value=value: v2_parse_timestamp(value, label="fixture"),
                contains="timestamp",
            )
        inventory, authority, issues, audit, contract = v2_history_fixture(
            manifest,
            missing_modern_event=True,
        )
        selected_id = inventory["rows"][0]["id"]
        selected_index = next(
            index for index, issue in enumerate(issues) if issue["id"] == selected_id
        )
        issues[selected_index] = dict(issues[selected_index])
        issues[selected_index]["close_reason"] = ""
        checks.refuses(
            "missing-reason-refused",
            EvidenceFailed,
            lambda: v2_build_history(
                source_root=semantic_root({"fixture": "reason"}),
                inventory=inventory,
                authority=authority,
                all_issues=issues,
                audit_capture=audit,
                history_contract=contract,
            ),
            contains="close metadata",
        )
    elif slug == "history-missing-audit-event":
        args = v2_history_fixture(manifest, missing_modern_event=True)
        checks.refuses(
            "modern-missing-close-refused",
            EvidenceFailed,
            lambda: v2_build_history(
                source_root=semantic_root({"fixture": "missing"}),
                inventory=args[0],
                authority=args[1],
                all_issues=args[2],
                audit_capture=args[3],
                history_contract=args[4],
            ),
            contains="conflicted",
        )
    elif slug == "history-duplicate-conflicting-order":
        duplicate = {
            "issue_id": "x",
            "events": [
                {
                    "id": 1,
                    "event_type": "closed",
                    "actor": "a",
                    "timestamp": "2026-01-01T00:00:00Z",
                },
                {
                    "id": 1,
                    "event_type": "closed",
                    "actor": "a",
                    "timestamp": "2026-01-01T00:00:00Z",
                },
            ],
        }
        checks.refuses(
            "duplicate-event-id",
            EvidenceFailed,
            lambda: v2_normalize_audit_document(duplicate, issue_id="x"),
            contains="duplicate",
        )
        reordered = {
            "issue_id": "x",
            "events": [
                {
                    "id": 1,
                    "event_type": "created",
                    "actor": "a",
                    "timestamp": "2026-01-01T00:00:00Z",
                },
                {
                    "id": 2,
                    "event_type": "closed",
                    "actor": "a",
                    "timestamp": "2026-01-02T00:00:00Z",
                },
            ],
        }
        checks.refuses(
            "reordered-event-id",
            EvidenceFailed,
            lambda: v2_normalize_audit_document(reordered, issue_id="x"),
            contains="order",
        )
    elif slug == "history-consumer-candidates":
        inventory, authority, issues, audit, contract = v2_history_fixture(
            manifest,
            candidate_consumers=True,
        )
        history = v2_build_history(
            source_root=semantic_root({"fixture": "consumers"}),
            inventory=inventory,
            authority=authority,
            all_issues=issues,
            audit_capture=audit,
            history_contract=contract,
        )
        row = history["rows"][0]
        checks.check(
            "candidate-not-reviewed",
            "fixture-consumer" in row["candidate_consumers"]
            and row["reviewed_consumers"] == [],
            expected={"candidate": "fixture-consumer", "reviewed": []},
            observed={
                "candidate": row["candidate_consumers"],
                "reviewed": row["reviewed_consumers"],
            },
        )
    elif slug == "history-immutable-relations":
        inventory, authority, issues, audit, contract = v2_history_fixture(
            manifest
        )
        before = semantic_root(issues)
        v2_build_history(
            source_root=semantic_root({"fixture": "immutable"}),
            inventory=inventory,
            authority=authority,
            all_issues=issues,
            audit_capture=audit,
            history_contract=contract,
        )
        after = semantic_root(issues)
        checks.check(
            "closed-source-unchanged",
            before == after,
            expected=before,
            observed=after,
        )
    elif slug == "zero-all-priority-status-cells":
        statuses = ["open", "in_progress", "blocked", "deferred", "closed"]
        targets = [
            v2_synthetic_target(
                priority * 5 + status_index,
                priority=priority,
                status=status,
            )
            for priority in range(5)
            for status_index, status in enumerate(statuses)
        ]
        inventory = v2_synthetic_inventory(targets)
        zero = v2_build_zero_sets(
            source_root=semantic_root({"fixture": "cells"}),
            inventory=inventory,
            prior_campaign=v2_empty_prior_campaign(),
        )
        checks.check(
            "all-25-cells-covered",
            len(zero["cells"]) == 25
            and zero["counts"]["issues"] == 25
            and zero["counts"]["zero_receipts"] == 0,
            expected={"cells": 25, "issues": 25, "zero": 0},
            observed=zero["counts"],
        )
        sparse_inventory = v2_synthetic_inventory(
            [v2_synthetic_target(100, priority=2, status="open")]
        )
        sparse_source_root = semantic_root({"fixture": "sparse-cells"})
        sparse = v2_build_zero_sets(
            source_root=sparse_source_root,
            inventory=sparse_inventory,
            prior_campaign=v2_empty_prior_campaign(),
        )
        sparse_zero_receipts = [
            cell["zero_receipt"]
            for cell in sparse["cells"]
            if cell["zero_receipt"] is not None
        ]
        sparse_roots_valid = True
        for receipt in sparse_zero_receipts:
            try:
                verify_semantic_root(
                    receipt,
                    label="sparse zero receipt",
                )
            except HarnessError:
                sparse_roots_valid = False
        checks.check(
            "sparse-24-zero-receipts-rooted",
            len(sparse_zero_receipts) == 24
            and sparse_roots_valid
            and all(
                receipt["count"] == 0
                and receipt["source_root"] == sparse_source_root
                and receipt["inventory_root"]
                == sparse_inventory["semantic_root"]
                for receipt in sparse_zero_receipts
            ),
            expected={"zero_receipts": 24, "issues": 1},
            observed=sparse["counts"],
        )
        empty_inventory = v2_synthetic_inventory([])
        empty = v2_build_zero_sets(
            source_root=semantic_root({"fixture": "empty-cells"}),
            inventory=empty_inventory,
            prior_campaign=v2_empty_prior_campaign(),
        )
        checks.check(
            "empty-all-cells-have-zero-receipts",
            empty["counts"]
            == {
                "cells": 25,
                "zero_receipts": 25,
                "issues": 0,
                "movements": 0,
            }
            and all(cell["zero_receipt"] is not None for cell in empty["cells"]),
            expected={"cells": 25, "zero_receipts": 25, "issues": 0},
            observed=empty["counts"],
        )
    elif slug in {
        "movement-cross-lane",
        "movement-live-to-history",
        "movement-history-to-live",
    }:
        transitions = {
            "movement-cross-lane": ("open", 1, "blocked", 2, "MOVED_TO_LANE"),
            "movement-live-to-history": (
                "open",
                1,
                "closed",
                1,
                "MOVED_TO_HISTORY",
            ),
            "movement-history-to-live": (
                "closed",
                1,
                "open",
                1,
                "MOVED_FROM_HISTORY",
            ),
        }
        before_status, before_priority, after_status, after_priority, expected_class = (
            transitions[slug]
        )
        current_target = v2_synthetic_target(
            0,
            status=after_status,
            priority=after_priority,
        )
        inventory = v2_synthetic_inventory([current_target])
        current = inventory["rows"][0]
        prior = v2_rooted(
            {
                "state": "PROVIDED",
                "terminal_root": semantic_root({"terminal": "prior"}),
                "rows": [
                    {
                        "id": current["id"],
                        "status": before_status,
                        "priority": before_priority,
                        "target_root": semantic_root({"before": current["id"]}),
                        "destination": (
                            "history-v2.json"
                            if before_status == "closed"
                            else "review-plan-v2.json"
                        ),
                    }
                ],
            }
        )
        zero = v2_build_zero_sets(
            source_root=semantic_root({"fixture": slug}),
            inventory=inventory,
            prior_campaign=prior,
        )
        movement = zero["movements"][0]
        checks.check(
            "movement-class-lineage",
            movement["class"] == expected_class
            and movement["successor_lineage_root"].startswith("sha256-v1:"),
            expected=expected_class,
            observed=movement,
        )
    else:
        raise EvidenceFailed(f"no history executor for {slug}")
    return checks.finish(parameters)


def v2_synthetic_payloads(
    manifest: Mapping[str, Any],
    *,
    oversize: bool = False,
) -> tuple[dict[str, bytes], tuple[Any, ...]]:
    components = v2_synthetic_bundle_components(
        manifest,
        target_count=1 if oversize else 2,
        oversize=oversize,
    )
    source, inventory, authority, review, history, zero, optional = components
    reproduction = [
        str(SCRIPT_REL),
        "--review-plan",
        "--artifact-root",
        "target/v2-fixture",
        "--artifact-dir",
        "run",
    ]
    payloads = v2_bundle_payloads(
        mode="review-plan",
        subject_mode="review-plan",
        artifact_root="target/v2-fixture",
        artifact_dir="run",
        source=source,
        inventory=inventory,
        authority=authority,
        review_plan=review,
        history=history,
        zero_sets=zero,
        optional_payloads=optional,
        reproduction=reproduction,
    )
    return payloads, components


def v2_execute_artifact_cases(
    case: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    slug = str(case["id"]).removeprefix("template-lint-v2.")
    checks = V2CheckCollector(case["id"], case["assertion_id"])
    parameters = {"family": "artifact-log-replay", "slug": slug}
    if slug == "artifact-review-base-set":
        payloads, components = v2_synthetic_payloads(manifest)
        checks.check(
            "review-nine-base-files",
            set(V2_RUN_ARTIFACTS).issubset(payloads)
            and len(set(V2_RUN_ARTIFACTS)) == 9
            and components[3]["state"] == "POPULATED"
            and components[4]["state"] == "NOT_REQUESTED",
            expected={"base": 9, "review": "POPULATED", "history": "NOT_REQUESTED"},
            observed={
                "base": len(set(V2_RUN_ARTIFACTS) & set(payloads)),
                "review": components[3]["state"],
                "history": components[4]["state"],
            },
        )
    elif slug == "artifact-history-base-set":
        inventory, authority, issues, audit, contract = v2_history_fixture(
            manifest
        )
        source = v2_minimal_source(
            manifest_root=manifest["semantic_root"],
            subject="history-bundle",
        )
        review = v2_not_requested(
            schema=V2_REVIEW_PLAN_SCHEMA,
            mode="history-plan",
            source_root=source["semantic_root"],
            inventory_root=inventory["semantic_root"],
            projection="review-plan",
        )
        history = v2_build_history(
            source_root=source["semantic_root"],
            inventory=inventory,
            authority=authority,
            all_issues=issues,
            audit_capture=audit,
            history_contract=contract,
        )
        zero = v2_build_zero_sets(
            source_root=source["semantic_root"],
            inventory=inventory,
            prior_campaign=v2_empty_prior_campaign(),
        )
        payloads = v2_bundle_payloads(
            mode="history-plan",
            subject_mode="history-plan",
            artifact_root="target/v2-fixture",
            artifact_dir="history-run",
            source=source,
            inventory=inventory,
            authority=authority,
            review_plan=review,
            history=history,
            zero_sets=zero,
            optional_payloads={},
            reproduction=["script", "--history-plan"],
        )
        checks.check(
            "history-nine-base-files",
            set(payloads) == set(V2_RUN_ARTIFACTS)
            and history["state"] == "POPULATED"
            and review["state"] == "NOT_REQUESTED",
            expected={"files": 9, "history": "POPULATED", "review": "NOT_REQUESTED"},
            observed={
                "files": len(payloads),
                "history": history["state"],
                "review": review["state"],
            },
        )
    elif slug == "artifact-not-requested-projection":
        document = v2_not_requested(
            schema=V2_HISTORY_SCHEMA,
            mode="review-plan",
            source_root=semantic_root({"source": 1}),
            inventory_root=semantic_root({"inventory": 1}),
            projection="history",
        )
        expected = {
            "schema",
            "state",
            "mode",
            "projection",
            "source_root",
            "inventory_root",
            "rows",
            "no_claim",
            "semantic_root",
        }
        checks.check(
            "not-requested-closed-schema",
            set(document) == expected
            and document["state"] == "NOT_REQUESTED"
            and document["rows"] == [],
            expected=sorted(expected),
            observed=sorted(document),
        )
        mutated = dict(document)
        mutated["state"] = "POPULATED"
        checks.check(
            "empty-populated-root-differs",
            semantic_root(
                {key: value for key, value in mutated.items() if key != "semantic_root"}
            )
            != document["semantic_root"],
            expected="different root",
            observed=semantic_root(
                {key: value for key, value in mutated.items() if key != "semantic_root"}
            ),
        )
    elif slug == "artifact-oversize-inventory":
        payloads, components = v2_synthetic_payloads(
            manifest,
            oversize=True,
        )
        review = components[3]
        entry = review["oversize_content"][0]
        optional_name = entry["relative_path"]
        checks.check(
            "registered-optional-exact",
            optional_name in payloads
            and entry["byte_length"] == len(payloads[optional_name])
            and entry["aggregate_cap_accounting"] == len(payloads[optional_name]),
            expected=entry["byte_length"],
            observed=len(payloads[optional_name]),
        )
        checks.check(
            "optional-roots-complete",
            bool(entry["target_roots"])
            and bool(entry["field_roots"])
            and bool(entry["clause_roots"]),
            expected=True,
            observed=entry,
        )
    elif slug == "artifact-no-overwrite-partial":
        checks.refuses(
            "existing-output-refused",
            EvidenceFailed,
            lambda: require_fresh_run_dir(
                REPO_ROOT / "target",
                label="fixture output",
            ),
            contains="already exists",
        )
        prefix = {"source-v2.json": b"{}\n"}
        checks.check(
            "partial-prefix-not-complete",
            "terminal.json" not in prefix,
            expected="terminal absent",
            observed=sorted(prefix),
        )
    elif slug == "artifact-membership-path-safety":
        for index, value in enumerate(
            ("../escape", "/absolute", "a//b", "a/./b", "a\\b", "e\u0301")
        ):
            checks.refuses(
                f"unsafe-member-{index}",
                InputRefused,
                lambda value=value: v2_safe_member(value, label="fixture"),
                contains="unsafe",
            )
        expected = sorted(V2_RUN_ARTIFACTS)
        extra = sorted([*V2_RUN_ARTIFACTS, "extra.json"])
        checks.check(
            "extra-membership-detected",
            expected != extra,
            expected=expected,
            observed=extra,
        )
        identities: dict[tuple[int, int], str] = {}
        v2_register_member_identity(
            identities,
            device=1,
            inode=100,
            link_count=1,
            relative="source-v2.json",
            kind="file",
        )
        identity_projection = lambda: [
            {
                "device": device,
                "inode": inode,
                "path": path,
            }
            for (device, inode), path in sorted(identities.items())
        ]
        checks.refuses(
            "inode-alias-refused",
            InputRefused,
            lambda: v2_register_member_identity(
                identities,
                device=1,
                inode=100,
                link_count=1,
                relative="inventory-v2.json",
                kind="file",
            ),
            contains="inode alias",
            projection=identity_projection,
            projection_label="BUNDLE_MEMBER_IDENTITIES",
        )
        checks.refuses(
            "external-hardlink-refused",
            InputRefused,
            lambda: v2_register_member_identity(
                identities,
                device=1,
                inode=101,
                link_count=2,
                relative="authority-v2.json",
                kind="file",
            ),
            contains="hard-link alias",
            projection=identity_projection,
            projection_label="BUNDLE_MEMBER_IDENTITIES",
        )
    elif slug == "log-assertion-argv-evidence":
        payloads, _ = v2_synthetic_payloads(manifest)
        events = v2_read_event_stream(payloads["events.jsonl"])
        checks.check(
            "events-have-required-evidence",
            all(
                event["case_id"]
                and event["assertion_id"] == event["executor_id"]
                and isinstance(event["argv"], list)
                and event["parameter_root"].startswith("sha256-v1:")
                for event in events
            ),
            expected=True,
            observed={
                "events": len(events),
                "root": semantic_root(events),
            },
        )
    elif slug == "log-raw-result-redaction":
        hostile = [
            "/Users/example/repo",
            "access_token=top-secret",
            "ordinary",
        ]
        sanitized = v2_sanitized_argv(hostile)
        checks.check(
            "sensitive-argv-redacted",
            sanitized[0].startswith("<redacted:")
            and sanitized[1].startswith("<redacted:")
            and sanitized[2] == "ordinary"
            and "top-secret" not in canonical_bytes(sanitized).decode("utf-8"),
            expected=["redacted", "redacted", "ordinary"],
            observed=sanitized,
        )
        body = b"secret stream body"
        receipt = {
            "bytes": len(body),
            "root": "sha256-v1:" + hashlib.sha256(body).hexdigest(),
            "body_retained": False,
        }
        checks.check(
            "raw-stream-root-only",
            "secret stream body" not in canonical_bytes(receipt).decode("utf-8")
            and receipt["body_retained"] is False,
            expected=False,
            observed=receipt,
        )
    elif slug == "log-order-terminal":
        payloads, _ = v2_synthetic_payloads(manifest)
        events = v2_read_event_stream(payloads["events.jsonl"])
        terminal = strict_json_loads(
            payloads["terminal.json"],
            label="synthetic terminal",
            require_canonical=True,
        )
        checks.check(
            "dense-unique-last-terminal",
            [event["sequence"] for event in events] == list(range(len(events)))
            and sum(event["terminal"] is not None for event in events) == 1
            and events[-1]["terminal"] == "Pass",
            expected=True,
            observed={
                "sequence": [event["sequence"] for event in events],
                "terminals": [
                    event["sequence"]
                    for event in events
                    if event["terminal"] is not None
                ],
            },
        )
        checks.check(
            "terminal-binds-event-roots",
            terminal["event_roots"] == [semantic_root(event) for event in events]
            and terminal["terminal_event_root"] == semantic_root(events[-1]),
            expected=terminal["event_roots"],
            observed=[semantic_root(event) for event in events],
        )
    elif slug in {"replay-review-offline", "replay-history-offline"}:
        before_receipts = len(_command_receipts)
        subject_mode = (
            "review-plan"
            if slug == "replay-review-offline"
            else "history-plan"
        )
        payloads, accepted_manifest, retained = v2_replayable_fixture_bundle(
            manifest,
            subject_mode=subject_mode,
        )
        terminal = strict_json_loads(
            payloads["terminal.json"],
            label="replayable fixture terminal",
            require_canonical=True,
        )
        events = v2_read_event_stream(payloads["events.jsonl"])
        reconstructed = v2_reconstruct_retained(
            artifact_root="target/v2-fixture",
            input_dir=f"{subject_mode}-run",
            payloads=payloads,
            terminal=terminal,
            events=events,
            accepted_manifest=accepted_manifest,
        )
        after_receipts = len(_command_receipts)
        checks.check(
            "artifact-only-reconstruction-exact",
            [row["semantic_root"] for row in reconstructed[:6]]
            == [row["semantic_root"] for row in retained[:6]]
            and reconstructed[6] == retained[6],
            expected=[row["semantic_root"] for row in retained[:6]],
            observed=[row["semantic_root"] for row in reconstructed[:6]],
        )
        checks.check(
            "selected-and-not-requested-projections",
            (
                reconstructed[3]["state"] == "POPULATED"
                and reconstructed[4]["state"] == "NOT_REQUESTED"
                if subject_mode == "review-plan"
                else reconstructed[3]["state"] == "NOT_REQUESTED"
                and reconstructed[4]["state"] == "POPULATED"
            ),
            expected=(
                ("POPULATED", "NOT_REQUESTED")
                if subject_mode == "review-plan"
                else ("NOT_REQUESTED", "POPULATED")
            ),
            observed=(
                reconstructed[3]["state"],
                reconstructed[4]["state"],
            ),
        )
        checks.check(
            "offline-no-command-calls",
            before_receipts == after_receipts,
            expected=before_receipts,
            observed=after_receipts,
        )
    elif slug == "replay-tamper-live-denied":
        payloads, _ = v2_synthetic_payloads(manifest)
        terminal = strict_json_loads(
            payloads["terminal.json"],
            label="synthetic terminal",
            require_canonical=True,
        )
        changed = payloads["inventory-v2.json"] + b" "
        observed_identity = v2_artifact_identity(
            relative_path="inventory-v2.json",
            schema_kind=V2_INVENTORY_SCHEMA,
            payload=changed,
        )
        checks.check(
            "changed-artifact-identity-detected",
            observed_identity
            != terminal["artifact_identities"]["inventory-v2.json"],
            expected=terminal["artifact_identities"]["inventory-v2.json"],
            observed=observed_identity,
        )
        checks.check(
            "mixed-v1-schema-refused",
            RUN_TERMINAL_SCHEMA != V2_TERMINAL_SCHEMA,
            expected=V2_TERMINAL_SCHEMA,
            observed=RUN_TERMINAL_SCHEMA,
        )
    else:
        raise EvidenceFailed(f"no artifact executor for {slug}")
    return checks.finish(parameters)


def v2_failure_evidence(
    *,
    stage: str,
    detail: str,
    expected: Any,
    observed: Any,
) -> dict[str, Any]:
    divergence = first_projection_divergence(expected, observed)
    return v2_rooted(
        {
            "schema": "frankensim.beads-template-hygiene.failure-evidence.v2",
            "stage": stage,
            "terminal": "EvidenceFailed",
            "detail_root": text_root(detail),
            "first_divergence": divergence or "$",
            "recovery": "MANUAL_RECOVERY_REQUIRED",
            "complete_green_prefix": False,
            "no_claim": (
                "failure evidence identifies recovery work and never seals a "
                "partial or ambiguous publication as green"
            ),
        }
    )


def v2_execute_fault_resource_cases(
    case: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    slug = str(case["id"]).removeprefix("template-lint-v2.")
    checks = V2CheckCollector(case["id"], case["assertion_id"])
    parameters = {"family": "fault-resource-mutation", "slug": slug}
    inventory = v2_synthetic_inventory([v2_synthetic_target(0)])

    if slug == "mutation-forged-approval":
        forged = v2_synthetic_receipts(inventory, declared=True)
        forged_row = dict(forged["receipts"][0])
        forged_row["target_root"] = semantic_root({"forged": "stale-target"})
        forged_row = v2_rooted(forged_row)
        forged_document = dict(forged)
        forged_document["receipts"] = [forged_row]
        forged_document = v2_rooted(forged_document)
        checks.refuses(
            "stale-root-forged-approval",
            InputRefused,
            lambda: v2_derive_authority(
                inventory,
                forged_document,
                current_br_version="0.3.0",
                allow_mechanical_fixture=True,
            ),
            contains="exact source roots drifted",
        )
        untrusted = v2_synthetic_authority(
            inventory,
            declared=True,
            manual=True,
            external_verdict="VALID",
            conditional_verdict="VALID",
            gate_verdict="VALID",
            version="0.3.0",
            allow_mechanical_fixture=False,
        )["decisions"][0]
        checks.check(
            "self-issued-receipts-never-mechanical",
            untrusted["readiness"] == "DECLARED_READY"
            and untrusted["remediation_route"] == "MANUAL_BR_REVIEW"
            and untrusted["external_authority"]["verified"] is False
            and untrusted["conditional_write_capability"]["verified"] is False,
            expected=("DECLARED_READY", "MANUAL_BR_REVIEW", False, False),
            observed=(
                untrusted["readiness"],
                untrusted["remediation_route"],
                untrusted["external_authority"]["verified"],
                untrusted["conditional_write_capability"]["verified"],
            ),
        )
        undeclared_receipts = v2_synthetic_receipts(
            inventory,
            manual=True,
        )
        undeclared = v2_derive_authority(
            inventory,
            undeclared_receipts,
            current_br_version="0.2.19",
        )["decisions"][0]
        checks.check(
            "undeclared-manual-route-refused",
            undeclared["remediation_route"] == "ANALYSIS_ONLY"
            and undeclared["implementation_owner"] == "UNRESOLVED"
            and undeclared["terminal_consumer"] == "UNRESOLVED",
            expected=("ANALYSIS_ONLY", "UNRESOLVED", "UNRESOLVED"),
            observed=(
                undeclared["remediation_route"],
                undeclared["implementation_owner"],
                undeclared["terminal_consumer"],
            ),
        )

    elif slug == "mutation-forged-capability-deferred":
        for verdict in ("LYING", "NON_ATOMIC", "VERSION_MISMATCHED"):
            decision = v2_synthetic_authority(
                inventory,
                declared=True,
                external_verdict="VALID",
                conditional_verdict=verdict,
                gate_verdict="VALID",
                version="0.3.0",
                allow_mechanical_fixture=True,
            )["decisions"][0]
            checks.check(
                f"{verdict.lower()}-cannot-upgrade",
                decision["readiness"] != "MECHANICALLY_APPLY_ELIGIBLE"
                and decision["remediation_route"] != "AUTOMATED_CONDITIONAL",
                expected="non-mechanical",
                observed=(decision["readiness"], decision["remediation_route"]),
            )
        deferred_inventory = v2_synthetic_inventory(
            [v2_synthetic_target(0, status="deferred")]
        )
        deferred = v2_synthetic_authority(
            deferred_inventory,
            declared=True,
            manual=True,
            external_verdict="VALID",
            conditional_verdict="VALID",
            gate_verdict="VALID",
            version="0.3.0",
            allow_mechanical_fixture=True,
        )["decisions"][0]
        checks.check(
            "deferred-never-upgrades",
            deferred["remediation_route"] == "ANALYSIS_ONLY"
            and deferred["deferred_apply_prohibition"] is True,
            expected=("ANALYSIS_ONLY", True),
            observed=(
                deferred["remediation_route"],
                deferred["deferred_apply_prohibition"],
            ),
        )

    elif slug == "mutation-missing-payload-edge":
        authority = v2_synthetic_authority(inventory)
        plan, _ = v2_build_review_plan(inventory, authority, max_targets=10)
        child = plan["children"][0]
        missing = dict(child)
        missing.pop("notes")
        missing = v2_rooted(missing)
        checks.refuses(
            "missing-required-payload-field",
            InputRefused,
            lambda: v2_validate_child_payload(missing),
            contains="non-closed schema",
        )
        membership = dict(child)
        membership["target_ids"] = [child["target_ids"][0], child["target_ids"][0]]
        membership["target_roots"] = [
            child["target_roots"][0],
            child["target_roots"][0],
        ]
        membership = v2_rooted(membership)
        checks.refuses(
            "duplicate-target-membership",
            EvidenceFailed,
            lambda: v2_validate_child_payload(membership),
            contains="target membership",
        )
        for label, edges in (
            (
                "missing-edge",
                child["intended_generated_edges"][:1],
            ),
            (
                "reversed-edge",
                [
                    {
                        **child["intended_generated_edges"][0],
                        "from": child["intended_generated_edges"][0]["to"],
                        "to": child["intended_generated_edges"][0]["from"],
                    },
                    child["intended_generated_edges"][1],
                ],
            ),
            (
                "duplicate-edge",
                [
                    child["intended_generated_edges"][0],
                    child["intended_generated_edges"][0],
                ],
            ),
            (
                "undeclared-edge",
                [
                    *child["intended_generated_edges"],
                    {
                        "from": "generated:extra",
                        "to": child["lane_parent"],
                        "type": "blocks",
                    },
                ],
            ),
        ):
            mutated = dict(child)
            mutated["intended_generated_edges"] = edges
            mutated = v2_rooted(mutated)
            checks.refuses(
                label,
                EvidenceFailed,
                lambda mutated=mutated: v2_validate_child_payload(mutated),
                contains="generated-edge",
            )

    elif slug == "mutation-relation-neighbor-drift":
        targets = [v2_synthetic_target(0), v2_synthetic_target(1)]
        before_inventory = v2_synthetic_inventory(targets)
        before_authority = v2_synthetic_authority(before_inventory)
        before_plan, _ = v2_build_review_plan(
            before_inventory,
            before_authority,
            max_targets=10,
        )
        changed_target = dict(targets[0])
        changed_target["dependencies"] = [
            {
                "id": "synthetic-neighbor",
                "type": "blocks",
                "status": "open",
                "priority": 1,
            }
        ]
        changed_target["dependency_neighborhood_root"] = semantic_root(
            {
                "dependencies": [
                    ("blocks", "synthetic-neighbor", "open", 1)
                ],
                "dependents": [],
            }
        )
        changed_target["target_root"] = semantic_root(
            {
                "prior": targets[0]["target_root"],
                "dependency_neighborhood_root": changed_target[
                    "dependency_neighborhood_root"
                ],
            }
        )
        changed_target = v2_rooted(changed_target)
        after_inventory = v2_synthetic_inventory([changed_target, targets[1]])
        after_authority = v2_synthetic_authority(after_inventory)
        after_plan, _ = v2_build_review_plan(
            after_inventory,
            after_authority,
            max_targets=10,
        )
        before_keys = {
            child["target_ids"][0]: child["child_key"]
            for child in before_plan["children"]
        }
        after_keys = {
            child["target_ids"][0]: child["child_key"]
            for child in after_plan["children"]
        }
        checks.check(
            "only-affected-logical-key-moves",
            before_keys[targets[0]["id"]] != after_keys[targets[0]["id"]]
            and before_keys[targets[1]["id"]] == after_keys[targets[1]["id"]],
            expected={"affected": "changed", "unrelated": "stable"},
            observed={"before": before_keys, "after": after_keys},
        )

    elif slug == "mutation-dropped-duplicate-target":
        targets = [v2_synthetic_target(0), v2_synthetic_target(1)]
        exact_inventory = v2_synthetic_inventory(targets)
        exact_authority = v2_synthetic_authority(exact_inventory)
        plan, _ = v2_build_review_plan(
            exact_inventory,
            exact_authority,
            max_targets=10,
        )
        mutations: list[tuple[str, list[dict[str, Any]]]] = [
            ("dropped", plan["children"][:1]),
            ("duplicated", [*plan["children"], plan["children"][0]]),
        ]
        for label, children in mutations:
            mutated = dict(plan)
            mutated["children"] = children
            mutated = v2_rooted(mutated)
            checks.refuses(
                f"{label}-target",
                EvidenceFailed,
                lambda mutated=mutated: v2_validate_review_plan(
                    mutated,
                    exact_inventory,
                    exact_authority,
                ),
                contains="target mapping",
            )
        cross_lane_child = dict(plan["children"][0])
        cross_lane_child["priority"] = 4
        cross_lane_child = v2_rooted(cross_lane_child)
        cross_lane_plan = dict(plan)
        cross_lane_plan["children"] = [
            cross_lane_child,
            *plan["children"][1:],
        ]
        cross_lane_plan = v2_rooted(cross_lane_plan)
        checks.refuses(
            "cross-lane-substitution",
            EvidenceFailed,
            lambda: v2_validate_review_plan(
                cross_lane_plan,
                exact_inventory,
                exact_authority,
            ),
            contains="hard partition",
        )

    elif slug == "fault-output-lock-export":
        checks.refuses(
            "existing-output-reservation",
            EvidenceFailed,
            lambda: require_fresh_run_dir(
                REPO_ROOT,
                label="fixture guaranteed-existing root",
            ),
            contains="already exists",
        )
        checks.refuses(
            "exclusive-writer-no-overwrite",
            EvidenceFailed,
            lambda: v2_write_exclusive(REPO_ROOT / SCRIPT_REL, b"forbidden"),
            contains="overwrite is forbidden",
        )
        failure = v2_failure_evidence(
            stage="output-reservation",
            detail="output parent identity or export state changed",
            expected={"reserved": True, "export": "stable"},
            observed={"reserved": False, "export": "drifted"},
        )
        checks.check(
            "fault-is-non-green",
            failure["terminal"] == "EvidenceFailed"
            and failure["complete_green_prefix"] is False
            and failure["recovery"] == "MANUAL_RECOVERY_REQUIRED",
            expected=("EvidenceFailed", False, "MANUAL_RECOVERY_REQUIRED"),
            observed=(
                failure["terminal"],
                failure["complete_green_prefix"],
                failure["recovery"],
            ),
        )

    elif slug == "cancel-signals-drained":
        global _cancel_requested
        prior_cancel = _cancel_requested
        try:
            for signum in (signal.SIGINT, signal.SIGTERM):
                request_cancel(signum, None)
                for boundary in (
                    "source",
                    "audit",
                    "packing",
                    "publication",
                    "replay",
                ):
                    checks.refuses(
                        f"{signal.Signals(signum).name.lower()}-{boundary}",
                        CancelledDrained,
                        check_cancel,
                        contains="cancellation requested",
                    )
                _cancel_requested = False
        finally:
            _cancel_requested = prior_cancel
        checks.check(
            "cancellation-state-restored",
            _cancel_requested == prior_cancel,
            expected=prior_cancel,
            observed=_cancel_requested,
        )

    elif slug == "restart-deterministic":
        first_plan, first_authority, first_optional = v2_plan_for_count(
            14,
            max_targets=10,
        )
        second_plan, second_authority, second_optional = v2_plan_for_count(
            14,
            max_targets=10,
        )
        checks.check(
            "fresh-restart-same-logical-state",
            first_plan["semantic_root"] == second_plan["semantic_root"]
            and first_authority["semantic_root"]
            == second_authority["semantic_root"]
            and first_optional == second_optional,
            expected=(
                first_plan["semantic_root"],
                first_authority["semantic_root"],
            ),
            observed=(
                second_plan["semantic_root"],
                second_authority["semantic_root"],
            ),
        )
        checks.refuses(
            "restart-still-no-overwrite",
            EvidenceFailed,
            lambda: require_fresh_run_dir(
                REPO_ROOT,
                label="fixture guaranteed-existing root",
            ),
            contains="already exists",
        )

    elif slug == "resource-all-caps":
        caps = {
            "inventory_rows": V2_INVENTORY_ROWS_CAP,
            "warning_rows": V2_WARNING_ROWS_CAP,
            "warnings_per_issue": V2_WARNINGS_PER_ISSUE_CAP,
            "description": V2_CHILD_DESCRIPTION_CAP,
            "acceptance": V2_CHILD_ACCEPTANCE_CAP,
            "design": V2_CHILD_DESIGN_CAP,
            "notes": V2_CHILD_NOTES_CAP,
            "retained_payload": V2_CHILD_PAYLOAD_CAP,
            "argv_count": V2_COMMAND_ARGUMENTS_CAP,
            "argv_bytes": V2_COMMAND_ARGUMENT_BYTES_CAP,
            "relative_path": CAPS["relative_path_bytes"],
            "path_component": CAPS["path_component_bytes"],
            "path_depth": CAPS["path_depth"],
            "artifact": RUN_ARTIFACT_CAP,
            "events": V2_LOG_EVENTS_CAP,
            "log_line": V2_LOG_LINE_BYTES_CAP,
            "synopsis": V2_SYNOPSIS_BYTES_CAP,
            "selected_ids": V2_SYNOPSIS_ID_PREVIEW_CAP,
        }
        checks.check(
            "manifest-cap-bindings",
            all(manifest["caps"].get(name) == cap for name, cap in {
                "max_inventory_rows": V2_INVENTORY_ROWS_CAP,
                "max_warning_rows": V2_WARNING_ROWS_CAP,
                "max_warnings_per_issue": V2_WARNINGS_PER_ISSUE_CAP,
                "max_description_bytes": V2_CHILD_DESCRIPTION_CAP,
                "max_acceptance_bytes": V2_CHILD_ACCEPTANCE_CAP,
                "max_design_bytes": V2_CHILD_DESIGN_CAP,
                "max_notes_bytes": V2_CHILD_NOTES_CAP,
                "max_retained_child_payload_bytes": V2_CHILD_PAYLOAD_CAP,
                "max_command_arguments": V2_COMMAND_ARGUMENTS_CAP,
                "max_command_argument_bytes": V2_COMMAND_ARGUMENT_BYTES_CAP,
                "max_relative_path_bytes": CAPS["relative_path_bytes"],
                "max_path_component_bytes": CAPS["path_component_bytes"],
                "max_path_depth": CAPS["path_depth"],
                "max_artifact_bytes": RUN_ARTIFACT_CAP,
                "max_log_events": V2_LOG_EVENTS_CAP,
                "max_log_line_bytes": V2_LOG_LINE_BYTES_CAP,
                "max_synopsis_bytes": V2_SYNOPSIS_BYTES_CAP,
                "max_synopsis_selected_ids": V2_SYNOPSIS_ID_PREVIEW_CAP,
            }.items()),
            expected="all manifest caps bound",
            observed=manifest["caps"],
        )
        for name, cap in caps.items():
            observed = [value <= cap for value in (cap - 1, cap, cap + 1)]
            checks.check(
                f"{name}-n-minus-one-n-plus-one",
                observed == [True, True, False],
                expected=[True, True, False],
                observed=observed,
            )
        v2_sanitized_argv(["x"] * V2_COMMAND_ARGUMENTS_CAP)
        checks.refuses(
            "argv-count-cap-plus-one",
            EvidenceFailed,
            lambda: v2_sanitized_argv(
                ["x"] * (V2_COMMAND_ARGUMENTS_CAP + 1)
            ),
            contains="argument-count",
        )
        v2_sanitized_argv(["x" * V2_COMMAND_ARGUMENT_BYTES_CAP])
        checks.refuses(
            "argv-byte-cap-plus-one",
            EvidenceFailed,
            lambda: v2_sanitized_argv(
                ["x" * (V2_COMMAND_ARGUMENT_BYTES_CAP + 1)]
            ),
            contains="argument-byte",
        )
        safe_relative(
            "x" * CAPS["path_component_bytes"],
            label="fixture component",
        )
        checks.refuses(
            "path-component-cap-plus-one",
            UsageRefused,
            lambda: safe_relative(
                "x" * (CAPS["path_component_bytes"] + 1),
                label="fixture component",
            ),
            contains="unsafe or over-cap",
        )

    elif slug == "mutation-truncation-oversize-root":
        components = v2_synthetic_bundle_components(
            manifest,
            target_count=1,
            oversize=True,
        )
        source, inv, authority, review, history, zero, optional = components
        entry = review["oversize_content"][0]
        altered_entry = dict(entry)
        altered_entry["byte_length"] = int(entry["byte_length"]) - 1
        altered_review = dict(review)
        altered_review["oversize_content"] = [altered_entry]
        altered_review = v2_rooted(altered_review)
        checks.refuses(
            "oversize-registry-byte-laundering",
            EvidenceFailed,
            lambda: v2_bundle_payloads(
                mode="review-plan",
                subject_mode="review-plan",
                artifact_root="target/v2-fixture",
                artifact_dir="run",
                source=source,
                inventory=inv,
                authority=authority,
                review_plan=altered_review,
                history=history,
                zero_sets=zero,
                optional_payloads=optional,
                reproduction=["script", "--review-plan"],
            ),
            contains="identity differs",
        )
        checks.refuses(
            "oversize-content-omission",
            EvidenceFailed,
            lambda: v2_bundle_payloads(
                mode="review-plan",
                subject_mode="review-plan",
                artifact_root="target/v2-fixture",
                artifact_dir="run",
                source=source,
                inventory=inv,
                authority=authority,
                review_plan=review,
                history=history,
                zero_sets=zero,
                optional_payloads={},
                reproduction=["script", "--review-plan"],
            ),
            contains="membership differs",
        )

    elif slug == "fault-partial-publication-conflict":
        payloads, _ = v2_synthetic_payloads(manifest)
        events = v2_read_event_stream(payloads["events.jsonl"])
        second_terminal = dict(events[-1])
        second_terminal["sequence"] = len(events)
        duplicate_terminal = b"".join(
            canonical_bytes(row) for row in [*events, second_terminal]
        )
        checks.refuses(
            "duplicate-terminal-refused",
            EvidenceFailed,
            lambda: v2_read_event_stream(duplicate_terminal),
            contains="unique and last",
        )
        prefix = {name: value for name, value in payloads.items() if name != "terminal.json"}
        recovery = v2_failure_evidence(
            stage="publication",
            detail="terminal missing after a partial publication",
            expected=sorted(payloads),
            observed=sorted(prefix),
        )
        checks.check(
            "partial-prefix-fail-sealed",
            "terminal.json" not in prefix
            and recovery["complete_green_prefix"] is False
            and recovery["first_divergence"] != "$",
            expected=True,
            observed=recovery,
        )

    elif slug == "fault-replay-manual-recovery":
        expected = {
            "source": "root-a",
            "audit": {"actor": "a", "event": 1},
            "publication": ["source-v2.json", "terminal.json"],
        }
        variants = {
            "replay-mismatch": {
                **expected,
                "source": "root-b",
            },
            "ambiguous-source": {
                **expected,
                "publication": ["source-v2.json", "source-copy.json"],
            },
            "audit-conflict": {
                **expected,
                "audit": {"actor": "b", "event": 1},
            },
            "publication-ambiguity": {
                **expected,
                "publication": [
                    "source-v2.json",
                    "terminal.json",
                    "terminal-copy.json",
                ],
            },
        }
        for label, observed in variants.items():
            evidence = v2_failure_evidence(
                stage=label,
                detail=f"{label} requires operator adjudication",
                expected=expected,
                observed=observed,
            )
            checks.check(
                label,
                evidence["first_divergence"].startswith("$")
                and evidence["recovery"] == "MANUAL_RECOVERY_REQUIRED"
                and evidence["terminal"] == "EvidenceFailed",
                expected=("first divergence", "MANUAL_RECOVERY_REQUIRED"),
                observed=evidence,
            )

    else:
        raise EvidenceFailed(f"no fault/resource executor for {slug}")
    return checks.finish(parameters)


def v2_nomock_bv_projection() -> dict[str, Any]:
    completed = run_command(("bv", "--robot-triage"))
    document = strict_json_loads(
        completed.stdout,
        label="bv robot triage",
    )
    if not isinstance(document, dict):
        raise InfrastructureFailed("bv robot triage is not an object")
    triage = document.get("triage")
    if not isinstance(triage, dict):
        raise InfrastructureFailed("bv robot triage lacks triage data")
    meta = triage.get("meta")
    health = triage.get("project_health")
    quick = triage.get("quick_ref")
    if not all(isinstance(value, dict) for value in (meta, health, quick)):
        raise InfrastructureFailed("bv robot triage stable sections are malformed")
    projection = {
        "data_hash": str(document.get("data_hash") or ""),
        "version": str(meta.get("version") or ""),
        "issue_count": meta.get("issue_count"),
        "phase2_ready": meta.get("phase2_ready"),
        "counts": health.get("counts"),
        "graph": health.get("graph"),
        "quick_ref": {
            "open_count": quick.get("open_count"),
            "actionable_count": quick.get("actionable_count"),
            "blocked_count": quick.get("blocked_count"),
            "in_progress_count": quick.get("in_progress_count"),
        },
        "robot_mode": "--robot-triage",
        "raw_generated_at_retained": False,
    }
    return v2_rooted(projection)


def v2_nomock_history_evidence(
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    global _v2_nomock_history_cache
    if _v2_nomock_history_cache is not None:
        return dict(_v2_nomock_history_cache)
    contract = manifest["history_contract"]
    anchor_id = str(contract["legacy_coverage_anchor_issue"])
    legacy_id = "frankensim-0aeh"
    first = normalize_document(
        br_read_json(
            "list",
            "--status",
            "closed",
            "--json",
            "--limit",
            "0",
        ),
        label="no-mock closed list first",
    )
    second = normalize_document(
        br_read_json(
            "list",
            "--status",
            "closed",
            "--json",
            "--limit",
            "0",
        ),
        label="no-mock closed list second",
    )
    first_projection = [
        {
            "id": str(row.get("id") or ""),
            "closed_at": str(row.get("closed_at") or ""),
            "close_reason": str(row.get("close_reason") or ""),
        }
        for row in first
    ]
    second_projection = [
        {
            "id": str(row.get("id") or ""),
            "closed_at": str(row.get("closed_at") or ""),
            "close_reason": str(row.get("close_reason") or ""),
        }
        for row in second
    ]
    first_projection.sort(key=lambda row: row["id"])
    second_projection.sort(key=lambda row: row["id"])
    if first_projection != second_projection:
        raise InputRefused("ConcurrentDrift: closed membership changed during audit")
    show_document = br_read_json("show", anchor_id, legacy_id, "--json")
    if not isinstance(show_document, list) or len(show_document) != 2:
        raise InfrastructureFailed("no-mock history show did not return two rows")
    show_by_id = {
        str(row.get("id") or ""): v2_full_issue_projection(row)
        for row in show_document
        if isinstance(row, dict)
    }
    if set(show_by_id) != {anchor_id, legacy_id}:
        raise InfrastructureFailed("no-mock history show membership differs")
    anchor_first = v2_capture_one_audit(anchor_id, capture_ordinal=1)[1]
    anchor_second = v2_capture_one_audit(anchor_id, capture_ordinal=2)[1]
    legacy_first = v2_capture_one_audit(legacy_id, capture_ordinal=1)[1]
    legacy_second = v2_capture_one_audit(legacy_id, capture_ordinal=2)[1]
    if anchor_first != anchor_second or legacy_first != legacy_second:
        raise InputRefused("ConcurrentDrift: selected audit evidence changed")
    anchor_timestamp = v2_parse_timestamp(
        contract["legacy_coverage_anchor_closed_at"],
        label="no-mock history anchor",
    )
    legacy_rows = sorted(
        [
            {
                "id": row["id"],
                "closed_at": row["closed_at"],
                "close_reason_root": text_root(row["close_reason"]),
            }
            for row in first_projection
            if row["closed_at"]
            and v2_parse_timestamp(
                row["closed_at"],
                label=f"no-mock closed {row['id']}",
            )
            < anchor_timestamp
        ],
        key=lambda row: row["id"],
    )
    status_id = int(contract["legacy_coverage_anchor_status_event_id"])
    close_id = int(contract["legacy_coverage_anchor_close_event_id"])
    anchor_events = {
        event["id"]: event for event in anchor_first["events"]
    }
    status_event = anchor_events.get(status_id)
    close_event = anchor_events.get(close_id)
    if (
        status_event is None
        or close_event is None
        or status_event.get("new_value") != "closed"
        or close_event.get("comment") != show_by_id[anchor_id]["close_reason"]
        or status_event.get("actor") != close_event.get("actor")
    ):
        raise EvidenceFailed("no-mock known closer pair differs")
    if legacy_first["events"]:
        raise EvidenceFailed("no-mock legacy issue unexpectedly has audit events")
    evidence = v2_rooted(
        {
            "anchor_id": anchor_id,
            "anchor_target_root": v2_target_root(show_by_id[anchor_id]),
            "anchor_audit_root": anchor_first["semantic_root"],
            "anchor_actor": close_event.get("actor"),
            "known_pair": [status_id, close_id],
            "legacy_id": legacy_id,
            "legacy_target_root": v2_target_root(show_by_id[legacy_id]),
            "legacy_audit_root": legacy_first["semantic_root"],
            "legacy_count": len(legacy_rows),
            "legacy_rows_root": semantic_root(legacy_rows),
            "closed_membership_count": len(first_projection),
            "closed_projection_root": semantic_root(first_projection),
            "capture_count": 2,
            "raw_stream_bodies_retained": False,
            "no_claim": (
                "released-br audit evidence proves only exact known and "
                "legacy-unavailable closure provenance"
            ),
        }
    )
    if (
        evidence["legacy_count"] != contract["legacy_coverage_count"]
        or evidence["legacy_rows_root"] != contract["legacy_coverage_rows_root"]
    ):
        raise EvidenceFailed("no-mock legacy coverage receipt differs")
    _v2_nomock_history_cache = evidence
    return dict(evidence)


def v2_reconstruct_payload_map(
    payloads: Mapping[str, bytes],
    accepted_manifest: Mapping[str, Any],
) -> tuple[
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, bytes],
]:
    terminal = strict_json_loads(
        payloads["terminal.json"],
        label="no-mock retained terminal",
        require_canonical=True,
    )
    events = v2_read_event_stream(payloads["events.jsonl"])
    return v2_reconstruct_retained(
        artifact_root=str(terminal["artifact_root"]),
        input_dir=str(terminal["artifact_dir"]),
        payloads=payloads,
        terminal=terminal,
        events=events,
        accepted_manifest=accepted_manifest,
    )


def v2_fixture_fingerprint() -> dict[str, Any]:
    context = fixture_br_context()
    issue_ids = sorted(
        str(row["id"]) for row in fixture_br_list()
    )
    shows = [
        fixture_br_semantic_projection(fixture_br_show(issue_id))
        for issue_id in issue_ids
    ]
    audits = [
        {
            "issue_id": issue_id,
            "root": semantic_root(
                fixture_br_json("audit", "log", issue_id)
            ),
        }
        for issue_id in issue_ids
    ]
    status = fixture_br_json("status", "--no-activity")
    sync_status = fixture_br_json("sync", "--status")
    witness = fixture_br_json("sync", "--witness")
    return v2_rooted(
        {
            "fixture_db": context["db"],
            "issue_ids": issue_ids,
            "shows": shows,
            "audits": audits,
            "status": status,
            "sync_status": sync_status,
            "export_witness": witness,
            "no_claim": (
                "fingerprint covers the isolated released-br fixture and "
                "does not claim wall-clock or append-only audit reversibility"
            ),
        }
    )


def v2_execute_nomock_cases(
    case: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    slug = str(case["id"]).removeprefix("template-lint-v2.")
    checks = V2CheckCollector(case["id"], case["assertion_id"])
    parameters = {"family": "released-br-no-mock-e2e", "slug": slug}

    if slug == "nomock-br-robot-cells":
        issues = normalize_document(
            br_read_json(
                "list",
                "--all",
                "--deferred",
                "--json",
                "--limit",
                "0",
            ),
            label="no-mock all-status list",
        )
        bv_projection = v2_nomock_bv_projection()
        statuses = Counter(str(row.get("status") or "") for row in issues)
        priorities = Counter(int(row.get("priority")) for row in issues)
        cells = {
            f"P{priority}/{status}": sum(
                row.get("priority") == priority
                and row.get("status") == status
                for row in issues
            )
            for priority in range(5)
            for status in STATUS_SCOPES
        }
        zero_receipts = [
            v2_rooted(
                {
                    "cell": cell,
                    "count": count,
                    "issue_ids_root": semantic_root(
                        sorted(
                            str(row.get("id") or "")
                            for row in issues
                            if f"P{row.get('priority')}/{row.get('status')}" == cell
                        )
                    ),
                }
            )
            for cell, count in sorted(cells.items())
            if count == 0
        ]
        checks.check(
            "released-br-all-state-coverage",
            set(statuses) == set(STATUS_SCOPES)
            and set(priorities) == set(range(5))
            and any(str(row.get("assignee") or "") for row in issues)
            and any(not str(row.get("assignee") or "") for row in issues),
            expected={
                "statuses": list(STATUS_SCOPES),
                "priorities": list(range(5)),
                "assigned_and_ownerless": True,
            },
            observed={
                "statuses": dict(statuses),
                "priorities": dict(priorities),
                "assigned": sum(bool(row.get("assignee")) for row in issues),
            },
        )
        checks.check(
            "robot-deferred-accounting",
            bv_projection["issue_count"]
            == len(issues) - statuses["deferred"]
            and bv_projection["robot_mode"] == "--robot-triage",
            expected=len(issues) - statuses["deferred"],
            observed=bv_projection["issue_count"],
        )
        checks.check(
            "exact-25-cells-and-zero-receipts",
            len(cells) == 25
            and all(
                receipt["count"] == 0
                and receipt["issue_ids_root"] == semantic_root([])
                for receipt in zero_receipts
            ),
            expected={"cells": 25, "zero_receipts_root_bound": True},
            observed={
                "cells": len(cells),
                "zero_receipts": len(zero_receipts),
            },
        )

    elif slug == "nomock-audit-known-legacy-conflict":
        evidence = v2_nomock_history_evidence(manifest)
        checks.check(
            "released-known-and-legacy",
            evidence["known_pair"]
            == [
                manifest["history_contract"][
                    "legacy_coverage_anchor_status_event_id"
                ],
                manifest["history_contract"][
                    "legacy_coverage_anchor_close_event_id"
                ],
            ]
            and evidence["legacy_count"]
            == manifest["history_contract"]["legacy_coverage_count"],
            expected="pinned known pair and legacy receipt",
            observed=evidence,
        )
        conflicting = {
            "issue_id": "fixture-conflict",
            "events": [
                {
                    "id": 3,
                    "event_type": "closed",
                    "actor": "other",
                    "timestamp": "2026-01-02T00:00:00.002000+00:00",
                },
                {
                    "id": 2,
                    "event_type": "closed",
                    "actor": "first",
                    "timestamp": "2026-01-02T00:00:00.001000+00:00",
                },
                {
                    "id": 1,
                    "event_type": "status_changed",
                    "actor": "first",
                    "timestamp": "2026-01-02T00:00:00+00:00",
                    "new_value": "closed",
                },
            ],
        }
        normalized = v2_normalize_audit_document(
            conflicting,
            issue_id="fixture-conflict",
        )
        checks.check(
            "conflicting-close-state-retained",
            sum(
                event["event_type"] == "closed"
                for event in normalized["events"]
            )
            == 2,
            expected=2,
            observed=normalized["events"],
        )

    elif slug == "nomock-current-br-two-gate-nodata":
        context = fixture_br_context()
        before = fixture_br_semantic_projection(
            fixture_br_show(context["stable_id"])
        )
        payloads, accepted, components = v2_replayable_fixture_bundle(
            manifest,
            subject_mode="review-plan",
            issue=fixture_br_show(context["stable_id"]),
        )
        inventory = components[1]
        receipts = v2_synthetic_receipts(
            inventory,
            declared=True,
            external_verdict="VALID",
            conditional_verdict="NODATA",
        )
        version = br_version()["version"]
        authority = v2_derive_authority(
            inventory,
            receipts,
            current_br_version=version,
        )
        after = fixture_br_semantic_projection(
            fixture_br_show(context["stable_id"])
        )
        decision = authority["decisions"][0]
        checks.check(
            "current-released-two-gate-nodata",
            version == "0.2.19"
            and decision["readiness"] == "DECLARED_READY"
            and authority["current_br_conditional_write_capability"]
            == "AUTOMATION_NODATA"
            and decision["remediation_route"] == "ANALYSIS_ONLY",
            expected=(
                "0.2.19",
                "DECLARED_READY",
                "AUTOMATION_NODATA",
                "ANALYSIS_ONLY",
            ),
            observed=(
                version,
                decision["readiness"],
                authority["current_br_conditional_write_capability"],
                decision["remediation_route"],
            ),
        )
        checks.check(
            "released-target-unchanged",
            before == after and bool(payloads) and bool(accepted),
            expected=semantic_root(before),
            observed=semantic_root(after),
        )

    elif slug in {"nomock-review-plan", "nomock-history-plan"}:
        subject_mode = (
            "review-plan"
            if slug == "nomock-review-plan"
            else "history-plan"
        )
        if subject_mode == "history-plan":
            live_history = v2_nomock_history_evidence(manifest)
        else:
            live_history = None
        context = fixture_br_context()
        authentic_issue = (
            fixture_br_show(context["stable_id"])
            if subject_mode == "review-plan"
            else None
        )
        payloads, accepted, retained = v2_replayable_fixture_bundle(
            manifest,
            subject_mode=subject_mode,
            issue=authentic_issue,
        )
        reconstructed = v2_reconstruct_payload_map(payloads, accepted)
        terminal = strict_json_loads(
            payloads["terminal.json"],
            label="no-mock plan terminal",
            require_canonical=True,
        )
        checks.check(
            "exact-nine-file-released-bundle",
            set(payloads) == set(V2_RUN_ARTIFACTS)
            and terminal["terminal"] == "Pass",
            expected=sorted(V2_RUN_ARTIFACTS),
            observed=sorted(payloads),
        )
        checks.check(
            "released-plan-offline-replay",
            [row["semantic_root"] for row in reconstructed[:6]]
            == [row["semantic_root"] for row in retained[:6]],
            expected=[row["semantic_root"] for row in retained[:6]],
            observed=[row["semantic_root"] for row in reconstructed[:6]],
        )
        checks.check(
            "mode-specific-projection-and-audit",
            (
                reconstructed[3]["state"] == "POPULATED"
                and reconstructed[4]["state"] == "NOT_REQUESTED"
                if subject_mode == "review-plan"
                else reconstructed[3]["state"] == "NOT_REQUESTED"
                and reconstructed[4]["state"] == "POPULATED"
                and live_history is not None
                and live_history["legacy_count"] > 0
            ),
            expected=subject_mode,
            observed=(
                reconstructed[3]["state"],
                reconstructed[4]["state"],
            ),
        )

    elif slug == "nomock-unrelated-drift-stability":
        context = fixture_br_context()
        stable_before = fixture_br_semantic_projection(
            fixture_br_show(context["stable_id"])
        )
        regression_before = fixture_br_show(context["regression_id"])
        old_description = str(regression_before.get("description") or "")
        changed_projection: dict[str, Any] | None = None
        try:
            fixture_br_set_description(
                context["regression_id"],
                old_description + "\n\nUnrelated coherent drift fixture.",
            )
            changed_projection = fixture_br_semantic_projection(
                fixture_br_show(context["regression_id"])
            )
            stable_during = fixture_br_semantic_projection(
                fixture_br_show(context["stable_id"])
            )
        finally:
            fixture_br_set_description(
                context["regression_id"],
                old_description,
            )
        stable_after = fixture_br_semantic_projection(
            fixture_br_show(context["stable_id"])
        )
        regression_after = fixture_br_semantic_projection(
            fixture_br_show(context["regression_id"])
        )
        checks.check(
            "unrelated-observation-moves-selected-key-stable",
            changed_projection is not None
            and changed_projection != fixture_br_semantic_projection(
                regression_before
            )
            and stable_before == stable_during == stable_after,
            expected="unrelated changed; selected stable",
            observed={
                "selected_root": semantic_root(stable_after),
                "unrelated_changed_root": semantic_root(changed_projection or {}),
            },
        )
        checks.check(
            "unrelated-fixture-restored",
            regression_after
            == fixture_br_semantic_projection(regression_before),
            expected=semantic_root(
                fixture_br_semantic_projection(regression_before)
            ),
            observed=semantic_root(regression_after),
        )

    elif slug == "nomock-selected-drift-refusal":
        result = fixture_br_description_round_trip(
            case["id"],
            assert_stale_guard=True,
        )
        checks.check(
            "released-selected-stale-root-refused-restored",
            result["semantic_restoration"]
            == "exact-excluding-append-only-audit-history"
            and result["checks"] >= 5,
            expected="exact restoration with stale guard",
            observed=result,
        )
        target = v2_synthetic_target(0)
        roots = {
            "field": target["target_root"],
            "dependency": target["dependency_neighborhood_root"],
            "owner": semantic_root(
                {
                    "assignee": target["tracker_assignee"],
                    "owner": target["tracker_owner"],
                }
            ),
            "receipt": semantic_root({"receipt": "before"}),
        }
        changed = {
            "field": semantic_root({"field": "changed"}),
            "dependency": semantic_root({"dependency": "changed"}),
            "owner": semantic_root({"owner": "changed"}),
            "receipt": semantic_root({"receipt": "after"}),
        }
        checks.check(
            "all-selected-logical-dimensions-invalidate",
            all(roots[key] != changed[key] for key in roots),
            expected="all affected roots changed",
            observed={"before": roots, "after": changed},
        )

    elif slug == "nomock-active-conflict-deferred":
        context = fixture_br_context()
        issue_id = context["stable_id"]
        fixture_br_set_status(issue_id, "open")
        before = fixture_br_semantic_projection(fixture_br_show(issue_id))
        try:
            fixture_br_set_status(issue_id, "in_progress")
            active = fixture_br_show(issue_id)
            payloads, _, components = v2_replayable_fixture_bundle(
                manifest,
                subject_mode="review-plan",
                issue=active,
            )
            decision = components[2]["decisions"][0]
            child = components[3]["children"][0]
            during = fixture_br_semantic_projection(fixture_br_show(issue_id))
        finally:
            fixture_br_set_status(issue_id, "open")
        after = fixture_br_semantic_projection(fixture_br_show(issue_id))
        checks.check(
            "active-conflict-defers-generated-child",
            decision["active_work_context"]["conflict"] is True
            and decision["remediation_route"] == "ANALYSIS_ONLY"
            and child["desired_status"] == "deferred"
            and bool(payloads),
            expected=(True, "ANALYSIS_ONLY", "deferred"),
            observed=(
                decision["active_work_context"]["conflict"],
                decision["remediation_route"],
                child["desired_status"],
            ),
        )
        checks.check(
            "planner-does-not-mutate-active-target",
            during["status"] == "in_progress"
            and after == before,
            expected=("in_progress", semantic_root(before)),
            observed=(during["status"], semantic_root(after)),
        )

    elif slug == "nomock-oversize":
        context = fixture_br_context()
        issue_id = context["stable_id"]
        fixture_br_set_status(issue_id, "open")
        before_issue = fixture_br_show(issue_id)
        before = fixture_br_semantic_projection(before_issue)
        old_description = str(before_issue.get("description") or "")
        try:
            fixture_br_set_description(
                issue_id,
                "x" * (V2_CHILD_DESCRIPTION_CAP + 1),
            )
            changed = fixture_br_show(issue_id)
            payloads, accepted, components = v2_replayable_fixture_bundle(
                manifest,
                subject_mode="review-plan",
                issue=changed,
            )
            review = components[3]
            reconstructed = v2_reconstruct_payload_map(payloads, accepted)
        finally:
            fixture_br_set_description(issue_id, old_description)
        after = fixture_br_semantic_projection(fixture_br_show(issue_id))
        checks.check(
            "released-oversize-complete-registry",
            len(review["oversize_content"]) == 1
            and review["children"][0]["disposition_workflow"]
            == "OVERSIZE_REVIEW_REQUIRED"
            and bool(components[6])
            and reconstructed[3]["semantic_root"] == review["semantic_root"],
            expected="one complete oversize registry and replay",
            observed={
                "registry": review["oversize_content"],
                "child": review["children"][0]["disposition_workflow"],
            },
        )
        checks.check(
            "released-oversize-target-restored",
            after == before,
            expected=semantic_root(before),
            observed=semantic_root(after),
        )

    elif slug == "nomock-offline-replay":
        before_receipts = len(_command_receipts)
        roots: dict[str, list[str]] = {}
        for subject_mode in ("review-plan", "history-plan"):
            payloads, accepted, retained = v2_replayable_fixture_bundle(
                manifest,
                subject_mode=subject_mode,
            )
            reconstructed = v2_reconstruct_payload_map(payloads, accepted)
            roots[subject_mode] = [
                row["semantic_root"] for row in reconstructed[:6]
            ]
            checks.check(
                f"{subject_mode}-offline-exact",
                roots[subject_mode]
                == [row["semantic_root"] for row in retained[:6]],
                expected=[row["semantic_root"] for row in retained[:6]],
                observed=roots[subject_mode],
            )
        checks.check(
            "offline-denies-live-reads",
            len(_command_receipts) == before_receipts,
            expected=before_receipts,
            observed=len(_command_receipts),
        )

    elif slug == "nomock-live-target-unchanged":
        before = v2_fixture_fingerprint()
        v2_replayable_fixture_bundle(
            manifest,
            subject_mode="review-plan",
        )
        v2_replayable_fixture_bundle(
            manifest,
            subject_mode="history-plan",
        )
        after = v2_fixture_fingerprint()
        checks.check(
            "read-only-modes-preserve-full-fixture-fingerprint",
            before["semantic_root"] == after["semantic_root"],
            expected=before["semantic_root"],
            observed=after["semantic_root"],
        )

    elif slug == "nomock-v1-regression":
        environment = os.environ.copy()
        environment["FS_TEMPLATE_HYGIENE_FORBID_V2_LOAD"] = "1"
        completed = subprocess.run(
            [str(REPO_ROOT / SCRIPT_REL), "--self-test"],
            cwd=REPO_ROOT,
            text=False,
            capture_output=True,
            env=environment,
            timeout=CAPS["subprocess_timeout_seconds"],
            check=False,
        )
        if (
            len(completed.stdout) > CAPS["subprocess_stdout_bytes"]
            or len(completed.stderr) > CAPS["subprocess_stdout_bytes"]
        ):
            raise InfrastructureFailed("frozen v1 regression output exceeds cap")
        rows = [
            strict_json_loads(
                line + b"\n",
                label=f"v1 regression line {index}",
                require_canonical=True,
            )
            for index, line in enumerate(completed.stdout.splitlines())
            if line
        ]
        summary = rows[-1] if rows else {}
        checks.check(
            "frozen-v1-suite-exact",
            completed.returncode == 0
            and summary.get("terminal") == "Pass"
            and summary.get("cases") == 21
            and summary.get("total_checks") == 148
            and summary.get("case_manifest_root")
            == manifest["compatibility_contract"][
                "v1_case_manifest_semantic_root"
            ],
            expected={
                "exit": 0,
                "cases": 21,
                "checks": 148,
                "root": manifest["compatibility_contract"][
                    "v1_case_manifest_semantic_root"
                ],
            },
            observed={
                "exit": completed.returncode,
                "cases": summary.get("cases"),
                "checks": summary.get("total_checks"),
                "root": summary.get("case_manifest_root"),
                "stderr_root": (
                    "sha256-v1:"
                    + hashlib.sha256(completed.stderr).hexdigest()
                ),
            },
        )
        checks.check(
            "v1-replay-and-v2-independence",
            "template-lint.artifact-replay"
            in summary.get("artifact_replay_cases", [])
            and "v2 manifest loading was forbidden"
            not in completed.stderr.decode("utf-8", errors="replace"),
            expected=True,
            observed={
                "artifact_replay_cases": summary.get(
                    "artifact_replay_cases", []
                ),
                "stderr_root": (
                    "sha256-v1:"
                    + hashlib.sha256(completed.stderr).hexdigest()
                ),
            },
        )

    else:
        raise EvidenceFailed(f"no no-mock executor for {slug}")
    return checks.finish(parameters)


def v2_build_executor_registry(
    manifest: Mapping[str, Any],
) -> dict[
    str,
    Callable[[Mapping[str, Any], Mapping[str, Any]], dict[str, Any]],
]:
    families = (
        (1, 12, v2_execute_schema_cli_ux),
        (13, 24, v2_execute_source_cases),
        (25, 36, v2_execute_authority_cases),
        (37, 48, v2_execute_packing_cases),
        (49, 60, v2_execute_history_cases),
        (61, 72, v2_execute_artifact_cases),
        (73, 84, v2_execute_fault_resource_cases),
        (85, 96, v2_execute_nomock_cases),
    )
    registry: dict[
        str,
        Callable[[Mapping[str, Any], Mapping[str, Any]], dict[str, Any]],
    ] = {}
    for case in manifest["case"]:
        ordinal = int(case["ordinal"])
        executor = next(
            (
                family
                for first, last, family in families
                if first <= ordinal <= last
            ),
            None,
        )
        if executor is None:
            raise InputRefused(
                f"{case['assertion_id']} has no compiled executor family"
            )
        assertion_id = str(case["assertion_id"])
        if assertion_id in registry:
            raise InputRefused(
                f"duplicate compiled executor {assertion_id}"
            )
        registry[assertion_id] = executor
    return registry


def v2_test_log_event(
    *,
    case: Mapping[str, Any],
    sequence: int,
    stage: str,
    result: Mapping[str, Any] | None = None,
    terminal: str | None = None,
) -> dict[str, Any]:
    assertion_id = str(case["assertion_id"])
    parameter = {
        "case_id": case["id"],
        "assertion_id": assertion_id,
        "ordinal": case["ordinal"],
        "fixture_engine": case["fixture_engine"],
        "stage": stage,
    }
    return v2_rooted(
        {
            "schema": V2_EVENT_SCHEMA,
            "case_id": case["id"],
            "assertion_id": assertion_id,
            "executor_id": assertion_id,
            "parameter_root": semantic_root(parameter),
            "stage": stage,
            "sequence": sequence,
            "argv": v2_sanitized_argv(
                [str(SCRIPT_REL), "--case-v2", str(case["id"])]
            ),
            "exit_code": (
                TERMINAL_EXIT.get(terminal, 0) if terminal else 0
            ),
            "result_category": (
                "PASS"
                if terminal in {None, "Pass"}
                else "EXPECTED_ASSERTION_FAILURE"
            ),
            "semantic_projection": (
                {
                    "result_root": result["semantic_root"],
                    "check_count": result["check_count"],
                    "ordered_check_roots": result["ordered_check_roots"],
                }
                if result is not None
                else {
                    "criterion_ids": list(case["criterion_ids"]),
                    "kinds": list(case["kinds"]),
                    "fixture_engine": case["fixture_engine"],
                }
            ),
            "stdout_byte_length": 0,
            "stdout_root": "sha256-v1:" + hashlib.sha256(b"").hexdigest(),
            "stderr_byte_length": 0,
            "stderr_root": "sha256-v1:" + hashlib.sha256(b"").hexdigest(),
            "first_divergence": None,
            "recovery": "NOT_REQUIRED",
            "terminal": terminal,
            "safe_relative_artifacts": {
                "base": [],
                "optional_count": 0,
                "optional_paths_root": semantic_root([]),
            },
            "no_claim": (
                "the test event records one compiled bounded assertion and "
                "mints no tracker, semantic, implementation, or release authority"
            ),
        }
    )


def run_v2_self_tests(
    manifest: Mapping[str, Any],
    *,
    selected: str | None = None,
) -> dict[str, Any]:
    registry = v2_build_executor_registry(manifest)
    cases = [
        row
        for row in manifest["case"]
        if selected is None or row["id"] == selected
    ]
    if selected is not None and not cases:
        raise UsageRefused(f"unknown v2 case: {selected}")
    results: list[dict[str, Any]] = []
    event_roots: list[str] = []
    sequence = 0
    for case in cases:
        check_cancel()
        start = v2_test_log_event(
            case=case,
            sequence=sequence,
            stage="assertion-start",
        )
        json_stdout(start)
        event_roots.append(start["semantic_root"])
        sequence += 1
        executor = registry.get(str(case["assertion_id"]))
        if executor is None:
            failure = EvidenceFailed(
                f"{case['assertion_id']} has no exact compiled executor"
            )
            terminal = failure.terminal
            failed = v2_test_log_event(
                case=case,
                sequence=sequence,
                stage="assertion-terminal",
                terminal=terminal,
            )
            json_stdout(failed)
            event_roots.append(failed["semantic_root"])
            suite_failure = v2_rooted(
                {
                    "schema": V2_TEST_TERMINAL_SCHEMA,
                    "terminal": terminal,
                    "exit_code": TERMINAL_EXIT[terminal],
                    "selected": selected,
                    "manifest_root": manifest["semantic_root"],
                    "completed_case_roots": [
                        row["semantic_root"] for row in results
                    ],
                    "ordered_event_roots": event_roots,
                    "first_divergence": str(case["assertion_id"]),
                    "recovery": "IMPLEMENT_COMPILED_EXECUTOR",
                    "no_claim": "the v2 suite did not pass",
                }
            )
            json_stdout(suite_failure)
            raise SystemExit(TERMINAL_EXIT[terminal])
        try:
            result = executor(case, manifest)
        except HarnessError as error:
            failed = v2_test_log_event(
                case=case,
                sequence=sequence,
                stage="assertion-terminal",
                terminal=error.terminal,
            )
            failed["semantic_projection"] = {
                "diagnostic_root": text_root(str(error)),
                "subject": case["subject"],
            }
            failed = v2_rooted(failed)
            json_stdout(failed)
            event_roots.append(failed["semantic_root"])
            suite_failure = v2_rooted(
                {
                    "schema": V2_TEST_TERMINAL_SCHEMA,
                    "terminal": error.terminal,
                    "exit_code": TERMINAL_EXIT[error.terminal],
                    "selected": selected,
                    "manifest_root": manifest["semantic_root"],
                    "completed_case_roots": [
                        row["semantic_root"] for row in results
                    ],
                    "ordered_event_roots": event_roots,
                    "first_divergence": str(case["assertion_id"]),
                    "diagnostic_root": text_root(str(error)),
                    "recovery": "FIX_ASSERTION_OR_SUBJECT",
                    "no_claim": "the v2 suite did not pass",
                }
            )
            json_stdout(suite_failure)
            raise SystemExit(TERMINAL_EXIT[error.terminal])
        if (
            result.get("schema") != V2_ASSERTION_RESULT_SCHEMA
            or result.get("case_id") != case["id"]
            or result.get("assertion_id") != case["assertion_id"]
            or result.get("executor_id") != case["assertion_id"]
            or result.get("terminal") != "Pass"
            or not result.get("ordered_checks")
            or result.get("check_count") != len(result["ordered_checks"])
        ):
            raise EvidenceFailed(
                f"{case['assertion_id']} returned malformed assertion evidence"
            )
        json_stdout(result)
        terminal_event = v2_test_log_event(
            case=case,
            sequence=sequence,
            stage="assertion-terminal",
            result=result,
            terminal="Pass",
        )
        json_stdout(terminal_event)
        event_roots.append(terminal_event["semantic_root"])
        sequence += 1
        results.append(result)
    terminal = v2_rooted(
        {
            "schema": V2_TEST_TERMINAL_SCHEMA,
            "terminal": "Pass",
            "exit_code": 0,
            "selected": selected,
            "manifest_root": manifest["semantic_root"],
            "manifest_content_identity": manifest["content_identity"],
            "case_count": len(results),
            "assertion_count": len(results),
            "check_count": sum(int(row["check_count"]) for row in results),
            "ordered_case_ids": [row["case_id"] for row in results],
            "ordered_assertion_ids": [
                row["assertion_id"] for row in results
            ],
            "ordered_case_roots": [
                row["semantic_root"] for row in results
            ],
            "ordered_event_roots": event_roots,
            "compiled_executor_ids_root": semantic_root(sorted(registry)),
            "skipped_assertions": [],
            "aggregate_count_only": False,
            "first_divergence": None,
            "recovery": "NOT_REQUIRED",
            "no_claim": (
                "the suite proves only its 96 bounded executable contracts; "
                "it does not mutate or complete live target Beads"
            ),
        }
    )
    json_stdout(terminal)
    return terminal


class HarnessArgumentParser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        raise UsageRefused(message)


def preflight_arguments(arguments: Sequence[str]) -> set[str]:
    seen: set[str] = set()
    for value in arguments:
        if not value.startswith("--"):
            continue
        option = value.split("=", 1)[0]
        if option == "--output":
            raise UsageRefused("--output is unsupported; use --artifact-root/--artifact-dir")
        if option in seen:
            raise UsageRefused(f"duplicate option is forbidden: {option}")
        seen.add(option)
    return seen


def parse_arguments(arguments: Sequence[str]) -> argparse.Namespace:
    provided = preflight_arguments(arguments)
    parser = HarnessArgumentParser(
        prog=str(SCRIPT_REL),
        description=(
            "Freeze, classify, plan, replay, and explicitly apply exact Beads "
            "template-hygiene work without inventing semantic criteria."
        ),
        allow_abbrev=False,
    )
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--list", action="store_true", help="list frozen case IDs")
    modes.add_argument("--check", action="store_true", help="validate live inventory")
    modes.add_argument("--self-test", action="store_true", help="run all 21 cases")
    modes.add_argument("--inventory", action="store_true", help="write live inventory")
    modes.add_argument("--plan", action="store_true", help="write live review plan")
    modes.add_argument("--apply-manifest", metavar="REL")
    modes.add_argument("--negative", metavar="CASE")
    modes.add_argument("--replay", metavar="REL")
    modes.add_argument("--closeout", action="store_true", help="require zero lint debt")
    modes.add_argument(
        "--review-plan",
        action="store_true",
        help="emit the normalized bounded v2 review plan",
    )
    modes.add_argument(
        "--history-plan",
        action="store_true",
        help="emit immutable v2 closed-history accounting",
    )
    modes.add_argument(
        "--self-test-v2",
        action="store_true",
        help="run all 96 executable v2 assertions",
    )
    modes.add_argument("--case-v2", metavar="CASE")
    parser.add_argument(
        "--artifact-root",
        default="target/beads-template-hygiene",
        help="repository-relative artifact root",
    )
    parser.add_argument(
        "--artifact-dir",
        help="repository-relative directory below artifact root for write modes",
    )
    parser.add_argument(
        "--review-receipts",
        metavar="REL",
        help="safe relative closed review-receipt JSON",
    )
    parser.add_argument(
        "--prior-campaign",
        metavar="REL",
        help="safe relative retained prior v2 campaign bundle",
    )
    parser.add_argument(
        "--priorities",
        metavar="P0,P1,...",
        help="exact comma-separated v2 priority filter",
    )
    parser.add_argument(
        "--statuses",
        metavar="open,closed,...",
        help="exact comma-separated v2 status filter",
    )
    parser.add_argument(
        "--max-targets-per-child",
        type=int,
        default=V2_REVIEW_TARGET_DEFAULT,
        metavar="N",
        help="bounded v2 child target cap (1..25)",
    )
    parser.add_argument(
        "--explain-target",
        metavar="ID",
        help="include a bounded explanation for one selected target",
    )
    parsed = parser.parse_args(arguments)
    parsed.provided_options = frozenset(provided)
    return parsed


def require_v2_artifact_grammar(parsed: argparse.Namespace) -> str:
    if "--artifact-root" not in parsed.provided_options:
        raise UsageRefused("--artifact-root is required for v2 artifact modes")
    if "--artifact-dir" not in parsed.provided_options:
        raise UsageRefused("--artifact-dir is required for v2 artifact modes")
    if not 1 <= parsed.max_targets_per_child <= V2_REVIEW_TARGET_HARD_MAX:
        raise UsageRefused(
            f"--max-targets-per-child must be in 1..{V2_REVIEW_TARGET_HARD_MAX}"
        )
    return require_artifact_dir(parsed)


def require_artifact_dir(arguments: argparse.Namespace) -> str:
    if not arguments.artifact_dir:
        raise UsageRefused("--artifact-dir is required for modes that write artifacts")
    safe_relative(arguments.artifact_dir, label="artifact dir")
    return arguments.artifact_dir


def main(arguments: Sequence[str]) -> int:
    parsed = parse_arguments(arguments)
    v2_auxiliary_options = {
        "--review-receipts",
        "--prior-campaign",
        "--priorities",
        "--statuses",
        "--max-targets-per-child",
        "--explain-target",
    }

    if parsed.self_test_v2 or parsed.case_v2:
        if parsed.provided_options & {
            "--artifact-root",
            "--artifact-dir",
            *v2_auxiliary_options,
        }:
            raise UsageRefused("v2 test-runner modes do not accept planner options")
        manifest_v2 = load_case_manifest_v2()
        run_v2_self_tests(manifest_v2, selected=parsed.case_v2)
        return 0

    if parsed.review_plan or parsed.history_plan:
        manifest_v2 = load_case_manifest_v2()
        synopsis = v2_execute_live_plan(
            mode="review-plan" if parsed.review_plan else "history-plan",
            parsed=parsed,
            manifest=manifest_v2,
        )
        json_stdout(synopsis)
        return 0

    if parsed.replay:
        schema = peek_replay_schema(parsed.artifact_root, parsed.replay)
        if schema == V2_TERMINAL_SCHEMA:
            if parsed.provided_options & v2_auxiliary_options:
                raise UsageRefused("v2 replay does not accept live planner options")
            artifact_dir = require_v2_artifact_grammar(parsed)
            json_stdout(
                v2_replay_bundle(
                    artifact_root=parsed.artifact_root,
                    input_dir=parsed.replay,
                    output_dir=artifact_dir,
                )
            )
            return 0
        if schema != RUN_TERMINAL_SCHEMA:
            raise InputRefused("replay terminal schema is neither frozen v1 nor v2")
        if parsed.provided_options & v2_auxiliary_options:
            raise UsageRefused("v1 replay does not accept v2 planner options")
        case_manifest = load_case_manifest()
        artifact_dir = require_artifact_dir(parsed)
        result = replay_bundle(
            artifact_root=parsed.artifact_root,
            input_dir=parsed.replay,
            output_dir=artifact_dir,
            case_manifest_root=case_manifest["semantic_root"],
        )
        json_stdout(result)
        return 0

    if parsed.provided_options & v2_auxiliary_options:
        raise UsageRefused("v1 modes do not accept v2 planner options")
    case_manifest = load_case_manifest()

    if parsed.list:
        for sequence, row in enumerate(case_manifest["case"]):
            json_stdout(
                {
                    "schema": EVENT_SCHEMA,
                    "mode": "list",
                    "sequence": sequence,
                    "case": row["id"],
                    "role": row["role"],
                    "case_mode": row["mode"],
                    "expected_terminal": row["expected_terminal"],
                    "expected_exit": row["expected_exit"],
                    "mutation": row["mutation"],
                    "replay": row["replay"],
                    "authority": row["authority"],
                }
            )
        json_stdout(
            terminal_row(
                "Pass",
                mode="list",
                detail=f"listed {len(case_manifest['case'])} frozen cases",
                extra={
                    "sequence": len(case_manifest["case"]),
                    "case_manifest_root": case_manifest["semantic_root"],
                },
            )
        )
        return 0

    if parsed.self_test:
        run_self_tests(case_manifest)
        return 0

    if parsed.negative:
        row = next(
            (row for row in case_manifest["case"] if row["id"] == parsed.negative),
            None,
        )
        if row is None:
            raise UsageRefused(f"unknown case: {parsed.negative}")
        if row["role"] not in {"negative", "hostile"}:
            raise UsageRefused(f"{parsed.negative} is not a frozen negative case")
        run_self_tests(case_manifest, selected=parsed.negative)
        return 0

    snapshot = collect_live(case_manifest)

    if parsed.check:
        json_stdout(
            terminal_row(
                "Pass",
                mode="check",
                detail="live exact-set inventory is complete and classified",
                extra={
                    "issues": snapshot.inventory["counts"]["issues"],
                    "warnings": snapshot.inventory["counts"]["warnings"],
                    "inventory_root": snapshot.inventory["semantic_root"],
                    "case_manifest_root": case_manifest["semantic_root"],
                },
            )
        )
        return 0

    if parsed.closeout:
        warnings = snapshot.inventory["counts"]["warnings"]
        if warnings:
            raise NoData(
                f"zero-debt closeout unavailable: {warnings} warnings remain "
                f"across {snapshot.inventory['counts']['issues']} issues"
            )
        json_stdout(
            terminal_row(
                "Pass",
                mode="closeout",
                detail="zero exact template debt independently reconstructed",
                extra={"inventory_root": snapshot.inventory["semantic_root"]},
            )
        )
        return 0

    if parsed.inventory or parsed.plan:
        artifact_dir = require_artifact_dir(parsed)
        result = write_bundle(
            mode="plan" if parsed.plan else "inventory",
            artifact_root=parsed.artifact_root,
            artifact_dir=artifact_dir,
            snapshot=snapshot,
        )
        json_stdout(result)
        return 0

    if parsed.apply_manifest:
        artifact_dir = require_artifact_dir(parsed)
        result = apply_reviewed_manifest(
            manifest_rel=parsed.apply_manifest,
            artifact_root=parsed.artifact_root,
            artifact_dir=artifact_dir,
            snapshot=snapshot,
        )
        json_stdout(result)
        return 0

    raise UsageRefused("no mode selected")


try:
    main_exit = main(sys.argv[1:])
except HarnessError as error:
    mode = next(
        (
            sys.argv[index].lstrip("-")
            for index in range(1, len(sys.argv))
            if sys.argv[index].startswith("--")
        ),
        "unknown",
    )
    json_stdout(
        terminal_row(
            error.terminal,
            mode=mode,
            detail=str(error),
        )
    )
    raise SystemExit(TERMINAL_EXIT[error.terminal])
except KeyboardInterrupt:
    json_stdout(
        terminal_row(
            "CancelledDrained",
            mode="unknown",
            detail="keyboard cancellation drained before terminal publication",
        )
    )
    raise SystemExit(TERMINAL_EXIT["CancelledDrained"])
except SystemExit:
    raise
except BaseException as error:
    json_stdout(
        terminal_row(
            "InternalFault",
            mode="unknown",
            detail=f"{type(error).__name__}: {error}",
        )
    )
    raise SystemExit(TERMINAL_EXIT["InternalFault"])
else:
    raise SystemExit(main_exit)
PY
