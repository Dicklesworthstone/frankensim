#!/usr/bin/env bash
# M8 data checks. Append each later item's named tests to EXPECTED_TESTS.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
readonly REPO_ROOT
readonly -a EXPECTED_TESTS=(
  "fs-vvreg:instrument_data:ledger_header_and_column_count_are_exact"
  "fs-vvreg:instrument_data:every_required_field_refuses_when_blank"
  "fs-vvreg:instrument_data:unknown_enumerated_values_refuse"
  "fs-vvreg:instrument_data:duplicate_source_id_refuses"
  "fs-vvreg:instrument_data:field_byte_cap_and_cap_plus_one"
  "fs-vvreg:instrument_data:decision_consistency_rules_refuse_each_violation"
  "fs-vvreg:instrument_data:figure_digitization_requires_error_for_acceptance_but_records_fixtures"
  "fs-vvreg:instrument_data:committed_graph_read_packs_are_fixture_only"
  "fs-vvreg:instrument_data:identifiers_dates_and_artifact_paths_are_checked"
  "fs-vvreg:instrument_data:row_permutation_is_lookup_invariant"
  "fs-vvreg:instrument_data:committed_ledger_parses_and_every_artifact_resolves"
  "fs-vvreg:instrument_data:committed_ledger_has_no_orphan_instrument_data_files"
  "fs-vvreg:instrument_data:committed_clarinet_grader_receipts_are_fixture_only"
  "fs-vvreg:instrument_data:committed_ledger_seed_covers_existing_music_sources"
  "fs-vvreg:instrument_data:falsifier_licence_mismatch_is_caught"
  "python:scripts/tests/test_instrument_references.py:ReferenceTests.test_schema_and_required_fields"
  "python:scripts/tests/test_instrument_references.py:ReferenceTests.test_ids_and_urls_unique"
  "python:scripts/tests/test_instrument_references.py:ReferenceTests.test_every_uiowa_link_is_manifested_or_excluded"
  "python:scripts/tests/test_instrument_references.py:ReferenceTests.test_mono_stereo_pairing_is_explicit"
  "python:scripts/tests/test_instrument_references.py:ReferenceTests.test_stereo_not_gradable_until_derivation_documented"
  "python:scripts/tests/test_instrument_references.py:ReferenceTests.test_chain_assignment_by_era"
  "python:scripts/tests/test_instrument_references.py:ReferenceTests.test_ambient_entry_or_absence_per_session"
  "python:scripts/tests/test_instrument_references.py:ReferenceTests.test_ledger_rows_resolve_both_ways"
  "python:scripts/tests/test_instrument_references.py:ReferenceTests.test_measured_content_defects_are_recorded"
  "python:scripts/tests/test_instrument_references.py:ReferenceTests.test_fetch_refuses_data_volume_destination"
  "python:scripts/tests/test_instrument_references.py:ReferenceTests.test_fetch_refuses_mismatched_existing_file_and_deletes_nothing"
  "python:scripts/tests/test_instrument_references.py:ReferenceTests.test_fetch_verify_mode_never_downloads_when_present"
  "python:scripts/tests/test_instrument_references.py:ReferenceTests.test_fetch_refuses_escaping_paths"
  "python:scripts/tests/test_instrument_references.py:ReferenceTests.test_fetch_checks_remote_bytes_before_retaining"
)

exec python3 -I - "$REPO_ROOT" "${1:---check}" "${EXPECTED_TESTS[@]}" <<'PY'
import datetime
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

root, mode, *expected = sys.argv[1:]
root = Path(root)
if mode == "--list":
    print("\n".join(expected))
    sys.exit(0)
if mode not in ("--check", "--run"):
    sys.exit("usage: instrument_data.sh [--list|--check|--run]")

# Group exact functions on one native target; module filters run separately.
groups = {}
for entry in expected:
    crate, target, test = entry.split(":", 2)
    if crate == "python":
        source = root / target
        sources = [source]
        needle = "def " + test.split(".")[-1] + "("
    else:
        source = root / "crates" / crate / "tests" / (target + ".rs")
        sources = [source, *source.with_suffix("").rglob("*.rs")]
        needle = ("mod " + test.rstrip(":").split("::")[-1]) if test.endswith("::") else ("fn " + test.split("::")[-1] + "(")
    if not source.is_file() or not any(path.is_file() and needle in path.read_text() for path in sources):
        sys.exit(f"expected test missing: {entry} ({source})")
    groups.setdefault((crate, target, test.endswith("::")), []).append(test)
if len(set(expected)) != len(expected):
    sys.exit("duplicate EXPECTED_TESTS entry")
if mode == "--check":
    print(f"instrument_data check: {len(expected)} expected tests present; no builds executed")
    sys.exit(0)
if any(crate != "python" for crate, _, _ in groups):
    if not shutil.which("rch") or not shutil.which("cargo"):
        sys.exit("rch and cargo are required for native --run")

def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

base = os.environ.get("IRPDATA_RECEIPT_BASE") or os.environ.get("RCH_TARGET_BASE") or tempfile.gettempdir()
output = Path(tempfile.mkdtemp(prefix="frankensim-instrument-data-", dir=base))
head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
host = subprocess.check_output(["uname", "-n"], text=True).strip()
started = utc()
summary = {
    "schema": "frankensim.instrument-data.e2e.v1", "started_utc": started,
    "host": host, "git_head": head, "expected_tests": expected,
    "command_file": "command.txt", "stdout_log": "stdout.log",
    "stderr_log": "stderr.log", "events_log": "events.jsonl",
    "verdict": "running", "producer_exit": None,
}
def save_summary():
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

save_summary()
print(f"instrument data E2E receipt (live logs): {output}", flush=True)
failures = []
exits = []
required = {"bead", "test", "source_ids", "artifacts", "rows_checked",
            "refusals_exercised", "ladder", "tolerance", "observed", "verdict", "elapsed_ms"}
env = dict(os.environ, RCH_REQUIRE_REMOTE="1")
target_base = env.get("RCH_TARGET_BASE") or env.get("TMPDIR") or "/tmp"
with (output / "command.txt").open("w") as commands, \
     (output / "stdout.log").open("w+b") as stdout, \
     (output / "stderr.log").open("w+b") as stderr, \
     (output / "events.jsonl").open("w") as events:
    def log(payload, command):
        events.write(json.dumps(dict(payload, utc=utc(), host=host, git_head=head, command=command)) + "\n")
        events.flush()

    for (crate, target, module), tests in groups.items():
        if crate == "python":
            command = ["python3", "-I", str(root / target), "-v", *tests]
        else:
            command = ["rch", "exec", "--no-self-healing", "--", "env",
                       "CARGO_BUILD_JOBS=12", f"CARGO_TARGET_DIR={target_base}/rch_target_frankensim_test",
                       "cargo", "test", "--locked", "-p", crate, "--test", target,
                       "--", "--test-threads=1", "--show-output"]
            if not module:
                command.append("--exact")
            command.extend(tests)
        rendered = shlex.join(command)
        commands.write(f"utc={utc()}\nRCH_REQUIRE_REMOTE=1\ncommand={rendered}\n")
        commands.flush()
        log({"event": "start", "tests": tests}, rendered)
        out_start, err_start = stdout.tell(), stderr.tell()
        try:
            result = subprocess.run(command, cwd=root, env=env, stdout=stdout, stderr=stderr, check=False)
            code = result.returncode
        except OSError as error:
            stderr.write((str(error) + "\n").encode())
            code = 127
        exits.append(code)
        stdout.flush()
        stderr.flush()
        stdout.seek(out_start)
        stderr.seek(err_start)
        transcript = stdout.read().decode(errors="replace") + "\n" + stderr.read().decode(errors="replace")
        stdout.seek(0, 2)
        stderr.seek(0, 2)
        log({"event": "producer-exit", "status": code}, rendered)
        if code != 0:
            failures.append(f"producer exit {code}: {crate}:{target}")
        if crate == "python":
            if not re.search(r"^Ran [1-9][0-9]* tests? in ", transcript, re.M) or not re.search(r"^OK$", transcript, re.M):
                failures.append(f"missing nonempty unittest success: {target}")
            for test in tests:
                method, cls = test.split(".")[-1], test.split(".")[0]
                if not re.search(r"^" + re.escape(method) + r" \([^\n]*\b" + re.escape(cls) + r"\." + re.escape(method) + r"\)(?:\n[^\n]*)? \.\.\. ok$", transcript, re.M):
                    failures.append(f"missing test success: {test}")
        else:
            if not re.search(r"^test result: ok\. [1-9][0-9]* passed", transcript, re.M):
                failures.append(f"missing nonempty libtest success: {crate}:{target}")
            for test in tests:
                pattern = re.escape(test) + (r"\S+" if module else "")
                if not re.search(r"^test " + pattern + r" \.\.\. ok$", transcript, re.M):
                    failures.append(f"missing test success: {test}")
        emitted = set()
        for line in transcript.splitlines():
            if not line.startswith("IRPDATA "):
                continue
            try:
                payload = json.loads(line[len("IRPDATA "):])
                if not required <= payload.keys():
                    raise ValueError("missing IRPDATA fields")
                emitted.add(payload["test"])
                if payload["verdict"] not in ("pass", "deferred-to-M8-W6a"):
                    failures.append(f"artifact check failed: {payload['test']}")
                log(payload, rendered)
            except (ValueError, TypeError, AttributeError) as error:
                failures.append(f"invalid IRPDATA: {error}")
        for test in tests:
            if not any(name.startswith(test) if module else name == test for name in emitted):
                failures.append(f"missing IRPDATA event: {test}")
    log({"event": "complete", "failures": failures}, "scripts/e2e/instrument_data.sh --run")
summary.update(finished_utc=utc(), producer_exits=exits,
               producer_exit=next((code for code in exits if code), 0),
               verdict="failed" if failures else "pass", failures=failures)
save_summary()
print(f"instrument data E2E {summary['verdict']}; retained receipt: {output}")
for failure in failures:
    print(failure, file=sys.stderr)
sys.exit(1 if failures else 0)
PY
