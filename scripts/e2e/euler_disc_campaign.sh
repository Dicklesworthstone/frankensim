#!/usr/bin/env bash
# Deterministic E2E runner for the committed Euler-disc campaign executable.
#
# This validates JSONL structure and encoded numerical/software receipts only.
# It is not experimental validation, physical validation, mechanism evidence,
# or a prediction/ranking of any observed Euler-disc outcome.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

readonly SCHEMA="euler-disc-campaign-jsonl-v1"
readonly DIGEST_DOMAIN="org.frankensim.euler-disc-campaign-jsonl.v1"
readonly DIGEST_SCOPE="preceding-data-records-LF-joined-no-trailing-LF"
readonly -a EXPECTED_SCENARIOS=(
  "geometry-sharp-squat-disc"
  "geometry-filleted-squat-disc"
  "conservative-steady-oracle"
  "dynamic-unilateral-contact"
  "reduced-flexible-base"
  "contour-only-decay"
  "boundary-layer-only-decay"
  "combined-decay"
  "reduced-exterior-wrench-passivity"
  "campaign-complete"
)

OUTPUT=""
STDERR_LOG=""
PREBUILT_BINARY=""
DRY_RUN=0

usage() {
  printf '%s\n' \
    'usage: scripts/e2e/euler_disc_campaign.sh --output PATH --stderr-log PATH [options]' \
    '' \
    'Runs the committed euler_disc_campaign executable through strict-remote RCH,' \
    'or runs --prebuilt-binary directly. JSONL is written only to --output and' \
    'all program/RCH stderr is retained at --stderr-log.' \
    '' \
    'Required:' \
    '  --output PATH        New JSONL evidence path; existing paths are refused.' \
    '  --stderr-log PATH    New stderr log path; existing paths are refused.' \
    '' \
    'Options:' \
    '  --prebuilt-binary PATH  Executable to run instead of RCH/Cargo.' \
    '  --dry-run               Validate inputs and print the resolved execution plan;' \
    '                          do not create files, invoke RCH, Cargo, or the binary.' \
    '  -h, --help              Show this help.' \
    '' \
    'The RCH path requires root Cargo.toml, Cargo.lock, and crates/ to be tracked' \
    'and clean, so no transitive dirty Cargo input can masquerade as HEAD. A' \
    'prebuilt run records the binary SHA-256 as a digest-only source identity;' \
    'it does not establish the source revision that produced the binary.' \
    '' \
    'This runner is numerical/software validation only. It performs no experimental' \
    'validation and makes no physical-validation, mechanism, or video-fit claim.'
}

die() {
  printf 'euler_disc_campaign result=refused elapsed_s=%s reason=%s\n' "$SECONDS" "$*" >&2
  exit 2
}

stage() {
  printf 'euler_disc_campaign stage=%s elapsed_s=%s %s\n' "$1" "$SECONDS" "$2" >&2
}

sha256_file() {
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 -- "$1" | awk '{print $1}'
  elif command -v sha256sum >/dev/null 2>&1; then
    sha256sum -- "$1" | awk '{print $1}'
  else
    die 'prebuilt provenance requires shasum or sha256sum for SHA-256'
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --output)
      [[ $# -ge 2 && -n "${2:-}" ]] || die '--output requires a nonempty path'
      OUTPUT="$2"
      shift 2
      ;;
    --stderr-log)
      [[ $# -ge 2 && -n "${2:-}" ]] || die '--stderr-log requires a nonempty path'
      STDERR_LOG="$2"
      shift 2
      ;;
    --prebuilt-binary)
      [[ $# -ge 2 && -n "${2:-}" ]] || die '--prebuilt-binary requires a path'
      PREBUILT_BINARY="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown argument: $1"
      ;;
  esac
done

[[ -n "$OUTPUT" ]] || die '--output is required'
[[ -n "$STDERR_LOG" ]] || die '--stderr-log is required'
[[ "$OUTPUT" != "$STDERR_LOG" ]] || die '--output and --stderr-log must differ'
[[ ! -e "$OUTPUT" && ! -L "$OUTPUT" ]] || die "refusing to overwrite existing output: $OUTPUT"
[[ ! -e "$STDERR_LOG" && ! -L "$STDERR_LOG" ]] || die "refusing to overwrite existing stderr log: $STDERR_LOG"

if [[ -n "$PREBUILT_BINARY" ]]; then
  [[ -f "$PREBUILT_BINARY" && -x "$PREBUILT_BINARY" ]] || \
    die "--prebuilt-binary must name an executable file: $PREBUILT_BINARY"
  BINARY_SHA256="$(sha256_file "$PREBUILT_BINARY")"
  [[ "$BINARY_SHA256" =~ ^[0-9a-fA-F]{64}$ ]] || \
    die "could not obtain a SHA-256 digest for prebuilt binary: $PREBUILT_BINARY"
  EXECUTION="prebuilt"
  SOURCE_IDENTITY="prebuilt-binary-sha256:$BINARY_SHA256"
else
  command -v rch >/dev/null 2>&1 || die 'rch is required unless --prebuilt-binary is supplied'
  command -v git >/dev/null 2>&1 || die 'git is required for committed-source admission'
  git ls-files --error-unmatch crates/fs-euler-disc-e2e/src/bin/euler_disc_campaign.rs \
    >/dev/null 2>&1 || die 'campaign binary source is not committed/tracked; use --prebuilt-binary'
  git ls-files --error-unmatch Cargo.toml Cargo.lock >/dev/null 2>&1 || \
    die 'root Cargo.toml and Cargo.lock must be committed/tracked for RCH admission'
  if [[ -n "$(git status --porcelain --untracked-files=all -- Cargo.toml Cargo.lock crates)" ]]; then
    die 'Cargo.toml, Cargo.lock, or crates/ has uncommitted or untracked content; use --prebuilt-binary'
  fi
  SOURCE_REVISION="$(git rev-parse --verify HEAD)"
  SOURCE_IDENTITY="git-head:$SOURCE_REVISION;cargo-inputs=clean:Cargo.toml,Cargo.lock,crates"
  BINARY_SHA256="not-applicable-rch-build"
  EXECUTION="rch"
fi
RUNNER_SHA256="$(sha256_file "${BASH_SOURCE[0]}")"
[[ "$RUNNER_SHA256" =~ ^[0-9a-fA-F]{64}$ ]] || \
  die 'could not obtain a SHA-256 digest for the campaign runner'

if [[ "$DRY_RUN" -eq 1 ]]; then
  stage plan "result=ready execution=$EXECUTION source_identity=$SOURCE_IDENTITY binary_sha256=$BINARY_SHA256 runner_sha256=$RUNNER_SHA256 output=$OUTPUT stderr_log=$STDERR_LOG schema=$SCHEMA"
  printf '%s\n' 'euler_disc_campaign result=dry-run numerical_software_validation_only=true' >&2
  exit 0
fi

mkdir -p "$(dirname "$OUTPUT")" "$(dirname "$STDERR_LOG")"
set -o noclobber
if ! exec 8>"$OUTPUT"; then
  die "could not exclusively create output: $OUTPUT"
fi
if ! exec 9>"$STDERR_LOG"; then
  die "could not exclusively create stderr log: $STDERR_LOG"
fi
set +o noclobber
printf 'euler_disc_campaign provenance execution=%s source_identity=%s binary_sha256=%s runner_sha256=%s numerical_software_validation_only=true\n' \
  "$EXECUTION" "$SOURCE_IDENTITY" "$BINARY_SHA256" "$RUNNER_SHA256" >&9

stage execute "execution=$EXECUTION source_identity=$SOURCE_IDENTITY binary_sha256=$BINARY_SHA256 runner_sha256=$RUNNER_SHA256 output=$OUTPUT stderr_log=$STDERR_LOG"
if [[ -n "$PREBUILT_BINARY" ]]; then
  if "$PREBUILT_BINARY" >&8 2>&9; then
    RUN_STATUS=0
  else
    RUN_STATUS=$?
  fi
else
  RCH_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_euler_disc_campaign"
  if RCH_REQUIRE_REMOTE=1 RCH_NO_SELF_HEALING=1 rch --no-self-healing exec -- \
    env CARGO_TARGET_DIR="$RCH_TARGET_DIR" \
    cargo run --locked -p fs-euler-disc-e2e --bin euler_disc_campaign -- \
    >&8 2>&9; then
    RUN_STATUS=0
  else
    RUN_STATUS=$?
  fi
fi

RCH_TRANSPORT_STATUS="not-applicable"
RCH_ARTIFACT_RETRIEVAL="not-applicable"
RCH_RECEIPT_AUTHORITY="not-applicable"
NUMERICAL_STDOUT_SOURCE="direct-binary-stdout"
if [[ "$EXECUTION" == "rch" ]]; then
  RCH_TRANSPORT_STATUS="$RUN_STATUS"
  if [[ "$RUN_STATUS" -ne 0 && "$RUN_STATUS" -ne 102 ]]; then
    printf 'euler_disc_campaign stage=execute elapsed_s=%s result=nonzero-exit status=%s output_retained=%s stderr_log_retained=%s partial_output_rejected=true\n' \
      "$SECONDS" "$RUN_STATUS" "$OUTPUT" "$STDERR_LOG" >&2
    exit "$RUN_STATUS"
  fi

  # RCH owns the remote process' stdout while it orchestrates the build.  Some
  # versions replay that stdout into the retained diagnostic transcript rather
  # than the caller's stdout.  Recover only byte-exact campaign records, and
  # only after the transcript proves that the remote command itself exited 0.
  # Exit 102 is admitted solely for RCH-E309: the numerical stdout is complete,
  # but retrieval of the separately built executable artifact is incomplete.
  if python3 - "$OUTPUT" "$STDERR_LOG" 8 9 "$SCHEMA" "$RUN_STATUS" 8>&8 9>&9 <<'PY'
import os
import re
import stat
import sys

output_path = sys.argv[1]
stderr_path = sys.argv[2]
output_fd = int(sys.argv[3])
stderr_fd = int(sys.argv[4])
schema = sys.argv[5].encode("utf-8")
status = int(sys.argv[6])


def open_bound_reader(path: str, retained_fd: int) -> int:
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        reader_fd = os.open(path, flags)
    except OSError as error:
        raise SystemExit(f"RCH stdout recovery refusal: cannot securely open {path!r}: {error}")
    reader = os.fstat(reader_fd)
    retained = os.fstat(retained_fd)
    if not stat.S_ISREG(reader.st_mode) or not stat.S_ISREG(retained.st_mode):
        raise SystemExit("RCH stdout recovery refusal: retained path is not a regular file")
    if (reader.st_dev, reader.st_ino) != (retained.st_dev, retained.st_ino):
        raise SystemExit("RCH stdout recovery refusal: retained path identity changed")
    return reader_fd


def read_fd(fd: int) -> bytes:
    try:
        size = os.fstat(fd).st_size
    except OSError as error:
        raise SystemExit(f"RCH stdout recovery refusal: retained fd {fd} is unavailable: {error}")
    chunks = []
    offset = 0
    while offset < size:
        chunk = os.pread(fd, min(1024 * 1024, size - offset), offset)
        if not chunk:
            raise SystemExit("RCH stdout recovery refusal: retained file was truncated while reading")
        chunks.append(chunk)
        offset += len(chunk)
    try:
        final_size = os.fstat(fd).st_size
    except OSError as error:
        raise SystemExit(f"RCH stdout recovery refusal: retained fd {fd} closed while reading: {error}")
    if final_size != size:
        raise SystemExit("RCH stdout recovery refusal: retained file changed while reading")
    return b"".join(chunks)


def write_fd(fd: int, payload: bytes) -> None:
    if os.fstat(fd).st_size != 0:
        raise SystemExit("RCH stdout recovery refusal: output became nonempty before recovery")
    os.ftruncate(fd, 0)
    offset = 0
    while offset < len(payload):
        written = os.pwrite(fd, payload[offset:], offset)
        if written <= 0:
            raise SystemExit("RCH stdout recovery refusal: could not write recovered output")
        offset += written
    os.fsync(fd)

output = read_fd(open_bound_reader(output_path, output_fd))
transcript = read_fd(open_bound_reader(stderr_path, stderr_fd))
prefix = b'{"schema":"' + schema + b'"'
lines = transcript.split(b"\n")
record_frames = []
record_indices = []
for index, line in enumerate(lines):
    if not line.startswith(prefix):
        continue
    if line.endswith(b"\r") or (index == len(lines) - 1 and not transcript.endswith(b"\n")):
        raise SystemExit("RCH stdout recovery refusal: campaign record is not LF-framed")
    record_frames.append(line + b"\n")
    record_indices.append(index)
recovered = b"".join(record_frames)

ansi_csi = re.compile(rb"\x1b\[[0-?]*[ -/]*[@-~]")
try:
    clean_lines = [ansi_csi.sub(b"", line).decode("utf-8") for line in lines]
except UnicodeDecodeError as error:
    raise SystemExit(f"RCH stdout recovery refusal: transcript is not UTF-8: {error}")

remote_receipts = [
    index
    for index, line in enumerate(clean_lines)
    if "rch::hook::transfer_orchestration" in line
    and re.search(r"Remote command finished: exit=0 in [0-9]+ms$", line)
]
remote_succeeded = len(remote_receipts) == 1

artifact_retrieval_failed = False
if status == 102 and record_indices and remote_succeeded:
    selected = []
    submitted = []
    retrievals = []
    summaries = []
    e309_receipts = []
    terminal_receipts = []
    for index, line in enumerate(clean_lines):
        match = re.search(r"Selected worker: ([A-Za-z0-9._-]+) at .+", line)
        if match:
            selected.append((index, match.group(1)))
        match = re.fullmatch(r"\[\*\] Job (j-[0-9]+) submitted to ([A-Za-z0-9._-]+) \(.+\)", line)
        if match:
            submitted.append((index, match.group(1), match.group(2)))
        match = re.search(r"Retrieving artifacts from .+ on ([A-Za-z0-9._-]+)$", line)
        if match:
            retrievals.append((index, match.group(1)))
        match = re.fullmatch(r"\[W\] Worker: ([A-Za-z0-9._-]+) \| Job: (j-[0-9]+)", line)
        if match:
            summaries.append((index, match.group(2), match.group(1)))
        match = re.fullmatch(
            r"\[RCH\] RCH-E309 remote compile on ([A-Za-z0-9._-]+) SUCCEEDED "
            r"but build artifacts could not be retrieved .+ \(exit 102\); re-run to "
            r"rebuild, or check connectivity to the worker\.",
            line,
        )
        if match:
            e309_receipts.append((index, match.group(1)))
        match = re.fullmatch(r"\[RCH\] remote ([A-Za-z0-9._-]+) failed \(exit 102\)", line)
        if match:
            terminal_receipts.append((index, match.group(1)))

    last_nonempty = max(index for index, line in enumerate(clean_lines) if line)
    running_index = record_indices[0] - 1
    running_bound = running_index >= 0 and re.fullmatch(
        r"\s+Running `[^`\n]*/euler_disc_campaign`", clean_lines[running_index]
    )
    if (
        len(selected)
        == len(submitted)
        == len(retrievals)
        == len(summaries)
        == len(e309_receipts)
        == len(terminal_receipts)
        == 1
        and running_bound
    ):
        selected_index, selected_worker = selected[0]
        submitted_index, submitted_job, submitted_worker = submitted[0]
        retrieval_index, retrieval_worker = retrievals[0]
        summary_index, summary_job, summary_worker = summaries[0]
        e309_index, e309_worker = e309_receipts[0]
        terminal_index, terminal_worker = terminal_receipts[0]
        remote_index = remote_receipts[0]
        artifact_retrieval_failed = (
            submitted_job == summary_job
            and selected_worker
            == submitted_worker
            == retrieval_worker
            == summary_worker
            == e309_worker
            == terminal_worker
            and selected_index < submitted_index < running_index < record_indices[0]
            and record_indices[-1]
            < remote_index
            < retrieval_index
            < summary_index
            < e309_index
            < terminal_index
            and terminal_index == last_nonempty
        )

if status == 102 and not (remote_succeeded and artifact_retrieval_failed and record_frames):
    raise SystemExit(
        "RCH stdout recovery refusal: exit 102 lacks the ordered job/worker, "
        "executed-binary, remote-exit-0, RCH-E309, terminal-exit-102 receipt "
        "chain or LF-framed campaign records"
    )
if status == 0 and not output:
    raise SystemExit(
        "RCH stdout recovery refusal: empty direct stdout is not recoverable "
        "for ordinary RCH status 0"
    )
if output and record_frames and output != recovered:
    raise SystemExit(
        "RCH stdout recovery refusal: direct stdout and transcript campaign "
        "records differ"
    )

if not output:
    write_fd(output_fd, recovered)
    if status == 102:
        raise SystemExit(0)
    raise SystemExit(1)
elif record_frames:
    raise SystemExit(11)
else:
    raise SystemExit(12)
PY
  then
    RCH_RECOVERY_STATUS=0
  else
    RCH_RECOVERY_STATUS=$?
  fi

  case "$RCH_RECOVERY_STATUS" in
    0)
      RCH_RECOVERY_RESULT="recovered-transcript-rch-e309"
      ;;
    11)
      RCH_RECOVERY_RESULT="direct-and-transcript-identical"
      ;;
    12)
      RCH_RECOVERY_RESULT="direct-stdout"
      ;;
    *)
      printf 'euler_disc_campaign stage=recover-rch-stdout elapsed_s=%s result=refused status=%s output_retained=%s stderr_log_retained=%s partial_output_rejected=true\n' \
        "$SECONDS" "$RUN_STATUS" "$OUTPUT" "$STDERR_LOG" >&2
      exit 2
      ;;
  esac

  case "$RCH_RECOVERY_RESULT" in
    recovered-transcript-rch-e309)
      RCH_ARTIFACT_RETRIEVAL="incomplete-rch-e309"
      RCH_RECEIPT_AUTHORITY="local-exclusive-transcript-not-cryptographic-attestation"
      NUMERICAL_STDOUT_SOURCE="verified-rch-transcript"
      ;;
    direct-and-transcript-identical)
      if [[ "$RUN_STATUS" -eq 102 ]]; then
        RCH_ARTIFACT_RETRIEVAL="incomplete-rch-e309"
        RCH_RECEIPT_AUTHORITY="local-exclusive-transcript-not-cryptographic-attestation"
      else
        RCH_ARTIFACT_RETRIEVAL="complete"
        RCH_RECEIPT_AUTHORITY="rch-process-exit-status"
      fi
      NUMERICAL_STDOUT_SOURCE="direct-and-transcript-identical"
      ;;
    direct-stdout)
      RCH_ARTIFACT_RETRIEVAL="complete"
      RCH_RECEIPT_AUTHORITY="rch-process-exit-status"
      NUMERICAL_STDOUT_SOURCE="direct-rch-stdout"
      ;;
    *)
      die "unexpected RCH stdout recovery result: $RCH_RECOVERY_RESULT"
      ;;
  esac
  printf 'euler_disc_campaign recovery_receipt result=accepted transport_status=%s numerical_stdout_source=%s artifact_retrieval=%s receipt_authority=%s runner_sha256=%s\n' \
    "$RCH_TRANSPORT_STATUS" "$NUMERICAL_STDOUT_SOURCE" "$RCH_ARTIFACT_RETRIEVAL" \
    "$RCH_RECEIPT_AUTHORITY" "$RUNNER_SHA256" >&9
  stage recover-rch-stdout \
    "result=accepted transport_status=$RCH_TRANSPORT_STATUS numerical_stdout_source=$NUMERICAL_STDOUT_SOURCE artifact_retrieval=$RCH_ARTIFACT_RETRIEVAL receipt_authority=$RCH_RECEIPT_AUTHORITY"
elif [[ "$RUN_STATUS" -ne 0 ]]; then
  printf 'euler_disc_campaign stage=execute elapsed_s=%s result=nonzero-exit status=%s output_retained=%s stderr_log_retained=%s partial_output_rejected=true\n' \
    "$SECONDS" "$RUN_STATUS" "$OUTPUT" "$STDERR_LOG" >&2
  exit "$RUN_STATUS"
fi

stage validate-jsonl "schema=$SCHEMA digest_domain=$DIGEST_DOMAIN digest_scope=$DIGEST_SCOPE"
python3 - \
  "$OUTPUT" 8 "$STDERR_LOG" 9 \
  "$EXECUTION" "$SOURCE_IDENTITY" "$BINARY_SHA256" "$RUNNER_SHA256" \
  "$NUMERICAL_STDOUT_SOURCE" "$RCH_TRANSPORT_STATUS" "$RCH_ARTIFACT_RETRIEVAL" \
  "$RCH_RECEIPT_AUTHORITY" "$SCHEMA" "$DIGEST_DOMAIN" "$DIGEST_SCOPE" \
  "${EXPECTED_SCENARIOS[@]}" 8>&8 9>&9 <<'PY'
import hashlib
import json
import os
import re
import stat
import sys

try:
    import blake3
except ImportError as error:
    raise SystemExit(
        "JSONL refusal: repository-standard Python interpreter lacks the "
        f"blake3 module needed to verify the manifest digest: {error}"
    )

(
    output_path,
    output_fd_text,
    stderr_path,
    stderr_fd_text,
    execution,
    source_identity,
    binary_sha256,
    runner_sha256,
    numerical_stdout_source,
    rch_transport_status,
    artifact_retrieval,
    receipt_authority,
    schema,
    digest_domain,
    digest_scope,
    *expected_scenarios,
) = sys.argv[1:]
output_fd = int(output_fd_text)
stderr_fd = int(stderr_fd_text)


def verify_path_identity(path_name: str, descriptor_fd: int, label: str):
    descriptor = os.fstat(descriptor_fd)
    path = os.stat(path_name, follow_symlinks=False)
    if not stat.S_ISREG(descriptor.st_mode) or not stat.S_ISREG(path.st_mode):
        raise SystemExit(f"JSONL refusal: {label} path or retained descriptor is not a regular file")
    if (descriptor.st_dev, descriptor.st_ino) != (path.st_dev, path.st_ino):
        raise SystemExit(f"JSONL refusal: {label} path no longer names the exclusively created file")
    return descriptor


def open_bound_reader() -> int:
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        reader_fd = os.open(output_path, flags)
    except OSError as error:
        raise SystemExit(f"JSONL refusal: cannot securely open output: {error}")
    reader = os.fstat(reader_fd)
    descriptor = os.fstat(output_fd)
    if (reader.st_dev, reader.st_ino) != (descriptor.st_dev, descriptor.st_ino):
        raise SystemExit("JSONL refusal: output reader does not match retained descriptor")
    return reader_fd


def read_output(reader_fd: int) -> bytes:
    size = os.fstat(reader_fd).st_size
    chunks = []
    offset = 0
    while offset < size:
        chunk = os.pread(reader_fd, min(1024 * 1024, size - offset), offset)
        if not chunk:
            raise SystemExit("JSONL refusal: output was truncated while reading")
        chunks.append(chunk)
        offset += len(chunk)
    if os.fstat(reader_fd).st_size != size:
        raise SystemExit("JSONL refusal: output changed while reading")
    return b"".join(chunks)


verify_path_identity(output_path, output_fd, "output")
raw = read_output(open_bound_reader())
if not raw:
    raise SystemExit("JSONL refusal: output is empty")
if not raw.endswith(b"\n"):
    raise SystemExit("JSONL refusal: output lacks a final newline and is treated as partial")
if b"\r" in raw:
    raise SystemExit("JSONL refusal: output must use exact LF framing, not CR or CRLF")
try:
    text = raw.decode("utf-8")
except UnicodeDecodeError as error:
    raise SystemExit(f"JSONL refusal: output is not UTF-8: {error}")

lines = text[:-1].split("\n")
if not lines or any(not line.strip() for line in lines):
    raise SystemExit("JSONL refusal: blank record line")
records = []
for number, line in enumerate(lines, start=1):
    try:
        record = json.loads(line)
    except json.JSONDecodeError as error:
        raise SystemExit(f"JSONL refusal: line {number} is invalid JSON: {error}")
    if not isinstance(record, dict):
        raise SystemExit(f"JSONL refusal: line {number} must be a JSON object")
    records.append(record)

actual_scenarios = [record.get("scenario") for record in records]
if actual_scenarios != expected_scenarios:
    raise SystemExit(
        "JSONL refusal: scenario sequence (the producer's record kind) differs; "
        f"expected={expected_scenarios!r} actual={actual_scenarios!r}"
    )
for number, record in enumerate(records, start=1):
    if record.get("schema") != schema:
        raise SystemExit(f"JSONL refusal: line {number} has wrong schema")
    for field in ("model", "source", "authority", "units", "budget", "terminal", "residual", "no_claim"):
        if field not in record:
            raise SystemExit(f"JSONL refusal: line {number} lacks required field {field}")
    if "numerical-slice-only:no-physical-validation-or-target-ranking" not in record["no_claim"]:
        raise SystemExit(f"JSONL refusal: line {number} lacks numerical-only no-claim receipt")

manifest = records[-1]
if manifest.get("scenario") != "campaign-complete" or manifest.get("terminal") != "completed":
    raise SystemExit("JSONL refusal: final record is not a completed campaign manifest")
if any(key in manifest for key in ("inputs", "powers_w", "work_j")) and (
    manifest.get("powers_w") != {} or manifest.get("work_j") != {}
):
    raise SystemExit("JSONL refusal: manifest power/work fields must be empty")
record_count = len(records) - 1
budget = manifest.get("budget")
residual = manifest.get("residual")
if not isinstance(budget, dict) or budget.get("record_count") != record_count:
    raise SystemExit("JSONL refusal: manifest budget.record_count does not match preceding records")
if not isinstance(residual, dict) or residual.get("record_count") != record_count:
    raise SystemExit("JSONL refusal: manifest residual.record_count does not match preceding records")
digest = manifest.get("digest_blake3")
if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
    raise SystemExit("JSONL refusal: manifest digest_blake3 must be 64 lowercase hexadecimal characters")
if manifest.get("digest_domain") != digest_domain:
    raise SystemExit(
        "JSONL refusal: manifest digest_domain differs; "
        f"expected={digest_domain!r} actual={manifest.get('digest_domain')!r}"
    )
if manifest.get("digest_scope") != digest_scope:
    raise SystemExit(
        "JSONL refusal: manifest digest_scope differs; "
        f"expected={digest_scope!r} actual={manifest.get('digest_scope')!r}"
    )
payload = b"\n".join(raw[:-1].split(b"\n")[:-1])
actual_digest = blake3.blake3(payload, derive_key_context=digest_domain).hexdigest()
if digest != actual_digest:
    raise SystemExit(
        "JSONL refusal: manifest digest mismatch "
        f"expected={digest} recomputed={actual_digest}"
    )

output_identity = verify_path_identity(output_path, output_fd, "output")
stderr_identity = verify_path_identity(stderr_path, stderr_fd, "stderr log")
final_receipt = {
    "schema": "euler-disc-campaign-runner-receipt-v1",
    "result": "passed",
    "execution": execution,
    "source_identity": source_identity,
    "binary_sha256": binary_sha256,
    "runner_sha256": runner_sha256,
    "output_sha256": hashlib.sha256(raw).hexdigest(),
    "output_identity": {"device": output_identity.st_dev, "inode": output_identity.st_ino},
    "stderr_log_identity": {"device": stderr_identity.st_dev, "inode": stderr_identity.st_ino},
    "record_count": len(records),
    "numerical_stdout_source": numerical_stdout_source,
    "rch_transport_status": rch_transport_status,
    "artifact_retrieval": artifact_retrieval,
    "receipt_authority": receipt_authority,
    "numerical_software_validation_only": True,
}
receipt_bytes = (
    json.dumps(final_receipt, sort_keys=True, separators=(",", ":")) + "\n"
).encode("utf-8")
written = 0
while written < len(receipt_bytes):
    count = os.write(stderr_fd, receipt_bytes[written:])
    if count <= 0:
        raise SystemExit("JSONL refusal: could not append the retained final receipt")
    written += count
os.fsync(stderr_fd)
verify_path_identity(output_path, output_fd, "output")
verify_path_identity(stderr_path, stderr_fd, "stderr log")

print(
    "JSONL validation passed: "
    f"records={len(records)} payload_records={record_count} digest={digest}",
    file=sys.stderr,
)
PY

stage complete "result=passed records=${#EXPECTED_SCENARIOS[@]} output=$OUTPUT stderr_log=$STDERR_LOG numerical_stdout_source=$NUMERICAL_STDOUT_SOURCE rch_transport_status=$RCH_TRANSPORT_STATUS artifact_retrieval=$RCH_ARTIFACT_RETRIEVAL receipt_authority=$RCH_RECEIPT_AUTHORITY numerical_software_validation_only=true"
