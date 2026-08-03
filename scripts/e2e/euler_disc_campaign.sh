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

if [[ "$DRY_RUN" -eq 1 ]]; then
  stage plan "result=ready execution=$EXECUTION source_identity=$SOURCE_IDENTITY binary_sha256=$BINARY_SHA256 output=$OUTPUT stderr_log=$STDERR_LOG schema=$SCHEMA"
  printf '%s\n' 'euler_disc_campaign result=dry-run numerical_software_validation_only=true' >&2
  exit 0
fi

mkdir -p "$(dirname "$OUTPUT")" "$(dirname "$STDERR_LOG")"
printf 'euler_disc_campaign provenance execution=%s source_identity=%s binary_sha256=%s numerical_software_validation_only=true\n' \
  "$EXECUTION" "$SOURCE_IDENTITY" "$BINARY_SHA256" >"$STDERR_LOG"

stage execute "execution=$EXECUTION source_identity=$SOURCE_IDENTITY binary_sha256=$BINARY_SHA256 output=$OUTPUT stderr_log=$STDERR_LOG"
if [[ -n "$PREBUILT_BINARY" ]]; then
  if "$PREBUILT_BINARY" >"$OUTPUT" 2>>"$STDERR_LOG"; then
    RUN_STATUS=0
  else
    RUN_STATUS=$?
  fi
else
  RCH_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_euler_disc_campaign"
  if RCH_REQUIRE_REMOTE=1 RCH_NO_SELF_HEALING=1 rch --no-self-healing exec -- \
    env CARGO_TARGET_DIR="$RCH_TARGET_DIR" \
    cargo run --locked -p fs-euler-disc-e2e --bin euler_disc_campaign -- \
    >"$OUTPUT" 2>>"$STDERR_LOG"; then
    RUN_STATUS=0
  else
    RUN_STATUS=$?
  fi
fi

if [[ "$RUN_STATUS" -ne 0 ]]; then
  printf 'euler_disc_campaign stage=execute elapsed_s=%s result=nonzero-exit status=%s output_retained=%s stderr_log_retained=%s partial_output_rejected=true\n' \
    "$SECONDS" "$RUN_STATUS" "$OUTPUT" "$STDERR_LOG" >&2
  exit "$RUN_STATUS"
fi

stage validate-jsonl "schema=$SCHEMA digest_domain=$DIGEST_DOMAIN digest_scope=$DIGEST_SCOPE"
python3 - "$OUTPUT" "$SCHEMA" "$DIGEST_DOMAIN" "$DIGEST_SCOPE" "${EXPECTED_SCENARIOS[@]}" <<'PY'
import json
import re
import sys

try:
    import blake3
except ImportError as error:
    raise SystemExit(
        "JSONL refusal: repository-standard Python interpreter lacks the "
        f"blake3 module needed to verify the manifest digest: {error}"
    )

output_path, schema, digest_domain, digest_scope, *expected_scenarios = sys.argv[1:]
raw = open(output_path, "rb").read()
if not raw:
    raise SystemExit("JSONL refusal: output is empty")
if not raw.endswith(b"\n"):
    raise SystemExit("JSONL refusal: output lacks a final newline and is treated as partial")
try:
    text = raw.decode("utf-8")
except UnicodeDecodeError as error:
    raise SystemExit(f"JSONL refusal: output is not UTF-8: {error}")

lines = text.splitlines()
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
payload = b"\n".join(raw.splitlines()[:-1])
actual_digest = blake3.blake3(payload, derive_key_context=digest_domain).hexdigest()
if digest != actual_digest:
    raise SystemExit(
        "JSONL refusal: manifest digest mismatch "
        f"expected={digest} recomputed={actual_digest}"
    )

print(
    "JSONL validation passed: "
    f"records={len(records)} payload_records={record_count} digest={digest}",
    file=sys.stderr,
)
PY

stage complete "result=passed records=${#EXPECTED_SCENARIOS[@]} output=$OUTPUT stderr_log=$STDERR_LOG numerical_software_validation_only=true"
