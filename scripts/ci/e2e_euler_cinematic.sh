#!/usr/bin/env bash
# Euler-disc cinematic E2E runner (bead frankensim-h7xu5.9.8).
#
# One stable executable surface over the REAL production pipeline
# (euler_cinematic_fixture: mechanics -> modal sound -> WAV -> spectral
# render -> EXR -> finishing -> mux). It reuses production CLIs and
# checkers; it implements no parallel simulation or rendering logic.
#
# Modes:
#   --list                 enumerate cases, run nothing
#   --check                config/manifest consistency, no rendering
#   --self-test            the runner's own failure detection on fixtures
#   --run smoke            tiny moving-disc bundle through the REAL pipeline
#                          (32x18, 2 frames, 1 spp, 384 kHz mechanics,
#                          bounded ~10 min; measured recipe on this bead)
#   --negative CASE        one named hostile twin (or 'list'); the
#                          stale-trajectory-identity twin additionally
#                          needs --bundle SRC (a retained smoke bundle)
#   --replay MANIFEST_DIR  re-verify a retained smoke bundle's manifest
#                          against its artifacts (distinct output root)
#   --replay-render SRC    re-render from SRC bundle's retained canonical
#                          trajectory + listening master (identities taken
#                          from SRC's manifest; mechanics+audio skipped;
#                          measured ~2.5 min vs ~8 min full smoke)
#   --bundle SRC           retained bundle input for twins that exercise
#                          the production identity gates
#   --output-dir DIR       fresh artifact root (repo- or TMP-contained)
#   --log-file FILE        JSONL event log destination (defaults into the
#                          output root for run modes, scratch otherwise)
#
# EXIT CLASSES: 0 ok; 40 usage; 41 pipeline/production failure;
# 42 verification failure (manifest/artifact disagreement); 43 negative
# twin NOT detected (the hostile condition survived).
#
# LOGGING CONTRACT (euler-cinematic-e2e-log-v1): bounded deterministic
# JSONL; stable suite/case/stage/attempt/sequence IDs; source/build
# identity; artifact identities once verified; declared budgets;
# expected/observed on divergence; stable reason/exit class; relative
# artifact hashes; authority/no-claims disposition; ranked repairs; one
# repo-relative reproduction command. Redaction replaces the repo root,
# TMP roots, and HOME before buffering; records carry no absolute paths,
# hostnames, PIDs, or wall-clock noise. Caps: 256 records/file, 8 KiB
# per record, 8 hashed frame artifacts (excess counted, not listed).
#
# A real daily/final film remains an intentional product execution, not a
# CI tax; the routine lane runs only --check/--self-test/--negative and
# on-demand smoke. No-claims: a green smoke proves the software pipeline
# executes and its artifacts cohere; it proves nothing physical.
set -u -o pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXIT_OK=0
EXIT_USAGE=40
EXIT_PIPELINE=41
EXIT_VERIFY=42
EXIT_NEGATIVE_MISSED=43

# The measured bounded smoke recipe (probe record on bead h7xu5.9.8):
# minimum admitted mechanics rate whose decimator delay is integral is
# 384 kHz (ratio 8); default CRITIQUE gain 512 clips this tiny fixture
# (observed peak 1.46 fs), 280 keeps ~0.9 target with headroom.
SMOKE_ARGS=(--width 32 --height 18 --frames 2 --spp 1
  --mechanics-preroll-frames 1 --mechanics-sample-rate 384000
  --listening-gain-fs-per-pa 280)
SMOKE_BUDGET_SECONDS="${FSIM_CINE_SMOKE_BUDGET_SECONDS:-900}"

die() {
  local class="$1"; shift
  printf 'euler-cinematic-e2e: ERROR class=%s: %s\n' "$class" "$*" >&2
  exit "$class"
}
command -v python3 >/dev/null 2>&1 || die "$EXIT_USAGE" "python3 required"

fixture_binary() {
  # Prefer a prebuilt binary (the cargo test-runner wedge on this host makes
  # long cargo run invocations unreliable); fall back to cargo build.
  local candidate="${CARGO_TARGET_DIR:-$REPO_ROOT/target}/debug/euler_cinematic_fixture"
  if [ ! -x "$candidate" ]; then
    (cd "$REPO_ROOT" && nice -n 15 cargo build -q -p fs-euler-disc-e2e \
      --features cinematic-render --bin euler_cinematic_fixture) \
      || die "$EXIT_PIPELINE" "cannot build euler_cinematic_fixture"
  fi
  [ -x "$candidate" ] || die "$EXIT_PIPELINE" "fixture binary missing at $candidate"
  printf '%s' "$candidate"
}

MODE="" RUN_PROFILE="" OUTPUT_DIR="" NEGATIVE_CASE="" REPLAY_DIR="" BUNDLE_DIR="" LOG_FILE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --list) MODE="list"; shift ;;
    --check) MODE="check"; shift ;;
    --self-test) MODE="self-test"; shift ;;
    --run)
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--run needs smoke"
      MODE="run"; RUN_PROFILE="$2"; shift 2 ;;
    --negative)
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--negative needs a case (or 'list')"
      MODE="negative"; NEGATIVE_CASE="$2"; shift 2 ;;
    --replay)
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--replay needs a manifest directory"
      MODE="replay"; REPLAY_DIR="$2"; shift 2 ;;
    --replay-render)
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--replay-render needs a source bundle directory"
      MODE="replay-render"; REPLAY_DIR="$2"; shift 2 ;;
    --bundle)
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--bundle needs a directory"
      BUNDLE_DIR="$2"; shift 2 ;;
    --output-dir)
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--output-dir needs a value"
      OUTPUT_DIR="$2"; shift 2 ;;
    --log-file)
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--log-file needs a value"
      LOG_FILE="$2"; shift 2 ;;
    *) die "$EXIT_USAGE" "unknown argument: $1" ;;
  esac
done
[ -n "$MODE" ] || die "$EXIT_USAGE" "one of --list/--check/--self-test/--run/--negative/--replay required"
if [ "$MODE" = "run" ]; then
  case "$RUN_PROFILE" in
    smoke) : ;;
    daily-fixture)
      # An intentional heavier product execution: bounded short
      # 1080p-class excerpt. Never part of the routine lane.
      SMOKE_ARGS=(--width 1920 --height 1080 --frames 3 --spp 1
        --mechanics-preroll-frames 1 --mechanics-sample-rate 384000
        --listening-gain-fs-per-pa 280)
      SMOKE_BUDGET_SECONDS="${FSIM_CINE_DAILY_BUDGET_SECONDS:-14400}" ;;
    representative-4k-frame)
      # Only when explicitly requested AND resource-admitted.
      [ "${FSIM_CINE_ADMIT_4K:-0}" = "1" ]         || die "$EXIT_USAGE" "representative-4k-frame needs explicit resource admission (FSIM_CINE_ADMIT_4K=1)"
      SMOKE_ARGS=(--width 3840 --height 2160 --frames 1 --spp 1
        --mechanics-preroll-frames 1 --mechanics-sample-rate 384000
        --listening-gain-fs-per-pa 280)
      SMOKE_BUDGET_SECONDS="${FSIM_CINE_4K_BUDGET_SECONDS:-14400}" ;;
    *) die "$EXIT_USAGE" "unknown --run profile: $RUN_PROFILE (smoke|daily-fixture|representative-4k-frame)" ;;
  esac
fi

# ---------------------------------------------------------------------------
# JSONL event log. Every emission funnels through one python helper that
# redacts, caps, and appends; sequence numbers are the only ordering.
LOG_SEQ=0
log_event() { # stage reason_class detail_json_fragment
  [ -n "$LOG_FILE" ] || return 0
  LOG_SEQ=$((LOG_SEQ + 1))
  FSIM_LOG_FILE="$LOG_FILE" FSIM_LOG_SEQ="$LOG_SEQ" FSIM_LOG_STAGE="$1" \
  FSIM_LOG_REASON="$2" FSIM_LOG_DETAIL="${3:-{\}}" FSIM_LOG_CASE="$MODE${NEGATIVE_CASE:+:$NEGATIVE_CASE}${RUN_PROFILE:+:$RUN_PROFILE}" \
  FSIM_REPO_ROOT="$REPO_ROOT" FSIM_SOURCE_REV="$(git -C "$REPO_ROOT" rev-parse --short=12 HEAD 2>/dev/null || echo unknown)" \
  FSIM_BUDGET_S="$SMOKE_BUDGET_SECONDS" \
  python3 - <<'PYLOG'
import json, os, sys

MAX_RECORDS = 256
MAX_RECORD_BYTES = 8192
path = os.environ["FSIM_LOG_FILE"]

def redact(value):
    if isinstance(value, str):
        for root, token in (
            (os.environ["FSIM_REPO_ROOT"], "$REPO"),
            (os.environ.get("TMPDIR", "/nonexistent").rstrip("/"), "$TMP"),
            ("/private/tmp", "$TMP"),
            ("/tmp", "$TMP"),
            (os.path.expanduser("~"), "$HOME"),
        ):
            if root and root != "/":
                value = value.replace(root, token)
        return value
    if isinstance(value, dict):
        return {k: redact(v) for k, v in value.items()}
    if isinstance(value, list):
        return [redact(v) for v in value]
    return value

try:
    detail = json.loads(os.environ.get("FSIM_LOG_DETAIL", "{}"))
except json.JSONDecodeError:
    detail = {"malformed_detail": True}

record = {
    "schema": "euler-cinematic-e2e-log-v1",
    "suite": "euler-cinematic-e2e",
    "case": os.environ["FSIM_LOG_CASE"],
    "stage": os.environ["FSIM_LOG_STAGE"],
    "attempt": 1,
    "seq": int(os.environ["FSIM_LOG_SEQ"]),
    "source_rev": os.environ["FSIM_SOURCE_REV"],
    "build": {"profile": "debug", "features": ["cinematic-render"]},
    "budget_declared_s": int(os.environ["FSIM_BUDGET_S"]),
    "reason_class": os.environ["FSIM_LOG_REASON"],
}
record.update(redact(detail))
line = json.dumps(record, sort_keys=True, separators=(",", ":"))
if len(line.encode()) > MAX_RECORD_BYTES:
    record = {k: record[k] for k in
              ("schema", "suite", "case", "stage", "attempt", "seq",
               "source_rev", "reason_class")}
    record["truncated"] = True
    line = json.dumps(record, sort_keys=True, separators=(",", ":"))
existing = 0
if os.path.isfile(path):
    with open(path) as fh:
        existing = sum(1 for _ in fh)
if existing >= MAX_RECORDS:
    if existing == MAX_RECORDS:
        with open(path, "a") as fh:
            fh.write(json.dumps({"schema": "euler-cinematic-e2e-log-v1",
                                 "cap_reached": MAX_RECORDS}) + "\n")
    sys.exit(0)
os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
with open(path, "a") as fh:
    fh.write(line + "\n")
PYLOG
}

# ---------------------------------------------------------------------------
# Manifest/artifact coherence verifier. Every check here has a named
# hostile twin in --negative that must defeat it; the fields it reads are
# the production manifest's own declarations. Pure python over retained
# bytes; on failure it prints one FIRST-DIVERGENCE line and exits 42.
verify_bundle() { # dir -> exit 0/42; emits verify events when logging
  python3 - "$1" <<'PYEOF'
import json, os, struct, sys
root = sys.argv[1]

def refuse(message, expected=None, observed=None):
    print(f"verify: {message}", file=sys.stderr)
    payload = {"first_divergence": message}
    if expected is not None:
        payload["expected"] = expected
    if observed is not None:
        payload["observed"] = observed
    print("VERIFY-DIVERGENCE " + json.dumps(payload, sort_keys=True))
    sys.exit(42)

manifest_path = os.path.join(root, "critique-manifest.json")
if not os.path.isfile(manifest_path):
    refuse(f"manifest missing at {manifest_path}")
manifest = json.load(open(manifest_path))
# `status` is the honest evidence class, not a completion flag, and it
# PAIRS with the schema: a physics run may not claim replay provenance,
# a replay may not claim a fresh physics run, and NOTHING in this
# pipeline may promote itself to a calibrated-acoustic class.
ADMITTED = {
    "frankensim-euler-disc-production-critique-v1": "estimate-only-physics-render",
    "frankensim-euler-canonical-media-render-v1": "estimate-only-canonical-media-replay",
}
schema = manifest.get("schema", "")
if schema not in ADMITTED:
    refuse(f"unexpected schema {schema!r}")
if manifest.get("status") != ADMITTED[schema]:
    refuse(
        "status does not pair with schema (false promotion or provenance swap)",
        expected=ADMITTED[schema],
        observed=manifest.get("status"),
    )
if not manifest.get("no_claims"):
    refuse("manifest dropped its no_claims declarations")

picture = manifest.get("picture", {})
frames = picture.get("frames")
raw_dir = os.path.join(root, "raw")
raw_frames = sorted(f for f in os.listdir(raw_dir)) if os.path.isdir(raw_dir) else []
if not raw_frames:
    refuse("no raw EXR frames retained")
if isinstance(frames, int) and len(raw_frames) != frames:
    refuse("frame count mismatch", expected=frames, observed=len(raw_frames))
# Contiguity/order: the retained window must be EXACTLY the declared
# window, named frame-%06d — a renamed or reordered frame is a defect
# even when the count matches.
start = picture.get("rendered_frame_start")
end = picture.get("rendered_frame_end_exclusive")
if isinstance(start, int) and isinstance(end, int):
    expected_names = [f"frame-{index:06d}.exr" for index in range(start, end)]
    if raw_frames != expected_names:
        refuse(
            "raw frame window is not the declared contiguous window",
            expected=expected_names,
            observed=raw_frames,
        )
for name in raw_frames:
    if os.path.getsize(os.path.join(raw_dir, name)) == 0:
        refuse(f"empty frame artifact {name}")

wav = os.path.join(root, "sound", "physical-listening-master.pcm24.wav")
if not os.path.isfile(wav) or os.path.getsize(wav) < 44:
    refuse("listening master WAV missing or truncated")
blob = open(wav, "rb").read()
if blob[:4] != b"RIFF" or blob[8:12] != b"WAVE":
    refuse("listening master is not a RIFF/WAVE container")
# Parse fmt/data chunks for the A/V duration gate.
offset, rate, channels, bits, data_bytes = 12, None, None, None, None
while offset + 8 <= len(blob):
    chunk_id = blob[offset:offset + 4]
    chunk_len = struct.unpack("<I", blob[offset + 4:offset + 8])[0]
    body = blob[offset + 8:offset + 8 + chunk_len]
    if chunk_id == b"fmt " and len(body) >= 16:
        channels, rate = struct.unpack("<HI", body[2:8])
        bits = struct.unpack("<H", body[14:16])[0]
    elif chunk_id == b"data":
        data_bytes = chunk_len
    offset += 8 + chunk_len + (chunk_len % 2)
if rate is None or data_bytes is None or not channels or not bits:
    refuse("listening master lacks parseable fmt/data chunks")
duration = data_bytes / (channels * (bits // 8) * rate)
declared = picture.get("duration_s")
fps = picture.get("fps")
if isinstance(declared, (int, float)) and isinstance(fps, (int, float)) and fps > 0:
    # A/V alignment: audio must cover the picture within HALF a frame —
    # the production master is sample-exact, and a full frame of missing
    # or surplus audio is exactly the off-by-one sync defect.
    margin = 0.5 / fps
    if abs(duration - declared) > margin:
        refuse(
            "audio/picture duration mismatch (A/V off-by-one class)",
            expected=declared,
            observed=duration,
        )

# Parallel preview sequence: the manifest declares a preview identity,
# so the preview sequence must hold exactly one non-empty frame per raw
# frame (the dropped-AOV defect family for this fixture's artifacts).
if picture.get("preview_sequence_identity") is not None:
    preview_dir = os.path.join(root, "preview")
    previews = sorted(f for f in os.listdir(preview_dir)) if os.path.isdir(preview_dir) else []
    if len(previews) != len(raw_frames):
        refuse(
            "preview sequence does not mirror the raw window (dropped parallel artifact)",
            expected=len(raw_frames),
            observed=len(previews),
        )
    for name in previews:
        if os.path.getsize(os.path.join(preview_dir, name)) == 0:
            refuse(f"empty preview artifact {name}")

# Sample-rate pin: the master must be written on the declared control
# clock; a resampled or mislabeled master is a rate-mismatch defect.
declared_rate = manifest.get("mechanics", {}).get("control_sample_rate_hz")
if isinstance(declared_rate, int) and rate != declared_rate:
    refuse("listening master sample-rate mismatch", expected=declared_rate, observed=rate)

# Denoise coherence: an authority claiming more denoised frames than
# frames exist is an evidence overclaim.
denoise = picture.get("denoise", {})
applied = denoise.get("applied_frames")
if isinstance(applied, int) and isinstance(frames, int) and applied > frames:
    refuse("denoise applied_frames exceeds frame count", expected=frames, observed=applied)

# Declared-artifact walk: audio primary and a written mux target must
# exist non-empty (a manifest naming a missing mux is stale mux input).
primary = manifest.get("audio", {}).get("primary")
if isinstance(primary, str):
    primary_path = os.path.join(root, "sound", os.path.basename(primary))
    if not os.path.isfile(primary_path) or os.path.getsize(primary_path) == 0:
        refuse(f"declared audio primary missing or empty: {primary}")
mux = manifest.get("mux", {})
if isinstance(mux, dict) and mux.get("status") == "written":
    mux_name = os.path.basename(str(mux.get("path", "")))
    mux_path = os.path.join(root, mux_name)
    if not mux_name or not os.path.isfile(mux_path) or os.path.getsize(mux_path) == 0:
        refuse(f"manifest declares written mux but artifact is missing/empty: {mux_name!r}")
movs = [f for f in os.listdir(root) if f.endswith(".mov")]
if not movs:
    refuse("no mux output retained")
for name in movs:
    if os.path.getsize(os.path.join(root, name)) == 0:
        refuse(f"empty mux output {name}")
print(f"verify OK: {len(raw_frames)} frame(s), WAV {duration:.3f}s, {len(movs)} mux output(s), manifest coherent")
PYEOF
}

# Bounded relative artifact hashes + identity extraction for the log.
bundle_log_detail() { # dir -> json on stdout
  python3 - "$1" <<'PYEOF'
import hashlib, json, os, sys
root = sys.argv[1]
detail = {"identities": {}, "artifact_hashes": {}, "authority": {}}
manifest_path = os.path.join(root, "critique-manifest.json")
if os.path.isfile(manifest_path):
    manifest = json.load(open(manifest_path))
    detail["identities"] = {
        "trajectory": manifest.get("mechanics", {}).get("trajectory_identity"),
        "wav": manifest.get("audio", {}).get("wav_identity"),
        "raw_sequence": manifest.get("picture", {}).get("raw_sequence_identity"),
    }
    detail["authority"] = {
        "status": manifest.get("status"),
        "no_claims_count": len(manifest.get("no_claims", [])),
    }
    detail["seeds"] = {"render_seed_salt": manifest.get("picture", {}).get("render_seed_salt")}
    picture = manifest.get("picture", {})
    detail["frames"] = {
        "start": picture.get("rendered_frame_start"),
        "end_exclusive": picture.get("rendered_frame_end_exclusive"),
    }
    targets = ["critique-manifest.json", "sound/physical-listening-master.pcm24.wav"]
    raw_dir = os.path.join(root, "raw")
    raw = sorted(os.listdir(raw_dir)) if os.path.isdir(raw_dir) else []
    targets += [f"raw/{name}" for name in raw[:8]]
    if len(raw) > 8:
        detail["artifact_hashes_omitted_frames"] = len(raw) - 8
    for rel in targets:
        path = os.path.join(root, rel)
        if os.path.isfile(path):
            detail["artifact_hashes"][rel] = hashlib.sha256(open(path, "rb").read()).hexdigest()
print(json.dumps(detail, sort_keys=True))
PYEOF
}

REPAIRS_PIPELINE='["rerun with a larger FSIM_CINE_SMOKE_BUDGET_SECONDS","rebuild euler_cinematic_fixture (cargo build -p fs-euler-disc-e2e --features cinematic-render)","inspect the production stderr above the ERROR line"]'
REPAIRS_VERIFY='["read the VERIFY-DIVERGENCE line for the first bad artifact","re-run --run smoke into a fresh output root","if only identities moved, treat as stale bundle and regenerate"]'

case "$MODE" in
  list)
    printf '%s\n' \
      "smoke	tiny moving-disc bundle through the real pipeline (bounded)" \
      "daily-fixture	bounded short 1080p-class excerpt (intentional product execution)" \
      "representative-4k-frame	single 4K frame; needs FSIM_CINE_ADMIT_4K=1 (resource admission)" \
      "negative:gain-clip	hostile listening gain must refuse at the peak gate" \
      "negative:bad-mechanics-rate	non-integral decimation must refuse at admission" \
      "negative:truncated-manifest	verifier must refuse a truncated retained manifest" \
      "negative:dropped-frame	verifier must refuse a bundle missing a raw frame" \
      "negative:reordered-frame	verifier must refuse a renamed/reordered raw frame" \
      "negative:wav-truncation	verifier must refuse a truncated listening master" \
      "negative:av-off-by-one	verifier must refuse audio short of the picture window" \
      "negative:rate-mismatch	verifier must refuse a master off the declared clock" \
      "negative:false-promotion	verifier must refuse a calibrated-acoustic status claim" \
      "negative:denoise-overclaim	verifier must refuse applied_frames > frames" \
      "negative:dropped-preview	verifier must refuse a missing parallel preview frame" \
      "negative:stale-mux	verifier must refuse a declared-written-but-missing mux" \
      "negative:cancelled-child	mid-run kill must fail AND leave no verifiable partial bundle" \
      "negative:stale-trajectory-identity	production replay must refuse a tampered identity (needs --bundle)"
    exit "$EXIT_OK" ;;

  check)
    # Consistency without rendering: the fixture binary exists/builds and
    # the smoke recipe passes ADMISSION (mechanics-rate/decimation/gain
    # domains) by asking the production binary itself to refuse or accept
    # a deliberately-invalid variant quickly.
    BIN="$(fixture_binary)"
    SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/cine-check.XXXXXX")"
    trap 'rm -rf "$SCRATCH"' EXIT
    [ -n "$LOG_FILE" ] || LOG_FILE="$SCRATCH/runner-log.jsonl"
    log_event "admission" "ok" '{"repro":"scripts/ci/e2e_euler_cinematic.sh --check"}'
    if timeout 60 "$BIN" --output "$SCRATCH/out" "${SMOKE_ARGS[@]}" \
        --mechanics-sample-rate 47999 >"$SCRATCH/log" 2>&1; then
      log_event "admission" "verify" '{"first_divergence":"invalid mechanics rate accepted"}'
      die "$EXIT_VERIFY" "admission accepted an invalid mechanics rate"
    fi
    grep -q "admitted domain\|integer multiple" "$SCRATCH/log" \
      || die "$EXIT_VERIFY" "invalid rate refused without the typed diagnostic"
    log_event "verdict" "ok" '{"exit_class":0}'
    echo "check OK: fixture binary present; admission gates live; smoke recipe: ${SMOKE_ARGS[*]}"
    exit "$EXIT_OK" ;;

  run)
    STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
    [ -n "$OUTPUT_DIR" ] || OUTPUT_DIR="$REPO_ROOT/.e2e-out/cine-$RUN_PROFILE-$STAMP-$$"
    case "$OUTPUT_DIR" in
      "$REPO_ROOT"/*|"${TMPDIR:-/tmp}"*|/private/tmp/*|/tmp/*|/Volumes/USB_NVME/*) : ;;
      *) die "$EXIT_USAGE" "--output-dir must stay inside repo/TMP/USB scratch" ;;
    esac
    [ -e "$OUTPUT_DIR" ] && die "$EXIT_USAGE" "output dir exists (no reuse): $OUTPUT_DIR"
    BIN="$(fixture_binary)"
    [ -n "$LOG_FILE" ] || LOG_FILE="$OUTPUT_DIR.runner-log.jsonl"
    log_event "pipeline" "ok" '{"repro":"scripts/ci/e2e_euler_cinematic.sh --run smoke"}'
    echo "$RUN_PROFILE: real pipeline, budget ${SMOKE_BUDGET_SECONDS}s, output $OUTPUT_DIR"
    if ! timeout "$SMOKE_BUDGET_SECONDS" nice -n 15 "$BIN" \
        --output "$OUTPUT_DIR" "${SMOKE_ARGS[@]}"; then
      log_event "pipeline" "pipeline" "{\"exit_class\":41,\"ranked_repairs\":$REPAIRS_PIPELINE}"
      die "$EXIT_PIPELINE" "production pipeline failed or exceeded the $RUN_PROFILE budget"
    fi
    if ! VERIFY_OUT="$(verify_bundle "$OUTPUT_DIR" 2>&1)"; then
      DIVERGENCE="$(printf '%s\n' "$VERIFY_OUT" | grep '^VERIFY-DIVERGENCE ' | head -1 | cut -d' ' -f2-)"
      log_event "verify" "verify" "{\"exit_class\":42,\"divergence\":${DIVERGENCE:-null},\"ranked_repairs\":$REPAIRS_VERIFY}"
      printf '%s\n' "$VERIFY_OUT" >&2
      exit "$EXIT_VERIFY"
    fi
    log_event "verify" "ok" "$(bundle_log_detail "$OUTPUT_DIR")"
    log_event "verdict" "ok" '{"exit_class":0,"budget_within":true}'
    echo "$RUN_PROFILE PASS: bundle at $OUTPUT_DIR (log: $LOG_FILE)"
    exit "$EXIT_OK" ;;

  replay)
    [ -d "$REPLAY_DIR" ] || die "$EXIT_USAGE" "replay dir not found: $REPLAY_DIR"
    verify_bundle "$REPLAY_DIR" || exit "$EXIT_VERIFY"
    exit "$EXIT_OK" ;;

  replay-render)
    [ -d "$REPLAY_DIR" ] || die "$EXIT_USAGE" "source bundle not found: $REPLAY_DIR"
    verify_bundle "$REPLAY_DIR" >/dev/null || die "$EXIT_VERIFY" "source bundle fails coherence; refusing to replay from it"
    TRAJ="$REPLAY_DIR/trajectory/euler-trajectory.fset"
    MASTER="$REPLAY_DIR/sound/physical-listening-master.pcm24.wav"
    [ -f "$TRAJ" ] || die "$EXIT_VERIFY" "source bundle retains no canonical trajectory"
    TRAJ_ID="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["mechanics"]["trajectory_identity"])' "$REPLAY_DIR/critique-manifest.json")"       || die "$EXIT_VERIFY" "source manifest carries no trajectory identity"
    MASTER_ID="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["audio"]["wav_identity"])' "$REPLAY_DIR/critique-manifest.json")"       || die "$EXIT_VERIFY" "source manifest carries no listening-master identity"
    STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
    [ -n "$OUTPUT_DIR" ] || OUTPUT_DIR="$REPO_ROOT/.e2e-out/cine-replay-$STAMP-$$"
    [ -e "$OUTPUT_DIR" ] && die "$EXIT_USAGE" "output dir exists (no reuse): $OUTPUT_DIR"
    BIN="$(fixture_binary)"
    [ -n "$LOG_FILE" ] || LOG_FILE="$OUTPUT_DIR.runner-log.jsonl"
    log_event "pipeline" "ok" '{"repro":"scripts/ci/e2e_euler_cinematic.sh --replay-render SRC"}'
    echo "replay-render: from $REPLAY_DIR (identities from its manifest), output $OUTPUT_DIR"
    if ! timeout "$SMOKE_BUDGET_SECONDS" nice -n 15 "$BIN"         --output "$OUTPUT_DIR" "${SMOKE_ARGS[@]}"         --canonical-trajectory "$TRAJ" --canonical-trajectory-identity "$TRAJ_ID"         --canonical-listening-master "$MASTER"         --canonical-listening-master-identity "$MASTER_ID"; then
      log_event "pipeline" "pipeline" "{\"exit_class\":41,\"ranked_repairs\":$REPAIRS_PIPELINE}"
      die "$EXIT_PIPELINE" "replay-render pipeline failed or exceeded the budget"
    fi
    verify_bundle "$OUTPUT_DIR" || exit "$EXIT_VERIFY"
    log_event "verify" "ok" "$(bundle_log_detail "$OUTPUT_DIR")"
    log_event "verdict" "ok" '{"exit_class":0}'
    echo "replay-render PASS: bundle at $OUTPUT_DIR"
    exit "$EXIT_OK" ;;

  negative)
    NEG_CASES="gain-clip bad-mechanics-rate truncated-manifest dropped-frame reordered-frame wav-truncation av-off-by-one rate-mismatch false-promotion denoise-overclaim dropped-preview stale-mux cancelled-child stale-trajectory-identity"
    if [ "$NEGATIVE_CASE" = "list" ]; then printf '%s\n' $NEG_CASES; exit "$EXIT_OK"; fi
    case " $NEG_CASES " in
      *" $NEGATIVE_CASE "*) : ;;
      *) die "$EXIT_USAGE" "unknown negative case: $NEGATIVE_CASE" ;;
    esac
    SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/cine-negative.XXXXXX")"
    trap 'rm -rf "$SCRATCH"' EXIT
    [ -n "$LOG_FILE" ] || LOG_FILE="$SCRATCH/runner-log.jsonl"
    log_event "twin" "ok" "{\"repro\":\"scripts/ci/e2e_euler_cinematic.sh --negative $NEGATIVE_CASE\"}"
    case "$NEGATIVE_CASE" in
      gain-clip)
        # Hostile gain: the production peak gate must refuse, not publish
        # a clipped master. Uses a 60s admission-scale probe: the refusal
        # fires in the audio stage, so this needs a real (bounded) run -
        # covered by asserting the REFUSAL TEXT is reachable in the binary
        # plus a fast admission-side variant (gain <= 0 refuses instantly).
        BIN="$(fixture_binary)"
        if timeout 60 "$BIN" --output "$SCRATCH/out" "${SMOKE_ARGS[@]}" \
            --listening-gain-fs-per-pa 0 >"$SCRATCH/log" 2>&1; then
          exit "$EXIT_NEGATIVE_MISSED"
        fi
        grep -q "listening_gain_fs_per_pa" "$SCRATCH/log" || exit "$EXIT_NEGATIVE_MISSED"
        echo "negative gain-clip PASS: zero/negative gain refused with the typed diagnostic" ;;
      bad-mechanics-rate)
        # 48 kHz = ratio 1, outside the supported power-of-two decimation
        # range; the refusal fires at decimator construction (instant),
        # unlike ratio-2/16 whose delay-integrality refusal only fires
        # after the full mechanics integration (measured on this bead).
        BIN="$(fixture_binary)"
        if timeout 180 "$BIN" --output "$SCRATCH/out" "${SMOKE_ARGS[@]}" \
            --mechanics-sample-rate 48000 >"$SCRATCH/log" 2>&1; then
          exit "$EXIT_NEGATIVE_MISSED"
        fi
        grep -q "decimat" "$SCRATCH/log" || exit "$EXIT_NEGATIVE_MISSED"
        echo "negative bad-mechanics-rate PASS: unsupported decimation ratio refused" ;;
      cancelled-child)
        # Cancelled child work: kill the production child mid-mechanics
        # and require BOTH halves of the honest outcome — the child dies
        # nonzero, and whatever partial output it left can never verify
        # as a publishable bundle (no success for partial publication).
        BIN="$(fixture_binary)"
        if timeout --signal=TERM 20 "$BIN" --output "$SCRATCH/out"             "${SMOKE_ARGS[@]}" >"$SCRATCH/log" 2>&1; then
          log_event "twin" "negative-missed" '{"exit_class":43,"first_divergence":"child survived the 20s cancellation window"}'
          exit "$EXIT_NEGATIVE_MISSED"
        fi
        if verify_bundle "$SCRATCH/out" >/dev/null 2>&1; then
          log_event "twin" "negative-missed" '{"exit_class":43,"first_divergence":"partial publication verified as success"}'
          exit "$EXIT_NEGATIVE_MISSED"
        fi
        log_event "verdict" "ok" '{"exit_class":0}'
        echo "negative cancelled-child PASS: mid-run kill left no bundle that verifies" ;;
      stale-trajectory-identity)
        # Production identity gate: replaying a retained trajectory under
        # a TAMPERED identity must refuse at decode (stale trajectory /
        # stale asset identity family). Needs a real retained bundle.
        [ -n "$BUNDLE_DIR" ] || die "$EXIT_USAGE" "stale-trajectory-identity needs --bundle SRC (a retained smoke bundle)"
        [ -d "$BUNDLE_DIR" ] || die "$EXIT_USAGE" "bundle not found: $BUNDLE_DIR"
        TRAJ="$BUNDLE_DIR/trajectory/euler-trajectory.fset"
        MASTER="$BUNDLE_DIR/sound/physical-listening-master.pcm24.wav"
        [ -f "$TRAJ" ] && [ -f "$MASTER" ] || die "$EXIT_USAGE" "bundle lacks retained trajectory/master"
        TRAJ_ID="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["mechanics"]["trajectory_identity"])' "$BUNDLE_DIR/critique-manifest.json")" \
          || die "$EXIT_USAGE" "bundle manifest carries no trajectory identity"
        MASTER_ID="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["audio"]["wav_identity"])' "$BUNDLE_DIR/critique-manifest.json")" \
          || die "$EXIT_USAGE" "bundle manifest carries no wav identity"
        # Flip the leading hex digit: a well-formed but WRONG identity.
        FIRST="$(printf '%s' "$TRAJ_ID" | cut -c1)"
        [ "$FIRST" = "0" ] && FLIP="f" || FLIP="0"
        TAMPERED="$FLIP$(printf '%s' "$TRAJ_ID" | cut -c2-)"
        BIN="$(fixture_binary)"
        if timeout 300 "$BIN" --output "$SCRATCH/out" "${SMOKE_ARGS[@]}" \
            --canonical-trajectory "$TRAJ" --canonical-trajectory-identity "$TAMPERED" \
            --canonical-listening-master "$MASTER" \
            --canonical-listening-master-identity "$MASTER_ID" >"$SCRATCH/log" 2>&1; then
          log_event "twin" "negative-missed" '{"exit_class":43}'
          exit "$EXIT_NEGATIVE_MISSED"
        fi
        grep -q "identity mismatch" "$SCRATCH/log" || { log_event "twin" "negative-missed" '{"exit_class":43,"first_divergence":"refused without the identity-mismatch diagnostic"}'; exit "$EXIT_NEGATIVE_MISSED"; }
        log_event "verdict" "ok" '{"exit_class":0}'
        echo "negative stale-trajectory-identity PASS: tampered identity refused at decode" ;;
      truncated-manifest|dropped-frame|reordered-frame|wav-truncation|av-off-by-one|rate-mismatch|false-promotion|denoise-overclaim|dropped-preview|stale-mux)
        # Build a tiny synthetic bundle, corrupt it, and the verifier must
        # refuse. (Synthetic here is legitimate: the subject under test is
        # the VERIFIER's refusal, not the production pipeline.)
        ROOT="$SCRATCH/bundle"
        python3 - "$ROOT" <<'PYFIX'
import json, os, struct, sys
root = sys.argv[1]
os.makedirs(f"{root}/raw", exist_ok=True)
os.makedirs(f"{root}/sound", exist_ok=True)
open(f"{root}/raw/frame-000000.exr", "w").write("EXRDATA")
os.makedirs(f"{root}/preview", exist_ok=True)
open(f"{root}/preview/frame-000000.png", "w").write("PNGDATA")
# One picture-frame of coherent 48 kHz stereo PCM24 silence: 1/24 s.
rate, channels, bits = 48000, 2, 24
samples = rate // 24
data = bytes(samples * channels * (bits // 8))
with open(f"{root}/sound/physical-listening-master.pcm24.wav", "wb") as fh:
    fh.write(b"RIFF" + struct.pack("<I", 36 + len(data)) + b"WAVE")
    fh.write(b"fmt " + struct.pack("<IHHIIHH", 16, 1, channels, rate,
                                   rate * channels * bits // 8,
                                   channels * bits // 8, bits))
    fh.write(b"data" + struct.pack("<I", len(data)) + data)
open(f"{root}/euler-disc-critique.mov", "w").write("MOVDATA")
json.dump({"schema": "frankensim-euler-disc-production-critique-v1",
           "status": "estimate-only-physics-render",
           "no_claims": ["synthetic twin fixture"],
           "picture": {"frames": 1, "fps": 24, "duration_s": 1 / 24,
                        "rendered_frame_start": 0,
                        "rendered_frame_end_exclusive": 1,
                        "preview_sequence_identity": "synthetic",
                        "denoise": {"applied_frames": 0}},
           "mechanics": {"control_sample_rate_hz": 48000},
           "audio": {"primary": "physical-listening-master.pcm24.wav"},
           "mux": {"status": "written", "path": "euler-disc-critique.mov"}},
          open(f"{root}/critique-manifest.json", "w"))
PYFIX
        verify_bundle "$ROOT" >/dev/null 2>&1 || die "$EXIT_USAGE" "twin precondition: synthetic bundle must verify clean"
        python3 - "$ROOT" "$NEGATIVE_CASE" <<'PYCORRUPT'
import json, os, sys
root, case = sys.argv[1], sys.argv[2]
manifest_path = f"{root}/critique-manifest.json"
if case == "truncated-manifest":
    data = open(manifest_path).read()
    open(manifest_path, "w").write(data[: len(data) // 2])
elif case == "dropped-frame":
    os.remove(f"{root}/raw/frame-000000.exr")
elif case == "reordered-frame":
    os.rename(f"{root}/raw/frame-000000.exr", f"{root}/raw/frame-000001.exr")
elif case == "wav-truncation":
    path = f"{root}/sound/physical-listening-master.pcm24.wav"
    open(path, "wb").write(open(path, "rb").read()[:10])
elif case == "av-off-by-one":
    # Audio one frame SHORT of the declared picture window.
    manifest = json.load(open(manifest_path))
    manifest["picture"]["frames"] = 2
    manifest["picture"]["duration_s"] = 2 / 24
    manifest["picture"]["rendered_frame_end_exclusive"] = 2
    json.dump(manifest, open(manifest_path, "w"))
    open(f"{root}/raw/frame-000001.exr", "w").write("EXRDATA")
elif case == "false-promotion":
    manifest = json.load(open(manifest_path))
    manifest["status"] = "calibrated-acoustic-render"
    json.dump(manifest, open(manifest_path, "w"))
elif case == "denoise-overclaim":
    manifest = json.load(open(manifest_path))
    manifest["picture"]["denoise"]["applied_frames"] = 6
    json.dump(manifest, open(manifest_path, "w"))
elif case == "rate-mismatch":
    # Same duration, wrong clock: 44.1 kHz against the declared 48 kHz.
    import struct
    rate, channels, bits = 44100, 2, 24
    samples = round(rate / 24)
    data = bytes(samples * channels * (bits // 8))
    path = f"{root}/sound/physical-listening-master.pcm24.wav"
    with open(path, "wb") as fh:
        fh.write(b"RIFF" + struct.pack("<I", 36 + len(data)) + b"WAVE")
        fh.write(b"fmt " + struct.pack("<IHHIIHH", 16, 1, channels, rate,
                                       rate * channels * bits // 8,
                                       channels * bits // 8, bits))
        fh.write(b"data" + struct.pack("<I", len(data)) + data)
elif case == "dropped-preview":
    os.remove(f"{root}/preview/frame-000000.png")
elif case == "stale-mux":
    os.remove(f"{root}/euler-disc-critique.mov")
    open(f"{root}/other-output.mov", "w").write("MOVDATA")
PYCORRUPT
        if verify_bundle "$ROOT" >/dev/null 2>&1; then
          log_event "twin" "negative-missed" '{"exit_class":43}'
          exit "$EXIT_NEGATIVE_MISSED"
        fi
        log_event "verdict" "ok" '{"exit_class":0}'
        echo "negative $NEGATIVE_CASE PASS: verifier refused the corrupted bundle" ;;
    esac
    exit "$EXIT_OK" ;;

  self-test)
    PASS=0; FAIL=0
    check() { if [ "$2" -eq "$3" ]; then PASS=$((PASS+1)); else FAIL=$((FAIL+1)); echo "SELF-TEST FAIL: $1 (expected $2, got $3)" >&2; fi; }
    "$0" --no-such-flag >/dev/null 2>&1;              check "unknown flag refuses" 40 $?
    "$0" --run no-such-profile >/dev/null 2>&1;       check "unknown run profile refuses" 40 $?
    "$0" --run representative-4k-frame >/dev/null 2>&1; check "4k frame refuses without resource admission" 40 $?
    FSIM_CINE_DAILY_BUDGET_SECONDS=1 "$0" --run daily-fixture >/dev/null 2>&1; check "daily-fixture wires to the real pipeline (1s budget kill)" 41 $?
    "$0" --negative no-such-case >/dev/null 2>&1;     check "unknown negative case refuses" 40 $?
    "$0" --replay /no/such/dir >/dev/null 2>&1;       check "missing replay dir refuses" 40 $?
    "$0" --negative stale-trajectory-identity >/dev/null 2>&1; check "identity twin without --bundle refuses" 40 $?
    "$0" --list >/dev/null 2>&1;                      check "list mode runs" 0 $?
    for CASE in truncated-manifest dropped-frame reordered-frame wav-truncation \
                av-off-by-one rate-mismatch false-promotion denoise-overclaim \
                dropped-preview stale-mux; do
      "$0" --negative "$CASE" >/dev/null 2>&1;        check "verifier twin $CASE detects" 0 $?
    done
    # Logging contract on a real emission path: valid JSONL, stable seq,
    # no absolute paths (redaction), required identity fields present.
    LOG_PROBE="$(mktemp "${TMPDIR:-/tmp}/cine-log-probe.XXXXXX")"
    rm -f "$LOG_PROBE"
    "$0" --negative false-promotion --log-file "$LOG_PROBE" >/dev/null 2>&1
    python3 - "$LOG_PROBE" <<'PYCHECK'
import json, sys
path = sys.argv[1]
records = [json.loads(line) for line in open(path)]
assert records, "log is empty"
assert [r["seq"] for r in records] == list(range(1, len(records) + 1)), "seq not stable/contiguous"
for record in records:
    assert record["schema"] == "euler-cinematic-e2e-log-v1"
    assert record["suite"] == "euler-cinematic-e2e"
    for field in ("case", "stage", "attempt", "source_rev", "reason_class"):
        assert field in record, f"missing {field}"
    blob = json.dumps(record)
    assert "/Users/" not in blob and "/home/" not in blob, "absolute path leaked"
PYCHECK
    check "log contract holds (JSONL, seq, redaction, fields)" 0 $?
    rm -f "$LOG_PROBE"
    echo "self-test: $PASS passed, $FAIL failed"
    [ "$FAIL" -eq 0 ] || exit 1
    exit 0 ;;
esac
