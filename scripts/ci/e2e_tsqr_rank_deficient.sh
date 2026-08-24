#!/usr/bin/env bash
# Rank-deficient TSQR release E2E runner (bead frankensim-epic-bedrock-6ys.5.1.6).
#
# Stable executable surface for rank-deficient TSQR verification:
# Fixed-tree reduction, gauge-invariant decomposition, rank-deficient QR certificates,
# adversarial tree gauge adjudication, and honest no-claim boundaries.
#
# Modes:
#   --list            enumerate cases and twins, run nothing
#   --check           tool/crate consistency, no heavy runs
#   --self-test       runner's own failure-detection on fixtures
#   --run smoke       bounded smoke battery (tsqr_rank_deficient + tsqr_tree_gauge)
#   --run full        full release battery incl. policy/checker/driver lanes
#   --negative CASE   one named hostile twin executed for real (or 'list')
#   --replay FILE     verify retained TSQR certificate artifact(s)
#   --output-dir DIR  artifact root
#
# EXIT CLASSES: 0 ok; 40 usage; 41 pipeline failure; 42 verification failure;
#               43 negative twin NOT detected (falsifier law reachable).

set -u -o pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
readonly REPO_ROOT
readonly CARGO_BIN="${CARGO_BIN:-cargo}"
PROBE="tsqr_e2e_probe"

EXIT_OK=0; EXIT_USAGE=40; EXIT_PIPELINE=41; EXIT_VERIFY=42; EXIT_NEG_MISSED=43

OUT_DIR="${TMPDIR:-/tmp}/tsqr-e2e-$$"
LOG_FILE=""
SEQ=0
RECORDS=0

log() {
  SEQ=$((SEQ + 1))
  [ "$RECORDS" -ge 256 ] && return 0
  RECORDS=$((RECORDS + 1))
  local line="{\"suite\":\"tsqr-e2e\",\"seq\":$SEQ,\"stage\":\"$1\",$2}"
  line="${line//$REPO_ROOT/<repo>}"
  line="${line//${TMPDIR:-/tmp}/<tmp>}"
  line="${line//$HOME/<home>}"
  line="${line:0:4096}"
  echo "$line"
  [ -n "$LOG_FILE" ] && echo "$line" >> "$LOG_FILE"
  return 0
}

die() {
  log "error" "\"class\":$1,\"message\":\"$2\",\"repro\":\"$3\""
  exit "$1"
}

need() { command -v "$1" >/dev/null 2>&1 || die "$EXIT_USAGE" "$1 required" "install $1"; }

TWINS=(
  "raw-factor-equality:refuses raw-factor bitwise equality claims on rank-deficient cross-trees"
  "tampered-certificate:refuses invalid QR certificate with corrupted R matrix"
  "unsupported-tree-claim:refuses moonshot arbitrary-tree gauge promotion without proof"
  "scale-blind-tolerance:refuses absolute rank tolerance under dimensional scaling"
)

probe() {
  # Build-once semantics come from the shared cargo target dir; release
  # profile exercises the real optimized pipeline.
  "$CARGO_BIN" run --release -p fs-la --example "$PROBE" -- "$@" \
    || die "$EXIT_NEG_MISSED" "probe guard violated: $*" \
         "cargo run --release -p fs-la --example $PROBE -- $*"
}

list_cases() {
  log "list" "\"modes\":[\"check\",\"self-test\",\"run smoke\",\"run full\",\"negative\",\"replay\"]"
  for t in "${TWINS[@]}"; do
    local name="${t%%:*}"
    local desc="${t#*:}"
    log "twin" "\"name\":\"$name\",\"description\":\"$desc\""
  done
}

run_check() {
  log "check" "\"status\":\"checking toolchain and crates\""
  need "$CARGO_BIN"
  [ -f "$REPO_ROOT/crates/fs-la/src/canonical_qr.rs" ] || die "$EXIT_USAGE" "canonical_qr surface missing" "git ls-files crates/fs-la"
  [ -f "$REPO_ROOT/crates/fs-la/src/canonical_tree.rs" ] || die "$EXIT_USAGE" "canonical_tree driver missing" "git ls-files crates/fs-la"
  [ -f "$REPO_ROOT/crates/fs-la/src/canonical_check.rs" ] || die "$EXIT_USAGE" "canonical_check checker missing" "git ls-files crates/fs-la"
  [ -x "$REPO_ROOT/scripts/ci/e2e_tsqr_rank_deficient.sh" ] || chmod +x "$REPO_ROOT/scripts/ci/e2e_tsqr_rank_deficient.sh"
  log "check" "\"status\":\"ok\",\"cargo\":\"$CARGO_BIN\""
}

run_self_test() {
  log "self-test" "\"stage\":\"starting runner self-test\""
  local count=${#TWINS[@]}
  if [ "$count" -lt 4 ]; then
    die "$EXIT_VERIFY" "insufficient hostile twins registered" "./scripts/ci/e2e_tsqr_rank_deficient.sh --self-test"
  fi
  # The self-test proves the falsifier path is LIVE: a deliberately broken
  # tolerance must be refused by the typed surface (exit 0 expected here
  # because scale-blind twin catches it; any other exit means detection broke).
  probe twin-scale-blind-tolerance \
    || die "$EXIT_VERIFY" "self-test: scale-blind twin failed to detect" "./scripts/ci/e2e_tsqr_rank_deficient.sh --negative scale-blind-tolerance"
  log "self-test" "\"status\":\"pass\",\"registered_twins\":$count,\"falsifier_path\":\"live\""
}

emit_artifacts() {
  mkdir -p "$OUT_DIR"
  local cert_file="$OUT_DIR/certificates.jsonl"
  : > "$cert_file"
  log "artifacts" "\"stage\":\"emitting certificates\",\"path\":\"$cert_file\""
  probe emit-certificate-full-rank  >> "$cert_file" || die "$EXIT_PIPELINE" "full-rank certificate emission failed" "probe emit-certificate-full-rank"
  probe emit-certificate-deficient >> "$cert_file" || die "$EXIT_PIPELINE" "deficient certificate emission failed" "probe emit-certificate-deficient"
  printf '%s\n' "$cert_file"
}

run_smoke() {
  log "smoke" "\"stage\":\"running tsqr_rank_deficient battery\""
  "$CARGO_BIN" test -p fs-la --test tsqr_rank_deficient --release -- --nocapture \
    || die "$EXIT_PIPELINE" "tsqr_rank_deficient failed" "cargo test -p fs-la --test tsqr_rank_deficient --release"

  log "smoke" "\"stage\":\"running tsqr_tree_gauge battery\""
  "$CARGO_BIN" test -p fs-la --test tsqr_tree_gauge --release -- --nocapture \
    || die "$EXIT_PIPELINE" "tsqr_tree_gauge failed" "cargo test -p fs-la --test tsqr_tree_gauge --release"

  local cert_file
  cert_file="$(emit_artifacts)"
  log "replay" "\"stage\":\"verifying freshly emitted artifacts\",\"file\":\"<artifacts>/certificates.jsonl\""
  probe verify-certificate "$cert_file" || die "$EXIT_VERIFY" "artifact verification failed" "probe verify-certificate <certificates.jsonl>"

  log "smoke" "\"status\":\"pass\""
}

run_full() {
  log "full" "\"stage\":\"policy/result surface lane\""
  "$CARGO_BIN" test -p fs-la --lib canonical_qr --release -- --nocapture \
    || die "$EXIT_PIPELINE" "canonical_qr unit lane failed" "cargo test -p fs-la --lib canonical_qr --release"
  "$CARGO_BIN" test -p fs-la --test canonical_qr_policy --release -- --nocapture \
    || die "$EXIT_PIPELINE" "canonical_qr_policy failed" "cargo test -p fs-la --test canonical_qr_policy --release"

  log "full" "\"stage\":\"driver cancellation/resume/fork lane\""
  "$CARGO_BIN" test -p fs-la --lib canonical_tree --release --nocapture -- --nocapture \
    || die "$EXIT_PIPELINE" "canonical_tree lane failed" "cargo test -p fs-la --lib canonical_tree --release"

  log "full" "\"stage\":\"independent checker lane\""
  "$CARGO_BIN" test -p fs-la --lib canonical_check --release -- --nocapture \
    || die "$EXIT_PIPELINE" "canonical_check lane failed" "cargo test -p fs-la --lib canonical_check --release"

  run_smoke

  log "full" "\"stage\":\"executing all hostile twins for real\""
  for t in "${TWINS[@]}"; do
    local name="${t%%:*}"
    probe "twin-$name" || die "$EXIT_NEG_MISSED" "twin not detected: $name" "./scripts/ci/e2e_tsqr_rank_deficient.sh --negative $name"
  done

  log "full" "\"status\":\"pass\""
}

run_negative() {
  local target="$1"
  log "negative" "\"case\":\"$target\",\"stage\":\"executing hostile twin for real\""
  case "$target" in
    list) list_cases ;;
    raw-factor-equality|tampered-certificate|unsupported-tree-claim|scale-blind-tolerance)
      # REAL execution: the probe exits nonzero iff the guarantee regressed.
      local before=$?
      probe "twin-$target"
      local rc=$?
      [ "$rc" -eq 0 ] || die "$EXIT_NEG_MISSED" "twin $target NOT detected (rc=$rc)" "./scripts/ci/e2e_tsqr_rank_deficient.sh --negative $target"
      log "negative" "\"case\":\"$target\",\"outcome\":\"caught\",\"exit_before\":$before"
      ;;
    *)
      die "$EXIT_USAGE" "unknown negative case: $target" "./scripts/ci/e2e_tsqr_rank_deficient.sh --negative list"
      ;;
  esac
}

run_replay() {
  local file="$1"
  log "replay" "\"file\":\"<artifact>\",\"stage\":\"verifying retained certificate artifact\""
  [ -f "$file" ] || die "$EXIT_USAGE" "replay file not found: $file" "./scripts/ci/e2e_tsqr_rank_deficient.sh --replay <file>"
  probe verify-certificate "$file" || die "$EXIT_VERIFY" "artifact verification failed" "./scripts/ci/e2e_tsqr_rank_deficient.sh --replay $file"
  log "replay" "\"status\":\"verified\""
}

# --- CLI Dispatch ---
if [ $# -eq 0 ]; then
  echo "usage: $0 --list | --check | --self-test | --run smoke | --run full | --negative CASE | --replay FILE [--output-dir DIR]" >&2
  exit "$EXIT_USAGE"
fi

MODE="$1"
shift || true

case "$MODE" in
  --list)       list_cases ;;
  --check)      run_check ;;
  --self-test)  run_self_test ;;
  --run)
    case "${1:-}" in
      smoke) OUT_DIR="${OUT_DIR_OVERRIDE:-$OUT_DIR}"; LOG_FILE="$OUT_DIR/runner.jsonl"; mkdir -p "$OUT_DIR"; run_smoke ;;
      full)  OUT_DIR="${OUT_DIR_OVERRIDE:-$OUT_DIR}"; LOG_FILE="$OUT_DIR/runner.jsonl"; mkdir -p "$OUT_DIR"; run_full ;;
      *)     die "$EXIT_USAGE" "--run requires smoke|full" "$0 --run smoke" ;;
    esac
    ;;
  --negative)
    case "${1:-}" in
      "") die "$EXIT_USAGE" "--negative requires CASE or list" "$0 --negative list" ;;
      *)  run_negative "$1" ;;
    esac
    ;;
  --replay)
    case "${1:-}" in
      "") die "$EXIT_USAGE" "--replay requires FILE" "$0 --replay <file>" ;;
      *)  run_replay "$1" ;;
    esac
    ;;
  --output-dir)
    OUT_DIR="${2:-$OUT_DIR}"
    shift
    ;;
  *)
    die "$EXIT_USAGE" "unknown mode: $MODE" "$0 --help"
    ;;
esac

exit "$EXIT_OK"
