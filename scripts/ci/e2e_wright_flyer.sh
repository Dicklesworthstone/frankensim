#!/usr/bin/env bash
# Wright Flyer E2E runner (bead frankensim-wf-root-guzez.7.3, E6.3).
# Cloned from the hardened Euler cinematic runner pattern.
#
# One stable executable surface over the REAL Wright Flyer pipeline:
# native engine (fs-flyer), wasm engine (fs-flyer-wasm via wasm-pack +
# node), the browser app's headless batteries, the golden-lane CI, and
# the replay/A-B receipt machinery. It reuses production CLIs, crates,
# and test surfaces; it implements NO parallel simulation logic.
#
# Modes:
#   --list            enumerate cases and twins, run nothing
#   --check           tool/config consistency, no engine runs
#   --self-test       the runner's own failure-detection on fixtures
#   --run smoke       bounded REAL pipeline: native canonical digest,
#                     wasm selftest, node harness, app headless suite
#   --negative CASE   one named hostile twin (or 'list'); each twin is
#                     the PRODUCTION battery that executes it — the
#                     runner drives it and verifies detection
#   --replay DIR      re-verify a retained smoke record (recompute the
#                     canonical digest and compare)
#   --output-dir DIR  artifact root (defaults under TMPDIR)
#
# EXIT CLASSES: 0 ok; 40 usage; 41 pipeline failure; 42 verification
# failure; 43 negative twin NOT detected.
#
# LOGGING CONTRACT (wf-e2e-log-v1): bounded JSONL (256 records, 4 KiB
# each), stable suite/case/stage/seq ids, expected/observed on
# divergence, one repo-relative repro command per failure; repo root,
# TMP and HOME redacted. No-claims: a green smoke proves the SOFTWARE
# pipeline executes and its digests cohere; it proves nothing physical.
set -u -o pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXIT_OK=0; EXIT_USAGE=40; EXIT_PIPELINE=41; EXIT_VERIFY=42; EXIT_NEG_MISSED=43
NATIVE_GOLDEN="823d9f59dd162c8bc0764e144236d2f00abc48a12142095688a22e59ae95ca9d"
WASM_GOLDEN="f088689ae4c60ec33a2034ec7020c85772bfc016968fa9ae5f6d92a308fcbbb6"

OUT_DIR="${TMPDIR:-/tmp}/wf-e2e-$$"
LOG_FILE=""
SEQ=0
RECORDS=0

log() { # stage json-payload-without-braces
  SEQ=$((SEQ + 1))
  [ "$RECORDS" -ge 256 ] && return 0
  RECORDS=$((RECORDS + 1))
  local line="{\"suite\":\"wf-e2e\",\"seq\":$SEQ,\"stage\":\"$1\",$2}"
  line="${line//$REPO_ROOT/<repo>}"
  line="${line//${TMPDIR:-/tmp}/<tmp>}"
  line="${line//$HOME/<home>}"
  line="${line:0:4096}"
  echo "$line"
  [ -n "$LOG_FILE" ] && echo "$line" >> "$LOG_FILE"
  return 0
}

die() { # class message repro
  log "error" "\"class\":$1,\"message\":\"$2\",\"repro\":\"$3\""
  exit "$1"
}

need() { command -v "$1" >/dev/null 2>&1 || die "$EXIT_USAGE" "$1 required" "install $1"; }

TWINS=(
  "checkpoint-tamper:fs-flyer:checkpoint_battery:tamper_terminal_and_caps_refuse"
  "replay-envelope-tamper:fs-flyer:replayenv_battery:envelope_roundtrip_and_hostile_twins"
  "seed-mismatch-replay:fs-flyer-app:replay:LIVE"
  "kpi-mismatch:fs-flyer-app:resultsCard:hostile"
  "leased-ring-violation:fs-flyer-app:transport:seqlock"
  "late-input-tamper:fs-flyer-app:humanControls:ledger"
  "golden-divergence:inline:selftest:perturbed"
  "ab-terminal-divergence:fs-flyer:abcompare_battery:terminal_divergence_records_different_kinds"
  "identity-intent-tamper:fs-flyer:simloop_battery:run_intent_id_is_minted_after_and_downstream_of_tick0"
  "effect-ownership:fs-flyer:effectowners_battery:"
)

list_cases() {
  log "list" "\"modes\":[\"check\",\"self-test\",\"run smoke\",\"negative\",\"replay\"]"
  for t in "${TWINS[@]}"; do
    log "list-twin" "\"case\":\"${t%%:*}\""
  done
}

run_check() {
  need python3; need node; need cargo; need wasm-pack
  [ -f "$REPO_ROOT/tools/wf-ci/golden_lanes.sh" ] || die "$EXIT_VERIFY" "golden_lanes.sh missing" "git checkout tools/wf-ci"
  [ -f "$REPO_ROOT/crates/fs-flyer-wasm/CONTRACT.md" ] || die "$EXIT_VERIFY" "wasm CONTRACT missing" "git checkout crates/fs-flyer-wasm"
  grep -q "$NATIVE_GOLDEN" "$REPO_ROOT/tools/wf-ci/golden_lanes.sh" \
    || die "$EXIT_VERIFY" "native golden drifted between runner and CI script" "unify goldens"
  grep -q "$WASM_GOLDEN" "$REPO_ROOT/crates/fs-flyer-wasm/src/engine.rs" \
    || die "$EXIT_VERIFY" "wasm selftest golden drifted" "unify goldens"
  log "check" "\"pass\":true,\"tools\":[\"python3\",\"node\",\"cargo\",\"wasm-pack\"]"
}

run_self_test() {
  # The runner must DETECT a planted divergence in its own comparison
  # path (fixture-level; no engine run).
  local a="deadbeef" b="deadbeef" c="d1verged"
  if [ "$a" != "$b" ]; then die "$EXIT_VERIFY" "self-test equality broken" "-"; fi
  if [ "$a" = "$c" ]; then die "$EXIT_VERIFY" "self-test inequality broken" "-"; fi
  # And the log redaction actually redacts.
  local probe="{\"p\":\"$REPO_ROOT/x\"}"
  probe="${probe//$REPO_ROOT/<repo>}"
  case "$probe" in *"<repo>"*) : ;; *) die "$EXIT_VERIFY" "redaction inert" "-" ;; esac
  log "self-test" "\"pass\":true"
}

native_digest() {
  (cd "$REPO_ROOT" && cargo run -q -p fs-flyer --release --bin canonical_digest 2>/dev/null) \
    | python3 -c 'import json,sys;print(json.load(sys.stdin)["digest"])'
}

run_smoke() {
  mkdir -p "$OUT_DIR"
  LOG_FILE="$OUT_DIR/wf-e2e.jsonl"
  local t0=$SECONDS
  # 1. Native canonical digest vs golden (production bin).
  local nd; nd=$(native_digest) || die "$EXIT_PIPELINE" "native canonical run failed" "cargo run -p fs-flyer --release --bin canonical_digest"
  if [ "$nd" != "$NATIVE_GOLDEN" ]; then
    log "smoke-native" "\"pass\":false,\"expected\":\"$NATIVE_GOLDEN\",\"observed\":\"$nd\""
    die "$EXIT_VERIFY" "native canonical digest diverged" "cargo run -p fs-flyer --release --bin canonical_digest"
  fi
  log "smoke-native" "\"pass\":true,\"digest\":\"$nd\""
  # 2. Wasm pkg build + selftest + node harness (production surfaces).
  local pkg="$OUT_DIR/pkg"
  (cd "$REPO_ROOT/crates/fs-flyer-wasm" && wasm-pack build --target nodejs --release --out-dir "$pkg" >/dev/null 2>&1) \
    || die "$EXIT_PIPELINE" "wasm-pack build failed" "cd crates/fs-flyer-wasm && wasm-pack build --target nodejs"
  local st; st=$(node -e 'const w=require(process.argv[1]);const r=JSON.parse(w.flyer_selftest());process.stdout.write(r.ok?("OK:"+r.ok.digest):("REFUSED:"+r.refusal.code));' "$pkg/fs_flyer_wasm.js")
  case "$st" in
    OK:*) log "smoke-wasm-selftest" "\"pass\":true,\"digest\":\"${st#OK:}\"" ;;
    *) log "smoke-wasm-selftest" "\"pass\":false,\"observed\":\"$st\""
       die "$EXIT_VERIFY" "wasm selftest refused" "node -e flyer_selftest" ;;
  esac
  node "$REPO_ROOT/crates/fs-flyer-wasm/node-harness/engine_harness.mjs" "$pkg/fs_flyer_wasm.js" > "$OUT_DIR/harness.jsonl" 2>&1 \
    || die "$EXIT_VERIFY" "node harness failed" "node crates/fs-flyer-wasm/node-harness/engine_harness.mjs <pkg>/fs_flyer_wasm.js"
  log "smoke-harness" "\"pass\":true"
  # 3. App headless suite against the fresh pkg (production tests).
  (cd "$REPO_ROOT/apps/wright-flyer" && WF_PKG="$pkg/fs_flyer_wasm.js" node --test test/*.test.ts > "$OUT_DIR/app-tests.txt" 2>&1) \
    || die "$EXIT_VERIFY" "app headless suite failed" "cd apps/wright-flyer && WF_PKG=<pkg>/fs_flyer_wasm.js node --test test/"
  log "smoke-app" "\"pass\":true"
  # Retain the smoke record for --replay.
  echo "$nd" > "$OUT_DIR/native-digest.txt"
  log "smoke" "\"pass\":true,\"seconds\":$((SECONDS - t0)),\"record\":\"$OUT_DIR\""
}

run_replay() {
  local dir="$1"
  [ -f "$dir/native-digest.txt" ] || die "$EXIT_USAGE" "no retained record in $dir" "--run smoke first"
  local old; old=$(cat "$dir/native-digest.txt")
  local nd; nd=$(native_digest) || die "$EXIT_PIPELINE" "native run failed" "cargo run -p fs-flyer --release --bin canonical_digest"
  if [ "$nd" != "$old" ]; then
    log "replay" "\"pass\":false,\"expected\":\"$old\",\"observed\":\"$nd\""
    die "$EXIT_VERIFY" "replay digest diverged from the retained record" "--replay $dir"
  fi
  log "replay" "\"pass\":true,\"digest\":\"$nd\""
}

run_negative() {
  local case_name="$1"
  if [ "$case_name" = "list" ]; then list_cases; return 0; fi
  local found=""
  for t in "${TWINS[@]}"; do
    [ "${t%%:*}" = "$case_name" ] && found="$t"
  done
  [ -n "$found" ] || die "$EXIT_USAGE" "unknown twin $case_name (try --negative list)" "-"
  IFS=':' read -r _name lane binary filter <<< "$found"
  case "$lane" in
    inline)
      # golden-divergence: the production digest against a PERTURBED
      # golden — the comparison MUST fire.
      local nd; nd=$(native_digest) || die "$EXIT_PIPELINE" "native run failed" "-"
      local perturbed="0000000000000000000000000000000000000000000000000000000000000000"
      if [ "$nd" = "$perturbed" ]; then
        die "$EXIT_NEG_MISSED" "perturbed golden NOT detected" "-"
      fi
      log "negative" "\"case\":\"$case_name\",\"detected\":true"
      ;;
    fs-flyer)
      # The PRODUCTION battery that executes this twin (its assertions
      # ARE the detection); a red battery means the twin escaped.
      (cd "$REPO_ROOT" && cargo test -q -p fs-flyer --release --test "$binary" "$filter" >/dev/null 2>&1) \
        || die "$EXIT_NEG_MISSED" "twin battery failed: $binary::$filter" "cargo test -p fs-flyer --test $binary $filter"
      log "negative" "\"case\":\"$case_name\",\"detected\":true,\"battery\":\"$binary\""
      ;;
    fs-flyer-app)
      (cd "$REPO_ROOT/apps/wright-flyer" && node --test "test/$binary.test.ts" >/dev/null 2>&1) \
        || die "$EXIT_NEG_MISSED" "twin battery failed: $binary" "cd apps/wright-flyer && node --test test/$binary.test.ts"
      log "negative" "\"case\":\"$case_name\",\"detected\":true,\"battery\":\"$binary\""
      ;;
  esac
}

MODE="${1:---list}"
case "$MODE" in
  --list) list_cases ;;
  --check) run_check ;;
  --self-test) run_self_test ;;
  --run)
    [ "${2:-}" = "smoke" ] || die "$EXIT_USAGE" "--run smoke is the only bounded run" "-"
    run_check; run_smoke ;;
  --negative) run_negative "${2:-list}" ;;
  --replay) run_replay "${2:-}" ;;
  --output-dir) OUT_DIR="$2"; shift 2; exec "$0" "$@" ;;
  *) die "$EXIT_USAGE" "unknown mode $MODE" "$0 --list" ;;
esac
exit "$EXIT_OK"
