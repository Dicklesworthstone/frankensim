#!/usr/bin/env bash
# new-domains shared E2E runner (bead frankensim-ext-epic-gov-rjoq.8).
#
# Owns deterministic discovery, validation, isolation, execution, logging,
# and exit semantics for the E0a..E8 phase batteries. Phase batteries own
# fixtures, oracles, QoIs, tolerances, and promotion decisions in their
# versioned manifests at tests/e2e/new_domains/<phase>.toml; this runner
# never weakens a phase's scientific acceptance criteria.
#
# CLI:
#   --phase <e0a|e0b|e0c|e0d|e1|e2|e3|e4|e5|e6|e7|e8>
#   --case <id>            (repeatable; selects semantic case IDs)
#   --manifest <path>      (explicit manifest instead of phase discovery)
#   --list                 (list selected cases, run nothing)
#   --validate-only        (validate manifests, run nothing)
#   --output-dir <path>    (default: a fresh directory under the repo-local
#                           .e2e-out; never reused, never escaped)
#   --seed <u64>           (overrides a case seed ONLY if the manifest
#                           declares seed_overridable = true)
#   --max-wall-seconds <n> (tighten-only global budget override)
#   --replay <receipt>     (re-validate a retained summary against its logs)
#
# EXIT CLASSES (stable):
#    0 every selected case reached its declared terminal authority/refusal
#   10 runner usage / manifest schema error
#   11 admission refusal mismatch (expected refusal did not match)
#   12 production failure (entry point failed where authority was expected)
#   13 scientific acceptance failure (checker command refused the result)
#   14 timeout / budget exhaustion
#   17 tamper / checker failure (replay disagreement, log truncation)
#   18 infrastructure unqualified (missing tool, non-file manifest, escape)
#
# A skipped, filtered, missing, or unqualified case NEVER counts as pass.
# Logging: bounded schema-versioned JSONL, one event per line, referencing
# the frozen fs-obs event-content identity domain (V.3.1) rather than
# inventing a local logging authority.
set -u -o pipefail

SCHEMA_VERSION="frankensim.new-domains.case-manifest.v1"
LOG_SCHEMA="frankensim.new-domains.runner-log.v1"
OBS_IDENTITY_DOMAIN="org.frankensim.fs-obs.event-content.v10"
PHASES="e0a e0b e0c e0d e1 e2 e3 e4 e5 e6 e7 e8"

EXIT_OK=0
EXIT_USAGE=10
EXIT_ADMISSION=11
EXIT_PRODUCTION=12
EXIT_ACCEPTANCE=13
EXIT_BUDGET=14
EXIT_TAMPER=17
EXIT_UNQUALIFIED=18

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
MANIFEST_DIR="$REPO_ROOT/tests/e2e/new_domains"

die() { # class message
  local class="$1"; shift
  printf 'new-domains-e2e: ERROR class=%s: %s\n' "$class" "$*" >&2
  exit "$class"
}

command -v python3 >/dev/null 2>&1 || die "$EXIT_UNQUALIFIED" "python3 (tomllib) is required"

PHASE="" MANIFEST="" LIST=0 VALIDATE_ONLY=0 OUTPUT_DIR="" SEED="" MAX_WALL="" REPLAY=""
declare -a CASES=()
while [ $# -gt 0 ]; do
  case "$1" in
    --phase) [ $# -ge 2 ] || die "$EXIT_USAGE" "--phase needs a value"; PHASE="$2"; shift 2 ;;
    --case) [ $# -ge 2 ] || die "$EXIT_USAGE" "--case needs a value"; CASES+=("$2"); shift 2 ;;
    --manifest) [ $# -ge 2 ] || die "$EXIT_USAGE" "--manifest needs a value"; MANIFEST="$2"; shift 2 ;;
    --list) LIST=1; shift ;;
    --validate-only) VALIDATE_ONLY=1; shift ;;
    --output-dir) [ $# -ge 2 ] || die "$EXIT_USAGE" "--output-dir needs a value"; OUTPUT_DIR="$2"; shift 2 ;;
    --seed) [ $# -ge 2 ] || die "$EXIT_USAGE" "--seed needs a value"; SEED="$2"; shift 2 ;;
    --max-wall-seconds) [ $# -ge 2 ] || die "$EXIT_USAGE" "--max-wall-seconds needs a value"; MAX_WALL="$2"; shift 2 ;;
    --replay) [ $# -ge 2 ] || die "$EXIT_USAGE" "--replay needs a value"; REPLAY="$2"; shift 2 ;;
    *) die "$EXIT_USAGE" "unknown argument: $1" ;;
  esac
done

if [ -n "$SEED" ] && ! printf '%s' "$SEED" | grep -Eq '^[0-9]{1,20}$'; then
  die "$EXIT_USAGE" "--seed must be an unsigned integer, got: $SEED"
fi
if [ -n "$MAX_WALL" ] && ! printf '%s' "$MAX_WALL" | grep -Eq '^[1-9][0-9]{0,8}$'; then
  die "$EXIT_USAGE" "--max-wall-seconds must be a positive integer"
fi
if [ -n "$PHASE" ]; then
  case " $PHASES " in
    *" $PHASE "*) : ;;
    *) die "$EXIT_USAGE" "unknown phase: $PHASE (admitted: $PHASES)" ;;
  esac
fi
if [ -n "$PHASE" ] && [ -n "$MANIFEST" ]; then
  die "$EXIT_USAGE" "--phase and --manifest are mutually exclusive"
fi

# ---------------------------------------------------------------- replay --
if [ -n "$REPLAY" ]; then
  [ -f "$REPLAY" ] || die "$EXIT_UNQUALIFIED" "replay receipt not found: $REPLAY"
  python3 - "$REPLAY" <<'PYEOF'
import json, sys, os, hashlib
receipt_path = sys.argv[1]
with open(receipt_path) as fh:
    summary = json.load(fh)
log_path = os.path.join(os.path.dirname(receipt_path), summary.get("log_file", ""))
if not os.path.isfile(log_path):
    print(f"replay: log file missing: {log_path}", file=sys.stderr); sys.exit(17)
digest = hashlib.sha256(open(log_path, "rb").read()).hexdigest()
if digest != summary.get("log_sha256"):
    print("replay: log digest mismatch (tamper or truncation)", file=sys.stderr); sys.exit(17)
seq = -1
terminals = {}
for line in open(log_path):
    event = json.loads(line)
    if event["seq"] <= seq:
        print("replay: non-monotonic sequence", file=sys.stderr); sys.exit(17)
    seq = event["seq"]
    if event["event"] == "case-terminal":
        terminals[event["case"]] = event["status"]
recorded = {row["case"]: row["status"] for row in summary["cases"]}
if terminals != recorded:
    print(f"replay: terminal statuses disagree: log={terminals} summary={recorded}", file=sys.stderr)
    sys.exit(17)
print(f"replay OK: {len(terminals)} case terminals agree with the retained summary")
PYEOF
  exit $?
fi

# ------------------------------------------------------------- discovery --
declare -a MANIFESTS=()
if [ -n "$MANIFEST" ]; then
  [ -f "$MANIFEST" ] || die "$EXIT_UNQUALIFIED" "manifest is not a file: $MANIFEST"
  MANIFESTS+=("$MANIFEST")
else
  [ -d "$MANIFEST_DIR" ] || die "$EXIT_UNQUALIFIED" "manifest directory missing: $MANIFEST_DIR"
  # Deterministic discovery: bytewise-sorted phase manifests.
  while IFS= read -r file; do
    base="$(basename "$file" .toml)"
    if [ -n "$PHASE" ] && [ "$base" != "$PHASE" ]; then continue; fi
    MANIFESTS+=("$file")
  done < <(LC_ALL=C ls "$MANIFEST_DIR"/*.toml 2>/dev/null | LC_ALL=C sort)
  if [ "${#MANIFESTS[@]}" -eq 0 ]; then
    die "$EXIT_USAGE" "no manifests selected (phase=${PHASE:-<all>}) in $MANIFEST_DIR"
  fi
fi

# ------------------------------------------- validation and case listing --
# Parses every selected manifest, applies the schema, and emits one JSON
# object per admitted case on stdout. Any violation is a typed refusal.
validate_and_project() {
  python3 - "$SCHEMA_VERSION" "$REPO_ROOT" "${MANIFESTS[@]}" <<'PYEOF'
import json, sys, os, re
try:
    import tomllib
except ModuleNotFoundError:
    print("python3 tomllib unavailable (need >=3.11)", file=sys.stderr); sys.exit(18)
schema, repo_root, paths = sys.argv[1], sys.argv[2], sys.argv[3:]

REQUIRED_CASE_KEYS = {
    "id": str, "version": int, "purpose": str, "owning_bead": str,
    "gauntlet_tier": str, "entry_command": list, "seed": int,
    "max_wall_seconds": int, "expected": str,
}
OPTIONAL_CASE_KEYS = {
    "checker_command": list, "expected_refusal_pattern": str,
    "seed_overridable": bool, "expected_output_pattern": str,
    "qoi_notes": str, "determinism_class": str,
}
ADMITTED_EXPECTED = {"authority", "refusal"}

def refuse(msg):
    print(f"manifest refusal: {msg}", file=sys.stderr)
    sys.exit(10)

seen_ids = set()
projected = []
for path in paths:
    with open(path, "rb") as fh:
        try:
            doc = tomllib.load(fh)
        except tomllib.TOMLDecodeError as error:
            refuse(f"{path}: not valid TOML: {error}")
    if doc.get("schema") != schema:
        refuse(f"{path}: schema is {doc.get('schema')!r}, runner admits only {schema!r}")
    phase = doc.get("phase")
    expected_phase = os.path.basename(path)[:-5]
    if phase != expected_phase:
        refuse(f"{path}: phase field {phase!r} must equal file stem {expected_phase!r}")
    cases = doc.get("case")
    if not isinstance(cases, list) or not cases:
        refuse(f"{path}: needs at least one [[case]]")
    unknown_top = set(doc) - {"schema", "phase", "case"}
    if unknown_top:
        refuse(f"{path}: unknown top-level fields {sorted(unknown_top)} (no silent semantics)")
    for case in cases:
        unknown = set(case) - set(REQUIRED_CASE_KEYS) - set(OPTIONAL_CASE_KEYS)
        if unknown:
            refuse(f"{path}: case {case.get('id')!r} has unknown fields {sorted(unknown)}")
        for key, kind in REQUIRED_CASE_KEYS.items():
            if not isinstance(case.get(key), kind):
                refuse(f"{path}: case {case.get('id')!r} field {key!r} missing or not {kind.__name__}")
        for key, kind in OPTIONAL_CASE_KEYS.items():
            if key in case and not isinstance(case[key], kind):
                refuse(f"{path}: case {case['id']!r} field {key!r} not {kind.__name__}")
        if not re.fullmatch(r"[a-z0-9][a-z0-9-]{2,63}", case["id"]):
            refuse(f"{path}: case id {case['id']!r} is not a stable slug")
        if case["id"] in seen_ids:
            refuse(f"duplicate semantic case id {case['id']!r}")
        seen_ids.add(case["id"])
        if case["expected"] not in ADMITTED_EXPECTED:
            refuse(f"{path}: case {case['id']!r} expected must be one of {sorted(ADMITTED_EXPECTED)}")
        if case["expected"] == "refusal" and "expected_refusal_pattern" not in case:
            refuse(f"{path}: refusal case {case['id']!r} must declare expected_refusal_pattern")
        if case["max_wall_seconds"] < 1 or case["max_wall_seconds"] > 7200:
            refuse(f"{path}: case {case['id']!r} max_wall_seconds outside [1, 7200]")
        for word in case["entry_command"] + case.get("checker_command", []):
            if not isinstance(word, str):
                refuse(f"{path}: case {case['id']!r} command words must be strings")
            if word.startswith("/") or ".." in word.split(os.sep):
                refuse(f"{path}: case {case['id']!r} command word {word!r} escapes the repo (absolute or ..)")
        case["_phase"] = phase
        case["_manifest"] = os.path.relpath(path, repo_root)
        projected.append(case)
print(json.dumps(projected))
PYEOF
}

PROJECTED="$(validate_and_project)" || exit $?

# Case selection (repeatable --case); unknown selections are usage errors.
if [ "${#CASES[@]}" -gt 0 ]; then
  PROJECTED="$(printf '%s' "$PROJECTED" | python3 -c '
import json, sys
selected = set(sys.argv[1:])
cases = json.load(sys.stdin)
known = {case["id"] for case in cases}
missing = selected - known
if missing:
    print(f"unknown case id(s): {sorted(missing)}", file=sys.stderr)
    sys.exit(10)
print(json.dumps([case for case in cases if case["id"] in selected]))
' "${CASES[@]}")" || exit $?
fi

COUNT="$(printf '%s' "$PROJECTED" | python3 -c 'import json,sys; print(len(json.load(sys.stdin)))')"

if [ "$VALIDATE_ONLY" -eq 1 ]; then
  echo "validate-only OK: ${#MANIFESTS[@]} manifest(s), $COUNT case(s) admitted"
  exit "$EXIT_OK"
fi
if [ "$LIST" -eq 1 ]; then
  printf '%s' "$PROJECTED" | python3 -c '
import json, sys
for case in json.load(sys.stdin):
    row = [case["_phase"], case["id"], "v%d" % case["version"], case["expected"], case["_manifest"]]
    print("\t".join(row))' || exit "$EXIT_UNQUALIFIED"
  exit "$EXIT_OK"
fi
[ "$COUNT" -gt 0 ] || die "$EXIT_USAGE" "no cases selected"

# -------------------------------------------------------------- execution --
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
if [ -z "$OUTPUT_DIR" ]; then
  OUTPUT_DIR="$REPO_ROOT/.e2e-out/new-domains-$STAMP-$$"
fi
case "$OUTPUT_DIR" in
  "$REPO_ROOT"/*) : ;;
  *) die "$EXIT_UNQUALIFIED" "--output-dir must stay inside the repository: $OUTPUT_DIR" ;;
esac
[ -e "$OUTPUT_DIR" ] && die "$EXIT_UNQUALIFIED" "output dir already exists (no reuse): $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR" || die "$EXIT_UNQUALIFIED" "cannot create output dir"

LOG="$OUTPUT_DIR/runner-log.jsonl"
SUMMARY="$OUTPUT_DIR/summary.json"
SEQ=0
HEAD_SHA="$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"

emit() { # event-name json-fields...
  SEQ=$((SEQ + 1))
  python3 - "$LOG" "$LOG_SCHEMA" "$OBS_IDENTITY_DOMAIN" "$SEQ" "$HEAD_SHA" "$@" <<'PYEOF'
import json, sys, datetime
log, schema, domain, seq, head = sys.argv[1:6]
event = {"schema": schema, "obs_identity_domain": domain, "seq": int(seq),
         "head_sha": head,
         "at": datetime.datetime.now(datetime.timezone.utc).isoformat()}
event["event"] = sys.argv[6]
for pair in sys.argv[7:]:
    key, _, value = pair.partition("=")
    event[key] = value
with open(log, "a") as fh:
    fh.write(json.dumps(event, sort_keys=True) + "\n")
PYEOF
}

emit run-start "manifest_count=${#MANIFESTS[@]}" "case_count=$COUNT" \
  "seed_override=${SEED:-none}" "max_wall_override=${MAX_WALL:-none}"

WORST="$EXIT_OK"
declare -a ROWS=()
while IFS= read -r case_json; do
  ID="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])')"
  EXPECTED="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin)["expected"])')"
  WALL="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin)["max_wall_seconds"])')"
  CASE_SEED="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin)["seed"])')"
  OVERRIDABLE="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(str(json.load(sys.stdin).get("seed_overridable", False)).lower())')"
  REFUSAL_PATTERN="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin).get("expected_refusal_pattern", ""))')"
  OUTPUT_PATTERN="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin).get("expected_output_pattern", ""))')"

  # Budget overrides only tighten the admitted manifest.
  if [ -n "$MAX_WALL" ] && [ "$MAX_WALL" -lt "$WALL" ]; then WALL="$MAX_WALL"; fi
  if [ -n "$SEED" ]; then
    if [ "$OVERRIDABLE" != "true" ]; then
      emit case-terminal "case=$ID" "status=unqualified" "reason=seed-not-overridable"
      ROWS+=("{\"case\":\"$ID\",\"status\":\"unqualified\"}")
      [ "$WORST" -lt "$EXIT_UNQUALIFIED" ] && WORST="$EXIT_UNQUALIFIED"
      continue
    fi
    CASE_SEED="$SEED"
  fi

  CASE_DIR="$OUTPUT_DIR/$ID"
  mkdir -p "$CASE_DIR"
  emit case-start "case=$ID" "expected=$EXPECTED" "seed=$CASE_SEED" "wall_budget=$WALL"

  set -- $(printf '%s' "$case_json" | python3 -c 'import json,sys
for word in json.load(sys.stdin)["entry_command"]:
    print(word)')
  START=$(date +%s)
  ( cd "$REPO_ROOT" && NEW_DOMAINS_SEED="$CASE_SEED" NEW_DOMAINS_CASE_DIR="$CASE_DIR" \
      timeout "$WALL" "$@" ) >"$CASE_DIR/stdout.txt" 2>"$CASE_DIR/stderr.txt"
  STATUS_CODE=$?
  ELAPSED=$(( $(date +%s) - START ))

  RESULT="" REASON=""
  if [ "$STATUS_CODE" -eq 124 ]; then
    RESULT="failed"; REASON="wall-budget-exhausted"; CLASS="$EXIT_BUDGET"
  elif [ "$EXPECTED" = "authority" ]; then
    if [ "$STATUS_CODE" -eq 0 ]; then
      if [ -n "$OUTPUT_PATTERN" ] && ! grep -Eq "$OUTPUT_PATTERN" "$CASE_DIR/stdout.txt" "$CASE_DIR/stderr.txt"; then
        RESULT="failed"; REASON="authority-output-pattern-unmatched"; CLASS="$EXIT_ACCEPTANCE"
      else
        RESULT="passed"; REASON="authority"; CLASS="$EXIT_OK"
      fi
    else
      RESULT="failed"; REASON="production-exit-$STATUS_CODE"; CLASS="$EXIT_PRODUCTION"
    fi
  else # expected refusal
    if [ "$STATUS_CODE" -ne 0 ] && grep -Eq "$REFUSAL_PATTERN" "$CASE_DIR/stdout.txt" "$CASE_DIR/stderr.txt"; then
      RESULT="expected-refused"; REASON="refusal-matched"; CLASS="$EXIT_OK"
    elif [ "$STATUS_CODE" -eq 0 ]; then
      RESULT="failed"; REASON="expected-refusal-but-authority"; CLASS="$EXIT_ADMISSION"
    else
      RESULT="failed"; REASON="refusal-pattern-unmatched-exit-$STATUS_CODE"; CLASS="$EXIT_ADMISSION"
    fi
  fi

  # Independent checker phase (only meaningful after a passing production).
  if [ "$RESULT" = "passed" ]; then
    CHECKER_LEN="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(len(json.load(sys.stdin).get("checker_command", [])))')"
    if [ "$CHECKER_LEN" -gt 0 ]; then
      set -- $(printf '%s' "$case_json" | python3 -c 'import json,sys
for word in json.load(sys.stdin).get("checker_command", []):
    print(word)')
      ( cd "$REPO_ROOT" && NEW_DOMAINS_CASE_DIR="$CASE_DIR" timeout "$WALL" "$@" ) \
        >"$CASE_DIR/checker-stdout.txt" 2>"$CASE_DIR/checker-stderr.txt"
      if [ $? -ne 0 ]; then
        RESULT="failed"; REASON="checker-refused"; CLASS="$EXIT_ACCEPTANCE"
      fi
    fi
  fi

  emit case-terminal "case=$ID" "status=$RESULT" "reason=$REASON" \
    "exit_code=$STATUS_CODE" "elapsed_seconds=$ELAPSED"
  ROWS+=("{\"case\":\"$ID\",\"status\":\"$RESULT\"}")
  [ "$WORST" -lt "$CLASS" ] && WORST="$CLASS"
done < <(printf '%s' "$PROJECTED" | python3 -c '
import json, sys
for case in json.load(sys.stdin):
    print(json.dumps(case))')

emit run-terminal "worst_class=$WORST"

python3 - "$SUMMARY" "$LOG" "$WORST" "${ROWS[@]}" <<'PYEOF'
import json, sys, hashlib, os
summary_path, log_path, worst = sys.argv[1], sys.argv[2], int(sys.argv[3])
rows = [json.loads(row) for row in sys.argv[4:]]
counts = {}
for row in rows:
    counts[row["status"]] = counts.get(row["status"], 0) + 1
summary = {
    "schema": "frankensim.new-domains.runner-summary.v1",
    "log_file": os.path.basename(log_path),
    "log_sha256": hashlib.sha256(open(log_path, "rb").read()).hexdigest(),
    "worst_exit_class": worst,
    "counts": counts,
    "cases": rows,
}
with open(summary_path, "w") as fh:
    json.dump(summary, fh, indent=2, sort_keys=True)
    fh.write("\n")
print(f"summary: {json.dumps(counts, sort_keys=True)} -> exit {worst}")
print(f"retained: {summary_path}")
PYEOF

exit "$WORST"
