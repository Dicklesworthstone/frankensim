#!/usr/bin/env bash
# Foundations P0-P6 retained-evidence revalidation driver
# (bead frankensim-epic-foundations-huq.25).
#
# Replays what remains reproducible of the historical Foundations closure
# evidence against the CURRENT revision and classifies every manifest row.
# It never rewrites or reopens closed history and never infers proof from a
# closure date, an aggregate pass count, or the absence of a replay surface.
#
# Classification (driver-owned, per run):
#   Current        every replay command passed at this revision
#   Stale          a replay command failed at this revision (the historical
#                  claim no longer holds here; this is not a reopening)
#   Blocked        the milestone bead never closed - nothing to revalidate
#   NoData         the retained closure carries no evidence content and no
#                  replay surface
#   HistoricalOnly prose evidence whose declared file roots all exist; the
#                  claim survives as history, not as current proof
#   Unsupported    a declared file root is missing or a replay command
#                  cannot execute in this environment
#
# Modes:
#   --list                 enumerate manifest rows, run nothing
#   --check                validate the manifest schema, run nothing
#   --self-test            exercise the driver's own classifiers on fixtures
#   --run smoke            file-root checks + xtask-gate replays only
#   --run full             smoke plus the cargo test batteries
#   --replay RECEIPT       re-verify a retained summary against its log
#   --negative CASE        run one named hostile twin (or 'list'); PASS iff
#                          the driver refuses with the exact expected class
#   --output-dir DIR       fresh, repo-contained artifact root
#
# EXIT CLASSES: 0 = every row reached a terminal classification and no row
# classified Stale or Unsupported; 20 usage/manifest error; 21 a row is
# Stale; 22 a row is Unsupported; 23 replay tamper/disagreement.
# (Blocked/NoData/HistoricalOnly are honest terminals, not failures.)
set -u -o pipefail

SCHEMA="frankensim.foundations.p0p6-revalidation-manifest.v1"
LOG_SCHEMA="frankensim.foundations.p0p6-revalidation-log.v1"
OBS_IDENTITY_DOMAIN="org.frankensim.fs-obs.event-content.v10"

EXIT_OK=0
EXIT_USAGE=20
EXIT_STALE=21
EXIT_UNSUPPORTED=22
EXIT_TAMPER=23

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
MANIFEST="$REPO_ROOT/tests/foundations/manifests/p0-p6-revalidation.toml"
PER_COMMAND_TIMEOUT="${FSIM_P0P6_COMMAND_TIMEOUT_SECONDS:-1200}"

die() {
  local class="$1"; shift
  printf 'revalidate-p0-p6: ERROR class=%s: %s\n' "$class" "$*" >&2
  exit "$class"
}

command -v python3 >/dev/null 2>&1 || die "$EXIT_USAGE" "python3 (tomllib) is required"

MODE="" RUN_PROFILE="" OUTPUT_DIR="" REPLAY="" NEGATIVE_CASE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --list) MODE="list"; shift ;;
    --check) MODE="check"; shift ;;
    --self-test) MODE="self-test"; shift ;;
    --run)
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--run needs smoke|full"
      MODE="run"; RUN_PROFILE="$2"; shift 2 ;;
    --replay)
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--replay needs a receipt path"
      MODE="replay"; REPLAY="$2"; shift 2 ;;
    --negative)
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--negative needs a case name (or 'list')"
      MODE="negative"; NEGATIVE_CASE="$2"; shift 2 ;;
    --output-dir)
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--output-dir needs a value"
      OUTPUT_DIR="$2"; shift 2 ;;
    *) die "$EXIT_USAGE" "unknown argument: $1" ;;
  esac
done
[ -n "$MODE" ] || die "$EXIT_USAGE" "one of --list/--check/--self-test/--run/--replay is required"
if [ "$MODE" = "run" ] && [ "$RUN_PROFILE" != "smoke" ] && [ "$RUN_PROFILE" != "full" ]; then
  die "$EXIT_USAGE" "--run admits smoke or full, got: $RUN_PROFILE"
fi

# --------------------------------------------------------------- manifest --
project_manifest() { # path -> JSON rows on stdout, typed refusal on stderr
  python3 - "$1" "$SCHEMA" <<'PYEOF'
import json, sys
try:
    import tomllib
except ModuleNotFoundError:
    print("python3 tomllib unavailable (need >=3.11)", file=sys.stderr); sys.exit(20)
path, schema = sys.argv[1], sys.argv[2]
try:
    with open(path, "rb") as fh:
        doc = tomllib.load(fh)
except (OSError, tomllib.TOMLDecodeError) as error:
    print(f"manifest refusal: {error}", file=sys.stderr); sys.exit(20)
if doc.get("schema") != schema:
    print(f"manifest refusal: schema {doc.get('schema')!r} != {schema!r}", file=sys.stderr)
    sys.exit(20)
rows = doc.get("milestone")
if not isinstance(rows, list) or not rows:
    print("manifest refusal: needs [[milestone]] rows", file=sys.stderr); sys.exit(20)
seen = set()
for row in rows:
    for key, kind in (("bead", str), ("title", str), ("closed_at", str),
                      ("evidence_grade", str), ("replay_commands", list)):
        if not isinstance(row.get(key), kind):
            print(f"manifest refusal: row {row.get('bead')!r} field {key} missing/mistyped",
                  file=sys.stderr)
            sys.exit(20)
    if row["bead"] in seen:
        print(f"manifest refusal: duplicate bead {row['bead']!r}", file=sys.stderr)
        sys.exit(20)
    seen.add(row["bead"])
    if row["evidence_grade"] not in ("executable", "prose", "none"):
        print(f"manifest refusal: row {row['bead']!r} grade {row['evidence_grade']!r}",
              file=sys.stderr)
        sys.exit(20)
    for command in row["replay_commands"]:
        if not isinstance(command, list) or not all(isinstance(w, str) for w in command):
            print(f"manifest refusal: row {row['bead']!r} malformed replay command",
                  file=sys.stderr)
            sys.exit(20)
        for word in command:
            if word.startswith("/") or ".." in word.split("/"):
                print(f"manifest refusal: row {row['bead']!r} command word {word!r} escapes",
                      file=sys.stderr)
                sys.exit(20)
print(json.dumps(rows))
PYEOF
}

# ----------------------------------------------------------------- replay --
if [ "$MODE" = "replay" ]; then
  [ -f "$REPLAY" ] || die "$EXIT_TAMPER" "receipt not found: $REPLAY"
  python3 - "$REPLAY" <<'PYEOF'
import hashlib, json, os, sys
receipt_path = sys.argv[1]
summary = json.load(open(receipt_path))
log_path = os.path.join(os.path.dirname(receipt_path), summary.get("log_file", ""))
if not os.path.isfile(log_path):
    print("replay: log missing", file=sys.stderr); sys.exit(23)
if hashlib.sha256(open(log_path, "rb").read()).hexdigest() != summary.get("log_sha256"):
    print("replay: log digest mismatch", file=sys.stderr); sys.exit(23)
seq, terminals = -1, {}
for line in open(log_path):
    event = json.loads(line)
    if event["seq"] <= seq:
        print("replay: non-monotonic sequence", file=sys.stderr); sys.exit(23)
    seq = event["seq"]
    if event["event"] == "row-terminal":
        terminals[event["bead"]] = event["classification"]
if terminals != {row["bead"]: row["classification"] for row in summary["rows"]}:
    print("replay: terminal classifications disagree", file=sys.stderr); sys.exit(23)
print(f"replay OK: {len(terminals)} row terminals agree")
PYEOF
  exit $?
fi

# --------------------------------------------------------------- negative --
# Named hostile twins. Each constructs the hostile condition in scratch,
# runs the REAL driver against it, and passes iff the driver refuses with
# the exact expected nonzero class - a twin the driver survives is a FAIL.
if [ "$MODE" = "negative" ]; then
  NEG_CASES="severed-file-root stale-schema failing-replay tampered-log truncated-log"
  if [ "$NEGATIVE_CASE" = "list" ]; then
    printf '%s\n' $NEG_CASES
    exit "$EXIT_OK"
  fi
  case " $NEG_CASES " in
    *" $NEGATIVE_CASE "*) : ;;
    *) die "$EXIT_USAGE" "unknown negative case: $NEGATIVE_CASE (admitted: $NEG_CASES)" ;;
  esac
  SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/p0p6-negative.XXXXXX")"
  trap 'rm -rf "$SCRATCH"' EXIT
  FIX="$SCRATCH/fixture.toml"
  EXPECTED=0
  case "$NEGATIVE_CASE" in
    severed-file-root)
      # Evidence severing: the retained closure names a root that is gone.
      printf '%s\n' "schema = \"$SCHEMA\"" '[[milestone]]' 'bead = "twin"' \
        'title = "t"' 'closed_at = "2026-01-01"' 'evidence_grade = "prose"' \
        'replay_commands = []' 'file_roots = ["no/such/severed-root.txt"]' > "$FIX"
      EXPECTED="$EXIT_UNSUPPORTED"
      FSIM_P0P6_MANIFEST="$FIX" "$0" --run smoke --output-dir "$SCRATCH/out" >/dev/null 2>&1
      GOT=$? ;;
    stale-schema)
      printf '%s\n' 'schema = "frankensim.foundations.p0p6-revalidation-manifest.v0-stale"' > "$FIX"
      EXPECTED="$EXIT_USAGE"
      FSIM_P0P6_MANIFEST="$FIX" "$0" --run smoke --output-dir "$SCRATCH/out" >/dev/null 2>&1
      GOT=$? ;;
    failing-replay)
      # Input perturbation: the historical claim fails at this revision.
      printf '%s\n' "schema = \"$SCHEMA\"" '[[milestone]]' 'bead = "twin"' \
        'title = "t"' 'closed_at = "2026-01-01"' 'evidence_grade = "executable"' \
        'replay_commands = [["false"]]' > "$FIX"
      EXPECTED="$EXIT_STALE"
      FSIM_P0P6_MANIFEST="$FIX" "$0" --run full --output-dir "$SCRATCH/out" >/dev/null 2>&1
      GOT=$? ;;
    tampered-log|truncated-log)
      printf '%s\n' "schema = \"$SCHEMA\"" '[[milestone]]' 'bead = "twin"' \
        'title = "t"' 'closed_at = "2026-01-01"' 'evidence_grade = "executable"' \
        'replay_commands = [["true"]]' > "$FIX"
      FSIM_P0P6_MANIFEST="$FIX" "$0" --run full --output-dir "$SCRATCH/out" >/dev/null 2>&1 \
        || die "$EXIT_USAGE" "twin precondition run failed"
      if [ "$NEGATIVE_CASE" = "tampered-log" ]; then
        printf '%s\n' '{"seq":999,"event":"row-terminal","bead":"twin","classification":"Current"}' \
          >> "$SCRATCH/out/runner-log.jsonl"
      else
        sed -i '' -e '$d' "$SCRATCH/out/runner-log.jsonl" 2>/dev/null \
          || sed -i -e '$d' "$SCRATCH/out/runner-log.jsonl"
      fi
      EXPECTED="$EXIT_TAMPER"
      "$0" --replay "$SCRATCH/out/summary.json" >/dev/null 2>&1
      GOT=$? ;;
  esac
  if [ "$GOT" -eq "$EXPECTED" ]; then
    echo "negative twin '$NEGATIVE_CASE' PASS: driver refused with class $GOT"
    exit "$EXIT_OK"
  fi
  echo "negative twin '$NEGATIVE_CASE' FAIL: expected class $EXPECTED, got $GOT" >&2
  exit 1
fi

# -------------------------------------------------------------- self-test --
if [ "$MODE" = "self-test" ]; then
  SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/p0p6-selftest.XXXXXX")"
  trap 'rm -rf "$SCRATCH"' EXIT
  PASS=0; FAIL=0
  check() {
    if [ "$2" -eq "$3" ]; then PASS=$((PASS + 1)); else
      FAIL=$((FAIL + 1)); echo "SELF-TEST FAIL: $1 (expected $2, got $3)" >&2
    fi
  }
  FIX="$SCRATCH/fixture.toml"
  printf '%s\n' 'schema = "wrong.v9"' > "$FIX"
  project_manifest "$FIX" >/dev/null 2>&1
  check "wrong schema refuses as 20" 20 $?
  printf '%s\n' "schema = \"$SCHEMA\"" '[[milestone]]' 'bead = "b-1"' \
    'title = "t"' 'closed_at = "2026-01-01"' 'evidence_grade = "bogus"' \
    'replay_commands = []' > "$FIX"
  project_manifest "$FIX" >/dev/null 2>&1
  check "unknown evidence grade refuses as 20" 20 $?
  printf '%s\n' "schema = \"$SCHEMA\"" '[[milestone]]' 'bead = "b-1"' \
    'title = "t"' 'closed_at = "2026-01-01"' 'evidence_grade = "executable"' \
    'replay_commands = [["/bin/true"]]' > "$FIX"
  project_manifest "$FIX" >/dev/null 2>&1
  check "absolute command word refuses as 20" 20 $?
  # Classifier truths through a real tiny run.
  printf '%s\n' "schema = \"$SCHEMA\"" \
    '[[milestone]]' 'bead = "open-ms"' 'title = "t"' 'closed_at = ""' \
    'evidence_grade = "none"' 'replay_commands = []' \
    '[[milestone]]' 'bead = "nodata-ms"' 'title = "t"' 'closed_at = "2026-01-01"' \
    'evidence_grade = "none"' 'replay_commands = []' \
    '[[milestone]]' 'bead = "current-ms"' 'title = "t"' 'closed_at = "2026-01-01"' \
    'evidence_grade = "executable"' 'replay_commands = [["true"]]' \
    '[[milestone]]' 'bead = "stale-ms"' 'title = "t"' 'closed_at = "2026-01-01"' \
    'evidence_grade = "executable"' 'replay_commands = [["false"]]' > "$FIX"
  FSIM_P0P6_MANIFEST="$FIX" "$0" --run full --output-dir "$SCRATCH/run-1" >/dev/null 2>&1
  check "a stale row drives exit 21" 21 $?
  # The same fixture under smoke must NOT claim Current: non-gate replays
  # are deferred, so the executable rows classify HistoricalOnly.
  FSIM_P0P6_MANIFEST="$FIX" "$0" --run smoke --output-dir "$SCRATCH/run-smoke" >/dev/null 2>&1
  check "smoke defers non-gate replays without claiming Current" 0 $?
  python3 - "$SCRATCH/run-smoke/summary.json" <<'PYSMOKE'
import json, sys
rows = {row["bead"]: row["classification"] for row in json.load(open(sys.argv[1]))["rows"]}
sys.exit(0 if rows["current-ms"] == "HistoricalOnly" and rows["stale-ms"] == "HistoricalOnly" else 1)
PYSMOKE
  check "smoke-deferred rows read HistoricalOnly, never Current" 0 $?
  python3 - "$SCRATCH/run-1/summary.json" <<'PYEOF'
import json, sys
rows = {row["bead"]: row["classification"] for row in json.load(open(sys.argv[1]))["rows"]}
expected = {"open-ms": "Blocked", "nodata-ms": "NoData",
            "current-ms": "Current", "stale-ms": "Stale"}
sys.exit(0 if rows == expected else 1)
PYEOF
  check "classifications are exact (Blocked/NoData/Current/Stale)" 0 $?
  "$0" --replay "$SCRATCH/run-1/summary.json" >/dev/null 2>&1
  check "replay of an untouched run agrees" 0 $?
  printf '%s\n' '{"seq":999,"event":"row-terminal","bead":"x","classification":"Current"}' \
    >> "$SCRATCH/run-1/runner-log.jsonl"
  "$0" --replay "$SCRATCH/run-1/summary.json" >/dev/null 2>&1
  check "a tampered log fails replay as 23" 23 $?
  # A missing declared file root is Unsupported and drives exit 22.
  printf '%s\n' "schema = \"$SCHEMA\"" '[[milestone]]' 'bead = "gone-root"' \
    'title = "t"' 'closed_at = "2026-01-01"' 'evidence_grade = "prose"' \
    'replay_commands = []' 'file_roots = ["no/such/file.txt"]' > "$FIX"
  FSIM_P0P6_MANIFEST="$FIX" "$0" --run smoke --output-dir "$SCRATCH/run-2" >/dev/null 2>&1
  check "a missing file root drives exit 22 (Unsupported)" 22 $?
  echo "self-test: $PASS passed, $FAIL failed"
  [ "$FAIL" -eq 0 ] || exit 1
  exit 0
fi

MANIFEST="${FSIM_P0P6_MANIFEST:-$MANIFEST}"
ROWS_JSON="$(project_manifest "$MANIFEST")" || exit $?

if [ "$MODE" = "check" ]; then
  COUNT="$(printf '%s' "$ROWS_JSON" | python3 -c 'import json,sys; print(len(json.load(sys.stdin)))')"
  echo "manifest OK: $COUNT milestone row(s)"
  exit "$EXIT_OK"
fi
if [ "$MODE" = "list" ]; then
  printf '%s' "$ROWS_JSON" | python3 -c '
import json, sys
for row in json.load(sys.stdin):
    state = "closed" if row["closed_at"] else "OPEN"
    commands = "%d replay cmd(s)" % len(row["replay_commands"])
    print("\t".join([row["bead"], state, row["evidence_grade"], commands]))' || exit "$EXIT_USAGE"
  exit "$EXIT_OK"
fi

# -------------------------------------------------------------------- run --
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
[ -n "$OUTPUT_DIR" ] || OUTPUT_DIR="$REPO_ROOT/.e2e-out/p0p6-$STAMP-$$"
case "$OUTPUT_DIR" in
  "$REPO_ROOT"/*|"${TMPDIR:-/tmp}"*|/private/tmp/*|/tmp/*) : ;;
  *) die "$EXIT_USAGE" "--output-dir must stay inside the repository or TMPDIR" ;;
esac
[ -e "$OUTPUT_DIR" ] && die "$EXIT_USAGE" "output dir already exists (no reuse): $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR" || die "$EXIT_USAGE" "cannot create output dir"

LOG="$OUTPUT_DIR/runner-log.jsonl"
SEQ=0
HEAD_SHA="$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
emit() {
  SEQ=$((SEQ + 1))
  python3 - "$LOG" "$LOG_SCHEMA" "$OBS_IDENTITY_DOMAIN" "$SEQ" "$HEAD_SHA" "$@" <<'PYEOF'
import datetime, json, sys
log, schema, domain, seq, head = sys.argv[1:6]
event = {"schema": schema, "obs_identity_domain": domain, "seq": int(seq),
         "head_sha": head,
         "at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
         "event": sys.argv[6]}
for pair in sys.argv[7:]:
    key, _, value = pair.partition("=")
    event[key] = value
with open(log, "a") as fh:
    fh.write(json.dumps(event, sort_keys=True) + "\n")
PYEOF
}

emit run-start "profile=$RUN_PROFILE" "manifest=$MANIFEST"
WORST="$EXIT_OK"
declare -a ROWS=()

while IFS= read -r row_json; do
  BEAD="$(printf '%s' "$row_json" | python3 -c 'import json,sys; print(json.load(sys.stdin)["bead"])')"
  CLOSED_AT="$(printf '%s' "$row_json" | python3 -c 'import json,sys; print(json.load(sys.stdin)["closed_at"])')"
  GRADE="$(printf '%s' "$row_json" | python3 -c 'import json,sys; print(json.load(sys.stdin)["evidence_grade"])')"
  N_CMDS="$(printf '%s' "$row_json" | python3 -c 'import json,sys; print(len(json.load(sys.stdin)["replay_commands"]))')"
  emit row-start "bead=$BEAD" "grade=$GRADE"

  CLASSIFICATION="" DETAIL=""
  if [ -z "$CLOSED_AT" ]; then
    CLASSIFICATION="Blocked"; DETAIL="milestone never closed; nothing to revalidate"
  else
    # Declared file roots must exist regardless of profile.
    MISSING_ROOT="$(printf '%s' "$row_json" | python3 -c '
import json, os, sys
row = json.load(sys.stdin)
root_dir = sys.argv[1]
for root in row.get("file_roots", []):
    if not os.path.exists(os.path.join(root_dir, root)):
        print(root)
        break
' "$REPO_ROOT")"
    if [ -n "$MISSING_ROOT" ]; then
      CLASSIFICATION="Unsupported"; DETAIL="declared file root missing: $MISSING_ROOT"
    elif [ "$N_CMDS" -eq 0 ]; then
      if [ "$GRADE" = "none" ]; then
        CLASSIFICATION="NoData"; DETAIL="closure retains no evidence content"
      else
        CLASSIFICATION="HistoricalOnly"; DETAIL="prose evidence; declared roots present"
      fi
    else
      # smoke replays only xtask gates; full replays everything. A cargo
      # test battery skipped by profile leaves prose-grade truth: the row
      # classifies HistoricalOnly, never Current, because nothing ran.
      CLASSIFICATION=""; INDEX=0; RAN=0
      while IFS= read -r command_json; do
        INDEX=$((INDEX + 1))
        IS_GATE="$(printf '%s' "$command_json" | python3 -c 'import json,sys
words = json.load(sys.stdin)
print(1 if "xtask" in words else 0)')"
        if [ "$RUN_PROFILE" = "smoke" ] && [ "$IS_GATE" -ne 1 ]; then
          emit replay-skip "bead=$BEAD" "index=$INDEX" "reason=smoke-profile"
          continue
        fi
        declare -a WORDS=()
        while IFS= read -r word; do WORDS+=("$word"); done < <(
          printf '%s' "$command_json" | python3 -c 'import json,sys
for word in json.load(sys.stdin):
    print(word)')
        RAN=$((RAN + 1))
        ( cd "$REPO_ROOT" && timeout "$PER_COMMAND_TIMEOUT" nice -n 15 "${WORDS[@]}" ) \
          >"$OUTPUT_DIR/$BEAD.cmd$INDEX.out" 2>&1
        STATUS=$?
        emit replay-command "bead=$BEAD" "index=$INDEX" "exit=$STATUS"
        if [ "$STATUS" -ne 0 ]; then
          CLASSIFICATION="Stale"
          DETAIL="replay command $INDEX exited $STATUS at this revision"
          break
        fi
      done < <(printf '%s' "$row_json" | python3 -c 'import json,sys
for command in json.load(sys.stdin)["replay_commands"]:
    print(json.dumps(command))')
      if [ -z "$CLASSIFICATION" ]; then
        if [ "$RAN" -gt 0 ]; then
          CLASSIFICATION="Current"; DETAIL="$RAN replay command(s) passed"
        else
          CLASSIFICATION="HistoricalOnly"; DETAIL="all replays deferred by smoke profile"
        fi
      fi
    fi
  fi

  emit row-terminal "bead=$BEAD" "classification=$CLASSIFICATION" "detail=$DETAIL"
  ROWS+=("{\"bead\":\"$BEAD\",\"classification\":\"$CLASSIFICATION\"}")
  case "$CLASSIFICATION" in
    Stale) [ "$WORST" -lt "$EXIT_STALE" ] && WORST="$EXIT_STALE" ;;
    Unsupported) [ "$WORST" -lt "$EXIT_UNSUPPORTED" ] && WORST="$EXIT_UNSUPPORTED" ;;
  esac
done < <(printf '%s' "$ROWS_JSON" | python3 -c 'import json,sys
for row in json.load(sys.stdin):
    print(json.dumps(row))')

emit run-terminal "worst_class=$WORST"
python3 - "$OUTPUT_DIR/summary.json" "$LOG" "$WORST" "${ROWS[@]}" <<'PYEOF'
import hashlib, json, os, sys
summary_path, log_path, worst = sys.argv[1], sys.argv[2], int(sys.argv[3])
rows = [json.loads(row) for row in sys.argv[4:]]
counts = {}
for row in rows:
    counts[row["classification"]] = counts.get(row["classification"], 0) + 1
summary = {
    "schema": "frankensim.foundations.p0p6-revalidation-summary.v1",
    "log_file": os.path.basename(log_path),
    "log_sha256": hashlib.sha256(open(log_path, "rb").read()).hexdigest(),
    "worst_exit_class": worst,
    "counts": counts,
    "rows": rows,
}
with open(summary_path, "w") as fh:
    json.dump(summary, fh, indent=2, sort_keys=True)
    fh.write("\n")
print(f"summary: {json.dumps(counts, sort_keys=True)} -> exit {worst}")
print(f"retained: {summary_path}")
PYEOF

exit "$WORST"
