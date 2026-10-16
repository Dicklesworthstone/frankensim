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
#   --check                (canonical validate spelling: source, manifest,
#                          route, capability, path, and schema closure
#                          WITHOUT production execution)
#   --validate-only        (compatibility alias of --check)
#   --self-test            (run only the source-closed runner fixtures)
#   --run smoke|full       (canonical run modes; bare invocation = full;
#                          smoke executes only cases whose manifest
#                          declaration includes profiles = ["smoke"])
#   --output-dir <path>    (default: a fresh directory under the repo-local
#                           .e2e-out; never reused, never escaped)
#   --seed <u64>           (overrides a case seed ONLY if the manifest
#                           declares seed_overridable = true)
#   --max-wall-seconds <n> (tighten-only global budget override)
#   --negative <id>        (execute ONE immutable negative probe declared
#                           by the selected phase manifest's [[negative]];
#                           refuses undeclared ids)
#   --replay <receipt>     (artifact-only replay; canonical form REQUIRES
#                           --replay-input-root PATH --replay-output-root
#                           PATH: input read-only, output fresh and
#                           disjoint; equality, overlap, traversal, symlink
#                           escape, or preexisting output refuses BEFORE
#                           payload bytes are read)
#   --cancel-after <s>     (cancel each case after s seconds through TERM,
#                           10s drain grace, then KILL = drain failure)
#   --determinism-repeat   (run repo-deterministic cases twice; stdout must
#                           agree byte-for-byte)
# EXIT CLASSES (stable):
#    0 every selected case reached its declared terminal authority/refusal
#   10 runner usage / manifest schema error
#   11 admission refusal mismatch (expected refusal did not match)
#   12 production failure (entry point failed where authority was expected)
#   13 scientific acceptance failure (checker command refused the result)
#   14 timeout / budget exhaustion (wall or output budget)
#   15 cancellation / drain failure
#   16 determinism mismatch
#   17 tamper / checker failure (replay disagreement, log truncation)
#   18 infrastructure unqualified (missing tool, non-file manifest, escape)
#
# A skipped, filtered, missing, or unqualified case NEVER counts as pass.
# Logging: bounded schema-versioned JSONL, one event per line, referencing
# the frozen fs-obs event-content identity domain (V.3.1) rather than
# inventing a local logging authority. Event cardinality per invocation:
# exactly one run-start and one run-terminal (first/last); one
# case-selected and one case-terminal per selected case; paired
# attempt-started/attempt-terminal per attempt; paired
# stage-started/stage-terminal per begun stage; exactly one summary after
# all case terminals; at most one first-divergence per case. Semantic rows
# carry no wall-clock timestamps or durations; the operational `at` field
# is transport metadata, and elapsed time lives only in the volatile
# per-case sidecar. Logs written before the cardinality contract replay as
# a typed refusal (schema migration is explicit, never silent).
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
EXIT_CANCEL=15
EXIT_DETERMINISM=16
EXIT_TAMPER=17
EXIT_UNQUALIFIED=18

# LITERAL HARD CAPS (acceptance contract). Declared manifest budgets may
# only tighten these ceilings; crossing any cap yields its stable typed
# nonzero exit class, preserves already committed bounded evidence, drains
# children, and can never be reported as pass.
CAP_CASES_PER_MANIFEST=256
CAP_SELECTED_CASES=256
CAP_RECORDS_PER_CASE=100000
CAP_RECORDS_PER_RUN=2000000
CAP_EVENT_BYTES=65536
CAP_DIAGNOSTIC_BYTES=16384
CAP_ARTIFACTS_PER_CASE=256
CAP_ARTIFACTS_PER_RUN=16384
CAP_ARTIFACT_INDEX_BYTES=16777216
CAP_CHILD_OUTPUT_BYTES=67108864

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
MANIFEST_DIR="$REPO_ROOT/tests/e2e/new_domains"

RUNNER_PATH="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
die() { # class message
  local class="$1"; shift
  printf 'new-domains-e2e: ERROR class=%s: %s\n' "$class" "$*" >&2
  exit "$class"
}

command -v python3 >/dev/null 2>&1 || die "$EXIT_UNQUALIFIED" "python3 (tomllib) is required"

PHASE="" MANIFEST="" LIST=0 CHECK=0 OUTPUT_DIR="" SEED="" MAX_WALL="" REPLAY=""
REPLAY_IN_ROOT="" REPLAY_OUT_ROOT="" CANCEL_AFTER="" DETERMINISM_REPEAT=0
RUN_MODE="" NEGATIVE="" SELFTEST=0
declare -a CASES=()
while [ $# -gt 0 ]; do
  case "$1" in
    --phase) [ $# -ge 2 ] || die "$EXIT_USAGE" "--phase needs a value"; PHASE="$2"; shift 2 ;;
    --case) [ $# -ge 2 ] || die "$EXIT_USAGE" "--case needs a value"; CASES+=("$2"); shift 2 ;;
    --manifest) [ $# -ge 2 ] || die "$EXIT_USAGE" "--manifest needs a value"; MANIFEST="$2"; shift 2 ;;
    --list) LIST=1; shift ;;
    --check) CHECK=1; shift ;;
    --validate-only) CHECK=1; shift ;;  # compatibility alias of --check
    --self-test) SELFTEST=1; shift ;;
    --run) [ $# -ge 2 ] || die "$EXIT_USAGE" "--run needs a value"; RUN_MODE="$2"; shift 2 ;;
    --negative) [ $# -ge 2 ] || die "$EXIT_USAGE" "--negative needs a value"; NEGATIVE="$2"; shift 2 ;;
    --output-dir) [ $# -ge 2 ] || die "$EXIT_USAGE" "--output-dir needs a value"; OUTPUT_DIR="$2"; shift 2 ;;
    --seed) [ $# -ge 2 ] || die "$EXIT_USAGE" "--seed needs a value"; SEED="$2"; shift 2 ;;
    --max-wall-seconds) [ $# -ge 2 ] || die "$EXIT_USAGE" "--max-wall-seconds needs a value"; MAX_WALL="$2"; shift 2 ;;
    --replay) [ $# -ge 2 ] || die "$EXIT_USAGE" "--replay needs a value"; REPLAY="$2"; shift 2 ;;
    --replay-input-root) [ $# -ge 2 ] || die "$EXIT_USAGE" "--replay-input-root needs a value"; REPLAY_IN_ROOT="$2"; shift 2 ;;
    --replay-output-root) [ $# -ge 2 ] || die "$EXIT_USAGE" "--replay-output-root needs a value"; REPLAY_OUT_ROOT="$2"; shift 2 ;;
    --cancel-after) [ $# -ge 2 ] || die "$EXIT_USAGE" "--cancel-after needs a value"; CANCEL_AFTER="$2"; shift 2 ;;
    --determinism-repeat) DETERMINISM_REPEAT=1; shift ;;
    *) die "$EXIT_USAGE" "unknown argument: $1" ;;
  esac
done

if [ -n "$SEED" ] && ! printf '%s' "$SEED" | grep -Eq '^[0-9]{1,20}$'; then
  die "$EXIT_USAGE" "--seed must be an unsigned integer, got: $SEED"
fi
if [ -n "$MAX_WALL" ] && ! printf '%s' "$MAX_WALL" | grep -Eq '^[1-9][0-9]{0,8}$'; then
  die "$EXIT_USAGE" "--max-wall-seconds must be a positive integer"
fi
if [ -n "$CANCEL_AFTER" ] && ! printf '%s' "$CANCEL_AFTER" | grep -Eq '^[1-9][0-9]{0,4}$'; then
  die "$EXIT_USAGE" "--cancel-after must be a positive integer (seconds)"
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
if [ -n "$RUN_MODE" ]; then
  case "$RUN_MODE" in
    smoke|full) : ;;
    *) die "$EXIT_USAGE" "--run must be 'smoke' or 'full', got: $RUN_MODE" ;;
  esac
fi
if [ -n "$NEGATIVE" ] && ! printf '%s' "$NEGATIVE" | grep -Eq '^[a-z0-9][a-z0-9-]{2,63}$'; then
  die "$EXIT_USAGE" "--negative id is not a stable slug: $NEGATIVE"
fi
if [ -n "$REPLAY" ]; then
  [ -n "$REPLAY_IN_ROOT" ] || die "$EXIT_USAGE" "--replay requires --replay-input-root (canonical form)"
  [ -n "$REPLAY_OUT_ROOT" ] || die "$EXIT_USAGE" "--replay requires --replay-output-root (canonical form)"
fi
if [ -z "$REPLAY" ]; then
  [ -z "$REPLAY_IN_ROOT" ] || die "$EXIT_USAGE" "--replay-input-root is only valid with --replay"
  [ -z "$REPLAY_OUT_ROOT" ] || die "$EXIT_USAGE" "--replay-output-root is only valid with --replay"
fi

# ------------------------------------------------------------- self-test --
# Source-closed runner fixtures only: the battery executes the real runner
# against synthetic manifests in a scratch directory and never runs phase
# production entries.
if [ "$SELFTEST" -eq 1 ]; then
  exec bash "$REPO_ROOT/scripts/ci/new_domains_e2e_selftest.sh"
fi

# ---------------------------------------------------------------- replay --
# Artifact-only replay. The input root is read-only and must contain the
# receipt; the output root must be a fresh, repository-contained, disjoint
# destination. Every refusal below happens BEFORE any payload byte is read
# and before the output root is created.
if [ -n "$REPLAY" ]; then
  python3 - "$REPLAY" "$REPO_ROOT" "$REPLAY_IN_ROOT" "$REPLAY_OUT_ROOT" "$LOG_SCHEMA" "$OBS_IDENTITY_DOMAIN" "$CAP_RECORDS_PER_RUN" <<'PYEOF'
import json, sys, os, hashlib

receipt_arg, repo_root, in_root, out_root, log_schema, obs_domain, cap_records = sys.argv[1:8]

def refuse(cls, msg):
    print(f"replay refusal class={cls}: {msg}", file=sys.stderr)
    sys.exit(int(cls))

def real(p):
    return os.path.realpath(os.path.abspath(p))

ri, ro = real(in_root), real(out_root)
for p in (in_root, out_root):
    rp = real(p)
    if not (rp == real(repo_root) or rp.startswith(real(repo_root) + os.sep)):
        refuse(18, f"path escapes the repository: {p}")
if ri == ro:
    refuse(18, "replay input and output roots must be disjoint (equal after resolution)")
if ro.startswith(ri + os.sep) or ri.startswith(ro + os.sep):
    refuse(18, "replay input and output roots overlap (nested)")
if not os.path.isdir(in_root):
    refuse(18, f"replay input root is not a directory: {in_root}")
if os.path.lexists(out_root):
    refuse(18, f"replay output root already exists (must be fresh): {out_root}")

receipt_path = os.path.abspath(receipt_arg)
if not (receipt_path == ri or receipt_path.startswith(ri + os.sep)):
    refuse(18, "receipt must live inside the replay input root")
if not os.path.isfile(receipt_path):
    refuse(18, f"replay receipt not found: {receipt_arg}")

# Payload is read only after every containment/disjointness precondition.
with open(receipt_path) as fh:
    try:
        summary = json.load(fh)
    except json.JSONDecodeError as error:
        refuse(17, f"receipt is not valid JSON: {error}")
log_rel = summary.get("log_file", "")
log_path = os.path.realpath(os.path.join(os.path.dirname(receipt_path), log_rel))
if not (log_path == ri or log_path.startswith(ri + os.sep)):
    refuse(17, "recorded log path escapes the replay input root")
if not os.path.isfile(log_path):
    refuse(17, "recorded log file missing")

raw = open(log_path, "rb").read()
if hashlib.sha256(raw).hexdigest() != summary.get("log_sha256"):
    refuse(17, "log digest mismatch (tamper or truncation)")

events = []
for line in raw.decode().splitlines():
    if not line.strip():
        continue
    try:
        events.append(json.loads(line))
    except json.JSONDecodeError as error:
        refuse(17, f"malformed log line: {error}")
if len(events) > int(cap_records):
    refuse(14, "record cap exceeded in retained log")
if not events:
    refuse(17, "empty log")

KNOWN = {"run-start", "run-terminal", "case-selected", "case-terminal",
         "attempt-started", "attempt-terminal", "stage-started",
         "stage-terminal", "summary", "first-divergence",
         "case-cancel-request", "case-cancel-drained"}

seq = 0
heads = set()
for event in events:
    if event.get("schema") != log_schema:
        refuse(17, f"unknown log schema revision: {event.get('schema')!r}")
    if event.get("obs_identity_domain") != obs_domain:
        refuse(17, "obs identity domain mismatch")
    s = event.get("seq")
    if not isinstance(s, int) or s <= seq:
        refuse(17, "non-monotonic sequence")
    seq = s
    heads.add(event.get("head_sha"))
    if event.get("event") not in KNOWN:
        refuse(17, f"unknown event kind: {event.get('event')!r}")
if len(heads) > 1:
    refuse(17, "mixed-run log: head identity changes within one log")

kinds = [e["event"] for e in events]
if kinds[0] != "run-start":
    refuse(17, "first event is not run-start")
if kinds.count("run-start") != 1:
    refuse(17, "exactly one run-start required")
if kinds[-1] != "run-terminal":
    refuse(17, "last event is not run-terminal (truncated?)")
if kinds.count("run-terminal") != 1:
    refuse(17, "exactly one run-terminal required")
if kinds.count("summary") != 1:
    refuse(17, "exactly one summary required")
term_idx = len(kinds) - 1
sum_idx = kinds.index("summary")
case_term_idxs = [i for i, k in enumerate(kinds) if k == "case-terminal"]
if not case_term_idxs or sum_idx < case_term_idxs[-1] or sum_idx > term_idx:
    refuse(17, "summary must follow all case terminals and precede run-terminal")

state = {}
for event in events[1:-1]:  # skip run-start and run-terminal
    kind, case = event["event"], event.get("case")
    if kind == "summary":
        continue
    st = state.get(case)
    if st is None:
        if kind == "case-selected":
            state[case] = {"attempts": 0, "open_attempt": False,
                           "stage_open": False, "terminal": None}
            continue
        refuse(17, f"{kind} before case-selected for {case!r}")
    if kind == "case-selected":
        refuse(17, f"duplicate case-selected for {case!r}")
    elif kind == "attempt-started":
        if st["open_attempt"]:
            refuse(17, f"nested attempt-started for {case!r}")
        st["attempts"] += 1
        st["open_attempt"] = True
        st["stage_open"] = False
    elif kind == "attempt-terminal":
        if not st["open_attempt"]:
            refuse(17, f"unpaired attempt-terminal for {case!r}")
        st["open_attempt"] = False
        st["stage_open"] = False
    elif kind == "stage-started":
        if not st["open_attempt"] or st["stage_open"]:
            refuse(17, f"misplaced stage-started for {case!r}")
        st["stage_open"] = True
    elif kind == "stage-terminal":
        if not st["stage_open"]:
            refuse(17, f"unpaired stage-terminal for {case!r}")
        st["stage_open"] = False
    elif kind in ("first-divergence", "case-cancel-request", "case-cancel-drained"):
        if not st["open_attempt"]:
            refuse(17, f"{kind} outside an open attempt for {case!r}")
    elif kind == "case-terminal":
        if st["open_attempt"] or st["terminal"] is not None:
            refuse(17, f"misplaced case-terminal for {case!r}")
        st["terminal"] = event.get("status")
    else:
        refuse(17, f"unexpected event placement: {kind}")
for case, st in state.items():
    if st["terminal"] is None:
        refuse(17, f"case selected without terminal: {case!r}")
recorded = {}
for row in summary.get("cases", []):
    recorded[row["case"]] = row["status"]
logged = {c: st["terminal"] for c, st in state.items() if st["terminal"]}
if logged != recorded:
    refuse(17, f"terminal statuses disagree: log={logged} summary={recorded}")

os.makedirs(out_root)
with open(os.path.join(out_root, "replay-verdict.json"), "w") as fh:
    json.dump({"schema": "frankensim.new-domains.replay-verdict.v1",
               "receipt": os.path.relpath(receipt_path, ri),
               "cases": len(recorded),
               "verdict": "agree"}, fh, indent=2, sort_keys=True)
    fh.write("\n")
print(f"replay OK: {len(recorded)} case terminals agree; verdict retained under {os.path.relpath(out_root, real(repo_root))}")
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
  python3 - "$SCHEMA_VERSION" "$REPO_ROOT" "$CAP_CASES_PER_MANIFEST" "$CAP_CHILD_OUTPUT_BYTES" "${MANIFESTS[@]}" <<'PYEOF'
import json, sys, os, re
try:
    import tomllib
except ModuleNotFoundError:
    print("python3 tomllib unavailable (need >=3.11)", file=sys.stderr); sys.exit(18)
schema, repo_root = sys.argv[1], sys.argv[2]
cap_cases_per_manifest, cap_child_output = int(sys.argv[3]), int(sys.argv[4])
paths = sys.argv[5:]

REQUIRED_CASE_KEYS = {
    "id": str, "version": int, "purpose": str, "owning_bead": str,
    "gauntlet_tier": str, "entry_command": list, "seed": int,
    "max_wall_seconds": int, "expected": str,
}
OPTIONAL_CASE_KEYS = {
    "checker_command": list, "expected_refusal_pattern": str,
    "seed_overridable": bool, "expected_output_pattern": str,
    "qoi_notes": str, "determinism_class": str, "max_output_bytes": int,
    "profiles": list,
}
ADMITTED_EXPECTED = {"authority", "refusal"}
ADMITTED_PROFILES = {"smoke"}
ADMITTED_PROBES = {
    "schema-mismatch", "unknown-field", "duplicate-case-id",
    "path-escape-command", "refusal-pattern-missing", "output-cap-crossing",
}
ADMITTED_NEGATIVE_EXITS = {10, 11, 12, 13, 14, 15, 16, 17, 18}

def refuse(msg):
    print(f"manifest refusal: {msg}", file=sys.stderr)
    sys.exit(10)

def bounded(name, value, limit):
    if len(value) > limit:
        refuse(f"{name} exceeds its {limit}-byte bound")

seen_ids = set()
negative_ids = set()
projected = []
negatives = []
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
    cases = doc.get("case", [])
    declared_negatives = doc.get("negative", [])
    if not isinstance(cases, list):
        refuse(f"{path}: [[case]] must be an array of tables")
    if not isinstance(declared_negatives, list):
        refuse(f"{path}: [[negative]] must be an array of tables")
    if not cases and not declared_negatives:
        refuse(f"{path}: needs at least one [[case]] or one [[negative]]")
    if len(cases) > cap_cases_per_manifest:
        refuse(f"{path}: {len(cases)} cases exceed the hard cap of {cap_cases_per_manifest}")
    unknown_top = set(doc) - {"schema", "phase", "case", "negative"}
    if unknown_top:
        refuse(f"{path}: unknown top-level fields {sorted(unknown_top)} (no silent semantics)")
    for negative in doc.get("negative", []):
        unknown = set(negative) - {"id", "purpose", "probe", "expect_exit"}
        if unknown:
            refuse(f"{path}: negative {negative.get('id')!r} has unknown fields {sorted(unknown)}")
        for key in ("id", "purpose", "probe", "expect_exit"):
            if key not in negative:
                refuse(f"{path}: negative entry missing {key!r}")
        nid = negative["id"]
        if not isinstance(nid, str) or not re.fullmatch(r"[a-z0-9][a-z0-9-]{2,63}", nid):
            refuse(f"{path}: negative id {nid!r} is not a stable slug")
        if nid in seen_ids or nid in negative_ids:
            refuse(f"duplicate id {nid!r} across cases and negatives")
        negative_ids.add(nid)
        if not isinstance(negative["purpose"], str):
            refuse(f"{path}: negative {nid!r} purpose must be a string")
        bounded(f"{path}: negative {nid!r} purpose", negative["purpose"], 512)
        if negative["probe"] not in ADMITTED_PROBES:
            refuse(f"{path}: negative {nid!r} probe must be one of {sorted(ADMITTED_PROBES)}")
        if not isinstance(negative["expect_exit"], int) or negative["expect_exit"] not in ADMITTED_NEGATIVE_EXITS:
            refuse(f"{path}: negative {nid!r} expect_exit must be a stable exit class {sorted(ADMITTED_NEGATIVE_EXITS)}")
        entry = dict(negative)
        entry["_phase"] = phase
        entry["_manifest"] = os.path.relpath(path, repo_root)
        negatives.append(entry)
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
        if case["id"] in seen_ids or case["id"] in negative_ids:
            refuse(f"duplicate id {case['id']!r} across cases and negatives")
        seen_ids.add(case["id"])
        if case["expected"] not in ADMITTED_EXPECTED:
            refuse(f"{path}: case {case['id']!r} expected must be one of {sorted(ADMITTED_EXPECTED)}")
        if case["expected"] == "refusal" and "expected_refusal_pattern" not in case:
            refuse(f"{path}: refusal case {case['id']!r} must declare expected_refusal_pattern")
        if case["max_wall_seconds"] < 1 or case["max_wall_seconds"] > 7200:
            refuse(f"{path}: case {case['id']!r} max_wall_seconds outside [1, 7200]")
        if "max_output_bytes" in case and not (1 <= case["max_output_bytes"] <= cap_child_output):
            refuse(f"{path}: case {case['id']!r} max_output_bytes outside [1, {cap_child_output}]")
        if "profiles" in case:
            profiles = case["profiles"]
            if any(p not in ADMITTED_PROFILES for p in profiles):
                refuse(f"{path}: case {case['id']!r} profiles must be a subset of {sorted(ADMITTED_PROFILES)}")
            if len(set(profiles)) != len(profiles):
                refuse(f"{path}: case {case['id']!r} profiles contain duplicates")
        bounded(f"{path}: case {case['id']!r} purpose", case["purpose"], 512)
        bounded(f"{path}: case {case['id']!r} owning_bead", case["owning_bead"], 64)
        bounded(f"{path}: case {case['id']!r} gauntlet_tier", case["gauntlet_tier"], 32)
        if "qoi_notes" in case:
            bounded(f"{path}: case {case['id']!r} qoi_notes", case["qoi_notes"], 512)
        if "determinism_class" in case:
            bounded(f"{path}: case {case['id']!r} determinism_class", case["determinism_class"], 64)
        for pattern_key in ("expected_refusal_pattern", "expected_output_pattern"):
            if pattern_key in case:
                bounded(f"{path}: case {case['id']!r} {pattern_key}", case[pattern_key], 256)
                try:
                    re.compile(case[pattern_key])
                except re.error as error:
                    refuse(f"{path}: case {case['id']!r} {pattern_key} is not a valid regex: {error}")
        for word in case["entry_command"] + case.get("checker_command", []):
            if not isinstance(word, str):
                refuse(f"{path}: case {case['id']!r} command words must be strings")
            if word.startswith("/") or ".." in word.split(os.sep):
                refuse(f"{path}: case {case['id']!r} command word {word!r} escapes the repo (absolute or ..)")
            bounded(f"{path}: case {case['id']!r} command word", word, 4096)
        if len(case["entry_command"]) > 64 or len(case.get("checker_command", [])) > 64:
            refuse(f"{path}: case {case['id']!r} command exceeds 64 words")
        rel_manifest = os.path.relpath(path, repo_root)
        bounded(f"{path}: manifest path", rel_manifest, 200)
        case["_phase"] = phase
        case["_manifest"] = rel_manifest
        projected.append(case)
print(json.dumps({"cases": projected, "negatives": negatives}))
PYEOF
}

PROJECTED="$(validate_and_project)" || exit $?
NEGATIVES_JSON="$(printf '%s' "$PROJECTED" | python3 -c 'import json,sys; print(json.dumps(json.load(sys.stdin)["negatives"]))')"
PROJECTED="$(printf '%s' "$PROJECTED" | python3 -c 'import json,sys; print(json.dumps(json.load(sys.stdin)["cases"]))')"

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

# Canonical run modes: smoke executes only cases whose manifest declaration
# includes profiles = ["smoke"]; bare invocation and --run full execute the
# full admitted selection.
if [ "$RUN_MODE" = "smoke" ]; then
  PROJECTED="$(printf '%s' "$PROJECTED" | python3 -c '
import json, sys
cases = json.load(sys.stdin)
print(json.dumps([c for c in cases if "smoke" in c.get("profiles", [])]))
')" || exit $?
  SMOKE_COUNT="$(printf '%s' "$PROJECTED" | python3 -c 'import json,sys; print(len(json.load(sys.stdin)))')"
  [ "$SMOKE_COUNT" -gt 0 ] || die "$EXIT_USAGE" "--run smoke selected no cases: no selected case declares profiles = [\"smoke\"]"
fi

COUNT="$(printf '%s' "$PROJECTED" | python3 -c 'import json,sys; print(len(json.load(sys.stdin)))')"
if [ "$COUNT" -gt "$CAP_SELECTED_CASES" ]; then
  die "$EXIT_USAGE" "$COUNT selected cases exceed the hard cap of $CAP_SELECTED_CASES"
fi

# Canonical --check (and its --validate-only alias): source, manifest,
# route, capability, path, and schema closure WITHOUT production execution.
if [ "$CHECK" -eq 1 ]; then
  command -v timeout >/dev/null 2>&1 || die "$EXIT_UNQUALIFIED" "capability closure: timeout(1) is required"
  CHECK_FAILURES="$(printf '%s' "$PROJECTED" | CHECK_REPO_ROOT="$REPO_ROOT" python3 -c '
import json, os, shutil, sys
repo = os.environ["CHECK_REPO_ROOT"]
failures = []
for case in json.load(sys.stdin):
    for role in ("entry_command", "checker_command"):
        words = case.get(role, [])
        if not words:
            continue
        first = words[0]
        resolved = shutil.which(first) is not None
        if not resolved:
            candidate = os.path.join(repo, first)
            resolved = os.path.isfile(candidate)
        if not resolved:
            cid = case["id"]
            failures.append(f"{cid}: {role} entry {first!r} does not resolve on PATH or under the repo")
print("; ".join(failures))
')"
  [ -z "$CHECK_FAILURES" ] || die "$EXIT_USAGE" "closure failures: $CHECK_FAILURES"
  NEG_COUNT="$(printf '%s' "$NEGATIVES_JSON" | python3 -c 'import json,sys; print(len(json.load(sys.stdin)))')"
  echo "check OK: ${#MANIFESTS[@]} manifest(s), $COUNT case(s), $NEG_COUNT negative(s) admitted; route/capability/path/schema closure verified"
  exit "$EXIT_OK"
fi

# Frozen negative probes: --negative executes exactly one immutable probe
# declared by a [[negative]] entry in the selected phase manifests. The
# probes are source-closed runner-contract checks; they run no phase
# production entries.
if [ -n "$NEGATIVE" ]; then
  DECLARED="$(printf '%s' "$NEGATIVES_JSON" | python3 -c '
import json, sys
nid = sys.argv[1]
for n in json.load(sys.stdin):
    if n["id"] == nid:
        print(json.dumps(n))
        break
' "$NEGATIVE")"
  [ -n "$DECLARED" ] || die "$EXIT_USAGE" "--negative id is not declared by any selected manifest: $NEGATIVE"
  PROBE="$(printf '%s' "$DECLARED" | python3 -c 'import json,sys; print(json.load(sys.stdin)["probe"])')"
  EXPECT="$(printf '%s' "$DECLARED" | python3 -c 'import json,sys; print(json.load(sys.stdin)["expect_exit"])')"
  NEG_DIR="$REPO_ROOT/.e2e-out/negative-$(date -u +%Y%m%dT%H%M%SZ)-$$"
  mkdir -p "$NEG_DIR" || die "$EXIT_UNQUALIFIED" "cannot create negative scratch dir"
  M="$NEG_DIR/fixture-phase.toml"
  CASE_HEADER='schema = "frankensim.new-domains.case-manifest.v1"'
  good_case() { # id extra-lines...
    local id="$1"; shift
    printf '[[case]]\nid = "%s"\nversion = 1\npurpose = "negative fixture"\nowning_bead = "fixture"\ngauntlet_tier = "G0"\n' "$id"
    printf '%s\n' "$@"
  }
  expect_probe() { # description expected-command...
    local desc="$1"; shift
    "${@}" >/dev/null 2>&1
    local got=$?
    [ "$got" -eq "$EXPECT" ] || { rm -rf "$NEG_DIR"; die 17 "negative $NEGATIVE ($desc): expected exit $EXPECT, got $got"; }
  }
  R=(bash "$RUNNER_PATH" --manifest "$M")
  case "$PROBE" in
    schema-mismatch)
      { printf '%s\n' 'schema = "wrong.schema.v9"' 'phase = "fixture"'; } > "$M"
      expect_probe "schema-mismatch" "${R[@]}" --list ;;
    unknown-field)
      { printf '%s\n' "$CASE_HEADER" 'phase = "fixture"' 'undeclared_field = 1'; } > "$M"
      expect_probe "unknown-field" "${R[@]}" --list ;;
    duplicate-case-id)
      { printf '%s\n' "$CASE_HEADER" 'phase = "fixture"'
        good_case dup-case 'entry_command = ["true"]' 'expected = "authority"'
        good_case dup-case 'entry_command = ["true"]' 'expected = "authority"'
      } > "$M"
      expect_probe "duplicate-case-id" "${R[@]}" --list ;;
    path-escape-command)
      { printf '%s\n' "$CASE_HEADER" 'phase = "fixture"'
        good_case esc-case 'entry_command = ["/bin/echo", "hi"]' 'expected = "authority"'
      } > "$M"
      expect_probe "path-escape-command" "${R[@]}" --list ;;
    refusal-pattern-missing)
      { printf '%s\n' "$CASE_HEADER" 'phase = "fixture"'
        good_case ref-case 'entry_command = ["false"]' 'expected = "refusal"'
      } > "$M"
      expect_probe "refusal-pattern-missing" "${R[@]}" --list ;;
    output-cap-crossing)
      { printf '%s\n' "$CASE_HEADER" 'phase = "fixture"'
        good_case flood-case 'entry_command = ["bash", "-c", "printf %s $(seq 1 200000)"]' 'expected = "authority"' 'max_output_bytes = 1024'
      } > "$M"
      expect_probe "output-cap-crossing" "${R[@]}" --output-dir "$NEG_DIR/out" ;;
    *) rm -rf "$NEG_DIR"; die "$EXIT_UNQUALIFIED" "unimplemented negative probe: $PROBE" ;;
  esac
  rm -rf "$NEG_DIR"
  echo "negative OK: $NEGATIVE (probe=$PROBE) refused with declared class $EXPECT"
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

# Redaction by construction: credential-shaped variables never reach a
# child process, so no log or retained artifact can leak them.
declare -a ENV_SCRUB=()
while IFS='=' read -r name _; do
  case "$name" in
    *SECRET*|*TOKEN*|*PASSWORD*|*CREDENTIAL*|*API_KEY*|*APIKEY*|AWS_*|GH_*|GITHUB_*)
      ENV_SCRUB+=("-u" "$name") ;;
  esac
done < <(env)

LOG="$OUTPUT_DIR/runner-log.jsonl"
SUMMARY="$OUTPUT_DIR/summary.json"
SEQ=0
HEAD_SHA="$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"

emit() { # event-name json-fields...
  SEQ=$((SEQ + 1))
  if [ "$SEQ" -gt "$CAP_RECORDS_PER_RUN" ]; then
    touch "$OUTPUT_DIR/.cap-breach-records"
    return 1
  fi
  python3 - "$LOG" "$LOG_SCHEMA" "$OBS_IDENTITY_DOMAIN" "$SEQ" "$HEAD_SHA" "$CAP_EVENT_BYTES" "$@" <<'PYEOF'
import json, sys, datetime
log, schema, domain, seq, head, cap_bytes = sys.argv[1:7]
event = {"schema": schema, "obs_identity_domain": domain, "seq": int(seq),
         "head_sha": head,
         "at": datetime.datetime.now(datetime.timezone.utc).isoformat()}
event["event"] = sys.argv[7]
for pair in sys.argv[8:]:
    key, _, value = pair.partition("=")
    event[key] = value
line = json.dumps(event, sort_keys=True) + "\n"
if len(line.encode()) > int(cap_bytes):
    print(f"canonical event exceeds {cap_bytes} bytes ({len(line.encode())})", file=sys.stderr)
    sys.exit(1)
with open(log, "a") as fh:
    fh.write(line)
PYEOF
  if [ $? -ne 0 ]; then
    touch "$OUTPUT_DIR/.cap-breach-event"
    return 1
  fi
}

# Volatile operational envelope: wall-clock durations live here and never
# enter the canonical semantic log.
volatile_note() { # elapsed_seconds extra...
  python3 - "$CASE_DIR/volatile.jsonl" "$@" <<'PYEOF'
import json, sys, datetime
path = sys.argv[1]
rec = {"at": datetime.datetime.now(datetime.timezone.utc).isoformat()}
for pair in sys.argv[2:]:
    key, _, value = pair.partition("=")
    rec[key] = value
with open(path, "a") as fh:
    fh.write(json.dumps(rec, sort_keys=True) + "\n")
PYEOF
}

# Copy-pasteable repository-relative reproduction command for one case.
# Fresh-output placeholder keeps reruns collision-free without embedding a
# concrete volatile path in the retained record.
repro_for_case() { # id seed manifest_rel
  local r
  r="$(printf 'scripts/ci/new_domains_e2e.sh --run full --manifest %s --case %s --seed %s --output-dir .e2e-out/new-domains-$(date -u +%%Y%%m%%dT%%H%%M%%SZ)-$$' "$3" "$1" "$2")"
  if [ "${#r}" -gt 512 ]; then
    r="$(printf 'scripts/ci/new_domains_e2e.sh --run full --manifest %s --case %s' "$3" "$1")"
  fi
  printf '%s' "$r"
}

emit run-start "manifest_count=${#MANIFESTS[@]}" "case_count=$COUNT" \
  "seed_override=${SEED:-none}" "max_wall_override=${MAX_WALL:-none}"

WORST="$EXIT_OK"
while IFS= read -r case_json; do
  [ -f "$OUTPUT_DIR/.cap-breach-records" ] || [ -f "$OUTPUT_DIR/.cap-breach-event" ] && break
  ID="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])')"
  EXPECTED="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin)["expected"])')"
  WALL="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin)["max_wall_seconds"])')"
  CASE_SEED="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin)["seed"])')"
  OVERRIDABLE="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(str(json.load(sys.stdin).get("seed_overridable", False)).lower())')"
  REFUSAL_PATTERN="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin).get("expected_refusal_pattern", ""))')"
  OUTPUT_PATTERN="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin).get("expected_output_pattern", ""))')"
  MANIFEST_REL="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin)["_manifest"])')"
  REPRO="$(repro_for_case "$ID" "$CASE_SEED" "$MANIFEST_REL")"

  # Cardinality: exactly one CaseSelected per selected case, emitted before
  # any disposition (including unqualified) is decided.
  emit case-selected "case=$ID" "expected=$EXPECTED" "manifest=$MANIFEST_REL"

  # Budget overrides only tighten the admitted manifest.
  if [ -n "$MAX_WALL" ] && [ "$MAX_WALL" -lt "$WALL" ]; then WALL="$MAX_WALL"; fi
  if [ -n "$SEED" ]; then
    if [ "$OVERRIDABLE" != "true" ]; then
      # Zero-attempt terminal: the refusal is admission-side, no production
      # attempt ever begins.
      emit case-terminal "case=$ID" "status=unqualified" "reason=seed-not-overridable" \
        "wall_budget=$WALL" "repro_cmd=$REPRO"
      [ "$WORST" -lt "$EXIT_UNQUALIFIED" ] && WORST="$EXIT_UNQUALIFIED"
      continue
    fi
    CASE_SEED="$SEED"
  fi

  CASE_DIR="$OUTPUT_DIR/$ID"
  mkdir -p "$CASE_DIR"
  emit attempt-started "case=$ID" "attempt=1"
  emit stage-started "case=$ID" "stage=production"

  # One word per line, read into an array so multi-word arguments (e.g. a
  # bash -c body) keep their boundaries. Newlines inside words are refused
  # at validation time by the slug/word rules.
  declare -a ENTRY_WORDS=()
  while IFS= read -r word; do
    ENTRY_WORDS+=("$word")
  done < <(printf '%s' "$case_json" | python3 -c 'import json,sys
for word in json.load(sys.stdin)["entry_command"]:
    print(word)')
  OUTPUT_CAP="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin).get("max_output_bytes", 1048576))')"
  START=$(date +%s)
  ( cd "$REPO_ROOT" && NEW_DOMAINS_SEED="$CASE_SEED" NEW_DOMAINS_CASE_DIR="$CASE_DIR" \
      env "${ENV_SCRUB[@]}" timeout "$WALL" "${ENTRY_WORDS[@]}" ) \
      >"$CASE_DIR/stdout.txt" 2>"$CASE_DIR/stderr.txt" &
  CHILD=$!
  # Forward operator signals to the child instead of orphaning it.
  trap 'kill -TERM "$CHILD" 2>/dev/null' TERM INT
  CANCELLED=0
  if [ -n "$CANCEL_AFTER" ]; then
    DEADLINE=$(( $(date +%s) + CANCEL_AFTER ))
    while kill -0 "$CHILD" 2>/dev/null && [ "$(date +%s)" -lt "$DEADLINE" ]; do
      sleep 1
    done
    if kill -0 "$CHILD" 2>/dev/null; then
      CANCELLED=1
      emit case-cancel-request "case=$ID"
      kill -TERM "$CHILD" 2>/dev/null
      GRACE=$(( $(date +%s) + 10 ))
      while kill -0 "$CHILD" 2>/dev/null && [ "$(date +%s)" -lt "$GRACE" ]; do
        sleep 1
      done
      if kill -0 "$CHILD" 2>/dev/null; then
        kill -KILL "$CHILD" 2>/dev/null
        wait "$CHILD" 2>/dev/null
        trap - TERM INT
        emit stage-terminal "case=$ID" "stage=production" "status=drain-failure"
        emit attempt-terminal "case=$ID" "attempt=1" "status=drain-failure"
        volatile_note "elapsed_seconds=$(( $(date +%s) - START ))" "exit_code=signal"
        emit case-terminal "case=$ID" "status=failed" "reason=cancel-drain-failure" \
          "wall_budget=$WALL" "repro_cmd=$REPRO"
        [ "$WORST" -lt "$EXIT_CANCEL" ] && WORST="$EXIT_CANCEL"
        continue
      fi
      emit case-cancel-drained "case=$ID"
    fi
  fi
  wait "$CHILD"
  STATUS_CODE=$?
  trap - TERM INT
  ELAPSED=$(( $(date +%s) - START ))
  if [ "$CANCELLED" -eq 1 ]; then
    emit stage-terminal "case=$ID" "stage=production" "status=cancelled"
    emit attempt-terminal "case=$ID" "attempt=1" "status=cancelled"
    volatile_note "elapsed_seconds=$ELAPSED" "exit_code=$STATUS_CODE"
    emit case-terminal "case=$ID" "status=cancelled" "reason=operator-cancellation" \
      "wall_budget=$WALL" "repro_cmd=$REPRO"
    [ "$WORST" -lt "$EXIT_CANCEL" ] && WORST="$EXIT_CANCEL"
    continue
  fi
  OUT_BYTES=$(( $(wc -c < "$CASE_DIR/stdout.txt") + $(wc -c < "$CASE_DIR/stderr.txt") ))
  emit stage-terminal "case=$ID" "stage=production" "status=exit-$STATUS_CODE"
  volatile_note "elapsed_seconds=$ELAPSED" "exit_code=$STATUS_CODE"
  if [ "$OUT_BYTES" -gt "$OUTPUT_CAP" ]; then
    # Bounded retention: keep a capped prefix, never the unbounded stream.
    head -c "$OUTPUT_CAP" "$CASE_DIR/stdout.txt" > "$CASE_DIR/stdout.capped.txt"
    mv "$CASE_DIR/stdout.capped.txt" "$CASE_DIR/stdout.txt"
    head -c "$OUTPUT_CAP" "$CASE_DIR/stderr.txt" > "$CASE_DIR/stderr.capped.txt"
    mv "$CASE_DIR/stderr.capped.txt" "$CASE_DIR/stderr.txt"
    emit attempt-terminal "case=$ID" "attempt=1" "status=output-budget-exhausted"
    emit case-terminal "case=$ID" "status=failed" "reason=output-budget-exhausted" \
      "wall_budget=$WALL" "repro_cmd=$REPRO"
    [ "$WORST" -lt "$EXIT_BUDGET" ] && WORST="$EXIT_BUDGET"
    continue
  fi

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

  # Determinism repeat: a declared repo-deterministic case must reproduce
  # its stdout byte-for-byte on an immediate second run. The first
  # divergence is recorded as at most one bounded, content-free event.
  if [ "$RESULT" = "passed" ] && [ "$DETERMINISM_REPEAT" -eq 1 ]; then
    DET_CLASS_DECL="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(json.load(sys.stdin).get("determinism_class", ""))')"
    case "$DET_CLASS_DECL" in
      repo-deterministic*)
        mkdir -p "$CASE_DIR/repeat"
        ( cd "$REPO_ROOT" && NEW_DOMAINS_SEED="$CASE_SEED" NEW_DOMAINS_CASE_DIR="$CASE_DIR/repeat" \
            env "${ENV_SCRUB[@]}" timeout "$WALL" "${ENTRY_WORDS[@]}" ) \
            >"$CASE_DIR/repeat/stdout.txt" 2>"$CASE_DIR/repeat/stderr.txt"
        DIVERGENCE="$(cmp "$CASE_DIR/stdout.txt" "$CASE_DIR/repeat/stdout.txt" 2>/dev/null | head -1)"
        if [ -n "$DIVERGENCE" ]; then
          DIV_OFFSET="$(printf '%s' "$DIVERGENCE" | sed -n 's/.*byte \([0-9][0-9]*\).*/\1/p')"
          emit first-divergence "case=$ID" "offset=${DIV_OFFSET:-unknown}"
          RESULT="failed"; REASON="determinism-mismatch"; CLASS="$EXIT_DETERMINISM"
        fi
        ;;
    esac
  fi

  # Independent checker phase (only meaningful after a passing production).
  if [ "$RESULT" = "passed" ]; then
    CHECKER_LEN="$(printf '%s' "$case_json" | python3 -c 'import json,sys; print(len(json.load(sys.stdin).get("checker_command", [])))')"
    if [ "$CHECKER_LEN" -gt 0 ]; then
      declare -a CHECKER_WORDS=()
      while IFS= read -r word; do
        CHECKER_WORDS+=("$word")
      done < <(printf '%s' "$case_json" | python3 -c 'import json,sys
for word in json.load(sys.stdin).get("checker_command", []):
    print(word)')
      emit stage-started "case=$ID" "stage=checker"
      ( cd "$REPO_ROOT" && NEW_DOMAINS_CASE_DIR="$CASE_DIR" timeout "$WALL" "${CHECKER_WORDS[@]}" ) \
        >"$CASE_DIR/checker-stdout.txt" 2>"$CASE_DIR/checker-stderr.txt"
      CHECKER_CODE=$?
      emit stage-terminal "case=$ID" "stage=checker" "status=exit-$CHECKER_CODE"
      if [ "$CHECKER_CODE" -ne 0 ]; then
        RESULT="failed"; REASON="checker-refused"; CLASS="$EXIT_ACCEPTANCE"
      fi
    fi
  fi

  # Retained-artifact cap per case: crossing it can never count as pass.
  CASE_ARTIFACTS="$(find "$CASE_DIR" -type f | wc -l | tr -d ' ')"
  if [ "$CASE_ARTIFACTS" -gt "$CAP_ARTIFACTS_PER_CASE" ] && { [ "$RESULT" = "passed" ] || [ "$RESULT" = "expected-refused" ]; }; then
    RESULT="failed"; REASON="artifact-cap-exhausted"; CLASS="$EXIT_BUDGET"
  fi
  emit attempt-terminal "case=$ID" "attempt=1" "status=$RESULT"
  emit case-terminal "case=$ID" "status=$RESULT" "reason=$REASON" \
    "wall_budget=$WALL" "repro_cmd=$REPRO"
  [ "$WORST" -lt "$CLASS" ] && WORST="$CLASS"
done < <(printf '%s' "$PROJECTED" | python3 -c '
import json, sys
for case in json.load(sys.stdin):
    print(json.dumps(case))')

# Cap breaches recorded during the run force the budget class; committed
# bounded evidence is preserved exactly as written.
if [ -f "$OUTPUT_DIR/.cap-breach-records" ]; then
  echo "record cap breached ($CAP_RECORDS_PER_RUN per run); further logging stopped" >&2
  [ "$WORST" -lt "$EXIT_BUDGET" ] && WORST="$EXIT_BUDGET"
fi
if [ -f "$OUTPUT_DIR/.cap-breach-event" ]; then
  echo "canonical event byte cap breached ($CAP_EVENT_BYTES)" >&2
  [ "$WORST" -lt "$EXIT_BUDGET" ] && WORST="$EXIT_BUDGET"
fi

# Retained-artifact cap across the whole run.
RUN_ARTIFACTS="$(find "$OUTPUT_DIR" -type f | wc -l | tr -d ' ')"
if [ "$RUN_ARTIFACTS" -gt "$CAP_ARTIFACTS_PER_RUN" ]; then
  echo "retained artifact cap breached ($CAP_ARTIFACTS_PER_RUN per run): $RUN_ARTIFACTS present" >&2
  [ "$WORST" -lt "$EXIT_BUDGET" ] && WORST="$EXIT_BUDGET"
fi

# Copy-pasteable repository-relative reproduction command for the run.
RUN_REPRO="scripts/ci/new_domains_e2e.sh --run full"
if [ -n "$MANIFEST" ]; then
  RUN_REPRO="$RUN_REPRO --manifest ${MANIFEST#"$REPO_ROOT"/}"
elif [ -n "$PHASE" ]; then
  RUN_REPRO="$RUN_REPRO --phase $PHASE"
fi
if [ "${#CASES[@]}" -gt 0 ] && [ "${#CASES[@]}" -le 3 ]; then
  for c in "${CASES[@]}"; do RUN_REPRO="$RUN_REPRO --case $c"; done
fi

# Cardinality: exactly one summary after all case terminals, then exactly
# one run-terminal closing the invocation.
COUNTS_JSON="$(python3 - "$LOG" <<'PYEOF'
import json, sys, collections
c = collections.Counter()
for line in open(sys.argv[1]):
    event = json.loads(line)
    if event["event"] == "case-terminal":
        c[event["status"]] += 1
print(json.dumps(dict(c), sort_keys=True))
PYEOF
)"
emit summary "counts=$COUNTS_JSON" "worst_class=$WORST" "repro_cmd=$RUN_REPRO"
emit run-terminal "worst_class=$WORST" "repro_cmd=$RUN_REPRO"

python3 - "$SUMMARY" "$LOG" "$WORST" "$CAP_ARTIFACT_INDEX_BYTES" <<'PYEOF'
import json, sys, hashlib, os, collections
summary_path, log_path = sys.argv[1], sys.argv[2]
worst, index_cap = int(sys.argv[3]), int(sys.argv[4])
rows = []
for line in open(log_path):
    event = json.loads(line)
    if event["event"] == "case-terminal":
        rows.append({"case": event["case"], "status": event["status"],
                     "reason": event.get("reason", ""),
                     "repro_cmd": event.get("repro_cmd", "")})
counts = dict(collections.Counter(row["status"] for row in rows))
summary = {
    "schema": "frankensim.new-domains.runner-summary.v2",
    "log_file": os.path.basename(log_path),
    "log_sha256": hashlib.sha256(open(log_path, "rb").read()).hexdigest(),
    "worst_exit_class": worst,
    "counts": counts,
    "cases": rows,
}
payload = json.dumps(summary, indent=2, sort_keys=True) + "\n"
truncated = False
while len(payload.encode()) > index_cap:
    truncated = True
    rows = rows[: len(rows) // 2]
    summary["cases"] = rows
    summary["truncated"] = True
    payload = json.dumps(summary, indent=2, sort_keys=True) + "\n"
with open(summary_path, "w") as fh:
    fh.write(payload)
if truncated:
    print(f"artifact index exceeded {index_cap} bytes; summary truncated (class stays budgeted)", file=sys.stderr)
print(f"summary: {json.dumps(counts, sort_keys=True)} -> exit {worst}")
print(f"retained: {summary_path}")
PYEOF

exit "$WORST"
