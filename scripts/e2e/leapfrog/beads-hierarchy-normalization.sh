#!/usr/bin/env bash
# Beads rollup-relation normalization tool
# (bead frankensim-leapfrog-2026-program-i94v.7.1.8).
#
# Corrects a graph-semantics defect - epics encoding rollups as
# epic-depends-on-dot-prefix-child BLOCKS edges - without losing or
# reordering any implementation scope. The controlling scope is a frozen,
# content-addressed review manifest derived from the exact live br DB
# identity; neither title prose nor a stale bv snapshot supplies expected
# membership, and every candidate row requires semantic review because some
# dot-prefix edges are genuine functional prerequisites.
#
# Modes (this revision is READ-ONLY: --apply intentionally does not exist
# yet and will refuse-by-absence; it lands only after reviewed rows exist):
#   --list             enumerate live candidate edges, run nothing
#   --freeze [OUT]     write the content-addressed review manifest (all
#                      rows verdict=pending) to OUT, defaulting to the
#                      tracked tests/leapfrog/manifests location
#   --check MANIFEST   validate a manifest's schema, counts, and DB binding
#   --plan MANIFEST    emit the exact two-step br command sequence for rows
#                      reviewed verdict=migrate, plus the inverse plan;
#                      refuses if ANY row is still pending
#   --self-test        exercise the tool's own refusals on fixtures
#
# Migration shape (never performed here): for each reviewed relation,
#   br dep remove <epic> <child>            # drop the blocks edge FIRST
#   br dep add <child> <epic> --type parent-child
# Never add before remove; never write JSONL or SQLite directly.
#
# No-claim boundary: this corrects issue-tracker hierarchy and planning
# truth only. It implements no product capability, closes no child, and
# does not make bv output the actionability source of truth; br ready
# remains canonical.
set -u -o pipefail

SCHEMA="frankensim.beads-hierarchy-normalization-manifest.v1"
REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BEADS_FILE="${FSIM_HIERNORM_BEADS:-$REPO_ROOT/.beads/issues.jsonl}"
DEFAULT_MANIFEST="$REPO_ROOT/tests/leapfrog/manifests/hierarchy-normalization-manifest.json"

EXIT_OK=0
EXIT_USAGE=30
EXIT_STALE=31
EXIT_UNREVIEWED=32

die() {
  local class="$1"; shift
  printf 'hierarchy-normalization: ERROR class=%s: %s\n' "$class" "$*" >&2
  exit "$class"
}
command -v python3 >/dev/null 2>&1 || die "$EXIT_USAGE" "python3 is required"

MODE="" MANIFEST_ARG=""
case "${1:-}" in
  --list) MODE="list" ;;
  --freeze) MODE="freeze"; MANIFEST_ARG="${2:-}" ;;
  --self-test) MODE="self-test" ;;
  --check|--plan)
    MODE="${1#--}"
    [ $# -ge 2 ] || die "$EXIT_USAGE" "$1 needs a manifest path"
    MANIFEST_ARG="$2" ;;
  --apply)
    die "$EXIT_USAGE" "--apply does not exist in this revision by design: it lands only after the frozen manifest carries reviewed verdicts (see the bead)" ;;
  *) die "$EXIT_USAGE" "usage: --list | --freeze | --check MANIFEST | --plan MANIFEST | --self-test" ;;
esac

tool() { # subcommand [manifest]
  python3 - "$1" "$BEADS_FILE" "$SCHEMA" "${2:-}" <<'PYEOF'
import hashlib
import json
import sys

mode, beads_path, schema, manifest_path = sys.argv[1:5]

EXIT_USAGE = 30
EXIT_STALE = 31
EXIT_UNREVIEWED = 32


def load(path):
    issues, deps = {}, []
    with open(path) as fh:
        text = fh.read()
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        if not isinstance(row, dict):
            continue
        if "issue_id" in row and "depends_on_id" in row:
            deps.append(row)
        elif "id" in row:
            issues[row["id"]] = row
            for dep in row.get("dependencies") or []:
                if isinstance(dep, dict) and "depends_on_id" in dep:
                    deps.append({"issue_id": row["id"], **dep})
    return text, issues, deps


def candidates(issues, deps):
    # An epic-side rollup encoded as blocks: the EPIC depends on its own
    # dot-prefix CHILD. (child depends-on parent edges are a different,
    # legitimate ordering pattern and are not candidates.)
    rows = []
    seen = set()
    for dep in deps:
        if dep.get("type") != "blocks":
            continue
        parent, child = dep["issue_id"], dep["depends_on_id"]
        if not child.startswith(parent + "."):
            continue
        key = (parent, child)
        if key in seen:
            continue
        seen.add(key)
        parent_row = issues.get(parent, {})
        child_row = issues.get(child, {})
        digest = hashlib.sha256(
            (
                (child_row.get("description") or "")
                + "\x00"
                + (child_row.get("acceptance_criteria") or "")
            ).encode()
        ).hexdigest()
        rows.append(
            {
                "parent": parent,
                "child": child,
                "old_type": "blocks",
                "intended_type": "parent-child",
                "parent_status": parent_row.get("status", "missing"),
                "child_status": child_row.get("status", "missing"),
                "child_priority": child_row.get("priority"),
                "child_issue_type": child_row.get("issue_type", child_row.get("type")),
                # Lossless binding: the beads DB remains the text authority;
                # this hash pins the child's exact description+acceptance
                # bytes at freeze time so review verdicts cannot silently
                # apply to different text.
                "child_text_sha256": digest,
                "verdict": "pending",
                "review_note": "",
            }
        )
    rows.sort(key=lambda row: (row["parent"], row["child"]))
    return rows


text, issues, deps = load(beads_path)
rows = candidates(issues, deps)


def candidate_set_digest(rows):
    """Binding scoped to what review verdicts actually depend on: the
    candidate edge set and each child's exact text. Unrelated churn
    (comments, other issues, other edges) must NOT invalidate a frozen
    review; a new/removed candidate edge or edited candidate text MUST."""
    triples = sorted(
        (row["parent"], row["child"], row["child_text_sha256"]) for row in rows
    )
    return hashlib.sha256(json.dumps(triples).encode()).hexdigest()

if mode == "list":
    for row in rows:
        print(
            "\t".join(
                [
                    row["parent"],
                    row["child"],
                    row["parent_status"],
                    row["child_status"],
                ]
            )
        )
    print(f"# {len(rows)} candidate edge(s)", file=sys.stderr)
    sys.exit(0)

if mode == "freeze":
    manifest = {
        "schema": schema,
        "authority_statement": (
            "REVIEW MANIFEST, not migration authority until every row "
            "carries a non-pending verdict. The beads DB remains the sole "
            "clause and hierarchy authority; rows bind child text by hash."
        ),
        "candidate_set_sha256": candidate_set_digest(rows),
        "expected_parent_count": len({row["parent"] for row in rows}),
        "expected_relation_count": len(rows),
        "rows": rows,
    }
    print(json.dumps(manifest, indent=2, sort_keys=True))
    sys.exit(0)

# check / plan need the manifest.
try:
    manifest = json.load(open(manifest_path))
except (OSError, json.JSONDecodeError) as error:
    print(f"manifest refusal: {error}", file=sys.stderr)
    sys.exit(EXIT_USAGE)
if manifest.get("schema") != schema:
    print(f"manifest refusal: schema {manifest.get('schema')!r}", file=sys.stderr)
    sys.exit(EXIT_USAGE)
mrows = manifest.get("rows")
if not isinstance(mrows, list) or not mrows:
    print("manifest refusal: rows missing", file=sys.stderr)
    sys.exit(EXIT_USAGE)
declared_parents = manifest.get("expected_parent_count")
declared_relations = manifest.get("expected_relation_count")
if declared_relations != len(mrows) or declared_parents != len(
    {row["parent"] for row in mrows}
):
    print(
        "manifest refusal: header counts disagree with the row set "
        f"(declared {declared_parents}/{declared_relations}, actual "
        f"{len({row['parent'] for row in mrows})}/{len(mrows)})",
        file=sys.stderr,
    )
    sys.exit(EXIT_USAGE)

if mode == "check":
    live_digest = candidate_set_digest(rows)
    if manifest.get("candidate_set_sha256") != live_digest:
        print(
            "manifest is STALE against the live candidate set: frozen "
            f"{manifest.get('candidate_set_sha256', '')[:16]} vs live {live_digest[:16]}; "
            "a candidate edge appeared/vanished or a candidate child's text "
            "moved - re-freeze and re-review the delta deliberately "
            "(unrelated DB churn does not trip this)",
            file=sys.stderr,
        )
        sys.exit(EXIT_STALE)
    verdicts = {}
    for row in mrows:
        verdicts[row.get("verdict", "?")] = verdicts.get(row.get("verdict", "?"), 0) + 1
    print(f"manifest OK: {len(mrows)} rows, verdicts {json.dumps(verdicts, sort_keys=True)}")
    sys.exit(0)

if mode == "plan":
    admitted = {"pending", "migrate", "keep-functional-prerequisite"}
    bad = [row for row in mrows if row.get("verdict") not in admitted]
    if bad:
        print(
            f"plan refusal: {len(bad)} row(s) carry verdicts outside {sorted(admitted)}",
            file=sys.stderr,
        )
        sys.exit(EXIT_USAGE)
    pending = [row for row in mrows if row.get("verdict") == "pending"]
    if pending:
        print(
            f"plan refusal: {len(pending)} row(s) still pending semantic review; "
            "a plan over an unreviewed manifest would launder the review step",
            file=sys.stderr,
        )
        sys.exit(EXIT_UNREVIEWED)
    live = {(row["parent"], row["child"]) for row in rows}
    plan, inverse = [], []
    for row in mrows:
        if row["verdict"] != "migrate":
            continue
        if (row["parent"], row["child"]) not in live:
            print(
                f"plan refusal: reviewed edge {row['parent']} -> {row['child']} "
                "no longer exists live; re-freeze",
                file=sys.stderr,
            )
            sys.exit(EXIT_STALE)
        # Remove FIRST, then add - never both edges at once.
        plan.append(["br", "dep", "remove", row["parent"], row["child"]])
        plan.append(
            ["br", "dep", "add", row["child"], row["parent"], "--type", "parent-child"]
        )
        inverse.append(["br", "dep", "remove", row["child"], row["parent"]])
        inverse.append(["br", "dep", "add", row["parent"], row["child"], "--type", "blocks"])
    print(
        json.dumps(
            {
                "schema": "frankensim.beads-hierarchy-normalization-plan.v1",
                "manifest_sha256": hashlib.sha256(
                    json.dumps(manifest, sort_keys=True).encode()
                ).hexdigest(),
                "migrate_rows": sum(1 for r in mrows if r["verdict"] == "migrate"),
                "kept_rows": sum(
                    1 for r in mrows if r["verdict"] == "keep-functional-prerequisite"
                ),
                "commands": plan,
                "inverse_commands": inverse,
                "note": "inverse plan is retained as explicit br commands and is never auto-applied",
            },
            indent=2,
            sort_keys=True,
        )
    )
    sys.exit(0)

print(f"unknown mode {mode}", file=sys.stderr)
sys.exit(EXIT_USAGE)
PYEOF
}

if [ "$MODE" = "self-test" ]; then
  SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/hiernorm-selftest.XXXXXX")"
  trap 'rm -rf "$SCRATCH"' EXIT
  PASS=0; FAIL=0
  check() {
    if [ "$2" -eq "$3" ]; then PASS=$((PASS + 1)); else
      FAIL=$((FAIL + 1)); echo "SELF-TEST FAIL: $1 (expected $2, got $3)" >&2
    fi
  }
  FIXDB="$SCRATCH/issues.jsonl"
  printf '%s\n' \
    '{"id":"x-epic","status":"open","description":"epic","dependencies":[{"depends_on_id":"x-epic.1","type":"blocks"},{"depends_on_id":"x-other","type":"blocks"}]}' \
    '{"id":"x-epic.1","status":"open","description":"child"}' \
    '{"id":"x-other","status":"open","description":"unrelated functional dep"}' \
    > "$FIXDB"
  LISTED="$(FSIM_HIERNORM_BEADS="$FIXDB" "$0" --list 2>/dev/null | wc -l | tr -d ' ')"
  [ "$LISTED" = "1" ]
  check "discovery finds exactly the dot-prefix blocks edge" 0 $?
  FSIM_HIERNORM_BEADS="$FIXDB" "$0" --freeze "$SCRATCH/m.json" >/dev/null 2>&1
  FSIM_HIERNORM_BEADS="$FIXDB" "$0" --check "$SCRATCH/m.json" >/dev/null 2>&1
  check "fresh manifest checks clean" 0 $?
  printf '%s\n' '{"id":"x-unrelated-new","status":"open"}' >> "$FIXDB"
  FSIM_HIERNORM_BEADS="$FIXDB" "$0" --check "$SCRATCH/m.json" >/dev/null 2>&1
  check "unrelated DB churn does NOT invalidate the freeze" 0 $?
  printf '%s\n' '{"id":"x-epic.2","status":"open","description":"new child"}' \
    '{"issue_id":"x-epic","depends_on_id":"x-epic.2","type":"blocks"}' >> "$FIXDB"
  FSIM_HIERNORM_BEADS="$FIXDB" "$0" --check "$SCRATCH/m.json" >/dev/null 2>&1
  check "a NEW candidate edge makes the manifest STALE (31)" 31 $?
  FSIM_HIERNORM_BEADS="$FIXDB" "$0" --plan "$SCRATCH/m.json" >/dev/null 2>&1
  check "a pending row refuses the plan (32)" 32 $?
  python3 - "$SCRATCH/m.json" <<'PYFIX'
import json, sys
manifest = json.load(open(sys.argv[1]))
manifest["rows"][0]["verdict"] = "migrate"
json.dump(manifest, open(sys.argv[1], "w"))
PYFIX
  FSIM_HIERNORM_BEADS="$FIXDB" "$0" --plan "$SCRATCH/m.json" > "$SCRATCH/plan.json" 2>/dev/null
  check "a fully reviewed manifest plans" 0 $?
  python3 - "$SCRATCH/plan.json" <<'PYPLAN'
import json, sys
plan = json.load(open(sys.argv[1]))
commands = plan["commands"]
assert commands[0][:3] == ["br", "dep", "remove"], commands
assert commands[1][:3] == ["br", "dep", "add"], commands
assert plan["inverse_commands"], "inverse plan retained"
PYPLAN
  check "plan removes FIRST then adds, inverse retained" 0 $?
  python3 - "$SCRATCH/m.json" <<'PYBAD'
import json, sys
manifest = json.load(open(sys.argv[1]))
manifest["rows"][0]["verdict"] = "looks-fine"
json.dump(manifest, open(sys.argv[1], "w"))
PYBAD
  FSIM_HIERNORM_BEADS="$FIXDB" "$0" --plan "$SCRATCH/m.json" >/dev/null 2>&1
  check "an out-of-vocabulary verdict refuses (30)" 30 $?
  "$0" --apply >/dev/null 2>&1
  check "--apply refuses by absence (30)" 30 $?
  echo "self-test: $PASS passed, $FAIL failed"
  [ "$FAIL" -eq 0 ] || exit 1
  exit 0
fi

case "$MODE" in
  list) tool list ;;
  freeze)
    TARGET="${MANIFEST_ARG:-$DEFAULT_MANIFEST}"
    mkdir -p "$(dirname "$TARGET")"
    tool freeze > "$TARGET" || exit $?
    COUNT="$(python3 -c "import json,sys;print(json.load(open(sys.argv[1]))['expected_relation_count'])" "$TARGET")"
    echo "frozen: $TARGET ($COUNT candidate rows, all pending review)" ;;
  check) tool check "$MANIFEST_ARG" ;;
  plan) tool plan "$MANIFEST_ARG" ;;
esac
