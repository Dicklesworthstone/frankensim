#!/usr/bin/env bash
# claim_integrity_inventory.sh (bead frankensim-extreal-program-f85xj.2.1) —
# the checked live inventory of claim-integrity defects.
#
# A claim-integrity defect is a public surface that can assert a STRONGER
# epistemic state than its evidence establishes (docs/CLAIM_INTEGRITY.md).
# Those defects only gate anything if they are countable, so this lane turns
# `br list -l claim-integrity` into a verified report:
#
#   * every entry is logged verbosely: status, priority, severity, crate scope,
#     assignee, and the taxonomy verdict for that entry;
#   * FAIL when an open P0 (severity:default-path) defect has no owner — an
#     unowned P0 is an inventory entry nobody is burning down;
#   * FAIL when an entry carries zero or several severity labels — the gate
#     reads severity to decide what to block, so an ambiguous severity is an
#     ungateable entry;
#   * FAIL CLOSED when the beads store cannot be read or parsed. An inventory
#     that cannot be read is NOT an empty inventory, and a lane that reported
#     "0 defects" because `br` was missing would itself be the exact defect
#     class this lane exists to count.
#
# Priority/severity disagreement and missing crate scope are reported as
# warnings: they degrade the gate's precision without making the entry
# invisible, and closed historical beads keep their original priority by
# policy (history is a record, not a thing to rewrite).
#
# Output is JSON-lines verdict rows plus a human summary, matching the repo's
# output convention. Deterministic: entries are sorted by id.
#
# Usage:
#   scripts/ci/claim_integrity_inventory.sh
#   FSIM_CLAIM_INTEGRITY_JSON=fixture.json scripts/ci/claim_integrity_inventory.sh
#
# The FSIM_CLAIM_INTEGRITY_JSON override feeds the lane a fixture snapshot
# instead of the live store; bead .2.4's seeded-fault drills use it to prove
# this lane fails closed without mutating the real beads database.
#
# Exit codes: 0 clean, 1 inventory defects found, 2 fail-closed read error.
set -euo pipefail

cd "$(dirname "$0")/../.."

LABEL="claim-integrity"
FIXTURE="${FSIM_CLAIM_INTEGRITY_JSON:-}"
RAW=""
RAW_CLOSED=""

if [[ -n "$FIXTURE" ]]; then
  if [[ ! -f "$FIXTURE" ]]; then
    printf '{"lane":"claim-integrity-inventory","verdict":"FAIL_CLOSED","reason":"fixture %s does not exist"}\n' "$FIXTURE"
    exit 2
  fi
  printf 'claim-integrity inventory: reading FIXTURE %s (not the live store)\n' "$FIXTURE" >&2
  RAW=$(cat "$FIXTURE")
else
  if ! command -v br >/dev/null 2>&1; then
    printf '{"lane":"claim-integrity-inventory","verdict":"FAIL_CLOSED","reason":"br not on PATH; an unreadable inventory is not an empty inventory"}\n'
    exit 2
  fi
  if ! RAW=$(CI=1 br list -l "$LABEL" --limit 4000 --json 2>/dev/null); then
    printf '{"lane":"claim-integrity-inventory","verdict":"FAIL_CLOSED","reason":"br list -l %s failed; refusing to report an empty inventory"}\n' "$LABEL"
    exit 2
  fi
  # `br list` omits closed beads, so the open query alone would let the report
  # print "0 closed" — an inventory that has not looked is not an inventory
  # that found nothing. The closed set is the sweep's known-answer set; it is
  # reported but never gates.
  if ! RAW_CLOSED=$(CI=1 br list -l "$LABEL" --status=closed --limit 4000 --json 2>/dev/null); then
    printf '{"lane":"claim-integrity-inventory","verdict":"FAIL_CLOSED","reason":"br list -l %s --status=closed failed; refusing to report an unknown closed history as empty"}\n' "$LABEL"
    exit 2
  fi
fi

# The payload travels by file, not by stdin: stdin belongs to the here-doc
# that carries the program, so a pipe here would be silently discarded and the
# lane would report an empty inventory — the exact fail-open shape this lane
# exists to catch.
PAYLOAD=$(mktemp "${TMPDIR:-/tmp}/claim-integrity-inventory.XXXXXX")
PAYLOAD_CLOSED=$(mktemp "${TMPDIR:-/tmp}/claim-integrity-closed.XXXXXX")
trap 'rm -f "$PAYLOAD" "$PAYLOAD_CLOSED"' EXIT
printf '%s' "$RAW" >"$PAYLOAD"
printf '%s' "$RAW_CLOSED" >"$PAYLOAD_CLOSED"

python3 - "$LABEL" "$PAYLOAD" "$PAYLOAD_CLOSED" <<'PY'
import json
import sys

label = sys.argv[1]
with open(sys.argv[2], "r", encoding="utf-8") as handle:
    raw = handle.read()
with open(sys.argv[3], "r", encoding="utf-8") as handle:
    raw_closed = handle.read()

def fail_closed(reason):
    print(json.dumps({
        "lane": "claim-integrity-inventory",
        "verdict": "FAIL_CLOSED",
        "reason": reason,
    }))
    raise SystemExit(2)

# Read-then-validate: another agent may be flushing the store mid-read, and a
# truncated document must refuse rather than parse as a short inventory.
try:
    document = json.loads(raw)
except json.JSONDecodeError as error:
    fail_closed(f"beads inventory is not valid JSON ({error}); the store may be mid-flush")

issues = document if isinstance(document, list) else document.get("issues")
if not isinstance(issues, list):
    fail_closed("beads inventory has no issues array")
if isinstance(document, dict) and document.get("has_more"):
    fail_closed("beads inventory was truncated by the query limit; a partial inventory undercounts")

# Merge the closed history (known-answer set). Absent for fixture runs, where
# the fixture is the whole world by construction.
if raw_closed.strip():
    try:
        closed_document = json.loads(raw_closed)
    except json.JSONDecodeError as error:
        fail_closed(f"closed claim-integrity history is not valid JSON ({error})")
    closed_issues = (
        closed_document if isinstance(closed_document, list)
        else closed_document.get("issues")
    )
    if not isinstance(closed_issues, list):
        fail_closed("closed claim-integrity history has no issues array")
    if isinstance(closed_document, dict) and closed_document.get("has_more"):
        fail_closed("closed claim-integrity history was truncated by the query limit")
    seen = {str(i.get("id")) for i in issues}
    issues = issues + [i for i in closed_issues if str(i.get("id")) not in seen]

SWEEP_SENTINEL = "Filed by the E02 claim-integrity audit sweep"

SEVERITY_PRIORITY = {
    "severity:default-path": 0,
    "severity:gated": 1,
    "severity:doc-only": 2,
}
CLOSED = {"closed"}

failures = []
warnings = []
rows = []

for issue in sorted(issues, key=lambda i: str(i.get("id", ""))):
    issue_id = str(issue.get("id", "<no-id>"))
    labels = issue.get("labels") or []
    if label not in labels:
        # br filtered by label; an entry without it means the query semantics
        # changed underneath us. Refuse rather than silently drop the entry.
        fail_closed(f"{issue_id} lacks the {label!r} label despite a label-filtered query")

    status = str(issue.get("status", "unknown"))
    priority = issue.get("priority")
    issue_type = str(issue.get("issue_type") or issue.get("type") or "")
    assignee = issue.get("assignee") or issue.get("owner") or ""
    severities = sorted(l for l in labels if l.startswith("severity:"))
    scopes = sorted(l for l in labels if l.startswith("crate:"))
    is_open = status not in CLOSED

    # The label also marks the E02 program (this doctrine, its sweep, its
    # gate). Those are work, not exposure; the bead type is the discriminator.
    # Counting them as defects would make the gate block on its own epic
    # forever, and a gate that can never go green teaches everyone to bypass it.
    if issue_type != "bug":
        rows.append({
            "lane": "claim-integrity-inventory",
            "id": issue_id,
            "verdict": "PROGRAM",
            "status": status,
            "priority": priority,
            "issue_type": issue_type,
            "severity": "<n/a: program bead>",
            "scopes": scopes,
            "assignee": assignee or "<unassigned>",
            "title": str(issue.get("title", ""))[:120],
            "failures": [],
            "warnings": [],
        })
        continue

    entry_failures = []
    entry_warnings = []

    # Taxonomy: exactly one severity label, drawn from the canonical set.
    if len(severities) != 1:
        entry_failures.append(
            f"carries {len(severities)} severity labels {severities}; the gate reads severity to "
            "decide what to block, so exactly one canonical label is required"
        )
        severity = None
    else:
        severity = severities[0]
        if severity not in SEVERITY_PRIORITY:
            entry_failures.append(
                f"severity {severity!r} is not canonical; expected one of "
                f"{sorted(SEVERITY_PRIORITY)}"
            )
            severity = None

    # An open P0 with nobody on it is an entry nobody is burning down.
    gating = severity == "severity:default-path" or priority == 0
    if is_open and gating and not assignee:
        entry_failures.append(
            "is an open gating (P0 / severity:default-path) defect with no owner; "
            "assign it with `br update <id> --assignee <agent>`"
        )

    # A finding without a reachable repro and a stated honest claim is not
    # actionable: the fix bead needs a destination, and the gate's count must
    # mean "verified exposure", not "someone's suspicion". Beads filed by the
    # sweep are held to this; retro-tagged history predates the schema and is
    # only warned about, since rewriting closed history is not the goal.
    body = str(issue.get("description") or "")
    swept = SWEEP_SENTINEL in body
    missing_fields = [
        name for name, marker in (("repro", "REPRO"), ("honest claim", "HONEST CLAIM"))
        if marker not in body
    ]
    if missing_fields:
        message = (
            f"is missing required field(s) {missing_fields}; every defect must carry a "
            "reachable repro and the honest claim the surface should make instead "
            "(docs/CLAIM_INTEGRITY.md)"
        )
        (entry_failures if swept else entry_warnings).append(message)

    if is_open and severity is not None:
        expected = SEVERITY_PRIORITY[severity]
        if priority != expected:
            entry_warnings.append(
                f"priority P{priority} disagrees with {severity} (expected P{expected})"
            )
    if is_open and not scopes:
        entry_warnings.append(
            "has no crate: scope; the promotion gate treats unscoped defects as GLOBAL "
            "(fail closed), which blocks every promotion"
        )

    verdict = "FAIL" if entry_failures else ("WARN" if entry_warnings else "OK")
    rows.append({
        "lane": "claim-integrity-inventory",
        "id": issue_id,
        "verdict": verdict,
        "status": status,
        "priority": priority,
        "issue_type": issue_type,
        "severity": severity or "<ambiguous>",
        "scopes": scopes,
        "assignee": assignee or "<unassigned>",
        "title": str(issue.get("title", ""))[:120],
        "failures": entry_failures,
        "warnings": entry_warnings,
    })
    failures.extend(f"{issue_id}: {message}" for message in entry_failures)
    warnings.extend(f"{issue_id}: {message}" for message in entry_warnings)

for row in rows:
    print(json.dumps(row, sort_keys=True))

defects = [r for r in rows if r["verdict"] != "PROGRAM"]
program = [r for r in rows if r["verdict"] == "PROGRAM"]
open_rows = [r for r in defects if r["status"] not in CLOSED]
gating_open = [
    r for r in open_rows
    if r["severity"] == "severity:default-path" or r["priority"] == 0
]

summary = {
    "lane": "claim-integrity-inventory",
    "verdict": "FAIL" if failures else "PASS",
    "labelled": len(rows),
    "program_beads": len(program),
    "defects": len(defects),
    "open": len(open_rows),
    "closed": len(defects) - len(open_rows),
    "open_gating": len(gating_open),
    "open_gating_ids": [r["id"] for r in gating_open],
    "failures": len(failures),
    "warnings": len(warnings),
}
print(json.dumps(summary, sort_keys=True))

print("", file=sys.stderr)
print("claim-integrity inventory", file=sys.stderr)
print(f"  labelled beads  : {len(rows)} ({len(program)} E02 program, {len(defects)} defects)", file=sys.stderr)
print(f"  defects         : {len(open_rows)} open, {summary['closed']} closed", file=sys.stderr)
print(f"  open gating (P0): {len(gating_open)}", file=sys.stderr)
for row in gating_open:
    print(f"      {row['id']}  [{row['assignee']}]  {row['title']}", file=sys.stderr)
if warnings:
    print(f"  warnings        : {len(warnings)}", file=sys.stderr)
    for message in warnings:
        print(f"      WARN {message}", file=sys.stderr)
if failures:
    print(f"  failures        : {len(failures)}", file=sys.stderr)
    for message in failures:
        print(f"      FAIL {message}", file=sys.stderr)
    print("", file=sys.stderr)
    print("Inventory is not gate-ready; see docs/CLAIM_INTEGRITY.md.", file=sys.stderr)
    raise SystemExit(1)

print("  verdict         : PASS (every entry is countable and every open P0 is owned)", file=sys.stderr)
PY
