#!/usr/bin/env bash
#
# Formula-bound Beads planning-hygiene inventory.
#
# Bead: frankensim-ug0yb (slice 1: read-only inventory core).
#
# Planning and ownership inventories must be formula-bound rather than
# copied from one filtered snapshot. This revision computes the exact row
# sets defined by the bead contract from a frozen live-DB identity:
#   active_status_set         = {open,in_progress,blocked}
#   raw_missing               = active rows whose estimate is null or <= 0
#   deferred_exempt           = membership of a REVIEWED manifest only
#   forecast_required         = raw_missing - deferred_exempt  (exact ids)
#   selected_now_unforecasted = br ready executable set ∩ forecast_required
#   stale_claims              = emitted rows of `br stale --days 7 --status
#                               in_progress` at the frozen basis
#   hard_priority_inversions  = direct type==blocks edges whose dependent's
#                               numeric priority is LOWER than its blocker's;
#                               related/parent-child/inactive excluded
# Counts are recomputed from emitted membership arrays and never replace
# them. Role classification is rule-based PROPOSAL with the matched rule
# recorded per row; it confers no authority. Mutation modes (--plan,
# --apply, --negative, --replay) intentionally do not exist yet;
# --owner-report exists as READ-ONLY evidence assembly: it emits typed
# staleness-evidence classes (ROSTER_ACTIVE, RECENT_PROOF, RECENTLY_UPDATED,
# IDLE_NO_EVIDENCE) per stale in_progress claim. Evidence classes are NOT
# dispositions; ABANDONED_CONFIRMED still requires owner response.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}"

BEADS_FILE="${FSIM_PLANHYG_BEADS:-$REPO_ROOT/.beads/issues.jsonl}"
DEFAULT_MANIFEST="$REPO_ROOT/tests/ci/beads_plan_hygiene/deferred-exempt-manifest.json"

EXIT_USAGE=30
EXIT_DRIFT=31
EXIT_UNIMPLEMENTED=32

die() {
  local class="$1"; shift
  printf 'beads-plan-hygiene: ERROR class=%s: %s\n' "$class" "$*" >&2
  exit "$class"
}

command -v python3 >/dev/null 2>&1 || die "$EXIT_USAGE" "python3 is required"
command -v br >/dev/null 2>&1 || die "$EXIT_USAGE" "br is required for ready/stale inputs"

MODE="${1:-}"
STALE_DAYS=7
case "$MODE" in
  --list)
    shift
    OUT=""
    if [ "${1:-}" = "--out" ]; then
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--out needs a path"
      OUT="$2"; shift 2
    fi
    [ $# -le 1 ] || die "$EXIT_USAGE" "usage: --list [--out PATH] [MANIFEST]"
    MANIFEST="${1:-$DEFAULT_MANIFEST}" ;;
  --check)
    [ $# -ge 2 ] || die "$EXIT_USAGE" "--check needs a deferred-exempt manifest path"
    MANIFEST="$2"; OUT="" ;;
  --self-test)
    # Extract the embedded python core to a temp module so fixture cases can
    # exercise the formulas directly without touching live tracker state.
    PLANHYG_PY_LIB="$(mktemp "${TMPDIR:-/tmp}/planhyg-core.XXXXXX")"
    trap 'rm -f "$PLANHYG_PY_LIB"' EXIT
    sed -n "/^# ---CORE-BEGIN---$/,/^# ---CORE-END---$/p" "$0" \
      | sed '1d;$d' > "$PLANHYG_PY_LIB"
    export PLANHYG_PY_LIB
    exec python3 - <<'SELFTEST'
from importlib.machinery import SourceFileLoader
import importlib.util
import os
import sys

PASS = 0
FAIL = 0

def check(name, expected, actual):
    global PASS, FAIL
    if expected == actual:
        PASS += 1
    else:
        FAIL += 1
        print(f"SELF-TEST FAIL: {name} (expected {expected!r}, got {actual!r})", file=sys.stderr)

_loader = SourceFileLoader("planhyg", os.environ["PLANHYG_PY_LIB"])
_spec = importlib.util.spec_from_loader("planhyg", _loader)
planhyg = importlib.util.module_from_spec(_spec)
_loader.exec_module(planhyg)

def issue(iid, status="open", est=None, prio=1, typ="task", deps=None):
    row = {"id": iid, "status": status, "priority": prio, "issue_type": typ}
    if est is not None:
        row["estimated_minutes"] = est
    if deps:
        row["dependencies"] = deps
    return row

# 1. raw_missing exact membership over status x estimate combinations.
rows = [
    issue("a.open.missing"),
    issue("b.open.zero", est=0),
    issue("c.open.neg", est=-5),
    issue("d.open.pos", est=30),
    issue("e.ip.missing", status="in_progress"),
    issue("f.blocked.missing", status="blocked"),
    issue("g.closed.missing", status="closed"),
]
issues = {r["id"]: r for r in rows}
raw = planhyg.raw_missing(issues)
check("raw_missing exact set",
      {"a.open.missing", "b.open.zero", "c.open.neg", "e.ip.missing", "f.blocked.missing"}, raw)

# 2. forecast_required is exact set subtraction (empty exempt boundary).
check("forecast subtraction empty-exempt",
      sorted(raw), sorted(planhyg.forecast_required(raw, set())))
deferred = {"b.open.zero"}
check("forecast subtraction one-exempt",
      sorted({"a.open.missing", "c.open.neg", "e.ip.missing", "f.blocked.missing"}),
      sorted(planhyg.forecast_required(raw, deferred)))

# 3. Inversions: direct blocks only; related/parent-child excluded; ties and
#    inactive endpoints excluded; duplicates collapsed.
edges = [
    {"issue_id": "p1.task", "depends_on_id": "p3.prereq", "type": "blocks"},
    {"issue_id": "p1.task", "depends_on_id": "p3.prereq", "type": "blocks"},
    {"issue_id": "p2.a", "depends_on_id": "p2.b", "type": "blocks"},
    {"issue_id": "p1.rel", "depends_on_id": "p3.other", "type": "related"},
    {"issue_id": "pc.child", "depends_on_id": "pc.parent", "type": "parent-child"},
    {"issue_id": "closed.dep", "depends_on_id": "p3.x", "type": "blocks"},
]
graph_issues = {
    "p1.task": issue("p1.task", prio=1),
    "p3.prereq": issue("p3.prereq", prio=3),
    "p2.a": issue("p2.a", prio=2),
    "p2.b": issue("p2.b", prio=2),
    "p1.rel": issue("p1.rel", prio=1),
    "p3.other": issue("p3.other", prio=3),
    "pc.child": issue("pc.child", prio=1),
    "pc.parent": issue("pc.parent", prio=3),
    "closed.dep": issue("closed.dep", status="closed", prio=1),
    "p3.x": issue("p3.x", prio=3),
}
inv = planhyg.hard_inversions(graph_issues, edges)
check("inversions detect exactly the direct hard blocks edge",
      [("p1.task", "p3.prereq")], inv)

# 4. Ready intersection join.
check("selected_now_unforecasted join",
      {"c.open.neg"},
      planhyg.intersect({"c.open.neg", "d.open.pos", "zz.ready"}, raw))

# 5. Input-order invariance.
shuffled = dict(reversed(list(issues.items())))
check("order invariance", sorted(raw), sorted(planhyg.raw_missing(shuffled)))

# 6. Zero-item and one-item boundaries.
check("zero missing boundary", set(), planhyg.raw_missing({"only": issue("only", est=10)}))
check("one missing boundary", {"solo"}, planhyg.raw_missing({"solo": issue("solo")}))

from datetime import datetime, timezone

NOW = datetime(2026, 8, 23, 12, 0, tzinfo=timezone.utc)

def claim(updated="2026-08-20T00:00:00Z", comments=None, assignee=None):
    row = {"updated_at": updated, "comments": comments or []}
    if assignee is not None:
        row["assignee"] = assignee
    return row

klass, _ = planhyg.classify_claim(
    claim(assignee="NobleLion"), NOW, 7,
    {"NobleLion": "2026-08-23T11:30:00Z"})
check("roster-active assignee classifies ROSTER_ACTIVE", "ROSTER_ACTIVE", klass)


klass, _ = planhyg.classify_claim(
    claim(comments=[{"author": "x", "created_at": "2026-08-20T00:00:00Z",
                     "text": "battery green, pushed"}]), NOW, 7, {})
check("recent proof comment classifies RECENT_PROOF", "RECENT_PROOF", klass)

klass, _ = planhyg.classify_claim(claim(), NOW, 7, {})
check("freshly touched no-comment classifies RECENTLY_UPDATED",
      "RECENTLY_UPDATED", klass)

klass, detail = planhyg.classify_claim(
    claim(updated="2026-07-14T00:00:00Z",
          comments=[{"author": "x", "created_at": "2026-07-15T00:00:00Z",
                     "text": "looking at it"}]), NOW, 7, {})
check("old idle claim classifies IDLE_NO_EVIDENCE", "IDLE_NO_EVIDENCE", klass)

klass, _ = planhyg.classify_claim(
    claim(comments=[{"author": "x", "created_at": "2026-08-20T00:00:00Z",
                     "text": "battery green, pushed"}],
          assignee="NobleLion"), NOW, 7,
    {"NobleLion": "2026-08-23T11:30:00Z"})
check("roster activity outranks proof comment", "ROSTER_ACTIVE", klass)
klass, _ = planhyg.classify_claim(
    claim(), NOW, 7, {}, last_commit_days=2)
check("recent bead-mentioning commit classifies COMMIT_ACTIVE",
      "COMMIT_ACTIVE", klass)

index = planhyg.commit_activity_index(
    ["2026-08-22T10:00:00+00:00\tfeat: closes frankensim-b (frankensim-a)",
     "2026-08-01T10:00:00+00:00\told frankensim-a work",
     "no tab line"],
    {"frankensim-a", "frankensim-b"}, NOW, 7)
check("commit index keeps freshest in-cutoff mention per id",
      {"frankensim-a": 1.08, "frankensim-b": 1.08},
      {k: round(v, 2) for k, v in index.items()})


print(f"self-test: {PASS} passed, {FAIL} failed")
sys.exit(0 if FAIL == 0 else 1)
SELFTEST
    ;;
  --owner-report)
    shift
    STALE_DAYS=7
    COMMIT_EVIDENCE=""
    if [ "${1:-}" = "--days" ]; then
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--days needs a number"
      STALE_DAYS_OVERRIDE="$2"; shift 2
    fi
    if [ "${1:-}" = "--commits" ]; then
      COMMIT_EVIDENCE="1"; shift 1
    fi
    OUT=""
    if [ "${1:-}" = "--out" ]; then
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--out needs a path"
      OUT="$2"; shift 2
    fi
    [ $# -le 0 ] || die "$EXIT_USAGE" "usage: --owner-report [--days N] [--commits] [--out PATH]"
    [ -z "${STALE_DAYS_OVERRIDE:-}" ] || STALE_DAYS="$STALE_DAYS_OVERRIDE"
    MANIFEST="" ;;
  --plan|--apply|--negative|--replay)
    die "$EXIT_UNIMPLEMENTED" \
      "$MODE lands only after a reviewed deferred-exempt manifest and owner-evidence workflow exist (bead frankensim-ug0yb slice gate); this revision is read-only" ;;
  *)
    die "$EXIT_USAGE" "usage: beads_plan_hygiene.sh --list [--out PATH] [MANIFEST] | --check MANIFEST | --owner-report [--days N] [--out PATH] | --self-test" ;;
esac

exec python3 - "$MODE" "$BEADS_FILE" "$DEFAULT_MANIFEST" "$OUT" "$MANIFEST" "$STALE_DAYS" "${COMMIT_EVIDENCE:-}" <<'PYEOF'
# ---CORE-BEGIN---
import hashlib
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone, timedelta

ACTIVE = {"open", "in_progress", "blocked"}
SCHEMA = "frankensim.beads-plan-hygiene-inventory.v1"

HORIZON_ID_PREFIXES = ("frankensim-epic-addendum-xpck",)
ROLE_RULES = (
    ("horizon_deferral", lambda row: any(row["id"].startswith(p) for p in HORIZON_ID_PREFIXES)
     or bool({"horizon", "deferral", "deferred"} & set(row.get("labels") or []))),
    ("human_study", lambda row: bool(re.search(
        r"study|cohort|participant|preregistration|interview", row.get("title", ""), re.I))),
    ("performance_campaign", lambda row: bool(re.search(
        r"\bbenchmark|\bcampaign\b|throughput|device matrix", row.get("title", ""), re.I))),
    ("evidence_adjudication_leaf", lambda row: bool(re.search(
        r"adjudicat|receipt|independent(ly)? (verify|audit)|closeout evidence",
        row.get("title", ""), re.I))),
)


def load_beads(path):
    issues, edges = {}, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not isinstance(rec, dict):
                continue
            if "issue_id" in rec and "depends_on_id" in rec:
                edges.append(rec)
            elif "id" in rec:
                issues[rec["id"]] = rec
                for dep in rec.get("dependencies") or []:
                    if isinstance(dep, dict) and "depends_on_id" in dep:
                        merged = dict(dep)
                        merged.setdefault("issue_id", rec["id"])
                        edges.append(merged)
    return issues, edges


def raw_missing(issues):
    """Active rows whose estimate is null or nonpositive."""
    out = set()
    for iid, rec in issues.items():
        if rec.get("status") not in ACTIVE:
            continue
        est = rec.get("estimated_minutes")
        if isinstance(est, bool) or not isinstance(est, (int, float)) or est <= 0:
            out.add(iid)
    return out


def forecast_required(missing, exempt):
    return set(missing) - set(exempt)


def intersect(left, right):
    return set(left) & set(right)


def hard_inversions(issues, edges):
    """Direct type==blocks edges whose DEPENDENT has numerically lower
    priority than its BLOCKER. Both endpoints must be active; related and
    parent-child edges are excluded by contract."""
    seen, out = set(), []
    for edge in edges:
        if edge.get("type") != "blocks":
            continue
        dependent, blocker = edge.get("issue_id"), edge.get("depends_on_id")
        key = (dependent, blocker)
        if key in seen or dependent not in issues or blocker not in issues:
            continue
        if issues[dependent].get("status") not in ACTIVE:
            continue
        if issues[blocker].get("status") not in ACTIVE:
            continue
        pd = issues[dependent].get("priority")
        pb = issues[blocker].get("priority")
        if isinstance(pd, int) and isinstance(pb, int) and pd < pb:
            seen.add(key)
            out.append(key)
    return sorted(out)


def child_counts(issues):
    counts = {}
    for iid in issues:
        parent, sep, _leaf = iid.rpartition(".")
        while sep:
            counts[parent] = counts.get(parent, 0) + 1
            parent, sep, _leaf = parent.rpartition(".")
    return counts


def classify(iid, rec, child_counts_map):
    title = rec.get("title", "")
    labels = rec.get("labels") or []
    probe = {"id": iid, "title": title, "labels": labels}
    if rec.get("issue_type") == "epic":
        n_children = child_counts_map.get(iid, 0)
        if n_children:
            return "exact_set_rollup", f"epic with {n_children} dot-prefix children"
        for name, pred in ROLE_RULES:
            if pred(probe):
                return name, f"epic without children; matched rule {name}"
        return "malformed_container", "epic without children; no deferral rule matched"
    for name, pred in ROLE_RULES:
        if pred(probe):
            return name, f"matched rule {name}"
    return "executable_leaf", "default classification"


def load_manifest(path):
    try:
        with open(path) as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        print(f"manifest refusal: {error}", file=sys.stderr)
        sys.exit(30)
    if data.get("schema") != "frankensim.beads-plan-deferred-exempt.v1":
        print(f"manifest refusal: schema {data.get('schema')!r}", file=sys.stderr)
        sys.exit(30)
    rows = data.get("rows")
    if not isinstance(rows, list):
        print("manifest refusal: rows missing", file=sys.stderr)
        sys.exit(30)
    return data, rows


PROOF_PATTERN = re.compile(
    r"test|receipt|proof|pushed|landed|green", re.I)

BEAD_ID_PATTERN = re.compile(r"frankensim-[A-Za-z0-9]+(?:[.\-][A-Za-z0-9]+)*")


def iso_age_days(iso_text, now):
    if not isinstance(iso_text, str) or not iso_text:
        return None
    try:
        then = datetime.fromisoformat(iso_text.replace("Z", "+00:00"))
    except ValueError:
        return None
    if then.tzinfo is None:
        then = then.replace(tzinfo=timezone.utc)
    return (now - then).total_seconds() / 86400.0


def commit_activity_index(commit_lines, ids, now, cutoff_days):
    """Map bead id -> freshest mentioning-commit age (days), for ids in the
    stale set only. Lines are 'ISO-DATE<TAB>SUBJECT' from git log."""
    best = {}
    for line in commit_lines:
        if "\t" not in line:
            continue
        date_text, subject = line.split("\t", 1)
        age = iso_age_days(date_text.strip(), now)
        if age is None or age > cutoff_days:
            continue
        for iid in BEAD_ID_PATTERN.findall(subject):
            if iid in ids and (iid not in best or age < best[iid]):
                best[iid] = age
    return best


def classify_claim(rec, now, stale_days, roster, last_commit_days=None):
    """Pure staleness-evidence classifier for one in_progress claim.
    Returns (evidence_class, detail). Classes in priority order:
    ROSTER_ACTIVE > COMMIT_ACTIVE > RECENT_PROOF > RECENTLY_UPDATED >
    IDLE_NO_EVIDENCE. These are EVIDENCE classes, never dispositions."""
    assignee = rec.get("assignee")
    roster_age = iso_age_days(roster.get(assignee), now) if assignee else None
    if roster_age is not None and roster_age <= 2 / 24:
        return "ROSTER_ACTIVE", {"basis": "assignee active on mail roster"}
    if last_commit_days is not None and last_commit_days <= stale_days:
        return "COMMIT_ACTIVE", {
            "basis": f"a commit referencing this bead landed "
                     f"{last_commit_days:.1f}d ago",
        }
    comments = rec.get("comments") or []
    last = max(
        comments,
        key=lambda c: c.get("created_at", "") if isinstance(c, dict) else "",
        default=None,
    )
    last_age = iso_age_days(
        last.get("created_at"), now) if isinstance(last, dict) else None
    last_text = (last.get("text") or "") if isinstance(last, dict) else ""
    proof = bool(last_text) and bool(PROOF_PATTERN.search(last_text))
    if last_age is not None and last_age <= stale_days and proof:
        return "RECENT_PROOF", {
            "basis": f"last comment {last_age:.1f}d ago cites work",
            "comment_author": last.get("author"),
            "comment_excerpt": last_text[:160],
        }
    updated_age = iso_age_days(rec.get("updated_at"), now)
    if updated_age is not None and updated_age <= stale_days:
        return "RECENTLY_UPDATED", {"basis": f"updated {updated_age:.1f}d ago"}
    return "IDLE_NO_EVIDENCE", {
        "basis": (
            f"updated {updated_age:.1f}d ago, no proof comment"
            if updated_age is not None
            else "no parsable timestamps"
        )
    }


def run_br(args):
    result = subprocess.run(["br", *args], capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        print(f"br {' '.join(args)} failed: {result.stderr[-200:]}", file=sys.stderr)
        sys.exit(30)
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError:
        print(f"br {' '.join(args)} returned unparsable output", file=sys.stderr)
        sys.exit(30)


def main(argv):
    mode = (argv[0] or "").lstrip("-")
    beads_path, default_manifest, out_arg, manifest_arg = argv[1:5]
    commit_evidence = bool(argv[6]) if len(argv) > 6 else False
    stale_days = int(argv[5]) if len(argv) > 5 and argv[5] else 7

    issues, edges = load_beads(beads_path)
    db_sha = hashlib.sha256(open(beads_path, "rb").read()).hexdigest()
    missing = raw_missing(issues)

    if mode == "list":
        if manifest_arg != default_manifest and not os.path.exists(manifest_arg):
            print(f"manifest refusal: not found: {manifest_arg}", file=sys.stderr)
            sys.exit(30)
        manifest_rows = []
        if os.path.exists(manifest_arg):
            _, manifest_rows = load_manifest(manifest_arg)
        exempt_ids = {row.get("id") for row in manifest_rows}
        unknown = sorted(e for e in exempt_ids if e not in missing)
        if unknown:
            print(
                "manifest drift refusal: reviewed exempt ids no longer in "
                f"raw_missing: {unknown[:10]}",
                file=sys.stderr,
            )
            sys.exit(31)
        required = forecast_required(missing, exempt_ids)
        ready_ids = {row.get("id") for row in run_br(["ready", "--json"])}
        stale_ids = {
            row.get("id")
            for row in run_br(["stale", "--days", "7", "--status", "in_progress", "--json"])
        }
        inversions = hard_inversions(issues, edges)
        child_counts_map = child_counts(issues)

        rows = [{
            "event": "inventory-header",
            "schema": SCHEMA,
            "bead": "frankensim-ug0yb",
            "formula_version": "ug0yb-slice1",
            "db_sha256": db_sha,
            "frozen_basis_utc": datetime.now(timezone.utc).isoformat(),
            "stale_cutoff_days": 7,
            "active_status_set": sorted(ACTIVE),
        }]
        for iid in sorted(required):
            rec = issues[iid]
            role, basis = classify(iid, rec, child_counts_map)
            rows.append({
                "event": "forecast_required_row",
                "id": iid,
                "status": rec.get("status"),
                "priority": rec.get("priority"),
                "issue_type": rec.get("issue_type"),
                "estimated_minutes": rec.get("estimated_minutes"),
                "role_proposed": role,
                "classification_basis": basis,
            })
        for row in sorted(manifest_rows, key=lambda r: r.get("id", "")):
            rows.append({"event": "deferred_exempt_row", **row})
        for sid in sorted(stale_ids):
            rows.append({"event": "stale_claim_row", "id": sid})
        for dependent, blocker in inversions:
            rows.append({
                "event": "hard_priority_inversion_row",
                "dependent": dependent,
                "blocker": blocker,
                "dependent_priority": issues[dependent].get("priority"),
                "blocker_priority": issues[blocker].get("priority"),
            })
        declared = {
            "raw_missing_count": len(missing),
            "deferred_exempt_count": len(exempt_ids),
            "forecast_required_count": len(required),
            "selected_now_unforecasted_count": len(intersect(ready_ids, required)),
            "stale_claims_count": len(stale_ids),
            "hard_priority_inversion_count": len(inversions),
        }
        derived = {
            "deferred_exempt_count": sum(1 for r in rows if r["event"] == "deferred_exempt_row"),
            "forecast_required_count": sum(1 for r in rows if r["event"] == "forecast_required_row"),
            "stale_claims_count": sum(1 for r in rows if r["event"] == "stale_claim_row"),
            "hard_priority_inversion_count": sum(
                1 for r in rows if r["event"] == "hard_priority_inversion_row"),
        }
        for field in ("forecast_required_count", "deferred_exempt_count",
                      "stale_claims_count", "hard_priority_inversion_count"):
            if declared[field] != derived[field]:
                print(
                    f"count/member drift on {field}: declared {declared[field]} "
                    f"vs emitted {derived[field]}",
                    file=sys.stderr,
                )
                sys.exit(33)
        rows.append({"event": "inventory-summary", **declared})
        rows.append({
            "event": "no-claim",
            "text": "forecasts are not deadlines; classifications are review proposals; "
                    "this inventory grants no implementation, scientific, performance, "
                    "release, or promotion authority",
        })
        rendered = "".join(json.dumps(r, sort_keys=True) + "\n" for r in rows)
        if out_arg:
            target = os.path.realpath(out_arg)
            if target == os.path.realpath(beads_path) or target.startswith(
                    os.path.realpath(beads_path) + os.sep):
                print(f"path-overlap refusal: {out_arg} overlaps tracker export",
                      file=sys.stderr)
                sys.exit(30)
            parent_dir = os.path.dirname(out_arg)
            if parent_dir:
                os.makedirs(parent_dir, exist_ok=True)
            with open(out_arg, "w") as handle:
                handle.write(rendered)
            print(f"inventory written: {out_arg}")
        else:
            sys.stdout.write(rendered)
        return 0

    if mode == "check":
        _, manifest_rows = load_manifest(manifest_arg)
        problems = []
        for row in manifest_rows:
            iid = row.get("id")
            rec = issues.get(iid)
            if rec is None:
                problems.append(("missing_record", iid))
            elif rec.get("status") not in ACTIVE:
                problems.append(("not_active", iid))
            elif iid not in missing:
                problems.append(("no_longer_missing_estimate", iid))
            for field in ("role", "reason"):
                if not row.get(field):
                    problems.append((f"absent_{field}", iid))
        if problems:
            print(
                f"manifest drift refusal: first divergence {problems[0]}; "
                f"{len(problems)} problem(s)",
                file=sys.stderr,
            )
            sys.exit(31)
        print(
            f"manifest OK: {len(manifest_rows)} reviewed exempt rows bound to live DB "
            f"identity {db_sha[:16]}"
        )
        return 0

    if mode == "owner-report":
        stale_rows = run_br([
            "stale", "--days", str(stale_days),
            "--status", "in_progress", "--json"])
        roster = {}
        roster_path = os.environ.get("FSIM_PLANHYG_MAIL_ROSTER", "")
        if roster_path and os.path.exists(roster_path):
            try:
                with open(roster_path) as handle:
                    roster_list = json.load(handle)
                roster = {
                    entry.get("name"): entry.get("last_active_ts")
                    for entry in roster_list
                    if isinstance(entry, dict)
                }
            except (OSError, json.JSONDecodeError):
                roster = {}
        now = datetime.now(timezone.utc)
        commit_index = {}
        if commit_evidence:
            repo_root = os.path.dirname(os.path.dirname(
                os.path.abspath(beads_path)))
            since = (now - timedelta(days=stale_days)).strftime("%Y-%m-%d")
            log_result = subprocess.run(
                ["git", "-C", repo_root, "log",
                 f"--since={since}", "--date=iso-strict",
                 "--format=%ad%x09%s"],
                capture_output=True, text=True, timeout=300)
            if log_result.returncode != 0:
                print(f"git log failed: {log_result.stderr[-200:]}",
                      file=sys.stderr)
                sys.exit(30)
            stale_ids = {row.get("id") for row in stale_rows}
            commit_index = commit_activity_index(
                log_result.stdout.splitlines(), stale_ids, now, stale_days)
        rows = [{
            "event": "owner-report-header",
            "schema": SCHEMA,
            "bead": "frankensim-ug0yb",
            "formula_version": "ug0yb-owner-report-slice3",
            "commits_bound": bool(commit_evidence),
            "db_sha256": db_sha,
            "frozen_basis_utc": now.isoformat(),
            "stale_cutoff_days": stale_days,
            "roster_bound": bool(roster),
            "roster_liveness_note": (
                "mail roster last_active_ts values are registration-time "
                "snapshots on this server build; ROSTER_ACTIVE is therefore "
                "conservative and may under-report live assignees"),
        }]
        counts = {}
        for claim in sorted(stale_rows, key=lambda r: r.get("id", "")):
            iid = claim.get("id")
            rec = issues.get(iid, {})
            klass, detail = classify_claim(
                rec, now, stale_days, roster,
                last_commit_days=commit_index.get(iid))
            counts[klass] = counts.get(klass, 0) + 1
            rows.append({
                "event": "owner_report_row",
                "id": iid,
                "title": (rec.get("title") or claim.get("title") or "")[:120],
                "priority": rec.get("priority"),
                "assignee": rec.get("assignee"),
                "evidence_class": klass,
                **detail,
            })
        derived_total = sum(counts.values())
        if derived_total != len(stale_rows):
            print(
                f"count/member drift: classified {derived_total} vs "
                f"stale rows {len(stale_rows)}",
                file=sys.stderr,
            )
            sys.exit(33)
        rows.append({"event": "owner-report-summary",
                     "stale_claims": len(stale_rows), **counts})
        rows.append({
            "event": "no-claim",
            "text": "evidence classes are not dispositions; "
                    "ABANDONED_CONFIRMED requires owner response per the "
                    "ug0yb contract; this report mutates nothing",
        })
        rendered = "".join(json.dumps(r, sort_keys=True) + "\n" for r in rows)
        if out_arg:
            target = os.path.realpath(out_arg)
            if target == os.path.realpath(beads_path) or target.startswith(
                    os.path.realpath(beads_path) + os.sep):
                print(f"path-overlap refusal: {out_arg} overlaps tracker export",
                      file=sys.stderr)
                sys.exit(30)
            parent_dir = os.path.dirname(out_arg)
            if parent_dir:
                os.makedirs(parent_dir, exist_ok=True)
            with open(out_arg, "w") as handle:
                handle.write(rendered)
            print(f"owner report written: {out_arg}")
        else:
            sys.stdout.write(rendered)
        return 0

    print(f"unknown mode {mode}", file=sys.stderr)
    return 30


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
# ---CORE-END---
PYEOF
