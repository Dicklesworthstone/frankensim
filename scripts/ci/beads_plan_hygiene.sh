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
# --apply, --negative, --replay, --owner-report) intentionally do not exist
# in this revision and refuse by absence; they land only after a reviewed
# deferred-exempt manifest and an owner-evidence workflow exist.
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

print(f"self-test: {PASS} passed, {FAIL} failed")
sys.exit(0 if FAIL == 0 else 1)
SELFTEST
    ;;
  --plan|--apply|--negative|--replay|--owner-report)
    die "$EXIT_UNIMPLEMENTED" \
      "$MODE lands only after a reviewed deferred-exempt manifest and owner-evidence workflow exist (bead frankensim-ug0yb slice gate); this revision is read-only" ;;
  *)
    die "$EXIT_USAGE" "usage: beads_plan_hygiene.sh --list [--out PATH] [MANIFEST] | --check MANIFEST | --self-test" ;;
esac

exec python3 - "$MODE" "$BEADS_FILE" "$DEFAULT_MANIFEST" "$OUT" "$MANIFEST" <<'PYEOF'
# ---CORE-BEGIN---
import hashlib
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone

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

    print(f"unknown mode {mode}", file=sys.stderr)
    return 30


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:6]))
# ---CORE-END---
PYEOF
