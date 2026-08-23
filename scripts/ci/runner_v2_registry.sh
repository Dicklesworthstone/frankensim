#!/usr/bin/env bash
#
# Runner V2 acceptance-corpus registry inventory.
#
# Bead: frankensim-epic-foundations-huq.24.6 (slice 1: read-only core).
#
# Builds a source-authoritative machine-readable inventory of the Runner V2
# acceptance corpus from a frozen live-DB identity and detects the
# mechanically decidable conflict classes before any prose-level clause
# normalization is attempted:
#   corpus            = huq-family issues with id under .24 plus any
#                       huq issue whose title contains "Runner V2"
#   contract_record   = one row per corpus producer (status, priority,
#                       type, owner, estimate)
#   conflict classes  = HIERARCHY_DUPLICATION (parent-child AND blocks
#                       edges on the same pair), CLOSED_CONTAINER_OPEN_CHILD,
#                       OPEN_CONTAINER_ALL_CHILDREN_CLOSED (informational),
#                       MISSING_PARENT (dot-prefix parent not in DB),
#                       OWNER_MISSING, ESTIMATE_MISSING
# Counts are recomputed from emitted rows and never substitute for
# membership. Mutation modes (--freeze-review, --apply) intentionally do
# not exist yet and refuse by absence until a reviewed normalization
# manifest exists.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}"

BEADS_FILE="${FSIM_RV2REG_BEADS:-$REPO_ROOT/.beads/issues.jsonl}"

EXIT_USAGE=30
EXIT_DRIFT=31
EXIT_UNIMPLEMENTED=32

die() {
  local class="$1"; shift
  printf 'runner-v2-registry: ERROR class=%s: %s\n' "$class" "$*" >&2
  exit "$class"
}

command -v python3 >/dev/null 2>&1 || die "$EXIT_USAGE" "python3 is required"

MODE="${1:-}"
case "$MODE" in
  --list)
    shift
    OUT=""
    if [ "${1:-}" = "--out" ]; then
      [ $# -ge 2 ] || die "$EXIT_USAGE" "--out needs a path"
      OUT="$2"; shift 2
    fi
    [ $# -le 0 ] || die "$EXIT_USAGE" "usage: --list [--out PATH]" ;;
  --check)
    [ $# -ge 2 ] || die "$EXIT_USAGE" "--check needs a registry JSONL path"
    OUT="$2" ;;
  --self-test)
    RV2REG_PY_LIB="$(mktemp "${TMPDIR:-/tmp}/rv2reg-core.XXXXXX")"
    trap 'rm -f "$RV2REG_PY_LIB"' EXIT
    sed -n "/^# ---CORE-BEGIN---$/,/^# ---CORE-END---$/p" "$0" \
      | sed '1d;$d' > "$RV2REG_PY_LIB"
    export RV2REG_PY_LIB
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

_loader = SourceFileLoader("rv2reg", os.environ["RV2REG_PY_LIB"])
_spec = importlib.util.spec_from_loader("rv2reg", _loader)
rv2reg = importlib.util.module_from_spec(_spec)
_loader.exec_module(rv2reg)

FAM = "frankensim-epic-foundations-huq"

def issue(iid, title="", status="open", prio=1, typ="task", owner=None, est=None):
    row = {"id": iid, "title": title, "status": status,
           "priority": prio, "issue_type": typ}
    if owner is not None:
        row["assignee"] = owner
    if est is not None:
        row["estimated_minutes"] = est
    return row

# Fixtures are keyed by FULL id exactly as production code sees them.
issues = {
    FAM: issue(FAM, "EPIC FOUNDATIONS"),
    FAM + ".24": issue(FAM + ".24", "RUNNER V2 PROGRAM"),
    FAM + ".24.1": issue(FAM + ".24.1", "Runner V2 framing"),
    FAM + ".31": issue(FAM + ".31", "unrelated foundations work"),
    FAM + ".77": issue(FAM + ".77", "Runner V2 side note"),
}

# 1. Discovery: .24 subtree by identity, title match anywhere in family,
#    unrelated family rows excluded.
corpus = rv2reg.corpus_ids(issues)
check("corpus discovery exact",
      {FAM + ".24", FAM + ".24.1", FAM + ".77"}, corpus)

# 2. Hierarchy duplication: pair carrying BOTH edge types is flagged;
#    single-type pairs (even duplicated within a type) are not.
both = [
    {"issue_id": FAM + ".24.1", "depends_on_id": FAM + ".24", "type": "parent-child"},
    {"issue_id": FAM + ".24.1", "depends_on_id": FAM + ".24", "type": "blocks"},
]
check("mixed-type duplicate pair detected",
      [(FAM + ".24.1", FAM + ".24")], rv2reg.hierarchy_duplicates(both))
single_type = [
    {"issue_id": FAM + ".24.1", "depends_on_id": FAM + ".24", "type": "parent-child"},
    {"issue_id": FAM + ".24.1", "depends_on_id": FAM + ".24", "type": "parent-child"},
]
check("same-type duplicate not flagged", [], rv2reg.hierarchy_duplicates(single_type))

# 3. Status contradiction: closed container with an active child; and the
#    informational stale-open inverse.
tree = {k: dict(v) for k, v in issues.items()}
tree[FAM + ".24"]["status"] = "closed"
containment = {(FAM + ".24.1", FAM + ".24"): True}
contras, stale_open = rv2reg.status_contradictions(tree, containment)
check("closed container open child flagged", [FAM + ".24"], contras)
check("stale-open empty here", [], stale_open)
tree[FAM + ".24"]["status"] = "open"
tree[FAM + ".24.1"]["status"] = "closed"
contras, stale_open = rv2reg.status_contradictions(tree, containment)
check("open container all children closed informational",
      [FAM + ".24"], stale_open)
check("no contradiction when container open", [], contras)

# 4. Missing parent detection: dot-prefix parent absent from DB.
orphans = rv2reg.missing_parents({FAM + ".24.1": tree[FAM + ".24.1"]})
check("missing parent listed once", [FAM + ".24"], orphans)

# 5. Owner / estimate flags on active producers only.
flags = rv2reg.owner_estimate_flags({
    FAM + ".a": issue(FAM + ".a", owner="NobleLion", est=30),
    FAM + ".b": issue(FAM + ".b"),
    FAM + ".c": issue(FAM + ".c", est=0),
    FAM + ".d": issue(FAM + ".d", status="closed"),
})
check("owner+estimate missing flagged", ["ESTIMATE_MISSING", "OWNER_MISSING"],
      sorted(n for n, _ in flags.get(FAM + ".b", [])))
check("zero estimate flagged alongside owner", ["ESTIMATE_MISSING", "OWNER_MISSING"],
      sorted(n for n, _ in flags.get(FAM + ".c", [])))
check("healthy producer unflagged", [], flags.get(FAM + ".a", []))
check("closed row unflagged", [], flags.get(FAM + ".d", []))

# 6. Input-order invariance.
shuffled = dict(reversed(list(issues.items())))
check("order invariance", sorted(corpus), sorted(rv2reg.corpus_ids(shuffled)))

print(f"self-test: {PASS} passed, {FAIL} failed")
sys.exit(0 if FAIL == 0 else 1)
SELFTEST
    ;;
  --freeze-review|--apply)
    die "$EXIT_UNIMPLEMENTED" \
      "$MODE lands only after a reviewed normalization manifest exists (bead frankensim-epic-foundations-huq.24.6 slice gate); this revision is read-only" ;;
  *)
    die "$EXIT_USAGE" "usage: runner_v2_registry.sh --list [--out PATH] | --check REGISTRY | --self-test" ;;
esac

exec python3 - "$MODE" "$BEADS_FILE" "$OUT" <<'PYEOF'
# ---CORE-BEGIN---
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

FAMILY_PREFIX = "frankensim-epic-foundations-huq"
PROGRAM_SUBTREE = FAMILY_PREFIX + ".24"
TITLE_MARKER = "Runner V2"
SCHEMA = "frankensim.runner-v2-registry-inventory.v1"
ACTIVE = {"open", "in_progress", "blocked"}


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


def corpus_ids(issues):
    out = set()
    for iid, rec in issues.items():
        if not iid.startswith(FAMILY_PREFIX):
            continue
        if iid.startswith(PROGRAM_SUBTREE) or TITLE_MARKER in rec.get("title", ""):
            out.add(iid)
    return out


def pair_types(edges):
    pairs = {}
    for edge in edges:
        key = (edge.get("issue_id"), edge.get("depends_on_id"))
        pairs.setdefault(key, set()).add(edge.get("type"))
    return pairs


def hierarchy_duplicates(edges):
    """Pairs carrying BOTH hierarchy semantics at once."""
    out = []
    for (dependent, blocker), types in sorted(pair_types(edges).items()):
        if "parent-child" in types and "blocks" in types:
            out.append((dependent, blocker))
    return out


def child_pairs(issues):
    """Dot-prefix containment restricted to pairs where both exist."""
    out = {}
    for iid in issues:
        parent, sep, _leaf = iid.rpartition(".")
        while sep:
            if parent in issues:
                out[(iid, parent)] = True
            parent, sep, _leaf = parent.rpartition(".")
    return out


def status_contradictions(issues, containment):
    closed_container_open_child, stale_open = [], []
    containers = {}
    for (child, parent) in containment:
        containers.setdefault(parent, []).append(child)
    for parent, children in sorted(containers.items()):
        prec = issues[parent]
        child_statuses = [issues[c].get("status") for c in children if c in issues]
        if not child_statuses:
            continue
        if prec.get("status") == "closed" and any(s in ACTIVE for s in child_statuses):
            closed_container_open_child.append(parent)
        elif prec.get("status") in ACTIVE and all(s == "closed" for s in child_statuses):
            stale_open.append(parent)
    return sorted(set(closed_container_open_child)), sorted(set(stale_open))


def missing_parents(issues):
    out = []
    for iid in issues:
        parent, sep, _leaf = iid.rpartition(".")
        if sep and parent not in issues:
            out.append(parent)
    return sorted(set(out))


def owner_estimate_flags(records):
    flags = {}
    for iid, rec in records.items():
        if rec.get("status") not in ACTIVE:
            continue
        owner = rec.get("assignee")
        est = rec.get("estimated_minutes")
        missing_est = isinstance(est, bool) or not isinstance(est, (int, float)) or est <= 0
        entry = []
        if not owner:
            entry.append(("OWNER_MISSING", True))
        if missing_est:
            entry.append(("ESTIMATE_MISSING", True))
        if entry:
            flags[iid] = entry
    return flags


def main(argv):
    mode = (argv[0] or "").lstrip("-")
    beads_path, out_arg = argv[1], argv[2]

    issues, edges = load_beads(beads_path)
    db_sha = hashlib.sha256(open(beads_path, "rb").read()).hexdigest()
    corpus = corpus_ids(issues)
    containment = child_pairs(issues)
    dup = hierarchy_duplicates(
        e for e in edges
        if e.get("issue_id") in corpus and e.get("depends_on_id") in corpus)
    contras, stale_open = status_contradictions(issues, containment)
    orphans = missing_parents(issues)
    records = {iid: issues[iid] for iid in sorted(corpus)}
    flags = owner_estimate_flags(records)

    if mode == "list":
        rows = [{
            "event": "registry-header",
            "schema": SCHEMA,
            "bead": "frankensim-epic-foundations-huq.24.6",
            "formula_version": "huq246-slice1",
            "db_sha256": db_sha,
            "frozen_basis_utc": datetime.now(timezone.utc).isoformat(),
            "discovery_rule": {
                "family_prefix": FAMILY_PREFIX,
                "program_subtree": PROGRAM_SUBTREE,
                "title_marker": TITLE_MARKER,
            },
        }]
        for iid in sorted(corpus):
            rec = issues[iid]
            rows.append({
                "event": "contract_record",
                "id": iid,
                "title": rec.get("title", ""),
                "status": rec.get("status"),
                "priority": rec.get("priority"),
                "issue_type": rec.get("issue_type"),
                "assignee": rec.get("assignee"),
                "estimated_minutes": rec.get("estimated_minutes"),
                "active_children": sum(
                    1 for (child, parent) in containment
                    if parent == iid and issues.get(child, {}).get("status") in ACTIVE),
            })
        for dependent, blocker in dup:
            rows.append({"event": "conflict_row", "class": "HIERARCHY_DUPLICATION",
                         "dependent": dependent, "blocker": blocker})
        for parent in contras:
            rows.append({"event": "conflict_row", "class": "CLOSED_CONTAINER_OPEN_CHILD",
                         "container": parent})
        for parent in stale_open:
            rows.append({"event": "info_row", "class": "OPEN_CONTAINER_ALL_CHILDREN_CLOSED",
                         "container": parent})
        for parent in orphans:
            rows.append({"event": "conflict_row", "class": "MISSING_PARENT",
                         "child_prefix": parent})
        for iid, entries in sorted(flags.items()):
            for name, blocking in entries:
                rows.append({"event": "conflict_row" if blocking else "info_row",
                             "class": name, "id": iid})
        declared = {
            "corpus_count": len(corpus),
            "hierarchy_duplication_count": len(dup),
            "closed_container_open_child_count": len(contras),
            "stale_open_container_count": len(stale_open),
            "missing_parent_count": len(orphans),
            "owner_missing_count": sum(1 for e in flags.values() if any(n == "OWNER_MISSING" for n, _ in e)),
            "estimate_missing_count": sum(1 for e in flags.values() if any(n == "ESTIMATE_MISSING" for n, _ in e)),
        }
        derived = {
            "corpus_count": sum(1 for r in rows if r["event"] == "contract_record"),
            "hierarchy_duplication_count": sum(
                1 for r in rows if r.get("class") == "HIERARCHY_DUPLICATION"),
            "closed_container_open_child_count": sum(
                1 for r in rows if r.get("class") == "CLOSED_CONTAINER_OPEN_CHILD"),
            "stale_open_container_count": sum(
                1 for r in rows if r.get("class") == "OPEN_CONTAINER_ALL_CHILDREN_CLOSED"),
            "missing_parent_count": sum(
                1 for r in rows if r.get("class") == "MISSING_PARENT"),
            "owner_missing_count": sum(
                1 for r in rows if r.get("class") == "OWNER_MISSING"),
            "estimate_missing_count": sum(
                1 for r in rows if r.get("class") == "ESTIMATE_MISSING"),
        }
        for field, expect in declared.items():
            if expect != derived[field]:
                print(
                    f"count/member drift on {field}: declared {expect} "
                    f"vs emitted {derived[field]}",
                    file=sys.stderr,
                )
                sys.exit(33)
        rows.append({"event": "registry-summary", **declared})
        rows.append({
            "event": "no-claim",
            "text": "structural conflicts are review findings, not verdicts; "
                    "this registry grants no execution, scientific, or "
                    "closeout authority and does not normalize any clause",
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
            print(f"registry written: {out_arg}")
        else:
            sys.stdout.write(rendered)
        return 0

    if mode == "check":
        try:
            with open(out_arg) as handle:
                prior = [json.loads(line) for line in handle if line.strip()]
        except (OSError, json.JSONDecodeError) as error:
            print(f"registry refusal: {error}", file=sys.stderr)
            sys.exit(30)
        header = next((r for r in prior if r.get("event") == "registry-header"), None)
        if not header or header.get("schema") != SCHEMA:
            print("registry refusal: missing or unknown schema header", file=sys.stderr)
            sys.exit(30)
        prior_corpus = {r["id"]: r for r in prior if r.get("event") == "contract_record"}
        live_fields = ("title", "status", "priority", "issue_type",
                       "assignee", "estimated_minutes")
        drift = []
        for iid in sorted(set(prior_corpus) & set(issues)):
            for field in live_fields:
                if prior_corpus[iid].get(field) != issues[iid].get(field):
                    drift.append((iid, field))
        gone = sorted(set(prior_corpus) - set(issues))
        new = sorted(corpus - set(prior_corpus))
        if drift or gone or new:
            print(
                f"registry STALE against live candidate set: "
                f"{len(drift)} field drift(s), {len(gone)} vanished, "
                f"{len(new)} new; first divergence {drift[:1] or gone[:1] or new[:1]}; "
                "re-freeze the registry deliberately",
                file=sys.stderr,
            )
            sys.exit(31)
        print(
            f"registry OK: {len(prior_corpus)} records bound to live DB "
            f"identity {db_sha[:16]}"
        )
        return 0

    print(f"unknown mode {mode}", file=sys.stderr)
    return 30


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:4]))
# ---CORE-END---
PYEOF
