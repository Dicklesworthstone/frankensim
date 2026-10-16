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
# Modes (--apply mutates ONLY through the exact br commands of a reviewed
# plan artifact; every other mode is read-only):
#   --list             enumerate live candidate edges, run nothing
#   --freeze [OUT]     write the content-addressed review manifest (all
#                      rows verdict=pending) to OUT, defaulting to the
#                      tracked tests/leapfrog/manifests location
#   --check MANIFEST   validate a manifest's schema, counts, and DB binding
#   --plan MANIFEST    emit the exact two-step br command sequence for rows
#                      reviewed verdict=migrate, plus the inverse plan;
#                      refuses if ANY row is still pending
#   --apply PLAN       execute ONLY the plan artifact's exact br dep
#                      commands (remove FIRST, then add), halting on the
#                      first unexpected failure; log under .e2e-out/
#   --negative         fixture battery proving every replay/refusal class
#   --self-test        exercise the tool's own refusals on fixtures
#
# Migration shape (per reviewed relation, only via --apply over a plan):
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
EXIT_MISMATCH=33

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
  --negative) MODE="negative" ;;
  --check|--plan)
    MODE="${1#--}"
    [ $# -ge 2 ] || die "$EXIT_USAGE" "$1 needs a manifest path"
    MANIFEST_ARG="$2" ;;
  --inverse-plan)
    MODE="inverse-plan"
    [ $# -ge 2 ] || die "$EXIT_USAGE" "--inverse-plan needs the plan artifact path"
    MANIFEST_ARG="$2"; INVERSE_OUT="${3:-}" ;;
  --replay)
    MODE="replay"
    [ $# -ge 4 ] || die "$EXIT_USAGE" "--replay needs PLAN PRE POST [OUT]"
    MANIFEST_ARG="$2"; REPLAY_PRE="$3"; REPLAY_POST="$4"; REPLAY_OUT="${5:-}" ;;
  *) die "$EXIT_USAGE" "usage: --list | --freeze | --check MANIFEST | --plan MANIFEST | --apply PLAN | --inverse-plan PLAN [OUT] | --replay PLAN PRE POST [OUT] | --self-test | --negative" ;;
esac

tool() { # subcommand [manifest]
  python3 - "$1" "$BEADS_FILE" "$SCHEMA" "${@:2}" <<'PYEOF'
import hashlib
import json
import sys
import os
from collections import Counter

mode, beads_path, schema = sys.argv[1:4]
manifest_path = sys.argv[4] if len(sys.argv) > 4 else ""
extra = sys.argv[5:]


def refuse_overlap(out_path, inputs):
    """Artifact-only guarantee: replay output must land outside the live
    tracker export and must never alias an immutable input."""
    out_real = os.path.realpath(out_path)
    beads_real = os.path.realpath(beads_path)
    if out_real == beads_real or out_real.startswith(beads_real + os.sep):
        print(
            f"path-overlap refusal: {out_path} overlaps the tracker export",
            file=sys.stderr,
        )
        sys.exit(EXIT_USAGE)
    for source in inputs:
        if source and out_real == os.path.realpath(source):
            print(
                f"path-overlap refusal: {out_path} is an input path",
                file=sys.stderr,
            )
            sys.exit(EXIT_USAGE)


def edge_multiset(deps):
    return Counter(
        (dep["issue_id"], dep["depends_on_id"], dep.get("type", ""))
        for dep in deps
        if isinstance(dep.get("issue_id"), str)
        and isinstance(dep.get("depends_on_id"), str)
    )


def load_plan_for_execution(plan_path, refusal_prefix):
    try:
        plan = json.load(open(plan_path))
    except (OSError, json.JSONDecodeError) as error:
        print(f"{refusal_prefix} refusal: plan unreadable: {error}", file=sys.stderr)
        sys.exit(EXIT_USAGE)
    if plan.get("schema") != "frankensim.beads-hierarchy-normalization-plan.v1":
        print(
            f"{refusal_prefix} refusal: plan schema {plan.get('schema')!r}",
            file=sys.stderr,
        )
        sys.exit(EXIT_USAGE)
    return plan


def validate_br_dep_commands(commands, refusal_prefix):
    for command in commands:
        if not (isinstance(command, list) and command[:2] == ["br", "dep"]
                and command[2] in ("remove", "add")):
            print(
                f"{refusal_prefix} refusal: non-br-dep command in plan: {command}",
                file=sys.stderr,
            )
            sys.exit(EXIT_USAGE)

EXIT_USAGE = 30
EXIT_STALE = 31
EXIT_UNREVIEWED = 32
EXIT_MISMATCH = 33


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

# check / plan need a reviewed MANIFEST; inverse-plan and replay consume
# PLAN artifacts instead and validate their own schema below.
if mode in ("check", "plan"):
    try:
        manifest = json.load(open(manifest_path))
    except (OSError, json.JSONDecodeError) as error:
        print(f"manifest refusal: {error}", file=sys.stderr)
        sys.exit(EXIT_USAGE)
    if manifest.get("schema") != schema:
        print(
            f"manifest refusal: schema {manifest.get('schema')!r}",
            file=sys.stderr,
        )
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
    admitted = {
        "pending",
        "migrate",
        "migrate-remove-only",
        "keep-functional-prerequisite",
    }
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
        if row["verdict"] not in ("migrate", "migrate-remove-only"):
            continue
        if (row["parent"], row["child"]) not in live:
            print(
                f"plan refusal: reviewed edge {row['parent']} -> {row['child']} "
                "no longer exists live; re-freeze",
                file=sys.stderr,
            )
            sys.exit(EXIT_STALE)
        # Remove FIRST, then add - never both edges at once. A row whose
        # child->parent parent-child edge ALREADY exists (double encoding)
        # is remove-only: adding again would duplicate the hierarchy edge.
        plan.append(["br", "dep", "remove", row["parent"], row["child"]])
        inverse_add = ["br", "dep", "add", row["parent"], row["child"], "--type", "blocks"]
        if row["verdict"] == "migrate":
            plan.append(
                ["br", "dep", "add", row["child"], row["parent"], "--type", "parent-child"]
            )
            inverse.append(["br", "dep", "remove", row["child"], row["parent"]])
        inverse.append(inverse_add)
    print(
        json.dumps(
            {
                "schema": "frankensim.beads-hierarchy-normalization-plan.v1",
                "manifest_sha256": hashlib.sha256(
                    json.dumps(manifest, sort_keys=True).encode()
                ).hexdigest(),
                "migrate_rows": sum(1 for r in mrows if r["verdict"] == "migrate"),
                "remove_only_rows": sum(
                    1 for r in mrows if r["verdict"] == "migrate-remove-only"
                ),
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

if mode == "inverse-plan":
    plan = load_plan_for_execution(manifest_path, "inverse-plan")
    inverse = plan.get("inverse_commands")
    if not isinstance(inverse, list) or not inverse:
        print(
            "inverse-plan refusal: plan carries no retained inverse_commands",
            file=sys.stderr,
        )
        sys.exit(EXIT_USAGE)
    validate_br_dep_commands(inverse, "inverse-plan")
    artifact = {
        "schema": "frankensim.beads-hierarchy-normalization-inverse-plan.v1",
        "source_plan_sha256": hashlib.sha256(
            open(manifest_path, "rb").read()
        ).hexdigest(),
        "command_count": len(inverse),
        "commands": inverse,
        "note": (
            "retained explicit br commands for rollback only; "
            "never auto-applied by this tool"
        ),
    }
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    if extra and extra[0]:
        refuse_overlap(extra[0], [manifest_path])
        with open(extra[0], "w") as handle:
            handle.write(rendered)
        print(f"inverse plan written: {extra[0]} ({len(inverse)} commands)")
    else:
        print(rendered, end="")
    sys.exit(0)

if mode == "replay":
    # Artifact-only adjudication: simulate the retained ledger over the
    # pre-state snapshot and compare against the post-state snapshot.
    # Reads immutable inputs, writes one disjoint receipt, never touches br.
    if len(extra) < 2:
        print("replay refusal: --replay needs PLAN PRE POST [OUT]", file=sys.stderr)
        sys.exit(EXIT_USAGE)
    pre_path, post_path = extra[0], extra[1]
    out_path = extra[2] if len(extra) > 2 and extra[2] else None
    for source in (pre_path, post_path):
        if not os.path.isfile(source):
            print(f"replay refusal: snapshot not found: {source}", file=sys.stderr)
            sys.exit(EXIT_USAGE)
    plan = load_plan_for_execution(manifest_path, "replay")
    commands = plan.get("commands") or []
    if not commands:
        print("replay refusal: plan has no commands", file=sys.stderr)
        sys.exit(EXIT_USAGE)
    validate_br_dep_commands(commands, "replay")
    _, pre_issues, pre_deps = load(pre_path)
    _, post_issues, post_deps = load(post_path)
    pre_edges = edge_multiset(pre_deps)
    post_edges = edge_multiset(post_deps)

    sim = Counter(pre_edges)
    removal_seen = set()
    ordering_ok = True
    first_divergence = None
    for index, command in enumerate(commands):
        actor, target, kind = command[3], command[4], command[2]
        if kind == "remove":
            key = (actor, target, "blocks")
            if sim.get(key, 0) <= 0:
                first_divergence = {
                    "index": index,
                    "class": "remove-missing-edge",
                    "edge": [actor, target],
                }
                break
            sim[key] -= 1
            if not sim[key]:
                del sim[key]
            removal_seen.add((actor, target))
        else:
            key = (actor, target, "parent-child")
            if sim.get(key, 0) > 0:
                first_divergence = {
                    "index": index,
                    "class": "add-duplicate-edge",
                    "edge": [actor, target],
                }
                break
            if (target, actor) not in removal_seen:
                ordering_ok = False
                first_divergence = {
                    "index": index,
                    "class": "add-before-remove",
                    "edge": [target, actor],
                }
                break
            sim[key] += 1

    planned_removals = sorted({(c[3], c[4]) for c in commands if c[2] == "remove"})
    planned_adds = sorted({(c[3], c[4]) for c in commands if c[2] == "add"})
    residue = [
        [parent, child]
        for parent, child in planned_removals
        if post_edges.get((parent, child, "blocks"), 0) > 0
    ]
    add_faults = [
        [child, parent, post_edges.get((child, parent, "parent-child"), 0)]
        for child, parent in planned_adds
        if post_edges.get((child, parent, "parent-child"), 0) != 1
    ]
    drift_added = Counter(post_edges) - Counter(sim)
    drift_removed = Counter(sim) - Counter(post_edges)

    def slim(row):
        return {
            k: v
            for k, v in row.items()
            if k not in ("dependencies", "updated_at")
        }

    changed = [
        issue_id
        for issue_id, row in post_issues.items()
        if issue_id in pre_issues and slim(pre_issues[issue_id]) != slim(row)
    ]
    added_ids = [i for i in post_issues if i not in pre_issues]
    removed_ids = [i for i in pre_issues if i not in post_issues]

    projection_pass = (
        first_divergence is None
        and ordering_ok
        and not residue
        and not add_faults
    )

    def bounded(items, cap=20):
        return sorted(list(item) if isinstance(item, tuple) else item for item in items)[:cap]

    rows = [
        {
            "event": "replay-header",
            "schema": "frankensim.beads-hierarchy-normalization-replay-receipt.v1",
            "plan_sha256": hashlib.sha256(
                open(manifest_path, "rb").read()
            ).hexdigest(),
            "pre_sha256": hashlib.sha256(open(pre_path, "rb").read()).hexdigest(),
            "post_sha256": hashlib.sha256(open(post_path, "rb").read()).hexdigest(),
            "commands": len(commands),
        },
        {
            "event": "stage",
            "name": "ledger-simulation",
            "verdict": "pass" if first_divergence is None else "fail",
            "ordering_remove_before_add": ordering_ok,
            "first_divergence": first_divergence,
        },
        {
            "event": "stage",
            "name": "target-old-edges-absent",
            "verdict": "pass" if not residue else "fail",
            "checked": len(planned_removals),
            "residue_sample": bounded(residue),
            "residue_count": len(residue),
        },
        {
            "event": "stage",
            "name": "target-new-edges-present-once",
            "verdict": "pass" if not add_faults else "fail",
            "checked": len(planned_adds),
            "fault_sample": bounded(add_faults),
            "fault_count": len(add_faults),
        },
        {
            "event": "stage",
            "name": "non-target-edge-drift",
            "verdict": "informational",
            "post_extra_edges": sum(drift_added.values()),
            "post_missing_edges": sum(drift_removed.values()),
            "post_extra_sample": bounded(drift_added.elements()),
            "post_missing_sample": bounded(drift_removed.elements()),
        },
        {
            "event": "stage",
            "name": "non-target-issue-fields",
            "verdict": "informational",
            "changed": len(changed),
            "added_issues": len(added_ids),
            "removed_issues": len(removed_ids),
            "changed_sample": sorted(changed)[:20],
            "added_sample": sorted(added_ids)[:20],
            "removed_sample": sorted(removed_ids)[:20],
        },
        {
            "event": "replay-verdict",
            "projection_verdict": "pass" if projection_pass else "fail",
            "exit_code": 0 if projection_pass else EXIT_MISMATCH,
            "reproduction_command": (
                f"scripts/e2e/leapfrog/beads-hierarchy-normalization.sh "
                f"--replay {manifest_path} {pre_path} {post_path}"
                + (f" {out_path}" if out_path else "")
            ),
        },
    ]
    rendered = "".join(json.dumps(r, sort_keys=True) + "\n" for r in rows)
    if out_path:
        refuse_overlap(out_path, [manifest_path, pre_path, post_path])
        with open(out_path, "w") as handle:
            handle.write(rendered)
        print(f"replay receipt written: {out_path}")
    else:
        print(rendered, end="")
    sys.exit(0 if projection_pass else EXIT_MISMATCH)

print(f"unknown mode {mode}", file=sys.stderr)
sys.exit(EXIT_USAGE)
PYEOF
}

run_negative() {
  SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/hiernorm-negative.XXXXXX")"
  trap 'rm -rf "$SCRATCH"' EXIT
  PASS=0; FAIL=0
  check() {
    if [ "$2" -eq "$3" ]; then PASS=$((PASS + 1)); else
      FAIL=$((FAIL + 1)); echo "NEGATIVE FAIL: $1 (expected $2, got $3)" >&2
    fi
  }
  # PRE: epic e blocks its own children e.1 and e.2; e.2 already carries
  # the child->parent parent-child edge (remove-only row); an unrelated
  # functional edge e->x must survive untouched.
  printf '%s\n' \
    '{"id":"e","status":"open","priority":1,"dependencies":[{"depends_on_id":"e.1","type":"blocks"},{"depends_on_id":"e.2","type":"blocks"},{"depends_on_id":"x","type":"blocks"}]}' \
    '{"id":"e.1","status":"open","priority":2,"description":"child one"}' \
    '{"id":"e.2","status":"open","priority":3,"description":"child two","dependencies":[{"depends_on_id":"e","type":"parent-child"}]}' \
    '{"id":"x","status":"open","description":"functional"}' \
    > "$SCRATCH/pre.jsonl"
  # POST: both old edges gone; e.1 gained the hierarchy edge; non-target
  # fields identical apart from legitimate updated_at churn.
  printf '%s\n' \
    '{"id":"e","status":"open","priority":1,"updated_at":"2026-08-14T18:00:00Z","dependencies":[{"depends_on_id":"x","type":"blocks"}]}' \
    '{"id":"e.1","status":"open","priority":2,"description":"child one","updated_at":"2026-08-14T18:00:00Z","dependencies":[{"depends_on_id":"e","type":"parent-child"}]}' \
    '{"id":"e.2","status":"open","priority":3,"description":"child two","dependencies":[{"depends_on_id":"e","type":"parent-child"}]}' \
    '{"id":"x","status":"open","description":"functional"}' \
    > "$SCRATCH/post.jsonl"
  # PRE missing one planned removal target; POST missing the planned add;
  # POST retaining one old edge.
  printf '%s\n' \
    '{"id":"e","status":"open","dependencies":[{"depends_on_id":"e.2","type":"blocks"},{"depends_on_id":"x","type":"blocks"}]}' \
    '{"id":"e.1","status":"open"}' \
    '{"id":"e.2","status":"open"}' \
    > "$SCRATCH/pre_missing.jsonl"
  printf '%s\n' \
    '{"id":"e","status":"open","dependencies":[{"depends_on_id":"x","type":"blocks"}]}' \
    '{"id":"e.1","status":"open"}' \
    '{"id":"e.2","status":"open","dependencies":[{"depends_on_id":"e","type":"parent-child"}]}' \
    > "$SCRATCH/post_missing_add.jsonl"
  printf '%s\n' \
    '{"id":"e","status":"open","dependencies":[{"depends_on_id":"e.2","type":"blocks"},{"depends_on_id":"x","type":"blocks"}]}' \
    '{"id":"e.1","status":"open","dependencies":[{"depends_on_id":"e","type":"parent-child"}]}' \
    '{"id":"e.2","status":"open","dependencies":[{"depends_on_id":"e","type":"parent-child"}]}' \
    > "$SCRATCH/post_residue.jsonl"
  python3 - "$SCRATCH" <<'PYFIX'
import json, sys
scratch = sys.argv[1]
base = {
    "schema": "frankensim.beads-hierarchy-normalization-plan.v1",
    "manifest_sha256": "fixture",
    "commands": [
        ["br", "dep", "remove", "e", "e.1"],
        ["br", "dep", "add", "e.1", "e", "--type", "parent-child"],
        ["br", "dep", "remove", "e", "e.2"],
    ],
    "inverse_commands": [
        ["br", "dep", "remove", "e.1", "e"],
        ["br", "dep", "add", "e", "e.1", "--type", "blocks"],
    ],
}
json.dump(base, open(scratch + "/plan.json", "w"))
bad_schema = dict(base); bad_schema["schema"] = "other.v1"
json.dump(bad_schema, open(scratch + "/plan_bad_schema.json", "w"))
non_br = dict(base); non_br["commands"] = [["cargo", "test"]] + base["commands"]
json.dump(non_br, open(scratch + "/plan_nonbr.json", "w"))
add_first = dict(base); add_first["commands"] = list(reversed(base["commands"]))
json.dump(add_first, open(scratch + "/plan_addfirst.json", "w"))
PYFIX
  ENV=(FSIM_HIERNORM_BEADS="$SCRATCH/pre.jsonl")
  env "${ENV[@]}" "$0" --replay "$SCRATCH/plan.json" "$SCRATCH/pre.jsonl" "$SCRATCH/post.jsonl" "$SCRATCH/receipt.jsonl" >/dev/null 2>&1
  check "clean synthetic migration replays PASS" 0 $?
  grep -q '"projection_verdict": "pass"' "$SCRATCH/receipt.jsonl" 2>/dev/null
  check "PASS receipt records projection_verdict=pass" 0 $?
  env "${ENV[@]}" "$0" --replay "$SCRATCH/plan_bad_schema.json" "$SCRATCH/pre.jsonl" "$SCRATCH/post.jsonl" >/dev/null 2>&1
  check "tampered plan schema refuses (30)" 30 $?
  env "${ENV[@]}" "$0" --replay "$SCRATCH/plan_nonbr.json" "$SCRATCH/pre.jsonl" "$SCRATCH/post.jsonl" >/dev/null 2>&1
  check "non-br-dep ledger command refuses (30)" 30 $?
  env "${ENV[@]}" "$0" --replay "$SCRATCH/plan_addfirst.json" "$SCRATCH/pre.jsonl" "$SCRATCH/post.jsonl" >/dev/null 2>&1
  check "add-before-remove ordering fails closed (33)" 33 $?
  env "${ENV[@]}" "$0" --replay "$SCRATCH/plan.json" "$SCRATCH/pre_missing.jsonl" "$SCRATCH/post.jsonl" >/dev/null 2>&1
  check "PRE missing a planned removal mismatches (33)" 33 $?
  env "${ENV[@]}" "$0" --replay "$SCRATCH/plan.json" "$SCRATCH/pre.jsonl" "$SCRATCH/post_missing_add.jsonl" >/dev/null 2>&1
  check "POST missing a planned add mismatches (33)" 33 $?
  env "${ENV[@]}" "$0" --replay "$SCRATCH/plan.json" "$SCRATCH/pre.jsonl" "$SCRATCH/post_residue.jsonl" >/dev/null 2>&1
  check "POST retaining an old edge mismatches (33)" 33 $?
  env "${ENV[@]}" "$0" --replay "$SCRATCH/plan.json" "$SCRATCH/pre.jsonl" "$SCRATCH/post.jsonl" "$SCRATCH/pre.jsonl" >/dev/null 2>&1
  check "output aliasing an input refuses (30)" 30 $?
  env "${ENV[@]}" "$0" --inverse-plan "$SCRATCH/plan_bad_schema.json" >/dev/null 2>&1
  check "inverse-plan refuses unknown schema (30)" 30 $?
  env "${ENV[@]}" "$0" --inverse-plan "$SCRATCH/plan.json" "$SCRATCH/inverse.json" >/dev/null 2>&1
  check "inverse-plan emits retained artifact" 0 $?
  python3 - "$SCRATCH/inverse.json" <<'PYINV'
import json, sys
artifact = json.load(open(sys.argv[1]))
assert artifact["schema"] == "frankensim.beads-hierarchy-normalization-inverse-plan.v1"
assert artifact["command_count"] == len(artifact["commands"]) == 2
PYINV
  check "inverse-plan artifact carries exact retained commands" 0 $?
  echo "negative: $PASS passed, $FAIL failed"
  [ "$FAIL" -eq 0 ] || exit 1
  exit 0
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
  apply)
    PLAN_FILE="$MANIFEST_ARG"
    [ -f "$PLAN_FILE" ] || die "$EXIT_USAGE" "plan artifact not found: $PLAN_FILE"
    command -v br >/dev/null 2>&1 || die "$EXIT_USAGE" "br is required for --apply"
    # The apply executes ONLY the plan artifact's exact commands - never a
    # recomputed set - and halts on the first unexpected failure with a
    # pointer to the retained inverse plan.
    STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
    APPLY_LOG="${FSIM_HIERNORM_APPLY_LOG:-$REPO_ROOT/.e2e-out/hiernorm-apply-$STAMP.jsonl}"
    mkdir -p "$(dirname "$APPLY_LOG")"
    python3 - "$PLAN_FILE" "$APPLY_LOG" <<'PYAPPLY'
import hashlib
import json
import subprocess
import sys
import datetime

plan_path, log_path = sys.argv[1], sys.argv[2]
import os
start = int(os.environ.get("FSIM_HIERNORM_APPLY_START", "0"))
stop = int(os.environ.get("FSIM_HIERNORM_APPLY_STOP", "0")) or None
plan = json.load(open(plan_path))
if plan.get("schema") != "frankensim.beads-hierarchy-normalization-plan.v1":
    print(f"apply refusal: plan schema {plan.get('schema')!r}", file=sys.stderr)
    sys.exit(30)
commands = plan.get("commands") or []
if not commands:
    print("apply refusal: plan has no commands", file=sys.stderr)
    sys.exit(30)
for command in commands:
    if not (isinstance(command, list) and command[:2] == ["br", "dep"]
            and command[2] in ("remove", "add")):
        print(f"apply refusal: non-br-dep command in plan: {command}", file=sys.stderr)
        sys.exit(30)

plan_digest = hashlib.sha256(open(plan_path, "rb").read()).hexdigest()
seq = 0


def emit(event, **fields):
    global seq
    seq += 1
    row = {"schema": "frankensim.beads-hierarchy-normalization-apply-log.v1",
           "seq": seq, "plan_sha256": plan_digest,
           "at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
           "event": event}
    row.update(fields)
    with open(log_path, "a") as fh:
        fh.write(json.dumps(row, sort_keys=True) + "\n")


emit("apply-start", command_count=len(commands), start=start, stop=stop or len(commands))
applied = 0
for index, command in enumerate(commands):
    if index < start or (stop is not None and index >= stop):
        continue
    result = subprocess.run(command, capture_output=True, text=True, timeout=120)
    ok = result.returncode == 0
    emit("command", index=index, argv=command, exit=result.returncode,
         stderr_tail=result.stderr.strip()[-200:] if not ok else "")
    if not ok:
        emit("apply-halt", applied=applied, failed_index=index)
        print(
            f"apply HALTED at command {index}/{len(commands)} "
            f"({' '.join(command)}): {result.stderr.strip()[-300:]}\n"
            f"{applied} command(s) applied; the inverse plan in the artifact "
            f"covers rollback; log: {log_path}",
            file=sys.stderr,
        )
        sys.exit(31)
    applied += 1
    if applied % 200 == 0:
        print(f"  progress: {applied}/{len(commands)}")
emit("apply-complete", applied=applied)
print(f"apply complete: {applied}/{len(commands)} commands; log: {log_path}")
PYAPPLY
    APPLY_STATUS=$?
    [ "$APPLY_STATUS" -eq 0 ] || exit "$APPLY_STATUS"
    # Post-condition: re-discover live candidates; the migrated set must be
    # gone. A nonzero residue is loud, not silent.
    RESIDUE="$(tool list 2>/dev/null | wc -l | tr -d ' ')"
    echo "post-apply live candidate edges: $RESIDUE (expected 0 if the whole plan applied)"
    ;;
  inverse-plan)
    PLAN_FILE="$MANIFEST_ARG"
    [ -f "$PLAN_FILE" ] || die "$EXIT_USAGE" "plan artifact not found: $PLAN_FILE"
    tool inverse-plan "$PLAN_FILE" "$INVERSE_OUT"
    ;;
  replay)
    PLAN_FILE="$MANIFEST_ARG"
    [ -f "$PLAN_FILE" ] || die "$EXIT_USAGE" "plan artifact not found: $PLAN_FILE"
    [ -f "$REPLAY_PRE" ] || die "$EXIT_USAGE" "pre-state snapshot not found: $REPLAY_PRE"
    [ -f "$REPLAY_POST" ] || die "$EXIT_USAGE" "post-state snapshot not found: $REPLAY_POST"
    tool replay "$PLAN_FILE" "$REPLAY_PRE" "$REPLAY_POST" "$REPLAY_OUT"
    ;;
  negative) run_negative ;;
esac
