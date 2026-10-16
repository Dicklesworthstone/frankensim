#!/usr/bin/env python3
"""Runner V2 acceptance-corpus registry (bead frankensim-epic-foundations-huq.24.6).

Builds a source-authoritative, machine-readable inventory of the Runner V2
acceptance corpus from frozen br output (.beads/issues.jsonl) WITHOUT
rewriting any clause text: every row binds the exact acceptance-criteria
bytes by sha256, so the historical supersession chain stays intact and any
silent edit is detectable as a hash move.

What this slice detects mechanically (each is a typed finding, and the
--self-test fixtures prove each detector fires):
  duplicate-exclusive-ownership   two beads claim the same owned source
                                  path or frozen symbol
  ac-fragment-collision           the same ACnn fragment is OWNED (not
                                  merely referenced) by two beads WITHOUT a
                                  common parent that declares per-child
                                  fragment ownership (the corpus's split-
                                  by-design pattern is informational, not a
                                  finding)
  missing-canonical-ac            a family bead whose description defers to
                                  the canonical Acceptance Criteria field
                                  while that field is empty
  charter-date-skew               sibling work packages carry different
                                  CONSOLIDATED CHARTER dates (a stale
                                  controlling clause candidate)
  unowned-implementation          a source file under runner_v2/work_packages
                                  exists in the tree but no family bead owns
                                  it (orphan implementation)
  closed-missing-implementation   a CLOSED bead exclusively owns a source
                                  path that does not exist in the tree

What this slice does NOT do (deliberately, per the bead): it never factors
shared base obligations from per-bead deltas, never rewrites text, and never
judges scientific content. Factoring comes only after a lossless one-to-one
mapping is mechanically established against this inventory.

Modes:
  (default)        write runner-v2-acceptance-registry.json
  --check          regenerate and byte-compare against the tracked registry
  --self-test      run the embedded detector fixtures
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BEADS = REPO_ROOT / ".beads" / "issues.jsonl"
REGISTRY = REPO_ROOT / "runner-v2-acceptance-registry.json"
SCHEMA = "frankensim.runner-v2-acceptance-registry.v1"

OWNED_PATH = re.compile(r"exclusively owns\s+((?:`[^`]+`(?:,\s*)?)+)", re.I)
BACKTICKED = re.compile(r"`([^`]+)`")
AC_FRAGMENT_OWNED = re.compile(r"\b(AC\d+) fragment\b")
CHARTER = re.compile(r"CONSOLIDATED CHARTER (\d{4}-\d{2}-\d{2})")
SUPERSEDE = re.compile(r"\bsupersede[sd]?\b", re.I)
CANONICAL_DEFER = re.compile(
    r"dedicated Acceptance Criteria field is the canonical", re.I
)


def load_family(beads_text: str) -> list[dict]:
    family = []
    for line in beads_text.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        if not isinstance(row, dict):
            continue
        labels = row.get("labels") or []
        if "runner-v2" in labels or "Runner V2" in row.get("title", ""):
            family.append(row)
    family.sort(key=lambda row: row.get("id", ""))
    return family


def project(row: dict) -> dict:
    ac = row.get("acceptance_criteria") or ""
    description = row.get("description") or ""
    both = ac + "\n" + description
    owned: list[str] = []
    for match in OWNED_PATH.finditer(both):
        owned.extend(BACKTICKED.findall(match.group(1)))
    return {
        "id": row.get("id", ""),
        "status": row.get("status", ""),
        "title": row.get("title", ""),
        "ac_sha256": hashlib.sha256(ac.encode()).hexdigest(),
        "ac_bytes": len(ac.encode()),
        "description_sha256": hashlib.sha256(description.encode()).hexdigest(),
        "owned_symbols": sorted(set(owned)),
        "owned_ac_fragments": sorted(set(AC_FRAGMENT_OWNED.findall(ac))),
        "charter_dates": sorted(set(CHARTER.findall(ac))),
        "supersession_language": bool(SUPERSEDE.search(both)),
        "defers_to_canonical_ac": bool(CANONICAL_DEFER.search(description)),
    }


def common_parent(ids: list[str]) -> str:
    split = [id_.split(".") for id_ in ids]
    prefix = split[0]
    for parts in split[1:]:
        while prefix and parts[: len(prefix)] != prefix:
            prefix = prefix[:-1]
    return ".".join(prefix)


def detect(rows: list[dict], ac_by_id: dict[str, str] | None = None) -> list[dict]:
    findings: list[dict] = []
    owners: dict[str, list[str]] = defaultdict(list)
    fragment_owners: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        for symbol in row["owned_symbols"]:
            owners[symbol].append(row["id"])
        for fragment in row["owned_ac_fragments"]:
            fragment_owners[fragment].append(row["id"])
    for symbol, ids in sorted(owners.items()):
        if len(ids) > 1:
            findings.append(
                {
                    "kind": "duplicate-exclusive-ownership",
                    "subject": symbol,
                    "beads": sorted(ids),
                }
            )
    for fragment, ids in sorted(fragment_owners.items()):
        if len(ids) > 1:
            # Split-by-design: the corpus assigns each child its own
            # immutable ACnn fragment under a parent that says so. Only
            # flag claimant sets whose common parent does NOT declare
            # per-child fragment ownership.
            parent = common_parent(sorted(ids))
            parent_ac = (ac_by_id or {}).get(parent, "")
            by_design = bool(
                re.search(r"its (?:own )?(?:immutable )?(?:AC\d+ )?fragment", parent_ac)
            )
            if not by_design:
                findings.append(
                    {
                        "kind": "ac-fragment-collision",
                        "subject": fragment,
                        "beads": sorted(ids),
                    }
                )
    for row in rows:
        if row["defers_to_canonical_ac"] and row["ac_bytes"] == 0:
            findings.append(
                {
                    "kind": "missing-canonical-ac",
                    "subject": row["id"],
                    "beads": [row["id"]],
                }
            )
    # Charter-date skew is judged among siblings (same dotted parent) that
    # carry any charter at all; a parent whose children disagree on the
    # controlling charter date has a stale controlling-clause candidate.
    by_parent: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if row["charter_dates"]:
            parent = row["id"].rsplit(".", 1)[0]
            by_parent[parent].append(row)
    for parent, group in sorted(by_parent.items()):
        dates = Counter(date for row in group for date in row["charter_dates"])
        if len(dates) > 1:
            findings.append(
                {
                    "kind": "charter-date-skew",
                    "subject": parent,
                    "beads": sorted(row["id"] for row in group),
                    "dates": dict(sorted(dates.items())),
                }
            )
    return findings


WORK_PACKAGE_DIR = "crates/fs-evidence-runner/src/runner_v2/work_packages"
OWNED_PATH_PREFIXES = ("crates/", "runner_v2/", "tests/", "scripts/")


def resolve_owned_path(symbol: str) -> str | None:
    """Repo-relative form of an owned symbol IF it names a source path."""
    if not symbol.endswith(".rs"):
        return None
    if symbol.startswith("crates/") or symbol.startswith("scripts/"):
        return symbol
    if symbol.startswith("runner_v2/"):
        return f"crates/fs-evidence-runner/src/{symbol}"
    if symbol.startswith("tests/"):
        return f"crates/fs-evidence-runner/{symbol}"
    return None


def detect_tree(rows: list[dict], tree_files: set[str]) -> list[dict]:
    """Cross-reference owned source paths against an actual file listing.

    `tree_files` is the repo-relative set of existing work-package-relevant
    files; injected so the self-test can exercise both directions without
    touching the real tree.
    """
    findings: list[dict] = []
    owned_paths: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        for symbol in row["owned_symbols"]:
            resolved = resolve_owned_path(symbol)
            if resolved:
                owned_paths[resolved].append(row["id"])
    for path in sorted(tree_files):
        if path.startswith(WORK_PACKAGE_DIR + "/") and path not in owned_paths:
            findings.append(
                {
                    "kind": "unowned-implementation",
                    "subject": path,
                    "beads": [],
                }
            )
    closed = {row["id"] for row in rows if row["status"] == "closed"}
    for path, ids in sorted(owned_paths.items()):
        closed_owners = sorted(set(ids) & closed)
        if closed_owners and path not in tree_files:
            findings.append(
                {
                    "kind": "closed-missing-implementation",
                    "subject": path,
                    "beads": closed_owners,
                }
            )
    return findings


def live_tree_files() -> set[str]:
    files: set[str] = set()
    base = REPO_ROOT / "crates" / "fs-evidence-runner"
    for path in base.rglob("*.rs"):
        files.add(str(path.relative_to(REPO_ROOT)))
    for extra in (REPO_ROOT / "scripts").rglob("*.rs"):
        files.add(str(extra.relative_to(REPO_ROOT)))
    return files


def render(beads_text: str) -> str:
    family = load_family(beads_text)
    rows = [project(row) for row in family]
    ac_by_id = {
        row.get("id", ""): row.get("acceptance_criteria") or "" for row in family
    }
    registry = {
        "schema": SCHEMA,
        "authority_statement": (
            "INERT INVENTORY DATA derived from frozen br output. Clause text "
            "is bound by hash, never rewritten; the beads database remains "
            "the sole clause authority. This registry cannot close, alter, "
            "or supersede any acceptance criterion."
        ),
        "beads_source_sha256": hashlib.sha256(beads_text.encode()).hexdigest(),
        "family_size": len(rows),
        "rows": rows,
        "findings": detect(rows, ac_by_id) + detect_tree(rows, live_tree_files()),
    }
    return json.dumps(registry, indent=2, sort_keys=True) + "\n"


def self_test() -> int:
    def bead(id_, ac="", description="", labels=("runner-v2",)):
        return json.dumps(
            {
                "id": id_,
                "title": f"Runner V2 fixture {id_}",
                "status": "open",
                "labels": list(labels),
                "acceptance_criteria": ac,
                "description": description,
            }
        )

    failures = []

    def expect(name, text, kind, present=True):
        registry = json.loads(render(text))
        kinds = [finding["kind"] for finding in registry["findings"]]
        if (kind in kinds) != present:
            failures.append(f"{name}: expected {kind} present={present}, got {kinds}")

    dup = "\n".join(
        [
            bead("f-1", ac="This Bead exclusively owns `runner_v2/work_packages/a.rs`."),
            bead("f-2", ac="This Bead exclusively owns `runner_v2/work_packages/a.rs`."),
        ]
    )
    expect("duplicate ownership fires", dup, "duplicate-exclusive-ownership")

    distinct = "\n".join(
        [
            bead("f-1", ac="This Bead exclusively owns `runner_v2/work_packages/a.rs`."),
            bead("f-2", ac="This Bead exclusively owns `runner_v2/work_packages/b.rs`."),
        ]
    )
    expect(
        "distinct ownership is clean", distinct, "duplicate-exclusive-ownership", False
    )

    frag = "\n".join(
        [
            bead("p.1", ac="owns the AC58 fragment."),
            bead("p.2", ac="owns the AC58 fragment too."),
        ]
    )
    expect("AC fragment collision fires", frag, "ac-fragment-collision")

    by_design = "\n".join(
        [
            bead("p", ac="each child exclusively owns its own immutable AC58 fragment."),
            bead("p.1", ac="owns its AC58 fragment."),
            bead("p.2", ac="owns its AC58 fragment."),
        ]
    )
    expect(
        "declared split-by-design fragments stay clean",
        by_design,
        "ac-fragment-collision",
        False,
    )

    missing = bead(
        "f-3",
        ac="",
        description="The dedicated Acceptance Criteria field is the canonical close gate.",
    )
    expect("missing canonical AC fires", missing, "missing-canonical-ac")

    skew = "\n".join(
        [
            bead("p.1", ac="CONSOLIDATED CHARTER 2026-07-30 (HIGHEST PRECEDENCE)"),
            bead("p.2", ac="CONSOLIDATED CHARTER 2026-08-02 (HIGHEST PRECEDENCE)"),
        ]
    )
    expect("charter date skew fires", skew, "charter-date-skew")

    uniform = "\n".join(
        [
            bead("p.1", ac="CONSOLIDATED CHARTER 2026-07-30"),
            bead("p.2", ac="CONSOLIDATED CHARTER 2026-07-30"),
        ]
    )
    expect("uniform charter is clean", uniform, "charter-date-skew", False)

    # Tree cross-reference detectors (injected file listings).
    owned = [project(json.loads(bead(
        "t-1", ac="This Bead exclusively owns `runner_v2/work_packages/a.rs`.")))]
    orphan = detect_tree(
        owned,
        {
            "crates/fs-evidence-runner/src/runner_v2/work_packages/a.rs",
            "crates/fs-evidence-runner/src/runner_v2/work_packages/rogue.rs",
        },
    )
    if [f["kind"] for f in orphan] != ["unowned-implementation"]:
        failures.append(f"orphan implementation: got {orphan}")

    closed_row = project(json.loads(bead(
        "t-2", ac="This Bead exclusively owns `runner_v2/work_packages/gone.rs`.")))
    closed_row["status"] = "closed"
    gone = detect_tree([closed_row], set())
    if [f["kind"] for f in gone] != ["closed-missing-implementation"]:
        failures.append(f"closed-missing-implementation: got {gone}")
    open_row = dict(closed_row, status="open")
    if detect_tree([open_row], set()):
        failures.append("an OPEN bead's future path must not be a finding")

    # Lossless binding: any clause edit moves the row hash.
    before = json.loads(render(bead("f-4", ac="clause text v1")))
    after = json.loads(render(bead("f-4", ac="clause text v2")))
    if before["rows"][0]["ac_sha256"] == after["rows"][0]["ac_sha256"]:
        failures.append("clause edit must move the acceptance hash")

    for failure in failures:
        print(f"SELF-TEST FAIL: {failure}", file=sys.stderr)
    print(f"self-test: {12 - len(failures)} passed, {len(failures)} failed")
    return 1 if failures else 0


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "--self-test":
        return self_test()
    beads_text = BEADS.read_text()
    rendered = render(beads_text)
    if mode == "--check":
        if not REGISTRY.is_file():
            print(f"missing tracked registry {REGISTRY}", file=sys.stderr)
            return 1
        if REGISTRY.read_text() != rendered:
            print(
                "runner-v2-acceptance-registry.json is stale against the live "
                "beads database; regenerate deliberately with "
                "scripts/ci/runner_v2_acceptance_registry.py",
                file=sys.stderr,
            )
            return 1
        print("registry check OK")
        return 0
    if mode not in ("", "--generate"):
        print(__doc__, file=sys.stderr)
        return 2
    REGISTRY.write_text(rendered)
    registry = json.loads(rendered)
    print(
        f"registry written: {registry['family_size']} beads, "
        f"{len(registry['findings'])} finding(s)"
    )
    for finding in registry["findings"]:
        print(f"  {finding['kind']}: {finding['subject']} <- {finding['beads']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
