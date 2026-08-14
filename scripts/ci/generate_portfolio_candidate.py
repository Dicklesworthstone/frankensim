#!/usr/bin/env python3
"""Emit the CandidatePortfolioAuthoritySnapshot (bead f85xj.16.11).

Produces `candidate-portfolio-snapshot.json`: the machine-readable candidate
projection of PORTFOLIO_DISPOSITION_CONTRACT_V1 (transcribed verbatim from
the owning bead's canonical Design field), bound to the exact live Beads
snapshot it was derived from. Every count is recomputed here from
`.beads/issues.jsonl` — no copied literal drives authority, per the bead's
own acceptance rule.

THE CANDIDATE IS INERT. It cannot change robot selection, public capability
projection, active budgets, or promotion authority. Only the independent
16.12 adjudication may reconstruct the decision from bounded inputs and
atomically install an ActivePortfolioAuthoritySnapshot; every non-Pass
preserves the exact prior active snapshot (today: none exists).

Producer-side validation (16.11's scope) runs inside this generator and its
results are embedded in the artifact: disposition-target existence, owner
existence, exactly-one-primary, role-vocabulary closure, and duplicate
detection. A validation failure refuses to emit rather than emitting a
partially trusted candidate.
"""

import datetime
import hashlib
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BEADS = REPO / ".beads" / "issues.jsonl"
OUT = REPO / "candidate-portfolio-snapshot.json"

PRIMARY_SPINE = "frankensim-extreal-program-f85xj.6"
PRIMARY_ACCEPTANCE = "frankensim-extreal-program-f85xj.6.11"
ACTIVATION_OWNER = "frankensim-extreal-program-f85xj.16.12"

# PORTFOLIO_DISPOSITION_CONTRACT_V1, transcribed from the canonical Design
# field of frankensim-extreal-program-f85xj.16.11. The bead is the source of
# truth for these rows; edit THERE first, then re-transcribe.
DISPOSITIONS = [
    {
        "issue_id": "frankensim-ext-e0d-machine-slice-7jxk",
        "spine_role": "shared_infrastructure",
        "disposition": "historical_evidence",
        "integration_owner": "frankensim-ext-epic-e0d-36db",
        "acceptance_owner": "frankensim-ext-epic-gov-rjoq.8",
        "attention_budget_class": "shared_20",
        "compute_budget_class": "shared_20",
        "wip_cap": 1,
        "falsifier": "machine_graph_admission_or_audit_or_replay_identity_fails_or_requires_duplicate_primary_authority",
        "promotion_gate": "shared_only_after_e0d_replay_pass_and_16_12_active_snapshot_pass",
        "forbidden_promotion": "primary_product_or_primary_acceptance",
    },
    {
        "issue_id": "frankensim-ext-flagship-coldplate-2kkb",
        "spine_role": "primary_enhancement",
        "disposition": "secondary_candidate",
        "integration_owner": "frankensim-ext-epic-e7-3nxa",
        "acceptance_owner": "frankensim-ext-epic-e7-3nxa.1",
        "attention_budget_class": "primary_70_after_cooling_l3",
        "compute_budget_class": "primary_70_after_cooling_l3",
        "wip_cap": 1,
        "falsifier": "conjugate_heat_pressure_energy_or_supplier_step_receipts_fail_the_frozen_reference_envelopes",
        "promotion_gate": "cooling_6_11_pass_then_e7_battery_pass_then_16_12_active_snapshot_pass",
        "forbidden_promotion": "l4_or_l5_by_inference",
    },
]

_FRONTIER = [
    ("frankensim-ext-flagship-geneva-b4hj", "frankensim-ext-epic-e2-q8wx",
     "frankensim-ext-e2-test-battery-xxgv", "e2",
     "intermittent_event_timing_retention_contact_or_wear_receipt_fails"),
    ("frankensim-ext-flagship-gear-ei5t", "frankensim-ext-epic-e3-e1c9",
     "frankensim-ext-e3-test-battery-ztso", "e3",
     "transmission_error_contact_life_or_acoustic_receipt_fails"),
    ("frankensim-ext-flagship-motor-m8n3", "frankensim-ext-epic-e4-tlca",
     "frankensim-ext-e4-test-battery-5kv0", "e4",
     "torque_power_loss_demagnetization_thermal_or_control_closure_fails"),
    ("frankensim-ext-flagship-wankel-jwx0", "frankensim-ext-epic-e5-3n3f",
     "frankensim-ext-e5-test-battery-u06c", "e5",
     "trochoid_seal_contact_mass_energy_chemistry_or_event_balance_fails"),
    ("frankensim-ext-flagship-genset-mr28", "frankensim-ext-epic-e7-3nxa",
     "frankensim-ext-epic-e7-3nxa.1", "e7",
     "coupled_ice_generator_energy_speed_control_or_interface_closure_fails"),
    ("frankensim-ext-flagship-digital-twin-pkfe", "frankensim-ext-epic-e7-3nxa",
     "frankensim-ext-epic-e7-3nxa.1", "e7",
     "posterior_coverage_drift_detection_escalation_or_replay_identity_fails"),
    ("frankensim-ext-flagship-turbo-efuel-7zyp", "frankensim-ext-epic-e7-3nxa",
     "frankensim-ext-epic-e7-3nxa.1", "e7",
     "boost_chemistry_knock_emissions_or_coupled_balance_receipt_fails"),
    ("frankensim-ext-flagship-pump-bearing-z65l", "frankensim-ext-epic-e7-3nxa",
     "frankensim-ext-epic-e7-3nxa.1", "e7",
     "head_efficiency_cavitation_fluid_film_stability_or_life_receipt_fails"),
    ("frankensim-ext-flagship-induction-i40z", "frankensim-ext-epic-e7-3nxa",
     "frankensim-ext-epic-e7-3nxa.1", "e7",
     "team_torque_loss_thermal_moving_conductor_or_control_closure_fails"),
    ("frankensim-ext-flagship-constant-width-cvng", "frankensim-ext-epic-e7-3nxa",
     "frankensim-ext-epic-e7-3nxa.1", "e7",
     "rolling_no_slip_contact_geometry_or_motion_falsifier_fails"),
]
for issue, integ, accept, battery, falsifier in _FRONTIER:
    DISPOSITIONS.append({
        "issue_id": issue,
        "spine_role": "secondary_candidate",
        "disposition": "capped_frontier",
        "integration_owner": integ,
        "acceptance_owner": accept,
        "attention_budget_class": "frontier_10_collective",
        "compute_budget_class": "frontier_10_collective",
        "wip_cap": 1,
        "falsifier": falsifier,
        "promotion_gate": (
            "cooling_6_11_pass_and_allocation_exception_then_"
            f"{battery}_battery_pass_then_16_12_active_snapshot_pass"
        ),
        "forbidden_promotion": "primary_or_l4_l5_by_inference",
    })

AGGREGATE_CONSTRAINTS = [
    "exactly_one_primary_product_and_acceptance_owner",
    "e0d_never_projects_as_product_authority",
    "coldplate_cannot_block_cooling_before_6_11_pass",
    "capped_frontier_rows_share_one_collective_10_percent_budget_and_one_active_integration_wip_per_phase",
    "every_promotion_requires_nonloss_inventory_and_16_12_activation",
    "no_row_implies_l4_or_l5",
]

ROLE_VOCABULARY = {
    "primary_product_spine", "primary_enhancement", "prerequisite_consumed",
    "shared_infrastructure", "capped_research_or_second_vertical",
    "secondary_candidate", "deferred_by_gate", "superseded_or_absorbed",
    "historical_evidence",
}


def load_issues():
    issues = {}
    for line in BEADS.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(row, dict) and row.get("id"):
            issues[row["id"]] = row
    return issues


def descendants(issues, root):
    """Parent-child closure under `root` (dependency type parent-child)."""
    children = {}
    for issue in issues.values():
        for dep in issue.get("dependencies") or []:
            if isinstance(dep, dict) and dep.get("type") == "parent-child":
                children.setdefault(dep.get("depends_on"), []).append(issue["id"])
    seen, stack = set(), [root]
    while stack:
        node = stack.pop()
        for child in children.get(node, []):
            if child not in seen:
                seen.add(child)
                stack.append(child)
    return seen


def main():
    issues = load_issues()
    failures = []
    checks = []

    def check(name, ok, detail):
        checks.append({"check": name, "pass": bool(ok), "detail": detail})
        if not ok:
            failures.append(f"{name}: {detail}")

    # Producer-side completeness validation (16.11 scope).
    ids = [d["issue_id"] for d in DISPOSITIONS]
    check("no-duplicate-disposition-rows", len(ids) == len(set(ids)),
          f"{len(ids)} rows, {len(set(ids))} unique")
    for d in DISPOSITIONS:
        for key in ("issue_id", "integration_owner", "acceptance_owner"):
            check(f"exists:{d[key]}", d[key] in issues,
                  f"{key} of {d['issue_id']} resolves in the live graph")
        check(f"role-vocabulary:{d['issue_id']}",
              d["spine_role"] in ROLE_VOCABULARY, d["spine_role"])
    check("exists:primary-spine", PRIMARY_SPINE in issues, PRIMARY_SPINE)
    check("exists:primary-acceptance", PRIMARY_ACCEPTANCE in issues, PRIMARY_ACCEPTANCE)
    check("exists:activation-owner", ACTIVATION_OWNER in issues, ACTIVATION_OWNER)
    primary_roles = [d for d in DISPOSITIONS if d["spine_role"] == "primary_product_spine"]
    check("exactly-one-primary",
          len(primary_roles) == 0,  # the primary is named at the top level, never a row
          "primary identity lives in primary_spine_id; no disposition row may claim it")

    if failures:
        sys.stderr.write("REFUSED: producer validation failed; no candidate emitted\n")
        for f in failures:
            sys.stderr.write(f"  {f}\n")
        return 1

    # Live-graph binding: every count recomputed from the snapshot.
    status = {}
    for issue in issues.values():
        status[issue.get("status", "?")] = status.get(issue.get("status", "?"), 0) + 1
    on_path = descendants(issues, PRIMARY_SPINE)
    on_path_open = sum(1 for i in on_path
                       if issues.get(i, {}).get("status") in ("open", "in_progress", "blocked"))
    head = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                          capture_output=True, text=True).stdout.strip()
    br_version = subprocess.run(["br", "--version"], capture_output=True,
                                text=True).stdout.strip() or "br (version unavailable)"
    data_root = hashlib.sha256(BEADS.read_bytes()).hexdigest()

    snapshot = {
        "schema": "extreal.portfolio-disposition.v1",
        "kind": "CandidatePortfolioAuthoritySnapshot",
        "authority_statement": (
            "INERT CANDIDATE DATA. This snapshot cannot change robot selection, "
            "public capability projection, active budgets, or promotion authority. "
            "Only frankensim-extreal-program-f85xj.16.12 may independently adjudicate "
            "non-loss, uniqueness, dependencies, budgets, and user projection, and "
            "alone may atomically install an ActivePortfolioAuthoritySnapshot; every "
            "non-Pass preserves the exact prior active snapshot (currently: none)."
        ),
        "primary_spine_id": PRIMARY_SPINE,
        "primary_acceptance_owner": PRIMARY_ACCEPTANCE,
        "activation_owner": ACTIVATION_OWNER,
        "allocation_basis": "separate_human_attention_and_compute_70_20_10",
        "default_transition_rule": "no_role_change_without_nonloss_inventory_and_16_12_pass",
        "dispositions": DISPOSITIONS,
        "aggregate_constraints": AGGREGATE_CONSTRAINTS,
        "live_graph_binding": {
            "captured_at": datetime.datetime.now(datetime.timezone.utc)
                .strftime("%Y-%m-%dT%H:%M:%SZ"),
            "beads_data_root_sha256": data_root,
            "head_sha": head,
            "br_version": br_version,
            "derivation": "all counts recomputed from .beads/issues.jsonl at the recorded root; no copied literal",
            "issue_total": len(issues),
            "status_counts": dict(sorted(status.items())),
            "primary_spine_descendants": len(on_path),
            "primary_spine_open_or_active": on_path_open,
        },
        "producer_validation": {
            "producer_bead": "frankensim-extreal-program-f85xj.16.11",
            "checks_run": len(checks),
            "checks_failed": 0,
            "checks": checks,
        },
        "history_preservation": (
            "No disposition deletes or collapses ambitions: absorbed and historical "
            "records retain their features, children, proofs, no-claims, and lineage "
            "in their own bead records, which this snapshot references but never rewrites."
        ),
    }
    OUT.write_text(json.dumps(snapshot, indent=2) + "\n")
    print(f"candidate emitted: {OUT.name} ({len(DISPOSITIONS)} dispositions, "
          f"{len(checks)} producer checks green, data root {data_root[:16]})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
