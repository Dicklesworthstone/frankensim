//! G0 battery for the generated requirement-to-evidence traceability ledger.

use std::collections::BTreeSet;

use fs_govern::{
    MAX_PROOF_OBLIGATION_OWNERS, MAX_REQUIREMENT_PO_LINKS, MAX_REQUIREMENT_ROWS,
    MAX_TRACEABILITY_FIELD_BYTES, PROOF_OBLIGATION_COUNT, ProofObligation, RequirementRow,
    TRACEABILITY_AUTHORITY, TRACEABILITY_SCHEMA, TraceabilityAudit, TraceabilityField,
    audit_traceability, generate_traceability_ledger, proof_obligations, requirements,
    traceability_ledger_json,
};

fn assert_field_failure(
    row: RequirementRow<'_>,
    expected_scope: &str,
    expected_field: TraceabilityField,
) -> TraceabilityAudit {
    let audit = generate_traceability_ledger(&[row], proof_obligations()).unwrap_err();
    assert!(
        audit.diagnostics.iter().any(|diagnostic| {
            diagnostic.requirement_id == expected_scope && diagnostic.field == expected_field
        }),
        "expected {expected_scope}.{} failure, got {:?}",
        expected_field.name(),
        audit.diagnostics
    );
    audit
}

#[test]
fn canonical_registry_has_the_exact_closed_requirement_and_po_inventories() {
    let expected_requirements: BTreeSet<&str> = [
        "B1",
        "B2",
        "B3",
        "B4",
        "B5",
        "B6",
        "B7",
        "B8",
        "B9",
        "B10",
        "B11",
        "B12",
        "B13",
        "B14",
        "RQ-ROLL",
        "RQ-GEAR",
        "RQ-FRICTION",
        "RQ-CONSTITUTIVE",
        "RQ-DENSITY",
        "RQ-MECHMAT",
        "RQ-ELEC",
        "RQ-MAG",
        "RQ-PHASE",
        "RQ-FLUID",
        "RQ-PERMEATE",
        "RQ-WET",
        "RQ-MOTORGEN",
        "RQ-ICE",
        "RQ-ACOUSTIC",
        "RQ-ACTIVE",
    ]
    .into_iter()
    .collect();
    let actual_requirements: BTreeSet<&str> = requirements()
        .iter()
        .map(|row| row.requirement_id)
        .collect();
    assert_eq!(requirements().len(), 30);
    assert_eq!(actual_requirements, expected_requirements);

    assert_eq!(proof_obligations().len(), PROOF_OBLIGATION_COUNT);
    let actual_pos: Vec<&str> = proof_obligations()
        .iter()
        .map(|obligation| obligation.id)
        .collect();
    let expected_pos: Vec<String> = (1..=PROOF_OBLIGATION_COUNT)
        .map(|number| format!("PO-{number}"))
        .collect();
    assert_eq!(
        actual_pos,
        expected_pos.iter().map(String::as_str).collect::<Vec<_>>()
    );
}

#[test]
fn canonical_registry_is_complete_and_every_requirement_links_known_pos() {
    let audit = audit_traceability(requirements(), proof_obligations());
    assert!(
        audit.ok(),
        "canonical registry gaps: {:?}",
        audit.diagnostics
    );
    assert_eq!(audit.total, 30);
    assert_eq!(audit.complete, 30);

    let known: BTreeSet<&str> = proof_obligations()
        .iter()
        .map(|obligation| obligation.id)
        .collect();
    for row in requirements() {
        assert!(
            !row.proof_obligations.is_empty(),
            "{} has no proof gate",
            row.requirement_id
        );
        for id in row.proof_obligations {
            assert!(
                known.contains(id),
                "{} links unknown {id}",
                row.requirement_id
            );
        }
    }
}

#[test]
#[allow(clippy::too_many_lines)] // one table pins every mandatory source field
fn every_required_row_field_has_a_named_fail_closed_diagnostic() {
    let base = requirements()[0];
    let cases = [
        (
            RequirementRow {
                requirement_id: "",
                ..base
            },
            "<orphaned-requirement>",
            TraceabilityField::RequirementId,
        ),
        (
            RequirementRow {
                capability_property: "",
                ..base
            },
            "B1",
            TraceabilityField::CapabilityProperty,
        ),
        (
            RequirementRow {
                blocker: "",
                ..base
            },
            "B1",
            TraceabilityField::Blocker,
        ),
        (
            RequirementRow {
                owner_artifact: "",
                ..base
            },
            "B1",
            TraceabilityField::OwnerArtifact,
        ),
        (
            RequirementRow {
                prerequisite_phase: "",
                ..base
            },
            "B1",
            TraceabilityField::PrerequisitePhase,
        ),
        (
            RequirementRow {
                milestone: "",
                ..base
            },
            "B1",
            TraceabilityField::Milestone,
        ),
        (
            RequirementRow {
                flagship: "",
                ..base
            },
            "B1",
            TraceabilityField::Flagship,
        ),
        (
            RequirementRow {
                benchmark_data: "",
                ..base
            },
            "B1",
            TraceabilityField::BenchmarkData,
        ),
        (
            RequirementRow {
                proof_obligations: &[],
                ..base
            },
            "B1",
            TraceabilityField::ProofObligation,
        ),
        (
            RequirementRow {
                claim_boundary: "",
                ..base
            },
            "B1",
            TraceabilityField::ClaimBoundary,
        ),
        (
            RequirementRow { status: "", ..base },
            "B1",
            TraceabilityField::Status,
        ),
    ];

    for (row, scope, field) in cases {
        let _ = assert_field_failure(row, scope, field);
    }
}

#[test]
fn deliberately_orphaned_requirement_names_the_id_field_and_reason() {
    let row = RequirementRow {
        requirement_id: "RQ-ORPHAN",
        owner_artifact: "",
        ..requirements()[0]
    };
    let audit = assert_field_failure(row, "RQ-ORPHAN", TraceabilityField::OwnerArtifact);
    let diagnostic = audit
        .diagnostics
        .iter()
        .find(|diagnostic| diagnostic.field == TraceabilityField::OwnerArtifact)
        .unwrap();
    assert!(diagnostic.reason.contains("no owner or artifact route"));
}

#[test]
fn tracker_closed_is_not_accepted_as_scientific_status() {
    let row = RequirementRow {
        status: "closed",
        ..requirements()[0]
    };
    let audit = assert_field_failure(row, "B1", TraceabilityField::Status);
    assert!(audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "B1"
            && diagnostic
                .reason
                .contains("cannot mint scientific proof status")
    }));

    for promoted in ["verified", "validated"] {
        let audit = assert_field_failure(
            RequirementRow {
                status: promoted,
                ..requirements()[0]
            },
            "B1",
            TraceabilityField::Status,
        );
        assert!(
            audit
                .diagnostics
                .iter()
                .any(|diagnostic| diagnostic.reason.contains("unbound registry")),
            "{promoted} must require a source-bound promotion receipt"
        );
    }
}

#[test]
fn duplicate_requirements_and_dangling_or_repeated_po_links_refuse() {
    let base = requirements()[0];
    let duplicate = generate_traceability_ledger(&[base, base], proof_obligations()).unwrap_err();
    assert!(duplicate.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "B1"
            && diagnostic.field == TraceabilityField::RequirementId
            && diagnostic.reason.contains("duplicate")
    }));

    let dangling = RequirementRow {
        proof_obligations: &["PO-26"],
        ..base
    };
    let dangling_audit = assert_field_failure(dangling, "B1", TraceabilityField::ProofObligation);
    assert!(
        dangling_audit.diagnostics[0]
            .reason
            .contains("unknown proof obligation")
    );

    let repeated = RequirementRow {
        proof_obligations: &["PO-5", "PO-5"],
        ..base
    };
    let repeated_audit = assert_field_failure(repeated, "B1", TraceabilityField::ProofObligation);
    assert!(
        repeated_audit
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.reason.contains("more than once"))
    );

    let blank = RequirementRow {
        proof_obligations: &[""],
        ..base
    };
    let blank_audit = assert_field_failure(blank, "B1", TraceabilityField::ProofObligation);
    assert!(blank_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "B1"
            && diagnostic.reason.contains("empty proof-obligation link")
    }));
}

#[test]
fn closed_requirement_inventory_rejects_missing_and_extra_ids() {
    let mut missing = requirements().to_vec();
    let removed = missing.pop().unwrap();
    let missing_audit = generate_traceability_ledger(&missing, proof_obligations()).unwrap_err();
    assert!(missing_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == removed.requirement_id
            && diagnostic.field == TraceabilityField::RequirementId
            && diagnostic.reason.contains("closed registry is missing")
    }));

    let mut extra = requirements().to_vec();
    extra.push(RequirementRow {
        requirement_id: "RQ-UNREGISTERED",
        ..requirements()[0]
    });
    let extra_audit = generate_traceability_ledger(&extra, proof_obligations()).unwrap_err();
    assert!(extra_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "RQ-UNREGISTERED"
            && diagnostic.field == TraceabilityField::RequirementId
            && diagnostic.reason.contains("outside the closed")
    }));
}

#[test]
#[allow(clippy::too_many_lines)] // one battery pins the complete PO-definition refusal surface
fn po_index_rejects_missing_alias_duplicate_and_ownerless_definitions() {
    let base = requirements()[0];
    let missing =
        generate_traceability_ledger(&[base], &proof_obligations()[..PROOF_OBLIGATION_COUNT - 1])
            .unwrap_err();
    assert!(
        missing
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.reason.contains("missing PO-25"))
    );

    let mut alias = proof_obligations().to_vec();
    alias[0] = ProofObligation {
        id: "PO-01",
        ..alias[0]
    };
    let alias_audit = generate_traceability_ledger(&[base], &alias).unwrap_err();
    assert!(
        alias_audit
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.reason.contains("outside the closed"))
    );
    assert!(
        alias_audit
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.reason.contains("missing PO-1"))
    );

    let mut duplicate_owner = proof_obligations().to_vec();
    duplicate_owner[0] = ProofObligation {
        owner_beads: &["owner-a", "owner-a"],
        ..duplicate_owner[0]
    };
    let duplicate_owner_audit =
        generate_traceability_ledger(&[base], &duplicate_owner).unwrap_err();
    assert!(duplicate_owner_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "PO-1"
            && diagnostic.field == TraceabilityField::OwnerArtifact
            && diagnostic.reason.contains("duplicate")
    }));

    let mut ownerless = proof_obligations().to_vec();
    ownerless[0] = ProofObligation {
        owner_beads: &[],
        ..ownerless[0]
    };
    let ownerless_audit = generate_traceability_ledger(&[base], &ownerless).unwrap_err();
    assert!(ownerless_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "PO-1"
            && diagnostic.field == TraceabilityField::OwnerArtifact
            && diagnostic.reason.contains("no owning Bead")
    }));

    let mut summaryless = proof_obligations().to_vec();
    summaryless[0] = ProofObligation {
        summary: "",
        ..summaryless[0]
    };
    let summaryless_audit = generate_traceability_ledger(&[base], &summaryless).unwrap_err();
    assert!(summaryless_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "PO-1"
            && diagnostic.field == TraceabilityField::ProofObligation
            && diagnostic.reason.contains("missing its executable summary")
    }));

    let mut duplicate_definition = proof_obligations().to_vec();
    duplicate_definition[1] = ProofObligation {
        id: "PO-1",
        ..duplicate_definition[1]
    };
    let duplicate_definition_audit =
        generate_traceability_ledger(requirements(), &duplicate_definition).unwrap_err();
    assert!(
        duplicate_definition_audit
            .diagnostics
            .iter()
            .any(|diagnostic| {
                diagnostic.requirement_id == "PO-1"
                    && diagnostic
                        .reason
                        .contains("duplicate proof-obligation definition")
            })
    );

    let mut idless = proof_obligations().to_vec();
    idless[0] = ProofObligation {
        id: "",
        ..idless[0]
    };
    let idless_audit = generate_traceability_ledger(requirements(), &idless).unwrap_err();
    assert!(idless_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "<proof-obligation-index>"
            && diagnostic.reason.contains("missing its id")
    }));

    let mut blank_owner = proof_obligations().to_vec();
    blank_owner[0] = ProofObligation {
        owner_beads: &[""],
        ..blank_owner[0]
    };
    let blank_owner_audit = generate_traceability_ledger(requirements(), &blank_owner).unwrap_err();
    assert!(blank_owner_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "PO-1"
            && diagnostic.field == TraceabilityField::OwnerArtifact
            && diagnostic.reason.contains("no owning Bead")
    }));
}

#[test]
#[allow(clippy::too_many_lines)] // exact and plus-one collection boundaries stay adjacent
fn empty_scope_and_collection_caps_refuse_without_crossing_the_boundary() {
    let empty = generate_traceability_ledger(&[], proof_obligations()).unwrap_err();
    assert!(empty.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "<ledger>"
            && diagnostic.field == TraceabilityField::Ledger
            && diagnostic.reason.contains("empty requirement scope")
    }));

    let at_row_cap = vec![requirements()[0]; MAX_REQUIREMENT_ROWS];
    let at_row_cap_audit = audit_traceability(&at_row_cap, proof_obligations());
    assert!(
        !at_row_cap_audit
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.reason.contains("source contains")),
        "the exact row cap must reach semantic validation"
    );
    let above_row_cap = vec![requirements()[0]; MAX_REQUIREMENT_ROWS + 1];
    let above_row_cap_audit = audit_traceability(&above_row_cap, proof_obligations());
    assert_eq!(above_row_cap_audit.complete, 0);
    assert!(above_row_cap_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "<ledger>" && diagnostic.reason.contains("source contains")
    }));

    let mut above_po_count = proof_obligations().to_vec();
    above_po_count.push(proof_obligations()[0]);
    let above_po_count_audit = audit_traceability(requirements(), &above_po_count);
    assert_eq!(above_po_count_audit.complete, 0);
    assert!(above_po_count_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "<proof-obligation-index>"
            && diagnostic.reason.contains("index contains")
    }));

    let all_po_links: Vec<String> = (1..=PROOF_OBLIGATION_COUNT)
        .map(|number| format!("PO-{number}"))
        .collect();
    let all_po_link_refs: Vec<&str> = all_po_links.iter().map(String::as_str).collect();
    assert_eq!(all_po_link_refs.len(), MAX_REQUIREMENT_PO_LINKS);
    let mut exact_link_rows = requirements().to_vec();
    let exact_link_base = exact_link_rows[0];
    exact_link_rows[0] = RequirementRow {
        proof_obligations: &all_po_link_refs,
        ..exact_link_base
    };
    assert!(
        generate_traceability_ledger(&exact_link_rows, proof_obligations()).is_ok(),
        "the exact per-row PO-link cap is valid"
    );
    let mut above_link_refs = all_po_link_refs.clone();
    above_link_refs.push("PO-1");
    let mut above_link_rows = requirements().to_vec();
    let above_link_base = above_link_rows[0];
    above_link_rows[0] = RequirementRow {
        proof_obligations: &above_link_refs,
        ..above_link_base
    };
    let above_link_audit =
        generate_traceability_ledger(&above_link_rows, proof_obligations()).unwrap_err();
    assert!(above_link_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "B1"
            && diagnostic.reason.contains("proof obligations; maximum")
    }));

    let owner_strings: Vec<String> = (0..MAX_PROOF_OBLIGATION_OWNERS)
        .map(|index| format!("owner-{index}"))
        .collect();
    let owner_refs: Vec<&str> = owner_strings.iter().map(String::as_str).collect();
    let mut exact_owner_index = proof_obligations().to_vec();
    let exact_owner_base = exact_owner_index[0];
    exact_owner_index[0] = ProofObligation {
        owner_beads: &owner_refs,
        ..exact_owner_base
    };
    assert!(
        generate_traceability_ledger(requirements(), &exact_owner_index).is_ok(),
        "the exact PO-owner cap is valid"
    );
    let extra_owner = String::from("owner-over-cap");
    let mut above_owner_refs = owner_refs.clone();
    above_owner_refs.push(&extra_owner);
    let mut above_owner_index = proof_obligations().to_vec();
    let above_owner_base = above_owner_index[0];
    above_owner_index[0] = ProofObligation {
        owner_beads: &above_owner_refs,
        ..above_owner_base
    };
    let above_owner_audit =
        generate_traceability_ledger(requirements(), &above_owner_index).unwrap_err();
    assert!(above_owner_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "PO-1" && diagnostic.reason.contains("owners; maximum")
    }));
}

#[test]
#[allow(clippy::too_many_lines)] // exact and plus-one byte boundaries stay adjacent
fn scalar_and_reference_byte_caps_refuse_at_limit_plus_one() {
    let exact_field = "x".repeat(MAX_TRACEABILITY_FIELD_BYTES);
    let above_field = "x".repeat(MAX_TRACEABILITY_FIELD_BYTES + 1);
    let mut exact_rows = requirements().to_vec();
    let exact_row_base = exact_rows[0];
    exact_rows[0] = RequirementRow {
        capability_property: &exact_field,
        ..exact_row_base
    };
    assert!(
        generate_traceability_ledger(&exact_rows, proof_obligations()).is_ok(),
        "the exact scalar-byte cap is valid"
    );
    let mut above_rows = requirements().to_vec();
    let above_row_base = above_rows[0];
    above_rows[0] = RequirementRow {
        capability_property: &above_field,
        ..above_row_base
    };
    let above_field_audit =
        generate_traceability_ledger(&above_rows, proof_obligations()).unwrap_err();
    assert!(above_field_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "B1"
            && diagnostic.field == TraceabilityField::CapabilityProperty
            && diagnostic.reason.contains("maximum")
    }));

    let exact_requirement_id = "R".repeat(MAX_TRACEABILITY_FIELD_BYTES);
    let above_requirement_id = "R".repeat(MAX_TRACEABILITY_FIELD_BYTES + 1);
    let mut exact_id_rows = requirements().to_vec();
    let exact_requirement_base = exact_id_rows[0];
    exact_id_rows[0] = RequirementRow {
        requirement_id: &exact_requirement_id,
        ..exact_requirement_base
    };
    let exact_requirement_audit = audit_traceability(&exact_id_rows, proof_obligations());
    assert!(
        !exact_requirement_audit
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.requirement_id == "<oversized-requirement-id>"),
        "the exact requirement-id cap must reach closed-id validation"
    );
    let mut above_id_rows = requirements().to_vec();
    let above_requirement_base = above_id_rows[0];
    above_id_rows[0] = RequirementRow {
        requirement_id: &above_requirement_id,
        ..above_requirement_base
    };
    let above_requirement_audit = audit_traceability(&above_id_rows, proof_obligations());
    assert!(
        above_requirement_audit
            .diagnostics
            .iter()
            .any(|diagnostic| {
                diagnostic.requirement_id == "<oversized-requirement-id>"
                    && diagnostic.field == TraceabilityField::RequirementId
                    && diagnostic.reason.contains("maximum")
            })
    );

    let exact_summary = "s".repeat(MAX_TRACEABILITY_FIELD_BYTES);
    let above_summary = "s".repeat(MAX_TRACEABILITY_FIELD_BYTES + 1);
    let mut exact_summary_index = proof_obligations().to_vec();
    let exact_summary_base = exact_summary_index[0];
    exact_summary_index[0] = ProofObligation {
        summary: &exact_summary,
        ..exact_summary_base
    };
    assert!(
        generate_traceability_ledger(requirements(), &exact_summary_index).is_ok(),
        "the exact PO-summary cap is valid"
    );
    let mut above_summary_index = proof_obligations().to_vec();
    let above_summary_base = above_summary_index[0];
    above_summary_index[0] = ProofObligation {
        summary: &above_summary,
        ..above_summary_base
    };
    let above_summary_audit =
        generate_traceability_ledger(requirements(), &above_summary_index).unwrap_err();
    assert!(above_summary_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "PO-1"
            && diagnostic.reason.contains("summary is")
            && diagnostic.reason.contains("maximum")
    }));

    let exact_owner = "o".repeat(MAX_TRACEABILITY_FIELD_BYTES);
    let above_owner = "o".repeat(MAX_TRACEABILITY_FIELD_BYTES + 1);
    let exact_owner_refs = [exact_owner.as_str()];
    let above_owner_refs = [above_owner.as_str()];
    let mut exact_owner_index = proof_obligations().to_vec();
    let exact_owner_base = exact_owner_index[0];
    exact_owner_index[0] = ProofObligation {
        owner_beads: &exact_owner_refs,
        ..exact_owner_base
    };
    assert!(
        generate_traceability_ledger(requirements(), &exact_owner_index).is_ok(),
        "the exact owner-id cap is valid"
    );
    let mut above_owner_index = proof_obligations().to_vec();
    let above_owner_base = above_owner_index[0];
    above_owner_index[0] = ProofObligation {
        owner_beads: &above_owner_refs,
        ..above_owner_base
    };
    let above_owner_audit =
        generate_traceability_ledger(requirements(), &above_owner_index).unwrap_err();
    assert!(above_owner_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "PO-1"
            && diagnostic.reason.contains("owner Bead is")
            && diagnostic.reason.contains("maximum")
    }));

    let exact_long_id = "I".repeat(MAX_TRACEABILITY_FIELD_BYTES);
    let above_long_id = "I".repeat(MAX_TRACEABILITY_FIELD_BYTES + 1);
    let mut exact_id_index = proof_obligations().to_vec();
    let exact_id_base = exact_id_index[0];
    exact_id_index[0] = ProofObligation {
        id: &exact_long_id,
        ..exact_id_base
    };
    let exact_id_audit = audit_traceability(requirements(), &exact_id_index);
    assert!(
        !exact_id_audit.diagnostics.iter().any(|diagnostic| {
            diagnostic.reason.contains("id is") && diagnostic.reason.contains("maximum")
        }),
        "the exact PO-id byte cap must reach canonical-id validation"
    );
    let mut above_id_index = proof_obligations().to_vec();
    let above_id_base = above_id_index[0];
    above_id_index[0] = ProofObligation {
        id: &above_long_id,
        ..above_id_base
    };
    let above_id_audit = audit_traceability(requirements(), &above_id_index);
    assert!(above_id_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "<oversized-proof-obligation-id>"
            && diagnostic.reason.contains("id is")
            && diagnostic.reason.contains("maximum")
    }));

    let exact_long_link = "L".repeat(MAX_TRACEABILITY_FIELD_BYTES);
    let above_long_link = "L".repeat(MAX_TRACEABILITY_FIELD_BYTES + 1);
    let exact_link_refs = [exact_long_link.as_str()];
    let above_link_refs = [above_long_link.as_str()];
    let mut exact_link_rows = requirements().to_vec();
    let exact_link_base = exact_link_rows[0];
    exact_link_rows[0] = RequirementRow {
        proof_obligations: &exact_link_refs,
        ..exact_link_base
    };
    let exact_link_audit = audit_traceability(&exact_link_rows, proof_obligations());
    assert!(
        !exact_link_audit.diagnostics.iter().any(|diagnostic| {
            diagnostic.requirement_id == "B1"
                && diagnostic.reason.contains("link is")
                && diagnostic.reason.contains("maximum")
        }),
        "the exact link byte cap must reach reference validation"
    );
    let mut above_link_rows = requirements().to_vec();
    let above_link_base = above_link_rows[0];
    above_link_rows[0] = RequirementRow {
        proof_obligations: &above_link_refs,
        ..above_link_base
    };
    let above_link_audit = audit_traceability(&above_link_rows, proof_obligations());
    assert!(above_link_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.requirement_id == "B1"
            && diagnostic.reason.contains("link is")
            && diagnostic.reason.contains("maximum")
    }));
}

#[test]
fn success_and_failure_outputs_are_input_permutation_invariant() {
    let canonical = traceability_ledger_json().unwrap();
    let mut rows = requirements().to_vec();
    rows.reverse();
    let mut obligations = proof_obligations().to_vec();
    obligations.reverse();
    assert_eq!(
        generate_traceability_ledger(&rows, &obligations).unwrap(),
        canonical
    );

    let missing_owner = RequirementRow {
        owner_artifact: "",
        ..requirements()[0]
    };
    let missing_gate = RequirementRow {
        milestone: "",
        ..requirements()[1]
    };
    let forward = audit_traceability(&[missing_owner, missing_gate], proof_obligations());
    let reverse = audit_traceability(&[missing_gate, missing_owner], proof_obligations());
    assert_eq!(forward, reverse);

    let canonical_b1 = requirements()[0];
    let invalid_duplicate_b1 = RequirementRow {
        owner_artifact: "",
        ..canonical_b1
    };
    let duplicate_forward =
        audit_traceability(&[canonical_b1, invalid_duplicate_b1], proof_obligations());
    let duplicate_reverse =
        audit_traceability(&[invalid_duplicate_b1, canonical_b1], proof_obligations());
    assert_eq!(duplicate_forward, duplicate_reverse);
}

#[test]
fn generated_json_contains_every_column_and_the_complete_po_index() {
    let json = traceability_ledger_json().unwrap();
    assert!(json.starts_with(&format!("{{\"schema\":\"{TRACEABILITY_SCHEMA}\"")));
    assert!(json.ends_with('}'));
    assert!(json.contains(&format!("\"authority\":\"{TRACEABILITY_AUTHORITY}\"")));
    assert!(json.contains("\"source_snapshot\":null"));
    assert!(!json.contains("\"status\":\"verified\""));
    assert!(!json.contains("\"status\":\"validated\""));
    for field in [
        "requirement_id",
        "capability_property",
        "blocker",
        "owner_artifact",
        "prerequisite_phase",
        "milestone",
        "flagship",
        "benchmark_data",
        "proof_obligations",
        "claim_boundary",
        "status",
    ] {
        assert!(json.contains(&format!("\"{field}\":")), "missing {field}");
    }
    assert_eq!(json.matches("\"requirement_id\":").count(), 30);
    assert_eq!(json.matches("\"summary\":").count(), 25);
    assert!(json.contains("\"id\":\"PO-1\""));
    assert!(json.contains("\"id\":\"PO-25\""));
}
