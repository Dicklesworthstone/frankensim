//! V.1.2 complete verification-inventory G0/G3 conformance.

use fs_blake3::ContentHash;
use fs_vmanifest::inventory::{
    INVENTORY_AUTHORITY_POLICY_VERSION, INVENTORY_RECONCILIATION_POLICY_VERSION,
    InventoryConflictKind, InventoryDiffKind, InventoryDraft, InventoryFactDraft, InventoryField,
    InventoryLimits, InventoryObservationDraft, InventoryRole, InventorySourceDraft,
    InventorySourceKind, ObservationAdjudication, ObservationCompleteness, ObservationExecution,
    ObservationIntegrity, ReconciliationDraft, ReconciliationKind, compile_inventory,
};
use fs_vmanifest::v1::{
    ClaimId, ClaimKind, ClaimRelationReceipt, ClaimRevision, MANIFEST_RECORD_FIELDS,
    QuantifierVariance, RelationKind, SourceAuthority, SourcePin,
};

fn digest(byte: u8) -> ContentHash {
    assert_ne!(byte, 0, "test fixture digests must be present");
    ContentHash([byte; 32])
}

fn revision(lineage: &str, statement: &str, supersedes: Option<ContentHash>) -> ClaimRevision {
    ClaimRevision {
        claim: ClaimId::new(lineage).expect("lineage id admits"),
        kind: ClaimKind::Behavioral,
        statement: statement.to_owned(),
        quantifiers: "for every admitted deterministic fixture".to_owned(),
        units_conventions: "SI; dimensionless predicates are explicit".to_owned(),
        hypotheses: "finite inputs; exact source snapshots available".to_owned(),
        domain: "the frozen miniature repository".to_owned(),
        surface: format!("fs-example::{lineage}/CONTRACT#claim"),
        no_claim: "inventory presence is not scientific proof".to_owned(),
        supersedes,
    }
}

fn revision_id(revision: &ClaimRevision) -> ContentHash {
    revision.revision_id().expect("fixture revision admits")
}

fn source(
    kind: InventorySourceKind,
    role: InventoryRole,
    locator: &str,
    snapshot: u8,
) -> InventorySourceDraft {
    source_with_identity(
        kind,
        role,
        locator,
        digest(snapshot),
        "fixture-adapter-v1",
        1,
    )
}

fn source_with_identity(
    kind: InventorySourceKind,
    role: InventoryRole,
    locator: &str,
    snapshot: ContentHash,
    adapter_version: &str,
    adapter_policy_version: u32,
) -> InventorySourceDraft {
    let source = match kind {
        InventorySourceKind::Beads => {
            InventorySourceDraft::beads(locator, snapshot, adapter_version, adapter_policy_version)
        }
        InventorySourceKind::Contract => InventorySourceDraft::contract(
            locator,
            snapshot,
            adapter_version,
            adapter_policy_version,
        ),
        InventorySourceKind::TypedRegistry => InventorySourceDraft::typed_registry(
            locator,
            snapshot,
            adapter_version,
            adapter_policy_version,
        ),
        InventorySourceKind::CodeRegistration => InventorySourceDraft::code_registration(
            locator,
            snapshot,
            adapter_version,
            adapter_policy_version,
        ),
        InventorySourceKind::TestRegistration => InventorySourceDraft::test_registration(
            locator,
            snapshot,
            adapter_version,
            adapter_policy_version,
        ),
        InventorySourceKind::VvArtifact => InventorySourceDraft::vv_artifact(
            locator,
            snapshot,
            adapter_version,
            adapter_policy_version,
        ),
        InventorySourceKind::BenchmarkRegistry => InventorySourceDraft::benchmark_registry(
            locator,
            snapshot,
            adapter_version,
            adapter_policy_version,
        ),
        InventorySourceKind::LedgerReceipt => InventorySourceDraft::ledger_receipt(
            locator,
            snapshot,
            adapter_version,
            adapter_policy_version,
        ),
        InventorySourceKind::FrozenInventory => {
            panic!("frozen replay sources come only from a sealed inventory")
        }
        _ => panic!("fixture must be updated for a new source kind"),
    };
    assert_eq!(source.role(), role);
    assert_eq!(source.kind(), kind);
    assert_eq!(source.pin().authority, kind.authority());
    source
}

fn base_sources() -> Vec<InventorySourceDraft> {
    vec![
        source(
            InventorySourceKind::Beads,
            InventoryRole::Obligation,
            ".beads/issues.jsonl",
            1,
        ),
        source(
            InventorySourceKind::Contract,
            InventoryRole::DeclaredSemantics,
            "crates/fs-example/CONTRACT.md",
            2,
        ),
        source(
            InventorySourceKind::CodeRegistration,
            InventoryRole::ExecutableRegistration,
            "crates/fs-example/src/lib.rs",
            3,
        ),
        source(
            InventorySourceKind::TestRegistration,
            InventoryRole::ValidationContext,
            "crates/fs-example/tests/conformance.rs",
            4,
        ),
        source(
            InventorySourceKind::VvArtifact,
            InventoryRole::ObservedEvidence,
            "artifacts/fs-example/vv.json",
            5,
        ),
    ]
}

fn source_index(authority: SourceAuthority) -> usize {
    match authority {
        SourceAuthority::BeadObligation => 0,
        SourceAuthority::Contract => 1,
        SourceAuthority::GeneratedArtifact => 2,
        SourceAuthority::TestSource => 3,
        SourceAuthority::FrozenSnapshot => {
            panic!("source-snapshots and claim-revision are compiler-derived")
        }
        _ => panic!("fixture must be updated for a new source authority"),
    }
}

fn values(field: InventoryField, ordinal: usize) -> Vec<String> {
    let suffix = ordinal.to_string();
    let value = match field {
        InventoryField::BeadObligation => format!("frankensim-fixture-{suffix}"),
        InventoryField::Stratum => "core".to_owned(),
        InventoryField::CampaignProfiles => "standard".to_owned(),
        InventoryField::Ambition => "S".to_owned(),
        InventoryField::PublicSurface => format!("fs-example::surface-{suffix}"),
        InventoryField::CaseIds => format!("case:fixture-{suffix}"),
        InventoryField::JourneyIds => format!("journey:fixture-{suffix}"),
        InventoryField::Ownership => "verification-team".to_owned(),
        InventoryField::FixtureIds => format!("fixture:{suffix}"),
        InventoryField::OracleIds => format!("oracle:{suffix}"),
        InventoryField::CheckerIds => format!("checker:{suffix}"),
        InventoryField::TcbOverlap => "production-and-checker-tcb-disjoint".to_owned(),
        InventoryField::ToleranceDerivation => "exact integer predicate".to_owned(),
        InventoryField::Budgets => "time=30s; memory=256MiB".to_owned(),
        InventoryField::Capabilities => "deterministic-cpu".to_owned(),
        InventoryField::EventKinds => "verification.completed".to_owned(),
        InventoryField::Retention => "retain-receipt-and-artifact".to_owned(),
        InventoryField::ReplayCommand => "cargo test -p fs-vmanifest --test inventory".to_owned(),
        InventoryField::DsrLane => "dsr quality --tool frankensim".to_owned(),
        InventoryField::ReceiptExpectations => "operation-metadata-only".to_owned(),
        InventoryField::SourceSnapshots | InventoryField::ClaimRevision => {
            panic!("derived fields do not accept facts")
        }
        _ => panic!("fixture must be updated for a new inventory field"),
    };
    vec![value]
}

fn complete_draft(revisions: Vec<ClaimRevision>) -> InventoryDraft {
    let mut facts = Vec::new();
    let mut observations = Vec::new();
    for (ordinal, revision) in revisions.iter().enumerate() {
        let id = revision_id(revision);
        for field in InventoryField::ALL {
            if field.is_derived() {
                continue;
            }
            facts.push(InventoryFactDraft {
                revision: id,
                field,
                values: values(field, ordinal),
                source_index: source_index(field.spec().authority),
            });
        }
        observations.push(InventoryObservationDraft {
            revision: id,
            observation_id: format!("observation:fixture-{ordinal}"),
            artifact_digest: digest(
                u8::try_from(20 + ordinal).expect("fixture observation digest fits u8"),
            ),
            execution: ObservationExecution::Completed,
            completeness: ObservationCompleteness::Complete,
            integrity: ObservationIntegrity::Verified,
            adjudication: ObservationAdjudication::Refuted,
            source_index: 4,
        });
    }
    InventoryDraft {
        sources: base_sources(),
        revisions,
        relations: Vec::new(),
        facts,
        observations,
        reconciliations: Vec::new(),
        authority_policy_version: INVENTORY_AUTHORITY_POLICY_VERSION,
        reconciliation_policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
    }
}

fn reverse_inputs(mut draft: InventoryDraft) -> InventoryDraft {
    let source_count = draft.sources.len();
    draft.sources.reverse();
    for fact in &mut draft.facts {
        fact.source_index = source_count - 1 - fact.source_index;
    }
    for observation in &mut draft.observations {
        observation.source_index = source_count - 1 - observation.source_index;
    }
    for reconciliation in &mut draft.reconciliations {
        let source_index = match reconciliation {
            ReconciliationDraft::Alias { source_index, .. }
            | ReconciliationDraft::Rename { source_index, .. }
            | ReconciliationDraft::Split { source_index, .. }
            | ReconciliationDraft::Merge { source_index, .. } => source_index,
        };
        *source_index = source_count - 1 - *source_index;
        match reconciliation {
            ReconciliationDraft::Split { successors, .. } => successors.reverse(),
            ReconciliationDraft::Merge { predecessors, .. } => predecessors.reverse(),
            ReconciliationDraft::Alias { .. } | ReconciliationDraft::Rename { .. } => {}
        }
    }
    draft.revisions.reverse();
    draft.relations.reverse();
    draft.facts.reverse();
    draft.observations.reverse();
    draft.reconciliations.reverse();
    draft
}

fn one_revision_draft() -> InventoryDraft {
    complete_draft(vec![revision(
        "claim/fixture",
        "the fixture behavior is deterministic",
        None,
    )])
}

fn assert_resource(
    error: &fs_vmanifest::inventory::InventoryRefusal,
    quantity: &'static str,
    required: u64,
    admitted: u64,
    unit: &'static str,
) {
    assert_eq!(error.rule(), "inventory-resource-limit");
    let resource = error
        .resource_refusal()
        .expect("resource refusal carries typed evidence");
    assert_eq!(resource.quantity, quantity);
    assert_eq!(resource.required, required);
    assert_eq!(resource.admitted, admitted);
    assert_eq!(resource.unit, unit);
}

fn expected_preconflict_semantic_bytes(draft: &InventoryDraft) -> u64 {
    let source_count = u64::try_from(draft.sources.len()).expect("source count");
    let revision_count = u64::try_from(draft.revisions.len()).expect("revision count");
    let derived_resolution_bytes = source_count
        .checked_mul(revision_count)
        .and_then(|entries| entries.checked_mul(128))
        .and_then(|bytes| bytes.checked_add(revision_count.checked_mul(64)?))
        .expect("fixture derived accounting fits");
    let source_and_fact_bytes = draft
        .sources
        .iter()
        .map(|source| source.pin().source.len() + source.adapter_version().len())
        .chain(
            draft
                .facts
                .iter()
                .flat_map(|fact| fact.values.iter().map(String::len)),
        );
    let observation_bytes = draft
        .observations
        .iter()
        .map(|observation| observation.observation_id.len());
    let reconciliation_bytes =
        draft
            .reconciliations
            .iter()
            .map(|reconciliation| match reconciliation {
                ReconciliationDraft::Alias {
                    alias,
                    canonical,
                    rationale,
                    ..
                } => rationale.len() + alias.as_str().len() + canonical.as_str().len(),
                ReconciliationDraft::Rename {
                    previous,
                    current,
                    rationale,
                    ..
                } => rationale.len() + previous.as_str().len() + current.as_str().len(),
                ReconciliationDraft::Split { rationale, .. }
                | ReconciliationDraft::Merge { rationale, .. } => rationale.len(),
            });
    let normalized_text_bytes = source_and_fact_bytes
        .chain(observation_bytes)
        .chain(reconciliation_bytes)
        .try_fold(0u64, |total, bytes| {
            total.checked_add(u64::try_from(bytes).expect("text length fits"))
        })
        .expect("fixture text accounting fits");
    derived_resolution_bytes
        .checked_add(normalized_text_bytes)
        .expect("fixture semantic accounting fits")
}

#[test]
fn g0_complete_registry_inventory_is_sealed_and_projection_complete() {
    let expected = [
        "source-snapshots",
        "bead-obligation",
        "claim-revision",
        "stratum",
        "campaign-profiles",
        "ambition",
        "public-surface",
        "case-ids",
        "journey-ids",
        "ownership",
        "fixture-ids",
        "oracle-ids",
        "checker-ids",
        "tcb-overlap",
        "tolerance-derivation",
        "budgets",
        "capabilities",
        "event-kinds",
        "retention",
        "replay-command",
        "dsr-lane",
        "receipt-expectations",
    ];
    assert_eq!(MANIFEST_RECORD_FIELDS.len(), InventoryField::ALL.len());
    for ((field, spec), expected_name) in InventoryField::ALL
        .into_iter()
        .zip(MANIFEST_RECORD_FIELDS)
        .zip(expected)
    {
        assert_eq!(field.as_str(), expected_name);
        assert_eq!(field.spec(), spec);
    }

    let inventory = compile_inventory(&one_revision_draft(), InventoryLimits::DEFAULT)
        .expect("complete inventory admits");
    assert!(inventory.conflict_free());
    assert!(!inventory.has_blocking_conflicts());
    assert_eq!(inventory.resolutions().len(), InventoryField::ALL.len());
    assert_eq!(inventory.observations().len(), 1);
    assert_eq!(
        inventory.observations()[0].adjudication(),
        ObservationAdjudication::Refuted
    );
    assert!(inventory.resolutions().iter().any(|resolution| {
        resolution.field == InventoryField::Ownership
            && resolution
                .values
                .as_ref()
                .is_some_and(|values| values == &["verification-team".to_owned()])
    }));

    let rows = inventory.semantic_rows().expect("semantic rows render");
    let human = inventory.render_human().expect("human projection");
    let json = inventory.render_json_lines().expect("json projection");
    let ledger = inventory.render_ledger_rows().expect("ledger projection");
    assert_eq!(human.lines().count(), rows.len());
    assert_eq!(json.lines().count(), rows.len());
    assert_eq!(ledger.lines().count(), rows.len());
    assert!(
        json.lines()
            .all(|line| line.starts_with("{\"schema_version\":"))
    );
    assert!(
        ledger
            .lines()
            .all(|line| line.starts_with("scope=operation outcome=inventory-metadata "))
    );
    assert!(!ledger.contains("scope=job"));
    let projection_digests = inventory
        .projection_digests()
        .expect("all projection digests derive from the sealed rows");
    assert_ne!(projection_digests.semantic, projection_digests.human);
    assert_ne!(projection_digests.semantic, projection_digests.json_lines);
    assert_ne!(projection_digests.semantic, projection_digests.ledger);
    assert_ne!(projection_digests.human, projection_digests.json_lines);
    assert_ne!(projection_digests.human, projection_digests.ledger);
    assert_ne!(projection_digests.json_lines, projection_digests.ledger);
    inventory
        .verify_replay_receipt(inventory.receipt())
        .expect("self receipt resolves exactly");
}

#[test]
fn g0_missing_required_authority_is_blocking_and_resolution_stays_explicit() {
    let mut draft = one_revision_draft();
    draft
        .facts
        .retain(|fact| fact.field != InventoryField::PublicSurface);
    let inventory = compile_inventory(&draft, InventoryLimits::DEFAULT)
        .expect("missing required authority remains inspectable");
    assert!(inventory.has_blocking_conflicts());
    assert!(inventory.conflicts().iter().any(|conflict| {
        conflict.kind == InventoryConflictKind::MissingRequiredField
            && conflict.field == Some(InventoryField::PublicSurface)
    }));
    assert_eq!(inventory.resolutions().len(), InventoryField::ALL.len());
    assert!(inventory.resolutions().iter().any(|resolution| {
        resolution.field == InventoryField::PublicSurface && resolution.values.is_none()
    }));
}

#[test]
fn g3_enumeration_order_and_tightened_limits_do_not_change_identity() {
    let draft = one_revision_draft();
    let canonical = compile_inventory(&draft, InventoryLimits::DEFAULT).expect("canonical admits");
    let reversed = compile_inventory(&reverse_inputs(draft.clone()), InventoryLimits::DEFAULT)
        .expect("reversed inputs admit");
    assert_eq!(canonical.digest(), reversed.digest());
    assert_eq!(canonical.source_set_digest(), reversed.source_set_digest());
    assert_eq!(
        canonical.projection_digests().expect("canonical digests"),
        reversed.projection_digests().expect("reversed digests")
    );
    assert_eq!(
        canonical.semantic_rows().expect("canonical rows"),
        reversed.semantic_rows().expect("reversed rows")
    );

    let row_count =
        u32::try_from(canonical.semantic_rows().expect("rows").len()).expect("row count fits u32");
    let tight = InventoryLimits {
        max_sources: u32::try_from(draft.sources.len()).expect("source count"),
        max_revisions: u32::try_from(draft.revisions.len()).expect("revision count"),
        max_relations: 0,
        max_facts: u32::try_from(draft.facts.len()).expect("fact count"),
        max_observations: u32::try_from(draft.observations.len()).expect("observation count"),
        max_reconciliations: 0,
        max_reconciliation_endpoints: 0,
        max_values_per_fact: 1,
        max_semantic_bytes: InventoryLimits::DEFAULT.max_semantic_bytes,
        max_projection_rows: row_count,
    };
    let tightened = compile_inventory(&draft, tight).expect("exact tight limits admit");
    assert_eq!(canonical.digest(), tightened.digest());
    assert_ne!(canonical.limits(), tightened.limits());
}

#[test]
fn g0_empty_and_preflight_resource_refusals_are_typed() {
    let empty = InventoryDraft {
        sources: Vec::new(),
        revisions: Vec::new(),
        relations: Vec::new(),
        facts: Vec::new(),
        observations: Vec::new(),
        reconciliations: Vec::new(),
        authority_policy_version: INVENTORY_AUTHORITY_POLICY_VERSION,
        reconciliation_policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
    };
    assert_eq!(
        compile_inventory(&empty, InventoryLimits::DEFAULT)
            .expect_err("empty inventory refuses")
            .rule(),
        "inventory-empty-sources"
    );

    let draft = one_revision_draft();
    let limits = InventoryLimits {
        max_sources: u32::try_from(draft.sources.len() - 1).expect("source cap"),
        ..InventoryLimits::DEFAULT
    };
    let error = compile_inventory(&draft, limits).expect_err("source cap+1 refuses first");
    assert_resource(
        &error,
        "source count",
        u64::try_from(draft.sources.len()).expect("required count"),
        u64::try_from(draft.sources.len() - 1).expect("admitted count"),
        "sources",
    );

    let limits = InventoryLimits {
        max_observations: 0,
        ..InventoryLimits::DEFAULT
    };
    let error = compile_inventory(&draft, limits).expect_err("observation cap+1 refuses");
    assert_resource(&error, "observation count", 1, 0, "observations");

    let required_semantic_bytes = expected_preconflict_semantic_bytes(&draft);
    let exact_semantic_limits = InventoryLimits {
        max_semantic_bytes: required_semantic_bytes,
        ..InventoryLimits::DEFAULT
    };
    compile_inventory(&draft, exact_semantic_limits).expect("exact semantic cap admits");
    let below_semantic_limits = InventoryLimits {
        max_semantic_bytes: required_semantic_bytes - 1,
        ..InventoryLimits::DEFAULT
    };
    let error = compile_inventory(&draft, below_semantic_limits)
        .expect_err("semantic exact-cap minus one refuses");
    assert_resource(
        &error,
        "aggregate inventory semantic bytes",
        required_semantic_bytes,
        required_semantic_bytes - 1,
        "bytes",
    );

    let mut conflict_draft = one_revision_draft();
    conflict_draft.facts.push(InventoryFactDraft {
        revision: revision_id(&conflict_draft.revisions[0]),
        field: InventoryField::PublicSurface,
        values: vec!["fs-example::stale-projection-row".to_owned()],
        source_index: 2,
    });
    let conflicted = compile_inventory(&conflict_draft, InventoryLimits::DEFAULT)
        .expect("retained conflict inventory admits");
    assert_eq!(conflicted.conflicts().len(), 1);
    let conflict_detail_bytes = conflicted
        .conflicts()
        .iter()
        .try_fold(0u64, |total, conflict| {
            total.checked_add(u64::try_from(conflict.detail.len()).expect("detail length fits"))
        })
        .expect("fixture conflict accounting fits");
    let conflict_semantic_bytes = expected_preconflict_semantic_bytes(&conflict_draft)
        .checked_add(conflict_detail_bytes)
        .expect("fixture conflict semantic accounting fits");
    compile_inventory(
        &conflict_draft,
        InventoryLimits {
            max_semantic_bytes: conflict_semantic_bytes,
            ..InventoryLimits::DEFAULT
        },
    )
    .expect("exact late conflict-detail semantic cap admits");
    let error = compile_inventory(
        &conflict_draft,
        InventoryLimits {
            max_semantic_bytes: conflict_semantic_bytes - 1,
            ..InventoryLimits::DEFAULT
        },
    )
    .expect_err("late conflict-detail semantic exact-cap minus one refuses");
    assert_resource(
        &error,
        "aggregate inventory semantic bytes",
        conflict_semantic_bytes,
        conflict_semantic_bytes - 1,
        "bytes",
    );
    let final_rows = u32::try_from(conflicted.semantic_rows().expect("rows").len())
        .expect("fixture rows fit u32");
    let projection_limits = InventoryLimits {
        max_projection_rows: final_rows - 1,
        ..InventoryLimits::DEFAULT
    };
    let error = compile_inventory(&conflict_draft, projection_limits)
        .expect_err("derived conflict row crosses the final projection cap");
    assert_resource(
        &error,
        "projection row count",
        u64::from(final_rows),
        u64::from(final_rows - 1),
        "rows",
    );
}

#[test]
fn g0_duplicate_equal_authority_and_cross_role_facts_are_distinct() {
    let mut duplicate = one_revision_draft();
    let revision = revision_id(&duplicate.revisions[0]);
    duplicate.facts.push(InventoryFactDraft {
        revision,
        field: InventoryField::PublicSurface,
        values: vec!["fs-example::different".to_owned()],
        source_index: 1,
    });
    assert_eq!(
        compile_inventory(&duplicate, InventoryLimits::DEFAULT)
            .expect_err("one source cannot assert one slot twice")
            .rule(),
        "inventory-duplicate-fact"
    );

    let mut equal_authority = one_revision_draft();
    equal_authority.sources.push(source(
        InventorySourceKind::Contract,
        InventoryRole::DeclaredSemantics,
        "crates/fs-example/SECOND_CONTRACT.md",
        31,
    ));
    equal_authority.facts.push(InventoryFactDraft {
        revision,
        field: InventoryField::PublicSurface,
        values: vec!["fs-example::contradiction".to_owned()],
        source_index: 5,
    });
    let inventory = compile_inventory(&equal_authority, InventoryLimits::DEFAULT)
        .expect("semantic disagreement is retained, not structurally rejected");
    assert!(inventory.has_blocking_conflicts());
    assert!(inventory.conflicts().iter().any(|conflict| {
        conflict.kind == InventoryConflictKind::EqualAuthorityDisagreement
            && conflict.field == Some(InventoryField::PublicSurface)
    }));
    assert!(inventory.resolutions().iter().any(|resolution| {
        resolution.field == InventoryField::PublicSurface && resolution.values.is_none()
    }));

    let mut corroborating_cross_role = one_revision_draft();
    corroborating_cross_role.facts.push(InventoryFactDraft {
        revision,
        field: InventoryField::PublicSurface,
        values: values(InventoryField::PublicSurface, 0),
        source_index: 2,
    });
    let inventory = compile_inventory(&corroborating_cross_role, InventoryLimits::DEFAULT)
        .expect("matching registration is retained");
    assert!(!inventory.has_blocking_conflicts());

    let mut stale_code = one_revision_draft();
    stale_code.facts.push(InventoryFactDraft {
        revision,
        field: InventoryField::PublicSurface,
        values: vec!["fs-example::stale-registration".to_owned()],
        source_index: 2,
    });
    let inventory = compile_inventory(&stale_code, InventoryLimits::DEFAULT)
        .expect("stale code is a retained semantic conflict");
    assert!(inventory.has_blocking_conflicts());
    assert!(inventory.conflicts().iter().any(|conflict| {
        conflict.kind == InventoryConflictKind::CrossRoleDisagreement
            && conflict.field == Some(InventoryField::PublicSurface)
    }));

    let mut optional_cross_role = one_revision_draft();
    optional_cross_role
        .facts
        .retain(|fact| fact.field != InventoryField::JourneyIds);
    optional_cross_role.facts.push(InventoryFactDraft {
        revision,
        field: InventoryField::JourneyIds,
        values: vec!["journey:code-only".to_owned()],
        source_index: 2,
    });
    let inventory = compile_inventory(&optional_cross_role, InventoryLimits::DEFAULT)
        .expect("optional cross-role value remains inspectable");
    assert!(inventory.has_blocking_conflicts());
    assert!(inventory.conflicts().iter().any(|conflict| {
        conflict.kind == InventoryConflictKind::CrossRoleDisagreement
            && conflict.field == Some(InventoryField::JourneyIds)
    }));

    let mut optional_empty_cross_role = one_revision_draft();
    optional_empty_cross_role
        .facts
        .retain(|fact| fact.field != InventoryField::JourneyIds);
    optional_empty_cross_role.facts.push(InventoryFactDraft {
        revision,
        field: InventoryField::JourneyIds,
        values: Vec::new(),
        source_index: 2,
    });
    let inventory = compile_inventory(&optional_empty_cross_role, InventoryLimits::DEFAULT)
        .expect("non-owning empty set agrees with authoritative optional absence");
    assert!(!inventory.has_blocking_conflicts());
}

#[test]
fn g0_source_snapshot_authority_resolves_only_the_same_locator() {
    let mut resolved = one_revision_draft();
    resolved.sources.push(source(
        InventorySourceKind::CodeRegistration,
        InventoryRole::ExecutableRegistration,
        ".beads/issues.jsonl",
        32,
    ));
    let inventory = compile_inventory(&resolved, InventoryLimits::DEFAULT)
        .expect("higher-authority same-locator pin is explicit");
    assert!(!inventory.has_blocking_conflicts());
    assert!(
        inventory
            .conflicts()
            .iter()
            .any(|conflict| { conflict.kind == InventoryConflictKind::SourceSnapshotResolved })
    );

    let mut equal = one_revision_draft();
    equal.sources.push(source(
        InventorySourceKind::Beads,
        InventoryRole::Obligation,
        ".beads/issues.jsonl",
        33,
    ));
    let inventory = compile_inventory(&equal, InventoryLimits::DEFAULT)
        .expect("equal-authority pin disagreement is retained");
    assert!(inventory.has_blocking_conflicts());
    assert!(
        inventory
            .conflicts()
            .iter()
            .any(|conflict| { conflict.kind == InventoryConflictKind::SourceSnapshotConflict })
    );

    let mut different_locator = one_revision_draft();
    different_locator.sources.push(source(
        InventorySourceKind::Beads,
        InventoryRole::Obligation,
        ".beads/archive.jsonl",
        34,
    ));
    let inventory = compile_inventory(&different_locator, InventoryLimits::DEFAULT)
        .expect("different logical sources do not conflict");
    assert!(!inventory.conflicts().iter().any(|conflict| {
        matches!(
            conflict.kind,
            InventoryConflictKind::SourceSnapshotConflict
                | InventoryConflictKind::SourceSnapshotResolved
        )
    }));
}

#[test]
fn g3_all_source_classes_and_adapter_identity_mutations_are_visible() {
    let mut draft = one_revision_draft();
    draft.sources.extend([
        source(
            InventorySourceKind::TypedRegistry,
            InventoryRole::ExecutableRegistration,
            "registries/fs-example/catalog.json",
            51,
        ),
        source(
            InventorySourceKind::BenchmarkRegistry,
            InventoryRole::ObservedEvidence,
            "artifacts/fs-example/benchmarks.json",
            52,
        ),
        source(
            InventorySourceKind::LedgerReceipt,
            InventoryRole::ObservedEvidence,
            "ledger:fs-example/receipt-set",
            53,
        ),
    ]);
    let baseline =
        compile_inventory(&draft, InventoryLimits::DEFAULT).expect("all source classes admit");
    for kind in [
        InventorySourceKind::Beads,
        InventorySourceKind::Contract,
        InventorySourceKind::TypedRegistry,
        InventorySourceKind::CodeRegistration,
        InventorySourceKind::TestRegistration,
        InventorySourceKind::VvArtifact,
        InventorySourceKind::BenchmarkRegistry,
        InventorySourceKind::LedgerReceipt,
    ] {
        assert!(
            baseline
                .sources()
                .iter()
                .any(|source| source.kind() == kind)
        );
    }

    for index in 0..draft.sources.len() {
        let original = &draft.sources[index];
        let kind = original.kind();
        let role = original.role();
        let locator = original.pin().source.clone();
        let snapshot = original.pin().snapshot;
        let policy = original.adapter_policy_version();
        let replacements = [
            source_with_identity(kind, role, &locator, snapshot, "fixture-adapter-v2", policy),
            source_with_identity(
                kind,
                role,
                &locator,
                digest(u8::try_from(100 + index).expect("source mutation digest fits u8")),
                original.adapter_version(),
                policy,
            ),
            source_with_identity(
                kind,
                role,
                &locator,
                snapshot,
                original.adapter_version(),
                policy.checked_add(1).expect("fixture policy increments"),
            ),
        ];
        for replacement in replacements {
            let mut mutated = draft.clone();
            mutated.sources[index] = replacement;
            let inventory = compile_inventory(&mutated, InventoryLimits::DEFAULT)
                .expect("adapter identity mutation admits as a new source set");
            assert_ne!(baseline.source_set_digest(), inventory.source_set_digest());
            assert_ne!(baseline.digest(), inventory.digest());
        }
    }

    let mut replay_draft = draft;
    replay_draft
        .sources
        .push(baseline.replay_source_draft("fixture-replay-adapter-v1", 1));
    let replay_inventory = compile_inventory(&replay_draft, InventoryLimits::DEFAULT)
        .expect("metadata-only frozen replay source admits");
    assert!(
        replay_inventory
            .sources()
            .iter()
            .any(|source| source.kind() == InventorySourceKind::FrozenInventory)
    );
}

#[test]
fn g0_frozen_replay_pointer_cannot_author_new_payload() {
    let prior = compile_inventory(&one_revision_draft(), InventoryLimits::DEFAULT)
        .expect("prior inventory seals replay identity");
    let replay = prior.replay_source_draft("fixture-replay-adapter-v1", 1);
    assert_eq!(replay.kind(), InventorySourceKind::FrozenInventory);
    assert_eq!(replay.role(), InventoryRole::FrozenReplay);
    assert_eq!(replay.pin().snapshot, prior.digest().content_hash());
    assert_eq!(replay.pin().source, format!("inventory:{}", prior.digest()));

    let mut forged = one_revision_draft();
    forged.sources.push(replay.clone());
    forged.facts.push(InventoryFactDraft {
        revision: revision_id(&forged.revisions[0]),
        field: InventoryField::PublicSurface,
        values: vec!["forged-frozen-override".to_owned()],
        source_index: 5,
    });
    assert_eq!(
        compile_inventory(&forged, InventoryLimits::DEFAULT)
            .expect_err("prior inventory pointer cannot mint new facts")
            .rule(),
        "inventory-frozen-replay-payload"
    );

    let mut forged_observation = one_revision_draft();
    forged_observation.sources.push(replay.clone());
    forged_observation.observations[0].source_index = 5;
    assert_eq!(
        compile_inventory(&forged_observation, InventoryLimits::DEFAULT)
            .expect_err("prior inventory pointer cannot mint observations")
            .rule(),
        "inventory-observation-source-role"
    );

    let mut forged_reconciliation = one_revision_draft();
    forged_reconciliation.sources.push(replay);
    forged_reconciliation
        .reconciliations
        .push(ReconciliationDraft::Alias {
            alias: ClaimId::new("claim/fixture").expect("fixture lineage"),
            canonical: ClaimId::new("claim/fixture").expect("fixture lineage"),
            source_index: 5,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "prior inventory tries to rewrite lineage".to_owned(),
        });
    assert_eq!(
        compile_inventory(&forged_reconciliation, InventoryLimits::DEFAULT)
            .expect_err("prior inventory pointer cannot mint reconciliation")
            .rule(),
        "inventory-reconciliation-source-role"
    );
}

#[test]
fn g0_observation_axes_and_source_roles_fail_closed() {
    let mut invalid = one_revision_draft();
    invalid.observations[0].execution = ObservationExecution::TimedOut;
    invalid.observations[0].adjudication = ObservationAdjudication::Supported;
    assert_eq!(
        compile_inventory(&invalid, InventoryLimits::DEFAULT)
            .expect_err("timeout cannot support a claim")
            .rule(),
        "inventory-observation-axes"
    );

    let mut partial = one_revision_draft();
    partial.observations[0].completeness = ObservationCompleteness::Partial;
    assert_eq!(
        compile_inventory(&partial, InventoryLimits::DEFAULT)
            .expect_err("partial evidence cannot terminally adjudicate")
            .rule(),
        "inventory-observation-axes"
    );

    let mut failed_integrity = one_revision_draft();
    failed_integrity.observations[0].integrity = ObservationIntegrity::Failed;
    assert_eq!(
        compile_inventory(&failed_integrity, InventoryLimits::DEFAULT)
            .expect_err("failed integrity cannot terminally adjudicate")
            .rule(),
        "inventory-observation-axes"
    );

    let mut nonterminal = one_revision_draft();
    nonterminal.observations[0].execution = ObservationExecution::TimedOut;
    nonterminal.observations[0].completeness = ObservationCompleteness::Partial;
    nonterminal.observations[0].integrity = ObservationIntegrity::Failed;
    nonterminal.observations[0].adjudication = ObservationAdjudication::Unknown;
    compile_inventory(&nonterminal, InventoryLimits::DEFAULT)
        .expect("unknown adjudication retains orthogonal unhealthy axes honestly");

    let mut wrong_source = one_revision_draft();
    wrong_source.observations[0].source_index = 2;
    assert_eq!(
        compile_inventory(&wrong_source, InventoryLimits::DEFAULT)
            .expect_err("code registration is not observed evidence")
            .rule(),
        "inventory-observation-source-role"
    );

    let mut unavailable = one_revision_draft();
    unavailable.observations[0].artifact_digest = ContentHash([0; 32]);
    assert_eq!(
        compile_inventory(&unavailable, InventoryLimits::DEFAULT)
            .expect_err("missing artifact digest refuses")
            .rule(),
        "inventory-observation-artifact"
    );

    let second_revision = revision(
        "claim/fixture-second",
        "the second fixture behavior is deterministic",
        None,
    );
    let mut duplicate_id = complete_draft(vec![
        revision(
            "claim/fixture-first",
            "the first fixture behavior is deterministic",
            None,
        ),
        second_revision,
    ]);
    duplicate_id.observations[1].observation_id =
        duplicate_id.observations[0].observation_id.clone();
    assert_eq!(
        compile_inventory(&duplicate_id, InventoryLimits::DEFAULT)
            .expect_err("source-local observation ids are globally unique")
            .rule(),
        "inventory-duplicate-observation"
    );

    let mut contradictory = one_revision_draft();
    contradictory.sources.push(source(
        InventorySourceKind::BenchmarkRegistry,
        InventoryRole::ObservedEvidence,
        "artifacts/fs-example/benchmark.json",
        41,
    ));
    contradictory.observations.push(InventoryObservationDraft {
        revision: revision_id(&contradictory.revisions[0]),
        observation_id: "observation:benchmark-supported".to_owned(),
        artifact_digest: digest(42),
        execution: ObservationExecution::Completed,
        completeness: ObservationCompleteness::Complete,
        integrity: ObservationIntegrity::Verified,
        adjudication: ObservationAdjudication::Supported,
        source_index: 5,
    });
    let inventory = compile_inventory(&contradictory, InventoryLimits::DEFAULT)
        .expect("contradictory observations remain inspectable");
    assert!(inventory.has_blocking_conflicts());
    assert!(inventory.conflicts().iter().any(|conflict| {
        conflict.kind == InventoryConflictKind::ObservationAdjudicationConflict
    }));
}

#[test]
fn g3_source_mutation_diff_and_historical_replay_are_exact() {
    let old_draft = one_revision_draft();
    let old = compile_inventory(&old_draft, InventoryLimits::DEFAULT).expect("old inventory");
    let old_rows = old.semantic_rows().expect("old rows");
    let old_receipt = old.receipt();

    let mut current_draft = old_draft.clone();
    current_draft.sources[1] = source(
        InventorySourceKind::Contract,
        InventoryRole::DeclaredSemantics,
        "crates/fs-example/CONTRACT.md",
        40,
    );
    let public_surface = current_draft
        .facts
        .iter_mut()
        .find(|fact| fact.field == InventoryField::PublicSurface)
        .expect("public surface fact");
    public_surface.values = vec!["fs-example::surface-mutated".to_owned()];
    let current = compile_inventory(&current_draft, InventoryLimits::DEFAULT)
        .expect("current source set compiles separately");
    assert_ne!(old.digest(), current.digest());
    assert_ne!(old.source_set_digest(), current.source_set_digest());
    assert!(current.verify_replay_receipt(old_receipt).is_err());

    let diff = old.diff(&current).expect("semantic diff");
    assert!(!diff.is_empty());
    assert_eq!(diff.from(), old.digest());
    assert_eq!(diff.to(), current.digest());
    assert!(diff.entries().iter().any(|entry| {
        entry.kind() == InventoryDiffKind::Changed
            && entry.row_kind() == "resolution"
            && entry.field() == "public-surface"
    }));
    assert!(
        diff.entries().iter().any(|entry| {
            entry.kind() == InventoryDiffKind::Added && entry.row_kind() == "source"
        })
    );
    assert!(diff.entries().iter().any(|entry| {
        entry.kind() == InventoryDiffKind::Removed && entry.row_kind() == "source"
    }));
    assert_eq!(diff.render_human().lines().count(), diff.entries().len());

    let exact_pins: Vec<SourcePin> = old_draft
        .sources
        .iter()
        .map(|source| source.pin().clone())
        .collect();
    old.verify_replay_source_availability(&exact_pins)
        .expect("exact old pins are available");
    let current_pins: Vec<SourcePin> = current_draft
        .sources
        .iter()
        .map(|source| source.pin().clone())
        .collect();
    let error = old
        .verify_replay_source_availability(&current_pins)
        .expect_err("current contract snapshot cannot replace historical snapshot");
    assert_eq!(error.rule(), "inventory-replay-source-unavailable");
    assert!(error.detail().contains("crates/fs-example/CONTRACT.md"));
    assert!(error.detail().contains(&digest(2).to_hex()));
    assert!(error.detail().contains(&digest(40).to_hex()));

    let mut caller_mutation = old_draft;
    caller_mutation.sources.clear();
    caller_mutation.facts.clear();
    assert_eq!(old.semantic_rows().expect("sealed rows remain"), old_rows);
}

#[test]
fn g3_projection_escaping_preserves_one_row_per_semantic_row() {
    let mut draft = one_revision_draft();
    let fact = draft
        .facts
        .iter_mut()
        .find(|fact| fact.field == InventoryField::PublicSurface)
        .expect("public surface fact");
    fact.values = vec!["fs-example::quoted\"\nline scope=job outcome=refuted".to_owned()];
    let inventory = compile_inventory(&draft, InventoryLimits::DEFAULT).expect("text admits");
    let rows = inventory.semantic_rows().expect("semantic rows");
    let json = inventory.render_json_lines().expect("json rows");
    let human = inventory.render_human().expect("human rows");
    let ledger = inventory.render_ledger_rows().expect("ledger rows");
    assert_eq!(json.lines().count(), rows.len());
    assert_eq!(human.lines().count(), rows.len());
    assert_eq!(ledger.lines().count(), rows.len());
    assert!(json.contains("quoted\\\"\\nline scope=job outcome=refuted"));
    assert!(!json.contains("quoted\"\nline scope=job outcome=refuted"));
    assert!(human.contains("quoted\"\\nline\\x20scope\\x3djob\\x20outcome\\x3drefuted"));
    assert!(ledger.contains("outcome=inventory-metadata"));
    assert!(ledger.contains("scope%3Djob%20outcome%3Drefuted"));
    assert!(!ledger.contains("scope=job"));
    assert!(!ledger.contains("outcome=refuted"));
}

fn topology_draft() -> (InventoryDraft, Vec<ContentHash>) {
    let root = revision("claim/root", "root historical statement", None);
    let root_id = revision_id(&root);
    let split_a = revision(
        "claim/split-a",
        "first split successor statement",
        Some(root_id),
    );
    let split_a_id = revision_id(&split_a);
    let split_b = revision(
        "claim/split-b",
        "second split successor statement",
        Some(root_id),
    );
    let split_b_id = revision_id(&split_b);
    let merged = revision(
        "claim/merged",
        "merged successor statement",
        Some(split_a_id),
    );
    let merged_id = revision_id(&merged);
    let alias = revision("claim/root-alias", "presentation alias statement", None);
    let alias_id = revision_id(&alias);
    let previous = revision(
        "claim/previous-name",
        "previous presentation statement",
        None,
    );
    let previous_id = revision_id(&previous);

    let mut draft = complete_draft(vec![root, split_a, split_b, merged, alias, previous]);
    draft
        .observations
        .retain(|observation| observation.revision == root_id);
    draft.reconciliations = vec![
        ReconciliationDraft::Merge {
            predecessors: vec![split_b_id, split_a_id],
            successor: merged_id,
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "two exact split branches converge".to_owned(),
        },
        ReconciliationDraft::Alias {
            alias: ClaimId::new("claim/root-alias").expect("alias id"),
            canonical: ClaimId::new("claim/root").expect("canonical id"),
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "historical alternate presentation".to_owned(),
        },
        ReconciliationDraft::Split {
            predecessor: root_id,
            successors: vec![split_b_id, split_a_id],
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "one claim decomposed into two exact lineages".to_owned(),
        },
        ReconciliationDraft::Rename {
            previous: ClaimId::new("claim/previous-name").expect("previous id"),
            current: ClaimId::new("claim/merged").expect("current id"),
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "explicit presentation rename".to_owned(),
        },
    ];
    (
        draft,
        vec![
            root_id,
            split_a_id,
            split_b_id,
            merged_id,
            alias_id,
            previous_id,
        ],
    )
}

#[test]
fn g3_typed_reconciliation_is_atomic_order_invariant_and_non_transferring() {
    let (draft, ids) = topology_draft();
    let inventory = compile_inventory(&draft, InventoryLimits::DEFAULT).expect("topology admits");
    assert!(!inventory.has_blocking_conflicts());
    assert_eq!(inventory.reconciliations().len(), 4);
    let kinds: Vec<ReconciliationKind> = inventory
        .reconciliations()
        .iter()
        .map(|receipt| receipt.kind())
        .collect();
    assert_eq!(
        kinds,
        vec![
            ReconciliationKind::Alias,
            ReconciliationKind::Rename,
            ReconciliationKind::Split,
            ReconciliationKind::Merge,
        ]
    );
    let split = inventory
        .reconciliations()
        .iter()
        .find_map(|receipt| receipt.split())
        .expect("split receipt");
    assert_eq!(split.0, ids[0]);
    assert_eq!(split.1, &[ids[1].min(ids[2]), ids[1].max(ids[2])]);
    let merge = inventory
        .reconciliations()
        .iter()
        .find_map(|receipt| receipt.merge())
        .expect("merge receipt");
    assert_eq!(merge.0, &[ids[1].min(ids[2]), ids[1].max(ids[2])]);
    assert_eq!(merge.1, ids[3]);
    assert!(inventory.graph().edges().is_empty());
    assert_eq!(inventory.observations().len(), 1);
    assert_eq!(inventory.observations()[0].revision(), ids[0]);
    assert_eq!(inventory.graph().representatives().len(), ids.len());

    let reversed = compile_inventory(&reverse_inputs(draft), InventoryLimits::DEFAULT)
        .expect("all enumeration and endpoint orders normalize");
    assert_eq!(inventory.digest(), reversed.digest());
    assert_eq!(
        inventory.semantic_rows().expect("canonical rows"),
        reversed.semantic_rows().expect("reversed rows")
    );
}

#[test]
fn g3_every_admitted_semantic_component_is_identity_forming() {
    let (mut draft, ids) = topology_draft();
    draft.relations.push(ClaimRelationReceipt {
        kind: RelationKind::Implication,
        from: ids[0],
        to: ids[1],
        checker: "checker:identity-baseline".to_owned(),
        tcb: "independent relation checker".to_owned(),
        variance: QuantifierVariance::Preserved,
        domain_note: "identical frozen miniature domain".to_owned(),
        policy_version: 1,
    });
    let baseline =
        compile_inventory(&draft, InventoryLimits::DEFAULT).expect("identity baseline admits");

    let mut graph_mutation = draft.clone();
    graph_mutation.relations[0].checker = "checker:identity-mutated".to_owned();
    let mut fact_mutation = draft.clone();
    fact_mutation
        .facts
        .iter_mut()
        .find(|fact| fact.revision == ids[0] && fact.field == InventoryField::PublicSurface)
        .expect("root public-surface fact")
        .values = vec!["fs-example::identity-mutated".to_owned()];
    let mut artifact_mutation = draft.clone();
    artifact_mutation.observations[0].artifact_digest = digest(70);
    let mut axis_mutation = draft.clone();
    axis_mutation.observations[0].adjudication = ObservationAdjudication::Failed;
    let mut reconciliation_mutation = draft.clone();
    match &mut reconciliation_mutation.reconciliations[0] {
        ReconciliationDraft::Alias { rationale, .. }
        | ReconciliationDraft::Rename { rationale, .. }
        | ReconciliationDraft::Split { rationale, .. }
        | ReconciliationDraft::Merge { rationale, .. } => {
            *rationale = "identity-mutated reconciliation rationale".to_owned();
        }
    }

    for (label, mutated) in [
        ("graph", graph_mutation),
        ("fact", fact_mutation),
        ("observation artifact", artifact_mutation),
        ("observation axis", axis_mutation),
        ("reconciliation", reconciliation_mutation),
    ] {
        let inventory = compile_inventory(&mutated, InventoryLimits::DEFAULT)
            .unwrap_or_else(|error| panic!("{label} mutation should admit: {error}"));
        assert_ne!(
            baseline.digest(),
            inventory.digest(),
            "{label} mutation must change inventory identity"
        );
    }

    let mut unsupported_authority = draft.clone();
    unsupported_authority.authority_policy_version = INVENTORY_AUTHORITY_POLICY_VERSION
        .checked_add(1)
        .expect("policy increment");
    assert_eq!(
        compile_inventory(&unsupported_authority, InventoryLimits::DEFAULT)
            .expect_err("unsupported authority policy refuses")
            .rule(),
        "inventory-authority-policy-version"
    );
    let mut unsupported_reconciliation = draft;
    unsupported_reconciliation.reconciliation_policy_version =
        INVENTORY_RECONCILIATION_POLICY_VERSION
            .checked_add(1)
            .expect("policy increment");
    assert_eq!(
        compile_inventory(&unsupported_reconciliation, InventoryLimits::DEFAULT)
            .expect_err("unsupported reconciliation policy refuses")
            .rule(),
        "inventory-reconciliation-policy-version"
    );
}

#[test]
fn g0_reconciliation_structure_and_authority_refuse_invalid_inputs() {
    let (mut arity, ids) = topology_draft();
    arity.reconciliations = vec![ReconciliationDraft::Split {
        predecessor: ids[0],
        successors: vec![ids[1]],
        source_index: 0,
        policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
        rationale: "invalid singleton split".to_owned(),
    }];
    assert_eq!(
        compile_inventory(&arity, InventoryLimits::DEFAULT)
            .expect_err("split requires two successors")
            .rule(),
        "inventory-reconciliation-arity"
    );

    let (mut duplicate_endpoint, ids) = topology_draft();
    duplicate_endpoint.reconciliations = vec![ReconciliationDraft::Split {
        predecessor: ids[0],
        successors: vec![ids[1], ids[1]],
        source_index: 0,
        policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
        rationale: "duplicate endpoint".to_owned(),
    }];
    assert_eq!(
        compile_inventory(&duplicate_endpoint, InventoryLimits::DEFAULT)
            .expect_err("duplicate exact endpoints refuse")
            .rule(),
        "inventory-reconciliation-duplicate-endpoint"
    );

    let split_root = revision("claim/history-split-root", "history split root", None);
    let split_root_id = revision_id(&split_root);
    let wrong_anchor = revision("claim/history-wrong-anchor", "wrong split anchor", None);
    let wrong_anchor_id = revision_id(&wrong_anchor);
    let wrongly_anchored_successor = revision(
        "claim/history-split-wrong",
        "wrongly anchored split successor",
        Some(wrong_anchor_id),
    );
    let wrongly_anchored_successor_id = revision_id(&wrongly_anchored_successor);
    let correctly_anchored_successor = revision(
        "claim/history-split-correct",
        "correctly anchored split successor",
        Some(split_root_id),
    );
    let correctly_anchored_successor_id = revision_id(&correctly_anchored_successor);
    let mut wrong_split_history = complete_draft(vec![
        split_root,
        wrong_anchor,
        wrongly_anchored_successor,
        correctly_anchored_successor,
    ]);
    wrong_split_history.reconciliations = vec![ReconciliationDraft::Split {
        predecessor: split_root_id,
        successors: vec![
            wrongly_anchored_successor_id,
            correctly_anchored_successor_id,
        ],
        source_index: 0,
        policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
        rationale: "one successor names the wrong immutable predecessor".to_owned(),
    }];
    assert_eq!(
        compile_inventory(&wrong_split_history, InventoryLimits::DEFAULT)
            .expect_err("split successor must anchor the exact predecessor")
            .rule(),
        "inventory-split-history"
    );

    let merge_a = revision("claim/history-merge-a", "history merge a", None);
    let merge_a_id = revision_id(&merge_a);
    let merge_b = revision("claim/history-merge-b", "history merge b", None);
    let merge_b_id = revision_id(&merge_b);
    let merge_outside = revision("claim/history-merge-outside", "outside merge anchor", None);
    let merge_outside_id = revision_id(&merge_outside);
    let merge_successor = revision(
        "claim/history-merge-successor",
        "merge successor with outside anchor",
        Some(merge_outside_id),
    );
    let merge_successor_id = revision_id(&merge_successor);
    let mut wrong_merge_history =
        complete_draft(vec![merge_a, merge_b, merge_outside, merge_successor]);
    wrong_merge_history.reconciliations = vec![ReconciliationDraft::Merge {
        predecessors: vec![merge_a_id, merge_b_id],
        successor: merge_successor_id,
        source_index: 0,
        policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
        rationale: "successor anchor is outside the declared predecessor set".to_owned(),
    }];
    assert_eq!(
        compile_inventory(&wrong_merge_history, InventoryLimits::DEFAULT)
            .expect_err("merge anchor must belong to the exact predecessor set")
            .rule(),
        "inventory-merge-history"
    );

    let (mut unauthorized, ids) = topology_draft();
    unauthorized.reconciliations = vec![ReconciliationDraft::Split {
        predecessor: ids[0],
        successors: vec![ids[1], ids[2]],
        source_index: 2,
        policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
        rationale: "code tries to rewrite lineage".to_owned(),
    }];
    assert_eq!(
        compile_inventory(&unauthorized, InventoryLimits::DEFAULT)
            .expect_err("code cannot author reconciliation")
            .rule(),
        "inventory-reconciliation-source-role"
    );

    let (mut duplicate_topology, ids) = topology_draft();
    duplicate_topology.reconciliations = vec![
        ReconciliationDraft::Split {
            predecessor: ids[0],
            successors: vec![ids[1], ids[2]],
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "first rationale".to_owned(),
        },
        ReconciliationDraft::Split {
            predecessor: ids[0],
            successors: vec![ids[2], ids[1]],
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "second rationale".to_owned(),
        },
    ];
    assert_eq!(
        compile_inventory(&duplicate_topology, InventoryLimits::DEFAULT)
            .expect_err("rationale cannot duplicate one topology")
            .rule(),
        "inventory-duplicate-reconciliation-topology"
    );

    let limits = InventoryLimits {
        max_reconciliation_endpoints: 3,
        ..InventoryLimits::DEFAULT
    };
    let (draft, _) = topology_draft();
    let error = compile_inventory(&draft, limits).expect_err("endpoint cap refuses preflight");
    assert_resource(&error, "reconciliation endpoint count", 10, 3, "endpoints");
}

#[test]
fn g0_reconciliation_forks_and_cycles_are_explicit() {
    let root = revision("claim/presentation-root", "root presentation", None);
    let target_a = revision("claim/presentation-a", "presentation a", None);
    let target_b = revision("claim/presentation-b", "presentation b", None);
    let mut fork = complete_draft(vec![root, target_a, target_b]);
    fork.reconciliations = vec![
        ReconciliationDraft::Alias {
            alias: ClaimId::new("claim/presentation-root").expect("root id"),
            canonical: ClaimId::new("claim/presentation-a").expect("target a"),
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "first presentation target".to_owned(),
        },
        ReconciliationDraft::Rename {
            previous: ClaimId::new("claim/presentation-root").expect("root id"),
            current: ClaimId::new("claim/presentation-b").expect("target b"),
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "second presentation target".to_owned(),
        },
    ];
    let inventory = compile_inventory(&fork, InventoryLimits::DEFAULT)
        .expect("ambiguous presentation remains inspectable");
    assert!(inventory.has_blocking_conflicts());
    assert!(
        inventory
            .conflicts()
            .iter()
            .any(|conflict| { conflict.kind == InventoryConflictKind::ReconciliationConflict })
    );

    let mut cycle = fork;
    cycle.reconciliations = vec![
        ReconciliationDraft::Alias {
            alias: ClaimId::new("claim/presentation-root").expect("root id"),
            canonical: ClaimId::new("claim/presentation-a").expect("target a"),
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "forward alias".to_owned(),
        },
        ReconciliationDraft::Rename {
            previous: ClaimId::new("claim/presentation-a").expect("target a"),
            current: ClaimId::new("claim/presentation-root").expect("root id"),
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "backward rename".to_owned(),
        },
    ];
    assert_eq!(
        compile_inventory(&cycle, InventoryLimits::DEFAULT)
            .expect_err("combined alias/rename cycles refuse")
            .rule(),
        "inventory-reconciliation-cycle"
    );

    let predecessor_a = revision("claim/merge-a", "merge predecessor a", None);
    let predecessor_a_id = revision_id(&predecessor_a);
    let predecessor_b = revision("claim/merge-b", "merge predecessor b", None);
    let predecessor_b_id = revision_id(&predecessor_b);
    let successor_a = revision(
        "claim/merge-successor-a",
        "first merge successor",
        Some(predecessor_a_id),
    );
    let successor_a_id = revision_id(&successor_a);
    let successor_b = revision(
        "claim/merge-successor-b",
        "second merge successor",
        Some(predecessor_a_id),
    );
    let successor_b_id = revision_id(&successor_b);
    let mut merge_fork =
        complete_draft(vec![predecessor_a, predecessor_b, successor_a, successor_b]);
    merge_fork.reconciliations = vec![
        ReconciliationDraft::Merge {
            predecessors: vec![predecessor_a_id, predecessor_b_id],
            successor: successor_a_id,
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "first merge target".to_owned(),
        },
        ReconciliationDraft::Merge {
            predecessors: vec![predecessor_b_id, predecessor_a_id],
            successor: successor_b_id,
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "second merge target".to_owned(),
        },
    ];
    let inventory = compile_inventory(&merge_fork, InventoryLimits::DEFAULT)
        .expect("merge target fork remains inspectable");
    assert!(inventory.has_blocking_conflicts());
    assert!(inventory.conflicts().iter().any(|conflict| {
        conflict.kind == InventoryConflictKind::ReconciliationConflict
            && conflict.detail.contains("multiple successors")
    }));

    let overlap_a = revision("claim/overlap-a", "overlap predecessor a", None);
    let overlap_a_id = revision_id(&overlap_a);
    let overlap_b = revision("claim/overlap-b", "overlap predecessor b", None);
    let overlap_b_id = revision_id(&overlap_b);
    let overlap_c = revision("claim/overlap-c", "overlap predecessor c", None);
    let overlap_c_id = revision_id(&overlap_c);
    let overlap_x = revision(
        "claim/overlap-x",
        "first partially overlapping merge successor",
        Some(overlap_a_id),
    );
    let overlap_x_id = revision_id(&overlap_x);
    let overlap_y = revision(
        "claim/overlap-y",
        "second partially overlapping merge successor",
        Some(overlap_a_id),
    );
    let overlap_y_id = revision_id(&overlap_y);
    let mut overlap = complete_draft(vec![overlap_a, overlap_b, overlap_c, overlap_x, overlap_y]);
    overlap.reconciliations = vec![
        ReconciliationDraft::Merge {
            predecessors: vec![overlap_a_id, overlap_b_id],
            successor: overlap_x_id,
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "first partially overlapping merge".to_owned(),
        },
        ReconciliationDraft::Merge {
            predecessors: vec![overlap_a_id, overlap_c_id],
            successor: overlap_y_id,
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "second partially overlapping merge".to_owned(),
        },
    ];
    let inventory = compile_inventory(&overlap, InventoryLimits::DEFAULT)
        .expect("partially overlapping merge fork remains inspectable");
    assert!(inventory.has_blocking_conflicts());
    assert!(inventory.conflicts().iter().any(|conflict| {
        conflict.kind == InventoryConflictKind::ReconciliationConflict
            && conflict
                .detail
                .contains("participates in 2 exact merge topologies")
            && conflict.detail.contains(&overlap_a_id.to_hex())
    }));
}

#[test]
fn g3_wide_merge_fork_diagnostics_scale_with_endpoint_count() {
    const WIDTH: usize = 64;
    let mut revisions = Vec::new();
    let mut predecessor_ids = Vec::new();
    for ordinal in 0..WIDTH {
        let predecessor = revision(
            &format!("claim/wide-merge-predecessor-{ordinal}"),
            &format!("wide merge predecessor {ordinal}"),
            None,
        );
        predecessor_ids.push(revision_id(&predecessor));
        revisions.push(predecessor);
    }
    let first_successor = revision(
        "claim/wide-merge-successor-a",
        "first wide merge successor",
        Some(predecessor_ids[0]),
    );
    let first_successor_id = revision_id(&first_successor);
    revisions.push(first_successor);
    let second_successor = revision(
        "claim/wide-merge-successor-b",
        "second wide merge successor",
        Some(predecessor_ids[0]),
    );
    let second_successor_id = revision_id(&second_successor);
    revisions.push(second_successor);

    let mut draft = complete_draft(revisions);
    let mut reversed_predecessors = predecessor_ids.clone();
    reversed_predecessors.reverse();
    draft.reconciliations = vec![
        ReconciliationDraft::Merge {
            predecessors: predecessor_ids,
            successor: first_successor_id,
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "first wide merge topology".to_owned(),
        },
        ReconciliationDraft::Merge {
            predecessors: reversed_predecessors,
            successor: second_successor_id,
            source_index: 0,
            policy_version: INVENTORY_RECONCILIATION_POLICY_VERSION,
            rationale: "second wide merge topology".to_owned(),
        },
    ];
    let inventory = compile_inventory(&draft, InventoryLimits::DEFAULT)
        .expect("wide merge fork remains bounded and inspectable");
    let predecessor_conflicts: Vec<_> = inventory
        .conflicts()
        .iter()
        .filter(|conflict| conflict.detail.contains("canonical merge ordinals"))
        .collect();
    assert_eq!(predecessor_conflicts.len(), WIDTH);
    let diagnostic_bytes: usize = predecessor_conflicts
        .iter()
        .map(|conflict| conflict.detail.len())
        .sum();
    assert!(
        diagnostic_bytes < WIDTH * 256,
        "per-predecessor diagnostics must not repeat complete wide topology sets"
    );
}

#[test]
fn g3_certified_equivalence_remains_v1_scientific_graph_only() {
    let left = revision("claim/equivalent-left", "left equivalent statement", None);
    let right = revision("claim/equivalent-right", "right equivalent statement", None);
    let left_id = revision_id(&left);
    let right_id = revision_id(&right);
    let mut draft = complete_draft(vec![left, right]);
    draft.relations.push(ClaimRelationReceipt {
        kind: RelationKind::CertifiedEquivalence,
        from: left_id,
        to: right_id,
        checker: "checker:equivalence".to_owned(),
        tcb: "independent symbolic checker".to_owned(),
        variance: QuantifierVariance::Preserved,
        domain_note: "identical frozen miniature domain".to_owned(),
        policy_version: 1,
    });
    let forward = compile_inventory(&draft, InventoryLimits::DEFAULT).expect("forward edge");
    assert!(forward.reconciliations().is_empty());
    assert_eq!(forward.graph().edges().len(), 1);
    assert_eq!(
        forward.graph().representative_of(&left_id),
        forward.graph().representative_of(&right_id)
    );

    let mut reverse = reverse_inputs(draft);
    reverse.relations[0].from = right_id;
    reverse.relations[0].to = left_id;
    let reverse = compile_inventory(&reverse, InventoryLimits::DEFAULT).expect("reverse edge");
    assert_eq!(forward.digest(), reverse.digest());
    assert_eq!(
        forward.semantic_rows().expect("forward rows"),
        reverse.semantic_rows().expect("reverse rows")
    );
}
