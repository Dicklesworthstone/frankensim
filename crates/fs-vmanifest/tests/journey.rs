//! V.4.1 Journey DSL and scoped-receipt conformance.

use fs_blake3::{ContentHash, hash_domain};
use fs_vmanifest::FiveExplicits;
use fs_vmanifest::journey::{
    ArtifactSandbox, AttemptReceipt, CampaignReceipt, ClaimAdjudication, ClaimRecord,
    DomainApplicability, EpistemicGrade, EvidenceCompleteness, EvidenceIntegrity, EvidenceMethod,
    EvidenceMethodSet, ExecutionDisposition, JobReceipt, JourneyCursor, JourneyDefaults,
    JourneyManifest, JourneyPhase, LEGACY_ASSESSMENT_MIGRATION_IDENTITY_DOMAIN,
    LEGACY_ASSESSMENT_MIGRATION_SCHEMA_VERSION, LegacyAssessmentMigrationLoss, LegacyClaimStateV0,
    LegacyCombinedEvidenceV0, OperationReceipt, OperationVerb, OperationalSupport, ProcessCode,
    PromotionEffect, PublicSurfaceIdentity, ReceiptOutcome, RequestedPredicateOutcome,
    ScientificAssessment, TypedSkip, WorkloadScale, migrate_legacy_assessment_v0,
};
use fs_vmanifest::v1::{
    ClaimId, ClaimKind, ClaimRelationReceipt, ClaimRevision, JourneyId, QuantifierVariance,
    RelationKind,
};
use fs_vmanifest::v1_selection::{BuiltinProfile, CompositeProfile, ProfileId, Stratum};

fn digest(label: &str) -> ContentHash {
    hash_domain("org.frankensim.fs-vmanifest.journey-test", label.as_bytes())
}

fn manifest() -> JourneyManifest {
    JourneyManifest {
        journey: JourneyId::new("journey/conformance").expect("journey id"),
        schema_version: 1,
        explicits: FiveExplicits {
            units: "SI units; nanoseconds for operation budgets",
            seeds: "Philox stream journey/conformance, seed 42",
            budgets: "time=60s; memory=1GiB; accuracy=1e-12",
            versions: "journey-v1; rust-2024; fixture-v3",
            capabilities: "cpu; deterministic; artifact-cas",
        },
        artifact_sandbox: ArtifactSandbox {
            relative_root: "artifacts/journey-conformance".to_owned(),
            max_artifacts: 128,
            max_bytes: 1 << 30,
            retention_policy: "retain receipts permanently; payloads 30 days".to_owned(),
        },
        public_surface: PublicSurfaceIdentity {
            catalog: "fs-cli/catalog".to_owned(),
            schema_version: 1,
            content_digest: digest("catalog"),
        },
        stratum: Stratum::Core,
        profile: ProfileId::Builtin(BuiltinProfile::Standard),
        selection_digest: digest("selection"),
        scale: WorkloadScale {
            logical_cases: 32,
            shards: 4,
            max_concurrency: 2,
        },
        defaults: JourneyDefaults {
            operation_timeout_ns: 60_000_000_000,
            drain_timeout_ns: 5_000_000_000,
            max_attempts: 2,
        },
    }
}

fn revision(statement: &str, domain: &str) -> ClaimRevision {
    ClaimRevision {
        claim: ClaimId::new("claim/energy-conservation").expect("claim id"),
        kind: ClaimKind::QuantitativeBound,
        statement: statement.to_owned(),
        quantifiers: "for every admitted deterministic fixture".to_owned(),
        units_conventions: "SI; relative error is dimensionless".to_owned(),
        hypotheses: "finite inputs; deterministic mode; declared tolerance".to_owned(),
        domain: domain.to_owned(),
        surface: "fs-solid::advance/CONTRACT#energy".to_owned(),
        no_claim: "no authority outside the declared fixture domain".to_owned(),
        supersedes: None,
    }
}

fn revision_id(revision: &ClaimRevision) -> ContentHash {
    revision.revision_id().expect("test revision admits")
}

fn claim_record() -> ClaimRecord {
    let revision = revision(
        "relative energy drift is at most 1e-10",
        "bounded elastic fixture corpus",
    );
    ClaimRecord {
        subject: revision_id(&revision),
        revisions: vec![revision],
        relations: Vec::new(),
    }
}

fn science(adjudication: ClaimAdjudication) -> ScientificAssessment {
    ScientificAssessment {
        claim: claim_record(),
        adjudication,
        methods: EvidenceMethodSet::new([
            EvidenceMethod::Property,
            EvidenceMethod::IndependentOracle,
        ]),
        grade: EpistemicGrade::Corroborated,
        domain: DomainApplicability::Admitted,
        support: OperationalSupport::Supported,
        completeness: EvidenceCompleteness::Complete,
        integrity: EvidenceIntegrity::Verified,
        promotion: if adjudication == ClaimAdjudication::Supported {
            PromotionEffect::Promote
        } else {
            PromotionEffect::Block
        },
    }
}

fn job(adjudication: ClaimAdjudication) -> JobReceipt {
    JobReceipt {
        job_id: "job/energy".to_owned(),
        journey: JourneyId::new("journey/conformance").expect("journey id"),
        requested_predicate: "claim/prove".to_owned(),
        phase: JourneyPhase::Verify,
        outcome: ReceiptOutcome {
            execution: ExecutionDisposition::Completed,
            predicate: if adjudication == ClaimAdjudication::Supported {
                RequestedPredicateOutcome::Satisfied
            } else {
                RequestedPredicateOutcome::Unsatisfied
            },
        },
        attempts: vec![digest("attempt")],
        science: science(adjudication),
        skips: Vec::new(),
    }
}

#[test]
fn production_phase_graph_admits_workflow_edges_and_refuses_jumps() {
    assert_eq!(JourneyPhase::ALL.len(), 23, "phase inventory is explicit");
    let mut cursor = JourneyCursor::new(&manifest()).expect("manifest admits");
    for phase in [
        JourneyPhase::Preflight,
        JourneyPhase::Author,
        JourneyPhase::Validate,
        JourneyPhase::Estimate,
        JourneyPhase::Plan,
        JourneyPhase::Submit,
        JourneyPhase::Admit,
        JourneyPhase::Queue,
        JourneyPhase::Execute,
        JourneyPhase::Observe,
        JourneyPhase::Checkpoint,
        JourneyPhase::Pause,
        JourneyPhase::Migrate,
        JourneyPhase::Resume,
        JourneyPhase::Execute,
        JourneyPhase::Cancel,
        JourneyPhase::Inspect,
        JourneyPhase::Verify,
        JourneyPhase::Report,
        JourneyPhase::Share,
        JourneyPhase::Replay,
    ] {
        cursor.transition(phase).expect("declared edge admits");
    }
    assert_eq!(cursor.phase(), JourneyPhase::Replay);
    assert_eq!(cursor.history().len(), 22);

    let before = cursor.clone();
    let error = cursor
        .transition(JourneyPhase::Queue)
        .expect_err("replay cannot jump to queue");
    assert_eq!(error.rule(), "journey-illegal-transition");
    assert_eq!(cursor, before, "refusal is non-mutating");
}

#[test]
fn import_fork_pause_and_replay_branches_are_explicit() {
    assert!(JourneyPhase::Preflight.allows(JourneyPhase::Import));
    assert!(JourneyPhase::Checkpoint.allows(JourneyPhase::Fork));
    assert!(JourneyPhase::Pause.allows(JourneyPhase::Migrate));
    assert!(JourneyPhase::Cancel.allows(JourneyPhase::Report));
    assert!(JourneyPhase::Replay.allows(JourneyPhase::Verify));
    assert!(!JourneyPhase::Author.allows(JourneyPhase::Execute));
    assert!(!JourneyPhase::Queue.allows(JourneyPhase::Share));
}

#[test]
fn manifest_requires_five_explicits_and_a_normalized_sandbox() {
    let mut missing = manifest();
    missing.explicits.seeds = "";
    assert_eq!(
        missing.validate().expect_err("missing explicit").rule(),
        "journey-five-explicits"
    );

    let mut escaping = manifest();
    escaping.artifact_sandbox.relative_root = "../outside".to_owned();
    assert_eq!(
        escaping.validate().expect_err("path escape").rule(),
        "journey-artifact-sandbox"
    );
}

#[test]
fn profile_is_atomic_and_scale_is_an_orthogonal_digest_field() {
    let base = manifest();
    base.validate().expect("base manifest");
    let mut larger = base.clone();
    larger.scale.logical_cases *= 10;
    assert_eq!(larger.profile, base.profile, "scale cannot change profile");
    assert_ne!(
        larger.digest().expect("larger digest"),
        base.digest().expect("base digest"),
        "scale remains identity-bearing"
    );

    let mut malformed = base;
    malformed.profile = ProfileId::Composite(CompositeProfile {
        id: "release-plus-release".to_owned(),
        version: 1,
        inputs: vec![BuiltinProfile::Release, BuiltinProfile::Release],
        precedence_rule: "first wins".to_owned(),
    });
    assert_eq!(
        malformed
            .validate()
            .expect_err("repeated composite input")
            .rule(),
        "v1-profile-composition"
    );
}

#[test]
fn process_projection_covers_the_exact_operation_code_table() {
    let table = [
        (RequestedPredicateOutcome::Satisfied, ProcessCode::Satisfied),
        (
            RequestedPredicateOutcome::Unsatisfied,
            ProcessCode::Unsatisfied,
        ),
        (
            RequestedPredicateOutcome::InvalidSchemaOrAdmission,
            ProcessCode::InvalidSchemaOrAdmission,
        ),
        (
            RequestedPredicateOutcome::IndeterminateOrIncomplete,
            ProcessCode::IndeterminateOrIncomplete,
        ),
        (
            RequestedPredicateOutcome::Unsupported,
            ProcessCode::Unsupported,
        ),
        (
            RequestedPredicateOutcome::CancelledAndDrained,
            ProcessCode::CancelledAndDrained,
        ),
        (
            RequestedPredicateOutcome::TimeoutFinalized,
            ProcessCode::TimeoutFinalized,
        ),
        (
            RequestedPredicateOutcome::InfrastructureError,
            ProcessCode::InfrastructureError,
        ),
        (
            RequestedPredicateOutcome::IntegrityOrSecurityFailure,
            ProcessCode::IntegrityOrSecurityFailure,
        ),
        (
            RequestedPredicateOutcome::BudgetExhaustedFinalized,
            ProcessCode::BudgetExhaustedFinalized,
        ),
    ];
    assert_eq!(
        table.map(|(_, code)| code.value()),
        [0, 10, 11, 12, 13, 14, 15, 16, 17, 18]
    );
    for (predicate, expected) in table {
        let receipt = OperationReceipt {
            operation_id: "operation/code-table".to_owned(),
            journey: JourneyId::new("journey/conformance").expect("journey id"),
            verb: OperationVerb::Execute,
            phase: JourneyPhase::Execute,
            outcome: ReceiptOutcome {
                execution: ExecutionDisposition::Completed,
                predicate,
            },
            referenced_receipt: Some(digest("opaque-job")),
            skips: Vec::new(),
        };
        assert_eq!(receipt.process_code(), expected);
    }
}

#[test]
fn status_and_cancel_acceptance_do_not_reinterpret_the_job() {
    let mut refuted = job(ClaimAdjudication::Refuted);
    refuted.outcome = ReceiptOutcome {
        execution: ExecutionDisposition::Failed,
        predicate: RequestedPredicateOutcome::Unsatisfied,
    };
    let snapshot = refuted.clone();

    let status = OperationReceipt::for_job(
        "operation/status",
        OperationVerb::Status,
        JourneyPhase::Inspect,
        &refuted,
    )
    .expect("status receipt");
    assert_eq!(status.process_code(), ProcessCode::Satisfied);
    assert_eq!(refuted, snapshot, "content reference never mutates job");
    assert_eq!(status.referenced_receipt, Some(refuted.digest().unwrap()));

    let mut running = job(ClaimAdjudication::Pending);
    running.outcome = ReceiptOutcome {
        execution: ExecutionDisposition::Running,
        predicate: RequestedPredicateOutcome::IndeterminateOrIncomplete,
    };
    running.science.completeness = EvidenceCompleteness::None;
    running.science.integrity = EvidenceIntegrity::Unknown;
    running.science.grade = EpistemicGrade::None;
    running.science.methods = EvidenceMethodSet::default();
    running.science.promotion = PromotionEffect::Hold;
    let cancel = OperationReceipt::for_job(
        "operation/cancel-accepted",
        OperationVerb::Cancel,
        JourneyPhase::Cancel,
        &running,
    )
    .expect("cancel acceptance receipt");
    assert_eq!(cancel.process_code(), ProcessCode::Satisfied);
    assert_eq!(running.outcome.execution, ExecutionDisposition::Running);
}

#[test]
fn refutation_satisfies_adjudicate_but_not_prove() {
    let refuted = job(ClaimAdjudication::Refuted);
    let prove = OperationReceipt::for_job(
        "operation/prove",
        OperationVerb::Prove,
        JourneyPhase::Verify,
        &refuted,
    )
    .expect("prove receipt");
    let adjudicate = OperationReceipt::for_job(
        "operation/adjudicate",
        OperationVerb::Adjudicate,
        JourneyPhase::Verify,
        &refuted,
    )
    .expect("adjudication receipt");
    assert_eq!(prove.process_code(), ProcessCode::Unsatisfied);
    assert_eq!(adjudicate.process_code(), ProcessCode::Satisfied);
}

#[test]
fn restricted_surviving_claim_is_a_new_exact_revision() {
    let broad = revision(
        "relative energy drift is at most 1e-10",
        "all admitted elastic fixtures",
    );
    let restricted = revision(
        "relative energy drift is at most 1e-10",
        "admitted elastic fixtures with condition number below 1e6",
    );
    assert_ne!(revision_id(&broad), revision_id(&restricted));
    let relation = ClaimRelationReceipt {
        kind: RelationKind::Restriction,
        from: revision_id(&broad),
        to: revision_id(&restricted),
        checker: "fs-vmanifest/restriction-v1".to_owned(),
        tcb: "rustc+fs-blake3".to_owned(),
        variance: QuantifierVariance::Weakened,
        domain_note: "target is an exact proper subdomain".to_owned(),
        policy_version: 1,
    };
    let record = ClaimRecord {
        subject: revision_id(&restricted),
        revisions: vec![restricted.clone(), broad],
        relations: vec![relation],
    };
    record.normalized_graph().expect("relation graph admits");
    assert_eq!(record.subject_revision(), Some(&restricted));
}

#[test]
fn completeness_and_integrity_remain_independent_axes() {
    let mut complete_corrupt = science(ClaimAdjudication::Supported);
    complete_corrupt.integrity = EvidenceIntegrity::Failed;
    complete_corrupt.promotion = PromotionEffect::Block;
    complete_corrupt
        .validate()
        .expect("complete but corrupt is representable");

    let mut partial_verified = science(ClaimAdjudication::Pending);
    partial_verified.completeness = EvidenceCompleteness::Partial;
    partial_verified.integrity = EvidenceIntegrity::Verified;
    partial_verified.promotion = PromotionEffect::Hold;
    partial_verified
        .validate()
        .expect("partial but verified is representable");

    complete_corrupt.promotion = PromotionEffect::Promote;
    assert_eq!(
        complete_corrupt
            .validate()
            .expect_err("corrupt evidence cannot promote")
            .rule(),
        "journey-promotion-laundering"
    );
}

#[test]
fn outside_domain_and_missing_hardware_are_distinct_but_both_unsupported() {
    let mut outside = job(ClaimAdjudication::Supported);
    outside.science.domain = DomainApplicability::OutsideDomain;
    outside.science.promotion = PromotionEffect::Block;
    let outside_operation = OperationReceipt::for_job(
        "operation/outside-domain",
        OperationVerb::Prove,
        JourneyPhase::Verify,
        &outside,
    )
    .expect("outside-domain receipt");

    let mut missing_hardware = job(ClaimAdjudication::Supported);
    missing_hardware.science.support = OperationalSupport::MissingCapability;
    missing_hardware.science.promotion = PromotionEffect::Block;
    let hardware_operation = OperationReceipt::for_job(
        "operation/missing-hardware",
        OperationVerb::Prove,
        JourneyPhase::Verify,
        &missing_hardware,
    )
    .expect("missing-hardware receipt");

    assert_ne!(
        outside.science.domain, missing_hardware.science.domain,
        "domain axis remains distinct"
    );
    assert_ne!(
        outside.science.support, missing_hardware.science.support,
        "support axis remains distinct"
    );
    assert_eq!(outside_operation.process_code(), ProcessCode::Unsupported);
    assert_eq!(hardware_operation.process_code(), ProcessCode::Unsupported);
}

#[test]
fn typed_skips_are_retained_and_identity_bearing() {
    let base = OperationReceipt {
        operation_id: "operation/verify".to_owned(),
        journey: JourneyId::new("journey/conformance").expect("journey id"),
        verb: OperationVerb::Verify,
        phase: JourneyPhase::Verify,
        outcome: ReceiptOutcome {
            execution: ExecutionDisposition::Completed,
            predicate: RequestedPredicateOutcome::IndeterminateOrIncomplete,
        },
        referenced_receipt: Some(digest("job")),
        skips: Vec::new(),
    };
    let mut skipped = base.clone();
    skipped.skips.push(TypedSkip {
        id: "skip/gpu-unavailable".to_owned(),
        predicate: "host exposes capability gpu".to_owned(),
        reason: "admitted host has no GPU".to_owned(),
        owner: "campaign/operator".to_owned(),
        promotion_effect: PromotionEffect::Block,
    });
    assert_ne!(
        base.digest().expect("base digest"),
        skipped.digest().expect("skip digest")
    );
    assert_eq!(skipped.skips.len(), 1);
}

#[test]
fn receipt_scopes_are_distinct_and_content_referenced() {
    let attempt = AttemptReceipt {
        attempt_id: "attempt/one".to_owned(),
        job_id: "job/energy".to_owned(),
        journey: JourneyId::new("journey/conformance").expect("journey id"),
        requested_predicate: "claim/prove".to_owned(),
        phase: JourneyPhase::Execute,
        outcome: ReceiptOutcome {
            execution: ExecutionDisposition::Completed,
            predicate: RequestedPredicateOutcome::Satisfied,
        },
        parent_attempt: None,
        artifacts: vec![digest("artifact")],
    };
    let mut job = job(ClaimAdjudication::Supported);
    job.attempts = vec![attempt.digest().expect("attempt digest")];
    let operation = OperationReceipt::for_job(
        "operation/status",
        OperationVerb::Status,
        JourneyPhase::Inspect,
        &job,
    )
    .expect("operation receipt");
    let campaign = CampaignReceipt {
        campaign_id: "campaign/standard".to_owned(),
        journey: job.journey.clone(),
        manifest_digest: manifest().digest().expect("manifest digest"),
        stratum: Stratum::Core,
        profile: ProfileId::Builtin(BuiltinProfile::Standard),
        selection_digest: digest("selection"),
        requested_predicate: "campaign/all-jobs-satisfied".to_owned(),
        phase: JourneyPhase::Report,
        outcome: ReceiptOutcome {
            execution: ExecutionDisposition::Completed,
            predicate: RequestedPredicateOutcome::Satisfied,
        },
        jobs: vec![job.digest().expect("job digest")],
        science: job.science.clone(),
        skips: Vec::new(),
    };
    let identities = [
        operation.digest().expect("operation digest"),
        attempt.digest().expect("attempt digest"),
        job.digest().expect("job digest"),
        campaign.digest().expect("campaign digest"),
    ];
    let unique = identities
        .into_iter()
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(unique.len(), 4, "receipt scope is identity-bearing");
}

#[test]
fn pretty_and_jsonl_surfaces_share_digest_and_process_code() {
    let receipt = OperationReceipt::for_job(
        "operation/prove",
        OperationVerb::Prove,
        JourneyPhase::Verify,
        &job(ClaimAdjudication::Supported),
    )
    .expect("operation receipt");
    let digest = receipt.digest().expect("digest").to_hex();
    let pretty = receipt.render_pretty().expect("pretty");
    let jsonl = receipt.render_json_line().expect("jsonl");
    assert!(pretty.contains(&digest));
    assert!(jsonl.contains(&digest));
    assert!(pretty.contains("code=0"));
    assert!(jsonl.contains("\"process_code\":0"));
}

#[test]
fn favorable_terminal_projection_mismatch_refuses() {
    let mismatch = ReceiptOutcome {
        execution: ExecutionDisposition::TimedOutFinalized,
        predicate: RequestedPredicateOutcome::Satisfied,
    };
    assert_eq!(
        mismatch
            .validate()
            .expect_err("timeout cannot project green")
            .rule(),
        "journey-outcome-inconsistent"
    );
}

#[test]
fn legacy_exact_claim_states_migrate_without_adjudication_laundering() {
    let exact = [
        (LegacyClaimStateV0::Pending, ClaimAdjudication::Pending),
        (LegacyClaimStateV0::Supported, ClaimAdjudication::Supported),
        (LegacyClaimStateV0::Failed, ClaimAdjudication::Failed),
        (LegacyClaimStateV0::Refuted, ClaimAdjudication::Refuted),
        (LegacyClaimStateV0::Unknown, ClaimAdjudication::Unknown),
    ];
    for (legacy, adjudication) in exact {
        let migration =
            migrate_legacy_assessment_v0(legacy, LegacyCombinedEvidenceV0::CompleteAndVerified);
        assert_eq!(migration.legacy_claim_state(), legacy);
        assert_eq!(migration.adjudication(), Some(adjudication));
        assert!(!migration.is_authority_complete());
        assert!(
            migration
                .losses()
                .any(|loss| loss == LegacyAssessmentMigrationLoss::EvidenceMethodsMissing)
        );
    }
}

#[test]
fn legacy_partial_and_unsupported_states_remain_typed_ambiguities() {
    let partial = migrate_legacy_assessment_v0(
        LegacyClaimStateV0::Partial,
        LegacyCombinedEvidenceV0::Partial,
    );
    assert_eq!(partial.adjudication(), None);
    assert_eq!(partial.completeness(), Some(EvidenceCompleteness::Partial));
    assert_eq!(partial.integrity(), None);
    assert!(
        partial
            .losses()
            .any(|loss| { loss == LegacyAssessmentMigrationLoss::PartialClaimMeaningCollapsed })
    );

    let unsupported = migrate_legacy_assessment_v0(
        LegacyClaimStateV0::Unsupported,
        LegacyCombinedEvidenceV0::Verified,
    );
    assert_eq!(unsupported.adjudication(), None);
    assert_eq!(unsupported.domain(), None);
    assert_eq!(unsupported.support(), None);
    assert_eq!(unsupported.completeness(), None);
    assert_eq!(unsupported.integrity(), Some(EvidenceIntegrity::Verified));
    assert!(
        unsupported
            .losses()
            .any(|loss| { loss == LegacyAssessmentMigrationLoss::UnsupportedCauseCollapsed })
    );
    assert_ne!(partial.digest(), unsupported.digest());
}

#[test]
fn legacy_combined_evidence_recovers_only_explicit_axes() {
    let table = [
        (
            LegacyCombinedEvidenceV0::Absent,
            Some(EvidenceCompleteness::None),
            None,
        ),
        (
            LegacyCombinedEvidenceV0::Partial,
            Some(EvidenceCompleteness::Partial),
            None,
        ),
        (
            LegacyCombinedEvidenceV0::Complete,
            Some(EvidenceCompleteness::Complete),
            None,
        ),
        (
            LegacyCombinedEvidenceV0::Verified,
            None,
            Some(EvidenceIntegrity::Verified),
        ),
        (
            LegacyCombinedEvidenceV0::CompleteAndVerified,
            Some(EvidenceCompleteness::Complete),
            Some(EvidenceIntegrity::Verified),
        ),
        (
            LegacyCombinedEvidenceV0::IntegrityFailed,
            None,
            Some(EvidenceIntegrity::Failed),
        ),
    ];
    for (legacy, completeness, integrity) in table {
        let first = migrate_legacy_assessment_v0(LegacyClaimStateV0::Supported, legacy);
        let replay = migrate_legacy_assessment_v0(LegacyClaimStateV0::Supported, legacy);
        assert_eq!(first.legacy_evidence(), legacy);
        assert_eq!(first.completeness(), completeness);
        assert_eq!(first.integrity(), integrity);
        assert_eq!(first.methods(), None);
        assert_eq!(first.grade(), None);
        assert_eq!(first.promotion(), None);
        assert_eq!(first, replay);
        assert_eq!(first.digest(), replay.digest());
        assert!(!first.is_authority_complete());
    }
}

#[test]
fn legacy_migration_identity_has_fixed_field_and_loss_order() {
    let migration = migrate_legacy_assessment_v0(
        LegacyClaimStateV0::Supported,
        LegacyCombinedEvidenceV0::CompleteAndVerified,
    );
    let mut expected = Vec::new();
    expected.extend_from_slice(&LEGACY_ASSESSMENT_MIGRATION_SCHEMA_VERSION.to_be_bytes());
    expected.extend_from_slice(&[
        1, // Legacy Supported.
        4, // Legacy CompleteAndVerified.
        1, 1, // Recovered Supported adjudication.
        0, // Methods absent.
        0, // Grade absent.
        0, // Domain absent.
        0, // Operational support absent.
        1, 0, // Recovered Complete evidence.
        1, 0, // Recovered Verified integrity.
        0, // Promotion effect absent.
    ]);
    expected.extend_from_slice(&5_u32.to_be_bytes());
    expected.extend_from_slice(&[4, 5, 6, 7, 8]);
    assert_eq!(
        migration.digest(),
        hash_domain(LEGACY_ASSESSMENT_MIGRATION_IDENTITY_DOMAIN, &expected)
    );
    assert_eq!(
        migration.losses().collect::<Vec<_>>(),
        vec![
            LegacyAssessmentMigrationLoss::EvidenceMethodsMissing,
            LegacyAssessmentMigrationLoss::EpistemicGradeMissing,
            LegacyAssessmentMigrationLoss::DomainApplicabilityMissing,
            LegacyAssessmentMigrationLoss::OperationalSupportMissing,
            LegacyAssessmentMigrationLoss::PromotionEffectMissing,
        ]
    );
}
