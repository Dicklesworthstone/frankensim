//! G0/G3 conformance for the ContextOfUse-to-ExperimentCampaignIR boundary.

use fs_blake3::ContentHash;
use fs_evidence::vv::{
    AcceptanceCriterion, ApplicabilityDomain, ApplicabilityPolicy, ArtifactHeader, ArtifactId,
    ContextOfUse, DeclaredBudget, QoiId, QoiSpec, SeedDeclaration, UnitId,
};
use fs_ir::campaign::*;

fn hash(label: &str) -> ContentHash {
    fs_blake3::hash_domain("org.frankensim.fs-ir.campaign-test.v1", label.as_bytes())
}

fn artifact_id(value: &str) -> ArtifactId {
    ArtifactId::try_new(value).expect("valid artifact id")
}

fn qoi_id(value: &str) -> QoiId {
    QoiId::try_new(value).expect("valid QoI id")
}

fn unit_id(value: &str) -> UnitId {
    UnitId::try_new(value).expect("valid unit id")
}

fn context() -> ContextOfUse {
    let unit = unit_id("k");
    let header = ArtifactHeader::try_new(
        artifact_id("thermal-release-context"),
        vec![unit.clone()],
        SeedDeclaration::Fixed(0x11_2026),
        DeclaredBudget::Limit(0.25),
        DeclaredBudget::Limit(60_000),
        DeclaredBudget::Limit(64 * 1_024 * 1_024),
        vec![
            ("fs-evidence".to_owned(), "1.0.0".to_owned()),
            ("fs-ir".to_owned(), "1.0.0".to_owned()),
        ],
        vec!["campaign-compile".to_owned()],
    )
    .expect("valid context header");
    ContextOfUse::try_new(
        header,
        "Decide whether the thermal response satisfies the release envelope.",
        vec![
            QoiSpec::try_new(
                qoi_id("peak-temperature"),
                "peak temperature",
                unit,
                AcceptanceCriterion::ClosedRange {
                    lo: 290.0,
                    hi: 330.0,
                },
            )
            .expect("valid QoI"),
        ],
        ApplicabilityDomain::unconstrained(),
        ApplicabilityPolicy::Refuse,
    )
    .expect("valid ContextOfUse")
}

fn claim_id(value: &str) -> CampaignClaimId {
    CampaignClaimId::try_new(value).expect("valid claim id")
}

fn channel_id(value: &str) -> MeasurementChannelId {
    MeasurementChannelId::try_new(value).expect("valid channel id")
}

fn specimen_id(value: &str) -> SpecimenId {
    SpecimenId::try_new(value).expect("valid specimen id")
}

fn run_id(value: &str) -> CampaignRunId {
    CampaignRunId::try_new(value).expect("valid run id")
}

#[allow(clippy::too_many_lines)]
fn base_draft() -> ExperimentCampaignDraft {
    let claim = claim_id("thermal-release");
    let qoi = qoi_id("peak-temperature");
    let unit = unit_id("k");
    let channel = channel_id("thermocouple-main");
    let assembly = AssemblyId::try_new("thermal-rig").expect("valid assembly");
    let factor = FactorId::try_new("ambient-temperature").expect("valid factor");
    let resource = ResourceId::try_new("daq-rig").expect("valid resource");
    let calibration_specimen = specimen_id("coupon-calibration");
    let validation_specimen = specimen_id("coupon-validation");

    ExperimentCampaignDraft {
        history: None,
        budget: CampaignBudget {
            max_runs: 8,
            max_specimens: 8,
            max_wall_time_ms: 20_000,
            max_memory_bytes: 8 * 1_024 * 1_024,
        },
        randomization: RandomizationPlan {
            seed: 0x5eed_11,
            algorithm: "philox-4x32-10/v1".to_owned(),
            blind_assignment_commitment: hash("blind-assignment"),
        },
        claims: vec![CampaignClaim {
            id: claim.clone(),
            qois: vec![qoi.clone()],
            hypothesis: "Peak temperature remains inside the declared release interval.".to_owned(),
            decision_consequence: "Release only if the preregistered validation analysis passes."
                .to_owned(),
            evidence_gaps: vec![EvidenceGap {
                id: EvidenceGapId::try_new("heldout-temperature-response").expect("valid gap id"),
                qoi: qoi.clone(),
                expected_evidence: "physical-validation-interval".to_owned(),
                description: "No held-out physical response currently closes this QoI.".to_owned(),
            }],
        }],
        dependencies: Vec::new(),
        specimens: vec![
            SpecimenSpec {
                id: validation_specimen.clone(),
                kind: "machined-coupon".to_owned(),
            },
            SpecimenSpec {
                id: calibration_specimen.clone(),
                kind: "machined-coupon".to_owned(),
            },
        ],
        assemblies: vec![AssemblySpec {
            id: assembly.clone(),
            specimens: vec![validation_specimen.clone(), calibration_specimen.clone()],
        }],
        factors: vec![FactorSpec {
            id: factor.clone(),
            unit: unit.clone(),
            levels: vec![310.0, 295.0],
        }],
        resources: vec![CampaignResource {
            id: resource.clone(),
            capabilities: vec!["sensor.read".to_owned(), "daq.capture".to_owned()],
            max_concurrent_runs: 1,
        }],
        channels: vec![MeasurementChannel {
            id: channel.clone(),
            claim: claim.clone(),
            qoi: qoi.clone(),
            unit,
            decision_consequence: "Supplies the peak-temperature acceptance statistic.".to_owned(),
        }],
        runs: vec![
            CampaignRun {
                id: run_id("run-validation"),
                specimen: validation_specimen,
                assembly: assembly.clone(),
                partition: CampaignPartition::Validation,
                claims: vec![claim.clone()],
                channels: vec![channel.clone()],
                factors: vec![FactorSetting {
                    factor: factor.clone(),
                    level: 310.0,
                }],
                randomization_slot: 7,
                blinded: true,
                resource: resource.clone(),
                wall_time_ms: 2_000,
                memory_bytes: 2 * 1_024 * 1_024,
            },
            CampaignRun {
                id: run_id("run-calibration"),
                specimen: calibration_specimen,
                assembly,
                partition: CampaignPartition::Calibration,
                claims: vec![claim.clone()],
                channels: vec![channel],
                factors: vec![FactorSetting {
                    factor,
                    level: 295.0,
                }],
                randomization_slot: 3,
                blinded: false,
                resource,
                wall_time_ms: 2_000,
                memory_bytes: 2 * 1_024 * 1_024,
            },
        ],
        analyses: vec![PreregisteredAnalysis {
            id: AnalysisId::try_new("heldout-interval-check").expect("valid analysis id"),
            claim,
            qois: vec![qoi],
            partition: CampaignPartition::Validation,
            preregistration_hash: hash("analysis-v1"),
            method: "exact-rank-interval/v1".to_owned(),
        }],
        rules: vec![
            CampaignRule {
                id: CampaignRuleId::try_new("stop-evidence-complete").expect("valid rule id"),
                kind: CampaignRuleKind::StopAfterDrain,
                predicate: "declared validation precision reached".to_owned(),
                action: "drain active acquisition and finalize a stop receipt".to_owned(),
            },
            CampaignRule {
                id: CampaignRuleId::try_new("abort-overtemperature").expect("valid rule id"),
                kind: CampaignRuleKind::AbortToSafeState,
                predicate: "temperature exceeds the hard safety envelope".to_owned(),
                action: "de-energize, drain acquisition, and finalize the safe-state receipt"
                    .to_owned(),
            },
        ],
    }
}

#[test]
fn canonical_round_trip_retains_exact_campaign_intent() {
    let admitted = ExperimentCampaignIr::compile(context(), base_draft()).expect("campaign admits");
    let bytes = admitted.canonical_bytes().to_vec();
    let decoded =
        ExperimentCampaignIr::from_canonical_bytes(&bytes).expect("canonical bytes readmit");

    assert_eq!(decoded.id(), admitted.id());
    assert_eq!(decoded.wire_hash(), admitted.wire_hash());
    assert_eq!(decoded.canonical_bytes(), bytes.as_slice());
    assert_eq!(decoded.context(), admitted.context());
    assert_eq!(decoded.claims(), admitted.claims());
    assert_eq!(decoded.runs(), admitted.runs());
    assert!(decoded.warnings().is_empty());
}

#[test]
fn caller_run_reordering_is_nonsemantic() {
    let baseline =
        ExperimentCampaignIr::compile(context(), base_draft()).expect("baseline campaign");
    let mut reordered = base_draft();
    reordered.runs.reverse();
    reordered.specimens.reverse();
    reordered.assemblies[0].specimens.reverse();
    reordered.factors[0].levels.reverse();
    reordered.resources[0].capabilities.reverse();

    let reordered =
        ExperimentCampaignIr::compile(context(), reordered).expect("reordered campaign");
    assert_eq!(reordered.id(), baseline.id());
    assert_eq!(reordered.canonical_bytes(), baseline.canonical_bytes());
}

#[test]
fn calibration_validation_specimen_leakage_refuses() {
    let mut draft = base_draft();
    let leaking_specimen = draft.runs[0].specimen.clone();
    draft.runs[1].specimen = leaking_specimen;
    let error = ExperimentCampaignIr::compile(context(), draft)
        .expect_err("exclusive partitions cannot share a specimen");
    assert!(matches!(error, CampaignError::PartitionLeakage { .. }));
    assert_eq!(error.code(), "CampaignPartitionLeakage");
}

#[test]
fn claim_without_context_acceptance_qoi_refuses() {
    let mut draft = base_draft();
    draft.claims[0].qois.push(qoi_id("undeclared-qoi"));
    let error = ExperimentCampaignIr::compile(context(), draft)
        .expect_err("claim cannot invent a QoI with no ContextOfUse acceptance");
    assert!(matches!(
        error,
        CampaignError::UnknownReference {
            field: "campaign.claim.qois",
            ..
        }
    ));
}

#[test]
fn orphan_measurement_is_flagged_without_inventing_use() {
    let mut draft = base_draft();
    draft.channels.push(MeasurementChannel {
        id: channel_id("thermocouple-orphan"),
        claim: claim_id("thermal-release"),
        qoi: qoi_id("peak-temperature"),
        unit: unit_id("k"),
        decision_consequence: "Reserved diagnostic channel with no admitted run.".to_owned(),
    });
    let admitted = ExperimentCampaignIr::compile(context(), draft).expect("orphan is diagnostic");
    assert_eq!(
        admitted.warnings(),
        &[CampaignWarning::UnusedMeasurement {
            channel: channel_id("thermocouple-orphan"),
        }]
    );
}

#[test]
fn duplicate_specimen_identity_refuses() {
    let mut draft = base_draft();
    let duplicate = draft.specimens[0].clone();
    draft.specimens.push(duplicate);
    let error =
        ExperimentCampaignIr::compile(context(), draft).expect_err("duplicate specimen refuses");
    assert!(matches!(
        error,
        CampaignError::Duplicate {
            field: "campaign.specimens",
            ..
        }
    ));
}

#[test]
fn campaign_budget_conflict_refuses_before_identity() {
    let mut draft = base_draft();
    draft.budget.max_runs = 1;
    let error =
        ExperimentCampaignIr::compile(context(), draft).expect_err("run count exceeds budget");
    assert_eq!(
        error,
        CampaignError::BudgetConflict {
            field: "campaign.budget.max-runs",
            required: 2,
            limit: 1,
        }
    );
}

#[test]
fn preregistration_mutation_moves_campaign_identity() {
    let baseline =
        ExperimentCampaignIr::compile(context(), base_draft()).expect("baseline campaign");
    let mut changed = base_draft();
    changed.analyses[0].preregistration_hash = hash("analysis-v2");
    let changed =
        ExperimentCampaignIr::compile(context(), changed).expect("changed preregistration");

    assert_ne!(changed.id(), baseline.id());
    assert_ne!(changed.canonical_bytes(), baseline.canonical_bytes());
}

#[test]
fn circular_calibration_validation_dependencies_refuse() {
    let mut draft = base_draft();
    let second_claim = claim_id("thermal-stability");
    draft.claims.push(CampaignClaim {
        id: second_claim.clone(),
        qois: vec![qoi_id("peak-temperature")],
        hypothesis: "Peak temperature is stable across the admitted run window.".to_owned(),
        decision_consequence: "A stability failure blocks release.".to_owned(),
        evidence_gaps: Vec::new(),
    });
    let second_channel = channel_id("thermocouple-stability");
    draft.channels.push(MeasurementChannel {
        id: second_channel.clone(),
        claim: second_claim.clone(),
        qoi: qoi_id("peak-temperature"),
        unit: unit_id("k"),
        decision_consequence: "Supplies the preregistered stability statistic.".to_owned(),
    });
    for run in &mut draft.runs {
        run.claims.push(second_claim.clone());
        run.channels.push(second_channel.clone());
    }
    draft.analyses.push(PreregisteredAnalysis {
        id: AnalysisId::try_new("heldout-stability-check").expect("valid analysis id"),
        claim: second_claim.clone(),
        qois: vec![qoi_id("peak-temperature")],
        partition: CampaignPartition::Validation,
        preregistration_hash: hash("stability-analysis"),
        method: "bounded-drift-check/v1".to_owned(),
    });
    draft.dependencies = vec![
        ClaimDependency {
            prerequisite: claim_id("thermal-release"),
            dependent: second_claim.clone(),
            use_kind: EvidenceUse::CalibrationInput,
        },
        ClaimDependency {
            prerequisite: second_claim,
            dependent: claim_id("thermal-release"),
            use_kind: EvidenceUse::ValidationInput,
        },
    ];

    let error = ExperimentCampaignIr::compile(context(), draft)
        .expect_err("circular evidence flow refuses");
    assert!(matches!(error, CampaignError::DependencyCycle { .. }));
}

#[test]
fn predecessor_intent_anchor_round_trips_without_equivalence_claim() {
    let mut draft = base_draft();
    draft.history = Some(CampaignHistoryAnchor {
        source_schema_version: 0,
        source_canonical_hash: hash("legacy-campaign-bytes"),
        source_intent_hash: hash("legacy-campaign-intent"),
    });
    let admitted = ExperimentCampaignIr::compile(context(), draft).expect("history anchor admits");
    let decoded = ExperimentCampaignIr::from_canonical_bytes(admitted.canonical_bytes())
        .expect("history anchor readmits");
    assert_eq!(decoded.history(), admitted.history());
    assert_eq!(decoded.id(), admitted.id());
}
