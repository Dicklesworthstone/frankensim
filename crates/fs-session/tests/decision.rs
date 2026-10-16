//! G0/G3 decision-projection binding, refusal, and replay-identity tests.

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::{
    action::ActionKind,
    uncertainty::{
        EngineeringUncertaintyBudget, EngineeringUncertaintyKind, EngineeringUncertaintyTerm,
        RequirementRelation, ScalarRequirement, TermValue, UncertaintyArtifactRef,
    },
    vv::{ArtifactId, ArtifactKind, ArtifactRef},
};
use fs_package::{EvidencePackage, Provenance, VerifiedPackage};
use fs_session::{
    AppliedSafetyFactor, DecisionAssessment, DecisionAssessmentError, DecisionRequirement,
    EvidenceRef, RequirementAuthority, RequirementAuthorityKind,
};
use fs_voi::{
    ActionValue, RecommendedEvidence, UnknownResolutionCandidate, recommend_unknown_resolutions,
};

#[derive(Debug, Clone, PartialEq, Eq)]
struct MaximumTemperature;

fn digest(label: &str) -> ContentHash {
    hash_domain(
        "org.frankensim.fs-session.test.decision.v1",
        label.as_bytes(),
    )
}

fn artifact(label: &str) -> UncertaintyArtifactRef {
    UncertaintyArtifactRef::new(label, digest(label)).expect("valid artifact fixture")
}

fn budget() -> EngineeringUncertaintyBudget {
    let terms = EngineeringUncertaintyKind::ALL
        .into_iter()
        .map(|kind| {
            let value = if kind == EngineeringUncertaintyKind::BoundaryConditions {
                TermValue::unknown("fan tolerance lacks a retained population authority")
                    .expect("named unknown")
            } else {
                TermValue::negligible(format!("{} is exact in this fixture", kind.name()))
                    .expect("named negligible term")
            };
            EngineeringUncertaintyTerm::try_new(kind, value, artifact(kind.name()))
                .expect("valid term")
        })
        .collect();
    EngineeringUncertaintyBudget::try_new("temperature:max", "kelvin", terms)
        .expect("complete budget")
}

fn scalar_requirement(limit: f64) -> ScalarRequirement {
    ScalarRequirement::try_new(
        "junction-temperature-limit",
        "temperature:max",
        "kelvin",
        RequirementRelation::AtMost,
        limit,
        artifact("requirement:thermal-safety"),
    )
    .expect("valid scalar requirement")
}

fn decision_requirement(policy: &str) -> DecisionRequirement {
    DecisionRequirement::try_new(
        scalar_requirement(100.0),
        RequirementAuthority::try_new(
            RequirementAuthorityKind::Datasheet,
            "cpu-thermal-specification",
            "rev-7",
            "table-5:tj-max",
        )
        .expect("valid requirement authority"),
        AppliedSafetyFactor::try_new(1.25, artifact(policy)).expect("valid safety factor"),
        RequirementAuthority::try_new(
            RequirementAuthorityKind::InternalPolicy,
            "thermal-derating-policy",
            "2026.1",
            "section-4.2",
        )
        .expect("valid factor authority"),
    )
    .expect("sourced decision requirement")
}

fn context(label: &str, hash: ContentHash) -> ArtifactRef {
    ArtifactRef::new(
        ArtifactKind::ContextOfUse,
        ArtifactId::try_new(label).expect("valid context id"),
        hash,
    )
}

fn replay_package(version: &str) -> VerifiedPackage {
    EvidencePackage::new(Provenance::new(version, "Cargo.lock:test"))
        .into_verified()
        .expect("empty deny-all package is structurally valid")
}

fn assemble(
    quantity_hash: ContentHash,
    context_ref: ArtifactRef,
    requirement: DecisionRequirement,
    package: &VerifiedPackage,
) -> Result<DecisionAssessment<MaximumTemperature>, DecisionAssessmentError> {
    let budget = budget();
    let compliance = budget
        .assess_requirement(90.0, requirement.scalar(), &[])
        .expect("valid compliance replay");
    let attribution = budget
        .attribute_requirement(90.0, requirement.scalar(), &[])
        .expect("valid attribution replay");
    let actions = recommend_unknown_resolutions(&compliance, &[]);
    DecisionAssessment::try_assemble(
        EvidenceRef::try_new(
            "temperature:max",
            "kelvin",
            "fs-evidence:certified-f64:v1",
            quantity_hash,
        )?,
        requirement,
        context_ref,
        compliance,
        budget,
        attribution,
        actions,
        package,
    )
}

#[test]
fn offline_reassembly_is_identical_and_explain_is_complete() {
    let package = replay_package("decision-test-a");
    let first = assemble(
        digest("quantity"),
        context("thermal-context", digest("context")),
        decision_requirement("safety-factor-policy"),
        &package,
    )
    .expect("complete decision projection");
    let replayed = assemble(
        digest("quantity"),
        context("thermal-context", digest("context")),
        decision_requirement("safety-factor-policy"),
        &package,
    )
    .expect("offline projection from the same artifacts");

    assert_eq!(first, replayed);
    assert_eq!(first.content_hash(), replayed.content_hash());
    assert!(first.validate_content_hash());
    assert_eq!(first.actions().len(), 1);
    assert_eq!(first.flip_conditions().len(), 1);
    assert!(first.largest_known_budget_link().is_some());
    assert!(first.strongest_decision_link().is_some());
    assert_eq!(first.replay_package(), package.report().merkle_root());
    let explain = first.render_explain();
    for expected in [
        "quantity=temperature:max unit=kelvin",
        "requirement=junction-temperature-limit",
        "safety-factor=1.25",
        "verdict=indeterminate",
        "budget-view ranked-known=",
        "decision-view ranked=",
        "next-actions:",
        "unpriced-kind=sensor-campaign",
    ] {
        assert!(
            explain.contains(expected),
            "missing {expected:?}: {explain}"
        );
    }
}

#[test]
fn every_top_level_authority_moves_identity() {
    let package_a = replay_package("decision-test-a");
    let package_b = replay_package("decision-test-b");
    let baseline = assemble(
        digest("quantity-a"),
        context("thermal-context", digest("context-a")),
        decision_requirement("factor-policy-a"),
        &package_a,
    )
    .expect("baseline");
    let baseline_id = baseline.content_hash();

    let changed_quantity = assemble(
        digest("quantity-b"),
        context("thermal-context", digest("context-a")),
        decision_requirement("factor-policy-a"),
        &package_a,
    )
    .expect("changed quantity");
    assert_ne!(baseline_id, changed_quantity.content_hash());

    let changed_context = assemble(
        digest("quantity-a"),
        context("thermal-context", digest("context-b")),
        decision_requirement("factor-policy-a"),
        &package_a,
    )
    .expect("changed context");
    assert_ne!(baseline_id, changed_context.content_hash());

    let changed_policy = assemble(
        digest("quantity-a"),
        context("thermal-context", digest("context-a")),
        decision_requirement("factor-policy-b"),
        &package_a,
    )
    .expect("changed factor policy");
    assert_ne!(baseline_id, changed_policy.content_hash());

    let changed_package = assemble(
        digest("quantity-a"),
        context("thermal-context", digest("context-a")),
        decision_requirement("factor-policy-a"),
        &package_b,
    )
    .expect("changed replay package");
    assert_ne!(baseline_id, changed_package.content_hash());
}

#[test]
fn missing_and_cross_bound_artifacts_refuse() {
    assert_eq!(
        EvidenceRef::<MaximumTemperature>::try_new(
            "temperature:max",
            "kelvin",
            "fs-evidence:certified-f64:v1",
            ContentHash([0; 32]),
        ),
        Err(DecisionAssessmentError::MissingArtifact {
            field: "quantity.artifact"
        })
    );
    assert_eq!(
        AppliedSafetyFactor::try_new(0.99, artifact("factor-policy")),
        Err(DecisionAssessmentError::InvalidSafetyFactor)
    );

    let package = replay_package("decision-test-a");
    let missing_context = assemble(
        digest("quantity"),
        context("thermal-context", ContentHash([0; 32])),
        decision_requirement("factor-policy"),
        &package,
    );
    assert_eq!(
        missing_context,
        Err(DecisionAssessmentError::MissingArtifact {
            field: "context.hash"
        })
    );

    let wrong_context = ArtifactRef::new(
        ArtifactKind::ValidationPlan,
        ArtifactId::try_new("not-a-context").expect("valid id"),
        digest("not-a-context"),
    );
    assert_eq!(
        assemble(
            digest("quantity"),
            wrong_context,
            decision_requirement("factor-policy"),
            &package,
        ),
        Err(DecisionAssessmentError::WrongContextKind)
    );

    let budget = budget();
    let requirement = decision_requirement("factor-policy");
    let compliance = budget
        .assess_requirement(90.0, requirement.scalar(), &[])
        .expect("compliance");
    let attribution = budget
        .attribute_requirement(90.0, requirement.scalar(), &[])
        .expect("attribution");
    let error = DecisionAssessment::<MaximumTemperature>::try_assemble(
        EvidenceRef::try_new(
            "temperature:max",
            "kelvin",
            "fs-evidence:certified-f64:v1",
            digest("quantity"),
        )
        .expect("quantity"),
        requirement,
        context("thermal-context", digest("context")),
        compliance,
        budget,
        attribution,
        Vec::new(),
        &package,
    )
    .expect_err("one flipping unknown requires one explicit action");
    assert_eq!(error, DecisionAssessmentError::ActionMismatch);
}

#[test]
fn quantity_and_requirement_cross_binding_refuses() {
    let package = replay_package("decision-test-a");
    let budget = budget();
    let requirement = decision_requirement("factor-policy");
    let compliance = budget
        .assess_requirement(90.0, requirement.scalar(), &[])
        .expect("compliance");
    let attribution = budget
        .attribute_requirement(90.0, requirement.scalar(), &[])
        .expect("attribution");
    let actions = recommend_unknown_resolutions(&compliance, &[]);
    let error = DecisionAssessment::<MaximumTemperature>::try_assemble(
        EvidenceRef::try_new(
            "temperature:mean",
            "kelvin",
            "fs-evidence:certified-f64:v1",
            digest("quantity"),
        )
        .expect("quantity ref is locally valid"),
        requirement,
        context("thermal-context", digest("context")),
        compliance,
        budget,
        attribution,
        actions,
        &package,
    )
    .expect_err("a different QoI must not project as this decision");
    assert_eq!(error, DecisionAssessmentError::QuantityMismatch);
}

#[test]
fn fs_voi_priced_zero_cost_action_preserves_infinite_ratio() {
    let package = replay_package("decision-test-a");
    let budget = budget();
    let requirement = decision_requirement("factor-policy");
    let compliance = budget
        .assess_requirement(90.0, requirement.scalar(), &[])
        .expect("compliance");
    let attribution = budget
        .attribute_requirement(90.0, requirement.scalar(), &[])
        .expect("attribution");
    let actions = recommend_unknown_resolutions(
        &compliance,
        &[UnknownResolutionCandidate::new(
            EngineeringUncertaintyKind::BoundaryConditions,
            ActionKind::SensorCampaign,
            ActionValue {
                action: "reuse-retained-fan-campaign".to_string(),
                value: 2.0,
                cost: 0.0,
                value_per_cost: f64::INFINITY,
            },
        )],
    );
    let assessment = DecisionAssessment::<MaximumTemperature>::try_assemble(
        EvidenceRef::try_new(
            "temperature:max",
            "kelvin",
            "fs-evidence:certified-f64:v1",
            digest("quantity"),
        )
        .expect("quantity"),
        requirement,
        context("thermal-context", digest("context")),
        compliance,
        budget,
        attribution,
        actions,
        &package,
    )
    .expect("fs-voi's eligible zero-cost action remains admissible");

    assert!(matches!(
        &assessment.actions()[0].recommended_evidence,
        RecommendedEvidence::Priced {
            value_per_cost,
            ..
        } if value_per_cost.is_infinite()
    ));
    assert!(assessment.render_explain().contains("value-per-cost=inf"));
}

#[test]
fn priced_action_refuses_a_mismatched_value_per_cost() {
    let package = replay_package("decision-test-a");
    let budget = budget();
    let requirement = decision_requirement("factor-policy");
    let smallest_subnormal = f64::from_bits(1);
    let compliance = budget
        .assess_requirement(90.0, requirement.scalar(), &[])
        .expect("compliance");
    let attribution = budget
        .attribute_requirement(90.0, requirement.scalar(), &[])
        .expect("attribution");
    let mut actions = recommend_unknown_resolutions(
        &compliance,
        &[UnknownResolutionCandidate::new(
            EngineeringUncertaintyKind::BoundaryConditions,
            ActionKind::SensorCampaign,
            ActionValue {
                action: "measure-fan-tolerance".to_string(),
                value: smallest_subnormal + smallest_subnormal,
                cost: smallest_subnormal,
                value_per_cost: 2.0,
            },
        )],
    );
    let RecommendedEvidence::Priced { value_per_cost, .. } = &mut actions[0].recommended_evidence
    else {
        panic!("positive priced action must remain recommended");
    };
    *value_per_cost = 1.0;

    let error = DecisionAssessment::<MaximumTemperature>::try_assemble(
        EvidenceRef::try_new(
            "temperature:max",
            "kelvin",
            "fs-evidence:certified-f64:v1",
            digest("quantity"),
        )
        .expect("quantity"),
        requirement,
        context("thermal-context", digest("context")),
        compliance,
        budget,
        attribution,
        actions,
        &package,
    )
    .expect_err("a retained ranking must equal its retained value and cost");
    assert_eq!(error, DecisionAssessmentError::InvalidActionValue);
}

#[test]
fn authority_lineage_is_identity_bearing_and_trim_canonical() {
    let package = replay_package("decision-test-lineage");
    let first = assemble(
        digest("quantity"),
        context("thermal-context", digest("context")),
        decision_requirement("safety-factor-policy"),
        &package,
    )
    .expect("complete decision projection");
    let changed = DecisionRequirement::try_new(
        scalar_requirement(100.0),
        RequirementAuthority::try_new(
            RequirementAuthorityKind::Datasheet,
            "cpu-thermal-specification",
            "rev-8",
            "table-5:tj-max",
        )
        .expect("changed source version"),
        AppliedSafetyFactor::try_new(1.25, artifact("safety-factor-policy"))
            .expect("valid safety factor"),
        RequirementAuthority::try_new(
            RequirementAuthorityKind::InternalPolicy,
            "thermal-derating-policy",
            "2026.1",
            "section-4.2",
        )
        .expect("valid factor authority"),
    )
    .expect("changed sourced requirement");
    let second = assemble(
        digest("quantity"),
        context("thermal-context", digest("context")),
        changed,
        &package,
    )
    .expect("complete changed projection");

    assert_ne!(first.content_hash(), second.content_hash());
    assert_eq!(first.requirement().source().version(), "rev-7");
    assert!(first.render_explain().contains("requirement-version=rev-7"));
    assert_eq!(
        RequirementAuthority::try_new(
            RequirementAuthorityKind::Standard,
            " standard ",
            "2026",
            "4.2",
        ),
        Err(DecisionAssessmentError::InvalidField {
            field: "requirement.source.document",
        })
    );
}
