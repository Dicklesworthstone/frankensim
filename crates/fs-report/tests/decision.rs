//! G0/G3 decision-headline projection and tri-state preservation tests.

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::{
    uncertainty::{
        EngineeringUncertaintyBudget, EngineeringUncertaintyKind, EngineeringUncertaintyTerm,
        RequirementRelation, ScalarRequirement, TermValue, UncertaintyArtifactRef,
    },
    vv::{ArtifactId, ArtifactKind, ArtifactRef},
};
use fs_package::{EvidencePackage, Provenance};
use fs_project::{
    ConsequenceClass, DecisionGate, Metadata, ProjectDecisionAuthority, RequirementDirection,
    RequirementSeverity, RequirementSource, RequirementSourceKind, SafetyFactorPolicy,
    ThermalLimit,
};
use fs_qty::QtyAny;
use fs_report::{decision_headline_markdown, project_decision_gate_markdown};
use fs_session::{AppliedSafetyFactor, DecisionAssessment, DecisionRequirement, EvidenceRef};
use fs_session::{RequirementAuthority, RequirementAuthorityKind};
use fs_voi::recommend_unknown_resolutions;

#[derive(Debug, Clone, PartialEq, Eq)]
struct MaximumTemperature;

fn digest(label: &str) -> ContentHash {
    hash_domain(
        "org.frankensim.fs-report.test.decision.v1",
        label.as_bytes(),
    )
}

fn artifact(label: &str) -> UncertaintyArtifactRef {
    UncertaintyArtifactRef::new(label, digest(label)).expect("valid artifact fixture")
}

fn budget(with_unknown: bool) -> EngineeringUncertaintyBudget {
    let terms = EngineeringUncertaintyKind::ALL
        .into_iter()
        .map(|kind| {
            let value = if with_unknown && kind == EngineeringUncertaintyKind::BoundaryConditions {
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

fn assessment(estimate: f64, with_unknown: bool) -> DecisionAssessment<MaximumTemperature> {
    let budget = budget(with_unknown);
    let scalar = ScalarRequirement::try_new(
        "junction-temperature-limit",
        "temperature:max",
        "kelvin",
        RequirementRelation::AtMost,
        100.0,
        artifact("requirement:thermal-safety"),
    )
    .expect("valid requirement");
    let requirement = DecisionRequirement::try_new(
        scalar,
        RequirementAuthority::try_new(
            RequirementAuthorityKind::Datasheet,
            "cpu-thermal-specification",
            "rev-7",
            "table-5:tj-max",
        )
        .expect("valid requirement source"),
        AppliedSafetyFactor::try_new(1.25, artifact("safety-factor-policy")).expect("valid factor"),
        RequirementAuthority::try_new(
            RequirementAuthorityKind::InternalPolicy,
            "thermal-derating-policy",
            "2026.1",
            "section-4.2",
        )
        .expect("valid factor source"),
    )
    .expect("sourced effective requirement");
    let compliance = budget
        .assess_requirement(estimate, requirement.scalar(), &[])
        .expect("valid compliance replay");
    let attribution = budget
        .attribute_requirement(estimate, requirement.scalar(), &[])
        .expect("valid attribution replay");
    let actions = recommend_unknown_resolutions(&compliance, &[]);
    let package = EvidencePackage::new(Provenance::new("decision-report-test", "Cargo.lock:test"))
        .into_verified()
        .expect("empty deny-all package is structurally valid");
    DecisionAssessment::try_assemble(
        EvidenceRef::try_new(
            "temperature:max",
            "kelvin",
            "fs-evidence:certified-f64:v1",
            digest("quantity"),
        )
        .expect("quantity evidence"),
        requirement,
        ArtifactRef::new(
            ArtifactKind::ContextOfUse,
            ArtifactId::try_new("thermal-context").expect("valid context id"),
            digest("context"),
        ),
        compliance,
        budget,
        attribution,
        actions,
        &package,
    )
    .expect("complete decision assessment")
}

fn project_authority(gate: DecisionGate) -> ProjectDecisionAuthority {
    ProjectDecisionAuthority::try_from_project_parts(
        &Metadata {
            name: "reference-cooling-v1".to_string(),
            created: "2026-07-23".to_string(),
            context_of_use: "thermal design review".to_string(),
            intended_decision: "release the cooling design".to_string(),
            decision_gate: gate,
            consequence: ConsequenceClass::Advisory,
        },
        &ThermalLimit {
            qoi: "temperature:max".to_string(),
            class: "junction".to_string(),
            region: "cpu".to_string(),
            direction: RequirementDirection::AtMost,
            limit: QtyAny::new(100.0, fs_project::spec::dims::TEMPERATURE),
            margin: QtyAny::new(10.0, fs_project::spec::dims::TEMPERATURE),
            source: RequirementSource {
                kind: RequirementSourceKind::Datasheet,
                document: "cpu-thermal-specification".to_string(),
                version: "rev-7".to_string(),
                locator: "table-5:tj-max".to_string(),
            },
            safety_factor: SafetyFactorPolicy {
                factor: 1.25,
                source: RequirementSource {
                    kind: RequirementSourceKind::InternalPolicy,
                    document: "thermal-derating-policy".to_string(),
                    version: "2026.1".to_string(),
                    locator: "section-4.2".to_string(),
                },
            },
            severity: RequirementSeverity::ReliabilityDerating,
        },
    )
    .expect("valid project decision authority")
}

#[test]
fn indeterminate_headline_retains_units_authorities_and_flip_action() {
    let decision = assessment(90.0, true);
    let markdown = decision_headline_markdown(&decision);

    for expected in [
        "**Verdict:** `indeterminate`",
        "known band `[90, 90] kelvin`",
        "junction-temperature-limit",
        "at most `100` `kelvin`",
        "Declared safety factor:** `1.25`",
        "already reflected in the effective limit",
        "requirement:thermal-safety",
        "document `cpu-thermal-specification` version `rev-7` locator `table-5:tj-max`",
        "safety-factor-policy",
        "document `thermal-derating-policy` version `2026.1` locator `section-4.2`",
        "fs-evidence:certified-f64:v1",
        "thermal-context",
        "boundary-conditions",
        "suggested evidence `sensor-campaign`",
        "The assessment retains 1 explicit evidence recommendation(s).",
        "### Exact audit projection",
        "    decision-assessment-v2",
        "Projection only:",
    ] {
        assert!(
            markdown.contains(expected),
            "missing {expected:?}:\n{markdown}"
        );
    }
    assert!(markdown.contains(&decision.content_hash().to_string()));
    assert!(markdown.contains(&decision.replay_package().to_string()));
    assert_eq!(markdown, decision_headline_markdown(&decision));
}

#[test]
fn binary_headlines_preserve_direction_and_do_not_invent_flip_actions() {
    let compliant = decision_headline_markdown(&assessment(90.0, false));
    assert!(compliant.contains("**Verdict:** `compliant` with residual margin `"));
    assert!(compliant.contains(" kelvin`"));
    assert!(compliant.contains("No admitted unknown is reported as verdict-flipping."));
    assert!(compliant.contains("next-actions=none-required-by-current-verdict"));

    let non_compliant = decision_headline_markdown(&assessment(110.0, false));
    assert!(non_compliant.contains("**Verdict:** `non-compliant` with residual shortfall `"));
    assert!(non_compliant.contains(" kelvin`"));
    assert!(non_compliant.contains("No admitted unknown is reported as verdict-flipping."));
    assert!(!non_compliant.contains("**Verdict:** `compliant`"));
}

#[test]
fn project_gate_reports_same_indeterminate_physics_differently_by_context() {
    let budget = budget(true);
    let scoping = project_authority(DecisionGate::ScopingEstimate);
    let compliance = budget
        .assess_requirement(90.0, scoping.requirement().scalar(), &[])
        .expect("valid indeterminate replay");
    let scoping_report = project_decision_gate_markdown(&scoping, &compliance);
    let signoff_report = project_decision_gate_markdown(
        &project_authority(DecisionGate::ComplianceSignoff),
        &compliance,
    );

    for report in [&scoping_report, &signoff_report] {
        assert!(report.contains("Lower-layer verdict:** `indeterminate`"));
        assert!(report.contains("document `cpu-thermal-specification` version `rev-7`"));
        assert!(report.contains("document `thermal-derating-policy` version `2026.1`"));
        assert!(report.contains("Projection only:"));
    }
    assert!(scoping_report.contains("Gate:** `scoping-estimate`"));
    assert!(scoping_report.contains("Gate outcome:** **admitted**"));
    assert!(signoff_report.contains("Gate:** `compliance-signoff`"));
    assert!(signoff_report.contains("refused: this context requires a determinate assessment"));
    assert_eq!(
        scoping_report,
        project_decision_gate_markdown(&scoping, &compliance)
    );
}
