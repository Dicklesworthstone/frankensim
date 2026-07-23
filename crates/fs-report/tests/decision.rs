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
use fs_report::decision_headline_markdown;
use fs_session::{AppliedSafetyFactor, DecisionAssessment, DecisionRequirement, EvidenceRef};
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
        AppliedSafetyFactor::try_new(1.25, artifact("safety-factor-policy")).expect("valid factor"),
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
        "safety-factor-policy",
        "fs-evidence:certified-f64:v1",
        "thermal-context",
        "boundary-conditions",
        "suggested evidence `sensor-campaign`",
        "The assessment retains 1 explicit evidence recommendation(s).",
        "### Exact audit projection",
        "    decision-assessment-v1",
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
