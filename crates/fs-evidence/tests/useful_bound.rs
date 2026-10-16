//! G0/G3 laws for contextual useful bounds and requirement absorption.

use std::cell::Cell;

use fs_blake3::hash_domain;
use fs_evidence::{
    BoundInterval, BoundOutcome, ClaimClass, NoUsefulBoundCause, UsefulBoundError,
    UsefulnessCriterion,
    uncertainty::{
        AttributionVerdictState, ComplianceVerdict, EngineeringUncertaintyBudget,
        EngineeringUncertaintyKind, EngineeringUncertaintyTerm, RequirementRelation,
        ScalarRequirement, TermValue, UncertaintyArtifactRef,
    },
};

fn artifact(label: &str) -> UncertaintyArtifactRef {
    UncertaintyArtifactRef::new(
        label,
        hash_domain("fs-evidence:test:useful-bound", label.as_bytes()),
    )
    .expect("bounded test artifact")
}

fn negligible_budget() -> EngineeringUncertaintyBudget {
    let terms = EngineeringUncertaintyKind::ALL
        .into_iter()
        .map(|kind| {
            EngineeringUncertaintyTerm::try_new(
                kind,
                TermValue::negligible(format!("{} is exact in this fixture", kind.name()))
                    .expect("nonempty explanation"),
                artifact(kind.name()),
            )
            .expect("complete term")
        })
        .collect();
    EngineeringUncertaintyBudget::try_new("temperature:max", "kelvin", terms)
        .expect("complete budget")
}

fn requirement() -> ScalarRequirement {
    ScalarRequirement::try_new(
        "junction-limit",
        "temperature:max",
        "kelvin",
        RequirementRelation::AtMost,
        100.0,
        artifact("requirement"),
    )
    .expect("sourced requirement")
}

fn criterion(max_width: f64) -> UsefulnessCriterion {
    UsefulnessCriterion::try_new("junction trip decision", "kelvin", max_width)
        .expect("valid criterion")
}

#[test]
fn same_enclosure_is_contextually_useful_or_refused() {
    let interval = BoundInterval::try_new(80.0, 90.0).expect("ordered interval");
    let loose = BoundOutcome::classify(
        interval,
        criterion(11.0),
        NoUsefulBoundCause::HorizonTooLong,
        ClaimClass::LongHorizonMeanLoad,
    );
    let tight = BoundOutcome::classify(
        interval,
        criterion(5.0),
        NoUsefulBoundCause::HorizonTooLong,
        ClaimClass::LongHorizonMeanLoad,
    );

    assert!(loose.bound().is_some());
    let refusal = tight.no_useful_bound().expect("tight criterion refuses");
    assert_eq!(refusal.cause(), NoUsefulBoundCause::HorizonTooLong);
    assert_eq!(
        refusal.suggested_reformulation(),
        ClaimClass::LongHorizonMeanLoad
    );
    assert!(refusal.width_achieved() > 10.0);
}

#[test]
fn constructors_and_context_mismatches_refuse() {
    assert!(matches!(
        BoundInterval::try_new(f64::NAN, 1.0),
        Err(UsefulBoundError::NaNInterval)
    ));
    assert!(matches!(
        BoundInterval::try_new(2.0, 1.0),
        Err(UsefulBoundError::InvertedInterval { .. })
    ));
    assert!(matches!(
        BoundInterval::try_new(f64::INFINITY, f64::INFINITY),
        Err(UsefulBoundError::PointAtInfinity)
    ));
    assert!(matches!(
        UsefulnessCriterion::try_new("", "kelvin", 1.0),
        Err(UsefulBoundError::InvalidField {
            field: "decision_context",
            ..
        })
    ));
    assert!(matches!(
        UsefulnessCriterion::try_new("thermal trip", "kelvin", 0.0),
        Err(UsefulBoundError::InvalidMaximumWidth { .. })
    ));

    let left = BoundOutcome::classify(
        BoundInterval::try_new(0.0, 1.0).expect("interval"),
        criterion(2.0),
        NoUsefulBoundCause::BudgetExhausted,
        ClaimClass::RootOrEventTime,
    );
    let right = BoundOutcome::classify(
        BoundInterval::try_new(0.0, 1.0).expect("interval"),
        UsefulnessCriterion::try_new("different decision", "kelvin", 2.0).expect("criterion"),
        NoUsefulBoundCause::BudgetExhausted,
        ClaimClass::RootOrEventTime,
    );
    assert!(matches!(
        left.compose_absorbing(
            &right,
            NoUsefulBoundCause::BudgetExhausted,
            ClaimClass::RootOrEventTime,
            |_, _| BoundInterval::try_new(0.0, 1.0),
        ),
        Err(UsefulBoundError::IncompatibleCriteria)
    ));
}

#[test]
fn refusal_is_absorbing_and_never_invokes_the_combiner() {
    let refusal = BoundOutcome::classify(
        BoundInterval::try_new(0.0, 10.0).expect("interval"),
        criterion(5.0),
        NoUsefulBoundCause::LipschitzBlowup,
        ClaimClass::BroadbandSpectrum,
    );
    let useful = BoundOutcome::classify(
        BoundInterval::try_new(1.0, 2.0).expect("interval"),
        criterion(5.0),
        NoUsefulBoundCause::LipschitzBlowup,
        ClaimClass::BroadbandSpectrum,
    );
    let called = Cell::new(false);
    let composed = refusal
        .compose_absorbing(
            &useful,
            NoUsefulBoundCause::LipschitzBlowup,
            ClaimClass::BroadbandSpectrum,
            |_, _| {
                called.set(true);
                BoundInterval::try_new(0.0, 1.0)
            },
        )
        .expect("absorbing composition");

    assert!(!called.get());
    assert_eq!(composed, refusal);
}

#[test]
fn explicit_semantic_failure_cannot_be_laundered_by_narrowness() {
    let outcome = BoundOutcome::refuse(
        BoundInterval::try_new(99.0, 99.01).expect("narrow interval"),
        criterion(1.0),
        NoUsefulBoundCause::BudgetExhausted,
        ClaimClass::RootOrEventTime,
    );
    let refusal = outcome.no_useful_bound().expect("explicit refusal");
    assert!(refusal.width_achieved() < refusal.criterion().max_width());
    assert_eq!(refusal.cause(), NoUsefulBoundCause::BudgetExhausted);
}

#[test]
fn requirement_gate_forces_indeterminate_with_cause() {
    let budget = negligible_budget();
    let requirement = requirement();
    let ordinary = budget
        .assess_requirement(90.0, &requirement, &[])
        .expect("ordinary assessment");
    assert!(matches!(ordinary, ComplianceVerdict::Compliant { .. }));

    let outcome = BoundOutcome::refuse(
        BoundInterval::try_new(89.0, 91.0).expect("interval"),
        criterion(5.0),
        NoUsefulBoundCause::DomainExit,
        ClaimClass::LongHorizonMeanLoad,
    );
    let verdict = budget
        .assess_requirement_with_bound(90.0, &outcome, &requirement, &[])
        .expect("typed gate");
    let refusal = verdict
        .no_useful_bound()
        .expect("forced indeterminate retains cause");
    assert_eq!(refusal.cause(), NoUsefulBoundCause::DomainExit);
    let known_band = if let ComplianceVerdict::Indeterminate {
        known_lower,
        known_upper,
        ..
    } = &verdict
    {
        Some((*known_lower, *known_upper))
    } else {
        None
    };
    assert!(known_band.is_some());
    assert_ne!(known_band, Some((89.0, 91.0)));
    assert!(matches!(verdict, ComplianceVerdict::Indeterminate { .. }));

    let attribution = budget
        .attribute_requirement_with_bound(90.0, &outcome, &requirement, &[])
        .expect("absorbing attribution");
    assert_eq!(attribution.baseline(), &verdict);
    assert!(attribution.decision_ranked().iter().all(|entry| {
        entry.baseline_state() == AttributionVerdictState::Indeterminate
            && entry.frozen_state() == AttributionVerdictState::Indeterminate
            && entry.influence() == 0.0
    }));

    let wrong_unit = BoundOutcome::refuse(
        BoundInterval::try_new(89.0, 91.0).expect("interval"),
        UsefulnessCriterion::try_new("junction trip decision", "seconds", 5.0).expect("criterion"),
        NoUsefulBoundCause::DomainExit,
        ClaimClass::LongHorizonMeanLoad,
    );
    assert!(
        budget
            .assess_requirement_with_bound(90.0, &wrong_unit, &requirement, &[])
            .is_err(),
        "a criterion in different units must never gate this requirement"
    );
}
