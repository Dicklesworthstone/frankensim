//! G0 governance battery for the expansion-program PR-001--PR-012 register.

use fs_govern::program_risks::{
    AssessmentStatus, MAX_OBSERVED_UNIT_PREVIEW_BYTES, ObservationDomain, ProgramRiskId,
    ProgramRiskObservation, TriggerComparison, assess_program_risks, program_risk_register_json,
    program_risks,
};

fn clear_observation(id: ProgramRiskId) -> ProgramRiskObservation<'static> {
    let risk = fs_govern::program_risks::program_risk(id);
    let value = match risk.trigger.comparison {
        TriggerComparison::GreaterThanOrEqual => risk.trigger.threshold - 1.0,
        TriggerComparison::LessThan => risk.trigger.threshold,
    };
    ProgramRiskObservation {
        id,
        value,
        unit: risk.trigger.unit,
        samples: risk.trigger.min_samples,
    }
}

#[test]
fn g0_program_register_has_exactly_twelve_ordered_unique_rows() {
    let register = program_risks();
    assert_eq!(register.len(), 12);
    for (row, id) in register.iter().zip(ProgramRiskId::ALL) {
        assert_eq!(row.id, id);
    }
    let codes: Vec<_> = register.iter().map(|row| row.id.code()).collect();
    assert_eq!(
        codes,
        vec![
            "PR-001", "PR-002", "PR-003", "PR-004", "PR-005", "PR-006", "PR-007", "PR-008",
            "PR-009", "PR-010", "PR-011", "PR-012",
        ]
    );
}

#[test]
#[allow(clippy::too_many_lines)] // The explicit tuple table is the governance mapping lock.
#[allow(clippy::float_cmp)] // Canonical threshold literals are intentionally exact artifact data.
fn g0_program_register_locks_owner_trigger_and_review_gate_mapping() {
    let expected = [
        (
            "PR-001",
            "validated-step",
            "frankensim-ext-time-validated-step-ow2o",
            ">=",
            10.0,
            "ratio-to-estimated-baseline",
            "non-negative-real",
            5,
            "E2",
        ),
        (
            "PR-002",
            "contact-detection",
            "frankensim-ext-contact-detection-ccd-tqag",
            ">=",
            100.0,
            "candidate-certificates-per-accepted-pair",
            "non-negative-real",
            10,
            "E2",
        ),
        (
            "PR-003",
            "hx-preconditioner",
            "frankensim-ext-solver-hx-preconditioner-12l2",
            ">=",
            200.0,
            "iterations",
            "non-negative-integer",
            12,
            "E4",
        ),
        (
            "PR-004",
            "chemistry-ladder",
            "frankensim-ext-gas-chemistry-ladder-paqh",
            ">=",
            1_000.0,
            "chemistry-substeps-per-flow-step",
            "non-negative-real",
            100,
            "E5",
        ),
        (
            "PR-005",
            "material-dataset",
            "frankensim-ext-matdb-seed-dataset-1sxe",
            "<",
            0.70,
            "fraction",
            "fraction-0-to-1",
            10,
            "E0c",
        ),
        (
            "PR-006",
            "adjoint-composition",
            "frankensim-ext-adjoint-composition-easb",
            "<",
            0.80,
            "fraction",
            "fraction-0-to-1",
            10,
            "E2",
        ),
        (
            "PR-007",
            "manifest-fixture",
            "frankensim-ext-manifest-fixture-r56j",
            ">=",
            1.0,
            "cycles",
            "non-negative-integer",
            1,
            "E0b",
        ),
        (
            "PR-008",
            "ledger-migration",
            "frankensim-ext-ledger-package-migration-h61n",
            ">=",
            1.0,
            "mismatches",
            "non-negative-integer",
            1,
            "E0a",
        ),
        (
            "PR-009",
            "theorem-foundry",
            "frankensim-ext-theorem-foundry-infra-zxob",
            ">=",
            1.0,
            "disagreements",
            "non-negative-integer",
            1,
            "E0d",
        ),
        (
            "PR-010",
            "scale-qualification",
            "frankensim-ext-scale-qualification-0h2j",
            ">=",
            2.0,
            "flagship-decks",
            "non-negative-integer",
            2,
            "E6",
        ),
        (
            "PR-011",
            "workflow-interop",
            "frankensim-ext-workflow-interop-lz8f",
            "<",
            0.25,
            "fraction",
            "fraction-0-to-1",
            8,
            "E7",
        ),
        (
            "PR-012",
            "safety-assurance",
            "frankensim-ext-safety-emc-assurance-te0w",
            ">=",
            1.0,
            "exported-reports",
            "non-negative-integer",
            1,
            "E7",
        ),
    ];

    for (risk, expected) in program_risks().iter().zip(expected) {
        assert_eq!(
            (
                risk.id.code(),
                risk.owner.role,
                risk.owner.bead_id,
                risk.trigger.comparison.symbol(),
                risk.trigger.threshold,
                risk.trigger.unit,
                risk.trigger.domain.code(),
                risk.trigger.min_samples,
                risk.review_gate.code(),
            ),
            expected,
            "{} governance mapping",
            risk.id.code(),
        );
    }
}

#[test]
fn g0_program_rows_are_schema_complete_and_quantitative() {
    for row in program_risks() {
        assert!(!row.name.trim().is_empty(), "{} name", row.id.code());
        assert!(!row.owner.role.trim().is_empty(), "{} role", row.id.code());
        assert!(
            row.owner.bead_id.starts_with("frankensim-"),
            "{} owner {}",
            row.id.code(),
            row.owner.bead_id
        );
        assert!(
            (1..=5).contains(&row.likelihood.score()),
            "{} likelihood",
            row.id.code()
        );
        assert!(
            (1..=5).contains(&row.impact.score()),
            "{} impact",
            row.id.code()
        );
        assert!(!row.leading_indicator.trim().is_empty());
        assert!(row.trigger.threshold.is_finite());
        assert!(!row.trigger.unit.trim().is_empty());
        assert!(!row.trigger.domain.code().trim().is_empty());
        assert!(row.trigger.min_samples > 0);
        assert!(!row.mitigation.trim().is_empty());
        assert!(!row.contingency.trim().is_empty());
        assert!((1..=5).contains(&row.residual_likelihood.score()));
        assert!((1..=5).contains(&row.residual_impact.score()));
        assert!(!row.review_gate.code().trim().is_empty());
    }
}

#[test]
fn g0_trigger_boundaries_follow_the_declared_comparator() {
    for risk in program_risks() {
        let threshold = risk.trigger.threshold;
        // The one-below probe must stay INSIDE the declared observation
        // domain: an epsilon nudge off an integer threshold is a
        // non-integral value that correctly assesses OutOfRange, not a
        // boundary sample (PR-003 exposed this).
        let just_below = match risk.trigger.domain {
            ObservationDomain::NonNegativeInteger => threshold - 1.0,
            _ => threshold - f64::EPSILON * threshold.max(1.0),
        };
        match risk.trigger.comparison {
            TriggerComparison::GreaterThanOrEqual => {
                assert_eq!(
                    risk.trigger
                        .assess(threshold, risk.trigger.min_samples, risk.trigger.unit),
                    AssessmentStatus::Triggered,
                    "{} boundary",
                    risk.id.code()
                );
                assert_eq!(
                    risk.trigger
                        .assess(just_below, risk.trigger.min_samples, risk.trigger.unit,),
                    AssessmentStatus::Clear,
                    "{} below boundary",
                    risk.id.code()
                );
            }
            TriggerComparison::LessThan => {
                assert_eq!(
                    risk.trigger
                        .assess(threshold, risk.trigger.min_samples, risk.trigger.unit),
                    AssessmentStatus::Clear,
                    "{} boundary",
                    risk.id.code()
                );
                assert_eq!(
                    risk.trigger
                        .assess(just_below, risk.trigger.min_samples, risk.trigger.unit,),
                    AssessmentStatus::Triggered,
                    "{} below boundary",
                    risk.id.code()
                );
            }
        }
    }
}

#[test]
fn g0_assessment_fails_closed_on_missing_duplicate_nonfinite_and_undersampled() {
    let observations = [
        ProgramRiskObservation {
            id: ProgramRiskId::Pr001,
            value: 1.0,
            unit: "ratio-to-estimated-baseline",
            samples: 5,
        },
        ProgramRiskObservation {
            id: ProgramRiskId::Pr001,
            value: 2.0,
            unit: "ratio-to-estimated-baseline",
            samples: 5,
        },
        ProgramRiskObservation {
            id: ProgramRiskId::Pr002,
            value: f64::NAN,
            unit: "candidate-certificates-per-accepted-pair",
            samples: 10,
        },
        ProgramRiskObservation {
            id: ProgramRiskId::Pr003,
            value: 20.0,
            unit: "iterations",
            samples: 11,
        },
        ProgramRiskObservation {
            id: ProgramRiskId::Pr004,
            value: 20.0,
            unit: "iterations",
            samples: 100,
        },
        ProgramRiskObservation {
            id: ProgramRiskId::Pr005,
            value: 1.5,
            unit: "fraction",
            samples: 10,
        },
    ];
    let assessment = assess_program_risks(&observations);
    assert_eq!(assessment.rows().len(), 12);
    assert_eq!(assessment.rows()[0].status, AssessmentStatus::Duplicate);
    assert_eq!(assessment.rows()[1].status, AssessmentStatus::NonFinite);
    assert_eq!(assessment.rows()[2].status, AssessmentStatus::UnderSampled);
    assert_eq!(assessment.rows()[3].status, AssessmentStatus::UnitMismatch);
    assert_eq!(assessment.rows()[4].status, AssessmentStatus::OutOfRange);
    assert_eq!(assessment.rows()[5].status, AssessmentStatus::Missing);
    assert_eq!(assessment.alert_count(), 12);
    assert!(!assessment.all_clear());
    let json = assessment.to_json();
    assert!(!json.contains("NaN"));
    assert!(json.contains("\"status\":\"duplicate\""));
    assert!(json.contains("\"status\":\"non-finite\""));
    assert!(json.contains("\"status\":\"under-sampled\""));
    assert!(json.contains("\"status\":\"unit-mismatch\""));
    assert!(json.contains("\"observed_unit\":\"iterations\""));
    assert!(json.contains("\"status\":\"out-of-range\""));
}

#[test]
fn g0_typed_numeric_domains_reject_negative_and_fractional_counts() {
    let cycle = fs_govern::program_risks::program_risk(ProgramRiskId::Pr007);
    assert_eq!(
        cycle
            .trigger
            .assess(-1.0, cycle.trigger.min_samples, cycle.trigger.unit),
        AssessmentStatus::OutOfRange
    );
    assert_eq!(
        cycle
            .trigger
            .assess(0.5, cycle.trigger.min_samples, cycle.trigger.unit),
        AssessmentStatus::OutOfRange
    );
    let fraction = fs_govern::program_risks::program_risk(ProgramRiskId::Pr011);
    assert_eq!(
        fraction
            .trigger
            .assess(1.01, fraction.trigger.min_samples, fraction.trigger.unit),
        AssessmentStatus::OutOfRange
    );
}

#[test]
fn g0_caller_units_have_a_utf8_safe_bounded_preview_and_exact_length() {
    let oversized_unit = format!("{}é", "x".repeat(63));
    let assessment = assess_program_risks(&[ProgramRiskObservation {
        id: ProgramRiskId::Pr001,
        value: 1.0,
        unit: &oversized_unit,
        samples: 5,
    }]);
    let row = &assessment.rows()[0];
    assert_eq!(row.status, AssessmentStatus::UnitMismatch);
    let preview = row.observed_unit.as_deref().expect("unit preview");
    assert_eq!(preview, "x".repeat(63));
    assert!(preview.len() <= MAX_OBSERVED_UNIT_PREVIEW_BYTES);
    assert_eq!(row.observed_unit_bytes, Some(oversized_unit.len()));
}

#[test]
fn g0_observation_order_does_not_change_assessment_artifact() {
    let mut forward: Vec<_> = ProgramRiskId::ALL
        .into_iter()
        .map(clear_observation)
        .collect();
    let expected = assess_program_risks(&forward);
    assert!(expected.all_clear());
    forward.reverse();
    let reversed = assess_program_risks(&forward);
    assert_eq!(expected.to_json(), reversed.to_json());
    assert_eq!(reversed.alert_count(), 0);
}

#[test]
fn g0_register_artifact_is_deterministic_and_uses_numeric_thresholds() {
    let first = program_risk_register_json();
    let second = program_risk_register_json();
    assert_eq!(first, second);
    assert!(first.starts_with("{\"schema\":\"frankensim.program-risk-register.v1\",\"risks\":["));
    assert_eq!(first.matches("\"id\":\"PR-").count(), 12);
    assert_eq!(first.matches("\"review_gate\":").count(), 12);
    assert_eq!(first.matches("\"contingency\":").count(), 12);
    assert_eq!(first.matches("\"threshold\":").count(), 12);
    assert_eq!(first.matches("\"domain\":").count(), 12);
    assert!(!first.contains("\"threshold\":\""));
}
