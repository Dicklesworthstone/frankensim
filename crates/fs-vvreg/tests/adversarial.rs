//! G0/G3 adversarial validation registry and honesty-tripwire tests.

use std::collections::BTreeSet;

use fs_vvreg::adversarial::{
    AdversarialCase, AdversarialEvidence, AdversarialEvidenceBasis, AdversarialOutcome,
    AdversarialRegistry, AdversarialRegistryError, AdversarialScorecardError, AttackedAssumption,
    DominantUncertainty, HonestyVerdict, adversarial_registry,
};

fn predicted(
    absolute_error: f64,
    allowed_error: f64,
    dominant: DominantUncertainty,
) -> AdversarialOutcome {
    AdversarialOutcome::Prediction {
        absolute_error,
        allowed_error,
        dominant,
    }
}

#[test]
fn seed_has_all_required_case_families_and_retained_bindings_resolve() {
    let registry = adversarial_registry();
    assert!(registry.cases().len() >= 6);
    let bases = registry
        .cases()
        .iter()
        .map(|case| case.evidence_basis)
        .collect::<BTreeSet<_>>();
    assert_eq!(
        bases,
        BTreeSet::from([
            AdversarialEvidenceBasis::Analytic,
            AdversarialEvidenceBasis::CrossCode,
            AdversarialEvidenceBasis::ControlledExperiment,
            AdversarialEvidenceBasis::InstrumentedRig,
        ])
    );
    assert!(
        registry
            .cases()
            .iter()
            .any(|case| case.evidence.is_retained())
    );
    assert!(
        registry
            .cases()
            .iter()
            .any(|case| !case.evidence.is_retained())
    );
}

#[test]
fn honest_prediction_refusal_and_demotion_pass() {
    let registry = adversarial_registry();
    let prediction = registry
        .assess(
            "contact-dominated-two-layer-stack",
            predicted(0.25, 0.25, DominantUncertainty::ContactResistance),
        )
        .unwrap();
    assert_eq!(prediction.verdict(), HonestyVerdict::Pass);
    assert!(!prediction.is_false_acceptance());

    let refusal = registry
        .assess(
            "fan-stall-multiple-operating-points",
            AdversarialOutcome::Refused {
                dominant: DominantUncertainty::FanOperatingPoint,
            },
        )
        .unwrap();
    assert_eq!(refusal.verdict(), HonestyVerdict::Pass);

    let demotion = registry
        .assess(
            "uncertain-blockable-vent-leakage",
            AdversarialOutcome::Demoted {
                dominant: DominantUncertainty::BoundaryCondition,
            },
        )
        .unwrap();
    assert_eq!(demotion.verdict(), HonestyVerdict::Pass);
}

#[test]
fn seeded_confident_wrong_prediction_trips_false_acceptance_counter() {
    let registry = adversarial_registry();
    let wrong = registry
        .assess(
            "radiation-dominated-low-flow-enclosure",
            predicted(
                0.500_000_000_000_000_1,
                0.5,
                DominantUncertainty::RadiationModel,
            ),
        )
        .unwrap();
    assert_eq!(wrong.verdict(), HonestyVerdict::Fail);
    assert!(wrong.is_false_acceptance());
    assert!(wrong.render_log().contains("false_acceptance=true"));

    let report = registry.render_regime_limitations(&[wrong]).unwrap();
    assert!(report.contains("false_acceptance_count: 1"));
    assert!(report.contains("| radiation-dominated-low-flow-enclosure |"));
}

#[test]
fn refusal_or_demotion_with_wrong_attribution_fails_without_false_acceptance() {
    let registry = adversarial_registry();
    for outcome in [
        AdversarialOutcome::Refused {
            dominant: DominantUncertainty::MaterialProperty,
        },
        AdversarialOutcome::Demoted {
            dominant: DominantUncertainty::FlowTopology,
        },
    ] {
        let receipt = registry
            .assess("fan-stall-multiple-operating-points", outcome)
            .unwrap();
        assert_eq!(receipt.verdict(), HonestyVerdict::Fail);
        assert!(!receipt.is_false_acceptance());
    }
}

#[test]
fn no_data_is_never_rendered_as_zero_or_pass() {
    let registry = adversarial_registry();
    let explicit = registry
        .assess(
            "natural-convection-cavity-reversal",
            AdversarialOutcome::NoData,
        )
        .unwrap();
    assert_eq!(explicit.verdict(), HonestyVerdict::NoData);
    let report = registry.render_regime_limitations(&[explicit]).unwrap();
    assert!(report.contains("NO-DATA:frankensim-extreal-program-f85xj.4.3"));
    assert!(report.contains("| natural-convection-cavity-reversal |"));
    assert!(!report.contains("false_acceptance_count: 1"));
}

#[test]
fn invalid_prediction_arithmetic_refuses_before_receipt() {
    let registry = adversarial_registry();
    for outcome in [
        predicted(f64::NAN, 1.0, DominantUncertainty::FlowTopology),
        predicted(1.0, f64::INFINITY, DominantUncertainty::FlowTopology),
        predicted(-0.0 - 1.0, 1.0, DominantUncertainty::FlowTopology),
        predicted(1.0, -0.0 - 1.0, DominantUncertainty::FlowTopology),
    ] {
        assert!(
            registry
                .assess("recirculation-behind-strip-fins", outcome)
                .is_err()
        );
    }
}

#[test]
fn registry_identity_is_order_stable_and_semantic_mutations_move_it() {
    let seed = adversarial_registry();
    let mut reversed = seed.cases().to_vec();
    reversed.reverse();
    assert_eq!(
        AdversarialRegistry::build(reversed).unwrap().identity(),
        seed.identity()
    );

    let mut changed = seed.cases().to_vec();
    changed[0].regime_limitation = "changed limitation";
    assert_ne!(
        AdversarialRegistry::build(changed).unwrap().identity(),
        seed.identity()
    );
}

#[test]
fn retained_evidence_basis_mismatch_refuses() {
    let case = AdversarialCase {
        id: "mismatched-coordinate",
        title: "Mismatched coordinate",
        regime: "synthetic",
        attacked_assumption: AttackedAssumption::AttachedFlow,
        expected_dominant_uncertainty: DominantUncertainty::FlowTopology,
        evidence_basis: AdversarialEvidenceBasis::CrossCode,
        evidence: AdversarialEvidence::Retained {
            dataset_id: "thermal-a-contact-series",
        },
        regime_limitation: "Must refuse rather than relabel analytic evidence as cross-code.",
    };
    assert!(matches!(
        AdversarialRegistry::build(vec![case]),
        Err(AdversarialRegistryError::EvidenceBasisMismatch { .. })
    ));
}

#[test]
fn scorecard_refuses_duplicate_and_foreign_receipts() {
    let registry = adversarial_registry();
    let receipt = registry
        .assess(
            "contact-dominated-two-layer-stack",
            AdversarialOutcome::Refused {
                dominant: DominantUncertainty::ContactResistance,
            },
        )
        .unwrap();
    assert!(matches!(
        registry.render_regime_limitations(&[receipt.clone(), receipt.clone()]),
        Err(AdversarialScorecardError::DuplicateAssessment { .. })
    ));

    let mut changed_cases = registry.cases().to_vec();
    changed_cases[0].regime_limitation = "foreign registry";
    let changed = AdversarialRegistry::build(changed_cases).unwrap();
    assert!(matches!(
        changed.render_regime_limitations(&[receipt]),
        Err(AdversarialScorecardError::ForeignAssessment { .. })
    ));
}

#[test]
fn report_is_byte_deterministic_and_exposes_every_regime_limitation() {
    let registry = adversarial_registry();
    let receipts = vec![
        registry
            .assess(
                "recirculation-behind-strip-fins",
                AdversarialOutcome::Demoted {
                    dominant: DominantUncertainty::FlowTopology,
                },
            )
            .unwrap(),
        registry
            .assess(
                "biot-extremes-lumped-breakdown",
                predicted(0.0, 0.0, DominantUncertainty::SpatialTemperature),
            )
            .unwrap(),
    ];
    let first = registry.render_regime_limitations(&receipts).unwrap();
    let second = registry.render_regime_limitations(&receipts).unwrap();
    assert_eq!(first, second);
    for case in registry.cases() {
        assert!(first.contains(case.id));
        assert!(first.contains(case.regime));
        assert!(first.contains(case.regime_limitation));
    }
}
