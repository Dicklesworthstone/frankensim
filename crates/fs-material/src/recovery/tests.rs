use super::*;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use std::collections::BTreeMap;

fn source(uniform: bool) -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId(
            if uniform {
                UNIFORM_RATE_LAW
            } else {
                SINGLE_RATE_LAW
            }
            .into(),
        ),
        law_version: 1,
        parameters: BTreeMap::from([
            (
                ULTIMATE_RECOVERY_PARAMETER.into(),
                LawParameter {
                    value: 0.9,
                    dims: Dims::NONE,
                },
            ),
            (
                if uniform {
                    MAXIMUM_RATE_PARAMETER
                } else {
                    RATE_PARAMETER
                }
                .into(),
                LawParameter {
                    value: 0.02,
                    dims: Frequency::DIMS,
                },
            ),
        ]),
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                ELAPSED_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                0.0,
                100_000.0,
            )
            .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 9.0, 11.0),
        sources: vec![hash_domain(
            "recovery-synthetic-test",
            b"authored fixture Rinf=0.9 k=0.02/s",
        )],
        provenance: Provenance {
            source: "synthetic curve; not an ore test".into(),
            license: "project fixture".into(),
            artifact: None,
        },
    }
}
fn conditions() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 10.0)
        .unwrap()
}
fn value(law: &RecoveryKinetics, t: f64) -> f64 {
    law.evaluate(Time::new(t), &conditions()).unwrap().value()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn single_rate_known_answer_and_limits() {
    let law = RecoveryKinetics::from_card(source(false)).unwrap();
    close(value(&law, 100.0), 0.7781982450870486, 2e-15);
    assert_eq!(value(&law, 0.0), 0.0);
    assert_eq!(value(&law, 100_000.0), 0.9);
    assert!(value(&law, 50.0) < value(&law, 100.0));
    close(value(&law, 1e-10), 0.9 * 0.02 * 1e-10, 1e-23);
    assert_eq!(law.source(), &source(false));
}

#[test]
fn uniform_rate_matches_independent_simpson_integral() {
    let law = RecoveryKinetics::from_card(source(true)).unwrap();
    for time in [1e-9, 0.5, 6.24, 6.25, 100.0, 500.0] {
        // Independent integration of individual rate populations across [0,kmax].
        let n = 2048;
        let mut sum = 0.0;
        for i in 0..=n {
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let k = 0.02 * f64::from(i) / f64::from(n);
            sum += weight * -(-k * time).exp_m1();
        }
        close(value(&law, time), 0.9 * sum / (3.0 * f64::from(n)), 3e-12);
    }
    close(value(&law, 100.0), 0.5109008774564757, 2e-15);
    assert!(
        value(&law, 100.0) < value(&RecoveryKinetics::from_card(source(false)).unwrap(), 100.0)
    );
    close(value(&law, 1e-10), 0.9 * 0.02 * 1e-10 / 2.0, 1e-23);
}

#[test]
fn zero_rate_or_zero_recoverable_population_is_explicit() {
    for uniform in [false, true] {
        for name in [
            ULTIMATE_RECOVERY_PARAMETER,
            if uniform {
                MAXIMUM_RATE_PARAMETER
            } else {
                RATE_PARAMETER
            },
        ] {
            let mut card = source(uniform);
            card.parameters.get_mut(name).unwrap().value = 0.0;
            assert_eq!(
                value(&RecoveryKinetics::from_card(card).unwrap(), 100.0),
                0.0
            );
        }
    }
}

#[test]
fn complete_typed_applicability_is_required() {
    let law = RecoveryKinetics::from_card(source(false)).unwrap();
    for point in [
        QueryPoint::new(),
        QueryPoint::new()
            .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 12.0)
            .unwrap(),
        QueryPoint::new()
            .with_quantity("pH", QuantitySpec::dimensional(Time::DIMS), 10.0)
            .unwrap(),
    ] {
        assert_eq!(
            law.evaluate(Time::new(100.0), &point),
            Err(RecoveryError::OutsideValidity)
        );
    }
    assert_eq!(
        law.evaluate(Time::new(100_001.0), &conditions()),
        Err(RecoveryError::OutsideValidity)
    );
    assert!(
        law.evaluate(
            Time::new(100.0),
            &conditions()
                .with_quantity(
                    ELAPSED_TIME_AXIS,
                    QuantitySpec::dimensional(Time::DIMS),
                    99.0
                )
                .unwrap()
        )
        .is_err()
    );
    assert!(law.evaluate(Time::new(f64::NAN), &conditions()).is_err());
}

#[test]
fn wrong_laws_dimensions_state_and_unbacked_cards_refuse() {
    for mutation in 0..8 {
        let mut card = source(false);
        match mutation {
            0 => card.law_version = 2,
            1 => card.sources.clear(),
            2 => card.provenance.source.clear(),
            3 => card.parameters.get_mut(RATE_PARAMETER).unwrap().dims = Time::DIMS,
            4 => {
                card.parameters
                    .get_mut(ULTIMATE_RECOVERY_PARAMETER)
                    .unwrap()
                    .value = 90.0
            }
            5 => card.parameters.get_mut(RATE_PARAMETER).unwrap().value = -0.01,
            6 => card.initial_state = InitialStatePolicy::RequiresDeclaredState,
            _ => {
                card.validity = ValidityDomain::unconstrained().with(ELAPSED_TIME_AXIS, 0.0, 100.0)
            }
        }
        assert!(
            RecoveryKinetics::from_card(card).is_err(),
            "mutation {mutation}"
        );
    }
}

#[test]
fn arithmetic_range_refuses_and_source_identity_changes() {
    let mut card = source(false);
    let original = card.content_hash();
    card.parameters.get_mut(RATE_PARAMETER).unwrap().value = f64::MAX;
    assert_ne!(original, card.content_hash());
    let law = RecoveryKinetics::from_card(card).unwrap();
    assert!(matches!(
        law.evaluate(Time::new(100.0), &conditions()),
        Err(RecoveryError::Arithmetic(_))
    ));
}
