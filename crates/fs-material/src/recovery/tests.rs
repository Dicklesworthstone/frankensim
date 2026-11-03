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

fn gamma_source(shape: f64, scale: f64) -> ConstitutiveModelCard {
    let mut card = source(false);
    card.law = LawId(GAMMA_RATE_LAW.into());
    card.parameters.remove(RATE_PARAMETER);
    card.parameters.insert(
        GAMMA_SCALE_PARAMETER.into(),
        LawParameter {
            value: scale,
            dims: Frequency::DIMS,
        },
    );
    card.parameters.insert(
        GAMMA_SHAPE_PARAMETER.into(),
        LawParameter {
            value: shape,
            dims: Dims::NONE,
        },
    );
    card
}

#[test]
fn gamma_batch_matches_exact_laplace_transform_and_retains_shape_scale() {
    for shape in [0.25_f64, 0.5, 1.0, 2.0, 12.0] {
        let original = gamma_source(shape, 0.02);
        let law = RecoveryKinetics::from_card(original.clone()).unwrap();
        assert_eq!(law.distribution(), RecoveryRateDistribution::Gamma);
        assert_eq!(law.gamma_shape().unwrap().value(), shape);
        assert_eq!(law.rate().value(), 0.02);
        assert_eq!(law.source(), &original);
        assert_eq!(value(&law, 0.0), 0.0);
        for t in [1e-10_f64, 1.0, 100.0, 1000.0] {
            let expected = 0.9 * (-(-shape * (0.02 * t).ln_1p()).exp_m1());
            close(value(&law, t), expected, 2e-14 * expected);
        }
    }
    // Shape one is the exponential-rate population, not a single rate.
    let law = RecoveryKinetics::from_card(gamma_source(1.0, 0.02)).unwrap();
    close(value(&law, 100.0), 0.6, 3e-16);
    assert!(
        value(&law, 100.0) < value(&RecoveryKinetics::from_card(source(false)).unwrap(), 100.0)
    );
}

#[test]
fn gamma_batch_refuses_invalid_shape_units_and_unrepresentable_recovery() {
    for shape in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(RecoveryKinetics::from_card(gamma_source(shape, 0.02)).is_err());
    }
    let mut wrong = gamma_source(2.0, 0.02);
    wrong
        .parameters
        .get_mut(GAMMA_SHAPE_PARAMETER)
        .unwrap()
        .dims = Time::DIMS;
    assert!(RecoveryKinetics::from_card(wrong).is_err());
    let mut wrong = gamma_source(2.0, 0.02);
    wrong
        .parameters
        .get_mut(GAMMA_SCALE_PARAMETER)
        .unwrap()
        .dims = Time::DIMS;
    assert!(RecoveryKinetics::from_card(wrong).is_err());
    let zero = RecoveryKinetics::from_card(gamma_source(2.0, 0.0)).unwrap();
    assert_eq!(value(&zero, 100.0), 0.0);
    let law = RecoveryKinetics::from_card(gamma_source(1e-300, 1e-300)).unwrap();
    assert!(matches!(
        law.evaluate(Time::new(1.0), &conditions()),
        Err(RecoveryError::Arithmetic(_))
    ));
    // Separate logarithms retain a representable response when scale*time
    // would overflow, or scale*time would underflow before shape is applied.
    let law = RecoveryKinetics::from_card(gamma_source(1e-300, f64::MAX)).unwrap();
    assert!(value(&law, 100.0) > 0.0);
    let law = RecoveryKinetics::from_card(gamma_source(1e300, 1e-300)).unwrap();
    close(value(&law, 1.0), 0.9 * (1.0 - (-1.0_f64).exp()), 4e-14);
}
