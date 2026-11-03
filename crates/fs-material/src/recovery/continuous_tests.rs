use super::*;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use std::collections::BTreeMap;

fn source(uniform: bool, rate: f64) -> ConstitutiveModelCard {
    let declaration = format!(
        "synthetic stationary local-rate fixture: uniform={uniform}; k={rate}; Rinf=0.9; equal ideal tanks; no measured calibration"
    );
    let hash = hash_domain("continuous-recovery-test", declaration.as_bytes());
    ConstitutiveModelCard {
        law: LawId(
            if uniform {
                CONTINUOUS_UNIFORM_RATE_LAW
            } else {
                CONTINUOUS_SINGLE_RATE_LAW
            }
            .into(),
        ),
        law_version: 1,
        parameters: BTreeMap::from([
            (
                super::super::ULTIMATE_RECOVERY_PARAMETER.into(),
                LawParameter {
                    value: 0.9,
                    dims: Dims::NONE,
                },
            ),
            (
                if uniform {
                    super::super::MAXIMUM_RATE_PARAMETER
                } else {
                    super::super::RATE_PARAMETER
                }
                .into(),
                LawParameter {
                    value: rate,
                    dims: Frequency::DIMS,
                },
            ),
        ]),
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                MEAN_RESIDENCE_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                1.0,
                100_000.0,
            )
            .with_quantity(
                TANK_COUNT_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
                64.0,
            )
            .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 9.0, 11.0),
        sources: vec![hash],
        provenance: Provenance {
            source: declaration,
            license: "project synthetic fixture".into(),
            artifact: Some(hash),
        },
    }
}
fn point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 10.0)
        .unwrap()
}
fn evaluate(card: &ContinuousRecoveryCard, n: usize, mean: f64, feed: f64) -> ContinuousTreatment {
    card.evaluate_train(
        n,
        Volume::new(mean),
        VolumetricFlowRate::new(1.0),
        MassFlowRate::new(feed),
        &point(),
    )
    .unwrap_or_else(|error| {
        panic!(
            "train n={n} mean={mean} parameters={:?}: {error:?}",
            card.source().parameters
        )
    })
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

fn gamma_source(shape: f64, scale: f64) -> ConstitutiveModelCard {
    let mut card = source(false, scale);
    card.law = LawId(CONTINUOUS_GAMMA_RATE_LAW.into());
    card.parameters.remove(super::super::RATE_PARAMETER);
    card.parameters.insert(
        super::super::GAMMA_SCALE_PARAMETER.into(),
        LawParameter {
            value: scale,
            dims: Frequency::DIMS,
        },
    );
    card.parameters.insert(
        super::super::GAMMA_SHAPE_PARAMETER.into(),
        LawParameter {
            value: shape,
            dims: Dims::NONE,
        },
    );
    card
}

#[test]
fn gamma_cstr_matches_independent_ratewise_balances_including_singular_shape() {
    // Integrate RATE populations, independently of production's log-residence
    // integral. For shape=1/2, k=theta*u² removes the rate-density singularity;
    // shape=2 has density u*exp(-u). Neither oracle uses candidate outputs.
    for shape in [0.5, 2.0] {
        let original = gamma_source(shape, 0.02);
        let card = ContinuousRecoveryCard::from_card(original.clone()).unwrap();
        for tanks in [1, 2, 8, 64] {
            let samples = 16384;
            let upper = if shape == 0.5 { 12.0 } else { 64.0 };
            let mut recovered = 0.0;
            let mut surviving = 0.0;
            for j in 0..=samples {
                let u = upper * f64::from(j) / f64::from(samples);
                let (rate, density) = if shape == 0.5 {
                    (
                        0.02 * u * u,
                        2.0 / std::f64::consts::PI.sqrt() * (-u * u).exp(),
                    )
                } else {
                    (0.02 * u, u * (-u).exp())
                };
                let mut survivor = 1.0;
                for _ in 0..tanks {
                    survivor /= 1.0 + rate * 100.0 / tanks as f64;
                }
                let weight = if j == 0 || j == samples {
                    1.0
                } else if j % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                recovered += weight * density * (1.0 - survivor);
                surviving += weight * density * survivor;
            }
            let result = evaluate(&card, tanks, 100.0, 10.0);
            close(
                result.recovery.value(),
                0.9 * recovered * upper / (3.0 * f64::from(samples)),
                1e-10,
            );
            close(
                result.unrecovered_fraction.value(),
                0.1 + 0.9 * surviving * upper / (3.0 * f64::from(samples)),
                1e-10,
            );
            assert_eq!(result.source_card, original);
            close(result.conservation_residual.value(), 0.0, 2e-15);
            let diagnostics = result.gamma_integration.unwrap();
            assert!(diagnostics.evaluations > 0 && diagnostics.evaluations <= 65_536);
            assert!(diagnostics.absolute_error_estimate.value() < 1e-11);
            assert!(diagnostics.omitted_tail_bound.value() < 1e-14);
        }
    }
}

#[test]
fn gamma_cstr_residence_distribution_and_narrow_rate_limit_remain_distinct() {
    let card = ContinuousRecoveryCard::from_card(gamma_source(0.5, 0.04)).unwrap();
    let one = evaluate(&card, 1, 100.0, 1.0);
    let eight = evaluate(&card, 8, 100.0, 1.0);
    let batch_at_mean = 0.9 * (1.0 - 1.0 / 5.0_f64.sqrt());
    assert!(one.recovery.value() < eight.recovery.value());
    assert!(eight.recovery.value() < batch_at_mean);
    assert!((eight.recovery.value() - batch_at_mean).abs() > 0.005);
    let concentrated = ContinuousRecoveryCard::from_card(gamma_source(1e6, 0.02 / 1e6)).unwrap();
    let single = ContinuousRecoveryCard::from_card(source(false, 0.02)).unwrap();
    for tanks in [1, 8, 64] {
        close(
            evaluate(&concentrated, tanks, 100.0, 1.0).recovery.value(),
            evaluate(&single, tanks, 100.0, 1.0).recovery.value(),
            4e-7,
        );
    }
    let twice_flow = card
        .evaluate_train(
            8,
            Volume::new(100.0),
            VolumetricFlowRate::new(2.0),
            MassFlowRate::new(1.0),
            &point(),
        )
        .unwrap();
    assert!(twice_flow.recovery.value() < eight.recovery.value());
}

#[test]
fn gamma_cstr_preserves_tiny_responses_and_refuses_unresolved_numerics() {
    let tiny = ContinuousRecoveryCard::from_card(gamma_source(2.0, 1e-200)).unwrap();
    for tanks in [1, 8] {
        let result = evaluate(&tiny, tanks, 1.0, 1e200);
        close(result.recovered_metal_flow.value(), 1.8, 2e-12);
        assert!(result.recovery.value() > 0.0);
    }
    let huge_shape = ContinuousRecoveryCard::from_card(gamma_source(1e300, 1e-300)).unwrap();
    close(
        evaluate(&huge_shape, 1, 1.0, 1.0).recovery.value(),
        0.45,
        2e-12,
    );
    let mut tail_source = gamma_source(0.1, f64::MAX);
    tail_source
        .parameters
        .get_mut(super::super::ULTIMATE_RECOVERY_PARAMETER)
        .unwrap()
        .value = 1.0;
    let tail = ContinuousRecoveryCard::from_card(tail_source).unwrap();
    let result = evaluate(&tail, 1, 100.0, 1.0);
    assert!(result.unrecovered_fraction.value() > 0.0);
    assert!(result.unrecovered_fraction.value() < 1e-30);
    // This positive survival is beyond binary64; it must not become zero.
    let impossible = ContinuousRecoveryCard::from_card(gamma_source(f64::MAX, f64::MAX)).unwrap();
    assert!(matches!(
        impossible.evaluate_train(
            64,
            Volume::new(1.0),
            VolumetricFlowRate::new(1.0),
            MassFlowRate::new(1.0),
            &point()
        ),
        Err(RecoveryError::Arithmetic(_))
    ));
}

#[test]
fn gamma_requires_its_own_stationary_source_and_complete_support() {
    for shape in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(ContinuousRecoveryCard::from_card(gamma_source(shape, 0.02)).is_err());
    }
    for mutate in 0..5 {
        let mut source = gamma_source(2.0, 0.02);
        match mutate {
            0 => source.law = LawId(super::super::GAMMA_RATE_LAW.into()),
            1 => {
                source
                    .parameters
                    .remove(super::super::GAMMA_SHAPE_PARAMETER);
            }
            2 => {
                source
                    .parameters
                    .get_mut(super::super::GAMMA_SHAPE_PARAMETER)
                    .unwrap()
                    .dims = Frequency::DIMS
            }
            3 => {
                source
                    .parameters
                    .get_mut(super::super::GAMMA_SCALE_PARAMETER)
                    .unwrap()
                    .dims = Time::DIMS
            }
            _ => source.sources.clear(),
        }
        assert!(ContinuousRecoveryCard::from_card(source).is_err());
    }
    let model = ContinuousRecoveryCard::from_card(gamma_source(2.0, 0.02)).unwrap();
    assert_eq!(
        model.evaluate_train(
            8,
            Volume::new(100.0),
            VolumetricFlowRate::new(1.0),
            MassFlowRate::new(1.0),
            &QueryPoint::new()
        ),
        Err(RecoveryError::OutsideValidity)
    );
    let zero = ContinuousRecoveryCard::from_card(gamma_source(2.0, 0.0)).unwrap();
    let result = evaluate(&zero, 8, 100.0, 10.0);
    assert_eq!(result.recovery.value(), 0.0);
    assert_eq!(result.gamma_integration, None);
}

#[test]
fn single_tank_closes_actual_reaction_balance_and_retains_source() {
    let original = source(false, 0.02);
    let card = ContinuousRecoveryCard::from_card(original.clone()).unwrap();
    let result = evaluate(&card, 1, 100.0, 10.0);
    close(result.recovery.value(), 0.6, 2e-16);
    close(result.recovered_metal_flow.value(), 6.0, 2e-15);
    close(result.tail_metal_flow.value(), 4.0, 2e-15);
    // Independent Q*C0-Q*C-k*V*C=0 for the recoverable population.
    let recoverable_exit_concentration = 9.0 / (1.0 + 0.02 * 100.0);
    close(
        result.recovered_metal_flow.value(),
        0.02 * 100.0 * recoverable_exit_concentration,
        2e-15,
    );
    close(result.conservation_residual.value(), 0.0, 2e-15);
    assert_eq!(result.source_card, original);
    assert_eq!(result.conditions.axes()[MEAN_RESIDENCE_TIME_AXIS], 100.0);
    assert_eq!(result.conditions.axes()[TANK_COUNT_AXIS], 1.0);
}

#[test]
fn staged_recovery_is_not_batch_at_the_mean() {
    let card = ContinuousRecoveryCard::from_card(source(false, 0.02)).unwrap();
    let one = evaluate(&card, 1, 100.0, 10.0);
    let eight = evaluate(&card, 8, 100.0, 10.0);
    close(eight.recovery.value(), 0.749005056, 3e-16);
    let batch_at_mean = 0.9 * (-(-2.0_f64).exp_m1());
    assert!(one.recovery.value() < eight.recovery.value());
    assert!(eight.recovery.value() < batch_at_mean);
    assert!((eight.recovery.value() - batch_at_mean).abs() > 0.02);
}

#[test]
fn uniform_population_matches_independent_ratewise_tank_balances() {
    // Equal mean rates: a single k=0.02/s versus uniform 0..0.04/s.
    let single = ContinuousRecoveryCard::from_card(source(false, 0.02)).unwrap();
    let uniform = ContinuousRecoveryCard::from_card(source(true, 0.04)).unwrap();
    for tanks in [1, 2, 8, 64] {
        let sample_count = 8192;
        let mut sum = 0.0;
        for j in 0..=sample_count {
            let k = 0.04 * f64::from(j) / f64::from(sample_count);
            let mut exiting = 1.0;
            for _ in 0..tanks {
                exiting /= 1.0 + k * 100.0 / tanks as f64;
            }
            let weight = if j == 0 || j == sample_count {
                1.0
            } else if j % 2 == 0 {
                2.0
            } else {
                4.0
            };
            sum += weight * (1.0 - exiting);
        }
        let result = evaluate(&uniform, tanks, 100.0, 10.0);
        close(
            result.recovery.value(),
            0.9 * sum / (3.0 * f64::from(sample_count)),
            6e-13,
        );
        assert!(result.recovery.value() < evaluate(&single, tanks, 100.0, 10.0).recovery.value());
        close(
            result.recovered_metal_flow.value() + result.tail_metal_flow.value(),
            10.0,
            4e-15,
        );
    }
    let uniform = ContinuousRecoveryCard::from_card(source(true, 0.02)).unwrap();
    close(
        evaluate(&uniform, 1, 100.0, 1.0).recovery.value(),
        0.4056244700993507,
        2e-16,
    );
    close(
        evaluate(&uniform, 2, 100.0, 1.0).recovery.value(),
        0.45,
        2e-16,
    );
    close(
        evaluate(&uniform, 8, 100.0, 1.0).recovery.value(),
        0.49356781714285713,
        3e-16,
    );
    for a in [0.01_f64, 0.12499999999999999, 0.125] {
        let card = ContinuousRecoveryCard::from_card(source(true, a)).unwrap();
        close(
            evaluate(&card, 1, 1.0, 1.0).recovery.value(),
            0.9 * (1.0 - a.ln_1p() / a),
            1e-15,
        );
    }
}

#[test]
fn tiny_conversion_and_large_rate_tails_remain_positive() {
    for uniform in [false, true] {
        let card = ContinuousRecoveryCard::from_card(source(uniform, 1e-200)).unwrap();
        for n in [1, 8] {
            let result = evaluate(&card, n, 1.0, 1e200);
            close(
                result.recovered_metal_flow.value(),
                if uniform { 0.45 } else { 0.9 },
                5e-16,
            );
            assert!(result.recovery.value() > 0.0);
        }
        let mut source = source(uniform, f64::MAX);
        source
            .parameters
            .get_mut(super::super::ULTIMATE_RECOVERY_PARAMETER)
            .unwrap()
            .value = 1.0;
        let card = ContinuousRecoveryCard::from_card(source).unwrap();
        let result = evaluate(&card, 1, 1.0, 1e308);
        assert!(result.tail_metal_flow.value() > 0.0);
        assert!(result.unrecovered_fraction.value() > 0.0);
    }
}

#[test]
fn actual_flow_controls_residence_and_stopped_flow_refuses() {
    let card = ContinuousRecoveryCard::from_card(source(false, 0.02)).unwrap();
    let original = evaluate(&card, 8, 100.0, 1.0);
    let half_flow = card
        .evaluate_train(
            8,
            Volume::new(100.0),
            VolumetricFlowRate::new(0.5),
            MassFlowRate::new(1.0),
            &point(),
        )
        .unwrap();
    assert_eq!(half_flow.mean_residence_time.value(), 200.0);
    assert!(half_flow.recovery.value() > original.recovery.value());
    let scaled = card
        .evaluate_train(
            8,
            Volume::new(200.0),
            VolumetricFlowRate::new(2.0),
            MassFlowRate::new(1.0),
            &point(),
        )
        .unwrap();
    assert_eq!(scaled.recovery, original.recovery);
    assert!(
        card.evaluate_train(
            8,
            Volume::new(100.0),
            VolumetricFlowRate::new(0.0),
            MassFlowRate::new(1.0),
            &point()
        )
        .is_err()
    );
    assert_eq!(
        evaluate(&card, 8, 100.0, 0.0).recovered_metal_flow.value(),
        0.0
    );
}

#[test]
fn complete_operating_validity_and_derived_coordinates_are_enforced() {
    let card = ContinuousRecoveryCard::from_card(source(false, 0.02)).unwrap();
    for conditions in [
        QueryPoint::new(),
        QueryPoint::new()
            .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 12.0)
            .unwrap(),
        QueryPoint::new()
            .with_quantity("pH", QuantitySpec::dimensional(Time::DIMS), 10.0)
            .unwrap(),
    ] {
        assert_eq!(
            card.evaluate_train(
                8,
                Volume::new(100.0),
                VolumetricFlowRate::new(1.0),
                MassFlowRate::new(1.0),
                &conditions
            ),
            Err(RecoveryError::OutsideValidity)
        );
    }
    assert_eq!(
        card.evaluate_train(
            8,
            Volume::new(100_001.0),
            VolumetricFlowRate::new(1.0),
            MassFlowRate::new(1.0),
            &point()
        ),
        Err(RecoveryError::OutsideValidity)
    );
    let conflict = point()
        .with_quantity(
            MEAN_RESIDENCE_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            99.0,
        )
        .unwrap();
    assert!(
        card.evaluate_train(
            8,
            Volume::new(100.0),
            VolumetricFlowRate::new(1.0),
            MassFlowRate::new(1.0),
            &conflict
        )
        .is_err()
    );
}

#[test]
fn batch_cards_and_unbacked_or_untyped_continuous_declarations_refuse() {
    for change in 0..9 {
        let mut card = source(false, 0.02);
        match change {
            0 => card.law = LawId(super::super::SINGLE_RATE_LAW.into()),
            1 => {
                card.validity = card.validity.with_quantity(
                    ELAPSED_TIME_AXIS,
                    QuantitySpec::dimensional(Time::DIMS),
                    0.0,
                    100.0,
                )
            }
            2 => card.sources.clear(),
            3 => card.provenance.license.clear(),
            4 => {
                card.parameters
                    .get_mut(super::super::RATE_PARAMETER)
                    .unwrap()
                    .dims = Time::DIMS
            }
            5 => card.law_version = 2,
            6 => {
                card.parameters
                    .get_mut(super::super::ULTIMATE_RECOVERY_PARAMETER)
                    .unwrap()
                    .value = 90.0
            }
            7 => {
                card.validity = ValidityDomain::unconstrained()
                    .with(MEAN_RESIDENCE_TIME_AXIS, 1.0, 100.0)
                    .with(TANK_COUNT_AXIS, 1.0, 64.0)
            }
            _ => card.initial_state = InitialStatePolicy::RequiresDeclaredState,
        }
        assert!(
            ContinuousRecoveryCard::from_card(card).is_err(),
            "change {change}"
        );
    }
}

#[test]
fn zero_rate_and_zero_recoverable_fraction_and_arithmetic_refusals() {
    for uniform in [false, true] {
        let card = ContinuousRecoveryCard::from_card(source(uniform, 0.0)).unwrap();
        let result = evaluate(&card, 8, 100.0, 10.0);
        assert_eq!(result.recovery.value(), 0.0);
        assert_eq!(result.tail_metal_flow.value(), 10.0);
        let mut no_population = source(uniform, 0.02);
        no_population
            .parameters
            .get_mut(super::super::ULTIMATE_RECOVERY_PARAMETER)
            .unwrap()
            .value = 0.0;
        assert_eq!(
            evaluate(
                &ContinuousRecoveryCard::from_card(no_population).unwrap(),
                8,
                100.0,
                10.0
            )
            .recovery
            .value(),
            0.0
        );
    }
    let card = ContinuousRecoveryCard::from_card(source(false, 0.02)).unwrap();
    for (n, v, q, m) in [
        (0, 1.0, 1.0, 1.0),
        (65, 1.0, 1.0, 1.0),
        (1, f64::NAN, 1.0, 1.0),
        (1, 1.0, f64::INFINITY, 1.0),
        (1, 1.0, 1.0, -1.0),
        (1, 1.0, 1.0, f64::NAN),
        (1, f64::MAX, f64::MIN_POSITIVE, 1.0),
    ] {
        assert!(
            card.evaluate_train(
                n,
                Volume::new(v),
                VolumetricFlowRate::new(q),
                MassFlowRate::new(m),
                &point()
            )
            .is_err()
        );
    }
    let mut extreme = source(false, f64::MAX);
    extreme
        .parameters
        .get_mut(super::super::ULTIMATE_RECOVERY_PARAMETER)
        .unwrap()
        .value = 1.0;
    let extreme = ContinuousRecoveryCard::from_card(extreme).unwrap();
    assert!(
        extreme
            .evaluate_train(
                8,
                Volume::new(1.0),
                VolumetricFlowRate::new(1.0),
                MassFlowRate::new(1.0),
                &point()
            )
            .is_err()
    );
    assert!(
        extreme
            .evaluate_train(
                1,
                Volume::new(1.0),
                VolumetricFlowRate::new(1.0),
                MassFlowRate::new(f64::MIN_POSITIVE),
                &point()
            )
            .is_err()
    );
}
