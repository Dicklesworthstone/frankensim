use super::*;
use fs_qty::Length;

fn context() -> RecoveryTestContext {
    RecoveryTestContext {
        specimen: "synthetic recoverable component".into(),
        method: "authored analytical fixture, not measured plant calibration".into(),
        provenance: Provenance {
            source: "independent mathematical known answers".into(),
            license: "project fixture".into(),
            artifact: None,
        },
    }
}

fn domain() -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity(
        RESIDENCE_TIME_AXIS,
        RESIDENCE_TIME_QUANTITY,
        0.0,
        f64::MAX,
    )
}

fn card(rates: RateDistribution, reactor: RecoveryReactor) -> RecoveryCard {
    RecoveryCard::new(Dimensionless::new(0.8), rates, reactor, context(), domain()).unwrap()
}

fn single() -> RateDistribution {
    RateDistribution::Single {
        rate: Frequency::new(0.01),
    }
}

fn recovery(card: &RecoveryCard, time: f64) -> RecoveryResponse {
    card.evaluate(Time::new(time), &QueryPoint::new()).unwrap()
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn batch_and_plug_flow_match_independent_known_answers() {
    for reactor in [RecoveryReactor::Batch, RecoveryReactor::PlugFlow] {
        let fixed = card(single(), reactor);
        let original = fixed.clone();
        close(
            recovery(&fixed, 100.0).recovered.value(),
            0.8 * (1.0 - (-1.0_f64).exp()),
            2e-15,
        );
        let uniform = card(
            RateDistribution::Rectangular {
                maximum_rate: Frequency::new(0.02),
            },
            reactor,
        );
        // At k_max*t=2: 1-(1-exp(-2))/2 = (1+exp(-2))/2.
        close(
            recovery(&uniform, 100.0).recovered.value(),
            0.4 * (1.0 + (-2.0_f64).exp()),
            2e-15,
        );
        let gamma = card(
            RateDistribution::Gamma {
                shape: Dimensionless::new(2.0),
                scale: Frequency::new(0.005),
            },
            reactor,
        );
        // Gamma shape 2, scale*t=1: recovery = R_inf * 3/4.
        close(recovery(&gamma, 200.0).recovered.value(), 0.6, 2e-15);
        assert_eq!(fixed, original);
        assert_eq!(fixed.test_context(), &context());
        assert_eq!(fixed.validity(), &domain());
    }
}

#[test]
fn mixed_tanks_have_distinct_recovery_and_converge_to_plug_flow() {
    let make = |tanks| card(single(), RecoveryReactor::EqualCstrCascade { tanks });
    let one = recovery(&make(1), 100.0).recovered.value();
    let two = recovery(&make(2), 100.0).recovered.value();
    close(one, 0.4, 2e-15); // 1 - 1/(1+1).
    close(two, 0.8 * (1.0 - 4.0 / 9.0), 2e-15); // TOTAL residence is 100 s.
    let many = recovery(&make(100_000), 100.0).recovered.value();
    let plug = recovery(&card(single(), RecoveryReactor::PlugFlow), 100.0)
        .recovered
        .value();
    assert!(one < two && two < many && many < plug);
    assert!(plug - many < 2e-6);
}

#[test]
fn discrete_populations_preserve_component_weights_through_the_cascade() {
    let make = |weight: f64| {
        card(
            RateDistribution::Discrete {
                populations: vec![
                    RatePopulation {
                        weight: Dimensionless::new(weight),
                        rate: Frequency::new(0.01),
                    },
                    RatePopulation {
                        weight: Dimensionless::new(3.0 * weight),
                        rate: Frequency::new(0.03),
                    },
                ],
            },
            RecoveryReactor::EqualCstrCascade { tanks: 2 },
        )
    };
    let expected = 0.8 * (0.25 * (1.0 - 4.0 / 9.0) + 0.75 * (1.0 - 4.0 / 25.0));
    close(
        recovery(&make(1.0), 100.0).recovered.value(),
        expected,
        3e-15,
    );
    close(
        recovery(&make(1e300), 100.0).recovered.value(),
        expected,
        3e-15,
    );
    let average_rate = card(
        RateDistribution::Single {
            rate: Frequency::new(0.025),
        },
        RecoveryReactor::EqualCstrCascade { tanks: 2 },
    );
    assert!((recovery(&average_rate, 100.0).recovered.value() - expected).abs() > 0.01);
}

#[test]
fn disparate_population_weights_preserve_partition_closure_or_refuse() {
    let mut populations = vec![RatePopulation {
        weight: Dimensionless::new(1.0),
        rate: Frequency::new(1e-30),
    }];
    populations.extend((0..255).map(|_| RatePopulation {
        weight: Dimensionless::new(1e-16),
        rate: Frequency::new(100.0),
    }));
    let model = card(
        RateDistribution::Discrete { populations },
        RecoveryReactor::Batch,
    );
    let response = recovery(&model, 1.0);
    let small_total = 255.0e-16;
    close(
        response.recovered.value(),
        0.8 * small_total / (1.0 + small_total),
        3e-15,
    );
    let residue =
        response.recoverable_remaining.value() + response.intrinsically_unrecoverable.value();
    assert!((response.recovered.value() + residue - 1.0).abs() <= 2.0 * f64::EPSILON);
    crate::separation::Partition::from_fractions(response.recovered, Dimensionless::new(residue))
        .unwrap();

    let extreme = RateDistribution::Discrete {
        populations: vec![
            RatePopulation {
                weight: Dimensionless::new(f64::MAX),
                rate: Frequency::new(1.0),
            },
            RatePopulation {
                weight: Dimensionless::new(f64::from_bits(1)),
                rate: Frequency::new(1.0),
            },
        ],
    };
    assert!(
        RecoveryCard::new(
            Dimensionless::new(0.8),
            extreme,
            RecoveryReactor::Batch,
            context(),
            domain()
        )
        .is_err()
    );
}

#[test]
fn recovery_is_monotone_bounded_and_component_accounts_close() {
    let cards = [
        card(single(), RecoveryReactor::Batch),
        card(
            RateDistribution::Rectangular {
                maximum_rate: Frequency::new(0.02),
            },
            RecoveryReactor::Batch,
        ),
        card(
            RateDistribution::Gamma {
                shape: Dimensionless::new(0.3),
                scale: Frequency::new(0.1),
            },
            RecoveryReactor::PlugFlow,
        ),
        card(single(), RecoveryReactor::EqualCstrCascade { tanks: 6 }),
    ];
    for card in cards {
        let mut previous = 0.0;
        for time in [0.0, 1e-15, 0.1, 1.0, 10.0, 100.0, 1e4, 1e100] {
            let response = recovery(&card, time);
            let recovered = response.recovered.value();
            assert!(recovered >= previous && recovered <= 0.8);
            assert!(response.recoverable_remaining.value() >= 0.0);
            close(
                recovered
                    + response.recoverable_remaining.value()
                    + response.intrinsically_unrecoverable.value(),
                1.0,
                4e-15,
            );
            previous = recovered;
        }
        assert_eq!(recovery(&card, 0.0).recovered.value(), 0.0);
        close(recovery(&card, 1e100).recovered.value(), 0.8, 2e-15);
    }
    for ultimate in [0.0, 1.0] {
        let card = RecoveryCard::new(
            Dimensionless::new(ultimate),
            single(),
            RecoveryReactor::Batch,
            context(),
            domain(),
        )
        .unwrap();
        assert_eq!(
            recovery(&card, 0.0).intrinsically_unrecoverable.value(),
            1.0 - ultimate
        );
    }
}

#[test]
fn small_recoveries_and_long_time_tails_remain_representable() {
    let uniform = card(
        RateDistribution::Rectangular {
            maximum_rate: Frequency::new(1.0),
        },
        RecoveryReactor::Batch,
    );
    close(recovery(&uniform, 1e-20).recovered.value(), 0.4e-20, 2e-15);
    let gamma = card(
        RateDistribution::Gamma {
            shape: Dimensionless::new(2.0),
            scale: Frequency::new(1.0),
        },
        RecoveryReactor::Batch,
    );
    close(recovery(&gamma, 1e-20).recovered.value(), 1.6e-20, 2e-15);
    let many = card(
        RateDistribution::Single {
            rate: Frequency::new(1.0),
        },
        RecoveryReactor::EqualCstrCascade { tanks: u32::MAX },
    );
    close(recovery(&many, 1e-20).recovered.value(), 0.8e-20, 2e-15);
    let batch = card(single(), RecoveryReactor::Batch);
    let tail = recovery(&batch, 4_000.0);
    assert_eq!(tail.recovered.value(), 0.8); // Recovery rounded to its plateau.
    close(
        tail.recoverable_remaining.value(),
        0.8 * (-40.0_f64).exp(),
        4e-15,
    );
    assert!(tail.recoverable_remaining.value() > 0.0);
    let huge = card(
        RateDistribution::Rectangular {
            maximum_rate: Frequency::new(1e200),
        },
        RecoveryReactor::Batch,
    );
    close(
        recovery(&huge, 1e109).recoverable_remaining.value(),
        0.8e-309,
        1e-13,
    );
}

#[test]
fn gamma_preserves_extreme_products_and_fixed_rate_limit() {
    let narrow = card(
        RateDistribution::Gamma {
            shape: Dimensionless::new(1e9),
            scale: Frequency::new(1e-11),
        },
        RecoveryReactor::Batch,
    );
    let fixed = recovery(&card(single(), RecoveryReactor::Batch), 100.0)
        .recovered
        .value();
    assert!((recovery(&narrow, 100.0).recovered.value() - fixed).abs() < 3e-10);
    let tiny = card(
        RateDistribution::Gamma {
            shape: Dimensionless::new(1e300),
            scale: Frequency::new(1e-200),
        },
        RecoveryReactor::Batch,
    );
    close(recovery(&tiny, 1e-200).recovered.value(), 0.8e-100, 3e-15);
    let small_shape = card(
        RateDistribution::Gamma {
            shape: Dimensionless::new(1e-10),
            scale: Frequency::new(1e10),
        },
        RecoveryReactor::Batch,
    );
    // Multiplying the small shape by time first would erase this result.
    assert!(
        (recovery(&small_shape, 1e-320).recovered.value() - 0.8e-320).abs()
            <= 2.0 * f64::from_bits(1)
    );
    let huge = card(
        RateDistribution::Gamma {
            shape: Dimensionless::new(0.001),
            scale: Frequency::new(1e300),
        },
        RecoveryReactor::Batch,
    );
    close(
        recovery(&huge, 1e300).recovered.value(),
        0.8 * (1.0 - 10.0_f64.powf(-0.6)),
        4e-15,
    );
}

#[test]
fn source_conditions_require_coordinates_and_exact_quantity_conventions() {
    let length = QuantitySpec::dimensional(Length::DIMS);
    let support = domain()
        .with_quantity("p80_m", length, 70e-6, 200e-6)
        .with("pH", 8.0, 11.0);
    let card = RecoveryCard::new(
        Dimensionless::new(0.8),
        single(),
        RecoveryReactor::Batch,
        context(),
        support,
    )
    .unwrap();
    let conditions = QueryPoint::new()
        .with_quantity("p80_m", length, 100e-6)
        .unwrap()
        .with("pH", 9.5)
        .unwrap();
    assert!(card.evaluate(Time::new(100.0), &conditions).is_ok());
    assert!(matches!(
        card.evaluate(Time::new(100.0), &QueryPoint::new()),
        Err(RecoveryError::MissingCondition { .. })
    ));
    let untyped = QueryPoint::new()
        .with("p80_m", 100e-6)
        .unwrap()
        .with("pH", 9.5)
        .unwrap();
    assert!(matches!(
        card.evaluate(Time::new(100.0), &untyped),
        Err(RecoveryError::QuantityMismatch { .. })
    ));
    let outside = conditions.clone().with("p80_m", 201e-6).unwrap();
    assert!(matches!(
        card.evaluate(Time::new(100.0), &outside),
        Err(RecoveryError::OutsideValidity { .. })
    ));
    let duplicate = conditions.with(RESIDENCE_TIME_AXIS, 100.0).unwrap();
    assert!(card.evaluate(Time::new(100.0), &duplicate).is_err());
    let limited = RecoveryCard::new(
        Dimensionless::new(0.8),
        single(),
        RecoveryReactor::Batch,
        context(),
        domain().with(RESIDENCE_TIME_AXIS, 1.0, 100.0),
    )
    .unwrap();
    assert!(
        limited
            .evaluate(Time::new(0.0), &QueryPoint::new())
            .is_err()
    );
    assert!(
        limited
            .evaluate(Time::new(100.0_f64.next_up()), &QueryPoint::new())
            .is_err()
    );
    for time in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(card.evaluate(Time::new(time), &QueryPoint::new()).is_err());
    }
}

#[test]
fn invalid_parameters_missing_source_and_unsupported_reactors_refuse() {
    let make = |ultimate, rates, reactor, context, domain| {
        RecoveryCard::new(
            Dimensionless::new(ultimate),
            rates,
            reactor,
            context,
            domain,
        )
    };
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            make(
                invalid,
                single(),
                RecoveryReactor::Batch,
                context(),
                domain()
            )
            .is_err()
        );
    }
    assert!(
        make(
            1.0_f64.next_up(),
            single(),
            RecoveryReactor::Batch,
            context(),
            domain()
        )
        .is_err()
    );
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            make(
                0.8,
                RateDistribution::Single {
                    rate: Frequency::new(invalid)
                },
                RecoveryReactor::Batch,
                context(),
                domain()
            )
            .is_err()
        );
        assert!(
            make(
                0.8,
                RateDistribution::Gamma {
                    shape: Dimensionless::new(invalid),
                    scale: Frequency::new(1.0)
                },
                RecoveryReactor::Batch,
                context(),
                domain()
            )
            .is_err()
        );
    }
    for rates in [
        RateDistribution::Rectangular {
            maximum_rate: Frequency::new(0.01),
        },
        RateDistribution::Gamma {
            shape: Dimensionless::new(1.0),
            scale: Frequency::new(0.01),
        },
    ] {
        assert_eq!(
            make(
                0.8,
                rates,
                RecoveryReactor::EqualCstrCascade { tanks: 1 },
                context(),
                domain()
            ),
            Err(RecoveryError::UnsupportedReactorDistribution)
        );
    }
    assert!(
        make(
            0.8,
            single(),
            RecoveryReactor::EqualCstrCascade { tanks: 0 },
            context(),
            domain()
        )
        .is_err()
    );
    assert!(
        make(
            0.8,
            RateDistribution::Discrete {
                populations: vec![]
            },
            RecoveryReactor::Batch,
            context(),
            domain()
        )
        .is_err()
    );
    let zero_weight = RateDistribution::Discrete {
        populations: vec![RatePopulation {
            weight: Dimensionless::new(0.0),
            rate: Frequency::new(1.0),
        }],
    };
    assert!(
        make(
            0.8,
            zero_weight,
            RecoveryReactor::Batch,
            context(),
            domain()
        )
        .is_err()
    );
    assert!(
        make(
            0.8,
            single(),
            RecoveryReactor::Batch,
            context(),
            ValidityDomain::unconstrained()
        )
        .is_err()
    );
    assert!(
        make(
            0.8,
            single(),
            RecoveryReactor::Batch,
            context(),
            ValidityDomain::unconstrained().with(RESIDENCE_TIME_AXIS, 0.0, 100.0)
        )
        .is_err()
    );
    let mut missing_source = context();
    missing_source.provenance.source.clear();
    assert!(
        make(
            0.8,
            single(),
            RecoveryReactor::Batch,
            missing_source,
            domain()
        )
        .is_err()
    );
}
