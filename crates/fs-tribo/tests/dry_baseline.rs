//! G0/G1/G3 public conformance for the dependency-independent dry baseline.
//!
//! Every coefficient and thermal property here is synthetic. These tests check
//! equations, refusal, and replay; they do not validate a material system.

use fs_tribo::{
    ApplicabilityRange, ArchardLaw, ContactFrame, DryFrictionApplicability, DryFrictionState,
    DryInterfaceSystemCard, FlashTemperatureEstimate, FlashTemperatureInput,
    FlashTemperatureUnknown, FrictionLaw, FrictionRegime, HeatPartition, HertzCylinderPlane,
    HertzSpherePlane, InputAuthority, InterfaceMedium, InterfaceSystemRef,
    SurfaceThermalProperties, TangentialSlip, TriboError, WearState, WorkLedger,
    flash_temperature_candidate,
};

fn interface() -> InterfaceSystemRef {
    InterfaceSystemRef::new(
        "fixture/steel-coated-a->ceramic-b",
        "fixture/dry-history-v1",
        "fixture/synthetic-dry-card-v1",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .expect("synthetic dry identity")
}

fn slip(speed_mps: f64) -> TangentialSlip {
    let frame = ContactFrame::new([0.0, 0.0, 1.0]).expect("unit normal");
    TangentialSlip::new(&frame, [speed_mps, 0.0, 0.0]).expect("tangent velocity")
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1.0e-11 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g0_card_requires_ordered_domain_and_matching_slip_speed() {
    let card = DryInterfaceSystemCard::new(
        interface(),
        FrictionLaw::Coulomb {
            static_mu: 0.6,
            kinetic_mu: 0.4,
        },
        DryFrictionApplicability::new(
            ApplicabilityRange::new(280.0, 320.0).unwrap(),
            ApplicabilityRange::new(1.0e5, 5.0e6).unwrap(),
            ApplicabilityRange::new(0.0, 2.0).unwrap(),
        )
        .unwrap(),
    )
    .unwrap();
    let state = DryFrictionState::new(300.0, 2.0e6, 1.5).unwrap();
    let receipt = card.query(state, 10.0, slip(1.5)).unwrap();
    assert_eq!(receipt.response().regime, FrictionRegime::Sliding);
    assert_eq!(
        receipt.response().provenance().authority(),
        InputAuthority::SyntheticFixture
    );
    assert!(matches!(
        card.query(
            DryFrictionState::new(321.0, 2.0e6, 1.5).unwrap(),
            10.0,
            slip(1.5)
        ),
        Err(TriboError::OutsideApplicability {
            field: "temperature_kelvin",
            ..
        })
    ));
    assert!(matches!(
        card.query(
            DryFrictionState::new(300.0, 2.0e6, 1.0).unwrap(),
            10.0,
            slip(1.5)
        ),
        Err(TriboError::InvalidInput {
            field: "state.slip_speed_mps"
        })
    ));
}

#[test]
fn g1_hertz_and_block_on_incline_match_closed_forms() {
    let sphere = HertzSpherePlane {
        effective_radius: 0.02,
        reduced_modulus: 200.0e9,
    }
    .response(&interface(), 2.0e-6)
    .unwrap();
    close(sphere.contact_radius_m, 2.0e-4);
    close(sphere.normal_force_n, 106.666_666_666_666_67);

    let cylinder = HertzCylinderPlane {
        effective_radius: 0.02,
        reduced_modulus: 200.0e9,
    }
    .response(&interface(), 1_000.0)
    .unwrap();
    close(cylinder.half_width_m, 1.128_379_167_095_512_6e-5);
    close(cylinder.peak_pressure_pa, 56_418_958.354_775_63);

    let static_mu = 0.5;
    let law = FrictionLaw::Coulomb {
        static_mu,
        kinetic_mu: 0.4,
    };
    let mass_kg = 3.0;
    let gravity_m_per_s2 = 9.806_65;
    let theta = static_mu.atan(); // det-ok: test-only incline oracle
    let normal_force = mass_kg * gravity_m_per_s2 * theta.cos(); // det-ok: test-only incline oracle
    let downslope_force = mass_kg * gravity_m_per_s2 * theta.sin(); // det-ok: test-only incline oracle
    let response = law.evaluate(&interface(), normal_force, slip(0.0)).unwrap();
    close(response.static_limit, downslope_force);
}

#[test]
fn g3_reversal_scaling_and_replay_are_deterministic() {
    let law = FrictionLaw::Stribeck {
        static_mu: 0.7,
        kinetic_mu: 0.3,
        characteristic_speed: 1.0,
        viscous_per_speed: 0.02,
    };
    let forward = law.evaluate(&interface(), 10.0, slip(0.5)).unwrap();
    let reverse = law.evaluate(&interface(), 10.0, slip(-0.5)).unwrap();
    let scaled = law.evaluate(&interface(), 20.0, slip(0.5)).unwrap();
    close(reverse.traction_n()[0], -forward.traction_n()[0]);
    close(reverse.dissipated_power_w(), forward.dissipated_power_w());
    close(scaled.traction_n()[0], 2.0 * forward.traction_n()[0]);

    let law = ArchardLaw {
        wear_coefficient: 2.0e-6,
    };
    let mut first = WearState::default();
    let mut second = WearState::default();
    for (force, distance) in [(100.0, 0.2), (200.0, 0.1), (50.0, 0.4)] {
        law.advance(&interface(), &mut first, force, distance, 2.0e9)
            .unwrap();
        law.advance(&interface(), &mut second, force, distance, 2.0e9)
            .unwrap();
    }
    assert_eq!(first, second);
}

#[test]
fn g0_heat_partition_ledger_and_flash_candidate_are_closed_or_unknown() {
    let partition = HeatPartition::new(0.5, 0.25, 0.25).unwrap();
    let step = fs_tribo::DissipationStep::from_power(100.0, 0.2, partition).unwrap();
    let mut ledger = WorkLedger::default();
    ledger.record(step).unwrap();
    close(ledger.dissipated_work_j(), 20.0);
    close(ledger.surface_a_heat_j(), 10.0);
    close(ledger.surface_b_heat_j(), 5.0);
    close(ledger.other_work_j(), 5.0);

    let surface_a = SurfaceThermalProperties::new(100.0, 1.0e-5).unwrap();
    let surface_b = SurfaceThermalProperties::new(50.0, 4.0e-6).unwrap();
    let input = FlashTemperatureInput::new(
        100.0,
        partition,
        0.01,
        0.002,
        0.5,
        Some(surface_a),
        Some(surface_b),
    )
    .unwrap();
    let estimate = flash_temperature_candidate(&interface(), input).unwrap();
    assert!(matches!(estimate, FlashTemperatureEstimate::Candidate(_)));
    let FlashTemperatureEstimate::Candidate(candidate) = estimate else {
        return;
    };
    close(candidate.dwell_time_s(), 0.004);
    close(candidate.surface_a_heat_flux_w_per_m2(), 5_000.0);
    close(candidate.surface_b_heat_flux_w_per_m2(), 2_500.0);
    assert!(candidate.surface_a_rise_k() > 0.0);
    assert!(candidate.surface_b_rise_k() > 0.0);
    assert_eq!(
        candidate.provenance().authority(),
        InputAuthority::SyntheticFixture
    );

    let missing_b =
        FlashTemperatureInput::new(100.0, partition, 0.01, 0.002, 0.5, Some(surface_a), None)
            .unwrap();
    assert_eq!(
        flash_temperature_candidate(&interface(), missing_b).unwrap(),
        FlashTemperatureEstimate::Unknown(
            FlashTemperatureUnknown::MissingSurfaceBThermalProperties
        )
    );
}

#[test]
fn g3_flash_candidate_scales_with_power() {
    let partition = HeatPartition::new(0.5, 0.5, 0.0).unwrap();
    let thermal = Some(SurfaceThermalProperties::new(60.0, 1.2e-5).unwrap());
    let candidate = |power| {
        let input =
            FlashTemperatureInput::new(power, partition, 0.01, 0.004, 1.0, thermal, thermal)
                .unwrap();
        let estimate = flash_temperature_candidate(&interface(), input).unwrap();
        assert!(matches!(estimate, FlashTemperatureEstimate::Candidate(_)));
        let FlashTemperatureEstimate::Candidate(value) = estimate else {
            return None;
        };
        Some(value)
    };
    let one = candidate(25.0).unwrap();
    let two = candidate(50.0).unwrap();
    close(two.surface_a_rise_k(), 2.0 * one.surface_a_rise_k());
    close(two.surface_b_rise_k(), 2.0 * one.surface_b_rise_k());
}

#[test]
fn g1_effusivity_partition_equalizes_incremental_half_space_temperature_rise() {
    let a = SurfaceThermalProperties::new(12.0, 4e-6).unwrap();
    let b = SurfaceThermalProperties::new(6.0, 1e-6).unwrap();
    let equal = HeatPartition::semi_infinite_surfaces(a, b).unwrap();
    close(equal.surface_a(), 0.5);
    close(equal.surface_b(), 0.5);
    let more_conductive = SurfaceThermalProperties::new(24.0, 4e-6).unwrap();
    let partition = HeatPartition::semi_infinite_surfaces(more_conductive, b).unwrap();
    close(partition.surface_a(), 2.0 / 3.0);
    close(partition.surface_b(), 1.0 / 3.0);
    let input = FlashTemperatureInput::new(
        100.0,
        partition,
        0.01,
        0.004,
        1.0,
        Some(more_conductive),
        Some(b),
    )
    .unwrap();
    let FlashTemperatureEstimate::Candidate(candidate) =
        flash_temperature_candidate(&interface(), input).unwrap()
    else {
        panic!("complete properties");
    };
    close(candidate.surface_a_rise_k(), candidate.surface_b_rise_k());
    let reversed = HeatPartition::semi_infinite_surfaces(b, more_conductive).unwrap();
    close(reversed.surface_a(), partition.surface_b());
    close(reversed.surface_b(), partition.surface_a());
    assert_eq!(partition.other(), 0.0);
    // Stable normalization must not form an overflowing e_a + e_b.
    let huge = SurfaceThermalProperties::new(f64::MAX, 1.0).unwrap();
    assert_eq!(
        HeatPartition::semi_infinite_surfaces(huge, huge).unwrap(),
        equal
    );
}

fn slope_laws(scale: f64) -> [FrictionLaw; 5] {
    [
        FrictionLaw::Coulomb {
            static_mu: 0.9,
            kinetic_mu: 0.3,
        },
        FrictionLaw::VelocityDependent {
            static_mu: 0.9,
            mu_zero: 0.7,
            slope_per_speed: -0.1 / scale,
        },
        FrictionLaw::Stribeck {
            static_mu: 0.9,
            kinetic_mu: 0.3,
            characteristic_speed: scale,
            viscous_per_speed: 0.03 / scale,
        },
        FrictionLaw::Hyperbolic {
            static_mu: 0.9,
            kinetic_mu: 0.3,
            characteristic_speed: scale,
        },
        FrictionLaw::ShiftedPower {
            static_mu: 0.9,
            kinetic_mu: 0.3,
            characteristic_speed: scale,
            exponent: 0.4,
        },
    ]
}

#[test]
fn g1_kinetic_slope_has_physical_right_limits_and_closed_form_values() {
    let expected = [
        (0.3, 0.0),
        (0.7, -0.5),
        (0.9, 0.15),
        (0.9, -3.0),
        (0.9, -1.2),
    ];
    for (law, (expected_mu, expected_slope)) in slope_laws(0.2).into_iter().zip(expected) {
        for zero in [0.0, -0.0] {
            let (mu, slope) = law.kinetic_coefficient_and_slope(zero).unwrap();
            assert_eq!(
                mu.to_bits(),
                law.zero_speed_kinetic_coefficient().unwrap().to_bits()
            );
            close(mu, expected_mu);
            close(slope, expected_slope);
        }
    }
    // det-ok: independent test-only exponential/power closed forms.
    let decay = (-1.0_f64).exp();
    let laws = slope_laws(0.2);
    let (mu, slope) = laws[2].kinetic_coefficient_and_slope(0.2).unwrap();
    close(mu, 0.3 + 0.6 * decay + 0.03);
    close(slope, 0.15 - 6.0 * decay);
    let (mu, slope) = laws[3].kinetic_coefficient_and_slope(0.2).unwrap();
    close(mu, 0.6);
    close(slope, -0.75);
    let (mu, slope) = laws[4].kinetic_coefficient_and_slope(0.2).unwrap();
    close(mu, 0.3 + 0.6 * 2.0_f64.powf(-0.4));
    close(slope, -1.2 * 2.0_f64.powf(-1.4));

    let clamped = FrictionLaw::VelocityDependent {
        static_mu: 0.9,
        mu_zero: 0.75,
        slope_per_speed: -0.25,
    };
    assert_eq!(
        clamped.kinetic_coefficient_and_slope(3.0).unwrap(),
        (0.0, 0.0)
    );
    assert_eq!(
        clamped.kinetic_coefficient_and_slope(4.0).unwrap(),
        (0.0, 0.0)
    );
    assert_eq!(
        clamped
            .kinetic_coefficient_and_slope(3.0_f64.next_down())
            .unwrap()
            .1,
        -0.25
    );
    for (growth, expected) in [(-0.25, 0.0), (0.25, 0.25)] {
        let law = FrictionLaw::VelocityDependent {
            static_mu: 0.9,
            mu_zero: 0.0,
            slope_per_speed: growth,
        };
        assert_eq!(
            law.kinetic_coefficient_and_slope(0.0).unwrap(),
            (0.0, expected)
        );
    }
}

#[test]
fn g3_kinetic_slope_matches_independent_differences_and_existing_value_bits() {
    for scale in [1e-12, 1e-4, 1.0, 1e4, 1e12] {
        for law in slope_laws(scale) {
            for ratio in [0.01, 0.3, 1.0, 3.0, 8.0] {
                let speed = scale * ratio;
                let (mu, slope) = law.kinetic_coefficient_and_slope(speed).unwrap();
                assert_eq!(
                    mu.to_bits(),
                    law.kinetic_coefficient(speed).unwrap().to_bits()
                );
                let step = speed * 1e-5;
                let difference = (law.kinetic_coefficient(speed + step).unwrap()
                    - law.kinetic_coefficient(speed - step).unwrap())
                    / (2.0 * step);
                assert!(
                    ((slope - difference) * scale).abs() < 2e-8,
                    "law={law:?}, speed={speed:e}, slope={slope:e}, oracle={difference:e}"
                );
            }
        }
        let hyperbolic = FrictionLaw::Hyperbolic {
            static_mu: 0.9,
            kinetic_mu: 0.3,
            characteristic_speed: scale,
        };
        let power = FrictionLaw::ShiftedPower {
            static_mu: 0.9,
            kinetic_mu: 0.3,
            characteristic_speed: scale,
            exponent: 1.0,
        };
        for ratio in [0.0, 0.01, 1.0, 10.0, 1e6] {
            let (hm, hs) = hyperbolic
                .kinetic_coefficient_and_slope(ratio * scale)
                .unwrap();
            let (pm, ps) = power.kinetic_coefficient_and_slope(ratio * scale).unwrap();
            close(hm, pm);
            assert!((hs - ps).abs() <= 8.0 * f64::EPSILON * hs.abs());
        }
    }
}

#[test]
fn g3_kinetic_slope_retains_rounded_away_tails_and_dimensional_scaling() {
    // det-ok: test-only independent closed forms for analytic tail slopes.
    let tail_cases = [
        (
            FrictionLaw::Stribeck {
                static_mu: 0.8,
                kinetic_mu: 0.3,
                characteristic_speed: 0.05,
                viscous_per_speed: 0.0,
            },
            0.4,
            -160.0 * (-64.0_f64).exp(),
        ),
        (
            FrictionLaw::Hyperbolic {
                static_mu: 0.8,
                kinetic_mu: 0.3,
                characteristic_speed: 0.05,
            },
            5e15,
            -1e-33,
        ),
        (
            FrictionLaw::ShiftedPower {
                static_mu: 0.8,
                kinetic_mu: 0.3,
                characteristic_speed: 0.05,
                exponent: 3.0,
            },
            5e4,
            -30.0 / 1_000_001.0_f64.powi(4),
        ),
    ];
    for (law, speed, expected) in tail_cases {
        let (mu, slope) = law.kinetic_coefficient_and_slope(speed).unwrap();
        assert_eq!(mu, 0.3);
        assert!(slope < 0.0);
        assert!(
            (slope - expected).abs() <= 2e-13 * expected.abs(),
            "{law:?}: {slope:e} vs {expected:e}"
        );
    }
    // Scaling Delta_mu, v and v0 together leaves the derivative unchanged,
    // even when dimensional products or v0^2 would overflow or underflow.
    for scale in [f64::from_bits(1), 1e-300, 1.0, 1e300] {
        for (law, expected) in [
            (
                FrictionLaw::Hyperbolic {
                    static_mu: scale,
                    kinetic_mu: 0.0,
                    characteristic_speed: scale,
                },
                -0.25,
            ),
            (
                FrictionLaw::Stribeck {
                    static_mu: scale,
                    kinetic_mu: 0.0,
                    characteristic_speed: scale,
                    viscous_per_speed: 0.0,
                },
                -2.0 * (-1.0_f64).exp(),
            ),
            (
                FrictionLaw::ShiftedPower {
                    static_mu: scale,
                    kinetic_mu: 0.0,
                    characteristic_speed: scale,
                    exponent: 2.0,
                },
                -0.25,
            ),
        ] {
            let (mu, slope) = law.kinetic_coefficient_and_slope(scale).unwrap();
            assert_eq!(
                mu.to_bits(),
                law.kinetic_coefficient(scale).unwrap().to_bits()
            );
            close(slope, expected);
        }
    }
}

#[test]
fn g0_kinetic_slope_refuses_invalid_inputs_and_unresolved_derivatives() {
    for law in slope_laws(0.2) {
        for speed in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(law.kinetic_coefficient_and_slope(speed).is_err());
        }
    }
    let invalid = FrictionLaw::Stribeck {
        static_mu: 0.3,
        kinetic_mu: 0.8,
        characteristic_speed: 0.2,
        viscous_per_speed: 0.0,
    };
    assert!(invalid.kinetic_coefficient_and_slope(0.0).is_err());
    for (law, speed) in [
        (
            FrictionLaw::Hyperbolic {
                static_mu: 1.0,
                kinetic_mu: 0.0,
                characteristic_speed: f64::from_bits(1),
            },
            f64::from_bits(1),
        ),
        (
            FrictionLaw::Stribeck {
                static_mu: 0.8,
                kinetic_mu: 0.3,
                characteristic_speed: 1e-300,
                viscous_per_speed: 0.0,
            },
            3e-299,
        ),
        (
            FrictionLaw::ShiftedPower {
                static_mu: 0.8,
                kinetic_mu: 0.3,
                characteristic_speed: 1e-300,
                exponent: 2000.0,
            },
            1e-300,
        ),
    ] {
        // Value-only consumers retain their original admission. A derivative
        // is evaluated only when requested; lost decay cannot be a zero slope.
        assert!(law.kinetic_coefficient(speed).unwrap().is_finite());
        assert!(law.kinetic_coefficient_and_slope(speed).is_err());
    }
    let far_tail = FrictionLaw::Stribeck {
        static_mu: 0.8,
        kinetic_mu: 0.3,
        characteristic_speed: 0.1,
        viscous_per_speed: 0.0,
    };
    assert_eq!(
        far_tail.kinetic_coefficient_and_slope(10.0).unwrap(),
        (0.3, 0.0)
    );
}
