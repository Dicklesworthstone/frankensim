//! G0/G1/G5 constitutive tests; all coefficients are synthetic benchmark data.
use fs_rand::StreamKey;
use fs_tribo::{
    ApplicabilityRange, DryFrictionApplicability, DryInterfaceSystemCard, FrictionLaw,
    FrictionRegime, HeatPartition, InputAuthority, InterfaceMedium, InterfaceSystemRef,
    TangentialLaw, TangentialSolveConfig, TangentialStepInput, TriboError, WorkLedger,
};

fn laws() -> [FrictionLaw; 3] {
    [
        FrictionLaw::Coulomb {
            static_mu: 0.8,
            kinetic_mu: 0.3,
        },
        FrictionLaw::VelocityDependent {
            static_mu: 0.8,
            mu_zero: 0.4,
            slope_per_speed: 0.1,
        },
        FrictionLaw::Stribeck {
            static_mu: 0.8,
            kinetic_mu: 0.3,
            characteristic_speed: 0.05,
            viscous_per_speed: 0.0,
        },
    ]
}

fn input() -> TangentialStepInput {
    TangentialStepInput {
        normal_n: 1.0,
        slip_increment_m: [0.03, 0.0],
        dt_s: 0.01,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: Some(2.0e5),
        compliance_m_per_n: 0.001,
        solve: TangentialSolveConfig {
            max_iterations: 80,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-14,
            force_relative_tolerance: 1e-14,
        },
    }
}

fn card() -> DryInterfaceSystemCard {
    DryInterfaceSystemCard::new(
        InterfaceSystemRef::new(
            "fixture/a->b",
            "fixture/dry-history",
            "fixture/tangential",
            InputAuthority::SyntheticFixture,
            InterfaceMedium::Dry,
        )
        .unwrap(),
        laws()[0].clone(),
        DryFrictionApplicability::new(
            ApplicabilityRange::new(280.0, 320.0).unwrap(),
            ApplicabilityRange::new(0.0, 5e5).unwrap(),
            ApplicabilityRange::new(0.0, 2.0).unwrap(),
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn g0_regularized_traction_continuous_at_regularization_speed() {
    let v_s = 0.05;
    for law in laws() {
        for sign in [-1.0, 1.0] {
            for normal in [0.1, 1.0, 100.0] {
                let left = law
                    .regularized_traction_1d(sign * v_s * (1.0 - 1e-12), normal, v_s)
                    .unwrap();
                let right = law
                    .regularized_traction_1d(sign * v_s * (1.0 + 1e-12), normal, v_s)
                    .unwrap();
                assert!(
                    (left - right).abs() <= 1e-12 * 0.8 * normal,
                    "{law:?}: {left} vs {right}"
                );
                let endpoint = law
                    .regularized_traction_1d(sign * v_s, normal, v_s)
                    .unwrap();
                let expected = -sign * law.kinetic_coefficient(v_s).unwrap() * normal;
                assert!((endpoint - expected).abs() <= f64::EPSILON * normal);
            }
        }
    }
}

#[test]
fn g0_step_work_closes_and_dissipation_nonnegative() {
    let mut stream = StreamKey {
        seed: 0xdef0_6000_0000_0001,
        kernel: 146,
        tile: 0,
    }
    .stream();
    for law in laws() {
        let mut state = law.initial_state();
        let mut ledger = WorkLedger::default();
        let partition = HeatPartition::new(0.5, 0.25, 0.25).unwrap();
        let mut sticking = 0;
        let mut sliding = 0;
        for _ in 0..10_000 {
            let x = TangentialStepInput {
                normal_n: 0.1 + 5.0 * stream.next_f64(),
                slip_increment_m: [
                    0.02 * (stream.next_f64() - 0.5),
                    0.02 * (stream.next_f64() - 0.5),
                ],
                dt_s: 0.01 + 0.01 * stream.next_f64(),
                compliance_m_per_n: 0.01 * stream.next_f64(),
                ..input()
            };
            let before = law.stored_energy_j(&state);
            let step = law.advance(&state, &x).unwrap();
            let work = step.traction_n[0] * step.slip_increment_m[0]
                + step.traction_n[1] * step.slip_increment_m[1];
            let delta = law.stored_energy_j(&step.state) - before;
            let tolerance = 1e-12 * work.abs().max(f64::EPSILON);
            assert!(step.dissipated_j >= 0.0);
            assert!((delta + step.dissipated_j - work).abs() <= tolerance);
            assert!((step.stored_change_j - delta).abs() <= tolerance);
            assert!(
                (step.closure_residual_j - (delta + step.dissipated_j - work)).abs() <= tolerance
            );
            for i in 0..2 {
                assert!(
                    (step.slip_increment_m[i]
                        - (x.slip_increment_m[i] - x.compliance_m_per_n * step.traction_n[i]))
                        .abs()
                        <= 2e-17
                );
            }
            match step.regime {
                FrictionRegime::Sticking => {
                    sticking += 1;
                    assert_eq!(step.slip_increment_m, [0.0; 2]);
                    assert!(
                        fs_math::det::hypot(step.traction_n[0], step.traction_n[1])
                            <= 0.8 * x.normal_n * (1.0 + 4.0 * f64::EPSILON)
                    );
                }
                FrictionRegime::Sliding => {
                    sliding += 1;
                    let slip =
                        fs_math::det::hypot(step.slip_increment_m[0], step.slip_increment_m[1]);
                    let force = fs_math::det::hypot(step.traction_n[0], step.traction_n[1]);
                    let expected = law.kinetic_coefficient(slip / x.dt_s).unwrap() * x.normal_n;
                    assert!((force - expected).abs() <= 4e-14 * (1.0 + expected));
                }
            }
            assert_eq!(step.state.contact_temperature_k(), x.contact_temperature_k);
            ledger
                .record(step.dissipation_step(partition).unwrap())
                .unwrap();
            state = step.state;
        }
        assert!(sticking > 0 && sliding > 0);
        assert!(ledger.dissipated_work_j() > 0.0);
        assert!(!law.model_id().is_empty());
    }
}

#[test]
fn g1_coulomb_trait_reproduces_block_on_incline() {
    let law = FrictionLaw::Coulomb {
        static_mu: 0.5,
        kinetic_mu: 0.4,
    };
    let (mass, gravity, dt) = (3.0, 9.806_65, 0.01);
    for theta in [0.2_f64, 0.7] {
        let normal = mass * gravity * theta.cos(); // det-ok: independent test oracle
        let free_acceleration = gravity * theta.sin(); // det-ok: independent test oracle
        // Exact held-force displacement of a free block initially at rest.
        let x = TangentialStepInput {
            normal_n: normal,
            dt_s: dt,
            slip_increment_m: [0.5 * free_acceleration * dt * dt, 0.0],
            compliance_m_per_n: 0.5 * dt * dt / mass,
            solve: TangentialSolveConfig {
                max_iterations: 1,
                ..input().solve
            },
            ..input()
        };
        let step = law.advance(&law.initial_state(), &x).unwrap();
        assert_eq!(
            step.iterations, 0,
            "Coulomb contact has an exact affine solution"
        );
        if theta.tan() < 0.5 {
            // det-ok: independent test oracle
            assert_eq!(step.regime, FrictionRegime::Sticking);
            assert_eq!(step.slip_increment_m, [0.0; 2]);
            assert!((step.traction_n[0] - mass * free_acceleration).abs() < 1e-12);
        } else {
            let expected_acceleration = gravity * (theta.sin() - 0.4 * theta.cos()); // det-ok: independent test oracle
            let acceleration = 2.0 * step.slip_increment_m[0] / (dt * dt);
            assert_eq!(step.regime, FrictionRegime::Sliding);
            assert!((acceleration - expected_acceleration).abs() <= 1e-12);
        }
    }
}

#[test]
fn g0_applicability_cap_and_cap_plus_one() {
    let law = card();
    let state = law.initial_state();
    let x = TangentialStepInput {
        dt_s: 1.0,
        slip_increment_m: [2.0, 0.0],
        compliance_m_per_n: 0.0,
        contact_temperature_k: Some(320.0),
        nominal_pressure_pa: Some(5e5),
        ..input()
    };
    assert!(law.advance(&state, &x).is_ok());
    for (field, bad) in [
        (
            "slip_speed_mps",
            TangentialStepInput {
                slip_increment_m: [2.0_f64.next_up(), 0.0],
                ..x
            },
        ),
        (
            "temperature_kelvin",
            TangentialStepInput {
                contact_temperature_k: Some(320.0_f64.next_up()),
                ..x
            },
        ),
        (
            "contact_pressure_pa",
            TangentialStepInput {
                nominal_pressure_pa: Some(5e5_f64.next_up()),
                ..x
            },
        ),
    ] {
        assert!(
            matches!(law.advance(&state, &bad), Err(TriboError::OutsideApplicability { field: got, .. }) if got == field)
        );
    }
    for normal in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            law.advance(
                &state,
                &TangentialStepInput {
                    normal_n: normal,
                    ..x
                }
            )
            .is_err()
        );
    }
    for bad in [
        TangentialStepInput {
            contact_temperature_k: None,
            ..x
        },
        TangentialStepInput {
            nominal_pressure_pa: None,
            ..x
        },
        TangentialStepInput { dt_s: 0.0, ..x },
        TangentialStepInput {
            compliance_m_per_n: -1.0,
            ..x
        },
        TangentialStepInput {
            slip_increment_m: [f64::NAN, 0.0],
            ..x
        },
    ] {
        assert!(law.advance(&state, &bad).is_err());
    }
    let mut at_cap = x;
    at_cap.solve.max_iterations = 128;
    assert!(law.advance(&state, &at_cap).is_ok());
    at_cap.solve.max_iterations = 129;
    assert!(matches!(
        law.advance(&state, &at_cap),
        Err(TriboError::InvalidInput {
            field: "max_iterations"
        })
    ));
    assert_eq!(state, law.initial_state());
}

#[test]
fn g0_force_budget_and_failed_root_refuse_without_publishing_state() {
    let law = laws()[0].clone();
    let state = law.initial_state();
    let mut x = TangentialStepInput {
        compliance_m_per_n: 0.0,
        ..input()
    };
    x.solve.maximum_force_n = 0.3;
    assert!(law.advance(&state, &x).is_ok());
    x.normal_n = 1.0_f64.next_up();
    assert!(law.advance(&state, &x).is_err());
    let mut underresolved = input();
    underresolved.solve.max_iterations = 1;
    assert!(matches!(
        laws()[2].advance(&state, &underresolved),
        Err(TriboError::TangentialSolve { iterations: 1, .. })
    ));
    let mut unbracketed = input();
    unbracketed.solve.maximum_force_n = 0.1;
    assert!(law.advance(&state, &unbracketed).is_err());
    for invalid in [
        FrictionLaw::Coulomb {
            static_mu: 0.1,
            kinetic_mu: 0.5,
        },
        FrictionLaw::Stribeck {
            static_mu: 0.5,
            kinetic_mu: 0.2,
            characteristic_speed: 0.0,
            viscous_per_speed: 0.0,
        },
    ] {
        assert!(invalid.advance(&state, &input()).is_err());
    }
    assert_eq!(state, law.initial_state());
}

/// G1: compare the affine fast path with an independent, tightly bisected
/// force/velocity intersection, including strengthening, weakening and the
/// zero-traction branch. No numerical stick regularization is introduced.
#[test]
fn g1_affine_sliding_matches_force_velocity_intersection_without_iterations() {
    for slope in [-2.0, -0.1, 0.0, 0.1, 20.0] {
        let law = FrictionLaw::VelocityDependent {
            static_mu: 0.8,
            mu_zero: 0.4,
            slope_per_speed: slope,
        };
        for normal in [0.1, 1.0, 10.0] {
            for direction in [[1.0, 0.0], [-1.0, 0.0], [0.6, 0.8]] {
                let x = TangentialStepInput {
                    normal_n: normal,
                    slip_increment_m: [0.03 * direction[0], 0.03 * direction[1]],
                    solve: TangentialSolveConfig {
                        max_iterations: 1,
                        ..input().solve
                    },
                    ..input()
                };
                let step = law.advance(&law.initial_state(), &x).unwrap();
                assert_eq!(step.regime, FrictionRegime::Sliding);
                assert_eq!(step.iterations, 0);
                let mut lo = 0.0;
                let mut hi = 0.03 / x.compliance_m_per_n;
                for _ in 0..160 {
                    let force = 0.5 * (lo + hi);
                    let velocity = (0.03 - x.compliance_m_per_n * force) / x.dt_s;
                    let kinetic = normal * (0.4 + slope * velocity).max(0.0);
                    if force < kinetic {
                        lo = force;
                    } else {
                        hi = force;
                    }
                }
                let reference = 0.5 * (lo + hi);
                let magnitude = fs_math::det::hypot(step.traction_n[0], step.traction_n[1]);
                assert!((magnitude - reference).abs() <= 1e-12 * (1.0 + reference));
                assert!(step.dissipated_j >= 0.0);
                assert_eq!(step.closure_residual_j, 0.0);
            }
        }
    }
}

#[test]
fn g0_static_capacity_endpoint_is_closed_but_zero_sliding_refuses() {
    let law = laws()[2].clone();
    let state = law.initial_state();
    let at_cap = TangentialStepInput {
        slip_increment_m: [0.8, 0.0],
        compliance_m_per_n: 1.0,
        ..input()
    };
    assert_eq!(
        law.advance(&state, &at_cap).unwrap().regime,
        FrictionRegime::Sticking
    );
    let over_cap = TangentialStepInput {
        slip_increment_m: [0.8_f64.next_up(), 0.0],
        ..at_cap
    };
    assert!(matches!(
        law.advance(&state, &over_cap),
        Err(TriboError::InvalidInput {
            field: "zero_sliding_increment"
        })
    ));
    assert_eq!(state, law.initial_state());
}

#[test]
fn g5_checkpoint_resume_bitwise() {
    let law = laws()[2].clone();
    let mut uninterrupted = law.initial_state();
    let mut resumed = law.initial_state();
    for i in 0..1000 {
        let x = TangentialStepInput {
            slip_increment_m: [if i % 2 == 0 { 0.01 } else { -0.02 }, 0.0],
            ..input()
        };
        let a = law.advance(&uninterrupted, &x).unwrap();
        let b = law.advance(&resumed, &x).unwrap();
        for (left, right) in a
            .traction_n
            .into_iter()
            .zip(b.traction_n)
            .chain(a.slip_increment_m.into_iter().zip(b.slip_increment_m))
            .chain(a.state.slip_m().into_iter().zip(b.state.slip_m()))
            .chain([
                (a.dissipated_j, b.dissipated_j),
                (a.stored_change_j, b.stored_change_j),
                (a.closure_residual_j, b.closure_residual_j),
            ])
        {
            assert_eq!(left.to_bits(), right.to_bits());
        }
        assert_eq!(a.regime, b.regime);
        assert_eq!(a.iterations, b.iterations);
        uninterrupted = a.state;
        resumed = if i == 499 { b.state.clone() } else { b.state };
    }
}
