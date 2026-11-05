use super::*;

fn budget() -> KineticsBudget {
    KineticsBudget {
        max_delayed_groups: 8,
        max_squarings: 40,
        max_taylor_terms: 24,
        population_tolerance: 1.0e-7,
        integral_tolerance_population_s: 1.0e-6,
    }
}

fn model(generation: f64, groups: &[(f64, f64)]) -> Arc<PointKineticsModel> {
    Arc::new(
        PointKineticsModel::new(
            "synthetic-equation-verification-v1".into(),
            generation,
            groups
                .iter()
                .map(|&(fraction, decay)| DelayedNeutronGroup {
                    delayed_fraction: fraction,
                    decay_constant_s_inv: decay,
                })
                .collect(),
        )
        .unwrap(),
    )
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

#[test]
fn g1_prompt_population_and_integral_match_affine_exponential() {
    for reactivity in [-0.4_f64, 0.0, 0.2] {
        let mut state = PointKineticsState::new(model(2.0, &[]), 3.0, vec![]).unwrap();
        let duration = 4.0;
        let source = 0.7;
        let rate = reactivity / 2.0;
        let (expected_n, expected_integral) = if rate == 0.0 {
            (
                3.0 + source * duration,
                3.0 * duration + 0.5 * source * duration * duration,
            )
        } else {
            let expm1 = (rate * duration).exp_m1();
            (
                3.0 * (rate * duration).exp() + source * expm1 / rate,
                3.0 * expm1 / rate + source * (expm1 / rate - duration) / rate,
            )
        };
        let step = state
            .advance(duration, reactivity, source, budget(), || Ok(()))
            .unwrap();
        near(state.neutron_population(), expected_n, 2.0e-13);
        near(step.integrated_population_s, expected_integral, 5.0e-13);
        near(step.external_neutrons, 2.8, 1.0e-15);
        assert!(step.population_balance_residual.abs() < 1.0e-12);
    }
}

#[test]
fn g1_one_group_matches_independent_two_eigenmode_solution() {
    let generation = 0.02;
    let beta = 0.1;
    let decay = 0.3;
    let n0 = 2.0;
    let c0 = 0.7;
    let duration = 2.0;
    for reactivity in [-0.02_f64, 0.01, 0.12] {
        let a = (reactivity - beta) / generation;
        let b = beta / generation;
        let separation = ((a + decay) * (a + decay) + 4.0 * decay * b).sqrt();
        let slow = 0.5 * (a - decay + separation);
        let fast = 0.5 * (a - decay - separation);
        let slow_coefficient = (a * n0 + decay * c0 - fast * n0) / separation;
        let fast_coefficient = n0 - slow_coefficient;
        let slow_n = slow_coefficient * (slow * duration).exp();
        let fast_n = fast_coefficient * (fast * duration).exp();
        let expected_n = slow_n + fast_n;
        let expected_c = (slow * slow_n + fast * fast_n - a * expected_n) / decay;
        let expected_integral = slow_coefficient * (slow * duration).exp_m1() / slow
            + fast_coefficient * (fast * duration).exp_m1() / fast;
        let mut state =
            PointKineticsState::new(model(generation, &[(beta, decay)]), n0, vec![c0]).unwrap();
        let step = state
            .advance(duration, reactivity, 0.0, budget(), || Ok(()))
            .unwrap();
        near(state.neutron_population(), expected_n, 5.0e-12);
        near(state.precursor_populations()[0], expected_c, 3.0e-11);
        near(step.integrated_population_s, expected_integral, 8.0e-12);
        assert!(state.neutron_population() >= 0.0 && state.precursor_populations()[0] >= 0.0);
    }
}

#[test]
fn g1_multigroup_critical_equilibrium_and_population_time_are_preserved() {
    let data = model(
        0.002,
        &[(0.001, 0.03), (0.002, 0.2), (0.001, 0.9), (0.0005, 3.0)],
    );
    let precursors = data.equilibrium_precursors(2.0).unwrap();
    let mut state = PointKineticsState::new(data, 2.0, precursors.clone()).unwrap();
    let step = state.advance(40.0, 0.0, 0.0, budget(), || Ok(())).unwrap();
    near(state.neutron_population(), 2.0, 2.0e-11);
    for (&actual, &expected) in state.precursor_populations().iter().zip(&precursors) {
        near(actual, expected, 3.0e-10);
    }
    near(step.integrated_population_s, 80.0, 1.0e-9);
    assert!(step.population_balance_residual.abs() < 1.0e-9);
}

#[test]
fn g1_stiff_prompt_relaxation_uses_logarithmic_scaling() {
    let generation = 1.0e-7;
    let beta = 0.1;
    let decay = 1.0;
    let rate = beta / generation + decay;
    let equilibrium_n = decay / rate;
    let mut state =
        PointKineticsState::new(model(generation, &[(beta, decay)]), 1.0, vec![0.0]).unwrap();
    let tolerance = KineticsBudget {
        population_tolerance: 1.0e-5,
        integral_tolerance_population_s: 1.0e-5,
        ..budget()
    };
    let step = state.advance(1.0, 0.0, 0.0, tolerance, || Ok(())).unwrap();
    let expected_n = equilibrium_n + (1.0 - equilibrium_n) * (-rate).exp();
    let expected_integral = equilibrium_n + (1.0 - equilibrium_n) * -(-rate).exp_m1() / rate;
    near(state.neutron_population(), expected_n, 2.0e-14);
    near(state.precursor_populations()[0], 1.0 - expected_n, 2.0e-8);
    near(step.integrated_population_s, expected_integral, 3.0e-14);
    assert!(step.squarings < 30 && step.taylor_terms <= 24);
}

#[test]
fn g1_source_builds_precursors_and_shutdown_retains_them() {
    let data = model(0.1, &[(0.1, 0.2)]);
    let mut state = PointKineticsState::new(data, 0.0, vec![0.0]).unwrap();
    let irradiation = state.advance(4.0, 0.0, 3.0, budget(), || Ok(())).unwrap();
    near(
        state.neutron_population() + state.precursor_populations()[0],
        12.0,
        3.0e-12,
    );
    assert!(irradiation.integrated_population_s > 0.0);
    let precursor_history = state.precursor_populations()[0];
    let cooled = state.advance(2.0, -0.4, 0.0, budget(), || Ok(())).unwrap();
    assert!(state.precursor_populations()[0] > 0.0);
    assert!(state.precursor_populations()[0] < precursor_history);
    assert!(state.neutron_population() > 0.0 && cooled.integrated_population_s > 0.0);
    near(cooled.population_balance_residual, 0.0, 1.0e-11);
}

#[test]
fn g3_partitioned_integrals_and_g5_checkpoint_replay_agree() {
    let data = model(0.07, &[(0.03, 0.2), (0.02, 0.9)]);
    let start = PointKineticsState::new(data, 1.7, vec![0.6, 0.2]).unwrap();
    let (whole, report) = start.trial(3.0, 0.01, 0.3, budget(), || Ok(())).unwrap();
    let replay = start.trial(3.0, 0.01, 0.3, budget(), || Ok(())).unwrap();
    assert_eq!((whole.clone(), report), replay);
    let mut partitioned = start;
    let mut integrated = 0.0;
    for _ in 0..30 {
        integrated += partitioned
            .advance(0.1, 0.01, 0.3, budget(), || Ok(()))
            .unwrap()
            .integrated_population_s;
    }
    near(
        partitioned.neutron_population(),
        whole.neutron_population(),
        2.0e-12,
    );
    for (&actual, &expected) in partitioned
        .precursor_populations()
        .iter()
        .zip(whole.precursor_populations())
    {
        near(actual, expected, 2.0e-12);
    }
    near(integrated, report.integrated_population_s, 4.0e-12);
}

#[test]
fn g0_zero_duration_and_vacuum_are_exact_but_validate_inputs() {
    let start = PointKineticsState::new(model(0.1, &[(0.2, 0.3)]), 2.0, vec![4.0]).unwrap();
    let no_work = KineticsBudget {
        max_squarings: 0,
        max_taylor_terms: 0,
        population_tolerance: 0.0,
        integral_tolerance_population_s: 0.0,
        ..budget()
    };
    let (unchanged, report) = start.trial(0.0, 0.1, 1.0, no_work, || Ok(())).unwrap();
    assert_eq!(unchanged, start);
    assert_eq!(report, KineticsStep::default());
    assert!(start.trial(0.0, f64::NAN, 0.0, no_work, || Ok(())).is_err());
    let vacuum = PointKineticsState::new(model(0.1, &[(0.2, 0.3)]), 0.0, vec![0.0]).unwrap();
    let (unchanged, report) = vacuum.trial(10.0, 0.5, 0.0, no_work, || Ok(())).unwrap();
    assert_eq!(unchanged, vacuum);
    assert_eq!(report.integrated_population_s, 0.0);
    assert_eq!(report.squarings, 0);
}

#[test]
fn g1_tiny_interval_retains_positive_population_time() {
    let mut state = PointKineticsState::new(model(0.001, &[(0.01, 0.1)]), 1.0, vec![0.0]).unwrap();
    let report = state
        .advance(1.0e-14, -0.02, 0.0, budget(), || Ok(()))
        .unwrap();
    assert!(report.integrated_population_s > 0.0);
    near(report.integrated_population_s, 1.0e-14, 2.0e-27);
    assert!(state.neutron_population() <= 1.0);
    assert!(state.precursor_populations()[0] > 0.0);
}

#[test]
fn g4_cancellation_work_accuracy_and_overflow_refuse_atomically() {
    let start = PointKineticsState::new(model(0.1, &[(0.2, 0.3)]), 2.0, vec![4.0]).unwrap();
    let mut state = start.clone();
    let mut polls = 0;
    let cancellation = state.advance(4.0, 0.1, 0.4, budget(), || {
        polls += 1;
        if polls == 30 {
            Err(KineticsError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert_eq!(cancellation, Err(KineticsError::Cancelled));
    assert_eq!(state, start);
    for limited in [
        KineticsBudget {
            max_delayed_groups: 0,
            ..budget()
        },
        KineticsBudget {
            max_squarings: 0,
            ..budget()
        },
        KineticsBudget {
            max_taylor_terms: 1,
            ..budget()
        },
        KineticsBudget {
            population_tolerance: 0.0,
            ..budget()
        },
        KineticsBudget {
            integral_tolerance_population_s: 0.0,
            ..budget()
        },
    ] {
        assert!(matches!(
            state.advance(4.0, 0.1, 0.4, limited, || Ok(())),
            Err(KineticsError::Budget(_))
        ));
        assert_eq!(state, start);
    }
    let mut exploding = PointKineticsState::new(model(1.0e-6, &[]), 1.0, vec![]).unwrap();
    let before = exploding.clone();
    assert_eq!(
        exploding.advance(0.01, 0.5, 0.0, budget(), || Ok(())),
        Err(KineticsError::Arithmetic)
    );
    assert_eq!(exploding, before);
}

#[test]
fn g0_invalid_data_and_mismatched_precursor_history_are_rejected() {
    assert!(PointKineticsModel::new(" ".into(), 1.0, vec![]).is_err());
    assert!(PointKineticsModel::new("x".into(), 0.0, vec![]).is_err());
    assert!(
        PointKineticsModel::new(
            "x".into(),
            1.0,
            vec![DelayedNeutronGroup {
                delayed_fraction: 1.0,
                decay_constant_s_inv: 0.2,
            }]
        )
        .is_err()
    );
    assert!(
        PointKineticsModel::new(
            "x".into(),
            1.0,
            vec![DelayedNeutronGroup {
                delayed_fraction: -0.1,
                decay_constant_s_inv: 0.2,
            }]
        )
        .is_err()
    );
    let data = model(0.1, &[(0.2, 0.3)]);
    assert!(PointKineticsState::new(Arc::clone(&data), 1.0, vec![]).is_err());
    assert!(PointKineticsState::new(data, 1.0, vec![f64::NAN]).is_err());
}

#[test]
fn g1_explicit_linear_and_log_temperature_feedback_obey_domains() {
    let linear = TemperatureFeedback::new(
        "synthetic-linear-fit".into(),
        300.0,
        250.0,
        1200.0,
        TemperatureFeedbackLaw::Linear {
            coefficient_per_kelvin: -2.0e-5,
        },
    )
    .unwrap();
    let logarithmic = TemperatureFeedback::new(
        "synthetic-doppler-fit".into(),
        600.0,
        300.0,
        1500.0,
        TemperatureFeedbackLaw::Logarithmic {
            coefficient: -0.008,
        },
    )
    .unwrap();
    assert_eq!(linear.reactivity_shift(300.0).unwrap(), 0.0);
    near(linear.reactivity_shift(800.0).unwrap(), -0.01, 1.0e-17);
    near(
        logarithmic.reactivity_shift(1200.0).unwrap(),
        -0.008 * 2.0_f64.ln(),
        1.0e-17,
    );
    assert!(linear.reactivity_shift(1201.0).is_err());
    assert!(logarithmic.reactivity_shift(0.0).is_err());
    assert!(logarithmic.reactivity_shift(f64::NAN).is_err());
    let data = model(0.1, &[(0.01, 0.2)]);
    let equilibrium = data.equilibrium_precursors(1.0).unwrap();
    let start = PointKineticsState::new(data, 1.0, equilibrium).unwrap();
    let cool = start
        .trial(
            1.0,
            0.005 + linear.reactivity_shift(300.0).unwrap(),
            0.0,
            budget(),
            || Ok(()),
        )
        .unwrap();
    let hot = start
        .trial(
            1.0,
            0.005 + linear.reactivity_shift(800.0).unwrap(),
            0.0,
            budget(),
            || Ok(()),
        )
        .unwrap();
    assert!(cool.0.neutron_population() > 1.0);
    assert!(hot.0.neutron_population() < 1.0);
    assert!(hot.1.integrated_population_s < cool.1.integrated_population_s);
}
