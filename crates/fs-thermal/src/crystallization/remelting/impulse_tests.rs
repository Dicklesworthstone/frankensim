use super::*;

fn config() -> CrystallizationConfig {
    CrystallizationConfig {
        identity: [0x73; 32],
        avrami_exponent: 2.0,
        maximum_crystalline_mass_fraction: 0.6,
        temperature_range_k: [250.0, 650.0],
        melting_temperature_k: 500.0,
    }
}

fn parent(t: f64) -> Result<f64, CrystallizationError> {
    Ok(2000.0 * (t - 300.0) + 100_000.0)
}

fn crystal(t: f64) -> Result<f64, CrystallizationError> {
    Ok(1000.0 * (t - 300.0))
}

fn rate(t: f64) -> Result<f64, CrystallizationError> {
    Ok(0.001 * (500.0 - t).max(0.0))
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_instantaneous_sensible_changes_freeze_history_even_with_positive_kinetics() {
    let queries = std::cell::Cell::new(0);
    let kinetic = |t| {
        queries.set(queries.get() + 1);
        rate(t)
    };
    let law = CrystallizationLaw::new(&parent, &crystal, &kinetic, config()).unwrap();
    let remelting = law.equivalent_age_remelting().unwrap();
    let initial = law.state_at_temperature(450.0, 0.7).unwrap();
    let solve = CrystallizationSolve::default();
    assert!(law.rate_s_inv(300.0).unwrap() > 0.0);
    for model in [None, Some(remelting)] {
        for temperature in [300.0, 490.0] {
            let h = law
                .specific_enthalpy_j_kg(temperature, initial.crystalline_mass_fraction())
                .unwrap();
            queries.set(0);
            let transition = law
                .apply_specific_enthalpy(initial, h, solve, model, || false)
                .unwrap();
            // No temperature-path kinetics are evaluated. Selecting remelting
            // only rechecks its prepared rate boundary at Tm.
            assert_eq!(queries.get(), usize::from(model.is_some()));
            assert_eq!(transition.initial, initial);
            assert_eq!(transition.state.exposure(), initial.exposure());
            assert_eq!(
                transition.state.crystalline_mass_fraction(),
                initial.crystalline_mass_fraction()
            );
            assert_eq!(transition.state.specific_enthalpy_j_kg(), h);
            assert_eq!(
                transition.specific_heat_j_kg,
                h - initial.specific_enthalpy_j_kg()
            );
            near(
                transition.state.temperature_k(),
                temperature,
                solve.temperature_tolerance_k,
            );
            law.validate_state(transition.state, solve).unwrap();
        }
    }
    let hot = law.state_at_temperature(550.0, 0.0).unwrap();
    let cooled = law
        .apply_specific_enthalpy(hot, parent(300.0).unwrap(), solve, Some(remelting), || {
            false
        })
        .unwrap();
    assert_eq!(cooled.state.exposure(), 0.0);
    assert_eq!(cooled.state.crystalline_mass_fraction(), 0.0);
    near(
        cooled.state.temperature_k(),
        300.0,
        solve.temperature_tolerance_k,
    );
}

#[test]
fn g1_instantaneous_partial_full_melting_and_superheat_use_one_original_account() {
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let model = law.equivalent_age_remelting().unwrap();
    let solve = CrystallizationSolve::default();
    let initial = law.state_at_temperature(500.0, 1.0).unwrap();
    let half_fraction = 0.5 * initial.crystalline_mass_fraction();
    let partial_h = 500_000.0 - 300_000.0 * half_fraction;
    assert!(matches!(
        law.apply_specific_enthalpy(initial, partial_h, solve, None, || false),
        Err(CrystallizationError::EnthalpyOutOfRange { .. })
    ));
    let partial = law
        .apply_specific_enthalpy(initial, partial_h, solve, Some(model), || false)
        .unwrap();
    assert_eq!(partial.state.temperature_k(), 500.0);
    near(
        partial.state.crystalline_mass_fraction(),
        half_fraction,
        5e-14,
    );
    near(
        partial.state.exposure(),
        (-(1.0 - half_fraction / 0.6).ln()).sqrt(),
        5e-13,
    );
    near(partial.specific_heat_j_kg, 300_000.0 * half_fraction, 1e-9);
    assert!(partial.state.exposure() < initial.exposure());
    let full = law
        .apply_specific_enthalpy(partial.state, 500_000.0, solve, Some(model), || false)
        .unwrap();
    near(
        full.state.temperature_k(),
        500.0,
        solve.temperature_tolerance_k,
    );
    assert_eq!(full.state.exposure(), 0.0);
    assert_eq!(full.state.crystalline_mass_fraction(), 0.0);
    let superheated = law
        .apply_specific_enthalpy(full.state, 600_000.0, solve, Some(model), || false)
        .unwrap();
    near(
        superheated.state.temperature_k(),
        550.0,
        solve.temperature_tolerance_k,
    );
    assert_eq!(superheated.state.exposure(), 0.0);
    assert_eq!(superheated.state.specific_enthalpy_j_kg(), 600_000.0);
    near(
        partial.specific_heat_j_kg + full.specific_heat_j_kg + superheated.specific_heat_j_kg,
        600_000.0 - initial.specific_enthalpy_j_kg(),
        1e-9,
    );
    for transition in [partial, full, superheated] {
        assert!(transition.caloric_residual_j_kg.abs() <= solve.specific_enthalpy_tolerance_j_kg);
        law.validate_state(transition.state, solve).unwrap();
    }
}

#[test]
fn g0_zero_heat_preserves_the_accepted_residual_and_positive_rate_state_bitwise() {
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let original = law.state_at_temperature(350.0, 0.4).unwrap();
    let rounded = CrystallizationState {
        specific_enthalpy_j_kg: f64::from_bits(original.specific_enthalpy_j_kg().to_bits() + 1),
        ..original
    };
    let solve = CrystallizationSolve::default();
    law.validate_state(rounded, solve).unwrap();
    assert!(law.rate_s_inv(rounded.temperature_k()).unwrap() > 0.0);
    let model = law.equivalent_age_remelting().unwrap();
    for remelting in [None, Some(model)] {
        let transition = law
            .apply_specific_enthalpy(
                rounded,
                rounded.specific_enthalpy_j_kg(),
                solve,
                remelting,
                || false,
            )
            .unwrap();
        assert_eq!(transition.initial, rounded);
        assert_eq!(transition.state, rounded);
        assert_eq!(transition.specific_heat_j_kg, 0.0);
        assert_eq!(
            transition.caloric_residual_j_kg,
            original.specific_enthalpy_j_kg() - rounded.specific_enthalpy_j_kg()
        );
    }
}

#[test]
fn g3_impulse_melting_retains_large_datums_and_saturated_finite_history() {
    let solve = CrystallizationSolve {
        specific_enthalpy_tolerance_j_kg: 0.001,
        temperature_tolerance_k: 1e-6,
        ..CrystallizationSolve::default()
    };
    let run = |offset| {
        let a = |t| parent(t).map(|h| h + offset);
        let c = |t| crystal(t).map(|h| h + offset);
        let law = CrystallizationLaw::new(&a, &c, &rate, config()).unwrap();
        let initial = law.state_at_temperature(500.0, f64::MAX).unwrap();
        let target = offset + 410_000.0;
        let transition = law
            .apply_specific_enthalpy(
                initial,
                target,
                solve,
                Some(law.equivalent_age_remelting().unwrap()),
                || false,
            )
            .unwrap();
        assert_eq!(transition.state.temperature_k(), 500.0);
        assert_eq!(transition.state.specific_enthalpy_j_kg(), target);
        near(transition.state.crystalline_mass_fraction(), 0.3, 4e-9);
        near(transition.specific_heat_j_kg, 90_000.0, 0.001);
        assert!(transition.state.exposure().is_finite() && transition.state.exposure() < f64::MAX);
        law.validate_state(transition.state, solve).unwrap();
        transition.state
    };
    let original = run(0.0);
    let shifted = run(1e12);
    near(
        original.crystalline_mass_fraction(),
        shifted.crystalline_mass_fraction(),
        4e-9,
    );
    near(original.exposure(), shifted.exposure(), 1e-8);
}

#[test]
fn g0_impulse_controls_domain_and_prepared_source_mismatch_refuse() {
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let initial = law.state_at_temperature(350.0, 0.4).unwrap();
    let solve = CrystallizationSolve::default();
    solve.validate().unwrap();
    for invalid in [
        CrystallizationSolve {
            maximum_iterations: 0,
            ..solve
        },
        CrystallizationSolve {
            maximum_iterations: 257,
            ..solve
        },
        CrystallizationSolve {
            specific_enthalpy_tolerance_j_kg: 0.0,
            ..solve
        },
        CrystallizationSolve {
            temperature_tolerance_k: f64::NAN,
            ..solve
        },
    ] {
        assert!(matches!(
            invalid.validate(),
            Err(CrystallizationError::Invalid(_))
        ));
        assert!(matches!(
            law.apply_specific_enthalpy(
                initial,
                initial.specific_enthalpy_j_kg(),
                invalid,
                None,
                || false,
            ),
            Err(CrystallizationError::Invalid(_))
        ));
    }
    let model = law.equivalent_age_remelting().unwrap();
    assert!(matches!(
        law.apply_specific_enthalpy(initial, f64::NAN, solve, Some(model), || false),
        Err(CrystallizationError::Invalid(_))
    ));
    assert!(matches!(
        law.apply_specific_enthalpy(initial, parent(651.0).unwrap(), solve, Some(model), || {
            false
        }),
        Err(CrystallizationError::EnthalpyOutOfRange { .. })
    ));
    let other = CrystallizationLaw::new(
        &parent,
        &crystal,
        &rate,
        CrystallizationConfig {
            identity: [0x74; 32],
            ..config()
        },
    )
    .unwrap();
    assert_eq!(
        law.apply_specific_enthalpy(
            initial,
            initial.specific_enthalpy_j_kg(),
            solve,
            Some(other.equivalent_age_remelting().unwrap()),
            || false,
        ),
        Err(CrystallizationError::StateMismatch)
    );
    assert_eq!(
        other.apply_specific_enthalpy(
            initial,
            initial.specific_enthalpy_j_kg(),
            solve,
            None,
            || false,
        ),
        Err(CrystallizationError::StateMismatch)
    );
    // Deliberately reusing a declared bundle ID must not permit cached plateau
    // calorics from a different actual source datum.
    let shifted_parent = |t| parent(t).map(|h| h + 100.0);
    let shifted_crystal = |t| crystal(t).map(|h| h + 100.0);
    let shifted =
        CrystallizationLaw::new(&shifted_parent, &shifted_crystal, &rate, config()).unwrap();
    let shifted_state = shifted.state_at_temperature(350.0, 0.4).unwrap();
    assert_eq!(
        shifted.apply_specific_enthalpy(
            shifted_state,
            shifted_state.specific_enthalpy_j_kg(),
            solve,
            Some(model),
            || false,
        ),
        Err(CrystallizationError::StateMismatch)
    );
}

#[test]
fn g4_impulse_cancellation_and_original_source_failures_leave_history_untouched() {
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let initial = law.state_at_temperature(450.0, 0.7).unwrap();
    let target = law
        .specific_enthalpy_j_kg(320.0, initial.crystalline_mass_fraction())
        .unwrap();
    let solve = CrystallizationSolve::default();
    assert_eq!(
        law.apply_specific_enthalpy(initial, target, solve, None, || true),
        Err(CrystallizationError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        law.apply_specific_enthalpy(initial, target, solve, None, || {
            polls += 1;
            polls == 7
        }),
        Err(CrystallizationError::Cancelled)
    );
    assert_eq!(polls, 7);
    assert_eq!(initial, law.state_at_temperature(450.0, 0.7).unwrap());
    let unavailable = |t| {
        if (319.0..321.0).contains(&t) {
            Err(CrystallizationError::Invalid(
                "original caloric source unavailable",
            ))
        } else {
            parent(t)
        }
    };
    let broken = CrystallizationLaw::new(&unavailable, &crystal, &rate, config()).unwrap();
    assert_eq!(
        broken.apply_specific_enthalpy(initial, target, solve, None, || false),
        Err(CrystallizationError::Invalid(
            "original caloric source unavailable"
        ))
    );
}
