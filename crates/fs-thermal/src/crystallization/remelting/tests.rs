use super::*;

fn config() -> CrystallizationConfig {
    CrystallizationConfig {
        identity: [0x52; 32],
        avrami_exponent: 2.0,
        maximum_crystalline_mass_fraction: 0.6,
        temperature_range_k: [250.0, 650.0],
        melting_temperature_k: 500.0,
    }
}

fn amorphous(t: f64) -> Result<f64, CrystallizationError> {
    Ok(2000.0 * (t - 300.0) + 100_000.0)
}

fn crystalline(t: f64) -> Result<f64, CrystallizationError> {
    Ok(1000.0 * (t - 300.0))
}

fn rate(t: f64) -> Result<f64, CrystallizationError> {
    Ok(0.001 * (500.0 - t).max(0.0))
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} within {tolerance}"
    );
}

#[test]
fn g1_partial_melting_uses_original_latent_gap_and_equivalent_surviving_age() {
    let law = CrystallizationLaw::new(&amorphous, &crystalline, &rate, config()).unwrap();
    let remelting = law.equivalent_age_remelting().unwrap();
    let solve = CrystallizationSolve::default();
    let initial = law.state_at_temperature(500.0, 1.0).unwrap();
    let expected_fraction = initial.crystalline_mass_fraction() * 0.5;
    let target_h = 500_000.0 - 300_000.0 * expected_fraction;
    let step = law
        .advance_specific_enthalpy_with_remelting(initial, target_h, 3.0, solve, remelting, || {
            false
        })
        .unwrap();

    assert_eq!(remelting.identity(), config().identity);
    assert_eq!(step.initial, initial);
    assert_eq!(step.state.temperature_k(), 500.0);
    assert_eq!(step.state.specific_enthalpy_j_kg(), target_h);
    near(
        step.state.crystalline_mass_fraction(),
        expected_fraction,
        5e-14,
    );
    let expected_age = (-(1.0 - expected_fraction / 0.6).ln()).sqrt();
    near(step.state.exposure(), expected_age, 5e-13);
    near(step.specific_heat_j_kg, 300_000.0 * expected_fraction, 1e-9);
    assert!(step.caloric_residual_j_kg.abs() <= solve.specific_enthalpy_tolerance_j_kg);
    assert_eq!(
        step.state.crystalline_mass_fraction().to_bits(),
        law.crystalline_mass_fraction(step.state.exposure())
            .unwrap()
            .to_bits()
    );
    law.validate_state(step.state, solve).unwrap();

    // The original irreversible API keeps its previous refusal contract.
    assert!(matches!(
        law.advance_specific_enthalpy(initial, target_h, 3.0, solve, || false),
        Err(CrystallizationError::EnthalpyOutOfRange { .. })
    ));
    let dwell = law
        .advance_specific_enthalpy_with_remelting(
            step.state,
            target_h,
            10.0,
            solve,
            remelting,
            || false,
        )
        .unwrap();
    assert_eq!(dwell.state, step.state);
    assert_eq!(dwell.specific_heat_j_kg, 0.0);

    // Plateau projection admits an absolute residual. Exercise a validated
    // accepted state one enthalpy ULP below its caloric reconstruction: a fresh
    // inverse would spuriously move below Tm during this zero-heat hold.
    let rounded = CrystallizationState {
        specific_enthalpy_j_kg: f64::from_bits(step.state.specific_enthalpy_j_kg().to_bits() - 1),
        ..step.state
    };
    law.validate_state(rounded, solve).unwrap();
    let dwell = law
        .advance_specific_enthalpy_with_remelting(
            rounded,
            rounded.specific_enthalpy_j_kg(),
            10.0,
            solve,
            remelting,
            || false,
        )
        .unwrap();
    assert_eq!(dwell.state, rounded);
    assert_eq!(dwell.specific_heat_j_kg, 0.0);
}

#[test]
fn g1_complete_melting_superheat_and_recrystallization_share_one_enthalpy_account() {
    let law = CrystallizationLaw::new(&amorphous, &crystalline, &rate, config()).unwrap();
    let remelting = law.equivalent_age_remelting().unwrap();
    let solve = CrystallizationSolve::default();
    let initial = law.state_at_temperature(500.0, 1.0).unwrap();
    let mut current = initial;
    let mut heat = 0.0;
    for temperature in [500.0, 550.0, 500.0] {
        let target_h = amorphous(temperature).unwrap();
        let step = law
            .advance_specific_enthalpy_with_remelting(
                current,
                target_h,
                1.0,
                solve,
                remelting,
                || false,
            )
            .unwrap();
        near(step.state.temperature_k(), temperature, 1e-9);
        assert_eq!(step.state.exposure(), 0.0);
        assert_eq!(step.state.crystalline_mass_fraction(), 0.0);
        assert_eq!(step.state.specific_enthalpy_j_kg(), target_h);
        assert!(step.caloric_residual_j_kg.abs() <= solve.specific_enthalpy_tolerance_j_kg);
        heat += step.specific_heat_j_kg;
        current = step.state;
    }

    // After a complete melt, cooling starts at zero effective exposure.
    // The amorphous 400 K midpoint gives K=0.1/s, hence z=0.1 at one second.
    let cooled = law
        .advance_specific_enthalpy_with_remelting(current, 100_000.0, 1.0, solve, remelting, || {
            false
        })
        .unwrap();
    let expected_fraction = 0.6 * (1.0 - (-0.01_f64).exp());
    let expected_t = 300.0 + 100_000.0 * expected_fraction / (2000.0 - 1000.0 * expected_fraction);
    near(cooled.state.exposure(), 0.1, 1e-11);
    near(
        cooled.state.crystalline_mass_fraction(),
        expected_fraction,
        1e-12,
    );
    near(cooled.state.temperature_k(), expected_t, 1e-9);
    assert!(cooled.state.temperature_k() > 300.0);
    law.validate_state(cooled.state, solve).unwrap();
    heat += cooled.specific_heat_j_kg;
    near(initial.specific_enthalpy_j_kg() + heat, 100_000.0, 1e-9);
}

#[test]
fn g0_adiabatic_recalescence_cannot_be_erased_by_a_projected_midpoint() {
    let fast_rate = |t: f64| Ok(0.1 * (500.0 - t).max(0.0));
    let law = CrystallizationLaw::new(&amorphous, &crystalline, &fast_rate, config()).unwrap();
    let remelting = law.equivalent_age_remelting().unwrap();
    let solve = CrystallizationSolve::default();
    let initial = law.state_at_temperature(480.0, 0.0).unwrap();
    let h = initial.specific_enthalpy_j_kg();

    // The oversized midpoint grows more crystals than this fixed energy can
    // support below Tm. Projecting to Tm sets Kmid=0; silently accepting old z
    // here would keep this undercooled state frozen on every repeated step.
    assert_eq!(
        law.advance_specific_enthalpy_with_remelting(initial, h, 2.0, solve, remelting, || false),
        Err(CrystallizationError::NonConverged)
    );
    assert_eq!(initial, law.state_at_temperature(480.0, 0.0).unwrap());

    let resolved = law
        .advance_specific_enthalpy_with_remelting(initial, h, 0.01, solve, remelting, || false)
        .unwrap();
    assert!(resolved.state.crystalline_mass_fraction() > 0.0);
    assert!(resolved.state.temperature_k() > 480.0);
    assert!(resolved.state.temperature_k() < 500.0);
    assert_eq!(resolved.state.specific_enthalpy_j_kg(), h);
    assert_eq!(resolved.specific_heat_j_kg, 0.0);
    law.validate_state(resolved.state, solve).unwrap();
}

#[test]
fn g0_partial_remelting_respects_large_datums_and_saturated_finite_ages() {
    let solve = CrystallizationSolve::default();
    for datum in [0.0, 1e12] {
        let a = |t: f64| amorphous(t).map(|h| datum + h);
        let c = |t: f64| crystalline(t).map(|h| datum + h);
        let law = CrystallizationLaw::new(&a, &c, &rate, config()).unwrap();
        let remelting = law.equivalent_age_remelting().unwrap();
        let initial = law.state_at_temperature(500.0, f64::MAX).unwrap();
        assert_eq!(initial.crystalline_mass_fraction(), 0.6);
        let target_h = datum + 425_000.0;
        let step = law
            .advance_specific_enthalpy_with_remelting(
                initial,
                target_h,
                1.0,
                solve,
                remelting,
                || false,
            )
            .unwrap();
        assert_eq!(step.state.temperature_k(), 500.0);
        assert_eq!(step.state.specific_enthalpy_j_kg(), target_h);
        near(step.state.crystalline_mass_fraction(), 0.25, 1e-12);
        near(
            step.state.exposure(),
            (-(1.0_f64 - 0.25 / 0.6).ln()).sqrt(),
            1e-11,
        );
        assert_eq!(step.specific_heat_j_kg, 105_000.0);
        // The common datum cannot enlarge the original absolute energy budget.
        assert!(step.caloric_residual_j_kg.abs() <= solve.specific_enthalpy_tolerance_j_kg);
        law.validate_state(step.state, solve).unwrap();
    }
}

#[test]
fn g4_remelting_cancellation_and_source_policy_refusals_publish_no_history() {
    let law = CrystallizationLaw::new(&amorphous, &crystalline, &rate, config()).unwrap();
    let remelting = law.equivalent_age_remelting().unwrap();
    let solve = CrystallizationSolve::default();
    let initial = law.state_at_temperature(500.0, 1.0).unwrap();
    let mut calls = 0;
    assert_eq!(
        law.advance_specific_enthalpy_with_remelting(
            initial,
            475_000.0,
            1.0,
            solve,
            remelting,
            || {
                calls += 1;
                calls == 4
            }
        ),
        Err(CrystallizationError::Cancelled)
    );
    assert_eq!(calls, 4);
    assert_eq!(initial, law.state_at_temperature(500.0, 1.0).unwrap());

    let other = CrystallizationLaw::new(
        &amorphous,
        &crystalline,
        &rate,
        CrystallizationConfig {
            identity: [0x53; 32],
            ..config()
        },
    )
    .unwrap();
    assert_eq!(
        law.advance_specific_enthalpy_with_remelting(
            initial,
            475_000.0,
            1.0,
            solve,
            other.equivalent_age_remelting().unwrap(),
            || false
        ),
        Err(CrystallizationError::StateMismatch)
    );
    let subcooled = CrystallizationLaw::new(
        &amorphous,
        &crystalline,
        &rate,
        CrystallizationConfig {
            temperature_range_k: [250.0, 499.0],
            ..config()
        },
    )
    .unwrap();
    assert!(matches!(
        subcooled.equivalent_age_remelting(),
        Err(CrystallizationError::TemperatureOutOfRange { .. })
    ));
    assert!(matches!(
        law.advance_specific_enthalpy_with_remelting(
            initial,
            amorphous(651.0).unwrap(),
            1.0,
            solve,
            remelting,
            || false
        ),
        Err(CrystallizationError::EnthalpyOutOfRange { .. })
    ));
}
