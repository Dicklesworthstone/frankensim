use super::*;

#[test]
fn prepared_third_body_isomerization_matches_density_dependent_relaxation() {
    use fs_thermochem::kinetics::{PressureDependenceV1, ThirdBodyEfficienciesV1};

    let reference_rate = 0.2;
    let reaction = reaction("A-B-M", &[("A", 1)], &[("B", 1)], reference_rate)
        .with_pressure_dependence(PressureDependenceV1::ThirdBody(
            ThirdBodyEfficienciesV1::default(),
        ));
    let mechanism = mechanism(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 4.0, 0.0, 0.0),
        ],
        vec![1, 1],
        vec![reaction],
    );
    let mut prepared = IsochoricReactor::new(mechanism, controls(1e-5), || false).unwrap();
    let duration = 0.05;
    let mut previous_product = 0.0;
    for volume in [0.02, 0.005] {
        // A <-> B conserves molecular count and has identical source calorics.
        // Unit efficiencies therefore give a constant collider concentration
        // of 1 mol / V, making the full pressure-dependent ODE analytic.
        let effective_rate = reference_rate / volume;
        let expected_product = 0.5 * (1.0 - (-2.0 * effective_rate * duration).exp());
        let endpoint = prepared
            .advance_adiabatic(
                &[Amount::new(1.0), Amount::new(0.0)],
                Volume::new(volume),
                Temperature::new(500.0),
                Time::new(duration),
                || false,
            )
            .unwrap();
        let actual_product = endpoint.amounts()[1].value();
        assert!((actual_product - expected_product).abs() < 4e-4);
        assert!(actual_product > previous_product);
        assert!((endpoint.temperature().value() - 500.0).abs() < 1e-8);
        assert!(endpoint.balance().energy_residual_j.abs() <= controls(1e-5).energy_tolerance_j);
        assert!(endpoint.balance().entropy_change_j_k > 0.0);
        previous_product = actual_product;
    }
}

#[test]
fn prepared_reactor_reuses_storage_for_independent_dense_cells() {
    let mechanism = isomers(0.0, 3.0);
    let solve = controls(1e-5);
    let mut prepared = IsochoricReactor::new(mechanism.clone(), solve, || false).unwrap();
    assert_eq!(prepared.mechanism(), &mechanism);
    assert_eq!(prepared.controls(), solve);
    let output_storage = prepared.amounts.as_ptr();
    let stage_storage = prepared.current.as_ptr();
    let concentration_storage = prepared.current_concentrations.as_ptr();
    for (a, b, volume, temperature, duration) in [
        (1.0, 0.0, 0.01, 500.0, 0.1),
        (0.0, 3.0, 0.027, 900.0, 0.13),
        (0.3, 0.7, 0.004, 600.0, 0.0),
        (2.0, 0.1, 0.015, 300.0, 0.2),
    ] {
        let input = [Amount::new(a), Amount::new(b)];
        let endpoint = prepared
            .advance_adiabatic(
                &input,
                Volume::new(volume),
                Temperature::new(temperature),
                Time::new(duration),
                || false,
            )
            .unwrap();
        let total = a + b;
        let expected_b = total / 2.0 + (b - total / 2.0) * (-6.0 * duration).exp();
        assert!((endpoint.amounts()[1].value() - expected_b).abs() < 4e-4 * total);
        assert!((endpoint.temperature().value() - temperature).abs() < 1e-8);
        assert_eq!(endpoint.amounts().as_ptr(), output_storage);
        assert_eq!(endpoint.elapsed_time(), Time::new(duration));
        assert_eq!(endpoint.volume(), Volume::new(volume));
        assert!(endpoint.balance().energy_residual_j.abs() <= solve.energy_tolerance_j);
        assert_eq!(input, [Amount::new(a), Amount::new(b)]);
        assert_eq!(prepared.current.as_ptr(), stage_storage);
        assert_eq!(
            prepared.current_concentrations.as_ptr(),
            concentration_storage
        );
    }
}

#[test]
fn prepared_explicit_energy_matches_owned_endpoint_and_keeps_source_identity() {
    let mechanism = isomers(7000.0, 8.0);
    let solve = controls(1e-5);
    let feed = [(id("A"), Amount::new(0.4)), (id("B"), Amount::new(0.6))];
    let source_energy = advance_adiabatic(
        &mechanism,
        &feed,
        Volume::new(0.03),
        Temperature::new(550.0),
        Time::new(0.0),
        solve,
        || false,
    )
    .unwrap()
    .internal_energy();
    let expected = advance_adiabatic_energy(
        &mechanism,
        &feed,
        source_energy,
        Volume::new(0.03),
        Temperature::new(800.0),
        Time::new(0.05),
        solve,
        || false,
    )
    .unwrap();
    let mut prepared = IsochoricReactor::new(mechanism.clone(), solve, || false).unwrap();
    let endpoint = prepared
        .advance_adiabatic_energy(
            &[Amount::new(0.4), Amount::new(0.6)],
            source_energy,
            Volume::new(0.03),
            Temperature::new(800.0),
            Time::new(0.05),
            || false,
        )
        .unwrap();
    assert_eq!(endpoint.mechanism(), &mechanism);
    assert_eq!(endpoint.target_internal_energy(), source_energy);
    assert_eq!(endpoint.to_owned(|| false).unwrap(), expected);
    assert_eq!(expected.sources().len(), 2);
    assert!(expected.temperature().value() > 550.0);
}

#[test]
fn prepared_workspace_recovers_after_late_cancellation_and_bad_energy() {
    let mechanism = isomers(7000.0, 8.0);
    let solve = controls(1e-5);
    let mut prepared = IsochoricReactor::new(mechanism.clone(), solve, || false).unwrap();
    let input = [Amount::new(1.0), Amount::new(0.0)];
    let mut polls = 0;
    let result = prepared.advance_adiabatic(
        &input,
        Volume::new(0.01),
        Temperature::new(500.0),
        Time::new(0.1),
        || {
            polls += 1;
            polls > 80
        },
    );
    assert!(matches!(result, Err(IsochoricKineticsError::Cancelled)));
    assert!(polls > 80);
    let result = prepared.advance_adiabatic_energy(
        &input,
        Energy::new(1e8),
        Volume::new(0.01),
        Temperature::new(500.0),
        Time::new(0.1),
        || false,
    );
    assert!(matches!(
        result,
        Err(IsochoricKineticsError::EnergyNotBracketed { .. })
    ));
    let expected = run(&mechanism, 500.0, 0.1, solve);
    let recovered = prepared
        .advance_adiabatic(
            &input,
            Volume::new(0.01),
            Temperature::new(500.0),
            Time::new(0.1),
            || false,
        )
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    assert_eq!(recovered, expected);
    assert_eq!(input, [Amount::new(1.0), Amount::new(0.0)]);
}

#[test]
fn prepared_attempt_limit_does_not_contaminate_the_next_cell() {
    let mechanism = isomers(0.0, 3.0);
    let solve = IsochoricKineticsSolve {
        maximum_steps: 1,
        ..controls(1e-5)
    };
    let mut prepared = IsochoricReactor::new(mechanism, solve, || false).unwrap();
    let result = prepared.advance_adiabatic(
        &[Amount::new(1.0), Amount::new(0.0)],
        Volume::new(0.01),
        Temperature::new(500.0),
        Time::new(0.1),
        || false,
    );
    assert!(matches!(
        result,
        Err(IsochoricKineticsError::WorkLimit { .. })
    ));
    let endpoint = prepared
        .advance_adiabatic(
            &[Amount::new(0.2), Amount::new(0.8)],
            Volume::new(0.025),
            Temperature::new(750.0),
            Time::new(0.0),
            || false,
        )
        .unwrap();
    assert_eq!(endpoint.amounts(), &[Amount::new(0.2), Amount::new(0.8)]);
    assert_eq!(endpoint.balance().accepted_steps, 0);
    assert_eq!(endpoint.balance().rejected_steps, 0);
    assert!((endpoint.temperature().value() - 750.0).abs() < 1e-10);
}
