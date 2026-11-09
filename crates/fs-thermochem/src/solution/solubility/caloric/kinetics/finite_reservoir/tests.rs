use super::super::tests::{component, controls, model, near, state};
use super::*;
use crate::solution::solubility::kinetics::PrecipitationPolicyV1;

fn sensible() -> SolubilityThermalKineticsV1 {
    model(
        vec![component("A", 600.0, 1000.0, 30.0, 90.0)],
        &[0.0],
        PrecipitationPolicyV1::ExistingSolidOnly,
    )
}

#[test]
fn two_sensible_bodies_refine_to_analytical_exchange_and_count_only_fine_heat() {
    let model = sensible();
    let initial = state(&model, &[(1.0, 0.25)]);
    let bath_capacity: f64 = 2545.0;
    let reservoir_capacity: f64 = 700.0;
    let conductance: f64 = 300.0;
    let duration: f64 = 8.0;
    let reservoir =
        SensibleHeatReservoirV1::new(Temperature::new(330.0), reservoir_capacity).unwrap();
    let equilibrium =
        (bath_capacity * 300.0 + reservoir_capacity * 330.0) / (bath_capacity + reservoir_capacity);
    let decay = (-conductance * (1.0 / bath_capacity + 1.0 / reservoir_capacity) * duration).exp();
    let exact = [
        equilibrium + (300.0 - equilibrium) * decay,
        equilibrium + (330.0 - equilibrium) * decay,
    ];
    let mut errors = Vec::new();
    for interval in [1.0, 0.25] {
        let mut ctl = controls(interval);
        ctl.relative_tolerance = 1.0;
        let result = model
            .advance_with_finite_heat_reservoir(
                &initial,
                Time::new(duration),
                &reservoir,
                conductance,
                ctl,
                || false,
            )
            .unwrap();
        assert_eq!(result.balance().rejected_steps(), 0);
        assert_eq!(
            result.balance().accepted_steps(),
            (duration / interval) as usize
        );
        let h = 0.5 * interval;
        let mut bath = 300.0;
        let mut store = 330.0;
        let mut heat = 0.0;
        // Independent two-capacity backward-Euler algebra. These are exactly
        // the accepted half stages; the coarse comparisons supply no heat.
        for _ in 0..2 * result.balance().accepted_steps() {
            let increment = (store - bath)
                / (1.0 / (conductance * h) + 1.0 / bath_capacity + 1.0 / reservoir_capacity);
            bath += increment / bath_capacity;
            store -= increment / reservoir_capacity;
            heat += increment;
        }
        let actual_bath = result.balance().final_state().temperature().value();
        let actual_store = result.reservoir().temperature().value();
        near(actual_bath, bath, 3e-9);
        near(actual_store, store, 3e-9);
        near(result.balance().heat_input_j(), heat, 3e-7);
        near(result.reservoir_enthalpy_change_j(), -heat, 3e-7);
        near(
            reservoir.enthalpy_change_j(result.reservoir()).unwrap(),
            -heat,
            3e-7,
        );
        assert!(result.closed_system_energy_residual_j().abs() <= 1e-7);
        assert_eq!(result.balance().final_state().phases(), initial.phases());
        assert!(actual_bath <= actual_store);
        errors.push((actual_bath - exact[0]).abs() + (actual_store - exact[1]).abs());
    }
    assert!(errors[1] < 0.3 * errors[0], "{errors:?}");
}

#[test]
fn stiff_exchange_reaches_common_temperature_without_exposure_overflow() {
    let model = sensible();
    let initial = state(&model, &[(1.0, 0.25)]);
    let reservoir = SensibleHeatReservoirV1::new(Temperature::new(330.0), 500.0).unwrap();
    let result = model
        .advance_with_finite_heat_reservoir(
            &initial,
            Time::new(100.0),
            &reservoir,
            1e308,
            controls(100.0),
            || false,
        )
        .unwrap();
    let equilibrium = (2545.0 * 300.0 + 500.0 * 330.0) / 3045.0;
    near(
        result.balance().final_state().temperature().value(),
        equilibrium,
        2e-10,
    );
    near(result.reservoir().temperature().value(), equilibrium, 2e-10);
    near(
        result.balance().heat_input_j(),
        2545.0 * (equilibrium - 300.0),
        2e-7,
    );
    assert!(result.closed_system_energy_residual_j().abs() <= 1e-7);
    assert_eq!(result.balance().accepted_steps(), 1);
}

#[test]
fn small_reservoir_error_cannot_hide_behind_a_large_bath_capacity() {
    let model = sensible();
    let initial = state(&model, &[(1.0, 0.25)]);
    let reservoir_capacity: f64 = 0.01;
    let reservoir =
        SensibleHeatReservoirV1::new(Temperature::new(330.0), reservoir_capacity).unwrap();
    let conductance: f64 = 0.01;
    let mut ctl = controls(1.0);
    ctl.relative_tolerance = 0.0;
    ctl.temperature_absolute_tolerance_k = 1e-3;
    let result = model
        .advance_with_finite_heat_reservoir(
            &initial,
            Time::new(1.0),
            &reservoir,
            conductance,
            ctl,
            || false,
        )
        .unwrap();
    let equilibrium = (2545.0 * 300.0 + reservoir_capacity * 330.0) / (2545.0 + reservoir_capacity);
    let exact_store = equilibrium
        + (330.0 - equilibrium) * (-conductance * (1.0 / 2545.0 + 1.0 / reservoir_capacity)).exp();
    assert!(result.balance().rejected_steps() > 0);
    assert!(result.balance().accepted_steps() > 8);
    // Bath-only control would accept the initial full interval: its change is
    // below0.001K, while two coarse half steps miss reservoir T by over2K.
    near(result.reservoir().temperature().value(), exact_store, 0.15);
    assert!((result.balance().final_state().temperature().value() - 300.0).abs() < 1e-3);
    assert!(result.closed_system_energy_residual_j().abs() <= 1e-7);
}

#[test]
fn repeated_windows_keep_energy_when_published_reservoir_temperature_is_unchanged() {
    let model = sensible();
    let mut bath = state(&model, &[(1.0, 0.25)]);
    let original = SensibleHeatReservoirV1::new(Temperature::new(310.0), 1e20).unwrap();
    let mut reservoir = original;
    let mut heat = 0.0;
    for _ in 0..8 {
        let result = model
            .advance_with_finite_heat_reservoir(
                &bath,
                Time::new(0.25),
                &reservoir,
                1.0,
                controls(0.25),
                || false,
            )
            .unwrap();
        heat += result.balance().heat_input_j();
        let (balance, next) = result.into_parts();
        bath = balance.final_state().clone();
        reservoir = next;
        assert_eq!(reservoir.temperature(), original.temperature());
    }
    assert!(heat > 0.0);
    near(reservoir.energy_change_j(), -heat, 3e-14 * heat);
    near(
        original.enthalpy_change_j(&reservoir).unwrap(),
        -heat,
        3e-14 * heat,
    );
    assert_eq!(
        reservoir.reference_temperature(),
        original.reference_temperature()
    );

    // Even a subsequent1J exchange beside a retained1e16J change survives,
    // though BOTH the rounded total-energy and temperature getters stay fixed.
    let large = original.after_heat(HeatReceipt::from_value(1e16)).unwrap();
    let trace = large.after_heat(HeatReceipt::from_value(1.0)).unwrap();
    assert_eq!(large.temperature(), trace.temperature());
    assert_eq!(large.energy_change_j(), trace.energy_change_j());
    assert_eq!(large.enthalpy_change_j(&trace).unwrap(), -1.0);
    let second = trace.after_heat(HeatReceipt::from_value(1.0)).unwrap();
    assert_eq!(large.enthalpy_change_j(&second).unwrap(), -2.0);
}

#[test]
fn finite_heat_drives_phase_exchange_and_closes_both_original_energy_accounts() {
    let model = model(
        vec![component("A", 600.0, 20_000.0, 30.0, 30.0)],
        &[1e-4],
        PrecipitationPolicyV1::ExistingSolidOnly,
    );
    let initial = state(&model, &[(2.0, 0.05)]);
    let reservoir = SensibleHeatReservoirV1::new(Temperature::new(330.0), 500.0).unwrap();
    let result = model
        .advance_with_finite_heat_reservoir(
            &initial,
            Time::new(5.0),
            &reservoir,
            100.0,
            controls(0.5),
            || false,
        )
        .unwrap();
    let bath = result.balance().final_state();
    let amount = bath.phases()[0].dissolved().value();
    let bath_energy = 2560.0 * (bath.temperature().value() - 300.0) + 20_000.0 * (amount - 0.05);
    let reservoir_energy = 500.0 * (result.reservoir().temperature().value() - 330.0);
    assert!(amount > 0.05);
    assert_eq!(
        bath.phases()[0].total_solute(),
        initial.phases()[0].total_solute()
    );
    assert_eq!(bath.solvent_volume(), initial.solvent_volume());
    assert!(result.balance().heat_input_j() > 0.0);
    assert!(result.reservoir().temperature().value() < 330.0);
    near(bath_energy, result.balance().heat_input_j(), 1e-7);
    near(reservoir_energy, result.reservoir_enthalpy_change_j(), 1e-7);
    near(bath_energy + reservoir_energy, 0.0, 1e-7);
    assert!(result.closed_system_energy_residual_j().abs() <= 1e-7);
}

#[test]
fn finite_store_validation_zero_limits_and_cancellation_preserve_both_inputs() {
    let model = sensible();
    let initial = state(&model, &[(1.0, 0.25)]);
    let snapshot = initial.clone();
    let reservoir = SensibleHeatReservoirV1::new(Temperature::new(330.0), 500.0).unwrap();
    for (temperature, capacity) in [
        (0.0, 1.0),
        (f64::NAN, 1.0),
        (300.0, 0.0),
        (300.0, -1.0),
        (300.0, f64::INFINITY),
    ] {
        assert!(SensibleHeatReservoirV1::new(Temperature::new(temperature), capacity).is_err());
    }
    let other = SensibleHeatReservoirV1::new(Temperature::new(330.0), 501.0).unwrap();
    assert!(reservoir.enthalpy_change_j(&other).is_err());
    let adiabatic = model
        .advance(&initial, Time::new(1.0), 0.0, controls(0.25), || false)
        .unwrap();
    let insulated = model
        .advance_with_finite_heat_reservoir(
            &initial,
            Time::new(1.0),
            &reservoir,
            0.0,
            controls(0.25),
            || false,
        )
        .unwrap();
    assert_eq!(insulated.balance(), &adiabatic);
    assert_eq!(insulated.reservoir(), &reservoir);
    let zero = model
        .advance_with_finite_heat_reservoir(
            &initial,
            Time::new(0.0),
            &reservoir,
            100.0,
            controls(0.25),
            || false,
        )
        .unwrap();
    assert_eq!(zero.balance().final_state(), &initial);
    assert_eq!(zero.reservoir(), &reservoir);
    assert_eq!(zero.closed_system_energy_residual_j(), 0.0);
    let mut polls = 0;
    assert_eq!(
        model.advance_with_finite_heat_reservoir(
            &initial,
            Time::new(5.0),
            &reservoir,
            100.0,
            controls(0.5),
            || {
                polls += 1;
                polls > 40
            },
        ),
        Err(SolubilityCaloricErrorV1::Cancelled)
    );
    let mut ctl = controls(0.5);
    ctl.maximum_attempts = 1;
    assert_eq!(
        model.advance_with_finite_heat_reservoir(
            &initial,
            Time::new(5.0),
            &reservoir,
            100.0,
            ctl,
            || false,
        ),
        Err(SolubilityCaloricErrorV1::IntegrationLimit)
    );
    assert_eq!(initial, snapshot);
    assert_eq!(reservoir.energy_change_j(), 0.0);
    let external = SensibleHeatReservoirV1::new(Temperature::new(1000.0), 1.0).unwrap();
    let supported = model
        .advance_with_finite_heat_reservoir(
            &initial,
            Time::new(1.0),
            &external,
            1000.0,
            controls(0.25),
            || false,
        )
        .unwrap();
    assert!(supported.balance().final_state().temperature().value() < 340.0);
    // An arbitrary absolute C*T may overflow while every actual exchange fits.
    let enormous = SensibleHeatReservoirV1::new(Temperature::new(310.0), 1e308).unwrap();
    let trace = model
        .advance_with_finite_heat_reservoir(
            &initial,
            Time::new(0.1),
            &enormous,
            1.0,
            controls(0.1),
            || false,
        )
        .unwrap();
    assert!(trace.reservoir_enthalpy_change_j() < 0.0);
    assert_eq!(trace.reservoir().temperature(), enormous.temperature());
}

#[test]
fn harmonic_product_preserves_heat_with_subnormal_coefficients() {
    let minimum = f64::from_bits(1);
    let exposure = minimum * 1e308;
    let source = FiniteReservoir {
        reservoir: SensibleHeatReservoirV1::new(Temperature::new(3.0), exposure).unwrap(),
        conductance_w_k: minimum,
    };
    let heat = source
        .heat(
            HeatReceipt::default(),
            Time::new(1e308),
            Temperature::new(1.0),
        )
        .unwrap();
    near(heat, exposure, 3e-15 * exposure);
    let source = FiniteReservoir {
        reservoir: SensibleHeatReservoirV1::new(Temperature::new(1e308), minimum).unwrap(),
        conductance_w_k: 1.0,
    };
    let heat = source
        .heat(
            HeatReceipt::default(),
            Time::new(2.0 * minimum),
            Temperature::new(1.0),
        )
        .unwrap();
    let expected = (minimum * 1e308) / 1.5;
    near(heat, expected, 3e-15 * expected);
}
