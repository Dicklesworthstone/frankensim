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
fn sensible_bath_refines_to_exponential_and_counts_only_accepted_endpoint_powers() {
    let model = sensible();
    let initial = state(&model, &[(1.0, 0.25)]);
    let capacity: f64 = 50.0 * 50.0 + 0.75 * 30.0 + 0.25 * 90.0;
    let conductance: f64 = 300.0;
    let duration: f64 = 8.0;
    for reservoir in [290.0, 330.0] {
        let exact = reservoir + (300.0 - reservoir) * (-conductance * duration / capacity).exp();
        let mut errors = Vec::new();
        for interval in [2.0, 0.5] {
            let mut ctl = controls(interval);
            ctl.relative_tolerance = 1.0;
            let result = model
                .advance_with_heat_transfer(
                    &initial,
                    Time::new(duration),
                    Temperature::new(reservoir),
                    conductance,
                    ctl,
                    || false,
                )
                .unwrap();
            assert_eq!(result.rejected_steps(), 0);
            assert_eq!(result.accepted_steps(), (duration / interval) as usize);
            // Independent backward-Euler recurrence for a constant sensible
            // capacity. Every accepted comparison supplies TWO physical half
            // powers; its coarse trial is not a third heating contribution.
            let h = 0.5 * interval;
            let mut expected_temperature = 300.0;
            let mut expected_heat = 0.0;
            for _ in 0..2 * result.accepted_steps() {
                expected_temperature = (capacity * expected_temperature
                    + conductance * h * reservoir)
                    / (capacity + conductance * h);
                expected_heat += conductance * h * (reservoir - expected_temperature);
            }
            near(
                result.final_state().temperature().value(),
                expected_temperature,
                2e-9,
            );
            near(result.heat_input_j(), expected_heat, 2e-7);
            near(
                result.enthalpy_change_j(),
                capacity * (expected_temperature - 300.0),
                3e-7,
            );
            assert!(result.energy_residual_j().abs() <= ctl.root.energy_tolerance_j);
            assert_eq!(result.final_state().phases(), initial.phases());
            assert!(
                (300.0_f64.min(reservoir)..=300.0_f64.max(reservoir))
                    .contains(&result.final_state().temperature().value())
            );
            errors.push((result.final_state().temperature().value() - exact).abs());
        }
        assert!(errors[1] < 0.3 * errors[0], "{errors:?}");
    }
}

#[test]
fn reservoir_phase_change_closes_original_energy_and_retains_real_inventories() {
    for enthalpy in [20_000.0, -20_000.0] {
        let model = model(
            vec![component("A", 600.0, enthalpy, 30.0, 30.0)],
            &[1e-4],
            PrecipitationPolicyV1::ExistingSolidOnly,
        );
        let initial = state(&model, &[(2.0, 0.05)]);
        let snapshot = initial.clone();
        let reservoir = 310.0;
        let conductance = 40.0;
        let result = model
            .advance_with_heat_transfer(
                &initial,
                Time::new(5.0),
                Temperature::new(reservoir),
                conductance,
                controls(1.0),
                || false,
            )
            .unwrap();
        let endpoint = result.final_state();
        let amount = endpoint.phases()[0].dissolved().value();
        let temperature = endpoint.temperature().value();
        let independent =
            (2500.0 + 2.0 * 30.0) * (temperature - 300.0) + enthalpy * (amount - 0.05);
        assert!(amount > 0.05);
        assert!(endpoint.phases()[0].solid().value() > 0.0);
        assert_eq!(
            endpoint.phases()[0].total_solute(),
            initial.phases()[0].total_solute()
        );
        assert_eq!(endpoint.solvent_volume(), initial.solvent_volume());
        near(independent, result.heat_input_j(), 1e-7);
        near(
            model.enthalpy_change_j(&initial, endpoint).unwrap(),
            result.heat_input_j(),
            1e-7,
        );
        assert!(result.energy_residual_j().abs() <= 1e-7);
        // Evolving bath temperature matters: replacing its integrated history
        // by the final endpoint power over the whole window gives other heat.
        assert!(
            (result.heat_input_j() - conductance * 5.0 * (reservoir - temperature)).abs() > 0.1
        );
        assert_eq!(initial, snapshot);
    }
}

#[test]
fn external_reservoir_may_exceed_source_support_but_the_bath_cannot() {
    let model = sensible();
    let initial = state(&model, &[(1.0, 0.25)]);
    let result = model
        .advance_with_heat_transfer(
            &initial,
            Time::new(1.0),
            Temperature::new(1000.0),
            1.0,
            controls(0.25),
            || false,
        )
        .unwrap();
    assert!(result.final_state().temperature().value() > 300.0);
    assert!(result.final_state().temperature().value() < 340.0);
    let mut ctl = controls(1.0);
    ctl.minimum_step = Time::new(0.5);
    assert!(matches!(
        model.advance_with_heat_transfer(
            &initial,
            Time::new(1.0),
            Temperature::new(1000.0),
            1e6,
            ctl,
            || false,
        ),
        Err(SolubilityCaloricErrorV1::HeatOutsideBracket { .. })
    ));
}

#[test]
fn zero_exchange_and_zero_duration_preserve_original_limits_and_validation() {
    let model = model(
        vec![component("A", 600.0, 20_000.0, 30.0, 30.0)],
        &[1e-4],
        PrecipitationPolicyV1::ExistingSolidOnly,
    );
    let initial = state(&model, &[(2.0, 0.05)]);
    let original = model
        .advance(&initial, Time::new(1.0), 0.0, controls(0.25), || false)
        .unwrap();
    let insulated = model
        .advance_with_heat_transfer(
            &initial,
            Time::new(1.0),
            Temperature::new(1000.0),
            0.0,
            controls(0.25),
            || false,
        )
        .unwrap();
    assert_eq!(insulated, original);
    let zero = model
        .advance_with_heat_transfer(
            &initial,
            Time::new(0.0),
            Temperature::new(310.0),
            100.0,
            controls(0.25),
            || false,
        )
        .unwrap();
    assert_eq!(zero.final_state(), &initial);
    assert_eq!(zero.heat_input_j(), 0.0);
    assert_eq!(zero.energy_residual_j(), 0.0);
    assert_eq!(zero.accepted_steps(), 0);
    for (temperature, conductance) in [
        (0.0, 1.0),
        (f64::NAN, 1.0),
        (f64::INFINITY, 1.0),
        (310.0, -1.0),
        (310.0, f64::NAN),
        (310.0, f64::INFINITY),
    ] {
        assert!(matches!(
            model.advance_with_heat_transfer(
                &initial,
                Time::new(1.0),
                Temperature::new(temperature),
                conductance,
                controls(0.25),
                || false,
            ),
            Err(SolubilityCaloricErrorV1::Invalid(_))
        ));
    }
}

#[test]
fn cancelled_or_work_limited_exchange_never_publishes_trial_heat_or_state() {
    let model = sensible();
    let initial = state(&model, &[(1.0, 0.25)]);
    let snapshot = initial.clone();
    let mut calls = 0;
    assert_eq!(
        model.advance_with_heat_transfer(
            &initial,
            Time::new(5.0),
            Temperature::new(330.0),
            300.0,
            controls(0.5),
            || {
                calls += 1;
                calls > 40
            },
        ),
        Err(SolubilityCaloricErrorV1::Cancelled)
    );
    assert!(calls > 40);
    let mut ctl = controls(0.5);
    ctl.maximum_attempts = 1;
    assert_eq!(
        model.advance_with_heat_transfer(
            &initial,
            Time::new(5.0),
            Temperature::new(330.0),
            300.0,
            ctl,
            || false,
        ),
        Err(SolubilityCaloricErrorV1::IntegrationLimit)
    );
    let mut ctl = controls(2.0);
    ctl.minimum_step = Time::new(2.0);
    ctl.relative_tolerance = 0.0;
    ctl.temperature_absolute_tolerance_k = 1e-8;
    assert!(matches!(
        model.advance_with_heat_transfer(
            &initial,
            Time::new(2.0),
            Temperature::new(330.0),
            300.0,
            ctl,
            || false,
        ),
        Err(SolubilityCaloricErrorV1::TemporalTolerance { .. })
    ));
    assert_eq!(initial, snapshot);
}

#[test]
fn heat_product_uses_representable_orders_and_refuses_true_overflow_or_underflow() {
    let endpoint = Temperature::new(1.0 - 2.0_f64.powi(-40));
    let source = Reservoir {
        temperature_k: 1.0,
        conductance_w_k: 1e308,
    };
    let heat = source.heat(Time::new(100.0), endpoint).unwrap();
    let expected = 1e308 * (100.0 * (1.0 - endpoint.value()));
    near(heat, expected, 2e-15 * expected);
    let source = Reservoir {
        temperature_k: 1e100,
        conductance_w_k: 1e-300,
    };
    let heat = source
        .heat(Time::new(1e-30), Temperature::new(1.0))
        .unwrap();
    near(heat, 1e-230, 2e-245);
    let source = Reservoir {
        temperature_k: 2.0,
        conductance_w_k: 1e308,
    };
    assert!(
        source
            .heat(Time::new(1e308), Temperature::new(1.0))
            .is_err()
    );
    let source = Reservoir {
        temperature_k: 1.5,
        conductance_w_k: f64::from_bits(1),
    };
    assert!(source.heat(Time::new(0.5), Temperature::new(1.0)).is_err());
    assert_eq!(
        source.heat(Time::new(0.5), Temperature::new(1.5)).unwrap(),
        0.0
    );
}
