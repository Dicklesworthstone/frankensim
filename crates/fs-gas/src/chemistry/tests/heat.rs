use super::*;
use crate::chemistry::isobaric::IsobaricReactor;

#[cfg(feature = "material")]
#[path = "heat/material.rs"]
mod finite_body;

fn thermal_mechanism(binding: f64, rate: f64, reference: f64) -> IdealGasKineticsV1 {
    let model = |name, formation_r| {
        let original = source(name, 0.01, 4.0, 0.0, 0.0);
        let mut card = original.regions()[0].card().clone();
        card.parameters.get_mut("a7").unwrap().value = formation_r;
        Nasa9StandardStateModelV1::new(
            original.species().clone(),
            original.molar_mass(),
            original.convention().clone(),
            vec![Nasa9RegionV1::from_card(card).unwrap()],
        )
        .unwrap()
    };
    mechanism(
        vec![model("A", reference), model("B", reference - binding)],
        vec![1, 1],
        vec![reaction("A-B", &[("A", 1)], &[("B", 1)], rate).irreversible()],
    )
}

fn bath(conductance: f64, temperature: f64) -> HeatBath {
    HeatBath::new(
        ThermalConductance::new(conductance),
        Temperature::new(temperature),
    )
    .unwrap()
}

fn thermal_controls(method: IsochoricKineticsMethod) -> IsochoricKineticsSolve {
    IsochoricKineticsSolve {
        method,
        maximum_step_s: 0.05,
        ..controls(2e-7)
    }
}

struct Endpoint {
    temperature: f64,
    a: f64,
    energy: f64,
    target: f64,
    internal_energy: f64,
    work: f64,
    balance: IsochoricKineticsBalance,
}

fn advance(
    mechanism: IdealGasKineticsV1,
    pressure: bool,
    feed: [Amount; 2],
    temperature: f64,
    duration: f64,
    bath: HeatBath,
    controls: IsochoricKineticsSolve,
) -> Endpoint {
    if pressure {
        let mut reactor = IsobaricReactor::new(mechanism, controls, || false).unwrap();
        let storage = reactor.workspace_bytes();
        let result = reactor
            .advance_with_heat_bath(
                &feed,
                Pressure::new(R * temperature),
                Temperature::new(temperature),
                Time::new(duration),
                bath,
                || false,
            )
            .unwrap()
            .to_owned(|| false)
            .unwrap();
        assert_eq!(storage, reactor.workspace_bytes());
        assert!(result.pressure_residual_pa().abs() < 1e-9);
        Endpoint {
            temperature: result.temperature().value(),
            a: result.amounts()[0].value(),
            energy: result.enthalpy().value(),
            target: result.target_enthalpy().value(),
            internal_energy: result.internal_energy().value(),
            work: result.boundary_work().value(),
            balance: result.balance(),
        }
    } else {
        let mut reactor = IsochoricReactor::new(mechanism, controls, || false).unwrap();
        let storage = reactor.workspace_bytes();
        let result = reactor
            .advance_with_heat_bath(
                &feed,
                Volume::new(1.0),
                Temperature::new(temperature),
                Time::new(duration),
                bath,
                || false,
            )
            .unwrap()
            .to_owned(|| false)
            .unwrap();
        assert_eq!(storage, reactor.workspace_bytes());
        Endpoint {
            temperature: result.temperature().value(),
            a: result.amounts()[0].value(),
            energy: result.internal_energy().value(),
            target: result.target_internal_energy().value(),
            internal_energy: result.internal_energy().value(),
            work: 0.0,
            balance: result.balance(),
        }
    }
}

fn check_balance(result: &Endpoint, initial_u: f64, bath_temperature: f64) {
    let balance = result.balance;
    assert!(balance.energy_residual_j.abs() <= 1e-8);
    assert!((result.energy - result.target - balance.heat_input_j).abs() < 3e-8);
    assert!((result.internal_energy - initial_u + result.work - balance.heat_input_j).abs() < 3e-8);
    assert_eq!(
        balance.entropy_transfer_j_k,
        balance.heat_input_j / bath_temperature
    );
    assert_eq!(
        balance.entropy_production_j_k,
        balance.entropy_change_j_k - balance.entropy_transfer_j_k
    );
    assert!(balance.entropy_production_j_k > 0.0);
    assert!(balance.minimum_step_entropy_production_j_k >= -1e-10);
    assert!(balance.maximum_relative_element_residual < 1e-13);
    assert!(balance.maximum_estimated_local_error <= 1.0);
}

#[test]
fn frozen_populations_heat_and_cool_with_actual_reservoir_entropy() {
    for pressure in [false, true] {
        let capacity = if pressure { 4.0 } else { 3.0 } * R;
        for method in [
            IsochoricKineticsMethod::BackwardEuler,
            IsochoricKineticsMethod::Sdirk2,
        ] {
            for (initial, reservoir) in [(800.0, 300.0), (300.0, 900.0)] {
                // No A is present: every reaction rate is exactly zero, but
                // the original implicit stages must still evolve energy.
                let result = advance(
                    thermal_mechanism(0.0, 2.0, 0.0),
                    pressure,
                    [Amount::new(0.0), Amount::new(1.0)],
                    initial,
                    0.4,
                    bath(40.0, reservoir),
                    thermal_controls(method),
                );
                let expected = reservoir + (initial - reservoir) * (-40.0 / capacity * 0.4).exp();
                let tolerance = if method == IsochoricKineticsMethod::BackwardEuler {
                    0.12
                } else {
                    0.002
                };
                assert!(
                    (result.temperature - expected).abs() < tolerance,
                    "{method:?}, p={pressure}: {} != {expected}",
                    result.temperature
                );
                assert_eq!(result.a, 0.0);
                assert!(
                    (result.balance.heat_input_j - capacity * (result.temperature - initial)).abs()
                        < 2e-8
                );
                assert!(
                    (result.balance.entropy_change_j_k
                        - capacity * (result.temperature / initial).ln())
                    .abs()
                        < 1e-11
                );
                if initial > reservoir {
                    assert!(result.balance.heat_input_j < 0.0);
                    assert!(result.balance.entropy_change_j_k < 0.0);
                    assert!(result.balance.minimum_step_entropy_change_j_k < 0.0);
                } else {
                    assert!(result.balance.heat_input_j > 0.0);
                }
                check_balance(&result, 3.0 * R * initial, reservoir);
            }
        }
    }
}

#[test]
fn simultaneous_reaction_and_heat_match_independent_closed_form() {
    let (rate, binding, duration) = (2.0, 512.0, 0.3);
    for pressure in [false, true] {
        let capacity = if pressure { 4.0 } else { 3.0 } * R;
        for method in [
            IsochoricKineticsMethod::BackwardEuler,
            IsochoricKineticsMethod::Sdirk2,
        ] {
            for (initial, reservoir, relaxation) in [(700.0, 300.0, 3.0), (300.0, 900.0, rate)] {
                let result = advance(
                    thermal_mechanism(binding, rate, 0.0),
                    pressure,
                    [Amount::new(1.0), Amount::new(0.0)],
                    initial,
                    duration,
                    bath(capacity * relaxation, reservoir),
                    thermal_controls(method),
                );
                let remaining = (-rate * duration).exp();
                let thermal = (-relaxation * duration).exp();
                let convolution = if relaxation == rate {
                    duration * thermal
                } else {
                    (remaining - thermal) / (relaxation - rate)
                };
                let expected = reservoir
                    + (initial - reservoir) * thermal
                    + (binding * R * rate / capacity) * convolution;
                let temperature_tolerance = if method == IsochoricKineticsMethod::BackwardEuler {
                    0.12
                } else {
                    0.002
                };
                assert!(
                    (result.temperature - expected).abs() < temperature_tolerance,
                    "{method:?}, p={pressure}: {} != {expected}",
                    result.temperature
                );
                assert!((result.a - remaining).abs() < 3e-4);
                let accounted_heat =
                    capacity * (result.temperature - initial) - binding * R * (1.0 - result.a);
                assert!((result.balance.heat_input_j - accounted_heat).abs() < 3e-8);
                check_balance(&result, 3.0 * R * initial, reservoir);
            }
        }
    }
}

#[test]
fn thermal_reaction_retains_large_elemental_energy_origins() {
    for pressure in [false, true] {
        let mut baseline: Option<(f64, f64, f64)> = None;
        for reference in [0.0, 2.0_f64.powi(40), -2.0_f64.powi(40)] {
            let result = advance(
                thermal_mechanism(512.0, 2.0, reference),
                pressure,
                [Amount::new(1.0), Amount::new(0.0)],
                500.0,
                0.3,
                bath(50.0, 300.0),
                thermal_controls(IsochoricKineticsMethod::Sdirk2),
            );
            assert!(result.balance.energy_residual_j.abs() <= 1e-8);
            if let Some((temperature, amount, heat)) = baseline {
                assert!((result.temperature - temperature).abs() < 1e-6);
                assert!((result.a - amount).abs() < 1e-9);
                assert!((result.balance.heat_input_j - heat).abs() < 1e-5);
            } else {
                baseline = Some((result.temperature, result.a, result.balance.heat_input_j));
            }
        }
    }
}

#[test]
fn thermal_constrained_jacobian_matches_temperature_sensitive_rate_differences() {
    let mechanism = mechanism(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 5.0, -2048.0 * R, 0.0),
        ],
        vec![1, 1],
        (0..5)
            .map(|index| {
                ElementaryReactionV1::new(
                    ReactionId::new(format!("activated-{index}")).unwrap(),
                    vec![(id("A"), 1)],
                    vec![(id("B"), 1)],
                    ArrheniusRateV1::new(
                        Temperature::new(500.0),
                        ReactionProgressRateQuantityV1::new(0.4),
                        1.5,
                        MolarEnergyQuantityV1::new(12000.0),
                    )
                    .unwrap(),
                )
                .unwrap()
                .irreversible()
            })
            .collect(),
    );
    for pressure in [false, true] {
        let mut reactor = solve::Reactor::new(
            mechanism.clone(),
            thermal_controls(IsochoricKineticsMethod::Sdirk2),
        )
        .unwrap();
        if pressure {
            reactor.prepare_pressure().unwrap();
        }
        reactor
            .initialize(
                &[Amount::new(1.0), Amount::new(0.0)],
                Volume::new(1.0),
                pressure.then_some(Pressure::new(R * 500.0)),
                Temperature::new(500.0),
                None,
                &mut || false,
            )
            .unwrap();
        let thermal_bath = bath(1000.0, 700.0);
        let extents = [0.016; 5];
        let (_, analytic) = reactor
            .thermal_rate_linearization(&extents, 500.0, thermal_bath, 0.1)
            .unwrap();
        for column in 0..5 {
            let mut plus = extents;
            let mut minus = extents;
            plus[column] += 1e-6;
            minus[column] -= 1e-6;
            let (plus, _) = reactor
                .thermal_rate_linearization(&plus, 500.0, thermal_bath, 0.1)
                .unwrap();
            let (minus, _) = reactor
                .thermal_rate_linearization(&minus, 500.0, thermal_bath, 0.1)
                .unwrap();
            for row in 0..5 {
                let numeric = (plus[row] - minus[row]) / 2e-6;
                let expected = analytic[row * 5 + column];
                assert!(
                    (numeric - expected).abs() < 2e-7 * expected.abs().max(1.0),
                    "p={pressure}, column {column}, row {row}: {numeric} != {expected}"
                );
            }
        }
    }
}

#[test]
fn insulated_bath_and_cancelled_calls_preserve_prepared_replay() {
    let mechanism = thermal_mechanism(512.0, 2.0, 0.0);
    let controls = thermal_controls(IsochoricKineticsMethod::Sdirk2);
    let feed = [Amount::new(1.0), Amount::new(0.0)];
    let temperature = Temperature::new(500.0);
    let duration = Time::new(0.3);
    let mut volume = IsochoricReactor::new(mechanism.clone(), controls, || false).unwrap();
    let adiabatic = volume
        .advance_adiabatic(&feed, Volume::new(1.0), temperature, duration, || false)
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    let insulated = volume
        .advance_with_heat_bath(
            &feed,
            Volume::new(1.0),
            temperature,
            duration,
            bath(0.0, 900.0),
            || false,
        )
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    assert_eq!(adiabatic, insulated);
    let heated = volume
        .advance_with_heat_bath(
            &feed,
            Volume::new(1.0),
            temperature,
            duration,
            bath(50.0, 300.0),
            || false,
        )
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    let mut polls = 0;
    let refused = volume.advance_with_heat_bath(
        &feed,
        Volume::new(1.0),
        temperature,
        duration,
        bath(50.0, 900.0),
        || {
            polls += 1;
            polls > 300
        },
    );
    assert!(matches!(refused, Err(IsochoricKineticsError::Cancelled)));
    let replay = volume
        .advance_with_heat_bath(
            &feed,
            Volume::new(1.0),
            temperature,
            duration,
            bath(50.0, 300.0),
            || false,
        )
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    assert_eq!(heated, replay);
    assert_eq!(feed, [Amount::new(1.0), Amount::new(0.0)]);

    let mut pressure = IsobaricReactor::new(mechanism, controls, || false).unwrap();
    let p = Pressure::new(R * 500.0);
    let adiabatic = pressure
        .advance_adiabatic(&feed, p, temperature, duration, || false)
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    let insulated = pressure
        .advance_with_heat_bath(&feed, p, temperature, duration, bath(0.0, 900.0), || false)
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    assert_eq!(adiabatic, insulated);
    let heated = pressure
        .advance_with_heat_bath(&feed, p, temperature, duration, bath(50.0, 300.0), || false)
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    let mut polls = 0;
    let refused =
        pressure.advance_with_heat_bath(&feed, p, temperature, duration, bath(50.0, 900.0), || {
            polls += 1;
            polls > 300
        });
    assert!(matches!(refused, Err(IsochoricKineticsError::Cancelled)));
    let replay = pressure
        .advance_with_heat_bath(&feed, p, temperature, duration, bath(50.0, 300.0), || false)
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    assert_eq!(heated, replay);
}

#[test]
fn explicit_initial_energy_precedes_heat_and_bath_inputs_are_admitted() {
    for conductance in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(
            HeatBath::new(
                ThermalConductance::new(conductance),
                Temperature::new(300.0)
            )
            .is_err()
        );
    }
    for temperature in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            HeatBath::new(ThermalConductance::new(1.0), Temperature::new(temperature)).is_err()
        );
    }
    let mechanism = thermal_mechanism(0.0, 2.0, 0.0);
    let controls = thermal_controls(IsochoricKineticsMethod::Sdirk2);
    let feed = [Amount::new(0.0), Amount::new(1.0)];
    let thermal_bath = bath(50.0, 300.0);
    let mut volume = IsochoricReactor::new(mechanism.clone(), controls, || false).unwrap();
    let result = volume
        .advance_energy_with_heat_bath(
            &feed,
            Energy::new(3.0 * R * 600.0),
            Volume::new(1.0),
            Temperature::new(500.0),
            Time::new(0.0),
            thermal_bath,
            || false,
        )
        .unwrap();
    assert!((result.temperature().value() - 600.0).abs() < 1e-9);
    assert_eq!(result.balance().heat_input_j, 0.0);
    assert_eq!(result.balance().entropy_production_j_k, 0.0);
    let mut pressure = IsobaricReactor::new(mechanism, controls, || false).unwrap();
    let result = pressure
        .advance_enthalpy_with_heat_bath(
            &feed,
            Energy::new(4.0 * R * 600.0),
            Pressure::new(R * 500.0),
            Temperature::new(500.0),
            Time::new(0.0),
            thermal_bath,
            || false,
        )
        .unwrap();
    assert!((result.temperature().value() - 600.0).abs() < 1e-9);
    assert!((result.initial_volume().value() - 1.2).abs() < 1e-13);
    assert_eq!(result.boundary_work().value(), 0.0);
    assert_eq!(result.balance().heat_input_j, 0.0);
    assert_eq!(result.balance().entropy_production_j_k, 0.0);
}
