use super::*;
use crate::chemistry::heat::gas::{FiniteGasHeatBody, FiniteGasHeatReport, GasHeatConstraint};
use crate::chemistry::heat::{HeatExchange, HeatStage};
use fs_qty::Area;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis};

const SIGMA: f64 = 5.670_374_419e-8;

fn coolant(variable: bool, datum: f64) -> FrozenIdealGasMixtureModelV1 {
    let model = |name, mass, capacity, slope| {
        let original = source(name, mass, capacity, datum * R, datum);
        let mut card = original.regions()[0].card().clone();
        card.parameters.get_mut("a3").unwrap().value = slope;
        Nasa9StandardStateModelV1::new(
            original.species().clone(),
            original.molar_mass(),
            original.convention().clone(),
            vec![Nasa9RegionV1::from_card(card).unwrap()],
        )
        .unwrap()
    };
    FrozenIdealGasMixtureModelV1::new(
        vec![
            model("coolant-A", 0.02, 3.0, if variable { 0.002 } else { 0.0 }),
            model("coolant-B", 0.04, 5.0, if variable { 0.004 } else { 0.0 }),
        ],
        Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.75]).unwrap(),
    )
    .unwrap()
}

fn body(
    model: &FrozenIdealGasMixtureModelV1,
    pressure: bool,
    initial: f64,
    amount: f64,
    g: f64,
    bounds: [f64; 2],
) -> FiniteGasHeatBody<'_> {
    let constraint = if pressure {
        GasHeatConstraint::ConstantPressure(Pressure::new(2e5))
    } else {
        GasHeatConstraint::ConstantVolume(Volume::new(amount * R * initial / 2e5))
    };
    FiniteGasHeatBody::new(
        model,
        Amount::new(amount),
        constraint,
        Temperature::new(initial),
        bounds.map(Temperature::new),
        ThermalConductance::new(g),
        || false,
    )
    .unwrap()
}

struct Pair<'m> {
    gas_t: f64,
    a: f64,
    gas_u: f64,
    gas_work: f64,
    balance: IsochoricKineticsBalance,
    body: FiniteGasHeatBody<'m>,
    report: FiniteGasHeatReport,
}

fn run<'m>(
    pressure: bool,
    method: IsochoricKineticsMethod,
    initial: f64,
    duration: f64,
    original: FiniteGasHeatBody<'m>,
    binding: f64,
    reacting: bool,
) -> Pair<'m> {
    let mechanism = thermal_mechanism(binding, 2.0, 0.0);
    let feed = if reacting {
        [Amount::new(1.0), Amount::new(0.0)]
    } else {
        [Amount::new(0.0), Amount::new(1.0)]
    };
    let controls = thermal_controls(method);
    if pressure {
        let mut reactor = IsobaricReactor::new(mechanism, controls, || false).unwrap();
        let bytes = reactor.workspace_bytes();
        let result = reactor
            .advance_with_finite_gas_body(
                &feed,
                Pressure::new(R * initial),
                Temperature::new(initial),
                Time::new(duration),
                original,
                || false,
            )
            .unwrap();
        let gas = result.gas();
        let pair = Pair {
            gas_t: gas.temperature().value(),
            a: gas.amounts()[0].value(),
            gas_u: gas.internal_energy().value(),
            gas_work: gas.boundary_work().value(),
            balance: gas.balance(),
            body: result.body(),
            report: *result.body_report(),
        };
        assert_eq!(bytes, reactor.workspace_bytes());
        pair
    } else {
        let mut reactor = IsochoricReactor::new(mechanism, controls, || false).unwrap();
        let bytes = reactor.workspace_bytes();
        let result = reactor
            .advance_with_finite_gas_body(
                &feed,
                Volume::new(1.0),
                Temperature::new(initial),
                Time::new(duration),
                original,
                || false,
            )
            .unwrap();
        let gas = result.gas();
        let pair = Pair {
            gas_t: gas.temperature().value(),
            a: gas.amounts()[0].value(),
            gas_u: gas.internal_energy().value(),
            gas_work: 0.0,
            balance: gas.balance(),
            body: result.body(),
            report: *result.body_report(),
        };
        assert_eq!(bytes, reactor.workspace_bytes());
        pair
    }
}

fn balances(result: &Pair<'_>, original: FiniteGasHeatBody<'_>, initial_gas: f64) {
    assert!(result.balance.energy_residual_j.abs() <= 1e-8);
    assert!(result.report.heat_residual_j.abs() <= 1e-8);
    assert_eq!(result.balance.heat_input_j, -result.report.heat_input_j);
    let first_law = result.gas_u - 3.0 * R * initial_gas
        + result.gas_work
        + result.report.internal_energy_change_j
        + result.report.boundary_work_j;
    assert!(first_law.abs() < 4e-8, "combined first law: {first_law}");
    assert_eq!(
        result.balance.entropy_transfer_j_k,
        -result.report.entropy_change_j_k
    );
    assert!(
        (result.balance.entropy_production_j_k
            - result.balance.entropy_change_j_k
            - result.report.entropy_change_j_k)
            .abs()
            < 1e-12
    );
    assert!(result.balance.entropy_production_j_k >= -1e-10);
    assert!(result.balance.minimum_step_entropy_production_j_k >= -1e-10);
    assert!(result.balance.maximum_estimated_local_error <= 1.0);
    assert!(std::ptr::eq(result.body.model(), original.model()));
    assert_eq!(result.body.amount(), original.amount());
    assert_eq!(result.body.mass(), original.mass());
    match original.constraint() {
        GasHeatConstraint::ConstantVolume(v) => {
            assert_eq!(result.body.volume(), v);
            assert_eq!(result.report.boundary_work_j, 0.0);
        }
        GasHeatConstraint::ConstantPressure(p) => {
            assert_eq!(result.body.pressure(), p);
            let work = original.amount().value()
                * R
                * (result.body.temperature().value() - original.temperature().value());
            assert!((result.report.boundary_work_j - work).abs() < 1e-10);
        }
    }
}

#[test]
fn finite_gas_couples_both_mechanical_constraints_to_the_reacting_parcel() {
    let model = coolant(false, 0.0);
    for pressure in [false, true] {
        let cg = if pressure { 4.0 } else { 3.0 } * R;
        for partner_pressure in [false, true] {
            let cb = 2.0 * if partner_pressure { 4.5 } else { 3.5 } * R;
            for method in [
                IsochoricKineticsMethod::BackwardEuler,
                IsochoricKineticsMethod::Sdirk2,
            ] {
                for (gas_initial, body_initial, binding, reacting) in
                    [(800.0, 300.0, 512.0, true), (300.0, 900.0, 0.0, false)]
                {
                    let original = body(
                        &model,
                        partner_pressure,
                        body_initial,
                        2.0,
                        40.0,
                        [200.0, 1500.0],
                    );
                    let result = run(
                        pressure,
                        method,
                        gas_initial,
                        0.3,
                        original,
                        binding,
                        reacting,
                    );
                    // Independent constant-C two-inventory solution with an
                    // analytic first-order reaction population and heat release.
                    let lambda = 40.0 * (1.0 / cg + 1.0 / cb);
                    let remaining = (-0.6_f64).exp();
                    let decay = (-0.3 * lambda).exp();
                    let difference = (gas_initial - body_initial) * decay
                        + binding * R * 2.0 / cg * (remaining - decay) / (lambda - 2.0);
                    let weighted =
                        cg * gas_initial + cb * body_initial + binding * R * (1.0 - remaining);
                    let expected_gas = (weighted + cb * difference) / (cg + cb);
                    let expected_body = (weighted - cg * difference) / (cg + cb);
                    let tolerance = if method == IsochoricKineticsMethod::BackwardEuler {
                        0.16
                    } else {
                        0.003
                    };
                    assert!(
                        (result.gas_t - expected_gas).abs() < tolerance,
                        "{method:?}, p={pressure}/{partner_pressure}: {} != {expected_gas}",
                        result.gas_t
                    );
                    assert!((result.body.temperature().value() - expected_body).abs() < tolerance);
                    assert!(
                        (result.report.entropy_change_j_k
                            - cb * (result.body.temperature().value() / body_initial).ln())
                        .abs()
                            < 2e-11
                    );
                    if reacting {
                        assert!((result.a - remaining).abs() < 3e-4);
                    }
                    balances(&result, original, gas_initial);
                }
            }
        }
    }
}

fn body_energy(pressure: bool, initial: f64, terminal: f64) -> f64 {
    let a = if pressure { 4.5 } else { 3.5 };
    2.0 * R * (terminal - initial) * (a + 0.5 * 0.0035 * (terminal + initial))
}

#[test]
fn finite_gas_keeps_original_variable_cp_radiation_and_analytic_feedback() {
    let model = coolant(true, 0.0);
    let controls = thermal_controls(IsochoricKineticsMethod::Sdirk2);
    for pressure in [false, true] {
        let cg = if pressure { 4.0 } else { 3.0 } * R;
        for partner_pressure in [false, true] {
            let original = body(&model, partner_pressure, 350.0, 2.0, 10.0, [250.0, 1500.0])
                .with_gray_radiation(Area::new(0.2), 0.8)
                .unwrap();
            let gas_t = |time: f64, body_t: f64| {
                1000.0
                    + (512.0 * R * (1.0 - (-2.0 * time).exp())
                        - body_energy(partner_pressure, 350.0, body_t))
                        / cg
            };
            let derivative = |time, body_t: f64| {
                let gas_t = gas_t(time, body_t);
                (10.0 * (gas_t - body_t) + SIGMA * 0.16 * (gas_t.powi(4) - body_t.powi(4)))
                    / (2.0 * R * (if partner_pressure { 4.5 } else { 3.5 } + 0.0035 * body_t))
            };
            let integrate = |steps: usize| {
                let mut t = 350.0;
                let dt = 0.3 / steps as f64;
                for index in 0..steps {
                    let time = index as f64 * dt;
                    let a = derivative(time, t);
                    let b = derivative(time + 0.5 * dt, t + 0.5 * dt * a);
                    let c = derivative(time + 0.5 * dt, t + 0.5 * dt * b);
                    let d = derivative(time + dt, t + dt * c);
                    t += dt * (a + 2.0 * b + 2.0 * c + d) / 6.0;
                }
                t
            };
            let expected_body = integrate(4096);
            assert!((expected_body - integrate(8192)).abs() < 1e-9);
            let result = run(
                pressure,
                IsochoricKineticsMethod::Sdirk2,
                1000.0,
                0.3,
                original,
                512.0,
                true,
            );
            assert!((result.body.temperature().value() - expected_body).abs() < 0.003);
            assert!((result.gas_t - gas_t(0.3, expected_body)).abs() < 0.004);
            balances(&result, original, 1000.0);

            let mut fixed = Sum::default();
            fixed.add(1500.0);
            let stage =
                HeatStage::implicit(fixed, Some(HeatExchange::Gas(&original)), 2.0, 0.2, 0.0)
                    .unwrap();
            let (heat, feedback) = stage.heat(850.0, controls, &mut || false).unwrap();
            let (mut lo, mut hi) = (250.0, 1500.0);
            for _ in 0..80 {
                let t = 0.5 * (lo + hi);
                let residual = body_energy(partner_pressure, 350.0, t)
                    + 3000.0
                    + 0.2 * (10.0 * (t - 850.0) + SIGMA * 0.16 * (t.powi(4) - 850.0_f64.powi(4)));
                if residual < 0.0 {
                    lo = t;
                } else {
                    hi = t;
                }
            }
            let body_t = 0.5 * (lo + hi);
            assert!(
                (heat.value() * 2.0 + body_energy(partner_pressure, 350.0, body_t)).abs() < 2e-8
            );
            let capacity = 2.0 * R * (if partner_pressure { 4.5 } else { 3.5 } + 0.0035 * body_t);
            let expected_feedback = 0.1 * (10.0 + 4.0 * SIGMA * 0.16 * 850.0_f64.powi(3))
                / (1.0 + 0.2 * (10.0 + 4.0 * SIGMA * 0.16 * body_t.powi(3)) / capacity);
            assert!((feedback - expected_feedback).abs() < 1e-9);
            let plus = stage
                .heat(850.001, controls, &mut || false)
                .unwrap()
                .0
                .value();
            let minus = stage
                .heat(849.999, controls, &mut || false)
                .unwrap()
                .0
                .value();
            assert!((feedback + (plus - minus) / 0.002).abs() < 3e-7 * feedback.max(1.0));
        }
    }
}

#[test]
fn finite_gas_source_energy_and_entropy_datums_cancel_before_evolution() {
    for pressure in [false, true] {
        for partner_pressure in [false, true] {
            let mut baseline: Option<[f64; 4]> = None;
            for datum in [0.0, 2.0_f64.powi(45), -2.0_f64.powi(45)] {
                let model = coolant(true, datum);
                let original = body(&model, partner_pressure, 400.0, 2.0, 20.0, [300.0, 700.0]);
                // Reactor begins outside the partner's declared T domain;
                // only the actual partner trajectory must stay inside it.
                let result = run(
                    pressure,
                    IsochoricKineticsMethod::Sdirk2,
                    1200.0,
                    0.2,
                    original,
                    512.0,
                    true,
                );
                let values = [
                    result.gas_t,
                    result.body.temperature().value(),
                    result.report.heat_input_j,
                    result.report.entropy_change_j_k,
                ];
                if let Some(expected) = baseline {
                    assert_eq!(values.map(f64::to_bits), expected.map(f64::to_bits));
                } else {
                    baseline = Some(values);
                }
                assert!(
                    result.body.temperature().value() > 400.0
                        && result.body.temperature().value() < 700.0
                );
                balances(&result, original, 1200.0);
            }
        }
    }
}

#[test]
fn tiny_finite_gas_inventory_retains_temperature_accuracy_below_global_joule_budget() {
    let model = coolant(false, 0.0);
    for partner_pressure in [false, true] {
        let capacity = 1e-13 * if partner_pressure { 4.5 } else { 3.5 } * R;
        let original = body(
            &model,
            partner_pressure,
            300.0,
            1e-13,
            2.0 * capacity,
            [200.0, 1200.0],
        );
        // The entire partner's admitted energy span is below the reactor's
        // absolute J tolerance. That does not permit a source-endpoint snap
        // or hiding its temporal error from the existing adaptive controller.
        assert!(capacity * 1000.0 < controls(2e-7).energy_tolerance_j);
        let result = run(
            false,
            IsochoricKineticsMethod::Sdirk2,
            800.0,
            0.3,
            original,
            0.0,
            false,
        );
        let expected = 800.0 - 500.0 * (-0.6_f64).exp();
        assert!((result.body.temperature().value() - expected).abs() < 0.003);
        assert!((result.gas_t - 800.0).abs() < 1e-7);
        assert!(result.report.heat_input_j > 0.0);
        balances(&result, original, 800.0);
    }
}

#[test]
fn finite_gas_source_exhaustion_cancellation_and_insulated_replay_are_atomic() {
    let model = coolant(false, 0.0);
    let original = body(&model, true, 400.0, 2.0, 40.0, [300.0, 700.0])
        .with_gray_radiation(Area::new(0.2), 0.8)
        .unwrap();
    let controls = thermal_controls(IsochoricKineticsMethod::Sdirk2);
    let mut reactor =
        IsochoricReactor::new(thermal_mechanism(512.0, 2.0, 0.0), controls, || false).unwrap();
    let feed = [Amount::new(1.0), Amount::new(0.0)];
    let expected = reactor
        .advance_with_finite_gas_body(
            &feed,
            Volume::new(1.0),
            Temperature::new(1000.0),
            Time::new(0.2),
            original,
            || false,
        )
        .unwrap();
    let expected_gas = expected.gas().to_owned(|| false).unwrap();
    let expected_body = expected.body();
    let mut polls = 0;
    assert!(matches!(
        reactor.advance_with_finite_gas_body(
            &feed,
            Volume::new(1.0),
            Temperature::new(1000.0),
            Time::new(0.2),
            original,
            || {
                polls += 1;
                polls > 100
            }
        ),
        Err(IsochoricKineticsError::Cancelled)
    ));
    let replay = reactor
        .advance_with_finite_gas_body(
            &feed,
            Volume::new(1.0),
            Temperature::new(1000.0),
            Time::new(0.2),
            original,
            || false,
        )
        .unwrap();
    assert_eq!(expected_gas, replay.gas().to_owned(|| false).unwrap());
    assert_eq!(expected_body, replay.body());
    assert_eq!(original.temperature(), Temperature::new(400.0));
    let insulated = body(&model, false, 400.0, 2.0, 0.0, [300.0, 700.0]);
    let with_body = reactor
        .advance_with_finite_gas_body(
            &feed,
            Volume::new(1.0),
            Temperature::new(1000.0),
            Time::new(0.2),
            insulated,
            || false,
        )
        .unwrap();
    assert_eq!(with_body.body(), insulated);
    let with_body = with_body.gas().to_owned(|| false).unwrap();
    let adiabatic = reactor
        .advance_adiabatic(
            &feed,
            Volume::new(1.0),
            Temperature::new(1000.0),
            Time::new(0.2),
            || false,
        )
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    assert_eq!(with_body, adiabatic);
    let edge = body(&model, true, 500.0, 2.0, 40.0, [300.0, 500.0]);
    assert!(
        reactor
            .advance_with_finite_gas_body(
                &feed,
                Volume::new(1.0),
                Temperature::new(1000.0),
                Time::new(0.02),
                edge,
                || false
            )
            .is_err()
    );
    assert_eq!(edge.temperature(), Temperature::new(500.0));
    let refused = FiniteGasHeatBody::new(
        &model,
        Amount::new(2.0),
        GasHeatConstraint::ConstantVolume(Volume::new(1.0)),
        Temperature::new(300.0),
        [Temperature::new(100.0), Temperature::new(500.0)],
        ThermalConductance::new(1.0),
        || false,
    );
    assert!(refused.is_err());
}
