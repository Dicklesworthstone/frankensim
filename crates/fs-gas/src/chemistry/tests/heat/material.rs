use super::*;
use crate::chemistry::heat::material::FiniteHeatBody;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};

fn curve(datum: f64, latent: bool) -> EquilibriumEnthalpyPhaseCurve {
    // Authored G1 calorimeter, not a measured material or source-validation deck.
    let identity = source("manufactured-calorimeter", 1.0, 4.0, 0.0, 0.0).regions()[0]
        .card()
        .content_hash();
    let knot = |h, t, f, rho| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: datum + h,
        temperature_k: t,
        liquid_mass_fraction: f,
        bulk_density_kg_m3: rho,
    };
    if latent {
        EquilibriumEnthalpyPhaseCurve::try_new(
            identity,
            vec![
                knot(0.0, 200.0, 0.0, 1000.0),
                knot(20_000.0, 400.0, 0.0, 950.0),
                knot(40_000.0, 400.0, 1.0, 800.0),
                knot(100_000.0, 1000.0, 1.0, 700.0),
            ],
        )
        .unwrap()
    } else {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            identity,
            SolidLiquidPhase::Solid,
            vec![
                knot(0.0, 200.0, 0.0, 1000.0),
                knot(20_000.0, 1200.0, 0.0, 500.0),
            ],
        )
        .unwrap()
    }
}

fn body<'a>(curve: &'a EquilibriumEnthalpyPhaseCurve, h: f64, g: f64) -> FiniteHeatBody<'a> {
    FiniteHeatBody::new(
        curve,
        Mass::new(1.0),
        Pressure::new(2e5),
        h,
        ThermalConductance::new(g),
    )
    .unwrap()
}

#[derive(Clone, Copy)]
struct Coupled {
    gas_t: f64,
    body_t: f64,
    body_h: f64,
    body_f: f64,
    a: f64,
    gas_u: f64,
    gas_work: f64,
    body_u_change: f64,
    body_work: f64,
    body_h_change: f64,
    body_entropy: f64,
    balance: IsochoricKineticsBalance,
}

fn run(
    pressure: bool,
    method: IsochoricKineticsMethod,
    initial_t: f64,
    duration: f64,
    original: FiniteHeatBody<'_>,
    binding: f64,
    reacting: bool,
) -> Coupled {
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
            .advance_with_finite_body(
                &feed,
                Pressure::new(R * initial_t),
                Temperature::new(initial_t),
                Time::new(duration),
                original,
                || false,
            )
            .unwrap();
        let gas = result.gas();
        let body = result.body().state();
        let wall = result.body_report();
        let out = Coupled {
            gas_t: gas.temperature().value(),
            body_t: body.temperature_k(),
            body_h: body.specific_enthalpy_j_kg(),
            body_f: body.liquid_mass_fraction(),
            a: gas.amounts()[0].value(),
            gas_u: gas.internal_energy().value(),
            gas_work: gas.boundary_work().value(),
            body_u_change: wall.internal_energy_change_j,
            body_work: wall.boundary_work_j,
            body_h_change: wall.enthalpy_change_j,
            body_entropy: result.body_entropy_change_j_k(),
            balance: gas.balance(),
        };
        assert!(std::ptr::eq(result.body().curve(), original.curve()));
        assert_eq!(bytes, reactor.workspace_bytes());
        out
    } else {
        let mut reactor = IsochoricReactor::new(mechanism, controls, || false).unwrap();
        let bytes = reactor.workspace_bytes();
        let result = reactor
            .advance_with_finite_body(
                &feed,
                Volume::new(1.0),
                Temperature::new(initial_t),
                Time::new(duration),
                original,
                || false,
            )
            .unwrap();
        let gas = result.gas();
        let body = result.body().state();
        let wall = result.body_report();
        let out = Coupled {
            gas_t: gas.temperature().value(),
            body_t: body.temperature_k(),
            body_h: body.specific_enthalpy_j_kg(),
            body_f: body.liquid_mass_fraction(),
            a: gas.amounts()[0].value(),
            gas_u: gas.internal_energy().value(),
            gas_work: 0.0,
            body_u_change: wall.internal_energy_change_j,
            body_work: wall.boundary_work_j,
            body_h_change: wall.enthalpy_change_j,
            body_entropy: result.body_entropy_change_j_k(),
            balance: gas.balance(),
        };
        assert!(std::ptr::eq(result.body().curve(), original.curve()));
        assert_eq!(bytes, reactor.workspace_bytes());
        out
    }
}

fn balances(result: Coupled, initial_u: f64) {
    assert!(result.balance.energy_residual_j.abs() < 1e-8);
    assert!((result.body_h_change + result.balance.heat_input_j).abs() < 1e-8);
    let first_law =
        result.gas_u - initial_u + result.body_u_change + result.gas_work + result.body_work;
    assert!(
        first_law.abs() < 4e-8,
        "combined first-law residual {first_law}"
    );
    assert_eq!(result.balance.entropy_transfer_j_k, -result.body_entropy);
    assert!(
        (result.balance.entropy_production_j_k
            - (result.balance.entropy_change_j_k + result.body_entropy))
            .abs()
            < 1e-12
    );
    assert!(result.balance.entropy_production_j_k >= -1e-10);
    assert!(result.balance.minimum_step_entropy_production_j_k >= -1e-10);
    assert!(result.balance.maximum_estimated_local_error <= 1.0);
}

#[test]
fn finite_sensible_body_and_reaction_match_independent_two_body_solution() {
    let chart = curve(0.0, false);
    let (duration, rate, g, cb) = (0.3, 2.0, 40.0, 20.0);
    for pressure in [false, true] {
        let cg = if pressure { 4.0 } else { 3.0 } * R;
        for method in [
            IsochoricKineticsMethod::BackwardEuler,
            IsochoricKineticsMethod::Sdirk2,
        ] {
            for (tg0, tb0, binding, reacting) in [
                (800.0, 300.0, 0.0, false),
                (300.0, 900.0, 0.0, false),
                (700.0, 300.0, 512.0, true),
            ] {
                let original = body(&chart, cb * (tb0 - 200.0), g);
                let result = run(pressure, method, tg0, duration, original, binding, reacting);
                let lambda = g * (1.0 / cg + 1.0 / cb);
                let chemical = binding * R;
                let remaining = (-rate * duration).exp();
                let decay = (-lambda * duration).exp();
                let difference = (tg0 - tb0) * decay
                    + chemical * rate / cg * (remaining - decay) / (lambda - rate);
                let weighted = cg * tg0 + cb * tb0 + chemical * (1.0 - remaining);
                let expected_tg = (weighted + cb * difference) / (cg + cb);
                let expected_tb = (weighted - cg * difference) / (cg + cb);
                let tolerance = if method == IsochoricKineticsMethod::BackwardEuler {
                    0.12
                } else {
                    0.002
                };
                assert!(
                    (result.gas_t - expected_tg).abs() < tolerance,
                    "{pressure} {method:?}: {} != {expected_tg}",
                    result.gas_t
                );
                assert!((result.body_t - expected_tb).abs() < tolerance);
                if reacting {
                    assert!((result.a - remaining).abs() < 3e-4);
                }
                assert!((result.body_entropy - cb * (result.body_t / tb0).ln()).abs() < 1e-11);
                balances(result, 3.0 * R * tg0);
            }
        }
    }
}

#[test]
fn latent_body_heats_cools_and_preserves_original_enthalpy_datum() {
    for pressure in [false, true] {
        let cg = if pressure { 4.0 } else { 3.0 } * R;
        for method in [
            IsochoricKineticsMethod::BackwardEuler,
            IsochoricKineticsMethod::Sdirk2,
        ] {
            for initial_t in [300.0, 800.0] {
                let mut reference: Option<Coupled> = None;
                for datum in [0.0, 1_048_576.0, -1_048_576.0] {
                    let chart = curve(datum, true);
                    let original = body(&chart, datum + 25_000.0, 40.0);
                    let result = run(pressure, method, initial_t, 0.3, original, 0.0, false);
                    let expected_tg = 400.0 + (initial_t - 400.0) * (-40.0 / cg * 0.3).exp();
                    let expected_h = 25_000.0 - cg * (expected_tg - initial_t);
                    assert_eq!(result.body_t, 400.0);
                    assert!((result.gas_t - expected_tg).abs() < 0.12);
                    assert!((result.body_h - datum - expected_h).abs() < 4.0);
                    assert!(
                        (result.body_f - (result.body_h - datum - 20_000.0) / 20_000.0).abs()
                            < 1e-14
                    );
                    balances(result, 3.0 * R * initial_t);
                    if let Some(first) = reference {
                        assert!((result.gas_t - first.gas_t).abs() < 1e-7);
                        assert!((result.body_h - datum - first.body_h).abs() < 1e-6);
                        assert!(
                            (result.balance.heat_input_j - first.balance.heat_input_j).abs() < 1e-6
                        );
                    } else {
                        reference = Some(result);
                    }
                }
            }
        }
    }
}

#[test]
fn sensible_to_latent_crossing_follows_heat_before_melting() {
    let chart = curve(0.0, true);
    for pressure in [false, true] {
        let cg = if pressure { 4.0 } else { 3.0 } * R;
        let (tg0, tb0, cb, g, duration) = (1000.0, 399.0, 100.0, 40.0, 0.2);
        let equilibrium = (cg * tg0 + cb * tb0) / (cg + cb);
        let ratio = (equilibrium - 400.0) * (cg + cb) / (cg * (tg0 - tb0));
        let crossing = -ratio.ln() / (g * (1.0 / cg + 1.0 / cb));
        let gas_at_crossing = tg0 - cb * (400.0 - tb0) / cg;
        let expected = 400.0 + (gas_at_crossing - 400.0) * (-g / cg * (duration - crossing)).exp();
        let result = run(
            pressure,
            IsochoricKineticsMethod::Sdirk2,
            tg0,
            duration,
            body(&chart, 19_900.0, g),
            0.0,
            false,
        );
        assert_eq!(result.body_t, 400.0);
        assert!(result.body_f > 0.0 && result.body_f < 1.0);
        assert!((result.gas_t - expected).abs() < 0.002);
        balances(result, 3.0 * R * tg0);
    }
}

#[test]
fn finite_body_refusals_cancel_replay_and_insulated_parity() {
    let chart = curve(0.0, true);
    let original = body(&chart, 25_000.0, 40.0);
    let controls = thermal_controls(IsochoricKineticsMethod::Sdirk2);
    let feed = [Amount::new(1.0), Amount::new(0.0)];
    let mut reactor =
        IsochoricReactor::new(thermal_mechanism(512.0, 2.0, 0.0), controls, || false).unwrap();
    let expected = reactor
        .advance_with_finite_body(
            &feed,
            Volume::new(1.0),
            Temperature::new(700.0),
            Time::new(0.3),
            original,
            || false,
        )
        .unwrap();
    let expected_gas = expected.gas().to_owned(|| false).unwrap();
    let expected_body = expected.body();
    let mut polls = 0;
    let refused = reactor.advance_with_finite_body(
        &feed,
        Volume::new(1.0),
        Temperature::new(700.0),
        Time::new(0.3),
        original,
        || {
            polls += 1;
            polls > 50
        },
    );
    assert!(matches!(refused, Err(IsochoricKineticsError::Cancelled)));
    let replay = reactor
        .advance_with_finite_body(
            &feed,
            Volume::new(1.0),
            Temperature::new(700.0),
            Time::new(0.3),
            original,
            || false,
        )
        .unwrap();
    assert_eq!(expected_gas, replay.gas().to_owned(|| false).unwrap());
    assert_eq!(expected_body, replay.body());
    assert_eq!(original.state().specific_enthalpy_j_kg(), 25_000.0);

    let insulated = body(&chart, 25_000.0, 0.0);
    let result = reactor
        .advance_with_finite_body(
            &feed,
            Volume::new(1.0),
            Temperature::new(700.0),
            Time::new(0.3),
            insulated,
            || false,
        )
        .unwrap();
    assert_eq!(result.body(), insulated);
    let with_body = result.gas().to_owned(|| false).unwrap();
    let adiabatic = reactor
        .advance_adiabatic(
            &feed,
            Volume::new(1.0),
            Temperature::new(700.0),
            Time::new(0.3),
            || false,
        )
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    assert_eq!(with_body, adiabatic);
    let mut limited = IsochoricReactor::new(
        thermal_mechanism(0.0, 2.0, 0.0),
        IsochoricKineticsSolve {
            maximum_steps: 100,
            ..controls
        },
        || false,
    )
    .unwrap();
    let edge = body(&chart, 100_000.0, 40.0);
    assert!(
        limited
            .advance_with_finite_body(
                &[Amount::new(0.0), Amount::new(1.0)],
                Volume::new(1.0),
                Temperature::new(1500.0),
                Time::new(0.1),
                edge,
                || false
            )
            .is_err()
    );
    assert_eq!(edge.state().specific_enthalpy_j_kg(), 100_000.0);
    assert!(
        FiniteHeatBody::new(
            &chart,
            Mass::new(0.0),
            Pressure::new(1e5),
            0.0,
            ThermalConductance::new(1.0)
        )
        .is_err()
    );
    assert!(
        FiniteHeatBody::new(
            &chart,
            Mass::new(1.0),
            Pressure::new(-1.0),
            0.0,
            ThermalConductance::new(1.0)
        )
        .is_err()
    );
}

#[test]
fn finite_body_constrained_jacobian_retains_both_caloric_feedbacks() {
    let mechanism = mechanism(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 5.0, -2048.0 * R, 0.0),
        ],
        vec![1, 1],
        vec![
            ElementaryReactionV1::new(
                ReactionId::new("activated-finite-body").unwrap(),
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
            .irreversible(),
        ],
    );
    for pressure in [false, true] {
        for latent in [false, true] {
            let chart = curve(0.0, latent);
            let original = body(&chart, if latent { 25_000.0 } else { 3000.0 }, 1000.0);
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
            let (_, derivative) = reactor
                .finite_body_rate_linearization(&[0.08], 500.0, original, 0.1)
                .unwrap();
            let (plus, _) = reactor
                .finite_body_rate_linearization(&[0.080001], 500.0, original, 0.1)
                .unwrap();
            let (minus, _) = reactor
                .finite_body_rate_linearization(&[0.079999], 500.0, original, 0.1)
                .unwrap();
            let numeric = (plus[0] - minus[0]) / 2e-6;
            assert!(
                (numeric - derivative[0]).abs() < 2e-7 * numeric.abs().max(1.0),
                "p={pressure}, latent={latent}: {numeric} != {}",
                derivative[0]
            );
        }
    }
}
