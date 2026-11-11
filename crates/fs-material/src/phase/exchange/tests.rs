//! Independent manufactured calorics and exchange times; no measured material claim.

#[path = "tests/driven.rs"]
mod driven;

use super::*;
use crate::phase::{EnthalpyPhaseKnot, SolidLiquidPhase};
use fs_blake3::ContentHash;

fn source_curve(datum: f64, split: bool) -> EquilibriumEnthalpyPhaseCurve {
    let knot = |h, t, liquid, rho| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: datum + h,
        temperature_k: t,
        liquid_mass_fraction: liquid,
        bulk_density_kg_m3: rho,
    };
    let mut knots = vec![knot(0.0, 300.0, 0.0, 1000.0)];
    let mut capacities = Vec::new();
    let mut previous_cp = 1000.0;
    let solid: &[f64] = if split { &[350.0, 400.0] } else { &[400.0] };
    for &t in solid {
        let x = t - 300.0;
        let cp = 1000.0 + 10.0 * x;
        knots.push(knot(1000.0 * x + 5.0 * x * x, t, 0.0, 1000.0 - x));
        capacities.push(Some([previous_cp, cp]));
        previous_cp = cp;
    }
    knots.push(knot(170_000.0, 400.0, 1.0, 800.0));
    capacities.push(None);
    previous_cp = 2000.0;
    let liquid: &[f64] = if split { &[450.0, 500.0] } else { &[500.0] };
    for &t in liquid {
        let x = t - 400.0;
        let cp = 2000.0 + 20.0 * x;
        knots.push(knot(
            170_000.0 + 2000.0 * x + 10.0 * x * x,
            t,
            1.0,
            800.0 - x,
        ));
        capacities.push(Some([previous_cp, cp]));
        previous_cp = cp;
    }
    let mut curve = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0xa7; 32]), knots).unwrap();
    curve.retain_sensible_density_law();
    curve.retain_sensible_heat_capacities(capacities).unwrap();
    curve
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

// Independent primitive for Cp(T)=intercept+slope*T, used to choose a target
// temperature and derive its exposure, rather than repeat the production inverse.
fn exposure(from: f64, to: f64, bath: f64, intercept: f64, slope: f64) -> f64 {
    (intercept + slope * bath) * ((bath - from) / (bath - to)).ln() - slope * (to - from)
}

#[test]
fn g1_source_cp_exchange_crosses_sensible_latent_and_liquid_in_both_directions() {
    let curve = source_curve(0.0, false);
    let heating = exposure(320.0, 400.0, 550.0, -2000.0, 10.0)
        + 20_000.0 / 150.0
        + exposure(400.0, 460.0, 550.0, -6000.0, 20.0);
    let cooling = exposure(500.0, 400.0, 250.0, -6000.0, 20.0)
        + 20_000.0 / 150.0
        + exposure(400.0, 330.0, 250.0, -2000.0, 10.0);
    for (h0, target_h, target_t, bath, x) in [
        (22_000.0, 326_000.0, 460.0, 550.0, heating),
        (470_000.0, 34_500.0, 330.0, 250.0, cooling),
    ] {
        let initial = curve.state_at_specific_enthalpy(h0).unwrap();
        let result = curve.linear_heat_exchange(initial, bath, x).unwrap();
        near(result.state.temperature_k(), target_t, 2e-10);
        near(result.state.specific_enthalpy_j_kg(), target_h, 2e-7);
        near(result.specific_heat_j_kg, target_h - h0, 2e-7);
        near(result.mean_temperature_k, bath - (target_h - h0) / x, 2e-10);
        near(result.heat_residual_j_kg, 0.0, 2e-9);
        assert_eq!(result.initial, initial);
        assert_eq!(result.state.phase_curve_identity(), curve.identity());
        assert_eq!(
            result.state.material_card_identity(),
            curve.material_card_identity()
        );
        // Original density is evaluated at the accepted state, not fixed at
        // inlet density or averaged over the thermal trajectory.
        near(
            result.state.bulk_density_kg_m3(),
            if target_t > 400.0 {
                1200.0 - target_t
            } else {
                1300.0 - target_t
            },
            2e-10,
        );
    }
}

#[test]
fn g1_latent_residence_uses_actual_driving_temperature_and_retains_equilibrium() {
    let curve = source_curve(0.0, false);
    let initial = curve.state_at_specific_enthalpy(155_000.0).unwrap();
    for (bath, x, h, fraction) in [
        (500.0, 50.0, 160_000.0, 0.5),
        (300.0, 50.0, 150_000.0, 0.0),
        (500.0, 150.0, 170_000.0, 1.0),
    ] {
        let result = curve.linear_heat_exchange(initial, bath, x).unwrap();
        assert_eq!(result.state.temperature_k(), 400.0);
        assert_eq!(result.mean_temperature_k, 400.0);
        assert_eq!(result.state.specific_enthalpy_j_kg(), h);
        assert_eq!(result.state.liquid_mass_fraction(), fraction);
        assert_eq!(result.specific_heat_j_kg, (bath - 400.0) * x);
        assert_eq!(result.heat_residual_j_kg, 0.0);
        near(
            result.state.bulk_density_kg_m3(),
            ((1.0 - fraction) / 900.0 + fraction / 800.0).recip(),
            2e-12,
        );
    }
    for (bath, x) in [(400.0, 1e250), (600.0, 0.0)] {
        let result = curve.linear_heat_exchange(initial, bath, x).unwrap();
        assert_eq!(result.state, initial);
        assert_eq!(result.specific_heat_j_kg, 0.0);
        assert_eq!(result.mean_temperature_k, 400.0);
    }
}

#[test]
fn g3_chart_subdivision_and_exchange_composition_preserve_original_calorics() {
    for datum in [0.0, 2.0_f64.powi(30)] {
        let original = source_curve(datum, false);
        let split = source_curve(datum, true);
        let h = datum + 22_000.0;
        let initial = original.state_at_specific_enthalpy(h).unwrap();
        let direct = original
            .linear_heat_exchange(initial, 550.0, 1200.0)
            .unwrap();
        let subdivided = split
            .linear_heat_exchange(split.state_at_specific_enthalpy(h).unwrap(), 550.0, 1200.0)
            .unwrap();
        near(
            direct.state.temperature_k(),
            subdivided.state.temperature_k(),
            2e-9,
        );
        near(
            direct.specific_heat_j_kg,
            subdivided.specific_heat_j_kg,
            2e-6,
        );
        let first = original
            .linear_heat_exchange(initial, 550.0, 600.0)
            .unwrap();
        let second = original
            .linear_heat_exchange(first.state, 550.0, 600.0)
            .unwrap();
        near(
            second.state.temperature_k(),
            direct.state.temperature_k(),
            2e-9,
        );
        near(
            first.specific_heat_j_kg + second.specific_heat_j_kg,
            direct.specific_heat_j_kg,
            2e-6,
        );
        near(
            0.5 * first.mean_temperature_k + 0.5 * second.mean_temperature_k,
            direct.mean_temperature_k,
            2e-9,
        );
        assert!(direct.heat_residual_j_kg.abs() < 5e-7);
    }
}

#[test]
fn g1_authored_mushy_intervals_retain_their_total_enthalpy_law() {
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xc7; 32]),
        [0.0, 0.5, 1.0]
            .into_iter()
            .map(|fraction| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 120_000.0 * fraction,
                temperature_k: 300.0 + 200.0 * fraction,
                liquid_mass_fraction: fraction,
                bulk_density_kg_m3: 1000.0 - 200.0 * fraction,
            })
            .collect(),
    )
    .unwrap();
    let initial = curve.state_at_specific_enthalpy(12_000.0).unwrap();
    let x = 600.0 * (280.0_f64 / 120.0).ln();
    let response = curve
        .linear_heat_exchange_response(initial, 600.0, x)
        .unwrap();
    let result = response.exchange;
    near(response.jacobian.inlet_enthalpy_gain, 3.0 / 7.0, 3e-14);
    near(
        response.jacobian.reservoir_enthalpy_derivative_j_kg_k,
        2400.0 / 7.0,
        3e-10,
    );
    near(
        response.jacobian.exposure_enthalpy_derivative_k,
        120.0,
        3e-11,
    );
    near(
        response.jacobian.mean_log_exposure_derivative_k,
        480.0 - 600.0 + 96_000.0 / x,
        3e-10,
    );
    near(result.state.temperature_k(), 480.0, 2e-10);
    near(result.state.specific_enthalpy_j_kg(), 108_000.0, 2e-8);
    near(result.state.liquid_mass_fraction(), 0.9, 2e-13);
    near(result.specific_heat_j_kg, 96_000.0, 2e-8);
}

#[test]
fn g1_weak_exchange_retains_relative_heat_and_reports_absolute_state_rounding() {
    let curve = source_curve(0.0, false);
    let initial = curve.state_at_specific_enthalpy(62_500.0).unwrap();
    let x = 1e-8;
    let result = curve.linear_heat_exchange(initial, 550.0, x).unwrap();
    // Independent local expansion of dh/dX=550-T, dT/dh=1/1500.
    let heat = 200.0 * x - 0.5 * (200.0 / 1500.0) * x * x;
    near(result.specific_heat_j_kg, heat, 1e-19);
    near(
        result.mean_temperature_k,
        350.0 + 0.5 * (200.0 / 1500.0) * x,
        2e-13,
    );
    assert_ne!(result.heat_residual_j_kg, 0.0);
    assert_eq!(
        result.heat_residual_j_kg,
        result.specific_heat_j_kg - result.represented_enthalpy_change_j_kg
    );
}

#[test]
fn g0_domain_identity_quantization_and_late_cancellation_refuse_transactionally() {
    let curve = source_curve(0.0, false);
    let before = curve.clone();
    let initial = curve.state_at_specific_enthalpy(22_000.0).unwrap();
    for (bath, x) in [
        (0.0, 1.0),
        (f64::NAN, 1.0),
        (550.0, -1.0),
        (550.0, f64::INFINITY),
        (550.0, f64::NAN),
        (600.0, 1e10),
        (200.0, 1e10),
    ] {
        assert!(curve.linear_heat_exchange(initial, bath, x).is_err());
    }
    let other = source_curve(1.0, false);
    assert!(matches!(
        curve.linear_heat_exchange(
            other.state_at_specific_enthalpy(22_001.0).unwrap(),
            550.0,
            1.0
        ),
        Err(PhaseHeatExchangeError::Phase(
            PhaseStateError::StateCurveMismatch
        ))
    ));
    let shifted = source_curve(2.0_f64.powi(59), false);
    let start = shifted
        .state_at_specific_enthalpy(2.0_f64.powi(59))
        .unwrap();
    assert!(matches!(
        shifted.linear_heat_exchange(start, 550.0, 1e-6),
        Err(PhaseHeatExchangeError::Phase(
            PhaseStateError::LostSpecificEnergyIncrement
        ))
    ));
    let mut checks = 0;
    let refused = curve.linear_heat_exchange_with_checkpoint(initial, 550.0, 1200.0, || {
        checks += 1;
        checks == 3
    });
    assert!(matches!(refused, Err(PhaseHeatExchangeError::Cancelled)));
    assert_eq!(checks, 3);
    assert_eq!(curve, before);
    assert_eq!(initial, curve.state_at_specific_enthalpy(22_000.0).unwrap());
}

// Independent derivative of the manufactured exposure primitive with respect
// to the bath, integrated over the prescribed path rather than perturbed marches.
fn bath_moment(from: f64, to: f64, bath: f64, intercept: f64, slope: f64) -> f64 {
    (intercept + slope * bath) * (1.0 / (bath - to) - 1.0 / (bath - from))
        + slope * ((bath - to) / (bath - from)).ln()
}

#[test]
fn g1_source_exact_response_carries_latent_residence_in_both_directions() {
    let curve = source_curve(0.0, false);
    for (h0, from, to, bath, first_a, first_b, last_a, last_b, latent) in [
        (
            22_000.0, 320.0, 460.0, 550.0, -2000.0, 10.0, -6000.0, 20.0, 20_000.0,
        ),
        (
            470_000.0, 500.0, 330.0, 250.0, -6000.0, 20.0, -2000.0, 10.0, -20_000.0,
        ),
    ] {
        let x = exposure(from, 400.0, bath, first_a, first_b)
            + latent / (bath - 400.0)
            + exposure(400.0, to, bath, last_a, last_b);
        let response = curve
            .linear_heat_exchange_response(curve.state_at_specific_enthalpy(h0).unwrap(), bath, x)
            .unwrap();
        let j = response.jacobian;
        let bath_derivative = (bath - to)
            * (bath_moment(from, 400.0, bath, first_a, first_b)
                + latent / (bath - 400.0).powi(2)
                + bath_moment(400.0, to, bath, last_a, last_b));
        near(j.inlet_enthalpy_gain, (bath - to) / (bath - from), 3e-14);
        near(
            j.inlet_enthalpy_loss,
            1.0 - (bath - to) / (bath - from),
            3e-14,
        );
        near(
            j.reservoir_enthalpy_derivative_j_kg_k,
            bath_derivative,
            3e-10,
        );
        near(j.exposure_enthalpy_derivative_k, bath - to, 3e-11);
        let target_h = if to > 400.0 {
            let delta = to - 400.0;
            170_000.0 + 2000.0 * delta + 10.0 * delta * delta
        } else {
            let delta = to - 300.0;
            1000.0 * delta + 5.0 * delta * delta
        };
        near(
            j.mean_log_exposure_derivative_k,
            to - bath + (target_h - h0) / x,
            3e-10,
        );
        near(
            j.mean_reservoir_derivative,
            1.0 - bath_derivative / x,
            3e-13,
        );
        assert_eq!(
            response.exchange,
            curve
                .linear_heat_exchange(response.exchange.initial, bath, x)
                .unwrap(),
        );
    }
}

#[test]
fn g1_weak_and_rounded_equilibrium_responses_retain_nonzero_derivatives() {
    let curve = source_curve(0.0, false);
    let initial = curve.state_at_specific_enthalpy(62_500.0).unwrap(); // T=350, Cp=1500.
    let x = 1e-200;
    let weak = curve
        .linear_heat_exchange_response(initial, 350.0, x)
        .unwrap();
    assert_eq!(weak.exchange.state, initial);
    assert_eq!(weak.jacobian.inlet_enthalpy_gain, 1.0);
    near(weak.jacobian.inlet_enthalpy_loss / (x / 1500.0), 1.0, 3e-15);
    near(
        weak.jacobian.reservoir_enthalpy_derivative_j_kg_k / x,
        1.0,
        3e-15,
    );
    // Forming X times this mean response first would underflow to zero.
    near(
        weak.jacobian.mean_reservoir_derivative / (x / 3000.0),
        1.0,
        3e-15,
    );
    assert_eq!(weak.jacobian.mean_log_exposure_derivative_k, 0.0);

    let constant = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0xd7; 32]),
        SolidLiquidPhase::Solid,
        [(0.0, 300.0), (100_000.0, 800.0)]
            .into_iter()
            .map(|(h, t)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            })
            .collect(),
    )
    .unwrap();
    for x in [8000.0, 200_000.0] {
        let response = constant
            .linear_heat_exchange_response(
                constant.state_at_specific_enthalpy(0.0).unwrap(),
                600.0,
                x,
            )
            .unwrap();
        let gain = (-x / 200.0).exp();
        assert_eq!(response.exchange.state.temperature_k(), 600.0);
        assert_eq!(response.jacobian.inlet_enthalpy_gain, gain);
        assert_eq!(
            response.jacobian.exposure_enthalpy_derivative_k,
            300.0 * gain
        );
        near(
            response.jacobian.reservoir_enthalpy_derivative_j_kg_k,
            200.0,
            3e-12,
        );
        near(
            response.jacobian.mean_reservoir_derivative,
            1.0 - 200.0 / x,
            3e-15,
        );
        near(
            response.jacobian.mean_log_exposure_derivative_k,
            300.0 * ((200.0 / x) * (1.0 - gain) - gain),
            3e-12,
        );
    }
}

#[test]
fn g0_equilibrium_latent_interior_has_a_response_but_nondifferentiable_corners_do_not() {
    let curve = source_curve(0.0, true);
    let plateau = curve
        .linear_heat_exchange_response(
            curve.state_at_specific_enthalpy(160_000.0).unwrap(),
            400.0,
            50.0,
        )
        .unwrap();
    assert_eq!(plateau.jacobian.inlet_enthalpy_gain, 1.0);
    assert_eq!(plateau.jacobian.inlet_enthalpy_loss, 0.0);
    assert_eq!(plateau.jacobian.reservoir_enthalpy_derivative_j_kg_k, 50.0);
    assert_eq!(plateau.jacobian.exposure_enthalpy_derivative_k, 0.0);
    assert_eq!(plateau.jacobian.mean_reservoir_derivative, 0.0);
    assert_eq!(plateau.jacobian.mean_log_exposure_derivative_k, 0.0);
    // An ordinary source subdivision with equal Cp on both sides is smooth.
    assert!(
        curve
            .linear_heat_exchange_response(
                curve.state_at_specific_enthalpy(62_500.0).unwrap(),
                350.0,
                50.0,
            )
            .is_ok()
    );
    for h in [150_000.0, 170_000.0] {
        let initial = curve.state_at_specific_enthalpy(h).unwrap();
        assert_eq!(
            curve
                .linear_heat_exchange(initial, 400.0, 50.0)
                .unwrap()
                .state,
            initial
        );
        assert!(matches!(
            curve.linear_heat_exchange_response(initial, 400.0, 50.0),
            Err(PhaseHeatExchangeError::Invalid(_)),
        ));
        let zero = curve
            .linear_heat_exchange_response(initial, 400.0, 0.0)
            .unwrap();
        assert_eq!(zero.jacobian.inlet_enthalpy_gain, 1.0);
        assert_eq!(zero.jacobian.reservoir_enthalpy_derivative_j_kg_k, 0.0);
        // A nonstationary crossing still has a unique enthalpy response.
        assert!(
            curve
                .linear_heat_exchange_response(initial, 450.0, 20.0)
                .is_ok()
        );
    }
    let initial = curve.state_at_specific_enthalpy(22_000.0).unwrap();
    let mut checkpoints = 0;
    assert!(matches!(
        curve.linear_heat_exchange_response_with_checkpoint(initial, 550.0, 1200.0, || {
            checkpoints += 1;
            checkpoints == 3
        }),
        Err(PhaseHeatExchangeError::Cancelled),
    ));
    assert_eq!(checkpoints, 3);
}

#[test]
fn g3_response_composes_and_preserves_the_original_source_under_subdivision() {
    for datum in [0.0, 2.0_f64.powi(30)] {
        let curve = source_curve(datum, false);
        let subdivided = source_curve(datum, true);
        let initial_h = datum + 22_000.0;
        let initial = curve.state_at_specific_enthalpy(initial_h).unwrap();
        let direct = curve
            .linear_heat_exchange_response(initial, 550.0, 1200.0)
            .unwrap();
        let split = subdivided
            .linear_heat_exchange_response(
                subdivided.state_at_specific_enthalpy(initial_h).unwrap(),
                550.0,
                1200.0,
            )
            .unwrap();
        near(
            direct.jacobian.inlet_enthalpy_gain,
            split.jacobian.inlet_enthalpy_gain,
            3e-14,
        );
        near(
            direct.jacobian.reservoir_enthalpy_derivative_j_kg_k,
            split.jacobian.reservoir_enthalpy_derivative_j_kg_k,
            3e-9,
        );
        near(
            direct.jacobian.mean_reservoir_derivative,
            split.jacobian.mean_reservoir_derivative,
            3e-14,
        );
        near(
            direct.jacobian.mean_log_exposure_derivative_k,
            split.jacobian.mean_log_exposure_derivative_k,
            3e-10,
        );
        let first = curve
            .linear_heat_exchange_response(initial, 550.0, 600.0)
            .unwrap();
        let second = curve
            .linear_heat_exchange_response(first.exchange.state, 550.0, 600.0)
            .unwrap();
        let (a, b, d) = (first.jacobian, second.jacobian, direct.jacobian);
        near(
            d.inlet_enthalpy_gain,
            b.inlet_enthalpy_gain * a.inlet_enthalpy_gain,
            3e-12,
        );
        near(
            d.reservoir_enthalpy_derivative_j_kg_k,
            b.inlet_enthalpy_gain * a.reservoir_enthalpy_derivative_j_kg_k
                + b.reservoir_enthalpy_derivative_j_kg_k,
            3e-8,
        );
        near(
            d.mean_reservoir_derivative,
            0.5 * (a.mean_reservoir_derivative + b.mean_reservoir_derivative)
                + b.inlet_enthalpy_loss * a.reservoir_enthalpy_derivative_j_kg_k / 1200.0,
            3e-12,
        );
    }
}

#[test]
fn g1_phase_mean_exposure_response_retains_sub_ulp_temperature_travel() {
    // Put h=0 at 350 K so a representable tiny heat input remains storable even
    // though both absolute material temperatures round to the same 350 K.
    let curve = source_curve(-62_500.0, false);
    let initial = curve.state_at_specific_enthalpy(0.0).unwrap();
    assert_eq!(initial.temperature_k(), 350.0);
    let x = 1e-200;
    for bath in [250.0, 550.0] {
        let response = curve.linear_heat_exchange_response(initial, bath, x).unwrap();
        assert_eq!(response.exchange.state.temperature_k(), 350.0);
        assert_eq!(response.exchange.mean_temperature_k, 350.0);
        assert_ne!(response.exchange.state.specific_enthalpy_j_kg(), 0.0);
        // Independent local ODE expansion: Cp(350)=1500, so
        // T_end - T_mean = (T_bath-350)*X/(2*1500) + O(X^2).
        let expected = (bath - 350.0) * (x / 3000.0);
        near(response.jacobian.mean_log_exposure_derivative_k / expected, 1.0, 3e-14);
    }
}

