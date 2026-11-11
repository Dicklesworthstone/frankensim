//! Original sensible/latent source oracles for signed applied-power exchange.
use super::*;

// Independent integral of T*Cp(T) for an affine sensible interval.
fn temperature_moment(from: f64, to: f64, intercept: f64, slope: f64) -> f64 {
    let width = to - from;
    let cp = intercept + slope * from;
    width * (from * cp + 0.5 * (cp + slope * from) * width + slope * width * width / 3.0)
}

#[test]
fn g1_phase_driven_crosses_latent_heat_with_signed_source_and_negative_virtual_equilibrium() {
    let curve = source_curve(0.0, false);
    for (h0, target_h, target_t, virtual_t) in [
        (22_000.0, 326_000.0, 460.0, 550.0),
        (470_000.0, 34_500.0, 330.0, -200.0),
    ] {
        let initial = curve.state_at_specific_enthalpy(h0).unwrap();
        let x = if target_h > h0 {
            exposure(320.0, 400.0, virtual_t, -2000.0, 10.0)
                + 20_000.0 / (virtual_t - 400.0)
                + exposure(400.0, target_t, virtual_t, -6000.0, 20.0)
        } else {
            exposure(500.0, 400.0, virtual_t, -6000.0, 20.0)
                + 20_000.0 / (400.0 - virtual_t)
                + exposure(400.0, target_t, virtual_t, -2000.0, 10.0)
        };
        let drive = x * (virtual_t - initial.temperature_k());
        let response = curve.driven_heat_exchange_response(initial, drive, x).unwrap();
        let (value, jacobian) = (response.exchange, response.jacobian);
        near(value.state.specific_enthalpy_j_kg(), target_h, 2e-7);
        near(value.state.temperature_k(), target_t, 2e-10);
        near(value.specific_heat_j_kg, target_h - h0, 2e-7);
        near(value.mean_temperature_change_k, (drive - (target_h - h0)) / x, 2e-10);
        near(value.temperature_change_k, target_t - initial.temperature_k(), 2e-10);
        near(jacobian.inlet_enthalpy_gain,
            (virtual_t - target_t) / (virtual_t - initial.temperature_k()), 2e-12);
        near(jacobian.remaining_drive_j_kg, x * (virtual_t - target_t), 2e-7);
        near(jacobian.mean_log_duration_derivative_k,
            value.temperature_change_k - value.mean_temperature_change_k, 2e-10);
        assert_eq!(value.state.phase_curve_identity(), curve.identity());
        assert_eq!(value.state.material_card_identity(), curve.material_card_identity());
    }
}

#[test]
fn g1_zero_conductance_retains_exact_source_heat_and_latent_weighted_time_mean() {
    let curve = source_curve(0.0, false);
    for (h0, target_h, target_t, moment) in [
        (22_000.0, 326_000.0, 460.0,
            temperature_moment(320.0, 400.0, -2000.0, 10.0) + 400.0 * 20_000.0
                + temperature_moment(400.0, 460.0, -6000.0, 20.0)),
        (470_000.0, 34_500.0, 330.0,
            temperature_moment(500.0, 400.0, -6000.0, 20.0) - 400.0 * 20_000.0
                + temperature_moment(400.0, 330.0, -2000.0, 10.0)),
    ] {
        let initial = curve.state_at_specific_enthalpy(h0).unwrap();
        let q = target_h - h0;
        let response = curve.driven_heat_exchange_response(initial, q, 0.0).unwrap();
        assert_eq!(response.exchange.specific_heat_j_kg, q);
        assert_eq!(response.jacobian.source_heat_gain, 1.0);
        assert_eq!(response.jacobian.inlet_enthalpy_gain, 1.0);
        near(response.exchange.state.specific_enthalpy_j_kg(), target_h, 2e-7);
        near(response.exchange.mean_temperature_k, moment / q, 2e-10);
        near(response.jacobian.mean_source_derivative_kg_k_j, (target_t - moment / q) / q, 2e-14);
        near(response.exchange.heat_residual_j_kg, 0.0, 2e-7);
    }
    let initial = curve.state_at_specific_enthalpy(155_000.0).unwrap();
    let latent = curve.driven_heat_exchange_response(initial, 5000.0, 1e6).unwrap();
    assert_eq!(latent.exchange.state.specific_enthalpy_j_kg(), 160_000.0);
    assert_eq!(latent.exchange.state.liquid_mass_fraction(), 0.5);
    assert_eq!(latent.exchange.mean_temperature_k, 400.0);
    assert_eq!(latent.exchange.temperature_change_k, 0.0);
    assert_eq!(latent.jacobian.source_heat_gain, 1.0);
    assert_eq!(latent.jacobian.mean_source_derivative_kg_k_j, 0.0);
    assert_eq!(latent.jacobian.remaining_drive_j_kg, 5000.0);
}

#[test]
fn g3_driven_phase_source_mean_inlet_and_duration_response_cross_latent_residence() {
    let curve = source_curve(0.0, false);
    let (h, drive, x) = (22_000.0, 500_000.0, 1200.0);
    let initial = curve.state_at_specific_enthalpy(h).unwrap();
    let response = curve.driven_heat_exchange_response(initial, drive, x).unwrap();
    assert!(response.exchange.state.liquid_mass_fraction() == 1.0);
    let [lo, hi] = [-0.1, 0.1].map(|d| curve.driven_heat_exchange(initial, drive + d, x).unwrap());
    near((hi.specific_heat_j_kg - lo.specific_heat_j_kg) / 0.2, response.jacobian.source_heat_gain, 2e-7);
    near((hi.mean_temperature_k - lo.mean_temperature_k) / 0.2,
        response.jacobian.mean_source_derivative_kg_k_j, 2e-9);
    let [lo, hi] = [-0.1, 0.1].map(|d| {
        let inlet = curve.state_at_specific_enthalpy(h + d).unwrap();
        curve.driven_heat_exchange(
            inlet, drive - x * (inlet.temperature_k() - initial.temperature_k()), x,
        ).unwrap().state.specific_enthalpy_j_kg()
    });
    near((hi - lo) / 0.2, response.jacobian.inlet_enthalpy_gain, 2e-7);
    let [lo, hi] = [-1e-4, 1e-4].map(|d| {
        curve.driven_heat_exchange(initial, drive * (1.0 + d), x * (1.0 + d)).unwrap()
    });
    near((hi.mean_temperature_k - lo.mean_temperature_k) / 2e-4,
        response.jacobian.mean_log_duration_derivative_k, 2e-6);
    let first = curve.driven_heat_exchange_response(initial, drive / 2.0, x / 2.0).unwrap();
    let second = curve.driven_heat_exchange(
        first.exchange.state, first.jacobian.remaining_drive_j_kg, x / 2.0,
    ).unwrap();
    near(first.exchange.specific_heat_j_kg + second.specific_heat_j_kg,
        response.exchange.specific_heat_j_kg, 2e-7);
    near(0.5 * (first.exchange.mean_temperature_change_k + first.exchange.temperature_change_k
        + second.mean_temperature_change_k), response.exchange.mean_temperature_change_k, 2e-10);
}

#[test]
fn g1_phase_driven_weak_heat_retains_both_temperature_moments() {
    let curve = source_curve(-62_500.0, false);
    let initial = curve.state_at_specific_enthalpy(0.0).unwrap();
    let x: f64 = 5.0;
    let a = x / 1500.0;
    let phi = -(-a).exp_m1() / a;
    let psi = (1.0 - phi) / a;
    for q in [-1e-18, 1e-18] {
        let response = curve.driven_heat_exchange_response(initial, q, x).unwrap();
        assert_eq!(response.exchange.state.temperature_k(), 350.0);
        assert_eq!(response.exchange.mean_temperature_k, 350.0);
        assert_ne!(response.exchange.state.specific_enthalpy_j_kg(), 0.0);
        near(response.exchange.specific_heat_j_kg / (q * phi), 1.0, 3e-14);
        near(response.exchange.temperature_change_k / (q * phi / 1500.0), 1.0, 3e-14);
        near(response.exchange.mean_temperature_change_k / (q * psi / 1500.0), 1.0, 1e-12);
        assert_eq!(response.exchange.heat_residual_j_kg, 0.0);
    }
}

#[test]
fn g0_phase_driven_corners_domain_identity_and_cancellation_preserve_original_contract() {
    let curve = source_curve(0.0, false);
    for h in [150_000.0, 170_000.0] {
        let initial = curve.state_at_specific_enthalpy(h).unwrap();
        for x in [0.0, 1000.0] {
            assert_eq!(curve.driven_heat_exchange(initial, 0.0, x).unwrap().state, initial);
            assert!(curve.driven_heat_exchange_response(initial, 0.0, x).is_err());
        }
    }
    let plateau = curve.state_at_specific_enthalpy(155_000.0).unwrap();
    let eq = curve.driven_heat_exchange_response(plateau, 0.0, 0.0).unwrap();
    assert_eq!(eq.jacobian.source_heat_gain, 1.0);
    assert_eq!(eq.jacobian.mean_source_derivative_kg_k_j, 0.0);
    let initial = curve.state_at_specific_enthalpy(22_000.0).unwrap();
    for q in [-1e9, 1e9, f64::NAN, f64::INFINITY] {
        assert!(curve.driven_heat_exchange(initial, q, 0.0).is_err());
    }
    let foreign = source_curve(1.0, false).state_at_specific_enthalpy(22_001.0).unwrap();
    assert!(curve.driven_heat_exchange(foreign, 0.0, 0.0).is_err());
    let shifted = source_curve(2.0_f64.powi(40), false);
    let large = shifted.state_at_specific_enthalpy(2.0_f64.powi(40) + 62_500.0).unwrap();
    assert!(matches!(
        shifted.driven_heat_exchange(large, 1e-20, 0.0),
        Err(PhaseHeatExchangeError::Phase(PhaseStateError::LostSpecificEnergyIncrement))
    ));
    let mut calls = 0;
    let expected = curve.driven_heat_exchange_with_checkpoint(initial, 304_000.0, 0.0, || {
        calls += 1; false
    }).unwrap();
    assert!(calls >= 5);
    let mut replay_calls = 0;
    assert!(matches!(
        curve.driven_heat_exchange_with_checkpoint(initial, 304_000.0, 0.0, || {
            replay_calls += 1; replay_calls == calls
        }),
        Err(PhaseHeatExchangeError::Cancelled)
    ));
    assert_eq!(curve.driven_heat_exchange(initial, 304_000.0, 0.0).unwrap(), expected);
}

#[test]
fn g1_reachable_latent_heat_survives_a_rounded_virtual_equilibrium_at_the_boundary() {
    let (initial_t, melting_t) = (0.5, 1.5 + f64::EPSILON);
    let width = melting_t - initial_t;
    let sensible_heat = 1000.0 * width;
    let mut curve = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x9d; 32]),
        vec![
            EnthalpyPhaseKnot { specific_enthalpy_j_kg: 0.0, temperature_k: initial_t,
                liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1000.0 },
            EnthalpyPhaseKnot { specific_enthalpy_j_kg: sensible_heat, temperature_k: melting_t,
                liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1000.0 },
            EnthalpyPhaseKnot { specific_enthalpy_j_kg: sensible_heat + 1.0, temperature_k: melting_t,
                liquid_mass_fraction: 1.0, bulk_density_kg_m3: 900.0 },
        ],
    ).unwrap();
    curve.retain_sensible_heat_capacities(vec![Some([1000.0; 2]), None]).unwrap();
    let drive = 2.0_f64.powi(100);
    let exposure = drive * (1.0 - f64::EPSILON);
    assert_eq!(drive / exposure, width);
    assert_eq!((-exposure).mul_add(width, drive), 0.0625);
    let initial = curve.state_at_specific_enthalpy(0.0).unwrap();
    let response = curve.driven_heat_exchange_response(initial, drive, exposure).unwrap();
    near(response.exchange.specific_heat_j_kg, sensible_heat + 0.0625, 3e-12);
    near(response.exchange.state.liquid_mass_fraction(), 0.0625, 3e-12);
    near(response.jacobian.remaining_drive_j_kg, 0.0625, 2e-15);
    assert_eq!(response.exchange.state.temperature_k(), melting_t);
}

#[test]
fn g1_source_only_phase_storage_is_independent_of_sensible_subdivisions() {
    let datum = 1e20;
    let curve = |split: bool| {
        let offsets: &[f64] = if split { &[0.0, 16_384.0, 49_152.0, 98_304.0] }
            else { &[0.0, 98_304.0] };
        EquilibriumEnthalpyPhaseCurve::try_single_phase(ContentHash([0xe3; 32]),
            SolidLiquidPhase::Solid,
            offsets.iter().map(|&offset| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: datum + offset,
                temperature_k: 300.0 + offset / 163.84,
                liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1000.0,
            }).collect()).unwrap()
    };
    let (direct, split) = (curve(false), curve(true));
    for (offset, heat) in [(0.0, 73_728.0), (98_304.0, -73_728.0)] {
        let initial = split.state_at_specific_enthalpy(datum + offset).unwrap();
        let expected = split.advance_specific_energy(initial, heat).unwrap();
        let result = split.driven_heat_exchange_response(initial, heat, 0.0).unwrap();
        let unsplit = direct.driven_heat_exchange(
            direct.state_at_specific_enthalpy(datum + offset).unwrap(), heat, 0.0,
        ).unwrap();
        // Native storage rounds once from the supplied input. Reconstructing
        // the last source interval gives a different tie direction here.
        assert_eq!(result.exchange.state, expected);
        assert_eq!(result.exchange.state.specific_enthalpy_j_kg(), unsplit.state.specific_enthalpy_j_kg());
        assert_eq!(result.exchange.specific_heat_j_kg, heat);
        assert_eq!(result.exchange.heat_residual_j_kg,
            heat - (expected.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg()));
        near(result.exchange.temperature_change_k, heat / 163.84, 2e-10);
        near(result.exchange.mean_temperature_k, initial.temperature_k() + heat / (2.0 * 163.84), 2e-10);
        assert_eq!(result.jacobian.source_heat_gain, 1.0);
        near(result.jacobian.mean_source_derivative_kg_k_j, 1.0 / (2.0 * 163.84), 2e-14);
    }
}

#[test]
fn g1_source_only_phase_exchange_reaches_finite_source_endpoints() {
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(ContentHash([0xe4; 32]),
        SolidLiquidPhase::Solid,
        [0.0, 1.0, 4.0, 7.0].into_iter().map(|h| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h, temperature_k: 300.0 + h,
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1000.0,
        }).collect()).unwrap();
    for (h0, heat) in [(0.0, 7.0), (7.0, -7.0)] {
        let initial = curve.state_at_specific_enthalpy(h0).unwrap();
        let result = curve.driven_heat_exchange_response(initial, heat, 0.0).unwrap();
        assert_eq!(result.exchange.state, curve.advance_specific_energy(initial, heat).unwrap());
        assert_eq!(result.exchange.specific_heat_j_kg, heat);
        assert_eq!(result.exchange.heat_residual_j_kg, 0.0);
        near(result.exchange.mean_temperature_k, 303.5, 2e-13);
        near(result.jacobian.mean_log_duration_derivative_k, heat / 2.0, 2e-13);
        assert_eq!(result.jacobian.source_heat_gain, 1.0);
    }
}
