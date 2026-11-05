//! G1/G3: carrier inertia, phase, heat and reversible kinetic storage.
use super::*;
use fs_math::c64::C64;

fn complex_relative(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn g1_original_carriers_have_drude_phase_and_separate_heat_from_kinetic_storage() {
    let gas = material(8000.0, 100_000.0);
    let prepared = gas.prepare_harmonic(|| false).unwrap();
    assert!(std::ptr::eq(prepared.source(), &gas));
    let field = [
        C64::from_re(0.01),
        C64::new(0.0, 0.02),
        C64::new(0.005, -0.003),
    ];
    let mean_e2 = 0.5 * field.iter().map(|e| e.norm_sq()).sum::<f64>();
    for omega in [1e-6, 1e9, core::f64::consts::TAU * 1e9, 1e10, 1e18] {
        let response = prepared.response(omega, field, || false).unwrap();
        assert!(std::ptr::eq(response.source, &gas));
        let mut conductivity = C64::ZERO;
        let mut kinetic = 0.0;
        let mut heat = [0.0; 2];
        for (i, carrier) in gas.carriers().iter().enumerate() {
            let charge = carrier.charge_number() as f64 * E;
            let coefficient = carrier.number_density_per_m3() * charge * charge / carrier.mass_kg();
            let denominator = C64::new(carrier.relaxation().rate_per_s(), omega);
            conductivity = conductivity + C64::from_re(coefficient) / denominator;
            let mobility = C64::from_re(charge / carrier.mass_kg()) / denominator;
            for (actual, electric) in response.carrier_velocity_m_s[i].iter().zip(field) {
                complex_relative(*actual, mobility * electric, 2e-13);
            }
            let mean_speed2 = mobility.norm_sq() * mean_e2;
            let inertia = carrier.number_density_per_m3() * carrier.mass_kg();
            kinetic += 0.5 * inertia * mean_speed2;
            heat[usize::from(carrier.is_electron())] +=
                inertia * carrier.relaxation().rate_per_s() * mean_speed2;
        }
        complex_relative(response.conductivity_s_m, conductivity, 2e-13);
        for (current, electric) in response.current_density_a_m2.iter().zip(field) {
            complex_relative(*current, conductivity * electric, 2e-13);
        }
        relative(response.mean_drift_kinetic_energy_j_m3, kinetic, 2e-13);
        relative(
            response.mean_heating_w_m3.electron_power_density_w_m3(),
            heat[1],
            2e-13,
        );
        relative(
            response.mean_heating_w_m3.heavy_power_density_w_m3(),
            heat[0],
            2e-13,
        );
        relative(
            response.mean_electrical_power_w_m3,
            conductivity.re * mean_e2,
            2e-13,
        );
        relative(response.reactive_power_var_m3, 2.0 * omega * kinetic, 2e-13);
        assert!(response.force_balance_relative_residual < 2e-13);
        assert!(response.power_balance_relative_residual < 2e-13);
    }
    let low = prepared.response(1e-6, field, || false).unwrap();
    relative(low.conductivity_s_m.re, gas.conductivity_s_m(), 2e-13);
    let high = prepared.response(1e18, field, || false).unwrap();
    let weight: f64 = gas
        .carriers()
        .iter()
        .map(|c| {
            c.number_density_per_m3() * (c.charge_number() as f64 * E).powi(2) / c.mass_kg()
        })
        .sum();
    relative(-high.conductivity_s_m.im * 1e18, weight, 2e-13);
    assert!(high.mean_electrical_power_w_m3 < 1e-14 * low.mean_electrical_power_w_m3);
}

#[test]
fn g1_coulomb_harmonics_match_the_independent_two_fluid_complex_inverse() {
    let state = frozen_state(0.001, 6000.0, 12000.0, 1.0);
    let th = state.heavy_temperature().value();
    let te = state.electron_temperature().value();
    let p = state.pressure().value();
    let gas = IonizedGasConductivity::from_state(
        state.clone(),
        vec![
            source_two_temperature("e-", th, te, p, 1e10),
            source_two_temperature("H+", th, te, p, 1e6),
        ],
        limits(),
        || false,
    )
    .unwrap()
    .with_coulomb_drag(coulomb_model(&state), || false)
    .unwrap();
    let prepared = gas.prepare_harmonic(|| false).unwrap();
    let electron = gas
        .carriers()
        .iter()
        .position(|c| c.is_electron())
        .unwrap();
    let ion = gas
        .carriers()
        .iter()
        .position(|c| !c.is_electron())
        .unwrap();
    let e = &gas.carriers()[electron];
    let i = &gas.carriers()[ion];
    let me = e.number_density_per_m3() * e.mass_kg();
    let mi = i.number_density_per_m3() * i.mass_kg();
    let ae = me * e.relaxation().rate_per_s();
    let ai = mi * i.relaxation().rate_per_s();
    let pair = &gas.coulomb_pairs()[0];
    let k = pair.friction_coefficient_kg_m3_s();
    let b = e.number_density_per_m3() * E;
    let field = [C64::from_re(0.01), C64::ZERO, C64::ZERO];
    for omega in [1e-6, 1e6, 1e10, 1e15] {
        let response = prepared.response(omega, field, || false).unwrap();
        let de = C64::new(ae, omega * me);
        let di = C64::new(ai, omega * mi);
        let determinant = de * di + (de + di).scale(k);
        let ue = di.scale(-b) / determinant;
        let ui = de.scale(b) / determinant;
        complex_relative(
            response.carrier_velocity_m_s[electron][0],
            ue.scale(0.01),
            3e-11,
        );
        complex_relative(
            response.carrier_velocity_m_s[ion][0],
            ui.scale(0.01),
            3e-11,
        );
        complex_relative(response.conductivity_s_m, (ui - ue).scale(b), 3e-11);
        let relative_speed2 = (ue - ui).norm_sq() * 0.5e-4;
        let total_mass = e.mass_kg() + i.mass_kg();
        let electron_heat =
            0.5e-4 * ae * ue.norm_sq() + k * relative_speed2 * i.mass_kg() / total_mass;
        let heavy_heat =
            0.5e-4 * ai * ui.norm_sq() + k * relative_speed2 * e.mass_kg() / total_mass;
        let correction = pair.thermal_conductance_w_m3_k() * (te - th) * relative_speed2
            / (2.0 * pair.relative_thermal_variance_m2_s2());
        relative(
            response.mean_heating_w_m3.electron_friction_power_density_w_m3(),
            electron_heat,
            5e-11,
        );
        relative(
            response.mean_heating_w_m3.heavy_friction_power_density_w_m3(),
            heavy_heat,
            5e-11,
        );
        relative(
            response
                .mean_heating_w_m3
                .thermal_exchange_correction_to_electrons_w_m3(),
            correction,
            5e-11,
        );
        relative(
            response.mean_drift_kinetic_energy_j_m3,
            0.25e-4 * (me * ue.norm_sq() + mi * ui.norm_sq()),
            5e-11,
        );
        assert!(response.power_balance_relative_residual < 1e-10);
    }
    let dc = gas.joule_heating_w_m3([0.01, 0.0]).unwrap();
    let low = prepared.response(1e-6, field, || false).unwrap();
    relative(
        low.mean_heating_w_m3.total_power_density_w_m3(),
        0.5 * dc.total_power_density_w_m3(),
        3e-11,
    );
    relative(low.conductivity_s_m.re, gas.conductivity_s_m(), 3e-11);
}

#[test]
fn g0_actual_polarized_peak_drift_and_frequency_semantics_are_enforced() {
    let gas = material(8000.0, 100_000.0);
    assert!(matches!(
        gas.prepare_harmonic(|| true),
        Err(IonizedGasError::Cancelled)
    ));
    let prepared = gas.prepare_harmonic(|| false).unwrap();
    let limit = gas.maximum_electric_field_v_m();
    // A circular vector has constant magnitude 0.8*limit, not sqrt(2)*0.8*limit.
    let circular = [
        C64::from_re(0.8 * limit),
        C64::new(0.0, 0.8 * limit),
        C64::ZERO,
    ];
    prepared.response(1e-6, circular, || false).unwrap();
    let excessive = [C64::from_re(1.1 * limit), C64::ZERO, C64::ZERO];
    assert!(prepared.response(1e-6, excessive, || false).is_err());
    // At high frequency inertia reduces actual drift; the old DC field bound
    // is not incorrectly reused as an AC amplitude bound.
    prepared.response(1e12, excessive, || false).unwrap();
    for frequency in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(prepared.response(frequency, circular, || false).is_err());
    }
    assert!(
        prepared
            .response(1e9, [C64::new(f64::NAN, 0.0); 3], || false)
            .is_err()
    );
    assert!(matches!(
        prepared.response(1e9, circular, || true),
        Err(IonizedGasError::Cancelled)
    ));
    let zero = prepared.response(1e9, [C64::ZERO; 3], || false).unwrap();
    assert_eq!(zero.current_density_a_m2, [C64::ZERO; 3]);
    assert_eq!(zero.mean_heating_w_m3.total_power_density_w_m3(), 0.0);
    assert_eq!(zero.mean_drift_kinetic_energy_j_m3, 0.0);
}
