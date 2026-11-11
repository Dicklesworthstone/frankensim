//! G1/G3/G4: signed gyromotion, original mutual drag, and physical work/heat.
use super::*;
use fs_material::ionized_gas::IonizedGasMagneticLimits;
use fs_math::c64::C64;

fn magnetic_limits() -> IonizedGasMagneticLimits {
    IonizedGasMagneticLimits {
        maximum_cyclotron_to_thermal_energy_ratio: 0.01,
    }
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    let scale = actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{actual:?} != {expected:?}"
    );
}
fn coupled() -> IonizedGasConductivity {
    let state = frozen_state(0.001, 6000.0, 12000.0, 1.0);
    let th = state.heavy_temperature().value();
    let te = state.electron_temperature().value();
    let p = state.pressure().value();
    IonizedGasConductivity::from_state(
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
    .unwrap()
}

#[test]
fn g1_zero_field_retains_original_harmonic_sources_velocities_and_thermal_pools() {
    let field = [
        C64::new(0.01, -0.005),
        C64::new(0.0, 0.02),
        C64::new(0.003, 0.004),
    ];
    for gas in [material(8000.0, 100_000.0), coupled()] {
        let original = gas.prepare_harmonic(|| false).unwrap();
        for omega in [1e-6, 1e6, 1e10, 1e15] {
            let reference = original.response(omega, field, || false).unwrap();
            let prepared = original
                .prepare_magnetized(omega, [0.0; 3], magnetic_limits(), || false)
                .unwrap();
            let response = prepared.response(field, || false).unwrap();
            assert!(std::ptr::eq(prepared.source(), &gas));
            assert!(std::ptr::eq(response.source, &gas));
            assert_eq!(prepared.maximum_cyclotron_to_thermal_energy_ratio(), 0.0);
            for row in 0..3 {
                close(
                    response.current_density_a_m2[row],
                    reference.current_density_a_m2[row],
                    3e-12,
                );
                for column in 0..3 {
                    let expected = if row == column {
                        reference.conductivity_s_m
                    } else {
                        C64::ZERO
                    };
                    close(
                        prepared.conductivity_tensor_s_m()[row][column],
                        expected,
                        2e-14,
                    );
                }
                for carrier in 0..gas.carriers().len() {
                    close(
                        response.carrier_velocity_m_s[carrier][row],
                        reference.carrier_velocity_m_s[carrier][row],
                        3e-12,
                    );
                }
            }
            relative(
                response.mean_heating_w_m3.electron_power_density_w_m3(),
                reference.mean_heating_w_m3.electron_power_density_w_m3(),
                3e-12,
            );
            relative(
                response.mean_heating_w_m3.heavy_power_density_w_m3(),
                reference.mean_heating_w_m3.heavy_power_density_w_m3(),
                3e-12,
            );
            relative(
                response.mean_drift_kinetic_energy_j_m3,
                reference.mean_drift_kinetic_energy_j_m3,
                3e-12,
            );
            relative(
                response.reactive_power_var_m3,
                reference.reactive_power_var_m3,
                3e-12,
            );
            assert_eq!(response.gyroscopic_reactive_power_var_m3, 0.0);
            assert_eq!(response.mean_magnetic_power_w_m3, 0.0);
            assert!(response.force_balance_relative_residual < 1e-10);
            assert!(response.power_balance_relative_residual < 1e-10);
        }
    }
}

#[test]
fn g1_signed_drude_hall_response_matches_independent_cartesian_inverse() {
    let gas = material(8000.0, 100_000.0);
    let original = gas.prepare_harmonic(|| false).unwrap();
    let field = [
        C64::new(0.01, 0.005),
        C64::new(-0.004, 0.003),
        C64::new(0.002, -0.001),
    ];
    for b in [-0.08_f64, 0.08] {
        for omega in [1e-6, 2e9, E * b.abs() / ME] {
            let prepared = original
                .prepare_magnetized(omega, [0.0, 0.0, b], magnetic_limits(), || false)
                .unwrap();
            let response = prepared.response(field, || false).unwrap();
            let mut current = [C64::ZERO; 3];
            for (i, carrier) in gas.carriers().iter().enumerate() {
                let charge_mass = carrier.charge_number() as f64 * E / carrier.mass_kg();
                let cyclotron = charge_mass * b;
                let d = C64::new(carrier.relaxation().rate_per_s(), omega);
                let determinant = d * d + C64::from_re(cyclotron * cyclotron);
                let expected = [
                    (d * field[0] + field[1].scale(cyclotron)).scale(charge_mass) / determinant,
                    (d * field[1] - field[0].scale(cyclotron)).scale(charge_mass) / determinant,
                    field[2].scale(charge_mass) / d,
                ];
                for dimension in 0..3 {
                    close(
                        response.carrier_velocity_m_s[i][dimension],
                        expected[dimension],
                        3e-13,
                    );
                    current[dimension] = current[dimension]
                        + expected[dimension].scale(
                            carrier.number_density_per_m3() * carrier.charge_number() as f64 * E,
                        );
                }
            }
            for dimension in 0..3 {
                close(
                    response.current_density_a_m2[dimension],
                    current[dimension],
                    3e-13,
                );
            }
            relative(
                response.mean_electrical_power_w_m3,
                response.mean_heating_w_m3.total_power_density_w_m3(),
                3e-13,
            );
            assert!(
                response.mean_magnetic_power_w_m3.abs()
                    < 1e-13 * response.mean_electrical_power_w_m3
            );
        }
    }
    let low = original
        .prepare_magnetized(1e-6, [0.0, 0.0, 0.08], magnetic_limits(), || false)
        .unwrap();
    // Electron-dominated Hall current points +y for E along +x and B along +z.
    assert!(low.conductivity_tensor_s_m()[1][0].re > 0.0);
    assert_eq!(
        low.conductivity_tensor_s_m()[1][0],
        -low.conductivity_tensor_s_m()[0][1]
    );
}

#[test]
fn g1_both_charge_signs_resonate_in_their_own_circular_polarization() {
    let b = 0.1;
    let owner = chemistry();
    let state = equilibrium(&owner, 8000.0, 100_000.0);
    let gas = IonizedGasConductivity::new(
        &owner,
        state,
        vec![
            source("e-", 8000.0, 100_000.0, E * b / ME / 20.0),
            source("H+", 8000.0, 100_000.0, E * b / (MH - ME) / 20.0),
        ],
        limits(),
        || false,
    )
    .unwrap();
    let original = gas.prepare_harmonic(|| false).unwrap();
    for (i, carrier) in gas.carriers().iter().enumerate() {
        let sign = (carrier.charge_number() as f64).signum();
        let omega = (carrier.charge_number() as f64 * E / carrier.mass_kg() * b).abs();
        let prepared = original
            .prepare_magnetized(omega, [0.0, 0.0, b], magnetic_limits(), || false)
            .unwrap();
        let field = [C64::from_re(0.01), C64::new(0.0, sign * 0.01), C64::ZERO];
        let on = prepared.response(field, || false).unwrap();
        let off = prepared.response(field.map(C64::conj), || false).unwrap();
        let resonant = carrier.charge_number() as f64 * E
            / carrier.mass_kg()
            / carrier.relaxation().rate_per_s();
        for d in 0..3 {
            close(
                on.carrier_velocity_m_s[i][d],
                field[d].scale(resonant),
                3e-13,
            );
        }
        assert!(on.carrier_velocity_m_s[i][0].abs() > 39.0 * off.carrier_velocity_m_s[i][0].abs());
        // The real magnetic work remains zero while its reactive term is material.
        assert!(
            on.gyroscopic_reactive_power_var_m3.abs()
                > 0.1 * 2.0 * omega * on.mean_drift_kinetic_energy_j_m3
        );
        let gross = (2.0 * omega * on.mean_drift_kinetic_energy_j_m3)
            .max(on.gyroscopic_reactive_power_var_m3.abs());
        assert!(
            (on.reactive_power_var_m3
                - 2.0 * omega * on.mean_drift_kinetic_energy_j_m3
                - on.gyroscopic_reactive_power_var_m3)
                .abs()
                < 2e-13 * gross
        );
        assert!(on.mean_magnetic_power_w_m3.abs() < 2e-13 * on.mean_electrical_power_w_m3);
    }
}

#[test]
fn g1_mutual_coulomb_drag_keeps_the_complete_circular_two_fluid_inverse_and_heat() {
    let gas = coupled();
    let original = gas.prepare_harmonic(|| false).unwrap();
    let electron = gas.carriers().iter().position(|c| c.is_electron()).unwrap();
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
    let charge = e.number_density_per_m3() * E;
    let b = 0.3;
    for omega in [1e6, 1e10, 5e10] {
        let prepared = original
            .prepare_magnetized(omega, [0.0, 0.0, b], magnetic_limits(), || false)
            .unwrap();
        for hand in [-1.0, 1.0] {
            let field = [C64::from_re(0.01), C64::new(0.0, hand * 0.01), C64::ZERO];
            let response = prepared.response(field, || false).unwrap();
            let de = C64::new(ae, me * (omega + hand * E * b / e.mass_kg()));
            let di = C64::new(ai, mi * (omega - hand * E * b / i.mass_kg()));
            let determinant = de * di + (de + di).scale(k);
            let ue = di.scale(-charge) / determinant;
            let ui = de.scale(charge) / determinant;
            for dimension in 0..3 {
                close(
                    response.carrier_velocity_m_s[electron][dimension],
                    ue * field[dimension],
                    3e-11,
                );
                close(
                    response.carrier_velocity_m_s[ion][dimension],
                    ui * field[dimension],
                    3e-11,
                );
                close(
                    response.current_density_a_m2[dimension],
                    (ui - ue).scale(charge) * field[dimension],
                    3e-11,
                );
            }
            // Circular velocity has constant speed, so no additional factor 1/2.
            let relative2 = (ue - ui).norm_sq() * 1e-4;
            let total_mass = e.mass_kg() + i.mass_kg();
            let electron_heat = 1e-4 * ae * ue.norm_sq() + k * relative2 * i.mass_kg() / total_mass;
            let heavy_heat = 1e-4 * ai * ui.norm_sq() + k * relative2 * e.mass_kg() / total_mass;
            let correction = pair.thermal_conductance_w_m3_k()
                * (e.temperature_k() - i.temperature_k())
                * relative2
                / (2.0 * pair.relative_thermal_variance_m2_s2());
            relative(
                response
                    .mean_heating_w_m3
                    .electron_friction_power_density_w_m3(),
                electron_heat,
                5e-11,
            );
            relative(
                response
                    .mean_heating_w_m3
                    .heavy_friction_power_density_w_m3(),
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
                0.5e-4 * (me * ue.norm_sq() + mi * ui.norm_sq()),
                5e-11,
            );
            assert!(response.force_balance_relative_residual < 1e-10);
            assert!(response.power_balance_relative_residual < 1e-10);
        }
    }
}

const ROTATION: [[f64; 3]; 3] = [
    [1.0 / 9.0, -4.0 / 9.0, 8.0 / 9.0],
    [8.0 / 9.0, 4.0 / 9.0, 1.0 / 9.0],
    [-4.0 / 9.0, 7.0 / 9.0, 4.0 / 9.0],
];
fn rotate(v: [C64; 3]) -> [C64; 3] {
    std::array::from_fn(|i| (0..3).fold(C64::ZERO, |a, j| a + v[j].scale(ROTATION[i][j])))
}

#[test]
fn g3_rotation_covariance_and_onsager_casimir_field_reversal_preserve_work() {
    let gas = coupled();
    let original = gas.prepare_harmonic(|| false).unwrap();
    let b = [0.05, -0.08, 0.1];
    let rotated_b = rotate(b.map(C64::from_re)).map(|z| z.re);
    let e = [
        C64::new(0.01, -0.005),
        C64::new(0.002, 0.003),
        C64::new(-0.004, 0.008),
    ];
    let baseline = original
        .prepare_magnetized(1e10, b, magnetic_limits(), || false)
        .unwrap();
    let rotated = original
        .prepare_magnetized(1e10, rotated_b, magnetic_limits(), || false)
        .unwrap();
    let reversed = original
        .prepare_magnetized(1e10, b.map(|v| -v), magnetic_limits(), || false)
        .unwrap();
    let first = baseline.response(e, || false).unwrap();
    let second = rotated.response(rotate(e), || false).unwrap();
    for i in 0..3 {
        close(
            second.current_density_a_m2[i],
            rotate(first.current_density_a_m2)[i],
            2e-12,
        );
        for j in 0..3 {
            close(
                reversed.conductivity_tensor_s_m()[i][j],
                baseline.conductivity_tensor_s_m()[j][i],
                2e-12,
            );
        }
        for carrier in 0..gas.carriers().len() {
            close(
                second.carrier_velocity_m_s[carrier][i],
                rotate(first.carrier_velocity_m_s[carrier])[i],
                2e-12,
            );
        }
    }
    relative(
        second.mean_electrical_power_w_m3,
        first.mean_electrical_power_w_m3,
        2e-12,
    );
    relative(
        second.mean_heating_w_m3.total_power_density_w_m3(),
        first.mean_heating_w_m3.total_power_density_w_m3(),
        2e-12,
    );
    relative(
        second.mean_drift_kinetic_energy_j_m3,
        first.mean_drift_kinetic_energy_j_m3,
        2e-12,
    );
    relative(
        second.gyroscopic_reactive_power_var_m3,
        first.gyroscopic_reactive_power_var_m3,
        2e-12,
    );
}

#[test]
fn g0_g4_resonant_peak_drift_classical_domains_and_cancellation_are_enforced() {
    let gas = material(8000.0, 100_000.0);
    let original = gas.prepare_harmonic(|| false).unwrap();
    let omega = 1e12;
    let b = omega
        * gas
            .carriers()
            .iter()
            .find(|c| c.is_electron())
            .unwrap()
            .mass_kg()
        / E;
    let prepared = original
        .prepare_magnetized(omega, [0.0, 0.0, b], magnetic_limits(), || false)
        .unwrap();
    let limit = gas.maximum_electric_field_v_m();
    let circular = [
        C64::from_re(0.8 * limit),
        C64::new(0.0, -0.8 * limit),
        C64::ZERO,
    ];
    prepared.response(circular, || false).unwrap();
    let excessive = circular.map(|v| v.scale(2.0));
    // Unmagnetized inertia admits this field; cyclotron resonance removes that
    // protection, so the actual magnetized velocity must be checked again.
    original.response(omega, excessive, || false).unwrap();
    assert!(prepared.response(excessive, || false).is_err());
    assert!(
        prepared
            .response([C64::new(f64::NAN, 0.0); 3], || false)
            .is_err()
    );
    assert!(matches!(
        prepared.response(circular, || true),
        Err(IonizedGasError::Cancelled)
    ));
    for w in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            original
                .prepare_magnetized(w, [0.0, 0.0, b], magnetic_limits(), || false)
                .is_err()
        );
    }
    for bad_b in [[f64::NAN; 3], [f64::INFINITY; 3], [0.0, 0.0, 1e6]] {
        assert!(
            original
                .prepare_magnetized(omega, bad_b, magnetic_limits(), || false)
                .is_err()
        );
    }
    for ceiling in [0.0, -1.0, 1.0, f64::NAN] {
        assert!(
            original
                .prepare_magnetized(
                    omega,
                    [0.0, 0.0, b],
                    IonizedGasMagneticLimits {
                        maximum_cyclotron_to_thermal_energy_ratio: ceiling
                    },
                    || false
                )
                .is_err()
        );
    }
    let zero = prepared.response([C64::ZERO; 3], || false).unwrap();
    assert_eq!(zero.current_density_a_m2, [C64::ZERO; 3]);
    assert_eq!(zero.mean_electrical_power_w_m3, 0.0);
    assert_eq!(zero.mean_heating_w_m3.total_power_density_w_m3(), 0.0);
    assert_eq!(zero.mean_drift_kinetic_energy_j_m3, 0.0);
    let coupled = coupled();
    let owner = coupled.prepare_harmonic(|| false).unwrap();
    let mut calls = 0;
    assert!(matches!(
        owner.prepare_magnetized(1e10, [0.0, 0.0, 0.3], magnetic_limits(), || {
            calls += 1;
            calls > 16
        }),
        Err(IonizedGasError::Cancelled)
    ));
    // A failed preparation has no mutable source or cached partial factor.
    owner
        .prepare_magnetized(1e10, [0.0, 0.0, 0.3], magnetic_limits(), || false)
        .unwrap();
}
