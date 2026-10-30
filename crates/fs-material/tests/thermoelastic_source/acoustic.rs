//! G1/G3: independent constant-entropy perturbations of the source potential.
use super::*;
use fs_blake3::ContentHash;
use fs_material::elastic::thermal::acoustic::AdiabaticThermoelasticInput;

fn input(strain: [f64; 6], capacity: f64) -> AdiabaticThermoelasticInput {
    AdiabaticThermoelasticInput {
        total_strain_mandel: strain,
        reference_density_kg_m3: 6.0,
        caloric_heat_capacity_j_m3_k: capacity,
        mechanical_material_identity: ContentHash([31; 32]),
        caloric_state_identity: ContentHash([41; 32]),
    }
}

#[test]
fn g1_unstressed_thermal_expansion_changes_compression_waves_but_not_shear() {
    for alpha in [-0.0005, 0.0, 0.0005] {
        let card = source_card(0.0, alpha, 0.0, 0.0, false);
        let at = law(&card).at_temperature(300.0).unwrap();
        let response = at.adiabatic_acoustics(input([0.0; 6], 1.2)).unwrap();
        let acoustic = response.acoustic();
        let bulk = 1000.0 / (3.0 * (1.0 - 2.0 * 0.2));
        let shear = 1000.0 / (2.0 * (1.0 + 0.2));
        // Independent isotropic thermodynamic identity K_S/K_T = C_sigma/C_epsilon.
        let fixed_strain_capacity = 1.2 - 9.0 * 300.0 * bulk * alpha * alpha;
        let adiabatic_bulk = bulk * 1.2 / fixed_strain_capacity;
        close(
            response.fixed_strain_heat_capacity_j_m3_k(),
            fixed_strain_capacity,
            3e-15,
        );
        close(acoustic.bulk_modulus_pa(), adiabatic_bulk, 2e-12);
        for direction in [[1.0, 0.0, 0.0], [0.3, -0.7, 0.2]] {
            let waves = acoustic.directional_waves(direction).unwrap();
            for mode in &waves.modes[..2] {
                close(mode.phase_speed_m_s, (shear / 6.0).sqrt(), 3e-13);
            }
            close(
                waves.modes[2].phase_speed_m_s,
                ((adiabatic_bulk + 4.0 * shear / 3.0) / 6.0).sqrt(),
                3e-13,
            );
        }
        let pressure = acoustic.hydrostatic_response(0.5).unwrap();
        close(pressure.volumetric_strain, -0.5 / adiabatic_bulk, 1e-17);
        close(
            pressure.stored_energy_density_j_m3,
            0.5 * 0.5_f64.powi(2) / adiabatic_bulk,
            1e-18,
        );
        assert_eq!(response.source().elastic, at.elastic);
        assert_eq!(response.source().expansion, at.expansion);
        if alpha == 0.0 {
            assert_eq!(
                *acoustic.stiffness_mandel_pa(),
                response.isothermal_stiffness_mandel_pa()
            );
            assert_eq!(*response.temperature_gradient_mandel_k(), [0.0; 6]);
        } else {
            assert!(acoustic.bulk_modulus_pa() > bulk);
            assert_eq!(
                response.temperature_gradient_mandel_k()[0].signum(),
                -alpha.signum()
            );
        }
    }
}

// Direct stress from independently authored E(T), nu(T) and integrated alpha(T).
fn authored_stress(t: f64, total: [f64; 6], parameters: [f64; 4]) -> [f64; 6] {
    let [young_slope, nu_slope, alpha, alpha_slope] = parameters;
    let dt = t - 300.0;
    let e = 1000.0 + young_slope * dt;
    let nu = 0.2 + nu_slope * dt;
    let free = alpha * dt + 0.5 * alpha_slope * dt * dt;
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let trace = total[..3].iter().sum::<f64>() - 3.0 * free;
    core::array::from_fn(|k| {
        if k < 3 {
            lambda * trace + 2.0 * mu * (total[k] - free)
        } else {
            2.0 * mu * total[k]
        }
    })
}

fn independent_entropy(t: f64, strain: [f64; 6], parameters: [f64; 4], capacity: f64) -> f64 {
    let [young_slope, nu_slope, alpha, alpha_slope] = parameters;
    let energy = |temperature| {
        varying_energy(
            temperature,
            strain,
            young_slope,
            nu_slope,
            alpha,
            alpha_slope,
            1.0,
        )
    };
    let step = 0.01;
    let derivative = (energy(t - 2.0 * step) - 8.0 * energy(t - step) + 8.0 * energy(t + step)
        - energy(t + 2.0 * step))
        / (12.0 * step);
    capacity * (t / 300.0).ln() - derivative
}

#[test]
fn g1_full_tangent_matches_stress_differences_on_an_independent_isentrope() {
    let strain = [0.013, -0.004, 0.002, 0.007, -0.003, 0.005];
    let temperature = 321.0;
    let capacity = 1.2;
    for parameters in [
        [0.0, 0.0, 0.0003, 0.0],
        [2.0, 0.0008, -0.0001, 0.000003],
        [-2.0, -0.0007, 0.0001, -0.000002],
        [2.0, 0.001, 0.0, 0.0],
    ] {
        let [young_slope, nu_slope, alpha, alpha_slope] = parameters;
        let card = source_card(young_slope, alpha, nu_slope, alpha_slope, false);
        let at = law(&card).at_temperature(temperature).unwrap();
        let response = at.adiabatic_acoustics(input(strain, capacity)).unwrap();
        let entropy = independent_entropy(temperature, strain, parameters, capacity);
        let step = 1e-5;
        for column in 0..6 {
            let mut stress = [[0.0; 6]; 2];
            let mut temperatures = [0.0; 2];
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                let mut displaced = strain;
                displaced[column] += sign * step;
                let mut low = temperature - 1.0;
                let mut high = temperature + 1.0;
                assert!(independent_entropy(low, displaced, parameters, capacity) < entropy);
                assert!(independent_entropy(high, displaced, parameters, capacity) > entropy);
                for _ in 0..64 {
                    let midpoint = 0.5 * (low + high);
                    if independent_entropy(midpoint, displaced, parameters, capacity) < entropy {
                        low = midpoint;
                    } else {
                        high = midpoint;
                    }
                }
                temperatures[side] = 0.5 * (low + high);
                stress[side] = authored_stress(temperatures[side], displaced, parameters);
            }
            close(
                (temperatures[1] - temperatures[0]) / (2.0 * step),
                response.temperature_gradient_mandel_k()[column],
                2e-5,
            );
            for row in 0..6 {
                close(
                    (stress[1][row] - stress[0][row]) / (2.0 * step),
                    response.acoustic().stiffness_mandel_pa()[row][column],
                    2e-4,
                );
            }
        }
        // The native mixed derivative must agree with stress at fixed TOTAL
        // strain. This also detects a wrong shear convention or wrong sign.
        let low = authored_stress(temperature - 0.001, strain, parameters);
        let high = authored_stress(temperature + 0.001, strain, parameters);
        for k in 0..6 {
            close(
                (high[k] - low[k]) / 0.002,
                response.stress_temperature_derivative_mandel_pa_k()[k],
                2e-10,
            );
        }
    }
}

#[test]
fn g3_caloric_capacity_and_thermal_source_refresh_change_native_wave_response() {
    let card = source_card(2.0, 0.0004, 0.0003, 0.000001, false);
    let law = law(&card);
    let strain = [0.01, -0.002, 0.004, 0.003, 0.0, -0.001];
    let cold = law
        .at_temperature(300.0)
        .unwrap()
        .adiabatic_acoustics(input(strain, 1.2))
        .unwrap();
    let hot = law
        .at_temperature(330.0)
        .unwrap()
        .adiabatic_acoustics(input(strain, 1.2))
        .unwrap();
    let large = law
        .at_temperature(300.0)
        .unwrap()
        .adiabatic_acoustics(input(strain, 2.4))
        .unwrap();
    assert_ne!(
        cold.acoustic().material_state_identity(),
        hot.acoustic().material_state_identity()
    );
    assert_ne!(
        cold.acoustic().material_state_identity(),
        large.acoustic().material_state_identity()
    );
    let direction = [1.0, 0.3, -0.1];
    let speed = |state: &fs_material::elastic::thermal::acoustic::AdiabaticThermoelasticState| {
        state.acoustic().directional_waves(direction).unwrap().modes[2].phase_speed_m_s
    };
    assert!((speed(&hot) - speed(&cold)).abs() > 0.01);
    assert!(speed(&large) < speed(&cold));
    for (a, b) in cold
        .source()
        .elastic
        .properties()
        .iter()
        .zip(hot.source().elastic.properties())
    {
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
}

#[test]
fn g4_unstable_capacity_and_invalid_reference_inputs_refuse() {
    let card = source_card(0.0, 0.0005, 0.0, 0.0, false);
    let at = law(&card).at_temperature(300.0).unwrap();
    let accepted = at.adiabatic_acoustics(input([0.0; 6], 1.2)).unwrap();
    let elastic_capacity = at.state([0.0; 6], 1.0).unwrap().heat_capacity_j_k;
    assert!(
        at.adiabatic_acoustics(input([0.0; 6], -elastic_capacity))
            .is_err()
    );
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut values = input([0.0; 6], 1.2);
        values.reference_density_kg_m3 = invalid;
        assert!(at.adiabatic_acoustics(values).is_err());
        values = input([0.0; 6], invalid);
        assert!(at.adiabatic_acoustics(values).is_err());
    }
    let mut values = input([0.0; 6], 1.2);
    values.caloric_state_identity = ContentHash([0; 32]);
    assert!(at.adiabatic_acoustics(values).is_err());
    values = input([0.0; 6], 1.2);
    values.mechanical_material_identity = ContentHash([0; 32]);
    assert!(at.adiabatic_acoustics(values).is_err());
    assert!(at.adiabatic_acoustics(input([0.051; 6], 1.2)).is_err());
    assert_eq!(
        at.adiabatic_acoustics(input([0.0; 6], 1.2)).unwrap(),
        accepted
    );
}
