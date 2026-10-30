use super::*;
use crate::elastic::{IsotropicElastic, OrthotropicElastic};

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()),
        "actual {actual:e}, expected {expected:e}"
    );
}

fn identity() -> ContentHash {
    ContentHash([7; 32])
}

#[test]
fn isotropic_waves_and_pressure_share_one_elastic_energy() {
    let law = IsotropicElastic::new(210e9, 0.3, 0.001).unwrap();
    let rho = 7800.0;
    let state = ElasticAcousticState::from_isotropic(&law, rho, identity()).unwrap();
    let bulk = law.youngs / (3.0 * (1.0 - 2.0 * law.poisson));
    let shear = law.youngs / (2.0 * (1.0 + law.poisson));
    assert_eq!(state.material_state_identity(), identity());
    assert_eq!(state.stiffness_mandel_pa(), &law.stiffness());
    assert_eq!(state.maximum_absolute_mandel_strain(), law.strain_limit);
    near(state.bulk_modulus_pa(), bulk, 1e-14);
    near(state.compressibility_pa_inv(), 1.0 / bulk, 1e-14);
    for direction in [
        [1.0, 0.0, 0.0],
        [1.0, 2.0, -3.0],
        [f64::MAX, 0.0, 0.0],
        [1e-300, 0.0, 0.0],
    ] {
        let waves = state.directional_waves(direction).unwrap();
        assert_eq!(waves.material_state_identity, identity());
        for (index, modulus) in [shear, shear, bulk + 4.0 * shear / 3.0]
            .into_iter()
            .enumerate()
        {
            let mode = waves.modes[index];
            near(mode.directional_modulus_pa, modulus, 1e-14);
            near(mode.phase_speed_m_s, (modulus / rho).sqrt(), 1e-14);
            near(mode.impedance_pa_s_m, rho * mode.phase_speed_m_s, 1e-14);
            let alignment: f64 = (0..3)
                .map(|i| mode.polarization[i] * waves.direction[i])
                .sum();
            if index < 2 {
                assert!(alignment.abs() < 1e-14);
            } else {
                near(alignment.abs(), 1.0, 1e-14);
            }
        }
    }
    for pressure in [-1e7, 0.0, 1e7] {
        let response = state.hydrostatic_response(pressure).unwrap();
        assert_eq!(response.material_state_identity, identity());
        for value in &response.strain_mandel[..3] {
            near(*value, -pressure / (3.0 * bulk), 1e-14);
        }
        assert_eq!(&response.strain_mandel[3..], &[0.0; 3]);
        near(response.volumetric_strain, -pressure / bulk, 1e-14);
        near(
            response.stored_energy_density_j_m3,
            pressure * pressure / (2.0 * bulk),
            1e-14,
        );
    }
    assert!(
        state
            .hydrostatic_response(4.0 * bulk * law.strain_limit)
            .is_err()
    );
}

#[test]
fn orthotropic_principal_and_oblique_wave_speeds_match_closed_form() {
    let law =
        OrthotropicElastic::new([120e9, 60e9, 30e9], [0.0; 3], [20e9, 10e9, 15e9], 0.002).unwrap();
    let state = ElasticAcousticState::from_orthotropic(&law, 5000.0, identity()).unwrap();
    near(
        state.compressibility_pa_inv(),
        law.e.iter().map(|e| 1.0 / e).sum(),
        1e-14,
    );
    let response = state.hydrostatic_response(1e6).unwrap();
    for i in 0..3 {
        near(response.strain_mandel[i], -1e6 / law.e[i], 1e-14);
    }
    let principal = state.directional_waves([1.0, 0.0, 0.0]).unwrap();
    for (mode, modulus) in principal.modes.iter().zip([15e9, 20e9, 120e9]) {
        near(mode.directional_modulus_pa, modulus, 1e-14);
    }
    let direction = [0.6, 0.8, 0.0];
    let waves = state.directional_waves(direction).unwrap();
    let [nx, ny, _] = direction;
    let a = 120e9 * nx * nx + 20e9 * ny * ny;
    let d = 20e9 * nx * nx + 60e9 * ny * ny;
    let b = 20e9 * nx * ny;
    let discriminant = ((a - d) * (a - d) + 4.0 * b * b).sqrt();
    let mut expected = [
        0.5 * (a + d - discriminant),
        0.5 * (a + d + discriminant),
        15e9 * nx * nx + 10e9 * ny * ny,
    ];
    expected.sort_by(f64::total_cmp);
    for (mode, modulus) in waves.modes.iter().zip(expected) {
        near(mode.directional_modulus_pa, modulus, 1e-14);
    }
    let heavier = ElasticAcousticState::from_orthotropic(&law, 20000.0, identity())
        .unwrap()
        .directional_waves(direction)
        .unwrap();
    for (light, heavy) in waves.modes.iter().zip(heavier.modes) {
        near(heavy.phase_speed_m_s, 0.5 * light.phase_speed_m_s, 1e-14);
        near(heavy.impedance_pa_s_m, 2.0 * light.impedance_pa_s_m, 1e-14);
    }
}

#[test]
fn full_tensor_rotation_preserves_hydrostatic_work_and_wave_modes() {
    let matrix: Tangent6 = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            if i == j {
                (10 + i) as f64 * 1e9
            } else {
                (i + j + 1) as f64 * 1e8
            }
        })
    });
    let q = crate::tensor::rotation([1.0 / 3.0_f64.sqrt(); 3], 0.47);
    let rotated = crate::tensor::rotate_mandel(matrix, q);
    let a = ElasticAcousticState::try_new_mandel(matrix, 3000.0, identity(), 0.001).unwrap();
    let b = ElasticAcousticState::try_new_mandel(rotated, 3000.0, identity(), 0.001).unwrap();
    near(a.bulk_modulus_pa(), b.bulk_modulus_pa(), 1e-13);
    let pressure_a = a.hydrostatic_response(1e5).unwrap();
    let pressure_b = b.hydrostatic_response(1e5).unwrap();
    assert!(
        pressure_a.strain_mandel[3..]
            .iter()
            .any(|value| value.abs() > 1e-8)
    );
    let mandel_rotation = crate::tensor::mandel_rotation(q);
    for (i, row) in mandel_rotation.iter().enumerate() {
        let expected = (0..6).map(|j| row[j] * pressure_a.strain_mandel[j]).sum();
        near(pressure_b.strain_mandel[i], expected, 1e-13);
    }
    near(
        pressure_a.stored_energy_density_j_m3,
        pressure_b.stored_energy_density_j_m3,
        1e-13,
    );
    let direction = [0.6, 0.0, 0.8];
    let rotated_direction = core::array::from_fn(|i| (0..3).map(|j| q[i][j] * direction[j]).sum());
    let waves_a = a.directional_waves(direction).unwrap();
    let waves_b = b.directional_waves(rotated_direction).unwrap();
    for (mode_a, mode_b) in waves_a.modes.iter().zip(waves_b.modes) {
        near(mode_a.phase_speed_m_s, mode_b.phase_speed_m_s, 1e-13);
    }
    for i in 0..3 {
        for j in 0..3 {
            let mut expected = 0.0;
            for k in 0..3 {
                for l in 0..3 {
                    expected += q[i][k] * waves_a.christoffel_pa[k][l] * q[j][l];
                }
            }
            near(waves_b.christoffel_pa[i][j], expected, 1e-13);
        }
    }
}

#[test]
fn invalid_source_tensors_directions_and_pressure_refuse() {
    let law = IsotropicElastic::new(100e9, 0.25, 0.001).unwrap();
    let state = ElasticAcousticState::from_isotropic(&law, 1000.0, identity()).unwrap();
    for direction in [[0.0; 3], [f64::NAN, 1.0, 0.0], [f64::INFINITY, 0.0, 0.0]] {
        assert!(state.directional_waves(direction).is_err());
    }
    for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(ElasticAcousticState::from_isotropic(&law, density, identity()).is_err());
    }
    for limit in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            ElasticAcousticState::try_new_mandel(law.stiffness(), 1000.0, identity(), limit)
                .is_err()
        );
    }
    assert!(ElasticAcousticState::from_isotropic(&law, 1000.0, ContentHash([0; 32])).is_err());
    let mut asymmetric = law.stiffness();
    asymmetric[0][1] += 1.0;
    assert!(ElasticAcousticState::try_new_mandel(asymmetric, 1000.0, identity(), 0.001).is_err());
    let mut indefinite = law.stiffness();
    indefinite[0][0] = -1.0;
    assert!(ElasticAcousticState::try_new_mandel(indefinite, 1000.0, identity(), 0.001).is_err());
    for pressure in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(state.hydrostatic_response(pressure).is_err());
    }
}
