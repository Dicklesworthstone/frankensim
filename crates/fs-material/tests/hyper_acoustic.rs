//! G0/G3: small waves retain the original finite-strain tangent and signed modes.

use fs_material::MaterialError;
use fs_material::hyper::acoustic::HyperelasticAcousticState;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::hyper::{FibreInvariants, Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};

const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

fn nh(mu: f64, lambda: f64, limit: f64) -> Hyperelastic {
    Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, limit).unwrap()
}

fn laws() -> Vec<Hyperelastic> {
    vec![
        nh(2e5, 8e5, 4.0),
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6e4,
                c01: 4e4,
                kappa: 9e5,
            },
            4.0,
        )
        .unwrap(),
        Hyperelastic::new(
            HyperelasticModel::Ogden {
                terms: vec![
                    OgdenTerm {
                        mu_pa: 2e5,
                        alpha: 1.3,
                    },
                    OgdenTerm {
                        mu_pa: -1e4,
                        alpha: -2.0,
                    },
                ],
                volumetric: Volumetric::LogSquared { kappa: 8e5 },
            },
            4.0,
        )
        .unwrap(),
        Hyperelastic::new(
            HyperelasticModel::HolzapfelGasserOgden {
                c10: 5e4,
                k1: 1.2e5,
                k2: 1.1,
                dispersion: 0.05,
                fibres: vec![[0.8, 0.0, 0.6], [0.0, 0.6, 0.8]],
                invariants: FibreInvariants::Full,
                volumetric: Volumetric::Quadratic { kappa: 8e5 },
            },
            4.0,
        )
        .unwrap(),
    ]
}

fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} vs {b:e}; scale {scale:e}, tolerance {tolerance:e}"
    );
}

fn det(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

#[test]
fn g0_reference_waves_match_small_strain_longitudinal_and_shear_speeds() {
    let (e, nu, rho) = (7e5, 0.3, 1050.0);
    let mu = e / (2.0 * (1.0 + nu));
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let state = nh(mu, lambda, 4.0).acoustic_state(&IDENTITY, rho).unwrap();
    assert_eq!(state.deformation_gradient(), &IDENTITY);
    assert_eq!(state.reference_density_kg_m3(), rho);
    assert_eq!(state.current_density_kg_m3(), rho);
    for direction in [[1.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-2e-300, 3e-300, 1e-300]] {
        let waves = state.directional_waves(direction).unwrap();
        for (mode, expected) in
            waves
                .modes
                .iter()
                .zip([mu / rho, mu / rho, (lambda + 2.0 * mu) / rho])
        {
            close(mode.squared_phase_speed_m2_s2, expected, expected, 2e-13);
            close(
                mode.phase_speed_m_s.unwrap(),
                expected.sqrt(),
                expected.sqrt(),
                2e-13,
            );
            close(
                mode.directional_modulus_pa,
                rho * expected,
                rho * expected,
                2e-13,
            );
            for (velocity, normal) in mode
                .energy_velocity_m_s
                .unwrap()
                .into_iter()
                .zip(waves.direction)
            {
                close(velocity, expected.sqrt() * normal, expected.sqrt(), 3e-13);
            }
        }
        for mode in &waves.modes[..2] {
            close(
                (0..3)
                    .map(|i| mode.polarization[i] * waves.direction[i])
                    .sum(),
                0.0,
                1.0,
                2e-13,
            );
        }
        close(
            (0..3)
                .map(|i| waves.modes[2].polarization[i] * waves.direction[i])
                .sum::<f64>()
                .abs(),
            1.0,
            1.0,
            2e-13,
        );
    }
}

#[test]
fn g0_material_factory_accepts_trait_objects_and_preserves_native_law_results() {
    let f = [1.2, 0.15, 0.03, -0.04, 0.93, 0.08, 0.02, -0.05, 1.08];
    for law in laws() {
        let material: &dyn HyperelasticMaterial = &law;
        let generic = HyperelasticAcousticState::from_material(material, &f, 1000.0).unwrap();
        let original = law.acoustic_state(&f, 1000.0).unwrap();
        assert_eq!(generic, original);
        assert_eq!(
            generic.directional_waves([0.4, -0.3, 0.8]).unwrap(),
            original.directional_waves([0.4, -0.3, 0.8]).unwrap(),
        );
        for invalid in [[0.0; 9], [f64::NAN; 9], [f64::MAX; 9]] {
            assert!(matches!(
                HyperelasticAcousticState::from_material(material, &invalid, 1000.0),
                Err(MaterialError::State { .. }),
            ));
        }
    }
}

#[test]
fn g0_finite_neo_hookean_tensor_matches_independent_closed_form_and_current_mass() {
    let (mu, lambda, rho0) = (2e5, 8e5, 1000.0);
    let law = nh(mu, lambda, 4.0);
    for f in [
        [1.2, 0.15, 0.03, -0.04, 0.93, 0.08, 0.02, -0.05, 1.08],
        [0.65, 0.0, 0.0, 0.0, 1.1, 0.0, 0.0, 0.0, 1.2],
        [1.7, 0.0, 0.0, 0.0, 1.7, 0.0, 0.0, 0.0, 1.7],
    ] {
        let state = law.acoustic_state(&f, rho0).unwrap();
        let waves = state.directional_waves([0.4, -0.3, 0.8]).unwrap();
        let n = waves.direction;
        let q: [f64; 3] = core::array::from_fn(|j| (0..3).map(|i| f[3 * i + j] * n[i]).sum());
        let nbn: f64 = q.iter().map(|v| v * v).sum();
        let j = det(&f);
        let beta = lambda + mu - lambda * j.ln();
        close(state.current_density_kg_m3(), rho0 / j, rho0 / j, 1e-14);
        for i in 0..3 {
            for k in 0..3 {
                let expected = (if i == k { mu * nbn } else { 0.0 }) / j + beta * n[i] * n[k] / j;
                close(waves.acoustic_tensor_pa[i][k], expected, 2e6, 5e-13);
            }
        }
        let mut expected = [mu * nbn / rho0, mu * nbn / rho0, (mu * nbn + beta) / rho0];
        expected.sort_by(f64::total_cmp);
        for (mode, value) in waves.modes.iter().zip(expected) {
            close(mode.squared_phase_speed_m2_s2, value, 2e3, 5e-13);
            close(
                mode.directional_modulus_pa,
                waves.current_density_kg_m3 * value,
                2e6,
                5e-13,
            );
            let bn: [f64; 3] = core::array::from_fn(|i| (0..3).map(|j| f[3 * i + j] * q[j]).sum());
            let longitudinal = (0..3)
                .map(|i| mode.polarization[i] * n[i])
                .sum::<f64>()
                .abs()
                > 0.5;
            let speed = mode.phase_speed_m_s.unwrap();
            for (i, velocity) in mode.energy_velocity_m_s.unwrap().into_iter().enumerate() {
                let expected =
                    (mu * bn[i] + if longitudinal { beta * n[i] } else { 0.0 }) / (rho0 * speed);
                close(velocity, expected, speed, 2e-12);
            }
        }
    }
}

#[test]
fn g3_all_original_laws_match_rank_one_stress_and_energy_variations() {
    let f = [1.2, 0.15, 0.03, -0.04, 0.93, 0.08, 0.02, -0.05, 1.08];
    for law in laws() {
        let state = law.acoustic_state(&f, 1000.0).unwrap();
        for direction in [[1.0, 0.2, -0.3], [-0.2, 0.7, 0.5]] {
            let waves = state.directional_waves(direction).unwrap();
            let q: [f64; 3] =
                core::array::from_fn(|j| (0..3).map(|i| f[3 * i + j] * waves.direction[i]).sum());
            for polarization in [[1.0, 0.0, 0.0], [0.2, -0.4, 0.7]] {
                let increment: [f64; 9] =
                    core::array::from_fn(|ij| polarization[ij / 3] * q[ij % 3]);
                let step = 2e-6;
                let plus = core::array::from_fn(|i| f[i] + step * increment[i]);
                let minus = core::array::from_fn(|i| f[i] - step * increment[i]);
                let pp = law.piola(&plus).unwrap();
                let pm = law.piola(&minus).unwrap();
                for i in 0..3 {
                    let derivative = (0..3)
                        .map(|j| (pp[3 * i + j] - pm[3 * i + j]) * q[j])
                        .sum::<f64>()
                        / (2.0 * step * det(&f));
                    let qa = (0..3)
                        .map(|k| waves.acoustic_tensor_pa[i][k] * polarization[k])
                        .sum();
                    close(derivative, qa, 2e6, 2e-7);
                }
                let step = 2e-4;
                let plus = core::array::from_fn(|i| f[i] + step * increment[i]);
                let minus = core::array::from_fn(|i| f[i] - step * increment[i]);
                let curvature = (law.energy(&plus) - 2.0 * law.energy(&f) + law.energy(&minus))
                    / (step * step * det(&f));
                let expected: f64 = (0..3)
                    .map(|i| {
                        (0..3)
                            .map(|k| {
                                polarization[i] * waves.acoustic_tensor_pa[i][k] * polarization[k]
                            })
                            .sum::<f64>()
                    })
                    .sum();
                close(curvature, expected, 2e6, 2e-6);
            }
            for mode in waves.modes {
                for i in 0..3 {
                    let qa: f64 = (0..3)
                        .map(|k| waves.acoustic_tensor_pa[i][k] * mode.polarization[k])
                        .sum();
                    close(
                        qa,
                        mode.directional_modulus_pa * mode.polarization[i],
                        2e6,
                        2e-12,
                    );
                }
            }
        }
    }
}

#[test]
fn g3_spatial_rotation_preserves_speeds_and_rotates_the_full_anisotropic_tensor() {
    let f = [1.2, 0.15, 0.03, -0.04, 0.93, 0.08, 0.02, -0.05, 1.08];
    let r = [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]];
    let rotate = |v: [f64; 3]| core::array::from_fn(|i| (0..3).map(|j| r[i][j] * v[j]).sum());
    let rf = core::array::from_fn(|ij| (0..3).map(|k| r[ij / 3][k] * f[3 * k + ij % 3]).sum());
    let n = [0.4, -0.3, 0.8];
    for law in laws() {
        let a = law
            .acoustic_state(&f, 1000.0)
            .unwrap()
            .directional_waves(n)
            .unwrap();
        let b = law
            .acoustic_state(&rf, 1000.0)
            .unwrap()
            .directional_waves(rotate(n))
            .unwrap();
        close(
            a.current_density_kg_m3,
            b.current_density_kg_m3,
            1000.0,
            1e-13,
        );
        for i in 0..3 {
            for j in 0..3 {
                let mut expected = 0.0;
                for k in 0..3 {
                    for l in 0..3 {
                        expected += r[i][k] * a.acoustic_tensor_pa[k][l] * r[j][l];
                    }
                }
                close(b.acoustic_tensor_pa[i][j], expected, 2e6, 2e-12);
            }
        }
        for i in 0..3 {
            close(
                a.modes[i].squared_phase_speed_m2_s2,
                b.modes[i].squared_phase_speed_m2_s2,
                2e3,
                2e-12,
            );
            let expected_velocity = rotate(a.modes[i].energy_velocity_m_s.unwrap());
            for (actual, expected) in b.modes[i]
                .energy_velocity_m_s
                .unwrap()
                .into_iter()
                .zip(expected_velocity)
            {
                close(actual, expected, a.modes[i].phase_speed_m_s.unwrap(), 3e-11);
            }
            // Projected eigenvectors remain valid even when a repeated root
            // permits the shared Jacobi solver to choose a different basis.
            let p = rotate(a.modes[i].polarization);
            for j in 0..3 {
                let qp: f64 = (0..3).map(|k| b.acoustic_tensor_pa[j][k] * p[k]).sum();
                close(qp, b.modes[i].directional_modulus_pa * p[j], 2e6, 2e-12);
            }
        }
    }
}

#[test]
fn g0_admitted_expansion_keeps_negative_longitudinal_squared_speed() {
    let law = nh(2.0, 40.0, 4.0);
    let f = [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0];
    let waves = law
        .acoustic_state(&f, 1000.0)
        .unwrap()
        .directional_waves([1.0, 0.0, 0.0])
        .unwrap();
    let expected = (10.0 + 40.0 * (1.0 - 8.0_f64.ln())) / 1000.0;
    assert!(expected < 0.0);
    close(
        waves.modes[0].squared_phase_speed_m2_s2,
        expected,
        0.1,
        2e-13,
    );
    assert_eq!(waves.modes[0].phase_speed_m_s, None);
    assert_eq!(waves.modes[0].energy_velocity_m_s, None);
    assert!(waves.modes[1].phase_speed_m_s.unwrap() > 0.0);
    assert!(waves.modes[2].phase_speed_m_s.unwrap() > 0.0);
    close(waves.modes[0].polarization[0].abs(), 1.0, 1.0, 1e-13);
    // Independent negative rank-one energy curvature rules out a global
    // polyconvexity declaration for this admitted material domain.
    let h = 1e-4;
    let mut plus = f;
    let mut minus = f;
    plus[0] += 2.0 * h;
    minus[0] -= 2.0 * h;
    let curvature = (law.energy(&plus) - 2.0 * law.energy(&f) + law.energy(&minus)) / (h * h);
    assert!(curvature < 0.0);
    close(curvature, expected * 1000.0, 100.0, 1e-6);
    assert_eq!(law.admissibility().polyconvex, None);
    assert_eq!(nh(2.0, 40.0, 1.1).admissibility().polyconvex, Some(true));
    assert_eq!(nh(2.0, 0.0, 10.0).admissibility().polyconvex, Some(true));
    for (c01, expected) in [(0.0, Some(true)), (1.0, None)] {
        let mr = Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 2.0,
                c01,
                kappa: 40.0,
            },
            4.0,
        )
        .unwrap();
        assert_eq!(mr.admissibility().polyconvex, expected);
    }
}

#[test]
fn g3_units_density_and_direction_scaling_preserve_physical_wave_predictions() {
    let f = [1.2, 0.15, 0.03, -0.04, 0.93, 0.08, 0.02, -0.05, 1.08];
    let reference = nh(2e5, 8e5, 4.0)
        .acoustic_state(&f, 1000.0)
        .unwrap()
        .directional_waves([1.0, -2.0, 3.0])
        .unwrap();
    for scale in [1e-12, 1e12] {
        let state = nh(2e5 * scale, 8e5 * scale, 4.0)
            .acoustic_state(&f, 1000.0 * scale)
            .unwrap();
        let result = state.directional_waves([1e300, -2e300, 3e300]).unwrap();
        for i in 0..3 {
            close(
                result.modes[i].squared_phase_speed_m2_s2,
                reference.modes[i].squared_phase_speed_m2_s2,
                2e3,
                2e-12,
            );
            close(
                result.modes[i].directional_modulus_pa / scale,
                reference.modes[i].directional_modulus_pa,
                2e6,
                2e-12,
            );
            for (actual, expected) in result.modes[i]
                .energy_velocity_m_s
                .unwrap()
                .into_iter()
                .zip(reference.modes[i].energy_velocity_m_s.unwrap())
            {
                close(
                    actual,
                    expected,
                    reference.modes[i].phase_speed_m_s.unwrap(),
                    3e-12,
                );
            }
        }
    }
    let doubled_density = nh(2e5, 8e5, 4.0)
        .acoustic_state(&f, 2000.0)
        .unwrap()
        .directional_waves([-1.0, 2.0, -3.0])
        .unwrap();
    for i in 0..3 {
        close(
            doubled_density.modes[i].squared_phase_speed_m2_s2,
            0.5 * reference.modes[i].squared_phase_speed_m2_s2,
            2e3,
            1e-13,
        );
        for (actual, expected) in doubled_density.modes[i]
            .energy_velocity_m_s
            .unwrap()
            .into_iter()
            .zip(reference.modes[i].energy_velocity_m_s.unwrap())
        {
            close(
                actual,
                -expected / 2.0_f64.sqrt(),
                reference.modes[i].phase_speed_m_s.unwrap(),
                3e-12,
            );
        }
    }
}

#[test]
fn g3_deformed_energy_transport_matches_frequency_gradients_for_all_native_laws() {
    let f = [1.2, 0.15, 0.03, -0.04, 0.93, 0.08, 0.02, -0.05, 1.08];
    let k = [0.4, -0.3, 0.8];
    for law in laws() {
        let state = law.acoustic_state(&f, 1000.0).unwrap();
        let waves = state.directional_waves(k).unwrap();
        for mode in 0..3 {
            let speed = waves.modes[mode].phase_speed_m_s.unwrap();
            let velocity = waves.modes[mode].energy_velocity_m_s.unwrap();
            close(
                (0..3).map(|i| waves.direction[i] * velocity[i]).sum(),
                speed,
                speed,
                3e-12,
            );
            for axis in 0..3 {
                let h = 1e-5;
                let mut plus = k;
                let mut minus = k;
                plus[axis] += h;
                minus[axis] -= h;
                let frequency = |k: [f64; 3]| {
                    k.iter().map(|v| v * v).sum::<f64>().sqrt()
                        * state.directional_waves(k).unwrap().modes[mode]
                            .phase_speed_m_s
                            .unwrap()
                };
                let derivative = (frequency(plus) - frequency(minus)) / (2.0 * h);
                close(velocity[axis], derivative, speed, 2e-8);
            }
        }
    }
}

#[test]
fn g0_invalid_mass_directions_material_mutations_and_out_of_domain_states_refuse() {
    let law = nh(2e5, 8e5, 4.0);
    for rho in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            law.acoustic_state(&IDENTITY, rho),
            Err(MaterialError::Parameters { .. })
        ));
    }
    let state = law.acoustic_state(&IDENTITY, 1000.0).unwrap();
    for direction in [[0.0; 3], [f64::NAN, 0.0, 1.0], [f64::INFINITY, 1.0, 0.0]] {
        assert!(matches!(
            state.directional_waves(direction),
            Err(MaterialError::Parameters { .. })
        ));
    }
    for f in [
        [-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        [5.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
    ] {
        assert!(matches!(
            law.acoustic_state(&f, 1000.0),
            Err(MaterialError::State { .. })
        ));
    }
    let mut changed = law;
    changed.model = HyperelasticModel::NeoHookean {
        mu: -1.0,
        lambda: 8e5,
    };
    assert!(matches!(
        changed.acoustic_state(&IDENTITY, 1000.0),
        Err(MaterialError::Parameters { .. })
    ));
    // The previously prepared state is an immutable snapshot of its admitted law.
    assert!(
        state
            .directional_waves([1.0, 0.0, 0.0])
            .unwrap()
            .modes
            .iter()
            .all(|m| m.phase_speed_m_s.is_some())
    );
}
