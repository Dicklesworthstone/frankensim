//! G0/G3: exact volume preservation and work-conjugate incompressible reduction.

use fs_material::MaterialError;
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::{FibreInvariants, Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};

fn nh(mu: f64, lambda: f64) -> Hyperelastic {
    Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 8.0).unwrap()
}

fn tilted_fibre(scale: f64, invariants: FibreInvariants) -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::HolzapfelGasserOgden {
            c10: 5e4 * scale,
            k1: 1.2e5 * scale,
            k2: 1.1,
            dispersion: 0.05,
            fibres: vec![[0.8, 0.0, 0.6], [0.0, 0.6, 0.8]],
            invariants,
            volumetric: Volumetric::Quadratic { kappa: 8e5 * scale },
        },
        8.0,
    )
    .unwrap()
}

fn laws() -> Vec<Hyperelastic> {
    vec![
        nh(2e5, 8e5),
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6e4,
                c01: 4e4,
                kappa: 9e5,
            },
            8.0,
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
            8.0,
        )
        .unwrap(),
        tilted_fibre(1.0, FibreInvariants::Full),
        tilted_fibre(1.0, FibreInvariants::Isochoric),
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
fn g0_neo_hookean_matches_independent_incompressible_energy_stress_and_tangent() {
    let mu = 2e5;
    let law = nh(mu, 3.0);
    for a in [
        [1.0, 0.0, 0.0, 1.0],
        [1.35, 0.2, -0.07, 0.9],
        [1.0, 0.6, 0.0, 1.0],
        [2.0, 0.0, 0.0, 2.0],
        [0.5, 0.1, 0.0, 0.5],
    ] {
        let response = law
            .incompressible_plane_stress(&a, PlaneStressSettings::default())
            .unwrap();
        let h = 1.0 / (a[0] * a[3] - a[1] * a[2]);
        let inverse_transpose = [a[3] * h, -a[2] * h, -a[1] * h, a[0] * h];
        let energy = 0.5 * mu * (a.iter().map(|v| v * v).sum::<f64>() + h * h - 3.0);
        assert_eq!(response.deformation_gradient[8], h);
        assert_eq!(response.iterations, 0);
        close(det(&response.deformation_gradient), 1.0, 1.0, 5e-15);
        close(
            response.energy_density_pa,
            energy,
            energy.abs().max(mu),
            1e-13,
        );
        // The native NH remainder contains -mu*ln(J); its multiplier is
        // shifted by -mu relative to the unshifted mu*I1/2 convention.
        close(
            response.pressure_pa,
            mu * (h * h - 1.0),
            mu * (h * h + 1.0),
            1e-13,
        );
        for i in 0..4 {
            let expected = mu * (a[i] - h * h * inverse_transpose[i]);
            close(response.piola[i], expected, expected.abs().max(mu), 2e-13);
            for j in 0..4 {
                let (r, c, k, l) = (i / 2, i % 2, j / 2, j % 2);
                let tangent = mu * (if i == j { 1.0 } else { 0.0 })
                    + mu * h
                        * h
                        * (2.0 * inverse_transpose[i] * inverse_transpose[j]
                            + inverse_transpose[2 * r + l] * inverse_transpose[2 * k + c]);
                close(
                    response.tangent[i][j],
                    tangent,
                    tangent.abs().max(mu),
                    5e-13,
                );
            }
        }
    }
}

#[test]
fn g3_all_original_models_retain_work_and_consistent_tangent_after_elimination() {
    let a = [1.3, 0.22, 0.12, 1.16];
    let settings = PlaneStressSettings {
        relative_tolerance: 1e-12,
        ..PlaneStressSettings::default()
    };
    for law in laws() {
        let response = law.incompressible_plane_stress(&a, settings).unwrap();
        close(det(&response.deformation_gradient), 1.0, 1.0, 5e-15);
        for j in 0..4 {
            let step = 2e-6;
            let mut plus = a;
            let mut minus = a;
            plus[j] += step;
            minus[j] -= step;
            let plus = law
                .incompressible_plane_stress_state(&plus, settings)
                .unwrap();
            let minus = law
                .incompressible_plane_stress_state(&minus, settings)
                .unwrap();
            close(
                (plus.energy_density_pa - minus.energy_density_pa) / (2.0 * step),
                response.piola[j],
                2e5,
                2e-7,
            );
            for i in 0..4 {
                close(
                    (plus.piola[i] - minus.piola[i]) / (2.0 * step),
                    response.tangent[i][j],
                    2e5,
                    2e-6,
                );
                close(response.tangent[i][j], response.tangent[j][i], 2e5, 5e-12);
            }
        }
    }
}

#[test]
fn g0_tilted_fibres_relax_shears_and_pressure_closes_all_physical_tractions() {
    let a = [1.3, 0.22, 0.12, 1.16];
    for invariants in [FibreInvariants::Full, FibreInvariants::Isochoric] {
        let law = tilted_fibre(1.0, invariants);
        let state = law
            .incompressible_plane_stress_state(&a, PlaneStressSettings::default())
            .unwrap();
        let f = state.deformation_gradient;
        assert!(state.iterations > 0);
        assert!(f[2].abs() > 1e-3 && f[5].abs() > 1e-3);
        // Reconstruct physical P directly from the original remainder and
        // independent 3-D determinant gradient, without using the reducer.
        let (_, remainder) = fs_ad::gradient::<9>(f, |x| law.isochoric_energy(&x));
        let (_, cofactor) = fs_ad::gradient::<9>(f, |x| {
            x[0] * (x[4] * x[8] - x[5] * x[7]) - x[1] * (x[3] * x[8] - x[5] * x[6])
                + x[2] * (x[3] * x[7] - x[4] * x[6])
        });
        let p: [f64; 9] = core::array::from_fn(|i| remainder[i] - state.pressure_pa * cofactor[i]);
        for i in 0..9 {
            close(state.full_piola[i], p[i], 2e5, 2e-12);
        }
        for i in [2, 5, 8] {
            close(p[i], 0.0, 2e5, 5e-11);
        }
        // Cauchy symmetry is an independent objectivity consequence; it also
        // checks the lower transverse Piola entries and pressure cofactors.
        for i in 0..3 {
            for j in 0..3 {
                let sigma_ij = (0..3).map(|k| p[3 * i + k] * f[3 * j + k]).sum::<f64>();
                let sigma_ji = (0..3).map(|k| p[3 * j + k] * f[3 * i + k]).sum::<f64>();
                close(sigma_ij, sigma_ji, 2e5, 5e-12);
                if i == 2 || j == 2 {
                    close(sigma_ij, 0.0, 2e5, 5e-11);
                }
            }
        }
        let mut unrelaxed = f;
        unrelaxed[2] = 0.0;
        unrelaxed[5] = 0.0;
        assert!(state.energy_density_pa < law.isochoric_energy(&unrelaxed));
    }
}

#[test]
fn g0_bulk_parameters_do_not_enter_exact_incompressibility_even_at_extreme_scales() {
    let a = [1.3, 0.22, 0.12, 1.16];
    for original in laws() {
        let reference = original
            .incompressible_plane_stress(&a, PlaneStressSettings::default())
            .unwrap();
        for value in [1e-15, 1.0, 1e300] {
            let mut model = original.model.clone();
            match &mut model {
                HyperelasticModel::NeoHookean { lambda, .. } => *lambda = value,
                HyperelasticModel::MooneyRivlin { kappa, .. } => *kappa = value,
                HyperelasticModel::Ogden { volumetric, .. }
                | HyperelasticModel::HolzapfelGasserOgden { volumetric, .. } => {
                    *volumetric = Volumetric::LogSquared { kappa: value };
                }
            }
            let law = Hyperelastic::new(model, original.stretch_limit).unwrap();
            let result = law
                .incompressible_plane_stress(&a, PlaneStressSettings::default())
                .unwrap();
            assert_eq!(result.deformation_gradient, reference.deformation_gradient);
            assert_eq!(
                result.energy_density_pa.to_bits(),
                reference.energy_density_pa.to_bits()
            );
            assert_eq!(
                result.pressure_pa.to_bits(),
                reference.pressure_pa.to_bits()
            );
            assert_eq!(result.piola, reference.piola);
            assert_eq!(result.full_piola, reference.full_piola);
            assert_eq!(result.tangent, reference.tangent);
            assert_eq!(result.iterations, reference.iterations);
        }
    }
}

#[test]
fn g3_superposed_rotation_transforms_physical_stress_and_reduced_tangent() {
    let a = [1.3, 0.22, 0.12, 1.16];
    let (c, s) = (0.8, 0.6);
    let rotate = |x: [f64; 4]| {
        [
            c * x[0] - s * x[2],
            c * x[1] - s * x[3],
            s * x[0] + c * x[2],
            s * x[1] + c * x[3],
        ]
    };
    let q = [
        [c, 0.0, -s, 0.0],
        [0.0, c, 0.0, -s],
        [s, 0.0, c, 0.0],
        [0.0, s, 0.0, c],
    ];
    for law in laws() {
        let reference = law
            .incompressible_plane_stress(&a, PlaneStressSettings::default())
            .unwrap();
        let result = law
            .incompressible_plane_stress(&rotate(a), PlaneStressSettings::default())
            .unwrap();
        close(
            result.energy_density_pa,
            reference.energy_density_pa,
            2e5,
            2e-11,
        );
        close(result.pressure_pa, reference.pressure_pa, 2e5, 2e-11);
        for (actual, expected) in result.piola.into_iter().zip(rotate(reference.piola)) {
            close(actual, expected, 2e5, 2e-11);
        }
        let f = reference.deformation_gradient;
        close(
            result.deformation_gradient[2],
            c * f[2] - s * f[5],
            1.0,
            2e-11,
        );
        close(
            result.deformation_gradient[5],
            s * f[2] + c * f[5],
            1.0,
            2e-11,
        );
        close(result.deformation_gradient[8], f[8], 1.0, 2e-14);
        for i in 0..4 {
            for j in 0..4 {
                let mut expected = 0.0;
                for k in 0..4 {
                    for l in 0..4 {
                        expected += q[i][k] * reference.tangent[k][l] * q[j][l];
                    }
                }
                close(result.tangent[i][j], expected, 2e5, 5e-10);
            }
        }
    }
}

#[test]
fn g3_modulus_scaling_preserves_thickness_and_relaxed_shears() {
    let a = [1.3, 0.22, 0.12, 1.16];
    let reference = tilted_fibre(1.0, FibreInvariants::Full)
        .incompressible_plane_stress(&a, PlaneStressSettings::default())
        .unwrap();
    for scale in [1e-12, 1e12] {
        let result = tilted_fibre(scale, FibreInvariants::Full)
            .incompressible_plane_stress(&a, PlaneStressSettings::default())
            .unwrap();
        for (actual, expected) in result
            .deformation_gradient
            .into_iter()
            .zip(reference.deformation_gradient)
        {
            close(actual, expected, 1.0, 1e-11);
        }
        close(
            result.energy_density_pa / scale,
            reference.energy_density_pa,
            2e5,
            2e-11,
        );
        close(
            result.pressure_pa / scale,
            reference.pressure_pa,
            2e5,
            2e-11,
        );
        for i in 0..9 {
            close(
                result.full_piola[i] / scale,
                reference.full_piola[i],
                2e5,
                2e-11,
            );
        }
        for i in 0..4 {
            for j in 0..4 {
                close(
                    result.tangent[i][j] / scale,
                    reference.tangent[i][j],
                    2e5,
                    2e-10,
                );
            }
        }
    }
}

#[test]
fn g0_state_and_tangent_paths_share_the_same_equilibrium() {
    for law in laws() {
        let a = [1.3, 0.22, 0.12, 1.16];
        let state = law
            .incompressible_plane_stress_state(&a, PlaneStressSettings::default())
            .unwrap();
        let full = law
            .incompressible_plane_stress(&a, PlaneStressSettings::default())
            .unwrap();
        assert_eq!(state.deformation_gradient, full.deformation_gradient);
        assert_eq!(
            state.energy_density_pa.to_bits(),
            full.energy_density_pa.to_bits()
        );
        assert_eq!(state.pressure_pa.to_bits(), full.pressure_pa.to_bits());
        assert_eq!(state.full_piola, full.full_piola);
        assert_eq!(state.piola, full.piola);
        assert_eq!(
            state.relative_transverse_residual.to_bits(),
            full.relative_transverse_residual.to_bits()
        );
        assert_eq!(state.iterations, full.iterations);
    }
}

#[test]
fn g0_inversion_exact_thickness_domain_and_exhausted_budgets_refuse() {
    let law = nh(2e5, 8e5);
    for a in [
        [-1.0, 0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0, 0.0],
        [f64::NAN, 0.0, 0.0, 1.0],
        [9.0, 0.0, 0.0, 1.0],
        [4.0, 0.0, 0.0, 4.0],
    ] {
        assert!(matches!(
            law.incompressible_plane_stress(&a, PlaneStressSettings::default()),
            Err(MaterialError::State { .. })
        ));
    }
    // The last state has admissible in-plane stretches, but exact J=1 would
    // require F22=1/16 outside [1/8,8]. Thickness cannot be clipped.
    assert!(
        law.piola(&[4.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.0])
            .is_ok()
    );
    for settings in [
        PlaneStressSettings {
            relative_tolerance: f64::NAN,
            ..PlaneStressSettings::default()
        },
        PlaneStressSettings {
            relative_tolerance: 0.0,
            ..PlaneStressSettings::default()
        },
        PlaneStressSettings {
            max_iterations: 0,
            ..PlaneStressSettings::default()
        },
        PlaneStressSettings {
            max_backtracks: 0,
            ..PlaneStressSettings::default()
        },
    ] {
        assert!(matches!(
            law.incompressible_plane_stress(&[1.0, 0.0, 0.0, 1.0], settings),
            Err(MaterialError::Parameters { .. })
        ));
    }
    assert!(
        tilted_fibre(1.0, FibreInvariants::Full)
            .incompressible_plane_stress(
                &[1.3, 0.22, 0.12, 1.16],
                PlaneStressSettings {
                    max_iterations: 1,
                    ..PlaneStressSettings::default()
                }
            )
            .is_err()
    );
    let mut invalid_domain = law;
    invalid_domain.stretch_limit = f64::NAN;
    assert!(matches!(
        invalid_domain
            .incompressible_plane_stress(&[1.0, 0.0, 0.0, 1.0], PlaneStressSettings::default()),
        Err(MaterialError::Parameters { .. })
    ));
}
