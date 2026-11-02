//! G0/G3: traction-free hyperelastic reduction and independent work derivatives.

use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::{FibreInvariants, Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};
use fs_material::{IsotropicElastic, MaterialError};

const IN_PLANE: [usize; 4] = [0, 1, 3, 4];
const TRANSVERSE: [usize; 3] = [2, 5, 8];

fn nh(youngs: f64, poisson: f64, limit: f64) -> Hyperelastic {
    let (lambda, mu) = IsotropicElastic::new(youngs, poisson, 1.0).unwrap().lame();
    Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, limit).unwrap()
}

fn tilted_fibre(scale: f64) -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::HolzapfelGasserOgden {
            c10: 5e4 * scale,
            k1: 1.2e5 * scale,
            k2: 1.1,
            dispersion: 0.05,
            fibres: vec![[0.8, 0.0, 0.6], [0.0, 0.6, 0.8]],
            invariants: FibreInvariants::Full,
            volumetric: Volumetric::Quadratic { kappa: 8e5 * scale },
        },
        4.0,
    )
    .unwrap()
}

fn laws() -> Vec<Hyperelastic> {
    vec![
        nh(5e5, 0.3, 4.0),
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
        tilted_fibre(1.0),
    ]
}

fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} vs {b:e}; scale {scale:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn g0_reference_tangent_is_isotropic_plane_stress_not_plane_strain() {
    let e = 3e5;
    for nu in [0.0, 0.25, 0.49] {
        let law = nh(e, nu, 4.0);
        let response = law
            .plane_stress(&[1.0, 0.0, 0.0, 1.0], PlaneStressSettings::default())
            .unwrap();
        let c = e / (1.0 - nu * nu);
        let mu = e / (2.0 * (1.0 + nu));
        let expected = [
            [c, 0.0, 0.0, c * nu],
            [0.0, mu, mu, 0.0],
            [0.0, mu, mu, 0.0],
            [c * nu, 0.0, 0.0, c],
        ];
        assert_eq!(response.iterations, 0);
        assert_eq!(response.deformation_gradient[8], 1.0);
        close(response.energy_density_pa, 0.0, e, 1e-14);
        for i in 0..4 {
            close(response.piola[i], 0.0, e, 1e-14);
            for j in 0..4 {
                close(response.tangent[i][j], expected[i][j], e, 2e-13);
            }
        }
    }
}

#[test]
fn g0_finite_neo_hookean_response_satisfies_the_independent_thickness_equation() {
    let (mu, lambda) = (2e5, 8e5);
    let law = Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0).unwrap();
    let f = [1.35, 0.2, -0.07, 0.9];
    let response = law
        .plane_stress(&f, PlaneStressSettings::default())
        .unwrap();
    let area = f[0] * f[3] - f[1] * f[2];
    let thickness = response.deformation_gradient[8];
    assert!(thickness < 1.0 && thickness > 0.0);
    close(response.deformation_gradient[2], 0.0, 1.0, 1e-14);
    close(response.deformation_gradient[5], 0.0, 1.0, 1e-14);
    // P33*t = mu*(t^2-1) + lambda*ln(det(F2)*t) = 0.
    close(
        mu * (thickness * thickness - 1.0) + lambda * (area * thickness).ln(),
        0.0,
        mu,
        2e-11,
    );
    // Substituting P33=0 into P gives mu*(F2-t^2*F2^-T).
    let inverse_transpose = [f[3] / area, -f[2] / area, -f[1] / area, f[0] / area];
    for i in 0..4 {
        close(
            response.piola[i],
            mu * (f[i] - thickness * thickness * inverse_transpose[i]),
            mu,
            3e-11,
        );
    }
    assert!(
        law.piola(&[f[0], f[1], 0.0, f[2], f[3], 0.0, 0.0, 0.0, 1.0])
            .unwrap()[8]
            .abs()
            > 1e4
    );
}

#[test]
fn g0_all_original_energies_keep_consistent_reduced_gradients_and_tangents() {
    let f = [1.23, 0.19, -0.04, 0.93];
    let settings = PlaneStressSettings::default();
    for law in laws() {
        let response = law.plane_stress(&f, settings).unwrap();
        let full_stress = law.piola(&response.deformation_gradient).unwrap();
        for i in TRANSVERSE {
            close(full_stress[i], 0.0, 2e5, 2e-10);
        }
        for i in 0..4 {
            assert_eq!(response.piola[i], full_stress[IN_PLANE[i]]);
        }
        for col in 0..4 {
            let (mut plus, mut minus) = (f, f);
            let h = 2e-6;
            plus[col] += h;
            minus[col] -= h;
            let up = law.plane_stress_state(&plus, settings).unwrap();
            let down = law.plane_stress_state(&minus, settings).unwrap();
            close(
                response.piola[col],
                (up.energy_density_pa - down.energy_density_pa) / (2.0 * h),
                2e5,
                2e-7,
            );
            for row in 0..4 {
                close(
                    response.tangent[row][col],
                    (up.piola[row] - down.piola[row]) / (2.0 * h),
                    2e5,
                    3e-7,
                );
                close(
                    response.tangent[row][col],
                    response.tangent[col][row],
                    2e5,
                    1e-12,
                );
            }
        }
    }
}

#[test]
fn g0_tilted_fibres_relax_transverse_shear_as_well_as_thickness() {
    let law = tilted_fibre(1.0);
    let f = [1.3, 0.2, 0.05, 0.94];
    let initial = [f[0], f[1], 0.0, f[2], f[3], 0.0, 0.0, 0.0, 1.0];
    let p0 = law.piola(&initial).unwrap();
    assert!(p0[2].abs() > 1e3);
    let response = law
        .plane_stress(&f, PlaneStressSettings::default())
        .unwrap();
    assert!(response.deformation_gradient[2].abs() > 1e-3);
    let p = law.piola(&response.deformation_gradient).unwrap();
    for i in TRANSVERSE {
        close(p[i], 0.0, 1e5, 1e-10);
    }
    assert!(response.energy_density_pa < law.energy(&initial));
}

fn rotate(a: [f64; 4], cos: f64, sin: f64) -> [f64; 4] {
    [
        cos * a[0] - sin * a[2],
        cos * a[1] - sin * a[3],
        sin * a[0] + cos * a[2],
        sin * a[1] + cos * a[3],
    ]
}

#[test]
fn g3_superposed_in_plane_rotation_preserves_energy_and_rotates_stress_and_relaxation() {
    let f = [1.23, 0.19, -0.04, 0.93];
    let (sin, cos) = 0.47_f64.sin_cos();
    for law in laws() {
        let reference = law
            .plane_stress_state(&f, PlaneStressSettings::default())
            .unwrap();
        let rotated = law
            .plane_stress_state(&rotate(f, cos, sin), PlaneStressSettings::default())
            .unwrap();
        close(
            rotated.energy_density_pa,
            reference.energy_density_pa,
            2e5,
            1e-11,
        );
        for (actual, expected) in rotated
            .piola
            .into_iter()
            .zip(rotate(reference.piola, cos, sin))
        {
            close(actual, expected, 2e5, 1e-10);
        }
        let a = reference.deformation_gradient;
        let b = rotated.deformation_gradient;
        close(b[2], cos * a[2] - sin * a[5], 1.0, 1e-10);
        close(b[5], sin * a[2] + cos * a[5], 1.0, 1e-10);
        close(b[8], a[8], 1.0, 1e-10);
    }
}

#[test]
fn g3_modulus_rescaling_cannot_change_the_relaxed_deformation() {
    let f = [1.3, 0.2, 0.05, 0.94];
    let settings = PlaneStressSettings::default();
    let reference = tilted_fibre(1.0).plane_stress(&f, settings).unwrap();
    for scale in [1e-12, 1e12] {
        let response = tilted_fibre(scale).plane_stress(&f, settings).unwrap();
        for (a, b) in response
            .deformation_gradient
            .into_iter()
            .zip(reference.deformation_gradient)
        {
            close(a, b, 1.0, 1e-10);
        }
        close(
            response.energy_density_pa / scale,
            reference.energy_density_pa,
            2e5,
            1e-10,
        );
        for i in 0..4 {
            close(response.piola[i] / scale, reference.piola[i], 2e5, 1e-10);
            for j in 0..4 {
                close(
                    response.tangent[i][j] / scale,
                    reference.tangent[i][j],
                    2e5,
                    1e-10,
                );
            }
        }
    }
}

#[test]
fn g0_large_bulk_modulus_cannot_hide_a_transverse_traction_residual() {
    let mu = 2e5;
    for bulk_ratio in [1.0, 1e2, 1e4] {
        let law = Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: mu / 2.0,
                c01: 0.0,
                kappa: mu * bulk_ratio,
            },
            4.0,
        )
        .unwrap();
        let response = law
            .plane_stress(&[1.3, 0.2, 0.05, 0.94], PlaneStressSettings::default())
            .unwrap();
        let p = law.piola(&response.deformation_gradient).unwrap();
        for i in TRANSVERSE {
            close(p[i], 0.0, mu, 2e-11);
        }
    }
}

#[test]
fn g0_state_only_and_tangent_paths_share_exactly_the_same_equilibrium() {
    for law in laws() {
        let f = [1.23, 0.19, -0.04, 0.93];
        let a = law
            .plane_stress_state(&f, PlaneStressSettings::default())
            .unwrap();
        let b = law
            .plane_stress(&f, PlaneStressSettings::default())
            .unwrap();
        assert_eq!(a.deformation_gradient, b.deformation_gradient);
        assert_eq!(a.energy_density_pa.to_bits(), b.energy_density_pa.to_bits());
        assert_eq!(a.piola, b.piola);
        assert_eq!(
            a.relative_transverse_residual.to_bits(),
            b.relative_transverse_residual.to_bits()
        );
        assert_eq!(a.iterations, b.iterations);
    }
}

#[test]
fn g0_inverted_states_material_domain_and_exhausted_budgets_refuse() {
    let law = nh(5e5, 0.3, 4.0);
    for f in [
        [-1.0, 0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0, 0.0],
        [f64::NAN, 0.0, 0.0, 1.0],
        [5.0, 0.0, 0.0, 1.0],
    ] {
        assert!(matches!(
            law.plane_stress(&f, PlaneStressSettings::default()),
            Err(MaterialError::State { .. })
        ));
    }
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
            law.plane_stress(&[1.0, 0.0, 0.0, 1.0], settings),
            Err(MaterialError::Parameters { .. })
        ));
    }
    assert!(
        law.plane_stress(
            &[1.35, 0.2, -0.07, 0.9],
            PlaneStressSettings {
                max_iterations: 1,
                ..PlaneStressSettings::default()
            }
        )
        .is_err()
    );
    // In-plane stretches are admissible; the required thickness contraction is
    // not. The reducer must refuse instead of clamping thickness at 1/limit.
    let bounded = nh(5e5, 0.45, 1.2);
    assert!(
        bounded
            .piola(&[1.2, 0.0, 0.0, 0.0, 1.2, 0.0, 0.0, 0.0, 1.0])
            .is_ok()
    );
    assert!(
        bounded
            .plane_stress(&[1.2, 0.0, 0.0, 1.2], PlaneStressSettings::default())
            .is_err()
    );
}
