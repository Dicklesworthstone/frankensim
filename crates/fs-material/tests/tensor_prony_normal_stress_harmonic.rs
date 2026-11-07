//! Independent one-normal constraints for the original anisotropic material.
//! These are manufactured constitutive fixtures, not measured material claims.

use fs_material::plane_stress_viscosity::MaterialPlane;
use fs_material::visco::{Mandel6, PronyBranch, TensorProny, ViscoError};
use fs_math::c64::C64;

fn matrix(diagonal: [f64; 6], vector: [f64; 6], scale: f64) -> Mandel6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            scale * (vector[i] * vector[j] + if i == j { diagonal[i] } else { 0.0 })
        })
    })
}

fn source() -> TensorProny {
    let equilibrium = matrix(
        [17.0, 13.0, 11.0, 7.0, 5.0, 9.0],
        [1.3, -0.7, 0.9, 0.4, -0.5, 0.6],
        1e6,
    );
    let kelvin = matrix(
        [0.7, 0.8, 0.9, 0.5, 0.3, 0.4],
        [0.3, 0.9, -0.4, 0.7, 0.8, -0.6],
        1e3,
    );
    let branch = matrix(
        [1.1, 0.9, 0.7, 0.5, 0.4, 0.3],
        [0.7, 0.4, -0.8, 0.9, -0.3, 0.5],
        1e6,
    );
    let second = matrix(
        [0.2, 0.1, 0.3, 0.4, 0.2, 0.3],
        [-0.2, 0.8, 0.3, -0.5, 0.9, 0.4],
        1e6,
    );
    TensorProny::new(
        equilibrium,
        Some(kelvin),
        vec![
            PronyBranch::new(branch, 0.014).unwrap(),
            PronyBranch::new(second, 0.23).unwrap(),
        ],
        (0.0, 500.0),
    )
    .unwrap()
}

const PLANES: [MaterialPlane; 4] = [
    MaterialPlane::Axes12,
    MaterialPlane::Axes21,
    MaterialPlane::Axes13,
    MaterialPlane::Axes23,
];

fn coordinates(plane: MaterialPlane) -> ([usize; 5], usize) {
    match plane {
        MaterialPlane::Axes12 => ([0, 1, 3, 5, 4], 2),
        MaterialPlane::Axes21 => ([1, 0, 3, 4, 5], 2),
        MaterialPlane::Axes13 => ([0, 2, 5, 3, 4], 1),
        MaterialPlane::Axes23 => ([1, 2, 4, 3, 5], 0),
    }
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual {actual:?}, expected {expected:?}"
    );
}

fn complex(storage: &Mandel6, loss: &Mandel6) -> [[C64; 6]; 6] {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| C64::new(storage[i][j], loss[i][j]))
    })
}

#[test]
fn independent_scalar_schur_retains_every_shear_and_cross_coefficient() {
    let law = source();
    let engineering = [1.0, 1.0, 0.5_f64.sqrt(), 0.5_f64.sqrt(), 0.5_f64.sqrt()];
    for plane in PLANES {
        let (inside, normal) = coordinates(plane);
        for omega in [0.0, 0.01, 37.0, 300.0] {
            let prepared = law.prepare_normal_stress_harmonic(plane, omega).unwrap();
            assert!(core::ptr::eq(prepared.source(), &law));
            assert_eq!(prepared.plane(), plane);
            assert_eq!(prepared.omega_rad_s(), omega);
            let (storage, loss) = law.complex_modulus(omega).unwrap();
            let c = complex(&storage, &loss);
            for i in 0..5 {
                for j in 0..5 {
                    let expected = (c[inside[i]][inside[j]]
                        - c[inside[i]][normal] * (c[normal][inside[j]] / c[normal][normal]))
                        .scale(engineering[i] * engineering[j]);
                    let actual = C64::new(prepared.storage_pa()[5 * i + j], prepared.loss_pa()[5 * i + j]);
                    close(actual, expected, 3e-13);
                    assert_eq!(prepared.storage_pa()[5 * i + j], prepared.storage_pa()[5 * j + i]);
                    assert_eq!(prepared.loss_pa()[5 * i + j], prepared.loss_pa()[5 * j + i]);
                }
            }
            // General anisotropy couples retained in-plane and transverse shear.
            assert!((0..3).any(|i| (3..5).any(|j| prepared.storage_pa()[5 * i + j] != 0.0)));
            assert!(prepared.storage_pa()[18] > 0.0);
            assert!(prepared.storage_pa()[24] > 0.0);
        }
    }
}

#[test]
fn physical_negative_time_phase_satisfies_only_the_original_normal_constraint() {
    let law = source();
    let omega = 37.0;
    let (storage, loss) = law.complex_modulus(omega).unwrap();
    let engineering = [1.0, 1.0, 0.5_f64.sqrt(), 0.5_f64.sqrt(), 0.5_f64.sqrt()];
    let input = [
        C64::new(2e-4, -3e-5),
        C64::new(-7e-5, 1e-4),
        C64::new(6e-5, 2e-5),
        C64::new(1e-4, -4e-5),
        C64::new(-9e-5, 7e-5),
    ];
    for plane in PLANES {
        let prepared = law.prepare_normal_stress_harmonic(plane, omega).unwrap();
        let strain = prepared.recover_strain(&input).unwrap();
        let (inside, normal) = coordinates(plane);
        let stress: [C64; 6] = core::array::from_fn(|i| {
            (0..6).fold(C64::ZERO, |sum, j| {
                sum + C64::new(storage[i][j], -loss[i][j]) * strain[j]
            })
        });
        let scale = stress.iter().map(|z| z.abs()).fold(1.0_f64, f64::max);
        assert!(stress[normal].abs() <= 3e-14 * scale);
        assert!(stress[inside[3]].abs() > 1.0);
        assert!(stress[inside[4]].abs() > 1.0);
        for i in 0..5 {
            close(strain[inside[i]], input[i].scale(engineering[i]), 1e-16);
            let reduced = (0..5).fold(C64::ZERO, |sum, j| {
                sum + C64::new(prepared.storage_pa()[5 * i + j], -prepared.loss_pa()[5 * i + j]) * input[j]
            });
            close(reduced, stress[inside[i]].scale(engineering[i]), 4e-13);
        }
        let full_work = (0..6).fold(C64::ZERO, |sum, i| sum + strain[i].conj() * stress[i]);
        let reduced_work = (0..5).fold(C64::ZERO, |sum, i| {
            let stress = (0..5).fold(C64::ZERO, |part, j| {
                part + C64::new(prepared.storage_pa()[5 * i + j], -prepared.loss_pa()[5 * i + j]) * input[j]
            });
            sum + input[i].conj() * stress
        });
        close(reduced_work, full_work, 4e-14);
        assert!(full_work.re > 0.0);
        assert!(full_work.im < 0.0);
    }
}

#[test]
fn exact_tangent_matches_independent_scalar_quotient_and_frequency_differences() {
    let law = source();
    let engineering = [1.0, 1.0, 0.5_f64.sqrt(), 0.5_f64.sqrt(), 0.5_f64.sqrt()];
    for plane in PLANES {
        let (inside, normal) = coordinates(plane);
        for omega in [0.0, 0.01, 37.0, 500.0] {
            let prepared = law.prepare_normal_stress_harmonic(plane, omega).unwrap();
            let (ds, dl) = prepared.frequency_derivative().unwrap();
            let (storage, loss) = law.complex_modulus(omega).unwrap();
            let derivative = law.complex_modulus_frequency_derivative(omega).unwrap();
            let c = complex(&storage, &loss);
            let dc = complex(&derivative.storage_pa_s, &derivative.loss_pa_s);
            for i in 0..5 {
                for j in 0..5 {
                    let (a, b) = (inside[i], inside[j]);
                    let ratio = c[normal][b] / c[normal][normal];
                    let ratio_rate = (dc[normal][b] - dc[normal][normal] * ratio) / c[normal][normal];
                    let expected = (dc[a][b] - dc[a][normal] * ratio - c[a][normal] * ratio_rate)
                        .scale(engineering[i] * engineering[j]);
                    close(C64::new(ds[5 * i + j], dl[5 * i + j]), expected, 2e-12);
                }
            }
        }
        let omega = 37.0;
        let h = 1e-3;
        let center = law.prepare_normal_stress_harmonic(plane, omega).unwrap();
        let (ds, dl) = center.frequency_derivative().unwrap();
        let left = law.prepare_normal_stress_harmonic(plane, omega - h).unwrap();
        let right = law.prepare_normal_stress_harmonic(plane, omega + h).unwrap();
        for k in 0..25 {
            let difference = C64::new(
                (right.storage_pa()[k] - left.storage_pa()[k]) / (2.0 * h),
                (right.loss_pa()[k] - left.loss_pa()[k]) / (2.0 * h),
            );
            close(C64::new(ds[k], dl[k]), difference, 2e-7);
        }
    }
}

#[test]
fn isotropic_shear_uses_the_original_shear_modulus_without_a_correction_factor() {
    let lambda = 3e9;
    let shear = 2e9;
    let mut c = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = lambda + if i == j { 2.0 * shear } else { 0.0 };
        }
    }
    for i in 3..6 {
        c[i][i] = 2.0 * shear;
    }
    let law = TensorProny::new(c, None, vec![], (0.0, 100.0)).unwrap();
    let material = law.prepare_normal_stress_harmonic(MaterialPlane::Axes12, 10.0).unwrap();
    for i in 2..5 {
        close(C64::from_re(material.storage_pa()[5 * i + i]), C64::from_re(shear), 3e-15);
    }
    assert_eq!(material.loss_pa(), &[0.0; 25]);
    assert_eq!(material.frequency_derivative().unwrap(), ([0.0; 25], [0.0; 25]));
    // The in-plane part agrees with the existing full plane-stress owner
    // for this uncoupled source, while transverse shear remains retained.
    let (plane, _) = law.plane_stress_complex(MaterialPlane::Axes12, 10.0).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            close(C64::from_re(material.storage_pa()[5 * i + j]), C64::from_re(plane[3 * i + j]), 3e-15);
        }
    }
}

#[test]
fn full_normal_and_shear_strain_cycles_retain_the_original_source_domain() {
    let mut c = core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
    c[0][0] = 5.0;
    c[0][2] = 2.0;
    c[2][0] = 2.0;
    let law = TensorProny::new_mandel(c, vec![], 0.01).unwrap().with_band((0.0, 2.0)).unwrap();
    let prepared = law.prepare_normal_stress_harmonic(MaterialPlane::Axes12, 1.0).unwrap();
    let mut strain = [C64::ZERO; 5];
    strain[0] = C64::from_re(0.004);
    assert!(prepared.recover_strain(&strain).is_ok());
    strain[0] = C64::from_re(0.006);
    assert!(prepared.recover_strain(&strain).is_err()); // eliminated normal exceeds full norm
    strain[0] = C64::new(0.004, 0.004);
    assert!(prepared.recover_strain(&strain).is_err()); // phase snapshots alone would admit
    strain[0] = C64::from_re(0.003);
    strain[1] = C64::new(0.0, 0.008);
    assert!(prepared.recover_strain(&strain).is_ok()); // actual ellipse, not sum of snapshot norms
    strain = [C64::ZERO; 5];
    strain[3] = C64::from_re(0.02);
    assert!(prepared.recover_strain(&strain).is_err()); // original physical transverse shear
    strain[3] = C64::new(f64::NAN, 0.0);
    assert!(prepared.recover_strain(&strain).is_err());
    assert!(matches!(
        law.prepare_normal_stress_harmonic(MaterialPlane::Axes12, 3.0),
        Err(ViscoError::OutOfBand { .. })
    ));
    let dc = law.prepare_normal_stress_harmonic(MaterialPlane::Axes12, 0.0).unwrap();
    let mut imaginary = [C64::ZERO; 5];
    imaginary[4] = C64::new(0.0, 1e-4);
    assert!(dc.recover_strain(&imaginary).is_err());
}

#[test]
fn unrelated_huge_normal_loss_preserves_weak_retained_dissipation_and_slope() {
    let equilibrium = core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
    let mut kelvin = [[0.0; 6]; 6];
    kelvin[0][0] = 1e-100;
    kelvin[2][2] = 1e300;
    let law = TensorProny::new(equilibrium, Some(kelvin), vec![], (0.0, 2.0)).unwrap();
    let prepared = law.prepare_normal_stress_harmonic(MaterialPlane::Axes12, 1.0).unwrap();
    assert_eq!(prepared.loss_pa()[0], 1e-100);
    let (storage_rate, loss_rate) = prepared.frequency_derivative().unwrap();
    assert_eq!(storage_rate, [0.0; 25]);
    assert_eq!(loss_rate[0], 1e-100);
    assert_eq!(loss_rate[24], 0.0);
}
