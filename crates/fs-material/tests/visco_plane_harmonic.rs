//! Original tensor plane-stress phasors, transverse constraint and source limits.
//! Synthetic constitutive fixtures; no measured material qualification claim.

use fs_material::plane_stress_viscosity::MaterialPlane;
use fs_material::visco::{Mandel6, PronyBranch, TensorProny, ViscoError};
use fs_math::{c64::C64, det};

fn outer(diagonal: f64, a: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| a[i] * a[j] + if i == j { diagonal } else { 0.0 })
    })
}

fn anisotropic() -> TensorProny {
    TensorProny::new_mandel(
        outer(5.0, [0.4, -0.5, 0.7, 0.2, -0.3, 0.6]),
        vec![
            PronyBranch::new(outer(0.6, [1.0, -2.0, 0.5, 0.25, -0.75, 1.5]), 0.125).unwrap(),
            PronyBranch::new(outer(0.9, [-0.5, 0.75, 1.0, 1.25, 0.5, -1.0]), 2.0).unwrap(),
        ],
        0.05,
    )
    .unwrap()
    .with_band((0.0, 100.0))
    .unwrap()
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn every_ordered_plane_satisfies_original_full_complex_stress_constraint() {
    let law = anisotropic();
    let e = [
        C64::new(0.006, -0.002),
        C64::new(-0.004, 0.003),
        C64::new(0.005, 0.002),
    ];
    for (plane, inside, outside) in [
        (MaterialPlane::Axes12, [0, 1, 3], [2, 4, 5]),
        (MaterialPlane::Axes21, [1, 0, 3], [2, 4, 5]),
        (MaterialPlane::Axes13, [0, 2, 5], [1, 3, 4]),
        (MaterialPlane::Axes23, [1, 2, 4], [0, 3, 5]),
    ] {
        for omega in [0.25, 1.0, 8.0, 50.0] {
            let prepared = law.prepare_plane_stress_harmonic(plane, omega).unwrap();
            assert_eq!(prepared.plane(), plane);
            assert_eq!(prepared.omega_rad_s(), omega);
            let strain = prepared.recover_strain(&e).unwrap();
            let (storage, loss) = law.complex_modulus(omega).unwrap();
            let stress: [C64; 6] = core::array::from_fn(|i| {
                (0..6).fold(C64::ZERO, |sum, j| {
                    sum + C64::new(storage[i][j], -loss[i][j]) * strain[j]
                })
            });
            let engineering = [1.0, 1.0, core::f64::consts::FRAC_1_SQRT_2];
            for i in 0..3 {
                close(strain[inside[i]], e[i].scale(engineering[i]), 1e-17);
                close(stress[outside[i]], C64::ZERO, 2e-15);
                let reduced = (0..3).fold(C64::ZERO, |sum, j| {
                    sum + C64::new(
                        prepared.storage_pa()[3 * i + j],
                        -prepared.loss_pa()[3 * i + j],
                    ) * e[j]
                });
                close(reduced, stress[inside[i]].scale(engineering[i]), 2e-15);
            }
        }
    }
}

fn isotropic(k: f64, g: f64) -> Mandel6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            if i < 3 && j < 3 {
                k - 2.0 * g / 3.0 + if i == j { 2.0 * g } else { 0.0 }
            } else if i == j {
                2.0 * g
            } else {
                0.0
            }
        })
    })
}

#[test]
fn isotropic_bulk_relaxation_has_the_correct_complex_poisson_strain() {
    let omega = 3.0;
    let tau = 0.4;
    let law = TensorProny::new_mandel(
        isotropic(12.0, 4.0),
        vec![PronyBranch::new(isotropic(9.0, 0.0), tau).unwrap()],
        0.05,
    )
    .unwrap();
    let e = [
        C64::new(0.004, -0.003),
        C64::new(-0.001, 0.002),
        C64::new(0.006, 0.001),
    ];
    let strain = law
        .prepare_plane_stress_harmonic(MaterialPlane::Axes12, omega)
        .unwrap()
        .recover_strain(&e)
        .unwrap();
    let x = omega * tau;
    let k = C64::new(12.0 + 9.0 * x * x / (1.0 + x * x), -9.0 * x / (1.0 + x * x));
    let expected =
        -((k - C64::from_re(8.0 / 3.0)) / (k + C64::from_re(16.0 / 3.0))) * (e[0] + e[1]);
    close(strain[2], expected, 2e-17);
    close(
        strain[3],
        e[2].scale(core::f64::consts::FRAC_1_SQRT_2),
        2e-17,
    );
    assert_eq!(strain[4], C64::ZERO);
    assert_eq!(strain[5], C64::ZERO);
}

#[test]
fn complete_cycle_admission_distinguishes_rotation_from_in_phase_strain() {
    let diagonal =
        core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 2.0 } else { 0.0 }));
    let law = TensorProny::new_mandel(diagonal, vec![], 0.01).unwrap();
    let prepared = law
        .prepare_plane_stress_harmonic(MaterialPlane::Axes12, 1.0)
        .unwrap();
    // Norm is 0.009 throughout a rotating ellipse; complex Frobenius norm
    // would incorrectly reject it. Neither coordinate is independently maximized.
    let circular = [C64::from_re(0.009), C64::new(0.0, 0.009), C64::ZERO];
    assert!(prepared.recover_strain(&circular).is_ok());
    // Both real/imaginary snapshots are below the bound, but their in-phase
    // combination exceeds it between the sampled phases.
    let diagonal_phase = [C64::new(0.008, 0.008), C64::ZERO, C64::ZERO];
    assert!(matches!(
        prepared.recover_strain(&diagonal_phase),
        Err(ViscoError::Parameters { .. })
    ));
    for angle in [0.0, 0.3, 1.0, 2.7] {
        let (sin, cos) = det::sin_cos(angle);
        let phase = C64::new(cos, sin);
        assert!(
            prepared
                .recover_strain(&circular.map(|v| v * phase))
                .is_ok()
        );
        assert!(
            prepared
                .recover_strain(&diagonal_phase.map(|v| v * phase))
                .is_err()
        );
    }
}

#[test]
fn recovered_transverse_strain_remains_inside_original_source_limit() {
    let law = TensorProny::new_mandel(isotropic(100.0, 1.0), vec![], 0.01).unwrap();
    let prepared = law
        .prepare_plane_stress_harmonic(MaterialPlane::Axes12, 1.0)
        .unwrap();
    let in_plane = [C64::from_re(0.006), C64::from_re(0.006), C64::ZERO];
    assert!(det::sqrt(2.0) * 0.006 < 0.01);
    assert!(prepared.recover_strain(&in_plane).is_err());
}

#[test]
fn original_band_real_dc_and_nonfinite_inputs_are_enforced() {
    let law = anisotropic().with_band((1.0, 10.0)).unwrap();
    for omega in [0.0, 0.5, 10.1] {
        assert!(matches!(
            law.prepare_plane_stress_harmonic(MaterialPlane::Axes12, omega),
            Err(ViscoError::OutOfBand { .. })
        ));
    }
    assert!(
        law.prepare_plane_stress_harmonic(MaterialPlane::Axes12, f64::NAN)
            .is_err()
    );
    let dc = anisotropic()
        .prepare_plane_stress_harmonic(MaterialPlane::Axes12, 0.0)
        .unwrap();
    assert_eq!(dc.recover_strain(&[C64::ZERO; 3]).unwrap(), [C64::ZERO; 6]);
    assert!(dc.recover_strain(&[C64::from_re(0.001); 3]).is_ok());
    assert!(dc.recover_strain(&[C64::new(0.0, 0.001); 3]).is_err());
    assert!(
        dc.recover_strain(&[C64::from_re(f64::INFINITY); 3])
            .is_err()
    );
}
