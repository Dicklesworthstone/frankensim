//! G0/G1/G3 checks of original scalar and tensor Maxwell dispersion.
//! Synthetic constitutive laws; no measured-material or transient-memory claim.

use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_material::visco::{GeneralizedMaxwell, ViscoError};

fn identity(scale: f64) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { scale } else { 0.0 }))
}

fn diagonal_plus_outer(diagonal: f64, v: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| v[i] * v[j] + if i == j { diagonal } else { 0.0 })
    })
}

fn tensor(branches: &[(Mandel6, f64)], kelvin: Option<Mandel6>, band: (f64, f64)) -> TensorProny {
    TensorProny::new(
        identity(5.0),
        kelvin,
        branches
            .iter()
            .map(|&(c, tau)| PronyBranch::new(c, tau).unwrap())
            .collect(),
        band,
    )
    .unwrap()
}

fn anisotropic_fixture() -> (TensorProny, Mandel6, Mandel6, Mandel6) {
    let a = diagonal_plus_outer(1.0, [1.0, -2.0, 0.5, 0.25, -0.75, 1.5]);
    let b = diagonal_plus_outer(2.0, [-0.5, 0.75, 1.0, 1.25, 0.5, -1.0]);
    let kelvin = diagonal_plus_outer(0.25, [0.125, -0.25, 0.375, -0.5, 0.625, -0.75]);
    (
        tensor(&[(a, 0.125), (b, 2.0)], Some(kelvin), (0.0, 100.0)),
        a,
        b,
        kelvin,
    )
}

fn close(actual: f64, expected: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    let tolerance = 3e-13 * expected.abs() + 8.0 * f64::from_bits(1);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance {tolerance:e}"
    );
}

// Independent direct formulas are used only at moderate rates and times.
fn branch_slopes(omega: f64, tau: f64) -> [f64; 3] {
    let x = omega * tau;
    let denominator = (1.0 + x * x).powi(2);
    [
        2.0 * tau * x / denominator,
        tau * (1.0 - x * x) / denominator,
        -2.0 * tau * tau * x / denominator,
    ]
}

#[test]
fn scalar_closed_form_slopes_have_signed_frequency_parity_and_dc_limit() {
    let law = GeneralizedMaxwell::new(5.0, vec![(3.0, 0.25), (7.0, 2.0)]).unwrap();
    for omega in [0.0, 0.125, 0.5, 1.0, 4.0, 64.0] {
        let mut expected = [0.0; 2];
        for &(c, tau) in &law.terms {
            let slope = branch_slopes(omega, tau);
            expected[0] += c * slope[0];
            expected[1] += c * slope[1];
        }
        for sign in [-1.0, 1.0] {
            let actual = law.modulus_frequency_derivative(sign * omega).unwrap();
            close(actual.0, sign * expected[0]);
            close(actual.1, expected[1]);
        }
        let positive = law.modulus(omega);
        let negative = law.modulus(-omega);
        close(negative.0, positive.0);
        close(negative.1, -positive.1);
    }
    assert_eq!(law.modulus_frequency_derivative(0.0).unwrap(), (0.0, 14.75));
}

#[test]
fn tensor_slopes_keep_original_anisotropy_kelvin_viscosity_and_symmetry() {
    let (law, a, b, kelvin) = anisotropic_fixture();
    for omega in [0.0, 0.25, 0.5, 2.0, 5.0] {
        let actual = law.complex_modulus_frequency_derivative(omega).unwrap();
        let sa = branch_slopes(omega, 0.125);
        let sb = branch_slopes(omega, 2.0);
        for i in 0..6 {
            for j in 0..6 {
                close(actual.storage_pa_s[i][j], a[i][j] * sa[0] + b[i][j] * sb[0]);
                close(
                    actual.loss_pa_s[i][j],
                    kelvin[i][j] + a[i][j] * sa[1] + b[i][j] * sb[1],
                );
                close(
                    actual.effective_viscosity_pa_s2[i][j],
                    a[i][j] * sa[2] + b[i][j] * sb[2],
                );
                assert_eq!(actual.storage_pa_s[i][j], actual.storage_pa_s[j][i]);
                assert_eq!(actual.loss_pa_s[i][j], actual.loss_pa_s[j][i]);
                assert_eq!(
                    actual.effective_viscosity_pa_s2[i][j],
                    actual.effective_viscosity_pa_s2[j][i]
                );
            }
        }
    }
    let kelvin_only = tensor(&[], Some(kelvin), (0.0, 100.0));
    for omega in [0.0, 3.0, 100.0] {
        let slope = kelvin_only
            .complex_modulus_frequency_derivative(omega)
            .unwrap();
        assert_eq!(slope.storage_pa_s, identity(0.0));
        assert_eq!(slope.loss_pa_s, kelvin);
        assert_eq!(slope.effective_viscosity_pa_s2, identity(0.0));
    }
}

#[test]
fn tensor_frequency_slopes_match_all_entries_of_native_primal_differences() {
    let (law, _, _, _) = anisotropic_fixture();
    for omega in [0.35, 1.7, 4.0] {
        let h = 2e-5 * omega;
        let (plus_storage, plus_loss) = law.complex_modulus(omega + h).unwrap();
        let (minus_storage, minus_loss) = law.complex_modulus(omega - h).unwrap();
        let slope = law.complex_modulus_frequency_derivative(omega).unwrap();
        for i in 0..6 {
            for j in 0..6 {
                let differences = [
                    (plus_storage[i][j] - minus_storage[i][j]) / (2.0 * h),
                    (plus_loss[i][j] - minus_loss[i][j]) / (2.0 * h),
                    (plus_loss[i][j] / (omega + h) - minus_loss[i][j] / (omega - h)) / (2.0 * h),
                ];
                let analytic = [
                    slope.storage_pa_s[i][j],
                    slope.loss_pa_s[i][j],
                    slope.effective_viscosity_pa_s2[i][j],
                ];
                for (difference, derivative) in differences.into_iter().zip(analytic) {
                    assert!(
                        (difference - derivative).abs() <= 2e-7 * (1.0 + derivative.abs()),
                        "omega={omega}, entry=({i},{j}): {difference:e} != {derivative:e}"
                    );
                }
            }
        }
    }
}

#[test]
fn tensor_frequency_admission_uses_the_original_nonnegative_band() {
    let law = tensor(&[(identity(3.0), 0.5)], None, (1.0, 8.0));
    for omega in [1.0, 8.0] {
        assert!(law.complex_modulus_frequency_derivative(omega).is_ok());
    }
    for omega in [-1.0, 0.0, 9.0] {
        assert!(matches!(
            law.complex_modulus_frequency_derivative(omega),
            Err(ViscoError::OutOfBand { .. })
        ));
    }
    for omega in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            law.complex_modulus_frequency_derivative(omega),
            Err(ViscoError::Parameters { .. })
        ));
    }
}

#[test]
fn scalar_frequency_derivative_refuses_invalid_public_parameters() {
    let law = GeneralizedMaxwell::new(1.0, vec![(2.0, 0.5)]).unwrap();
    for omega in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(law.modulus_frequency_derivative(omega).is_err());
    }
    for equilibrium in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut invalid = law.clone();
        invalid.e_inf = equilibrium;
        assert!(invalid.modulus_frequency_derivative(1.0).is_err());
    }
    for term in [
        (-1.0, 0.5),
        (f64::NAN, 0.5),
        (f64::INFINITY, 0.5),
        (2.0, 0.0),
        (2.0, -1.0),
        (2.0, f64::NAN),
        (2.0, f64::INFINITY),
    ] {
        let mut invalid = law.clone();
        invalid.terms[0] = term;
        assert!(invalid.modulus_frequency_derivative(1.0).is_err());
    }
    let zero_branch = GeneralizedMaxwell::new(1.0, vec![(0.0, 1e300)]).unwrap();
    assert_eq!(
        zero_branch.modulus_frequency_derivative(0.0).unwrap(),
        (0.0, 0.0)
    );
}

#[test]
fn finite_primal_and_modulus_slopes_do_not_require_a_representable_viscosity_slope() {
    let scalar = GeneralizedMaxwell::new(1.0, vec![(1e150, 1e200)]).unwrap();
    assert_eq!(scalar.modulus(0.0), (1.0, 0.0));
    assert!(scalar.modulus_frequency_derivative(0.0).is_err());
    let law = tensor(&[(identity(1e150), 1e200)], None, (0.0, 1.0));
    assert_eq!(
        law.complex_modulus(0.0).unwrap(),
        (identity(5.0), identity(0.0))
    );
    assert!(law.complex_modulus_frequency_derivative(0.0).is_err());

    // At omega*tau=1, modulus slopes fit but d(C''/omega)/domega=-tau²/2 does not.
    let tau = 2.0_f64.powi(600);
    let omega = 2.0_f64.powi(-600);
    let scalar = GeneralizedMaxwell::new(1.0, vec![(1.0, tau)]).unwrap();
    assert_eq!(
        scalar.modulus_frequency_derivative(omega).unwrap(),
        (0.5 * tau, 0.0)
    );
    let law = tensor(&[(identity(1.0), tau)], None, (0.0, 1.0));
    assert!(law.complex_modulus(omega).is_ok());
    assert!(law.complex_modulus_frequency_derivative(omega).is_err());
}

#[test]
fn overflowing_rate_time_preserves_physical_loss_and_weak_mutual_slopes() {
    let mut c = identity(1e150);
    c[0][1] = -1e140;
    c[1][0] = c[0][1];
    let law = tensor(&[(c, 1e250)], None, (0.0, 1e100));
    let (storage, loss) = law.complex_modulus(1e100).unwrap();
    let slope = law.complex_modulus_frequency_derivative(1e100).unwrap();
    close(storage[0][1], -1e140);
    close(loss[0][0], 1e-200);
    close(loss[0][1], -1e-210);
    close(slope.loss_pa_s[0][0], -1e-300);
    close(slope.loss_pa_s[0][1], 1e-310);
    assert_eq!(slope.loss_pa_s[0][1], slope.loss_pa_s[1][0]);
    let scalar = GeneralizedMaxwell::new(1.0, vec![(1e150, 1e250)]).unwrap();
    close(scalar.modulus(1e100).1, 1e-200);
    close(
        scalar.modulus_frequency_derivative(1e100).unwrap().1,
        -1e-300,
    );
}

#[test]
fn underflowing_squared_rate_time_preserves_weak_mutual_storage_and_slopes() {
    let mut c = identity(1e150);
    c[0][1] = -1e140;
    c[1][0] = c[0][1];
    let law = tensor(&[(c, 1e-100)], None, (0.0, 1.0));
    let (storage, loss) = law.complex_modulus(1e-100).unwrap();
    let slope = law.complex_modulus_frequency_derivative(1e-100).unwrap();
    close(storage[0][1], -1e-260);
    close(loss[0][1], -1e-60);
    close(slope.storage_pa_s[0][1], -2e-160);
    close(slope.loss_pa_s[0][1], -1e40);
    close(slope.effective_viscosity_pa_s2[0][1], 2e-260);
    let scalar = GeneralizedMaxwell::new(1e-260, vec![(1e150, 1e-100)]).unwrap();
    close(scalar.modulus(1e-100).0, 1e-260 + 1e-250);
    close(
        scalar.modulus_frequency_derivative(1e-100).unwrap().0,
        2e-150,
    );
}

#[test]
fn minimum_subnormal_relaxation_time_retains_admitted_tensor_response() {
    let omega = 2.0_f64.powi(500);
    let tau = f64::from_bits(1);
    let mut c = identity(2.0_f64.powi(500));
    c[0][1] = -2.0_f64.powi(498);
    c[1][0] = c[0][1];
    let law = tensor(&[(c, tau)], None, (0.0, omega));
    let (storage, loss) = law.complex_modulus(omega).unwrap();
    let slope = law.complex_modulus_frequency_derivative(omega).unwrap();
    close(storage[0][1], -2.0_f64.powi(-650));
    close(loss[0][1], -2.0_f64.powi(-76));
    close(slope.loss_pa_s[0][1], -2.0_f64.powi(-576));
    let scalar = GeneralizedMaxwell::new(1.0, vec![(2.0_f64.powi(500), tau)]).unwrap();
    close(
        scalar.modulus_frequency_derivative(omega).unwrap().1,
        2.0_f64.powi(-574),
    );
}

#[test]
fn compensated_relaxation_slopes_retain_a_small_net_loss_derivative() {
    // The first and last loss slopes cancel exactly; their magnitude hides
    // the middle contribution in an ordinary left-to-right f64 sum.
    let large = 2.0_f64.powi(60);
    let scalar =
        GeneralizedMaxwell::new(1.0, vec![(large, 0.5), (1.0, 0.5), (large, 2.0)]).unwrap();
    close(scalar.modulus_frequency_derivative(1.0).unwrap().1, 0.24);
    let law = tensor(
        &[
            (identity(large), 0.5),
            (identity(1.0), 0.5),
            (identity(large), 2.0),
        ],
        None,
        (0.0, 2.0),
    );
    let slope = law.complex_modulus_frequency_derivative(1.0).unwrap();
    for i in 0..6 {
        close(slope.loss_pa_s[i][i], 0.24);
    }
}
