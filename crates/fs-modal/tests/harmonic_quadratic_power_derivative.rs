//! Original first-order work for resonance refinement, without unused curvature.
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicResponseError, harmonic_quadratic_power_curvature_with_cancel as curvature,
    harmonic_quadratic_power_derivative_with_cancel as observe,
};

#[test]
fn physical_first_order_work_survives_an_unrepresentable_unused_curvature() {
    let v = [C64::from_re(2.0_f64.powi(-500))];
    let first = [C64::from_re(2.0_f64.powi(600))];
    let result = observe(&v, &first, &[1.], || false).unwrap();
    assert_eq!(result.average_power_w, 2.0_f64.powi(-1001));
    assert_eq!(result.angular_derivative, 2.0_f64.powi(100));
    // P_omegaomega=|v_omega|^2 is outside f64, but the two requested
    // physical observations above are finite and must remain available.
    assert!(curvature(&v, &first, &[C64::ZERO], &[1.], || false).is_err());
}

#[test]
fn complex_work_and_original_small_coefficients_preserve_gyro_separation() {
    let v = [C64::new(1., 2.), C64::new(2., -1.)];
    let first = [C64::new(3., -1.), C64::new(-1., 2.)];
    let result = observe(&v, &first, &[2., 8., -6., 5.], || false).unwrap();
    assert_eq!(result.average_power_w, 17.5);
    assert_eq!(result.angular_derivative, -8.);
    assert_eq!(
        result,
        observe(&v, &first, &[2., 1., 1., 5.], || false).unwrap()
    );

    let field = 2.0_f64.powi(500);
    let weak = 2.0_f64.powi(-900);
    let gyro = 2.0_f64.powi(1000);
    let v = [C64::from_re(field); 2];
    let first = [C64::from_re(3. * field); 2];
    let result = observe(&v, &first, &[weak, gyro, -gyro, 2. * weak], || false).unwrap();
    assert_eq!(result.average_power_w, 1.5 * 2.0_f64.powi(100));
    assert_eq!(result.angular_derivative, 9. * 2.0_f64.powi(100));

    let tiny = f64::from_bits(1);
    let result = observe(&v, &first, &[0., tiny, 0., 0.], || false).unwrap();
    assert_eq!(result.average_power_w, 2.0_f64.powi(-75));
    assert_eq!(result.angular_derivative, 3. * 2.0_f64.powi(-74));

    let low = 2.0_f64.powi(-60);
    let v = [C64::ONE, -C64::ONE];
    let result = observe(&v, &v, &[1., 2., low, 1.], || false).unwrap();
    assert_eq!(result.average_power_w, -0.5 * low);
    assert_eq!(result.angular_derivative, -low);
}

#[test]
fn malformed_or_cancelled_first_order_work_is_atomic_and_reusable() {
    let v = [C64::ONE, C64::new(0.5, 0.25)];
    let first = [C64::ZERO; 2];
    let matrix = [2., 0.5, 0.5, 3.];
    let baseline = observe(&v, &first, &matrix, || false).unwrap();
    assert!(observe(&[], &[], &[], || false).is_err());
    assert!(observe(&v, &first[..1], &matrix, || false).is_err());
    assert!(observe(&v, &first, &matrix[..3], || false).is_err());
    assert!(observe(&v, &[C64::ZERO, C64::new(f64::NAN, 0.)], &matrix, || false).is_err());
    assert!(observe(&v, &first, &[2., f64::INFINITY, 0.5, 3.], || false).is_err());
    let mut polls = 0;
    assert!(matches!(
        observe(&v, &first, &matrix, || {
            polls += 1;
            polls == 3
        }),
        Err(HarmonicResponseError::Cancelled)
    ));
    assert_eq!(polls, 3);
    assert_eq!(baseline, observe(&v, &first, &matrix, || false).unwrap());
}
